import argparse
from email import policy
from email.parser import BytesParser
import json
import math
import re
import tempfile
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, quote, urlsplit
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parent
CSV_FILE = ROOT / "RCSB_MASTER_DATA_EXPANDED_ANNOTATED.csv"
MAX_PDB_UPLOAD_BYTES = 25 * 1024 * 1024
COHERENT_PAIR_METRICS = (
    "diffusion_length_angstrom",
    "total_pairs",
    "coherent_pairs",
    "coherent_pair_fraction",
)
NUMERICAL_FIELDS = {
    "FRET_Diffusion_Length_Angstroms",
    "DET_Diffusion_Length_Angstroms",
    "TaxID",
}


def _download_structure(source, identifier):
    identifier = identifier.strip().upper()
    if source == "rcsb":
        if not re.fullmatch(r"[A-Z0-9]{4}", identifier):
            raise ValueError("Enter a valid four-character RCSB PDB ID.")
        structure_url = f"https://files.rcsb.org/download/{identifier}.pdb"
        filename = f"{identifier}.pdb"
    elif source == "alphafold":
        accession_match = re.fullmatch(
            r"(?:AF-)?([A-Z][0-9][A-Z0-9]{3,8}[0-9])(?:-F([0-9]+))?"
            r"(?:-MODEL_V[0-9]+)?(?:\.PDB)?",
            identifier,
        )
        if not accession_match:
            raise ValueError(
                "Enter a UniProt accession or AlphaFold ID such as AF-P05067-F1."
            )
        accession = accession_match.group(1)
        metadata_request = Request(
            f"https://alphafold.ebi.ac.uk/api/prediction/{quote(accession)}",
            headers={"User-Agent": "Plexcal/1.0"},
        )
        with urlopen(metadata_request, timeout=30) as response:
            metadata = json.loads(response.read(1024 * 1024))
        if not isinstance(metadata, list) or not metadata:
            raise ValueError(f"No AlphaFold structure found for {accession}.")
        requested_fragment = accession_match.group(2)
        prediction = next(
            (
                item
                for item in metadata
                if isinstance(item, dict)
                and (
                    requested_fragment is None
                    or item.get("entryId", "").upper().endswith(
                        f"-F{requested_fragment}"
                    )
                )
            ),
            None,
        )
        if prediction is None:
            raise ValueError(f"No AlphaFold structure found for {identifier}.")
        structure_url = prediction.get("pdbUrl")
        parsed_structure_url = (
            urlsplit(structure_url) if isinstance(structure_url, str) else None
        )
        if (
            parsed_structure_url is None
            or parsed_structure_url.scheme != "https"
            or parsed_structure_url.hostname != "alphafold.ebi.ac.uk"
            or not parsed_structure_url.path.startswith("/files/")
            or not parsed_structure_url.path.lower().endswith(".pdb")
        ):
            raise ValueError(f"AlphaFold returned no valid PDB file for {accession}.")
        prediction_id = prediction.get("entryId")
        if not isinstance(prediction_id, str):
            raise ValueError(f"AlphaFold returned an invalid entry ID for {accession}.")
        identifier = prediction_id.upper()
        filename = f"{identifier}.pdb"
    else:
        raise ValueError("Choose upload, RCSB PDB ID, or AlphaFold structure ID.")

    structure_request = Request(
        structure_url, headers={"User-Agent": "Plexcal/1.0"}
    )
    with urlopen(structure_request, timeout=60) as response:
        pdb_bytes = response.read(MAX_PDB_UPLOAD_BYTES + 1)
    if len(pdb_bytes) > MAX_PDB_UPLOAD_BYTES:
        raise ValueError("The downloaded PDB file exceeds the 25 MiB limit.")
    if not any(
        line.startswith((b"ATOM  ", b"HETATM")) for line in pdb_bytes.splitlines()
    ):
        raise ValueError(f"No PDB atom records were downloaded for {identifier}.")
    return pdb_bytes, identifier, filename


def _request_json(url, payload=None):
    request_data = (
        json.dumps(payload).encode("utf-8") if payload is not None else None
    )
    headers = {"User-Agent": "Plexcal/1.0", "Accept": "application/json"}
    if request_data is not None:
        headers["Content-Type"] = "application/json"
    request = Request(url, data=request_data, headers=headers)
    with urlopen(request, timeout=30) as response:
        response_body = response.read(5 * 1024 * 1024 + 1)
    if len(response_body) > 5 * 1024 * 1024:
        raise ValueError("The structure metadata response exceeded 5 MiB.")
    return json.loads(response_body)


def _get_structure_metadata(source, identifier):
    identifier = identifier.strip().upper()
    if source == "rcsb":
        if not re.fullmatch(r"[A-Z0-9]{4}", identifier):
            raise ValueError("Enter a valid four-character RCSB PDB ID.")
        query = (
            'query { entry(entry_id: "%s") { '
            "struct { title } "
            "rcsb_entry_info { experimental_method resolution_combined "
            "polymer_entity_count deposited_atom_count } "
            "rcsb_primary_citation { title rcsb_authors journal_abbrev year "
            "pdbx_database_id_DOI pdbx_database_id_PubMed } "
            "polymer_entities { rcsb_entity_source_organism { "
            "ncbi_scientific_name ncbi_common_names } } } }"
        ) % identifier
        response = _request_json(
            "https://data.rcsb.org/graphql", {"query": query}
        )
        entry = response.get("data", {}).get("entry")
        if not isinstance(entry, dict):
            raise ValueError(f"No RCSB structure metadata found for {identifier}.")

        info = entry.get("rcsb_entry_info") or {}
        citation = entry.get("rcsb_primary_citation") or {}
        organisms = []
        for entity in entry.get("polymer_entities") or []:
            for organism in (entity or {}).get("rcsb_entity_source_organism") or []:
                name = organism.get("ncbi_scientific_name")
                if name and name not in organisms:
                    organisms.append(name)

        doi = citation.get("pdbx_database_id_DOI")
        pubmed_id = citation.get("pdbx_database_id_PubMed")
        resolution = info.get("resolution_combined") or []
        return {
            "source": "RCSB PDB",
            "identifier": identifier,
            "title": (entry.get("struct") or {}).get("title"),
            "organism": ", ".join(organisms) or None,
            "method": info.get("experimental_method"),
            "resolution_angstrom": resolution[0] if resolution else None,
            "polymer_entity_count": info.get("polymer_entity_count"),
            "atom_count": info.get("deposited_atom_count"),
            "publication_title": citation.get("title"),
            "authors": citation.get("rcsb_authors") or [],
            "journal": citation.get("journal_abbrev"),
            "publication_year": citation.get("year"),
            "doi": doi,
            "pubmed_id": pubmed_id,
            "record_url": f"https://www.rcsb.org/structure/{identifier}",
            "publication_url": (
                f"https://doi.org/{quote(str(doi), safe='/')}" if doi else None
            ),
            "pubmed_url": (
                f"https://pubmed.ncbi.nlm.nih.gov/{quote(str(pubmed_id), safe='')}/"
                if pubmed_id
                else None
            ),
        }

    if source == "alphafold":
        accession_match = re.fullmatch(
            r"(?:AF-)?([A-Z][0-9][A-Z0-9]{3,8}[0-9])(?:-F([0-9]+))?"
            r"(?:-MODEL_V[0-9]+)?(?:\.PDB)?",
            identifier,
        )
        if not accession_match:
            raise ValueError(
                "Enter a UniProt accession or AlphaFold ID such as AF-P05067-F1."
            )
        accession = accession_match.group(1)
        predictions = _request_json(
            f"https://alphafold.ebi.ac.uk/api/prediction/{quote(accession)}"
        )
        requested_fragment = accession_match.group(2)
        prediction = next(
            (
                item
                for item in predictions
                if isinstance(item, dict)
                and (
                    requested_fragment is None
                    or item.get("entryId", "").upper().endswith(
                        f"-F{requested_fragment}"
                    )
                )
            ),
            None,
        ) if isinstance(predictions, list) else None
        if prediction is None:
            raise ValueError(f"No AlphaFold structure found for {identifier}.")

        entry_id = prediction.get("entryId")
        if not isinstance(entry_id, str):
            raise ValueError(f"AlphaFold returned an invalid entry ID for {identifier}.")
        uniprot_accession = prediction.get("uniprotAccession") or accession
        return {
            "source": "AlphaFold Protein Structure Database",
            "identifier": entry_id,
            "title": prediction.get("uniprotDescription"),
            "organism": prediction.get("organismScientificName"),
            "method": "AlphaFold predicted model",
            "resolution_angstrom": None,
            "mean_plddt": prediction.get("globalMetricValue"),
            "plddt_very_low_fraction": prediction.get("fractionPlddtVeryLow"),
            "plddt_low_fraction": prediction.get("fractionPlddtLow"),
            "plddt_confident_fraction": prediction.get("fractionPlddtConfident"),
            "plddt_very_high_fraction": prediction.get("fractionPlddtVeryHigh"),
            "model_updated": prediction.get("modelCreatedDate"),
            "sequence_version_date": prediction.get("sequenceVersionDate"),
            "uniprot_accession": uniprot_accession,
            "record_url": (
                "https://alphafold.ebi.ac.uk/entry/"
                + quote(entry_id, safe="-")
            ),
            "uniprot_url": (
                "https://www.uniprot.org/uniprotkb/"
                + quote(uniprot_accession, safe="")
            ),
            "publication_title": (
                "The AlphaFold Protein Structure Database in 2021"
            ),
            "authors": ["Varadi et al."],
            "journal": "Nucleic Acids Research",
            "publication_year": 2022,
            "doi": "10.1093/nar/gkab1061",
            "publication_url": "https://doi.org/10.1093/nar/gkab1061",
        }

    raise ValueError("Choose an RCSB or AlphaFold structure ID.")


def _finite_metric(result, field):
    value = result.get(field)
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    if field in {"total_pairs", "coherent_pairs"} and number.is_integer():
        return int(number)
    return number


def clean_record(row):
    record = {}
    for key, raw_value in row.items():
        value = (raw_value or "").strip()
        if key not in NUMERICAL_FIELDS:
            record[key] = value
        elif not value:
            record[key] = None
        else:
            try:
                number = float(value)
            except ValueError:
                record[key] = value
            else:
                if not math.isfinite(number):
                    record[key] = value
                elif key == "TaxID":
                    record[key] = int(number)
                else:
                    record[key] = number
    return record


class PDBRequestHandler(SimpleHTTPRequestHandler):
    def send_json(self, status, payload):
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        request = urlsplit(self.path)
        if request.path == "/api/status":
            self.send_json(
                200,
                {
                    "fretCalculationAvailable": True,
                    "detCalculationAvailable": True,
                },
            )
            return

        if request.path == "/api/structure":
            query = parse_qs(request.query)
            source = query.get("source", [""])[0].strip().lower()
            identifier = query.get("id", [""])[0].strip()
            try:
                pdb_bytes, _, filename = _download_structure(
                    source, identifier
                )
            except ValueError as error:
                self.send_json(400, {"error": str(error)})
                return
            except (
                HTTPError,
                URLError,
                TimeoutError,
                OSError,
                json.JSONDecodeError,
            ) as error:
                self.send_json(
                    502,
                    {"error": f"Could not download structure {identifier}: {error}"},
                )
                return

            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-pdb")
            self.send_header("Content-Length", str(len(pdb_bytes)))
            self.send_header(
                "Content-Disposition", f'inline; filename="{filename}"'
            )
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(pdb_bytes)
            return

        if request.path == "/api/structure-metadata":
            query = parse_qs(request.query)
            source = query.get("source", [""])[0].strip().lower()
            identifier = query.get("id", [""])[0].strip()
            try:
                metadata = _get_structure_metadata(source, identifier)
            except ValueError as error:
                self.send_json(400, {"error": str(error)})
                return
            except (
                HTTPError,
                URLError,
                TimeoutError,
                OSError,
                json.JSONDecodeError,
            ) as error:
                self.send_json(
                    502,
                    {"error": f"Could not retrieve metadata for {identifier}: {error}"},
                )
                return
            self.send_json(200, metadata)
            return

        super().do_GET()

    def do_POST(self):
        request = urlsplit(self.path)
        if request.path not in {"/api/calculate", "/api/calculate-fret"}:
            self.send_json(404, {"error": "API endpoint not found."})
            return
        legacy_fret_endpoint = request.path == "/api/calculate-fret"

        try:
            content_length = int(self.headers.get("Content-Length", ""))
        except ValueError:
            self.send_json(411, {"error": "A valid Content-Length is required."})
            return

        if content_length <= 0:
            self.send_json(400, {"error": "Upload a non-empty PDB file."})
            return
        if content_length > MAX_PDB_UPLOAD_BYTES + 1024 * 1024:
            self.send_json(413, {"error": "The upload exceeds the 25 MiB limit."})
            return
        if not self.headers.get("Content-Type", "").lower().startswith(
            "multipart/form-data"
        ):
            self.send_json(
                415, {"error": "Upload the PDB file as multipart form data."}
            )
            return

        body = self.rfile.read(content_length)
        if len(body) != content_length:
            self.send_json(400, {"error": "The upload was incomplete."})
            return

        try:
            message = BytesParser(policy=policy.default).parsebytes(
                b"Content-Type: "
                + self.headers["Content-Type"].encode("ascii")
                + b"\r\nMIME-Version: 1.0\r\n\r\n"
                + body
            )
        except (UnicodeEncodeError, ValueError) as error:
            self.send_json(400, {"error": f"Could not read the upload: {error}"})
            return
        if not message.is_multipart():
            self.send_json(400, {"error": "The upload form is malformed."})
            return

        upload_part = None
        source = "upload"
        identifier = ""
        for part in message.iter_parts():
            field_name = part.get_param("name", header="content-disposition")
            if field_name == "pdbFile":
                upload_part = part
            elif field_name == "inputSource":
                source = (part.get_payload(decode=True) or b"").decode(
                    "ascii", errors="replace"
                ).strip().lower()
            elif field_name == "structureId":
                identifier = (part.get_payload(decode=True) or b"").decode(
                    "ascii", errors="replace"
                ).strip()

        if legacy_fret_endpoint and source != "upload":
            self.send_json(400, {"error": "The legacy FRET endpoint accepts uploads only."})
            return

        if source == "upload":
            if upload_part is None:
                self.send_json(400, {"error": "Choose a PDB file to calculate."})
                return

            pdb_bytes = upload_part.get_payload(decode=True) or b""
            if not pdb_bytes:
                self.send_json(400, {"error": "The selected file is empty."})
                return
            if len(pdb_bytes) > MAX_PDB_UPLOAD_BYTES:
                self.send_json(
                    413, {"error": "The PDB file must be 25 MiB or smaller."}
                )
                return
            if not any(
                line.startswith((b"ATOM  ", b"HETATM"))
                for line in pdb_bytes.splitlines()
            ):
                self.send_json(
                    400, {"error": "The uploaded file has no PDB atom records."}
                )
                return
            original_name = Path(
                upload_part.get_filename() or "uploaded.pdb"
            ).name
            identifier = original_name
        elif source in {"rcsb", "alphafold"}:
            if upload_part is not None:
                self.send_json(
                    400, {"error": "Provide either a PDB upload or a structure ID."}
                )
                return
            if not identifier:
                self.send_json(400, {"error": "Enter a structure ID to calculate."})
                return
            try:
                pdb_bytes, identifier, original_name = _download_structure(
                    source, identifier
                )
            except ValueError as error:
                self.send_json(400, {"error": str(error)})
                return
            except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
                self.send_json(
                    502, {"error": f"Could not download the requested PDB structure: {error}"}
                )
                return
        else:
            self.send_json(400, {"error": "Choose upload, RCSB, or AlphaFold input."})
            return

        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(
                prefix="pdb-calculation-", suffix=".pdb", delete=False
            ) as temporary_file:
                temporary_file.write(pdb_bytes)
                temporary_path = Path(temporary_file.name)

            from process_single_pdb import process_single_pdb

            fret_result = process_single_pdb(temporary_path)
            if not isinstance(fret_result, dict):
                raise TypeError("The FRET processor must return a metric dictionary.")

            if legacy_fret_endpoint:
                fret_length = _finite_metric(fret_result, "diffusion_length_angstrom")
                if fret_length is None:
                    self.send_json(
                        422,
                        {
                            "error": (
                                "The FRET calculation failed for this PDB. Check "
                                "that aromatic residue ring atoms and coordinates are present."
                            )
                        },
                    )
                else:
                    self.send_json(
                        200,
                        {
                            "fileName": original_name,
                            "diffusionLengthAngstroms": fret_length,
                        },
                    )
                return

            from det_code_clean import process_single_pdb as process_det_pdb

            det_result = process_det_pdb(temporary_path)
            if not isinstance(det_result, dict):
                raise TypeError("The DET processor must return a metric dictionary.")

            def metrics_for(result, method):
                metrics = {
                    field: _finite_metric(result, field)
                    for field in COHERENT_PAIR_METRICS
                }
                error = None
                if metrics["diffusion_length_angstrom"] is None:
                    error = (
                        "Check that aromatic residue ring atoms and coordinates "
                        "are present."
                        if method == "fret"
                        else "Check that the structure contains valid TRP, TYR, "
                        "or PHE residue coordinates."
                    )
                if error:
                    metrics["error"] = error
                return metrics

            fret_metrics = metrics_for(fret_result, "fret")
            det_metrics = metrics_for(det_result, "det")
            response = {
                "input": {
                    "source": source,
                    "identifier": identifier,
                    "filename": original_name,
                },
                "fret": fret_metrics,
                "det": det_metrics,
            }
            failed = (
                fret_metrics["diffusion_length_angstrom"] is None
                or det_metrics["diffusion_length_angstrom"] is None
            )
            self.send_json(422 if failed else 200, response)
        except (OSError, ImportError, TypeError, ValueError) as error:
            self.send_json(
                500, {"error": f"Could not run the requested calculation: {error}"}
            )
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(
        description="Serve the local FRET/DET structure analysis application."
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()

    handler = partial(PDBRequestHandler, directory=str(ROOT))
    server = ThreadingHTTPServer((args.host, args.port), handler)
    print(f"Serving Plexcal structure analysis at http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
