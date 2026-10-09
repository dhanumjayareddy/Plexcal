const form = document.getElementById("fretUploadForm");
const statusMessage = document.getElementById("fretUploadStatus");
const resultSection = document.getElementById("fretUploadResult");
const resultError = document.getElementById("fretUploadResultError");
const fileInput = document.getElementById("pdbFile");
const sourceInput = document.getElementById("structureSource");
const identifierInput = document.getElementById("structureIdentifier");
const identifierLabel = document.getElementById("structureIdentifierLabel");
const uploadFields = document.getElementById("uploadSourceFields");
const identifierFields = document.getElementById("identifierSourceFields");
const calculateButton = document.getElementById("calculateButton");
const previewButton = document.getElementById("loadStructurePreview");
const viewer = document.getElementById("structureViewer");
const viewerStatus = document.getElementById("structureViewerStatus");
const expandButton = document.getElementById("expandStructureViewer");
const visibilityControls = document.getElementById("structureVisibilityControls");
const metadataStatus = document.getElementById("structureMetadataStatus");
const metadataList = document.getElementById("structureMetadata");
const metadataLinks = document.getElementById("structureMetadataLinks");
const metadataSourceBadge = document.getElementById("metadataSourceBadge");
const metadataHint = document.getElementById("structureIdentifierHint");
const analysisEmptyState = document.getElementById("analysisEmptyState");
const maximumPdbBytes = 25 * 1024 * 1024;
const resultStorageKey = "plexcalLatestAnalysis";
let molstarScriptPromise;
let viewerInstance;
let activeUploadUrl;
let viewerLoadSequence = 0;
let activeStructureKey;
let metadataLoadSequence = 0;

function setStatus(message, isError = false) {
    statusMessage.textContent = message;
    statusMessage.classList.toggle("is-error", isError);
}

function updateInputSource() {
    const isUpload = sourceInput.value === "upload";
    uploadFields.classList.toggle("hidden", !isUpload);
    identifierFields.classList.toggle("hidden", isUpload);
    fileInput.required = isUpload;
    fileInput.disabled = !isUpload;
    identifierInput.required = !isUpload;
    identifierInput.disabled = isUpload;
    previewButton.disabled = isUpload;
    identifierLabel.textContent =
        sourceInput.value === "rcsb" ? "RCSB PDB ID" : "AlphaFold structure ID";
    identifierInput.placeholder = sourceInput.value === "rcsb"
        ? "e.g. 1TUB"
        : "e.g. P05067 or AF-P05067-F1";
    metadataHint.textContent = sourceInput.value === "rcsb"
        ? "Enter a four-character PDB accession."
        : "Enter a UniProt accession or AlphaFold model ID.";
}

function finiteValue(value) {
    const number = Number(value);
    return value !== null && value !== "" && Number.isFinite(number)
        ? number
        : null;
}

function formatMetric(value, digits = 3) {
    const number = finiteValue(value);
    return number === null ? "—" : number.toLocaleString(undefined, {
        minimumFractionDigits: digits,
        maximumFractionDigits: digits,
    });
}

function metadataField(list, label, value) {
    if (value === null || value === undefined || value === "") {
        return;
    }
    const row = document.createElement("div");
    const term = document.createElement("dt");
    term.textContent = label;
    const description = document.createElement("dd");
    description.textContent = String(value);
    row.append(term, description);
    list.append(row);
}

function metadataLink(label, url) {
    if (!url) {
        return;
    }
    const link = document.createElement("a");
    link.href = url;
    link.target = "_blank";
    link.rel = "noopener noreferrer";
    link.textContent = label;
    metadataLinks.append(link);
}

function renderStructureMetadata(metadata) {
    metadataList.replaceChildren();
    metadataLinks.replaceChildren();
    metadataSourceBadge.textContent = metadata.source || "Structure";

    metadataField(metadataList, "Structure / protein", metadata.title);
    metadataField(metadataList, "Organism", metadata.organism);
    metadataField(metadataList, "Experimental method / model", metadata.method);
    metadataField(
        metadataList,
        "Resolution",
        finiteValue(metadata.resolution_angstrom) === null
            ? null
            : `${formatMetric(metadata.resolution_angstrom, 2)} Å`,
    );
    metadataField(metadataList, "Polymer entities", metadata.polymer_entity_count);
    metadataField(
        metadataList,
        "Deposited atoms",
        metadata.atom_count === null || metadata.atom_count === undefined
            ? null
            : formatMetric(metadata.atom_count, 0),
    );
    metadataField(
        metadataList,
        "Mean pLDDT",
        finiteValue(metadata.mean_plddt) === null
            ? null
            : `${formatMetric(metadata.mean_plddt, 1)} / 100`,
    );
    metadataField(
        metadataList,
        "Model updated",
        metadata.model_updated
            ? new Date(metadata.model_updated).toLocaleDateString()
            : null,
    );
    metadataField(
        metadataList,
        "Sequence version",
        metadata.sequence_version_date
            ? new Date(metadata.sequence_version_date).toLocaleDateString()
            : null,
    );
    metadataField(metadataList, "Publication", metadata.publication_title);
    metadataField(
        metadataList,
        "Authors",
        Array.isArray(metadata.authors) ? metadata.authors.join(", ") : metadata.authors,
    );
    metadataField(
        metadataList,
        "Journal",
        [metadata.journal, metadata.publication_year].filter(Boolean).join(" · "),
    );
    if (metadata.doi) {
        metadataField(metadataList, "DOI", metadata.doi);
    }

    for (const [label, fraction] of [
        ["Very low confidence", metadata.plddt_very_low_fraction],
        ["Low confidence", metadata.plddt_low_fraction],
        ["Confident", metadata.plddt_confident_fraction],
        ["Very high confidence", metadata.plddt_very_high_fraction],
    ]) {
        const number = finiteValue(fraction);
        if (number !== null) {
            metadataField(metadataList, label, `${formatMetric(number * 100, 1)}%`);
        }
    }

    metadataLink("Open structure record", metadata.record_url);
    metadataLink("Read primary publication", metadata.publication_url);
    metadataLink("PubMed record", metadata.pubmed_url);
    metadataLink("UniProt protein record", metadata.uniprot_url);
    metadataStatus.textContent = metadata.publication_title
        ? "Database annotations and the associated publication are shown below."
        : "Structure metadata loaded.";
}

function parseUploadedPdbMetadata(pdbText, filename) {
    const values = {
        source: "Uploaded PDB header",
        identifier: filename,
        record_url: null,
        publication_url: null,
    };
    const titles = [];
    const citationAuthors = [];
    const journalTitles = [];
    const journalReferences = [];
    const organisms = [];
    let method;
    let resolution;
    let doi;
    let pubmed;

    for (const line of pdbText.split(/\r?\n/)) {
        const record = line.slice(0, 6).trim();
        const content = line.slice(10).trim();
        if (record === "TITLE") {
            titles.push(content);
        } else if (record === "SOURCE" && /ORGANISM_SCIENTIFIC:/i.test(line)) {
            const name = line.split(/ORGANISM_SCIENTIFIC:/i).slice(1).join(":")
                .trim().replace(/;$/, "");
            if (name && !organisms.includes(name)) {
                organisms.push(name);
            }
        } else if (record === "EXPDTA") {
            method = content;
        } else if (record === "REMARK" && /RESOLUTION\./i.test(line)) {
            resolution = Number(line.match(/RESOLUTION\.\s+([0-9.]+)/i)?.[1]);
        } else if (record === "JRNL" && /^\s*AUTH\b/i.test(content)) {
            citationAuthors.push(content.replace(/^\s*AUTH\s+/i, ""));
        } else if (record === "JRNL" && /\bTITL\b/i.test(line)) {
            journalTitles.push(content.replace(/^\s*TITL\s*(?:\d)?\s+/i, ""));
        } else if (record === "JRNL" && /\bREF\b/i.test(line)) {
            journalReferences.push(content.replace(/^\s*REF\s*(?:\d)?\s+/i, ""));
        } else if (record === "JRNL" && /\bDOI\b/i.test(line)) {
            doi = content.replace(/^\s*DOI\s+/i, "").trim();
        } else if (record === "JRNL" && /\bPMID\b/i.test(line)) {
            pubmed = content.match(/\d+/)?.[0];
        }
    }

    values.title = titles.join(" ").replace(/\s+/g, " ") || filename;
    values.organism = organisms.join(", ") || null;
    values.method = method;
    values.resolution_angstrom = Number.isFinite(resolution) ? resolution : null;
    values.authors = citationAuthors.join(" ").split(/\s*,\s*/).filter(Boolean);
    values.publication_title = journalTitles.join(" ").replace(/\s+/g, " ");
    values.journal = journalReferences.join(" ").replace(/\s+/g, " ");
    values.doi = doi;
    values.publication_url = doi ? `https://doi.org/${encodeURI(doi)}` : null;
    values.pubmed_url = pubmed
        ? `https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(pubmed)}/`
        : null;
    return values;
}

async function loadStructureMetadata(source, identifier, pdbText, filename) {
    const sequence = ++metadataLoadSequence;
    metadataStatus.textContent = "Retrieving structure annotations and publication…";
    metadataSourceBadge.textContent = source === "rcsb" ? "RCSB PDB" : "AlphaFold DB";
    if (source === "upload") {
        renderStructureMetadata(parseUploadedPdbMetadata(pdbText, filename));
        return;
    }

    const query = new URLSearchParams({ source, id: identifier });
    try {
        const response = await fetch(`/api/structure-metadata?${query}`, {
            cache: "no-store",
        });
        const metadata = await response.json();
        if (!response.ok) {
            throw new Error(metadata.error || `Metadata request returned HTTP ${response.status}.`);
        }
        if (sequence === metadataLoadSequence) {
            renderStructureMetadata(metadata);
        }
    } catch (error) {
        if (sequence !== metadataLoadSequence) {
            return;
        }
        metadataList.replaceChildren();
        metadataLinks.replaceChildren();
        metadataStatus.textContent =
            `Structure loaded; metadata could not be retrieved: ${error.message || error}`;
        console.error("Structure metadata error:", error);
    }
}

function addMetric(list, name, value) {
    const row = document.createElement("div");
    const label = document.createElement("dt");
    label.textContent = name;
    const result = document.createElement("dd");
    result.textContent = value;
    row.append(label, result);
    list.append(row);
}

function displayModelMetrics(method, metrics) {
    const list = document.getElementById(`${method}Metrics`);
    list.replaceChildren();
    const fraction = finiteValue(metrics.coherent_pair_fraction);
    addMetric(
        list,
        "Diffusion length",
        metrics.diffusion_length_angstrom === null
            ? "—"
            : `${formatMetric(metrics.diffusion_length_angstrom)} Å`,
    );
    addMetric(list, "Total pairs", formatMetric(metrics.total_pairs, 0));
    addMetric(list, "Coherent pairs", formatMetric(metrics.coherent_pairs, 0));
    addMetric(
        list,
        "Coherent pair fraction",
        fraction === null ? "—" : `${formatMetric(fraction * 100, 2)}%`,
    );
}

function createBar(value, maximum, colorClass, valueLabel) {
    const column = document.createElement("div");
    column.className = "chart-column";
    const barArea = document.createElement("div");
    barArea.className = "chart-bar-pair";
    const unit = document.createElement("div");
    unit.className = "chart-bar-unit";
    const valueElement = document.createElement("span");
    valueElement.className = "chart-value";
    valueElement.textContent = value === null ? "—" : valueLabel(value);
    const bar = document.createElement("div");
    bar.className = `chart-bar ${colorClass}`;
    const barHeight =
        value === null || maximum <= 0 ? 0 : (value / maximum) * 140;
    bar.style.height = `${barHeight}px`;
    if (value !== null && value > 0) {
        bar.style.minHeight = "2px";
    }
    bar.setAttribute("aria-hidden", "true");
    unit.append(valueElement, bar);
    barArea.append(unit);
    const modelLabel = document.createElement("span");
    modelLabel.className = "chart-label";
    column.append(barArea, modelLabel);
    return { column, modelLabel };
}

function addLegend(chart, entries) {
    chart.nextElementSibling?.remove();
    const legend = document.createElement("div");
    legend.className = "chart-legend";
    for (const { label, className } of entries) {
        const item = document.createElement("span");
        const swatch = document.createElement("i");
        if (className) {
            swatch.classList.add(className);
        }
        item.append(swatch, document.createTextNode(label));
        legend.append(item);
    }
    chart.after(legend);
}

function renderComparisonChart(chart, fret, det, field, valueFormatter) {
    chart.replaceChildren();
    const values = [finiteValue(fret[field]), finiteValue(det[field])];
    const maximum = Math.max(0, ...values.filter((value) => value !== null));
    for (const [index, label] of ["FRET", "DET"].entries()) {
        const bar = createBar(
            values[index],
            maximum,
            index === 0 ? "fret" : "det",
            valueFormatter,
        );
        bar.modelLabel.textContent = label;
        chart.append(bar.column);
    }
}

function renderPairCountChart(chart, fret, det) {
    chart.replaceChildren();
    const values = [
        finiteValue(fret.total_pairs),
        finiteValue(fret.coherent_pairs),
        finiteValue(det.total_pairs),
        finiteValue(det.coherent_pairs),
    ];
    const maximum = Math.max(0, ...values.filter((value) => value !== null));
    const groups = [
        { name: "FRET", total: values[0], coherent: values[1], color: "fret" },
        { name: "DET", total: values[2], coherent: values[3], color: "det" },
    ];
    for (const group of groups) {
        const column = document.createElement("div");
        column.className = "chart-column";
        const bars = document.createElement("div");
        bars.className = "chart-bar-pair";
        for (const [value, kind] of [
            [group.total, "total"],
            [group.coherent, "coherent"],
        ]) {
            const unit = document.createElement("div");
            unit.className = "chart-bar-unit";
            const valueLabel = document.createElement("span");
            valueLabel.className = "chart-value";
            valueLabel.textContent = formatMetric(value, 0);
            const bar = document.createElement("div");
            bar.className = `chart-bar ${group.color} ${kind}`;
            const barHeight =
                value === null || maximum <= 0 ? 0 : (value / maximum) * 140;
            bar.style.height = `${barHeight}px`;
            if (value !== null && value > 0) {
                bar.style.minHeight = "2px";
            }
            bar.setAttribute(
                "aria-label",
                `${group.name} ${kind}: ${formatMetric(value, 0)}`,
            );
            const kindLabel = document.createElement("span");
            kindLabel.className = "chart-kind-label";
            kindLabel.textContent = kind === "total" ? "Total" : "Coherent";
            unit.append(valueLabel, bar, kindLabel);
            bars.append(unit);
        }
        const modelLabel = document.createElement("span");
        modelLabel.className = "chart-label";
        modelLabel.textContent = group.name;
        column.append(bars, modelLabel);
        chart.append(column);
    }
    addLegend(chart, [
        { label: "FRET total", className: "" },
        { label: "FRET coherent", className: "coherent-swatch" },
        { label: "DET total", className: "det-swatch" },
        { label: "DET coherent", className: "det-coherent-swatch" },
    ]);
}

function renderCharts(fret, det) {
    const diffusion = document.getElementById("diffusionLengthChart");
    renderComparisonChart(
        diffusion,
        fret,
        det,
        "diffusion_length_angstrom",
        (value) => `${formatMetric(value)} Å`,
    );
    addLegend(diffusion, [
        { label: "FRET", className: "" },
        { label: "DET", className: "det-swatch" },
    ]);

    renderPairCountChart(document.getElementById("pairCountsChart"), fret, det);

    const fractionChart = document.getElementById("coherentPairFractionChart");
    const fractionFret = finiteValue(fret.coherent_pair_fraction);
    const fractionDet = finiteValue(det.coherent_pair_fraction);
    renderComparisonChart(
        fractionChart,
        { fraction: fractionFret === null ? null : fractionFret * 100 },
        { fraction: fractionDet === null ? null : fractionDet * 100 },
        "fraction",
        (value) => `${formatMetric(value, 2)}%`,
    );
    addLegend(fractionChart, [
        { label: "FRET", className: "" },
        { label: "DET", className: "det-swatch" },
    ]);
}

function displayAnalysis(payload, persist = true) {
    displayModelMetrics("fret", payload.fret);
    displayModelMetrics("det", payload.det);
    document.getElementById("analysisInputSummary").textContent =
        `${payload.input.filename} · ${payload.input.source.toUpperCase()}`;
    renderCharts(payload.fret, payload.det);

    const errors = [
        payload.fret.error && `FRET: ${payload.fret.error}`,
        payload.det.error && `DET: ${payload.det.error}`,
    ].filter(Boolean);
    resultError.textContent = errors.join(" ");
    resultError.classList.toggle("hidden", errors.length === 0);
    resultSection.classList.remove("hidden");
    analysisEmptyState.classList.add("hidden");

    if (persist) {
        try {
            sessionStorage.setItem(resultStorageKey, JSON.stringify(payload));
        } catch (error) {
            console.warn("Could not preserve analysis results for refresh.", error);
        }
    }
}

function restoreAnalysis() {
    try {
        const stored = sessionStorage.getItem(resultStorageKey);
        if (!stored) {
            return;
        }
        const payload = JSON.parse(stored);
        if (
            payload?.input
            && payload?.fret
            && payload?.det
            && typeof payload.input.filename === "string"
            && typeof payload.input.source === "string"
        ) {
            displayAnalysis(payload, false);
            if (payload.input.source !== "upload") {
                void openStructureViewer(
                    payload.input.source,
                    payload.input.identifier,
                );
            }
        } else {
            throw new TypeError("Stored analysis data has an invalid shape.");
        }
    } catch (error) {
        try {
            sessionStorage.removeItem(resultStorageKey);
        } catch (storageError) {
            console.warn("Could not remove invalid stored analysis data.", storageError);
        }
        console.warn("Could not restore previous analysis results.", error);
    }
}

function loadMolstar() {
    if (window.PDBeMolstarPlugin) {
        return Promise.resolve();
    }
    if (!molstarScriptPromise) {
        molstarScriptPromise = (async () => {
            const sources = [
                "https://cdn.jsdelivr.net/npm/pdbe-molstar@3.12.0/build/pdbe-molstar-plugin.js",
                "https://unpkg.com/pdbe-molstar@3.12.0/build/pdbe-molstar-plugin.js",
            ];
            for (const source of sources) {
                try {
                    await new Promise((resolve, reject) => {
                        const script = document.createElement("script");
                        script.src = source;
                        script.onload = resolve;
                        script.onerror = () => {
                            script.remove();
                            reject(new Error(`Could not load PDBe Mol* from ${source}.`));
                        };
                        document.head.append(script);
                    });
                    if (window.PDBeMolstarPlugin) {
                        return;
                    }
                    throw new Error(`The PDBe Mol* bundle at ${source} did not initialize.`);
                } catch (error) {
                    if (source === sources[sources.length - 1]) {
                        throw new Error(
                            `The interactive Mol* viewer could not be loaded. ${error.message}`,
                        );
                    }
                }
            }
        })();
    }
    return molstarScriptPromise;
}

function structureViewerUrl(source, identifier) {
    if (source === "upload") {
        return activeUploadUrl;
    }
    const query = new URLSearchParams({ source, id: identifier });
    return `/api/structure?${query}`;
}

async function openStructureViewer(source, identifier, uploadedFile = null) {
    const sequence = ++viewerLoadSequence;
    const structureKey = `${source}:${source === "upload" ? activeUploadUrl : identifier}`;
    viewerStatus.textContent = "Loading structure viewer…";
    try {
        const pdbUrl = structureViewerUrl(source, identifier);
        if (!pdbUrl) {
            throw new Error("The selected PDB file is not available to the viewer.");
        }
        const pdbText = uploadedFile
            ? await uploadedFile.text()
            : await fetchStructureText(pdbUrl);
        if (sequence !== viewerLoadSequence) {
            return;
        }
        void loadStructureMetadata(
            source,
            identifier,
            pdbText,
            uploadedFile?.name || identifier,
        );
        const isAlphaFold = source === "alphafold";
        const options = {
            customData: { url: pdbUrl, format: "pdb", binary: false },
            alphafoldView: isAlphaFold,
            sequencePanel: true,
            hideControls: false,
            loadingOverlay: true,
            subscribeEvents: false,
            bgColor: { r: 255, g: 255, b: 255 },
        };
        await loadMolstar();
        if (sequence !== viewerLoadSequence) {
            return;
        }
        if (viewerInstance) {
            await viewerInstance.visual.update(options, true);
        } else {
            viewerInstance = new window.PDBeMolstarPlugin();
            await viewerInstance.render(viewer, options);
        }
        if (sequence !== viewerLoadSequence) {
            return;
        }
        activeStructureKey = structureKey;
        visibilityControls.disabled = false;
        document.getElementById("viewerPlaceholder").classList.add("hidden");
        viewerStatus.textContent = isAlphaFold
            ? "AlphaFold model loaded. Full Mol* sequence, representation, selection, and display controls are available."
            : "Structure loaded. Full Mol* sequence, representation, selection, and display controls are available.";
        await applyInitialStructureVisibility();
    } catch (error) {
        if (sequence === viewerLoadSequence) {
            viewerStatus.textContent = error.message || "Could not load this structure.";
            document.getElementById("viewerPlaceholder").classList.remove("hidden");
            console.error("3D structure viewer error:", error);
        }
    }
}

async function applyInitialStructureVisibility() {
    for (const checkbox of visibilityControls.querySelectorAll(
        "[data-structure-visibility]",
    )) {
        try {
            await viewerInstance.visual.visibility({
                [checkbox.dataset.structureVisibility]: checkbox.checked,
            });
        } catch (error) {
            viewerStatus.textContent =
                `Could not set ${checkbox.parentElement.textContent.trim()} visibility: ${error.message || error}`;
            console.error("Mol* initial visibility error:", error);
            return;
        }
    }
}

async function fetchStructureText(url) {
    const response = await fetch(url, { cache: "no-store" });
    if (!response.ok) {
        const payload = await response.json().catch(() => ({}));
        throw new Error(
            payload.error || `Structure download returned HTTP ${response.status}.`,
        );
    }
    return response.text();
}

function validateSelectedUpload(file) {
    if (!file) {
        throw new Error("Choose a PDB file first.");
    }
    if (file.size === 0) {
        throw new Error("The selected file is empty.");
    }
    if (file.size > maximumPdbBytes) {
        throw new Error("The PDB file must be 25 MiB or smaller.");
    }
}

async function previewSelectedStructure() {
    const source = sourceInput.value;
    const identifier = identifierInput.value.trim();
    if (!identifier) {
        setStatus("Enter a structure ID to preview.", true);
        identifierInput.focus();
        return;
    }
    viewerStatus.textContent = `Loading ${identifier} from ${source === "rcsb" ? "RCSB PDB" : "AlphaFold DB"}…`;
    await openStructureViewer(source, identifier);
}

async function calculateStructure(event) {
    event.preventDefault();
    resultError.classList.add("hidden");
    setStatus("Calculating FRET and DET locally. Larger structures may take a while.");

    const source = sourceInput.value;
    const file = fileInput.files[0];
    const identifier = identifierInput.value.trim();
    try {
        if (source === "upload") {
            validateSelectedUpload(file);
        } else if (!identifier) {
            throw new Error("Enter a structure ID to calculate.");
        }

        calculateButton.disabled = true;
        calculateButton.textContent = "Calculating…";
        const formData = new FormData();
        formData.append("inputSource", source);
        if (source === "upload") {
            formData.append("pdbFile", file, file.name);
        } else {
            formData.append("structureId", identifier);
        }
        const response = await fetch("/api/calculate", {
            method: "POST",
            body: formData,
        });
        const payload = await response.json();
        if (!response.ok || !payload.fret || !payload.det || !payload.input) {
            throw new Error(
                payload.error || `Calculation failed with HTTP ${response.status}.`,
            );
        }

        if (source === "upload") {
            if (!activeUploadUrl) {
                activeUploadUrl = URL.createObjectURL(file);
            }
        }
        displayAnalysis(payload);
        setStatus(`FRET and DET diffusion-length analysis complete for ${payload.input.filename}.`);
        const structureKey =
            `${source}:${source === "upload" ? activeUploadUrl : identifier}`;
        if (activeStructureKey !== structureKey) {
            void openStructureViewer(source, identifier);
        }
    } catch (error) {
        setStatus(error.message || "Structure calculation failed.", true);
        console.error("Structure calculation error:", error);
    } finally {
        calculateButton.disabled = false;
        calculateButton.textContent = "Calculate FRET and DET";
    }
}

async function checkFretCalculator() {
    try {
        const response = await fetch("/api/status");
        if (!response.ok) {
            throw new Error(`Local server returned HTTP ${response.status}.`);
        }
        const status = await response.json();
        if (
            status.fretCalculationAvailable !== true
            || status.detCalculationAvailable !== true
        ) {
            throw new Error("FRET and DET calculations are not available on this server.");
        }
        calculateButton.disabled = false;
        setStatus("Local FRET and DET calculators ready.");
    } catch (error) {
        setStatus(
            "Structure calculations require the local app. Run `python server.py` and open http://127.0.0.1:8000.",
            true,
        );
        console.info("Local structure calculator unavailable:", error);
    }
}

expandButton.addEventListener("click", () => {
    viewerInstance?.canvas.toggleExpanded();
});

visibilityControls.addEventListener("change", async (event) => {
    const target = event.target;
    if (!(target instanceof HTMLInputElement) || !target.dataset.structureVisibility) {
        return;
    }
    try {
        await viewerInstance?.visual.visibility({
            [target.dataset.structureVisibility]: target.checked,
        });
    } catch (error) {
        viewerStatus.textContent =
            `Could not update structure visibility: ${error.message || error}`;
        console.error("Mol* visibility update error:", error);
    }
});

sourceInput.addEventListener("change", updateInputSource);
fileInput.addEventListener("change", () => {
    const file = fileInput.files[0];
    if (!file) {
        return;
    }
    try {
        validateSelectedUpload(file);
        if (activeUploadUrl) {
            URL.revokeObjectURL(activeUploadUrl);
        }
        activeUploadUrl = URL.createObjectURL(file);
        void openStructureViewer("upload", file.name, file);
    } catch (error) {
        viewerStatus.textContent = error.message;
        setStatus(error.message, true);
    }
});
previewButton.addEventListener("click", previewSelectedStructure);
form.addEventListener("submit", calculateStructure);

updateInputSource();
restoreAnalysis();
void checkFretCalculator();
