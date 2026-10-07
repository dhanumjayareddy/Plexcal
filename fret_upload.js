const fretUploadForm = document.getElementById("fretUploadForm");
const fretUploadButton = document.getElementById("fretUploadButton");
const fretUploadStatus = document.getElementById("fretUploadStatus");
const fretUploadResult = document.getElementById("fretUploadResult");
const fretUploadValue = document.getElementById("fretUploadValue");
const detUploadValue = document.getElementById("detUploadValue");
const fretResultRow = document.getElementById("fretResultRow");
const detResultRow = document.getElementById("detResultRow");
const fretUploadResultError = document.getElementById("fretUploadResultError");
const pdbFileInput = document.getElementById("pdbFile");
const calculationOptions = document.querySelectorAll('input[name="calculation"]');
const calculationButtonLabel = fretUploadButton.querySelector("span:first-child");
const maximumPdbBytes = 25 * 1024 * 1024;

function selectedCalculation() {
    return document.querySelector('input[name="calculation"]:checked').value;
}

function updateCalculationButton() {
    const labels = {
        fret: "Calculate FRET",
        det: "Calculate DET",
        both: "Calculate both",
    };
    calculationButtonLabel.textContent = labels[selectedCalculation()];
}

async function checkFretCalculator() {
    try {
        const response = await fetch("/api/status");
        if (!response.ok) {
            throw new Error(`Local server returned HTTP ${response.status}`);
        }

        const status = await response.json();
        if (
            status.fretCalculationAvailable !== true
            || status.detCalculationAvailable !== true
        ) {
            throw new Error("FRET and DET calculations are not available on this server.");
        }

        fretUploadButton.disabled = false;
        fretUploadStatus.textContent = "Local FRET and DET calculators ready.";
    } catch (error) {
        fretUploadStatus.textContent =
            "PDB upload calculations require the local app. Run `python server.py` and open http://127.0.0.1:8000.";
        fretUploadStatus.classList.add("is-error");
        console.info("Local FRET calculator unavailable:", error);
    }
}

async function calculateUploadedStructure(event) {
    event.preventDefault();
    const file = pdbFileInput.files[0];
    fretUploadResult.classList.add("hidden");
    fretUploadResultError.classList.add("hidden");
    fretUploadStatus.classList.remove("is-error");

    if (!file) {
        fretUploadStatus.textContent = "Choose a PDB file to calculate.";
        fretUploadStatus.classList.add("is-error");
        return;
    }
    if (file.size === 0) {
        fretUploadStatus.textContent = "The selected file is empty.";
        fretUploadStatus.classList.add("is-error");
        return;
    }
    if (file.size > maximumPdbBytes) {
        fretUploadStatus.textContent = "The PDB file must be 25 MiB or smaller.";
        fretUploadStatus.classList.add("is-error");
        return;
    }

    const calculation = selectedCalculation();
    fretUploadButton.disabled = true;
    calculationButtonLabel.textContent = "Calculating…";
    const calculationDescription = calculation === "both"
        ? "FRET and DET"
        : calculation.toUpperCase();
    fretUploadStatus.textContent =
        `Calculating ${calculationDescription} locally. Larger structures may take a while.`;

    try {
        const formData = new FormData();
        formData.append("pdbFile", file, file.name);
        formData.append("calculation", calculation);
        const response = await fetch("/api/calculate", {
            method: "POST",
            body: formData,
        });
        const payload = await response.json();
        const results = payload.calculations;
        if (!results) {
            throw new Error(payload.error || `Calculation failed with HTTP ${response.status}`);
        }

        fretResultRow.classList.toggle("hidden", !results.fret);
        detResultRow.classList.toggle("hidden", !results.det);
        const errors = [];
        if (results.fret) {
            if (results.fret.error) {
                fretUploadValue.textContent = "Unavailable";
                errors.push(`FRET: ${results.fret.error}`);
            } else {
                fretUploadValue.textContent =
                    `${Number(results.fret.diffusionLengthAngstroms).toFixed(3)} Å`;
            }
        }
        if (results.det) {
            if (results.det.error) {
                detUploadValue.textContent = "Unavailable";
                errors.push(`DET: ${results.det.error}`);
            } else {
                detUploadValue.textContent =
                    `${Number(results.det.diffusionLengthAngstroms).toFixed(3)} Å`;
            }
        }
        fretUploadResult.classList.remove("hidden");
        if (errors.length) {
            fretUploadResultError.textContent = errors.join(" ");
            fretUploadResultError.classList.remove("hidden");
            fretUploadStatus.textContent =
                `Calculation finished with an error for ${payload.fileName}.`;
            fretUploadStatus.classList.add("is-error");
        } else {
            fretUploadStatus.textContent =
                `${calculationDescription} calculation complete for ${payload.fileName}.`;
        }
    } catch (error) {
        fretUploadStatus.textContent = error.message || "FRET calculation failed.";
        fretUploadStatus.classList.add("is-error");
        console.error("Uploaded PDB calculation error:", error);
    } finally {
        fretUploadButton.disabled = false;
        updateCalculationButton();
    }
}

calculationOptions.forEach((option) => {
    option.addEventListener("change", updateCalculationButton);
});
fretUploadForm.addEventListener("submit", calculateUploadedStructure);
updateCalculationButton();
checkFretCalculator();
