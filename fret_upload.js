const fretUploadForm = document.getElementById("fretUploadForm");
const fretUploadButton = document.getElementById("fretUploadButton");
const fretUploadStatus = document.getElementById("fretUploadStatus");
const fretUploadResult = document.getElementById("fretUploadResult");
const fretUploadResultError = document.getElementById("fretUploadResultError");
const pdbFileInput = document.getElementById("pdbFile");
const inputSource = document.getElementById("inputSource");
const structureIdInput = document.getElementById("structureId");
const uploadSourceFields = document.getElementById("uploadSourceFields");
const idSourceFields = document.getElementById("idSourceFields");
const calculationButtonLabel = fretUploadButton.querySelector("span:first-child");
const maximumPdbBytes = 25 * 1024 * 1024;

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
            "Structure calculations require the local app. Run `python server.py` and open http://127.0.0.1:8000.";
        fretUploadStatus.classList.add("is-error");
        console.info("Local structure calculator unavailable:", error);
    }
}

function updateInputSource() {
    const isUpload = inputSource.value === "upload";
    uploadSourceFields.classList.toggle("hidden", !isUpload);
    idSourceFields.classList.toggle("hidden", isUpload);
    pdbFileInput.required = isUpload;
    pdbFileInput.disabled = !isUpload;
    structureIdInput.required = !isUpload;
    structureIdInput.disabled = isUpload;
    structureIdInput.placeholder = inputSource.value === "rcsb"
        ? "e.g. 1TUB"
        : "e.g. P05067 or AF-P05067-F1";
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

function displayModelMetrics(modelName, metrics) {
    const prefix = modelName === "fret" ? "fretUpload" : "detUpload";
    document.getElementById(`${prefix}Length`).textContent =
        metrics.diffusion_length_angstrom === null
            ? "—"
            : `${formatMetric(metrics.diffusion_length_angstrom)} Å`;
    document.getElementById(`${prefix}TotalPairs`).textContent =
        formatMetric(metrics.total_pairs, 0);
    document.getElementById(`${prefix}CoherentPairs`).textContent =
        formatMetric(metrics.coherent_pairs, 0);
    const fraction = finiteValue(metrics.coherent_pair_fraction);
    document.getElementById(`${prefix}Fraction`).textContent =
        fraction === null ? "—" : `${formatMetric(fraction * 100, 2)}%`;
}

function makeBar(value, maximum, className, valueFormatter) {
    const column = document.createElement("div");
    column.className = "chart-column";
    const shownValue = document.createElement("span");
    shownValue.className = "chart-value";
    shownValue.textContent = value === null ? "—" : valueFormatter(value);
    const pair = document.createElement("div");
    pair.className = "chart-bar-pair";
    const bar = document.createElement("div");
    bar.className = `chart-bar ${className}`;
    bar.style.height = `${value === null || maximum <= 0 ? 0 : (value / maximum) * 140}px`;
    bar.setAttribute("aria-hidden", "true");
    pair.append(bar);
    const label = document.createElement("span");
    label.className = "chart-label";
    column.append(shownValue, pair, label);
    return { column, label };
}

function addLegend(container, entries) {
    const legend = document.createElement("div");
    legend.className = "chart-legend";
    for (const entry of entries) {
        const item = document.createElement("span");
        const swatch = document.createElement("i");
        if (entry.className) {
            swatch.classList.add(entry.className);
        }
        item.append(swatch, document.createTextNode(entry.label));
        legend.append(item);
    }
    container.after(legend);
}

function renderModelChart(host, fret, det, field, valueFormatter) {
    host.replaceChildren();
    host.nextElementSibling?.remove();
    const values = [finiteValue(fret[field]), finiteValue(det[field])];
    const available = values.filter((value) => value !== null);
    const maximum = Math.max(0, ...available);
    for (const [index, model] of ["FRET", "DET"].entries()) {
        const bar = makeBar(values[index], maximum, index === 0 ? "fret" : "det", valueFormatter);
        bar.label.textContent = model;
        host.append(bar.column);
    }
}

function renderPairCountChart(host, fret, det) {
    host.replaceChildren();
    host.nextElementSibling?.remove();
    const values = [
        finiteValue(fret.total_pairs),
        finiteValue(fret.coherent_pairs),
        finiteValue(det.total_pairs),
        finiteValue(det.coherent_pairs),
    ];
    const maximum = Math.max(0, ...values.filter((value) => value !== null));
    const groups = [
        { model: "FRET", total: values[0], coherent: values[1], className: "fret" },
        { model: "DET", total: values[2], coherent: values[3], className: "det" },
    ];
    for (const group of groups) {
        const column = document.createElement("div");
        column.className = "chart-column";
        const bars = document.createElement("div");
        bars.className = "chart-bar-pair";
        for (const [value, type] of [
            [group.total, "total"],
            [group.coherent, "coherent"],
        ]) {
            const bar = document.createElement("div");
            bar.className = `chart-bar ${group.className} ${type}`;
            bar.style.height =
                `${value === null || maximum <= 0 ? 0 : (value / maximum) * 140}px`;
            bar.setAttribute("aria-label", `${group.model} ${type}: ${formatMetric(value, 0)}`);
            bars.append(bar);
        }
        const valueLabel = document.createElement("span");
        valueLabel.className = "chart-value";
        valueLabel.textContent = `Total ${formatMetric(group.total, 0)} · Coherent ${formatMetric(group.coherent, 0)}`;
        const label = document.createElement("span");
        label.className = "chart-label";
        label.textContent = group.model;
        column.append(valueLabel, bars, label);
        host.append(column);
    }
    addLegend(host, [
        { label: "FRET total", className: "" },
        { label: "FRET coherent", className: "coherent-swatch" },
        { label: "DET total", className: "det-swatch" },
        { label: "DET coherent", className: "det-coherent-swatch" },
    ]);
}

function renderCharts(fret, det) {
    const diffusionChart = document.getElementById("diffusionChart");
    renderModelChart(
        diffusionChart,
        fret,
        det,
        "diffusion_length_angstrom",
        (value) => `${formatMetric(value)} Å`,
    );
    addLegend(diffusionChart, [
        { label: "FRET", className: "" },
        { label: "DET", className: "det-swatch" },
    ]);
    renderPairCountChart(
        document.getElementById("pairCountsChart"),
        fret,
        det,
    );
    const fractionChart = document.getElementById("fractionChart");
    const fretFraction = finiteValue(fret.coherent_pair_fraction);
    const detFraction = finiteValue(det.coherent_pair_fraction);
    const percentageFret = fretFraction === null ? null : fretFraction * 100;
    const percentageDet = detFraction === null ? null : detFraction * 100;
    renderModelChart(
        fractionChart,
        { coherent_pair_fraction: percentageFret },
        { coherent_pair_fraction: percentageDet },
        "coherent_pair_fraction",
        (value) => `${formatMetric(value, 2)}%`,
    );
    addLegend(fractionChart, [
        { label: "FRET", className: "" },
        { label: "DET", className: "det-swatch" },
    ]);
}

function displayAnalysis(payload) {
    const { fret, det } = payload;
    displayModelMetrics("fret", fret);
    displayModelMetrics("det", det);
    document.getElementById("analysisTitle").textContent =
        `${payload.input.identifier} · ${payload.input.filename}`;
    renderCharts(fret, det);
    const errors = [fret.error && `FRET: ${fret.error}`, det.error && `DET: ${det.error}`]
        .filter(Boolean);
    if (errors.length) {
        fretUploadResultError.textContent = errors.join(" ");
        fretUploadResultError.classList.remove("hidden");
    } else {
        fretUploadResultError.classList.add("hidden");
    }
    fretUploadResult.classList.remove("hidden");
}

async function calculateStructure(event) {
    event.preventDefault();
    fretUploadResult.classList.add("hidden");
    fretUploadResultError.classList.add("hidden");
    fretUploadStatus.classList.remove("is-error");

    const source = inputSource.value;
    const file = pdbFileInput.files[0];
    const identifier = structureIdInput.value.trim();
    if (source === "upload") {
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
    } else if (!identifier) {
        fretUploadStatus.textContent = "Enter a structure ID to calculate.";
        fretUploadStatus.classList.add("is-error");
        structureIdInput.focus();
        return;
    }

    fretUploadButton.disabled = true;
    calculationButtonLabel.textContent = "Calculating…";
    fretUploadStatus.textContent =
        "Calculating FRET and DET locally. Larger structures may take a while.";

    try {
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
        if (!payload.fret || !payload.det || !payload.input) {
            throw new Error(payload.error || `Calculation failed with HTTP ${response.status}`);
        }
        displayAnalysis(payload);
        if (!response.ok) {
            fretUploadStatus.textContent =
                `Calculation finished with an error for ${payload.input.filename}.`;
            fretUploadStatus.classList.add("is-error");
        } else {
            fretUploadStatus.textContent =
                `FRET and DET calculations complete for ${payload.input.filename}.`;
        }
    } catch (error) {
        fretUploadStatus.textContent = error.message || "Structure calculation failed.";
        fretUploadStatus.classList.add("is-error");
        console.error("Structure calculation error:", error);
    } finally {
        fretUploadButton.disabled = false;
        calculationButtonLabel.textContent = "Calculate FRET + DET";
    }
}

inputSource.addEventListener("change", updateInputSource);
fretUploadForm.addEventListener("submit", calculateStructure);
updateInputSource();
checkFretCalculator();
