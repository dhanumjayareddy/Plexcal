const searchForm = document.getElementById("searchForm");
const searchButton = document.getElementById("searchButton");
const searchInput = document.getElementById("pdbInput");
const searchStatus = document.getElementById("searchStatus");
const connectionStatus = document.getElementById("connectionStatus");
const recordCount = document.getElementById("recordCount");

const loadedShards = new Map();
let manifest = null;
let manifestPromise;

async function loadManifest() {
    try {
        const response = await fetch("data/index.json");
        if (!response.ok) {
            throw new Error(`Static index returned HTTP ${response.status}`);
        }

        manifest = await response.json();
        recordCount.textContent = `${manifest.recordCount.toLocaleString()} structures indexed · only matching data is loaded`;
        connectionStatus.classList.add("is-connected");
        connectionStatus.innerHTML = '<span class="status-dot"></span> Static database ready';
    } catch (error) {
        recordCount.textContent = "Static data unavailable";
        connectionStatus.classList.add("is-disconnected");
        connectionStatus.innerHTML = '<span class="status-dot"></span> Static data unavailable';
        searchStatus.textContent = "Could not load the data index. Refresh the page to try again.";
        console.error("Static data index error:", error);
    }
}

function normalizePDBID(value) {
    return value.trim().toUpperCase();
}

async function loadShard(prefix) {
    if (loadedShards.has(prefix)) {
        return loadedShards.get(prefix);
    }

    const response = await fetch(`data/shards/${encodeURIComponent(prefix)}.json`);
    if (!response.ok) {
        throw new Error(`Data shard "${prefix}" returned HTTP ${response.status}`);
    }

    const shard = await response.json();
    loadedShards.set(prefix, shard);
    return shard;
}

async function searchDatabase(event) {
    event.preventDefault();
    const pdbID = normalizePDBID(searchInput.value);
    searchInput.value = pdbID;

    if (!pdbID) {
        showInitialMessage();
        searchStatus.textContent = "Enter a four-character PDB ID to begin.";
        searchInput.focus();
        return;
    }

    if (!/^[A-Z0-9]{4}$/.test(pdbID)) {
        showInitialMessage();
        searchStatus.textContent = "PDB IDs contain exactly four letters or numbers.";
        searchInput.focus();
        return;
    }

    if (!manifest && manifestPromise) {
        searchStatus.textContent = "Preparing the static data index…";
        await manifestPromise;
    }
    if (!manifest) {
        searchStatus.textContent = "The static data index is not available. Refresh the page to try again.";
        return;
    }

    const prefix = pdbID.slice(0, 2).toLowerCase();
    if (!Object.hasOwn(manifest.shards, prefix)) {
        showNoResult(pdbID);
        return;
    }

    setView("initialMessage");
    searchButton.disabled = true;
    searchButton.classList.add("is-loading");
    searchButton.querySelector("span:first-child").textContent = "Searching";
    searchStatus.textContent = `Loading the small data shard for ${pdbID}…`;

    try {
        const shard = await loadShard(prefix);
        const record = shard[pdbID];
        if (!record) {
            showNoResult(pdbID);
            return;
        }
        displayRecord(record, pdbID);
    } catch (error) {
        setView("initialMessage");
        searchStatus.textContent = "Could not retrieve this data shard. Check your connection and try again.";
        console.error("Static PDB lookup error:", error);
    } finally {
        searchButton.disabled = false;
        searchButton.classList.remove("is-loading");
        searchButton.querySelector("span:first-child").textContent = "Find structure";
    }
}

function displayRecord(record, requestedID) {
    setView("resultSection");
    setText("pdbIdDisplay", record.PDB_ID || requestedID);
    setText("fretValue", formatNumber(record.FRET_Diffusion_Length_Angstroms));
    setText("detValue", formatNumber(record.DET_Diffusion_Length_Angstroms));
    setText("species", record.Species_x);
    setText("scientificName", record.Scientific_Name);
    setText("taxId", record.TaxID);
    setText("kingdom", record.Kingdom);
    setText("phylum", record.Phylum);
    setText("className", record.Class);
    setText("order", record.Order);
    setText("family", record.Family);
    setText("genus", record.Genus);
    setText("speciesY", record.Species_y);
    displayLineage(record.Full_Lineage);
    displayRawRecord(record);
    searchStatus.textContent = `Showing database record ${record.PDB_ID || requestedID}`;
}

function formatNumber(value) {
    if (value === null || value === undefined || value === "") {
        return "—";
    }
    const number = Number(value);
    return Number.isFinite(number) ? number.toFixed(3) : String(value);
}

function setText(id, value) {
    const element = document.getElementById(id);
    if (!element) {
        return;
    }
    element.textContent = value === null || value === undefined || value === ""
        ? "—"
        : String(value);
}

function displayLineage(lineage) {
    const parts = typeof lineage === "string"
        ? lineage.split(";").map((item) => item.trim()).filter(Boolean)
        : [];
    document.getElementById("lineage").textContent = parts.length
        ? parts.join("  →  ")
        : "—";
}

function displayRawRecord(record) {
    const container = document.getElementById("rawRecord");
    container.replaceChildren();

    for (const [key, value] of Object.entries(record)) {
        const field = document.createElement("div");
        field.className = "raw-field";

        const fieldName = document.createElement("span");
        fieldName.className = "raw-field-name";
        fieldName.textContent = key.replaceAll("_", " ");

        const fieldValue = document.createElement("span");
        fieldValue.className = "raw-field-value";
        fieldValue.textContent = value === null || value === undefined || value === ""
            ? "—"
            : String(value);

        field.append(fieldName, fieldValue);
        container.append(field);
    }
}

function setView(visibleId) {
    for (const id of ["resultSection", "noResult", "initialMessage"]) {
        document.getElementById(id).classList.toggle("hidden", id !== visibleId);
    }
}

function showNoResult(pdbID) {
    setView("noResult");
    setText("noResultId", pdbID);
    searchStatus.textContent = `No record found for ${pdbID}`;
}

function showInitialMessage() {
    setView("initialMessage");
    searchStatus.textContent = "";
}

searchForm.addEventListener("submit", searchDatabase);
manifestPromise = loadManifest();
