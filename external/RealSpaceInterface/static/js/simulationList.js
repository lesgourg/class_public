/**
 * Number of decimal digits to round values to in the simulation list dialog.
 * @const
 * @type {number}
 */
const SIMU_LIST_ROUNDING = 3;

/**
 * Will be called when a simulation is deleted.
 *
 * @callback SimuTable~deleteCallback
 * @param {number} index - index of simulation to be deleted
 */

/**
 * Will be called when a simulation is activated.
 *
 * @callback SimuTable~activateCallback
 * @param {number} index - index of simulation to be activated
 */

/**
 * Will be called when a simulation is deactivated.
 * @callback SimuTable~deactivateCallback
 * @param {number} index - index of simulation to be deactivated
 */

/**
 * Will be called when a .gif needs to be created for a simulation.
 * @callback SimuTable~gifCallback
 * @param {number} index - index of simulation to create a .gif of
 */

/**
 * Will be called when an image needs to be created for a simulation at the current frame.
 * @callback SimuTable~imageCallback
 * @param {number} index - index of simulation to create an image of
 */

/**
 * Will be called when user requests simulation to be displayed in 'single' mode
 * @callback SimuTable~singleModeCallback
 * @param {number} index - index of simulation to display in single mode
 */

/**
 * Will be called when user requests camera to be focused on a simulation.
 * @callback SimuTable~focusCameraCallback
 * @param {number} index - index of simulation to focus in view
 */

/**
 * Will be called when user requests to load parameters of simulation
 * into control panel.
 * @callback SimuTable~loadParamCallback
 * @param {number} index - index of simulation of which to load parameters
 */

/**
 *
 * Implementation of a list dialog which displays previously run simulations
 * along with their respective set of cosmological parameters.
 * Exposes additional functionality such as hiding and showing individual
 * simulations, entering 'single' mode, in which only a single simulation
 * (but all quantities of that simulation) are shown, creating a .gif and
 * deleting a simulation from memory.
 * The passed callbacks are bound to the appropriate buttons and actions.
 *
 * @constructor
 *
 * @param {SimuTable~deleteCallback} deleteCallback
 * @param {SimuTable~activateCallback} activateCallback
 * @param {SimuTable~deactivateCallback} deactivateCallback
 * @param {SimuTable~gifCallback} gifCallback
 * @param {SimuTable~imageCallback} imageCallback
 * @param {SimuTable~singleModeCallback} singleModeCallback
 * @param {SimuTable~focusCameraCallback} focusCameraCallback
 * @param {SimuTable~loadParamCallback} loadParamCallback
 */
function SimuTable(deleteCallback, activateCallback, deactivateCallback,
    gifCallback, imageCallback, singleModeCallback, focusCameraCallback,
    loadParamCallback) {
    this.deleteCallback = deleteCallback;
    this.activateCallback = activateCallback;
    this.deactivateCallback = deactivateCallback;
    this.gifCallback = gifCallback;
    this.imageCallback = imageCallback;
    this.singleModeCallback = singleModeCallback;
    this.focusCameraCallback = focusCameraCallback;
    this.loadParamCallback = loadParamCallback;

    this.table = document.getElementById("simulationTable");
    this.tableHead = document.getElementById("simulationTableHead");
    this.tableBody = document.getElementById("simulationTableBody");
    this.rowTemplate = document.getElementById("simulationTableRowTemplate");

    this.activeSet = 0;

    this.displayParams = ["omega_b", "omega_m", "Omega_k", "N_ur", "w0_fld", "wa_fld"];
    this.displayLabels = ["&omega;<sub>b</sub>", "&omega;<sub>m</sub>", "&Omega;<sub>k</sub>",
                            "N<sub>ur</sub>", "w<sub>0,fld</sub>", "w<sub>a,fld</sub>"];
}

/**
 * Creates the table header for the modal.
 * Reads table column headings from {@link COSMOLOGICAL_PARAMETER_LIST}.
 *
 */
SimuTable.prototype.createHeader = function() {
    var headRow = document.createElement("tr");

    var activeTh = document.createElement("th");
    activeTh.textContent = "Active";
    headRow.appendChild(activeTh);

    var idxTh = document.createElement("th");
    idxTh.textContent = "#";
    idxTh.setAttribute("scope", "col");

    headRow.appendChild(idxTh);

    for (var entry of COSMOLOGICAL_PARAMETER_LIST) {
        var th = document.createElement("th");
        th.innerHTML = entry.displayName;
        headRow.appendChild(th);
    }

    var actionsTh = document.createElement("th");
    actionsTh.textContent = "Actions";

    var colorTh = document.createElement("th");
    colorTh.textContent = "Color";

    headRow.appendChild(actionsTh);
    headRow.appendChild(colorTh);

    this.tableHead.appendChild(headRow);
}

/**
 * Given a list of simulations and a list of active simulations,
 * populates the dialog with their data.
 *
 * @param {Simulation[]} simulations - list of all simulations
 * @param {number[]} activeList - array of indices of active simulations.
 *                                required to determine whether to show
 *                                individual simulations as active or
 *                                inactive.
 */
SimuTable.prototype.populate = function(simulations, activeList) {
    var self = this;
    // Data
    simulations.forEach(function(sim, i) {
        var isActive = activeList.includes(i);

        var sim = simulations[i];
        var row = document.createElement("tr");
        row.setAttribute("scope", "row");
        row.classList.toggle("table-primary", isActive);

        // Active checkbox
        var checkbox = document.createElement("input");
        checkbox.setAttribute("type", "checkbox");
        checkbox.classList.add("visible-checkbox");

        checkbox.onclick = function() {
            row.classList.toggle("table-primary");
            if (checkbox.checked) {
                self.activateCallback(i);
            } else {
                self.deactivateCallback(i);
            }
        };
        var checkboxTd = document.createElement("td");
        checkbox.checked = isActive;
        checkboxTd.appendChild(checkbox);
        row.appendChild(checkboxTd);

        // Index
        var idxCol = document.createElement("td");
        idxCol.textContent = i;
        row.appendChild(idxCol);

        for (var entry of COSMOLOGICAL_PARAMETER_LIST) {
            var td = document.createElement("td");
            td.textContent = round(sim.params[entry.name], 3);
            row.appendChild(td);
        }

        /* BUTTONS */
        var buttonGroup = document.createElement("div");
        buttonGroup.classList.add("btn-group", "btn-group-sm");

        // Expand in Single Mode Button
        var singleModeButton = document.createElement("button");
        singleModeButton.classList.add("btn", "btn-primary", "simulationTableSingleModeBtn");
        var singleIconSpan = document.createElement("span");
        singleIconSpan.classList.add("oi", "oi-layers");
        singleModeButton.appendChild(singleIconSpan);
        singleModeButton.onclick = function() {
            self.singleModeCallback(i);
        };
        singleModeButton.setAttribute("data-toggle", "tooltip");
        singleModeButton.setAttribute("data-placement", "top");
        singleModeButton.setAttribute("title", "Show all quantities side by side");
        buttonGroup.appendChild(singleModeButton);

        // Image Button
        var imgButton = document.createElement("button");
        imgButton.classList.add("btn", "btn-success", "simulationTableImgBtn");
        var imgIconSpan = document.createElement("span");
        imgIconSpan.classList.add("oi", "oi-image");
        imgButton.appendChild(imgIconSpan);
        imgButton.onclick = function() {
            self.imageCallback(i);
        };
        imgButton.setAttribute("data-toggle", "tooltip");
        imgButton.setAttribute("data-placement", "top");
        imgButton.setAttribute("title", "Take Snapshot");
        buttonGroup.appendChild(imgButton);

        // GIF Button
        var gifButton = document.createElement("button");
        gifButton.classList.add("btn", "btn-success", "simulationTableGifBtn");
        var gifIconSpan = document.createElement("span");
        gifIconSpan.classList.add("oi", "oi-video");
        gifButton.appendChild(gifIconSpan);
        gifButton.onclick = function() {
            self.gifCallback(i);
        };
        gifButton.setAttribute("data-toggle", "tooltip");
        gifButton.setAttribute("data-placement", "top");
        gifButton.setAttribute("title", "Create .gif");
        buttonGroup.appendChild(gifButton);

        // Load Parameters Button
        var loadParamButton = document.createElement("button");
        loadParamButton.classList.add("btn", "btn-secondary",
            "simulationTableLoadParamBtn");
        var paramIconSpan = document.createElement("span");
        paramIconSpan.classList.add("oi", "oi-list-rich");
        loadParamButton.appendChild(paramIconSpan);
        loadParamButton.onclick = function() {
            self.loadParamCallback(i);
        };
        loadParamButton.setAttribute("data-toggle", "tooltip");
        loadParamButton.setAttribute("data-placement", "top");
        loadParamButton.setAttribute("title", "Load Parameters in Control Panel");
        buttonGroup.appendChild(loadParamButton);

        // Focus Camera Button
        var focusCameraButton = document.createElement("button");
        focusCameraButton.classList.add("btn", "btn-secondary", "simulationTableFocusCameraBtn");
        var focusIconSpan = document.createElement("span");
        focusIconSpan.classList.add("oi", "oi-aperture");
        focusCameraButton.appendChild(focusIconSpan);
        focusCameraButton.onclick = function() {
            self.focusCameraCallback(i);
        };
        focusCameraButton.setAttribute("data-toggle", "tooltip");
        focusCameraButton.setAttribute("data-placement", "top");
        focusCameraButton.setAttribute("title", "Focus in View");
        buttonGroup.appendChild(focusCameraButton);

        // Delete Button
        var deleteButton = document.createElement("button");
        deleteButton.classList.add("btn", "btn-danger", "simulationTableDelBtn");
        var deleteIconSpan = document.createElement("span");
        deleteIconSpan.classList.add("oi", "oi-x");
        deleteButton.appendChild(deleteIconSpan);
        deleteButton.onclick = function() {
            self.deleteCallback(i);
        };
        deleteButton.setAttribute("data-toggle", "tooltip");
        deleteButton.setAttribute("data-placement", "top");
        deleteButton.setAttribute("title", "Delete (cannot be undone!)");
        buttonGroup.appendChild(deleteButton);


        var buttonTd = document.createElement("td");
        buttonTd.appendChild(buttonGroup);
        row.appendChild(buttonTd);


        // Color
        var colorCircleTd = document.createElement("td");
        colorCircleTd.classList.add("color-circle-container");
        var colorCircle = document.createElement("div");
        colorCircleTd.appendChild(colorCircle);
        colorCircle.classList.add("color-circle");
        colorCircle.style.backgroundColor = colors[i % colors.length];

        row.appendChild(colorCircleTd);


        self.tableBody.appendChild(row);
    });

    $(function () {
      $('[data-toggle="tooltip"]').tooltip()
    });

};

/**
 * Activate a given simulation.
 *
 * @param {number} idx - index of simulation to active
 * @deprecated
 * @ignore
 */
SimuTable.prototype.setActive = function(idx) {
    if (this.activeSet >= 0) {
        this.tableBody.children[this.activeSet].classList.remove("table-primary");
        this.activeSet = idx;
        this.tableBody.children[this.activeSet].classList.add("table-primary");
    }
}

/**
 * Clears all rows.
 */
SimuTable.prototype.clear = function() {
    while (this.tableBody.firstChild) {
        this.tableBody.removeChild(this.tableBody.firstChild);
    }
}

/**
 * Refreshes (i.e. clears and repopulates).
 * @deprecated
 * @ignore
 */
SimuTable.prototype.refresh = function() {
    this.clear();
    this.populate();
}
