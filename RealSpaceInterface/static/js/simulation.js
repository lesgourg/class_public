/** @const {string[]} */
var QUANTITIES = ["d_g", "d_b", "d_ur", "d_g/4 + psi", "d_cdm"];
/** @const {string[]} */
var QUANTITY_LABELS =  [
    "&delta;<sub>&gamma;</sub>",
    "&delta;<sub>b</sub>",
    "&delta;<sub>&nu;</sub>",
    "&delta;<sub>&gamma;</sub>/4 + &psi;",
    "&delta;<sub>cdm</sub>"
];

/**
 * @typedef {Object} Simulation
 * @typedef {Object} SimulationManager
 * @typedef {Object} QuantityObject
 */

/**
 * Represents a single simulation.
 *
 * @constructor
 */
function Simulation() {
    /*
    These objects (once filled) can be accessed like this:
    transfer[QUANTITY][FRAME][DATAPOINT]
    */
    this.real = this.createQuantityObject();
    this.transfer = this.createQuantityObject();
    this.fourier = this.createQuantityObject();

    this.params = null;

    // Cl's
    this.l = null;
    this.tCl = null;

    this.mPk = null;

    this.frameCount = 0;
}

/**
 * Create a new quantity object, i.e. an object whose
 * keys are the quantities of {@link QUANTITIES} and values
 * empty arrays.
 */
Simulation.prototype.createQuantityObject = function() {
    var obj = {};
    for (var qty of QUANTITIES) {
        obj[qty] = [];
    }
    return obj;
}

/**
 * Set the cosmological parameters associated with this simulation.
 *
 * @param {Object} params - The parameter object
 */
Simulation.prototype.setParameters = function(params) {
    this.params = params;
};

/**
 * Add a frame to this simulation.
 *
 * @param {QuantityObject} real - real data
 * @param {QuantityObject} transfer - transfer function data
 * @param {QuantityObject} fourier - fourier data (not used)
 */
Simulation.prototype.addFrame = function(real, transfer, fourier) {
    for (var qty of QUANTITIES) {
        this.real[qty].push(real[qty]);
        this.transfer[qty].push(transfer[qty]);
        this.fourier[qty].push(fourier[qty]);
    }
    this.frameCount++;
};

/**
 * Add a frame to this simulation.
 *
 * @param {Float32Array} real - real data
 * @param {Float32Array} transfer - transfer function data
 * @param {Float32Array} fourier - fourier data (not used)
 */
Simulation.prototype.addInitial = function(real, transfer, fourier) {
    for (var qty of QUANTITIES) {
        this.real[qty].push(real);
        this.transfer[qty].push(transfer);
        this.fourier[qty].push(fourier);
    }
    this.frameCount++;
}

/**
 * Get the real data of `quantity` at frame `index`.
 * @param {string} quantity
 * @param {number} index - frame number
 * @return {Float32Array} real data
 */
Simulation.prototype.getReal = function(quantity, index) {
    return this.real[quantity][index];
};

/**
 * Get the transfer function data of `quantity` at frame `index`.
 * @param {string} quantity
 * @param {number} index - frame number
 * @return {Float32Array} transfer function array
 */
Simulation.prototype.getTransfer = function(quantity, index) {
    console.assert(index < this.transfer.length);
    return this.transfer[quantity][index];
};

/**
 * Get the Fourier data of `quantity` at frame `index`.
 * @deprecated Unused
 * @param {string} quantity
 * @param {number} index - frame number
 * @return {Float32Array} Fourier data
 */
Simulation.prototype.getFourier = function(quantity, index) {
    console.assert(index < this.fourier.length);
    return this.fourier[quantity][index];
};

/**
 * Get the number of frames in this simulation.
 * @return {number} frame count
 */
Simulation.prototype.getFrameCount = function() {
    return this.frameCount;
};

/*************************************************************************************************
*************************************************************************************************/

/**
 * A container class which is capable of managing multiple simulations.
 * Provides an abstraction over dealing with individual simulations manually.
 *
 * @constructor
 * @param {THREE.Texture} colorMap - A texture containing the initial color map
 *
 * @property {Simulation[]} simulations - An array of all stored simulations
 * @property {number} frame - The currently active frame
 * @property {number[]} active - The indices of active simulations
 * @property {THREE.Group} planeGroup - The object group of planes in 'regular' mode
 * @property {THREE.Group} collisionPlaneGroup - The group containing the meshes used
 *                                               for mouse picking
 * @property {number} resolution - mesh resolution
 * @property {number} PLANE_SIZE - size of planes in world space
 * @property {number} PLANE_MARGIN - margin between planes in world units
 * @property {number} GRID_DIVISIONS - number of grid divisions
 *
 * @property {THREE.Texture} colorMap - active color map
 * @property {Simulation} targetSimulation - reference to the (unfinished) simulation object that is still receiving data
 *
 * @property {Float32Array} realInit - initial state for real data, that is the same for all quantities and passed to every newly created instance of {@link Simulation}
 * @property {GridLayout} layout - layout manager which calculates the position of the individual planes
 * @property {THREE.Group} activeBoxes - group of color outlines implemented as BoxHelpers
 *
 * @property {bool} singleMode - true if application is in single mode, i.e. all quantities of a _single_ simulation are shown side by side. false otherwise.
 * @property {number} singleSimulationIndex - index of simulation that is shown in single mode (if any)
 * @property {THREE.Group} singleGroup - object group of planes and grids in single mode
 * @property {THREE.Mesh[]} singleMeshes - array of plane meshes in single mode
 * @property {THREE.Mesh[]} singleCollisionMeshes -
 *      array of low resolution plane meshes in single mode used for mouse picking
 * @property {THREE.Mesh[]} singleGrids - array of grids in single mode
 */
function SimulationManager(colorMap) {
    this.simulations = [];

    this.frame = 0;

    this.active = []; // List of active (i.e. visible) simulations
    this.planeGroup = new THREE.Group();
    this.collisionPlaneGroup = new THREE.Group();

    // Definition of constants
    this.resolution = 200;
    this.PLANE_SIZE = 500; // Size of Ground Plane in world units
    this.PLANE_MARGIN = 50;
    this.GRID_DIVISIONS = 20;

    this.colorMap = colorMap;

    // active simulation, i.e. simulation that is (still) receiving data
    this.targetSimulation = null;

    this.realInit = null;
    this.transferInit = null;
    this.fourierInit = null;

    this.layout = new GridLayout(1, this.PLANE_SIZE, this.PLANE_MARGIN);

    // Group of highlight boxes
    this.activeBoxes = new THREE.Group();

    // Single Mode = Showing only 1 simulation, but all 4 quantities at once
    this.singleMode = false;
    this.singleSimulationIndex = 0;
    this.singleGroup = new THREE.Group();
    this.singleMeshes = [];
    this.singleCollisionMeshes = [];
    this.singleGrids = [];
}

/**
 * Test whether the user hovers over any of the simulations and react accordingly,
 * i.e. create a colored outline in the appropriate color.
 * If in single mode, also display the quantity of the currently hovered plane
 * in the status bar.
 *
 * @param {THREE.Raycaster} raycaster - raycaster to be used for intersection checking
 */
SimulationManager.prototype.mousePick = function(raycaster) {
    var self = this;
    this.activeBoxes.children.forEach(function(box) {
        self.activeBoxes.remove(box);
    });
    $("#displayHover").text("");

    if (this.singleMode) {
        var intersects = raycaster.intersectObjects(this.singleCollisionMeshes);
        for (var i = 0; i < intersects.length; i++) {
            // Find corresponding simulation
            var qtyIdx = intersects[i].object.index;
            var box = new THREE.BoxHelper(this.singleCollisionMeshes[qtyIdx]);
            box.material.linewidth = 5;
            box.material.color = new THREE.Color(colors[qtyIdx % colors.length]);
            this.activeBoxes.add(box);
            $("#displayHover").html(" | " + QUANTITY_LABELS[qtyIdx]);
        }
    } else {
        var intersects = raycaster.intersectObjects(this.getCollisionGroup().children);
        for (var i = 0; i < intersects.length; i++) {
            // Find corresponding simulation
            var simulationIndex = intersects[i].object.index;
            var simulation = simulationManager.getSimulations()[simulationIndex];
            var box = new THREE.BoxHelper(simulation.collisionMesh);
            box.material.linewidth = 5;
            box.material.color = new THREE.Color(colors[simulationIndex % colors.length]);
            this.activeBoxes.add(box);
        }
    }
};

/**
 * @return {THREE.Mesh[]} array of active boxes (colored outlines).
 */
SimulationManager.prototype.getActiveBoxes = function() {
    return this.activeBoxes;
};

/**
 * Create a new set of uniforms for the shader material for a new {@link Simulation} instance.
 * @private
 * @return {Object} uniform object
 */
SimulationManager.prototype.createUniforms = function() {
    return {
        notFlat: {
            type: 'i',
            value: parameterGui.config.notFlat
        },
        amplitude: {
            type: 'f',
            value: parameterGui.config.amplitude
        },
        cmap: {
            type: "t",
            value: this.colorMap,
        },
    };
}

/**
 * Set the resolution for future new simulations.
 * @parameter {number} resolution - resolution (side length of mesh as number of vertices)
 */
SimulationManager.prototype.setResolution = function(resolution) {
    // If updating, also update this.getResolution
    this.resolution = resolution;
    if (this.singleMode && this.singleGroup.children.length > 0 ||
        !this.singleMode && this.planeGroup.children.length  > 0) {
            console.warn("Resolution set to " + resolution + ", but there are existing planes!");
    }
};

/**
 * @return {number} resolution
 */
SimulationManager.prototype.getResolution = function() {
    return this.resolution;
};

/**
 * Set a new color map.
 *
 * @param {THREE.Texture} texture - new color map to be set
 */
SimulationManager.prototype.setColorMap = function(texture) {
    this.colorMap = texture;

    for (var plane of this.planeGroup.children) {
        if (plane.material.uniforms)
            plane.material.uniforms.cmap.value = texture;
    }
    for (var plane of this.singleMeshes) {
        plane.material.uniforms.cmap.value = texture;
    }
};

/**
 * Create a new material for a {@link Simulation} instance.
 * @private
 *
 * @return {THREE.ShaderMaterial} ShaderMaterial for simulation
 */
SimulationManager.prototype.createMeshMaterial = function() {
    var uniforms = this.createUniforms();
    return new THREE.ShaderMaterial({
        uniforms: uniforms,
        vertexShader: document.getElementById('vertexshader').textContent,
        fragmentShader: document.getElementById('fragmentshader').textContent,
    });
};

/**
 * Create a new simulation and initialize it with the stored initial state.
 */
SimulationManager.prototype.createSimulation = function() {
    var sim = new Simulation();
    this.simulations.push(sim);
    this.targetSimulation = sim;
    sim.addInitial(this.realInit, this.transferInit, this.fourierInit);
    this.activateSimulation(sim);
};


/**
 * @return {Simulation} the current target simulation
 */
SimulationManager.prototype.getTarget = function() {
    return this.targetSimulation;
};

/**
 * @return {bool} whether or not there exists a target simulation.
 */
SimulationManager.prototype.hasTarget = function() {
    return this.targetSimulation !== null;
};

/**
 * Finalize target (should be called once all data for a simulation has been received).
 */
SimulationManager.prototype.finalizeTarget = function() {
    this.targetSimulation = null;
};

/**
 * Destroys the current target. This should be called in case an error occurs during
 * the data transfer from the server.
 */
SimulationManager.prototype.destroyTarget = function() {
    this.delete(this.simulations.indexOf(this.targetSimulation));
    this.targetSimulation = null;
};

/**
 * Update the grids after the user toggles grid visibility in the control panel.
 */
SimulationManager.prototype.updateGrids = function() {
    var self = this;
    if (this.singleMode) {
        this.singleGrids.forEach(function(grid) {
            grid.visible = parameterGui.config.showGrid;
        });
    } else {
        this.active.map(function(idx) { return self.simulations[idx]; }).forEach(function(sim) {
            sim.grid.visible = parameterGui.config.showGrid;
        });
    }
};

/**
 * Enter single mode for the given simulation index, i.e. show all quantities
 * for that simulation side by side.
 *
 * @param {number} simulationIdx - index of simulation to view in single mode
 */
SimulationManager.prototype.enableSingleMode = function(simulationIdx) {
    this.singleMode = true;
    this.singleSimulationIndex = simulationIdx;

    var self = this;
    QUANTITIES.forEach(function(qty, qtyIdx) {
        var plane = self.createPlane();
        var collisionPlane = self.createCollisionPlane(qtyIdx);
        var grid = new THREE.GridHelper(self.PLANE_SIZE, self.GRID_DIVISIONS);
        self.singleMeshes.push(plane);
        self.singleCollisionMeshes.push(collisionPlane);
        self.singleGrids.push(grid);

        self.singleGroup.add(plane);
        self.singleGroup.add(grid);
    });

    this.refresh();
    this.loadFrame("d_g", this.frame); // Quantity is actually unused here and merely serves as a dummy
};

/**
 * Leave single mode.
 */
SimulationManager.prototype.disableSingleMode = function() {
    this.singleMode = false;

    for (var i = 0; i < this.singleGroup.children.length; ++i) {
        this.singleGroup.remove(this.singleGroup.children[i]);
    }


    for (var i = 0; i < QUANTITIES.length; ++i) {
        this.singleMeshes[i].geometry.dispose();
        this.singleCollisionMeshes[i].geometry.dispose();
        this.singleGrids[i].geometry.dispose();
    }

    this.singleMeshes = [];
    this.singleCollisionMeshes = [];
    this.singleGrids = [];

    this.refresh();
};

/**
 * @return {THREE.Group} {@link SimulationManager#singleGroup}
 */
SimulationManager.prototype.getSingleGroup = function() {
    return this.singleGroup;
};

/**
 * Activate the simulation for the given index.
 *
 * @param {number} index - index of the simulation to activate
 */
SimulationManager.prototype.activate = function(index) {
    console.log("Activating simulation #" + index);

    this.active.push(index);

    this.simulations[index].mesh = this.createPlane();
    this.simulations[index].collisionMesh = this.createCollisionPlane(index);
    this.simulations[index].grid = new THREE.GridHelper(this.PLANE_SIZE, this.GRID_DIVISIONS);
    this.simulations[index].grid.visible = parameterGui.config.showGrid;

    this.planeGroup.add(this.simulations[index].mesh);
    this.planeGroup.add(this.simulations[index].grid);
    this.collisionPlaneGroup.add(this.simulations[index].collisionMesh);

    this.refresh();
};

/**
 * Activate the given simulation.
 *
 * @param {Simulation} sim - simulation to activate
 */
SimulationManager.prototype.activateSimulation = function(sim) {
    this.activate(this.simulations.indexOf(sim));
};

/**
 * Deactivate the simulation for the given index.
 *
 * @param {number} index - index of simulation to deactivate
 */
SimulationManager.prototype.deactivate = function(index) {
    if (!this.active.includes(index)) {
        return;
    }
    this.active.splice(this.active.indexOf(index), 1);
    var simulation = this.simulations[index];

    // Remove meshes from rendering groups
    this.planeGroup.remove(simulation.mesh);
    this.planeGroup.remove(simulation.grid);
    this.collisionPlaneGroup.remove(simulation.collisionMesh);

    // Dispose of mesh geometry
    if (simulation.mesh) {
        simulation.mesh.geometry.dispose();
        simulation.mesh.material.dispose();
        simulation.mesh = null;
    }
    if (simulation.collisionMesh) {
        simulation.collisionMesh.geometry.dispose();
        simulation.collisionMesh = null;
    }
    if (simulation.grid) {
        simulation.grid.geometry.dispose();
        simulation.grid = null;
    }

    this.refresh();
};

/**
 * Deactivate given simulation.
 *
 * @param {Simulation} sim - simulation to deactivate
 */
SimulationManager.prototype.deactivateSimulation = function(sim) {
    this.deactivate(this.simulations.indexOf(sim));
};

/**
 * Delete all simulations.
 */
SimulationManager.prototype.deleteAll = function() {
    // for (var activeIdx of this.active) {
    //     this.deactivate(activeIdx);
    // }

    // Evaluate count outside for loop since loop body modifies this.simulations
    var count = this.simulations.length;
    for (var i = 0; i < count; ++i) {
        // Always delete simulation at index #0, since by deleting simulation #0,
        // what was previously simulation #1 now becomes #0.
        this.delete(0);
    }

    this.active = [];
    this.simulations = [];
    console.assert(this.planeGroup.children.length == 0);
    this.refresh();
};

/**
 * Delete a single simulation.
 *
 * @param {number} index of simulation to delete.
 */
SimulationManager.prototype.delete = function(index) {
    if (index >= 0 && index < this.simulations.length) {
        this.deactivate(index);

        // Since a simulation was removed from the list,
        // the indices of this.active aren't accurate anymore.
        // Every active index > `index` (i.e. the index of the simulation
        // to delete) needs to be shifted down by one so that they continue
        // to refer to the correct simulation.
        this.active = this.active.map(function(idx) {
            return (idx > index) ? idx - 1 : idx;
        });

        // After this, no more reference the simulation should exist
        this.simulations.splice(index, 1);

        // Update simulationIndex on collision meshes
        this.simulations.forEach(function(simulation, index) {
            if (simulation.collisionMesh) {
                simulation.collisionMesh.index = index;
            }
        });
    }
};

/**
 * Needs to be called after simulations have been activated/deactivated
 * to re-layout the remaining (activated) planes.
 */
SimulationManager.prototype.refresh = function() {
    if (this.singleMode) {
        this.layout.count = QUANTITIES.length;
        this.layout.recalculate();

        for (var i = 0; i < QUANTITIES.length; ++i) {
            var position = this.layout.getWorldPosition(i);
            this.singleMeshes[i].position.x = position[0];
            this.singleMeshes[i].position.z = position[1];
            this.singleCollisionMeshes[i].position.x = position[0];
            this.singleCollisionMeshes[i].position.z = position[1];
            this.singleGrids[i].position.x = position[0];
            this.singleGrids[i].position.z = position[1];
        }
    } else {
        this.layout.count = this.active.length;
        this.layout.recalculate();

        var self = this;
        var objects = ["mesh", "collisionMesh", "grid"];
        this.active.forEach(function(active, gridIdx) {
            var position = self.layout.getWorldPosition(gridIdx);
            var x = position[0];
            var y = position[1];

            var sim = self.simulations[active];
            objects.forEach(function(obj) {
                sim[obj].position.x = x;
                sim[obj].position.z = y;
            });
        });
    }
};

/**
 * Create a new mesh for a simulation.
 * @private
 * @return {THREE.Mesh} mesh
 */
SimulationManager.prototype.createPlane = function() {

    var bufferGeometry = new THREE.PlaneBufferGeometry(
        this.PLANE_SIZE, this.PLANE_SIZE,
        this.resolution - 1, this.resolution - 1
    );
    var result = new THREE.Mesh(bufferGeometry, this.createMeshMaterial());

    var displacementBuffer = new Float32Array(this.resolution * this.resolution);
    var displacementAttribute = new THREE.BufferAttribute(displacementBuffer, 1);
    result.geometry.addAttribute("displacement", displacementAttribute);
    result.rotation.x = -Math.PI/2;

    return result;
};

/**
 * Create a plane for mouse hover checking for the simulation at given index.
 *
 * @param {number} index - index of simulation
 */
SimulationManager.prototype.createCollisionPlane = function(index) {
    // index is required, as it is stored as mesh attribute to allow the
    // mouse picker to identify the simulation that is associated with the collision mesh.
    var bufferGeometry = new THREE.PlaneBufferGeometry(this.PLANE_SIZE, this.PLANE_SIZE);
    var result = new THREE.Mesh(bufferGeometry);
    result.rotation.x = -Math.PI/2;
    // Apply transformation such that ray caster detects plane in correct orientation
    // even if it is not added to the scene
    result.updateMatrix();
    result.geometry.applyMatrix(result.matrix);
    result.rotation.set(0, 0, 0);
    result.updateMatrix();
    result.index = index;
    return result;
};

/**
 * @return {Simulation[]} {@link SimulationManager#simulations}
 */
SimulationManager.prototype.getSimulations = function() {
    return this.simulations;
};

/**
 * @return {number[]} {@link SimulationManager#active}
 */
SimulationManager.prototype.getActive = function() {
    return this.active;
};

/**
 * Return a list of active simulations as tuples, each of which consists of the
 * index of the simulation and the simulation itself.
 * @return {Object[]}
 */
SimulationManager.prototype.getActiveSimulations = function() {
    var self = this;
    return this.active.map(function(idx) { return [idx, self.simulations[idx]]; });
};

/**
 * @callback SimulationManager~activeCallback
 * @param {number} activeIndex
 * @param {Simulation} activeSimulation
 */

/**
 * Iterate over active simulations.
 *
 * @param {SimulationManager~activeCallback} - a function called for each active simulation
 */
SimulationManager.prototype.forEachActive = function(callback) {
    var self = this;
    this.active.forEach(function(activeIndex) {
        callback(activeIndex, self.simulations[activeIndex]);
    });
}

/**
 * Check whether a given simulation specified by index is active.
 *
 * @param {number} index - index of simulation
 */
SimulationManager.prototype.isActive = function(index) {
    return this.active.includes(index);
};

/**
 * Update meshes to show the given frame for the given quantity.
 *
 * @param {string} quantity - quantity to display
 * @param {number} frame - frame to display
 */
SimulationManager.prototype.loadFrame = function(quantity, frame) {
    this.frame = frame;

    if (this.singleMode) {
        this.loadSingleFrame(frame);
    } else {
        var self = this;
        this.active.forEach(function(index) {
            self.loadFrameForIndex(quantity, index, frame);
        });
    }
};

/**
 * Similar to {@link SimulationManager#loadFrame}, but for single mode (i.e.
 * all quantities).
 *
 * @param {number} frame - frame to load
 */
SimulationManager.prototype.loadSingleFrame = function(frame) {
    var self = this;
    QUANTITIES.forEach(function(qty, qtyIdx) {
        var disp = self.singleMeshes[qtyIdx].geometry.attributes.displacement;
        disp.array = self.simulations[self.singleSimulationIndex].getReal(qty, frame);
        disp.needsUpdate = true;
    });
};

/**
 * Load the given frame for the given quantity for the given simulation index.
 *
 * @private
 * @param {string} quantity - quantity
 * @param {number} index - simulation index
 * @param {number} frame - frame
 */
SimulationManager.prototype.loadFrameForIndex = function(quantity, index, frame) {
    var simulation = this.simulations[index];
    var plane = simulation.mesh;

    // Update plane
    var displacementArray = plane.geometry.attributes.displacement.array;
    var newArray = simulation.getReal(quantity, frame);
    console.assert(newArray.length === plane.geometry.attributes.displacement.array.length);
    plane.geometry.attributes.displacement.array = newArray;
    plane.geometry.attributes.displacement.needsUpdate = true;
};

/**
 * Check whether there already exists a simulation for the given set of cosmological parameters.
 *
 * @param {Object} cosmoParams - object of cosmological parameters
 * @return {bool}
 */
SimulationManager.prototype.wasAlreadySimulated = function(cosmoParams) {
    for (var simulation of this.simulations) {
        // Require that cosmo params have already been set.
        // This is not the case for the first simulation.
        if (simulation.params) {
            var sProps = Object.getOwnPropertyNames(simulation.params);
            var nProps = Object.getOwnPropertyNames(cosmoParams);

            var lengthMatch = sProps.length == nProps.length;
            var propMatch = true;

            for (var i = 0; i < sProps.length; i++) {
                var pName = sProps[i];
                if (simulation.params[pName] !== cosmoParams[pName]) {
                    propMatch = false;
                    break;
                }
            }

            if (lengthMatch && propMatch) {
                return true;
            }
        }
    }
    return false;
};

/**
 * Get number of frames.
 * @return {number} frame count
 */
SimulationManager.prototype.getFrameCount = function() {
    if (this.simulations.length > 0)
        return this.simulations[0].getFrameCount();
    else
        return 0;
};

/**
 * @return {THREE.Group} {@link SimulationManager#planeGroup}
 */
SimulationManager.prototype.getGroup = function() {
    return this.planeGroup;
};

/**
 * @return {THREE.Group} {@link SimulationManager#collisionPlaneGroup}
 */
SimulationManager.prototype.getCollisionGroup = function() {
    return this.collisionPlaneGroup;
};

/**
 * Set a uniform value for all materials.
 *
 * @param {string} key - uniform key
 * @param {Object} value - value for key
 */
SimulationManager.prototype.setUniform = function(key, value) {
    var self = this;
    if (this.singleMode) {
        this.singleMeshes.forEach(function(mesh) {
            mesh.material.uniforms[key].value = value;
        });
    } else {
        this.active.forEach(function(idx) {
            var simulation = self.simulations[idx];
            if (simulation.mesh) {
                simulation.mesh.material.uniforms[key].value = value;
            }
        });
    }
};

/**
 * Get transfer function arrays of all simulations for given quantity and frame.
 *
 * @param {string} quantity
 * @param {number} frame
 *
 * @return {Float32Array[]} array of transfer function data for all simulations
 */
SimulationManager.prototype.getTransferData = function(quantity, frame) {
    var self = this;
    return this.active.map(function(idx) { return self.simulations[idx].transfer[quantity][frame]; });
};

/**
 * Get transfer function data for given quantity and frame of active single-mode simulation.
 *
 * @param {string} quantity
 * @param {number} frame
 *
 * @return {Float32Array}
 */
SimulationManager.prototype.getTransferDataOfSingle = function(quantity, frame) {
    var self = this;
    return self.simulations[self.singleSimulationIndex].transfer[quantity][frame];
};

/**
 * Get index and simulation itself of currently active single mode simulation.
 *
 * @return {Object[]}
 */
SimulationManager.prototype.getSingleSimulation = function() {
    return [this.singleSimulationIndex, this.simulations[this.singleSimulationIndex]];
};
