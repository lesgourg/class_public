/** @module paramGui */

/**
 * Initial Condition Type. Either Gaussian (0) or Scale Invariant (1).
 *
 * @typedef {number} IC_TYPE
 */
const IC_TYPE = {
    GAUSSIAN: 0,
    SCALE_INVARIANT: 1,
};

/**
 * Set of initial parameters providing the data structure for the ParameterGui.
 * Reads default parameters from {@link DEFAULT_INITIAL} and default redshifts
 * from {@link DEFAULT_REDSHIFTS}.
 *
 * @typedef {Object} ParamSet
 * @constructor
 *
 * @property {number} scale - side length of 2D cross section of universe in Mpc
 * @property {module:paramGui~IC_TYPE} initialType - either
 * {@link IC_TYPE.GAUSSIAN} or
 * {@link IC_TYPE.SCALE_INVARIANT}
 * @property {bool} limitModes - true if modes are to be limited, false otherwise
 * @property {number} minMode - minimum mode wave number in 1 / Mpc
 * @property {number} minMode - maximum mode wave number in 1 / Mpc
 * @property {number} ns - primordial tilt n_s
 * @property {number} resolution - the side length of the mesh in number of vertices
 *
 * @param {module:RSI~initialConditionSetCallback} setInitialConditions - callback
 */
function ParamSet(setInitialConditions) {
    this.scale = DEFAULT_INITIAL.scale;
    this.initialType = IC_TYPE.GAUSSIAN;
    this.sigmaGauss = DEFAULT_INITIAL.sigmaGauss;
    this.limitModes = false;
    this.minMode = 1;
    this.maxMode = 10;
    this.ns = DEFAULT_INITIAL.n_s;
    this.resolution = DEFAULT_INITIAL.resolution;

    var self = this;

    var initialCallback = function() {
        setInitialConditions(self);
    };

    this.setInitialSI = initialCallback;
    this.setInitialGauss = initialCallback;

    this.redshift = Object.assign([], DEFAULT_REDSHIFTS);
    this.modal = new RedshiftModal(
        this.redshift,
        function(redshift) {
            self.redshift = redshift;
        }
    );

    this.configureRedshift = function() {
        self.modal.show();
    };
};

/**
 * Returns an object according to the communication protocol
 * to be sent as JSON object to the server specifying the initial
 * conditions.
 *
 * @return {Object} the object to be sent to the server
 */
ParamSet.prototype.serializeInitial = function() {
    var dataType = (this.initialType === IC_TYPE.GAUSSIAN) ? "Gaussian" : "SI";
    var SILimit;
    if (this.initialType === IC_TYPE.GAUSSIAN || !this.limitModes)
        SILimit = "None";
    else
        SILimit = [this.minMode, this.maxMode];

    var result = {
        type: "Initial",
        params: {
            xScale: Math.floor(this.scale),
            resolution: this.resolution,
            initialDataType: dataType,
            sigma: this.sigmaGauss,
            SILimit: SILimit,
            n_s: this.ns,
            redshift: this.redshift,
        }
    };

    return result;
};

/***************************************************************************************/

/**
 * Represents the set of cosmological parameters (that are
 * not initial state parameters). Reads the set of parameters
 * exposed to the user from {@link COSMOLOGICAL_PARAMETER_LIST}.
 *
 * @typedef {Object} CosmoParams
 * @constructor
 * @param {module:RSI~cosmologicalParametersSetCallback}
 */
function CosmoParams(cosmoCallback) {
    this.config = COSMOLOGICAL_PARAMETER_LIST;
    this.additionalParameters = ADDITIONAL_CLASS_PARAMETERS;

    this.parameters = this.buildParameterObject();
    this.cosmoCallback = cosmoCallback;
}

/**
 * Builds the object which is passed to dat.gui.
 */
CosmoParams.prototype.buildParameterObject = function() {
    var result = {};
    for (var entry of this.config) {
        result[entry.name] = entry.default;
    }
    return result;
};

/**
 * Resets all parameters to their default value.
 */
CosmoParams.prototype.resetParameters = function() {
    for (var entry of this.config) {
        this.parameters[entry.name] = entry.default;
    }
};

/**
 * This is called once the user presses the button to set the cosmological parameters
 * (which coincides with the button to start the calculation).
 * This will call the {@link module:RSI~cosmologicalParametersSetCallback} instance
 * passed in the {@link CosmoParams} constructor.
 */
CosmoParams.prototype.setCosmological = function() {
    var result = {
        type: "Cosmo",
        params: this.serialize(),
    };
    this.cosmoCallback(this, result);
};

/**
 * Converts the user specified cosmological parameters into a representation
 * that can be sent to the server.
 *
 * @return {Object} object to be sent to the server
 */
CosmoParams.prototype.serialize = function() {
    var serialized = {};
    for (var entry of this.config) {
        if (entry.hasOwnProperty("replaceBy")) {
            var replacement = entry.replaceBy(this.parameters);
            for (var key in replacement) {
                serialized[key] = replacement[key];
            }
        }
        else {
            serialized[entry.name] = this.parameters[entry.name];
        }
    }
    for (var key in this.additionalParameters) {
        serialized[key] = this.additionalParameters[key];
    }

    return serialized;
};

CosmoParams.prototype.loadParameters = function(parameters) {
    var self = this;
    Object.keys(parameters).forEach(function(key) {
        self.parameters[key] = parameters[key];
    });
};


/**
 * This class provides the graphical user interface that controls the inital
 * + cosmological parameters of the simulation as well as other properties,
 * such as the animation rate, whether or not grids are shown, etc.
 *
 * @typedef {Object} ParameterGui
 * @constructor
 * @param {module:RSI~initialConditionSetCallback} setInitialConditions
 * @param {module:RSI~cosmologicalParametersSetCallback} setCosmologicalParams
 */
function ParameterGui(setInitialConditions, setCosmologicalParams) {
    this.gui = new dat.GUI({ width: 300 });
    this.gui.domElement.id = "datgui";
    this.paramSet = new ParamSet(setInitialConditions);
    this.cosmoParams = new CosmoParams(setCosmologicalParams);

    this.config = {
        notFlat: true,
        amplitude: 25.0,
        resetCamera: function() {controls.reset()},
        orthographicMode: false,
        showAxes: true,
        dashedAxes: true,
        showGrid: true,
        animationSpeed: 25,
        backgroundColor: "#4a4a4a"
    };


    this.buildFolders();
}

/**
 * Create folders for initial condition and cosmological parameters
 * and populate them.
 * Called via the constructor of {@link ParameterGui}.
 */
ParameterGui.prototype.buildFolders = function() {
    this.initialFolder = this.gui.addFolder("Initial Condition");
    this.cosmoFolder = this.gui.addFolder("Cosm. Parameters");

    this.buildInitialFolder();
    this.buildCosmoFolder();
    this.buildDisplayFolder();
    this.buildControlFolder();

    this.initialFolder.open();
    this.siFolder.open();
    this.switchToScaleInvariant();
    this.cosmoFolder.open();
};

/**
 * Build folder containing initial condition parameters.
 */
ParameterGui.prototype.buildInitialFolder = function() {
    this.initialFolder.add(this.paramSet, "scale")
        .min(100).max(1600).name("Scale [Mpc]");
    this.initialFolder.add(this.paramSet, "resolution")
        .min(64).max(1024).step(1).name("Resolution");
    this.initialFolder.add(this.paramSet, "configureRedshift")
        .name("Configure Redshifts");

    var self = this;
    this.gaussianFolder = this.initialFolder.addFolder("Gaussian Initial Conditions");
    this.sigmaGaussCtl = this.gaussianFolder
        .add(this.paramSet, "sigmaGauss").min(1).max(30).name("&sigma; [Mpc]");
    this.gaussianFolder
        .add(this.paramSet, "setInitialGauss").name("Set Initial Condition");
    this.gaussianFolder.domElement.onclick = function() {
        self.switchToGaussian();
    };

    var siButtonString = "Set Initial Condition";

    // Create SI folder
    this.siFolder = this.initialFolder.addFolder("Scale Invariant Initial Conditions");
    // Add fields
    this.siFolder
        .add(this.paramSet, "ns").name("n<sub>s</sub>");
    this.limitModesCtl = this.siFolder
        .add(this.paramSet, "limitModes").name("Limit Modes");
    this.limitModesCtl.onChange(function (limitModes) {
        if (limitModes) {
            self.siFolder.remove(self.set_initial_si);
            self.minModeCtl = self.siFolder
                .add(self.paramSet, "minMode").name("Min. Mode [1/Mpc]");
            self.maxModeCtl = self.siFolder
                .add(self.paramSet, "maxMode").name("Max. Mode [1/Mpc]");
            self.set_initial_si = self.siFolder
                .add(self.paramSet, "setInitialSI")
                .name(siButtonString);
        } else {
            self.siFolder.remove(self.minModeCtl);
            self.siFolder.remove(self.maxModeCtl);
        }
    });
    this.set_initial_si = this.siFolder
        .add(this.paramSet, "setInitialSI").name(siButtonString);
    this.siFolder.domElement.onclick = function() {
        self.switchToScaleInvariant();
    };
};

/**
 * Build folder containing cosmological parameters.
 */
ParameterGui.prototype.buildCosmoFolder = function() {
    for (var param of this.cosmoParams.config) {
        var step = param.hasOwnProperty("step") ? param.step : 0.001;
        this.cosmoFolder
        .add(this.cosmoParams.parameters, param.name)
        .step(step)
        .min(param.min)
        .max(param.max)
        .name(param.displayName);
    }

    // Rather ugly solution but necessary since calling `listen()` on the
    // config attributes above to auto-update them on value change
    // prevents user from manually entering a value into the text boxes.
    // NOTE: This relies on the internal implementation of dat.gui and is
    // therefore not guaranteed to be compatible with future versions of dat.gui.
    var self = this;
    this.resetConfig = {
        resetCosmological: function() {
            self.cosmoParams.resetParameters();
            self.refreshCosmoFolder();
        }
    };

    this.cosmoFolder.add(this.resetConfig, "resetCosmological").name("Reset Parameters");
    this.cosmoFolder.add(this.cosmoParams, "setCosmological").name("Start Calculation");
};

ParameterGui.prototype.refreshCosmoFolder = function() {
    this.gui.__folders["Cosm. Parameters"].__controllers.forEach(function(ctrl) {
        ctrl.updateDisplay();
    });
};

ParameterGui.prototype.loadCosmoParameters = function(parameters) {
    this.cosmoParams.loadParameters(parameters);
    this.refreshCosmoFolder();
};

ParameterGui.prototype.buildDisplayFolder = function() {
    displaySettings = this.gui.addFolder("Display Settings");

    var self = this;
    displaySettings
        .add(this.config, 'notFlat')
        .name("Show in 3D")
        .onChange(function() {
            simulationManager.setUniform("notFlat", self.config.notFlat);
        });

    displaySettings
        .add(this.config, "showGrid")
        .name("Show Grid")
        .onChange(function() {
            simulationManager.updateGrids();
        });

    displaySettings
        .add(this.config , 'amplitude', 0, 50)
        .name("Amplitude")
        .onChange(function() {
            simulationManager.setUniform("amplitude", self.config.amplitude);
        });

};

/**
 * Build the control folder, containing properties like the animation rate,
 * background color etc.
 */
ParameterGui.prototype.buildControlFolder = function() {
    var self = this;

    this.gui.add(this.config, "animationSpeed", 1, 60).name("FPS").listen().onChange(function() {
        fps = self.config.animationSpeed;
    });

    // this.gui.add(this.config, "orthographicMode")
    //     .name("Orthographic Camera")
    //     .onChange(toggleCameraMode);

    this.gui.add(this.config, "resetCamera"). name("Reset Camera");

    this.gui.addColor(this.config, "backgroundColor").name("Background Color")
    .onChange(function (value) {
        scene.background = new THREE.Color(value);
    });
};

/**
 * Switch to the menu for Gaussian initial conditions
 */
ParameterGui.prototype.switchToGaussian = function() {
    if (this.paramSet.initialType === IC_TYPE.SCALE_INVARIANT) {
        this.siFolder.close();
        this.paramSet.initialType = IC_TYPE.GAUSSIAN;
    }
};

/**
 * Switch to the menu for scale invariant initial conditions
 */
ParameterGui.prototype.switchToScaleInvariant = function() {
    if (this.paramSet.initialType === IC_TYPE.GAUSSIAN) {
        this.gaussianFolder.close();
        this.paramSet.initialType = IC_TYPE.SCALE_INVARIANT;
    }
};
