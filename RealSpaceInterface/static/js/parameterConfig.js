/*
The parameters the application exposes to the user are configured using this list.
Each parameter is represented by an entry in that list.
An entry is an object containing some required and some optional fields.

The required fields are:
------------------------
- `name`: 			The `name` property is the one that will be used to pass the
                    parameter value to CLASS.

- `displayName`: 	The `displayName` property will be displayed as the label
                    of the property in both the control panel (on the left of the
                    application) and in the simulation list.
                    HTML is allowed. Subscripts use the <sub> tag.

- `min`: 			Minimum parameter value.

- `max`: 			Maximum parameter value.

- `default`: 		Default parameter value.

The optional fields are:
------------------------
- `step`: 			The increment by which the parameter will be increased using
                    the sliders in the control panel. Defaults to 0.001.

- `replaceBy`: 		It might be necessary to pass parameters to CLASS
                    which are functions of other parameters, e.g. the user might
                    set a value of `omega_m` but CLASS requires `omega_cdm` to be
                    passed. In this case, a function accepting a dictionary of
                    all the parameters and returning a dictionary containining
                    the key-value pairs that are supposed to be passed to CLASS
                    (see example for `omega_m` below).
*/

/** @global */
var COSMOLOGICAL_PARAMETER_LIST = [
    {
        name: "h",
        displayName: "h",
        min: 0.0,
        max: 2.0,
        default: 0.67556,
    },
    {
        name: "omega_b",
        displayName: "&omega;<sub>b</sub>",
        min: 0.0,
        max: 1.0,
        default: 0.022,
    },
    {
        name: "omega_m",
        displayName: "&omega;<sub>m</sub>",
        min: 0.0,
        max: 1.0,
        default: 0.142,
        replaceBy: function(parameters) {
            return {
                "omega_cdm":  parameters.omega_m - parameters.omega_b
            };
        },
    },
    {
        name: "Omega_k",
        displayName: "&Omega;<sub>k</sub>",
        min: -0.2,
        max: 0.2,
        default: 0.0,
    },
    {
        name: "N_ur",
        displayName: "N<sub>eff</sub>",
        min: 0.0,
        max: 30.0,
        default: 3.046
    },
    {
        name: "w0_fld",
        displayName: "w<sub>0,fld</sub>",
        min: -2.0,
        max: 0.0,
        default: -1.0,
    },
    {
        name: "wa_fld",
        displayName: "w<sub>a,fld</sub>",
        min: -1.0,
        max: 1.0,
        default: 0.0
    },
];


/**
 * It is possible to pass additional parameters to CLASS using
 * the following object.
 *
 * @global
 */
var ADDITIONAL_CLASS_PARAMETERS = {
    "Omega_Lambda": 0,
    "YHe": 0.25,
};


/**
 * Default redshifts.
 * Depending on the value of log, the individual arrays will be created either
 * as
 * <pre>
 * numpy.logspace(from, to, points) (for log = true)
 *  or
 * numpy.linspace(from, to, points) (for  log = false).
 * </pre>
 *
 * The resulting individual arrays will be concatenated and duplicates removed.
 *
 * @global
 */
var DEFAULT_REDSHIFTS = [
    {
        from: 1000000,
        to:     10000,
        log:     true,
        points:    60,
    },
    {
        from:   10000,
        to:      1089,
        log:     true,
        points:    80,
    },
    {
        from:    1089,
        to:      0.01,
        log:    false,
        points:    40,
    },
];

/**
 * Default Initial Condition Values.
 *
 * @global
 */
var DEFAULT_INITIAL = {
    scale: 400,
    resolution: 200,
    sigmaGauss: 10.0,
    n_s: 0.96
};

