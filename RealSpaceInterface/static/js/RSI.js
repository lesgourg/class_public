/** @module RSI */

/**
 * Display a warning message if web browser does not support WebGL.
 */
if (!Detector.webgl) {
    Detector.addGetWebGLMessage();
    document.getElementById('container').innerHTML = "";
}

//////////////////////
// Global Variables //
//////////////////////

/**
 * Indicates whether initial condition has been set.
 * @type {boolean}
 */
var initialSet = false; //initialConditions are not set
/**
 * Indicates whether cosmological parameters have been set.
 * @type {boolean}
 */
var cosmoSet = false; //cosmological parameters are not set
/**
 * indicates whether calculation is running (i.e. initial or cosmological)
 * @type {boolean}
 */
var calculationRunning = false;
/**
 * indicates whether initial state has been received.
 * @type {boolean}
 */
var receivedInitial = false;

/**
 * array containing redshift values.
 * @type {number[]}
 */
var redshift = ["Initial"];
/**
 * True if simulation is to be run until today, else until decoupling.
 * @type {boolean}
 */
var stopAtDecoupling = true;
/**
 * The last frame (inclusive) to be displayed when simulation is run until decoupling.
 * @type {number}
 */
var decouplingFrame = 0;

/**
 * Active quantity. Must be from {@link QUANTITIES}.
 *
 * @type {string}
 */
var activeQuantity = "d_g";

var container, stats;

/**
 * @type {THREE.Camera}
 */
var camera;
/**
 * @type {THREE.Scene}
 */
var scene;
/**
 * @type {THREE.OrbitControls}
 */
var controls;
/**
 * @type {THREE.WebGLRenderer}
 */
var renderer;

/**
 * @type {THREE.Clock}
 */
var clock = new THREE.Clock();

/**
 * Size of grid divisions.
 * Set in {@link module:RSI~onInitialConditionSet}.
 * @type {string}
 */
var realScale = "0 Mpc";

/*
 * @type {SimulationManager}
 */
var simulationManager;

/**
 * List of color map objects:
 * <pre><tt>
 * [{name: &ltname&gt;, texture: &lttexture&gt;}, ...]
 * </tt><pre>
 *
 * @type {Array.<Object.<string, THREE.Texture>>}
 * */
var colorMaps = [];

/**
 * GUI responsible for controlling inital + cosmological parameters.
 * @type {module:paramGui~ParamGui}
 */
var parameterGui;

/**
 * modal dialog containing table of previous simulations and
 * their previous parameters.
 *
 * @type {SimuTable}
 * @see {SimuTable}
 */
var simuTable;

/**
 * Player controls.
 *
 * @type {PlayerPanel}
 */
var playerPanel;

/**
 * flot.js instance for transfer function plot.
 */
var transferFunctionPlot;

/**
 * flot.js instance for Cl plot.
 */
var tClPlot;

/**
 * Indicates whether plots panel is visible.
 * @type {boolean}
 */
var plotWindowVisible = true;

/**
 * Holds the list of k's for the transfer function plot.
 * @type {number[]}
 */
var kRange;

/**
 * Indicates whether websocket connection is closed.
 * @type {boolean}
 */
var closed = true;
/**
 * Current frame.
 * @type {number}
 */
var frame = 0;
/**
 * Indicates whether animation is currently running.
 * @ŧype {boolean}
 */
var animationRunning = false;

/**
 * Field of View (in degrees) of perspective camera
 * @type {number}
 */
var FoV = 45;

/**
 * true if scene is to be rendered using an orthographic camera.
 * @type {boolean}
 * @see {@link module:RSI~toggleCameraMode}
 */
var orthographicMode = false;

/**
 * Renderer used to create images / .gifs.
 * @type {THREE.WebGLRenderer}
 */
var imgRenderer = new THREE.WebGLRenderer();

/**
 * Orthographic camera used for image / .gif creation.
 * @type {THREE.OrthographicCamera}
 */
var orthoCam;

/**
 * Number of decimal digits to round redshift values to.
 * @type {number}
 */
var REDSHIFT_ROUND = 1;

/**
 * List of visually distinct colors for use in plots,
 * taken from {@link https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/}
 */
var colors = [
    "#e6194b",
    "#3cb44b",
    // "#ffe119", // yellow, hard to see on white ground
    "#0082c8",
    "#f58231",
    "#911eb4",
    // "#46f0f0", // cyan, hard to see on white ground
    "#f032e6",
    "#d2f53c",
    "#fabebe",
    "#008080",
    "#e6beff",
    "#aa6e28",
    "#fffac8",
    "#800000",
    "#aaffc3",
    "#808000",
    "#ffd8b1",
    "#000080",
    "#808080",
    "#FFFFFF",
    "#000000"
];

/**
 * Plot options for Cl's
 */
var tClPlotOptions = {
    legend: {
        show: false,
    },
    xaxis: {
        axisLabel: "&#8467;",
        // transform: Math.log,
        // inverseTransform: Math.exp,
    },
    yaxis: {
        tickFormatter: function(num, obj) {
            if (num > 0 && num < 1e-3) {
                var power = Math.floor(Math.log10(num));
                var prefix = num / Math.pow(10, power);
                return round(prefix, 2) + "&middot;10<sup>" + power + "</sup>";
            }
            return num;
        },
        axisLabel: "&#8467; (&#8467;+1) C<sub>&#8467;</sub><sup>TT</sup> / (2&pi;)",
    },
};

/**
 * Plot options for matter spectrum
 */
var mPkPlotOptions = {
    legend: {
        show: false,
    },
    xaxis: {
        transform: Math.log,
        inverseTransform: Math.exp,
        axisLabel: "k [h/Mpc]",
        ticks: [[1e-3, "10<sup>-3</sup>"], [1e-2, "10<sup>-2</sup>"], [1e-1, "10<sup>-1</sup>"], 1, 10],
    },
    yaxis: {
        transform: Math.log,
        inverseTransform: Math.exp,
        tickFormatter: function(num, obj) {
            var power = Math.floor(Math.log10(num));
            var prefix = num / Math.pow(10, power);
            var result = "";
            if (prefix != 1)
                result += round(prefix, 2) + "&middot;";
            return result + "10<sup>" + power + "</sup>";
        },
        axisLabel: "P(k) [(Mpc/h)<sup>3</sup>]",
        ticks: [1, 10, 100, 1000, 10000, 100000, 1000000],
        min: 0.1,
    },
};

/**
 * Plot options for transfer function(s)
 */
var transferPlotOptions = {
    legend: {
        show: false,
    },
    xaxis: {
        axisLabel: "k [Mpc<sup>-1</sup>]",
        transform: Math.log,
        inverseTransform: Math.exp,
        ticks: [[1e-3, "10<sup>-3</sup>"], [1e-2, "10<sup>-2</sup>"], [1e-1, "10<sup>-1</sup>"], 1],
    },
    yaxis: {
        axisLabel: "transfer function",
        min: -5,
        max: 5,
    },
};


/**
 * Raycaster for mouse picking
 * @type {THREE.Raycaster}
 */
var raycaster = new THREE.Raycaster();
/**
 * 2d vector containing mouse position.
 * @type {THREE.Vector2}
 */
var mouse = new THREE.Vector2();


/* ENTRY POINT */
$(document).ready(connect);

/**
 * Mouse Event Callback which is called on every mouse move.
 * It is required to track the mouse's position since it's required
 * to feed its position to the raycaster to determine over which mesh
 * the mouse is currently hovering.
 */
function onMouseMove(ev) {
    var rect = container.getBoundingClientRect();
    mouse.x = ((ev.clientX - rect.left) / window.innerWidth) * 2 - 1;
    mouse.y = -((ev.clientY - rect.top) / window.innerHeight) * 2 + 1;
}


/**
 * Once the websocket in {@link module:RSI~connect} has established a connection,
 * this function will be called, which initializes the application.
 */
function start() {
    init();
    step();
    requestAnimationFrame(animate);
}

/**
 * Logging function
 *
 * @param {Object} msg
 */
function log(msg) {
    console.log("D: " + msg);
}

/**
 * Entry point for the whole application.
 *
 * Establishes a connection to the WebSocket provided by the backend.
 */
function connect() {
    // connection = new SockJS('http://' + window.location.host + '/datasocket');
    connection = new WebSocket("ws://" + window.location.host + "/datasocket");

    log("Connecting...");

    connection.onopen = function() {
        log("Connected!");
        start();
    };

    connection.onmessage = messageCallback;

    connection.onclose = function() {
        log('Disconnected.');
        connection = null;
        onClose();
    };
}

/**
 * Message callback for the WebSocket connection which gets called for each
 * received message.
 *
 * @param {Object} e - message event object
 */
function messageCallback(e) {
    var message = JSON.parse(e.data);
    /**
     * Once data has been received, update the progress in the progress
     * in the progress modal and call {@link module:RSI~onData} to process
     * the received data.
     * Also, flag {@link receivedInitial} as true if no initial data has been
     * received prior.
     */
    if (message.type == 'data') {
        if (message.hasOwnProperty("progress")) {
            var percentage = (message.progress * 100).toFixed(1) + "%";
            $("#simulationProgressBar").css({width: percentage});
            $("#simulationProgressBar").text(percentage);
        }
        onData(message, !receivedInitial);
        if (!receivedInitial) {
            receivedInitial = true;
        }
    /**
     * The server sends a list of k values (at which the transfer function is sampled)
     * immediately after the client connected.
     */
    } else if (message.type === 'krange') {
        log("Received k range from server.");
        kRange = message.k;
    }
    /**
     * Handles the temperature Cl's sent by the server.
     * The Cl's (and l's) sent by the server conventionally do not contain the samples for l = 0 and l = 1.
     * They are, however, not already in the conventional form of l * (l + 1) Cl / (2 * Pi),
     * so they are processed below.
     * After that, they are added to the active simulation exposed by
     * {@link SimulationManager#getTarget}.
     */
    else if (message.type === "Cl") {
        log("Received Cl's from server.");
        var l = message.l;
        var tCl = message.tCl;

        message.tCl = l.map(function(ll, i) { return tCl[i] * ll * (ll + 1) / (2 * Math.PI); });
        simulationManager.getTarget().Cls = message;
    }
    /**
     * Handles matter power spectrum, similar to Cl's above.
     */
    else if (message.type == "mPk") {
        log("Received matter spectrum from server");
        var kh = message.kh;
        var Pkh = message.Pkh;

        simulationManager.getTarget().mPk = message;
        plotStaticIfVisible();
    }
    /**
     * The server signals end of calculations or data transfers via a "success" message
     * of one of several types.
     * Upon reception of that message the activity indicator is hidden.
     */
    else if (message.type == "success") {
        calculationRunning = false;
        hideActivityIndicator();

        if (message.sort == 'Initial') {
            initialSet = true;
            loadFrame(0);
        }
        /**
         * Once the server has confirmed that it has received cosmological parameters,
         * immediately request the simulation results.
         */
        else if (message.sort == 'Cosmo') {
            cosmoSet = true;
            $("#simulationProgressInfo").text("Receiving simulation data...");
            requestSimulationData();
        }
        /**
         * Handle the server signaling end of data transfer (for all frames).
         * The progress modal is hidden, and the list of simulations is rebuild from scratch.
         * Also, the current frame is reloaded to display the newly received data
         * at the currently active frame.
         */
        else if (message.sort == "Data") {
            log("Success receiving simulation.");
            $("#simulationProgressModal").modal("hide");
            simulationManager.finalizeTarget();
            // Update the simulation list/table
            simuTable.clear();
            var simulations = simulationManager.getSimulations();
            var active = simulationManager.getActive();
            simuTable.populate(simulationManager.getSimulations(), simulationManager.getActive());
            loadFrame(frame, true);
        }
    }
    /**
     * Server sends redshifts immediately after sending initial condition data.
     */
    else if (message.type == "redshift") {
        redshift = ["Initial"].concat(message.redshift);
    }
    /**
     * Server also sends exact z of decoupling extracted from CLASS's background parameters.
     * Also contains the index of the closest frame *before* decoupling, which this
     * application requires to stop the simulation playback at the appropriate frame.
     */
    else if (message.type == 'decoupling') {
        log("Decoupling at z = " + message.z + " (frame #" + message.frame + ")");
        // Add +1 because redshift[0] === "Initial"
        decouplingFrame = message.frame + 1;
    }
    /**
     * The server notifies the client about the resolution (i.e. the side length of the square
     * data matrix) it is about to receive.
     * This is done before an initial state is sent.
     */
    else if (message.type == "resolution") {
        log("Setting resolution to " + message.value);
        simulationManager.setResolution(message.value);
    }
    /**
     * Currently unused.
     */
    else if (message.type == "extrema") {
    }
    /**
     * In case of an exception occuring on the server side, the client is notified of that
     * and displays the transmitted stack trace to the user.
     * Also, if one exists, the target simulation (the one, which is currently receiving data)
     * is destroyed in case of an exception, to guarantee the continued functionality of the
     * application.
     */
    else if (message.type == "exception") {
        $("#exceptionModalMessage").text(message.exception);
        if (calculationRunning) {
            $("#simulationProgressModal").modal("hide");
            calculationRunning = false;
        }
        simulationManager.destroyTarget();
        $("#exceptionModal").modal("show");
    }
}

function sendParams(params,type) {
    connection.send(JSON.stringify({type:type,params:params}));
}

/**
 * Disconnect from the server.
 */
function disconnect() {
    log("disconnecting");
    connection.close();
    connection = null;
}


/**
 * Requests the transfer of the actual data (for all the redshifts) from the server.
 */
function requestSimulationData() {
    log("Requesting calculation result from server...");
    if (!calculationRunning) {
        connection.send(JSON.stringify({
            type: "Start",
            params: [],
        }));
        calculationRunning = true;
        showActivityIndicator();
    }
}

/**
 * A simple rounding implementation since JavaScript doesn't provide a built-in
 * solution.
 * This will round the given number to a maximum number of decimal places (no trailing zeros).
 *
 * @param {number} num - The number to round
 * @param {number} decimalPlaces - Maximum number of decimal places
 * @return {number} rounding result
 */
function round(num, decimalPlaces) {
    return num.toFixed(decimalPlaces).replace(/(\.)?0+$/, "");
}

/**
 * Since the data transfer between server and client (in that direction) is the main
 * bottleneck of the application, the data is mostly (except for plots, redshifts, ...)
 * encoded in base64.
 * This function serves the purpose of decing that data.
 * The data is assumed to consist of 32-bit floats. Other data formats (e.g. 64-bit floats)
 * would require a modification of this function.
 *
 * @param {string} b64 - base64 encoded data received from the server
 * @return {Float32Array} The decoded data array
 */
function b64ToFloatArray(b64) {
    var binary = window.atob(b64);
    var dv = new DataView(new ArrayBuffer(4));
    var resultLength = binary.length / 4;
    var result = new Float32Array(resultLength);

    for (var i = 0; i < resultLength; i++) {
        var idx = 4 * i;
        for (var offset = 0; offset < 4; offset++) {
            dv.setUint8(offset, binary.charCodeAt(idx + offset));
        }
        result[i] = dv.getFloat32(0, true);
    }

    return result;
}

/**
 * A data array (for a single frame) sent by the server has the following structure:
 * <pre><tt>
 * {
 *  'd_g': &lt;base64 encoded string&gt;,
 *  'd_b': &lt;base64 encoded string&gt;,
 *  ...
 * }
 * </tt></pre>
 * Hence, each field must be decoded individually,
 * which is done using {@link module:RSI~b64ToFloatArray}.
 *
 * <pre><tt>
 * {
 *  'd_g': &lt;Float32Array&gt;,
 *  'd_b': &lt;Float32Array&gt;,
 *  ...
 * }
 * </tt></pre>
 * @param {Object} obj - base64 object as described above
 * @return {Object} decoded object, as described above
 */
function decodeArray(obj) {
    return Object.keys(obj).reduce(function(acc, key) {
        acc[key] = b64ToFloatArray(obj[key]);
        return acc;
    }, {});
}

/**
 * Toggles the camera mode between an orthographic and a perspective camera
 */
function toggleCameraMode() {
    var aspect = $("#container").width() / $("#container").height();
    if (orthographicMode) {
        console.log("Switching to perspective mode");
        camera = new THREE.PerspectiveCamera(FoV, aspect, 0.1, 10000);
        controls.object = camera;
    }
    else {
        console.log("Switching to orthographic mode");
        camera = new THREE.OrthographicCamera(
            -500 * aspect, 500 * aspect,
            500, -500,
            -10000, 10000
        );
        controls.object = camera;
    }

    orthographicMode = !orthographicMode;
}

/**
 * Switches to so-called 'single mode', which, instead of displaying <i>one</i> quantity
 * for <i>multiple</i> simulations side by side, displays <i>all</i> quantities for a
 * <i>single</i> simulation side by side.
 * In order to enter single mode, following steps will occur in this function:
 * <ol>
 *  <li>Enable single mode in the {@link SimulationManager} instance</li>
 *  <li>Remove the 'normal' plane group from the scene</li>
 *  <li>Add the single plane group to the scene</li>
 *  <li>Adapt the top navigation bar, i.e. add a return button and remove quantity selector</li>
 *  <li>Replot</li>
 * </ol>
 *
 * @param {number} simulationIndex - the index of the simulation to display in single mode
 */
function enterSingleMode(simulationIndex) {
    console.assert(simulationIndex >= 0);
    console.assert(simulationIndex < simulationManager.getSimulations().length);

    simulationManager.enableSingleMode(simulationIndex);
    scene.remove(simulationManager.getGroup());
    scene.add(simulationManager.getSingleGroup());
    $("#exitSingleMode").show();
    $("#quantityDropdownWrapper").hide();
    $("#simulationListButton").hide();

    plotTransferFunction(frame);
    plotStatic();
}

/**
 * The opposite of {@link module:RSI~enterSingleMode}.
 * In order to exit single mode, following steps will occur in this function:
 * <ol>
 *  <li>Disable single mode in the {@link SimulationManager} instance</li>
 *  <li>Remove the single plane group from the scene</li>
 *  <li>Add the 'normal' plane group to the scene</li>
 *  <li>Adapt the top navigation bar back to its original state</li>
 *  <li>Replot</li>
 * </ol>
 */
function exitSingleMode() {
    simulationManager.disableSingleMode();
    scene.remove(simulationManager.getSingleGroup());
    scene.add(simulationManager.getGroup());
    $("#exitSingleMode").hide();
    $("#quantityDropdownWrapper").show();
    $("#simulationListButton").show();

    plotStatic();
    loadFrame(frame, true);
}

/**
 * Once data is received, this method is called.
 * This function is responsible for decoding the received data using
 * {@link module:RSI~b64ToFloatArray} and store the result using  the
 * {@link SimulationManager} instance.
 *
 * @param {Object} data - A data object as described in {@link module:RSI~b64ToFloatArray}
 * @param {bool} initial - true if initial data, false otherwise
 */
function onData(data, initial) {
    /**
     * If initial data, store it in the simulationManager (instead of simulationManager.getTarget()
     * as below in the else branch) and create new Simulation instance using simulationManager
     */
    if (initial) {
        var realArray = b64ToFloatArray(data.real);
        var transferArray = b64ToFloatArray(data.transfer);
        var fourierArray = new Float32Array(data.fourier);

        simulationManager.realInit = realArray;
        simulationManager.transferInit = transferArray;
        simulationManager.fourierInit = fourierArray;

        simulationManager.createSimulation();
    } else {
        simulationManager.getTarget().addFrame(
            decodeArray(data.real),
            decodeArray(data.transfer),
            decodeArray(data.fourier)
        );
    }
}


/**
 * Render an image of the specified simulation at a resolution specified by size.
 *
 * @param {number} simulationIndex - the index (with respect to the simulations array stored
 * in the {@link SimulationManager} instance) of the simulation of which to render an image
 * @param {number} size - side length in pixels of the image to render
 *
 * @return {Object} Canvas DOM element
 */
function renderImage(simulationIndex, size) {
    var simulation = simulationManager.getSimulations()[simulationIndex];
    console.assert(simulation);
    console.assert(simulation.mesh);

    var size = simulationManager.PLANE_SIZE / 2;
    orthoCam = new THREE.OrthographicCamera(-size, size, size, -size, 0, 1000);
    orthoCam.position.copy(simulation.mesh.position);
    orthoCam.position.y = 500; // Lift camera above mesh
    orthoCam.lookAt(simulation.mesh.position);

    imgRenderer.setSize(size, size);
    imgRenderer.render(scene, orthoCam);

    // var image = new Image();
    // image.src = imgRenderer.domElement.toDataURL();

    return imgRenderer.domElement;
}

/**
 * @callback gifProgressCallback
 * @param {number} frame - number of currently rendered frame
 * @param {number} frameCount - total number of all frames (rendered + not yet rendered)
 */

/**
 * @callback gifFinishCallback
 * @param {blob} blob - blob of .gif
 */

/**
 * Renders a .gif of the specified simulation.
 *
 * @param {number} simulationIndex - the index (with respect to the simulations array stored
 * in the {@link SimulationManager} instance) of the simulation of which to render the .gif
 * @param {number} size - side length in pixels of the .gif to render
 * @param {number} fps - frames per second
 * @param {number} quality - quality of resulting .gif, corresponds to pixel sampling interval.
 *                           Lower values are better, but also produce bigger files.
 * @param {module:RSI~gifProgressCallback} progressCallback - progress callback
 * @param {module:RSI~gifFinishCallback} finishCallback - finish callback
 */
function renderGIF(simulationIndex, size, fps, quality, progressCallback, finishCallback) {
    var gif = new GIF({
        workers: 2,
        quality: quality,
        workerScript: "/static/js/gif.worker.js",
    });

    var previousFrame = frame;

    var frameCount = getLastFrame() + 1;
    var delay = 1000.0 / fps;

    for (var i = 0; i < frameCount; ++i) {
        loadFrame(i);
        gif.addFrame(renderImage(simulationIndex, size), {
            copy: true,
            delay: delay
        });
        if (progressCallback) {
            progressCallback(i, frameCount);
        }
    }
    progressCallback(frameCount, frameCount);

    if (finishCallback) {
        gif.on('finished', finishCallback);
    }

    gif.render();

    loadFrame(previousFrame);
}


/**
 * Gets called on close / loss of connection and displays a corresponding message to the user.
 */
function onClose() {
    closed = true;
    alert("Connection to server lost; Reload page to start a new session.");
}

/**
 * Initializes the quantity dropdown by assigning an event handler, which
 * updates the currently active quantity and reloads the current frame in order
 * for the change to take effect.
 */
function initQuantityDropdown() {
    var dropDown = $("#quantity-dropdown");
    dropDown.children().click(function(e) {
        dropDown.find(".active").removeClass("active");
        $(this).addClass("active");
        $("#quantity-dropdown-btn").html($(this).html());
        activeQuantity = $(this).attr("data-quantity");
        loadFrame(frame, true);

        e.preventDefault();
    });
}

/**
 * Initializes the dropup color map selector by assigning an event listener
 * which loads the appropriate texture using {@link SimulationManager#setColorMap}.
 */
function initColorMapDropdown() {
    var menu = $("#colormap-selector-menu");
    menu.find(".colormap-selector-item").click(function(e) {
        menu.find(".colormap-selector-item.active").removeClass("active");
        $(this).addClass("active");
        $("#colormap-selector-button") .css({
            "background": "url(" + $(this).find("img").attr("src") + ")",
            "background-size": "contain"
        });

        var texture = new THREE.Texture($(this).find("img").get(0));
        texture.needsUpdate = true;
        simulationManager.setColorMap(texture);

        e.preventDefault();
    });

    menu.find(":contains(default)").click();
}

/**
 * Initializes the checkbox which toggles between stopping the playback
 * at decoupling and today.
 */
function initRedshiftEndToggle() {
    $("#redshift-end-toggle").change(function(e) {
        stopAtDecoupling = !$(this).get(0).checked;
        console.log("Stop at decoupling? " + stopAtDecoupling);
        var lastFrame = getLastFrame();
        if (frame > lastFrame) {
            loadFrame(lastFrame);
        }
        updateTimeline(lastFrame);
    });
}

/**
 * Initializes the .gif creation dialog
 */
function initGifDialog() {
    log("Initializing .gif dialog");
    var modal = $("#gifExportModal");
    modal.on("show.bs.modal", function() {
        $(this).find("#gifExportProgressContainer").hide();
        $(this).find("#gifExportResultContainer").hide();
    });
    $("#gifExportCreateBtn").click(function() {
        onCreateGif(modal.data("simulationIndex"));
    });
}

/**
 * Initializes the image creation dialog
 */
function initImgDialog() {
    var modal = $("#imgExportModal");
    modal.on("show.bs.modal", function() {
        $(this).find("#imgExportResultContainer").hide();
    });
    $("#imgExportCreateBtn").click(function() {
        $("#imgExportResultContainer").show();

        // Hide grids if necessary
        var gridShown = parameterGui.config.showGrid;
        parameterGui.config.showGrid = $("#imgExportGrid").is(":checked");
        simulationManager.updateGrids();

        var canvas = renderImage(modal.data("simulationIndex"), $("#imgExportSize").val());
        $("#imgExportResultImage").attr("src", canvas.toDataURL());

        // Re-show grids if necessary
        parameterGui.config.showGrid = gridShown;
        simulationManager.updateGrids();
    });
}

/**
 * Callback function which is called from the {@link SimulationList} instance once
 * the user presses the 'Create .gif' button there.
 *
 * @param {number} simulationIndex - index of simulation to render as .gif
 */
function onCreateGif(simulationIndex) {
    log("Creating gif for simulation #" + simulationIndex);
    var size = $("#gifExportSize").val();
    var fps = $("#gifExportFPS").val();
    parameterGui.config.showGrid = $("#gifExportGrid").is(":checked");
    simulationManager.updateGrids();

    var quality = $("#gifExportQuality").val();

    $("#gifExportProgressContainer").show();
    $("#gifExportResultContainer").hide();

    renderGIF(simulationIndex, size, fps, quality, gifProgressCallback, gifFinishCallback);
}

/**
 * Implements {@link module:RSI~gifProgressCallback} and is called after every frame
 * rendered to the .gif.
 * Uses this data to update the progress bar in the .gif creation dialog.
 *
 * @param {number} frame - number of currently rendered frame
 * @param {number} frameCount - total number of all frames (rendered + not yet rendered)
 */
function gifProgressCallback(frame, frameCount) {
    var progressPercent = Math.floor(frame / frameCount * 100) + "%";
    $("#gifExportProgressBar").css({width: progressPercent});
    if (frame == frameCount) {
        $("#gifExportProgressBar").html("Finalizing .gif&hellip;");
    } else {
        $("#gifExportProgressBar").text("Rendering Frame " + frame + "/" + frameCount);
    }

}

/**
 * Implements {@link module:RSI~gifFinishCallback} and is called once the .gif has been rendered.
 * Proceeds to reset the progress bar in the .gif creation dialog and displays the resulting
 * .gif.
 *
 * @param {blob} blob - data blob of .gif
 */
function gifFinishCallback(blob) {
    $("#gifExportProgressBar") .css({width: 0})
    $("#gifExportProgressContainer").hide();
    $("#gifExportResultContainer").show();
    $("#gifExportResultImage").attr("src", URL.createObjectURL(blob));
    $("#gifExportResultMessage").alert();
}

/**
 * Toggles the visibility of the plot panel.
 */
function togglePlotWindow() {
    plotWindowVisible = !plotWindowVisible;
    $("#plotWindowWrapper").toggleClass("plotWindowWrapperVisible").toggleClass("plotWindowWrapperHidden");
    $("#plotWindowWrapper").find("div").toggleClass("plotWindowVisible").toggleClass("plotWindowHidden");
    $("#plotWindowToggleIcon").toggleClass("oi-caret-top").toggleClass("oi-caret-bottom");

    if (plotWindowVisible) {
        if (simulationManager.getActive().length == 0) {
            initClPlot();
            initTransferPlot();
        } else {
            plotStatic();
            plotTransferFunction(frame, true);
        }
    }
}

/**
 * Initializes the button used to collapse the plot panel.
 */
function initPlotCollapseButton() {
    $("#plotWindowToggle").click(function(e) {
        togglePlotWindow();
    });
}


/* SECTION: PLOT INITIALIZATION */

/**
 * Initialize transfer function plot.
 */
function initTransferPlot() {
    var initialOptions = Object.assign({}, transferPlotOptions);
    initialOptions.xaxis.min = 1e-4;
    initialOptions.xaxis.max = 10;
    transferFunctionPlot = $.plot($("#transferFunctionPlot"), [{shadowSize: 0, color: colors[0], data: []}], initialOptions);
}

/**
 * Initialize Cl's plot.
 */
function initClPlot() {
    var initialOptions = Object.assign({}, tClPlotOptions);
    initialOptions.xaxis.min = 2;
    initialOptions.xaxis.max = 2500;
    $.plot($("#tClPlot"), [{shadowSize: 0, color: colors[0], data: []}], initialOptions);
}

/**
 * Initialize matter spectrum plot.
 */
function initmPkPlot() {
    var initialOptions = Object.assign({}, mPkPlotOptions);
    initialOptions.xaxis.min = 1e-3;
    initialOptions.xaxis.max = 10;
    $.plot($("#mPkPlot"), [{shadowSize: 0, data: []}], initialOptions);
}
/* END OF PLOT INITIALIZATION */

/**
 * Callback to be called once cosmological parameters have been set
 * @callback cosmologicalParametersSetCallback
 * @param {module:paramGui~CosmoParams} cosmoParamsInstance
 * @param {Object} serializedMessage - message to be sent to the server
*/

/**
 * Callback to be called once initial condition has been set
 * @callback initialConditionSetCallback
 * @param {module:paramGui~ParamSet} paramSet
*/

/**
 * Called once user presses button to generate initial state.
 * Obtains the JSON representation and sends it to the server.
 * Since it isn't particularly enlightening to compare simulations (resulting from different
 * cosmological parameters) for different initial state seeds, this also cleans up
 * any existing simulations.
 *
 * @param {module:paramGui~ParamSet} paramSet - initial condition parameter set
 */
function onInitialConditionSet(paramSet) {
    if (!calculationRunning) {
        if (simulationManager.singleMode) {
            exitSingleMode();
        }
        log("Sending initial parameters to server.");
        var serializedParamSet = paramSet.serializeInitial();
        var paramString = JSON.stringify(serializedParamSet);
        log(paramString);
        calculationRunning = true;
        showActivityIndicator();
        connection.send(paramString);

        simulationManager.deleteAll();
        simuTable.clear();

        var realScale = serializedParamSet.params.xScale / simulationManager.GRID_DIVISIONS;
        $("#displayRealScale").text(round(realScale, 1) + " Mpc");

        receivedInitial = false;
    }
}

/**
 * Called once the user presses the button which sets the cosmological parameters.
 * <br>
 * Implements {@link module:RSI~cosmologicalParametersSetCallback}.
 * <br>
 * If an initial condition has been set, the following steps will occur:
 * First, this function checks using {@link SimulationManager~wasAlreadySimulated}
 * whether a simulation set for that particular set of cosmological parameters has already
 * been computed.
 * If so, a message indicating that will be displayed.
 * If not, the parameters (which have already been encoded in the appropriate structure
 * by {@link module:paramGui~CosmoParams}) will be sent to the server and the progress
 * modal will be displayed.
 * Also, a new simulation will be created and a copy of the parameters will be stored
 * in it.
 *
 * @param {module:paramGui~CosmoParams} cosmoParamsInstance
 * @param {Object} serializedMessage - message to be sent to the server
 */
function onCosmologicalParamsSet(cosmoParamsInstance, serializedMessage) {
    if (!initialSet) {
        alert("Initial condition has to be set!");
        return;
    }
    if (!calculationRunning) {
        log("Sending cosmological parameters to server.");

        var alreadySimulated = simulationManager.wasAlreadySimulated(cosmoParamsInstance.parameters);
        log("Already simulated? " + alreadySimulated);
        if (!alreadySimulated) {
            calculationRunning = true;
            showActivityIndicator();
            connection.send(JSON.stringify(serializedMessage));

            $("#simulationProgressModal").modal("show");
            $("#simulationProgressInfo").text("Running calculation... (this may take some time depending on the choice of parameters!)");
            $("#simulationProgressBar").css({width: "0%"});
            $("#simulationProgressBar").text("0%");

            if (!simulationManager.hasTarget()) {
                simulationManager.createSimulation();
            }
            /**
             * At this point, it is absolutely necessary to create a copy of the
             * parameters. If not, any change by the user in the parameters in the
             * control panel will be reflected by the parameters stored in the simulation
             * object, since both references refer to the same object.
             * This results in simulationManager.wasAlreadySimulated always returning
             * true and hence preventing any further simulation runs.
             */
            var simParams = Object.assign({}, cosmoParamsInstance.parameters);
            simulationManager.getTarget().setParameters(simParams);
        } else {
            $("#alreadySimulatedModal").modal("show");
        }
    }
}

/**
 * Loads the given frame.
 *
 * @param {number} f - frame to load
 * @param {bool} replot - if true, fully replot (which includes recalculating
 *                        data limits, axis labels, etc.). otherwise, just swap
 *                        out the data (faster, but this effect decreases with data
 *                        size, since then most of the CPU time will be spent on plotting
 *                        the curve anyways)
 */
function loadFrame(f, replot) {
    frame = f;
    simulationManager.loadFrame(activeQuantity, f);

    plotTransferFunctionIfVisible(f, replot);

    if (f > 0)
        $("#DisplayRedshift").text(round(redshift[f], REDSHIFT_ROUND));
    else if (f == 0)
        $("#DisplayRedshift").text(redshift[0]);

    var frameCount = simulationManager.getFrameCount();
    // var percentage = (frameCount > 1) ? frame / (frameCount - 1) : frame / frameCount;
    updateTimeline(getLastFrame());
}

/**
 * Plots the transfer function for the given frame if and only if the plot window
 * is visible.
 *
 * @param {number} f - frame number
 * @param {bool} replot - same as in {@link module:RSI~loadFrame}
 */
function plotTransferFunctionIfVisible(f, replot) {
    if (plotWindowVisible) {
        plotTransferFunction(f, replot);
    }
}

/**
 * Finds the minimum and maximum of <i>all</i> transfer functions
 * over <i>all</i> frames.
 *
 * @return {number[]} tuple containing minimum and maximum, in that order
 */
function findTransferFunctionMinMax() {
    var min = Infinity, max = -Infinity;
    var frameCount = simulationManager.getFrameCount();
    for (var _f = 0; _f < frameCount; ++_f) {
        var transferFunctions = simulationManager.getTransferData(activeQuantity, _f);
        for (var i = 0; i < transferFunctions.length; ++i) {
            var transferFunction = transferFunctions[i];
            var _min = Math.min.apply(null, transferFunction);
            var _max = Math.max.apply(null, transferFunction);

            min = Math.min(min, _min);
            max = Math.max(max, _max);
        }
    }

    return [min, max];
}

/**
 * Plot transfer function in non-single mode, i.e. one quantity but for
 * multiple simulations.
 *
 * @param {number} f - frame number
 * @param {bool} replot - same as in {@link module:RSI~loadFrame}
 */
function plotTransferFunctionMulti(f, replot) {
    // Transfer function
    var tData = simulationManager.getTransferData(activeQuantity, f);
    var activeIndices = simulationManager.getActive();
    var datas = [];
    for (var i = 0; i < activeIndices.length; ++i) {
        var color = colors[activeIndices[i] % colors.length];
        var data = zip(kRange, tData[i]);
        datas.push({
            shadowSize: 0,
            color: color,
            data: data,
        });
    }

    transferFunctionPlot = $.plot($("#transferFunctionPlot"), datas, transferPlotOptions);
}

/**
 * Plot transfer function in single mode.
 *
 * @param {number} f - frame number
 */
function plotTransferFunctionSingle(f) {
    var replot = true;
    var datas = [];
    // Transfer function
    for (var i = 0; i < QUANTITIES.length; ++i) {
        var tData = simulationManager.getTransferDataOfSingle(QUANTITIES[i], f);
        var color = colors[i % colors.length];
        var data = zip(kRange, tData);
        datas.push({
            shadowSize: 0,
            color: color,
            data: data,
        });
    }

    transferFunctionPlot = $.plot($("#transferFunctionPlot"), datas, transferPlotOptions);
    return;

    if (replot) {
    }
    else {
        transferFunctionPlot.setData(datas);
        transferFunctionPlot.draw();
    }
}

/**
 * Plot the transfer function, automatically in the correct mode (single vs. non-single).
 *
 * @param {number} f - frame number
 * @param {bool} replot - same as in {@link module:RSI~loadFrame}
 */
function plotTransferFunction(f, replot) {
    if (simulationManager.singleMode) {
        plotTransferFunctionSingle(f);
    } else {
        plotTransferFunctionMulti(f, replot);
    }
}

/*******************************************************************************
* Cl plotting
*******************************************************************************/

/**
 * Plot Cl spectrum in non-single mode.
 */
function plotClsMulti() {
    var plotDatas = [];

    simulationManager.forEachActive(function(index, simulation) {
        var dataItem = {
            shadowSize: 0,
            color: colors[index % colors.length],
            data: zip(simulation.Cls.l, simulation.Cls.tCl),
        }
        plotDatas.push(dataItem);
    });

    $.plot($("#tClPlot"), plotDatas, tClPlotOptions);
}

/**
 * Plot Cl spectrum in single mode.
 */
function plotClsSingle() {
    var tuple = simulationManager.getSingleSimulation();
    var singleIndex = tuple[0];
    var singleSimulation = tuple[1];

    var plotData = {
        shadowSize: 0,
        color: colors[singleIndex % colors.length],
        data: zip(singleSimulation.Cls.l, singleSimulation.Cls.tCl)
    };
    $.plot($("#tClPlot"), [plotData], tClPlotOptions);
}

/*******************************************************************************
* mPk plotting
*******************************************************************************/

/**
 * Plot matter spectrum in non-single mode.
 */
function plotmPkMulti() {
    var plotDatas = [];

    simulationManager.forEachActive(function(index, simulation) {
        var dataItem = {
            shadowSize: 0,
            color: colors[index % colors.length],
            data: zip(simulation.mPk.kh, simulation.mPk.Pkh),
        }
        plotDatas.push(dataItem);
    });

    $.plot($("#mPkPlot"), plotDatas, mPkPlotOptions);
}

/**
 * Plot matter spectrum in single mode.
 */
function plotmPkSingle() {
    var tuple = simulationManager.getSingleSimulation();
    var singleIndex = tuple[0];
    var singleSimulation = tuple[1];

    var plotData = {
        shadowSize: 0,
        color: colors[singleIndex % colors.length],
        data: zip(singleSimulation.mPk.kh, singleSimulation.mPk.Pkh)
    };
    $.plot($("#mPkPlot"), [plotData], mPkPlotOptions);
}

/*******************************************************************************
* Static plotting (i.e. Cl's and mPk)
*******************************************************************************/

/**
 * Plot the static plots (Cl's and matter spectrum) in the appropriate mode
 */
function plotStatic() {
    if (simulationManager.singleMode) {
        plotStaticSingle();
    } else {
        plotStaticMulti();
    }
}

/**
 * Plot the static plots (Cl's and matter spectrum) only if the plot panel
 * is visible.
 */
function plotStaticIfVisible() {
    if (plotWindowVisible)
        plotStatic();
}

/**
 * Plot Cl's and matter spectrum in non-single mode.
 */
function plotStaticMulti() {
    plotClsMulti();
    plotmPkMulti();
}

/**
 * Plot Cl's and matter spectrum in single mode.
 */
function plotStaticSingle() {
    plotClsSingle();
    plotmPkSingle();
}


/**
 * Analogon of the zip function in python. Takes two arrays [a1, a2, a3, ...] and
 * [b1, b2, b3, ...] and returns [[a1, b1], [a2, b2], [a3, b3], ...].
 * Required for plotting, as flot.js requires its data to be specified in a 'zipped'
 * format.
 *
 * @param {Object[]} a - first array
 * @param {Object[]} b - second array
 *
 * @return {Array.<Array.<Object>>}} the zipped array
 */
function zip(a, b) {
    var result = [];
    var limit = Math.min(a.length, b.length);
    for (var i = 0; i < limit; ++i) {
        result.push([a[i], b[i]]);
    }
    return result;
}

/*
 * Player Panel Callbacks
 */

/**
 * Implements {@link PlayerPanel~playPauseCallback}.
 */
function onPlayPause(playing) {
    animationRunning = playing;
    if (playing) {
        if (frame == getLastFrame()) {
            frame = 0;
        }
    } else {
    }
}

/**
 * Implements {@link PlayerPanel~scrubCallback}.
 */
function scrubCallback(progress) {
    var f = Math.floor((getLastFrame() + 1) * progress);
    f = Math.min(getLastFrame(), f);
    loadFrame(f);
}

/**
 * Implements {@link PlayerPanel~backFrameCallback}.
 */
function onBackFrame() {
    if (frame > 0 && receivedInitial) {
        loadFrame(frame - 1);
    }
}

/**
 * Implements {@link PlayerPanel~forwardFrameCallback}.
 */
function onForwardFrame() {
    if (frame < getLastFrame())
        loadFrame(frame + 1);
}

/**
 * Implements {@link PlayerPanel~toStartCallback}.
 */
function onToStart() {
    frame = 0;
    loadFrame(frame);
}

/**
 * Implements {@link PlayerPanel~toEndCallback}.
 */
function onToEnd() {
    if (receivedInitial) {
        loadFrame(getLastFrame());
    }
}

/**
 * Returns the number of the last frame to be displayed.
 * If {@link module:RSI~stopAtDecoupling} is true, return the last frame before
 * decoupling, otherwise the actual last frame stored.
 *
 * @return {number} number of last frame
 */
function getLastFrame() {
    var frameCount = simulationManager.getFrameCount();
    return lastFrame = (stopAtDecoupling) ? Math.min(decouplingFrame, frameCount - 1) : frameCount - 1;
}


/**
 * Initializes all components.
 */
function init() {
    hideActivityIndicator();

    window.addEventListener("mousemove", onMouseMove, false);

    parameterGui = new ParameterGui(onInitialConditionSet, onCosmologicalParamsSet);

    simulationManager = new SimulationManager(cmapTexture);
    playerPanel = new PlayerPanel(
        onPlayPause,
        onBackFrame,
        onForwardFrame,
        onToStart,
        onToEnd,
        function(e) {},
        scrubCallback
    );

    initQuantityDropdown();
    initColorMapDropdown();
    initGifDialog();
    initImgDialog();

    initTransferPlot();
    initClPlot();
    initmPkPlot();

    initPlotCollapseButton();

    initRedshiftEndToggle();

    ///////////
    // Scene //
    ///////////

    container = document.getElementById("container");

    scene = new THREE.Scene();
    scene.add(simulationManager.getGroup());
    scene.add(simulationManager.getActiveBoxes());
    // scene.add(simulationManager.getCollisionGroup());
    scene.background = new THREE.Color("#4a4a4a");

    var width = window.innerWidth;
    var height = window.innerHeight;
    var aspect = width / height; // view aspect ratio

    ////////////
    // Camera //
    ////////////

    camera = new THREE.PerspectiveCamera(FoV, aspect, 0.1, 10000);

    camera.position.y = 300;
    camera.position.z = 500;
    camera.position.x = -50;
    camera.lookAt(scene.position);


    //////////////
    // Renderer //
    //////////////
    renderer = new THREE.WebGLRenderer({antialias: true});
    renderer.setClearColor(0x53504d);
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.autoClear = false;

    container.innerHTML = "";
    container.appendChild(renderer.domElement);

    //////////////
    // Controls //
    //////////////
    controls = new THREE.OrbitControls(camera, renderer.domElement);

    controls.rotateSpeed = 1.0;
    controls.zoomSpeed = 3;
    controls.panSpeed = 0.8;

    controls.noZoom = false;
    controls.noPan = false;

    controls.maxPolarAngle = Math.PI/2;
    controls.minDistance = 100;
    controls.maxDistance = 10000;

    controls.staticMoving = true;
    controls.dynamicDampingFactor = 0.3;

    controls.keys = [65, 83, 68];

    controls.addEventListener('change', render);

    /******************************************************************
     * SIMULATION TABLE
     ******************************************************************/

    var deleteCallback = function(index) {
        simulationManager.delete(index);
        simuTable.clear();
        simuTable.populate(simulationManager.getSimulations(), simulationManager.getActive());
        plotStaticIfVisible();
        loadFrame(frame, true);
    };

    var activateCallback = function(index) {
        simulationManager.activate(index);
        loadFrame(frame, true);
        plotStaticIfVisible();
    };

    var deactivateCallback =  function(index) {
        simulationManager.deactivate(index);
        loadFrame(frame, true);
        plotStaticIfVisible();
    };

    var imageCallback = function(index) {
        if (simulationManager.isActive(index)) {
            $("#simulationListModal").modal("hide");
            $("#imgExportModal").data("simulationIndex", index).modal("show");
        }
        else {
            alert("To create an image, the selected simulation needs to be active.");
        }
    };

    var gifCallback = function(index) {
        if (simulationManager.isActive(index)) {
            $("#simulationListModal").modal("hide");
            $("#gifExportModal").data("simulationIndex", index).modal("show");
        } else {
            alert("To create a .gif, the selected simulation needs to be active.");
        }
        // renderGIF(index);
    };

    var singleModeCallback = function(index) {
        $("#simulationListModal").modal("hide");
        enterSingleMode(index);
    };

    var focusCameraCallback = function(index) {
        if (simulationManager.isActive(index)) {
            controls.center.copy(simulationManager.getSimulations()[index].mesh.position);
            controls.rotateUp(controls.getPolarAngle());
            controls.rotateLeft(controls.getAzimuthalAngle());
            $("#simulationListModal").modal("hide");
        } else {
            alert("To focus the camera on the selected simulation, it needs to be active.");
        }
    };

    var loadParamCallback = function(index) {
        var parameters = simulationManager.getSimulations()[index].params;
        parameterGui.loadCosmoParameters(parameters);
    };

    simuTable = new SimuTable(deleteCallback, activateCallback, deactivateCallback,
        gifCallback, imageCallback, singleModeCallback, focusCameraCallback,
        loadParamCallback);

    simuTable.createHeader();
    /******************************************************************
     * END OF SIMULATION TABLE
     ******************************************************************/

    // stats = new Stats();
    // stats.domElement.style.position = 'absolute';
    // stats.domElement.style.bottom = '0px';
    // stats.domElement.style.zIndex = 100;
    // container.appendChild(stats.domElement);

    window.addEventListener('resize', onWindowResize, false);
}

/**
 * Shows the activity indicator.
 */
function showActivityIndicator() {
    $("#activityIndicator").removeClass("activityIndicator-inactive").addClass("activityIndicator-active");
}

/**
 * Hides the activity indicator.
 */
function hideActivityIndicator() {
    $("#activityIndicator").addClass("activityIndicator-inactive").removeClass("activityIndicator-active");
}

/**
 * Called when the window is resized to update the renderer (which requires
 * the window size to render) and the camera (which requires the aspect ratio
 * of the window to construct its projection matrix).
 */
function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();

    renderer.setSize(window.innerWidth, window.innerHeight);

    // Re-read the width of the time line once window refreshes
    playerPanel.timeline.refresh();
}

/**
 * This function is called every time an update (such as an advance in frame, etc.)
 * is needed.
 * It handles the logic required to play back the animation, such as advancing,
 * repeating (if repeat is enabled), pausing and continuing the playback, etc.
 * The rate at which this function is called is determined by the FPS.
 */
function step() {
    var totalFrames = simulationManager.getFrameCount();
    var lastFrame = getLastFrame();
    if(animationRunning) {
        if(frame < lastFrame) {
            loadFrame(++frame);
        }
        else if (frame == lastFrame) {
            if (playerPanel.repeat) {
                frame = 1;
                loadFrame(frame);
            } else
                playerPanel.pause();
        }
    }
    setTimeout(step, 1000 / parameterGui.config.animationSpeed);
}

/**
 * Positions the 'scrubber' of the playback timeline at {@link module:RSI~frame}
 * given the last frame position.
 *
 * @param {number} lastFrame - index of last frame
 */
function updateTimeline(lastFrame) {
    playerPanel.timeline.scrubTo(frame / lastFrame);
}


/**
 * Function that is called on every frame of the main loop.
 */
function animate() {
    requestAnimationFrame(animate);
    controls.update();

    render();
}

/**
 * Rendering function, called on every frame of the main loop via
 * {@link module:RSI~animate}.
 * On each frame, this renders the scene and updates the raycaster.
 * Then determines over which plane (if at all) the mouse hovers
 * and instructs the {@link SimulationManager} instance to highlight it.
 */
function render() {
    raycaster.setFromCamera(mouse, camera);
    simulationManager.mousePick(raycaster);

    // renderer.clear();
    renderer.render(scene, camera);
    // renderer.clearDepth();
    // stats.update();
}
