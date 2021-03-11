/**
 * This function will be called when the user presses the play/pause button.
 * @callback PlayerPanel~playPauseCallback
 * @param {bool} playing - true if movie should be played, false if movie should be paused.
 */

/**
 * This function will be called when the user presses the "back one frame" button.
 * @callback PlayerPanel~backFrameCallback
 */

/**
 * This function will be called when the user presses the "forward one frame" button.
 * @callback PlayerPanel~forwardFrameCallback
 */

/**
 * This function will be called when the user presses the "Go to Start" button.
 * @callback PlayerPanel~toStartCallback
 */

/**
 * This function will be called when the user presses the "Go to End" button.
 * @callback PlayerPanel~toEndCallback
 */

/**
 * This function will be called when the user presses the play/pause button.
 * @callback PlayerPanel~repeatCallback
 * @param {bool} repeating - true if movie should be played in repeat, else false.
 */

/**
 * This function will be called when the user uses the 'scrubber' to move to a specific frame.
 * @callback PlayerPanel~scrubCallback
 * @param {number} progress - number between 0.0 and 1.0 (inclusive).
 */

/**
 * Represents the player panel.
 * @constructor
 * @param {PlayerPanel~playPauseCallback}
 * @param {PlayerPanel~backFrameCallback}
 * @param {PlayerPanel~forwardFrameCallback}
 * @param {PlayerPanel~toStartCallback}
 * @param {PlayerPanel~toEndCallback}
 * @param {PlayerPanel~repeatCallback}
 * @param {PlayerPanel~scrubCallback}
 */
function PlayerPanel(playPauseCallback, backFrameCallback, forwardFrameCallback, toStartCallback, toEndCallback, repeatCallback,
    scrubCallback) {
    this.playPauseCallback = playPauseCallback;
    this.backFrameCallback = backFrameCallback;
    this.forwardFrameCallback = forwardFrameCallback;
    this.toStartCallback = toStartCallback;
    this.toEndCallback = toEndCallback;
    this.repeatCallback = repeatCallback;

    this.playPauseButton = document.getElementById("playerPlayPauseBtn");
    this.backFrameButton = document.getElementById("playerBackFrameBtn");
    this.forwardFrameButton = document.getElementById("playerForwardFrameBtn");
    this.toStartButton = document.getElementById("playerToStartBtn");
    this.toEndButton = document.getElementById("playerToEndBtn");
    this.repeatButton = document.getElementById("playerRepeatBtn");

    this.timeline = new Timeline("timeline", "scrubber", scrubCallback);
    this.timeline.scrubTo(0);

    this.playing = false;
    this.repeat = false;

    this.connectCallbacks();
}

PlayerPanel.prototype.connectCallbacks = function() {
    var self = this;
    this.playPauseButton.onclick = function() {
        // Additionally, toggle button
        var buttonSpan = self.playPauseButton.firstElementChild;
        if (self.playing) {
            buttonSpan.classList.remove("oi-media-pause");
            buttonSpan.classList.add("oi-media-play");
        } else {
            buttonSpan.classList.remove("oi-media-play");
            buttonSpan.classList.add("oi-media-pause");
        }
        self.playing = !self.playing;
        self.playPauseCallback(self.playing);
    };
    this.backFrameButton.onclick = this.backFrameCallback;
    this.forwardFrameButton.onclick = this.forwardFrameCallback;
    this.toStartButton.onclick = this.toStartCallback;
    this.toEndButton.onclick = this.toEndCallback;

    this.repeatButton.onclick = function() {
        if (self.repeat) {
            self.repeatButton.classList.remove("btn-success");
            self.repeatButton.classList.add("btn-outline-success");
        } else {
            self.repeatButton.classList.remove("btn-outline-success");
            self.repeatButton.classList.add("btn-success");
        }
        self.repeat = !self.repeat;
        if (self.repeatCallback) {
            self.repeatCallback(self.repeat);
        }
    };
};

PlayerPanel.prototype.pause = function() {
    if (this.playing) {
        this.playPauseButton.onclick();
    }
};

PlayerPanel.prototype.play = function() {
    if (!this.playing) {
        this.playPauseButton.onclick();
    }
};




/* TIMELINE */
function Timeline(timelineId, scrubberId, scrubCallback) {
    this.timeline = document.getElementById(timelineId);
    this.scrubber = document.getElementById(scrubberId);

    this.scrubCallback = scrubCallback;

    this.timelineWidth = this.timeline.offsetWidth;
    this.scrubberRadius = this.scrubber.offsetWidth / 2;


    // Center
    var dy = this.scrubber.offsetHeight / 2 - this.timeline.offsetHeight / 2;
    this.scrubber.style.marginTop = "-" + Math.floor(dy) + "px";

    // Event handling
    this.mouseDown = false;
    var self = this;
    this.scrubber.addEventListener("mousedown", function(e) {
        self.mouseDown = true;
        e.stopPropagation();
    });

    this.timeline.addEventListener("mousedown", function(e) {
        self.handleClick(e);
    });

    document.addEventListener("mouseup", function() {
        self.mouseDown = false;
    });

    document.addEventListener("mousemove", function(e) {
        if (self.mouseDown) {
            self.handleClick(e);
        }
    });
}

Timeline.prototype.refresh = function() {
    this.timelineWidth = this.timeline.offsetWidth;
};

Timeline.prototype.handleClick = function(e) {
    var leftBoundary = this.timeline.getBoundingClientRect().left;
    var rightBoundary = leftBoundary + this.timelineWidth;
    if (e.clientX >= leftBoundary && e.clientX <= rightBoundary) {
        var percentage = (e.clientX - leftBoundary) / this.timelineWidth;
        this.scrubTo(percentage);
        this.scrubCallback(percentage);
    }
};

Timeline.prototype.scrubTo = function(percentage) {
    var pos = percentage * this.timelineWidth - this.scrubberRadius;
    this.scrubber.style.marginLeft = pos + "px";

    this.timeline.style.background = this.generateGradientCSS(percentage);
};

Timeline.prototype.generateGradientCSS = function(percent) {
    percent *= 100;
    // var result = "linear-gradient(to right, rgba(255, 0, 0, 0.75) 0%, rgba(255, 0, 0, 0.75) " + percent + "%, #222 " + percent + "%)";
    var result = "linear-gradient(to right, transparent 0%, transparent " + percent + "%, #222 " + percent + "%)";
    result = result + ", repeating-linear-gradient(45deg, hsl(134, 61%, 35%) 0%, hsl(134, 61%, 35%) 2%, hsl(134, 61%, 41%) 2%, hsl(134, 61%, 41%) 4%)";
    return result;
};
