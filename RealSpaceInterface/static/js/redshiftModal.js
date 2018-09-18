/**
 * @callback RedshiftModal~saveCallback
 * @param {Object} config - redshift configuration;
 *      for documentation on format, see {@link DEFAULT_REDSHIFTS}
 */

/**
 * Redshift Modal that allows user control over redshift sampling.
 * @typedef {Object} RedshiftModal
 *
 * @constructor
 * @param {Object} initial state, probably {@link DEFAULT_REDSHIFTS}
 * @param {RedshiftModal~saveCallback} callback - callback to be called once user saves
 */
function RedshiftModal(initial, callback) {
    this.modal = $("#redshift-modal");
    this.grid = $("#redshift-modal-grid");
    this.initial = initial;

    this.init(initial);
    this.callback = callback;
}

/**
 * Initialize the modal with the initial configuration passed to the {@link RedshiftModal} constructor.
 *
 * @param {Object} initial - initial configuration
 */
RedshiftModal.prototype.init = function(initial) {
    this.load(initial);

    // Deactivate first delete button
    this.grid.find(".redshift-modal-btn-delete").first().prop("disabled", true);

    var self  = this;
    this.modal.find("#redshift-modal-save").click(function(e) {
        self.clearAlerts();
        if (self.validate()) {
            self.callback(self.serialize());
            self.hide();
        }
    });

    this.modal.find("#redshift-modal-cancel").click(function(e) {
        self.clear();
        self.load(self.initial);
    });
};

/**
 * Removes all entries
 */
RedshiftModal.prototype.clear = function() {
    this.grid.children().remove();
};

/**
 * Show the modal
 */
RedshiftModal.prototype.show = function() {
    this.modal.modal("show");
};

/**
 * Hide the modal
 */
RedshiftModal.prototype.hide = function() {
    this.modal.modal("hide");
};

/**
 * Create an empty row template
 * @return {Object} row
 */
RedshiftModal.prototype.createRow = function() {
    return $("<div></div").addClass("row").addClass("mb-2");
};

/**
 * Creates a 'left row' by reading the corresponding template from RedshiftModal.html.
 * A 'left row' contains an input field for a z value.
 *
 * @return {Object} row
 */
RedshiftModal.prototype.createLeftRow = function() {
    return this.createRow().append($("#redshift-modal-left").html());
};

/**
 * Creates a 'right row' by reading the corresponding template from RedshiftModal.html.
 * A 'right row' contains an input field for the number of z samples, a dropdown for the
 * mode (either log or linear) and the two add/delete buttons.
 *
 * @return {Object} row
 */
RedshiftModal.prototype.createRightRow = function() {
    var row = this.createRow().append($("#redshift-modal-right").html());
    var self = this;
    row.find(".redshift-modal-btn-add").click(function() {
        var rowIndex = ($(this).parent().parent().parent().index() - 1) / 2;
        self.addRow(rowIndex);
    });
    row.find(".redshift-modal-btn-delete").click(function() {
        var rowIndex = $(this).parent().parent().parent().index();
        self.deleteRow(rowIndex);
    });
    return row;
};

/**
 * Deletes a row.
 *
 * @param {number} index - row number
 */
RedshiftModal.prototype.deleteRow = function(index) {
    // this.grid.find(".row").eq(index).remove();
    this.grid.find(".row").eq(index - 1).remove();
    this.grid.find(".row").eq(index - 1).remove();
};

/**
 * Adds a row after the given index.
 * Note that this index is interpreted as index of groups consisting
 * of z input row and sample input row, i.e. the absolute index in terms
 * of rows is 2 * index + 1.
 *
 * @param {number} after - the group index after which to add a row
 */
RedshiftModal.prototype.addRow = function(after) {
    this.grid.find(".row").eq(2 * after + 1).after(
        this.createLeftRow(),
        this.createRightRow()
    );
};

/**
 * Display a message in the modal. Used for signaling invalid input values.
 *
 * @param {string} message - message to display
 */
RedshiftModal.prototype.alert = function(message) {
    var element = $("<div></div>").addClass("alert").addClass("alert-warning").text(message);
    this.modal.find(".modal-body").append(element);
};

/**
 * Clear all alerts.
 */
RedshiftModal.prototype.clearAlerts = function() {
    this.modal.find(".alert").remove();
};

/**
 * Validate input fields and display warnings if input fields
 * contain illegal values.
 *
 * @return {bool} true if all inputs are valid; false otherwise.
 */
RedshiftModal.prototype.validate = function() {
    var self = this;
    var return_ = false;

    var values = this.modal.find(".redshift-modal-z").map(function() {
        if (isNaN(this.value)) {
            self.alert("Invalid value encountered: " + this.value);
            return_ = true;
        }
        var value = parseFloat(this.value);
        return value;
    }).get();

    var points = this.modal.find(".redshift-modal-samples").map(function() {
        if (isNaN(this.value)) {
            self.alert("Invalid point count encountered: " + this.value);
            return_ = true;
        }
        var value = parseInt(this.value);
        if (value <= 0) {
            self.alert("Invalid point count encountered: " + value);
            return_ = true;
        }
        return value;
    }).get();

    if (return_) {
        return false;
    }

    for (var i = 0; i < values.length - 1; ++i) {
        if (values[i] <= values[i + 1]) {
            self.alert("Redshifts must be decreasing!");
            return false;
        }
    }

    return true;
};

/**
 * Constructs an object representing the user's choice of z values.
 * For an example of the format, see {@link DEFAULT_REDSHIFTS}.
 */
RedshiftModal.prototype.serialize = function() {
    var z = this.modal.find(".redshift-modal-z").map(function() {
        return parseFloat(this.value);
    });

    var points = this.modal.find(".redshift-modal-samples").map(function() {
        return parseInt(this.value);
    });

    var logarithmic = this.modal.find(".redshift-modal-spacing").map(function() {
        return this.value == 0;
    });

    var result = [];
    for (var i = 0; i < z.length - 1; ++i) {
        result.push({
            from: z[i],
            to: z[i + 1],
            points: points[i],
            log: logarithmic[i],
        });
    }
    return result;
};

/**
 * The 'inverse' operation of {@link RedshiftModal#serialize}.
 * Takes an object of the same format and populates the modal
 * accordingly. Note that this does NOT clear the modal but merely
 * appends to the end of it.
 * To clear, call {@link RedshiftModal#clear}.
 *
 * @param {Object[]} data - the data to populate the modal with
 */
RedshiftModal.prototype.load = function(data) {
    for (var i = 0; i < data.length; ++i) {
        if (i == 0) {
            var lr = this.createLeftRow();
            lr.find(".redshift-modal-z").val(data[i].from);
            this.grid.append(lr);
        }

        var rr = this.createRightRow();
        rr.find(".redshift-modal-samples").val(data[i].points);
        this.grid.append(rr);

        var lr2 = this.createLeftRow();
        lr2.find(".redshift-modal-z").val(data[i].to);
        this.grid.append(lr2);
    }
};
