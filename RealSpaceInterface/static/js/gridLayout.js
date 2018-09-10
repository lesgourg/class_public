/**
 * @class GridLayout
 *
 * @property {number} count - Maximum number of items in this layout
 * @property {number} margin - space between planes
 * @property {number} size - size of each plane
 */
function GridLayout(count, size, margin) {
    this.count = count;
    this.margin = margin;
    this.size = size;

    this.recalculate();
}

/**
 * Recalculate the number of rows of the layout.
 *
 * Must be called after each change of either columns or count.
 */
GridLayout.prototype.recalculate = function() {
    this.columns = Math.ceil(Math.sqrt(this.count));
    this.rows = Math.ceil(this.count / this.columns);
}

/**
 * Get the position of the item defined by index.
 */
GridLayout.prototype.getPosition = function(index) {
    if (index >= this.count) {
        throw "GridLayout error: index higher than maximum number of items in layout";
    }

    row = Math.floor(index / this.columns);
    column = index % this.columns;

    return [row, column];
}

/**
 * Get the actual position (in correct units of size) of
 * the item defined by index.
 */
GridLayout.prototype.getWorldPosition = function(index) {
    var pos = this.getPosition(index);

    var row = pos[0];
    var col = pos[1];

    var x = 0, y = 0;
    var dx = Math.floor(this.columns / 2);
    var dy = Math.floor(this.rows / 2);
    if (this.columns % 2 == 0)
        dx -= 1 / 2;
    if (this.rows % 2 == 0)
        dy -= 1 / 2;
    x = (col - dx) * (this.size + this.margin);
    y = (row - dy) * (this.size + this.margin);

    return [x, y];
}
