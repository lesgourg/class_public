import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import os

OUTPUT_DIR = os.path.join("static", "images", "colormaps")
WIDTH = 512

def create_image(cmap, width):
    values = np.linspace(0, 1, width)
    colors = cmap(values).reshape((1, width, 4))
    image = Image.fromarray(np.uint8(255 * colors))
    return image

cmap_names = {}
cmap_names['Uniform'] = [
            'viridis', 'plasma', 'inferno', 'magma']
cmap_names['Diverging'] = [
            'seismic', 'RdYlBu', 'Spectral'
            ]
cmap_names['Miscellaneous'] = ['jet']

if __name__ == "__main__":
    for category in cmap_names:
        category_dir = os.path.join(OUTPUT_DIR, category)
        if not os.path.exists(category_dir):
            os.mkdir(category_dir)
        for name in cmap_names[category]:
            result = create_image(plt.get_cmap(name), width=WIDTH)
            output_path = os.path.join(category_dir, "{}.png".format(name))
            print(output_path)
            result.save(output_path)

