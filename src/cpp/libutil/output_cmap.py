import numpy as np
import matplotlib.pyplot as plt

def colormap_to_colorizer(colormap_name, n=256):
    """
    Generate a Colorizer-style string for a given matplotlib colormap.
    
    Returns a string like: "0:000000,0.003921:000003,...,1:ffffff"
    """
    cmap = plt.get_cmap(colormap_name, n)
    colors = (cmap(np.linspace(0, 1, n))[:, :3] * 255).astype(int)
    
    entries = []
    for i, (r, g, b) in enumerate(colors):
        pos = i / (n - 1)  # normalize position to [0,1]
        hexcolor = f"{r:02X}{g:02X}{b:02X}"
        entries.append(f"{pos}:{hexcolor}")
    
    return ",".join(entries)

# Example usage
if __name__ == "__main__":
    for cmap_name in ['afmhot']:
        colorizer_string = colormap_to_colorizer(cmap_name)
        print(f'Colorizer Colorizer::{cmap_name}("{colorizer_string}");\n')

