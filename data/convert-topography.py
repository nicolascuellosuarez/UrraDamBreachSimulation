import rasterio
import numpy as np
from scipy.ndimage import zoom

with rasterio.open('/home/nicolas/UrraDamBreachSimulation/data/n08_w077_1arc_v3.tif') as src:
    zb_original = src.read(1).astype(np.float64)
    print(f"Tamaño original: {zb_original.shape}")

nx, ny = 600, 520
factor_x = nx / zb_original.shape[1]
factor_y = ny / zb_original.shape[0]
zb_redim = zoom(zb_original, (factor_y, factor_x), order=1)

zb_redim.tofile('urra_topography.bin')
print(f"Guardado: {zb_redim.shape} -> urra_topography.bin")
