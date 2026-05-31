import rasterio
import numpy as np
from scipy.ndimage import zoom

with rasterio.open('/home/nicolas/UrraDamBreachSimulation/data/urra_only.tif') as src:
    zb = src.read(1).astype(np.float64)

print(f"Tamaño original: {zb.shape}")

nx, ny = 600, 520
factor_x = nx / zb.shape[1]
factor_y = ny / zb.shape[0]
zb_redim = zoom(zb, (factor_y, factor_x), order=1)

zb_redim.tofile('urra_specific_topography.bin')
print(f"Guardado: {zb_redim.shape} -> data/urra_specific_topography.bin")