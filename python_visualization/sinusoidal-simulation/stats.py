import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob('/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'))

h0 = np.fromfile(files[0], dtype=np.float64).reshape(520, 600)
hf = np.fromfile(files[-1], dtype=np.float64).reshape(520, 600)

print(f"Inicial: min={h0.min():.3f}, max={h0.max():.3f}, mean={h0.mean():.3f}")
print(f"Final:   min={hf.min():.3f}, max={hf.max():.3f}, mean={hf.mean():.3f}")