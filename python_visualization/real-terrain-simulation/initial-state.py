import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob('/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'))
h0 = np.fromfile(files[0], dtype=np.float64).reshape(520, 600)

plt.figure(figsize=(15,6))

plt.imshow(
    h0,
    cmap='Blues',
    vmin=np.min(h0),
    vmax=np.max(h0)
)

plt.colorbar(label='Profundidad (m)')
plt.axvline(x=314, color='red', linestyle='--', label='Presa')
plt.title('Estado inicial')
plt.legend()

plt.savefig('initial_visualizer.png', dpi=200)
plt.show()