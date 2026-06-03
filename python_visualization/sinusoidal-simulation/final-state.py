import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob('/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'))
hf = np.fromfile(files[-1], dtype=np.float64).reshape(520, 600)

plt.figure(figsize=(15,6))

plt.imshow(
    hf,
    cmap='Blues',
    vmin=np.min(hf),
    vmax=np.percentile(hf,99)
)

plt.colorbar(label='Profundidad (m)')
plt.axvline(x=314, color='red', linestyle='--', label='Presa')
plt.title('Estado final')
plt.legend()

plt.savefig('real_final_state.png', dpi=200)
plt.show()