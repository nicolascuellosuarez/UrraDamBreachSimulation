import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob(
    '/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'
))

print(f"Snapshots encontrados: {len(files)}")
print(f"Primer snapshot: {files[0]}")

h = np.fromfile(files[0], dtype=np.float64).reshape(520, 600)

print(f"min={h.min():.3f}, max={h.max():.3f}, mean={h.mean():.3f}")

plt.figure(figsize=(15,6))

im = plt.imshow(
    h,
    cmap='Blues',
    aspect='auto'
)

plt.colorbar(im, label='Profundidad (m)')
plt.axvline(x=314, color='red', linestyle='--', linewidth=2, label='Presa')

plt.title('Primer snapshot de la simulación')
plt.xlabel('Distancia (celdas)')
plt.ylabel('Ancho (celdas)')
plt.legend()

plt.tight_layout()
plt.savefig('primer_snapshot.png', dpi=200)
plt.show()