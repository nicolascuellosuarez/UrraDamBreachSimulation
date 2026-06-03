import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob('/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'))

h0 = np.fromfile(files[0], dtype=np.float64).reshape(520, 600)
hf = np.fromfile(files[-1], dtype=np.float64).reshape(520, 600)

delta_h = hf - h0

print("Percentiles Δh:")
for p in [0, 1, 5, 25, 50, 75, 95, 99, 100]:
    print(f"{p:3d}%: {np.percentile(delta_h, p):8.3f} m")

vmin = np.percentile(delta_h, 5)
vmax = np.percentile(delta_h, 95)

plt.figure(figsize=(15, 6))

plt.imshow(
    delta_h,
    cmap='RdBu_r',
    vmin=vmin,
    vmax=vmax
)

plt.colorbar(label='Δh (m)')
plt.axvline(x=314, color='black', linestyle='--', linewidth=2, label='Presa')

plt.title('Cambio de profundidad (Final - Inicial)\nEscala basada en percentiles 5%-95%')
plt.xlabel('Distancia (celdas)')
plt.ylabel('Ancho (celdas)')
plt.legend()

plt.tight_layout()
plt.savefig('delta_h_percentiles.png', dpi=300)
plt.show()

print(f"\nΔh mínimo real: {delta_h.min():.3f} m")
print(f"Δh máximo real: {delta_h.max():.3f} m")
print(f"Δh promedio:    {delta_h.mean():.3f} m")