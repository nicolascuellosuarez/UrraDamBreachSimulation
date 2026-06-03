import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob('/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'))

if not files:
    print("ERROR: No hay snapshots")
    exit(1)

h_final = np.fromfile(files[-1], dtype=np.float64).reshape(520, 600)

plt.figure(figsize=(15, 6))

plt.imshow(
    h_final,
    cmap='Blues',
    vmin=0,
    vmax=np.percentile(h_final, 95)
)

plt.colorbar(label='Profundidad (m)')
plt.axvline(x=314, color='red', linestyle='--', linewidth=2, label='Presa')
plt.title('Estado final - Rotura de Presa Urrá')
plt.xlabel('Distancia (celdas)')
plt.ylabel('Ancho (celdas)')
plt.legend()

plt.text(0.02, 0.95, f'Prof. máxima: {h_final.max():.2f} m\nProf. media: {h_final.mean():.2f} m',
         transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.8))

plt.savefig('2d_final_state.png', dpi=200)
plt.show()

print(f"Archivos procesados: {len(files)}")
print(f"Profundidad máxima: {h_final.max():.2f} m")
print(f"Profundidad media: {h_final.mean():.2f} m")