# visualizacion_ultimo_snapshot.py
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

NX, NY = 600, 520
DX, DY = 50.0, 50.0
DAM_X_POS = 314

zb_path = '/home/nicolas/UrraDamBreachSimulation/data/urra_specific_topography.bin'
zb = np.fromfile(zb_path, dtype=np.float64).reshape(NY, NX)

snapshot_dir = '/home/nicolas/UrraDamBreachSimulation/snapshots'
files = sorted(glob.glob(os.path.join(snapshot_dir, 'h_*.bin')))

if not files:
    print("ERROR: No hay snapshots")
    print("Ejecuta primero la simulación")
    exit(1)

ultimo_file = files[-1]
h = np.fromfile(ultimo_file, dtype=np.float64).reshape(NY, NX)

print(f"Último snapshot: {os.path.basename(ultimo_file)}")
print(f"Profundidad máxima: {h.max():.2f} m")
print(f"Profundidad media: {h.mean():.2f} m")

fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111, projection='3d')

X = np.arange(NX) * DX
Y = np.arange(NY) * DY
X, Y = np.meshgrid(X, Y)

surf = ax.plot_surface(X, Y, zb, cmap='terrain', alpha=0.7, linewidth=0, antialiased=True)

h_filtrado = np.where(h > 0.5, h, 0)
water_surface = zb + h_filtrado
water_masked = np.where(h_filtrado > 0, water_surface, np.nan)

water = ax.plot_surface(X, Y, water_masked, cmap='Blues', alpha=0.8, linewidth=0, antialiased=True)

presa_x = DAM_X_POS * DX
ax.plot([presa_x, presa_x], [0, NY*DY], [0, 200], color='red', linewidth=3, label='Presa')

ax.set_xlabel('Distancia (m)')
ax.set_ylabel('Ancho (m)')
ax.set_zlabel('Elevación (m)')
ax.set_title(f'Estado final - Rotura de Presa Urrá\nTiempo simulado: último snapshot')
ax.view_init(elev=25, azim=-60)
ax.legend()

plt.savefig('estado_final_3d.png', dpi=200)
plt.show()

print("Visualización guardada como 'estado_final_3d.png'")