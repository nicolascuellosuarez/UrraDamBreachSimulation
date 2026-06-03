# video_3d_simple.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import glob
import os

NX, NY = 600, 520
DX, DY = 50.0, 50.0
dam_x = 314
UMBRAL_AGUA = 0.5
fps = 10
skip = 4 

x = np.arange(0, NX, skip) * DX / 1000.0  # km
y = np.arange(0, NY, skip) * DY / 1000.0  # km
X, Y = np.meshgrid(x, y)
nx_red, ny_red = len(x), len(y)

zb_path = '/home/nicolas/UrraDamBreachSimulation/data/urra_specific_topography.bin'
zb_full = np.fromfile(zb_path, dtype=np.float64).reshape(NY, NX)
zb_red = zb_full[::skip, ::skip]

snapshot_dir = '/home/nicolas/UrraDamBreachSimulation/snapshots'
files = sorted(glob.glob(os.path.join(snapshot_dir, 'h_*.bin')))
print(f"Encontrados {len(files)} snapshots")

if not files:
    print("ERROR: No hay snapshots")
    exit(1)

h_frames = []
for f in files:
    h_full = np.fromfile(f, dtype=np.float64).reshape(NY, NX)
    h_red = h_full[::skip, ::skip]
    h_frames.append(h_red)

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

terreno = ax.plot_surface(X, Y, zb_red, cmap='terrain', alpha=0.6, linewidth=0, antialiased=True)

h0 = np.where(h_frames[0] > UMBRAL_AGUA, h_frames[0], 0)
water_surface = zb_red + h0
agua = ax.plot_surface(X, Y, water_surface, cmap='Blues', alpha=0.8, linewidth=0, antialiased=True)

# Línea de la presa
presa_x_km = dam_x * DX / 1000.0
ax.plot([presa_x_km, presa_x_km], [0, max(y)], [0, 150], color='red', linewidth=3, label='Presa')

ax.set_xlabel('Distancia (km)')
ax.set_ylabel('Ancho (km)')
ax.set_zlabel('Elevación (m)')
ax.set_title('Simulación Rotura de Presa Urrá')
ax.view_init(elev=25, azim=-60)
ax.legend()

time_text = ax.text2D(0.02, 0.95, '', transform=ax.transAxes, color='white',
                      bbox=dict(facecolor='black', alpha=0.7))

def update(frame):
    ax.clear()
    ax.plot_surface(X, Y, zb_red, cmap='terrain', alpha=0.6, linewidth=0, antialiased=True)
    
    h_actual = np.where(h_frames[frame] > UMBRAL_AGUA, h_frames[frame], 0)
    water_surface = zb_red + h_actual
    ax.plot_surface(X, Y, water_surface, cmap='Blues', alpha=0.8, linewidth=0, antialiased=True)
    
    ax.plot([presa_x_km, presa_x_km], [0, max(y)], [0, 150], color='red', linewidth=3, label='Presa')
    
    tiempo = frame * 30
    ax.text2D(0.02, 0.95, f'Tiempo: {tiempo:.0f} s ({tiempo/60:.1f} min)', 
              transform=ax.transAxes, color='white',
              bbox=dict(facecolor='black', alpha=0.7))
    
    ax.set_xlabel('Distancia (km)')
    ax.set_ylabel('Ancho (km)')
    ax.set_zlabel('Elevación (m)')
    ax.set_title('Simulación Rotura de Presa Urrá')
    ax.view_init(elev=25, azim=-60)
    ax.legend()
    
    return []

print("Generando video 3D... Esto puede tomar varios minutos...")
ani = FuncAnimation(fig, update, frames=len(h_frames), interval=1000/fps, blit=False)

output_path = '/home/nicolas/UrraDamBreachSimulation/real_topography_urra_3d.mp4'
ani.save(output_path, writer='ffmpeg', fps=fps, dpi=150)
print(f"Video 3D guardado: {output_path}")

plt.close()
print("¡Video 3D generado con éxito!")