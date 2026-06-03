import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os

dam_x = 314
fps = 10

snapshot_dir = '/home/nicolas/UrraDamBreachSimulation/snapshots'
files = sorted(glob.glob(os.path.join(snapshot_dir, 'h_*.bin')))
print(f"Encontrados {len(files)} snapshots")

if not files:
    print("ERROR: No hay snapshots")
    exit(1)

fig, ax = plt.subplots(figsize=(12, 6))

h0 = np.fromfile(files[0], dtype=np.float64).reshape(520, 600)
im = ax.imshow(h0, cmap='Blues', vmin=0, vmax=None)
plt.colorbar(im, ax=ax, label='Profundidad (m)')
ax.axvline(x=dam_x, color='red', linestyle='--', linewidth=2, label='Presa')
ax.set_xlabel('Distancia (celdas)')
ax.set_ylabel('Ancho (celdas)')
ax.set_title('Simulación Rotura de Presa Urrá')
ax.legend()

time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='white',
                    bbox=dict(facecolor='black', alpha=0.7))

def update(frame):
    h = np.fromfile(files[frame], dtype=np.float64).reshape(520, 600)
    im.set_array(h)
    
    tiempo = frame * 30
    time_text.set_text(f'Tiempo: {tiempo:.0f} s ({tiempo/60:.1f} min)')
    
    im.set_clim(vmin=0, vmax=np.percentile(h, 95))
    
    return [im, time_text]

ani = FuncAnimation(fig, update, frames=len(files), interval=1000/fps, blit=True)

output_path = '/home/nicolas/UrraDamBreachSimulation/real_topography_urra_2d.mp4'
ani.save(output_path, writer='ffmpeg', fps=fps, dpi=150)
print(f"Video guardado: {output_path}")

gif_path = '/home/nicolas/UrraDamBreachSimulation/real_topography_urra_2d.gif'
ani.save(gif_path, writer='pillow', fps=fps)
print(f"GIF guardado: {gif_path}")

plt.close()
print("¡Video 2D generado con éxito!")