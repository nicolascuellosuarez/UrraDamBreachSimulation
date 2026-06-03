import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob

files = sorted(glob.glob(
    '/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'
))

h = np.fromfile(files[-1], dtype=np.float64).reshape(520, 600)

x = np.arange(h.shape[1])
y = np.arange(h.shape[0])

X, Y = np.meshgrid(x, y)

fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(
    X,
    Y,
    h,
    cmap='Blues',
    linewidth=0,
    antialiased=True
)

ax.set_title('Estado Final - Profundidad del Agua')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Profundidad (m)')

fig.colorbar(surf, shrink=0.5, aspect=12)

plt.tight_layout()
plt.savefig('final_state_3d.png', dpi=300)
plt.show()

print(f"Min: {h.min():.2f}")
print(f"Max: {h.max():.2f}")
print(f"Mean: {h.mean():.2f}")