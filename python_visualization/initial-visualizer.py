import numpy as np
import matplotlib.pyplot as plt

NX, NY = 600, 520
DX, DY = 50.0, 50.0

try:
    zb = np.fromfile('/home/nicolas/UrraDamBreachSimulation/data/urra_specific_topography.bin', dtype=np.float64).reshape(NY, NX)
    print("Topografía real cargada")
except:
    x = np.arange(NX) * DX
    y = np.arange(NY) * DY
    X, Y = np.meshgrid(x, y)
    zb = 50.0 - 0.00015 * X + 2.0 * np.sin(0.001 * Y)
    print("Topografía sintética generada")

h = np.zeros((NY, NX))
H_UPSTREAM, H_DOWNSTREAM = 130.5, 68.0
DAM_X_POS = 314  

for j in range(NY):
    for i in range(NX):
        if i < DAM_X_POS and zb[j, i] < H_UPSTREAM:
            h[j, i] = H_UPSTREAM - zb[j, i]
        else:
            h[j, i] = 0.01
        
        if i >= DAM_X_POS and zb[j, i] < H_DOWNSTREAM + 10.0:
            rio_depth = H_DOWNSTREAM - zb[j, i]
            if rio_depth > h[j, i]:
                h[j, i] = rio_depth
        
        if h[j, i] < 0.01:
            h[j, i] = 0.01

presa_x = DAM_X_POS * DX
presa_y = 0

fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111, projection='3d')

X = np.arange(NX) * DX
Y = np.arange(NY) * DY
X, Y = np.meshgrid(X, Y)

surf = ax.plot_surface(X, Y, zb, cmap='terrain', alpha=0.7)
water_surface = zb + h
water = ax.plot_surface(X, Y, water_surface, cmap='Blues', alpha=0.6)

ax.axvline(x=presa_x, color='red', linestyle='--', linewidth=2, label='Presa')

ax.set_xlabel('Distancia (m)')
ax.set_ylabel('Ancho (m)')
ax.set_zlabel('Elevación (m)')
ax.set_title('Estado inicial - Presa Urrá')
ax.legend()

plt.savefig('initial_state_3d.png')
plt.show()
print(f"Visualización guardada. Presa en X={presa_x} m (celda {DAM_X_POS})")