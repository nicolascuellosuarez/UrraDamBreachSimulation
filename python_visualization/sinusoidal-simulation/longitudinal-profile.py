import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob('/home/nicolas/UrraDamBreachSimulation/snapshots/h_*.bin'))

h0 = np.fromfile(files[0], dtype=np.float64).reshape(520, 600)
hf = np.fromfile(files[-1], dtype=np.float64).reshape(520, 600)

perfil_ini = h0.mean(axis=0)
perfil_fin = hf.mean(axis=0)

plt.figure(figsize=(12,5))

plt.plot(perfil_ini, label='Inicial')
plt.plot(perfil_fin, label='Final')

plt.axvline(314, color='red', linestyle='--', label='Presa')

plt.xlabel('Distancia (celdas)')
plt.ylabel('Profundidad media (m)')
plt.title('Perfil longitudinal promedio')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('longitudinal_profile.png', dpi=300)
plt.show()