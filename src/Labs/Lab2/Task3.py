import numpy as np
import matplotlib.pyplot as plt
from time import time

xmin, xmax, step = -2, 2, 0.004
N = int((xmax - xmin) / step) + 1
max_iter = 30

roots = np.array([
    1 + 0j,
    np.exp(2j * np.pi / 3),
    np.exp(-2j * np.pi / 3)
], dtype=np.complex64)

x = np.linspace(xmin, xmax, N, dtype=np.float32)
y = np.linspace(xmin, xmax, N, dtype=np.float32)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y
Z = Z.astype(np.complex64)

start = time()
Zcur = Z.copy()
for i in range(max_iter):
    F = Zcur**3 - 1
    dF = 3 * Zcur**2
    # безопасное деление
    mask = np.abs(dF) > 1e-12
    Zcur[mask] = Zcur[mask] - F[mask] / dF[mask]
    # "плохие" точки игнорируются в дальнейшем шаге

# определяем к ближайшему корню принадлежит каждая точка
dist = np.abs(Zcur[..., None] - roots)
root_ind = np.argmin(dist, axis=2)

end = time()
plt.figure(figsize=(8, 8))
plt.imshow(root_ind, extent=(xmin, xmax, xmin, xmax), origin='lower', cmap='tab10')
plt.xlabel('Re(z)')
plt.ylabel('Im(z)')
plt.title(f'Newton fractal ($z^3-1=0$), solved in {end-start:.2f} s')
plt.tight_layout()
plt.show()
