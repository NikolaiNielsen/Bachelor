import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
np.set_printoptions(threshold=np.nan)

N = 1001
L = 10
xarr = np.linspace(-L, L, N)
dx = xarr[1] - xarr[0]
V = np.zeros(N)
V0 = 0.5
V[[np.floor(N / 6).astype(int), np.floor(N / 2).astype(int), np.floor(5 * N / 6).astype(int)]] = 1 / dx

D1 = np.diag(np.ones(N), 0)
D2 = -2 * np.diag(np.ones(N - 1), 1)
D3 = np.diag(np.ones(N - 1), -1)
extra = np.zeros((N, N))
extra[[0, -1], [-1, 0]] = 1

H = (D1 + D2 + D3 + extra) / (2 * dx**2) + V0 * np.diag(V, 0)

E, psi = np.linalg.eigh(H)
prob = psi * np.conj(psi)

fig = plt.figure()
ax = fig.gca()
ax.plot(xarr, prob[:, 0])

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.plot(xarr, V)

plt.show()
