import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
np.set_printoptions(threshold=np.nan)

N = 10
t = 1
E0 = 5
Delta = 1

D1 = E0 * np.diag(np.ones(N), 0)
D2 = -t * np.diag(np.ones(N - 1), 1)
D3 = -t * np.diag(np.ones(N - 1), -1)

H_HW = D1 + D2 + D3
H_PBC = np.copy(H_HW)
H_PBC[0, -1] = -t
H_PBC[-1, 0] = -t

# E_HW = np.linalg.eigvalsh(H_HW)
# E_PBC = np.linalg.eigvalsh(H_PBC)
# fig = plt.figure()
# ax = fig.gca()
# ax.plot(E_HW)
# ax.plot(E_PBC)
# plt.show()

blueprint = np.zeros((N, N))
blueprint[0, 0] = E0
blueprint[(1, -1, 0, 0), (0, 0, 1, -1)] = -t
H = np.zeros((N**2, N**2))
counter = 0
for i in range(N):
    for j in range(N):
        Hi = np.roll(blueprint, i, axis=0)
        Hj = np.roll(Hi, j, axis=1)
        H[:, counter] = Hj.flatten()
        counter += 1

E_H = np.linalg.eigvalsh(H)
E_H2 = np.linalg.eigvals(H)
x, y = np.meshgrid(np.arange(N), np.arange(N))
fig = plt.figure()
#ax = fig.gca(projection="3d")
#ax.plot_surface(x, y, E_H.reshape((N, N)))
ax = fig.gca()
ax.plot(E_H)
ax.plot(E_H2)

plt.show()
