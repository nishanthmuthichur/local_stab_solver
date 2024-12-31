import numpy as np

import num_sch as ns

import matplotlib.pyplot as plt

PI = np.pi

N_pts = 3

N = N_pts - 1

xo_vec = np.cos(np.linspace(0, N, N_pts) * PI / N)

x_vec = (xo_vec + 1) / 2

func_1 = np.cos(2 * PI * 10 * x_vec)

D_mat = ns.get_cheb_mat(xo_vec)

dfunc_1n = 2 * D_mat @ func_1
dfunc_1n = np.asarray(dfunc_1n)
dfunc_1n = dfunc_1n.flatten()

dfunc_1a = -2 * PI * 10 * np.sin(2 * PI * 10 * x_vec) 

fig, ax = plt.subplots(1, 1)
fig.set_dpi(300)

ax.plot(x_vec, func_1, 'b-')

ax.plot(x_vec, dfunc_1a, 'r-', label = 'Num')
ax.plot(x_vec, dfunc_1n, 'k--', label = 'Ana')
ax.grid()

