"""
    Lorentz Attractor


"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# Parameters
rho = 28
sigma = 10.0
beta = 8.0 / 3.0

# Initial states
init_states = [[1.0, 1.0, 1.0], [1.0, 0, 0]]
colors = ['#1C91A6', '#AD6C2E']


# Attractor function
def f(state, t):
    x, y, z = state  # Unpack the state vector
    return sigma * (y - x), x * (rho - z) - y, x * y - beta * z  # Derivatives


# Create time array
t = np.arange(0.0, 100.0, 0.001)

# Compute states
print('Solving equations...')
states = []
for ix, init_state in enumerate(init_states):
    states.append(odeint(f, init_state, t))

# Create figure
fig = plt.figure(figsize=(50, 50), dpi=1200)
ax = fig.gca(projection='3d')

# Plot main line
print('Plotting...')
prim_col = '#1C91A6'
for ix, state in enumerate(states):
    ax.plot(state[:, 0], state[:, 1], state[:, 2], c=colors[ix], lw=0.05, solid_capstyle='butt', zorder=0)

    # Plot neon lines
    neon_n = 2
    for n in range(1, neon_n + 1):
        ax.plot(state[:, 0], state[:, 1], state[:, 2], c=colors[ix], alpha=0.03, lw=0.05 + (1 * 1.05 * n), solid_capstyle='butt', zorder=1)

# Style plot
back_col = '#1E1E1E'
ax.set_facecolor(back_col)
ax.axis('off')

# Show
plt.show()

