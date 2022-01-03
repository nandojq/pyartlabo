"""
    Doble pendulum equation solver and plot creation

"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap


def solve_equations(l1, l2, m1, m2, g, sim_time, theta1_init, theta2_init):
    def deriv(y, t, l1, l2, m1, m2, g):
        """ Returns the first derivatives of y = theta1, theta1_dot, theta2, theta2_dot """

        # Unpack state variables
        theta1, theta1_dot, theta2, theta2_dot = y

        # Cosinus and sinus shortcut
        c, s = np.cos(theta1 - theta2), np.sin(theta1 - theta2)

        # First order differential equations for acceleration
        theta1_dotdot = (m2 * g * np.sin(theta2) * c - m2 * s * (l1 * (theta1_dot ** 2) * c + l2 * (theta2_dot ** 2)) - (m1 + m2) * g * np.sin(theta1)) / (l1 * (m1 + m2 * (s ** 2)))
        theta2_dotdot = ((m1 + m2) * (l1 * (theta1_dot ** 2) * s - g * np.sin(theta2) + g * np.sin(theta1) * c) + m2 * l2 * (theta2_dot ** 2) * s * c) / (l2 * (m1 + m2 * (s ** 2)))
        return theta1_dot, theta1_dotdot, theta2_dot, theta2_dotdot

    # time step and time grid (all in s)
    dt = 0.01
    t = np.arange(0, sim_time + dt, dt)

    # Initial conditions for speeds
    d_theta1_init = 0
    d_theta2_init = 0

    y0 = np.array([np.radians(theta1_init), d_theta1_init, np.radians(theta2_init), d_theta2_init])

    # Make an image every di time points, corresponding to a frame rate of fps
    fps = 10
    di = int(1 / fps / dt)

    # EQUATION SOLVING
    print('Solving differential equations...')
    # Do the numerical integration of the equations of motion
    y = odeint(deriv, y0, t, args=(l1, l2, m1, m2, g))

    # Unpack z and theta as a function of time
    theta1, theta1_dot, theta2, theta2_dot = y[:, 0], y[:, 1], y[:, 2], y[:, 3]

    # Convert to Cartesian coordinates of the two bob positions.
    x1 = l1 * np.sin(theta1)
    y1 = -l1 * np.cos(theta1)
    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)
    speed = (theta1_dot * l1) + (theta2_dot * l2)
    # st.write('{} points'.format(len(x2)))
    return x2, y2, speed


def generate_plot(x2, y2, back_col, prim_col):
    """ Make final plot """

    # Create figure
    fig = plt.figure(figsize=(20, 20), dpi=1200)
    ax = fig.add_subplot(111)

    # Styling
    fig.set_facecolor(back_col)

    lwidth = 1
    ax.plot(x2, y2, c=prim_col, lw=lwidth, solid_capstyle='butt', zorder=0)
    neon_n = 8
    for n in range(1, neon_n + 1):
        ax.plot(x2, y2, c=prim_col, alpha=0.03, lw=2 + (lwidth * 1.05 * n), solid_capstyle='butt', zorder=1)

    fig.set_facecolor(back_col)
    plt.axis('off')
    # plt.tight_layout()
    return fig


