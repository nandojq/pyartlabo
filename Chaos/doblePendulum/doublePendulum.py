"""
    CHAOS COLLECTION: Double pendulum

"""

import sys
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap


# Functions
def deriv(y, t, l1, l2, m1, m2):
    """ Returns the first derivatives of y = theta1, theta1_dot, theta2, theta2_dot """

    # Unpack state variables
    theta1, theta1_dot, theta2, theta2_dot = y

    # Cosinus and sinus shortcut
    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    # First order differential equations for acceleration
    theta1_dotdot = (m2*g*np.sin(theta2)*c - m2*s*(l1*(theta1_dot**2)*c + l2*(theta2_dot**2)) - (m1+m2)*g*np.sin(theta1)) / (l1 * (m1 + m2*(s**2)))
    theta2_dotdot = ((m1+m2)*(l1*(theta1_dot**2)*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + m2*l2*(theta2_dot**2)*s*c) / (l2 * (m1 + m2*(s**2)))
    return theta1_dot, theta1_dotdot, theta2_dot, theta2_dotdot


def complement_color(color):
    """ Returns complementary color """
    r, g, b = color
    def hilo(a, b, c):
        if c < b: b, c = c, b
        if b < a: a, b = b, a
        if c < b: b, c = c, b
        return a + c

    k = hilo(r, g, b)
    return tuple(k - u for u in (r, g, b))

def CustomCmap(from_rgb,to_rgb):

    # from color r,g,b
    r1,g1,b1 = from_rgb

    # to color r,g,b
    r2,g2,b2 = to_rgb

    cdict = {'red': ((0, r1, r1),
                   (1, r2, r2)),
           'green': ((0, g1, g1),
                    (1, g2, g2)),
           'blue': ((0, b1, b1),
                   (1, b2, b2))}

    cmap = LinearSegmentedColormap('custom_cmap', cdict)
    return cmap


def generate_plot(x2, y2, speed, filename, type='raw'):
    """ Make final plot """

    # Create figure
    fig = plt.figure(figsize=(8.3333, 6.25), dpi=180)
    ax = fig.add_subplot(111)

    if type == 'raw':
        # x2 = x2[::5]
        # y2 = y2[::5]
        ax.plot(x2, y2)

    elif type == 'color':
        # Styling
        back_color = (np.random.random(), np.random.random(), np.random.random())
        shade_factor = 0.80
        shaded_back_color = (back_color[0] * (1 - shade_factor), back_color[1] * (1 - shade_factor), back_color[2] * (1 - shade_factor))
        corrected_shaded_back_color = [0, 0, 0]
        for ix, color_elem in enumerate(shaded_back_color):
            if color_elem < 0:
                corrected_shaded_back_color[ix] = 0
            elif color_elem > 1:
                corrected_shaded_back_color[ix] = 1
            else:
                corrected_shaded_back_color[ix] = color_elem
        fig.set_facecolor(tuple(corrected_shaded_back_color))

        lwidth = 1
        lcolor = complement_color(back_color)
        tint_factor = 0.5
        tinted_back_color = (lcolor[0] + (1 - lcolor[0]) * tint_factor, lcolor[1] + (1 - lcolor[0]) * tint_factor, lcolor[2] + (1 - lcolor[0]) * tint_factor)
        corrected_tinted_back_color = [0, 0, 0]
        for ix, color_elem in enumerate(tinted_back_color):
            if color_elem < 0:
                corrected_tinted_back_color[ix] = 0
            elif color_elem > 1:
                corrected_tinted_back_color[ix] = 1
            else:
                corrected_tinted_back_color[ix] = color_elem
        ax.plot(x2, y2, c=tuple(corrected_tinted_back_color), lw=lwidth, solid_capstyle='butt', zorder=0)
        neon_n = 10
        for n in range(1, neon_n+1):
            ax.plot(x2, y2, c=tuple(corrected_tinted_back_color), alpha=0.03, lw=2+(lwidth*1.05*n), solid_capstyle='butt', zorder=1)

    elif type == 'multicolor_cortado':
        # Styling
        back_color = (np.random.random(), np.random.random(), np.random.random())
        shade_factor = 0.80
        shaded_back_color = (back_color[0] * (1 - shade_factor), back_color[1] * (1 - shade_factor), back_color[2] * (1 - shade_factor))
        corrected_shaded_back_color = [0, 0, 0]
        for ix, color_elem in enumerate(shaded_back_color):
            if color_elem < 0:
                corrected_shaded_back_color[ix] = 0
            elif color_elem > 1:
                corrected_shaded_back_color[ix] = 1
            else:
                corrected_shaded_back_color[ix] = color_elem
        fig.set_facecolor(tuple(corrected_shaded_back_color))

        lwidth = 1
        lcolor = complement_color(back_color)
        tint_factor = 0.5
        tinted_back_color = (lcolor[0] + (1 - lcolor[0]) * tint_factor, lcolor[1] + (1 - lcolor[0]) * tint_factor, lcolor[2] + (1 - lcolor[0]) * tint_factor)
        corrected_tinted_back_color = [0, 0, 0]
        for ix, color_elem in enumerate(tinted_back_color):
            if color_elem < 0:
                corrected_tinted_back_color[ix] = 0
            elif color_elem > 1:
                corrected_tinted_back_color[ix] = 1
            else:
                corrected_tinted_back_color[ix] = color_elem

        cmap = CustomCmap(corrected_shaded_back_color, corrected_tinted_back_color)
        # norm = np.linalg.norm(speed)
        # normal_speed = speed / norm
        s = 2
        ax.scatter(x2, y2, c=cmap(speed), norm=speed, s=s, zorder=0)
        neon_n = 10
        for n in range(1, neon_n + 1):
            ax.scatter(x2, y2, c=cmap(speed), alpha=0.03, zorder=1, s=5 + (s * 1.05 * n))

    else:
        print('Nope')

    # Centre the image on the fixed anchor point, and ensure the axes are equal
    ax.set_aspect('equal', adjustable='box')
    # fig.set_facecolor('black')
    plt.axis('off')
    plt.savefig('{}.png'.format(filename), facecolor=fig.get_facecolor(), dpi=180)
    # plt.show()
    plt.cla()


# Run
if __name__ == '__main__':

    print('##### Double pendulum canvas initialized!')

    n = 10
    for paint in range(n):

        print('Painting {} of {}'.format(paint, n))

        # INITIAL CONDITIONS
        print('Randomizing initial conditions...')
        l1 = np.random.random()*10
        l2 = np.random.random()*10
        m1 = np.random.random()*10
        m2 = np.random.random()*10
        theta1_init = np.random.random()*180
        theta2_init = np.random.random()*180
        tmax = 5
        print(' L1 = {} m\n'.format(l1),
              'L2 = {} m\n'.format(l2),
              'm1 = {} kg\n'.format(m1),
              'm2 = {} kg\n'.format(m2),
              'Theta1 = {} deg\n'.format(theta1_init),
              'Theta2 = {} deg\n'.format(theta2_init),
              'Simulation time = {} s'.format(tmax))

        # SIMULATION PARAMETERS
        # The gravitational acceleration (m.s-2)
        g = 9.81
        # time step and time grid (all in s)
        dt = 0.01
        t = np.arange(0, tmax+dt, dt)
        # Initial conditions for speeds
        d_theta1_init = 0
        d_theta2_init = 0
        y0 = np.array([np.radians(theta1_init), d_theta1_init, np.radians(theta2_init), theta2_init])
        # Make an image every di time points, corresponding to a frame rate of fps
        fps = 10
        di = int(1 / fps / dt)

        # EQUATION SOLVING
        print('Solving differential equations...')
        # Do the numerical integration of the equations of motion
        y = odeint(deriv, y0, t, args=(l1, l2, m1, m2))

        # Unpack z and theta as a function of time
        theta1, theta1_dot, theta2, theta2_dot = y[:, 0], y[:, 1], y[:, 2], y[:, 3]

        # Convert to Cartesian coordinates of the two bob positions.
        x1 = l1 * np.sin(theta1)
        y1 = -l1 * np.cos(theta1)
        x2 = x1 + l2 * np.sin(theta2)
        y2 = y1 - l2 * np.cos(theta2)
        print('{} points'.format(len(x2)))

        print('Generating plot...')
        speed = (theta1_dot * l1) + (theta2_dot * l2)
        generate_plot(x2, y2, speed, 'DoublePendulum_{}'.format(paint), type='multicolor_cortado')
        print('Done!')
