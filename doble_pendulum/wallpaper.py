

import os
from DoblePendulum import solve_equations, generate_plot
import matplotlib.pyplot as plt


# Figures
figures = {
    'Saco': [1.553, 4.281, 8.043, 1.167, 9.81, 200, 89.725, 17.49, '#1C91A6'],
    'Aro': [4.777, 0.966, 5.011, 1.532, 9.81, 140, 162.488, 180.0, '#AD6C2E'],
    'Red': [1.553, 6.881, 8.043, 2.8, 9.81, 194, 73.872, 17.79, '#25B777'],
}
back_color = '#1E1E1E'

# Create figure
fig, axs = plt.subplots(1, 3, figsize=(54, 18), facecolor=back_color)

for ix, figure in enumerate(figures):
    # Extract figure data
    fig_data = figures[figure]

    # Compute trajectory
    x2, y2, speed = solve_equations(fig_data[0],
                                    fig_data[1],
                                    fig_data[2],
                                    fig_data[3],
                                    fig_data[4],
                                    fig_data[5],
                                    fig_data[6],
                                    fig_data[7])

    # Plot figure
    lwidth = 1
    axs[ix].plot(x2, y2, c=fig_data[8], lw=lwidth, solid_capstyle='butt', zorder=0)
    neon_n = 8
    for n in range(1, neon_n + 1):
        axs[ix].plot(x2, y2, c=fig_data[8], alpha=0.03, lw=2 + (lwidth * 1.05 * n), solid_capstyle='butt', zorder=1)
    axs[ix].axis('off')

# Create figure
plt.tight_layout()

# Save fig
print('Saving figure...')
fig.savefig(os.path.join(r'C:\Users\fjorg\Desktop\Git\pyartlabo\doble_pendulum\Favs', 'wallpaper.png'), format='png', dpi=100)
print('Done!')
