"""
    Doble Pendulum WebApp (Interactive)
    --------------------------
    Made with Streamlit
    Creator: Fernando Jorge

"""

import streamlit as st
import os
from DoblePendulum import solve_equations, generate_plot


class PyArtLaboApplication:

    def __init__(self):

        # Use the full page instead of a narrow central column
        st.set_page_config(page_title="Double Pendulum", layout="wide")

        # Define parameter widgets
        col1, col2 = st.columns([1, 3])
        self.l1 = col1.number_input(
            'L1', step=0.1, format="%.3f", min_value=0.0, max_value=10.0, value=2.234
        )
        self.l2 = col1.number_input(
            'L2', step=0.1, format="%.3f", min_value=0.0, max_value=10.0, value=4.325
        )
        self.m1 = col1.number_input(
            'M1', step=0.1, format="%.3f", min_value=0.0, max_value=10.0, value=5.425
        )
        self.m2 = col1.number_input(
            'M2', step=0.1, format="%.3f", min_value=0.0, max_value=10.0, value=2.846
        )
        self.theta1_init = col1.number_input(
            'Theta 1 init', step=0.1, format="%.3f", min_value=0.0, max_value=180.0, value=0.0
        )
        self.theta2_init = col1.number_input(
            'Theta 2 init', step=0.1, format="%.3f", min_value=0.0, max_value=180.0, value=180.0
        )
        self.g = col1.number_input(
            'Gravity', step=0.1, format="%.3f", min_value=0.0, max_value=100.0, value=9.81
        )
        self.sim_time = col1.slider(
            'Simulation time', step=1, format="%.3f", min_value=0, max_value=500, value=5
        )
        self.btn_save = col1.button('Save', on_click=self.save_figure)

        self.prim_col = col1.color_picker('Pick a primary color', value='#1C91A6')
        self.back_col = col1.color_picker('Pick a background color', value='#202020')

        # Compute trajectory
        x2, y2, speed = solve_equations(self.l1,
                                        self.l2,
                                        self.m1,
                                        self.m2,
                                        self.g,
                                        self.sim_time,
                                        self.theta1_init,
                                        self.theta2_init)

        # Create figure
        self.fig = generate_plot(x2, y2, self.back_col, self.prim_col)
        col2.pyplot(self.fig, use_column_width=True)

    def save_figure(self):
        # Build file name with parameters
        filename = '{}_{}_{}_{}_{}_{}_{}_{}.svg'.format(self.l1,
                                                        self.l2,
                                                        self.m1,
                                                        self.m2,
                                                        self.g,
                                                        self.sim_time,
                                                        self.theta1_init,
                                                        self.theta2_init)

        # Save figure
        self.fig.savefig(os.path.join(r'C:\Users\fjorg\Desktop\Git\pyartlabo\doble_pendulum\Favs', filename), format='svg', dpi=1200)
        st.write('File Saved!')


if __name__ == '__main__':
    app = PyArtLaboApplication()
