"""
    Doble Pendulum WebApp (Random Parameters)
    --------------------------
    Made with Streamlit
    Creator: Fernando Jorge

"""

import streamlit as st
import os
import numpy as np
from DoblePendulum import solve_equations, generate_plot


class PyArtLaboApplication:

    def __init__(self):

        # Use the full page instead of a narrow central column
        st.set_page_config(page_title="Double Pendulum", layout="wide")

        # Define parameter widgets
        col1, col2 = st.columns([1, 3])
        self.l1 = col1.text_input('L1', value=round(np.random.random()*10, 3))
        self.l2 = col1.text_input('L2', value=round(np.random.random()*10, 3))
        self.m1 = col1.text_input('M1', value=round(np.random.random()*10, 3))
        self.m2 = col1.text_input('M2', value=round(np.random.random()*10, 3))
        self.theta1_init = col1.text_input('Theta 1', value=round(np.random.random()*180, 3))
        self.theta2_init = col1.text_input('Theta 2', value=round(np.random.random()*180, 3))
        self.g = 9.81
        self.sim_time = 200

        self.btn_save = col1.button('Save', on_click=self.save_figure)

        self.prim_col = col1.color_picker('Pick a primary color', value='#1C91A6')
        self.back_col = col1.color_picker('Pick a background color', value='#202020')

        # Compute trajectory
        x2, y2, speed = solve_equations(float(self.l1),
                                        float(self.l2),
                                        float(self.m1),
                                        float(self.m2),
                                        float(self.g),
                                        float(self.sim_time),
                                        float(self.theta1_init),
                                        float(self.theta2_init))

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
