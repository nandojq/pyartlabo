"""
    CHAOS COLLECTION: Mandelbrot set

"""

import numpy as np
from tqdm import tqdm
from keras.preprocessing.image import array_to_img
import matplotlib.pyplot as plt


def create_img_array(px: int, py: int):
    # Check pixel ratio
    ratio = px/py
    if ratio != 4/3:
        py = int(px * 3 / 4)
        print('Pixel ratio inserted is not 4/3. Image will be resized to ({}, {})'.format(px, py))
    # Create image array
    array = np.zeros((py, px, 1))
    return array


def mandelbrot(array, portion_dt, max_it, thres):
    # Get width and height
    w = array.shape[1]
    h = array.shape[0]
    # Get limits from portion data [Re_center, Im_center, zoom]
    re_min = portion_dt[0] - (portion_dt[2] * 0.004)
    re_max = portion_dt[0] + (portion_dt[2] * 0.004)
    im_min = portion_dt[1] - (portion_dt[2] * 0.003)
    im_max = portion_dt[1] + (portion_dt[2] * 0.003)
    print(re_min, re_max, im_min, im_max)
    # Iterate over pixels
    for x in tqdm(range(w)):
        for y in range(h):
            # Map a and b relative to x and y
            a = ((re_max - re_min) / w) * x + re_min
            b = ((im_max - im_min) / h) * y + im_min
            # Save initial point parameters
            ca = a
            cb = b
            # Start Mandelbrot iteration
            it = 0
            conj = 0
            color = 0
            while it < max_it:
                # Compute this iteration parameters
                aa = a * a - b * b + ca
                bb = 2 * a * b + cb
                # Save parameters for next iteration
                a = aa
                b = bb
                # Break if point went to infinity
                if conj > thres:
                    color = ((255 / max_it) * it)
                    break
                conj = (a * a) + (b * b)
                it += 1
            img_array[y, x, :] = color
    # Return image
    return img_array


def create_figure(array, dpi=None, cmap='magma_r'):
    # Convert img data array to image
    print('Converting data array to image...')
    img_pil = array_to_img(array)

    # Choose DPI
    if dpi is None:
        dpi = matplotlib.rcParams['savefig.dpi']
    dpi = float(dpi)

    # Create figure and add axes
    h = 4000 / dpi
    w = 3000 / dpi
    fig = plt.figure(figsize=(h, w), dpi=dpi, frameon=False)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # Plot image from data array
    ax.imshow(img_pil, cmap=cmap, vmin=0, vmax=255)

    # Remove axes from plot
    ax.axis('off')

    # Return image
    return fig


if __name__ == '__main__':

    type_run = 'Create Data'
    data_to_load = r'data/Julia Island_1920_1440_500it.npy'

    if type_run == 'Create Data':

        # Choose Mandelbrot set portion to render. [Re_center, Im_center, zoom]
        portions = {'Mandelbrot': [-0.737442923, 0.02283105, 438],
                    'Julia Island': [-1.768778822, -0.001738915, 2192091967]
                    }
        chosen_portion = 'Julia Island'
        portion = portions[chosen_portion]

        # Create image array
        print('Initializing image array...')
        x_pixels = 1920
        y_pixels = 1440
        img_array = create_img_array(x_pixels, y_pixels)
        print('Done! Image size (height, width, channels): {}'.format(img_array.shape))

        # Generate Mandelbrot set in image
        max_iterations = 500
        threshold = 4
        print('Generating Mandelbrot set in image size {}'.format(img_array.shape))
        mdt = mandelbrot(img_array, portion, max_iterations, threshold)
        print('Done!')

        # Save data array
        data_filename = r'data/{}_{}_{}_{}it.npy'.format(chosen_portion, x_pixels, y_pixels, max_iterations)
        np.save(data_filename, mdt)
        print('Data saved in file {}'.format(data_filename))

    else:

        # Load data
        mdt = np.load(data_to_load)

        # Extract parameters
        sptd = data_to_load.split('/')[1].split('it')[0].split('_')
        chosen_portion = sptd[0]
        x_pixels = sptd[1]
        y_pixels = sptd[2]
        max_iterations = sptd[3]

    # Create figure
    print('Creating image...')
    dpi = 120
    colormap = 'inferno'
    fig = create_figure(mdt, dpi, colormap)

    # Save figure
    fig.savefig(r'imgs/{}_{}_{}_{}it.png'.format(chosen_portion, x_pixels, y_pixels, max_iterations), format='png', dpi=100)
    print('Done!')
