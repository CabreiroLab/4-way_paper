#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 09:16:20 2017

@author: Povilas
"""

"""
===============
Rain simulation
===============

Simulates rain drops on a surface by animating the scale and opacity
of 50 scatter points.

Author: Nicolas P. Rougier
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



replicate='Rep4'
plate='PM1'
indx=13
#for indx in range(1,13):
#rfile="{}/{}/NGM_Control_Rep1.tif".format(sourceloc,replicate)
tfile="{}/{}/{}/{}_{}.tif".format(sourceloc,replicate,plate,plate,str(indx).zfill(3))

if indx in alldeletions[replicate][plate].keys():
    delworms=alldeletions[replicate][plate][indx]['Worms']
else:
    delworms=[]
    

image = tiff.imread(tfile)
image_clean=wormdel(image,delworms)
imghsv=color.rgb2hsv(image)
imgrgb=color.hsv2rgb(imghsv)

h=imghsv[:,:,0]
s=imghsv[:,:,1]
v=imghsv[:,:,2]

plt.imshow(v)

hmin=np.percentile(h,95)
print hmin


labeled_worms,hmax,hthres=labeling2(imghsv,3.4757,-0.6754,0.02,hmin=0.28)
plt.imshow(labeled_worms)
print hmax,hthres





def plot_comparison(figs, labels):
    nfig=len(figs)
    nlab=len(labels)

    figure, axes = plt.subplots(ncols=nfig, figsize=(4*nfig, 4), sharex=True,
                                   sharey=True)
    for fid,fplot in enumerate(figs):
        ax = axes[fid]
        ax.imshow(fplot) #, cmap=plt.cm.gray
        ax.set_title(labels[fid])
        ax.axis('off')
        ax.set_adjustable('box-forced')

    return figure,axes
    


nfig=2
figure, axes = plt.subplots(ncols=nfig, figsize=(8*nfig, 8), sharex=True,
                                   sharey=True)

ax = axes[0]
ax.imshow(labeled_worms)
ax = axes[1]
ax.imshow(imgrgb)

ax = axes[1]
ax.imshow(imghsv)

figure.show()



# Create new Figure and an Axes which fills it.
fig = plt.figure(figsize=(7, 7))
ax = fig.add_axes([0, 0, 1, 1], frameon=False)
ax.set_xlim(0, 1), ax.set_xticks([])
ax.set_ylim(0, 1), ax.set_yticks([])

# Create rain data
n_drops = 50
rain_drops = np.zeros(n_drops, dtype=[('position', float, 2),
                                      ('size',     float, 1),
                                      ('growth',   float, 1),
                                      ('color',    float, 4)])

# Initialize the raindrops in random positions and with
# random growth rates.
rain_drops['position'] = np.random.uniform(0, 1, (n_drops, 2))
rain_drops['growth'] = np.random.uniform(50, 200, n_drops)

# Construct the scatter which we will update during animation
# as the raindrops develop.
scat = ax.scatter(rain_drops['position'][:, 0], rain_drops['position'][:, 1],
                  s=rain_drops['size'], lw=0.5, edgecolors=rain_drops['color'],
                  facecolors='none')


def update(frame_number):
    # Get an index which we can use to re-spawn the oldest raindrop.
    current_index = frame_number % n_drops

    # Make all colors more transparent as time progresses.
    rain_drops['color'][:, 3] -= 1.0/len(rain_drops)
    rain_drops['color'][:, 3] = np.clip(rain_drops['color'][:, 3], 0, 1)

    # Make all circles bigger.
    rain_drops['size'] += rain_drops['growth']

    # Pick a new position for oldest rain drop, resetting its size,
    # color and growth factor.
    rain_drops['position'][current_index] = np.random.uniform(0, 1, 2)
    rain_drops['size'][current_index] = 5
    rain_drops['color'][current_index] = (0, 0, 0, 1)
    rain_drops['growth'][current_index] = np.random.uniform(50, 200)

    # Update the scatter collection, with the new colors, sizes and positions.
    scat.set_edgecolors(rain_drops['color'])
    scat.set_sizes(rain_drops['size'])
    scat.set_offsets(rain_drops['position'])


# Construct the animation, using the update function as the animation
# director.
animation = FuncAnimation(fig, update, interval=10)
plt.show()




def labeling(v, thr_1, thr_2, elevation='sobel'):
	markers = np.zeros_like(v)
	markers[v < thr_1] = 1
	markers[v > thr_2] = 2

	if elevation == 'canny':
		elevation_map = canny(v)

	else:
		elevation_map = sobel(v)

	segmentation = watershed(elevation_map, markers)
	segmentation = ndi.binary_fill_holes(segmentation - 1)
	labeled_worms, _ = ndi.label(segmentation)

	return labeled_worms



def thresholding2(imghsv, thresholds, sizes, delworms):
	# file = '20170503/image_011.tif'
	# Changed from skio.imread, which could not handle 16bit TIFF

	if len(sizes) == 2:
		vc_size, hc_size = sizes
	else:
		print 'Fix your sizes parameters!'
		vc_size = sizes[0]
		hc_size = sizes[1]

	vmin, vmax, hmin, hmax = thresholds
	# image.shape
	# image.dtype


	h = imghsv[:, :, 0]
	# s =imghsv[:,:,1]
	v = imghsv[:, :, 2]

	vmask = (v > vmin).astype(np.float64);  # Thresholding in Brightness
	vmask = (vmask < vmax).astype(np.float64)

	hmask = (h > hmin).astype(np.float64);  # Thresholding in Hue
	hmask = (hmask < hmax).astype(np.float64)

	# Works reallt well
	# Morphological filtering
	v_closing = closing(vmask, selem=disk(vc_size))
	h_closing = closing(hmask, selem=disk(hc_size))
	# plt.imshow(hopened); plt.title('Hue mask')
	# plt.imshow(hdilated); plt.title('Hue mask')

	img2 = imghsv.copy()
	img2[h_closing.astype(bool), :] = 0;  # Set the pixels to zero by Hue mask
	img2[v_closing.astype(bool), :] = 0;  # Set the pixels to zero by Lightness mask
	for worm in delworms:
		x1, y1, x2, y2 = worm
		xs = [x1, x2]
		ys = [y1, y2]

		img2[min(ys):max(ys), min(xs):max(xs), :] = 0

	img2rgb = color.hsv2rgb(img2)

	return v_closing, h_closing, img2, img2rgb