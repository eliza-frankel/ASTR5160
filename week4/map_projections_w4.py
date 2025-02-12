#!/usr/local/Anaconda2024/bin/python3.11

#############################################
#											#
#   Week 4 in class work: Map Projections	#
#       Author: Eliza Frankel				#
#											#
#############################################

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
from astropy.coordinates import SkyCoord
from astropy import wcs
import astropy.units as u
import dustmaps
from dustmaps.config import config
from dustmaps.sfd import SFDQuery


#### Task 1

# EAF - randomly creating 10000 points on the surface of a sphere equally in area
ra = 2 * np.pi * (random(10000) - 0.5)
dec = np.arcsin(1. - random(10000) * 2)
plt.plot(ra, dec, marker='.', color='mediumpurple', linestyle='none')

plt.xlabel('RA [rad]')
plt.ylabel('Dec [rad]')
plt.title('Random RA and Dec on a Cartesian Grid')

plt.savefig('random_cartesian.png')
plt.show()

### Task 2

fig = plt.figure(figsize=(10,7))
# EAF - the following line changes the projection style and makes subplots
ax = fig.add_subplot(111, projection='aitoff')
# EAF - s = size of the marker, alpha = opacity
ax.scatter(ra, dec, marker='o', color='dodgerblue', s=1, alpha=0.6)

# EAF - changing the axes from degrees to hours
xlabels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
ax.set_xticklabels(xlabels, weight=800)

# EAF - creates a blue grid over the projection
ax.grid(color='b', linestyle='dashed', linewidth=1)
ax.set_title('Aitoff Projection\n')
plt.savefig('aitoff_projection.png')

plt.show()
plt.close()

# Lambert projection

# EAF - changing the figure size so that the axis labels are not overlapping
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='lambert')
ax.scatter(ra, dec, marker='o', color='dodgerblue', s=1, alpha=0.6)

xlabels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
ax.set_xticklabels(xlabels, weight=800)

ax.grid(color='b', linestyle='dashed', linewidth=1)
ax.set_title('Lambert Projection\n')
plt.savefig('lambert_projection.png')

plt.show()



#### TASK 4

dustdir = '/d/scratch/ASTR5160/data/dust/v0_1/maps'
config["data_dir"] = dustdir
sfd = SFDQuery()

# EAF - creating RA and Dec values, binned by 1 degree
RA = np.linspace(0.5, 359.5, 360)
Dec = np.linspace(-89.5, 89.5, 180)

# EAF - creates meshgrid using RA and Dec values
RA_mesh, Dec_mesh = np.meshgrid(RA, Dec)
coordinates = SkyCoord(ra=RA_mesh*u.deg, dec= Dec_mesh*u.deg) 
ebmv = sfd(coordinates)

w = wcs.WCS(naxis=2)
w.wcs.ctype = ["RA---AIT", "DEC--AIT"]
x, y = w.wcs_world2pix(RA_mesh, Dec_mesh, 1)


fig = plt.figure(figsize=(16,9))

# EAF - left subplot
ax = fig.add_subplot(121, projection='aitoff')
ax.scatter(RA_mesh, Dec_mesh, marker='.', color='lightskyblue', s=0.7, alpha=0.5)
xlabels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
ax.set_xticklabels(xlabels, weight=800)
ax.set_title('Degrees\n')
ax.grid(color='b', linestyle='dashed', linewidth=0.75)

# EAF - right subplot
ax = fig.add_subplot(122, projection='aitoff')
ax.scatter(x, y, marker='.', color='lightskyblue', s=0.7, alpha=0.5)
xlabels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
ax.set_xticklabels(xlabels, weight=800)
ax.set_title('Cartesian\n')
ax.grid(color='b', linestyle='dashed', linewidth=0.75)

plt.savefig('degrees_vs_cartesian.png')
plt.show()

# EAF - Yes, the difference is very notable. In the degrees subplot, there are 
#		notable horizontal stripes, while in the cartesian plot it is generally uniform,
#		with more concentration around the poles.


# Contour Plot

fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111, projection='aitoff')
ax.set_title('Dust Contour Plot')

xlabels = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
ax.set_xticklabels(xlabels, weight=800)

ax.grid(color='b', linestyle='dashed', linewidth=0.75)

plt.contourf(x, y, ebmv, levels=100, vmax=55, vmin=11, cmap='plasma')
plt.colorbar()
plt.savefig('contour_plot.png')
plt.show()

# EAF - This contour plot code takes a bit of time to run. I also was not able to get a
#		contour using the levels. After 500, it runs even slower and has no improvement
#		visually. I plan on working on this code more to try and figure it out,
#		but I wanted to submit what I have so far.