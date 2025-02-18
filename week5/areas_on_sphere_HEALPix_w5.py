########################################################
#						                               #
#    In Class Work: Areas on the Sphere and HEALPix    #
#		               2.18.2025			           #
#	             Author: Eliza Frankel                 #
#                                                      #                        
########################################################

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
import healpy as hp

pi = np.pi

###  Task 1  ###

# EAF - randomly generating 1,000,000 points on the surface of a sphere.
ra = 360. * (random(1000000))
dec=(180 / pi) * np.arcsin(1. - random(1000000) * 2.)

###  Task 2  ###

# EAF - nside = # of divisions along the side of a base-res. pixel
#	needed to reach a desired high-res partition (JPL)
nside = 1
# EAF - lonlat=True means that inputted ra and dec are the RA and Dec
#       This is the conversion from sky coordinates to pixels
inside_region = hp.ang2pix(nside, ra, dec, lonlat=True)

# EAF - the following gives the area of an nside. Returns square radians (or degrees)
area = hp.nside2pixarea(nside, degrees=False)
print("The area of an nside=1 HEALPix is {} radians^2.\n".format(round(area, 5)))

###  TASK 3  ###

# EAF - provides the number of points in each pixel.
#       Returns the indices of the pixels [0] and the number of counts per pixel [1]
uniqueness = np.unique(inside_region, return_counts=True)
print("The number of points in each pixel have been found to be:\n{}".format(uniqueness[1]))
print("Though these are not all the same, there are {} total points in the {} pixels.".format(len(inside_region), len(uniqueness[1])))
print("The average number of points per pixel is {}.".format(round(np.average(uniqueness[1]), 2)))
print("The values of points are close enough that it can be considered consistent with equal area pixels.\n")

###  TASK 4  ###

# EAF - creates Boolean mask for the points in specific pixels
i2 = inside_region == 2
i5 = inside_region == 5
i8 = inside_region == 8

# EAF - plotting the pixels on top of all points. For some reason, all RA and Dec were always 
#       plotting on top of the pixels, regardless of the order of my statements.
#       The zorder command allows me to choose the order things are plotted.
plt.figure(figsize=(16, 9))
plt.scatter(ra, dec, marker='o', s=50, color='lightskyblue', label='All Points', zorder=1)
plt.scatter(ra[i2], dec[i2], marker='o', s=50, color='goldenrod', label='Pixel 2', zorder=4)
plt.scatter(ra[i5], dec[i5], marker='o', s=50, color='deeppink', label='Pixel 5', zorder=3)
plt.scatter(ra[i8], dec[i8], marker='o', s=50, color='forestgreen', label='Pixel 8' ,zorder=3)

plt.legend(loc='upper right', fontsize=15)
plt.title('nside = 1', size=20)
plt.xticks(size=15)
plt.xlabel('Right Ascension [degrees]', size=15)
plt.yticks(size=15)
plt.ylabel('Declination [degrees]', size=15)
plt.tight_layout()
plt.savefig('nside_1.png', dpi=300)
plt.show()
plt.close()


###  TASK 5  ###
nside = 2
inside_nside2 = hp.ang2pix(nside, ra, dec, lonlat=True)

# EAF - determining things in Pixel 5
pix = inside_nside2[i5]
uniqueness2 = np.unique(pix, return_counts=True)
print('The pixels from nside = 2 that are in Pixel 5 for nside = 1 are {}.'.format(uniqueness2[0]))

i14 = inside_nside2 == 14
i21 = inside_nside2 == 21
i22 = inside_nside2 == 22
i30 = inside_nside2 == 30

plt.figure(figsize=(16, 9))
plt.scatter(ra, dec, marker='o', s=50, color='lightskyblue', label='All Points', zorder=1)
plt.scatter(ra[i14], dec[i14], marker='o', s=50, color='coral', label='Pixel 14', zorder=3)
plt.scatter(ra[i21], dec[i21], marker='o', s=50, color='palegreen', label='Pixel 21', zorder=4)
plt.scatter(ra[i22], dec[i22], marker='o', s=50, color='darkorchid', label='Pixel 22', zorder=5)
plt.scatter(ra[i30], dec[i30], marker='o', s=50, color='cornflowerblue', label='Pixel 30', zorder=6)

plt.legend(loc='upper right', fontsize=15)
plt.title('nside = 2', size=20)
plt.xticks(size=15)
plt.xlabel('Right Ascension [degrees]', size=15)
plt.yticks(size=15)
plt.ylabel('Declination [degrees]', size=15)
plt.tight_layout()
plt.savefig('nside_2.png', dpi=300)
plt.show()
plt.close()

# EAF - I am unsure if this is the correct way to determine the same pixels.
#       I didn't want to plot all 48 pixels because I didn't know an efficient way to do so.
