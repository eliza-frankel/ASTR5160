#!/usr/local/Anaconda2024/bin/python3.11

#########################################
#										#
#   Week 3 in class work: Dust Maps		#
#       Author: Eliza Frankel			#
#										#
#########################################

import dustmaps
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from dustmaps.config import config

### Following the Lecture Notes

# Converting to Galactic Coordinates
RA = '12h42m30s'
Dec = '+41d12m00s'
c = SkyCoord(RA, Dec).galactic

dustdir = '/d/scratch/ASTR5160/data/dust/v0_1/maps'
config["data_dir"] = dustdir
from dustmaps.sfd import SFDQuery
sfd = SFDQuery()
print('With Galactic coordinates:', sfd(c))

# Not converting to Galactic coords.
RA = '12h42m30s'
Dec = '+41d12m00s'
c = SkyCoord(RA, Dec)
ebv = sfd(c)
print('Without Galactic coordinates:', ebv)

# Extinction for SDSS ugriz filters
ugriz = np.array([4.239, 3.303, 2.285, 1.698, 1.263])
A = ebv*ugriz
print('Extinction for the SDSS ugriz filters:', A)



### TASK 1

ra1 = 246.933  # degrees
dec1 = 40.795  # degrees
z1 = 1.08  # redshift
# EAF - using SDSS Navigator tool,
g1 = 18.81
r1 = 18.73
i1 = 18.82

coord1 = SkyCoord(ra1*u.deg, dec1*u.deg)
ebv1 = sfd(coord1)  # EAF - accounting for extinction
A_1 = ebv1 * ugriz  # EAF - gets the extinction for each filter. The extinction is for (u, g, r, i, z)

ra2 = 236.562  # degrees
dec2 = 2.440  #degrees
z2 = 1.08  # redshift
# EAF - again using SDSS Navigator
g2 = 19.10
r2 = 18.79
i2 = 18.73

coord2 = SkyCoord(ra2*u.deg, dec2*u.deg)
ebv2 = sfd(coord2)
A_2 = ebv2 * ugriz 

# EAF - for the g-r vs r-i plots
g_r1 = g1 - r1
r_i1 = r1 - i1
g_r1_ext = (g1 - A_1[1]) - (r1 - A_1[2]) # EAF - the following color_ext variables are accounting for extinction
r_i1_ext = (r1 - A_1[2]) - (i1 - A_1[3])

g_r2 = g2 - r2
r_i2 = r2 - i2
g_r2_ext = (g2 - A_2[1]) - (r2 - A_2[2])
r_i2_ext = (r2 - A_2[2]) - (i2 - A_2[3])

# Plotting both quasars
plt.plot(g_r1, r_i1, marker='o', markersize=10, linestyle='None', label='Quasar 1')
plt.plot(g_r2, r_i2, marker='o', markersize=10, linestyle='None', label='Quasar 2')
plt.xlim(-0.03, 0.35)  #EAF - resizing the figure
plt.ylim(-0.1, 0.07)

plt.xlabel('g - r')
plt.ylabel('r - i')
plt.title('g-r vs r-i')
plt.legend()
plt.savefig('color_color.png')
plt.show()
plt.close()


# EAF - the two quasars are similar colors (with a difference of ~0.1 in r-i). They're much farther apart in the g-r (~0.25 apart)
#		Because they're at a similar redshift, I think that it makes sense they have similar colors. However, I'd expect them to be closer in color

plt.plot(g_r1_ext, r_i1_ext, marker='o', markersize=10, linestyle='None', label='Corrected Quasar 1')
plt.plot(g_r2_ext, r_i2_ext, marker='o', markersize=10, linestyle='None', label='Corrected Quasar 2')
plt.xlim(-0.03, 0.35)  # EAF - resizing the figure to have the same axes as the uncorrected values
plt.ylim(-0.1, 0.07)

plt.xlabel('g - r')
plt.ylabel('r - i')
plt.title('Corrected g-r vs r-i')
plt.legend()
plt.savefig('corrected_color_color.png')
plt.show()
plt.close()

# EAF - After correcting for extinction, the colors of the two quasars are closer together, although not very similar.
#		I would expect them to be even closer together than they are.


### TASK 2 (and Task 3)

# For quasar 2:
x2 = np.linspace(231.6, 241.6, 100)  # EAF - makes an array centered at 236.6 with 100 steps, each a bin of 0.1
y2 = np.linspace(-2.6, 7.4, 100)  # EAF - centered at 2.4 degrees

xm2, ym2 = np.meshgrid(x2, y2)  # EAF - X Mesh 2, Y Mesh 2
cm2 = SkyCoord(ra=xm2*u.deg, dec=ym2*u.deg)  # EAF - making SkyCoord object for each 0.1 bin
ebv_mesh2 = sfd(cm2)


plt.contourf(x2, y2, ebv_mesh2)
plt.colorbar()
plt.title('Quasar 2 Dust Map')
plt.xlabel(u'Right Ascension (bins of 0.1\N{DEGREE SIGN})')
plt.ylabel(u'Declination (bins of 0.1\N{DEGREE SIGN})')
plt.savefig('quasar2_dustmap.png')
plt.show()

# For quasar 1:
x1 = np.linspace(240.4, 253.4, 100)  # EAF - makes array with 100 bins of 0.13 degrees, centered about 246.9
y1 = np.linspace(35.8, 45.8, 100)  # EAF - centered about 40.8 deg

xm1, ym1 = np.meshgrid(x1, y1)  # EAF - X Mesh 1, Y Mesh 1
cm1 = SkyCoord(ra=xm1*u.deg, dec=ym1*u.deg)  # EAF - making SkyCoord object for each 0.1 bin
ebv_mesh1 = sfd(cm1)


plt.contourf(x1, y1, ebv_mesh1)  # EAF - plots contours of the dust
plt.colorbar()
plt.title('Quasar 1 Dust Map')
plt.xlabel(u'Right Ascension (bins of 0.13\N{DEGREE SIGN})')
plt.ylabel(u'Declination (bins of 0.1\N{DEGREE SIGN})')
plt.savefig('quasar1_dustmap.png')
plt.show()