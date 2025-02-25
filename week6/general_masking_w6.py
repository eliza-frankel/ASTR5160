#!/usr/local/Anaconda2024/bin/python3

#######################################
#						              
#       In Class Work: General Masking    
#		    Due: 2.27.2025			  
#	      Author: Eliza Frankel       
#                                                            
#######################################

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

# EAF - importing my functions from week 5
import sys
sys.path.insert(0, '/d/users/eliza/Classes/Astro_Techniques_II/week5/')
from spherical_caps_w5 import *

import pymangle

import time

# EAF - creating a variable for convertine between radians and degrees
rad_to_deg = 180 / np.pi

###  TASK 1  ###

# RA in hourangles
ra1 = '5h'
ra2 = '6h'

# RA in degrees 
# EAF - Converts from hourangles to degrees (1 hrangle = 15 deg)
ra1_deg = float(ra1[:-1]) * 15
ra2_deg = float(ra2[:-1]) * 15

ra1_cap = RA_cap_bound(ra1, 1)
ra2_cap = RA_cap_bound(ra2, -1)

# Declination in degrees
dec1 = 50
dec2 = 60

dec1_cap = Dec_cap_bound(dec1, 1)
dec2_cap = Dec_cap_bound(dec2, -1)

steradians = ((ra2_deg / rad_to_deg) - (ra1_deg / rad_to_deg)) \
	 * (np.sin(dec2 / rad_to_deg) - np.sin(dec1 / rad_to_deg))

# EAF - calling in my function spherical_cap_output from spherical_caps_w5
caps = [ra1_cap, ra2_cap, dec1_cap, dec2_cap]
spherical_cap_output([caps], 1, [steradians], "pythontask.ply")


###  TASK 2  ###

ra1_t2 = '11h'
ra2_t2 = '12h'

ra1_deg_t2 = float(ra1_t2[:-1]) * 15
ra2_deg_t2 = float(ra2_t2[:-1]) * 15

ra1_cap_t2 = RA_cap_bound(ra1_t2, 1)
ra2_cap_t2 = RA_cap_bound(ra2_t2, -1)

# Declination in degrees
dec1_t2 = 60
dec2_t2 = 70

dec1_cap_t2 = Dec_cap_bound(dec1_t2,1 )
dec2_cap_t2 = Dec_cap_bound(dec2_t2, -1)

steradians_t2 = ((ra2_deg_t2 / rad_to_deg) - (ra1_deg_t2 / rad_to_deg)) \
	 * (np.sin(dec2_t2 / rad_to_deg) - np.sin(dec1_t2 / rad_to_deg))

caps_t2 = [ra1_cap_t2, ra2_cap_t2, dec1_cap_t2, dec2_cap_t2]

total_caps = [caps, caps_t2]
total_steradians = [steradians, steradians_t2]
spherical_cap_output(total_caps, 2, total_steradians, "pythontask.ply")


###  TASK 3  ###

# EAF - generating 1 million objects distributed over the sphere

x = time.time()
ra = 360. * (random(1000000))
dec = (180 / np.pi) * np.arcsin(1. - random(1000000) * 2.)

###  TASK 4  ###

mask = pymangle.Mangle("pythontask.ply")
within_polygons = mask.contains(ra, dec)


plt.figure(figsize = (16,9))
plt.scatter(ra, dec, marker = 'o', color = "skyblue")
plt.scatter(ra[within_polygons], dec[within_polygons], marker = '*', color = 'darkblue', label = "Points within polygons")

plt.title("Points in the lat-lon Rectangular Polygons", fontsize=18)
plt.xlabel("Right Ascension [degrees]", fontsize=15)
plt.xticks(fontsize=12)
plt.ylabel("Declination [degrees]", fontsize=15)
plt.yticks(fontsize=12)
plt.legend(loc='upper right', fontsize=15)

plt.savefig("contains.png", dpi=300)

y = time.time()
print("Time task 4:", y-x)

plt.show()


### TASK 5  ###

x = time.time()

ra_random, dec_random = mask.genrand(10000)

plt.figure(figsize = (16,9))
plt.scatter(ra, dec, marker = 'o', color = "palegoldenrod")
plt.scatter(ra_random, dec_random, marker = '*', color = 'goldenrod', label = "Randomly Generated points")

plt.title("GenRand", fontsize=18)
plt.xlabel("Right Ascension [degrees]", fontsize=15)
plt.xticks(fontsize=12)
plt.ylabel("Declination [degrees]", fontsize=15)
plt.yticks(fontsize=12)
plt.legend(loc='upper right', fontsize=15)

plt.savefig("genrand.png", dpi=300)

y = time.time()
print("Time task 5:", y-x)

plt.show()

# EAF - Generating random numbers takes slightly longer than using mask.contains. 
#		I used the time.time() function to determine this (as seen by my print statement),
#		However visually it feels longer to calculate too.
#
#		The relative area of the randomly generated lat-lon rectangle is consistent
#		with the area of the .contains() points. Visually, they look identical.