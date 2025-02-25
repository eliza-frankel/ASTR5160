#!/usr/local/Anaconda2024/bin/python3

#######################################
#						              
#       In Class Work: Mangle    
#		       2.25.2025			  
#	      Author: Eliza Frankel       
#                                                            
#######################################

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

# EAF - importing my functions from the previous class
import sys
sys.path.insert(0, '/d/users/eliza/Classes/Astro_Techniques_II/week5/')
from spherical_caps_w5 import *

import pymangle



###  CLASS NOTES  ###

# polygons - shapes represented by the intersection of multiple caps

###  TASK 1  ###

theta = 5  # unit - degrees

cap1_coordinates = SkyCoord(ra=76, dec=36, frame='icrs', unit=(u.degree, u.degree))
cap2_coordinates = SkyCoord(ra=75, dec=35, frame='icrs', unit=(u.degree, u.degree))

cap1 = spherical_cap(cap1_coordinates, theta, 1)
cap2 = spherical_cap(cap2_coordinates, theta, 1)


###  TASK 2  ###
# Create a file containing these two caps as a single polygon in the Mangle format. 
# The file should be called intersection.ply.

# spherical_cap_output takes caps (list/np.array), polygon (int), and output_file (string) as args

caps = [cap1, cap2]
spherical_cap_output([caps], 1, [0], "intersection.ply")

#  Create a second Mangle file containing these two caps as two different polygons (i.e. polygon 1
# should contain cap 1 and polygon 2 should contain cap 2). This second file should be called bothcaps.ply

caps_both = [[cap1], [cap2]]
steradians = [0, 0]
spherical_cap_output(caps_both, 2, steradians, "bothcaps.ply")


###  TASK 3  ###

m_inter = pymangle.Mangle("intersection.ply")
m_both = pymangle.Mangle("bothcaps.ply")

inter_ra, inter_dec = m_inter.genrand(10000)
both_ra, both_dec = m_both.genrand(10000)

plt.figure(figsize=(16,9))
plt.scatter(both_ra, both_dec, marker = 'o', color = 'peru', alpha = 0.5, label = "Both Caps")
plt.scatter(inter_ra, inter_dec, marker = '*', color = 'rebeccapurple', alpha = 0.5, label = "Intersection")

plt.title("Intersection vs Both Caps", fontsize=18)
plt.xlabel("Right Ascension [degrees]", fontsize=15)
plt.xticks(fontsize=12)
plt.ylabel("Declination [degrees]", fontsize=15)
plt.yticks(fontsize=12)
plt.legend(loc='upper right', fontsize=15)
plt.savefig("intersection_v_bothcaps.png", dpi=300)
plt.show()

# EAF - I think I understand why this figure looks the way it does. The masks are 
#		different due to the number of polygons. Polygons are the shapes represented 
#		by the intersection of multiple caps. When we have two separate polygon, we see all
# 		of the area that is covered by each individual cap.
#		With just one polygon, we only see where the two caps intersect, removing anything
#		outside of the shared RA and Decs.


###  TASK 4  ###

cap1 = spherical_cap(cap1_coordinates, theta, -1)
cap2 = spherical_cap(cap2_coordinates, theta, 1)
caps = [cap1, cap2]
spherical_cap_output([caps], 1, [0], "intersection_flip1.ply")
m_flip1 = pymangle.Mangle("intersection_flip1.ply")
flip1_ra, flip1_dec = m_flip1.genrand(10000)

plt.figure(figsize=(16, 9))
plt.scatter(flip1_ra, flip1_dec, marker ='o', color = 'cornflowerblue', alpha = 0.5, label = "cap1 Flip")
plt.scatter(inter_ra, inter_dec, marker = '*', color = 'rebeccapurple', alpha = 0.5, label = "Intersection")

plt.title("Intersection vs Cap1 Flip", fontsize=18)
plt.xlabel("Right Ascension [degrees]", fontsize=15)
plt.xticks(fontsize=12)
plt.ylabel("Declination [degrees]", fontsize=15)
plt.yticks(fontsize=12)
plt.legend(loc='upper right', fontsize=15)
plt.savefig("intersection_v_cap1_flip.png", dpi=300)
plt.show()


###  TASK 5  ###

cap1 = spherical_cap(cap1_coordinates, theta, 1)
cap2 = spherical_cap(cap2_coordinates, theta, -1)
caps = [cap1, cap2]

spherical_cap_output([caps], 1, [0], "intersection_flip2.ply")
m_flip2 = pymangle.Mangle("intersection_flip2.ply")
flip2_ra, flip2_dec = m_flip2.genrand(10000)

plt.figure(figsize=(16, 9))
plt.scatter(flip1_ra, flip1_dec, marker ='o', color = 'cornflowerblue', alpha = 0.5, label = "cap1 Flip")
plt.scatter(flip2_ra, flip2_dec, marker ='v', color = 'pink', alpha = 0.5, label = "cap2 Flip")
plt.scatter(inter_ra, inter_dec, marker = '*', color = 'rebeccapurple', alpha = 0.5, label = "Intersection")

plt.title("Intersection vs Cap1 & Cap2 Flip", fontsize=18)
plt.xlabel("Right Ascension [degrees]", fontsize=15)
plt.xticks(fontsize=12)
plt.ylabel("Declination [degrees]", fontsize=15)
plt.yticks(fontsize=12)
plt.legend(loc='upper right', fontsize=15)
plt.savefig("intersection_v_cap2_flip.png", dpi=300)
plt.show()


###  TASK 6  ###

cap1 = spherical_cap(cap1_coordinates, theta, -1)
cap2 = spherical_cap(cap2_coordinates, theta, -1)

caps = [cap1, cap2]
spherical_cap_output([caps], 1, [0], "intersection_flip_both.ply")
m_flip_both = pymangle.Mangle("intersection_flip_both.ply")
flip_both_ra, flip_both_dec = m_flip_both.genrand(1000000)

plt.figure(figsize=(16, 9))
plt.scatter(flip_both_ra, flip_both_dec, marker ='o', color = 'goldenrod', alpha = 0.5, label = "cap1 & cap2 flip")

plt.title("Cap 1 & Cap 2 Flip", fontsize=18)
plt.xlabel("Right Ascension [degrees]", fontsize=15)
plt.xticks(fontsize=12)
plt.ylabel("Declination [degrees]", fontsize=15)
plt.yticks(fontsize=12)
# plt.legend(loc='upper right', fontsize=15)
plt.savefig("intersection_both_flip.png", dpi=300)
plt.show()
