#!/usr/local/Anaconda2024/bin/python3

#######################################
#						              
#    In Class Work: Coding Week
#		       3.6.2025			  
#	      Author: Eliza Frankel       
#                                                             
#######################################

from astropy.coordinates import SkyCoord
import astropy.units as u
from week5.spherical_caps_w5 import spherical_cap


# EAF - taking this from week 6 mangle_w6.py file

theta = 5  # unit - degrees
cap1_coordinates = SkyCoord(ra=76, dec=36, frame='icrs', unit=(u.degree, u.degree))

print(cap1_coordinates)