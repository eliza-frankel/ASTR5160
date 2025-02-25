#!/usr/local/Anaconda2024/bin/python3.11

#######################################
#						              
#    In Class Work: Spherical Caps    
#		       2.20.2025			  
#	      Author: Eliza Frankel       
#                                                             
#######################################

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
from astropy.coordinates import SkyCoord
from astropy import units as u

###  TASK 1  ###

def RA_cap_bound(RA):
    """ Converts RA from hms to cartesian coordinates
    
    Parameters:
    ra (string): the right ascension, specifically in units of hour angles or degrees

    Returns:
    np array consisting of:
	    x (float or int): x component of RA
	    y (float or int): y component of RA
	    z (float or int): z component of RA
	    l_cos0 (float or int): Width of RA component (1-h)
    """

    declination='0d0m0s'

    # EAF - Determining whether RA is in hourangles or degrees
    if 'h' in RA:
    	coord = SkyCoord(RA, declination, frame='icrs', unit=(u.hourangle, u.degree))
    else:
    	coord = SkyCoord(RA, declination, frame='icrs', unit=(u.degree, u.degree))
    

	# EAF - I can't think of a better way to add the 90 degrees to the RA than this 
    coord_degrees = SkyCoord(ra = coord.ra.degree + 90, dec = coord.dec.degree, unit=(u.degree, u.degree))
    coord_degrees.representation_type = 'cartesian'

    x = float(coord_degrees.x)
    y = float(coord_degrees.y)
    z = float(coord_degrees.z)
    l_cos0 = 1

    return np.array([x, y, z, l_cos0])


###  TASK 2  ###

def Dec_cap_bound(declination):
	""" Converts Declination from hms to cartesian coordinates
    
    Parameters:
    declination (string): the declination, specifically in units of degrees

    Returns:
    np array consisting of:
	    x (float or int): x component of Dec
	    y (float or int): y component of Dec
	    z (float or int): z component of Dec
	    l_sin0 (float or int): Width of Dec component (1-h)
    """
	coord = SkyCoord(ra=0, dec=90, frame='icrs', unit=(u.degree, u.degree))
	coord.representation_type = 'cartesian'

	x = float(coord.x)
	y = float(coord.y)
	z = float(coord.z)
	# EAF - convert declination to radians
	l_sin0 = 1 - np.sin(declination * (np.pi / 180))

	return np.array([x, y, z, l_sin0])




###  TASK 3  ###

def spherical_cap(coordinates, theta):
	""" Creates the vector 4-array for the spherical cap at (ra, dec)

	Parameters:
	coordinates (SkyCoord object): the coordinates of the ra and dec
	theta (float or int): should be in degrees

	Returns:
	np.array consisting of:
		x (float): x component
		y (float): y component
		z (float): z component
		l-cos0 (float): height of the cap
	"""

	coordinates.representation_type = 'cartesian'

	x = float(coordinates.x)
	y = float(coordinates.y)
	z = float(coordinates.z)
	l_cos0 = 1 - np.cos(theta * (np.pi / 180))

	return np.array([x, y, z, l_cos0])


###  TASK 4  ###

# EAF - the following is a function to write the results of the above 3 functions
#		to a text file. 
# 		It has been modified from the previous submission to allow for differing numbers
#		of polygons and caps to be used.

def spherical_cap_output(caps, polygons, outfile_name):
	""" Writes output of RA_cap_bound(), Dec_cap_bound(), and spherical_cap() to a file

	Parameters:
	caps (list or np.array): list of arrays from the ra, dec, and spherical cap functions
	polygons (int): the number of polygons in the output file
	outfile_name (string): a string containing the name you'd like to save your file as
		- Should also include the suffix (ie .txt, .ply, etc.)

	Notes:
		Creates a text file with the desired results
		* number of caps must be divisible by number of polygons in this iteration of the function
			- Currently thinking of something better to do this.
	"""

	with open(outfile_name, 'w') as file:
		index = 0
		file.write(f'{polygons} polygons\n')
		# EAF - creates the number of polygons in the file
		for i in range(polygons):
			file.write(f'polygon {i+1} ( {len(caps) // polygons} caps, 1 weight, 0 pixel, 0 str):\n')	
			
			# EAF - For each cap, write the data in this specified format
			for j in range(len(caps) // polygons):
				file.write(f'\t{caps[index][0]:19.16f} {caps[index][1]:19.16f} {caps[index][2]:19.16f} {caps[index][3]:19.16f}\n')
				index += 1


if __name__ == "__main__":
	# Task 1
	RA = RA_cap_bound('5h')
	print(RA)

	# Task 2
	Dec = Dec_cap_bound(36)
	print(Dec)

	# Task 3
	coords = SkyCoord(ra='5h', dec='36d', frame='icrs', unit=(u.hourangle, u.degree))
	sphere = spherical_cap(coords, 1)
	print(sphere)

	# Task 4
	caps = np.array([RA, Dec, sphere])
	spherical_cap_output(caps, 1, 'cap_output.txt')
