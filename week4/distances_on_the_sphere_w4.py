#!/usr/local/Anaconda2024/bin/python3.11

#####################################################
#													#
#   Week 4 in class work: Distance on the Sphere	#
#       Author: Eliza Frankel						#
#													#
#####################################################

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

### Notes

# EAF - X.separation(Y) finds the separation between respective x, y, z points.
#		Can also take a single point and find all objects within a specified 
#		distance to it.

ra1, dec1 = [8., 9., 10.]*u.degree, [0., 0., 0.]*u.degree
ra2, dec2 = [18., 19., 21.]*u.degree, [0., 0., 0.]*u.degree
c1 = SkyCoord(ra1, dec1, frame='icrs')
c2 = SkyCoord(ra2, dec2, frame='icrs')
# EAF - prints in the form 00d00m00s instead of __ deg
# print(c1.separation(c2)) 

# EAF - separation of all points to one specified point (within 1.5 deg here)
#		Returns 'True' if coordinates are within the specified range,
#		Returns 'False' if coordinates are outside of it.
ra1, dec1 = [20.]*u.degree, [0.]*u.degree
ra2, dec2 = [18., 19., 21.]*u.degree, [0., 0., 0.]*u.degree
c1 = SkyCoord(ra1, dec1, frame='icrs')
c2 = SkyCoord(ra2, dec2, frame='icrs')
ii = c1.separation(c2) < 1.5*u.degree



# EAF - Y.search_around_sky(X) finds the separation between EVERY point
# same ra1, dec1 and ra2, dec2 as above example

id1, id2, d2, d3 = c2.search_around_sky(c1, 9.1*u.degree)
# EAF - id1 and id2 are indexes in the ra and dec arrays.
#		Tells matches bw 2 sets of points
# print(id1, id2)
# print(ra1[id1], dec1[id1])


### TASK 1

def separation_angle(obj1, obj2):
	""" Returns the dot product of two coordinates
	
	Parameters:
	obj1 - class: SkyCoord object
		must be in cartesian coordinates
	obj2 - class: SkyCoord object
		must be in cartesian coordinates

	Returns:
	z_angle - class: np.array
		The angle between the object and the zenith

	Note:
		The equation for the dot product is a.b = |a||b|cos(z_angle),
		so rearranging this gives the equation for z_angle
	"""

	a = np.array([obj1.x, obj1.y, obj1.z])
	b = np.array([obj2.x, obj2.y, obj2.z])

	a_dot_b = np.sum(a*b)
	abs_a = (np.sum(a**2)) ** (1/2)
	abs_b = (np.sum(b**2)) ** (1/2)

	denom = abs_a * abs_b
	
	z_angle = np.arccos(a_dot_b / denom)
	# EAF - converting from radians to degrees
	z_angle *= (180 / np.pi)

	return z_angle

ra1 = 263.75*u.degree
dec1 = -17.9*u.degree
coordinates1 = SkyCoord(ra1, dec1, frame='icrs')


ra2 = '20h24m59.9s'
dec2 = '10d6m0s'
coordinates2 = SkyCoord(ra2, dec2, frame='icrs')


separation_cartesian = coordinates1.separation(coordinates2)
separation_cartesian.deg
print("The separation using '.separation()' is {}".format(separation_cartesian))

coordinates1 = coordinates1.cartesian
coordinates2 = coordinates2.cartesian

separation_dot_product = separation_angle(coordinates1, coordinates2)
print("The separation using a dot product is {} deg".format(separation_dot_product))
print("These are the same values.\n")


### TASK 2

# EAF - creates an array of 100 random numbers between 0 and 1, twice for RA and twice for Dec
rand_ra_1 = np.random.random(100)
# # EAF - brings the random RA values between 2h and 3h
rand_ra_1 += 2

rand_ra_2 = np.random.random(100)
rand_ra_2 += 2

rand_dec_1 = np.random.random(100)
# EAF - makes the rand_dec span 4 degrees
rand_dec_1 *= 4
# EAF - moves the dec between the range -2 to 2
rand_dec_1 -= 2

rand_dec_2 = np.random.random(100)
rand_dec_2 *= 4
rand_dec_2 -= 2

# EAF - creating SkyCoord objects of my random RA and Decs
rand_obj1 = SkyCoord(rand_ra_1, rand_dec_1, unit=(u.hourangle, u.degree))
rand_obj2 = SkyCoord(rand_ra_2, rand_dec_2, unit=(u.hourangle, u.degree))

plt.figure(figsize=(16,9))
plt.scatter(rand_ra_1, rand_dec_1, color='mediumslateblue', marker='*', s=100)
plt.scatter(rand_ra_2, rand_dec_2, color='lightpink', marker='p', s=100)
# plt.tight_layout()
plt.xlabel('Right Ascension [h]')
plt.ylabel(u'Declination [\N{DEGREE SIGN}]')
plt.title('Randomly Generaly RA and Decs')
plt.savefig('random_ra_dec.png', dpi=300)
# plt.show()

### TASK 3

# EAF - setting the range of points near each other
range = (10*u.degree) / 60
id1, id2, d2, d3 = rand_obj2.search_around_sky(rand_obj1, range)
# EAF - id1 and id2 are indexes in the ra and dec arrays.
#		Tells matches bw 2 sets of points

plt.scatter(rand_ra_1[id1], rand_dec_1[id1], color='goldenrod', label="Within 10'")
plt.scatter(rand_ra_2[id2], rand_dec_2[id2], color='goldenrod')

plt.legend(loc='upper right')
plt.savefig('random_ra_dec_within_10arcmin.png', dpi=300)
plt.show()
plt.close()


### Task 4 and Task 5

all_ra = np.concatenate([rand_ra_1, rand_ra_2])
all_dec = np.concatenate([rand_dec_1, rand_dec_2])
all = SkyCoord(all_ra, all_dec, unit=(u.hourangle, u.deg))

# EAF - creating SkyCoord object for the plate centered at a specific point
plate_radius = 1.8*u.degree
plate_ra = '2h20m5s'
plate_dec = '-0d6m12s'
plate = SkyCoord(plate_ra, plate_dec, frame='icrs', unit=(u.hourangle, u.deg))

ii = plate.separation(all) < 1.8*u.degree

plt.figure(figsize=(16,9))
plt.scatter(all_ra, all_dec, color='darkolivegreen', label='All Random Points')
plt.scatter(all_ra[ii], all_dec[ii], color='goldenrod', label=u"< 1.8\N{DEGREE SIGN} from (2h20m5s, -0\N{DEGREE SIGN}6'12\")")

# plt.tight_layout()
plt.title(u'Points overlapping with the 1.8\N{DEGREE SIGN} Plate')
plt.xlabel('Right Ascension [h]')
plt.ylabel(u'Declination [\N{DEGREE SIGN}]')
plt.legend(loc='upper right')
plt.savefig('overlapping_separation.png', dpi=300)
plt.show() 


### ELIZA - CHECK TASK 3 AGAIN BECAUSE OVERLAPPING GUYS ARE WEIRD???
