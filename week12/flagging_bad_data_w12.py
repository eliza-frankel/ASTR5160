#!/usr/local/Anaconda2024/bin/python3

################################################
#						              
#    In Class Work: Flagging Bad Data
#		          4.8.2025		  
#	      Author: Eliza Frankel       
#                                                             
################################################

import numpy as np
from week8.cross_matching_survey_w8 import find_sweep
from week10.magnitude_systems_w10 import magnitude
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table, vstack
import astropy.units as u


def good_objects(sweep_file, BITMASK, flag_exp):
	""" Finds which sources are not bad data

	Parameters:
	-----------
	sweep_file (astropy Table) - Table of the sweep file in question
	BITMASK (string) - the bitmask desired.
	flag_exp (int) - integer corresponding to the flag exponent. Is determined
		by what you're looking for (ie saturated is 1, interpolated is 2).
		found: https://www.legacysurvey.org/dr9/bitmasks/#allmask-x-anymask-x

	Returns:
	--------
	goodobjs (list) - only the good objects

	Notes:
	------
	This comes from the 'Flagging Bad Data.pdf' file from class
	"""

	flag = 2**flag_exp
	ii = (sweep_file[BITMASK] & flag) == 0
	goodobjs = sweep_file[ii]

	return goodobjs



def mask_data(desired_object, goodobjs, BITMASK):
	""" Updates the data table using good_objects()

	Parameters:
	-----------
	desired_object (SkyCoord object) - the SkyCoord obj of desired target
	goodobjs (Table) - results of good_objects()
	BITMASK (string) - desired BTIMASK

	Returns:
	---------
	info (int) - the BITMASK column for the desired object
	"""
	good_sweep_objects = SkyCoord(goodobjs["RA"], goodobjs["DEC"],\
		unit=(u.degree, u.degree))

	separation = desired_object.separation(good_sweep_objects)

	good_object_index = separation.argmin()
	good_object_information = goodobjs[good_object_index]

	info = good_object_information[BITMASK]

	return good_object_information




if __name__ == "__main__":

	### TASK 1 ###
	""" Find the closest object in the sweep files to the RA/Dec
	location (α, δ) = (188.53667 o, 21.04572o )

		- What is the TYPE of this object? Is its light profile
		“stellar” (like a star) or extended (like a galaxy)?
	"""

	RA = 188.53667
	DEC = 21.04572

	obj = SkyCoord(RA, DEC, unit=(u.degree, u.degree))

	path = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/'

	# EAF - determining which sweep file this object is in!
	sweep = find_sweep(path, RA, DEC)

	# EAF - reading in the data from the sweep file.
	data = Table.read(path + sweep[0], unit_parse_strict='silent')

	sweep_RA = data["RA"]
	sweep_DEC = data["DEC"]
	sweep_TYPE = data["TYPE"]

	sweep_objects = SkyCoord(sweep_RA, sweep_DEC, unit=(u.degree, u.degree))

	separation = obj.separation(sweep_objects)

	object_index = separation.argmin()

	# EAF - determining the object's type.

	obj_type = sweep_TYPE[object_index]

	print("The object at RA: {}, Dec: {} is of type {}, which means it is exponential"\
		.format(RA, DEC, obj_type))

	# EAF - I believe that the EXP morphological model is extended, like a galaxy.



	### TASK 2 ###
	"""
	Consider the ALLMASK_G , ALLMASK_R and ALLMASK_Z bitmasks for this object.
	In how many bands out of the three (g, r, z) is the object saturated in
	all exposures?

		- Find this object in the Legacy Surveys Sky Viewer (see
		the links page). Is it saturated? Is it even a galaxy?
	"""
	desired_object = data[object_index]


	ALLMASK_G = desired_object["ALLMASK_G"]
	ALLMASK_R = desired_object["ALLMASK_R"]
	ALLMASK_Z = desired_object["ALLMASK_Z"]

	# print(ALLMASK_G, ALLMASK_R, ALLMASK_Z)

	# EAF - determining which masks are saturated in this object

	goodobjs_G = good_objects(data, "ALLMASK_G", 1)
	info_G = mask_data(obj, goodobjs_G, "ALLMASK_G")

	goodobjs_R = good_objects(data, "ALLMASK_R", 1)
	info_R = mask_data(obj, goodobjs_R, "ALLMASK_R")

	goodobjs_Z = good_objects(data, "ALLMASK_Z", 1)
	info_Z = mask_data(obj, goodobjs_Z, "ALLMASK_Z")

	# EAF - checking to see if my object and the object found via mask_data
	#		are the same:

	# print(desired_object["RA"], desired_object["DEC"])
	# print(info_G["RA"], info_G["DEC"])
	
	obj_new = SkyCoord(188.53400825178812, 21.045324773607916, unit=(u.degree, u.degree))
	separation = obj_new.separation(sweep_objects)
	i = separation.argmin()

	print("Are the desired object and the object found via mask_data the same?",\
		object_index == i)


	print("\nThe band is saturated in all three exposures:")
	print("\t g:", desired_object["ALLMASK_G"])
	print("\t r:", desired_object["ALLMASK_R"])
	print("\t z:", desired_object["ALLMASK_Z"])

	# Looking at the Legacy Sky Viewer, I don't believe that this object
	# is very saturated, although there is some saturation. - This leads me
	# to the question of what did I do wrong to get no saturation?
	# It also doesn't look like a galalaxy. According to the 'Explore' tab, 
	# it has the flag 'PSF', so it would be a star.

	

	### TASK 3 ###

	directory = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/'

	central_object = SkyCoord(180, 30, unit=(u.degree, u.degree))

	# # EAF - find what sweep file this object is in and confirm the RA and DEC ranges
	# ra_range = [177, 180, 183]
	# dec_range = [27, 30, 33]

	# sweeps = find_sweep(directory, ra_range, dec_range)

	# EAF - using the same sweep files from week10 classification_in_imaging_w10.py
	# EAF - I tried using my find_sweep() function to figure out the files, but
	# 		it was giving me only two sweep files, not the 4 that are necessary.

	sweeps = ['sweep-170p025-180p030.fits', 'sweep-170p030-180p035.fits',\
	'sweep-180p025-190p030.fits', 'sweep-180p030-190p035.fits']

	# EAF - making one large sweep file with the data from all (both) files
	#		containing these objects.
	super_sweep = vstack([Table.read(directory + sweep, unit_parse_strict = 'silent')\
		for sweep in sweeps])


	# EAF - determining which objects in the sweep files are PSFs
	psf_mask = (super_sweep["TYPE"] == b"PSF")
	sweep_psf = super_sweep[psf_mask]

	# EAF - converting flux to magnitude and to narrow down filters
	R_flux = sweep_psf["FLUX_R"]
	R_mag = magnitude(R_flux)

	mag_filter = (R_mag < 20)

	sweep_psf_mag = sweep_psf[mag_filter]

	super_sweep_SkyCoords = SkyCoord(sweep_psf_mag["RA"], sweep_psf_mag["DEC"],\
		unit=(u.deg, u.deg))
	
	# super_sweep_SkyCoords = SkyCoord(sweep_psf["RA"], sweep_psf["DEC"],\
	# 	unit=(u.deg, u.deg))


	# EAF - narrowing down objects within 3 degrees of our coordinates.
	separations_pt3 = central_object.separation(super_sweep_SkyCoords)
	i_less_than_3deg = separations_pt3 < (3 * u.degree)

	psfobjs = super_sweep_SkyCoords[i_less_than_3deg]

	print('num psf objects:', len(psfobjs))


	quasar_file = '/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits'

	quasar_data = Table.read(quasar_file, unit_parse_strict = 'silent')

	quasars = SkyCoord(quasar_data["RA"], quasar_data["DEC"], unit=(u.deg, u.deg))

	# EAF matching my objects to the quasars:
	index = psfobjs.search_around_sky(quasars, u.arcsecond * 1)[0]

	# EAF - all quasars that match with the psf objects
	qsos = quasar_data[index]

	print("There are {} objects that are quasars and have r magnitude < 20"\
		.format(len(qsos)))

