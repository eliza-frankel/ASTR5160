#!/usr/local/Anaconda2024/bin/python

#########################################
#                                                             
#    Homework 3 - Data Mining    
#            4.11.2025                        
#      Author: Eliza Frankel       
#                                                             
#########################################

import numpy as np
import matplotlib.pyplot as plt
from time import time, sleep
import argparse

from astropy.coordinates import search_around_sky, SkyCoord
from astropy.io import fits
import astropy.units as u

from week8.cross_matching_survey_w8 import find_sweep
from week10.magnitude_systems_w10 import magnitude
import sdssDR9query



def read_file(path, parameters):
	""" Reads a FITS file and returns an array of coordinates

	Parameters:
	-----------
	path (string) - path to directory with FIRST FITS files
	parameters (list of strings) - list of parameters of the files to return
		i.e. "RA", "DEC", "OBJID", "FLUX_R", etc.

	Returns:
	--------
	information (list of lists) - all objects in the specified columns
	"""

	data = fits.open(path, memmap=True)[1].data

	information = [data[param] for param in parameters]

	return information



def determine_sources(RAs, DECs, region_RA, region_DEC, radius, units=u.degree):
	""" Determines what objects lie in a specified radius
	
	Parameters:
	-----------
	RAs (np.array) - array of RAs to compare to an object/region
	DECs (np.array) - array of Decs to compare to an object/region
	region_RA (float) - RA at the center of a specified region
	region_DEC (float) - Dec at the center of a specified region
	radius (float) - Radius of the circular region
	units (astropy unit) - unit of the RA and Decs. If it is not specified,
		units will default to u.degree.

	Returns:
	--------
	survey_objects (SkyCoord) - SkyCoord objects within the desired range
	"""

	objects = SkyCoord(RAs, DECs, unit=(units, units))
	region_center = SkyCoord(region_RA, region_DEC, unit=(units, units))

	# EAF - finding separation between the FIRST objects and the region_center.
	#		In a list says True if it is within the radius and False otherwise.
	mask = region_center.separation(objects) < (radius * units)
	
	# EAF - Creates SkyCoord objects within the radius of the survey
	survey_objects = objects[mask]

	return survey_objects



def fits_criteria(path, objects):
	""" Checks that the r-band magnitude is < 22 and W1 - W2 mag is > 0.5.
		Also gets and returns all Legacy fluxes.

	Parameters:
	-----------
	path (string) - the path to a sweep file
	objects (np.array of SkyCoord) - list of SkyCoord objects to search for
	
	Returns:
	--------
	final_RA (list) - list of RAs for objects that fit criteria
	final_DEC (list) - list of RAs for objects that fit criteria
	"""

	# EAF - first reads the sweep file and gets necessary columns

	data = read_file(path, ["RA", "DEC", "FLUX_R", "FLUX_W1", "FLUX_W2"])
	RA = np.array(data[0])
	DEC = np.array(data[1])

	all_sources = SkyCoord(RA, DEC, unit=(u.degree, u.degree))

	# EAF - determine if object is in this sweep file
	index = objects.search_around_sky(all_sources, seplimit = u.arcsecond * 1.0)[0]

	updated_RA = RA[index]
	updated_DEC = DEC[index]

	# EAF - separating fluxes
	R = data[2][index]
	W1 = data[3][index]
	W2 = data[4][index]

	# EAF - converting from fluxes to magnitudes. 
	R_mag = magnitude(R)
	W1_mag = magnitude(W1)
	W2_mag = magnitude(W2)

	# EAF - checking the criteria (r-band magnitude < 22 and W1 - W2 > 0.5)
	mask = (R_mag < 22) & ((W1_mag - W2_mag) > 0.5)
	print()

	final_RA = updated_RA[mask]
	final_DEC = updated_DEC[mask]

	return final_RA, final_DEC



def query(RA, DEC):
	""" Makes an SDSS query for desired obects. Based on sdssDR9query.py 
		provided by Adam, but modified to write information to a list, rather
		than printing to the commandline or an external file.

	Parameters:
	-----------
	RA (np.array) - list of RAs of object to query
	DEC (np.array) - list of Decs of object to query

	Returns:
	--------
	ras (list) - a list of source RAs
	decs (list) - a list of source Decs
	u_mags (list) - a list of source u-band magnitudes
	i_mags (list) - a list of source i-band magnitudes
	"""

	print("Please hold! Currently querying your sources.\n")

	# ADM initialize the query.
	qry = sdssDR9query.sdssQuery()

    # ADM the query to be executed. You can substitute any query, here!
    # EAF - I have modified this line to just get ra, dec, u, and i mag.
	query = ["""SELECT top 1 ra,dec,u,i,GNOE.distance*60 FROM PhotoObj as PT
	JOIN dbo.fGetNearbyObjEq({},{},0.02) as GNOE
	on PT.objID = GNOE.objID ORDER BY GNOE.distance""".format(str(RA[i]), str(DEC[i]))\
		for i in range(len(RA))]

	objects = []

	for i in range(len(query)):
	    # ADM execute the query.
		qry.query = query[i]
		for line in qry.executeQuery():
			result = line.strip()

		# ADM NEVER remove this line! It won't speed up your code, it will
	    # ADM merely overwhelm the SDSS server (a denial-of-service attack)!
		sleep(1)
		
		# ADM the server returns a byte-type string. Convert it to a string.
		objects.append(result.decode())


	# EAF - determining which objects were found in the query versus not.
	objects = np.array(objects)
	valid_objects = objects != 'No objects have been found'
	final_list = objects[valid_objects]

	# EAF - splitting the results of the query into RA, DEC, u-magnitude 
	#		and i-magnitude
	split = np.array([o.split(',') for o in final_list])
	
	ras = split[:,0]
	decs = split[:,1]
	u_mags = split[:,2]
	i_mags = split[:,3]

	return ras, decs, u_mags, i_mags



def source_fluxes(RAs, DECs, u_mags, i_mags):
	""" Queries the FIRST objects, determines the brightest object, 
		and finds the nine fluxes for the brightest object.

	Parameters:
	-----------
	RAs (np.array) - list of RAs
	DECs (np.array) - list of Decs
	u_mags (np.array) - list of u-band magnitudes
	i_mags (np.array) - list of i-band magnitudes

	Returns:
	--------
	ubrite1 (int) - index of brightest object (in u-band)
	fluxes (list) - list of the 9 available fluxes for the brightest source
			SDSS u, i and Legacy g, r, z, W1, W2, W3, W4
	"""

	# EAF - determine the index of the brightest source
	ubrite1 = np.argmin(u_mags)

	# EAF - convert from magnitude to flux
	u_flux = 10 ** ((22.5 - float(u_mags[ubrite1])) / 2.5)
	i_flux = 10 ** ((22.5 - float(i_mags[ubrite1])) / 2.5)

	# EAF - determining the Legacy Survey fluxes!
	# EAF - first finding which sweep file:
	path = '/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/'
	sweep = find_sweep(path, float(RAs[ubrite1]), float(DECs[ubrite1]))

	data = read_file(path + sweep[0], ["RA", "DEC", "FLUX_G", "FLUX_R", "FLUX_Z", "FLUX_W1",\
		"FLUX_W2", "FLUX_W3", "FLUX_W4"])

	# EAF - determining which object in the sweep file is the brightest.
	brightest_object = SkyCoord(RAs[ubrite1], DECs[ubrite1], unit=(u.degree, u.degree))
	data_objs = SkyCoord(data[0], data[1], unit=(u.degree, u.degree))
	separations = brightest_object.separation(data_objs)
	bright_obj_index = separations.argmin()

	# EAF - making a list of the fluxes for this bright object.
	fluxes = [u_flux, i_flux, data[2][bright_obj_index], data[3][bright_obj_index], \
		data[4][bright_obj_index], data[5][bright_obj_index], data[6][bright_obj_index], \
		data[7][bright_obj_index], data[8][bright_obj_index]]

	return ubrite1, fluxes



def plotting(fluxes):
	""" Plots the desired flux vs wavelength graph

	Parameters:
	-----------
	fluxes (list) - list of the 9 fluxes for the object ubrite1
		the fluxes are 'u, i, g, r, z, W1, W2, W3, W4' in that order

	Notes:
	------
	By default, this function shows a figure on the screen.
	Optionally, it can be saved to a figure, the command currently
		being commented out.
	"""

	plt.figure(figsize=(16,9))
	plt.scatter(3543E-10, fluxes[0], marker = 'o', s=50, label = 'u-band')
	plt.scatter(4770E-10, fluxes[2], marker = 'o', s=50, label = 'g-band')
	plt.scatter(6231E-10, fluxes[3], marker = 'o', s=50, label = 'r-band')
	plt.scatter(7625E-10, fluxes[1], marker = 'o', s=50, label = 'i-band')
	plt.scatter(9134E-10, fluxes[4], marker = 'o', s=50, label = 'z-band')

	plt.scatter(3.4E-6, fluxes[5], marker = 'o', s=50, label = 'W1')
	plt.scatter(4.6E-6, fluxes[6], marker = 'o', s=50, label = 'W2')
	plt.scatter(12E-6, fluxes[7], marker = 'o', s=50, label = 'W3')
	plt.scatter(22E-6, fluxes[8], marker = 'o', s=50, label = "W4")

	plt.title("Fluxes of ubrite1 as a function of wavelength")
	plt.xlabel("Wavelength (meters)")
	plt.ylabel("Fluxes")

	plt.legend()

	# plt.savefig("flux_v_wavelength.png", dpi=300)
	plt.show()



def main(RA, DEC, radius, sweep_path):
	""" Joins all of the above functions into one function with the goal of
	making "if __name__ == '__main__'" less bulky

	Parameters:
	-----------
	RA (float or int) - central RA of the object/region desired *
	DEC (float or int) - dentral DEC of the object/region desired *
	radius (float or int) - radius to look for objects *
	sweep_path (str) - Path to sweep files

	Returns:
	--------
	time (float) - time (in seconds) it takes to run this function 
		(AKA the bulk of the code)

	Notes:
	------
	* RA, DEC, and radius must all be in degrees.
	"""

	start = time()

	path_FIRST = '/d/scratch/ASTR5160/data/first/first_08jul16.fits'
	FIRST_RA, FIRST_DEC = read_file(path_FIRST, ["RA", "DEC"])

	# EAF - determining which FIRST sources are in the survey of radius 3 degrees
	#		and separating their RAs and Decs - in degrees. 

	survey_FIRST_objects = determine_sources(FIRST_RA, FIRST_DEC, RA,\
		DEC, radius)
	survey_RA = survey_FIRST_objects.ra.degree
	survey_DEC = survey_FIRST_objects.dec.degree

	# EAF - Crossmatching to north sweep survey and determining which fit
	#		the criteria for fluxes, too.

	sweeps = find_sweep(sweep_path, survey_RA, survey_DEC)
	coord = [fits_criteria(sweep_path + sweep, survey_FIRST_objects) for sweep in sweeps]

	ra = np.concatenate([coord[i][0] for i in range(len(sweeps))])
	dec = np.concatenate([coord[i][1] for i in range(len(sweeps))])

	print("The total number of sources found is {}\n".format(len(ra)))

	# EAF - Finding U and I band mags from the queries
	source_RAs, source_DECs, u_mags, i_mags = query(ra, dec)
	print("The number of sources with a match in the SDSS database is {}\n"\
		.format(len(source_RAs)))

	ubrite1, fluxes = source_fluxes(source_RAs, source_DECs, u_mags, i_mags)
	# fluxes = u, i, g, r, z, W1, W2, W3, W4

	end = time()

	plotting(fluxes)

	print("The brightest object can be found at RA: {}, Dec: {} \n".\
	format(source_RAs[ubrite1], source_DECs[ubrite1]))

	return end - start



if __name__ == "__main__":
	
	help_message_path = """ 
	Please input the path to the sweep files you'd like to search.
	These files will be used to search the desired area of the sky.
	"""

	help_message_coordinates = """
	Please input the coordinates (RA, Dec, search radius) all in degrees.
	If nothing is inputted, these values will default to
	RA = 163, Dec = 50, radius = 3.
	"""

	# EAF - setting up argparse, ensuring there is a help message
	parser = argparse.ArgumentParser(add_help=True)
	parser.add_argument("path", help=help_message_path)
	parser.add_argument("-coordinates", help=help_message_coordinates, default='163 50 3.0')

	args = parser.parse_args()
	sweep_path = args.path
	coordinates = args.coordinates

	coordinates = coordinates.split()
	center_RA = float(coordinates[0])
	center_Dec = float(coordinates[1])
	radius = float(coordinates[2])
	

	# EAF - calling main() which runs all above functions
	time = main(center_RA, center_Dec, radius, sweep_path)

	text = """
	ubrite is found at a redshift of 1.035. According to SDSS Viewer,
	this object has a type 'STAR'. However, via the explore tab,
	the sourceType and class are both 'QSO'. 
	ubrite1 was determined by finding the lowest u-band magnitude, which
	would be the brightest object from the list. ubrite1 being so bright at
	a redshift of z = 1.035 means the classification of 'QSO' makes a lot
	of sense.
	The spectrum plotted matches that of a quasar (plotted flux vs 
	wavelength), again confirming that the identification of ubrite1 
	makes sense. Additionally,using the magnitudes of ubrite1 creates a 
	spectra very similar to that on the SDSS Viewer.
	Based on these results, the classification of quasar makes sense for 
	this object.
	"""

	print(text)

	print("This code took {} seconds to run.".format(round(time, 3)))

	