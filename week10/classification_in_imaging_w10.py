#!/usr/local/Anaconda2024/bin/python3

#################################################
#						              
#    In Class Work: Clasification in Imaging
#		          3.27.2025		  
#	      Author: Eliza Frankel!!!       
#                                                             
#################################################

from astropy.table import Table, vstack
from week8.cross_matching_survey_w8 import find_sweep
from astropy.coordinates import SkyCoord
import astropy.units as u
from week10.magnitude_systems_w10 import magnitude
import matplotlib.pyplot as plt
import numpy as np


def is_star(g_z, r_W1):
    """
    Parameters:
    ------------
    g_z (np.array) - g - z magnitude of object
    r_W1 (np.array) - r - W1 magnitude of object

    Returns:
    ---------
    (boolean) - True if star, False if quasar

    Note:
    -----
    This is not perfect, there are some outliers on both sides.
    However, this does account for a majority of the stars and quasar
    classifications.
    """

    b = -1
    m = 1
    x = np.linspace(-2, 6, 1000)
    boundary = m * x + b

    if r_W1 > boundary.all() and g_z < x.all():
        return False
    else: 
        return True



if __name__ == "__main__":

	###  TASK 1  ###
	directory1 = '/d/scratch/ASTR5160/week10/'
	stars = Table.read(directory1 + "stars-ra180-dec30-rad3.fits", unit_parse_strict='silent')
	quasars = Table.read(directory1 + "qsos-ra180-dec30-rad3.fits", unit_parse_strict='silent')

	# EAF - getting RA and Decs
	staRA = stars["RA"]
	starDec = stars["DEC"]

	quasaRA = quasars["RA"]
	quasarDec = quasars["DEC"]

	star_obj = SkyCoord(staRA, starDec, unit=(u.degree, u.degree))
	quasar_obj = SkyCoord(quasaRA, quasarDec, unit=(u.degree, u.degree))


	directory2 = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/'
	star_sweeps = find_sweep(directory2, staRA, starDec)
	quasar_sweeps = find_sweep(directory2, quasaRA, quasarDec)

	sweeps = star_sweeps

	# EAF - stacking all 4 sweep files into one Table and making them SkyCoord Objects
	super_sweep = vstack([Table.read(directory2 + sweep, unit_parse_strict = 'silent') for sweep in sweeps])

	super_coordinates = SkyCoord(super_sweep["RA"], super_sweep["DEC"], unit=(u.degree, u.degree))

	# EAF - taking only the star object ID for objects less than 0.5 arcsec
	star_index = star_obj.search_around_sky(super_coordinates, u.arcsecond * 0.5)[0]
	quasar_index = quasar_obj.search_around_sky(super_coordinates, u.arcsecond * 0.5)[0]

	# EAF - making a mask of the stars and quasars within the 0.5 arcsec range
	#       then getting their fluxes and the MW transmissiopn fluxes for task 2
	close_stars = super_sweep[star_index]

	star_flux_g = close_stars["FLUX_G"]
	star_flux_r = close_stars["FLUX_R"]
	star_flux_z = close_stars["FLUX_Z"]
	star_flux_W1 = close_stars["FLUX_W1"]
	star_flux_W2 = close_stars["FLUX_W2"]

	mw_g_star =  close_stars["MW_TRANSMISSION_G"]
	mw_r_star = close_stars["MW_TRANSMISSION_R"]
	mw_z_star = close_stars["MW_TRANSMISSION_Z"]
	mw_W1_star = close_stars["MW_TRANSMISSION_W1"]
	mw_W2_star =  close_stars["MW_TRANSMISSION_W2"]


	close_quasars = super_sweep[quasar_index]

	quasar_flux_g = close_quasars["FLUX_G"]
	quasar_flux_r = close_quasars["FLUX_R"]
	quasar_flux_z = close_quasars["FLUX_Z"]
	quasar_flux_W1 = close_quasars["FLUX_W1"]
	quasar_flux_W2 = close_quasars["FLUX_W2"]

	mw_g_quasar = close_quasars["MW_TRANSMISSION_G"]
	mw_r_quasar = close_quasars["MW_TRANSMISSION_R"]
	mw_z_quasar = close_quasars["MW_TRANSMISSION_Z"]
	mw_W1_quasar = close_quasars["MW_TRANSMISSION_W1"]
	mw_W2_quasar = close_quasars["MW_TRANSMISSION_W2"]


	###  TASK 2  ###

	# EAF - correcting for Galactic dust
	corrected_star_g = star_flux_g / mw_g_star
	corrected_star_r = star_flux_r / mw_r_star
	corrected_star_z = star_flux_z / mw_z_star
	corrected_star_W1 = star_flux_W1 / mw_W1_star
	corrected_star_W2 = star_flux_W2 / mw_W2_star

	corrected_quasar_g = quasar_flux_g / mw_g_quasar
	corrected_quasar_r = quasar_flux_r / mw_r_quasar
	corrected_quasar_z = quasar_flux_z / mw_z_quasar
	corrected_quasar_W1 = quasar_flux_W1 / mw_W1_quasar
	corrected_quasar_W2 = quasar_flux_W2 / mw_W2_quasar

	# EAF - converting to magnitudes, nusing magnitude() from week 10 magnitude_systems

	mag_g_star = magnitude(corrected_star_g)
	mag_r_star = magnitude(corrected_star_r)
	mag_z_star = magnitude(corrected_star_z)
	mag_W1_star = magnitude(corrected_star_W1)
	mag_W2_star = magnitude(corrected_star_W2)

	mag_g_quasar = magnitude(corrected_quasar_g)
	mag_r_quasar = magnitude(corrected_quasar_r)
	mag_z_quasar = magnitude(corrected_quasar_z)
	mag_W1_quasar = magnitude(corrected_quasar_W1)
	mag_W2_quasar = magnitude(corrected_quasar_W2)

	# # EAF - the following is a series of print statements to show all of the fluxes. I have it commented out currently.
	# print('Star Magnitudes')
	# print(mag_g_star)
	# print(mag_r_star)
	# print(mag_z_star)
	# print(mag_W1_star)
	# print(mag_W2_star)

	# print('\nQuasar Magnitudes')
	# print(mag_g_quasar)
	# print(mag_r_quasar)
	# print(mag_z_quasar)
	# print(mag_W1_quasar)
	# print(mag_W2_quasar)


	###  TASK 3  ###

	star_r_W1 = mag_r_star - mag_W1_star
	star_g_z = mag_g_star - mag_z_star

	quasar_r_W1 = mag_r_quasar - mag_W1_quasar
	quasar_g_z = mag_g_quasar - mag_z_quasar

	b = -1
	m = 1
	x = np.linspace(-2, 6, 1000)
	y = m * x + b

	plt.figure(figsize=(16,9))
	plt.scatter(star_g_z, star_r_W1, marker='*', c='goldenrod', label='Stars')
	plt.scatter(quasar_g_z, quasar_r_W1, c='slateblue', label='Quasars')
	plt.plot(x, y)

	plt.xlabel("g - z", fontsize=15)
	plt.ylabel("r - W1", fontsize=15)

	plt.title("Color Cuts", fontsize=20)
	plt.legend(fontsize=15)

	plt.show()



	# EAF - the function is is_star().

	print(is_star(quasar_g_z[0], quasar_r_W1[0]))
	print(is_star(star_g_z[0], star_r_W1[0]))