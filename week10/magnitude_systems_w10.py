#!/usr/local/Anaconda2024/bin/python3

################################################
#						              
#    In Class Work: Cross Matching Surveys
#		          3.25.2025			  
#	      Author: Eliza Frankel       
#                                                             
################################################

from week8.cross_matching_survey_w8 import find_sweep
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np


def magnitude(flux):
    """
    Parameters:
    ------------
    flux (float) - flux from sweep file, in nanomaggies

    Returns:
    --------
    magnitude
    """

    return 22.5 - 2.5 * np.log10(flux)


###  NOTES  ###

"""
- Flux in the sweep files are in nanomaggies (3631*10^9) Jy or 10^9 f_0
- Flux (in sweep files) to magnitudes - m = 22.5 - 2.5log10(flux)
"""

if __name__ == "__main__":

    ###  TASK 1  ###
    """
    Consider UBVRI for the (standard) star PG1633+099A,
    the conversion between UBVRI and ugriz, and the SDSS
    Navigate Tool — all on the WyoCourses links page
    • Use the UBVRI to ugriz transformations to show that the
      g and z magnitudes displayed for PG1633+099A in the
      SDSS Navigate Tool are near the expected values
    """

    #  According to the links, the following is the UBVRI for PG1633+099A
    V = 15.256
    B_V = 0.873
    U_B = 0.320
    V_R = 0.505
    R_I = 0.511

    # Using the final section for stars with Rc-Ic < 1.15:

    u_g = 1.28 * U_B + 1.13
    g_r = 1.02 * B_V - 0.22
    r_i = 0.91 * R_I - 0.20
    r_z = 1.72 * R_I - 0.41
    g = V + 0.60 * B_V - 0.12
    r = V - 0.42 * B_V + 0.11

    sdss_ugriz = ['17.30', '15.70', '15.19', '14.71', '14.55']
    calc_ugriz = [u_g + g, g, r, r - r_i, r - r_z]

    print(" --------------------------------------")
    print("\t\tugriz Values \n --------------------------------------")
    print("\tCalculated \t SDSS Navigate \n --------------------------------------")
    print("   u\t  {} \t     {}".format(calc_ugriz[0], sdss_ugriz[0]))
    print("   g\t  {} \t     {}".format(calc_ugriz[1], sdss_ugriz[1]))
    print("   r\t  {} \t     {}".format(calc_ugriz[2], sdss_ugriz[2]))
    print("   i\t  {} \t     {}".format(calc_ugriz[3], sdss_ugriz[3]))
    print("   z\t  {} \t     {}".format(calc_ugriz[4], sdss_ugriz[4]))
    print(" --------------------------------------")

    print("These magnitudes are similar to the magnitudes shown in the SDSS Navigate Tool.\n")


    ###  TASK 2  ###

    path = "/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/"
    ra = 248.85833  # 16:35:26
    dec = 9.79806  # +09:47:53

    # EAF - using function from last week to find the sweep file!
    sweep = find_sweep(path, [ra], [dec])

    data = Table.read(path + sweep[0], unit_parse_strict='silent')

    dataRA = data["RA"]
    dataDec = data["DEC"]

    # EAF - creating SkyCoord objects for PG1633+099A and for all objects in sweep file
    PG1633_099A = SkyCoord('16:35:26', '+09:47:53', unit=(u.hourangle, u.degree))
    sweep_objects = SkyCoord(dataRA, dataDec, unit=(u.degree, u.degree))

    # EAF - finding the index of PG1633+099A from the sweep file
    #       separation() shows the difference between the coordinates of PG1633+099A
    #       and all the sweep object coordinates
    separations = PG1633_099A.separation(sweep_objects)
    index = separations.argmin()

    # EAF - Getting grz fluxes, and WISE magnitudes!

    flux_g = data["FLUX_G"][index]
    flux_r = data["FLUX_R"][index]
    flux_z = data["FLUX_Z"][index]
    flux_W1 = data["FLUX_W1"][index]
    flux_W2 = data["FLUX_W2"][index]
    flux_W3 = data["FLUX_W3"][index]
    flux_W4 = data["FLUX_W4"][index]

    mag_g = magnitude(flux_g)
    mag_r = magnitude(flux_r)
    mag_z = magnitude(flux_z)
    mag_W1 = magnitude(flux_W1)
    mag_W2 = magnitude(flux_W2)
    mag_W3 = magnitude(flux_W3)
    mag_W4 = magnitude(flux_W4)

    print(" --------------------------------------")
    print("\t\tugriz Values \n --------------------------------------")
    print("\tFrom Sweep \t SDSS Navigate \n --------------------------------------")
    print("   g\t  {} \t     {}".format(round(mag_g, 2), sdss_ugriz[1]))
    print("   r\t  {} \t     {}".format(round(mag_r, 2), sdss_ugriz[2]))
    print("   z\t  {} \t     {}".format(round(mag_z, 2), sdss_ugriz[4]))
    print(" --------------------------------------")

    print("These magnitudes from sweep files are similar to the magnitudes shown in the SDSS Navigate Tool.\n")


    print(" ------------------------------")
    print("\tWISE Magnitudes \n ------------------------------")
    print("   W1\t  {} ".format(mag_W1))
    print("   W2\t  {} ".format(mag_W2))
    print("   W3\t  {} ".format(mag_W3))
    print("   W4\t  {} ".format(mag_W4))

    print(" ------------------------------")

    print("PG1633+099A was detected in W1, W2, and W3 bands. It was not detected in W4 (as seen above)")