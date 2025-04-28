#!/usr/local/Anaconda2024/bin/python

#############################################
#                                                             
#    Homework 3 - Object Classification  
#           	 	5.2.2025                        
#    		  Author: Eliza Frankel       
#   	                                                          
#############################################

import numpy as np
import matplotlib.pyplot as plt
from time import time, sleep
import argparse

from astropy.coordinates import search_around_sky, SkyCoord
from astropy.table import Table, vstack
from astropy.io import fits
import astropy.units as u

from week8.cross_matching_survey_w8 import find_sweep
from week10.magnitude_systems_w10 import magnitude
from homework3 import read_file
import warnings



"""
EAF
    The following 3 functions (crossmatching, determining_colorcuts, and testing)
    are functions I wrote to determine the colorcuts used in splendid_function,
    and to test the functionality of said function. I am leaving these here to show
    my process, but they are not needed for the actual assignment grading.
"""



def crossmatching(sweep_path, quasar_RAs, quasar_DECs):
    """ Determines which objects are in the sweep files, and gets all necessary 
    information about them.

    Parameters:
    -----------
    sweep_path (string) - path to sweep files
    quasar_RA (list or np.array) - list of RAs for desired objects
    quasar_DEC (list or np.array) - list of DECs for desired objects
    
    Returns:
    --------
    quasar_sweeps (Astropy Table) - table containing all necessary information about the 
    quasar objects

    Notes:
    ------
    This function is used to determine the initial color cuts
    """

    # EAF - determining which sweep files the quasars are located in
    sweeps = find_sweep(path_to_sweeps, quasar_RAs, quasar_Decs)
	
    super_sweep = vstack([Table.read(sweep_path + sweep,\
        unit_parse_strict = 'silent') for sweep in sweeps])
    
    # EAF - first checking to see if these sources have a magnitude of r < 19
    super_sweep_R_flux = magnitude(super_sweep["FLUX_R"]) < 19
    super_sweep_criteria = super_sweep[super_sweep_R_flux]
    
    super_SkyCoords = SkyCoord(super_sweep_criteria["RA"], super_sweep_criteria["DEC"],\
        unit=(u.degree, u.degree))

    
    quasars = SkyCoord(quasar_RAs, quasar_DECs, unit=(u.degree, u.degree))

    # EAF - search_around_sky[] returns the index for these objects
    indices = quasars.search_around_sky(super_SkyCoords,\
        seplimit = u.arcsecond * 0.5)[0]

    print(len(indices))
    
    # EAF - creating the table of only quasars with r-magnitude < 19
    quasar_sweeps = super_sweep[indices]

    return quasar_sweeps



def determining_colorcuts(sweep_path, quasar_RAs, quasar_DECs):
    """ Determining the colorcuts to use in the homework

    Parameters:
    -----------
    sweep_path (string) - path to sweep files
    quasar_RA (list or np.array) - list of RAs for desired objects
    quasar_DEC (list or np.array) - list of DECs for desired objects

    Note:
    -----
    Only used to determine the initial color-color cuts!
    """

    sweep = crossmatching(path_to_sweeps, quasaRAs, quasar_Decs)

    g_mag = magnitude(sweep["FLUX_G"])
    r_mag = magnitude(sweep["FLUX_R"])
    z_mag = magnitude(sweep["FLUX_Z"])
    w1_mag = magnitude(sweep["FLUX_W1"])
    w2_mag = magnitude(sweep["FLUX_W2"])
    w3_mag = magnitude(sweep["FLUX_W3"])
    w4_mag = magnitude(sweep["FLUX_W4"])

    # EAF - finding the boundaries for the combinations I'll use to hardcode

    cc_cut_name = ['g-r', 'g-z', 'g-w2', 'g-w4', 'r-z', 'r-w1', 'r-w3',\
        'r-w4', 'z-w2', 'z-w4','w1-w4', 'w2-w3']
    cc_cut = [g_mag - r_mag, g_mag - z_mag, g_mag - w2_mag, g_mag - w4_mag,\
        r_mag - w1_mag, r_mag - w3_mag, r_mag - w4_mag, z_mag - w2_mag,\
        z_mag - w4_mag, w1_mag - w4_mag, w2_mag - w3_mag]

    with open('hw4_colorcuts.txt', 'w') as f:

        for i in range(len(cc_cut_name)):
            f.write("{} max: {} min: {}\n"\
                .format(cc_cut_name[i], max(cc_cut[i]),min(cc_cut[i])))



def testing():
    """ Testing function to ensure splendid_function works.
    The following code is copied from my script classification_in_imaging_w10.py.
    It is part of the in class tasks for March 27th.

    Returns:
    --------
    (astropy Table of star and quasar objects)
    """

    directory1 = '/d/scratch/ASTR5160/week10/'
    stars = Table.read(directory1 + "stars-ra180-dec30-rad3.fits", unit_parse_strict='silent')
    quasars = Table.read(directory1 + "qsos-ra180-dec30-rad3.fits", unit_parse_strict='silent')


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
    close_quasars = super_sweep[quasar_index]


    #EAF - This is the only thing I have changed, in order to create a Table of
    #       known star and quasar objects.

    objects = vstack([close_stars, close_quasars])

    return objects


##################################################################################

def splendid_function(data_table):
    """ Uses predetermined colorcuts to decide if an object is a quasar or not

    Parameters:
    -----------
    data_table (Astropy Table) - Table containing information from sweep files
        about objects

    Returns:
    --------
    (np.array of Booleans) - an array with an entry for each row of data_table;
        will be True if the object is determined to be a quasar, false otherwise

    Notes:
    ------
     -  Uses colorcuts predetermined (and hardcoded in) from using 
        crossmatching() and determining_colorcuts() above. For the sake of
        timing, these color cuts have been hardcoded in, however in a more
        optimized/generic code with more possible training targets, these may
        be softcoded.
     -  The predetermined colorcuts using the above functions are also located
        in the file "hw4_colorcuts.txt" which is created when 
        determining_colorcuts() is run.
    """

    # EAF - ignoring warnings, specifically for 
    #       'RuntimeWarning: invalid value encountered in subtract'
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    # EAF - creating arrays of the magnitudes of each object in data_table
    g_mag = magnitude(data_table["FLUX_G"])
    r_mag = magnitude(data_table["FLUX_R"])
    z_mag = magnitude(data_table["FLUX_Z"])
    w1_mag = magnitude(data_table["FLUX_W1"])
    w2_mag = magnitude(data_table["FLUX_W2"])
    w3_mag = magnitude(data_table["FLUX_W3"])
    w4_mag = magnitude(data_table["FLUX_W4"])

    # EAF - all combos of colorcuts 
    cc_cut = [g_mag - r_mag, g_mag - z_mag, g_mag - w2_mag, g_mag - w4_mag,\
        r_mag - w1_mag, r_mag - w3_mag, r_mag - w4_mag, z_mag - w2_mag,\
        z_mag - w4_mag, w1_mag - w4_mag, w2_mag - w3_mag]

    # EAF - the maximum and minimum values were found from the function
    #       determining_colorcuts. They are copied into lists for making
    #       the if statement more compact.
    # EAF - these are in the same order as the color color cuts in cc_cut

    cc_max = [5.9665985107421875, 7.963407516479492, 10.792217254638672,\
        12.070793151855469, 5.220605850219727, 8.312885284423828,\
        11.722053527832031, 3.648937225341797, 11.858993530273438,\
        9.923791885375977, 8.087688446044922]

    cc_min = [-1.0673351287841797, -2.6758460998535156, -3.4638748168945312,\
        1.1107330322265625, -3.699237823486328, 0.2900714874267578,\
        0.4204483032226562, -4.749082565307617, -0.1328983306884765,\
        -0.4219303131103515, -1.5249576568603516]

    # EAF - create a Boolean array. Entries will be changed to True if a quasar.
    object_type_quasar = np.full(len(g_mag), False)

    # EAF - long if statement to determine if object is considered quasar or not
    for i in range(len(g_mag)):
        if ( (cc_cut[0][i] < cc_max[0] and cc_cut[0][i] > cc_min[0])\
                and (cc_cut[1][i] < cc_max[1] and cc_cut[1][i] > cc_min[1])\
                and (cc_cut[2][i] < cc_max[2] and cc_cut[2][i] > cc_min[2])\
                and (cc_cut[3][i] < cc_max[3] and cc_cut[3][i] > cc_min[3])\
                and (cc_cut[4][i] < cc_max[4] and cc_cut[4][i] > cc_min[4])\
                and (cc_cut[5][i] < cc_max[5] and cc_cut[5][i] > cc_min[5])\
                and (cc_cut[6][i] < cc_max[6] and cc_cut[6][i] > cc_min[6])\
                and (cc_cut[7][i] < cc_max[7] and cc_cut[7][i] > cc_min[7])\
                and (cc_cut[8][i] < cc_max[8] and cc_cut[8][i] > cc_min[8])\
                and (cc_cut[9][i] < cc_max[9] and cc_cut[9][i] > cc_min[9])\
                and (cc_cut[10][i] < cc_max[10] and cc_cut[10][i] > cc_min[10])\
                ):

            object_type_quasar[i] = True
        # else:
        #     object_type_quasar[i] = False

    return object_type_quasar


    



if __name__ == "__main__":

    # # EAF - Code to determine the colorcuts & test splendid_function
    # #       Not necessary for homework, but included!

    # quasar_file = '/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits'

    # # EAF - reading in the quasar file as an astropy table
    # data = Table.read(quasar_file, unit_parse_strict='silent')
    # quasaRAs = data["RA"]
    # quasar_Decs = data["DEC"]

    # path_to_sweeps = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/'
    # determining_colorcuts(path_to_sweeps, quasaRAs, quasar_Decs)

    # # EAF - this is the code to test if splendid_function works
    data = testing()

    ###############################################################################

    # EAF - setting up argparse to accept directory to sweepfile
    help_msg = """
    Please input the path to the directory of the sweep file you wish to look at
    and determine if objects are quasars or not.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=str, help=help_msg)

    args = parser.parse_args()
    directory = args.directory

    data = Table.read(directory, unit_parse_strict='silent')

    start = time()

    # EAF - calling splending_results to compare data to the color cuts
    splendid_results = splendid_function(data)

    # EAF - counting the number of objects that are quasars. array.sum counts
    #       all 'True' values in the array.
    count = splendid_results.sum()

    print(f"The number of objects classified as quasars is {count}")

    end = time()

    print("It took {} seconds to find your quasars!".format(round(end - start, 3)))

