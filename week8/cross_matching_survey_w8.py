#!/usr/local/Anaconda2024/bin/python3

################################################
#						              
#    In Class Work: Cross Matching Surveys
#		       3.13.2025			  
#	      Author: Eliza Frankel       
#                                                             
###############################################

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import glob
import os

###  FUNCTIONS  ###

def find_sweep(path, ra, dec):
    """ Uses decode_sweep_name() and a slightly modified version of is_in_box() 
    to determine what sweeps a list of objects are located in

    Parameters:
    path (str) - path to directory containing sweep fits files.
    ra (np.array or list) - list of object RAs, must be in degrees
    dec (np.array or list) - list of object Decs, must be in degrees

    Returns:
    sweeps (list) - List containing all sweeps these objects would be in.
    """

    name_list = glob.glob(path + "*fits")
    objects = np.array([ra, dec])
    sweeps = []

    for i in range(len(name_list)):
        ramin, ramax, decmin, decmax = decode_sweep_name(name_list[i])
        location = is_in_box(objects, [ramin, ramax, decmin, decmax])
        if location.any() == True:
            pathways = name_list[i].split('/')
            sweeps.append(pathways[-1])
            # EAF - the above 2 lines get just the filename, ignoring the path.

    return sweeps


# EAF - The following 2 functions (is_in_box() and decode_sweep_name() come from
#	DESI codes, https://github.com/desihub/desitarget/blob/2.4.0/py/desitarget/geomask.py#L1338-L1369 and https://github.com/desihub/desitarget/blob/2.4.0/py/desitarget/io.py#L2354-L2397, respectively.
#	This code was provided to us by Adam).
#	I wasn't sure the best way to include them in this assignment, so I opted to
#	copy them into this file.

def decode_sweep_name(sweepname, nside=None, inclusive=True, fact=4):
    """Retrieve RA/Dec edges from a full directory path to a sweep file

    Parameters
    ----------
    sweepname : :class:`str`
        Full path to a sweep file, e.g., /a/b/c/sweep-350m005-360p005.fits
    nside : :class:`int`, optional, defaults to None
        (NESTED) HEALPixel nside
    inclusive : :class:`book`, optional, defaults to ``True``
        see documentation for `healpy.query_polygon()`
    fact : :class:`int`, optional defaults to 4
        see documentation for `healpy.query_polygon()`

    Returns
    -------
    :class:`list` (if nside is None)
        A 4-entry list of the edges of the region covered by the sweeps file
        in the form [RAmin, RAmax, DECmin, DECmax]
        For the above example this would be [350., 360., -5., 5.]
    :class:`list` (if nside is not None)
        A list of HEALPixels that touch the  files at the passed `nside`
        For the above example this would be [16, 17, 18, 19]
    """
    # ADM extract just the file part of the name.
    sweepname = os.path.basename(sweepname)

    # ADM the RA/Dec edges.
    ramin, ramax = float(sweepname[6:9]), float(sweepname[14:17])
    decmin, decmax = float(sweepname[10:13]), float(sweepname[18:21])

    # ADM flip the signs on the DECs, if needed.
    if sweepname[9] == 'm':
        decmin *= -1
    if sweepname[17] == 'm':
        decmax *= -1

    if nside is None:
        return [ramin, ramax, decmin, decmax]

    pixnum = hp_in_box(nside, [ramin, ramax, decmin, decmax],
                       inclusive=inclusive, fact=fact)

    return pixnum


def is_in_box(objs, radecbox):
    """Determine which of an array of objects are inside an RA, Dec box.

    NOTE - EAF Modified this to change objs["RA"] and objs["DEC"] to 0 and 1, the
    position of these arrays in a list.

    Parameters
    ----------
    objs : :class:`~numpy.ndarray`
        An array of objects. Must include at least the columns "RA" and "DEC".
    radecbox : :class:`list`
        4-entry list of coordinates [ramin, ramax, decmin, decmax] forming the
        edges of a box in RA/Dec (degrees).

    Returns
    -------
    :class:`~numpy.ndarray`
        ``True`` for objects in the box, ``False`` for objects outside of the box.

    Notes
    -----
        - Tests the minimum RA/Dec with >= and the maximum with <
    """
    ramin, ramax, decmin, decmax = radecbox

    # ADM check for some common mistakes.
    if decmin < -90. or decmax > 90. or decmax <= decmin or ramax <= ramin:
        msg = "Strange input: [ramin, ramax, decmin, decmax] = {}".format(radecbox)
        log.critical(msg)
        raise ValueError(msg)

    # ii = ((objs["RA"] >= ramin) & (objs["RA"] < ramax)
    #       & (objs["DEC"] >= decmin) & (objs["DEC"] < decmax))

    ii = ((objs[0] >= ramin) & (objs[0] < ramax)
          & (objs[1] >= decmin) & (objs[1] < decmax))

    return ii

if __name__ == '__main__':

    ###  TASK 1  ###

    file = '/d/scratch/ASTR5160/data/first/first_08jul16.fits'
    data = Table.read(file)

    RA = data["RA"]
    Dec = data["DEC"]

    # plt.figure(figsize=(16,9))
    # plt.scatter(RA, Dec, marker='*')
    # plt.show()


    ###  TASK 2  ###

    # I chose an object at RA: 130.02403 and Dec: 39.99085
    # sdssDR9query.py returned: 130.02402626,39.99084523,15.743527,14.060556,13.161932,13.347169,14.092394,0.02004858


    ###  TASK 3  ###

    first_100_RA = RA[:100]
    first_100_Dec = Dec[:100]


    # EAF - running CLI code using os.system and saving to 'query.txt'
    #for i in range(len(first_100_RA)):
    #	ra = first_100_RA[i]
    #	dec = first_100_Dec[i]

    #	os.system(f"python sdssDR9query.py {ra} {dec} >> query.txt")
    	

    ###  TASK 4  ###

    directory = "/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/"

    #sweep1 = Table.read(directory + "sweep-350p030-360p035.fits")
    #sweep2 = Table.read(directory + "sweep-190p060-200p065.fits")

    # EAF - to read headers:
    # print(sweep1[0])
    # print(sweep2[0])

    ###  TASK 6  ###

    sweep_files = find_sweep(directory, first_100_RA, first_100_Dec)

    print(sweep_files)


