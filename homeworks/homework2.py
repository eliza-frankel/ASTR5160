#!/usr/local/Anaconda2024/bin/python

#########################################
#                                                             
#    Homework 2 - Areas on the Sphere    
#                3.28.2025                        
#        Author: Eliza Frankel       
#                                                             
#########################################

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
import argparse


# EAF - defining constants and conversion factors

pi = np.pi
to_degrees = 180 / pi
to_radians = 1 / to_degrees  # EAF - this is redundant, but makes code more readable


def area(ra_min, ra_max, dec_min, dec_max):
    """ Finds the area of a field in the sky
    Parameters:
    -----------
    ra_min (float) - minimum RA value, in degrees
    ra_max (float) - maximum RA value, in degrees
    dec_min (float) - minimum Dec value, in degrees
    dec_max (float) - maximum Dec value, in degrees

    Returns:
    --------
    area (float) - area of the field, in square degrees
    """

    dec_diff = np.sin(dec_max * to_radians) - np.sin(dec_min * to_radians)
    ra_diff = ra_max - ra_min

    area = to_degrees * ra_diff * dec_diff

    return area



def aitoff_plot(ra_min, ra_max, directory, overlap=None):

    """ Finds the area of a field in the sky
    Parameters:
    -----------
    ra_min (float) - minimum RA value, in degrees
    ra_max (float) - maximum RA value, in degrees
    directory (string) - location to save outputted figure
    overlap (integer) - chooses which field region to plot the overlapping, 
                        randomly generated points. If no parameter, defaults to None.

    """

    # EAF - Arbitrarily determine 4 regions of declinations

    dec_mins = [60, 15, -30, -75]
    dec_maxs = [90, 45, 0, -45]


    area1 = area(ra_min, ra_max, dec_mins[0], dec_maxs[0])
    area2 = area(ra_min, ra_max, dec_mins[1], dec_maxs[1])
    area3 = area(ra_min, ra_max, dec_mins[2], dec_maxs[2])
    area4 = area(ra_min, ra_max, dec_mins[3], dec_maxs[3])

    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111, projection='aitoff')


    # EAF - plotting overlapping randomly generated points
    if overlap != None:

        regions = populate(ra_min, ra_max, dec_mins[overlap - 1], dec_maxs[overlap - 1])
        ras = np.array([coord[0] for coord in regions])
        decs = np.array([coord[1] for coord in regions])

        ax.scatter(ras * to_radians, decs * to_radians)

    # EAF - plotting the different spots on an aitoff projection.

    ax.axhspan(dec_mins[3] * to_radians, dec_maxs[3] * to_radians, (ra_min + 180) / 360.,\
         (ra_max + 180) / 360., edgecolor='goldenrod', facecolor='None', linewidth=4, label=area4)    
    ax.axhspan(dec_mins[2] * to_radians, dec_maxs[2] * to_radians, (ra_min + 180) / 360.,\
         (ra_max + 180) / 360., edgecolor='slateblue', facecolor='None', linewidth=4, label=area3)
    ax.axhspan(dec_mins[1] * to_radians, dec_maxs[1] * to_radians, (ra_min + 180) / 360.,\
         (ra_max + 180) / 360., edgecolor='forestgreen', facecolor='None', linewidth=4, label=area2)
    ax.axhspan(dec_mins[0] * to_radians, dec_maxs[0] * to_radians, (ra_min + 180) / 360.,\
         (ra_max + 180) / 360., edgecolor='pink', facecolor='None', linewidth=4, label=area1)


    # rand_ra = np.array(360. * (random(1000000)))
    # rand_dec= np.array((180 / pi) * np.arcsin(1. - random(1000000) * 2.))
    # ax.scatter(rand_ra, rand_dec)


    ax.grid(True)
    plt.title('Aitoff Projection\n')
    plt.legend(title="Area")
    plt.savefig(directory + '/aitoff_plot.png', dpi=300)
    # plt.show()



def populate(ra_min, ra_max, dec_min, dec_max):
    """ Finds the area of a field in the sky
    Parameters:
    -----------
    ra_min (float) - minimum RA value of field, in degrees
    ra_max (float) - maximum RA value of field, in degrees
    dec_min (float) - minimum Dec value of field, in degrees
    dec_max (float) - maximum Dec value of field, in degrees

    Returns:
    --------
    coordinates (np.array or list) - list of coordinates that are located in the field provided
                                     Will be in the format of [[ra0, dec0], [ra1, dec1], ...]
    """

    # EAF - randomly generate set of 1 million points on surface of sphere
    rand_ra = np.array(360. * (random(1000000)) - 180)
    rand_dec= np.array((180 / pi) * np.arcsin(1. - random(1000000) * 2.))

    overlapping = [[rand_ra[i], rand_dec[i]] for i in range(1000000) if rand_ra[i] \
        < ra_max and rand_ra[i] > ra_min and rand_dec[i] < dec_max and rand_dec[i] > dec_min]

    return overlapping





if __name__ == "__main__":

    # # EAF - All checks from the homework assignment
    #         Currently commented out. Feel free to uncomment to check!


    # check = area(0, 360, 0, 90)
    # halfsphere = 41252.96/2

    # if round(check, 2) == halfsphere:
    #     print("area() is working properly. The area of this region is {}".format(round(check, 2)))
    # else:
    #     print("Keep working!")


    # # EAF - checking that random generation has the right number of random points 
    # #       relative to entire sphere

    # overlap = populate(-10, 25, 60, 90)
    # fraction_overlap = len(overlap) / 1000000

    # area_of_sphere = 4 * pi * 180 * 180 / pi / pi
    # area_of_region = area(-10, 25, 60, 90)
    # fraction_area = area_of_region / area_of_sphere

    # print("The number of points in the region is {}.".format(len(overlap)))
    # print("The fraction of these points compared to the total 1,000,000 is {}."\
    #     .format(fraction_overlap))
    # print("The fraction of area in the region compared to the area of the sphere is {}."\
    #     .format(round(fraction_area, 6)))
    # print("These are practically the same fraction!")



    # EAF - calling argparse and running the functions

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=str,\
         help="Please input the path to the directory you'd like to save your figure.")

    args = parser.parse_args()
    directory = args.directory

    aitoff_plot(-10, 25, directory)

    print("Done! Check {} for the final plot.".format(directory))

    
