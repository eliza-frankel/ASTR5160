#!/usr/local/Anaconda2024/bin/python3.11

#########################################
#										#
#   Homework 1: Survey Observations		#
#       Author: Eliza Frankel			#
#		Due: February 28, 2025			#
#										#
#########################################

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy import units as u
from astropy.time import Time
import argparse

# EAF - reading in the data file, and writing it to the list 'quasarfile'

filename = '/d/scratch/ASTR5160/week4/HW1quasarfile.txt'
quasarfile = np.genfromtxt(filename, dtype=str)

# EAF - separating RA and Dec
ra = [i[0:9] for i in quasarfile]
dec = [i[9:] for i in quasarfile]

# EAF - putting the RA and Dec into proper formatting
ra_hourangle = np.array([f'{i[0:2]}h{i[2:4]}m{i[4:9]}' for i in quasarfile])
dec_degree = np.array([f'{i[9:12]}d{i[12:14]}m{i[14:]}' for i in quasarfile])



# EAF - makes SkyCoord object for each quasar
all_coords = SkyCoord(ra=ra_hourangle, dec=dec_degree, frame='icrs')#, unit=(u.hourangle, u.degree))

# # Making time offset for MST to UTC
utc_offset = -7 * u.hour

# EAF - starting times at Jan 1, 2025, and the months all through 2025
# EAF - UTC is 7 hours ahead of MST, so 11:00 pm MST is 6:00 am UTC
jan1 = Time('2025-1-1 23:00:00') - utc_offset
year = str(jan1)[0:4]

mjd_jan1 = jan1.mjd
times = np.arange(mjd_jan1, mjd_jan1+365, 1)
times_iso = [Time(times[i], format = 'mjd').iso for i in range(len(times))]
times_iso = np.array(times_iso)



# EAF - The following chunk of code gets the desired month of observation from the command line
parser = argparse.ArgumentParser()
# EAF - choice gives argparse a set of valid numbers. Anything outside of this range is invalid.
parser.add_argument("month", type=int, help="What month would you like to observe in? (Between 1-12)", \
	choices=range(1, 13))
args = parser.parse_args()

MONTH = args.month

# EAF - the following is a list with the indices when the months begin in the quasarfile
if int(year) % 4 == 0:
	# EAF - This is a leap year. Beginning indices must account for this
	beginning_index = [-1, 30, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
elif int(year) % 4 != 0:
	beginning_index = [-1, 30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333]

# EAF - determining the indices (in quasarfile) for the desired observation month 
# Months with 31 days
if MONTH in [1, 3, 5, 7, 8, 10, 12]:
	# EAF - these are the months with 31 days (index should be + 1 at the start to account for time difference
	indices = np.arange(beginning_index[MONTH - 1] + 1, beginning_index[MONTH - 1] + 32, 1)
	
elif MONTH in [4, 6, 9, 11]:
	# EAF - the months with 30 days
	indices = np.arange(beginning_index[MONTH - 1] + 1, beginning_index[MONTH - 1] + 31, 1)
else:
	# EAF - Determining if it's a leap year or not
	if int(year) % 4 == 0:
		# EAF - Leap Year	
		indices = np.arange(beginning_index[MONTH - 1] + 1, beginning_index[MONTH - 1] + 30, 1)
	elif int(year) % 4 != 0:
		indices = np.arange(beginning_index[MONTH - 1] + 1, beginning_index[MONTH - 1] + 29, 1)



# EAF - list of the lowest airmass for each day 
#		and the index of the quasar with lowest airmass for that day
KPNO = EarthLocation.of_site("kpno")
lowest_airmass = []
quasar_index = []

for i in range(len(times_iso)):
# for i in range(1):
	alt_az = AltAz(location=KPNO, obstime = times_iso[i])
	# EAF - the following creates an airmass for every single quasar on just one day
	coords = all_coords.transform_to(alt_az)

	# EAF - finding the  maximum altitude (which will have the minimum airmass)
	max_alt_min_airmass_index = np.where(coords.alt.deg == max(coords.alt.deg))[0][0]
	
	# EAF - gets the airmass of the maximum altitude quasar and appends to lowest_airmass
	airmass = coords[max_alt_min_airmass_index].secz
	lowest_airmass.append(airmass)

	# EAF - recording coordinates of lowest airmass quasar
	quasar_index.append(max_alt_min_airmass_index)

lowest_airmass = np.array(lowest_airmass)
quasar_name = np.array(quasarfile[quasar_index])

# EAF - gets the SkyCoord object for only the desired quasars
quasar_coords = all_coords[quasar_index]

# EAF - choosing the airmass, coordinates, ra, and dec of only the desired quasars
desired_airmass = lowest_airmass[indices]
desired_quasar = quasar_name[indices]
desired_ra = quasar_coords[indices].ra.deg
desired_dec = quasar_coords[indices].dec.deg
desired_dates = times_iso[indices]
desired_dates = Time(desired_dates) + utc_offset


# EAF - printing the results
print('-------------------------------------------------------------------------------------------------------------')
print('\tDate\t\t\tQuasar Coordinates \t RA (°)\t\t Dec(°)\t\t Airmass')
print('\t\t\t\t    (hms.ss°\'") \t\t\t')
print('-------------------------------------------------------------------------------------------------------------')

for i in range(len(indices)):
	print('{} \t{}\t{} \t{} \t{}'.format(Time(desired_dates[i]), desired_quasar[i], \
		round(desired_ra[i], 5), round(desired_dec[i], 5), desired_airmass[i]))

print('-------------------------------------------------------------------------------------------------------------')