from astropy import units as u
from astropy.coordinates import SkyCoord   # EAF - for converting RA and Dec to decimals
from astropy.coordinates import EarthLocation  # EAF - used for Earth lat and long, not RA and Dec
from astropy.coordinates import AltAz
from astropy.time import Time  # EAF - obtaining MJD and JD
import numpy as np

def dec_to_degrees(dec):
	""" Converts declination from ° ' " to degrees

	Parameters:
	dec - class: 'str'
		must be in the format (_d_m_s)

	Returns:
	degrees - class: 'float'
	"""
	d, temp = dec.split('d')  # EAF - separates the degrees from the rest of the dec. Degrees saved in 'd', everything else in 'temp'
	d_sign = d[0]  # EAF - gets the sign of the degree
	d = d[1:]  # EAF - gets the value of degree
	m, s = temp.split('m')  # EAF - separates arcmin and arcsec, placing them in 'm' and 's', respectively
	s = s[:-1]  # EAF - removes the last character from the string, s.

	if d_sign == '-':  # EAF - this if else statement assigns the sign of the first term in the product, depending on the sign of the given declination
		sign = -1
	else:
		sign = 1

	degrees = sign * (int(d) + (int(m) / int(60)) + (int(s) / 3600))
	return degrees

def ra_to_degrees(ra):
	""" Converts RA from h, m, s to degrees

	Parameters:
	ra - class: 'str'
		must be in the format (_h_m_s)

	Returns:
	degrees - class: 'float'
	"""

	h, temp = ra.split('h')  # # EAF - separates hours from the rest of RA. Hours stored in 'h', everything else in 'temp'
	m, s = temp.split('m')  # EAF - separates arcmin and arcsec, placing them in 'm' and 's', respectively
	s = s[:-1]  # EAF - removes the last character from the string, s.

	degrees = 15 * (int(h) + (int(m) / int(60)) + (int(s) / 3600))
	return degrees

# Task 2

coord_RA = '23h12m11s'
coord_Dec = '-40d12m13s'  # EAF - I'm not sure if this is the format you are requesting, but through documentation couldn't determine what to do

coordinates = SkyCoord(ra=coord_RA, dec=coord_Dec, frame='icrs')  # EAF makes SkyCoord obj using data 
ra_deg = coordinates.ra.degree  # EAF - converts RA to degrees
dec_deg = coordinates.dec.degree  # EAF - converts dec to degrees

# print(coordinates)  # EAF - This prints the (ra, dec) in units of degrees

print('Using SkyCoords, the RA and Declination were determined to be:\n (' + str(ra_deg) + ', ' + str(dec_deg))
# EAF - the following line uses the functions dec_to_degree() and ra_to_degree() written above
print('Using the equations from class, the RA and Declination were determined to be:\n (' + str(ra_to_degrees(coord_RA)) + ', ' + str(dec_to_degrees(coord_Dec)))
print('These are the same.\n')

# Task 3

now = Time.now()
MJD = now.mjd  # EAF - determines modified Julian date based on current time
JD = now.jd  # EAF - determines Julian date based on current time

calc_MJD = JD - 2400000.5  # EAF - calculates modified JD based on equation from notes

print('Using Time.now(), the Julian date was determined to be ' + str(JD) +' and modifid Julian date to be ' + str(MJD))
print('Using the equation from class, the modified JD ws determined to be', calc_MJD)
print('These are the same.\n')

# Task 4

similar_dates = np.arange(MJD-10, MJD+10, 1)  # EAF - creates an array of Julian dates +/- today's MJD
print("20 dates around today's MJD are:")
print(similar_dates)
print()

# Task 5

# EAF - Longitude is used as RA and latitude as Dec in the equitorial coordinate system
WIRO_longitude = '+105h58m33s'  # EAF - the format is wrong bc degrees but for my function it works
WIRO_long_degree = ra_to_degrees(WIRO_longitude)  # EAF - this is because the factor of 15 is not in the degrees to decimal degrees eqn for longitude
	 # EAF - The above converts longitude using RA to Degrees function above
WIRO_long_degree /= -15

WIRO_latitude = '+41d5m49s'
WIRO_lat_degree = dec_to_degrees(WIRO_latitude) 

WIRO_height = 2943 # [m]

WIRO = EarthLocation(lat=WIRO_lat_degree*u.deg, lon=WIRO_long_degree*u.deg, height=WIRO_height)

obj_RA = '12h0m0s'
obj_RA_deg = ra_to_degrees(obj_RA)

obj_dec = '+30d0m0s'
obj_dec_deg = dec_to_degrees(obj_dec)

obj_coordinates = SkyCoord(ra=obj_RA_deg*u.deg, dec=obj_dec_deg*u.deg, frame='icrs')  # EAF - makes SkyCoord obj using data 

utc_offset = -7 * u.hour   # EAF - Mountain Standard Time
t_tonight = Time('2025-1-30 23:00:00') - utc_offset  # EAF - setting the time to night 
t_month = Time('2025-3-2 23:00:00') - utc_offset

obj_altaz = obj_coordinates.transform_to(AltAz(obstime=[t_tonight, t_month], location=WIRO))

airmass = obj_altaz.secz
print('The airmass today, 1.30.2025 at 11:00 pm is', airmass[0])
print('The airmass in 1 month, 3.2.2025', airmass[1])