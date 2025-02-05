from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
import numpy as np
import matplotlib.pyplot as plt

# EAF - writing functions to convert from RA and dec to cartesian coordinates

def spherical_to_cartesian(coords):
    """ Converts RA and DEC to cartesian coordinates
    Parameters:
    coords - class: SkyCoords

    Returns:
    x - class: float
      x coordinate in cartesian
    y - class: float
      y coordinate in cartesian
    z - class: float
      z coordinate in cartesian
    """
    spherical = coords.to_string('decimal')  # EAF - this creates a string with the coordinates in the form 'ra   dec'
    ra, dec = spherical.split()  # EAF - splitting the string into ra and dec, respectively
    ra = float(ra)  # EAF - converting the strings to floats so they are usable mathematically
    dec = float(dec)

    # EAF - the following three lines calculate cartesian coordinates using cosine and sin functions. np.sin and np.cos take radians, so np.radians converts from degrees to radians
    x = np.cos(np.radians(ra)) * np.cos(np.radians(dec))
    y = np.sin(np.radians(ra)) * np.cos(np.radians(dec))
    z = np.sin(np.radians(dec))
    return x, y, z

### TASK 1

RA = '23h12m11s'  # EAF - the RA and dec are the same that I used in last week' in class assignment
DEC = '-40d12m13s'

coordinates = SkyCoord(ra=RA, dec=DEC, frame='icrs')  # EAF - ICRS is the reference system
c = spherical_to_cartesian(coordinates)  # EAF - using equations from class to calculate x, y, x
coordinates.representation_type = 'cartesian' # EAF - can be written in one line (using SkyCoord), or separated like this.

print("Using SkyCoord representation_type, the cartesian coordinates are: x = {}, y = {}, and z = {}".format(round(float(coordinates.x),8), round(float(coordinates.y),8), round(float(coordinates.z),8)))

print("Using the equations from class, the cartesian corrdinates are: x = {}, y = {}, z = {}".format(round(c[0], 8), round(c[1],8), round(c[2],8)))

print("These are the same values.\n")


### TASK 2

coord_galactic = SkyCoord(frame=Galactic, l=0*u.deg, b=0*u.deg)  # EAF - creates SkyCoord of the galactic center, (0deg, 0deg), using the Galactic frame keyword
coord_icrs = coord_galactic.icrs  # EAF - converts from galactic coordinates to RA and Dec. Can also use coord_galactic.transform_to('icrs')
print("The coordinate for the center of the galaxy are,")
print(coord_galactic)
print(coord_icrs)

print("In hours and minutes the RA of the center of the galaxy is {}".format(coord_icrs.ra.hms))
print("In degrees and minutes, the Dec of the center of the galaxy is {}".format(coord_icrs.dec.dms))

# EAF - The center of the galaxy falls in the constellation Sagittarius. It falls near the edge of the constellation, with Sag's RA being +19 05.94 and the galactic center has a RA of +17 45.37.


### TASK 3

lar_dec = 40
lar_ra = np.arange(0, 360, 360/365.25)  # EAF - creates an array of 0 to 360 degrees in steps of the degrees of 1 day (360 / 365.25)

laramie_icrs = SkyCoord(ra=lar_ra, dec=lar_dec, frame='icrs', unit='deg')  # EAF - makes SkyCoord object using the Laramie RA array
laramie_galactic = laramie_icrs.galactic

lar_l = laramie_galactic.l
lar_b = laramie_galactic.b

plt.plot(lar_l, lar_b, 'g*')
plt.xlabel('Galactic Longitude, l')
plt.ylabel('Galactic Latitude, b')
plt.title('Galactic Coordinates Over Laramie Over One Year')
plt.savefig('galactic_coordinates_laramie.png')
plt.show()
plt.close()


### Task 4

# EAF - the following 20 lines of code are creating ICRS SkyCoord objects for the sun, moon, and planets (+ pluto),
#       and then converting them to heliocentrictrueecliptic objects.
#       Latitude and Longitude from 2/4/2025, distance from 2/5/2025.

sun = SkyCoord(ra='21h14m42s', dec='-15d58m24s', distance=0.9860501*u.au, frame='icrs')
sun_ecliptic = sun.heliocentrictrueecliptic

moon = SkyCoord(ra='2h39m22s', dec='+19d8m58s', distance=372093*u.km, frame='icrs')
moon_ecliptic = moon.heliocentrictrueecliptic

mercury = SkyCoord(ra='21h4m28.4s', dec='-18d53m48s',distance=1.407*u.au , frame='icrs')
mercury_ecliptic = mercury.heliocentrictrueecliptic

venus = SkyCoord(ra='23h56m48.4s', dec='+2d29m9s', distance=0.489*u.au, frame='icrs')
venus_ecliptic = venus.heliocentrictrueecliptic

mars = SkyCoord(ra='7h24m49s', dec='+26d14m11s', distance=0.705*u.au, frame='icrs')
mars_ecliptic = mars.heliocentrictrueecliptic

jupiter = SkyCoord(ra='4h37m41.4s', dec='+21d36m24s', distance=4.611*u.au, frame='icrs')
jupiter_ecliptic = jupiter.heliocentrictrueecliptic

saturn = SkyCoord(ra='23h17m8.7s', dec='-6d41m39s', distance=10.453*u.au, frame='icrs')
saturn_ecliptic = saturn.heliocentrictrueecliptic

uranus = SkyCoord(ra='3h22m21.3s', dec='+18d16m4s', distance=19.417*u.au, frame='icrs')
uranus_ecliptic = uranus.heliocentrictrueecliptic

neptune = SkyCoord(ra='23h53m44.0s', dec='-2d3m40s', distance=30.630*u.au, frame='icrs')
neptune_ecliptic = neptune.heliocentrictrueecliptic

pluto = SkyCoord(ra='20h19m34.4s', dec='-22d58m54s', distance=36.146*u.au, frame='icrs')  # EAF - rest in peace Pluto the planet
pluto_ecliptic = pluto.heliocentrictrueecliptic


plt.figure(figsize=(10,7))
# EAF - the following 10 lines of code are plotting the latitude and longitued of each object
plt.plot(sun_ecliptic.lon, sun_ecliptic.lat, marker='*', markersize=10, color='gold', linestyle='None', label='Sun')
plt.plot(moon_ecliptic.lon, moon_ecliptic.lat, marker='o', markersize=8, color='red', linestyle='None', label='Moon')
plt.plot(mercury_ecliptic.lon, mercury_ecliptic.lat, marker='o', markersize=8, color='orangered', linestyle='None', label='Mercury')
plt.plot(venus_ecliptic.lon, venus_ecliptic.lat, marker='o', markersize=8, color='orange', linestyle='None', label='Venus')
plt.plot(mars_ecliptic.lon, mars_ecliptic.lat, marker='o', markersize=8, color='yellowgreen', linestyle='None', label='Mars')
plt.plot(jupiter_ecliptic.lon, jupiter_ecliptic.lat, marker='o', markersize=8, color='green', linestyle='None', label='Jupiter')
plt.plot(saturn_ecliptic.lon, saturn_ecliptic.lat, marker='o', markersize=8, color='turquoise', linestyle='None', label='Saturn')
plt.plot(uranus_ecliptic.lon, uranus_ecliptic.lat, marker='o', markersize=8, color='dodgerblue', linestyle='None', label='Uranus')
plt.plot(neptune_ecliptic.lon, neptune_ecliptic.lat, marker='o', markersize=8, color='mediumslateblue', linestyle='None', label='Neptune')
plt.plot(pluto_ecliptic.lon, pluto_ecliptic.lat, marker='o', markersize=8, color='orchid', linestyle='None', label='Pluto')


plt.legend(loc='upper right')
plt.xticks(fontsize=10)  # EAF - the following 2 lines change the x and y axes sizes
plt.yticks(fontsize=10)
plt.xlabel('Right Ascension', fontsize=12)
plt.ylabel('Declination', fontsize=12)
plt.title('Planets, Sun, and Moon in Ecliptic Coordinates', fontsize=15)

plt.savefig('solarsystem_ecliptic.png')
plt.show()

