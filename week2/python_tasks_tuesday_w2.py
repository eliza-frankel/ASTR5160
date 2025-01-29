from astropy.table import QTable
from astropy.table import Table  # EAF - for python tasks. This is for reading in fits files
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

##### EAF - Following the astropy documentation to create numpy rec arrays as tables

a = np.array([1, 4, 5], dtype=np.int32)  # EAF - numpy array with specified type int32
b = [2.0, 5.0, 8.5]  # EAF - list of float values
c = ['x', 'y', 'z']  # EAF - list of strings
d = [10, 20, 30] * u.m / u.s  # EAF - Quantity array (a Quantity is a number with associated units)

t = QTable([a, b, c, d], names=('a', 'b', 'c', 'd'), meta={'name': 'first table'})

# t.show_in_browser(jsviewer=True)  # EAF - displays interactive table in browser


##### EAF - PYTHON TASKS DAY 3

struc = Table.read("/d/scratch/ASTR5160/week2/struc.fits")  # reads in the struc.fits file from the network
RA = struc["RA"]   # EAF - creates a variable for the RA column
DEC = struc["DEC"]  # EAF - creates a variable for the Dec column
extinction = struc["EXTINCTION"][:,0]  # EAF - creates a variable of the extinction column

# EAF - Boolean indexing struc.fits where extinction is greater than 0.22
big_extinction = extinction > 0.22  # EAF - creates a mask. True if extinction > 0.22, False otherwise
extinction_RA = RA[big_extinction]  # EAF - puts the values of RA where extinction > 0.22 in an array
extinction_DEC = DEC[big_extinction]  # EAF - puts the values of declination where extinction > 0.22 in an array

plt.plot(RA, DEC, 'ko', label="All objects")  # EAF - plots all values of RA and Dec
plt.plot(extinction_RA, extinction_DEC, 'r*', label="Extinction > 0.22")  # EAF - only plots RA and Dec for objs with extinction > 0.22

plt.xlabel(r'Right Ascension, $\alpha$')  # EAF - plots the declination  # r allows me to use LaTeX formatting
plt.ylabel(r'Declination, $\delta$')  # EAF - plots RA
plt.legend()
plt.tight_layout()  # EAF - this line removes excess whitespace in the final product
plt.savefig("ra_vs_dec_extinction.png")
plt.show()



##### EAF - Generating recarray & writing to a FITS file

rand_1 = np.random.uniform(0,1000,100)  # EAF - randomly generates 100 random numbers between 0 and 1000
rand_2 = np.random.uniform(0,1000,100)
rand_3 = np.random.uniform(0,1000,100)

a = np.hstack((rand_1, rand_2, rand_3))  # EAF - stacks the 3 randomly generated arrays into one
randomnum = np.reshape(a, (100, 3))  # EAF - reshapes it into an array with 3 columns and 100 rows. I will be honest, I'm not sure why we are using this because it needs an array as an argument.

rec_array = QTable([RA, DEC, randomnum], names=('RA', 'DEC', 'randomnum'))  # EAF - creates a rec array using provided information.

FITS = rec_array.write('recarray.fits', format='fits', overwrite=True)   # EAF - writes the recarray to a fits file


# EAF - I'm not sure if I did the final task (generating a recarray) correctly