#!/usr/local/Anaconda2024/bin/python3

#######################################
#						              
#    In Class Work: SQL queries
#		       3.10.2025			  
#	      Author: Eliza Frankel       
#                                                             
#######################################

import numpy as np
import matplotlib.pyplot as plt


###  TASK 2  ###

#  EAF - file comes from SDSS SQL Query I made in the following format:
"""
select
    p.ra, p.dec, p.g
from
    photoObj p, dbo.fGetNearbyObjEq(300,-1,2) n
where
    p.objID = n.objID
"""

file = 'result.csv'

#  EAF - reading in csv, specifying to split at commas, and ignoring first row.
data = np.loadtxt(file, skiprows=1, delimiter=',', dtype=float)

#  EAF - separating RA, Dec, and g magnitudes
ra = data[:, 0]
dec = data[:, 1]
g = data[:, 2]

#  EAF - Plotting data!

plt.figure(figsize=(16,9))
plt.scatter(ra, dec, marker='o', c='chocolate', s=50)

# EAF - changing the x ticks to not be in scientific notation.
x = [299.97, 299.98, 299.99, 300.0, 300.01, 300.02, 300.03]
plt.xticks(x, x, fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Right Ascension', fontsize=18)
plt.ylabel('Declination', fontsize=18)
plt.title('SDSS Objects: RA vs Dec', fontsize=20)

plt.savefig('ra_vs_dec.png', dpi=300)
plt.show()
plt.close()


###  TASK 3  ###

# s = [20*2**n for n in range(len(x))]

#  EAF - getting the indices sorted by g-band magnitude from lowest to highest
g_sort_indices = np.argsort(g)
length = len(g_sort_indices)

ra_sort = ra[g_sort_indices]
dec_sort = dec[g_sort_indices]
g_sort = g[g_sort_indices]

# creating sizes that decrease with increasing g
s = [20*0.4*(length - n) for n in range(length)]

plt.figure(figsize=(16,9))
plt.scatter(ra_sort, dec_sort, s=s, marker='o', c=g_sort, cmap='magma', alpha=0.8)

x = [299.97, 299.98, 299.99, 300.0, 300.01, 300.02, 300.03]
plt.xticks(x, x, fontsize=15)

plt.colorbar(label='g-band magnitude')
plt.xlabel('Right Ascension', fontsize=18)
plt.ylabel('Declination', fontsize=18)
plt.title('SDSS Objects: RA vs Dec\nSize and Color Varying with g-band Magnitude', fontsize=20)

plt.savefig('g_band_mag_ordering.png', dpi=300)
plt.show()
plt.close()


###  TASK 4  ###

# I think that my plot looks reasonably like the SDSS Navigate Tool image. 
# Obviously it isn't a perfect match but I think it's relatively similar.