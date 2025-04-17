#!/usr/local/Anaconda2024/bin/python3

################################################
#						              
#    In Class Work: Fitting a Line
#		          4.15.2025		  
#	      Author: Eliza Frankel       
#                                                             
################################################

import matplotlib.pyplot as plt
from scipy.stats import chi2
import numpy as np
from itertools import product


def E(x, m_b):
	""" Determines the E/y value for a line at a specified xvalue and slope and y intercept

	Parameters:
	-----------
	x (float or int) - whatever the x value is
	m_b (tuple) - has both slope and y-intercept value in the form (m, b)

	Returns:
	-------
	(float) - E value
	"""

	return (m_b[0] * x) + m_b[1]



### NOTES ###

# chi2.sf(chi^2 val, degrees of freedom)
# print(chi2.sf(2.2958, 2))
# print(chi2.sf (1, 1) , chi2.sf (2.3, 2) , chi2.sf (3.5, 3))


### TASK 1 ###
file = '/d/scratch/ASTR5160/week13/line.data'

# EAF - unpack=True transposes the data so data is arrays of columns instead
#		of rows!
data = np.loadtxt(file, dtype=float, unpack=True)
data_range = range(len(data))

# EAF - determining the means for all y in each x range
means = [np.mean(data[i]) for i in data_range]

# EAF - determining variance for all y in each x range
variances = [np.var(data[i], ddof=1) for i in data_range]



### TASK 2 ###
plt.figure(figsize=(16, 9))
plt.grid(True, zorder=1)

# EAF - Making the x values - same as the header in file but np.loadtxt doesn't include
x = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])

# EAF - determining possible slope lines using multiple ms and bs
m_range = np.arange((22/9), (32/9), 0.5)
b_range = np.arange(3, 9.5, 1)

# EAF - itertools.product() takes 2 lists and finds all possible combinations of items
options = list(product(m_range, b_range))

E_all = [options[i][0] * x + options[i][1] for i in range(len(options))]

[plt.plot(x, E_all[i], zorder=1) for i in range(len(E_all))]
# EAF - plotting the data with corresponding x values!
for i in data_range:
	[plt.scatter(x[i], data[i][j], c='saddlebrown', s=50, zorder=3) for j in range(len(data[i]))]

plt.title('line.data Data & Possible Slopes', fontsize=20)
plt.xlabel('x-bin', fontsize=18)
plt.xticks(fontsize=15)

plt.ylabel('y measurement', fontsize=18)
plt.yticks(fontsize=15)

plt.savefig('data.png', dpi=300)



### TASKS 3 AND 4 ###

# EAF - notes on chi squared:
# χ^2 = Σ_i (O_i - E_i)^2/σ_i^2
# where E = y = mx+b and O = observed y values  
# instead of using all y values for O, I am using the mean values (and same for variance)


plt.figure(figsize=(16, 9))
plt.grid(True, zorder=1)

chi_squared = []
for num in range(len(options)):
	temp = []
	for i, j, k in zip(E_all[num], means, variances):
		val = (i - j) **2 / k**2
		temp.append(val)
	chi_squared.append(np.sum(temp))


chi_squared = np.array(chi_squared)

# EAF - finding the minimum value
minimum = np.argmin(chi_squared)

print("The best fit model parameters are m = {} and b = {}"\
	.format(round(options[minimum][0], 3), options[minimum][1]))

# EAF - The following figure looks weird with the chi^2 plot (it isn't very pretty).
#		I tried with the same values as Caleb and got the same resulting plot, so I think it
#		just has something to do with my m and b ranges.

m = np.linspace(22/9, 32/9, len(chi_squared))
[plt.scatter(m[i], chi_squared[i], zorder=2) for i in range(len(chi_squared))]

plt.title("m and b vs. $\chi^{2}$")
plt.xticks(ticks=m, labels=['({},{})'.format(round(i[0], 3), i[1]) for i in options], rotation=300)


plt.savefig('chisquare.png', dpi=300)



### TASK 5 ###

""" For each pair of parameters in your grid of m and b
determine the 68% (α = 0.32 as I defined it) and 95%
(α = 0.05) confidence limits for your parameters from ∆χ2
• remember that you’re fitting ∆χ2 for 2 parameters
"""

# ∆χ2 = χ2 - χ2min and for a 2 parameter fit, ∆χ2 = 2.3
# print(chi2.sf(2.3, 2))

delta_chi2 = [chi_squared[i] - chi_squared[minimum] for i in range(len(chi_squared))]
print(delta_chi2)


# EAF - I am not quite sure how to complete this section. I don't think I 
#		understand how exactly to use chi2.sf() or use this ∆χ2 to get the
#		alpha = 0.32 or alpha = 0.05.


