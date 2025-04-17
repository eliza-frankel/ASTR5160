#!/usr/local/Anaconda2024/bin/python3

#########################################################
#						              
#    In Class Work: Fitting a Line with Correlated Data
#		       		   4.17.2025		  
#	   		   Author: Eliza Frankel       
#                                                             
#########################################################

import matplotlib.pyplot as plt
from scipy.stats import chi2
import numpy as np


### TASK 1 ###
file = '/d/scratch/ASTR5160/week13/line.data'

data = np.loadtxt(file, dtype=float, unpack=True)

# EAF - according to numpy documentation, np.cov() takes the parameter m
#		m is an array containing multiple variables and observations.
#		each row of m represents a variable and each column is a single observation
#		of those variables

covariance = np.cov(data)
print("The covariance matrix:")
print(covariance)

# EAF - The resulting covariance matrix should be a 10x10 matrix. This is because 
#		in line.data, there are 10 variables (the x values), which means there will
#		be 10 variances in the leading diagonal of the covariance matrix. As a result,
#		we end with a 10x10 matrix because there are 10 diagonal variables.

# EAF - axis=1 says to calcualte the variance of the rows
variances = np.var(data, axis=1, ddof=1)
print("\nVariances:")
print(variances)
print("The variances agree with the diagonal of the covariance matrix")

# EAF - Visualizing the matrix
plt.imshow(data, cmap='magma')
plt.show()


### TASK 2 ###

# EAF - making correlation matrix.

correlation = np.corrcoef(data)
print(correlation)

# EAF - determining which column is most correlated and which is least
#		Not sure if finding the sum of the values is the best way to do this, 
#		but I can't think of another way.


s = [np.sum(correlation[i]) for i in range(len(correlation))]

anti = np.argmin(s)
most = np.argmax(s)

print('The most anti-correlated column is column', anti)
print('This column is:', correlation[anti])

print('The most correlated column is column', most)
print('This column is:', correlation[most])


# EAF - ignoring the perfectly correlated diagonal

s2 = [np.sum(correlation[i]) - correlation[i][i] for i in range(len(correlation))]

most = np.argmax(s2)

print('The most correlated column ignoring the perfect diagonal is still column', most)

# EAF - I am not sure if this is the correct way to find most correlated columns/
#		if this is what the question is asking