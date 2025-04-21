#!/usr/local/Anaconda2024/bin/python3

######################################################
#						              
#    In Class Work: Likelihood Functions and MCMC
#		        	  4.22.2025		  
#	     		 Author: Eliza Frankel       
#                                                             
######################################################

import numpy as np
from numpy import random

### TASK 1 ###

"""
In /d/scratch/ASTR5160/week13/ is a file of (x, y)
data called line.data . The columns correspond to x bins
of 0 < x < 1 (bin x0), 1 < x < 2 (x1) ... up to 9 < x < 10
(x 9). Each of the 20 rows is a y measurement in that x bin
	• Read the file and find the variance ( σi; remember to pass
	ddof=1 ) and mean (y i) of the y data in each bin (i) of x
"""

file = '/d/scratch/ASTR5160/week13/line.data'

data = np.loadtxt(file, dtype=float, unpack=True)
data_range = range(len(data))

# EAF - determining the means & variances for all y in each x range
means = [np.mean(data[i]) for i in data_range]
variances = [np.var(data[i], ddof=1) for i in data_range]


def likelihood(m, b, x, y):
	""" Calculates the ln(Likelihood)

	Parameters:
	-----------
	m (int or float) - Slope of the linear fit
	b (int or float) - y-intercept of the linear fit
	x (np.array) - array of all bin values
	y (np.array) - y = mx + b

	Returns:
	--------
	ln_L (float) - the log of the likelihood of these values
	"""

	# EAF - Equation 1 - makes a np array and gets the sum of each term
	sum_term = np.sum(np.array([((means[i] - (m*x[i] + b))**2)/variances[i]\
		+ np.log(2 * np.pi * variances[i]) for i in range(len(x))]))
	ln_L = (-1/2) * sum_term

	return ln_L



def posterior_probability(m, b, x, y):
	""" Calculates the (ln) posterior probability for a linear fit of data

	Parameters:
	-----------
	m (int or float) - Slope of the linear fit
	b (int or float) - y-intercept of the linear fit
	x (np.array) - array of all bin values
	y (np.array) - y = mx + b

	Returns:
	--------
	(float) - posterior probability

	Note:
	-----
	Posterior probability ∝ Likelihood (L) x Prior
	Because this function finds ln(Posterior probability), this equation becomes:
		ln(PP) ∝ ln(L x Prior) = ln(L) + ln(Prior)
	"""

	ln_L = likelihood(m, b, x, y)

	if (0 < b < 8) and (0 < m < 5):
		ln_Prior = 0
	else:
		# EAF - if (b > 8) or (b < 0):
		# EAF - Prior = 0
		ln_Prior = - np.inf

	return ln_L + ln_Prior



def mcmc(m, b, iterations=5000):
	""" Uses Metropolis-Hastings algorith to create MCMC chain

	Parameters:
	-----------
	m (int or float) - Slope of the linear fit
	b (int or float) - y-intercept of the linear fit
	iterations (int) - number of times to iterate through the M-H algorithm.
		If not specified, it will default to 5000 iterations. This value was
		chosen because of the emcee documentation for fitting a line 
		(https://emcee.readthedocs.io/en/stable/tutorials/line/)

	Returns:
	-------
	slopes (list) - list of slopes. Found using a Gaussian proposal function
	intercepts (list) - list of intercepts. Found using a Gaussian proposal function
	posterior (list) - list of the posterior probability
	likelihoods (list) - list of the likelihoods for each parameter
	acceptance_rate (float) - The acceptance rate of R (should be around 30%)
	"""

	x = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
	y = m * x + b

	slopes = [m]
	intercepts = [b]
	posterior = [posterior_probability(m, b, x, y)]
	likelihoods = [likelihood(m, b, x, y)]

	accepted = 0

	for i in range(iterations):
		# EAF - making the shifted slope and intercept variables
		# EAF - indexing using -1 to account for not accepting new parameters.
		#		Takes the last item appended to the lists.

		G_prop_slopes = random.normal(slopes[-1], 0.35)
		G_prop_intercepts = random.normal(intercepts[-1], 0.35)

		# EAF - defining variables for current m and b (for simplicity sake)
		m_old = slopes[-1]
		b_old = intercepts[-1]

		# lnR = lnPnew - lnPold 
		ln_R = posterior_probability(G_prop_slopes, G_prop_intercepts, x, y)\
			- posterior_probability(m_old, b_old, x, y)
		R = np.exp(ln_R)

		# EAF - Determining if R is 
		if R > 1:
			# EAF - of R > 1, take the new pro
			slopes.append(G_prop_slopes)
			intercepts.append(G_prop_intercepts)
			posterior.append(posterior_probability(G_prop_slopes, G_prop_intercepts, x, y))
			likelihoods.append(likelihood(G_prop_slopes, G_prop_intercepts, x, y))
			accepted += 1
		else:
			# EAF - if R < 1, generate a random number between 0 and 1
			rand_number = random.random()

			if R < rand_number:
				# EAF - if R is less than rand_number, reject these new parameters
				# EAF - 'continue' skips the rest of the loop and moves to next iteration
				continue
			elif R > rand_number:
				# EAF - if R is greater than rand_number accept these new parameters
				slopes.append(G_prop_slopes)
				intercepts.append(G_prop_intercepts)
				posterior.append(posterior_probability(G_prop_slopes, G_prop_intercepts, x, y))
				likelihoods.append(likelihood(G_prop_slopes, G_prop_intercepts, x, y))
				accepted += 1

	acceptance_rate = (accepted / iterations) * 100

	return slopes, intercepts, posterior, likelihoods, round(acceptance_rate, 3)



if __name__ == "__main__":


	### TASK 2 ###

	"""
	The data have been drawn from a line of the form y = mx
	+ b and scattered according to a Gaussian.
		• Write a function to calculate (ln) posterior probability
		(Eqn. 1) for a linear fit to the data when passed m and b

		• Use “flat” priors that reflect the extent of the m, b space
		you need to sample. Say, 0 < b < 8; beyond this range
		Prior = 0, so ln(Prior) = a very -ve number, e.g. -np.inf 
	"""

	# EAF - From fitting_a_line_w13.py, I found that the best fit parameters are m = 2.944 and b = 5.0

	m = 2.944
	b = 5.0

	x = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
	y = m * x + b


	print(posterior_probability(m, b, x, y))


	### TASK 3 ###

	"""
	3.Using the Metropolis-Hastings algorithm walk through
	the parameter space and create an MCMC chain
		• Start at values of m and b that seem reasonable

		• Use a Gaussian proposal function with a step size (the
		standard deviation of the Gaussian) of 0.1

		• Note that “accepting new parameters with probability
		R” can be done by generating a random number
		between 0 and 1 and checking that it’s less than R
	"""

	# EAF - Calling mcmc function above (I made it a function to be callable for 
	#		Thursday's tasts)

	slopes, iters, post, like, acceptance_rate = mcmc(m, b, iterations=5000)


	### TASK 4 ###

	# EAF - I have updated the step size to 0.35 to get an acceptance rate of ~ 30%

	print("Acceptance rate of {} %".format(acceptance_rate))


	### TASK 5 ###
	# print(np.argmax(like), like[np.argmax(like)])
