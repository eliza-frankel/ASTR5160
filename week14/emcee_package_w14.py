#!/usr/local/Anaconda2024/bin/python3

######################################################
#						              
#    In Class Work: The emcee package
#		         4.24.2025		  
#	     	Author: Eliza Frankel       
#                                                             
######################################################

import matplotlib.pyplot as plt
import numpy as np
import emcee
import likelihood_functions_w14 as lf
from scipy.optimize import minimize
import corner


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


### TASK 2 ###

"""
2.Follow the tutorial on fitting a linear model to data at
https://emcee.readthedocs.io/en/stable/tutorials/line/ and
adapt it to write your own code that uses emcee to fit a
linear model to the data in line.data
	• Note that, in the final part of the line-fitting tutorial,
	display(Math(txt)) can be replaced by print(txt)
	or print(labels[i], mcmc[1], q[0], q[1])
"""

m = 2.944
b = 5.0

x = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
y = m*x + b

slopes, inters, post, like, acceptance_rate = lf.mcmc(m, b, iterations=5000)

# EAF - determining maximum likelihood

max_index = np.argmax(like)
max_m, max_b, max_like = slopes[max_index], inters[max_index], like[max_index]



# EAF - plotting the data

plt.figure(figsize=(16, 9))
plt.grid(True, zorder=1)


for i in data_range:
	[plt.scatter(x[i], data[i][j], c='saddlebrown', s=50, zorder=2) for j in range(len(data[i]))]

# EAF - Plotting the maximum liklihood estimates
max_likelihood_eqn = max_m * x + max_b
plt.plot(x, max_likelihood_eqn, color = 'blueviolet', zorder=3, label="Maximum Likelihood Estimate")


plt.title('line.data Data', fontsize=20)
plt.xlabel('x-bin', fontsize=18)
plt.xticks(fontsize=15)

plt.ylabel('y measurement', fontsize=18)
plt.yticks(fontsize=15)
plt.legend()

# plt.savefig('data.png', dpi=300)
plt.show()
plt.close()

# EAF - following https://emcee.readthedocs.io/en/stable/tutorials/line/
# EAF - initializing the walkers

# EAF - The functions that I wrote were giving me errors regardless of 

# sampler = emcee.EnsembleSampler(nwalkers, ndim, lf.posterior_probability, args=(slopes, inters, x, y))
# sampler.run_mcmc(pos, 5000, progress=True)

# sampler = emcee.EnsembleSampler(nwalkers, ndim, lf.posterior_probability(np.array(slopes), np.array(inters), x, y))
# sampler.run_mcmc(pos, 5000, progress=True)

# EAF - this is the error I am getting:
#		posterior_probability() takes 4 positional arguments but 5 were given
#		even though the input above lf.posterior_probability, args=(slopes, inters, x, y)
#		or about inputting an array. 

# EAF - I haven't figured out how to mimick the tutorial using my functions, so I have copied
#		the functions from https://emcee.readthedocs.io/en/stable/tutorials/line/

soln = [max_m, max_b]

def log_likelihood(theta, x, y):
    m, b = theta
    model = m * x + b
    # sigma2 = yerr**2 + model**2 * np.exp(2 * log_f)
    sum_term = np.sum(np.array([((means[i] - (m*x[i] + b))**2)/variances[i]\
        + np.log(2 * np.pi * variances[i]) for i in range(len(x))]))
    ln_L = (-1/2) * sum_term
    return ln_L

def log_prior(theta):
    m, b = theta
    if (0 < b < 8) and (0 < m < 5):
        return 0.0
    return -np.inf

def log_probability(theta, x, y):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y)


pos = soln + 1e-4 * np.random.randn(32,2)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y))
sampler.run_mcmc(pos, 5000, progress=True);


fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")

plt.show()
plt.close()


tau = sampler.get_autocorr_time()
# print(tau)

# EAF - creating a flat list of samples
flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
# print(flat_samples.shape)

# EAF - creating a corner plot of posterior probability distributions
fig = corner.corner(flat_samples, labels=labels, truths=[m, b])
plt.show()
plt.close()

# EAF - plotting some of the results over real data
inds = np.random.randint(len(flat_samples), size=100)
for i in data_range:
	[plt.scatter(x[i], data[i][j], c='saddlebrown', s=50, zorder=2) for j in range(len(data[i]))]

for ind in inds:
    sample = flat_samples[ind]
    plt.plot(x, np.dot(np.vander(x, 2), sample[:2]), "C1", alpha=0.1)
plt.plot(x, m * x + b, "k", label="truth")
plt.legend(fontsize=14)
plt.xlim(0, 10)
plt.xlabel("x")
plt.ylabel("y")
plt.show()
plt.close()

for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    txt = "{} = {} -{}/+{}"
    txt = txt.format(labels[i], round(mcmc[1],3), round(q[0],3), round(q[1],3))
    print(txt)