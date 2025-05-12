#!/usr/local/Anaconda2024/bin/python3

######################################################
#						              
#    Final Project: Fitting Models to Data
#		             5.16.2025		  
#	     	   Author: Eliza Frankel       
#                                                             
######################################################

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

import emcee
from scipy.optimize import minimize
import corner



def max_likelihood(parameters, x, y, yerr, model_type):
    """ This is modified code from the emcee fitting model to data website,
    https://emcee.readthedocs.io/en/stable/tutorials/line/
    Finds the most parameters with the maximum liklihood.

    Parameters:
    -----------
    parameters (tuple or list) - contains the parameters for the desired
        equation: m, b for linear and a2, a1, a0 for quadratic
    x (np.array) - x values
    y (np.array) - y values
    yerr (np.array) - list of y error
    model_type (str) - name of type of model desired

    Returns:
    --------
    max_parameters (np.array) - most maximized parameters

    Notes:
    ------
    The variable names for this function come from the emcee documentation
        for the simplicity of my understanding
    """
    nll = lambda *args: -likelihood(*args)

    # EAF - minimizing the function. The -likelihood_function is so that
    #		whatever is found to be most minimum is actually the maximum.
    initial = np.array(parameters) + 0.1 * np.random.randn(len(parameters))

    soln = minimize(nll, initial, args=(x, y, yerr, model_type))

    # max_parameters = soln.x

    # return max_parameters
    return soln.x



def likelihood(parameters, x, y, yerr, model_type):
    """ Calculates the ln(Likelihood) for the both equations
    *** This is my code from week14/likelihood_functions_w14.py, although
    	I have modified it to actually work with the emcee module.

    Parameters:
    -----------
    parameters (tuple or list) - contains the parameters for the desired
    	equation: m, b for linear and a2, a1, a0 for quadratic
    x (np.array) - x values
    y (np.array) - y values
    yerr (np.array) - list of y error
    model_type (str) - name of type of model desired

    Returns:
    --------
    ln_L (float) - the log of the likelihood of these values
    """
    
    # EAF - getting the proper equation for the model depending on 
    #		specified equation type

    if model_type == "linear":
        m, b = parameters
        model = m * x + b
    elif model_type == "quadratic":
        a2, a1, a0 = parameters
        model = (a2 * (x**2)) + (a1 * x) + a0

    # EAF - calculating the likelihood
    ln_L = np.sum((((y - model)**2) / (yerr**2)) + np.log(2 * np.pi * yerr**2))

    return (-1/2) * ln_L



def posterior_probability(parameters, x, y, yerr, model_type):
    """ Calculates the (ln) posterior probability
    *** This is my code from week14/likelihood_functions_w14.py, although
        I have modified it to actually work with the emcee module and to
        be able to choose which model used

    Parameters:
    -----------
    parameters (tuple or list) - contains the parameters for the desired
        equation: m, b for linear and a2, a1, a0 for quadratic
    x (np.array) - x values
    y (np.array) - y values
    yerr (np.array) - list of y error
    model_type (str) - name of type of model desired

    Returns:
    --------
    (float) - posterior probability

    Note:
    -----
    Posterior probability ∝ Likelihood (L) x Prior
    Because this function finds ln(Posterior probability), this equation becomes:
        ln(PP) ∝ ln(L x Prior) = ln(L) + ln(Prior)
    """

    ln_L = likelihood(parameters, x, y, yerr, model_type)

    # EAF - Determining the priors depending on which type of function
    #       is being used
    if model_type == "linear":
        m, b = parameters
        if -5.0 < m < 0.5 and 0 < b < 7:
            ln_prior = 0.0
        else:
            ln_prior = -np.inf

    elif model_type == "quadratic":
        a2, a1, a0 = parameters
        if 0.0 < a2 < 0.12 and -5 < a1 < 5 and 0 < a0 < 15:
            ln_prior = 0.0
        else:
            ln_prior = -np.inf

    return ln_L + ln_prior



def using_emcee(parameters, x, y, yerr, model_type, num_steps=5000):
    """This is modified code from the emcee fitting model to data website,
    https://emcee.readthedocs.io/en/stable/tutorials/line/
    This function uses the emcee package to find projections and plots them.

    Parameters:
    -----------
    parameters (tuple or list) - contains the parameters for the desired
        equation: m, b for linear and a2, a1, a0 for quadratic
    x (np.array) - x values
    y (np.array) - y values
    yerr (np.array) - list of y error
    model_type (str) - name of type of model desired
    num_steps (int) - number of steps to do in MCMC. If not specified,
        defaults to 5000 steps

    Notes:
    ------
    Most of the  variable names for this function come from the emcee 
        documentation for the simplicity of my understanding
    """

    max_parameters = max_likelihood(parameters, x, y, yerr, model_type)

    # EAF - Initializing walkers in a tiny Gaussian ball around max 
    #       likelihood result
    pos = max_parameters + 1e-4 * np.random.randn(32, len(parameters))

    nwalkers, ndim = pos.shape

    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior_probability,\
        args=(x, y, yerr, model_type))
    sampler.run_mcmc(pos, num_steps, progress=True)

    # EAF - creating a flat list of samples
    # EAF - There are 10432 samples total
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

    mcmc, negative, positive = calculate_percentile(ndim, flat_samples, model_type)

    # EAF - plotting the cornerplot
    plotting(flat_samples, mcmc, negative, positive, parameters, x, y, yerr,\
        model_type, cornerplot=True, projection=False)

    # EAF - plotting the regular scatter plot with best fit model overlaid
    plotting(flat_samples, mcmc, negative, positive, parameters, x, y, yerr,\
        model_type, cornerplot=False, projection=True)



def plotting(flat_samples, max_params, neg_err, pos_err, parameters,\
    x, y, yerr, model_type, cornerplot=False, projection=False):
    """ Creates figures. Can either make a corner plot, or a scatter plot.
    Has the option to add projections if desired, however it defaults to False.

    Parameters:
    -----------
    flat_samples - flat list of samples from MCMC
    max_parameters (np.array) - most maximized parameters
    neg_err (list) - list of negative error for most maximized params
    pos_err (list) - list of positive error for most maximized params
    parameters (tuple or list) - contains the parameters for the desired
        equation: m, b for linear and a2, a1, a0 for quadratic
    x (np.array) - x values
    y (np.array) - y values
    yerr (np.array) - list of y error
    model_type (str) - name of type of model desired
    corner (Boolean) - Determines which type of plot to make.
        Defaults to a regular scatter plot, but if True
        makes a corner plot
    """

    # EAF - setting labels, parameters, and models depending on 
    #       type of model desired.

    if model_type == "linear":
        m, b = max_params
        param_labels = ["m", "b"]
        model = m * x + b

        # EAF - creating text for the legend.
        # EAF - The following formatting for m_txt and b_txt was written with 
        #       help from Lexi Boone.
        # EAF - the 'f' is for formatting and 'r' is to allow for LaTeX syntax
        m_txt = fr"{max_params[0]:.3f}$\pm$ $^{{{pos_err[0]:.3f}}}_{{{neg_err[0]:.3f}}}$"
        b_txt = fr"{max_params[1]:.3f}$\pm$ $^{{{pos_err[1]:.3f}}}_{{{neg_err[1]:.3f}}}$"
        eqn = m_txt + ' x + ' + b_txt

    elif model_type == "quadratic":
        a2, a1, a0 = max_params
        param_labels = ["a2", "a1", "a0"]
        model = (a2 * (x**2)) + (a1 * x) + a0

        # EAF - creating text for the legend
        a2_txt = fr"{max_params[0]:.3f}$\pm$ $^{{{pos_err[0]:.3f}}}_{{{neg_err[0]:.3f}}}$"
        a1_txt = fr"{max_params[1]:.3f}$\pm$ $^{{{pos_err[1]:.3f}}}_{{{neg_err[1]:.3f}}}$"
        a0_txt = fr"{max_params[2]:.3f}$\pm$ $^{{{pos_err[2]:.3f}}}_{{{neg_err[2]:.3f}}}$"

        eqn = a2_txt + r" x$^2$ + " + a1_txt + " x + " + a0_txt

    # EAF - plotting either corner or regular x vs y plot:

    if projection == True:
        # EAF - plotting the projection of results, plotting a few sample
        #       chains

        # EAF - I am choosing the first 150 samples to plot
        inds = np.random.randint(len(flat_samples), size=150)
        for ind in inds:
            sample = flat_samples[ind]
            plt.plot(x, np.dot(np.vander(x, len(parameters)),\
                sample[:len(parameters)]), "skyblue", alpha=0.1)


    if cornerplot==True:
        # EAF - creating a cornerplot using the sampled data
        fig = corner.corner(flat_samples, labels=param_labels)
        plt.savefig(model_type + "_corner_mcmc.png", dpi=300)
        plt.show()
        plt.close()
    else:
        plt.errorbar(x, y, yerr=yerr, ms=15, fmt='.k', capsize=0)
        plt.plot(x, model, color='k', label = "Best-fitting parameters:\n" + eqn)
        plt.title("The {} best-fit model".format(model_type))
        plt.legend()
        plt.savefig(model_type + "_mcmc.png", dpi=300)
        plt.show()   
        plt.close() 



def calculate_percentile(ndim, flat_samples, model_type):
    """ This is modified code from the emcee fitting model to data website,
    https://emcee.readthedocs.io/en/stable/tutorials/line/

    Parameters:
    -----------
    ndim (int) - number of dimensions of the best fit parameters
    flat_samples (np.array) - array of a flat set of samples
    model_type (str) - name of type of model desired. For print statement.

    Returns:
    --------
    mcmc_value (list) - list of best mcmc values
    negative_error (list) - list of negative error for mcmc value
    positive_error (list) - list of positive error for mcmc value
    """
    if model_type == "linear":
        param_labels = ["m", "b"]
    elif model_type == "quadratic":
        param_labels = ["a2", "a1", "a0"]

    print("\nThe {} model!\n".format(model_type))

    mcmc_value = []
    negative_error = []
    positive_error = []

    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:,i], [16, 50, 84])
        q  = np.diff(mcmc)
        mcmc_value.append(mcmc[1])
        negative_error.append(q[0])
        positive_error.append(q[1])

        print("For {}, the best-fitting parameters are:"\
            .format(param_labels[i]))
        print("16th percentile: {0:.3f} \u00B1 ({1:.3f}, -{2:.3f})"\
            .format(mcmc[0], q[1], q[0]))
        print("50th percentile: {0:.3f} \u00B1 ({1:.3f}, -{2:.3f})"\
            .format(mcmc[1], q[1], q[0]))
        print("84th percentile: {0:.3f} \u00B1 ({1:.3f}, -{2:.3f})\n"\
            .format(mcmc[2], q[1], q[0]))

    return mcmc_value, negative_error, positive_error



if __name__ == "__main__":

    # EAF - filepath to the experimental data
    file_path = '/d/scratch/ASTR5160/final/dataxy.fits'

    # EAF - reading in the data
    data = Table.read(file_path)

    # EAF - separating the columns by data type
    x = data["x"]
    y = data["y"]
    yerr = data["yerr"]

	# EAF - These are the starting values for m and b for the linear case
    # EAF - I found these by eye.
    m = -0.8
    b = 3

    np.random.seed(123456)


    m_max, b_max = max_likelihood((m, b), x, y, yerr, 'linear')
    using_emcee((m_max, b_max), x, y, yerr, 'linear', num_steps=5000)



    # EAF - These are the starting values for a2, a1, a0 for quadratic case
    # EAF - I found these by eye.
    
    # EAF - a2 determines the tightness of the parabola and direction
    a2 = 0.04
    # EAF - a1 determines which part of parabola intersects with the y-axis 
    a1 = -1.8
    # EAF - a0 is like the y intercept parameter
    a0 = 7

    a2_max, a1_max, a0_max = max_likelihood((a2, a1, a0), x, y, yerr, 'quadratic')
    using_emcee((a2_max, a1_max, a0_max), x, y, yerr, 'quadratic')


    which_model = """
    Looking at the probability distribution for a_2, it can be seen that there 
    is a correlation between the variables. On this plot, the resulting 
    contours are very elliptical - if they were uniformly spherical (and 
    small), it would mean that they are not correlated with one another.
    Looking at the best fit equation found using the MCMC method, the value for 
    a_2 was found to be 0.060 \u00B1 0.026. This value is very 
    close to zero. If this value were to be zero, the resulting equation would 
    still likely fit the points - however it would be a linear fit rather than 
    quadratic.
    In conclusion, a quadratic model can be used to represent the data, but
    it is not necessary. A linear model is sufficient with this data.
    """

    print(which_model)




