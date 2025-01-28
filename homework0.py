# Homework 0 - ASTR 5160
# Eliza Frankel

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot():
    """ Compiles all defined functions into one; Plots a line of 10 random data points with two slopes 
        - one inputted by the user, the other determined to be the optimal fit using scipy.curve_fit

    User input:
    m - class: 'int' or 'float' - converted to float
    b - class: 'int' or 'float' - converted to float

    Notes:

    """
    m = float(input("Please input a slope, m = "))
    b = float(input("Please input a y-intercept, b = "))
    plot_data(m, b)

def create_data():
    """ Creates all components of a line to eventually graph

    Returns:
    x - class: 'np.array' - 10 random floats
    y - class: 'np.array' - the values of x scattered by a random value
    y_error - class: 'np.array' - array of error on y

    """

    x = np.random.uniform(0, 10.0, 10)  # EAF randomly chooses 10 numbers between 0 and 10.0
    y = scatter(x)   
    y_error = np.ones(10)*(0.5)  # EAF - Creates array of error; in this case, y_error is always 0.5  

    return x, y, y_error

def eqn(m, x, b):
    """ Returns values of y found through the equation y = mx + b
    
    Parameters:
	m - class: 'int' or 'float'
	x - class: 'np.array'
	b - class: 'int' or 'float'

    Returns:
	- class: 'np.array'
    """

    y = m * x + b
    return y

def scatter(y):
    """ Scatters values of y using np.random.normal

    Parameters:
    y - class: 'np.array'

    Returns:
    - class: 'np.array'
    
    Notes:
    Once np.random.normal has been found, add these values to the parameter y
    """

    offset = np.random.normal(0, 0.5, 10)  # EAF - the gaussian is centered at 0, providing an offset 
    y_adjusted = y + offset  # EAF - randomly changes the value of y by the determined offset value
    return y_adjusted

    
def line_fit(x, y):
    """ Uses curve_fit to determine optimal slope and y intercept values

    Parameters:
    x - class: 'np.array' or list
    y - class: 'np.array' or 'list'

    Returns:
    - class: 'float'
    - class: 'float'

    Notes: 
    Ignoring the covariance array given through scipy.curve_fit
    The results of curve_fit are the modified m & b values.
    This function is slightly redudant, but nice to have. 
    """

    optimal, cov = curve_fit(eqn, x, y)  # EAF - curve_fit returns 2 arrays; the optimal values for parameters, and the covariance of the optimal values. Ignore cov in this case. 
    return optimal[0], optimal[1]

def plot_data(m, b):
    """ Plots the data and best fit lines, both inputted values and those determined using curve_fit
    
    Parameters:
    m - class: 'int' or 'float'
    b - class: 'int' or 'float'

    Returns:
    Nothing

    Notes:
    Creates a PNG file of the plot
    """

    x, y, y_error = create_data()  # EAF - this creates the data and error.
    plt.plot(x, y, marker='o', color='k', linestyle='none', markersize=5, label='Randomly Generated Data')
    plt.errorbar(x, y, y_error, capsize=4, ecolor='k', linestyle='none')  # EAF - plots black errorbars.

    # EAF - The following plots the slope given by the user
    plt.plot(x, eqn(m, x, b), color='tab:olive', linewidth=2, label='User Slope')

    #EAF - The following plots the optimized slope, found using scipy.curve_fit
    m2, b2 = line_fit(x, y)
    plt.plot(x, eqn(m2, x, b2), color='b', linestyle='dotted', linewidth=3, label='Optimized Slope')
 
    plt.legend()  # EAF - displays the legend
    plt.savefig('homework0_plot.png')
    plt.show()

if __name__ == '__main__':
    plot()
