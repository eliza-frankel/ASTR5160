import matplotlib.pyplot as plt
import numpy as np

def plot_data():
    data = np.loadtxt('w1_PythonTasks.txt')

    plt.plot(data[:,0], data[:,1])
    plt.savefig('line.png')
    plt.show()
    plt.close()

    plt.plot(data[:,0], data[:,1], color='yellow', marker='o', linestyle='none')
    plt.savefig('yellow_dots.png')
    plt.show()
    plt.close()

# if __name__ == "__main__":
#    plot_data()

def quad(x):
    """  Returns the answer to a quadratic equation

    Parameters:
        x : type: 'int' or 'float' or np.array
    
    Returns:
        y : the answer to the given quadratic equation
    """

    y = x**2 + 3*x + 8
    return y

def plot(x):
    """
    Produces a plot of the quadratic from the function quad(x)
    """
    import numpy as np
    import matplotlib.pyplot as plt

    plt.plot(x, quad(x))
    plt.savefig('plotting_quad.png')
    plt.show()

if __name__ == "__main__":
    x = np.arange(10)-4.5
    quad(x)
    plot(x)

