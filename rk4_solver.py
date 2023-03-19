import numpy as np
import pandas as pd
import matplotlib as plt
import scipy.special as sp

def tumor_cell_population_function():
    """
    Creates a differential equation to model a tumor cell population over time
    
    Parameters: TODO
    ----------
    num1 : int
        First number to add.
    num2 : int
        Second number to add.

    Returns: TODO
    -------
    int
        The sum of ``num1`` and ``num2``.
    """
    return 0

def effector_cell_population_function():
    """
    Creates a differential equation to model an effector cell population over time
    
    Parameters: TODO
    ----------
    num1 : int
        First number to add.
    num2 : int
        Second number to add.

    Returns: TODO
    -------
    int
        The sum of ``num1`` and ``num2``.
    """
    return 0

def newton_raphson(function, derivative, initial_guess, x, initial_y, step_size, error):
    """
    Solves a (set of) differential equations through the Newton-Raphson (Newton's) method. Simple method for root finding
    
    Parameters:
    ----------
    function : func
        Function definition for which we are numerically solving
    derivative : func
        Derivative of function input
    initial_guess : float
        Initial guess for root
    x : float
        Initial x value
    initial_y: float
        Initial y value
    step_size: int
        Step size to determine how far we are looking for root
    error: float
        Newton's method requires a tolerance for how close the root is to predicted
    
    Returns:
    -------
    float
        Approximation of root of derivative given starting x and y values.
    """
    guess = initial_guess - (initial_guess - function(x, initial_guess) * step_size - initial_y) / (1 + derivative(x, initial_guess) * step_size)

    if np.abs((guess - initial_guess) / initial_guess) < error:
        return guess
    else:
        newton_raphson(function, derivative, guess, x, initial_y, step_size, error)

def rk4(x_start, x_finish, init_condition, num_steps, function):
    """
    Solves a (set of) differential equations through the Runge-Kutta (RK4) method
    
    Parameters:
    ----------
    start : float
        Starting point of the numerical sim
    finish : float
        End point of the numerical sim
    num_steps : int
        Number of steps to take
    init_condition : float
        Initial condition for the IVP of the ODE
    function : func
        Function definition for which we are numerically solving

    Returns:
    -------
    int
        Approximation w to function (y) at the (num_steps + 1) value of x.
    """
    step_size = abs(x_finish - x_start) / num_steps # Defines the step size to take when solving rk4 numerically
    x = [x_start] # Store x values as we compute them, with x_start being first
    y = [init_condition] # Store y values as we compute them, with init_condition being first

    for i in range(num_steps): # Loop through the number of steps to approximate solution
        # https://mathworld.wolfram.com/Runge-KuttaMethod.html
        k1 = step_size * function(x[i], y[i])
        k2 = step_size * function(x[i] + (step_size / 2), y[i] + (k1 / 2))
        k3 = step_size * function(x[i] + (step_size / 2), y[i] + (k2 / 2))
        k4 = step_size * function(x[i], y[i] + k3)

        x.append(x_start + (i + 1) * step_size) # Store new computed x
        y.append(y[i] + (k1 + (2 * k2) + (2 * k3) + k4) / 6) # Store new computed y
    
    return (x, y) # Return final (x, y) lists

if __name__ == "__main__":
    print("hello")