import numpy as np
import pandas as pd
import matplotlib as plt
import scipy.special as sp


def tumor_cell_population_function(T,E,M):
    """
    Creates a differential equation to model a tumor cell population over time
    
    Parameters: 
    ----------
    r : float
        rate of tumor growth
    T : float
        number of tumor cells for time t
    b : float
        capacity of the tumor cell
    a : float
        paramter of cancer cleanup
    E : float
        number of effector cells at time t
    g : float
        half-saturation for cancer cleanup
    K_t : float
        rate of tumor cell death by chemotherapy drug
    M : float
        concentration of chemotherapy drug at time t
    Returns: 
    -------
    float
        change in tumor cells after one timestep
    """
    r = 4.31 * 10**(-3)
    b = 10**(-9)
    a = 3.41 * 10**(-10)
    g = 10**5
    K_t = 1 # not sure what this is yet
    dT = r*T*(1-(b*T))-a*((E*T)/(T+g))-(K_t*M*T)
    return dT

def effector_cell_population_function(T,E,M):
    """
    Creates a differential equation to model an effector cell population over time
    
    Parameters: 
    ----------
    s : float
        Growth rate of normal / effector cells
    mu : float
        Rate of natural demise of effector cells
    E : int
        Number of effector cells at time t
    p : float
        degree of recruitment of maximum immune-effector cells
    h : int
        steepness coefficient of the recruitment curve of effectorimmune cells 
    T : int
        Number of tumor cells at time t
    m : float
        degree of inactivation of effector cells by tumor cells
    M : int
        concentration fo chemotherapy drug
    K_e : float
        rate at which chemotherapy drugs kill effector cells
    Returns: 
    -------
    float
        The change in efffector cells after one timestep
    """
    s = 1.2 * 10**4
    mu = 4.12 * 10**(-2)
    p = 0.015
    h = 20.2
    m = 2 * 10**(-11)
    K_e = 1 # not sure what this is yet
    dE = s-(mu*E) + p*((E*T)/(h+T)) - (m*E*T) - (K_e*M*E)
    return dE

def chemotherapy_drug_concentration_function(M):
    """
    Creates a differential equation to model the concentration of chemotherapy drug over time
    
    Parameters: 
    ----------
    gamma : float
        rate of decrease in concentration of chemotherapy drug
    M : int
        concentration fo chemotherapy drug
    V_m : float
        outside addition of drug, in the litterature, this is a function, not sure how it is a function of time. 
    Returns: 
    -------
    float
        The change in chemotherapy drug concentration
    """
    gamma = 0.9
    V_m = 1 # not sure what this is yet
    dM = -gamma*M+V_m
    return dM

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