import math
import numpy as np
import matplotlib.pyplot as plt

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

def tumor_effector_chemo_forward_integrate(initial_conditions, params, t):
    """
    Forward integrates the tumor, effector, chemotherapy model
    
    Parameters: 
    ----------
    initial_conditions : (int, int, int)
        Tuple of initial conditions for the IVP
    params : floats
        parameters of the model
    t : [ints]
        an array of timepoints, ASSUMED TO BE EQUALLY SPACED. 
    Returns: 
    -------
    np.array, np.array, np.array
        Returns a numpy array of T, E, and M arrays.
    """
    T0, E0, M0 = initial_conditions
    T, E, M = [T0], [E0], [M0]
    p, r, b, a, g, s, m, mu, gamma, h, K_t, K_e, V_m = params
    dt = t[1] - t[0]

    for _ in t[1:]:
        Tt = T[-1] + (r*T[-1]*(1-(b*T[-1]))-a*((E[-1]*T[-1])/(T[-1]+g))-(K_t*M[-1]*T[-1]))*dt
        Et = E[-1] + (s-(mu*E[-1]) + p*((E[-1]*T[-1])/(h+T[-1])) - (m*E[-1]*T[-1]) - (K_e*M[-1]*E[-1]))*dt
        Mt = M[-1] + (-gamma*M[-1]+V_m)*dt
        T.append(Tt)
        E.append(Et)
        M.append(Mt)
    return T, E, M

if __name__ == "__main__":
    T0 = 40000
    E0 = 30000
    M0 = 0
    initial_conditions = (T0, E0, M0)

    p = 0.015
    r = 0.00431
    b = 10
    a = 3.41
    g = 10e5
    s = 1.2e4
    m = 2
    mu = 4.12
    gamma = 0.9
    h = 2.02
    K_t = 10
    K_e = 5
    V_m = 2
    params = p, r, b, a, g, s, m, mu, gamma, h, K_t, K_e, V_m

    # Timesteps in days
    t_max = 5
    dt = 1
    t = np.linspace(0, t_max, int(t_max/dt) + 1)

    results = tumor_effector_chemo_forward_integrate(initial_conditions, params, t)

    # Plot Results
    colors = {
    'T':'#22223B',
    'E':'#F13030',
    'M':'#FFA737'}

    # Set up the axes
    T,E,M = results
    plt.plot(t,T,label='T',color=colors['T'])
    plt.plot(t,E,label='E',color=colors['E'],linewidth=2)
    plt.plot(t,M,label='M',color=colors['M'])

    # Make the plot attractive
    plt.legend(loc='best')
    plt.xlabel('time (days)')
    plt.ylabel('Cell Population')


# Constants for tumor cell proliferation IVP => T(0) = 40000 cells
T_0 = 0
Y_0 = 40000

T_FINAL = 50 # This is the final time/timestep the simulation should end at. This was the final time in the numerical simulation graphs from Lestari et al.
DT = 1 # Timestep size used in Lestari et al. simulations
DT_EXACT = 0.1 # Higher resolution timestep...May not be needed

def rk4_step(t_i, y_i, dt, f):
    k1 = f(t_i, y_i)
    k2 = f(t_i + dt/2, y_i + k1*dt/2)
    k3 = f(t_i + dt/2, y_i + k2*dt/2)
    k4 = f(t_i + dt, y_i + k3*dt)
    return y_i + dt/6*(k1 + 2*k2 + 2*k3 + k4)
