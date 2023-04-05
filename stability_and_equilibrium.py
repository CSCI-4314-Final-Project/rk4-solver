class Stability_Equilibirum_Analysis():
    def __init__(self, s, r, b, m, p, h, g, a, gamma, mu, K_e, K_t, V_m):
        """
        Initalize equilibrium analaysis with necessary parameter values 

        Parameters:
        ----------
        s: float
            growth rate of normal/effector cells
        r: float
            rate of tumor growth
        b: float
            capacity of tumor cell
        m: float
            degree of inactivation of effector
        p: float
            degree of recruitment of maximum immune-effector cells
        h: float
            steepness coefficient of the recruitment curve of effector immune cells
        g: float
            half saturation for cancer cleanup
        a: float
            parameter of cancer cleanup
        gamma: float
            rate of decrease in concentration of chemotherapy drug
        mu: float
            rate of natural demise of effector cells
        K_e: float
            rate of chemotherapy drugs causing demise of effector cells
        K_t: float
            rate of decrease in population of tumor cells due to influence of interaction between tumor cells and chemotherapy drug
        V_m: float
            the amount of concentration of chemotherapy drug that increases due to occurence of a drug outside the body
        """
        
        self.s = s
        self.r = r
        self.b = b
        self.m = m
        self.p = p
        self.h = h
        self.g = g
        self.a = a
        self.gamma = gamma
        self.mu = mu 
        self.K_e = K_e
        self.K_t = K_t
        self.V_m = V_m

    def compute_equilibrium_point(self):
        """
        Compute the equilibrium point for the 3-ODE system. To do this we will calculate 3 different points: when E=T=0, E=0, and T=0

        Parameters:
        ----------
        s: float
            growth rate of normal/effector cells
        r: float
            rate of tumor growth
        b: float
            capacity of tumor cell
        gamma: float
            rate of decrease in concentration of chemotherapy drug
        mu: float
            rate of natural demise of effector cells
        K_e: float
            rate of chemotherapy drugs causing demise of effector cells
        K_t: float
            rate of decrease in population of tumor cells due to influence of interaction between tumor cells and chemotherapy drug
        V_m: float
            the amount of concentration of chemotherapy drug that increases due to occurence of a drug outside the body

        Returns: 
        -------
        list(list(float, float, float), list(float, float, float), list(float, float, float))
            The equilibirum points for the systems when [E=0, T=0, E,T=0]
        """
        
        E_eq = [0.0, (self.r - (self.K_t*(self.V_m/self.gamma)))/(self.r * self.b), self.V_m/self.gamma]
        T_eq = [(self.s * self.gamma) / ((self.mu * self.gamma)+(self.K_e * self.V_m)), 0.0,  self.V_m / self.gamma]
        ET_eq = [0.0, 0.0, self.V_m/self.gamma]

        return [ET_eq, E_eq, T_eq]
    
    def compute_stability_analysis(self):
        """
        Compute the stability analaysis of our equilibrium point for the 3 equilibrium solutions of our system. Local asymptotic analysis for P_0 = equilbrium solution when E,T = 0, 
        P_1 = equilibrium solution when E= 0, P_2 = equilibrium solution when T=0

        Parameters:
        ----------
        s: float
            growth rate of normal/effector cells
        r: float
            rate of tumor growth
        b: float
            capacity of tumor cell
        m: float
            degree of inactivation of effector
        p: float
            degree of recruitment of maximum immune-effector cells
        h: float
            steepness coefficient of the recruitment curve of effector immune cells
        g: float
            half saturation for cancer cleanup
        a: float
            parameter of cancer cleanup
        gamma: float
            rate of decrease in concentration of chemotherapy drug
        mu: float
            rate of natural demise of effector cells
        K_e: float
            rate of chemotherapy drugs causing demise of effector cells
        K_t: float
            rate of decrease in population of tumor cells due to influence of interaction between tumor cells and chemotherapy drug
        V_m: float
            the amount of concentration of chemotherapy drug that increases due to occurence of a drug outside the body

        Returns: 
        -------
        list(tuple(list(float, float, float), int), tuple(list(float, float, float), int), tuple(list(float, float, float), int))
            This is a list of tuples where each index in the list corresponds to a tuple of (equilibrium points, stability analysis) for when E and T=0, E= 0, T=0. The stability 
            analysis will result in either a 0 or a 1, where 0=unstable and 1=local asympotically stable
        """


        p_0, p_1, p_2 = self.compute_equilibrium_point()

        # flags for each equilibrium point
        stableE, stableT, stableET = "asympotically unstable", "asympotically unstable", "asympotically unstable"

        # case one: stability analysis when E,T = 0
        if(self.V_m > ((self.r * self.gamma)/self.K_t)):
            flag_0 = 1
            stableET = "asympotically stable"
        
        # case two: stability analysis when E = 0
        T_hat = p_1[1]
        comparison1 = (((self.m*T_hat) + self.mu + self.K_e*(self.V_m/self.gamma)) * ((self.h + T_hat)/T_hat))
        comparison2 = (self.r * (1.0 - (2.0 * self.b * T_hat))) * (self.gamma/self.K_t)
        if(self.p > comparison1 and self.V_m < comparison2):
            flag_1 = 1
            stableE = "asympotically stable"
        
        # case three: stability analysis when T = 0
        E_hat = p_2[0]
        comparison = (self.r - (self.K_t*(self.V_m/self.gamma))) * (self.g/E_hat)
        if(self.a > comparison):
            flag_2 = 1
            stableT = "asympotically stable"
        
        with open("StabilityLog.txt", "w") as f:
            f.write('Equilibrium solution when E,T=0 %s' % p_0)
            f.write('\n')
            f.write('Stability analysis: %s' % stableET)
            f.write('\n')
            f.write('\n')

            f.write('Equilibrium solution when E=0 %s' % p_1)
            f.write('\n')
            f.write('Stability analysis: %s' % stableE)
            f.write('\n')
            f.write('\n')

            f.write('Equilibrium solution when T=0 %s' % p_2)
            f.write('\n')
            f.write('Stability analysis: %s' % stableT)
            f.write('\n')
            f.write('\n')

            f.close()
        
        return [(p_0, stableET), (p_1, stableE), (p_2, stableT)]




####### To test 
# st = Stability_Equilibirum_Analysis(1.2e4, 4.31 * 10**(-3), 10**(-9), 2e-11, 0.015, 2.02e1, 10**5, 3.41 * 10**(-10), 0.9, 4.12e-2, 0.6, 1, 0.5)
# st.compute_stability_analysis()