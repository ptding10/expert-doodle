import math
import numpy as np
from functools import reduce
import pandas as pd

class Integration:
    def __init__(self, a, b ,n,func = None, rule='simpson', tol=None):
        self.a = a # start of interval
        self.b = b # end of interval
        self.n = n # number of partitions
        self.func = func
        self.rule = rule
        self.tol = tol

    def choose_rule(self):
        if self.rule == 'midpoint':
            return Integration.midpoint_rule
        elif self.rule == 'trapezoidal':
            return Integration.trapezoidal_rule
        else:
            return Integration.simpson_rule

    def midpoint_rule(self):
        h = (self.b - self.a) / self.n
        I_midpoint = 0.0
        for i in range(1, self.n+1):
            I_midpoint = I_midpoint + self.func(self.a + (i - 0.5)*h)
        return h * I_midpoint

    def trapezoidal_rule(self):
        h = (self.b - self.a) / self.n
        I_trap = self.func(self.a)/2 + self.func(self.b)/2
        for i in range(1, self.n):
            I_trap = I_trap + self.func(self.a + i * h)
        return h * I_trap

    def simpson_rule(self):
        h = (self.b - self.a) / self.n
        I_simpson = self.func(self.a)/6 + self.func(self.b)/6
        for i in range(1, self.n):
            I_simpson = I_simpson + self.func(self.a + i*h)/3
        for i in range(1,self.n + 1):
            I_simpson = I_simpson + 2 * self.func(self.a + (i - 0.5)* h)/3
        return h * I_simpson

    def given_tolerance(self):
        I_numerical = Integration.choose_rule(self)
        I_old = I_numerical(self)
        self.n = self.n * 2
        I_new = I_numerical(self)

        while (abs(I_new - I_old) > self.tol):
            I_old = I_new
            self.n = 2*self.n
            I_new = I_numerical(self)
        return I_new

######### RUN HERE #######################################################################################
# Change inputs: start, end, # periods, function
# print("Midpoint Rule:", Integration(0,2,4, lambda x: math.exp(-x**2)).midpoint_rule())
# print("Trapezoidal Rule:",Integration(0,2,4, lambda x: math.exp(-x**2)).trapezoidal_rule())
# print("Simpson Rule:",Integration(0,2,32, lambda x: math.exp(-x**2)).simpson_rule())
# print(Integration(0,2,4, lambda x: math.exp(-x**2), rule='trapezoidal',tol=5*10**-7).given_tolerance())
##########################################################################################################


class Bond_math(Integration):
    def __init__(self, t_cashflow=None , v_cashflow= None, r_zero = None,
                 r_inst = None, rule = 'simpson', vec_tol = None, y=None, bond_price = 100):
        self.n = len(t_cashflow)
        self.t_cashflow = t_cashflow
        self.v_cashflow = v_cashflow
        self.r_zero = r_zero
        self.r_inst = r_inst
        self.rule = rule
        self.vec_tol = vec_tol
        self.y = y
        self.bond_price = bond_price

    def bond_price_zero_rates(self):
        B = 0.0
        for i in range(0,self.n):
            disc_i = math.exp(-self.t_cashflow[i]*self.r_zero(self.t_cashflow[i]))
            B = B + self.v_cashflow[i] * disc_i

        # single line way to do it
        # B = reduce(lambda x,y: x+y, list(map(lambda v, t: v* math.exp(-t*self.r_zero(t)),self.v_cashflow,self.t_cashflow)))
        return B

    def bond_price_inst_rt_curve(self):
        B = 0.0
        for i in range(0,self.n):
            I_numerical_i = Integration(a = 0, b = self.t_cashflow[i], n = self.n, func= self.r_inst,
                                        rule='simpson', tol = self.vec_tol[i]).given_tolerance()
            disc_i = math.exp(-I_numerical_i)
            B += self.v_cashflow[i] * disc_i
            print(i, I_numerical_i, disc_i, B, self.t_cashflow[i], self.v_cashflow[i])
        return B

    def bond_dur_convexity(self):
        B, D, C = 0.0, 0.0, 0.0
        for i in range(0,self.n):
            disc_i = math.exp(-self.t_cashflow[i] * self.y)
            B += self.v_cashflow[i] * disc_i
            D += self.t_cashflow[i] * self.v_cashflow[i] * disc_i
            C += (self.t_cashflow[i]**2) * self.v_cashflow[i] * disc_i
        return B, D/B, C/B

    def bond_yield(self):
        x_0 = 0.1 # initial guess 10% yield
        x_new = x_0
        x_old = x_0 - 1
        tol = 1e-06
        while (abs(x_new - x_old) > tol):
            x_old = x_new
            func =sum([v*math.exp(-x_old*t) for v, t in zip(self.v_cashflow, self.t_cashflow)])-self.bond_price
            func_prime = sum([t*v*math.exp(-x_old*t) for v, t in zip(self.v_cashflow, self.t_cashflow)])
            x_new = x_old + func/func_prime
        return x_new



######### RUN HERE #######################################################################################
# print(Bond_math(t_cashflow = [2/12, 8/12, 14/12, 20/12],
#                 v_cashflow = [3,3,3,103],
#                 r_zero = lambda t : 0.0525 + np.log(1+2*t)/200,
#                 ).bond_price_zero_rates())

# print(Bond_math(t_cashflow = [2/12, 8/12, 14/12, 20/12],
#                 v_cashflow = [3, 3, 3, 103],
#                 r_inst = lambda t: 0.0525 + np.log(1+2*t)/200 + t/(100*(1+2*t)),
#                 vec_tol = [1e-04,1e-04,1e-04,1e-06],
#                 rule='simpson'
#                 ).bond_price_inst_rt_curve())

# print("(Bond price, duration, convexity): \n",Bond_math(t_cashflow = [2/12, 8/12, 14/12, 20/12],
#                 v_cashflow = [3,3,3,103],
#                 y = 0.065
#                 ).bond_dur_convexity())

# print(Bond_math(t_cashflow = [4/12, 10/12, 16/12, 22/12, 28/12, 34/12],
#                 v_cashflow = [4,4,4,4,4,104],
#                 bond_price = 105).bond_yield())

# lambda t: 0.0525 + 1/(100*(1+1+np.exp(-t**2)))
# lambda t: 0.0525 + np.log(1+2*t)/200 + t/(100*(1+2t))
##########################################################################################################

class BSM:
    def __init__(self,t=0.0, S=0.0, K=0.0, T=0.0, sigma=0.25, r=0.0, q=0.0, C = 0.0):
        self.t = t
        self.S = S
        self.K = K
        self.T = T
        self.sigma = sigma
        self.r = r
        self.q = q
        self.C = C

    @property
    def d1(self):
        return (np.log(self.S / self.K) + (self.r - self.q + (self.sigma ** 2) / 2) * (self.T - self.t)) / (self.sigma * math.sqrt(self.T - self.t))

    @property
    def d2(self):
        return self.d1 - self.sigma * math.sqrt(self.T - self.t)


    @staticmethod
    def cum_dist_normal(x):
        # computes approximation N(t) for cumulative distribution of Z
        z = abs(x)
        y = 1/(1+0.2316419*z)
        a1 = 0.319381530
        a2 = -0.356563782
        a3 = 1.781477937
        a4 = -1.821255978
        a5 = 1.330274429
        m = 1 - np.exp(-x**2/2)*(a1*y + a2*y**2 + a3*y**3 +
                                 a4*y**4 + a5*y**5)/math.sqrt(2*math.pi)
        if x > 0:
            return m
        else:
            return 1-m


    def call_price(self):
        return self.S*np.exp(-self.q*(self.T-self.t))*BSM.cum_dist_normal(self.d1) - \
               self.K*np.exp(-self.r*(self.T-self.t))*BSM.cum_dist_normal(self.d2)

    def put_price(self):
        return self.K*np.exp(-self.r*(self.T-self.t))*BSM.cum_dist_normal(-self.d2) - \
               self.S*np.exp(-self.q*(self.T-self.t))*BSM.cum_dist_normal(-self.d1)

    # call price's sensitivity to changes in underlying price
    def call_delta(self):
        return math.exp(-self.q*(self.T-self.t))*BSM.cum_dist_normal(self.d1)

    def put_delta(self):
        return -math.exp(-self.q*(self.T-self.t))*BSM.cum_dist_normal(-self.d1)

    # call delta sensitivity to changes in underlying price
    def call_gamma(self):
        return math.exp(-self.q*(self.T-self.t))/(self.S*self.sigma*math.sqrt(self.T-self.t))* \
               (math.exp(-(self.d1**2)/2)/math.sqrt(2*math.pi))

    def put_gamma(self):
        return BSM.call_gamma(self)

    # call price sensitivity to underlying price volatility
    def call_vega(self):
        return self.S*math.exp(-self.q * (self.T - self.t))*math.sqrt(self.T-self.t)* \
               (math.exp(-(self.d1**2)/2)/math.sqrt(2*math.pi))

    def put_vega(self):
        return BSM.call_vega(self)

    # change in call value to time
    def call_theta(self):
        return (-self.S*self.sigma*math.exp(-self.q*(self.T-self.t))*math.exp(-(self.d1**2)/2))/ \
               (2*math.sqrt(2*math.pi*(self.T-self.t))) + \
               self.q*self.S*math.exp(-self.q*(self.T-self.t))*BSM.cum_dist_normal(self.d1) - \
               self.r* self.K*math.exp(-self.r*(self.T-self.t))*BSM.cum_dist_normal(self.d2)

    def put_theta(self):
        return (-self.S*self.sigma*math.exp(-self.q*(self.T-self.t))*math.exp(-(self.d1**2)/2))/ \
               (2*math.sqrt(2*math.pi*(self.T-self.t))) - \
               self.q*self.S*math.exp(-self.q*(self.T-self.t))*BSM.cum_dist_normal(-self.d1) + \
               self.r* self.K*math.exp(-self.r*(self.T-self.t))*BSM.cum_dist_normal(-self.d2)

    # change in call value to interest rates
    def call_rho(self):
        return self.K*(self.T-self.t)*math.exp(-self.r*(self.T-self.t))*BSM.cum_dist_normal(self.d2)

    def put_rho(self):
        return -self.K*(self.T-self.t)*math.exp(-self.r*(self.T-self.t))*BSM.cum_dist_normal(-self.d2)

    def bsm_table(self):
        df = pd.DataFrame([[BSM.call_price(self),BSM.put_price(self)],
                           [BSM.call_delta(self),BSM.put_delta(self)],
                           [BSM.call_gamma(self),BSM.put_gamma(self)],
                           [BSM.call_vega(self),BSM.put_vega(self)],
                           [BSM.call_theta(self),BSM.put_theta(self)],
                           [BSM.call_rho(self),BSM.put_rho(self)]],
                          columns=['Call','Put'],
                          index=['Price','Delta','Gamma','Vega','Theta','Rho'])
        print(df)

    def implied_vol(self):
        x_0 = 0.25
        x_new = x_0
        x_old = x_0 -1
        tol = 1e-06
        while (abs(x_new - x_old) > tol):
            x_old = x_new
            self.sigma = x_new
            x_new = x_new - (BSM.call_price(self) - self.C)/BSM.call_vega(self)
        return x_new


######### RUN HERE ##########################################################################
# print(BSM.cum_dist_normal(x=.0040123546))
# print(BSM(S=139, K=135, T=0.5, t=0, sigma=0.3, q=0.0377, r=0.01).bsm_table())
# print(BSM(S=139, K=135, T=0.5, t=0, sigma=0.3, q=0.0377, r=0.01).call_price())
print(BSM(C=7 ,S=25, K=20, T=1, t=0, q = 0.00, r=0.05).implied_vol())
#############################################################################################

class Numerical_methods:
    def __init__(self, x_0 = None, x_neg1 = None, left= None, right = None, func = None, func_prime = None,
                 tol_int = None, tol_approx = None, tol_consec = None):
        self.x_0 = x_0
        self.x_neg1 = x_neg1
        self.left = left
        self.right = right
        self.func = func
        self.func_prime = func_prime
        self.tol_int = tol_int
        self.tol_approx = tol_approx
        self.tol_consec = tol_consec

    def bisection_method(self):
        """cuts the interval in half and then interatively sees if there's a sign change to solve
        for x for f(x)=0"""
        x_l = self.left
        x_r = self.right
        while (max(abs(self.func(x_l)),abs(self.func(x_r))) > self.tol_approx) or (x_r - x_l > self.tol_int):
            x_m = (x_l + x_r)/2
            if (self.func(x_l)*self.func(x_m) < 0):
                x_r = x_m
            else:
                x_l = x_m
        return x_m

    def newtons_method(self):
        """x intercept of the tangent line (to func_prime) of the graph of f(x) at the point x_k"""
        x_new = self.x_0
        x_old = self.x_0-1
        while (abs(self.func(x_new)) > self.tol_approx) or (abs(x_new - x_old) > self.tol_consec):
            x_old = x_new
            x_new = x_old - self.func(x_old)/self.func_prime(x_old)
        return x_new

    def secant_method(self):
        """use over newton's method when you can't solve for func_prime directly"""
        x_new = self.x_0
        x_old = self.x_neg1
        while (abs(self.func(x_new)) > self.tol_approx) or (abs(x_new - x_old) > self.tol_consec):
            x_oldest = x_old
            x_old = x_new
            x_new = x_old - self.func(x_old)*((x_old-x_oldest)/(self.func(x_old)-self.func(x_oldest)))
        return x_new


######### RUN HERE ##########################################################################
# print(Numerical_methods(left=-2,
#                         right=3,
#                         func= lambda x: x**4 - 5*x**2+4-(1/(1+math.exp(x**3))),
#                         tol_int = 1e-06,
#                         tol_approx=1e-09).bisection_method())
#
# print(Numerical_methods(x_0 = -0.5, # initial guess
#                         func= lambda x: x**4 - 5*x**2+4-(1/(1+math.exp(x**3))),
#                         func_prime = lambda x: 4*x**3-10*x+(3*x**3*math.exp(x**3)/(1+math.exp(x**3))**2),
#                         tol_consec = 1e-06,
#                         tol_approx=1e-09).newtons_method())
#
# print(Numerical_methods(x_0 = -3, # initial guess
#                         x_neg1 = -3-.01, # something close to initial guess
#                         func= lambda x: x**4 - 5*x**2+4-(1/(1+math.exp(x**3))),
#                         tol_consec = 1e-06,
#                         tol_approx=1e-09).secant_method())
##############################################################################################