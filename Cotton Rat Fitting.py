# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:20:28 2017

@author: sams club
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 17:37:11 2017

@author: sams club
"""


import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

params = np.ones(5)
# Parameter array

# data = np.loadtxt("Lung 3 Day Old.csv", delimiter = ",")
# data = np.loadtxt("Lungs_14_Day_Old.csv", delimiter = ",")
data = np.loadtxt("Lung 28 Day Old.csv", delimiter = ",")
# data = np.loadtxt("Lung Adult .csv", delimiter = ",")

# data = np.loadtxt("Trachea 3 Day Old.csv", delimiter = ",")
# data = np.loadtxt("Trachea 14 Day Old.csv", delimiter = ",")
# data = np.loadtxt("Trachea 28 Day Old.csv", delimiter = ",")
# data = np.loadtxt("Trachea Adult Infection.csv", delimiter = ",")

# data = np.loadtxt("Nose 3 Day Old Infection.csv", delimiter = ",")
# data = np.loadtxt("Nose 14 Day Old.csv", delimiter = ",")
# data = np.loadtxt("Nose 28 Day Old.csv", delimiter = ",")
# data = np.loadtxt("Nose Adult.csv", delimiter = ",")

# Importing data from a spreadsheet into a variable to used to call
# the data set. 


x_data = data[:,][:,0]
# Initializing x variable to store 1st column of data, time.

y_data = data[:,][:,1]
# Initializing y variable to store 2nd column of data, viral count.

#model for reference
# y = T,E,I,V
# dydt[0] = dT
# dydt[1] = dE
# dydt[2] = dI
# dydt[3] = dV
# y[0] = T
# y[1] = E
# y[2] = I
# y[3] = V
# z[0] = B
# z[1] = k
# z[2] = d
# z[3] = p
# z[4] = c

def F(y,t):
    """
        Defining function which uses simple viral model where:
            dT/t = -BTV
            dE/t = BTV - kE
            dI/t = kE - dI
            dV/t = pI- cV
        The dydt array contains the values for the variables in order:
        (dT/t,dE/t,dI/t,and dV/t)
        The y array contains the values for the parameters in order:
        (B,k,d,p, and c)
        The z array contains values for the parameters in order:
        (T,E,I, and V)
        Needed for use in odeint command.
    """
    global params

    dydt = np.zeros_like(y)
    dydt[0] = -params[0] * y[0] * y[3]
    dydt[1] = params[0] * y[0] * y[3] - params[1] * y[1]
    dydt[2] = params[1] * y[1] - params[2] * y[2]
    dydt[3] = params[3] * y[2] - params[4] * y[3]
    return dydt


t = np.linspace(0, 7, 600)
# Initializing time variable for a time span of 20 days beginning on day 1.

guess = np.array( [1e-5, 4, 4, 4e6, 5, 1e4] )
# Initializing guesses for parameters (pulled from handout on simple viral model)
# in order: B, k, d, p, c, and v0.

def SSR(guess):
    """
        SSR(sum of squared residuals) is the objective function that will
        later be minimized. The log of the integrated model and the original
        data set are taken so to account for the logarithmic nature of the
        data and the optimizer's tendency to overcompensate for large 
        differences.
        SSR Formula = sum(y_prdicted - y_data)^2
    """
    global params
    
    params = guess[0:5]
    y0 = [1, 0, 0, guess[5]]
    # Initializing initial conditions- using 1st data point of RS virus data set 
    # as the initial value of virus. Initial conditions in order: T, E, I, V.
    
    y_prediction = odeint(F,y0,x_data)
    virus_prediction = (y_prediction[:,][:,3])
    virus_prediction_log = np.log10(virus_prediction)
    virus_data_log = np.log10(y_data)
    #plt.plot(virus_prediction_log)
    #plt.plot(virus_data_log,'r')
    #plt.show(True)
    #Prints plotting out in real time. 
    diff = virus_prediction_log - virus_data_log
    return sum((diff)**2)

# result = minimize(SSR,guess, method = 'Nelder-Mead')
# result = minimize(SSR,guess, method = 'Powell')
# result = minimize(SSR,guess, method = 'CG')
# result = minimize(SSR,guess, method = 'BFGS')
# result = minimize(SSR,guess, method = 'Newton-CG')
result = minimize(SSR,guess, method = 'L-BFGS-B')
# result = minimize(SSR,guess, method = 'TNC')
# result = minimize(SSR,guess, method = 'COBYLA')
# result = minimize(SSR,guess, method = 'SLSQP')
# result = minimize(SSR,guess, method = 'dogleg')
# result = minimize(SSR,guess, method = 'trust-ncg')
print (result)
print (result.success)
# Finding the smallest SSR value and the corresponding parameter values.
# Printed for later reference.

params = result.x[0:5]
# Initializing an array to hold results of parameters obtained
# by optimiziation in order of : B, k, d, p, and c.

y0 = [1, 0, 0, result.x[5]]
y_fitted = odeint(F,y0,t)
virus_best_fit = y_fitted[:,][:,3]

# Obtaining predicted y values when the optimized parameters are used
# in the model. Virus_best_fit var is used to store the virus model as
# its in the 3rd index and listed 4th in T, E,I, and V.

optimal_B_value = result.x[0]
optimal_k_value = result.x[1]
optimal_d_value = result.x[2]
optimal_p_value = result.x[3]
optimal_c_value = result.x[4]
# Setting the parameters to the the returns of the optimized SSR.

print ("Optimal B value = ", optimal_B_value)
print ("Optimal k value = ", optimal_k_value)
print ("Optimal d value = ", optimal_d_value)
print ("Optimal p value = ", optimal_p_value)
print ("Optimal c value = ", optimal_c_value)

plt.plot(x_data,y_data, "o", label = "Dataset")
plt.plot(t,virus_best_fit, label = "Line of Best Fit")
plt.legend()
plt.xlabel("Time (days)")
ax = plt.gca()
ax.set_yscale('log')
plt.ylabel("Virus")
plt.title("RSV Count in Lungs of 28 Day Old Rats")
