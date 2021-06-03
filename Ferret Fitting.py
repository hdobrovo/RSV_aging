# -*- coding: utf-8 -*-
"""
Created on Wed May 23 12:21:17 2018

@author: sams club
"""

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
from matplotlib import rcParams

"""
HOW TO USE:
    Select the two datafiles you wish to use
    Choose the correct guess
    For the open circle graphs, copy and paste the name of the file that has
    the virus count for the day that you want. 
"""

params = np.ones(5)
# Parameter array

# data = np.loadtxt("Nose Ferrets Day 0.csv", delimiter = ",")
# data = np.loadtxt("New Ferret Nose 3 Closed.csv", delimiter = ",")
# data = np.loadtxt("New Ferret Nose 7 Closed.csv", delimiter = ",")
#data = np.loadtxt("Nose Ferrets Day 14 New.csv", delimiter = ",") 
data = np.loadtxt("New Ferret 28 Closed.csv", delimiter = ",")

# data2 = np.loadtxt("Nose Open Circle 0 Days.csv", delimiter = ",")
# data2 = np.loadtxt("New Ferret Nose 3 Open.csv", delimiter = ",")
# data2 = np.loadtxt("New Ferret Nose 7 Open.csv", delimiter = ",")
# data2 = np.loadtxt("Nose Open Circle 14 Days.csv", delimiter = ",")
# data2 = np.loadtxt("Nose Ferrets Open Circle Day 14 New.csv", delimiter = ",")
data2 = np.loadtxt("New Ferret Nose 28 Open.csv", delimiter = ",")

# data = np.loadtxt("Lung Ferrets Day 0.csv", delimiter = ",")
# data = np.loadtxt("Lung Ferrets Day 3.csv", delimiter = ",")
# data = np.loadtxt("Lung Ferrets Day 7.csv", delimiter = ",")
# data = np.loadtxt("Lung Ferrets Day 14.csv", delimiter = ",")
# data = np.loadtxt("Lung Ferrets Day 28.csv", delimiter = ",")

# data = np.loadtxt("Nose Ferrets Day 3.csv", delimiter = ",")
# data = np.loadtxt("Nose Ferrets Day 7.csv", delimiter = ",")
# data = np.loadtxt("Nose Ferrets Day 14.csv", delimiter = ",")
# data = np.loadtxt("Nose Ferrets Day 28.csv", delimiter = ",")
# data = np.loadtxt("Nose Open Circle 3 Days.csv", delimiter = ",")
# data = np.loadtxt("Nose Open Circle 7 Days.csv", delimiter = ",")
# data = np.loadtxt("Nose Open Circle 28 Days.csv", delimiter = ",")

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

t = np.linspace(0, 10, 1000)
# for day 0, day 7, day 14, and day 28

# t = np.linspace(0, 8, 600)
# for day 3 

guess = np.array( [1e-6, 4, 4, 1e8, 4, 15] ) 
# Initializing guesses for parameters (pulled from handout on simple viral model)
# in order: B, k, d, p, c, and v0.
##################### NOSE Ferret 1- Closed Circles #####################
# Nose 0 Day
# guess = np.array( [5e-7, 2, 3, 1e8, 2, 1e1] ) 
# SSR = 7.9987638534

# Nose 3 Day
# guess = np.array( [5e-7, 7, 2, 1e8, 10, 1e1] ) 
# NEW ONE
# guess = np.array( [5e-7, 2, 7, 1e8, 10, 1e1] ) 
# SSR = 3.74794799735

# Nose 7 Day
# guess = np.array( [5e-7, 2, 3, 1e8, 2, 1e1] )
# NEW ONE
# guess = np.array( [1e-6, 4, 4, 1e8, 4, 15] ) 
# SSR = 9.18460274775

# Nose 14 Day
# guess = np.array( [1e-6, 4, 1, 1e8, 4, 15] ) 
# NEW
# guess = np.array( [1e-6, 4, 4, 1e8, 4, 15] ) 
# SSR = 2.13112773996

# Nose 28 Day
# guess = np.array( [1e-6, 4, 4, 1e8, 4, 15] ) 
# SSR = 4.63185820437
############################################################################

##################### NOSE Ferret 2- Open Circles #####################
# Nose Open 0 Day
# guess = np.array( [1e-6, 10, 4, 1e8, 5, 1e1] )  
# SSR = 2.82472488319

# Nose Open 3 Day
# guess = np.array( [1e-6, 7, 7, 1e8, 7, 15] )  #NEW
# guess = np.array( [1e-6, 1, 1, 1e6, 1, 15] ) 
# SSR = 5.03112868245

# Nose Open 7 Day
# guess = np.array( [1e-6, 7, 7, 1e8, 7, 15] ) 
# NEW
# guess = np.array( [1e-7, 10, 4, 9e8, 12, 15] ) 
# SSR = 13.0907003303

# Nose Open 14 Day
# guess = np.array( [1e-6, 4, 1, 1e8, 4, 15] ) 
# NEW
# guess = np.array( [1e-7, 10, 4, 9e8, 12, 15] ) 
# SSR = 3.07872273134

# Nose Open 28 Day
# guess = np.array( [1e-6, 4, 4, 1e8, 4, 15] ) 
# SSR = 
#######################################################################
################### LUNG BELOW ###################
# Lung 0 Day
# guess = np.array( [1e-6, 7, 4, 1e8, 4, 15] )  
# SSR = 1.66096821826
    
# Lung 3 Day
# guess = np.array( [1e-6, 4, 2, 1e8, 4, 15] ) 
# SSR = 0.780053215282

# Lung 7 Day
# guess = np.array( [1e-6, 4, 4, 1e8, 4, 15] ) 
# SSR = 

# Lung 14 Day
# guess = np.array( [1e-6, 4, 2, 1e8, 4, 15] )   
# SSR = 7.22095865039

# Lung 28 Day
# guess = np.array( [1e-7, 10, 4, 1e8, 12, 15] ) 
# SSR = 7.86697095542

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
    
    times=[0.]
    params = guess[0:5]
    y0 = [1, 0, 0, guess[5]]
    # Initializing initial conditions- using 1st data point of RS virus data set 
    # as the initial value of virus. Initial conditions in order: T, E, I, V.
    times=np.append(times,x_data)
    y_prediction = odeint(F,y0,times)
    virus_prediction = (y_prediction[:,][1:,3])
    virus_prediction_log = np.log10(virus_prediction)
    virus_data_log = np.log10(y_data)
    #plt.plot(virus_prediction_log)
    #plt.plot(virus_data_log,'r')
    #plt.show(True)
    #Prints plotting out in real time. 
    diff = virus_prediction_log - virus_data_log
    return sum((diff)**2)

result = minimize(SSR,guess, method = 'Nelder-Mead', options={'maxfev':10000})
# result = minimize(SSR,guess, method = 'Powell')
# result = minimize(SSR,guess, method = 'CG')
# result = minimize(SSR,guess, method = 'BFGS')
# result = minimize(SSR,guess, method = 'Newton-CG')
# result = minimize(SSR,guess, method = 'L-BFGS-B')
# result = minimize(SSR,guess, method = 'TNC')
# result = minimize(SSR,guess, method = 'COBYLA')
# result = minimize(SSR,guess, method = 'SLSQP', bounds = bnds)
# result = minimize(SSR,guess, method = 'dogleg')
# result = minimize(SSR,guess, method = 'trust-ncg')
print (result)
#print ( result.success)
#print (result.fun)
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
optimal_v0_value = result.x[5]

# Setting the parameters to the the returns of the optimized SSR.

print ("Optimal B value = ", optimal_B_value)
print ("Optimal k value = ", optimal_k_value)
print ("Optimal d value = ", optimal_d_value)
print ("Optimal p value = ", optimal_p_value)
print ("Optimal c value = ", optimal_c_value)
print ("Optimal v0 value = ", optimal_v0_value)
print (result.fun)
print (result.success)

rfont = {'fontname':'Times New Roman'}

plt.plot(x_data,y_data, "o", label = "Dataset", color = "turquoise")
plt.plot(data2[:,][:,0], data2[:,][:,1], "o", label = "Dataset", color = "magenta")

# plt.plot(t,virus_best_fit, label = "Line of Best Fit", color="purple")
# Nose Rat 1
Open_Cirle_Virus = np.loadtxt("Virus Best Fit Open.csv", delimiter = ",")

virus_best_fit_0_day = Open_Cirle_Virus[:,][:,0]
virus_best_fit_3_day = Open_Cirle_Virus[:,][:,1]
virus_best_fit_7_day = Open_Cirle_Virus[:,][:,2]
virus_best_fit_14_day = Open_Cirle_Virus[:,][:,3]
virus_best_fit_28_day = Open_Cirle_Virus[:,][:,4]

plt.plot(t,virus_best_fit, label = "Line of Best Fit", color="turquoise")
plt.plot(t,virus_best_fit_28_day, label = "Line of Best Fit", color="magenta")
# Nose Rat 2

ax = plt.gca()
plt.xlabel("Time (days)",fontsize = 25,**rfont,)
ax.set_yscale('log',fontsize = 25,)
plt.ylabel("Virus Titer (PFU/g)", fontsize = 25,**rfont,) #family = 'serif', variant = 'normal')
plt.xticks(fontsize = 25,**rfont,)
plt.yticks(fontsize = 25,**rfont,)
rcParams.update({'figure.autolayout': True})
# ax.set_xticklabels(ax.get_xticks(), fontsize = 25, family = 'serif', variant = 'normal')
# ax.set_yticklabels(ax.get_yticks(), fontsize = 25, family = 'serif', variant = 'normal')

# plt.title("0 Day Old",fontsize = 30,**rfont,)
# plt.title("3 Day Old",fontsize = 30, **rfont,)
# plt.title("7 Day Old",fontsize = 30,**rfont,)
# plt.title("14 Day Old",fontsize = 30,**rfont,)
plt.title("28 Day Old",fontsize = 30,**rfont,)

# plt.savefig('0_day_nose_f.pdf', bbox_inches='tight')
# plt.savefig('3_day_nose_f.pdf', bbox_inches='tight')
# plt.savefig('7_day_nose_f.pdf', bbox_inches='tight')
# plt.savefig('14_day_nose_f.pdf', bbox_inches='tight')
plt.savefig('28_day_nose_f.pdf', bbox_inches='tight')

# plt.savefig('Nose_0_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_3_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_7_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_14_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_28_Day_Ferret_Fit.pdf', bbox_inches='tight')

# plt.savefig('New_Nose_0_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_3_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_7_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_14_Day_Ferret_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_28_Day_Ferret_Fit.pdf', bbox_inches='tight')

# plt.savefig('Nose_0_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_3_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_3_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_7_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_7_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_14_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_14_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_28_Day_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_28_Day_Fit.pdf', bbox_inches='tight')

# plt.savefig('Nose_Open_Circle_0_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_Open_Circle_3_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_Open_Circle_3_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_Open_Circle_7_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_Open_Circle_7_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_Open_Circle_14_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_Open_Circle_14_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('Nose_Open_Circle_28_Days_Fit.pdf', bbox_inches='tight')
# plt.savefig('New_Nose_Open_Circle_28_Days_Fit.pdf', bbox_inches='tight')

"""
plt.plot(x_data,y_data, "o", label = "Dataset", color = "black")

# plt.plot(t,virus_best_fit, label = "Line of Best Fit", color="red")
#Lung

# plt.plot(t,virus_best_fit, label = "Line of Best Fit", color="blue")
#Trachea

plt.plot(t,virus_best_fit, label = "Line of Best Fit", color="green")
#Nose

plt.xlabel("Time (days)",fontsize = 25,)
ax = plt.gca()
ax.set_yscale('log',fontsize = 25,)
plt.ylabel("Virus Titer (PFU/g)", fontsize = 25,) #family = 'serif', variant = 'normal')
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)

#Lung_3_Day_OldFile.close()
df = pd.DataFrame({'B Value' : boot_optimal_B_value_data_storage, 
                   'k Value' : boot_optimal_k_value_data_storage,
                   'd Value' : boot_optimal_d_value_data_storage,
                   'p Value' : boot_optimal_p_value_data_storage,
                   'c Value' : boot_optimal_c_value_data_storage,
                   'v0 Value' : boot_optimal_v0_value_data_storage,
                   'SSR Value': boot_optimal_SSR_value_data_storage})
    
# df.to_csv('Lung 3 Day Old.txt', sep='\t')
# df.to_csv('Lung 14 Day Old.txt', sep='\t')
# df.to_csv('Lung 28 Day Old.txt', sep='\t')
# df.to_csv('Lung Adult.txt', sep='\t')

# df.to_csv('Trachea 3 Day Old.txt', sep='\t')
# df.to_csv('Trachea 14 Day Old.txt', sep='\t')
# df.to_csv('Trachea 28 Day Old.txt', sep='\t')
# df.to_csv('Trachea Adult.txt', sep='\t')

# df.to_csv('Nose 3 Day Old.txt', sep='\t')
# df.to_csv('Nose 14 Day Old.txt', sep='\t')
# df.to_csv('Nose 28 Day Old.txt', sep='\t')
# df.to_csv('Nose Adult.txt', sep='\t')
"""
