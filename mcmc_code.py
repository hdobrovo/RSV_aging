# -*- coding: utf-8 -*-
#from __future__ import print_function, division
#rom IPython.display import display, Math

#import os
#import sys
import numpy as np
import corner
from scipy.integrate import odeint
from matplotlib import pyplot as plt

import emcee
global par
global data

#time[0], old_a[1], young_a[2], old_b[3], young_b[4]
global par
data_file = "1a- Cotton Rat Experimental Data Files/Nose 3 Day Old.dat" 
#data_file = "1b- Ferret Experimental Data Files/Open Circle Nose Ferret Day 28.dat"
tfile = np.loadtxt(data_file)
t = tfile[:,0]
#olda = tfile[:,1]
#younga = tfile[:,2]
#oldb = tfile[:,3]
#youngb = tfile[:,4]
vdata = np.log10(tfile[:,1])
print(vdata)

# ODE function describing the virus and cell equations
def virus_dt(Z,t):
    global par
    b = par[0]
    p = par[1]
    c = par[2]
    k = par[3]
    d = par[4]

    [T, I1, I2, V] = [0, 1, 2, 3]
    dT = -b*Z[T]*Z[V]  #target
    dI1 = b*Z[T]*Z[V]- k*Z[I1] #infected
    dI2 = k*Z[I1] - d*Z[I2] #productive
    dV = p*Z[I2] - c*Z[V] #virus

    return [dT, dI1, dI2, dV]

# Integrate the ODE function and return the virus
def solveeqs(pp):
    global par
    global t

    fpar=np.power(10,pp)
    par = []
    par.extend(fpar[:])
    upd = []
 #  [dT, dI1, dI2, dV]
    y = [1.0, 0.0, 0.0, fpar[-1]]
    sample_time=np.arange(8.)
#    sample_time=[0., 2., 3., 4., 5., 6., 7., 8., 10.]
    soln = odeint(virus_dt, y, sample_time)
    Vm = soln[:,-1]
    return Vm

# This is the initial guess 
# parameter models b, p, c, k, d (1/15) v0 
result=[2.15*10**6,5.39,5.39,4.40*10**7,5.40,26.6]

# Define the SSR
def lnlike(theta, t, vdata):
    b,p,c,k,d,v0 = theta
    model = solveeqs(theta)
    model = np.log10(model[1:])
    x = np.where(np.isnan(model))
    model[x] = 0
    print(np.sum((vdata-model)**2))
    return -np.sum((vdata-model)**2)

# Define the range for all the parameters
def lnprior(theta):
    b,p,c,k,d,v0 = theta
    if -8<b< -2 and  4< p < 10 and -3 < c< 2 and -3 < k < 2 and -3 < d < 2 and -3 < v0 < 5:
        return 0.0
    return -np.inf

# Defines the new parameters
def lnprob(theta, t, vdata):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, t, vdata)

print('Start')
# You can change the number of walkers (second number)
ndim, nwalkers = 6, 50
# Sets up all the walkers in a ball around the initial guess
pos = [np.log10(result) + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# Runs the Markov chain Monte Carlo sampling
move = emcee.moves.StretchMove(a=6.0)
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, vdata), moves=[move])
sampler.run_mcmc(pos, 1000); # The second argument is the number of steps taken, you can change this

res=plt.plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
plt.show()
res=plt.plot(sampler.chain[:,:,1].T, '-', color='k', alpha=0.3)
plt.show()
res=plt.plot(sampler.chain[:,:,2].T, '-', color='k', alpha=0.3)
plt.show()
res=plt.plot(sampler.chain[:,:,3].T, '-', color='k', alpha=0.3)
plt.show()
res=plt.plot(sampler.chain[:,:,4].T, '-', color='k', alpha=0.3)
plt.show()
res=plt.plot(sampler.chain[:,:,5].T, '-', color='k', alpha=0.3)
plt.show()
#es.savefig("chain1.png")
af = np.mean(sampler.acceptance_fraction)
print(af)

# Creates the corner plot
samples =sampler.chain[:, 100:, :].reshape((-1, ndim))
print('Done', samples.shape)
samples=samples[0::10]
fig = corner.corner(np.power(10,samples), labels=["$b$", "$p$","$c$", "$k$", "$d$", "$V_0$"], truths=result, plot_contours="False")
fig.savefig("test.pdf")
plt.show()

print('Done')

# Save the results
#fd=open('lung14.dat','w')
#fd.write(str(np.power(10,samples)).tolist())#,str(samples))
#fd.close()
#np.savetxt("nose3.dat",np.power(10,samples))

x=sampler.chain
print('Done', samples.shape, x.shape)
