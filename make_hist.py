import numpy as np
import matplotlib.pyplot as plt

decay1 = np.empty([4500])
R01 = np.empty([4500])
tinf1 = np.empty([4500])
decay2 = np.empty([4500])
R02 = np.empty([4500])
tinf2 = np.empty([4500])
decay3 = np.empty([4500])
R03 = np.empty([4500])
tinf3 = np.empty([4500])
decay4 = np.empty([4500])
R04 = np.empty([4500])
tinf4 = np.empty([4500])
decay5 = np.empty([4500])
R05 = np.empty([4500])
tinf5 = np.empty([4500])

a = np.loadtxt('ferret2_0.dat')
for i in np.arange(4500):
    decay1[i] = min(a[i,2:4])
    R01[i] = a[i,0]*a[i,1]/(a[i,2]*a[i,4])
    tinf1[i] = 24.*(2./(a[i,0]*a[i,1]))**(0.5)

a = np.loadtxt('ferret2_3.dat')
for i in np.arange(4500):
    decay2[i] = min(a[i,2:4])
    R02[i] = a[i,0]*a[i,1]/(a[i,2]*a[i,4])
    tinf2[i] = 24.*(2./(a[i,0]*a[i,1]))**(0.5)

a = np.loadtxt('ferret2_7.dat')
for i in np.arange(4500):
    decay3[i] = min(a[i,2:4])
    R03[i] = a[i,0]*a[i,1]/(a[i,2]*a[i,4])
    tinf3[i] = 24.*(2./(a[i,0]*a[i,1]))**(0.5)

a = np.loadtxt('ferret2_14.dat')
for i in np.arange(4500):
    decay4[i] = min(a[i,2:4])
    R04[i] = a[i,0]*a[i,1]/(a[i,2]*a[i,4])
    tinf4[i] = 24.*(2./(a[i,0]*a[i,1]))**(0.5)

a = np.loadtxt('ferret2_28.dat')
for i in np.arange(4500):
    decay5[i] = min(a[i,2:4])
    R05[i] = a[i,0]*a[i,1]/(a[i,2]*a[i,4])
    tinf5[i] = 24.*(2./(a[i,0]*a[i,1]))**(0.5)

fig, ax = plt.subplots()
ax.hist(decay1, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'red', alpha = 0.3)
ax.hist(decay2, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'blue', alpha = 0.3)
ax.hist(decay3, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'green', alpha = 0.3)
ax.hist(decay4, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'turquoise', alpha = 0.3)
ax.hist(decay5, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'orange', alpha = 0.3)
ax.legend(['Day 0','Day 3', 'Day 7', 'Day 14', 'Day 28'], fontsize=22)
ax.set_xlabel(r'Decay rate (\d)', fontsize=22)
plt.yticks(fontsize=18 )
plt.xticks(np.arange(0,21,5), ('0','5','10','15','20'),fontsize=18 )
plt.savefig("ferret2_decay.pdf", bbox_inches='tight', pad_inches=0.2)

fig, ax = plt.subplots()
ax.hist(R01, bins=100, range = (0,100), density=True, histtype='barstacked', color = 'red', alpha = 0.3)
ax.hist(R02, bins=100, range = (0,100), density=True, histtype='barstacked', color = 'blue', alpha = 0.3)
ax.hist(R03, bins=100, range = (0,100), density=True, histtype='barstacked', color = 'green', alpha = 0.3)
ax.hist(R04, bins=100, range = (0,100), density=True, histtype='barstacked', color = 'turquoise', alpha = 0.3)
ax.hist(R05, bins=100, range = (0,100), density=True, histtype='barstacked', color = 'orange', alpha = 0.3)
ax.legend(['Day 0','Day 3', 'Day 7', 'Day 14', 'Day 28'], fontsize=22)
ax.set_xlabel(r'$R_0$', fontsize=22)
plt.yticks(fontsize=18 )
plt.xticks(np.arange(0,101,20), ('0','20','40','60','80','100'),fontsize=18 )
plt.savefig("ferret2_R0.pdf", bbox_inches='tight', pad_inches=0.2)

fig, ax = plt.subplots()
ax.hist(tinf1, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'red', alpha = 0.3)
ax.hist(tinf2, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'blue', alpha = 0.3)
ax.hist(tinf3, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'green', alpha = 0.3)
ax.hist(tinf4, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'turquoise', alpha = 0.3)
ax.hist(tinf5, bins=100, range = (0,20), density=True, histtype='barstacked', color = 'orange', alpha = 0.3)
ax.legend(['Day 0','Day 3', 'Day 7', 'Day 14', 'Day 28'], fontsize=22)
ax.set_xlabel(r'$t_{inf}$ (h)', fontsize=22)
plt.yticks(fontsize=18 )
plt.xticks(np.arange(0,21,5), ('0','5','10','15','20'),fontsize=18 )
plt.savefig("ferret2_tinf.pdf", bbox_inches='tight', pad_inches=0.2)
