#import sys
#sys.path.append('/Users/stakeru/.python/versions/3.11.1/lib/python3.11/site-packages')
import numpy as np
import matplotlib.pyplot as plt

dt = np.dtype([("t","<d"), ("var","<15010d"), ("varave","<19513f")]) # data format per recl
fd = open("printvars32.dat","r")
tframe = 601
chunk = np.fromfile(fd, dtype=dt, count=tframe) #count = the number of time slices

y1ave = np.zeros(1501) # creat 0 array
y2ave = np.zeros(1501) # creat 0 array
#for j in range (0,tframe): 
for j in range (301,400):
    print ("t:",chunk[j]["t"])
    var = chunk[j]["var"].reshape((1501,10),order="F")
    #print "var(1,1):",var[0][0],var[3][0]
    varave = chunk[j]["varave"].reshape((1501,13),order="F")

    x = np.array([r[0] for r in var]) # r=R/R_star
    x = x - 1 # r => r-1
    #x = np.log10(x)
    y1 = np.array([r[2] for r in var]) # density
    y1 = y1*2.5e-7 # simulation unit => g/cm^3
    y1ave = y1ave + y1
    y2 = np.array([r[3] for r in var]) # velocity
    y2 = y2*7.505 # simulation unit => km/s
    y2ave = y2ave + y2
    #y = np.log10(y)
y1ave = y1ave * 0.01
y2ave = y2ave * 0.01
fig = plt.figure(figsize = (10,6), facecolor='lightblue')
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
#density plot----
ax1.plot(x, y1ave)
ax1.loglog(base=10)
ax1.set_xlabel("r-1",fontsize=20)
ax1.set_ylabel("density(g/cm^3)",fontsize=20)
ax1.set_xlim([1.e-4,100])
ax1.set_ylim([1.e-22,5.e-7])
# velocity plot----
ax2.plot(x, y2ave)
#ax2.semilogx(basex=10)
ax2.loglog(base=10)
ax2.set_xlabel("r-1",fontsize=20)
ax2.set_ylabel("velocity(km/s)",fontsize=20)
ax2.set_xlim([1.e-4,100])
ax2.set_ylim([1.,1000.])
#plt.show()
plt.pause(10)
