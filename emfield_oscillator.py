#usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def f(y,t,params):
    u,v = y
    a = 1.0
    delta, d_delta, tau, Omega0 = params
    derivs = [v, (a*delta - 1/tau)*v - (0.25*(Omega0**2) + a*(d_delta
        - delta/tau))*u]
    return derivs

def plotter(time, solution, title):
    fig = plt.figure(1, figsize=(8,8))
    #fig.add_title(title)

    ax1 = fig.add_subplot(311)
    ax1.plot(time, solution[:,0]**2)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("u")
#    ax1.set_xticklabels([])
#    ax1.set_yticklabels([])

    ax2 = fig.add_subplot(312)
    ax2.plot(time, solution[:,1])
    ax2.set_xlabel("Time")
    ax2.set_ylabel("v")
#    ax2.set_xticklabels([])
#    ax2.set_yticklabels([])

    ax3 = fig.add_subplot(313)
    ax3.plot(solution[:,0], solution[:,1], "-", ms=1)
    ax3.set_xlabel("u")
    ax3.set_ylabel("v")
#    ax3.set_xticklabels([])
#    ax3.set_yticklabels([])
#    ax3.axis('off')

    plt.tight_layout()
    plt.show()

#Parameters
delta = 0.0
d_delta = 0.0
tau = 10000
Omega0 = 0.5

#Initial values
u0 = 1.0
v0 = 0.0

#Parameter bundle
params = [delta, d_delta, tau, Omega0]

y0 = [u0, v0]

tStop = 200
tInc = 0.05
t = np.arange(0.0, tStop, tInc)

psoln = odeint(f, y0, t, args=(params,))
plotter(t, psoln, "some string")
