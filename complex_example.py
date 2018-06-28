#!/usr/bin/python

import matplotlib.pyplot as plt
from scipy.integrate import ode

def f(t, y, arg1, arg2):
    return -arg1*y -arg2*y

#Parameters
delta = 0.0
d_delta = 0.0
tau = 1000
Omega0 = 1.0j

f_arg = (1.0j * delta) - 1/tau
jac_arg = ((abs(Omega0**2))/4) + 1.0j*(d_delta - (delta/tau))

r = ode(f).set_integrator('zvode')
r.set_initial_value(1.0,1.0).set_f_params(f_arg,jac_arg)

tl = 20
dt = 0.001

t_vals = []
u_vals = []

while r.successful() and r.t < tl:
    r.integrate(r.t+dt)
    t_vals.append(r.t)
    u_vals.append(r.y)
    #print("%g %g" % (r.t, r.y))

plt.plot(t_vals, u_vals)
plt.show()
