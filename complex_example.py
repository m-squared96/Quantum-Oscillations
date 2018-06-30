#!/usr/bin/python

from numpy import pi
import matplotlib.pyplot as plt
from scipy.integrate import ode

def f(t, y, arg1, arg2):
    u, v = y
    return [v, arg1*v - arg2*u]

#Parameters
delta = 0.5
d_delta = 0.0
tau = 1e9
Omega0 = 1.0j/2*pi

f_arg = complex(-1/tau, delta)
jac_arg = complex(abs(Omega0**2)/4, d_delta - (delta/tau))

#f_arg = (1.0j * delta) - 1/tau
#jac_arg = ((abs(Omega0**2))/4) + 1.0j*(d_delta - (delta/tau))

r = ode(f).set_integrator('zvode')
r.set_initial_value([0.0, 1.0],0.0).set_f_params(f_arg,jac_arg)

tl = 100
dt = 0.01

t_vals = []
u_vals = []
usquare_vals = []
v_vals = []
vsquare_vals = []

while r.successful() and r.t < tl:
    r.integrate(r.t+dt)
    t_vals.append(r.t)

    u_vals.append(r.y[1])
    usquare_vals.append((abs(r.y[1]))**2)
    
    v_vals.append(r.y[0])
    vsquare_vals.append((abs(r.y[0]))**2)

plt.figure()
plt.plot(t_vals, u_vals, label="u(t)")
#plt.plot(t_vals, v_vals, label="u'(t)")
plt.xlabel("Time")
plt.legend()
plt.title("u(t) and u'(t)")

plt.figure()
plt.plot(t_vals, usquare_vals)
#plt.xlabel("Time")
#plt.title("abs(u(t)**2) and abs(v(t)**2)")

plt.show()
