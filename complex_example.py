#!/usr/bin/python

import matplotlib.pyplot as plt
from scipy.integrate import ode

def f(t, y, arg1, arg2):
    u, v = y
    return [v, arg1*v -arg2*u]

#Parameters
delta = 0.0
d_delta = 0.0
tau = 1000
Omega0 = 1.0j

f_arg = (1.0j * delta) - 1/tau
jac_arg = ((abs(Omega0**2))/4) + 1.0j*(d_delta - (delta/tau))

r = ode(f).set_integrator('zvode')
r.set_initial_value([1.0, 0.0],1.0).set_f_params(f_arg,jac_arg)

tl = 200
dt = 0.01

t_vals = []
u_vals = []
usquare_vals = []

while r.successful() and r.t < tl:
    r.integrate(r.t+dt)
    t_vals.append(r.t)
    u_vals.append(r.y)
    usquare_vals.append((r.y)**2)

plt.figure()
plt.plot(t_vals, u_vals)
plt.xlabel("Time")
plt.title("u(t) and v(t)")

plt.figure()
plt.plot(t_vals, usquare_vals)
plt.xlabel("Time")
plt.title("abs(u(t)**2) and abs(v(t)**2)")

plt.show()
