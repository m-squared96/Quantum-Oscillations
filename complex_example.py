#!/usr/bin/python

from numpy import pi,sqrt,cos,sin
import matplotlib.pyplot as plt
from scipy.integrate import ode

def u1squared(t, Omega_tilda, delta):
    return (cos((Omega_tilda*t)/2))**2 + ((delta/Omega_tilda)*sin(Omega_tilda*t/2))**2

def u2(t, y, arg1, arg2):
    u, v = y
    return [v, arg1*v - arg2*u]

def u2_solver(delta,d_delta,tau,Omega0):

    f_arg = complex(-1/tau, delta)
    jac_arg = complex(abs(Omega0**2)/4, d_delta - (delta/tau))

    #f_arg = (1.0j * delta) - 1/tau
    #jac_arg = ((abs(Omega0**2))/4) + 1.0j*(d_delta - (delta/tau))

    r = ode(u2).set_integrator('zvode')
    r.set_initial_value([0.0, 1.0],0.0).set_f_params(f_arg,jac_arg)

    tl = 20
    dt = 0.001

    t_vals = []
    u2_vals = []
    u2square_vals = []
    v2_vals = []
    v2square_vals = []

    while r.successful() and r.t < tl:
        r.integrate(r.t+dt)
        t_vals.append(r.t)

        u2_vals.append(r.y[0])
        u2square_vals.append((abs(r.y[0]))**2)

        v2_vals.append(r.y[1])
        v2square_vals.append((abs(r.y[1]))**2)

    return t_vals, u2_vals, u2square_vals, v2_vals, v2square_vals

def u1_solver(time,Omega_tilda,delta):

    u1square_vals = []

    for t in time:
        u1square_vals.append(u1squared(t,Omega_tilda,delta))

    return u1square_vals


def main():

    #Parameters
    delta = 0.5
    d_delta = 0.0
    tau = 1e9
    #a = float(input("Enter the numerator used for Omega0:  "))
    a = 1.274

    Omega0 = a/2*pi
    Omega_tilda = sqrt(delta**2 + abs(Omega0**2))

    t_vals, u2_vals, u2square_vals, v2_vals, v2square_vals = u2_solver(delta,d_delta,tau,Omega0)
    u1square_vals = u1_solver(t_vals,Omega_tilda,delta)

    diff = []
    for i,j in zip(u1square_vals,u2square_vals):
        diff.append(i+j)

    plt.figure()
    plt.plot(t_vals, u2_vals, label=r"$u_2(t)$")
    plt.plot(t_vals, v2_vals, label=r"$u_2'(t)$")
    plt.xlabel("Time")
    plt.legend()
    plt.title(r"$u_2(t)$ and $u_2'(t)$:$ \Omega_0 = \frac{%5.2f}{2 \pi}$" % (a))

    plt.figure()
    plt.plot(t_vals, u1square_vals, label=r"$|u_1(t)|^2$")
    plt.plot(t_vals, u2square_vals, label=r"$|u_2(t)|^2$")
    plt.plot(t_vals, diff, label=r"$|u_1(t)|^{2} + |u_2(t)|^{2}$", ls="--")
    plt.xlabel("Time")
    plt.ylabel("Population of states")
    plt.legend()
    plt.title(r"$|u_{2}(t)|^2$ and $|u_{1}(t)|^2$ : $\Omega_0 = \frac{%5.2f}{2 \pi}$" % (a))

    plt.show()

main()
