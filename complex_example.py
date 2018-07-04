#!/usr/bin/python

from numpy import pi
import matplotlib.pyplot as plt
from scipy.integrate import ode

def u1(t, y, arg1, arg2):
    return (arg1/2.0j)*arg2*y

def u2(t, y, arg1, arg2):
    u, v = y
    return [v, arg1*v - arg2*u]

def main():
    
    #Parameters
    delta = 0.5
    d_delta = 0.0
    tau = 1e9
    a = float(input("Enter the numerator used for Omega0:  "))

    Omega0 = a/2*pi

    t_vals, u2_vals, u2square_vals, v2_vals, v2square_vals = u2_solver(delta,d_delta,tau,Omega0)

    u1_vals = []
    u1square_vals = []
    v1_vals = []
    v1square_vals = []

    iterator = 0

    for i,j in zip(t_vals, u2_vals):
        u1_result, u1square_result = u1_solver(i,j,Omega0,iterator)

        u1_vals.append(u1_result)
        u1square_vals.append(u1square_result)

        iterator += 1

    print("u1(t) calculated successfully")
    
    plt.figure()
    plt.plot(t_vals, u2_vals, label=r"$u_2(t)$")
    plt.plot(t_vals, v2_vals, label=r"$u_2'(t)$")
    plt.xlabel("Time")
    plt.legend()
    plt.title(r"$u_2(t)$ and $u_2'(t)$:$ \Omega = \frac{%5.2f}{2 \pi}$" % (a))

    plt.figure()
    plt.plot(t_vals, u2square_vals)
    plt.xlabel("Time")
    plt.ylabel(r"$u^{2}_{2}(t)$: Population of excited state")
    plt.title(r"$|u_{2}(t)|^2$: $\Omega = \frac{%5.2f}{2 \pi}$" % (a))

    plt.figure()
    plt.plot(t_vals, u1_vals)

    plt.show()

def u2_solver(delta,d_delta,tau,Omega0):

    f_arg = complex(-1/tau, delta)
    jac_arg = complex(abs(Omega0**2)/4, d_delta - (delta/tau))

    #f_arg = (1.0j * delta) - 1/tau
    #jac_arg = ((abs(Omega0**2))/4) + 1.0j*(d_delta - (delta/tau))

    r = ode(u2).set_integrator('zvode')
    r.set_initial_value([0.0, 1.0],0.0).set_f_params(f_arg,jac_arg)

    tl = 20
    dt = 0.01

    t_vals = []
    u2_vals = []
    u2square_vals = []
    v2_vals = []
    v2square_vals = []

    while r.successful() and r.t < tl:
        r.integrate(r.t+dt)
        t_vals.append(r.t)

        u2_vals.append(r.y[1])
        u2square_vals.append((abs(r.y[1]))**2)
        
        v2_vals.append(r.y[0])
        v2square_vals.append((abs(r.y[0]))**2)

    print("u2(t) calculated successfully")

    return t_vals, u2_vals, u2square_vals, v2_vals, v2square_vals

def u1_solver(time,u2,Omega0,iterator):
    q = ode(u1).set_integrator('zvode')

    q.set_initial_value(time,0.0).set_f_params(Omega0,u2)
    q.integrate(q.t+0.01)

    if q.successful():

        u1_result = q.y
        u1square_result = (abs(q.y))**2 
        return u1_result, u1square_result


if __name__ == "__main__":
    main()
