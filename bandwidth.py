#!/usr/bin/python

from numpy import pi,sqrt,cos,sin,exp,conj,real,linspace,imag
import matplotlib.pyplot as plt
from scipy.integrate import ode

def ug(s,gamma,t):
    return exp(-1.0j*t*(s - 1.0j*gamma)/2)

def sg(t):
    return sin(t/2*pi) + cos(t/2*pi)

def main():
    time = linspace(0.0,10.0,1000)

    gamma = 0.5

    sg_vals = []
    ug_vals = []

    ug_reals = []
    ug_imag = []
    ug_squares = []

    for t in time:
        sg_vals.append(sg(t))

    for t,s in zip(time, sg_vals):
        ug_vals.append(ug(s,gamma,t))

    for u in ug_vals:
        ug_reals.append(real(u))
        ug_imag.append(imag(u))
        ug_squares.append(real(u * u.conj()))

    plt.figure()
    plt.plot(time, ug_reals, label=r"Real $u_g(t)$")
    plt.plot(time, ug_imag, label=r"Complex $u_g(t)$")
    plt.xlabel("Time")
    plt.legend()

    plt.figure()
    plt.plot(time, ug_squares)
    plt.xlabel("Time")
    plt.ylabel(r"$|u_g(t)^2|$")

    plt.show()

main()
