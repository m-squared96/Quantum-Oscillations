import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def wrapper(t, y):
	return f(t, y, params)           # solve_ivp doesn't have an args parameter so we use a wrapper

def f(t, y, params):
    u1, u2 = y                       # unpack current values of y
    d, Omega = params                # unpack parameters
    derivs = [(Omega/(2*1j)) * u2,   # list of dy/dt=f functions
             ((d * u2)/1j) + (np.conj(Omega)/(2*1j))*u1]
    return derivs

def pos(x):
	return (abs(x))**2

# Parameters
Omega = 1/(2*np.pi)              # Rabi parameter
d = 0.25                         # detuning

# Initial values
u1_0 = complex(1.0)              # u1 value at t = 0
u2_0 = complex(0.0)              # u2 value at t = 0

# Bundle parameters for ODE solver
params = [d, Omega]

# Bundle initial conditions for ODE solver
y0 = [u1_0, u2_0]

# Time interval for solution
t = [0, 100]

# Call the ODE solver
res = solve_ivp(wrapper, t, y0, t_eval=np.linspace(0, 100, 20000))


plt.figure()
plt.plot(res.t, pos(res.y[0]), label=r"$|u_1(t)|^2$")
plt.plot(res.t, pos(res.y[1]), label=r"$|u_2(t)|^2$")
plt.plot(res.t, sum([pos(res.y[0]), pos(res.y[1])]), label=r"$|u_1(t)|^{2} + |u_2(t)|^{2}$", ls="--")
plt.xlabel("Time")
plt.ylabel("Population of states")
plt.legend()
plt.title(r"$|u_{1}(t)|^2$ and $|u_{2}(t)|^2$ : $\Omega_0 = \frac{%5.2f}{2 \pi}$" % (Omega*2*np.pi))

plt.show()