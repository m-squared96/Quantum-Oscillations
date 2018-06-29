from scipy.integrate import complex_ode
import numpy as np
def model(t,y,params):
  
  omega_0 = params[0]
  tau = params[1]

  x = y[0]
  dx = y[1]  
  K = (omega_0**2)/4 - ((i*detuning)/tau)
  xdot = [[],[]]
  xdot[0] = dx
  xdot[1] = (-i*detuning + (1/tau)) * dx - K * x

  return xdot

omega_0 = 1/(2*np.pi)
tau = 1000
detuning = 0.01
i = 1j

params = [omega_0, tau]

#t = np.arange(0,200,0.05)
z = complex_ode(model)
z.set_initial_value(np.array([100,1,1], dtype=np.complex128),0)
print(z)

t1 = 10
dt = 1

sol = np.array([], dtype = np.complex128)
t   = np.array([], dtype = np.complex128)

while z.successful() and z.t < t1:
    t=np.append(t,z.t+dt)
    sol=np.append(sol, z.integrate(z.t+dt))

# plot results
A1 = sol[:, 0]
A2 = sol[:, 1]
A3 = sol[:, 2]

plt.figure()
plt.plot(t, abs(A1), label='A1')
plt.plot(t, abs(A2), label='A2')
plt.plot(t, abs(A3), label='A3')
plt.xlabel('t')
plt.ylabel('A')
plt.title('A-t')
plt.legend(loc=0)