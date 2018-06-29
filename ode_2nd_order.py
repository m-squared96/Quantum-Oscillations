from scipy.integrate import odeint
import numpy as np
def model(y,t,params):
  
  omega_0 = params[0]
  tau = params[1]

  x = y[0]
  dx = y[1]  
  K = (omega_0**2)/4
  xdot = [[],[]]
  xdot[0] = dx
  xdot[1] = (1/tau) * dx - K * x

  return xdot

omega_0 = 1/(2*np.pi)
tau = 1000

params = [omega_0, tau]

t = np.arange(0,200,0.05)
z = odeint(model,[1.0,7.0],t, args=(params,))
print(z)


# plot results
import matplotlib.pyplot as plt
plt.plot(t,z[:,0],'g:')
plt.plot(t,z[:,1],'k-.')
plt.xlabel('Time')
plt.show()
