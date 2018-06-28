from scipy.integrate import odeint
import numpy as np
def model(x,t):
  y = x[0]
  dy = x[1]  
  K = (omega_0**2)/4
  xdot = [[],[]]
  xdot[0] = dy
  xdot[1] = (1/tau) * dy - K * y

  return xdot

omega_0 = 1/(2*np.pi)
tau = 1000

time = np.linspace(0,200,0.05)
z = odeint(model,[2,-1],time)

# plot results
import matplotlib.pyplot as plt
plt.plot(time,z[:,0],'g:')
plt.plot(time,z[:,1],'k-.')
plt.xlabel('Time')
plt.show()
