import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.animation as an



V_spread = 100
dt = .0000075    # timestep
tNumSteps = 1000
nParticles = 100

radius = .02#(.01)*(2**(1/6))

limit = 1


x = np.zeros([tNumSteps,nParticles,2])
r = np.ones([nParticles])*radius
v = np.zeros([nParticles,2])


for i in range(nParticles):
    x[0,i,:] = [limit*np.random.random(),limit*np.random.random()]
    v[i,:] = [V_spread*np.random.randn(),V_spread*np.random.randn()]

# print(v)
message = 
