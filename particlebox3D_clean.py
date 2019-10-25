import numpy as np
import numpy.linalg as lina
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
SHOW=1
g = -0 #downward acc
k = 1.38064852*(10**-23) # (m^2 kg)/(s^2 K)
limit = 1

N = 50      # number of particles
dt = .00001    # timestep
RUNTIME = .01 # seconds
num_iters = int(RUNTIME/dt)

Argon = False
if Argon == True:
    epsilon = 1.65*(10**-21)
    radius = (3.4*(10**-10))*(2**(1/6)) # sigma *(2**(1/6))
    mass = .039948/(6.022*10**23) #(kg/mol)/(N_a)
else: # random
    epsilon = 1.65*(10**-21)
    radius = (.01)*(2**(1/6))
    mass = .039948/(6.022*10**23)

def Dot(a,b):
    return (a[0]*b[0])+(a[1]*b[1])+(a[2]*b[2])
def Mag(a):
    return math.sqrt((a[0]**2)+(a[1]**2)+(a[2]**2))

def sumMomentum(x): #x is particle Array
    totalPvec = [0,0,0]
    for i in range(len(x)):
        totalPvec = [totalPvec[0] + (x[i].m*x[i].v[0]),totalPvec[1] + (x[i].m*x[i].v[1]),totalPvec[2] + (x[i].m*x[i].v[2])]
    return(Mag(totalPvec))

def particleE(x): # x is a particle
    return 0.5*x.m*(Mag(x.v)**2)

def sumEnergy(x): #x is particle Array
    totalE = 0
    for i in range(len(x)):
        totalE = totalE + particleE(x[i])
    return totalE

class particle:
    def __init__(self,m,r,x,y,z,u,v,w):
        self.m = m
        self.r = r
        self.x = [x,y,z]
        self.v = [u,v,w]
        self.a = [0,0,g]
        self.F = [0,0,0]

    def hitwall(self,x):
        if (self.x[1]+self.r) >= x.top and self.v[1] >0: # time step must be small
            self.v[1] = -self.v[1] #elastic cols w/ wall
        if (self.x[1]-self.r) <= x.bot and self.v[1] <0:
            self.v[1] = -self.v[1]
        if (self.x[0]-self.r) <= x.lef and self.v[0] <0:
            self.v[0] = -self.v[0]
        if (self.x[0]+self.r) >= x.rig and self.v[0] >0:
            self.v[0] = -self.v[0]
        if (self.x[2]-self.r) <= x.bac and self.v[2] <0:
            self.v[2] = -self.v[2]
        if (self.x[2]+self.r) >= x.fro and self.v[2] >0:
            self.v[2] = -self.v[2]

    def accelerate(self):
        self.v[0] = self.v[0] + (self.a[0]*dt)
        self.v[1] = self.v[1] + (self.a[1]*dt)
        self.v[2] = self.v[2] + (self.a[2]*dt)

    def lennardJones(self,distance,sigma,epsilon): # x is particles, return force
        return (24*epsilon*(sigma**6) *(distance**6 -(2*sigma**6)))/(distance**13) # sigma and epsilon need refinement

    def vCollision(self,x,dx):
        V = [0.5*(x.v[0]+self.v[0]),0.5*(x.v[1]+self.v[1]),0.5*(x.v[2]+self.v[2])]# COM velocity
        g = [0.5*(self.v[0]-x.v[0]),0.5*(self.v[1]-x.v[1]),0.5*(self.v[2]-x.v[2])]# relative velocity
        gprime = [g[0] - ((2*Dot(g,dx)*dx[0])/(Mag(dx)**2)),g[1] - ((2*Dot(g,dx)*dx[1])/(Mag(dx)**2)) ,g[2] - ((2*Dot(g,dx)*dx[2])/(Mag(dx)**2))]
        self.v[0] = V[0] + gprime[0]
        self.v[1] = V[1] + gprime[1]
        self.v[2] = V[2] + gprime[2]
        x.v[0] = V[0] - gprime[0]
        x.v[1] = V[1] - gprime[1]
        x.v[2] = V[2] - gprime[2]

    def hitparticle(self,x,pNum): #x is particle array #maybe do range input
        for i in range(len(x)):
            if i != pNum:
                dx = [x[i].x[0] - self.x[0] , x[i].x[1] - self.x[1] , x[i].x[2] - self.x[2]]
                dv = [x[i].v[0] - self.v[0] , x[i].v[1] - self.v[1] , x[i].v[2] - self.v[2]]
                distance = Mag(dx)
                vDist = Dot(dx,dv) / distance
                if distance<(self.r + x[i].r):
    				## go if approaching each other
                    if vDist<0.:
    					## update
                        self.vCollision(x[i],dx)

    def timestep(self,x,pArray,pNum):
        self.accelerate()
        self.hitparticle(pArray,pNum)
        self.x[0] = self.x[0] + (self.v[0] * dt)
        self.x[1] = self.x[1] + (self.v[1] * dt)
        self.x[2] = self.x[2] + (self.v[2] * dt)
        self.hitwall(x)
        # print(sumMomentum(pArray))

class box:
    def __init__(self,top,bottom,left,right,front,back):
        self.top = top
        self.bot = bottom
        self.lef = left
        self.rig = right
        self.fro = front
        self.bac = back

#####################
box1 = box(limit,0,0,limit,limit,0)

prtcls = []
for i in range(N):
    prtcls.append(particle(mass,radius*1,limit*random.random(),limit*random.random(),limit*random.random(),limit*500*random.random()*(-1)**random.randint(0,10),limit*5000*random.random()*(-1)**random.randint(0,10),limit*500*random.random()*(-1)**random.randint(0,10)))

if SHOW:
    plt.ion()
    fig = plt.figure()
    ax = Axes3D(fig)#fig.add_subplot(111, projection='3d')
    datax = []
    datay = []
    dataz = []
    for i in range(N):
        datax.append(prtcls[i].x[0])
        datay.append(prtcls[i].x[1])
        dataz.append(prtcls[i].x[2])
    # line1, = ax.plot(datax, datay, dataz, 'ko')
    ax.scatter(datax, datay, dataz,c='r')
    # plt.show
    for time in range(num_iters):
        for i in range(N):
            prtcls[i].timestep(box1,prtcls,i)
            datax[i]= (prtcls[i].x[0])
            datay[i]= (prtcls[i].x[1])
            dataz[i]= (prtcls[i].x[2])

        ax.scatter(datax, datay, dataz,c='r')
        plt.axis([0, limit, 0, limit])
        ax.set_zlim(0, limit)
        plt.show()
        plt.pause(.00000000000001)
        plt.cla()

if abs(SHOW-1):
    # plt.ion()
    Edistr = []
    Vdistr = []
    fig = plt.figure()
    num_bins = 200
    for time in range(num_iters):
        for i in range(N):
            prtcls[i].timestep(box1,prtcls,i)
            Edistr.append(particleE(prtcls[i]))
            Vdistr.append(prtcls[i].v[0])

        # plt.clf()
        # n, bins, patches = plt.hist(Vdistr, num_bins, facecolor='blue', alpha=1)
        # plt.show()
        # plt.pause(.0000000000001)
    # plt.figure(0)
    n, bins, patches = plt.hist(Vdistr, num_bins, facecolor='blue', alpha=1)
    plt.title('Velocity Distribution')
    plt.show()
    # plt.figure(1)
    n, bins, patches = plt.hist(Edistr, num_bins, facecolor='blue', alpha=1)
    plt.title('Energy Distribution')
    plt.show()
    # plt.figure(1)
