import numpy as np
import numpy.linalg as lina
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import random
SHOW=1
g = -10000000 #downward acceleration
k = 1.38064852*(10**-23) # (m^2 kg)/(s^2 K)
limit = 1


N =100       # number of particles
dt = .0000075    # timestep
RUNTIME = 1 # seconds
num_iters = int(RUNTIME/dt)


radius = (.01)*(2**(1/6))
mass = 1#.039948/(6.022*10**23)

def Dot(a,b):
    return (a[0]*b[0])+(a[1]*b[1])
def Mag(a):
    return math.sqrt((a[0]**2)+(a[1]**2))

def sumMomentum(x): #x is particle Array
    totalPvec = [0,0]
    for i in range(len(x)):
        totalPvec = [totalPvec[0] + (x[i].m*x[i].v[0]),totalPvec[1] + (x[i].m*x[i].v[1])]
    return(Mag(totalPvec))

def particleE(x): # x is a particle
    return 0.5*x.m*(Mag(x.v)**2)

def sumEnergy(x): #x is particle Array
    totalE = 0
    for i in range(len(x)):
        totalE = totalE + particleE(x[i])
    return totalE

class particle:
    def __init__(self,m,r,u,v,x =limit*random.random(),y = limit*random.random()):
        self.m = m
        self.r = r
        self.vi = [u,v]
        self.v = [self.vi[0],self.vi[1]]
        self.a = [0,g]
        self.F = [0,0]
        self.xi = [x,y] #initial starting posistion
        self.x = [self.xi[0],self.xi[1]]

    def hitwall(self,x):
        if (self.x[1]+self.r) >= x.top and self.v[1] >0: # time step must be small
            self.v[1] = -self.v[1] #elastic cols w/ wall
        if (self.x[1]-self.r) <= x.bot and self.v[1] <0:
            self.v[1] = -self.v[1]
        if (self.x[0]-self.r) <= x.lef and self.v[0] <0:
            self.v[0] = -self.v[0]
        if (self.x[0]+self.r) >= x.rig and self.v[0] >0:
            self.v[0] = -self.v[0]
    def pacmanwall(self,x):
        if (self.x[1]) >= x.top and self.v[1] >0: # time step must be small
            self.x[1] = self.x[1]-1 #elastic cols w/ wall
        if (self.x[1]) <= x.bot and self.v[1] <0:
            self.x[1] = 1+self.x[1]
        if (self.x[0]) <= x.lef and self.v[0] <0:
            self.x[0] = 1+self.x[0]
        if (self.x[0]) >= x.rig and self.v[0] >0:
            self.x[0] = self.x[0]-1
    def specificStart(self,x):
        vel = [3000,50,1000]
        if (self.x[1]) >= x.top and self.v[1] >0: # time step must be small
            self.x[0] = 0
            self.x[1] = .5*limit+(.1*random.random()-.1*random.random())
            self.v[0] = vel[0]+vel[2]*random.random()
            self.v[1] = vel[1]*random.random()-vel[1]*random.random()
        if (self.x[1]) <= x.bot and self.v[1] <0:
            self.x[0] = 0
            self.x[1] = .5*limit
            self.v[0] = vel[0]+vel[2]*random.random()
            self.v[1] = vel[1]*random.random()-vel[1]*random.random()
        if (self.x[0]) <= x.lef and self.v[0] <0:
            self.x[0] = 0
            self.x[1] = .5*limit
            self.v[0] = vel[0]+vel[2]*random.random()
            self.v[1] = vel[1]*random.random()-vel[1]*random.random()
        if (self.x[0]) >= x.rig and self.v[0] >0:
            # self.x[0] = 0
            # self.x[1] = .5*limit
            # self.v[0] = vel[0]
            # self.v[1] = vel[1]*random.random()
            self.v[0] = -self.v[0] #wall
    def waveStart(self,x):
        vel = [3000,50,1000]
        if (self.x[1]) >= x.top and self.v[1] >0: # time step must be small
            self.x[0] = 0
            self.x[1] = self.xi[1]
            self.v[0] = self.vi[0]
            self.v[1] = self.vi[1]+vel[1]*random.random()-vel[1]*random.random()
        if (self.x[1]) <= x.bot and self.v[1] <0:
            self.x[0] = 0
            self.x[1] = self.xi[1]
            self.v[0] = self.vi[0]
            self.v[1] = self.vi[1]+vel[1]*random.random()-vel[1]*random.random()
        if (self.x[0]) <= x.lef and self.v[0] <0:
            self.x[0] = 0
            self.x[1] = self.xi[1]
            self.v[0] = self.vi[0]
            self.v[1] = self.vi[1]+vel[1]*random.random()-vel[1]*random.random()
        if (self.x[0]) >= x.rig and self.v[0] >0:
            self.x[0] = 0
            self.x[1] = self.xi[1]
            self.v[0] = self.vi[0]
            self.v[1] = self.vi[1]+vel[1]*random.random()-vel[1]*random.random()
            # self.v[0] = -self.v[0] #wall

    def accelerate(self):
        self.v[0] = self.v[0] + (self.a[0]*dt)
        self.v[1] = self.v[1] + (self.a[1]*dt)

    def lennardJones(self,distance,sigma,epsilon): # x is particles, return force
        return (24*epsilon*(sigma**6) *(distance**6 -(2*sigma**6)))/(distance**13) # sigma and epsilon need refinement

    def vCollision(self,x,dx):
        V = [0.5*(x.v[0]+self.v[0]),0.5*(x.v[1]+self.v[1])]# COM velocity
        g = [0.5*(self.v[0]-x.v[0]),0.5*(self.v[1]-x.v[1])]# relative velocity
        gprime = [g[0] - ((2*Dot(g,dx)*dx[0])/(Mag(dx)**2)),g[1] - ((2*Dot(g,dx)*dx[1])/(Mag(dx)**2))]
        self.v[0] = V[0] + gprime[0]
        self.v[1] = V[1] + gprime[1]
        x.v[0] = V[0] - gprime[0]
        x.v[1] = V[1] - gprime[1]

    def hitparticle(self,x,pNum): #x is particle array #maybe do range input
        for i in range(len(x)):
            if i != pNum:
                dx = [x[i].x[0] - self.x[0] , x[i].x[1] - self.x[1]]
                dv = [x[i].v[0] - self.v[0] , x[i].v[1] - self.v[1]]
                distance = math.sqrt((dx[0]**2)+(dx[1]**2))
                vDist = ((dx[0]*dv[0])+(dx[1]*dv[1]))# / distance
                if distance<(self.r + x[i].r):
    				## go if approaching each other
                    if vDist<0.:
    					## update
                        self.vCollision(x[i],dx)

    def triangleCheck(self,tri): #x is particle
        if self.x[1]-math.sqrt(((self.r**2)+(self.r*math.tan((math.pi/2)-tri.theta))**2)) < ((tri.side1[0]*self.x[0])+tri.side1[1]):
            if Dot(self.v,tri.unitn)*tri.unitn[1]<0:#/Dot(self.v,tri.unitn)*tri.unitn[1]
                self.v = self.triangleVelocity(tri)
                # print(self.x)
    def triangleVelocity(self,tri):
        vpara = [Dot(self.v,tri.unitx)*tri.unitx[0],Dot(self.v,tri.unitx)*tri.unitx[1]]
        vperp = [-Dot(self.v,tri.unitn)*tri.unitn[0],-Dot(self.v,tri.unitn)*tri.unitn[1]]
        # vpara = (self.v[0]/math.cos(tri.theta)) + (self.v[1]/math.cos(tri.theta))
        # vperp = -((self.v[0]/math.sin(tri.theta)) + (self.v[1]/math.sin(tri.theta)))
        return [vpara[0]+vperp[0],vpara[1]+vperp[1]]


    def timestep(self,box,pArray,pNum,tri):
        self.accelerate()
        self.hitparticle(pArray,pNum)
        self.x[0] = self.x[0] + (self.v[0] * dt)
        self.x[1] = self.x[1] + (self.v[1] * dt)
        # self.triangleCheck(triangle1)
        self.hitwall(box)

        # self.pacmanwall(box)
        # self.waveStart(box)
        # self.specificStart(box)
        # print(sumEnergy(pArray))

###########
class triangle:
    def __init__(self, pt11, pt12, pt21, pt22): # (1-2top)
        self.side1 = [(pt12-pt22)/(pt11-pt21), -((pt12-pt22)/(pt11-pt21) *(pt11))+pt12 ]   #[m,b] (y = mx+b)

        self.unitx = [pt21-pt11,pt22-pt12]
        self.magn = Mag(self.unitx)
        self.unitx = [self.unitx[0]/self.magn,self.unitx[1]/self.magn]
        self.unitn = [-self.unitx[1],self.unitx[0]]
        self.theta = math.atan(self.side1[0])


###########
class box:
    def __init__(self,top,bot,lef,rig):
        self.top = top
        self.bot = bot
        self.lef = lef
        self.rig = rig

#####################
box1 = box(limit,0,0,limit)
triangle1 = triangle(.25*limit,0*limit,1*limit,.3*limit)
# print(triangle1.theta,triangle1.side1)

prtcls = []
for i in range(N):
    prtcls.append(particle(mass,radius*1,limit*1000,limit*500,limit*random.random(),limit*random.random())) # *random.random()*(-1)**random.randint(0,10)

# starting point
# b = 10
# a  = int(N/b)
# for i in range(a):
#     for j in range(b):
#         prtcls.append(particle(mass,radius*1,1000,10*random.random()-10*random.random(),0-(.1*j),.04+int(i % a)/a))


if SHOW:
    fig = plt.figure(1, figsize=(5,5))
    ax = fig.add_subplot(111)
    datax = []
    datay = []
    for i in range(N):
        datax.append(prtcls[i].x[0])
        datay.append(prtcls[i].x[1])
    line1, = ax.plot(datax, datay, 'ko')
    # plt.show()
    for time in range(num_iters):
        for i in range(N):
            prtcls[i].timestep(box1,prtcls,i,triangle1)
            datax[i]= (prtcls[i].x[0])
            datay[i]= (prtcls[i].x[1])
        line1.set_xdata(datax)
        line1.set_ydata(datay)
        # print(sumEnergy(prtcls))
        # print(sumMomentum(prtcls))
        # print()
        plt.axis([0, limit, 0, limit], aspect=1.)
        fig.canvas.draw()
        plt.pause(0.000000001)
        fig.canvas.flush_events()
####ax = plt.axes([.1,.1,.8,.8], aspect=1.)plt.figure(1, figsize=(5,5))

if abs(SHOW-1):
    plt.ion()
    Vdistr = []
    fig = plt.figure()
    num_bins = 400
    for time in range(num_iters):
        if time < 100:
            Vdistr = []
        for i in range(N):
            prtcls[i].timestep(box1,prtcls,i,triangle1)
            Vdistr.append(prtcls[i].v[0])
        plt.clf()
        n, bins, patches = plt.hist(Vdistr, num_bins, facecolor='blue', alpha=1)
        plt.show()
        plt.pause(.000000001)
