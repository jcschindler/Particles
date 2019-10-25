import numpy as np
import numpy.linalg as lina
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.animation as an
import random
SHOW=1
g = -000000 #downward acceleration
k = 1.38064852*(10**-23) # (m^2 kg)/(s^2 K)
limit = 1

rows = 12
N =rows*200       # number of particles
dt = .000005    # timestep
RUNTIME = .004 # seconds
num_iters = int(RUNTIME/dt)
Tr = 5 #number of dt's

velocity = 1500
spread = 150
counter = 0
sty = dict(markersize=2.5)

interval = 20
fps = 30.

Argon = False
if Argon == True:
    epsilon = 1.65*(10**-21)
    radius = (3.4*(10**-10))*(2**(1/6)) # sigma *(2**(1/6))
    mass = .039948/(6.022*10**23) #(kg/mol)/(N_a)
else: # random
    epsilon = 1.65*(10**-21)
    radius = (.01)*(2**(1/6))
    radius = .007
    mass = .039948/(6.022*10**23)

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
    def __init__(self,m,r,u,v,x =limit*random.random(),y = limit*random.random(),active = False):
        self.m = m
        self.r = r
        self.vi = [u,v]
        self.v = [self.vi[0],self.vi[1]]
        self.a = [0,g]
        self.F = [0,0]
        self.xi = [x,y] #initial starting posistion
        self.x = [self.xi[0],self.xi[1]]
        self.active = active

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
    def waveBounds(self,x):
        if (self.x[1]) > x.top and self.v[1] >0: # time step must be small
            self.active = False
            self.x = [10,10]
            self.v = [0,0]
        if (self.x[1]) < x.bot and self.v[1] <0:
            self.active = False
            self.x = [10,10]
            self.v = [0,0]
        if (self.x[0]) < x.lef and self.v[0] <0:
            self.active = False
            self.x = [10,10]
            self.v = [0,0]
        if (self.x[0]) > x.rig and self.v[0] >0:
            self.active = False
            self.x = [10,10]
            self.v = [0,0]
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

    def timestep(self,box,pArray,pNum,tri,a,counter):
        if self.active == True:
            # self.accelerate()
            self.hitparticle(pArray,pNum)
            self.x[0] = self.x[0] + (self.v[0] * dt)
            self.x[1] = self.x[1] + (self.v[1] * dt)
            self.triangleCheck(tri)
            # self.hitwall(box)
            # self.pacmanwall(box)
            self.waveBounds(box)
            # self.specificStart(box)
            # print(sumMomentum(pArray))
def which_on(passiveArray):
    # x = [passiveArray.index(i) for i in passiveArray if i == True]
    return [i for i, x in enumerate(passiveArray) if x]



def addWave(pArray,a,counter):
    if counter/1. >= Tr:
        indices = which_on(queue(pArray))
        # print(indices)
        for i in range(a):
            pArray[indices[i]].x = [0,.04+int(i % a)/a]
            pArray[indices[i]].v = [pArray[indices[i]].vi[0],pArray[indices[i]].vi[1]]#pArray[indices[i]].vi
            pArray[indices[i]].active = True
        return 1 #reset count
    else:
        return 0



def queue(pArray):
    passiveArray = []
    for i in range(N):
        passiveArray.append(not pArray[i].active)
    return(passiveArray)
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
triangle1 = triangle(.25*limit,0*limit,1*limit,.45*limit)
# print(triangle1.side1)

prtcls = []
# nums = []
# for i in range(N):
#     prtcls.append(particle(mass,radius*1,limit*1000*random.random()*(-1)**random.randint(0,10),limit*500*random.random()*(-1)**random.randint(0,10),limit*random.random(),limit*random.random()))

## starting point
b = int(N/rows)
a  = int(N/b)
for i in range(a):
    for j in range(b):
        if j == 0:
            acti = True
        else:
            acti = False
        prtcls.append(particle(m = mass,r = radius*1,u = velocity+spread*random.random()-spread*random.random(),v = spread*random.random()-spread*random.random(),x = 0 - (100*j),y = .04+int(i % a)/a, active = acti))
        # nums.append(j+(i*b))

def gendata(box,prtcls,tri,a,counter,num_iters,N):
    datax = np.zeros((N,2,num_iters))
    for i in range(N):
        datax[i,:,0]=[prtcls[i].x[0],prtcls[i].x[1]]
    # plt.show()
    for time in range(num_iters):
        for i in range(N):
            prtcls[i].timestep(box1,prtcls,i,triangle1,a,counter)
            datax[i,:,time]=[prtcls[i].x[0],prtcls[i].x[1]]
        reset =(prtcls,a,counter)
        print(time)
        if reset == True:
            counter = 0
        counter = counter +1
    return datax

def run(i, data):
    ## status
    print(i)
    ## total energy
    #print("P = %s"%(np.round(data[5][:,i],2)))
    #print("E = %6.2f"%(data[4][i]))
    ## update
    # x = 1.*data[1]
    x = data
    ## plot
    points.set_data(x[:,0,i], x[:,1,i])
    ## return
    return [points,]
## plot
def init():
	plt.xticks([0,1])
	plt.yticks([0,1])
	plt.xlim(0,1)
	plt.ylim(0,1)
	#plt.tight_layout()
	return [points,]

# plt.plot([triangle1.side1[1],triangle1.side1[0]+triangle1.side1[1]],'k')


########
fig = plt.figure(1, figsize=(5,5))
ax = plt.axes([.1,.1,.8,.8], aspect=1.)
global points
points, = plt.plot([],[],'ko', **sty)
## make data
print("gendata")
data = gendata(box1,prtcls,triangle1,a,counter,num_iters,N)
## animate
print("animate")
ani = an.FuncAnimation(fig, run, fargs=(data,), interval=interval, init_func=init, frames=range(num_iters))
## save
print("save")
if False:
    ani.save("temp.gif", writer="imagemagick", dpi=100, fps=fps )
## show
if True:
    plt.show()
