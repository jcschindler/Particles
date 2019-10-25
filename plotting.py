import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.animation as an



sty = dict(markersize=2.5)

interval = 20
fps = 30.

xdata = np.loadtxt('/Users/Jeff/Desktop/python/Helpful Code/Particle Box/Fortran/xData.txt')
ydata = np.loadtxt('/Users/Jeff/Desktop/python/Helpful Code/Particle Box/Fortran/yData.txt')

num_iters = np.shape(xdata)[0]

def run(i, xdata,ydata):
    ## status
    print(i)
    ## total energy
    #print("P = %s"%(np.round(data[5][:,i],2)))
    #print("E = %6.2f"%(data[4][i]))
    ## update
    # x = 1.*data[1]
    # x = data
    ## plot
    points.set_data(xdata[i,:], ydata[i,:])
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


fig = plt.figure(1, figsize=(5,5))
ax = plt.axes([.1,.1,.8,.8], aspect=1.)
global points
points, = plt.plot([],[],'ko', **sty)
## make data
# print("gendata")
# data = gendata(box1,prtcls,triangle1,a,counter,num_iters,N)
## animate
print("animate")
ani = an.FuncAnimation(fig, run, fargs=(xdata,ydata,), interval=interval, init_func=init, frames=range(num_iters))
## save
print("save")
if False:
    ani.save("temp.gif", writer="imagemagick", dpi=100, fps=fps )
## show
if True:
    plt.show()
