import numpy as np
import matplotlib
import math
matplotlib.use('TkAgg')
import pylab
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

class Square():
    def __init__(self,poscm,radius,vcm,phi,k,d,eta,kt,dt):
        self.dt=0.005
        self.pcm=np.array([0,0])
        self.vcm=vcm
        self.acm=np.array([0,0])
        self.k=k
        self.d=d
        self.eta=eta
        self.kt=kt
        self.r=0.1
        self.m=0.1
        self.M=self.m*100
        self.a=0
        self.theta=0
        self.thetadot=0
        self.thetadotdot=0
        self.box=np.ndarray((10/d+10,10/d+10))
        self.box=np.zeros_like(self.box)
        self.pa=np.ndarray((100,2))
        self.pa=np.zeros_like(self.pa)
        self.boxx=np.zeros_like(np.ndarray((10/d+10,10/d+10)))
        self.boxy=np.zeros_like(np.ndarray((10/d+10,10/d+10)))


    def createsquare(self):
        self.pa=np.ndarray((100,2))
        startx=7
        starty=5
        i=0
        for y in range(0,10):
            for x in range(0,10):
                self.pa[i,0]=startx+0.2*x
                self.pa[i,1]=starty+0.2*y
                i+=1
                if i==100:
                    break
        for i in range(0,100):
            self.pcm=self.pcm+self.m*self.pa[i]/self.M
        self.pacm=np.ndarray((100,2))
        self.pav=np.ndarray((100,2))
        self.pavcm=np.zeros_like(np.ndarray((100,2)))
        for l in range(0,100):
            self.pacm[l]=self.pa[l]-self.pcm
            self.pav[l]=self.vcm

        

    def createboundary(self):
        self.ba=np.ndarray((1000,2))
        self.ba=np.zeros_like(self.ba)
        for j in range(0,5):
            for i in range(0,51):
                self.ba[i+j*51,0]=self.d*i
                self.ba[i+j*51,1]=self.d*j

        for k in range(0,5):
            for l in range(5,51):
                self.ba[255+l+46*k,0]=self.d*k
                self.ba[255+l+46*k,1]=self.d*l

        for m in range(0,5):
            for n in range(5,51):
                self.ba[485+n+m*46,0]=10-m*self.d
                self.ba[485+n+m*46,1]=self.d*n

        for o in range(0,5):
            for p in range(5,46):
                self.ba[715+p+o*41,0]=self.d*p
                self.ba[715+p+o*41,1]=10-self.d*o
                                       
        
        

    def neighbour_box_co(self,arrayij): # Takes coordinates of a box and returns its neighbouring boxes
        neighbours=[[arrayij[0],arrayij[1]],[arrayij[0]-1,arrayij[1]+1],[arrayij[0],arrayij[1]+1],[arrayij[0]+1,arrayij[1]+1],[arrayij[0]+1,arrayij[1]],[arrayij[0]-1,arrayij[1]],[arrayij[0]-1,arrayij[1]-1],[arrayij[0],arrayij[1]-1],[arrayij[0]+1,arrayij[1]-1]]
        return neighbours


    def boxlist_update(self): #assigns box number to each particle at the end of each time step
        d=self.d
        self.boxx,self.boxy=np.zeros_like(np.ndarray((10/d+10,10/d+10))),np.zeros_like(np.ndarray((10/d+10,10/d+10)))
        for i in range(0,1000):
            x,y=self.box_co(self.ba[i])
            self.boxx[x,y]=self.ba[i,0]
            self.boxy[x,y]=self.ba[i,1]
            

    def box_co(self,array):# takes particle's location and returns box's coordinates i and j
        x=int(array[0]/self.d)
        y=int(array[1]/self.d)
        return(x,y)


    def move(self):
        dt=self.dt
        self.theta=self.thetadot*dt+self.theta
        self.pcm=self.pcm+self.vcm*dt
        for i in range(0,100):
            self.pavcm[i,0],self.pavcm[i,1]=np.cross([0,0,self.thetadot],[self.pacm[i,0],self.pacm[i,1],0])[0],np.cross([0,0,self.thetadot],[self.pacm[i,0],self.pacm[i,1],0])[1]
            self.pav[i]=self.vcm+self.pavcm[i]
        self.pacm=self.pavcm*dt+self.pacm
        self.pa=self.pa+self.pav*dt


    def check_collision(self):
        dt=self.dt
        neighbours=[]
        pairs=[]
        for i in range(0,100):
            neighbours=self.neighbour_box_co(self.box_co(self.pa[i]))
            for k in range(0,len(neighbours)):
                if self.boxx[neighbours[k][0],neighbours[k][1]]!=0: #Checks if the neighbouring box has a particle
                    x,y=self.boxx[neighbours[k][0],neighbours[k][1]],self.boxy[neighbours[k][0],neighbours[k][1]]
                    if np.sqrt(np.dot(self.pa[i]-np.array([x,y]),self.pa[i]-np.array([x,y])))<self.d: #checks for collision
                        pairs.append([i,x,y])
        self.collided(pairs)
        return

    def fs(self,rij):#takes the difference in dispacements of the 2 colliding particles
        modrij=np.sqrt(np.dot(rij,rij))
        if modrij==0:
            return(np.array([0,0]))
        return(-self.k*(self.d-modrij)*rij/modrij)

    def fd(self,vij):# takes the difference in vel of 2 colliding particles
        return(self.eta*vij)

    def ft(self,vij,rij):# takes the difference in vel and displacement of 2 colliding particles
        modrij=np.sqrt(np.dot(rij,rij))
        if modrij==0:
            return(np.array([0,0]))
        vijt=vij-(np.dot(vij,rij/modrij))*rij/modrij
        return(self.kt*vijt)

    def MOI(self):
        self.I=np.ndarray((3,3))
        self.I=np.zeros_like(self.I)
        for i in range(0,100):
            self.I[0,0]=self.I[0,0]+self.m*self.pacm[i,1]**2
            self.I[1,1]=self.I[1,1]+self.m*self.pacm[i,0]**2
            self.I[2,2]=self.I[2,2]+self.m*self.pacm[i,0]**2+self.m*self.pacm[i,1]**2
            self.I[0,1]=self.I[1,2]-self.m*self.pacm[i,0]*self.pacm[i,1]
        self.I[1,0]=self.I[0,1]
        return(self.I)
        
        

    def collided(self,pairs):
        dt=self.dt
        F=np.array([0,0]) #total force
        T=0 #total torque
        for i in range(0,len(pairs)):
            rij=np.array(pairs[i][1]-self.pa[pairs[i][0],0],pairs[i][2]-self.pa[pairs[i][0],1])
            vij=-self.pav[pairs[i][0]]
            F=F+self.fs(rij)+self.fd(vij)+self.ft(vij,rij)
            T=T+np.cross(self.pacm[pairs[i][0]],F)
        self.thetadotdot=np.dot(np.linalg.inv(self.MOI()),np.array([0,0,T]))[2]
        self.thetadot=self.thetadot+self.thetadotdot*dt
        self.acm=F/self.M
        self.vcm=self.vcm+self.acm*dt
        return
            
            

            
            
def main():
    dt=0.001
    s=Square(np.array([5.9,5.9]),0.1,np.array([1,2]),1,200,0.2,10,1,0.005)
    s.createsquare()
    s.createboundary()
    fig=pylab.figure(1)
    ax=fig.add_subplot(111,xlim=(0,10),ylim=(0,10))
    scatter,=plt.plot([],[],'o',ms=3)
    bound,=plt.plot([],[],'o',ms=3)

    def init():
        scatter.set_data([],[])
        bound.set_data([s.ba[:,0]],[s.ba[:,1]])
        return(scatter,bound,)

    def animate(i):
        global dt
        scatter.set_data(s.pa[:,0],s.pa[:,1])
        s.move()
        s.check_collision()
        s.boxlist_update()
        return(scatter,)

    anim=FuncAnimation(fig,animate,frames=200,init_func=init,interval=dt*1000,blit=True)

    pylab.show()

if __name__=='__main__':
    main()


        
        
