import numpy as np
import math
import Untitled


def Evolution(s):
    s.move()
    s.check_collision()
    s.boxlist_update()

def Energy(s,dt,timearray,t):
    i,j=0,0
    print timearray[-1]
    while j<=len(timearray)-1:
        I=s.MOI()[2,2]
        E=I*s.thetadot**2/2+s.M*np.dot(s.vcm,s.vcm)/2
        Evolution(s)
        if i==timearray[j]:
            print 'Energy(t='+str(t)+')= ' + str(E)
            j=j+1
        i=i+1
        t=t+dt

    
        
def main():
    t=0.0
    dt=0.0005
    s=Untitled.Square(np.array([5.9,5.9]),0.1,np.array([6,3]),1,1000,0.2,30,2,dt)
    s.createsquare()
    s.createboundary()
    timearray=[0,5,15,100,200] # enter time in units of dt
    Energy(s,dt,timearray,t)

    
if __name__=='__main__':
    main()        
        
        
        
            
    
