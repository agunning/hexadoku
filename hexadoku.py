import math as math
import numpy as np
from time import perf_counter

''' USAGE: just run python3 hexadoku.py on command line '''

exclusion_grid=np.zeros((49,49),dtype=int)
#exclusion_grid[i,j]=1 <=> cell i and cell j must be distinct

#to build this operate in eisenstein integers omega is such that omega*2=-1-omega

for i in range(7):
    exclusion_grid[7*i:7*(i+1),7*i:7*(i+1)]=1


def eisenmult(A,B):
    a1,a2=A
    b1,b2=B
    return (a1*b1-a2*b2,a1*b2+b1*a2-b2*a2)

indices=[(0,0)]
X=(-1,-1)
for i in range(6):
    X=eisenmult(X,(1,1))
    indices.append(X)

locs = np.zeros((49,2),dtype=int)

for i in range(49):
    locs[i]= np.array(indices[i%7])+np.array(eisenmult([3,1],indices[i//7]))

X=np.array(
    [[2,-1],[0,3]])
Y=np.array([7,12])


W=np.zeros(65)
for i in range(7):
    W[2**i]=i

W[0]=-1

def print_soln(soln,binconvert=True):
    if binconvert:
        soln=W[soln]
    printthings=-np.ones((15,25),dtype=int)
    for i in range(49):
        t=X@locs[i]+Y
        printthings[t[0],t[1]]=soln[i]
    for i in range(15):
        s=''
        for j in range(25):
            if printthings[i,j]==-1:
                s+=" "
            else:
                s+=str(printthings[i,j])
        print(s)

exclusion_grid+=(locs[:,0]==locs[:,0].reshape(-1,1))
exclusion_grid+=(locs[:,1]==locs[:,1].reshape(-1,1))
exclusion_grid+=((locs[:,0]-locs[:,1])==(locs[:,0]-locs[:,1]).reshape(-1,1))

exclusion_grid=np.minimum(exclusion_grid,1)
if False:
    for i in range(6):
        for (j,k) in np.ndindex((2,2)):
            exclusion_grid[7*(i+1)+(i+j)%6+1,7*((i+1)%6+1)+(i+k)%6+1]=1
            exclusion_grid[7*((i+1)%6+1)+(i+k)%6+1,7*(i+1)+(i+j)%6+1]=1


l=2**np.arange(7,dtype=int)
solmatrix=np.zeros((1,49),dtype=int)
solmatrix[0,:7]=l

U=np.zeros(129,dtype=int)

for i in range(7):
    U[2**i:2**(i+1)]=1+U[:2**i]
U[-1]=-1


Z=np.bitwise_or.reduce(solmatrix.reshape(-1,49,1)*exclusion_grid,axis=1)
Z=np.maximum(128*(solmatrix>0),Z)
def do_thing(solmatrix,Y):
    m=solmatrix.shape[0]

    T=np.argmax(U[Y],axis=1)
    vals = np.take_along_axis(Y,T.reshape(-1,1),axis=1).flatten()
    newm=(7-U[vals]).sum()
    newsolmatrix=np.repeat(solmatrix,7-U[vals],axis=0).flatten()
    newY=np.repeat(Y,7-U[vals],axis=0).flatten()
    indices=np.repeat(T,7-U[vals].flatten(),axis=0)+np.arange(newm)*49
    reps=np.tile(l,m)
    reps= reps[(vals.repeat(7)&reps==0)]
    newY|=(reps.reshape(-1,1)*exclusion_grid[np.repeat(T,7-U[vals])]).flatten()
    newY[indices]=128
    newY=np.minimum(128,newY)
    newsolmatrix[indices]=reps
    return newsolmatrix.reshape(-1,49),newY.reshape(-1,49)
    print_soln(U[Y[0,:]])
    print(T)

def run_loop(iterations=42):
    t1_start=perf_counter()
    Xiter=solmatrix
    Yiter=Z
    for i in range(iterations):
        Xiter,Yiter=do_thing(Xiter,Yiter)
    
        #print("Cell "+str(i+8)+", "+str(Xiter.shape[0])+" possibilities")\
    for i in range(Xiter.shape[0]):
        print_soln(Xiter[i])
    t1_stop= perf_counter()

    print("time taken: "+str(t1_stop-t1_start)+" seconds")
    return Xiter
r=run_loop()
