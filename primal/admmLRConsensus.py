### Solve a linear regression, i.e.,
###            Ax = b
### minimize  (1/2)||Ax - b||_2^2.


from mpi4py import MPI
from numpy import *
import numpy as np
import time
import sys

iterMax = 4000
ABSTOL = 1e-10

def testErr(A,b,x):
    y = A*x
    return np.linalg.norm(y - b)

def objectfuction(A,x,b):
    return 0.5 * np.linalg.norm(A*x - b)

def admmLR(A,b,rho,alpha,comm,rank,size):
    matA = mat(A);matb = mat(b).T
    solution = mat(np.loadtxt("../data/solution.dat")).T
    n = shape(matA)[1]
    meanxold = zeros((n,1),dtype = np.float)
    meanx = zeros((n,1),dtype = np.float)
    x = zeros((n,1),dtype = np.float)
    y = zeros((n,1),dtype = np.float)
    lemon = (matA.T * matA + 2 * size * rho * eye((n))).I
    Atb = matA.T * matb
    itercount = 0
    commuTimeSum = 0.0
    while itercount < iterMax:
        #update x:
        x = lemon * (Atb + rho * size * x +rho * meanx - y)
        commuTime0 = time.time()
        comm.Allreduce(x,meanxold,op = MPI.SUM)
        meanxold = meanxold - x
        commuTime1 = time.time()
	
        #update meanx:
        meanx = meanxold/1.0
        meanxold.fill(0)

        #update y:
        y = y + rho * (x * size - meanx)
        
        #tol
        tol = np.linalg.norm(x - solution)/np.linalg.norm(solution)
        alltol = zeros((1,1))
        comm.Allreduce(tol,alltol,op = MPI.MIN)
        if rank == 0:
			logLine =  "itercout:%d primal %.15f objectionfuction --- >%.15f"%(itercount,alltol,objectfuction(matA,meanx/size,matb))
			#print logLine
        if alltol < ABSTOL:
			break
        itercount = itercount + 1
        commuTimeSum += commuTime1 - commuTime0
    return x,commuTimeSum,itercount

if __name__ == "__main__":
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	Adir = "../data/A%d.dat"%(rank%20);bdir = "../data/b%d.dat"%(rank%20)
	A = np.loadtxt(Adir);b = np.loadtxt(bdir)
	if rank == 0:
		time0 = time.time()
	#8 2.5
	#16 2.3
	x,ct,it = admmLR(A,b,float(sys.argv[1]),1.0,comm,rank,size-1)
	if rank == 0:
		time1 = time.time()
		print "itercout:",it,"rank %d"%(rank),"Elapsed time is","%4.3f"%(time1 - time0),\
		"Communication time is %4.9f"%ct,"seconds","average Communication time %4.9f"%(ct/it)
