### Solve a linear regression, i.e.,
###            Ax = b
### minimize  (1/2)||Ax - b||_2^2.


from mpi4py import MPI
from numpy import *
import numpy as np
import mpiNode
import time
import math
import sys

iterMax = 4000
ABSTOL = 1e-10

def objectfuction(A,x,b):
    return 0.5 * np.linalg.norm(A*x - b)

def admmLR(A,b,rho,localcomm,mastercomm,comm,rank,size):
	
	matA = mat(A);matb = mat(b).T
	n = shape(matA)[1]
	solution = mat(np.loadtxt("../data/solution.dat")).T
	meanxold = zeros((n,1),dtype = np.float)
	meanxtmp = zeros((n,1),dtype = np.float)
	meanxtmp1 = zeros((n,1),dtype = np.float)
	meanx = zeros((n,1),dtype = np.float)
	x = zeros((n,1),dtype = np.float)
	y = zeros((n,1),dtype = np.float)
	lemon = (matA.T * matA + 2 * size * rho * eye((n))).I
	Atb = matA.T * matb
	
	itercount = 0
	commuTimeSum = 0.0
	
	if mastercomm != None:
		numberofNodeLog = int(math.log(mastercomm.Get_size(),2))
		curRank = mastercomm.Get_rank()

	while itercount < iterMax:
		#update x:
		x = lemon * (Atb + rho * size * x +rho * meanx - y)
		commuTime0 = time.time()
		localcomm.Allreduce(x,meanxold,op = MPI.SUM)
		meanxold = meanxold - x
		
		if mastercomm != None:
			#butterfly algorithm
			k = itercount%numberofNodeLog
			if(k == 0):
				tempRank = curRank
			flag = tempRank & 1;
			tempRank = tempRank >> 1
			if flag == 0:
				swapRank = curRank | 1<<(k)
			else:
				swapRank = curRank & (~(1<<(k)))
			mastercomm.Sendrecv(meanxold + x,swapRank,0,meanxtmp,swapRank,0)
			meanxold = meanxold + meanxtmp
			
		#update meanx:
		commuTime1 = time.time()
		meanx = meanxold/1.0
		meanxold.fill(0)
		meanxtmp.fill(0)
		
		#update y:
		y = y + rho * (x * size - meanx)
		
		#tol
		tol = np.linalg.norm(x - solution)/np.linalg.norm(solution)
		alltol = zeros((1,1))
		comm.Allreduce(tol,alltol,MPI.MIN)
		
		if alltol < ABSTOL:
			break
		if rank == 0:
			logLine = "itercount:%d primal %.15f objectionfuction --- >%0.15f"%(itercount,alltol,objectfuction(matA,meanx/size,matb))
			print logLine
		itercount = itercount + 1
		commuTimeSum += commuTime1 - commuTime0
	return x,commuTimeSum,itercount

if __name__ == "__main__":
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	mpiNode.MPIN_init(comm,rank)
	if len(sys.argv) < 2:
		print "need two argument!!"
		comm.Abort()
	localcomm = mpiNode.MPIN_get_local_comm(comm,rank)
	mastercomm = mpiNode.MPIN_get_master_comm(comm,rank)
	nodeid = mpiNode.MPIN_get_node_by_rank(rank)
	masterRank = mpiNode.MPIN_get_master_rank(nodeid)
	if rank == masterRank:
		#butterfly network topo
		localsize = 2 * localcomm.Get_size() -1
	else:
		localsize = localcomm.Get_size()-1
		mastercomm = None
	Adir = "../data/A%d.dat"%(rank%20);bdir = "../data/b%d.dat"%(rank%20)
	A = np.loadtxt(Adir);b = np.loadtxt(bdir)
	time0 = time.time()
	#connect 8 11.67
	#cycle 8 14.2
	#butterfly 8 18.5
	x,ct,it = admmLR(A,b,float(sys.argv[1]),localcomm,mastercomm,comm,rank,localsize)
	time1 = time.time()
	if rank == 0:
		print "itercount:",it,"rank %d"%rank," Elapsed time is","%4.3f"%(time1 - time0),"seconds",\
		"Communication time is %4.9f"%ct,"seconds","average Communication time %4.9f"%(ct/it)
	
	
