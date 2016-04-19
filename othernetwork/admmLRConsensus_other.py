### Solve a linear regression, i.e.,
###            Ax = b
### minimize  (1/2)||Ax - b||_2^2.


from mpi4py import MPI
from numpy import *
import numpy as np
import mpiNode
import time
import sys

iterMax = 4000
ABSTOL = 1e-10
ALLCONNECT = 0
CYCLE = 1

def testErr(A,b,x):
    y = A*x
    return np.linalg.norm(y - b)

def objectfuction(A,x,b):
    return 0.5 * np.linalg.norm(A*x - b)

def admmLR(A,b,rho,alpha,localcomm,mastercomm,comm,rank,size,connectTOPO):
	tolList = []
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
	while itercount < iterMax:
		#update x:
		x = lemon * (Atb + rho * size * x +rho * meanx - y)
		commuTime0 = time.time()
		#localcomm.Allreduce(x,meanxold,op = MPI.SUM)
		#meanxold = meanxold - x
		localcomm.Reduce(x,meanxold,op=MPI.SUM,root=0)
		
		if mastercomm != None:
			#local all connect topo
			if connectTOPO == ALLCONNECT:
				mastercomm.Allreduce(meanxold,meanxtmp,op = MPI.SUM)
				meanxold = meanxtmp
				
			#cycle network topo
			if connectTOPO == CYCLE:
				curRank = mastercomm.Get_rank()
				preRank = (curRank + mastercomm.Get_size() - 1)%mastercomm.Get_size()
				nextRank = (curRank + mastercomm.Get_size() + 1)%mastercomm.Get_size()
				mastercomm.Sendrecv(meanxold,preRank,0,meanxtmp,nextRank,0)
				mastercomm.Sendrecv(meanxold,nextRank,1,meanxtmp1,preRank,1)
				meanxold = meanxold + meanxtmp1 + meanxtmp
				
		meanxold = localcomm.bcast(meanxold,root = 0)
		meanxold = meanxold - x
		
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
		tolList.append(alltol[0,0])
		if alltol < ABSTOL:
			break
		if rank == 0:
			logLine = "itercount:%d primal %.15f objectionfuction --- >%0.15f"%(itercount,alltol,objectfuction(matA,meanx/size,matb))
			#print logLine
		itercount = itercount + 1
		commuTimeSum += commuTime1 - commuTime0
	return x,commuTimeSum,tolList

if __name__ == "__main__":
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	mpiNode.MPIN_init(comm,rank)
	if len(sys.argv) < 3:
		print "need two argument!!"
		comm.Abort()
	connectTOPO = int(sys.argv[1])
	localcomm = mpiNode.MPIN_get_local_comm(comm,rank)
	mastercomm = mpiNode.MPIN_get_master_comm(comm,rank)
	nodeid = mpiNode.MPIN_get_node_by_rank(rank)
	masterRank = mpiNode.MPIN_get_master_rank(nodeid)
	if connectTOPO == ALLCONNECT:
		localsize = size - 1
	if connectTOPO == CYCLE:
			localsize = 3 * localcomm.Get_size() - 1
	if rank != masterRank:
		mastercomm = None
	Adir = "../data/A%d.dat"%(rank%20);bdir = "../data/b%d.dat"%(rank%20)
	A = np.loadtxt(Adir);b = np.loadtxt(bdir)
	time0 = time.time()
	#connect 8 11.67
	#cycle 8 14.2
	#butterfly 8 18.5
	x,ct,tolList = admmLR(A,b,float(sys.argv[2]),1.0,localcomm,mastercomm,comm,rank,localsize,connectTOPO)
	if rank == 0:
		time1 = time.time()
		print "itercount",len(tolList),"rank %d"%rank," Elapsed time is","%4.3f"%(time1 - time0),"seconds","Communication time is %4.9f"%ct,"seconds","average Communication Time is %4.9f"%(ct/len(tolList))
