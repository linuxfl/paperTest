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
BUTTERFLY = 2

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
	itercout = 0
	commuTimeSum = 0.0
	while itercount < iterMax:
		#update x:
		x = lemon * (Atb + rho * size * x +rho * meanx - y)
		commuTime0 = time.time()
		localcomm.Allreduce(x,meanxold,op = MPI.SUM)
		meanxold = meanxold - x
		if mastercomm != None:
			print "rank:",rank
			#local all connect topo
			if connectTOPO == ALLCONNECT:
				mastercomm.Allreduce((meanxold+x),meanxtmp,op = MPI.SUM)
				meanxold = meanxtmp - x
				
			#butterfly network topo
			if connectTOPO == BUTTERFLY:
				curRank = mastercomm.Get_rank()	
				if itercout % 2 == 0:
					if curRank % 2 == 0:
						swapRank = (curRank + mastercomm.Get_size() + 1)%mastercomm.Get_size()
					else:
						swapRank = (curRank + mastercomm.Get_size() - 1)%mastercomm.Get_size()
				else:
					if curRank % 2 == 0:
						swapRank = (curRank + mastercomm.Get_size() + 2)%mastercomm.Get_size()
					else:
						swapRank = (curRank + mastercomm.Get_size() - 2)%mastercomm.Get_size()
				mastercomm.Sendrecv(meanxold + x,swapRank,0,meanxtmp,swapRank,0)
				meanxold = meanxold + meanxtmp		
					
			#cycle network topo
			if connectTOPO == CYCLE:
				curRank = mastercomm.Get_rank()
				preRank = (curRank + mastercomm.Get_size() - 1)%mastercomm.Get_size()
				nextRank = (curRank + mastercomm.Get_size() + 1)%mastercomm.Get_size()
				mastercomm.Sendrecv((meanxold+x),preRank,0,meanxtmp,nextRank,0)
				mastercomm.Sendrecv((meanxold+x),nextRank,1,meanxtmp1,preRank,1)
				meanxold = meanxold + meanxtmp1 + meanxtmp
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
			print logLine
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
	if rank == masterRank:
		#all nodes  connect
		if connectTOPO == ALLCONNECT:
			#localsize = localcomm.Get_size()+mpiNode.MPIN_get_node_size()-2
			localsize = size - 1
		#butterfly network topo
		if connectTOPO == BUTTERFLY:
			localsize = 2 * localcomm.Get_size() -1
		#cycle network topo
		if connectTOPO == CYCLE:
			localsize = 3 * localcomm.Get_size() - 1
	else:
		localsize = localcomm.Get_size()-1
		mastercomm = None
	Adir = "../data/A%d.dat"%(rank%20);bdir = "../data/b%d.dat"%(rank%20)
	A = np.loadtxt(Adir);b = np.loadtxt(bdir)
	if rank == 0:
		time0 = time.time()
	#connect 8 11.67
	#cycle 8 14.2
	#butterfly 8 18.5
	x,ct,tolList = admmLR(A,b,float(sys.argv[2]),1.0,localcomm,mastercomm,comm,rank,localsize,connectTOPO)
	if rank == 0:
		time1 = time.time()
		print "itercout",len(tolList),"rank %d"%rank," Elapsed time is","%4.3f"%(time1 - time0),"seconds","Communication time is %4.9f"%ct,"seconds","average Communication Time is %4.9f"%(ct/len(tolList))
