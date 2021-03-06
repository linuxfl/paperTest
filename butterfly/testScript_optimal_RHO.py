import os
import sys
import socket

USER = "root"
numberofprocessor = 16

if __name__ == "__main__":
	#if len(sys.argv) < 3:
	#	print "need two paraments!!"
	#	exit()
	
	hostname = socket.gethostname()
	
	nodelist = []
	fp = open("hostfile")
	for node in fp.readlines():
		node = node.split(":")
		nodelist.append(node[0])
		
	for node in nodelist:
		if node != hostname:
			codeCpyCmd1 = "scp admmLRConsensus_butterfly.py %s@%s:/%s/paperTest/butterfly/"%(USER,node,USER)
			codeCpyCmd2 = "scp mpiNode.py %s@%s:/%s/paperTest/butterfly/"%(USER,node,USER)
			os.system(codeCpyCmd1)
			os.system(codeCpyCmd2)

	for j in range(1,2):
		for i in range(1,11):
			print "rho = ",float(i*0.01+j*0.1)
			codeRunCmd = "mpirun -f hostfile -np %d python admmLRConsensus_butterfly.py %f"%(numberofprocessor,float(i*0.01+j*0.1))
			os.system(codeRunCmd)
