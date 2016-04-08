import os
import sys
USER = "fangling"
rho = 1

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "need two paraments!!"
		exit()

	nodelist = []
	fp = open("hostfile")
	for node in fp.readlines():
		node = node.split(":")
		nodelist.append(node[0])
		
	for node in nodelist:
		codeCpyCmd1 = "scp admmLRConsensus_butterfly.py %s@%s:/home/%s/paperTest/butterfly/"%(USER,node,USER)
		codeCpyCmd2 = "scp mpiNode.py %s@%s:/home/%s/paperTest/butterfly/"%(USER,node,USER)
		os.system(codeCpyCmd1)
		os.system(codeCpyCmd2)

	codeRunCmd = "mpirun -f hostfile -np %d python admmLRConsensus_butterfly.py %f"%(int(sys.argv[1]),float(sys.argv[2]))
	os.system(codeRunCmd)
