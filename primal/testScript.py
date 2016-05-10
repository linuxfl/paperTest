import os
import sys
import socket
USER = "root"

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "need two paraments!!"
		exit()
	hostname = socket.gethostname()
	nodelist = []
	fp = open("hostfile")
	for node in fp.readlines():
		node = node.split(":")
		nodelist.append(node[0])
		
	for node in nodelist:
		if node != hostname:
			codeCpyCmd = "scp admmLRConsensus.py %s@%s:/%s/paperTest/primal/"%(USER,node,USER)
			os.system(codeCpyCmd)

	codeRunCmd = "mpirun -f hostfile -np %d python admmLRConsensus.py %f"%(int(sys.argv[1]),float(sys.argv[2]))
	os.system(codeRunCmd)
