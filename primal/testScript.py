import os
import sys
USER = "fangling"

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
		codeCpyCmd = "scp admmLRConsensus.py %s@%s:/home/%s/paperTest/primal/"%(USER,node,USER)
		os.system(codeCpyCmd1)

	codeRunCmd = "mpirun -f hostfile -np %d python admmLRConsensus.py %f"%(int(sys.argv[1]),float(sys.argv[2]))
	os.system(codeRunCmd)
