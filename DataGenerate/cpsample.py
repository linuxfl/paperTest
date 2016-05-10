import os
import sys
import socket

USER = "root"

if __name__ == "__main__":
	
	hostname = socket.gethostname()
	nodelist = []
	fp = open("../primal/hostfile")
	for node in fp.readlines():
		node = node.split(":")
		nodelist.append(node[0])
		
	for node in nodelist:
		if node != hostname:
			print "copy data to %s"%(node)
			codeCpyCmd = "scp ../data/* %s@%s:/%s/paperTest/data"%(USER,node,USER)
			os.system(codeCpyCmd)
