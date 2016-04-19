import os

base = 0.001
for i in range(1,100):
	#print i*base
	command = "./admmLRConsensus_cycle %f"%(i*base)
	os.system(command)
