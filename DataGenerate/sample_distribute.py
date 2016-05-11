from numpy import *
import numpy as np
import os

#user and sample information
USER = "root"
SAMPLE_FILE_DIR = "/paperTest/data"
SOURCE_SAMPLE_DIR_A = "/root/paperTest/DataGenerate/lineardataA.dat"
SOURCE_SAMPLE_DIR_b = "/root/paperTest/DataGenerate/lineardatab.dat"
COL = 500
HOSTFILE = "../hostfile"
METHOD = "process_type" #method is 'process_type' or 'node_type'

numberofprocess = int(input("please input the number of the process:"))
fileline = int(input("please input the number line of the source data file:"))
eachline = int(fileline/numberofprocess)
numberofnode = 0
hostdict = {}
hostlist = []
hpmap = []
strline = ""

class getoutofloop(Exception): pass

#get the host information
fp_host = open(HOSTFILE,"r")
strline = fp_host.readline()
while strline:
    pos = strline.find(":")
    host = strline[0:pos]
    if strline[(len(strline)-1)] == '\n':
        nop = strline[(pos+1):(len(strline)-1)]
    else:
        nop = strline[(pos+1):len(strline)]
    hostdict[host] = nop
    hostlist.append(host)
    strline = fp_host.readline()
fp_host.close()

#ceate the map of the node and process
#init the the map of host and process
for i in range(0,int(len(hostlist))):
    l = []
    hpmap.append(l)
    
rank = 0
filter_round = 0
try:
    while True:
        for i in range(0,int(len(hostlist))):
            j = 0
            while j < int(hostdict[hostlist[i]]):
                if rank >= numberofprocess:
					if j == 0:
						i = i - 1
					raise getoutofloop()
                j += 1
                hpmap[i].append(rank)
                rank += 1
	filter_round = 1
except getoutofloop:
    pass

#get the number of arranging node
if filter_round == 1:
	numberofnode = int(len(hostlist))
elif filter_round == 0:
	numberofnode = i+1
else:
	print("error flag value!")
print(numberofnode)

#clean up the sample direction
for i in range(0,numberofnode):
    command = "ssh "+hostlist[i]+" rm "+SAMPLE_FILE_DIR+"/* -rf"
    print(command)
    os.system(command)
    

#read data file and split it by number of process

dataA = np.loadtxt(SOURCE_SAMPLE_DIR_A)
dataA = dataA.reshape((fileline,COL))
datab = np.loadtxt(SOURCE_SAMPLE_DIR_b)
datab = datab.reshape((fileline,1))

for i in range(0,numberofprocess):
	desA = "A%d.dat"%i;desb = "b%d.dat"%i
	np.savetxt(desA,dataA[i*eachline:(i+1)*eachline,:])
	np.savetxt(desb,datab[i*eachline:(i+1)*eachline,:])
	
for i in range(0,numberofnode):
	dessolu = "solution.dat"
	commands = "scp "+ dessolu +" "+USER+"@"+str(hostlist[i])+":"+SAMPLE_FILE_DIR
	os.system(commands)
	for j in hpmap[i]:
		#example :scp 1 fangling@node1:/home/fangling/xxxx/
		desA = "A%d.dat"%j;desb = "b%d.dat"%j
		commandA = "scp "+ desA +" "+USER+"@"+str(hostlist[i])+":"+SAMPLE_FILE_DIR
		commandb = "scp "+ desb +" "+USER+"@"+str(hostlist[i])+":"+SAMPLE_FILE_DIR
		os.system(commandA)
		os.system(commandb)
	
#delete the every file in local
for i in range(0,numberofprocess):
    desA = "A%d.dat"%i;desb = "b%d.dat"%i
    commandA = "rm "+desA+" -rf"
    commandb = "rm "+desb+" -rf"
    os.system(commandA) 
    os.system(commandb)
print("success...")
