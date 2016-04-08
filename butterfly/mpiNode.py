import socket
from numpy import *
import numpy as np

hpmap = []
hostlist = []
numberofnode = 0
numberofprocess = 0
process = []

def MPIN_init(comm,rank):
	global numberofnode
	global numberofprocess
	global process
	
	hostname = socket.gethostname()
	
	dt = dtype([('hostname','S10'),('rank',np.int)])
	localhost = np.array([(hostname,rank)],dtype = dt)
	localhost = comm.gather(localhost,root = 0)
	localhost = comm.bcast(localhost,root = 0)
	numberofprocess = len(localhost)
	
	for i in range(numberofprocess):
		if localhost[i]['hostname'] not in hostlist:
			hostlist.append(localhost[i]['hostname'])
			
	numberofnode = len(hostlist)
	
	for i in range(numberofnode):
		l = []
		hpmap.append(l)
		
	process = range(numberofprocess)
	
	for i in range(numberofprocess):
		k = hostlist.index(localhost[i]['hostname'],0,numberofnode)
		hpmap[k].append(int(localhost[i]['rank']))
		process[i] = k

def MPIN_get_node_size():
	return numberofnode
	
def MPIN_get_node_process_size(nodeid):
	if nodeid >= numberofnode or nodeid < 0:
		print "MPIN_get_node_process_size:nodeid error!"
		return -1
	else:
		return len(hpmap[nodeid])
		 
def MPIN_get_master_rank(nodeid):
	if nodeid >= numberofnode or nodeid < 0:
		print "MPIN_get_master_rank:nodeid error!"
		return -1
	else:
		return hpmap[nodeid][0]

def MPIN_get_master_rank_list(nodeid,degree):
	if nodeid >= numberofnode or nodeid < 0:
		print "MPIN_get_master_rank:nodeid error!"
		return None
	else:
		numberp = MPIN_get_node_process_size(nodeid)
		if numberp != -1:
			masterlist = []
			if degree <= numberp:				
				for i in range(degree):
					masterlist.append(hpmap[nodeid][i])
			else:
				masterlist.append(hpmap[nodeid][0])
			return masterlist
		else:
			return None

def MPIN_get_node_process_rank(nodeid,loc):
	if nodeid >= numberofnode or nodeid < 0:
		print "MPIN_get_node_process_rank:nodeid error!"
		return -1
	elif loc >= len(hpmap[nodeid]) or loc < 0:
		print "MPIN_get_node_process_rank:loc error!"
		return -1
	else:
		return hpmap[nodeid][loc]
	
def MPIN_get_node_by_rank(rank):
	if rank >= numberofprocess or rank < 0:
		print "MPIN_get_node_by_rank:rank error!"
		return -1
	else:
		return process[rank]
	
def MPIN_get_local_comm(comm,rank):
	color = MPIN_get_node_by_rank(rank)
	localComm = comm.Split(color, 0)
	return localComm

def MPIN_get_master_comm(comm,rank):
	if numberofnode > 1:
		nodeid = MPIN_get_node_by_rank(rank)
		masterRank = MPIN_get_master_rank(nodeid)
		glist = []
		if masterRank != -1:
			group = comm.Get_group()
			for i in range(numberofnode):	
				mr = MPIN_get_master_rank(i)
				if mr != -1:
					glist.append(mr)
			newgroup = group.Incl(glist)
			masterComm = comm.Create(newgroup)
			return masterComm
		else:
			return None
	else:
		return None

def MPIN_get_master_comm(comm,rank):
	if numberofnode > 1:
		nodeid = MPIN_get_node_by_rank(rank)
		masterRank = MPIN_get_master_rank(nodeid)
		glist = []
		if masterRank != -1:
			group = comm.Get_group()
			for i in range(numberofnode):	
				mr = MPIN_get_master_rank(i)
				if mr != -1:
					glist.append(mr)
			newgroup = group.Incl(glist)
			masterComm = comm.Create(newgroup)
			return masterComm
		else:
			return None
	else:
		return None

def MPIN_get_master_comm_degree(comm,rank,degree):
	if numberofnode > 1:
		glist = []
		group = comm.Get_group()
		nodeid = MPIN_get_node_by_rank(rank)
		for i in range(numberofnode):	
			mrlist = MPIN_get_master_rank_list(i,degree)
			if mrlist != None:
				glist = glist + mrlist	
		#print glist
		newgroup = group.Incl(glist)
		masterComm = comm.Create(newgroup)
		return masterComm
	else:
		return None
