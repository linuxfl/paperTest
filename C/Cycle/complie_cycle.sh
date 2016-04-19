#!/bin/bash

gcc -I/usr/local/include/ -c admmLRConsensus_cycle.c 
gcc -L/usr/local/include/  admmLRConsensus_cycle.o -o admmLRConsensus_cycle -lgsl -lgslcblas -lm
