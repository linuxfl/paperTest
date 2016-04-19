#!/bin/bash

gcc -I/usr/local/include/ -c admmLRConsensus_butterfly.c 
gcc -L/usr/local/include/  admmLRConsensus_butterfly.o -o admmLRConsensus_butterfly -lgsl -lgslcblas -lm
