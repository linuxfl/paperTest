#
#   y = x1*a1 + x2*a2 + x3*a3 +...+ xn*an
#   Ax = b
#
from numpy import *
import numpy as np
import sys

def generateData(sampleNumber,paraNumber):
	
	s = np.random.uniform(-1.9061842247,1.14401391414,paraOfNumber)
	s = mat(s).T
	np.savetxt("solution.dat",s)
	Astr = "lineardataA.dat";bstr = "lineardatab.dat"
	a = np.random.uniform(-0.375111355247,0.322563975491,paraOfNumber * sampleOfNumber)
	a = a.reshape((sampleOfNumber,paraOfNumber))
	np.savetxt(Astr,a)
	b = a * s
	np.savetxt(bstr,b)
		
if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "argv must two praraments likes (samepleNubmer,paraNumber)!"
	sampleOfNumber = int(sys.argv[1])
	paraOfNumber = int(sys.argv[2])
	generateData(sampleOfNumber,paraOfNumber)
