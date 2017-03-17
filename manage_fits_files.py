import numpy as np
import random
def balanced_sample(SPT,band1,num_elements):
	x = np.zeros(10)
	
	c = list(zip(SPT,band1))
	random.shuffle(c)
	SPT, band1 = zip(*c)
	
	SPTs = []
	band1s = []
	for i in range(0,len(SPT)):
		if(x[int(SPT[i])] != num_elements and int(SPT[i]) != 9):
			SPTs.append(SPT[i])
			band1s.append(band1[i])	
			x[int(SPT[i])] = x[int(SPT[i])] + 1
	
	return SPTs,band1s
