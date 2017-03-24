import numpy as np
count = 1
print "250 2 0 1"
for t in np.linspace(-np.pi/2,np.pi/2,250):
	print count," ",50*np.cos(t)," ",50*np.sin(t)," ",1
	count += 1
print "250 0 1"
for i in range(1,250):
	print i," ",i," ",i+1," 1"
print "250 1 250 1"
print "0"
