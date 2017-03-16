import numpy as np
count = 1
print "25 2 0 1"
for t in np.linspace(-np.pi/2,np.pi/2,25):
	print count," ",0.05*np.cos(t)," ",0.05*np.sin(t)," ",1
	count += 1
print "25 0 1"
for i in range(1,25):
	print i," ",i," ",i+1," 1"
print "25 1 25 1"
print "0"
