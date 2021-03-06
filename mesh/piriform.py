import numpy as np
from math import sin,cos

def f_y(t):
    return 100-(1+sin(t))

def f_x(t):
    return -40*(cos(t)*(1+sin(t)))

count = 1
print "250 2 0 1"
for t in np.linspace(np.pi/2,np.pi*1.5,250):
	print count," ",f_x(t)," ",f_y(t)," ",1
	count += 1
print "250 0 1"
for i in range(1,250):
	print i," ",i," ",i+1," 1"
print "250 1 250 1"
print "0"
