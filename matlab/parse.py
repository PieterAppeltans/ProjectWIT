# Script to store mesh data in three .dat files for (manual, not gui) Matlab fem.m

filename = "../mesh/circle.1"
vi = open(filename+'.node','r')
vo = open('vertices.dat','w')
ti = open(filename+'.ele','r')
to = open('triangles.dat','w')
pi = open(filename+'.poly','r')
po = open('boundaries.dat','w')

vertices = []
nb_vertices = int(next(vi).split()[0])
for i in range(nb_vertices):
	arr = next(vi).split()
	t = arr[1:3]
	for l in range(0,2):
		t[l] = str(float(t[l])/1000)
	vertices.append(t)
	s = " ".join(t)
	vo.write(s)
	vo.write("\n")

vertices = []
nb_triangles = int(next(ti).split()[0])
for i in range(nb_triangles):
	arr = next(ti).split()
	s = " ".join(arr[1:])
	to.write(s)
	to.write("\n")

vertices = []
nb_lines = int(next(pi).split()[0])
for i in range(nb_lines):
	arr = next(pi).split()
	print (float(vertices[int(arr[1])-1][0]) <= 10**-15) and (float(vertices[int(arr[2])-1][0]) <= 10**-15)
	if not (float(vertices[int(arr[1])-1][0]) <= 10**-15 and float(vertices[int(arr[2])-1][0]) <= 10**-15):
		s = " ".join(arr[1:3])
		po.write(s)
		po.write("\n")

vi.close()
vo.close()
ti.close()
to.close()
pi.close()
po.close()
