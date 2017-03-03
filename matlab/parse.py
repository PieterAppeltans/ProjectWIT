vi = open('circle.1.node','r')
vo = open('vertices.dat','w')
ti = open('circle.1.ele','r')
to = open('triangles.dat','w')
pi = open('circle.1.poly','r')
po = open('boundary.dat','w')
vertices = []
for line in vi:
	arr = line.split()
	t = arr[1:3]
	for l in range(0,2):
		t[l] = str(float(t[l]))
	vertices.append(t)
	s = " ".join(t)
	vo.write(s)
	vo.write("\n")

for line in ti:
	arr = line.split()
	s = " ".join(arr[1:])
	to.write(s)
	to.write("\n")
for line in pi:
	arr = line.split()
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
