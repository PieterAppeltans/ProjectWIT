
def parse_input(filename):
	vi = open(filename+'.node','r')
	vo = open('vertices.dat','w')
	ti = open(filename+'.ele','r')
	to = open('triangles.dat','w')
	pi = open(filename+'.poly','r')
	po = open('boundary.dat','w')
	vertices = []
	elements = []
	for line in vi:
		arr = line.split()
		t = arr[1:3]
		for l in range(0,2):
			t[l] = float(t[l])/100.
		vertices.append(t)

	for line in ti:
		arr = line.split()
		s = arr[1:]
		for i in range(0,len(s)):
			s[i] = int(s[i])-1
		elements.append(s)

	vi.close()
	vo.close()
	ti.close()
	to.close()
	pi.close()
	po.close()
	return vertices,elements
