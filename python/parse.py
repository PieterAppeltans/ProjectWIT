
def parse_input(filename):
	vi = open('../triangle/'+ filename + '.node','r')
	ti = open('../triangle/'+ filename + '.ele','r')
	pi = open('../triangle/'+ filename + '.poly','r')
	vertices = []
	elements = []
	vi_lines = vi.readlines()
	for line in vi_lines[1:-1]:
		arr = line.split()
		t = arr[1:3]
		for l in range(0,2):
			t[l] = float(t[l])/100.
		vertices.append(t)
	ti_lines = ti.readlines()
	for line in ti_lines[1:-1]:
		arr = line.split()
		s = arr[1:]
		for i in range(0,len(s)):
			s[i] = int(s[i])-1
		elements.append(s)

	vi.close()
	ti.close()
	pi.close()

	return vertices,elements
