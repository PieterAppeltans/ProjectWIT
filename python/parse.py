import numpy as np
def parse_u_v():
	file_u = open('../triangle/result_u.out')
	file_v = open('../triangle/result_v.out')
	u = np.array([])
	v = np.array([])
	for line in file_u:
		u = np.append(u,float(line))
	for line in file_v:
		v= np.append(v,float(line))
	return u,v
def parse_input(filename):
	vi = open('../triangle/'+ filename + '.node','r')
	ti = open('../triangle/'+ filename + '.ele','r')
	pi = open('../triangle/'+ filename + '.poly','r')
	vertices = np.empty((0,2), float)
	elements = np.empty((0,2), float)
	vi_lines = vi.readlines()
	for line in vi_lines[1:-1]:
		arr = line.split()
		t = arr[1:3]
		for l in range(0,2):
			t[l] = float(t[l])

		vertices = np.append(vertices,np.array([t]),axis=0)
	ti_lines = ti.readlines()
	for line in ti_lines[1:-1]:
		arr = line.split()
		s = arr[1:3]
		for i in range(0,len(s)):
			s[i] = int(s[i])-1
		elements = np.append(elements,np.array([s]),axis=0)

	vi.close()
	ti.close()
	pi.close()

	return vertices,elements
