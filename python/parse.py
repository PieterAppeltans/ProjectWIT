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
	vertices = np.empty((0,2), float)
	elements = np.empty((0,3), int)
	vi_lines = vi.readlines()
	for line in vi_lines[1:-1]:
		arr = line.split()
		t = arr[1:3]
		for l in range(0,2):
			t[l] = float(t[l])/1000

		vertices = np.append(vertices,np.array([t]),axis=0)
	ti_lines = ti.readlines()
	for line in ti_lines[1:-1]:
		arr = line.split()
		s = arr[1:]
		for i in range(0,len(s)):
			s[i] = int(s[i])-1
		elements = np.append(elements,np.array([s]),axis=0)
	print elements
	vi.close()
	ti.close()


	return vertices,elements
def write_matlab(filename,vertices,elements):
	pi = open('../triangle/'+ filename + '.1.poly','r')
	vo = open('../matlab/vertices.dat','w')
	to = open('../matlab/triangles.dat','w')
	po = open('../matlab/boundaries.dat','w')
	pi_lines = pi.readlines()
	for v in vertices:
		v_str = np.char.mod('%f', v)
		s = " ".join(v_str)
		vo.write(s)
		vo.write("\n")
	for e in elements:
		e_str = np.char.mod('%i',e+1)
		s = " ".join(e_str)
		to.write(s)
		to.write("\n")
	for line in pi_lines[2:-2]:
		arr = line.split()
		if not (float(vertices[int(arr[1])-1][0]) <= 10**-15 and float(vertices[int(arr[2])-1][0]) <= 10**-15):
			s = " ".join(arr[1:3])
			po.write(s)
			po.write("\n")
	pi.close()
	vo.close()
	to.close()
	po.close()
