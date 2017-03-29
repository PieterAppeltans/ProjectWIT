import numpy as np
vertices = [[0,0],[0,1],[1,0]]
area = abs(vertices[0][0]*(vertices[1][1] -vertices[2][1])+vertices[1][0]*(vertices[2][1]-vertices[0][1])+vertices[2][0]*(vertices[0][1]-vertices[1][1] ))

print area

G = np.array([vertices[1][1]

A = 1./(2.*area)*(vertices[0][0]+vertices[1][0]+vertices[2][0])

