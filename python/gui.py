import matplotlib
import numpy as np
import scipy.interpolate
from parse import *
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import sys
import Tkinter as Tk

def destroy(e):
    sys.exit()

root = Tk.Tk()
root.wm_title("The Pear Project")

vertices,elements = parse_input("circle.1")
u,v = parse_u_v()

xmin = min(vertices[:,0])
xmax = max(vertices[:,0])
ymin = min(vertices[:,1])
ymax = max(vertices[:,1])


xlin = np.linspace(xmin,xmax,300)
ylin = np.linspace(ymin,ymax,300)

X,Y = np.meshgrid(xlin,ylin)
U = scipy.interpolate.griddata(vertices,u,(X,Y),'linear');
V = scipy.interpolate.griddata(vertices,v,(X,Y),'linear');


f = Figure(figsize=(5, 4), dpi=100)
a_u = f.add_subplot(121)

a_u.contourf(X,Y,U,10)
a_u.set_title('Contour u(r,z)')
a_u.set_xlabel('r(m)')
a_u.set_ylabel('z(m)')

a_v = f.add_subplot(122)

a_v.contourf(X,Y,V,10)
a_v.set_title('Contour v(r,z)')
a_v.set_xlabel('r(m)')
a_v.set_ylabel('z(m)')


# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

button = Tk.Button(master=root, text='Quit', command=sys.exit)
button.pack(side=Tk.BOTTOM)

Tk.mainloop()
