import matplotlib
from parse import *
matplotlib.use('TkAgg')

from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import sys
import Tkinter as Tk

def destroy(e):
    sys.exit()

root = Tk.Tk()
root.wm_title("Embedding in TK")

vertices,elements = parse_input("circle.1")

f = Figure(figsize=(5, 4), dpi=100)
a = f.add_subplot(111)

for elem in elements:
    a.plot([vertices[elem[0]][0],vertices[elem[1]][0],vertices[elem[2]][0],vertices[elem[0]][0]], [vertices[elem[0]][1],vertices[elem[1]][1],vertices[elem[2]][1],vertices[elem[0]][1]],'k')
a.set_title('CirclElements')
a.set_xlabel('X axis label')
a.set_ylabel('Y label')


# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

button = Tk.Button(master=root, text='Quit', command=sys.exit)
button.pack(side=Tk.BOTTOM)

Tk.mainloop()
