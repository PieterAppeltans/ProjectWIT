from Tkinter import *
import subprocess
import tkMessageBox
import numpy as np
import scipy.interpolate
from parse import *
import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
simulation_options = {"Orchard":[25,20.8,0.04],"Shelf life":[20,20.8,0],"Refrigerator":[7,20.8,0]}
TEMP = 0.
NU = 0.
NV = 0.
AREA = 0.
ANGLE = 0.
FILE_ = ""
class Header(Frame):

    def createWidgets(self):
        self.title = Label(self,text="Pear Project",font="Helvetica 16 bold")
        self.title.pack({"side": "top"})

        self.subtitle = Label(self,text="Set simulation conditions",font="Helvetica 10")
        self.subtitle.pack({"side": "top"})


    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

def next1(temp,nu,nv):
    print "Next 1"
    global TEMP,NU,NV,input_field,mesh_field,root
    try:
        TEMP = float(temp)
        NU = float(nu)
        NV = float(nv)
        input_field.grid_forget()
        input_field.destroy()
        mesh_field = MeshField(master=root)
    except:
        tkMessageBox.showerror("Error", "Please enter a floating point number")
def next2(file_,area,angle):
    global AREA,ANGLE,FILE
    try:
        AREA = float(area)
        ANGLE = float(angle)
        FILE = file_
        mesh_field.grid_forget()
        mesh_field.destroy()
    except:
        tkMessageBox.showerror("Error", "Please enter a floating point number")
    loading_screen = LoadingScreen(master=root)
    #subprocess.call("pwd")
    try:
        subprocess.call(["bash","run.sh",str(AREA),str(ANGLE),FILE],cwd="./")
    except:
        tkMessageBox.showerror("Error", "An error occured during mesh generation")
    try:
        subprocess.call(["bash","run2.sh",FILE],cwd="./")
    except:
        tkMessageBox.showerror("Error","An error occured during calculation")
    print "Done"
    loading_screen.destroy()
    plot = ResultPlot(master=root)
    #subprocess.call(["exec triangle -p -a{} -q{} ../triangle/{}".format(AREA,ANGLE,FILE)],cwd="../triangle")
    #subprocess.call("cd ../code/; ./test_eigen.o")
class LoadingScreen(Frame):
    def __init__(self,master=None):
        Frame.__init__(self, master)
        self.pack()
        self.create_text()

    def create_text(self):
        self.loading_text = Label(self,text="Waiting on results (May take 1-5 minutes)",font="Helvetica 10")
        self.loading_text.pack()
        print "Creating text"
class ResultPlot(Frame):
        global FILE
        def __init__(self, master=None):
            Frame.__init__(self, master)
            self.pack()
            self.create_plot()
        def create_plot(self):
            vertices,elements = parse_input(FILE+".1")
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
            canvas = FigureCanvasTkAgg(f, master=self)
            canvas.show()
            canvas.get_tk_widget().pack(side="top", fill="both", expand=1)

            canvas._tkcanvas.pack(side="top", fill="both", expand=1)
            print "done"


class MeshField(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()
    def createWidgets(self):
        self.label_option = Label(self, text="Choose grid:")
        self.label_option.grid(row=0)
        self.mesh = StringVar(self)
        self.mesh.set("pear") # default value

        self.mesh_menu = OptionMenu(self, self.mesh,"pear", "circle")
        self.mesh_menu.grid(row=0,column=1)

        self.label_area = Label(self, text="Area:")
        self.label_angle = Label(self, text="Angle:")

        self.label_area.grid(row=1)
        self.label_angle.grid(row=2)

        self.area = StringVar()
        self.angle = StringVar()

        self.area.set("0.00005")
        self.angle.set("30")

        self.entry_area = Entry(self,textvariable=self.area)
        self.entry_angle = Entry(self,textvariable=self.angle)

        self.entry_area.grid(row=1,column=1)
        self.entry_angle.grid(row=2,column=1)

        self.button = Button(self, text='Next',command=lambda: next2(self.mesh.get(),self.area.get(),self.angle.get()))
        self.button.grid(row=4)

class InputField(Frame):

    def preset_changed(self,arg):
        self.temp.set(self.simulation_options.get(arg)[0])
        self.nu.set(self.simulation_options.get(arg)[1])
        self.nv.set(self.simulation_options.get(arg)[2])

    def createWidgets(self):
        self.label_option = Label(self, text="Choose preset:")
        self.label_option.grid(row=0)
        self.option = StringVar(self)
        self.option.set("Custom preset") # default value

        self.option_menu = OptionMenu(self, self.option,"Custom preset","Orchard", "Shelf life", "Refrigerator",command=self.preset_changed)
        self.option_menu.grid(row=0,column=1)

        self.label_temp = Label(self, text="Temprature:")
        self.label_nu = Label(self, text="Ambient oxigen concentration:")
        self.label_nv = Label(self, text="Ambient oxigen concentration:")

        self.label_temp.grid(row=1)
        self.label_nu.grid(row=2)
        self.label_nv.grid(row=3)

        self.temp = StringVar()
        self.nu = StringVar()
        self.nv = StringVar()

        self.entry_temp = Entry(self,textvariable=self.temp)
        self.entry_nu = Entry(self,textvariable=self.nu)
        self.entry_nv = Entry(self,textvariable=self.nv)

        self.entry_temp.grid(row=1,column=1)
        self.entry_nu.grid(row=2,column=1)
        self.entry_nv.grid(row=3,column=1)

        self.button = Button(self, text='Next',command=lambda: next1(self.temp.get(),self.nu.get(),self.nv.get()))
        self.button.grid(row=4)

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.simulation_options = {"Custom preset":[0,0,0],"Orchard":[25,20.8,0.04],"Shelf life":[20,20.8,0],"Refrigerator":[7,20.8,0]}
        self.createWidgets()

class Footer(Frame):

    def createWidgets(self):
        self.footertext = Label(self,text="Pieter Appeltans & Lennart Bulteel",font="Helvetica 10")
        self.footertext.pack()


    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack({"side": "bottom"})
        self.createWidgets()

root = Tk()
root.wm_title("Pear project")
root.geometry("500x500")
header = Header(master=root)
input_field = InputField(master=root)
footer = Footer(master=root)
root.mainloop()
