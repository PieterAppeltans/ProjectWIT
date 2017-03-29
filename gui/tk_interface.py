from Tkinter import *
import subprocess
import tkMessageBox
import numpy as np
import time
import scipy.interpolate
from parse import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

simulation_options = {"Custom preset":[0,0,0],"Orchard":[25,20.8,0.04],"Shelf life":[20,20.8,0],
    "Refrigerator":[7,20.8,0],"Precooling":[-1,20.8,0],"Disorder inducing":[-1,2,5],"Optimal CA":[-1,2,0.7]}
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

        self.subtitle = Label(self,text="Simulate for various conditions, with various code versions",font="Helvetica 10")
        self.subtitle.pack({"side": "top"})

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()


def next1(temp,nu,nv):
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


def next2(file_,area,angle,compile,version):
    global AREA,ANGLE,FILE,TEMP,NU,NV,plot
    try:
        AREA = float(area)
        ANGLE = float(angle)
        FILE = file_
        mesh_field.grid_forget()
        mesh_field.destroy()
    except:
        tkMessageBox.showerror("Error", "Please enter a floating point number")
    mesh_plot = False
    try:
        print "=== CREATING MESH ====\n"
        subprocess.call(["bash","create_mesh.sh",str(AREA),str(ANGLE),FILE],cwd=None)
        vertices,elements = parse_input(FILE+".1")
        #mesh_plot = MeshPlot(vertices,elements,master=root)
        #root.update()
    except:
        tkMessageBox.showerror("Error", "An error occured during mesh generation")
    else:
        if compile:
            if version == "dense":
                subprocess.call(["bash","compile_cpp_dense.sh"],cwd=None)
                print "\n=== DENSE VERSION COMPILED ==="
            elif version == "sparse" or version == "sparse quasi":
                subprocess.call(["bash","compile_cpp_sparse.sh"],cwd=None)
                print "\n=== SPARSE VERSION COMPILED ==="
        try:
            if version == "dense":
                print "\n=== EXECUTING DENSE VERSION CPP ==="
                subprocess.call(["bash","execute_cpp_dense.sh",FILE,str(TEMP),str(NU),str(NV)],cwd=None)
            elif version == "sparse":
                print "\n=== EXECUTING SPARSE VERSION CPP ==="
                subprocess.call(["bash","execute_cpp_sparse.sh",FILE,str(TEMP),str(NU),str(NV),str(0)],cwd=None)
            elif version == "sparse quasi":
                print "\n=== EXECUTING SPARSE VERSION CPP WITH QUASI NEWTON ==="
                subprocess.call(["bash","execute_cpp_sparse.sh",FILE,str(TEMP),str(NU),str(NV),str(1)],cwd=None)
            elif version == "matlab":
                print "\n=== EXECUTING MATLAB (DENSE) ==="
                write_matlab(FILE,vertices,elements)
                subprocess.call(["matlab","-nojvm -nodisplay",'-r "GUI_FEM(%f,%f,%f);exit"'%(TEMP,NU,NV)],cwd="../matlab")
        except:
            tkMessageBox.showerror("Error","An error occured during calculation")
    if mesh_plot:
        mesh_plot.destroy()
    plot = ResultPlot(vertices,elements,master=root)

def next3():
        global plot,input_field
        plot.destroy()
        input_field = InputField(master=root)


class LoadingScreen(Frame):
    def __init__(self,master=None):
        Frame.__init__(self, master)
        self.pack()
        self.create_text()

    def create_text(self):
        self.loading_text = Label(self,text="Waiting on results (May take 1-5 minutes)",font="Helvetica 10")
        self.loading_text.pack()
        print "Creating text"

class MeshPlot(Frame):
    def __init__(self,vertices,elements,master=None):
        Frame.__init__(self,master)
        self.pack()
        self.create_mesh(vertices,elements)
    def create_mesh(self,vertices,elements):
        f = Figure()
        a = f.add_subplot(111)
        for elem in elements:
            a.plot([vertices[int(elem[0])][0],vertices[int(elem[1])][0],vertices[int(elem[2])][0],vertices[int(elem[0])][0]],
                [vertices[int(elem[0])][1],vertices[int(elem[1])][1],vertices[int(elem[2])][1],vertices[int(elem[0])][1]],'k')
        a.set_title('Mesh')
        a.set_xlabel('r(m)')
        a.set_ylabel('z(m)')
        canvas = FigureCanvasTkAgg(f, master=self)
        canvas.show()
        canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
        canvas._tkcanvas.pack(side="top", fill="both", expand=1)

class ResultPlot(Frame):
        global FILE
        def __init__(self,vertices,elements, master=None):
            Frame.__init__(self,master)
            self.pack(fill=BOTH, expand=YES)
            self.create_plot(vertices,elements)
        def create_plot(self,vertices,elements):
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


            f = Figure()
            a_u = f.add_subplot(121)

            a_u_cont = a_u.contourf(X,Y,U,10)
            a_u.set_title('Contour u(r,z)')
            a_u.set_xlabel('r(m)')
            a_u.set_ylabel('z(m)')

            cbar = f.colorbar(a_u_cont)

            a_v = f.add_subplot(122)

            a_v_cont = a_v.contourf(X,Y,V,10)
            a_v.set_title('Contour v(r,z)')
            a_v.set_xlabel('r(m)')
            a_v.set_ylabel('z(m)')

            cbar = f.colorbar(a_v_cont)

            # a tk.DrawingArea
            canvas = FigureCanvasTkAgg(f, master=self)
            canvas.show()
            canvas.get_tk_widget().pack(side="top", fill="both", expand=1)

            canvas._tkcanvas.pack(side="top", fill="both", expand=1)
            self.button = Button(self, text='Back to start',command=lambda: next3())
            self.button.pack({"side": "bottom"})
            print "\n=== DONE ==="


class MeshField(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        self.label_version = Label(self, text="Choose cpp version:")
        self.label_version.grid(row=0)
        self.version = StringVar(self)
        self.version.set("sparse quasi") # default value

        self.version_menu = OptionMenu(self, self.version,"sparse quasi","sparse","dense","matlab")
        self.version_menu.grid(row=0,column=1)


        self.label_option = Label(self, text="Choose grid:")
        self.label_option.grid(row=1)
        self.mesh = StringVar(self)
        self.mesh.set("pear") # default value

        self.mesh_menu = OptionMenu(self, self.mesh,"pear", "circle")
        self.mesh_menu.grid(row=1,column=1)

        self.label_area = Label(self, text="Area:")
        self.label_angle = Label(self, text="Angle:")

        self.label_area.grid(row=2)
        self.label_angle.grid(row=3)

        self.compile = BooleanVar()
        self.compile_menu = Checkbutton(self, text="Compile c++", variable=self.compile)
        self.compile_menu.grid(row=4,column=1)

        self.area = DoubleVar()
        self.angle = DoubleVar()

        self.area.set(5)
        self.angle.set(30)

        self.entry_area = Entry(self,textvariable=self.area)
        self.entry_angle = Entry(self,textvariable=self.angle)

        self.entry_area.grid(row=2,column=1)
        self.entry_angle.grid(row=3,column=1)

        self.button = Button(self, text='Next',command=lambda: next2(self.mesh.get(),self.area.get(),
            self.angle.get(),self.compile.get(),self.version.get()))
        self.button.grid(row=5)


class InputField(Frame):

    global simulation_options

    def preset_changed(self,arg):
        self.temp.set(simulation_options.get(arg)[0])
        self.nu.set(simulation_options.get(arg)[1])
        self.nv.set(simulation_options.get(arg)[2])

    def createWidgets(self):
        self.label_option = Label(self, text="Choose preset:")
        self.label_option.grid(row=0)
        self.option = StringVar(self)
        self.option.set("Custom preset") # default value

        self.option_menu = OptionMenu(self, self.option,*sorted(simulation_options.keys()),command=self.preset_changed)
        self.option_menu.grid(row=0,column=1)

        self.label_temp = Label(self, text="Temprature(C):")
        self.label_nu = Label(self, text="Ambient oxigen concentration(%):")
        self.label_nv = Label(self, text="Ambient oxigen concentration(%):")

        self.label_temp.grid(row=1)
        self.label_nu.grid(row=2)
        self.label_nv.grid(row=3)

        self.temp = DoubleVar()
        self.nu = DoubleVar()
        self.nv = DoubleVar()

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
root.geometry("800x600")
header = Header(master=root)
input_field = InputField(master=root)
footer = Footer(master=root)
root.mainloop()
