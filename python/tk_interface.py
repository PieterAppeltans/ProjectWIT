from Tkinter import *

class Header(Frame):    

    def createWidgets(self):
        self.title = Label(self,text="Pear Project",font="Helvetica 16 bold")
        self.title.pack({"side": "top"})

        self.subtitle = Label(self,text="Subtitle",font="Helvetica 10")
        self.subtitle.pack({"side": "top"})
        

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()
class InputField(Frame):
    
    def createWidgets(self):
        self.label_temp = Label(self, text="Temprature")
        self.label_nu = Label(self, text="Ambient oxigen concentration")
        self.label_nv = Label(self, text="Ambient oxigen concentration")
        
        self.label_temp.grid(row=0)
        self.label_nu.grid(row=1)
        self.label_nv.grid(row=2)
        
        self.entry_temp = Entry(self)
        self.entry_nu = Entry(self)
        self.entry_nv = Entry(self)

        self.entry_temp.grid(row=0,column=1)
        self.entry_nu.grid(row=1,column=1)
        self.entry_nv.grid(row=2,column=1)

        self.button = Button(self, text='Next')
        self.button.grid(row=3)

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

class Footer(Frame):    

    def createWidgets(self):

        self.footertext = Label(self,text="Footer",font="Helvetica 10")
        self.footertext.pack({"side": "bottom"})
        

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
root.destroy()
