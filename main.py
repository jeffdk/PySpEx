from numpy import *
#from Tkinter import *
import Tkinter as Tk
from dataFunction import *

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


#Never got these things working 

# import os

# remoteServerDirectory="zwicky:/home/jeff/tovStars/stdTOVFiles/tenMSsy*"

# localDir="/home/jeff/work/PySEx"

# os.execl('/bin/pwd', './')

# os.mkdir('data')
 
# os.execl('/usr/bin/rsync', '-e', ,'-r', '--rsh=ssh', "zwicky:/home/jeff/tovStars/stdTOVFiles/tenMSsy*/Run/DensestPoint.dat", 'data/')

# print os.listdir("zwicky:/home/jeff/tovStars/stdTOVFiles/tenMSsy*/Run/")


dirBase='data/'
fileList=['DensestPoint']
Prefix='Lev'
Levs=[0,1,2]




exit

class DatFile():
    def parseFile(self, filename):
        infile = open(filename, 'r')
        if(infile):
            print infile.read()
        else:
            print "ERROR BAD FILENAME"
            exit
            
        
            
    def __init__(self):
        print "inittetd"
    


file = DatFile()

#file.parseFile("DensestPoint.dat")

octsym=dataFunction()
octfile=open("DensestPoint.dat",'r')
octsym.readFuncDataFromFile(octfile,[0],[4])
mpl.plot(octsym.data,octsym.points)
mpl.legend(["Densest Point"])
mpl.title("Densests Point")
mpl.grid()
#mpl.show()

#root = Tk.Tk()
#root.wm_title("Embedding in TK")
#root.bind("<Destroy>", destroy)


f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)
t = arange(0.0,3.0,0.01)
s = sin(2*pi*t)

b = f.add_subplot(212)
b.plot(octsym.data,octsym.points)
b.legend(["Densest Point"])
#b.title("Densests Point")
b.grid()

a.plot(t,s)




class Application(Tk.Frame):
    def say_hi(self):
        print "hi there, everyone!"

    def createWidgets(self):
        # a tk.DrawingArea
        canvas = FigureCanvasTkAgg(f, master=root)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        
        toolbar = NavigationToolbar2TkAgg( canvas, root )
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
      
        self.QUIT = Tk.Button(self)
        self.QUIT["text"] = "QUIT"
        self.QUIT["fg"]   = "red"
        self.QUIT["command"] =  self.quit

        self.QUIT.pack({"side": "left"})

        self.hi_there = Tk.Button(self)
        self.hi_there["text"] = "Hello",
        self.hi_there["command"] = self.say_hi

        self.hi_there.pack({"side": "left"})

    def __init__(self, master=None):
        Tk.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

root = Tk.Tk()
app = Application(master=root)
app.mainloop()
root.destroy()



print "HEWWO"


