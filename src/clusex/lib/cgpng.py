import warnings
import pandas as pd
import tkinter as tk
import tkinter.font as tkFont
import matplotlib.pyplot as plt
from io import BytesIO
from tkinter import Tk
from tkinter import messagebox as tkMessageBox
from PIL import ImageTk, Image
warnings.filterwarnings("ignore")



def cgPNG(pngfile_in, pngfile_out):
    #names of input and output files needed to run the GUI
    with open(pngfile_in) as f:
        pnglist = f.readlines()
    
    png_paths = [s.replace("\n", "") for s in pnglist]
    
    print("Ready to classify")

    root = Tk()
    gui = GUI(root, png_paths, pnglist, pngfile_out)
    root.mainloop()

class GUI:
    def __init__(self, master, png_paths, pnglist, pngfile_out):
        self.master = master

        self.index = 0
        self.png_paths = png_paths
        self.n_paths = len(png_paths)

        self.pngfile_out = pngfile_out

        self.table  = pd.DataFrame({'Name': pnglist})
        self.obj_id = []

        master.title("Classify Galaxies with Python and PNGs")
        
        #setting title
        #setting window size
        width=543
        height=520
        screenwidth = master.winfo_screenwidth()
        screenheight = master.winfo_screenheight()
        alignstr = '%dx%d+%d+%d' % (width, height, (screenwidth - width) / 2, (screenheight - height) / 2)
        master.geometry(alignstr)
        master.resizable(width=False, height=False)
        
        global ImgTitle
        ImgTitle=tk.Label(master)
        ft = tkFont.Font(family='Arial',size=10)
        ImgTitle["font"] = ft
        ImgTitle["fg"] = "#333333"
        ImgTitle["justify"] = "left"
        ImgTitle["text"] = "%s" % (png_paths[0])
        ImgTitle.place(x=10,y=10,width=200,height=25)
        
        image_gal = plt.imread(png_paths[0])
        plt.figure()
        #plt.title('%s' % (png_paths[0]))
        plt.axis('off')
        plt.imshow(image_gal) 
        img_data = BytesIO()
        plt.tight_layout()
        plt.savefig(img_data)
        
        load = Image.open(img_data)
        resized = load.resize((520, 380),Image.ANTIALIAS)
        render = ImageTk.PhotoImage(resized)
        Image_Plot = tk.Label(master, image = render)
        Image_Plot.image = render 
        Image_Plot.place(x=10,y=35,width=520,height=380)
    
    
        Class_label=tk.Label(master)
        ft = tkFont.Font(family='Arial',size=10)
        Class_label["font"] = ft
        Class_label["fg"] = "#333333"
        Class_label["justify"] = "left"
        Class_label["text"] = "Select a class:"
        Class_label.place(x=10,y=415,width=100,height=25)  
      
        Button_E=tk.Button(master)
        Button_E["bg"] = "#efefef"
        ft = tkFont.Font(family='Arial',size=10)
        Button_E["font"] = ft
        Button_E["fg"] = "#000000"
        Button_E["justify"] = "center"
        Button_E["text"] = "E"
        Button_E.place(x=10,y=440,width=70,height=25)
        Button_E["command"] = lambda:[self.E()]

        Button_S0=tk.Button(master)
        Button_S0["bg"] = "#efefef"
        ft = tkFont.Font(family='Arial',size=10)
        Button_S0["font"] = ft
        Button_S0["fg"] = "#000000"
        Button_S0["justify"] = "center"
        Button_S0["text"] = "S0"
        Button_S0.place(x=100,y=440,width=70,height=25)
        Button_S0["command"] = self.S0

        Button_S0B=tk.Button(master)
        Button_S0B["bg"] = "#efefef"
        ft = tkFont.Font(family='Arial',size=10)
        Button_S0B["font"] = ft
        Button_S0B["fg"] = "#000000"
        Button_S0B["justify"] = "center"
        Button_S0B["text"] = "S0B"
        Button_S0B.place(x=190,y=440,width=70,height=25)
        Button_S0B["command"] = self.S0B

        Button_S=tk.Button(master)
        Button_S["bg"] = "#efefef"
        ft = tkFont.Font(family='Arial',size=10)
        Button_S["font"] = ft
        Button_S["fg"] = "#000000"
        Button_S["justify"] = "center"
        Button_S["text"] = "S"
        Button_S.place(x=280,y=440,width=70,height=25)
        Button_S["command"] = self.S

        Button_SB=tk.Button(master)
        Button_SB["bg"] = "#efefef"
        ft = tkFont.Font(family='Arial',size=10)
        Button_SB["font"] = ft
        Button_SB["fg"] = "#000000"
        Button_SB["justify"] = "center"
        Button_SB["text"] = "SB"
        Button_SB.place(x=370,y=440,width=70,height=25)
        Button_SB["command"] = self.SB
        
        Button_Other=tk.Button(master)
        Button_Other["bg"] = "#efefef"
        ft = tkFont.Font(family='Arial',size=10)
        Button_Other["font"] = ft
        Button_Other["fg"] = "#000000"
        Button_Other["justify"] = "center"
        Button_Other["text"] = "Other"
        Button_Other.place(x=460,y=440, width=70,height=25)
        Button_Other["command"] = self.Other
        
        Button_Quit=tk.Button(master)
        Button_Quit["bg"] = "#efefef"
        ft = tkFont.Font(family='Arial',size=10)
        Button_Quit["font"] = ft
        Button_Quit["fg"] = "#000000"
        Button_Quit["justify"] = "center"
        Button_Quit["text"] = "Quit"
        Button_Quit.place(x=460,y=480,width=70,height=25)
        Button_Quit["command"] = self.ask_quit
        
    def show_next_image(self):
        '''displays next image, or saves data when classification is over'''

        if len(self.obj_id) < self.n_paths: # opens next image
            self.index += 1
            
            # Changes file-name's label's text
            ImgTitle["text"] = "%s" % (self.png_paths[self.index])
            
            # Plots the next image, same code 
            image_gal = plt.imread(self.png_paths[self.index])
            plt.figure()
            #plt.title('%s' % (self.png_paths[self.index]))
            plt.axis('off')
            plt.imshow(image_gal) 
            img_data = BytesIO()
            plt.tight_layout()
            plt.savefig(img_data)
            
            load = Image.open(img_data)
            resized = load.resize((520, 380),Image.ANTIALIAS)
            render = ImageTk.PhotoImage(resized)
            Image_Plot = tk.Label(self.master, image = render)
            Image_Plot.image = render 
            Image_Plot.place(x=10,y=35,width=520,height=380)
            
        if len(self.obj_id) == self.n_paths: # saves data into file

            # creates pandas dataframe with file names column and class selections
            df = pd.DataFrame({'Name': self.png_paths,'Class': self.obj_id})

            # saves newly merged dataframe
            df.to_csv(path_or_buf=self.pngfile_out, sep='\t', index=False, header=False, mode='w')    
        
            print('Classification successful.')
            print('Data has been saved to %s.' %(self.pngfile_out))
            tkMessageBox.showinfo('Classification successful', 'Data has been saved to %s.' %(self.pngfile_out))

        if len(self.obj_id) > (self.n_paths):
            print('Everything is already classified.')
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def E(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('E')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def S0(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('S0')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def S0B(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('S0B')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def S(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('S')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def SB(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('SB')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")
    def Other(self):
        if len(self.obj_id) < self.n_paths:
            self.obj_id.append('Other')
            self.show_next_image()
        else:
            tkMessageBox.showinfo("Classification successful", "Everything is already classified.")

    def ask_quit(self):
        if tkMessageBox.askokcancel("Quit", "Are you sure you want to quit?\nData will be lost if you're not finished."):
            self.master.destroy()
            
