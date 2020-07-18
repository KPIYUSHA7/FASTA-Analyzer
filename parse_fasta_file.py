#Import the libraries

import pandas as pd
import tkinter as tk
from tkinter import filedialog
from tkinter.ttk import *

#Tkinter GUI

window = tk.Tk()
window.geometry("890x400")
window.configure(bg = "gray26")

frame1 = tk.Label(window, text = "FASTA SEQUENCE ANALYZER", bg = 'gray26', font = ('courier new', 20, 'bold'), fg = 'azure')
frame1.grid(row=0, column=1)

frame2 = tk.Frame(window)
frame2.grid(row=6, column=1)

frame3 = tk.Frame(window)
frame3.grid(row=9, column=1)

frame4 = tk.Frame(window)
frame4.grid(row=15,column=1)

frame5 = tk.Frame(window)
frame5.grid(row=12,column=1)

frame6 = tk.Frame(window)
frame6.grid(row=14, column=1)

l1 = tk.Label(window, text = "Upload the file (fasta format): ", fg = "azure", bg = 'gray26')
l1.configure(font = ('times new roman', 14))
l1.grid(row=3, column=0)

e1 = tk.Entry(window)
e1.grid(row=3, column=1)
e1.configure(width=65)

def browsefunc():  
    global info
    filename = filedialog.askopenfilename()
    e1.insert(tk.END, filename)
    f = open(filename, 'r')
    info = f.read()

b1 = tk.Button(window, text = "Browse", command=browsefunc, bg = "steelblue")
b1.configure(font =  ('courier new', 12))
b1.grid(row=3, column=2)

#------------------------------------------------------------------------------------------------------------------------------

#General Sequence Information

def general_statistics():
       
    without_header = info.split("\n", 1)[1]
    without_header = without_header.replace("\n", "")
       
    frame3.grid(row=9, column=1)

    counter_A = 0
    counter_G = 0
    counter_T = 0
    counter_C = 0
    counter_N = 0
    Unknown = 0

    for i in without_header:
        if(i == "A"):
            counter_A += 1
        elif(i == "G"):
            counter_G += 1
        elif(i == "T"):
            counter_T += 1
        elif(i == "C"):
            counter_C += 1
        elif(i == "N"):
            counter_N +=1
        else:
            Unknown += 1


    GC_content = ((counter_G + counter_C) / (counter_A + counter_G + counter_C + counter_T))  * 100
    
    gs1 = tk.Label(frame3, text = "Adenine Content: {} bp \n".format(counter_A))
    gs1.grid(row = 12, column = 1)
    gs1.configure(width=40, font = ('times new roman', 11))
    gs2 = tk.Label(frame3, text = "Guanine Content: {} bp \n".format(counter_G))
    gs2.configure(width=40, font = ('times new roman', 11))
    gs2.grid(row = 13, column = 1)
    gs3 = tk.Label(frame3, text = "Thymine Content: {} bp \n".format(counter_T))
    gs3.configure(width=40, font = ('times new roman', 11))
    gs3.grid(row = 14, column = 1)
    gs4 = tk.Label(frame3, text = "Cytosine Content: {} bp \n".format(counter_C))
    gs4.configure(width=40, font = ('times new roman', 11))
    gs4.grid(row = 15, column = 1)
    gs5 = tk.Label(frame3, text = "GC content:  {} %".format(GC_content))
    gs5.configure(width=40, font = ('times new roman', 11))
    gs5.grid(row = 16, column = 1)
      
    frame4.grid_forget()
    frame5.grid_forget()
    frame6.grid_forget()
    
b2 = tk.Button(frame2, text = "Show General Statistics", command = general_statistics, bg = 'steelblue')
b2.grid(row=8, column=0)
b2.configure(width=25, font =  ('courier new', 11))

#------------------------------------------------------------------------------------------------------------------------------

#Transcription Sequence - DNA to RNA Conversion

def transcription():
        
    frame4.grid(row=15,column=1)
    
       
    list_of_characters = []
    for i in info:
        list_of_characters.append(i)
        
    transcribed_RNA_sequence = ["U" if x == "T" else x for x in list_of_characters]
    s = ""
    transcribed = s.join(transcribed_RNA_sequence)
    write_file = open("DNA_to_RNA_transcription.txt", "wt")
    n = write_file.write(transcribed)
    write_file.close()
    result = tk.Label(frame4, text = "Completed! \n Transcribed Sequence is stored in the working directory as \n 'DNA_to_RNA_transcription.txt'")
    result.grid(row=12,column=1)
    result.configure(font = ("times new roman", 14))
       
    frame3.grid_forget()
    frame5.grid_forget()
    frame6.grid_forget()
    
    
b3 = tk.Button(frame2, text = "Perform Transcription", command = transcription, bg = 'steelblue')
b3.grid(row=8, column=1)
b3.configure(width=25, font =  ('courier new', 11))

#------------------------------------------------------------------------------------------------------------------------------

def reverse_complement():
         
    frame6.grid(row=14, column=1)
    
    without_header = info.split("\n", 1)[1]
    without_header = without_header.replace("\n", "")
  
    list_of_characters = []
    for i in without_header:
        if (i=="A"):
            list_of_characters.append("T")
        elif(i=="T"):
            list_of_characters.append("A")
        elif(i=="C"):
            list_of_characters.append("G")
        elif(i=="G"):
            list_of_characters.append("C")
        else:
            pass

    list_of_characters.reverse()
    s = ''
    complement = s.join(list_of_characters)
    file = open("Reverse_complement.txt", "w")
    file.writelines(list_of_characters)
    file.close()   
    rc = tk.Label(frame6, text="Completed! \n Saved in the working directory as 'Reverse_complement.txt'")
    rc.grid(row=8,column=0)
    rc.configure(font = ("times new roman", 14))
        
    frame3.grid_forget()
    frame4.grid_forget()
    frame5.grid_forget()
    
b4 = tk.Button(frame2, text = "Reverse Complement", command = reverse_complement, bg = 'steelblue')
b4.grid(row=12, column=0)
b4.configure(width=25,  font = ('courier new', 11))
    

#------------------------------------------------------------------------------------------------------------------------------

#RNA to Protein conversion

def translation():
            
    frame5.grid(row=12,column=1)
    
    # Scrape data from ttp://www.hgmd.cf.ac.uk/docs/cd_amino.html for the codon table using pandas
    
    codon_table = pd.read_html("http://www.hgmd.cf.ac.uk/docs/cd_amino.html")
    codon_table = codon_table[0]
    without_header = info.split("\n", 1)[1]
    without_header = without_header.replace("\n", "")
      
    n = 3
    amino_acids_list = [without_header[i:i+n] for i in range(0, len(without_header), n)]
    aa = tk.Label(frame5, text = "Completed! \n The total number of amino acids in this genome is: {}".format(len(amino_acids_list)))
    aa.grid(row=12,column=1)
    aa.configure(font = ("times new roman", 14))

    aa_present = []
    for i in amino_acids_list:
        if i in list(codon_table.iloc[:,0]):
            aa_present.append(codon_table[codon_table[0]== i])
    df = (pd.concat(aa_present))
    
    df.to_csv("AA_present.csv")
    save_file = tk.Label(frame5, text = "CSV file stored in working directory as AA_present.csv")
    save_file.grid(row=13,column=1)
    save_file.configure(font = ("times new roman", 14))
    
    frame3.grid_forget()
    frame4.grid_forget()
    frame6.grid_forget()
    
b5 = tk.Button(frame2, text = "Translation", command = translation, bg = 'steelblue')
b5.grid(row=12, column=1)
b5.configure(width=25, font =  ('courier new', 11))

#------------------------------------------------------------------------------------------------------------------------------

def clear_output():
    
    #Remove all frames
    
    frame3.grid_forget()
    frame4.grid_forget()
    frame5.grid_forget()
    frame6.grid_forget()


b6 = tk.Button(window, text = "Clear Output", command=clear_output, bg = 'steelblue')
b6.grid(row=8, column=2)
b6.configure(font =  ('courier new', 11))

pl = tk.Label(window, text = "** Pls press 'Clear Output' \n after each analysis **", bg = 'gray26', fg = 'red')
pl.configure(font = ('times new roman', 13))
pl.grid(row=4, column=0)


window.mainloop()

#------------------------------------------------------------------------------------------------------------------------------
