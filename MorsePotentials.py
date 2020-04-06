import matplotlib.pyplot as plt
import numpy as np
import math 
import scipy
import sympy as sym
import Tkinter as tk
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mtick
#import pandas as pd

#Data Frames to be exported to Excel
#MorsePotentials = pd.DataFrame(columns = ["Internuclear Distance", "Potential Energy"])
#FCF = pd.DataFrame(columns = ["Wavenumbers (cm-1)", "Intensity"])
#lines = False

#Quantum Methods with Mathematica: https://books.google.com/books?id=tXIukmOFqscC&pg=PA189&lpg=PA189&dq=morse+oscillator+mathematica&source=bl&ots=74Ru4dMOLB&sig=b5Xg4UFTiakcIv36ik3AN2yYvAQ&hl=en&sa=X&ved=0ahUKEwiM3sStopHVAhUIzGMKHZ-vA9EQ6AEIRTAG#v=onepage&q=morse%20oscillator%20mathematica&f=false
#Simulating the Physical World: https://books.google.com/books?id=6pzEKsDkEVAC&pg=PA87&lpg=PA87&dq=morse+oscillator+python&source=bl&ots=0pIPCftJfM&sig=4kTEuoDVs_ZvzidW0zhMMS7Bs74&hl=en&sa=X&ved=0ahUKEwjpwv-jqJHVAhUN3mMKHbs3CzI4ChDoAQguMAI#v=onepage&q=morse%20oscillator%20python&f=false

"""Constants"""
pi = np.pi                  #Pi
c = 3e10                    #Speed of Light (cm/s)
h = 6.626e-34               #Planck's Constant
hBar = h / (2 * pi)         #Reduced Planck's Constant
m = 0.1269 / (6.022e23)     #Mass
u = m * m / (m + m)         #Reduced Mass
re1 = 3.025e-10             #Excited State Internuclear Distance 
re2 = 2.666e-10             #Ground State Internuclear Distance


"""Student input for Morse potential and wave function calculation"""
#Enter values for D'0 and D''0 in wavenumbers:
D02 = 4357.91054
D01 = 12302.85736

#Enter a value for V'e and V''e in wavenumbers:
Ve1 = 133.8434639
Ve2 = 213.36

#Enter a value for V'eX'e and V''eX''e in wavenumbers:
VeXe1 = 1.01219076
VeXe2 = 0.14

#Enter values for D'e and D''e in wavenumbers:
De1 = 4424.579224
De2 = 12409.50236

#Enter a value for Eel in wavenumbers:
Eel = 15587.92


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
################ DO NOT EDIT BELOW THIS LINE ##################
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#Calculation of force constant k:
#k1 = u*((2*pi*c*Ve1)**2)
#k2 = u*((2*pi*c*Ve2)**2)
k1 = 4 * De1 / (Ve1)
k2 = 4 * De2 / (Ve2)
     


#Calculation of alpha:
#alpha1 = np.sqrt(k1*u)/(hBar*6.022*10**23)
#alpha2 = np.sqrt(k2*u)/(hBar*6.022*10**23)
#alpha1 = np.sqrt(k1/(2*De1))
#alpha2 = np.sqrt(k2/(2*De2))
alpha1 = Ve1 * np.sqrt(pi * c * u/(hBar * De1))
alpha2 = Ve2 * np.sqrt(pi * c * u/(hBar * De2))

#Calculation of beta:
beta1 = np.sqrt(k1 / (2 * h * c * De1))
beta2 = np.sqrt(k2 / (2 * h * c * De2))

#Data ranges for x y plot of Morse Potential
xRange = np.arange(0, 7e-9, 0.000000000001)
yExcitedMorse = np.zeros(len(xRange))
yGroundMorse = np.zeros(len(xRange))

#Calculation of Potential Energy
for i in np.arange(len(xRange)):
    yExcitedMorse[i] = (Eel + De1 * (1 - np.exp(-alpha1 * (xRange[i] - re1)))**2)
    yGroundMorse[i] = (De2 * (1 - np.exp(-alpha2 * (xRange[i] - re2)))**2)





"""
end
Used when 'OK' is pressed on tkinter dialogue box
grabs user input values for ground state to excited 
state transition and calculates appropriate values
for plotting the wave functions
"""
def end():
    global master, var1, var2, v1, v2, z1, z2, b1, b2, Laguerre1, Laguerre2, Nn1, Nn2, x
    #obtain ground state and excited state vibrational levels from dialogue box
    v1 = int(var1.get())
    v2 = int(var2.get())

#    v1 = 3
#    v2 = 14

    #calculate ground state wave function variables z2, b2, Nn2 and Laguerre2
    x = np.arange (-20, 100, 0.25)
    z2=[]
    for i in np.arange(len(x)):
        z2.append(k2 * np.exp(-alpha2 * (x[i])))
    
    b2 = k2 - 2 * v2 - 1
    
    Nn2 = (alpha2 * b2 * float(math.factorial(v2)) / sym.gamma(k2-v2))**(0.5)
    
    Laguerre2 = np.zeros(len(z2))
    for i in np.arange(len(z2)): 
        Laguerre2[i] = sym.assoc_laguerre(v2,b2,z2[i])
        
#    wfn2 = np.zeros(len(x))
#    for j in np.arange(len(z2)):    
#        wfn2[j] = Nn2 * np.exp(-z2[j]/2) * z2[j] ** (b2/2) * Laguerre2[j]

    #calculate excited state wave function variables z1, b1, Nn1 and Laguerre1
    z1=[]
    for i in np.arange(len(x)):
        z1.append(k1 * np.exp(-alpha1 * (x[i])))
    
    b1 = k1 - 2 * v1 - 1
    
    Nn1 = (alpha1 * b1 * float(math.factorial(v1)) / sym.gamma(k1-v1))**(0.5)
    
    Laguerre1 = np.zeros(len(z1))
    for i in np.arange(len(z1)): 
        Laguerre1[i] = sym.assoc_laguerre(v1,b1,z1[i])
    
#    wfn1 = np.zeros(len(x))
#    for j in np.arange(len(z1)):    
#        wfn1[j] = Nn1 * np.exp(-z1[j]/2) * z1[j] ** (b1/2) * Laguerre1[j]
   
#    fig,ax = plt.subplots()
#    plt.plot(x,wfn2,'k')     
#    plt.plot(x,wfn1,'r')

    #call plot function and kill dialogue window
    plot()    
    master.quit()
    master.destroy()
    
y1=0
y2=0

"""
plot
Handles plotting of Morse potentials, vibrational states and wave functions
"""
def plot():
    global n,s,x,Eel,Ve1,Ve2,alpha1,alpha2,xRange,yExcitedMorse,yGroundMorse,VeXe1,VeXe2,re1,re2,y1,y2,wfn2,wfn1,z1,b1,Laguerre1,z2,b2,Laguerre2
    
    #create a new figure to display the Morse potentials and vibrational states
    fig,ax = plt.subplots()
    plt.title("Morse Potentials with Vibrational States")
    plt.xlabel("Internuclear Distance")
    plt.ylabel("Potential Energy")
    ax.set_ylim(0, 30000)
    ax.set_xlim(1e-10, 7e-10)
    ax.axvline(2.666e-10)
    ax.axvline(3.025e-10)
     
    #draw vibrational states on figure and calculate wave function
    for i in np.arange(60):
        #draw excited vibrational states
        if (Eel + Ve1 * i - (i * VeXe1 * (i + 1))) <= yExcitedMorse[-1] and (Eel + Ve2 * i - (i * VeXe2 * (i + 1))) >= min(yExcitedMorse):
            ax.axhline((Eel + Ve1 * i - (i * VeXe1 * (i + 1))))
        #draw ground vibrational states
        if i != 0 and (Ve2 * i - (i * VeXe2 * (i + 1))) <= yGroundMorse[-1]:
            ax.axhline((Ve2 * i - (i * VeXe2 * (i + 1))))
        #calculate ground state wave function
        if i == v2+1:
            ax.axhline((Ve2 * i - (i * VeXe2 * (i + 1))),color = 'y')   

            wfn2 = np.zeros(len(x))
            for j in np.arange(len(z2)):    
                wfn2[j] = Nn2 * np.exp((-z2[j]) / 2) * (z2[j]) ** (b2/2) * Laguerre2[j] + (Ve2 * i - (i * VeXe2 * (i + 1)))
        #calculate excited state wave function
        if i == v1+1:
            ax.axhline((Eel + Ve1 * i-(i * VeXe1 * (i + 1))),color = 'y')
            
            wfn1 = np.zeros(len(x))
            for j in np.arange(len(z1)):
                wfn1[j] = Nn1 * np.exp(-z1[j] / 2) * z1[j] ** (b1 / 2) * Laguerre1[j] + (Eel + Ve1 * i - (i * VeXe1 * (i + 1)))
        
        """Scaled(?) anharmonic wave functions""" 
#        if i == v2+1:
#            ax.axhline((Ve2*i-(i*VeXe2*(i+1))),color = 'y')   
#
#            wfn2 = np.zeros(len(x))
#            for j in np.arange(len(z2)):    
#                wfn2[j] = Nn2 * np.exp((-z2[j])/2) * (z2[j]) ** (b2/2) * Laguerre2[j] + (Ve2*i-(i*VeXe2*(i+1)))
#            
#            wfn2 = wfn2*200
#            yAdjust = wfn2[0] - (Ve2*i-(i*VeXe2*(i+1)))
#            wfn2 = wfn2 - yAdjust
#            
#            x = np.arange (-20,100,.5)
#            x = x*10**-11.35 + re2
#            
#            plt.plot(x,wfn2,'r')
#
#            y1 = Ve2*i-(i*VeXe2*(i+1))
#            
#        if i == v1+1:
#            ax.axhline((Eel + Ve1*i-(i*VeXe1*(i+1))),color = 'y')
#            
#            
#            wfn1 = np.zeros(len(x))
#            for j in np.arange(len(z1)):    
#                wfn1[j] = Nn1 * np.exp(-z1[j]/2) * z1[j] ** (b1/2) * Laguerre1[j] + (Eel + Ve1*i-(i*VeXe1*(i+1)))
#            
#            wfn1 = 200*wfn1
#            yAdjust = wfn1[0] - (Eel + Ve1*i-(i*VeXe1*(i+1)))
#            wfn1 = wfn1 - yAdjust
#            
#            x = np.arange (-20,100,.5)
##            x = x*10**-11.65 + re1-.15e-10    works for 26
#            if v1 < 5:
#                x = x*10**-11.35 + re1
#            elif v1 < 10:
#                x = x*10**-n[v1-5] + re1
#            elif v1 >= 10:
#                x = x*10**-n[v1-5] + re1 - s[v1-10]
#                            
#            plt.plot(x,wfn1,'r')
#            
#            y2 = Eel + Ve1*i-(i*VeXe1*(i+1))

#    xScale = 1e-10
#    xTicks = ticker.FuncFormatter(lambda x,pos: '{0:g}'.format(x/xScale))
#    ax.xaxis.set_major_formatter(xTicks)
    
#    dy = [y2,y1]
#     xLin = [2.666*10**-10,2.666*10**-10]
#    fig,ax = plt.subplots()   

    #draw Morse potentials on figure
    ax.plot(xRange, yExcitedMorse, 'g')
    ax.plot(xRange, yGroundMorse, 'b')
    
###Arrow showing transition from ground state to excited state
#    plt.plot(xLin,dy,'k')
#    plt.annotate('', xy=(2.741*10**-10,y1), xycoords = 'data', xytext=(2.741*10**-10,y2),textcoords='data',arrowprops={'arrowstyle': '<->'})
#####    plt.annotate('',(2.666*10**-10,y2),(2.666*10**-10,y1+15000), arrowprops={'arrowstyle': '->'})


"""Button for showing vibrational states (doesn't work yet)"""
#vibLinesAx = plt.axes([0.69,0.05,0.1,0.075])
#vibLinesButton = Button(vibLinesAx, "Show Vibrational States")
#vibLinesButton.on_clicked(vibStates)


"""n and s used in calculating FCF for scaled(?) anharmonic wave functions"""
#n = []
#for i in np.arange(25):
#    if i <= 5:
#        n.append(11.4 + (.1/5)*i)
#    elif i <= 10:
#        n.append(11.5 + (.02/5)*(i-5))
#    elif i <= 15:
#        n.append(11.52 + (.07/5)*(i-10))
#    elif i <= 20:
#        n.append(11.57 + (.05/5)*(i-15))
#    elif i <= 25:
#        n.append(11.62 + (.15/5)*(i-20))
#s = []
#for i in np.arange(30):
#    if i <= 5:
#        s.append(.05e-10)
#    elif i <= 10:
#        s.append(.06e-10+((.03e-10)/5)*(i-9))
#    elif i <= 15:
#        s.append(.11e-10+((.03e-10)/5)*(i-14))
#    elif i <= 20:
#        s.append(.12e-10+((.07e-10)/5)*(i-14))
          

#Dialogue box for user to select a ground state to excited state transition
master = tk.Tk()
var1 = tk.IntVar()
var2 = tk.IntVar()
tk.Label(master, text = "Choose a Transition:", font=14).grid(row=0, column=0, padx = (10,10), pady = (10,10), sticky=tk.W)
tk.Label(master,text = u"\u03c8 v'':" ).grid(row=1, column=0, padx = (60,10), sticky=tk.W)
tk.Label(master,text = u"\u03c8 v':" ).grid(row=1, column=2, padx = (60,10), sticky=tk.W)
tk.Label(master,text = "to" ).grid(row=1, column=1, padx = (0,20), sticky=tk.W)
tk.Spinbox(master, from_=0, to=20, textvariable = var2).grid(row=2, column=0, padx = (10,10), sticky=tk.W)
tk.Spinbox(master, from_=0, to=26, textvariable = var1).grid(row=2, column=2, padx = (10,10), sticky=tk.W)
tk.Button(master, text = "OK", command = end, width = 15).grid(row=7, column=2, padx = (10,10), pady = (10,10), sticky=tk.W)
tk.mainloop()




#"""For loop used for calculating FCF from one ground state to a range of excited states"""
#for i in np.arange(1): #v''
#    v2 = 0
#    
#    offset = 20
#
#    x2 = np.arange (-20, 100 + offset, 0.5)
##    x2 = np.arange (-20,100,.5)
##    x2 = x2 + re2
#    
#    z2=[]
#    for k in np.arange(len(x2)):
#        z2.append(k2 * np.exp(-alpha2 * (x2[k])))
#    
#    b2 = k2 - 2 * v2 - 1
#    
#    Nn2 = (alpha2 * b2 * float(math.factorial(v2)) / sym.gamma(k2-v2))**(0.5)
#    
#    Laguerre2 = np.zeros(len(z2))
#    for k in np.arange(len(z2)): 
#        Laguerre2[k] = sym.assoc_laguerre(v2,b2,z2[k])
#        
#    wfn2 = np.zeros(len(x2))
#    for k in np.arange(len(z2)):    
#        wfn2[k] = Nn2 * np.exp(-z2[k]/2) * z2[k] ** (b2/2) * Laguerre2[k]
#    
#    FCF = []
#    fig,ax = plt.subplots()
#    for j in np.arange(5,25): #v'
#        v1 = j
#        x1 = np.arange (-20 - offset, 100, 0.5)
##        x1 = np.arange (-20,100,.5)
# 
#        z1=[]
#        for k in np.arange(len(x1)):
#            z1.append(k1 * np.exp(-alpha1 * (x1[k])))
#        
#        b1 = k1 - 2 * v1 - 1
#        
#        Nn1 = (alpha1 * b1 * float(math.factorial(v1)) / sym.gamma(k1-v1))**(0.5)
#        
#        Laguerre1 = np.zeros(len(z1))
#        for k in np.arange(len(z1)): 
#            Laguerre1[k] = sym.assoc_laguerre(v1,b1,z1[k])
#        
#        wfn1 = np.zeros(len(x1))
#        for k in np.arange(len(z1)):    
#            wfn1[k] = Nn1 * np.exp(-z1[k]/2) * z1[k] ** (b1/2) * Laguerre1[k]
#        
#        x1 = x1 + offset    
#        wfn3 = wfn1*wfn2
#        
#        FCF.append(scipy.integrate.simps(wfn3) ** 2)
#                
#        plt.plot(x2, wfn2, 'b')     
#        plt.plot(x1, wfn1, 'g')
        
        
#    fig,ax= plt.subplots()
#    ax.bar(np.arange(14,30),FCF,align='center')
#    #ax.plot(FCF,'o')        
#    ax.set_title(r"$[\int\psi_g \cdot \psi_e dx]^2$")




"""Calculating FCF from one ground state to another excited state and plotting overlap"""
alpha1 = np.sqrt(k1 * u) / (hBar * 6.022e23)
alpha2 = np.sqrt(k2 * u)/(hBar * 6.022e23)
k1 = u*((2 * pi * c * Ve1)**2)
k2 = u*((2 * pi * c * Ve2)**2)
#v2 = 0
#v1 = 14

offset = 20

x2 = np.arange (-20, 100 + offset, 0.5)
#x2 = np.arange (-20,100,.5)

z2=[]
for k in np.arange(len(x2)):
    z2.append(k2 * np.exp(-alpha2 * (x2[k])))

b2 = k2 - 2 * v2 - 1

Nn2 = (alpha2 * b2 * float(math.factorial(v2)) / sym.gamma(k2 - v2))**(0.5)

Laguerre2 = np.zeros(len(z2))
for k in np.arange(len(z2)): 
    Laguerre2[k] = sym.assoc_laguerre(v2, b2, z2[k])

wfn2 = np.zeros(len(x2))
for k in np.arange(len(z2)):    
    wfn2[k] = Nn2 * np.exp(-z2[k] / 2) * z2[k]**(b2/2) * Laguerre2[k]

#x2 = x2*10**-11.35 + re2
   
x1 = np.arange (-20 - offset, 100, 0.5)
#x1 = np.arange (-20,100,.5)

z1=[]
for k in np.arange(len(x1)):
    z1.append(k1 * np.exp(-alpha1 * (x1[k])))

b1 = k1 - 2 * v1 - 1

Nn1 = (alpha1 * b1 * float(math.factorial(v1)) / sym.gamma(k1-v1))**(0.5)

Laguerre1 = np.zeros(len(z1))
for k in np.arange(len(z1)): 
    Laguerre1[k] = sym.assoc_laguerre(v1, b1, z1[k])

wfn1 = np.zeros(len(x1))
for k in np.arange(len(z1)):    
    wfn1[k] = Nn1 * np.exp(-z1[k] / 2) * z1[k]**(b1 / 2) * Laguerre1[k]
    
x1 = x1 + offset    
    
wfn3 = wfn1 * wfn2


FCF = scipy.integrate.simps(wfn3)**2               

#if v1 < 5:
#    x1 = x1*10**-11.35 + re1
#elif v1 < 10:
#    x1 = x1*10**-n[v1-5] + re1
#elif v1 >= 10:
#    x1 = x1*10**-n[v1-5] + re1 - s[v1-10]       


fig, (intAx0,intAx1) = plt.subplots(nrows=2)
intAx0.set_title("Ground State and Excited State Wave Functions")
#intAx0.set_yticks([])
intAx0.plot(x1,wfn1,'g')
intAx0.plot(x1,wfn2,'b')
intAx0.set_ylim(-.2, 0.7)
intAx0.set_xlim(-10, 110)


intAx1.set_title("Product Wave Function")
intAx1.set_title("Internuclear Distance")
#intAx1.set_yticks([])
intAx1.plot(x1,wfn3, 'k')
intAx1.yaxis.set_major_formatter(ScalarFormatter())
#intAx1.set_ylim(0, 30000)
intAx1.set_xlim(-10, 110)
intAx1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))




