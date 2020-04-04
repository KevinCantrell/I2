import sys
import pandas as pd
import numpy as np
import Tkinter as tk
import tkFileDialog
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.widgets import Button
from scipy.signal import savgol_filter
from numpy import NaN, Inf, arange, isscalar, asarray, array

#Instance Variables
PeakTable = pd.DataFrame(columns=['Wavelength (nm)', 'Absorbance', "v''", '',"v'","v'' = 0 (nm)","v'","v'' = 1 (nm)","v'","v'' = 2 (nm)",])
WavenumberTable = pd.DataFrame(columns = ["v'+0.5","v'' = 0 (cm-1)","v'+0.5","v'' = 1 (cm-1)","v'+0.5","v'' = 2 (cm-1)",])
IntensityFrame = pd.DataFrame(columns = ["Wavelength (nm)","Intensity"])
WavenumberTable["v'+0.5"] = np.arange(50)
testData = []
coords = []
line = []
intensities = []
waveln = []
absorb=[]
index = -2
v0Waveln = []
v1Waveln = []
v2Waveln = []
v0Intensity = []
v1Intensity = []
v2Intensity = []

"""
peakdet
Finds peaks and valleys of a data set
"""
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

"""
on_pick_transition
When user clicks a point on the transitions graph, the point
will change color and be added to the list coords
"""
def on_pick_transition(event):
    global coords, testData,coll,colors,index

    #check to see if user clicked on more than one point at a time, if so choose the first point
    if (len(event.ind) > 1):
        index = np.array([event.ind[0]])
    else:
        index = event.ind

    double = False
    doubleIndex = -1
    #search and see if clicked point has been clicked before, if so remove from coords
    if len(coords) > 0:
        for i in np.arange(len(coords)):
            if (testData[index][0][0] == coords[i][0][0]):
                double = True
                doubleIndex = i
                break
        #if its a double, delete it from coords and change it back to blue
        if double:
            coll._facecolors[index,:] = (0, 0, 1, 1)
            coll._edgecolors[index,:] = (0, 0, 1, 1)
            del(coords[doubleIndex])
        #if it wasnt a double, add to coords and change to red
        else:
            coll._facecolors[index,:] = (1, 0, 0, 1)
            coll._edgecolors[index,:] = (1, 0, 0, 1)
            coords.append(testData[index])
    #if its the first point, it cant be a double so add it to coords and change to red
    else:
        coll._facecolors[index,:] = (1, 0, 0, 1)
        coll._edgecolors[index,:] = (1, 0, 0, 1)
        coords.append(testData[index])
            
    fig1.canvas.draw()

"""
on_pick_intensity
When user clicks a point on the intensity graph, the point
will change color and be added to the list line
"""
def on_pick_intensity(event):
    global  lineData,peaks,valleys,line,index
   
    #check to see if user clicked on more than one point at a time, if so choose the first point
    if (len(event.ind) > 1):
        index = np.array([event.ind[0]])
    else:
        index = event.ind
          
    double = False
    doubleIndex = -1
    #search and see if clicked point has been clicked before, if so remove from line
    if len(line) > 0:
        for i in np.arange(len(line)):
            if (lineData[index][0][0] == line[i][0][0]):
                double = True
                doubleIndex = i
                break        
        #if its a double, delete it from line and change it back to blue
        if double:
            valleys._facecolors[index,:] = (0, 0, 1, 1)
            valleys._edgecolors[index,:] = (0, 0, 1, 1)
            del(line[doubleIndex])
        #if it wasnt a double, add to line and change to red
        else:
            valleys._facecolors[index,:] = (1, 0, 0, 1)
            valleys._edgecolors[index,:] = (1, 0, 0, 1)
            line.append(lineData[index])
    #if its the first point, it cant be a double so add it to coords and change to red
    else:
        valleys._facecolors[index,:] = (1, 0, 0, 1)
        valleys._edgecolors[index,:] = (1, 0, 0, 1)
        line.append(lineData[index])

    fig2.canvas.draw()

"""
axButton Class
contains functions 
v0, v1, v2, finishPeak, plotFit, and finishValley
to be used for buttons on Transitions and Intensity Figures
"""
class axButton:
    def __init__(self,PeakTable):
        self.PeakTable = PeakTable
    #arrays used to hold and manipulate peak data         
    a = []
    d = []
    c = []
    e = []
    f = []
    g = []
    
    """Transition Functions"""    
    
    """
    v0
    When user clicks v'' = 0 button, the points previously clicked will be
    added to a corresponding dataframe and redrawn on the graph in green.
    A new graph will display representing v'' vs v'+1/2
    """
    def v0(self,event):
        global coords
        #store x and y coordinates of clicked points 
        a = array(coords)[:,0]
        c = array (a)[:,0]
        d = array(a)[:,1]
    
        #draw clicked points in green for v''=0 transitions
        ax1.scatter(c,d,color=["green"])
        plt.draw()
        
        #data set to plot v'' vs v' + 1/2
        e = arange(14.5,14.5+len(coords))
        f = pd.DataFrame((1*10**7/c))
        
        #new figure for v'' vs v'+1/2
        figV,axV = plt.subplots()
        axV.scatter(e,f)
        plt.title("v vs. v' + 0.5 for v'' = 0")
        
        #add to data frame
        WavenumberTable["v'' = 0 (cm-1)"] = f
        #index on data frame
        for index in arange(len(coords)):
            PeakTable.loc[(coords[index][0,0]==PeakTable['Wavelength (nm)']),"v''"] = "0"
            
        #clear coords for next set of peaks
        coords = []
        
    """
    v1
    When user clicks v'' = 1 button, the points previously clicked will be
    added to a corresponding dataframe and redrawn on the graph in purple.
    A new graph will display representing v'' vs v'+1/2
    """
    def v1(self,event):
        global coords
        #store x and y coordinates of clicked points 
        a = array(coords)[:,0]
        c = array (a)[:,0]
        d = array(a)[:,1]
        
        #draw clicked points in purple for v''=1 transitions
        ax1.scatter(c,d,color=["purple"])
        plt.draw()
        
        #data set to plot v'' vs v' + 1/2
        e = arange(14.5,14.5+len(coords))
        f = pd.DataFrame((1*10**7/c))
        
        #new figure for v'' vs v'+1/2
        figV,axV = plt.subplots()
        axV.scatter(e,f)
        plt.title("v vs. v' + 0.5 for v'' = 1")
        
        #add to data frame
        WavenumberTable["v'' = 1 (cm-1)"] = f
        #index on data frame
        for index in arange(len(coords)):
            PeakTable.loc[(coords[index][0,0]==PeakTable['Wavelength (nm)']),"v''"] = "1"
        
        #clear coords for next set of peaks
        coords = []

    """
    v2
    When user clicks v'' = 2 button, the points previously clicked will be
    added to a corresponding dataframe and redrawn on the graph in yellow.
    A new graph will display representing v'' vs v'+1/2
    """
    def v2(self,event):
        global coords
        #store x and y coordinates of clicked points 
        a = array(coords)[:,0]
        c = array (a)[:,0]
        d = array(a)[:,1]
        
        #draw clicked points in yellow for v''=2 transitions
        ax1.scatter(c,d,color=["yellow"])
        plt.draw()
        
        #data set to plot v'' vs v' + 1/2
        e = arange(14.5,14.5+len(coords))
        f = pd.DataFrame((1*10**7/c))
        
        #new figure for v'' vs v'+1/2
        figV,axV = plt.subplots()
        axV.scatter(e,f)
        plt.title("v vs. v' + 0.5 for v'' = 2")
        
        #add to data frame
        WavenumberTable["v'' = 2 (cm-1)"] = f
        #index on data frame
        for index in arange(len(coords)):
            PeakTable.loc[(coords[index][0,0]==PeakTable['Wavelength (nm)']),"v''"] = "2"
        
        #clear coords for next set of peaks
        coords = []
        
    """
    finishPeak
    When the 'Finish' button is clicked on the Transitions plot, three columns 
    for v'' = 0, 1 and 2 are created based on user input. The data is then exported
    to an excel file    
    """    
    def finishPeak(self,event):
        #find peaks for v''=0 and add them to a new column in the original data frame
        p = PeakTable['Wavelength (nm)'][(PeakTable["v''"]=="0")]
        q = p.to_frame(name="")
        r = pd.DataFrame(q.values)
        PeakTable["v'' = 0 (nm)"] = r
        #find peaks for v''=1 and add them to a new column in the original data frame        
        p = PeakTable['Wavelength (nm)'][(PeakTable["v''"]=="1")]
        q = p.to_frame(name="")
        r = pd.DataFrame(q.values)
        PeakTable["v'' = 1 (nm)"] = r
        #find peaks for v''=2 and add them to a new column in the original data frame       
        p = PeakTable['Wavelength (nm)'][(PeakTable["v''"]=="2")]
        q = p.to_frame(name="")
        r = pd.DataFrame(q.values)
        PeakTable["v'' = 2 (nm)"] = r
        #add blank spaces between columns for user to fill out in Excel
        WavenumberTable["v'+0.5"] = ""
                  
        #export data to an Excel file
        writer = pd.ExcelWriter(file_path+'_Processed_Transition_Data.xlsx')
        PeakTable.to_excel(writer,'Peak Data',index=False)
        WavenumberTable.to_excel(writer,"v'' vs v' + 0.5",index=False)
        writer.save()
        
        
    """Intensity Functions"""

    """
    plotFit
    User clicked points are made into a line of best fit through linear regression.
    Peaks above this line are stored and colored in yellow, the line segment is 
    drawn in green.
    """
    def plotFit(self,event):
        global line,testData,intensities,waveln,absorb
        #arrays used to manipulate and store data
        xFit = []
        yFit = []
        Rline= []
        peakPoint = []
        
        #iterate through clicked points, store x-y coordinates
        for i in arange(len(line)):
            xFit.append(line[i][0][0])
            yFit.append(line[i][0][1])
        
        #run a linear regression on the clicked points
        slope,intercept,r,p,stdErr = stats.linregress(xFit,yFit)
        
        #create a line from the regression data
        for j in arange(len(xFit)):
            Rline.append(slope*xFit[j]+intercept)
            
        #create points to plot from the line
        p1 = array([line[0][0][0],line[0][0][1]])
        p2= array([line[-1][0][0],line[-1][0][1]]) 
         
        for k in arange(len(testData)):
            if testData[k][0] > line[-1][0][0] and testData[k][0] < line[0][0][0]:
                peakPoint.append(testData[k])
        
        #iterate through peaks, find distance from line of fit to get relative intensity
        for l in arange(len(peakPoint)):
            p3= array([peakPoint[l][0],peakPoint[l][1]])
            intensity = np.linalg.norm(np.cross(p2-p1,p1-p3))/np.linalg.norm(p2[0]-p1[0])
            
            #store clicked point
            waveln.append(peakPoint[l][0])
            absorb.append(peakPoint[l][1])
            
            #store intensities
            intensities.append(intensity)

        #graph line of fit and peaks above it
        ax2.plot(xFit,Rline,'g')
        ax2.plot(waveln,absorb,'oc')
        
        #clear line, peakPoint and Rline for next set of points
        line=[]
        peakPoint=[]
        Rline=[]
    
    """
    finishValley
    When the user clicks the 'Finish' button on the intensity graph, the wavelength 
    and the intensity of each peak are added to a dataframe and a new figure is created
    depicting the intensity of v'' = 0, 1 and 2 to v' transitions
    """
    def finishValley(self,event):
        global intensities,waveln,v0Waveln,v1Waveln,v2Waveln,v0Intensity,v1Intensity,v2Intensity,IntensityFrame 
        #store total wavelength and intensities 
        IntensityFrame["Wavelength (nm)"] = waveln
        IntensityFrame["Intensity"] = intensities             
        
        #iterate through points and known states to separate v'' states by intensity
        for i in np.arange(len(waveln)):
            for j in np.arange(len(PeakTable["Wavelength (nm)"])):
                #find v''=0 peaks
                if (waveln[i] == PeakTable["Wavelength (nm)"][j]) and (PeakTable["v''"][j] == "0"):
                    v0Waveln.append(waveln[i])
                    v0Intensity.append(intensities[i])
                #find v''=1 peaks
                elif (waveln[i] == PeakTable["Wavelength (nm)"][j]) and (PeakTable["v''"][j] == "1"):
                    v1Waveln.append(waveln[i])
                    v1Intensity.append(intensities[i])
                #find v''=2 peaks
                elif (waveln[i] == PeakTable["Wavelength (nm)"][j]) and (PeakTable["v''"][j] == "2"):
                    v2Waveln.append(waveln[i])
                    v2Intensity.append(intensities[i])
        
        #create bar graph showing all three ground state to excited state Franck-Condon Factors
        fig, (intAx0,intAx1,intAx2) = plt.subplots(nrows=3)
        #v''=0
        intAx0.invert_xaxis()
        intAx0.bar(v0Waveln,v0Intensity,align='center',color = 'g')
        intAx0.set_title("v'' = 0 to v' Transition Intensities")
        #v''=1
        intAx1.bar(v1Waveln,v1Intensity,align='center',color = 'm')
        intAx1.invert_xaxis()
        intAx1.set_title("v'' = 1 to v' Transition Intensities")
        intAx1.set_ylabel("Intensity")
        #v''=2
        intAx2.bar(v2Waveln,v2Intensity,align='center',color = 'y')
        intAx2.invert_xaxis()
        intAx2.set_title("v'' = 2 to v' Transition Intensities")
        intAx2.set_xlabel("Wavelength (nm)")
        plt.tight_layout(pad=.4,w_pad=.5,h_pad=1.0)
        
        #add Franck-Condon Factors to data frame
        v0Transitions= pd.DataFrame({"v'' = 0 Peak (nm)": v0Waveln,"v'' = 0 Intensity": v0Intensity})
        v1Transitions= pd.DataFrame({"v'' = 1 Peak (nm)": v1Waveln,"v'' = 1 Intensity": v1Intensity})
        v2Transitions= pd.DataFrame({"v'' = 2 Peak (nm)": v2Waveln,"v'' = 2 Intensity": v2Intensity})
        IntensityTable = pd.concat([IntensityFrame,v0Transitions,v1Transitions,v2Transitions],axis=1)

        #export data to Excel
        writer = pd.ExcelWriter(file_path+'_Processed_Intensity_Data.xlsx')
        IntensityTable.to_excel(writer,'Intensity Data',index=False)
        writer.save()
        
        
"""Import data from csv file"""
#setting up interface for dialogue box 
root = tk.Tk()
root.withdraw()
file_path = tkFileDialog.askopenfilename()
df=pd.read_csv(file_path)
#create a dataframe of selected data
df=pd.read_csv(file_path,header=1,names=['nm','absorbance'],index_col=False)


"""Range of consideration in nm (CAN EDIT)"""
#lower wavenumber
lw=505
#higher wavenumber
hw=620

"""Smooth data"""
#addressing on data frame, look at absorbance and wavenumbers between selected range for x and y 
SegmentY=np.array(df['absorbance'][(df['nm']>=lw) & (df['nm']<=hw)])
SegmentX=np.array(df['nm'][(df['nm']>=lw) & (df['nm']<=hw)])

#saviski-golay smoothing of data
ySmoothed = savgol_filter(SegmentY, 15, 2)

#finding peaks and highlighting them on a graph
series = ySmoothed
maxtab, mintab = peakdet(series,.0003,x=SegmentX) #adjust sensitivity

#store peaks and valleys                    
v = array(mintab)[:,0]      #v - valley x positon
w = array(mintab)[:,1]      #w - valley y position
x = array(maxtab)[:,0]      #x - peak x position
y = array(maxtab)[:,1]      #y - peak y position
testData = maxtab
lineData = mintab

#Place peak data into PeakTable dataframe
PeakTable['Wavelength (nm)'] = x
PeakTable['Absorbance'] = y

#Plot for finding transitions
fig1, ax1 = plt.subplots()
ax1.invert_xaxis()
plt.subplots_adjust(bottom=0.2)
plt.title("Select the v'' to v' Transition Peaks")
plt.xlabel("Wavelength (nm)",position=(.2,.5))
plt.ylabel("Intensity")

#Vertical lines depicting the v'' to v'=14 transitions to help with labeling and orienting peaks
line0 = plt.axvline(x=577.4,color = 'g',label = "v'' = 0 to v'= 14")
line1 = plt.axvline(x=584.7,color = 'purple',label ="v'' = 1 to v'= 14")
line2 = plt.axvline(x=592.0,color = 'y',label ="v'' = 2 to v'= 14")

ax1.legend(handles=[line0,line1,line2])
ax1.plot(SegmentX,ySmoothed)

#Buttons v''=0, v''=1, v''=2 and Finished to help users choose appropriate peaks
button2 = axButton(PeakTable)
axV0 = plt.axes([0.47,0.05,0.1,0.075])
axFinished = plt.axes([0.8,0.05,0.1,0.075])
axV1 = plt.axes([0.58,0.05,0.1,0.075])
axV2 = plt.axes([0.69,0.05,0.1,0.075])
bV0 = Button(axV0, "v''=0")
bV1 = Button(axV1, "v''=1")
bV2 = Button(axV2, "v''=2")
bfinished = Button(axFinished, "Finished")
bV0.on_clicked(button2.v0)
bV1.on_clicked(button2.v1)
bV2.on_clicked(button2.v2)
bfinished.on_clicked(button2.finishPeak)

#scatter plot of peaks that are clickable
coll = ax1.scatter(x, y, color = ['blue']*len(testData), picker = True, s=[50]*len(testData))
fig1.canvas.mpl_connect('pick_event', on_pick_transition)
plt.show()

#Plot for finding intensities
fig2, ax2 = plt.subplots()
ax2.invert_xaxis()
plt.subplots_adjust(bottom=0.2)
plt.title("Select the Valleys to Get Fit Lines Along the Bottom")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")

#Buttons Show Fit and Finished to find bottom line to find intensities from
button1 = axButton(PeakTable)
fitAx = plt.axes([0.69,0.05,0.1,0.075])
doneAx = plt.axes([0.8,0.05,0.1,0.075])
fitButton = Button(fitAx, "Show Fit")
doneButton = Button(doneAx, "Finished")
fitButton.on_clicked(button1.plotFit)
doneButton.on_clicked(button1.finishValley)

ax2.plot(SegmentX,ySmoothed,'-k')
peaks = ax2.scatter(x,y, color=['green'])

#clickable valley points
valleys = ax2.scatter(v,w,color=['blue']*len(lineData), picker = True, s=[50]*len(lineData))
fig2.canvas.mpl_connect('pick_event', on_pick_intensity)
plt.show()
