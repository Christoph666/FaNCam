import ROOT
from ROOT import TFile, TTree
from ROOT import TCanvas, TGraph
from ROOT import gROOT

import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
import sys
from array import array

dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/"
measurement = "Calibration_ROOT_20190418"

file = ROOT.TFile(dir+"/"+measurement+"/UNFILTERED/compass_"+measurement+".root");

energy, PSD = [], []
n=0
for event in file.Data:
    if event.Board == 0 and event.Channel == 0:
        energy.append(event.Energy)
        PSD.append((event.Energy-event.EnergyShort)/float(event.Energy))
    n+=1
    # if n==1e5:
    #     break
print "Number of events:", n

fScatter = plt.figure("Scatter")
axSc = fScatter.add_subplot(111)
# axSc.hist(energy)
# axSc.hist(PSD)
axSc.hist2d(energy, PSD, range=[[0,2e3],[0,1]], bins=[150, 150], norm=LogNorm())
axSc.set_title(measurement)
axSc.set_xlabel("Energy / ADC")
axSc.set_ylabel("PSD")
plt.show()
