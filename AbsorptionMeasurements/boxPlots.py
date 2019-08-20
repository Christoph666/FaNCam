# plot the absorption of different materials as a box plot together in one graph
# insert material names to x axis

import ROOT
from ROOT import TFile, TTree
from ROOT import TCanvas, TGraph
from ROOT import gROOT

import numpy as np
from array import array
# import matplotlib.colors as colors
# import matplotlib
# matplotlib.use('Agg')
# from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
# from scipy import interpolate
# from scipy.optimize import curve_fit
# from scipy.optimize import curve_fit
import glob, os, sys


dir = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/";
vMeasurements = ["Aluminium_20190418", "Copper_20190418", "Paraffin_20190314", "Stone_20190426", "Lead_20190423"]
widths = ["1cm", "1cm", "5cm", "10cm", "5cm"]
colors=["r", "g", "orange", "blue", "pink", "purple", "turquoise"]


fBox = plt.figure("Absorptions")
# axBox = fBox.add_subplot(111)
axBox = fBox.add_axes((.1,.3,.8,.6))
axRatio=fBox.add_axes((.1,.1,.8,.2))
iBox = 0

for measurement in vMeasurements:

    os.chdir(dir+measurement+"/Images/")
    thisFile = glob.glob("*.root")[0]

    file1 = ROOT.TFile.Open(thisFile)
    neutronRates = np.array(file1.Get("neutronRates_subs"))
    neutronRates_err = np.array(file1.Get("neutronRates_subs_err"))
    gammaRates = np.array(file1.Get("gammaRates_subs"))
    gammaRates_err = np.array(file1.Get("gammaRates_subs_err"))

    avNeutronRate = np.mean(neutronRates)
    avNeutronRate_err = np.std(neutronRates)
    avGammaRate = np.mean(gammaRates)
    avGammaRate_err = np.std(gammaRates)

    axBox.bar(iBox-0.15, avNeutronRate, width=0.2, yerr=avNeutronRate_err, color=colors[iBox], label=measurement[:-9]+" "+widths[iBox])
    axBox.bar(iBox+0.15, avGammaRate, width=0.2, yerr=avGammaRate_err, color=colors[iBox])

    axBox.text(iBox-0.19, avNeutronRate-4, "N")
    axBox.text(iBox+0.1, avGammaRate-4, "G")

    axRatio.plot(iBox, avNeutronRate / avGammaRate, marker="o", ls="--", color=colors[iBox])
    iBox+=1

axBox.set_ylabel("$\Delta_{Rate}/$%")
axRatio.set_ylabel("Ratio N/G")
axRatio.set_ylim(0, 4)
axRatio.hlines(1, 0, iBox-1, color="red", linestyle="--")
axBox.legend()



fBox.savefig("Absorptions.png")
plt.show()
