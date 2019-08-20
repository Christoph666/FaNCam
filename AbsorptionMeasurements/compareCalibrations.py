import ROOT
from ROOT import TFile, TTree
from ROOT import TCanvas, TGraph
from ROOT import gROOT

import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colorbar as cbar
from scipy.optimize import curve_fit
from array import array
import glob, os, sys
sys.path.insert(0, "/home/christoph/Documents/MasterThesis/Analyse/Utlities/")
import python_utilities as util

save_output = 0
mapping = ("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4")

dir = "/home/christoph/Documents/MasterThesis/Data/Absorption/"
measurements = ["Al/Calibration_20190529", "Al/Calibration_20190603", "Al/Calibration_20190604_0", "Al/Calibration_20190604_1",  "Al/Calibration_20190611", "Cu/Calibration_20190701",
 "Fe/Calibration_20190624", "H2O/Calibration_20190617"]

allNeutronRates=[]
allNeutronRatesErr=[]
allGammaRates=[]
allGammaRatesErr=[]

for measurement in measurements:
    thisfile = dir + measurement + "/Results/" + measurement.split("/")[1] + "_Results.root";
    print "reading:", thisfile
    file1 = ROOT.TFile.Open(thisfile)

    neutronRates = (np.array(file1.Get("neutronRates")))
    neutronRates_err = np.array( file1.Get("neutronRates_err"))
    gammaRates = (np.array(file1.Get("gammaRates")))
    gammaRates_err = np.array( file1.Get("gammaRates_err"))

    # print "abs"
    # neutronRates_abs = np.array(file1.Get("neutronRates_abs"))
    # neutronRates_abs_err = np.array(file1.Get("neutronRates_abs_err"))
    # print neutronRates_abs_err[0] / neutronRates_abs[0] * 100.
    # print "calib"
    # neutronRates_calib = np.array(file1.Get("neutronRates_calib"))
    # neutronRates_calib_err = np.array(file1.Get("neutronRates_calib_err"))
    # print neutronRates_calib_err[0] / neutronRates_calib[0] * 100.
    # print "total: ", np.sqrt( (neutronRates_abs_err[0] / neutronRates_abs[0])**2 + (neutronRates_calib_err[0] / neutronRates_calib[0])**2 ) * 100.

    allNeutronRates.append(neutronRates)
    allNeutronRatesErr.append(neutronRates_err)
    allGammaRates.append(gammaRates)
    allGammaRatesErr.append(gammaRates_err)

# allNeutronRates /= np.mean(allNeutronRates[0])
# allNeutronRatesErr /= np.mean(allNeutronRates[0])
# allGammaRates /= np.mean(allGammaRates[0])
# allGammaRatesErr /= np.mean(allGammaRates[0])

f1 = plt.figure("f1")
axN = f1.add_subplot(211)
axG = f1.add_subplot(212)
axN.set_title("Neutrons")
axG.set_title("Gammas")
axG.set_xlabel("Pixel")
axN.set_ylabel("Relative rate")
axG.set_ylabel("Relative rate")
axN.set_xticks([])
axG.set_xticks(range(16))
for i in range(len(measurements)):
    axN.errorbar(range(16), allNeutronRates[i], yerr=allNeutronRatesErr[i], label=measurements[i], marker="o", markersize=5, linestyle="--", linewidth=0.5)
    axG.errorbar(range(16), allGammaRates[i], yerr=allGammaRatesErr[i], label=measurements[i], marker="o", markersize=5, linestyle="--", linewidth=0.5)
# axN.legend(loc='upper center', bbox_to_anchor=(0.5, -1.4), fancybox=True, shadow=True, ncol=4)
# f1.tight_layout()
axG.legend()

# print np.mean(allNeutronRates, axis=1)

# f2 = plt.figure("f2")
# axN_mean = f2.add_subplot(111)
# # axN_mean = f2.add_subplot(211)
# # axG_mean = f2.add_subplot(212)
# x = np.arange(1,len(measurements)+1)
# # axN_mean.errorbar(x, allNeutronRates[0], yerr=allNeutronRatesErr[0], marker="o", ls="", label="Neutrons")
# axN_mean.errorbar(x, np.mean(allNeutronRates, axis=1), yerr=np.std(allNeutronRates, axis=1), marker="o", ls="", label="Neutrons")
# axN_mean.errorbar(x, np.mean(allGammaRates, axis=1), yerr=np.std(allGammaRates, axis=1), marker="o", ls="", label="Gammas")
#
# if not fitWithOffset:
#     poptN, pcovN = curve_fit(exponential, x, np.mean(allNeutronRates, axis=1), sigma=np.std(allNeutronRates, axis=1), p0=[0.9,8.])
#     poptG, pcovG = curve_fit(exponential, x, np.mean(allGammaRates, axis=1), sigma=np.std(allGammaRates, axis=1), p0=[0.9,8.])
#     axN_mean.plot(x, exponential(x, *poptN), 'r-', label="Neutrons:\nA=%5.1f\n$x_{1/2}$=%5.1f cm" % tuple([poptN[0], poptN[1]]))
#     axN_mean.plot(x, exponential(x, *poptG), 'r-', label="Gammas:\nA=%5.1f\n$x_{1/2}$=%5.1f cm" % tuple([poptG[0], poptG[1]]))
#     print chisq(x, np.mean(allNeutronRates, axis=1), np.std(allNeutronRates, axis=1), exponential, poptN)
#     print chisq(x, np.mean(allGammaRates, axis=1), np.std(allGammaRates, axis=1), exponential, poptG)
#     axN_mean.text(6, 0.7, r"$f(x)=A\cdot \exp(-x/x_0)$")
#
# if fitWithOffset:
#     poptN, pcovN = curve_fit(exponential_offset, x, np.mean(allNeutronRates, axis=1), sigma=np.std(allNeutronRates, axis=1), p0=[0.9,8.,0.])
#     poptG, pcovG = curve_fit(exponential_offset, x, np.mean(allGammaRates, axis=1), sigma=np.std(allGammaRates, axis=1), p0=[0.9,8.,0.])
#     axN_mean.plot(x, exponential_offset(x, *poptN), 'r-', label="Neutrons:\nA=%5.1f\n$x_{1/2}$=%5.1f cm\nA_0=%5.1f" % tuple([poptN[0], poptN[1]*np.log(2.), poptN[2]]))
#     axN_mean.plot(x, exponential_offset(x, *poptG), 'r-', label="Gammas:\nA=%5.1f\n$x_{1/2}$=%5.1f cm\nA_0=%5.1f" % tuple([poptG[0], poptG[1]*np.log(2.), poptG[2]]))
#     print "chi square n", chisq(x, np.mean(allNeutronRates, axis=1), np.std(allNeutronRates, axis=1), exponential_offset, poptN)
#     print "chi square g", chisq(x, np.mean(allGammaRates, axis=1), np.std(allGammaRates, axis=1), exponential_offset, poptG)
#     axN_mean.text(6, 0.7, r"$f(x)=A\cdot \exp(-x/x_0) + A_0$")
#
#
#
# axN_mean.set_title("")
# axN_mean.set_xlabel("Thickness / cm")
# # axN_mean.set_ylabel("Decrease in rate / %")
# axN_mean.set_ylabel("Relative rate")
# axN_mean.set_xticks(range(1,len(measurements)+1))
# axN_mean.legend()
# # axG_mean.set_title("Neutrons")


plt.show()
