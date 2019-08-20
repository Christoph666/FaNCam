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

# sys.path.append("/home/christoph/Documents/Python/usefull_stuff/")
# from tempfile import TemporaryFile
# from read import *
# from write import *


savePlots = 0
do_other_stuff = 0

def exponential(x, N, tau):
    return N * np.exp(-x/tau)


# 1px 54.9 overvoltage
# filename = "Messung1"
# file1 = ROOT.TFile.Open("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung1/Plots/FOM/Messung1_25_270/Messung1_25_270_FOM_results.root")
# file2 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung1/Plots/Messung1_25_270_PSD.root")
# file3 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung1/"+filename+".root")

# 1px 53.9 overvoltage
# filename = "Messung1_50cm_539V"
# file1 = ROOT.TFile.Open("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung1_50cm_539V/Plots/FOM/Messung1_50cm_539V_25_270/Messung1_50cm_539V_25_270_FOM_results.root")
# file2 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung1_50cm_539V/Plots/Messung1_50cm_539V_25_270_PSD.root")
# file3 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung1_50cm_539V/"+filename+".root")

# 16px 54.9 overvoltage
# filename = "Messung3_0_0"
# file1 = ROOT.TFile.Open("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung3/Plots/FOM/Messung3_0_0/Messung3_0_0_FOM_results.root")
# file2 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung3/Plots/Messung3_0_0_25_270_PSD_full.root");
# file3 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung3/"+filename+".root");

# 16px 53.9 overvoltage
filename = "Messung3_539V_0_0"
file1 = ROOT.TFile.Open("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung3_539V/Plots/FOM/Messung3_539V_0_0/Messung3_539V_0_0_FOM_results.root")
file2 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung3_539V/Plots/Messung3_539V_0_0_25_270_PSD.root");
file3 = ROOT.TFile("/home/christoph/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung3_539V/"+filename+".root");

vSlices_n = np.array(file1.Get("vSlices_n"))
vSlices_g = np.array(file1.Get("vSlices_g"))
vSlices_FOM = np.array(file1.Get("vSlices_FOM"))
mu_n = np.array(file1.Get("mu_n"))
sig_n = np.array(file1.Get("sig_n"))
mu_g = np.array(file1.Get("mu_g"))
sig_g = np.array(file1.Get("sig_g"))

vPeak = np.array(file2.Get("vPeak"))
vPSD = np.array(file2.Get("vPSD"))
vValid = np.array(file2.Get("vValid"))

f_mu_n = interpolate.interp1d(vSlices_n, mu_n, kind='cubic')
f_sig_n = interpolate.interp1d(vSlices_n, sig_n, kind='cubic')
f_mu_g = interpolate.interp1d(vSlices_g, mu_g, kind='cubic')
f_sig_g = interpolate.interp1d(vSlices_g, sig_g, kind='cubic')

f1 = plt.figure("f1")
ax1 = f1.add_subplot(111)
ax1.set_title(filename)
ax1.set_xlabel("PH / mV")
ax1.set_ylabel("PSD")
ax1.plot(vSlices_n, mu_n + 3*sig_n, color="green")
ax1.plot(vSlices_n, mu_n - 3*sig_n, color="green")
ax1.plot(vSlices_g, mu_g + 3*sig_g, color="red")
ax1.plot(vSlices_g, mu_g - 3*sig_g, color="red")
plt.hist2d(vPeak, vPSD, range=[[0,1000],[0,1]], bins=[500, 500], norm=LogNorm())
plt.colorbar()
if savePlots:
    f1.savefig("Plots/"+filename+"_psd2d.pdf")

iAccEvent=0
iEvent=0
nCount=0
gCount=0
Vbase=0
Nbase=15

if filename == "Messung1_50cm_539V" or filename=="Messung1":
    nSamples=999
    iStart=int(0.13*nSamples)
    iEnd=int(0.4*nSamples)
if filename == "Messung3_0_0" or filename=="Messung3_539V_0_0":
    nSamples=743
    iStart=int(0.15*nSamples)
    iEnd=int(0.55*nSamples)

x = np.linspace(0, nSamples, nSamples) * 2 - 300
avNeutronPulse = np.zeros(nSamples)
avGammaPulse = np.zeros(nSamples)

# average pulses
f2 = plt.figure("f2")
ax2 = f2.add_subplot(111)

# plot valid stamps
# f_valid = plt.figure("f_valid")
# ax_valid = f_valid.add_subplot(111)
# nslice = -1
# ax_valid.plot(range(len(vValid))[0:nslice], vValid[0:nslice], lw=0.2)
# f_valid.savefig("Plots/valid.png")

# delta t
f3 = plt.figure("f3")
ax3 = f3.add_subplot(211)
ax3_1 = f3.add_subplot(212)
timeTags = []
deltaT = []

if do_other_stuff:

    # time tags
    f4 = plt.figure("f4")
    ax4 = f4.add_subplot(111)

    # psd plot with caen data
    f5 = plt.figure("f5")
    ax5 = f5.add_subplot(111)

    # channel display
    f6 = plt.figure("f6")
    ax6 = f6.add_subplot(111)

    endTimes = []
    vQS = []
    vQL = []
    vPSD_HW = []
    vPH_HW = []
    vCHN = []

ntot = 0
Nmax=1e3

for event in file3.Events:
    if vValid[iEvent] == 0:
        print "iEvent: ", iEvent
        print "iAccEvent: ", iAccEvent
    if vValid[iEvent] == 1:
        if iEvent<Nmax:
            Vbase=0
            samples = np.array(event.samples, dtype="float64") * 2e3/2**14

            # if do_other_stuff:
            timeTag = np.uint32(event.triggerTimeTag)*1./1e6
            # vPSD_HW.append(float((event.energy)-(event.energyShort)) / (event.energy))
            # vPH_HW.append(vPeak[iAccEvent])
            if iEvent == 0:
                t0=timeTag
            timeTag-=t0
            timeTags.append(timeTag)
            # endTimes.append(timeTag+768.*2./1e6)
            if iEvent > 0:
                deltaT.append(timeTags[iAccEvent] - timeTags[iAccEvent-1])

            for j in xrange(Nbase):
                Vbase += samples[j]/Nbase
            samples -= Vbase
            samples = -samples
            peak = vPeak[iAccEvent]
            psd = vPSD[iAccEvent]

            if peak<vSlices_FOM[-1] and peak>vSlices_FOM[0]:
                if psd<f_mu_n(peak) + 3*f_sig_n(peak) and psd>f_mu_n(peak) - 3*f_sig_n(peak) and psd>f_mu_g(peak) + 3*f_sig_g(peak):
                    # print "neutron"
                    if max(samples) < 940.:
                        avNeutronPulse += samples
                        # print samples
                        # print ""
                    nCount+=1
            if peak<vSlices_FOM[-1] and peak>vSlices_FOM[0]:
                if psd<f_mu_g(peak)+3*f_sig_g(peak) and psd>f_mu_g(peak)-3*f_sig_g(peak) and psd<f_mu_n(peak)-3*f_sig_n(peak):
                    # print "gamma"
                    if max(samples) < 940.:
                        avGammaPulse += samples
                    gCount += 1
        iAccEvent+=1
    iEvent+=1
    if iEvent==Nmax:
        print "Breaking."
        break
    ntot+=1

print "ntot: ", ntot
print "nacc: ", iAccEvent
print "neutrons: ", nCount
print "gammas: ", gCount

if do_other_stuff:
    print "mean delta t: ", np.mean(deltaT), " ms"
    print "stddev delta t: ", np.std(deltaT), " ms"
    print "mean rate: ", 1./np.mean(deltaT), "kHz"

avNeutronPulse /= float(nCount)
avGammaPulse /= float(gCount)

avNeutronPulse /= max(avNeutronPulse)
avGammaPulse /= max(avGammaPulse)

# average pulses
ax2.plot(x[iStart:iEnd], avNeutronPulse[iStart:iEnd], color="green", label="Neutron Pulse")
ax2.plot(x[iStart:iEnd], avGammaPulse[iStart:iEnd], color="red", label="Gamma Pulse")

# plt.title("Average Pulses")
ax2.set_xlabel("t / ns")
ax2.set_ylabel("normed pulse height")
ax2.legend(loc="best")
if savePlots:
    f2.savefig("Plots/"+filename+"_average_pulse.png")
    f2.savefig("Plots/"+filename+"_average_pulse.pdf")
# write("av_pulse.txt", [x, avNeutronPulse, avGammaPulse], 3)

xArr, yArrN, yArrG = array('d'), array('d'), array('d')
for i in range(len(avNeutronPulse)):
    xArr.append(i*2)
    yArrN.append(avNeutronPulse[i])
    yArrG.append(avGammaPulse[i])

grN = TGraph(len(avNeutronPulse), xArr, yArrN)
grN.SetLineColor(2)
grN.SetLineWidth(4)
grN.SetMarkerColor(4)
grN.SetMarkerStyle(21)
grN.SetTitle('Average Neutron Pulse')
grN.GetXaxis().SetTitle('t / ns')
grN.GetYaxis().SetTitle('Normed Pulse Height')
grN.Draw("ACP")

grG = TGraph(len(avNeutronPulse), xArr, yArrG)
grG.SetLineColor(2)
grG.SetLineWidth(4)
grG.SetMarkerColor(4)
grG.SetMarkerStyle(21)
grG.SetTitle('Average Gamma Pulse')
grG.GetXaxis().SetTitle('t / ns')
grG.GetYaxis().SetTitle('Normed Pulse Height')
grG.Draw("ACP")

f_out = TFile("Plots/"+filename+".root", "recreate")
grN.Write()
grG.Write()




# hist of delta t with exp fit
range=[0,5e3]
nBins=100
counts, bin_edges = np.histogram(deltaT, range=range, bins=nBins)
bins = (bin_edges[:-1] + bin_edges[1:]) / 2
popt, pcov = curve_fit(exponential, bins, counts, p0=[50.,500.])
ax3.plot(bins, exponential(bins, *popt), 'r-', label="exp fit:\nN=%5.1f\ntau=%5.1f mu s" % tuple(popt))

# ax3.yscale('log')
ax3.set_title(r"$\Delta_{t}$ between pulses")
ax3.hist(deltaT, range=range, bins=nBins)
ax3_1.hist(deltaT, range=[0, 1e2], bins=30)
ax3_1.set_xlabel("mu s")
ax3.set_ylabel("")
ax3.legend(loc="best")
if savePlots:
    f3.savefig("Plots/"+filename+"_deltaT.png")



if do_other_stuff:

    # time tags
    ax4.vlines(timeTags, 0,1, lw=0.02, color="black", alpha=0.7)
    # ax4.vlines(endTimes, 0,1, lw=0.02, color="red", alpha=0.7)
    # ax4.set_xlim([0, 1e2])
    ax4.set_xlabel("ms")
    if savePlots:
        f4.savefig("Plots/timeTags.pdf")

    # psd plot caen
    # print vPH_HW
    ax5.hist2d(vPH_HW, vPSD_HW, range=[[0,1000],[0,1]], bins=[300, 300], norm=LogNorm())
    # ax5.plot(range(len(vPH_HW)), vPH_HW)
    # ax5.plot(range(len(vPSD_HW)), vPSD_HW)
    if savePlots:
        f5.savefig("Plots/"+filename+"_psd_caen.pdf")

    # vCHN.append(7)
    # print len(vCHN)
    # ax6.hist(vCHN, range=[0,8], bins=8)
    # # ax6.set_xlim(0,8)
    # f6.savefig("channels.pdf")


file1.Close()
file2.Close()
file3.Close()

plt.show()
