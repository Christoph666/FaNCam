import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import time
from tqdm import tqdm

# dir = "~/Documents/MasterThesis/Data/CoMPASS/run_Cs137/"
# dir = "~/Documents/MasterThesis/Data/CoMPASS/run_AmBe/"
# dir = "~/Documents/MasterThesis/Data/CoMPASS/testrun_long_reclen/"
# dir = "~/Documents/MasterThesis/Data/CoMPASS/testrun_10/"
# dir = "~/Documents/MasterThesis/Data/CoMPASS/Messung2/"
dir = "~/Documents/MasterThesis/Data/CoMPASS/Stolberg/Messung1/"


# open root file
# file1 = ROOT.TFile.Open(dir+"data.root")
file1 = ROOT.TFile.Open(dir+"Messung1.root")

plotTimetags=0
plotWaveform=1
plotPSD=0

if plotWaveform:
    f1 = plt.figure("Waveform", figsize=(9,6))
    ax1 = f1.add_subplot(111)
    ax1.set_xlabel("t / ns")
    ax1.set_ylabel("ADC counts")

n=0
nMax=1e1
timeTags=[]
vEnergy=[]
vEnergyShort=[]
vPSD=[]
allTimeTags=[[],[]]
vWaveMean=[[],[]]
for i in xrange(8):
    for j in xrange(2):
        allTimeTags[j].append([])


# loop through events
for event in tqdm(file1.Events):

    # read CAEN device number
    device = event.device

    # read CAEN channel number
    channel = event.channel

    # compute PSD
    if plotPSD:
        vPSD.append((event.energy-event.energyShort)/event.energy)
        vEnergy.append(event.energy)
        vEnergyShort.append(event.energyShort)

    if plotTimetags:
        # measure time of first event
        if n==0:
            print "sample length:", len(np.array(event.samples, dtype="float64"))
            # t0=np.uint64(event.triggerTimeTag)/1e9
            t0=np.float64(event.triggerTimeTag)/1e9
            print "\nt0=", t0

        # timeTag = (np.uint64(event.triggerTimeTag)/1e9)#*2./1e6 - t0) # mus
        timeTag = (np.float64(event.triggerTimeTag)/1e9)#*2./1e6 - t0) # mus
        # timeTag = ((event.triggerTimeTagMuS))#*2./1e6 - t0) # ms
        # print timeTag
        allTimeTags[device][channel].append(timeTag)

    if plotWaveform:
        samples = np.array(event.samples, dtype="float64")
        ax1.plot(np.linspace(0, len(samples), len(samples))*2., samples)

    n+=1
    if n==nMax:
        print "breaking at n =", n
        break
####
# time tags
if plotTimetags:
    lw=0.3
    f2 = plt.figure("Time Tags", figsize=(9,6))
    ax2 = f2.add_subplot(111)
    for i in range(2):
        for j in range(8):
            # ax2.vlines(allTimeTags[i][j], i*8+j, i*8+j+1, lw=lw, color="black", alpha=1.)#, label=channelMapping[i][j])
            ax2.plot(range(len(allTimeTags[i][j])), allTimeTags[i][j], ".", ls="--", label=str(i)+" "+str(j))
            # ax2.hist(allTimeTags[i][j], bins=100, ls="-", fill=0, histtype='step', label=str(i)+" "+str(j))
    # leg = plt.legend(bbox_to_anchor=(0., -.2, 1., -0.95), loc=3, ncol=4, mode="expand", borderaxespad=0.)
    # for line in leg.get_lines():
    #     line.set_linewidth(4.)
    ax2.set_title("Time Tags")
    ax2.set_xlabel("index")
    ax2.set_ylabel("time in mu s")
    ax2.legend()
    # ax2.set_yticks(np.arange(0.5,15.5,step=1), tuple(channelMapping[0]+channelMapping[1]))
    plt.tight_layout()

    # compute trigger rate
    # avRate=0
    # nRate=0
    # vDeltaT=[]
    # for i in range(len(timeTags)-1):
    #     # print timeTags[i]
    #     if abs(timeTags[i+1]-timeTags[i]) > 0:
    #         avRate += 1./abs(timeTags[i+1]-timeTags[i])
    #         vDeltaT.append(abs(timeTags[i+1]-timeTags[i]))
    #         nRate+=1
    # avRate /= float(nRate)
    # print "min delta t:", min(vDeltaT)*1e3, "mu sec"
    #
    # print "average rate: ", avRate, "kHz"
###

if plotPSD:
    f3 = plt.figure("PSD", figsize=(9,6))
    ax3 = f3.add_subplot(111)
    # plt.hist2d(vEnergy, vPSD, bins=[300,300], range=[[0,1e4],[-1,1]], norm=LogNorm())
    plt.hist2d(vEnergy, vEnergyShort)#, bins=[300,300], range=[[0,1e4],[-1,1]], norm=LogNorm())
    plt.colorbar()
    ax3.set_title("PSD")
    ax3.set_xlabel("energy [a.u.]")
    ax3.set_ylabel("PSD")

plt.show()
