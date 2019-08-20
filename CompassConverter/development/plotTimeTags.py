################################################################################
### CAEN CoMPASS binary file to ROOT file converter                          ###
### Version 1.0                                                              ###
### AUTHOR: Christoph Guenther                                               ###
################################################################################

import sys

import matplotlib.pyplot as plt
import numpy as np

from ROOT import TFile, TTree
from array import array

# file directories
# Cs137
# dir_in = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/run_Cs137/UNFILTERED/";
# file_in = "run_Cs137"

# AmBe
# dir_in = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/run_AmBe/UNFILTERED/";
# file_in = "run_AmBe"

# testrun
dir_in = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/run1/UNFILTERED/";
file_in = "run1"

# Stolberg
# file_in = "Messung2"
# dir_in = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/"+file_in+"/UNFILTERED/";

# test
# dir_in = "/home/christoph/Documents/MasterThesis/Data/CoMPASS/test/DAQ/run_bin/UNFILTERED/";
# file_in = "run_bin"

# data storage
vTT = [[],[]]

saveMusTimeTag=0
plotWaveform=1

firstBlock = True # check if first block is read (need some data for branching)

# device names (serial number)
vDEV=["50", "58"]

breakAfterN = 1e2
eventCounter=[0,0]
waveforms = []

if plotWaveform:
    fWave = plt.figure("fWave", figsize=(7,6))
    axWave = fWave.add_subplot(111)
    iWaveform=3


# loop over devices and channels
for iDEV in range(0,1):
    DEV = vDEV[iDEV]
    print "**********************"
    print "Reading data of device", DEV, "\b."

    for iCHN in range(0,2):
        filename = str(iCHN)+"@N6730B #3-11-"+DEV+"_Data_"+file_in+".bin";
        print "Reading file:", filename

        try:
            # input binary file
            f_in = open(dir_in+filename, "rb")
            # print "Openend file."

            # get number of samples in a waveform and bytes per block
            bytes = f_in.read(22)
            b = bytearray()
            b.extend(bytes)
            # print "bytes:", bytes
            # print "b:", b
            # print "len b:", len(b)

            num_of_samples = 0
            for i in range(4):
                # print "b[", 18+i, "] =", b[18+i]
                num_of_samples += b[18+i]*2**(i*8)
            # 16bit integers take 2 bytes, 32bit take 4 bytes etc...
            bytes_per_block = 2+2+8+2+2+2+4+num_of_samples*2

            # go back to start of file
            f_in.seek(0, 0)
            # print "Seeking to 0."

            print "Samples in one waveform:", num_of_samples
            print "Bytes per block:", bytes_per_block
            blockCount=0


            while(1):

                sys.stdout.flush()
                sys.stdout.write("block: "+str(blockCount)+"\r")

                # read one block of data
                bytes = f_in.read(bytes_per_block)
                b = bytearray()
                b.extend(bytes)

                # print "len b:", len(b)
                # for i in range(len(b)):
                #     print b[i]

                # check if end of file is reached
                if len(b) < bytes_per_block:
                    print "Detected end of file."
                    break

                # interpret bytes
                board = b[0] + b[1]*2**8
                channel = b[2] + b[3]*2**8
                timeStamp = 0
                for i in range(8):
                    timeStamp += b[4+i]*2**(i*8)
                if saveMusTimeTag:
                    timeStampMuS = timeStamp*2./1e3
                energy = b[12] + b[13]*2**8
                energyShort = b[14] + b[15]*2**8
                flags = b[16] + b[17]*2**8
                nSamples = 0
                for i in range(4):
                    nSamples += b[18+i]*2**(i*8)
                vSamples=[]
                for i in range(num_of_samples):
                    sample = b[22+i*2] + b[23+i*2]*2**8
                    vSamples.append(sample)

                # print "\nreading block", blockCount
                # print "\t- board:", board
                # print "\t- channel:", channel
                # print "\t- timeStamp:", timeStamp
                # print "\t- energy:", energy
                # print "\t- energyShort:", energyShort
                # print "\t- flags:", flags
                # print "\t- nSamples:", nSamples
                # for i in range(10):
                #     print "\t- sample", i, ":", vSamples[i]

                # vTT[iCHN].append(timeStamp)
                if iDEV==0 and iCHN==0:
                    vTT[0].append(timeStamp)
                if iDEV==0 and iCHN==1:
                    vTT[1].append(timeStamp)
                if plotWaveform:
                    if eventCounter[iCHN]==iWaveform:
                        # axWave.plot(np.linspace(0, len(vSamples)*2., len(vSamples)), vSamples, label=str(iCHN))
                        waveforms.append(vSamples)
                eventCounter[iCHN]+=1


                # print timeStamp
                # print timeStampMuS

                if breakAfterN > 0:
                    if blockCount==breakAfterN:
                        break
                blockCount+=1

            print "Found", blockCount, "blocks."
            print "Finished reading file", filename, "\b."
            print ""

            f_in.close()

        except:
            print "File", filename, "not found or empty.\n"

if plotWaveform:
    axWave.plot(np.linspace(0, len(waveforms[0]), len(waveforms[0]))*2, waveforms[0], label=str(0))
    axWave.plot(np.linspace(0, len(waveforms[1]), len(waveforms[1]))*2.-np.ones(len(waveforms[1]))*(vTT[1][iWaveform]-vTT[0][iWaveform])/1e6, waveforms[1], label=str(1))
    axWave.set_xlabel("t / ns")
    axWave.set_ylabel("ADC")
    axWave.legend()

f1 = plt.figure("f1", figsize=(7,6))
ax1 = f1.add_subplot(111)
ax1.plot(range(len(vTT[0])), vTT[0], ".", color="blue")
ax1.plot(range(len(vTT[1])), vTT[1], ".", color="orange")
# ax1.plot(range(len(vTT[1][1:])), np.array(vTT[0][1:])-np.array(vTT[1][1:]), ".", color="orange")

f2 = plt.figure("f2", figsize=(7,6))
ax2 = f2.add_subplot(111)
ax2.hist((np.array(vTT[1][1:])-np.array(vTT[0][1:]))/1e6, bins=20)
ax2.set_xlabel("t / ns")
ax2.set_ylabel("Counts")

plt.show()
