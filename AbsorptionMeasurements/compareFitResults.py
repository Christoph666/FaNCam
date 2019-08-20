import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

mapping = ("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4")


dir = "/home/christoph/Documents/MasterThesis/Data/"
folder = sys.argv[-2]
measurement = sys.argv[-1]
dir = dir + folder + "/"
thisfile = dir + measurement + "/Fit_Results/Results_off.txt";
thisfile_sum = dir + measurement + "/Fit_Results/Results_sum_off.txt";
print "reading:", thisfile

data = pd.read_csv(thisfile, sep="\t", header=0)
data_sum = pd.read_csv(thisfile_sum, sep="\t", header=0)
print data
print ""
print data_sum

fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(range(len(data.xn)), data.xn, yerr=data.xn_err, label="Neutrons", marker="o", color="green")
ax.errorbar(range(len(data.xg)), data.xg, yerr=data.xg_err, label="Gammas", marker="o", color="red")
ax.hlines(data_sum.xn, 0, 15, color="green", linestyle="--", linewidth=2, label="sum result")
ax.hlines(data_sum.xn - data_sum.xn_err, 0, 15, color="green", linestyle="--", linewidth=0.5, label=None)
ax.hlines(data_sum.xn + data_sum.xn_err, 0, 15, color="green", linestyle="--", linewidth=0.5, label=None)
ax.fill_between(range(16), data_sum.xn - data_sum.xn_err, data_sum.xn + data_sum.xn_err, facecolor="green", alpha=0.2, label="sum error")

ax.hlines(data_sum.xg, 0, 15, color="red", linestyle="--", linewidth=2, label="sum result")
ax.hlines(data_sum.xg - data_sum.xg_err, 0, 15, color="red", linestyle="--", linewidth=0.5, label=None)
ax.hlines(data_sum.xg + data_sum.xg_err, 0, 15, color="red", linestyle="--", linewidth=0.5, label=None)
ax.fill_between(range(16), data_sum.xg - data_sum.xg_err, data_sum.xg + data_sum.xg_err, facecolor="red", alpha=0.2, label="sum error")
ax.legend(loc="best")

ax.set_title(measurement)
ax.set_xlabel("Pixel")
ax.set_ylabel(r"$x$ / cm")
ax.set_xticks(range(16))
ax.set_xticklabels(mapping)

plt.show()
