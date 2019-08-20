import numpy as np
import matplotlib.pyplot as plt

x12_Lead = 6.6
x12_Paraffin = 6.8

def relAbsorption(x, x12):
    return np.exp(-x/x12*np.log(2.)) - 1.

print "Absorption after 5cm Paraffin:", relAbsorption(5, x12_Paraffin)
print "Absorption after 6mm Lead:", relAbsorption(0.6, x12_Lead)
print "Total Absorption:", np.exp(-5./x12_Paraffin*np.log(2.)) * np.exp(-0.6/x12_Lead*np.log(2.)) - 1.

# x = np.linspace(0, 30., 100)
# plt.plot(x, relAbsorption(x, x12_Paraffin))
# plt.show()
