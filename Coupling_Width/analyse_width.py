import numpy as np
import matplotlib.pyplot as plt

vWidth_small = [0.55, 0.4, 0.51, 0.6, 0.56, 0.65, 0.65, 0.56]
vWidth_large = [0.59, 0.59, 0.62, 0.67]

f1 = plt.figure('f1')
ax = f1.add_subplot(111)
plt.title('Optical Coupling Width')
plt.xlabel('Width / mm')
plt.ylabel('Counts')
nBins=6
plt.hist(vWidth_small, bins=nBins, range=[0.4, 0.7], ls='-', lw=2, alpha=0.5, color='orange', edgecolor='orange', hatch="/", label='small')
plt.hist(vWidth_large, bins=nBins, range=[0.4, 0.7], ls='-', lw=2, alpha=0.5, color='blue', edgecolor='blue', hatch="\\", label='large')
plt.legend(loc='best')
plt.text(0.41, 2.7, 'mean small: {:.2f}'.format(np.mean(vWidth_small)), color='orange')
plt.text(0.41, 2.5, 'mean large: {:.2f}'.format(np.mean(vWidth_large)), color='blue')
f1.savefig('width_plot.pdf')
# plt.show()

print'mean small: ', np.mean(vWidth_small)
print'mean large: ', np.mean(vWidth_large)
