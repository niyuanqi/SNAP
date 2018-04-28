#################################################################
# Name:     plotkasen.py                                        #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 27, 2017                                       #
# Function: Program plots output of KasenMC.py                  #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt

#essential files
from ObjData import *

a13s, out1, out2, out2, out3, out4, out5 = np.loadtxt(kasfile, unpack=True)
outangles = [out1, out2, out3, out4, out5]
style = ['k-', 'k--']
for n, conf in enumerate(confs):
    #Plot angles ruled out for each a13 at this conf
    plt.plot(a13s, outangles[n], style[n])

plt.xlim(0,10.0)
plt.ylim(0,185)
#plot positions of 1RG, 6MS, 2MS
plt.plot([0.05,0.05], [175,185], 'k:', linewidth=1)
plt.plot([0.2,0.2], [175,185], 'k:', linewidth=1)
plt.plot([2.0,2.0], [175,185], 'k:', linewidth=1)
plt.ylabel("Unacceptable viewing angles (deg)", fontsize=16)
plt.xlabel("Separation Distance ($10^{13}$ cm)", fontsize=16)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.show()
