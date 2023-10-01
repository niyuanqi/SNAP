#################################################################
# Name:     plotkasen.py                                        #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 27, 2017                                       #
# Function: Program plots output of KasenMC.py                  #
#################################################################

#essential modules
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

#essential files
from ObjData import *

outs = np.loadtxt(kasfile, unpack=True)
a13s = outs[0]
outangles = outs[1:]
for n, sn in enumerate(limSNs):
    if sn in plot:
        print "Sig and conf:", sn, norm.cdf(sn)
        #Plot angles ruled out for each a13 at this conf
        plt.plot(a13s, outangles[n], style[plot[sn]])

plt.xlim(0,10.0)
plt.ylim(0,185)
#plot positions of 1RG, 6MS, 2MS
plt.plot([0.05,0.05], [175,185], 'k:', linewidth=1)
plt.plot([0.2,0.2], [175,185], 'k:', linewidth=1)
plt.plot([2.0,2.0], [175,185], 'k:', linewidth=1)
plt.ylabel("Lower Limit of Acceptable \nViewing Angles (Degrees)", fontsize=16)
plt.xlabel("Separation Distance ($10^{13}$ cm)", fontsize=16)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.show()
