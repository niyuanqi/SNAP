#essential modules
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'axes.linewidth': 1.})
#plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams['text.latex.preamble']=[r"\usepackage{amssymb}"]

#essential files
from SNAP.Analysis.ShockMod import *
from SNAP.Analysis.Cosmology import *
from SNAP.Analysis.Ni56Mod import *

#example SN parameters
m_ej = 1.0
m_c = 1.0/1.4
e_ej = 0.8
lim = -8

#shock breakout model
x = [0, 10]
t0 = x[0] #days
R8 = x[1] #10^8 cm
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_sh = np.array([ShockCoolingFit(ti, wave_0[1], 0, 0, m_ej/1.4, R8, 0.1, t0) for ti in tT])
#predictions
tcT = tT[tT > 0]
McT = Flux_toMag('B',B_sh[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'g-', zorder=8)
plt.fill_between(tcT-t0, McT, -8*np.ones(len(tcT)), color='g', alpha=0.5)
plt.text(0.3, -8.15, "SBO", color='w', fontsize=18)

#CSM-disk interaction model
x = [0, 1.0, 0.0001]
t0 = x[0] 
Mext = x[1] #0.01 Msun
Rext = x[2] #10^13 cm
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_csm = np.array([CSMFit(ti, wave_0[1], 0, 0, m_ej, e_ej, Mext, Rext, t0) for ti in tT])
#predictions
tcT = tT[tT > 0]
McT1 = Flux_toMag('B',B_csm[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'r--', zorder=8)
#CSM-disk interaction model
x = [0, 1.0, 0.1]
t0 = x[0] #days
Mext = x[1] #0.01 Msun
Rext = x[2] #10^13 cm
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_csm = np.array([CSMFit(ti, wave_0[1], 0, 0, m_ej, e_ej, Mext, Rext, t0) for ti in tT])
#predictions
tcT = tT[tT > 0]
McT2 = Flux_toMag('B',B_csm[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'r--', zorder=8)
plt.fill_between(tcT-t0, McT2, McT1, color='r', alpha=0.4)
plt.text(1.15, -8.5, "CSM", color='w', fontsize=22)

#WD/He-star Companion interaction model
x = [0, 0.0001, 0]
t0 = x[0] #days
a13 = x[1] #10^13 cm
theta = x[2] %180
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_kas = np.array([KasenFit(ti, a13, 1.0, wave_0[1], 0, m_c, e_ej, 0, t0)
                  for ti in tT]) * Kasen_isocorr(theta)
#predictions
tcT = tT[tT > 0]
McT1 = Flux_toMag('B',B_kas[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'b:', zorder=8)
#Companion interaction model
x = [0, 0.0001, 179.9]
t0 = x[0] #days
a13 = x[1] #10^13 cm
theta = x[2] %180
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_kas = np.array([KasenFit(ti, a13, 1.0, wave_0[1], 0, m_c, e_ej, 0, t0)
                  for ti in tT]) * Kasen_isocorr(theta)
#predictions
tcT = tT[tT > 0]
McT2 = Flux_toMag('B',B_kas[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'c:', zorder=8)
plt.fill_between(tcT-t0, McT2, McT1, color='b', alpha=0.45, label=r"WD/He$\bigstar$")

#1MS Companion interaction model
x = [0, 0.015, 0]
t0 = x[0] #days
a13 = x[1] #10^13 cm
theta = x[2] %180
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_kas = np.array([KasenFit(ti, a13, 1.0, wave_0[1], 0, m_c, e_ej, 0, t0)
                  for ti in tT]) * Kasen_isocorr(theta)
#predictions
tcT = tT[tT > 0]
McT1 = Flux_toMag('B',B_kas[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'b:', zorder=8)
#Companion interaction model
x = [0, 0.015, 179.9]
t0 = x[0] #days
a13 = x[1] #10^13 cm
theta = x[2] %180
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_kas = np.array([KasenFit(ti, a13, 1.0, wave_0[1], 0, m_c, e_ej, 0, t0)
                  for ti in tT]) * Kasen_isocorr(theta)
#predictions
tcT = tT[tT > 0]
McT2 = Flux_toMag('B',B_kas[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'c:', zorder=8)
plt.fill_between(tcT-t0, McT2, McT1, color='c', alpha=0.5)
plt.text(2, -10.8, "1M$_{\odot}$ MS", color='w', fontsize=22)

#1RG Companion interaction model
x = [0, 2, 0]
t0 = x[0] #days
a13 = x[1] #10^13 cm
theta = x[2] %180
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_kas = np.array([KasenFit(ti, a13, 1.0, wave_0[1], 0, m_c, e_ej, 0, t0)
                  for ti in tT]) * Kasen_isocorr(theta)
#predictions
tcT = tT[tT > 0]
McT1 = Flux_toMag('B',B_kas[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'b:', zorder=8)
#Companion interaction model
x = [0, 2, 179.9]
t0 = x[0] #days
a13 = x[1] #10^13 cm
theta = x[2] %180
tT = np.linspace(x[0]+0.0001, 8, 3000)
B_kas = np.array([KasenFit(ti, a13, 1.0, wave_0[1], 0, m_c, e_ej, 0, t0)
                  for ti in tT]) * Kasen_isocorr(theta)
#predictions
tcT = tT[tT > 0]
McT2 = Flux_toMag('B',B_kas[tT>0]*1.e-6)
#plt.plot(tcT-t0, McT, 'c:', zorder=8)
plt.fill_between(tcT-t0, McT2, McT1, color='m', alpha=0.3)
plt.text(2, -15.5, "1M$_{\odot}$ RG", color='w', fontsize=22)

#Ni56 logistic model
td, Ld, Mc, Ek, beta, x_2 = 20.670, 1.150e+43, 0.585, 0.642, 8.442, 0.317
t_pred = np.arange(td/1000,10,0.01)
x = [0.013, 0.25]
x_s, a_s = x[0], x[1]
#estimate optical depth c/v
popt = [20.22,-0.40] #i-band subtracted background star
#perr = [0.48, 0.68]
#power 2 decay
def vpow(t, a, b):
    return a*np.power(t-b, -0.22)
#optical depth by band
def tau_pow(t, K=1):
    return K*300.0/vpow((t), *popt)
K = 0.05523718
taus_pred = tau_pow(t_pred, K)
L_o = predNi56mod(t_pred, wave_0[1], 0, 0, taus_pred, 
                  td, Ld, Mc, Ek, beta, x_2)
M_o = Flux_toMag('B', L_o*1e-6)
#plt.plot(t_pred, M_o, color='k', linestyle=':')
#Ni56 clump model
L_pred = predShallowNimod(t_pred, wave_0[1], 0, 0, taus_pred, 
                          td, Ld, Mc, Ek, beta, x_2, x_s, a_s)
M_pred = Flux_toMag('B', L_pred*1e-6)
plt.fill_between(t_pred, M_pred, M_o, color='k', alpha=0.2, label="$^{56}$Ni Clump")

#KMTNet detection limits
plt.text(2.6, -8.1, "  8 Mpc")
plt.plot([0, 3], [-10, -10], color='k', linestyle='--')
plt.text(2.6, -10.1, "20 Mpc")
plt.plot([0, 3], [-12, -12], color='k', linestyle='--')
plt.text(2.6, -12.1, "50 Mpc")
plt.plot([0, 3], [-13, -13], color='k', linestyle='--')
plt.text(2.6, -13.1, "75 Mpc")
    
#plot axes
plt.ylabel("B-band Absolute Magnitude", fontsize=14)
plt.xlabel("Days since explosion", fontsize=14)
plt.ylim([lim, -19])
plt.xlim([0, 3])
plt.legend(loc='upper left', fontsize=14)
plt.tight_layout()
plt.show()
