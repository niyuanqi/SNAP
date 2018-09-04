#################################################################
# Name:     SNMCcalib.py                                        #
# Author:   Yuan Qi Ni                                          #
# Date:     Apr. 26, 2018                                       #
# Function: Program performs Monte Carlo calculation of SN      #
#           parameters using SNooPy software to calculate       #
#           K corrections, fit for dm15/sBV. Gets redshift.     #
#           Update snpy light curve and ObjData.py first.       #
#################################################################

#essential modules
from snpy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#essential imports
from SNAP.Analysis.LCFitting import *
from SNAP.Analysis.Cosmology import *
from ObjData import *

N = 50 #number of Monte Carlo trials to run
plot = False #plot polynomial fits to light curves
#window of phillips relation in which intersection is estimated
sBV1 = 0.75
sBV2 = 0.90
#polynomial order with which to fit phillips section
npoly = 2
#dispersion intervals away from fit to consider point an outlier (for removal)
Nsig = 1.0
#window of redshifts in which SN redshift is estimated
z1 = 0.004
z2 = 0.007

print "Loading SN File"
s = get_sn(sn_file)

print "Loading Calibration File"
SNz, SNsBV, SNsBV_err, SNdm15, SNdm15_err, SNH, SNH_err, SNJ, SNJ_err, SNY, SNY_err, SNi, SNi_err, SNV, SNV_err, SNB, SNB_err = np.loadtxt(ph_file, comments='#', usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), unpack=True)
SNsBV = [SNsBV[SNB<90],SNsBV[SNV<90],SNsBV[SNi<90]]
SNsBV_err = [SNsBV_err[SNB<90],SNsBV_err[SNV<90],SNsBV_err[SNi<90]]
SNdm15 = [SNdm15[SNB<90],SNdm15[SNV<90],SNdm15[SNi<90]]
SNdm15_err = [SNdm15_err[SNB<90],SNdm15_err[SNV<90],SNdm15_err[SNi<90]]
SNM = [SNB[SNB<90],SNV[SNV<90],SNi[SNi<90]]
SNM_err = [SNB_err[SNB<90],SNV_err[SNV<90],SNi_err[SNi<90]]

#crop window of sBV
for i in range(len(band)):
    mask = np.logical_and(SNsBV[i]>sBV1,SNsBV[i]<sBV2)
    SNsBV_err[i] = SNsBV_err[i][mask]
    SNdm15[i] = SNdm15[i][mask]
    SNdm15_err[i] = SNdm15_err[i][mask]
    SNM[i] = SNM[i][mask]
    SNM_err[i] = SNM_err[i][mask]
    SNsBV[i] = SNsBV[i][mask]

#don't plot fit
s.replot = 0
#arrays to store Monte Carlo values
print "Performing Monte Carlo Calculations"
z = np.linspace(z1,z2,N)
#z = np.array([0.063])
Mdm = np.zeros([N,len(band)])
Mdm_err = np.zeros([N,len(band)])
dm, dm_err = np.zeros(N), np.zeros(N)
Mst = np.zeros([N,len(band)])
Mst_err = np.zeros([N,len(band)])
st, st_err = np.zeros(N), np.zeros(N)
#for each trial redshift
for i in range(N):
    print ""
    print "-----------------------------"
    print "Trial redshift", z[i]
    print "-----------------------------"
    s.z = z[i]
    #fit SN dm15
    s.choose_model("EBV_model2", stype="dm15")
    s.fit(snband)
    s.kcorr()
    print s.Tmax, s.e_Tmax
    #t, m, e, b = s.get_max(band, restframe=1, deredden=1, use_model=1)
    dm[i] = s.dm15
    dm_err[i] = s.e_dm15
    #get max magnitude of light curve from model fit
    for j in range(len(band)):
        #get absolute magnitude in rest frame
        M = s.data[band[j]].mag
        M_err = s.data[band[j]].e_mag
        kcorr = s.ks[band[j]]
        #deredden
        M = M - Coefs[j]*EBVgal
        t = s.data[band[j]].MJD
        M_abs = absMag(M, z[i], Kcorr=kcorr)
        t_abs = absTime(t, z[i])
        if j != 2:
            fit, fit_err, params, params_err = LCpolyFit(t_abs, M_abs, M_err, order=6, N=30, plot=plot)
        else:
            fit, fit_err, params, params_err = LCpolyFit(t_abs, M_abs, M_err, order=8, N=30, plot=plot)
        print params, params_err
        #maximum absolute magnitude
        Mdm[i][j] = params[1]
        Mdm_err[i][j] = params_err[1]
    print dm[i], dm_err[i]
    print Mdm[i][0], Mdm_err[i][0]

    #fit SN sBV
    s.choose_model("EBV_model2", stype="st")
    s.fit(snband)
    s.kcorr()
    print s.Tmax, s.e_Tmax
    #t, m, e, b = s.get_max(band, restframe=1, deredden=1, use_model=1)
    st[i] = s.st
    st_err[i] = s.e_st
    #get max magnitude of light curve from model fit
    for j in range(len(band)):
        #get absolute magnitude in rest frame
        M = s.data[band[j]].mag
        M_err = s.data[band[j]].e_mag
        kcorr = s.ks[band[j]]
        #deredden
        M = M - Coefs[j]*EBVgal
        t = s.data[band[j]].MJD
        M_abs = absMag(M, z[i], Kcorr=kcorr)
        t_abs = absTime(t, z[i])
        if j != 2:
            fit, fit_err, params, params_err = LCpolyFit(t_abs, M_abs, M_err, order=6, N=30, plot=plot)
        else:
            fit, fit_err, params, params_err = LCpolyFit(t_abs, M_abs, M_err, order=8, N=30, plot=plot)
        print params, params_err
        #maximum absolute magnitude
        Mst[i][j] = params[1]
        Mst_err[i][j] = params_err[1]
    print st[i], st_err[i]
    print Mst[i][0], Mst_err[i][0]

print ""
poly = []
disp = []
for i in range(len(band)):
    #fit polynomial to phillips relation, and plot
    print "Fitting phillips relation with polynomial."
    poly.append(np.polyfit(SNsBV[i], SNM[i], npoly))
    fitvals = np.polyval(poly[i], SNsBV[i])
    plt.errorbar(SNsBV[i],SNM[i], xerr=SNsBV_err[i],yerr=SNM_err[i], fmt='r+')
    plt.scatter(SNsBV[i], fitvals)
    disp.append(np.std(SNM[i]-fitvals))
    print disp[i]
    plt.ylim(-16,-20)
    plt.xlim(1.4,0.2)
    plt.show()

    #remove some phillips relation outliers, and plot again
    print "Kick out ouliers (>Nsig*dispersion from fit) and fit again."
    sBVcur = SNsBV[i][np.absolute(SNM[i]-fitvals)<Nsig*disp[i]]
    SNMcur = SNM[i][np.absolute(SNM[i]-fitvals)<Nsig*disp[i]]
    sBVcur_err = SNsBV_err[i][np.absolute(SNM[i]-fitvals)<Nsig*disp[i]]
    SNMcur_err = SNM_err[i][np.absolute(SNM[i]-fitvals)<Nsig*disp[i]]
    poly[i] = np.polyfit(sBVcur, SNMcur, npoly)
    fitvals = np.polyval(poly[i], sBVcur)
    plt.errorbar(sBVcur,SNMcur, xerr=sBVcur_err,yerr=SNMcur_err, fmt='r+')
    plt.scatter(sBVcur, fitvals)
    disp[i] = np.sqrt(np.std(SNMcur-fitvals)**2+np.mean(SNMcur_err**2))
    chi2dof = np.sum(np.square((SNMcur-fitvals)/np.sqrt(SNMcur_err**2 + sBVcur_err**2)))/(len(fitvals)-npoly-1)#
    corr = np.corrcoef(SNMcur, sBVcur)
    print "Fit quality"
    print "Dispersion:",disp[i]
    print "X2/dof:",chi2dof
    print "Pearson Correlation:",corr
    print "Number of phillips data points used:",len(sBVcur)
    plt.ylim(-16,-20)
    plt.xlim(1.4,0.2)
    plt.show()
print ""
print "Phillips Polynomial fits in B, V, I:",poly
print "Dispersion in B, V, I:",disp
print ""

zs = []
zs_err = []
for i in range(len(band)):
    print "Determine best fit redshift."
    curst = st[i]
    curst_err = st_err[i]
    curM = Mst.T[i]
    curM_err = Mst_err.T[i]

    fitvals = np.polyval(poly[i], curst)
    disps = np.sqrt(np.absolute(curM-fitvals)**2+curM_err**2)
    print z
    zs.append(z[np.argmin(disps)])
    zs_err.append(np.absolute(z[disps<disp[i]][0] - z[disps<disp[i]][-1]))
    print ""
w = 1/np.square(zs_err)
z_mean = np.sum(zs*w)/np.sum(w)
z_err = np.sqrt(1/np.sum(w))
print ""
print "---------------------------------------"
print "Best fit redshifts by band"
print zs, zs_err
print "Best fit mean z:", z_mean, z_err
print "---------------------------------------"
print ""

#colorbar parameterized by redshift
norm = mpl.colors.Normalize(vmin=z[0],vmax=z[-1])
cmap = mpl.cm.jet
smap = mpl.cm.ScalarMappable(cmap=cmap,norm=norm)
smap.set_array([])
scatter_kwargs = {"zorder":100}
error_kwargs = {"zorder":0}

#generate plots
f, ax = plt.subplots(3,2)
ax[1][0].set_ylabel(r"$M_{\lambda}-A_{\lambda}$", fontsize = 14)
ax[2][1].set_xlabel(r"$\Delta M_{15}(B)$", fontsize = 14)
ax[2][0].set_xlabel(r"$s_{BV}$", fontsize = 16)
ax[0][0].text(0.3,-19.2,"B", fontsize = 14, fontstyle='italic', fontweight='bold')
ax[1][0].text(0.3,-19.2,"V", fontsize = 14, fontstyle='italic', fontweight='bold')
ax[2][0].text(0.3,-19.2,"I", fontsize = 14, fontstyle='italic', fontweight='bold')
ax[0][1].text(1.95,-19.2,"B", fontsize = 14, fontstyle='italic', fontweight='bold')
ax[1][1].text(1.95,-19.2,"V", fontsize = 14, fontstyle='italic', fontweight='bold')
ax[2][1].text(1.95,-19.2,"I", fontsize = 14, fontstyle='italic', fontweight='bold')

print "Plotting Phillips Data"
    
for i in range(len(band)):
    ax[i][1].errorbar(SNdm15[i],SNM[i], xerr=SNdm15_err[i],yerr=SNM_err[i], fmt='r+')
    ax[i][1].errorbar(dm,Mdm.T[i], xerr=dm_err,yerr=Mdm_err.T[i], fmt='g+')
    ax[i][1].scatter(dm,Mdm.T[i],c=smap.to_rgba(z), **scatter_kwargs)

    ax[i][0].errorbar(SNsBV[i],SNM[i], xerr=SNsBV_err[i],yerr=SNM_err[i], fmt='r+')
    ax[i][0].errorbar(st,Mst.T[i], xerr=st_err,yerr=Mst_err.T[i], fmt='g+')
    ax[i][0].scatter(st,Mst.T[i],c=smap.to_rgba(z), **scatter_kwargs)

    ax[i][0].yaxis.set_ticks([-19,-18,-17])
    ax[i][1].yaxis.set_ticks([])
    ax[i][0].xaxis.set_ticks([0.4,0.8,1.2])
    ax[i][1].xaxis.set_ticks([1.0,1.5,2.0])
    ax[i][0].tick_params(labelsize=12)
    ax[i][1].tick_params(labelsize=12)
    ax[i][0].set_ylim(-16,-20)
    ax[i][1].set_ylim(-16,-20)
    ax[i][0].set_xlim(1.4,0.2)
    ax[i][1].set_xlim(0.6,2.1)

f.subplots_adjust(wspace=0)
f.subplots_adjust(hspace=0)

cbar_ax = f.add_axes([0.91, 0.15, 0.015, 0.7])
cbar_ax.tick_params(labelsize=12)
cbar_ax.set_yticks([0.035,0.045,0.055,0.065,0.075])
cbar_ax.set_yticklabels(z,rotation=90)
cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap,norm=norm, orientation='vertical')
cb1.set_label('Trial Redshifts', fontsize = 14)

plt.show()

