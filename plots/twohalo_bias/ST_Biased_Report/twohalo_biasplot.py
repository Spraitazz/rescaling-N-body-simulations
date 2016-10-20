import matplotlib.pylab        as plt
import matplotlib.pyplot as pl
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import ScalarFormatter
from   matplotlib.ticker       import FixedFormatter
import pylab                   as pl
import numpy                   as np
import math, os, sys
import glob, pickle
from scipy import interpolate

formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3,3))

fig_width_pt  = 504.0 # Get this from LaTex using \the\columnwidth                                               
                                                                                                                                                                                                                                            
inches_per_pt = 1.0/72.27

golden_mean  = 1.0
golden_mean *= (np.sqrt(5)-1.0)/2.0

fig_width  = fig_width_pt*inches_per_pt # width in inches                                              
                                                                                                                                                                                                                                            
fig_height = fig_width*golden_mean # height in inches                                          
                                                                                                                                                                                                                                            
fig_size   = [fig_width, fig_height]

params     = {'axes.labelsize':12,
              'font.size':12,
              'legend.fontsize':9,
              'xtick.labelsize':10.0,
              'ytick.labelsize':10.0,
              'figure.figsize':fig_size,
              'font.family': 'serif'}
              
pl.rcParams.update(params)
pl.clf()
pl.figure()

fig = pl.figure()

axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)

axes.xaxis.set_major_formatter(formatter)
axes.yaxis.set_major_formatter(formatter)

#END FANCY PLOT ----------------------------------------------


startind = int(sys.argv[1])


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

simno = 0
while os.path.isfile("./ZA_debiased_pk_" + str(startind + simno) + ".dat"):
	simno+=1

rows = file_len("./ZA_debiased_pk_" + str(startind)+ ".dat")

data_debiased = np.loadtxt("./ZA_debiased_pk_" + str(startind) + ".dat")
data_debiased_mmin = np.loadtxt("./ZA_debiased_1e13_pk_" + str(startind) + ".dat")
data_debiased_mmax = np.loadtxt("./ZA_debiased_1e14_pk_" + str(startind) + ".dat")

kmin = 0.05
kmax = 0.3
minInd = 0
maxInd = 14
ks = data_debiased[minInd:maxInd,0]

data_mono_debiased = data_debiased[minInd:maxInd,1]
data_mono_debiased_mmin = data_debiased_mmin[minInd:maxInd,1]
data_mono_debiased_mmax = data_debiased_mmax[minInd:maxInd,1]

for i in xrange(1, simno):

	this_data_debiased = np.loadtxt("./ZA_debiased_pk_" + str(i + startind) + ".dat")
	this_data_debiased_mmin = np.loadtxt("./ZA_debiased_1e13_pk_" + str(i+startind) + ".dat")
	this_data_debiased_mmax = np.loadtxt("./ZA_debiased_1e14_pk_" + str(i+startind) + ".dat")
	
	data_mono_debiased = np.vstack((this_data_debiased[minInd:maxInd,1], data_mono_debiased))
	data_mono_debiased_mmin = np.vstack((this_data_debiased_mmin[minInd:maxInd,1], data_mono_debiased_mmin))
	data_mono_debiased_mmax = np.vstack((this_data_debiased_mmax[minInd:maxInd,1], data_mono_debiased_mmax))
	
	
data_model_bveff = np.loadtxt("./model_pk_beff.dat")
data_model_bvmin = np.loadtxt("./model_pk_bvmin.dat")
data_model_bvmax = np.loadtxt("./model_pk_bvmax.dat")

mean_mono_debiased = np.mean(data_mono_debiased, axis = 0)
mean_mono_debiased_mmin = np.mean(data_mono_debiased_mmin, axis = 0)
mean_mono_debiased_mmax = np.mean(data_mono_debiased_mmax, axis = 0)

var_mono_debiased = np.var(data_mono_debiased, axis = 0)
var_mono_debiased_mmin = np.var(data_mono_debiased_mmin, axis = 0)
var_mono_debiased_mmax = np.var(data_mono_debiased_mmax, axis = 0)

err_on_mean_mono_debiased = np.sqrt(var_mono_debiased) / np.sqrt(simno)
err_on_mean_mono_debiased_mmin = np.sqrt(var_mono_debiased_mmin) / np.sqrt(simno)
err_on_mean_mono_debiased_mmax = np.sqrt(var_mono_debiased_mmax) / np.sqrt(simno)

#END data input, error calc----------------------------------------

ax = pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")

pl.plot(data_model_bveff[:,0], data_model_bveff[:,1], color="Aqua", label=r"model $P_{0}$ (all haloes)")
pl.errorbar(ks, mean_mono_debiased, yerr=err_on_mean_mono_debiased, linestyle="None", color="Blue", label=r"$P_{0}$ (all haloes)")

pl.plot(data_model_bvmin[:,0], data_model_bvmin[:,1], color="DeepPink", label=r"model $P_{0}$ ($10^{12} M_{\odot}$ haloes)")
pl.errorbar(ks, mean_mono_debiased_mmin, yerr=err_on_mean_mono_debiased_mmin, linestyle="None", color="Red", label=r"$P_{0}$ ($10^{12} M_{\odot}$ haloes)")

pl.plot(data_model_bvmax[:,0], data_model_bvmax[:,1], color="GreenYellow", label=r"model $P_{0}$ ($10^{13} M_{\odot}$ haloes)")
pl.errorbar(ks, mean_mono_debiased_mmax, yerr=err_on_mean_mono_debiased_mmax, linestyle="None", color="Green", label=r"$P_{0}$ ($10^{13} M_{\odot}$ haloes)")


pl.legend(loc="best", ncol=2)

pl.ylim(2e2,3e3)
pl.xlim(0.05, 0.2)

pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")
pl.xlabel(r"$k, \ [h \, \rm{Mpc} ^{-1}]$")


pl.savefig("./twohalo_bias.pdf", bbox_inches='tight')


