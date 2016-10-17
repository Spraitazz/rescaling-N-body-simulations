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


#data = np.loadtxt("./ZA_pk_" + str(startind) + ".dat")
#data_4fold = np.loadtxt("./ZA_pk_4fold_" + str(startind) + ".dat")
#data = np.vstack((data[0:16,:], data_4fold[3:,:]))



data_debiased = np.loadtxt("./ZA_debiased_pk_" + str(startind) + ".dat")
data_debiased_4fold = np.loadtxt("./ZA_debiased_pk_4fold_" + str(startind) + ".dat")
data_debiased = np.vstack((data_debiased[0:16,:], data_debiased_4fold[3:,:]))

ks = data_debiased[:,0]
#data_mono = data[:,1]
#data_quadro = data[:,2]

data_mono_debiased = data_debiased[:,1]

for i in xrange(1, simno):
	#this_data = np.loadtxt("./ZA_pk_" + str(i + startind) + ".dat")
	this_data_debiased = np.loadtxt("./ZA_debiased_pk_" + str(i + startind) + ".dat")
	
	#this_data_4fold = np.loadtxt("./ZA_pk_4fold_" + str(i + startind) + ".dat")
	#this_data = np.vstack((this_data[0:16,:], this_data_4fold[3:,:]))
	
	this_data_debiased_4fold = np.loadtxt("./ZA_debiased_pk_4fold_" + str(i + startind) + ".dat")
	this_data_debiased = np.vstack((this_data_debiased[0:16,:], this_data_debiased_4fold[3:,:]))
	
	#data_mono = np.vstack((this_data[:,1], data_mono))
	#data_quadro = np.vstack((this_data[:,2], data_quadro))	
	
	data_mono_debiased = np.vstack((this_data_debiased[:,1], data_mono_debiased))
	
	
data_model_debiased = np.loadtxt("./model_pk_debiased.dat")
data_model = np.loadtxt("./model_pk.dat")

#mean_mono = np.mean(data_mono, axis = 0)
#mean_quadro = np.mean(data_quadro, axis = 0)

mean_mono_debiased = np.mean(data_mono_debiased, axis = 0)

#var_mono = np.var(data_mono, axis = 0)
#var_quadro = np.var(data_quadro, axis = 0)

var_mono_debiased = np.var(data_mono_debiased, axis = 0)

#err_on_mean_mono = np.sqrt(var_mono) / np.sqrt(simno)
#err_on_mean_quadro = np.sqrt(var_quadro) / np.sqrt(simno)

err_on_mean_mono_debiased = np.sqrt(var_mono_debiased) / np.sqrt(simno)

#END data input, error calc----------------------------------------

ax = pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")

#pl.errorbar(ks, mean_mono, yerr=err_on_mean_mono, label="monopole")
#pl.plot(ks, mean_mono, label="monopole")
pl.errorbar(ks, mean_mono_debiased, yerr=err_on_mean_mono_debiased, label="debiased monopole")
pl.plot(data_model[:,0], data_model[:,1], label="model debiased (b average))")
pl.plot(data_model_debiased[:,0], data_model_debiased[:,1], label="model debiased (b_eff S&T HMF)")

pl.legend(loc="best")

pl.ylim(2e2,9e2)
pl.xlim(0.03, 0.4)

pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")
pl.xlabel(r"$k, \ [h \, \rm{Mpc} ^{-1}]$")


pl.savefig("./twohalo_bias.pdf", bbox_inches='tight')


