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


startind = 0

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

simno = 0
while os.path.isfile("./ZA_RSD_pk_" + str(startind + simno) + ".dat"):
	simno+=1

rows = file_len("./ZA_RSD_pk_" + str(startind)+ ".dat")

data = np.loadtxt("./ZA_RSD_pk_" + str(startind) + ".dat")
ks = data[0:26,0]
data_mono = data[0:26,1]
data_quadro = abs(data[0:26,2])
data_hex = abs(data[0:26,3])

for i in xrange(1, simno):
	#print i
	this_data = np.loadtxt("./ZA_RSD_pk_" + str(i + startind) + ".dat")
	data_mono = np.vstack((this_data[0:26,1], data_mono))
	data_quadro = np.vstack((abs(this_data[0:26,2]), data_quadro))
	data_hex = np.vstack((abs(this_data[0:26,3]), data_hex))	
	
	
data_model = np.loadtxt("./model_pk.dat")

#kmin = ks[0]
#kmax = ks[rows-1]

mean_mono = np.mean(data_mono, axis = 0)
mean_quadro = np.mean(data_quadro, axis = 0)
mean_hex = np.mean(data_hex, axis = 0)

var_mono = np.var(data_mono, axis = 0)
var_quadro = np.var(data_quadro, axis = 0)
var_hex = np.var(data_hex, axis = 0)

err_on_mean_mono = np.sqrt(var_mono) / np.sqrt(simno)
err_on_mean_quadro = np.sqrt(var_quadro) / np.sqrt(simno)
err_on_mean_hex = np.sqrt(var_hex) / np.sqrt(simno)
print err_on_mean_hex[0]
print err_on_mean_hex[16]

#END data input, error calc----------------------------------------

mono_over_hex = np.divide(mean_mono, mean_hex)

frac_err_mono = np.divide(err_on_mean_mono, mean_mono)
frac_err_hex = np.divide(err_on_mean_hex, mean_hex)

err_on_mean_mono_over_hex = np.multiply(np.add(frac_err_mono, frac_err_hex), mono_over_hex)

f_GR = pow(0.24, 0.545)
ratio = f_GR*f_GR*8.0/35.0
ratio = (1.0 + (2.0/3.0)*f_GR + (1.0/5.0)*f_GR*f_GR)/((4.0/3.0)*f_GR + (4.0/7.0)*f_GR*f_GR)
ratio = (1.0 + (2.0/3.0)*f_GR + (1.0/5.0)*f_GR*f_GR)/((8.0/35.0)*f_GR*f_GR)
print ratio

f, axarr = pl.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
axarr[0].set_xscale("log")
axarr[0].set_yscale("log")

axarr[0].plot(data_model[:,0], data_model[:,1], color="Aqua", label=r"model $P_{0}(k)$")
axarr[0].plot(data_model[:,0], data_model[:,2], color="Yellow", label=r"model $P_{2}(k)$")
axarr[0].plot(data_model[:,0], data_model[:,3], label=r"model $P_{4}(k)$")

axarr[0].errorbar(ks, mean_mono, yerr=err_on_mean_mono, linestyle="None", color="Blue", label=r"measured $P_{0}(k)$")
axarr[0].errorbar(ks, mean_quadro, yerr=err_on_mean_quadro, linestyle="None", label=r"measured $P_{2}(k)$")
axarr[0].errorbar(ks, mean_hex, yerr=err_on_mean_hex, linestyle="None", label=r"measured $P_{4}(k)$")


axarr[0].legend(loc="best", ncol=2)
axarr[0].set_ylim(1.0,2e3)
#axarr[0].set_xlim(kmin, kmax)
axarr[0].set_ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

axarr[1].set_xscale("log")
axarr[1].errorbar(ks, mono_over_hex, yerr=err_on_mean_mono_over_hex, label=r"$\left(\frac{P_{0}(k)}{P_{4}(k)}\right)", fmt="o", ms = 0.5)
#axarr[1].set_ylim(0.8*ratio, 1.2*ratio)
axarr[1].axhline(y=ratio, ls="--", lw=0.5, color="k")
axarr[1].set_ylabel(r"$\frac{P_{0}}{P_{4}} \,$",rotation=0)

pl.xlim(0.05, 0.4);

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")

pl.savefig("./ZA_RSD.pdf", bbox_inches='tight')


