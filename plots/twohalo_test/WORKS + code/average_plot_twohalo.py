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
while os.path.isfile("./measured_pk_" + str(startind + simno) + ".dat"):
	simno+=1

rows = file_len("./measured_pk_" + str(startind)+ ".dat")
print "simulation number: " + str(simno) + ", rows: " + str(rows)



data = np.loadtxt("./measured_pk_" + str(startind) + ".dat")
ks = data[:,0]
data_mono = data[:,1]
data_quadro = data[:,2]

for i in xrange(1, simno):
	this_data = np.loadtxt("./measured_pk_" + str(i + startind) + ".dat")
	data_mono = np.vstack((this_data[:,1], data_mono))
	data_quadro = np.vstack((this_data[:,2], data_quadro))	
	
	
data_model = np.loadtxt("./model_pk.dat")

kmin = ks[0]
kmax = ks[rows-1]

mean_mono = np.mean(data_mono, axis = 0)
mean_quadro = np.mean(data_quadro, axis = 0)

var_mono = np.var(data_mono, axis = 0)
var_quadro = np.var(data_quadro, axis = 0)

err_on_mean_mono = np.sqrt(var_mono) / np.sqrt(simno)
err_on_mean_quadro = np.sqrt(var_quadro) / np.sqrt(simno)

#END data input, error calc----------------------------------------

f_model = interpolate.interp1d(data_model[:,0], data_model[:,1])
spectrum_model_interpolated = []
errorbars_residuals = []
for i in xrange(0, len(ks)):
	if (ks[i] >= kmin and ks[i] <= kmax):	
		model_current = float(f_model(ks[i]))	
		spectrum_model_interpolated.append(model_current)
		errorbars_residuals.append(err_on_mean_mono[i] / model_current)
				

frac_diffs = abs(np.divide(mean_mono, spectrum_model_interpolated)) - 1.0

f, axarr = pl.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
axarr[0].set_xscale("log")
axarr[0].set_yscale("log")

axarr[0].errorbar(ks, mean_mono, yerr=err_on_mean_mono, label="monopole")
axarr[0].plot(data_model[:,0], data_model[:,1], label="model monopole")

axarr[0].errorbar(ks, mean_quadro, yerr=err_on_mean_quadro, label="quadrupole")
axarr[0].plot(data_model[:,0], data_model[:,2], label="model quadrupole")

axarr[0].legend(loc="best")
axarr[0].set_ylim(0.0,2.3e3)
axarr[0].set_xlim(kmin, kmax)
axarr[0].set_ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

axarr[1].set_xscale("log")
axarr[1].errorbar(ks, frac_diffs, yerr=errorbars_residuals, label="absolute fractional difference", fmt="o", ms = 0.5)
axarr[1].set_ylim(-0.1, 0.1)
axarr[1].axhline(y=0.0, ls="--", lw=0.5, color="k")
axarr[1].set_ylabel("residuals")

'''
ax = pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")

pl.errorbar(ks, mean, yerr=err_on_mean, label="monopole")
pl.plot(data_model[:,0], data_model[:,1], label="model monopole")
#pl.plot(data_model_nospline[:,0], data_model_nospline[:,1] + inv_dens, label="model without cubic spline")

pl.legend(loc="best")


pl.ylim(5e2,2e4)
pl.xlim(0.015, 0.65)
'''

pl.xlabel(r"$k, \ [h \, \rm{Mpc} ^{-1}]$")


pl.savefig("./ZA_test.pdf", bbox_inches='tight')


