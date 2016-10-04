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

	
data_current = np.loadtxt("./delsq_current.dat")
data_target = np.loadtxt("./delsq_target.dat")
data_current_z = np.loadtxt("./delsq_current_z.dat")
data_current_zs = np.loadtxt("./delsq_current_zs.dat")



f, axarr = pl.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
axarr[0].set_xscale("log")
axarr[0].set_yscale("log")


axarr[0].plot(data_current[:,0], data_current[:,1], label="current")

axarr[0].plot(data_target[:,0], data_target[:,1], label="target")

axarr[0].plot(data_current_z[:,0], data_current_z[:,1], label="current z")

axarr[0].plot(data_current_zs[:,0], data_current_zs[:,1], label="current zs")

residuals = np.divide(data_current_zs[:,1], data_target[:,1]) - 1.0

axarr[0].legend(loc="best")
#axarr[0].set_ylim(0.0,2.3e3)
#axarr[0].set_xlim(kmin, kmax)
axarr[0].set_ylabel(r"$\Delta^{2}(k)$")

axarr[1].set_xscale("log")
axarr[1].plot(data_target[:,0], residuals, label="absolute fractional difference")

axarr[1].set_ylim(-0.2, 0.2)
axarr[0].set_xlim(0.01, 0.1)
axarr[1].axhline(y=0.0, ls="-", lw=0.5, color="k")

axarr[1].set_ylabel("Residual")

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


pl.savefig("./delsq_plot.pdf", bbox_inches='tight')


