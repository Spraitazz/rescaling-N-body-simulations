import numpy as np
import matplotlib.pyplot as pl
import matplotlib.pylab        as plt
import matplotlib.pyplot
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import ScalarFormatter
from   matplotlib.ticker       import FixedFormatter
import pylab                   as pl
import numpy                   as np
import math, os, sys
import glob, pickle
from matplotlib.lines import Line2D

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
              'text.fontsize':12,
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

data = np.loadtxt("./sigmas_TEST_NEW.dat")	


f, axarr = pl.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
axarr[0].set_xscale("log")
#axarr[0].set_yscale("log")

axarr[0].plot(data[:,0], data[:,1], label=r"original")
axarr[0].plot(data[:,0], data[:,3], label=r"original, $z$ scaled")
axarr[0].plot(data[:,0], data[:,4], ls="--", color="k", lw = 1.3, label=r"original, $z + L$ scaled")
axarr[0].plot(data[:,0], data[:,2], color="r", lw = 0.7, label=r"target")


axarr[0].legend(loc="best")
axarr[0].set_ylim(0.8, 1.9)

axarr[0].set_ylabel(r"$\sigma(R)$")

axarr[1].set_xscale("log")
axarr[1].scatter(data[:,0], np.divide(data[:,2], data[:,4]) - 1.0, label="fractional difference")
axarr[1].set_ylim(-0.05, 0.05)
axarr[1].axhline(y=0.0, ls="--", lw=0.5, color="k")
axarr[1].set_ylabel("residuals")

pl.xlim(0.87, 1.48)
pl.xlabel(r"$R,\, [h^{-1} \rm{Mpc]}$")

pl.savefig("./rescaled_test_sigmas_withz_cur.pdf", bbox_inches='tight')
pl.clf()



