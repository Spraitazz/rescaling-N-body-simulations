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




pl.semilogx(data[:,0], data[:,1], label=r"original")
pl.semilogx(data[:,0], data[:,3], label=r"original, $z$ scaled")
pl.semilogx(data[:,0], data[:,4], ls="--", color="k", lw = 1.3, label=r"original, $z + L$ scaled")
pl.semilogx(data[:,0], data[:,2], color="r", lw = 0.7, label=r"target")


pl.ylim(0.1, 1.2)
pl.xlim(2.77, 21.1)

pl.legend(loc="best")
pl.xlabel(r"$R,\, [h^{-1} \rm{Mpc]}$")
pl.ylabel(r"$\sigma(R)$")

pl.savefig("./rescaled_test_sigmas_withz_cur.pdf", bbox_inches='tight')
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


