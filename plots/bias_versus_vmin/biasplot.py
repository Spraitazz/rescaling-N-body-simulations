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



data = np.loadtxt("./bias_plot.dat")

vmin = 0.207563
vmax = 2.031407	


#data = np.loadtxt("GR_1fold.txt")
#pl.title(r"$b_{eff}$ vs $\nu_{min}$")
pl.plot(data[:,0], data[:,1])
#pl.loglog(data_F5[:,0], data_F5[:,1], label="F5 monopole")
pl.axvline(x=vmin, linestyle="--", color="black")
pl.axvline(x=vmax, linestyle="--", color="black")
#pl.ylim(3.9e3,2.5e5)
pl.xlim(0.0, 6.0)


pl.xlabel(r"$\nu_{\rm{min}}$")
pl.ylabel(r"$b_{\rm{eff}}$")

pl.savefig("./bias_test_GR1.pdf", bbox_inches='tight')
pl.clf()


