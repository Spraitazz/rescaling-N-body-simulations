import matplotlib.pylab        as plt
import matplotlib.pyplot as pl
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

#-----------------------------------------------------------------------

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

data = np.loadtxt("./test_Pk_model_spline.dat")

	
ax = pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")


pl.loglog(data[:,0], data[:,1], label="z = 0")
pl.loglog(data[:,0], data[:,2], label="z = 2")
pl.loglog(data[:,0], data[:,3], label="z = 4")
pl.loglog(data[:,0], data[:,4], label="z = 6")
pl.loglog(data[:,0], data[:,5], label="z = 8")

pl.legend(loc="best")

#pl.ylim(7.9e3,8.2e3)
#pl.xlim(0.025, 0.8)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./Pk_diffz_comp_model.pdf", bbox_inches='tight')


