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



data = np.loadtxt("./cdm.dat")
data = np.log10(data)

pl.plot(data[:,0], data[:,2], label=r"$\Omega_m = 0.3, z=0.6$")
pl.plot(data[:,0], data[:,4], label=r"$\Omega_m = 0.3, z=1.2$")
pl.plot(data[:,0], data[:,1], label=r"$\Omega_m = 0.1, z=0.6$")
pl.plot(data[:,0], data[:,3], label=r"$\Omega_m = 0.1, z=1.2$")



pl.legend(loc="best")

pl.ylabel(r"$c_{\rm{dm}}$")

#pl.xlim(5e12, 1e15)
#pl.xlim(0.026, 0.3);

pl.xlabel(r"$M, [M_{\odot}]$")

pl.savefig("./ZA_RSD.pdf", bbox_inches='tight')


