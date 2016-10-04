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



data_GR = np.loadtxt("./GR_Dplus.dat")
data_F4_1 = np.loadtxt("./F4_Dplus_k_1e+0.dat")
data_F4_2 = np.loadtxt("./F4_Dplus_k_1e-1.dat")
data_F4_3 = np.loadtxt("./F4_Dplus_k_1e-2.dat")
'''
data_F5_1 = np.loadtxt("./F5_Dplus_k_1e+0.dat")
data_F5_2 = np.loadtxt("./F5_Dplus_k_1e-1.dat")
data_F5_3 = np.loadtxt("./F5_Dplus_k_1e-2.dat")

data_F6_1 = np.loadtxt("./F6_Dplus_k_1e+0.dat")
data_F6_2 = np.loadtxt("./F6_Dplus_k_1e-1.dat")
data_F6_3 = np.loadtxt("./F6_Dplus_k_1e-2.dat")
'''


#pl.title(r"$b_{eff}$ vs $\nu_{min}$")
lw = 0.7



pl.plot(data_F4_1[:,1], data_F4_1[:,2], label=r"F4, $k = 1\, h\rm{Mpc^{-1}}$", color="k", linestyle="solid", linewidth=lw)
pl.plot(data_F4_2[:,1], data_F4_2[:,2], label=r"F4, $k = 0.1\, h\rm{Mpc^{-1}}$", color="k", linestyle="dashed", linewidth=lw)
pl.plot(data_F4_3[:,1], data_F4_3[:,2], label=r"F4, $k = 0.01\, h\rm{Mpc^{-1}}$", color="k", linestyle="-.", linewidth=lw)

'''
pl.plot(data_F5_1[:,1], data_F5_1[:,2], label=r"F5, $k = 1\, h\rm{Mpc^{-1}}$", color="b", linestyle="solid", linewidth=lw)
pl.plot(data_F5_2[:,1], data_F5_2[:,2], label=r"F5, $k = 0.1\, h\rm{Mpc^{-1}}$", color="b", linestyle="dashed", linewidth=lw)
pl.plot(data_F5_3[:,1], data_F5_3[:,2], label=r"F5, $k = 0.01\, h\rm{Mpc^{-1}}$", color="b", linestyle="-.", linewidth=lw)
'''
pl.plot(data_GR[:,1], data_GR[:,2], label="GR", color="r", linewidth=lw)
'''
pl.plot(data_F6_1[:,1], data_F6_1[:,2], label=r"F6, $k = 1 h\rm{Mpc^{-1}}$", color="b", linestyle="solid", linewidth=lw)
pl.plot(data_F6_2[:,1], data_F6_2[:,2], label=r"F6, $k = 0.1 h\rm{Mpc^{-1}}$", color="b", linestyle="dashed", linewidth=lw)
pl.plot(data_F6_3[:,1], data_F6_3[:,2], label=r"F6, $k = 0.01 h\rm{Mpc^{-1}}$", color="b", linestyle="-.", linewidth=lw)
'''
pl.ylim(0.3, 0.65)
pl.xlim(0.6, 1.2)
pl.legend(loc="best")

pl.xlabel(r"$z$")
pl.ylabel(r"$D_{+}(k,\, z)$")

pl.savefig("./Dplus_comparison.pdf", bbox_inches='tight')
pl.clf()


