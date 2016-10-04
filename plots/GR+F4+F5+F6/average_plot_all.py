import matplotlib.pylab        as plt
import matplotlib.pyplot as pl
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import ScalarFormatter
from   matplotlib.ticker       import FixedFormatter
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

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

simno = int(sys.argv[1])
rows = file_len("./GR_combined_1.dat")
#print rows
data_GR = []
data_F4 = []
data_F5 = []
data_F6 = []
monos_GR = np.zeros(rows)
quadros = np.zeros(rows)
monos_F4 = np.zeros(rows)
monos_F5 = np.zeros(rows)
monos_F6 = np.zeros(rows)


for i in xrange(0, simno):
	this_data_GR = np.loadtxt("./GR_combined_" + str(i+1) + ".dat")
	this_data_F4 = np.loadtxt("./F4_combined_" + str(i+1) + ".dat")
	this_data_F5 = np.loadtxt("./F5_combined_" + str(i+1) + ".dat")
	this_data_F6 = np.loadtxt("./F6_combined_" + str(i+1) + ".dat")

	monos_GR += this_data_GR[:,1]
	monos_F4 += this_data_F4[:,1]
	monos_F5 += this_data_F5[:,1]
	monos_F6 += this_data_F6[:,1]

	data_GR.append(this_data_GR)
	data_F4.append(this_data_F4)
	data_F5.append(this_data_F5)
	data_F6.append(this_data_F6)

ks = data_GR[0][:,0]
monos_GR /= simno
monos_F4 /= simno
monos_F5 /= simno
monos_F6 /= simno

ax = pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")

pl.errorbar(ks, monos_GR, label="GR")
pl.errorbar(ks, monos_F4, label="F4")
pl.errorbar(ks, monos_F5, label="F5")
pl.errorbar(ks, monos_F6, label="F6")

pl.legend(loc="best")


pl.ylim(5e3,2e5)
pl.xlim(0.003, 1.3)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./Different_gravities_realspace.pdf", bbox_inches="tight")
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


