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



data = np.loadtxt("./Pk_current_0.dat")
ks = data[0:33,0]
mono = data[0:33,1]
quadro = abs(data[0:33,2])


for i in xrange(1, 30):

	this_data = np.loadtxt("./Pk_current_" + str(i) + ".dat")
	mono = np.vstack((this_data[0:33,1], mono))
	quadro = np.vstack((abs(this_data[0:33,2]), quadro))
	
	
	
model = np.loadtxt("Pk_current_model.dat")

mean_mono = np.mean(mono, axis = 0)
mean_quadro = np.mean(quadro, axis = 0)



#-----

var_mono = np.var(mono, axis = 0)
var_quadro = np.var(quadro, axis = 0)

err_on_mean_mono = np.sqrt(var_mono) / np.sqrt(20)
err_on_mean_quadro = np.sqrt(var_quadro) / np.sqrt(20)


#END data input, error calc----------------------------------------



pl.loglog(model[:,0], model[:,1], color="Orange", label=r"model redshift-space $P_0$")
pl.loglog(model[:,0], model[:,2], color="BurlyWood", label=r"model redshift-space $P_2$")


#measured

pl.errorbar(ks, mean_mono, yerr=err_on_mean_mono, linestyle="None", color="OrangeRed", label=r"measured redshift-space $P_0$")

pl.errorbar(ks, mean_quadro, yerr=err_on_mean_quadro, linestyle="None", color="Chocolate", label=r"measured $P_{2}$")

#axarr[0].errorbar(ks, mean_hex, yerr=err_on_mean_hex, label=r"measured $P_{4}(k)$")

#models_legend = axarr[0].legend(handles = [line1,line2,line3,line4,line5], loc = 1, ncol = 2)
#measured_legend = axarr[0].legend(handles = [line6,line7,line8,line9,line10], loc = 4, ncol = 2)

pl.legend(loc="best", ncol=2)

pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")



pl.xlim(0.05, 0.5)
pl.ylim(3e1, 3e3)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")

pl.savefig("./ZA_RSD.pdf", bbox_inches='tight')


