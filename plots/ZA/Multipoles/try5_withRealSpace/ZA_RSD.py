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


startind = 0

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    


simno = 0
while os.path.isfile("./ZA_redshift_pk_" + str(startind + simno) + ".dat"):
	simno+=1

rows = file_len("./ZA_redshift_pk_" + str(startind)+ ".dat")

print "simulation number: " + str(simno) + ", file rows: " + str(rows)

data = np.loadtxt("./ZA_redshift_pk_" + str(startind) + ".dat")
ks = data[:,0]
data_mono = data[:,1]
data_quadro = abs(data[:,2])
data_hex = abs(data[:,3])

data_f4 = np.loadtxt("./ZA_redshift_f4_pk_" + str(startind)+ ".dat")
mono_f4 = data_f4[:,1]
quadro_f4 = abs(data_f4[:,2])

data_real = np.loadtxt("./ZA_real_pk_" + str(startind) + ".dat")
mono_real = data_real[:,1]

for i in xrange(1, simno):
	this_data = np.loadtxt("./ZA_redshift_pk_" + str(i + startind) + ".dat")
	data_mono = np.vstack((this_data[:,1], data_mono))
	data_quadro = np.vstack((abs(this_data[:,2]), data_quadro))
	data_hex = np.vstack((abs(this_data[:,3]), data_hex))
	
	this_data_real = np.loadtxt("./ZA_real_pk_" + str(i + startind) + ".dat")
	mono_real = np.vstack((this_data_real[:,1], mono_real))
	
	this_data_f4 = np.loadtxt("./ZA_redshift_f4_pk_" + str(i + startind) + ".dat")
	mono_f4 = np.vstack((this_data_f4[:,1], mono_f4))
	quadro_f4 = np.vstack((this_data_f4[:,2], quadro_f4))
	
	
model_real = np.loadtxt("./model_real_pk.dat")
model_redshift = np.loadtxt("./model_redshift_pk.dat")
model_f4 = np.loadtxt("./model_redshift_f4_pk.dat")

mean_mono = np.mean(data_mono, axis = 0)
mean_quadro = np.mean(data_quadro, axis = 0)
mean_hex = np.mean(data_hex, axis = 0)

mean_mono_real = np.mean(mono_real, axis = 0)

mean_mono_f4 = np.mean(mono_f4, axis = 0)
mean_quadro_f4 = np.mean(quadro_f4, axis = 0)

#-----

var_mono = np.var(data_mono, axis = 0)
var_quadro = np.var(data_quadro, axis = 0)
var_hex = np.var(data_hex, axis = 0)

var_mono_real = np.var(mono_real, axis = 0)

var_mono_f4 = np.var(mono_f4, axis = 0)
var_quadro_f4 = np.var(quadro_f4, axis = 0)

err_on_mean_mono = np.sqrt(var_mono) / np.sqrt(simno)
err_on_mean_quadro = np.sqrt(var_quadro) / np.sqrt(simno)
err_on_mean_hex = np.sqrt(var_hex) / np.sqrt(simno)

err_on_mean_mono_real = np.sqrt(var_mono_real) / np.sqrt(simno)

err_on_mean_mono_f4 = np.sqrt(var_mono_f4) / np.sqrt(simno)
err_on_mean_quadro_f4 = np.sqrt(var_quadro_f4) / np.sqrt(simno)

#END data input, error calc----------------------------------------

mono_over_quad = np.divide(mean_mono, mean_quadro)
mono_over_quad_f4 = np.divide(mean_mono_f4, mean_quadro_f4)

frac_err_mono = np.divide(err_on_mean_mono, mean_mono)
frac_err_quad = np.divide(err_on_mean_quadro, mean_quadro)

frac_err_mono_f4 = np.divide(err_on_mean_mono_f4, mean_mono_f4)
frac_err_quad_f4 = np.divide(err_on_mean_quadro_f4, mean_quadro_f4)

err_on_mean_mono_over_quad = np.multiply(np.add(frac_err_mono, frac_err_quad), mono_over_quad)
err_on_mean_mono_over_quad_f4 = np.multiply(np.add(frac_err_mono, frac_err_quad), mono_over_quad)

f_GR = pow(0.24, 0.545)
#ratio = f_GR*f_GR*8.0/35.0
ratio = (1.0 + (2.0/3.0)*f_GR + (1.0/5.0)*f_GR*f_GR) / ((4.0/3.0)*f_GR + (4.0/7.0)*f_GR*f_GR)

f_f4 = pow(0.24, 0.4)
ratio_f4 = (1.0 + (2.0/3.0)*f_f4 + (1.0/5.0)*f_f4*f_f4) / ((4.0/3.0)*f_f4 + (4.0/7.0)*f_f4*f_f4)

f, axarr = pl.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
axarr[0].set_xscale("log")
axarr[0].set_yscale("log")

#models

line1, = axarr[0].plot(model_real[:,0], model_real[:,1], color="Orange", label=r"model real-space $P(k)$")

line2, = axarr[0].plot(model_redshift[:,0], model_redshift[:,1], color="BurlyWood", label=r"model $P_{0}(k) (\gamma = 0.545)$")
line3, = axarr[0].plot(model_redshift[:,0], model_redshift[:,2], color="DarkCyan", label=r"model $P_{2}(k) (\gamma = 0.545)$")

line4, = axarr[0].plot(model_f4[:,0], model_f4[:,1], color="PaleGreen", label=r"model $P_{0}(k) (\gamma = 0.4)$")
line5, = axarr[0].plot(model_f4[:,0], model_f4[:,2], color="PeachPuff", label=r"model $P_{2}(k) (\gamma = 0.4)$")

#measured

axarr[0].errorbar(ks, mean_mono_real, yerr=err_on_mean_mono_real, linestyle="None", color="OrangeRed", label=r"measured real-space $P(k)$")

axarr[0].errorbar(ks, mean_mono, yerr=err_on_mean_mono, linestyle="None", color="Chocolate", label=r"measured $P_{0}(k) (\gamma = 0.545)$")
axarr[0].errorbar(ks, mean_quadro, yerr=err_on_mean_quadro, linestyle="None", color="DarkBlue", label=r"measured $P_{2}(k) (\gamma = 0.545)$")

axarr[0].errorbar(ks, mean_mono_f4, yerr=err_on_mean_mono_f4, linestyle="None", color="OliveDrab", label=r"measured $P_{0}(k) (\gamma = 0.4)$")
axarr[0].errorbar(ks, mean_quadro_f4, yerr=err_on_mean_quadro_f4, linestyle="None", color="Salmon", label=r"measured $P_{2}(k) (\gamma = 0.4)$")
#axarr[0].errorbar(ks, mean_hex, yerr=err_on_mean_hex, label=r"measured $P_{4}(k)$")

#models_legend = axarr[0].legend(handles = [line1,line2,line3,line4,line5], loc = 1, ncol = 2)
#measured_legend = axarr[0].legend(handles = [line6,line7,line8,line9,line10], loc = 4, ncol = 2)

axarr[0].legend(loc="best", ncol=2)
axarr[0].set_ylim(6e1,2e3)
#axarr[0].set_xlim(kmin, kmax)
axarr[0].set_ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

axarr[1].set_xscale("log")
axarr[1].errorbar(ks, mono_over_quad, yerr=err_on_mean_mono_over_quad, label=r"$\gamma = 0.54$", fmt="o", ms = 0.5)
axarr[1].errorbar(ks, mono_over_quad_f4, yerr=err_on_mean_mono_over_quad_f4, label=r"$\gamma = 0.4$", fmt="o", ms = 0.5)
axarr[1].set_ylim(1.0, 3.0)
axarr[1].axhline(y=ratio, ls="--", lw=0.5, color="k")
axarr[1].axhline(y=ratio_f4, ls="--", lw=0.5, color="k")
axarr[1].set_ylabel(r"$\frac{P_{0}(k)}{P_{2}(k)} \,\,\,\,\,$", rotation=0)
axarr[1].legend(loc="best")

pl.xlim(0.05, 0.2);

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")

pl.savefig("./ZA_RSD.pdf", bbox_inches='tight')


