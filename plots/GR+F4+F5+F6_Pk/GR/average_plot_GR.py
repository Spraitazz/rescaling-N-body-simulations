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
rows0 = file_len("./GR_redshift_1_0fold.dat")
rows2 = file_len("./GR_redshift_1_2fold.dat")
rows4 = file_len("./GR_redshift_1_4fold.dat")
rows8 = file_len("./GR_redshift_1_8fold.dat")
ind0 = 11
ind11 = 6
ind12 = 13
ind21 = 8
ind22 = 13
ind31 = 8
ind32 = 17
ind41 = 12
ind42 = 13

length = ind0 + (ind12-ind11) + (ind22 - ind21) + (ind32-ind31)
#print rows
data_0 = []
data_2 = []
data_4 = []
data_8 = []
data_16 = []
monos_0 = np.zeros(ind0)
monos_2 = np.zeros(ind12-ind11)
monos_4 = np.zeros(ind22-ind21)
monos_8 = np.zeros(ind32-ind31)
monos_16 = np.zeros(ind42-ind41)

monos_0_rs = np.zeros(ind0)
monos_2_rs = np.zeros(ind12-ind11)
monos_4_rs = np.zeros(ind22-ind21)
monos_8_rs = np.zeros(ind32-ind31)
monos_16_rs = np.zeros(ind42-ind41)

quadros_0 = np.zeros(ind0)
quadros_2 = np.zeros(ind12-ind11)
quadros_4 = np.zeros(ind22-ind21)
quadros_8 = np.zeros(ind32-ind31)
quadros_16 = np.zeros(ind42-ind41)
#quadros_rs = np.zeros(rows)
#monos_rs = np.zeros(rows)

this_data_rs_0 = np.loadtxt("./GR_redshift_" + str(1) + "_0fold.dat")[0:ind0,:]
this_data_rs_2 = np.loadtxt("./GR_redshift_" + str(1) + "_2fold.dat")[ind11:ind12,:]
this_data_rs_4 = np.loadtxt("./GR_redshift_" + str(1) + "_4fold.dat")[ind21:ind22,:]
this_data_rs_8 = np.loadtxt("./GR_redshift_" + str(1) + "_8fold.dat")[ind31:ind32:,:]
this_data_rs_16 = np.loadtxt("./GR_redshift_" + str(1) + "_16fold.dat")[ind41:ind42:,:]
#this_data_rs_32 = np.loadtxt("./GR_redshift_" + str(1) + "_16fold.dat")[ind41:ind42:,:]

this_data_0 = np.loadtxt("./GR_" + str(1) + "_1fold.dat")[0:ind0,:]
this_data_2 = np.loadtxt("./GR_" + str(1) + "_2fold.dat")[ind11:ind12,:]
this_data_4 = np.loadtxt("./GR_" + str(1) + "_4fold.dat")[ind21:ind22,:]
this_data_8 = np.loadtxt("./GR_" + str(1) + "_8fold.dat")[ind31:ind32:,:]
this_data_16 = np.loadtxt("./GR_" + str(1) + "_16fold.dat")[ind41:ind42:,:]
#this_data_32 = np.loadtxt("./GR_" + str(1) + "_16fold.dat")[ind41:ind42:,:]

data_joint = np.vstack((this_data_0, this_data_2))
data_joint = np.vstack((data_joint, this_data_4))
data_joint = np.vstack((data_joint, this_data_8))
data_joint = np.vstack((data_joint, this_data_16))
#data_joint = np.vstack((data_joint, this_data_32))

data_joint_rs = np.vstack((this_data_rs_0, this_data_rs_2))
data_joint_rs = np.vstack((data_joint_rs, this_data_rs_4))
data_joint_rs = np.vstack((data_joint_rs, this_data_rs_8))
data_joint_rs = np.vstack((data_joint_rs, this_data_rs_16))

data_mono = data_joint[:,1]
data_mono_rs = data_joint_rs[:,1]
data_quadro_rs = data_joint_rs[:,2]


for i in xrange(1, simno):
	#this_data = np.loadtxt("./GR_redshift_" + str(i+1) + ".dat")
	this_data_rs_0 = np.loadtxt("./GR_redshift_" + str(i+1) + "_0fold.dat")[0:ind0,:]
	this_data_rs_2 = np.loadtxt("./GR_redshift_" + str(i+1) + "_2fold.dat")[ind11:ind12,:]
	this_data_rs_4 = np.loadtxt("./GR_redshift_" + str(i+1) + "_4fold.dat")[ind21:ind22,:]
	this_data_rs_8 = np.loadtxt("./GR_redshift_" + str(i+1) + "_8fold.dat")[ind31:ind32:,:]
	this_data_rs_16 = np.loadtxt("./GR_redshift_" + str(i+1) + "_16fold.dat")[ind41:ind42:,:]
	
	this_data_0 = np.loadtxt("./GR_" + str(i+1) + "_1fold.dat")[0:ind0,:]
	this_data_2 = np.loadtxt("./GR_" + str(i+1) + "_2fold.dat")[ind11:ind12,:]
	this_data_4 = np.loadtxt("./GR_" + str(i+1) + "_4fold.dat")[ind21:ind22,:]
	this_data_8 = np.loadtxt("./GR_" + str(i+1) + "_8fold.dat")[ind31:ind32:,:]
	this_data_16 = np.loadtxt("./GR_" + str(i+1) + "_16fold.dat")[ind41:ind42:,:]
	
	data_joint = np.vstack((this_data_0, this_data_2))
	data_joint = np.vstack((data_joint, this_data_4))
	data_joint = np.vstack((data_joint, this_data_8))
	data_joint = np.vstack((data_joint, this_data_16))
	
	data_joint_rs = np.vstack((this_data_rs_0, this_data_rs_2))
	data_joint_rs = np.vstack((data_joint_rs, this_data_rs_4))
	data_joint_rs = np.vstack((data_joint_rs, this_data_rs_8))
	data_joint_rs = np.vstack((data_joint_rs, this_data_rs_16))
	
	data_mono = np.vstack((data_mono, data_joint[:,1]))
	data_mono_rs = np.vstack((data_mono_rs, data_joint_rs[:,1]))
	data_quadro_rs = np.vstack((data_quadro_rs, data_joint_rs[:,2]))
	
	
	data_0.append(this_data_rs_0)
	data_2.append(this_data_rs_2)
	data_4.append(this_data_rs_4)
	data_8.append(this_data_rs_8)
	data_16.append(this_data_rs_16)
	


mean_mono = np.mean(data_mono, axis = 0) 
mean_mono_rs = np.mean(data_mono_rs, axis = 0)
mean_quadro_rs = np.mean(data_quadro_rs, axis = 0)


var_mono = np.sqrt(np.var(data_mono, axis = 0)) / np.sqrt(simno)
var_mono_rs = np.sqrt(np.var(data_mono_rs, axis = 0)) / np.sqrt(simno)
var_quadro_rs = np.sqrt(np.var(data_quadro_rs, axis = 0)) / np.sqrt(simno)


mean_quadro_rs = abs(mean_quadro_rs)

ks_0 = data_0[0][:,0]
ks_2 = data_2[0][:,0]
ks_4 = data_4[0][:,0]
ks_8 = data_8[0][:,0]
ks_16 = data_16[0][:,0]

ks = np.concatenate((ks_0, ks_2))
ks = np.concatenate((ks, ks_4))
ks = np.concatenate((ks, ks_8))
ks = np.concatenate((ks, ks_16))
#monos /= simno
#monos_rs /= simno
#quadros_rs /= simno

#err_monos = np.zeros(rows)
#err_monos_rs = np.zeros(rows)
#err_quadros_rs = np.zeros(rows)
#err_zeros = np.zeros(rows)

#variances
'''
for i in xrange(0, simno):
	err_monos += ((data[i][:,1] - monos)**2.0)/simno
	err_monos_rs += ((data_rs[i][:,1] - monos_rs)**2.0)/simno
	err_quadros_rs += ((data_rs[i][:,2] - quadros_rs)**2.0)/simno

#want stdevs
err_monos = err_monos**0.5
err_monos_rs = err_monos_rs**0.5
err_quadros_rs = err_quadros_rs**0.5

#error on mean
err_monos /= np.sqrt(simno)
err_monos_rs /= np.sqrt(simno)
err_quadros_rs /= np.sqrt(simno)
	'''
ax = pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")

pl.errorbar(ks, mean_mono, yerr=var_mono, label=r"real-space $P_0(k)$")
pl.errorbar(ks, mean_mono_rs, yerr=var_mono_rs, label=r"redshift-space $P_0(k)$")
pl.errorbar(ks[2:len(ks)-2], mean_quadro_rs[2:len(ks)-2], yerr=var_quadro_rs[2:len(ks)-2], label=r"redshift-space $P_2(k)$")
#pl.plot(ks_2, monos_2, label="2folded")
#pl.plot(ks_4, monos_4, label="4folded")
#pl.plot(ks_8, monos_8, label="8folded")

#pl.errorbar(ks, monos, yerr=err_monos, label=r"real-space $P_0(k)$")
#pl.errorbar(ks, monos_rs, yerr=err_monos_rs, label=r"redshift-space $P_0(k)$")
#pl.errorbar(ks[6:], quadros_rs[6:], yerr=err_quadros_rs[6:], label=r"redshift-space $P_2(k)$")

pl.legend(loc="best")


pl.ylim(4e1,2.2e5)
pl.xlim(0.01, 1.0)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./GR_ALL.pdf", bbox_inches="tight")
pl.clf()


