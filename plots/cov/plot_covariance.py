import matplotlib.pylab        as plt
import matplotlib.pyplot as pl
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import ScalarFormatter
from   matplotlib.ticker       import FixedFormatter
import numpy                   as np
import math, os, sys
import glob, pickle
from mpl_toolkits.axes_grid1 import AxesGrid

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

#-------------------------------------------------------

f, ax = pl.subplots()

data_p0p0 = np.loadtxt("./305_mocks_cov_0.dat")
data_p0p2 = np.loadtxt("./305_mocks_cov_1.dat")
data_p2p2 = np.loadtxt("./305_mocks_cov_2.dat")

data_1 = np.hstack((data_p0p2, data_p2p2))
data_2 = np.hstack((data_p0p0, data_p0p2))

data = np.vstack((data_2, data_1))

data = np.loadtxt("./305_mocks_cov.dat")

pl.imshow(data, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.5, vmax = 1.)

pl.xticks([0, 10, 20, 30], [str(0.018), str(0.59),  str(0.19), str(0.6)])
pl.yticks([0, 10, 20, 30], [str(0.018), str(0.59),  str(0.19), str(0.6)])

pl.xlabel(r'$k \, [h$Mpc$^{-1}$]')
pl.ylabel(r'$k \, [h$Mpc$^{-1}$]')

# vertical 
f = 1./6. 

# horizontal
g = 1.1/7.

pl.text(1.2*g,  f/2., r'$P_0 \times \ P_0$', transform=ax.transAxes)
pl.text(4.5*g, 5.5*f, r'$P_2 \times \ P_2$',  transform=ax.transAxes)

pl.text(1.2*g,  5.5*f, r'$P_0 \times \ P_2$', transform=ax.transAxes)
pl.text(4.5*g, f/2., r'$P_2 \times \ P_0$',  transform=ax.transAxes)

pl.colorbar()


'''
if data_p0p0.shape[0] != data_p0p0.shape[1]:
	print "wtf matrix?"
	sys.exit(0)
	
dim = data_p0p0.shape[0]
kmin_hod = 1.756204e-02
kmax_hod = 8.446513e-01
k30_hod = 5.329911e-01
kdiff_hod = k30_hod - kmin_hod
kstep_hod = kdiff_hod / 30

xticks = [6*i for i in xrange(0, 6)]
#xticks.append(33)
yticks = [6*i for i in xrange(0, 6)]
#yticks.append(33)
xlabels = [round(kmin_hod + i*kstep_hod,2) for i in xrange(0, 6)]
xlabels.append(round(kmax_hod, 2))
ylabels = [round(kmin_hod + i*kstep_hod,2) for i in xrange(0, 6)]
ylabels.append(round(kmax_hod, 2))
'''
'''
fig = pl.figure()
grid = AxesGrid(fig, 111,
                nrows_ncols=(2, 2),
                axes_pad=0.05,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )
                
pl.imshow(data, origin='lower', interpolation='nearest', vmin = -0.3, vmax = 1.)
                
grid[0].imshow(data, origin='lower', interpolation='nearest', vmin = -0.3, vmax = 1.)
grid[1].imshow(data_p0p2)
grid[2].imshow(data_p0p2)
im = grid[3].imshow(data_p2p2)
grid.cbar_axes[0].colorbar(im)

grid.axes_llc.set_xticks(xticks, xlabels)
#grid.axes_llc.set_xlabels(xlabels)

#grid.cbar_axes[0].xticks(xticks, xlabels)
#fig.yticks(yticks, ylabels)
	

'''
'''
f, axarr = pl.subplots(2, 2)
axarr[1,0].imshow(data_p0p0, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.5, vmax = 1.)
axarr[1,1].imshow(data_p0p2, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.5, vmax = 1.)
axarr[0,0].imshow(data_p0p2, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.5, vmax = 1.)
im = axarr[0,1].imshow(data_p2p2, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.5, vmax = 1.)

#f.subplots_adjust(hspace=0.0, wspace=1.0)
pl.tight_layout(pad=0, w_pad=0.3, h_pad=0.3)

#axarr[0,0].axis("off")
#axarr[1].axis("off")
#axarr[1,1].axis("off")

axarr[1,0].title.set_text(r"$P_0 \times P_0$")
#axarr[1,1].title.set_text(r"$P_0 \times P_2$")
#axarr[0,0].title.set_text(r"$P_0 \times P_2$")
axarr[0,1].title.set_text(r"$P_2 \times P_2$")

axarr[1,0].set_xticks(xticks)
axarr[1,0].set_xticklabels(xlabels)

axarr[0,1].set_yticks(yticks)
axarr[0,1].set_yticklabels(ylabels)

#pl.colorbar(im, use_gridspec=True)
f.colorbar(im, ax=axarr.ravel().tolist())
#pl.imshow(data_p2p2, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.3, vmax = 1.)


#pl.xticks(xticks, xlabels)
#pl.yticks(yticks, ylabels)

pl.xlabel(r'$k \, [h$Mpc$^{-1}$]')
pl.ylabel(r'$k \, [h$Mpc$^{-1}$]')

#pl.subplots_adjust(wspace=0, hspace=0)

#pl.colorbar()
'''

pl.savefig("./HOD_305_ALL_SUBPLOTS.pdf", bbox_inches='tight')

		
