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

no = sys.argv[1]
data = np.loadtxt("./cov_" + str(no) + "_sim0.dat")
data_p0p2 = np.loadtxt("./cov_" + str(no) + "_sim1.dat")
data_p2p2 = np.loadtxt("./cov_" + str(no) + "_sim2.dat")





if data.shape[0] != data.shape[1]:
	print "wtf matrix?"
	sys.exit(0)
	
dim = data.shape[0]
kmin_hod = 1.756204e-02
kmax_hod = 8.446513e-01
k30_hod = 5.329911e-01
kdiff_hod = k30_hod - kmin_hod
kstep_hod = kdiff_hod / 30

xticks = [6*i for i in xrange(0, 6)]
xticks.append(33)
yticks = [6*i for i in xrange(0, 6)]
yticks.append(33)
xlabels = [round(kmin_hod + i*kstep_hod,2) for i in xrange(0, 6)]
xlabels.append(round(kmax_hod, 2))
ylabels = [round(kmin_hod + i*kstep_hod,2) for i in xrange(0, 6)]
ylabels.append(round(kmax_hod, 2))

fig = pl.figure()
grid = AxesGrid(fig, 111,
                nrows_ncols=(2, 2),
                axes_pad=0.05,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )
                
#fig.imshow(data, origin='lower', interpolation='nearest', vmin = -0.3, vmax = 1.)
                
grid[0].imshow(data, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.3, vmax = 1.)
grid[1].imshow(data_p0p2, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.3, vmax = 1.)
grid[2].imshow(data_p0p2, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.3, vmax = 1.)
im = grid[3].imshow(data_p2p2, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.3, vmax = 1.)
grid.cbar_axes[0].colorbar(im)

grid.axes_llc.set_xticks(xticks, xlabels)
#grid.axes_llc.set_xlabels(xlabels)

#grid.cbar_axes[0].xticks(xticks, xlabels)
#fig.yticks(yticks, ylabels)	


pl.xlabel(r'$k \, [h$Mpc$^{-1}$]')
pl.ylabel(r'$k \, [h$Mpc$^{-1}$]')


pl.savefig("./HOD_" + str(no) + "ALL.pdf", bbox_inches='tight')

		
