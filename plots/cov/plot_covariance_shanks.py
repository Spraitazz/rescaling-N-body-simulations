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

#-------------------------------------------------------

no = sys.argv[1]
data = np.loadtxt("./cov_" + str(no) + "_sim0.dat")

matrices = []
expected = np.loadtxt("./cov_305_sim0.dat")


for i in xrange(3, 30):
	tmp = np.loadtxt("./cov_" + str(i) + "_sim0.dat")
	matrices.append(tmp)
	
def partialSum(n):
	return matrices[n]	
	
def shanks(n):
	if n < 1 or n > 25:
		return 0
	else:
		psn = partialSum(n)
		psn_minus = partialSum(n-1)
		psn_plus = partialSum(n+1)
		top = np.multiply(psn_plus,psn_minus) - np.multiply(psn, psn)
		bottom = psn_plus - 2*psn + psn_minus
		np.fill_diagonal(top, 1)		
		np.fill_diagonal(bottom, 1)
		inter = top/bottom
		return inter
	
def discrep(matr):
	inter = (matr-expected)/expected
	return abs(np.sum(inter))
	
data = shanks(int(no))
#data = np.flipud(data)

if data.shape[0] != data.shape[1]:
	print "wtf matrix?"
	sys.exit(0)
	
dim = data.shape[0]
kmin_hod = 1.756204e-02
kmax_hod = 8.446513e-01
k30_hod = 5.329911e-01
kdiff_hod = k30_hod - kmin_hod
kstep_hod = kdiff_hod / 30



fig, ax = pl.subplots()
pl.imshow(data, origin='lower', cmap='YlOrRd', interpolation='nearest', vmin = -0.3, vmax = 1.)
#ax.xaxis.tick_top()

xticks = [6*i for i in xrange(0, 6)]
xticks.append(33)
yticks = [6*i for i in xrange(0, 6)]
yticks.append(33)
xlabels = [round(kmin_hod + i*kstep_hod,2) for i in xrange(0, 6)]
xlabels.append(round(kmax_hod, 2))
ylabels = [round(kmin_hod + i*kstep_hod,2) for i in xrange(0, 6)]
ylabels.append(round(kmax_hod, 2))



#plt.gca().invert_yaxis()
#plt.gca().tick_bottom()
#ylabels = ylabels[::-1]
#yticks = yticks[::-1]

#pl.xticks([0, 6, 13, 20, 24], [str(0.019), str(0.052),  str(0.13), str(0.34), str(0.88)])
#pl.yticks([0, 6, 13, 20, 24], [str(0.019), str(0.052),  str(0.13), str(0.34), str(0.88)])

pl.xticks(xticks, xlabels)
pl.yticks(yticks, ylabels)

pl.xlabel(r'$k \, [h$Mpc$^{-1}$]')
pl.ylabel(r'$k \, [h$Mpc$^{-1}$]')

pl.colorbar()

#pl.show()
pl.savefig("./HOD_" + str(no) + "_shanks.pdf", bbox_inches='tight')
		
