import numpy as np
import matplotlib.pyplot as pl
import os
import sys



data_GR = np.loadtxt("./GR_monopole.dat")
data_F5 = np.loadtxt("./F5_monopole.dat")
data_GR_resc = np.loadtxt("./GR_rescaled_test.dat")

	


#data = np.loadtxt("GR_1fold.txt")
pl.title("GR -> F5 attempt")
pl.loglog(data_GR[:,0], data_GR[:,1], label="GR monopole")
pl.loglog(data_F5[:,0], data_F5[:,1], label="F5 monopole")
pl.loglog(data_GR_resc[:,0], data_GR_resc[:,1], label="GR rescaled monopole")
volume = 1500.0 * 1500.0 * 1500.0
number = 50000.0
inv_no_dens = volume / number
#pl.axhline(y=inv_no_dens, linestyle="--")

pl.ylim(3.9e3,2.5e5)
pl.xlim(0.018, 1.2)

pl.legend(loc="best")
pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./rescaled_test.pdf")
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


