import numpy as np
import matplotlib.pyplot as pl
import os
import sys



data_GR = np.loadtxt("./linear_matter_pk_sig8_0.593_z_0.75.dat")
data_F5 = np.loadtxt("./linear_matter_pk_0.76.dat")
mult = 0.36**3.0

	


#data = np.loadtxt("GR_1fold.txt")
pl.title("Box scaling test model->model")
pl.loglog(data_GR[:,0], data_GR[:,1], label="lin pk z 0.75 current")
pl.loglog(data_F5[:,0], data_F5[:,1], label="lin pk z 0.76 target")
pl.loglog(data_GR[:,0], data_GR[:,1]*mult, label="after volume change")
volume = 1500.0 * 1500.0 * 1500.0
number = 50000.0
inv_no_dens = volume / number
#pl.axhline(y=inv_no_dens, linestyle="--")

#pl.ylim(3.9e3,2.5e5)
#pl.xlim(0.018, 1.2)

pl.legend(loc="best")
pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./rescaled_test.pdf")
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


