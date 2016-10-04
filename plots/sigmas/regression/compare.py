import numpy as np
import matplotlib.pyplot as pl
import os
import sys

#count = 0
#while os.path.isfile("./power_spectra/power_spectrum" + str(count) + ".pdf"):
	#count+=1
#if os.path.isfile("./power_spectrum.pdf"):


	


model = np.loadtxt("./linear_matter_pk_sig8_0.593_z_0.75.dat")
extended = np.loadtxt("./Pk_regression_test.dat")
pl.loglog(model[:,0], model[:,1], label="model pk") 
pl.loglog(extended[:,0], extended[:,1], label="model extended") 


pl.legend(loc="best")
#pl.ylim(30000.0,90000.0)
#pl.xlim(0.02, 0.9)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./result.pdf")

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


