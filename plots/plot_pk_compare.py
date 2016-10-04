import numpy as np
import matplotlib.pyplot as pl
import os
import sys

#count = 0
#while os.path.isfile("./power_spectra/power_spectrum" + str(count) + ".pdf"):
	#count+=1
#if os.path.isfile("./power_spectrum.pdf"):


input_len = len(sys.argv)
fileNames = sys.argv

for i in range(1, input_len-1):
	data = np.loadtxt(fileNames[i])
	pl.loglog(data[:,0], data[:,1], label=fileNames[i])
	
pl.legend(loc="best")

#data = np.loadtxt("GR_1fold.txt")
#pl.loglog(data[:,0], data[:,1]/8.) 
volume = 1500.0 * 1500.0 * 1500.0
number = 50000.0
inv_no_dens = volume / number
#pl.axhline(y=inv_no_dens, linestyle="--")

#pl.ylim(30000.0,90000.0)
#pl.xlim(0.02, 0.9)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./" + fileNames[input_len - 1] +  ".pdf")
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


