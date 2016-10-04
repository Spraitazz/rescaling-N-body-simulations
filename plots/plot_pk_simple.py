import numpy as np
import matplotlib.pyplot as pl
import os
import sys

'''
count = 0
while os.path.isfile("./power_spectra/power_spectrum" + str(count) + ".pdf"):
	count+=1
#if os.path.isfile("./power_spectrum.pdf"):
'''	

data = np.loadtxt(sys.argv[1])
'''
for i in range(0, len(data)):
	if data[i,1] == 0.:
		del data[i,:]
		print("zeros found")
'''
monopole = pl.loglog(data[:,0], data[:,1], label="monopole")
quadrupole = pl.loglog(data[:,0], data[:,2], label="quadrupole")
pl.legend()

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.ylim(200,120000)
#pl.xlim(0.005, 5)

pl.savefig("./" + sys.argv[2] + ".png")
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


