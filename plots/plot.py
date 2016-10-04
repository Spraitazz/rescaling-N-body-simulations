import numpy as np
import matplotlib.pyplot as pl
import os
import sys

#count = 0
#while os.path.isfile("./power_spectra/power_spectrum" + str(count) + ".pdf"):
	#count+=1
#if os.path.isfile("./power_spectrum.pdf"):


#randoms = np.loadtxt("randoms_pk_joined.dat")
#two_halo = np.loadtxt("randoms_2halo_nonoise.dat")
#final = np.loadtxt("randoms_combined_pk_joined.dat")
model = np.loadtxt("./model_pk_onehalo.dat")
one_halo_ngp = np.loadtxt("./randoms_onehalo_NGP.dat")
one_halo_cic = np.loadtxt("./randoms_onehalo_CIC.dat")
one_halo = np.loadtxt("./randoms_onehalo_0.dat")
#pl.loglog(two_halo[:,0], two_halo[:,1], label="2-halo term minus noise")
pl.loglog(model[:,0], model[:,1], label="model monopole")
#pl.loglog(one_halo_ngp[:,0], one_halo_ngp[:,1], label="1halo using NGP")
#pl.loglog(one_halo_cic[:,0], one_halo_cic[:,1], label="1halo using CIC")
pl.loglog(one_halo[:,0], one_halo[:,1], label="1halo")
	
pl.legend(loc="best")

#data = np.loadtxt("GR_1fold.txt")
#pl.loglog(data[:,0], data[:,1]/8.) 
#volume = 1500.0 * 1500.0 * 1500.0
#number = 50000.0
#inv_no_dens = volume / number
#pl.axhline(y=inv_no_dens, linestyle="--")

pl.ylim(50.0, 1000.0)
#pl.xlim(0.02, 0.9)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./1halo_NGP_CIC.pdf")
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


