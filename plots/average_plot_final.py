import numpy as np
import matplotlib.pyplot as pl
import os
import sys

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#count = 0
#while os.path.isfile("./power_spectra/power_spectrum" + str(count) + ".pdf"):
	#count+=1
#if os.path.isfile("./power_spectrum.pdf"):

simno = int(sys.argv[1])
rows = file_len("./final/randoms_final_0.dat")
#print rows
data = []
monos = np.zeros(rows)
quadros = np.zeros(rows)

for i in xrange(0, simno):
	this_data = np.loadtxt("./final/randoms_final_" + str(i) + ".dat")
	#print data[:,1]
	monos += this_data[:,1]
	quadros += this_data[:,2]
	data.append(this_data)

ks = data[0][:,0]
monos /= simno
quadros /= simno
data_model = np.loadtxt("./model_pk_final.dat")
#quadros = abs(quadros)

err_monos = np.zeros(rows)
err_quadros = np.zeros(rows)
err_zeros = np.zeros(rows)

#variances
for i in xrange(0, simno):
	err_monos += ((data[i][:,1] - monos)**2.0)/simno
	err_quadros += ((data[i][:,2] - quadros)**2.0)/simno

#want stdevs
err_monos = err_monos**0.5
err_quadros = err_quadros**0.5
	
ax = pl.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_title("Two-halo + one-halo terms for randoms")
#pl.
pl.errorbar(ks, monos, yerr=err_monos, label="one+two halo term monopole")
pl.loglog(data_model[:,0], data_model[:,1], label="model monopole")
#pl.errorbar(ks, quadros, xerr=err_zeros, yerr=err_quadros, label="quadrupole")
#data = np.loadtxt("GR_1fold.txt")
#pl.loglog(ks, monos, label="monopole") 
#pl.loglog(ks, quadros, label="quadrupole, interesting")
pl.legend(loc="best")
#volume = 1500.0 * 1500.0 * 1500.0
#number = 50000.0
#inv_no_dens = volume / number
#pl.axhline(y=inv_no_dens, linestyle="--")

pl.ylim(60.0,7500.0)
pl.xlim(0.02, 4.5)

pl.xlabel(r"$k, \ [h \rm{Mpc} ^{-1}]$")
pl.ylabel(r"$P(k), \ [(h^{-1} \rm{Mpc})^3]$")

pl.savefig("./final/final_average_" + str(simno) + ".pdf")
pl.clf()

#data2 = np.loadtxt("coef.txt")

#pl.plot(data2[:,0], data2[:,1])
#pl.savefig("coefficients.pdf")


