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
rows = file_len("./toaverage/randoms_onehalo_0.dat")
#print rows
data = []
monos = np.zeros(rows)


for i in xrange(0, simno):
	this_data = np.loadtxt("./toaverage/randoms_onehalo_" + str(i) + ".dat")
	#print data[:,1]
	monos += this_data[:,1]
	data.append(this_data)

ks = data[0][:,0]
monos /= simno

outFile = open("./toaverage/randoms_average.dat", "w")
for i in xrange(0, rows):
	outFile.write(str(ks[i]) + "\t" + str(monos[i]) + "\n")
outFile.close()
#np.savetxt("./toaverage/randoms_average.dat", (ks, monos), delimiter="\t", fmt="%1.4e")


