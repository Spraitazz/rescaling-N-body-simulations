import numpy as np
import matplotlib.pyplot as pl
import os
import sys

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


simno = int(sys.argv[1])
rows = file_len("./F5_combined_1.dat")

data = []
data_rs = []
monos = np.zeros(rows)
quadros = np.zeros(rows)
monos_rs = np.zeros(rows)


for i in xrange(0, simno):
	this_data = np.loadtxt("./F5_combined_" + str(i+1) + ".dat")
	this_data_rs = np.loadtxt("./F5_combined_redshift_" + str(i+1) + ".dat")
	#print data[:,1]
	monos += this_data[:,1]
	monos_rs += this_data_rs[:,1]
	quadros += this_data_rs[:,2]
	data.append(this_data)
	data_rs.append(this_data_rs)

ks = data[0][:,0]
monos /= simno
monos_rs /= simno
quadros /= simno
#data_model = np.loadtxt("./model_pk_onehalo.dat")
#quadros = abs(quadros)

err_monos = np.zeros(rows)
err_monos_rs = np.zeros(rows)
err_quadros = np.zeros(rows)
err_zeros = np.zeros(rows)

#variances
for i in xrange(0, simno):
	err_monos += ((data[i][:,1] - monos)**2.0)/simno
	err_monos_rs += ((data_rs[i][:,1] - monos_rs)**2.0)/simno
	err_quadros += ((data_rs[i][:,2] - quadros)**2.0)/simno

#want stdevs
err_monos = err_monos**0.5
err_monos_rs = err_monos_rs**0.5
err_quadros = err_quadros**0.5

realspace_monos_out = open("./F5_realspace_mono_avg.dat", "w")
redshift_monos_out = open("./F5_redshift_mono_avg.dat", "w")
redshift_quadros_out = open("./F5_redshift_quadro_avg.dat", "w")

for i in xrange(0, rows):
	realspace_monos_out.write(str(ks[i]) + "\t" + str(monos[i]) + "\t" + str(err_monos[i]) + "\n")
	redshift_monos_out.write(str(ks[i]) + "\t" + str(monos_rs[i]) + "\t" + str(err_monos_rs[i]) + "\n")
	redshift_quadros_out.write(str(ks[i]) + "\t" + str(quadros[i]) + "\t" + str(err_quadros[i]) + "\n")	

realspace_monos_out.close()
redshift_monos_out.close()
redshift_quadros_out.close()



