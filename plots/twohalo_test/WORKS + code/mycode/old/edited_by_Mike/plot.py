import numpy as np
import matplotlib.pyplot as pl

x=[]
y=[]

x, y = np.loadtxt("pk.txt", delimiter = ' ', unpack = True);
pl.plot(x, y);
pl.show();


