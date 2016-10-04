import numpy as np

matrices = []
expected = np.loadtxt("./cov_305_sim0.dat")

for i in xrange(3, 30):
	data = np.loadtxt("./cov_" + str(i) + "_sim0.dat")
	matrices.append(data)
	
def partialSum(n):
	return matrices[n]
	
def simpleCov(n):
	return partialSum(n)	
	
def shanks(n):
	if n < 1 or n > 25:
		return 0
	else:
		psn = partialSum(n)
		psn_minus = partialSum(n-1)
		psn_plus = partialSum(n+1)
		top = np.multiply(psn_plus,psn_minus) - np.multiply(psn, psn)
		bottom = psn_plus - 2*psn + psn_minus
		np.fill_diagonal(top, 1)		
		np.fill_diagonal(bottom, 1)
		inter = top/bottom
		return inter
		
def shanksCov(n):	
	return shanks(n)
	
def discrep(matr):
	inter = (matr-expected)/expected
	return abs(np.sum(inter))
		

for i in xrange(1, 26):
	print "simple calc: " + str(discrep(simpleCov(i)))
	print "shanks: " + str(discrep(shanksCov(i)))
	print "\n"
	
print shanksCov(25)

'''
print shanksCov(2)
print "\n"
print shanksCov(4)
'''





