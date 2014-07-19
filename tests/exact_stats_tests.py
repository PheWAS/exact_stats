#!/usr/bin/env python

'''
File: exact_stats_tests.py
Created: 19-July-2014 Joseph P. Bochenek
1.0
Description: A python module for asssociation studies with rare phenomena. 
'''


import exact_stats

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


def test_basic():
	'''
	Do some tests
	'''

	a = [8, 5, 2, 4]
	
	mod = 10.
	a = [mod*a_ for a_ in a]
	a = [29403, 297, 297, 90]
	a = [29404, 298, 298, 7]

	a = [8, 5, 2, 4]

	bincenters = []

	print a
	for i in range(0,35):
		phi = -2+float(i)*0.2
		bincenters.append(phi)
#     	print "\n{}, {}/{}, {}".format(phi, exact_post(phi, a ),   exact_post_num(phi, a, 40 ), exact_pvalue_num(phi, a , 10 ))
	
	y = [exact_stats.exact_post( phi, a, True) for phi in [bin for bin in bincenters]]
 	z = [exact_stats.exact_pvalue_num( phi, a, 500) for phi in [bin for bin in bincenters]]
	x = [exact_stats.exact_post_num( phi, a, 500) for phi in [bin for bin in bincenters]]
	w = [integrate.quad( exact_post, -2, phi,  args=(a)) for phi in [bin for bin in bincenters]]
		
	for i in range(len(bincenters)):
		print "{}, {:.2f}/{:.2f}, cum.:{:.2f}/{:.2f}".format(bincenters[i],  y[i], x[i],  z[i], float(w[i][0]))

	total = np.sum([0.2*x1 for x1 in x ])
	print "Integral:", total

	import matplotlib.pyplot as plt

	fig = plt.figure()
	ax = fig.add_subplot(111)
 	l = ax.plot(bincenters, z, 'b--', linewidth=1)
	m = ax.plot(bincenters, y, 'g--', linewidth=1)
	n = ax.plot(bincenters, x, 'k--', linewidth=1)
	n = ax.plot(bincenters, w, 'r--', linewidth=1)

	ax.set_xlabel('')
	ax.set_ylabel('Probability')

	ax.grid(True)
	plt.savefig("test_cumulative.png")




