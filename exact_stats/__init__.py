import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpmath import *

import scipy.integrate as integrate
import scipy.special as special


# ---------------------------------------------------------------------------------------
# Latorre 1982
# ---------------------------------------------------------------------------------------

def B_(a, b):
# 	return np.prod(  [math.gamma(a_i) for a_i in a]  ) / math.gamma(np.sum( a )  )
# 	res =  np.sum([math.lgamma(a_i) for a_i in a]) - math.lgamma(np.sum( a )) 
# 	print np.sum([math.lgamma(a_i) for a_i in a]),  math.lgamma(np.sum( a )) , res, math.exp(res)
	res = special.betaln(a, b)
	return res



def F_num(a, b, c, x, k, verbose=False):
#  	print a, b, c, x
 	
 	coeff = math.lgamma(c)- (math.lgamma(a) + math.lgamma(b))
	total = 0
	prev = 100000
	i = 0
	tracker = 0
	
	if x < 10e-10:
		return 1.
	
	for j in range(k+1):
	
		i+=1
		tracker += 1
		liter = (math.lgamma(a+i) + math.lgamma(b+i) + i*math.log(x) ) - ( math.lgamma(c+i) + math.log(math.factorial(i))  ) + coeff 


		if (i > 2) and (total == 0.) and (tracker > 10):
			dif = liter - liter_prev
			if verbose:
				print "skipping forward", liter, dif, int((liter)/(2*dif))
			i += int((liter)/(2*dif))
			tracker = 0

		if verbose:
			print i, "liter", liter, "prev", prev, "\tboth", liter, "\ttotal:", total
# 		if total and prev and (total - prev < 0.0000000001):
# 			print "CONVERGENCE"
# 			break
		total = total + math.exp(liter) 
		prev = total
		liter_prev = liter
	return total



def F_int_approx(a, b, c, A, x, k, verbose=False):
#  	print a, b, c, x 	
	total = 0

	if x < 10e-10:
		return 0
	
	
	for i in range(1, k+1):
		if not( special.betainc(A + 1, i+1, x) ):
			beta_inc = 10e-320
		else:
			beta_inc = special.betainc(A + 1, i+1, x)
# 		print A+1, i+1, x, (math.lgamma(a+i) + math.lgamma(b+i) + math.lgamma(c) ) - (math.lgamma(a) + math.lgamma(b) + math.lgamma(c+i) + math.log(math.factorial(i))  ), special.betaln(A + 1, i+1)
		total += \
		math.exp( \
		(math.lgamma(a+i) + math.lgamma(b+i) + math.lgamma(c) ) - \
		(math.lgamma(a) + math.lgamma(b) + math.lgamma(c+i) + math.log(math.factorial(i))  ) +  \
		math.log( beta_inc ) + \
		special.betaln(A + 1, i+1) 
		)
	return total
	

def F_int(a, b, c, A, x, k, verbose=False):
#  	print a, b, c, x 	
	total = 0

	prev = 100000
	i = 0
	tracker = 0
	liter_prev = 0
	final_prec = 0
	l_conv_prev = 0.1

	if x < 10e-10:
		return 0
	
	liter_zero = 0
	
	PREC = 10e-5
	
	for j in range(1, k+1):

		i+=1
		tracker += 1
		
		# if the integral is too small it causes a domain error for math.log
		if not( special.betainc(A + 1, i+1, x) ):
			beta_inc = 10e-320
		else:
			beta_inc = special.betainc(A + 1, i+1, x)
		
		liter = \
		(math.lgamma(a+i) + math.lgamma(b+i) + math.lgamma(c) ) - \
		(math.lgamma(a) + math.lgamma(b) + math.lgamma(c+i) + math.log(math.factorial(i))  ) +  \
		math.log( beta_inc ) + \
		special.betaln(A + 1, i+1) 

		if i == 1:
			liter_zero = liter

		# skip the iterator ahead to save time
		dif =  liter - liter_prev
		
		
		total += math.exp( liter )

		if not(i%1000):
			print "i = {}".format(i)

		# Check for convergence
		if not(i%10):
# 			print i, liter, total, l_conv_prev, math.fabs((total - l_conv_prev)/l_conv_prev), PREC, math.fabs((total - l_conv_prev)/l_conv_prev) < PREC
			if (math.fabs((total - l_conv_prev)/l_conv_prev) < PREC) and (liter > -300):
 				print "Precision of {} reached after {} iterations.\t".format(PREC, i)
				return total
			l_conv_prev = total
					
		
		if (i > 2) and (tracker > 2) and (dif < 0) and (liter < -100):
#  			print "skipping forward", liter, dif, log((liter)/(2*dif)), 2*int(math.log((liter)/(2*dif)))
			i += min(100, 2*int(math.log((liter)/(dif))))
			tracker = 0




		liter_prev = liter
		prev = total

		if i > k:
			break
	
				
 	print "Precision of {} was reached after {} iterations.  \n".format(final_prec, i)
	return total



def K_(a):
	res = B_(a[0]+a[1], a[2]+a[3]) - (B_(a[0], a[2]) + B_(a[1],a[3]))
	return math.exp(res)




def F_(a, b, c, x): 
	maxterms = 100000
	try:
		hyp = mp.hyp2f1(a,b,c,x, maxterms=maxterms)
		return hyp
	except:
		pass
	try:
		hyp = mp.hyp2f1(a,b,c,x, maxterms=maxterms)
		return float(hyp)
	except:
		print "Warning: Both numpy and sympy fail to converge, returning zero for hypergeometric function, mpmath.hyp2f1 ({:.2f}, {:.2f}, {:.2f}, {}).".format(a,b,c,x, maxterms)
		return 0.


def exact_post(phi, a, verbose=False):
	a_1 = a[0] + a[2]
	a1_ = a[0] + a[1]
	a__ = np.sum(a)

	a_2 = a[1] + a[3]
	a2_ = a[2] + a[3]

	t = math.exp(phi)
	
	if phi > 0:
# 		res = K_(a) * t**(-1.*(a[1]+1)) *      F_(a1_, a_2, a__, 1-t**(-1.))
		res = K_(a) * t**(-1.*(a[1])) *      F_(a1_, a_2, a__, 1-t**(-1.))

#  		res = K_(a) * math.exp(-1.*phi*a[1]) * F_(a[0]+a[1], a[1]+a[3], np.sum(a), 1-math.exp(-1.*phi))
		if verbose:
			print "F({};{:.2e})={:.2e}".format(a, phi,  float(res))
		return float(res)
	else:
#		res = K_(a) * t**(a[0]-1.) * F_(a_1, a1_, a__, 1-t)
		res = K_(a) * t**(a[0]) * F_(a_1, a1_, a__, 1-t)

# 		res = K_(a) * exp(phi*a[0]) * F_(a[0]+a[2], a[0]+a[1], np.sum(a), 1-math.exp(phi))
		if verbose:
			print "F({};{:.2e})={:.2e}".format(a, float(phi),  float(res))
		return float(res)



def exact_post_num(phi, a, k, verbose=False):
	a_1 = a[0] + a[2]
	a1_ = a[0] + a[1]
	a__ = np.sum(a)

	a_2 = a[1] + a[3]
	a2_ = a[2] + a[3]

	t = math.exp(phi)
	
	if phi > 0:
# 		res = K_(a) * t**(-1.*(a[1]+1)) *      F_(a1_, a_2, a__, 1-t**(-1.))
		res = K_(a) * t**(-1.*(a[1])) * F_num(a1_, a_2, a__, 1-t**(-1.), k)

#  		res = K_(a) * math.exp(-1.*phi*a[1]) * F_(a[0]+a[1], a[1]+a[3], np.sum(a), 1-math.exp(-1.*phi))
		if verbose:
			print "F({};{:.2e})={:.2e}".format(a, phi,  float(res))
		return float(res)
	else:
#		res = K_(a) * t**(a[0]-1.) * F_(a_1, a1_, a__, 1-t)
		res = K_(a) * t**(a[0]) * F_num(a_1, a1_, a__, 1-t, k)

# 		res = K_(a) * exp(phi*a[0]) * F_(a[0]+a[2], a[0]+a[1], np.sum(a), 1-math.exp(phi))
		if verbose:
			print "F({};{:.2e})={:.2e}".format(a, float(phi),  float(res))
		return float(res)
		
		
		

def exact_pvalue_num(phi, a, k, verbose=False):
	a_1 = a[0] + a[2]
	a1_ = a[0] + a[1]
	a__ = np.sum(a)

	a_2 = a[1] + a[3]
	a2_ = a[2] + a[3]

	t = math.exp(phi)
	
	if phi > 0:
		res_sum = F_int(a1_, a_2, a__, a[1]-1, t**(-1.), k) 
		res =  1 - K_(a) * (  t**(-1.*(a[1]-1)) / (a[1]-1)  +  res_sum ) 
		print phi, a[1], K_(a), t**(-1.*(a[1]-1)) / (a[1]-1), res_sum
		if verbose:
			print "F({};{:.2e})={:.2e}".format(a, phi,  float(res))
		return 1.-float(res)
	# if phi < 0
	else:
		res_sum = F_int(a_1, a1_, a__, a[0]-1, t, k) 
		res = K_(a) * (t**(a[0]+1)/(a[0]+1) +  res_sum)

		if verbose:
			print "F({};{:.2e})={:.2e}".format(a, float(phi),  float(res))
		return float(res)
		





def exact_pvalue(phi, a, verbose=False):
	return integrate.quad( exact_post, phi, 4.*phi, args=(a))
	

def exact_pvalue2(phi, a, k, verbose=False):
	'''
	Return the exact p-value up to order k (Latorre 1982)
	'''
	t = math.exp(phi)

	a_1 = a[0] + a[2]
	a1_ = a[0] + a[1]
	a__ = np.sum(a)

	a_2 = a[1] + a[3]
	a2_ = a[2] + a[3]

	if t > 1.:
		F_0 = -1*t**(-1.*(a[1]-1.))/(a[1]-1)
		F_sum = 0
		
		for l in range(1, k+1):
			
			F_i =  ( math.lgamma(a1_ + l) + math.lgamma(a_2 + l) + math.lgamma(a__) ) - \
				   ( math.lgamma(a1_) + math.lgamma(a_2) + math.lgamma(a__ + l) + math.log(math.factorial(l)) ) 
			
			A = -1 * (-1)**(l) * 1./(l + a[1] - 1) * t**(1-l-a[1]) * F_(-l, -l-a[1]+1, -l-a[1]+2, t)
			
			F_sum += A * math.exp(F_i)
			
			if verbose:
				print "\n", l, phi, a[1], t
				print "\t (", math.lgamma(a1_ + l) , math.lgamma(a_2 + l) ,  math.lgamma(a__) , ") [", (math.lgamma(a1_ + l) + math.lgamma(a_2 + l) +  math.lgamma(a__) ) ,"] - \
					  \n\t(", math.lgamma(a1_)     , math.lgamma(a_2) ,      math.lgamma(a__ + l) ,    math.log(math.factorial(l)) , ") [", (math.lgamma(a1_) + math.lgamma(a_2) + math.lgamma(a__ + l) + math.log(math.factorial(l))),"] + \n\t(", B_(a[1], l+1) , math.log(special.betainc(a[1], l+1, 1./t)), ")"
				print "\t F_{} = {} \t exp(F_i) = {} \t F_sum = {}".format(l, F_i, math.exp(F_i) , F_sum)
				
				print "\tI_{} = {}".format(l, 1 + K_(a) * ( F_0 + F_sum  ))

		return 1 + K_(a) * ( F_0 + F_sum )



	else:
		F_0 = t**(a[0]+1.)/(a[0]+1)
		F_sum = 0

		for i in range(1, k+1):
			F_i = (math.lgamma(a_1 + i) + math.lgamma(a1_ + i) + math.lgamma(a__) ) - \
	   			  (math.lgamma(a_1) + math.lgamma(a1_) + math.lgamma(a__ + i) + math.log(math.factorial(i)) ) + \
			      B_(a[0]+1, i+1) 
			F_sum +=  math.exp(F_i)

		return K_(a) * (F_0 + F_sum)



def test_run():
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
	
	y = [exact_post( phi, a, True) for phi in [bin for bin in bincenters]]
 	z = [exact_pvalue_num( phi, a, 500) for phi in [bin for bin in bincenters]]
	x = [exact_post_num( phi, a, 500) for phi in [bin for bin in bincenters]]
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




