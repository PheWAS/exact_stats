#!/usr/bin/env python

'''
File: __init__.py
Created: 19-July-2014 Joseph P. Bochenek
1.0
Description: A python module for asssociation studies with rare phenomena. 
'''


import math
import numpy as np
import scipy.stats as stats

import scipy.integrate as integrate
import scipy.special as special
from prettytable import PrettyTable
import mpmath
mpmath.mp.dps = 100

def b_l(l, a):
	'''
	Compute the l'th cumulant for the posterior distribution of the log-odds ratio with contingency table a
	
	Reference:
	'A new approximation of the posterior distribution of the log-odds ratio' - Fredette and Angers, 2002
	'''
	a = [int(a_i) for a_i in a]
	C = [1.,-1.,-1.,1.]
	if l == 1:
		return np.sum( [ c * special.psi(a_j) for c, a_j in zip(C, a) ] )
	else:
		return math.factorial(l-1) * np.sum( [ (-1. * c)**l * special.zeta(l, a_j) for c, a_j in zip(C, a) ]  ) 

def B_(a, b):
	'''
	log Beta function
	'''

	res = special.betaln(a, b)
	return res


def F_func(i, a, b, c, x):
	liter = math.lgamma(a+i) - math.lgamma(a) \
		  + math.lgamma(b+i) - math.lgamma(b) \
		  + math.lgamma(c) - math.lgamma(c+i) \
		  + i*math.log(x) - math.lgamma(i+1)
	
	liter_pre = liter
	if liter > 709.:
		liter = 709.
	return math.exp(liter)



def F_num(a, b, c, x, verbose=False):
	if verbose:
		print '\nComputing F_num():'
		print 'F_num({}, {}, {}, {}, verbose={})'.format(a, b, c, x, verbose)

	ans, err = integrate.quad( F_func, 0, 1e3, (a,b,c,x))

	# if the integrand is greater than 1., we need to modify the upper limit of integration
	relerr = 1.
	ind = 1.
	prev = ans

	while (relerr > 10e-3):
		upper = 10e3 * (ind)
		if verbose:
			print 'raising integration limit to', upper, ans, err
		if ind > 20:
			if verbose:
				print 'F_num(): ans > 1000, no convergence'
			break 	
		ans, err = integrate.quad( F_func, 0, upper, (a,b,c,x))
		ind += 1

		if prev == 0:
			relerr = 1.
		else:
			relerr = math.fabs(ans-prev)/prev
		prev = ans

	# if the answer is small (<1.
	else:
		if err/ans > 10e-2:
			if verbose: print 'No convergence, increasing precision to 10e4 iters', ans, err
			ans, err = integrate.quad( F_func, 0, 10e4, (a,b,c,x))
		if ans==0:
				if verbose: print 'No convergence, increasing precision to 20e5 iters', ans, err
				ans, err = integrate.quad( F_func, 0, 10e5, (a,b,c,x))
		else:
			if (err/ans > 10e-2):
				if verbose: print 'No convergence, increasing precision to 20e5 iters', ans, err
				ans, err = integrate.quad( F_func, 0, 10e5, (a,b,c,x))

		if not(ans==0):
			if err/ans > 10e-2 or math.isnan(ans):
				if verbose: print 'Numerical integration did not converge.  Relative error:', err/ans

	if verbose:
		print 'F=', ans, '+-', err
		
		

	return ans


def F_num_exp(a, b, c, x, k, prec=10e-4, verbose=False):
	'''
	Performs the work for exact_post_num
	'''
	if verbose:
		print 'F_int({}, {}, {}, {}, {}, prec=10e-4, verbose=False):'.format(a, b, c, x, k)

	coeff = math.lgamma(c)- (math.lgamma(a) + math.lgamma(b))
	total = 0
	prev = 100000
	i = 0
	tracker = 0
	l_conv_prev = 0.1

	if x < 10e-10:
		return 0.
	
	for j in range(k+1):
	
		i+=1
		tracker += 1
		liter = (math.lgamma(a+i) + math.lgamma(b+i) + i*math.log(x) ) - ( math.lgamma(c+i) + math.log(math.factorial(i))  ) + coeff 

		# Check for convergence
		if not(i%10):
			final_prec = math.fabs((total - l_conv_prev)/l_conv_prev)
			if (final_prec < prec) and (liter > -300):
				if verbose:
	 				print "Precision of {} was reached after {} iterations.\t".format(prec, i)
				return (total, final_prec)
# 			print i, liter, total, l_conv_prev, final_prec, dif
			l_conv_prev = total

		if (i > 2) and (total == 0.) and (tracker > 10):
			dif = liter - liter_prev
			if verbose:
				print "skipping forward", liter, dif, int((liter)/(2*dif))
			i += int((liter)/(2*dif))
			tracker = 0

		if verbose:
			if not(i%1000):
				print "i = {}, prec = {}/{}, liter = {}".format(i, final_prec, prec, total)

		total = total + math.exp(liter) 
		prev = total
		liter_prev = liter
		
	if verbose:		
	 	print "Reached the maximum interations.  Precision of {} was reached after {} iterations.  \n".format(final_prec, i)

	return (total, final_prec)


def F_int_func2(i, a, b, c, A, x, verbose=False):
	
# 	beta_inc = special.betainc(A, i+1, x)

# 	if not( special.betainc(A, i+1, x) ):
# # 		if verbose:
# 		print 'parameters', i, a, b, c, A, x
# 		print 'rounding error 1', A, i+1, x, beta_inc, mpmath.betainc(A, i+1, x)
# 		return 0.

# 	beta_inc = special.betainc(A, i+1, x)

	try:
		lbeta_inc = mpmath.mp.log( mpmath.mp.betainc(A, i+1, 0, x))
	except:
		return 0
# 	print lbeta_inc
	liter = mpmath.mp.loggamma(a+i) - mpmath.mp.loggamma(a) \
		  + mpmath.mp.loggamma(b+i) - mpmath.mp.loggamma(b) \
		  + mpmath.mp.loggamma(c) - mpmath.mp.loggamma(c+i) \
		  - mpmath.mp.loggamma(i+1) \
		  + lbeta_inc \
# 		  + special.betaln(A, i+1) 
		  
# 	print mpmath.mp.loggamma(a+i), '-', mpmath.mp.loggamma(a), '/', math.lgamma(a+i), '-', math.lgamma(a)
# 	print '+', mpmath.mp.loggamma(b+i), '-', mpmath.mp.loggamma(b) 
# 	print '+', mpmath.mp.loggamma(c), '-', mpmath.mp.loggamma(c+i) 
# 	print '-', mpmath.mp.loggamma(i+1) 
# 	print '+lbeta', lbeta_inc , '/', math.log( beta_inc ) + special.betaln(A, i+1) 
# 	
# 	print liter
#  	print i, x, beta_inc, liter, mpmath.mp.exp(liter).real, float(mpmath.mp.exp(liter).real)
	return mpmath.mp.exp(liter)


def F_int_func(i, a, b, c, A, x, verbose=False):
	beta_inc = special.betainc(A, i+1, x)
	if not(beta_inc ):
		alt = F_int_func2(i, a, b, c, A, x, verbose=False)
		if verbose:
			print 'rounding error 1', i, beta_inc, special.betaln(A, i+1), mpmath.betainc(A, i+1, x), alt
		return alt

		beta_inc = 10e-320

	liter = math.lgamma(a+i) - math.lgamma(a) \
		  + math.lgamma(b+i) - math.lgamma(b) \
		  + math.lgamma(c) - math.lgamma(c+i) \
		  - math.lgamma(i+1) \
		  + math.log( beta_inc ) \
		  + special.betaln(A, i+1) 
		  
	if liter > 709.:
# 		if verbose:
		print 'rounding error 2', i, 
		liter = 709.
	return math.exp(liter)



def F_int2(a, b, c, A, x, max=mpmath.mp.inf, prec=100, verbose=False):	
 	mpmath.mp.dps = prec
	if verbose:
		print '\nComputing F_num():'
		print 'F_int({}, {}, {}, {}, {}, verbose={}, max={}, prec={})'.format(a, b, c, A, x, verbose, max, mpmath.mp.dps)

	f = lambda i: F_int_func2(i, a, b, c, A, x)
# 	ans, err = integrate.quad( F_int_func, 1, 1e3, (a, b, c, A, x))
	ans, err =  mpmath.mp.quad(f, [1, max], error=True, verbose=verbose) #mpmath.mp.inf
	if verbose:
		print "Result", ans, err
	return ans, err


def F_int(a, b, c, A, x,  verbose=False):	
	if verbose:
		print '\nComputing F_num():'
		print 'F_int({}, {}, {}, {}, {}, verbose={})'.format(a, b, c, A, x, verbose)

	ans, err = integrate.quad( F_int_func, 1, 1e3, (a, b, c, A, x))

	# if the integrand is greater than 1., we need to modify the upper limit of integration
	relerr = 1.
	abserr = 1.

	ind = 1.
	prev = ans

	while ((relerr > 10e-3) or (abserr > 10e-2)):
		
		upper = 2e2 * 2**(ind)

		if ind > 9:
			if verbose:
				print 'ans > 1000, no convergence'
			break 	
		ans, err = integrate.quad( F_int_func, 1, upper, (a,b,c,A,x))
		ind += 1

		if ans == 0:
			abserr = 1.
		else:
			abserr = err/ans

		if prev == 0:
			relerr = 1.
		else:
			relerr = math.fabs(ans-prev)/prev
		prev = ans

		if verbose:
			print 'raising integration limit to {} -- relerr/abserr:{:.2e}/{:.2e}'.format( upper, relerr, abserr )

	

	# Try second method (using mpmath for high precision integration, much slower)
# 	if not(ans):
# 		try:
# 			print 'trying mpmath'
# 			ans,err = F_int2(a, b, c, A, x,  verbose=False)
# 		except:
# 			print 'mpmath didn\'t work either'
# 	elif err/ans > 10e-2:
# 		try:
# 			print 'trying mpmath'
# 			ans,err = F_int2(a, b, c, A, x,  verbose=False)
# 		except:
# 			print 'mpmath didn\'t work either'

	if verbose:
		print 	'Result', ans , '+=', err 
	if ans:
		if verbose:
			print ' error margin:', err/ans, '/', 10e-3
		if err:
			if err/ans > 0.1:
				print 'Numerical integration did not converge.  Relative error:', err/ans
				return float('nan')
		else:
			return float('nan')
	print "Answer:", ans
	return ans


def F_int_exp(a, b, c, A, x, k, prec=10e-4, verbose=False):
	'''
	Performs the work for pdf2.
	
	Reference:
	Fredette/Angers 2002
	'''
	
	

	total = 0

	prev = 100000
	i = 0
	tracker = 0
	liter_prev = 0
	final_prec = 0
	l_conv_prev = 0.1

	if x < 10e-10:
		return 0
	
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

		# Skip the iterator ahead to save time
		dif =  liter - liter_prev
		total += math.exp( liter )

		if verbose:
			if not(i%1000):
				print "i = {}, prec = {}/{}".format(i, final_prec, prec)

		# Check for convergence
		if not(i%10):
			final_prec = math.fabs((total - l_conv_prev)/l_conv_prev)
			if (final_prec < prec) and (liter > -300):
				if verbose:
	 				print "Precision of {} was reached after {} iterations.\t".format(prec, i)
 				return (total, final_prec)
# 			print i, liter, total, l_conv_prev, final_prec, dif
			l_conv_prev = total
					
		
		if (i > 2) and (tracker > 2) and (dif < 0) and (liter < -100):
			i += min(100, 2*int(math.log((liter)/(dif))))
# 			print '\t', i, 'skipping ahead by ', min(100, 2*int(math.log((liter)/(dif))))
			tracker = 0

		if (final_prec > 1.) and (tracker > 2):
			i += 100
# 			print '\t', i, 'skipping ahead by ', 100
			tracker = 0


		liter_prev = liter
		prev = total

		if i > k:
			break

	if verbose:		
	 	print "Reached the maximum interations.  Precision of {} was reached after {} iterations.  \n".format(final_prec, i)

	return (total, final_prec)

def K_(a):
	'''
	Compute the K function
	'''
# 	res = B_(a[0]+a[1], a[2]+a[3]) - (B_(a[0], a[2]) + B_(a[1],a[3]))
	
	return mpmath.beta(a[0]+a[1], a[2]+a[3])/(mpmath.beta(a[0], a[2]) * mpmath.beta(a[1],a[3]))
# 	return math.exp(res)

def F_(a, b, c, x, verbose=False): 
	'''
	Compute the F function using two different algorithms which have different regions of convergence.
	'''
	if verbose:
		print '\nComputing F_():'
		print 'F_({}, {}, {}, {}, verbose={})'.format(a, b, c, x, verbose)

	from mpmath import hyp2f1 as mp_hyp2f1
	from scipy.special import hyp2f1 as np_hyp2f1
	maxterms = 20000

	hyp = 0.
	try:
		hyp = np_hyp2f1(a,b,c,x)
		if verbose:
			print 'trying Numpy, F=', hyp
	except:
		if verbose:
			print 'Warning: np doesnt converge'
		try:
			hyp = mp_hyp2f1(a,b,c,x, maxterms=maxterms)
			if verbose:
				print 'trying MP, F=', hyp
			
		except:
			print "Warning: Both numpy and sympy fail to converge, returning zero for hypergeometric function, mpmath.hyp2f1 ({:.2f}, {:.2f}, {:.2f}, {}).".format(a,b,c,x, maxterms)
	
	if math.isnan(hyp) or math.isinf(hyp):
		hyp = F_num(a, b, c, x, verbose=verbose)
		if verbose:
			print 'trying numerical calculation, F=', hyp


	return hyp



def pi_Z(m, a, z, verbose=False):
	'''
	This implements the edgeworth expansion recusion relation from Fredette/Angers 2002
	'''
	D = [0]*(m+1)
	M = [0]*(m+1)

	D[0] = 1./math.sqrt(2.*math.pi)
	D[1] = 0
	D[2] = 0
	
	sigma = math.sqrt(b_l(2, a))


	if m > 2:
		for k in range(3,m+1):
			D_sum = 0
			for l in range(3, k+1):
				coeff = 1. / ( k * math.factorial(l-1) * (-1.*sigma)**(l) )
				D_iter = D[k-l] * b_l(l, a)
				D_sum += coeff * D_iter
			D[k] = D_sum

	M[0] = math.exp(-0.5 * z**2)
	M[1] = -1. * z * math.exp(-0.5 * z**2)

	if m > 1:
		for k in range(2,m+1):
			M[k] = -1. * z * M[k-1] - (k-1) * M[k-2]

	if verbose:
		print ""
		for k in range(m+1):
			print "D_{}(a) = {:.3e} \t M_{}({:.2f}) = {:.2f} \t DM = {:.2e}".format(k, D[k], k, z, M[k], D[k]*M[k])

	return D, M

def pdf(phi, a, verbose=False):
	'''
	Compute the posterior of the log-odds-ratio
	
	Description: Based on approximation in [Latorre, 1982]
	
	Input:
	phi - log-OR to calculate the p-value
	a - 2x2 contingency table

	Keyword Arguments:
	verbose - be loud
	
	Returns: The posterior value at phi
	'''



	t = math.exp(phi)
	
	if phi <= 0:
		a = [a[3], a[2], a[1], a[0]]

		a_1 = a[0] + a[2]
		a1_ = a[0] + a[1]
		a__ = np.sum(a)

		F_res = F_(a_1, a1_, a__, 1-t, verbose=verbose)
		res = K_(a) * t**(a[0]) * F_res
		if verbose:
			print 'phi:{:.2e}\na:{}'.format( phi, a )
			print "K() F() t**(-1.*({}))= {:.3e} * {:.3e} * {:.3e} = {:.3e}".format(a[1], K_(a), F_res, t**((a[0])), float(res))
		return float(res)

	else:
		a1_ = a[0] + a[1]
		a__ = np.sum(a)
		a_2 = a[1] + a[3]

		F_res = F_(a1_, a_2, a__, 1-t**(-1.), verbose=verbose)

		# Algorithm diverges for large arguments, but it is symmetric under this transformation:
		if math.isnan(F_res):
			a = [a[3], a[2], a[1], a[0]]
			a1_ = a[0] + a[1]
			a__ = np.sum(a)
			a_2 = a[1] + a[3]
			F_res = F_(a1_, a_2, a__, 1-t**(-1.), verbose=verbose)

		res = K_(a) * t**(-1.*(a[1])) * F_res
		if verbose:
			print '\nphi:{}\na:{}\nt={}'.format( phi, a, t )
			print "K() F() t**(-1.*({}))= {:.3e} * {:.3e} * {:.3e} = {:.3e}".format(a[1], K_(a), F_res, t**(-1.*(a[1])), float(res))
		return float(res)



def pval_norm(a, verbose=False):

	a = normalize_table(a)

	prior = [1., 1., 1., 1.]
	prior = [0.5, 0.5, 0.5, 0.5]
	prior = [0., 0., 0., 0.]

	odds_ratio = (float(a[1][1])/float(a[1][0]))/(float(a[0][1])/float(a[0][0]))
	lOR = math.log(odds_ratio)

	# Parameters for the posterior of the null hypothesis 	
	a = [a[0][0], a[0][1], a[1][0], a[1][1]]
	a_null = [a[0]+prior[0], a[1]+prior[1], a[2]+prior[2], float(a[2])*float(a[1])/float(a[0]+prior[3])]

	mu = b_l(1, a_null)
	sigma = math.sqrt(b_l(2, a_null))

	C = [1.,-1.,-1.,1.]

	mu    = np.sum(  [  c*(math.log(j) - 1./(2.*j)) for c, j in zip(C, a_null) ]  )
	sigma = math.sqrt(	np.sum(  [  1./(j) for c, j in zip(C, a_null) ]  )	)

	mu    = np.sum(  [  c*(math.log(j) - 1./(2.*j)) for c, j in zip(C, a) ]  )
	sigma = math.sqrt(	np.sum(  [  1./(j) for c, j in zip(C, a) ]  )	)
	
	mpmath.mp.dps=100

	if verbose:
		print 'log(OR)', lOR
		print 'mu, sigma', mu, sigma
		print 'numpy:', 1.-stats.norm.cdf(0., mu, sigma)
		print 'mpmath:', 1.-mpmath.mp.ncdf(0., mu=mu, sigma=sigma)

	if lOR > 0:
		pval = float(str(mpmath.mp.ncdf(0., mu=mu, sigma=sigma)))
		return odds_ratio, pval
	else:
		pval = float(str(1. - mpmath.mp.ncdf(0., mu=mu, sigma=sigma)))
		return odds_ratio, pval


def cdf(phi, a, k=10000, prec=10e-4, verbose=False):
	'''
	Compute the p-value using the posterior distribution for the log-odds ratio
	
	Based on approximation in [Latorre, 1982] (corresponds to pdf())
	
	Input:
	phi - log-OR to calculate the p-value
	a - 2x2 contingency table

	Keyword Arguments:
	k - maximum number of terms for the numerical calculation
	prec - desired precision to terminate the loop
	
	Returns: The p-value
	'''

	t = math.exp(phi)
	PREC = -1
	
	if phi < 0:
		a_r = [a[3], a[2], a[1], a[0]]
		a_r = [a[0], a[1], a[2], a[3]]

		a_1 = a_r[0] + a_r[2]
		a1_ = a_r[0] + a_r[1]
		a__ = np.sum(a_r)

		F_res = F_int(a_1, a1_, a__, a_r[0], t, verbose=verbose) 
		res = K_(a_r) * (  t**(a_r[0]+1.)/(a_r[0]+1.) + F_res )

		if verbose:
			print '\nphi:{}\na:{}\nt={}\nF_res={}'.format( phi, a_r, t, F_res )

			print "K_(a_r) * (  t**(a_r[0]+1.)/(a_r[0]+1.) + F_res ) = {} * ({} + {}) ".format(float(K_(a_r)), 			(  t**(a_r[0]+1.)/(a_r[0]+1.)), 			float(F_res) )
		print phi, res
		return float(res)

	else:
		a1_ = a[0] + a[1]
		a__ = np.sum(a)
		a_2 = a[1] + a[3]
		F_res, err = F_int2(a1_, a_2, a__, a[1], t**(-1.), verbose=verbose, prec=100)

		# Algorithm diverges for large arguments, but it is symmetric under this transformation:
		if math.isnan(F_res):
			if verbose:
				print "Inverting matrix"
			a = [a[3], a[2], a[1], a[0]]
			a1_ = a[0] + a[1]
			a__ = np.sum(a)
			a_2 = a[1] + a[3]
			F_res, err = F_int2(a1_, a_2, a__, a[1], t**(-1.), verbose=verbose)

		if math.isnan(F_res):
			return F_res
		
		res =  1 - K_(a) *  (  mpmath.mpf(t**(-1.*(a[1])) / (a[1])) + F_res )
		if verbose:
			print '\nphi:{}\na:{}\nt={}\tF_res={}'.format( phi, a, t, F_res )
			print phi, a[1], K_(a), t**(-1.*(a[1])) / (a[1]), F_res
			print a[1], mpmath.mpf(t**(-1.*(a[1])) / (a[1]))
			print "K_(a_r) * (  t**(a_r[0]+1.)/(a_r[0]+1.) + F_res ) = 1. - {} * ({} + {}) = 1.-{}".format(float(K_(a)),  mpmath.mpf(t**(-1.*(a[1])) / (a[1])), F_res, ( K_(a)* (t**(-1.*(a[1])) / (a[1])+F_res  )))
		return float(res)


def pdf2(phi, a, k=10000, prec=10e-4,  verbose=False):
	''' 
	Same as pdf but uses F_num() to compute the hypergeometic function instead of the 
	built-in functions in scipy.  Can be useful for debugging and evaluating the approximation.
	'''

	t = math.exp(phi)
	
	if phi > 0:
		a1_ = a[0] + a[1]
		a__ = np.sum(a)
		a_2 = a[1] + a[3]

		F_res = F_num(a1_, a_2, a__, 1-t**(-1.), verbose=verbose)
		res = K_(a) * t**(-1.*(a[1])) * F_res
		if verbose:
			print '\nphi:{}\na:{}'.format( phi, a )
			print "K() F() t**(-1.*(a[1]))= {:.3e} * {:.3e} * {:.3e} = {:.3e}".format(K_(a), F_res, t**(-1.*(a[1])), float(res))
		return float(res)

	else:
		a_r = [a[3], a[2], a[1], a[0]]

		a_1 = a_r[0] + a_r[2]
		a1_ = a_r[0] + a_r[1]
		a__ = np.sum(a_r)

		F_res = F_num(a_1, a1_, a__, 1-t, verbose=verbose)
		res = K_(a_r) * t**(a_r[0]) * F_res
		if verbose:
			print 'phi:{:.2e}\na:{}'.format( phi, a )
			print "K() F() t({})**(-1.*(a[1]))= {:.3e} * {:.3e} * {:.3e} = {:.3e}".format(t, K_(a_r), F_res, t**((a_r[0])), float(res))
		return float(res)

def pdf3(phi, a, m=4, verbose=False):
	'''
	This implementation uses the Edgeworth expansion to compute the pdf

	Reference:
	'A new approximation of the posterior distribution of the log-odds ratio' - Fredette and Angers, Statistica Neerlandica 2002
	'''

	for a_i in a:
		assert (float(a_i).is_integer()), "Error: this implementation only works with integers, {}".format(a)

	a = [int(a_i) for a_i in a]

	if verbose:
		print "mu: {:.2f}, sigma: {:.2f}".format(b_l(1, a), b_l(2, a))
	z =  (phi - b_l(1,a)) / b_l(2, a)**(0.5)
	D, M = pi_Z(m, a, z, verbose)
	res = np.sum( [m_*d_ for m_,d_ in zip(M,D)] )/b_l(2, a)**0.5

	return res

def pdf4(phi, a, m=4, verbose=False):
	'''
	This implementation expands the posterior distribution in terms of the standard normal
	distribution using the Gram-Charlier series (http://en.wikipedia.org/wiki/Edgeworth_series).
	This series generally doesn't converge, but can be really fast.  Generally, don't use more 
	than four or five terms (m=4).  This makes the distribution accurate to the fourth moment, 
	the kurtosis and is much better than using a normal distribution (which is the first two moments). 
	
	Actually, just use pdf3.
	
	References:
	https://indico.cern.ch/event/74919/session/16/contribution/88/material/slides/1.pdf
	'''

	sigma = math.sqrt(b_l(2,a))
	mu = b_l(1,a)
	z = (phi - mu)/sigma

	C = [1.,-1.,-1.,1.]

	normpart = stats.norm.pdf(phi, loc=mu, scale=sigma)
	seriespart = [  b_l(k,a)/(math.factorial(k) * sigma**(k) )  *  special.eval_hermitenorm(k, z)[0] for k in range(3,m+1) ]

	return normpart * (1 + np.sum(seriespart))


def positive_OR(a, verbose=True):
	A = a
	a = np.array(A)
	
	assert (a.shape == (2,3)) or (a.shape == (2,2)), "Error: input must be a 2x2 or 2x3 contingency table, {}".format(a)

	a_flip = []
	if a.shape == (2,3):
		a_flip = [ [a[0][2], a[0][1], a[0][0]], [a[1][2], a[1][1], a[1][0]] ]
	if a.shape == (2,2):
		a_flip = [ [ a[0][1], a[0][0]], [ a[1][1], a[1][0]] ]

	return a_flip

def normalize_table(a, verbose=False):
	A = a
	a = np.array(A)
	
	assert (a.shape == (2,3)) or (a.shape == (2,2)), "Error: input must be a 2x2 or 2x3 contingency table, {}".format(a)
		
	# Convert the contingency table from (people x exposure) to (chromosome x exposure), convert 2x3 to 2x2
	if a.shape == (2,3):
		x = PrettyTable(["Minor Allele", "X=0", "X=1", "X=2"])
		x.add_row(["Control", a[0][0], a[0][1], a[0][2]] )
		x.add_row(["Cases",  a[1][0],  a[1][1], a[1][2]] )
		if verbose:
			print "\nConverting 3x2 to 2x2 contingency table."
			print x
			print "--->"

		a = np.array([[2.*a[0][0]+a[0][1], a[0][1]+2.*a[0][2]],[2.*a[1][0]+a[1][1], a[1][1]+2.*a[1][2]]])
		x = PrettyTable(["Allele", "X=a", "X=A"])
		x.add_row(["Control", a[0][0], a[0][1]] )
		x.add_row(["Cases",  a[1][0],  a[1][1]] )
		if verbose:
			print x
			print ""
		
	else:
		if verbose:
			print "\nContingency Table:"
		x = PrettyTable(["Allele", "X=a", "X=A"])
		x.add_row(["Control", a[0][0], a[0][1]] )
		x.add_row(["Cases",  a[1][0],  a[1][1]] )
		if verbose:	
			print x
			print ""

	
	return a
	

def pval_null(a, k=10000, prec=10e-4, verbose=False):

	'''
	Compute the p-value for the null hypothesis given a contingency table
	
	Input:
	a - a 2x2 or 2x3 contingency table (assuming an additive model)

	Keyword Arguments:
	k - maximum number of terms for the numerical calculation
	prec - desired precision to terminate the loop	

	Return Value:
	a tuple with the OR and p-value
	'''

	if len(a[0])<2:
		return (-1, -1)
	if int(a[0][1]) == 0 or int(a[1][1]) == 0:
		return (-1, -1)


	a = normalize_table(a)
	
	

	odds_ratio = (float(a[1][1])/float(a[1][0]))/(float(a[0][1])/float(a[0][0]))


	if not(odds_ratio):
		return (0., 0.)
	phi = math.log(odds_ratio)
	if verbose:
		print "Odds Ratio: {:.3f}".format(odds_ratio)
		print "log(OR): {:.3f}".format(phi)
	
	# Jeffries Prior
	prior = 0.5
	
	a = [item +prior for sublist in a for item in sublist]	
	
	# Jeffrey's Prior
	prior = [1., 1., 1., 1.]
	prior = [0.5, 0.5, 0.5, 0.5]
	prior = [0., 0., 0., 0.]
	
	
	# Parameters for the posterior of the null hypothesis 	
	a_null = [a[0]+prior[0], a[1]+prior[1], a[2]+prior[2], float(a[2])*float(a[1])/float(a[0])]

	if verbose:
		print '\nPrior parameter table:\t', 
		print '\nNull model parameter table:', a_null
		print '\nMean: {:.3e}, Sigma: {:.3e}'.format(b_l(1, a_null), b_l(1, a_null)**0.5 )

	pval = 1.
	if phi > 0:
		pval = cdf(0., a, k=k, prec=prec, verbose=verbose)
	else:
		pval = 1.-cdf(0., a, k=k, prec=prec, verbose=verbose)

	if verbose:
		print "Null hypothesis P-value = {:.3e}\n".format(pval)

	return (odds_ratio, pval)


def odds_ratio(a):
	a = normalize_table(a)
	odds_ratio = (float(a[1][1])/float(a[1][0]))/(float(a[0][1])/float(a[0][0]))
	return (odds_ratio)


def Bayes_factor(a, k=10000, prec=10e-4, verbose=False):

	'''
	Compute the p-value for the null hypothesis given a contingency table
	
	Input:
	a - a 2x2 or 2x3 contingency table (assuming an additive model)

	Keyword Arguments:
	k - maximum number of terms for the numerical calculation
	prec - desired precision to terminate the loop	

	Return Value:
	a tuple with the OR and p-value
	'''
	

	A = a
	a = np.array(A)
	
	assert (a.shape == (2,3)) or (a.shape == (2,2)), "Error: input must be a 2x2 or 2x3 contingency table, {}".format(a)
		
	# Convert the contingency table from (people x exposure) to (chromosome x exposure), convert 2x3 to 2x2
	if a.shape == (2,3):
		x = PrettyTable(["Minor Allele", "X=0", "X=1", "X=2"])
		x.add_row(["Control", a[0][0], a[0][1], a[0][2]] )
		x.add_row(["Cases",  a[1][0],  a[1][1], a[1][2]] )
		if verbose:
			print "\nConverting 3x2 to 2x2 contingency table."
			print x
			print "--->"

		a = np.array([[2.*a[0][0]+a[0][1], a[0][1]+2.*a[0][2]],[2.*a[1][0]+a[1][1], a[1][1]+2.*a[1][2]]])
		x = PrettyTable(["Allele", "X=a", "X=A"])
		x.add_row(["Control", a[0][0], a[1][0]] )
		x.add_row(["Cases",  a[0][1],  a[1][1]] )
		if verbose:
			print x
			print ""
		
	else:
		if verbose:
			print "\nContingency Table:"
		x = PrettyTable(["Allele", "X=a", "X=A"])
		x.add_row(["Control", a[0][0], a[1][0]] )
		x.add_row(["Cases",  a[0][1],  a[1][1]] )
		if verbose:	
			print x
			print ""
				

	odds_ratio = (float(a[1][1])/float(a[1][0]))/(float(a[0][1])/float(a[0][0]))


	if not(odds_ratio):
		return (0., 0., 0.)
	phi = math.log(odds_ratio)
	if verbose:
		print "Odds Ratio: {:.3f}".format(odds_ratio)
		print "log(OR): {:.3f}".format(phi)
	
	a = [item for sublist in a for item in sublist]	
	
	# Jeffrey's Prior
	prior = [1., 1., 1., 1.]
	prior = [0.5, 0.5, 0.5, 0.5]
	
	# Parameters for the posterior of the null hypothesis 	
	a_null = [a[0]+prior[0], a[1]+prior[1], a[2]+prior[2], float(a[2])*float(a[1])/float(a[0])+prior[3]]
	a_alt =  [a[0]+prior[0], a[1]+prior[1], a[2]+prior[2], a[3]+prior[3]]

	if verbose:
		print '\nNull model parameter table:\n', a_null
		print 'Mean: {:.3e}, Sigma: {:.3e}'.format(b_l(1, a_null), b_l(1, a_null)**0.5 )

		print '\nAlternat model parameter table:\n{}'.format(a_alt)
		print 'Mean: {:.3e}, Sigma: {:.3e}\n'.format(b_l(1, a_alt), b_l(1, a_alt)**0.5 )

	pval = 1.
	
	H0 = pdf(phi, a_null, verbose=verbose)
	print 'H0', H0
	H1 = pdf(phi, a_alt, verbose=verbose)
	print 'H1', H0

	if verbose:
		print "H1/H0 posterior value = {:.3e} / {:.3e}\n".format(H1, H0)
		print 

	bayes_factor = H1/H0

		
	return (bayes_factor, odds_ratio)

def cdf_norm(phi, a):
	sigma = math.sqrt(b_l(2,a))
	mu = b_l(1,a)
	z = (phi - mu)/sigma
	cdf = stats.norm.cdf(z)
	return cdf

def cdf2(phi, a, m=4, verbose=False):
	'''
	Computes the p-value corresponding to the method of pdf3()
	
	[1] eq. 10
	'''
	sigma = math.sqrt(b_l(2,a))
	mu = b_l(1,a)

	z = (phi - mu)/sigma

	D, M = pi_Z(m, a, z, verbose)

	cdf = stats.norm.cdf(z) + np.sum([ D[k] * M[k-1] for k in range(1,m+1)])

	return cdf

def cdf_num(phi, a, verbose=False):
	'''
	Compute the p-value by numerical integration -- for validation
	'''
 	return integrate.quad( pdf, phi, 30., args=(a))
	mpmath.mp.dps = 100

	f = lambda i: pdf(i, a)
# 	ans, err = integrate.quad( F_int_func, 1, 1e3, (a, b, c, A, x))
	ans, err =  mpmath.mp.quad(f, [phi, mpmath.mp.inf], error=True, verbose=verbose) #mpmath.mp.inf

	print ans, err
	return (ans, err)

def fisher_exact(a):
	if len(a[0])<2:
		return (-1, -1)
	if int(a[0][1]) == 0 or int(a[1][1]) == 0:
		return (-1, -1)
	N = normalize_table(a)
# 	odds_ratio = (float(a[1][1])/float(a[1][0]))/(float(a[0][1])/float(a[0][0]))
	fisher_res = stats.fisher_exact( N )
	return fisher_res

def test_plot_with_toys(a, filename='logOR_test'):
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab

	N = normalize_table(a)
	
	# Manipulate the table
	N_null = [[N[0][0], N[0][1]], [N[1][0], float(N[1][0])*float(N[0][1])/float(N[0][0])]]
	N_flat = [item for sublist in N for item in sublist]
	N_flat_null = [item for sublist in N_null for item in sublist]
	
	C = [1.,-1.,-1.,1.]


	OR =  float(N[0][0]) * float(N[1][1]) /  float(N[1][0] * N[0][1])
	lOR = math.log(OR)

	print 'Null model:', N_null
	print "OR:", OR
	print "log(OR):", lOR
	print 
	a_prior = 0.5
	C_prior = [0.5,0.5,0.5,0.5]
	C_prior_null = [0.5,0.5,0.5,0.]

	tot = 30000.

	a_ = [n+a_prior for n in N_flat]
	a_null = [n+c for n, c in zip(N_flat_null, C_prior_null)]
	
	# Compute the approximations of the first two moments (or mean and sigma at least)
	mu    = np.sum(  [  c*(math.log(j-0.5+1./(24.*j))) for c, j in zip(C, a_) ]  )
	sigma = math.sqrt(	np.sum(  [  1./(j-0.5+1./(24.*j)) for c, j in zip(C, a_) ]  )	)
	print "Bloch/Watson  {:.4e}, sigma {:.4e}\n".format(mu, sigma)


	mu    = np.sum(  [  c*(math.log(j)) for c, j in zip(C, a_) ]  )
	sigma = math.sqrt(	np.sum(  [  1./(j) for c, j in zip(C, a_) ]  )	)
	print "Naive mu  {:.4e}, sigma {:.4e}\n".format(mu, sigma)


	mu    = np.sum(  [  c*(math.log(j) - 1./(2.*j)) for c, j in zip(C, a_) ]  )
	sigma = math.sqrt(	np.sum(  [  1./(j) for c, j in zip(C, a_) ]  )	)
	print "Normal mu  {:.4e}, sigma {:.4e}\n".format(mu, sigma)

	mu = b_l(1,a_)
	sigma = math.sqrt(b_l(2,a_))
	print "Exact mu (b_l)  {:.4e}, sigma {:.4e}\n".format(mu, sigma )
	
	# This is an estimate of the inaccuracy of using a normal distribution 	
	print "O(min_a) = ", 1./np.min(N)**3

	print
	print "Posterior Prob. of Null Hypothesis:", pdf( 0, a_)
	print "Posterior P-val of Null Hypothesis:", cdf( 0, a_)

	print

	size = 100000

	fig = plt.figure()
	ax = fig.add_subplot(111)

	# Generate toys to compare the approximate distributions:
	S = np.random.dirichlet(  a_, size=size)
	A = [tot*x[0] for x in S]
	B = [tot*x[3] for x in S]
	C = [tot*x[2] for x in S]
	D = [tot*x[1] for x in S]

	import collections
	dict_ = collections.Counter()

	x = []
	for i in range(len(A)):
			val = 0
			a = A[i]
			b = B[i]
			c = C[i] 
			d = D[i]
			if  B[i] == 0:
				val = 0
			if  D[i] == 0:
				val = 0
			if (A[i] * b * C[i] * d > 0):	
				val = math.log(  float(A[i])*float(B[i]) / (float(C[i])*float(D[i]))  )   #/    (sigma*math.sqrt(tot))
	# 		if (val < 0.005):
	# 			print a, b, c, d, val
			dict_[b]+=1
			x.append(val)

	n, bins, patches = ax.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

	print "Kurtosis:", stats.kurtosis(x)


	# Same for the null model
	S = np.random.dirichlet( a_null[:3], size=size )
	A = [tot*x[0] for x in S]
	B = [a_null[3] for x in S]
	C = [tot*x[2] for x in S]
	D = [tot*x[1] for x in S]

	import collections
	dict_ = collections.Counter()

	y = []
	for i in range(len(A)):
			val = 0
			a = A[i]
			b = B[i]
			c = C[i] 
			d = D[i]
			if  B[i] == 0:
				val = 0
			if  D[i] == 0:
				val = 0
			if (A[i] * b * C[i] * d > 0):	
				val = math.log(  float(A[i])*float(B[i]) / (float(C[i])*float(D[i]))  )   #/    (sigma*math.sqrt(tot))
	# 		if (val < 0.005):
	# 			print a, b, c, d, val
			dict_[b]+=1
			y.append(val)

 	n_null, bins_null, patches_null = ax.hist(y, 50, normed=1, facecolor='orange', alpha=0.75)

	bincenters = 0.5 * (bins[1:]+bins[:-1])
 	bincenters = 0.5 * (bins_null[1:]+bins_null[:-1])
	
	
	print a_null
	print a_
	

	print "Posterior model:"

	print "pdf"	
	q = [pdf( phi, a_) for phi in [bin for bin in bincenters]]
	l = ax.plot(bincenters, q, 'r--', linewidth=1)

	print "cdf"
# 	w = [cdf( phi, a_, verbose=True) for phi in [bin for bin in bincenters]]
# 	l = ax.plot(bincenters, w, 'k.', linewidth=1)


# 	z = stats.norm.cdf( bincenters, loc=mu, scale=sigma) 
# 	l = ax.plot(bincenters, z, 'b--', linewidth=1)
# 
# 	z = stats.norm.pdf( bincenters, loc=mu, scale=sigma) 
# 	l = ax.plot(bincenters, z, 'y.', linewidth=1)
# 

	print "Null model:"

# Plot the distributions
	print "pdf"
	u = [pdf( phi, a_null) for phi in [bin for bin in bincenters]]
	l = ax.plot(bincenters, u, 'r', linewidth=1)
# # 
# 	print "cdf"
# 	v = [cdf( phi, a_null, verbose=True) for phi in [bin for bin in bincenters]]
# 	l = ax.plot(bincenters, v, 'c.', linewidth=1)


	C_null = [1.,1.,1.,0.]

	mu    = 0# np.sum(  [  c*(math.log(j) - 1./(2.*j)) for c, j in zip(C_null, a_null) ]  )
	sigma = math.sqrt(	np.sum(  [  c/(j) for c, j in zip(C_null, a_null) ]  )	)

# 	mu = b_l(1, a_null)
# 	sigma = math.sqrt(b_l(2, a_null))

	print "a_null", a_null
	print C_null
	print "Null mu: {:.2f}, sigma: {:.2f}".format(mu, sigma)

	z = stats.norm.cdf( bincenters, loc=0., scale=sigma) 
	l = ax.plot(bincenters, z, 'b', linewidth=1)

	z = stats.norm.pdf( bincenters, loc=0., scale=sigma) 
	l = ax.plot(bincenters, z, 'k', linewidth=1)
	
	
	# Mark the lOR and the mean of the normal	
	ax.axvline(x=b_l(1,a_), ymin=0, ymax=100, color='r')
	ax.axvline(x=lOR, ymin=0, ymax=100, color='b')
	
	# Compute the fisher also for fun
	fisher_res = stats.fisher_exact(N)
	bayes_res = pval_null(N)
	norm_res = pval_norm(N)
# 	logit_res = gm.logistic_regression(df_pheno, df_snp, loud=False)	

	print "Fisher p-value:{}".format(fisher_res)
	print "Bayesian p-value:{}".format(bayes_res)
	print "Norm p-value:{}".format(norm_res)

	fig.text(0.15, 1.-0.15, r'Cont. table: {}'.format(a_), fontsize=15)
	fig.text(0.15, 1.-0.2, r'Null Cont. table: {}'.format(a_null), fontsize=15)
	fig.text(0.15, 1.-0.25, r'OR: {}'.format(OR), fontsize=15)
	fig.text(0.15, 1.-0.30, r'p-val: {:.2e}'.format(bayes_res[1]), fontsize=15)

	plt.savefig("{}.png".format(filename))
	plt.savefig("{}.pdf".format(filename))
	
	


