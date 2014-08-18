import poly
import numpy as np
import cmath
# naive multiplication algorithm
# complexity O(n^2)
def naive(p1,p2):
	coeffs=[0 for _ in range(p1.degree+p2.degree+1)]
	for i,a in enumerate(p1.coeffs):
		for j,b in enumerate(p2.coeffs):
			coeffs[i+j] += a*b
	return poly.Poly(coeffs)

# get the inverse vandermonde matrix
# from precomputed values
def __tk_inverse_vandermonde(N):
	if N==3:
		return np.matrix([[ 0.5, -1. ,  0.5],[-1.5,  2. , -0.5],[ 1. ,  0. ,  0. ]])
	# use the points 0,1,2,..,N-1 as evaluation points
	x = __tk_get_evaluationpoints(N)
	return np.matrix(np.vander(x,N)).I
def __tk_get_evaluationpoints(N):
	if N==3:
		return np.array([0,1,2])
	return np.arange(N)
# gets a list of coefficients representing a polynomial
# p_0 + p_1*Y + ... + p_n*Y^n where p_i are polynomials over X
# and substitutes X^degree_subst for Y
def __tk_restore_polynomial(l, degree_subst, max_degree):
	coeff_list = [0 for _ in range(max_degree+1)]
	for i,p in enumerate(l):
		for j,a in enumerate(p.coeffs):
			if j+degree_subst*i <= max_degree:
				coeff_list[j+degree_subst*i] += a
	return poly.Poly(coeff_list)
# multiplies 2 polynomials of 
# polynomials by using Lagrangian interpolation
# The multiplication is done by recursively calling 
# The toomCook algorithm of degree k
def __tk_multiply_by_interpolation(p1,p2, k):
	deg = p1.degree+p2.degree
	# number of evaluation points
	N = deg-1
	# evaluation points
	x = __tk_get_evaluationpoints(N)
	# Inverse van der Monde Matrix 
	V = __tk_inverse_vandermonde(N)
	# evaluate the polynomials at N points and multiply the results
	# compute the product polynomial by recursively calling toom-cook
	y = np.array([toomCook(k,p1(xi),p2(xi)) for xi in x])
	# solve system of linear equations to get coefficients
	# NOTE: This could be made faster by making use of the matrix structure 
	a = np.dot(V,y)
	# make polynomial from the calculated coefficients
	list_a = a.tolist()
	list_a = list_a[0]
	list_a.reverse()
	if not isinstance(p1.coeffs[0],poly.Poly):
		#list_a = [round(x) for x in list_a]
		return poly.Poly(list_a)
	# The substitution degree is one more than the highest coefficient degree in the first polynomial
	degree_subst = len(p1.coeffs[0].coeffs)
	p = __tk_restore_polynomial(list_a, degree_subst, deg*degree_subst)
	return p

# splits the polynomial into k smaller polynomials
# s.t. p = p_0+X^n/k*p_1+...+X^(m*n/k)p_m
def __tk_split_polynomial(p, k):
	if p.degree <= k:
		return p
	polys = []
	# number of coefficients in one polynomial
	n = p.degree/k
	for i in range(0,p.degree, n):
		q = poly.Poly(p.coeffs[i:i+n])
		polys.append(q)
	if len(polys) == 1:
		return p
	return poly.Poly(polys)

# generic Toom-Cook multiplication 
# Special cases: 
# k=2 : Karatsuba algorithm 
# So far only works when the degree is k^m
# complexity O(n^(lok_k (2k+1)))
def toomCook(k,p1,p2):
	if (not isinstance(p1,poly.Poly) and not isinstance(p2, poly.Poly)) or (p1.degree <= 1 and p2.degree <= 1):
		return p1*p2
	# split into k polynomials

	p1_split = __tk_split_polynomial(p1,k)
	p2_split = __tk_split_polynomial(p2,k)
	#print "________________________________________"
	#print p1,p2
	#print __tk_multiply_by_interpolation(p1_split, p2_split, k)
	#print "________________________________________"
	return __tk_multiply_by_interpolation(p1_split, p2_split, k)

# The discrete fourier transform algorithm
# input is a list A=[a_0,..a_n] of some field and an N-th root of unity of this field
# applies the fast fourier transform algorithm to this vector
def __fft_dft(A, omega):
	N = len(A)
	if N==1:
		return A
	# split into even and odd list entries
	A_even = A[::2]
	A_odd = A[1::2]
	#omega^2 is an N/2-th root of unity
	omega2 = omega*omega
	# recursively apply fft on even and odd list elements
	DFT_even = __fft_dft(A_even, omega2) 
	DFT_odd = __fft_dft(A_odd, omega2)
	# make list for result
	A_transformed = [0 for _ in range(N)]
	# omega_pow will contain omega^i in the i-th call of the loop
	omega_pow = 1
	for i in range(N/2):
		# The first half of DFT(N)*A is given by DFT(N/2)*A_even+diag(1,omega,..., omega^(N/2-1))*A_odd
		A_transformed[i] = DFT_even[i]+omega_pow*DFT_odd[i]
		# The second half of DFT(N)*A is given by DFT(N/2)*A_even-diag(1,omega,..., omega^(N/2-1))*A_odd
		A_transformed[N/2+i] = DFT_even[i]-omega_pow*DFT_odd[i]
		# omega value for next iteration
		omega_pow *= omega
	return A_transformed

# computes the inverse of the discrete fourier transform
# We want to compute A=vander(omega,N)^-1*A_trans = 1/N*vander(omega^-1,N)*A_trans
# Thus DFT^-1(A,omega) = DFT(A,omega^-1)/N
def __fft_inverse_dft(A,omega):
	N = len(A)
	# inverse of omega
	omega_inverse = 1/omega
	# Apply DFT
	return [int(round((x/N).real)) for x in __fft_dft(A,omega_inverse)]

# multiplication using fast fourier transform
def fft(p1,p2):
	# assume both polynomials have the same degree which is a power of 2
	#degree of the resulting polynomial
	N = len(p1.coeffs)+len(p2.coeffs)
	# Fill coefficient list for fast fourier transform so that both vectors have length N
	A1 = p1.coeffs+[0 for _ in range(N-len(p1.coeffs))]
	A2 = p2.coeffs+[0 for _ in range(N-len(p2.coeffs))]
	# assume complex polynomials thus root of unity is e^((2pi/N)i)
	omega = np.exp(-2j * np.pi  / N)
	# apply discrete fourier transformation to both polynomials
	A1_transformed = __fft_dft(A1, omega)
	A2_transformed =  __fft_dft(A2, omega)
	# combine the result by multipilying them to get the result of the FFT of p1*p2
	A_transformed = [A1_transformed[i]*A2_transformed[i] for i in range(N)]
	# apply inverse fourier transform to obtain the coefficients
	A = __fft_inverse_dft(A_transformed, omega)
	# return polynomial with coefficients from A
	return poly.Poly(A)



