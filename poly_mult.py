import poly
import numpy as np
# naive multiplication algorithm
# complexity O(n^2)
def naive(p1,p2):
	coeffs=[0 for _ in range(p1.degree+p2.degree+1)]
	for i,a in enumerate(p1.coeffs):
		for j,b in enumerate(p2.coeffs):
			coeffs[i+j] += a*b
	return poly.Poly(coeffs)

# multiplies 2 polynomials of degree <= deg 
# by using Lagrangian interpolation
def __tk_multiply_by_interpolation(deg,p1,p2):
	# number of evaluation points
	N = 2*deg-1
	# use the points 0,1,2,..,N-1 as evaluation points
	x = np.arange(N)
	# Van der Monde Matrix 
	V = np.vander(x,N)
	# evaluate the polynomials at N points and multiply the results
	y = np.array([p1(xi)*p2(xi) for xi in x])
	# solve system of linear equations to get coefficients
	# NOTE: This could be made faster by making use of the matrix structure 
	a = np.linalg.solve(V,y)
	# make polynomial from the calculated coefficients
	list_a = a.tolist()
	list_a.reverse()
	p = poly.Poly(list_a)
	return p

# splits the polynomial into smaller polynomials s.t.
# each of the polynomials has degree < k
# i.e. returns polynomials p_0,..p_m of deg < k
# s.t. p = p_0+X^k*p_1+...+X^(m*k)p_m
def __tk_split_polynomial(p, k):
	polys = []
	for i in range(0,p.degree, k):
		q = poly.Poly(p.coeffs[i:i+k])
		polys.append(q)
	return polys

# generic Toom-Cook multiplication 
# Special cases: 
# k=1 : naive multiplication
# k=2 : Karatsuba algorithm 
# complexity O(n^(lok_k (2k+1)))
def toomCook(k,p1,p2):
	# split into polynomials of degree k
	if p1.degree <= k and p2.degree <= k:
		return __tk_multiply_by_interpolation(k, p1,p2)
	p1_split = __tk_split_polynomial(p1,k)

