import poly
# naive multiplication algorithm
# complexity O(n^2)
def naive_multiplication(p1,p2):
	coeffs={}
	for i,a in enumerate(p1.coeffs):
		for j,b in enumerate(p2.coeffs):
			coeffs[i+j] = a*b
	return poly.Poly(coeffs)
