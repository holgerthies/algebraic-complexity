import poly
import numpy as np

def random_polynomial(deg,max_val=1000):
	coeffs = np.random.randint(max_val, size=deg+1).tolist()
	return poly.Poly(coeffs)

def poly_map(f,p):
	coeffs = [f(x) for x in p.coeffs]
	return poly.Poly(coeffs)