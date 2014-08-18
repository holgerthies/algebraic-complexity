""" 
A basic class for polynomials with generic coefficient type
"""
import types
class Poly:
	coeffs = []
	degree = 0
	# construct from coefficient list or a string in the form a_nX^n+...+a0
	def __init__(self, coeffs):
		if isinstance(coeffs, list):
			# init from list of coefficients
			self.coeffs = coeffs
		elif isinstance(coeffs, dict):
			# init from dict of coefficients
			self.coeffs = [0 for _ in range(max(coeffs.keys())+1)]
			for i in coeffs.keys(): self.coeffs[i] = coeffs[i]
		else:
			# init from string
			self.__read(coeffs)
		self.degree = len(self.coeffs)


	# Print polynomial in the form a_nX^n+...+a0
	def __repr__(self):
		if len(self.coeffs) == 0:
			return ""
		l = []
		for i in range(len(self.coeffs)):
			if self.coeffs[i] != 0:
				if i == 0:
					l.append(str(self.coeffs[i]))
				elif i==1:
					l.append(str(self.coeffs[i])+"X")
				else:
					l.append(str(self.coeffs[i])+"X^"+str(i))
		l.reverse()
		return ' + '.join(l)

	# scalar multiplication
	def __mul__(self, other):
		new_coeffs = [x*other for x in self.coeffs]
		return Poly(new_coeffs)

	def __rmul__(self, other):
		return self*other

	# add two polynomials
	def __add__(self,other):
		new_coeffs = [0 for i in range(max(len(self.coeffs), len(other.coeffs)))]
		for i,x in enumerate(self.coeffs):
			new_coeffs[i] = x
		for i,x in enumerate(other.coeffs):
			new_coeffs[i] += x
		return Poly(new_coeffs)

	# Read a polynomial in the form the __str__  method returns
	def __read(self, str):
		# remove white spaces
		str.replace(" ", "")
		# get strings of form a_kX^k
		coeff_list = str.split("+")
		coeff = {}
		for s in coeff_list:
			# split string on X
			s_parsed = s.split("X")
			# set coefficient dict
			if len(s_parsed) == 1:
				coeff[0] = int(s_parsed[0])
			elif s_parsed[1] == "":
				coeff[1] = int(s_parsed[0])
			else:
				coeff[int(s_parsed[1].replace("^", ""))] = int(s_parsed[0])
		# set coefficient list with coefficients from dict
		self.coeffs = [0 for _ in range(max(coeff.keys())+1)]
		for i in coeff.keys(): self.coeffs[i] = coeff[i]

	# evaluation using Horner Schema
	def __call__(self,x):
		a = self.coeffs
		# Horner scheme
		n = len(a)-1
		# b_n := a_n
		b = [a[n]] # Note b is in reversed order so b[i] = b[n-i]
		for i in range(1,n+1):
			b.append(a[n-i]+b[i-1]*x)
		return b[n]




