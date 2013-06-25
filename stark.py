# -*- coding: iso-8859-1 -*-
# Copyright (C) 2013 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

class stark(object):
	def __compute_tau_xi(self):
		from mpmath import sqrt, asin, ellipf, mpc, atan
		if self.bound:
			# Gradshtein 3.131.3.
			u = self.__init_coordinates[0]**2 / 2
			c,b,a = self.__roots_xi
			gamma = asin(sqrt((u - c)/(b - c)))
			q = sqrt((b - c)/(a - c))
			# NOTE: here it's q**2 instead of q because of the difference in
			# convention between G and mpmath :(
			# NOTE: the external factor 1 / sqrt(8 * self.eps) comes from the fact
			# that G gives the formula for a polynomial normalised by its leading coefficient.
			retval = 1 / sqrt(8 * self.eps) * (2 / sqrt(a - c)) * ellipf(gamma,q**2)
		else:
			# Here we will need two cases: one for when the other two roots are real
			# but negative, one for when the other two roots are imaginary.
			if isinstance(self.__roots_xi[1],mpc):
				# G 3.138.7.
				m = self.__roots_xi[1].real
				n = abs(self.__roots_xi[1].imag)
				# Only real root is the first one.
				a = self.__roots_xi[0]
				u = self.__init_coordinates[0]**2 / 2
				p = sqrt((m - a)**2 + n**2)
				retval = 1 / sqrt(8 * self.eps) * (1 / sqrt(p)) * ellipf(2 * atan(sqrt((u - a) / p)),(p + m - a) / (2*p))
			else:
				# G 3.131.7.
				u = self.__init_coordinates[0]**2 / 2
				c,b,a = self.__roots_xi
				mu = asin(sqrt((u - a)/(u - b)))
				q = sqrt((b - c)/(a - c))
				retval = 1 / sqrt(8 * self.eps) * (2 / sqrt(a - c)) * ellipf(mu,q**2)
		# Fix the sign according to the initial momentum.
		if self.__init_momenta[0] >= 0:
			return -abs(retval)
		else:
			return abs(retval)
	def __compute_tau_eta(self):
		# Gradshtein 3.131.5.
		from mpmath import sqrt, asin, ellipf
		u = self.__init_coordinates[1]**2 / 2
		c,b,a = self.__roots_eta
		k = asin(sqrt(((a - c) * (u - b)) / ((a - b) * (u - c))))
		p = sqrt((a - b) / (a - c))
		retval = 1 / sqrt(8 * self.eps) * (2 / sqrt(a - c)) * ellipf(k,p**2)
		if self.__init_momenta[1] >= 0:
			return -abs(retval)
		else:
			return abs(retval)
	def __init__(self,eps,x0,v0):
		from numpy import dot
		from mpmath import polyroots, mpf, mpc, sqrt, atan2, polyval
		from weierstrass_elliptic import weierstrass_elliptic as we
		if eps <= 0:
			raise ValueError('thrust must be strictly positive')
		eps = mpf(eps)
		# Unitary grav. parameter.
		mu = mpf(1.)
		x,y,z = [mpf(x) for x in x0]
		vx,vy,vz = [mpf(v) for v in v0]
		r = sqrt(x**2 + y**2 + z**2)
		xi = sqrt(r+z)
		eta = sqrt(r-z)
		phi = atan2(y,x)
		vr = dot(v0,x0) / r
		vxi = (vr + vz) / (2 * sqrt(r+z))
		veta = (vr - vz) / (2 * sqrt(r-z))
		vphi = (vy*x - vx*y) / (x**2 + y**2)
		pxi = (xi**2 + eta**2) * vxi
		peta = (xi**2 + eta**2) * veta
		pphi = xi**2*eta**2*vphi
		if pphi == 0:
			raise ValueError('bidimensional case')
		# Energy constant.
		h = (pxi**2 + peta**2) / (2*(xi**2 + eta**2)) + pphi**2 / (2*xi**2*eta**2) - (2 * mu) / (xi**2+eta**2) - eps * (xi**2 - eta**2) / 2
		# Alpha constants.
		alpha1 = -eps * xi**4 / 2 - h * xi**2 + pxi**2 / 2 + pphi**2 / (2*xi**2)
		alpha2 = eps * eta**4 / 2 - h * eta**2 + peta**2 / 2 + pphi**2 / (2*eta**2)
		# Analysis of the cubic polynomials in the equations for pxi and peta.
		roots_xi, _ = polyroots([8*eps,8*h,4*alpha1,-pphi**2],error=True,maxsteps=100)
		roots_eta, _ = polyroots([-8*eps,8*h,4*alpha2,-pphi**2],error=True,maxsteps=100)
		# NOTE: these are all paranoia checks that could go away if we used the exact cubic formula.
		if not (all([isinstance(x,mpf) for x in roots_xi]) or (isinstance(roots_xi[0],mpf) and isinstance(roots_xi[1],mpc) and isinstance(roots_xi[2],mpc))):
			raise ValueError('invalid xi roots detected: ' + str(roots_xi))
		if not (all([isinstance(x,mpf) for x in roots_eta]) or (isinstance(roots_eta[0],mpf) and isinstance(roots_eta[1],mpc) and isinstance(roots_eta[2],mpc))):
			raise ValueError('invalid eta roots detected: ' + str(roots_eta))
		# For xi we need to understand which of the real positive roots will be or was reached
		# given the initial condition.
		rp_roots_extract = lambda x: isinstance(x,mpf) and x > 0
		rp_roots_xi = [sqrt(2 * _) for _ in filter(rp_roots_extract,roots_xi)]
		rp_roots_eta = [sqrt(2. * _) for _ in filter(rp_roots_extract,roots_eta)]
		# Paranoia.
		if not len(rp_roots_xi) in [1,3]:
			raise ValueError('invalid xi roots detected: ' + str(roots_xi))
		if len(rp_roots_eta) != 2:
			raise ValueError('invalid eta roots detected: ' + str(roots_eta))
		# We choose as reachable/reached roots always those corresponding to the "pericentre"
		# for the two coordinates.
		if len(rp_roots_xi) == 1:
			# Here there's really no choice, only 1 root available.
			rr_xi = rp_roots_xi[0]
		else:
			# If motion is unbound, take the only root, otherwise take the smallest of the
			# two roots of the bound motion.
			rr_xi = rp_roots_xi[-1] if xi >= rp_roots_xi[-1] else rp_roots_xi[0]
		# No choice to be made here.
		rr_eta = rp_roots_eta[0]
		# Store parameters and constants.
		self.__init_coordinates = [xi,eta,phi]
		self.__init_momenta = [pxi,peta,pphi]
		self.__eps = eps
		self.__h = h
		self.__alpha1 = alpha1
		self.__alpha2 = alpha2
		self.__rp_roots_xi = rp_roots_xi
		self.__rp_roots_eta = rp_roots_eta
		self.__rr_xi = rr_xi
		self.__rr_eta = rr_eta
		self.__roots_xi = roots_xi
		self.__roots_eta = roots_eta
		# Create the Weierstrass objects for xi and eta.
		a1, a2, a3, a4 = 2*eps, (4 * h)/3, alpha1, -pphi**2
		g2 = -4 * a1 * a3 + 3 * a2**2
		g3 = 2 * a1 * a2 * a3 - a2**3 - a1**2*a4
		self.__f_xi = [4*a1,6*a2,4*a3,a4]
		self.__fp_xi = [12*a1,12*a2,4*a3]
		self.__fpp_xi = [24*a1,12*a2]
		self.__w_xi = we(g2,g3)
		# Eta.
		a1,a3 = -a1,alpha2
		g2 = -4 * a1 * a3 + 3 * a2**2
		g3 = 2 * a1 * a2 * a3 - a2**3 - a1**2*a4
		self.__f_eta = [4*a1,6*a2,4*a3,a4]
		self.__fp_eta = [12*a1,12*a2,4*a3]
		self.__fpp_eta = [24*a1,12*a2]
		self.__w_eta = we(g2,g3)
		# Compute the taus.
		tau_xi = self.__compute_tau_xi()
		tau_eta = self.__compute_tau_eta()
		self.__tau_xi = tau_xi
		self.__tau_eta = tau_eta
		# Store the real periods.
		self.__T_xi = self.__w_xi.periods[0]
		self.__T_eta = self.__w_eta.periods[0]
		# Delta bound (for debugging).
		xi_roots = self.__w_xi.roots
		# Determine the root corresponding to the real half-period.
		e_R = min(xi_roots,key = lambda x: abs(self.__w_xi.P(self.__T_xi/2) - x))
		self.__Dbound = e_R - polyval(self.__fpp_xi,xi**2/2) / 24
	def xi(self,tau):
		from mpmath import sqrt, polyval
		xi_r = self.__rr_xi
		retval = xi_r**2 + polyval(self.__fp_xi,xi_r**2/2) / (2 * (self.__w_xi.P(tau - self.tau_xi) - polyval(self.__fpp_xi,xi_r**2/2) / 24))
		return sqrt(retval)
	def eta(self,tau):
		from mpmath import sqrt, polyval
		eta_r = self.__rr_eta
		retval = eta_r**2 + polyval(self.__fp_eta,eta_r**2/2) / (2 * (self.__w_eta.P(tau - self.tau_eta) - polyval(self.__fpp_eta,eta_r**2/2) / 24))
		return sqrt(retval)
	@staticmethod
	def __ln_expansion(wp,u,n_iter = 25):
		from mpmath import ln, pi, sin, exp, mpc
		# TODO: test for pathological cases?
		om1, om3 = wp.periods[0]/2, wp.periods[1]/2
		assert(u.real < 2*om1 and u.real >= 0)
		assert(u.imag < 2*om3.imag and u.imag >= 0)
		tau = om3/om1
		q = exp(mpc(0,1)*pi()*tau)
		eta1 = wp.zeta(om1)
		om1_2 = om1 * 2
		retval = ln(om1_2/pi()) + eta1*u**2/(om1_2) + ln(sin(pi()*u/(om1_2)))
		for r in range(1,n_iter + 1):
			q_2 = q**(2*r)
			retval += q_2/(r * (1 - q_2))*(2. * sin(r*pi()*u/(om1_2)))**2
		return retval
	@staticmethod
	def __reduce_to_frp(z,om1):
		from mpmath import floor, mpc
		N = int(floor(z.real/(2*om1)))
		xs = z.real - N*2*om1
		return mpc(xs,z.imag),N
	@staticmethod
	def __ln_sigma(wp,u):
		from mpmath import mpc, pi
		# TODO: test for pathological cases?
		om1, om3 = wp.periods[0]/2, wp.periods[1]/2
		assert(u.imag >= 0 and u.imag < 2*om3.imag)
		further_reduction = False
		if u.imag/(2*om3.imag) > .5:
			further_reduction = True
			u = -u + 2*om1 + 2*om3
		u_F,N = stark.__reduce_to_frp(u,om1)
		eta1 = wp.zeta(om1)
		retval = stark.__ln_expansion(wp,u_F) + 2*N*eta1*(u_F + N*om1) - mpc(0,1)*N*pi()
		if further_reduction:
			retval -= (u - om1 - om3) * (2*eta1 + 2*wp.zeta(om3))
		return retval
	def phi(self,tau):
		from mpmath import polyval
		phi0 = self.__init_coordinates[2]
		pphi = self.__init_momenta[2]
		fp_xi, fp_eta = self.__fp_xi, self.__fp_eta
		fpp_xi, fpp_eta = self.__fpp_xi, self.__fpp_eta
		xi_r, eta_r = self.__rr_xi, self.__rr_eta
		w_xi, w_eta = self.__w_xi, self.__w_eta
		tau_xi, tau_eta = self.__tau_xi, self.__tau_eta
		beta_xi, beta_eta = -polyval(fpp_xi,xi_r**2/2) / 24, -polyval(fpp_eta,eta_r**2/2) / 24
		gamma_xi, gamma_eta = 2*xi_r**2, 2*eta_r**2
		delta_xi, delta_eta = polyval(fp_xi,xi_r**2/2) + 2*xi_r**2*beta_xi, polyval(fp_eta,eta_r**2/2) + 2*eta_r**2*beta_eta
		u_xi, u_eta = w_xi.Pinv(-delta_xi/gamma_xi), w_eta.Pinv(-delta_eta/gamma_eta)
		retval = phi0 + 2*pphi*(tau * (1 / gamma_xi + 1 / gamma_eta) + \
			(delta_xi - beta_xi*gamma_xi) / (gamma_xi**2 * w_xi.Pprime(u_xi)) * \
			(stark.__ln_sigma(w_xi,tau-tau_xi+u_xi) - stark.__ln_sigma(w_xi,-tau+tau_xi+u_xi) - stark.__ln_sigma(w_xi,-tau_xi+u_xi) + stark.__ln_sigma(w_xi,tau_xi+u_xi) - 2*tau*w_xi.zeta(u_xi)) + \
			(delta_eta - beta_eta*gamma_eta) / (gamma_eta**2 * w_eta.Pprime(u_eta)) * \
			(stark.__ln_sigma(w_eta,tau-tau_eta+u_eta) - stark.__ln_sigma(w_eta,-tau+tau_eta+u_eta) - stark.__ln_sigma(w_eta,-tau_eta+u_eta) + stark.__ln_sigma(w_eta,tau_eta+u_eta) - 2*tau*w_eta.zeta(u_eta)))
		return retval
	def t(self,tau):
		from mpmath import polyval
		xi_r, eta_r = self.__rr_xi, self.__rr_eta
		fp_xi, fp_eta = self.__fp_xi, self.__fp_eta
		fpp_xi, fpp_eta = self.__fpp_xi, self.__fpp_eta
		w_xi, w_eta = self.__w_xi, self.__w_eta
		tau_xi, tau_eta = self.__tau_xi, self.__tau_eta
		# These are the roots at the denominator.
		e_xi, e_eta = polyval(fpp_xi,xi_r**2/2)/24, polyval(fpp_eta,eta_r**2/2)/24
		# Link the roots to the half-periods.
		def e_to_omega(wp,e):
			return min([p/2 for p in list(wp.periods) + [sum(wp.periods)]],key = lambda p: abs(wp.P(p) - e))
		om_xi, om_eta = e_to_omega(w_xi,e_xi), e_to_omega(w_eta,e_eta)
		g2_xi, g2_eta = w_xi.invariants[0], w_eta.invariants[0]
		print(abs(w_xi.P(om_xi)-e_xi))
		print(abs(w_eta.P(om_eta)-e_eta))
		print(om_xi,om_eta)
		retval = tau * (xi_r**2 + eta_r**2) + \
			polyval(fp_xi,xi_r**2/2) / (2 * (g2_xi/4 - 3*e_xi**2)) * (tau*e_xi+w_xi.zeta(tau-tau_xi-om_xi)-w_xi.zeta(-tau_xi-om_xi)) + \
			polyval(fp_eta,eta_r**2/2) / (2 * (g2_eta/4 - 3*e_eta**2)) * (tau*e_eta+w_eta.zeta(tau-tau_eta-om_eta)-w_eta.zeta(-tau_eta-om_eta))
		return retval
	def state(self,tau):
		return [self.xi(tau),self.eta(tau),self.phi(tau)]
	def cart_state(self,tau):
		from mpmath import cos, sin
		xi,eta,phi = self.state(tau)
		return xi*eta*cos(phi),xi*eta*sin(phi),(xi**2 - eta**2)/2
	def __repr__(self):
		retval = ''
		retval += 'initial_coordinates: ' + repr(self.initial_coordinates) + '\n'
		retval += 'initial_momenta: ' + repr(self.initial_momenta) + '\n'
		retval += 'eps: ' + repr(self.eps) + '\n'
		retval += 'h: ' + repr(self.h) + '\n'
		retval += 'alpha1: ' + repr(self.alpha1) + '\n'
		retval += 'alpha2: ' + repr(self.alpha2) + '\n'
		retval += 'rp_roots_xi: ' + repr(self.rp_roots_xi) + '\n'
		retval += 'rp_roots_eta: ' + repr(self.rp_roots_eta) + '\n'
		retval += 'rr_xi: ' + repr(self.rr_xi) + '\n'
		retval += 'rr_eta: ' + repr(self.rr_eta) + '\n'
		retval += 'bound: ' + repr(self.bound) + '\n'
		retval += 'Dbound: ' + repr(self.__Dbound) + '\n'
		retval += 'tau_xi: ' + repr(self.tau_xi) + '\n';
		retval += 'tau_eta: ' + repr(self.tau_eta) + '\n';
		retval += 'T_xi: ' + repr(self.T_xi) + '\n';
		retval += 'T_eta: ' + repr(self.T_eta) + '\n';
		return retval
	@property
	def initial_coordinates(self):
		return self.__init_coordinates
	@property
	def initial_momenta(self):
		return self.__init_momenta
	@property
	def eps(self):
		return self.__eps
	@property
	def h(self):
		return self.__h
	@property
	def alpha1(self):
		return self.__alpha1
	@property
	def alpha2(self):
		return self.__alpha2
	@property
	def rp_roots_xi(self):
		return self.__rp_roots_xi
	@property
	def rp_roots_eta(self):
		return self.__rp_roots_eta
	@property
	def rr_xi(self):
		return self.__rr_xi
	@property
	def rr_eta(self):
		return self.__rr_eta
	@property
	def bound(self):
		return len(self.rp_roots_xi) == 3 and self.rr_xi != self.rp_roots_xi[-1]
	@property
	def tau_xi(self):
		return self.__tau_xi
	@property
	def tau_eta(self):
		return self.__tau_eta
	@property
	def T_xi(self):
		return self.__T_xi
	@property
	def T_eta(self):
		return self.__T_eta

def ln_expansion(wp,u,n_iter = 50):
	from mpmath import ln, pi, sin, exp, mpc
	# TODO: test for pathological cases?
	if wp.periods[1].real == 0:
		om1, om3 = wp.periods[0]/2, wp.periods[1]/2
	else:
		om1, om3 = wp.periods[1].conjugate()/2, wp.periods[1]/2
	tau = om3/om1
	q = exp(mpc(0,1)*pi()*tau)
	eta1 = wp.zeta(om1)
	om1_2 = om1 * 2
	retval = ln(om1_2/pi()) + eta1*u**2/(om1_2) + ln(sin(pi()*u/(om1_2)))
	for r in range(1,n_iter + 1):
		q_2 = q**(2*r)
		retval += q_2/(r * (1 - q_2))*(2. * sin(r*pi()*u/(om1_2)))**2
	return retval

def ln_expansion2(wp,u,n_iter = 25):
	from mpmath import ln, pi, sin, exp, mpc
	# TODO: test for pathological cases?
	om1, om3 = wp.periods[0]/2, wp.periods[1]/2
	assert(u.real < 2*om1 and u.real >= 0)
	assert(abs(u.imag) < 2*om3.imag)
	tau = om3/om1
	q = exp(mpc(0,1)*pi()*tau)
	eta1 = wp.zeta(om1)
	om1_2 = om1 * 2
	retval = ln(om1_2/pi()) + eta1*u**2/(om1_2) + ln(sin(pi()*u/(om1_2)))
	for r in range(1,n_iter + 1):
		q_2 = q**(2*r)
		retval += q_2/(r * (1 - q_2))*(2. * sin(r*pi()*u/(om1_2)))**2
	return retval
