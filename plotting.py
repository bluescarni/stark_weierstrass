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

def plot_poly(var,eps,x0,v0):
	from math import sqrt, atan2
	from pylab import plot, xticks, yticks, grid
	from numpy import dot, linspace
	x,y,z = x0
	vx,vy,vz = v0
	r = sqrt(x**2 + y**2 + z**2)
	xi = sqrt(r+z)
	eta = sqrt(r-z)
	phi = atan2(y,x)
	vr = dot(v0,x0) / r
	vxi = (vr + vz) / (2. * sqrt(r+z))
	veta = (vr - vz) / (2. * sqrt(r-z))
	vphi = (vy*x - vx*y) / (x**2 + y**2)
	pxi = (xi**2 + eta**2) * vxi
	peta = (xi**2 + eta**2) * veta
	pphi = xi**2*eta**2*vphi
	# Energy constant.
	E = 0.5 * (pxi**2 + peta**2) / (xi**2 + eta**2) + 0.5 * pphi**2 / (xi**2*eta**2) - 2. / (xi**2+eta**2) - eps * (xi**2 - eta**2) / 2.
	# Alpha constants.
	alpha1 = -eps*xi**4/2. - E*xi**2 + pxi**2/2. + 0.5*pphi**2/xi**2
	alpha2 = eps*eta**4/2. - E*eta**2 + peta**2/2. + 0.5*pphi**2/eta**2
	def func_xi(xi):
		try:
			return 1./abs(xi) * sqrt(eps*xi**6+2.*E*xi**4+2.*alpha1*xi**2-pphi**2)
		except ValueError:
			return float('nan')
	def func_eta(eta):
		try:
			return 1./abs(eta) * sqrt(-eps*eta**6+2.*E*eta**4+2.*alpha2*eta**2-pphi**2)
		except ValueError:
			return float('nan')
	if var == 'xi':
		func = func_xi
	else:
		func = func_eta
	rng = linspace(-10,10,10000)
	plot(rng,[func(xi) for xi in rng],'k-',linewidth=2)
	plot(rng,[-func(xi) for xi in rng],'k-',linewidth=2)
	xticks([0],[''])
	yticks([0],[''])
	grid()

def plot_all():
	from pylab import xlabel, ylabel, xlim, ylim, subplot, title, rc
	rc('text', usetex=True)
	subplot(141)
	title(r'(a)')
	plot_poly('xi',1,[.1,.2,.3],[.1,-10,-.5])
	xlabel(r'$\xi$')
	ylabel(r'$p_\xi$')
	subplot(142)
	title(r'(b)')
	plot_poly('xi',1,[.6,.2,.3],[.1,.5,-.5])
	xlim(-3,3)
	ylim(-12,12)
	xlabel(r'$\xi$')
	ylabel(r'$p_\xi$')
	subplot(143)
	title(r'(c)')
	plot_poly('xi',1,[.8,.2,.3],[.1,.5,-.5])
	xlim(-3,3)
	ylim(-12,12)
	xlabel(r'$\xi$')
	ylabel(r'$p_\xi$')
	subplot(144)
	title(r'(d)')
	plot_poly('eta',1,[70.,.2,.3],[.1,.5,-.5])
	xlim(-10,10)
	ylim(-200,200)
	xlabel(r'$\eta$')
	ylabel(r'$p_\eta$')

def plot_all_2d():
	from pylab import xlabel, ylabel, xlim, ylim, subplot, title, rc
	rc('text', usetex=True)
	subplot(151)
	title(r'(a)')
	plot_poly('xi',1,[.1,.2,.3],[.1,.2,1.])
	xlabel(r'$\xi$')
	ylabel(r'$p_\xi$')
	xlim(-4,4)
	ylim(-10,10)
	subplot(152)
	title(r'(b)')
	plot_poly('xi',1,[.1,.2,.3],[.1,.2,3.])
	xlabel(r'$\xi$')
	ylabel(r'$p_\xi$')
	xlim(-5,5)
	ylim(-25,25)
	subplot(153)
	title(r'(c)')
	plot_poly('xi',80,[.1,.2,.3],[.1,.2,.5])
	xlabel(r'$\xi$')
	ylabel(r'$p_\xi$')
	xlim(-6,6)
	ylim(-600,600)
	subplot(154)
	title(r'(d)')
	plot_poly('eta',4,[.1,.2,1],[.1,.2,5.5])
	xlabel(r'$\eta$')
	ylabel(r'$p_\eta$')
	xlim(-3,3)
	ylim(-6,6)
	subplot(155)
	title(r'(e)')
	plot_poly('eta',10,[.1,.2,.3],[.1,.2,.1])
	xlabel(r'$\eta$')
	ylabel(r'$p_\eta$')
	xlim(-.7,.7)
	ylim(-3,3)

class dynsys(object):
	# ODE function.
	@staticmethod
	def __func(y,_,E,eps,pphi):
		pxi, peta, xi, eta = y
		pxi_tau = 2.*eps*xi**3+2.*E*xi+pphi**2/xi**3
		peta_tau = -2.*eps*eta**3+2.*E*eta+pphi**2/eta**3
		return [pxi_tau,peta_tau,pxi,peta]
	def __init__(self,eps,x0,v0):
		from math import sqrt, atan2
		from numpy import dot
		x,y,z = x0
		vx,vy,vz = v0
		r = sqrt(x**2 + y**2 + z**2)
		xi = sqrt(r+z)
		eta = sqrt(r-z)
		phi = atan2(y,x)
		vr = dot(v0,x0) / r
		vxi = (vr + vz) / (2. * sqrt(r+z))
		veta = (vr - vz) / (2. * sqrt(r-z))
		vphi = (vy*x - vx*y) / (x**2 + y**2)
		pxi = (xi**2 + eta**2) * vxi
		peta = (xi**2 + eta**2) * veta
		pphi = xi**2*eta**2*vphi
		# Energy constant.
		E = 0.5 * (pxi**2 + peta**2) / (xi**2 + eta**2) + 0.5 * pphi**2 / (xi**2*eta**2) - 2. / (xi**2+eta**2) - eps * (xi**2 - eta**2) / 2.
		self.__E = E
		self.__eps = eps
		self.__pphi = pphi
		self.__y0 = [pxi,peta,xi,eta]
	@property
	def eps(self):
		return self.__eps
	@property
	def E(self):
		return self.__E
	@property
	def y0(self):
		return self.__y0
	@property
	def pphi(self):
		return self.__pphi
	def integrate(self,t_fin,steps):
		from numpy import linspace
		from scipy.integrate import odeint
		t = linspace(0.,t_fin,steps)
		return odeint(self.__func,self.__y0,t,args=(self.__E,self.__eps,self.__pphi))

class dynsys_cart(object):
	@staticmethod
	def __func(y,_,eps):
		from math import sqrt
		x,y,z,vx,vy,vz = y
		r = sqrt(x**2 + y**2 + z**2)
		return vx,vy,vz,-x/r**3,-y/r**3,-z/r**3+eps
	def __init__(self,eps,x0,v0):
		self.__eps = eps
		self.__y0 = [x0[0],x0[1],x0[2],v0[0],v0[1],v0[2]]
	def integrate(self,t_fin,steps):
		from numpy import linspace
		from scipy.integrate import odeint
		t = linspace(0.,t_fin,steps)
		return odeint(self.__func,self.__y0,t,args=(self.__eps,))

def plot_bound_vs_unbound():
	from pylab import xlabel, ylabel, xlim, ylim, subplot, title, rc, linspace, plot, xticks, yticks, grid
	from stark import stark
	rc('text', usetex=True)
	npoints = 1000
	s = stark(.07,[1.,1.,.3],[.1,.4,.6])
	r = linspace(-s.T_xi*2,s.T_xi*2,npoints)
	subplot(241)
	title(r'(a)')
	plot(r,[s.xi(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$\xi\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()
	subplot(242)
	title(r'(b)')
	plot(r,[s.eta(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$\eta\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()
	subplot(243)
	title(r'(c)')
	plot(r,[s.phi(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$\phi\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()
	subplot(244)
	title(r'(d)')
	plot(r,[s.t(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$t\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()
	s = stark(.1,[1.,1.,.3],[.1,.4,.6])
	r = linspace(-s.T_xi*2,s.T_xi*2,npoints)
	subplot(245)
	title(r'(e)')
	plot(r,[s.xi(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$\xi\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()
	subplot(246)
	title(r'(f)')
	plot(r,[s.eta(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$\eta\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()
	subplot(247)
	title(r'(g)')
	plot(r,[s.phi(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$\phi\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()
	subplot(248)
	title(r'(h)')
	plot(r,[s.t(tau).real for tau in r],'k-',linewidth=2)
	xlabel(r'$\tau$')
	ylabel(r'$t\left(\tau\right)$')
	xticks([0],[''])
	yticks([0],[''])
	xlim(float(r[0]),float(r[-1]))
	#ylim(-.1,4)
	grid()

def plot_3d():
	from numpy import array
	from stark import stark
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from pylab import xlabel, ylabel, xlim, ylim, subplot, title, rc, linspace, plot, xticks, yticks, grid
	rc('text', usetex=True)
	npoints = 1000
	fig = plt.figure()
	ax = fig.add_subplot(121, projection='3d')
	s = stark(.08,[1.,0.,.05],[0,1,.02])
	title(r'(a)')
	r = linspace(0,s.T_xi*2,npoints)
	states = array([s.cart_state(tau) for tau in r])
	ax.plot([_.real for _ in states[:,0]],[_.real for _ in states[:,1]],[_.real for _ in states[:,2]],'k-',linewidth=2)
	ax.set_zlim(0,0.3)
	ax = fig.add_subplot(122, projection='3d')
	s = stark(.1793,[1.,0.,.05],[0,1,.02])
	title(r'(b)')
	r = linspace(0,s.T_xi*2,npoints)
	states = array([s.cart_state(tau) for tau in r])
	ax.plot([_.real for _ in states[:,0]],[_.real for _ in states[:,1]],[_.real for _ in states[:,2]],'k-',linewidth=2)

def plot_quasi_periodic():
	from numpy import array
	from stark import stark
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from pylab import xlabel, ylabel, xlim, ylim, subplot, title, rc, linspace, plot, xticks, yticks, grid
	rc('text', usetex=True)
	npoints = 1000
	x0,v0,eps = [1.0, 0, 6.123233995736766e-17],[-0.2833431329044651, 0.8928727330044925, 0.34957175101847515],0.0612722970737
	s = stark(eps,x0,v0)
	T = 3.4589043219148699*5
	r = linspace(0,T*3,npoints)
	states = array([s.cart_state(_) for _ in r])
	fig = plt.figure()
	ax = fig.add_subplot(121, projection='3d')
	title(r'(a)')
	ax.plot([_.real for _ in states[:,0]],[_.real for _ in states[:,1]],[_.real for _ in states[:,2]],'k-',linewidth=2)
	ax = fig.add_subplot(122, projection='3d')
	title(r'(b)')
	ax.plot([_.real for _ in states[:,0]],[_.real for _ in states[:,1]],[_.real for _ in states[:,2]],'k-',linewidth=2)


def plot_periodic():
	from mpmath import mpf
	from numpy import array
	from stark import stark
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from pylab import xlabel, ylabel, xlim, ylim, subplot, title, rc, linspace, plot, xticks, yticks, grid
	rc('text', usetex=True)
	npoints = 2000
	x0,v0,eps =  [1.0, 0, 6.123233995736766e-17], [0.020495878242533995, 1.1180786341556543, 0.05717527804838409], 0.080969489037
	s = stark(eps,x0,v0)
	T =  43.335885307166
	r = linspace(0,T,npoints)
	states = array([s.cart_state(_) for _ in r])
	fig = plt.figure()
	ax = fig.add_subplot(121, projection='3d')
	title(r'(a)')
	ax.plot([_.real for _ in states[:,0]],[_.real for _ in states[:,1]],[_.real for _ in states[:,2]],'k-',linewidth=2)
	ax = fig.add_subplot(122, projection='3d')
	title(r'(b)')
	ax.plot([_.real for _ in states[:,0]],[_.real for _ in states[:,1]],[_.real for _ in states[:,2]],'k-',linewidth=2)
