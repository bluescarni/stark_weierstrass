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

#NOTE: this will need sympy 0.7.2 in order to work out of the box.
def check():
	from numpy import array, dot
	from sympy import Integer, Rational, Symbol, sqrt, simplify, diff
	mu,z,eps,s = [Symbol(_) for _ in ['mu','z','eps','s']]
	r_vec = array([sqrt((z*mu/eps)**Rational(2,3)-z**2),Integer(0),z])
	x,y,_ = r_vec
	r = sqrt(x**2 + y**2 + z**2)
	print(r)
	v_vec = array([Integer(0),sqrt(mu*((z*mu/eps)**Rational(2,3)-z**2)/(r**3)),Integer(0)])
	print(v_vec[1])
	v_x, v_y, v_z = v_vec
	xi, eta = sqrt(r+z), sqrt(r-z)
	rp = dot(r_vec,v_vec) / r
	v_xi, v_eta, v_phi = (rp + v_z) / (2 * sqrt(r + z)), (rp - v_z) / (2 * sqrt(r - z)), (v_y*x - v_x*y)/(x**2 + y**2)
	p_xi, p_eta, p_phi = (xi**2 + eta**2)*v_xi, (xi**2 + eta**2)*v_eta, xi**2*eta**2*v_phi
	h = (p_xi**2 + p_eta**2) / (2*(xi**2 + eta**2)) + p_phi**2/(2*xi**2*eta**2) - 2*mu / (xi**2 + eta**2) - eps * (xi**2 - eta**2) / 2
	alpha_1 = -eps * xi**4 / 2 - h*xi**2 + p_xi**2/2 + p_phi**2/(2*xi**2)
	alpha_2 =  eps * eta**4 / 2 - h*eta**2 + p_eta**2/2 + p_phi**2/(2*eta**2)
	f_xi = 8*eps*s**3 + 8*h*s**2 + 4*alpha_1*s - p_phi**2
	f_xi_p = diff(f_xi,'s')
	f_eta = -8*eps*s**3 + 8*h*s**2 + 4*alpha_2*s - p_phi**2
	f_eta_p = diff(f_eta,'s')
	retval = simplify(f_xi.subs('s',(xi**2)/2)) == 0 and simplify(f_eta.subs('s',(eta**2)/2)) == 0 and \
		simplify(f_xi_p.subs('s',(xi**2)/2)) == 0 and simplify(f_eta_p.subs('s',(eta**2)/2)) == 0
	return retval
