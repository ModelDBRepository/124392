"""
 * Copyright (C) 2006 Evan Thomas
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
"""

"""
/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

This file implements:
   +----------------------------------------------------------------------+
   | Vijayalakshmi Santhakumar, Ildiko Aradi and Ivan Soltesz             |
   | Role of Mossy Fiber Sprouting and Mossy Cell Loss in                 |
   | Hyperexcitability:A Network Model of the Dentate Gyrus Incorporating |
   | Cell Types and Axonal Topography                                     |
   | J Neurophysiol 93: 437-453, 2005.                                    |
   +----------------------------------------------------------------------+

Mostly copied from the source on ModelDB

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
"""

import sys
sys.path = ['C:\Documents and Settings\evan\My Documents\Visual Studio Projects\parplex']+sys.path

from p3Build import VGCDynamicsBuilder, CaDynamicsBuilder, Gate, build

######################
# Regular Na current #
######################
h = Gate()
alpha_h = '0.23/exp((V+60+5)/20)'
beta_h  = '3.33/(exp(-(V+60-47.5)/10.0)+1)'
h.name = 'h'
h.exponent = 1
h.alpha = alpha_h
h.beta  = beta_h

m = Gate()
alpha_m = '-0.3*(V+60-17)/(exp((V+60-17)/-5) - 1)'
beta_m  = '0.3*(V+60-45)/(exp((V+60-45)/5) - 1)'
m.name = 'm'
m.exponent = 3
m.alpha = alpha_m
m.beta  = beta_m

Na = VGCDynamicsBuilder()
Na.gates = [h, m]
Na.description = 'Na: Na model'
Na.name = 'Na'

#########################################
# Na current with slow inactivation and #
# allowance for 'mutation' of gating    #
#########################################
alpha_h = '0.23/exp((V-Vhalfh+60+5)/20)'
beta_h  = '3.33/(exp(-(V-Vhalfh+60-47.5)/10.0)+1)'
h = Gate()
h.name = 'h'
h.exponent = 1
h.inf = '(%s)/((%s)+(%s))' % (alpha_h, alpha_h, beta_h)
h.tau = 'Ah/((%s)+(%s))' % (alpha_h, beta_h)

alpha_m = '-0.3*(V-Vhalfm+60-17)/(exp((V-Vhalfm+60-17)/-5) - 1)'
beta_m  = '0.3*(V-Vhalfm+60-45)/(exp((V-Vhalfm+60-45)/5) - 1)'
m = Gate()
m.name = 'm'
m.exponent = 3
m.inf = '(%s)/((%s)+(%s))' % (alpha_m, alpha_m, beta_m)
m.tau = 'Am/((%s)+(%s))' % (alpha_m, beta_m)

s = Gate()
s.name = 's'
s.exponent = 1
s.inf = '(%s)/((%s)+(%s))' % (alpha_h, alpha_h, beta_h)
s.tau = 'As*cAs*exp(-cBs*V)'

Naslow = VGCDynamicsBuilder()
Naslow.gates = [h, m, s]
Naslow.variable('Vhalfm', 0, ' activation voltage shift')
Naslow.variable('Vhalfh', 0, 'fast inactivation voltage shift')
Naslow.variable('Vhalfs', 0, 'slow inactivation voltage shift')
Naslow.variable('Am', 1, 'activation rate change')
Naslow.variable('Ah', 1, 'fast inactivation rate change')
Naslow.variable('As', 1, 'slow inactivation rate change')
Naslow.constant('cAs', 550)
Naslow.constant('cBs', 0.046)
Naslow.description = 'Naslow: Na model with slow kinetics'
Naslow.name = 'Naslow'

####################################
# Fast delayed rectifier K current #
####################################
n = Gate()
alpha_n = '-0.07*(V-Vhalfn+65-47)/(exp((V-Vhalfn+65-47)/-6)-1)'
beta_n  = '0.264/exp((V-Vhalfn+65-22.0)/40.0)'
n.name = 'n'
n.exponent = 4
n.alpha = alpha_n
n.beta  = beta_n

Kdrf = VGCDynamicsBuilder()
Kdrf.gates = [n]
Kdrf.variable('Vhalfn', 0, ' activation voltage shift')
Kdrf.description = 'K: fast delayed rectifier'
Kdrf.name = 'Kdrf'

##########################
# Slow delayed rectifier #
##########################
n = Gate()
n.alpha = '-0.028*(V-Vhalfn+65-35)/(exp((V-Vhalfn+65-35)/-6)-1)'
n.beta  = '0.1056/exp((V-Vhalfn+65-10)/40)'
n.name  = 'n'
n.exponent = 4

Kdrs = VGCDynamicsBuilder()
Kdrs.gates = [n]
Kdrs.variable('Vhalfn', 0, ' activation voltage shift')
Kdrs.description = 'K: slow delayed rectifier'
Kdrs.name = 'Kdrs'

#############
# A current #
#############
k = Gate()
alpha = 'exp(1e-3*zetan*(V-vhalfn)*9.648e4/(8.315*(273.16+celsius)))'
beta  = 'exp(1e-3*zetan*gmn*(V-vhalfn)*9.648e4/(8.315*(273.16+celsius)))'
k.inf = '1/(' + alpha + '+1)'
k.tau = beta + '/(q10*a0n*(1+' + alpha + '))'
k.name  = 'k'
k.exponent = 1
k.description = 'activation'
l = Gate()
alpha = 'exp(1e-3*zetal*(V-vhalfl)*9.648e4/(8.315*(273.16+celsius)))'
beta  = 'exp(1e-3*zetal*gml*(V-vhalfl)*9.648e4/(8.315*(273.16+celsius)))'
l.inf = '1/(' + alpha + '+1)'
l.tau = beta + '/(q10*a0l*(1 + ' + alpha + '))'
l.name  = 'l'
l.exponent = 1
l.description = 'inactivation'

Ka = VGCDynamicsBuilder()
Ka.gates = [k, l]
Ka.description = 'Ka: A current'
Ka.name = 'Ka'
Ka.variable('vhalfn', -33.6, 'Vhalf of activation')
Ka.variable('vhalfl', -83, 'Vhalf of inactivation')
Ka.constant('a0l', 0.08)
Ka.constant('a0n', 0.02)
Ka.constant('zetan', -3)
Ka.constant('zetal', 4)
Ka.constant('gmn', 0.6)
Ka.constant('gml', 1)
Ka.constant('celsius', 6.3)
Ka.constant('q10', 0.0740) # This is different for different channels! */

##################
# fast Ih system #
##################
h = Gate()
h.inf = '1 / (1 + exp( (V-Vhalfhf+91)/10 ))'
h.tau = '14.9 + 14.1 / (1+exp(-(V+95.2)/0.5))'
h.exponent = 2
h.name = 'h'
h.description = 'activation'

Ihf = VGCDynamicsBuilder()
Ihf.gates = [h]
Ihf.description = 'fast Ih model'
Ihf.name = 'Ihf'
Ihf.variable('Vhalfhf', 0, ' activation voltage shift')

##################
# slow Ih system #
##################
h = Gate()
h.inf = '1 / (1 + exp( (V-Vhalfhs+91)/10 ))'
h.tau = '80 + 172.7 / (1+exp(-(V+59.3)/-0.83))'
h.exponent = 2
h.name = 'h'
h.description = 'activation'

Ihs = VGCDynamicsBuilder()
Ihs.gates = [h]
Ihs.description = 'slow Ih model'
Ihs.name = 'Ihs'
Ihs.variable('Vhalfhs', 0, ' activation voltage shift')

#######
# CaT #
#######
a = Gate()
a.alpha = '(0.2*(19.26-(V-Vhalfmt))/(exp((19.26-(V-Vhalfmt))/10)-1))'
a.beta  = '(0.009*exp(-(V-Vhalfmt)/22.03))'
a.exponent = 2
a.description = 'activation'
a.name = 'a'
b = Gate()
b.alpha = '(1e-6*exp(-V/16.26))'
b.beta  = '(1.0/(exp((29.79-V)/10.0)+1.0))'
b.exponent = 1
b.description = 'inactivation'
b.name = 'b'

CaT = VGCDynamicsBuilder()
CaT.description = 'T type calcium channels'
CaT.gates = [a, b]
CaT.variable('Vhalfmt', 0, ' activation voltage shift')
CaT.name = 'CaT'

#######
# CaN #
#######
c = Gate()
c.alpha = '(0.19 * (19.88-(V-Vhalfmn)) / (exp((19.88-(V-Vhalfmn))/10) - 1))'
c.beta  = '(0.046*exp(-(V-Vhalfmn)/20.73))'
c.exponent = 2
c.description = 'activation'
c.name = 'c'
d = Gate()
d.alpha = '(1.6e-4/exp(-V/48.4))'
d.beta  = '(1.0/(exp((39.0-V)/10.0)+1.0))'
d.exponent = 1
d.description = 'inactivation'
d.name = 'd'

CaN = VGCDynamicsBuilder()
CaN.description = 'N type calcium channels'
CaN.gates = [c, d]
CaN.variable('Vhalfmn', 0, ' activation voltage shift')
CaN.name = 'CaN'

#######
# CaL #
#######
e = Gate()
e.alpha = '(15.69*(81.5-(V-Vhalfml))/(exp((81.5-(V-Vhalfml))/10)-1.0))'
e.beta  = '(0.29*exp(-(V-Vhalfml)/10.86))'
e.exponent = 2
e.description = 'activation'
e.name = 'e'

CaL = VGCDynamicsBuilder()
CaL.description = 'L type calcium channels'
CaL.gates = [e]
CaL.variable('Vhalfml', 0, ' activation voltage shift')
CaL.name = 'CaL'

#######
# Kbk #
#######
o = Gate()
o.alpha = 'Ca*abar/(Ca + (k1*exp(-2*d1*FARADAY*V/R/(273.15 + celsius))))'
o.beta  = 'bbar/(1 + Ca/(k2*exp(-2*d2*FARADAY*V/R/(273.15 + celsius))))'
o.exponent = 1
o.description = 'activation'
o.name = 'o'

Kbk = VGCDynamicsBuilder()
Kbk.description = 'Ca and voltage activated BK channels'
Kbk.gates = [o]
Kbk.name = 'Kbk'
Kbk.constant('d1', 0.84)
Kbk.constant('d2', 1.0)
Kbk.constant('k1', 0.48e-3)
Kbk.constant('k2', 0.13e-6)
Kbk.constant('abar', 0.28)
Kbk.constant('bbar', 0.48)
Kbk.constant('FARADAY', 96.4853)
Kbk.constant('R', 8.313424)
Kbk.constant('celsius', 6.3)

#######
# Ksk #
#######
q = Gate()
q.alpha = '(1.25e1 * Ca * Ca)'
q.beta  = '(0.00025)'
q.exponent = 2
q.description = 'activation'
q.name = 'q'

Ksk = VGCDynamicsBuilder()
Ksk.description = 'Ca activated SK channels'
Ksk.gates = [q]
Ksk.name = 'Ksk'


build('aradi', 'The Santhakumar et al dentate gyrus model')
