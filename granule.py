"""
/*
 * Copyright (C) 2007 Evan Thomas
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
 */
"""

#
# Implements the granule cell of Santhakumar et al 2005
#
from p3 import *
from math import pi, log
old = False
if old:
  from aradi_old import *
  method     = rk32
  timeSample = 0.01
else:
  from aradi import *
  method     = hines
  timeSample = 0
from Exp2Syn import *


# Units
Ohm    = 1
Siemen = 1
coul   = 1
mole   = 1e3
meter  = 1
litre  = 1
cm     = meter*1e-2
um     = meter*1e-6
nm     = meter*1e-9
mS     = Siemen*1e-3
ms     = 1  # Time is in ms!!!
uF     = 1e-6
mM     = mole*1e-3/litre
# Current is in units of coul/ms !
nA     = 1e-9 * 1e3
pA     = 1e-12 * 1e3

# Constants
FARADAY = 96520*coul
depth = 200*nm

class DGcmpt(Compartment):
	
  def __init__(self, cell, parent, Gleak, Ra, Cm, length, dia, Vleak=-70):
    Compartment.__init__(self, cell)
    		
    # unit assignment
    Cm     = Cm*uF/cm/cm
    dia    = dia*um
    length = length*um
    Ra     = Ra*Ohm*cm
    Gleak  = Gleak*mS/cm/cm
    
    self.Em          = Vleak
    # Assuming an open cylinder - why?
    self.area        = pi*dia*length
    self.capacitance = Cm * self.area * 1e3
    self.Ra          = Ra
    self.dia         = dia
    self.length      = length
    self.leak        = Leak(self)
    self.leak.Er     = self.Em
    self.leak.Gmax   = self.area*Gleak

    if parent:
      r1 = 1/(pi*dia*dia/2.0/Ra/length)
      r2 = 1/(pi*parent.dia*parent.dia/2.0/Ra/parent.length)
      self.connect(parent, 1/(r1+r2))
      
  def addG(self, useSlow=False, gNa=0, gKdrf=0, gKdrs=0, gKA=0,
           gCaT=0, gCaN=0, gCaL=0, gBK=0, gSK=0, gIh=0):
    area = self.area*mS/cm/cm
    ENa = 55
    ECa = 130
    EK  = -90
    Eh  = -40
    if gNa>0 and not useSlow:
      self.Na        = Na(self)
      self.Na.Er     = ENa
      self.Na.Gmax   = gNa*area
      self.Na.m      = self.Na.m_inf(self.Em)
      self.Na.h      = self.Na.h_inf(self.Em)
    if gNa>0 and useSlow:
      self.Na        = Naslow(self)
      self.Na.Er     = ENa
      self.Na.Am     = 1
      self.Na.Ah     = 1
      self.Na.As     = 1
      self.Na.Vhalfm = 0
      self.Na.Vhalfh = 0
      self.Na.Vhalfs = 0
      self.Na.m      = self.Na.m_inf(self.Em)
      self.Na.h      = self.Na.h_inf(self.Em)
      self.Na.s      = self.Na.s_inf(self.Em)
      self.Na.Gmax   = gNa*area/self.Na.s
    if gKdrf>0:
      self.Kdrf      = Kdrf(self)
      self.Kdrf.Er   = EK
      self.Kdrf.Gmax = gKdrf*area
      self.Kdrf.n    = self.Kdrf.n_inf(self.Em)
    if gKdrs>0:
      self.Kdrs      = Kdrs(self)
      self.Kdrs.Er   = EK
      self.Kdrs.Gmax = gKdrs*area
      self.Kdrs.n    = self.Kdrs.n_inf(self.Em)
    if gKA>0:
      self.Ka        = Ka(self)
      self.Ka.Er     = EK
      self.Ka.Gmax   = gKA*area
      self.Ka.k      = self.Ka.k_inf(self.Em)
      self.Ka.l      = self.Ka.l_inf(self.Em)

    # IntCa is always present!
    self.IntCa        = CaDynamics(self)
    self.IntCa.charge = 2
    self.IntCa.B      = 1/depth/FARADAY/area/1e5
    self.IntCa.Catau  = 10
    self.IntCa.Ca0    = 5.e-6 * mM
    try:
      self.IntCa.Ca     = self.IntCa.Ca0
    except AttributeError:
      self.IntCa.CaLi   = self.IntCa.Ca0/3
      self.IntCa.CaNi   = self.IntCa.Ca0/3
      self.IntCa.CaTi   = self.IntCa.Ca0/3
    self.IntCa.co     = 2*mM
    
    if gCaN>0:
      self.CaN       = CaN(self)
      self.CaN.Er    = ECa
      self.CaN.Gmax  = gCaN*area
      self.IntCa.CaProviders.append(self.CaN)
      self.CaN.useGHKrectification = True
      self.CaN.c     = self.CaN.c_inf(self.Em)
      self.CaN.d     = self.CaN.d_inf(self.Em)
      self.CaN.Cadynamics = self.IntCa
    if gCaL>0:
      self.CaL            = CaL(self)
      self.CaL.Er         = ECa
      self.CaL.Gmax       = gCaL*area
      self.CaL.e          = self.CaL.e_inf(self.Em)
      self.CaL.ki         = 1*mM
      self.CaL.useGHKrectification = True
      self.IntCa.CaProviders.append(self.CaL)
      self.CaL.Cadynamics = self.IntCa
    if gCaT>0:
      self.CaT            = CaT(self)
      self.CaT.Er         = ECa
      self.CaT.Gmax       = gCaT*area
      self.CaT.a          = self.CaT.a_inf(self.Em)
      self.CaT.b          = self.CaT.b_inf(self.Em)
      self.CaT.useGHKrectification = True
      self.IntCa.CaProviders.append(self.CaT)
      self.CaT.Cadynamics = self.IntCa
    if gBK>0:
      self.Kbk       = Kbk(self)
      self.Kbk.Er    = EK
      self.Kbk.Gmax  = gBK*area
      self.Kbk.o     = self.Kbk.o_inf(self.Em, self.IntCa.Ca0)
      self.Kbk.Cadynamics = self.IntCa
    if gSK>0:
      self.Ksk       = Ksk(self)
      self.Ksk.Er    = EK
      self.Ksk.Gmax  = gSK*area
      self.Ksk.q     = 0
      self.Ksk.Cadynamics = self.IntCa
    if gIh>0:
      self.Ihs       = Ihs(self)
      self.Ihs.Er    = Eh
      self.Ihs.Gmax  = gIh*area
      self.Ihs.h     = self.Ihs.h_inf(self.Em)
      self.Ihf       = Ihf(self)
      self.Ihf.Er    = Eh
      self.Ihf.Gmax  = gIh*area
      self.Ihf.h     = self.Ihf.h_inf(self.Em)


class SquareWaveCurrent(PyDynamics):
  def __init__(self, cell):
    PyDynamics.__init__(self, cell)
    self.HinesIntegrable = True
    self.Iinject = 0
    self.dur = 40
    
  def current(self, t):
    if t%self.dur<=self.dur/2:
      return 0
    else:
      return self.Iinject

class Granule(mCell):
  P = [0.2409,    1.0000,    0.4832,    0.0000,    1.0000]

  def __str__(self):
    return 'Granule cell'

  def __init__(self, useSlow=False):
    mCell.__init__(self)

    if old:
      self.CaLscale = 1
      self.CaNscale = 1
      self.CaTscale = 1
      self.BKscale  = 1
      self.SKscale  = 1
    else:
      self.CaLscale = self.P[0]
      self.CaNscale = self.P[1]
      self.CaTscale = self.P[2]
      self.BKscale  = self.P[3]
      self.SKscale  = self.P[4]

    self.method = method
    self.sparseLinearAlgebra = False
    self.timeSample = timeSample

    # Soma (table 1 in the paper indicates 2 somatic compartments ??)
    self.soma    = DGcmpt(self, None, 0.04, 210, 1, 16.8, 16.8)
    self.soma.addG(useSlow=useSlow, gNa=120, gKdrf=16, gKdrs=6,
                   gKA=12, gCaN=2*self.CaNscale, gCaL=5*self.CaLscale,
                   gCaT=0.037*self.CaTscale,
                   gSK=1*self.SKscale, gBK=0.6*self.BKscale)

    # Inject current into the soma
    #d = SquareWaveCurrent(self.soma)
    d = CurrentInjection(self.soma)
    d.start = 0
    d.end   = 1000000000
    self.current = d

    self.soma.APthreshold = 10

    self.synlist = []
    self.PPlist  = []

    # dendrites
    self.dendrite = []
    for i in range(2):
      self.dendrite.append([])
      # GCL
      d = DGcmpt(self, self.soma, 0.04, 210, 1, 50, 3)
      self.dendrite[i].append(d)
      d.addG(useSlow=useSlow, gNa=18, gKdrf=4, gKdrs=6,
             gCaN=3*self.CaNscale, gCaL=7.5*self.CaLscale,
             gCaT=0.075*self.CaLscale,
             gSK=0.4*self.SKscale, gBK=0.6*self.BKscale)
      
      # proximal
      d = DGcmpt(self, d, 0.063, 210, 1.6, 150, 3)
      self.dendrite[i].append(d)
      d.addG(useSlow=useSlow, gNa=13, gKdrf=4, gKdrs=6,
             gCaN=1*self.CaNscale, gCaL=7.5*self.CaLscale,
             gCaT=0.25*self.CaTscale,
             gSK=0.2*self.SKscale, gBK=1*self.BKscale)

      # medial
      d = DGcmpt(self, d, 0.063, 210, 1.6, 150, 3)
      self.dendrite[i].append(d)
      d.addG(useSlow=useSlow, gNa=8, gKdrf=1, gKdrs=6,
             gCaN=1*self.CaNscale, gCaL=0.5*self.CaLscale,
             gCaT=0.5*self.CaTscale,
             gSK=0*self.SKscale, gBK=2.4*self.BKscale)

      # distal
      d = DGcmpt(self, d, 0.063, 210, 1.6, 150, 3)
      self.dendrite[i].append(d)
      d.addG(useSlow=useSlow, gNa=0, gKdrf=1, gKdrs=8,
             gCaN=1*self.CaNscale, gCaL=0*self.CaLscale,
             gCaT=1*self.CaTscale,
             gSK=0*self.SKscale, gBK=2.4*self.BKscale)

  def make1synapse(self, cmpt, tau1, tau2, Er):
    syn = Exp2Syn(cmpt)
    syn.tau1 = tau1; syn.tau2 = tau2; syn.Er = Er; syn.Gmax = nA
    self.synlist.append(syn)
    return syn
  
  def makesynapses(self):
    # PP syn based on data from Greg Hollrigel and Kevin Staley
    d = self.make1synapse(self.dendrite[0][3], 1.5, 5.5, 0)
    d.Gmax = d.Gmax * 0.02
    self.PPlist.append(d)
    d = self.make1synapse(self.dendrite[1][3], 1.5, 5.5, 0)
    d.Gmax = d.Gmax * 0.02
    self.PPlist.append(d)
    # MC syn *** Estimated
    self.make1synapse(self.dendrite[0][1], 1.5, 5.5, 0)
    self.make1synapse(self.dendrite[1][1], 1.5, 5.5, 0)
    # HIPP  syn based on Harney and Jones corrected for temp
    self.make1synapse(self.dendrite[0][3], 0.5, 6.0, -70)
    self.make1synapse(self.dendrite[0][3], 0.5, 6.0, -70)
    # BC  syn syn based on Bartos
    self.make1synapse(self.soma, 0.26, 5.5, -70)
    # NOTE: SPROUTED SYNAPSE based on Molnar and Nadler
    self.make1synapse(self.dendrite[0][1], 1.5, 5.5, 0)
    self.make1synapse(self.dendrite[1][1], 1.5, 5.5, 0)


class Basket(mCell):
  P=[0.25, 1, 0, 1]
  
  def addcmpt(self, parent, useSlow=False, gNa=0, gKdrf=0, length=0, dia=0):
    cmpt = DGcmpt(self, parent, 0.18, 100, 1.4, length, dia, Vleak=-60.06)
    cmpt.addG(useSlow=useSlow, gNa=gNa, gKdrf=gKdrf, gKA=0.15,
              gCaN=0.8*self.CaNscale, gCaL=5*self.CaLscale,
              gSK=0.002*self.SKscale, gBK=0.2*self.BKscale)
    return cmpt

  def __init__(self, useSlow=False):
    mCell.__init__(self)

    if old:
      self.CaLscale = 1
      self.CaNscale = 1
      self.BKscale  = 1
      self.SKscale  = 1
    else:
      self.CaLscale = self.P[0]
      self.CaNscale = self.P[1]
      self.BKscale  = self.P[2]
      self.SKscale  = self.P[3]
    
    self.method = method
    self.sparseLinearAlgebra = False
    self.timeSample = timeSample

    # Soma (table 1 in the paper indicates 2 somatic compartments ??)
    self.soma = self.addcmpt(None, useSlow=useSlow,
                             gNa=120, gKdrf=13, length=20, dia=15)
    
    # Inject current into the soma
    #d = SquareWaveCurrent(self.soma)
    #d = CurrentInjection(self.soma)
    #d.start = 0
    #d.end   = 1000000000
    #self.current = d

    self.soma.APthreshold = 10

    self.synlist  = []
    self.PPlist   = []
    self.dendrite = []

    napdend  = 2
    nbasdend = 2
    
    # apical dendrites
    for i in range(napdend):
      self.dendrite.append([])
      # proximal
      d = self.addcmpt(self.soma, useSlow=useSlow,
                       gNa=120, gKdrf=13, length=75, dia=4)
      self.dendrite[i].append(d)
      
      # mid 1
      d = self.addcmpt(d, useSlow=useSlow,
                       gNa=0, gKdrf=0, length=75, dia=3)
      self.dendrite[i].append(d)
      
      # mid 2
      d = self.addcmpt(d, useSlow=useSlow,
                       gNa=0, gKdrf=0, length=75, dia=2)
      self.dendrite[i].append(d)

      # distal
      d = self.addcmpt(d, useSlow=useSlow,
                       gNa=0, gKdrf=0, length=75, dia=1)
      self.dendrite[i].append(d)
    
    # basal dendrites
    for i in range(nbasdend):
      self.dendrite.append([])
      # proximal
      d = self.addcmpt(self.soma, useSlow=useSlow,
                       gNa=120, gKdrf=13, length=50, dia=4)
      self.dendrite[i+2].append(d)
      
      # mid 1
      d = self.addcmpt(d, useSlow=useSlow,
                       gNa=0, gKdrf=0, length=50, dia=3)
      self.dendrite[i+2].append(d)
      
      # mid 2
      d = self.addcmpt(d, useSlow=useSlow,
                       gNa=0, gKdrf=0, length=50, dia=2)
      self.dendrite[i+2].append(d)

      # distal
      d = self.addcmpt(d, useSlow=useSlow,
                       gNa=0, gKdrf=0, length=50, dia=1)
      self.dendrite[i+2].append(d)

  def __str__(self):
    return 'Basket cell'

  def make1synapse(self, cmpt, tau1, tau2, Er):
    syn = Exp2Syn(cmpt)
    syn.tau1 = tau1; syn.tau2 = tau2; syn.Er = Er; syn.Gmax = nA
    self.synlist.append(syn)
    return syn

  def makesynapses(self):
    # PP(AMPA) syn to apical dist dend Dingledine '95
    d = self.make1synapse(self.dendrite[0][3], 2.0 , 6.3,   0)
    d.Gmax = d.Gmax * 0.01
    self.PPlist.append(d)
    d = self.make1synapse(self.dendrite[1][3], 2.0 , 6.3,   0)
    d.Gmax = d.Gmax * 0.01
    self.PPlist.append(d)
    # GC(AMPA) syn to prox dend Geiger '97
    self.make1synapse(self.dendrite[0][0], 0.3 , 0.6,   0)
    self.make1synapse(self.dendrite[1][0], 0.3 , 0.6,   0)
    self.make1synapse(self.dendrite[2][0], 0.3 , 0.6,   0)
    self.make1synapse(self.dendrite[3][0], 0.3 , 0.6,   0)
    # MC(AMPA) syn to apical IML dend
    # *** Estimated based on CA3>BC min stim Dingledine '95
    self.make1synapse(self.dendrite[0][1], 0.9 , 3.6,   0)
    self.make1synapse(self.dendrite[1][1], 0.9 , 3.6,   0)
    # BC(GABA) syn to apical IML dend Bartos
    self.make1synapse(self.dendrite[0][1], 0.16, 1.8, -70)
    self.make1synapse(self.dendrite[1][1], 0.16, 1.8, -70)
    # HIPP(GABA) syn to apical distal dend
    # *** Estimated as HIPP>GC
    self.make1synapse(self.dendrite[0][3], 0.4 , 5.8, -70)
    self.make1synapse(self.dendrite[1][3], 0.4 , 5.8, -70)

class Mossy(mCell):
  P = [0, 0, 1, 0.25]

  def addcmpt(self, parent, useSlow=False, gNa=0, gKdrf=0, length=0, dia=0, Cm=1, gl=1):
    cmpt = DGcmpt(self, parent, gl, 100, Cm, length, dia, Vleak=-59)
    cmpt.addG(useSlow=useSlow, gNa=gNa, gKdrf=gKdrf, gKA=0.01,
              gCaN=0.08*self.CaNscale, gCaL=0.6*self.CaLscale,
              gSK=16*self.SKscale, gBK=16.5*self.BKscale, gIh=0.005)
              
    return cmpt

  def __init__(self, useSlow=False):
    mCell.__init__(self)

    if old:
      self.CaLscale = 1
      self.CaNscale = 1
      self.BKscale  = 1
      self.SKscale  = 1
    else:
      self.CaLscale = self.P[0]
      self.CaNscale = self.P[1]
      self.BKscale  = self.P[2]
      self.SKscale  = self.P[3]

    self.method = method
    self.sparseLinearAlgebra = False
    self.timeSample = timeSample

    # Soma (table 1 in the paper indicates 2 somatic compartments ??)
    self.soma = self.addcmpt(None, useSlow=useSlow,
                             gNa=120, gKdrf=0.5, length=20, dia=20,
                             Cm=0.6, gl=0.011)

    # Inject current into the soma
    d = CurrentInjection(self.soma)
    d.start = 0
    d.end   = 1000000000
    self.current = d

    self.soma.APthreshold = 10

    self.synlist = []
    self.PPlist  = []

    self.dendrite = []
    # dendrites
    for i in range(4):
      self.dendrite.append([])
      # proximal
      d = self.addcmpt(self.soma, useSlow=useSlow,
                       gNa=120, gKdrf=0.5, length=50, dia=5.78,
                       Cm=2.4, gl=0.044)
      self.dendrite[i].append(d)

      dia = [4, 2.5, 1]
      for j in range(3):
        # mid 1
        d = self.addcmpt(d, useSlow=useSlow,
                         gNa=0, gKdrf=0, length=50, dia=dia[j],
                         Cm=2.4, gl=0.044)
        self.dendrite[i].append(d)

  def __str__(self):
    return 'Mossy cell'

  def make1synapse(self, cmpt, tau1, tau2, Er):
    syn = Exp2Syn(cmpt)
    syn.tau1 = tau1; syn.tau2 = tau2; syn.Er = Er; syn.Gmax = nA
    self.synlist.append(syn)
    return syn

  def makesynapses(self):
    # PP(AMPA) syn to dist dend similar to PP to GC
    d = self.make1synapse(self.dendrite[0][3], 1.5 , 5.5,   0)
    d.Gmax = d.Gmax * 0.005
    self.PPlist.append(d)
    d = self.make1synapse(self.dendrite[1][3], 1.5 , 5.5,   0)
    d.Gmax = d.Gmax * 0.005
    self.PPlist.append(d)
    d = self.make1synapse(self.dendrite[2][3], 1.5 , 5.5,   0)
    d.Gmax = d.Gmax * 0.005
    self.PPlist.append(d)
    d = self.make1synapse(self.dendrite[3][3], 1.5 , 5.5,   0)
    d.Gmax = d.Gmax * 0.005
    self.PPlist.append(d)
    # GC(AMPA) syn to prox dend similar to GC>CA3 Jonas '93
    self.make1synapse(self.dendrite[0][0], 0.5 , 6.2,   0)
    self.make1synapse(self.dendrite[1][0], 0.5 , 6.2,   0)
    self.make1synapse(self.dendrite[2][0], 0.5 , 6.2,   0)
    self.make1synapse(self.dendrite[3][0], 0.5 , 6.2,   0)
    # MC(AMPA) syn to prox dend similar to CA#>CA3 Aaron
    self.make1synapse(self.dendrite[0][0], 0.45, 2.2,   0)
    self.make1synapse(self.dendrite[1][0], 0.45, 2.2,   0)
    self.make1synapse(self.dendrite[2][0], 0.45, 2.2,   0)
    self.make1synapse(self.dendrite[3][0], 0.45, 2.2,   0)
    # BC(GABA) syn to prox dend based on BC>CA3 Bartos PNAS (mice)
    self.make1synapse(self.soma,           0.3 , 3.3, -70)
    # HIPP(GABA) syn to prox dend based on Hilar>GC Harney&Jones
    self.make1synapse(self.dendrite[0][2], 0.5 , 6.0, -70)
    self.make1synapse(self.dendrite[1][2], 0.5 , 6.0, -70)
    self.make1synapse(self.dendrite[2][2], 0.5 , 6.0, -70)
    self.make1synapse(self.dendrite[3][2], 0.5 , 6.0, -70)
    
    
class Hipp(mCell):
  P=[0, 1, 0.7]

  def addcmpt(self, parent, useSlow=False, gNa=0, gKdrf=0,
              length=0, dia=0, Cm=1.1, gl=0.036):
    cmpt = DGcmpt(self, parent, gl, 100, Cm, length, dia, Vleak=-70.45)
    cmpt.addG(useSlow=useSlow, gNa=gNa, gKdrf=gKdrf, gKA=0.8,
              gCaL=1.5*self.CaLscale,
              gSK=3*self.SKscale, gBK=3*self.BKscale, gIh=0.015)
              
    return cmpt

  def __init__(self, useSlow=False):
    mCell.__init__(self)

    if old:
      self.CaLscale = 1
      self.BKscale  = 1
      self.SKscale  = 1
    else:
      self.CaLscale = self.P[0]
      self.BKscale  = self.P[1]
      self.SKscale  = self.P[2]
      
    self.method = method
    self.sparseLinearAlgebra = False
    self.timeSample = timeSample

    # Soma (table 1 in the paper indicates 2 somatic compartments ??)
    self.soma = self.addcmpt(None, useSlow=useSlow,
                             gNa=200, gKdrf=6, length=20, dia=10)

    # Inject current into the soma
    #d = SquareWaveCurrent(self.soma)
    d = CurrentInjection(self.soma)
    d.start = 0
    d.end   = 1000000000
    self.current = d

    self.PPlist  = []

    self.soma.APthreshold = 10

    self.dendrite = []
    length = [75, 75, 50, 50]
    dia = [3, 2, 1]
    gNa = [200, 0, 0]
    gK  = [6, 0, 0]
    # dendrites
    for i in range(4):
      self.dendrite.append([])
      d = self.soma
      for j in range(3):
        d = self.addcmpt(d, useSlow=useSlow, gNa=gNa[j], gKdrf=gK[j],
                         length=length[i], dia=dia[j])
        self.dendrite[i].append(d)

  def __str__(self):
    return 'HIPP cell'

  def make1synapse(self, cmpt, tau1, tau2, Er):
    syn = Exp2Syn(cmpt)
    syn.tau1 = tau1; syn.tau2 = tau2; syn.Er = Er; syn.Gmax = nA
    self.synlist.append(syn)
    return syn

  def makesynapses(self):
    self.synlist = []
    # GC(AMPA) syn to prox dend similar to GC>BC
    self.make1synapse(self.dendrite[0][0], 0.3 , 0.6,   0)
    self.make1synapse(self.dendrite[1][0], 0.3 , 0.6,   0)
    self.make1synapse(self.dendrite[2][0], 0.3 , 0.6,   0)
    self.make1synapse(self.dendrite[3][0], 0.3 , 0.6,   0)
    # MC(AMPA) syn to mid dend similar to CA3>int Aaron
    # *** Assumed data at physio temp
    self.make1synapse(self.dendrite[0][1], 0.9 , 3.6,   0)
    self.make1synapse(self.dendrite[1][1], 0.9 , 3.6,   0)
    self.make1synapse(self.dendrite[2][1], 0.9 , 3.6,   0)
    self.make1synapse(self.dendrite[3][1], 0.9 , 3.6,   0)

def modifyAll(Network=[], \
    Vhalfnf=0, Vhalfns=0, \
    Vhalfm=0, Vhalfh=0, Vhalfs=0, \
    Am=1, Ah=1, As=1, \
    Ggaba=1, \
    Vhalfma=0, Vhalfha=0, \
    Vhalfmt=0, \
    Vhalfml=0, \
    Vhalfhf=0, \
    Vhalfhs=0, \
    Vhalfmn=0, \
    ):
  for cell in Network:
    for cmpt in cell.compartments:
      try:
        cmpt.Na.Vhalfm = Vhalfm
        cmpt.Na.Vhalfh = Vhalfh
        cmpt.Na.Vhalfs = Vhalfs
        cmpt.Na.Am     = Am
        cmpt.Na.Ah     = Ah
        cmpt.Na.As     = As
      except AttributeError:
        pass
      try:
        cmpt.Ka.vhalfn = cmpt.Ka.vhalfn-Vhalfma
        cmpt.Ka.vhalfl = cmpt.Ka.vhalfl-Vhalfha
      except AttributeError:
        pass
      try:
        cmpt.Ihf.vhalfhf = Vhalfhf
      except AttributeError:
        pass
      try:
        cmpt.Ihs.vhalfhs = Vhalfhs
      except AttributeError:
        pass
      try:
        cmpt.Kdrf.Vhalfn = Vhalfnf
      except AttributeError:
        pass
      try:
        cmpt.CaL.Vhalfml = Vhalfml
      except AttributeError:
        pass
      try:
        cmpt.CaN.Vhalfmn = Vhalfmn
      except AttributeError:
        pass
      try:
        cmpt.CaT.Vhalfmt = Vhalfmt
      except AttributeError:
        pass
      try:
        cmpt.Kdrs.Vhalfn = Vhalfns
      except AttributeError:
        pass
    for syn in cell.synlist:
        if syn.Er==-70:
            syn.Gmax = syn.Gmax*Ggaba

def nc_append(src, tgt, synapse, weight, delay, threshold, cnxtype=''):
    from __main__ import Network, ngcell, nbcell, nhcell, nmcell, sprout
    from random import random
    Ncells = len(Network)

    tgtcell = Network[tgt]
    tgtdyn  = tgtcell.synlist[synapse]
    cmpt    = tgtdyn.owner
    cmpt.APthreshold = threshold
    if src<Ncells:
        srccell = Network[src]
        if isinstance(tgtcell, Granule) and isinstance(srccell, Granule):
            if sprout<=random():
                return
        s = Synapse(srccell, tgtdyn)
        s.trans_time = delay
        s.nominal_strength = weight
    else:
        if isinstance(tgtcell, Granule):
            tgtcell = Network[tgt+200]
            tgtdyn  = tgtcell.synlist[synapse]
        tgtdyn.enq(delay)
