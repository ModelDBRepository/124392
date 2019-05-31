"""
 * Copyright (C) 2004 Evan Thomas
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
"""

"""
Reproduces Fig 6 from Thomas, et al Mossy fiber sprouting interacts with 
sodium channel mutations to increase dentate gyrus excitability, 
Epilepsia (2009)
"""

import sys, os

from p3 import *
from time import clock
from granule import *
from support import *

comment = ''
for c in sys.argv[1:]:
    if c=='' or c==' ': continue
    comment = comment + ' %s' % c
apfn   = 'ap%s' % comment
setAPfilename(apfn)

#// define network size
ngcell = 500
nbcell = 6
nmcell = 15
nhcell = 6

#####################################
#// NETWORK SPECIFICATION INTERFACE #
#####################################
Nstate = 0
Ncmpt  = 0
Gcell  = []
for i in range(ngcell):
    c = Granule(useSlow=True)
    c.makesynapses()
    Nstate = Nstate + len(c.Y)
    Ncmpt  = Ncmpt  + len(c.compartments)
    Gcell.append(c)
Bcell = []
for i in range(nbcell):
    c = Basket(useSlow=True)
    c.makesynapses()
    Nstate = Nstate + len(c.Y)
    Ncmpt  = Ncmpt  + len(c.compartments)
    Bcell.append(c)
Mcell = []
for i in range(nmcell):
    c = Mossy(useSlow=True)
    c.makesynapses()
    Nstate = Nstate + len(c.Y)
    Ncmpt  = Ncmpt  + len(c.compartments)
    Mcell.append(c)
Hcell = []
for i in range(nhcell):
    c = Hipp(useSlow=True)
    c.makesynapses()
    Nstate = Nstate + len(c.Y)
    Ncmpt  = Ncmpt  + len(c.compartments)
    Hcell.append(c)

duration = 400
Vhalfm  = 0
Vhalfnf = 0
Vhalfhf = 0
Vhalfhs = 0
Vhalfns = 0
Vhalfha = 0
Vhalfma = 0
Vhalfmt = 0
Vhalfmn = 0
Vhalfml = 0
Am = 1
Ah = 1
sprout = 0
Ggaba = 1
for s in sys.argv[1:]:
    exec(s)
    
Network = Gcell+Bcell+Mcell+Hcell
s = 'Network with %d cells, %d compartments and %d state variables\n'
message_print(info, s % (len(Network), Ncmpt, Nstate))

modifyAll(Network=Network, Vhalfm=Vhalfm, Am=Am, Ah=Ah, \
    Vhalfns=Vhalfns, Vhalfnf=Vhalfnf, \
    Vhalfmt=Vhalfmt, \
    Vhalfml=Vhalfml, \
    Vhalfmn=Vhalfmn, \
    Vhalfhf=Vhalfhf, \
    Vhalfhs=Vhalfhs, \
    Ggaba=Ggaba, \
    Vhalfma=Vhalfma, Vhalfha=Vhalfha, \
    )

message_print(info, 'Making connections.\n')
execfile('conx_sprout.py')
            
###############
# Run options #
###############
gd = GD()
gd.duration   = duration
gd.tolerance  = 1e-3
gd.minStep    = 0.05
gd.network    = Network
gd.ap_handler = ap_print
gd.trace_handler     = trace_print
gd.endWindow_handler = TimeTicker
gd.stepTrace_handler = None
gd.dumpCell_handler  = dumpcell

#######
# Run #
#######
try:
    start = clock()
    message_print(info, 'Starting run.\n')
    parplex(gd)
    message_print(info, 'made it in %fs\n' % (clock()-start))
except ParplexRuntimeError:
    (e, v) = sys.exc_info()[:2]
    message_print(fatal, 'Caught ParplexRuntimeError: %s\n' % v)
    dumpcell(e.currentCell, fatal)
    message_print(fatal, 'Didn\'t make it\n')
    sys.exit(0)

#try:
#    from py2mat import Matwrap
#    m = Matwrap()
#    m.closeOnDel = False
#    m.write('cd \'' + os.getcwd() + '\'')
#    m.write('figure;traceplot(\'%s.dat\');\n' % trfn)
#    m.write("legend('Granule', 'Basket', 'Mossy', 'HIPP');\n")
#    m.write('title(\'%s\')\n' % comment)
#    m.write('figure;raster(\'%s.dat\');\n' % apfn)
#    m.write('title(\'%s\')\n' % comment)
#except:
#    pass
