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

"""-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
     This script is used to build and install user C language extensions
-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"""

import sys, re
from distutils.core import setup, Extension

def guess_modname(files):
    exp = '\s*MAKE_P3_MODEL\s*\(\s*(.+)\s*,'
    for fn in files:
        f = open(fn)
        ls = f.readlines()
        for l in ls:
            o = re.search(exp, l)
            if o:
                return o.group(1)
        

if sys.argv[1]=='-d':
    debug = True
    del sys.argv[1]
else:
    debug = False

srcfiles = sys.argv[1:]

name = guess_modname(srcfiles)

if debug:
    sys.argv[1:] = ['build_ext', '--inplace', '--debug']
else:
    sys.argv[1:] = ['build_ext', '--inplace']


defmac = [('MPI', 1)]
E = Extension(name,  srcfiles, define_macros=defmac)
setup(ext_modules=[E])

if debug:
    import os
    # When using the debug flag the load module will be
    # called name_d.pyd. Rename it to name.pyd, but don't
    # rename the debug symbol thingy name_d.pdb
    os.system('mv %s_d.pyd %s.pyd' % (name, name))
