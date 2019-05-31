from p3 import *

set_message_option('nodebug')

import os
if os.name == 'nt':
    # Whadya know, windows still can't multitask properly
    try:
        from win32process import SetPriorityClass, GetCurrentProcess, \
             HIGH_PRIORITY_CLASS, NORMAL_PRIORITY_CLASS, \
             REALTIME_PRIORITY_CLASS, IDLE_PRIORITY_CLASS
        # Drop priority so that we can still use the computer
        SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS)
    except:
        pass


##############################################
# Special trace capture so that neurons with #
# differing numbers of compartments can be   #
# handled                                    #
##############################################
tracex = None
trfilename = 'trace'
def setTRfilename(s):
    global trfilename
    trfilename = s

def trace_print(cell):
    global tracex
    tracex = openTraceFile(tracex)

    for j in range(len(cell.compartments)):
        cmpt = cell.compartments[j]
        if not cmpt.emtrace: continue
        for i in range(len(cmpt.traceTimes)):
            time = cmpt.traceTimes[i]
            data = cmpt.traceData[i]
            tracex.write('%d %d %g %g\n' % (cell.id, j, time, data))

    if cell.synlist and cell.synlist[0].trace:
        d = cell.synlist[0]
        cmpt = d.owner
        for i in range(len(cmpt.traceTimes)):
            time = cmpt.traceTimes[i]
            A = d.ATrace[i]
            B = d.BTrace[i]
            tracex.write('%d %d %g %g\n' % (cell.id+1000, 0, time, A))
            tracex.write('%d %d %g %g\n' % (cell.id+1001, 0, time, B))
            
    tracex.flush()

def multi_trace_print(cell):

    if not cell.soma.emtrace: return
    
    if not cell.traceFileHandle:
        cell.traceFileHandle = open(cell.traceFileName, 'w')
        message_print(debug, 'Open of trace file %s successful.\n'
                      % cell.traceFileName)

    cmpt = cell.soma
    for i in range(len(cmpt.traceTimes)):
        time = cmpt.traceTimes[i]
        data = cmpt.traceData[i]
        cell.traceFileHandle.write('%d %d %g %g\n' % (cell.id, 0, time, data))

def multi_ca_trace_print(cell):

    if not cell.soma.IntCa.trace: return
    
    if not cell.CatraceFileHandle:
        cell.CatraceFileHandle = open(cell.CatraceFileName, 'w')
        message_print(debug, 'Open of trace file %s successful.\n'
                      % cell.CatraceFileName)
        
    cmpt  = cell.soma
    IntCa = cell.soma.IntCa
    for i in range(len(cmpt.traceTimes)):
        time = cmpt.traceTimes[i]
        try:
            data = IntCa.CaTrace[i]
        except AttributeError: 
            data = IntCa.CaLiTrace[i]+IntCa.CaNiTrace[i]+IntCa.CaTiTrace[i]
        cell.CatraceFileHandle.write('%d %d %g %g\n' % (cell.id, 0, time, data))


def openTraceFile(tracex):
    if tracex==None:
        if mpi_size == 1:
            fn = trfilename + '.dat'
        else:
            fn = '%s_%d_%d.dat' % \
                 (trfilename, mpi_size, mpi_rank)
        tracex = open(fn, 'w')
        message_print(debug, 'Open of trace file %s successful.\n' % fn)
    return tracex

def stepTrace2(cell):
    """Dump the entire state vector"""
    global tracex
    tracex = openTraceFile(tracex)

    print cell.time
    
    tracex.write('%d %d %g %g' % \
                 (cell.id, cell.laststepaccept, cell.time, cell.step))
    #for x in cell.Y:
    tracex.write(' %g' % cell.Y[0])

    tracex.write('\n')
    tracex.flush()

def stepTrace(cell):
    if cell.id==0 and cell.time==0:
        dumpcell(cell, 'info')

def stepTicker(cell):
    """Tick the integration steps"""
    global stepx

    if cell.laststepaccept:
        s = 'succeess'
    else:
        s = 'failure'

    msg = 'Cell.id=%d soma.Em=%g time=%g step=%g %s\n' % \
          (cell.id, cell.soma.Em, cell.time, cell.step, s)
    
    message_print(debug, msg) 
    stepx.write(msg)
    stepx.flush()

def dumpcell(cell, msglevel):
    """Diagnostic dump of a cell"""
    message_print(msglevel, 'Dumping %s, id=%d, time=%g, step=%g:\n' % \
                  (cell, cell.id, cell.time, cell.step))
    for c in cell.compartments:
        message_print(msglevel, ': %s\n' % c)
        for s in c.synapse:
            message_print(msglevel, '::  %s\n' % s)
    for i in range(len(cell.Y)):
        message_print(msglevel, ':::   Y[%d]=%g DYDT[%d]=%g\n' % (i, cell.Y[i], i, cell.DYDT[i]))
        

def dumpSolverStats(gd, cell, msglevel):
    message_print(msglevel, 'cellid=%d method=%s tol=%g stepTotal=%d stepAccepts=%d (%g%%) functionCnt=%d jacobianCnt=%d newtonCnt=%d sparseLinearAlgebra=%s\n' % \
                  (cell.id, cell.method, gd.tolerance, cell.stepTotal, cell.stepAccepts, \
                   float(cell.stepAccepts)/float(cell.stepTotal)*100.0, \
                   cell.functionCnt, cell.jacobianCnt, cell.newtonCnt, \
                   cell.sparseLinearAlgebra))

def TimeTicker(gd):
    time = (gd.windowID+1)*gd.window
    if time%1000 == 0:
        message_print(info, '%d (ms) elapsed\n' % time)

def stepTicker(cell):
    print cell.time, cell.step, cell.soma.Em, cell.laststepaccept


def CaDumper(cell):
    """Dump the Ca concentration"""
    print 't=%g  Em=%g' % (cell.time, cell.soma.Em)
    for i in [0]:#range(len(cell.compartments)):
      cmpt = cell.compartments[i]
      ICa = 0
      t = cell.time
      Ca = 0
      q = 0
      for dyn in cmpt.dynamics:
        if isinstance(dyn, CaDynamics):
          Ca = dyn.Ca
        if isinstance(dyn, CaL):
          pass
        if isinstance(dyn, CaN):
          pass
        if isinstance(dyn, CaT):
          pass
        if isinstance(dyn, Ksk):
          q = dyn.q
      step.write('%g %g %g %g\n' % (t, cmpt.Em, Ca, q))


apx = None
apfilename = 'ap'
def setAPfilename(s):
    global apfilename
    apfilename = s

def openAPFile(apx):
    if apx==None:
        if mpi_size == 1:
            fn = apfilename + '.dat'
        else:
            fn = '%s_%d_%d.dat' % \
                 (apfilename, mpi_size, mpi_rank)
        apx = open(fn, 'w')
        message_print(debug, 'Open of AP file %s successful.\n' % fn)
    return apx

def ap_print(cell):
    global apx, apfilename
    apx = openAPFile(apx)

    for tm in cell.soma.APtimes:
        apx.write('%d %g "%s"\n' % (cell.id, tm, cell))

    apx.flush()
