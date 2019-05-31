#!/usr/bin/python

import sys, string, os

def writecmd(comment):
    print 'Running: ', comment
    os.system('python sprout_run.py %s\n' % comment)
    
for sprout in [0.4, 0.5, 0.6, 0.7, 1]:
    cmd = 'sprout=%g' % sprout
    writecmd(cmd)
    cmd = 'Vhalfmn=2 sprout=%g' % sprout
    writecmd(cmd)
    cmd = 'Vhalfm=-2 Vhalfmn=2 sprout=%g' % sprout
    writecmd(cmd)
    cmd = 'Vhalfm=-2 Vhalfns=2 Vhalfmn=2 sprout=%g' % sprout
    writecmd(cmd)

