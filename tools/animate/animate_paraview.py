#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='Animate Paraview data')
parser.add_argument('-l','--layout',  help='Paraview State file (.pvsm-file)')
parser.add_argument('-r','--reader',  help='Path to custom reader plugin')
parser.add_argument('-n','--nomovie', help='Only generate pictures', action='store_true')
parser.add_argument('-f','--fps', type=int, default=10, help='Frames per second')
parser.add_argument('-b','--bitrate', type=int, default=10000, help='Bitrate of movie')
parser.add_argument('-m','--mpi', type=int, default=1, help='Number of MPI procs for rendering')
parser.add_argument('-s','--scale', type=int, default=1, help='Magnification of rendering')
parser.add_argument('-o','--output', default='',  help='Appendix for filenames')
parser.add_argument('plotfiles', nargs='+', help='Files to animate (.vtu/.pvtu-files)')

args = parser.parse_args()

plotfiles = [f for f in args.plotfiles if (os.path.splitext(f)[1] in ['.pvtu', '.vtu', '.plt', '.vtm', '.h5']) ]

i = 0
for p in plotfiles : 
    i = i+1
    sys.stdout.write('\r%05.2f %% Animate: %s' % (100.0 * i / len(plotfiles), p))
    sys.stdout.flush()
    fn = os.path.splitext(p)[0]+'.py'
    f = open(fn, 'w')
    f.write("""from paraview.simple import * 
import os

paraview.simple._DisableFirstRenderCameraReset()
""")
    if args.reader :
        f.write("""servermanager.LoadPlugin('%s')
""" % (args.reader))

    f.write("""servermanager.LoadState('%s')
statefilename = GetSources() 
plotfilename = None
for k in statefilename.keys() :
    if os.path.splitext(k[0])[1] in ['.pvtu', '.vtu', '.plt', '.vtm', '.h5'] :
        plotfilename = k[0]
        break

if not plotfilename : exit(1)

reader = FindSource(plotfilename)
reader.FileName = ['%s'] 
reader.FileNameChanged() 
SetActiveView(GetRenderView())
Render()
WriteImage('%s',  Magnification=%d)
""" % (args.layout, p, os.path.splitext(p)[0] + args.output + '.png', args.scale))
    f.close()
    if args.mpi > 1 :
        cmd = ['mpirun', '-np', str(args.mpi), 'pvbatch', '--use-offscreen-rendering', fn]
    else :
        cmd = ['pvbatch', '--use-offscreen-rendering', fn]
    p = subprocess.Popen(cmd)
    p.wait()
    os.remove(fn)
sys.stdout.write('\n')


if not args.nomovie :
    print 'Generate movie ....'
    cmd = ['mencoder']
    cmd.append('mf://%s/*%s.png' % (os.path.dirname(os.path.abspath(args.layout)),args.output ))
    cmd.append('-mf')
    cmd.append('fps=%d' % args.fps)
    cmd.append('-o')
    cmd.append('%s' % os.path.splitext(args.layout)[0] + '.avi')
    cmd.append('-ovc')
    cmd.append('lavc')
    cmd.append('-lavcopts')
    cmd.append('vcodec=msmpeg4v2:vbitrate=%d' % args.bitrate)
    p = subprocess.Popen(cmd)
    p.wait()
    print 'done'

