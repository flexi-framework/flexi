#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='Animate Paraview data')
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument('-l','--layout',  help="Paraview State file (.pvsm-file); gained via 'Save State...' in ParaView", required=True)
optional.add_argument('-r','--reader',  help="Path to custom reader plugin (e.g. 'path/to/libvisuReader.so'), required for .h5 files")
optional.add_argument('-n','--nomovie', help='Only generate pictures', action='store_true')
optional.add_argument('-f','--fps',     type=int,   default=10,    help='Frames per second')
optional.add_argument('-b','--bitrate', type=int,   default=10000, help='Bitrate of movie')
optional.add_argument('-m','--mpi',     type=int,   default=1,     help='Number of MPI procs for rendering')
optional.add_argument('-s','--scale',   type=int,   default=1,     help='Magnification of rendering')
optional.add_argument('-o','--output',  default='', help='Appendix for filenames')
optional.add_argument('-x','--folder',  default='', help='relative output folder for images and movies')
required.add_argument('plotfiles',      nargs='+',  help='Files to animate (.vtu/.pvtu-files)')
parser._action_groups.append(optional)


args = parser.parse_args()

# if output folder does not exist, create it
#cwd = os.getcwd()
fp=os.path.join(os.getcwd(), args.folder)
if not os.path.exists(fp):
    os.makedirs(fp)

plotfiles = [f for f in args.plotfiles if (os.path.splitext(f)[1] in ['.pvtu', '.vtu', '.plt', '.vtm', '.h5']) ]

has_h5_plotfiles = any([(os.path.splitext(f)[1] == '.h5') for f in plotfiles])
if has_h5_plotfiles and not args.reader :
   sys.exit("Please specifiy path to reader plugin (e.g. '-r path/to/libvisuReader.so')  if input is HDF5!")

i = 0
for p in plotfiles :
    i = i+1
    sys.stdout.write('\r%05.2f %% Animate: %s' % (100.0 * i / len(plotfiles), p))
    sys.stdout.flush()
    fn = os.path.splitext(p)[0]+ args.output + '.py'
    f = open(fn, 'w')
    # get filename
    of = os.path.splitext(p)[0] + args.output + '.png'
    of2=os.path.basename(of)
    # create output filename
    of=os.path.join(os.getcwd(), args.folder, of2)
    f.write("""from paraview.simple import *
import os

paraview.simple._DisableFirstRenderCameraReset()
""")
    if args.reader :
        f.write("servermanager.LoadPlugin('%s')\n" % (args.reader))

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
RenderView1 = GetRenderView()
if RenderView1.InteractionMode == "2D" :
    RenderView1.CameraParallelProjection=1
SetActiveView(RenderView1)
Render()
WriteImage('%s',  Magnification=%d)
""" % (args.layout, p, of, args.scale))
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
    cmd.append('mf://%s/*%s.png' % (fp,args.output ))
    cmd.append('-mf')
    cmd.append('fps=%d' % args.fps)
    cmd.append('-o')
    cmd.append('%s' % os.path.join(fp, os.path.basename(os.path.splitext(args.layout)[0])) + '.avi')
    cmd.append('-ovc')
    cmd.append('lavc')
    cmd.append('-lavcopts')
    cmd.append('vcodec=msmpeg4v2:vbitrate=%d' % args.bitrate)
    p = subprocess.Popen(cmd)
    p.wait()
    print 'done'

