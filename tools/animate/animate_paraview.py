#!/usr/bin/env python3
# -*- coding: utf8 -*-

import os
import sys
import shutil
import argparse
import subprocess

FILETYPES = ['.pvtu', '.vtu', '.plt', '.vtm', '.h5']
"""Supported filetypes for rendering with ParaView."""


def parse_commandline_args():
    """Adds commandline arguments and parses them."""
    parser = argparse.ArgumentParser(description='Animate Paraview data')
    parser.add_argument('-l','--layout',  help="Paraview State file (.pvsm-file); gained via 'Save State...' in ParaView", required=True)
    parser.add_argument('-r','--reader',  help="Path to custom reader plugin (e.g. 'path/to/libvisuReader.so'), required for .h5 files")
    parser.add_argument('-n','--nomovie', help='Only generate pictures', action='store_true')
    parser.add_argument('-t','--trim'   , help='Trim surrounding whitespace', action='store_true')
    parser.add_argument('-f','--fps',     type=int,   default=10,    help='Frames per second')
    parser.add_argument('-b','--bitrate', type=int,   default=10000, help='Bitrate of movie')
    parser.add_argument('-m','--mpi',     type=int,   default=1,     help='Number of MPI procs for rendering')
    parser.add_argument('-s','--scale',   type=int,   default=1,     help='Magnification of rendering')
    parser.add_argument('-o','--output',  default='', help='Appendix for filenames')
    parser.add_argument('-x','--folder',  default='', help='relative output folder for images and movies')
    parser.add_argument('plotfiles',      nargs='+',  help='Files to animate (.vtu/.pvtu-files)')
    return parser.parse_args()


def command_available(command):
    """Check if command is available."""
    if shutil.which(command):
        return True
    return False


def run_command(cmd, abort_on_fail=False):
    """Run command with error handling."""
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(e)
        if abort_on_fail:
            sys.exit(1)


def main(args):
    # Check if all required commands are available
    if not command_available('pvbatch'):
        sys.exit("Please install 'pvbatch' for batch processing with ParaView!")
    if not args.nomovie and not command_available('mencoder'):
        sys.exit("Please install 'mencoder' for generating movies!")
    if args.trim and not command_available('convert'):
        sys.exit("Please install 'convert' to enable trimming of whitespace!")

    # if output folder does not exist, create it
    fp=os.path.join(os.getcwd(), args.folder)
    if not os.path.exists(fp):
        os.makedirs(fp)

    plotfiles = [f for f in args.plotfiles if (os.path.splitext(f)[1] in FILETYPES) ]

    has_h5_plotfiles = any([(os.path.splitext(f)[1] == '.h5') for f in plotfiles])
    if has_h5_plotfiles and not args.reader :
        sys.exit("Please specifiy path to reader plugin (e.g. '-r path/to/libvisuReader.so') if input is HDF5!")

    for i, p in enumerate(plotfiles):
        # print progress
        sys.stdout.write('\r%05.2f %% Animate: %s' % (100. * i / len(plotfiles), p))
        sys.stdout.flush()
        # get output filename
        of = os.path.splitext(os.path.basename(p))[0]    # get filename, remove extension
        of = of + args.output + '.png'                   # add appendix, add png extension
        of = os.path.join(os.getcwd(), args.folder, of)  # add output folder
        # write file
        fn = os.path.splitext(p)[0]+ args.output + '.py'
        with open(fn, 'w', encoding='utf-8') as f:
            f.write(f"""
import os
import sys
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()
servermanager.LoadState('{args.layout}')
statefilename = GetSources()

plotfilename = None
for k in statefilename.keys() :
    if os.path.splitext(k[0])[1] in {FILETYPES}:
        plotfilename = k[0]
        break

if not plotfilename:
    print('No plotfile found in state file!')
    sys.exit(1)

reader = FindSource(plotfilename)
reader.FileName = ['{p}']
reader.FileNameChanged()
RenderView1 = GetRenderView()
if RenderView1.InteractionMode == "2D":
    RenderView1.CameraParallelProjection=1
SetActiveView(RenderView1)
Render()
SaveScreenshot('{of}', magnification={int(args.scale)})
""")

        # create pvbatch command
        cmd = []
        if args.reader:
            cmd.extend(['env', 'PV_PLUGIN_PATH='+str(args.reader)])
        if args.mpi > 1:
            cmd.extend(['mpirun', '-np', str(args.mpi)])
        cmd.extend(['pvbatch', '--force-offscreen-rendering', fn])
        run_command(cmd)

        # trim whitespace
        if args.trim:
            cmd = ['convert', of, '-trim', of]
            run_command(cmd)

        # Remove temporary file
        os.remove(fn)

    if not args.nomovie:
        print('\nGenerate movie ....')
        cmd = ['mencoder']
        cmd.append(f'mf://{fp}/*{args.output}.png')
        cmd.append('-mf')
        cmd.append(f'fps={args.fps}')
        cmd.append('-o')
        cmd.append(f'{os.path.join(fp, os.path.basename(os.path.splitext(args.layout)[0])) + ".avi"}')
        cmd.append('-ovc')
        cmd.append('lavc')
        cmd.append('-lavcopts')
        cmd.append(f'vcodec=msmpeg4v2:vbitrate={args.bitrate}')
        run_command(cmd, abort_on_fail=True)


if __name__ == '__main__':
    args = parse_commandline_args()
    main(args)
