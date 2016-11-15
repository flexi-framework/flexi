#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

parser = argparse.ArgumentParser(description='Merge pictures to movie')
parser.add_argument('-t','--trim', help='Trim pictures before merge', action='store_true')
parser.add_argument('-f','--fps', type=int, default=10, help='Frames per second')
parser.add_argument('-b','--bitrate', type=int, default=10000, help='Bitrate of movie')
parser.add_argument('-o','--output', type=str, default='output.avi', help='Output-Filename')
parser.add_argument('pictures', nargs='+', help='Picture to animate')

args = parser.parse_args()

tmpdir = tempfile.mkdtemp()

try :
   ext = os.path.splitext(args.pictures[0])[1]
except :
   exit(1)

cmd = ['identify']
cmd.append('-format')
cmd.append('%@')
cmd.append(args.pictures[0])
p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
p.wait()
dimension = p.stdout.read().strip()

i = 0
for p in args.pictures : 
   i = i+1
   sys.stdout.write('\r%05.2f %% Process: %s' % (100.0 * i / len(args.pictures), p))
   sys.stdout.flush()
   if args.trim :
      cmd = ['convert']
      cmd.append(p)
      cmd.append('-crop')
      cmd.append(dimension)
   else :
      cmd = ['cp']
      cmd.append(p)
   cmd.append(os.path.join(tmpdir, os.path.basename(p)))
   p = subprocess.Popen(cmd)
   p.wait()

sys.stdout.write('\n')
    

print 'Generate movie ....'
cmd = ['mencoder']
cmd.append('mf://%s/*%s' % (tmpdir ,ext))
cmd.append('-mf')
cmd.append('fps=%d' % args.fps)
cmd.append('-o')
cmd.append('%s' % args.output)
cmd.append('-ovc')
cmd.append('lavc')
cmd.append('-lavcopts')
cmd.append('vcodec=msmpeg4v2:vbitrate=%d' % args.bitrate)
p = subprocess.Popen(cmd)
p.wait()
print 'done'

shutil.rmtree(tmpdir)

