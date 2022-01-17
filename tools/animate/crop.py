#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

parser = argparse.ArgumentParser(description='Crop pictures to same size')
parser.add_argument('pictures', nargs='+', help='Picture to crop')

args = parser.parse_args()


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
   cmd = ['convert']
   cmd.append(p)
   cmd.append('-crop')
   cmd.append(dimension)
   tmp = os.path.splitext(p)
   cmd.append(tmp[0] + '_crop' + tmp[1])
   p = subprocess.Popen(cmd)
   p.wait()

sys.stdout.write('\n')
