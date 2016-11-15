#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

parser = argparse.ArgumentParser(description='Concatenate pictures')
parser.add_argument('-d','--direction', choices=['n','e','s','w'], default='e', help='Direction where to append picture (north, east, south, west)')
parser.add_argument('-p','--pics', nargs='+', help='Base pictures (to each of this pictures, the corresponding of "--appends" will be appended)')
parser.add_argument('-a','--appends', nargs='+', help='Append pictures')

args = parser.parse_args()


no = 0
for a,b in zip(sorted(args.pics), sorted(args.appends)) :
   first = None
   last = None
   for i in range(len(a)) :
      if not first  :
         if a[i] != b[i] :
            first = i
            break

   for i in range(len(a)) :
      if not last  :
         if a[len(a)-1-i] != b[len(b)-1-i] :
            last = i
            break

   filename = a[0:len(a)-last] + '_' + b[first:len(b)-last] + a[len(a)-last:]
   no = no+1
   sys.stdout.write('\r%05.2f %% Process: %s and %s' % (100.0 * no / len(args.pics), a, b))
   sys.stdout.flush()

   cmd = ['convert']
   if args.direction in 'ns' :
      cmd.append('-append')
   else :
      cmd.append('+append')

   if args.direction in 'es' : 
      cmd.append(a)
      cmd.append(b)
   else : 
      cmd.append(b)
      cmd.append(a)

   cmd.append(filename)
   p = subprocess.Popen(cmd)
   p.wait()

sys.stdout.write('\n')
