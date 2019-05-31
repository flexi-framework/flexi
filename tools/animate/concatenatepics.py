#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

parser = argparse.ArgumentParser(description='Concatenate pictures')
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional.add_argument('-d','--direction', choices=['n','e','s','w'], default='e', help='Direction where to append picture (north, east, south, west)')
required.add_argument('-p','--pics', nargs='+', help='Base pictures (to each of this pictures, the corresponding of "--appends" will be appended)',required=True)
required.add_argument('-a','--appends', nargs='+', help='Append pictures',required=True)
parser._action_groups.append(optional)

args = parser.parse_args()


no = 0
for a,b in zip(sorted(args.pics), sorted(args.appends)) :
   first = None
   last = None
   af=os.path.basename(a)
   bf=os.path.basename(b)
   for i in range(len(af)) :
      if not first  :
         if af[i] != bf[i] :
            first = i
            break

   for i in range(len(af)) :
      if not last  :
         if af[len(af)-1-i] != bf[len(bf)-1-i] :
            last = i
            break
   filename =os.path.join(os.getcwd() ,af[0:len(af)-last] + '_' + bf[first:len(bf)-last] + af[len(af)-last:])
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
