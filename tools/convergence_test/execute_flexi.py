#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import subprocess

def execute_flexi(flexi_path, prm_path, projectname, analyze_fcts = None, log = True, \
      mpi_procs = 1) :
   if mpi_procs == 1 :
      cmd = []
   else :
      cmd = ["mpirun", "-np", "%d" % mpi_procs]
   cmd.append(flexi_path)
   cmd.append(prm_path)
   #print cmd
   p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
   lines = []
   while p.poll() is None :
      line = p.stdout.readlines()
      if line :
         lines.extend(line)
   if p.wait() != 0 :
       for line in lines :
           print line,
       print "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       print "!!     Flexi crashed!    !!"
       print "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       return None

   #for line in lines :
      #print line,
   if log :
      #log_path = os.path.splitext(os.path.basename(prm_path))[0] + ".log"
      log_path = projectname + ".log"
      f = open(log_path, 'a')
      for line in lines :
         f.write(line)
      f.close()

   if analyze_fcts :
      results = []
      if type(analyze_fcts) != list :
         return analyze_fcts(lines)
      for analyze_fct in analyze_fcts :
         results.append(analyze_fct(lines))

      return results
