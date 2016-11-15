#!/usr/bin/python
# -*- coding: utf8 -*-

import argparse
import os
import subprocess
import sys
import re

parser = argparse.ArgumentParser(description='Create DefineParemeters subrountine for new prm-parser')
parser.add_argument('file', help='File to work on')

args = parser.parse_args()

type = {"STR" : "String", "INT": "Int", "REAL": "Real", "LOGICAL": "Logical"}

options = []
module = ""

lines = open(args.file, 'r').readlines()
p  = re.compile("GET(STR|INT|REAL|LOGICAL)\s*\(([^)]*)\)")
pa = re.compile("GET(STR|INT|REAL|LOGICAL)ARRAY\s*\(([^)]*)\)")
for l in lines : 
    if l.startswith("MODULE") :
        module = l.split()[1].split("MOD_")[1]
    if "DefineParameters" in l : 
        print "File already processed"
        exit(1)
    m = p.search(l)
    if m: # non-array option 
        func = m.group(1)
        tmp = m.group(2).split(',')
        name = tmp[0]
        if len(tmp) > 1 :
            proposal = tmp[1]
            options.append("CALL prms%%Create%sOption(%s, \"TODO\", %s)" % (type[func], name, proposal))
        else :
            options.append("CALL prms%%Create%sOption(%s, \"TODO\")" % (type[func], name))

        continue
    m = pa.search(l)
    if m: # array option
        func = m.group(1)
        tmp = m.group(2).split(',')
        name = tmp[0]
        no = tmp[1]
        if len(tmp) > 2 :
            proposal = tmp[2]
            options.append("CALL prms%%Create%sArrayOption(%s, \"TODO\", %s)" % (type[func], name, proposal))
        else :
            options.append("CALL prms%%Create%sArrayOption(%s, \"TODO\")" % (type[func], name))
        continue

if len(options) == 0 : 
    print "No parameters in this file"
    exit(0)

s = """SUBROUTINE DefineParameters%s()
!==================================================================================================================================
! Define parameters 
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Parameters ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================\n""" % (module)
s = s+ "CALL prms%%SetSection(\"%s\")\n" %(module)


for o in options :
    s = s + o + "\n"

s = s + "END SUBROUTINE DefineParameters%s\n" % (module)

f = open(args.file, 'w')

for l in lines : 
    if l.strip() == "CONTAINS" :
        f.write("PUBLIC::DefineParameters%s\n" % (module)) 
    f.write(l)
    if l.strip() == "CONTAINS" :
        f.write(s)

        
