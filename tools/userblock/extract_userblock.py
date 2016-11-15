#!/usr/bin/python
#************************************************************************************
#
# Author:       Matthias Sonntag
# Institution:  Inst. of Aero- and Gasdynamics, University of Stuttgart
# Date:         14.07.2016
#
# Description:  This script contains routines to extract userblock data from an HDF5
#               file containing git revision information, eventuelly a patch to
#               the Git remote and the ini (parameter) file used for the computation.
#
#************************************************************************************

import argparse
import subprocess

# extract userblock from HDF5 state file
def get_userblock(filename,userblock) :
    linesread = 0
    HDFfound = False
    f = open(filename, 'r')
    while True :
      line = f.readline()
      linesread = linesread + 1

      if line == "" :
        break

      if linesread == 1 and not line.startswith('{[(') :
        print 'Error: HDF5 state file %s contains no userblock.' % filename
        exit(1)

      if line.startswith('{[( END USERBLOCK )]}') :
        break

      # check for HDF identifier (here the original HDF state file begins)
      for i in range(len(line)) :
          c = line[i]
          if ord(c) == 0 : continue
          if ord(c) == 137 : 
              if line[i+1:i+4] == 'HDF' : HDFfound = True
      if HDFfound :
        break

      # check for compressed data
      if line.startswith('{[( COMPRESSED )]}') :
        try:
          filenamec=f.readline()      #size in bytes
          filenametar=filenamec.strip() + ".tar.xz"
          filesize=int(f.readline()) #size in bytes
          userblock_compressed=f.read(filesize)
        except:
          print 'Error: Could not extract compressed data.'
          exit(1)
        fc = open(filenametar, 'w')
        fc.write(userblock_compressed)
        fc.close()
        try:
          p = subprocess.call("tar -xJf " + filenametar,shell=True)
        except:
          print 'Error while extracting userblock data.'
          exit(1)

        try:
          userblock = get_userblock(filenamec.strip(),userblock)
        except:
          print 'Error while reading compressed userblock data.'
          exit(1)
        continue

      # everything ok
      userblock = userblock + line

    return userblock

# show all parts of the userblock
def print_all_parts(userblock) :
    for line in userblock.split('\n') :
        # try if line contains a part identifier: {[( IDENTIFIER )]}
        try :
            if not line.startswith("{[(") : continue
            identifier = line.split("{[(")[1].split(")]}")[0] 
            if "END USERBLOCK" in identifier : break
            # print identifier
            print identifier
        except :
            continue

def get_part(userblock,part) :
    ret = "" 
    output = False
    for line in userblock.split('\n') :
        # try if line contains a part identifier: {[( IDENTIFIER )]}
        try :
            if line.startswith("{[(") :
                identifier = line.split("{[(")[1].split(")]}")[0] 
                # if identifier is found -> start output
                if part in identifier :
                    output = True
                    continue
                else :
                    # we found another identifier -> stop output
                    if output : break
                if "END USERBLOCK" in identifier : break
        except :
            continue
        if output : ret = ret + line + "\n"
    return ret


# output a specific part of the userblock
def print_part(userblock,part) :
    tmp = get_part(userblock,part)
    if tmp :
       print tmp,

# output whole userblock
def print_all(userblock) :
    for line in userblock.split('\n') :
        print line


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract information from userblock of a HDF5-state file.')
    parser.add_argument('-s', '--show', action='store_true', help='show all available parts of the userblock')
    parser.add_argument('-p', '--part', help='indicates which part of the userblock should be extracted')
    parser.add_argument('-a', '--all', help='extract complete userblock', action='store_true')
    parser.add_argument('filename', type=str, help='filename of hdf5 file')

    args = parser.parse_args()

    userblock = get_userblock(args.filename,'') 
    if args.show :
        print_all_parts(userblock)
    if args.part :
        print_part(userblock,args.part)
    if args.all :
       print_all(userblock)
