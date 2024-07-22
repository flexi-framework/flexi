#!/usr/bin/python

''' Compute error norm between Flexi solution and exact solution
    error = L2(Flexi - exact)
	for each variable
	'''

def parseFile(filename):

  from numpy import loadtxt, genfromtxt

  with open(filename, 'r') as infile:
    with open(r'output.txt', 'w') as outfile:
      data = infile.read()
      data = data.replace(",", " ")
      outfile.write(data)
  infile.close()
  outfile.close()

  #with open(filename) as file:
  with open(r'output.txt') as file:

      A=np.genfromtxt(file, comments='#', skip_header=0, unpack=True)

      file.seek(0, 0)

      NM = np.genfromtxt(file, comments='#', dtype='str', unpack=True)
      names = NM[:,0]
  file.close()

  #print("names: ", names)

  na = [n.strip('"') for n in names[:]]
  #print("names: \n", na)

  d = {k: v for k, v in zip(na, A[:,1:])}
  #print("d{}:", d,"\n")

  return na, d

import sys, getopt
import numpy as np
from numpy.linalg import norm

def main(argv):

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:t:h", ["help", "infile=","outfile=","tsv="])

        if len(opts) == 0:
            print ('usage: multiXYplot.py -i <first xy plot file name> -i <second xy file> ... -t <tab separated data file> -o <plot file name>')
            sys.exit(2)

    except getopt.GetoptError as err:
        # print help information and exit:
        sys.stderr.write("%s: %s\n" % (sys.argv[0], "Unknown option"))
        sys.stderr.write("Usage: `%s -help' for more information\n" % sys.argv[0])
        sys.exit(2)

    infile = ''
    outfile = ''
    var = ''

    RHOE = ['EnergyStagnationDensity','rhoE','Erho']
    MOMX = ['MomentumX','momx', 'mx','mom']
    DEN = ['dens','Density','rho']
    PRES = ['P','pres','Pres', 'PF']
    EINT = ['e','eint','eF']
    VEL = ['u','velx','uF','V']
    #CSP = ['cspd','Cspd','c','C','cspdF']

    #group = [DEN, PRES, EINT, VEL]
    group = [DEN, RHOE, MOMX, PRES, EINT, VEL]

    den = []
    pres = []
    eint = []
    vel = []
    rhoe = []
    momx = []

    #VARS = [den,pres,eint,vel]
    VARS = [den,rhoe,momx,pres,eint,vel]
    X = []
    vname = []
    ID = []

    for o, a in opts:
        if o in ("-h", "-help"):
            print( 'multiXY.py -i <inputfile> -t <tab separated data file> -o <outputfile>')
            sys.exit(0)
        elif o in ("-i", "-infile"):
            # For exactpack output files
            '''            infile = a
            fig1 = plt.figure()
            ax = plt.subplot(111)
            varnameParts = a.split('.')
            var = varnameParts[0]
            print( 'Input file is: ', infile)

            data = np.loadtxt(infile, skiprows=1)

            ax.plot(data[:,0], data[:,1], label=var)
            plt.xlabel(varnameParts[0])
            '''
        elif o in ("-t", "-tsv"):
            # Read Paraview exported tsv file; first line has var names
            # Also, toro_exact output file

            names, d = parseFile(a)

            for idx, gr in enumerate(group, start=0):
                for  n in gr:
                    #print("n: ", n)
                    if n in d:
                        ID.append(idx)
                        #print("idx: ", idx, " d[", n, "]: ", d[n],"\n")
                        #print("x: ", d[names[0]])
                        #print("len(d[", n, "]: ", len(d[n]))
                        #print("VARS: ", VARS)
                        #print("n: ", n)
                        VARS[idx].append(d[n])
                        vname.append(n)
                        #ax[I,J].plot(d[names[0]], d[n], label=n)
        elif o in ("-o", "-outfile"):
            outfile = a
            print( 'Output file is: ', outfile)

    #print("VARS: ", VARS)
    #print("vname: ", vname)

    if outfile != '':
        with open(outfile, 'w') as f:
            i = 0
            for var in VARS:
                if(var == []):
                   continue
                #print(var)
                #print("vname[", i, "]: ", vname[i])
                #print("var[0]: ", var[0])
                #print("var[1]: ", var[1])
                #print("max(|var[0]-var[1]|): ", np.max(np.abs(var[0]-var[1])))
                L0 = len(var[0])
                L1 = len(var[1])
                if(L0 == L1):
                   #L2norm = np.sqrt(norm(var[0] - var[1]))/np.sqrt(norm(var[1]))
                   L2norm = np.sqrt(norm(var[0] - var[1]))/L0
                   print(f"{vname[i]},  L2 norm: {L2norm}")
                   print(vname[i], L2norm, file=f)
                   i += 1
                elif(L0 != L1):
                   print("ERROR: Input files have different lengths!\n")
                   exit()

if __name__ == "__main__":
   main(sys.argv[1:])








