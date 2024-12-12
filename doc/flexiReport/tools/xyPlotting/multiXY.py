#!/usr/bin/python

''' Plots each variable on separate plot'''

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
import matplotlib.pyplot as plt

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

    L = ['dens','eint','Density','rho']
    DEN = ['dens','Density','rho']
    PRES = ['P','Pres', 'pres','PF']
    EINT = ['u','eint','eF']
    VEL = ['v','velx','V', 'uF', '|v_i|']
    CSP = ['cspd','Cspd','c','C','cspdF']

    group = [DEN, PRES, EINT, VEL, CSP]

    IDX = [[0,0],[1,0],[2,0],[0,1],[1,1]]

    fig, ax = plt.subplots(3, 2)

    for o, a in opts:
        if o in ("-h", "-help"):
            print( 'multiXY.py -i <inputfile> -t <tab separated data file> -o <outputfile>')
            sys.exit(0)
        elif o in ("-i", "-infile"):
            # For exactpack output files
            infile = a

            fig1 = plt.figure()
            ax = plt.subplot(111)
            varnameParts = a.split('.')
            var = varnameParts[0]
            print( 'Input file is: ', infile)

            data = np.loadtxt(infile, skiprows=1)

            ax.plot(data[:,0], data[:,1], label=var)
            plt.xlabel(varnameParts[0])  
            
        elif o in ("-t", "-tsv"):
            # Read Paraview exported tsv file; first line has var names
            # Also, toro_exact output file
            
            names, d = parseFile(a)

            for idx, gr in enumerate(group, start=0):

                I = IDX[idx][0]
                J = IDX[idx][1]
            
                for  n in gr:
                    #print("n: ", n)
                    if n in d:
                        #print(n, d[n],"\n")
                        ax[I,J].plot(d[names[0]], d[n], label=n)
                        ax[I,J].set_ylabel(n)

                if n in EINT:
                    #print("n: ", n)
                    ax[I][J].set_yscale('log')

                # Shrink current axis by 20%
                #box = ax[I,J].get_position()
                #ax[I,J].set_position([box.x0, box.y0, box.width * 0.8, box.height])

                # Put a legend to the right of the current axis
                #ax[I,J].legend(loc='center left', bbox_to_anchor=(1, 0.5))
                #ax[I,J].legend(loc='right')
            #ux = mx/rho
            #sie = E/rho - ux*ux/2
            #P = (5./3. - 1.)*(E - rho*ux*ux/2.)
            
            ax[2,1].set_axis_off()


           #plt.legend()
        elif o in ("-o", "-outfile"):
            outfile = a
            print( 'Output file is: ', outfile)

    # Shrink current axis by 20%
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #plt.legend()
    plt.autoscale('y',tight=True)    
    fig.tight_layout()
    if outfile != '':
        plt.savefig(outfile)  

    plt.show()

if __name__ == "__main__":
   main(sys.argv[1:])









