#!/usr/bin/python

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

     #print("data: ", data)
     #A=np.genfromtxt(data, comments='#', skip_header=1, unpack=True)
      A=np.genfromtxt(file, comments='#', skip_header=1, unpack=True)
      X = A[0,1:]
      file.seek(0, 0)
     #NM = np.loadtxt(data, comments='#', dtype='str', unpack=True)
      NM = np.genfromtxt(file, comments='#', dtype='str', unpack=True)
      names = NM[:,0]
  file.close()

  #print("names: ", names)
  #print("X: ", X)
  #print("data[:,1:]: ", A[:,1:])

  return names, X, A

import sys, getopt
import numpy as np
import matplotlib.pyplot as plt

def main(argv):

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:t:h", ["help", "infile=","outfile=","tsv="])

        if len(opts) == 0:
            print ('usage: simpleXYplot.py -i <first xy plot file name> -i <second xy file> ... -t <tab separated data file> -o <plot file name>')
            sys.exit(2)
        
    except getopt.GetoptError as err:
        # print help information and exit:
        sys.stderr.write("%s: %s\n" % (sys.argv[0], "Unknown option"))
        sys.stderr.write("Usage: `%s -help' for more information\n" % sys.argv[0])
        sys.exit(2)

    infile = ''
    outfile = ''
    var = ''

    Exact = []

    fig = plt.figure()
    ax = plt.subplot(111)

    for o, a in opts:
        if o in ("-h", "-help"):
            print( 'simpleXYplot.py -i <inputfile> -t <tab separated data file> -o <outputfile>')
            sys.exit(0)
        elif o in ("-i", "-infile"):
            infile = a
            varnameParts = a.split('.')
            var = varnameParts[0]
            print( 'Input file is: ', infile)

            data = np.loadtxt(infile, skiprows=1)

            ax.plot(data[:,0], data[:,1], label=var)
            plt.xlabel(varnameParts[0])  
            
        elif o in ("-t", "-tsv"):
            # Read Paraview exported tsv file; first line has var names
            
            names, X, data = parseFile(a)
            
            #print("X: ", X)
            #print("data: ", data.T)
           
            ctr = 1
            for col in data[ctr:,1:]:
                #print("ctr: ", ctr)
                #print("col: ", col)
                #print("name: ", names[ctr])

                ax.plot(X, col, label=names[ctr])
                ctr += 1
            plt.xlabel(names[0])  

            #plt.legend(loc='best')
            #ux = mx/rho
            #sie = E/rho - ux*ux/2
            #P = (5./3. - 1.)*(E - rho*ux*ux/2.)
            
        elif o in ("-o", "-outfile"):
            outfile = a
            print( 'Output file is: ', outfile)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #plt.legend()
    
    if outfile != '':
        plt.savefig(outfile)  

    plt.show()
    
if __name__ == "__main__":
   main(sys.argv[1:])









