#!/usr/bin/env python3

import sys
import os
import argparse

#============================================
# setupParser: Setup command line argument parser
#============================================
def setupParser():
	'''Setup the argument parser and return the parser object'''

	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--inFilePrefix',dest='vtmFile',
		                required=True,
		                help='prefix name of Solution .vtm file: e.g.,  testname from testname_Solution_0000000.000000000.vtm')
	parser.add_argument('-o','--out',dest='filename',
						default='vtm2pdv.pdv',
		                help='enter the output pvd file name')
	parser.add_argument('-t','--dt',dest='dt',
						required=True,
		                help='enter the time step between solution output files')
	parser.add_argument('-e','--extension',dest='ext',
		                default='vtm',
		                help='File extension: e.g.,  pvd from testname_Solution_0000000.000000000.pvd')
	return parser

#============================================
#                   MAIN
#============================================
if(__name__ == '__main__'):

	#--------------------------------
	# Handle arguments
	#--------------------------------
	parser = setupParser()
	args = parser.parse_args()

	# Extract parameters relevant to all problems
	inPrefix = args.vtmFile
	outExt = args.ext
	outFile = args.filename
	print("inPrefix: ", inPrefix)
	print("outFile: ", outFile)

	sys.stdout.write('  Writing to file: {:s}\n'.format(outFile))
	file = open(outFile, 'w+')

	file.write("<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n")
	file.write("  <Collection>\n")

	exvar = True
	time = 0.
	dt = float(args.dt) 

#==============================================
#   Setup input file name search
#==============================================

	while exvar:

		fname = '{:s}_Solution_{:0>17.9f}.{:s}'.format(inPrefix, time, outExt)

		if os.path.exists(fname):
			print("Writing file: {:s}\n".format(fname))

			# write file name into the pvd file 
			file.write("     <DataSet timestep=\"{:f}\" part=\"0\" file=\"{:s}\"/>\n".format( time, fname))
			time += dt;
		else:
			exvar = False

	file.write("  </Collection>\n</VTKFile>\n")
	file.close()

