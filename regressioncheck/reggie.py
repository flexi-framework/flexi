import argparse
import os

import logger
import combinations 
import compile 

parser = argparse.ArgumentParser(description='Regression checker for NRG codes.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--mode', choices=['build', 'run'], default='build', help='''  --mode build : compile code for all binary-combinations and for all binaries run all examples with all run-combinations
  --mode run   : run all binaries for all examples with all run-combinations, compile missing binary-combinations.''')
parser.add_argument('-c', '--carryon', action='store_true', help='''Continue build/run process. 
  --carryon --mode build : build non-existing binary-combinations and run all examples for thoses builds
  --carryon --mode run   : run all failed examples''')
parser.add_argument('-e', '--exe', help='Path to executable of code that should be tested.')
parser.add_argument('-d', '--debug', type=int, default=0, help='Debug level.')
parser.add_argument('-j', '--buildprocs', type=int, default=1, help='Number of processors used for compiling (make -j XXX).')
parser.add_argument('-b', '--basedir', help='Path to basedir of code that should be tested (contains CMakeLists.txt).')
parser.add_argument('check', help='Path to check-/example-directory.')

args = parser.parse_args()

# setup logger for printing information, debug messages to stdout
logger.setup(args.debug)

# remove following lines!!!
builds = combinations.getCombinations(os.path.join(args.check, 'builds.ini'))


for build in builds :
    print "BUILD:", build
    for example in ["freestream_3D"] :
        print " EXAMPLE:", example
        example_path = os.path.join(args.check, example)
        flexis   = combinations.getCombinations(os.path.join(example_path,'flexi.ini')) # mesh= mesh1, mesh2 
        reggies  = combinations.getCombinations(os.path.join(example_path,'reggie.ini')) # MPI=1,2,3
        excludes = combinations.getCombinations(os.path.join(example_path,'excludes.ini'))
        if combinations.anyIsSubset(excludes, build) : 
            print "Example is excluded"
            continue
        for reggie in reggies :
            print "  REGGIE:", reggie
            for flexi in flexis :
                print "   FLEXI:", flexi
                
                #for sample in samples : 



            #analyze(build, reggieconf)


#print type(c), c
#basedir = compile.find_basedir()
#compile.cmake('build', c[0], basedir)
