import argparse
import os
import logging

import logger
import compile 
import check

parser = argparse.ArgumentParser(description='Regression checker for NRG codes.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--mode', choices=['build', 'run'], default='build', help='''  --mode build : compile code for all binary-combinations and for all binaries run all examples with all run-combinations
  --mode run   : run all binaries for all examples with all run-combinations, compile missing binary-combinations.''')
parser.add_argument('-c', '--carryon', action='store_true', help='''Continue build/run process. 
  --carryon --mode build : build non-existing binary-combinations and run all examples for thoses builds
  --carryon --mode run   : run all failed examples''')
parser.add_argument('-e', '--exe', help='Path to executable of code that should be tested.')
parser.add_argument('-d', '--debug', type=int, default=0, help='Debug level.')
parser.add_argument('-j', '--buildprocs', type=int, default=0, help='Number of processors used for compiling (make -j XXX).')
parser.add_argument('-b', '--basedir', help='Path to basedir of code that should be tested (contains CMakeLists.txt).')
parser.add_argument('check', help='Path to check-/example-directory.')

args = parser.parse_args()

# setup logger for printing information, debug messages to stdout
logger.setup(args.debug)
log = logging.getLogger('logger')

basedir = os.path.abspath('dummy_basedir') #compile.find_basedir()

builds = check.getBuilds(basedir, os.path.join(args.check, 'builds.ini'))

for build in builds :
    log.info(str(build))
    build.compile(args.buildprocs)
    build.examples = check.getExamples(args.check, build)
    for example in build.examples :
        log.info(str(example))
        example.reggies = check.getReggies(os.path.join(example.path,'reggie.ini'), example) # MPI=1,2,3
        for reggie in example.reggies :
            log.info(str(reggie))
            reggie.runs = check.getRuns(os.path.join(example.path,'flexi.ini' ), reggie) # mesh= mesh1, mesh2 
            for run in reggie.runs :
                log.info(str(run))
                #run.execute()

print "=========================="

for build in builds :
    if not build.successful : 
        print "BUILD: failed", build.configuration
    for example in build.examples :
        if not example.successful : 
            print "  EXAMPLE: failed", example.path
        for reggie in example.reggies :
            if not reggie.successful : 
                print "    REGGIE: failed", reggie.parameters
            for run in reggie.runs :
                if not run.successful : 
                    print "      RUN: failed", run.parameters

