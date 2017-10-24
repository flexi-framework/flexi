import argparse
import numpy as np
import os
import logging
import tools
import check
from timeit import default_timer as timer

print('='*132)
print "reggie2.0, add nice ASCII art here"
print('='*132)
                                
start = timer()

parser = argparse.ArgumentParser(description='Regression checker for NRG codes.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', '--carryon', action='store_true', help='''Continue build/run process. 
  --carryon         : build non-existing binary-combinations and run all examples for thoses builds
  --carryon --run   : run all failed examples''')
parser.add_argument('-e', '--exe', help='Path to executable of code that should be tested.')
parser.add_argument('-d', '--debug', type=int, default=0, help='Debug level.')
parser.add_argument('-j', '--buildprocs', type=int, default=0, help='Number of processors used for compiling (make -j XXX).')
parser.add_argument('-b', '--basedir', help='Path to basedir of code that should be tested (contains CMakeLists.txt).')
parser.add_argument('-y', '--dummy', action='store_true',help='use dummy_basedir and dummy_checks for fast testing on dummy code')
parser.add_argument('-r', '--run', action='store_true' ,help='run all binaries for all examples with all run-combinations for all existing binaries')
parser.add_argument('check', help='Path to check-/example-directory.')

args = parser.parse_args() # reggie command line arguments
cwd = os.getcwd()                                          # start with current working directory
found = os.path.exists(os.path.join(cwd,args.check)) # check if directory exists
if not found :
    print "Check directory not found: ",os.path.join(cwd,args.check)
    exit(1)

# setup logger for printing information, debug messages to stdout
tools.setup_logger(args.debug)
log = logging.getLogger('logger')

# setup basedir (search upward from staring point of reggie)
if args.dummy : # 
    args.check = 'dummy_checks/test'
    print "Check directoryswitched to ",args.check
    args.basedir = os.path.abspath('dummy_basedir')
    #print "basedir = ".ljust(25)+basedir
else :
    try :
        args.basedir = tools.find_basedir()
        #print "basedir = [".ljust(15)+basedir,"]"
    except Exception,ex :
        args.basedir = os.path.abspath('dummy_basedir')

# delete the building directory when [carryon = False] and [run = False] before getBuilds is called
if not args.carryon and not args.run : tools.clean_folder("reggie_outdir")

# get builds from checks directory if no executable is supplied
if args.exe is None : # if not exe is supplied, get builds
    # read build combinations from checks/XX/builds.ini
    builds = check.getBuilds(args.basedir, args.check)
else :
    found = os.path.exists(args.exe) # check if executable exists
    if not found :
        print tools.red("no executable found under ")
        exit(1)
    else :
        builds = [check.Standalone(args.exe,args.check)] # set builds list to contain only the supplied executable
        args.run = True
        args.basedir = None


if args.run :
    print "args.run -> skip building"
    # remove all build from builds when build.binary_exists() = False
    if builds[0].binary_exists() :
        builds = [build for build in builds]
        print builds
    else :
        print tools.red("no binary found under "+builds[0].binary_path)
        print tools.red(str(builds))
        exit(1)



# display all command line arguments
print "Running with the following command line options"
for arg in args.__dict__ :
    print arg.ljust(15)," = [",getattr(args,arg),"]"
print('='*132)




# perform build/run/analyze routines
check.PerformCheck(start,builds,args,log)




# print table with summary of errors if all builds where successful
check.SummaryOfErrors(start,builds)

