import argparse
import os
import logging
import tools
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
parser.add_argument('-y', '--dummy', action='store_true',help='use dummy_basedir and dummy_checks for fast testing on dummy code')
parser.add_argument('check', help='Path to check-/example-directory.')

args = parser.parse_args() # reggie command line arguments

print('='*132)
print "reggie2.0, add nice ASCII art here"
print('='*132)
print "Running with the following command line options"
#print "args=",args
for arg in args.__dict__ :
    print arg.ljust(25)," = [",getattr(args,arg),"]"
print('='*132)
cwd = os.getcwd()                                          # start with current working directory
found = os.path.exists(os.path.join(cwd,args.check)) # check if directory exists
if not found :
    print "Check directory not found: ",os.path.join(cwd,args.check)
    exit(1)

# setup logger for printing information, debug messages to stdout
tools.setup_logger(args.debug)
log = logging.getLogger('logger')

# setup basedir (search upward from staring point of reggie)
if args.dummy :
    args.check = 'dummy_checks/test'
    print "Check directoryswitched to ",args.check
    basedir = os.path.abspath('dummy_basedir')
    print "basedir = ".ljust(25),basedir
else :
    try :
        basedir = tools.find_basedir()
        print "basedir = [".ljust(15),basedir,"]"
    except Exception,ex :
        basedir = os.path.abspath('dummy_basedir')
        print "basedir = ".ljust(25),basedir

# delete the building directory
tools.clean_folder()

# get builds from checks directory
builds = check.getBuilds(basedir, os.path.join(args.check, 'builds.ini'))

print " "

try : # if compiling fails -> go to exception
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
                    run.execute(build,reggie)
except check.BuildFailedException,ex:
    print ex
    print ex.build.make_cmd
    print " Build failed, see: ",ex.build.stdout_filename
    print " Build failed, see: ",ex.build.stderr_filename
    for line in ex.build.stderr[-20:] :
        print line,
    exit(1)



print('='*132)

for build in builds :
    if not build.successful : 
        print "BUILD: failed", build.configuration
    else :
        print "BUILD: successful", build.configuration
    for example in build.examples :
        if not example.successful : 
            print "  EXAMPLE: failed", example.path
        for reggie in example.reggies :
            if not reggie.successful : 
                print "    REGGIE: failed", reggie.parameters
            for run in reggie.runs :
                if not run.successful : 
                    print "      RUN: failed", run.parameters

