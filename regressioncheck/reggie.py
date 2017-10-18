import argparse
import numpy as np
import os
import logging
import tools
import check
from timeit import default_timer as timer
#import ast
#import re
import analyze_functions
                                
start = timer()
global_run_number=0
global_errors=0
build_number=0

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

# get builds from checks directory
if args.exe is None : # if not exe is supplied, get builds
    builds = check.getBuilds(args.basedir, os.path.join(args.check, 'builds.ini')) # read build combinations from checks/XX/builds.ini
else :
    found = os.path.exists(args.exe) # check if executable exists
    if not found :
        print tools.bcolors.FAIL+"no executable found under ",args.exe+tools.bcolors.ENDC
        exit(1)
    else :
        builds = [args.exe]
        args.mode = 'run'
        args.basedir = None






print "Running with the following command line options"
for arg in args.__dict__ :
    print arg.ljust(15)," = [",getattr(args,arg),"]"
print('='*132)




# delete the building directory when [carryon = False] and [mode = 'build']
if not args.carryon and args.mode=='build' : tools.clean_folder("reggie_outdir")

try : # if compiling fails -> go to exception
    for build in builds :
        build_number+=1
        print "Build Cmake Configuration ",build_number," of ",len(builds)," ...",
        log.info(str(build))
        
        build.compile(args.buildprocs)
        if not args.carryon and args.mode=='run' :
            tools.clean_folder(os.path.join(build.target_directory,"examples"))
        # get examples: run_basic/example1, run_basic/example2
        build.examples = check.getExamples(args.check, build)
        for example in build.examples :
            log.info(str(example))
            # get MPI=1,2,3
            example.command_lines = \
                    check.getCommand_Lines(os.path.join(example.source_directory,'command_line.ini'), example)
            # get analyze: L2, convtest, line integral
            example.analyzes = \
                    check.getAnalyzes(os.path.join(example.source_directory,'analyze.ini'), example)
            for command_line in example.command_lines :
                log.info(str(command_line))
                command_line.runs = \
                        check.getRuns(os.path.join(example.source_directory,'parameter.ini' ), command_line)
                for run in command_line.runs :
                    log.info(str(run))
                    if run.skip :
                        continue
                    run.execute(build,command_line)
                    global_run_number+=1
                    run.globalnumber=global_run_number

                    #run.analyze_results = []
                    if run.successful : # only do analysis if the run was successful
                        for analyze in example.analyzes :
                            log.info(str(analyze))
                            # L2 error check
                            L2_tolerance = float(analyze.parameters.get('analyze_L2',-1.))
                            if L2_tolerance > 0 :
                                L2_errors = np.array(analyze_functions.get_last_L2_error(run.stdout))
                                if (L2_errors > L2_tolerance).any() :
                                    run.analyze_results.append("analysis failed: L2 error >"+str(L2_tolerance))
                                    global_errors+=1
                                    analyze.successful=False
                    else :
                        global_errors+=1
        print('='*132)
except check.BuildFailedException,ex:
    print tools.bcolors.WARNING+""
    print ex # display error msg
    print tools.indent(" ".join(build.cmake_cmd),1)
    print tools.indent(" ".join(ex.build.make_cmd),1)
    print tools.indent("Build failed, see: "+ex.build.stdout_filename,1)
    print tools.indent("                   "+ex.build.stderr_filename,1)+tools.bcolors.ENDC
    print tools.bcolors.FAIL
    for line in ex.build.stderr[-20:] :
        print tools.indent(line,4),
    print tools.bcolors.ENDC
    tools.finalize(start,"FAILED!")
    exit(1)









print('='*132)
param_str_old=""
print " Summary of Errors"+"\n"
#print tools.indent("Compile flags",2)+" build/example/reggie/run".rjust(25)
d = ' '
d2 = '.'
#invalid_keys = {"MPI", "binary", "analyze*"}
#parameters_removed = tools.without_keys(command_line.parameters, invalid_keys)

print "#run".center(8,d),"options".center(37,d2),"path".center(44,d),"MPI".center(9,d2),"runtime".rjust(10,d),"Information".rjust(15,d2)
for build in builds :
    #print('-'*132)
    print " "
    print " ".join(build.cmake_cmd)
    for example in build.examples :
        for command_line in example.command_lines :
            #line=", ".join(["%s=%s"%item for item in command_line.parameters.items()])
            #print tools.yellow(tools.indent(line,4," "))
            for run in command_line.runs :
                if run.skip :
                    continue
                line=", ".join(["%s=%s"%item for item in run.parameters.items()[1:]]) # skip first index
                if line != param_str_old : # only print when the parameter set changes
                    print tools.yellow(tools.indent(line,5))
                param_str_old=line
                line=str(run.globalnumber).center(9,d)+" "*3 # global run number

                line+= tools.yellow("%s=%s"%(run.parameters.items()[0])) # only use first index
                line=line.ljust(55,d) # inner most run variable (e.g. TimeDiscMethod)

                # build/example/reggie/run info
                line+=os.path.relpath(run.target_directory,"reggie_outdir").ljust(25,d2)

                line+=command_line.parameters.get('MPI','-').center(9,d)
                line+="%2.2f".rjust(12,d2) % (run.execution_time)
                line+=run.result.rjust(25,d) # add result (successful or failed)
                print line
                for result in run.analyze_results :
                    print tools.red(result).rjust(137)
        print ""

if global_errors > 0 :
    tools.finalize(start,"Failed! Number of errors: "+str(global_errors))
else :
    tools.finalize(start,"successful")



