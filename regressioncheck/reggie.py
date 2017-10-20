import argparse
import numpy as np
import os
import logging
import tools
import check
import analyze
from timeit import default_timer as timer

print('='*132)
print "reggie2.0, add nice ASCII art here"
print('='*132)
                                
start = timer()
global_run_number=0
global_errors=0
build_number=0

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
        print tools.red("no binary found"+str(builds))
        exit(1)



# display all command line arguments
print "Running with the following command line options"
for arg in args.__dict__ :
    print arg.ljust(15)," = [",getattr(args,arg),"]"
print('='*132)




# General workflow:
# 1.   loop over alls builds
# 1.1    compile the build if args.run is false and the binary is non-existent
# 1.1    read all example directories in the check directory
# 2.   loop over all example directories
# 2.1    read the command line options in 'command_line.ini' for binary execution (e.g. number of threads for mpirun) 
# 2.2    read the analyze options in 'analyze.ini' within each example directory (e.g. L2 error analyze)
# 3.   loop over all command_line options
# 3.1    read the executable parameter file 'parameter.ini' (e.g. flexi.ini with which flexi will be started)
# 4.   loop over all parameter combinations supplied in the parameter file 'parameter.ini'
# 4.1    execute the binary file for one combination of parameters
# 5.   loop over all successfully executed binary results and perform analyze tests
# 6.   rename all run directories for which the analyze step has failed for at least one test

# compile and run loop
try : # if compiling fails -> go to exception

    # 1.   loop over alls builds
    for build in builds :
        build_number+=1 # count number of builds
        print "Build Cmake Configuration ",build_number," of ",len(builds)," ...",
        log.info(str(build))

        # 1.1    compile the build if args.run is false and the binary is non-existent
        build.compile(args.buildprocs)
        if not args.carryon :
            tools.clean_folder(os.path.join(build.target_directory,"examples"))
        
        # 1.1    read the example directories
        # get example folders: run_basic/example1, run_basic/example2 from check folder
        print args.check
        print build
        build.examples = check.getExamples(args.check, build)
        log.info("build.examples"+str(build.examples))

        # 2.   loop over all example directories
        for example in build.examples :
            log.info(str(example))
            
            # 2.1    read the command line options in 'command_line.ini' for binary execution 
            #        (e.g. number of threads for mpirun)
            example.command_lines = \
                    check.getCommand_Lines(os.path.join(example.source_directory,'command_line.ini'), example)
            
            # 2.2    read the analyze options in 'analyze.ini' within each example directory (e.g. L2 error analyze)
            example.analyzes = \
                    analyze.getAnalyzes(os.path.join(example.source_directory,'analyze.ini'), example)

            # 3.   loop over all command_line options
            for command_line in example.command_lines :
                log.info(str(command_line))

                # 3.1    read the executable parameter file 'parameter.ini' (e.g. flexi.ini with which 
                #        flexi will be started), N=, mesh=, etc.
                command_line.runs = \
                        check.getRuns(os.path.join(example.source_directory,'parameter.ini' ), command_line)

                # 4.   loop over all parameter combinations supplied in the parameter file 'parameter.ini'
                for run in command_line.runs :
                    log.info(str(run))

                    # 4.1    execute the binary file for one combination of parameters
                    run.execute(build,command_line)
                    global_run_number+=1
                    run.globalnumber=global_run_number
                    if not run.successful :
                        global_errors+=1

                # 5.   loop over all successfully executed binary results and perform analyze tests
                runs_successful = [run for run in command_line.runs if run.successful]
                if runs_successful : # do analyzes only if runs_successful is not emtpy
                    for analyze in example.analyzes :
                        #print tools.blue(">>>>>>>>>>>>>> ANALYZE <<<<<<<<<<<<<<<")
                        print tools.blue(str(analyze))
                        analyze.perform(runs_successful)

                # 6.   rename all run directories for which the analyze step has failed for at least one test
                for run in runs_successful : # all successful runs (failed runs are already renamed)
                    if not run.analyze_successful : # if 1 of N analyzes fails: rename
                        #print tools.blue(">>>>>>>>>>>>>> RENAME <<<<<<<<<<<<<<<")
                        #print run.target_directory
                        run.rename_failed()
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
d  = ' '
d2 = ' '
d3 = ' '
d4 = ' '
#invalid_keys = {"MPI", "binary", "analyze*"} # define keys to be removed from a dict
#parameters_removed = tools.without_keys(command_line.parameters, invalid_keys) # remove keys from dict

print "#run".center(5,d4)+"options".center(51,d4)+"path".center(65,d4)+"MPI".center(3,d4)+"time".rjust(5,d4)+"Information".rjust(12,d4)
for build in builds :
    #print('-'*132)
    print " "
    if type(build) is check.Build : print " ".join(build.cmake_cmd)
    for example in build.examples :
        for command_line in example.command_lines :
            #line=", ".join(["%s=%s"%item for item in command_line.parameters.items()])
            #print tools.yellow(tools.indent(line,4," "))
            for run in command_line.runs :
                line=", ".join(["%s=%s"%item for item in run.parameters.items()[1:]]) # skip first index
                if line != param_str_old : # only print when the parameter set changes
                    print tools.yellow(tools.indent(line,3))
                param_str_old=line
                line=str(run.globalnumber).rjust(4,d3)+" "*3 # global run number

                line+= tools.yellow("%s=%s"%(run.parameters.items()[0])) # only use first index
                line=line.ljust(65,d3) # inner most run variable (e.g. TimeDiscMethod)

                # build/example/reggie/run info
                line+=os.path.relpath(run.target_directory,"reggie_outdir").ljust(65,d2)

                line+=command_line.parameters.get('MPI','-').center(4,d3)
                #run.execution_time=21000.20
                #print "%2.1f".rjust(10,d2) % (run.execution_time)
                line+="%2.1f".rjust(5,d2) % (run.execution_time)
                line+=run.result.center(21,d3) # add result (successful or failed)
                print line
                for result in run.analyze_results :
                    print tools.red(result).rjust(148)
        print ""

if global_errors > 0 :
    tools.finalize(start,"Failed! Number of errors: "+str(global_errors))
else :
    tools.finalize(start,"successful")



