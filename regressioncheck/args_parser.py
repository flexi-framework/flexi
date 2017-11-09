import argparse
import os
import tools
import check
from outputdirectory import OutputDirectory

def getArgsAndBuilds() :
    """get command line arguments and builds in check directory from 'builds.ini'"""
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
    
    # get reggie command line arguments
    args = parser.parse_args()
    
    # setup basedir
    if args.dummy : 
        # For testing reggie during reggie-developement: 
        # Overwrite basedir and check directory with dummy directories.
        reggieDir = os.path.dirname(os.path.realpath(__file__))
        args.basedir = os.path.join(reggieDir, 'dummy_basedir')
        args.check =   os.path.join(reggieDir, 'dummy_checks/test')
        print "Basedir directory switched to '%s'" % args.basedir
        print "Check   directory switched to '%s'" % args.check
    else :
        # For real reggie-execution:
        # Setup basedir (containing CMakeLists.txt) by searching upward from current working directory 
        try :
            args.basedir = tools.find_basedir()
        except Exception,ex :
            print tools.red("Basedir (containing 'CMakeLists.txt') not found!\nEither specify the basedir on the command line or execute reggie within a project with a 'CMakeLists.txt'.")
            exit(1)
    
        if not os.path.exists(args.check) : # check if directory exists
            print tools.red("Check directory not found: '%s'" % args.check)
            exit(1)
    
    
    # delete the building directory when [carryon = False] and [run = False] before getBuilds is called
    if not args.carryon and not args.run : tools.remove_folder(OutputDirectory.output_dir)
    
    # get builds from checks directory if no executable is supplied
    if args.exe is None : # if not exe is supplied, get builds
        # read build combinations from checks/XX/builds.ini
        builds = check.getBuilds(args.basedir, args.check)
    else :
        if not os.path.exists(args.exe) : # check if executable exists
            print tools.red("No executable found under '%s'" % args.exe)
            exit(1)
        else :
            builds = [check.Standalone(args.exe,args.check)] # set builds list to contain only the supplied executable
            args.run = True      # set 'run-mode' do not compile the code
            args.basedir = None  # since code will not be compiled, the basedir is not needed
    
    if args.run :
        print "args.run -> skip building"
        # in 'run-mode' remove all build from list of builds if their binaries do not exist (build.binary_exists() == False)
        builds = [build for build in builds if build.binary_exists()]
    
    if len(builds) == 0 :
        print tools.red("List of 'builds' is empty! Maybe switch off '--run'.")
        exit(1)
    
    # display all command line arguments
    print "Running with the following command line options"
    for arg in args.__dict__ :
        print arg.ljust(15)," = [",getattr(args,arg),"]"
    print('='*132)
    
    
    return args, builds
