from timeit import default_timer as timer
import logging
import tools
import check
import args_parser
"""
General workflow:
1.  get the command line arguments 'args' and all valid build combinations in the check directory from 'builds.ini'
2.  set the logger 'log' with the debug level from 'args' to determine the level of logging which displays output to the user
3.  perform the regression check by a) building executables
                                    b) running the code
                                    c) performing the defined analyzes
4.  display the summary table with information for each build, run and analysis step
5.  display if regression check was successful or not and return the corresponding error code
"""

print(132*'='+"\n"+"reggie2.0, add nice ASCII art here"+"\n"+132*'=')
start = timer()

# 1.  get the command line arguments 'args' and all valid build combinations in the check directory from 'builds.ini'
args, builds = args_parser.getArgsAndBuilds()

# 2.  set the logger 'log' with the debug level from 'args' to determine the level of logging which displays output to the user
tools.setup_logger(args.debug)
log = logging.getLogger('logger')

# 3.  perform the regression check by a) building executables
#                                     b) running the code
#                                     c) performing the defined analyzes
check.PerformCheck(start,builds,args,log)

# 4.  display the summary table with information for each build, run and analysis step
check.SummaryOfErrors(builds)

# 5.  display if regression check was successful or not and return the corresponding error code
tools.finalize(start, 0, check.Run.total_errors, check.Analyze.total_errors)
