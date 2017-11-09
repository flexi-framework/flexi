import logging
import shutil
import os
from timeit import default_timer as timer
import re
import tools

class bcolors :
    """color and font style definitions for changing output appearance"""
    # Reset (user after applying a color to return to normal coloring)
    ENDC   ='\033[0m'    

    # Regular Colors
    BLACK  ='\033[0;30m' 
    RED    ='\033[0;31m' 
    GREEN  ='\033[0;32m' 
    YELLOW ='\033[0;33m' 
    BLUE   ='\033[0;34m' 
    PURPLE ='\033[0;35m' 
    CYAN   ='\033[0;36m' 
    WHITE  ='\033[0;37m' 

    # Text Style
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def red(text) :
    return bcolors.RED+text+bcolors.ENDC

def green(text) :
    return bcolors.GREEN+text+bcolors.ENDC

def blue(text) :
    return bcolors.BLUE+text+bcolors.ENDC

def yellow(text) :
    return bcolors.YELLOW+text+bcolors.ENDC

def indent(text, amount, ch=' '):
    """Indent text line by amount times a white space """
    padding = amount * 2 * ch
    return ''.join(padding+line for line in text.splitlines(True))

def setup_logger(debug_level):
    """Setups a global logger with the name 'logger'. 
    This logger can accessed in any function by "log = logging.getLogger('logger')".
    Three different logging levels:
        0 : print no logging messages
        1 : print information messages (i.e. print all messages invoked with "log.info(message)")
        2 : print debug + information messages (i.e. print all messages invoked with "log.info(message)" or "log.debug(message)")
    """ 

    if debug_level == 0   : # no logging
        formatter = logging.Formatter()  
    elif debug_level == 1 : # info 
        formatter = logging.Formatter(fmt='%(message)s')
    elif debug_level == 2 : # debug
        formatter = logging.Formatter(fmt='%(levelname)s - %(module)s: %(message)s')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger('logger')
    if debug_level == 0 :   # no logging
        logger.setLevel(0)
    elif debug_level == 1 : # info
        logger.setLevel(logging.INFO)
    elif debug_level == 2 : # debug
        logger.setLevel(logging.DEBUG)

    logger.addHandler(handler)
    return logger

def find_basedir() :
    """Search 'CMakeLists.txt' in directories above current working directory.
    The directory containing the 'CMakeLists.txt' is the 'basedir'."""
    basedir = os.getcwd()                                           # start with current working directory
    found = os.path.exists(os.path.join(basedir, "CMakeLists.txt")) # check if actual directory is the basedir
    while not found :                                               # look upwards until basedir found
        basedir = os.path.dirname(basedir)                              # basedir = basedir/..
        found = os.path.exists(os.path.join(basedir, "CMakeLists.txt")) # check if actual directory is the basedir
        if basedir == "/" : break                                       # check if root of filesystem is reached

    if not found :
        raise Exception("No basedir found. Started searching for 'CMakeLists.txt' in '%s'" % os.getcwd())

    return basedir


def remove_folder(path) :
    print tools.yellow("[remove_folder]: deleting folder '%s'" % path)
    shutil.rmtree(path,ignore_errors=True)
    #shutil.rmtree(path)

def finalize(start, build_errors, run_errors, analyze_errors) :
    """Display if regression check was successful or not and return the corresponding error code"""
    if build_errors + run_errors + analyze_errors > 0 :
        print bcolors.RED + 132*'='
        print "reggie 2.0  FAILED!",
        return_code = 1
    else :
        print bcolors.BLUE + 132*'='
        print "reggie 2.0  successful!",
        return_code = 0

    if start > 0 : # only calculate run time and display output when start > 0
        end = timer()
        print "in [%2.2f sec]" % (end - start)
    else :
        print ""

    print "Number of build   errors: %d" % build_errors
    print "Number of run     errors: %d" % run_errors
    print "Number of analyze errors: %d" % analyze_errors

    print '='*132 + bcolors.ENDC
    exit(return_code)

def diff_lists(x,x_ref,tol,tol_type) :
    """
    determine diff of two lists of floats, either relative of absolute 
    (if the reference value is zero, use absolute comparison)
    x        : vector of real values
    x_ref    : vector of real values (reference)
    tol      : tolerance value
    tol_type : tolerance type, relative or absolute
    """

    # check tolerance type: absolute/relative (is the reference value is zero, absolute comparison is used)
    if tol_type == 'absolute' :
        diff = [abs(a-b) for (a,b) in zip(x,x_ref)]
    else : # relative comparison
        # if the reference value is zero, use absolute comparison
        diff = [abs(a/b-1.0) if abs(b) > 0.0 else a for (a,b) in zip(x,x_ref) ]

    # determie success logical list for return variable
    success = [d <= tol for d in diff]

    # display information when a diff is not successful, display value+reference+difference
    if not all(success) :
        print "Differences in vector comparison:"
        print "%13s   %13s   %13s" % ("x","x_ref","diff")
        for i in range(len(diff)) :
            if diff[i] <= tol : continue
            print "%13.6e   %13.6e   %13.6e" % (x[i],x_ref[i],diff[i])

    return success

