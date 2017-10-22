import logging
import shutil
import os
from timeit import default_timer as timer
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def indent(text, amount, ch=' '):
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


def clean_folder(path) :
    print "clean_folder: deleting folder ",path
    shutil.rmtree(path,ignore_errors=True)
    #shutil.rmtree(path)

def red(text) :
    return bcolors.FAIL+text+bcolors.ENDC

def green(text) :
    return bcolors.OKGREEN+text+bcolors.ENDC

def blue(text) :
    return bcolors.OKBLUE+text+bcolors.ENDC

def yellow(text) :
    return bcolors.WARNING+text+bcolors.ENDC

def finalize(start,text,global_errors=0) :
    text+=str(global_errors)
    if global_errors > 0:
        print bcolors.FAIL+""
    else :
        print bcolors.OKBLUE+""
    print('='*132)
    if start > 0 : # only calculate run time and display output when start > 0
        end = timer()
        print "reggie2.0 ",text," [%2.2f sec]" % (end - start)
        print('='*132)
    print ""+bcolors.ENDC


#invalid_keys = {"MPI", "binary", "analyze*"} # define keys to be removed from a dict
#parameters_removed = tools.without_keys(command_line.parameters, invalid_keys) # remove keys from dict

#def without_keys(d, keys) :
#    # remove keys from a dict and return a dict
#    return {x: d[x] for x in d if x not in keys}


