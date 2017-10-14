import os
import shutil

from errorhandler import ErrorHandler
from execute_cmd  import execute_cmd
import combinations 

def indent(text, amount, ch=' '):
    padding = amount * 2 * ch
    return ''.join(padding+line for line in text.splitlines(True))

class Build(ErrorHandler) :
    def __init__(self, basedir, number, configuration) :
        self.basedir = basedir
        self.number = number
        self.configuration = configuration
        self.build_directory = os.path.join('reggie_outdir', 'build_%04d' % self.number)

        self.cmake_cmd = ["cmake"]                        # start composing cmake command
        for (key, value) in self.configuration.items() :  # add configuration to the cmake command
            self.cmake_cmd.append("-D%s=%s" % (key, value))    
        self.cmake_cmd.append(self.basedir)               # add basedir to the cmake command

    def compile(self, buildprocs) :
        # skip compiling if build_directory already exists
        if os.path.exists(self.build_directory) :
            return

        # CMAKE
        os.makedirs(self.build_directory)                 # create build directory
        # execute cmd in build_directory
        self.cmake_return_code, self.cmake_stdout, self.cmake_stderr = execute_cmd(self.cmake_cmd, self.build_directory)
        if self.cmake_return_code != 0 :
            shutil.rmtree(self.build_directory)
            raise Exception("CMAKE failed")

        # MAKE
        cmd = ["make", "-j"]
        if buildprocs > 0 : cmd.append(str(buildprocs))

        # execute cmd in build_directory
        self.make_return_code, self.make_stdout, self.make_stderr = execute_cmd(cmd, self.build_directory)
        if self.make_return_code != 0 :
            shutil.rmtree(self.build_directory)
            raise Exception("MAKE failed")

    def __str__(self) :
        s = "BUILD in: " + self.build_directory + "\n"
        s += " ".join(self.cmake_cmd)
        return s


def getBuilds(basedir, path) :
    builds = []
    i = 0
    for b in combinations.getCombinations(path) :
        builds.append(Build(basedir, i, b))
        i += 1
    return builds

#==================================================================================================

class Example(ErrorHandler) :
    def __init__(self, path) :
        self.path = path

    def __str__(self) :
        s = "EXAMPLE in: " + self.path
        return indent(s,1)

def getExamples(path, build_configuration) :
    example_paths = [os.path.join(path,p) for p in sorted(os.listdir(path)) if os.path.isdir(os.path.join(path,p))]
    examples = []
    # iterate over all example paths (directories of the examples)
    for p in example_paths :
        # check if example should be excluded for the build_configuration
        exlcude_path = os.path.join(p, 'excludeBuild.ini')
        if os.path.exists(exlcude_path) :
            excludes = combinations.getCombinations(exlcude_path) 
            if combinations.anyIsSubset(excludes, build_configuration) :
                continue # any of the excludes matches the build_configuration. Skip this example for the build_configuration

        # append example to the return list
        examples.append(Example(p))
    return  examples


#==================================================================================================

class Reggie(ErrorHandler) :
    def __init__(self, parameters) :
        self.parameters = parameters

    def __str__(self) :
        s = "REGGIE parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return indent(s,2)

def getReggies(path) :
    return [Reggie(r) for r in combinations.getCombinations(path)]




#==================================================================================================

class Run(ErrorHandler) :
    def __init__(self, parameters) :
        self.parameters = parameters

    def __str__(self) :
        s = "RUN parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return indent(s,3)

def getRuns(path) :
    return [Run(r) for r in combinations.getCombinations(path)]



