import os
import shutil

from errorhandler import ErrorHandler
from execute_cmd  import execute_cmd
import combinations 

def indent(text, amount, ch=' '):
    padding = amount * 2 * ch
    return ''.join(padding+line for line in text.splitlines(True))

class Build(ErrorHandler) :
    def __init__(self, basedir, configuration, number) :
        self.basedir = basedir
        self.configuration = configuration
        ErrorHandler.__init__(self, None, 'build', number, mkdir=False)   

        self.cmake_cmd = ["cmake"]                        # start composing cmake command
        for (key, value) in self.configuration.items() :  # add configuration to the cmake command
            self.cmake_cmd.append("-D%s=%s" % (key, value))    
        self.cmake_cmd.append(self.basedir)               # add basedir to the cmake command

    def compile(self, buildprocs) :
        # skip compiling if build directory already exists
        if os.path.exists(self.directory) :
            return

        # CMAKE
        os.makedirs(self.directory)                 # create build directory
        # execute cmd in build directory
        self.cmake_return_code, self.cmake_stdout, self.cmake_stderr = execute_cmd(self.cmake_cmd, self.directory)
        if self.cmake_return_code != 0 :
            shutil.rmtree(self.directory)
            raise Exception("CMAKE failed")

        # MAKE
        cmd = ["make", "-j"]
        if buildprocs > 0 : cmd.append(str(buildprocs))

        # execute cmd in build directory
        self.make_return_code, self.make_stdout, self.make_stderr = execute_cmd(cmd, self.directory)
        if self.make_return_code != 0 :
            shutil.rmtree(self.directory)
            raise Exception("MAKE failed")

    def __str__(self) :
        s = "BUILD in: " + self.directory + "\n"
        s += " ".join(self.cmake_cmd)
        return s


def getBuilds(basedir, path) :
    builds = []
    i = 0
    for b in combinations.getCombinations(path) :
        builds.append(Build(basedir, b, i))
        i += 1
    return builds

#==================================================================================================

class Example(ErrorHandler) :
    def __init__(self, path, build) :
        self.path = path
        ErrorHandler.__init__(self, build, os.path.basename(self.path))   

    def __str__(self) :
        s = "EXAMPLE in: " + self.path
        return indent(s,1)

def getExamples(path, build) :
    example_paths = [os.path.join(path,p) for p in sorted(os.listdir(path)) if os.path.isdir(os.path.join(path,p))]
    examples = []
    # iterate over all example paths (directories of the examples)
    for p in example_paths :
        # check if example should be excluded for the build.configuration
        exlcude_path = os.path.join(p, 'excludeBuild.ini')
        if os.path.exists(exlcude_path) :
            excludes = combinations.getCombinations(exlcude_path) 
            if combinations.anyIsSubset(excludes, build.configuration) :
                continue # any of the excludes matches the build.configuration. Skip this example for the build.configuration

        # append example to the return list
        examples.append(Example(p, build))
    return  examples


#==================================================================================================

class Reggie(ErrorHandler) :
    def __init__(self, parameters, example, number) :
        self.parameters = parameters
        ErrorHandler.__init__(self, example, 'reggie', number)


    def __str__(self) :
        s = "REGGIE parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return indent(s,2)


def getReggies(path, example) :
    reggies = []
    i = 0
    for r in combinations.getCombinations(path) :
        reggies.append(Reggie(r, example, i))
        i += 1
    return reggies




#==================================================================================================

class Run(ErrorHandler) :
    def __init__(self, parameters, reggie, number) :
        self.parameters = parameters
        ErrorHandler.__init__(self, reggie, 'run', number)

    def __str__(self) :
        s = "RUN parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return indent(s,3)

def getRuns(path, reggie) :
    runs = []
    i = 0
    for r in combinations.getCombinations(path) :
        runs.append(Run(r, reggie, i))
        i += 1
    return runs



