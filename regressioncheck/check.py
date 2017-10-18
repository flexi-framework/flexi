import os
import shutil
import combinations 
from loop import Loop
import tools
from timeit import default_timer as timer

class Build(Loop) :
    def __init__(self, basedir, configuration, number) :
        self.basedir = basedir
        self.configuration = configuration
        Loop.__init__(self, None, 'build', number, mkdir=False)   
        
        # move 'binary' from 'configuration' dict to 'parameters' dict
        self.parameters = {'binary':self.configuration.get('binary','no binary supplied')}
        self.configuration.pop('binary', None) # remove binary from config dict

        # set path to binary/executable
        self.binary_path = os.path.abspath(self.target_directory+'/'+self.parameters['binary'])
        self.skip = os.path.exists(self.binary_path)

        # set cmake command
        self.cmake_cmd = ["cmake"]                        # start composing cmake command
        for (key, value) in self.configuration.items() :  # add configuration to the cmake command
            self.cmake_cmd.append("-D%s=%s" % (key, value))    
        self.cmake_cmd.append(self.basedir)               # add basedir to the cmake command

    def compile(self, buildprocs) :
        # skip compiling if build directory already exists
        if self.skip :  # if the binary exists, return
            print "skipping"
            return
        else :
            print "removing folder, ",
            shutil.rmtree(self.target_directory,ignore_errors=True)
        print "building"

        # CMAKE
        os.makedirs(self.target_directory)          # create build directory
        # execute cmd in build directory
        print "C-making with ["," ".join(self.cmake_cmd),"] ...",
        self.execute_cmd(self.cmake_cmd)
        if self.return_code != 0 :
            #shutil.rmtree(self.target_directory)
            raise BuildFailedException(self) # "CMAKE failed"

        # MAKE
        self.make_cmd = ["make", "-j"]
        if buildprocs > 0 : self.make_cmd.append(str(buildprocs))
        # execute cmd in build directory
        print "Building with ["," ".join(self.make_cmd),"] ...",
        self.execute_cmd(self.make_cmd)
        if self.return_code != 0 :
            #shutil.rmtree(self.target_directory) # remove reggie_outdir/build_0000
            raise BuildFailedException(self) # "MAKE failed"
        print('-'*132)

    def __str__(self) :
        s = "BUILD in: " + self.target_directory + "\n"
        s += " ".join(self.cmake_cmd)
        return s


def getBuilds(basedir, path) :
    builds = []
    i = 0
    for b in combinations.getCombinations(path) :
        builds.append(Build(basedir, b, i))
        i += 1
    print "Total number of valid builds: ",i
    return builds

class BuildFailedException(Exception) :
    def __init__(self, build):
        self.build = build
    def __str__(self):
        return "build.compile failed in directory '%s'." % (self.build.target_directory)

#==================================================================================================

class Example(Loop) :
    def __init__(self, source_directory, build) :
        self.source_directory = source_directory
        Loop.__init__(self, build, os.path.join("examples",os.path.basename(self.source_directory)))

    def __str__(self) :
        s = "EXAMPLE in: " + self.source_directory
        return tools.indent(s,1)

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
class Command_Lines(Loop) :
    def __init__(self, parameters, example, number) :
        self.parameters = parameters
        Loop.__init__(self, example, 'command_line', number)
    def __str__(self) :
        s = "command_line parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return tools.indent(s,2)

def getCommand_Lines(path, example) :
    #print path
    command_lines = []
    i = 0
    for r in combinations.getCombinations(path) :
        command_lines.append(Command_Lines(r, example, i))
        i += 1
    return command_lines


#==================================================================================================
class Analyze(Loop) :
    def __init__(self, parameters, example, number) :
        self.parameters = parameters
        Loop.__init__(self, example, 'analyze', number, mkdir=False)
    def __str__(self) :
        s = "analyze parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return tools.indent(s,2)

def getAnalyzes(path, example) :
    #print path
    analyze = []
    i = 0
    for r in combinations.getCombinations(path) :
        analyze.append(Analyze(r, example, i))
        i += 1
    return analyze

#==================================================================================================
class Run(Loop) :
    def __init__(self, parameters, path, command_line, number) :
        self.analyze_results = []
        self.parameters = parameters
        self.source_directory = os.path.dirname(path)
        Loop.__init__(self, command_line, 'run', number)
        if self.skip :
            return
        for f in os.listdir(self.source_directory) : # copy all files in the source directory (example) to the target directory
          src = os.path.abspath(os.path.join(self.source_directory,f))
          dst = os.path.abspath(os.path.join(self.target_directory,f))
          shutil.copyfile(src, dst)

    def rename_failed(self) :
        shutil.rmtree(self.target_directory+"_failed",ignore_errors=True)  # remove if exists
        shutil.move(self.target_directory,self.target_directory+"_failed") # rename folder

    def execute(self, build, command_line) :

        # set path to parameter file (single combination of values for execution "parameter.ini" for example)
        self.parameter_path = os.path.join(self.target_directory, "parameter.ini")

        # create parameter file with one set of combinations
        combinations.writeCombinationsToFile(self.parameters, self.parameter_path)

        # check MPI threads for mpirun
        MPIthreads = command_line.parameters.get('MPI')
        if MPIthreads :
            cmd = ["mpirun","-np",MPIthreads]
        else :
            cmd = []
        
        cmd.append(build.binary_path)
        cmd.append("parameter.ini")

        # append suffix commands, e.g., a second parameter file 'DSMC.ini' or '-N 12'
        cmd_suffix = command_line.parameters.get('cmd_suffix')
        if cmd_suffix :
            cmd.append(cmd_suffix)

        # execute the command 'cmd'
        start = timer()
        print "Running ["," ".join(cmd),"]",
        self.execute_cmd(cmd) # run the code
        end = timer()
        self.execution_time = end - start

        if self.return_code != 0 :
            self.successful = False
            self.rename_failed()

    def __str__(self) :
        s = "RUN parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return tools.indent(s,3)

def getRuns(path, command_line) :
    runs = []
    i = 0
    for r in combinations.getCombinations(path) : # path to parameter.ini (source)
        runs.append(Run(r, path, command_line, i))
        i += 1
    return runs



