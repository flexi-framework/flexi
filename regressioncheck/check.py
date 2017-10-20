import os
import shutil
import combinations 
from loop import Loop
import tools
from timeit import default_timer as timer

class Build(Loop) :
    def __init__(self, basedir, source_directory,configuration, number, name='build') :
        self.basedir = basedir
        self.source_directory = source_directory
        self.configuration = configuration
        Loop.__init__(self, None, name, number)  
        
        # move 'binary' from 'configuration' dict to 'parameters' dict
        self.parameters = {'binary':self.configuration.get('binary','no binary supplied')}
        self.configuration.pop('binary', None) # remove binary from config dict

        # set path to binary/executable
        self.binary_path = os.path.abspath(self.target_directory+'/'+self.parameters['binary'])

        # set cmake command
        self.cmake_cmd = ["cmake"]                        # start composing cmake command
        for (key, value) in self.configuration.items() :  # add configuration to the cmake command
            self.cmake_cmd.append("-D%s=%s" % (key, value))    
        self.cmake_cmd.append(self.basedir)               # add basedir to the cmake command

    def compile(self, buildprocs) :
        # don't compile if build directory already exists
        if self.binary_exists() :  # if the binary exists, return
            print "skipping"
            return
        else : # for build carryon: when a binary is missing remove all examples (re-run all examples)
            print "removing folder, ",
            shutil.rmtree(self.target_directory,ignore_errors=True)
            os.makedirs(self.target_directory)
        print "building"

        # CMAKE
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

    def binary_exists(self) :
        return os.path.exists(self.binary_path)

class Standalone(Build) :
    def __init__(self,binary_path,source_directory) :
        Build.__init__(self, None, source_directory, {}, -1, "standalone")
        self.binary_path = binary_path

    def compile(self, buildprocs) :
        pass

    def __str__(self) :
        s = "standalone :       binary_path= " + self.binary_path + "\n"
        s+= "              target_directory= " + self.target_directory
        return s

def getBuilds(basedir, source_directory) :
    builds = []
    i = 0
    for b in combinations.getCombinations(os.path.join(source_directory, 'builds.ini')) :
        builds.append(Build(basedir, source_directory,b, i))
        i += 1
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
    # checks directory with 'builds.ini'
    if os.path.exists(os.path.join(build.source_directory, 'builds.ini')) :
        example_paths = [os.path.join(path,p) for p in sorted(os.listdir(path)) if os.path.isdir(os.path.join(path,p))]
    else :
        example_paths = [path]

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
class Run(Loop) :
    def __init__(self, parameters, path, command_line, number) :
        self.analyze_results = []
        self.analyze_successful = True
        self.parameters = parameters
        self.source_directory = os.path.dirname(path)
        Loop.__init__(self, command_line, 'run', number, mkdir=False)

        self.skip = os.path.exists(self.target_directory)
        if self.skip :
            return

        os.makedirs(self.target_directory)

        # copy all files in the source directory (example) to the target directory: always overwrite
        for f in os.listdir(self.source_directory) :
          src = os.path.abspath(os.path.join(self.source_directory,f))
          dst = os.path.abspath(os.path.join(self.target_directory,f))
          shutil.copyfile(src, dst)

    def rename_failed(self) :
        shutil.rmtree(self.target_directory+"_failed",ignore_errors=True)  # remove if exists
        shutil.move(self.target_directory,self.target_directory+"_failed") # rename folder (non-existent folder fails)
        self.target_directory = self.target_directory+"_failed" # set new name for summary of errors

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
        run = Run(r, path, command_line, i)
        if not run.skip :
            runs.append(run)
        i += 1
    return runs



