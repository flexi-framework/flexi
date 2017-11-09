import os
import shutil
import collections
import combinations 
from outputdirectory import OutputDirectory
from externalcommand import ExternalCommand
import tools
from analysis import Analyze, getAnalyzes
import collections

class Build(OutputDirectory,ExternalCommand) :

    def __init__(self, basedir, source_directory,configuration, number, name='build', binary_path=None) :
        self.basedir          = basedir
        self.source_directory = source_directory
        self.configuration    = configuration
        OutputDirectory.__init__(self, None, name, number)  
        ExternalCommand.__init__(self)

        # initialize result as empty list
        self.result = tools.yellow("skipped building")

        # initialize examples as empty list
        self.examples = []
        
        # set path to binary/executable
        if binary_path :
            self.binary_path = binary_path
        else :
            # get 'binary' from 'configuration' dict and remove it 
            try :
                binary_name = self.configuration["binary"]
            except :
                print tools.red("No 'binary'-option with the name of the binary specified in 'builds.ini'")
                exit(1)
            self.configuration.pop('binary', None) # remove binary from config dict
            self.binary_path = os.path.abspath(os.path.join(self.target_directory, binary_name))

        # set cmake command
        self.cmake_cmd = ["cmake"]                        # start composing cmake command
        self.cmake_cmd_color = ["cmake"]                  # start composing cmake command with colors
        for (key, value) in self.configuration.items() :  # add configuration to the cmake command
            self.cmake_cmd.append("-D%s=%s" % (key, value))    
            self.cmake_cmd_color.append(tools.blue("-D")+"%s=%s" % (key, value))    
        self.cmake_cmd.append(self.basedir)               # add basedir to the cmake command
        self.cmake_cmd_color.append(self.basedir)               # add basedir to the cmake command

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

        # CMAKE: execute cmd in build directory
        print "C-making with ["," ".join(self.cmake_cmd_color),"] ...",   # print colored string
        if self.execute_cmd(self.cmake_cmd, self.target_directory) != 0 : # use unclolored string for cmake
            raise BuildFailedException(self) # "CMAKE failed"

        # MAKE: default with '-j'
        self.make_cmd = ["make", "-j"]
        if buildprocs > 0 : self.make_cmd.append(str(buildprocs))
        # execute cmd in build directory
        print "Building with ["," ".join(self.make_cmd),"] ...",
        if self.execute_cmd(self.make_cmd, self.target_directory) != 0 :
            raise BuildFailedException(self) # "MAKE failed"
        print('-'*132)

    def __str__(self) :
        s = "BUILD in: " + self.target_directory
        return s

    def binary_exists(self) :
        return os.path.exists(self.binary_path)

class Standalone(Build) :
    def __init__(self,binary_path,source_directory) :
        Build.__init__(self, None, source_directory, {}, -1, "standalone", os.path.abspath(binary_path))

    def compile(self, buildprocs) :
        pass

    def __str__(self) :
        s = "standalone :       binary_path= " + self.binary_path + "\n"
        s+= "              target_directory= " + self.target_directory
        return s

def getBuilds(basedir, source_directory) :
    builds = []
    i = 1
    combis, digits = combinations.getCombinations(os.path.join(source_directory, 'builds.ini')) 
    for b in combis :
        builds.append(Build(basedir, source_directory,b, i))
        i += 1
    return builds

class BuildFailedException(Exception) :
    def __init__(self, build):
        self.build = build
    def __str__(self):
        return "build.compile failed in directory '%s'." % (self.build.target_directory)

#==================================================================================================

class Example(OutputDirectory) :
    def __init__(self, source_directory, build) :
        self.source_directory = source_directory
        OutputDirectory.__init__(self, build, os.path.join("examples",os.path.basename(self.source_directory)))

    def __str__(self) :
        s = "EXAMPLE in: " + self.source_directory
        return tools.indent(s,1)

def getExamples(path, build, log) :
    # checks directory with 'builds.ini'
    if os.path.exists(os.path.join(build.source_directory, 'builds.ini')) :
        example_paths = [os.path.join(path,p) for p in sorted(os.listdir(path)) \
                                              if os.path.isdir(os.path.join(path,p))]
    else :
        example_paths = [path]

    examples = [] # list of examples for each build
    # iterate over all example paths (directories of the examples)
    for p in example_paths :
        log.info('-'*132)
        log.info(tools.blue("example "+str(p)))
        # check if example should be excluded for the build.configuration
        exclude_path = os.path.join(p, 'excludeBuild.ini')
        if os.path.exists(exclude_path) :
            log.info(tools.blue("excludes under "+str(exclude_path)))
            # get all keys+values in 'excludeBuild.ini'
            options, _, _ = combinations.readKeyValueFile(exclude_path)
            # list of all excludes for comparison with 'build.configuration'
            excludes = [ { option.name : value } for option in options for value in option.values ]
            if combinations.anyIsSubset(excludes, build.configuration) :
                log.info(tools.red("  skipping example"))
                continue # any of the excludes matches the build.configuration. 
                         # Skip this example for the build.configuration
            else :
                log.info(tools.yellow("  not skipping"))
        examples.append(Example(p, build))
    return  examples


#==================================================================================================
class Command_Lines(OutputDirectory) :
    def __init__(self, parameters, example, number) :
        self.parameters = parameters
        OutputDirectory.__init__(self, example, 'cmd', number)
    def __str__(self) :
        s = "command_line parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return tools.indent(s,2)

def getCommand_Lines(path, example) :
    command_lines = []
    i = 1
    combis, digits = combinations.getCombinations(path) 
    for r in combis :
        command_lines.append(Command_Lines(r, example, i))
        i += 1
    return command_lines


#==================================================================================================
class Run(OutputDirectory, ExternalCommand) :
    total_errors = 0
    total_number_of_runs = 0

    def __init__(self, parameters, path, command_line, number, digits) :
        self.successful = True
        self.globalnumber = -1
        self.analyze_results = []
        self.analyze_successful = True
        self.parameters = parameters
        self.digits = digits
        self.source_directory = os.path.dirname(path)
        OutputDirectory.__init__(self, command_line, 'run', number, mkdir=False)
        ExternalCommand.__init__(self)

        self.skip = os.path.exists(self.target_directory)
        if self.skip :
            return

        os.makedirs(self.target_directory)

        # copy all files in the source directory (example) to the target directory: always overwrite
        for f in os.listdir(self.source_directory) :
          src = os.path.abspath(os.path.join(self.source_directory,f))
          dst = os.path.abspath(os.path.join(self.target_directory,f))
          if os.path.isdir(src) : # check if file or directory needs to be copied
              shutil.copytree(src, dst) # copy tree
          else :
              shutil.copyfile(src, dst) # copy file

    def rename_failed(self) :
        """Rename failed run directories in order to repeat the run when the regression check is repeated.
        This routine is called if either the execution fails or an analysis."""
        shutil.rmtree(self.target_directory+"_failed",ignore_errors=True)  # remove if exists
        shutil.move(self.target_directory,self.target_directory+"_failed") # rename folder (non-existent folder fails)
        self.target_directory = self.target_directory+"_failed" # set new name for summary of errors

    def execute(self, build, command_line) :
        Run.total_number_of_runs += 1
        self.globalnumber = Run.total_number_of_runs

        # set path to parameter file (single combination of values for execution "parameter.ini" for example)
        self.parameter_path = os.path.join(self.target_directory, "parameter.ini")

        # create parameter file with one set of combinations
        combinations.writeCombinationsToFile(self.parameters, self.parameter_path)

        # check MPI threads for mpirun
        MPIthreads = command_line.parameters.get('MPI')

        # check MPI built binary (only for reggie-compiled binaries possible)
        MPI_built_flag=os.path.basename(build.binary_path).upper()+"_MPI"
        MPIbuilt = build.configuration.get(MPI_built_flag,'ON')

        if MPIthreads :
            if MPIbuilt == "ON" :
                cmd = ["mpirun","-np",MPIthreads]
            else :
                print tools.indent(tools.yellow("Found %s=%s (binary has been built with MPI=OFF) with command_line setting MPIthreads=%s, running case in single (without 'mpirun -np')" % (MPI_built_flag,MPIbuilt,MPIthreads)),3)
                cmd = []
        else :
            cmd = []
        
        cmd.append(build.binary_path)
        cmd.append("parameter.ini")

        # append suffix commands, e.g., a second parameter file 'DSMC.ini' or '-N 12'
        cmd_suffix = command_line.parameters.get('cmd_suffix')
        if cmd_suffix :
            cmd.append(cmd_suffix)

        # append restart file name
        cmd_restart_file = command_line.parameters.get('restart_file')
        if cmd_restart_file :
            cmd.append(cmd_restart_file)

        # execute the command 'cmd'
        print tools.indent("Running [%s]" % (" ".join(cmd)), 2),
        self.execute_cmd(cmd, self.target_directory) # run the code

        if self.return_code != 0 :
            self.successful = False
            self.rename_failed()

    def __str__(self) :
        s = "RUN parameters:\n"
        s += ",".join(["%s: %s" % (k,v) for k,v in self.parameters.items()])    
        return tools.indent(s,3)

def getRuns(path, command_line) :
    """Get all combinations in 'parameter.ini'"""
    runs = []
    i = 1
    combis, digits = combinations.getCombinations(path)  # path to parameter.ini (source)
    for r in combis :
        run = Run(r, path, command_line, i, digits)
        if not run.skip :
            runs.append(run)
        i += 1
    return runs


def PerformCheck(start,builds,args,log) :
    """
    General workflow:
    1.   loop over alls builds
    1.1    compile the build if args.run is false and the binary is non-existent
    1.1    read all example directories in the check directory
    2.   loop over all example directories
    2.1    read the command line options in 'command_line.ini' for binary execution (e.g. number of threads for mpirun) 
    2.2    read the analyze options in 'analyze.ini' within each example directory (e.g. L2 error analyze)
    3.   loop over all command_line options
    3.1    read the executable parameter file 'parameter.ini' (e.g. flexi.ini with which flexi will be started)
    4.   loop over all parameter combinations supplied in the parameter file 'parameter.ini'
    4.1    execute the binary file for one combination of parameters
    5.   loop over all successfully executed binary results and perform analyze tests
    6.   rename all run directories for which the analyze step has failed for at least one test
    """

    build_number=0
    
    # compile and run loop
    try : # if compiling fails -> go to exception
    
        # 1.   loop over alls builds
        for build in builds :
            remove_build_when_successful = True
            build_number+=1 # count number of builds
            print "Build Cmake Configuration ",build_number," of ",len(builds)," ...",
            log.info(str(build))
    
            # 1.1    compile the build if args.run is false and the binary is non-existent
            build.compile(args.buildprocs)
            if not args.carryon : # remove examples folder if not carryon, in order to re-run all examples
                tools.remove_folder(os.path.join(build.target_directory,"examples"))
            
            # 1.1    read the example directories
            # get example folders: run_basic/example1, run_basic/example2 from check folder
            print build
            build.examples = getExamples(args.check, build,log)
            log.info("build.examples"+str(build.examples))
    
            if len(build.examples) == 0 :
                s = tools.yellow("No matching examples found for this build!")
                build.result += ", " + s
                print s
            # 2.   loop over all example directories
            for example in build.examples :
                log.info(str(example))
                print str(example)
                
                # 2.1    read the command line options in 'command_line.ini' for binary execution 
                #        (e.g. number of threads for mpirun)
                example.command_lines = \
                        getCommand_Lines(os.path.join(example.source_directory,'command_line.ini'), example)
                
                # 2.2    read the analyze options in 'analyze.ini' within each example directory (e.g. L2 error analyze)
                example.analyzes = \
                        getAnalyzes(os.path.join(example.source_directory,'analyze.ini'), example)
    
                # 3.   loop over all command_line options
                for command_line in example.command_lines :
                    log.info(str(command_line))
    
                    # 3.1    read the executable parameter file 'parameter.ini' (e.g. flexi.ini with which 
                    #        flexi will be started), N=, mesh=, etc.
                    command_line.runs = \
                            getRuns(os.path.join(example.source_directory,'parameter.ini' ), command_line)
    
                    # 4.   loop over all parameter combinations supplied in the parameter file 'parameter.ini'
                    for run in command_line.runs :
                        log.info(str(run))
    
                        # 4.1    execute the binary file for one combination of parameters
                        run.execute(build,command_line)
                        if not run.successful :
                            Run.total_errors+=1 # add error if run fails
    
                    # 5.   loop over all successfully executed binary results and perform analyze tests
                    runs_successful = [run for run in command_line.runs if run.successful]
                    if runs_successful : # do analyzes only if runs_successful is not emtpy
                        for analyze in example.analyzes :
                            print tools.indent(tools.blue(str(analyze)),2)
                            analyze.perform(runs_successful)
                    else : # don't delete build folder after all examples/runs
                        remove_build_when_successful = False
    
                    # 6.   rename all run directories for which the analyze step has failed for at least one test
                    for run in runs_successful :         # all successful runs (failed runs are already renamed)
                        if not run.analyze_successful :  # if 1 of N analyzes fails: rename
                            run.rename_failed()
                    if not any([run.analyze_successful for run in runs_successful]) : remove_build_when_successful = False # don't delete build folder after all examples/runs

            if remove_build_when_successful :
                tools.remove_folder(build.target_directory)
            print('='*132)

    # catch exception if bulding fails
    except BuildFailedException,ex:
        # print table with summary of errors
        SummaryOfErrors(builds)
    
        # display error message
        print ex # display error msg
        print tools.indent(" ".join(ex.build.cmake_cmd),1)
        print tools.indent(" ".join(ex.build.make_cmd),1)
        print tools.indent("Build failed, see: "+ex.build.stdout_filename,1)
        print tools.indent("                   "+ex.build.stderr_filename,1)
        print tools.bcolors.RED
        for line in ex.build.stderr[-20:] :
            print tools.indent(line,4),
        print tools.bcolors.ENDC

        print "run 'reggie' with the command line option '-c/--carryon' to skip successful builds."
        tools.finalize(start, 1, Run.total_errors, Analyze.total_errors)
        exit(1)


def SummaryOfErrors(builds) :
    """
    General workflow:
    1. loop over all builds, examples, command_lines, runs and for every run set the output strings 
       and get the maximal lengths of those strings
    2. print header 
    3. loop over all builds 
    3.1  print some information of the build
    3.2  within each build loop over all examples, command_lines, runs and for every run print some information:
    3.2.1  print an empty separation line if number of MPI threads changes
    3.2.2  print (only if changes) a line with all run parameters except the inner most, which is printed in 3.2.3
    3.2.3  print a line with following information:
             run.globalnumber, run.parameters[0] (the one not printed in 3.2.2), run.target_directory, MPI, run.walltime, run.result 
    3.2.4  print the analyze results line by line
    """
    
    param_str_old = ""
    str_MPI_old   = "-"

    # 1. loop over all runs and set output strings
    max_lens = collections.OrderedDict([ ("#run",4) , ("options",7) , ("path",4) , ("MPI",3), ("time",4) , ("Info",4) ])
    for build in builds :
        for example in build.examples :
            for command_line in example.command_lines :
                for run in command_line.runs :
                    run.output_strings = {}
                    run.output_strings['#run']    = str(run.globalnumber)
                    run.output_strings['options'] = ""
                    if run.digits.items()[0][1] > 0 :
                        run.output_strings['options'] += "%s=%s"%(run.parameters.items()[0])
                    run.output_strings['path']    = os.path.relpath(run.target_directory,OutputDirectory.output_dir)
                    run.output_strings['MPI']     = command_line.parameters.get('MPI', '-') 
                    run.output_strings['time']    = "%2.1f" % run.walltime
                    run.output_strings['Info']    = run.result
                    for key in run.output_strings.keys() :
                        max_lens[key] = max(max_lens[key], len(run.output_strings[key])) # set max column widths for summary table
    
    # 2. print header
    print 132*"="
    print " Summary of Errors"+"\n"
    spacing = 1
    for key, value in max_lens.items() :
        print key.ljust(value),spacing*' ',
    print ""
    
    # 3. loop over alls builds
    for build in builds :
    
        # 3.1 print cmake flags if no external binary was used for execution
        print('-'*132)
        if isinstance(build, Standalone) :
            print "Binary supplied externally under ",build.binary_path
        elif isinstance(build, Build) : 
            print "Build %d of %d (%s) compiled with in [%.2f sec]:" % (build.number, len(builds), build.result, build.walltime)
            print " ".join(build.cmake_cmd_color)
            if build.return_code != 0 : break # stop output as soon as a failed build in encountered
    
        # 3.2 loop over all examples, command_lines and runs
        for example in build.examples :
            for command_line in example.command_lines :
                for run in command_line.runs :
                    # 3.2.1 print separation line if MPI threads change
                    if run.output_strings["MPI"] != str_MPI_old :
                        print ""
                        str_MPI_old = run.output_strings["MPI"]
    
                    # 3.2.2 print the run parameters, execpt the inner most (this one is displayed in # 3.2.3)
                    paramsWithMultipleValues = [item for item in run.parameters.items()[1:] if run.digits[item[0]]>0 ]
                    param_str =", ".join(["%s=%s"%item for item in paramsWithMultipleValues]) # skip first index
                    if param_str  != param_str_old : # only print when the parameter set changes
                        print "".ljust(max_lens["#run"]), spacing*' ', tools.yellow(param_str)
                    param_str_old=param_str

                    # 3.2.3 print all output_strings
                    for key,value in max_lens.items() :
                        if key == "options" :
                            print tools.yellow(run.output_strings[key].ljust(value)),
                        else :
                            print run.output_strings[key].ljust(value),
                        print spacing*' ',
                    print ""

                    # 3.2.4  print the analyze results line by line
                    for result in run.analyze_results :
                        print tools.red(result).rjust(150)
    


