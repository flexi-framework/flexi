import os
import shutil
import combinations 
from loop import Loop
import tools
from timeit import default_timer as timer
import analysis

class Build(Loop) :
    def __init__(self, basedir, source_directory,configuration, number, name='build') :
        self.basedir          = basedir
        self.source_directory = source_directory
        self.configuration    = configuration
        Loop.__init__(self, None, name, number)  

        # initialize result as empty list
        self.result = tools.yellow("not built")

        # initialize examples as empty list
        #self.examples = []
        
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

        # CMAKE: execute cmd in build directory
        print "C-making with ["," ".join(self.cmake_cmd),"] ...",
        self.execute_cmd(self.cmake_cmd)
        if self.return_code != 0 :
            raise BuildFailedException(self) # "CMAKE failed"

        # MAKE: default with '-j'
        self.make_cmd = ["make", "-j"]
        if buildprocs > 0 : self.make_cmd.append(str(buildprocs))
        # execute cmd in build directory
        print "Building with ["," ".join(self.make_cmd),"] ...",
        self.execute_cmd(self.make_cmd)
        if self.return_code != 0 :
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
        self.binary_path = os.path.abspath(binary_path)

    def compile(self, buildprocs) :
        pass

    def __str__(self) :
        s = "standalone :       binary_path= " + self.binary_path + "\n"
        s+= "              target_directory= " + self.target_directory
        return s

def getBuilds(basedir, source_directory) :
    builds = []
    i = 1
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
    i = 1
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
    i = 1
    for r in combinations.getCombinations(path) : # path to parameter.ini (source)
        run = Run(r, path, command_line, i)
        if not run.skip :
            runs.append(run)
        i += 1
    return runs


def PerformCheck(start,builds,args,log) :
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
    build_number=0
    global_run_number=0
    
    # compile and run loop
    try : # if compiling fails -> go to exception
    
        # 1.   loop over alls builds
        for build in builds :
            build_number+=1 # count number of builds
            print "Build Cmake Configuration ",build_number," of ",len(builds)," ...",
            log.info(str(build))
    
            # 1.1    compile the build if args.run is false and the binary is non-existent
            build.compile(args.buildprocs)
            if not args.carryon : # remove examples folder if not carryon, in order to re-run all examples
                tools.clean_folder(os.path.join(build.target_directory,"examples"))
            
            # 1.1    read the example directories
            # get example folders: run_basic/example1, run_basic/example2 from check folder
            print args.check
            print build
            build.examples = getExamples(args.check, build)
            log.info("build.examples"+str(build.examples))
    
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
                        analysis.getAnalyzes(os.path.join(example.source_directory,'analyze.ini'), example)
    
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
                        global_run_number+=1
                        run.globalnumber=global_run_number
                        if not run.successful :
                            build.total_errors+=1 # add error if run fails
    
                    # 5.   loop over all successfully executed binary results and perform analyze tests
                    runs_successful = [run for run in command_line.runs if run.successful]
                    if runs_successful : # do analyzes only if runs_successful is not emtpy
                        for analyze in example.analyzes :
                            print tools.blue(str(analyze))
                            analyze.perform(runs_successful)
                    # add errors after analyze
                    build.total_errors+=sum(run.total_errors for run in command_line.runs if run.successful)
    
                    # 6.   rename all run directories for which the analyze step has failed for at least one test
                    for run in runs_successful :         # all successful runs (failed runs are already renamed)
                        if not run.analyze_successful :  # if 1 of N analyzes fails: rename
                            run.rename_failed()
            print('='*132)

    # catch exception if bulding fails
    except BuildFailedException,ex:
        # print table with summary of errors
        SummaryOfErrors(-1.,builds)
    
        # display error message
        print tools.bcolors.WARNING+"" # activate yellow text color
        print ex # display error msg
        print tools.indent(" ".join(ex.build.cmake_cmd),1)
        print tools.indent(" ".join(ex.build.make_cmd),1)
        print tools.indent("Build failed, see: "+ex.build.stdout_filename,1)
        print tools.indent("                   "+ex.build.stderr_filename,1)+tools.bcolors.ENDC # de-activate yellow
        print tools.bcolors.FAIL
        for line in ex.build.stderr[-20:] :
            print tools.indent(line,4),
        print tools.bcolors.ENDC
        global_errors = sum([build.total_errors for build in builds]) # sum up all errors from running and analyzing
        tools.finalize(start,"FAILED!",min(1,global_errors))
        exit(1)


def SummaryOfErrors(start,builds) :
    # General workflow:
    # 1.   loop over alls builds
    # 1.1    display cmake flags if no external binary was used for execution
    # 2.   loop over all examples
    # 3.   loop over all command lines
    # 4.   loop over all runs
    # 4.1.   display the run parameters, execpt the inner most (this one is displayed in # 4.2.2)
    # 4.2.   display one line of interformation for each run
    # 4.2.1    start line with the global run number
    # 4.2.2    add the inner most run variable (this one was skipped in # 4.1)
    # 4.2.3    add the path to directory where the run is executed (build/example/reggie/run info)
    # 4.2.4    add the number of MPI threads that where used for execution
    # 4.2.5    add the time that was needed for the execution of the binary and print the line
    # 4.3.     display the analyze results line by line
    # 5.   display the number of errors encountered during build/execution/analyze
    
    print('='*132)
    param_str_old = ""
    str_MPI_old   = "-"
    print " Summary of Errors"+"\n"
    d  = ' '
    d1 = ' '
    
    print "#run".center(5,d1)+"options".center(51,d1)+"path".center(65,d1)+"MPI".center(3,d1)+"time".rjust(8,d1)+"Information".rjust(12,d1),"\n"
    
    # 1.   loop over alls builds
    for build in builds :
    
        # 1.1    display cmake flags if no external binary was used for execution
        print('-'*132)
        if build.__class__ is Build : 
            print "Build  ",build.number," of ",len(builds),"compiled with ","\n"," ".join(build.cmake_cmd),build.result
            if build.result == tools.red("Failed") : break # stop output as soon as a failed build in encountered
        elif build.__class__ is Standalone :
            print "Binary supplied externally under ",build.binary_path
    
        # 2.   loop over all examples
        for example in build.examples :
            print ""
            
            # 3.   loop over all command lines
            for command_line in example.command_lines :
                str_MPI = command_line.parameters.get('MPI','-')
                if str_MPI != str_MPI_old : print "" # print separation line if MPI threads change
                str_MPI_old = str_MPI
    
                # 4.   loop over all runs
                for run in command_line.runs :
    
                    # 4.1    display the run parameters, execpt the inner most (this one is displayed in # 4.2.2)
                    line=", ".join(["%s=%s"%item for item in run.parameters.items()[1:]]) # skip first index
                    if line != param_str_old : # only print when the parameter set changes
                        print tools.yellow(tools.indent(line,3))
                    param_str_old=line
    
                    # 4.1.1    start line with the global run number
                    line=str(run.globalnumber).rjust(4,d)+" "*3
    
                    # 4.2.2    add the inner most run variable (this one was skipped in # 4.1)
                    line+= tools.yellow("%s=%s"%(run.parameters.items()[0]))
                    line=line.ljust(65,d)
    
                    # 4.2.3    add the path to directory where the run is executed (build/example/reggie/run info)
                    line+=os.path.relpath(run.target_directory,"reggie_outdir").ljust(65,d)
    
                    # 4.2.4    add the number of MPI threads that where used for execution
                    line+=str_MPI.center(4,d)
    
                    # 4.2.5    add the time that was needed for the execution of the binary
                    line+="%7.1f" % (run.execution_time)
                    #line+="%10.1f".rjust(9,d) % (run.execution_time)
                    line+=run.result.center(21,d) # add result (successful or failed)
                    print line
    
                    # 4.3.     display the analyze results line by line
                    for result in run.analyze_results :
                        print tools.red(result).rjust(150)
                #print ""
            #print ""
    
    # 5.   display the number of errors encountered during build/execution/analyze
    global_errors = sum([build.total_errors for build in builds]) # sum up all errors from running and analyzing
    if global_errors > 0 :
        tools.finalize(start,"Failed! Number of errors: ",global_errors)
    else :
        tools.finalize(start,"successful")


