import os
import subprocess
import logging
import tools
from timeit import default_timer as timer

class ExternalCommand() :
    def __init__(self) :
        self.stdout = []
        self.stderr = []
        self.stdout_filename = None
        self.stderr_filename = None
        self.return_code = 0
        self.result = ""
        self.walltime = 0

    def execute_cmd(self, cmd, target_directory, name="std"):
        """Execute an external program specified by 'cmd'. The working directory of this program is set to target_directory.
        Returns the return_code of the external program.
        """
        log = logging.getLogger('logger')

        workingDir = os.path.abspath(target_directory)
        log.debug(workingDir)
        log.debug(cmd)
        start = timer()
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, \
                                        stderr=subprocess.PIPE, \
                                        universal_newlines=True, cwd=workingDir)
        self.stdout = []
        self.stderr = []
        for line in iter(process.stdout.readline, '') :
            log.info(line.rstrip())
            self.stdout.append(line)
        process.stdout.close()
        
        for line in iter(process.stderr.readline, '') :
            self.stderr.append(line)
        process.stderr.close()

        self.return_code = process.wait()

        end = timer()
        self.walltime = end - start

        # write std.out and err.out to disk
        self.stdout_filename = os.path.join(target_directory,name+".out")
        f = open(self.stdout_filename, 'w')
        for line in self.stdout :
            f.write(line)
        f.close()
        if self.return_code != 0 :
            self.result=tools.red("Failed")
            self.stderr_filename = os.path.join(target_directory,name+".err")
            f = open(self.stderr_filename, 'w')
            for line in self.stderr :
                f.write(line)
            f.close()
        else :
            self.result=tools.blue("Successful")
        print self.result+" [%.2f sec]" % self.walltime

        return self.return_code
    
