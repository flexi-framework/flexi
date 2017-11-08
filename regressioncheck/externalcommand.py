import os
import subprocess
import logging
import tools
import select
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
        (pipeOut_r, pipeOut_w) = os.pipe()
        (pipeErr_r, pipeErr_w) = os.pipe()
        process = subprocess.Popen(cmd, stdout=pipeOut_w, \
                                        stderr=pipeErr_w, \
                                        universal_newlines=True, cwd=workingDir)

        self.stdout = []
        self.stderr = []

        bufOut = ""
        bufErr = ""
        while process.poll() is None:
            # Loop long as the selct mechanism indicates there
            # is data to be read from the buffer
            while len(select.select([pipeOut_r], [], [], 0)[0]) == 1:
                # Read up to a 1 KB chunk of data
                bufOut = bufOut + os.read(pipeOut_r, 1024)
                tmp = bufOut.split('\n') 
                for line in tmp[:-1] :
                    self.stdout.append(line+'\n')
                    log.info(line)
                bufOut = tmp[-1]
            while len(select.select([pipeErr_r], [], [], 0)[0]) == 1:
                # Read up to a 1 KB chunk of data
                bufErr = bufErr + os.read(pipeErr_r, 1024)
                tmp = bufErr.split('\n') 
                for line in tmp[:-1] :
                    self.stderr.append(line+'\n')
                    log.info(line)
                bufErr = tmp[-1]

        os.close(pipeOut_w)
        os.close(pipeOut_r)
        os.close(pipeErr_w)
        os.close(pipeErr_r)


        self.return_code = process.returncode

        end = timer()
        self.walltime = end - start

        # write std.out and err.out to disk
        self.stdout_filename = os.path.join(target_directory,name+".out")
        with open(self.stdout_filename, 'w') as f :
            for line in self.stdout :
                f.write(line)
        if self.return_code != 0 :
            self.result=tools.red("Failed")
            self.stderr_filename = os.path.join(target_directory,name+".err")
            with open(self.stderr_filename, 'w') as f :
                for line in self.stderr :
                    f.write(line)
        else :
            self.result=tools.blue("Successful")
        print self.result+" [%.2f sec]" % self.walltime

        return self.return_code
    
