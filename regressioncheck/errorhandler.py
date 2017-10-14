import os
import subprocess
import logging

class ErrorHandler() :
    successful = True
    errors = []

    def __init__(self, parent, name, number = -1, mkdir=True) :
        self.number = number
        self.parent = parent
        if self.parent :
            parent_dir = self.parent.directory
        else :
            parent_dir = "reggie_outdir"
        if number >= 0 :
            self.directory = os.path.join(parent_dir, "%s_%04d" %(name, number))
        else :
            self.directory = os.path.join(parent_dir, name)

        if mkdir :
            if not os.path.exists(self.directory) :
                os.mkdir(self.directory)  # create example directory

    def execute_cmd(self, cmd):
        """Execute an external program specified by 'cmd'. The working directory of this program is set to self.directory.
        Returns return_code, stdout, stderr of the external program.
        """
        log = logging.getLogger('logger')

        workingDir = os.path.abspath(self.directory)
        log.debug(workingDir)
        log.debug(cmd)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, cwd=workingDir)
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

        return self.return_code, self.stdout, self.stderr

