import subprocess
import logging


def execute_cmd(cmd, workingDir):
    """Execute an external program specified by 'cmd'. The working directory of this program is set with 'workingDir'.
    Returns return_code, stdout, stderr of the external program.
    """
    log = logging.getLogger('logger')

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, cwd=workingDir)
    stdout = []
    stderr = []
    for line in iter(process.stdout.readline, '') :
        log.info(line.rstrip())
        stdout.append(line)
    process.stdout.close()
    
    for line in iter(process.stderr.readline, '') :
        stderr.append(line)
    process.stderr.close()

    return_code = process.wait()

    return return_code, stdout, stderr

