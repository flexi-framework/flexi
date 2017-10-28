import sys
import os
from combinations import getCombinations, splitValues

buildCombinations = getCombinations('examples/build_check/configurations.reggie')
runCombinations   = getCombinations('tmp.ini') 

iBuild = 0
for build in buildCombinations :
    iBuild +=1
    print 30*'='
    print "Build %d/%d:" % (iBuild, len(buildCombinations))
    print 30*'='

    cmake = "cmake"
    for name in build.keys() :
        cmake += " -D"+name+"="+build[name]
    cmake += " .."

    print cmake

    iRun = 0
    for run in runCombinations :
        iRun += 1
        print 30*'-'
        print "Run %d/%d:" % (iRun, len(runCombinations))
        print 30*'-'

        prm = ""
        for name in run.keys() :
            prm += name + "=" + run[name] + '\n'

        print prm
