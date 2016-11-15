#!/bin/bash
# Script that autogenerates a qsub runsript
# checks the given parameters and launches the computation
#
# Set parameters here or use run.sh 'runargs' numcores runtime
# e.g. run.sh './flexi parameter.ini restartxy.h5' 4096 4:00:00

JOBNAME=phill
PPN=24

LOCALDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
JOBSCRIPT=$JOBNAME.pbs


if [[ "$1" == '-test' ]] ; then
  QUEUE='#PBS -q test'
  shift
fi

# Check args
if [ $# -lt 3 ] ; then
  echo "usage: run.sh numnodes time executable param1 param2 ..."
  echo "e.g. run.sh 123 4:00:0 ./flexi parameter.ini restartxy.h5"
  exit 1
fi


if [ $# -ne 0 ] ; then
  # set cores and runtime, then assemble RUNSTRING 
  if ! [[ "$1" =~ ^[0-9]+$ ]] ; then
    echo "Number of cores needs to be numeric!"; exit 1
  fi
  if ! [[ "$2" =~ ^[0-9]+(([:][0-9]+)+)?$ ]] ; then
    echo "Time needs to be numeric, either just ss or hh:mm:ss!"; exit 1
  fi
  NUMNODES=$1
  RUNTIME=$2
  EXECUTABLE=$3
  shift
  shift
  shift

  if [ ! -x $EXECUTABLE ] ; then
    echo "======"
    echo "== WARNING: Programm $EXECUTABLE is not an executable !"
    echo "== Make sure that the run command is correct !"
    echo "======"
  fi

  PARAMETERS=""
  while [[ $# > 0 ]] ; do
    if [ ! -e $1 ] ; then
      echo "======"
      echo "== WARNING: Parameter $1 is not a file in current directory !"
      echo "== Make sure that the run command is correct !               "
      echo "======"
    else
      if [[ $1 == *".ini" ]] ; then
        PROJECT=$(grep -i '^\s*projectname' $1 | rev | cut -d'=' -f1 | sed -e 's/ //g' | rev)
        if [[ ! -z $PROJECT ]] ; then
          JOBNAME=$PROJECT
        fi
      fi
    fi
    PARAMETERS+=$1' '
    shift
  done
fi

if [ ! -n "$PARAMETERS" ]; then
  echo "======"
  echo "== WARNING: No parameters specified ! "
  echo "== Make sure that the run command is correct !               "
  echo "======"
fi

# override testqueue time to max time on hornet
if [ -n "$QUEUE" ]; then
  RUNTIME='00:25:00'
fi
NUMCORES=`echo "$NUMNODES*$PPN" | bc`

echo "+++ Starting job with parameters +++"
echo JOBNAME     =$JOBNAME
echo EXECUTABLE  =$EXECUTABLE
echo PARAMETERS  =$PARAMETERS
echo NUMNODES    =$NUMNODES
echo PPN         =$PPN
echo NUMCORES    =$NUMCORES
echo RUNTIME     =$RUNTIME
echo LOCALDIR    =$LOCALDIR

# backup old file
if [ -e $JOBSCRIPT ] ; then
  mv $JOBSCRIPT $JOBSCRIPT.bak
  echo "Old jobscript moved to $JOBSCRIPT.bak"
fi

# autogenerate file
# PBS -m a|b|e: send mail at abort, beginning, end of job
cat > $JOBSCRIPT << EOF
#!/bin/bash
# Set job name, total number of cores, cores per node and runtime
#PBS -N $JOBNAME
#PBS -l nodes=$NUMNODES:ppn=$PPN
#PBS -l walltime=$RUNTIME
#PBS -m ab
#PBS -M $USERMAIL
#PBS -o $JOBNAME.out
#PBS -e $JOBNAME.err
$QUEUE

# Change directory and launch parallel job on the allocated compute nodes
cd $LOCALDIR
aprun -n $NUMCORES -N $PPN $EXECUTABLE $PARAMETERS
EOF

chmod u+x $JOBSCRIPT

# get job
qsub $JOBSCRIPT
