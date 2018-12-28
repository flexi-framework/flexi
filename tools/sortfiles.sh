#!/bin/bash

PROJECTNAME=INSERT_PROJECTNAME_HERE

if [ ! -n "$1" ]; then
  echo You have to specify a directory name!
  exit 1
fi
DIR=$1

if [ ! -d "$DIR" ]; then
  mkdir $DIR
fi
if [ ! -d "$DIR/state/" ]; then
  mkdir $DIR/state
fi
if [ ! -d "$DIR/timeavg/" ]; then
  mkdir $DIR/timeavg
fi
if [ ! -d "$DIR/baseflow/" ]; then
  mkdir $DIR/baseflow
fi
if [ ! -d "$DIR/rp/" ]; then
  mkdir $DIR/rp
fi
if [ ! -d "$DIR/logs/" ]; then
  mkdir $DIR/logs
fi

# remove flexi error files and system error and logfiles
rm -f "$PROJECTNAME"_ERRORS_[0-9]*.out
rm -f "$PROJECTNAME"_ERROR_State_[0-9]*.[0-9]*.h5
rm -f "$PROJECTNAME".e[0-9]* "$PROJECTNAME".o[0-9]*
rm -f "$PROJECTNAME".err

   STATE=$(ls "$PROJECTNAME"_State_[0-9]*.[0-9]*.h5    | sort -nr | tail -n +2)
BASEFLOW=$(ls "$PROJECTNAME"_BaseFlow_[0-9]*.[0-9]*.h5 | sort -nr | tail -n +2)
 TIMEAVG=$(ls "$PROJECTNAME"_TimeAvg_[0-9]*.[0-9]*.h5  | sort -nr | tail -n +2)
      RP=$(ls "$PROJECTNAME"_RP_[0-9]*.[0-9]*.h5       | sort -nr | tail -n +2)
# move all but last
for file in $STATE;    do mv --backup=simple $file $DIR/state/.    ; done
for file in $BASEFLOW; do mv --backup=simple $file $DIR/baseflow/.    ; done
for file in $TIMEAVG;  do mv --backup=simple $file $DIR/timeavg/. ; done
for file in $RP;       do mv --backup=simple $file $DIR/rp/.  ; done
# copy last file to keep last set
cp --backup=simple "$PROJECTNAME"_State_[0-9]*.[0-9]*.h5    $DIR/state/.   
cp --backup=simple "$PROJECTNAME"_BaseFlow_[0-9]*.[0-9]*.h5 $DIR/baseflow/.   
cp --backup=simple "$PROJECTNAME"_TimeAvg_[0-9]*.[0-9]*.h5  $DIR/timeavg/.
cp --backup=simple "$PROJECTNAME"_RP_[0-9]*.[0-9]*.h5       $DIR/rp/. 


# copy data to files to logdir
mv --backup=t Log.*.sdb $DIR/logs/.
mv --backup=t $PROJECTNAME*.out $DIR/logs/.
cp --backup=t $PROJECTNAME*.dat $DIR/logs/.
