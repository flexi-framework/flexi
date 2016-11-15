#!/bin/bash
#-------------------------------------------------------------------------------------------------------------------------------------
# User Input 
#-------------------------------------------------------------------------------------------------------------------------------------
# Queue name
SIM_NAME=phill
# Queue name
INI_FILE=parameter.ini
# Number of cores
NUM_CORES=4096
# Running time for each subjob
SIM_TIME='4:00:00'
# File to restart from (leave empty if no restart)
RESTART_FILE=phill_State_0000700.0000000.h5
# Gap time before the end of running job when subjob will be submitted
GAP_TIME='0:30:00'
# Number of subjobs (total time=SUBJOBSxSIM_TIME)
NUM_SUBJOBS=4
#-------------------------------------------------------------------------------------------------------------------------------------
# Do not modify anything below this line
#-------------------------------------------------------------------------------------------------------------------------------------

WORKDIR=`pwd`
echo $WORKDIR
cd $WORKDIR

#-- Get seconds --
IFS=":"
Time=($SIM_TIME)
SIMULATION_TIME=$((${Time[0]}*3600 + ${Time[1]}*60 + ${Time[2]}))
Time=($GAP_TIME)
WAIT_TIME=$((${SIMULATION_TIME}-(${Time[0]}*3600 + ${Time[1]}*60 + ${Time[2]})))
IFS=""

#-- Update qsub file
sed -i 's/#PBS -N .*/#PBS -N '$SIM_NAME'/g' xjob.qsub
sed -i 's/walltime=.*/walltime='$SIM_TIME'/g' xjob.qsub
sed -i 's/mppwidth=.*/mppwidth='$NUM_CORES'/g' xjob.qsub
sed -i 's/NUM_CORES=.*/NUM_CORES='$NUM_CORES'/g' xjob.qsub
sed -i 's/INI_FILE=.*/INI_FILE='$INI_FILE'/g' xjob.qsub
sed -i 's/RESTART_FILE=.*/RESTART_FILE='$RESTART_FILE'/g' xjob.qsub

echo "----------------------------------------------------------------------------------------------------------------------" >> xjobs.logfile
echo "Simulation start"
echo "----------------------------------------------------------------------------------------------------------------------" >> xjobs.logfile
echo "---->>>>" >> xjobs.logfile
echo "-->> First job / Simulation restart from restart file" >> xjobs.logfile
qsub xjob.qsub >job_id_first
echo $(cat job_id_first)
JOB_ID=$(cat job_id_first | awk -F. '{print $1}')
echo "---->>>>" >> xjobs.logfile
echo -e "----------------------------------------------------------------------------------------------------------------------\n\n" >> xjobs.logfile

echo "----------------------------------------------------------------------------------------------------------------------" >> xjobs.logfile
echo  "-->>Loop over $NUM_SUBJOBS subjobs" >>xjobs.logfile
echo "----------------------------------------------------------------------------------------------------------------------" >> xjobs.logfile
# ------------------------------------------------------------------------------------------------------------
# loop over all subjobs 
# ------------------------------------------------------------------------------------------------------------
echo $NUM_SUBJOBS
# Update Time 
UPDATE_TIME=60
for ((i=1; i<=$NUM_SUBJOBS; i++));  do

  #-----------------------------------------------------------------------------------------------------------
  # Get Information from job 
  #-----------------------------------------------------------------------------------------------------------
  if [ -f job_id ]; then
  JOB_ID=$(cat job_id | awk -F. '{print $1}')
   
  fi
  STATUS=$(showq | grep $JOB_ID | awk '{print $3}')

  #-----------------------------------------------------------------------------------------------------------
  # Check status of current job
  #-----------------------------------------------------------------------------------------------------------

  echo "---->>>>" >> xjobs.logfile
  until [ "$STATUS" == 'Running' ]; do
    CHECK_OUTPUT=$(echo "$SIM_NAME.e$JOB_ID")
    if [ -f $CHECK_OUTPUT ]; then
      if  grep -q "FLEXI FINISHED" Log.$JOB_ID.sdb || ! grep -q "walltime" $CHECK_OUTPUT  ; then
        echo "********************************************************************************"  >>xjobs.logfile
        echo "$(date  '+%b %e %T'): Error occured or Flexi finished!! Please check Error-Filei " >>xjobs.logfile
        echo "********************************************************************************"   >>xjobs.logfile
        exit
      else
        echo "-->> do nothing !!!" >>xjobs.logfile
      fi
    fi
    STATUS=$(showq | grep $JOB_ID | awk '{print $3}')
    echo "-->> $(date  '+%b %e %T'): Submitted Job ($JOB_ID) is $STATUS" >>xjobs.logfile
    echo "-->> $(date  '+%b %e %T'): Check status ... until submitted job is running" >>xjobs.logfile
    sleep $UPDATE_TIME
  done
  echo "-->> $(date  '+%b %e %T'): Job ($JOB_ID) is finally $STATUS !! " >>xjobs.logfile
  echo "---->>>>" >> xjobs.logfile
echo -e "----------------------------------------------------------------------------------------------------------------------\n\n" >> xjobs.logfile
  #-----------------------------------------------------------------------------------------------------------
  # Check status of current job
  #-----------------------------------------------------------------------------------------------------------
 
  if [ -f job_id_first ]; then
  mv job_id_first job_id  
  fi


  
  if [ "$STATUS" == 'Running' ]; then
    #---------------------------------------------------------------------------------------------------------
    # Start time 
    #---------------------------------------------------------------------------------------------------------
    
    START_TIME=$(showq | grep $JOB_ID | awk '{print $7,$8,$9;}' )
    echo "$START_TIME " >>xjobs.logfile
    LogStat=0
    
    #---------------------------------------------------------------------------------------------------------
    # Loop until a new job will be submit for this case
    #---------------------------------------------------------------------------------------------------------
    until [ "$LogStat" == "1" ]; do
      CURRENT_TIME=$(date  '+%b %e %T')

      #-------------------------------------------------------------------------------------------------------      
      # Calculate time difference between running job and current time [time in sec]
      #-------------------------------------------------------------------------------------------------------      
      diff () {
              printf '%s' $(( $(date  -d"$CURRENT_TIME" +%s) -
                              $(date  -d"$START_TIME" +%s)))
      }
      DIFF=$(diff)
      # Check the time
      LogStat=$(echo "$DIFF > $WAIT_TIME" | bc)

      SLEEP_TIME=`echo  $WAIT_TIME - $DIFF | bc`

      #-------------------------------------------------------------------------------------------------------      
      # Wait until gap time is reached 
      #-------------------------------------------------------------------------------------------------------      

      echo "------>>>>" >>xjobs.logfile
      if [ "$LogStat" == "0" ]; then
        if [ "$SLEEP_TIME" > "1" ]; then
        echo "-->> $(date  '+%b %e %T'): waiting ... until a subjob can be submitted" >>xjobs.logfile
        fi
        sleep $SLEEP_TIME
      fi
      echo -e "------>>>>\n\n" >>xjobs.logfile
    done
    #---------------------------------------------------------------------------------------------------------

    echo "------>>>>" >>xjobs.logfile
    echo "-->> $(date  '+%b %e %T'): Time gap is over -->> start a new subjob" >>xjobs.logfile
    
    #-------------------------------------------------------------------------------------------------------      
    # Before submitting new job  -  check whether a simulation error is occured 
    #                               or T_END is reached 
    #-------------------------------------------------------------------------------------------------------      
    LOGQUEU=0
    CHECK_OUTPUT=$(echo "$SIM_NAME.e$JOB_ID")
    if [ -f $CHECK_OUTPUT ]; then
      if grep -q "FLEXI FINISHED" Log.$JOB_ID.sdb || ! grep -q "walltime" $CHECK_OUTPUT  ; then
        echo "***************************************************************************************"  >>xjobs.logfile
        echo "$(date  '+%b %e %T'): Error occured or Flexi finished!! Please check Error-Filei" >>xjobs.logfile
        echo "***************************************************************************************"   >>xjobs.logfile
        exit
      else
        qsub xjob.qsub >job_id_new
        JOB_ID_NEW=$(cat job_id_new | awk -F. '{print $1}')
        echo "-->> $(date  '+%b %e %T'): New subjob( $JOB_ID_NEW) is submitted, after previous job finished!!!" >>xjobs.logfile
        until [ "$LOGQUEU" == "1" ]; do
           CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.e$JOB_ID_NEW")
           if showq | grep -q $JOB_ID_NEW || [ -f $CHECK_OUTPUT_NEW ]; then         
             LOGQUEU=1
             if [ -f $CHECK_OUTPUT_NEW ]; then
               STATUS_NEW="ABORTED"
               echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
             else
               STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
               echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
             fi
           fi
           sleep 5 
        done           
      LOGQUEU=0
      fi
    else
      qsub xjob.qsub >job_id_new
      JOB_ID_NEW=$(cat job_id_new | awk -F. '{print $1}')
      echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is submitted!" >>xjobs.logfile
      until [ "$LOGQUEU" == "1" ]; do
         CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.e$JOB_ID_NEW")
         if showq | grep -q $JOB_ID_NEW || [ -f $CHECK_OUTPUT_NEW ]; then         
           LOGQUEU=1
           if [ -f $CHECK_OUTPUT_NEW ]; then
             STATUS_NEW="ABORTED"
             echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
           else
             STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
             echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
           fi
         fi
         sleep 5
      done           
      LOGQUEU=0
    fi
    echo -e "------>>>>" >>xjobs.logfile
    echo -e "----------------------------------------------------------------------------------------------------------------------\n\n" >> xjobs.logfile

    #---------------------------------------------------------------------------------------------------------      
    # Get the rest of the simulation time 
    #---------------------------------------------------------------------------------------------------------      
    CURRENT_TIME=$(date  '+%b %e %T')
    diff_1 () {
            printf '%s' $(( $(date  -d"$CURRENT_TIME" +%s) -
                            $(date  -d"$START_TIME" +%s)))
    }
    echo "CT = $CURRENT_TIME , ST = $START_TIME"    
    REST_SIMULATION_TIME=`echo $SIMULATION_TIME - $(diff_1) | bc`
    echo "REST_SIMULATION_TIME: $REST_SIMULATION_TIME"
    UPDATE_TIME_2=$(echo "scale=2; $REST_SIMULATION_TIME / 4.0" | bc)
    echo "----------------------------------------------------------------------------------------------------------------------">>xjobs.logfile
    echo "------>>>>" >>xjobs.logfile
    echo "-->> $(date  '+%b %e %T'): Check status of the running job" >>xjobs.logfile
    echo "-->> $(date  '+%b %e %T'): Rest of the runnig time: $REST_SIMULATION_TIME" >>xjobs.logfile
    echo "-->> $(date  '+%b %e %T'): Recycle time =$UPDATE_TIME_2 " >>xjobs.logfile
    echo -e "------>>>>\n" >>xjobs.logfile
    
     #Check whether current running job done or aborted (Error file exist only if one of this two action are happen)  
    LOGQUEU=0
    CHECK_OUTPUT=$(echo "$SIM_NAME.e$JOB_ID")
    if [ -f $CHECK_OUTPUT ]; then 
      if  grep -q "FLEXI FINISHED" Log.$JOB_ID.sdb || ! grep -q "walltime" $CHECK_OUTPUT  ; then
        echo "********************************************************************************"  >>xjobs.logfile
        echo "$(date  '+%b %e %T'): Error occured or Flexi finished!! Please check Error-File" >>xjobs.logfile
        echo "*********************************************************************************"   >>xjobs.logfile
        exit
      else 
        #Check status of the actual subjob 
        if ! showq | grep -q $JOB_ID_NEW; then
          CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.o$JOB_ID_NEW")
          if grep -q "subjob started too early" $CHECK_OUTPUT_NEW  ; then
             echo "------>>>>" >>xjobs.logfile
             echo "-->> $(date  '+%b %e %T'): Previous subjob ($JOB_ID_NEW) started too early!" >>xjobs.logfile
             qsub xjob.qsub >job_id_new
             JOB_ID_NEW=$(cat job_id_new | awk -F. '{print $1}')
             echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is submitted!" >>xjobs.logfile
             until [ "$LOGQUEU" == "1" ]; do
                CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.e$JOB_ID_NEW")
                if showq | grep -q $JOB_ID_NEW || [ -f $CHECK_OUTPUT_NEW ]; then         
                  LOGQUEU=1
                  if [ -f $CHECK_OUTPUT_NEW ]; then
                    STATUS_NEW="ABORTED"
                    echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                  else
                    STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
                    echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                  fi
                fi
                sleep 5 
             done           
             LOGQUEU=0
             #echo "-->> $(date  '+%b %e %T'): New subjob($JOB_ID_NEW) listed in the QUEU and current state of the subjob is $STATUS_NEW" >>xjobs.logfile
             echo -e "------>>>>\n\n" >>xjobs.logfile
          else
            echo "********************************************************************"  >>xjobs.logfile
            echo " Error occured !! Please check Error-File(2)" >>xjobs.logfile
            echo "********************************************************************"   >>xjobs.logfile
            exit
          fi
        fi
      fi
    else 
      #-------------------------------------------------------------------------------------------------------      
      # Wait until the current job is finished  
      #-------------------------------------------------------------------------------------------------------      
      STATUS=$(showq | grep $JOB_ID | awk '{print $3}')
      until [ "$STATUS" != "Running" ]; do
        #-----------------------------------------------------------------------------------------------------      
        # Wait a littel bit ....
        #-----------------------------------------------------------------------------------------------------      
        STATUS=$(showq | grep $JOB_ID | awk '{print $3}')
        if [ "$STATUS" == "Running" ]; then 
          sleep $UPDATE_TIME_2
          echo "------>>>" >>xjobs.logfile
          echo "-->> wait ... until the runnig job is finished " >>xjobs.logfile  
        fi

        #-----------------------------------------------------------------------------------------------------      
        # If previous SUBJOB is started too early, another one is submitted
        #-----------------------------------------------------------------------------------------------------      
        if ! showq | grep -q $JOB_ID_NEW ; then
          CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.o$JOB_ID_NEW")
          if grep -q "subjob started too early " $CHECK_OUTPUT_NEW  ; then
            echo "-->> $(date  '+%b %e %T'): Previous subjob ($JOB_ID_NEW) started to early!" >>xjobs.logfile
            qsub xjob.qsub >job_id_new
            JOB_ID_NEW=$(cat job_id_new | awk -F. '{print $1}')
            echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is submitted!" >>xjobs.logfile
            until [ "$LOGQUEU" == "1" ]; do
              sleep 5
              CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.e$JOB_ID_NEW")
              if showq | grep -q $JOB_ID_NEW || [ -f $CHECK_OUTPUT_NEW ]; then         
                LOGQUEU=1
                if [ -f $CHECK_OUTPUT_NEW ]; then
                  STATUS_NEW="ABORTED"
                  echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                else
                  STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
                  echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                fi
              fi
            done           
            LOGQUEU=0
            echo "------>>>" >>xjobs.logfile
          else
            echo "********************************************************************"  >>xjobs.logfile
            echo " Error occured !! Please check Error-File(3)" >>xjobs.logfile
            echo "********************************************************************"   >>xjobs.logfile
            exit
          fi  
        fi
        STATUS=$(showq | grep $JOB_ID | awk '{print $3}')
        CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.e$JOB_ID_NEW")
        if [ -f $CHECK_OUTPUT_NEW ]; then
          STATUS_NEW="ABORTED"
        else
          STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
        fi
        CURRENT_TIME=$(date  '+%b %e %T')
        echo "-->> $CURRENT_TIME: status of running job ($JOB_ID): $STATUS" >>xjobs.logfile
        echo "-->> $CURRENT_TIME: subjob ($JOB_ID_NEW) is: $STATUS_NEW " >>xjobs.logfile
        echo -e "------>>>\n" >>xjobs.logfile
      done
   fi
         
   
   echo "$(date  '+%b %e %T'): current job is successfully finished ($CURRENT_TIME) !! " >>xjobs.logfile
   echo -e "----------------------------------------------------------------------------------------------------------------------\n\n " >>xjobs.logfile
   echo "---------------------------------------------------------------------------------------------------------------------" >>xjobs.logfile
        

   #-----------------------------------------------------------------------------------------------------      
   # If previous SUBJOB is started too early, another one is submitted
   #-----------------------------------------------------------------------------------------------------      
   echo "-->> $(date  '+%b %e %T'): Check status of the subjob" >>xjobs.logfile
   if [ -f $CHECK_OUTPUT ]; then
     if grep -q "FLEXI FINISHED" Log.$JOB_ID.sdb || ! grep -q "walltime" $CHECK_OUTPUT  ; then
        echo "***************************************************************************************"  >>xjobs.logfile
        echo "$(date  '+%b %e %T'): Error occured or Flexi finished!! Please check Error-Filei" >>xjobs.logfile
        echo "***************************************************************************************"   >>xjobs.logfile
        exit
     else 
       CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.o$JOB_ID_NEW")
       if ! showq |  grep -q $JOB_ID_NEW; then
         if  grep -q " $(date  '+%b %e %T'): Subjob started too early " $CHECK_OUTPUT_NEW  ; then
            echo "------>>>" >>xjobs.logfile
            echo "-->> $(date  '+%b %e %T'): Previous subjob ($JOB_ID_NEW) started to early!" >>xjobs.logfile
            qsub xjob.qsub >job_id_new
            JOB_ID_NEW=$(cat job_id_new | awk -F. '{print $1}')
            echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is submitted!" >>xjobs.logfile
            until [ "$LOGQUEU" == "1" ]; do
               CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.e$JOB_ID_NEW")
               if showq | grep -q $JOB_ID_NEW || [ -f $CHECK_OUTPUT_NEW ]; then         
                 LOGQUEU=1
                 if [ -f $CHECK_OUTPUT_NEW ]; then
                   STATUS_NEW="ABORTED"
                   echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                 else
                   STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
                   echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                 fi
               fi
               sleep 5 
            done           
            LOGQUEU=0
            echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) listed in the QUEU and current state of the subjob is $STATUS_NEW" >>xjobs.logfile
            echo -e "------>>>>\n\n" >>xjobs.logfile
         else
           echo "********************************************************************"  >>xjobs.logfile
           echo " Error occured !! Please check Error-File" >>xjobs.logfile
           echo "********************************************************************"   >>xjobs.logfile
           exit
         fi
       else
         if [ STATUS_NEW == "ABORTED" ]; then
            qsub xjob.qsub >job_id_new
            JOB_ID_NEW=$(cat job_id_new | awk -F. '{print $1}')
            echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is submitted!" >>xjobs.logfile
            until [ "$LOGQUEU" == "1" ]; do
               CHECK_OUTPUT_NEW=$(echo "$SIM_NAME.e$JOB_ID_NEW")
               if showq | grep -q $JOB_ID_NEW || [ -f $CHECK_OUTPUT_NEW ]; then         
                 LOGQUEU=1
                 if [ -f $CHECK_OUTPUT_NEW ]; then
                   STATUS_NEW="ABORTED"
                   echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                 else
                   STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
                   echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) is $STATUS_NEW!" >>xjobs.logfile
                 fi
               fi
               sleep 5 
            done           
            LOGQUEU=0
            # echo "-->> $(date  '+%b %e %T'): New subjob ($JOB_ID_NEW) listed in the QUEU and current state of the subjob is $STATUS_NEW" >>xjobs.logfile
         fi
       fi
     fi
   else 
     echo "********************************************************************"  >>xjobs.logfile
     echo " Error occured !! Please check Error-File" >>xjobs.logfile
     echo "********************************************************************"   >>xjobs.logfile
     exit
   fi
  
    
  until [ "$STATUS_NEW" == 'Running' ]; do
    echo "------>>>>" >>xjobs.logfile
    STATUS_NEW=$(showq | grep $JOB_ID_NEW | awk '{print $3}')
    echo "-->> $(date  '+%b %e %T'): Submitted Job ($JOB_ID_NEW) is $STATUS_NEW" >>xjobs.logfile
    echo "-->> Check status  ... until current subjob is running" >>xjobs.logfile
    echo -e "------>>>>\n" >>xjobs.logfile
    sleep $UPDATE_TIME
  done
 
  mv job_id_new job_id
  echo "-->>... subjob is finally submitted and is running now !!" >>xjobs.logfile
  echo -e "------------------------------------------------------------------------------------------------------------------\n" >>xjobs.logfile
  echo "*********************************************************************************************************************" >>xjobs.logfile
  echo "...  loop $i of $NUM_SUBJOBS  is finished!!" >>xjobs.logfile
  echo "*********************************************************************************************************************\n\n" >>xjobs.logfile

fi
done


