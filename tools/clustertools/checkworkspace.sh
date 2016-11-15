#!/bin/bash
##############################################################################################################################
##
## Author: Thomas Bolemann, IAG, Stuttgart
##
## This script will login onto the specified host via ssh and get the list of workspaces with their expiry date.
## Sends weekly mails with calendar messeges with workspaces expiry dates and warns if and are soon to expire and 
## will try to automatically extend them if possible.a Built for Stuttgart HLRS Cray clusters
##
## Usage: 
## Place the script on a (reliable 24x7) server where you have access and create a cronjob (e.g. check daily at 6:00):
##
##  00 6 * * * ~/checkworkspace.sh iagxyzuser                    hazelhen.hww.de            bolemann@iag.uni-stuttgart.de
##
##                                 $1: your cluster login name,  $2: the clusters address,  $3: your email address
##
## You can set the number of remaining days, when to start email warnings or extend the workspace automatically 
## and the weekday when to send the reports.
##
##############################################################################################################################

fail () {
  SUBJECT="WARNING: WORKSPACE CHECK FOR ${USERNAME} ON ${HOST} HAS FAILED !!!";
  echo $TODAY > $EMAILMESSAGE;
  echo "Please check workspaces manually!" >> $EMAILMESSAGE;
  mail -s "$SUBJECT" "$EMAILADRESS" < $EMAILMESSAGE;
  exit 1;
}

reportday="Mon" # weekday (Mon,Tue,Wed,Thu,Fri,Sat,Sun) for workspace status notifications
ndayswarn=5     # number of remaining days, when to start emailing warnings
ndaysextend=3   # number of remaining days, when to auto extend workspace

if [ $# -ne 3 ] ; then
  echo "Usage: wsmail <username> <host> <email>"
  exit 1
fi

TODAY=`date`
USERNAME=$1
HOST=$2
EMAILADRESS=$3

WSDIR=`mktemp -d`
WSINFO=${WSDIR}/ws_list.info
EMAILMESSAGE=${WSDIR}/"emailmessage.txt"

ssh $USERNAME@$HOST ws_list > ${WSINFO}

# Check if getting wsinfo was successfull
if [ $? -ne 0 ] ; then
  fail
fi

wsdays=(`cat ${WSINFO} | grep days | sed -e 's/ days/#/' |cut -d '#' -f1 |rev |cut -d ' ' -f1|rev |sed ':a;N;$!ba;s/\n/ /g'`)
wsnames=(`cat ${WSINFO} | grep days | cut -d ' ' -f1`)
wsextensions=(`cat ${WSINFO} | grep extensions | cut -d ':' -f2-`)

if [ "${#wsnames[@]}" -le 0 ] || [ "${#wsnames[@]}" -ne "${#wsdays[@]}" ] || [ "${#wsnames[@]}" -ne "${#wsextensions[@]}" ]; then
  fail
fi
for ((i=0;i<${#wsdays[@]};i++)); do
  if ! [[ ${wsdays[$i]} =~ ^-?[0-9]+ ]] || ! [[ ${wsextensions[$i]} =~ ^-?[0-9]+ ]] ; then
    fail
  fi
done

#find minimum number of remaining days
expire=180
for ((i=0;i<${#wsdays[@]};i++)); do
  if [ "${wsdays[$i]}" -lt "$expire" ]; then
    expire=${wsdays[$i]}
  fi
done

#if remaining days <= ndaysextend send email!
for ((i=0;i<"${#wsdays[@]}";i++)); do
  if [ "${wsdays[$i]}" -le "$ndaysextend" ]; then
    if [ $(echo "${wsextensions[$i]} > 0"| bc) -eq 1 ]; then
      ssh $USERNAME@$HOST ws_extend "${wsnames[$i]}" "31"
      SUBJECT=" CRAY WORKSPACE ${wsnames[$i]} FOR $USERNAME HAS BEEN EXTENDED"
    else
      SUBJECT=" NO MORE EXTENSIONS LEFT FOR CRAY WORKSPACE ${wsnames[$i]} FOR $USERNAME !!!"
    fi
    echo $TODAY > $EMAILMESSAGE
    echo "" >> $EMAILMESSAGE
    mail -s "$SUBJECT" "$EMAILADRESS" < $EMAILMESSAGE
  fi
done

# Send mail once a week or if workspaces with short expiry date
weekday=`LANG=en_US date +%a`
if [ "$expire" -gt "$ndayswarn" ] && [ "$weekday" != "$reportday" ] ; then
  exit 0
fi

# Create ical notifications (mostly copied from HLRS script)
for ((i=0;i<${#wsnames[@]};i++)); do
  L=`cat ${WSINFO} | grep ^${wsnames[$i]}`
  DURATION=`echo $L | awk '{printf ( "%d %s %d %s", $(NF - 3), $(NF - 2),  $(NF - 1), $NF) }'`
  END_I_CAL_STR=`date --date="${DURATION}"  +%Y%m%dT%H%M00Z`
  END_STR=`date --date="${DURATION}"  +%c`
  CREATED=`date +%Y%m%dT%H%M00Z`
  I_CAL_FILE[$i]="$WSDIR/WS_${wsnames[$i]}.ics"

  cat > ${I_CAL_FILE[$i]} <<EOF
BEGIN:VCALENDAR
PRODID:-//HLRS Cluster Team//Workspace V2.1//EN
VERSION:2.0
METHOD:REQUEST
X-MS-OLK-FORCEINSPECTOROPEN:TRUE
BEGIN:VEVENT
CLASS:PUBLIC
CREATED:${CREATED}
DESCRIPTION:Workspace ${wsnames[$i]} will be deleted\non host ${HOST}
DTEND:${END_I_CAL_STR}
DTSTAMP:${CREATED}
DTSTART:${END_I_CAL_STR}
LAST-MODIFIED:${CREATED}
LOCATION:${HOST}
SEQUENCE:0
SUMMARY:Workspace ${wsnames[$i]}
TRANSP:OPAQUE
UID:587a1aa6-$$
X-MICROSOFT-CDO-BUSYSTATUS:BUSY
X-MICROSOFT-DISALLOW-COUNTER:TRUE
END:VEVENT
END:VCALENDAR
EOF
done

SUBJECT=" CRAY WORKSPACES FOR $USERNAME EXPIRE IN $expire DAYS !!!"
echo "Dear User," > $EMAILMESSAGE
echo "" >> $EMAILMESSAGE
echo "this is a list of your workspaces on host ${HOST}" >> $EMAILMESSAGE
cat ${WSINFO} >> $EMAILMESSAGE
echo "Date is: $TODAY" >> $EMAILMESSAGE
echo "" >> $EMAILMESSAGE
echo "Please find calendar files attached.." >> $EMAILMESSAGE
ATTACH=""
for ((i=0;i<${#wsnames[@]};i++)); do
  ATTACH+="-A ${I_CAL_FILE[$i]} "
done
mail -s "$SUBJECT" $ATTACH $EMAILADRESS < $EMAILMESSAGE

rm -r $WSDIR/*.ical

