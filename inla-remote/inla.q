#!/usr/bin/env bash

CYGWIN=0

## This is the inla.q-script for handling the remote inla-queue.  This
## version works both for Linux/Mac and Windows under CYGWIN. We just
## have to add the environment

maketemp () {
    mktemp -t inla.q.XXXXXXXX
}
makedtemp () {
    mktemp -d -t inla.q.XXXXXXXX
}

compress_put="z"
compress_get="z"
Logfile="Logfile.txt"
## this one is fixed, do not change
SSHDefaultPort=22

if [ "${INLA_DEBUG}XX" != "XX" ]; then
    set -vx
    echo "INLA_PATH $INLA_PATH"
    echo "INLA_HOME $INLA_HOME"
    echo "INLA_SSH_AUTH_SOCK $INLA_SSH_AUTH_SOCK"
    echo "INLA_OS $INLA_OS"
    echo "INLA_CYGWIN_HOME $INLA_CYGWIN_HOME"
    echo "INLA_HGVERSION $INLA_HGVERSION"
fi

if [ "$INLA_OS" == "windows" ]; then
    INIFILE=$INLA_HOME/.inlarc
    export SSH_AUTH_SOCK="$INLA_SSH_AUTH_SOCK"
    if [ "$INLA_CYGWIN_HOME"XX != "XX" ]; then
	export HOME=$INLA_CYGWIN_HOME
    fi
else
    INIFILE=~/.inlarc
    if [ "$SSH_AUTH_SOCK"XX == "XX" ]; then
	## this is not set, weird, but it might be that a manual setup is used
	export SSH_AUTH_SOCK="$INLA_SSH_AUTH_SOCK"
    fi   
fi

## these are the variables to be set by the user
RemoteINLA="/usr/local/bin/inla"
RemoteHost=inla.math.ntnu.no
RemoteUser=$USER
Port=$SSHDefaultPort
sshArguments="-x -o BatchMode=yes -o TCPKeepAlive=yes -e none"

if [ -r "$INIFILE" ]; then
    if [ "$INLA_OS" == "windows" ]; then
	TMP=$(maketemp)
	tr -d \\015 < "$INIFILE" > $TMP
	source $TMP
	rm -rf $TMP
    else
	source "$INIFILE"
    fi
else
    echo -e "\n\n\n$0: No such file $INIFILE"
    echo -e "Run command in R: inla.remote()"
    echo -e "\n\n"
    exit 1
fi

cmd="$1"
id="$2"
remove="$3"
rdir=tmp/.inla.remote
tarfile=results$RANDOM$RANDOM.tar
tarfile_to=$(maketemp).tar
logfile_to=$(maketemp).log
dirto=$(makedtemp)
if [ "$id" = "NULL" ]; then
    no=0
else
    no=$(echo "$id" | awk '{print int($1)}')
    if [ $no -eq 0 ]; then
	no=-1
    fi
fi  

if [ "$cmd" = "stat" ]; then

   ssh -p$Port $sshArguments $RemoteUser@$RemoteHost "\
        mkdir -p $rdir; \
        cd $rdir; \
	nno=0; \
    	for d in \$(ls -1 .); do \
            if [ -d \$d -a -f \$d/jobid -a -f \$d/.inla.pid -a \! -f \$d/working ]; then \
	        nno=\$[ \$nno + 1 ]; \
		siz=\$(du -sm \$d | cut	-f1); \
	        if [ $no -eq 0 -o $no -eq \$nno -o \$(cat \$d/jobid) = "$id" ]; then \
       	            if [ -f \$d/done ]; then \
	                status="\""Finished"\""; \
                    elif [ \$(echo \$(ps -p \$(cat \$d/.inla.pid) -o comm=) | grep -s inla) ]; then \
	                status="\""Running(\$(ps -p \$(cat \$d/.inla.pid) -o time=))"\""; \
	            else	\
	                status="\""Aborted"\""; \
                    fi; \
	            echo \"\$(cat \$d/jobid) \$nno  \$(cat \$d/.inla.pid) \$status \${siz}Mb \"; \
                fi; \
            fi; \
        done"

elif [ "$cmd" = "get" ]; then    

    ssh -p$Port $sshArguments $RemoteUser@$RemoteHost "
            cd $rdir; \
	    nno=0; \
    	    for d in \$(ls -1 .); do \
                if [ -d \$d -a -f \$d/jobid -a -f \$d/.inla.pid -a \! -f \$d/working ]; then \
	            nno=\$[ \$nno + 1 ]; \
	            if [ $no -eq \$nno -o \$(cat \$d/jobid) = "$id" ]; then \
	                 cd \$d; \
			 if [ -f done ]; then \
			     if \[ -f $Logfile \]; then cp $Logfile results.files; fi; \
			     tar cf ../$tarfile results.files 2>/dev/null ; \
			     if [ $remove -eq 1 ]; then \
                                 cd ..; rm -rf \$d; \
                             fi; \
			 fi; \
	            fi; \
                fi; \
             done"
    scp -P$Port -B -C -p -q $RemoteUser@$RemoteHost:$rdir/$tarfile "$tarfile_to" >/dev/null 2>&1 || \
        {  echo "ERROR Job is not yet finished or does not exist; try 'inla.qstat()'"; exit; }
    ssh -p$Port $sshArguments $RemoteUser@$RemoteHost "rm -f $rdir/$tarfile" 
    tar xfm  "$tarfile_to" -C $dirto 2>/dev/null 
    rm -f "$tarfile_to"
    rm -f "${tarfile_to%.tar}"
    echo "$dirto"

elif [ "$cmd" = "log" ]; then    

    ssh -p$Port $sshArguments $RemoteUser@$RemoteHost "
            cd $rdir; \
	    nno=0; \
    	    for d in \$(ls -1 .); do \
                if [ -d \$d ]; then \
	            nno=\$[ \$nno + 1 ]; \
	            if [ $no -eq \$nno -o \$(cat \$d/jobid) = "$id" ]; then \
		       if [ -f \$d/$Logfile ]; then \
                           cp \$d/$Logfile .; \
                       else \
		           touch $Logfile; \
                       fi; \
	            fi; \
                fi; \
             done"
    scp -P$Port -B -C -p -q $RemoteUser@$RemoteHost:$rdir/$Logfile "$logfile_to" >/dev/null 2>&1 
    ssh -p$Port $sshArguments $RemoteUser@$RemoteHost "rm -f $rdir/$Logfile" 
    echo "LOG $logfile_to" 

elif [ "$cmd" = "del" ]; then
    
    ssh -p$Port $sshArguments $RemoteUser@$RemoteHost "
            cd $rdir; \
	    nno=0; \
    	    for d in \$(ls -1 .); do \
                if [ -d \$d -a -f \$d/jobid ]; then \
	            nno=\$[ \$nno + 1 ]; \
	            if [ $no -eq \$nno -o \$(cat \$d/jobid) = "$id" ]; then \
	                if [ -f \$d/working ]; then \
                            kill \$(cat \$d/.inla.pid); \
                        fi; \
	                rm -rf \$d; \
                    fi; \
                fi; \
             done" > /dev/null 2>&1

    if [ $no -gt 0 ]; then
	echo "DELETE $no"
    else
	echo "DELETE $id"
    fi

elif [ "$cmd" = "nuke" ]; then
    
    del="$no"
    ssh -p$Port $sshArguments $RemoteUser@$RemoteHost "
            cd $rdir; \
    	    for d in \$(ls -1 .); do \
                if [ -d \$d -a -f \$d/jobid ]; then \
	            if [ -f \$d/working ]; then \
                        kill \$(cat \$d/.inla.pid); \
                    fi; \
	            rm -rf \$d; \
                fi; \
             done" > /dev/null 2>&1
    echo "NUKE"

fi
