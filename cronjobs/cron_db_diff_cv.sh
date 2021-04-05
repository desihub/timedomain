#!/bin/bash
## This is similar to cron_test.sh but is not called with sbatch, 
## Because it runs sbatch itself
## It interacts with the sqlite3 db called 
## /global/cfs/cdirs/desi/science/td/daily-search/transients_search.db
## Blame Antonella Palmese version Jan 2021 for the ugliness of this code

echo `date` Running daily time domain pipeline on `hostname`
#- Configure desi environment if needed
if [ -z "$DESI_ROOT" ]; then
    echo "Loading DESI modules"
    module use /global/common/software/desi/$NERSC_HOST/desiconda/startup/modulefiles
    echo "module load"
    module load desimodules/master
fi

tiles_path="/global/project/projectdirs/desi/spectro/redux/daily/tiles"
run_path="/global/cscratch1/sd/akim/project/timedomain/cronjobs/"
td_path="/global/cfs/cdirs/desi/science/td/daily-search/"


##################
#Now double check that we have run over all new exposures. If not, we send new runs.
#The database is updated in the next lines starting where it calls python ${run_path}exposure_db.py
#That script updates the exposures table, which is later compared to the processed exposures in 
#desitrip_exposures to find unprocessed exposures. 
#desitrip_exposures is then updated with new exposures that went through the classifier in the classifier script
#ATM this only works for desitrip outputs - needs to be added to specdiff

echo "Looking for new exposures"


python ${run_path}exposure_db.py daily

# query="select distinct obsdate,tileid from exposures
# where (tileid,obsdate) not in (select tileid,obsdate from desidiff_cv_spectra_exposures);"
query="select distinct obsdate,tileid from exposures
where (tileid,obsdate) not in (select tileid,obsdate from desidiff_cv_coadd_exposures);"

mapfile -t -d $'\n' obsdates_tileids < <( sqlite3 ${td_path}transients_search.db "$query" )


#Prepare sbatch file
now=$( date -d "today" '+%Y%m%d_%T' )
logfile="${td_path}desitrip/log/${now}.log"
echo "Putting log into "$logfile

#I think this is not clever at the moment: it sends different jobs for each obsdate,
#but most of the time is spent importing tensorflow. 
#It would be more clever to run multiple obsdate together

echo "#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=120
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=$logfile
#SBATCH --constraint=haswell
">${run_path}sbatch_file.sh


Nobsdates_tileids=${#obsdates_tileids[@]}
if [ $Nobsdates_tileids -eq 0 ]; then
    echo "No new observations found today"
else
    echo "$Nobsdates_tileids new observations found"

    echo "---------- Starting coadd differencing ----------"


    
    run_path_diff="/global/cscratch1/sd/akim/project/timedomain/timedomain/bin/"
    logfile="${td_path}/desitrip/log/${now}.log"

    # instead of doing all the date/tile pairs at once, split into pieces
    # the purpose is to have more things get processed/saved in case of
    # any kind of error
    nper=1
    nloop=$(((Nobsdates_tileids+nper-1)/nper))

    for ((i=0;i<$nloop;i++)); 
        do 
            subarr=("${obsdates_tileids[@]:$(($i*$nper)):$nper}")
            echo "${subarr[@]}"

    #     echo "${run_path_diff}/diff-db.py $lastnite CVLogic Date_SpectraPairs_Iterator daily
    #     coadd"
    #     srun -o ${logfile} ${run_path_diff}/diff.py $lastnite CVLogic
    #     Date_SpectraPairs_Iterator daily coadd
    
#             python ${run_path_diff}diff-db.py TileDate_TargetPairs_Iterator CVLogic daily spectra --obsdates_tilenumbers ${subarr[@]}
            python ${run_path_diff}diff-db.py TileDate_SpectraPairs_Iterator CVLogic daily coadd --obsdates_tilenumbers ${subarr[@]}
            if [ $? -eq 0 ]
            then
                echo "Successfully executed script"
                #Now add this tile info to the sqlite db
                for t in ${subarr[@]}; do
                    arrt=(${t//|/ })
#                     query="INSERT OR IGNORE INTO desidiff_cv_spectra_exposures(obsdate, tileid) VALUES(${arrt[0]},${arrt[1]});"
                    query="INSERT OR IGNORE INTO desidiff_cv_coadd_exposures(obsdate, tileid) VALUES(${arrt[0]},${arrt[1]});"
                    echo $query
                    sqlite3 ${td_path}transients_search.db "$query"
                done
            else
              # Redirect stdout from echo command to stderr.
                echo "Script encountered error." >&2
#               echo "Failure in $query" |  mail -s 'Failure: cron_db_diff.sh' agkim@lbl.gov
                exit 1
            fi
        done

fi



##################

#I think we should do this on a separate sbatch to avoid hitting walltime

#echo "---------- Starting coadd differencing ----------"

#run_path_diff="/global/u2/p/palmese/desi/timedomain/timedomain/bin/"
#logfile="${td_path}/desitrip/log/${now}.log"

#echo "${run_path_diff}/diff.py $lastnite CVLogic Date_SpectraPairs_Iterator daily coadd"

#srun -o ${logfile} ${run_path_diff}/diff.py $lastnite CVLogic Date_SpectraPairs_Iterator daily coadd

#echo "Jobs for coadd differencing sent"

