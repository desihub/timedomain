#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=20
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --constraint=haswell

#echo $HOSTNAME
#echo ${NERSC_HOST}
#source /global/common/software/desi/desi_environment.sh 20.8
#source ${DESICONDA}/etc/profile.d/conda.sh
#conda activate
#echo $0

### set -e
echo `date` Running daily time domain pipeline on `hostname`
#- Configure desi environment if needed
if [ -z "$DESI_ROOT" ]; then
    echo "Loading DESI modules"
    module use /global/common/software/desi/$NERSC_HOST/desiconda/startup/modulefiles
    echo "module load"
    module load desimodules/master
fi

#Look for new files from last night

# tiles_path="/global/project/projectdirs/desi/spectro/redux/daily/tiles"
# run_path="/global/u2/p/palmese/desi/timedomain/cronjobs/"
run_path="/global/cscratch1/sd/akim/project/timedomain/timedomain/bin/"
tiles=$( ls $tiles_path )
lastnite=$( date -d "yesterday 13:00" '+%Y%m%d' )
mapfile -d $'\0' newobs < <(find $tiles_path -type d -name $lastnite -print0)

#If there are new files, and if they have not been analyzed before, send a python script
Nfolders=${#newobs[@]}
if [ $Nfolders -eq 0 ]; then
    echo "No folders found for last night observations"
else
    echo "$Nfolders Folders found for $lastnite"
    tilenums=()
    for thisobs in ${newobs[@]}; do
        tilenums+=(`echo $thisobs | rev | cut -d\/ -f2 | rev`)
    done
fi

echo "Tile numbers: ${tilenums[@]}"

#Here we will want to make sure that we send only BGS tiles to desitrip - reading through petals takes time

# echo "python ${run_path}/cnn_classify_data.py --tilenum ${tilenums[@]} --obsdate $lastnite"
#python ${run_path}/cnn_classify_data.py --tilenum ${tilenums[@]} --obsdate $lastnite
# srun python ${run_path}/cnn_classify_data.py --tilenum ${tilenums[@]} --obsdate $lastnite >& ${run_path}/test.log
srun ${run_path}/diff.py $lastnite CVLogic Date_SpectraPairs_Iterator daily coadd
srun ${run_path}/diff.py $lastnite CVLogic Date_SpectraPairs_Iterator daily coadd

#If a run has successfully finished, we want to add the tile+date in a file to keep track 
