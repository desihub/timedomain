#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=240
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
    echo Loading DESI modules
    module use /global/common/software/desi/$NERSC_HOST/desiconda/startup/modulefiles
    echo module load
    module load desimodules/master
fi

#Look for new files from last night

run_path="/global/cscratch1/sd/akim/project/timedomain/timedomain/bin/"
lastnite=$( date -d "yesterday 13:00" '+%Y%m%d' )

echo "${run_path}/diff.py $lastnite CVLogic Date_SpectraPairs_Iterator daily coadd"

srun ${run_path}/diff.py $lastnite CVLogic Date_TargetPairs_Iterator daily spectra

echo "end"
#If a run has successfully finished, we want to add the tile+date in a file to keep track 
