
export CONDA_DIR=/cvmfs/des.opensciencegrid.org/fnal/anaconda2

source $CONDA_DIR/etc/profile.d/conda.sh

conda activate des18a

echo "Environment is ready. To start jupyter notebook server, run"

echo "jupyter notebook --no-browser --port=8889"


