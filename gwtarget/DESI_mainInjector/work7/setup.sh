unset PYTHONPATH
export DESGW_DIR=/global/homes/p/portmanm/timedomain/gwtarget/DESI_mainInjector/Main-Injector-master/python/
#/data/des60.b/data/palmese/DESI/DESI_mainInjector/Main-Injector-master/python
export DESGW_DATA_DIR=${DESGW_DIR}../data/
export PYTHONPATH=$DESGW_DIR:$PYTHONPATH
export PATH=$DESGW_DIR:${DESGW_DIR}equalArea:$PATH
source /opt/intel/bin/compilervars.sh intel64

