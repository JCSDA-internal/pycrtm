#!/bin/csh
set myshell=`echo $shell | rev | cut -d"/" -f1 | rev`
source $MODULESHOME/init/$myshell
module purge
setenv OPT /discover/swdev/jcsda/modules
module use $OPT/modulefiles
module use $OPT/modulefiles/apps
module use $OPT/modulefiles/core
module load jedi/intel-impi/19.1.0.166-v0.4
module list
