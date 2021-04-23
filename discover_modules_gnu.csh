#!/bin/csh
set myshell=`echo $shell | rev | cut -d"/" -f1 | rev`
source $MODULESHOME/init/$myshell
module purge
setenv OPT /discover/swdev/jcsda/modules
module use $OPT/modulefiles
module use $OPT/modulefiles/apps
module use $OPT/modulefiles/core
module load jedi/gnu-impi/9.2.0
module list
