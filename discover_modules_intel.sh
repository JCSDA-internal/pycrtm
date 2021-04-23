#!/bin/bash
myshell=`echo $SHELL | rev | cut -d"/" -f1 | rev`
source $MODULESHOME/init/$myshell
module purge
export OPT='/discover/swdev/jcsda/modules'
module use $OPT/modulefiles
module use $OPT/modulefiles/apps
module use $OPT/modulefiles/core
module load jedi/intel-impi/19.1.0.166-v0.4
module list
