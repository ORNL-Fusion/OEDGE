#!/bin/bash

#if [ -z "$PYTHON" ]; then
  export PYTHON=python3
#fi

rm *.o *.mod *.f90

rm *.pdf
rm -r errors
mkdir errors

#PARALLEL ADVECTION
cp  ../*/*/*.f90 .
cp  ../*/*.f90 .
cp  ../*.f90 .
make -f Makefiles/Make_paradif_test
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

#rm  *.mod *.o *.f90
./test_paradif.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
rm test_para_diffusion

rm *.aux *.log 
exit 0
