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
make -f Makefiles/Make_perpdif_test
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

#rm  *.mod *.o *.f90
./test_perpdif.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
rm test_perp_adv

rm *.aux *.log 
exit 0
