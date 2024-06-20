#!/bin/bash

if [ -z "$PYTHON" ]; then
  export PYTHON=python
fi

rm *.o *.mod *.f90

rm *.pdf
rm -r errors
mkdir errors

#cp Makefiles/Make_library_gfortran.inc Make_library.inc

#METRIC
cp  ../*/*/*.f90 .
cp  ../*/*.f90 .
cp  ../*.f90 .
make -f Makefiles/Make_metric_test
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

#rm  *.mod *.o *.f90
./test_metrics.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

rm test_metric

#PARALLEL DIFFUSION
#cp  ../*/*/*.f90 .
#cp  ../*/*.f90 .
#cp  ../*.f90 .
make -f Makefiles/Make_paradif_test
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

#rm  *.mod *.o *.f90
./test_paradif.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
rm test_para_diffusion

#PARALLEL ADVECTION
#cp  ../*/*/*.f90 .
#cp  ../*/*.f90 .
#cp  ../*.f90 .
make -f Makefiles/Make_advection_test
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

#rm  *.mod *.o *.f90
./test_advection.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
rm test_para_adv

#PERPENDICULER DIFFUSION
#cp  ../*/*/*.f90 .
#cp  ../*/*.f90 .
#cp  ../*.f90 .
make -f Makefiles/Make_perpdif_test
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
# rm  *.mod *.o *.f90
./test_perpdif.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
rm test_perp_dif

#COLLISIONS
./run_tests_col.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

pdflatex regression_report.tex
rm *.aux *.log Make_library.inc
rm -r error
exit 0
