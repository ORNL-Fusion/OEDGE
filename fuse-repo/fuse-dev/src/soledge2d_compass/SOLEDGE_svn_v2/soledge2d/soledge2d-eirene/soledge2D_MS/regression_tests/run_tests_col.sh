rm *.pdf
rm -r Results Results1 Results2 Results3
mkdir errors

#cp Makefiles/Make_library_gfortran.inc Make_library.inc

#COLLISIONS GRID1
# cp  ../*/*/*.f90 .
# cp  ../*/*.f90 .
# cp  ../*.f90 .
rm -f test_collisions
make -f Makefiles/Make_collisions_test1
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

# rm  *.mod *.o *.f90
./test_collisions1.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

cp -r Results Results1

# cp  ../*/*/*.f90 .
# cp  ../*/*.f90 .
# cp  ../*.f90 .
rm -f test_collisions
make -f Makefiles/Make_collisions_test2
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

# rm  *.mod *.o *.f90
./test_collisions1.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

cp -r Results Results2

# cp  ../*/*/*.f90 .
# cp  ../*/*.f90 .
# cp  ../*.f90 .
rm -f test_collisions
make -f Makefiles/Make_collisions_test3
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

# rm  *.mod *.o *.f90
./test_collisions1.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

cp -r Results Results3

${PYTHON} python/plot_err_collisions_grid1v.py
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

${PYTHON} python/plot_err_collisions_grid1T.py
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi

rm -r Results1 Results2 Results3 Results

#COLLISIONS GRID2
# cp  ../*/*/*.f90 .
# cp  ../*/*.f90 .
# cp  ../*.f90 .
rm -f test_collisions
make -f Makefiles/Make_collisions_test1
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
# rm  *.mod *.o *.f90
./test_collisions2.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
cp -r Results Results1

# cp  ../*/*/*.f90 .
# cp  ../*/*.f90 .
# cp  ../*.f90 .
rm -f test_collisions
make -f Makefiles/Make_collisions_test2
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
# rm  *.mod *.o *.f90
./test_collisions2.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
cp -r Results Results2

# cp  ../*/*/*.f90 .
# cp  ../*/*.f90 .
# cp  ../*.f90 .
rm -f test_collisions
make -f Makefiles/Make_collisions_test3
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
# rm  *.mod *.o *.f90
./test_collisions2.sh
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi
cp -r Results Results3

${PYTHON} python/plot_err_collisions_grid2v.py
${PYTHON} python/plot_err_collisions_grid2T.py
STAT=$?;if [ $STAT -ne 0 ]; then exit 1; fi


rm -r Results1 Results2 Results3 Results
rm test_collisions
exit 0