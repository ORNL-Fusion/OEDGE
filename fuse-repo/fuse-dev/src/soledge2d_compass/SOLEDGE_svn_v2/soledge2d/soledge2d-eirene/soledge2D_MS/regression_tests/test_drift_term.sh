rm -r Results*

ulimit -s unlimited

./test_drift_term 2 1
./test_drift_term 3 1
./test_drift_term 4 1
./test_drift_term 5 1
./test_drift_term 6 1
./test_drift_term 7 1
./test_drift_term 8 1
./test_drift_term 9 1

python plot_err_drift_grid1.py

rm -r Results

./test_drift_term 2 3
./test_drift_term 3 3
./test_drift_term 4 3
./test_drift_term 5 3
./test_drift_term 6 3
./test_drift_term 7 3
./test_drift_term 8 3
./test_drift_term 9 3

python plot_err_drift_grid2.py

rm -r Results

