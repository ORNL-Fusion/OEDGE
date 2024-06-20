rm -r Results

ulimit -s unlimited

./test_metric 2 1
./test_metric 3 1
./test_metric 4 1
./test_metric 5 1
./test_metric 6 1
./test_metric 7 1
./test_metric 8 1
./test_metric 9 1

$PYTHON python/plot_err_metric_grid1.py

./test_metric 2 3
./test_metric 3 3
./test_metric 4 3
./test_metric 5 3
./test_metric 6 3
./test_metric 7 3
./test_metric 8 3
./test_metric 9 3

$PYTHON python/plot_err_metric_grid3.py

#rm -r Results
