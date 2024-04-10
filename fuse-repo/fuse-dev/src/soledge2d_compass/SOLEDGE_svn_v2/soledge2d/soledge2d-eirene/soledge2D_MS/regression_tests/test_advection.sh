rm -r Results

ulimit -s unlimited

./test_para_adv 2 1
./test_para_adv 3 1
./test_para_adv 4 1
./test_para_adv 5 1
./test_para_adv 6 1
./test_para_adv 7 1
#./test_para_adv 8 1
#./test_para_adv 9 1

${PYTHON} python/plot_err_advection_grid1.py

./test_para_adv 2 2
./test_para_adv 3 2
./test_para_adv 4 2
./test_para_adv 5 2
./test_para_adv 6 2
./test_para_adv 7 2
#./test_para_adv 8 2
#./test_para_adv 9 2

${PYTHON} python/plot_err_advection_grid2.py

rm -r Results
