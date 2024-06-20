rm -r Results

ulimit -s unlimited

./test_perp_dif 2 1
./test_perp_dif 3 1
./test_perp_dif 4 1
./test_perp_dif 5 1
./test_perp_dif 6 1
./test_perp_dif 7 1
#./test_perp_dif 8 1
#./test_perp_dif 9 1

${PYTHON} python/plot_err_perpdif_grid1.py

./test_perp_dif 2 2
./test_perp_dif 3 2
./test_perp_dif 4 2
./test_perp_dif 5 2
./test_perp_dif 6 2
./test_perp_dif 7 2
#./test_perp_dif 8 2
#./test_perp_dif 9 2

${PYTHON} python/plot_err_perpdif_grid2.py

rm -r Results
