rm -r Results

ulimit -s unlimited

./test_para_diffusion 2 1
./test_para_diffusion 3 1
./test_para_diffusion 4 1
./test_para_diffusion 5 1
./test_para_diffusion 6 1
./test_para_diffusion 7 1
#./test_para_diffusion 8 1
#./test_para_diffusion 9 1

${PYTHON} python/plot_err_paradif_grid1.py

./test_para_diffusion 2 2
./test_para_diffusion 3 2
./test_para_diffusion 4 2
./test_para_diffusion 5 2
./test_para_diffusion 6 2
./test_para_diffusion 7 2
#./test_para_diffusion 8 2
#./test_para_diffusion 9 2

${PYTHON} python/plot_err_paradif_grid2.py

rm -r Results
