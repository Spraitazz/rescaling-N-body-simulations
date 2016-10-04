gcc -Wall driver.c -o driver_data.out -lm -lfftw3 -lgsl -lgslcblas 
GSL_RNG_TYPE="taus" GSL_RNG_SEED=123  ./driver_data.out
#./driver_data.out
#python3 plot.py
