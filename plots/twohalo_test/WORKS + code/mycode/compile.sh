#gcc -std=c11 main.c -o main.out -lm -L/shome/jonasp/usr/lib/ -lfftw3 -lgsl -lgslcblas
gcc -o main.out main.c -lm -lgsl -lgslcblas -lfftw3 -std=c11
#-I/home/jonas/libs/gsl/include -I/home/jonas/libs/gsl/include/gsl  -L/home/jonas/libs/gsl/lib

