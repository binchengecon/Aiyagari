# Builder using GCC ---
DIR="cd "$( dirname "$0" )""
DIR2=""$( dirname "$0" )""
$DIR

rm -rf II_model 
export OMP_STACKSIZE=1000000000; 
ulimit -s 65532
ulimit -a
xcrun g++-10 -g main.cpp -Ofast -fopenmp -o II_model -I/usr/local/include  -L/usr/local/lib -lgsl -lgslcblas;time ./II_model


