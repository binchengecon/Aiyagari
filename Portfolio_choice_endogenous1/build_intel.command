# Builder using Intel compiler (ICC) ---
DIR="cd "$( dirname "$0" )""
DIR2=""$( dirname "$0" )""
$DIR

rm -Rf II_model

# export path for missing link.
export CPATH=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include
export LIBRARY_PATH=$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/usr/lib

export OMP_STACKSIZE=1000000000; 
ulimit -s 65532
ulimit -a
icpc main.cpp -fast -fopenmp -o II_model  -lgsl -lgslcblas;time ./II_model

