#!/bin/bash
# script automotion for FOODiE's OpenMP benchmark using provided tests
cd ../../../
echo "Building 1D Euler OpenMP test executable"
FoBiS.py build -mode openmp-benchmark --build_dir papers/announcement/openmp_benchmark/tests-openmp > buils.log
cd -
echo "1D Euler test benchmark"
rm -rf results-1D-euler-opnemp
mkdir results-1D-euler-opnemp
./tests-openmp/euler-1D-openmp --Ni 1000 --omp_threads 1  > results-1D-euler-opnemp/summary-01thds.log
./tests-openmp/euler-1D-openmp --Ni 1000 --omp_threads 2  > results-1D-euler-opnemp/summary-02thds.log
./tests-openmp/euler-1D-openmp --Ni 1000 --omp_threads 4  > results-1D-euler-opnemp/summary-04thds.log
./tests-openmp/euler-1D-openmp --Ni 1000 --omp_threads 8  > results-1D-euler-opnemp/summary-08thds.log
./tests-openmp/euler-1D-openmp --Ni 1000 --omp_threads 16 > results-1D-euler-opnemp/summary-16thds.log
exit 0
