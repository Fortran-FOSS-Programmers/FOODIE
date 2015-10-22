#!/bin/bash
# script automotion for FOODIE's OpenMP benchmark using provided tests

# defaults
code='all'
bench='all'
ripetitions=1
Ni=100

# functions
function print_usage {
  # print help message for correct CLI usage
  echo
  echo "`basename $0`"
  echo "script automotion for FOODIE OpenMP benchmark"
  echo "Usage: `basename $0` [opts [args]]"
  echo "  [ -c #code_type -b #benchmark_type -r #ripetitions -Ni #cells ]"
  echo
  echo "Defaults of optional arguments:"
  echo "  -c all # type of code: valid value are 'all', 'foodie', 'no-foodie'"
  echo "  -b all # type of benchmark: valid value are 'all', 'strong', 'weak'"
  echo "  -r 1 # number of ripetitions for averaging the benchmarks results"
  echo "  -Ni 1000 # test size for strong scaling or minimum test size for weak scaling benchamrk"
  echo
  echo "Examples:"
  echo "  `basename $0` # perform $ripetitions ripetitions with default test size of $Ni cells"
  echo "  `basename $0` -r 5 # perform 5 ripetitions and compute the average"
  echo "  `basename $0` -Ni 2000 -r 10 # use a test size of 2000 cells"
  echo
}
function run_exe {
  # perform the code execution
  exe=$1
  thd=$2
  size=$3
  average_cpu_time=0.
  for r in $( seq $ripetitions ); do
    cpu_time=$($exe --Ni $size --omp_threads $thd | awk '{print $2}')
    average_cpu_time=`echo "scale=10; $average_cpu_time+$cpu_time" | bc -l`
  done
  average_cpu_time=`echo "scale=10; $average_cpu_time/$ripetitions" | bc -l`
  echo $average_cpu_time
}
function strong_scaling {
  # perform strong scaling analysis
  foodie=$1
  if [ "$foodie" == "yes" ]; then
    exe='../tests-openmp/euler-1D-openmp'
    output='strong-scaling.dat'
  else
    exe='../tests-openmp/euler-1D-openmp-no-foodie'
    output='strong-scaling-no-foodie.dat'
  fi
  echo 'TITLE="Strong scaling analysis, cells '$Ni'"' > $output
  echo 'VARIABLES="OpenMP threads number" "CPU time"' >> $output
  echo "Benchmarking without OpenMP threads"
  average_cpu_time=$(run_exe $exe'-off' 1 $Ni)
  echo 0 $average_cpu_time >> $output
  c=1
  while [ "$c" -lt $ncores ]; do
    echo "Benchmarking with "$c" OpenMP threads"
    average_cpu_time=$(run_exe $exe'-on' $c $Ni)
    echo $c $average_cpu_time >> $output
    c=$((c*2))
    if [ "$c" -gt $ncores ]; then
      # doing the last benchmark
      echo "Benchmarking with "$ncores" OpenMP threads"
      average_cpu_time=$(run_exe $exe'-on' $ncores $Ni)
      echo $ncores $average_cpu_time >> $output
    fi
  done
}
function weak_scaling {
  # perform weak scaling analysis
  foodie=$1
  if [ "$foodie" == "yes" ]; then
    exe='../tests-openmp/euler-1D-openmp'
    output='weak-scaling.dat'
  else
    exe='../tests-openmp/euler-1D-openmp-no-foodie'
    output='weak-scaling-no-foodie.dat'
  fi
  echo 'TITLE="Weak scaling analysis"' > $output
  echo 'VARIABLES="Size" "OpenMP threads number" "CPU time"' >> $output
  echo "Benchmarking without OpenMP threads and size "$Ni
  average_cpu_time=$(run_exe $exe'-off' 1 $Ni)
  echo $Ni 0 $average_cpu_time >> $output
  c=1
  s=$Ni
  while [ "$c" -lt $ncores ]; do
    echo "Benchmarking with "$c" OpenMP threads and size "$s
    average_cpu_time=$(run_exe $exe'-on' $c $s)
    echo $s $c $average_cpu_time >> $output
    c=$((c*2))
    s=$((c*Ni))
    if [ "$c" -gt $ncores ]; then
      # doing the last benchmark
      s=$((ncores*Ni))
      echo "Benchmarking with "$ncores" OpenMP threads and size "$s
      average_cpu_time=$(run_exe $exe'-on' $ncores $s)
      echo $s $ncores $average_cpu_time >> $output
    fi
  done
}

# parsing CLI
while [ $# -gt 0 ]; do
  case "$1" in
  "-c")
    shift; code=$1
    if [ "$bench" != "all" ] && [ "$bench" != "foodie" ] && [ "$bench" != "no-foodie" ] ; then
      echo; echo "invalid value $1"; print_usage; exit 1
    fi
    ;;
  "-b")
    shift; bench=$1
    if [ "$bench" != "all" ] && [ "$bench" != "strong" ] && [ "$bench" != "weak" ] ; then
      echo; echo "invalid value $1"; print_usage; exit 1
    fi
    ;;
  "-r")
    shift; ripetitions=$1
    ;;
  "-Ni")
    shift; Ni=$1
    ;;
  *)
    echo; echo "Unknown switch $1"; print_usage; exit 1
    ;;
  esac
  shift
done

# building codes
echo "Building 1D Euler OpenMP test executable"
date > builds.log
cd ../../../
if [ "$code" == "all" ] || [ "$code" == "foodie" ]; then
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/obj
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/mod
  FoBiS.py build -f src/tests/parallel/euler-1D-openmp/fobos -mode benchmark-openmp-on --build_dir papers/announcement/openmp_benchmark/tests-openmp >> papers/announcement/openmp_benchmark/builds.log
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/obj
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/mod
  FoBiS.py build -f src/tests/parallel/euler-1D-openmp/fobos -mode benchmark-openmp-off --build_dir papers/announcement/openmp_benchmark/tests-openmp >> papers/announcement/openmp_benchmark/builds.log
fi
if [ "$code" == "all" ] || [ "$code" == "no-foodie" ]; then
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/obj
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/mod
  FoBiS.py build -f src/tests/parallel/euler-1D-openmp-no-foodie/fobos -mode benchmark-openmp-on --build_dir papers/announcement/openmp_benchmark/tests-openmp >> papers/announcement/openmp_benchmark/builds.log
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/obj
  rm -rf papers/announcement/openmp_benchmark/tests-openmp/mod
  FoBiS.py build -f src/tests/parallel/euler-1D-openmp-no-foodie/fobos -mode benchmark-openmp-off --build_dir papers/announcement/openmp_benchmark/tests-openmp >> papers/announcement/openmp_benchmark/builds.log
fi
cd -
# benchmarks
tecplot=$(which tec360)
ncpus=`grep "^physical id" /proc/cpuinfo | sort -u | wc -l`
ncores_cpu=`grep "^core id" /proc/cpuinfo | sort -u | wc -l`
ncores=$((ncpus*ncores_cpu))
echo "Available physical cores "$ncores

echo "1D Euler test benchmark"
mkdir -p results-euler-1D-openmp
cd results-euler-1D-openmp
if [ "$code" == "all" ] || [ "$code" == "foodie" ]; then
  if [ "$bench" == "all" ] || [ "$bench" == "strong" ]; then
    strong_scaling yes
  fi
  if [ "$bench" == "all" ] || [ "$bench" == "weak" ]; then
    weak_scaling yes
  fi
fi
if [ "$code" == "all" ] || [ "$code" == "no-foodie" ]; then
  if [ "$bench" == "all" ] || [ "$bench" == "strong" ]; then
    strong_scaling no
  fi
  if [ "$bench" == "all" ] || [ "$bench" == "weak" ]; then
    weak_scaling no
  fi
fi
if [ -x "$tecplot" ] ; then
  ln -fs ../utilities-euler-1D-openmp/*lay .
  ln -fs ../../utilities/lay_export* .
  for file in $( ls *lay ); do
  ./lay_export_all_f $file
  done
fi
cd -

exit 0
