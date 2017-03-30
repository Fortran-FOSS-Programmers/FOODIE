#!/bin/bash
# script automotion for FOODIE's errors analysis

# solvers enabling flags (1 => enabled, disabled otherwise)
ab=0
abm=0
am=0
bdf=0
emd_rk=0
euler=0
lf=0
lf_raw=0
lmm_ssp=0
ls_rk=0
tvd_rk=0

function print_usage {
  echo
  echo "`basename $0`"
  echo "script automotion for FOODIE's errors analysis"
  echo "Usage: `basename $0` #solver"
  echo "Valid solver switches are:"
  echo "  -ab      # => adams-bashforth"
  echo "  -abm     # => adams-bashforth-moulton"
  echo "  -am      # => adams-moulton"
  echo "  -bdf     # => backward differentiation formula"
  echo "  -emd-rk  # => embeded RK"
  echo "  -euler   # => euler"
  echo "  -lf      # => leapfrog"
  echo "  -lf-raw  # => leapfrog raw filtered"
  echo "  -lmm-ssp # => linear multistep methods SSP"
  echo "  -ls-rk   # => low storage RK"
  echo "  -tvd-rk  # => TVD RK"
  echo "  -all     # => all solvers!"
}

#parsing command line
while [ $# -gt 0 ]; do
  case "$1" in
    "-ab")
      ab=1
      ;;
    "-abm")
      abm=1
      ;;
    "-am")
      am=1
      ;;
    "-bdf")
      bdf=1
      ;;
    "-emd-rk")
      emd_rk=1
      ;;
    "-euler")
      euler=1
      ;;
    "-lf")
      lf=1
      ;;
    "-lf-raw")
      lf_raw=1
      ;;
    "-lmm-ssp")
      lmm_ssp=1
      ;;
    "-ls-rk")
      ls_rk=1
      ;;
    "-tvd-rk")
      tvd_rk=1
      ;;
    "-all")
      ab=1
      abm=1
      am=1
      bdf=1
      emd_rk=1
      euler=1
      lf=1
      lf_raw=1
      ls_rk=1
      tvd_rk=1
      ;;
    *)
      echo; echo "Unknown switch $1"; print_usage; exit 1
      ;;
  esac
  shift
done

echo "Building tests executables"
date > builds.log
cd ../../../
FoBiS.py build -f src/tests/accuracy/oscillation/fobos -mode errors-analysis --build_dir papers/announcement/errors_analysis/tests >> papers/announcement/errors_analysis/builds.log
cd -

echo "Oscillation test errors analysis"
# rm -rf results-oscillation
mkdir -p results-oscillation
cd results-oscillation
if [ $ab -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver adams-bashforth --ss 1 4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis-adams-bashforth.log
  ../tests/oscillation -r -tf 1e6 --solver adams-bashforth --ss 1 6 -Dt 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis-adams-bashforth.log
  ../tests/oscillation -r -tf 1e7 --solver adams-bashforth --ss 7 -Dt 400.d0 250.d0 200.d0 --errors_analysis >> analysis-adams-bashforth.log
  ../tests/oscillation -r -tf 1e6 --solver adams-bashforth --ss 8 -Dt 250.d0 200.d0 100.d0 --errors_analysis >> analysis-adams-bashforth.log
  ../tests/oscillation -r -tf 1e5 --solver adams-bashforth --ss 9 -Dt 160.d0 80.d0 --errors_analysis >> analysis-adams-bashforth.log
  ../tests/oscillation -r -tf 1e5 --solver adams-bashforth --ss 10 -Dt 50.d0 1.d0 --errors_analysis >> analysis-adams-bashforth.log
  ../tests/oscillation -r -tf 1e3 --solver adams-bashforth --ss 12 -Dt 1.d0 0.01d0 --errors_analysis >> analysis-adams-bashforth.log
  ../tests/oscillation -r -tf 1e2 --solver adams-bashforth --ss 13 -Dt 0.2d0 0.001d0 --errors_analysis >> analysis-adams-bashforth.log
fi
if [ $abm -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver adams-bashforth-moulton --ss 1 4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e6 --solver adams-bashforth-moulton --ss 1 6 -Dt 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e7 --solver adams-bashforth-moulton --ss 7 -Dt 400.d0 250.d0 200.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e6 --solver adams-bashforth-moulton --ss 8 -Dt 250.d0 200.d0 100.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e5 --solver adams-bashforth-moulton --ss 9 -Dt 160.d0 80.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e5 --solver adams-bashforth-moulton --ss 10 -Dt 50.d0 1.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e3 --solver adams-bashforth-moulton --ss 12 -Dt 1.d0 0.01d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e2 --solver adams-bashforth-moulton --ss 13 -Dt 0.2d0 0.001d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
fi
if [ $am -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver adams-moulton --ss 1 4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e6 --solver adams-moulton --ss 1 6 -Dt 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e7 --solver adams-moulton --ss 7 -Dt 400.d0 250.d0 200.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e6 --solver adams-moulton --ss 8 -Dt 250.d0 200.d0 100.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e5 --solver adams-moulton --ss 9 -Dt 160.d0 80.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e5 --solver adams-moulton --ss 10 -Dt 50.d0 1.d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e3 --solver adams-moulton --ss 12 -Dt 1.d0 0.01d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
  ../tests/oscillation -r -tf 1e2 --solver adams-moulton --ss 13 -Dt 0.2d0 0.001d0 --errors_analysis >> analysis-adams-bashforth-moulton.log
fi
if [ $bdf -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver backward-diff-formula --ss 1 6 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis-backward-diff-formula.log
fi
if [ $emd_rk -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver emd-runge-kutta --ss 6 7 --tolerance 1e-8 1e-9 --errors_analysis > analysis-emd-runge-kutta.log
  ../tests/oscillation -r -tf 1e6 --solver emd-runge-kutta --ss 9 --tolerance 1e-8 1e-9 --errors_analysis >> analysis-emd-runge-kutta.log
  ../tests/oscillation -r -tf 1e6 --solver emd-runge-kutta --ss 17 --tolerance 1e-8 1e-9 --errors_analysis >> analysis-emd-runge-kutta.log
fi
if [ $euler -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver euler --ss 1 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis-euler.log
fi
if [ $lf -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver leapfrog --ss 2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis-leapfrog.log
fi
if [ $lf_raw -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver leapfrog-raw --ss 2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis #> analysis-leapfrog.log
fi
if [ $lmm_ssp -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver lmm-ssp --ss 3 5 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis > analysis-lmm-ssp.log
fi
if [ $ls_rk -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver ls-runge-kutta --ss 1 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis #> analysis-leapfrog.log
  ../tests/oscillation -r -tf 1e6 --solver ls-runge-kutta --ss 5 7 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis #> analysis-leapfrog.log
  ../tests/oscillation -r -tf 1e6 --solver ls-runge-kutta --ss 12 14 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis #> analysis-leapfrog.log
fi
if [ $tvd_rk -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --solver tvd-runge-kutta --ss 1 3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis #> analysis-leapfrog.log
  ../tests/oscillation -r -tf 1e6 --solver tvd-runge-kutta --ss 5 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis #> analysis-leapfrog.log
fi
cd -
exit 0
