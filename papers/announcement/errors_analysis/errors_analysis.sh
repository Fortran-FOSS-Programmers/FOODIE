#!/bin/bash
# script automotion for FOODiE's errors analysis of the provided tests
cd ../../../
echo "Building tests executables"
FoBiS.py build -mode errors-analysis --build_dir papers/announcement/errors_analysis/tests > buils.log
cd -
echo "Oscillation test errors analysis"
rm -rf results-oscillation
mkdir results-oscillation
cd results-oscillation
../tests/oscillation --errors_analysis --output errors_analysis-oscillation > summary.log
tecplot=$(which tec360)
if [ -x "$tecplot" ] ; then
  ln -fs ../utilities-oscillation/*lay .
  ln -fs ../utilities/lay_export* .
  ./lay_export_all_f errors_analysis-oscillation-adams-bashforth-2.lay
  ./lay_export_all_f errors_analysis-oscillation-adams-bashforth-3.lay
  ./lay_export_all_f errors_analysis-oscillation-adams-bashforth-4.lay
  ./lay_export_all_f errors_analysis-oscillation-leapfrog.lay
  ./lay_export_all_f errors_analysis-oscillation-ls-runge-kutta-5.lay
  ./lay_export_all_f errors_analysis-oscillation-tvd-runge-kutta-2.lay
  ./lay_export_all_f errors_analysis-oscillation-tvd-runge-kutta-3.lay
  ./lay_export_all_f errors_analysis-oscillation-tvd-runge-kutta-5.lay
fi
cd -
exit 0
