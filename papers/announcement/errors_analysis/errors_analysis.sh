#!/bin/bash
# script automotion for FOODiE's errors analysis of the provided tests
cd ../../../
echo "Building tests executables"
FoBiS.py build -mode error-analysis --build_dir papers/announcement/error_analysis/tests > buils.log
cd -
echo "Oscillation test analysis"
rm -rf results-oscillation
echo "  Analize leapfrog solver"
./utilities-oscillation/error_analysis.py -Dts 5000 2500 1250 625 320 100 -ex -s leapfrog > leapfrog.log
echo "  Analize Adams-Bashforth solver"
./utilities-oscillation/error_analysis.py -Dts 5000 2500 1250 625 320 100 -ex -s adams-bashforth --ss 2 > adams-bashforth-2.log
./utilities-oscillation/error_analysis.py -Dts 5000 2500 1250 625 320 100 -ex -s adams-bashforth --ss 3 > adams-bashforth-3.log
echo "  Analize low storage Runge-Kutta solver"
./utilities-oscillation/error_analysis.py -Dts 5000 2500 1250 625 320 100 -ex -s ls-runge-kutta --ss 5 > ls-runge-kutta-5.log
echo "  Analize TVD/SSP Runge-Kutta solver"
./utilities-oscillation/error_analysis.py -Dts 5000 2500 1250 625 320 100 -ex -s tvd-runge-kutta --ss 2 > tvd-runge-kutta-2.log
./utilities-oscillation/error_analysis.py -Dts 5000 2500 1250 625 320 100 -ex -s tvd-runge-kutta --ss 3 > tvd-runge-kutta-3.log
./utilities-oscillation/error_analysis.py -Dts 5000 2500 1250 625 320 100 -ex -s tvd-runge-kutta --ss 5 > tvd-runge-kutta-5.log
mkdir -p results-oscillation
mv *log *dat results-oscillation
cd results-oscillation
cp ../utilities-oscillation/*lay .
cp ../utilities/lay_export* .
lay_export_all_f error_analysis-oscillation-adams-bashforth-2.lay
lay_export_all_f error_analysis-oscillation-adams-bashforth-3.lay
lay_export_all_f error_analysis-oscillation-leapfrog.lay
lay_export_all_f error_analysis-oscillation-ls-runge-kutta-5.lay
lay_export_all_f error_analysis-oscillation-tvd-runge-kutta-2.lay
lay_export_all_f error_analysis-oscillation-tvd-runge-kutta-3.lay
lay_export_all_f error_analysis-oscillation-tvd-runge-kutta-5.lay
cd -
exit 0
