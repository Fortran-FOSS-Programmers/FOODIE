#!/bin/bash
# script automotion for FOODIE's errors analysis

# defaults
clean=0
# schemes enabling flags (1 => enabled, disabled otherwise)
ab=0
abm=0
am=0
bdf=0
euler=0
lf=0
lf_raw=0
lmm_ssp=0
lmm_ssp_vss=0
rk_emd=0
rk_ls=0
rk_ssp=0

function print_usage {
  echo
  echo "`basename $0`"
  echo "script automotion for FOODIE's errors analysis"
  echo "Usage: `basename $0` #scheme [-clean]"
  echo "  -clean         clean old analyis results"
  echo "Valid scheme switches are:"
  echo "  -ab            adams-bashforth"
  echo "  -abm           adams-bashforth-moulton"
  echo "  -am            adams-moulton"
  echo "  -bdf           backward differentiation formula"
  echo "  -euler         euler"
  echo "  -lf            leapfrog"
  echo "  -lf_raw        leapfrog raw filtered"
  echo "  -lmm_ssp       linear multistep methods SSP"
  echo "  -lmm_ssp-vss   linear multistep methods SSP with Variable Step Size"
  echo "  -rk_emd        embeded RK"
  echo "  -rk_ls         low storage RK"
  echo "  -rk_ssp        TVD RK"
  echo "  -all           all schemes!"
}

#parsing command line
while [ $# -gt 0 ]; do
  case "$1" in
    "-clean")
      clean=1
      ;;
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
    "-euler")
      euler=1
      ;;
    "-lf")
      lf=1
      ;;
    "-lf_raw")
      lf_raw=1
      ;;
    "-lmm_ssp")
      lmm_ssp=1
      ;;
    "-lmm_ssp-vss")
      lmm_ssp_vss=1
      ;;
    "-rk_emd")
      rk_emd=1
      ;;
    "-rk_ls")
      rk_ls=1
      ;;
    "-rk_ssp")
      rk_ssp=1
      ;;
    "-all")
      ab=1
      abm=1
      am=1
      bdf=1
      euler=1
      lf=1
      lf_raw=1
      lmm_ssp=1
      lmm_ssp_vss=1
      rk_emd=1
      rk_ls=1
      rk_ssp=1
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
if [ $clean -eq 1 ] ; then
  rm -rf results-oscillation
fi
mkdir -p results-oscillation
cd results-oscillation
if [ $ab -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_1 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >  analysis_adams_bashforth.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_bashforth.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_bashforth.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_bashforth.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_5 -Dt 1250.d0 625.d0 320.d0 100.d0                 --errors_analysis >> analysis_adams_bashforth.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_6 -Dt 1250.d0 625.d0 320.d0 100.d0                 --errors_analysis >> analysis_adams_bashforth.log
  ../tests/oscillation -r -tf 1e5 --scheme adams_bashforth_7 -Dt 10.d0 5.d0 1.d0                              --errors_analysis >> analysis_adams_bashforth.log
  ../tests/oscillation -r -tf 1e4 --scheme adams_bashforth_8 -Dt 1.d0 0.5d0 0.1d0                             --errors_analysis >> analysis_adams_bashforth.log
fi
if [ $abm -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_1 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >  analysis_adams_bashforth_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_bashforth_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_bashforth_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_bashforth_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_5 -Dt 1250.d0 625.d0 320.d0 100.d0                 --errors_analysis >> analysis_adams_bashforth_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_6 -Dt 1250.d0 625.d0 320.d0 100.d0                 --errors_analysis >> analysis_adams_bashforth_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_7 -Dt 400.d0 250.d0 200.d0                         --errors_analysis >> analysis_adams_bashforth_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_bashforth_moulton_8 -Dt 250.d0 200.d0 100.d0                         --errors_analysis >> analysis_adams_bashforth_moulton.log
fi
if [ $am -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_1 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >  analysis_adams_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_adams_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_5 -Dt 1250.d0 625.d0 320.d0 100.d0                 --errors_analysis >> analysis_adams_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_6 -Dt 1250.d0 625.d0 320.d0 100.d0                 --errors_analysis >> analysis_adams_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_7 -Dt 400.d0 250.d0 200.d0                         --errors_analysis >> analysis_adams_moulton.log
  ../tests/oscillation -r -tf 1e6 --scheme adams_moulton_8 -Dt 250.d0 200.d0 100.d0                         --errors_analysis >> analysis_adams_moulton.log
fi
if [ $bdf -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme back_df_1 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >  analysis_back_df.log
  ../tests/oscillation -r -tf 1e6 --scheme back_df_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_back_df.log
  ../tests/oscillation -r -tf 1e6 --scheme back_df_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_back_df.log
  ../tests/oscillation -r -tf 1e6 --scheme back_df_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_back_df.log
  ../tests/oscillation -r -tf 1e6 --scheme back_df_5 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_back_df.log
  ../tests/oscillation -r -tf 1e6 --scheme back_df_6 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_back_df.log
fi
if [ $euler -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme euler_explicit -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis_euler_explicit.log
fi
if [ $lf -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme leapfrog -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis_leapfrog.log
fi
if [ $lf_raw -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme leapfrog_raw -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis > analysis_leapfrog_raw.log
fi
if [ $lmm_ssp -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_steps_3_order_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >  analysis_lmm_ssp.log
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_steps_4_order_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >> analysis_lmm_ssp.log
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_steps_5_order_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >> analysis_lmm_ssp.log
fi
if [ $lmm_ssp_vss -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_vss_steps_2_order_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >  analysis_lmm_ssp_vss.log
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_vss_steps_3_order_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >> analysis_lmm_ssp_vss.log
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_vss_steps_3_order_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >> analysis_lmm_ssp_vss.log
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_vss_steps_4_order_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >> analysis_lmm_ssp_vss.log
  ../tests/oscillation -r -tf 1e6 --scheme lmm_ssp_vss_steps_5_order_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 50.d0 10.d0 --errors_analysis >> analysis_lmm_ssp_vss.log
fi
if [ $rk_emd -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_emd_stages_6_order_5   --tolerance 1e-8 1e-9 --errors_analysis >  analysis_runge_kutta_emd.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_emd_stages_7_order_4   --tolerance 1e-8 1e-9 --errors_analysis >> analysis_runge_kutta_emd.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_emd_stages_9_order_6   --tolerance 1e-8 1e-9 --errors_analysis >> analysis_runge_kutta_emd.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_emd_stages_17_order_10 --tolerance 1e-8 1e-9 --errors_analysis >> analysis_runge_kutta_emd.log
fi
if [ $rk_ls -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ls_stages_1_order_1  -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >  analysis_runge_kutta_ls.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ls_stages_5_order_4  -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ls.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ls_stages_6_order_4  -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ls.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ls_stages_7_order_4  -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ls.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ls_stages_12_order_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ls.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ls_stages_13_order_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ls.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ls_stages_14_order_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ls.log
fi
if [ $rk_ssp -eq 1 ] ; then
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ssp_stages_1_order_1 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >  analysis_runge_kutta_ssp.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ssp_stages_2_order_2 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ssp.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ssp_stages_3_order_3 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ssp.log
  ../tests/oscillation -r -tf 1e6 --scheme runge_kutta_ssp_stages_5_order_4 -Dt 5000.d0 2500.d0 1250.d0 625.d0 320.d0 100.d0 --errors_analysis >> analysis_runge_kutta_ssp.log
fi
cd -
exit 0
