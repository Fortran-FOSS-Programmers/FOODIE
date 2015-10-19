#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import argparse
import math
# import os
import numpy as np
import subprocess
import sys


def syswork(cmd):
  """
  Execute system command 'cmd'.
  """
  work_error = 0
  try:
    work_output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
  except subprocess.CalledProcessError as err:
    work_error = err.returncode
    work_output = err.output
  if sys.version_info[0] > 2:
    work_output = work_output.decode("ascii")
  return [work_error, str(work_output)]


def compute_error(s_coarse, s_fine, Dx):
  """
  Compute (relative) error between two solutions.

  Parameters
  ----------
  s_coarse : numpy array
    solution on coarse grid overaged on average grid
  s_fine : numpy array
    solution on fine grid overaged on average grid
  Dx : real
    resolution of average grid

  Returns
  error_estimate : real
    Norm L2 of error estimate
  """
  error_estimate = 0.
  for index, sol_coarse in np.ndenumerate(s_coarse):
    error_estimate += Dx * (s_fine[index] - sol_coarse) ** 2
  error_estimate = math.sqrt(error_estimate)
  return error_estimate


if __name__ == '__main__':
  cliparser = argparse.ArgumentParser(prog='error_analysis.py', description='Perform error analysis of 1D Euler test of FOODIE library')
  cliparser.add_argument('-ex', '--execute', required=False, action='store_true', default=False, help='Execute the euler-1D program (avoided if results files already generated)')
  cliparser.add_argument('-p', '--problem', required=False, action='store', default='sod', choices=['sod', 'smooth'], help='Problem to be solved')
  cliparser.add_argument('--Ni_grids', required=False, action='store', nargs='+', default=[100, 200, 400, 800, 1600, 3200], help='Grids number of cells')
  cliparser.add_argument('--av_Ni', required=False, action='store', default=100, help='Number of cells of grid over averages are performed')
  cliparser.add_argument('-s', '--solver', required=False, action='store', default='ls-runge-kutta',
                         choices=['adams-bashforth', 'euler', 'leapfrog', 'ls-runge-kutta', 'tvd-runge-kutta'], help='ODE solver')
  cliparser.add_argument('--ss', required=False, action='store', default=5, help='Stages/steps used')
  cliargs = cliparser.parse_args()
  Ni_grids = cliargs.Ni_grids
  av_Ni = cliargs.av_Ni
  Dx_grids = [1. / float(Ni) for Ni in Ni_grids]
  av_Dx = 1. / float(av_Ni)
  rho = []
  for Ni in Ni_grids:
    if cliargs.execute:
      [error, output] = syswork('./euler-1D --solver ' + cliargs.solver +
                                ' --ss ' + str(cliargs.ss) + ' -r --problem ' + cliargs.problem +
                                ' --av_Ni ' + str(cliargs.av_Ni) + ' --Ni ' + str(Ni) +
                                ' --output error_analysis')
      if error != 0:
        sys.stderr.write(output)
    outputfile = 'error_analysis-' + cliargs.solver
    if cliargs.solver in ['adams-bashforth', 'ls-runge-kutta', 'tvd-runge-kutta']:
      outputfile += '-' + str(cliargs.ss)
    outputfile += '-' + str(Ni) + '_cells.dat'
    with open(outputfile) as foodie_results:
      rho_av = np.fromstring(foodie_results.readlines()[-2], dtype=np.float64, sep=' ')
      rho.append(rho_av)
  errors = []
  for index in range(len(Ni_grids) - 1):
    errors.append(compute_error(s_coarse=rho[index], s_fine=rho[-1], Dx=av_Dx))
  observed_orders = []
  for index in range(len(errors) - 1):
    observed_orders.append(math.log10(errors[index] / errors[index + 1]) / math.log10(Dx_grids[index] / Dx_grids[index + 1]))
  observed_order = 0.
  for index in range(len(errors) - 1):
    print("Ni: " + str(Ni_grids[index]) + " Dx: " + str(Dx_grids[index]) + " error: " + str(errors[index]) + " order: " + str(observed_orders[index]))
    observed_order += observed_orders[index]
  observed_order = observed_order / max(1, (len(errors) - 1))
  print("Average observed order: " + str(observed_order))
