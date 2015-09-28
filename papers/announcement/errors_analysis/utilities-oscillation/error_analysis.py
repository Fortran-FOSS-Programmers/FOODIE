#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import argparse
import numpy as np
import re
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


def compute_error(exact, numerical):
  """
  Compute error of numerical solution once the exact one is provided.

  Parameters
  ----------
  exact : numpy array
    exact solution
  numerical : numpy array
    numerical solution

  Returns
  error_estimate : real
    Norm L2 of error estimate
  """
  error_estimate = 0.
  for i in range(numerical.shape[0]):
    error_estimate += (numerical[i][0] - exact[i][0])**2 + (numerical[i][1] - exact[i][1])**2
  return np.sqrt(error_estimate)


def load_oscillation_results(filename):
  """Load Osciallation test results.

  Parameters
  ----------
  filename : str
    output file name

  Returns
  -------
  numerical : list
    list of numerical results
  exact : list
    list of exact results
  """
  numerical = []
  exact = []
  with open(filename) as foodie_results:
    zones = re.split('ZONE.*\n', foodie_results.read())
    sol_steps = zones[1].split('\n')
    for sol in sol_steps:
      if sol.strip() != '':
        variables = re.split(' *', sol.strip())
        numerical.append([variables[1], variables[2]])
    sol_steps = zones[2].split('\n')
    for sol in sol_steps:
      if sol.strip() != '':
        variables = re.split(' *', sol.strip())
        exact.append([variables[1], variables[2]])
  return numerical, exact


def collect_errors(cliargs):
  """Collect the errors on all simulations performed.

  Parameters
  ----------
  cliargs : argpargse parsed arguments

  Returns
  -------
  errors : list
    list of float64 errors
  """
  time_steps = cliargs.time_steps
  t_final = cliargs.t_final
  errors = []
  for Dt in time_steps:
    Dt = np.float64(Dt)
    steps = int(t_final / Dt)
    if cliargs.execute:
      command = './tests/oscillation --solver ' + cliargs.solver + ' --ss ' + str(cliargs.ss) + ' -r -Dt ' + str(Dt) + ' -tf ' + str(t_final) + ' --output error_analysis-oscillation'
      [error, output] = syswork(command)
      if error != 0:
        sys.stderr.write(command)
        sys.stderr.write(output)
    outputfile = 'error_analysis-oscillation-' + "{0:010d}".format(steps) + '-time_steps-' + cliargs.solver
    if cliargs.solver in ['adams-bashforth', 'ls-runge-kutta', 'tvd-runge-kutta']:
      outputfile += '-' + str(cliargs.ss)
    outputfile += '.dat'
    sol_numerical, sol_exact = load_oscillation_results(filename=outputfile)
    sol_numerical_np = np.zeros((len(sol_numerical), 2), dtype='float64')
    for t, sol in enumerate(sol_numerical):
      sol_numerical_np[t][0] = np.float64(sol[0])
      sol_numerical_np[t][1] = np.float64(sol[1])
    sol_exact_np = np.zeros((len(sol_numerical), 2), dtype='float64')
    for t, sol in enumerate(sol_exact):
      sol_exact_np[t][0] = np.float64(sol[0])
      sol_exact_np[t][1] = np.float64(sol[1])
    errors.append(compute_error(exact=sol_exact_np, numerical=sol_numerical_np))
  return errors


def estimate_orders(time_steps, errors):
  """Estimate the observed order of each integration pair (fine/coarse time step).

  Parameters
  ----------
  time_steps : list
    list of time steps used
  errors: list
    list of errors associated to each integration

  Returns
  -------
  observed_orders : list
    list of orders estimations
  """
  orders = []
  for index in range(len(errors) - 1):
    orders.append(np.log10(errors[index] / errors[index + 1]) /
                  np.log10(np.float64(time_steps[index]) / np.float64(time_steps[index + 1])))
  return orders


def main():
  cliparser = argparse.ArgumentParser(prog='error_analysis.py', description='Perform error analysis of Oscillation test of FOODiE library')
  cliparser.add_argument('-ex', '--execute', required=False, action='store_true', default=False, help='Execute the oscillation program (avoided if results files already generated)')
  cliparser.add_argument('-s', '--solver', required=False, action='store', default='ls-runge-kutta',
                         choices=['adams-bashforth', 'euler', 'leapfrog', 'ls-runge-kutta', 'tvd-runge-kutta'], help='ODE solver')
  cliparser.add_argument('--ss', required=False, action='store', default=5, help='Stages/steps used')
  cliparser.add_argument('-Dts', '--time_steps', required=False, action='store', nargs='+', default=[100, 200, 400, 800, 1600], help='Time steps used')
  cliparser.add_argument('-tf', '--t_final', required=False, default=1000000., action='store', help='Final integration time')
  cliargs = cliparser.parse_args()
  errors = collect_errors(cliargs=cliargs)
  observed_orders = estimate_orders(time_steps=cliargs.time_steps, errors=errors)
  observed_order = 0.
  for index in range(len(errors)):
    if index < len(errors) - 1:
      print("Dt: " + cliargs.time_steps[index] + " error: " + str(errors[index]) + " order: " + str(observed_orders[index]))
      observed_order += observed_orders[index]
    else:
      print("Dt: " + cliargs.time_steps[index] + " error: " + str(errors[index]))
  observed_order = observed_order / max(1, (len(errors) - 1))
  print("Average observed order: " + str(observed_order))

if __name__ == '__main__':
  main()
