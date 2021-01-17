#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 10:51:45 2021

@author: zc252
"""

# =============================================================================
#
# Imports
#
# =============================================================================

import numpy as np

import os

# =============================================================================
#
# Function definitions
#
# =============================================================================

''' ===========================================================================

Integrate meshgrid laminar data with PDFs

=========================================================================== '''
def integrate(cbDict):

  work_dir = cbDict['work_dir']

  cwd = os.getcwd()

  # set up environment
  print("\n Changing directory to:")
  os.chdir(work_dir)
  os.system("pwd")
  #
  os.system("pkill 'flare'")
  #
  if os.path.isfile('unit01_h01.dat'):
    os.system("rm unit*")

  # execute the Fortran code
  os.system("mpirun -np " + str(cbDict['n_procs']) + " "
            + cbDict['exe_path'])

  # change directory back
  print("\n Changing directory to:")
  os.chdir(cwd)
  os.system("pwd")

''' ===========================================================================

Assemble integrated MPI unit files

=========================================================================== '''
def assemble(cbDict):

  work_dir = cbDict['work_dir']

  # read integrated units
  fln = work_dir
  #
  for i in range(1,cbDict['n_points_h']+1):
    #
    print('Reading unit' + '%02d' % 1 + '_h' + '%02d' % i + ' ... \n')
    M = np.loadtxt(fln + 'unit01_h' + '%02d' % i + '.dat')
    #
    for j in range(2,cbDict['int_pts_z']+1):
      #
      print('Reading unit' + '%02d' % j + '_h' + '%02d' % i + ' ... \n')
      #
      tmp = np.loadtxt(fln + 'unit' + '%02d' % j + '_h' + '%02d' % i
                     + '.dat')
      #
      M = np.insert(tmp,0,M,axis=0)

  # remove unwanted columns - h(12),qdot(13),Yc_max(14)
  if(cbDict['scaled_PV']):
    rm_list = [12,13,14]
  else:
    rm_list = [12,13]
  MM = np.delete(M,rm_list,axis=1)

  # separate scalars and Yis
  if cbDict['nYis'] > 0:
    ind = -cbDict['nYis']
    MS = MM[:,:ind]
    YM = M[:,ind:]

  # write assembled table
  fln = work_dir + cbDict['output_fln']
  print('Writing assembled table ...')
  with open(fln,'a') as strfile:
    #
    strfile.write(str(cbDict['int_pts_z']) + '\t' +
                  str(cbDict['int_pts_c']) + '\t' +
                  str(cbDict['int_pts_gz']) + '\t' +
                  str(cbDict['int_pts_gc']) + '\t' +
                  str(cbDict['int_pts_gcz']) + '\t' +
                  str(MS.shape[1]) + '\n' )
    #
    np.savetxt(strfile,MS,fmt='%.5E',delimiter='\t')
    #
  strfile.close()

  #append Yis to the table
  if cbDict['nYis'] > 0:
    #
    with open(fln,'a') as strfile:
      #
      strfile.write(str(cbDict['nYis']) + '\n')
      #
      np.savetxt(strfile,YM,fmt='%.5E',delimiter='\t')
      #
    strfile.close()


  print('\nDone.')


