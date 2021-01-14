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
def integrate(work_dir,exe_path,n_procs):

  # execute the Fortran code
  print("\n Changing directory to:")
  os.chdir(work_dir)
  os.system("pwd")
  os.system("mpirun -np " + str(n_procs) + " " + exe_path)


''' ===========================================================================

Assemble integrated MPI unit files

=========================================================================== '''
def assemble(work_dir,cbDict):

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

  # remove unwanted columns - last two for qdot and h / Yis
  MM = np.delete(M,[12,13],axis=1)
  if cbDict['nYis'] > 0:
    ind = -cbDict['nYis']
    MM = MM[:,:ind]
    YM = M[:,ind:]

  # write assembled table
  fln = work_dir + 'flare.tbl'
  print('Writing assembled table ...')
  with open(fln,'a') as strfile:
    #
    strfile.write(str(cbDict['int_pts_z']) + '\t' +
                  str(cbDict['int_pts_c']) + '\t' +
                  str(cbDict['int_pts_gz']) + '\t' +
                  str(cbDict['int_pts_gc']) + '\t' +
                  str(cbDict['int_pts_gcz']) + '\n' )
    #
    np.savetxt(strfile,MM[:,0:12],fmt='%.5E',delimiter='\t')
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


