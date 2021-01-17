#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 12:26:49 2020

@author: zc252
"""

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Imports
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

import cantera as ct
import numpy as np
import os

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Global parameters
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def load_dict():

    chemMech = 'gri30.xml'

    commDict =  {
                     'work_dir' : ('/home/zc252/OpenFOAM_run/' +
                                   'TNF/new-case_reacting_2p6/' +
                                   'flareTable_unscaledPV/'),

                     'output_fln' : 'flare.tbl',

                    'chemMech' : chemMech,
                    'nSpeMech' : ct.Solution(chemMech).n_total_species,
                    'transModel' : 'Mix',
                    'scaled_PV' : False,

                    'f_min' : 2.85e-02,
                    'f_max' : 15.0e-02,
                    'nchemfile' : 50,
                    'nVarCant' : 10,
                    'cUnifPts' : 1001,
                    'n_points_z' : 501,
                    'n_points_c' : 401,
                    'n_points_h' : 1,

                    'int_pts_z' : 80,
                    'int_pts_c' : 51,
                    'int_pts_gz' : 15,
                    'int_pts_gc' : 21,
                    'int_pts_gcz' : 1,
                    'nScalars' : 12,

                    'nYis' : 3,
                    'spc_names' : ['H2O','CO','CO2'],

                    'small' : 1.0e-4,
                    'n_procs' : 16,
                    'exe_path' : '~/Dropbox/Codes/flare/bin/flare'
                }

    return commDict


def create_meshgrid(cbDict):

  z_space = np.linspace(0,1,cbDict['n_points_z'])
  nn = int(cbDict['n_points_z']/20*19)
  for i in range(nn):
      z_space[i] = cbDict['f_max']*1.2/(nn-1)*i
  z_space[nn-1:] = np.linspace(cbDict['f_max']*1.2,1.0,
                               cbDict['n_points_z']-nn+1)
  #
  c_space = np.linspace(0,1,cbDict['n_points_c'])
  #
  meshgrid = { 'z_space' : z_space,
               'c_space' : c_space, }

  return meshgrid

def create_manifold(cbDict):

  work_dir = cbDict['work_dir']

  # assign the manifold control parameters
  z_int = np.linspace(0,1,cbDict['int_pts_z'])
  nn = int(cbDict['int_pts_z']/20*19)
  for i in range(nn):
    z_int[i] = cbDict['f_max']*1.2/(nn-1)*i
  z_int[nn-1:] = np.linspace(cbDict['f_max']*1.2,1.0,
                             cbDict['int_pts_z']-nn+1)
  #
  c_int = np.linspace(0,1,cbDict['int_pts_c'])
  #
  gz_int = np.zeros(cbDict['int_pts_gz'])
  gz_int[1:] = np.logspace(-4,-1,cbDict['int_pts_gz']-1)
  #
  gc_int = np.linspace(0,1-cbDict['small'],cbDict['int_pts_gc'])
  #
  gcz_int = np.linspace(0,0,cbDict['int_pts_gcz'])

  if not os.path.isdir(work_dir): os.mkdir(work_dir)

  # write the manifold control parameters
  fln = work_dir + 'integrate.inp'
  print('Writing integrate.inp ...')
  with open(fln,'w') as strfile:
    #
    strfile.write(str(cbDict['nSpeMech']) + '\n')
    #
    strfile.write('%.5E' % cbDict['f_min'] + '\t' +
                    '%.5E' % cbDict['f_max'] + '\n')
    #
    strfile.write(str(cbDict['nchemfile']) + '\n')
    #
    strfile.write(str(cbDict['n_points_z']) + '\t' +
                  str(cbDict['n_points_c']) + '\t' +
                  str(cbDict['n_points_h']) + '\n')
    #
    strfile.write(str(cbDict['int_pts_z']) + '\t' +
                  str(cbDict['int_pts_c']) + '\t' +
                  str(cbDict['int_pts_gz']) + '\t' +
                  str(cbDict['int_pts_gc']) + '\t' +
                  str(cbDict['int_pts_gcz']) + '\n' )
    #
    np.savetxt(strfile,np.expand_dims(z_int,axis=0),
               fmt='%.5E',delimiter='\t')
    np.savetxt(strfile,np.expand_dims(c_int,axis=0),
               fmt='%.5E',delimiter='\t')
    np.savetxt(strfile,np.expand_dims(gz_int,axis=0),
               fmt='%.5E',delimiter='\t')
    np.savetxt(strfile,np.expand_dims(gc_int,axis=0),
               fmt='%.5E',delimiter='\t')
    np.savetxt(strfile,np.expand_dims(gcz_int,axis=0),
               fmt='%.5E',delimiter='\t')
    #
    strfile.write(str(cbDict['nScalars']) + '\t' +
                  str(cbDict['nYis']) + '\n')
    #
  strfile.close()

  print('\nDone.')

  return