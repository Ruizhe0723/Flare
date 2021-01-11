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

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Global parameters
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def load_dict():

    chemMech = 'gri30.xml'
    
    commDict =  {
                    'chemMech' : chemMech,
                    'nSpeMech' : ct.Solution(chemMech).n_total_species,
                    'transModel' : 'Mix',
                        
                    'f_min' : 2.85e-02,
                    'f_max' : 8.5e-02,
                    'nchemfile' : 20,
                    'nVarCant' : 10,
                    'cUnifPts' : 1001,
                    'n_points_z' : 501,
                    'n_points_c' : 401,
                    'n_points_h' : 1,
                    
                    'int_pts_z' : 20,
                    'int_pts_c' : 21,
                    'int_pts_gz' : 15,
                    'int_pts_gc' : 11,
                    'int_pts_gcz' : 1,
                    'nScalars' : 11,
                    'nYis' : 1
                }
    
    return commDict


def create_meshgrid(cbDict):
    
    z_space = np.zeros([cbDict['n_points_z'],])
    nn = int((cbDict['n_points_z'] - 1)/5*4 + 1)
    for i in range(len(z_space)):
        if i < nn: 
            z_space[i] = cbDict['f_max']*1.2/(nn-1)*i
        else:
            z_space[i] = cbDict['f_max']*1.2 + (1.-cbDict['f_max']*1.2)/(cbDict['n_points_z']-nn)*(i-nn+1)
    c_space = np.linspace(0.,1.,cbDict['n_points_c'])

    meshgrid = {
                'z_space' : z_space,
                'c_space' : c_space,
           }
    
    return meshgrid

def create_manifold(cbDict):
    
    z_int = np.zeros([cbDict['n_points_z'],])
    nn = int((cbDict['n_points_z'] - 1)/5*4 + 1)
    for i in range(len(z_space)):
        if i < nn: 
            z_space[i] = cbDict['f_max']*1.2/(nn-1)*i
        else:
            z_space[i] = cbDict['f_max']*1.2 + (1.-cbDict['f_max']*1.2)/(cbDict['n_points_z']-nn)*(i-nn+1)
    c_space = np.linspace(0.,1.,cbDict['n_points_c'])

    meshgrid = {
                'z_space' : z_space,
                'c_space' : c_space,
           }
    
    return meshgrid