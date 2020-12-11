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
#                    'solIdx' : 1,
                    'cUnifPts' : 1001,
                    # 'nScalars' : 7,
                    'n_points_z' : 501,
                    'n_points_c' : 401,
                }
    
    return commDict


def create_meshgrid(cbDict):
    
    Z_space = np.zeros([cbDict['n_points_z'],])
    nn = int((cbDict['n_points_z'] - 1)/5*4 + 1)
    for i in range(len(Z_space)):
        if i < nn: 
            Z_space[i] = cbDict['f_max']*1.2/(nn-1)*i
        else:
            Z_space[i] = cbDict['f_max']*1.2 + (1.-cbDict['f_max']*1.2)/(cbDict['n_points_z']-nn)*(i-nn+1)
    c_space = np.linspace(0.,1.,cbDict['n_points_c'])

    meshgrid = {
                'Z_space' : Z_space,
                'c_space' : c_space,
           }
    
    return meshgrid
