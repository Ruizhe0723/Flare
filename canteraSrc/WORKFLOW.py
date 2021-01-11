#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 14:53:00 2020

@author: zc252
"""

# %% ==========================================================================
#
# Imports
#
# =============================================================================

import commonBlock 

import canteraSim 

import interpToMeshgrid

import os

#%% ==========================================================================
#
# Function definitions
#
# =============================================================================

cbDict = commonBlock.load_dict()

meshgrid = commonBlock.create_meshgrid(cbDict)

# %% ==========================================================================
#
# 
#
# =============================================================================

solIdx = 1

work_dir = ('/home/zc252/OpenFOAM/zc252-7/run/TNF/new-case_reacting_2p6/' 
            + 'flareTable/')
if not os.path.isdir(work_dir): os.mkdir(work_dir)
        
canteraSim.adiabaticFlame(work_dir,solIdx,cbDict)

interpToMeshgrid.interpLamFlame(work_dir,solIdx,cbDict,meshgrid)

