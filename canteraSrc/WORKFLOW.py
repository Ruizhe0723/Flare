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

import pdf

import os

#%% ==========================================================================
#
# Initial setup
#
# =============================================================================

cbDict = commonBlock.load_dict()

meshgrid = commonBlock.create_meshgrid(cbDict)

work_dir = ('/home/zc252/OpenFOAM_run/TNF/new-case_reacting_2p6/'
            + 'flareTable/')
if not os.path.isdir(work_dir): os.mkdir(work_dir)

commonBlock.create_manifold(work_dir,cbDict)

# %% ==========================================================================
#
#
#
# =============================================================================

solIdx = 1

canteraSim.adiabaticFlame(work_dir,solIdx,cbDict)

interpToMeshgrid.interpLamFlame(work_dir,solIdx,cbDict,meshgrid)

n_procs = 16
exe_path = '/home/zc252/Dropbox/Codes/flare/bin/flare'
pdf.integrate(work_dir,exe_path,n_procs)

pdf.assemble(work_dir,cbDict)

