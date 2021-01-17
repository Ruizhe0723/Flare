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

#%% ==========================================================================
#
# Initial setup
#
# =============================================================================

cbDict = commonBlock.load_dict()

meshgrid = commonBlock.create_meshgrid(cbDict)

commonBlock.create_manifold(cbDict)

# %% ==========================================================================
#
#
#
# =============================================================================

solIdx = 1

canteraSim.adiabaticFlame(solIdx,cbDict)

interpToMeshgrid.interpLamFlame(solIdx,cbDict,meshgrid)

pdf.integrate(cbDict)

pdf.assemble(cbDict)

