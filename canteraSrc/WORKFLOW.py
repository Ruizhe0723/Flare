#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 14:53:00 2020

@author: zc252
"""

# %% ======================================================================
#
# Imports
#
# =========================================================================

import commonBlock

import canteraSim

import interpToMeshgrid

import pdf

# %% ======================================================================
#
# Load common dictionary
#
# =========================================================================

cbDict = commonBlock.load_dict()

contVarDict = commonBlock.create_manifold(cbDict)

solIdx = 1

# %% ======================================================================
#
# Run cantera simulations
#
# =========================================================================

canteraSim.adiabaticFlame(solIdx,cbDict)

# %% ======================================================================
#
# Interpolate to meshgrid and write first small part of table
#
# =========================================================================

meshgrid = commonBlock.create_meshgrid(cbDict)

d2Yeq_table = interpToMeshgrid.interpLamFlame(solIdx,cbDict,contVarDict,
                                              meshgrid)

# %% ======================================================================
#
# Integrate pdf and write big part of table
#
# =========================================================================

pdf.integrate(cbDict)

pdf.assemble(cbDict,contVarDict,d2Yeq_table)

