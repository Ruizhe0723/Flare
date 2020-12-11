#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 14:53:00 2020

@author: zc252
"""
# =============================================================================
#
# Imports
#
# =============================================================================

import commonBlock 

import canteraSim 

import interpToMeshgrid

# =============================================================================
#
# Function definitions
#
# =============================================================================

cbDict = commonBlock.load_dict()

meshgrid = commonBlock.create_meshgrid(cbDict)

# =============================================================================
#
# 
#
# =============================================================================

solIdx = 1

canteraSim.adiabaticFlame(solIdx,cbDict)

interpToMeshgrid.interpLamFlame(solIdx,cbDict,meshgrid)

