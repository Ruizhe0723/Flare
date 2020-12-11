#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 12:04:25 2020

@author: zc252
"""

# =============================================================================
#
# Imports
#
# =============================================================================

import cantera as ct
import numpy as np
import sys
# =============================================================================
#
# Function definitions
#
# =============================================================================

def adiabaticFlame(solIdx,cbDict):
    solFln = ('canteraData/solution_' 
              + str('{:02d}'.format(solIdx-1)) + '/')
    CASENAME = 'CH4' # case name
    p =  101325  # pressure [Pa]
    Lx = 0.05 # Domain size for the simulation [m] - has to be larger than flame thickness
    chemMech = cbDict['chemMech'] # chemical mechanism
    
    ## Fuel characteristics
    fuel_species = 'CH4' # Fuel is assumed to be of the form CxHy
    fuel_C = 1. # number of C atoms in the fuel
    fuel_H = 4. # number of H atoms in the fuel
    stoich_O2 = fuel_C+fuel_H/4. # DO NOT CHANGE - stoichiometric air mole fraction
    
    W_fuel = fuel_C * 12. + fuel_H * 1.0 # DO NOT CHANGE - fuel molar weight
    
    T_fuel = 200. # Fuel temperature [K]
    X_fuel = 'CH4:1' # Fuel composition (in mole fraction)
    
    ## Oxidiser characteristics
    W_O2 = 2. * 16. # DO NOT CHANGE - molar weight of O2
    W_N2 = 2. * 14. # DO NOT CHANGE - molar weight of N2
    
    T_ox = 333. # oxidiser temperature [K]
    X_ox = 'O2:0.21, N2:0.79' # oxidiser composition (in mole fraction)
    
    ## Mixture properties
    Zst = (W_fuel) / (W_fuel + stoich_O2 * ( W_O2 + 3.76 * W_N2) ) # DO NOT CHANGE - stoichiometric mixture fraction
    Z = np.linspace(cbDict['f_min'],cbDict['f_max'],cbDict['nchemfile']) # DO NOT CHANGE - array of mixture fraction of interest
    
    # DO NOT CHANGE BELOW THIS LINE
    phi = Z*(1.0 - Zst) / (Zst*(1.0 - Z))
    phi_tab = np.zeros((len(phi),8))
    sL = []
    
    for i in range(len(phi)):
    #  Z = (phi[i]*Zst)/(1-Zst+Zst*phi[i])
    
      reactants = {fuel_species: phi[i] / stoich_O2, 'O2': 1.0, 'N2': 3.76}
    
      ## Load chemical mechanism
      gas = ct.Solution(chemMech)
      
      # Stream A (air)
      A = ct.Quantity(gas, constant='HP')
      A.TPX = T_ox, p, X_ox # Define oxidiser
      
      # Stream B (methane)
      B = ct.Quantity(gas, constant='HP')
      B.TPX = T_fuel, p, X_fuel # Define fuel
    
      # Set the molar flow rates corresponding to stoichiometric reaction,
      # CH4 + 2 O2 -> CO2 + 2 H2O
      A.moles = 1
      nO2 = A.X[A.species_index('O2')]
      B.moles = nO2 / stoich_O2 * phi[i]
      
      # Compute the mixed state
      M = A + B
    
      # unburned gas temperature [K]
      Tin = M.T
      print('Flamelet No. ',i,'|',phi[i],Z[i],Tin,p,reactants)
    
      # set reactants state  
      gas.TPX = Tin, p, reactants
    
      # Flame object
      f = ct.FreeFlame(gas,width=Lx)
    
      transModel = cbDict['transModel']
      if transModel == 'Mix' or transModel == 'Multi':
          # Solve with the energy equation disabled
          f.energy_enabled = False
          f.transport_model = transModel
          f.solve(loglevel=1, refine_grid=False)
        
          # Solve with the energy equation enabled
          f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14) 
          f.energy_enabled = True
          f.solve(loglevel=1,auto=True)
          print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))
      
      elif transModel == 'Multi':
          # Solve with multi-component transport properties
          f.transport_model = 'Multi'
          f.solve(loglevel=1, auto=True)
          print('multicomponent flamespeed = {0:7f} m/s'.format(f.u[0]))
        
      else: 
          sys.exit(' !! Error - incorrect transport model specified !!')
          
      sL.append(f.u[0]) # store flame speed in array sL
      
      # store useful data for future simulation
      iCO = gas.species_index('CO')
      iCO2 = gas.species_index('CO2')
      data = np.zeros((len(f.grid),cbDict['nSpeMech']+cbDict['nVarCant']))
      
      # unscaled progress variable
      data[:,0] = f.Y[iCO] + f.Y[iCO2]
      
      # Reaction rate of unscaled progress variable
      data[:,1] = f.gas.molecular_weights[iCO]*f.net_production_rates[iCO,:] \
          + f.gas.molecular_weights[iCO2]*f.net_production_rates[iCO2,:]
    
      data[:,2] = f.T
      data[:,3] = f.density_mass
      data[:,4] = f.cp_mass
      data[:,5] = np.dot(np.transpose(f.X),f.gas.molecular_weights) 
    
      # formation enthalpy
      for j in range(len(f.grid)):
          dumGas = ct.Solution(chemMech) # dummy working variable
          dumGas.TPY = 298.15,p,f.Y[:,j]
          data[j,6] = dumGas.enthalpy_mass
    
      data[:,7] = f.viscosity/f.density
      data[:,8] = f.enthalpy_mass
      data[:,9] = f.heat_release_rate
    
      data[:,cbDict['nVarCant']:
             cbDict['nSpeMech']+cbDict['nVarCant']] = np.transpose(f.Y) # store species mass fractions
    
      # save flamelet data
      fln = solFln + CASENAME + '_' + '{:03d}'.format(i) + '.csv'
      np.savetxt(fln,data)
      
      ###### Compute global flame properties #######
      ###### calculate flame thickness
      DT = np.gradient(f.T,f.grid)
      dl = (f.T[-1]-f.T[0])/max(DT)  
     
      phi_tab[i,0] = phi[i]                      #equivalence ratio
      phi_tab[i,1] = Z[i]                        #mixture fraction
      phi_tab[i,2] = len(f.grid)                 #number of grid points
      phi_tab[i,3] = f.u[0]                      #flame speed
      phi_tab[i,4] = dl                          #flame thickness
      phi_tab[i,5] = (f.T[-1]-f.T[0])/f.T[0]     #heat release parameter
    
      c = (f.Y[iCO,:] + f.Y[iCO2,:]) / max(f.Y[iCO,:] + f.Y[iCO2,:]) 
      alpha = f.thermal_conductivity/f.density/f.cp_mass
      Dc = np.gradient(c,f.grid)
      Nc = alpha*Dc*Dc
      PDF_c = Dc*f.viscosity/f.density/f.u[0]  
      integ_1 = np.trapz(f.density*Nc*np.gradient(f.u,f.grid)*PDF_c,c)
      integ_2 = np.trapz(f.density*Nc*PDF_c,c)  
      phi_tab[i,6] = dl/f.u[0] *integ_1/integ_2/phi_tab[i,5]   #KcStar  
    
      ###### calculate integral of cp in T space from 298.15 to Tin
      gasCP = ct.Solution(chemMech)
      gasCP.TPX = 298.15,p,reactants
      cp_0 = gasCP.cp_mass
      if abs(Tin - 298.15) < 0.1:
        phi_tab[i,7] = 0.
      else:
        sum_CpdT = 0.
        dT = (Tin-298.15)/(int(100*(Tin-298.15))-1)
        for kk in range(1,int(100*abs(Tin-298.15))):
          gasCP.TPX = (298.15 + kk*dT),p,reactants
          cp_1 = gasCP.cp_mass
          sum_CpdT = sum_CpdT + 0.5*(cp_0 + cp_1)*dT
          cp_0 = cp_1 
        phi_tab[i,7] = sum_CpdT
    
    ###### calculate boundary conditions for pure fuel and oxidiser
    BCdata = np.zeros((2,7))
    gas_fuel = ct.Solution(chemMech,'gri30_multi')
    gas_fuel.TPX = T_fuel,p,X_fuel
    BCdata[0,0] = gas_fuel.T
    BCdata[0,1] = gas_fuel.density_mass
    
    gas_fuelCP = ct.Solution(chemMech)
    gas_fuelCP.TPX = 298.15,p,X_fuel
    cp_0 = gas_fuelCP.cp_mass
    if abs(T_fuel - 298.15) < 0.1:
      BCdata[0,2] = cp_0
    else:
      sum_CpdT = 0.
      dT = (T_fuel-298.15)/(int(100*(T_fuel-298.15))-1)
      for kk in range(1,int(100*abs(T_fuel-298.15))):
        gas_fuelCP.TPX = (298.15 + kk*dT),p,X_fuel
        cp_1 = gas_fuelCP.cp_mass
        sum_CpdT = sum_CpdT + 0.5*(cp_0 + cp_1)*dT
        cp_0 = cp_1
      BCdata[0,2] = sum_CpdT / abs(T_fuel-298.15)
      
    BCdata[0,3] = np.dot(np.transpose(gas_fuel.X),gas_fuel.molecular_weights)
    
    gas_fuelHf = ct.Solution(chemMech)
    gas_fuelHf.TPX = 298.15,p,X_fuel
    BCdata[0,4] = gas_fuelHf.enthalpy_mass
    
    BCdata[0,5] = gas_fuel.viscosity/gas_fuel.density_mass
    BCdata[0,6] = gas_fuel.enthalpy_mass
    
    print('T_fuel_approx: {0:7f}'.format((BCdata[0,6]-BCdata[0,4])/BCdata[0,2]+298.15))
    
    gas_ox = ct.Solution(chemMech,'gri30_multi')
    gas_ox.TPX = T_ox,p,X_ox
    BCdata[1,0] = gas_ox.T
    BCdata[1,1] = gas_ox.density_mass
    
    gas_oxCP = ct.Solution(chemMech)
    gas_oxCP.TPX = 298.15,p,X_ox
    cp_0 = gas_oxCP.cp_mass
    if abs(T_ox - 298.15) < 0.1:
      BCdata[1,2] = cp_0
    else:
      sum_CpdT = 0.0
      dT = (T_ox-298.15)/(int(100*(T_ox-298.15))-1)
      for kk in range(1,int(100*abs(T_ox-298.15))):
        gas_oxCP.TPX = (298.15 + kk*dT),p,X_ox
        cp_1 = gas_oxCP.cp_mass
        sum_CpdT = sum_CpdT + 0.5*(cp_0 + cp_1)*dT
        cp_0 = cp_1
      BCdata[1,2] = sum_CpdT / abs(T_ox-298.15)
      
    BCdata[1,3] = np.dot(np.transpose(gas_ox.X),gas_ox.molecular_weights)
    
    gas_oxHf = ct.Solution(chemMech)
    gas_oxHf.TPX = 298.15,p,X_ox
    BCdata[1,4] = gas_oxHf.enthalpy_mass
    
    BCdata[1,5] = gas_ox.viscosity
    BCdata[1,6] = gas_ox.enthalpy_mass/gas_ox.density_mass
    
    print('T_ox_approx: {0:7f}'.format((BCdata[1,6]-BCdata[1,4])/BCdata[1,2]+298.15))
    
    # save the laminar parameters of all the flamelets 
    fln_phi_tab = solFln + 'lamParameters.txt'
    with open(fln_phi_tab,'w') as strfile:
      strfile.write(CASENAME + '\n')
      np.savetxt(strfile,phi_tab,fmt='%.5e %.5e %04d %.5e %.5e %.5e %.5e %.5e')
      np.savetxt(strfile,BCdata,fmt='%.5e %.5e %.5e %.5e %.5e %.5e %.5e')
    strfile.close()
