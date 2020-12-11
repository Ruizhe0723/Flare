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

import numpy as np
import numpy.matlib 
#from scipy.interpolate import griddata
from scipy import interpolate
#from scipy.interpolate import RectBivariateSpline
# =============================================================================
#
# Function definitions
#
# =============================================================================

''' ===========================================================================

Load Cantera data and interpolate to meshgrid 

=========================================================================== '''
def interpLamFlame(solIdx,cbDict,meshgrid):

    # read laminar flame global parameters & write Table_5 
    solFln = ('canteraData/solution_' 
              + str('{:02d}'.format(solIdx-1)) + '/')
    lamArr = []
    with open(solFln + 'lamParameters.txt') as f:
        ll = 0
        for line in f:
            if ll == 0: 
                casename = line.strip()
            elif ll > cbDict['nchemfile']:
                break
            else:
                line = line.strip()
                lamArr.append(line.split(' '))
            ll += 1
    lamArr = np.array(lamArr,dtype='float')    

    streamBC = np.loadtxt(solFln + 'lamParameters.txt', 
                          skiprows=1+cbDict['nchemfile'])
    
    np.savetxt('Table_5_lamParameters_' 
               + str('{:02d}'.format(solIdx)) + '.dat', 
               np.transpose([lamArr[:,1],lamArr[:,3],lamArr[:,4],
                             lamArr[:,5],lamArr[:,6]]),
               fmt='%.5e %.5e %.5e %.5e %.5e'
              )
    
    # read cantera solutions & calculate c,omega_c
    nScalCant = cbDict['nSpeMech']+cbDict['nVarCant']
    mainData = np.zeros((cbDict['nchemfile'],int(max(lamArr[:,2])),nScalCant))
    cIn = np.zeros(np.shape(mainData[:,:,0]))
    omega_cIn = np.zeros(np.shape(mainData[:,:,0]))
    Yc_eq = np.zeros((cbDict['nchemfile']))
    for i in range(cbDict['nchemfile']):
        fln = (solFln + casename + '_' + str('{:03d}'.format(i)) + '.csv')
        print('\nReading --> ' + fln)
        
        len_grid = int(lamArr[i,2])
        with open(fln) as f:
            j = 0
            for line in f:
                if j >= len_grid:
                    break
                line = line.strip()
                mainData[i,j,:] = line.split(' ')
                j += 1
                
            # for j in range(int(max(lamArr[:,2]))):
            #         if j < int(lamArr[i,2]):
            #             mainData[i,j,:] = np.loadtxt(fln,skiprows=j,max_rows=1)
            #             end = j
            # else: mainData[i,j,:] = np.loadtxt(fln,skiprows=end,max_rows=1)
            
        imax = np.argmax(mainData[i,:,0])
        cIn[i,:] = mainData[i,:,0] / mainData[i,imax,0]
        if mainData[i,imax,0]/mainData[i,len_grid-1,0] > 1.0:
            print('c_max/c_end =',mainData[i,imax,0]/mainData[i,len_grid-1,0],
                  ' --> overshooting')
        
        omega_cIn[i,:] = mainData[i,:,1]/ mainData[i,imax,0]
        
        # Yc_eq[i] = mainData[i,-1,0]
        # # calculate d2Yc_eq/dZ2 & write Table 2
        # generateTable2(lamArr[:,1],Yc_eq)
                
    # interpolate in c space & write out for each flamelet               
    MatScl_c = np.zeros((cbDict['nchemfile'],cbDict['cUnifPts'],nScalCant))
    for i in range(cbDict['nchemfile']):
        len_grid = int(lamArr[i,2])
        ctrim = cIn[i,:len_grid-1]
        MatScl_c[i,:,0] = np.linspace(0.,1.,cbDict['cUnifPts']) 
        MatScl_c[i,:,1] = np.matlib.interp(MatScl_c[i,:,0],ctrim, 
                                         omega_cIn[i,:len_grid-1])
        for k in range(2,nScalCant):
            MatScl_c[i,:,k] = np.matlib.interp(MatScl_c[i,:,0],ctrim,
                                             mainData[i,:len_grid-1,k])
        MatScl_c[i,:,4] = calculateCp_eff(MatScl_c[i,:,2],MatScl_c[i,:,4],
                                        lamArr[i,7])
        fln = (solFln+'Unf_'+casename+'_'+str('{:03d}'.format(i))+'.dat')
        print('\nWriting --> ' + fln)
        np.savetxt(fln,MatScl_c[i,:,:],fmt='%.5e')
        
    oriSclMat = np.zeros([cbDict['nchemfile']+2,cbDict['cUnifPts'],nScalCant])
    for i in range(len(streamBC[:,0])):
        if i == 0: j = len(oriSclMat[:,0,0]) - 1
        else: j = 0
        for k in range(len(streamBC[0,:])):
            oriSclMat[j,:,k+2] = streamBC[i,k] 
        # INSERT Yi BCs HERE
    oriSclMat[1:cbDict['nchemfile']+1,:,:] = MatScl_c
    
    intpSclMat = np.zeros([len(meshgrid['Z_space']),len(meshgrid['c_space']),
                           cbDict['nVarCant']+1])
    Yi_vals = np.zeros([len(meshgrid['Z_space']),len(meshgrid['c_space']),
                        cbDict['nSpeMech']])
    Z_pts = np.insert(lamArr[:,1],0,[0.],axis=0)
    Z_pts = np.insert(Z_pts,len(Z_pts),[1.],axis=0)
    c_pts = MatScl_c[0,:,0]
    print('\nInterpolating...')
    intpSclMat[:,:,0] = interp2D(oriSclMat[:,:,3],Z_pts,c_pts,meshgrid) # rho
    for k in [1,4,5,6]: # omega_c,cp_e,mw,hf
        intpSclMat[:,:,k] = interp2D(oriSclMat[:,:,k],Z_pts,c_pts,meshgrid)
    intpSclMat[:,:,7] = interp2D(oriSclMat[:,:,2],Z_pts,c_pts,meshgrid) # T
    intpSclMat[:,:,8] = interp2D(oriSclMat[:,:,7],Z_pts,c_pts,meshgrid) # nu
    intpSclMat[:,:,9] = interp2D(oriSclMat[:,:,8],Z_pts,c_pts,meshgrid) # h
    intpSclMat[:,:,10] = interp2D(oriSclMat[:,:,9],Z_pts,c_pts,meshgrid) # qdot
    for y in range(cbDict['nSpeMech']):
        Yi_vals[:,:,y] = interp2D(oriSclMat[:,:,y+cbDict['nVarCant']],
                                            Z_pts,c_pts,meshgrid)
    print('\nInterpolation done. ')
    if np.sum(np.isnan(intpSclMat)) > 0: 
        print('\nNumber of Nans detected: ', np.sum(np.isnan(intpSclMat)))
    else: print('\nNo Nans detected. Well done!')
    
    print('\nwriting chemTab file...')
    fln = './chemTab_' + str('{:02d}'.format(solIdx)) + '.dat'
    Arr_c,Arr_Z = np.meshgrid(meshgrid['c_space'],meshgrid['Z_space'])
    Idx_outLmt = np.hstack([(Arr_Z>cbDict['f_max']).nonzero(),
                             (Arr_Z<cbDict['f_min']).nonzero()])
    intpSclMat[:,:,1][Idx_outLmt[0],Idx_outLmt[1]] = 0.
    chemMat = np.insert(intpSclMat,0,Arr_c,axis=2)
    chemMat = np.insert(chemMat,0,Arr_Z,axis=2)
    stackMat = np.reshape(chemMat,[np.shape(chemMat)[0]*np.shape(chemMat)[1],
                              np.shape(chemMat)[2]])
    np.savetxt(fln,stackMat,fmt='%.5E')
    print('\nDone.')
    
''' ===========================================================================

Subroutine functions

=========================================================================== '''

def calculateCp_eff(T,cp_m,lam_sumCpdT):
    cp_e = np.zeros(np.shape(cp_m))
    T_0 = 298.15
    if abs(T[0] - T_0) > 0.1: cp_e[0] = lam_sumCpdT / (T[0] - T_0)
    else: cp_e[0] = cp_m[0]
    for ii in range(1,len(T)):
            tmp_sum = 0.
            for j in range(1,ii+1):
                tmp_sum = (tmp_sum + 0.5*(cp_m[j]+cp_m[j-1])*(T[j]-T[j-1]))
            cp_e[ii] = (tmp_sum + lam_sumCpdT) / abs(T[ii] - T_0)
    return cp_e        
    
def interp2D(M_Zc,Z_pts,c_pts,meshgrid):
    f = interpolate.interp2d(c_pts,Z_pts,M_Zc)
    intpM_Zc = f(meshgrid['c_space'],meshgrid['Z_space'])
    return intpM_Zc     

def generateTable2(Z0,Yc_eq0):
    Z0 = lamArr[:,1]
    Yc_eq0 = Yc_eq
    from scipy.interpolate import UnivariateSpline
    sp = UnivariateSpline(Z0,Yc_eq0,s=0)
    Z1 = np.linspace(min(Z0),0.085,101);
    Yc_eq1 = sp(Z1)
#    plt.plot(Z0,Yc_eq0,label='original')
#    plt.plot(Z1,Yc_eq1,label='spline')
#    plt.legend()
#    plt.show()
    d2 = sp.derivative(n=2)
    d2Yc_eq1 = d2(Z1)
    sp = UnivariateSpline(Z1,d2Yc_eq1)
    d2Yc_eq2 = sp(Z1)
    from scipy.signal import savgol_filter
    d2Yc_eq3 = savgol_filter(d2Yc_eq2, 11, 3)
    plt.plot(Z1,d2Yc_eq1,label='original')
    plt.plot(Z1,d2Yc_eq3,label='spline')
    plt.legend()
    plt.show()
    '''
        INCOMPLETE
                    '''
    
#def scatterInterp(M_Zc):
#    Zcoor_1D = np.reshape(Z_grid,[np.shape(Z_grid)[0]*np.shape(Z_grid)[1],1])
#    ccoor_1D = np.reshape(c_grid,[np.shape(c_grid)[0]*np.shape(c_grid)[1],1])
#    M_Zc1D = np.reshape(M_Zc,[np.shape(M_Zc)[0]*np.shape(M_Zc)[1],1])    
#    M_Zc_intp = np.squeeze(griddata(np.hstack((ccoor_1D,Zcoor_1D)),
#                                    M_Zc1D,(grid_Z, grid_c),
#                                    method='linear'
#                                    )
#                            )
#    return M_Zc_intp
