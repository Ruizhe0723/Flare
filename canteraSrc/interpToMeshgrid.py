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

    work_dir = cbDict['work_dir']

    # read laminar flame global parameters & write Table_5
    solFln = (work_dir + 'canteraData/solution_'
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

    fln = work_dir + cbDict['output_fln']
    with open(fln,'w') as strfile:
      #
      strfile.write('%.5E' % streamBC[0,6] + '\t' +
                    '%.5E' % streamBC[1,6] + '\n')
      #
      strfile.write(str(cbDict['nchemfile']) + '\n')
      #
      np.savetxt(strfile,
                 np.transpose([lamArr[:,1],lamArr[:,3],lamArr[:,4],
                               lamArr[:,5],lamArr[:,6]]),
                 fmt='%.5E',delimiter='\t')
      # strfile.write(str(cbDict['int_pts_z']) + '\t' +
      #             str(cbDict['int_pts_gz']) + '\n' )
      #
    strfile.close()

    # read cantera solutions & calculate c,omega_c
    nScalCant = cbDict['nSpeMech'] + cbDict['nVarCant']
    mainData = np.zeros((cbDict['nchemfile'],int(max(lamArr[:,2])),nScalCant))
    cIn = np.zeros(np.shape(mainData[:,:,0]))
    omega_cIn = np.zeros(np.shape(mainData[:,:,0]))
    # Yc_eq = np.zeros((cbDict['nchemfile']))
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

        if(cbDict['scaled_PV']):
          omega_cIn[i,:] = mainData[i,:,1] / mainData[i,imax,0]
        else:
          omega_cIn[i,:] = mainData[i,:,1]

        # Yc_eq[i] = mainData[i,-1,0]
        # # calculate d2Yc_eq/dZ2 & write Table 2
        # generateTable2(lamArr[:,1],Yc_eq)

    # interpolate in c space & write out for each flamelet
    MatScl_c = np.zeros((cbDict['nchemfile'],cbDict['cUnifPts'],
                         nScalCant+1))
    for i in range(cbDict['nchemfile']):
        len_grid = int(lamArr[i,2])
        ctrim = cIn[i,:len_grid-1]
        # 0:c|1:omg_c
        MatScl_c[i,:,0] = np.linspace(0.,1.,cbDict['cUnifPts'])
        MatScl_c[i,:,1] = np.matlib.interp(MatScl_c[i,:,0],ctrim,
                                         omega_cIn[i,:len_grid-1])
        # 2:T|3:rho|4:cp|5:mw|6:hf_0|7:nu|8:h|9:qdot
        for k in range(2,nScalCant):
            MatScl_c[i,:,k] = np.matlib.interp(MatScl_c[i,:,0],ctrim,
                                             mainData[i,:len_grid-1,k])
        # cp-->cp_e
        MatScl_c[i,:,4] = calculateCp_eff(MatScl_c[i,:,2],MatScl_c[i,:,4],
                                          lamArr[i,7])
        # Yc_max
        MatScl_c[i,:,-1] = mainData[i,len_grid-1,0]

        # write inpterpolated 1D profiles
        fln = (solFln+'Unf_'+casename+'_'+str('{:03d}'.format(i))+'.dat')
        print('\nWriting --> ' + fln)
        np.savetxt(fln,MatScl_c[i,:,:],fmt='%.5e')

    oriSclMat = np.zeros([cbDict['nchemfile']+2,cbDict['cUnifPts'],
                          nScalCant+1])
    ind_list_Yis = []
    for i in range(len(streamBC[:,0])):

        if i == 0: j = len(oriSclMat[:,0,0]) - 1
        else: j = 0
        for k in range(len(streamBC[0,:])):
            # for thermo scalars
            if k < len(streamBC[0,:]) - 2*cbDict['nYis']:
                oriSclMat[j,:,k+2] = streamBC[i,k]
            else:
                nk = k
                break

        # for selected species
        for s in range(cbDict['nYis']):
            ispc=int(streamBC[i,nk+s*2])
            if i == 0: ind_list_Yis.append(ispc)
            iscal = cbDict['nVarCant']+ispc
            oriSclMat[j,:,iscal]=streamBC[i,nk+s*2+1]

    oriSclMat[1:cbDict['nchemfile']+1,:,:] = MatScl_c

    intpSclMat = np.zeros([len(meshgrid['z_space']),len(meshgrid['c_space']),
                           cbDict['nVarCant']+2])
    intpYiMat = np.zeros([len(meshgrid['z_space']),len(meshgrid['c_space']),
                        cbDict['nYis']])
    Z_pts = np.insert(lamArr[:,1],0,[0.],axis=0)
    Z_pts = np.insert(Z_pts,len(Z_pts),[1.],axis=0)
    c_pts = MatScl_c[0,:,0]

    np.array(ind_list_Yis)

    print('\nInterpolating...')
    intpSclMat[:,:,0] = interp2D(oriSclMat[:,:,3],Z_pts,c_pts,meshgrid) # rho
    for k in [1,4,5,6]: # omega_c,cp_e,mw,hf
        intpSclMat[:,:,k] = interp2D(oriSclMat[:,:,k],Z_pts,c_pts,meshgrid)
    intpSclMat[:,:,7] = interp2D(oriSclMat[:,:,2],Z_pts,c_pts,meshgrid) # T
    intpSclMat[:,:,8] = interp2D(oriSclMat[:,:,7],Z_pts,c_pts,meshgrid) # nu
    intpSclMat[:,:,9] = interp2D(oriSclMat[:,:,8],Z_pts,c_pts,meshgrid) # h
    intpSclMat[:,:,10] = interp2D(oriSclMat[:,:,9],Z_pts,c_pts,meshgrid) # qdot
    for y in range(cbDict['nYis']):
        iy = ind_list_Yis[y]
        intpYiMat[:,:,y] = interp2D(oriSclMat[:,:,iy+cbDict['nVarCant']],
                                            Z_pts,c_pts,meshgrid)
    intpSclMat[:,:,11] = interp2D(oriSclMat[:,:,-1],Z_pts,c_pts,meshgrid) # qdot

    print('\nInterpolation done. ')

    if np.sum(np.isnan(intpSclMat)) > 0:
        print('\nNumber of Nans detected: ', np.sum(np.isnan(intpSclMat)))
    else: print('\nNo Nans detected. Well done!')

    print('\nwriting chemTab file...')
    Arr_c,Arr_Z = np.meshgrid(meshgrid['c_space'],meshgrid['z_space'])
    Idx_outLmt = np.hstack([(Arr_Z>cbDict['f_max']).nonzero(),
                             (Arr_Z<cbDict['f_min']).nonzero()])
    ind_rates=[1,9]
    intpSclMat[:,:,ind_rates][Idx_outLmt[0],Idx_outLmt[1]] = 0.
    chemMat = np.append(intpSclMat,intpYiMat,axis=2)
    chemMat = np.insert(chemMat,0,Arr_c,axis=2)
    chemMat = np.insert(chemMat,0,Arr_Z,axis=2)
    stackMat = np.reshape(chemMat,[np.shape(chemMat)[0]*np.shape(chemMat)[1],
                              np.shape(chemMat)[2]])

    fln = work_dir + 'chemTab_' + str('{:02d}'.format(solIdx)) + '.dat'
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
    intpM_Zc = f(meshgrid['c_space'],meshgrid['z_space'])
    return intpM_Zc

# def generateTable2(Z0,Yc_eq0):
#     Z0 = lamArr[:,1]
#     Yc_eq0 = Yc_eq
#     from scipy.interpolate import UnivariateSpline
#     sp = UnivariateSpline(Z0,Yc_eq0,s=0)
#     Z1 = np.linspace(min(Z0),0.085,101);
#     Yc_eq1 = sp(Z1)
# #    plt.plot(Z0,Yc_eq0,label='original')
# #    plt.plot(Z1,Yc_eq1,label='spline')
# #    plt.legend()
# #    plt.show()
#     d2 = sp.derivative(n=2)
#     d2Yc_eq1 = d2(Z1)
#     sp = UnivariateSpline(Z1,d2Yc_eq1)
#     d2Yc_eq2 = sp(Z1)nScalars
#     from scipy.signal import savgol_filter
#     d2Yc_eq3 = savgol_filter(d2Yc_eq2, 11, 3)
#     plt.plot(Z1,d2Yc_eq1,label='original')
#     plt.plot(Z1,d2Yc_eq3,label='spline')
#     plt.legend()
#     plt.show()
#     '''
#         INCOMPLETE
#                     '''

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
