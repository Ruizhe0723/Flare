! ----------------------------------------------------------------------------------
!                                              |
!   =======     //        \\      F laRe       |  University of Cambridge
!   ||         // \\    // \\     T urbulent   |  Department of Engineering
!   ||        //   \\  //   \\    F lame       |  Hopkinson Laboratory
!   =======  //     \\//     \\   M anifold    |
!                                              |
! CAMBRIDGE-MHI Combustion Instability Project |  Copyright (C) 2016-2019
! ------------------------------------------------------------------------------
! Copyright
!     This file is distributed within the Camrbidge-MHI framework only.
!
! Application
!     Cambridge FlaRe Combustion Model - Turbulent flame manifold integration
!
! Description
!     This Fortran code describes the flamelet manifold integration
!     to be used in the Cambridge FlaRe Combustion Model.
! ------------------------------------------------------------------------------

      program riemannCZ
      implicit none

      include 'mpif.h'
!      include 'data.inc'
      include 'integrate.inc'
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                          Declaration
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer ierr, nproc, my_id
      integer iUnit_low, iUnit_high

      integer i,j,k,l,m,i1,iii,iiii,jjjj,z_loc,c_loc,nn,solIdx
      double precision c_space(n_points_c), z_space(n_points_z)
     &     ,c_int(int_pts_c),z_int(int_pts_z)
     &     ,gz_int(int_pts_gz),gc_int(int_pts_gc),gcz_int(int_gcz)
     &     ,gc,gz,gcz!,dum,c_mean
     &     ,yint(nScalars),Yi_int(nYis)
     &     ,Src_vals(n_points_z,n_points_c,nScalars)
     &     ,Yi_vals(n_points_z,n_points_c,nYis)
      double precision theta,expnts(int_pts_gz-1)
      parameter(theta = 1.0D0)

      integer load_distro(1000),load_idx(1000)
      integer nwork,iproc

      character*2 strUnit2,strDouble
      ! character*50 dumStr

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                          MPI setup
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
! Load distribution across cores
       nwork = int_pts_z - 2!ZC: fist and last are instant
       do iproc=1,nproc
         load_distro(iproc) = floor(dble(nwork)/(nproc - iproc +1))
         nwork = nwork - load_distro(iproc)
       enddo
       load_idx(:) = 0
       load_distro(:) = 0
       load_idx(1) = load_distro(1) + 1
       do iproc = 1,nproc-1
         load_idx(iproc+1) = load_idx(iproc)+load_distro(iproc+1)
       enddo

       if (my_id.eq.0) then
         iUnit_low = 1
         iUnit_high = load_idx(1)
       elseif (my_id.eq.(nproc-1)) then
         iUnit_low = load_idx(my_id)+1
         iUnit_high = int_pts_z
       else
         iUnit_low = load_idx(my_id) + 1
         iUnit_high = load_idx(my_id+1)
       endif
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +           Interpolation from laminar flame table to chemTab
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO solIdx = 1,n_points_h

!        call interpLamFlame(z_space,c_space,Src_vals,Yi_vals,solIdx)
!        goto 10

        write(strDouble,'(I2.2)') solIdx
        write(*,*) 'Reading chemTab...'
        open(unit=20, file='chemTab_'//strDouble//'.dat', status='old')
!        read(20,*) dumStr
!        read(20,*) dumStr
        do i=1,n_points_z
            do j=1,n_points_c
            read(20,*) z_space(i),c_space(j)
     &                ,(Src_vals(i,j,iii),iii=1,nScalars)
!     &                ,(Yi_vals(i,j,jjj),jjj=1,nYis)
            enddo
        enddo
        close(20)
       write(*,*) 'Done reading chemTab.'

! 	open(unit=32,file='canteraData/solution_00/lamParameters.txt'
!      &     ,status='old')
! 	read(32,*) dumStr
! 	do i = 1,nchemfile
! 	  read(32,*) dum,z_lam(i),(dum,j=1,6)
! 	enddo
! 	close(32)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                          MPI allocation
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(*,*) my_id, ' integrating from ', iUnit_low, '-',
     &  iUnit_high
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +       Assign discretisation number for 5 control parameters:
! + +                    Z_Til,c_Til,gZ,gc,gZc
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!*!! Z space: non-linear, refined near Zst
      z_int(1) = 0.0
      z_int(int_pts_z) =  1.0
      nn = int_pts_z/20*19
      do i = 2,nn
        z_int(i) = fmix_max*1.2/(nn-1)*(i-1)
      enddo
      do i = nn+1,int_pts_z-1
        z_int(i) = fmix_max*1.2+(1.-fmix_max*1.2)/(int_pts_z-nn)*(i-nn)
      enddo


!     c space: linear between 0 and 1

      c_int(1) = 0.0
      c_int(int_pts_c) = 1.0

      do i = 2,int_pts_c-1
         c_int(i) = c_int(1) +(i - 1)*(c_int(int_pts_c)-c_int(1))
     &                               /(int_pts_c -1)
      enddo

!     gZ space: non-linear
!!*!! the upper limit depends on the cold flow solution usually <0.2
!      gz_int(int_pts_gz) = 0.05
!      gz_int(1) = 0.
!      gz_int(int_pts_gz-1) = 0.03
!      gz_int(int_pts_gz-2) = 0.02
!      do i = 2, int_pts_gz-3
!        gz_int(i) = gz_int(1) +(i - 1)*(gz_int(int_pts_gz-2)-gz_int(1))
!     &        /(int_pts_gz -3)
!      enddo

      gz_int(1) = 0.
      do i = 1,int_pts_gz-1
        expnts(i) = -4. + (i-1)*(-1.-(-4.))/(int_pts_gz-2)
        gz_int(i+1) = 10.**expnts(i)
      enddo

!     gc space: linear between 0 and 1
      gc_int(1) = 0.0
      gc_int(int_pts_gc) = 1. - smaller

      do i=2,int_pts_gc-1
          gc_int(i) = gc_int(1) +(i - 1)*(gc_int(int_pts_gc)-gc_int(1))
     &                                  /(int_pts_gc -1)
      enddo

!      gZc space: linear between -1 and 1
!      do i = 1,int_gcz
!        gcz_int(i) = -1.0 + 2./(int_gcz-1)*(i-1)
!      enddo
      gcz_int(1) = 0.

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                       Final integration
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     looping in Z space
         do i1 = 3,3!iUnit_low,iUnit_high
         write(strUnit2,'(I2.2)') i1
         write(strDouble,'(I2.2)') solIdx
         open(unit=7, file='unit'//strUnit2//'_h'//strDouble//'.dat')
!          open(unit=77, file='Yiunit'//strUnit2//'_'//strDouble//'.dat')

!     looping in c space
         do j=1,int_pts_c
!     looping in gZ space
           do k=1,int_pts_gz
!     looping in gc space
             do l=1,int_pts_gc
!     looping in gZc space
               do m=1,int_gcz       !!*!! gZc=0
!     convert to variances

           call locate2(z_space,n_points_z,z_int(i1),z_loc)
           call locate2(c_space,n_points_c,c_int(j),c_loc)
           gc=gc_int(l)*(c_int(j)*(1.0-c_int(j)))
           gz=gz_int(k)*(z_int(i1)*(1.0-z_int(i1)))

           !! Z=0 or Z=1
           IF((i1.eq.1).or.(i1.eq.int_pts_z))then
             yint(2) = 0.
             yint(3) = 0.
             yint(4) = 0.

             if(i1.eq.int_pts_z)then
               z_loc = z_loc + 1
             endif
             do jjjj=5,nScalars
              yint(jjjj) = Src_vals(z_loc,1,jjjj)
             enddo
             do iiii=1,nYis
              Yi_int(iiii) = Yi_vals(z_loc,1,iiii)
             enddo

           !! gZ=0 and gc=0
           ELSEIF(  ((k.eq.1).and.(l.eq.1))
     &          .or.((k.eq.1).and.(j.eq.1))
     &          .or.((k.eq.1).and.(j.eq.int_pts_c))  )then
             call delta(z_int(i1),c_int(j),z_space,c_space,yint
     &    ,Src_vals,Yi_int,Yi_vals)

           !! gZ=0 and gc>0
           ELSEIF(  (k.eq.1).and.(l.gt.1)  )then
             call cbeta(c_int(j),gc,c_space,z_int(i1),z_space,yint
     &    ,Src_vals,Yi_int,Yi_vals)

           !! gZ>0 and gc=0
           ELSEIF(  (k.gt.1).and.(l.eq.1)  )then
             call zbeta(z_int(i1),gz,z_space,c_int(j),c_space,yint
     &    ,Src_vals,Yi_int,Yi_vals)

           !! gZ>0 and gc>0
           ELSE
!             if((j.eq.1))then
!             c_mean = small
!             elseif((j.eq.int_pts_c))then
!             c_mean = 1.-small
!             else
!             c_mean = c_int(j)
!             endif
!             gc=gc_int(l)*(c_mean*(1.0-c_mean))
             if((j.eq.1).or.(j.eq.int_pts_c))then
               goto 99
             endif

             gcz = gcz_int(m)*sqrt(gc)*sqrt(gz)*0.98
               call int_point(z_int(i1),c_int(j),gc
     &             ,gz,gcz,yint,z_space,c_space,Src_vals
     &             ,Yi_int,Yi_vals,theta)

           ENDIF

!!*!! write out the scalars needed in cgs units
!     1:rho 2:omegac 3:c*omegac 4:Z*omegac 5:Cp 6:MW 7:formEnthal 8:T

 99       write(7,200)z_int(i1),c_int(j),gz_int(k),gc_int(l),gcz_int(m)
     &           ,(yint(iii),iii=2,3),(yint(iii),iii=5,11)

!        write(77,200)z_int(i1),c_int(j),gz_int(k),gc_int(l),gcz_int(m)
!      &       ,(Yi_int(iii),iii=1,nYis)

               enddo
             enddo
           enddo
         enddo

       close(7)
!        close(77)
       enddo

      ENDDO
       call MPI_FINALIZE ( ierr ) ! nakd2

! 100  format(a20,e9.3,a5,e9.3,a6,e9.3,a6,e9.3,a7,16e13.6)
 200  format(25e15.5)

      write(*,*) "Done"

      end
