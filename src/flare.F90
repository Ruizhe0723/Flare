! ----------------------------------------------------------------------
! Copyright
!
! Application
!   FlaRe Combustion Model - Turbulent flame manifold
!
! Description
!   This Fortran code describes the flamelet manifold integration
!     to be used in the FlaRe Combustion Model.
! ======================================================================
! Dec 2020 -- Rewritten by Z.X. Chen @ Camrbidge
! ----------------------------------------------------------------------
program main
  !
  use mpi
  use integrate
  use pdf
  use func, only: locate
  !
  implicit none
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                          Declaration
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer ::  ierr,nproc,my_id,load_distro(0:9999),load_idx(0:9999)
      integer ::  iUnit_low,iUnit_high,nwork,iproc
      integer ::  i,j,k,l,m,i1,iii,iiii,jjjj,z_loc,c_loc,solIdx
      real(8),allocatable :: c_space(:),z_space(:),yint(:),Yi_int(:), &
                             Src_vals(:,:,:),Yi_vals(:,:,:)
      real(8) :: gc,gz,gcz
      real(8),parameter :: theta = 1.0d0
      character(len=2) :: strUnit2,strDouble

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                          MPI setup
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +           Interpolation from laminar flame table to chemTab
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    call read_integrate_inp(my_id)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !
    allocate( &
      c_space(n_points_c),z_space(n_points_z),yint(nScalars), &
      Yi_int(nYis),Src_vals(n_points_z,n_points_c,nScalars), &
      Yi_vals(n_points_z,n_points_c,nYis),)

      DO solIdx = 1,n_points_h

        write(strDouble,'(I2.2)') solIdx
        if(my_id==1) write(*,*) 'Reading chemTab...'
        open(unit=20, file='chemTab_'//strDouble//'.dat', status='old')
!        read(20,*) dumStr
!        read(20,*) dumStr
        do i=1,n_points_z
            do j=1,n_points_c
            read(20,*) z_space(i),c_space(j) &
                     ,(Src_vals(i,j,iii),iii=1,nScalars) &
                     ,(Yi_vals(i,j,iii),iii=1,nYis)
            enddo
        enddo
        close(20)
        if(my_id==1) write(*,*) 'Done reading chemTab.'
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                          MPI allocation
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Load distribution across cores
   nwork = int_pts_z
   do iproc=1,nproc
     load_distro(iproc) = floor(dble(nwork)/(nproc - iproc +1))
     nwork = nwork - load_distro(iproc)
   enddo
   load_idx(1) = load_distro(1)
   do iproc = 1,nproc-1
     load_idx(iproc+1) = load_idx(iproc) + load_distro(iproc+1)
   enddo
!        write(*,*) load_idx(1:nproc)
   if (my_id.eq.0) then
     iUnit_low = 1
     iUnit_high = load_idx(1)
   elseif (my_id.eq.(nproc-1)) then
     iUnit_low = load_idx(my_id) +1
     iUnit_high = int_pts_z
   else
     iUnit_low = load_idx(my_id) + 1
     iUnit_high = load_idx(my_id+1)
   endif

    write(*,*) my_id,' integrating from ',iUnit_low, '-',iUnit_high

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + +
! + +                       Final integration
! + +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     looping in Z space
         do i1 = iUnit_low,iUnit_high
         write(strUnit2,'(I2.2)') i1
         write(strDouble,'(I2.2)') solIdx
         open(unit=7, file='unit'//strUnit2//'_h' &
              //strDouble//'.dat')
!          open(unit=77, file='Yiunit'//strUnit2//'_'//strDouble//'.dat')

!     looping in c space
         do j=1,int_pts_c
!     looping in gZ space
           do k=1,int_pts_gz
!     looping in gc space
             do l=1,int_pts_gc
!     looping in gZc space
               do m=1,int_pts_gcz       !!*!! gZc=0
!     convert to variances

           z_loc=locate(z_space,n_points_z,z_int(i1))
           c_loc=locate(c_space,n_points_c,c_int(j))
           gc=gc_int(l)*(c_int(j)*(1.0-c_int(j)))
           gz=gz_int(k)*(z_int(i1)*(1.0-z_int(i1)))

           !! Z=0 or Z=1
           IF((i1==1).or.(i1==int_pts_z))then
             yint(2) = 0.
             yint(3) = 0.
             yint(4) = 0.

             if(i1==int_pts_z)then
               z_loc = z_loc + 1
             endif
             do jjjj=5,nScalars
              yint(jjjj) = Src_vals(z_loc,1,jjjj)
             enddo
             do iiii=1,nYis
              Yi_int(iiii) = Yi_vals(z_loc,1,iiii)
             enddo

           !! gZ=0 and gc=0
         ELSEIF(  ((k==1).and.(l==1)) .or. ((k==1).and.(j==1)) &
              .or.((k==1).and.(j==int_pts_c))  )then
             call delta(z_int(i1),c_int(j),z_space,c_space,yint &
         ,Src_vals,Yi_int,Yi_vals)

           !! gZ=0 and gc>0
           ELSEIF(  (k==1).and.(l.gt.1)  )then
             call cbeta(c_int(j),gc,c_space,z_int(i1),z_space,yint &
         ,Src_vals,Yi_int,Yi_vals)

           !! gZ>0 and gc=0
           ELSEIF(  (k.gt.1).and.(l==1)  )then
             call zbeta(z_int(i1),gz,z_space,c_int(j),c_space,yint &
         ,Src_vals,Yi_int,Yi_vals)

           !! gZ>0 and gc>0
           ELSE
!             if((j==1))then
!             c_mean = small
!             elseif((j==int_pts_c))then
!             c_mean = 1.-small
!             else
!             c_mean = c_int(j)
!             endif
!             gc=gc_int(l)*(c_mean*(1.0-c_mean))
             if((j==1).or.(j==int_pts_c))then
               goto 99
             endif

             gcz = gcz_int(m)*sqrt(gc)*sqrt(gz)*0.98
               call int_point(z_int(i1),c_int(j),gc &
                  ,gz,gcz,yint,z_space,c_space,Src_vals &
                  ,Yi_int,Yi_vals,theta)

           ENDIF

!!*!! write out the scalars needed in cgs units
!     1:rho 2:omegac 3:c*omegac 4:Cp 5:MW 6:Hf0 7:T 8:nu 9:h 10:qdot

 99       write(7,200)z_int(i1),c_int(j),gz_int(k),gc_int(l),gcz_int(m)&
              ,(yint(iii),iii=2,3),(yint(iii),iii=5,nScalars) &
              ,(Yi_int(iii),iii=1,nYis)

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
       call MPI_FINALIZE ( ierr )

       if(my_id==1) write(*,*) "Done"

 200  format(25e15.5)

    end program main
