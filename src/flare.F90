! ----------------------------------------------------------------------
! Copyright
!
! Application
!   FlaRe Combustion Model-Turbulent flame manifold
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
  !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++
! ++                          Declaration
! ++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer ::  ierr,nproc,my_id,load_distro(0:9999),load_idx(0:9999)
  integer ::  iUnit_low,iUnit_high,nwork,iproc
  integer ::  i,j,iz,ic,igz,igc,igcz,jj,z_loc,c_loc,solidx
  real(8) :: c_var,z_var,co_var
  character(len=2) :: str_iz,str_ih

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++
! ++                          MPI setup
! ++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  !
  call read_integrate_inp(my_id)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !
  nwork=int_pts_z
  !
  do iproc=1,nproc
    !
    load_distro(iproc)=floor(dble(nwork)/(nproc-iproc+1))
    nwork=nwork-load_distro(iproc)
    !
  enddo
  !
  load_idx(1)=load_distro(1)
  !
  do iproc=1,nproc-1
   load_idx(iproc+1)=load_idx(iproc)+load_distro(iproc+1)
  enddo
  !
  if(my_id==0) then
    !
    iUnit_low=1
    iUnit_high=load_idx(1)
    !
  elseif(my_id.eq.(nproc-1)) then
    !
    iUnit_low=load_idx(my_id)+1
    iUnit_high=int_pts_z
    !
  else
    !
   iUnit_low=load_idx(my_id)+1
   iUnit_high=load_idx(my_id+1)
   !
  endif

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++
! ++           Interpolation from laminar flame table to chemTab
! ++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  do solidx=1,n_points_h
    !
    write(str_ih,'(I2.2)') solidx
    !
    if(my_id==1) print*, 'Reading chemTab...'
    !
    open(unit=20, file='chemTab_'//str_ih//'.dat', status='old')
    !
    do i=1,n_points_z
      !
      do j=1,n_points_c
        !
        read(20,*) z_space(i),c_space(j), &
                   (Src_vals(i,j,jj),jj=1,nScalars), &
                   (Yi_vals(i,j,jj),jj=1,nYis)
      enddo
      !
    enddo
    !
    close(20)
    !
    if(my_id==1) print*, 'Done reading chemTab.'
    !
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++
! ++                       Final integration
! ++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    print*,'Proc',my_id,' integrate',iUnit_low, '-',iUnit_high
    !
    do iz=iUnit_low,iUnit_high
      !
      write(str_iz,'(I2.2)') iz
      !
      write(str_ih,'(I2.2)') solidx
      !
      open(unit=7, file='unit'//str_iz//'_h'//str_ih//'.dat')
      !
      z_loc=locate(z_space,n_points_z,z_int(iz))
      !
      do ic=1,int_pts_c
        !
        c_loc=locate(c_space,n_points_c,c_int(ic))
        !
        do igz=1,int_pts_gz
          !
          do igc=1,int_pts_gc
            !
            do igcz=1,int_pts_gcz
              !
              ! z = 0 or 1
              if(iz==1.or.iz==int_pts_z) then
                !
                yint(2)=0.d0
                yint(3)=0.d0
                yint(4)=0.d0
                !
                if(iz==1) z_loc=1
                if(iz==int_pts_z) z_loc=n_points_z
                !
                yint(5:nScalars)=Src_vals(z_loc,1,5:nScalars)
                Yi_int(:)=Yi_vals(z_loc,1,:)
                !
              ! gZ=0 and gc=0
              elseif((igz==1.and.igc==1).or.(igz==1.and.ic==1) &
                .or.(igz==1.and.ic==int_pts_c)) then
                !
                call delta(z_int(iz),c_int(ic),z_space,c_space,yint, &
                  Src_vals,Yi_int,Yi_vals)
                  !
              ! gZ=0 and gc>0
              elseif(igz==1.and.igc>1) then
                !
                call cbeta(c_int(ic),gc_int(ii),c_space,z_int(iz), &
                  z_space,yint,Src_vals,Yi_int,Yi_vals)
                  !
              ! gZ>0 and gc=0
              elseif(igz>1.and.igc==1) then
                !
                call zbeta(z_int(iz),gz_int(ii),z_space,c_int(ic), &
                  c_space,yint,Src_vals,Yi_int,Yi_vals)
                  !
              ! gZ>0 and gc>0
              else
                !
                if((ic==1).or.(ic==int_pts_c)) goto 99
                !
                c_var=gc_int(igc)*(c_int(ic)*(1.0-c_int(ic)))
                z_var=gz_int(igz)*(z_int(iz)*(1.0-z_int(iz)))
                co_var=gcz_int(igcz)*sqrt(c_var)*sqrt(z_var)*0.98d0
                !
                call int_point(z_int(iz),c_int(ic),c_var,z_var,co_var, &
                  yint,z_space,c_space,Src_vals,Yi_int,Yi_vals,theta)
                  !
              endif
              ! 1:rho|2:omg_c|3:coc|4:zoc|5:cp|6:mw|7:Hf0
              ! 8:T|9:nu|10:Ycmax
99            write(7,200) &
            z_int(iz),c_int(ic),gz_int(igz),gc_int(igc),gcz_int(igcz), &
            (yint(jj),jj=2,3),(yint(jj),jj=5,nScalars), &
            (Yi_int(jj),jj=1,nYis)
              !
            enddo !igcz
          enddo !igc
        enddo !igz
      enddo !ic
      !
      close(7)
      !
    enddo !iz

  enddo !solidx
  !
  call MPI_FINALIZE(ierr)
  !
  if(my_id==1) print*, "Done"
  !
200  format(25e15.5)
  !
end program main
