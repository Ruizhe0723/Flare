!+---------------------------------------------------------------------+
!| The module is used to define global variables and arraries.         |
!+---------------------------------------------------------------------+
module integrate
  !
  implicit none
  !
  integer :: &
    n_points_z,n_points_c,int_pts_z,int_pts_c,int_pts_gcz,int_pts_gz, &
    int_pts_gc,nYis,nScalars,nchemfile,nSpeMech,n_points_h,ii
    !
  real(8),allocatable :: &
    c_int(:),z_int(:),gz_int(:),gc_int(:),gcz_int(:)
  real(8) :: fmix_min,fmix_max
  ! parameter (smaller = 1.0d-08)
  !
  contains
  !
  ! Z.Chen-----------interpolation points for laminar flame data------------------
  subroutine read_integrate_inp(my_id)
    !
    ! arguments
    integer, intent(in) :: my_id
    !
    open(unit=25, file='./integrate.inp', status='old')
    !
    if(my_id==0) print*, '>> reading from ./integrate.inp ...'
    read(25,*) nSpeMech
    if(my_id==0) print*, '>> read nSpeMech ---->', nSpeMech
    !
    read(25,*) fmix_min,fmix_max
    if(my_id==0) print*, '>> read fmix_min ---->',fmix_min
    if(my_id==0) print*, '>> read fmix_max ---->',fmix_max
    !
    read(25,*) nchemfile
    if(my_id==0) print*, '>> read nchemfile ---->', nchemfile
    !
    read(25,*) n_points_z,n_points_c,n_points_h
    if(my_id==0) print*, '>> read n_points_z ---->',n_points_z
    if(my_id==0) print*, '>> read n_points_c ---->',n_points_c
    if(my_id==0) print*, '>> read n_points_h ---->',n_points_h
    !
    ! Z.Chen-----------turbulent table dimension------------------------------------
    !
    read(25,*) int_pts_z,int_pts_c,int_pts_gz,int_pts_gc,int_pts_gcz
    if(my_id==0) print*, '>> read int_pts_z ---->',int_pts_z
    if(my_id==0) print*, '>> read int_pts_c ---->',int_pts_c
    if(my_id==0) print*, '>> read int_pts_gz ---->',int_pts_gz
    if(my_id==0) print*, '>> read int_pts_gc ---->',int_pts_gc
    if(my_id==0) print*, '>> read int_pts_gcz ---->',int_pts_gcz
    !
    allocate(z_int(int_pts_z),c_int(int_pts_c), &
      gz_int(int_pts_gz),gc_int(int_pts_gc),gcz_int(int_pts_gcz))
      !
    read(25,*) (z_int(ii),ii=1,int_pts_z)
    read(25,*) (c_int(ii),ii=1,int_pts_c)
    read(25,*) (gz_int(ii),ii=1,int_pts_gz)
    read(25,*) (gc_int(ii),ii=1,int_pts_gc)
    read(25,*) (gcz_int(ii),ii=1,int_pts_gcz)
    !
    read(25,*) nScalars,nYis
    if(my_id==0) print*, '>> read nScalars ---->',nScalars
    if(my_id==0) print*, '>> read nYis ---->',nYis
    !
  end subroutine read_integrate_inp
  !
end module integrate
