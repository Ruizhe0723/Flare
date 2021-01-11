!+---------------------------------------------------------------------+
!| The module is used to define global variables and arraries.         |
!+---------------------------------------------------------------------+
module integrate
  !
  implicit none
  !
  integer :: &
    n_points_z,n_points_c,int_pts_z,int_pts_c,int_gcz,int_pts_gz, &
    int_pts_gc,nYis,nScalars,nchemfile,nSpeMech,n_points_h
    !
  real(8) :: small,clip,smaller,fmix_min,fmix_max
  parameter (small = 1.0d-04)
  parameter (smaller = 1.0d-08)
  parameter (clip = 2.0d-04)
  !
  contains
  !
  ! Z.Chen-----------interpolation points for laminar flame data------------------
  parameter(nSpeMech = 53)

  parameter(fmix_min = 2.85e-02,fmix_max = 9.50e-02)

  parameter(nchemfile = 20)

  parameter (n_points_z = 501)

  parameter (n_points_c = 401)

  parameter (n_points_h = 1)

  ! Z.Chen-----------turbulent table dimension------------------------------------

      parameter (int_pts_z = 20)  !best if mod(int_pts_z-2)/Nproc = 0
      parameter (int_pts_gz = 15)
      parameter (int_gcz = 1)
      parameter (int_pts_c = 21)
      parameter (int_pts_gc = 11)

      parameter (nScalars = 11)   !number of scalars to be integrated
      ! parameter (scalars = nScalars + nSpeMech)   !number of scalars in chemTab
      parameter (nYis = 3)      !number of Yis to be integrated

end module integrate
