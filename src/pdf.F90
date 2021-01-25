!+---------------------------------------------------------------------+
!| The module is used to compute the PDFs                              |
!+---------------------------------------------------------------------+
module pdf
  !
  implicit none
  !
  contains
  !
  ! ********************************************************************
  ! ********************************************************************
  subroutine delta(z_mean,c_mean,z_space,c_space,y_int,sc_vals, &
                   Yi_int,Yi_vals)
  ! ********************************************************************
  ! ********************************************************************
    use integrate, only: n_points_z,n_points_c,nScalars,nYis
    use func, only: locate,intfac
    !
    ! arguments
    real(8), intent(in) :: &
      z_mean,c_mean,z_space(:),c_space(:),sc_vals(:,:,:),Yi_vals(:,:,:)
    real(8), intent(out) :: y_int(:),Yi_int(:)
    ! local data
    integer :: i,z_loc,c_loc
    real(8) :: z_fac,c_fac
    !
    z_loc=locate(z_space,n_points_z,z_mean)
    z_fac=intfac(z_mean,z_space,z_loc)
    !
    c_loc=locate(c_space,n_points_c,c_mean)
    c_fac=intfac(c_mean,c_space,c_loc)
    !
    y_int(:)=0.d0
    Yi_int(:)=0.d0
    !
    do i=1,2     !    scalars !number of integral scalars
      y_int(i)= &
        (1.d0-c_fac)*(z_fac*sc_vals(z_loc+1,c_loc,i)+ &
                      (1.d0-z_fac)*(sc_vals(z_loc,c_loc,i)))+ &
        c_fac*(z_fac*sc_vals(z_loc+1,c_loc+1,i)+ &
               (1.d0-z_fac)*sc_vals(z_loc,c_loc+1,i))
    enddo
    !
    do i=5,nScalars     !    scalars !number of integral scalars
      y_int(i)= &
        (1.d0-c_fac)*(z_fac*sc_vals(z_loc+1,c_loc,i)+ &
                      (1.d0-z_fac)*(sc_vals(z_loc,c_loc,i)))+ &
        c_fac*(z_fac*sc_vals(z_loc+1,c_loc+1,i)+ &
               (1.d0-z_fac)*sc_vals(z_loc,c_loc+1,i))
    enddo
    !
    y_int(2)=y_int(2)/y_int(1) !Zhi:omega_c/rho
    y_int(3)=y_int(2)*c_mean*y_int(nScalars) ! gc_source
    y_int(4)=y_int(2)*z_mean ! gz_source
    !
    do i=1,nYis     !    scalars !number of integral scalars Yi
      Yi_int(i)= (1.-c_fac)*(z_fac*Yi_vals(z_loc+1,c_loc,i)+ &
          (1.-z_fac)*(Yi_vals(z_loc,c_loc,i)))+ &
          c_fac*(z_fac*Yi_vals(z_loc+1,c_loc+1,i)+ &
          (1.-z_fac)*Yi_vals(z_loc,c_loc+1,i))
    enddo
    !
    do i=1,nScalars
      if(abs(y_int(i)).lt.1.e-20) y_int(i)=0.d0
    enddo
    do i=1,nYis
      if(Yi_int(i).lt.1.e-20) Yi_int(i)=0.d0
    enddo
    !
  end subroutine delta
  !
  ! ********************************************************************
  ! ********************************************************************
  subroutine cbeta(mean,g_var,space,z_mean,z_space,y_int,sc_vals, &
                   Yi_int,Yi_vals)
  ! ********************************************************************
  ! ********************************************************************
    !
    use integrate, only: n_points_z,n_points_c,nScalars,nYis
    use func, only: betai,locate,intfac
    !
    ! arguments
    real(8), intent(in) :: &
      mean,g_var,space(:),z_mean,z_space(:),sc_vals(:,:,:), &
      Yi_vals(:,:,:)
    real(8), intent(out) :: y_int(:),Yi_int(:)
    ! local data
    integer :: i,loc,j
    real(8) :: fac,alpha_c,beta_c,Yc_0,Yc_1
    real(8) :: dYdc(nYis,n_points_c+1),cdf001(n_points_c), &
      dsdc(nScalars,n_points_c+1),sc_vals_int(n_points_c,nScalars), &
      Yi_vals_int(n_points_c,nScalars)
      !
    loc=locate(z_space,n_points_z,z_mean)
    fac=intfac(z_mean,z_space,loc)
    !
    do i=1,n_points_c
      !
      do j=1,nScalars
        sc_vals_int(i,j)=fac*sc_vals(loc+1,i,j) &
            +((1.-fac)*sc_vals(loc,i,j))
      enddo
      !
      do j=1,nYis
        Yi_vals_int(i,j)=fac*Yi_vals(loc+1,i,j) &
            +((1.-fac)*Yi_vals(loc,i,j))
      enddo
      !
    enddo
    !
    y_int(:)=0.d0
    Yi_int(:)=0.d0
    !
    alpha_c=mean*((1.d0/g_var)-1.d0)
    beta_c=(1.d0-mean)*((1.d0/g_var)-1.d0)
    !
    !     ===== integration by part ====
    do i=1,n_points_c-1
      cdf001(i) =betai(alpha_c,beta_c,(space(i)+space(i+1))/2.0)
    enddo
    cdf001(n_points_c)=1.d0
    !
    do j=1,nScalars
      !
      ! compute derivatives
      do i=1,n_points_c-1
        !
        if(j==1 .or. j>=5) then !1:rho 5:cpe 6:mwt 7:sum_dhfi
          dsdc(j,i)=(sc_vals_int(i+1,j)-sc_vals_int(i,j)) &
              /(space(i+1)-space(i))
        elseif(j==2) then !2:omegac
          dsdc(j,i)=(sc_vals_int(i+1,2)/sc_vals_int(i+1,1) &
                -sc_vals_int(i,2)/sc_vals_int(i,1)) &
              /(space(i+1)-space(i))
        elseif(j==3) then !3:c*omegac
          Yc_0=space(i)*sc_vals_int(i,nScalars)
          Yc_1=space(i+1)*sc_vals_int(i+1,nScalars)
          dsdc(j,i)=(Yc_1*sc_vals_int(i+1,2)/sc_vals_int(i+1,1) &
                -Yc_0*sc_vals_int(i,2)/sc_vals_int(i,1)) &
              /(space(i+1)-space(i))
        else !4:Z*omegac
          dsdc(j,i)=(z_mean*sc_vals_int(i+1,2)/sc_vals_int(i+1,1) &
                -z_mean*sc_vals_int(i,2)/sc_vals_int(i,1)) &
              /(space(i+1)-space(i))
        endif
        !
      enddo !n_points_c-1
      !
      dsdc(j,n_points_c)=dsdc(j,n_points_c-1)
      dsdc(j,n_points_c+1)=dsdc(j,1)
      !
      ! start integration
      do i=1,n_points_c-2 !n_points_c-1 points, so n_points_c-2 sections
        y_int(j)= &
          y_int(j)-0.5*(dsdc(j,i)*cdf001(i)+dsdc(j,i+1)*cdf001(i+1))  &
            *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
      enddo
      !
      y_int(j)=y_int(j) &
        - dsdc(j,n_points_c)*cdf001(n_points_c) &!last half rectangle
          *(space(n_points_c)-space(n_points_c-1))/2.0 &!first half
        - dsdc(j,n_points_c+1)*cdf001(1) &!first half rectangle
          *(space(2)-space(1))/2.0+sc_vals_int(n_points_c,j)
          !
    enddo !nScalars
    !
    do j=1,nYis
      ! compute derivatives
      do i=1,n_points_c-1
         dYdc(j,i)=(Yi_vals_int(i+1,j)-Yi_vals_int(i,j)) &
                /(space(i+1)-space(i))
      enddo
      !
      dYdc(j,n_points_c)=dYdc(j,n_points_c-1)
      dYdc(j,n_points_c+1)=dYdc(j,1)
      !
      ! start integration
      do i=1,n_points_c-2 !n_points_c-1 points, so n_points_c-2 sections
        Yi_int(j)=Yi_int(j) &
             - 0.5*(dYdc(j,i)*cdf001(i)+dYdc(j,i+1)*cdf001(i+1)) &
           *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
      enddo
      !
      Yi_int(j)=Yi_int(j) &
        - dYdc(j,n_points_c)*cdf001(n_points_c) &!last half rectangle
          *(space(n_points_c)-space(n_points_c-1))/2.0 &!first half
        - dYdc(j,n_points_c+1)*cdf001(1) &!first half rectangle
          *(space(2)-space(1))/2.0+Yi_vals_int(n_points_c,j)
          !
    enddo !nYis
    !
    do i=1,nScalars
      if(abs(y_int(i)).lt.1.e-20) y_int(i)=0.d0
    enddo
    !
    do i=1,nYis
      if(Yi_int(i).lt.1.e-20) Yi_int(i)=0.d0
    enddo
    !
  end subroutine cbeta
  !
  ! ********************************************************************
  ! ********************************************************************
  subroutine zbeta(mean,g_var,space,c_mean,c_space,y_int,sc_vals, &
    Yi_int,Yi_vals)
  ! ********************************************************************
  ! ********************************************************************
    !
    use integrate, only: n_points_z,n_points_c,nScalars,nYis!,smaller
    use func, only: betai,locate,intfac
    !
    ! arguments
    real(8), intent(in) :: mean,g_var,space(:),c_mean,c_space(:), &
      sc_vals(:,:,:), Yi_vals(:,:,:)
    real(8), intent(out) :: y_int(:),Yi_int(:)
    ! local data
    integer :: i,loc,j
    real(8) :: fac,alpha_z,beta_z,Yc_0,Yc_1
    real(8) :: dYdc(nYis,n_points_z+1),cdf001(n_points_z), &
      dsdc(nScalars,n_points_z+1),sc_vals_int(n_points_z,nScalars), &
      Yi_vals_int(n_points_z,nScalars)
      !
    loc=locate(c_space,n_points_c,c_mean)
    fac=intfac(c_mean,c_space,loc)
    !
    !     interp sc_vals
    do i=1,n_points_z
      !
      do j=1,nScalars
        sc_vals_int(i,j)=fac*sc_vals(i,loc+1,j) &
            +((1.-fac)*sc_vals(i,loc,j))
      enddo
      !
      do j=1,nYis
        Yi_vals_int(i,j)=fac*Yi_vals(i,loc+1,j) &
            +((1.-fac)*Yi_vals(i,loc,j))
      enddo
      !
    enddo
    !
    y_int(:)=0.d0
    Yi_int(:)=0.d0
    !
    alpha_z=mean*((1.d0/g_var)-1.d0)
    beta_z=(1.d0-mean)*((1.d0/g_var)-1.d0)
    !
    !      ===== integration by part ====
    do i=1,n_points_z-1
        cdf001(i) =betai(alpha_z,beta_z,(space(i)+space(i+1))/2.0)
    enddo
    cdf001(n_points_z)=1.d0
    !
    do j=1,nScalars
      !! get derivative
      do i=1,n_points_z-1
        !
        if((j==1) .or. (j>=5))then !1:density 5:cpe 6:mwt 7:sum_dhfi
          dsdc(j,i)=(sc_vals_int(i+1,j)-sc_vals_int(i,j)) &
               /(space(i+1)-space(i))
               !
        elseif(j==2) then !2:omegac
          dsdc(j,i)= (sc_vals_int(i+1,2)/sc_vals_int(i+1,1) &
                 -sc_vals_int(i,2)/sc_vals_int(i,1)) &
               /(space(i+1)-space(i))
               !
        elseif(j==3) then !3:c*omegac
          Yc_0=c_mean*sc_vals_int(i,nScalars)
          Yc_1=c_mean*sc_vals_int(i+1,nScalars)
          dsdc(j,i)=(Yc_1*sc_vals_int(i+1,2)/sc_vals_int(i+1,1) &
                 -Yc_0*sc_vals_int(i,2)/sc_vals_int(i,1)) &
               /(space(i+1)-space(i))
               !
        else !4:Z*omegac
          dsdc(j,i)=(space(i+1)*sc_vals_int(i+1,2)/sc_vals_int(i+1,1) &
                 -space(i)*sc_vals_int(i,2)/sc_vals_int(i,1)) &
               /(space(i+1)-space(i))
        endif
        !
      enddo !n_points_z-1
      !
      dsdc(j,n_points_z)=dsdc(j,n_points_z-1)
      dsdc(j,n_points_z+1)=dsdc(j,1)
      !
      ! start integration
      do i=1,n_points_z-2 !n_points_z-1 points, so n_points_z-2 sections
         y_int(j)=y_int(j) &
               - 0.5*(dsdc(j,i)*cdf001(i)+dsdc(j,i+1)*cdf001(i+1)) &
          *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
      enddo
      y_int(j)=y_int(j) &
       - dsdc(j,n_points_z)*cdf001(n_points_z) &!last half rectangle
         *(space(n_points_z)-space(n_points_z-1))/2.0 &!first half
       - dsdc(j,n_points_z+1)*cdf001(1) &!first half rectangle
         *(space(2)-space(1))/2.0+sc_vals_int(n_points_z,j)
    enddo !nScalars
    !

    do j=1,nYis
      ! compute derivatives
      do i=1,n_points_z-1
         dYdc(j,i)= (Yi_vals_int(i+1,j)-Yi_vals_int(i,j)) &
                   /(space(i+1)-space(i))
      enddo
      !
      dYdc(j,n_points_z)=dYdc(j,n_points_z-1)
      dYdc(j,n_points_z+1)=dYdc(j,1)
      ! start integration
      do i=1,n_points_z-2 !n_points_z-1 points, so n_points_z-2 sections
         Yi_int(j)=Yi_int(j) &
               - 0.5*(dYdc(j,i)*cdf001(i)+dYdc(j,i+1)*cdf001(i+1)) &
               *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
      enddo
         Yi_int(j)=Yi_int(j) &
           - dYdc(j,n_points_z)*cdf001(n_points_z) &!last half rectangle
             *(space(n_points_z)-space(n_points_z-1))/2.0 &!first half
           - dYdc(j,n_points_z+1)*cdf001(1) &!first half rectangle
             *(space(2)-space(1))/2.0 &
                   +Yi_vals_int(n_points_z,j)
    enddo
    !
    do i=1,nScalars
      if(abs(y_int(i)).lt.1.e-20) y_int(i)=0.d0
    enddo
    !
    do i=1,nYis
      if(Yi_int(i).lt.1.e-20) Yi_int(i)=0.d0
    enddo
    !
  end subroutine zbeta

! **********************************************************************
! **********************************************************************
  subroutine int_point(z_mean,c_mean,c_var,z_var,co_var,y_int &
    ,z_space,c_space,sc_vals,Yi_int,Yi_vals,theta)
! **********************************************************************
! **********************************************************************
    !
    use integrate, only: n_points_z,n_points_c,nScalars,nYis
    !
    ! arguments
    real(8), intent(in) :: &
      z_mean,c_mean,c_var,z_var,co_var,z_space(:),c_space(:), &
      sc_vals(:,:,:),Yi_vals(:,:,:),theta
    real(8), intent(out) :: y_int(:),Yi_int(:)
    ! local data
    integer :: i,j,k
    real(8) :: alpha_z,beta_z,alpha_c,beta_c,rho
    real(8) :: &
      CDF_C(n_points_c+1),CDF_Z(n_points_z+1), &
      plaF(n_points_z+1,n_points_c+1), &
      dPsidc(n_points_z+1,n_points_c+1,nScalars), &
      Psi(n_points_z+1,n_points_c+1,nScalars), &
      Q_int(n_points_z+1,nScalars),dQdz(n_points_z+1,nScalars), &
      YiPsi(n_points_z+1,n_points_c+1,nYis),YiQ_int(n_points_z+1,nYis),&
      dYiQdz(n_points_z+1,nYis),dYiPsidc(n_points_z+1,n_points_c+1,nYis)
      !
    rho=co_var/(dsqrt(z_var)*dsqrt(c_var))
    ! if(abs(rho).ge.1.0e03) then
    !    theta=1.
    ! endif
    !
    alpha_z=z_mean*(((z_mean*(1.d0-z_mean))/z_var)-1.d0)
    alpha_c=c_mean*(((c_mean*(1.d0-c_mean))/c_var)-1.d0)
    beta_z=(1.d0-z_mean)*(((z_mean*(1.d0-z_mean))/z_var)-1.d0)
    beta_c=(1.d0-c_mean)*(((c_mean*(1.d0-c_mean))/c_var)-1.d0)
    !
    y_int(:)=0.d0
    Q_int(:,:)=0.d0
    Yi_int(:)=0.0
    YiQ_int(:,:)=0.0
    !
    ! -----assign values for Psi and obtain CDFs-----------------------
    ! main Z array
    do i=1,n_points_z-1
      !!!for Psi(1:n_points_z-1,1:n_points_c-1)
      do j=1,n_points_c-1
        call plackettFunc((c_space(j)+c_space(j+1))/2.0, &
                       (z_space(i)+z_space(i+1))/2.0, &
                       alpha_z,beta_z,alpha_c,beta_c, &
                       theta,CDF_C(j),CDF_Z(i), &
                       plaF(i,j))
        do k=5,nScalars
          Psi(i,j,k)=sc_vals(i,j,k)*plaF(i,j)
        enddo
        Psi(i,j,2)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
        Psi(i,j,3)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*c_space(j) &
                   *sc_vals(i,j,nScalars)
        Psi(i,j,4)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*z_space(i)
        do k=1,nYis
         YiPsi(i,j,k)=Yi_vals(i,j,k)*plaF(i,j)
        enddo
      enddo

      !for Psi(1:n_points_z-1,n_points_c)
      j=n_points_c
        call plackettFunc((c_space(j-1)+3.0*c_space(j))/4.0, &
                         (z_space(i)+z_space(i+1))/2.0, &
                         alpha_z,beta_z,alpha_c,beta_c, &
                         theta,CDF_C(j),CDF_Z(i), &
                         plaF(i,j))
      do k=5,nScalars
        Psi(i,j,k)=sc_vals(i,j,k)*plaF(i,j)
      enddo
      Psi(i,j,2)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
      Psi(i,j,3)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*c_space(j) &
                 *sc_vals(i,j,nScalars)
      Psi(i,j,4)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*z_space(i)
      do k=1,nYis
      YiPsi(i,j,k)=Yi_vals(i,j,k)*plaF(i,j)
      enddo

      !!for Psi(1:n_points_z-1,n_points_c+1) first point in c space
      j=n_points_c+1
      call plackettFunc((3.0*c_space(1)+c_space(2))/4.0, &
                       (z_space(i)+z_space(i+1))/2.0, &
                       alpha_z,beta_z,alpha_c,beta_c, &
                       theta,CDF_C(j),CDF_Z(i), &
                       plaF(i,j))
    enddo !n_points_z-1
    !
    !++! last point Z array
    i=n_points_z
      !!!for Psi(n_points_z,1:n_points_c-1)
      do j=1,n_points_c-1
        call plackettFunc((c_space(j)+c_space(j+1))/2.0,(z_space(i-1) &
                          +3.0*z_space(i))/4.0, &
                         alpha_z,beta_z,alpha_c,beta_c, &
                         theta,CDF_C(j),CDF_Z(i), &
                         plaF(i,j))
        do k=5,nScalars
          Psi(i,j,k)=sc_vals(i,j,k)*plaF(i,j)
        enddo
        Psi(i,j,2)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
        Psi(i,j,3)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*c_space(j) &
                   *sc_vals(i,j,nScalars)
        Psi(i,j,4)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*z_space(i)
        do k=1,nYis
          YiPsi(i,j,k)=Yi_vals(i,j,k)*plaF(i,j)
        enddo
      enddo
      !
      !!!!for Psi(n_points_z,1:n_points_c)
      j=n_points_c
      call plackettFunc((c_space(j-1)+3.0*c_space(j))/4.0, &
                       (z_space(i-1) &
                        +3.0*z_space(i))/4.0, &
                       alpha_z,beta_z,alpha_c,beta_c, &
                       theta,CDF_C(j),CDF_Z(i), &
                       plaF(i,j))
      do k=5,nScalars
        Psi(i,j,k)=sc_vals(i,j,k)*plaF(i,j)
      enddo
      Psi(i,j,2)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
      Psi(i,j,3)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*c_space(j) &
                 *sc_vals(i,j,nScalars)
      Psi(i,j,4)=sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)*z_space(i)
      do k=1,nYis
        YiPsi(i,j,k)=Yi_vals(i,j,k)*plaF(i,j)
      enddo
      !
      j=n_points_c+1
      !!!!for Psi(n_points_z,n_points_c+1) first point in c space
      call plackettFunc((3.0*c_space(1)+c_space(2))/4.0, &
                       (z_space(i-1)+z_space(i))/4.0, &
                       alpha_z,beta_z,alpha_c,beta_c, &
                       theta,CDF_C(j),CDF_Z(i), &
                       plaF(i,j))
    !++! first point Z array
    i=n_points_z+1
      !!!for Psi(n_points_z+1,1:n_points_c-1)
      do j=1,n_points_c-1
        call plackettFunc((c_space(j)+c_space(j+1))/2.0, &
                          (3.0*z_space(1)+z_space(2))/4.0, &
                          alpha_z,beta_z,alpha_c,beta_c, &
                          theta,CDF_C(j),CDF_Z(i), &
                          plaF(i,j))
      enddo
      !!!for Psi(n_points_z+1,n_points_c)
      j=n_points_c
      call plackettFunc((c_space(j-1)+3.0*c_space(j))/4.0, &
                        (3.0*z_space(1)+z_space(2))/4.0, &
                        alpha_z,beta_z,alpha_c,beta_c, &
                        theta,CDF_C(j),CDF_Z(i), &
                        plaF(i,j))
      !!!!for Psi(n_points_z+1,n_points_c+1) first point in c space
      j=n_points_c+1
      call plackettFunc((3.0*c_space(1)+c_space(2))/4.0, &
                        (3.0*z_space(1)+z_space(2))/4.0, &
                        alpha_z,beta_z,alpha_c,beta_c, &
                        theta,CDF_C(j),CDF_Z(i), &
                        plaF(i,j))

! -----2:omegac 3:c*omegac 4:z*omegac---------------------------------
    do k=2,3
      !calculate dPsidc(1:n_points_z,1:n_points_c-1)
      do i=1,n_points_z       ! looping on z

      do j=1,n_points_c-1    ! looping on c
       dPsidc(i,j,k)=(Psi(i,j+1,k)-Psi(i,j,k)+1.e-15) &
                     /(c_space(j+1) - c_space(j))
      enddo
      !last point in c space, dPsidc(1:n_points_z,n_points_c)
      j=n_points_c
        dPsidc(i,j,k)=dPsidc(i,j-1,k) !CZ:half rectangle
        !
      !first point in c space, dPsidc(1:n_points_z,n_points_c+1)
      j=n_points_c+1
        dPsidc(i,j,k)=dPsidc(i,1,k) !CZ:half rectangle
        !
      !start integration
      do j=1,n_points_c-2
       Q_int(i,k)=Q_int(i,k) &
         -0.5*(dPsidc(i,j,k)*CDF_C(j)+dPsidc(i,j+1,k)*CDF_C(j+1)) &
             *((c_space(j+2)+c_space(j+1)) &
                *0.5-(c_space(j+1)+c_space(j))*0.5)
      enddo
      !
      Q_int(i,k)=Q_int(i,k) &
        - dPsidc(i,n_points_c,k)*CDF_C(n_points_c) & !last half rect
          *(c_space(n_points_c)-c_space(n_points_c-1))/2.0 & !first half
        - dPsidc(i,n_points_c+1,k)*CDF_C(n_points_c+1) & !first half
          *(c_space(2)-c_space(1))/2.0 &
        +Psi(i,n_points_c,k)
      enddo

      !calculate dQdz(1:n_points_z-1)
      do i=1,n_points_z-1
        dQdz(i,k)=(Q_int(i+1,k) &
                 -Q_int(i,k)+1.e-15)/(z_space(i+1)-z_space(i))
      enddo
      !last point in z space, dQdz(n_points_z)
      i=n_points_z
        dQdz(i,k)=dQdz(i-1,k) !CZ:half rectangle
      !first point in z space, dQdz(n_points_z+1)
      i=n_points_z+1
        dQdz(i,k)=dQdz(1,k) !CZ:half rectangle
      !
      !start integration
      do i=1,n_points_z-2       ! looping on z
      y_int(k)=y_int(k) &
        -0.5*(dQdz(i,k)*CDF_Z(i)+dQdz(i+1,k)*CDF_Z(i+1)) &
        *((z_space(i+2)+z_space(i+1))*0.5-(z_space(i+1)+z_space(i))*0.5)
      enddo
      !
      y_int(k)=y_int(k) &
        +dQdz(n_points_z,k)*CDF_Z(n_points_z) &!last half
          *(z_space(n_points_z)-z_space(n_points_z-1))/2.0 & !first half
        +dQdz(n_points_z+1,k)*CDF_Z(n_points_z+1) & !last half rectangle
          *(z_space(2)-z_space(1))/2.0 & !first half rectangle ignored
        +Q_int(n_points_z,k)
    enddo !k=2,3

! c-----5:Cp 6:MW 7:formEnthal 8:T-----------------------------------
    do k=5,nScalars
      !
      !calculate dPsidc(1:n_points_z,1:n_points_c-1)
      do i=1,n_points_z       ! looping on z
        !
        do j=1,n_points_c-1    ! looping on c
          dPsidc(i,j,k)=(Psi(i,j+1,k)-Psi(i,j,k)) &
            /(c_space(j+1)-c_space(j))
        enddo
        !
        !last point in c space, dPsidc(1:n_points_z,n_points_c)
        j=n_points_c
          dPsidc(i,j,k)=dPsidc(i,j-1,k)
        !first point in c space, dPsidc(1:n_points_z,n_points_c+1)
        j=n_points_c+1
          dPsidc(i,j,k)=dPsidc(i,1,k)
          !
        !start integration
        do j=1,n_points_c-2
          Q_int(i,k)=Q_int(i,k) &
            -0.5*(dPsidc(i,j,k)*CDF_C(j)+dPsidc(i,j+1,k)*CDF_C(j+1)) &
               *((c_space(j+2)+c_space(j+1)) &
                  *0.5-(c_space(j+1)+c_space(j))*0.5)
        enddo
        !
        Q_int(i,k)=Q_int(i,k) &
          - dPsidc(i,n_points_c,k)*CDF_C(n_points_c) & !last half rect
            *(c_space(n_points_c)-c_space(n_points_c-1))/2.0 &
          - dPsidc(i,n_points_c+1,k)*CDF_C(n_points_c+1) & !first half
            *(c_space(2)-c_space(1))/2.0 &
          +Psi(i,n_points_c,k)
      enddo !i=1,n_points_z
      !
      !ccalculate dQdz(1:n_points_z-1)
      do i=1,n_points_z-1
        dQdz(i,k)=(Q_int(i+1,k)-Q_int(i,k))/(z_space(i+1)-z_space(i))
      enddo
      !last point in z space, dQdz(n_points_z)
      i=n_points_z
        dQdz(i,k)=dQdz(i-1,k) !CZ:half rectangle
      !first point in z space, dQdz(n_points_z+1)
      i=n_points_z+1
        dQdz(i,k)=dQdz(1,k) !CZ:half rectangle
        !
      !start integration
      do i=1,n_points_z-2       ! looping on z
        y_int(k)=y_int(k) &
          -0.5*(dQdz(i,k)*CDF_Z(i)+dQdz(i+1,k)*CDF_Z(i+1)) &
          *((z_space(i+2)+z_space(i+1)) &
          *0.5-(z_space(i+1)+z_space(i))*0.5)
      enddo
      !
      y_int(k)=y_int(k) &
        +dQdz(n_points_z,k)*CDF_Z(n_points_z) & !last half rect
           *(z_space(n_points_z)-z_space(n_points_z-1))/2.0 &
        +dQdz(n_points_z+1,k)*CDF_Z(n_points_z+1) & !first half rect
           *(z_space(2)-z_space(1))/2.0 &
        +Q_int(n_points_z,k)
    enddo !k=5,nScalars

! c-----Yis-----------------------------------
    do k=1,nYis
      !calculate dYiPsidc(1:n_points_z,1:n_points_c-1)
      do i=1,n_points_z       ! looping on z
        !
        do j=1,n_points_c-1    ! looping on c
          dYiPsidc(i,j,k)=(YiPsi(i,j+1,k)-YiPsi(i,j,k)) &
                       /(c_space(j+1) - c_space(j))
        enddo
        !
        !last point in c space, dYiPsidc(1:n_points_z,n_points_c)
        j=n_points_c
          dYiPsidc(i,j,k)=dYiPsidc(i,j-1,k) !CZ:half rectangle
        !first point in c space, dYiPsidc(1:n_points_z,n_points_c+1)
        j=n_points_c+1
          dYiPsidc(i,j,k)=dYiPsidc(i,1,k) !CZ:half rectangle
         !
        !start integration
        do j=1,n_points_c-2
          YiQ_int(i,k)=YiQ_int(i,k)-0.5*(dYiPsidc(i,j,k)*CDF_C(j) &
            +dYiPsidc(i,j+1,k)*CDF_C(j+1)) &
            *((c_space(j+2)+c_space(j+1)) &
            *0.5-(c_space(j+1)+c_space(j))*0.5)
        enddo
        !
        YiQ_int(i,k)=YiQ_int(i,k) &
          - dYiPsidc(i,n_points_c,k)*CDF_C(n_points_c) &
            *(c_space(n_points_c)-c_space(n_points_c-1))/2.0 &
          - dYiPsidc(i,n_points_c+1,k)*CDF_C(n_points_c+1) &
            *(c_space(2)-c_space(1))/2.0 &
          +YiPsi(i,n_points_c,k)
      enddo
      !
      !calculate dYiQdz(1:n_points_z-1)
      do i=1,n_points_z-1
        dYiQdz(i,k)=(YiQ_int(i+1,k) &
                -YiQ_int(i,k))/(z_space(i+1)-z_space(i))
      enddo
      !last point in z space, dYiQdz(n_points_z)
      i=n_points_z
        dYiQdz(i,k)=dYiQdz(i-1,k) !CZ:half rectangle
      !first point in z space, dYiQdz(n_points_z+1)
      i=n_points_z+1
        dYiQdz(i,k)=dYiQdz(1,k) !CZ:half rectangle

      !start integration
      do i=1,n_points_z-2       ! looping on z
        Yi_int(k)=Yi_int(k) &
              -0.5*(dYiQdz(i,k)*CDF_Z(i)+dYiQdz(i+1,k)*CDF_Z(i+1)) &
        *((z_space(i+2)+z_space(i+1))*0.5-(z_space(i+1)+z_space(i))*0.5)
      enddo

      Yi_int(k)=Yi_int(k) &
        +dYiQdz(n_points_z,k)*CDF_Z(n_points_z) &
           *(z_space(n_points_z)-z_space(n_points_z-1))/2.0 &
        +dYiQdz(n_points_z+1,k)*CDF_Z(n_points_z+1) &
           *(z_space(2)-z_space(1))/2.0 &
        +YiQ_int(n_points_z,k)
    enddo !k=1,nYis
    !
    do i=1,nScalars
      if(abs(y_int(i)).lt.1.e-20) y_int(i)=0.d0
    enddo
    !
    do i=1,nYis
      if(Yi_int(i).lt.1.e-20) Yi_int(i)=0.d0
    enddo
    !
  end subroutine int_point

! **********************************************************************
! **********************************************************************
  subroutine plackettFunc(c_point,z_point,alpha_z,beta_z,alpha_c, &
                          beta_c,theta,CDF_C,CDF_Z,plaF)
! **********************************************************************
! **********************************************************************
  !
  use func, only: betai
  !
  ! arguments
  real(8), intent(in) :: &
    c_point,z_point,alpha_z,beta_z,alpha_c,beta_c,theta
  real(8), intent(out) :: CDF_C,CDF_Z,plaF
  ! local data
  real(8) :: S,th
  !
  CDF_Z=betai(alpha_z,beta_z,z_point)
  CDF_C=betai(alpha_c,beta_c,c_point)
  !
  S=1.+(theta- 1.d0)*(CDF_Z+CDF_C)
  !
  if(theta.lt.1.e-5) then
    th=1.e-05
  else
    th=theta
  endif
  !
  if(abs(1.d0-th).lt.1.e-5) then
    plaF=1.d0
  else
    plaF=th*(1.d0+(th-1.d0)*(CDF_Z+CDF_C-2.0*(CDF_Z*CDF_C))) &
          /((S*S-(4.0*th*(th - 1.d0)*CDF_Z* CDF_C))**(3./2.))
  endif
  !
  end subroutine plackettFunc
  !
end module pdf
