************************************************************************
************************************************************************
      subroutine delta(z_mean,c_mean,z_space,c_space,y_int,sc_vals
     &                 ,Yi_int,Yi_vals)

************************************************************************
************************************************************************
      implicit none

      include 'integrate.inc'

      integer i,z_loc,c_loc

      double precision y_int(nScalars)
     &     ,sc_vals(n_points_z,n_points_c,nScalars)
     &     ,Yi_int(nYis),Yi_vals(n_points_z,n_points_c,nYis)
     &     ,c_mean,z_mean,c_space(n_points_c),z_space(n_points_z)
     &     ,z_fac,c_fac

      z_loc=1
      z_fac=0.
      do i=1,n_points_z-1
         if(z_mean.le.z_space(1))then
         elseif(z_mean.ge.z_space(i).and.z_mean.lt.z_space(i+1))then
            z_loc =i
            z_fac = (z_mean - z_space(i))/(z_space(i+1)-z_space(i))
            continue
         elseif(z_mean.ge.z_space(n_points_z))then
            z_loc = n_points_z-1
            z_fac=1.
            continue
         endif
      enddo

       c_loc=1
       c_fac=0.
      do i=1,n_points_c-1
         if(c_mean.le.c_space(1))then
         elseif(c_mean.ge.c_space(i).and.c_mean.lt.c_space(i+1))then
            c_loc =i
            c_fac = (c_mean - c_space(i))/(c_space(i+1)-c_space(i))
            continue
         elseif(c_mean.ge.c_space(n_points_c))then
            c_loc = n_points_c-1
            c_fac=1.
            continue
         endif
      enddo

      do i=1,nScalars
         y_int(i) = 0.
      enddo
      do i=1,nYis
         Yi_int(i) = 0.
      enddo

      if(c_loc.eq.0)c_loc=1
      if(z_loc.eq.0)z_loc=1

      do i=1,2     !    scalars !number of integral scalars
         y_int(i) = (1.-c_fac)*(z_fac*sc_vals(z_loc+1,c_loc,i)+
     &        (1.-z_fac)*(sc_vals(z_loc,c_loc,i))) +
     &        c_fac*(z_fac*sc_vals(z_loc+1,c_loc+1,i) +
     &        (1.-z_fac)*sc_vals(z_loc,c_loc+1,i))
      enddo

      do i=5,nScalars     !    scalars !number of integral scalars
         y_int(i) = (1.-c_fac)*(z_fac*sc_vals(z_loc+1,c_loc,i)+
     &        (1.-z_fac)*(sc_vals(z_loc,c_loc,i))) +
     &        c_fac*(z_fac*sc_vals(z_loc+1,c_loc+1,i) +
     &        (1.-z_fac)*sc_vals(z_loc,c_loc+1,i))
      enddo

      y_int(2)=y_int(2)/y_int(1) !Zhi:omega_c/rho
      y_int(3) = y_int(2)*c_mean ! gc_source
      y_int(4) = y_int(2)*z_mean ! gz_source

      do i=1,nYis     !    scalars !number of integral scalars Yi
         Yi_int(i) = (1.-c_fac)*(z_fac*Yi_vals(z_loc+1,c_loc,i)+
     &        (1.-z_fac)*(Yi_vals(z_loc,c_loc,i))) +
     &        c_fac*(z_fac*Yi_vals(z_loc+1,c_loc+1,i) +
     &        (1.-z_fac)*Yi_vals(z_loc,c_loc+1,i))
      enddo

      do i=1,nScalars
         if(abs(y_int(i)).lt.1.e-10)then
         y_int(i) = 0.
         endif
      enddo
      do i=1,nYis
         if(Yi_int(i).lt.1.e-10)then
         Yi_int(i) = 0.
         endif
      enddo

      return
      end

************************************************************************
************************************************************************
      subroutine cbeta(mean,g_var,space,z_mean,z_space,
     &       y_int,sc_vals,Yi_int,Yi_vals)
************************************************************************
************************************************************************

      implicit none

      include 'integrate.inc'

      integer i,loc,j

      double precision mean,g_var,y_int(nScalars)
     &     ,Yi_int(nYis),Yi_vals(n_points_z,n_points_c,nYis)
     &     ,sc_vals(n_points_z,n_points_c,nScalars)
     &     ,space(n_points_c),z_mean,z_space(n_points_z),fac
     &     ,sc_vals_int(n_points_c,nScalars)
     &     ,Yi_vals_int(n_points_c,nYis)

      double precision alpha_coeff,beta_coeff,alpha_c, beta_c,mean1,var1
      double precision betai,dYdc(nYis,n_points_c+1)
     &          ,cdf001(n_points_c),dsdc(nScalars,n_points_c+1)

      external betai  !CDF

      alpha_coeff(mean1,var1) = mean1*(((mean1*(1.0d0-mean1))/var1)-1.0)
      beta_coeff(mean1,var1)= (1.0d0-mean1)
     &                       *(((mean1*(1.0d0-mean1))/var1)-1.0)

      alpha_c=alpha_coeff(mean,g_var)
      beta_c =beta_coeff(mean,g_var)

      loc=1
      fac=0.
      do i=1,n_points_z-1
         if(z_mean.le.z_space(1))then
         elseif(z_mean.ge.z_space(i).and.z_mean.lt.z_space(i+1))then
            loc =i
            fac = (z_mean-z_space(i))/(z_space(i+1)-z_space(i))
            continue
         elseif(z_mean.ge.z_space(n_points_z))then
            loc = n_points_z-1
            fac=1.0
            continue
         endif
      enddo

c     interp sc_vals

      do i=1,n_points_c
         do j=1,nScalars
            sc_vals_int(i,j) = fac*sc_vals(loc+1,i,j)
     &           +((1.-fac)*sc_vals(loc,i,j))
         enddo
         do j=1,nYis
            Yi_vals_int(i,j) = fac*Yi_vals(loc+1,i,j)
     &           +((1.-fac)*Yi_vals(loc,i,j))
         enddo
      enddo


      do i=1,nScalars
         y_int(i) = 0.
      enddo
      do i=1,nYis
         Yi_int(i) = 0.
      enddo

c     ===== integration by part ====

       do i=1,n_points_c-1
            cdf001(i)  = betai(alpha_c,beta_c,(space(i)+space(i+1))/2.0)
       enddo

       DO j = 1,nScalars
        !! get derivative
        do i=1,n_points_c-1
           if((j.eq.1) .or. (j.ge.5))then !1:density 5:cpe 6:mwt 7:sum_dhfi
           dsdc(j,i)= (sc_vals_int(i+1,j)-sc_vals_int(i,j))
     &             /(space(i+1)-space(i))
           elseif(j.eq.2)then !2:omegac
           dsdc(j,i)= (sc_vals_int(i+1,2)/sc_vals_int(i+1,1)
     &               -sc_vals_int(i,2)/sc_vals_int(i,1))
     &             /(space(i+1)-space(i))
           elseif(j.eq.3)then !3:c*omegac
           dsdc(j,i)= (space(i+1)*sc_vals_int(i+1,2)/sc_vals_int(i+1,1)
     &               -space(i)*sc_vals_int(i,2)/sc_vals_int(i,1))
     &             /(space(i+1)-space(i))
           else !4:Z*omegac
           dsdc(j,i)= (z_mean*sc_vals_int(i+1,2)/sc_vals_int(i+1,1)
     &               -z_mean*sc_vals_int(i,2)/sc_vals_int(i,1))
     &             /(space(i+1)-space(i))
           endif
        enddo
        i=n_points_c
        dsdc(j,i) = dsdc(j,i-1)
        i=n_points_c+1
        dsdc(j,i) = dsdc(j,1)

        !! start integration
        do i=1,n_points_c-2 !n_points_c-1 points, so n_points_c-2 sections
           y_int(j) = y_int(j)
     &               - 0.5*(dsdc(j,i)*cdf001(i)+dsdc(j,i+1)*cdf001(i+1))
     &          *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
        enddo
           y_int(j) = y_int(j)
     &     - dsdc(j,n_points_c)*cdf001(n_points_c) !last half rectangle
     &       *(space(n_points_c)-space(n_points_c-1))/2.0 !first half rectangle ignored
     &     - dsdc(j,n_points_c+1)*cdf001(1) !first half rectangle
     &       *(space(2)-space(1))/2.0
     &              + sc_vals_int(n_points_c,j)
       ENDDO


       DO j = 1,nYis
        !! get derivative
        do i=1,n_points_c-1
           dYdc(j,i)= (Yi_vals_int(i+1,j)-Yi_vals_int(i,j))
     &             /(space(i+1)-space(i))
        enddo
        i=n_points_c
        dYdc(j,i) = dYdc(j,i-1)
        i=n_points_c+1
        dYdc(j,i) = dYdc(j,1)

        !! start integration
        do i=1,n_points_c-2 !n_points_c-1 points, so n_points_c-2 sections
           Yi_int(j) = Yi_int(j)
     &               - 0.5*(dYdc(j,i)*cdf001(i)+dYdc(j,i+1)*cdf001(i+1))
     &          *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
        enddo
           Yi_int(j) = Yi_int(j)
     &     - dYdc(j,n_points_c)*cdf001(n_points_c) !last half rectangle
     &       *(space(n_points_c)-space(n_points_c-1))/2.0 !first half rectangle ignored
     &     - dYdc(j,n_points_c+1)*cdf001(1) !first half rectangle
     &       *(space(2)-space(1))/2.0
     &              + Yi_vals_int(n_points_c,j)
       ENDDO

       do i=1,nScalars
         if(abs(y_int(i)).lt.1.e-10)then
         y_int(i) = 0.
         endif
      enddo
      do i=1,nYis
         if(Yi_int(i).lt.1.e-10)then
         Yi_int(i) = 0.
         endif
      enddo

      return
      end

************************************************************************
************************************************************************
      subroutine zbeta(mean,g_var,space,c_mean,c_space,
     &       y_int,sc_vals,Yi_int,Yi_vals)
************************************************************************
************************************************************************

      implicit none

      include 'integrate.inc'

      integer i,loc,j

      double precision mean,g_var,y_int(nScalars)
     &     ,Yi_int(nYis),Yi_vals(n_points_z,n_points_c,nYis)
     &     ,sc_vals(n_points_z,n_points_c,nScalars)
     &     ,space(n_points_z),c_mean,c_space(n_points_z),fac
     &     ,sc_vals_int(n_points_z,nScalars)
     &     ,Yi_vals_int(n_points_z,nYis)

      double precision alpha_coeff,beta_coeff,alpha_z, beta_z,mean1,var1
      double precision betai,dYdc(nYis,n_points_z+1)
     &          ,cdf001(n_points_z),dsdc(nScalars,n_points_z+1)

      external betai  !CDF

      alpha_coeff(mean1,var1) = mean1*(((mean1*(1.0d0-mean1))/var1)-1.0)
      beta_coeff(mean1,var1)= (1.0d0-mean1)
     &                       *(((mean1*(1.0d0-mean1))/var1)-1.0)

      alpha_z=alpha_coeff(mean,g_var)
      beta_z =beta_coeff(mean,g_var)

      loc=1
      fac=0.
      do i=1,n_points_c-1
         if(c_mean.le.c_space(1))then
         elseif(c_mean.ge.c_space(i).and.c_mean.lt.c_space(i+1))then
            loc =i
            fac = (c_mean-c_space(i))/(c_space(i+1)-c_space(i))
            continue
         elseif(c_mean.ge.c_space(n_points_c))then
            loc = n_points_c-1
            fac=1.0
            continue
         endif
      enddo

c     interp sc_vals

      do i=1,n_points_z
         do j=1,nScalars
            sc_vals_int(i,j) = fac*sc_vals(i,loc+1,j)
     &           +((1.-fac)*sc_vals(i,loc,j))
         enddo
         do j=1,nYis
            Yi_vals_int(i,j) = fac*Yi_vals(i,loc+1,j)
     &           +((1.-fac)*Yi_vals(i,loc,j))
         enddo
      enddo


      do i=1,nScalars
         y_int(i) = 0.
      enddo
      do i=1,nYis
         Yi_int(i) = 0.
      enddo

c     ===== integration by part ====

       open(unit=29, file='tmp.txt')
       space(1) = smaller
       do i=1,n_points_z-1
            cdf001(i)  = betai(alpha_z,beta_z,(space(i)+space(i+1))/2.0)
       enddo

       DO j = 1,nScalars
        !! get derivative
        do i=1,n_points_z-1
           if((j.eq.1) .or. (j.ge.5))then !1:density 5:cpe 6:mwt 7:sum_dhfi
           dsdc(j,i)= (sc_vals_int(i+1,j)-sc_vals_int(i,j))
     &             /(space(i+1)-space(i))
           elseif(j.eq.2)then !2:omegac
           dsdc(j,i)= (sc_vals_int(i+1,2)/sc_vals_int(i+1,1)
     &               -sc_vals_int(i,2)/sc_vals_int(i,1))
     &             /(space(i+1)-space(i))
           elseif(j.eq.3)then !3:c*omegac
           dsdc(j,i)= (c_mean*sc_vals_int(i+1,2)/sc_vals_int(i+1,1)
     &               -c_mean*sc_vals_int(i,2)/sc_vals_int(i,1))
     &             /(space(i+1)-space(i))
           else !4:Z*omegac
           dsdc(j,i)= (space(i+1)*sc_vals_int(i+1,2)/sc_vals_int(i+1,1)
     &               -space(i)*sc_vals_int(i,2)/sc_vals_int(i,1))
     &             /(space(i+1)-space(i))
           write(29,*) c_mean,dsdc(1,i),sc_vals_int(i,1),dsdc(5,i)
     & ,sc_vals_int(i,5)
           endif
        enddo
        i=n_points_z
        dsdc(j,i) = dsdc(j,i-1)
        i=n_points_z+1
        dsdc(j,i) = dsdc(j,1)

        !! start integration
        if((j.le.nScalars))then
        do i=1,n_points_z-2 !n_points_z-1 points, so n_points_z-2 sections
           y_int(j) = y_int(j)
     &               - 0.5*(dsdc(j,i)*cdf001(i)+dsdc(j,i+1)*cdf001(i+1))
     &          *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
        enddo
           y_int(j) = y_int(j)
     &     - dsdc(j,n_points_z)*cdf001(n_points_z) !last half rectangle
     &       *(space(n_points_z)-space(n_points_z-1))/2.0 !first half rectangle ignored
     &     - dsdc(j,n_points_z+1)*cdf001(1) !first half rectangle
     &       *(space(2)-space(1))/2.0
     &              + sc_vals_int(n_points_z,j)
        else !for Yis
        do i=1,n_points_z-2 !n_points_z-1 points, so n_points_z-2 sections
           Yi_int(j) = y_int(j)
     &               - 0.5*(dsdc(j,i)*cdf001(i)+dsdc(j,i+1)*cdf001(i+1))
     &          *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
        enddo
           Yi_int(j) = y_int(j)
     &     - dsdc(j,n_points_z)*cdf001(n_points_z) !last half rectangle
     &       *(space(n_points_z)-space(n_points_z-1))/2.0 !first half rectangle ignored
     &     - dsdc(j,n_points_z+1)*cdf001(1) !first half rectangle
     &       *(space(2)-space(1))/2.0
     &              + sc_vals_int(n_points_z,j)
        endif
       ENDDO


       DO j = 1,nYis
        !! get derivative
        do i=1,n_points_z-1
           dYdc(j,i)= (Yi_vals_int(i+1,j)-Yi_vals_int(i,j))
     &             /(space(i+1)-space(i))
        enddo
        i=n_points_z
        dYdc(j,i) = dYdc(j,i-1)
        i=n_points_z+1
        dYdc(j,i) = dYdc(j,1)

        !! start integration
        do i=1,n_points_z-2 !n_points_z-1 points, so n_points_z-2 sections
           Yi_int(j) = Yi_int(j)
     &               - 0.5*(dYdc(j,i)*cdf001(i)+dYdc(j,i+1)*cdf001(i+1))
     &          *((space(i+2)+space(i+1))*0.5-(space(i+1)+space(i))*0.5)
        enddo
           Yi_int(j) = Yi_int(j)
     &     - dYdc(j,n_points_z)*cdf001(n_points_z) !last half rectangle
     &       *(space(n_points_z)-space(n_points_z-1))/2.0 !first half rectangle ignored
     &     - dYdc(j,n_points_z+1)*cdf001(1) !first half rectangle
     &       *(space(2)-space(1))/2.0
     &              + Yi_vals_int(n_points_z,j)
       ENDDO

      do i=1,nScalars
         if(abs(y_int(i)).lt.1.e-10)then
         y_int(i) = 0.
         endif
      enddo
      do i=1,nYis
         if(Yi_int(i).lt.1.e-10)then
         Yi_int(i) = 0.
         endif
      enddo

      return
      end

************************************************************************
************************************************************************
      subroutine int_point(z_mean,c_mean,c_var,z_var,co_var,y_int
     &     ,z_space,c_space,sc_vals,Yi_int,Yi_vals,theta)

************************************************************************
************************************************************************
      implicit none
!      include 'data.inc'
      include 'integrate.inc'

      integer i,j,k,dum

      double precision c_mean, z_mean, c_var, z_var, co_var
     &     ,c_space(n_points_c), z_space(n_points_z)
     &     ,alpha_c,alpha_z,beta_c,beta_z, rho,theta
     &     ,alpha_coeff,beta_coeff
     &     ,mean,var,y_int(nScalars)

      double precision sc_vals(n_points_z,n_points_c,nScalars)
     &     ,Yi_int(nYis),Yi_vals(n_points_z,n_points_c,nYis)
     &  ,CDF_C(n_points_c+1),CDF_Z(n_points_z+1)
     &  ,plaF(n_points_z+1,n_points_c+1)
     &  ,dPsidc(n_points_z+1,n_points_c+1,nScalars)
     &  ,Psi(n_points_z+1,n_points_c+1,nScalars)
     &  ,Q_int(n_points_z+1,nScalars)
     &  ,dQdz(n_points_z+1,nScalars)
     &  ,YiPsi(n_points_z+1,n_points_c+1,nYis)
     &  ,YiQ_int(n_points_z+1,nYis)
     &  !,dYiQdz(n_points_z+1,nYis)
     &  !,dYiPsidc(n_points_z+1,n_points_c+1,nYis)

      alpha_coeff(mean,var) = mean*(((mean*(1.0d0-mean))/var)-1.0)
      beta_coeff(mean,var)=(1.0d0-mean)*(((mean*(1.0d0-mean))/var)-1.0)

      rho = co_var/(dsqrt(z_var)*dsqrt(c_var))
      if(abs(rho).ge.1.0e03)then
         theta=1.
      endif

      alpha_z = alpha_coeff(z_mean,z_var)
      alpha_c = alpha_coeff(c_mean,c_var)
      beta_z = beta_coeff(z_mean,z_var)
      beta_c = beta_coeff(c_mean,c_var)
!      write(*,*) alpha_z,alpha_c,beta_z,beta_c,c_mean,c_var

      do k=1,nScalars
      y_int(k) = 0.0
        do i=1,n_points_z
        Q_int(i,k) = 0.0
        enddo
      enddo

      do k=1,nYis
      Yi_int(k) = 0.0
        do i=1,n_points_z
        YiQ_int(i,k) = 0.0
        enddo
      enddo

c-----assign values for Psi and obtain CDFs------------------------
            do dum=1,1
      !++! main Z array
      do i=1,n_points_z-1

        !!!for Psi(1:n_points_z-1,1:n_points_c-1)
        do j=1,n_points_c-1
        call plackettFunc((c_space(j)+c_space(j+1))/2.0,
     &                    (z_space(i)+z_space(i+1))/2.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))
          do k=5,nScalars
           Psi(i,j,k) = sc_vals(i,j,k)*plaF(i,j)
          enddo
           Psi(i,j,2) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
           Psi(i,j,3) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *c_space(j)
           Psi(i,j,4) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *z_space(i)
          do k=1,nYis
           YiPsi(i,j,k) = Yi_vals(i,j,k)*plaF(i,j)
          enddo
        enddo

        !for Psi(1:n_points_z-1,n_points_c)
        j=n_points_c
        call plackettFunc((c_space(j-1)
     &                     +3.0*c_space(j))/4.0, !CZ:half trapz
     &                    (z_space(i)+z_space(i+1))/2.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))
        do k=5,nScalars
        Psi(i,j,k)=sc_vals(i,j,k)*plaF(i,j)
        enddo
        Psi(i,j,2) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
        Psi(i,j,3) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *c_space(j)
        Psi(i,j,4) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *z_space(i)
        do k=1,nYis
        YiPsi(i,j,k)=Yi_vals(i,j,k)*plaF(i,j)
        enddo

        !!for Psi(1:n_points_z-1,n_points_c+1) first point in c space
        j=n_points_c+1
        call plackettFunc((3.0*c_space(1)
     &                     +c_space(2))/4.0, !CZ:half trapz
     &                    (z_space(i)+z_space(i+1))/2.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))
      enddo

      !++! last point Z array
      i=n_points_z

        !!!for Psi(n_points_z,1:n_points_c-1)
        do j=1,n_points_c-1
        call plackettFunc((c_space(j)+c_space(j+1))/2.0,
     &                    (z_space(i-1)
     &                     +3.0*z_space(i))/4.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))
          do k=5,nScalars
          Psi(i,j,k) = sc_vals(i,j,k)*plaF(i,j)
          enddo
        Psi(i,j,2) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
        Psi(i,j,3) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *c_space(j)
        Psi(i,j,4) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *z_space(i)
          do k=1,nYis
          YiPsi(i,j,k) = Yi_vals(i,j,k)*plaF(i,j)
          enddo
        enddo

        !!!!for Psi(n_points_z,1:n_points_c)
        j=n_points_c
        call plackettFunc((c_space(j-1)
     &                     +3.0*c_space(j))/4.0, !CZ:half trapz
     &                    (z_space(i-1)
     &                     +3.0*z_space(i))/4.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))
        do k=5,nScalars
        Psi(i,j,k)=sc_vals(i,j,k)*plaF(i,j)
        enddo
        Psi(i,j,2) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
        Psi(i,j,3) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *c_space(j)
        Psi(i,j,4) = sc_vals(i,j,2)*plaF(i,j)/sc_vals(i,j,1)
     &                 *z_space(i)
        do k=1,nYis
        YiPsi(i,j,k)=Yi_vals(i,j,k)*plaF(i,j)
        enddo

        j=n_points_c+1
        !!!!for Psi(n_points_z,n_points_c+1) first point in c space
        call plackettFunc((3.0*c_space(1)
     &                     +c_space(2))/4.0, !CZ:half trapz
     &                    (z_space(i-1)
     &                     +z_space(i))/4.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))

      !++! first point Z array
      i = n_points_z+1
        !!!for Psi(n_points_z+1,1:n_points_c-1)
        do j=1,n_points_c-1
        call plackettFunc((c_space(j)+c_space(j+1))/2.0,
     &                    (3.0*z_space(1)
     &                     +z_space(2))/4.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))
         enddo
        !!!for Psi(n_points_z+1,n_points_c)
        j=n_points_c
        call plackettFunc((c_space(j-1)
     &                     +3.0*c_space(j))/4.0, !CZ:half trapz
     &                    (3.0*z_space(1)
     &                     +z_space(2))/4.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))
        !!!!for Psi(n_points_z+1,n_points_c+1) first point in c space
        j=n_points_c+1
        call plackettFunc((3.0*c_space(1)
     &                     +c_space(2))/4.0, !CZ:half trapz
     &                    (3.0*z_space(1)
     &                     +z_space(2))/4.0,
     &                    alpha_z,beta_z,alpha_c,beta_c,
     &                    theta,CDF_C(j),CDF_Z(i),
     &                    plaF(i,j))

              enddo

c-----2:omegac 3:c*omegac 4:z*omegac-----------------------------------

              do k = 2,3
      !calculate dPsidc(1:n_points_z,1:n_points_c-1)
      do i = 1,n_points_z       ! looping on z

         do j = 1,n_points_c-1    ! looping on c
           dPsidc(i,j,k) = (Psi(i,j+1,k)-Psi(i,j,k)+1.e-15)
     &                    /(c_space(j+1) - c_space(j))
         enddo
         !last point in c space, dPsidc(1:n_points_z,n_points_c)
         j = n_points_c
         dPsidc(i,j,k)=dPsidc(i,j-1,k) !CZ:half rectangle
         !first point in c space, dPsidc(1:n_points_z,n_points_c+1)
         j = n_points_c+1
         dPsidc(i,j,k)=dPsidc(i,1,k) !CZ:half rectangle

         !start integration
         do j = 1,n_points_c-2
           Q_int(i,k) = Q_int(i,k)
     &        -0.5*(dPsidc(i,j,k)*CDF_C(j)+dPsidc(i,j+1,k)*CDF_C(j+1))
     &            *((c_space(j+2)+c_space(j+1))
     &               *0.5-(c_space(j+1)+c_space(j))*0.5)
         enddo

         Q_int(i,k) = Q_int(i,k)
     &     - dPsidc(i,n_points_c,k)*CDF_C(n_points_c) !last half rectangle
     &       *(c_space(n_points_c)-c_space(n_points_c-1))/2.0 !first half rectangle ignored
     &     - dPsidc(i,n_points_c+1,k)*CDF_C(n_points_c+1) !first half rectangle
     &       *(c_space(2)-c_space(1))/2.0
     &              + Psi(i,n_points_c,k)
      enddo

      !calculate dQdz(1:n_points_z-1)
      do i = 1,n_points_z-1
      dQdz(i,k) = (Q_int(i+1,k)
     &            -Q_int(i,k)+1.e-15)/(z_space(i+1)-z_space(i))
      enddo
      !last point in z space, dQdz(n_points_z)
      i = n_points_z
      dQdz(i,k) = dQdz(i-1,k) !CZ:half rectangle
      !first point in z space, dQdz(n_points_z+1)
      i = n_points_z+1
      dQdz(i,k) = dQdz(1,k) !CZ:half rectangle

      !start integration
      do i = 1,n_points_z-2       ! looping on z
      y_int(k) = y_int(k)
     &         -0.5*(dQdz(i,k)*CDF_Z(i)+dQdz(i+1,k)*CDF_Z(i+1))
     &*((z_space(i+2)+z_space(i+1))*0.5-(z_space(i+1)+z_space(i))*0.5)
      enddo

      y_int(k) = y_int(k)
     &    + dQdz(n_points_z,k)*CDF_Z(n_points_z) !last half rectangle
     &      *(z_space(n_points_z)-z_space(n_points_z-1))/2.0 !first half rectangle ignored
     &    + dQdz(n_points_z+1,k)*CDF_Z(n_points_z+1) !last half rectangle
     &      *(z_space(2)-z_space(1))/2.0 !first half rectangle ignored
     &          + Q_int(n_points_z,k)

              enddo

c-----5:Cp 6:MW 7:formEnthal 8:T-----------------------------------
              do k = 5,nScalars
      !calculate dPsidc(1:n_points_z,1:n_points_c-1)
      do i = 1,n_points_z       ! looping on z

         do j = 1,n_points_c-1    ! looping on c
           dPsidc(i,j,k) = (Psi(i,j+1,k)-Psi(i,j,k))
     &                    /(c_space(j+1) - c_space(j))
         enddo
         !last point in c space, dPsidc(1:n_points_z,n_points_c)
         j = n_points_c
         dPsidc(i,j,k)=dPsidc(i,j-1,k) !CZ:half rectangle
         !first point in c space, dPsidc(1:n_points_z,n_points_c+1)
         j = n_points_c+1
         dPsidc(i,j,k)=dPsidc(i,1,k) !CZ:half rectangle

         !start integration
         do j = 1,n_points_c-2
           Q_int(i,k) = Q_int(i,k)
     &        -0.5*(dPsidc(i,j,k)*CDF_C(j)+dPsidc(i,j+1,k)*CDF_C(j+1))
     &            *((c_space(j+2)+c_space(j+1))
     &               *0.5-(c_space(j+1)+c_space(j))*0.5)
         enddo

         Q_int(i,k) = Q_int(i,k)
     &     - dPsidc(i,n_points_c,k)*CDF_C(n_points_c) !last half rectangle
     &       *(c_space(n_points_c)-c_space(n_points_c-1))/2.0
     &     - dPsidc(i,n_points_c+1,k)*CDF_C(n_points_c+1) !first half rectangle
     &       *(c_space(2)-c_space(1))/2.0
     &              + Psi(i,n_points_c,k)
      enddo

      !ccalculate dQdz(1:n_points_z-1)
      do i = 1,n_points_z-1
      dQdz(i,k) = (Q_int(i+1,k)
     &            -Q_int(i,k))/(z_space(i+1)-z_space(i))
      enddo
      !last point in z space, dQdz(n_points_z)
      i = n_points_z
      dQdz(i,k) = dQdz(i-1,k) !CZ:half rectangle
      !first point in z space, dQdz(n_points_z+1)
      i = n_points_z+1
      dQdz(i,k) = dQdz(1,k) !CZ:half rectangle

      !start integration
      do i = 1,n_points_z-2       ! looping on z
      y_int(k) = y_int(k)
     &         -0.5*(dQdz(i,k)*CDF_Z(i)+dQdz(i+1,k)*CDF_Z(i+1))
     &*((z_space(i+2)+z_space(i+1))*0.5-(z_space(i+1)+z_space(i))*0.5)
      enddo

      y_int(k) = y_int(k)
     &    + dQdz(n_points_z,k)*CDF_Z(n_points_z) !last half rectangle
     &      *(z_space(n_points_z)-z_space(n_points_z-1))/2.0
     &    + dQdz(n_points_z+1,k)*CDF_Z(n_points_z+1) !first half rectangle
     &      *(z_space(2)-z_space(1))/2.0
     &          + Q_int(n_points_z,k)

              enddo

c-----Yis-----------------------------------
!              do k = 1,nYis
!       !calculate dYiPsidc(1:n_points_z,1:n_points_c-1)
!       do i = 1,n_points_z       ! looping on z
!
!          do j = 1,n_points_c-1    ! looping on c
!            dYiPsidc(i,j,k) = (YiPsi(i,j+1,k)-YiPsi(i,j,k))
!      &                    /(c_space(j+1) - c_space(j))
!          enddo
!          !last point in c space, dYiPsidc(1:n_points_z,n_points_c)
!          j = n_points_c
!          dYiPsidc(i,j,k)=dYiPsidc(i,j-1,k) !CZ:half rectangle
!          !first point in c space, dYiPsidc(1:n_points_z,n_points_c+1)
!          j = n_points_c+1
!          dYiPsidc(i,j,k)=dYiPsidc(i,1,k) !CZ:half rectangle
!
!          !start integration
!          do j = 1,n_points_c-2
!            YiQ_int(i,k) = YiQ_int(i,k)
!      &    -0.5*(dYiPsidc(i,j,k)*CDF_C(j)+dYiPsidc(i,j+1,k)*CDF_C(j+1))
!      &            *((c_space(j+2)+c_space(j+1))
!      &               *0.5-(c_space(j+1)+c_space(j))*0.5)
!          enddo
!
!          YiQ_int(i,k) = YiQ_int(i,k)
!      &     - dYiPsidc(i,n_points_c,k)*CDF_C(n_points_c) !last half rectangle
!      &       *(c_space(n_points_c)-c_space(n_points_c-1))/2.0 !first half rectangle ignored
!      &     - dYiPsidc(i,n_points_c+1,k)*CDF_C(n_points_c+1) !first half rectangle
!      &       *(c_space(2)-c_space(1))/2.0
!      &              + YiPsi(i,n_points_c,k)
!       enddo
!
!       !calculate dYiQdz(1:n_points_z-1)
!       do i = 1,n_points_z-1
!       dYiQdz(i,k) = (YiQ_int(i+1,k)
!      &            -YiQ_int(i,k))/(z_space(i+1)-z_space(i))
!       enddo
!       !last point in z space, dYiQdz(n_points_z)
!       i = n_points_z
!       dYiQdz(i,k) = dYiQdz(i-1,k) !CZ:half rectangle
!       !first point in z space, dYiQdz(n_points_z+1)
!       i = n_points_z+1
!       dYiQdz(i,k) = dYiQdz(1,k) !CZ:half rectangle
!
!       !start integration
!       do i = 1,n_points_z-2       ! looping on z
!       Yi_int(k) = Yi_int(k)
!      &         -0.5*(dYiQdz(i,k)*CDF_Z(i)+dYiQdz(i+1,k)*CDF_Z(i+1))
!      &*((z_space(i+2)+z_space(i+1))*0.5-(z_space(i+1)+z_space(i))*0.5)
!       enddo
!
!       Yi_int(k) = Yi_int(k)
!      &    + dYiQdz(n_points_z,k)*CDF_Z(n_points_z) !last half rectangle
!      &      *(z_space(n_points_z)-z_space(n_points_z-1))/2.0 !first half rectangle ignored
!      &    + dYiQdz(n_points_z+1,k)*CDF_Z(n_points_z+1) !last half rectangle
!      &      *(z_space(2)-z_space(1))/2.0 !first half rectangle ignored
!      &          + YiQ_int(n_points_z,k)
!
!               enddo

      do i=1,nScalars
         if(abs(y_int(i)).lt.1.e-10)then
         y_int(i) = 0.
         endif
      enddo
      do i=1,nYis
         if(Yi_int(i).lt.1.e-10)then
         Yi_int(i) = 0.
         endif
      enddo

      return
      end

************************************************************************
c
      subroutine plackettFunc(c_point,z_point,alpha_z,beta_z,alpha_c
     &                        ,beta_c,theta,CDF_C,CDF_Z,plaF)
************************************************************************
c     subroutine to calculate the plackett function

      implicit none
      include 'integrate.inc'

      double precision c_point,z_point,alpha_z,beta_z,alpha_c,
     &                 beta_c,theta,CDF_C,plaF,CDF_Z,S,betai

      external betai

      CDF_Z = betai(alpha_z,beta_z,z_point)
      CDF_C = betai(alpha_c,beta_c,c_point)
!       write(6,*) CDF_Z(z_point),CDF_C(c_point)

      S = 1. + (theta - 1.) * (CDF_Z + CDF_C)

      if(theta.lt.1.0e-05)then
      theta=1.0e-05
      endif

      if(abs(1.0-theta).lt.1.0e-05)then
      plaF = 1.0

      else
      plaF = theta *
     &    (1.0+(theta-1.0)*(CDF_Z+CDF_C
     &    -2*(CDF_Z*CDF_C)))
     &       /((S*S -
     &     (4.0 * theta * (theta - 1.0) * CDF_Z
     &    * CDF_C))**(3./2.))
      endif

      return
      end

c************************************************************************

      subroutine locate2(xx,n,x,j)

c************************************************************************

      implicit none
      integer n,j,i1
      double precision xx(n),x

c$$$      if(flag.eq.2)then
c$$$      write(6,*) xx
c$$$      stop
c$$$      endif

      if(x.le.xx(1))then
         j=1
         goto 10
      elseif(x.ge.xx(n))then
         j=n-1
         goto 10
      else
         do i1= 1,n-1
            if(x.gt.xx(i1).and.x.le.xx(i1+1))then
               j=i1
               goto 10
            endif
         enddo
      endif


 10   return
      end

c************************************************************************

      subroutine intfac(xx,low,high,fac)

c************************************************************************

      implicit none

      double precision xx,low,high,fac

      if(xx.le.low)then
         fac=0.
      elseif(xx.ge.high)then
         fac=1.
      else
         fac = (xx -low)/(high-low)
      endif

      return
      end

c************************************************************************
