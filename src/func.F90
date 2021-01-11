!+---------------------------------------------------------------------+
! Functions neccesary for computing Beta- and Gauss distribution
! Source: Numerical Recipes for FORTRAN Handbook
! http://www.library.cornell.edu/nr/cbookfpdf.html
!+---------------------------------------------------------------------+
module func
  !
  implicit none
  !
  contains
  !
  real(8) function betai(a,b,x)
    !
    ! arguments
    real(8),intent(in) :: a,b,x
    ! local data
    real(8) :: bt

    IF (x.EQ.0.0d0) THEN
      betai=0.0d0
      RETURN
    ELSEIF (x.EQ.1.0d0) THEN
      betai=1.0d0
      RETURN
    ELSE
      bt=exp(gammln(a+b)-gammln(a)-gammln(b) &
              +a*log(x)+b*log(1.0d0-x))
      IF (x.LT.(a+1.0d0)/(a+b+2.0d0)) THEN
         betai=bt*betacf(a,b,x)/a
         RETURN
      ELSE
         betai=1.0d0-bt*betacf(b,a,1.0d0-x)/b
         RETURN
      ENDIF
    ENDIF

  end function betai

! c***************************************************
  pure real(8) function  betacf(a,b,x)
! c---- Continued Fraction approach Beta function ----
    !
    ! arguments
    real(8),intent(in) :: a,b,x
    ! local data
    integer :: maxit,m,m2
    real(8) :: eps,fpmin,aa,c,d,del,h,qab,qap,qam
    parameter(maxit=100,eps=3.0e-7,fpmin=1.0e-30)

    qab=a+b
    qap=a+1.0d0
    qam=a-1.0d0
    c=1.0d0
    d=1.0d0-qab*x/qap
    IF (abs(d).LT.fpmin) d=fpmin
    d=1.0d0/d
    h=d
    DO m=1,maxit
      m2=2*m
      aa=m*(b-m)*x/((qam+m2)*(a+m2))
      d=1.0d0+aa*d
      IF(abs(d).LT.fpmin) d=fpmin
      c=1.0d0+aa/c
      IF(abs(c).LT.fpmin) c=fpmin
      d=1.0d0/d
      h=h*d*c
      aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
      d=1.0d0+aa*d
      IF(abs(d).LT.fpmin) d=fpmin
      c=1.0d0+aa/c
      IF(abs(c).LT.fpmin) c=fpmin
      d=1.0d0/d
      del=d*c
      h=h*del
      IF(abs(del-1.0d0).LT.eps) exit
    ENDDO
    betacf=h

  end function betacf

  ! c**************************************************
    real(8) function gammln(xx)
  ! c---- Logarithmic of Gamma function ---------------

      real(8), intent(in) :: xx

      INTEGER ::         j
      real(8) :: ser,stp,tmp,x,y,cof(6)

      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
        24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
        -.5395239384953d-5,2.5066282746310005d0/

      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      DO j=1,6
        y=y+1.0d0
        ser=ser+cof(j)/y
      ENDDO
      gammln=tmp+dlog(stp*ser/x)

    end FUNCTION gammln

! ! c**************************************************
  real(8) function gammp(a,x)
! c---- Series or Continued Fraction approach -------
    !
    ! arguments
    real(8),intent(in) :: a,x
    ! local data
    real(8) :: gamser,gammcf

       IF (x.LT.a+1.0d0) THEN
          CALL gser(gamser,a,x)
          gammp=gamser
       ELSE
          CALL gcf(gammcf,a,x)
          gammp=1.0d0-gammcf
       ENDIF

  end function gammp

! c**************************************************
  SUBROUTINE gser(gamser,a,x)
! c---- Series expansion approach -------------------
    real(8), intent(in) ::a,x
    real(8), intent(out) ::gamser

    INTEGER ::         maxit,n
    real(8) :: eps,ap,del,sm
    PARAMETER        (maxit=100,eps=3.0e-7)

    IF (x.LE.0.0d0) THEN
      gamser=0.0d0
      RETURN
    ENDIF
    ap=a
    sm=1.0d0/a
    del=sm
    DO n=1,maxit
      ap=ap+1
      del=del*x/ap
      sm=sm+del
      IF (abs(del).LT.abs(sm)*eps) exit
    ENDDO
    gamser=sm*exp(-x+a*log(x)-gammln(a))

  end SUBROUTINE gser

! c**************************************************
  SUBROUTINE gcf(gammcf,a,x)
! c---- Continued Fraction approach Gamma function --
    real(8), intent(in) ::a,x
    real(8), intent(out) ::gammcf

    INTEGER ::         maxit,i
    real(8) :: eps,del,fpmin,an,b,c,d,h
    PARAMETER        (maxit=100,eps=3.0e-7,fpmin=1.0e-30)

    b=x+1.0d0-a
    c=1.0d0/fpmin
    d=1.0d0/b
    h=d
    DO i=1,maxit
      an=-i*(i-a)
      b=b+2.0d0
      d=b+an*d
      IF (abs(d).LT.fpmin) d=fpmin
      c=b+an/c
      IF (abs(c).LT.fpmin) c=fpmin
      d=1.0d0/d
      del=d*c
      h=h*del
      IF (abs(del-1.0d0).LT.eps) exit
    ENDDO
    gammcf=exp(-x+a*log(x)-gammln(a))*h

  end SUBROUTINE gcf


  function locate(xx,n,x) result(j)

    integer, intent(in) :: n
    real(8), intent(in) :: xx(n),x

    integer :: i1,j

    j=1
    if(x.le.xx(1))then
      continue
    elseif(x.ge.xx(n))then
      j=n-1
    else
      do i1= 1,n-1
        if(x.gt.xx(i1).and.x.le.xx(i1+1))then
           j=i1
           exit
        endif
      enddo
    endif

  end function locate

! c************************************************************************

  function intfac(xx,xarray,loc_low) result(fac)

! c************************************************************************
    integer, intent(in) :: loc_low
    real(8), intent(in) :: xx,xarray(:)

    real(8) :: fac

    fac=0.
    if(xx.le.xarray(loc_low))then
       continue
    elseif(xx.ge.xarray(loc_low+1))then
       fac=1.
    else
       fac = (xx -xarray(loc_low))/(xarray(loc_low+1)-xarray(loc_low))
    endif

  end function intfac
!
! c************************************************************************
end module func
