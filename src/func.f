c--- Functions neccesary for computing Beta- and Gauss distribution --------
c    Source: Numerical Recipes for FORTRAN Handbook
c    http://www.library.cornell.edu/nr/cbookfpdf.html
c---------------------------------------------------------------------------
       FUNCTION betai(a,b,x)

c---- Incomplete Beta function ---------------------
       DOUBLE PRECISION   betai,a,b,x
       DOUBLE PRECISION   bt,betacf,gammln


       IF (x.EQ.0.0d0) THEN
          betai=0.0d0
          RETURN
       ELSEIF (x.EQ.1.0d0) THEN
          betai=1.0d0
          RETURN
       ELSE
          bt=exp(gammln(a+b)-gammln(a)-gammln(b)
     v          +a*log(x)+b*log(1.0d0-x))
          IF (x.LT.(a+1.0d0)/(a+b+2.0d0)) THEN
             betai=bt*betacf(a,b,x)/a
             RETURN
          ELSE
             betai=1.0d0-bt*betacf(b,a,1.0d0-x)/b
             RETURN
          ENDIF
       ENDIF

       END

c***************************************************
       FUNCTION  betacf(a,b,x)

c---- Continued Fraction approach Beta function ----
       INTEGER           maxit,m,m2
       DOUBLE PRECISION  betacf,a,b,x,eps,fpmin
       DOUBLE PRECISION  aa,c,d,del,h,qab,qap,qam
       PARAMETER         (maxit=100,eps=3.0e-7,fpmin=1.0e-30)

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
          IF(abs(del-1.0d0).LT.eps) GOTO 1
       ENDDO
1      betacf=h

       RETURN
       END

c**************************************************
       FUNCTION gammp(a,x)

c---- Series or Continued Fraction approach -------
       DOUBLE PRECISION gammp,a,x,gammcf,gamser

       IF (x.LT.a+1.0d0) THEN
          CALL gser(gamser,a,x)
          gammp=gamser
       ELSE
          CALL gcf(gammcf,a,x)
          gammp=1.0d0-gammcf
       ENDIF

       RETURN
       END

c**************************************************
       SUBROUTINE gser(gamser,a,x)

c---- Series expansion approach -------------------
       INTEGER          maxit,n
       DOUBLE PRECISION a,x,gamser,eps,ap,del,sm,gammln
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
          IF (abs(del).LT.abs(sm)*eps) GOTO 1
       ENDDO
1      gamser=sm*exp(-x+a*log(x)-gammln(a))

       RETURN
       END

c**************************************************
       SUBROUTINE gcf(gammcf,a,x)

c---- Continued Fraction approach Gamma function --
       INTEGER          maxit,i
       DOUBLE PRECISION a,x,gammcf,eps,del,gammln,fpmin
       DOUBLE PRECISION an,b,c,d,h
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
          IF (abs(del-1.0d0).LT.eps) GOTO 1
       ENDDO
1      gammcf=exp(-x+a*log(x)-gammln(a))*h

       RETURN
       END

c************************************************************************
c$$$      function betapdf(a,b,x)
c$$$
c$$$      double precision betapdf, a, b, x, fac
c$$$
c$$$      fac = gammln(a+b) - gammln(a) - gammln(b)
c$$$
c$$$      betapdf = exp((a-1.)*log(x)+(b-1.)*log(1.-x)+fac)
c$$$
c$$$      return
c$$$      end


c**************************************************
       FUNCTION gammln(xx)

c---- Logarithmic of Gamma function ---------------
       INTEGER          j
       DOUBLE PRECISION gammln,xx,ser,stp,tmp,x,y,cof(6)

       SAVE cof,stp
       DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     v  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     v  -.5395239384953d-5,2.5066282746310005d0/
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

       RETURN
       END

c************************************************************************
c$$$      function alpha_coeff(mean,var)
c$$$
c$$$c     return the alpha coefficient for a beta function
c$$$      double precision alpha_coeff,mean,var
c$$$
c$$$      alpha_coeff = mean*(((mean*(1.-mean))/var)-1.0)
c$$$
c$$$      return
c$$$      end
c$$$
c$$$c************************************************************************
c$$$      function beta_coeff(mean,var)
c$$$
c$$$c     return the beta coefficient for a beta function
c$$$      double precision beta_coeff,mean,var
c$$$
c$$$      beta_coeff = (1.-mean)*(((mean*(1.-mean))/var)-1.0)
c$$$
c$$$      return
c$$$      end

c************************************************************************
