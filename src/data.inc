c--- Initializing paramaters for sizes matrices ------------------
      INTEGER           Numflam,Npoint,Nspec,Nv,Ne
      INTEGER           Numflami,Npointi,Npointj,Nspeci,Nvi
      PARAMETER         (Numflami=150,Npointi=2000,Npointj=500)
      PARAMETER         (Nspeci=200,Nvi=Nspeci+20,Ne=4)
      
c--- Integers ----------------------------------------------------
      INTEGER           icv1,icv2,ivcv1,ivcv2,isvcv2,iCO,iNO,isNO
      INTEGER           iC2H2,iC6H6,iO2,iO,iOH,iH2O,ich4,in2,ico2
      INTEGER		ih2,ic3h8,isch4,isco2,ih,ihcn,isc3h8
      INTEGER           idens,ipress,itemp,imeanmass,iscv2,
     v                  imassflow,icp,isco,ienth
      INTEGER           Ksi(Ne,Nspeci),ic12h26
      COMMON /lamflamdatai/Numflam,Npoint,Nspec,Nv,
     v   icv1,icv2,ivcv1,ivcv2,idens,ipress,itemp,imassflow,
     v   icp,imeanmass,iscv2,isvcv2,Ksi,iCO,iNO,isNO,
     v   iC2H2,iC6H6,iO2,iO,iOH,iH2O,iCH4,iC3H8,iCO2,
     v	 ih2,in2,isch4,isco2,ih,ihcn,isc3h8,ic12h26,isco,
     v   ienth

c--- Reals -------------------------------------------------------      
      DOUBLE PRECISION  MaxMF,MinMF,Runiv,Zstoich,Pi
      DOUBLE PRECISION  Ele(Ne,Numflami),Ycvmax(Numflami)
      DOUBLE PRECISION  Mixf(Numflami),Atm(Nspeci)
      DOUBLE PRECISION  Purefuel(Nvi+Nspeci+3),
     v                  Pureoxi(Nvi+Nspeci+3)
      DOUBLE PRECISION  Yread(Nvi+Nspeci+3,Numflami*Npointj)
      DOUBLE PRECISION  Sootrates(6,Numflami*Npointj), 
     v                  Ymax(Nvi+Nspeci+3),max_co2,fmix_cut
      DOUBLE PRECISION co2_max(Numflami)


      PARAMETER         (Runiv=8.314d0,Pi=3.14159265358979d0)
      COMMON /lamflamdatar/MaxMF,MinMF,yread,Ele,Ycvmax, 
     v   Atm,Mixf,purefuel,pureoxi,Zstoich,Sootrates,
     v	 co2_max

c--- Characters --------------------------------------------------      
      CHARACTER*10      CV1,CV2,VCV1,VCV2,SVCV2
      CHARACTER*15      Name(Nvi+Nspeci+3)
      COMMON /lamflamdatac/CV1,CV2,VCV1,VCV2,SVCV2,Name
c---Logicals----
      Logical Premixed
      COMMON /lamflamdatalog/Premixed      

      Logical Soot
      Integer nsoot
      COMMON /SOOTORNOSOOT/Soot,nSoot
    
