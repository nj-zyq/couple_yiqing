      PROGRAM COUPLE
c
c     This program is the driver for a range dependent coupled
c     normal mode model for acoustic propagation in a sequence of
c     environmental regions with depth dependent variations of
c     the sound speed, attenuation, and density. The environmental
c     regions consist of a water layer and a bottom sediment layer,
c     separated by a linearly sloping water-sediment interface.
c
c     A pressure release boundary is imposed at the surface of the
c     water and at the greatest depth in the bottom sediment.
c     the water depth may vary with range but the total thickness
c     of water and bottom sediment layers is fixed.
c
c     The linearly slopping bottom is approximated by a sequence
c     of locally flat computational regions that give rise to
c     the name stepwise coupled modes. Mode coupling occurs,
c     discretely, at the vertical interfaces between the regions.
c
c     The depth variations of the sound speed and attenuation are
c     assumed to be such that the complex wave number squared is
c     continuous and piecewise linear, in both layers. The density
c     is assumed to be continuous and piecewise exponential, in both
c     layers. There is usually a discontinuity at the water-sediment
c     interface, i.e., at the bottom. This environment is referred 
c     to as the complex foreground model. The normal modes are
c     needed, for the complex foreground model, in each of the
c     locally flat computational regions. These normal modes, or
c     eigenfunctions, are approximated with the Galerkin method.  
c
c     The Galerkin method is employed, using a basis of real
c     eigenfunctions from a lossless background model consisting
c     of a homogeneous water layer over a homogeneous bottom
c     sediment layer. The density discontinuity at the water-
c     sediment interface in the background model is chosen to
c     match the density discontinuity, at the water-sediment
c     interface, in the complex foreground model. The 
c     eigenfunctions of the complex foreground model are
c     approximated as a finite sum of the real eigenfunctions
c     of the background model, with complex coefficients.
c     consequently, the differential eigenvalue problem, for
c     the eigenfunctions, is replace by a finite dimensional
c     generalized matrix eigenvalue probelm. The complex eigenvalues
c     are approximated by the matrix eigenvalues. The complex
c     eigenfunctons, or modes, are constructed from the matrix
c     eigenvectors. The Galerkin approximation procedure is
c     carried out in subroutine GALRKN
c
c     In the case that the fully elliptic two way solution is
c     computed (IOUT=0), the decoupling algorithm is used to
c     solve the matrix two-point boundary value problem that
c     occurs in the stepwise coupled mode formulation. This
c     eliminates the numerical problem inherent in the shooting
c     method that was initially used.
c
c     The decoupling algorithm was modified to compute a
c     fundamental matrix solution on 9/20/04. 
c
c     Accounted for polarity change of the modes in the adiabatic
c     approximation on 8/3/05.
c
c     Modified subroutine CPRES to allow a long first region 8/18/06.
c
c     Converted to the free FORTRAN compiler g95 on Jan. 31, 2007,
c     which can write scratch files in excess of 3gb.
c     Wrote T on SCRACH41.DAT - SCRACH44.DAT to reduce file size.
c     Wrote UC on SCRACH23.DAT to reduce file size.
c     Wrote UB on SCRACH24.DAT to reduce file size.
c     Wrote Y1 on SCRACH31.DAT to reduce file size.
c     Wrote Y3 on SCRACH33.DAT to reduce file size.
c     Wrote Y4 on SCRACH34.DAT to reduce file size. 2/1/07
c     Removed MAMP and IREIG(J) from input file 2/2/07
c
c     Eliminated the search mode - use only Galerkin 9/21/07
c     Corrected an old (4/29/97) error in subroutine AMAT that
c     existed for non-constant density profile in the water
c     layer and passed FKBL and FKBR to subroutine CRSCUP 9/23/07
c     Implemented automatic choice of CW and CB 9/25/07
c
c     Required profile input for range zero and allowed
c     slopes between all input range depth pairs 10/3/07
c
c     Eliminated multiple source depth array 10/10/07
c     Printed out the CW and CB 11/3/07
c
c     The main variables are as follows.
c
c        Input-
c
c          HB=Total thickness (m) of water and bottom sediment layers
c          FREQ=Frequency (Hz)
c          ZS=Source depth (m)
c
c
c          M=Number of contributing modes (max.=MX=400)
c          IGEOM=Source-geometry flag.
c                    IGEOM=0: Use a point source in cylindrical
c                             geometry.
c                    IGEOM=1: Use a line source in plane
c                             geometry.
c          NEWSBC=Flag for the source region boundary condition at
c                 range zero.
c                    NEWSBC=0: Use the correct boundary condition
c                              for perfect cylindrically symmetric  
c                              geometry that allows standing waves
c                              in the source region.
c                    NEWSBC=1: Use a radiation boundary condition
c                              at the source that is consistent 
c                              with the single scatter approximation
c                              in cylindrical geometry and is
c                              always used with plane geometry. 
c
c          NDEP=Number of receiver depths for pressure calculation 
c                       (max.=MD=400)
c          ZMIN=Min. receiver depth (m) for pressure calculation
c          ZINC=Inc. in receiver depth (m) for pressure calculation 
c
c          RMIN=Min. range (km) for pressure calculation 
c          RMAX=Max. range (km) for pressure calculation
c          RINC=Inc. in range (km) for pressure calculation 
c
c          N=Number of input range-depth pairs (RANG,DPTH) to follow.
c          IPRT=Print flag for details of calculation such as eigen-
c               values and coefficients of outgoing and backscattered
c               waves. If IPRT >0 the results are printed out
c               for the first and the last region and every IPRT
c               region in between.
c          IOUT=Flag which controls what sort of mode coupling
c                      is used:
c               IOUT= -1: Compute an approximation to the single
c                         scatter approximation, due to Porter 
c                         and Jensen (see subroutine PROPM).
c               IOUT= 0: Compute the fully elliptic two way 
c                        solution that matches both  pressure 
c                        and radial particle velocity.
c               IOUT= 1: Compute the one-way pressure 
c                        matching solution.
c               IOUT= 2: Compute the adiabatic approximation
c
c               ** All the remaining inputs are repeated for J=1,N **
c
c          RANG(J)=Range (km) at the start of the J th environmental
c               region from RANG(J) to RANG(J+1) determined by the
c               following profile and/or new linearly sloping region.
c               RANG(1)=0.0 is required. RANG(N+1) is automatically
c               assumed to be 10,000 km and is not input.
c          IRLIN(J)=Number of additional regions generated to rep-
c                resent a linearly sloping bottom between RANG(J) and 
c                and RANG(J+1) for J=1,N. IRLIN(N) is ignored since
c                the bathymetry is assumed to be flat beyond RANG(N). 
c          DPTH(J)=Depth (m) of water at RANG(J) 
c          NPW(J)=Number of points in the water sound speed,
c             attenuation and density profile in the region between
c             RANG(J) and RANG(J+1). NPW(1) >1 for the 
c             first environmental region and a profile must be input.
c             If NPW(J) <0, in subsequent regions, then no profile
c             follows and the previously input profile is reused with
c             interpolation or extrapolation to accommodate bathymetry.
c
c          NPB(J)=Number of points in the bottom sound speed,
c             attenuation and density profile in the region between
c             RANG(J) and RANG(J+1). If NPW(J) >1 then 
c             NPB(J)>1 and a bottom profile must be provided.
c             If NPW(J) <0 then NPB(J)<0, no profile follows,
c             and the previously input bottom profile is reused.
c             The previous bottom profile is translated up and 
c             down to follow the bathymetry and retain same sound
c             speed, attenuation and density at the top of the
c             sediment. The translation forces the interpolation
c             or extrapolation, of the bottom profiles, to occur
c             at HB.
c
c               * The following inputs are repeated for L=1,NWP(J) * 
c
c               DEPW(L,J)=Depth (m) of sound speed in the water with
c                           DEPW(1,J)=0.0 and DEPW(NWP(J),J)=DPTH(J)
c               SVPW(L,J)=Sound speed (m/sec) in the water 
c               DBPWLW(L,J)=Attenuation (dB/wave length) in the
c                           water.
c               RHOW(L,J)=Density (gm/cm**3) in the water 
c
c               * The following inputs are repeated for L=1,NWB(J) * 
c
c               DEPB(L,J)=Depth (m) of sound speed in the bottom
c                         with DEPB(1,J)=DEPW(NPW(J),J) and
c                         DEPB(NWB(J),J)=HB. Note that DEPB(L,J)
c                         is measured from the ocean surface.
c               SVPB(L,J)=Sound speed (m/sec) in the bottom.
c               DBPWLB(L,J)=Attenuation (db/wave length) in the
c                           bottom.
c               RHOB(L,J)=Density in the bottom (g/cm**3).
c
c        Intermediate-
c
c          CW= Background (average) sound speed in the water
c              saved in cw_bck on a region by region basis
c          CB= Background (average) sound speed in the bottom.
c              CB is used compute the index of refraction.
c              saved in cb_bck on a region by region basis
c          ROHW= Background density in the water layer
c          ROHB= Background density in the bottom sediment layer
c
c          NREG=Number of computational regions on return from
c               subroutine BATHY .
c
c               * The following are defined for JJ=1,NREG *
c
c           RANGE(JJ)=Ending ranges of locally flat computational 
c                regions generated in subroutine BATHY to represent
c                a sloping bottom. RANGE(NREG)=10,000 km.
c           DEPTH(JJ)=Ending depths (and starting depths) of locally 
c                flat computational regions generated in subroutine
c                BATHY. DEPTH(NREG)=DEPTH(NREG-1).
c           NEWSVP(JJ)=Character array equal to 'NEW' when a new
c                 environmental profile needs to be incorporated
c                 into the calculation and equal to 'OLD' when
c                 a previously input environmental profile should
c                 be reused.
c
c          BSOUT=Integer flag which is set to abs(IOUT) in INPUT. 
c                If BSOUT =0 the two-way calculation is done.
c                If BSOUT >0 it turns off calculation of backscatter
c                and yields an outgoing solution that is either the
c                the one-way, the single scatter approximation or the
c                adiabatic approximation depending on IOUT. 
c          LHW=The water depth HW is between DEPW(LHW) and 
c                DEPW(LHW+1). If LHW=NPW extrapolation is done.
c          LHB=The greatest bottom depth HB is between DEPB(LHB) 
c                and DEPB(LHB+1). If LHB=NPB extrapolation is done.
c          cwnsw=Intercepts in the representation of the complex  
c                wave number squared in the water. Stored on LUSVP.
c          cgwnsw=Gradients in the representation of the complex 
c                 wave number squared in the water. Stored on LUSVP.
c          CINTW=Intercepts in the representation of the complex  
c                index of refraction squared in the water.
c          CGRADW=Gradients in the representation of the complex 
c                 index of refraction squared in the water.
c          GRHOW=Logarithmic gradients of the density in the water
c                (1/m). Computed in subroutine INPUT and stored on
c                LUSVP.
c          cwnsb=Intercepts in the representation of the complex 
c                wave number squared in the bottom. Stored on LUSVP.
c          cgwnsb=Gradients in the representation of the complex  
c                 wave number squared in the bottom. Stored on LUSVP.
c          CGRADB=Gradients in the representation of the complex
c                 index of refraction squared in the bottom. 
c          CINTB=Intercepts in the representation of the complex 
c                index of refraction squared in the bottom
c          GRHOB=Logarithmic gradients of the density in the bottom
c                (1/m). computed in subroutine INPUT and stored on
c                LUSVP.
c          ROHWL=Density of water on left side of vertical interface
c                in subroutine CRSCUP that was used to find the basic
c                eigenvalues and eigenfunctions.
c          ROHWR=Density of water on right side of vertical interface
c                in subroutine CRSCUP that was used to find the basic
c                eigenvalues and eigenfunctions.
c          ROHBL=Density of bottom on left side of a vertical
c                interface in subroutine CRSCUP that was used to
c                find the basic eigenvalues and eigenfunctions.
c          ROHBR=Density of bottom on right side of a vertical
c                interface in subroutine CRSCUP that was used to
c                find the basic eigenvalues and eigenfunctions.
c          BEGVL=Basic set of eigenvalues on left
c          BEGVR=Basic set of eigenvalues on right
c          BMATL=Matrix of normalized eigenvectors corresponding
c                to BEGVL (coefficients for the depth functions)
c          BMATR=Matrix of normalized eigenvectors corresponding
c                to BEGVR (coefficients for the depth functions)
c          EGVL=Eigenvalue on the left
c          EGVR=Eigenvalue on the right
c          FNORM=Normalization factor
c          CLR=Cross coupling matrix (uses density profile on right)
c          CRL=Cross coupling matrix (uses density profile on left)
c          R=Propagator matrix   
c          T=Unitary transformation matrix, R(J)*T(J)=T(J+1)*U(J)
c          T3=Sub block 3 of T in last region
c          T4=Sub block 4 of T in last region
c          T3Y1=T3*Y1 in last region
c          U=Upper triangular matrix
c          Y1=Decoupled outgoing solution (block 1)
c          Y3,Y4=Decoupled ingoing solution (block 3,4)
c          C=Vector (2*M x 1)that generates solution to matrix 
c                boundary value problem.
c          C2,C1=Top and bottom half of C
c          BA=Wave vector (2*M x 1)
c          B=Coeficients of ingoing waves, top half of BA
c          A=Coeficients of outgoing waves, bottom half of BA
c          CDF=Depth functions at output depths
c          ZR=Reciever depths
c          ZZ=Real eigenvalues
c
c        Computed (in CPRES)-
c
c          P=Complex pressure
c          P0=Scattered complex pressure in region 1
c             if BSOUT =0
c          P1=Outgoing complex pressure if BSOUT=0
c          P2=Ingoing complex pressure if BSOUT=0
c          TL=Transmision loss
c          TL0=Scattered transmision loss in region 1
c              if BSOUT =0
c          TL1=Outgoing transmision loss if BSOUT=0
c          TL2=Ingoing transmision loss if BSOUT=0
c
c     The files used are as follows.
c
c        Input-
c
c          UNIT5=Input
c
c        Intermediate-
c
c          UNIT24=Block matrix UB or U4 in upper 
c                 triangular matrix U
c          UNIT23=Block matrix UC or U3 in upper 
c                 triangular matrix U
c          UNIT31=Decoupled outgoing solution Y1
c          UNIT33=Decoupled ingoing solutions Y3
c          UNIT34=Decoupled ingoing solutions Y4
c          UNIT44=Unitary transformation matrix T4
c          UNIT43=Unitary transformation matrix T3
c          UNIT42=Unitary transformation matrix T2
c          UNIT41=Unitary transformation matrix T1
c          UNIT8=Eigenvalues and normalization factors.
c                Also basic eigenvalues and normalized
c                eigenvectors       
c          LUSVP=30, Tables of complex wave number squared and 
c                density in the water and bottom. The linear and
c                logarithmic gradients of the wave number squared
c                and density, respectively, are also tabulated.
c                Written in subroutine INPUT and read in the main
c                program COUPLE. 
c
c        Output-
c
c          UNIT6= Log file
c          UNIT9= Printed output
c          UNIT7= Real and imaginary parts of complex pressure
c          UNIT10= Transmission loss
c          UNIT12= Outgoing transmission loss if BSOUT=0
c          UNIT13= Ingoing transmission loss if BSOUT=0
c          UNIT14= Outgoing complex pressure if BSOUT=0
c          UNIT15= Ingoing complex pressure if BSOUT=0
c          UNIT16= Scattered complex pressure in region 1
c                  if BSOUT=0
c          UNIT17= Scattered transmission loss in region 1
c                  BSOUT=0
c
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MX=400,MXT2=2*MX,MD=400,NRX=2002,MXSD2=MX*(MX+1)/2)
C
      DIMENSION IPS(MX),INTH(MX)
      DIMENSION RANG(NRX),DPTH(NRX),RANGE(NRX),DEPTH(NRX)
      DIMENSION IRLIN(NRX)
      DIMENSION DEPW(MX),SVPW(MX),DBPWLW(MX),RHOW(MX),GRHOW(MX)
      DIMENSION DEPB(MX),SVPB(MX),DBPWLB(MX),RHOB(MX),GRHOB(MX)
      DIMENSION DEPWL(MX),RHOWL(MX),GRHOWL(MX)
      DIMENSION DEPWR(MX),RHOWR(MX),GRHOWR(MX)
      DIMENSION DEPBL(MX),RHOBL(MX),GRHOBL(MX)
      DIMENSION DEPBR(MX),RHOBR(MX),GRHOBR(MX)
      DIMENSION cw_bck(NRX),cb_bck(NRX)
      DIMENSION RHOWBOT(NRX),RHOBTOP(NRX)
      DIMENSION ZR(MD),TL(MD),TL0(MD),TL1(MD),TL2(MD)
C
      REAL*8 ZZ(MX)
      REAL*8 RSPD(MX,MX),VMAT(MX,MX)
      REAL*8 FA(MXSD2),FB(MXSD2),FU(MX,MX),FV(MX,MX)
      REAL*8 FCS(MX),FZ1(MX),FZ2(MX),FD(MX),FE(MX)
      COMPLEX*16 CINTW(NRX),CGRADW(NRX)
      COMPLEX*16 CINTB(NRX),CGRADB(NRX)
      COMPLEX*16 cwnsw(NRX),cgwnsw(NRX)
      COMPLEX*16 cwnsb(NRX),cgwnsb(NRX)
      COMPLEX*16 FNORM(MX)
      COMPLEX*16 EGVL(MX),EGVR(MX),BEGVL(MX),BEGVR(MX)
      COMPLEX*16 CDSRT,BDF
      COMPLEX*16 CRL(MX,MX),CLR(MX,MX),BMATL(MX,MX),BMATR(MX,MX)
      COMPLEX*16 R(MXT2,MXT2),T(MXT2,MXT2),U(MXT2,MXT2)
      COMPLEX*16 WRK2(MXT2,MXT2),WRK1(MX,MX),WRK3(MX,MX)
      COMPLEX*16 Y1(MX,MX),Y3(MX,MX),Y4(MX,MX)
      COMPLEX*16 T3(MX,MX),T4(MX,MX),T3Y1(MX,MX)
      COMPLEX*16 UB(MX,MX),UC(MX,MX),UE(MX,MX)
      COMPLEX*16 C1(MX),C2(MX),A(MX),B(MX),WRK(MX)
      COMPLEX*16 C(MXT2),BA(MXT2),WRK12(MXT2)
      COMPLEX*16 CDF(MX,MD)
      COMPLEX*16 P(MD),P0(MD),P1(MD),P2(MD)
      COMPLEX*16 CI,ONE,FKH
      COMPLEX*16 CXP,ZERO
      COMPLEX*16 HK01,SUM,SUMA
      COMPLEX*16 ALPHA(MX), FAC(MX), RAT(MX),
     1           BETA(MX),ALPHAR(MX),BETAR(MX),FNORMR(MX),
     2           ALPHAL(MX),BETAL(MX),FNORML(MX)
C
      INTEGER BSOUT
      CHARACTER*3 NEWSVP(NRX)
      CHARACTER*10 TITLE(8)
      LOGICAL READP(NRX)
C
      COMMON /ACCUR/ CONST1
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
      COMMON /BLKCON/ CI,ONE,ZERO
      COMMON /FLAGS/ IGEOM,NEWSBC
C
      EQUIVALENCE (CRL(1,1),UB(1,1)),(CLR(1,1),UC(1,1))
      EQUIVALENCE (BMATL(1,1),UE(1,1))
      EQUIVALENCE (R(1,1),U(1,1))
C
C     START UNIX OR g95 TIMER
C
c      ISTARTS=time()
C
C     START LAHEY TIMER
C
      CALL TIMER(IHU)
      STARTS=.01*IHU
C
      PI=2.*ACOS(0.0)
      CI=DCMPLX(0.0,1.0)
      ONE=DCMPLX(1.0,0.0)
      ZERO=DCMPLX(0.0,0.0)
C
C     OPEN INPUT AND PRINTED OUTPUT FILES. Unit6 is the command
c     window when using g95 and need not be opened. It may be
c     opened to divert this output to couple.log
C
      OPEN(5,FILE='COUPLE.DAT',FORM='FORMATTED') 
c      OPEN(6,FILE='COUPLE.LOG',FORM='FORMATTED') 
      OPEN(9,FILE='COUPLE.PRT',FORM='FORMATTED') 
      LUSVP=30
      OPEN(LUSVP,FILE='LUSVP.DAT',FORM='FORMATTED')
C
      WRITE(6,10)
      WRITE(9,10)
   10 FORMAT(1X,/,'  ********  ',/,
     & ' THIS PC VERSION OF COUPLE IS RANGE DEPENDENT ',/,
     & ' HAS PLANE OR CYLINDRICAL GEOMETRY OPTION ',/,
     & ' AND STANDARD OR NEW INITIAL SOURCE BOUNDARY CONDITIONS ',/,
     & '  ********  ',/)
C
C     GET INPUT DATA 
C
      CALL INPUT(TITLE,FREQ,M,N,IPRT,IOUT,BSOUT,NDEP,
     1           RANGE,DEPTH,NREG,
     2           ZS,RMIN,RMAX,RINC,ZR,NDINC,MD,NRX,READP,
     3           RANG,DPTH,IRLIN,MX,
     4           DEPW,SVPW,DBPWLW,cgwnsw,cwnsw,RHOW,GRHOW,
     5           DEPB,SVPB,DBPWLB,cgwnsb,cwnsb,RHOB,GRHOB,
     6           NEWSVP,LUSVP)
C
C
c      istop=1
c      if(istop .eq. 1) stop
C
      H1=DEPTH(1)
C
C      RKM1=RANGE(1)
C
C     NEWSVP(1)='NEW' is always true so NEWSVP(1) is not referenced.
C
      REWIND LUSVP
      READ(LUSVP,*) NPW,NPB
      DO 11 KK=1,NPW
      READ(LUSVP,*) DEPW(KK),cwnsw(KK),cgwnsw(KK),RHOW(KK),GRHOW(KK)
   11 CONTINUE
      DO 12 KK=1,NPB
      READ(LUSVP,*) DEPB(KK),cwnsb(KK),cgwnsb(KK),RHOB(KK),GRHOB(KK)
   12 CONTINUE
c
c
      CALL BOUNDS(DEPW,cgwnsw,cwnsw,NPW,MX,H1,fminw,favgw,fmaxw,LHW)
      CALL BOTDEPS(H1,DEPB,cgwnsb,cwnsb,NPB,LHB,MX)
      CALL BOUNDS(DEPB,cgwnsb,cwnsb,NPB,MX,HB,fminb,favgb,fmaxb,LHB)
c
c  Use the average water wave number and the average bottom wave
c  number to define the sound speeds used in the background two-
c  layer problem.
c
      FKW=SQRT(favgw)
      CW=2.0*PI*FREQ/FKW
      FKB=SQRT(favgb)
      CB=2.0*PI*FREQ/FKB
c
      write(6,*) ' '
      write(6,*) '1, RANGE(1) (km)= ',1,RANGE(1)
      write(6,*) 'CW and CB (m/sec)= ',CW,CB
      write(6,*) ' '
      write(9,*) ' '
      write(9,*) '1, RANGE(1) (km)= ',1,RANGE(1)
      write(9,*) 'CW and CB (m/sec)= ',CW,CB
      write(9,*) ' '
c
c      FKWmin=SQRT(fminw)
c      FKWavg=SQRT(favgw)
c      FKWmax=SQRT(fmaxw)
c      CWmin=2.0*PI*FREQ/FKWmax
c      CWavg=2.0*PI*FREQ/FKWavg
c      CWmax=2.0*PI*FREQ/FKWmin
c      write(*,*) 'CWmin, CWavg, CWmax= ',CWmin, CWavg, CWmax
c
c      FKBmin=SQRT(fminb)
c      FKBavg=SQRT(favgb)
c      FKBmax=SQRT(fmaxb)
c      CBmin=2.0*PI*FREQ/FKBmax
c      CBavg=2.0*PI*FREQ/FKBavg
c      CBmax=2.0*PI*FREQ/FKBmin
c      write(*,*) 'CBmin, CBavg, CBmax= ',CBmin, CBavg, CBmax
C
      cw_bck(1)=CW
      cb_bck(1)=CB
      RHOWBOT(1)=RHOW(LHW)*DEXP(GRHOW(LHW)*(H1-DEPW(LHW)))
      ROHW=RHOWBOT(1)
      RHOBTOP(1)=RHOB(1)
      ROHB=RHOBTOP(1)
C
C     CHECK FOR MODE CUTOFF, AVOID IF NECESSARY BY INCREASING
C     THE WATER DEPTH BY 1/1000 OF A WAVE LENGTH, AND REDEFINE
C     THE maximum WAVE NUMBER IN THE WATER AND THE DEPTHS
C     AND INTERCEPTS IN THE BOTTOM
C
        CUTOFF=RCF(H1,1.0D0)
        IF(DABS(CUTOFF) .LT. 1.0D-12) THEN
        DEPINC=2.0D0*PI/(1000.0D0*FKW)
        WRITE(6,*) ' *** WARNING *** '
        WRITE(6,*) 'MODE CUTOFF AT RANGE=',RANGE(1),' KM'
        WRITE(6,*) 'DEPTH INCREASED BY ',DEPINC,' M'
        WRITE(9,*) ' '
        WRITE(9,*) ' *** WARNING *** '
        WRITE(9,*) 'MODE CUTOFF AT RANGE=',RANGE(1),' KM'
        WRITE(9,*) 'DEPTH INCREASED BY ',DEPINC,' M'
        DEPTH(1)=DEPTH(1)+DEPINC
        H1=DEPTH(1)
c
        CALL BOUNDS(DEPW,cgwnsw,cwnsw,NPW,MX,H1,fminw,favgw,fmaxw,LHW)
        CALL BOTDEPS(H1,DEPB,cgwnsb,cwnsb,NPB,LHB,MX)
        CALL BOUNDS(DEPB,cgwnsb,cwnsb,NPB,MX,HB,fminb,favgb,fmaxb,LHB)
c
        FKW=SQRT(favgw)
        CW=2.0*PI*FREQ/FKW
        FKB=SQRT(favgb)
        CB=2.0*PI*FREQ/FKB
c
        cw_bck(1)=CW
        cb_bck(1)=CB
c
        RHOWBOT(1)=RHOW(LHW)*DEXP(GRHOW(LHW)*(H1-DEPW(LHW)))
        ROHW=RHOWBOT(1)
        END IF
C
C
      NPWL=NPW
      LHWL=LHW
      DO 14 KK=1,NPWL
      DEPWL(KK)=DEPW(KK)
      RHOWL(KK)=RHOW(KK)
      GRHOWL(KK)=GRHOW(KK)
   14 CONTINUE
C
      NPBL=NPB
      LHBL=LHB
      DO 16 KK=1,NPBL
      DEPBL(KK)=DEPB(KK)
      RHOBL(KK)=RHOB(KK)
      GRHOBL(KK)=GRHOB(KK)
   16 CONTINUE
C
C      WRITE(9,*) 'NPBL,LHBL= ',NPBL,LHBL
C      WRITE(6,*) 'NPBL,LHBL= ',NPBL,LHBL
C
CCC
      FKWL=FKW
      FKBL=FKB
      ROHBL=ROHB
      ROHWL=ROHW
C
C  Define the complex index of refraction squared for use in GALRKN
C
         DO 25 KK=1,NPW
         CINTW(KK)=cwnsw(kk)/fkb**2
         CGRADW(KK)=cgwnsw(kk)/fkb**2
   25    CONTINUE
C
         DO 26 KK=1,NPB
         CINTB(KK)=cwnsb(kk)/fkb**2
         CGRADB(KK)=cgwnsb(kk)/fkb**2
   26    CONTINUE
c
c         WRITE(9,*) 'CW,CB= ',CW,CB
c         WRITE(9,*) 'NPW,LHW= ',NPW,LHW
c         DO 27 KK=1,NPW
c         WRITE(9,28) DEPW(KK),CINTW(KK)+CGRADW(KK)*DEPW(KK),CGRADW(KK)
c   27    CONTINUE
c         WRITE(9,*) 'NPB,LHB= ',NPB,LHB
c         DO 29 KK=1,NPB
c         WRITE(9,26) DEPB(KK),CINTB(KK)+CGRADB(KK)*DEPB(KK),CGRADB(KK)
c   28    FORMAT(1X,F10.4,2X,F10.4,2X,E12.4,2X,E12.4,2X,E12.4)
c   29    CONTINUE
C
C     OPEN SCRATCH AND OTHER FILES
C
      CALL FILE(BSOUT,MX,MXT2)
C
C     WRITE HEADERS ON THE OUTPUT FILES
C
      CALL HEADER(TITLE,FREQ,ZS,NDEP,NDINC,ZR,NREG,RANGE,DEPTH,MD,
     1            NRX,BSOUT)
C
C     FIND THE FIRST SET OF EIGENVALUES AND NORMALIZATION
C     FACTORS
C
      CALL GALRKN(H1,M,ZZ,EGVL,FNORM,
     1            BEGVL,BMATL,VMAT,MX,WRK,WRK1,RSPD,CLR,CRL,IPS,INTH,
     2            MXSD2,FA,FB,FU,FV,FCS,FZ1,FZ2,FD,FE,
     3            DEPW,CGRADW,CINTW,RHOW,GRHOW,NPW,LHW,
     4            DEPB,CGRADB,CINTB,RHOB,GRHOB,NPB,LHB)
C
C
c      WRITE(*,35) 1,RANGE(1),DEPTH(1)
C
      IF(IPRT .LE. 0) GO TO 60
C
      WRITE(9,30)
   30 FORMAT(/)
      WRITE(6,35) 1,RANGE(1),DEPTH(1)
      WRITE(9,35) 1,RANGE(1),DEPTH(1)
   35 FORMAT(1X,26HTHE EIGENVALUES IN REGION ,I5,11H ENDING AT ,F10.3,
     1       16H  KM WITH DEPTH ,F10.2,18H  M ARE AS FOLLOWS//)
      WRITE(9,37)
   37 FORMAT(2X,3HNO.,6X,10HEIGENVALUE,9X,20HNORMALIZATION FACTOR,
     1       8X,15HHORIZ. WAVE NO./)
C
      DO 50 I=1,M
      FKH=FKB*CDSRT(ONE-EGVL(I))
      IF(REAL(FKH) .LT. 0.0)   THEN
      WRITE(9,*) ' ERROR IN FIRST REGION :'
      WRITE(9,*) ' NEGATIVE WAVE NUMBER FOR MODE NO. ',I
      WRITE(9,*) ' WAVE NUMBER : ',FKH
      END IF
      WRITE(9,45) I,EGVL(I),FNORM(I),FKH
   45 FORMAT(1X,I3,1P,2(1X,2E11.3),1X,2E13.5)
   50 CONTINUE
C
C
   60 CONTINUE
C
C
      IF(BSOUT .NE. 0) GO TO 65
C
C     FOR THE DECOUPLING ALGORITHM ONLY
C
      WRITE(8) (BEGVL(I),I=1,M)
      WRITE(8) ((BMATL(II,I),II=1,M),I=1,M)
      WRITE(8) (EGVL(I),I=1,M)
      WRITE(8) (FNORM(I),I=1,M)
   65 CONTINUE
C
      RNG1=1000.*RANGE(1)
      DO 67 I=1,M
      C(I)=ZERO
      C(M+I)=ZERO
   67 CONTINUE
C
C     STORE THE VECTOR d determined by the source depth ZS
C     IN BOTTOM HALF OF C.
C
C     USE DEPTH FUNCTIONS CONSTRUCTED FROM THE BASIC
C     EIGENFUNCTIONS AND THE MATRIX EIGENVECTORS IN
C     THE ARRAY BMATL 
C
      DO 80 I=1,M
      ALPHAL(I)=FKB*CDSRT((FKW/FKB)**2+BEGVL(I)-ONE)
      BETAL(I)=FKB*CDSRT(BEGVL(I))
      FKH=FKB*CDSRT(ONE-EGVL(I))
      IF(IGEOM .EQ. 1) THEN
      FAC(I)=CDEXP(CI*RNG1*FKH)*CI/(2.0D0*FKH)
      ELSE
      FAC(I)=HK01(RNG1*FKH)*CI/4.0D0
      END IF
   80 CONTINUE
C
      IF(ZS .GT. HB) STOP 'SOURCE DEPTHS MUST BE .LE. HB'
C
C     FIND THE DENSITY AT ZS
C
      IF(ZS .LE. H1) THEN
      DO 82 L=1,LHW
      IF(ZS .GE. DEPW(L)) THEN
      RHOZSD=RHOW(L)*DEXP(GRHOW(L)*(ZS-DEPW(L)))
      END IF
   82 CONTINUE
      ELSE
      DO 86 L=1,LHB
      IF(ZS .GE. DEPB(L)) THEN
      RHOZSD=RHOB(L)*DEXP(GRHOB(L)*(ZS-DEPB(L)))
      END IF
   86 CONTINUE
      END IF
C
      DO 100 I=1,M
C
C     SUM OVER BASIC EIGENFUNCTIONS
C
      SUM=ZERO
      DO 87 II=1,M
      SUM=SUM+BMATL(II,I)*BDF(ZS,H1,ALPHAL(II),
     1                        BETAL(II),FNORM(II))
   87 CONTINUE
C
      C(M+I)=SUM*FAC(I)/RHOZSD
C
  100 CONTINUE
C
      IF(BSOUT .EQ. 0) GO TO 200
C
C     FIND THE OUTGOING COEFFICIENTS IN THE FIRST REGION (BSOUT>0)
C
      WRITE(9,105)
  105 FORMAT(//)
      WRITE(6,110) BSOUT
      WRITE(9,110) BSOUT
  110 FORMAT(1X,'BACKSCATTER WILL BE NEGLECTED AND ONLY AN',
     1       ' OUTGOING SOLUTION IS FOUND,','BSOUT= ',I5,/)   
C    
      IF (IOUT .EQ. 1) THEN  
      WRITE(6,112) IOUT
      WRITE(9,112) IOUT
  112 FORMAT(1X,'ONE-WAY SOLUTION IS USED: IOUT= ',I5,/)   
      ELSE IF (IOUT .EQ. 2) THEN
      WRITE(6,113) IOUT
      WRITE(9,113) IOUT
  113 FORMAT(1X,'ADIABATIC APPROXIMATION IS USED: IOUT= ',I5,///)
      ELSE IF (IOUT .EQ. -1) THEN
      WRITE(6,114) IOUT
      WRITE(9,114) IOUT
  114 FORMAT(1X,'SINGLE SCATTER APPROXIMATION IS USED: IOUT= ',I5,///)
      END IF
C
      DO 120 I=1,M
      A(I)=C(M+I)
      B(I)=ZERO
  120 CONTINUE
C
C      UNORMALIZED IN AND OUTGOING COEFICIENTS ARE PRINTED
C      BELOW AND IN ALL SUBSQUENT PRINTS.
C
c      WRITE(*,130) 1,RANGE(1),DEPTH(1)
C
      IF(IPRT .LE. 0) GO TO 150
C
      WRITE(6,130) 1,RANGE(1),DEPTH(1)
      WRITE(9,130) 1,RANGE(1),DEPTH(1)
  130 FORMAT(1X,27HTHE COEFFICIENTS IN REGION ,I5,11H ENDING AT ,F10.3,
     1       16H  KM WITH DEPTH ,F10.2,18H  M ARE AS FOLLOWS///
     2       10H  MODE NO.,11X,8HOUTGOING,12X,13HOUTGOING (DB)//)
C
      DO 140 I=1,M
C
      ATL=200.
      AABS=CDABS(A(I))
      IF(AABS .GT. 1.0D-10) ATL=-20.*DLOG10(AABS)
C
      WRITE(9,135) I,A(I),ATL
  135 FORMAT(1X,I5,5X,1X,2E12.4,5X,F8.2)
  140 CONTINUE
C
  150 CONTINUE
C
C     THIS INITIALIZES THE RANGE LOOP FOR THE PRESSURE CALCULATION
C     (BSOUT>0)
C
      RKM=RMIN
C
      IF(RANGE(1) .LT. RMIN) GO TO 200
C
      R2KM=RANGE(1)
      R1KM=RANGE(1)
C
      JBOT=0
C
C     COMPUTE THE DEPTH FUNCTIONS AND PRESSURE IN THE FIRST REGION
C     (BSOUT>0)
C
      CALL CDFUN(M,EGVL,FNORM,H1,ZR,NDEP,CDF,MX,MD,
     1           BEGVL,BMATL,ALPHA,FAC,RAT,BETA)
C
C
      CALL CPRES(M,RINC,RMAX,R1KM,R2KM,RKM,EGVL,A,B,
     1            NDEP,NDINC,CDF,JBOT,MX,MD,P,P0,P1,P2,
     2            TL,TL0,TL1,TL2,BSOUT)
C
C
  200 CONTINUE
C
C
C ******** START OUTGOING LOOP TO FIND EIGENVALUES,COMPUTE *********
C ******** PRESSURE (BSOUT>0) AND GENERATE R, T, U AND Y1   *********
C
      NREM1=NREG-1
      DO 1000 J=1,NREM1
      H1=DEPTH(J)
      H2=DEPTH(J+1)
C
C      RKMJP1=RANGE(J+1)
C
C     FIND EIGENVALUES AND NORMALIZATION FACTORS
C
       IF(NEWSVP(J+1) .EQ. 'NEW')   THEN
         READ(LUSVP,*) NPW,NPB
         DO 215 KK=1,NPW
         READ(LUSVP,*) DEPW(KK),cwnsw(KK),cgwnsw(KK),RHOW(KK),GRHOW(KK)
  215    CONTINUE
         DO 216 KK=1,NPB
         READ(LUSVP,*) DEPB(KK),cwnsb(KK),cgwnsb(KK),RHOB(KK),GRHOB(KK)
  216    CONTINUE
       END IF
      CALL BOUNDS(DEPW,cgwnsw,cwnsw,NPW,MX,H2,fminw,favgw,fmaxw,LHW)
      CALL BOTDEPS(H2,DEPB,cgwnsb,cwnsb,NPB,LHB,MX)
      CALL BOUNDS(DEPB,cgwnsb,cwnsb,NPB,MX,HB,fminb,favgb,fmaxb,LHB)
C
      FKW=SQRT(favgw)
      CW=2.0*PI*FREQ/FKW
      FKB=SQRT(favgb)
      CB=2.0*PI*FREQ/Fkb
C
      cw_bck(j+1)=CW
      cb_bck(j+1)=CB
c
      if(mod(j+1,iprt) .eq. 0) then
      write(6,*) ' '
      write(6,*) 'J+1, RANGE(J+1) (km)= ',J+1,RANGE(J+1)
      write(6,*) 'CW and CB (m/sec)= ',CW,CB
      write(6,*) ' '
      write(9,*) ' '
      write(9,*) 'J+1, RANGE(J+1) (km)= ',J+1,RANGE(J+1)
      write(9,*) 'CW and CB (m/sec)= ',CW,CB
      write(9,*) ' '
      end if
c
       RHOWBOT(J+1)=RHOW(LHW)*DEXP(GRHOW(LHW)*(H2-DEPW(LHW)))
       ROHW=RHOWBOT(J+1)
       RHOBTOP(J+1)=RHOB(1)
       ROHB=RHOBTOP(J+1)
C
C     CHECK FOR MODE CUTOFF, AVOID IF NECESSARY BY INCREASING
C     THE WATER DEPTH BY 1/1000 OF A WAVE LENGTH, AND REDEFINE
C     THE maximum WAVE NUMBER IN THE WATER AND THE DEPTHS
C     AND INTERCEPTS IN THE BOTTOM
C
        CUTOFF=RCF(H2,1.0D0)
        IF(DABS(CUTOFF) .LT. 1.0D-12) THEN
        DEPINC=2.0D0*PI/(1000.0D0*FKW)
        WRITE(6,*) ' *** WARNING *** '
        WRITE(6,*) 'MODE CUTOFF AT RANGE=',RANGE(J+1),' KM'
        WRITE(6,*) 'DEPTH INCREASED BY ',DEPINC,' M'
        WRITE(9,*) ' '
        WRITE(9,*) ' *** WARNING *** '
        WRITE(9,*) 'MODE CUTOFF AT RANGE=',RANGE(J+1),' KM'
        WRITE(9,*) 'DEPTH INCREASED BY ',DEPINC,' M'
        DEPTH(J+1)=DEPTH(J+1)+DEPINC
        H2=DEPTH(J+1)
        CALL BOUNDS(DEPW,cgwnsw,cwnsw,NPW,MX,H2,fminw,favgw,fmaxw,LHW)
        CALL BOTDEPS(H2,DEPB,cgwnsb,cwnsb,NPB,LHB,MX)
        CALL BOUNDS(DEPB,cgwnsb,cwnsb,NPB,MX,HB,fminb,favgb,fmaxb,LHB)
C
        FKW=SQRT(favgw)
        CW=2.0*PI*FREQ/FKW
        FKB=SQRT(favgb)
        CB=2.0*PI*FREQ/FKB
C
        cw_bck(j+1)=CW
        cb_bck(j+1)=CB
c
        RHOWBOT(J+1)=RHOW(LHW)*DEXP(GRHOW(LHW)*(H2-DEPW(LHW)))
        ROHW=RHOWBOT(J+1)
        END IF
C
C  Define the complex index of refraction squared for use in GALRKN
C
         DO 218 KK=1,NPW
         CINTW(KK)=cwnsw(kk)/fkb**2
         CGRADW(KK)=cgwnsw(kk)/fkb**2
  218    CONTINUE
C
         DO 220 KK=1,NPB
         CINTB(KK)=cwnsb(kk)/fkb**2
         CGRADB(KK)=cgwnsb(kk)/fkb**2
  220    CONTINUE
C
c         WRITE(9,*) 'CW,CB= ',CW,CB
c         WRITE(9,*) 'NPW,LHW= ',NPW,LHW
c         DO 225 KK=1,NPW
c         WRITE(9,227) DEPW(KK),CINTW(KK)+CGRADW(KK)*DEPW(KK),CGRADW(KK)
c  225    CONTINUE
c         WRITE(9,*) 'NPB,LHB= ',NPB,LHB
c         DO 229 KK=1,NPB
c         WRITE(9,227) DEPB(KK),CINTB(KK)+CGRADB(KK)*DEPB(KK),CGRADB(KK)
c  227    FORMAT(1X,F10.4,2X,F10.4,2X,E12.4,2X,E12.4,2X,E12.4)
c  229    CONTINUE
C
C
        CALL GALRKN(H2,M,ZZ,EGVR,FNORM,
     1             BEGVR,BMATR,VMAT,MX,WRK,WRK1,RSPD,CLR,CRL,IPS,INTH,
     2             MXSD2,FA,FB,FU,FV,FCS,FZ1,FZ2,FD,FE,
     3             DEPW,CGRADW,CINTW,RHOW,GRHOW,NPW,LHW,
     4             DEPB,CGRADB,CINTB,RHOB,GRHOB,NPB,LHB)
C
C
      NPWR=NPW
      LHWR=LHW
      DO 250 KK=1,NPWR
      DEPWR(KK)=DEPW(KK)
      RHOWR(KK)=RHOW(KK)
      GRHOWR(KK)=GRHOW(KK)
  250 CONTINUE
C
      NPBR=NPB
      LHBR=LHB
      DO 253 KK=1,NPBR
      DEPBR(KK)=DEPB(KK)
      RHOBR(KK)=RHOB(KK)
      GRHOBR(KK)=GRHOB(KK)
  253 CONTINUE
C
C
      FKWR=FKW
      FKBR=FKB
      ROHBR=ROHB
      ROHWR=ROHW
C
C
      DO 280 I=1,M
      FKH=FKB*CDSRT(ONE-EGVR(I))
      IF(REAL(FKH) .LT. 0.0)   THEN
      WRITE(9,*) ' ERROR IN REGION NO :',J+1
      WRITE(9,*) ' NEGATIVE WAVE NUMBER FOR MODE NO. ',I
      WRITE(9,*) ' WAVE NUMBER : ',FKH
      END IF
  280 CONTINUE
      JP1=J+1
C
c      WRITE(*,35) JP1,RANGE(JP1),DEPTH(JP1)
C
      IF(IPRT .LE. 0) GO TO 330
      IF(JP1 .EQ. NREG) GO TO 290
      IF(MOD(JP1,IPRT) .NE. 0) GO TO 330
  290 CONTINUE
C
      WRITE(9,30)
      WRITE(6,35) JP1,RANGE(JP1),DEPTH(JP1)
      WRITE(9,35) JP1,RANGE(JP1),DEPTH(JP1)
      WRITE(9,37)
C
      DO 300 I=1,M
      FKH=FKB*CDSRT(ONE-EGVR(I))
      WRITE(9,45) I,EGVR(I),FNORM(I),FKH
  300 CONTINUE
C
  330 CONTINUE
C
      IF(BSOUT .NE. 0) GO TO 350
      WRITE(8) (BEGVR(I),I=1,M)
      WRITE(8) ((BMATR(II,I),II=1,M),I=1,M)
      WRITE(8) (EGVR(I),I=1,M)
      WRITE(8) (FNORM(I),I=1,M)
  350 CONTINUE
C
C     COMPUTE THE CROSS COUPLING MATRICES, EVEN IN THE CASE
C     OF THE ADIABATIC APPROXIMATION WHERE THEY ARE NEEDED
C     TO CHECK FOR A POLARITY CHANGE OF THE MODES (SEE:
C     "COMPUTAIONAL ACOUSTICS," BY JENSEN, KUPERMAN, PORTER
C     AND SCHMIDT, P. 323 (AIP PRESS, 1994)
C
      CALL CRSCUP(H1,H2,M,BEGVL,BEGVR,CLR,CRL,MX,
     1            BMATL,BMATR,WRK1,ALPHAR,BETAR,FNORMR,
     2            ALPHAL,BETAL,FNORML,FKWL,FKWR,FKBL,FKBR,
     3            ROHWL,NPWL,DEPWL,RHOWL,GRHOWL,LHWL,
     4            ROHWR,NPWR,DEPWR,RHOWR,GRHOWR,LHWR,
     5            ROHBL,NPBL,DEPBL,RHOBL,GRHOBL,LHBL,
     6            ROHBR,NPBR,DEPBR,RHOBR,GRHOBR,LHBR)
C      
C
      IF(IPRT .EQ. 1) CALL PRTMAT(M,M,'CLR ',CLR,MX) 
      IF(IPRT .EQ. 1) CALL PRTMAT(M,M,'CRL ',CRL,MX) 
C
C      WRITE(9,*) 'ROHBL,ROHBR= ',ROHBL,ROHBR
C      WRITE(9,*) 'H1,H2= ',H1,H2
C      WRITE(*,*) 'H1,H2= ',H1,H2
C
      RNG1=RANGE(J)*1000.
      IF(J .GT. 1)  RNG1=RANGE(J-1)*1000.
      RNG2=RANGE(J)*1000.
      MT2=M*2
C
C     COMPUTE THE PROPAGATOR MATRIX R
C
      CALL PROPM(M,RNG1,RNG2,EGVL,EGVR,CLR,CRL,R,MX,MXT2,IOUT)
C
C
      IF(BSOUT .EQ. 0) GO TO 500
C
C
C     GENERATE THE OUTGOING COEFFICIENTS (BSOUT=1 OR 2)
C
C     COMPUTE THE 1-WAY OUTGOING SOLUTION IN THE CASE IOUT=1,
C     THE ADIABATIC APPROXIMATION IN THE CASE IOUT=2, OR THE
C     APPROXIMATE SINGLE SCATTER APPROXIMATION IN THE CASE IOUT=-1.
C     THESE CASES ARE DISTINGUISHED BY THE WAY THE PROPAGATOR
C     MATRIC R IS DEFINED IN PROPM.
C
      DO 380 II=1,M
      SUMA=ZERO
      DO 375 I=1,M
      SUMA=SUMA+R(II+M,I+M)*A(I)
  375 CONTINUE
      WRK(II)=SUMA
  380 CONTINUE
C
      DO 385 I=1,M
      A(I)=WRK(I)
      B(I)=ZERO
  385 CONTINUE  
C
c      WRITE(*,130) JP1,RANGE(JP1),DEPTH(JP1)
C
      IF(IPRT .LE. 0) GO TO 400
      IF(JP1 .EQ. NREG) GO TO 390
      IF(MOD(JP1,IPRT) .NE. 0) GO TO 400
  390 CONTINUE
C
      WRITE(9,105)
      WRITE(6,130) JP1,RANGE(JP1),DEPTH(JP1)
      WRITE(9,130) JP1,RANGE(JP1),DEPTH(JP1)
C
      DO 395 I=1,M
C
      ATL=200.
      AABS=CDABS(A(I))
      IF(AABS .GT. 1.0D-10) ATL=-20.*DLOG10(AABS)
      WRITE(9,135) I,A(I),ATL
  395 CONTINUE
C
  400 CONTINUE
C
      IF((RANGE(J+1) .LT. RMIN) .OR. (RANGE(J) .GE. RMAX)) GO TO 900
      IF(RKM .GT. RANGE(J+1)) GO TO 900
C
      R1KM=RANGE(J)
      R2KM=RANGE(J+1)
      JBOT=0
C
C     CALCULATE DEPTH FUNCTION AND PRESSURE (BSOUT=1 OR 2)
C
C
      CALL CDFUN(M,EGVR,FNORM,H2,ZR,NDEP,CDF,MX,MD,
     1           BEGVR,BMATR,ALPHA,FAC,RAT,BETA)
C
C
      CALL CPRES(M,RINC,RMAX,R1KM,R2KM,RKM,EGVR,A,B,
     1            NDEP,NDINC,CDF,JBOT,MX,MD,P,P0,P1,P2,
     2            TL,TL0,TL1,TL2,BSOUT)
C
C
      GO TO 900
C 
  500 CONTINUE  
C
      IF(J .GT. 1) GO TO 580
C
      DO 550 I=1,MT2
      DO 540 II=1,MT2
      WRK2(II,I)=R(II,I)
  540 CONTINUE
  550 CONTINUE
C
      GO TO 600
C
  580 CONTINUE
C
C     STORE THE PRODUCT R*T IN WRK2 (WRK2=R*T)
C
      CALL MMULT(R,T,WRK2,MT2,MXT2) 
C
C
  600 CONTINUE
C
C     DO MODIFIED GRAM-SCHMIDT DECOMPOSITION OF (WRK2=T*U) 
C
C
      CALL MGS(WRK2,T,U,MT2,MXT2)
C
C     THE MATRICES T AND  U ARE PARTITIONED AS FOLLOWS:
C
C         T4 T3              UB  UC
C     T =         AND   U = 
C         T2 T1              0   UE
C
C     DEFINE UE
C
      DO 650 I=1,M
      DO 640 II=1,M
      UE(II,I)=U(M+II,M+I)
  640 CONTINUE
  650 CONTINUE
C
C     STORE UB, UC AND T
C
      WRITE(24) ((U(II,I),II=1,I),I=1,M)
      WRITE(23) ((U(II,M+I),II=1,M),I=1,M)
C
      WRITE(44) ((T(II,I),II=1,M),I=1,M)
      WRITE(43) ((T(II,M+I),II=1,M),I=1,M)
      WRITE(42) ((T(II+M,I),II=1,M),I=1,M)
      WRITE(41) ((T(II+M,I+M),II=1,M),I=1,M)
C
      IF(J .GT. 1) GO TO 820
C
      DO 800 I=1,M
      DO 700 II=1,M
C
      Y1(II,I)=UE(II,I)
C
  700 CONTINUE
  800 CONTINUE
C
      GO TO 850
C
  820 CONTINUE
C
C     ADVANCE Y1 (Y1=UE*Y1)
C
      CALL MMULT(UE,Y1,WRK1,M,MX) 
C
C
      DO 840 I=1,M
      DO 830 II=1,I
      Y1(II,I)=WRK1(II,I)
  830 CONTINUE
  840 CONTINUE
C
  850 CONTINUE
C
      WRITE(31) ((Y1(II,I),II=1,I),I=1,M)
C
  900 CONTINUE
C
      DO 920 I=1,M
      BEGVL(I)=BEGVR(I)
      DO 910 II=1,M
      BMATL(II,I)=BMATR(II,I)
  910 CONTINUE
  920 CONTINUE
C
      DO 950 I=1,M
      EGVL(I)=EGVR(I)
  950 CONTINUE
C
      NPWL=NPWR
      LHWL=LHWR
      DO 970 KK=1,NPWL
      DEPWL(KK)=DEPWR(KK)
      RHOWL(KK)=RHOWR(KK)
      GRHOWL(KK)=GRHOWR(KK)
  970 CONTINUE
C
      NPBL=NPBR
      LHBL=LHBR
      DO 975 KK=1,NPBL
      DEPBL(KK)=DEPBR(KK)
      RHOBL(KK)=RHOBR(KK)
      GRHOBL(KK)=GRHOBR(KK)
  975 CONTINUE
CCC
      FKWL=FKWR
      FKBL=FKBR
      ROHBL=ROHBR
      ROHWL=ROHWR
CCC
C     END OF OUTGOING LOOP
C
 1000 CONTINUE
C
C ****** IF BSOUT=1 OR 2,THE OUTGOING SOLUTION HAS BEEN  ******
C ******        OBTAINED THE CALCULATION IN OVER         ******
C
      IF(BSOUT .NE. 0) GO TO 8000
C
C     THE MATRIX Y1, IN THE LAST REGION, IS NOT READ FROM UNIT31
C     IN THE LOOP TO GENERATE Y3 & Y4, BUT IT IS MULTIPLLIED BY T3
C     AND SAVED IN T3Y1. BACKSPACE OVER END OF FILE AND LAST Y1.
C
      ENDFILE 31
      BACKSPACE 31
      BACKSPACE 31
C
C     END FILES AND BACKSPACE OVER END OF FILE
C
      ENDFILE 24
      BACKSPACE 24
      ENDFILE 23
      BACKSPACE 23
C
      ENDFILE 44
      REWIND 44
      ENDFILE 43
      REWIND 43
      ENDFILE 42
      REWIND 42
      ENDFILE 41
      REWIND 41
C
      ENDFILE 8
      REWIND 8
C
C     INITIALIZE Y4 & Y3 IN PREPARATION FOR FINDING THE FUNDAMENTAL
C     MATRIX SOLUTION. SAVE T4 AND T3Y1 FOR LATER USE.
C
      DO 1200 I=1,M
      DO 1100 II=1,M
      T4(II,I)=T(II,I)
      T3(II,I)=T(II,M+I)
      Y4(II,I)=ZERO
      IF(II .EQ. I) Y4(II,I)=ONE
      Y3(II,I)=ZERO
 1100 CONTINUE
 1200 CONTINUE
C
C     T3Y1 = T3 * Y1
C
      CALL MMULT(T3,Y1,T3Y1,M,MX)
C
C
      WRITE(33) ((Y3(II,I),II=1,M),I=1,M)  
      WRITE(34) ((Y4(II,I),II=1,M),I=1,M)  
C
C********* START IN GOING LOOP TO GENERATE Y3 & Y4 **************
C
C
      DO 2000 JJ=1,NREM1
      J=NREM1-JJ+1
C
      BACKSPACE 24
      READ(24) ((UB(II,I),II=1,I),I=1,M)
      BACKSPACE 24
C
      BACKSPACE 23
      READ(23) ((UC(II,I),II=1,M),I=1,M)
      BACKSPACE 23
C
      IF(J .GT. 1) GO TO 1500
C
      DO 1400 I=1,M
      DO 1300 II=1,M
      WRK1(II,I)=Y3(II,I)-UC(II,I)
      WRK3(II,I)=Y4(II,I)
 1300 CONTINUE
 1400 CONTINUE
      GO TO 1900
C
 1500 CONTINUE
C
      BACKSPACE 31
      READ(31) ((Y1(II,I),II=1,I),I=1,M)
      BACKSPACE 31
C
C
      CALL MMULT(UC,Y1,WRK1,M,MX)
C
C
      DO 1800 I=1,M
      DO 1700 II=1,M
      WRK1(II,I)=Y3(II,I)-WRK1(II,I)
      WRK3(II,I)=Y4(II,I)
 1700 CONTINUE
 1800 CONTINUE
C
C
 1900 CONTINUE
C
C
      CALL BACKSUB(UB,Y3,WRK1,M,MX)
      CALL BACKSUB(UB,Y4,WRK3,M,MX)
C
C
      IF(J .GT. 1) THEN
      WRITE(33) ((Y3(II,I),II=1,M),I=1,M)
      WRITE(34) ((Y4(II,I),II=1,M),I=1,M)
      END IF
C
C     END OF INGOING LOOP
C
 2000 CONTINUE
C
C     END FILES AND BACKSPACE OVER END OF FILES
C
      ENDFILE 33
      BACKSPACE 33
C
      ENDFILE 34
      BACKSPACE 34
C 
      REWIND 31
C
C     SET UP THE LINEAR SYSSTEM OF EQUATIONS FOR THE VECTOR C
C
      READ(8) (BEGVL(I),I=1,M)
      READ(8) ((BMATL(II,I),II=1,M),I=1,M)
      READ(8) (EGVL(I),I=1,M)
      READ(8) (FNORM(I),I=1,M)
C
      H1=DEPTH(1)
C
      ROHW=RHOWBOT(1)
      ROHB=RHOBTOP(1)
      CW=cw_bck(1)
      FKW=2.0*PI*FREQ/CW
      CB=cb_bck(1)
      FKB=2.0*PI*FREQ/CB
C
      RNG1=RANGE(1)*1000.
C
      IF(NEWSBC .GT. 0)   THEN
      WRITE(6,2100)
      WRITE(9,2100)
 2100 FORMAT(1X,/,
     & ' *** NEW INITIAL SOURCE BOUNDARY CONDITIONS ARE COMPUTED ***',/)
      END IF
C
      DO 2500 II=1,M
C
      FKH=FKB*CDSRT(ONE-EGVL(II))
      IF(-DIMAG(2.*FKH*RNG1) .LT. -88.) GO TO 2150
      IF(IGEOM .EQ. 1) THEN
      CXP=CDEXP(2.*CI*FKH*RNG1)
      ELSE
      CXP=-CI*CDEXP(2.*CI*FKH*RNG1)
      END IF
      GO TO 2200
 2150 CONTINUE
      CXP=ZERO
 2200 CONTINUE
      IF(NEWSBC .NE. 0) CXP=ZERO
      DO 2400 I=1,M
      T(II,I)=T4(II,I)
      T(II,M+I)=T3Y1(II,I)
      T(M+II,I)=-CXP*Y4(II,I)
      T(M+II,M+I)=-CXP*Y3(II,I)
      IF(I .EQ. II) T(M+II,M+I)=ONE+T(M+II,M+I)
 2400 CONTINUE
C
 2500 CONTINUE
C
C
C ********** SOLVE MATRIX EQUATION FOR C *****************
C
      CALL CLINEQS(MT2,MXT2,1,1,T,C,WRK2,WRK12,BA,IPS,IER)
C     
C
      IFLAG=3
C
C     ABORT NOT COMMENTED OUT
C
      IF(IER .NE. 0) CALL ABORTC(IFLAG)
C
      WRITE(9,105)
      IF(IER .EQ. 1)   THEN
      WRITE(6,1261) 
      WRITE(9,1261) 
      ELSE
      WRITE(6,1260) IER
      WRITE(9,1260) IER
      END IF
 1260 FORMAT(1X,4HIER=,I5///)
 1261 FORMAT(1X,///,
     & ' ***************************************************** ',
     & //,' WARNING :  FAILURE TO IMPROVE BY IMPRUV IN CLINEQS ',/, 
     &    '            AS ITERATION DID NOT CONVERGE IN SOLVE ',//,
     & ' *********************************************** ',///)
C
c      WRITE(*,3690) 1,RANGE(1),DEPTH(1)
      WRITE(6,3690) 1,RANGE(1),DEPTH(1)
      WRITE(9,3690) 1,RANGE(1),DEPTH(1)
 3690 FORMAT(1X,27HTHE COEFFICIENTS IN REGION ,I5,11H ENDING AT ,F10.3,
     1       16H  KM WITH DEPTH ,F10.2,18H  M ARE AS FOLLOWS///
     2       10H  MODE NO.,11X,8HOUTGOING,14X,11HBACKSCATTER,
     3       12X,13HOUTGOING (DB),4X,16HBACKSCATTER (DB),3X,
     4       10HTOTAL (DB)//)
C
C     GENERATE THE COEFFICIENTS IN THE FIRST REGION
C
      DO 4000 II=1,M
      A(II)=C(M+II)
      C1(II)=C(M+II)
      C2(II)=C(II)
 4000 CONTINUE
C
C     B=Y3*C1+Y4+C2
C
      CALL CVMULT(M,Y3,C1,MX,WRK) 
      CALL CVMULT(M,Y4,C2,MX,WRK) 
C
      DO 4025 II=1,M
      B(II)=C1(II)+C2(II)
 4025 CONTINUE
C
      DO 4100 I=1,M
C
      ATL=200.
      BTL=200.
      TTL=200.
      AABS=CDABS(A(I))
      BABS=CDABS(B(I))
      TABS=CDABS(A(I)+B(I))
      IF(AABS .GT. 1.0D-10) ATL=-20.*DLOG10(AABS)
      IF(BABS .GT. 1.0D-10) BTL=-20.*DLOG10(BABS)
      IF(TABS .GT. 1.0D-10) TTL=-20.*DLOG10(TABS)
C
      WRITE(9,4050) I,A(I),B(I),ATL,BTL,TTL
 4050 FORMAT(1X,I5,5X,2(1X,2E12.4),5X,F8.2,2(10X,F8.2))
 4100 CONTINUE
C
C     THIS INITIALIZES THE RANGE LOOP FOR THE PRESSURE CALCULATION
C     (BSOUT .EQ. 0)
C
      RKM=RMIN
C
      IF(RANGE(1) .LT. RMIN) GO TO 4200
C
      R2KM=RANGE(1)
      R1KM=RANGE(1)
C
      JBOT=0
C
C     COMPUTE THE DEPTH FUNCTIONS AND PRESSURE IN THE FIRST REGION
C     (BSOUT .EQ. 0)
C
      CALL CDFUN(M,EGVL,FNORM,H1,ZR,NDEP,CDF,MX,MD,
     1           BEGVL,BMATL,ALPHA,FAC,RAT,BETA)
C
C
      CALL CPRES(M,RINC,RMAX,R1KM,R2KM,RKM,EGVL,A,B,
     1            NDEP,NDINC,CDF,JBOT,MX,MD,P,P0,P1,P2,
     2            TL,TL0,TL1,TL2,BSOUT)
C
C
 4200 CONTINUE
C
C ************** LOOP FOR PRESSURE CALCULATION *********************
C **************         (BSOUT .EQ. 0)         *********************
C
C
      DO 7000 J=1,NREM1
C
      READ(8) (BEGVR(I),I=1,M)
      READ(8) ((BMATR(II,I),II=1,M),I=1,M)
      READ(8) (EGVR(I),I=1,M)
      READ(8) (FNORM(I),I=1,M)
C
      H1=DEPTH(J)
      H2=DEPTH(J+1)
C
      ROHW=RHOWBOT(J+1)
      ROHB=RHOBTOP(J+1)
      CW=cw_bck(j+1)
      FKW=2.0*PI*FREQ/CW
      CB=cb_bck(j+1)
      FKB=2.0*PI*FREQ/CB
C
      READ(44) ((T(II,I),II=1,M),I=1,M)
      READ(43) ((T(II,M+I),II=1,M),I=1,M)
      READ(42) ((T(II+M,I),II=1,M),I=1,M)
      READ(41) ((T(II+M,I+M),II=1,M),I=1,M)
C
 6500 CONTINUE
C 
      READ(31) ((Y1(II,I),II=1,I),I=1,M)
C
      BACKSPACE 33
      READ(33) ((Y3(II,I),II=1,M),I=1,M)
      BACKSPACE 33
C
      BACKSPACE 34
      READ(34) ((Y4(II,I),II=1,M),I=1,M)
      BACKSPACE 34
C
C     GENERATE COEFFICIENTS A AND B
C
      CALL GENABC(T,Y1,Y3,Y4,R,C,BA,A,B,WRK12,WRK2,M,MT2,MX,MXT2)
C
C
      JP1=J+1
c      WRITE(*,3690) JP1,RANGE(JP1),DEPTH(JP1)
      IF(IPRT .LE. 0) GO TO 6850
      IF(JP1 .EQ. NREG) GO TO 6780
      IF(MOD(JP1,IPRT) .NE. 0) GO TO 6850
 6780 CONTINUE
C
      WRITE(9,105)
      WRITE(6,3690) JP1,RANGE(JP1),DEPTH(JP1)
      WRITE(9,3690) JP1,RANGE(JP1),DEPTH(JP1)
C
      DO 6800 I=1,M
C
      ATL=200.
      BTL=200.
      TTL=200.
      AABS=CDABS(A(I))
      BABS=CDABS(B(I))
      TABS=CDABS(A(I)+B(I))
      IF(AABS .GT. 1.0D-10) ATL=-20.*DLOG10(AABS)
      IF(BABS .GT. 1.0D-10) BTL=-20.*DLOG10(BABS)
      IF(TABS .GT. 1.0D-10) TTL=-20.*DLOG10(TABS)
C
      WRITE(9,4050) I,A(I),B(I),ATL,BTL,TTL
 6800 CONTINUE
      WRITE(9,105)
C
 6850 CONTINUE
C
      IF(J .LT. NREM1) GO TO 6880
C
      BMAX=0.
      DO 6860 I=1,M
      BMAX=DMAX1(BMAX,CDABS(B(I)))
      B(I)=ZERO
 6860 CONTINUE
C
      WRITE (6,105)
      WRITE(6,6870) BMAX
c      WRITE(*,6870) BMAX
      WRITE(9,6870) BMAX
 6870 FORMAT(1X,'THE MAXIMUM MAGNITUDE OF THE BACKSCATTER COEFFICIENTS',
     1 ' IN THE LAST REGION WAS',/,E12.4,25H BEFORE BEING SET TO ZERO/)
C
 6880 CONTINUE
C
      IF((RANGE(J+1) .LT. RMIN) .OR. (RANGE(J) .GE. RMAX)) GO TO 6900
      IF(RKM .GT. RANGE(J+1)) GO TO 6900
C
      R1KM=RANGE(J)
      R2KM=RANGE(J+1)
      JBOT=0
C
C     CALCULATE DEPTH FUNCTION AND PRESSURE (BSOUT .EQ. 0)
C
      CALL CDFUN(M,EGVR,FNORM,H2,ZR,NDEP,CDF,MX,MD,
     1           BEGVR,BMATR,ALPHA,FAC,RAT,BETA)
C
C
      CALL CPRES(M,RINC,RMAX,R1KM,R2KM,RKM,EGVR,A,B,
     1            NDEP,NDINC,CDF,JBOT,MX,MD,P,P0,P1,P2,
     2            TL,TL0,TL1,TL2,BSOUT)
C
C
 6900 CONTINUE
C
      DO 6920 I=1,M
      BEGVL(I)=BEGVR(I)
      DO 6910 II=1,M
      BMATL(II,I)=BMATR(II,I)
 6910 CONTINUE
 6920 CONTINUE
C
      DO 6950 I=1,M
      EGVL(I)=EGVR(I)
 6950 CONTINUE
C
C     END OF PRESSURE CALCULATION LOOP BSOUT=0
C
 7000 CONTINUE
C
C     END OF PRESSURE CALCULATION BSOUT>0
C
 8000 CONTINUE
C
      ENDFILE 7
      ENDFILE 10
      IF(BSOUT .EQ. 0) THEN
      ENDFILE 12
      ENDFILE 13
      ENDFILE 14
      ENDFILE 15
      ENDFILE 16
      ENDFILE 17
      END IF
C
      REWIND 7
      REWIND 10
      REWIND LUSVP
      IF(BSOUT .EQ. 0) THEN
      REWIND 12
      REWIND 13
      REWIND 14
      REWIND 15
      REWIND 16
      REWIND 17
      END IF
C
C     PRINT TRANSMISSION LOSS
C
      WRITE(6,*) 'CALL PRTTL'
      WRITE(9,*) 'CALL PRTTL'
      CALL PRTTL(TITLE,FREQ,ZS,ZR,TL)
C
C     CLOSE THE SCRATCH FILES
C
      CLOSE(24,STATUS='DELETE')
      CLOSE(23,STATUS='DELETE')
      CLOSE(34,STATUS='DELETE')
      CLOSE(33,STATUS='DELETE')
      CLOSE(31,STATUS='DELETE')
      CLOSE(44,STATUS='DELETE')
      CLOSE(43,STATUS='DELETE')
      CLOSE(42,STATUS='DELETE')
      CLOSE(41,STATUS='DELETE')
      CLOSE(8,STATUS='DELETE')
      CLOSE(18,STATUS='DELETE')
      CLOSE(LUSVP,STATUS='DELETE')
C 
C     FIND ELAPSED TIME IN UNIX or g95
C
c      IENDS=time()
c      ELAPSE=IENDS-ISTARTS
C 
C     FIND ELAPSED TIME IN LAHEY
C
      CALL TIMER(IHU)
      ENDS=.01*IHU
      ELAPSE=ENDS-STARTS
C
c      WRITE(*,8100) ELAPSE
      WRITE(6,8100) ELAPSE
      WRITE(9,8100) ELAPSE
 8100 FORMAT(1X,'TOTAL ELAPSED TIME= ',F10.3,'  SECONDS')    
C
 9000 CONTINUE
C
      ENDFILE 9
c
c      ENDFILE 6
C
      STOP
      END
