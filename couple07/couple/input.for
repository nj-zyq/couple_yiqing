      SUBROUTINE INPUT(TITLE,FREQ,M,N,IPRT,IOUT,BSOUT,NDEP,
     1                 RANGE,DEPTH,NREG,
     2                 ZS,RMIN,RMAX,RINC,ZR,NDINC,MD,NRX,READP,
     3                 RANG,DPTH,IRLIN,MX,  
     4                 DEPW,SVPW,DBPWLW,cgwnsw,cwnsw,RHOW,GRHOW,
     5                 DEPB,SVPB,DBPWLB,cgwnsb,cwnsb,RHOB,GRHOB,
     6                 NEWSVP,LUSVP)
C
C     THIS SUBROUTINE READS THE INPUTS THAT CONTROL THE STEP-
C     WISE COUPLED MODE CALCULATION, SUCH AS THE FREQUENCY (FREQ),
C     THE SOURCE DEPTH (ZS), THE NUMBER OF MODES (M), THE NUMBER
C     OF INPUT RANGE-DEPTH PAIRS (N), AND THE TOTAL THICKNESS OF
C     THE WAVEGUIDE (HB).
C
C     IT CALLS THE SUBROUTINE RWPROF TO READ THE ENVIRONMENTAL
C     PROFILES (DEPTH, SOUND SPEED, ATTENUATION AND DENSITY),
C     GENERATE THE WAVE NUMBER SQUARED TABLES AND WRITE THE
C     TABLES ON DISK (UNIT=LUSVP).
C
C     SUBROUTINE INPUT CALLS SUBROUTINE BATHY TO INTERPOLATE
C     THE BATHYMETRY OVER LINEARLY SLOPING BOTTOMS BASED ON
C     THE INPUT RANGE-DEPTH PAIRS AND IRLIN. SUBROUTINE BATHY
C     ALSO GENRATES THE CHARACTER VARIABLE NEWSVP THE CONTROL
C     WHEN A NEW ENVIRONMENTAL PROFILE IS INCORPORATED INTO
C     THE CALCULATION.
C
C     SUBROUTINE INPUT ALSO GENERATES THE RECEIVER DEPTHS
C     AND RANGES, AND SETS THE FLAGS IOUT AND BSOUT THAT
C     CONTROL THE TYPE OF MODE COUPLING, WHICH IS USED.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER BSOUT
      DIMENSION RANGE(NRX),DEPTH(NRX)
      DIMENSION ZR(MD)
      DIMENSION DEPW(MX),SVPW(MX),DBPWLW(MX),RHOW(MX),GRHOW(MX)
      DIMENSION DEPB(MX),SVPB(MX),DBPWLB(MX),RHOB(MX),GRHOB(MX)
      DIMENSION RANG(NRX),DPTH(NRX),IRLIN(NRX)
C
      COMPLEX*16 cgwnsw(MX),cwnsw(MX)
      COMPLEX*16 cgwnsb(MX),cwnsb(MX)
C
      CHARACTER*10 TITLE(8)
      CHARACTER*3 NEWSVP(NRX)
      LOGICAL READP(NRX)
C
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
      COMMON/FLAGS/IGEOM,NEWSBC
C
      NREG=0
C
C     READ INPUT 
C
      WRITE(9,50)
   50 FORMAT(1H1)
      READ(5,100) (TITLE(I),I=1,8)
  100 FORMAT(8A10)
      WRITE(6,110) (TITLE(I),I=1,8)
      WRITE(9,110) (TITLE(I),I=1,8)
  110 FORMAT(1X,8A10/)
C
      READ(5,*) HB,FREQ,ZS
      WRITE(6,120) HB,FREQ,ZS
      WRITE(9,120) HB,FREQ,ZS
  120 FORMAT(1X,'HB=',F10.2,2X,'FREQ=',F10.2,2X,'ZS=',F10.2,/)
C
      WRITE(9,125)
  125 FORMAT(/)
C
      READ(5,*) M,IGEOM,NEWSBC
      IF(IGEOM .EQ. 0)   THEN
      WRITE(6,155)
      WRITE(9,155)
      ELSE
      WRITE(6,160)
      WRITE(9,160)
      IGEOM=1
      NEWSBC=1
      END IF
  155 FORMAT(1X,//,' ********  CYLINDRICAL GEOMETRY PROBLEM *******',/)   
  160 FORMAT(1X,//,' ********  PLANE GEOMETRY PROBLEM ******* ',/)   
C
      IF(NEWSBC .EQ. 0)   THEN
      WRITE(6,165)
      WRITE(9,165)
      ELSE
      WRITE(6,170)
      WRITE(9,170)
      NEWSBC=1
      END IF
  165 FORMAT(1X,/,'  ********* ',/,
     & '  OLD INITIAL SOURCE BOUNDARY CONDITIONS PROBLEM',
     & /,'  ********* ',/)
  170 FORMAT(1X,/,' ********',
     & '  NEW INITIAL SOURCE BOUNDARY CONDITIONS PROBLEM',
     & '  ******** ',/)
C
      WRITE(9,125)
      WRITE(6,180) M
      WRITE(9,180) M
  180 FORMAT(1X,2HM=,I5,/)
C
C     READ RECEIVER SPECIFICATIONS
C
      WRITE (6,125)
      READ(5,*) NDEP,ZMIN,ZINC
      WRITE(6,210) NDEP,ZMIN,ZINC
      WRITE(9,210) NDEP,ZMIN,ZINC
  210 FORMAT(1X,'NDEP= ',I5,5X,'ZMIN=',F10.2,5X,'ZINC= ',F10.2)
C
      IF(NDEP .GT. MD) NDEP=MD
      NDINC= 1 
      DO 240 J=1,NDEP
      ZR(J)=ZMIN+(J-1)*ZINC
  240 CONTINUE
C
      WRITE (6,125)
      READ(5,*) RMIN,RMAX,RINC
      WRITE(6,250) RMIN,RMAX,RINC
      WRITE(9,250) RMIN,RMAX,RINC
  250 FORMAT(1X,'RMIN= ',F10.3,5X,'RMAX= ',F10.3,5X,'RINC= ',F10.3)
C
      READ(5,*) N,IPRT,IOUT
C
C     SET IOUT AND BSOUT
C
      IF(IOUT .EQ. 0) THEN
C
C     TWO-WAY SOLUTION
C
      BSOUT=0
      ELSE IF(IOUT .EQ. 1) THEN
C
C     ONE-WAY SOLUTION
C
      BSOUT=1
      ELSE IF(IOUT .EQ. 2) THEN
C
C     ADIABATIC APPROXIMATION
C
      BSOUT=2
      ELSE
C
C     SINGLE SCATTER APPROXIMATION
C
      IOUT=-1
      BSOUT=1
      END IF
C
      WRITE (6,125)
      WRITE(6,270) N,IPRT,IOUT,BSOUT
      WRITE(9,270) N,IPRT,IOUT,BSOUT
  270 FORMAT(1X,'N=',I5,3X,'IPRT=',I5,3X,'IOUT=',I5,3x,'BSOUT=',I5)
C
       IF(N .GT. NRX) THEN
       NRXX=NRX
       WRITE(9,280) N,NRXX
  280  FORMAT(1X,//,'  TOO MANY INPUT REGIONS  ',
     & /,' NUMBER OF REGIONS  : ',I5,
     & /,' MAX ALLOWED NUMBER : ',I5,
     & /,' EXECUTION IS TERMINATED ')
       STOP
       END IF 
C     WRITE(9,290) N
  290 FORMAT(1X,/,' ******',/,
     & ' THE BATHYMETRY AND SOUND VELOCITY PROFILES',
     & ' FOR THE ',I5,' INPUT REGION(S) ARE:',/,' ******',/)
C
      NRTOT=1
C
      DO 500 J=1,N
C
      WRITE(6,310) J
      WRITE(9,310) J
  310 FORMAT(1X,///,'   REGION NO. ',I4,/)
C
      READ(5,*) RANG(J),IRLIN(J),DPTH(J),NPW,NPB
C
      IF( (J .EQ. 1) .AND. (RANG(J) .GT. 0.0E0) ) THEN
      RANG(1)=0.0E0
      WRITE(9,320)
  320 FORMAT(1X,/,' WARNING : RANGE FOR FIRST PROFILE HAS BEEN',
     & ' RESET TO ZERO KM ')
      END IF
C
       IF( (NPW .GT. MX) .OR. (NPB .GT. MX) ) THEN
       MXP=MX
       WRITE(9,330) J,NPW,NPB,MXP
  330  FORMAT(1X,//,'  TOO MANY POINTS IN PROFILE NO. ',I5,
     & /,' NUMBER OF WATER PROFILE POINTS : ',I5,
     & /,' NUMBER OF BOTTOM PROFILE POINTS : ',I5,
     & /,' MAX ALLOWED NUMBER : ',I5,
     & /,' EXECUTION IS TERMINATED ')
       STOP
       END IF 
C
C     SET READP(J) AND CHECK FOR INCONSISTENT INPUTS
C
      IF((NPW .GE. 2) .AND. (NPB .GE. 2)) THEN
C
C     PROFILES WILL BE READ
C
      READP(J)=.TRUE.
C
      ELSE IF( (NPW .LT. 0) .AND. (NPB .LT. 0) .AND. (J .GT. 1) ) THEN
C
C     PROFILES WILL BE REUSED 
C
      READP(J)=.FALSE.
      ELSE
C
      WRITE(6,*) 'INCONSISTENT NPW AND NPB IN REGION= ',J
      WRITE(9,*) 'INCONSISTENT NPW AND NPB IN REGION= ',J
      STOP 'INCONSISTENT NPW AND NPB'
C
      END IF
C
      IF(J.LT.N)   THEN
      NRTOT=NRTOT+MAX(1,IRLIN(J))
      IF(NRTOT+N-J.GT.NRX)   THEN
      WRITE(9,*) ' TOO MANY REGIONS. REVISE INPUT RUN STREAM '
      STOP
      END IF
      END IF
C
      WRITE(6,490) RANG(J),DPTH(J),IRLIN(J)
      WRITE(9,490) RANG(J),DPTH(J),IRLIN(J)
  490 FORMAT(1X,/,' ****** ',/,
     & ' STARTING RANGE AND DEPTH, AND THE NUMBER OF',
     & ' SUBDIVISIONS BEORE NEXT RANGE: ',/,
     & 1X,F10.2,3X,F10.2,5X,I5,/,' ****** ',/)
C     
      IF(READP(J)) THEN
C
C     READ PROFILES, GENERATE INDEX OF REFRACTION TABLES AND
C     WRITE TO LUSVP.
C
      CALL RWPROF(FREQ,NPW,NPB,
     1                 DEPW,SVPW,DBPWLW,cgwnsw,cwnsw,RHOW,GRHOW,
     2                 DEPB,SVPB,DBPWLB,cgwnsb,cwnsb,RHOB,GRHOB,
     3                 LUSVP,MX)
C
      END IF
C
C
      IF(J.GT.1)   THEN
       IF(RANG(J).LE.RANG(J-1))   THEN
       WRITE(9,*) ' ERROR : REVISE INPUT RANGE ',J
       STOP
       END IF
      END IF
C
      IF(DPTH(J).EQ.0.0)   THEN
      DPTH(J)=1.0E-6
      WRITE(9,*) ' WARNING : '
      WRITE(9,*)' WATER DEPTH IS INCREASED TO 1.E-6 M IN SECTOR ',J
      END IF
C     
  500 CONTINUE
C
      ENDFILE LUSVP
C
      WRITE(9,125)
C
C     GENERATE REPRESENTATION OF A SLOPING BOTTOM
C
      CALL BATHY(N,RANG,DPTH,IRLIN,RANGE,DEPTH,NREG,NRX,
     &           NEWSVP,READP)
C
C      PRINT RANGE,DEPTH PAIRS
C
      WRITE(6,510)
      WRITE(9,510)
  510 FORMAT(1X,'THE RANGE DEPTH PAIRS ARE AS FOLLOWS:'/)
C
      NREG=MAX0(N,NREG)
C
C     The profile at range zero is always read ,
C     by default. NEWSVP(1) is not used in COUPLE 
C
      DO 600 J=1,NREG
      WRITE(6,550) J,RANGE(J),DEPTH(J),NEWSVP(J)
      WRITE(9,550) J,RANGE(J),DEPTH(J),NEWSVP(J)
  550 FORMAT(1X,I5,5X,F10.3,5X,F10.2,5X,A3)
C
  600 CONTINUE
C
C
      RETURN
      END