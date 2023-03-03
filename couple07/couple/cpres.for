      SUBROUTINE CPRES(M,RINC,RMAX,R1KM,R2KM,RKM,EGVL,
     1                  A,B,NDEP,NDINC,CDF,JBOT,MX,MD,P,P0,P1,P2,
     2                  TL,TL0,TL1,TL2,BSOUT)
C
C     THIS SUBROUTINE COMPUTES THE COMPLEX PRESSURE,INTENSITY AND 
C     TRANSMISSION LOSS AT EACH RANGE RKM IN THE INTERVAL R1KM
C     .LT. RKM .LE. R2KM AT EACH OF THE NDEP OUTPUT DEPTHS ZR. THE 
C     PRESSURE (REAL,AIMAG) AND TL ARE WRITTEN ON UNIT7 AND UNIT10,
C     RESPECTIVELY. IF BSOUT .EQ. 0, THE OUTGOING, INGOING AND
C     SCATTERED COMPONENTS OF THE TL ARE WRITTEN ON UNIT12, UNIT13
C     AND UNIT17 RESPECTIVELY. IF BSOUT.EQ. 0, THE OUTGOING,INGOING
C     AND SCATTERED COMPONENTS OF THE COMPLEX PRESSURE ARE WRITTEN
C     ON UNIT14 UNIT15 AND UNIT16 RESPECTIVELY.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TL(MD),TL0(MD),TL1(MD),TL2(MD)
      INTEGER BSOUT
C
      REAL*8 PI
C
      COMPLEX*16 P(MD),P0(MD),P1(MD),P2(MD)
      COMPLEX*16 EGVL(MX)
      COMPLEX*16 A(MX),B(MX),CDF(MX,MD)
      COMPLEX*16 CHDH,CJDJ,CRH,CRJ
      COMPLEX*16 CDJ0,HK02,CDSRT,FKH
      COMPLEX*16 CARG,CARG1,CARGP,CARGM
      COMPLEX*16 CI,ZERO,ONE,EGV,SCALE
      COMPLEX*16 RCOEF,RCOEF0,RCOEF1,RCOEF2
      COMPLEX*16 CMODE,CMODE0,CMODE1,CMODE2
C
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
      COMMON/FLAGS/IGEOM,NEWSBC
C
      PI=2.*DACOS(DBLE(0.0))
      ZERO=DCMPLX(0.0,0.0)
      ONE=DCMPLX(1.0,0.0)
      CI=DCMPLX(0.0,1.0)
C
C     FIND NORMALIZATION FOR CYLINDRICAL GEOMETRY
C
      SCALE=DCMPLX(4.0D0*PI,0.0D0)
C
      IF(IGEOM.EQ.1) THEN
C
C     THE NORMALIZATION IN PLANE GEOMETRY DEPENDS ON FKW
C     DIVIDE BY (I/4)*HK01(FKW*R) AT R=1 M, IGNORING THE I.
C
C      SCALE=4.0D0/DCMPLX(BESSJ0(FKW),BESSY0(FKW))
C
C     USE THE NEAR FIELD FORM OF THE HANKEL FUNCTION.
C
C      SCALE=DCMPLX(-2.0D0*PI/DLOG(FKW),0.0D0)
C
C     USE THE FAR FIELD FORM OF THE HANKEL FUNCTION.
C
      SCALE=DCMPLX(DSQRT(8.0D0*PI*FKW),0.0D0)
C
      END IF
C
      RNG1=R1KM*1000.
C
  100 CONTINUE
C
      RNG=RKM*1000.
C
      DO 150 J=1,NDEP
      P(J)=ZERO
      IF(BSOUT .NE. 0) GO TO 150
      P0(J)=ZERO
      P1(J)=ZERO
      P2(J)=ZERO
  150 CONTINUE
C
      DO 300 I=1,M
C
      EGV=EGVL(I)
C
C     ALLOW LONG FIRST REGION
C
      IF((RNG .LT. RNG1) .AND. (CDABS(A(I)) .LT. 1.0D-30)) THEN
      CRH=ZERO
      ELSE
      CRH=CHDH(EGV,RNG1,DMAX1(RNG,1.0D-6))
      END IF
C
      IF(BSOUT .NE. 0)   THEN
      RCOEF=A(I)*CRH
      ELSE
C
      IF(CDABS(B(I)) .GT. 1.0D-30) THEN
      CRJ=CJDJ(EGV,RNG1,DMAX1(RNG,1.0D-6))
      RCOEF=A(I)*CRH+B(I)*CRJ
      ELSE
      RCOEF=A(I)*CRH
      END IF
C
      IF((RNG .LE. RNG1) .AND. (IGEOM .NE. 1)) THEN
C
C     IN CYLINDRICAL GEOMTERY THE SCATTERED FIELD IN THE
C     FIRST REGION IS IN TERMS OF THE BESSEL FUNCTION J0
C     DIVIDED BY THE HANKEL FUNCTIONS OF ORDER ZERO TYPE
C     TWO EVALUATED AT RNG1
C
      FKH=FKB*CDSRT(ONE-EGV)
      CARG=FKH*RNG
      CARG1=FKH*RNG1
C
      IF(CDABS(CARG) .LT. 2.0D0) THEN
C
C     THE POWER SERIES IS USED FOR J0 WHEN CARG IS SMALL
C     AVOID OVER. FLOW OF HK02(CARG1)
C
      IF(DIMAG(CARG1) .GT. 675.0D0) THEN
      RCOEF0=ZERO
      ELSE
      RCOEF0=B(I)*(2.0D0*CDJ0(CARG))/HK02(CARG1)
      END IF
C
      ELSE
C
C     OTHERWISE J0 IS WRITTEN IN TERMS OF HANKEL FUNCTIONS
C     AND THE LARGE ARGUMENT APPROXIMATIONS ARE USED
C
      CARGP=CI*FKH*(RNG1+RNG)
      CARGM=CI*FKH*(RNG1-RNG)
      RCOEF0=B(I)*(CDEXP(CARGP)+CDEXP(CARGM))*DSQRT(RNG1/RNG)
C
      END IF
C
      END IF
C
      RCOEF1=A(I)*CRH
      RCOEF2=B(I)*CRJ
C
      IF((RNG .LE. RNG1) .AND. (IGEOM .EQ. 1)) THEN
C
C     IN PLANE GEOMTERY THE SCATTERED FIELD IN THE
C     FIRST REGION IS JUST THE INGOING PLANE WAVE
C
      RCOEF0=RCOEF2
      END IF
C
      END IF
C
      DO 200 J=1,NDEP
C
C     THIS IS THE CRUX OF THE WHOLE CALCULATION.
C
      IF(BSOUT .NE. 0)   THEN
      CMODE=RCOEF*CDF(I,J)
      P(J)=P(J)+CMODE
      ELSE
      CMODE=RCOEF*CDF(I,J)
      CMODE1=RCOEF1*CDF(I,J)
      CMODE2=RCOEF2*CDF(I,J)
      P(J)=P(J)+CMODE
      P1(J)=P1(J)+CMODE1
      P2(J)=P2(J)+CMODE2
C
      IF(RNG .LE. RNG1) THEN
      CMODE0=RCOEF0*CDF(I,J)
      P0(J)=P0(J)+CMODE0
      END IF
C
      END IF
C
  200 CONTINUE
C
  300 CONTINUE
C
C
      IF(BSOUT .NE. 0) THEN
C
      DO 400 J=1,NDEP
      P(J)=P(J)*SCALE
      XI=CDABS(P(J))
      IF( XI .LE. 1.0E-10)   THEN
      TL(J)=200.
      ELSE IF( XI .GE. 1.0E10 )   THEN
      TL(J)=-200.
      ELSE
      TL(J)=-20.*DLOG10(XI)
      END IF
  400 CONTINUE
C
      WRITE(7) SNGL(RKM),
     1         (SNGL(DREAL(P(J))),SNGL(DIMAG(P(J))),J=1,NDEP)
      WRITE(10,450) RKM,(TL(J),J=1,NDEP,NDINC)
  450 FORMAT(F13.7,20F13.7)
C
      ELSE
C
      DO 500 J=1,NDEP
C
      IF(RNG .LE. RNG1) THEN
      P0(J)=P0(J)*SCALE
      XI0=CDABS(P0(J))
      IF( XI0 .LE. 1.0E-10)   THEN
      TL0(J)=200.
      ELSE IF( XI0 .GE. 1.0E10 )   THEN
      TL0(J)=-200.
      ELSE
      TL0(J)=-20.*DLOG10(XI0)
      END IF
      END IF
C
      P(J)=P(J)*SCALE
      P1(J)=P1(J)*SCALE
      P2(J)=P2(J)*SCALE
      XI=CDABS(P(J))
      XI1=CDABS(P1(J))
      XI2=CDABS(P2(J))
C
      IF( XI .LE. 1.0E-10)   THEN
      TL(J)=200.
      ELSE IF( XI .GE. 1.0E10 )   THEN
      TL(J)=-200.
      ELSE
      TL(J)=-20.*DLOG10(XI)
      END IF
C
      IF( XI1 .LE. 1.0E-10)   THEN
      TL1(J)=200.
      ELSE IF( XI1 .GE. 1.0E10 )   THEN
      TL1(J)=-200.
      ELSE
      TL1(J)=-20.*DLOG10(XI1)
      END IF
C
      IF( XI2 .LE. 1.0E-10)   THEN
      TL2(J)=200.
      ELSE IF( XI2 .GE. 1.0E10 )   THEN
      TL2(J)=-200.
      ELSE
      TL2(J)=-20.*DLOG10(XI2)
      END IF
C
  500 CONTINUE
C
      WRITE(7) SNGL(RKM),
     1         (SNGL(DREAL(P(J))),SNGL(DIMAG(P(J))),J=1,NDEP)
      WRITE(14) SNGL(RKM),
     1         (SNGL(DREAL(P1(J))),SNGL(DIMAG(P1(J))),J=1,NDEP)
      WRITE(15) SNGL(RKM),
     1         (SNGL(DREAL(P2(J))),SNGL(DIMAG(P2(J))),J=1,NDEP)
C     
      WRITE(10,450) RKM,(TL(J),J=1,NDEP,NDINC)
      WRITE(12,450) RKM,(TL1(J),J=1,NDEP,NDINC)
      WRITE(13,450) RKM,(TL2(J),J=1,NDEP,NDINC)
C
      IF(RNG .LE. RNG1) THEN
      WRITE(17,450) RKM,(TL0(J),J=1,NDEP,NDINC)
      WRITE(16) SNGL(RKM),
     1         (SNGL(DREAL(P0(J))),SNGL(DIMAG(P0(J))),J=1,NDEP)
      END IF
C
      END IF
C
      RKM=RKM+RINC
      IF((RKM .GT.R2KM) .OR. (RKM .GT.RMAX)) GO TO 900
C
      GO TO 100
C
  900 CONTINUE
C
      RETURN
      END