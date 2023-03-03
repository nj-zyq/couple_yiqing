      SUBROUTINE BATHY(N,RANG,DPTH,IRLIN,RANGE,DEPTH,NREG,NRX,
     &                  NEWSVP,READP)
C
C     THIS SUBROUTINE GENERATES A STEPWISE REPRESENTATION OF A
C     SLOPING BOTTOM. THE LOGIC DEPENDS ON IRLIN (IRLIN(N) IS
C     NOT USED). RANG,DPTH are the input range, depth pairs.
C     The pairs RANGE, DEPTH are generated to be consistent with
C     the COUPLE logic: RANGE(J) is the end range of the Jth
C     interval (RANG(1) is zero but RANGE(1) is positive).
C      
C
C     IRLIN(J) .LE. 1 - The depth in the interval RANG(J) to
C                       RANG(J+1) is DPTH(J) and changes to
C                       DPTH(J+1) at RANG(J+1).
C
C     IRLIN(J) .GT. 1 - The interval between RANG(J) and RANG(J+1)
C                       with length DRNG is divided into IRLIN(J)+1
C                       subintervals. The first and last are DRNG/2
C                       long with depth DPTH(J) and DPTH(J+1)
C                       respectively. The remaining IRLIN(J)-1
C                       subintervals are DRNG long. THE MID-
C                       POINTS OF THE TOPS OF THESE STEPS
C                       LIE ON THE LINE BETWEEN (DPTH(J),RANG(J))
C                       AND (DPTH(J+1),RANG(J+1)).
C
C     ADDED logical READP TO ALLOW REUSE OF PROFILES R. B. EVANS, 
C     SAIC, NEW LONDON, CT ON 4/20/98.
C
C     This version assumes that RANG(1) is zero.
C     (R. B. Evans, N. Stonington, CT, 10/10/07).
C     
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RANG(NRX),DPTH(NRX),IRLIN(NRX)
      DIMENSION RANGE(NRX),DEPTH(NRX)
      CHARACTER*3 NEWSVP(NRX)
      LOGICAL READP(NRX)
C
      DO 50 NR=1,NRX
      NEWSVP(NR)='OLD'
   50 CONTINUE
C
      JJ=0
C
C     NEWSVP determines when to start using a new profile.
C     NEWSVP(1) is always 'NEW' but it is not referenced in
C     COUPLE because the first profile is always read,
C     by default.
C
      DO 150 J=1,N-1
C
      IF(READP(J)) NEWSVP(JJ+1)='NEW'
C
      IF(IRLIN(J) .GT. 1 )   THEN
C
C     GENERATE LINEAR REPRESENTATION OF THE SLOPE
C
      DRNG=(RANG(J+1)-RANG(J))/IRLIN(J)
      DDEP=(DPTH(J+1)-DPTH(J))/IRLIN(J)
C
      DO 100 I=1,IRLIN(J)
      JJ=JJ+1
      RANGE(JJ)=RANG(J)+(FLOAT(I-1)+0.5)*DRNG
      DEPTH(JJ)=DPTH(J)+FLOAT(I-1)*DDEP
  100 CONTINUE
      JJ=JJ+1
      RANGE(JJ)=RANG(J+1)
      DEPTH(JJ)=DPTH(J+1)
C
      ELSE
C
      JJ=JJ+1
      RANGE(JJ)=RANG(J+1)
      DEPTH(JJ)=DPTH(J)
C
      END IF
C
  150 CONTINUE
C
      JJ=JJ+1
      IF(READP(N)) NEWSVP(JJ)='NEW'
      DEPTH(JJ)=DPTH(N)
      RANGE(JJ)=1.0E4
C
      NREG=JJ
C
      RETURN
      END
