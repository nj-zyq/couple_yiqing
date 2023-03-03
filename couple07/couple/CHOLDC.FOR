      SUBROUTINE CHOLDC(A,G,M,MX)
C
C     CHOLESKY DECOMPOSITON OF A REAL SYMETRIC POSITIVE DEFINITE
C     MATRIX A IN THE FORM:
C                                    T
C                         A =  G   G
C
C     WHERE G IS LOWER TRIANGULAR. TAKEN FROM "FUNDEMENTALS OF 
C     MATRIX COMPUTATIONS" BY D. S. WATKINS, PG. 21 (JOHN WILEY,
C     NEW YORK, 1991)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(MX,MX),G(MX,MX)
      REAL*8 SUM,ZERO
C
      ZERO=0.0D0
C
      DO 100 J=1,M
      JM1=J-1
      JP1=J+1
      SUM=0.0D0
      IF(J .EQ. 1) GO TO 30
      DO 20 K=1,JM1
      SUM=SUM+G(J,K)**2
   20 CONTINUE
   30 CONTINUE
      IF(A(J,J) .LE. SUM) THEN
      WRITE(6,*) '*** NOT POSITIVE DEFINITE IN CHOLDC ***'
      WRITE(9,*) '*** NOT POSITIVE DEFINITE IN CHOLDC ***'
      STOP 'FAILURE: NOT POSITIVE DEFINITE IN CHOLDC'
      END IF
      G(J,J)=DSQRT(A(J,J)-SUM)
C
      IF(J .EQ. M) GO TO 100
C
      DO 60 I=JP1,M
      SUM=0.0D0
      IF(J .EQ. 1) GO TO 50
      DO 40 K=1,JM1
      SUM=SUM+G(I,K)*G(J,K)
   40 CONTINUE
   50 CONTINUE
      G(I,J)=(A(I,J)-SUM)/G(J,J)
      G(J,I)=ZERO
   60 CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
