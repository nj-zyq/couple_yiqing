      SUBROUTINE CLINEQS(N,NDIM,M,MDIM,A,B,CLU,X,R,IPS,IER) 
C     
C     THIS SUBROUTINE SOLVES THE COMPLEX MATRIX EQUATION A*X=B USING
C     COMPLEX*16 VERSIONS OF THE REAL NCAR SUBROUTINES IN LINEQSV.
C     A IS A NXN COMPLEX MATRIX AND B IS NXM COMPLEX MATRIX. THE RESULT
C     IS RETURNED IN B.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPS(NDIM)
C     
      REAL*8 DIGITS,XMIN,XMAX,XMAG
C
      COMPLEX*16 B(NDIM,MDIM),X(NDIM),R(NDIM)
      COMPLEX*16 A(NDIM,NDIM),CLU(NDIM,NDIM),CDET
C
C     SCALE THE MATRIX A AND B
C
      XMIN=DBLE(.1E+38)
      XMAX=DBLE(0.0)
C
      DO 20 I=1,N
      DO 20 II=1,N
C
      XMIN=DMIN1(XMIN,CDABS(A(II,I)))
      XMAX=DMAX1(XMAX,CDABS(A(II,I)))
C
   20 CONTINUE
C
      XMAG=(XMIN+XMAX)/2.
      IFLAG=2
      IF(XMAG .EQ. 0.0) CALL ABORTC(IFLAG)
C
      DO 40 I=1,N
      DO 40 II=1,N
      A(II,I)=A(II,I)/XMAG
   40 CONTINUE
C
      DO 50 I=1,M
      DO 50 II=1,N
      B(II,I)=B(II,I)/XMAG
   50 CONTINUE
C
C 
      CALL DECOMP(N,NDIM,A,CLU,CDET,IPS,IER)
C
C
      IF(CDABS(CDET) .EQ. 0.0) CALL ABORTC(IFLAG)
C
      DO 500 I=1,M
C
C
      CALL SOLVE(N,NDIM,CLU,B(1,I),IPS,X)
C
C
      CALL IMPRUV(N,NDIM,A,CLU,B(1,I),IPS,X,R,DIGITS,IER) 
C
C
      DO 100 II=1,N
      B(II,I)=X(II)
  100 CONTINUE
C
  500 CONTINUE
C
      RETURN
      END
C LINEQSV    FROM NSSL                                     08/15/79    
      SUBROUTINE DECOMP (N,NDIM,A,CLU,D,IPS,IER)    
C    
C PACKAGE LINEQSV    
C    
C    
C LATEST REVISION        JANUARY 1979    
C    
C PURPOSE                LINEQSV IS A PACKAGE OF THREE    
C                        SUBROUTINES--DECOMP, SOLVE AND IMPRUV--WHICH   
C                        MAY BE USED TO SOLVE A SYSTEM OF LINEAR    
C                        EQUATIONS BY GAUSSIAN ELIMINATION WITH PARTIAL 
C                        PIVOTING AND (OPTIONALLY) WITH ITERATIVE    
C                        IMPROVEMENT.    
C    
C ACCESS CARDS           *FORTRAN,S=ULIB,N=LINEQSV    
C    
C    
C USAGE                  TO SOLVE A SYSTEM OF EQUATIONS, AX = B, THE    
C                        USER FIRST CALLS DECOMP BY    
C                            CALL DECOMP (N,NDIM,A,CLU,D,IPS,IER)    
C                        WHICH DECOMPOSES A INTO THE PRODUCT OF A UNIT  
C                        LOWER TRIANGULAR MATRIX AND AN UPPER TRIANGULAR
C                        MATRIX.  AFTER CHECKING D TO INSURE THAT THE   
C                        DETERMINANT OF A IS NOT ZERO, THE USER THEN    
C                        FOLLOWS THIS CALL WITH A CALL TO SOLVE BY    
C                            CALL SOLVE (N,NDIM,CLU,B,IPS,X)    
C                        WHICH USES THE DECOMPOSITION TO FIND THE    
C                        SOLUTION X.  IF DESIRED, THE USER MAY THEN TRY 
C                        TO IMPROVE THE ACCURACY OF THE SOLUTION BY    
C                        CALLING IMPRUV BY    
C                            CALL IMPRUV (N,NDIM,A,CLU,B,IPS,X,R,DIGITS, 
C                                         IER)    
C                        IF THE USER WANTS TO SOLVE SEVERAL LINEAR    
C                        SYSTEMS OF EQUATIONS WITH THE SAME COEFFICIENT 
C                        MATRIX A, ONLY ONE CALL TO DECOMP IS NECESSARY.
C                        THIS IS THEN FOLLOWED BY A CALL TO SOLVE FOR   
C                        EACH DIFFERENT B VECTOR.    
C                        . EACH OF THE ABOVE ROUTINES IS WRITTEN UP IN  
C                          DETAIL BELOW.    
C    
C ENTRY POINTS           DECOMP, SOLVE, IMPRUV    
C    
C SPACE REQUIRED         723 (OCTAL) = 467 (DECIMAL) LOCATIONS    
C    
C SUBROUTINE DECOMP (N,NDIM,A,CLU,D,IPS,IER)    
C    
C    
C DIMENSION OF           A(NDIM,N),CLU(NDIM,N),IPS(N)    
C ARGUMENTS    
C    
C PURPOSE                DECOMPOSE MATRIX INTO THE PRODUCT OF UPPER AND 
C                        UNIT LOWER TRIANGULAR MATRICES USING GAUSSIAN  
C                        ELIMINATION WITH PARTIAL PIVOTING.    
C    
C USAGE                  CALL DECOMP (N,NDIM,A,CLU,D,IPS,IER)    
C    
C ARGUMENTS    
C    
C ON INPUT               N    
C                          ORDER OF MATRIX A.    
C    
C                        NDIM    
C                          DECLARED FIRST (ROW) DIMENSIONS OF A AND CLU  
C                          IN THE CALLING PROGRAM.    
C    
C                        A    
C                          MATRIX TO BE DECOMPOSED.    
C    
C ON OUTPUT              CLU    
C                          LOWER AND UPPER TRIANGULAR MATRICES OF A,    
C                          WHERE    
C                              CLU(I,J), FOR I .LE. J,    
C                          IS THE UPPER TRIANGULAR MATRIX WITH ROW    
C                          INTERCHANGES COMPLETED, AND    
C                              CLU(I,J), FOR I .GT. J,    
C                          IS THE NEGATIVE OF LOWER TRIANGULAR MATRIX   
C                          EXCLUDING DIAGONAL ELEMENTS, BUT WITH ROW    
C                          INTERCHANGES PARTIALLY COMPLETED.    
C    
C                          IF MATRIX A NEED NOT BE SAVED, THEN THE    
C                          MATRIX A MAY BE USED AS THE DECOMPOSITION    
C                          MATRIX CLU IN THE CALL, I.E.,    
C                              CALL DECOMP (N,NDIM,A,A,D,IPS,IER)    
C                          NOTE:  IF SUBROUTINE IMPRUV WILL BE USED, A  
C                                 MUST BE RETAINED.    
C    
C                        D    
C                          DETERMINANT OF A.    
C    
C                        IPS    
C                          VECTOR (OF DIMENSION AT LEAST N) STORING ROW 
C                          PERMUTATIONS RESULTING FROM PIVOTING, WHERE  
C                          IPS(K) = KTH PIVOT ROW AFTER ROW INTERCHANGES
C                          FROM PREVIOUS PIVOTS.    
C    
C                        IER    
C                          ERROR FLAG USED BY STANDARD ERROR MESSAGE    
C                          ROUTINE, ULIBER, AND RETURNED TO USER.    
C                          =  0,  NO ERROR.    
C                          = 32,  DETERMINANT OF A EQUALS ZERO, HENCE   
C                                 DECOMPOSITION IS IMPOSSIBLE.  A    
C                                 MESSAGE IS PRINTED.    
C    
C    
C    
C    
C    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NONE
C    
C PRECISION              COMPLEX*16
C    
C REQUIRED ULIB          NONE    
C ROUTINES    
C    
C SPECIALIST             R. K. SATO, NCAR, BOULDER, COLORADO  80303    
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                TRANSCRIBED FROM A PAPER BY C. B. MOLER,    
C                        LINEAR EQUATION SOLVER, COMM. ACM 15    
C                        (APRIL, 1972), 74.    
C    
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPS(NDIM)
C    
      COMPLEX*16 A(NDIM,NDIM),CLU(NDIM,NDIM)
      COMPLEX*16 T,D
C    
      IER = 0    
C    
C INITIALIZE ARRAYS    
C    
C    
      DO 101 J=1,N    
         IPS(J) = J    
  101 CONTINUE    
C    
      DO 103 J=1,N    
         DO 102 I=1,N    
            CLU(I,J) = A(I,J)    
  102    CONTINUE    
  103 CONTINUE    
  104 D = 1.    
C    
C LU DECOMPOSITION USING GAUSSIAN ELIMINATION WITH PARTIAL PIVIOTING.   
C    
      DO 110 K=1,N    
         IF (K .EQ. N) GO TO 109    
         KP1 = K+1    
         M = K    
C    
C FIND PIVOT ELEMENT.    
C    
         DO 105 I=KP1,N    
            IF (CDABS(CLU(I,K)) .GT. CDABS(CLU(M,K))) M = I    
  105    CONTINUE    
C    
C SAVE ROW PERMUTATIONS.    
C    
         IPS(K) = M    
         IF (M .NE. K) D = -D    
         T = CLU(M,K)    
         CLU(M,K) = CLU(K,K)    
         CLU(K,K) = T    
         IF (T .EQ. 0.0) GO TO 112    
C    
C CALCULATE AND STORE ELEMENTS OF FACTORS, L AND U.    
C    
         DO 106 I=KP1,N    
            CLU(I,K) = -CLU(I,K)/T    
  106    CONTINUE    
         DO 108 J=KP1,N    
            T = CLU(M,J)    
            CLU(M,J) = CLU(K,J)    
            CLU(K,J) = T    
            IF (T .EQ. 0.0) GO TO 108    
            DO 107 I=KP1,N    
               CLU(I,J) = CLU(I,J)+CLU(I,K)*T    
  107       CONTINUE    
  108    CONTINUE    
  109    IF (CLU(K,K) .EQ. 0.0) GO TO 112    
  110 CONTINUE    
C    
C COMPUTE DETERMINANT OF MATRIX A.    
C    
      DO 111 K=1,N    
         D = D*CLU(K,K)    
  111 CONTINUE    
      RETURN    
  112 IER = 32    
      D = 0.0    
      RETURN    
      END    
      SUBROUTINE SOLVE (N,NDIM,CLU,B,IPS,X)    
C    
C    
C DIMENSION OF           CLU(NDIM,N),B(N),IPS(N),X(N)    
C ARGUMENTS    
C    
C PURPOSE                SOLVES FOR VECTOR X OF UNKNOWNS OF THE LINEAR  
C                        SYSTEM AX = B USING THE LU DECOMPOSITION OF A  
C                        FROM SUBROUTINE DECOMP.    
C    
C USAGE                  CALL SOLVE (N,NDIM,CLU,B,IPS,X)    
C    
C ARGUMENTS    
C    
C ON INPUT               N    
C                          ORDER OF MATRIX A (NUMBER OF EQUATIONS).    
C    
C                        NDIM    
C                          DECLARED FIRST (ROW) DIMENSION OF CLU IN THE  
C                          CALLING PROGRAM.    
C    
C                        CLU    
C                          LU DECOMPOSITION MATRIX FROM SUBROUTINE    
C                          DECOMP.    
C    
C                        B    
C                          RIGHT HAND SIDE VECTOR (OF DIMENSION AT    
C                          LEAST N) OF CONSTANTS.    
C    
C                        IPS    
C                          VECTOR (OF DIMENSION AT LEAST N) FROM    
C                          SUBROUTINE DECOMP STORING ROW PERMUTATIONS   
C                          RESULTING FROM PIVOTING, WHERE IPS(K) = KTH  
C                          PIVOT ROW AFTER ROW INTERCHANGES FROM    
C                          PREVIOUS PIVOTS.    
C    
C ON OUTPUT              X    
C                          SOLUTION VECTOR (OF DIMENSION AT LEAST N).   
C                          IF VECTOR B NEED NOT BE SAVED, THEN THE    
C                          VECTOR B MAY BE USED AS THE SOLUTION VECTOR  
C                          IN THE CALL, I.E.,    
C                              CALL SOLVE (N,NDIM,CLU,B,IPS,B)    
C                          NOTE:  IF SUBROUTINE IMPRUV WILL BE USED, B  
C                                 MUST BE RETAINED.    
C    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NONE    
C    
C PRECISION              COMPLEX*16
C    
C REQUIRED ULIB          NONE    
C ROUTINES    
C    
C SPECIALIST             R. K. SATO, NCAR, BOULDER, COLORADO  80303    
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                TRANSCRIBED FROM A PAPER BY C. B. MOLER,    
C                        LINEAR EQUATION SOLVER, COMM. ACM 15    
C                        (APRIL, 1972), 74.    
C    
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPS(NDIM)
C   
      COMPLEX*16 CLU(NDIM,NDIM)
      COMPLEX*16 B(NDIM),X(NDIM)  
      COMPLEX*16 T  
C    
C    
C    
C    
C    
      DO 101 K=1,N    
         X(K) = B(K)    
  101 CONTINUE    
  102 IF (N .EQ. 1) GO TO 107    
      NM1 = N-1    
      DO 104 K=1,NM1    
         KP1 = K+1    
         M = IPS(K)    
         T = X(M)    
         X(M) = X(K)    
         X(K) = T    
         DO 103 I=KP1,N    
            X(I) = X(I)+CLU(I,K)*T    
  103    CONTINUE    
  104 CONTINUE    
      DO 106 KB=1,NM1    
         KM1 = N-KB    
         K = KM1+1    
         X(K) = X(K)/CLU(K,K)    
         T = -X(K)    
         DO 105 I=1,KM1    
            X(I) = X(I)+CLU(I,K)*T    
  105    CONTINUE    
  106 CONTINUE    
  107 X(1) = X(1)/CLU(1,1)    
      RETURN    
      END    
      SUBROUTINE IMPRUV (N,NDIM,A,CLU,B,IPS,X,R,DIGITS,IER)    
C    
C    
C DIMENSION OF           A(NDIM,N),CLU(NDIM,N),B(N),IPS(N),X(N),R(N)    
C ARGUMENTS    
C    
C PURPOSE                ITERATIVE IMPROVEMENT OF THE SOLUTION OBTAINED 
C                        FROM SUBROUTINE SOLVE.    
C    
C USAGE                  CALL IMPRUV (N,NDIM,A,CLU,B,IPS,X,R,DIGITS,IER) 
C    
C ARGUMENTS    
C    
C ON INPUT               N    
C                          ORDER OF MATRIX A (NUMBER OF EQUATIONS).    
C    
C                        NDIM    
C                          DECLARED FIRST (ROW) DIMENSIONS OF A AND CLU  
C                          IN THE CALLING PROGRAM.    
C    
C                        A    
C                          ORIGINAL COEFFICIENT MATRIX.    
C    
C                        CLU    
C                          LU DECOMPOSITION MATRIX OF A FROM SUBROUTINE 
C                          DECOMP.    
C    
C                        B    
C                          ORIGINAL RIGHT HAND SIDE VECTOR (OF DIMENSION
C                          AT LEAST N) OF CONSTANTS.    
C    
C                        IPS    
C                          VECTOR (OF DIMENSION AT LEAST N) FROM    
C                          SUBROUTINE DECOMP STORING ROW PERMUTATIONS   
C                          RESULTING FROM PIVOTING, WHERE IPS(K) = KTH  
C                          PIVOT ROW AFTER ROW INTERCHANGES FROM    
C                          PREVIOUS PIVOTS.    
C    
C                        X    
C                          SOLUTION VECTOR (OF DIMENSION AT LEAST N)    
C                          FROM SUBROUTINE SOLVE.    
C    
C                        R    
C                          WORK ARRAY (OF DIMENSION AT LEAST N) USED IN 
C                          THE COMPUTATION OF RESIDUALS.    
C    
C ON OUTPUT              X    
C                          SOLUTION VECTOR AFTER ITERATIVE IMPROVEMENT. 
C    
C                        DIGITS    
C                          APPROXIMATE NUMBER OF CORRECT DIGITS IN X.   
C    
C                        IER    
C                          ERROR FLAG USED BY STANDARD ERROR MESSAGE    
C                          ROUTINE ULIBER AND RETURNED TO USER.    
C                          = 0,  NO ERROR.    
C                          = 1,  NO CONVERGENCE IN THE ITERATION.    
C    
C COMMON BLOCKS          NONE    
C    
C I/O                    NOME
C    
C PRECISION              COMPLEX*16
C    
C REQUIRED ULIB          NONE    
C ROUTINES    
C    
C SPECIALIST             R. K. SATO, NCAR, BOULDER, COLORADO  80303    
C    
C LANGUAGE               FORTRAN    
C    
C HISTORY                TRANSCRIBED FROM FORSYTHE AND MOLER, COMPUTER  
C                        SOLUTIONS OF LINEAR ALGEBRAIC SYSTEMS, 1967.   
C    
C    
C    
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPS(NDIM)
C
      REAL*8 EPS,DIGITS,XNORM,DXNORM
C
      COMPLEX*16 A(NDIM,NDIM),CLU(NDIM,NDIM)    
      COMPLEX*16 B(NDIM),X(NDIM),R(NDIM)
      COMPLEX*16 SUM,T   
C
C      
C EPS AND ITMAX ARE MACHINE DEPENDENT CONSTANTS.    
C EPS = LARGEST FLOATING POINT NUMBER SUCH THAT 1.0 + EPS = 1.0    
C ITMAX = NUMBER OF BINARY BITS IN THE MANTISSA OF A FLOATING POINT    
C         NUMBER AND IS THE MAXIMUM NUMBER OF ITERATIONS.    
      DATA EPS,ITMAX/1.E-15,48/    
C    
      IER = 0    
      XNORM = 0.0    
      DO 101 I=1,N    
         XNORM = DMAX1(XNORM,CDABS(X(I)))    
  101 CONTINUE    
      IF (XNORM .GT. 0.0) GO TO 102    
      DIGITS = -DLOG10(EPS)    
      GO TO 108    
  102 DO 107 ITER=1,ITMAX    
C    
C CALCULATION OF DOUBLE PRECISION RESIDUALS.    
C    
         DO 104 I=1,N    
            SUM = 0.0    
            DO 103 J=1,N    
               SUM = SUM+A(I,J)*X(J)    
  103       CONTINUE    
            SUM = B(I)-SUM    
            R(I) = SUM    
  104    CONTINUE    
         CALL SOLVE (N,NDIM,CLU,R,IPS,R)    
         DXNORM = 0.0    
         DO 105 I=1,N    
            T = X(I)    
            X(I) = X(I)+R(I)    
            DXNORM = DMAX1(DXNORM,CDABS(X(I)-T))    
  105    CONTINUE    
         IF (ITER .GT. 1) GO TO 106    
         DIGITS = -DLOG10(DMAX1(DXNORM/XNORM,EPS))    
  106    IF (DXNORM .LE. EPS*XNORM) GO TO 108    
  107 CONTINUE    
C    
C ITERATION DID NOT CONVERGE    
C    
      IER = 1    
  108 CONTINUE    
      RETURN    
C    
C REVISION HISTORY---    
C    
C JUNE 1977        REMOVED REFERENCE TO THE  LOC  FUNCTION    
C                  TO ENHANCE PORTABILITY.    
C    
C JANUARY 1978     DELETED REFERENCES TO THE  *COSY  CARDS, MOVED    
C                  THE REVISION HISTORIES TO APPEAR BEFORE THE    
C                  FINAL END CARD, AND MOVED THE INITIAL COMMENT    
C                  CARDS TO APPEAR AFTER THE FIRST SUBROUTINE CARD    
C JANUARY 1979     REMOVED COMMENTS ABOUT LOC FUNCTION.    
C-----------------------------------------------------------------------
      END    
