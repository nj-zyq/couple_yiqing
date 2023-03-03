      SUBROUTINE ABORTC(IFLAG)       
C    
C     THIS SUBROUTINE ABORTS COUPLE FOR THE FOLLOWING REASONS:
C
C      1. Non-Physical Eigenvalue found in GALRKN   
C      2. ZERO DETERMINANT RETURNED BY DECOMP IN CLINEQS
C      3. FAILURE TO IMPROVE BY IMPRUV IN CLINEQS ON 
C         RETURN TO COUPLE
C    
      IMPLICIT REAL*8 (A-H,O-Z)
      WRITE(9,100) IFLAG
  100 FORMAT(1X,///,' ***** RUN ABORTED FOR REASON ',I5,
     1       '  *****',////)
C
      WRITE(9,200)
  200 FORMAT(5X,'1. Non-Physical Eigenvalue found in GALRKN',/
     1       5X,'2. ZERO DETERMINANT RETURNED BY DECOMP',/
     2       5X,'3. FAILURE TO INPROVE BY IMPRUV',/)
C              
      STOP
      END
