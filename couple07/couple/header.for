      SUBROUTINE HEADER(TITLE,FREQ,ZS,NDEP,NDINC,ZR,N,RANGE,DEPTH,MD,
     1                  NRX,BSOUT)
C
C     THIS SUBROUTINE WRITES THE HEADERS ON THE OUTPUT FILES.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZR(MD)
      DIMENSION RANGE(NRX),DEPTH(NRX)
      INTEGER BSOUT
C
      CHARACTER*10 TITLE(8),DA
      CHARACTER*24 date_time
C
C    FIND DATE IN G95 or Unix
C               
c      call fdate(date_time)
c      da=date_time(5:7)//date_time(9:10)//','//date_time(21:24)
C
C    FIND DATE IN LAHEY
C               
C
C YIQING: 
C THERE IS NO 'DATE' IN GFORTRAN
C USE FDATE INSTEAD
C https://gcc.gnu.org/onlinedocs/gfortran/FDATE.html
C      CALL DATE(DA) 
      CALL FDATE(DA)
C
C     FIND DATE IN SVS
C
C      CALL GETDAT(IY,IM,ID)
C     WRITE(DA,50) IM,ID,IY
C   50 FORMAT(I2,'-',I2,'-',I4)    
C
      NR=0
C
C     UNFORMATTED COMPLEX PRESSURE: TOTAL,OUTGOING,INGOING
C     AND SCATTERED.
C     
      WRITE(7) (TITLE(I),I=1,8),DA,SNGL(FREQ),SNGL(ZS),NR,NDEP,
     1         (SNGL(ZR(I)),I=1,NDEP)
C
      IF(BSOUT .eq. 0) THEN
      WRITE(14) (TITLE(I),I=1,8),DA,SNGL(FREQ),SNGL(ZS),NR,NDEP,
     1         (SNGL(ZR(I)),I=1,NDEP)
      WRITE(15) (TITLE(I),I=1,8),DA,SNGL(FREQ),SNGL(ZS),NR,NDEP,
     1         (SNGL(ZR(I)),I=1,NDEP)
      WRITE(16) (TITLE(I),I=1,8),DA,SNGL(FREQ),SNGL(ZS),NR,NDEP,
     1         (SNGL(ZR(I)),I=1,NDEP)
      END IF
C
C      THE FOLLOWING RECORD WAS USED FOR COUTOUR PLOTS AND IS
C      CURRENTLY INACTIVE
C      WRITE(7) N,(SNGL(RANGE(I)),I=1,N),(SNGL(DEPTH(I)),I=1,N)
C
C     FORMATTED TRANSMISSION LOSS
C
      XND=FLOAT(NDEP)
      XNI=FLOAT(NDINC)
      ND=XND/XNI+.999
      WRITE(10,100) (TITLE(I),I=1,8)
  100 FORMAT(8A10)
      WRITE(10,200) DA,FREQ,ZS,NR,ND
  200 FORMAT(A10,2F13.7,2I5)
      WRITE(10,300) (ZR(I),I=1,NDEP,NDINC)
  300 FORMAT(20F13.7)
C
      IF(BSOUT .eq. 0) THEN
C
C     OUTGOING TL
C
      WRITE(12,100) (TITLE(I),I=1,8)
      WRITE(12,200) DA,FREQ,ZS,NR,ND
      WRITE(12,300) (ZR(I),I=1,NDEP,NDINC)
C
C     INGOING TL
C
      WRITE(13,100) (TITLE(I),I=1,8)
      WRITE(13,200) DA,FREQ,ZS,NR,ND
      WRITE(13,300) (ZR(I),I=1,NDEP,NDINC)
C
C     SCATTERED TL
C
      WRITE(17,100) (TITLE(I),I=1,8)
      WRITE(17,200) DA,FREQ,ZS,NR,ND
      WRITE(17,300) (ZR(I),I=1,NDEP,NDINC)
C
      END IF
C
      RETURN
      END
