      COMPLEX*16 FUNCTION CEXIB(H1,H2,Z1,Z2,DEPB,RHOB,GRHOB,
     1                          ALPHA1,ALPHA2,BETA1,BETA2)
C
C     THIS FUNCTION EVALUATES THE INTEGRAL OF TWO BASIC DEPTH 
C     FUNCTIONS OVER THE INTERVAL Z1 TO Z2 IN THE BOTTOM, 
C     WEIGHTED BY AN EXPONENTIAL DENSITY PROFILE. IT CAN BE 
C     USED WITHIN ONE SET OF BASIC DEPTH FUNCTIONS OR IN THE
C     COMPUTATION OF CROSS COUPLING INTEGRALS BETWEEN DIFFERENT
C     SETS OF BASIC DEPTH FUNCTIONS.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMPLEX*16 ALPHA1,ALPHA2,BETA1,BETA2
      COMPLEX*16 B1,B2,CFAC,CGRHOB,CI
      COMPLEX*16 U,V,W,X
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
      CI=DCMPLX(0.0D0,1.0D0)
      CGRHOB=DCMPLX(GRHOB,0.0D0)
      THBMH1=2.0D0*HB-H1
      THBMH2=2.0D0*HB-H2
C
      CEXIB=CMPLX(0.0D0,0.0D0)
      IF(Z2 .LE. Z1) GO TO 200
C
      CFAC=FKB**4*CDSIN(ALPHA1*H1)*CDSIN(ALPHA2*H2)/
     1            (RHOB*ALPHA1*ALPHA2)
C
      B1=2.0D0*BETA1*(HB-H1)
      IF(DIMAG(B1) .LT. 88.0D0) THEN
      CFAC=CFAC/(1.0D0-CDEXP(CI*B1))
      END IF
C
      B2=2.0D0*BETA2*(HB-H2)
      IF(DIMAG(B2) .LT. 88.0D0) THEN
      CFAC=CFAC/(1.0D0-CDEXP(CI*B2))
      END IF
C
      W=CI*(BETA1+BETA2)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA1*(Z2-H1)+BETA2*(Z2-H2))-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA1*(Z1-H1)+BETA2*(Z1-H2))-CGRHOB*(Z1-DEPB))
      CEXIB=CEXIB+(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA1*(DEPB-H1)+BETA2*(DEPB-H2)))
      CEXIB=CEXIB+X*(Z2-Z1)
      END IF
C
      W=CI*(BETA1-BETA2)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA1*(Z2-H1)+BETA2*(THBMH2-Z2))-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA1*(Z1-H1)+BETA2*(THBMH2-Z1))-CGRHOB*(Z1-DEPB))
      CEXIB=CEXIB-(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA1*(DEPB-H1)+BETA2*(THBMH2-DEPB)))
      CEXIB=CEXIB-X*(Z2-Z1)
      END IF
C
      W=CI*(BETA2-BETA1)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA2*(Z2-H2)+BETA1*(THBMH1-Z2))-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA2*(Z1-H2)+BETA1*(THBMH1-Z1))-CGRHOB*(Z1-DEPB))
      CEXIB=CEXIB-(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA2*(DEPB-H2)+BETA1*(THBMH1-DEPB)))
      CEXIB=CEXIB-X*(Z2-Z1)
      END IF
C
      W=-CI*(BETA1+BETA2)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA1*(THBMH1-Z2)+BETA2*(THBMH2-Z2))-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA1*(THBMH1-Z1)+BETA2*(THBMH2-Z1))-CGRHOB*(Z1-DEPB))
      CEXIB=CEXIB+(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA1*(THBMH1-DEPB)+BETA2*(THBMH2-DEPB)))
      CEXIB=CEXIB+X*(Z2-Z1)
      END IF
C
      CEXIB=CFAC*CEXIB
C
  200 CONTINUE
C
      RETURN
      END