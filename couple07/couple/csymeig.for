      SUBROUTINE CSYMEIG(a,b,u,v,cs,z1,z2,d,e,n,eivec,nx,nxsd2)
c
c From ALGOL program in Handbook Series Linear Algebra: "A Jacobi type
c method for complex symmetric matrices" by P. J. Anderson and G. Loizou,
c Numer. Math. 25, 347-363 (1976). Converted to FORTRAN by K. Santinello
c and R. B. Evans, ODSI, North Stonington, CT (1987).  
c
c Solves by a Jacobi-like method the eigenproblem for the n x n complex
c symmetric matrix A + Bi, whose diagonal and super-diagonal elements
c are stored by column in arrays a,b[1:n(n+1)/2]. If eivec is TRUE then
c the eigenvectors are generated and the arrays u,v[1:n,1:n] must be
c provided to receive them. If eivec is FALSE then u,v are not required
c and can be trival in the procedure call. Iterations are limited to
c 15. Eps is set to 10**-16 in the program. It is used to test for 
c convergence and should be no smaller than the machine accurracy.
c The arrays are dimensioned for nx, which should be greater than or 
c equal to n, and nxsd2=nx*(nx+1)/2.
c
c The calling sequence to find the n complex eigenvalues CEIG (and the 
c n complex eigenvectors, stored as colunms in CVECT) of the symmetric
c complex matrix C is as follows:
c
c      REAL*8 a(nxsd2),b(nxsd2),u(nx,nx),v(nx,nx)
c      REAL*8 cs(nx),z1(nx),z2(nx),d(nx),e(nx)
c      LOGICAL eivec
c      ....
c
c      DO 30 j=1,n
c         DO 20 i=1,j
c            jji = (j*(j-1)/2)+i
c            a(jji) = DREAL(C(i,j))
c            b(jji) = DIMAG(C(i,j))
c   20    CONTINUE
c   30 CONTINUE
c
c      eivec = .TRUE.
c      CALL CSYMEIG(a,b,u,v,cs,z1,z2,d,e,n,eivec,nx,nxsd2) 
c
c      DO 50 i=1,n
c         ii=i*(i+1)/2
c         CEIG(i)=DCMPLX(a(ii),b(ii))
c   50 CONTINUE
c
c      DO 70 i=1,n
c         DO 60 ii=1,n
c            CVECT(ii,i)=DCMPLX(u(ii,i),v(ii,i))
c   60    CONTINUE
c   70 CONTINUE
c
      INTEGER n,nx,nxsd2
      INTEGER i,j,k,p,q,r,it,pq,pp,qq,ip,iq,ij,ii,jj,n1,ip1
      LOGICAL eivec
      LOGICAL notrans,wlarge,wsmall
      REAL*8 a(nxsd2),b(nxsd2),u(nx,nx),v(nx,nx)
      REAL*8 cs(nx),z1(nx),z2(nx),d(nx),e(nx)
      REAL*8 eps,sqrt3,cos45,sk,t1,max,trhold,app,bpp,aqq,bqq,apq,bpq,
     >       a1,a2,a3,a4,a5,aip,bip,aiq,biq,d1,d2,b1,b2,w,delta,fx,den,
     >       xr1,xr2,cy,sy,c2y,s2y,cx,sx,c2x,s2x,cot4x,cot2x,csec2x,
     >       cotx,t2,t3,uip,vip,uiq,viq
c
c
      eps = 1.0D-16
      n1 = n*(n-1)/2
      sqrt3 = DSQRT(3.0D0)
      cos45 = 1.0D0/DSQRT(2.0D0)
c
c Initialize u and v
c
      IF (eivec) THEN
	DO 10 i=1,n
	   u(i,i) = 1.0D0
	   v(i,i) = 0.0D0
	   ip1=i+1
	   DO 20 j=ip1,n
	      u(i,j) = 0.0D0
	      u(j,i) = 0.0D0
	      v(i,j) = 0.0D0
	      v(j,i) = 0.0D0
20         CONTINUE 
10      CONTINUE
      END IF
      notrans = .FALSE.
c 
c Iteration loop
c
      DO 500 it=1,15
	 IF (notrans) GOTO 900
	 sk = 0.0D0
	 DO 40 i=1,n
	    t1 = 0.0D0
	    DO 50 j=1,n
	       IF (i .NE. j) THEN 
		  IF (i .LT. j) THEN 
		     ij = j*(j-1)/2+i
		  ELSE 
		     ij = i*(i-1)/2+j
		  END IF
		  t1 = DABS(a(ij)) + DABS(b(ij)) + t1
	       END IF
50          CONTINUE               
	       sk = sk + t1
	       ii = i*(i+1)/2
	       cs(i) = t1 + DABS(a(ii)) + DABS(b(ii))
40       CONTINUE
c
c Interchange rows and columns
c
	 nm1 = n-1
	 DO 60 p=1,nm1
	    max = cs(p)
	    j = p
	    ip1 = p+1
	    DO 70 i=ip1,n
	       IF (cs(i) .GT. max) THEN
		  max = cs(i)
		  j = i
	       END IF
70          CONTINUE
	    IF (j .NE. p) THEN
	       cs(j) = cs(p)
	       DO 80 i=1,n
		  IF ((i .NE. j) .AND. (i .NE. p)) THEN
		     IF (i .LT. p) THEN
			ip=p*(p-1)/2+i
		     ELSE 
			ip=i*(i-1)/2+p
		     END IF  
		     IF (i .LT. j) THEN 
			ij=j*(j-1)/2+i
		     ELSE 
			ij=i*(i-1)/2+j
		     END IF
		     t1 = a(ip)
		     a(ip) = a(ij)
		     a(ij) = t1
		     t1 = b(ip)
		     b(ip) = b(ij)
		     b(ij) = t1
		  END IF
		  IF (eivec) THEN
		     t1 = u(i,p)
		     u(i,p) = u(i,j)
		     u(i,j) = t1
		     t1 = v(i,p)
		     v(i,p) = v(i,j)
		     v(i,j) = t1
		  END IF
80             CONTINUE
	       pp = p*(p+1)/2
	       jj = j*(j+1)/2
	       t1 = a(pp)
	       a(pp) = a(jj)
	       a(jj) = t1
	       t1 = b(pp)
	       b(pp) = b(jj)
	       b(jj) = t1                  
	    END IF
60       CONTINUE
C
c
c Decide whether off-diagonal elements sufficiently small
c
	 IF (sk .LT. (10000.0D0*eps)) GO TO 900
c
c Diagonal elements are in d,e
c
	 DO 90 i=1,n
	    ii = i*(i+1)/2
	    d(i) = a(ii)
	    e(i) = b(ii)
	    z1(i) = 0.0D0
	    z2(i) = 0.0D0
90       CONTINUE
c
c Sweep starts
c
	 notrans = .TRUE.        
	 k = 0
	 r = 0
	 trhold = sk/(1.0D1*n1)
	 IF (it .GT. 5) trhold = 1.0D-2*trhold
	 p = 1
	 q = 2
c
c (r < n1) loop
c
100       CONTINUE
	    pq = q*(q-1)/2+p
	    IF (it .LT. 6) k=k+1
	    IF ((DABS(a(pq)) + DABS(b(pq))) .GE. trhold) THEN
c 
c Pivoting on (p,q)
c       
	       r = r + 1
	       pp = p*(p+1)/2
	       qq = q*(q+1)/2
	       pq = q*(q-1)/2+p
	       app = a(pp)
	       bpp = b(pp)
	       aqq = a(qq)
	       bqq = b(qq)
	       apq = a(pq)
	       bpq = b(pq)
	       d1 = app-aqq
	       d2 = bpp-bqq
	       b1 = 2.0D0*apq
	       b2 = 2.0D0*bpq
	       a2 = 0.0D0
	       a4 = 0.0D0
	       DO 105 i=1,n
		  IF ((i .NE. p) .AND. (i .NE. q)) THEN
		     IF (i .LT. p) THEN 
			ip = p*(p-1)/2+i
		     ELSE 
			ip = i*(i-1)/2+p
		     END IF
		     IF (i .LT. q) THEN
			iq = q*(q-1)/2+i
		     ELSE 
			iq = i*(i-1)/2+q
		     END IF
		     aip = a(ip)
		     bip = b(ip)
		     aiq = a(iq)
		     biq = b(iq)
		     a2 = a2 + (aip-biq)**2 + (bip+aiq)**2
		     a4 = a4 + bip*aiq-aip*biq
		  END IF
105            CONTINUE
		  a1 = (d1-b2)**2 + (d2+b1)**2
		  a5 = 4.0D0*(d2*b1-d1*b2) + 8.0D0*a4
		  a3 = 6.0D0*(a1+a2)
		  a4 = 4.0D0*(a1+a2) + 8.0D0*a4
		  a2 = 4.0D0*a1+2.0D0*a2
		  t3 = a1*16.0D0+a2*8.0D0+a3*4.0D0+a4*2.0D0+a5
		  wlarge = .TRUE.
		  IF (t3 .GT. 0.0D0) wlarge =  .FALSE.
		  t3 = a1*16.0D0-a2*24.0D0+a3*36.0D0-a4*54.0D0+a5*81.0D0
		  wsmall = .TRUE.
		  IF (t3 .LT. 0.0D0) wsmall = .FALSE.
		  IF (wlarge .or. wsmall) THEN
		     cy = 2.0D0/sqrt3
		     sy = -cy/2.0D0
		     c2y = 5.0D0/3.0D0
		     s2y = -4.0D0/3.0D0
		     IF (wlarge) THEN 
			sy = -sy
			s2y = -s2y
		     END IF
		  ELSE
c
c Solve for w by Newton-Raphson
c
		     i = 0
		     w = 2.0D0
		     delta = 1.0D0
115                     CONTINUE
			i = i + 1
			fx = (((a1*w+a2)*w+a3)*w+a4)*w+a5
			den = ((4.0D0*a1*w+3.0D0*a2)*w+2.0D0*a3)*w+a4
			delta = -fx/den
			w = w + delta
			IF ((DABS(delta) .GT. (1.0D1*eps*DABS(w)))
     >                        .AND.  (i .LT. 30)) GO TO 115
		     IF (DABS(w) .LT. eps) THEN
			sy = 0.0D0
			s2y = 0.0D0
			cy = 1.0D0
			c2y = 1.0D0
		     ELSE
			t1 = 2.0D0*DSQRT(1.0D0+W)
			sy = w/t1
			cy = (2.0D0+w)/t1
			t1 = 2.0D0*w+w*w
			t2 = 2.0D0*(1.0D0+w)
			s2y = t1/t2 
			c2y = (2.0D0+t1)/t2
		     END IF
		  END IF
		  xr1 = 2.0D0*(b1*d1+b2*d2)
		  xr2 = d1*d1+d2*d2-b1*b1-b2*b2
		  IF (DABS(xr1) .LE. eps*DABS(xr2)) THEN
		     IF (xr2 .GE. 0.0D0) THEN
c 
c x is small
c
			IF (sy .EQ. 0.0D0) THEN
			   GO TO 121
			ELSE
			   sx = 0.0D0
			   s2x = 0.0D0
			   cx = 1.0D0
			   c2x = 1.0D0
			END IF
		     ELSE
c
c x is near pi/4 or -pi/4
c
			c2x = 0.0D0
			IF (xr1 .GE. 0.0D0) THEN 
			   s2x = 1.0D0
			ELSE 
			   s2x = -1.0D0
			END IF
			cx = cos45
			sx = s2x*cos45
		     END IF
		  ELSE
		     cot4x = xr2/xr1
		     cot2x = DABS(cot4x) + DSQRT(1.0D0+cot4x*cot4x)
		     IF (xr2 .LT. 0.0D0) cot2x = 1.0D0/cot2x
		     csec2x = DSQRT(1.0D0+cot2x*cot2x)
		     s2x = 1.0D0/csec2x
		     c2x = s2x*cot2x
		     cotx = cot2x+csec2x
		     sx = 1.0D0/DSQRT(1.0D0+cotx*cotx)
		     cx = sx*cotx
		     IF (xr1 .LT. 0.0D0) THEN
			s2x = -s2x
			sx = -sx
		     END IF
		  END IF  
c
c transformation section
c
		  notrans = .FALSE.
		  k = 0
		  DO 110 i=1,n
		     IF ((i .NE. p) .AND. (i .NE. q)) THEN
			IF (i .LT. p) THEN
			   ip = p*(p-1)/2+i
			ELSE 
			   ip = i*(i-1)/2+p
			END IF
			IF (i .LT. q) THEN
			   iq = q*(q-1)/2+i
			ELSE 
			   iq = i*(i-1)/2+q
			END IF
			aip = a(ip)
			bip = b(ip)
			aiq = a(iq)
			biq = b(iq)
			a(ip) = (aip*cy-biq*sy)*cx+(aiq*cy+bip*sy)*sx
			b(ip) = (bip*cy+aiq*sy)*cx+(biq*cy-aip*sy)*sx
			a(iq) = (aiq*cy+bip*sy)*cx-(aip*cy-biq*sy)*sx
			b(iq) = (biq*cy-aip*sy)*cx-(bip*cy+aiq*sy)*sx
		     END IF
110               CONTINUE
		     a(pq) =((b1*c2y+d2*s2y)*c2x-(d1*c2y-b2*s2y)*s2x)/2.0D0
		     b(pq) =((b2*c2y-d1*s2y)*c2x-(d2*c2y+b1*s2y)*s2x)/2.0D0
		     t1 =(-2.0D0*d1*sx*sx+(d1*2.0D0*sy*sy-b2*s2y)*c2x+
     >                   (b1*c2y+d2*s2y)*s2x)/2.0D0
		     t2 =(-2.0D0*d2*sx*sx+(d2*2.0D0*sy*sy+b1*s2y)*c2x+
     >                   (b2*c2y-d1*s2y)*s2x)/2.0D0
		     a(pp) = app + t1
		     b(pp) = bpp + t2
		     a(qq) = aqq - t1
		     b(qq) = bqq - t2
		     z1(p) = z1(p) + t1
		     z1(q) = z1(q) - t1
		     z2(p) = z2(p) + t2
		     z2(q) = z2(q) - t2
		     IF (eivec) GO TO 119
		     GO TO 121
119                  CONTINUE
		     DO 120 i=1,n
			uip = u(i,p)
			uiq = u(i,q)
			vip = v(i,p)
			viq = v(i,q)
			u(i,p)=(uip*cy-viq*sy)*cx+(uiq*cy+vip*sy)*sx
			v(i,p)=(vip*cy+uiq*sy)*cx+(viq*cy-uip*sy)*sx
			u(i,q)=(uiq*cy+vip*sy)*cx-(uip*cy-viq*sy)*sx
			v(i,q)=(viq*cy-uip*sy)*cx-(vip*cy+uiq*sy)*sx
120                  CONTINUE
121         CONTINUE
	    END IF
	    IF (k .eq. n1) THEN  
	       trhold = 1.0D-2*trhold
	       k = 0
	    END IF
	    IF (q .LT. n) THEN
		q = q + 1
	    ELSE 
	       IF (p .LT. (n-1)) THEN
		  p = p + 1
		  q = p + 1
	       ELSE
		  p = 1
		  q = 2
		  IF (it .GT. 5) trhold = 1.0D-4*trhold
	       END IF
	    END IF
	 IF (r .LT. n1) GO TO 100     
c
c diagonal elements are corrected
c
	 DO 130 i=1,n
	    ii = i*(i+1)/2
	    a(ii) = d(i) + z1(i)
	    b(ii) = e(i) + z2(i)
130      CONTINUE
C
C
500   CONTINUE         
900   CONTINUE
      RETURN
      END
