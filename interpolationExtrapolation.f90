SUBROUTINE polint(xa,ya,x,y,dy)
USE nrtype; USE nrutil, ONLY: assert_eq,iminloc,nrerror
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(IN):: xa,ya
REAL(SP),INTENT(IN):: x
REAL(SP),INTENT(OUT):: y,dy
! Given arrays xa and ya of length N, and given a value x, this
! routine returns a value y, and an error estimate dy. If P(x)
! is the polynomial of degree N-1 such that P(xa_i)=ya_i,i=1,..,N,
! then the returned value y=P(x)
INTEGER(I4B):: m,n,ns
REAL(SP),DIMENSION(size(xa)):: c,d,den,ho
n=assert_eq(size(xa),size(ya),'polint')
c=ya			! Initialize the tableau of c's and d's
d=ya
ho=xa-x
ns=iminloc(abs(x-xa))	! Find index ns of closest table entry.
y=ya(ns)				! This is the initial approximation to y.
ns=ns-1
do m=1,n-1			! For each column of the tableau
	den(1:n-m)=ho(1:n-m)-ho(1+m:n)
	if (any(den(1:n-m) == 0.0)) &
		call nrerror('polint: calculation failure')
		! This error can occur only if two input xa's are (to within roundoff) identical.
	den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
	d(1:n-m)=ho(1+m:n)*den(1:n-m)	! Here the c's and d's are updated
	c(1:n-m)=ho(1:n-m)*den(1:n-m)
	if (2*ns < n-m) then	! After each column in the tableau is completed,
		dy=c(ns+1)			! we decide which correction, c or d, we want to
	else					! add to our accumulating value of y, i.e., which
		dy=d(ns)			! path to take through the tableau - forking up or 
		ns=ns-1				! down. We do this in such a way as to take the most
	end if					! "straight line" route through the tableau to its apex,
	y=y+dy					! updating ns accordingly to keep track of where we are.
end do						! This route keeps the partial approximations centered 
END SUBROUTINE polint		! (insofar as possible) on the target x. The last dy added is thus the error indication.


SUBROUTINE ratint(xa,ya,x,y,dy)
USE nrtype; USE nrutil, ONLY: assert_eq,iminloc,nrerror
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(IN):: xa,ya
REAL(SP),INTENT(IN):: x
REAL(SP),INTENT(OUT):: y,dy
! Given array xa and ya of length N, and given a value of x, this
! routine returns a value of y and an accuracy estimate dy. The value
! returned is that of the diagonal rational function, evaluated at x,
! that passes through the N points (xa_i,ya_i),i=1,..,N
INTEGER(I4B):: m,n,ns
REAL(SP),DIMENSION(size(xa)):: c,d,dd,h,t
REAL(SP),PARAMETER:: TINY=1.0e-25_sp	! A small number.
n=assert_eq(size(xa),size(ya),'ratint')
h=xa-x
ns=iminloc(abs(h))
y=ya(ns)
if (x == xa(ns)) then
	dy=0.0
	RETURN
end if
c=ya
d=ya+TINY
ns=ns-1					! The TINY part is needed to prevent a
do m=1,n-1				! rare zero-over-zero condition.
	t(1:n-m)=(xa(1:n-m)-x)*d(1:n-m)/h(1+m:n)	! h will never be zero
	dd(1:n-m)=t(1:n-m)-c(2:n-m+1)	! since this was tested in the initializing loop
	if (any(dd(1:n-m) == 0.0)) &
		call nrerror('failure in ratint')
	dd(1:n-m)=(c(2:n-m+1)-d(1:n-m))/dd(1:n-m)	! The error condition indicates
	d(1:n-m)=c(2:n-m+1)*dd(1:n-m)	! that the interpolating function has a pole
	c(1:n-m)=t(1:n-m)*dd(1:n-m)		! at the requested value of x.
	if (2*ns < n-m) then
		dy=c(ns+1)
	else
		dy=d(ns)
		ns=ns-1
	end if
	y=y+dy
end do
END SUBROUTINE ratint

