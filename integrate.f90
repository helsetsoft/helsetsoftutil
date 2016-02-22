SUBROUTINE trapzd(func,a,b,s,n)
USE nrtype; USE nrutil, ONLY: arth
IMPLICIT NONE
REAL(SP),INTENT(IN):: a,b
REAL(SP),INTENT(INOUT):: s
INTEGER(I4B),INTENT(IN):: n
INTERFACE
	FUNCTION func(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	END FUNCTION func
END INTERFACE
! This routine computes the nth stage of refinement of an extended trapeziodal
! rule. func is input as the name of the function to be integrated between
! limits a and b, also input. When called with n=1, the routine returns as s 
! the crudest estimate of int_a^b f(x)dx. Subsequent calls with n=2,3,.. (in
! that sequential order) will improve the accuracy of s by adding 2^(n-2) additional
! interior points, a should not be modified between sequential calls.
REAL(SP):: del,fsum
INTEGER(I4B):: it
if(n==1) then
	s=0.5_sp*(b-a)*sum(func((/ a,b /)))
else
	it=2**(n-2)
	del=(b-a)/it	! This is the spacing of the points to be added.
	fsum=sum(func(arth(a+0.5_sp*del,del,it)))
	s=0.5_sp*(s+del*fsum)	! This replaces s by its refined value.
end if
END SUBROUTINE trapzd


FUNCTION qtrap(func,a,b)
USE nrtype; USE nrutil, ONLY: nrerror
IMPLICIT NONE
REAL(SP),INTENT(IN):: a,b
REAL(SP):: qtrap
INTERFACE
	FUNCTION func(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	END FUNCTION func
END INTERFACE
INTEGER(I4B),PARAMETER:: JMAX=20
REAL(SP),PARAMETER:: EPS=1.0e-6_sp, UNLIKELY=-1.0e30_sp
! Returns the integral of the function func from a to b. The
! parameter EPS should be set to the desired fractional accuracy
! and JMAX so that 2 to the power JMAX-1 is the maximum allowed
! number of steps. Integration is performed by the trapezoidal rule.
REAL(SP):: olds
INTEGER(I4B):: j	! Any number that is unlikely to be the average
olds=UNLIKELY		! of the function at its endpoints will do here.
do j=1,JMAX
	call trapzd(func,a,b,qtrap,j)
	if (j>5) then	! Avoid spurious early convergence.
		if (abs(qtrap-olds) < EPS*abs(olds) .or. &
			(qtrap == 0.0 .and. olds == 0.0)) RETURN
	end if
	olds=qtrap
end do
call nrerror('qtrap: too many steps')
END FUNCTION qtrap


FUNCTION qsimp(func,a,b)
USE nrtype; USE nrutil, ONLY: nrerror
IMPLICIT NONE
REAL(SP),INTENT(IN):: a,b
REAL(SP):: qsimp
INTERFACE
	FUNCTION func(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	END FUNCTION func
END INTERFACE
INTEGER(I4B),PARAMETER:: JMAX=20
REAL(SP),PARAMETER:: EPS=1.0e-6_sp, UNLIKELY=-1.0e30_sp
! Returns the integral of the function func from a to b. The
! parameter EPS should be set to the desired fractional accuracy
! and JMAX so that 2 to the power JMAX-1 is the maximum allowed
! number of steps. Integration is performed by Simpson's rule.
INTEGER(I4B):: j
REAL(SP):: os,ost,st
ost=UNLIKELY
os=UNLIKELY
do j=1,JMAX
	call trapzd(func,a,b,st,j)
	qsimp=(4.0_sp*st-ost)/3.0_sp	! Compare equation (4.2.4)
	if (j>5) then					! Avoid spurious early convergence
		if (abs(qsimp-os) < EPS*abs(os) .or. &
			(qsimp == 0.0 .and. os == 0.0)) RETURN
	end if
	os=qsimp
	ost=os
end do
call nrerror('qsimp: too many steps')
END FUNCTION qsimp


FUNCTION qromb(func,a,b)
USE nrtype; USE nrutil, ONLY: nrerror
IMPLICIT NONE
REAL(SP),INTENT(IN):: a,b
REAL(SP):: qromb
INTERFACE
	FUNCTION func(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	END FUNCTION func
END INTERFACE
INTEGER(I4B),PARAMETER:: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
REAL(SP),PARAMETER:: EPS=1.0e-6_sp
! Returns the integral of the function func from a to b.
! Integration is performed by Romberg's method of order
! 2K, where, e.g., K=2 is Simpson's rule.
! Parameters: EPS is the fractional accuracy desired, as
! determined by the extrapolation error estimate; JMAX
! limits the total number of steps; K is the number of points
! used in the extrapolation.
REAL(SP),DIMENSION(JMAXP):: h,s		! These store the successive trapezoidal
REAL(SP):: dqromb					! approximations and their relative stepsizes.
INTEGER(I4B):: j
h(1)=1.0
do j=1,JMAX
	call trapzd(func,a,b,s(j),j)
	if (j >= K) then
		call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb,dqromb)
		if (abs(dqromb) <= EPS*abs(qromb)) RETURN
	end if
	s(j+1)=s(j)
	h(j+1)=0.25_sp*h(j)	! This is a key step: The factor is 0.25 even
end do					! though the stepsize is decreased by only 0.5
call nrerror('qromb: too many steps')	! This makes the extrapolation 
END FUNCTION qromb		! a polynomial is h^2 as allowed by equation (4.2.1)
						! not just a polynomial in h.


SUBROUTINE midpnt(func,a,b,s,n)
USE nrtype; USE nrutil, ONLY: arth
IMPLICIT NONE
REAL(SP),INTENT(IN):: a,b
REAL(SP),INTENT(INOUT):: s
INTEGER(I4B),INTENT(IN):: n
INTERFACE
	FUNCTION func(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	END FUNCTION func
END INTERFACE
! This routine computes the nth stage of refinement of an
! extended midpoint rule. func is input as the name of the 
! function to be integrated between limits a and b, also 
! input. When called with n=1, the routine as s the crudest
! estimate of int^b_a f(x)dx. Subsequent calls with n=2,3,..
! (in that sequential order) will improve the accuracy of s
! by adding (2/3)x3^(n-1) additional interior points. s should
! not be modified between sequential calls.
REAL(SP):: del
INTEGER(I4B):: it
REAL(SP),DIMENSION(2*3**(n-2)):: x
if (n == 1) then
	s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
else
	it=3**(n-2)
	del=(b-a)/(3.0_sp*it)	! The added points alternate in spacing between
	x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)	! del and 2*del.
	x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
	s=s/3.0_sp+del*sum(func(x))	! The new sum is combined with the old integral
end if							! to give a refined integral.
END SUBROUTINE midpnt


FUNCTION qromo(func,a,b,choose)
USE nrtype; USE nrutil, ONLY: nrerror
IMPLICIT NONE
REAL(SP),INTENT(IN):: a,b
REAL(SP):: qromo
INTERFACE
	FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	END FUNCTION func
	
	SUBROUTINE choose(funk,aa,bb,s,n)
	USE nrtype
	IMPLICIT NONE
	REAL(SP),INTENT(IN):: aa,bb
	REAL(SP),INTENT(INOUT):: s
	INTEGER(I4B),INTENT(IN):: n
	INTERFACE
		FUNCTION funk(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP),DIMENSION(:),INTENT(IN):: x
		REAL(SP),DIMENSION(size(x)):: funk
		END FUNCTION funk
	END INTERFACE
	END SUBROUTINE choose
END INTERFACE
INTEGER(I4B),PARAMETER:: JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1
REAL(SP),PARAMETER:: EPS=1.0e-6
! Romberg integration on an open interval. Returns the integral
! of the function func from a to b, using any specified integrating
! subroutine choose and Romberg's method. Normally choose will be
! an open formula, not evaluating the function at the endpoints. 
! It is assumed that choose triples the number of steps on each call,
! and that its error series contains only even powers of the number of
! steps. The routines midpnt, midinf, midsql, midsqu and midexp are
! possible choices for choose. The parameters have the same meaning as
! in qromb.
REAL(SP),DIMENSION(JMAXP):: h,s
REAL(SP):: dqromo
INTEGER(I4B):: j
h(1)=1.0
do j=1,JMAX
	call choose(func,a,b,s(j),j)
	if (j >= K) then
		call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromo,dqromo)
		if (abs(dqromo) <= EPS*abs(qromo)) RETURN
	end if
	s(j+1)=s(j)
	h(j+1)=h(j)/9.0_sp	! This is where the assumption of step
end do					! tripling and an even error series is used.
call nrerror('qromo: too many steps')
END FUNCTION qromo


SUBROUTINE midinf(funk,aa,bb,s,n)
USE nrtype; USE nrutil, ONLY: arth, assert
IMPLICIT NONE
REAL(SP),INTENT(IN):: aa,bb
REAL(SP),INTENT(INOUT):: s
INTEGER(I4B),INTENT(IN):: n
INTERFACE
	FUNCTION funk(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: funk
	END FUNCTION funk
END INTERFACE
! This routine is an exact replacement for midpnt, i.e.,
! returns as s the nth stage of refinement of the integral
! of funk from aa to bb, except that the function is 
! evaluated at evenly spaced points in 1/x rather than in
! x. This allows the upper limit bb to be as large and positive
! as the computer allows, or the lower limit aa to be as large
! and negative, but not both. aa and bb must have the same sign.
REAL(SP):: a,b,del
INTEGER(I4B):: it
REAL(SP),DIMENSION(2*3**(n-2)):: x
call assert(aa*bb > 0.0, 'midinf args')
b=1.0_sp/aa		! These two statements change the limits of 
a=1.0_sp/bb		! integration accordingly.
if (n==1) then	! From this point on, the routine is exactly
	s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))	! identical to midpnt.
else
	it=3**(n-2)
	del=(b-a)/(3.0_sp*it)
	x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
	x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
	s=s/3.0_sp+del*sum(func(x))
end if
CONTAINS
	FUNCTION func(x)	! This internal effects the change of variable.
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	func=funk(1.0_sp/x)/x**2
	END FUNCTION func
END SUBROUTINE midinf


SUBROUTINE midsql(funk,aa,bb,s,n)
USE nrtype; USE nrutil, ONLY: arth
IMPLICIT NONE
REAL(SP),INTENT(IN):: aa,bb
REAL(SP),INTENT(INOUT):: s
INTEGER(I4B),INTENT(IN):: n
INTERFACE
	FUNCTION funk(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: funk
	END FUNCTION funk
END INTERFACE
! This routine is an exact replacement for midpnt, i.e.,
! returns as s the nth stage of refinement of the integral
! of funk from aa to bb, except that it allows for an inverse
! square-root singularity in the integrand at the lower limit aa.
REAL(SP):: a,b,del
INTEGER(I4B):: it
REAL(SP),DIMENSION(2*3**(n-2)):: x
b=sqrt(bb-aa)	! These two statements change the limits of 
a=0.0			! integration accordingly.
if (n==1) then	! From this point on, the routine is exactly identical to midpnt.
	s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
else
	it=3**(n-2)
	del=(b-a)/(3.0_sp*it)
	x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
	x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
	s=s/3.0_sp+del*sum(func(x))
end if
CONTAINS
	FUNCTION func(x)	! This internal function effects the change of variable.
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	func=2.0_sp*x*funk(aa+x**2)
	END FUNCTION func
END SUBROUTINE midsql


SUBROUTINE midsqu(funk,aa,bb,s,n)
USE nrtype; USE nrutil, ONLY: arth
IMPLICIT NONE
REAL(SP),INTENT(IN):: aa,bb
REAL(SP),INTENT(INOUT):: s
INTEGER(I4B),INTENT(IN):: n
INTERFACE
	FUNCTION funk(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: funk
	END FUNCTION funk
END INTERFACE
! This routine is an exact replacement for midpnt, i.e.,
! returns as s the nth stage of refinement of the integral
! of funk from aa to bb, except that is allows for an inverse
! square-root singularity in the integrand at the upper limit bb.
REAL(SP):: a,b,del
INTEGER(I4B):: it
REAL(SP),DIMENSION(2*3**(n-2)):: x
b=sqrt(bb-aa)	! These two statements change the limits of integration
a=0.0			! accordingly
if (n==1) then	! From this point on, the routine is exactly identical to midpnt.
	s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
else
	it=3**(n-2)
	del=(b-a)/(3.0_sp*it)
	x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
	x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
	s=s/3.0_sp+del*sum(func(x))
end if
CONTAINS
	FUNCTION func(x)	! This internal function effects the change of variable.
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	func=2.0_sp*x*funk(bb-x**2)
	END FUNCTION func
END SUBROUTINE midsqu

SUBROUTINE midexp(funk,aa,bb,s,n)
USE nrtype; USE nrutil, ONLY: arth
IMPLICIT NONE
REAL(SP),INTENT(IN):: aa,bb
REAL(SP),INTENT(INOUT):: s
INTEGER(I4B),INTENT(IN):: n
INTERFACE
	FUNCTION funk(x)
	USE nrtype
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: funk
	END FUNCTION funk
END INTERFACE
! This routine is an exact replacement for midpnt, i.e.,
! returns as s the nth stage of refinement of the integral
! of funk from aa to bb, except that bb is assumed to be
! infinite (value passed not actually used). It is assumed
! that the function funk decreases exponentially rapidly at
! infinity.
REAL(SP):: a,b,del
INTEGER(I4B):: it
REAL(SP),DIMENSION(2*3**(n-2)):: x
b=exp(-aa)	! These two statements change the limits of
a=0.0		! integration accordingly.
if (n==1) then 	! From this point on, the routine is exactly identical to midpnt.
	s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
else
	it=3**(n-2)
	del=(b-a)/(3.0_sp*it)
	x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
	x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
	s=s/3.0_sp+del*sum(func(x))
end if
CONTAINS
	FUNCTION func(x)	! This internal function effects the change of variable.
	REAL(SP),DIMENSION(:),INTENT(IN):: x
	REAL(SP),DIMENSION(size(x)):: func
	func=funk(-log(x))/x
	END FUNCTION func
END SUBROUTINE midexp
