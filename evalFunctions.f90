SUBROUTINE eulsum(sum,term,jterm)
USE nrtype; USE nrutil, ONLY: poly_term,reallocate
IMPLICIT NONE
REAL(SP),INTENT(INOUT):: sum
REAL(SP),INTENT(IN):: term
INTEGER(I4B),INTENT(IN):: jterm
! Incorporates into sum the jterm'th term, with value term, of an
! alternating series. sum is input as the previous partial sum, and
! is output as the new partial sum. The first call to this routine,
! with the first term in the series, should be with jterm=1. On the 
! second call, term should be set to the second term of the series,
! with sign opposite to that of the first call, and jterm should be
! 2. And so on.
REAL(SP),DIMENSION(:),POINTER,SAVE:: wksp
INTEGER(I4B),SAVE:: nterm		! Number of saved differences in wksp
LOGICAL(LGT),SAVE:: init=.true.
if (init) then		! Initialize.
	init=.false.
	nullify(wksp)
end if
if (jterm == 1) then
	nterm=1
	wksp => reallocate(wksp,100)
	wksp(1)=term
	sum=0.5_sp*term	! Return first estimate.
else	
	if (nterm+1 > size(wksp)) wksp => reallocate(wksp,2*size(wksp))
	wksp(2:nterm+1)=0.5_sp*wksp(1:nterm)	! Update saved quantities by van 
	wksp(1)=term							! Wijn-gaarden's algorithm
	wksp(1:nterm+1)=poly_term(wksp(1:nterm+1),0.5_sp)
	if (abs(wksp(nterm+1)) <= abs(wksp(nterm))) then ! Favorable to increase p
		sum=sum+0.5_sp*wksp(nterm+1)
		nterm=nterm+1						! and the table becomes longer
	else									! Favorable to increase n,
		sum=sum+wksp(nterm+1)				! the table doesn't becomes longer.
	end if
end if
END SUBROUTINE eulsum


SUBROUTINE ddpoly(c,x,pd)
USE nrtype; USE nrutil, ONLY: arth,cumprod,poly_term
IMPLICIT NONE
REAL(SP),INTENT(IN):: x
REAL(SP),DIMENSION(:),INTENT(IN):: c
REAL(SP),DIMENSION(:),INTENT(OUT):: pd
! Given the coefficients of a polynomial of degree N_c-1as an array
! c(1:N_c) with c(1) being the constant term, and given a value x,
! this routine returns the polynomial evaluated at x as pd(1) and
! N_d-1 derivatives as pd(2:N_d).
INTEGER(I4B):: i,nc,nd
REAL(SP),DIMENSION(size(pd)):: fac
REAL(SP),DIMENSION(size(c)):: d
nc=size(c)
nd=size(pd)
d(nc:1:-1)=poly_term(c(nc:1:-1),x)
do i=2,min(nd,nc)
	d(nc:i:-1)=poly_term(d(nc:i:-1),x)
end do
pd=d(1:nd)
fac=cumprod(arth(1.0_sp,1.0_sp,nd))	! After the first derivative, 
pd(3:nd)=fac(2:nd-1)*pd(3:nd)		! factorial constants come in
END SUBROUTINE ddpoly
