MODULE random

INTERFACE ran
	MODULE PROCEDURE ran
END INTERFACE

INTERFACE ran0
	MODULE PROCEDURE ran0_s,ran0_v
END INTERFACE

INTERFACE ran1
	MODULE PROCEDURE ran1_s,ran1_v
END INTERFACE

INTERFACE ran2
	MODULE PROCEDURE ran2_s,ran2_v
END INTERFACE

INTERFACE gasdev
	MODULE PROCEDURE gasdev_s,gasdev_v
END INTERFACE

CONTAINS

FUNCTION ran(idum)
IMPLICIT NONE
INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL :: ran
! "Minimal" random number generator of Park and Miller combined
! with a Marsaglia shift sequence. Returns a uniform random deviate
! between 0.0 and 1.0 (exclusive of the endpoint values). This fully
! portable, scalar generator has the "traditional" (not Fortran 90) calling
! sequence with a random deviate as the return function value: call with
! idum a negative integer to initialize; thereafter, do not alter idum except
! to reinitialize. The period of this generator is about 3.1x10^18
INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
REAL, SAVE :: am
INTEGER(K4B), SAVE :: ix=-1, iy=-1,k
if (idum <= 0 .or. iy < 0) then ! Initialize
	am=nearest(1.0,-1.0)/IM
	iy=ior(ieor(888889999,abs(idum)),1)
	ix=ieor(777755555,abs(idum))
	idum=abs(idum)+1	! Set idum positive.
end if
ix=ieor(ix,ishft(ix,13))	! Marsaglia shift sequence with period 2^32-1
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/IQ					! Park-Miller sequence by Schrage's mmethod. period 2^31-2
iy=IA*(iy-k*IQ)-IR*k
if (iy < 0) iy=iy+IM
ran=am*ior(iand(IM,ieor(ix,iy)),1)  ! Combine the two generators with masking to
END FUNCTION ran 					! ensure nonzero value.


SUBROUTINE ran0_s(harvest)
USE nrtype
USE ran_state, ONLY: amm,lenran,ran_init,iran0,jran0,kran0,nran0,rans
IMPLICIT NONE
REAL(SP),INTENT(OUT) ::  harvest
! Lagged Fibonacci generator combined with a Marsaglia shift sequence.
! Returns as harvest a uniform random deviate between 0.0 and 1.0 (exclusive
! of the endpoint values). This generator has the same calling and initialization
! conventions as Fortran 90's random_number routine. Use ran_seed to initialize
! or reinitialize to a particular sequence. The period of this geenrator is about
! 2.0x10^28, and it fully vectorizes. Validity of the integer model assumed by this
! generator is tested at initialization.
if (lenran < 1) call ran_init(1) ! Initialization routine in ran_state
rans = iran0-kran0 ! Update Fibonacci generator, which has period p^2+p+1,p=2^31-69
iran0=jran0
jran0=kran0
kran0=rans
nran0=ieor(nran0,ishft(nran0,13)) ! Update Marsaglia shift sequence with period 2^32-1
nran0=ieor(nran0,ishft(nran0,-17))
nran0=ieor(nran0,ishft(nran0,5))
rans=ieor(nran0,rans)
rans=ieor(nran0,rans)		! Combine the generators.
harvest=amm*merge(rans,not(rans),rans<0) ! Make the result positive define (not that amm is negative).
END SUBROUTINE ran0_s


SUBROUTINE ran0_v(harvest)
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,ran_init,iran,jran,kran,nran,ranv
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(OUT) ::  harvest
INTEGER(K4B) :: n
n=size(harvest)
if (lenran < n+1) call ran_init(n+1)
ranv(1:n)=iran(1:n)-kran(1:n)
where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
iran(1:n)=jran(1:n)
jran(1:n)=kran(1:n)
kran(1:n)=ranv(1:n)
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
ranv(1:n)=ieor(nran(1:n),ranv(1:n))
harvest=amm*merge(ranv(1:n),not(ranv(1:n)),ranv(1:n)<0)
END SUBROUTINE ran0_v


SUBROUTINE ran1_s(harvest)
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,ran_init,&
	iran0,jran0,kran0,nran0,mran0,rans
IMPLICIT NONE
REAL(SP),INTENT(OUT) :: harvest
! Lagged Fibonacci generator combined with two Marsglia shift sequence.
! On output, returns as harvest a uniform random deviate between 0.0 and
! 1.0 (exclusive of the endpoint values). This generator has the same calling
! and initialization conventions as Fortran 90's random_number routine. Use
! ran_seed to initialize or reinitialize to a particular sequence. The period
! of this generator is about 8.5x10^37, and it fully vectorizes. Validity of the 
! integer model assumed by this generator is tested at initialization.
if (lenran < 1) call ran_init(1) ! Initialization routine in ran_state.
rans=iran0-kran0	! Update Fibonacci generator, which has period p^2+p+1,p=2^31-69
iran0=jran0
jran0=kran0
kran0=rans
nran0=ieor(nran0,ishft(nran0,13)) ! Update Marsaglia shift sequence
nran0=ieor(nran0,ishft(nran0,-17))
nran0=ieor(nran0,ishft(nran0,5))
! Once only per cycle, advance sequence by 1, shortening its period 2^32-2
if (nran0 == 1) nran0=270369_k4b
mran0=ieor(mran0,ishft(mran0,5)) ! Update Marsaglia shift sequence with 2^32-1
mran0=ieor(mran0,ishft(mran0,-13))
mran0=ieor(mran0,ishft(mran0,6))
rans=ieor(nran0,rans)+mran0
! Combine the generators. The above statement has wrap-around addition
harvest=amm*merge(rans,not(rans),rans<0) ! Make the result positive definite(note that amm is negative)
END SUBROUTINE ran1_s


SUBROUTINE ran1_v(harvest)
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,&
	iran,jran,kran,nran,mran,ranv
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(OUT) :: harvest
INTEGER(K4B) :: n
n=size(harvest)
if(lenran<n+1) ranv(1:n)=ranv(1:n)+2147483579_k4b
iran(1:n)=jran(1:n)
jran(1:n)=kran(1:n)
kran(1:n)=ranv(1:n)
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
where (nran(1:n) == 1) nran(1:n)=270369_k4b
mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
harvest=amm*merge(ranv(1:n),not(ranv(1:n)),ranv(1:n)<0)
END SUBROUTINE ran1_v


SUBROUTINE ran2_s(harvest)
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,ran_init,&
	iran0,jran0,kran0,nran0,mran0,rans
IMPLICIT NONE
REAL(SP),INTENT(OUT) :: harvest
! Lagged Fibonacci generator combined with a Marsaglia shift
! sequence and a linear congruential generator. Returns a harvest
! a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
! values). This generator has the same calling and initalization 
! conventions as Fortran 90's random_number routine. Use ran_seed to
! initialize or reinitialize to a particular sequence. The period of
! this generator is about 8.5x10^37, and it fully vectorizes. Validity
! of the integer model assumed by this generator is tested at initialization.
if (lenran < 1) call ran_init(1) ! Initialization routine in ran_state
rans=iran0-kran0 	! Update Fibonacci generator, which has period p^2+p+1,p=2^31-69
iran0=jran0
jran0=kran0
kran0=rans
nran0=ieor(nran0,ishft(nran0,13))	! Update Marsaglia shift sequence with period 2^32-1
nran0=ieor(nran0,ishft(nran0,-17))
nran0=ieor(nran0,ishft(nran0,5))
rans=iand(mran0,65535)
! Update the sequence m <- 69069m+820265819 mod 2^32 using shifts unstead of multiplies.
! Wrap-around addition (tested at initialization) is used.
mran0=ishft(3533*ishft(mran0,-16)+rans,16)+&
	3533*rans+820265819_k4b
rans=ieor(nran0,kran0)+mran0 ! Combine the generators
harvest=amm*merge(rans,not(rans),rans<0) ! Make the result positive definite (note that amm is negative).
END SUBROUTINE ran2_s


SUBROUTINE ran2_v(harvest)
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,ran_init,&
	iran,jran,kran,nran,mran,ranv
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(OUT) :: harvest
INTEGER(K4B) :: n
n=size(harvest)
if (lenran < n+1) call ran_init(n+1)
ranv(1:n)=iran(1:n)-kran(1:n)
where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
iran(1:n)=jran(1:n)
jran(1:n)=kran(1:n)
kran(1:n)=ranv(1:n)
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
ranv(1:n)=iand(mran(1:n),65535)
mran(1:n)=ishft(3533*ishft(mran(1:n),-16)+ranv(1:n),16)+&
	3533*ranv(1:n)+820265819_k4b
ranv(1:n)=ieor(nran(1:n),kran(1:n))+mran(1:n)
harvest=amm*merge(ranv(1:n),not(ranv(1:n)),ranv(1:n)<0)
END SUBROUTINE ran2_v


SUBROUTINE gasdev_s(harvest)
USE nrtype
IMPLICIT NONE
REAL(SP),INTENT(OUT):: harvest
! Returns in harvest a normally distributed deviate with zero mean
! and unit variane, using ran1 as the source of uniform deviates.
REAL(SP):: rsq,v1,v2
REAL(SP),SAVE:: g
LOGICAL,SAVE:: gaus_stored=.false.
if (gaus_stored) then	! We have an extra deviate handy
	harvest=g			! so return it,
	gaus_stored=.false.	! and unset the flag
else					! We don't have an extra deviate handy, so
	do
		call ran1(v1)	! pick two uniform numbers in the square
		call ran1(v2)	! extending from -1 to +1 in each direction,
		v1=2.0_sp*v1-1.0_sp
		v2=2.0_sp*v2-1.0_sp
		rsq=v1**2+v2**2	! see if they are in unit circle,
		if (rsq > 0.0 .and. rsq < 1.0) exit
	end do				! otherwise try again.
	rsq=sqrt(-2.0_sp*log(rsq)/rsq)	! Noe make the Box-Muller transformation
	harvest=v1*rsq		! to get two normal deviates. Return one and 
	g=v2*rsq			! save the other for next time
	gaus_stored=.false.
end if
END SUBROUTINE gasdev_s

SUBROUTINE gasdev_v(harvest)
USE nrtype; USE nrutil, ONLY: array_copy
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(OUT):: harvest
REAL(SP),DIMENSION(size(harvest)):: rsq,v1,v2
REAL(SP),ALLOCATABLE,DIMENSION(:),SAVE:: g
INTEGER(I4B):: n,ng,nn,m
INTEGER(I4B),SAVE:: last_allocated=0
LOGICAL,SAVE:: gaus_stored=.false.
LOGICAL,DIMENSION(size(harvest)):: mask
n=size(harvest)
if (n /= last_allocated) then
	if (last_allocated /= 0) deallocate(g)
	allocate(g(n))
	last_allocated=n
	gaus_stored=.false.
end if
if (gaus_stored) then
	harvest=g
	gaus_stored=.false.
else 
	ng=1
	do
		if (ng > n) exit
		call ran1(v1(ng:n))
		call ran1(v2(ng:n))
		v1(ng:n)=2.0_sp*v1(ng:n)-1.0_sp
		v2(ng:n)=2.0_sp*v2(ng:n)-1.0_sp
		rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
		mask(ng:n)=(rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
		call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
		v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
		rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
		ng=ng+nn
	end do
	rsq=sqrt(-2.0_sp*log(rsq)/rsq)
	harvest=v1*rsq
	g=v2*rsq
	gaus_stored=.true.
end if
END SUBROUTINE gasdev_v


END MODULE random
