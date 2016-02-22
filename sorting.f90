SUBROUTINE sort_pick(arr)
USE nrtype
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
! Sort an array arr into ascending numerical order, by straight
! insertion. arr is replaced on output by its sorted rearrangement.
INTEGER(I4B):: i,j,n
REAL(SP):: a
n=size(arr)
do j=2,n		! Pick out each element in turn.
	a=arr(j)
	do i=j-1,1,-1	! Look for the place to insert it.
		if (arr(i) <= a) exit
		arr(i+1)=arr(i)
	end do
	arr(i+1)=a		! Insert it.
end do
END SUBROUTINE sort_pick


SUBROUTINE sort_shell(arr)
USE nrtype
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
! Sort an array arr into ascending numerical order by Shell's method
! (diminishing increment sort). arr is replaced on output by its
! sorted rearrangement.
INTEGER(I4B):: i,j,inc,n
REAL(SP):: v
n=size(arr)
inc=1
do				! Determine the starting increment.
	inc=3*inc+1
	if (inc > n) exit
end do
do				! Loop over the partial sorts.
	inc=inc/3
	do i=inc+1,n	! Outer loop of straight insertion.
		v=arr(i)
		j=i
		do		! Inner loop of straight insertion
			if (arr(j-inc) <= v) exit
			arr(j)=arr(j-inc)
			j=j-inc
			if (j <= inc) exit
		end do
		arr(j)=v
	end do
	if (inc <= 1) exit
end do
END SUBROUTINE sort_shell


SUBROUTINE sort_byreshape(arr)
USE nrtype; USE nrutil, ONLY: swap
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
! Sort an array arr by bubble sorting a succession of reshapings into
! array slices. The method is similar to Shell sort, but allows parallelization
! within the vectorized masked swap calls.
REAL(SP),DIMENSION(:,:),ALLOCATABLE:: tab
REAL(SP),PARAMETER:: big=huge(arr)
INTEGER(I4B):: inc,n,m
n=size(arr)
inc=1
do				! Find the largest increment that fits.
	inc=2*inc+1
	if (inc > n) exit
end do
do				! Loop over the different shapes for the reshaped array.
	inc=inc/2
	m=(n+inc-1)/inc
	allocate(tab(inc,m))	! Allocate space and reshape the array. big
	tab=reshape(arr,(/inc,m/),(/big/))	! ensures that fill elements stay
	do						! at the end.
		! Bubble sort all the rows in parallel.
		call swap(tab(:,1:m-1:2),tab(:,2:m:2),&
			tab(:,1:m-1:2)>tab(:,2:m:2))
		call swap(tab(:,2:m-1:2),tab(:,3:m:2),&
			tab(:,2:m-1:2)>tab(:,3:m:2))
		if (all(tab(:,1:m-1) <= tab(:,2:m))) exit
	end do
	arr=reshape(tab,shape(arr))	! Put the array back together for the next shape
	deallocate(tab)
	if (inc <= 1) exit
	end do
END SUBROUTINE sort_byreshape
		
		
SUBROUTINE sort(arr)
USE nrtype; USE nrutil, ONLY: swap,nrerror
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
INTEGER(I4B),PARAMETER:: NN=15,NSTACK=50
! Sorts an array arr into ascending numerical order using the 
! Quicksort algorithm. arr is replaced on output by its sorted
! rearrangement. Parameters: NN is the size of subarrays sorted
! by straight insertion and NSTACK is the required auxilliary storage.
REAL(SP):: a
INTEGER(I4B):: n,k,i,j,jstack,l,r
INTEGER(I4B),DIMENSION(NSTACK):: istack
n=size(arr)
jstack=0
l=1
r=n
do
	if (r-1 < NN) then	! Insertion sort when subarray small enough
		do j=l+1,r
			a=arr(j)
			do i=j-1,l,-1
				if (arr(i) <= a) exit
				arr(i+1)=arr(i)
			end do
			arr(i+1)=a
		end do
		if (jstack == 0) RETURN
		r=istack(jstack)	! Pop stack and begin a new round of partitioning.
		l=istack(jstack-1)
		jstack=jstack-2
	else					! Choose median of left, center and right elements
		k=(l+r)/2			! as partitioning element a. Also rearrange so that
		call swap(arr(k),arr(l+1))	! a(l) <= a(l+1) <= a(r)
		call swap(arr(l),arr(r),arr(l)>arr(r))
		call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
		call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
		i=l+1				! Initialize pointers for partitioning
		j=r
		a=arr(l+1)			! Partitioning element
		do					! Here is the meat
			do				! Scan up to find element >= a
				i=i+1
				if (arr(i) >= a) exit
			end do
			do				! Scan down to find element <= a
				j=j-1
				if (arr(j) <= a) exit
			end do
			if (j < i) exit	! Pointers crossed. Exit with partitioning complete.
			call swap(arr(i),arr(j))	! Exchange elements
		end do
		arr(l+1)=arr(j)		! Insert partitioning element
		arr(j)=a
		jstack=jstack+2
		! Push pointers to larger subarray on stack; process smaller subarray immediately
		if (jstack > NSTACK) call nrerror('sort: NSTACK too small')
		if (j-i+1 >= j-l) then
			istack(jstack)=r
			istack(jstack-1)=i
			r=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
	end do
END SUBROUTINE sort
	

SUBROUTINE sort2(arr,slave)
USE nrtype; USE nrutil, ONLY: assert_eq
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr,slave
! Sorts an array arr into ascending order using Quicksor, while
! making the corresponding rearrangement of the same-size array
! slave. The sorting and rearrangement are performed by means
! of an index array
INTEGER(I4B):: ndum
INTEGER(I4B),DIMENSION(size(arr)):: index
ndum=assert_eq(size(arr),size(slave),'sort2')
call indexx(arr,index)		! Make the index array
arr=arr(index)				! Sort arr.
slave=slave(index)			! Rearrange slave.
END SUBROUTINE sort2


RECURSIVE SUBROUTINE sort_bypack(arr)
USE nrtype; USE nrutil, ONLY: array_copy,swap
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
! Sort an array arr by recursively applying the Fortran 90 pack
! intrinsic. The method is similar to Quicksort, but this variant
! allows parallelization by the Fortran 90 compiler.
REAL(SP):: a
INTEGER(I4B):: n,k,nl,nerr
INTEGER(I4B),SAVE:: level=0
LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE:: mask
REAL(SP),DIMENSION(:),ALLOCATABLE,SAVE:: temp
n=size(arr)
if (n <= 1) RETURN
k=(1+n)/2
call swap(arr(1),arr(k),arr(1)>arr(k))	! Pivot element is median
call swap(arr(k),arr(n),arr(k)>arr(n))	! of first, middle and last
call swap(arr(1),arr(k),arr(1)>arr(k))
if (n <= 3) RETURN
level=level+1	! Keep track of recursion level to avoid allocation overhead
if (level == 1) allocate(mask(n),temp(n))
a=arr(k)
mask(1:n) = (arr <= a)		! Which elements move to left?
mask(k) = .false.
call array_copy(pack(arr,mask(1:n)),temp,nl,nerr)	! Move them.
mask(k) = .true.
temp(nl+2:n)=pack(arr,.not. mask(1:n))	! Move others to right
temp(nl+1)=a
arr=temp(1:n)
call sort_bypack(arr(1:nl))		! And recurse
call sort_bypack(arr(nl+2:n))
if (level == 1) deallocate(mask,temp)
level=level-1
END SUBROUTINE sort_bypack


SUBROUTINE sort_heap(arr)
USE nrtype; USE nrutil, ONLY: swap
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
! Sorts an array arr into ascending numerical order using the 
! Heapsort algorithm. arr is replaced on output by its sorted
! rearrangement
INTEGER(I4B):: i,n
n=size(arr)
do i=n/2,1,-1
	! The index i, which here determines the "left" range of the
	! sift-down, i.e., the element to be sifted down, is decremented
	! from n/2 down to 1 during the "hiring" (heap creation) phase.
	call sift_down(i,n)
end do
do i=n,2,-1
	! Here the "right" range of the sift-down is decremented from
	! n-1 down to 1 during the "retirement-and-promotion" (heap
	! selection) phase
	call swap(arr(1),arr(i))	! Call a space at the end of the array,
	call sift_down(1,i-1)		! and retire the top of the heap into it.
end do
CONTAINS

SUBROUTINE sift_down(l,r)
INTEGER(I4B),INTENT(IN):: l,r
! Carry out the sift-down on element arr(l) to maintain the heap structure.
INTEGER(I4B):: j,jold
REAL(SP):: a
a=arr(l)
jold=l
j=l+l
do		! "Do while j <= r."
	if (j > r) exit
	if (j < r) then
		if (arr(j) < arr(j+1)) j=j+1	! Compare to the better underling
	end if
	if (a >= arr(j)) exit	! Found a's level. Terminate the sift-down.
	arr(jold)=arr(j)		! Otherwise, demote a and continue
	jold=j
	j=j+j
end do
arr(jold)=a			! Put a into its slot
END SUBROUTINE sift_down
END SUBROUTINE sort_heap
	

SUBROUTINE sort_radix(arr)
USE nrtype; USE nrutil, ONLY: array_copy,nrerror
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
! Sort an array arr by radix sort on its bits.
INTEGER(I4B),DIMENSION(size(arr)):: narr,temp
LOGICAL,DIMENSION(size(arr)):: msk
INTEGER(I4B):: k,negm,ib,ia,n,nl,nerr
! Because we are going to transfer reals to integers, we must check
! that the number of bits is the same in each:
ib=bit_size(narr)
ia=ceiling(log(real(maxexponent(arr)-minexponent(arr),sp))/log(2.0_sp))&
	+ digits(arr)
if (ib /= ia) call nrerror('sort_radix: bit sizes not compatible')
negm=not(ishftc(1,-1))		! Mask for all bits except sign bit.
n=size(arr)
narr=transfer(arr,narr,n)
where (btest(narr,ib-1)) narr=ieor(narr,negm)	! Flip all bits on neg. numbers.
do k=0,ib-2
	! Work from low- to high-order bits, and partition the array
	! according to the value of the bit.
	msk=btest(narr,k)
	call array_copy(pack(narr,.not. msk),temp,nl,nerr)
	temp(nl+1:n)=pack(narr,msk)
	narr=temp
end do
msk=btest(narr,ib-1)		! The sign bit gets separate treatment,
call array_copy(pack(narr,msk),temp,nl,nerr)	! since here 1 comes
temp(nl+1:n)=pack(narr,.not. msk)	! before 0
narr=temp
where(btest(narr,ib-1)) narr=ieor(narr,negm)	! Unflip all bits on neg. numbers.
arr=transfer(narr,arr,n)
END SUBROUTINE sort_radix


SUBROUTINE indexx_sp(arr,index)
USE nrtype; USE nrutil, ONLY: arth,assert_eq,nrerror,swap
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(IN):: arr
INTEGER(I4B),DIMENSION(:),INTENT(OUT):: index
INTEGER(I4B),PARAMETER:: NN=15,NSTACK=50
! Indexes an array arr, i.e., outputs the array index of length N 
! such that arr(index(j)) is in ascending order for j=1,..,N.
! The input quantity arr is not changed.
REAL(SP):: a
INTEGER(I4B):: n,k,i,j,indext,jstack,l,r
INTEGER(I4B),DIMENSION(NSTACK):: istack
n=assert_eq(size(index),size(arr),'indexx_sp')
index=arth(1,1,n)
jstack=0
l=1
r=n
do
	if (r-l < NN) then
		do j=l+1,r
			indext=index(j)
			a=arr(indext)
			do i=j-1,l,-1
				if (arr(index(i)) <= a) exit
				index(i+1)=index(i)
			end do
			index(i+1)=index(i)
		end do
		if (jstack == 0) RETURN
		r=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+r)/2
		call swap(index(k),index(l+1))
		call icomp_xchg(index(l),index(r))
		call icom_xchg(index(l+1),index(r))
		call icom_xchg(index(l),index(l+1))
		i=l+1
		j=r
		indext=index(l+1)
		a=arr(indext)
		do
			do
				i=i+1
				if (arr(index(i)) >= a) exit
			end do
			do
				j=j-1
				if (arr(index(j)) <= a) exit
			end do
			if (j < i) exit
			call swap(index(i),index(j))
		end do
		index(l+1)=index(j)
		index(j)=indext
		jstack=jstack+2
		if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
		if (r-i+1 >= j-1) then
			istack(jstack)=r
			istack(jstack-1)=i
			r=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
end do
CONTAINS

SUBROUTINE icomp_xchg(i,j)
INTEGER(I4B),INTENT(INOUT):: i,j
INTEGER(I4B):: swp
if (arr(j) < arr(i)) then
	swp=i
	i=j
	j=swp
end if
END SUBROUTINE icomp_xchg
END SUBROUTINE indexx_sp

SUBROUTINE indexx_i4b(iarr,index)
USE nrtype; USE nrutil, ONLY: arth,assert_eq,nrerror,swap
IMPLICIT NONE
INTEGER(I4B),DIMENSION(:),INTENT(IN):: iarr
INTEGER(I4B),DIMENSION(:),INTENT(OUT):: index
INTEGER(I4B),PARAMETER:: NN=15,NSTACK=50
INTEGER(I4B):: a
INTEGER(I4B):: n,k,i,j,indext,jstack,l,r
INTEGER(I4B),DIMENSION(NSTACK):: istack
n=assert_eq(size(index),size(iarr),'indexx_sp')
index=arth(1,1,n)
jstack=0
l=1
r=n
do
	if (r-l < NN) then
		do j=l+1,r
			indext=index(j)
			a=iarr(indext)
			do i=j-1,1,-1
				if (iarr(index(i)) <= a) exit
				index(i+1)=index(i)
			end do
			index(i+1)=indext
		end do
		if (jstack == 0) RETURN
		r=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+r)/2
		call swap(index(k),index(l+1))
		call icomp_xchg(index(l),index(r))
		call icomp_xchg(index(l+1),index(r))
		call icomp_xchg(index(l),index(l+1))
		i=l+1
		j=r
		indext=index(l+1)
		a=iarr(indext)
		do
			do
				i=i+1
				if (iarr(index(i)) >= a) exit
			end do
			do
				j=j-1
				if (iarr(index(j)) <= a) exit
			end do
			if (j < i) exit
			call swap(index(i),index(j))
		end do
		index(l+1)=index(j)
		index(j)=indext
		jstack=jstack+2
		if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
		if (r-i+1 >= j-1) then
			istack(jstack)=r
			istack(jstack-1)=i
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
		end if
	end if
end do
CONTAINS

SUBROUTINE icomp_xchg(i,j)
INTEGER(I4B),INTENT(INOUT):: i,j
INTEGER(I4B):: swp
if (iarr(j) < iarr(i)) then
	swp=i
	i=j
	j=swp
end if
END SUBROUTINE icomp_xchg
END SUBROUTINE indexx_i4b


SUBROUTINE sort3(arr,slave1,slave2)
USE nrtype; USE nrutil, ONLY: assert_eq
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr,slave1,slave2
! Sorts an array arr into ascending order using Quicksor, while
! making the corresponding rearrangement of the same-size arrays
! slave1 and slave2. The sorting and rearrangement are performed
! by means of an index array.
INTEGER(I4B):: ndum
INTEGER(I4B),DIMENSION(size(arr)):: index
ndum=assert_eq(size(arr),size(slave1),size(slave2),'sort3')
call indexx(arr,index)		! Make the index array
arr=arr(index)				! Sort arr.
slave1=slave1(index)		! Rearrange slave1
slave2=slave2(index)		! and slave2
END SUBROUTINE sort3


FUNCTION rank(index)
USE nrtype; USE nrutil, ONLY: arth
IMPLICIT NONE
INTEGER(I4B),DIMENSION(:),INTENT(IN):: index
INTEGER(I4B),DIMENSION(size(index)):: rank
! Given index as output from the routine indexx, this routine
! returns a same-size array rank, the corresponding table of ranks
rank(index(:))=arth(1,1,size(index))
END FUNCTION rank


RECURSIVE SUBROUTINE index_bypack(arr,index,partial)
USE nrtype; USE nrutil, ONLY: array_copy,arth,assert_eq
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(IN):: arr
INTEGER(I4B),DIMENSION(:),INTENT(INOUT):: index
INTEGER,OPTIONAL,INTENT(IN):: partial
! Indexes an array arr, i.e., outputs the array index of length N
! such that arr(index(j)) is in ascending order for j=1,..,N. The
! method is to apply recursively the Fortran 90 pack intrinsic.
! This is similar to Quicksort, but allows parallelization by the
! Fortran 90 compiler. partial is an optional argument that is used
! only internally on the recursive calls.
REAL(SP):: a
INTEGER(I4B):: n,k,nl,indext,nerr
INTEGER(I4B),SAVE:: level=0
LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE:: mask
INTEGER(I4B),DIMENSION(:),ALLOCATABLE,SAVE:: temp
if (present(partial)) then
	n=size(index)
else
	n=assert_eq(size(index),size(arr),'indexx_bypack')
	index=arth(1,1,n)
end if
if (n <= 1) RETURN
k = (1+n)/2
call icomp_xchg(index(1),index(k))		! Pivot element is median
call icomp_xchg(index(k),index(n))		! of first, middle and last
call icomp_xchg(index(1),index(k))
if (n <= 3) RETURN
level=level+1	! Keep track of recursion level to avoid allocation overhead
if (level == 1) allocate(mask(n),temp(n))
indext=index(k)
a=arr(indext)
mask(1:n)=(arr(index) <= a)		! Which elements move to left?
mask(k)=.false.
call array_copy(pack(index,mask(1:n)),temp,nl,nerr)	! Move them.
mask(k) = .true.
temp(nl+2:n)=pack(index,.not. mask(1:n))	! Move others to right
temp(nl+1)=indext
index=temp(1:n)
call index_bypack(arr,index(1:nl),partial=1)	! And recurse
call index_bypack(arr,index(nl+2:n),partial=1)
if (level == 1) deallocate(mask,temp)
level=level-1
CONTAINS

SUBROUTINE icomp_xchg(i,j)
IMPLICIT NONE
INTEGER(I4B),INTENT(INOUT):: i,j
! Swap or don'r swap integer arguments, depending on the ordering
! of their corresponding elements in an array arr
INTEGER(I4B):: swp
if (arr(j) < arr(i)) then
	swp=i
	i=j
	j=swp
end if
END SUBROUTINE icomp_xchg
END SUBROUTINE index_bypack


FUNCTION select(k,arr)
USE nrtype; USE nrutil, ONLY: assert,swap
IMPLICIT NONE
INTEGER(I4B),INTENT(IN):: k
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
REAL(SP):: select
! Returns the kth smallest value in the array arr. The input array
! will be rearranged to have this value in location arr(k), with all
! smaller elements moved to arr(1:k-1) (in arbitrary order) and all
! larger elements in arr(k+1:) (also in arbitrary order).
INTEGER(I4B):: i,r,j,l,n
REAL(SP):: a
n=size(arr)
call assert(k >= 1, k <= n, 'select args')
l=1
r=n
do
	if (r-l <= 1) then	! Active partition contains 1 or 2 elements
		if (r-1 == 1) call swap(arr(l),arr(r),arr(l)>arr(r)) ! Active partition contains 2 elements
		select=arr(k)
		RETURN
	else
		i=(l+r)/2	! Choose median of left, center and right elements
		call swap(arr(i),arr(l+1))	! as partitioning element a. Also
		call swap(arr(l),arr(r),arr(l)>arr(r))	! rearrange so that arr(l)<=arr(l+1)<=arr(r)
		call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
		call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
		i=l+1		! Initialize pointers for partitioning
		j=r
		a=arr(l+1)	! Partitioning element.
		do			! Here is the meat
			do 		! Scan up to find element > a
				i=i+1
				if(arr(i) >= a) exit
			end do
			do		! Scan down to find element < a
				j=j-1
				if (arr(j) <= a) exit
			end do
			if (j < i) exit	! Pointers crossed. Exith with partitioning complete.
			call swap(arr(i),arr(j))	! Exchange elements.
		end do
		arr(l+1)=arr(j)		! Insert partitioning element
		arr(j)=a
		if (j >= k) r=j-1	! Keep active the partition that contains the kth element
		if (j <= k) l=i
	end if
end do
END FUNCTION select


FUNCTION select_inplace(k,arr)
USE nrtype
IMPLICIT NONE
INTEGER(I4B),INTENT(IN):: k
REAL(SP),DIMENSION(:),INTENT(IN):: arr
REAL(SP):: select_inplace
! Returns the kth smallest value in the array arr, without altering
! the input array. In Fortran 90's assumed memory-rich environment,
! we just call select in scratch space.
REAL(SP),DIMENSION(size(arr)):: tarr
tarr=arr
select_inplace=select(k,tarr)
END FUNCTION select_inplace


FUNCTION select_bypack(k,arr)
USE nrtype; USE nrutil, ONLY: array_copy,assert,swap
IMPLICIT NONE
INTEGER(I4B),INTENT(IN):: k
REAL(SP),DIMENSION(:),INTENT(INOUT):: arr
REAL(SP):: select_bypack
! Returns the kth smallest value in the array arr. The input array
! will be rearranged to have this value in location arr(k), with
! all smaller elements moved to arr(1:k-1) (in arbitrary order) and
! all larger elements in arr(k+1:) (also in arbitrary order). This
! implementation allows parallelization in the Fortran 90 pack intrinsic.
LOGICAL,DIMENSION(size(arr)):: mask
REAL(SP),DIMENSION(size(arr)):: temp
INTEGER(I4B):: i,r,j,l,n,nl,nerr
REAL(SP):: a
n=size(arr)
call assert(k >= 1, k <= n, 'select_bypack args')
l=1		! Initial left and right bounds
r=n
do		! Keep partitioning until desired element is found
	if (r-l <= 1) exit
	i=(l+r)/2
	call swap(arr(l),arr(i),arr(l)>arr(i))	! Pivot element is median
	call swap(arr(i),arr(r),arr(i)>arr(r))	! of first, middle and last
	call swap(arr(l),arr(i),arr(l)>arr(i))
	a=arr(i)
	mask(l:r)=(arr(l:r) <= a)
	mask(i) = .false.
	call array_copy(pack(arr(l:r),mask(l:r)),temp(l:n),nl,nerr) ! Move them
	j=l+nl
	mask(i)= .true.
	temp(j+1:r)=pack(arr(l:r),.not. mask(l:r))	! Move others to right
	temp(j)=a
	arr(l:r)=temp(l:r)
	if (k > j) then		! Reset bounds to whichever side has the desired element
		l=j+1
	else if (k < j) then
		r=j-1
	else
		l=j
		r=j
	end if
end do
if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r))	! Case of only two left
select_bypack=arr(k)
END FUNCTION select_bypack


SUBROUTINE select_heap(arr,heap)
USE nrtype; USE nrutil, ONLY: nrerror,swap
IMPLICIT NONE
REAL(SP),DIMENSION(:),INTENT(IN):: arr
REAL(SP),DIMENSION(:),INTENT(OUT):: heap
! Returns in heap, an array of length M, the largest M elements
! of the array of length N, with heap(1) guaranteed to be the
! Mth largest element. The array arr is not altered. For efficiency,
! this routine should be used only when M << N.
INTEGER(I4B):: i,j,k,m,n
m=size(heap)
n=size(arr)
if (m > n/2 .or. m < 1) call nrerror('probable misuse of select_heap')
heap=arr(1:m)
call sort(heap)		! Create innitial heap by overkill! We assume m << n
do i=m+1,n			! For each remaining element..
	if (arr(i) > heap(1)) then	! Put it on the heap?
		heap(1)=arr(i)
		j=1
		do			! Sift down
			k=2*j
			if (k > m) exit
			if (k /= m) then
				if (heap(k) > heap(k+1)) k=k+1
			end if
			if (heap(j) <= heap(k)) exit
			call swap(heap(k),heap(j))
			j=k
		end do
	end if
end do
END SUBROUTINE select_heap


FUNCTION eclass(lista,listb,n)
USE nrtype; USE nrutil, ONLY: arth,assert_eq
IMPLICIT NONE
INTEGER(I4B),DIMENSION(:),INTENT(IN):: lista,listb
INTEGER(I4B),INTENT(IN):: n
INTEGER(I4B),DIMENSION(n):: eclass
! Given M equivalences between pairs of n individual elements in
! the form of the input arrays lista and listb of length, this 
! routine returns in an array of length n the number of the 
! equivalence class of each of the n elements, integers between
! 1 and n (not all such integers used).
INTEGER:: j,k,l,m
m=assert_eq(size(lista),size(listb),'eclass')
eclass(1:n)=arth(1,1,n)		! initialize each element its own class
do l=1,m					! For each piece of input informatino..
	j=lista(l)
	do
		if (eclass(j) == j) exit
		j=eclass(j)
	end do
	k=listb(l)
	do			! Track second element up to its ancestor
		if (eclass(k) == k) exit
		k=eclass(k)
	end do
	if (j /= k) eclass(j)=k	! If they are not already related, make them so.
end do
do j=1,n					! Final sweep up to highest ancestors
	do
		if (eclass(j) == eclass(eclass(j))) exit
		eclass(j)=eclass(eclass(j))
	end do
end do
END FUNCTION eclass

FUNCTION eclazz(equiv,n)
USE nrtype; USE nrutil, ONLY: arth
IMPLICIT NONE
INTERFACE
	FUNCTION equiv(i,j)
	USE nrtype
	IMPLICIT NONE
	LOGICAL(LGT):: equiv
	INTEGER(I4B),INTENT(IN):: i,j
	END FUNCTION equiv
END INTERFACE
INTEGER(I4B),INTENT(IN):: n
INTEGER(I4B),DIMENSION(n):: eclazz
! Given a user-supplied logical function equiv that tells whether
! a pair of elements, each in the range 1..n, are related, return
! in an array of length n equivalence class numbers for each element.
INTEGER:: i,j
eclazz(1:n)=arth(1,1,n)
do i=2,n			! Loop over first element of all pairs.
	do j=1,i-1		! Loop over second element of all pairs.
		eclazz(j)=eclazz(eclazz(j))	! Sweep it up this much
		if (equiv(i,j)) eclazz(eclazz(eclazz(j)))=i
		! Good exercise for the reader to figure out why this much ancestry is necessary
	end do
end do
do i=1,n			! Only this much sweeping is needed finally
	eclazz(i)=eclazz(eclazz(i))
end do
END FUNCTION eclazz
