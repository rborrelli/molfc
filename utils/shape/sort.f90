module sort

    ! Implementa alcuni algoritmi di sorting
    ! bubble sort, heap sort e quicksort
    
    use param
    
    implicit none
    
    private
    
    public :: qsort2, heapsort
    
    integer, parameter :: NPAR_ARTH=16, NPAR2_ARTH=8
    
    interface swap
        module procedure swap_i
        module procedure swap_dp
    end interface swap
    
    interface arth
        module procedure arth_i
    end interface arth
    
    contains
    
    !==========================================================================!
	
    subroutine heapsort(arr)
    
    implicit none
    real(DP), dimension(:), intent(inout) :: arr
    integer :: i,n
    n=size(arr)
    do i=n/2,1,-1
        call sift_down(i,n)
    end do
    do i=n,2,-1
        call swap(arr(1),arr(i))
        call sift_down(1,i-1)
    end do
    contains
!BL
    subroutine sift_down(l,r)
    integer, intent(in) :: l,r
    integer :: j,jold
    real(DP) :: a
    a=arr(l)
    jold=l
    j=l+l
    do
        if (j > r) exit
        if (j < r) then
            if (arr(j) < arr(j+1)) j=j+1
        end if
        if (a >= arr(j)) exit
        arr(jold)=arr(j)
        jold=j
        j=j+j
    end do
    arr(jold)=a
    end subroutine sift_down    
    
    end subroutine heapsort
    
    !==========================================================================!
    
    subroutine swap_dp(a,b)
    real(DP), intent(inout) :: a,b
    real(DP) :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_dp

    subroutine swap_i(a,b)
    integer, intent(inout) :: a,b
    integer :: dum
    dum=a
    a=b
    b=dum
    end subroutine swap_i

    !==========================================================================!
    
    subroutine qsort2(arr,slave)
    
    implicit none
    real(DP), dimension(:), intent(inout) :: arr,slave
    integer :: i, alloc_err
    !INTEGER, DIMENSION(size(arr)) :: index  ! This causes problem on the stack
    integer, dimension(:), allocatable :: index
    real(DP), allocatable :: temp(:)
    
    allocate(index(1:size(arr)),stat=alloc_err)
    if (alloc_err /= 0) then
		print *, 'Cannot allocate INDEX array for reordering.'
		stop
	end if

    call indexx(arr,index)
    
    allocate(temp(1:size(arr)),stat=alloc_err)
    if (alloc_err /= 0) then
		print *, 'Cannot allocate temporary array for reordering.'
		stop
	end if
    ! the simple declaration arr = arr(index) can give a segfault signal if the array are
    ! too large. The solution below avoid this problem with an array allocated on the heap.
    ! The forall contruct is almost as fast as arr = arr(index).
    forall (i=1:size(arr)) temp(i) = arr(index(i))
    arr = temp
    forall (i=1:size(arr)) temp(i) = slave(index(i))
    slave = temp

    deallocate(temp,index)
    
    return    
    end subroutine qsort2
    
    !==========================================================================!
    
    subroutine indexx(arr,index)
    implicit none
    real(DP), dimension(:), intent(in) :: arr
    integer, dimension(size(arr)), intent(out) :: index
    integer, parameter :: NN=15, NSTACK=50
    real(DP) :: a
    integer :: n,k,k2,i,j,indext,jstack,l,r, temp
    integer, dimension(NSTACK) :: istack
    integer :: first, increment 

    n = size(arr)
!    index=arth(1,1,n)
!   The procedure arth_i has been inlined to avoid allocation of large temporary arrays.
!	the following is much less efficient but safer.
	forall (i=1:n) index(i) = i
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
                index(i+1)=indext
            end do
            if (jstack == 0) return
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
            if (jstack > NSTACK) then 
                print *, 'indexx: NSTACK too small'
            end if
            if (r-i+1 >= j-l) then
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
    contains
!BL
    subroutine icomp_xchg(i,j)
    integer, intent(inout) :: i,j
    integer :: swp
    if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
    end if
    end subroutine icomp_xchg
    end subroutine indexx
    
    function arth_i(first,increment,n)
    integer, intent(in) :: first,increment,n
    integer, dimension(n) :: arth_i
    integer :: k,k2,temp

    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
                arth_i(k)=arth_i(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
                arth_i(k)=arth_i(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
                if (k >= n) exit
                k2=k+k
                arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
                temp=temp+temp
                k=k2
        end do
    end if
    return
    end function arth_i
    
end module sort
