module matfun

  use parameters
  use errors
  
  implicit none

  private

  public :: expm, sqrtm, det, eig, eigv, invm, vecsort, kron, vsort, diag, zinvmv

  interface eigv
    module procedure eigv
    module procedure eigav
  end interface eigv

  interface vsort
    module procedure vsort_dp
    module procedure vsort_int
  end interface vsort
  
  interface det
    module procedure det
    module procedure zdet
  end interface det

  interface invm
    module procedure invm
    module procedure zinvm
  end interface invm

  interface diag
    module procedure diagm
    module procedure diagv
    module procedure zdiagv
    module procedure zdiagm
  end interface diag

	interface sqrtm
		module procedure zsqrtm
		module procedure sqrtm
	end interface sqrtm

  contains

  ! ------------------------------------------------------

  function sqrtm(A)

  ! Al termine la matrice A non viene distrutta.
  ! E' la stessa di quella iniziale.
  real(kind=dp), target :: A(:,:)
  real(kind=dp) sqrtm(1:size(A,1),1:size(A,1))
  real(kind=dp), allocatable :: eig(:)
  real(kind=dp), allocatable :: work(:)
  real(kind=dp) temp
  real(kind=dp), allocatable :: Ap(:,:)

  integer S1, S2, info, lwork, i, j ,k
  
  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) stop
  
  lwork = 3*S1 - 1
  allocate (work(1:lwork), eig(1:S1))        

  allocate(Ap(1:S1,1:S1))
  Ap = A
  call dsyev ('V', 'U', S1, Ap, S1, eig, work, lwork, info)

  sqrtm = zero
  do i = 1, S1
    do j = 1, S1
      do k = 1, S1
        temp = Ap(i,k)*sqrt(eig(k))
        sqrtm(i,j) = sqrtm(i,j) + temp*Ap(j,k)
      end do
    end do
  end do
  
  deallocate(work,eig,Ap)

  return
  end function sqrtm

  ! ------------------------------------------------------

  function zsqrtm(A)

  ! Al termine la matrice A non viene distrutta.
  ! E' la stessa di quella iniziale.
  complex(kind=dpc) :: A(:,:)
  complex(kind=dpc) zsqrtm(1:size(A,1),1:size(A,1))
  real(kind=dp), allocatable :: eig(:)
  real(kind=dp), allocatable :: rwork(:)
  complex(kind=dpc), allocatable :: work(:)
  complex(kind=dpc) :: temp
  complex(kind=dpc), allocatable :: Ap(:,:)

  integer S1, S2, info, lwork, i, j ,k

  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) stop

  lwork = 3*S1 - 1
  allocate (work(1:lwork), eig(1:S1))
  allocate (rwork(1:3*S1-2))

  allocate(Ap(1:S1,1:S1))
  Ap = A
  call zheev ('V', 'U', S1, Ap, S1, eig, work, lwork, rwork, info)

  zsqrtm = zero
  do i = 1, S1
    do j = 1, S1
      do k = 1, S1
        temp = Ap(i,k)*sqrt(eig(k))
        zsqrtm(i,j) = zsqrtm(i,j) + temp*Ap(j,k)
      end do
    end do
  end do

  deallocate(work,eig,Ap)

  return
  end function zsqrtm

  ! ------------------------------------------------------

  function expm(A)

  ! Al termine la matrice A non viene distrutta.
  ! E' la stessa di quella iniziale.
  real(kind=dp), target :: A(:,:)
  real(kind=dp) expm(1:size(A,1),1:size(A,1))
  real(kind=dp), allocatable :: eig(:)
  real(kind=dp), allocatable :: work(:)
  real(kind=dp) temp
  real(kind=dp), allocatable :: Ap(:,:)

  integer S1, S2, info, lwork, i, j ,k
  
  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) stop
  
  lwork = 3*S1 - 1
  allocate (work(1:lwork), eig(1:S1))        

  allocate(Ap(1:S1,1:S1))
  Ap = A
  call dsyev ('V', 'U', S1, Ap, S1, eig, work, lwork, info)

  expm = zero
  do i = 1, S1
    do j = 1, S1
      do k = 1, S1
        temp = Ap(i,k)*exp(eig(k))
        expm(i,j) = expm(i,j) + temp*Ap(j,k)
      end do
    end do
  end do
  
  deallocate(work,eig,Ap)

  return
  end function expm

  ! ------------------------------------------------------

  function invm(A)

  ! Invert a non singular square matrix.
  ! Al termine la matrice A non viene distrutta.
  ! E' la stessa di quella iniziale.
  real(kind=dp), intent(in) :: A(:,:)
  real(kind=dp) invm(1:size(A,1),1:size(A,1))
  real(kind=dp), allocatable :: Ap(:,:), work(:)
  integer, allocatable :: ipiv(:)

  integer ns, S2, info, lwork, ierr
  
  !print *, 'matrix A:'
  !print *, A, size(A,1)

  ns = size(A,1)
  S2 = size(A,2)

  if (ns /= S2) ierr = error(0,"Cannot invert a rectangular matrix.")
  
  lwork = max(1,2*ns)  ! TODO: calcolare le dimensioni ottimali per WORK.
  allocate (work(1:lwork))
  allocate (ipiv(1:ns))
  allocate (Ap(1:ns,1:ns))

  Ap = A ! Copia A in Ap: altrimenti A verrebbe distrutta dalla routine.
  call dgetrf (ns, ns, Ap, ns, ipiv, info)

  if (info /= 0) ierr = error(0,"Error in LU factorization. Cannot invert matrix.")

  call dgetri (ns, Ap, ns, ipiv, work, lwork, info)

  if (info /= 0) ierr = error(0,"Error in DGETRI(). Cannot invert matrix.")

  invm = Ap

  deallocate(work, ipiv, Ap)

  return
  end function invm

  ! ------------------------------------------------------
   
  function zinvm(A)

  ! Invert a non singular square matrix.
  ! Al termine la matrice A non viene distrutta.
  ! E' la stessa di quella iniziale.
  complex(kind=dpc), intent(in) :: A(:,:)
  complex(kind=dpc) zinvm(1:size(A,1),1:size(A,1))
  complex(kind=dpc), allocatable :: Ap(:,:), work(:)
  integer, allocatable :: ipiv(:)

  integer ns, S2, info, lwork, ierr

  ns = size(A,1)
  S2 = size(A,2)

  if (ns /= S2) ierr = error(0,"Cannot invert a rectangular matrix.")

  lwork = max(1,2*ns)  ! TODO: calcolare le dimensioni ottimali per WORK.
  allocate (work(1:lwork))
  allocate (ipiv(1:ns))
  allocate (Ap(1:ns,1:ns))

  Ap = A ! Copia A in Ap: altrimenti A verrebbe distrutta dalla routine.
  call zgetrf (ns, ns, Ap, ns, ipiv, info)

  if (info /= 0) ierr = error(0,"Error in ZGETRF (LU factorization). Cannot invert matrix.")

  call zgetri (ns, Ap, ns, ipiv, work, lwork, info)

  if (info /= 0) ierr = error(0,"Error in ZGETRI. Cannot invert matrix.")

  zinvm = Ap

  deallocate(work, ipiv, Ap)

  return
  end function zinvm

  ! ------------------------------------------------------

  function zinvmv(A,b) result(zvec)

  ! Compute the vector (A^-1)*b
  ! by solving the linear system: A*x = b
  ! Al termine la matrice A non viene distrutta
  ! e' la stessa di quella iniziale.
  complex(kind=dpc), intent(in) :: A(:,:), b(:)
  complex(kind=dpc) zvec(1:size(A,1))
  complex(kind=dpc), allocatable :: Ac(:,:), bc(:), work(:)
  real(kind=dp), allocatable :: anormv(:)
  real(kind=dp) :: anorm, rcond
  integer :: i
  integer, allocatable :: ipiv(:)

  integer ns, info, lwork, ierr

  ns = size(A,1)
  if (ns == 1) then
	  zvec = b/A(1,1)
	  return
  end if

  if (A(1,2) /= A(2,1)) ierr = error(0,'Error in ZINVMV. Matrix non symmetric')

  allocate (bc(1:ns))
  allocate (Ac(1:ns,1:ns),source=A)
  bc = b ! b cannot be destroyed!

  if (ns /= size(A,2)) ierr = error(0,"Error in ZINVMV. Cannot invert a rectangular matrix.")

  lwork = max(1,2*ns)  ! TODO: calcolare le dimensioni ottimali per WORK.
  allocate (work(1:lwork))
  allocate (ipiv(1:ns))

  ! A viene distrutta
  call zsytrf ('U',ns, Ac, ns, ipiv, work, lwork, info)

  if (info /= 0) ierr = error(0,"Error in ZSYTRF. Problem in LU factor.")

  allocate(anormv(1:ns))
  forall (i=1:ns) anormv(i) = sum(abs(A(:,i)))
  anorm = maxval(anormv)
  call zsycon ('U', ns, Ac, ns, ipiv, anorm, rcond, work, info)
  !print *, 'RCOND =', RCOND
  if (rcond < 1e-14) print *, 'BADLY CONDITIONED MATRIX, RCOND =', RCOND

  call zsytrs ('U', ns, 1, Ac, ns, ipiv, bc, ns, info)

  if (info /= 0) ierr = error(0,"Error in ZSYTRS. Cannot solve linear system.")

  zvec = bc

  deallocate(ipiv,work,bc,anormv)

  return
  end function zinvmv

  ! ------------------------------------------------------

  function det(A)
  
  !----------------------------------------------!
  ! Calcola il determinante di una matrice reale !
  !----------------------------------------------!

  real(kind=dp) det
  real(kind=dp) perm
  real(kind=dp) A(:,:)
  real(kind=dp), allocatable :: AC(:,:)
  integer, allocatable :: ipiv(:)

  integer nrow, ncol, info, cp, i, ierr
  
  nrow = size(A,1)
  ncol = size(A,2)

  if (nrow /= ncol) ierr = error(0,"Erron in DET. Cannot compute determinant of a rectangular matrix.")
  if (nrow == 1) then
	  det = A(1,1)
  end if
  
  allocate(AC(1:nrow,1:ncol))
  AC = A
  allocate (ipiv(1:min(nrow,ncol)))

  call dgetrf (nrow, ncol, AC, nrow, ipiv, info)

  if (info /= 0) ierr = error(0,"Error in DET(): LU factorization. Cannot compute determinant")

  ! note that ipiv(n) = n always!
  perm = one
  do i = 1, nrow-1
    if (ipiv(i) /= i) perm = -perm
  end do
  
  det = perm
  do i = 1, nrow
    det = det*AC(i,i)
  end do

  deallocate(ipiv, AC)

  end function det

  ! ------------------------------------------------------

  function zdet(A)

  !--------------------------------------------------!
  ! Calcola il determinante di una matrice complessa !
  !--------------------------------------------------!

  complex(kind=dpc) zdet
  real(kind=dp) perm
  complex(kind=dpc) A(:,:)
  complex(kind=dpc), allocatable :: AC(:,:)
  integer, allocatable :: ipiv(:)

  integer nrow, ncol, info, cp, i, ierr

  nrow = size(A,1)
  ncol = size(A,2)

  if (nrow /= ncol) ierr = error(0,"Erron in ZDET. Cannot compute determinant of a rectangular matrix.")
  if (nrow == 1) then
	  zdet = A(1,1)
	  return
  end if

  allocate(AC(1:nrow,1:ncol))
  AC = A
  allocate (ipiv(1:min(nrow,ncol)))

  call zgetrf (nrow, ncol, AC, nrow, ipiv, info)

  if (info /= 0) then
	  ierr = error(0,"Error in ZDET: LU factorization. Cannot compute determinant")
  end if

  ! note that ipiv(n) = n always!
  perm = one
  do i = 1, nrow-1
    if (ipiv(i) /= i) perm = -perm
  end do

  zdet = perm
  do i = 1, nrow
    zdet = zdet*AC(i,i)
  end do

  deallocate(ipiv, AC)

  end function zdet

  ! ------------------------------------------------------

  function eig(A)

  real(kind=dp) A(:,:)
  real(kind=dp) eig(1:size(A,1))
  real(kind=dp), allocatable :: work(:)

  integer S1, S2, info, lwork, NS
  
  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) stop
  NS = S1
  
  lwork = 3*NS - 1
  allocate (work(1:lwork))        

  call dsyev ('N', 'U', NS, A, NS, eig, work, lwork, info)

  deallocate(work)

  end function eig

  ! ------------------------------------------------------

  subroutine eigv(A,W)

  real(kind=dp), target :: A(:,:)
  real(kind=dp)  w(1:size(A,1))
  real(kind=dp), allocatable :: work(:)
  real(kind=dp), allocatable :: Ap(:,:)

  integer S1, S2, info, lwork, i, j ,k
  
  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) stop
  
  lwork = 3*S1 - 1
  allocate (work(1:lwork))        

  allocate(Ap(1:S1,1:S1))
  Ap = A
  call dsyev ('V', 'U', S1, Ap, S1, w, work, lwork, info)

  A = Ap
  return
  end subroutine eigv

  ! ------------------------------------------------------

  subroutine eigav(A,W,U)

  real(kind=dp), intent(in) :: A(:,:)
  real(kind=dp)  w(1:size(A,1)), U(1:size(A,1),1:size(A,1))
  real(kind=dp), allocatable :: work(:)
  real(kind=dp), allocatable :: Ap(:,:)

  integer S1, S2, info, lwork, i, j ,k
  
  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) stop
  
  lwork = 3*S1 - 1
  allocate (work(1:lwork))        

  allocate(Ap(1:S1,1:S1))
  Ap = A ! copio A in modo da non distruggerla dopo dyev.

  call dsyev ('V', 'U', S1, Ap, S1, w, work, lwork, info)

  U = Ap ! metto gli autovettori (che ora sono in Ap) in U.
  
  deallocate(Ap,work)
  return
  end subroutine eigav

  ! ------------------------------------------------------

  function kron(A,B)

  ! Kronecker product (tensor product) of two matrices

  real(kind=dp), intent(in) :: A(:,:), B(:,:)
  real(kind=dp)  :: kron(1:size(A,1)*size(B,1),1:size(A,2)*size(B,2))
  integer :: na, ma, nb, mb, i, j

  na = size(A,1); ma = size(A,2)
  nb = size(B,1); mb = size(B,2)

  forall (i=1:na, j=1:ma)
      kron(1+(i-1)*nb:i*nb,1+(j-1)*mb:j*mb) = A(i,j)*B
  end forall

  end function kron

  !--------------------------------------------------------

  function vecsort(vec) result(vsr)

  ! Ascending order sort. Bubble sort algorithm.

  integer vec(:)
  integer vsr(1:size(vec))

  integer i, j, vt
  
  vsr = vec

  do i = 1, size(vec) - 1
    do j = i+1, size(vec)
      if (vsr(i) > vsr(j)) then
        vt = vsr(i)
        vsr(i) = vsr(j)
        vsr(j) = vt
      end if
    end do
  end do

  return
  end function vecsort

  !----------------------------------------------------------------------------!
  
  subroutine vsort_dp(vec,ids,order)

  ! Ascending order sort. Bubble sort algorithm.

  real(kind=dp), intent(inout) :: vec(:)
  integer, optional, intent(out) :: ids(1:size(vec))
  character(len=1), optional :: order  
  character(len=1) :: actord
  
  integer i, j, id, JMAX
  real(kind=dp) :: vt
  logical sort, cond
  
  ! Initialize ids
  forall (i=1:size(vec)) ids(i) = i
  
!  print *, 'before', vec
  
  actord = "D"
  if (present(order)) actord = order

    JMAX =  size(vec) - 1
    do i = 1, size(vec) - 1
        sort = .false.
        do j = 1, JMAX
          if (actord == "D") cond = vec(j) < vec(j+1) ! sort with ascending order
          if (actord == "A") cond = vec(j) > vec(j+1) ! sort with descending order
            if (cond) then
                sort = .true.
                ! Swap vector elements
                vt = vec(j)
                vec(j) = vec(j+1)
                vec(j+1) = vt
                ! Store their position accordingly
                if (present(ids)) then
                    id = ids(j)
                    ids(j) = ids(j+1)
                    ids(j+1) = id
                end if
            end if
        end do
        if (.not.sort) exit
        JMAX = JMAX - 1
    end do  		
  
  return
  end subroutine vsort_dp

  !----------------------------------------------------------------------------!
  
  subroutine vsort_int(vec,ids,order)

  ! Ascending order sort. Bubble sort algorithm.

  integer, intent(inout) :: vec(:)
  integer, optional, intent(out) :: ids(1:size(vec))
  character(len=1), optional :: order  
  character(len=1) :: actord
  
  integer i, j, id, JMAX
  real(kind=dp) :: vt
  logical sort, cond
    
  ! Initialize ids
  forall (i=1:size(vec)) ids(i) = i
  
  actord = "D"
  if (present(order)) actord = order

          JMAX = size(vec) - 1		
		  do i = 1, size(vec) - 1
		    sort = .false.
		    do j = 1, JMAX
      	if (actord == "D") cond = vec(j) < vec(j+1) ! sort with descending order
      	if (actord == "A") cond = vec(j) > vec(j+1) ! sort with ascending order
        if (cond) then
		        sort = .true.
				! Swap vector elements
		        vt = vec(j)
		        vec(j) = vec(j+1)
		        vec(j+1) = vt
		        ! Store their position accordingly
		        if (present(ids)) then
			        id = ids(j)
			        ids(j) = ids(j+1)
			        ids(j+1) = id
		        end if
		      end if
		    end do
		    if (.not.sort) exit
		    JMAX = JMAX - 1
		  end do  		

  return
  end subroutine vsort_int

  !==================================================================!

  function diagv(A)

  ! extract the diagonal of a square matrix
  real(kind=dp) A(:,:)
  real(kind=dp) diagv(1:size(A,1))

  integer S1, S2, i, ierr

  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) ierr = error(0,"error in DIAG")

  forall (i=1:S1) diagv(i) = A(i,i)

  end function diagv

  !==================================================================!

  function diagm(A)

	! Create a square diagonal matrix from a vector
  real(kind=dp) A(:)
  real(kind=dp) diagm(1:size(A,1),1:size(A,1))

  integer S1, S2, i

  S1 = size(A)

  diagm = zero
  forall (i=1:S1) diagm(i,i) = A(i)

  end function diagm

  !==================================================================!

  function zdiagv(A)

  complex(kind=dpc) A(:,:)
  complex(kind=dpc) zdiagv(1:size(A,1))

  integer S1, S2, i

  S1 = size(A,1)
  S2 = size(A,2)

  if (S1 /= S2) stop

  forall (i=1:S1) zdiagv(i) = A(i,i)

  end function zdiagv

  !==================================================================!

  function zdiagm(A)

  complex(kind=dpc) A(:)
  complex(kind=dpc) zdiagm(1:size(A,1),1:size(A,1))

  integer S1, S2, i

  S1 = size(A)

  zdiagm = cmplx(zero,zero)
  forall (i=1:S1) zdiagm(i,i) = A(i)

  end function zdiagm

end module matfun
