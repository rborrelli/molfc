module sysop

  use system_type
  use errors
  use util, only : isequal, isempty
  use matfun, only : eig
  use vecmath, only : cross, norm

  implicit none

  private

  public :: get_frequencies, get_geometry, set_geometry, get_atmass, get_symbol, &
            get_vibrations,  isLinear, get_atnum, get_bmat, get_normal_modes, &
            get_equilibrium_intc, get_geometry_vec, get_order, set_normal_modes, &
            read_intc, get_nat_bmat
            !set_equilibrium_intc,  get_equilibrium_intc, get_geometry_vec, get_order

  public :: get_state_from_id

  public :: operator(==)

  interface set_geometry
    module procedure set_geometry_mol
    module procedure set_geometry_state
  end interface set_geometry

  interface get_geometry
    module procedure get_geometry
    module procedure get_geometry_molecule
  end interface get_geometry

  interface get_geometry_vec
    module procedure get_geometry_molecule_vec
  end interface get_geometry_vec

  interface get_frequencies
    module procedure get_frequencies
    module procedure get_mol_frequencies
  end interface get_frequencies

  interface read_intc
    module procedure read_intc1
    module procedure read_intc2
  end interface read_intc

!  INTERFACE get_maxqn
!    MODULE PROCEDURE get_max_qn
!    MODULE PROCEDURE get_molecule_max_qn
!    MODULE PROCEDURE get_mol_max_qn
!  END INTERFACE get_maxqn

  interface get_atmass
    module procedure get_atmass
    module procedure get_molecule_atmass
  end interface get_atmass
  
!  INTERFACE set_excited_vibrations
!    MODULE PROCEDURE set_state_excited_vibrations
!    MODULE PROCEDURE set_mol_excited_vibrations
!  END INTERFACE set_excited_vibrations
  
  interface operator(==)
    module procedure mol_eq_mol
  end interface

!  INTERFACE write
!    MODULE PROCEDURE write_state
!  END INTERFACE write
 
  contains

!==============================================================================!

    function get_frequencies(estate) result(freq)

    type(state_t), intent(in) :: estate
    real(kind=dp) freq(1:estate%nvib)
    integer i, j, k

    k = 0
    do j = 1, estate%molecule%nvib
      k = k + 1
      freq(k) = estate%molecule%normodes%vibration(j)%freq
    end do

    return
    end function get_frequencies

!==============================================================================!

    function get_mol_frequencies(molecule) result(freq)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) freq(1:molecule%nvib)
    integer j

    do j = 1, molecule%nvib
      freq(j) = molecule%normodes%vibration(j)%freq
    end do

    return
    end function get_mol_frequencies

!==============================================================================!

    function get_geometry(estate) result(geometry)

    type(state_t), intent(in) :: estate
    real(kind=dp) geometry(1:3,1:estate%molecule%structure%numat)
    integer i
 
    do i = 1, estate%molecule%structure%numat
      geometry(1:3,i) = estate%molecule%structure%atom(i)%coord(1:3)
    end do
    
    return
    end function get_geometry

!==============================================================================!

    function get_geometry_molecule(molecule) result(geometry)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) geometry(1:3,1:molecule%structure%numat)
    integer i
 
    do i = 1, molecule%structure%numat
      geometry(1:3,i) = molecule%structure%atom(i)%coord(1:3)
    end do
    
    return
    end function get_geometry_molecule
    
!==============================================================================!

    function get_geometry_molecule_vec(molecule) result(geometry)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) geometry(1:3*molecule%structure%numat)
    integer i
 
    do i = 1, molecule%structure%numat
      geometry(3*(i-1)+1:3*i) = molecule%structure%atom(i)%coord(1:3)
    end do
    
    return
    end function get_geometry_molecule_vec
    
!==============================================================================!
    
    function get_vibrations(estate) result(vibrations)

    type(state_t), intent(in) :: estate
    real(kind=dp) vibrations(1:3*estate%molecule%structure%numat,1:estate%molecule%nvib)
    integer i, nm
 
    do nm = 1, estate%molecule%nvib
      do i = 1, estate%molecule%structure%numat
        vibrations(3*(i-1)+1:3*i,nm) = estate%molecule%normodes%vibration(nm)%atom(i)%d(1:3)
      end do
    end do
    
    return
    end function get_vibrations

!==============================================================================!

    function get_atmass(estate) result(atmass)

    type(state_t), intent(in) :: estate
    real(kind=dp) atmass(1:estate%molecule%structure%numat)
    integer i
 
    do i = 1, estate%molecule%structure%numat
      atmass(i) = estate%molecule%structure%atom(i)%elem%AM
    end do
    
    return
    end function get_atmass

!==============================================================================!

    function get_molecule_atmass(molecule) result(atmass)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) atmass(1:molecule%structure%numat)
    integer i
 
    do i = 1, molecule%structure%numat
      atmass(i) = molecule%structure%atom(i)%elem%AM
    end do

    return
    end function get_molecule_atmass
    
!==============================================================================!

    function get_symbol(estate) result(symbol)

    type(state_t), intent(in) :: estate
    character(len=3) symbol(1:estate%molecule%structure%numat)
    integer i
 
    do i = 1, estate%molecule%structure%numat
      symbol(i) = estate%molecule%structure%atom(i)%elem%sym
    end do
    
    return
    end function get_symbol

!==============================================================================!

    function get_atnum(estate) result(atnum)

    type(state_t), intent(in) :: estate
    integer atnum(1:estate%molecule%structure%numat)
    integer i
 
    do i = 1, estate%molecule%structure%numat
      atnum(i) = estate%molecule%structure%atom(i)%elem%Z
    end do
    
    return
    end function get_atnum

!==============================================================================!

    function get_bmat(molecule) result(bmat)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) bmat(1:molecule%intcoord%nintc,1:3*molecule%structure%numat)
    integer i, j, k
 
    bmat = zero
    do i = 1, molecule%intcoord%nintc
      do j = 1, size(molecule%intcoord%coord(i)%bmat(:,:),dim=1)
        k = molecule%intcoord%coord(i)%list(j)
        bmat(i,3*(k-1)+1:3*k) = molecule%intcoord%coord(i)%bmat(j,1:3)
      end do
    end do

    return
    end function get_bmat

!==============================================================================!

    function get_nat_bmat(molecule) result(bmat)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) bmat(1:molecule%intcoord%nintc,1:3*molecule%structure%numat)
    real(kind=dp) :: norm, ncoeff
    integer i, j, k, ij
 
    bmat = zero
    do i = 1, size(molecule%intcoord%ncoord)
      norm = 1.0
      !if (molecule%intcoord%ncoord(i)%type == "d") then
      !  norm = size(molecule%intcoord%ncoord(i)%coord)
      !else
      !  norm = sqrt(sum(molecule%intcoord%ncoord(i)%coord(:)%c**2)) !
      !end if
      !print *, 'norm coord ', i, norm, molecule%intcoord%ncoord(i)%type
      do ij = 1, size(molecule%intcoord%ncoord(i)%coord)
        !print *, 'c: ',molecule%intcoord%ncoord(i)%coord(ij)%c
        do j = 1, size(molecule%intcoord%ncoord(i)%coord(ij)%bmat(:,:),dim=1)
          k = molecule%intcoord%ncoord(i)%coord(ij)%list(j)
          ncoeff = molecule%intcoord%ncoord(i)%coord(ij)%c/norm 
          !print *, 'ncoeff', i, j, ncoeff
          bmat(i,3*(k-1)+1:3*k) = bmat(i,3*(k-1)+1:3*k) + &
                                  ncoeff*molecule%intcoord%ncoord(i)%coord(ij)%bmat(j,1:3)
        end do
      end do
    end do

    return
    end function get_nat_bmat

!==============================================================================!

    function get_normal_modes(molecule) result(t)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) t(1:3*molecule%structure%numat,1:molecule%nvib)
    integer i, j, k
 
    !---------------------------------------------------------!
    ! Ritorna una matrice 3N x 3N-6 le cui colonne contengono !
    ! i modi normali della molecola specificata.              !
    !---------------------------------------------------------!
    t = zero

    do j = 1, molecule%nvib
      do k = 1, molecule%structure%numat
        t(3*(k-1)+1:3*k,j) = molecule%normodes%vibration(j)%atom(k)%d(1:3)
      end do
    end do

    return
    end function get_normal_modes

!==============================================================================!

    subroutine set_normal_modes(molecule,t)

    type(molecule_t), intent(inout) :: molecule
    real(kind=dp), intent(in) :: t(1:3*molecule%structure%numat,1:molecule%nvib)
    integer i, j, k
 
    !---------------------------------------------------------!
    ! Assegna i modi normali a partire da una matrice         !
	! 3N x 3N-6 le cui colonne contengono                     !
    ! i modi normali della molecola specificata.              !
    !---------------------------------------------------------!

    do j = 1, molecule%nvib
      do k = 1, molecule%structure%numat
        molecule%normodes%vibration(j)%atom(k)%d(1:3) = t(3*(k-1)+1:3*k,j)  
      end do
    end do

    return
    end subroutine set_normal_modes

!==============================================================================!

    function get_state_from_id(system, stateid) result(id)
    
    type(system_t), intent(in) :: system
    character(len=*), intent(in) :: stateid
    integer :: id
    integer i, ierr
    character(len=80) errmesg

    id = 0
    do i = 1, system%nstate
        if (isequal(system%state(i)%id,stateid)) id = i
    end do

    if (id == 0) then 
      errmesg = "Cannot find state: "//stateid//" in your system."
      !do i = 1, system%nstate
      !   errmesg = errmesg//" Found state: "//system%state(i)%id//"."
      !end do
      ierr = error(0,errmesg)
    end if

    end function get_state_from_id 

!==============================================================================!

    function get_equilibrium_intc(molecule) result(q)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp) q(1:molecule%intcoord%nintc)
    integer nico
 
    nico = molecule%intcoord%nintc
    q(1:nico) = molecule%intcoord%coord(1:nico)%val

    end function get_equilibrium_intc

!==============================================================================!
    
!    SUBROUTINE set_equilibrium_intc(molecule)
!
!    type(molecule_t), intent(in out) :: molecule
!    real(kind=dp) intc(1:molecule%intcoord%nintc)
!    real(kind=dp) b1(3), b2(3), b3(3), b4(3)
!    real(kind=dp) n1(3), n2(3), phi, cosang
!    integer i, ierr
!    integer, pointer :: list(:)
! 
!    do i = 1, molecule%intcoord%nintc
!      list => molecule%intcoord%coord(i)%list 
!      select case(molecule%intcoord%coord(i)%type) 
!        case("s","S")
!          b1 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(1))%coord(:)  
!          molecule%intcoord%coord(i)%val = sqrt(dot_product(b1,b1))
!        case("b","B")
!          b1 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(1))%coord(:)  
!          b1 = b1/sqrt(dot_product(b1,b1)) 
!          b2 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(3))%coord(:)
!          b2 = b2/sqrt(dot_product(b2,b2)) 
!          molecule%intcoord%coord(i)%val = acos(dot_product(b1,b2))
!        case("d","D")
!          !if (molecule%intcoord%coord(i)%assigned) cycle
!          b1 = molecule%structure%atom(list(1))%coord(:) - molecule%structure%atom(list(2))%coord(:)  
!          b1 = b1/norm(b1)
!          b2 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(2))%coord(:)
!          b2 = b2/norm(b2)
!          b3 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(4))%coord(:)
!          b3 = b3/norm(b3)
!          n1 = cross(b1,b2) ! normale al piano 1 2 3
!          n1 = n1/norm(n1)
!          n2 = cross(b2,b3) ! normale al piano 2 3 4
!          n2 = n2/norm(n2)
!          phi = atan2 (dot_product(b2,cross(n1,n2)), dot_product(n1,n2)) 
!          ! ATAN2 restituisce phi in [-pi,pi], io lo rimetto in [0,2pi] 
!          ! mi serve per definire il verso di persorrenza della rotazione: vedi itransf.f90
!          if (phi < (-epsilon(one))) phi = phi + 2*PI
!          molecule%intcoord%coord(i)%val = phi
!        case("l","L")
!          molecule%intcoord%coord(i)%val = zero
!        case("w","W")
!          !if (molecule%intcoord%coord(i)%assigned) cycle
!!          write(fout,*)"Out-of-plane coordinate value read from input"
!          b1 = molecule%structure%atom(list(1))%coord(:) - molecule%structure%atom(list(4))%coord(:)  
!          b1 = b1/norm(b1)
!          b2 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(4))%coord(:)
!          b2 = b2/norm(b2)
!          b3 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(4))%coord(:)  
!          b3 = b3/norm(b3)
!          n2 = cross(b2,b3) ! normale al piano 2 3 4
!          n2 = n2/norm(n2)  ! nu = sin(phi1), see WDC pag 59.
!          molecule%intcoord%coord(i)%val = asin(dot_product(b1,n2))  ! Devo controllare se il segno va bene.
!        case("t","T")
!          write(fout,*)"Torsion coordinate value read from input"
!      end select
!    end do
!    
!    return
!    END SUBROUTINE set_equilibrium_intc
    
!==============================================================================!

    subroutine set_geometry_state(estate,coord)
    
    type(state_t), intent(inout) :: estate
    real(kind=dp), intent(in) :: coord(:,:)
    integer i
 
    do i = 1, estate%molecule%structure%numat
      estate%molecule%structure%atom(i)%coord(1:3) = coord(1:3,i)
    end do

    return
    end subroutine set_geometry_state

!==============================================================================!

    subroutine set_geometry_mol(molecule,coord)
    
    type(molecule_t), intent(inout) :: molecule
    real(kind=dp), intent(in) :: coord(:,:)
    integer i
 
    do i = 1, molecule%structure%numat
      molecule%structure%atom(i)%coord(1:3) = coord(1:3,i)
    end do

    return
    end subroutine set_geometry_mol

!==============================================================================!

    function isLinear(molecule) result(linear)

    type(molecule_t), intent(in) :: molecule
    real(kind=dp), allocatable :: crd(:,:), ti0(:,:), IABC(:), wm(:)
    real(kind=dp) tmd0, cdm
    integer i, j, k
    logical linear
    
    allocate(crd(1:3,1:molecule%structure%numat), IABC(1:3), ti0(1:3,1:3), &
             wm(1:molecule%structure%numat)) 
    
    linear = .false.
    
    if (molecule%structure%numat == 2) then
      linear = .true.
      return
    end if

    wm = zero; crd = zero
    crd = get_geometry(molecule)
    wm = get_atmass(molecule)
    
!    
! Put the molecule in the mass center.
!
    do i = 1, 3
      cdm = sum(wm(:)*crd(i,:))/sum(wm(:))
      crd(i,:) = crd(i,:) - cdm
    end do
          
! Per stabilire se la molecola e' lineare calola il tensore
! di inerzia.
! ---
! --- Calcola i tensori di inerzia delle due strutture
! ---
    tmd0 = zero
    do i = 1,3
      do j = 1, molecule%structure%numat
        tmd0 = tmd0 + wm(j)*(crd(i,j)**2)
      end do
    end do
   
    ti0 = zero
    do i=1,3
      do j=1,i
        do k=1, molecule%structure%numat
          ti0(i,j) = ti0(i,j) - wm(k)*crd(i,k)*crd(j,k)
          ti0(j,i) = ti0(i,j)
        end do
      end do
      ti0(i,i) = tmd0 + ti0(i,i)
    end do
   
    ti0 = ti0 / nav

    IABC = eig(ti0)
    
    if (any(abs(IABC(1:3)) <= epslin)) linear = .true.
    
    deallocate(crd,wm,IABC,ti0)
    
    return
    end function isLinear

!==============================================================================!

    function mol_eq_mol (m1, m2)
    
    type(molecule_t), intent(in) :: m1, m2
    logical :: mol_eq_mol
    real(kind=dp), allocatable :: geo1(:,:), geo2(:,:), f1(:), f2(:)
    real(DP) :: eps = 1.0e-6

    mol_eq_mol = .false.

    if (m1%structure%numat == m1%structure%numat) then
        allocate(geo1(1:3,1:m1%structure%numat),geo2(1:3,1:m1%structure%numat))
        geo1 = get_geometry(m1); geo2 = get_geometry(m2)
        allocate(f1(1:m1%nvib),f2(1:m1%nvib))
        mol_eq_mol = all((geo1 - geo2) < eps) .and. all((f1 - f2) < eps)
    end if

    deallocate(f1,f2,geo1,geo2)

    end function

!==============================================================================!

    function get_order(system,ilock,is) result(order)

    type(system_t), intent(in) :: system
    integer, intent(in) :: is, ilock
    integer order(1:system%state(is)%molecule%structure%numat)
    integer i, ia, ja, nato, iwarn, l, ierr
    real(kind=dp), allocatable :: dr(:)
    real(kind=dp), parameter :: LONG = 10000.0
    integer, external :: idamin

	! Init order array
	order = 0
    write(fout,'(2x,A,i3)') 'Finding order of atoms in state ', is
                
    nato = system%state(is)%molecule%structure%numat
    allocate(dr(1:nato))

NA: do ia = 1, nato
        ! Prende un atomo dello stato 1 e calcola la distanza da tutti quelli dello
        ! stesso tipo degli altri stati elettronici.
        dr = LONG
        do ja = 1, nato
			! skip atoms if already ordered
			if (any(ja == order(1:ia-1))) cycle
			if (system%state(ilock)%molecule%structure%atom(ia)%elem%Z /= &
				system%state(is)%molecule%structure%atom(ja)%elem%Z ) cycle
			dr(ja) = ZERO
			do l = 1, 3
				dr(ja) = dr(ja) + (system%state(ilock)%molecule%structure%atom(ia)%coord(l) - &
                         system%state(is)%molecule%structure%atom(ja)%coord(l))**2 
			end do
			dr(ja) = abs(sqrt(dr(ja)))
        end do
        !------------------------------------------------------------------------------------------
        ! Controlla se gli atomi sono troppo lontani tra loro. Se la distanza e' maggiore di 1 Ang 
        ! c'e' un warning.
        !------------------------------------------------------------------------------------------
        !if (control%printlevel >= 3) then
        if (minval(dr) > 1.0_dp) then
			write(fout,*) system%state(is)%molecule%id
			iwarn = warning(1)
			write(fout,'(a14,1x,a2,a,i2,a,x,a3,x,a2,a,i2,a,a3,f5.2)')'Distance ', &
				system%state(is-1)%molecule%structure%atom(ia)%elem%Sym,'(',ia,')',' - ', &
				system%state(is)%molecule%structure%atom(idamin(nato,dr,1))%elem%Sym,'(', &
				idamin(nato,dr,1),')', ' = ', dr(idamin(nato,dr,1))
        end if
        !end if

        ! Nuovo ordine.
        order(ia) = idamin(nato,dr,1)
    end do NA

	do i = 1, nato
		if (any(order(1:i-1) == order(i)) .or. any(order(i+1:nato) == order(i))) & 
			ierr = error(0,"Cannot find a proper ordering of the atoms.")
	end do

    end function get_order

	!=====================================================================!

	subroutine read_intc1(molecule,filename)

	! Thi subroutine read a file named INTVAL, placed into the working directory
	! which contains the equilibrium values of the internal coordinates.
	! The file INTVAL is structured as follows:
	! LINE 1: MOLECULE (1 or 2)
	! LINE 2: Numer of the coordinate, Its value
	! Example:
	! 1
	! 10 1.4
	! This set the internal coordinate 10 of molecule 1 equal to 1.4 Ang.
	
    type(molecule_t), intent(inout) :: molecule
	character(len=80), intent(in) :: filename
	integer, parameter :: ivc = 88
	integer :: i, istate, imol, icoo, ierr, ios
	real(kind=dp) :: rval
	character(len=80) :: line

	open(ivc,file=adjustl(filename),status='old')

	do	
		read(ivc,'(a)',end=10) line
		if (isempty(line)) exit
		read(line,*,iostat=ios) icoo, rval
		if (ios /= 0) ierr=error(0,"Cannot read internal coordinates from INTVAL file.")
		select case(molecule%intcoord%coord(icoo)%type)
			case("b","B","w","W","l","L","d","D","t","T")
				molecule%intcoord%coord(icoo)%val = (rval/180.0)*PI
			case("s","S")
				molecule%intcoord%coord(icoo)%val = rval
		end select
	end do

10  return
	end subroutine read_intc1

	!=====================================================================!

	subroutine read_intc2(molecule1,molecule2,filename)

	! Thi subroutine read a file named INTVAL, placed into the working directory
	! which contains the equilibrium values of the internal coordinates.
	! Each line of the file INTVAL is structured as follows:
	! MOLECULE (1 or 2),  Numer of the coordinate, value
	! Example:
	! 1 10 1.4
	! This set the internal coordinate 10 of molecule1 equal to 1.4 Ang.
	
    type(molecule_t), intent(inout) :: molecule1, molecule2
	character(len=80), intent(in) :: filename
	integer, parameter :: ivc = 88
	integer :: i, istate, imol, icoo, ierr, ios
	real(kind=dp) :: rval
	character(len=80) :: line

	open(ivc,file=adjustl(filename),status='old')

        print *, 'reading internal coordinates from external file'
	do	
		read(ivc,'(a)',end=10) line
                print *, line
		if (isempty(line)) exit
		read(line,*,iostat=ios) istate, icoo, rval
                print *, 'parse ', istate, icoo, rval
		if (ios /= 0) ierr=error(0,"Cannot read internal coordinates from INTVAL file.")
		print *, 'coordinate type',molecule1%intcoord%coord(icoo)%type
		select case(molecule1%intcoord%coord(icoo)%type)
			case("b","B","w","W","l","L","d","D","t","T")
                                print *, 'b '
                                if (istate == 1) then 
				  molecule1%intcoord%coord(icoo)%val = (rval/180.0)*PI
                                else if (istate == 2) then
				  molecule2%intcoord%coord(icoo)%val = (rval/180.0)*PI
                                end if
			case("s","S")
                                print *, 's '
                                if (istate == 1) then 
				  molecule1%intcoord%coord(icoo)%val = rval
                                else if (istate == 2) then
				  molecule2%intcoord%coord(icoo)%val = rval
                                end if
		end select
	end do

10  return
	end subroutine read_intc2

end module sysop
