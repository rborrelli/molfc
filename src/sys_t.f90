module system_type

	use parameters
	use iomatrix

	implicit none

	integer, parameter, private :: flen = 80, plen = 80

	type, public :: file_t
		integer :: iu = 0 ! The logical unit of the fortran file (must be set!!! during the program execution)
		character(len = flen) :: name = BLANK
		character(len = plen) :: path = BLANK
		character(len = flen) :: type
	end type file_t

	type, public :: element
		real(kind = dp) AM ! Atomic mass
		integer N ! Number of neutrons
		integer Z ! Atomic charge
		integer AN ! Atomic number
		integer val ! Valence (maximum number of bonds)
		real(kind = sp) rad ! Atomic Radius (used to check connectivity) (Angstrom)
		character(len = SYMLEN) Sym ! Symbol
	end type element

	type, public :: atom_t
		real(kind = dp) coord(1:3)
		type(element) elem
	contains
		procedure :: bond
	end type atom_t

	type, public :: displacement_t
		real(kind = dp) d(1:3)
		type(element) elem
	end type displacement_t

	type, public :: vibration_t
		real(kind = dp) freq
		real(kind = dp) :: amp = zero, & ! --- Parameters for gaussian type coupling.
		pos = zero, &
		coeff = zero
		real(kind = dp) :: fcpl = zero, & ! --- First order counpling
		scpl = zero ! --- Second order coupling
		integer id
		integer :: nq = 0
		type(displacement_t), allocatable :: atom(:)
                character(len=3) :: symlab
	end type vibration_t

	type, public :: curvc_t
		character(len = 4) :: type
		real(kind = dp) :: c = one ! coefficient of the coordinate if it is a linear combination in a natural coordinate
		real(kind = dp) :: val = zero ! Value of the internal coordinates
		integer, pointer :: list(:) ! Indices of the atoms specifying internal coordinate
		real(kind = dp), allocatable :: bmat(:,:) ! Wilson B matrix for the coordinate.
	end type curvc_t

    ! A natural coordinate is a linear combination of curvilinear coordinates
    type, public :: natcoord_t
        character(len = 4) :: type
	real(kind = dp) :: val = zero ! Value of the natural internal coordinate
        !real(kind = dp), allocatable :: bmat(:,:) ! Wilson B matrix for the coordinate.
        type(curvc_t), allocatable :: coord(:)
    end type natcoord_t

	type, public :: intcoord_t
		integer :: nintc = 0 ! Number of internal coordinates.
		type(curvc_t), allocatable :: coord(:)
		type(natcoord_t), allocatable :: ncoord(:)
	end type intcoord_t

	type, public :: structure_t
		integer :: numat = 0
		logical :: linear = .false.
		logical :: fromfile = .false.
		real(kind = dp) :: axsgn(1:3) = 1.d0
		character(len = IDLEN) id
		character(len = 4) :: units = "ang", &
		coord = "xyz"
		type(atom_t), allocatable :: atom(:)
		type(file_t) :: file
		integer, allocatable :: order(:) ! Ordine degli atomi nella molecola.
	end type structure_t

	type normodes_t
		logical :: massw = .false.
		logical :: sort = .false. ! if .true. vibrations are sorted assordin to symmetry labels. Implemented only for g09hp and TM.
		logical :: fromfile = .false.
		real(kind = dp) :: scfreq = one
		real(kind=dp) :: fscale = one ! scaling factor for the vibrational frequencies.
		character(len = 4) :: units = "ang"
		type(vibration_t), allocatable :: vibration(:)
		type(file_t) :: file
		integer, allocatable :: order(:) ! Ordinamento delle vibrazioni
	end type normodes_t

	type, public :: molecule_t
		integer :: nvib = 0
		logical :: reorient = .false. ! Di default la moecola non viene riorientata.
		logical :: reord = .true. ! Di default la molecola viene riordinata.
		integer :: orient(1:3) = one ! Matrice di orientazione degli assi.
		type(structure_t) structure
		type(normodes_t) normodes
		type(intcoord_t) intcoord
		integer, allocatable :: acm(:,:) ! Atom connectivity matrix
		character(len = IDLEN) id
	contains
		procedure :: setacm
        procedure :: getSubset
	end type molecule_t

	type, public :: state_t
		real(kind = dp) :: energy = zero ! --- Electronic energy of the state.
		real(kind = dp) :: zcpl = zero ! --- Zero order coupling between electronic states.
		integer :: nmol = 0 ! Number of molecules.
		integer :: nvib = 0 ! Total number of normal modes (sum of normal modes of all molecules).
		type(molecule_t) :: molecule
		character(len = IDLEN) id
	end type state_t

    type, public :: tm_t
        real(kind=dp) :: mu0(1:3) = 0.0
        real(kind=dp), allocatable :: dmudQ(:,:)
        real(kind=dp), allocatable :: dmuds(:,:) ! this is the derivative with respect to dimensionless normal coordinates
    end type tm_t

	type, public :: system_t
		logical :: model = .false.
		integer :: nstate = 0
		type(state_t), allocatable :: state(:)
		type(tm_t) :: tm
	end type system_t

!	type, public :: beta_t
!		real(kind = dp) :: zcpl = zero
!		real(kind = dp), allocatable :: fcpl(:)
!	end type beta_t

contains

	!========================================================================= !

	subroutine setacm(molecule)
		! This subroutine sets the Atom Connectivity Matrix.
		class(molecule_t), intent(inout) :: molecule
		integer, allocatable :: acm(:,:)
		integer i, j, numat

		numat = molecule%structure%numat

		allocate(acm(1:molecule%structure%numat, 1:molecule%structure%numat))
		acm = 0
		do i = 1, molecule%structure%numat
			do j = i + 1, molecule%structure%numat
				!acm(i,j) = bond(molecule%structure%atom(i),molecule%structure%atom(j))
				acm(i, j) = molecule%structure%atom(i)%bond(molecule%structure%atom(j))
				acm(j, i) = acm(i, j)
			end do
		end do

		! The diagonal elemet define the cardinality of the atom, i.e. the
		! number of bonds that it forms with all the other atoms.
		do i = 1, molecule%structure%numat
			acm(i, i) = sum(acm(i,:))
		end do

		allocate(molecule%acm(1:molecule%structure%numat, 1:molecule%structure%numat))

		! Assign acm matrix to the molecule
		molecule%acm = acm

		if (.true.) then
			write(fout, '(//)')
			write(fout, '(2x,a)') ' -------------------------'
			write(fout, '(2x,a,2x,a)') ' Atom Connectivity Matrix:', molecule%id
			write(fout, '(2x,a)') ' -------------------------'
			call layout(acm, numat, numat, numat, numat)
		end if

		deallocate(acm)

		return
	end subroutine setacm

	!========================================================================= !

	integer function bond(atom_i, atom_j)

		! This function return 1 if the two atoms A, B, are bonded.
		! The criterion for the existence of a bond is the standard used in most
		! crystallographic programs. Two atoms are bonded if their distance (d)
		! is in the range (Rcov(A) + Rcov(B) +/- t) where t is a tolerance which is
		! is set to 0.4 according to the Cambridge Structural Database.

		class(atom_t), intent(in) :: atom_i, atom_j
		real(kind = dp) dist(3), distance
		real(kind = dp), parameter :: toler = 0.4

		bond = 0

		dist = atom_i%coord - atom_j%coord
		distance = sqrt(sum(dist**2))

		if (abs(distance - toler) <= (atom_i%elem%rad + atom_j%elem%rad)) bond = 1

		return
	end function

	!========================================================================= !

        subroutine getSubset(molecule,incvib)
                
          class(molecule_t), intent(in out) :: molecule
          type(molecule_t)  :: newMolecule 
          integer, intent(in) :: incvib(:) ! vibrations defining the subset
          integer :: i, j, alloc_err, numat
          
	      newMolecule%nvib = size(incvib)
          numat = molecule%structure%numat

          allocate(newMolecule%normodes%vibration(1:newMolecule%nvib))
          do i = 1, newMolecule%nvib
            allocate(newMolecule%normodes%vibration(i)%atom(1:numat),stat=alloc_err)
            ! This assign all the elements of the vibration
            newMolecule%normodes%vibration(i) = molecule%normodes%vibration(incvib(i))
            newMolecule%normodes%vibration(i)%atom(1:numat) = molecule%normodes%vibration(incvib(i))%atom(1:numat)
          end do

          ! deallocate atom arrays
          do i=1,newMolecule%nvib
            deallocate(molecule%normodes%vibration(i)%atom)
          end do
          ! deallocate vibration array
          deallocate(molecule%normodes%vibration)
          ! re-allocate vibration array
          allocate(molecule%normodes%vibration(1:newMolecule%nvib))
          ! re-allocate atom arrays
          do i = 1, newMolecule%nvib
            allocate(molecule%normodes%vibration(i)%atom(1:numat),stat=alloc_err)
          end do

          ! Reassign vibrations
          do i = 1, newMolecule%nvib
              molecule%normodes%vibration(i) = newMolecule%normodes%vibration(i)
           ! print *, 'freq redef', molecule%normodes%vibration(i)%freq
          end do

          molecule%nvib = newMolecule%nvib

          write(fout,*) 'Molecule cropped to subset: ', incvib 

          return
        end subroutine getSubset

end module system_type
