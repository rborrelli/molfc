module fc_type

  use parameters
  use sysop, only : get_state_from_id
  use system_type

  implicit none

  private

!  public :: set_active_space, get_fill_active_space

  type, public :: spectrum_t
    character(len=80)  :: file = "spectrum.dat", &
                          shape = "gaussian",  &
                          form = "formatted"
    character(len=4) :: type = "fc" ! available : fc, emi, abs 
    real(kind=dp) ::  fwhm = -1.D0  ! Default value is < 0. This means: no convolution!
    real(kind=dp) ::  tol = 1.0D-5
    real(kind=dp) ::  Temp = 298.15 ! Default temperature for spectrum simulation.
    real(kind=dp) ::  dE = 0.5  ! Lo step di energia per la convoluzione in cm-1.
    real(kind=dp) ::  emin = 0.0
    real(kind=dp) ::  emax = -1.0 ! Se negativa allora viene automaticamente messa al valore massimo.
  end type spectrum_t

  type, public :: kubo_t
	  character(len=STRLEN) :: file = "spectrum.dat"
	  integer :: pow = 13 ! 2^pow is the number fo points used in the FFT
	  real(kind=dp) :: Temp = 298.15 ! Default temperature for spectrum simulation
	  real(kind=dp) :: OmR = 100000.0 ! energy range (E = [-OmR/2, +OmR/2])
	  real(kind=dp) :: fwhm = -1.0 ! full width at half maximum of the spectrum peaks
	  logical :: on = .false. ! if ON is TRUE the GF method is used
  end type kubo_t

  type, public :: activespace_t
    integer ::  nact = 0
    integer, allocatable :: vibid(:)
    integer, allocatable :: nqmax(:)
    real(kind=dp), allocatable :: freq(:)
  end type activespace_t

  type, public :: activem_t
    integer :: id ! Numero identificativo del modo.
    integer :: nq = 1 ! Numero di quanti sul modo MID. Il valore minimo per essere un modo attivo e' 1.
    real(kind=dp) :: freq = zero
  end type activem_t

  type, public :: active_t
    integer :: nact = 0
	logical :: autex = .false. ! automatic determination of excitation level
    character(len=STRLEN) state
    character(len=STRLEN) molecule
    type(activem_t), allocatable :: mode(:) ! ha dimensione nact
  end type active_t

  type, public :: incvib_t
    character(len=STRLEN) :: mol
    integer, allocatable :: id(:) ! ids of the group vibrations.
  end type incvib_t

  type, public :: group_t
    integer :: id = 0 ! identificative number of the group
    integer :: nvib = 0 ! number of vibrations in the group.
    type(incvib_t) :: incvib
    type(active_t), allocatable :: active(:) ! active space of the group
  end type group_t

  type, public :: fc_t
    character(len=STRLEN) :: bra = ""
    character(len=STRLEN) :: ket = ""
    integer :: ngroup = 0 ! ngroup deve essere almeno pari al numero di molecole.
    integer :: order = 2 ! default to second order perturbation formulae
    integer :: nclasses = 3
    real(kind=dp) :: ftol = 0.01, &
                     tol = 0.1
    real(kind=dp) :: Temp = 0.0 ! This is used in the fcclass Boltzmann
    integer :: printlevel = 1
    logical :: printfc = .true., &
               fcspec = .false., &
               debug = .false.,  &
               force = .false.,  &
               pert = .false.,  &
               recur = .false.,  &
               berk = .false.,   &  ! activate Berkowitz algorithm
               fcht = .false.    ! if true, triggers the calculation of Herzberg-Teller effects.
	logical :: class = .false. ! use Santoro algorithm for automatic spectrum calculation
!    type(active_t), allocatable :: active(:)
    type(group_t), allocatable :: group(:)
    type(spectrum_t) :: spectrum
    type(kubo_t) :: kubo
  end type fc_t

  type, public :: fci_t
    real(kind=dp), allocatable :: FC(:,:)
    real(kind=dp), allocatable :: FCV(:)
  end type fci_t
  
end module fc_type
