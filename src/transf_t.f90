module transf_type

  use parameters
  
  implicit none


  type, public :: mod_shift
    character(len=STRLEN) :: file = ""
    character(len=STRLEN) :: data = ""
    character(len=STRLEN) :: type = "DIM"
    character(len=STRLEN) :: coord = ""! tipo di coordinate (specifica se adimensionali o massa pesate)
  end type mod_shift

  type, public :: mod_rot
    character(len=STRLEN) :: file = ""
    character(len=STRLEN) :: data = ""
    character(len=STRLEN) :: coord = ""! tipo di coordinate (specifica se adimensionali o massa pesate)
  end type mod_rot
  
  type, public :: axsw_t
	  logical :: on = .true.
	  integer :: solution = -1
  end type axsw_t

  type, public :: transf_t
    character(len=80) :: icfile = "INTVAL"
    character(len=80) molecule
    character(len=10) :: coord = "xyz"
    real(kind=dp), allocatable :: KM(:)  ! normal modes shift vector
    real(kind=dp), allocatable :: JM(:,:) ! normal modes rotation matrix
    real(kind=dp), allocatable :: T0(:,:) ! axis switching matrix
    logical :: debug = .false., &
               model = .false., &
			   setic = .false., &
			   dusch = .true. ! if false set Duschinsky  matrix to identity
    logical :: cartesian = .true., &
               internal = .false., &
               intauto  = .true., & ! If internal coordinates are chose: default is automatic.
               nonlin   = .false., &
               natint = .false.
    logical :: lstre = .true., &
               lbend = .true., &
               lopb  = .false., & ! Out of plane bendings are not included by default.
               ldih  = .true.
    logical :: tm_rotate = .true. ! transform the transition moment to the ground state coordinates
    integer :: printlevel = 1
    type(mod_shift) :: modk
    type(mod_rot) :: modr
	type(axsw_t) :: axsw
  end type transf_t

end module transf_type

