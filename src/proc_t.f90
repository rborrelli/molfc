module proc_type

  use parameters
  
  implicit none

  type, public :: reord_t
    character(len=STRLEN) :: data=""
    character(len=STRLEN) :: state
    character(len=STRLEN) :: reference=""
    character(len=STRLEN) :: molecule=""
    character(len=STRLEN) :: ord=""
    !logical :: all = .false. ! Reorder all the molecules. not yet implemented.
  end type reord_t

  type, public :: frame_t
    character(len=STRLEN) :: state
    character(len=STRLEN) :: molecule
    integer :: axis(1:3)
  end type frame_t

  type, public :: subset_t
    character(len=IDLEN) :: state
    integer :: nvib = 0
    integer, allocatable :: incvib(:)
  end type subset_t

  type, public :: proc_t
    type(reord_t), allocatable :: reorder(:)
    type(frame_t), allocatable :: frame(:)
    type(subset_t), allocatable :: subset(:)
  end type proc_t

end module proc_type
