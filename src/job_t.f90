module job_type

  use parameters
  use transf_type
  use fc_type
  use dos_type

  implicit none

  private

  !========================================
  !  Derived  Data Type Definitions for Job
  !========================================

  type, public :: job_t
    type(fc_t), allocatable :: fc(:)
    character(len=10), allocatable :: method(:)
    type(transf_t), allocatable :: trns
    type(dos_t) :: dos
    integer :: nmeth
  end type job_t

end module job_type
