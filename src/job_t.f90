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
    integer :: nmeth
    character(len=10), allocatable :: method(:)
    type(transf_t), allocatable :: trns
    type(fc_t), allocatable :: fc(:)
    type(dos_t) :: dos
  end type job_t

end module job_type
