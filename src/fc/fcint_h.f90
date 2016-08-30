module fcint_h

    use parameters
    use errors
    use xmvar, only : system
    use vibut, only : energy
    use system_type
    use sysop
    !use job_type
    use fc_type

    implicit none

    type, public :: fccl_t
        integer :: n  ! size of the class
        integer, allocatable :: ivib(:) ! the size is n
        integer, allocatable :: gvib(:) ! ground state excited vibrations
        integer, allocatable :: zmax(:) ! the size is n
        real(kind=dp), allocatable :: fc(:) ! the size is product(wmax+1)
        real(kind=dp), allocatable :: fcht(:) ! the size is product(wmax+1)
    end type fccl_t

    real(kind=dp), public, allocatable    :: MDFC(:)
    real(kind=dp), public, allocatable    :: FCHT(:)

end module fcint_h
