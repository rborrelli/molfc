!---------------------------------------------------------------------
! This module contains the main variables used throughout the program.
!---------------------------------------------------------------------
module xmvar

    use job_type
    use proc_type
    use system_type
  
    implicit none

    private

    type(system_t), public            :: system ! Whole system (states and molecules)
    type(proc_t), public              :: proc   ! Processing instructions for reference frames
    type(job_t), public, allocatable  :: job(:) ! Actual jobs (Franck-Condon, Dynamics)
  
end module xmvar
