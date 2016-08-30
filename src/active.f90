!     
! File:   active.f90
! Author: lello
!
! Created on May 25, 2007, 5:21 PM
!
! NOTE:
!----------------------------------------------------------------------------
! For a correct definition of the active space it is important to know a few
! things. In the fc calculation job one has to specify the bar and ket
! electronic states. Then the active space of each is given.
! The numbering of normal modes is fundamental. In the fcjob this numbering
! starts from bra and continues with the ket. Normal modes of bra are
! numbered from 1 to N, and normal modes of ket are numbered from N+1 to 2N,
! however this is not the case for the transformation of normal coordinates
! and hence for the numbering of normal noordinates in the matrices M and Q.
! In those cases indeed the numeration obtained from the input sequence is
! applyed. This causes some problems. A possible solution would be to give
! the same enumeration to both the <transformation> job and <fc> job. This
! however would require additional keywords. On the other hand I'd prefer to
! keep the two jobs as separate as possible.
! The solution adopted here is a temporary one and simply enumerates normal
! modes according to the M and Q matrix in order to have a correct FC
! calculation.
!----------------------------------------------------------------------------

module active

  use parameters
  use system_type
  use sysop, only : get_state_from_id
  use fc_type

  implicit none

  private

  public :: get_active_space, get_group_active_space_st
            
  public :: set_active_space, set_no_active_space

  interface get_active_space
      module procedure get_group_active_space_st
  end interface get_active_space

  contains

    !==========================================================================!

    function get_group_active_space_st(system,group,istate) result (acs)
    !---------------------------------------------------------------------!
    ! This functions retrieve the active space of a group of vibrations.  
    ! In order to properly number the vibrations the state they belong to 
    ! is also specified.
    !---------------------------------------------------------------------!
    
    type(system_t), intent(in) :: system
    type(group_t), intent(in) :: group
    integer, intent(in) :: istate
    !character(len=*), intent(in) :: molid
    type(activespace_t) :: acs

    integer i, j, i0, i1, iact1, is, im, ioff
    integer idtemp, nqtemp
    real(kind=dp) frtemp

    if (allocated(acs%vibid)) deallocate(acs%vibid)
    if (allocated(acs%nqmax)) deallocate(acs%nqmax)
    if (allocated(acs%freq)) deallocate(acs%freq)

    ! Find size of active space of the group
    do i = 1, size(group%active(:))
      is = get_state_from_id(system,group%active(i)%state)
      !print *, 'istate ', is, istate
      if (is == istate) acs%nact = acs%nact + group%active(i)%nact
    end do

    !print *, 'acsnact ', acs%nact
    if (acs%nact == 0) return ! no active space found

    allocate(acs%vibid(1:acs%nact),acs%nqmax(1:acs%nact),acs%freq(1:acs%nact))

    j = 0
    i0 = 0;
    ! Set the active space for a single group
    do i = 1, size(group%active(:))
      is = get_state_from_id(system,group%active(i)%state)
      if (is == istate) then
        acs%vibid(1:group%active(i)%nact) = group%active(i)%mode(:)%id
        acs%nqmax(1:group%active(i)%nact) = group%active(i)%mode(:)%nq
        acs%freq(1:group%active(i)%nact) = group%active(i)%mode(:)%freq
      end if
    end do

    ! Sort active space modes in ascending order according to their ID
    ! Bubble sort. Not very efficient but very easy.
    do i = 1, size(acs%vibid)
      do j = i + 1, size(acs%vibid)
        if (acs%vibid(i) > acs%vibid(j)) then
          idtemp = acs%vibid(i)
          acs%vibid(i) = acs%vibid(j)
          acs%vibid(j) = idtemp
          nqtemp = acs%nqmax(i)
          acs%nqmax(i) = acs%nqmax(j)
          acs%nqmax(j) = nqtemp
          frtemp = acs%freq(i)
          acs%freq(i) = acs%freq(j)
          acs%freq(j) = frtemp
        end if
      end do
    end do

    ! Renumbering of VIBIDs. According to their position with respect to the 
    ! included vibrations. This is fundamental!!!!!
    ! Esempio:
    ! 1 5 7 9 10 11 12 13 : vibrazioni incluse <include>
    !     7 9       12    : vibrazioni attive <active> (vibid iniziali)
    !---------------------
    !     3 4       7     : posizioni                  (vibid finali)
    do i = 1, acs%nact 
      acs%vibid(i) = count(group%incvib%id(:) < acs%vibid(i)) + 1
    end do

    end function get_group_active_space_st

    !==========================================================================!

    function set_active_space(system,nqmax,istate) result (acs)
    !---------------------------------------------------------------------!
    ! This functions retrieve the active space of a group of vibrations.  
    ! In order to properly number the vibrations the state they belong to 
    ! is also specified.
    !---------------------------------------------------------------------!
    
    type(system_t), intent(in) :: system
    integer, intent(in) :: istate
    integer, intent(in) :: nqmax(1:system%state(istate)%molecule%nvib)
    type(activespace_t) :: acs

    integer i, j, k, i0

    if (allocated(acs%vibid)) deallocate(acs%vibid)
    if (allocated(acs%nqmax)) deallocate(acs%nqmax)
    if (allocated(acs%freq)) deallocate(acs%freq)

    ! Find size of active space of the group
    acs%nact = system%state(istate)%molecule%nvib

    !print *, 'acsnact ', acs%nact
    if (acs%nact == 0) return ! no active space found

    allocate(acs%vibid(1:acs%nact),acs%nqmax(1:acs%nact),acs%freq(1:acs%nact))

    ! Set the active space for a single group
	do j = 1, system%state(istate)%molecule%nvib
        acs%vibid(j) = j
        acs%freq(j) = system%state(istate)%molecule%normodes%vibration(j)%freq        
        acs%nqmax(j) = nqmax(j) ! this will be set in fcwd_Metz
    end do

	return
    end function set_active_space

	!==========================================================================!
	
    function set_no_active_space(system,istate) result (acs)
    !---------------------------------------------------------------------!
    ! This functions retrieve the active space of a group of vibrations.  
    ! In order to properly number the vibrations the state they belong to 
    ! is also specified.
    !---------------------------------------------------------------------!
    
    type(system_t), intent(in) :: system
    integer, intent(in) :: istate
    !character(len=*), intent(in) :: molid
    type(activespace_t) :: acs

    integer i, j, k, i0

    if (allocated(acs%vibid)) deallocate(acs%vibid)
    if (allocated(acs%nqmax)) deallocate(acs%nqmax)
    if (allocated(acs%freq)) deallocate(acs%freq)

    ! Find size of active space of the group
    acs%nact = 0

	return
    end function set_no_active_space

end module active
