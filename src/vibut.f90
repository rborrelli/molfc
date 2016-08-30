module vibut

  use parameters
  use xmvar, only : system

  public :: next_quantum_numbers, get_quantum_numbers, get_index, Energy

  interface next_quantum_numbers
    module procedure next_quantum_numbers
    module procedure next_q
  end interface
  
  contains

  ! ---------------------------------------------------------        

  subroutine next_q ( v, vmax, nvex, exvib )

  integer, intent(in out)  :: v(:)
  integer, intent(in)  :: vmax(:)
  integer, intent(in)  :: nvex
  integer, intent(in)  :: exvib(:)

  integer i

  if (any(v < vmax)) then

    do i = 1, nvex
      if (v(exvib(i)) < vmax(exvib(i))) then 
        v(exvib(i)) = v(exvib(i)) + 1; exit
      else
        if (v(exvib(i+1)) < vmax(exvib(i+1))) then
          v(exvib(i+1)) = v(exvib(i+1)) + 1
          v(1:exvib(i)) = 0
          exit
        end if
      end if
    end do

    return
  end if

  return
  end subroutine next_q

  ! ---------------------------------------------------------        

  subroutine next_quantum_numbers ( v, vmax )

  integer, intent(in out)  :: v(:)
  integer, intent(in)  :: vmax(:)

  integer i

  if (any(v < vmax)) then

    do i = 1, nmd - 1
      if (v(i) < vmax(i)) then 
        v(i) = v(i) + 1; exit
      else
        if (v(i+1) < vmax(i+1)) then
          v(i+1) = v(i+1) + 1
          v(1:i) = 0
          exit
        end if
      end if
    end do

    return
  end if

  return
  end subroutine next_quantum_numbers

  ! --------------------------------------------------------        
  subroutine get_quantum_numbers ( kq, id, v )
  
  ! kq = indice dello stato
  ! id(i) = vmax(i) + 1
  integer, intent(in)   :: kq, id(:)
  integer, intent(out)  :: v(:)
  integer  kkq, pp, i

  kkq = kq - 1

  v(1) = mod(kkq, id(1))
  pp = kkq
  do i = 2, nmd
    pp = int(pp/id(i-1))
    v(i) = mod(pp, id(i))
  end do

  return
  end subroutine get_quantum_numbers

  !-------------------------------------------------------------        
!
!  function get_index (v, w, ipdv, ipdw) result(ind)
!
!  integer ind
!  integer, intent(in)   ::  w(1:nmd), v(1:nmd)
!  integer, intent(in)   ::  ipdw(1:nmd), ipdv(1:nmd)
!
!  integer i, ind1, ind2, idim, ipp, idiff
!
!  ind1 = v(state(1)%exvib(1))
!  do i = 2, state(1)%nex
!    ipp = v(state(1)%exvib(i)) * ipdv(state(1)%exvib(i))
!    ind1 = ind1 + ipp
!  end do
!
!  ind2 = w(state(2)%exvib(1))
!  do i = 2, state(2)%nex
!    ipp = w(state(2)%exvib(i))*ipdw(state(2)%exvib(i))
!    ind2 = ind2 + ipp
!  end do
!  
!  ind = ind1 + ind2*product(state(1:2)%nvibstate)
!  ind = ind + 1
!  
!  return
!  end function get_index 

  ! -------------------------------------------------------------        

  function Energy (omega, v) result (MDHOE)
  
  real(kind=dp) MDHOE
  real(kind=dp), intent(in)  :: omega(:)
  integer, intent(in)  :: v(:)
  integer i

  ! I don't check the size of the two arrays so be careful in calling this
  ! function.

  MDHOE = zero
  MDHOE = sum(dble(v)*omega)

  return
  end function Energy 

end module  vibut
