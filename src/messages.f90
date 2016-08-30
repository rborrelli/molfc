module messages

  use parameters

  implicit none
  
  public :: message
  
  contains

    function message(which) result (written)

    integer, intent(in) :: which
    integer written

    select case(which)
      
      case (1)
        write(fout,*) '   Diagonalizing Hamiltonian...This can take a while'
        written = 0

      case default
        write(fout,*) 'No message associated.'
        written = 1

    end select

    return
    end function message

end module messages
