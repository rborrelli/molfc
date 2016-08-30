module vecmath

  use parameters
  use errors

  implicit none

  private

  public :: cross, norm

  contains

    !-----------------------------------------------------------!

    function norm(v) result (nr)

    real(kind=dp), intent(in) :: v(:)
    real(kind=dp) nr

    nr = dot_product(v,v)
    nr = sqrt(nr)

    return
    end function norm

    !-----------------------------------------------------------!

    function cross(u,v) result (cp)

    ! Cross product U x V.
    real(kind=dp), intent(in) :: u(3), v(3)
    real(kind=dp) cp(3), s
    integer ierr

    cp(1) = u(2)*v(3) - u(3)*v(2)
    cp(2) = u(3)*v(1) - u(1)*v(3)
    cp(3) = u(1)*v(2) - u(2)*v(1)

    !s = norm(cp)
    !if(s.le.0.0001d0) ierr= error(0,"Cross product error!")
    !cp = cp/s

    return
    end function cross

end module
