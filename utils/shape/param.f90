module param

  implicit none

  public

  integer, parameter :: dp = kind(1.0D0)
  integer, parameter :: sp = kind(1.0)


  real(kind=dp), parameter ::         &    
                        PI=3.141592653589793238462643383279502884197_dp,     &
                        PIO2=1.57079632679489661923132169163975144209858_dp, &
                        TWOPI=6.283185307179586476925286766559005768394_dp,  &
                        one      = 1.0_dp,            & 
                        zero     = 0.0_dp            

  integer, parameter ::  fout = 6, & 
                         finp = 10, &
                         flog = 9

end module param
