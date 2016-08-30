!     
! File:   hypfun.f90
! Author: lello
!
! Created on February 17, 2009, 2:24 PM
!

module hypfun
	
	use parameters

    implicit none

	private

	public :: sinh, cosh, tanh, coth
	
	interface sinh
		module procedure sinh
		module procedure sinhr
		!module procedure vsinhr
		!module procedure vsinh
	end interface sinh

	interface cosh
		module procedure cosh
		module procedure coshr
		!module procedure vcosh
		!module procedure vcoshr
	end interface cosh

	interface tanh
		module procedure tanh
		!module procedure vtanh
	end interface tanh

	interface coth
		module procedure coth
		!module procedure vcoth
	end interface coth

	contains

	!==========================================================================!

	elemental function sinhr(r) result(c)

	real(kind=dp), intent(in) :: r
	real(kind=dp) :: c

	c = 0.5_dp*(exp(r) - exp(-r))

	end function sinhr

	!==========================================================================!

	function vsinhr(r) result(c)

	real(kind=dp) :: r(:)
	real(kind=dp) :: c(1:size(r))

	c = 0.5_dp*(exp(r) - exp(-r))

	end function vsinhr

	!==========================================================================!

	elemental function sinh(z) result(c)

	complex(kind=dpc), intent(in) :: z
	complex(kind=dpc) :: c

	c = 0.5_dp*(exp(z) - exp(-z))

	end function sinh

	!==========================================================================!

	elemental function cosh(z) result(c)

	complex(kind=dpc), intent(in) :: z
	complex(kind=dpc) :: c

	c = 0.5_dp*(exp(z) + exp(-z))

	end function cosh

	!==========================================================================!

	elemental function coshr(z) result(c)

	real(kind=dp), intent(in) :: z
	real(kind=dpc) :: c

	c = 0.5_dp*(exp(z) + exp(-z))

	end function coshr

	!==========================================================================!

	elemental function tanh(z) result(c)

	complex(kind=dpc), intent(in) :: z
	complex(kind=dpc) :: z2
	complex(kind=dpc) :: c
	complex(kind=dpc) :: r1

	r1 = cmplx(1.0_dp,0.0_dp)

	z2 = 2.0_dp*z
	if (real(z) > zero) then
		c = (r1 - exp(-z2))/(r1 + exp(-z2))
	else
		c = (exp(z2) - r1)/(exp(z2) + r1)
	end if

	end function tanh

	!==========================================================================!

	elemental function coth(z) result(c)

	complex(kind=dpc), intent(in) :: z
	complex(kind=dpc) :: z2
	complex(kind=dpc) :: c

	z2 = 2.0_dp*z

	if (real(z) > zero) then
		c = (one + exp(-z2))/(one - exp(-z2))
	else
		c = (exp(z2) + one)/(exp(z2) - one)
	end if

	end function coth

	!==========================================================================!

	function vsinh(z) result(c)

	complex(kind=dpc) :: z(:)
	complex(kind=dpc) :: c(1:size(z))

	c = 0.5_dp*(exp(z) - exp(-z))

	end function vsinh

	!==========================================================================!

	function vcosh(z) result(c)

	complex(kind=dpc) :: z(:)
	complex(kind=dpc) :: c(1:size(z))

	c = 0.5_dp*(exp(z) + exp(-z))

	end function vcosh

	!==========================================================================!

	function vcoshr(z) result(c)

	real(kind=dpc) :: z(:)
	real(kind=dpc) :: c(1:size(z))

	c = 0.5_dp*(exp(z) + exp(-z))

	end function vcoshr

	!==========================================================================!

	function vtanh(z) result(c)

	complex(kind=dpc) :: z(:)
	complex(kind=dpc) :: c(1:size(z))
	integer :: i

	forall (i=1:size(z)) c(i) = tanh(z(i))

	end function vtanh

	!==========================================================================!

	function vcoth(z) result(c)

	complex(kind=dpc) :: z(:)
	complex(kind=dpc) :: c(1:size(z))
    integer :: i

	forall (i=1:size(z)) c(i) = coth(z(i))

	end function vcoth


end module hypfun
