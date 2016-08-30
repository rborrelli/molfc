module fc_h

  use parameters
  use errors
  use iomatrix
  use xmvar, only : job
  use matfun
  use system_type
  use sysop, only : get_bmat, get_normal_modes, get_nat_bmat
  use intc, only : set_equilibrium_intc
  use transf_type
  use fc_type

  implicit none

  real(kind=dp),  allocatable :: y(:), p(:,:), xm1(:,:), xm12(:,:)

  real(kind=dp), public, allocatable :: PR(:,:), TR(:), RR(:,:), PIR(:,:), TIR(:)
  !-----------------------------------------------+
  ! The super-matrices for the recursion formula. |
  !-----------------------------------------------+
  real(kind=dp), public, allocatable :: M(:,:), Q(:) 
  
  ! L'opzione "target" serve per il modulo fcint_doktorov.
  !real(kind=dp), public, allocatable, target :: M(:,:), Q(:) 
  
end module fc_h
  
