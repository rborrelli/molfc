program boltz

implicit none
real(8), parameter :: cfac = 
real(8), allocatable :: qvib(:), freq(:)

read(*,+) n
allocate(qvib(1:n))
allocate(freq(1:n))

do i = 1, n
  read(*,*) freq(i)
end do

do i = 1, n
  theta = freq(i)*cfac
  qvib(i) = exp(-theta/2.0)/(1-exp(-theta))
end do

qvib_tot = product(qvib(:))


end program boltz

