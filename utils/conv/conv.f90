program convolution


   implicit none

   type :: partition_t
     real(8) :: domega
     real(8), allocatable :: dens(:)
     real(8), allocatable :: energy(:)
   end type partition_t

   integer, parameter :: iuf = 99, fout = 12

   integer :: count_arg, stat_arg
   integer :: i, j, k, nj, iemin, iemax, NGR, nfile
   integer :: N, N1, N2, ifile(2)
   real(8), allocatable  :: t0(:), t1(:), dens(:), eminj(:)
   real(8) :: emin, emax, egrain, FCF, eshift, e0
   character(len=50) :: buffer, fname, foutname
   type(partition_t), allocatable :: part(:)
 
  integer, external :: nargs

 !  read input data
  count_arg = nargs()
  if (count_arg < 2) stop 100 ! error

  ! get the thresold for FC factors
  call getarg(1,buffer,stat_arg)
  print *, buffer

  ! get the files to convolute
  nfile = 0; i = 1
  do while (i <= count_arg - 1)
    print *, i
    call getarg(i,buffer,stat_arg)
  print *, 'buffer ', buffer
    if (stat_arg .lt. 0) stop 99 ! error
    if (buffer == "-o") then
      i = i + 1
      call getarg(i,foutname,stat_arg)
      if (stat_arg .lt. 0) stop 98 ! error
      i = i + 1
      print *, foutname
    else ! get input file names
      nfile = nfile + 1
      ifile(i)=nfile+iuf
      fname = buffer
      print *, 'file to read ', i, nfile, fname, ifile(i)
      !call getarg(i,fname,stat_arg)
      open(unit=ifile(i),file=fname,status='old')
      i = i + 1
      ! update number of input files
    end if
  end do

  !end do
  ! Apre il file di output 
  open (unit=fout,file=foutname,status='unknown')
 
  ! Only two densities can be convoluted at the moment but the algorithm can be extended to any number
  NGR = 2

  allocate(part(1:NGR))

  do k = 1, NGR
    print *, 'Reading file ',  ifile(k)
    read(ifile(k),*) N
    print *, 'N = ', N
    allocate(part(k)%dens(1:N),part(k)%energy(1:N)) 
    do i = 1, N
      read(ifile(k),*) part(k)%energy(i),part(k)%dens(i)
    end do
    part(k)%domega=part(k)%energy(2)-part(k)%energy(1) ! assuming it is a constant through the interval
    print *, k, part(k)%domega
  end do


 ! Set energy range for convolution
    emin = -1000.0
    emax = 20000.0
    egrain = 1.0
    iemin = int(emin/egrain)
    iemax = int(emax/egrain)

    print *, 'iemax: ', iemax, 'egrain: ', egrain
 	! initialize t: t(i) stores the actual number of states at energy e(i)
    allocate(t0(0:iemax),t1(0:iemax))
	t0 = 0.0; t1 = 0.0
	t0(0) = 1.0

    allocate(eminj(1:NGR))
    forall (i=1:NGR) eminj(i) = part(i)%energy(1) 
    eshift = minval(eminj)
    if (eshift >= 0.0) eshift = 0.0 ! se le densità sono zero per valori di energia negativi non c'è bisogno di alcuno shift.
    e0 = sum(eminj)

    print *, 'Applying shift to energy levels ', eshift

    ! Ciclo sui gruppi di modi vibrazionali
        print *, 'Starting convolution...'
        do j = 1, NGR
                ! Ciclo sull'intero rane di energia nel quale ci interessa la FCWD
                do i = iemax, 0, -1
                        ! Ciclo sui livelli energetici della densità j
                        do k = 1, size(part(j)%energy)
                                nj = nint((part(j)%energy(k)-eshift)/egrain) ! shift energies to zero
                                !nj = nint((k-1)/egrain) ! shift energies to zero
                                ! The exit method can be used only if energy levels are sorted in ascending order
                                if (nj > i) exit
                                !if (nj > i) cycle
                                FCF = part(j)%dens(k)*part(j)%domega ! mean value of FCF in that specific interval
                                !print *, i, nj, i-nj, t0(i-nj), FCF
                                t1(i) = t1(i) + t0(i-nj)*FCF
                        end do
                end do
                t0 = t1
                t1 = 0.0
        end do

	allocate(dens(0:iemax))

	dens = 0.0

    dens(0) = t0(0)/egrain
    do i = 1 , iemax
       dens(i) = t0(i)/egrain
    end do

    !Here I should put the energy levels back to their original value
    do i = 0, iemax
      !write(fout,*) egrain*i+e0, dens(i)
      write(fout,*) egrain*i+e0+eshift, dens(i)
    end do

end program convolution

