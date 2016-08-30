program mconvolution

  implicit none
  
  type llist_t
    real(8) fc
    real(8) en
    character(len=80) :: desc
    type(llist_t), pointer :: node => null()
  end type llist_t

  integer, parameter :: NMAXFILE = 10, &
                        NFCMAX = 50000000
  integer, parameter :: iuf = 99, &
                        fout = 12
  real(8), parameter :: fctol = 1.0d-8
  real(8)               fcthr 
  integer i, j, ifile, nfile, nf, nfcsel, NFC1, NFC2
  real(8), allocatable :: e1(:), fc1(:), e2(:), fc2(:)
  real(8) :: fctmp
  integer :: count_arg, stat_arg
  character(len=40) :: finpname(1:NMAXFILE)
  character(len=40) :: foutname
  character(len=50) :: buffer, fname
  character(len=150) :: line
  character(len=80), allocatable :: desc1(:), desc2(:)
  character(len=20) :: st1, st2

  type(llist_t), pointer :: fclist, current

  integer, external :: nargs

  count_arg = nargs()
  if (count_arg < 2) stop 100 ! error

  ! get the thresold for FC factors
  call getarg(1,buffer,stat_arg)
  print *, buffer
  read(buffer,*)fcthr

  ! get the files to convolute
  nfile = 0; i = 2
  do while (i <= count_arg - 1)
    call getarg(i,buffer,stat_arg)
  print *, buffer
    if (stat_arg .lt. 0) stop 99 ! error
    if (buffer == "-o") then
      i = i + 1
      call getarg(i,foutname,stat_arg)
      if (stat_arg .lt. 0) stop 98 ! error
      i = i + 1
      print *, foutname
    else ! get input file names
      i = i + 1
      nfile = nfile + 1
      fname = buffer
      print *, i, nfile, fname
      !call getarg(i,fname,stat_arg)
      open(unit=nfile+iuf,file=fname,status='old')
      ! update number of input files
    end if
  end do


  ! Apre i files di input
  !do i = 1, nfile
  !  fname = finpname(i)
  !  open(i+iuf,fname,status='old')
  !end do
  ! Apre il file di output 
  open (unit=fout,file=foutname,status='unknown')
 
  ! legge i fc nel file 1
  ifile = 1+iuf
  do 
    read(ifile,*,end=11)NFC1
  end do  

11  rewind ifile
  print *, NFC1

  allocate(e1(1:NFC1),fc1(1:NFC1),desc1(1:NFC1))
  do i = 1, NFC1
    !read(ifile,*) e1(i),fc1(i)
    read(ifile,'(a)') line
    read(line,*) e1(i),fc1(i)
    desc1(i)=line(index(line,'<'):index(line,'>'))
    !print *, desc1(i)
  end do

  ! legge i fc nel file i-simo e li moltiplica con quelli del file precedente.
  ! quindi aggiorna i fc selezionati e procede con i file successivi.
  do nf = 2, nfile

    ifile = nf+iuf
    do 
      read(nf+iuf,*,end=12)NFC2
    end do  

12  rewind ifile

    print *, NFC2

    allocate(e2(1:NFC2),fc2(1:NFC2),desc2(1:NFC2))
    do i = 1, NFC2
!      read(ifile,*) e2(i),fc2(i)
      read(ifile,'(a)') line
      read(line,*) e2(i),fc2(i)
      desc2(i)=line(index(line,'<'):index(line,'>'))
      !print *, desc2(i)
    end do

    ! inizializza la linked list
    nullify(fclist)

    nfcsel = 0
    do i = 1, NFC1
      do j = 1, NFC2
        fctmp = fc1(i)*fc2(j)
        ! metto i FC selezionati in una linked list
        !if (fctmp >= fctol) then
        if (fctmp >= fcthr) then
          nfcsel = nfcsel + 1
          if (nfcsel > nfcmax) exit
          allocate(current)
          current%fc = fctmp
          current%en = e1(i)+e2(j)
          st1 = desc1(i); st2=desc2(j)
          current%desc = st1(1:len_trim(st1))//st2(1:len_trim(st2))
          current%node => fclist 
          fclist => current
        end if
      end do
    end do

    !sposto i fc e le rispettive energie dalla linked list agli array fc1 ed e1

    NFC1 = nfcsel
    deallocate(fc1,e1,desc1)
    allocate(fc1(1:NFC1),e1(1:NFC1),desc1(1:NFC1))

    current => fclist ! point to the head of list
    i = 0
    do
      if (.not.associated(current)) exit
      i = i + 1
      fc1(i) = current%fc
      e1(i) = current%en
      desc1(i) = current%desc
      current=> current%node
    end do

    ! svuoto la linkes list
    current => fclist ! point to the head of list
    do 
      if (.not.associated(current)) exit
      fclist => current%node
      deallocate(current)
      current => fclist
    end do

    ! fc2 ed e2 devono contenere i fc del nuovo file
    deallocate(fc2,e2,desc2)

  end do

  print *, NFC1
  do i = 1, NFC1
      !write(fout,*) e1(i), fc1(i), desc1(i)
      write(fout,'(2x,f10.2,2x,g14.5,3x,a)') e1(i), fc1(i), desc1(i)
  end do

end program mconvolution
