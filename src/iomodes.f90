module iomodes

    use parameters
    use isotopes
    use errors
    use xmvar, only : system
    use system_type, only : vibration_t

    implicit none
  
    public :: read_normal_modes_gau, read_normal_modes_gaunohp, &
                read_normal_modes_tm, read_normal_modes_molden, &
                read_normal_modes_mopac, read_normal_modes_gamess, &
                read_normal_modes_demon_nohp
  
contains

    ! ============================================================================ !

    subroutine read_normal_modes_gau (iu, istate, imol, stat)

        ! In uscita stat = .true. se i modi normali sono stati letti
        !           stat = .false. se i modi normali non sono stati letti
        integer, intent(in) :: iu, istate, imol
        logical, intent(out) :: stat

        integer, parameter :: MAXCOL = 5
        integer icol, nr, nc, inmd, i, j, k, ios, i1, i2, i3
        real(kind=dp), allocatable  ::  aus(:), g(:)
        integer, allocatable :: mkc(:)
       
        character (len=80) line
        character (len=15) str
        character (len=3) crj(1:3) 

        rewind iu

        ! Solo per specie non lineari!! 
        ! Devo trovare una soluzione per specie lineari.
        nc = 3*system%state(istate)%molecule%structure%numat - 6
        nr = 3*system%state(istate)%molecule%structure%numat

        allocate(g(1:MAXCOL), aus(1:MAXCOL), mkc(1:MAXCOL))

        aus = ZERO

        stat = .false.

        ! --- Trova la stringa 'normal coordinates nell'input

        line = BLANK
        do while (index(line,'normal coordinates') == 0)
            read (iu,'(a80)',iostat=k,end=222)line
        end do

        RDWM:   do

            mkc = 0
            read (iu,'(a80)',end=222)line
            read (line,*,iostat=ios) (mkc(i),i=1,MAXCOL)
             
            icol=count(mkc > 0)

            if (icol > 0) then

                ! Read symmetry labels
                read (iu,'(a80)',end=222)line
                read (line,*,iostat=ios) (system%state(istate)%molecule%normodes%vibration(mkc(i))%symlab, i=1,icol)
!                do i = 1, icol
!                print *, 'SYMLAB : ', i, system%state(istate)%molecule%normodes%vibration(mkc(i))%symlab
!                end do

                ! Assegna le ID alle vibrazioni di ciscuna molecola.
                system%state(istate)%molecule%normodes%vibration(mkc(1):mkc(icol))%id = mkc(1:icol)

                ! --- Le frequenze che leggo sono positive e quindi le pongo uguali a -1
                ! --- per sapere se sono state lette o meno. Se la lettura e' corretta
                ! --- devono essere tutte positive

                g = -ONE
                do
                    read (iu,'(a80)',iostat=k)line
                    if (k /= 0) exit
                    read (line,*,iostat=ios) str, str, (g(i),i=1,icol)
                    if (ios == 0) then
                        system%state(istate)%molecule%normodes%vibration(mkc(1):mkc(icol))%freq = g(1:icol)
                    end if
                    ! Controlla che siano state effettivamente lette le frequenze.
                    if (count(g(1:icol) > ZERO) == icol) exit
                end do
            
                inmd = 0
                do
                    read (iu,'(a80)',iostat=k,end=222)line
                    if (k /= 0) stop "Cannot read line"
                    read(line,*,iostat=ios)i1,i2,i3,(aus(j),j=1,icol)
                    if (ios == 0) then
                        inmd = inmd + 1
                        do j = 1, icol
                            system%state(istate)%molecule%normodes%vibration(mkc(j))%atom(i2)%elem = elmnts(mapel(i3))
                            system%state(istate)%molecule%normodes%vibration(mkc(j))%atom(i2)%d(i1) = aus(j)
                        end do
                    end if
                    if ( inmd == nr ) exit
                end do

            end if

            if ( inmd == nr .and. mkc(icol) == nc) exit

        end do RDWM

        stat = .true.

        system%state(istate)%molecule%nvib = mkc(icol)

        deallocate(g,aus,mkc)

        !   do i = 1, nc
        !     write(18,'(2x,a17,f7.2,a3)')'<vibration freq="', &
        !   	       system%state(istate)%molecule%normodes%vibration(i)%freq,'" >'
        !     do j = 1, system%state(istate)%molecule%structure%numat
        !       write(18,'(2x,3(f10.5,3x))') (system%state(istate)%molecule%normodes%vibration(i)%atom(j)%d(k),k=1,3)
        !     end do
        !     write(18,'(2x,a12)')'</vibration>'
        !   end do

	    ! Here we sort the modes according to their symmetry labels

        if (system%state(istate)%molecule%normodes%sort) call sort_modes_symmetry(istate)

        call scale_frequencies(istate)

        return

        ! Exit on error
222     stat = .false.

        return
    end subroutine read_normal_modes_gau

    ! ============================================================================ !

    subroutine read_normal_modes_gaunohp (iu, istate, imol, stat)

        integer, intent(in) :: iu, istate, imol
        logical, intent(out) :: stat

        integer, parameter :: MAXCOL = 9
        integer icol, nr, nc, inmd, i, j, k, ios, i1, i2, i3
        real(kind=dp), allocatable  ::  aus(:), g(:)
        integer, allocatable :: mkc(:)

        character (len=80) line
        character (len=15) str
        character (len=3) crj(1:3) 

        rewind iu

        ! Solo per specie non lineari!!
        nc = 3*system%state(istate)%molecule%structure%numat - 6
        nr = system%state(istate)%molecule%structure%numat

        allocate(g(1:MAXCOL), aus(1:MAXCOL), mkc(1:3))

        aus = ZERO

        stat = .false.

        ! --- Trova la stringa 'normal coordinates nell'input

        line = BLANK
        do while (index(line,'normal coordinates') == 0)
            read (iu,'(a80)',iostat=k,end=222)line
        end do

        RDWM:   do

            mkc = 0
            read (iu,'(a80)',end=222)line
            !read (line,*,iostat=ios) (mkc(i),i=1,MAXCOL)
            read (line,*,iostat=ios) (mkc(i),i=1,3)
            icol=count(mkc > 0)

            if (icol > 0) then

                ! Assegna le ID alle vibrazioni di ciscuna molecola.
                system%state(istate)%molecule%normodes%vibration(mkc(1):mkc(icol))%id = mkc(1:icol)

                ! --- Le frequenze che leggo sono positive e quindi le pongo uguali a -1
                ! --- per sapere se sono state lette o meno. Se la lettura e' corretta
                ! --- devono essere tutte positive

                g = -ONE
                do
                    read (iu,'(a80)',iostat=k)line
                    if (k /= 0) exit
                    read (line,*,iostat=ios) str, str, (g(i),i=1,icol)
                    if (ios == 0) then
                        system%state(istate)%molecule%normodes%vibration(mkc(1):mkc(icol))%freq = g(1:icol)
                    end if
                    ! Controlla che siano state effettivamente lette le frequenze.
                    if (count(g(1:icol) > ZERO) == icol) exit
                end do
            
                inmd = 0
                do
                    read (iu,'(a80)',iostat=k,end=222)line
                    if (k /= 0) stop "Cannot read Gaussian file."
                    read(line,*,iostat=ios)i2,i3,(aus(j),j=1,3*icol)
                    if (ios == 0) then
                        inmd = inmd + 1
                        do j = 1, icol
                            !system%state(istate)%molecule%normodes%vibration(mkc(j))%atom(i2)%elem = elmnts(mapel(i3))
                            system%state(istate)%molecule%normodes%vibration(mkc(j))%atom(i2)%elem = &
                            system%state(istate)%molecule%structure%atom(i2)%elem
                            system%state(istate)%molecule%normodes%vibration(mkc(j))%atom(i2)%d(1:3) = aus(3*(j-1)+1:3*j)
                        end do
                    end if
                    if ( inmd == nr ) exit
                end do

            end if

            if ( inmd == nr .and. mkc(icol) == nc) exit

        end do RDWM

        stat = .true.

        system%state(istate)%molecule%nvib = mkc(icol)

        deallocate(g,aus,mkc)

        call scale_frequencies(istate)

        return

        ! Exit on error
222     stat = .false.

        return
    end subroutine read_normal_modes_gaunohp

    !==============================================================================!

    subroutine read_normal_modes_tm (iu, istate, imol, stat)

        ! In uscita stat = .true. se i modi normali sono stati letti 
        !           stat = .false. se i modi normali non sono stati letti 
        integer, intent(in) :: iu, istate, imol
        logical, intent(out) :: stat

        integer, parameter :: MAXCOL = 5
        integer ICOL, NR, NROW, NVBR, NMAT, NTBR, i, j, k, l, ios
        integer i1, i2, i3, irv, iat, ico, icoord, inum, ivib
        real(kind=dp) freq
        real(kind=dp), allocatable  ::  aus(:)
       
        character (len=80) line
        character (len=15) str
        character (len=3) crj(1:3) 

        rewind iu

        ! Solo per specie non lineari!! 
        ! Devo trovare una soluzione per specie lineari.
        NVBR = 3*system%state(istate)%molecule%structure%numat - 6
        NMAT = system%state(istate)%molecule%structure%numat
        NR = 3*NMAT

        allocate(aus(1:NR))

        NROW = int(NR/MAXCOL)
        if (modulo(NR,MAXCOL) /= 0) NROW = NROW + 1
        aus = ZERO

        stat = .false.

        ! Assegna gli elementi alle singole componenti delle vibrazioni.
        ! potrebbe anche essere fatto in automatico dopo la lettura della struttura.
        do k = 1, NVBR
            system%state(istate)%molecule%normodes%vibration(k)%id = k
            do l = 1, NMAT
                system%state(istate)%molecule%normodes%vibration(k)%atom(l)%elem = &
                system%state(istate)%molecule%structure%atom(l)%elem
            end do
        end do

        ! --- Trova la stringa '$vibrational normal modes nell'input

        line = BLANK
        do while (index(line,'$vibrational normal modes') == 0)
            read (iu,'(a80)',iostat=k,end=222)line
        end do

        RDWM:   do i = 1, NR

            ICOL = 0
            NTBR = MAXCOL
            do j = 1, NROW
                if (j == NROW) then
                    NTBR = modulo(NR,MAXCOL)
                    if (NTBR == 0) NTBR = MAXCOL
                end if
                ICOL = ICOL + NTBR
                read(iu,'(i3,i2)',advance="no") icoord, inum
                read (iu,*,iostat=ios)(aus(k),k=(j-1)*MAXCOL+1,ICOL)
              !read (iu,*,iostat=ios) icoord, inum, (aus(k),k=(j-1)*MAXCOL+1,ICOL)
            end do

            iat = int((icoord-1)/3)+1
            ico = modulo(icoord-1,3)+1
            do k = 1, NVBR
                system%state(istate)%molecule%normodes%vibration(k)%atom(iat)%d(ico) = aus(k+6)
            end do

        end do RDWM

        deallocate(aus)

        ! --- Trova la stringa '$vibrational spectrum' nell'input
        ! in some output the two sections may be inverted thu we have to rewind the file
        rewind(iu)

        line = BLANK
        do while (index(line,'$vibrational spectrum') == 0)
            read (iu,'(a80)',iostat=k)line
        end do

        irv = 0
        do
            read (iu,'(a80)',iostat=k)line
            if (k /= 0) exit
            read (line,*,iostat=ios) ivib
            if (ios == 0 .and. (ivib > 6)) then
                read (line,*,iostat=ios) ivib, str, freq
                irv = irv + 1
                system%state(istate)%molecule%normodes%vibration(irv)%freq = freq
                system%state(istate)%molecule%normodes%vibration(irv)%symlab = adjustl(str)
            end if
            if (irv == NVBR) exit
        end do

        stat = .true.

        system%state(istate)%molecule%nvib = NVBR

        if (system%state(istate)%molecule%normodes%sort) call sort_modes_symmetry(istate)

        ! do i = 1, NVBR
        !   write(18,'(2x,a17,f7.2,a3)')'<vibration freq="', &
        ! 	       system%state(istate)%molecule%normodes%vibration(i)%freq,'" >'
        !   do j = 1, system%state(istate)%molecule%structure%numat
        !     write(18,'(2x,3(f10.5,3x))') (system%state(istate)%molecule%normodes%vibration(i)%atom(j)%d(k),k=1,3)
        !   end do
        !   write(18,'(2x,a12)')'</vibration>'
        ! end do

        call scale_frequencies(istate)

        return

        ! Exit on error
222     stat = .false.

        return
    end subroutine read_normal_modes_tm

    !==============================================================================!

    subroutine read_normal_modes_molden (iu, istate, imol, stat)

        ! Legge i modi normali da un input in formato MOLDEN:
        ! Per lo standard di MOLDEN vedere il sito web del programma.
        ! Ricordo solo che gli spostamenti sono espressi in Unità Atomiche (Bohr)

        ! In uscita stat = .true. se i modi normali sono stati letti 
        !           stat = .false. se i modi normali non sono stati letti 
        integer, intent(in) :: iu, istate, imol
        logical, intent(out) :: stat

        integer ICOL, NR, NROW, NVBR, NMAT, NTBR, i, j, k, l, ios
        integer i1, i2, i3, irv, iat, ico, icoord, inum, ivib, inmd
        real(kind=dp) freq
        real(kind=dp) aus(1:3)
       
        character (len=80) line
        character (len=15) str
        character (len=3) crj(1:3) 

        rewind iu

        ! Solo per specie non lineari!! 
        ! Devo trovare una soluzione per specie lineari.
        NVBR = 3*system%state(istate)%molecule%structure%numat - 6
        NMAT = system%state(istate)%molecule%structure%numat
        NR = 3*NMAT

        stat = .false.

        ! Assegna gli elementi alle singole componenti delle vibrazioni.
        ! potrebbe anche essere fatto in automatico dopo la lettura della struttura.
        do k = 1, NVBR
            system%state(istate)%molecule%normodes%vibration(k)%id = k
            do l = 1, NMAT
                system%state(istate)%molecule%normodes%vibration(k)%atom(l)%elem = &
                system%state(istate)%molecule%structure%atom(l)%elem
            end do
        end do

        ! --- Trova la stringa '[fr-norm-coord] ell'input

        line = BLANK
        do while (index(line,'[FR-NORM-COORD]') == 0)
            read (iu,'(a80)',iostat=k,end=222)line
        end do

        RDWM:   do i = 1, NVBR

            read (iu,'(a80)',iostat=k,end=222)line
            if (k /= 0) stop "Cannot read MOLDEN file."
            read(line,*,iostat=ios)str, i2
            if (ios == 0) then
                if (index(str,'vibration') /= 0) then
                    do j = 1, NMAT
                        read(iu,*) aus(1:3)
                        system%state(istate)%molecule%normodes%vibration(i)%atom(j)%d(1:3) = aus(1:3)
                    end do
                end if
            end if

        end do RDWM


        ! Since normal modes have been transformed to bohr they must be renormalized
        !do j = 1, NVBR
        !    ! Renormalize the j-th normal coordinate.
        !    norm = zero
        !    do k = 1, system%state(is)%molecule%structure%numat
        !      x(1:3) = system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3)
        !      norm = norm + dot_product(x,x)
        !    end do
        !    do k = 1, system%state(is)%molecule%structure%numat
        !      system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3) = &
        !      system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3)/sqrt(norm)
        !    end do
        !end do

        ! --- Trova la stringa '[FREQ]' nell'input

        rewind iu

        line = BLANK
        do while (index(line,'[FREQ]') == 0)
            read (iu,'(a80)',iostat=k)line
        end do

        irv = 0
        do
            read (iu,'(a80)',iostat=k)line
            if (k /= 0) exit
            irv = irv + 1
            read (line,*,iostat=ios) freq
            system%state(istate)%molecule%normodes%vibration(irv)%freq = freq
            if (irv == NVBR) exit
        end do

        stat = .true.

        system%state(istate)%molecule%nvib = NVBR

        call scale_frequencies(istate)

        return

        ! Exit on error
222     stat = .false.

        return
    end subroutine read_normal_modes_molden


    !==============================================================================!

    subroutine read_normal_modes_mopac (iu, istate, imol, stat)

        ! In uscita stat = .true. se i modi normali sono stati letti
        !           stat = .false. se i modi normali non sono stati letti
        integer, intent(in) :: iu, istate, imol
        logical, intent(out) :: stat

        integer, parameter :: MAXCOL = 8
        integer :: RCOL
        integer icol, nr, nc, inmd, i, j, k, ios, i1, i2, iat
        real(kind=dp), allocatable  ::  aus(:), g(:)
        integer, allocatable :: mkc(:)

        character (len=100) line
        character (len=15) str, str1, str2
        character (len=3) crj(1:3)

        rewind iu

        ! Solo per specie non lineari!!
        ! Devo trovare una soluzione per specie lineari.
        nc = 3*system%state(istate)%molecule%structure%numat - 6
        nr = 3*system%state(istate)%molecule%structure%numat

        allocate(g(1:MAXCOL), aus(1:MAXCOL), mkc(1:MAXCOL))

        aus = ZERO

        stat = .false.

        ! --- Trova la stringa 'normal coordinates nell'input

        line = BLANK
        do while (index(line,'MASS-WEIGHTED COORD') == 0)
            read (iu,'(a100)',iostat=k,end=222)line
        end do

        RDWM:   do

            RCOL = 8
            mkc = 0
            read (iu,'(a)',end=222)line
            do
                read (line,*,iostat=ios) str1, str2, (mkc(i),i=1,RCOL)
                if (RCOL == 1) then
                    !print *, 'rejected row'
                    exit
                end if
                if (ios == 0) exit
                if (ios /= 0) then
                    !print *, 'error line ', line
                    RCOL = RCOL - 1
                  !print *, 'maxcol ', rcol
                end if
            end do

            if (ios == 0) then

                icol=count(mkc > 0)

                if (icol > 0) then
                    !write(*,*) str1, str2, (mkc(i),i=1,RCOL)


                    ! Assegna le ID alle vibrazioni di ciscuna molecola.
                    system%state(istate)%molecule%normodes%vibration(mkc(1):mkc(icol))%id = mkc(1:icol)

                    ! --- Le frequenze che leggo sono positive e quindi le pongo uguali a -1
                    ! --- per sapere se sono state lette o meno. Se la lettura e' corretta
                    ! --- devono essere tutte positive

                    g = -ONE
                    do
                        read (iu,'(a100)',iostat=k)line
                        if (k /= 0) exit
                        read (line,*,iostat=ios) (g(i),i=1,icol)
                        if (ios == 0) then
                            !print *, icol, g(1:icol)
                            system%state(istate)%molecule%normodes%vibration(mkc(1):mkc(icol))%freq = g(1:icol)
                        end if
                        ! Controlla che siano state effettivamente lette le frequenze.
                        if (count(g(1:icol) > ZERO) == icol) exit
                    end do

                    inmd = 0
                    do
                        read (iu,'(a100)',iostat=k,end=222)line
                        !print *, 'reading displ ', line
                        if (k /= 0) stop "Cannot read line"
                        read(line,*,iostat=ios)i1,(aus(j),j=1,icol)
                        if (ios == 0) then
                            iat = int((i1-1)/3)+1
                            i2 = mod(i1-1,3)+1
                            !  print *, 'i1: ', i1, 'iat: ', iat, 'i2: ', i2
                            inmd = inmd + 1
                            do j = 1, icol
                                system%state(istate)%molecule%normodes%vibration(mkc(j))%atom(iat)%elem = &
                                system%state(istate)%molecule%structure%atom(iat)%elem
                                system%state(istate)%molecule%normodes%vibration(mkc(j))%atom(iat)%d(i2) = aus(j)
                            end do
                        end if
                        if ( inmd == nr ) exit
                    end do

                end if

                if ( inmd == nr .and. mkc(icol) == nc) exit

            end if

        end do RDWM

        stat = .true.

        system%state(istate)%molecule%nvib = mkc(icol)

        deallocate(g,aus,mkc)

        !   do i = 1, nc
        !     write(18,'(2x,a17,f7.2,a3)')'<vibration freq="', &
        !   	       system%state(istate)%molecule%normodes%vibration(i)%freq,'" >'
        !     do j = 1, system%state(istate)%molecule%structure%numat
        !       write(18,'(2x,3(f10.5,3x))') (system%state(istate)%molecule%normodes%vibration(i)%atom(j)%d(k),k=1,3)
        !     end do
        !     write(18,'(2x,a12)')'</vibration>'
        !   end do

        call scale_frequencies(istate)

        return

        ! Exit on error
222     stat = .false.

        return
    end subroutine read_normal_modes_mopac

    !==============================================================================!

    subroutine read_normal_modes_gamess (iu, istate, imol, stat)

        ! In uscita stat = .true. se i modi normali sono stati letti
        !           stat = .false. se i modi normali non sono stati letti
        integer, intent(in) :: iu, istate, imol
        logical, intent(out) :: stat

        integer NVBR, NUMAT, inmd, i, j, k, ios, idk, iat, ivib
        real(kind=dp)  ::  disp(1:3), freq, orth, av(1:3), bv(1:3)

        character (len=100) line
        character (len=5) str

        rewind iu

        ! Solo per specie non lineari!!
        ! Devo trovare una soluzione per specie lineari.
        NVBR = 3*system%state(istate)%molecule%structure%numat - 6
        NUMAT = system%state(istate)%molecule%structure%numat
        system%state(istate)%molecule%nvib = NVBR

        stat = .false.

        ! --- Trova la stringa 'normal coordinates nell'input

        line = BLANK

        idk = 0
        do
                if (idk == NVBR) exit
                read (iu,'(a100)',iostat=k)line
                if (k /= 0) stop "Cannot read line"
                if (index(line,'MODE') /= 0) then
                    idk = idk + 1
                    read(line,*,iostat=ios)str,ivib, freq
                    !print *, 'reading ', str, ivib, freq
                    do iat = 1, NUMAT
                            !read(iu,'(a100)') line
                            !print *, 'disp line ', line
                            read(iu,*) disp(1:3)
                            system%state(istate)%molecule%normodes%vibration(idk)%id = idk
                            system%state(istate)%molecule%normodes%vibration(idk)%freq = freq
                            system%state(istate)%molecule%normodes%vibration(idk)%atom(iat)%elem = &
                            system%state(istate)%molecule%structure%atom(iat)%elem
                            system%state(istate)%molecule%normodes%vibration(idk)%atom(iat)%d(1:3) = disp(1:3)
                    end do
                end if
        end do

        ! check ortho
!        do i = 1, NVBR
!          do j = 1, NVBR
!            orth = 0.0_dp
! 	    do k = 1, NUMAT 
!	      av(1:3) = system%state(istate)%molecule%normodes%vibration(i)%atom(k)%d(1:3)
!	      bv(1:3) = system%state(istate)%molecule%normodes%vibration(j)%atom(k)%d(1:3)
!	      orth = orth + dot_product(av,bv)
!            end do
!            print *, 'ORTHO: ', i, j,  orth
!          end do
!        end do

        stat = .true.

        call scale_frequencies(istate)

        return

        ! Exit on error
222     stat = .false.

        return
    end subroutine read_normal_modes_gamess

! ============================================================================ !

        SUBROUTINE read_normal_modes_demon_nohp (iu, istate, imol, stat)

        integer, intent(in) :: iu, istate, imol
        logical, intent(out) :: stat

        integer, parameter :: MAXCOL = 9
        integer icol, nr, nc, nf, inmd, i, j, k, ios, i1, i2, i3, iam
        real(kind=dp) :: freq, dx, dy, dz
	integer :: mkc
        
        character (len=80) line
        character (len=10) str 
	character (len=3) asym
        character (len=3) crj(1:3)

        rewind iu

        ! Solo per specie non lineari!!
        nf = 6
        if (system%state(istate)%molecule%structure%numat == 2) nf = 5
        nc = 3*system%state(istate)%molecule%structure%numat - nf
        nr = system%state(istate)%molecule%structure%numat

        stat = .false.

! --- Trova la stringa 'normal coordinates nell'input

        line = BLANK
        do while (index(line,'VIBRATIONS') == 0)
          read (iu,'(a80)',iostat=k)line
        end do

            mkc=1

RDWM:    do
	    ! --- Trova la stringa 'MODE' nell'input.
	    line = BLANK
	    do while (index(line,'MODE:') == 0)
	      read (iu,'(a80)',iostat=k)line
	    end do

		system%state(istate)%molecule%normodes%vibration(mkc)%id = mkc

	    freq = -ONE

	    read (iu,'(a80)',iostat=k)line
	    read (line,*, iostat=ios) str, freq, str
	
	      if (ios == 0) then
	         system%state(istate)%molecule%normodes%vibration(mkc)%freq = freq 
	      end if

              ! Controlla che siano state effettivamente lette le frequenze.
              if (freq < 0) exit
	
	      iam=0
	    do
	      read (iu,'(a80)',iostat=k, end=222)line
               	if (k /= 0) stop "Cannot read deMon file."
	      read (line,*, iostat=ios)  asym , dx, dy, dz
	      if (ios == 0) then
	        iam = iam + 1
                !  qstate(istate)%molecule(imol)%normodes%vibration(mkc)%atom(iam)%elem%Sym = asym 
                  system%state(istate)%molecule%normodes%vibration(mkc)%atom(iam)%elem = elmnts(mapel(asym))
                  system%state(istate)%molecule%normodes%vibration(mkc)%atom(iam)%d(1:3) = (/dx, dy, dz/)
 	      end if
              
              if  ( iam == nr ) exit
	    end do
	  
	if ( iam == nr .and. mkc == nc) exit
	mkc = mkc+1
	end do RDWM 

        stat = .true.

        system%state(istate)%molecule%nvib = mkc

222     if (iam == 0) stat = .false.


        return
        END SUBROUTINE read_normal_modes_demon_nohp


!=====================================================================================

    !Sort modes according to their symmetry is it is available
    subroutine sort_modes_symmetry(istate)

    integer, intent(in) :: istate
    integer :: nvib, i
    type(vibration_t) :: vtmp
    logical exch

    !print *, 'sorting modes ofstate: ', istate

    nvib = system%state(istate)%molecule%nvib

    exch = .true.
        do while (exch)
            exch = .false.
            do i = 1, nvib-1
            ! compare labels
            if (system%state(istate)%molecule%normodes%vibration(i)%symlab.gt.&
                system%state(istate)%molecule%normodes%vibration(i+1)%symlab) then
                vtmp = system%state(istate)%molecule%normodes%vibration(i)
                system%state(istate)%molecule%normodes%vibration(i) = &
                system%state(istate)%molecule%normodes%vibration(i+1)
                system%state(istate)%molecule%normodes%vibration(i+1) = vtmp
                !print *, 'swap ', i,i+1
                exch = .true.
            end if
            end do
        end do

    return
    end subroutine sort_modes_symmetry

    !=====================================================================================

    subroutine scale_frequencies(istate)

    integer, intent(in) :: istate
    real(kind=dp) :: fscale

    fscale = system%state(istate)%molecule%normodes%fscale

    write(fout,'(2x,a,i3,a,f8.4)') 'Scaling frequencies of state ', istate, ' by a factor ', fscale
    system%state(istate)%molecule%normodes%vibration(:)%freq = &
    system%state(istate)%molecule%normodes%vibration(:)%freq * fscale

    return
    end subroutine  scale_frequencies

end module iomodes

