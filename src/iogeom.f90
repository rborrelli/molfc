module iogeom

    use parameters
    use isotopes
    use errors
    use xmvar, only : system

    implicit none

    public :: read_geometry_gau, read_geometry_tm, read_geometry_xyz, read_geometry_molden, &
    read_geometry_mopac, read_geometry_demon
contains
  
    !=======================================================================!

    subroutine read_geometry_gau(iu, istate, imol, iam)
  
        integer, intent(in) :: iu, istate, imol
        integer, intent(out) :: iam
        integer ios, cn, anum, at, numat, k, ngeo
        real(kind=dp) ax, ay, az
        character(len=80) line
        logical next_geometry
  
        numat = system%state(istate)%molecule%structure%numat
  
        rewind iu
  
        ! Legge tutte le geomentrie presenti nel file e ne memorizza solo l'ultima.
  
        next_geometry = .true.
        ngeo = 0
        do
            ! --- Trova la stringa 'Standard orientation' nell'input.
            line = BLANK
            do while (index(line,'Standard orientation') == 0)
                read (iu,'(a80)',iostat=k,end=222)line
            end do

            iam = 0
            do
                read (iu,'(a80)',iostat=k,end=222)line
                read (line,*, iostat=ios) cn, anum, at, ax, ay, az
                if (ios == 0) then
                    iam = iam + 1
                    system%state(istate)%molecule%structure%atom(iam)%elem = elmnts(mapel(anum))
                    system%state(istate)%molecule%structure%atom(iam)%coord(1:3) = (/ax, ay, az/)
                end if
                if (iam == numat) then
                    ngeo = ngeo + 1
                    exit
                end if
            end do

        end do

        ! Se non ha letto nessuna struttura esce con errore.
222     if (ngeo == 0) stop "Cannot read geometry from input file"

        return
    end subroutine read_geometry_gau

    !=======================================================================!

    subroutine read_geometry_tm(iu, istate, imol, iam)
 
        ! Legge la geometria dal file coord di Turbomole o dal file
        integer, intent(in) :: iu, istate, imol
        integer, intent(out) :: iam
        integer ios, cn, at, numat, k, ngeo, ierr
        real(kind=dp) ax, ay, az
        character(len=80) line
        character(len=2) asym
        logical next_geometry
  
        numat = system%state(istate)%molecule%structure%numat
  
        rewind iu
  
        ! --- Trova la stringa '$coord' nell'input.
        line = BLANK
        do while (index(line,'$coord') == 0)
            read (iu,'(a80)',iostat=k,end=222)line
        end do
  
        iam = 0
        do
            read (iu,'(a80)',iostat=k,end=222)line
            read (line,*, iostat=ios) ax, ay, az, asym
            if (ios == 0) then
                iam = iam + 1
                system%state(istate)%molecule%structure%atom(iam)%elem = elmnts(mapel(asym))
                ! Geometry is automatically converter from bohr to angstrom.
                system%state(istate)%molecule%structure%atom(iam)%coord(1:3) = (/ax, ay, az/)*bohr
            end if
            if (iam == numat) return
        end do

        ! Se non ha letto nessuna struttura esce con errore.
222     ierr = error(0,"Cannot read geometry from TURBOMOLE input file")

        return
    end subroutine read_geometry_tm

    !=======================================================================!

    subroutine read_geometry_xyz(iu, istate, imol, iam)
 
        ! Legge la geometria dal file coord di Turbomole o dal file
        integer, intent(in) :: iu, istate, imol
        integer, intent(out) :: iam
        integer ios, cn, at, numat, k, ngeo, ierr
        real(kind=dp) ax, ay, az
        character(len=80) line
        character(len=200) title
        character(len=2) asym
        logical next_geometry
  
        rewind iu
  
        !skip beginning empty lines
        line = BLANK
        do while (len_trim(line) == 0)
            read (iu,'(a80)',iostat=k,end=222)line
        end do

        read(line,*,iostat=k) numat
        if (numat /= system%state(istate)%molecule%structure%numat) ierr = error(0,'Cannot read atom number in XYZ file.')
        ! read title
        read(iu,'(a200)',iostat=k,end=222)title
        iam = 0
        do
            read (iu,'(a80)',iostat=k,end=222)line
            read (line,*, iostat=ios) asym, ax, ay, az
            if (ios == 0) then
                iam = iam + 1
                system%state(istate)%molecule%structure%atom(iam)%elem = elmnts(mapel(asym))
                system%state(istate)%molecule%structure%atom(iam)%coord(1:3) = (/ax, ay, az/)
            end if
            if (iam == numat) return
        end do

        ! Se non ha letto nessuna struttura esce con errore.
222     ierr = error(0,"Cannot read geometry from XYZ type input file")

        return
    end subroutine read_geometry_xyz

    !=======================================================================!

    subroutine read_geometry_molden(iu, istate, imol, iam)
 
        ! Legge la geometria da un input in formato MOLDEN:
        ! Per lo standard di MOLDEN vedere il sito web del programma.
        ! Ricordo solo che le coordinate in [FR-COORD] sono espresse in Unità Atomiche (Bohr)

        integer, intent(in) :: iu, istate, imol
        integer, intent(out) :: iam
        integer ios, cn, at, numat, k, ngeo, ierr
        real(kind=dp) ax, ay, az
        character(len=80) line
        character(len=2) asym
        logical next_geometry
  
        rewind iu
  
        numat = system%state(istate)%molecule%structure%numat

        line = BLANK
        do while (index(line,'[FR-COORD]') == 0)
            read (iu,'(a80)',iostat=k,end=222)line
        end do

        iam = 0
        do
            read (iu,'(a80)',iostat=k,end=222)line
            read (line,*, iostat=ios) asym, ax, ay, az
            if (ios == 0) then
                iam = iam + 1
                system%state(istate)%molecule%structure%atom(iam)%elem = elmnts(mapel(asym))
                system%state(istate)%molecule%structure%atom(iam)%coord(1:3) = (/ax, ay, az/)*bohr
            end if
            if (iam == numat) return
        end do

        ! Se non ha letto nessuna struttura esce con errore.
222     ierr = error(0,"Cannot read geometry from MOLDEN type input file")

        return
    end subroutine read_geometry_molden

    !=======================================================================!

    subroutine read_geometry_mopac(iu, istate, imol, iam)

        integer, intent(in) :: iu, istate, imol
        integer, intent(out) :: iam
        integer ios, cn, anum, at, numat, k, ngeo
        real(kind=dp) ax, ay, az
        character(len=80) line
        character(len=2) asym
        logical next_geometry

        numat = system%state(istate)%molecule%structure%numat

        rewind iu

        ! Legge tutte le geomentrie presenti nel file e ne memorizza solo l'ultima.

        next_geometry = .true.
        ngeo = 0
        do
            ! --- Trova la stringa 'Standard orientation' nell'input.
            line = BLANK
            do while (index(line,'ORIENTATION') == 0)
                read (iu,'(a80)',iostat=k,end=222)line
            end do

            iam = 0
            do
                read (iu,'(a80)',iostat=k,end=222)line
                read (line,*, iostat=ios) cn, asym, ax, ay, az
                if (ios == 0) then
                    iam = iam + 1
                    system%state(istate)%molecule%structure%atom(iam)%elem = elmnts(mapel(asym))
                    system%state(istate)%molecule%structure%atom(iam)%coord(1:3) = (/ax, ay, az/)
                end if
                if (iam == numat) then
                    ngeo = ngeo + 1
                    exit
                end if
            end do

        end do

        ! Se non ha letto nessuna struttura esce con errore.
222     if (ngeo == 0) stop "Cannot read geometry from input file"

        return
    end subroutine read_geometry_mopac

    !=======================================================================!

  SUBROUTINE read_geometry_demon(iu, istate, imol, iam)

  integer, intent(in) :: iu, istate, imol
  integer, intent(out) :: iam
  integer ios, cn, anum, at, numat, k, ngeo,  amass
  real(kind=dp) ax, ay, az
  character(len=80) line
 character(len=5) asym
  logical next_geometry

  numat = system%state(istate)%molecule%structure%numat

  rewind iu

  ! Legge tutte le geomentrie presenti nel file e ne memorizza solo l'ultima.

  next_geometry = .true.
  ngeo = 0
  do
    ! --- Trova la stringa 'Standard orientation' nell'input.
    line = BLANK
    do while (index(line,'STANDARD ORIENTATION IN ANGSTROM') == 0)
      read (iu,'(a80)',iostat=k,end=222)line
    end do

    iam = 0
    do
      read (iu,'(a80)',iostat=k)line
      read (line,*, iostat=ios) cn, asym, ax, ay, az, anum , amass
      if (ios == 0) then
        iam = iam + 1
        system%state(istate)%molecule%structure%atom(iam)%elem = elmnts(mapel(anum))
        system%state(istate)%molecule%structure%atom(iam)%coord(1:3) = (/ax, ay, az/)
      end if
      if (iam == numat) then
        ngeo = ngeo + 1
        exit
      end if
    end do

  end do

  ! Se non ha letto nessuna struttura esce con errore.
222 if (ngeo == 0) stop "Cannot read geometry from input file"

  return
  END SUBROUTINE read_geometry_demon


end module iogeom
