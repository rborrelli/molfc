module FCIO

    use fcint_h

    implicit none

    private

    public :: fcint_write

contains

    !==============================================================================!

    subroutine fcint_write(acs_initial,acs_final,fcjob,gid)

        type(activespace_t), intent(in) :: acs_initial, acs_final
        type(fc_t), intent(in) :: fcjob
        integer, intent(in) :: gid

        integer, allocatable :: ivib(:), zmax(:), z(:)
        integer :: i, j, k, iss, r, kx, jj, is, ix, zk, rmu
        integer :: nmt, nmd, NCBRA, NSZIV, NFC
        integer :: mu, nexvib, iscratch

        real(kind=dp) :: coeff, fcsum, enfc, ebra, eket, WFC
        real(kind=dp), allocatable :: fbra(:), fket(:)


        nmd = fcjob%group(gid)%nvib
        nmt = 2*nmd

        nexvib = acs_initial%nact + acs_final%nact
        if (acs_initial%nact > 0 .and. acs_final%nact > 0) NFC = product(acs_initial%nqmax+1)*product(acs_final%nqmax+1)
        if (acs_initial%nact == 0 .and. acs_final%nact > 0) NFC = product(acs_final%nqmax+1)
        if (acs_initial%nact > 0 .and. acs_final%nact == 0) NFC = product(acs_initial%nqmax+1)
        if (acs_initial%nact == 0 .and. acs_final%nact == 0) NFC = 1

        allocate (z(1:nmt), zmax(1:nmt))

        z = 0; zmax = 0

        if (nexvib > 0) allocate (ivib(1:nexvib))

        if (.not.allocated(ivib)) then
            write(fout,'(2x,a)',advance='no')'< 0 | 0 >'
            enfc = system%state(FINAL_STATE)%energy - system%state(INITIAL_STATE)%energy
            write(fout,'(2x,g14.6,4x,f10.2,4x)') MDFC(1), enfc
            return
        end if

        ivib = 0
        if (acs_initial%nact > 0) then
            ivib(1:acs_initial%nact) = acs_initial%vibid  ! Contiene le ID delle vibrazioni dello spazio attivo.
            zmax(ivib(1:acs_initial%nact)) = acs_initial%nqmax
        end if

        if (acs_final%nact > 0) then
            ivib(acs_initial%nact+1:nexvib) = acs_final%vibid + nmd ! Contiene le ID delle vibrazioni dello spazio attivo.
            zmax(ivib(acs_initial%nact+1:nexvib)) = acs_final%nqmax
        end if

        ! --- Stampa i FC ---

        write(fout,299)
  
        allocate (fbra(1:fcjob%group(gid)%nvib), fket(1:fcjob%group(gid)%nvib))
        fbra = system%state(INITIAL_STATE)%molecule%normodes%vibration(fcjob%group(gid)%incvib%id(:))%freq
        fket = system%state(FINAL_STATE)%molecule%normodes%vibration(fcjob%group(gid)%incvib%id(:))%freq

        NCBRA = count(ivib <= nmd)
        NSZIV = size(ivib)
  
        iscratch = 44
        open(iscratch,file="FCI.dat",status="unknown")
        do r = 1, NFC

            if (fcjob%fcht) then
                WFC = FCHT(r)
                else
                WFC = MDFC(r)
            end if

            write(iscratch,*) WFC

            if (abs(WFC) >= fcjob%ftol) then
                write(fout,'(2x,a3)',advance='no')'  <'
                do k = 1, NCBRA
                    write(fout,'(i3)', advance='no') z(ivib(k))
                end do
                write(fout,'(a1)',advance='no')'|'
                do k = 1+NCBRA, NSZIV
                    write(fout,'(i3)', advance='no') z(ivib(k))
                end do
                write(fout,'(a3)',advance='no')'> '
                eket = Energy(fket,z(nmd+1:nmt))+ system%state(FINAL_STATE)%energy
                ebra = Energy(fbra,z(1:nmd)) + system%state(INITIAL_STATE)%energy
                enfc = eket - ebra
                write(fout,'(2x,g14.6,4x,f10.2,4x)') WFC, enfc
            end if

            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0
                    exit
                end if
            end do

        end do

        write(fout,300) sum(MDFC**2)
        write(fout,301) sum(MDFC**2) - MDFC(1)**2
  
        deallocate (fbra, fket)
        deallocate (z,zmax,ivib)

        include 'formats'

        return
    end subroutine fcint_write
  
end module FCIO
