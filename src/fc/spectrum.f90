module spectrum

    use fcint_h
    use matfun, only : vecsort

    implicit none

    private

    public :: fcspec

contains

    subroutine fcspec(acs_initial,acs_final,fcjob,gid)

        !-------------------------------------------------------------------------------------!
        ! Per ora questa subroutine funziona solo se nell'input viene usato un solo           !
        ! gruppo di vibrazioni. Per usare diversi gruppi di vibrazioni dovrei specificare un file !
        ! per ogni gruppo e questo richiede una modifica del formato di input.                !
        !-------------------------------------------------------------------------------------!

        type(activespace_t), intent(in) :: acs_initial, acs_final
        type(fc_t), intent(in) :: fcjob
        integer, intent(in) :: gid
        ! character(len=*), intent(in) :: molid

        integer, allocatable :: ivib(:), zmax(:), z(:)
        integer :: i, j, k, r, is, ix, zk, rmu, npt, NFC, i2, i3, ia, ntot
        integer :: nmt, nmd, initial_vibrations_number, total_vibrations_number
        integer :: mu, nexvib

        real(kind=dp) :: efc, Ej, EM, Emin, Emax, Temp, fwhm, bw, dex, ebra, eket, coeff
        real(kind=dp), allocatable :: freq_final(:), freq_initial(:), spec(:)

        nexvib = acs_initial%nact + acs_final%nact
        nmd = fcjob%group(gid)%nvib
        nmt = 2*nmd
  
        if (acs_initial%nact > 0 .and. acs_final%nact > 0) NFC = product(acs_initial%nqmax+1)*product(acs_final%nqmax+1)
        if (acs_initial%nact == 0 .and. acs_final%nact > 0) NFC = product(acs_final%nqmax+1)
        if (acs_initial%nact > 0 .and. acs_final%nact == 0) NFC = product(acs_initial%nqmax+1)
        if (acs_initial%nact == 0 .and. acs_final%nact == 0) NFC = 1
 
        allocate(z(1:nmt),zmax(1:nmt));   zmax = 0

        if (nexvib > 0) allocate (ivib(1:nexvib))
        ivib = 0
  
        if (acs_initial%nact > 0) then
            ivib(1:acs_initial%nact) = acs_initial%vibid ! Contiene le ID delle vINITIAL_STATEzioni dello spazio attivo.
            zmax(ivib(1:acs_initial%nact)) = acs_initial%nqmax
        end if

        if (acs_final%nact > 0) then
            ivib(acs_initial%nact+1:nexvib) = acs_final%vibid + int(nmt/2) ! Contiene le ID delle vINITIAL_STATEzioni dello spazio attivo.
            zmax(ivib(acs_initial%nact+1:nexvib)) = acs_final%nqmax
        end if
 
        allocate (freq_initial(1:fcjob%group(gid)%nvib), freq_final(1:fcjob%group(gid)%nvib))
        freq_initial = system%state(INITIAL_STATE)%molecule%normodes%vibration(fcjob%group(gid)%incvib%id(:))%freq
        freq_final = system%state(FINAL_STATE)%molecule%normodes%vibration(fcjob%group(gid)%incvib%id(:))%freq

        ! apre il file su cui scrivere lo spettro
        i2 = int(gid/10)
        i3 = gid-i2*10
        i2 = i2+ichar('0')
        i3 = i3+ichar('0')

        if (fcjob%spectrum%form == "formatted") then
            open(unit=iuspec,file=fcjob%spectrum%file(1:len_trim(fcjob%spectrum%file))//char(i2)//char(i3), &
            status="unknown")
        else
            open(unit=iuspec,file=fcjob%spectrum%file(1:len_trim(fcjob%spectrum%file))//char(i2)//char(i3), &
            form="unformatted",status="unknown")
        end if
  
        ! -- Definisce i parametri della convoluzione.
        Emin = fcjob%spectrum%emin
        Emax = fcjob%spectrum%emax
        if (Emax < 0) Emax = Energy(freq_final,zmax(nmd+1:nmt)) + system%state(FINAL_STATE)%energy
        EM = fcjob%spectrum%dE
        npt = int((Emax - Emin)/EM) + 1
        fwhm = fcjob%spectrum%fwhm
        Temp = fcjob%spectrum%Temp

        write(fout,*) ' Simulated Spectrum written in file ', fcjob%spectrum%file(1:len_trim(fcjob%spectrum%file))
        write(fout,'(2x,a,f10.2)') 'Emin: ', Emin
        write(fout,'(2x,a,f10.2)') 'Emax: ', Emax
        write(fout,'(2x,a,f10.2)') 'DE: ', EM
        write(fout,'(2x,a,i7)') 'NPT: ', NPT
        write(fout,'(2x,a,f10.3)') 'FWHM: ', fwhm
        write(fout,'(2x,a,f10.3)') 'Temp: ', Temp
  
        ! --- Stampa lo spettro FC ---

        initial_vibrations_number = count(ivib <= nmd)
        total_vibrations_number = size(ivib)
  
        z = 0
        if (fwhm < zero) then
            !-------------------------------------------
            ! Scrive i FC senza fare alcuna convoluzione.
            ! Scrive i FC pesati con una Boltzmann sull'energia dello stato del ket
            ! (quello di partenza).
            ntot = 0
            do r = 1, NFC
                if (abs(MDFC(r)) >= fcjob%spectrum%tol) then
                    eket = Energy(freq_final,z(nmd+1:nmt))+ system%state(FINAL_STATE)%energy
                    ebra = Energy(freq_initial,z(1:nmd)) + system%state(INITIAL_STATE)%energy
                    dex = ebra/(bkcm*Temp)
                    efc = eket - ebra
                    if (efc >= emin .and. efc <= emax ) then
                        ntot = ntot + 1
                        coeff = exp(-dex)
                        if (fcjob%spectrum%type == "abs") coeff = efc*exp(-dex)
                        if (fcjob%spectrum%type == "ems") coeff = (efc**3)*exp(-dex)
						!print *, 'RR ', MDFC(r)
                        if (abs(MDFC(r)) <= 1.0D-90) then
                            ! or we have formatting problems...
                            write(iuspec,'(2x,f10.2,E14.7)',advance='no') efc, zero
                        else
                            if (fcjob%fcht) then
                                write(iuspec,'(2x,f10.2,E14.7,2x,E14.7)',advance='no') efc, (MDFC(r)**2)*coeff, (FCHT(r))*coeff
                            else
                                write(iuspec,'(2x,f10.2,E14.7)',advance='no') efc, (MDFC(r)**2)!*coeff
                            end if
                        end if
                            ! write(iuspec,'(2x,f10.2,g14.5,g14.5,f10.6)',advance='no') efc, (MDFC(r)**2)*coeff,(MDFC(r)**2), exp(-dex)
                            ! writes transition assignment
                        write(iuspec,'(2x,a3)',advance='no')'  <'
                        do k = 1, initial_vibrations_number
                            write(iuspec,'(i3)', advance='no') z(ivib(k))
                        end do
                        write(iuspec,'(a1)',advance='no')'|'
                        do k = 1+initial_vibrations_number, total_vibrations_number
                            write(iuspec,'(i3)', advance='no') z(ivib(k))
                        end do
                        write(iuspec,'(a3)')'> '
                    end if
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

        else ! (se fwhm > 0.0)
            ! Scrive i FC facendo la convoluzione con una lorentziana in base ai parametri
            ! specificati in input in <spectrum>.
            allocate(spec(1:npt))
            spec = 0.D0
            FC:   do r = 1, NFC

                if (abs(MDFC(r)) >= fcjob%spectrum%tol) then

                    efc = Energy(freq_final,z(nmd+1:nmt))-Energy(freq_initial,z(1:nmd)) + & 
                          system%state(FINAL_STATE)%energy - system%state(INITIAL_STATE)%energy
                    dex = Energy(freq_initial,z(1:nmd))/(bkcm*Temp)
                    bw = exp(-dex)

                    do j = 1, npt
                        Ej = Emin+(j-1)*EM
                        if (abs(Ej - efc) < 20*fwhm) then
                            spec(j) = spec(j) + bw * (MDFC(r)**2)*(1/pi)*(0.5*fwhm/(((Ej - efc)**2 + ((0.5*fwhm)**2))))
                        end if
                    end do

                end if

                do i = 1, nexvib
                    mu = ivib(i)
                    if (z(mu) < zmax(mu)) then
                        z(mu) = z(mu) + 1
                        if (i > 1) z(1:ivib(i-1)) = 0
                        exit
                    end if
                end do

            end do FC

            do j = 1, npt
                write(iuspec,*) Emin+(j-1)*EM, spec(j)
            end do

            deallocate (spec)
        end if
  
        deallocate (freq_final, freq_initial)
        deallocate (z, zmax)
        if (allocated(ivib)) deallocate(ivib)

        ! scrive sul file il numero di FC selezionati
        !rewind iuspec
        write(iuspec,'(a,1x,f8.2,a,f10.2,a,i7,a,f5.2,a)')'#<shape emin="',emin,'" emax="',emax,'" nfc="',ntot,'" w="',fwhm,'" />'
        ! chiude il file su cui ha scritto lo spettro
        close(iuspec)
  
        return
    end subroutine fcspec
  
end module spectrum
