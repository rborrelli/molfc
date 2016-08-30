module FCintegrals

    use fcint_h
    use fc_h, only : M, Q
    use fc, only : d
    use matfun, only : vsort

    implicit none

    private

    public :: fcint_init, fcint, fchtint, fcint_free, fcclasses, fchtclasses, &
               fchtint_old, fcht00, fcclasses_boltz, berkowitz, fcclasses_twostate

    integer :: NFC, NMT, nexvib
    integer, allocatable :: ivib(:), zmax(:), s(:)
    logical, save :: fc00written = .false.
    logical, save :: fcht00written = .false.
    real(kind=dp) :: sumfc, sumfcht, MDFCHT00
    real(kind=dp) :: FCHT0(1:3)

contains

    !==============================================================================!

    subroutine fcint_init (acs_initial,acs_final,FC00)

        !--------------------------------------------------------------------------------+
        ! ON INPUT                                                                       |
        ! ACS: is a structure containing NEX, IVIB and ZMAX, where                       |
        !
        ! NEX     Number of excited vibrations                                           |
        ! IVIB    Array di dimensioni NEX, che contiene le ID delle vibrazioni eccitate. |
        !         Esempio: modi eccitati = 1, 7 18                                       |
        !         IVIB = (/1, 7, 18/)                                                    |
        ! ZMAX    Maximum quanta on excited vibrations                                   |
        !--------------------------------------------------------------------------------|
        ! The matrices M and Q, used in the recurrence formulae are made accessible via  |
        ! a "use" statement. They are defined in the module fc_h.                        |
        !--------------------------------------------------------------------------------+

        ! N.B.: Le vibrazioni definite nello spazio attivo (acs),  come lette in input non sono in genere ordinate,
        ! a meno che  l'utente non le abbia ordinate a mano nell'input. D'altro canto e' necessario che vi sia
        ! corrispondenza tra l'ordine delle vibrazioni in M, Q ed s(i) e quello in ivib quindi
        ! IVIB deve essere riordinato.  (Il sorting di ivib e' fondamentale altrimenti si dovrebbero definire
        ! gli s(i) in altro modo.)
        ! Nella nuova versione di set_active_space() le vibid sono gia' riordinate
        type(activespace_t), intent(in) :: acs_initial, acs_final
        ! FC00    Franck-Condon <0|0>                                                    |
        real(kind=dp), intent(in) :: FC00

        ! local variables
        integer :: i, ierr, alloc_err, nfci, nfcf

        print *, 'FC00 ', FC00
        sumfc = FC00**2

        nexvib = acs_initial%nact + acs_final%nact
        !if (acs_initial%nact > 0 .and. acs_final%nact > 0) NFC = product(acs_initial%nqmax+1)*product(acs_final%nqmax+1)
        !if (acs_initial%nact == 0 .and. acs_final%nact > 0) NFC = product(acs_final%nqmax+1)
        !if (acs_initial%nact > 0 .and. acs_final%nact == 0) NFC = product(acs_initial%nqmax+1)
        !if (acs_initial%nact == 0 .and. acs_final%nact == 0) NFC = 1

        nfci = 1; nfcf = 1
        if (acs_initial%nact > 0) nfci = product(acs_initial%nqmax+1)
        if (acs_final%nact > 0) nfcf = product(acs_final%nqmax+1)
        NFC = nfci*nfcf

        allocate (MDFC(1:NFC),stat=alloc_err)
        if (alloc_err /= 0) ierr = error(0,'Cannot allocate FC integrals array.')
        MDFC = zero;
        !-------------------------------------------!
        ! FC00 e' il FC tra i ground states (<0|0>) !
        !-------------------------------------------!
        MDFC(1) = FC00

        ! The overall dimensionality of the problem which is twice the number of
        ! vibrations.
        nmt = size(M,2)

        allocate (zmax(1:nmt), s(1:nmt))
        zmax = 0

        if (nexvib > 0) then
            allocate (ivib(1:nexvib))
            ivib = 0
        else
            return ! no active space. Probably automatic calculation using fcclasses.
        end if

        if (acs_initial%nact > 0) then
            ivib(1:acs_initial%nact) = acs_initial%vibid ! Contiene le ID delle vibrazioni dello spazio attivo.
            zmax(ivib(1:acs_initial%nact)) = acs_initial%nqmax
        end if

        if (acs_final%nact > 0) then
            ivib(acs_initial%nact+1:nexvib) = acs_final%vibid + int(nmt/2) ! Contiene le ID delle vibrazioni dello spazio attivo.
            zmax(ivib(acs_initial%nact+1:nexvib)) = acs_final%nqmax
        end if

        s(1) = 1
        do i = 2, nmt
            s(i) = product(zmax(1:i-1)+1)
        end do

        return
    end subroutine fcint_init

    !==============================================================================!

    subroutine fcint ()

        ! Local variables
        real(kind=dp) :: coeff

        integer :: mu   ! Index of the excited vibration.
        integer :: i, j, k, iss, r, kx, zk, rmu

        ! z is used to store the 2N quantum numbers of each FC integrals
        integer, allocatable :: z(:)

        ! If the active space is empty just exit
        if (nexvib == 0) return

        !--------------------------!
        !     Inizio calcolo FC    !
        !--------------------------!
        allocate (z(1:nmt))
        z = 0
        do r = 2, NFC

            !--------------------------------------!
            ! generate next set of quantum numbers !
            !--------------------------------------!
            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0
                    exit
                end if
            end do

            rmu = r - s(mu)
    
            MDFC(r) = Q(mu) * MDFC(rmu)

            do k = 1, nexvib
                kx = ivib(k)
                if ( rmu <= s(kx) ) cycle
                iss = rmu - s(kx)
                zk = z(kx)
                if ( kx == mu ) zk = zk - 1
                coeff = sqrt(real(zk,kind(one)))*M(mu,kx)
                MDFC(r) = MDFC(r) + coeff * MDFC(iss)
            end do
    
            MDFC(r) = MDFC(r)/sqrt(real(z(mu),kind(one)))

        end do
        !-----------------!
        ! Fine calcolo FC |
        !-----------------!

        ! Free local variables
        deallocate (z)

        return
    end subroutine fcint

    !==============================================================================!

    subroutine fchtint_old()

        integer :: mu, i, j, k, kx, iss, r, zk, rpk, rmk, ixi, nvib

        real(kind=dp) :: cp, cm, coeff, FCHT_xyz, FCHT1_xyz, FCHT2_xyz, FCHT3_xyz, FCHT4
        real(kind=dp), parameter :: c2 = sqrt(one/2.0_dp)

        integer, allocatable :: z(:)

        ! If the active space is empty just exit
        if (nexvib == 0) return

        nvib = int(nmt/2) ! this must be an integer since nmt is even but I prefer casting

        allocate(FCHT(1:NFC))
        FCHT = zero

        allocate(z(1:nmt))
        !--------------------------!
        !     Inizio calcolo FCHT  !
        !--------------------------!
        z = 0

        ! This is for the <0|mu|0> term
        FCHT_xyz = zero
        do ixi = 1, 3
            FCHT_xyz = FCHT_xyz + MDFC(1)*system%tm%mu0(ixi)
            do i = 1, nvib
                FCHT_xyz = FCHT_xyz + c2*MDFC(1+s(i))*system%tm%dmuds(ixi,i)
            end do
            FCHT(1) = FCHT_xyz**2
        end do

        do r = 2, NFC

            !--------------------------------------!
            ! generate next set of quantum numbers !
            !--------------------------------------!
            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0
                    exit
                end if
            end do

            if (any(z(nvib+1:nmt) == zmax(nvib+1:nmt))) cycle !cannot compute term with z(i) -> z(i)+1

            ! Sum HT terms
XYZ:        do ixi = 1, 3

                FCHT1_xyz = system%tm%mu0(ixi)*MDFC(r)

                FCHT2_xyz = zero
                do i = nvib+1, nmt
                    if (z(i) == zmax(i)) cycle ! should be redundant
                    cp = sqrt(real(z(i)+1,kind(one))/2.0_dp)
                    FCHT2_xyz = FCHT2_xyz + system%tm%dmuds(ixi,i-nvib)*cp*MDFC(r+s(i))
                end do

                do i = nvib+1, nmt
                    if (z(i) == 0) cycle
                    cm = sqrt(real(z(i),kind(one))/2.0_dp)
                    FCHT2_xyz = FCHT2_xyz + system%tm%dmuds(ixi,i-nvib)*cm*MDFC(r-s(i))
               end do

                FCHT(r) = FCHT(r) + (FCHT1_xyz+FCHT2_xyz)**2
               ! print *, 'components ', ixi, FCHT1_xyz, FCHT2_xyz, FCHT(r)
            end do XYZ
        end do
        !--------------------!
        ! Fine calcolo FCHT  |
        !--------------------!

        deallocate(z)

        return
    end subroutine fchtint_old

    !==============================================================================!

    subroutine fchtint()

        integer :: mu, i, j, k, kx, iss, r, zk, rpk, rmk, ixi, nvib

        real(kind=dp) :: cm, coeff, FCHT_xyz, FCHT1_xyz, FCHT2_xyz, FCHT3_xyz, FCHT4
        real(kind=dp), parameter :: cp = sqrt(one/2.0_dp)

        integer, allocatable :: z(:)

        ! If the active space is empty just exit
        if (nexvib == 0) return

        nvib = int(nmt/2) ! this must be an integer since nmt is even but I prefer casting

        allocate(FCHT(1:NFC))
        FCHT = zero

        allocate(z(1:nmt))
        !--------------------------!
        !     Inizio calcolo FCHT  !
        !--------------------------!
        z = 0

        ! This is for the <0|mu|0> term
        do ixi = 1, 3
            ! Add FC term
            FCHT0(ixi) = sqrt(2.0_dp)*system%tm%mu0(ixi)
            ! Add pure HT terms
            do i = 1, nvib
                FCHT0(ixi) = FCHT0(ixi) + system%tm%dmudQ(ixi,i)*Q(i)
            end do
            FCHT1_xyz = FCHT0(ixi)*cp*MDFC(1)
            ! Here we only store the intensities, not the ampitudes
            ! We would have three sets of amplitudes to store (x,y, and z).
            !print *, ixi, FCHT_xyz
            FCHT(1) = FCHT(1) + FCHT1_xyz**2
        end do

        do r = 2, NFC

            !--------------------------------------!
            ! generate next set of quantum numbers !
            !--------------------------------------!
            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0
                    exit
                end if
            end do

            ! Sum HT terms
XYZ:        do ixi = 1, 3

                FCHT1_xyz = FCHT0(ixi)*MDFC(r)

                FCHT2_xyz = zero
                do i = 1, nvib
                    if (z(i) == 0) cycle
                    coeff = sqrt(real(z(i),kind(one)))
                    FCHT2_xyz = FCHT2_xyz + system%tm%dmudQ(ixi,i)*coeff*MDFC(r-s(i))
                end do

                FCHT3_xyz = zero
                do i = 1, nvib
                    FCHT4 = zero
                    do k = 1, nexvib
                        j = ivib(k)
                        if ( r <= s(j) ) cycle
                        iss = r - s(j)
                        zk = z(j)
                        if (j == i) zk = zk - 1
                        if (zk < 0) cycle
                        coeff = sqrt(real(zk,kind(one)))*M(i,j)
                        !print *, 'iss coeff ', iss, coeff, zk
                        FCHT4 = FCHT4 + coeff*MDFC(iss)
                    end do
                    FCHT3_xyz = FCHT3_xyz + system%tm%dmudQ(ixi,i)*FCHT4
                end do

                FCHT_xyz = (FCHT1_xyz + FCHT2_xyz+FCHT3_xyz)*cp
                ! Here we only store the intensities, not the ampitudes
                ! We would have three sets of amplitudes x,y, and z.
                FCHT(r) = FCHT(r) + FCHT_xyz**2
            end do XYZ
        end do
        !--------------------!
        ! Fine calcolo FCHT  |
        !--------------------!

        deallocate(z)

        return
    end subroutine fchtint

    !==============================================================================!

    subroutine fcint_cl (fccl)

        type(fccl_t), intent(inout) :: fccl

        ! Local variables
        real(kind=dp) :: coeff
        integer :: mu   ! Index of the excited vibration.
        integer :: i, j, k, iss, r, kx, zk, rmu, nvib
        real(kind=dp), allocatable :: freq(:)
        integer, save :: kfc = 1

        ! z is used to store the 2N quantum numbers of each FC integrals
        integer, allocatable :: z(:)

        nvib = int(nmt/2) ! this must be an integer since nmt is even but I prefer casting
        nexvib = fccl%n
        NFC  = product(fccl%zmax+1)

        ! Redefine ivib
        if (allocated(ivib)) deallocate(ivib)
        allocate(ivib(1:fccl%n))
        ivib = fccl%ivib

        ! Redefine zmax for all dofs
        zmax = 0
        zmax(fccl%ivib) = fccl%zmax

        ! redefine s
        s(1) = 1
        do i = 2, nmt
            s(i) = product(zmax(1:i-1)+1)
        end do

        ! Set the value of the <0|0> FC
        fccl%FC(1) = MDFC(1)
        !print *, fccl%FC(1)

        allocate(freq(1:nexvib))
        freq = system%state(1)%molecule%normodes%vibration(ivib-nvib)%freq
        !freq = system%state(1)%molecule%normodes%vibration(ivib)%freq

        !--------------------------!
        !     Inizio calcolo FC    !
        !--------------------------!
        allocate (z(1:nmt))
        z = 0
        do r = 2, NFC

            !--------------------------------------!
            ! generate next set of quantum numbers !
            !--------------------------------------!
            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0 ! I could insert here a boolean to check if an integral has to be writen or not.
                    exit
                end if
            end do

            rmu = r - s(mu)

            fccl%FC(r) = Q(mu) * fccl%FC(rmu)

            do k = 1, nexvib
                kx = ivib(k)
                if ( rmu <= s(kx) ) cycle
                iss = rmu - s(kx)
                zk = z(kx)
                if ( kx == mu ) zk = zk - 1
                coeff = sqrt(real(zk,kind(one)))*M(mu,kx)
                fccl%FC(r) = fccl%FC(r) + coeff * fccl%FC(iss)
            end do

            fccl%FC(r) = fccl%FC(r)/sqrt(real(z(mu),kind(one)))

            ! if any of the z is zero we are falling in a previous classs
            ! which we have already written to file.
            if (all(z(ivib) > 0)) then
                kfc = kfc + 1
                sumfc = sumfc + fccl%FC(r)**2
                !if (fccl%FC(r)**2 > 1.0d-14) then
                write(777,*) dot_product(freq,z(ivib)), fccl%FC(r)!, sumfc
                !end if
            end if
        end do
        !-----------------!
        ! Fine calcolo FC |
        !-----------------!

        !write(*,*) 'SUMFC ', sumfc
        ! Free local variables
        deallocate (z)

        return
    end subroutine fcint_cl

    !==============================================================================!

    subroutine fcint_cl_boltz (fccl,fe,qvib,tv,fg)

        type(fccl_t), intent(inout) :: fccl
        real(kind=dp), intent(in) :: qvib, fe(:)
        real(kind=dp), intent(in), optional ::fg(:),tv(:)  

        ! Local variables
        real(kind=dp) :: coeff, pv
        integer :: mu   ! Index of the excited vibration.
        integer :: i, j, k, iss, r, kx, zk, rmu, nvib
        real(kind=dp), allocatable :: freq(:)

        ! z is used to store the 2N quantum numbers of each FC integrals
        integer, allocatable :: z(:)

        nvib = int(nmt/2) ! this must be an integer since nmt is even but I prefer casting
        nexvib = fccl%n
        NFC  = product(fccl%zmax+1)

        ! Redefine ivib
        if (allocated(ivib)) deallocate(ivib)
        allocate(ivib(1:fccl%n))
        ivib = fccl%ivib

        ! Redefine zmax for all dofs
        zmax = 0
        zmax(fccl%ivib) = fccl%zmax

        ! redefine s
        s(1) = 1
        do i = 2, nmt
            s(i) = product(zmax(1:i-1)+1)
        end do

        ! Set the value of the <0|0> FC
        fccl%FC(1) = MDFC(1)
        !print *, fccl%FC(1)

        !--------------------------!
        !     Inizio calcolo FC    !
        !--------------------------!
        allocate (z(1:nmt))
        z = 0
        do r = 2, NFC

            !--------------------------------------!
            ! generate next set of quantum numbers !
            !--------------------------------------!
            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0 ! I could insert here a boolean to check if an integral has to be writen or not.
                    exit
                end if
            end do

            rmu = r - s(mu)

            fccl%FC(r) = Q(mu) * fccl%FC(rmu)

            do k = 1, nexvib
                kx = ivib(k)
                if ( rmu <= s(kx) ) cycle
                iss = rmu - s(kx)
                zk = z(kx)
                if ( kx == mu ) zk = zk - 1
                coeff = sqrt(real(zk,kind(one)))*M(mu,kx)
                fccl%FC(r) = fccl%FC(r) + coeff * fccl%FC(iss)
            end do

            fccl%FC(r) = fccl%FC(r)/sqrt(real(z(mu),kind(one)))

            ! if any of the z is zero we are falling in a previous class
            ! which we have already written to file.
            if (all(z(ivib) > 0)) then
                pv = 1/qvib
                if (present(fg)) pv = product(exp(-tv*z(fccl%gvib)))*pv
                sumfc = sumfc + pv*(fccl%FC(r)**2)
                !if (fccl%FC(r)**2 > 1.0d-2) then
                  !print *, 'IVIB', ivib
                  !print *, 'ZVIB ', z(ivib)
                  if (present(fg)) then
                    !print *, fg(1), z(1)
                    !pv = product(exp(-tv*z(fccl%gvib)))/qvib
                    write(777,*) dot_product((/fg, fe/),z(ivib)), pv*(fccl%FC(r)**2)
                  else
                    write(777,*) dot_product(fe,z(ivib)), pv*(fccl%FC(r)**2)
                  end if
                !end if
            end if
        end do
        !-----------------!
        ! Fine calcolo FC |
        !-----------------!

        !write(*,*) 'SUMFC ', sumfc
        ! Free local variables
        deallocate (z)

        return
    end subroutine fcint_cl_boltz

    !==============================================================================!

    subroutine fcint_cl_twostate (fccl,fe,fg)

        type(fccl_t), intent(inout) :: fccl
        real(kind=dp), intent(in) :: fe(:)
        real(kind=dp), intent(in), optional ::fg(:)

        ! Local variables
        real(kind=dp) :: coeff, pv
        integer :: mu   ! Index of the excited vibration.
        integer :: i, j, k, iss, r, kx, zk, rmu, nvib
        real(kind=dp), allocatable :: freq(:)

        ! z is used to store the 2N quantum numbers of each FC integrals
        integer, allocatable :: z(:)

        nvib = int(nmt/2) ! this must be an integer since nmt is even but I prefer casting
        nexvib = fccl%n
        NFC  = product(fccl%zmax+1)

        ! Redefine ivib
        if (allocated(ivib)) deallocate(ivib)
        allocate(ivib(1:fccl%n))
        ivib = fccl%ivib

        ! Redefine zmax for all dofs
        zmax = 0
        zmax(fccl%ivib) = fccl%zmax

        ! redefine s
        s(1) = 1
        do i = 2, nmt
            s(i) = product(zmax(1:i-1)+1)
        end do

        ! Set the value of the <0|0> FC
        fccl%FC(1) = MDFC(1)
        !print *, fccl%FC(1)

        !--------------------------!
        !     Inizio calcolo FC    !
        !--------------------------!
        allocate (z(1:nmt))
        z = 0
        do r = 2, NFC

            !--------------------------------------!
            ! generate next set of quantum numbers !
            !--------------------------------------!
            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0 ! I could insert here a boolean to check if an integral has to be writen or not.
                    exit
                end if
            end do

            rmu = r - s(mu)

            fccl%FC(r) = Q(mu) * fccl%FC(rmu)

            do k = 1, nexvib
                kx = ivib(k)
                if ( rmu <= s(kx) ) cycle
                iss = rmu - s(kx)
                zk = z(kx)
                if ( kx == mu ) zk = zk - 1
                coeff = sqrt(real(zk,kind(one)))*M(mu,kx)
                fccl%FC(r) = fccl%FC(r) + coeff * fccl%FC(iss)
            end do

            fccl%FC(r) = fccl%FC(r)/sqrt(real(z(mu),kind(one)))

            ! if any of the z is zero we are falling in a previous class
            ! which we have already written to file.
            if (all(z(ivib) > 0)) then
                sumfc = sumfc + (fccl%FC(r)**2)
                !if (fccl%FC(r)**2 > 1.0d-2) then
                  !print *, 'IVIB', ivib
                  !print *, 'ZVIB ', z(ivib)
                  if (present(fg)) then
                    !print *, fg(1), z(1)
                    !pv = product(exp(-tv*z(fccl%gvib)))/qvib
                    write(777,*) dot_product((/fg, fe/),z(ivib)), fccl%FC(r)
                  else
                    write(777,*) dot_product(fe,z(ivib)), fccl%FC(r)
                  end if
                !end if
            end if
        end do
        !-----------------!
        ! Fine calcolo FC |
        !-----------------!

        !write(*,*) 'SUMFC ', sumfc
        ! Free local variables
        deallocate (z)

        return
    end subroutine fcint_cl_twostate

    !==============================================================================!

    subroutine fchtint_cl(fccl)

        type(fccl_t), intent(inout) :: fccl

        real(kind=dp), allocatable :: freq(:)
        integer :: mu, i, j, k, kx, iss, r, zk, rpk, rmk, ixi, nvib

        real(kind=dp) :: cm, coeff, FCHT_xyz, FCHT1_xyz, FCHT2_xyz, FCHT3_xyz, FCHT4
        real(kind=dp), parameter :: cp = sqrt(one/2.0_dp)

        integer, allocatable :: z(:)

        nvib = int(nmt/2) ! this must be an integer since nmt is even but I prefer casting
        nexvib = fccl%n
        NFC  = product(fccl%zmax+1)

       ! Redefine ivib
        if (allocated(ivib)) deallocate(ivib)
        allocate(ivib(1:fccl%n))
        ivib = fccl%ivib
        allocate(freq(1:nexvib))
        freq = system%state(1)%molecule%normodes%vibration(ivib-nvib)%freq

        ! Redefine zmax for all dofs
        zmax = 0
        zmax(fccl%ivib) = fccl%zmax

        ! redefine s
        s(1) = 1
        do i = 2, nmt
            s(i) = product(zmax(1:i-1)+1)
        end do

        allocate(fccl%FCHT(1:NFC))
        fccl%FCHT = zero
        fccl%FCHT(1) = MDFCHT00

        allocate(z(1:nmt))
        !--------------------------!
        !     Inizio calcolo FCHT  !
        !--------------------------!
        z = 0

        do r = 2, NFC

            !--------------------------------------!
            ! generate next set of quantum numbers !
            !--------------------------------------!
            do i = 1, nexvib
                mu = ivib(i)
                if (z(mu) < zmax(mu)) then
                    z(mu) = z(mu) + 1
                    if (i > 1) z(1:ivib(i-1)) = 0
                    exit
                end if
            end do

            ! Sum HT terms
XYZ:        do ixi = 1, 3

                FCHT1_xyz = FCHT0(ixi)*fccl%FC(r)

                FCHT2_xyz = zero
                do i = 1, nvib
                    if (z(i) == 0) cycle
                    coeff = sqrt(real(z(i),kind(one)))
                    FCHT2_xyz = FCHT2_xyz + system%tm%dmuds(ixi,i)*coeff*fccl%FC(r-s(i))
                end do

                FCHT3_xyz = zero
                do i = 1, nvib
                    FCHT4 = zero
                    do k = 1, nexvib
                        kx = ivib(k)
                        if ( r <= s(kx) ) cycle
                        iss = r - s(kx)
                        zk = z(kx)
                        if (kx == i) zk = zk - 1
                        if (zk < 0) cycle
                        coeff = sqrt(real(zk,kind(one)))*M(i,kx)
                        !print *, 'iss coeff ', iss, coeff, zk
                        FCHT4 = FCHT4 + coeff*fccl%FC(iss)
                    end do
                    FCHT3_xyz = FCHT3_xyz + system%tm%dmuds(ixi,i)*FCHT4
                end do

                FCHT_xyz = (FCHT1_xyz + FCHT2_xyz+FCHT3_xyz)*cp
                ! Here we only store the intensities, not the ampitudes
                ! We would have three sets of amplitudes x,y, and z.
                fccl%FCHT(r) = fccl%FCHT(r) + FCHT_xyz**2
            end do XYZ


            ! if any of the z is zero we are falling in a previous classs
            ! which we have already written to file.
            if (all(z(ivib) > 0)) then
                sumfcht = sumfcht + fccl%FCHT(r)
                write(778,'(g12.6,2x,E16.6,<nexvib>(i4,2x),<nexvib>(i4,2x))') dot_product(freq,z(ivib)), fccl%FCHT(r), z(ivib), ivib-nvib
            end if

        end do
        !--------------------!
        ! Fine calcolo FCHT  |
        !--------------------!

        deallocate(z)

        return
    end subroutine fchtint_cl

    !==============================================================================!

    subroutine fcht00()

        integer :: ixi, i, nvib
        real(kind=dp) :: FCHT1_xyz
        real(kind=dp), parameter :: cp = sqrt(one/2.0_dp)

        nvib = int(nmt/2)

        ! This is for the <0|mu|0> term
        MDFCHT00 = zero
        do ixi = 1, 3
            ! Add FC term
            FCHT0(ixi) = sqrt(2.0_dp)*system%tm%mu0(ixi)
            ! Add pure HT terms
            do i = 1, nvib
                FCHT0(ixi) = FCHT0(ixi) + system%tm%dmuds(ixi,i)*Q(i)
            end do
            FCHT1_xyz = FCHT0(ixi)*cp*MDFC(1)
            ! Here we only store the intensities, not the ampitudes
            ! We would have three sets of amplitudes to store (x,y, and z).
            MDFCHT00 = MDFCHT00 + FCHT1_xyz**2
        end do

        ! Initialize sumfcht (this is used to check the convergence of the calculation
        ! of FCHT integrals.
        sumfcht = MDFCHT00


        return
    end subroutine fcht00

    !==============================================================================!

    ! This subroutine implements the algorithm proposed by Santoro (JCP 2008)

    subroutine fcclasses(nclasses)

        ! We divide the N vibrations into classes and then compute the FCs
        ! for each class. We can go up to 4 class
        integer, intent(in) :: NCLASSES
        integer :: i, nvib, nfccl, nsmx(1:7)
        integer, allocatable :: cvib(:), umax(:)
        real(kind=dp), parameter :: epsfc1 = 0.0010_dp
        ! The main variable which stores the FC integrals for each class
        type(fccl_t), allocatable :: fccl(:)

        ! Define overall zmax
        nvib = int(nmt/2)

        ! WRITE THE <0|0> FC
        write(777,*) zero, MDFC(1)

        !nsmx = (/30, 10, 5, 4, 2, 2, 1 /)
        nsmx = (/30, 10, 5, 2, 1, 1, 1 /)
        allocate(umax(1:nmt))
        !umax(1:nvib) = 60 
        umax(nvib+1:nmt) = 10 ! this is just to test the algorithm.
                              ! Excite all vibrations of the excited state with N quanta.
        allocate(cvib(1:NCLASSES))
        allocate(fccl(1:NCLASSES))
        do i = 1, NCLASSES
            umax(nvib+1:nmt) = nsmx(i)
            write(fout,*) 'Computing class :', i,' integrals.'
            fccl(i)%n = i ! class order
            allocate(fccl(i)%ivib(1:i))
            allocate(fccl(i)%zmax(1:i))
            cvib = 0
            do
                call combi(cvib(1:i),nvib,i)
                if (cvib(1) == 0) exit
                !print *, cvib(1:i)
                fccl(i)%ivib = cvib(1:i)+nvib !correct?
                !fccl(i)%ivib = cvib(1:i)       !wrong?
                ! set active space zmax
                fccl(i)%zmax = umax(fccl(i)%ivib)
                nfccl = product(fccl(i)%zmax+1)
                !print *, 'nfccl ', nfccl
                allocate(fccl(i)%fc(1:nfccl))
                ! compute FC
                call fcint_cl(fccl(i))
                ! deallocate FC since we are changing class
                deallocate(fccl(i)%fc)
                ! Here the main problem is to avoid the computation of the FC integrals
                ! which have already been computed in the previous class. For this
                ! reason my vectorial approach is not suitable and a binary tree description
                ! of the FC integrals is necessary.

                ! FC integrals are written in the fcint_cl() subroutine.
            end do
            write(fout,*) 'FCSUM: ', sumfc
        end do

        stop
        return
    end subroutine fcclasses

    !==============================================================================!

    subroutine fcclasses_boltz(nclasses,Temp)

        ! We divide the N vibrations into classes and then compute the FCs
        ! for each class. We can go up to 4 class
        integer, intent(in) :: NCLASSES
        real(kind=dp), intent(in) :: Temp
        integer :: i, j, i1, k, nvib, ngvib, nfccl, NGCLASSES
        integer, allocatable :: cvib(:), gvib(:), umax(:)
        real(kind=dp), parameter :: epsfc1 = 0.0010_dp, eps_boltz = 1.0E-3
        ! The main variable which stores the FC integrals for each class
        !type(fccl_t), allocatable :: fccl(:)
        type(fccl_t) :: fccl
        real(kind=dp), allocatable :: freq(:), fg(:), fe(:), tv(:)
        real(kind=dp) :: qvib, cv, thetav, p0, ptot, pv

        ! Define overall zmax
        nvib = int(nmt/2)

        allocate(freq(1:nvib))
        freq = system%state(2)%molecule%normodes%vibration(1:nvib)%freq

        ! this is just to test the algorithm.
        allocate(umax(1:nmt))
        umax(1:nvib) = 3 
        ! Excite all vibrations of the excited state with N quanta.
        umax(nvib+1:nmt) = 3 
       

        !allocate(nmax(0:nvib)) ! we need nmax(0) = 1 to have one cycle for the vibrational ground state.
        !nmax(0) = 1
        !nmax(1:nvib) = 10

        ! Compute vibrational partition function of the electronic ground state
        qvib = 1.0_dp
        do k = 1, nvib 
            thetav = hplanck*freq(k)*clgt/bk
            cv = exp(-thetav/Temp)
            qvib = qvib/(1.0_dp-cv)
        end do

        ! ground state population
        pv =  1.0_dp/qvib
        ! Modify the initial value of FC sum.
        sumfc = sumfc*pv
        ! WRITE THE <0|0> intensity
        write(777,*) zero, pv*MDFC(1)**2

        NGCLASSES=4
        allocate(cvib(1:NCLASSES),gvib(1:NGCLASSES))
        allocate(fg(1:NGCLASSES),fe(1:NCLASSES),tv(1:NGCLASSES))
        !allocate(fccl(1:NCLASSES))

        ngvib=nvib
NG:     do k = 0, NGCLASSES
          gvib = 0
CG:       do
            if (k > 0) then
              if (allocated(fccl%gvib)) deallocate(fccl%gvib)
              allocate(fccl%gvib(1:k))
              call combi(gvib(1:k),ngvib,k)
              !print *, 'gvib ', gvib(1:k)
              if (gvib(1) == 0) exit
              tv(1:k) = hplanck*freq(gvib(1:k))*clgt/(bk*Temp)
              !print *, 'B ', k, tv(1:k), exp(-tv(1:k)/Temp)/qvib, 1.0/qvib
              fccl%gvib = gvib(1:k)
              !cv = exp(-tv/Temp) 
              !pv = cv/qvib
            end if
            ! Avoid including high energy frequency 
            if (sum(tv(1:k)) > 1.8) cycle
            !ptot = ptot + pv
            !print *, 'Boltz: ', k
            ! skip if the population is less than a threshold
            !if (pv < eps_boltz) cycle
            !print *, 'nclasses ', nclasses
NC:         do i = 1, NCLASSES
              !allocate(fe(1:i)) ! At the moment fg is a scalar
              i1 = i + k !
              fccl%n = i1 ! excited state class order + ground state excited modes
              if (k > 0) then 
                fg = 0
                fg(1:k) = -system%state(2)%molecule%normodes%vibration(gvib(1:k))%freq
                !print *, 'gvib ', gvib(1:k), fg(1:k)
              end if
            ! Loop over state with a significant Boltzman population
              write(fout,*) 'Computing class :', i,' integrals.'
              print *, 'class ', i
              allocate(fccl%ivib(1:i1))
              allocate(fccl%zmax(1:i1))
              cvib = 0
              !print *, 'cvib', cvib
CE:           do
                  call combi(cvib(1:i),nvib,i)
                  fe(1:i) = system%state(1)%molecule%normodes%vibration(cvib)%freq
                  !print *, 'cvib', cvib
                  if (cvib(1) == 0) exit
                  if (k == 0) then
                    fccl%ivib = cvib(1:i)+nvib !correct?
                    fccl%zmax = umax(cvib(1:i)+nvib)
                  else
                    fccl%ivib = (/gvib(1:k), cvib(1:i)+nvib/) !correct?
                    !print *, 'ivib ', fccl(i)%ivib
                    !print *, 'umax ', umax(nvib+1:nmt)
                    fccl%zmax = (/umax(gvib(1:k)), umax(cvib(1:i)+nvib)/)
                  end if
                  nfccl = product(fccl%zmax+1)
                  allocate(fccl%fc(1:nfccl))
                  ! compute FC
                  if (k == 0) then
                    call fcint_cl_boltz(fccl,fe(1:i),qvib)
                  else 
                    call fcint_cl_boltz(fccl,fe(1:i),qvib,tv(1:k),fg(1:k))
                  end if
                  ! deallocate FC since we are changing mode combination
                  deallocate(fccl%fc)
                  ! FC integrals are written to file in the fcint_cl() subroutine.
              end do CE
              deallocate(fccl%ivib,fccl%zmax)
              !deallocate(fe)
              write(fout,*) 'FCSUM: ', sumfc
            end do NC
            if (k == 0) exit
          end do CG
        end do NG

        stop
        return
    end subroutine fcclasses_boltz

    !==============================================================================!

    subroutine fcclasses_twostate(nclasses1,nclasses2)

        ! We divide the N vibrations into classes and then compute the FCs
        ! for each class. We can go up to 4 class
        integer, intent(in) :: NCLASSES1,NCLASSES2
        integer :: i, j, i1, k, nvib, ngvib, nfccl
        integer, allocatable :: cvib(:), gvib(:), umax(:)
        real(kind=dp), parameter :: epsfc1 = 0.0010_dp, eps_boltz = 1.0E-3
        ! The main variable which stores the FC integrals for each class
        !type(fccl_t), allocatable :: fccl(:)
        type(fccl_t) :: fccl
        real(kind=dp), allocatable :: freq(:), fg(:), fe(:)
        real(kind=dp) :: p0, ptot

        ! Define overall zmax
        nvib = int(nmt/2)

        allocate(freq(1:nvib))
        freq = system%state(2)%molecule%normodes%vibration(1:nvib)%freq

        ! this is just to test the algorithm.
        allocate(umax(1:nmt))
        umax(1:nvib) = 3
        ! Excite all vibrations of the excited state with N quanta.
        umax(nvib+1:nmt) = 3


        ! WRITE THE <0|0> intensity
        write(777,*) zero, MDFC(1)

        allocate(cvib(1:NCLASSES2),gvib(1:NCLASSES1))
        allocate(fg(1:NCLASSES1),fe(1:NCLASSES2))

        ngvib=nvib
NG:     do k = 0, NCLASSES1
          gvib = 0
CG:       do
            if (k > 0) then
              if (allocated(fccl%gvib)) deallocate(fccl%gvib)
              allocate(fccl%gvib(1:k))
              call combi(gvib(1:k),ngvib,k)
              !print *, 'gvib ', gvib(1:k)
              if (gvib(1) == 0) exit
              fccl%gvib = gvib(1:k)
            end if
NC:         do i = 1, NCLASSES2
              i1 = i + k !
              fccl%n = i1 ! excited state class order + ground state excited modes
              if (k > 0) then
                fg = 0
                fg(1:k) = -system%state(2)%molecule%normodes%vibration(gvib(1:k))%freq
                !print *, 'gvib ', gvib(1:k), fg(1:k)
              end if
            ! Loop over state with a significant Boltzman population
              write(fout,*) 'Computing class :', i,' integrals.'
              print *, 'class ', i
              allocate(fccl%ivib(1:i1))
              allocate(fccl%zmax(1:i1))
              cvib = 0
              !print *, 'cvib', cvib
CE:           do
                  call combi(cvib(1:i),nvib,i)
                  fe(1:i) = system%state(1)%molecule%normodes%vibration(cvib)%freq
                  !print *, 'cvib', cvib
                  if (cvib(1) == 0) exit
                  if (k == 0) then
                    fccl%ivib = cvib(1:i)+nvib !correct?
                    fccl%zmax = umax(cvib(1:i)+nvib)
                  else
                    fccl%ivib = (/gvib(1:k), cvib(1:i)+nvib/) !correct?
                    !print *, 'ivib ', fccl(i)%ivib
                    !print *, 'umax ', umax(nvib+1:nmt)
                    fccl%zmax = (/umax(gvib(1:k)), umax(cvib(1:i)+nvib)/)
                  end if
                  nfccl = product(fccl%zmax+1)
                  allocate(fccl%fc(1:nfccl))
                  ! compute FC
                  if (k == 0) then
                    call fcint_cl_twostate(fccl,fe(1:i))
                  else
                    call fcint_cl_twostate(fccl,fe(1:i),fg(1:k))
                  end if
                  ! deallocate FC since we are changing mode combination
                  deallocate(fccl%fc)
                  ! FC integrals are written to file in the fcint_cl() subroutine.
              end do CE
              deallocate(fccl%ivib,fccl%zmax)
              !deallocate(fe)
              write(fout,*) 'FCSUM: ', sumfc
            end do NC
            if (k == 0) exit
          end do CG
        end do NG

        stop
        return
    end subroutine fcclasses_twostate

    !==============================================================================!

    ! This subroutine implements the algorithm proposed by Santoro (JCP 2008).
    ! The main problem of the entire algorithm is that it is not possible to avoid the computation
    ! of the FC integrals which have already been computed in the previous class.
    ! My vectorial approach is not suitable and a binary tree description
    ! of the FC integrals might help (though the final outcoume would probably be inefficient).

    subroutine fchtclasses(nclasses)

        ! We divide the N vibrations into classes and then compute the FCs
        ! for each class. We can go up to 4 class
        integer, intent(in) :: NCLASSES
        integer :: i, nvib, nfccl
        integer, allocatable :: cvib(:), umax(:)
        real(kind=dp), parameter :: epsfc1 = 0.0010_dp
        ! The main variable which stores the FC integrals for each class
        type(fccl_t), allocatable :: fccl(:)

        ! Define overall zmax
        nvib = int(nmt/2)

        ! Write the <0|\mu|0> intensity
        write(778,*) zero, MDFCHT00

        allocate(umax(1:nmt))
        umax = 0
        umax(nvib+1:nmt) = 4 ! this is just to test the algorithm.

        allocate(cvib(1:NCLASSES))
        allocate(fccl(1:NCLASSES))
        do i = 1, NCLASSES
            write(fout,*) 'Computing class :', i,' integrals.'
            fccl(i)%n = i ! class order
            allocate(fccl(i)%ivib(1:i))
            allocate(fccl(i)%zmax(1:i))
            cvib = 0
            do
                call combi(cvib(1:i),nvib,i)
                if (cvib(1) == 0) exit
                !print *, cvib(1:i)
                fccl(i)%ivib = cvib(1:i)+nvib
                ! set active space zmax
                fccl(i)%zmax = umax(fccl(i)%ivib)
                nfccl = product(fccl(i)%zmax+1)
                !print *, 'nfccl ', nfccl
                allocate(fccl(i)%fc(1:nfccl))
                ! compute FC
                call fcint_cl(fccl(i))
                ! compute FCHT
                call fchtint_cl(fccl(i))
                ! deallocate FC since we are changing class
                deallocate(fccl(i)%fc)
                deallocate(fccl(i)%fcht)
           end do
           write(fout,*) 'FCHTSUM: ', sumfcht
        end do

        stop
        return
    end subroutine fchtclasses

    !==============================================================================!

    subroutine fcint_free

        integer alloc_err, ierr

        write(fout,'(//,2x,a)') 'Freeing memory allocated for FC Integrals.'

        if (allocated(zmax)) deallocate(zmax)
        if (allocated(s)) deallocate(s)
        if (allocated(ivib)) deallocate(ivib)

        if (allocated(MDFC)) then
            deallocate(MDFC,stat=alloc_err)
        else
            write(fout,'(//,2x,a)') 'Not using a large array to store FC integrals.'
        end if
        if (alloc_err /= 0) ierr = error(0,"Cannot deallocate array used for FC integrals.")
 
        return
    end subroutine fcint_free

    !==============================================================================!

    SUBROUTINE combi(ia,n,j)

        integer :: ia(:)
        integer i, i1, j, n

        if(n .le. 0 .or. j .le. 0) return

        if(j .gt. n) then
            write(fout,103) j,n
        else if(ia(1) .eq. 0) then
            do i = 1,j
                ia(i)=i
            end do
            ia(j+1)=0
        else
            do i1 = 1,n
                i=i1
                if(ia(i+1) .ne. ia(i)+1) exit
                ia(i)=i
            end do
            ia(i)=ia(i)+1
            if(ia(j) .eq. n+1) ia(1)=0
        endif

        return

103     format('j = ',i5,' > ',i5,' = n is not permitted')
    END SUBROUTINE combi

    !==============================================================================!

    SUBROUTINE Berkowitz() 

    integer :: nvib
    integer :: i, mz(1:10), s(1:10)
    real(kind=dp) :: fci

    ! This is just a trial version of the routine
    ! We read the sequence of integrals to be computed
    do i = 1, 10
      mz = (/1, 1, 1, 1, i, 0, 2, 0, 0, 2/)
      ! Sort vector mz(i) in descending order ( Es: [0 2 1 3 4] -> [4 3 2 1 0]) and 
      ! store its pivot sequence s(i). Using the reordered vector ensures that the 
      ! number of FC to be computed is the minimum possible
      print *, 'initial vector :', mz
      call vsort(mz,s)
      print *, 'ordered vector :', mz
      print *, 'pivot sequence :', s 
      fci = FC_BK(mz,s)
      print *, 'FC :',  fci
    end do

    END SUBROUTINE Berkowitz

    !==============================================================================!

    RECURSIVE FUNCTION FC_BK(mz,s) RESULT(fci)

    ! Algoritmo di S. Berkowitz e F. J. Garner per il calcolo ricorsivo dei polinomi 
    ! di Hermite multidimensionali. (Math. of Comput. Vol 24, pag 537, 1970)

    integer, intent(in) :: mz(:)
    integer, intent(in) :: s(:)
    integer :: v(1:size(mz)), mzp(1:size(mz)), mzjp(1:size(mz))
    integer :: i, j, k, p, NVP, f, n
    real(kind=dp) :: fci, fc1, fc2

    if (all(mz == 0)) then 
      fci = MDFC(1)
      return
    end if

    !========= Here we must find which recurrence relation to apply ================!
    ! Define class V 

    print *, 'call fcbk'

    NVP = 1
    n = count(mz > 0)
PV: do p = 1, N

      ! Define new class V(p) with NVP elements
      ! NVP = m(s(1))*m(s(2))*... : numero di elementi nella classe Vp

      NVP = NVP*mz(s(p))
      ! Initialize the first vector of the class
      v = 0 
CL:   do i = 1, NVP
        ! next vector in class Vp
        v = next_sequence(v,mz,s,p)
        ! find pivot p and apply the p-th recursion
        f = sum(mz(s(p+1:n)))+(mz(s(p))-v(s(p)))-sum(mz(s(1:p-1))-v(s(1:p-1)))
        print *, 'berk index: ', f
        if (f >= 0) then
          ! use s(p)-th recurrence relation
          mzp(s(p)) =  mz(s(p)) - 1
          fc1 = Q(s(p))*FC_BK(mzp,s) 
          fc2 = 0.0_dp
          do j = 1, N 
            mzjp = mz
            mzjp(s(p)) =  mzjp(s(p)) - 1
            mzjp(s(j)) =  mzjp(s(j)) - 1
            fc2 = fc2 + M(s(p),s(j))*mz(s(j))*FC_BK(mzjp,s) 
          end do
          fci = fc1 + fc2
        end if
      end do CL

    end do PV

    return 
    END FUNCTION FC_BK

    !==============================================================================!

    FUNCTION next_sequence(v,mz,s,p) RESULT(nv)

    integer, intent(in) :: v(:), mz(:), s(:), p
    integer :: nv(1:size(v))
    integer :: k, mu

        print *, 'input v: ', v
        print *, 'input s: ', s
        print *, 'input mz: ', mz
        print *, 'input p: ', p

        nv = v
        do k = 1, p
            if (v(s(k)) < mz(s(k))) then
                nv(s(k)) = v(s(k)) + 1
                if (k > 1) nv(1:s((k-1))) = 0
                exit
            end if
        end do

        print *, 'next v: ', nv
    END FUNCTION next_sequence

end module FCintegrals
