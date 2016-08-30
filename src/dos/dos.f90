module dos
    
    use parameters
    use errors
    use matfun, only : vsort
    use xmvar, only : system
    use dos_type
    use sysop
      
    implicit none

    private

    public:: kdb, brsw
            
CONTAINS
     
    !===========================================================================!
      
    SUBROUTINE kdb(dosjob)
    
        !-------------------------------------------------------------------------------!
        ! This subroutine utilized the Kemper - Van Dijk - Buck backtracking algorithm  !
        ! for state search.                                                             !
        ! See: Chem. Phys. Lett. vol 53, pag 121 (1978)                                 !
        !-------------------------------------------------------------------------------!
    
        type(dos_t), intent(in) :: dosjob
    
        integer, parameter :: nqmax = 3 ! maximum quantum numbers allowed per mode
        integer, allocatable :: a(:), amax(:), ids(:), iv(:)
        real(kind=dp), allocatable :: freq(:)
        real(kind=dp) emin, emax, energy, eplus
        integer i, j, k, l, m, n, ng, istate, imol, nmd, nv, nv0, ic, nc
    
        ! number of dof:
        istate = get_state_from_id(system, dosjob%state)
        nmd = system%state(istate)%molecule%nvib

        allocate(freq(1:nmd))
        allocate(amax(1:nmd),a(1:nmd),ids(1:nmd))
        amax = nqmax    ! maximum quantum number allowed per mode
    
        ! frequencies in descending order
        freq = system%state(istate)%molecule%normodes%vibration(:)%freq
        
        call vsort(freq, ids, order='D') ! sort frequencies in descending order es: 3000, 2000, 200, 12, 2
        ! Set energy range for state search.
        emin = dosjob%emin
        emax = dosjob%emax

        ! open file to write the states of the search result.
        open(unit=fnst,file=dosjob%file//'.dos',status="unknown")
        ! open file to write the indices of vibrations excited in file FNST.
        open(unit=fivb,file=dosjob%file//'.ivb',status="unknown")
        ! Controllo che le frequenze siano state ordinate correttamente
        !print *, freq
        !print *, emin, emax
    
        energy = zero
        a = 0
        i = 1
        ng = 0
        do
            a(i) = a(i) + 1    ! eccita il modo i-esimo
            energy = energy + freq(i)
            if ((energy >= emin) .and. (energy <= emax) .and. (a(i) <= amax(i))) then
                ! Il ciclo for che segue serve solo per la stampa degli indici, si può rimuovere ai fini dell'algoritmo.
                !-----
                nc = count(a > 0)
                if (allocated(iv)) deallocate(iv)
                allocate(iv(1:nc))
                ic = 0    					! trova gli indici dei modi eccitati
                do l = 1, nmd    			! non so se questo ciclo for verrà mantenuto durante lo sviluppo
                   if (a(l) > 0) then		! dipende sa quanto rallenta l'algoritmo di ricerca.
                   	ic = ic + 1			! Devo fare un profiling...
                   	iv(ic) = l
                   end if
                end do
                write(fivb,'(<nc>(i3,1x))')(iv(l),l=1,nc) ! scrive gli indici dei modi eccitati
                write(fnst,'(<nc>i3,2x,f10.2)')(a(iv(l)),l=1,nc), energy ! scrive lo stato trovato e la sua energia
                ng = ng + 1    ! aggiorna il numero di stati eccitati
            else if ((energy > emax) .or.(a(i) > amax(i))) then     ! se l'energia è troppo alta o se si è superato
                a(i) = a(i) - 1            						! il numero max di quanti sul modo i-esimo
                energy = energy - freq(i)
                if (i < nmd) then ! passa al modo successivo
                    i = i + 1
                    cycle
                end if
                if (i == nmd) then    ! se si è raggiunto l'ultimo modo
                    j = -1
                    do k = nmd-1, 1, -1    ! indice del penultimo modo eccitato
                        if ((a(k) > 0)) then
                            j = k
                            exit
                        end if
                    end do
                    if (j < 0) exit    ! se solo il modo nmd è eccitato esce
                    a(j) = a(j) - 1    ! scala un quanto sul penultimo modo eccitato
                    energy = energy - freq(j)    ! scala anche l'energia di un quanto sul modo j
                    do l=j+1,nmd        ! ricalcola l'energia
                        energy = energy - a(l)*freq(l)
                    end do
                    a(j+1:nmd) = 0    ! mette a zero tutti i numeri quantici dal penultimo in poi
                    i = j + 1    ! i = penultimo eccitato + 1
                end if
            end if
        end do

        write(fout,99021)dosjob%state,emin, emax
        write(fout,99023) ng + 1
        write(fout,99024)dosjob%file

        include 'formats_dos'
    
        return
    END SUBROUTINE kdb

    !============================================================================!

    SUBROUTINE brsw(dosjob)

        ! Beyer and Swinehart algorithm applied to the density of states

        integer, allocatable :: t(:) , ifr(:)
        type(dos_t), intent(in) :: dosjob
        real(kind=dp), allocatable :: freq(:), dens(:), ssum(:)
        real(kind=dp) emin, emax, egrain, e
        integer i, j, nmd, imax, nv, nv0, imol, istate, nvic, ic, iarg
	
        ! number of dof:
        istate = get_state_from_id(system, dosjob%state)
        nmd = system%state(istate)%molecule%nvib

        allocate(freq(1:nmd))
    
        ! frequencies in descending order
        freq = system%state(istate)%molecule%normodes%vibration(:)%freq
        
        ! Set energy range for state search.
        emin = dosjob%emin
        emax = dosjob%emax
        egrain = dosjob%egrain
        imax = 1 + int(emax/egrain)
    	
        ! Inizializza l'array ir che contiene le frequenze arrotondate all'intero più vicino
        nvic = count(freq <= emax) ! number of included vibrations. Vibration with freq(i) > emax.
						           ! will not contribute to the number of states states with energy below emax.
        allocate(ifr(1:nvic))
        ic = 0
        do i = 1, nmd
            if (freq(i) <= emax) then
                ic = ic + 1
                ifr(ic) = nint(freq(i)/egrain)
            end if
        end do

        ! initialize t
        allocate(t(1:imax))
        t = zero
        ! this actually counts the number of multiply partitions of an integer n
        ! The algorithm is the one originally proposed by Beyer and Swinehart.
        call countmrp(ifr,t,imax)
	
        allocate(ssum(1:imax),dens(1:imax))
	
        ssum = zero; dens = zero
	
        ssum(1) = t(1)                    ! energy corresponds to top of energy grain
        dens(1) = t(1)/egrain
        do i = 2 , imax
            ssum(i) = t(i) + ssum(i-1)      ! sum of states at top of energy grain
            dens(i) = t(i)/egrain
        end do

        if ( ssum(1) .lt. 1.0 ) then      ! force semiclassical approximation
            ssum(1) = 1.0_dp
            dens(1) = 1.0_dp/egrain
        endif

        open(unit=fdos,file=dosjob%file//'.dos',status="unknown")
	
        do i = 1 , imax
            e = (i-1)*egrain
            write (fdos,*)  e, dens(i) , ssum(i)
        end do
	
        return
    END SUBROUTINE brsw

    !==========================================================================!

    subroutine countmrp(C, P, N)
        ! 	COUNTING MULTIPLY RESTRICTED PARTITIONS OF AN INTEGER
        !   ALGORITHM 448, COLLECTED ALGORITHMS FROM ACM.
            !   THIS WORK PUBLISHED IN COMMUNICATIONS OF THE ACM
            !   VOL. 16, NO. 6, June, 1973, P.379.
        integer, intent(in) :: N
        integer, intent(in) :: C(:)
        integer, intent(out) :: P(1:N)
        INTEGER I, J, K, JP1, M, MMJ
        !  COUNT COMPUTES THE NUMBER OF PARTITIONS OF AN INTEGER
        !  RESTRICTED TO C FOR INTEGERS IN THE RANGE 1 TO N.
        !  INPUT:  K -- A POSITIVE INTEGER.
        !          C -- AN ARRAY OF K POSITIVE INTEGERS.
        !          N -- AN INTEGER LARGER THAN THE MAXIMUM VALUE IN C.
        !  OUTPUT: P -- AN ARRAY OF N INTEGERS, WHERE P(M) IS THE
        !               NUMBER OF PARTITIONS OF M RESTRICTED TO C.
        !  INITIALIZE P
        P = 0
        !  EACH PASS THROUGH THE OUTER LOOP BELOW TRANSFORMS P FROM
        !  PARTITIONS RESTRICTED TO C(1), ..., C(I-1) TO
        !  PARTITIONS RESTRICTED TO C(1), ..., C(I).
        K = size(C)
        do I = 1, K
            J = C(I)
            JP1 = J + 1
            P(J) = P(J) + 1
            do M = JP1, N
                MMJ = M - J
                P(M) = P(M) + P(MMJ)
            end do
        end do

        return
    end subroutine countmrp

end module dos
