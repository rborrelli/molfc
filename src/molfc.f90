!+---------------------------------------------------------------+ 
!|                                                               |
!|  Franck-Condon Integrals Calculation and Quantum Dynamics     |
!|  Software.                                                    |
!|                                                               |
!|  References:                                                  |
!|  1) R. Borrelli, A. Peluso, J. Chem. Phys. 119, 8437 (2003).  |
!|                                                               |
!|  author
!|  Raffaele Borrelli                                            |
!|                                                               |
!+---------------------------------------------------------------+ 

program MOLFC

    use parameters
    use isotopes, only : init_isotopes, free_isotopes
    use input, only : read_input
    use fc_type
    use proc_type
    use errors
    use iomatrix
    use xmvar, only : job, system, proc
    use sysop
    use axsw
    use cartesian_transformation
    use internal_transformation
    use natural_internal_transformation
    use model_transformation
    use fc
    use fcintegrals
    use fcint_pert
    use spectrum
    use agric
    use intc
    use active
    use output
    use procsys
    use fcio
    use kubo
    use dos
    use tm_derivative
    use reorg

    implicit none
   
    character(len=40) :: finpname ! --- Input file name.
    character(len=40) :: foutname ! --- Output file name.
    character(len=50) :: buffer
    integer :: i, ierr, count_arg, stat_arg
    integer, allocatable :: zm(:)
    logical :: testonly = .false.

    integer :: is, im, ij, j, ic, id, k, l, p
    integer :: nvib
    real(kind=dp), allocatable :: freq_final_state(:), freq_initial_state(:)
    real(kind=dp) :: mdfc00
    type(molecule_t) :: molecola(1:2)
    type(activespace_t) :: acs, acs_initial, acs_final

    integer, external :: nargs
    count_arg = nargs()
    if (count_arg < 2) ierr = error(17)

    foutname = BLANK
    finpname = BLANK
    i = 1
    do while (i <= count_arg - 1)
        call getarg(i,buffer,stat_arg)
        if (stat_arg .lt. 0) ierr = error(17)
        if (buffer == "-o") then
            i = i + 1
            call getarg(i,foutname,stat_arg)
            if (stat_arg .lt. 0) ierr = error(19)
            i = i + 1
        else if (finpname == BLANK) then
            call getarg(i,finpname,stat_arg)
            i = i + 1
        end if
        if (buffer == "-t") testonly = .true.
    end do

    if (foutname == finpname) ierr = error(18)

    !-----------------------------+
    ! Nome di default dell output |
    !-----------------------------+
    if (foutname == BLANK) then
        i = index(finpname,'.',back=.true.)
        if (i /= 0) then
            foutname = finpname(1:i-1)//'.out'
        else
            foutname = finpname(1:len_trim(finpname))//'.out'
        end if
    end if

    !-----------------------------------!
    ! Apre il file di output principale !
    !-----------------------------------!
    open (unit=fout,file=foutname,status='unknown')
    ! ---------------------------------------------------------+
    ! Scrive l'intestazione del programma, il PID e lo USERID. |
    ! ---------------------------------------------------------+
    call job_headers
    !------------------------------------------------------+
    ! Richiama le masse di alcuni elementi piu' comuni     |
    !------------------------------------------------------+
    call init_isotopes
    !-----------------+
    ! Legge l'input.  |
    !-----------------+
    call read_input(finpname)
    !-----------------+
    ! Free memory...  |
    !-----------------+
    call free_isotopes
    !-----------------------------------------------------------------------+
    !  Controlla che gli stati letti in input siano tra loro compatibili,   |
    !  e applica altre trasformazioni a geometrie e modi normali:           |
    !  reorder, orient.                                                     |
    !-----------------------------------------------------------------------+
    if (.not.system%model) call process_system(system, proc)
    !-----------------------------------------------+
    ! Scrive le informazioni sul sistema, e sui job |
    !-----------------------------------------------+
    call print_system(system)

       !--------------------------------------!
       !  Lancia i job specificati in input   !
       !--------------------------------------!
    NJ: do ij = 1, size(job)
        !-------------------------------------------------------------------+ 
        ! Determina la trasformazione tra le coordinate normali del sistema.|
        !-------------------------------------------------------------------+ 
        !---------------------------------------------------------------------------------+
        ! Per ogni molecola calcolo la trasformazione tra coordinate normali
        !                              q(FINAL_STATE) = J*q(INITIAL_STATE) + d
        ! Here: FINAL_STATE = 1
        !       INITIAL_STATE = 2 
        ! Check the parameters module for the actual value of FINAL_STATE/INITIAL_STATE
        !---------------------------------------------------------------------------------+
        molecola = (/ system%state(FINAL_STATE)%molecule, system%state(INITIAL_STATE)%molecule /)
        if (job(ij)%trns%axsw%on) then
            call axis_switching(molecola,job(ij)%trns)
            ! reassign the molecule...I don't remember why...(?)
            system%state(FINAL_STATE)%molecule = molecola(1)
            system%state(INITIAL_STATE)%molecule = molecola(2)
            call print_structure(system%state(FINAL_STATE),fout)
            call print_structure(system%state(INITIAL_STATE),fout)
            call print_normal_modes(system%state(FINAL_STATE),fout)
            call print_normal_modes(system%state(INITIAL_STATE),fout)
        end if
        !+-------------------------------------------------------+
        !| Chose cartesian or internal coordinate representation |
        !| of normal modes.                                      |
        !+-------------------------------------------------------+
        if (job(ij)%trns%cartesian) then
            !------------------------------------------------------+
            ! Determine Transformation using Cartesian coordinates |
            !------------------------------------------------------+
            call cart_transf(molecola,job(ij)%trns)
            call print_cart_transf()
            call print_cart_transf_file()
        else if (job(ij)%trns%internal) then
            !----------------------------------------------------------+
            ! Define Atom connectivities and determine Wilson B-matrix |
            !----------------------------------------------------------+
            do im = 1, 2
                if (job(ij)%trns%intauto) call gric(molecola(im),job(ij)%trns)
                call set_equilibrium_intc(molecola(im))
                call bmat(molecola(im))
            end do
            ! Set internal coordinate values as given in an extergnal file
            !if (job(ij)%trns%setic) call read_intc(molecola(1),job(ij)%trns%icfile)
            !if (job(ij)%trns%setic) call read_intc(molecola(2),job(ij)%trns%icfile)
            if (job(ij)%trns%setic) call read_intc(molecola(1), molecola(2), job(ij)%trns%icfile)
            !----------------------------------------------------------------+
            ! Determine Transformation using linearized Internal coordinates |
            !----------------------------------------------------------------+
            call intc_transf(molecola,job(ij)%trns)
            call print_intc_transf()
            call print_intc_transf_file()
            call print_intc_ham_file()
        else if (job(ij)%trns%natint) then
            do im = 1, 2
                call set_equilibrium_natint(molecola(im))
                call nat_bmat(molecola(im))
            end do
            call nat_intc_transf(molecola,job(ij)%trns)
            call print_nat_intc_transf()
        else if (job(ij)%trns%model) then
            ! read transformation from file
            call model_transf(molecola,job(ij)%trns)
            call print_model_transf()
            call print_model_ham_file()
        end if
       
        ! Compute reorganizatione energies
        call reorganization_energies(molecola(1),job(ij)%trns)
 
NM:   do im = 1, job(ij)%nmeth

JOBS:   select case(job(ij)%method(im))

                case("fc")
                    ic = 1
                    FCJ:  if (allocated(job(ij)%fc)) then
                        call print_job(job(ij)%fc(ic))

                        !---------------------------------!
                        ! Compute Herzberg-Teller terms   !
                        !---------------------------------!
                        if (job(ij)%fc(ic)%fcht) call ht(job(ij)%trns)

                        !-----------------------------------------------------------+
                        ! Perform a FC calculation separately for each group of     |
                        ! vibrations defined in the <group> tag.                    |
                        !-----------------------------------------------------------+
                GROUPS: do i = 1, job(ij)%fc(ic)%ngroup
                            ! trova i modi normali che definiscono il gruppo
                            ! Il programma calcola i FC assumento la trasformazione nella forma:
                            !                   q(IS1) = J*q(IS2) + d
                            ! quindi i primi n numeri quantici sono del ket e gli altri sono del bra...
                            molecola = (/ system%state(FINAL_STATE)%molecule, system%state(INITIAL_STATE)%molecule /)

                            ! If we have used <subset> then we must redefine the included vibrations.
                            ! It is NOT POSSIBLE to combine <subset> and <include> in the same job!
                            if (allocated(proc%subset)) then
                              deallocate(job(ij)%fc(ic)%group(i)%incvib%id) 
                              allocate(job(ij)%fc(ic)%group(i)%incvib%id(1:molecola(1)%nvib)) 
                              ! redefine included vibrations if we have used the subset option.
                              !print *, 'redefine included vibrations'
                              job(ij)%fc(ic)%group(i)%nvib = molecola(1)%nvib 
                              forall (k=1:molecola(1)%nvib) job(ij)%fc(ic)%group(i)%incvib%id(k) = k 
                             ! do k = 1, molecola(1)%nvib
                             !   print *, 'job ', i, job(ij)%fc(ic)%group(i)%incvib%id(k) 
                             ! end do
                            end if
 
                            if (job(ij)%fc(ic)%kubo%on) then
                                nvib = molecola(1)%nvib
                                allocate(freq_final_state(1:nvib), freq_initial_state(1:nvib))
                                freq_final_state = molecola(1)%normodes%vibration(:)%freq  ! frequenze dello stato eccitato
                                freq_initial_state = molecola(2)%normodes%vibration(:)%freq  ! frequenze dello stato fondamentale
                                ! Devo trovare un modo per definire lo stato fondamentale ed eccitato.
                                ! al momento f1 sono le freq di quello eccitato e f2 quelle del fondamentale.
                                !call kubo_lineshape(job(ij)%trns,job(ij)%fc(ic)%kubo,f1,f2)
                                if (job(ij)%fc(ic)%fcht) then
                                    call kubo_lineshape(job(ij)%fc(ic)%group(i)%incvib,job(ij)%trns,job(ij)%fc(ic)%kubo,freq_initial_state,freq_final_state,system%tm)
                                else
                                    call kubo_lineshape(job(ij)%fc(ic)%group(i)%incvib,job(ij)%trns,job(ij)%fc(ic)%kubo,freq_initial_state,freq_final_state)
                                end if
                                ! skip any subsequent calculation
                                !cycle
                            end if
                            !-------------------------------------------------------------------------------------
                            ! IMPORTANT: Only normal modes defined in the vector INCVIB are included in the model.
                            ! Default: Include all normal modes!
                            !-------------------------------------------------------------------------------------
                            ! FC_INIT calculates the matrices M and Q used by FCINT
                            ! using only the vibrations included in the group.
                            !-------------------------------------------------------------------------------------
                            call fc_init(molecola,job(ij)%trns,job(ij)%fc(ic)%group(i)%incvib)

                            ! Get active space for the i-th molecule in the two electronic states
                            acs_initial = get_active_space(system,job(ij)%fc(ic)%group(i),INITIAL_STATE)
                            acs_final = get_active_space(system,job(ij)%fc(ic)%group(i),FINAL_STATE)

                            if (job(ij)%fc(ic)%pert) then
                                ! Perturbative FC calculation
                                call fcint_pt(acs_initial,acs_final,MDFC_00(),job(ij)%fc(ic),i)
                                exit
                            end if

                            ! initialize fc integrals defining the active space
                            ! MDFC_00() ! <- <0|0> FC integral: Watson formula is more stable than mine
                            mdfc00 = MDFC_00()
                            call fcint_init(acs_initial,acs_final,mdfc00)

                          !  if (job(ij)%fc(ic)%class) then
                          !      call fcclasses_twostate(2,2)
                          !  end if

                            ! Use the class approach
                            if (job(ij)%fc(ic)%class) then
                                !print *, 'Class algorithm'
                                if (job(ij)%fc(ic)%fcht) then
                                !print *, 'Class fcht algorithm'
                                    call fcht00()
                                    call fchtclasses(job(ij)%fc(ic)%nclasses)
                                else
                                    if (job(ij)%fc(ic)%Temp <= 10.0) then ! if temperature is less than 10 K, neglect Boltzmann population.
                                      call fcclasses(job(ij)%fc(ic)%nclasses)
                                    else
                                      call fcclasses_boltz(job(ij)%fc(ic)%nclasses,job(ij)%fc(ic)%Temp)
                                    end if
                                end if
                            else
                                ! compute Franck-Condon integrals using an explicitely defined active space
                                call fcint()
                                ! compute FCHT intensities
                                if (job(ij)%fc(ic)%fcht) then
                                    call fcht00()
                                    call fchtint()
                                end if
                            end if

                            !if (job(ij)%fc(ic)%berk) then
                            !  call Berkowitz()
                            !end if

                            if (job(ij)%fc(ic)%printfc) then
                                ! Print Franck-Condon integrals.
                                call fcint_write(acs_initial,acs_final,job(ij)%fc(ic),i)
                                if (job(ij)%fc(ic)%fcspec) then
                                    ! Calculate Franck-Condon spectrum
                                    call fcspec(acs_initial,acs_final,job(ij)%fc(ic),i)
                                end if
                            end if
                            call fcint_free()
                            call fc_end() ! deallocate arrays M and Q
                        end do GROUPS
                    end if FCJ

                case("dos")

                    if (job(ij)%dos%method == "kdb") then
                        call kdb(job(ij)%dos)
                    else if (job(ij)%dos%method == "brsw") then
                        call brsw(job(ij)%dos)
                    else
                        ierr = error(0,"Wrong method for Density of State calculation.")
                    end if

            end select JOBS
        end do NM
    end do NJ

    call job_endings

end program MOLFC

