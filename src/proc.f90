module procsys

    use parameters
    use proc_type
    use system_type
    use sysop
    use errors
    use util, only : isempty

    implicit none

    public :: process_system

    contains

	!==============================================================================

    subroutine process_system(system, proc)
    
    type(system_t), intent(inout) :: system
    type(proc_t), intent(in) :: proc

    real(kind=dp) norm, x(1:3), tmp_coord(1:3)
    real(kind=dp) conv, scf
    !type(vibration_t), pointer :: vibr
    type(vibration_t), allocatable :: tmp_vib(:)
    type(vibration_t) :: vibr
    type(atom_t), pointer :: atom
    type(structure_t) :: struct

    integer i, j, k, l, jj, kj, is, ist, im, ierr, ios, kk, kk1, in, io, nvb, iv
    integer ilock, alloc_err
    integer, allocatable :: nori(:), nv_ord(:)
    character(len=80) or 
    integer ia, ja, nato, iwarn
    integer, allocatable :: n_ord(:)
    real(kind=dp), allocatable :: dr(:)
    real(kind=dp), parameter :: LONG = 10000.0
    real(kind=dp) :: av(1:3), bv(1:3), orth
    integer, external :: idamin
    logical, allocatable :: reo(:)

    !------------------------------------------------|
    ! Controlla numero di atomi e numero di molecole |
    !------------------------------------------------|
    if ( system%state(1)%molecule%structure%numat /=  &
               system%state(2)%molecule%structure%numat) ierr = error(67)

    !-----------------------------------------------------------+
    ! If necessary trasform to mass-weighted normal coordinates |
    ! and renormalize them. This has to be done for Gaussian and
    ! Turbomole                                                 |
    !-----------------------------------------------------------+
    do is = 1, system%nstate
        if (.not.system%state(is)%molecule%normodes%massw) then
          ! mass-weight the j-th normal coordinate
          do j = 1, system%state(is)%molecule%nvib
            do k = 1, system%state(is)%molecule%structure%numat
              system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3) = &
              system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3) * &
              sqrt(system%state(is)%molecule%normodes%vibration(j)%atom(k)%amass)
            end do
            ! Renormalize the j-th normal coordinate.
            norm = zero
            do k = 1, system%state(is)%molecule%structure%numat
              x(1:3) = system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3)
              norm = norm + dot_product(x,x)
            end do
			norm = sqrt(norm)
            do k = 1, system%state(is)%molecule%structure%numat
              system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3) = &
              system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3)/norm
            end do
          end do
        end if
    end do

        ! check ortho
!    do is = 1, system%nstate
!        do i = 1, system%state(is)%molecule%nvib
!          do j = 1, system%state(is)%molecule%nvib
!            orth = 0.0_dp
!            do k = 1, system%state(is)%molecule%structure%numat
!	      av(1:3) = system%state(is)%molecule%normodes%vibration(i)%atom(k)%d(1:3)
!	      bv(1:3) = system%state(is)%molecule%normodes%vibration(j)%atom(k)%d(1:3)
!	      orth = orth + dot_product(av,bv)
!            end do
!            print *, 'ORTHO: ', i, j,  orth
!          end do
!        end do
!    end do

    !--------------------------------------------------------------------+
    ! Converte la geometria da Angstrom a Bohr (se specificato in input).| 
    !--------------------------------------------------------------------+
    do is = 1, system%nstate
        conv = ONE
        if (trim(system%state(is)%molecule%structure%units) == "bohr") conv = BOHR
        if (conv == ONE) exit
        write(fout,'(a,i2,a,i2,a)') '  Geometry of molecule ', i, ' state ', is ,' converted to Angstrom.'
        do j = 1, system%state(is)%molecule%structure%numat
          system%state(is)%molecule%structure%atom(j)%coord = &
          system%state(is)%molecule%structure%atom(j)%coord * conv
        end do

        conv = ONE
        if (trim(system%state(is)%molecule%normodes%units) == "bohr") conv = BOHR
        if (conv == ONE) exit

        do iv= 1, system%state(is)%molecule%nvib
          do j = 1, system%state(is)%molecule%structure%numat
            system%state(is)%molecule%normodes%vibration(iv)%atom(j)%d = & 
            system%state(is)%molecule%normodes%vibration(iv)%atom(j)%d * conv
          end do
        end do
    end do

    !---------------------------------------------------------------------+
    ! Scala le frequenze di uno stato col fattore scfreq (dato in input). | 
    !---------------------------------------------------------------------+
    do is = 1, system%nstate
        scf = system%state(is)%molecule%normodes%scfreq
        if (scf == ONE) exit
        write(fout,*) '  Scaling factor for frequencies of system%state ', is, ' molecule  ', j, ' : ', scf 
        do j = 1, system%state(is)%molecule%nvib
          system%state(is)%molecule%normodes%vibration(j)%freq = & 
          system%state(is)%molecule%normodes%vibration(j)%freq * scf
        end do
    end do

    !------------------------------------------------------------------------------------------+
    !  Riorienta gli assi del sistema di riferimento della molecola secondo quanto specificato |
    !  dal vettore <frame>.                                                                    |
    !------------------------------------------------------------------------------------------+
    allocate(nori(1:3))
    do i = 1, size(proc%frame(:))
        is = get_state_from_id(system,proc%frame(i)%state)
        system%state(is)%molecule%reorient = .true.
        !jj = 0
        !call build_data_array(proc%frame(i)%axis,system%state(is)%molecule%orient,jj)
        system%state(is)%molecule%orient = proc%frame(i)%axis
        nori = abs(system%state(is)%molecule%orient)
        !if (ios /= 0) ierr = error(0,"Cannot read new frame orientation. Check <frame> tag.")
        !------------------------+
        ! Riorienta la struttura |
        !------------------------+
        do kk = 1, system%state(is)%molecule%structure%numat 
          tmp_coord = system%state(is)%molecule%structure%atom(kk)%coord(:)
          system%state(is)%molecule%structure%atom(kk)%coord(:) = tmp_coord(nori(:))
          forall (io=1:3) system%state(is)%molecule%structure%atom(kk)%coord(io) = & 
                          sign(abs(one),real(system%state(is)%molecule%orient(io),kind(one)))* &
                          system%state(is)%molecule%structure%atom(kk)%coord(io)
        end do
        !-------------------------+ 
        ! Riorienta le vibrazioni |
        !-------------------------+
		!if (.true.) then
        do j = 1, system%state(is)%molecule%nvib
          do kk = 1, system%state(is)%molecule%structure%numat
             tmp_coord = system%state(is)%molecule%normodes%vibration(j)%atom(kk)%d(:)
             do kj = 1, 3
                system%state(is)%molecule%normodes%vibration(j)%atom(kk)%d(kj) = tmp_coord(nori(kj))
             end do
          end do
          do kk = 1, system%state(is)%molecule%structure%numat
            forall (io=1:3) system%state(is)%molecule%normodes%vibration(j)%atom(kk)%d(io) = &
                            sign(abs(one),real(system%state(is)%molecule%orient(io),kind(one)))* & 
                            system%state(is)%molecule%normodes%vibration(j)%atom(kk)%d(io) 
          end do
        end do
		!end if
    end do
    deallocate(nori)

    ! ----------------------------------------------------------------+ 
    ! N.B.: Le strutture ed i modi normali si possono riordinare SOLO |
    ! dopo aver riorientato gli assi del sistema.                     |
    ! ----------------------------------------------------------------+

    !-----------------------------------------------------------------+ 
    ! Riordina le geometrie ed i modi normali di una coppia di stati. |
    !-----------------------------------------------------------------+

    !-----------------------------------------------------------------+
    ! Definisce un array di interi che e' il nuovo ordine degli atomi |
    ! usando la sequenza specificata in <order> (letta in input).     |
    !-----------------------------------------------------------------+
    if (allocated(proc%reorder)) then
        allocate(reo(1:system%nstate))
        reo = .false.
RD:     do i = 1, size(proc%reorder(:))

            if (proc%reorder(i)%data == "atoms") then
                if (.not.isempty(proc%reorder(i)%ord)) then
                !print *, proc%reorder(i)%ord
                !print *, 'reordering not empty '

		            is = get_state_from_id(system,proc%reorder(i)%state)
                    nato = system%state(is)%molecule%structure%numat
                    allocate(system%state(is)%molecule%structure%order(1:nato),stat=alloc_err)
                    !ierr = error(0,"NOT IMPLEMENTED ANYMORE. TELL THE AUTHOR TO DO IT")
                    read(proc%reorder(i)%ord,*)(system%state(is)%molecule%structure%order(k),k=1,nato)
                    !call build_data_array(proc%reorder(i)%ord,system%state(is)%molecule%structure%order(1:nato))
                    allocate(n_ord(1:nato),stat=alloc_err)
                    !-----------------------+
                    ! Riordina le strutture |
                    !-----------------------+
                    n_ord = system%state(is)%molecule%structure%order(1:nato)
					allocate(struct%atom(1:nato))
                    struct%atom = system%state(is)%molecule%structure%atom
                    system%state(is)%molecule%structure%atom(1:nato) = struct%atom(n_ord(1:nato))
					deallocate(struct%atom)
                    !---------------------------------------+
                    ! Riordina gli atomi nei i modi normali |
                    !---------------------------------------+
					allocate(vibr%atom(1:nato))
                    do j = 1, system%state(is)%molecule%nvib
						vibr = system%state(is)%molecule%normodes%vibration(j)
                        system%state(is)%molecule%normodes%vibration(j)%atom(1:nato) = vibr%atom(n_ord(1:nato))
                    end do
					deallocate(vibr%atom)
                    deallocate(n_ord)
				else
				!print *, 'reorder ilock'
					ilock = 1 
					if (.not.isempty(proc%reorder(i)%reference)) then 
						ilock = get_state_from_id(system,proc%reorder(i)%reference)
					end if
					! Lock the selected state and reorder atoms in the other states accordingly
                    do ist = 1, system%nstate
						if (ist == ilock) cycle
                        is = ist
						nato = system%state(is)%molecule%structure%numat
                        nvb = system%state(is)%molecule%nvib
                        allocate(system%state(is)%molecule%structure%order(1:nato),stat=alloc_err)
                        if (alloc_err /= 0) ierr = error(0,"Cannot allocate order() array in sysop.")
                        allocate(n_ord(1:nato),stat=alloc_err)
                        system%state(is)%molecule%structure%order(1:nato) = get_order(system,ilock,is)
                        !-----------------------+
                        ! Riordina le strutture |
                        !-----------------------+
                        n_ord = system%state(is)%molecule%structure%order(1:nato)
				 		allocate(struct%atom(1:nato))
				 		struct%atom = system%state(is)%molecule%structure%atom
				 		system%state(is)%molecule%structure%atom(1:nato) = struct%atom(n_ord(1:nato))
				 		deallocate(struct%atom)
                        !---------------------------------------+
                        ! Riordina gli atomi nei i modi normali |
                        !---------------------------------------+
				 		allocate(vibr%atom(1:nato))
                        do j = 1, system%state(is)%molecule%nvib
							vibr = system%state(is)%molecule%normodes%vibration(j)
							system%state(is)%molecule%normodes%vibration(j)%atom(1:nato) = vibr%atom(n_ord(1:nato))
                        end do
				 		deallocate(vibr%atom)
                        deallocate(n_ord)
                    end do
				end if
            else if (proc%reorder(i)%data == "vibrations") then 
                is = get_state_from_id(system,proc%reorder(i)%state)
                nvb = system%state(is)%molecule%nvib
                ! L'ordine delle vibrazioni deve essere specificato in input.  Non esiste un  modo
                ! per definirlo in maniera automatica. (o meglio...non ne ho implementato nessuno).
                allocate(system%state(is)%molecule%normodes%order(1:nvb))
                allocate(n_ord(1:nvb))
                !ierr = error(0,"NOT IMPLEMENTED. TELL THE AUTHOR TO DO IT!!")
                !call build_data_array(proc%reorder(i)%ord,system%state(is)%molecule%normodes%order(1:nvb))
                print *, proc%reorder(i)%ord
                read(proc%reorder(i)%ord,*)(system%state(is)%molecule%normodes%order(k),k=1,nvb)
                n_ord = system%state(is)%molecule%normodes%order(1:nvb)
				allocate(tmp_vib(1:nvb))
				tmp_vib = system%state(is)%molecule%normodes%vibration(1:nvb)
                system%state(is)%molecule%normodes%vibration(1:nvb) = tmp_vib(n_ord(1:nvb))
				deallocate(tmp_vib)
                deallocate(n_ord)
                print *, 'OK..'
            end if
        end do RD
    end if

    !------------------------------------------------------------------+
    ! Aassegna le ID alle vibrazioni. E' fondamentale per definire le  |
    ! vibrazioni eccitate.                                             |
    !------------------------------------------------------------------+
    do is = 1, size(system%state)
        do j = 1, system%state(is)%molecule%nvib
          system%state(is)%molecule%normodes%vibration(j)%id = j
        end do
    end do

	!-----------------------------------------------------------------------------------+
	! Calcola le connectivity matrix delle molecole nei diversi stati elettronici       |
	! e controlla che siano le stesse. In caso contraro il reordering non ha funzionato |
	! e bisogna farlo a mano.														    |
	!-----------------------------------------------------------------------------------+

	do is = 1, system%nstate
		call system%state(is)%molecule%setacm()
	end do
	
	do i = 1, system%nstate
		do j = i+1, system%nstate
		if (any(system%state(i)%molecule%acm /= system%state(j)%molecule%acm))  &
		ierr = warning(0,"Detected different connectivities of the two electronic states. Please check your data!")
		end do
	end do

	if (allocated(proc%subset)) then
	    do i = 1, system%nstate
	        do j = 1, system%nstate
	            print *, 'SUBSET ID ',proc%subset(i)%state
	            print *, 'STATE ID', system%state(i)%id
	            if (proc%subset(j)%state == system%state(i)%id)  then
	                write(fout,*) 'Subset for STATE ', system%state(i)%id
	                call system%state(i)%molecule%getSubset(proc%subset(j)%incvib)
	                exit
	            end if
	        end do
	    end do
	end if

    include 'formats'

    return
    end subroutine process_system
    
end module procsys
