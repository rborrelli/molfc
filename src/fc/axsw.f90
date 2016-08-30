module axsw

  use parameters
  use iomatrix
  use xmvar, only : job
  use system_type
  use sysop
  use errors
  use job_type
  use transf_type

  implicit none

  real(kind=dp), allocatable :: axrot0(:,:)

  real(kind=dp), allocatable :: crd1(:,:), crd0(:,:), crd(:,:), wm(:)

  real(kind=dp) c(1:3,1:3), ti0(1:3,1:3), ti1(1:3,1:3), &
                ctc(1:3,1:3), work(1:8), e(1:3), lbd(1:3,1:8), rm(1:3,1:3), &
                cm1(1:3,1:3), t0tmp(1:3,1:3), cdm(1:3), scat(1:8)
  real(kind=dp) cond1, cond2, cond3, scat0, det, t011, t012, t013, &
                t0ij, tmd0, tmd1

  real(kind=dp), parameter :: toler = 1.5D0 ! Questo valore deve essere controllato: e' solo un trial.
  logical passed(1:8)

  integer i, j, k, isol, ipass, info, nato, natoi, natof, nm


  private

  public :: axis_switching

  contains

!====================================================================================

        subroutine axis_switching (molecule, jobtrn)

        type(molecule_t), intent(inout) :: molecule(1:2)
        type(transf_t), intent(in) :: jobtrn

        write(fout,799)

        call axis_switching_geom(molecule, jobtrn)

        call axis_switching_vibr(molecule(2), jobtrn)
        
        if (allocated(axrot0)) deallocate(axrot0)

799     format(/,2x,'+---------------------------+',/, &
                 2x,'| Axis Switching Procedure  |',/, &
                 2x,'+---------------------------+',/)

        return
        end subroutine axis_switching

!====================================================================================

        subroutine axis_switching_geom (molecule, jobtrn)

        !---------------------------------------------------------------------------+ 
        ! Subroutine utilizzata per la determinazione della matrice                 |
        ! di axis switchig. E' stata aggiornata ed ora funziona anche per molecole  |
        ! planari.                                                                  |
        !                                                                           |
        ! N.B.:                                                                     |
        ! Affinche' la procedura di axis switching funzioni correttamente e'        |
        ! necessario che gli atomi delle due geometrie siano correttamente          |
        ! ordinati. Nel programma viene quindi chiamata prima la subroutine         |
        ! reorder_state quindi axis_switching.                                      |
        !---------------------------------------------------------------------------+

        type(molecule_t), intent(inout) :: molecule(1:2)
        type(transf_t), intent(in) :: jobtrn

        real(kind=dp), allocatable :: axsw0(:,:,:)
        integer ierr, iwarn, fsol(1:1), imx(1:1)
        real(kind=dp) molw, sq(1:8), vdiff(1:3)
		integer(kind=dp), pointer :: mloc
		
        allocate(axrot0(1:3,1:3))
        allocate(axsw0(1:8,1:3,1:3))

        axrot0 = zero
        axsw0 = zero

! If both are linear do not perform axis switching

          if (isLinear(molecule(1)) .and. isLinear(molecule(2))) return

! If there is a linear -> non-linear transition exits.

          if (isLinear(molecule(1)) .and. (.not.isLinear(molecule(2))) &
              .or. &
              isLinear(molecule(2)) .and. (.not.isLinear(molecule(1)))) &
          ierr = error(0,"Cannot handle Linear -> Non Linear transitions")
    
          nato = molecule(1)%structure%numat

          allocate (wm(1:nato), crd(1:3,1:nato), crd0(1:3,1:nato), crd1(1:3,1:nato))

          crd0 = get_geometry(molecule(1))
          crd  = get_geometry(molecule(2))
          wm = get_atmass(molecule(1))
         
         !--------------------------------------------------------+
         ! Calcola i tensori di inerzia delle due strutture       |
         ! prima dell'axis switching e con le geometrie di input. |
         !--------------------------------------------------------+
ITEN:    if (jobtrn%printlevel >= 5) then
            write(fout,'(2x,a,/)')'Inertia tensors before axis switching'
            tmd0 = zero; tmd1 = zero
            do i = 1,3
              do j = 1, nato
                tmd0 = tmd0 + wm(j)*(crd0(i,j)**2)
                tmd1 = tmd1 + wm(j)*(crd(i,j)**2)
              end do
            end do
           
            ti0 = zero; ti1 = zero
            do i=1,3
              do j=1,i
                do k=1,nato
                  ti0(i,j) = ti0(i,j) - wm(k)*crd0(i,k)*crd0(j,k)
                  ti1(i,j) = ti1(i,j) - wm(k)*crd(i,k)*crd(j,k)
                  ti0(j,i) = ti0(i,j)   
                  ti1(j,i) = ti1(i,j)   
                end do
              end do
            ti0(i,i) = tmd0 + ti0(i,i)
            ti1(i,i) = tmd1 + ti1(i,i)
            end do
           
            ti0 = ti0 / nav
            ti1 = ti1 / nav

            write(fout,820) molecule(1)%id
            call matout (ti0,3,3,3)
            write(fout,820) molecule(2)%id
            call matout (ti1,3,3,3)

          end if ITEN

          !------------------------------------------------------------+
          ! Calcola il centro di massa e sposta le geometrie nel CM.   |
          ! Serve come controllo perche' le geometrie di input devono  |
          ! gia' essere in standard orientation.                       |
          !------------------------------------------------------------+

          molw = sum(wm(1:nato)) 
          cdm = zero
          do i = 1, 3
            cdm(i) = sum(wm(1:nato)*crd0(i,1:nato))/molw
            crd0(i,1:nato) = crd0(i,1:nato) - cdm(i)
          end do

          if (jobtrn%debug) then
          write(fout,'(/,2x,''Center of mass of state 1'')')
          write(fout,'(2x,3(f7.4,2x))') cdm
          if (any(abs(cdm) > 1.0d-5)) iwarn = warning(0,' Molecule 1 not in center of mass.')
          end if
         
          cdm = zero
          do i = 1, 3
            cdm(i) = sum(wm(1:nato)*crd(i,1:nato))/molw
            crd(i,1:nato) = crd(i,1:nato) - cdm(i)
          end do

          if (jobtrn%debug) then
          write(fout,'(/,2x,''Center of mass of state 2'')')
          write(fout,'(2x,3(f7.4,2x))') cdm
          if (any(abs(cdm) > 1.0d-5)) iwarn = warning(0,' Molecule 2 not in center of mass.')
          end if

          !----------------------------------------------------------+
          ! Valuta le condizioni di Eckhart alle geometrie iniziali. |
          !----------------------------------------------------------+
          cond1 = zero; cond2 = zero; cond3 = zero
         
          do i = 1, nato
            cond1 = cond1+wm(i)*(crd(2,i)*crd0(1,i)-crd(1,i)*crd0(2,i))
            cond2 = cond2+wm(i)*(crd(3,i)*crd0(1,i)-crd(1,i)*crd0(3,i))
            cond3 = cond3+wm(i)*(crd(3,i)*crd0(2,i)-crd(2,i)*crd0(3,i))
          end do
         
          scat0 = sqrt(cond1**2+cond2**2+cond3**2)
         
          if (jobtrn%printlevel >= 3) write(fout,823) molecule(1)%id, scat0
         
          !---------------------------------------------------------+
          ! Definisce tre matrici per il calcolo, c, ctc e lbd      |
          ! Il sistema che risolvo ha 8 soluzioni di cui solo 4 con |
          ! determinante 1 e solo 1  e' quella esatta.              |
          ! Trovo quella esatta paragonando le geometrie ruotate.   |
          !---------------------------------------------------------+
         
          lbd = reshape((/ one,  one,  one,  &
                          -one, -one, -one,  & 
                          -one,  one,  one,  &
                          -one, -one,  one,  &
                           one, -one,  one,  &
                           one, -one, -one,  & 
                           one,  one, -one,  &
                          -one,  one, -one /),(/3,8/))

          c = zero
          do i = 1, 3
            do j = 1, 3
              do k = 1, nato
                c(i,j) = c(i,j) + wm(k)*crd(i,k)*crd0(j,k) 
              end do
            end do
            !-----------------------------------------------------------------+
            ! Se la molecola e' planare (...o quasi).                         |
            ! Il valore 1.0d-4 e' "empirico", determinato cioe' a partire da  |
            ! prove su alcune molecole planari....                            |
            !-----------------------------------------------------------------+
            if (abs(c(i,i)) <= 1.0d-4) c(i,i) = one
          end do
         
          if (jobtrn%debug) then
            write(fout,'(/)')
            write(fout,831)
            call matout (c,3,3,3)
          end if

          !-----------------------+
          ! Calcola CTC = C(tr)*C |
          !-----------------------+
          ctc = zero
          do i = 1, 3
            do j = 1, 3
              do k = 1, 3
                ctc(i,j) = ctc(i,j) + c(k,i)*c(k,j)
              end do
            end do
          end do

          !--------------------------------------------------+
          ! Diagonalizza CTC: in CTC ci sono gli autovettori |
          ! e in e(k) gli autovalori corrispondenti.         |
          ! Serve per calcolare (C(tr)*C)^(1/2)              |
          !--------------------------------------------------+

          call dsyev ('v', 'U', 3, ctc, 3, e, work, 8, info)

          !-------------------------------------------------+
          ! Calcola c^-1 tramite il determinante e i minori |
          !-------------------------------------------------+

          rm = zero
         
          rm(1,1) = c(2,2)*c(3,3)-c(2,3)*c(3,2)
          rm(1,2) = c(2,1)*c(3,3)-c(2,3)*c(3,1)
          rm(1,3) = c(2,1)*c(3,2)-c(2,2)*c(3,1)
         
          rm(2,1) = c(1,2)*c(3,3)-c(1,3)*c(3,2)
          rm(2,2) = c(1,1)*c(3,3)-c(1,3)*c(3,1)
          rm(2,3) = c(1,1)*c(3,2)-c(1,2)*c(3,1)
         
          rm(3,1) = c(1,2)*c(2,3)-c(1,3)*c(2,2)
          rm(3,2) = c(1,1)*c(2,3)-c(1,3)*c(2,1)
          rm(3,3) = c(1,1)*c(2,2)-c(1,2)*c(2,1)
          
          det = c(1,1)*rm(1,1)-c(1,2)*rm(1,2)+c(1,3)*rm(1,3)
         
          cm1 = zero
          do i = 1, 3
            do j = 1, 3
              cm1(i,j) = ((-one)**(i+j))*rm(j,i)/det
            end do
          end do
         
		! --------------------- Inizia la soluzione del sistema ----------------

          !-------------------------------------------+
          ! La matrice di rotazione e' definita come: |
          !        T0= (C(tr)*C)^(1/2)*C^(-1)         |
          !-------------------------------------------+

          ipass = 0 
          passed = .false.

		  sq = -1
SOL:      do isol = 1, 8

            !-------------------------+
            ! Calcolo (C(tr)*C)^(1/2) |
            !-------------------------+

            !---------------------------------------------------------------+
            ! Il segno degli autovettori e' arbitrario per questo ci sono 8 |
            ! possibili soluzioni al problema.
            !---------------------------------------------------------------+
            do i = 1, 3
              do j = 1, 3
                t0ij = zero
                do k = 1, 3
                  t0ij = t0ij + ctc(i,k)*dble(lbd(k,isol))*sqrt(e(k))*ctc(j,k)
                end do
                ! Metto (C(tr)*C)^(1/2) in axsw0 |
                !axsw0(isol,i,j) = t0ij
                t0tmp(i,j) = t0ij
              end do
            end do
           
            axsw0(isol,:,:) = matmul(t0tmp,cm1)

            ! Calcola il determinante di axsw0
            t011 = axsw0(isol,2,2)*axsw0(isol,3,3)-axsw0(isol,2,3)*axsw0(isol,3,2)
            t012 = axsw0(isol,2,1)*axsw0(isol,3,3)-axsw0(isol,2,3)*axsw0(isol,3,1)
            t013 = axsw0(isol,2,1)*axsw0(isol,3,2)-axsw0(isol,2,2)*axsw0(isol,3,1)
            
            det = axsw0(isol,1,1)*t011-axsw0(isol,1,2)*t012+axsw0(isol,1,3)*t013
         
            !--------------------------------------------------------------------------------
            !  Ruota il sistema di riferimento
            !  Considera sia le rotazioni proprie, cioe' quelle con determinante 1
            !  che quelle con determinante -1 (rotazioni improprie). Queste ultime possono
            !  essere necessarie quando gli assi sono orientati in modo diverso.
            !  Il codice  e' modificato in modo da tener conto solo delle rotazioni
            !  proprie con un if (det > 0) then...
            !--------------------------------------------------------------------------------
                 
            !------------------------------------------------------------------------------
            !  Questa parte del programma presenta una serie di bug per cui e' necessario 
            !  controllare l'output per bene!!!  
            !------------------------------------------------------------------------------
         
ROT0:       if (det > 0) then
          
              crd1 = zero
              !-----------------------------------------------------+
              ! Ruota la geometria crd con axsw0 e la mette in crd1 |
              !-----------------------------------------------------+
              crd1 = matmul(axsw0(isol,1:3,1:3),crd)

              !---------------------------------------------------------------------------+         
              ! Controlla se la geometria trasformata con axsw0 e' simile alla geometria  |
              ! del primo stato elettronico                                               |
              !---------------------------------------------------------------------------+
			  sq(isol) = zero
			  do i = 1, nato
				  vdiff = crd1(1:3,i)-crd0(1:3,i)
				  sq(isol) = sq(isol) + dot_product(vdiff,vdiff)
			  end do
			  sq(isol) = sqrt(sq(isol))

              !-----------------------------------------------------+
              ! Verifica se la condizione di Eckart e' soddisfatta: |
			  ! questo è solo un check perché la soluzione è esatta |
              !-----------------------------------------------------+
              cond1 = zero; cond2 = zero; cond3 = zero
            
              do i = 1, nato
                cond1 = cond1+wm(i)*(crd1(2,i)*crd0(1,i)-crd1(1,i)*crd0(2,i))
                cond2 = cond2+wm(i)*(crd1(3,i)*crd0(1,i)-crd1(1,i)*crd0(3,i))
                cond3 = cond3+wm(i)*(crd1(3,i)*crd0(2,i)-crd1(2,i)*crd0(3,i)) 
              end do
             
              scat(isol) = sqrt(cond1**2+cond2**2+cond3**2)

              !-------------------------------------+
              ! Stampa la matrice di axis-switching |
              !-------------------------------------+
              if (jobtrn%debug) then
                write(fout,825) isol
                call matout (axsw0(isol,1:3,1:3),3,3,3)
                write(fout,'(/,2x,''Matrix determinant:'',f8.4)') det
                write(fout,'(/,2x,''Zero order Eckart Conditions after rotation: '',e12.5)') scat(isol)
                write(fout,'(/,2x,''Sum of atomic distances : '',e12.5,/)') sq(isol)
              end if

            end if ROT0
           
          end do SOL

          !----------------------------------------------------------------+
          ! Dopo il ciclo sulle 8 soluzioni possibili utilizza la migliore |
          ! delle soluzioni che verifica le condizioni sopra imposte.      |
          !----------------------------------------------------------------+
		  if (jobtrn%axsw%solution == -1) then
			!fsol(1:1) = minloc(scat,mask=passed)
			fsol(1:1) = minloc(sq,(sq > 0))
			!if (.not.any(passed)) ierr = error(0,"Cannot find any axis switching matrix.")
		  else
			fsol(1:1) = jobtrn%axsw%solution
		  end if

			if (jobtrn%debug) then
				write(fout,'(2x,a,i4)') 'Axis Switching Matrix selected among the four solutions:', fsol(1)
			end if
			axrot0(1:3,1:3) = axsw0(fsol(1),1:3,1:3)

          if (jobtrn%printlevel >= 1) then
            write(fout,825) fsol(1)
            call matout (axrot0(1:3,1:3),3,3,3)
            write(fout,'(/,2x,''Matrix determinant:'',f8.4)') det
            write(fout,'(/,2x,''Zero order Eckart Conditions after rotation: '',e12.5)') scat(fsol(1))
          end if

          !-------------------------------------------------------+          
          ! Aggiorna i valori delle coordinate della II struttura |
          !-------------------------------------------------------+          
          crd1 = matmul(axrot0(1:3,1:3),crd)
          crd = crd1

          !-----------------------+
          ! Aggiorna le strutture |
          !-----------------------+
          call set_geometry(molecule(1),crd0)
          call set_geometry(molecule(2),crd)

          !--------------------------------------------------------+
          ! Calcola i tensori di inerzia delle due nuove strutture |
          !--------------------------------------------------------+
          tmd0 = zero; tmd1 = zero
          do i = 1,3
            do j = 1, nato
              tmd0 = tmd0 + wm(j)*(crd0(i,j)**2)
              tmd1 = tmd1 + wm(j)*(crd1(i,j)**2)
            end do
          end do
         
          ti0 = zero; ti1 = zero
          do i=1,3
            do j=1,i
              do k=1,nato
                ti0(i,j) = ti0(i,j) - wm(k)*crd0(i,k)*crd0(j,k)
                ti1(i,j) = ti1(i,j) - wm(k)*crd1(i,k)*crd1(j,k)
                ti0(j,i) = ti0(i,j)   
                ti1(j,i) = ti1(i,j)   
              end do
            end do
            ti0(i,i) = tmd0 + ti0(i,i)
            ti1(i,i) = tmd1 + ti1(i,i)
          end do
         
          ti0 = ti0 / nav
          ti1 = ti1 / nav

          if (jobtrn%printlevel >= 3) then
            write(fout,820) molecule(1)%id
            call matout (ti0,3,3,3)
            write(fout,820) molecule(1)%id
            call matout (ti1,3,3,3)
          end if
         
          deallocate (wm, crd, crd0, crd1)

        !if (job%control%printlevel >= 4) write(fout,305)

        !write(fout,306)

        include 'formats'

        return 
        end subroutine axis_switching_geom

!====================================================================================

        subroutine axis_switching_vibr(molecule, jobtrn)

        !---------------------------------------------------------------------+
        ! Ruota le matrici dei modi normali con la matrice di Axis Switching. |
        !---------------------------------------------------------------------+

        type(molecule_t), intent(inout) :: molecule
        type(transf_t), intent(in) :: jobtrn

        integer i, j, k, l
        real(kind=dp) disp(1:3), t_disp(1:3)

        if (jobtrn%axsw%on) then
          if (isLinear(molecule)) return
          do j = 1, molecule%nvib
            do k = 1, molecule%structure%numat
              t_disp = molecule%normodes%vibration(j)%atom(k)%d(1:3)
              disp = matmul(axrot0(:,:),t_disp)
              molecule%normodes%vibration(j)%atom(k)%d(1:3) = disp(1:3)
            end do
          end do
        end if

        return
        end subroutine axis_switching_vibr

end module axsw
