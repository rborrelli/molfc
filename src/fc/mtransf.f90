module model_transformation

  use fc_h
  use job_type
  use util
  use fox_dom
  
  implicit none

  private
  
  real(kind=dp), private, allocatable :: KMT(:)
  
  type(molecule_t), pointer, private :: molecule(:)
  type(transf_t), pointer, private :: trns

  integer nmv

  public :: model_transf, print_model_transf, print_model_transf_file, print_model_ham_file


  contains

    !==============================================================================!

    subroutine model_transf (molecola, transf)
    
    type(molecule_t), target, intent(in) :: molecola(1:2)
    type(transf_t), target, intent(inout) :: transf

    integer i, j, ij, ierr, jj

    ! molecule e trns sono accessibili a tutto il modulo.
    molecule => molecola
    trns => transf

    ! Numero totale di vibrazioni per stato
    nmv = molecule(1)%nvib

    !----------------------------------------+
    ! Alloca le variabili del modulo FC.     !
    ! KMT viene deallocatao in print_transf  !
    !----------------------------------------+
    allocate (trns%JM(1:nmv,1:nmv), trns%KM(1:nmv), KMT(1:nmv))
    !----------------------------------------------------------------+
    ! A questo punto le geometrie sono quelle finali e le matrici
    ! dei modi normali sono state riordinate in base alle geometrie.
    !
    ! Calcola la matrice di Duschinsky GM
    !
    call nmrot()
    !
    ! Calcola il vettore di spostamenti dei modi normali GK.
    ! La geometria del secondo stato elettronico e' stata trasformata da AXSMT. 
    !
    call nmshifts() 
    !
    ! Calcola e stampa il determinante della trasformazione.    
    !
    call metric_determinant()

    include 'formats'

    return
    end subroutine model_transf

    !==============================================================================!

    subroutine nmshifts 
   
    real(kind=dp) deltax, sh
    real(kind=dp), allocatable :: KSH(:)

    integer  i, j, k, ierr
   
    !---------------------------------+
    ! Legge gli spostamenti da file   |
    !---------------------------------+     
    if (.not.isempty(trns%modk%file)) then
        open(unit=fkm,file=trns%modk%file,status="old")
        do j = 1, nmv
            read(fkm,*)trns%KM(j)
        end do
    end if

    ! Legge gli spostamenti dalla stringa interna modk%data
    if (.not.isempty(trns%modk%data)) then
        !ierr = error(0,"NOT IMPLEMENTED. PLEASE USE AN EXTERNAL FILE TO DEFINE THE SHIFT.")
       ! call build_data_array(trns%modk%data,trns%KM)
        read(trns%modk%data,*)(trns%KM(j),j=1,nmv)
    end if

    ! normal modes displacements set to zero is not given
    if (isempty(trns%modk%data) .and. isempty(trns%modk%file)) trns%KM = zero
     
    KMT = -matmul(trns%JM,trns%KM)
    
    return
    end subroutine nmshifts 
    
    !==============================================================================!

    subroutine nmrot 

    !------------------------------------------------------------------+  
    ! Calcola la matrice di trasfomazione tra i modi normali.          |
    ! Se si e' fatto l'axis switching T2 e' stata gia' ruotata con S0  |
    ! da axsmt.                                                        |
    !------------------------------------------------------------------+  
 
    integer i, j, k, ierr, ndusch, id, jd
    character(len=200) data
    real(kind=dp), allocatable :: vec(:)
    
    ! Legge la matrice di Duschinsky da file
    if (.not.isempty(trns%modr%file)) then
        open(unit=fjm,file=trns%modr%file,status="old")
        read(fjm,*) ndusch
        do i = 1, ndusch
            read(fjm,*)id, jd, trns%JM(id,jd)
        end do
    end if
    ! Legge la matrice di Duschinsky dalla stringa interna modr%data
    if (.not.isempty(trns%modr%data)) then
        if (.not.isempty(trns%modr%file)) ierr = error(0,'Duschinsky matrix already read from file. Remove data from input.')
        allocate(vec(1:nmv**2))
        !ierr = error(0,"NOT IMPLEMENTED. PLEASE USE AN EXTERNAL FILE TO DEFINE THE DUSCHINSKY MATRIX.")
        !call build_data_array(trns%modr%data,vec)
        read(trns%modr%data,*) (vec(id),id=1,nmv**2)
        do i = 1, nmv
            trns%JM(i,1:nmv) = vec(1+nmv*(i-1):i*nmv)
        end do
        deallocate(vec)
    end if

    if (isempty(trns%modr%data) .and. isempty(trns%modr%file)) then
        ! Duschinsky matrix set to identity.
        trns%JM = zero
        forall (i=1:nmv) trns%JM(i,i) = one
    end if
    
    write(fout,261)
    ! Scrive la matrice di Duschinsky in output.
    call layout(trns%JM,nmv,nmv,nmv,nmv)

    include 'formats'
    return
    end subroutine nmrot 
  
    !==============================================================================!
   
    subroutine metric_determinant 
   
    real(kind=dp) :: detgm
    detgm = det(trns%JM)
   
    write(fout,10) detgm
   
    include 'formats'

    return        
    end subroutine metric_determinant 

    !==============================================================================!        
   
    subroutine print_model_transf
    
    real(kind=dp), allocatable  :: vaus(:)
    real(kind=dp) conv, conv1
    integer i, j, l, k, icorr, icorr2, icorr3
    integer, allocatable  :: nmixed(:), mmixed(:,:), irec(:)
    integer Kcount, isp   
!    real(kind=dp), parameter :: gtol = 0.96
!    real(kind=dp), parameter :: gtolb = 0.10
   
    write(fout,140)
    
    write(fout,2401) molecule(1)%id
    write(fout,'(/)')
    write(fout,2413)
    ! New output version. Normal modes mixing is still missing.
    do k = 1, molecule(1)%nvib
        if (trns%modk%type == "dimless") then
            conv = one
            conv1 = one/sqrt(cfac*molecule(1)%normodes%vibration(k)%freq)
        else
            conv = sqrt(cfac*molecule(1)%normodes%vibration(k)%freq)
            conv1 = one
        end if
        write(fout,161) k, molecule(1)%normodes%vibration(k)%freq, &
                         molecule(2)%normodes%vibration(k)%freq, &
                         trns%KM(k)*conv, KMT(k)*conv, &
                         trns%KM(k)*conv1, KMT(k)*conv1
    end do
    
!    if (job%transf%printlevel >= 3) then
      write(fout,163)
      write(fout,2401) molecule(1)%id
      do k = 1, molecule(1)%nvib
        write(fout,162) k, molecule(1)%normodes%vibration(k)%freq*cfac, &
                           molecule(2)%normodes%vibration(k)%freq*cfac
      end do
 !   end if

    include 'formats'
158 format(8x,2(f10.3,3x))
!158 format(<isp+2>x,2(f10.3,3x))
   
    !------------------------------------------------------------------------------!
    ! After printing the transformation we can deallocate the auxilyary array KMT  !
    ! This is important for successive calls to cart_transf()                      !
    !------------------------------------------------------------------------------!
    deallocate(KMT)

    return
    end subroutine print_model_transf
   
    !==============================================================================!

	subroutine print_model_transf_file()

	real(kind=dp) :: tmp, sc, f, f1, f2, del, lij
	real(kind=dp), allocatable :: dd(:)
	integer i, j, k

	allocate(dd(1:nmv))
	dd = zero

	! Scrive la matrice di Duschinsky per colonne nel file DUSCH
	open(fdusch,file="DUSCH",status="unknown")
	do j = 1, nmv
		do i = 1, nmv
			write(fdusch,*)trns%JM(i,j)
		end do
	end do

	! Scrive la matrice di Duschinsky Adimensionale per colonne nel file TRASF_DIM
	open(ftransf,file="HAM_DIMLESS",status="unknown")

	tmp = zero
	do k = 1, nmv
		f = molecule(1)%normodes%vibration(k)%freq
		sc = sqrt(cfac*f)
		dd(k) = trns%KM(k)*sc
		tmp = tmp + (dd(k)**2)*f
	end do

	write(ftransf,'(a)')'KAPPA/2'
	write(ftransf,*)tmp/2.0_dp

	write(ftransf,'(a)')'DELTA'
	do j = 1, nmv
		del = zero
		f2 = molecule(2)%normodes%vibration(j)%freq
		do k = 1, nmv
			f1 = molecule(1)%normodes%vibration(k)%freq
			del = del + trns%JM(k,j)*sqrt(f1/f2)*f1*dd(k)
		end do
		write(ftransf,'(a,i2,a,f18.12,a)') '<g dof="',j,'">',del,'</g>'
	end do

	write(ftransf,'(a)')'LAMBDA/2'
	do j = 1, nmv
		f1 = molecule(2)%normodes%vibration(j)%freq
		do i = j, nmv
			f2 = molecule(2)%normodes%vibration(i)%freq
		    !tmp = sqrt(f1/f2)
			!write(ftransf,*)trns%JM(i,j)*tmp
			lij = zero
			do k = 1, nmv
				f = molecule(1)%normodes%vibration(k)%freq
				lij = lij + trns%JM(k,i)*trns%JM(k,j)*(f**2)
			end do
			lij = lij/sqrt(f1*f2)
			if (i == j) lij = lij*0.5_dp
			write(ftransf,'(a,i2,a,i2,a,f18.12,a)') '<g dof="',j,',',i,'">',lij,'</g>'
		end do
	end do

	end subroutine print_model_transf_file

    !==============================================================================!

    subroutine print_model_ham_file()

    real(kind=dp) :: tmp, sc, f, f1, f2, del, delc, lij
    real(kind=dp), allocatable :: DDL(:), KDL(:), JDL(:,:), GDL(:,:), FDL(:,:), wn(:), wo(:)
    integer i, j, k

    allocate(DDL(1:nmv), KDL(1:nmv), wn(1:nmv), wo(1:nmv))
    allocate(JDL(1:nmv,1:nmv))
    allocate(GDL(1:nmv,1:nmv))
    allocate(FDL(1:nmv,1:nmv))

    wn = molecule(1)%normodes%vibration(:)%freq
    wo = molecule(2)%normodes%vibration(:)%freq

    ! Dimensionless Duschinsky matrix
    do i = 1, nmv
        do j = 1, nmv
            JDL(i,j) = trns%JM(i,j)*sqrt(wn(i)/wo(j))
        end do
    end do


    ! Write the Hamiltonian model using dimensionless coordinates
    open(ftransf,file="HAM_DIMLESS",status="unknown")

    write(ftransf,'(a)')'INITIAL STATE HAMILTONIAN'
    write(ftransf,*)0.0_dp,  '#reorganization energy'
    write(ftransf,'(a)')'#POTENTIAL ENERGY: LINEAR TERMS'
    write(ftransf,*) 0
    write(ftransf,'(a)')'#POTENTIAL ENERGY: QUADRATIC TERMS'
    write(ftransf,*)nmv
    do j = 1, nmv
        write(ftransf,'(i3,2x,i3,2x,f18.1)') j, j, wo(j)
    end do
    write(ftransf,'(a)')'#KINETIC ENERGY'
    write(ftransf,*)nmv
    do j = 1, nmv
        write(ftransf,'(i3,2x,i3,2x,f18.1)') j, j, wo(j)
    end do

    write(ftransf,'(a)')'FINAL STATE HAMILTONIAN'

    ! Compute dimensionless displacements and total reorganization energy
    tmp = zero
    do k = 1, nmv
        sc = sqrt(cfac*wn(k))
        DDL(k) = trns%KM(k)*sc
        tmp = tmp + (ddl(k)**2)*wn(k)
    end do

    write(ftransf,*)tmp/2.0_dp, '#reorganization energy'

    ! Write the hamiltonian
    write(ftransf,'(a)')'#POTENTIAL ENERGY: LINEAR TERMS'

    ! Linear couplings: K = 2*J*tr(d*w)
    KDL = matmul(JDL,DDL*wn)
    write(ftransf,*)nmv
    do i = 1, nmv
        write(ftransf,*) i, KDL(i)
    end do

    ! Second order potential energy matrix: F = tr(J)*w*J
    write(ftransf,'(a)')'#POTENTIAL ENERGY: QUADRATIC TERMS'
    write(ftransf,*)nmv*nmv

    FDL = matmul(transpose(JDL),matmul(diag(wn),JDL))
    do i = 1, nmv
        do j = i, nmv
            write(ftransf,'(i3,2x,i3,2x,f18.1)') i, j, FDL(i,j)
        end do
    end do

    ! Kinetic energy term for the second electronic state. G = J^-1*w*tr(J)^-1
    GDL = matmul(invm(JDL),matmul(diag(wn),invm(transpose(JDL))))
    write(ftransf,'(a)')'#KINETIC ENERGY'
    write(ftransf,*)nmv*nmv
    do i = 1, nmv
        do j = i, nmv
            write(ftransf,'(i3,2x,i3,2x,f18.1)') i, j, GDL(i,j)
        end do
    end do

    return
    end subroutine print_model_ham_file

end module model_transformation
