module internal_transformation

  use fc_h
  use job_type
  use iofiles, only : newfile

  implicit none

  private
  
  real(kind=dp), private, allocatable :: KMT(:), ishift(:)
  
  type(molecule_t), pointer, private :: molecule(:)
  type(transf_t), pointer, private :: trns
  integer nmv

  public :: intc_transf, print_intc_transf, print_intc_transf_file, print_intc_ham_file

  contains

!==============================================================================!

    subroutine intc_transf (molecola, transf)
    
    type(molecule_t), target, intent(in) :: molecola(1:2)
    type(transf_t), target, intent(inout) :: transf

    integer i, j, ij, ierr, jj

    ! molecule e trns sono accessibili a tutto il modulo.
    molecule => molecola
    trns => transf

    ! Numero totale di vibrazioni per stato
    nmv = molecule(1)%nvib
    !------------------------------------+
    ! Alloca le variabili del modulo FC.
    !------------------------------------+
    allocate (trns%JM(1:nmv,1:nmv), trns%KM(1:nmv), KMT(1:nmv))
    !----------------------------------------------------------------+
    ! A questo punto le geometrie sono quelle finali e le matrici
    ! dei modi normali sono state riordinate in base alle geometrie.
    !
    ! Calcola la matrice di Duschinsky JM
    !
	if (transf%dusch) then
		call nmrot()
	else
		! set duschinsky matrix to identity
		trns%JM = zero
		forall (i=1:nmv) trns%JM(i,i) = one
	end if
    !
    ! Calcola il determinante della trasformazione.
    !
    call metric_determinant()
    !
    ! Calcola il vettore di spostamenti dei modi normali GK.
    ! La geometria del secondo stato elettronico e' stata trasformata da AXSMT. 
    !
    call nmshifts()

    include 'formats'

    return
    end subroutine intc_transf

!==============================================================================!

    subroutine nmshifts
   
    real(kind=dp) deltax
    real(kind=dp), allocatable :: B1(:,:), LN1(:,:), ds(:), c1(:,:)
    real(kind=dp), allocatable :: W1(:), V1(:,:), U1(:,:), LN1_inv(:,:)
    real(kind=dp), allocatable :: work(:)
    integer :: lwork, info

    integer  i, j, k, ierr, jk, il, ik, nnc, nna, NICO, NNCO
    external dgesvd
   
    NICO = molecule(1)%intcoord%nintc
    NNCO = molecule(1)%nvib
    nna = molecule(1)%structure%numat

    ! Questo ciclo definisce le differenze tra i valori di eq. delle
    ! coordinate interne. Alcuni valori devono essere aggiustati a mano.
    allocate(ds(1:NICO))
    do i = 1, molecule(1)%intcoord%nintc
      ds(i) = molecule(2)%intcoord%coord(i)%val - molecule(1)%intcoord%coord(i)%val 
      if (molecule(1)%intcoord%coord(i)%type == "b" .or. &
          molecule(1)%intcoord%coord(i)%type == "d") then
          if (ds(i) > pi) ds(i) = ds(i) - 2.D0*pi
          if (ds(i) < -pi) ds(i) = ds(i) + 2.D0*pi
      end if
      !if (molecule(1)%intcoord%coord(i)%type == "d") ds(i) = 0.0_dp ! no dihedrals
    end do

    allocate(B1(1:nico,1:3*nna))
    B1 = get_bmat(molecule(1))

    allocate(LN1(1:NICO,1:NNCO)) ! NICO x NNCO
    allocate(c1(1:3*molecule(1)%structure%numat,1:molecule(1)%nvib))

    c1 = get_normal_modes(molecule(1))
    do i = 1, molecule(1)%structure%numat
      c1(3*(i-1)+1:3*i,:) = c1(3*(i-1)+1:3*i,:)/sqrt(molecule(1)%structure%atom(i)%elem%AM)
    end do
    !----------------------+
    ! LN1 = B1*M^(-1/2)*C1 |
    !----------------------+
    LN1 = zero
    LN1 = matmul(B1,c1)

	write(fout,'(//,2x,''Normal modes in internal coordinates for state |1>'')')
    call layout(LN1,nico,nnco,nico,nnco)

    allocate(LN1_inv(1:NNCO,1:NICO))
    allocate(U1(1:NICO,1:NICO)) ! NICO x NICO
    allocate(V1(1:NNCO,1:NNCO)) ! NNCO x NNCO
    allocate(W1(1:NNCO))
 
    !-------------------------------------+
    ! Calcola LN^(-1) con una SVD         |
    !      LN1      = U * W * V'          |
    !      LN1^(-1) = V * 1/W * U'        |
    !-------------------------------------+

    ! Query size of WORK array
    allocate(work(1:NICO))
    call dgesvd('A','A',NICO,NNCO,LN1,NICO,W1,U1,NICO,V1,NNCO,WORK,-1,INFO)
    lwork = work(1)
    deallocate(work)
    allocate(work(1:lwork))
    call dgesvd('A','A',NICO,NNCO,LN1,NICO,W1,U1,NICO,V1,NNCO,WORK,LWORK,INFO)
    !call svdcmp(LN1,W1,V1)
     
    ! LN1 ora contiene U.
    do i = 1, NNCO
      !LN1_inv(i,:) = LN1(:,i) / W1(i)
      LN1_inv(i,:) = U1(:,i) / W1(i)
    end do

    LN1_inv = matmul(transpose(V1),LN1_inv)

    do i = 1, NNCO
      do j = 1, NICO
        write(444,*) i, j, LN1_inv(i,j)*ds(j)
      end do
    end do
    trns%KM = matmul(LN1_inv,ds)

    allocate(ishift(1:NICO))
    ishift = ds
    deallocate(LN1,LN1_inv,B1,c1,ds)
   
    KMT = -matmul(invm(trns%JM),trns%KM)

    return
    end subroutine nmshifts
    
!==============================================================================!

    subroutine nmrot

    !------------------------------------------------------------------+  
    ! Calcola la matrice di trasfomazione tra i modi normali.          |
    ! Se si e' fatto l'axis switching T2 e' stata gia' ruotata con S0  |
    ! da axsmt.                                                        |
    !------------------------------------------------------------------+  
 
    real(kind=dp), allocatable :: B1(:,:), B2(:,:), B1p(:,:), B2p(:,:)
    real(kind=dp), allocatable :: lw(:,:), ds(:), c1(:,:), c2(:,:), gm1(:,:)
    real(kind=dp), allocatable :: LN1(:,:), LN2(:,:)
    real(kind=dp), allocatable :: LN1_inv(:,:), V1(:,:), W1(:), U1(:,:)
    real(kind=dp), allocatable :: LN2_inv(:,:), V2(:,:), W2(:), U2(:,:)
    real(kind=dp), allocatable :: work(:)
    integer i, j, k, nnc,nna, NICO, NNCO, lwork, info
   
    external dgesvd

    NNCO = molecule(1)%nvib
    NICO = molecule(1)%intcoord%nintc
    nna = molecule(1)%structure%numat

    allocate(B1(1:nico,1:3*nna),B2(1:nico,1:3*nna))

    ! Store B matrices into two new arrays.
    ! The B matrices have been calculated by the call to bmat() in MolFC main.
    B1 = get_bmat(molecule(1))
    B2 = get_bmat(molecule(2))
    
    ! Write B matrices according to the Octave ascii format.
    write(22,*)'# name: B1'
    write(22,*)'# type: matrix'
    write(22,*)'# rows:', nico
    write(22,*)'# columns:', 3*nna
    do i = 1, nico
      write(22,*)(B1(i,j),j=1,3*nna)
    end do
    write(22,*)'# name: B2'
    write(22,*)'# type: matrix'
    write(22,*)'# rows:', nico
    write(22,*)'# columns:', 3*nna
    do i = 1, nico
      write(22,*)(B2(i,j),j=1,3*nna)
    end do


    allocate(c1(1:3*nna,1:molecule(1)%nvib))
    allocate(c2(1:3*nna,1:molecule(2)%nvib))

    c1 = get_normal_modes(molecule(1))
    c2 = get_normal_modes(molecule(2))

    !--------------------+
    ! C2 <= M^(-1/2)*C2  |
    !--------------------+
    do i = 1, molecule(2)%structure%numat
      c2(3*(i-1)+1:3*i,:) = c2(3*(i-1)+1:3*i,:)/sqrt(molecule(2)%structure%atom(i)%elem%AM)
    end do
    !--------------------+
    ! C1 <= M^(-1/2)*C1  |
    !--------------------+
    do i = 1, molecule(2)%structure%numat
      c1(3*(i-1)+1:3*i,:) = c1(3*(i-1)+1:3*i,:)/sqrt(molecule(1)%structure%atom(i)%elem%AM)
    end do

    allocate(LN1(1:NICO,1:NNCO))
    allocate(LN2(1:NICO,1:NNCO))
    !----------------------+b
    ! LN2 = B2*M^(-1/2)*C2 |
    !----------------------+
    LN2 = zero
    LN2 = matmul(B2,c2)
    !----------------------+
    ! LN1 = B1*M^(-1/2)*C1 |
    !----------------------+
    LN1 = zero
    LN1 = matmul(B1,c1)

	write(fout,'(//,2x,''Normal modes in internal coordinates for state |1>'')')
    call layout(LN1,nico,nnco,nico,nnco)
	write(fout,'(//,2x,''Normal modes in internal coordinates for state |2>'')')
    call layout(LN2,nico,nnco,nico,nnco)

    allocate(LN1_inv(1:NNCO,1:NICO))
    allocate(V1(1:NNCO,1:NNCO))
    !allocate(U1(1:NICO,1:NNCO))
    allocate(U1(1:NICO,1:NICO))
    allocate(W1(1:NNCO))

    !-------------------------------------+
    ! Calcola LN^(-1) con una SVD         |
    !      LN1      = U * W * V'          |
    !      LN1^(-1) = V * 1/W * U'        |
    !-------------------------------------+

    ! Query size of WORK array
    allocate(work(1:NICO))
    call dgesvd('A','A',NICO,NNCO,LN1,NICO,W1,U1,NICO,V1,NNCO,WORK,-1,INFO)
    lwork = work(1)
    deallocate(work)
    allocate(work(1:lwork))
    call dgesvd('A','A',NICO,NNCO,LN1,NICO,W1,U1,NICO,V1,NNCO,WORK,LWORK,INFO)

    ! LN1 ora contiene U.
    ! Qui devo controllare che W1(i) non sia nullo, in tal caso meglio usare il metodo 2 (vedi sopra).
    do i = 1, NNCO
      !LN1_inv(i,:) = LN1(:,i) / W1(i)
      LN1_inv(i,:) = U1(:,i) / W1(i)
    end do

    LN1_inv = matmul(transpose(V1),LN1_inv)

    !--------------------------+
    !      LN1^(-1) * LN2      |
    !--------------------------+
    trns%JM = matmul(LN1_inv,LN2) 

    write(fout,261)
    
    call layout(trns%JM,nmv,nmv,nmv,nmv)

    include 'formats'
    return
    end subroutine nmrot
   
    !==============================================================================!

	subroutine print_intc_transf_file()

	real(kind=dp) :: tmp, sc, f, f1, f2, del, delc, lij
	real(kind=dp), allocatable :: dd(:)
	integer i, j, k

	allocate(dd(1:nmv))
	dd = zero

	! Scrive la matrice di Dischinsky per colonne nel file DUSCH
	open(fdusch,file="DUSCH",status="unknown")
	do j = 1, nmv
		do i = 1, nmv
			write(fdusch,*)i, j, trns%JM(i,j)
		end do
	end do

	! Scrive la matrice di Duschinsky adimensionale per colonne nel file DUSCH_DIM
	open(ftransf,file="HAM_DIMLESS",status="unknown")

	tmp = zero
	do k = 1, nmv
		f = molecule(1)%normodes%vibration(k)%freq
		sc = sqrt(cfac*f)
		dd(k) = trns%KM(k)*sc
		tmp = tmp + (dd(k)**2)*f
	end do

        write(ftransf,'(a)')'STATE 1'
	do j = 1, nmv
		delc = molecule(2)%normodes%vibration(j)%freq
		write(ftransf,'(a)') '<g>'
		write(ftransf,'(a)') '<!--g*(p^2+q^2)/2 -->'
		write(ftransf,'(a,f18.12,a)') '<value> ',delc,' </value>'
		write(ftransf,'(a,i2,a)') '<dof="',j,'" optype="9" />'
		write(ftransf,'(a)') '</g>'
	end do

        write(ftransf,'(a)')'STATE 2'

!	write(ftransf,'(a)')'KAPPA/2'
!	write(ftransf,*)tmp/2.0_dp
!	write(ftransf,'(a)')'DELTA'

        ! Kinetic energy term for the second electronic state. I have to check if the transformation
        ! of the Hamiltonian is correct!!!! There should be some correction due to the non-orthogonality of the Duschinsky matrix.
	do j = 1, nmv
		delc = molecule(2)%normodes%vibration(j)%freq
		write(ftransf,'(a)') '<g>'
		write(ftransf,'(a)') '<!--g*(p^2) -->'
		write(ftransf,'(a,f18.12,a)') '<value> ',delc,' </value>'
		write(ftransf,'(a,i2,a)') '<dof="',j,'" optype="2" />'
		write(ftransf,'(a)') '</g>'
	end do

	do j = 1, nmv
		delc = zero
		f2 = molecule(2)%normodes%vibration(j)%freq
		do k = 1, nmv
			f = molecule(1)%normodes%vibration(k)%freq
		    ! the line below is Wrong but I use for some results already published.
			!del = del + trns%JM(j,k)*f*dd(k) !JM is not for dimensionless coordinates!!
			! Correct (this is the correct implementation)
			delc = delc + trns%JM(k,j)*sqrt(f2/f)*f*dd(k)
		end do
		write(ftransf,'(a)') '<g>'
		write(ftransf,'(a)') '<!--g*(q) -->'
		write(ftransf,'(a,f18.12,a)') '<value> ',delc,' </value>'
		write(ftransf,'(a,i2,a)') '<dof="',j,'" optype="0" pow="1" />'
		write(ftransf,'(a)') '</g>'
	end do

	write(ftransf,'(a)')'Diagonal elements:  Hessian/2; off diagonal: Hessian'
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
			if (i == j) lij=lij*0.5_dp
                        if (i == j) then
			write(ftransf,'(a)') '<g>'
		        write(ftransf,'(a)') '<!--g*(q^2) -->'
		        write(ftransf,'(a,f18.12,a)') '<value> ',lij,' </value>'
			write(ftransf,'(a,i2,a)') '<dof="',i,'" optype="0" pow="2" />'
			write(ftransf,'(a)') '</g>'
                        else
			write(ftransf,'(a)') '<g>'
		        write(ftransf,'(a)') '<!--g*(q_i q_j) -->'
		        write(ftransf,'(a,f18.12,a)') '<value> ',lij,' </value>'
			write(ftransf,'(a,i2,a)') '<dof="',j,'" optype="0" pow="1" />'
			write(ftransf,'(a,i2,a)') '<dof="',i,'" optype="0" pow="1" />'
			write(ftransf,'(a)') '</g>'
                        end if
		end do
	end do

	end subroutine print_intc_transf_file

    !==============================================================================!

    subroutine print_intc_ham_file()

    real(kind=dp) :: tmp, sc, f, f1, f2, del, delc, lij
    real(kind=dp), allocatable :: DDL(:), KDL(:), JDL(:,:), GDL(:,:), FDL(:,:), wn(:), wo(:)
    integer i, j, k, ioff

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

    ! ioff: is used as an offset number to print the final Hamiltonian.
    ! It can be used when we have two or more molecules and we have to number all
    ! the vibrational modes.
    ioff = 0
    ! Write the Hamiltonian model using dimensionless coordinates
    open(ftransf,file="HAM_DIMLESS",status="unknown")

    write(ftransf,'(a)')'INITIAL STATE HAMILTONIAN'
    write(ftransf,*)0.0_dp,  '#reorganization energy'
    write(ftransf,'(a)')'#POTENTIAL ENERGY: LINEAR TERMS'
    write(ftransf,*) 0
    write(ftransf,'(a)')'#POTENTIAL ENERGY: QUADRATIC TERMS'
    write(ftransf,*)nmv
    do j = 1, nmv
        write(ftransf,'(i3,2x,i3,2x,f18.1)') j+ioff, j+ioff, wo(j)
    end do
    write(ftransf,'(a)')'#KINETIC ENERGY'
    write(ftransf,*)nmv
    do j = 1, nmv
        write(ftransf,'(i3,2x,i3,2x,f18.1)') j+ioff, j+ioff, wo(j)
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

    ! Linear couplings: K = 2*J*(d*w)^tr
    KDL = matmul(JDL,DDL*wn)
    write(ftransf,*)nmv
    do i = 1, nmv
        write(ftransf,*) i+ioff, KDL(i)
    end do

    ! Second order potential energy matrix: F = J^(tr)*w*J
    write(ftransf,'(a)')'#POTENTIAL ENERGY: QUADRATIC TERMS'
    write(ftransf,*)nmv*(nmv+1)/2

    FDL = matmul(transpose(JDL),matmul(diag(wn),JDL))
    do i = 1, nmv
        do j = i, nmv
            write(ftransf,'(i3,2x,i3,2x,f18.1)') i+ioff, j+ioff, FDL(i,j)
        end do
    end do

    ! Kinetic energy term for the second electronic state. G = J^-1*w*tr(J)^-1
    GDL = matmul(invm(JDL),matmul(diag(wn),invm(transpose(JDL))))
    write(ftransf,'(a)')'#KINETIC ENERGY'
    write(ftransf,*)nmv*(nmv+1)/2
    do i = 1, nmv
        do j = i, nmv
            write(ftransf,'(i3,2x,i3,2x,f18.1)') i+ioff, j+ioff, GDL(i,j)
        end do
    end do

    return
    end subroutine print_intc_ham_file

!==============================================================================!
   
    subroutine metric_determinant 
   
    real(kind=dp) :: detgm
    detgm = det(trns%JM)
   
    write(fout,10) detgm
   
    include 'formats'

    return        
    end subroutine metric_determinant 

!==============================================================================!
   
    subroutine print_intc_transf
    
    real(kind=dp), allocatable  :: vaus(:)
    real(kind=dp) conv
    integer i, j, l, k, icorr, icorr2, icorr3
    integer, allocatable  :: nmixed(:), mmixed(:,:), irec(:)
    integer Kcount, isp
   
    real(kind=dp), parameter :: gtol = 0.96
    real(kind=dp), parameter :: gtolb = 0.10
   
    write(fout,'(/,2x,a)') "Writing normal modes shift into file:"

    call newfile(233,"shift",extension=".dat")

    write(fout,140)
    
    write(fout,2401) molecule(1)%id
    write(fout,'(/)')
    write(fout,2413)
    ! New output version. Normal modes mixing is still missing.
    do k = 1, molecule(1)%nvib
      write(fout,161) k, molecule(1)%normodes%vibration(k)%freq, &
                         molecule(2)%normodes%vibration(k)%freq, &
                         trns%KM(k)*sqrt(cfac*molecule(1)%normodes%vibration(k)%freq), & 
                         KMT(k)*sqrt(cfac*molecule(2)%normodes%vibration(k)%freq), & 
                         trns%KM(k), KMT(k)
      write(233,'(i3,2x,5(g15.8,2x))') k, molecule(1)%normodes%vibration(k)%freq, molecule(2)%normodes%vibration(k)%freq,  &
      (trns%KM(k)*sqrt(cfac*molecule(1)%normodes%vibration(k)%freq)), &
      molecule(1)%normodes%vibration(k)%freq*(trns%KM(k)*sqrt(cfac*molecule(1)%normodes%vibration(k)%freq)), &
      molecule(1)%normodes%vibration(k)%freq*(trns%KM(k)*sqrt(cfac*molecule(1)%normodes%vibration(k)%freq))**2/2.0_dp
      !write(233,*)  trns%KM(k)
    end do
    
    if (molecule(1)%intcoord%nintc > 0) then
      write(fout,159)
      write(fout,*)  '          FINAL STATE    INITIAL STATE '
      do k = 1, molecule(1)%intcoord%nintc
        select case(molecule(1)%intcoord%coord(k)%type)
          case("s")
            isp = 6; conv = one
            write(fout,157,advance="no") k, molecule(1)%intcoord%coord(k)%type, molecule(2)%intcoord%coord(k)%list(1:2)
          case("t")
            isp = 6; conv = 180.d0/pi
            write(fout,157,advance="no") k, molecule(1)%intcoord%coord(k)%type, molecule(2)%intcoord%coord(k)%list(1:2)
          case("b","l")
            isp = 3; conv = 180.d0/pi
            write(fout,157,advance="no") k, molecule(1)%intcoord%coord(k)%type, molecule(2)%intcoord%coord(k)%list
          case("w","d")
            isp = 0; conv = 180.d0/pi
            write(fout,157,advance="no") k, molecule(1)%intcoord%coord(k)%type, molecule(2)%intcoord%coord(k)%list
        end select
        write(fout,158) conv*molecule(1)%intcoord%coord(k)%val,conv*molecule(2)%intcoord%coord(k)%val,ishift(k)*conv
      end do
    end if
    deallocate (ishift)

    write(fout,163)
    write(fout,2401) molecule(1)%id
    do k = 1, molecule(1)%nvib
        write(fout,162) k, molecule(1)%normodes%vibration(k)%freq*cfac, &
                           molecule(2)%normodes%vibration(k)%freq*cfac
    end do

    include 'formats'
158 format(10(x),3(f10.3,3x))
   
	! Free KMT 
	deallocate(KMT)

    return
    end subroutine print_intc_transf

end module internal_transformation
