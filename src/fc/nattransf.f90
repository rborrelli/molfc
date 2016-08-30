module natural_internal_transformation

  use fc_h
  use job_type
  use lapack95, only : gesvd

  implicit none

  private
  
  real(kind=dp), private, allocatable :: KMT(:), ishift(:)
  
  type(molecule_t), pointer, private :: molecule(:)
  type(transf_t), pointer, private :: trns
  integer nmv

  public :: nat_intc_transf, print_nat_intc_transf

  contains

!==============================================================================!

    subroutine nat_intc_transf (molecola, transf)
    
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
    end subroutine nat_intc_transf

!==============================================================================!

    subroutine nmshifts
   
    real(kind=dp) deltac, dsb, coeff
    real(kind=dp), allocatable :: B1(:,:), LN1(:,:), ds(:), c1(:,:)
    real(kind=dp), allocatable :: W1(:), V1(:,:), U1(:,:), LN1_inv(:,:)

    integer  i, j, k, ierr, jk, il, ik, nnc, nna, NICO, NNCO
   
    NICO = molecule(1)%intcoord%nintc
    NNCO = molecule(1)%nvib
    !print *, 'NICOxNNCO', nico, nnco
    nna = molecule(1)%structure%numat

    ! Questo ciclo definisce le differenze tra i valori di eq. delle
    ! coordinate interne. Alcuni valori devono essere aggiustati a mano.
    allocate(ds(1:NICO))
    do i = 1, molecule(1)%intcoord%nintc
      ds(i) = molecule(2)%intcoord%ncoord(i)%val - molecule(1)%intcoord%ncoord(i)%val 
      if (molecule(1)%intcoord%ncoord(i)%type == "b" .or. &
            molecule(1)%intcoord%ncoord(i)%type == "d") then
        deltac = zero
        do j = 1, size(molecule(1)%intcoord%ncoord(i)%coord)
            dsb = molecule(2)%intcoord%ncoord(i)%coord(j)%val - molecule(1)%intcoord%ncoord(i)%coord(j)%val 
            if (dsb > pi) dsb = dsb - 2.D0*pi
            if (dsb < -pi) dsb = dsb + 2.D0*pi
            coeff = molecule(1)%intcoord%ncoord(i)%coord(j)%c
            deltac = deltac + coeff*dsb
        end do
        ds(i) = deltac
      end if
    end do

    allocate(B1(1:nico,1:3*nna))
    B1 = get_nat_bmat(molecule(1))

    allocate(LN1(1:NICO,1:NNCO))
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
    print *, 'B singular? :', det(matmul(B1,transpose(B1)))

    write(fout,'(//,2x,''Normal modes in internal coordinates for state |1>'')')
    call layout(LN1,nnco,nnco,nnco,nnco)

    if (NNCO /= NICO) then
        allocate(LN1_inv(1:NNCO,1:NICO))
        allocate(U1(1:NICO,1:NICO)) ! NICO x NICO
        allocate(V1(1:NNCO,1:NNCO)) ! NNCO x NNCO
        allocate(W1(1:NNCO))
 
        !-------------------------------------+
        ! Calcola LN^(-1) con una SVD         |
        !      LN1      = U * W * V'          |
        !      LN1^(-1) = V * 1/W * U'        |
        !-------------------------------------+

        call gesvd(LN1,W1,U1,V1)
     
        ! LN1 ora contiene U.
        do i = 1, NNCO
            LN1_inv(i,:) = U1(:,i) / W1(i)
        end do
        LN1_inv = matmul(transpose(V1),LN1_inv)
        trns%KM = matmul(LN1_inv,ds)
    else
        trns%KM = matmul(invm(LN1),ds)
    end if

    allocate(ishift(1:NICO))
    ishift = ds
    do i = 1, NICO
    write(995,*) ishift(i)
    end do
    deallocate(LN1,B1,c1,ds)
   
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
    real(kind=dp), allocatable :: LNR1(:,:), NR1(:,:)
    real(kind=dp), allocatable :: LNR2(:,:), NR2(:,:)
    integer i, j, k, nnc,nna, NICO, NNCO
   
    NNCO = molecule(1)%nvib
    NICO = molecule(1)%intcoord%nintc
    nna = molecule(1)%structure%numat

    allocate(B1(1:nico,1:3*nna),B2(1:nico,1:3*nna))

    ! Store B matrices into two new arrays.
    ! The B matrices have been calculated by the call to bmat() in MolFC main.
    B1 = get_nat_bmat(molecule(1))
    B2 = get_nat_bmat(molecule(2))
    
    ! Write B matrices according to the Octave ascii format.
    write(22,*)'# name: B1'
    write(22,*)'# type: matrix'
    write(22,*)'# rows:', nico
    write(22,*)'# columns:', 3*nna
    do i = 1, nico
        do j = 1, 3*nna
            write(22,*)B1(i,j)
        end do
    end do
    write(22,*)'# name: B2'
    write(22,*)'# type: matrix'
    write(22,*)'# rows:', nico
    write(22,*)'# columns:', 3*nna
    do i = 1, nico
        do j = 1, 3*nna
            write(22,*)B2(i,j)
        end do
    end do


    allocate(c1(1:3*nna,1:molecule(1)%nvib))
    allocate(c2(1:3*nna,1:molecule(2)%nvib))

    ! C1 e C2 sono coordinate cartesian massa pesate ossia x * sqrt(m)
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

    ! At thi spoint we have: C1 * x = Q1 ; C2 * x = Q2

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

    ! Write L matrices according to the Octave ascii format.
    write(22,*)'# name: L1'
    write(22,*)'# type: matrix'
    write(22,*)'# rows:', nico
    write(22,*)'# columns:', nico
    do i = 1, nico
      write(22,*)(LN1(i,j),j=1,nico)
    end do
    write(22,*)'# name: L2'
    write(22,*)'# type: matrix'
    write(22,*)'# rows:', nico
    write(22,*)'# columns:', nico
    do i = 1, nico
      write(22,*)(LN2(i,j),j=1,nico)
    end do


    !--------------------------+
    !      LN1^(-1) * LN2      |
    !--------------------------+
    trns%JM = matmul(invm(LN1),LN2)

    write(fout,261)
    
    call layout(trns%JM,nmv,nmv,nmv,nmv)

	! Scrive la matrice di Dischinsky per colonne nel file DUSCH
	open(fdusch,file="DUSCH",status="unknown")
	do j = 1, nmv
		do i = 1, nmv
			write(fdusch,*)i, j, trns%JM(i,j)
		end do
	end do

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
   
    subroutine print_nat_intc_transf
    
    real(kind=dp), allocatable  :: vaus(:)
    real(kind=dp) conv
    integer i, j, l, k, icorr, icorr2, icorr3
    integer, allocatable  :: nmixed(:), mmixed(:,:), irec(:)
    integer Kcount, isp
   
    real(kind=dp), parameter :: gtol = 0.96
    real(kind=dp), parameter :: gtolb = 0.10
   
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
    end do
    
    if (molecule(1)%intcoord%nintc > 0) then
      write(fout,159)
      do k = 1, molecule(1)%intcoord%nintc
        do j = 1, size(molecule(1)%intcoord%ncoord(k)%coord)
          select case(molecule(1)%intcoord%ncoord(k)%type)
            case("s")
              isp = 6; conv = one
              write(fout,156) k, molecule(1)%intcoord%ncoord(k)%type, molecule(1)%intcoord%ncoord(k)%coord(j)%list(1:2), &
                              molecule(1)%intcoord%ncoord(k)%coord(j)%val*conv, molecule(2)%intcoord%ncoord(k)%coord(j)%val*conv
            case("t")
              isp = 6; conv = 180.d0/pi
              write(fout,156) k, molecule(1)%intcoord%ncoord(k)%type, molecule(1)%intcoord%ncoord(k)%coord(j)%list(1:2), &
                              molecule(1)%intcoord%ncoord(k)%coord(j)%val*conv, molecule(2)%intcoord%ncoord(k)%coord(j)%val*conv
            case("b","l")
              isp = 9; conv = 180.d0/pi
              write(fout,1561) k, molecule(1)%intcoord%ncoord(k)%type, molecule(1)%intcoord%ncoord(k)%coord(j)%list, &
                              molecule(1)%intcoord%ncoord(k)%coord(j)%val*conv, molecule(2)%intcoord%ncoord(k)%coord(j)%val*conv
            case("w","d")
              isp = 12; conv = 180.d0/pi
              !write(fout,157,advance="no") k, molecule(1)%intcoord%ncoord(k)%type, molecule(1)%intcoord%ncoord(k)%coord(j)%list
              write(fout,1562) k, molecule(1)%intcoord%ncoord(k)%type, molecule(1)%intcoord%ncoord(k)%coord(j)%list, &
                              molecule(1)%intcoord%ncoord(k)%coord(j)%val*conv, molecule(2)%intcoord%ncoord(k)%coord(j)%val*conv
          end select
        end do
        write(fout,158) conv*molecule(1)%intcoord%ncoord(k)%val,conv*molecule(2)%intcoord%ncoord(k)%val,ishift(k)*conv
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
156 format(2x,i3,2x,a,3x,2(i2,1x),2(f10.3,3x))
1561 format(2x,i3,2x,a,3x,3(i2,1x),2(f10.3,3x))
1562 format(2x,i3,2x,a,3x,4(i2,1x),2(f10.3,3x))
158 format(<isp+11>x,3(f10.3,3x))
   
	! Free KMT 
	deallocate(KMT)

    return
    end subroutine print_nat_intc_transf

end module natural_internal_transformation
