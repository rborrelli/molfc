module cartesian_transformation

  use fc_h
  use job_type
  
  implicit none

  private
  
  real(kind=dp), private, allocatable :: KMT(:)
  
  type(molecule_t), pointer, private :: molecule(:)
  type(transf_t), pointer, private :: trns
  integer nmv

  public :: cart_transf, print_cart_transf, print_cart_transf_file


  contains

    !==============================================================================!

    subroutine cart_transf (molecola, transf)
    
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
    allocate (transf%JM(1:nmv,1:nmv), transf%KM(1:nmv), KMT(1:nmv))
    !----------------------------------------------------------------+
    ! A questo punto le geometrie sono quelle finali e le matrici
    ! dei modi normali sono state riordinate in base alle geometrie.
    !
    ! Calcola la matrice di Duschinsky GM
    !
	if (trns%dusch)  then
		call nmrot()
	else
		! set duschinsky matrix to identity
		trns%JM = zero
		forall (i=1:nmv) trns%JM(i,i) = one
	end if

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
    end subroutine cart_transf

    !==============================================================================!

    subroutine nmshifts 
   
    real(kind=dp) deltax, sh
    real(kind=dp), allocatable :: KSH(:)

    integer  i, j, k
   
    !----------------------------------------------------+
    ! Calcolo degli spostamenti in coordinate cartesiane |
    ! KSH = tr(T1)(x02 - x01)
    !----------------------------------------------------+
    !allocate(KSH(1:nmv))
    !KSH = zero
    deltax = zero
    trns%KM = zero
    do i = 1, nmv
        sh = zero
        do j = 1, molecule(1)%structure%numat
	        do k = 1, 3
	            ! deltax is the mass-weighted displacement of atom j along direction k
	           deltax = (molecule(2)%structure%atom(j)%coord(k) - &
	                     molecule(1)%structure%atom(j)%coord(k)) * &
	                     sqrt(molecule(1)%structure%atom(j)%elem%AM)
				trns%KM(i) = trns%KM(i) + molecule(1)%normodes%vibration(i)%atom(j)%d(k)*deltax
	        end do
        end do
    end do


    !----------------------------------------------------+
    ! Calcolo degli spostamenti in coordinate cartesiane |
    ! KMT = tr(T2)(x01 - x02)                            |
    !----------------------------------------------------+
    deltax = zero
    KMT = zero
    do i = 1, nmv
      do j = 1, molecule(1)%structure%numat 
        do k = 1, 3
          deltax = (molecule(1)%structure%atom(j)%coord(k) - &
                    molecule(2)%structure%atom(j)%coord(k)) * &
                    sqrt(molecule(2)%structure%atom(j)%elem%AM)
          KMT(i) = KMT(i) + molecule(2)%normodes%vibration(i)%atom(j)%d(k)*deltax
        end do
      end do
    end do

    !-------------------------------------------------------+
    ! Calcolo degli spostamenti in coordinate adimensionali |
    !-------------------------------------------------------+
    !f1 = molecule(1)%normodes%vibration(:)%freq * cfac
    !f2 = molecule(2)%normodes%vibration(:)%freq * cfac
	!trnsdl%KM = trns%KM*sqrt(f1)
    !forall(i=1:nmv,j=1:nmv) trnsdl%JM(i,j) = trns%JM(i,j)*sqrt(f1(i)/f2(j))
    !forall(i=1:nmv,j=1:nmv) trnsdl%JM(i,j) = trns%JM(i,j)*sqrt(f1(i)/f2(j))

    return
    end subroutine nmshifts 
    
    !==============================================================================!

    subroutine nmrot 

    !------------------------------------------------------------------+  
    ! Calcola la matrice di trasfomazione tra i modi normali.          |
    ! Se si e' fatto l'axis switching T2 e' stata gia' ruotata con S0  |
    ! da axsmt.                                                        |
    !------------------------------------------------------------------+  
 
	real(kind=dp) :: tmp
    integer i, j, k, nnc,nna
   
    !----------------------------------------------------+
    ! Calcolo della matrice di rotazione in coord. cart. |
    ! Viene calcolata:                                   | 
    !     Q(1)_mu = sum_nu JM_{mu,mu} Q(2)_nu            | 
    ! con                                                |
    !     JM = T1^+ T2                                   |
    !----------------------------------------------------+
    trns%JM = zero
    do i = 1, nmv
      do j = 1, nmv
        do k = 1, molecule(1)%structure%numat
          trns%JM(i,j) = trns%JM(i,j) + dot_product(molecule(1)%normodes%vibration(i)%atom(k)%d(1:3),&
                              molecule(2)%normodes%vibration(j)%atom(k)%d(1:3))
        end do
      end do
    end do

    ! Arbitrario: mette a zero tutti i numeri minori di 1e-5
    !where(abs(trns%JM) < 1e-5) trns%JM = zero

    write(fout,261)
    ! Scrive la matrice di Duschinsky in output.
    call layout(trns%JM,nmv,nmv,nmv,nmv)

    include 'formats'
    return
    end subroutine nmrot

    !==============================================================================!

	subroutine print_cart_transf_file()

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

	write(ftransf,'(a)')'KAPPA/2'
	write(ftransf,*)tmp/2.0_dp

	write(ftransf,'(a)')'DELTA'
	do j = 1, nmv
		del = zero
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
		write(ftransf,'(a,f18.12,a)') '<value> ',delc,' </value>'
		write(ftransf,'(a,i2,a)') '<dof="',j,'" optype="9" />'
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
		        write(ftransf,'(a,f18.12,a)') '<value> ',lij,' </value>'
			write(ftransf,'(a,i2,a)') '<dof="',i,'" optype="0" pow="2" />'
			write(ftransf,'(a)') '</g>'
                        else
			write(ftransf,'(a)') '<g>'
		        write(ftransf,'(a,f18.12,a)') '<value> ',lij,' </value>'
			write(ftransf,'(a,i2,a)') '<dof="',j,'" optype="0" pow="1" />'
			write(ftransf,'(a,i2,a)') '<dof="',i,'" optype="0" pow="1" />'
			write(ftransf,'(a)') '</g>'
                        end if
		end do
	end do

	end subroutine print_cart_transf_file

    !==============================================================================!
   
    subroutine metric_determinant 
   
    real(kind=dp) :: detgm
    detgm = det(trns%JM)
   
    write(fout,10) detgm
   
    include 'formats'

    return        
    end subroutine metric_determinant 

    !==============================================================================!        
   
    subroutine print_cart_transf
    
    real(kind=dp), allocatable  :: vaus(:)
    real(kind=dp) conv
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
      write(fout,161) k, molecule(1)%normodes%vibration(k)%freq, &
                         molecule(2)%normodes%vibration(k)%freq, &
                         trns%KM(k)*sqrt(cfac*molecule(1)%normodes%vibration(k)%freq), & 
                         KMT(k)*sqrt(cfac*molecule(2)%normodes%vibration(k)%freq), & 
                         trns%KM(k), KMT(k)
    end do
    
    if (molecule(1)%intcoord%nintc > 0) then
      write(fout,159)
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
        write(fout,158) conv*molecule(1)%intcoord%coord(k)%val,conv*molecule(2)%intcoord%coord(k)%val 
      end do
    end if


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
   
    !------------------------------------------------------------------------------!
    ! After printing the transformation we can deallocate the auxilyary array KMT  !
    ! This is important for successive calls to cart_transf()                      !
    !------------------------------------------------------------------------------!
    deallocate(KMT)

    return
    end subroutine print_cart_transf
   

end module cartesian_transformation

