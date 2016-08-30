module fc

  use fc_h
 
  implicit none

  private
 
  public :: fc_init, fc_end, MDFC_00

  integer, private :: nmv

  type(molecule_t), pointer, private :: molecule(:)
  type(transf_t), pointer, private :: trns
  real(kind=dp), allocatable, private :: gm(:,:), gk(:)
  real(kind=dp), allocatable :: B(:,:), A(:,:), X(:,:)
  real(kind=dp), allocatable, public :: d(:)
  real(kind=dp), allocatable :: f1(:), f2(:)
  type(incvib_t) :: incvib

  contains

	!==============================================================================!

    subroutine fc_init (molecola,transf,incvibr)
    
    type(molecule_t), target, intent(in) :: molecola(1:2)
    type(transf_t), target, intent(inout) :: transf
    type(incvib_t), intent(in) :: incvibr
    integer i, j, ij, ierr, jj, alloc_err
	real(kind=dp) :: theta, deg
	real(kind=dp), allocatable :: YY(:,:)
	
    ! molecule e' accessibile a tutto il modulo.
    molecule => molecola
    trns => transf
    !    allocate (incvib%id(1:molecola(1)%nvib))
  	! incvib: accessibile a tutto il modulo.
   	incvib = incvibr ! hidden automatic allocation: beware!

    ! Numero totale di vibrazioni. Con questa scelta includo tutte le vibrazioni
    ! di una molecola nel calcolo delle matrici M e Q (usate nelle formule di ricorrenza).
    !nmv = molecule(1)%nvib

    ! Con questa soluzione posso selezionare le vibrazioni incluse nel calcolo fc.
    ! Queste sono quelle incluse nell'array incvib definito in input dal tag <include>.
    nmv = size(incvib%id)
    !print *, 'NMV ', nmv, incvib%id
    if (allocated(gm)) deallocate(gm)
    if (allocated(gk)) deallocate(gk)
    allocate(gm(1:nmv,1:nmv), gk(1:nmv), stat=alloc_err)
    if (alloc_err /= 0) ierr = error(0,"Cannot allocate Duschinksy matrix")

    ! NOTE: It would be probably better to use gm and gk as pointers and set 
    ! gm => transf%JM(incvib,incvib)
    ! gk => transf%KM(incvib)
    ! However the solution adopted below  is not a serious memory drawback since
    ! gm an gk are not so big and are deallocated after the determination of M
    ! and Q.
    gm = transf%JM(incvib%id,incvib%id)
    gk = transf%KM(incvib%id)
    ! actually used frequencies.
    allocate(f1(1:nmv), f2(1:nmv))
    f1 = molecule(1)%normodes%vibration(incvib%id(:))%freq * cfac ! FINAL_STATE frequencies
    f2 = molecule(2)%normodes%vibration(incvib%id(:))%freq * cfac ! INITIAL_STATE frequencies

    ! Trasforma la mia G nella matrice di trasformazione in coordinate adimensionali
    !  J = W1^{-1/2}*G*W2^{1/2} 
    ! dove W1 e W2 sono le matrici diagonali con le frequenze.
    forall(i=1:nmv,j=1:nmv) gm(i,j) = gm(i,j)*sqrt(f1(i)/f2(j))
	if (.false.) then
    	write(fout,'(2x,a)')'Duschinsky matrix'
    	call layout(gm,nmv,nmv,nmv,nmv)
    	write(fout,'(/)')
	end if
    ! Trasforma GK in delta di Watson e lo mette in d: 
    ! d: lo spostamento delle coordinate adimensionali
    allocate(d(1:nmv))
    d = sqrt(f1)*GK
    if (transf%model) then
        if (transf%modk%type == "dimless") then
            d = GK
        else
            d = sqrt(f1)*GK
        end if
    end if
	!-----------------------------------------------------+
	! Calcola le matrici A, B e X e le rispettive inverse |
	!-----------------------------------------------------+	
    call basic_matrices
    !------------------------------------------------------------+
    ! Calcola le matrici utilizzate nelle formule di ricorrenza. |
    !------------------------------------------------------------+
    ! Con questa chiamata vengono alcolate le matrici M e Q che
    ! verranno utilizzate in molti altri moduli del programma.  
    !-----------------------------------------------------------+    
    call recursion_matrices

    include 'formats'

    return
    end subroutine fc_init

	!==============================================================================!

    subroutine basic_matrices 

    ! Calcola B, A e X
    ! utilizza come input GM, ossia la matrice di Duschinsky.
    ! NOTE: A and B are symmetric matrices.

    real(kind=dp) :: temp
    integer i, j, k
	
    ! A = 0.5*(I+tr(J)*J)
    allocate(A(1:nmv,1:nmv))
    A = zero
    do i = 1, nmv
      do j = i, nmv
        A(i,j) = dot_product(GM(1:nmv,i),GM(1:nmv,j))
        if (i /= j) A(j,i) = A(i,j)
      end do
    end do
    forall (i=1:nmv) A(i,i) = A(i,i) + 1.0_dp
    A = 0.5_dp*A

    ! B = 0.5*(I+J*tr(J))
    allocate(B(1:nmv,1:nmv))
    B = zero
    do i = 1, nmv
      do j = i, nmv
        B(i,j) = dot_product(GM(i,1:nmv),GM(j,1:nmv))
        if (i /= j) B(j,i) = B(i,j)
      end do
    end do
    forall (i=1:nmv) B(i,i) = B(i,i) + 1.0_dp
    B = 0.5_dp*B

    ! We need the inverse of A and B
    ! B <-- B^{-1}; A <-- A^{-1}
    A = invm(A) ! Inversa di A
    B = invm(B) ! Inversa di B

	! X = 0.5*(tr(J)+J^{-1})
    allocate(X(1:nmv,1:nmv))
    X = zero
    X = 0.5_dp*(transpose(GM) + invm(GM))
  	! X <-- X^{-1}
    X = invm(X)

    return
    end subroutine basic_matrices 

	!==============================================================================!        

    subroutine recursion_matrices 

    ! NOTA BENE: M E Q SONO LE MATRICI DA CUI DIPENDONO LE FORMULE DI
    ! RICORRENZA DEI FC. LA DEFINIZIONE ADOTTATA IN QUESTO PROGRAMMA PREVEDE
    ! CHE LE PRIME N VIBRAZIONI SIANO ASSOCIATE ALLO STATO INIZIALE E LE SECONDE
    ! N SIANO ASSOCIATE ALLO STATO FINALE. WATSON, NEL SUO ARTICOLO SUL CAN J.
    ! FA L'OPPOSTO, QUINDI LE SUE DEFINIZIONI SONO DIVERSE.
    ! ---------------------------------------------------------------------------
    ! Trasformazione di Duschinsky
    !  Q_e = J Q_g + K
    !  dove Q_e sono le vibrazioni dello stato eccitato e Q_g quelle
    !  dello stato fondamenta..
   
    real(kind=dp), allocatable :: aus(:,:)
    integer  i, j, k, nmt

	!-------------------------------------------------------------------------
	! M e Q sono una super-matrice ed un super-vettore contenenti le matrici
	! necessarie per applicare le formule di ricorrenza in modo compatto.
	!-------------------------------------------------------------------------
   
   	! aus = GM*GM'
    allocate(aus(1:nmv,1:nmv))
    aus = zero
    do i = 1, nmv
      do j = i, nmv
        aus(i,j) = dot_product(GM(i,1:nmv),GM(j,1:nmv))
        if (i /= j) aus(j,i) = aus(i,j)
      end do
    end do
    !----------------------------------------------------
    ! Dall'articolo sul JCP (2008) 129, 64116:
    !
    ! M : equivale a -A dell'articolo
    ! A : C dell'articolo
    ! X^{-1} : -A^{-1}*tr(J)
    ! B^{-1}*J*J' = J*A*J'
    !
    ! M = [A^{-1}-I    tr(X^-1)       ]
    !     [X^-1        J*A^{-1}tr(J)-I]
    !----------------------------------------------------
    nmt = 2*nmv
    allocate (M(1:nmt,1:nmt),Q(1:nmt))
    M = zero; Q = zero
    M(1:nmv,1:nmv) = A(1:nmv,1:nmv)
    M(nmv+1:nmt,nmv+1:nmt) = matmul(B,aus)
    M(1:nmv,nmv+1:nmt) = transpose(X(1:nmv,1:nmv))
    M(nmv+1:nmt,1:nmv) = X(1:nmv,1:nmv)
        
    forall (i=1:nmt) M(i,i) = M(i,i) - one

    !-----------------------------------------
    ! Nell'articolo
    ! y = Q =[-tr(d)*X^{-1} tr(d)*B^{-1} ]*sqrt(2)/2
    ! Qui
    ! Q =[-tr(d)*X^{-1} tr(d)*B^{-1} ]/sqrt(2)
    !
    ! Dall'articolo sul JCP:
    !
    !  -tr(d)*C^{-1}*J/2 = -tr(d)*X^{-1}/2
    ! I-J*C^{-1}*tr(J)/2 = tr(d)*[(I+J*tr(J))]^{-1} = tr(d)*B^{-1}/2
    !-----------------------------------------
    allocate ( TR(1:nmv), TIR(1:nmv) )
    TIR = -matmul(d,X)
    TR = matmul(d,B)
    Q = (/TIR, TR /)
    Q = Q/sqrt(2.0_dp)

    deallocate(TR,TIR,aus)

    return
    end subroutine recursion_matrices  
     
	!==============================================================================!

    function MDFC_00()  result(WFC00)
    !------------------------------------+
    ! Calculate <0|0> FC integral as     !
    ! <0|0> = det(X)^{-1/2}*exp(-d'*B*d) !
    !------------------------------------+   
    !---------------------------------------------------------------------------+   
    ! Questa function usa la formula di Watson che mi sembra numericamente      |
    ! stabile in quanto non calcola il termine 2^N that is incluso nel                 |
    ! determinante.                                                             |
    !---------------------------------------------------------------------------+    
    real(kind=dp) WFC00
    real(kind=dp) detx, dbd
    integer i, j
   
    ! N.B.: Nel mio caso X contiene X^{-1}, e det(X^-1)=1/det(X)   
    detx = det(X)

    dbd = zero
    do i = 1, nmv
      do j = 1, nmv
        dbd = dbd + d(i)*d(j)*B(i,j) ! dall'articolo di Watson
      end do
    end do

    dbd = -0.25_dp*dbd
    WFC00 = exp(dbd)*sqrt(abs(detx))  

    print *, detx, dbd
    print *, WFC00
    return
    end function MDFC_00 

	!==============================================================================!

    subroutine fc_end 
    !-----------------------------------------------------------------------------!
    ! Deallocate all arrayes used in the calculation of FC integrals by recursion !
    ! formulae                                                                    !
    !-----------------------------------------------------------------------------!    
	integer :: alloc_err, ierr
	
    deallocate(M,Q,stat=alloc_err)
    if (alloc_err /= 0) ierr = error(0,"Cannot deallocate M and Q in fc_end")
    deallocate(f1,f2,stat=alloc_err) 
    if (alloc_err /= 0) ierr = error(0,"Cannot deallocate f1 and f2 in fc_end")
    deallocate(A,B,X,d,stat=alloc_err)
    if (alloc_err /= 0) ierr = error(0,"Cannot deallocate A,B,X,d  in fc_end")
    
    return
    end subroutine fc_end
    
end module fc

