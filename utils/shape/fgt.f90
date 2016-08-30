module fgt

    use param
    use sort
        
	implicit none
	
	private
	
	real(kind=dp), allocatable :: w(:), s(:), g(:), x(:), sb(:)
	real(kind=dp), allocatable ::  tc(:,:), hc(:,:)
	integer, allocatable :: NB(:), nsi(:), nsf(:), nxi(:), nxf(:) 
	integer :: npts
	integer :: NBOXES, nterms, N, NF, MC, ML
	real(kind=dp) :: tol, boxdim, boxrad, delta, dsq
	! stores factorial from 0 to 12 to avoid recalculation
	real(kind=dp) :: fact(0:14)
	parameter (fact = (/1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, & 
	                    3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0 /))
	
	public :: gausst
	
	
	contains
      
	!==========================================================================!
	
    subroutine gausst(wg,sg,dx,delt,pot)
    !
	!     this subroutine computes the fast gauss transform, that is
	!     the sum of gaussians of the form
	!
	!                m
	!     pot(j) =  sum  w(i) * exp[(t(j) - s(i))^2/delta]
	!               i=1           
	!
	!     the source positions and targets must lie within
	!     the unit box [0,1].
    !
    real(kind=dp), intent(inout), target :: wg(:), sg(:)
    real(kind=dp), intent(in) :: dx, delt
    real(kind=dp), intent(out), allocatable :: pot(:) 
    real(kind=dp) :: h
    integer, allocatable :: ids(:)
    integer i, j, jlasts, jlastx, alloc_err
    logical, allocatable :: mask(:), mask2(:)
    
    ! w, s e dx sono accessibili all'intero modulo
    allocate(w(1:size(wg)),s(1:size(sg)),stat=alloc_err)
    if (alloc_err /= 0) then 
        print *, "Cannot allocate w and s"
        stop
    end if
    
    w = wg; s = sg;
    
    tol = zero
    tol= max(1.d-12,tol)
    tol= min(1.d-3,tol)

    ! npts è il numero di punti tra 0 e 1 (dx è l'intervallo, già scalato)
    npts = nint(1.0_dp/dx)
    ! alloca l'array che contiene la transformata finale
    allocate(g(1:npts),x(1:npts),stat=alloc_err)
    if (alloc_err /= 0) then 
        print *, "Cannot allocate g and x"
        stop
    end if
    ! Initialize g and x     
    g = zero
    forall (i = 1:npts) x(i) = (i-1)*dx
    
    ! Box dimension
    ! io passo delt=FWHM^2 = (2*sqrt(ln(2))*sigma)^2 = 4*ln(2)*delta (delta=2*sigma^2)
    ! quindi ricavo delta... 
    delta = delt/(4.0_dp*log(2.0_dp))
    dsq = 1.0_dp/sqrt(delta) 
    boxdim = sqrt(0.5d0*delta) ! equivale a fissare r = 1/2 nell'algoritmo
    boxrad = boxdim/2.0_dp
    
    ! Parametri dell'algoritmo
    NBOXES = nint(1.0_dp/boxdim)    
    nterms = 8 ! grado dell'espansione di hermite: 8 equivale a tol = 1e-6
    allocate(tc(0:nterms,1:NBOXES),stat=alloc_err) ! cumulative Taylor coefficients
    if (alloc_err /= 0) then 
        print *, "Cannot allocate tc (Tayloc coefficients)"
        stop
    end if
    allocate(hc(0:nterms,1:NBOXES),stat=alloc_err) ! cumulative Hermite coefficients
    if (alloc_err /= 0) then 
        print *, "Cannot allocate hc (Hermite coefficients)"
        stop
    end if
	! Initialize Taylor and Hermite coefficients arrays
    tc = zero; hc = zero
    ! 2N+1 è il numero di box che interagiscono con la box i-esima
    N = nint(sqrt(abs(2.0_dp*log(tol))))    
    NF = nterms + 2 ! in base ad NF e ML decido che algoritmo applicare.
    ML = nterms + 2
    ! MC: number of targets in each interval [sB(i)-rad, sB(i) + rad)    
    MC = int(npts/NBOXES)

    allocate(NB(1:NBOXES), nsi(1:NBOXES), nsf(1:NBOXES))
    allocate(nxi(1:NBOXES), nxf(1:NBOXES))
    allocate(sb(1:NBOXES)) ! stores the centers of the boxes
    ! Sort s  (and w accordingly) with a quick sort algorithm
    call qsort2(s,w)
    ! NB(i) = number of sources in i-th box B.
    ! defines NB(i): number of sources belonging to interval [sB(i)-rad, sB(i) + rad)
    ! N.B.: right open interval
    write(*,*)'# input data sorted'
    
    !allocate(mask(1:size(s)))
    !mask = .false.
    !allocate(mask2(1:size(x)))
    !mask2 = .false.

    ! Define centers of boxes   
    sb(1) = boxrad
    do i = 1, NBOXES-1
        sb(i+1) = sb(i) + boxdim
    end do
    
	NB = 0; nsi = 0; nsf = 0; nxi = 0; nxf = 0
	jlasts = 0; jlastx = 0
    ! questo ciclo è troppo lento...deve valutare troppi array logici: da migliorare!!
    do i = 1, NBOXES
		do j = jlasts+1, size(s)
			if (s(j) >= (sb(i) + boxrad)) exit ! gli elementi di s sono ordinati è inutile proseguire!
			if (s(j) >= (sb(i) - boxrad)) then
				if (nsi(i) == 0) nsi(i) = j
				if (s(j) < (sb(i+1) + boxrad)) then
					NB(i) = NB(i)+1
					nsf(i) = j
					jlasts = j
				end if
			end if
		end do
        !mask = (s >= (sb(i+1) - boxrad)) .and. (s < (sb(i+1) + boxrad))
        !nsi(i+1:i+1) = minloc(s,mask) ! posizione della prima sorgente nella box i
        !nsf(i+1:i+1) = maxloc(s,mask) ! posizione dell'ultima sorgente nella box i
		!end do
		do j = jlastx+1, size(x)
			if (x(j) >= (sb(i) + boxrad)) exit ! gli elementi di s sono ordinati è inutile proseguire!
			if (x(j) >= (sb(i) - boxrad)) then
				if(nxi(i) == 0) nxi(i) = j
				if (x(j) < (sb(i) + boxrad)) then
					nxf(i) = j
					jlastx = j
				end if
			end if
		end do
		!mask2 = (x >= (sb(i+1) - boxrad)) .and. (x < (sb(i+1) + boxrad))
        !nxi(i+1:i+1) = minloc(x,mask2) ! primo punto della cella i+1
        !nxf(i+1:i+1) = maxloc(x,mask2) ! ultimo punto della cella i+1
    end do

    ! Fa in modo che l'ultima box includa tutti i punti...non è proprio ottimale ma è semplice
    !nsf(NBOXES) = size(s)
    !nxf(NBOXES) = npts
        
    write(flog,*) 'delta ', delta
    write(flog,*) 'boxdim ', boxdim
    write(flog,*) 'boxrad ', boxrad
    write(flog,*) 'NBOXES ', nboxes
    write(flog,*) 'interacting boxes',2*N+1
    write(flog,*) 'NF, ML, MC', NF, ML, MC
    write(flog,*) 'Sum NB', sum(NB) ! should be equal to the number of FC
    !write(flog,*) 'indici sorgeti nsi, nsf'
    !do i = 1, nboxes
    !    write(flog,*) i, sb(i), NB(i), nsi(i), nsf(i)
    !end do
    !write(flog,*) 'indici target nxi, nxf'
    !do i = 1, nboxes
    !    write(flog,*) nxi(i), nxf(i)
    !end do
    !-----------------------------------------------------
    ! gafexp: create all expansion coefficients on grid,
    ! evaluate all appropriate far field expansions,
    ! evaluate all appropriate direct interactions.
    ! then calls tlreval to evaluate cumulated Taylor expansions
    !-----------------------------------------------------    
    call gafexp()
    
    allocate(pot(1:size(g)))
    pot = g
    
    !deallocate(w,s,g)
    !deallocate(nb,nsi,nsf)
    
    return
    end subroutine gausst

    !==========================================================================!
    
    subroutine gafexp()

    integer i, j, jboxi, jboxf, idir, itlr, ihrd, ihrtr
	
    idir = 0; itlr = 0; ihrd = 0; ihrtr = 0 ! set four counters to zero
    ! Main subroutine of the Fast Gauss Transform Algorithm    	
	do i = 1, NBOXES
	    if (NB(i) == 0) cycle ! se non ci sono sorgenti nella cella passa alla successiva	    
	    !Form the interaction list of (2n + 1) target boxes C within range of B.
	    ! Nel nostro caso le box che interaciscono con la i-esima sono N a sinistra e N a destra.
	    ! Se N a sinistra è negativo allora si una come inizio la prima box.
        jboxi = max(1,i-N)
        jboxf = min(i+N,NBOXES)
		if (NB(i) <= NF) then
            ! The Gaussians of sources in box I can be evaluated in box J either directly
            ! or by Taylor expansion about the center of box J.		      
			do j = jboxi, jboxf  ! ciclo sulle target box nel range di interazione con la source box
				if (MC <= ML) then
				    !Compute source/target interactions by direct evaluation of Gaussians.
				    ! Valuta direttamente la somma delle Gaussiane della box i nei punti
				    ! della box j
				    idir = idir + 1
				    call gdirect(i,j)
				else
                    !Convert each of the NB sources into a Taylor series about the center of
                    !box C via (20) and add to Taylor series for box C.
                    ! Calcola solo i coefficienti degli sviluppi in serie delle varie sorgenti in B 
                    ! che contribuiscono al campo nella box C. Dopo aver accumulato tutti i coefficienti li valuta
                    ! chiamando tlreval()
				    itlr = itlr + 1
				    call tlrexp(i,j)
				end if
			end do
		else !(NB(i) > NF)
			!Form Hermite expansion about center of box B due to NB sources via equation (20).
			call hrmexp(i)
            ! The Hermite expansion can be evaluated in the 2N+1 interacting boxes in two different ways.
			do j = jboxi, jboxf
				if (MC <= ML) then
					!Evaluate Hermite expansion at each target location and add to accumulated
					!potential.
					ihrd = ihrd + 1
					call hrmeval(i,j)
				else
				    ihrtr = ihrtr + 1
!					!Convert Hermite expansion into a Taylor series about the center of box C
!					!by means of (12) and add to Taylor series for box C.
					call hrm2tlr(i,j)	
				end if
			end do
		end if
	end do
	
	! Evaluate taylor expansions (we have stored all the coefficients in the array tc(i,j)
	call tlreval()

    print *, '# Gaussians - Direct ', idir
    print *, '# Gaussian - Taylor ', itlr
    print *, '# Hermite series - Direct ', ihrd
    print *, '# Hermite series - Taylor ', ihrtr
    		
	return
    end subroutine gafexp

    !==========================================================================!
    
    subroutine gdirect(ibox,jbox)
    !
    ! Direct evaluation of Gaussians 
    ! Valuta le Gaussiane che si trovano nella box IBOX nei punti della box JBOX
    !
    integer, intent(in) :: ibox, jbox
    integer :: i, j
 
    do j = nxi(jbox), nxf(jbox)        
	    do i = nsi(ibox), nsf(ibox) ! bsi ed nfs sono gli array determinati sopra
            g(j) = g(j) + w(i)*exp(-(x(j) - s(i))**2/delta)
	    end do   
    end do
                    
    return
    end subroutine gdirect
    
    !==========================================================================!
    
    subroutine tlrexp(ibox,jbox)
    !
    ! Valuta i coefficienti della serie di Taylor del campo gaussian dovuto alla box IBOX
    ! rispetto al centro della box JBOX.
    ! I coefficienti per ogni IBOX vengono poi sommati nell'array TC(i,j)
    ! Dopo che l'intero array è stato calcolato si chiama gameval per il calcolo.
    ! delle serie.
    !
    integer, intent(in) :: ibox, jbox
    integer i, j, k
    real(kind=dp), allocatable :: b(:), hx(:)
    real(kind=dp) :: xd, fac, sig
    
    allocate(b(0:nterms), hx(0:nterms))
    
    do k = nsi(ibox), nsf(ibox)

        b = zero; hx = zero
        ! Store Hermite polynomials of variable -(s(k) - centro(jbox))/sqrt(delta)     
		! In tutte le versioni precedenti c'è un errore nel segno di xd.
        xd = -(s(k) - sb(jbox))*dsq
	    hx(0) = exp(-(xd**2))*w(k)
	    hx(1) = 2.0_dp*xd*hx(0)
	    do i = 1, nterms - 1
	        hx(i+1) = 2.0_dp*xd*hx(i) - 2.0_dp*dble(i)*hx(i-1)   
	    end do
	  
		sig = -1.0_dp
	    do j = 0, nterms
			sig = -1.0_dp*sig
	        b(j) = hx(j)
	        fac = sig/dble(fact(j))
            b(j) = fac*b(j)	        
	    end do	    

        ! sum Taylor series coefficients        
        tc(0:nterms,jbox) = tc(0:nterms,jbox) + b(0:nterms)
        
    end do
    
    deallocate(b,hx)
    
    return
    end subroutine tlrexp

    !==========================================================================!
    
    subroutine tlreval()
    !
    ! Valuta le serie di taylor nella box jbox con coefficienti t(i,jbox), i = 0,...,nterms
    !
    real(kind=dp), allocatable :: tys(:) ! temporary taylor sum
    real(kind=dp) xp
    integer i, j, k
    
    allocate(tys(1:npts))
    tys = zero
        
    do i = 1, NBOXES
        if (MC > ML) then
	        do j = nxi(i), nxf(i)
	            xp = (x(j) - sb(i))*dsq
	            tys(j) = zero	           
                do k = nterms, 1, -1
                    tys(j) = (tys(j) + tc(k,i))*xp
                end do
	        end do
	        tys(nxi(i):nxf(i)) = tys(nxi(i):nxf(i)) + tc(0,i) ! add costant term to i-ith series 
        end if
    end do

    ! Add taylor series to final transformation g.
    g = g + tys
    
    deallocate(tys)
    
    return
    end subroutine tlreval
    
    !==========================================================================!
    
    subroutine hrmexp(ibox)
    
    ! Create Hermite expansion of the Gaussian field of source box IBOX
    ! Cioè calcola i coefficienti dell'espansione di Hermite che vengono
    ! memorizzati nell'array hc(i,j)
    ! Uso la relazione A(k) = (1/k!)*sum_i w(i)*d(i)^k con d = (s(i) - centro(ibox))/sqrt(delta)
    ! ossia:
    ! A(0) = sum_i w(i)
    ! A(1) = sum_i w(i)*d(i)
    ! A(2) = (1/2!)sum_i w(i)*d(i)^2
    ! ...
    integer, intent(in) :: ibox
    integer :: i, j, k
    real(kind=dp) :: xd
    
    ! E' più efficiente se prendo una sorgente e calcolo i contributi di questa a tutti i coefficienti
    do i = nsi(ibox), nsf(ibox)
    ! Calculate term A(0) it is independent from the sources s(i)
        hc(0,ibox) = hc(0,ibox) + w(i)
        xd = (s(i) - sb(ibox))*dsq
        do k = 1, nterms
            hc(k,ibox) = hc(k,ibox) + w(i)*xd
            xd = xd*xd
        end do        
    end do
    
    ! divide by k! (N.B.: fact(0:1) = 1)    
    do k = 2, nterms
        hc(k,ibox) = hc(k,ibox)/fact(k)
    end do
    
    return
    end subroutine hrmexp
    
    !==========================================================================!
    
    subroutine hrmeval(ibox,jbox)
    
    ! Direct evaluation of Hermite expansion of box IBOX in the targets of box JBOX
    ! add the result to the final  g array.
    integer, intent(in) :: ibox, jbox
    integer:: i, j, k
    real(kind=dp) :: xd
    real(kind=dp), allocatable :: hx(:) 
    
    allocate(hx(0:nterms))
    hx = zero
    
    ! For each target in box JBOX...
    do k = nxi(jbox), nxf(jbox)

        ! Store Hermite polynomials of variable (x(k) - sb(ibox))/sqrt(delta)         
	    xd = (x(k) - sb(ibox))*dsq
	    hx(0) = exp(-xd**2)
	    hx(1) = 2.0_dp*xd*hx(0)
	    do i = 1, nterms - 1
	        hx(i+1) = 2.0_dp*xd*hx(i) - 2.0_dp*i*hx(i-1)   
	    end do
     
        do j = 0, nterms   
            g(k) = g(k) + hc(j,ibox)*hx(j)
        end do
         
    end do
    
    deallocate(hx)
    
    return
    end subroutine hrmeval
                  
    !==========================================================================!
    
    subroutine hrm2tlr(ibox,jbox)
    
    ! Transform the Hermite expansion of the Gaussian field of box IBOX 
    ! to a Taylor expansion abot the center of box JBOX; calculate the coefficients of the expansions
    ! and add to the coefficents of other Taylor expansion in JBOX 
    
    integer, intent(in) :: ibox, jbox
    integer :: i, j, k
    real(kind=dp) :: xp, fac, lxp
    real(kind=dp), allocatable :: hx(:)
    
    allocate(hx(0:2*nterms))
    hx = zero
    !xp = (sb(ibox) - sb(jbox))*dsq ! controllare che non sia l'inverso...cambiarebbero segno i polinomi!
    xp = (sb(jbox) - sb(ibox))*dsq ! le formule sull'articolo riportano l'opposto di questo valore
     
    ! Store Hermite polynomials
    hx(0) = exp(-(xp**2))
    hx(1) = 2.0_dp*xp*hx(0)
    do i = 1, 2*nterms - 1
        hx(i+1) = 2.0_dp*xp*hx(i) - 2.0_dp*i*hx(i-1)     
    end do
    
    ! Add expansion to previously calculated Taylor coefficients tc(i,j) (see tlrexp)
    do i = 0, nterms
        fac = ((-1)**i)/fact(i)
        lxp = zero
        do k = 0, nterms
            lxp = lxp + hc(k,ibox)*hx(i+k)
        end do
        tc(i,jbox) = tc(i,jbox) + fac*lxp
    end do
     
    deallocate(hx)
        
    return
    end subroutine hrm2tlr
                      
end module fgt
