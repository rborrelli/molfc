MODULE fcint_pert

	USE fcint_h
	USE fc_h, only : M, Q
	USE fc, only : d
	
	implicit none
	
	public :: fcint_pt
	
	private
	
	integer :: NFC
	
	CONTAINS
		
	!==========================================================================!
	
	SUBROUTINE fcint_pt (acs_initial,acs_final,FC00,fcjob,gid)
	
	!---------------------------------------------------------------+
    ! The matrix M  used is made accessible via  a "use" statement. | 
    !---------------------------------------------------------------+
	
  	type(activespace_t), intent(in) :: acs_initial,acs_final
  	real(kind=dp), intent(in) :: FC00
  	type(fc_t), intent(in) :: fcjob
    integer, intent(in) :: gid
  	
	! Local variables  	  
	integer, pointer :: v(:), w(:)
	integer, allocatable :: wo(:), wp(:)
	type(fci_t), allocatable :: fci(:)
	
    real(kind=dp) :: coeff, time1, time2, fcsum, fac, fct, pt, cfs
    real(kind=dp) :: PT1, PT2, PT3, PT1E, PTRS, PTRR, fnum, half, pr
    real(kind=dp) :: MFC0, MFC2, MPFC
    real(kind=dp) :: QQ(1:2), A0(1:2,1:2)
    real(kind=dp), allocatable :: MFC(:), XP(:,:)
    real(kind=dp), parameter :: duschthr = zero ! threshold value for Duschinsky effect
    real(kind=dp) :: efc, Ej, EM, Emin, Emax, Temp, fwhm, bw, dex
    
    integer, allocatable :: ivib(:), zmax(:), s(:), NINITIAL(:), NFINAL(:)
    integer, allocatable, target ::  z(:), z0(:)
    integer :: zi, zj, zk, zl, zp, zt
        		
    integer :: mu   ! Index of the excited vibration.
    integer :: nexvib

    integer :: ierr, alloc_err
    integer :: i, j, k, l, p, r, t, nmt, nmd, iss, kx, ii, jj, is, ix, rmu
    integer :: ntot, npt, i2, i3
    integer :: NFC1, si, sj, iv, jv, kv, lv, kk, NVIB_INCLUDE, NSZIV
    logical :: lc1, lc2
    
    real(kind=dp) :: enfc, ebra, eket
    real(kind=dp), allocatable :: freq_initial(:), freq_finale(:) 

    write(fout,'(2x,a,/)')'Perturbative Franck-Condon calculation'
       
    nexvib = acs_initial%nact + acs_final%nact
    !print *, 'NEXVIB ', nexvib
    if (acs_initial%nact > 0 .and. acs_final%nact > 0) NFC = product(acs_initial%nqmax+1)*product(acs_final%nqmax+1)
    if (acs_initial%nact == 0 .and. acs_final%nact > 0) NFC = product(acs_final%nqmax+1)
    if (acs_initial%nact > 0 .and. acs_final%nact == 0) NFC = product(acs_initial%nqmax+1)
    if (acs_initial%nact == 0 .and. acs_final%nact == 0) NFC = 1

    !allocate (MDFC(1:NFC),stat=alloc_err)
    !if (alloc_err /= 0) ierr = error(0,'Cannot allocate MDFC integrals array.')

    nmt = size(M,2)
    nmd = int(nmt/2)
    
	allocate (z(1:nmt), zmax(1:nmt), s(1:nmt), z0(1:nmt))
	allocate (fci(1:nmd))
	allocate (NINITIAL(1:nmd),NFINAL(1:nmd))
	allocate(XP(1:nmt,1:nmt))
	XP = zero
	! N.B.: Le vibrazioni definite nello spazio attivo (acs), come lette in input non sono in genere ordinate,
	! a meno che  l'utente non le abbia ordinate a mano con la sezione <proc>.
	! D'altro canto e' necessario che vi sia
	! corrispondenza tra l'ordine delle vibrazioni in M, Q ed s(i) e quello in ivib quindi
	! IVIB deve essere riordinato.  (Il sorting di ivib e' fondamentale altrimenti si dovrebbero definire
	! gli s(i) in altro modo.)
	! Nella nuova versione di set_active_space() le vibid sono gia' riordinate

	if (nexvib > 0) then
	    allocate (ivib(1:nexvib))
	    ivib = 0
	else
        write(fout,'(2x,a,/)')'No active space selected.'
        return ! no active space.
	end if
    !if (nexvib > 0) allocate (ivib(1:nexvib))
    !ivib = 0

    zmax = 0
    if (acs_initial%nact > 0) then
        ivib(1:acs_initial%nact) = acs_initial%vibid ! Contiene le ID delle vibrazioni dello spazio attivo.
        zmax(ivib(1:acs_initial%nact)) = acs_initial%nqmax
    end if

    if (acs_final%nact > 0) then
        ivib(acs_initial%nact+1:nexvib) = nmd + acs_final%vibid ! Contiene le ID delle vibrazioni dello spazio attivo.
        zmax(ivib(acs_initial%nact+1:nexvib)) = acs_final%nqmax
    end if

    s(1) = 1
	do i = 2, nmt
	  s(i) = product(zmax(1:i-1)+1)
	end do 

    allocate (freq_initial(1:fcjob%group(gid)%nvib), freq_finale(1:fcjob%group(gid)%nvib))
	freq_initial = system%state(INITIAL_STATE)%molecule%normodes%vibration(fcjob%group(gid)%incvib%id(:))%freq
	freq_finale = system%state(FINAL_STATE)%molecule%normodes%vibration(fcjob%group(gid)%incvib%id(:))%freq
	
	NVIB_INCLUDE = count(ivib <= nmd)
	NSZIV = size(ivib)

	NINITIAL = 0; NFINAL = 0
    NINITIAL = zmax(1:nmd)
    NFINAL = zmax(nmd+1:nmt)

    ! Calculate zeroth order integrals: for each mode calculate FC of an 
    ! effective displaced oscillator		
	do i = 1, nmd
		! Zero order: A0 = M0
		QQ = (/Q(i), Q(i+nmd)/)
		A0(1,1) = M(i,i); 	A0(1,2) = M(i,i+nmd)
		A0(2,1) = A0(1,2);	A0(2,2) = M(i+nmd,i+nmd)
		allocate(fci(i)%FC(0:NINITIAL(i),0:NFINAL(i)))		
		call oned_fc(NINITIAL(i),NFINAL(i),fci(i)%FC,A0,QQ)
	end do
	
    !-------------------------------------------------------------!
	! Write fci on file. This is used in the fcwd_Metz subroutine !
	!-------------------------------------------------------------! 
	open(unit=ffcpt,file="fc.pt",status="unknown")
	do i = 1, nmd
		write(ffcpt,'(3(i5,2x))') i, NINITIAL(i), NFINAL(i)
		do j = 0, NINITIAL(i)
			do k = 0, NFINAL(i)
				write(ffcpt,'(3(i5,2x),g0.14)')i, j, k, fci(i)%FC(j,k)
			end do	
		end do
	end do	
	
    ! If spectrum output is required...
	if (fcjob%fcspec) then
	    ! apre il file su cui scrivere lo spettro
        i2 = int(gid/10)
        i3 = gid-i2*10
        i2 = i2+ichar('0')
        i3 = i3+ichar('0')	
	    ! Open file for spectrum
        open(unit=iuspec,file=fcjob%spectrum%file(1:len_trim(fcjob%spectrum%file))//char(i2)//char(i3), &
             status="unknown")
	    ! -- Definisce i parametri per il calcolo dello spettro 
	    Emin = fcjob%spectrum%emin
	    Emax = fcjob%spectrum%emax
	    if (Emax < 0) Emax = Energy(freq_finale,zmax(nmd+1:nmt)) + system%state(FINAL_STATE)%energy
	    EM = fcjob%spectrum%dE
	    npt = int((Emax - Emin)/EM) + 1
	    Temp = fcjob%spectrum%Temp	              
    end if
  
	
	!---------------------------------+	
    ! Start perturbative calculation  |
    !---------------------------------+
   	z = 0
    z0 = z ! copy z to z0 and change v, w pointers
	v => z0(1:nmd)
	w => z0(nmd+1:2*nmd) 
	fcsum = zero
	ntot = 0
	
	if (fcjob%order >= 1) then

	do r = 1, NFC

		z0 = z
		!----------------------------+
		! Zero order Franck-Condon   |
		!----------------------------+
		MFC0=fci(1)%FC(v(1),w(1))
    	do j = 2, NMD
    		MFC0 = MFC0*fci(j)%FC(v(j),w(j))
    	end do
    	    	    	
		!----------------------------+
		! First order perturbation   |
		!----------------------------+						
		pt1 = zero
PRT1:	if (fcjob%order >= 1) then
			do i = 1, nexvib
				do j = i+1, nexvib
					z0 = z
					if (ivib(j) == (ivib(i) + nmd)) cycle ! only if A0 is the diagonal part of M
					!if (abs(M(ivib(i),ivib(j))) <= duschthr) cycle 
					if (z(ivib(i)) >= 1 .and. z(ivib(j)) >= 1) then
						if ((r-s(ivib(i))-s(ivib(j))) > 0) then
							coeff = sqrt(dble(z(ivib(i))*z(ivib(j))))
							z0(ivib(i)) = z(ivib(i)) - 1
							z0(ivib(j)) = z(ivib(j)) - 1
							MPFC = fci(1)%FC(v(1),w(1))
							do p = 2, nmd
								MPFC = MPFC*fci(p)%FC(v(p),w(p))
							end do
							!pt = pt + coeff*M(i,j)*MDFC(r-s(i)-s(j)) !...OK with arrays
							pt1 = pt1 + coeff*M(ivib(i),ivib(j))*MPFC !...OK
							!pt1 = pt1 + coeff*XP(ivib(i),ivib(j))*MPFC !...OK
						end if
					end if
				end do
			end do
        end if PRT1
        
        pt2 = zero
PRT2:   if (fcjob%order >= 2) then
		!----------------------------+
		! Second order perturbation  |
		!----------------------------+
			do i = 1, nexvib
				zi = z(ivib(i))
				do j = i+1 , nexvib
					zj = z(ivib(j))
					if (ivib(j) == (ivib(i) + nmd)) cycle 
					!if (abs(M(ivib(i),ivib(j))) <= duschthr) cycle
					do k = 1, nexvib
						do l = k+1 , nexvib
						    zk = z(ivib(k))
                            zl = z(ivib(l))                         						
							z0 = z ! reset z0
							if (ivib(l) == (ivib(k) + nmd)) cycle 
							!if (abs(M(ivib(k),ivib(l))) <= duschthr) cycle
							if ((r-s(ivib(i))-s(ivib(j))-s(ivib(k))-s(ivib(l))) > 0) then
								if (ivib(i) == ivib(k) .or. ivib(j) == ivib(k)) zk = zk -1
								if (ivib(i) == ivib(l) .or. ivib(j) == ivib(l)) zl = zl -1
								if (zi*zj*zk*zl > 0 ) then ! istruzione superfula, ma per ora la tengo...
									coeff = sqrt(dble(zi*zj*zk*zl))
									z0(ivib(i)) = zi - 1
									z0(ivib(j)) = zj - 1
									z0(ivib(k)) = zk - 1
									z0(ivib(l)) = zl - 1
									MPFC = fci(1)%FC(v(1),w(1))
									do p = 2, nmd
										MPFC = MPFC*fci(p)%FC(v(p),w(p))
									end do
									!pt = pt + coeff*M(i,j)*M(k,l)*MDFC(r-s(i)-s(j)-s(k)-s(l)) !...OK
									pt2 = pt2 + coeff*M(ivib(i),ivib(j))*M(ivib(k),ivib(l))*MPFC !...OK
									!pt2 = pt2 + coeff*XP(ivib(i),ivib(j))*XP(ivib(k),ivib(l))*MPFC
								end if
							end if
						end do
					end do
				end do
			end do
		end if PRT2 
        
		MFC2 = (MFC0 + pt1 + 0.50_dp*pt2)*FC00
	
		fcsum = fcsum + MFC2**2

        ! This is needed to calculate all energies of states in both spectrum and output.
        ! N.B.: The two output can have different threshold values for printing
        if (fcjob%printfc .and. fcjob%fcspec) then
            if (abs(mfc2) >= min(fcjob%ftol,fcjob%spectrum%tol)) then
              eket = Energy(freq_finale,z(nmd+1:nmt))+ system%state(FINAL_STATE)%energy
              ebra = Energy(freq_initial,z(1:nmd)) + system%state(INITIAL_STATE)%energy
              enfc = eket - ebra                 
            end if
        end if
        		
		if (fcjob%printfc) then
			if (abs(MFC2) >= fcjob%ftol) then
			  write(fout,'(2x,a3)',advance='no')'  <'
			  do k = 1, NVIB_INCLUDE
			      write(fout,'(i3)', advance='no') z(ivib(k))
			  end do
			  write(fout,'(a1)',advance='no')'|'
			  do k = 1+NVIB_INCLUDE, NSZIV
			      write(fout,'(i3)', advance='no') z(ivib(k))
			  end do
			  write(fout,'(a3)',advance='no')'> '
			  write(fout,'(2x,g14.6,4x,f10.2,4x)') MFC2, enfc
		    end if
		end if
		
		
		if (fcjob%fcspec) then		    
		    if (abs(MFC2) >= fcjob%spectrum%tol) then
              dex = ebra/(bkcm*Temp) 
              if (enfc >= emin .and. enfc <= emax ) then
                ntot = ntot + 1
                cfs = exp(-dex)
                if (fcjob%spectrum%type == "abs") cfs = enfc*exp(-dex)
                if (fcjob%spectrum%type == "ems") cfs = (enfc**3)*exp(-dex)                
                write(iuspec,'(2x,f10.2,g14.5)',advance='no') enfc, (MFC2**2)*cfs
                !write(iuspec,'(2x,f10.2,g14.5,g14.5,f10.6)',advance='no') enfc, enfc*(MFC2**2)*exp(-dex),(MFC2**2), exp(-dex)                
                ! writes transition assignment
                  write(iuspec,'(2x,a3)',advance='no')'  <'
                  do k = 1, NVIB_INCLUDE
                      write(iuspec,'(i3)', advance='no') z(ivib(k))
                  end do
                  write(iuspec,'(a1)',advance='no')'|'
                  do k = 1+NVIB_INCLUDE, NSZIV
                      write(iuspec,'(i3)', advance='no') z(ivib(k))
                  end do
                  write(iuspec,'(a3)')'> '                  
              end if
            end if         				
		end if
				
		! change z...as usual
    	do i = 1, nexvib
          mu = ivib(i)
          if (z(mu) < zmax(mu)) then
            z(mu) = z(mu) + 1
            if (i > 1) z(1:ivib(i-1)) = 0
            exit
          end if
      	end do
		
	end do
    end if

    ! write FC factors sum
    write(fout,300) fcsum

    if (fcjob%fcspec) then
	    ! Write spectrum info on output
	    write(fout,*) ' Simulated Spectrum written in file ', fcjob%spectrum%file(1:len_trim(fcjob%spectrum%file))
	    write(fout,'(2x,a,f10.2)') 'Emin: ', Emin
	    write(fout,'(2x,a,f10.2)') 'Emax: ', Emax
	    write(fout,'(2x,a,f10.2)') 'DE: ', EM
	    write(fout,'(2x,a,i7)') 'NPT: ', NPT
	    write(fout,'(2x,a,f10.3)') 'FWHM: ', fwhm
	    write(fout,'(2x,a,f10.3)') 'Temp: ', Temp
    end if
    
    ! scrive sul file il numero di FC selezionati
	write(iuspec,'(''#'',1x,i7)')ntot
	! chiude il file su cui ha scritto lo spettro
	close(iuspec)
    	
	deallocate (freq_initial, freq_finale)
    deallocate (z,zmax,ivib,s,z0)
    deallocate(fci)	
	
	include 'formats'
	
	return
	END SUBROUTINE fcint_pt

	!==========================================================================!
	
	SUBROUTINE oned_fc(NBRA, NKET, FCI, A0, QQ)
	
	!----------------------------------------------------------!
	! This subroutine calculates one dimensional FC integrals  !
	! using a recursive formula.                               !
	!----------------------------------------------------------!
		
	integer, intent(in) :: NBRA, NKET
	real(kind=dp), intent(in) :: QQ(1:2), A0(1:2,1:2)
	real(kind=dp), intent(out) :: FCI(0:NBRA,0:NKET)
	real(kind=dp), allocatable :: ODFC(:)
!	real(kind=dp), intent(out) :: ODFC(1:(NBRA+1)*(NKET+1))
			
	! Local variables
	integer :: i, j, k, r, mu, rmu, iss, kx, zk, NFC1, nexvib, nabra, naket
	integer :: z(1:2), zmax(1:2), s(1:2)
	integer, allocatable :: ivib(:)
	real(kind=dp) :: coeff
    	
	NFC1 = (NBRA+1)*(NKET+1)
	
	zmax(1) = NBRA; zmax(2) = NKET

    nexvib = count(zmax > 0)
    
    if (nexvib > 0) then 
        allocate(ivib(1:nexvib))
        ivib = 0
    else
        FCI(0,0) = one
        return
    end if
    
    nabra = 0
    if (zmax(1) > 0) then 
        nabra = 1
        ivib(1) = 1
    end if
    
    naket = 0
    if (zmax(2) > 0) then
        naket = 1 
        ivib(1+nabra:nexvib) = 2
    end if    
   
    s(1) = 1
    s(2) = zmax(1) + 1
	
	allocate (ODFC(1:NFC1))
	
	ODFC = zero
	
	ODFC(1) = one
	z = 0
	
	do r = 2, NFC1

      do i = 1, nexvib
          mu = ivib(i)
          if (z(mu) < zmax(mu)) then
            z(mu) = z(mu) + 1
            if (i > 1) z(1:ivib(i-1)) = 0
            exit
          end if
      end do

      rmu = r - s(mu)
    
      ODFC(r) = QQ(mu) * ODFC(rmu)
	    
      do k = 1, nexvib
          kx = ivib(k)
          if ( rmu <= s(kx) ) cycle
          iss = rmu - s(kx)
          zk = z(kx)
          if ( kx == mu ) zk = zk - 1
          coeff = sqrt(real(zk,kind(one)))*A0(mu,kx)
          ODFC(r) = ODFC(r) + coeff * ODFC(iss)
      end do
        
      ODFC(r) = ODFC(r)/sqrt(real(z(mu),kind(one)))

    end do	
        
    FCI = reshape(ODFC,(/NBRA+1,NKET+1/))

    deallocate(ODFC,ivib)
    
    return
    END SUBROUTINE oned_fc
          	
END MODULE fcint_pert
