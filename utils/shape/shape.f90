program gshape

    use param
    use fgt
    
    implicit none

    integer i, j, r, NFC, NPT
    real(kind=dp) emin, emax, dE, fwhm, Ej, efc, tol, des, scalef, shift
    real(kind=dp) q, t1, t2, fwhm2, gnorm
    real(kind=dp), allocatable :: en(:), MDFC(:), spec(:)

    read(*,*) emin, emax, dE, fwhm!, tol
    read(*,*) NFC

    allocate(en(1:NFC),MDFC(1:NFC))
    do i = 1, NFC
        read(*,*)en(i), MDFC(i)        
    end do

    npt = nint((Emax - Emin)/dE)! + 1
    allocate(spec(1:npt))
    spec = 0.D0

    ! shift spectrum position
    en = en - emin
    ! scale position
    scalef = 1.0_dp/((emax-emin)+epsilon(one)) 
    en = en*scalef ! scale spectrum in [0 1)
    des = dE*scalef   ! scale gaussian amplitude
    fwhm2 = (scalef*fwhm)**2
    !normalize spectrum
    q = sum(MDFC)
    MDFC = MDFC/q
    
    print *, '# size spec', npt
    print *, '# spec scale factor', scalef
    call cpu_time(t1)
    call gausst(MDFC,en,des,fwhm2,spec)
    call cpu_time(t2)
    print *, '# total time ', t2 -t1
    
    ! Write transformed spectrum	
	gnorm = (2.0_dp*sqrt(log(2.0_dp)))/(fwhm*sqrt(PI))
	! This normalization implies: sum(spec)*dE = 1
    !spec = spec*gnorm
	! This normalization implies sum(spec) = 1
    spec = spec*gnorm*dE
	do j = 1, npt
	   write(*,*) Emin+(j-1)*dE, spec(j)
	end do

end program gshape
  
