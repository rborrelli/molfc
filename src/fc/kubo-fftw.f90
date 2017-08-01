!     
! File:   kubo.f90
! Author: lello
!
! Created on November 18, 2010, 5:33 PM
!

module kubo

    use parameters
    use hypfun
    use matfun
    use transf_type
    use fc_type
    use system_type
    use omp_lib
    use, intrinsic :: ieee_arithmetic

    implicit none

    private

    public :: kubo_lineshape

    real(kind=dp), allocatable :: x(:)
    complex(kind=dpc), allocatable :: f(:)

contains

 	!===========================================================================!

    subroutine genfun(D,S,wg,we,T,ntau,OmR,fwhm)

        ! wg: frequencies of the INITIAL STATE (groud)
        ! we: frequencies of the FINAL STATE (excited)
        ! I/O variables
        real(kind=dp), intent(in) :: D(:), S(:,:), wg(:), we(:), T, OmR, fwhm
        integer, intent(in) :: ntau

        ! Local variables
        complex(kind=dpc) :: z
        real(kind=dp) :: beta, tau, dtau, OmegaRange, domega
        real(kind=dp) :: TauRange, si, at, atm, atc
        real(kind=dp) :: wcot, wtan, StS
        real(kind=dp) :: t1, t2, xi(1:4), yi(1:4), fk_r, fk_i
        
        real(kind=dp), allocatable, dimension(:) :: tan2t, st2, pf!, thz1
        real(kind=dp), allocatable :: SR(:)
        complex(kind=dpc), allocatable, dimension(:) :: Tg, Cg, DC
        complex(kind=dpc), allocatable, dimension(:,:) :: Xm, Te, TH
        complex(kind=dpc), allocatable, dimension(:,:) :: Ym, Ce, TT
        complex(kind=dpc), allocatable, dimension(:) :: Y, phi, thz1
        complex(kind=dpc), allocatable, dimension(:) :: G, sdphi, detphi
        integer :: i, k, m, p, q, N, nvib

        ! The dimension of the problem (the number of vibrations)
        nvib = size(D)
        ! beta Temperature
        beta = 1.0_dp/(bkcm*T)
        ! Length of the Fourier transform N = 2^(p+1)
        N = 2*ntau;
        OmegaRange = 0.5_dp*OmR;
        domega = OmR/real(N,dp)
        dtau = TWOPI/OmR
        TauRange = ntau*dtau

        ! Write some information on a log file
        open(flog,file='dft.log',status='unknown')
        write(flog,*) 'N ', N
        write(flog,*) 'Omega range ', "[",-OmegaRange,OmegaRange,"]"
        write(flog,*) 'domega ', domega
        write(flog,*) 'dtau ', dtau
        write(flog,*) 'Tau range ', TauRange

        ! Write FFT information in the output file
        write(fout,*) '====================================== '
        write(fout,*) 'Discrete Fourier Transform parameters: '
        write(fout,*) 'N ', N
        write(fout,*) 'Omega range ', "[",-OmegaRange,OmegaRange,"]"
        write(fout,*) 'domega ', domega
        write(fout,*) 'dtau ', dtau
        write(fout,*) 'Tau range ', TauRange
        write(fout,*) '====================================== '

        allocate(SR(1:nvib,1:nvib))
        SR = invm(S)
        allocate(f(1:N))
        allocate(G(1:ntau),sdphi(1:ntau),detphi(1:ntau))
        allocate(pf(1:nvib))
        allocate(DC(1:nvib))

        DC = cmplx(D,zero,dpc)

        pf = (one-exp(-beta*wg))**2
        call cpu_time(t1)
        !$OMP PARALLEL  DEFAULT(none) PRIVATE(k,tau,z,p,q,m,StS,wtan,wcot,at,atm,atc,tan2t,thz1,Xm,Te,Tg,Ym,Ce,Cg,TH,TT,st2,phi,Y) SHARED(sdphi,G,S,SR,D,DC,dtau,wg,we,nvib,beta,pf,ntau,detphi)
        allocate(Xm(1:nvib,1:nvib),Tg(1:nvib),Te(1:nvib,1:nvib))
        allocate(Ym(1:nvib,1:nvib),Cg(1:nvib),Ce(1:nvib,1:nvib))
        allocate(TH(1:nvib,1:nvib),TT(1:nvib,1:nvib))
        allocate(Y(1:nvib),phi(1:nvib))
        allocate(tan2t(1:nvib), st2(1:nvib), thz1(1:nvib))
        !$OMP DO
        do k = 1, ntau-1
            tau = (real(k,dp))*dtau
            z = cmplx(beta,-tau,dpc)

            tan2t = tan((tau/2.0_dp)*we)
            thz1 = tanh((z/2.0_dp)*wg)
            Tg = wg*thz1
            Cg = wg/thz1

            do p = 1, nvib
                do q = p, nvib
                    at = zero; atm = zero; atc = zero
                    do m = 1, nvib
                        StS = S(m,p)*S(m,q)
                        wtan = we(m)*tan2t(m)
                        wcot = we(m)/tan2t(m)
                        at = at + StS*wtan
                        atc = atc + StS*wcot
                        atm = atm + SR(p,m)*SR(q,m)/wtan
                    end do
                    Te(p,q) = cmplx(zero,at,dpc)
                    Ce(p,q) = cmplx(zero,-atc,dpc)
                    TT(p,q) = cmplx(zero,-atm,dpc)
                    if (p /= q) then
                        Te(q,p) = Te(p,q)
                        Ce(q,p) = Ce(p,q)
                        TT(q,p) = TT(p,q)
                    end if
                end do
            end do
            ! Xm = Te + Tg
            Xm = Te
            forall(p=1:nvib)
                Xm(p,p) = Xm(p,p) + Tg(p)
            end forall
            ! Ym = Ce+Cg
            Ym = Ce
            forall(p=1:nvib)
                Ym(p,p) = Ym(p,p) + Cg(p)
            end forall
            ! TT = Te^{-1} + Tg^{-1}
            forall(p=1:nvib)
                TT(p,p) = TT(p,p) + one/Tg(p)
            end forall
            ! Start calculation of Phi and its determinant
            TH = matmul(Xm,Ym)
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)/(wg(p)*we(q))
            end forall
            st2 = sin(tau*we)
            phi = 0.5_dp*cmplx(st2*sin(tau*wg),st2*cos(tau*wg),dpc)* &
            (one-exp(-2.0_dp*z*wg))/pf
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)*phi(q)
            end forall
            detphi(k) = det(TH)
            sdphi(k) = 1.0_dp/sqrt(detphi(k))

            Y = zinvmv(TT,DC)
            G(k) = cfac*dot_product(D,Y)
        end do
        !$OMP END DO
        deallocate(Xm,Tg,Te,Ym,Cg,Ce)
        deallocate(TH,TT,Y,phi)
        deallocate(tan2t, st2, thz1)
        !$OMP END PARALLEL
        call cpu_time(t2)

        write(flog,*) 'Time ', t2-t1
        !------------------------------------------------------------------!
        ! Check the sign of the square root: this is fundamental           !
        ! we must always follow the same sheet of the complex square root  !
        !------------------------------------------------------------------!
        si = one
        do k = 2, ntau-1
            if (real(detphi(k)) < 0 .and. real(detphi(k-1)) < 0) then
                if (aimag(detphi(k-1))*aimag(detphi(k)) < 0) si = -si
            end if
            sdphi(k) = si*sdphi(k) ! Set the sign of the square root
			!NEW version 23/7/2014
			!print *, tau, real(detphi(k)), imag(detphi(k))
			!sdphi(k) = cmplx(real(sdphi(k),dp),si*imag(sdphi(k)),dpc)	
			!print *, 'ATAN2 ', k, atan2(imag(detphi(k)),real(detphi(k)))
        end do

        !------------------------------------------------------------------!
        ! Compute Generating Function and apply a window (apodiziation)
        !------------------------------------------------------------------!
        f(1) = cmplx(one,zero,dpc)*BH4(1,N)
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k,tau) SHARED(f,sdphi,G,ntau,N,dtau,detphi)
        do k = 1, ntau-1
            tau = (real(k,dp))*dtau
            !do k = 1, N
            ! Use Blackman-harris window
            f(k+1) = sdphi(k)*exp(-G(k))*BH4(k+1,N)
			!print *, 'DETPHI ', k, detphi(k), sdphi(k)
	        ! Use Kaiser-Bessel window
	        !f(ntau+k) = sdphi(k)*exp(-G(k))*KBESS(k,N)
	        !f(k+1) = sdphi(k)*exp(-G(k))*KBESS(k,N)
	        ! No window
	        !f(ntau+k) = sdphi(k)*exp(-G(k))
	        !f(k+1) = sdphi(k)*exp(-G(k))
        end do
        !$OMP END PARALLEL DO

        ! This is a workaround to remove singularities. It works only if a NaN exception
        ! is catched by the compiler.
        xi(1:4) = (/ 0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp/)
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k,yi,fk_r,fk_i) SHARED(f,ntau,xi)
        do k = 3, ntau-2
            if (ieee_is_nan(real(f(k))) .or. ieee_is_nan(aimag(f(k)))) then
                yi(1:4) = (/real(f(k-2)), real(f(k-1)), real(f(k+1)), real(f(k+2))/)
                fk_r = cubic_interp(2.0_dp, xi, yi)
                yi(1:4) = (/aimag(f(k-2)), imag(f(k-1)), imag(f(k+1)), imag(f(k+2))/)
                fk_i = cubic_interp(2.0_dp, xi, yi)
                f(k) = cmplx(fk_r, fk_i,dpc)
            end if
        end do
        !$OMP END PARALLEL DO

        if (fwhm > 0.0) then
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k) SHARED(f,ntau,N,fwhm)
            do k = 1, ntau
                ! Add Lorentzian broadening
                f(k) = f(k)*expl(k,N,fwhm)
            end do
        !$OMP END PARALLEL DO
        end if

        deallocate(G,sdphi,detphi)

        ! Use complex conjugate symmetry
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(f,ntau,N)
        do i= 1, ntau
            f(modulo(N-i+1,N)+1) = conjg(f(i))
        end do
        !$OMP END PARALLEL DO

        allocate(x(1:N))

        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(x,OmegaRange,domega,N)
        do i = 1, N
            x(i) = -OmegaRange +real(i,dp)*domega
        end do
        !$OMP END PARALLEL DO

        ! Shift the origin
        x = x + (sum(wg-we))/2.0_dp


        write(97,*)'origin shift', sum(wg-we)/2.0_dp
        return
    end subroutine genfun

 	    !===========================================================================!

    subroutine genfun_zero(D,S,wg,we,ntau,OmR,fwhm)

        ! I/O varibles
        real(kind=dp), intent(in) :: D(:), S(:,:), wg(:), we(:), OmR, fwhm
        integer, intent(in) :: ntau

        ! Local variables
        
        real(kind=dp) :: tau, dtau, OmegaRange, domega
        real(kind=dp) :: TauRange, si, at, atm, atc
        real(kind=dp) :: wcot, wtan, StS
        real(kind=dp) :: t1, t2, xi(1:4), yi(1:4), fk_r, fk_i
        real(kind=dp), allocatable, dimension(:) :: tan2t, st2
        complex(kind=dpc), allocatable, dimension(:) :: Tg, Cg, DC
        complex(kind=dpc), allocatable, dimension(:,:) :: TH1, Te, TH
        complex(kind=dpc), allocatable, dimension(:,:) :: TH2, Ce, TT
        complex(kind=dpc), allocatable, dimension(:) :: Y, phi
        complex(kind=dpc), allocatable, dimension(:) :: G, sdphi, detphi
        integer :: i, k, m, p, q, N, nvib

        ! The dimension of the problem (the number of vibrations)
        nvib = size(D)
        ! Length of the Fourier transform N = 2^(p+1)
        N = 2*ntau;
        OmegaRange = 0.5_dp*OmR;
        domega = OmR/real(N,dp)
        dtau = TWOPI/OmR
        TauRange = ntau*dtau

        ! Write some information on a log file
        open(flog,file='dft.log',status='unknown')
        write(flog,*) 'N ', N
        write(flog,*) 'Omega range ', OmegaRange
        write(flog,*) 'domega ', domega
        write(flog,*) 'dtau ', dtau
        write(flog,*) 'Tau range ', TauRange

        allocate(f(1:N))
        allocate(G(1:ntau),sdphi(1:ntau),detphi(1:ntau))
        allocate(DC(1:nvib))

        DC = cmplx(D,zero,dpc)

        allocate(Tg(1:nvib),Cg(1:nvib))
        Tg = wg
        Cg = wg

        call cpu_time(t1)
        !$OMP PARALLEL  DEFAULT(none) PRIVATE(i,k,tau,p,q,m,StS,wtan,wcot,at,atm,atc,tan2t,TH1,Te,TH2,Ce,TH,TT,st2,phi,Y) SHARED(sdphi,Tg,Cg,G,S,D,DC,dtau,wg,we,nvib,ntau,detphi)
        allocate(TH1(1:nvib,1:nvib),Te(1:nvib,1:nvib))
        allocate(TH2(1:nvib,1:nvib),Ce(1:nvib,1:nvib))
        allocate(TH(1:nvib,1:nvib),TT(1:nvib,1:nvib))
        allocate(Y(1:nvib),phi(1:nvib))
        allocate(tan2t(1:nvib), st2(1:nvib))
        !$OMP DO
        do k = 1, ntau-1
            tau = (real(k,dp))*dtau

            tan2t = tan((tau/2.0_dp)*we)

            do p = 1, nvib
                do q = p, nvib
                    at = zero; atm = zero; atc = zero
                    do m = 1, nvib
                        StS = S(m,p)*S(m,q)
                        wtan = we(m)*tan2t(m)
                        wcot = we(m)/tan2t(m)
                        at = at + StS*wtan
                        atc = atc + StS*wcot
                        atm = atm + StS/wtan
                    end do
                    !write(777,*) k, p, q, at, atc, atm
                    Te(p,q) = cmplx(zero,at,dpc)
                    Ce(p,q) = cmplx(zero,-atc,dpc)
                    TT(p,q) = cmplx(zero,-atm,dpc)
                    if (p /= q) then
                        Te(q,p) = Te(p,q)
                        Ce(q,p) = Ce(p,q)
                        TT(q,p) = TT(p,q)
                    end if
                end do
            end do
            TH1 = Te
            forall(p=1:nvib)
                TH1(p,p) = Te(p,p) + Tg(p)
            end forall
            TH2 = Ce
            forall(p=1:nvib)
                TH2(p,p) = TH2(p,p) + Cg(p)
            end forall
            forall(p=1:nvib)
                TT(p,p) = TT(p,p) + one/Tg(p)
            end forall
            ! Start calculation of Phi and its determinant
            TH = matmul(TH1,TH2)
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)/(wg(p)*we(q))
            end forall
            st2 = sin(tau*we)
            phi = 0.5_dp*cmplx(st2*sin(tau*wg),st2*cos(tau*wg),dpc)
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)*phi(q)
            end forall
            detphi(k) = det(TH)
            sdphi(k) = 1.0_dp/sqrt(detphi(k))

            ! compute the exponent part
            Y = zinvmv(TT,DC)
            G(k) = cfac*dot_product(D,Y)
        end do
        !$OMP END DO
        deallocate(TH1,Te,TH2,Ce)
        deallocate(TH,TT,Y,phi)
        deallocate(tan2t, st2)
        !$OMP END PARALLEL
        call cpu_time(t2)

        write(flog,*) 'Time ', t2-t1
        !------------------------------------------------------------------!
        ! Check the sign of the square root: this is fundamental           !
        ! we must always follow the same sheet of the complex square root  !
        !------------------------------------------------------------------!
        si = one
        do k = 2, ntau-1
            if (real(detphi(k)) < 0 .and. real(detphi(k-1)) < 0) then
                if (aimag(detphi(k-1))*aimag(detphi(k)) < 0) si = -si
            end if
            sdphi(k) = si*sdphi(k) ! Set the sign of the square root
        end do

        !------------------------------------------------------------------!
        ! Compute Generating Function and apply a window (apodiziation)
        !------------------------------------------------------------------!
        f(1) = cmplx(one,zero,dpc)
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k) SHARED(f,sdphi,G,ntau,N)
        do k = 1, ntau-1
            !do k = 1, N
            ! Use Lorentzian window
            !f(ntau+k) = sdphi(k)*exp(-G(k))*EXPL(k,N)
            ! Use Blackman-harris window
            !f(ntau+k) = sdphi(k)*exp(-G(k))*BH4(k,N)
            f(k+1) = sdphi(k)*exp(-G(k))*BH4(k+1,N)
	        ! Use Kaiser-Bessel window
	        !f(ntau+k) = sdphi(k)*exp(-G(k))*KBESS(k,N)
	        ! No window
	        !f(ntau+k) = sdphi(k)*exp(-G(k))
	        !f(k+1) = sdphi(k)*exp(-G(k))
        end do
        !$OMP END PARALLEL DO

        ! This is a workaround to remove singularities. It works only if a NaN exception
        ! is catched by the compiler.
        xi(1:4) = (/ 0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp/)
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k,yi,fk_r,fk_i) SHARED(f,ntau,xi)
        do k = 3, ntau-2
            if (ieee_is_nan(real(f(k))) .or. ieee_is_nan(aimag(f(k)))) then
                yi(1:4) = (/real(f(k-2)), real(f(k-1)), real(f(k+1)), real(f(k+2))/)
                fk_r = cubic_interp(2.0_dp, xi, yi)
                yi(1:4) = (/aimag(f(k-2)), imag(f(k-1)), imag(f(k+1)), imag(f(k+2))/)
                fk_i = cubic_interp(2.0_dp, xi, yi)
                f(k) = cmplx(fk_r, fk_i,dpc)
            end if
        end do
        !$OMP END PARALLEL DO

        if (fwhm > 0.0) then
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k) SHARED(f,ntau,N,fwhm)
            do k = 1, ntau
                ! Add Lorentzian broadening
                f(k) = f(k)*expl(k,N,fwhm)
            end do
        !$OMP END PARALLEL DO
        end if

        deallocate(G,sdphi,detphi)

        ! Use complex conjugate symmetry
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(f,ntau,N)
        do i= 1, ntau
            f(modulo(N-i+1,N)+1) = conjg(f(i))
        end do
        !$OMP END PARALLEL DO

        allocate(x(1:N))

        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(x,OmegaRange,domega,N)
        do i = 1, N
            x(i) = -OmegaRange +(real(i,dp))*domega
        end do
        !$OMP END PARALLEL DO

        ! Shift the origin
        x = x + (sum(wg-we))/2.0_dp

        write(97,*)'origin shift', sum(wg-we)/2.0_dp
        return
    end subroutine genfun_zero

 	    !===========================================================================!

    subroutine genfun_ht(D,S,wg,we,T,ntau,OmR,fwhm,tm)

        ! I/O variables
        real(kind=dp), intent(in) :: D(:), S(:,:), wg(:), we(:), T, OmR, fwhm
        type(tm_t), intent(in), optional :: tm
        integer, intent(in) :: ntau

        ! Local variables
        complex(kind=dpc) :: z
        real(kind=dp) :: beta, tau, dtau, OmegaRange, domega
        real(kind=dp) :: TauRange, si, at, atm, atc
        real(kind=dp) :: wcot, wtan, StS
        real(kind=dp) :: t1, t2, xi(1:4), yi(1:4), fk_r, fk_i
        
        real(kind=dp), allocatable, dimension(:) :: tan2t, st2, pf, mu1
        real(kind=dp) :: mu0(1:3)
        real(kind=dp), allocatable :: dmudQ(:,:)
        complex(kind=dpc), allocatable, dimension(:) :: Tg, Cg, DC
        complex(kind=dpc), allocatable, dimension(:,:) :: Xm, Te, TH
        complex(kind=dpc), allocatable, dimension(:,:) :: Ym, Ce, TT
        complex(kind=dpc), allocatable, dimension(:) :: Xmu1, Ymu1, Zmu1, XTeD, gmu, gmu1, gmu2
        complex(kind=dpc), allocatable, dimension(:) :: Y, phi, thz1
        complex(kind=dpc), allocatable, dimension(:) :: G, sdphi, detphi
        integer :: i, k, m, p, q, N, nvib, imu

        ! The dimension of the problem (the number of vibrations)
        nvib = size(D)
        ! beta Temperature
        beta = 1.0_dp/(bkcm*T)
        ! Length of the Fourier transform N = 2^(p+1); NTAU = 2^p
        N = 2*ntau;
        OmegaRange = 0.5_dp*OmR;
        domega = OmR/real(N,dp)
        dtau = TWOPI/OmR
        TauRange = ntau*dtau

        ! Write some information on a log file
        open(flog,file='dft.log',status='unknown')
        write(flog,*) 'N ', N
        write(flog,*) 'Omega range ', OmegaRange
        write(flog,*) 'domega ', domega
        write(flog,*) 'dtau ', dtau
        write(flog,*) 'Tau range ', TauRange

        allocate(f(1:N),gmu(1:ntau),gmu1(1:ntau),gmu2(1:ntau))
        allocate(G(1:ntau),sdphi(1:ntau),detphi(1:ntau))
        allocate(pf(1:nvib))
        allocate(DC(1:nvib))
        allocate(mu1(1:nvib),dmudQ(1:3,1:nvib))
        gmu = zero
        gmu1 = zero; gmu2 = zero

        ! Assign transition moment (tm) to local variable to be used in OMP directives
        mu0 = tm%mu0
        dmudQ = tm%dmudQ
        ! mu1 ian auxiliary vector
        mu1 = zero
        do imu = 1, 3
            mu1(1:nvib) = mu1(1:nvib) + mu0(imu)*dmudQ(imu,1:nvib)
        end do

        ! Normal mode shift as a complex number
        DC = cmplx(D,zero,dpc)
        pf = (one-exp(-beta*wg))**2

        gmu = zero
        call cpu_time(t1)
        !$OMP PARALLEL  DEFAULT(none) PRIVATE(k,tau,z,p,q,m,StS,wtan,wcot,at,atm,atc, &
        !$OMP& tan2t,thz1,Xm,Te,Tg,Ym,Ce,Cg,TH,TT,st2,phi,Y,Xmu1,Ymu1,Zmu1,XTeD)  &
        !$OMP& SHARED(sdphi,G,S,D,DC,dtau,wg,we,nvib,beta,pf,ntau,detphi,gmu, &
        !$OMP& gmu1,gmu2,dmudQ,mu0,mu1)
        allocate(Xm(1:nvib,1:nvib),Tg(1:nvib),Te(1:nvib,1:nvib))
        allocate(Ym(1:nvib,1:nvib),Cg(1:nvib),Ce(1:nvib,1:nvib))
        allocate(TH(1:nvib,1:nvib),TT(1:nvib,1:nvib))
        allocate(Y(1:nvib),phi(1:nvib))
        allocate(tan2t(1:nvib), st2(1:nvib), thz1(1:nvib))
        allocate(Xmu1(1:nvib),Ymu1(1:nvib),Zmu1(1:nvib),XTeD(1:nvib))
        !$OMP DO
        do k = 1, ntau-1
            tau = (real(k,dp))*dtau
            z = cmplx(beta,-tau,dpc)

            tan2t = tan((tau/2.0_dp)*we)
            thz1 = tanh((z/2.0_dp)*wg)
            Tg = wg*thz1
            Cg = wg/thz1

            do p = 1, nvib
                do q = p, nvib
                    at = zero; atm = zero; atc = zero
                    do m = 1, nvib
                        StS = S(m,p)*S(m,q)
                        wtan = we(m)*tan2t(m)
                        wcot = we(m)/tan2t(m)
                        at = at + StS*wtan
                        atc = atc + StS*wcot
                        atm = atm + StS/wtan
                    end do
                    Te(p,q) = cmplx(zero,at,dpc)
                    Ce(p,q) = cmplx(zero,-atc,dpc)
                    TT(p,q) = cmplx(zero,-atm,dpc)
                    if (p /= q) then
                        Te(q,p) = Te(p,q)
                        Ce(q,p) = Ce(p,q)
                        TT(q,p) = TT(p,q)
                    end if
                end do
            end do
            ! Xm = Te + Tg
            Xm = Te
            forall(p=1:nvib)
                Xm(p,p) = Xm(p,p) + Tg(p)
            end forall
            ! Ym = Ce + Cg
            Ym = Ce
            forall(p=1:nvib)
                Ym(p,p) = Ym(p,p) + Cg(p)
            end forall
            ! TT = Te^{-1} + Tg^{-1}
            forall(p=1:nvib)
                TT(p,p) = TT(p,p) + one/Tg(p)
            end forall
            ! Start calculation of Phi and its determinant
            TH = matmul(Xm,Ym)
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)/(wg(p)*we(q))
            end forall
            st2 = sin(tau*we)
            phi = 0.5_dp*cmplx(st2*sin(tau*wg),st2*cos(tau*wg),dpc)* &
            (one-exp(-2.0_dp*z*wg))/pf
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)*phi(q)
            end forall
            detphi(k) = det(TH)
            sdphi(k) = 1.0_dp/sqrt(detphi(k))

            Y = zinvmv(TT,DC)
            G(k) = cfac*dot_product(D,Y)

            XTeD = zinvmv(Xm,matmul(Te,D))
            gmu1(k) = gmu1(k) + 2.0_dp * dot_product(mu1,XTeD)

            do imu = 1, 3
                Xmu1 = zinvmv(Xm,cmplx(dmudQ(imu,:),zero,dpc))
                Ymu1 = zinvmv(Ym,cmplx(dmudQ(imu,:),zero,dpc))
                Zmu1 = Xmu1 - Ymu1
                gmu2(k) = gmu2(k) + 0.5_dp*dot_product(dmudQ(imu,:),Zmu1)/cfac
                do p=1,nvib
                    do q=1,nvib
                        gmu2(k) = gmu2(k)+ dmudQ(imu,p)*dmudQ(imu,q)*XTeD(p)*XTeD(q)
                    end do
                end do
            end do
            !print *, sum(mu0**2), gmu1(k), gmu2(k)
            gmu(k) = cmplx(sum(mu0**2),zero) + gmu1(k) + gmu2(k)

        end do
        !$OMP END DO
        deallocate(Xm,Tg,Te,Ym,Cg,Ce)
        deallocate(TH,TT,Y,phi)
        deallocate(tan2t, st2, thz1)
        !$OMP END PARALLEL
        call cpu_time(t2)

        write(flog,*) 'Time ', t2-t1
        !------------------------------------------------------------------!
        ! Check the sign of the square root: this is fundamental           !
        ! we must always follow the same sheet of the complex square root  !
        !------------------------------------------------------------------!
        si = one
        do k = 2, ntau-1
            if (real(detphi(k)) < 0 .and. real(detphi(k-1)) < 0) then
                if (aimag(detphi(k-1))*aimag(detphi(k)) < 0) si = -si
            end if
            ! Set the sign of the square root
            !sdphi(k) = sign(sdphi(k),si) not working with complex numbers
            sdphi(k) = si*sdphi(k) ! Set the sign of the square root
        end do

        !------------------------------------------------------------------!
        ! Compute Generating Function and apply a window (apodiziation)
        !------------------------------------------------------------------!
        ! f(1) is the generating function at tau = 0.
        f(1) = cmplx(sum(mu0**2),zero,dpc)
        allocate(Tg(1:nvib))
        Tg = (1.0/wg)*cosh((beta/2.0_dp)*wg)/sinh((beta/2.0_dp)*wg)/cfac
        do i = 1, 3
            f(1) = f(1) + 0.5_dp*dot_product(dmudQ(i,:),Tg*dmudQ(i,:))
        end do

        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k) SHARED(f,sdphi,G,ntau,N,gmu)
        do k = 1, ntau-1
            !do k = 1, N
            ! Use Lorentzian window
            !f(ntau+k) = sdphi(k)*exp(-G(k))*EXPL(k,N)
            ! Use Blackman-harris window
            !f(k+1) = sdphi(k)*exp(-G(k))*gmu(k)*BH4(k+1,N) ! this is the working expression for HT effects
            f(k+1) = sdphi(k)*exp(-G(k))*BH4(k+1,N)  ! this can be used when we have no HT effects
	        ! Use Kaiser-Bessel window
	        !f(ntau+k) = sdphi(k)*exp(-G(k))*KBESS(k,N)
	        ! No window
	        !f(k+1) = sdphi(k)*exp(-G(k))*gmu(k)
        end do
        !$OMP END PARALLEL DO

        ! This is a workaround to remove singularities. It works only if a NaN exception
        ! is catched by the compiler.
        xi(1:4) = (/ 0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp/)
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k,yi,fk_r,fk_i) SHARED(f,ntau,xi)
        do k = 3, ntau-2
            if (ieee_is_nan(real(f(k))) .or. ieee_is_nan(aimag(f(k)))) then
                print *, 'interpolating at', k
                yi(1:4) = (/real(f(k-2)), real(f(k-1)), real(f(k+1)), real(f(k+2))/)
                fk_r = cubic_interp(2.0_dp, xi, yi)
                yi(1:4) = (/aimag(f(k-2)), imag(f(k-1)), imag(f(k+1)), imag(f(k+2))/)
                fk_i = cubic_interp(2.0_dp, xi, yi)
                f(k) = cmplx(fk_r, fk_i,dpc)
            end if
        end do
        !$OMP END PARALLEL DO

        if (fwhm > 0.0) then
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k) SHARED(f,ntau,N,fwhm)
            do k = 1, ntau
                ! Add Lorentzian broadening
                f(k) = f(k)*expl(k,N,fwhm)
            end do
        !$OMP END PARALLEL DO
        end if

        deallocate(G,sdphi,detphi)

        do k = 1, ntau
           write(44,*) real(f(k)), imag(f(k))
        end do

        ! Use complex conjugate symmetry
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(f,ntau,N)
        do i= 1, ntau
            f(modulo(N-i+1,N)+1) = conjg(f(i))
        end do
        !$OMP END PARALLEL DO

        allocate(x(1:N))

        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(x,OmegaRange,domega,N)
        do i = 1, N
            x(i) = -OmegaRange +real(i,dp)*domega
        end do
        !$OMP END PARALLEL DO

        ! Shift the origin
        x = x + (sum(wg-we))/2.0_dp

        write(97,*)'origin shift', sum(wg-we)/2.0_dp
        return
    end subroutine genfun_ht

    !==========================================================================!

    function expl(k,N,gamma) result(pw)

        ! Lorentzian window
        integer, intent(in) :: k, N
        real(kind=dp), intent(in) :: gamma
        real(kind=dp) :: rk, rN
        real(kind=dp) :: pw

        rN = real(N,dp)
        rk = real(k,dp)-one-rN/2.0_dp

        pw = exp(gamma*rk/rN)

    end function expl

    !==========================================================================!

    function BH4(k,N) result(bh)

        ! Four term Blackman-Harris window
        ! The window has the maximum at the origin (k=0).

        integer, intent(in) :: k, N
        real(kind=dp) :: bh, rk, rN
        real(kind=dp) :: a(0:3)
        parameter (a=(/ 0.35875_dp, 0.48829_dp, 0.14128_dp, 0.01168_dp/))
        !parameter (b=(/ 0.40217_dp, 0.49703_dp, 0.09392_dp, 0.00183_dp/))

        rN = real(N,dp)
        rk = real(k,dp)-one-rN/2.0_dp

        bh = a(0) - a(1)*cos(TWOPI*rk/rN)+a(2)*cos(4.0_dp*PI*rk/rN)-a(3)*cos(6.0_dp*PI*rk/rN)

    end function BH4

    !==========================================================================!

    function kbess(k,N) result(bh)

        ! Kaiser-Bessel Window
        ! The window has the maximum at the origin (k=0).

        integer, intent(in) :: k, N
        real(kind=dp) :: bh, rk, rN, alpha, ax
        parameter (alpha = 3.0_dp)
        real(kind=dp), save :: bsnorm
        logical, save :: first = .true.
        real(kind=dp), external :: dbesi0

        !$OMP SINGLE
        if (first) then
            bsnorm = dbesi0(pi*alpha)
            first = .false.
        end if
        !$OMP END SINGLE

        rN = real(N,dp)
        rk = real(k,dp)

        ax = alpha*sqrt(one-(2.0_dp*rk/rN)**2)
        bh = dbesi0(pi*ax)/bsnorm

    end function kbess

    !==========================================================================!

    function cubic_interp(xi, x, y) result(yi)

        real(dp), intent(in) :: xi, x(1:4), y(1:4)
        real(dp) :: yi, a(1:4)
        real(dp) :: V(1:4,1:4)
        integer :: i

        do i = 1, 4
            V(i,1:4) = (/one, x(i), x(i)**2, x(i)**3/)
        end do

        a = matmul(invm(V),y)

        yi = a(1)+a(2)*xi+a(3)*(xi**2)+a(4)*(xi**3)

    end function cubic_interp

    !==========================================================================!

    subroutine kubo_lineshape(incvib,transf,kopt,freq_initial_state,freq_final_state,tm)

        use, intrinsic :: iso_c_binding 
        include 'fftw3.f03'
        !---------------------------------------------------------------------------------------------------
        ! NOTE: the density of states obtained by FFT is normalized in such a way that \sum_i \spec_i = N
        ! however we need it to be \Delta E \sum_i \rho_i = 1 (in the program Delta E = domega)
        ! Thus the real density \rho_i is obtained by  \rho_i = \spec_i/omr. (domega = omr/N)
        !---------------------------------------------------------------------------------------------------
        !I/O variables
        type(incvib_t), intent(in) :: incvib
        type(transf_t), intent(in) :: transf
        type(kubo_t), intent(in) :: kopt
        type(tm_t), intent(in), optional :: tm
        real(kind=dp), intent(in out) ::  freq_initial_state(:), freq_final_state(:)

        !Local variables
        complex(kind=dpc), allocatable :: spec(:)
        real(kind=dp), allocatable :: S(:,:), D(:), freq_ini(:), freq_fin(:)
        real(kind=dp) :: Temp, OmR, fwhm, srs, domega
        integer :: ntau, N, nvib, i, status, k
        type(C_PTR) :: plan

        ntau = 2**kopt%pow
        N = 2*ntau
        nvib = size(incvib%id)

        Temp = kopt%Temp
        OmR = kopt%OmR
        fwhm = kopt%fwhm
        domega = OmR/real(N,dp)

        allocate(D(1:nvib), S(1:nvib,1:nvib))
        S = zero; D = zero
        S = transf%JM(incvib%id,incvib%id)
        !D = -matmul(transpose(S),transf%KM(incvib%id)) ! note the sign
        ! Above I should be the inverse of S, not the transpose. Why am I using the
        ! transpose???? 'cause I am an idiot!!
        ! It is more correct to use:
        D = -matmul(invm(S),transf%KM(incvib%id)) ! note the sign

        allocate(freq_ini(1:nvib),freq_fin(1:nvib))
        ! f1: frequencies of the ground electronic state
        freq_ini = freq_initial_state(incvib%id)
        ! f2: frequencies of the excited electronic state
        freq_fin = freq_final_state(incvib%id)

		!print *, ' D', D
		!print *, 'F INI ', freq_ini
		!print *, 'F FIN ', freq_fin
		
        ! Compute the Generating Function of the spectral lineshape.
        if (Temp <= zero+epsilon(one)) then
            call genfun_zero(D,S,freq_ini,freq_fin,ntau,OmR,fwhm)
        else
            !if (allocated(tm%dmudQ)) then
            if (present(tm)) then
                !print *, "fcht..."
                call genfun_ht(D,S,freq_ini,freq_fin,Temp,ntau,OmR,fwhm,tm)
            else
                !print *, "fc..."
                call genfun(D,S,freq_ini,freq_fin,Temp,ntau,OmR,fwhm)
                !call genfun_complex(D,S,freq_ini,freq_fin,Temp,ntau,OmR,fwhm)
            end if
        end if

        ! This is equivalent to a circular shift of the inverse FFT
        forall (i=1:N) f(i) = f(i)*(-1.0_dp)**(i-1)

        ! Allocate the array where the spectrum is actually stored
        allocate(spec(1:N))
        spec = cmplx(zero,zero,dpc)

        ! Backward Fourier Transform
        call dfftw_plan_dft_1d(plan,N,f,spec,FFTW_BACKWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,f,spec)
        call dfftw_destroy_plan(plan)

        open(unit=iuspec,file=kopt%file(1:len_trim(kopt%file)),form="formatted",status="unknown")

        ! Write output
        srs = sum(real(spec))
        do k = 1, N
            write(222,*) k, real(f(k)), aimag(f(k))
            ! Write normalized spectrum
            !write(iuspec,*) x(k), real(spec(k))/(srs*real(N,dp))
            write(iuspec,*) x(k), real(spec(k))/omr
            ! Write un-normalized spectrum
            !write(iuspec,*) x(k), real(spec(k))/real(N,dp)
            write(223,*) x(k), real(spec(k)), aimag(spec(k))
        end do

        ! Le due grandezze calcolate sotto sono equivalenti, infatti: domega = omr/N -> domega/omr = 1/N
        write(flog,*) 'Norm (should be 1):', sum(real(spec))/real(N,dp)
        write(flog,*) 'Real spec Integral:', domega*sum(real(spec)/omr)

        return
    end subroutine kubo_lineshape

!------------------------------------------------------------------------------

    subroutine genfun_complex(D,S,wg,we,T,ntau,OmR,fwhm)

        ! wg: frequencies of the INITIAL STATE (groud)
        ! we: frequencies of the FINAL STATE (excited)
        ! S is the Duschinsky matrix
        ! D is the vector of normal mode displacements
        ! T is the temperature
        ! I/O variables
        real(kind=dp), intent(in) :: D(:), S(:,:), wg(:), we(:), T, OmR, fwhm
        integer, intent(in) :: ntau

        ! Local variables
        complex(kind=dpc) :: z, wtanz, wcotz, att, atm, atc, itau
        real(kind=dp) :: beta, tau, dtau, OmegaRange, domega, ge
        real(kind=dp) :: TauRange, si
        real(kind=dp) :: StS
        real(kind=dp) :: t1, t2, xi(1:4), yi(1:4), fk_r, fk_i
        
        real(kind=dp), allocatable, dimension(:) :: tan2t, pf
        complex(kind=dpc), allocatable, dimension(:) :: Tg, Cg, DC, wez, tanh2t, st2
        complex(kind=dpc), allocatable, dimension(:,:) :: Xm, Te, TH
        complex(kind=dpc), allocatable, dimension(:,:) :: Ym, Ce, TT
        complex(kind=dpc), allocatable, dimension(:) :: Y, phi, thz1
        complex(kind=dpc), allocatable, dimension(:) :: G, sdphi, detphi
        integer :: i, k, m, p, q, N, nvib

        write(fout,*)"Using complex plane shift for the calculation of the inverse matrices in KUBO."
        ! The dimension of the problem (the number of vibrations)
        nvib = size(D)
        ! beta Temperature
        beta = 1.0_dp/(bkcm*T)
        ! Length of the Fourier transform N = 2^(p+1)
        N = 2*ntau;
        OmegaRange = 0.5_dp*OmR;
        domega = OmR/real(N,dp)
        dtau = TWOPI/OmR
        TauRange = ntau*dtau

        ! Write some information on a log file
        open(flog,file='dft.log',status='unknown')
        write(flog,*) 'N ', N
        write(flog,*) 'Omega range ', OmegaRange
        write(flog,*) 'domega ', domega
        write(flog,*) 'dtau ', dtau
        write(flog,*) 'Tau range ', TauRange

        allocate(f(1:N))
        allocate(G(1:ntau),sdphi(1:ntau),detphi(1:ntau))
        allocate(pf(1:nvib))
        allocate(DC(1:nvib))
        allocate(wez(1:nvib))

        DC = cmplx(D,zero,dpc)
        pf = (one-exp(-beta*wg))**2
        ! complex constant (natural line width in cm-1)
        ge = 0.1_dp
        ! Shift the energies in the complex plane with an imaginary component "ge"
        wez =  cmplx(we,-ge)

        call cpu_time(t1)
        !$OMP PARALLEL  DEFAULT(none) PRIVATE(k,itau,tau,z,p,q,m,StS,wtanz,wcotz,att,atm,atc,tanh2t,thz1,Xm,Te,Tg,Ym,Ce,Cg,TH,TT,st2,phi,Y) SHARED(sdphi,G,S,D,DC,dtau,wg,we,ge,wez,nvib,beta,pf,ntau,detphi)
        allocate(Xm(1:nvib,1:nvib),Tg(1:nvib),Te(1:nvib,1:nvib))
        allocate(Ym(1:nvib,1:nvib),Cg(1:nvib),Ce(1:nvib,1:nvib))
        allocate(TH(1:nvib,1:nvib),TT(1:nvib,1:nvib))
        allocate(Y(1:nvib),phi(1:nvib))
        allocate(st2(1:nvib),thz1(1:nvib),tanh2t(1:nvib))
        !$OMP DO
        do k = 1, ntau-1
            tau = (real(k,dp))*dtau
            z = cmplx(beta,-tau,dpc)
            itau = cmplx(zero,tau,dpc)

            tanh2t = tanh((itau/2.0_dp)*wez)
            thz1 = tanh((z/2.0_dp)*wg)
            Tg = wg*thz1
            Cg = wg/thz1

            do p = 1, nvib
                do q = p, nvib
                    att = cmplx(zero,zero,dpc); atm = cmplx(zero,zero,dpc); atc = cmplx(zero,zero,dpc)
                    do m = 1, nvib
                        StS = S(m,p)*S(m,q)
                        wtanz = wez(m)*tanh2t(m)
                        wcotz = wez(m)/tanh2t(m)
                        att = att + StS*wtanz
                        atc = atc + StS*wcotz
                        atm = atm + StS/wtanz
                    end do
                    Te(p,q) = att
                    Ce(p,q) = atc
                    TT(p,q) = atm  ! TT = Te^-1
                    if (p /= q) then
                        Te(q,p) = Te(p,q)
                        Ce(q,p) = Ce(p,q)
                        TT(q,p) = TT(p,q)
                    end if
                end do
            end do
            ! Xm = Te + Tg
            Xm = Te
            forall(p=1:nvib)
                Xm(p,p) = Xm(p,p) + Tg(p)
            end forall
            ! Ym = Ce + Cg
            Ym = Ce
            forall(p=1:nvib)
                Ym(p,p) = Ym(p,p) + Cg(p)
            end forall
            ! TT = Te^{-1} + Tg^{-1}
            forall(p=1:nvib)
                TT(p,p) = TT(p,p) + one/Tg(p)
            end forall
            ! Start calculation of Phi and its determinant
            TH = matmul(Xm,Ym)
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)/(wg(p)*wez(q))
            end forall
            st2 = sinh(itau*wez)
            !phi = 0.5_dp*cmplx(st2*sin(tau*wg),st2*cos(tau*wg),dpc)*(one-exp(-2.0_dp*z*wg))/pf
            phi = 0.5_dp*sinh(itau*wez)*exp(-itau*wg)*(one-exp(-2.0_dp*z*wg))/pf  ! check the sign: (it should be exp(itau*wg)...)
            forall(p=1:nvib,q=1:nvib)
                TH(p,q)=TH(p,q)*phi(q)
            end forall
            detphi(k) = det(TH)
            sdphi(k) = 1.0_dp/sqrt(detphi(k))

            Y = zinvmv(TT,DC)
            G(k) = cfac*dot_product(D,Y)
        end do
        !$OMP END DO
        deallocate(Xm,Tg,Te,Ym,Cg,Ce)
        deallocate(TH,TT,Y,phi)
        deallocate(tanh2t, st2, thz1)
        !$OMP END PARALLEL
        call cpu_time(t2)

        write(flog,*) 'Time ', t2-t1
        !------------------------------------------------------------------!
        ! Check the sign of the square root: this is fundamental           !
        ! we must always follow the same sheet of the complex square root  !
        !------------------------------------------------------------------!
        si = one
        do k = 2, ntau-1
            if (real(detphi(k)) < 0 .and. real(detphi(k-1)) < 0) then
                if (aimag(detphi(k-1))*aimag(detphi(k)) < 0) si = -si
            end if
            !if (si < 0) sdphi(k) = si*sdphi(k) ! Set the sign of the square root
            sdphi(k) = si*sdphi(k) ! Set the sign of the square root
        end do

        !------------------------------------------------------------------!
        ! Compute Generating Function and apply a window (apodiziation)
        !------------------------------------------------------------------!
        f(1) = cmplx(one,zero,dpc)
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k) SHARED(f,sdphi,G,ntau,N)
        do k = 1, ntau-1
            ! Use Blackman-harris window
            !f(ntau+k) = sdphi(k)*exp(-G(k))*BH4(k,N)
            !f(k+1) = sdphi(k)*exp(-G(k))*BH4(k+1,N)
            ! No window-use with caution
            f(k+1) = sdphi(k)*exp(-G(k))
        end do
        !$OMP END PARALLEL DO

        ! This is a workaround to remove singularities. It works only if a NaN exception
        ! is catched by the compiler.
        xi(1:4) = (/ 0.0_dp, 1.0_dp, 3.0_dp, 4.0_dp/)
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(k,yi,fk_r,fk_i) SHARED(f,ntau,xi)
        do k = 3, ntau-2
            if (ieee_is_nan(real(f(k))) .or. ieee_is_nan(aimag(f(k)))) then
                yi(1:4) = (/real(f(k-2)), real(f(k-1)), real(f(k+1)), real(f(k+2))/)
                fk_r = cubic_interp(2.0_dp, xi, yi)
                yi(1:4) = (/aimag(f(k-2)), imag(f(k-1)), imag(f(k+1)), imag(f(k+2))/)
                fk_i = cubic_interp(2.0_dp, xi, yi)
                f(k) = cmplx(fk_r, fk_i,dpc)
            end if
        end do
        !$OMP END PARALLEL DO

        deallocate(G,sdphi,detphi)

        ! Use complex conjugate symmetry
        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(f,ntau,N)
        do i= 1, ntau
            f(modulo(N-i+1,N)+1) = conjg(f(i))
        end do
        !$OMP END PARALLEL DO
        do i = 1, N
          write(771,*)f(i)
        end do

        allocate(x(1:N))

        !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(x,OmegaRange,domega,N)
        do i = 1, N
            x(i) = -OmegaRange +real(i,dp)*domega
        end do
        !$OMP END PARALLEL DO

        ! Shift the origin
        x = x + (sum(wg-we))/2.0_dp


        write(97,*)'origin shift', sum(wg-we)/2.0_dp
        return
    end subroutine genfun_complex

end module kubo
