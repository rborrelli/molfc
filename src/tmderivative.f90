module tm_derivative

    use parameters
    use transf_type
    use xmvar, only : system

    implicit none

    private

    real(kind=dp), public :: mu0(1:3)
    real(kind=dp), public, allocatable :: dmudQ(:,:)
    public :: ht

contains

    !------------------------------------------------------------------!
    ! The input transition moments are defined with respect to the     !
    ! excited state equilibrium position and excited state vibrational !
    ! coordinate. We transform it to the ground state coordinates to   !
    ! perform the calculation of the generating function.              !
    !------------------------------------------------------------------!

    subroutine ht(transformation)

        type(transf_t), intent(in) :: transformation
        integer :: i, j, iat, k, N, iatom, ixi, istep, nvib
        logical :: logex
        integer, parameter :: futm = 19
        real(kind=dp) :: hstep, mass, l_iat_k_j, mu1, mu2, mu3, conv
        real(kind=dp), allocatable :: dmudxi(:,:,:), mu(:,:,:,:)

        ! Number of atoms and vibrations
        N = system%state(1)%molecule%structure%numat
        nvib = system%state(1)%molecule%nvib

        if (system%model) then

            allocate(dmudQ(1:3,1:nvib))
            allocate(system%tm%dmuds(1:3,1:nvib))

            dmudQ(1:3,1:nvib) = system%tm%dmudQ(1:3,1:nvib)
            mu0 = system%tm%mu0

            ! Transform the transition dipole moment derivative to the ground state coordinate
            ! and assign the values to the globale variable system%tm
            if (transformation%tm_rotate) then
                write(fout,*) 'TM rotated to the ground state coordinates.'
                do ixi = 1, 3
                    system%tm%dmudQ(ixi,1:nvib) = matmul(transpose(transformation%JM), dmudQ(ixi,1:nvib))
                    system%tm%mu0(ixi) = mu0(ixi) + dot_product(dmudQ(ixi,1:nvib),transformation%KM)
                end do
                ! TM w.r.t. dimensionless coordinates
                do i = 1, nvib
                    conv = system%state(2)%molecule%normodes%vibration(i)%freq*cfac
                    system%tm%dmuds(1:3,i) = system%tm%dmudQ(1:3,i)/sqrt(conv)
                end do
            else
                write(fout,*) 'TM defined w.r.t. the excited state coordinates.'
                ! TM w.r.t. dimensionless coordinates
                do i = 1, nvib
                    conv = system%state(1)%molecule%normodes%vibration(i)%freq*cfac
                    system%tm%dmuds(1:3,i) = system%tm%dmudQ(1:3,i)/sqrt(conv)
                end do
            end if

        else

            inquire (file="tm.xyz", exist=logex)
            if (.not.logex) then
                write(*,*)'Cannot find tm.xyx file. Cannot compute Herzberg-Teller effects.'
                write(*,*)'Execution aborted.'
                stop
            end if

            allocate(mu(1:3,1:N,1:3,1:2))
            allocate(dmudxi(1:3,1:N,1:3))
            dmudxi = zero

            ! Open the transition dipole moment
            open(unit=futm,file='tm.xyz',status='old')

            read(futm,*)hstep
            read(futm,*)mu0(1),mu0(2),mu0(3)
            mu = zero
            do i = 1, 6*N
                read(futm,*)iatom, ixi, istep, mu1, mu2, mu3
                mu(1,iatom,ixi,istep) = mu1
                mu(2,iatom,ixi,istep) = mu2
                mu(3,iatom,ixi,istep) = mu3
            end do

            do ixi = 1, 3
                do iat = 1, n
                    mass =  system%state(1)%molecule%structure%atom(iat)%elem%AM
                    do k = 1, 3
                        dmudxi(ixi,iat,k) = (mu(ixi,iat,k,1) - mu(ixi,iat,k,2))/(2.0_dp*hstep)
                        dmudxi(ixi,iat,k) = dmudxi(ixi,iat,k) / sqrt(mass)
                    end do
                end do
            end do

            allocate(dmudQ(1:3,1:nvib))
            dmudQ = zero

            ! Compute the transition dipole derivatives with respect to the excited state
            ! normal coordinates.
            do ixi = 1, 3
                do j = 1, nvib
                    do iat = 1, N
                        do k = 1, 3
                            l_iat_k_j = system%state(1)%molecule%normodes%vibration(j)%atom(iat)%d(k)
                            dmudQ(ixi,j) = dmudQ(ixi,j) + l_iat_k_j*dmudxi(ixi,iat,k)
                        end do
                    end do
                end do
            end do
            deallocate(dmudxi)

            allocate(system%tm%dmudQ(1:3,1:nvib))
            allocate(system%tm%dmuds(1:3,1:nvib))
            system%tm%dmudQ = zero
            if (transformation%tm_rotate) then
                write(fout,*) 'TM rotated to the ground state coordinates.'
                ! Transform the transition dipole moment derivative to the ground state coordinate
                ! and assign the values to the globale variable system%tm
                do ixi = 1, 3
                    system%tm%dmudQ(ixi,1:nvib) = matmul(transpose(transformation%JM), dmudQ(ixi,1:nvib))
                    system%tm%mu0(ixi) = mu0(ixi) + dot_product(dmudQ(ixi,1:nvib),transformation%KM)
                    !print *, 'dmudQ', ixi
                    !print *,  system%tm%dmudQ(ixi,1:nvib)
                    !print *, 'mu0', ixi, system%tm%mu0(ixi)
                end do
                ! TM w.r.t. dimensionless coordinates
                do i = 1, nvib
                    conv = system%state(2)%molecule%normodes%vibration(i)%freq*cfac
                    system%tm%dmuds(:,i) = system%tm%dmudQ(:,i)/sqrt(conv)
                end do
            else
                write(fout,*) 'TM defined w.r.t. the excited state coordinates.'
                do ixi = 1, 3
                    system%tm%dmudQ(ixi,1:nvib) = dmudQ(ixi,1:nvib)
                    system%tm%mu0(ixi) = mu0(ixi)
                end do
                ! TM w.r.t. dimensionless coordinates
                do i = 1, nvib
                    conv = system%state(1)%molecule%normodes%vibration(i)%freq*cfac
                    system%tm%dmuds(:,i) = system%tm%dmudQ(:,i)/sqrt(conv)
                end do
            end if

        end if

        ! Write tm derivatives
        write(fout,'(//,6x,a,i1,a)')'Transition Moment d/dQ_i mu'
        write(fout,'(4x,a)')'===================================================='
        write(fout,'(6x,a)')'          X               Y                 Z    '
        write(fout,'(4x,a)')'===================================================='
        write(fout,'(2x,i3,2x,3(f14.4,2x))')0,(system%tm%mu0(ixi),ixi=1,3)
        do i = 1, nvib
            write(fout,'(2x,i3,2x,3(f14.4,2x))')i,(system%tm%dmudQ(ixi,i),ixi=1,3)
        end do

    end subroutine ht

end module tm_derivative
