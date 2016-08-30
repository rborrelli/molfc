module reorg

    use parameters
    use system_type
    use transf_type

    implicit none

    private

    public :: reorganization_energies

    CONTAINS 

    !==============================================================================!

    SUBROUTINE reorganization_energies(molecule,transf)
    
    !---------------------------------------+
    ! Calcola l'energia di riorganizzazione |
    !---------------------------------------+

    ! Energia dello stato finale nel punto di minimo dello stato iniziale.
    ! nota che: stato finale = lo stato 1, stato iniziale = stato 2.

    type(molecule_t) :: molecule
    type(transf_t) :: transf
    real(kind=dp), allocatable :: f1(:), d(:)
    real(kind=dp) reorg, reorgm, eph
    integer i
    
    allocate(f1(1:size(transf%KM)))
    allocate(d(1:size(transf%KM)))
    f1 = molecule%normodes%vibration(:)%freq ! FINAL_STATE frequencies
    d = transf%KM*sqrt(cfac*f1) ! dimensionless displacements

    reorg = 0.d0; reorgm = 0.0_dp
    
    reorg = 0.5_dp*sum(cfac*(transf%KM*f1)**2)
    write(fout,2111)
    do i = 1, size(transf%KM)
      reorgm = 0.5_dp*(d(i)**2)*f1(i)
      eph = d(i)*f1(i)
      write(fout,2110) i, f1(i), reorgm, reorgm*1000/ev, eph !/sqrt(2.0_dp)
    end do
    write(fout,2113) reorg, reorg*1000/ev
    
    
2111 format(/,2x,'-----------------------',&
            /,2x,'Reorganization Energies',&
            /,2x,'-----------------------',&
            /,20x,'cm-1      meV       e-ph ')
2110 format(2x,i3,f8.2,2x,f10.3,2x,f10.5,2x,f8.1)
2113 format(/,2x,'-----------------------',&
            /,2x,'Total Reorganization Energy',&
            /,2x,'---------------------------',&
            /,2x,'Lambda ',f10.3,3x,f10.3/)

    return
    END SUBROUTINE reorganization_energies

end module reorg


