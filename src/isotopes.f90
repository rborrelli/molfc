module isotopes

  use parameters
  use util
  use system_type, only : element
  use errors

  implicit none

  type (element), allocatable, public :: elmnts (:)

  character (len=Symlen), private :: PeriodicTable(1:NELE)
  data PeriodicTable(1:NELE) /'H', 'D',  'He', &
                              'Li', 'Be', 'B',  'C', 'N', 'O', 'F', 'Ne', &
                              'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', &
                              'K', 'Ca', 'Fe' ,'Co', 'Ni', 'Cu', 'Zn', 'Br', 'Br81', 'Ge'/

  public :: init_isotopes, & 
            free_isotopes, & 
            mapel
  
  interface mapel
    module procedure mapel_sym
    module procedure mapel_z
  end interface mapel

  contains
  
    subroutine init_isotopes
    
    integer iel
    integer, parameter :: X = 6
    
    allocate(elmnts(1:NELE)) ! Allocate an array of elements.
    
    iel = mapel('H ')
    elmnts(iel)%Z   =   1                   
    elmnts(iel)%N   =   0                   
    elmnts(iel)%AN  =   1                   
    elmnts(iel)%Sym = 'H '                  
    elmnts(iel)%AM  =     1.007825032
    elmnts(iel)%val  =     1
    elmnts(iel)%rad  =     0.23
                                  
    iel = mapel('D ')
    elmnts(iel)%Z   =   1                   
    elmnts(iel)%N   =   1                   
    elmnts(iel)%AN  =   2                   
    elmnts(iel)%Sym = 'D '                  
    elmnts(iel)%AM  =     2.014101778       
    elmnts(iel)%val  =     1
    elmnts(iel)%rad  =     0.23
    
    iel = mapel('He')
    elmnts(iel)%AN  =   4                   
    elmnts(iel)%N   =   2                   
    elmnts(iel)%Z   =   2                   
    elmnts(iel)%Sym = 'He'                  
    elmnts(iel)%AM  =     4.002603250       
    elmnts(iel)%val  =     0    ! Never bonded
    elmnts(iel)%rad  =     1.50 ! Never bonded
    
    iel = mapel('Li')
    elmnts(iel)%Z   =   3                   
    elmnts(iel)%N   =   3                   
    elmnts(iel)%AN  =   6                   
    elmnts(iel)%Sym = 'Li'                  
    elmnts(iel)%AM  =     6.015122281       
    elmnts(iel)%val  =     3   
    elmnts(iel)%rad  =     0.68 
    
    iel = mapel('Be')
    elmnts(iel)%Z   =   4                   
    elmnts(iel)%N   =   3                   
    elmnts(iel)%AN  =   7                   
    elmnts(iel)%Sym = 'Li'                  
    elmnts(iel)%AM  =     9.012 !...very approximate...
    elmnts(iel)%val  =     3
    elmnts(iel)%rad  =     0.35
    
    iel = mapel('B ')
    elmnts(iel)%Z   =   5                   
    elmnts(iel)%N   =   5                   
    elmnts(iel)%AN  =   10                   
    elmnts(iel)%Sym = 'B '                  
    elmnts(iel)%AM  =     10.811 !...very approximate...
    elmnts(iel)%val  =     3
    elmnts(iel)%rad  =     0.83

    iel = mapel('C ')
    elmnts(iel)%Z   =   6                   
    elmnts(iel)%N   =   6                   
    elmnts(iel)%AN  =  12                   
    elmnts(iel)%Sym = 'C '                  
    elmnts(iel)%AM  =    12.000000000       
    elmnts(iel)%val  =     4
    elmnts(iel)%rad  =     0.68
    
    iel = mapel('N ')
    elmnts(iel)%Z   =   7                   
    elmnts(iel)%N   =   7                   
    elmnts(iel)%AN  =  14                   
    elmnts(iel)%Sym = 'N '                  
    elmnts(iel)%AM  =    14.003074005       
    elmnts(iel)%val  =     5
    elmnts(iel)%rad  =     0.68
    
    iel = mapel('O ')
    elmnts(iel)%Z   =   8                   
    elmnts(iel)%N   =   8                   
    elmnts(iel)%AN  =  16                   
    elmnts(iel)%Sym = 'O '                  
    elmnts(iel)%AM  =    15.994914622       
    elmnts(iel)%val  =     2
    elmnts(iel)%rad  =     0.68
    
    iel = mapel('F ')
    elmnts(iel)%Z   =   9                   
    elmnts(iel)%N   =  10 
    elmnts(iel)%AN  =  19                 
    elmnts(iel)%Sym = 'F '                  
    elmnts(iel)%AM  =  18.998403205       
    elmnts(iel)%val  =     1
    elmnts(iel)%rad  =     0.64

    iel = mapel('Ne')
    elmnts(iel)%Z   =  10                   
    elmnts(iel)%N   =  10                   
    elmnts(iel)%AN  =  20                   
    elmnts(iel)%Sym = 'Ne'                  
    elmnts(iel)%AM  =  19.9924401754      
    elmnts(iel)%val  =     0
    elmnts(iel)%rad  =     1.50

    iel = mapel('Na')
    elmnts(iel)%Z   =  11                   
    elmnts(iel)%N   =  12                   
    elmnts(iel)%AN  =  23                   
    elmnts(iel)%Sym = 'Na'                  
    elmnts(iel)%AM  =  22.989769281       
    elmnts(iel)%val  =     1
!    elmnts(iel)%rad  =     0.97
    elmnts(iel)%rad  =     1.90  ! the stadard value 0.97 gives some problem with Na+ in DNA

    iel = mapel('Mg')
    elmnts(iel)%Z   =  12
    elmnts(iel)%N   =  12
    elmnts(iel)%AN  =  24
    elmnts(iel)%Sym = 'Mg'
    elmnts(iel)%AM  =  23.985041699       
    elmnts(iel)%val  =     2
    elmnts(iel)%rad  =     1.20 ! the stadard value 1.10 gives some problem with Mg++ in chlorophylls

    iel = mapel('Al')
    elmnts(iel)%Z   =  13
    elmnts(iel)%N   =  14
    elmnts(iel)%AN  =  27
    elmnts(iel)%Sym = 'Al'
    elmnts(iel)%AM  =  26.981538631      
    elmnts(iel)%val  =     3
    elmnts(iel)%rad  =     1.35

    iel = mapel('Si')
    elmnts(iel)%Z   =  14
    elmnts(iel)%N   =  14
    elmnts(iel)%AN  =  28
    elmnts(iel)%Sym = 'Si'
    elmnts(iel)%AM  =  27.976926532
    elmnts(iel)%val  =     4
    elmnts(iel)%rad  =     1.20

    iel = mapel('P ')
    elmnts(iel)%Z   =  15
    elmnts(iel)%N   =  16
    elmnts(iel)%AN  =  31
    elmnts(iel)%Sym = 'P '
    elmnts(iel)%AM  =  30.973761632
    elmnts(iel)%val  =     5
    elmnts(iel)%rad  =     1.05

    iel = mapel('S ')
    elmnts(iel)%Z   =  16
    elmnts(iel)%N   =  16
    elmnts(iel)%AN  =  32
    elmnts(iel)%Sym = 'S '
    elmnts(iel)%AM  =  31.972071001
    elmnts(iel)%val  =     2
    elmnts(iel)%rad  =     1.02

    iel = mapel('Cl')
    elmnts(iel)%Z   =  17                   
    elmnts(iel)%N   =  18                   
    elmnts(iel)%AN  =  35                   
    elmnts(iel)%Sym = 'Cl'
    elmnts(iel)%AM  =  34.978852684
    elmnts(iel)%val  =     1
    elmnts(iel)%rad  =     0.99

    iel = mapel('Ar')
    elmnts(iel)%Z   =  18
    elmnts(iel)%N   =  22
    elmnts(iel)%AN  =  40
    elmnts(iel)%Sym = 'Ar'
    elmnts(iel)%AM  =  39.962383123       
    elmnts(iel)%val  =     0
    elmnts(iel)%rad  =     1.51

    iel = mapel('K ')
    elmnts(iel)%Z   =  19                   
    elmnts(iel)%N   =  10                   
    elmnts(iel)%AN  =  19                   
    elmnts(iel)%Sym = 'F '                  
    elmnts(iel)%AM  =    18.998403205       
    elmnts(iel)%val  =     1
    elmnts(iel)%rad  =     1.33

    iel = mapel('Ca')
    elmnts(iel)%Z   =  20                   
    elmnts(iel)%N   =  10                   
    elmnts(iel)%AN  =  19                   
    elmnts(iel)%Sym = 'F '                  
    elmnts(iel)%AM  =    18.998403205       
    elmnts(iel)%val  =     2
    elmnts(iel)%rad  =     0.99

    iel = mapel('Fe')                          
    elmnts(iel)%Z   =  26                   
    elmnts(iel)%N   =  30                   
    elmnts(iel)%AN  =  56                   
    elmnts(iel)%Sym = 'Fe'                  
    elmnts(iel)%AM  =    55.934942133       
    elmnts(iel)%val  =     X
    elmnts(iel)%rad  =     1.34

    iel = mapel('Co')                          
    elmnts(iel)%Z   =  27
    elmnts(iel)%N   =  32
    elmnts(iel)%AN  =  59
    elmnts(iel)%Sym = 'Co'
    elmnts(iel)%AM  =  58.933195070
    elmnts(iel)%val  =     X
    elmnts(iel)%rad  =     1.33

    iel = mapel('Ni')                          
    elmnts(iel)%Z   =  28
    elmnts(iel)%N   =  30
    elmnts(iel)%AN  =  58
    elmnts(iel)%Sym = 'Ni'
    elmnts(iel)%AM  =  57.935342973
    elmnts(iel)%val  =     X
    elmnts(iel)%rad  =     1.50

    iel = mapel('Cu')
    elmnts(iel)%Z   =  29
    elmnts(iel)%N   =  34
    elmnts(iel)%AN  =  63
    elmnts(iel)%Sym = 'Fe'
    elmnts(iel)%AM  =  62.92959756
    elmnts(iel)%val  =     X
    elmnts(iel)%rad  =     1.52

    iel = mapel('Zn')
    elmnts(iel)%Z   =  30
    elmnts(iel)%N   =  34
    elmnts(iel)%AN  =  64
    elmnts(iel)%Sym = 'Zn'
    elmnts(iel)%AM  =  63.92914482
    elmnts(iel)%val  =     X
    elmnts(iel)%rad  =     1.45

	iel = mapel('Br')
    elmnts(iel)%Z   =  35
    elmnts(iel)%N   =  44
    elmnts(iel)%AN  =  79
    elmnts(iel)%Sym = 'Br'
    elmnts(iel)%AM  =  78.918337122

    iel = mapel('Br81')
    elmnts(iel)%Z   =  35
    elmnts(iel)%N   =  46
    elmnts(iel)%AN  =  81
    elmnts(iel)%Sym = 'Br'
    elmnts(iel)%AM  =  80.916290621

    iel = mapel('Ge')
    elmnts(iel)%Z   =  32
    elmnts(iel)%N   =  42
    elmnts(iel)%AN  =  74
    elmnts(iel)%Sym = 'Ge'
    elmnts(iel)%AM  =  73.9211788 
    elmnts(iel)%rad  =     1.45


    return
    end subroutine init_isotopes 
    
    !------------------------------------------------------
    
    subroutine free_isotopes
        
    integer alloc_err, ierr

    deallocate(elmnts,stat=alloc_err)
    if (alloc_err /= 0) ierr=error(0,"Cannot free memory from element array.")
    
    return
    end subroutine free_isotopes

    !------------------------------------------------------

    integer function mapel_sym(Sym)
    
    character (len = *), intent(in) :: Sym
    character (len = Symlen) TSym
    character c1, c2
    integer i, ierr 
    logical found_element
    
    
    if (len(Sym) == 0) return 
    
    TSym = adjustl(Sym)
    c1 = TSym(1:1)
    TSym(1:1) = ToUpper (c1)

    if (len(Sym) > 1 ) then 
      c2 = TSym(2:2)
      TSym(2:2) = ToLower (c2)
    end if

    found_element = .false.
    do i = 1, NELE
      if (PeriodicTable(i) == TSym) then
        mapel_sym = i 
        found_element = .true.
        exit
      end if
    end do
    
    if (.not.found_element) mapel_sym = error(25,error_string=TSym)

    return
    end function mapel_sym

!!! =============================================================

    integer function mapel_z(Z)
    
    integer, intent(in) :: Z
    integer i, ierr 
    logical found_element
    
    
    if (Z == 0) ierr = error(0,'No element with atomic number: 0')
    
    found_element = .false.
    do i = 1, NELE
      if (elmnts(i)%Z == Z) then
        mapel_z = i 
        found_element = .true.
        exit
      end if
    end do
    
    if (.not.found_element) then 
		mapel_z = 0
		ierr = error(0,'Cannot found element with atomic number')
	end if

    return
    end function mapel_z

end module isotopes
