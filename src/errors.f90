module errors

  use parameters, only : fout, MAXFC

  implicit none

  public :: error, warning

  contains

    integer function error(error_type,error_string,error_iflag)

    integer, intent(in) :: error_type
    integer, optional :: error_iflag
    character(len=*), optional :: error_string

    select case (error_type)
      case (0)
        write(fout,*) '==================================='
        write(fout,*) 'ERROR:                             '
        write(fout,*) error_string
        write(fout,*) '==================================='
        stop
      case (1)
        write(fout,*) '==================================='
        write(fout,*) 'ERROR                              '
        write(fout,*) 'Unable to open file ', error_string,' for structure ', error_iflag
        write(fout,*) 'Check if it exist.'
        write(fout,*) '==================================='
        stop
      case (4)
        write(fout,*) '===================================='
        write(fout,*) 'ERROR                               '
        write(fout,*) 'Unable to open file 4 for vibrations' 
        write(fout,*) 'Check if it exist.                  '
        write(fout,*) '===================================='
        stop
      case (5)
        write(fout,*) '========================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'The program has read in different numbers'
        write(fout,*) 'of atoms for the two structures.         '
        write(fout,*) '========================================='
        stop
      case (6)
        write(fout,*) '===================================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'Unable to open the file containing the initial state'
        write(fout,*) 'Check if it exist.'
        write(fout,*) '===================================================='
        stop
      case (7)
        write(fout,*) '=========================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'Unable to re-allocate memory for molecular'
        write(fout,*) 'structure.  Probably this is a BUG!    '
        write(fout,*) '=========================================='
        stop
      case (8)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'Unable to allocate memory for FC matrices'
        write(fout,*) 'structure.                                  '
        write(fout,*) '============================================'
        stop
      case (9)
        write(fout,*) '================================================'
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'Unable to allocate memory for T1 and T2 matrices'
        write(fout,*) 'structure.  Probably this is a BUG!         '
        write(fout,*) '================================================'
        stop
      case (10)
        write(fout,*) '====================================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'Unable to allocate memory for GEO1 and GEO2 matrices'
        write(fout,*) 'structure.  Probably this is a BUG!         '
        write(fout,*) '====================================================='
        stop
      case (15)
        write(fout,*) '===================================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'The two structures have a different number of atoms'
        write(fout,*) 'Check the input.  '
        write(fout,*) '===================================================='
        stop
      case (16)
        write(fout,*) '===================================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'The number of normal modes is non consistent with the'
        write(fout,*) 'number of atoms in the molecule.'
        write(fout,*) 'Check the input.  '
        write(fout,*) '===================================================='
        stop
      case (17)
        write(fout,*) '===================================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'Not enough argument on command line.'
        write(fout,*) 'getarg error.'
        write(fout,*) '===================================================='
        stop
      case (18)
        write(fout,*) '===================================================='
        write(fout,*) 'Input and output files have the same name - exiting' 
        write(fout,*) '===================================================='
        stop
      case (19)
        write(fout,*) '========================================================='
        write(fout,*) 'You must specify the output file name after the -o option' 
        write(fout,*) '========================================================='
        stop
     case (168)
        write(fout,*) '========================================================='
        write(fout,*) 'Cannot allocate arrays in exact_evolution.               ' 
        write(fout,*) '========================================================='
        stop
      case (20)
        write(fout,*) '==================================='
        write(fout,*) 'ERROR                            '
        write(fout,*) 'Cannot find ',error_string,' keyword' 
        write(fout,*) 'Check the input file.'
        write(fout,*) '==================================='
        stop
      case (21)
        write(fout,*) '==================================='
        write(fout,*) 'ERROR parsing input stream.        '
        write(fout,*) '==================================='
        stop
      case (22)
        write(fout,*) '===================================================='
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'The ranges in <probability> are not correct.'
        write(fout,*) 'Check the input.  '
        write(fout,*) '===================================================='
        stop
      case (25)
        write(fout,*) '========================================================'
        write(fout,*) 'ERROR                                    '
        write(fout,*) 'Cannot find element ', error_string ,' in isotopes module'
        write(fout,*) 'You may want to update your library. Check the internal '
        write(fout,*) 'documentation.'
        write(fout,*) '========================================================'
        stop
      case (65)
        write(fout,*) '==================================='
        write(fout,*) 'ERROR reading vibrations.          '
        write(fout,*) '==================================='
        stop
      case (66)
        write(fout,*) '===================================='
        write(fout,*) 'ERROR: number of molecules in <state>'
        write(fout,*) ' must be the same.                  '
        write(fout,*) '===================================='
        stop
      case (67)
        write(fout,*) '==================================='
        write(fout,*) 'ERROR reorder structures in input. '
        write(fout,*) '==================================='
        stop
      case (68)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: increase the number of states for the'
        write(fout,*) 'mode ', error_iflag
        write(fout,*) '============================================'
        stop
      case (69)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: '
        write(fout,*) "No element input with the name: ", error_string
        write(fout,*) "is available. Double check your input!"
        write(fout,*) '============================================'
        stop
      case (70)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: '
        write(fout,*) "No attributes named ", error_string
        write(fout,*) "or wrog position. Check the input keywords!"
        write(fout,*) '============================================'
        stop
      case (71)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: '
        write(fout,*) 'Cannot allocate  ', error_string
        write(fout,*) '============================================'
        stop
      case (72)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: '
        write(fout,*) 'Cannot deallocate  ', error_string
        write(fout,*) '============================================'
        stop
      case (88)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: '
        write(fout,*) 'Cannot open file ', error_string,'!'
        write(fout,*) '============================================'
        stop
      case (190)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: '
        write(fout,*) 'Array ', error_string,' not allocated. '
        write(fout,*) 'THIS IS A BUG.'
        write(fout,*) '============================================'
        stop
      case (191)
        write(fout,*) '============================================'
        write(fout,*) 'ERROR: '
        write(fout,*) 'Maximum number of FC allowed:', MAXFC
        write(fout,*) '============================================'
        stop

      case default
        write(fout,*) 'UNKNOWN ERROR'
        stop

    end select 
    
    error = 1 

    return
    end function error

    ! --------------------------------------------------------------------------

    integer function warning(warning_type,warning_string,warning_iflag)
    integer, intent(in) :: warning_type
    integer, optional :: warning_iflag
    character(len=*), optional :: warning_string

    select case (warning_type)
      case (0)
        write(fout,*) '======================================'
        write(fout,*) '              WARNING!!!!'
        write(fout,*) warning_string
        write(fout,*) 
        write(fout,*) '======================================'

      case (1)
        write(fout,*) '======================================'
        write(fout,*) '              WARNING!!!!'
        write(fout,*) 'Found distances  > 1.0 Ang.'
        write(fout,*) 'Double Check your Molecular Structures'
        write(fout,*) '======================================'

      case (2)
        write(fout,*) '======================================'
        write(fout,*) '              WARNING!!!!'
        write(fout,*) 'Hamiltonian is huge.       '
        write(fout,*) 'If you really want to write it add '
        write(fout,*) 'force=".true." in control             '
        write(fout,*) '======================================'

      case (3)
        write(fout,*) '======================================'
        write(fout,*) '              WARNING!!!!'
        write(fout,*) 'The inversion of P matrix gives '
        write(fout,*) 'complex eigenvalues. '
        write(fout,*) '======================================'

      case (4)
        write(fout,*) '======================================'
        write(fout,*) 'You are skipping some elements of the'
        write(fout,*) 'Hamiltonian matrix. '
        write(fout,*) '======================================'

      case default
        write(fout,*) 'UNKNOWN WARNING'

    end select 
    
    warning = 1 

    return
    end function warning

    !==============================================================================

    subroutine die(str)
    character(len=*), intent(in), optional   :: str
    if (present(str)) then
       write(unit=0,fmt="(a)") trim(str)
    endif
    write(unit=0,fmt="(a)") "Stopping Program"
    stop
    end subroutine die

end module errors
