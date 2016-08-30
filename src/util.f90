module util

  use parameters
  use errors

  implicit none

  private

  public :: ToUpper, ToLower, openfileprob, isempty, isequal, extract
  
  interface extract
    module procedure extract_c
  end interface extract

  interface isempty
    module procedure isempty_char
  end interface isempty

  interface ToUpper
!    MODULE PROCEDURE ToUpper
    module procedure ToUpper_s
  end interface ToUpper

  contains 

!==============================================================================!

!    function ToUpper (Lchar) result (Uchar)
!    
!    character Lchar, Uchar
!    integer iu1, iu2, il1, il2, isu
!
!    iu1 = ichar('a'); iu2 = ichar('z')
!    il1 = ichar('A'); il2 = ichar('Z')
! 
!    isu = ichar(Lchar)
!
!    Uchar = Lchar
!
!    if ( isu >= iu1 .and. isu <= iu2 ) then
!       isu = isu - iu1
!       Uchar = char(il1+isu)
!    end if
!
!    return
!    end function ToUpper
    
!==============================================================================!

    function ToUpper_s (Lchar) result (Uchar)
    
    character(len=*) Lchar
    character(len=len(Lchar)) Uchar
    integer i, iu1, iu2, il1, il2, isu

    iu1 = ichar('a'); iu2 = ichar('z')
    il1 = ichar('A'); il2 = ichar('Z')
 
    do i = 1, len(Lchar)
    
    isu = ichar(Lchar(i:i))

    Uchar(i:i) = Lchar(i:i)

    if ( isu >= iu1 .and. isu <= iu2 ) then
       isu = isu - iu1
       Uchar(i:i) = char(il1+isu)
    end if

    end do
    
    return
    end function ToUpper_s

!==============================================================================!

    function ToLower (Uchar) result (Lchar)
   
    character Lchar, Uchar
    integer iu1, iu2, il1, il2, isu

    iu1 = ichar('a'); iu2 = ichar('z')
    il1 = ichar('A'); il2 = ichar('Z')
    isu = ichar(Uchar)

    Lchar = Uchar

    if ( isu >= il1 .and. isu <= il2 ) then
       isu = isu - il1
       Lchar = char(iu1+isu)
    end if

    return
    end function ToLower

!==============================================================================

    subroutine openfileprob ( k, file, filedir ) 

    character (len = *)              file
    character (len = 72)             filenew 
    character (len = *), optional :: filedir 
    integer i1, i2, i3, ia, k, kk
    
    i1 = int(k/100)
    ia = k-i1*100
    i2 = int(ia/10)
    i3 = ia-i2*10
    i1 = i1+48
    i2 = i2+48
    i3 = i3+48
    
    filenew = 'prob_'//char(i1)//char(i2)//char(i3)//'_'//&
             &file(1:len_trim(file))//'.dat'

    if (present(filedir)) then
     filenew=filedir(1:len_trim(filedir))//filenew
    end if

    kk = k +15 
    open (unit=kk,file=filenew(1:len_trim(filenew)),status='unknown')

    return
    end subroutine openfileprob 

!==============================================================================!
    
    function isempty_char(string)

    character (len = *) string
    logical isempty_char

    isempty_char = .false.
    if (len_trim(string) == 0) isempty_char = .true.

    return
    end function isempty_char

!==============================================================================!

    function isequal(string1, string2)

    character (len = *) string1, string2
    logical isequal
    integer l1, l2, i

    isequal = .true.

    string1 = adjustl(string1)
    string2 = adjustl(string2)

    l1 = len_trim(string1); l2 = len_trim(string2)

    if (l1 /= l2) then ! le stringhe sono diverse 
      isequal = .false. 
      return  
    end if

    do i = 1, l1  ! l1 e l2 sono uguali!
      if (string1(i:i) /= string2(i:i)) then 
        isequal = .false.
        return
      end if
    end do

    return
    end function isequal

    pure function extract_c(c,start,finish)

      implicit none
      character(*), intent(in)                  :: c
      integer, intent(in)                       :: start,finish
      character(len_extract_c(c,start,finish))  :: extract_c
      integer                                   :: is,if


      is = max(1,start)
      if = min(len(c),finish)
      if (if < is) then
          extract_c = ''
      else
          extract_c(1:if-is+1) = c(is:if)
      endif

      end function extract_c

      elemental function len_extract_c(c,start,finish)

      implicit none
      character(*), intent(in)  :: c
      integer, intent(in)       :: start,finish
      integer                   :: len_extract_c
      integer                   :: is,if


      is = max(1,start)
      if = min(len(c),finish)
      if (if < is) then
          len_extract_c = 0
      else
          len_extract_c = max(0,if-is) + 1
      endif

      end function len_extract_c

end module util
