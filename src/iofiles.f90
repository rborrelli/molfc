module iofiles
  
  use parameters
  use util
  use errors
  
  implicit none

  private
  
  public :: openfile, newfile

  contains
  
!==============================================================================!

    subroutine openfile ( iunit, prefix, suffix, extension, filedir ) 

    implicit none
    
! --- Variabile statica che conta il numero di files aperti.
    integer, save :: nopfile = 0

! --- Variabili di input-output della subroutine
    integer, intent(in)       :: iunit
    character (len = *)           :: prefix
    character (len = *), optional :: suffix
    character (len = *), optional :: extension
    character (len = *), optional :: filedir
    
! --- Variabili interne
    integer i, ierr
    character (len = 80) filenew, gfile
    character newchar
    logical logex

    logex = .false.
    
    if (present(suffix)) then
      filenew = prefix(1:len_trim(prefix))//'_'//suffix(1:len_trim(suffix))
    else
      filenew = prefix(1:len_trim(prefix))
    end if

    if (present(filedir)) then
      if (.not.isempty(filedir)) then
        if (filedir(len_trim(filedir):len_trim(filedir)) /= '/') & 
          filedir=filedir(1:len_trim(filedir))//'/'
        filenew=adjustl(filedir(1:len_trim(filedir)))//adjustl(filenew(1:len_trim(filenew)))
      end if
    end if

    if (present(extension)) then
      filenew=filenew(1:len_trim(filenew))//extension(1:len_trim(extension))
    end if

    inquire (file=filenew, exist=logex)

    if (.not.logex) ierr = error(88,filenew)
    
    filenew = filenew(1:len_trim(filenew))
    filenew = adjustl(filenew)
    open (unit=iunit,file=filenew(1:len_trim(filenew)),status='old')    

    return
    end subroutine openfile

!==============================================================================!

    subroutine newfile ( iunit, prefix, suffix, extension, filedir ) 

    implicit none
    
! --- Variabile statica che conta il numero di files aperti.
    integer, save :: nopfile = 0

! --- Variabili di input-output della subroutine
    integer, intent(in)       :: iunit
    character (len = *)           :: prefix
    character (len = *), optional :: suffix
    character (len = *), optional :: extension
    character (len = *), optional :: filedir
    
! --- Variabili interne
    integer i
    character (len = 80) filenew, gfile
    character newchar
    logical logex, logex2

    logex = .false.; logex2 = .false.
    
    if (present(suffix)) then
      filenew = prefix(1:len_trim(prefix))//'_'//suffix(1:len_trim(suffix))
    else
      filenew = prefix(1:len_trim(prefix))
    end if

    if (present(filedir)) then
      if (.not.isempty(filedir)) then
        if (filedir(len_trim(filedir):len_trim(filedir)) /= '/') & 
          filedir=filedir(1:len_trim(filedir))//'/'
        filenew=adjustl(filedir(1:len_trim(filedir)))//adjustl(filenew(1:len_trim(filenew)))
      end if
    end if

    if (present(extension)) then
      filenew=filenew(1:len_trim(filenew))//extension(1:len_trim(extension))
    end if

    inquire (file=filenew, exist=logex)

    ! Se il file esiste e lo si vuole aprire come nuovo allora viene rinominato.
    if (logex) then
     do i = 65, 90 
      newchar=char(i)
      gfile = filenew(1:len_trim(filenew))//'_'//newchar(1:len_trim(newchar))
      inquire (file=gfile, exist=logex2)
      if (.not.logex2) then
        open (unit=iunit,file=gfile(1:len_trim(gfile)),status='new')
        exit
      end if
     end do
    else
      gfile = filenew(1:len_trim(filenew))
      open (unit=iunit,file=gfile(1:len_trim(gfile)),status='new')
    end if

    nopfile = nopfile + 1
    write(fout,'(2x,i3,''.'',x,a)') nopfile, gfile
    
    if (logex2) stop 'Cannot open new files in your directory.'

    include 'formats'

    return
    end subroutine newfile
  
end module iofiles
