! Transform Turbomole input of internal coordinates to MolFC style.
program int2fc

implicit none
integer :: ios, ivib, i,j,k,l,ind,icoord, ircoord
character(len=100) ::  line, tline
character(len=4) :: akt,ctyp
real(8) :: c

! IRCOORD counts the total redundants internal coordinates
! ICOORD counts only the non-redundant coordinates
icoord = 0
ircoord = 0
do 
 
 read(*,'(a100)',end=2)line
 read(line,*,iostat=ios)ivib,akt
 if (ios == 0) then
   !print *, 'ios', ivib, akt
   ! IF we cannot find the string 'K' or 'F' the coordinate is part
   ! of a linear combination
   if (index(akt,'K') /= 0 ) then
     if (icoord > 0) write(19,'(2x,a)')"</ncoord>"
     read(line,*,iostat=ios)ivib,akt,c,c,ctyp
     if (ios == 0) then
       write(19,'(2x,a)')"<ncoord type=""d"">"
     else
       read(line,*,iostat=ios)ivib,akt,c,ctyp
       select case (ctyp)
         case("STRE") 
           write(19,'(2x,a)')"<ncoord type=""s"">"
         case("BEND") 
           write(19,'(2x,a)')"<ncoord type=""b"">"
         case("TORS") 
           write(19,'(2x,a)')"<ncoord type=""d"">"
         case("OUT","OUTP") 
           write(19,'(2x,a)')"<ncoord type=""w"">"
       end select
     end if
     icoord = icoord + 1 
   !print *, 'y-k:', c,ctyp
   else
     read(line,*,iostat=ios)c,ctyp
   !print *, 'n-k', c,ctyp
   end if
 ind = index(line,ctyp)+4
 tline=adjustl(line(ind:))
 !print *, ind
 !print *, 'tline', tline
 ! now we have "akt","c",and "ctype" hence we know what to read and how
 select case (ctyp)
   case("STRE") 
      read(tline,*)i,j
      ircoord = ircoord + 1 
      !write(9,'(2(i3,","))')ircoord,icoord
      write(19,'(4x,a,f15.12,a,i3,i3,a)') "<coord type=""s"" c=""",c,""">",i,j," </coord>" 
      !write(19,'(f15.12,",")') c
      !print *, 's ', i,j
   case("BEND") 
      read(tline,*)i,j,k
      ircoord = ircoord + 1 
      write(*,'(4(i3,","))')2,i,k,j
      !write(9,'(2(i3,","))')ircoord,icoord
      write(19,'(4x,a,f15.12,a,i3,i3,i3,a)') "<coord type=""b"" c=""",c,""">",i,k,j," </coord>" 
      !print *, 'b ', i,j,k
   case("TORS") 
      read(tline,*)i,j,k,l
      ircoord = ircoord + 1 
      write(*,'(5(i3,","))')3,i,j,k,l
      write(19,'(4x,a,f15.12,a,i3,i3,i3,i3,a)') "<coord type=""d"" c=""",c,""">",i,j,k,l," </coord>" 
      !print *, 't ', i,j,k,l
   case("OUTP","OUT") 
      read(tline,*)i,j,k,l
      ircoord = ircoord + 1 
      write(*,'(5(i3,","))')4,i,l,k,j
      write(19,'(4x,a,f15.12,a,i3,i3,i3,i3,a)') "<coord type=""w"" c=""",c,""">",i,j,k,l," </coord>" 
      !print *, 'p ', i,j,k,l
 end select
 end if

end do

2 if (icoord > 0) write(19,'(2x,a)')"</ncoord>"

!print *, 'COORD ', ircoord, icoord

end program int2fc
 
