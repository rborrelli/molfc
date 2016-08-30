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
   if (index(akt,'k') /= 0 ) then
     if (icoord > 0) write(19,'(2x,a)')"</ncoord>"
     read(line,*,iostat=ios)ivib,akt,c,ctyp
     select case (ctyp)
       case("stre") 
         write(19,'(2x,a)')"<ncoord type=""s"">"
       case("bend") 
         write(19,'(2x,a)')"<ncoord type=""b"">"
       case("tors") 
         write(19,'(2x,a)')"<ncoord type=""d"">"
       case("out","outp") 
         write(19,'(2x,a)')"<ncoord type=""w"">"
       case("linc") 
         write(19,'(2x,a)')"<ncoord type=""linc"">"
       case("linp") 
         write(19,'(2x,a)')"<ncoord type=""linp"">"
     end select
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
   case("stre") 
      read(tline,*)i,j
      ircoord = ircoord + 1 
      !write(9,'(2(i3,","))')ircoord,icoord
      write(19,'(4x,a,f15.12,a,i3,i3,a)') "<coord type=""s"" c=""",c,""">",i,j," </coord>" 
      !write(19,'(f15.12,",")') c
      !print *, 's ', i,j
   case("bend") 
      read(tline,*)i,j,k
      ircoord = ircoord + 1 
      write(*,'(4(i3,","))')2,i,k,j
      !write(9,'(2(i3,","))')ircoord,icoord
      write(19,'(4x,a,f15.12,a,i3,i3,i3,a)') "<coord type=""b"" c=""",c,""">",i,k,j," </coord>" 
      !print *, 'b ', i,j,k
   case("linc") 
      read(tline,*)i,j,k,l
      ircoord = ircoord + 1 
      write(*,'(4(i3,","))')2,i,k,j,l
      !write(9,'(2(i3,","))')ircoord,icoord
      write(19,'(4x,a,f15.12,a,4i3,a)') "<coord type=""linc"" c=""",c,""">",i,k,j,l," </coord>" 
      !print *, 'linc ', i,j,k
   case("linp") 
      read(tline,*)i,j,k,l
      ircoord = ircoord + 1 
      write(*,'(4(i3,","))')2,i,k,j,l
      !write(9,'(2(i3,","))')ircoord,icoord
      write(19,'(4x,a,f15.12,a,4i3,a)') "<coord type=""linp"" c=""",c,""">",i,k,j,l," </coord>" 
      !print *, 'linp ', i,j,k,l
   case("tors") 
      read(tline,*)i,j,k,l
      ircoord = ircoord + 1 
      write(*,'(5(i3,","))')3,i,j,k,l
      write(19,'(4x,a,f15.12,a,i3,i3,i3,i3,a)') "<coord type=""d"" c=""",c,""">",i,j,k,l," </coord>" 
      !print *, 't ', i,j,k,l
   case("outp","out") 
      read(tline,*)i,j,k,l
      ircoord = ircoord + 1 
      write(*,'(5(i3,","))')4,i,l,k,j
      write(19,'(4x,a,f15.12,a,i3,i3,i3,i3,a)') "<coord type=""w"" c=""",c,""">",i,j,k,l," </coord>" 
      !print *, 'p ', i,j,k,l
 end select
 end if

end do

2 if (icoord > 0) write(19,'(2x,a)')"</ncoord>"

!print *, 'coord ', ircoord, icoord

end program int2fc
 
