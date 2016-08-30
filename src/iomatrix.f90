module iomatrix

  use parameters
  
  private

  public :: matout, layout, write_octave_data

  interface write_octave_data
    module procedure write_octave_data
    module procedure write_octave_data_vec
  end interface write_octave_data

  interface layout
    module procedure layout
    module procedure layout_int
  end interface layout

  interface matout
    module procedure matout
    module procedure matout2
    module procedure vecout
  end interface matout
  
  contains

      subroutine matout(aa,nc,nr,ndim) 

      implicit real(kind=dp)(a-h,o-z)
      dimension aa(ndim,nc)

      ka=1
      kc=6
   10 kb=min0(kc,nc)
      write(fout,50) (i,i=ka,kb)
      write(fout,70)
      la=1
      lc=40
   20 lb=min0(lc,nr)
      n=0
      do 30 i=la,lb
         write(fout,80) i,(aa(i,j),j=ka,kb)
         n=n+1
         if (n.lt.10) go to 30
!         write(fout,70)
         n=0
   30 continue
      if (lb.eq.nr) go to 40
      la=lc+1
      lc=lc+40
!      write(fout,90)
      go to 20
   40 if (kb.eq.nc) return
      ka=kc+1
      kc=kc+6
!      if (nr.gt.25) write(fout,90)
      go to 10
!
   50 format (/,2x,'Mode',3x,i5,9i12)
   60 format (/4x,10f12.6)
   70 format (/)
   80 format (i5,10f12.5)
   90 format (1h1)

      return
      end subroutine matout

!-------------------------------------------------------------
      subroutine matout2(aa,bb,nc,nr,ndim) 

      implicit real(kind=dp)(a-h,o-z)
      dimension aa(ndim,nc), bb(nc)


      ka=1
      kc=6
   10 kb=min0(kc,nc)
      write(fout,50) (i,i=ka,kb)
      write(fout,70)
      write(fout,79) (bb(j),j=ka,kb)
      write(fout,70)
      la=1
      lc=40
   20 lb=min0(lc,nr)
      n=0
      do 30 i=la,lb
         write(fout,80) i,(aa(i,j),j=ka,kb)
         n=n+1
         if (n.lt.10) go to 30
!         write(fout,70)
         n=0
   30 continue
      if (lb.eq.nr) go to 40
      la=lc+1
      lc=lc+40
!      write(fout,90)
      go to 20
   40 if (kb.eq.nc) return
      ka=kc+1
      kc=kc+6
!      if (nr.gt.25) write(fout,90)
      go to 10
!
   50 format (/,2x,'Root',3x,i5,9i12)
   60 format (/4x,10f12.6)
   70 format (/)
   79 format (5x,10f12.1,//)
   80 format (i5,10f12.5)
   90 format (1h1)

      return
      end subroutine matout2

!===========================================================================

        subroutine vecout(a,n)

        real(kind=dp), intent(in) :: a(1:n)
        integer i

        write(fout,'(2x,a)') 'Vector :'
        do i = 1, size(a)
          write(fout,'(2x,f12.6)') a(i)
        end do

        return
        end subroutine vecout

!===========================================================================

        subroutine layout(a,km,kn,m,n)
!
!       layout prints the lower triangle of a square matrix a,
!       or the entire matrix if a is rectangular.
!
        implicit real(kind=dp)(a-h,o-z)
        real(kind=dp), intent(in) :: a(km,kn)

        lhi=0
!        if(m.eq.n) then
!1          low=lhi+1
!           write(fout,200) (jcol,jcol=low,min0(low+8,n))
!           do irow=low,n
!             lhi=min0(irow,low+8)
!             write(fout,201) irow,(a(jcol,irow),jcol=low,min0(low+8,irow))
!           end do
!           if(lhi.lt.n) go to 1
!        else
           do
!3            low=lhi+1
3            low=lhi+1
             lhi=min0(low+5,n)
             write(fout,200) (jcol,jcol=low,lhi)
             write(fout,'(1h )')
             do irow=1,m
               write(fout,201) irow,(a(irow,jcol),jcol=low,lhi)
             end do
!             if(lhi.lt.n) go to 3
             if(lhi.ge.n) exit
           end do
!        endif
        return

    200 format(/,8x,9i12)
    201 format(i5,3x,9f12.5)
        end subroutine layout

!===========================================================================

        subroutine layout_int(a,km,kn,m,n)
!
!       layout prints the lower triangle of a square matrix a,
!       or the entire matrix if a is rectangular.
!
        implicit real(kind=dp)(a-h,o-z)
        integer, intent(in) :: a(km,kn)

        integer, parameter :: nmxc = 15
        lhi=0
        if(m.eq.n) then
1          low=lhi+1
           write(fout,200) (jcol,jcol=low,min0(low+nmxc,n))
           do irow=low,n
             lhi=min0(irow,low+nmxc)
             write(fout,201) irow,(a(jcol,irow),jcol=low,min0(low+nmxc,irow))
           end do
           if(lhi.lt.n) go to 1
        else
           do
!3            low=lhi+1
3            low=lhi+1
             lhi=min0(low+nmxc,n)
             write(fout,200) (jcol,jcol=low,lhi)
             do irow=1,m
               write(fout,201) irow,(a(irow,jcol),jcol=low,lhi)
             end do
!             if(lhi.lt.n) go to 3
             if(lhi.gt.n) exit
           end do
        endif
        return

    200 format(//8x,30(3x,i3)//)
    201 format(i5,3x,30i6)
        end subroutine layout_int

    !==========================================================================!

    subroutine write_octave_data(A,imat)
    
    real(8), intent(in) :: A(:,:)
    integer, intent(in) :: imat
    integer i, j
    
    write(foct,'(a,i1)')'# name: A',imat
    write(foct,'(a)')'# type: matrix'
    write(foct,'(a,1x,i5)')'# rows:',size(A,dim=1)
    write(foct,'(a,1x,i5)')'# columns:',size(A,dim=2)
    
    do i = 1, size(A,dim=1)
       write(foct,*) (A(i,j),j=1,size(A,dim=2))
    end do

    return

    end subroutine write_octave_data

    !==========================================================================!

    subroutine write_octave_data_vec(A,imat)

    real(8), intent(in) :: A(:)
    integer, intent(in) :: imat
    integer i, j

    write(foct,'(a,i1)')'# name: E',imat
    write(foct,'(a)')'# type: matrix'
    write(foct,'(a,1x,i5)')'# rows:',size(A)
    write(foct,'(a)')'# columns: 1'

    do i = 1, size(A)
       write(foct,*) A(i)
    end do

    return
    end subroutine write_octave_data_vec
    
end module iomatrix

