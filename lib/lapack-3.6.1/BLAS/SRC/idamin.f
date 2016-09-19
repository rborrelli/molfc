      integer function idamin(n,dx,incx)
c
c     finds the index of element having min. absolute value.
c
      double precision dx(*),dmin 
      integer i,incx,ix,n

      idamin = 0
      if( n .lt. 1 ) return
      idamin = 1
      if(n.eq.1) return
      if(incx.eq.1) go to 30
c
c        code for increment not equal to 1
c
      ix = 1
      dmin = abs(dx(1))
      ix = ix + incx
      do 20 i = 2,n
         if( abs(dx(ix)).ge.dmin) go to 10
         idamin = i
         dmin = abs(dx(ix))
   10    ix = ix + incx
   20 continue
      return
c
c        code for increment equal to 1
c
   30 dmin = abs(dx(1))
      do 40 i = 2,n
         if( abs(dx(i)).ge.dmin) go to 40
         idamin = i
         dmin = abs(dx(i))
   40 continue
      return
      end 
