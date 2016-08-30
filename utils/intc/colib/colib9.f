
ctm  stems from cisrt2.f

      subroutine code3x (
     & buffer, intblk, stn121, stn2,
     & stn3,   nintkl, ksym,   k1s,
     & k1f,    l1,     iscr )
c
c  pack the 3-external header info into two working precision words.
c
c  03-jul-01 modified to allow large basis sets (>170). -rls
c
      implicit logical(a-z)
c     # dummy:
      integer intblk, stn121, stn2, stn3, nintkl,
     & ksym, k1s, k1f, l1, iscr
      integer i,j
      logical first

c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
      data first /.true./
      save first
c
*@ifdef int64
*       integer buffer(2)
*       integer oder,links
*@ifdef cray T3E64
*       oder(i,j)=or(i,j)
*       links(i,j)=shiftl(i,j)
*@else
*       oder(i,j)=ior(i,j)
*       links(i,j)=ishft(i,j)
*@endif
*c tm additional change for up to 1023 basis functions
*c  instead k1s 8bit k1f 8bit l1 16bit
*c          k1s 10bit k1f 10bit l1 12bit
*       if (first) then
*       write(0,*) 'more than 255 basis functions not supported in
*     & 64bit'
*       first=.false.
*       endif
*      buffer(1)=oder(links(intblk,48), oder(links(stn121,32),
*     &          oder(links(stn2,  16), stn3)))
*      buffer(2)=oder(links(nintkl,48), oder(links(ksym,  40),
*     &          oder(links(k1s,   32), oder(links(k1f,   24),
*     &          oder(links(l1,     8), iscr)))))
*@elif defined  (oldoldsun) || defined ( fps)
*      integer iword(2)
*      real*8 dword, buffer(2)
*      equivalence (iword(1),dword)
*      integer   lshift, or
*      intrinsic lshift, or
*      iword(1)=or(lshift(intblk,16), stn121)
*      iword(2)=or(lshift(stn2,  16), stn3)
*      buffer(1)=dword
*      iword(1)=or(lshift(nintkl,16), or(lshift(ksym,8), k1s))
*      iword(2)=or(lshift(k1f,   24), or(lshift(l1,  8), iscr))
*      buffer(2)=dword
*@elif defined  obsolete
*c     # 32-bit integer machines with f90 bit operators.
*c     # 16-bit limits on the buffer lengths.
*      integer iword(2)
*      real*8 dword, buffer(2)
*      equivalence (iword(1),dword)
*      integer   ishft, ior
*      intrinsic ishft, ior
*      iword(1)  = ior(ishft(intblk,16), stn121)
*      iword(2)  = ior(ishft(stn2,  16), stn3)
*      buffer(1) = dword
*      iword(1)  = ior(ishft(nintkl,16), ior(ishft(ksym,8), k1s))
*      iword(2)  = ior(ishft(k1f,   24), ior(ishft(l1,  8), iscr))
*      buffer(2) = dword
*@else
c     # 32-bit integer machines with f90 bit operators.
c     # 18-bit limits on the buffer lengths.
      integer iword(2)
      real*8 dword, buffer(2)
      equivalence (iword(1),dword)
c
      integer ib1hi,  ib2hi,  ib3hi,  ib4hi,  ib5hi
      integer ib1low, ib2low, ib3low, ib4low, ib5low
c
      integer   maxsiz,             maxbfn,       mask2
      parameter(maxsiz = 2**18 - 1, maxbfn = 2**10-1, mask2=3)
c
      integer   ishft, iand, ior
      intrinsic ishft, iand, ior
c
c     # check the sizes to make sure everything fits.
c
*@ifdef debug
* 1120  format('code3x:',10i6)
*       write(6,1120) intblk, stn121, stn2,
*     & stn3,   nintkl, ksym,   k1s, k1f, l1, iscr
*@endif

      if ( intblk .gt. maxsiz .or. intblk .lt. 0 ) then
         call bummer(': intblk=', intblk, faterr)
      elseif ( stn121 .gt. maxsiz .or. stn121 .lt. 0 ) then
         call bummer(': stn121=', stn121, faterr)
      elseif ( stn2 .gt. maxsiz .or. stn2 .lt. 0 ) then
         call bummer(': stn2=', stn2, faterr)
      elseif ( stn3 .gt. maxsiz .or. stn3 .lt. 0 ) then
         call bummer(': stn3=', stn3, faterr)
      elseif ( nintkl .gt. maxsiz .or. nintkl .lt. 0 ) then
         call bummer(': nintkl=', nintkl, faterr)
      elseif ( ksym .gt. 8 .or. ksym .lt. 0 ) then
         call bummer(': ksym=', ksym, faterr)
      elseif ( k1s .gt. maxbfn .or. k1s .lt. 0 ) then
         call bummer(': k1s=', k1s, faterr)
      elseif ( k1f .gt. maxbfn .or. k1f .lt. 0 ) then
         call bummer(': k1f=', k1f, faterr)
      elseif ( l1 .gt. maxbfn .or. l1 .lt. 0 ) then
         call bummer(': l1=', l1, faterr)
      elseif ( iscr .gt. 3 .or. iscr .lt. 0 ) then
         call bummer(': iscr=', iscr, faterr)
      endif
c
      ib1hi  = ishft(intblk,-2)
      ib1low = iand(intblk,mask2)
      ib2hi  = ishft(stn121,-2)
      ib2low = iand(stn121,mask2)
      ib3hi  = ishft(stn2,-2)
      ib3low = iand(stn2,mask2)
      ib4hi  = ishft(stn3,-2)
      ib4low = iand(stn3,mask2)
      ib5hi  = ishft(nintkl,-2)
      ib5low = iand(nintkl,mask2)
c
      iword(1)  = ior(ishft(ib1hi,16), ib2hi)
      iword(2)  = ior(ishft(ib3hi,  16), ib4hi)
      buffer(1) = dword
c     iword(1)  = ior(ishft(ib5hi,16), ior(ishft(k1s,8), k1f))
c     iword(2)  = ior(ishft(l1,16), ior(ishft(ksym,12),
c    &            ior(ishft(ib1low,10), ior(ishft(ib2low,8),
c    &            ior(ishft(ib3low,6), ior(ishft(ib4low,4),
c    &            ior(ishft(ib5low,2), iscr)))))))
c
c tm additional change for up to 1023 basis functions
c  instead k1s 8bit k1f 8bit l1 16bit
c          k1s 10bit k1f 10bit l1 12bit
c   k1f is distributed among 2 dwords !
c
      iword(1)  = ior(ishft(ib5hi,16), ior(ishft(k1s,6),
     .                ishft(k1f,-4)))
      iword(2)  = ior(ishft(iand(k1f,15),28),
     &            ior(ishft(l1,16), ior(ishft(ksym,12),
     &            ior(ishft(ib1low,10), ior(ishft(ib2low,8),
     &            ior(ishft(ib3low,6), ior(ishft(ib4low,4),
     &            ior(ishft(ib5low,2), iscr))))))))

      buffer(2) = dword
*@endif
      return
      end

ctm code4x stems from cisrt2.f


      subroutine code4x(
     & buffer, strtnf, stn121, stn3,
     & nintkl, lsym,   ksym,   l1s,
     & l1f,    k1s,    k1f,    k1l1 )
c
c  pack the 4-external header info into two working precision words.
c
c  03-jul-01 modified to allow large basis sets (>170). -rls
c
      implicit logical(a-z)
c     # dummy:
      integer strtnf, stn121, stn3, nintkl, lsym, ksym,
     & l1s, l1f, k1s, k1f, k1l1
      integer i,j
      logical first
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
       data first /.false./
       save first
c
*@ifdef int64
*       integer buffer(2)
*       integer oder,links
*@ifdef cray T3E64
*       oder(i,j)=or(i,j)
*       links(i,j)=shiftl(i,j)
*@else
*       oder(i,j)=ior(i,j)
*       links(i,j)=ishft(i,j)
*@endif
*       if (first) then
*       write(0,*) 'more than 255 basis functions not supported in
*     &  64bit'
*       first=.false.
*       endif
*      buffer(1)=oder(links(strtnf,48), oder(links(stn121,32),
*     &          oder(links(stn3,  16), nintkl)))
*      buffer(2)=oder(links(lsym,  56), oder(links(ksym,  48),
*     &          oder(links(l1s,   40), oder(links(l1f,   32),
*     &          oder(links(k1s,   24),
*     &          oder(links(k1f,   16), k1l1))))))
*@elif defined  (oldoldsun) || defined ( fps)
*      integer iword(2)
*      real*8 dword , buffer(2)
*      equivalence (iword(1),dword)
*      integer   lshift, or
*      intrinsic lshift, or
*      iword(1)=or(lshift(strtnf,16), stn121)
*      iword(2)=or(lshift(stn3,  16), nintkl)
*      buffer(1)=dword
*      iword(1)=or(lshift(lsym,  24), or(lshift(ksym,16),
*     &         or(lshift(l1s,    8), l1f)))
*      iword(2)=or(lshift(k1s,   24), or(lshift(k1f, 16), k1l1))
*      buffer(2)=dword
*@elif defined  obsolete
*c     # 32-bit integer machines with f90 bit operators.
*c     # record length parameters are limited to 16-bit fields.
*      integer iword(2)
*      real*8 dword, buffer(2)
*      equivalence (iword(1),dword)
*      integer   ior, ishft
*      intrinsic ior, ishft
*      iword(1)  = ior(ishft(strtnf,16), stn121)
*      iword(2)  = ior(ishft(stn3,  16), nintkl)
*      buffer(1) = dword
*      iword(1)  = ior(ishft(lsym,  24), ior(ishft(ksym,16),
*     &            ior(ishft(l1s,    8), l1f)))
*      iword(2)  = ior(ishft(k1s,   24), ior(ishft(k1f, 16), k1l1))
*      buffer(2) = dword
*@else
c     # 32-bit integer machines with f90 bit operators.
c     # record length parameters are limited to 18-bit fields.
c
c     # dummy:
      real*8 buffer(3)
c
c     # local:
      integer iword(2)
      real*8 dword
      equivalence (iword(1),dword)
c
      integer ib1hi,  ib2hi,  ib3hi,  ib4hi
      integer ib1low, ib2low, ib3low, ib4low
c
      integer   maxsiz,             maxbfn,       mask2
      parameter(maxsiz = 2**18 - 1, maxbfn = 511, mask2=3)
c
      integer   ior, ishft, iand
      intrinsic ior, ishft, iand
c
      if ( strtnf .gt. maxsiz .or. strtnf .lt. 0 ) then
         call bummer('code4x: strtnf=', strtnf, faterr)
      elseif ( stn121 .gt. maxsiz .or. stn121 .lt. 0 ) then
         call bummer('code4x: stn121=', stn121, faterr)
      elseif ( stn3 .gt. maxsiz .or. stn3 .lt. 0 ) then
         call bummer('code4x: stn3=', stn3, faterr)
      elseif ( nintkl .gt. maxsiz .or. nintkl .lt. 0 ) then
         call bummer('code4x: nintkl=', nintkl, faterr)
      elseif ( lsym .gt. 8 .or. lsym .lt. 1 ) then
         call bummer('code4x: lsym=', lsym, faterr)
      elseif ( ksym .gt. 8 .or. ksym .lt. 1 ) then
         call bummer('code4x: ksym=', ksym, faterr)
      elseif ( l1s .gt. maxbfn .or. l1s .lt. 0 ) then
         call bummer('code4x: l1s=', l1s, faterr)
      elseif ( l1f .gt. maxbfn .or. l1f .lt. 0 ) then
         call bummer('code4x: l1f=', l1f, faterr)
      elseif ( k1s .gt. maxbfn .or. k1s .lt. 0 ) then
         call bummer('code4x: k1s=', k1s, faterr)
      elseif ( k1f .gt. maxbfn .or. k1f .lt. 0 ) then
         call bummer('code4x: k1f=', k1f, faterr)
      elseif ( k1l1 .gt. maxsiz .or. k1l1 .lt. 0 ) then
         call bummer('code4x: k1l1=', k1l1, faterr)
      endif
c
      ib1hi  = ishft(strtnf,-2)
      ib1low = iand(strtnf,mask2)
      ib2hi  = ishft(stn121,-2)
      ib2low = iand(stn121,mask2)
      ib3hi  = ishft(stn3,-2)
      ib3low = iand(stn3,mask2)
      ib4hi  = ishft(nintkl,-2)
      ib4low = iand(nintkl,mask2)
c
      iword(1)  = ior(ishft(ib1hi,16), ib2hi)
      iword(2)  = ior(ishft(ib3hi,16), ib4hi)
      buffer(1) = dword
cold
c     iword(1)  = ior(ishft(l1s,  24), ior(ishft(l1f,16),
c    &            ior(ishft(k1s,    8),          k1f)))
      iword(1)  = 0
      iword(2)  = ior(ishft(lsym-1,29), ior(ishft(ksym-1,26),
     &            ior(ishft(ib1low,24), ior(ishft(ib2low,22),
     &            ior(ishft(ib3low,20), ior(ishft(ib4low,18),
     &            k1l1))))))
      buffer(2) = dword
cnew
       iword(1)  = ior(ishft(l1s,  16), l1f)
       iword(2)  = ior(ishft(k1s,  16), k1f)
       buffer(3) = dword
c      write(0,*) 'code4x:',l1s,l1f,k1s,k1f
*@endif
      return
      end



ctm dcod3x stems from ciudg7.f

      subroutine dcod3x(
     & rdbuf,  iblsym, ksym,   k1s,
     & k1f,    l1bl,   intblk, stn121,
     & stn3,   stn1,   stn2,   nintkl )
c
c  decodes integral labels prepared by cisrt.
c  see subroutine () for the packing convention.
c
c  03-jul-01 modified to allow large basis sets (>170). -rls
c
      implicit logical(a-z)
c     # dummy:
      integer iblsym, ksym, k1s, k1f, l1bl, intblk, stn121, stn3,
     & stn1, stn2, nintkl
       integer i,j
       logical first
       data first /.false./
       save first
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
*@ifdef int64
*       integer m8,m16
*       parameter(m8=2**8-1, m16=2**16-1)
*        integer wkword
*        integer rdbuf(*)
*        integer und,rechts
*@ifdef cray T3E64 T3D
*        und(i,j) =and(i,j)
*        rechts(i,j) = shiftr(i,j)
*@else
*        und(i,j) = iand(i,j)
*        rechts(i,j) = ishft(i,-j)
*@endif
*       if (first) then
*       write(0,*) 'more than 255 basis functions not supported in
*     & 64bit'
*       first=.false.
*       endif
*@ifdef cray or T3E64
*      wkword  =  rdbuf(2)
*      iblsym  =  und( wkword, mask(120) )
*      ksym  =  und( rechts(wkword, 40 ),mask(120) )
*      k1s   =  und( rechts( wkword, 32), mask(120) )
*      k1f   =  und( rechts(wkword, 24),mask(120) )
*      l1bl  =  und( rechts(wkword, 8),mask(120) )
*c
*      wkword  =  rdbuf(1)
*      intblk  =  und( rechts(wkword, 48 ), mask(112) )
*      stn121  =  und( rechts( wkword, 32), mask(112) )
*      stn3    =  und( wkword, mask(112) )
*c
*      if ( iblsym.ne.1 ) then
*          stn1    =  stn121
*          stn2    =  und( rechts(wkword, 16), mask(112) )
*          wkword  =  rdbuf(2)
*          nintkl  =  und( rechts(wkword, 48 ), mask(112) )
*      endif
*@else
*      wkword  =  rdbuf(2)
*      iblsym  =   und ( wkword, m8       )
*      ksym  =    und( rechts(wkword, 40 ),m8        )
*      k1s   =    und( rechts( wkword, 32), m8        )
*      k1f   =    und( rechts(wkword, 24),m8        )
*      l1bl  =    und( rechts(wkword,8), m8       )
*c
*      wkword  =  rdbuf(1)
*      intblk  =    und( rechts(wkword,48 ), m16       )
*      stn121  =    und( rechts( wkword, 32), m16      )
*      stn3    =    und( wkword, m16  )
*c
*      if ( iblsym.ne.1 ) then
*          stn1    =  stn121
*          stn2    =    und( rechts(wkword, 16), m16       )
*          wkword  =  rdbuf(2)
*          nintkl  =    und( rechts(wkword, 48 ), m16       )
*      endif
*@endif
*@elif defined  (oldoldsun) || defined ( fps)
*c     # local:
*      real*8 rdbuf(*)
*      real*8 wkword
*      integer iword(2),iword1,iword2
*      equivalence (iword(1),iword1,wkword), (iword(2),iword2)
*      integer   mask8,        mask16
*      parameter(mask8=2**8-1, mask16=2**16-1)
*c
*      integer   and, rshift
*      intrinsic and, rshift
*c
*      wkword  =  rdbuf(2)
*      ksym  =  and( rshift(iword1, 8), mask8)
*      k1s   =  and(        iword1,     mask8)
*      k1f   =  and( rshift(iword2,24), mask8)
*      l1bl  =  and( rshift(iword2, 8), mask8)
*      iblsym=  and(        iword2,     mask8)
*c
*      wkword  =  rdbuf(1)
*      intblk  =  and( rshift(iword1,16), mask16)
*      stn121  =  and(        iword1,     mask16)
*      stn3    =  and(        iword2,     mask16)
*c
*      if ( iblsym.ne.1 ) then
*          stn1   = stn121
*          stn2   = and( rshift(iword2,16), mask16)
*          wkword = rdbuf(2)
*          nintkl = and( rshift(iword1,16), mask16)
*      endif
*@elif defined  obsolete
*c  32-bit integer machines with f90 bit operators.
*c
*c  note: the 3-way equivalence version breaks the titan optimizer.
*c        use 2-way equivalences and array references instead. -rls
*c
*c     # local:
*      real*8 rdbuf(*)
*      real*8 wkword
*      integer iword(2)
*      equivalence (iword(1),wkword)
*      integer   mask8,        mask16
*      parameter(mask8=2**8-1, mask16=2**16-1)
*c
*      integer   iand, ishft
*      intrinsic iand, ishft
*c
*      wkword  =  rdbuf(2)
*      ksym    =  iand( ishft(iword(1), -8), mask8)
*      k1s     =  iand(       iword(1),      mask8)
*      k1f     =  iand( ishft(iword(2),-24), mask8)
*      l1bl    =  iand( ishft(iword(2), -8), mask8)
*      iblsym  =  iand(       iword(2),      mask8)
*c
*      wkword  =  rdbuf(1)
*      intblk  =  iand( ishft(iword(1),-16), mask16)
*
*      stn121  =  iand(       iword(1),      mask16)
*      stn3    =  iand(       iword(2),      mask16)
*c
*      if ( iblsym .ne. 1 ) then
*          stn1   = stn121
*          stn2   = iand( ishft(iword(2),-16), mask16)
*          wkword = rdbuf(2)
*          nintkl = iand( ishft(iword(1),-16), mask16)
*      endif
*@else
c     # 32-bit integer machines with f90 bit operators.
c     # 18-bit limits on the buffer lengths.
c
c  note: the 3-way equivalence version breaks the titan optimizer.
c        use 2-way equivalences and array references instead. -rls
c
      real*8 rdbuf(*)
c
c     # local:
      real*8 wkword
      integer iword(2)
      equivalence (iword(1),wkword)
c
      integer ib1hi,  ib2hi,  ib3hi,  ib4hi,  ib5hi
      integer ib1low, ib2low, ib3low, ib4low, ib5low
c
      integer   mask2,   mask4,    mask8,        mask16
      parameter(mask2=3, mask4=15, mask8=2**8-1, mask16=2**16-1)

      integer   mask10,mask12,mask6
      parameter(mask10=2**10-1, mask12=2**12-1,mask6=2**6-1)
c
      integer   iand, ishft
      intrinsic iand, ishft
c
      wkword  =  rdbuf(2)
      ib5hi   =        ishft(iword(1),-16)
      k1s     =  iand( ishft(iword(1), -6), mask10)
      k1f     =  iand(       iword(1),      mask6)
      k1f     =  ior(ishft(k1f,4),ishft(iword(2),-28))
      l1bl    =  iand( ishft(iword(2),-16), mask12)
      ksym    =  iand( ishft(iword(2),-12), mask4)
      ib1low  =  iand( ishft(iword(2),-10), mask2)
      ib2low  =  iand( ishft(iword(2), -8), mask2)
      ib3low  =  iand( ishft(iword(2), -6), mask2)
      ib4low  =  iand( ishft(iword(2), -4), mask2)
      ib5low  =  iand( ishft(iword(2), -2), mask2)
      iblsym  =  iand(       iword(2),      mask2)
c
      wkword  =  rdbuf(1)
      ib1hi   =        ishft(iword(1),-16)
      ib2hi   =  iand(       iword(1),      mask16)
      ib3hi   =        ishft(iword(2),-16)
      ib4hi   =  iand(       iword(2),      mask16)

      intblk  =  ior(ishft(ib1hi,2),ib1low)
      stn121  =  ior(ishft(ib2hi,2),ib2low)
      stn3    =  ior(ishft(ib4hi,2),ib4low)
c
      if ( iblsym .eq. 0 ) then
          stn1   = stn121
          stn2   = ior(ishft(ib3hi,2),ib3low)
          nintkl = ior(ishft(ib5hi,2),ib5low)
      endif
*@endif
      return
      end

ctm dcod4x stems from ciudg8.f

      subroutine dcod4x(
     & rdbuf,  intblk, stn121, stn3,
     & nintkl, lsym,   ksym,   l1s,
     & l1f,    k1s,    k1f,    k1l1 )
c
c  decodes integral labels prepared by cisrt
c  see subroutine code4x() for the packing convention.
c
c  03-jul-01 modified to allow large basis sets (>170). -rls
c
      implicit logical(a-z)
c     # dummy:
      integer intblk, stn121, stn3, nintkl, lsym, ksym, l1s, l1f,
     & k1s, k1f, k1l1
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
       logical first
       data first /.false./
       save first
c
*@ifdef int64
*       integer m8,m16
*       parameter (m8=2**8-1, m16=2**16-1)
*       integer wkword
*       integer rdbuf(*)
*       integer i,j,rechts, und
*@ifdef cray or T3E64
*       rechts(i,j) = shiftr(i,j)
*       und(i,j) = and(i,j)
*@else
*       rechts(i,j) =ishft(i,-j)
*       und(i,j)    =and(i,j)
*@endif
*       if (first) then
*       write(0,*) 'more than 255 basis functions not supported in
*     &  64bit'
*       first=.false.
*       endif
*@ifdef cray or T3E64
*      integer   and, shiftr
*      intrinsic and, shiftr
*      wkword  =  rdbuf(1)
*      intblk  =  and( shiftr(wkword,48), mask(112) )
*      stn121  =  and( shiftr(wkword,32), mask(112) )
*      stn3    =  and( shiftr(wkword,16), mask(112) )
*      nintkl  =  and( wkword, mask(112) )
*c
*      wkword  =  rdbuf(2)
*      lsym  =  and( shiftr(wkword,56), mask(120) )
*      ksym  =  and( shiftr(wkword,48), mask(120) )
*      l1s   =  and( shiftr(wkword,40), mask(120) )
*      l1f   =  and( shiftr(wkword,32), mask(120))
*      k1s   =  and( shiftr(wkword,24), mask(120))
*      k1f   =  and( shiftr(wkword,16), mask(120))
*      k1l1  =  and( wkword, mask(112))
*@else
*      wkword  =  rdbuf(1)
*      intblk  =    und( rechts(wkword,48), m16)
*      stn121  =    und( rechts(wkword,32), m16       )
*      stn3    =    und( rechts(wkword,16), m16       )
*      nintkl  =    und( wkword, m16       )
*c
*      wkword  =  rdbuf(2)
*      lsym  =    und( rechts(wkword,56), m8        )
*      ksym  =    und( rechts(wkword,48), m8        )
*      l1s   =    und( rechts(wkword,40), m8        )
*      l1f   =    und( rechts(wkword,32), m8       )
*      k1s   =    und( rechts(wkword,24), m8       )
*      k1f   =    und( rechts(wkword,16), m8       )
*      k1l1  =    und( wkword, m8)
*@endif
*@elif defined  (oldoldsun) || defined ( fps)
*c     # local:
*      real*8 wkword
*      real*8 rdbuf(*)
*      integer iword(2),iword1,iword2
*      equivalence (iword(1),iword1,wkword)
*      equivalence (iword(2),iword2)
*      integer mask4,mask8,mask15,mask16
*      parameter(mask4  = 2**4-1,  mask8  = 2**8-1)
*      parameter(mask15 = 2**15-1, mask16 = 2**16-1)
*c
*      integer   and, rshift
*      intrinsic and, rshift
*c
*      wkword  =  rdbuf(1)
*      intblk  =  and( rshift(iword1,16), mask15)
*      stn121  =  and(        iword1,     mask16)
*      stn3    =  and( rshift(iword2,16), mask16)
*      nintkl  =  and(        iword2,     mask16)
*      wkword  =  rdbuf(2)
*      lsym    =  and( rshift(iword1,24), mask4 )
*      ksym    =  and( rshift(iword1,16), mask8 )
*      l1s     =  and( rshift(iword1, 8), mask8 )
*      l1f     =  and(        iword1,     mask8 )
*      k1s     =  and( rshift(iword2,24), mask8 )
*      k1f     =  and( rshift(iword2,16), mask8 )
*      k1l1    =  and(        iword2,     mask16)
*@elif defined  obsolete
*c  32-bit integer machines with f90 bit operators.
*c
*c  note: the 3-way equivalence version breaks the titan optimizer.
*c        use 2-way equivalences and array references instead. -rls
*c
*c     # local:
*      real*8 wkword
*      real*8 rdbuf(*)
*      integer iword(2)
*      equivalence (iword(1),wkword)
*      integer mask8,mask16
*      parameter(mask8  = 2**8-1, mask16 = 2**16-1)
*c
*      integer   iand, ishft, ior
*      intrinsic iand, ishft, ior
*c
*      wkword  =  rdbuf(1)
*      intblk  =  iand( ishft(iword(1),-16), mask16)
*      stn121  =  iand(       iword(1),      mask16)
*      stn3    =  iand( ishft(iword(2),-16), mask16)
*      nintkl  =  iand(       iword(2),      mask16)
*c
*      wkword  =  rdbuf(2)
*      lsym    =  iand( ishft(iword(1),-24), mask8 )
*      ksym    =  iand( ishft(iword(1),-16), mask8 )
*      l1s     =  iand( ishft(iword(1), -8), mask8 )
*      l1f     =  iand(       iword(1),      mask8 )
*      k1s     =  iand( ishft(iword(2),-24), mask8 )
*      k1f     =  iand( ishft(iword(2),-16), mask8 )
*      k1l1    =  iand(       iword(2),      mask16)
*@else
c     # 32-bit integer machines with f90 bit operators.
c     # 18-bit limits on the buffer lengths.
c
c  note: the 3-way equivalence version breaks the titan optimizer.
c        use 2-way equivalences and array references instead. -rls
c
c     # dummy:
      real*8 rdbuf(*)
c
c     # local:
      real*8 wkword
      integer iword(2)
      equivalence (iword(1),wkword)
c
      integer ib1hi,  ib2hi,  ib3hi,  ib4hi,  ib5hi
      integer ib1low, ib2low, ib3low, ib4low, ib5low
c
      integer   mask2,   mask3,   mask8,        mask16
      parameter(mask2=3, mask3=7, mask8=2**8-1, mask16=2**16-1)
      integer   mask18
      parameter(mask18=2**18-1)
c
      integer   iand, ishft, ior
      intrinsic iand, ishft, ior
c
      wkword  =  rdbuf(1)
      ib1hi   =        ishft(iword(1),-16)
      ib2hi   =  iand(       iword(1),      mask16)
      ib3hi   =        ishft(iword(2),-16)
      ib4hi   =  iand(       iword(2),      mask16)
      intblk  =  iand( ishft(iword(1),-16), mask16)
      stn121  =  iand(       iword(1),      mask16)
      stn3    =  iand( ishft(iword(2),-16), mask16)
      nintkl  =  iand(       iword(2),      mask16)
c
      wkword  =  rdbuf(2)
cold
c     l1s     =        ishft(iword(1),-24)
c     l1f     =  iand( ishft(iword(1),-16), mask8 )
c     k1s     =  iand( ishft(iword(1), -8), mask8 )
c     k1f     =  iand(       iword(1),      mask8 )
      lsym    =        ishft(iword(2),-29) + 1
      ksym    =  iand( ishft(iword(2),-26), mask3 ) + 1
      ib1low  =  iand( ishft(iword(2),-24), mask2 )
      ib2low  =  iand( ishft(iword(2),-22), mask2 )
      ib3low  =  iand( ishft(iword(2),-20), mask2 )
      ib4low  =  iand( ishft(iword(2),-18), mask2 )
      k1l1    =  iand(       iword(2),      mask18)
c
      intblk  = ior(ishft(ib1hi,2),ib1low)
      stn121  = ior(ishft(ib2hi,2),ib2low)
      stn3    = ior(ishft(ib3hi,2),ib3low)
      nintkl  = ior(ishft(ib4hi,2),ib4low)
c     return
cnew

      wkword  =  rdbuf(3)
      l1s     =        ishft(iword(1),-16)
      l1f     =  iand(  iword(1), mask16 )
      k1s     =        ishft(iword(2),-16)
      k1f     =  iand(  iword(2), mask16 )
c      write(0,*) 'dcod4x:',l1s,l1f,k1s,k1f
*@endif
      return
      end

      subroutine encodf(
     & nv,     nsv,    symw,   head,
     & tail,   icode,  val,    isv,
     & yb,     yk,     nw,     wword,
     & isv2,   imode )
c**********************************************************************
c
c  pack formula info for a loop into a structure.
c
c  input:
c  nv = number of loop values.
c  nsv = number of symmetry versions of the loop.
c  symw = loop symmetry.
c  head = head within the drt.
c  tail = tail within the drt.
c  icode = loop code.
c  val(1:nv) = loop values.
c  isv(1:nsv) = contributing symmetry versions (bra).
c  isv(1:nsv) = contributing symmetry versions (ket).
c  imode = 0: isv(bra) = isv(ket) (mcscf)
c  imode = 1: isv(bra) <> isv(ket) (transft,transmom)
c  yb(1:nsv) = bra loop weights for each symmetry.
c  yk(1:nsv) = ket loop weights for each symmetry.
c
c  output:
c  nw = number of packed entries.
c  wword(1:nw) = working precision buffer of packed entries.
c
      implicit logical(a-z)
c     # dummy:
      integer nv, nsv, symw, head, tail, icode, nw
      integer isv(8), yb(8), yk(8) , isv2(8), imode
C
c     # local:
      integer i, j, i1, i2, i3, i4, i5, i6, i7, i8, mask, ipt
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # set up any required bit operators...
c
*@ifdef int64
*      integer val(3), wword(8), wwords
*      integer oder,und,links,rechts
*      integer iword
*      integer  dword
*      equivalence (dword,iword)
*c
*      integer   m4,        m8,        m16
*      parameter(m4=2**4-1, m8=2**8-1, m16=2**16-1)
*      integer  ip4, ip8, ip16,  ipack8, ipack6, sand
*@ifdef cray T3E64
*       oder(i,j)=or(i,j)
*       und(i,j)=and(i,j)
*       links(i,j) =shiftl(i,j)
*       rechts(i,j) = shiftr(i,j)
*@else
*       oder(i,j)=ior(i,j)
*       und(i,j)=iand(i,j)
*       links(i,j) =ishft(i,j)
*       rechts(i,j) = ishft(i,-j)
*@endif
*c
*      ip4(i,j)  = oder(links(i,4),  und(j,m4))
*      ip8(i,j)  = oder(links(i,8),  und(j,m8))
*      ip16(i,j) = oder(links(i,16),  und(j,m16))
*      ipack8(i1,i2,i3,i4,i5,i6,i7,i8) = ip16(
*     & ip8( ip4(i1,i2), ip4(i3,i4) ),
*     & ip8( ip4(i5,i6), ip4(i7,i8) )  )
*      ipack6(i1,i2,i3,i4,i5,i6) = ip16(
*     & ip8( ip4(i1,i2),i3), ip8(i4,ip4(i5,i6)) )
*      sand(i,j,mask) = und( rechts(i,j), mask )
*
*@elif defined  oldsun
*       real*8 val(3), wword(8), wwords
*       integer   m4,        m8,        m16
*       parameter(m4=2**4-1, m8=2**8-1, m16=2**16-1)
*       integer iword(2)
*       real*8 dword
*       equivalence (iword(1),dword)
*c
*       integer ip4, ip8, ip16, ipack8, ipack6, sand
*       ip4(i,j)  = or(ishft(i,4),and(j,m4))
*       ip8(i,j)  = or(ishft(i,8),and(j,m8))
*       ip16(i,j) = or(ishft(i,16),and(j,m16))
*       ipack8(i1,i2,i3,i4,i5,i6,i7,i8) = ip16(
*     & ip8( ip4(i1,i2), ip4(i3,i4) ),
*     & ip8( ip4(i5,i6), ip4(i7,i8) )  )
*       ipack6(i1,i2,i3,i4,i5,i6) = ip16(
*     & ip8( ip4(i1,i2),i3), ip8(i4,ip4(i5,i6)) )
*       sand(i,j,mask) = and( ishft(i,-j), mask )
*@elif defined ( milstd1753) || defined ( titan) || defined ( alliant) || defined ( f90 )
      real*8 val(3), wword(8), wwords
c     # 32-bit integer machines using f90 bit operators.
      integer   m4,        m8,        m16
      parameter(m4=2**4-1, m8=2**8-1, m16=2**16-1)
      integer iword(2)
      real*8 dword
      equivalence (iword(1),dword)
c
      integer   ior, iand, ishft
      intrinsic ior, iand, ishft
c
      integer  ip4, ip8, ip16, ipack8, ipack6, sand
      ip4(i,j)  = ior(ishft(i,4),iand(j,m4))
      ip8(i,j)  = ior(ishft(i,8),iand(j,m8))
      ip16(i,j) = ior(ishft(i,16),iand(j,m16))
      ipack8(i1,i2,i3,i4,i5,i6,i7,i8) = ip16(
     & ip8( ip4(i1,i2), ip4(i3,i4) ),
     & ip8( ip4(i5,i6), ip4(i7,i8) )  )
      ipack6(i1,i2,i3,i4,i5,i6) = ip16(
     & ip8( ip4(i1,i2),i3), ip8(i4,ip4(i5,i6)) )
      sand(i,j,mask) = iand( ishft(i,-j), mask )
*@else
*      real*8 val(3), wword(8), wwords
*c     # machine-independent version.
*c     # packing/unpacking uses arithmetic and external library calls.
*      integer ibit8, nsvx
*      integer unpack(16), ibit(8)
*      integer iword(2)
*      real*8 wlocal(1)
*c
*      data ibit / 1, 2, 4, 8, 16, 32, 64, 128 /
*@endif

c  error check

        if (imode.ne.0 .and. imode.ne.1)
     .      call bummer('encodf invalid imode, imode=',imode,2)
c


c
c     # pack the formula entries...
c
*@ifdef int64
*      nw = (nv+1) + (nsv+1)/2 + imode
*      wword(1) = oder(links(
*     & ipack8(isv(1),isv(2),isv(3),isv(4),
*     & isv(5),isv(6),isv(7),isv(8)),32),
*     & ipack6(nsv,symw,head,tail,nw,icode))
*
*      if (imode.eq.1) then
*      wword(2) = ipack8(isv2(1),isv2(2),isv2(3),isv2(4),
*     .                  isv2(5),isv2(6),isv2(7),isv2(8))
*      endif
*
*c
*      do 10 i = 1, nv
*         wword(i+1+imode) = val(i)
*10    continue
*c
*      ipt = 1
*      do 20 i = 1, ((nsv+1)/2)
*         wword(1+nv+i+imode) = oder(links(
*     &    ip16(yb(isv(ipt  )),yk(isv2(ipt  ))), 32),
*     &    ip16(yb(isv(ipt+1)),yk(isv2(ipt+1)))      )
*         ipt = ipt + 2
*20    continue
*@elif defined ( milstd1753) || defined ( titan) || defined ( alliant) || defined ( f90 )
c     # 32-bit integer machines using f90 bit-operators.
c     # wword(1)             : (isv(i),i=1,nsv),
c     #                        nsv,symw,head,tail,nw,icode
c     # wword(2)-wword(1+nv) : (val(i),i=1,nv)
c     # wword(2+nv)...       : (yb(i),yk(i),i=1,nsv)
c     # resulting in 3 to 8 words used for each loop.
c
c     real*8 val(3), wword(8), wwords
c     real*8 dword
c     equivalence (iword(1),dword)

      nw = (nv+1) + (nsv+1)/2+imode
      iword(1) = ipack8(
     & isv(1), isv(2), isv(3), isv(4),
     & isv(5), isv(6), isv(7), isv(8) )
      iword(2) = ipack6( nsv, symw, head, tail, nw, icode )
      wword(1) = dword
      if (imode.eq.1) then
      iword(1) = ipack8(
     & isv2(1), isv2(2), isv2(3), isv2(4),
     & isv2(5), isv2(6), isv2(7), isv2(8) )
      iword(2) = 0
      wword(2) = dword
      endif
c
      do 10 i = 1, nv
         wword(i+1+imode) = val(i)
10    continue
c
      ipt = 1
      do 20 i = 1, ((nsv+1)/2)
         iword(1) = ip16( yb(isv(ipt  )), yk(isv2(ipt  )) )
         iword(2) = ip16( yb(isv(ipt+1)), yk(isv2(ipt+1)) )
         wword(1+nv+i+imode) = dword
         ipt = ipt + 2
20    continue
*@else
*c     # machine-independent version.
*c     # wword(1)             : ibit8,nsv,symw,head,tail,nw,icode
*c     # wword(2)-wword(1+nv) : (val(i),i=1,nv)
*c     # wword(2+nv)...       : (yb(i),yk(i),i=1,nsv)
*c     # resulting in 3 to 8 words used for each loop.
*c
*      nw = nv + 1 + (nsv+1)/2 + imode
*c
*c     # bit-pack the isv(*) vector.
*c     # isv(*) entries are monotonic.
*c     # this requires at most nsym (<=8) bits.
*      ibit8 = 0
*      do 10 i = 1, nsv
*         ibit8 = ibit8 + ibit( isv(i) )
*10    continue
*
*      unpack(1) = ibit8 * 256 + nv * 16 + nsv
*      unpack(2) = head * 256 + tail
*      unpack(3) = nw * 16 + symw
*      unpack(4) = icode
*      call plab16( wword(1), unpack, 4 )
*c
*      if (imode.eq.1) then
*      ibit8=0
*      do 11 i=1,nsv
*         ibit8= ibit8+ibit(isv2(i))
*11    continue
*      unpack(1)=ibit8*256
*      call plab16(wword(2),unpack,4)
*      endif
*
*      do 20 i = 1, nv
*         wword(i+1+imode) = val(i)
*20    continue
*c
*      do 30 i = 1, nsv
*         unpack( 2 * (i-1) + 1 ) = yb( isv(i) )
*         unpack( 2 * (i-1) + 2 ) = yk( isv2(i) )
*30    continue
*      call plab16( wword(nv+2+imode), unpack, 2*nsv )
*c
*@endif
c
c
       return
c**********************************************************************
c deck pcodef
      entry pcodef( icode, nw, wwords )
c**********************************************************************
c
c     # pack a code value into the formula structure.
c     # this is used for icode=0 and icode=1.
c
*@ifdef int64
*      nw     = 1
*      IWORD  = ip4(nw,icode)
*      wwords = dword
*@elif defined ( milstd1753) || defined (  sun) || defined ( f90)
c     # 32-bit machines using f90 bit operators.
       nw       = 1
       iword(2) = ip4(nw,icode)
       wwords   = dword
*@else
*c     # machine-independent version.
*      nw = 1
*      unpack(1) = 0
*      unpack(2) = 0
*      unpack(3) = 1 * 16 + 0
*      unpack(4) = icode
*      call plab16( wlocal, unpack, 4 )
*      wwords = wlocal(1)
*@endif
      return
c
c**********************************************************************
c deck decodf
      entry decodf(
     & nv,     nsv,    symw,   head,
     & tail,   icode,  val,    isv,
     & yb,     yk,     nw,     wword ,
     & isv2, imode )
c**********************************************************************
c
c     # unpack a packed formula structure.
c     #
c     # input:
c     # wword(*) = packed structure.
c     #
c     # output:
c     # icode = loop code.
c     # nw = number of packed words.
c     # nv,nsv,symw,head,tail,val(*),isv(*),yb(*),yk(*) = formula
c     #                           entries. (computed only if icode>1 )
c
*@ifdef int64
*      dword = wword(1)
*      icode = sand(iword,0,m4)
*      nw    = sand(iword,4,m4)
*      if ( icode .gt. 1 ) then
*         tail   = sand( iword,  8, m8 )
*         head   = sand( iword, 16, m8 )
*         symw   = sand( iword, 24, m4 )
*         nsv    = sand( iword, 28, m4 )
*         isv(1) = sand( iword, 60, m4 )
*         isv(2) = sand( iword, 56, m4 )
*         isv(3) = sand( iword, 52, m4 )
*         isv(4) = sand( iword, 48, m4 )
*         isv(5) = sand( iword, 44, m4 )
*         isv(6) = sand( iword, 40, m4 )
*         isv(7) = sand( iword, 36, m4 )
*         isv(8) = sand( iword, 32, m4 )
*c
*         if (imode.eq.1) then
*         dword = wword(2)
*         isv2(1) = sand( iword, 28, m4 )
*         isv2(2) = sand( iword, 24, m4 )
*         isv2(3) = sand( iword, 20, m4 )
*         isv2(4) = sand( iword, 16, m4 )
*         isv2(5) = sand( iword, 12, m4 )
*         isv2(6) = sand( iword,  8, m4 )
*         isv2(7) = sand( iword,  4, m4 )
*         isv2(8) = sand( iword,  2, m4 )
*         endif
*
*
*         nv = nw - 1 - (nsv + 1) / 2 -imode
*         do 210 i = 1, nv
*            val(i) = wword(i+1+imode)
*210      continue
*         ipt = 1
*CDIR$ NEXTSCALAR
*         do 220 i = 1, (nsv+1)/2
*            dword     = wword(1+nv+i+imode)
*            yb(ipt)   = sand( iword, 48, m16 )
*            yk(ipt)   = sand( iword, 32, m16 )
*            yb(ipt+1) = sand( iword, 16, m16 )
*            yk(ipt+1) = sand( iword,  0, m16 )
*            ipt       = ipt + 2
*220      continue
*      endif
*@elif defined ( milstd1753) || defined ( titan) || defined ( alliant) || defined ( f90 )
c     # 32-bit integer machines using f90 bit operators.
      dword = wword(1)
      icode = sand(iword(2),0,m4)
      nw    = sand(iword(2),4,m4)
      if ( icode .gt. 1 ) then
         tail   = sand( iword(2),  8, m8 )
         head   = sand( iword(2), 16, m8 )
         symw   = sand( iword(2), 24, m4 )
         nsv    = sand( iword(2), 28, m4 )
         isv(1) = sand( iword(1), 28, m4 )
         isv(2) = sand( iword(1), 24, m4 )
         isv(3) = sand( iword(1), 20, m4 )
         isv(4) = sand( iword(1), 16, m4 )
         isv(5) = sand( iword(1), 12, m4 )
         isv(6) = sand( iword(1),  8, m4 )
         isv(7) = sand( iword(1),  4, m4 )
         isv(8) = sand( iword(1),  0, m4 )
c
         if (imode.eq.1) then
           dword = wword(2)
           isv2(1) = sand( iword(1), 28, m4 )
           isv2(2) = sand( iword(1), 24, m4 )
           isv2(3) = sand( iword(1), 20, m4 )
           isv2(4) = sand( iword(1), 16, m4 )
           isv2(5) = sand( iword(1), 12, m4 )
           isv2(6) = sand( iword(1),  8, m4 )
           isv2(7) = sand( iword(1),  4, m4 )
           isv2(8) = sand( iword(1),  0, m4 )
         endif

         nv = nw - 1 - (nsv + 1) / 2 -imode
         do 210 i = 1, nv
            val(i) = wword(i+1+imode)
210      continue
         ipt = 1
         do 220 i = 1, (nsv+1)/2
            dword     = wword(1+nv+i+imode)
            yb(ipt)   = sand( iword(1), 16, m16 )
            yk(ipt)   = sand( iword(1),  0, m16 )
            yb(ipt+1) = sand( iword(2), 16, m16 )
            yk(ipt+1) = sand( iword(2),  0, m16 )
            ipt       = ipt + 2
220      continue
      endif
*@else
*c     # machine-independent version.
*c
*      call ulab16( wword(1), unpack, 4 )
*      icode = unpack(4)
*      nw    = unpack(3) / 16
*      if ( icode .gt. 1 ) then
*         ibit8 =      unpack(1) / 256
*         nv    = mod( unpack(1) / 16,   16 )
*         nsv   = mod( unpack(1),        16 )
*         head  =      unpack(2) / 256
*         tail  = mod( unpack(2),       256 )
*         symw  = mod( unpack(3),        16 )
*c
*c        # unpack isv(*).  bit-packed method is used in this version as
*c        # an example in case this approach is needed for future work.
*         nsvx = 0
*         do 210 i = 1, 8
*            if ( mod( ibit8, 2 ) .eq. 1 ) then
*               nsvx      = nsvx + 1
*               isv(nsvx) = i
*            endif
*            ibit8 = ibit8 / 2
*210      continue
*         if ( nsvx .ne. nsv ) then
*            call bummer('decodf: (nsvx-nsv)=', (nsvx-nsv), faterr )
*         endif
*
*         if (imode.eq.1) then
*          call ulab16( wword(2), unpack, 4 )
*          ibit8 =      unpack(1) / 256
*         nsvx = 0
*         do 211 i = 1, 8
*            if ( mod( ibit8, 2 ) .eq. 1 ) then
*               nsvx      = nsvx + 1
*               isv2(nsvx) = i
*            endif
*            ibit8 = ibit8 / 2
*211      continue
*         if ( nsvx .ne. nsv ) then
*            call bummer('decodf: (nsvx-nsv)=', (nsvx-nsv), faterr )
*         endif
*         endif
*
*
*c
*         do 220 i = 1, nv
*            val(i) = wword(i+1+imode)
*220      continue
*c
*         call ulab16( wword(nv+2), unpack, 2*nsv )
*         do 230 i = 1, nsv
*            yb(i) = unpack( 2 * (i-1) + 1 )
*            yk(i) = unpack( 2 * (i-1) + 2 )
*230      continue
*      endif
*@endif
      return
      end

      function getstepventry(nbfn,occ,bfn,ref)

      integer getstepventry,nbfn,occ(*),bfn,ref
      integer ibit,iword,itmp
      integer wordlen
*@ifdef int64
*        parameter(wordlen=64)
*@else
      parameter(wordlen=32)
*@endif


      ibit=(ref-1)*nbfn*2 + (bfn)*2
      iword=(ibit-1)/wordlen
c     write(*,*) 'getstv: getting ibit,iword=',ibit,iword+1
      ibit = ibit-iword*wordlen
      iword=occ(iword+1)

c       extract the two bits no ibit and ibit-1 from iword


      itmp=iand(3,ishft(iword,-ibit+2))
c   write(*,*) 'getstepventry: extracting bfn,ref',bfn,ref,'=',itmp
      getstepventry=itmp

      return
      end

ctm in2 stems from ciden9.f


      function in2(i)
c     this function gives the number of full real words for a given
c     number of integer words
*@ifdef int64
*      in2 = i
*@else
       in2=(i-1)/2+1
*@endif
      return
      end

ctm in stems from ciden9.f

      function in(i)
c     this function ensures that for ibm machines an integer array
c     ends at a double word boundary
*@ifdef int64
*      in = i
*@else
       in=((i-1)/2)*2+2
*@endif
      return
      end

ctm labget stems from ciden2.f

      subroutine labget(word,lab,ipos)
c
c  this routine gets the integer label 'lab' from the 'ipos' position
c  of the working precision word 'word'.
c
       implicit integer(a-z)
*@ifdef int64
*      integer iword(2),lab,ipos
*      integer m32
*      parameter(m32=2**32-1)
*@ifdef T3E64 cray
*c
*      iword(1)=shiftr(word,32)
*      iword(2)=and(word,m32)
*      lab=iword(ipos)
*@else
*       iword(1)=ishft(word,-32)
*       iword(2)=iand(word,m32)
*       lab=iword(ipos)
*@endif
*@else
      integer word(2)
      lab=word(ipos)
*@endif
      return
      end

ctm labgt1 stems from ciden2.f

      subroutine labgt1(word,lab,ipos)
c
c  this routine gets the integer label 'lab' from the 'ipos' position
c  of the working precision word 'word'.
c
      implicit integer(a-z)
*@ifdef int64
*      integer iw16, word
*      parameter(m16=2**16-1)
*@ifdef T3E64 cray
*      iw16=and(shiftr(word,16*(4-ipos)),m16)
*      lab=iw16
*@else
*       iw16=iand(ishft(word,16*(ipos-4)),m16)
*       lab=iw16
*@endif
*@else
       real*8 word
       real*8 dpword
      integer*2 iw16(4)
      equivalence (iw16(1),dpword)

      dpword=word
      lab=iw16(ipos)
*@endif
      return
      end

ctm labpt1 stems from ciden9.f

      subroutine labpt1(word,lab,ipos)
c
c  this routine puts the integer label 'lab' into the 'ipos' position
c  of the working precision word 'word'.
c
      implicit integer(a-z)
*@ifdef int64
*      integer iw16(4),word
*      integer m16
*      parameter(m16=2**16-1)
*@ifdef T3E64 cray
*      integer dpack4
*      dpack4(i1,i2,i3,i4)=
*     +  or(shiftl(or(shiftl(i1,16),i2),32),or(shiftl(i3,16),i4))
*c
*      iw16(1)=shiftr(word,48)
*      iw16(2)=and(shiftr(word,32),m16)
*      iw16(3)=and(shiftr(word,16),m16)
*      iw16(4)=and(word,m16)
*      iw16(ipos)=and(lab,m16)
*      word=dpack4(iw16(1),iw16(2),iw16(3),iw16(4))
*@else
*      integer dpack4
*      dpack4(i1,i2,i3,i4)=
*     +  ior(ishft(ior(ishft(i1,16),i2),32),ior(ishft(i3,16),i4))
*c
*      iw16(1)=ishft(word,-48)
*      iw16(2)=iand(ishft(word,-32),m16)
*      iw16(3)=iand(ishft(word,-16),m16)
*      iw16(4)=iand(word,m16)
*      iw16(ipos)=iand(lab,m16)
*      word=dpack4(iw16(1),iw16(2),iw16(3),iw16(4))
*@endif
*@else
      real*8 word
      integer*2 iw16(4)
      real*8 dpword
      equivalence (iw16(1),dpword)
      dpword=word
      iw16(ipos)=lab
      word=dpword
*@endif
      return
      end

ctm labput stems from ciden9.f

      subroutine labput(word,lab,ipos)
c
c  this routine puts the integer label 'lab' into the 'ipos' position
c  of the working precision word 'word'.
c
      implicit integer(a-z)
*@ifdef int64
*       integer word
*       integer iword(2)
*       integer m32
*      parameter(m32=2**32-1)
*c
*@ifdef T3E64 cray
*       iword(1)=shiftr(word,32)
*       iword(2)=and(word,m32)
*       iword(ipos)=and(lab,m32)
*       word= or(shiftl(iword(1),32),iword(2))
*@else
*       iword(1)=ishft(word,-32)
*       iword(2)=iand(word,m32)
*       iword(ipos)=iand(lab,m32)
*       word= ior(ishft(iword(1),32),iword(2))
*@endif
*@else
       integer word(2)
       integer lab,ipos
       word(ipos)=lab
*@endif
      return
      end

ctm low32 stems from ciden9.f

      integer function low32( wkword )
c
c     returns the "lowest" 32 bits of the working precision word wkword
c
*@ifdef int64
*@ifdef T3E64 cray
*      integer wkword
*      low32 = and( wkword, mask(96) )
*@else
*       integer wkword,m32
*       parameter (m32=2**32-1)
*       low32=iand(wkword,m32)
*@endif
*@else
      integer*4 wkword(2)
      low32=wkword(2)
*@endif
      return
      end

c changed 13.1.98
      subroutine packlp( word, code, head, bver, ybt, ykt)
ctm
ctm old version
ctm   subroutine packlp( word, code, head, xbt, ybt, ykt)
ctm
c
c  28-jun-90 wlkxmx->xbarmx, wlkmx->nwlkmx changes. -rls
c  09-apr-90 rwlphd->head change completed. -rls
c  17-apr-89 ft packing check moved here from wloop. -rls
c  modification of packing in  2-12-87 (p. szalay)
c                       (also in ci program!!!)
c
c  note:  if bits become scarce, then xbt may be eliminated in the
c         ft since it depends on head which is now explicitly
c         specified in each ft entry.
c         of course, the programs which read the ft will require
c         modification consistent with this change also. -rls
c
      implicit integer(a-z)
c
      integer   nrowmx,        xbarmx,         nwlkmx
      parameter(nrowmx=1023    , nwlkmx=2**20-1)

c      code  head  bver  ykt  ybt
c        4    10    2    24   24             64bit
c        4    10    2   16/ 8  24            2*32bit

      integer code, head, ybt, ykt
      integer bver

c     # for enforcement of ft packing limits. note that xbarmx
c     # may be 2 to 4 times smaller than nwlkmx since xbt is only
c     # addressed above the first level and only for a single chaining
c     # scheme at any one time.  Specifically xbt=xbar(head,ichain)
c     # whereas nwlkmx = sum( xbar(1:4,1:3) ). -rls
c
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
*@ifdef int64
*      integer word,word1
*       integer oder,links
*@ifdef cray T3E64
*       oder(i,j)=or(i,j)
*       links(i,j) = shiftl(i,j)
*@else
*       oder(i,j) = ior(i,j)
*       links(i,j) = ishft(i,j)
*@endif
*@else
      real*8 word, word1
c     # for 32-bit integer machines...
      integer iword(2)
      equivalence (word1,iword(1))
*@endif

c
c
c     # check items to be packed...
*@ifdef debug
*      if(head.gt.nrowmx)then
*         call bummer('packlp: head=',head,faterr)
*      elseif(ybt.gt.nwlkmx)then
*         call bummer('packlp: ybt=',ybt,faterr)
*      elseif(ykt.gt.nwlkmx)then
*         call bummer('packlp: ykt=',ykt,faterr)
*      endif
*@endif
c
*@ifdef int64
*      word1=oder(oder(oder(oder(links(ybt,24),ykt),links(bver,48)),
*     + links(head,50)),links(code,60))
*@else
*@if defined  (oldoldsun) || defined ( fps)
*      iword(1)=or(or(or(lshift(code,28),lshift(head,18)),
*     + lshift(bver,16)),rshift(ybt, 8))
*      iword(2)=or(lshift(ybt,24),ykt)
*@else
c  32-bit integer machines with f90 bit operators...
      iword(1)=ior(ior(ior(ishft(code,28),ishft(head,18)),
     + ishft(bver,16)),ishft(ybt,-8))
      iword(2)=ior(ishft(ybt,24),ykt)
*@endif
*@endif
c
      word = word1
c
      return
      end


      subroutine packstep( nintorb, istepv, output,stepno )
c
c ==============================================================
c compress each step vector istepv into nintorb*2 bits
c starting from (stepno-1)*nintorb*2+1 ... stepno*nintorb*2
c
c  input: nintorb   number of internal orbitals
c         istepv    step vector
c         stepno    number of the valid z walk
c  output: output   packed step vector
c
c          1 ... nintorb (ivalz=1)
c          1 ... nintorb (ivalz=2) ....
c
c tm october, 1st 1999
c===============================================================
      implicit logical(a-z)
       integer stepno, nintorb,istepv(nintorb),output(*)

      integer firstbit,firstwordm1,ii,intorb,ibit
      integer wordlen
*@ifdef int64
*        parameter(wordlen=64)
*@else
      parameter(wordlen=32)
*@endif


      firstbit=(stepno-1)*nintorb*2
      firstwordm1= (firstbit)/wordlen

      ibit=firstbit-firstwordm1*wordlen

      intorb=1
      ii=firstwordm1
      ii=ii+1

c  output must have been initialized with zeros!
100   output(ii) = ior(output(ii),ishft(istepv(intorb),ibit))
      ibit=ibit+2
      intorb=intorb+1
      if (intorb.gt.nintorb) goto 999
      if (ibit.gt.(wordlen-2)) then
         ibit=0
         ii=ii+1
      endif
      goto 100
999   return
      end
c
ctm This file contains a collection of various strongly machine
ctm dependent packing routines - on a long-term basis they
ctm have to generalized and to be included into colib3.f
ctm
c deck zdummy
*@ifdef flag .and. (.not. flag)
*      subroutine zdummy
*c
*c     # dummy routine to keep this file from being empty.
*c
*c     # bummer error types.
*      integer   wrnerr,  nfterr,  faterr
*      parameter(wrnerr=0,nfterr=1,faterr=2)
*c
*c     # this routine should never be called.
*c
*      call bummer( 'zdummy:', 0, faterr )
*c
*      return
*      end
*@endif
c deck rl
      integer function rl(i)
c
c     this function maps real word lengths on integer word lengths
c
c**********************************************************************
c  the use of this function is dangerous, due to machine-dependent data
c  alignment constraints, and should be eliminated in the calling
c  programs. -rls
c**********************************************************************
c
      integer i
c
c      integer first
c      save    first
c      data first / 1 /
c
c      if ( first .eq. 1 ) then
cc        # print an annoying message the first time through. -rls
c         call bummer(' rl(): obsolete function called.'
c     &    //' the calling program should be updated', 0, 0)
c         first = 0
c      endif
c
*@ifdef int64
*      rl= i
*@else
c     # assuming i*4 declarations
      rl=2*i
*@endif
      return
      end

      subroutine unplp(word, code, rwlphd, wtlphd, mlp2, mlp1, xbar)
c
c     modification of packing   2-12-87 (p. szalay)
c                       (also in ftape program!!!)
c
      implicit logical(a-z)
      integer nrowmx
      parameter (nrowmx=1023)
c
      integer code, rwlphd, wtlphd, mlp2, mlp1
      integer bver, xbar(nrowmx,3)
c
      integer ybta, ybtb,mask10,mask24
      real*8 word1
      integer mask4, mask8, mask14, mask7,mask16, mask18, maskex
      parameter (mask4=2**4-1, mask8=2**8-1,mask10=2**10-1)
      parameter (mask14=2**14-1, mask16=2**16-1, mask7=2**7-1)
      parameter (mask24=2**24-1, maskex=2**24-2**8 )
*@ifdef int64
*       integer word
*@ifdef cray or T3E64
*c mask(i) returns  for  0<=i<=64  i left adjusted bits of 1
*c                  for  65<=i<=128 128-i righ adjusted bits of 1
*      integer   and, shiftr
*      intrinsic and, shiftr
*      code   = and(shiftr(word,60),mask(124))
*      rwlphd = and(shiftr(word,50),mask(118))
*      bver   = and(shiftr(word,48),mask(126))
*      wtlphd = xbar(rwlphd,bver)
*      mlp2   = and(shiftr(word,24),mask(104))
*      mlp1   = and(word,mask(104))
*@else
*       integer m4,m8,m16
*       parameter (m4=2**4-1, m8=2**8-1, m16=2**16-1)
*      code   = iand(ishft (word,-60),mask4)
*      rwlphd = iand(ishft (word,-50),mask10 )
*      bver   = iand(ishft (word,-48),3)
*      wtlphd = xbar(rwlphd,bver)
*      mlp2   = iand(ishft (word,-24),mask24  )
*      mlp1   = iand(word,mask24   )
*@endif
*@elif defined  (oldoldsun) || defined ( fps)
*      integer   and, or, rshift, lshift
*      intrinsic and, or, rshift, lshift
*      real*8 word
*      integer iword(2)
*      equivalence (word1,iword(1))
*      word1 = word
*      code   = and(rshift(iword(1),28),mask4)
*      rwlphd = and(rshift(iword(1),18),mask10)
*      bver   = and(rshift(iword(1),16),3)
*c     wtlphd = and(rshift(iword(1), 4),mask16)
*      wtlphd = xbar(rwlphd,bver)
*      ybta   = and(lshift(iword(1), 8),maskex)
*      ybtb   = and(rshift(iword(2),24),mask8 )
*      mlp2   = or(ybta,ybtb)
*      mlp1   = and(iword(2),mask24)
*@else
c     # 32-bit machines with vax/ibm/mil-std-1753/f90 bit operators.
      integer   iand, ishft, ior
      intrinsic iand, ishft, ior
      real*8 word
      integer iword(2)
      equivalence (word1,iword(1))
      word1 = word
      code   = iand(ishft(iword(1),-28),mask4)
      rwlphd = iand(ishft(iword(1),-18),mask10)
      bver   = iand(ishft(iword(1),-16),3)
      wtlphd= xbar(rwlphd,bver)
      ybta   = iand(ishft(iword(1), 8),maskex)
      ybtb   = iand(ishft(iword(2),-24),mask8 )
      mlp2   = ior(ybta,ybtb)
      mlp1   = iand(iword(2),mask24)
*@endif
      return
      end
c
