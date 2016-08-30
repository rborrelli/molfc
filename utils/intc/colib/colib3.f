ccolib3.f
ccolib part=3 of 9.  general utility library routines.
cversion 5.0
c deck plblks
      subroutine plblks( title, z, nblk, nrow, labr, ifmt, nlist )
c
c  print a lower-triangular blocked matrix.
c
c  input:
c  title  = character title to be printed before each block.
c  z(*)   = blocked matrix to be printed.
c  nblk   = number of blocks in the matrix z(*).
c  nrow(*)= number of rows in each block.
c  labr   = character*4 row label.
c  ifmt   = format type.
c  nlist  = output unit number.
c
c  ron shepard 17-may-84.
c
      implicit logical(a-z)
      integer nblk, ifmt, nlist
      character*(*) title
      character*(*) labr
      integer nrow(nblk)
      real*8 z(*)
c
      integer i, nr, r0, zpt
c
      r0=0
      zpt=1
      do 100 i=1,nblk
         write(nlist,6010)title,i
         nr=nrow(i)
         if(nr.gt.0)call plblk(z(zpt),nr,r0,labr,ifmt,nlist)
         r0=r0+nr
         zpt=zpt+(nr*(nr+1))/2
100   continue
      return
6010  format(/10x,a,' block',i4)
      end
c deck plblkt
      subroutine plblkt( title, z, nr, labr, ifmt, nlist )
c
c  print a titled lower-triangular packed matrix.
c
c  ron shepard 17-may-84.
c
      implicit logical(a-z)
      integer nr, ifmt, nlist
      character*(*) title
      character*(*) labr
      real*8 z(*)
c
      write(nlist,6010)title
      call plblk(z,nr,0,labr,ifmt,nlist)
      return
6010  format(/10x,a)
      end
c deck plblk
      subroutine plblk( z, nr, r0, labr, ifmt, nlist )
c
c  print a lower-triangular packed matrix.
c  this version prints eight columns at a time with three formats.
c  parameter ncol and formats 10 and 1-3 should be modified to print
c  a different number of columns or to use different formats.
c
c  input:
c  z(*) = matrix to be printed.
c  nr   = row and column dimension.
c  r0   = row number offset.
c  labr = character row and column label.
c  ifmt   = format type (1:f, 2:e, 3:g).
c  nlist= output unit nubmer.
c
c  format statement assignment version 5-jun-87 (rls).
c  ron shepard 17-may-84.
c
      implicit logical(a-z)
      integer nr, ifmt, nlist
      character*(*) labr
      real*8 z(*)
c
      integer    ncol
      parameter (ncol=8)
10    format(/8x,8(6x,a4,i4,1x))
1     format(1x,a4,i4,8f15.8)
2     format(1x,a4,i4,1p,8e15.6)
3     format(1x,a4,i4,1p,8g15.6)
c
      real*8    zero
      parameter(zero=0d0)
c
      integer fmtz, i, ij0, ilab, j, j2, jlab1, jlab2, jlast, jstrt,
     & jt, r0
c
      if(ifmt.le.1)then
         assign 1 to fmtz
      elseif(ifmt.eq.2)then
         assign 2 to fmtz
      else
         assign 3 to fmtz
      endif
c
      jlast=0
      do 400 jstrt=1,nr,ncol
         jlast=min(nr,jlast+ncol)
c
         jlab1=jstrt+r0
         jlab2=jlast+r0
         write(nlist,10)(labr,j,j=jlab1,jlab2)
c
         ij0=(jstrt*(jstrt-1))/2
         do 300 i=jstrt,nr
            ilab=i+r0
            j2=min(i,jlast)
c
c  print the row if a nonzero element is found.
c
            do 100 j=jstrt,j2
               if(z(ij0+j).ne.zero)then
                  write(nlist,fmtz)labr,ilab,(z(ij0+jt),jt=jstrt,j2)
                  go to 101
               endif
100         continue
101         ij0=ij0+i
300      continue
400   continue
c
      return
      end
c deck prblks
      subroutine prblks(
     & title, z, nblk, nrow, ncol, labr, labc, ifmt, nlist )
c
c  print a rectangular-packed blocked matrix.
c
c  input:
c  title  = character title to be printed before each block.
c  z(*)   = blocked rectangular matrix to be printed.
c  nblk   = number of blocks in the matrix z(*).
c  nrow(*)= number of rows in each block.
c  ncol(*)= number of columns in each block.
c  labr   = character*4 row label.
c  labc   = character*4 column label.
c  ifmt   = format type.
c  nlist  = output unit number.
c
c  ron shepard 17-may-84.
c
      implicit logical(a-z)
      integer nblk, ifmt, nlist
      character*(*) title
      character*(*) labr,labc
      integer nrow(nblk),ncol(nblk)
      real*8 z(*)
c
      integer c0, i, nc, nr, nrnc, r0, zpt
c
      r0=0
      c0=0
      zpt=1
      do 100 i=1,nblk
         write(nlist,6010)title,i
         nr=nrow(i)
         nc=ncol(i)
         nrnc=nr*nc
         if(nrnc.gt.0)call prblk
     &    (z(zpt),nr,nr,nc,r0,c0,labr,labc,ifmt,nlist)
         r0=r0+nr
         c0=c0+nc
         zpt=zpt+nrnc
100   continue
      return
6010  format(/10x,a,' block',i4)
      end
c deck prblkt
      subroutine prblkt(
     & title, z, nrd, nr, nc, labr, labc, ifmt, nlist )
c
c  print a titled rectangular matrix.
c
c  ron shepard 17-may-84.
c
      implicit logical(a-z)
      integer nrd, nr, nc, ifmt, nlist
      character*(*) title
      character*(*) labr,labc
      real*8 z(*)
c
      write(nlist,6010)title
      call prblk(z,nrd,nr,nc,0,0,labr,labc,ifmt,nlist)
      return
6010  format(/10x,a)
      end
c deck prblk
      subroutine prblk(
     & z, nrd, nr, nc, r0, c0, labr, labc, ifmt, nlist )
c
c  print a sub-block of a rectangular matrix.
c  this version prints eight columns at a time with three formats.
c  parameter ncol and formats 10 and 1-3 should be modified to print
c  a different number of columns or to use different formats.
c
c  input:
c  z(*) = matrix to be printed.
c  nrd  = row dimension.
c  nr   = number of rows to print.
c  nc   = column dimension.
c  r0   = row number offset.
c  c0   = column number offset.
c  labr = character row label.
c  labc = character column label.
c  ifmt   = format type (1:f, 2:e, 3:g).
c  nlist= output unit number.
c
c  format statement assignment version 5-jun-87 (rls).
c  ron shepard 17-may-84.
c
      implicit logical(a-z)
      integer nrd, nr, nc, r0, c0, ifmt, nlist
      character*(*) labr,labc
      real*8 z(nrd,nc)
c
      integer    ncol
      parameter (ncol=8)
10    format(/8x,8(6x,a4,i4,1x))
1     format(1x,a4,i4,8f15.8)
2     format(1x,a4,i4,1p,8e15.6)
3     format(1x,a4,i4,1p,8g15.6)
c
      integer fmtz, i, ilab, j, jlab1, jlab2, jlast, jstrt, jt
c
      real*8     zero
      parameter (zero=0d0)
c
      if(ifmt.le.1)then
         assign 1 to fmtz
      elseif(ifmt.eq.2)then
         assign 2 to fmtz
      else
         assign 3 to fmtz
      endif
c
      jlast=0
      do 400 jstrt=1,nc,ncol
         jlast=min(nc,jlast+ncol)
c
         jlab1=jstrt+c0
         jlab2=jlast+c0
         write(nlist,10)(labc,j,j=jlab1,jlab2)
c
         do 300 i=1,nr
            ilab=i+r0
c
c  print the row if a nonzero element is found.
c
            do 100 j=jstrt,jlast
               if(z(i,j).ne.zero)then
                  write(nlist,fmtz)labr,ilab,(z(i,jt),jt=jstrt,jlast)
                  go to 300
               endif
100         continue
300      continue
400   continue
c
      return
      end
c deck prvblk
      subroutine prvblk(
     & z, v, nrd, nr, nc, r0, c0, labr, labc, labv, ifmt, nlist )
c
c  print a subblock of a rectangular matrix and a corresponding vector.
c  this version prints 8 columns at a time with one of three formats.
c  parameter ncol and the appropriate formats should be modified to
c  print a different number of columns or to use different formats.
c
c  input:
c  z(*,*)= matrix to be printed.
c  v(*)  = vector to be printed.
c  nrd   = row dimension.
c  nr    = number of rows to print.
c  nc    = column dimension.
c  r0    = row number offset.
c  c0    = column number offset.
c  labr  = character row label.
c  labc  = character column label.
c  labv  = character vector label.
c  ifmt  = format type (1:f, 2:e, 3:g).
c  nlist = output unit number.
c
c  format statement assignment version 5-jun-87 (rls).
c  ron shepard 17-may-84.
c
      implicit logical(a-z)
      integer nrd, nr, nc, r0, c0, ifmt, nlist
      character*(*) labr,labc
      character*(*) labv
      real*8 z(nrd,nc),v(nc)
c
      integer    ncol
      parameter (ncol=8)
10    format(/8x,8(6x,a4,i4,1x))
1     format(1x,a4,i4,8f15.8)
2     format(1x,a4,i4,1p,8e15.6)
3     format(1x,a4,i4,1p,8g15.6)
11    format(/1x,a8,8f15.8)
12    format(/1x,a8,1p,8e15.6)
13    format(/1x,a8,1p,8g15.6)
c
      integer fmtv, fmtz, i, ilab, j, jlab1, jlab2, jlast, jstrt, jt
c
      real*8     zero
      parameter (zero=0d0)
c
      if(ifmt.le.1)then
         assign 1 to fmtz
         assign 11 to fmtv
      elseif(ifmt.eq.2)then
         assign 2 to fmtz
         assign 12 to fmtv
      else
         assign 3 to fmtz
         assign 13 to fmtv
      endif
c
      jlast=0
      do 400 jstrt=1,nc,ncol
         jlast=min(nc,jlast+ncol)
c
         jlab1=jstrt+c0
         jlab2=jlast+c0
         write(nlist,10)(labc,j,j=jlab1,jlab2)
         write(nlist,fmtv)labv,(v(j),j=jstrt,jlast)
         write(nlist,*)
c
         do 300 i=1,nr
            ilab=i+r0
c
c  print the row if a nonzero element is found.
c
            do 100 j=jstrt,jlast
               if(z(i,j).ne.zero)then
                  write(nlist,fmtz)labr,ilab,(z(i,jt),jt=jstrt,jlast)
                  go to 300
               endif
100         continue
300      continue
400   continue
c
      return
      end
c deck getima
c deck igetim
      subroutine getima( numva, timesj )
c
c  time/statistics interface routine.
c  return the current timing statistics for this job in array timesj(*).
c
c  usage:
c      real*8      times1(numva), times2(numva)
c      character*8 ctimes(numva)
c      call igetim( numva, ctimes, numret ) # initialization
c      call getima( numva, times1 )
c      ...code to be timed...
c      call getima( numva, times2 )
c      diff(i) = abs(times2(i)-times1(i))  # i = 1, numret
c
c  input:
c  numva = number of values allocated in the calling program. this is
c          the dimension of timesj(*) and ctimes(*) in the calling
c          program. no more than numva elements of timesj(*) are
c          referenced in this routine.
c          (numva .ge. 1) must be satisfied.
c
c  output:
c  timesj(*) = array for returned values.  the individual entries
c              may be either increasing or decreasing.  the first
c              entry must be the cpu time in seconds (see below).
c              other entries are machine dependent and are defined
c              by entry igetim().  in general, these should be
c              ordered by decreasing importance to account for cases
c              in which the calling program cannot store all of the
c              computed values.
c  ctimes(*) = character array defining the timesj(*) entries.  these
c              array entries are no more than 8 characters in length.
c
c  19-oct-01 f90 code added. -rls
c  24-apr-92 ibm rs6000 code added, code merged from getime(). -rls
c  09-apr-92 fujitsu_vp code added. (Ross Nobes, Roger Edberg) -rls
c  13-mar-91 posix version. -rls
c  16-mar-90 written by ron shepard with suggestions by don comeau,
c            eric stahlberg, and robert harrison.
c
c======================================================================
c  the design of this routine is based on the observation that it is
c  impossible to determine ap priori a set of timing statistics that
c  are both important and available on all machines.  this routine
c  allows these timing values to be defined on a machine-by-machine
c  basis.  these values, along with the definitions encoded in
c  character strings, are then returned to the calling program, which
c  may be written in a machine-independent manner, leaving all of the
c  machine dependence localized within this routine.
c
c  with the proper choice of defining character strings, this also
c  removes any ambiguity associated with an imprecise definition of the
c  returned values.  for example, is the returned cpu time on a
c  parallel machine the total over all processors, the average over
c  all active processors, or the time for the slowest thread on any
c  single processor?  different values may be appropriate on different
c  machines, and this routine allows any, or all, of these values to
c  be returned as necessary.
c
c  the higher-level routine timer() may be examined for examples of
c  machine-independent code that correctly handles an arbitrary number
c  of return values.
c
c  for simplicity, the first return value is defined to be the cpu
c  time in seconds.  this value should be suitable for the computation
c  of the mflops rate, using whatever conventions are appropriate for
c  that machine.  this appears to be the only important value that is
c  common to all machines, and for this reason it is treated in this
c  special manner.  this convention allows this routine to be called
c  as
c
c             real*8 flops, timera(1), timerb(1)
c             ....
c             call getima(1, timera(1) )
c             call work
c             call getima(1, timerb(1) )
c             flops = nflop / abs( timerb(1) - timera(1) )
c
c  in order to compute mflops rates.  note that timera(*) and
c  timerb(*) should be declared as arrays in order to satisfy standard
c  fortran calling conventions even though only a single element is
c  referenced.  note also the use of the abs() function in the flops
c  computation; the returned values from this routine may be either
c  increasing or decreasing.
c======================================================================
c
      implicit logical(a-z)
      integer numva, numret
      real*8 timesj(*)
      character*(*) ctimes(*)
c
      integer i
c
      integer   min
      intrinsic min
c
c  in the following sections of machine-dependent code, the following
c  local variables must be defined and/or determined:
c      numvl = number of local timesl(*) elements computed.  this
c              should be an integer parameter.
c      timesl(1:numvl) = local running timing statistics values.  this
c                        should be declared as a real*8 array.
c      ctimel(1:numvl) = local character*8 definitions of the timesj(*)
c                        entries.  these should be defined in a data
c                        statement. ctimel(*) is referenced in igetim().
*@ifdef f95
c
c     # f95 version returns cpu time and wallclock time.
c
      integer, parameter :: numvl = 2
      real*8 timesl(numvl)
      character(len=8), parameter ::
     & ctimel(numvl) = (/ 'cpu_time', 'walltime' /)
c
      integer count, count_rate
      real*8 cpu
      intrinsic cpu_time, system_clock
c
      call cpu_time( cpu )
      timesl(1) = cpu
c
      call system_clock( count, count_rate )  !wrap around is ignored
      timesl(2) = (count)
      timesl(2) = timesl(2) / (count_rate)
*@elif defined  posix
*c
*c     # return the user_time, user+system_time, child_usertime,
*c     # child_usertime+chile_system_time, and wall_elapsed_time.
*c     # pxftimes() returns values in clock-ticks.
*c     # pxftime() is accurate to only 1. sec.
*c
*      integer jtms, ierr, itime
*      logical uninit
*      integer   numvl
*      parameter(numvl=5)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*
*      save uninit, jtms
*c
*      data ctimel /
*     & '    user', 'user+sys', 'walltime', ' ch_user', 'ch_(u+s)'/
*c
*      data uninit / .true. /
*c
*c     # make sure the jtms structure has been initialized.
*      if ( uninit ) then
*         call pxfstructcreate( 'tms', jtms, ierr )
*         if ( ierr .ne. 0 ) return
*         uninit = .false.
*      endif
*c
*c     # update the jtms structure values.
*      call pxftimes( jtms, itime, ierr )
*      if ( ierr .ne. 0 ) then
*         return
*      endif
*c
*c     # extract the components.
*      call pxfintget( jtms, 'tms_utime',  itime, ierr )
*      if ( ierr .ne. 0 ) return
*      timesl(1) = (itime)
*      call pxfintget( jtms, 'tms_stime',  itime, ierr )
*      if ( ierr .ne. 0 ) return
*      timesl(2) = (itime) + timesl(1)
*      call pxfintget( jtms, 'tms_cutime', itime, ierr )
*      if ( ierr .ne. 0 ) return
*      timesl(4) = (itime)
*      call pxfintget( jtms, 'tms_cstime', itime, ierr )
*      if ( ierr .ne. 0 ) return
*      timesl(5) = (itime) + timesl(4)
*c
*c     # use pxftime() for wall_time to avoid wraparound in pxftimes().
*      call pxftime( itime, ierr )
*      if ( ierr .ne. 0 ) return
*      timesl(3) = (itime)
*c
*@elif defined  rs6000
*c
*c     # return user time plus user and system times for subprocesses.
*c     # see also $COLUMBUS/special/unix/* for the library routine
*c     # fwtime().
*c
*      integer   numvl
*      parameter(numvl=2)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel / '    user', 'walltime' /
*c
*      real*8     hundth
*      parameter( hundth=1d-2 )
*c
*      real*8   fwtime
*      external fwtime
*c
*      integer  mclock
*      external mclock
*c
*      timesl(1) = ( mclock() ) * hundth
*      timesl(2) = fwtime()
*@elif defined  ibm
*c
*c     # version log:
*c     # 04-may-92 code moved from getime(). -rls
*c     # 28-mar-90 use standard subroutine calls (rmp)
*c     # 07-may-88 added vm/cms code (rmp)
*c     # 12-oct-87 added wall clock time (dcc)
*c
*      integer nret
*      real*8 cpu, elap
*c
*      real*8     micro
*      parameter( micro=1d-6 )
*c
*      integer   numvl
*      parameter(numvl=2)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel/'     cpu', 'walltime'/
*c
*      call cputime( cpu, nret )
*      if( nret .ne. 0 ) write (*,*) 'cpu clock reset'
*      call clockx( elap )
*c
*      timesl(1) = cpu  * micro
*      timesl(2) = elap * micro
*@elif defined  hp
*c
*c     # return user time plus user and system times for subprocesses.
*c     # see also $COLUMBUS/special/unix/* for the library routine
*c     # fwtime().
*c
*      integer   numvl
*      parameter(numvl=2)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel / '    user', 'walltime' /
*c
*c
*      real*8   fwtime
*      external fwtime
*c
*      real*8   runsec
*      external runsec
*c
*      timesl(1) = runsec()
*      timesl(2) = fwtime()
*@elif defined  delta
*      integer   numvl
*      parameter(numvl=1)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel / 'walltime' /
*c
*      real*8   dclock
*      external dclock
*c
*      timesl(1) = dclock()
*@elif defined  fujitsu
*c
*c     # unix version specific for fujitsu vp.
*c     # return the cpu_time and vu_time.
*c
*      integer   numvl
*      parameter(numvl=2)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel / '     cpu', '      vu' /
*c
*      call clockv( timesl(2), timesl(1), 0, 2 )
*@elif defined  vax
*c
*c     # vax returns cpu_time, process_elapsed_time, nio, nfault.
*c     # code taken from getime() written by ray bair.
*c     # getime() history:
*c     #     vax timing routine - written by ray bair 6/11/79
*c     #     page fault statistics added by ray bair -- 5/16/82
*c
*      integer*4 ihour, imin, isec, ihun, itim
*      integer*4 jpi$_dirio, jpi$_cputim, jpi$_pageflts
*      integer*4 long(10)
*      integer*2 list(20)
*      equivalence(list(1),long(1))
*      integer*4 dirio, cputim, faults
*      integer*2 time(7), hour, minute, sec, hun, last, day
*      save last
*      equivalence (hour,time(4)),(minute,time(5)),(sec,time(6)),
*     * (hun,time(7)),(day,time(3))
*c
*      real*8     hundth
*      parameter( hundth=0.01d0 )
*c
*      integer   numvl
*      parameter(numvl=4)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel / '     cpu', '    elap', '     nio', '  nfault' /
*c
*      data last/-1/
*      data jpi$_dirio, jpi$_cputim, jpi$_pageflts
*     & /   z0000040b,  z00000407,   z0000040a      /
*      data long/10*0/
*c
*      list(1)  = 4
*      list(2)  = jpi$_dirio
*      long(2)  = %loc(dirio)
*      list(7)  = 4
*      list(8)  = jpi$_cputim
*      long(5)  = %loc(cputim)
*      list(13) = 4
*      list(14) = jpi$_pageflts
*      long(8)  = %loc(faults)
*c
*      call sys$getjpi(,,,list,,,)
*      call sys$numtim(time,)
*c
*      if ( last .lt. 0 )   last = day
*      if ( last .ne. day ) hour = hour + 24
*      ihour  = hour
*      imin   = minute
*      isec   = sec
*      ihun   = hun
*      itim   = ihun + 100 * (isec + 60 * (imin + 60 * ihour) )
*c
*c     # convert times to standard units (seconds).
*c
*      timesl(1) = hundth * (cputim)
*      timesl(2) = hundth * (itim)
*      timesl(3) = (dirio)
*      timesl(4) = (faults)
*c
*@elif defined  (T3E64) || defined ( T3D )
*c   Since T3D runs essentially in dedicated mode no cpu clock timer
*c   is available; instead wall clock time is used
*C    SECONDR returns walltime
*C    timef available on T3D ???
*      real*8   timef
*      external timef
*      integer   numvl
*      parameter(numvl=2)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel / 'cpu_time', 'walltime' /
*c
*      call secondr(timesl(1))
*      timesl(2) = timef() / 1000
*@elif defined  unicos
*c
*c     # cray version returns cpu_time and wall_time.
*c     # timef() value is returned in milliseconds.
*c
*      real*8   timef
*      external timef
*      integer   numvl
*      parameter(numvl=2)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel / 'cpu_time', 'walltime' /
*c
*      call second( timesl(1) )
*      timesl(2) = timef() / 1000
*@elif defined (  decalpha) || defined ( unix) || defined ( linux) || defined ( convex)
*c
*c     # this is a "generic" bsd  version:
*c     # return the user_time, user+system_time, and wall_elapsed_time.
*c     # on most machines etime() is accurate to only 0.01 sec.
*c     # time() is accurate to only 1. sec.
*c
*      real*4 tarray(2), etime
*      real*8 fwtime
*      external etime, fwtime
*c
*      integer   numvl
*c     parameter(numvl=3)
*      parameter(numvl=2)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*c     data ctimel / '    user', 'user+sys', 'walltime' /
*      data ctimel / 'user+sys', 'walltime' /
*c
*c     timesl(2) = ( etime(tarray) )
*c     timesl(2) = tarray(2)
*c     timesl(1) = tarray(1)
*      timesl(1) = ( etime(tarray))
*      timesl(2) = fwtime()
*@else
*c     # can't do anything else, so just return a counter.
*      integer ncount
*      save    ncount
*c
*      integer   numvl
*      parameter(numvl=1)
*      real*8 timesl(numvl)
*      character*8 ctimel(numvl)
*      data ctimel/'  ncount'/
*c
*      data ncount / 0 /
*c
*      ncount = ncount + 1
*      timesl(1) = (ncount)
*@endif
c
c     # copy the appropriate values to timesj(*)
c     # note that the number of returned values may be less
c     # than the number of computed local values, numvl.
c
      do 100 i = 1, min(numva, numvl)
         timesj(i) = timesl(i)
100   continue
c
      return
c
      entry igetim( numva, ctimes, numret )
c
c  this entry point is to determine numret and to return the appropriate
c  ctimes(*) definitions to the calling program.
c  the use of this entry point is optional in the sense that values
c  returned by getima() do not depend on this initialization.
c  however, the meaning of these values, which are defined by the
c  ctimes(*) array, is machine dependent.  therefore, it is recommended
c  that the calling program print these character strings along with
c  results determined from the returned values for definitness.
c
      numret = min( numva, numvl )
      do 200 i = 1, numret
         ctimes(i) = ctimel(i)
200   continue
c
      return
      end
c deck timer
      subroutine timer( title, oper, event, nlist )
c
c  event timer.  this timer keeps separate timing statistics for
c  up to nevtmx events, identified by their event number.
c
c  input:
c  title  = character string (.le.30 characters).
c  oper   = timing operation to be performed.
c         = tmin0 =0 initialize everything and activate the first event.
c         = tminit=1 initialize and activate a new event.
c         = tmprt =2 print the running timing statistics of the event.
c         = tmrein=3 print the timing and reinitialize the event.
c         = tmclr =4 print the timing and clear the event.
c         = tmsusp=5 suspend event timing without printing.
c         = tmresm=6 resume event timing without printing.
c  event  = event number (1 to nevtmx) for 2<=oper<=6.
c  nlist  = listing file unit number.
c
c  output:
c  event  = new event number (1 to nevtmx) for oper=0,1.
c
c  version log:
c  16-mar-90  getima() version written. -rls
c  15-mar-90  /timex0/ reappeard and was removed again.
c             some redundant and inconsistent oper values were
c             removed.  suspend and resume were added. -rls
c  17-aug-88  options by p.g. szalay merged (eas)
c  09-may-88  bummer() call added (rls).
c  07-may-88  harris and ibm vm/cms versions (russ pitzer).
c  12-oct-87  /timex0/ removed, no ibm call, error checking,
c             flexible output (don comeau).
c  15-jun-87  written by ron shepard
c
      implicit logical(a-z)
      character*(*) title
      integer oper, event, nlist
c
c  local...
c     # nevtmx = the maximum number of events that can be simultaneously
c                timed.
c     # numvmx = the maximum number of times(*) values that can be
c                accomodated.  the number that are actually printed
c                is machine dependent and is determined by getima().
      integer   nevtmx,    numvmx
      parameter(nevtmx=32, numvmx=10)
c
      integer ievent, j, numval
c
      integer status(nevtmx)
c
c     # times(*)  = local copy of job times.
c     # timest(*) = local event totals.
c     # times0(*) = local event offsets.
c     # ctimes(*) = character array defining the times(*) values.
c
      real*8 times(numvmx), timest(numvmx,nevtmx), times0(numvmx,nevtmx)
      character*8 ctimes(numvmx)
c
      save numval, status, timest, times0, ctimes
c
      real*8    zero
      parameter(zero=0d0)
c
c     # note: it may be useful to copy the following parameter
c     #       definitions to the calling program, and use them
c     #       instead of the equivalent integer constants when
c     #       making calls to this routine. -rls
c
c     # timer operation parameters...
      integer   tmin0
      parameter(tmin0=0)
      integer   tminit,  tmprt,  tmrein,  tmclr,  tmsusp,  tmresm
      parameter(tminit=1,tmprt=2,tmrein=3,tmclr=4,tmsusp=5,tmresm=6)
c
c     # internal event status values
      integer   uninit,    ready,   active,   suspnd
      parameter(uninit=-1, ready=0, active=1, suspnd=2)
c
c     # ierror(oper,status)=1 if the combination is not allowed.
c     # note: in this version it is illegal to suspend a suspended event
c     #       or to resume an active event.
      integer ierror(tminit:tmresm, uninit:suspnd)
c
c     # iaccum(oper)=1 if accumulation is to be done.
      integer iaccum(tminit:tmresm)
c
c     # iprint(oper)=1 if cumulative timing is to be printed.
      integer iprint(tminit:tmresm)
c
c     # ireint(oper)=1 if cumulative values are to be reinitialized.
      integer ireint(tminit:tmresm)
c
c     # outst(oper,status) = output status after oper is performed.
      integer outst(tminit:tmresm, active:suspnd)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      data ierror/
     & 1,1,1,1,1,1,
     & 0,1,1,1,1,1,
     & 0,0,0,0,0,1,
     & 0,0,0,0,1,0/
c
      data iaccum/0,1,1,1,1,0/
c
      data iprint/0,1,1,1,0,0/
c
      data ireint/1,0,1,0,0,0/
c
      data outst/
     & active, active, active, ready, suspnd, active,
     & active, suspnd, active, ready, suspnd, active/
c     # set the status of all events to an uninitialized state.
c     # the first call to timer() should be with oper=tmin0.
      data status/nevtmx*uninit/
      data times/numvmx*zero/
c
c     # check the validity of oper.
      if ((oper.lt.tmin0).or.(oper.gt.tmresm)) then
         call bummer('timer: oper out of range, oper=', oper, wrnerr)
         return
      endif
c
c     # initialize the ctimes(*) array and the numval variable.
      if ( oper .eq. tmin0 ) then
         call igetim(numvmx, ctimes, numval)
      endif
c
c  the following call returns the cumulative times(*) values for this
c  job.  individual entries may be either increasing or decreasing.
c  e.g. the cpu time may be the time left for the job.  this routine
c  computes the absolute values of the incremental quantities and
c  accumulates the appropriate totals. the first entry must be the cpu
c  time in seconds suitable for computing mflops rates, but the other
c  entries are machine-dependent and are defined by the ctimes(*) array.
c
      call getima( numval, times )
c
      if ( oper .eq. tmin0 ) then
c        # initialize everything and assign the first event.
         do 140 ievent = 1, nevtmx
            status(ievent) = ready
c
c           # times0(*) values are for determining increments.
c           # timest(*) values are for event totals.
            do 120 j = 1, numvmx
               times0(j,ievent) = times(j)
               timest(j,ievent) = zero
120         continue
140      continue
         event         = 1
         status(event) = active
         return
      endif
c
      if ( oper .eq. tminit ) then
c        # determine a new event number.
         do 200 ievent = 1, nevtmx
            event = ievent
            if ( status(event) .eq. ready ) goto 220
200      continue
c        # loop exit means no event numbers are ready for use.  this
c        # means either that initialization is required or that all
c        # events have been used.
         if ( status(event) .eq. uninit ) then
            call bummer('timer: timer() has not been initialized',
     &       status(event),wrnerr)
            return
         else
c           # all events have been used, print a warning and reuse the
c           # last event number.
            if (nlist.gt.0) write(nlist,6020)event
            status(event) = ready
         endif
220      continue
      endif
c
c     # check for input consistency errors.
      if ( event .le. 0 .or. event .gt. nevtmx ) then
         call bummer('timer: illegal event, event=',event,wrnerr)
         return
      elseif ( ierror( oper, status(event) ) .eq. 1 ) then
         if (nlist.gt.0) write(nlist,6030) oper, status(event), event
         call bummer('timer: inconsistent operation/status combination'
     &    //' requested, oper=',oper,wrnerr)
         return
      endif
c
      if ( iaccum(oper) .eq. 1 .and. status(event) .eq. active ) then
c        # accumulate the latest increments.
         do 320 j = 1, numval
            timest(j,event) = timest(j,event)
     &       +                abs( times(j) - times0(j,event) )
320      continue
      endif
c
      if ( iprint(oper) .eq. 1 ) then
c        # print out the cumulative timings.
c        # either active or suspended events may be printed.
         if (nlist.gt.0)
     .     write(nlist,6010)title,(ctimes(j),timest(j,event),j=1,numval)
      endif
c
      if ( ireint(oper) .eq. 1 ) then
c        # activate the event and reinitialize the cumulative stats.
c        # ready, active, or suspended events may be reinitialized.
         status(event) = active
         do 420 j = 1, numval
            timest(j,event) = zero
420      continue
      endif
c
c     # reset the output status of the current event.
c     # at this time, the event is either active or suspended.
      status(event) = outst( oper, status(event) )
c
      if ( status(event) .eq. active ) then
c        # before exit, the incremental offsets should be
c        # updated for active events.
         do 520 j = 1, numval
            times0(j,event) = times(j)
520      continue
      endif
c
      return
c
6010  format(' !timer: ',a,t41,4(1x,a,'=',f10.3):/(6(1x,a,'=',f10.3)))
6020  format(' timer: all events are active.',
     & ' reusing the last event number.',i4)
6030  format(' timer: illegal or inconsistent timer operation. oper=',
     & i3,' status=',i3,' event=',i3)
      end

c deck forbyt
      integer function forbyt( i )
c
c  compute the number of working precision words required to hold
c  an (unpacked) integer array.
c
      integer i
c
*@ifdef int64 cray
*      forbyt = i
*@else
      forbyt = (i+1)/2
*@endif
      return
      end
c deck atebyt
      integer function atebyt( i )
c
c  compute the number of working precision words required to hold
c  a real (working precision) array.
c
      integer i
c
      atebyt = i
      return
      end
c deck unique
      subroutine unique( nline, line )
c
c  eliminate redundant lines.
c  input:
c  nline = number of initial lines.
c  line(*) = initial lines.
c
c  output:
c  nline = the number of unique lines.
c  line(1:nline) = unique lines.
c
      integer nline
      character*80 line(*)
c
      integer nmax, i, j
      external streq
      logical  streq
c
      nmax  = nline
      nline = 1
      do 10 i = 2, nmax
         do 20 j = 1, nline
            if ( streq( line(j), line(i) ) ) goto 10
20       continue
c        ...loop exit means a new unique line.
         nline = nline + 1
         if ( nline .ne. i ) line(nline) = line(i)
10    continue
      return
      end


      subroutine plab8( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+7)/8)*8 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
       implicit logical(a-z)
       integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw8
*c
*      nuw8=((nuw+7)/8)*8
*      if ( nuw8 .ne. 0 ) call pack(p,8,u,nuw8)
*@elif defined  int64
*      integer p(*),m8,m32
*      integer u(8,*)
*      integer nuw8
*c
*      parameter( m32=2**32-1)
*      parameter( m8=2**8-1)
*c
*C
*      integer jscr8,jscr16,jscr16b
*      integer moder,oder,links,und
*        integer i,j
*@if defined  (T3E64) || defined ( T3D )
*      moder(i,j)=or(i,and(j,m8))
*      oder (i,j)=or(i,j)
*      und  (i,j)=and(i,j)
*      links (i,j)= shiftl(i,j)
*@else 
*ctm: this is compiler specific!
*ctm   moder(i,j)=ior(i,kiand(j,m8))
*      moder(i,j)=ior(i,iand(j,m8))
*      oder (i,j)=ior(i,j)
*      und  (i,j)=iand(i,j)
*      links (i,j)=ishft(i,j)
*@endif
*c
*      nuw8=(nuw+7)/8
*       do  i=9-nuw8*8+nuw, 8
*       u(i,nuw8)=0
*       enddo
*      do 10 i=1,nuw8
*       jscr8=moder( links( u(1,i),8),u(2,i))
*       jscr16=oder (links(jscr8,16), oder(links(u(3,i),8),u(4,i) ))
*      jscr8= moder(links( u(5,i),8),u(6,i))
*       jscr16b= oder(links(jscr8,16), moder(links(u(7,i),8),u(8,i)) )
*       p(i) = oder(links(jscr16,32), und(jscr16b,m32) )
*
*10    continue
*@elif defined  oldoldsun
*c  17-apr-91 mask added to fix garbage overwrite bug. -rls/tom kovar
*      integer p(2,*)
*      integer u(8,*)
*c
*      integer nuw8
*c
*      integer   or, lshift, and
*      intrinsic or, lshift, and
*c
*      integer    m8
*      parameter( m8=2**8-1 )
*c
*      integer i, j
*      integer mor
*      mor(i,j) = or( i, and(j,m8) )
*c
*      nuw8=(nuw+7)/8
*      do 10 i=1,nuw8
*         p(1,i)=mor(lshift(mor(lshift(mor(lshift(
*     &    u(1,i),8),u(2,i)),8),u(3,i)),8),u(4,i))
*         p(2,i)=mor(lshift(mor(lshift(mor(lshift(
*     &    u(5,i),8),u(6,i)),8),u(7,i)),8),u(8,i))
*10    continue
*@else
c  general byte-addressable 32-bit integer machines.
      logical*1 p(8,*)
      logical*1 u(4,8,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+7)/8)
            p(8,i)=u(1,1,i)
            p(7,i)=u(1,2,i)
            p(6,i)=u(1,3,i)
            p(5,i)=u(1,4,i)
            p(4,i)=u(1,5,i)
            p(3,i)=u(1,6,i)
            p(2,i)=u(1,7,i)
            p(1,i)=u(1,8,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+7)/8)
            p(1,i)=u(4,1,i)
            p(2,i)=u(4,2,i)
            p(3,i)=u(4,3,i)
            p(4,i)=u(4,4,i)
            p(5,i)=u(4,5,i)
            p(6,i)=u(4,6,i)
            p(7,i)=u(4,7,i)
            p(8,i)=u(4,8,i)
20       continue
      endif
*@endif
      return
      end

      subroutine ulab8( p, u, nuw )
c
c  unpack integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+7)/8)*8 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  20-jul-91 cray nuw check added. -galen gawboy/rls
c  written by ron shepard.
c
       implicit logical(a-z)
       integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw8
*c
*      nuw8=((nuw+7)/8)*8
*      if ( nuw8 .ne. 0 ) call unpack(p,8,u,nuw8)
*@elif defined  int64
*      integer p(*)
*      integer u(8,*)
*c
*      integer nuw8
*c
*      integer   m8
*      parameter(m8=2**8-1)
*c
*       integer und,rechts
*        integer i,j
*@if defined  (T3E64) || defined ( T3D )
*       und(i,j) = and(i,j)
*       rechts(i,j) = shiftr(i,j)
*@else 
*       und(i,j) = iand(i,j)
*       rechts(i,j)=rshift(i,j)
*@endif
*
*c
*c     call izero(nuw,u,1)
*      nuw8=(nuw+7)/8
*      do 10 i=1,nuw8
*         u(1,i)=und(rechts(p(i),56),m8)
*         u(2,i)=und(rechts(p(i),48),m8)
*         u(3,i)=und(rechts(p(i),40),m8)
*         u(4,i)=und(rechts(p(i),32),m8)
*         u(5,i)=und(rechts(p(i),24),m8)
*         u(6,i)=und(rechts(p(i),16),m8)
*         u(7,i)=und(rechts(p(i),8),m8)
*         u(8,i)=und(       p(i),    m8)
*10    continue
*@elif defined  oldoldsun
*      integer p(2,*)
*      integer u(8,*)
*c
*      integer i, nuw8
*c
*      integer   m8
*      parameter(m8=2**8-1)
*c
*      integer   and, rshift
*      intrinsic and, rshift
*c
*      nuw8=(nuw+7)/8
*      do 10 i=1,nuw8
*         u(1,i)=and(rshift(p(1,i),24),m8)
*         u(2,i)=and(rshift(p(1,i),16),m8)
*         u(3,i)=and(rshift(p(1,i), 8),m8)
*         u(4,i)=and(       p(1,i),    m8)
*         u(5,i)=and(rshift(p(2,i),24),m8)
*         u(6,i)=and(rshift(p(2,i),16),m8)
*         u(7,i)=and(rshift(p(2,i), 8),m8)
*         u(8,i)=and(       p(2,i),    m8)
*10    continue
*@else
c  general byte-addressable 32-bit integer machines.
      logical*1 p(8,*)
      logical*1 u(4,8,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      call izero(nuw,u,1)
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+7)/8)
            u(1,1,i)=p(8,i)
            u(1,2,i)=p(7,i)
            u(1,3,i)=p(6,i)
            u(1,4,i)=p(5,i)
            u(1,5,i)=p(4,i)
            u(1,6,i)=p(3,i)
            u(1,7,i)=p(2,i)
            u(1,8,i)=p(1,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+7)/8)
            u(4,1,i)=p(1,i)
            u(4,2,i)=p(2,i)
            u(4,3,i)=p(3,i)
            u(4,4,i)=p(4,i)
            u(4,5,i)=p(5,i)
            u(4,6,i)=p(6,i)
            u(4,7,i)=p(7,i)
            u(4,8,i)=p(8,i)
20       continue
      endif
*@endif
      return
      end


      subroutine plab16( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
       implicit logical(a-z)
        integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw16
*c
*      nuw16=((nuw+3)/4)*4
*      if ( nuw16 .ne. 0 ) call pack(p,16,u,nuw16)
*@elif defined  int64
*       integer p(*), u(4,*)
*       integer m16,m32
*       parameter (m16=2**16-1, m32=2**32-1)
*       integer i,j,jscr32a, jscr32b,nuw16
*      integer moder,oder,links,und
*@if defined  (T3E64) || defined ( T3D )
*      moder(i,j)=or(i,and(j,m16))
*      oder (i,j)=or(i,j)
*      und  (i,j)=and(i,j)
*      links (i,j)= shiftl(i,j)
*@else 
*ctm      moder(i,j)=ior(i,kiand(j,m16))
*      moder(i,j)=ior(i,iand(j,m16))
*      oder (i,j)=ior(i,j)
*      und  (i,j)=iand(i,j)
*      links (i,j)=ishft(i,j)
*@endif
*
*       nuw16=(nuw+3)/4
*        do i=5 - nuw16*4 +nuw,4
*        u(i,nuw16)=0
*         enddo
*       do 10 i=1,nuw16
*       jscr32a= moder(links(u(1,i),16),u(2,i))
*       jscr32b= moder(links(u(3,i),16),u(4,i))
*       p(i)= oder(links(jscr32a,32),und(jscr32b,m32) )
*  10   continue
*@else
c     # general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer*2 u(2,4,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+3)/4)
            p(4,i)=u(1,1,i)
            p(3,i)=u(1,2,i)
            p(2,i)=u(1,3,i)
            p(1,i)=u(1,4,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+3)/4)
            p(1,i)=u(2,1,i)
            p(2,i)=u(2,2,i)
            p(3,i)=u(2,3,i)
            p(4,i)=u(2,4,i)
20       continue
      endif
*@endif
      return
      end



      subroutine ulab16( p, u, nuw )
c
c  unpack one-electron integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+3)/4)*4 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  20-jul-91 cray nuw check added. -galen gawboy/rls
c  written by ron shepard.
c
       implicit logical(a-z)
       integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw16
*c
*      nuw16=((nuw+3)/4)*4
*      if ( nuw16 .ne. 0 ) call unpack(p,16,u,nuw16)
*@elif defined  int64
*       integer p(*), u(4,*)
*       integer m16,nuw16
*       parameter (m16 = 2**16-1)
*       integer und,rechts
*        integer i,j
*@if defined  (T3E64) || defined ( T3D )
*       und(i,j) = and(i,j)
*       rechts(i,j) = shiftr(i,j)
*@else 
*       und(i,j) = iand(i,j)
*       rechts(i,j)=ishft(i,-j)
*@endif
*
*       nuw16=((nuw+3)/4)
*c      call izero (nuw , u, 1)
*       do 10 i=1, nuw16
*        u(1,i)=und(rechts(p(i),48),m16)
*        u(2,i)=und(rechts(p(i),32),m16)
*        u(3,i)=und(rechts(p(i),16),m16)
*        u(4,i)=und(p(i),m16)
* 10    continue
*@else
c     # general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer*2 u(2,4,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      call izero ( nuw, u, 1 )
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+3)/4)
            u(1,1,i)=p(4,i)
            u(1,2,i)=p(3,i)
            u(1,3,i)=p(2,i)
            u(1,4,i)=p(1,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+3)/4)
            u(2,1,i)=p(1,i)
            u(2,2,i)=p(2,i)
            u(2,3,i)=p(3,i)
            u(2,4,i)=p(4,i)
20       continue
      endif
*@endif
      return
      end


      subroutine plab32( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+1)/2)*2 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
       implicit logical(a-z)
       integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw32
*c
*      nuw32=((nuw+1)/2)*2
*      if ( nuw32 .ne. 0 ) call pack(p,32,u,nuw32)
*@elif defined  int64
*      integer p(*), u(2,*)
*      integer j,m32,i,nuw32
*      parameter(m32=2**32-1)
*      integer moder,oder,links,und
*@if defined  (T3E64) || defined ( T3D )
*      moder(i,j)=or(i,and(j,m32))
*      oder (i,j)=or(i,j)
*      und  (i,j)=and(i,j)
*      links (i,j)= shiftl(i,j)
*@else 
*c     moder(i,j)=ior(i,kiand(j,m32))
*      moder(i,j)=ior(i,iand(j,m32))
*      oder (i,j)=ior(i,j)
*      und  (i,j)=iand(i,j)
*      links (i,j)=ishft(i,j)
*@endif
*
*
*      nuw32=(nuw+1)/2
*      do 10 i=1,nuw32
*      p(i)=moder(links(u(1,i),32), und(u(2,i),m32))
* 10   continue
*@else
      integer p(2,*)
      integer u(2,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+1)/2)
            p(2,i)=u(1,i)
            p(1,i)=u(2,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+1)/2)
            p(1,i)=u(1,i)
            p(2,i)=u(2,i)
20       continue
      endif
*@endif
      return
      end


      subroutine ulab32( p, u, nuw )
c
c  unpack integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+1)/2)*2 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  20-jul-91 cray nuw check added. -galen gawboy/rls
c  written by ron shepard.
c
       implicit logical(a-z)
       integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw32
*c
*      nuw32=((nuw+1)/2)*2
*      if ( nuw32 .ne. 0 ) call unpack(p,32,u,nuw32)
*@elif defined  int64
*       integer  p(*), u(2,*)
*       integer m32,nuw32
*       parameter (m32=2**32-1)
*       integer und,rechts
*        integer i,j
*@if defined  (T3E64) || defined ( T3D )
*       und(i,j) = and(i,j)
*       rechts(i,j) = shiftr(i,j)
*@else 
*       und(i,j) = iand(i,j)
*       rechts(i,j)=ishft(i,-j)
*@endif
*       nuw32=((nuw+1)/2)
*       call izero(nuw,u,1)
*       do 10 i=1, nuw32
*         u(1,i) = und(rechts(p(i),32),m32)
*         u(2,i) = und(p(i),m32)
* 10   continue
*@else
      integer p(2,*)
      integer u(2,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+1)/2)
            u(2,i)=p(1,i)
            u(1,i)=p(2,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+1)/2)
            u(1,i)=p(1,i)
            u(2,i)=p(2,i)
20       continue
      endif
*@endif
      return
      end


      subroutine plab1( p, u, nuw )
c
c  pack a bit vector.
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+63)/64)*64 ) are referenced.
c  nuw  = number of unpacked array elements.
c
c  written by ron shepard.
c  version: 5-jul-89
c
       implicit logical(a-z)
       integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw64
*c
*      nuw64=((nuw+63)/64)*64
*      if ( nuw64 .ne. 0 ) call pack(p,1,u,nuw64)
*@elif defined  int64
*        integer p(*), u(64,*)
*        integer oder,local,i,j,mpack,nuw64
*@if defined  (T3E64) || defined ( T3D )
*        mpack(i,j)=shiftl(and(i,1),j)
*        oder(i,j) = or(i,j)
*@else 
*        mpack(i,j) =ishft(iand(i,1),j)
*        oder(i,j)  = ior(i,j)
*@endif
*        nuw64= (nuw+63)/64
*        do 10 i=1,nuw64
*         p(i)=0
*        do 11 j=63,0,-1
*         p(i)= or(p(i), mpack(u(64-j,i),j))
*   11   continue
*   10   continue
*@elif defined  vax
*c
*c  the following code may be used for 32-bit integer, little-endian
*c  machines that use vax-type bit operators.
*c  30-jul-90  written by ron shepard.
*c
*      integer p(2,*)
*      integer u(*)
*c
*      integer local1, local2, u0
*      integer sor, i, j
*      sor(i,j) = ior(ishft(i,1), iand(j,1) )
*c
*      do 20 i=1,((nuw+63)/64)
*c
*         u0 = (i-1)*64
*         local1 = 0
*         do 11 j = 1, 32
*            local1 = sor(local1,u(u0+j))
*11       continue
*         p(2,i) = local1
*c
*         u0 = (i-1)*64 + 32
*         local2 = 0
*         do 12 j = 1, 32
*            local2 = sor(local2,u(u0+j))
*12       continue
*         p(1,i) = local2
*c
*20    continue
*@elif defined   fujitsu_vp
*c
*c  the following code may be used for 32-bit integer, big-endian
*c  machines that use vax-type bit operators.
*c  30-nov-89  written by ron shepard.
*c
*c  the general code breaks the titan compiler. use this version
*c  instead. -rls
*      integer p(2,*)
*      integer u(*)
*c
*      integer local1, local2, u0
*      integer sor, i, j
*      sor(i,j) = ior(ishft(i,1), iand(j,1) )
*c
*      do 20 i=1,((nuw+63)/64)
*c
*         u0 = (i-1)*64
*         local1 = 0
*         do 11 j = 1, 32
*            local1 = sor(local1,u(u0+j))
*11       continue
*         p(1,i) = local1
*c
*         u0 = (i-1)*64 + 32
*         local2 = 0
*         do 12 j = 1, 32
*            local2 = sor(local2,u(u0+j))
*12       continue
*         p(2,i) = local2
*c
*20    continue
*@else
c  general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer u(*)
c
      integer i, u0
      integer local1, local2, local3, local4
      integer*2 llocl1(2), llocl2(2), llocl3(2), llocl4(2)
      equivalence (llocl1,local1),(llocl2,local2)
      equivalence (llocl3,local3),(llocl4,local4)
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
c  define bit-packing functions recursively for maximum pipelining...
      intrinsic mod
      integer lp2,lp4,lp8,lp16
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16
      lp2(i1,i2) = mod(i1,2) * 2 + mod(i2,2)
      lp4(i1,i2,i3,i4) = (lp2(i1,i2) * 4 + lp2(i3,i4))
      lp8(i1,i2,i3,i4,i5,i6,i7,i8)=
     & (lp4(i1,i2,i3,i4) * 16 + lp4(i5,i6,i7,i8))
      lp16(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16)=
     & ( lp8(i1,i2,i3,i4,i5,i6,i7,i8) * 256 +
     & lp8(i9,i10,i11,i12,i13,i14,i15,i16) )
c
      longw = 1
      if(shortw(1).eq.1)then
c        # ...little-endian.
         do 10 i=1,((nuw+63)/64)
c
            u0=(i-1)*64
            local1=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(4,i)=llocl1(1)
c
            u0=(i-1)*64+16
            local2=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(3,i)=llocl2(1)
c
            u0=(i-1)*64+32
            local3=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(2,i)=llocl3(1)
c
            u0=(i-1)*64+48
            local4=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(1,i)=llocl4(1)
c
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+63)/64)
c
            u0=(i-1)*64
            local1=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(1,i)=llocl1(2)
c
            u0=(i-1)*64+16
            local2=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(2,i)=llocl2(2)
c
            u0=(i-1)*64+32
            local3=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(3,i)=llocl3(2)
c
            u0=(i-1)*64+48
            local4=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(4,i)=llocl4(2)
c
20       continue
      endif
*@endif
      return
      end


      subroutine ulab1( p, u, nuw )
c
c  unpack a bit vector.
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+63)/64)*64 ) are referenced.
c  nuw  = number of unpacked array elements.
c
c  written by ron shepard.
c  version: 5-jul-89
c
       implicit logical(a-z)
       integer nuw
*@ifdef cray
*      real*8 p(*)
*      integer u(nuw)
*c
*      integer nuw64
*c
*      nuw64=((nuw+63)/64)*64
*      if ( nuw64 .ne. 0 ) call unpack(p,1,u,nuw64)
*@elif defined  int64
*        integer p(*), mpack,nuw64,u(64,*),j,l,i,k,plocal ,m2
*        parameter (m2=2)
*@if defined  (T3E64) || defined ( T3D )
*         mpack(l,k)= and(shiftr(l,k),1)
*@else 
*         mpack(l,k)= iand(ishft(l,-k),1)
*@endif
*        nuw64= ((nuw+63)/64)
*        call izero(nuw,u,1)
*        do 10 i=1, nuw64
*         do 11 j=0,63
*          u(j+1,i) =mpack(p(i),63-j)
*   11    continue
*   10    continue
*@elif defined  vax
*c
*c  the following code may be used for 32-bit integer, little-endian
*c  machines that use vax-type bit operators.
*c  30-jul-90  written by ron shepard.
*c
*      integer p(2,*)
*      integer u(*)
*c
*      integer local1, local2, u0
*      integer sand, i, j
*      sand(i,j)=iand(ishft(i,-j),1)
*c
*      do 20 i=1,((nuw+63)/64)
*c
*         u0 = (i-1)*64
*         local1=p(2,i)
*         do 11 j=31,0,-1
*            u0=u0+1
*            u(u0)=sand(local1,j)
*11       continue
*c
*         u0 = (i-1)*64 + 32
*         local2=p(1,i)
*         do 12 j=31,0,-1
*            u0=u0+1
*            u(u0)=sand(local2,j)
*12       continue
*c
*20    continue
*@elif defined   fujitsu_vp
*c
*c  the following code may be used for 32-bit integer, big-endian
*c   machines that use vax-type bit operators.
*c  11-jun-89  written by ron shepard.
*c
*c  the general code breaks the alliant compiler. use this instead. -rls
*      integer p(2,*)
*      integer u(*)
*c
*      integer local1, local2, u0
*      integer sand, i, j
*      sand(i,j)=iand(ishft(i,-j),1)
*c
*      do 20 i=1,((nuw+63)/64)
*c
*         u0 = (i-1)*64
*         local1=p(1,i)
*         do 11 j=31,0,-1
*            u0=u0+1
*            u(u0)=sand(local1,j)
*11       continue
*c
*         u0 = (i-1)*64 + 32
*         local2=p(2,i)
*         do 12 j=31,0,-1
*            u0=u0+1
*            u(u0)=sand(local2,j)
*12       continue
*c
*20    continue
*@else
c  general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer u(*)
c
      integer i, ipower, j, u0
      integer local1, local2, local3, local4
      integer*2 llocl1(2), llocl2(2), llocl3(2), llocl4(2)
      equivalence (llocl1,local1),(llocl2,local2)
      equivalence (llocl3,local3),(llocl4,local4)
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      intrinsic mod
c
      local1=0
      local2=0
      local3=0
      local4=0
c
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         u0=0
         do 10 i=1,((nuw+63)/64)
c
            llocl1(1)=p(4,i)
            ipower=2**15
            do 1 j=1,16
               u0=u0+1
               u(u0)=mod(local1/ipower,2)
               ipower=ipower/2
1           continue
c
            llocl2(1)=p(3,i)
            ipower=2**15
            do 2 j=1,16
               u0=u0+1
               u(u0)=mod(local2/ipower,2)
               ipower=ipower/2
2           continue
c
            llocl3(1)=p(2,i)
            ipower=2**15
            do 3 j=1,16
               u0=u0+1
               u(u0)=mod(local3/ipower,2)
               ipower=ipower/2
3           continue
c
            llocl4(1)=p(1,i)
            ipower=2**15
            do 4 j=1,16
               u0=u0+1
               u(u0)=mod(local4/ipower,2)
               ipower=ipower/2
4           continue
c
10       continue
      else
c        # ...big-endian.
         u0=0
         do 20 i=1,((nuw+63)/64)
c
            llocl1(2)=p(1,i)
            ipower=2**15
            do 11 j=1,16
               u0=u0+1
               u(u0)=mod(local1/ipower,2)
               ipower=ipower/2
11          continue
c
            llocl2(2)=p(2,i)
            ipower=2**15
            do 12 j=1,16
               u0=u0+1
               u(u0)=mod(local2/ipower,2)
               ipower=ipower/2
12           continue
c
            llocl3(2)=p(3,i)
            ipower=2**15
            do 13 j=1,16
               u0=u0+1
               u(u0)=mod(local3/ipower,2)
               ipower=ipower/2
13           continue
c
            llocl4(2)=p(4,i)
            ipower=2**15
            do 14 j=1,16
               u0=u0+1
               u(u0)=mod(local4/ipower,2)
               ipower=ipower/2
14           continue
c
20       continue
      endif
*@endif
      return
      end

c deck moerr
c *** this routine is incremental ***
      subroutine moerr( filerr, syserr, errtyp )
c
c subroutine moerr: evaluate errors produced by moread/mowrit
c    written by: eric stahlberg
c                mar 1990
c
c  input parameters:
c     filerr : mo error number
c     syserr : system error number (if filerr eq 1)
c     errtyp : type of error handling requested
c     (0 = wrnerr, 1 = nfaterr, 2 = faterr)
c
      integer filerr, syserr, errtyp
c
c local variables
c
      integer    lowerr,    higerr
      parameter (lowerr=-14,higerr=1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      if (filerr.lt.lowerr.or.filerr.gt.higerr) then
         call bummer('moerr: bad error number',filerr,faterr)
      endif
      if (errtyp.gt.faterr.or.errtyp.lt.wrnerr) then
         call bummer('moerr: bad error type',errtyp,faterr)
      endif
      if (filerr.eq.-1) then
         call bummer('moread: clen too small',filerr,errtyp)
      elseif (filerr.eq.-2) then
         call bummer('moread: ntitle > mxtitl',filerr,errtyp)
      elseif (filerr.eq.-3) then
         call bummer('moread: nsym > mxsym or nsym < 1',filerr,errtyp)
      elseif (filerr.eq.-4) then
         call bummer('moread: incompatible file version',filerr,errtyp)
      elseif (filerr.eq.-5) then
         call bummer('moread: invalid file opcode',filerr,errtyp)
      elseif (filerr.eq.-6) then
         call bummer('moread: eof before field match',filerr,errtyp)
      elseif (filerr.eq.-7) then
         call bummer('moread: symblk out of range',filerr,errtyp)
      elseif (filerr.eq.-8) then
         call bummer('moread: symblk > nsym',filerr,errtyp)
      elseif (filerr.eq.-9) then
         call bummer('mowrit: invalid file opcode',filerr,errtyp)
      elseif (filerr.eq.-10) then
         call bummer('mowrit: invalid write previous blk',filerr,errtyp)
      elseif (filerr.eq.-11) then
         call bummer('mowrit: blocks out of order',filerr,errtyp)
      elseif (filerr.eq.-12) then
         call bummer('mowrit: invalid symmetry dimension',filerr,errtyp)
      elseif (filerr.eq.-13) then
         call bummer('mowrit: invalid block specification',
     &    filerr,errtyp)
      elseif (filerr.eq.-14) then
         call bummer('mowrit: header is not first field',filerr,errtyp)
      elseif (filerr.eq.1) then
         call bummer('moread: system i/o error',syserr,errtyp)
      endif
      return
      end
c deck moinqu
c *** this routine is incremental ***
      integer function moinqu( unit, flcod )
c
c subroutine: moinqu
c purpose   : inquire existence of a mo coefficient file field
c             in the standard mocoeff file format
c author    : eric stahlberg, osu chemistry (mar 1990)
c
c version   : 2.1  (mar 1990)
c
c variable definitions
c
c *** in variables ***
c unit  : unit number of mo file
c flcod: file read operation code
c           10=read header only
c           20=read after header mocoef field
c           30=read after header energy field
c           40=read after header orbocc field
c           + block number to be read (0 means all blocks)
c         fopcod lt 10 implies all blocks and that fopcod
c
      integer unit, flcod
c
c     # local variables.
c
      integer iversn, filerr, syserr, ifld
      integer symblk, fopcod
      character*6 inname, field(4), fldfmt, flname, inline
c
c     # d2h and subgroups only.
c
      integer    maxsym
      parameter (maxsym=8)
c
c     # latest file version
c
      integer    lstver
      parameter (lstver=2)
c
c     # last valid fop code
c
      integer    lstfop
      parameter (lstfop=4)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # initialize some static constants.
c
      flname = 'mocoef'
      call allcap(flname)
c
      field(1) = 'header'
      field(2) = 'mocoef'
      field(3) = 'energy'
      field(4) = 'orbocc'
      do 10 ifld = 1, 4
         call allcap(field(ifld))
10    continue
      fldfmt  =' (a6) '
c
c     # parse out symmetry block and real fopcod from fopcod.
c
      if (flcod.ge.10) then
         symblk = flcod - (flcod/10)*10
         fopcod = (flcod-symblk) / 10
      else
         fopcod = flcod
         symblk = 0
      endif
c
c     # check fopcod for valid value.
c
      if ((fopcod.gt.lstfop).or.(fopcod.lt.1)) then
         filerr=-5
         goto 1000
      endif
c
c     # initialize state (assume success).
c
      filerr=0
      syserr=0
c
c     # rewind file.
c
      rewind (unit)
c
c     # read initial line for parsing.
c
      read (unit,5001,iostat=syserr) iversn, inname
      if ( syserr .ne. 0 ) goto 1100
5001  format(i2,a6)
      call allcap(inname)
c
c     # (optional checks on file version number go here)
c
      if ( (iversn .gt. lstver)
     & .or. (inname .ne. flname) ) then
         filerr=-4
         goto 1000
      endif
c
c     # loop until to get to next field.
c     # drop out of loop if field is reached or eof detected.
c
130   read (unit,fldfmt,end=137,iostat=syserr) inline
      if (syserr.ne.0) goto 1100
      call allcap(inline)
      if (inline.eq.field(fopcod)) goto 140
      goto 130
c
c     # end of file detected before field is matched.
c
137   continue
      filerr=-6
      goto 1000
c
c     # end of until loop
140   continue
c
      moinqu=0
      return
c
c     # error handling portion.
c
c     # system errors.
1100  continue
      call bummer('moread: system file error',syserr,wrnerr)
      moinqu=1
      return
c
c     # file errors.
1000  continue
      moinqu=filerr
      return
c
      end
c deck moread
c *** this routine is incremental ***
      subroutine moread(
     & unit,   flcod,  filerr, syserr,
     & mxtitl, ntitle, titles, afmt,
     & nsym,   nbfpsy, nmopsy, labels,
     & clen,   c  )
c
c subroutine: moread
c purpose   : universal interface read of mo coefficient file
c author    : eric stahlberg, osu chemistry (may 1989)
c
c version   : 2.1  (feb 1990)
c
c variable definitions
c
c *** in variables ***
c unit  : unit number of mo file
c flcod: file read operation code
c           10=read header only
c           20=read after header mocoef field
c           30=read after header energy field
c           40=read after header orbocc field
c           + block number to be read (0 means all blocks)
c         fopcod lt 10 implies all blocks and that fopcod
c mxtitl: maximum number of title cards which may be read
c clen  : available space in array c
c *** out variables ***
c filerr: mo file error number
c syserr: i/o system error - if filerr eq 1
c ntitle: number of title cards read
c titles: actual title cards
c afmt  : format with which values were read
c nsym  : number of irreps in point group
c nbfpsy: number of basis functions per mo of each irrep
c nmopsy: number of mos per irrep
c labels: symmetry labels for each irrep
c c     : coefficient array space
c
       implicit logical(a-z)
      integer unit, flcod, filerr, syserr, mxtitl, ntitle, nsym
      integer nbfpsy(*), nmopsy(*), clen
      real*8 c(*)
      character*80 titles(*)
      character*80 afmt
      character*4 labels(*)
c
c     # local variables.
c
      integer totspc, index, iversn, i, isym, imo, ifld
      integer symblk, fopcod
      logical ldio
      real*8 cdummy
      character*6 inname, field(4), fldfmt, flname, inline
c
c     # abelian point groups only.
c
      integer    maxsym
      parameter (maxsym=8)
c
c     # latest file version
c
      integer    lstver
      parameter (lstver=2)
c
c     # last valid fop code.
c
      integer    lstfop
      parameter (lstfop=4)
c
      real*8     zero
      parameter( zero=0d0 )
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # initialize some static constants.
c
      flname = 'mocoef'
      call allcap(flname)
c
      field(1) = 'header'
      field(2) = 'mocoef'
      field(3) = 'energy'
      field(4) = 'orbocc'
      do 10 ifld = 1, 4
         call allcap(field(ifld))
10    continue
      fldfmt  = ' (a6) '
c
c     # parse out symmetry block and real fopcod from fopcod.
c
      if (flcod.ge.10) then
         symblk = flcod - (flcod/10)*10
         fopcod = (flcod-symblk) / 10
      else
         fopcod = flcod
         symblk = 0
      endif
c
c     # check fopcod for valid value.
c
      if ((fopcod.gt.lstfop).or.(fopcod.lt.1)) then
         filerr=-5
         goto 1000
      endif
c
c     # check symblk for valid value.
c
      if ((symblk.gt.maxsym).or.(symblk.lt.0)) then
         filerr=-7
         goto 1000
      endif
c
c
c     # initialize state (assume success).
c
      filerr = 0
      syserr = 0
      ldio   = .false.
c
c     # address to begin filling retrieved data into c(:).
c
      index=0
c
c     # rewind file.
c
      rewind (unit)
c
c     # read initial line for parsing.
c
      read (unit,5001,iostat=syserr) iversn, inname
      if (syserr.ne.0) goto 1100
5001  format(i2,a6)
      call allcap(inname)
c
c     # (optional checks on file version number go here)
c
      if ( (iversn .gt. lstver)
     & .or. (inname .ne. flname) ) then
         filerr=-4
         goto 1000
      endif
c
c     # read field 1, header.
c
      read (unit,fmt=fldfmt,iostat=syserr) inline
      if (syserr.ne.0) goto 1100
      call allcap(inline)
      if (inline.ne.field(1)) then
         filerr=-14
         goto 1000
      endif
c
c     # read in title information.
c
      read (unit,*,iostat=syserr) ntitle
      if (syserr.ne.0) goto 1100
      if (ntitle.gt.mxtitl) then
         filerr=-2
         goto 1000
      endif
      do 100 i=1,ntitle
         read (unit,5002,iostat=syserr) titles(i)
100   continue
      if (syserr.ne.0) goto 1100
5002  format (a80)
c
c     # read in symmetry information (dimensions).
c
      read (unit,*,iostat=syserr) nsym
      if (syserr.ne.0) goto 1100
      if (symblk.gt.nsym) then
         filerr=-8
         goto 1000
      endif
      if ((nsym.gt.maxsym).or.(nsym.lt.1)) then
         filerr=-3
         goto 1000
      endif
c
c     # read in nbfpsy() and nmopsy().
c
      read (unit,*,iostat=syserr) (nbfpsy(i),i=1,nsym),
     + (nmopsy(i),i=1,nsym)
      if (syserr.ne.0) goto 1100
c
c     # check for adequate space.
c
      totspc = 0
      do 140 isym=1,nsym
c
c        # check to see if block is to be read.
c
         if (symblk.eq.0.or.isym.eq.symblk) then
            if (fopcod.le.2) then
               totspc=totspc+nmopsy(isym)*nbfpsy(isym)
            else
               totspc=totspc+nmopsy(isym)
            endif
         endif
140   continue
      if ((totspc.gt.clen).and.(fopcod.gt.1)) then
         filerr=-1
         goto 1000
      endif
c
c     # clear out c vector before reading data.
c
      do 141 i=1,clen
         c(i) = zero
141   continue
c
c     # read in label information
c
      read (unit,5003,iostat=syserr)(labels(i),i=1,nsym)
      if (syserr.ne.0) goto 1100
5003  format(8(a4,1x))
c
c     # test for symmetry and title read only and return if necessary.
c
      if (fopcod.eq.1) return
c
c     # loop until to get to next field
c     # drop out of loop if field is reached or eof detected
c
130   read (unit,fldfmt,end=137,iostat=syserr) inline
      if (syserr.ne.0) goto 1100
      call allcap(inline)
      if (inline.eq.field(fopcod)) goto 135
      goto 130
c
c     # end of file detected before field is matched.
c
137   continue
      filerr=-6
      goto 1000
c
c     # end of until loop.
135   continue
c
c     # field is matched, begin reading in c() at position index.
c
      read (unit,5004,iostat=syserr) afmt
5004  format(a80)
      if (syserr.ne.0) goto 1100
c
c     # parse format designation.
c
      if (afmt.eq.'(*)') then
         ldio=.true.
      else
         ldio=.false.
      endif
c
c     # read as coefficients only if fopcod eq 2.
c
      if (fopcod.eq.2) then
c
c        # read in coefficients
c
         do 110 isym=1,nsym
            do 120 imo=1,nmopsy(isym)
c
c              # read block if it is to be read.
c
               if (isym.eq.symblk.or.symblk.eq.0) then
                  if (ldio) then
                     read (unit,*,iostat=syserr)
     +                (c(index+i),i=1,nbfpsy(isym))
c                  write(*,*) 'MOREAD1'
                  else
                     read (unit,fmt=afmt,iostat=syserr)
     +                (c(index+i),i=1,nbfpsy(isym))
c                  write(*,*) 'MOREAD2'
                  endif
                  if (syserr.ne.0) goto 1100
                  index=index+nbfpsy(isym)
c                  write(*,*) 'MOREAD3'
               else
c
c                 # skip values which are not to be read.
c
                  if (ldio) then
                     read (unit,*,iostat=syserr)
     +                (cdummy,i=1,nbfpsy(isym))
c                  write(*,*) 'MOREAD4'
                  else
                     read (unit,fmt=afmt,iostat=syserr)
     +                (cdummy,i=1,nbfpsy(isym))
c                  write(*,*) 'MOREAD5'
                  endif
                  if (syserr.ne.0) goto 1100
c                  write(*,*) 'MOREAD6'
               endif
120         continue
110      continue
c
      elseif ((fopcod.eq.3).or.(fopcod.eq.4)) then
c
c        # read in other fields for fopcod 3,4.
c
         do 150 isym=1,nsym
            if (symblk.eq.0.or.symblk.eq.isym) then
c
c              # read if block is to be read.
c
               if (ldio) then
                  read (unit,*,iostat=syserr)
     +             (c(index+i),i=1,nmopsy(isym))
c                  write(*,*) 'MOREAD7'
               else
                  read (unit,fmt=afmt,iostat=syserr)
     +             (c(index+i),i=1,nmopsy(isym))
c                  write(*,*) 'MOREAD8'
               endif
               if (syserr.ne.0) goto 1100
               index=index+nmopsy(isym)
c                  write(*,*) 'MOREAD9'
            else
c
c              # skip values.
c
               if (ldio) then
                  read (unit,*,iostat=syserr)
     +             (cdummy,i=1,nmopsy(isym))
c                  write(*,*) 'MOREAD10'
               else
                  read (unit,fmt=afmt,iostat=syserr)
     +             (cdummy,i=1,nmopsy(isym))
c                  write(*,*) 'MOREAD11'
               endif
               if (syserr.ne.0) goto 1100
c                  write(*,*) 'MOREAD12'
            endif
150      continue
      endif
c
      return
c
c     # error handling portion.
c
c     # system errors.
1100  continue
      call bummer('moread: system file error',syserr,wrnerr)
      filerr=1
      return
c
c     # file errors.
1000  continue
      return
c
      end
c deck mowrit
c *** this routine is incremental ***
      subroutine mowrit(
     & unit,   flcod,  filerr, syserr,
     & ntitle, titles, afmt,   nsym,
     & nbfpsy, nmopsy, labels, c  )
c
c subroutine: mowrit
c purpose   : universal interface write of mo coefficient file
c author    : eric stahlberg, osu chemistry (may 1989)
c
c version   : 2.1
c 17-apr-91 if();do 130 bug fixed. -eas
c
c variable definitions
c
c *** input variables ***
c unit  : unit number of mo file
c flcod: file operation code
c         10   = write header to new file
c         20   = append mo coefficients to file
c         30   = append orbital energies to file
c         40   = append orbital occupations to file
c          + symmetry block to be written (0 means all blocks)
c         flcod lt 10 implies all blocks and that fopcod
c ntitle: number of title cards to write
c titles: actual title cards
c nsym  : number of irreps in point group
c nbfpsy: number of basis functions per mo of each irrep
c nmopsy: number of mos per irrep
c labels: symmetry labels for each irrep
c c     : coefficient array space
c afmt  : format to write coefficients
c
c output values
c filerr : mowrit file error number
c syserr : mowrit system i/o file error number
c
       implicit logical(a-z)
      integer unit, ntitle, nsym, flcod, filerr, syserr
      integer nbfpsy(*), nmopsy(*)
      real*8 c(*)
      character*80 titles(*)
      character*4 labels(*)
      character*(*) afmt
c
c     # local variables.
c
      integer index, i, isym, imo, symblk, fopcod
      character*6 fname,field(4),fldfmt
c
      integer    iversn,  lstfop,  maxsym
      parameter (iversn=2,lstfop=4,maxsym=8)
c
c     # saved variables.
c
      integer    mxunit
      parameter (mxunit=99)
      integer prvfop(mxunit),prvblk(mxunit),nmblk(mxunit)
      save prvfop,prvblk,nmblk
c
c     # initialize state.
c
      syserr   = 0
      filerr   = 0
      fname    = 'mocoef'
      field(1) = 'header'
      field(2) = 'mocoef'
      field(3) = 'energy'
      field(4) = 'orbocc'
      fldfmt   = ' (a6) '
c
c     # parse out block number and real fopcod.
c
      if (flcod.ge.10) then
         symblk= flcod - (flcod/10)*10
         fopcod= flcod/10
      else
         fopcod = flcod
         symblk = 0
      endif
c
c     # check that passed fopcod is valid.
c
      if (fopcod.gt.lstfop.or.fopcod.lt.1) then
         filerr = -9
         goto 1000
      endif
c
c     # check that valid nsym is passed.
c
      if (nsym.gt.maxsym.or.nsym.lt.1) then
         filerr = -12
         goto 1000
      endif
c
c     # check that valid symblk exists.
c
      if (symblk.lt.0.or.symblk.gt.nsym) then
         filerr = -13
         goto 1000
      endif
c
c     # check to see if unit was used previously.
c
      if (prvfop(unit).lt.1.and.fopcod.ne.1) then
         filerr = -14
         goto 1000
      endif
c
c     # check that last block was written completely.
c
      if ((fopcod.eq.1).or.(prvfop(unit).ne.fopcod)) then
         if (prvblk(unit).ne.nmblk(unit)) then
            filerr = -10
            goto 1000
         endif
      endif
c
c     # check that blocks are in order for fields other than header.
c
      if (fopcod.ne.1) then
         if (symblk.eq.1.or.symblk.eq.0) then
            if (prvblk(unit).ne.nmblk(unit)) then
               filerr = -11
               goto 1000
            endif
         elseif  (prvblk(unit).ne.(symblk-1)) then
            filerr = -11
            goto 1000
         endif
      endif
c
c     # save requests status.
c
      prvfop(unit) = fopcod
      nmblk(unit)  = nsym
c
c     # create a new file iff fopcod eq 1.
c
      if (fopcod.eq.1) then
c
c        # write file header.
c
         write (unit,5001) iversn,fname
5001     format (i2,a6)
c
c        # write header field delimiter.
c
         write (unit,fmt=fldfmt) field(1)
c
c        # write number of title cards.
c
         write (unit,*) ntitle
c
c        # write title cards.
c
         do 110 i = 1, ntitle
            write (unit,5002) titles(i)
110      continue
5002     format (a80)
c
c        # write number of irreps.
c
         write (unit,*) nsym
c
c        # write nmopsy and nbfpsy.
c
         write (unit,*) (nbfpsy(i),i=1,nsym),(nmopsy(i),i=1,nsym)
c
c        # write labels.
c
         write (unit,5003) (labels(i),i=1,nsym)
5003     format (8(a4,1x))
c
c        # set prvblk = nsym.
c
         prvblk(unit) = nsym
c
c        # write mo coefficient field.
c
      elseif (fopcod.eq.2) then
         if (symblk.eq.0.or.symblk.eq.1) then
c
c           # write field header.
c
            write (unit,fmt=fldfmt) field(2)
c
c           # write format.
c
            write (unit,5004) afmt
5004        format (a)
         endif
c
c        # begin writing coefficients.
c
         index=0
         do 120 isym=1,nsym
            if (isym.eq.symblk.or.symblk.eq.0) then
               do 130 imo=1,nmopsy(isym)
                  if (afmt.eq.'(*)') then
                     write(unit,*) (c(index+i),i=1,nbfpsy(isym))
                  else
                     write (unit,fmt=afmt) (c(index+i),i=1,nbfpsy(isym))
                  endif
                  index=index+nbfpsy(isym)
130            continue
               prvblk(unit)=isym
            endif
120      continue
c
      elseif ((fopcod.eq.3).or.(fopcod.eq.4)) then
c
c        # write out new field.
c
         index = 0
         if (symblk.eq.0.or.symblk.eq.1) then
            write (unit,fmt=fldfmt) field(fopcod)
            write (unit,5004) afmt
         endif
         do 160 isym=1,nsym
            if (symblk.eq.0.or.symblk.eq.isym) then
               if (afmt.eq.'(*)') then
                  write(unit,*) (c(index+i),i=1,nmopsy(isym))
               else
                  write (unit,fmt=afmt) (c(index+i),i=1,nmopsy(isym))
               endif
               index=index+nmopsy(isym)
               prvblk(unit)=isym
            endif
160      continue
      endif
c
1000  continue
      return
      end
c deck ufvin
      subroutine ufvin(
     & unit,   reclen, vector, ciuvfl,
     & lenci,  ninfo,  info,   nenrgy,
     & ietype, energy, ntitle, title,
     & heap,   totfre, mxtitl, mxinfo,
     & mxnrgy  )
c
c  this subroutine reads the CI vector
c  in a standard sequential unformatted vector format,
c  and writes the vector to the (local) da file using putvec().
c
c  24-oct-91 energy(*) read fixed. -rls
c
c  author: eric stahlberg
c  date  : 10-oct-90
c  version: 2.0
c
       implicit logical(a-z)
      integer unit, vector, lenci, totfre, ciuvfl, reclen
      integer mxtitl, mxinfo, mxnrgy
      real*8 heap(*)
      integer ntitle, nenrgy, ninfo
      character*80 title(*)
      real*8 info(*), energy(*)
      integer ietype(*)
c
c     # local variables.
      integer blksiz, recamt, icsf, icode, mxread, i
      integer versn
c
c     # curver is the current version number compatible with this
c     # reading routine.
c
      integer    curver
      parameter (curver = 1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # read the first file header.
      read (ciuvfl) versn, blksiz, lenci
      if ( versn .ne. curver ) then
         call bummer('ufvin: disallowed version number',versn,faterr)
      elseif ( blksiz .gt. totfre ) then
         call bummer('ufvin: blocksize too large',blksiz,faterr)
      endif
c     # general information.
      read (ciuvfl ) ninfo, nenrgy, ntitle
c
      mxread = min (ninfo, mxinfo)
      read (ciuvfl ) (info(i), i=1,mxread)
c
      mxread = min (ntitle, mxtitl)
      read (ciuvfl ) (title(i), i=1,mxread)
c
      if (mxnrgy .le. 0 ) then
         read (ciuvfl)
      else
         mxread = min (nenrgy, mxnrgy)
         read (ciuvfl ) (ietype(i),i=1,mxread),(energy(i),i=1,mxread)
      endif
c
c     # list of coefficients.
      if (blksiz.ne.reclen) then
         icode = 2
      else
         icode = 0
      endif
      do 10 icsf = 1, lenci, blksiz
         recamt = min (blksiz, (lenci - icsf + 1))
         call seqrbf( ciuvfl, heap, recamt)
         call putvec(
     &    unit,   reclen, vector, lenci,
     &    heap,   icsf,   recamt, icode )
10    continue
c
      return
      end
c deck ufvout
      subroutine ufvout(
     & unit,   reclen, ivec,   ciuvfl,
     & lenci,  ninfo,  info,   nenrgy,
     & ietype, energy, ntitle, title,
     & heap,   totfre, loglen  )

ctm
c     lenci : number of elements to be written to sequential file
c     loglen: logical length of ci vector
ctm
c
c  this subroutine writes the specified CI vector
c  in a standard sequential unformatted vector format.
c
c  author: eric stahlberg
c  date  : 10-oct-90
c  version: 2.0
c
       implicit logical(a-z)
      integer unit, ivec, lenci, totfre, ciuvfl, reclen
      integer loglen
      real*8 heap(*)
      integer ntitle, nenrgy, ninfo
      character*80 title(*)
      real*8 info(*), energy(*)
      integer ietype(*)
c
c     # local variables.
c
      integer i, blksiz, recamt, icsf
      integer versn
c
c     # curver is the current version number for this writing routine.
c
      integer    curver
      parameter (curver = 1)
c
      blksiz = min( totfre, lenci, reclen )
      versn = curver
c
c     # file header.
      write(ciuvfl) versn, blksiz, lenci
c     # general information.
      write(ciuvfl ) ninfo, nenrgy, ntitle
      write(ciuvfl) (info(i),i=1,ninfo)
      write(ciuvfl) (title(i),i=1,ntitle)
      write(ciuvfl) (ietype(i),i=1,nenrgy),(energy(i),i=1,nenrgy)
c     # list of coefficients.
      do 10 icsf = 1, lenci, blksiz
         recamt = min (blksiz, (lenci - icsf + 1))
         call getvec( unit, reclen, ivec, loglen, heap, icsf, recamt )
         call seqwbf(ciuvfl,heap,recamt)
10    continue
c
      return
      end
c deck writex
      subroutine writex( a, n )
c writex:
c this routine writes out a formatted real*8 vector
c a : the vector to be listed
c n : the number of elements to list
c
c  10-dec-90 iwritx() entry added for unit initialization. -rls
c
      real*8 a(*)
      integer n,m, iunit
c
      integer i,j, nlist
      save nlist
      data nlist /6/
c
      m=1
      do j =1, (n/100)
      write(nlist,6011)m,m+99
      write(nlist,6010) (a(i),i=m,m+99)
      m=m+100
      enddo
      write(nlist,6011)m,n
      write(nlist,6010) (a(i),i=m,n)
      return
c
      entry iwritx( iunit )
      nlist = iunit
c
      return
6010  format(1x,10f10.6)
6011  format(1x, 4('----------'),'elements ',i4,' to ', i4,
     .           4('----------'))
      end
c deck getlcl
      subroutine getlcl( lcored, lcore, ierr )
c
c  get the lcore parameter from the command line
c  for workspace allocation.
c
c  input:
c  lcored = default value.  this is used if no replacement value can
c           be found on the command line.
c
c  output:
c  lcore = amount of workspace to be allocated in the calling program.
c  ierr  = return code. 0 for normal return, <>0 for error.
c
c  19-oct-01 f2kcli interface added. -rls
c  24-apr-92 variable format used in the internal read statement.
c            this is a workaround for an ibm rs6000 bug. -rls
c  13-mar-91 posix version added. -rls
c  30-nov-90 vms code added. -rls
c  29-nov-90 ierr returned. explicit error processing removed. -rls
c  20-mar-90 written by ron shepard.
c
       implicit logical(a-z)
      integer lcored, lcore, ierr
c
*@ifdef f2kcli
*c
*c     # based on the (proposed) fortran 2000 command line interface
*c     #
*c     # usage:  % program_name [ -M lcore ] [ other ignored options]
*c     #
*c     # where lcore is a positive integer used for workspace allocation
*c     # in the calling program.
*
*      integer   lenc
*      parameter(lenc=9)
*      character*(lenc) carg
*c
*      integer i, last
*      character*4 cfmt
*c
*      logical  streq
*      external streq
*c
*      integer  command_argument_count
*      external command_argument_count
*c
*      external get_command_argument
*c
*      ierr  = 0
*      lcore = lcored
*      cfmt  = '(i9)'
*c
*c     # loop over the arguments, looking for -M nnnnnn
*c     # if -M is the last argument, then just ignore it.
*c
*      do 10 i = 1, ( command_argument_count() - 1 )
*c
*         call get_command_argument( i, carg, last, ierr )
*         if ( ierr .ne. 0 ) return
*c
*c        # string comparisons are folded to upper case.
*         if ( streq( '-m', carg) ) then
*c
*c           # found the option, now get the value.
*c
*            call get_command_argument( (i+1), carg, last, ierr )
*            if ( ierr .ne. 0 ) return
*c
*c           # put the value of last into the cfmt variable.
*c           # added 24-apr-92. rs6000 does not allow fixed '(i9)'. -rls
*            write( cfmt(3:3), '(i1)' ) last
*c
*c           # in this version, the argument must be a nonnegative
*c           # integer;  character strings may be parsed in the future.
*c           # fmt should agree with lenc.
*c
*            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
*            if ( ierr .ne. 0 ) then
*               return
*            else
*c              # lcore was integer, check for validity.
*               if ( lcore .le. 0 ) then
*c                 # special meaning may be given to negative values
*c                 # later.  for now, treat this as an error.
*                  ierr = -2
*               endif
*               return
*            endif
*         endif
*10    continue
*c
*@elif defined  (T3E64) || defined ( T3D ) || defined (posix)
*      integer   lenc
*      parameter(lenc=9)
*      character*(lenc) carg
*c
*      integer i, last
*      character*4 cfmt
*c
*      integer  ipxfargc, fstrlen
*      logical                    streq
*      external ipxfargc, fstrlen, streq
*c
*      ierr  = 0
*      lcore = lcored
*      cfmt  = '(i9)'
*c
*c     # loop over the arguments, looking for -M nnnnnn
*c     # if -M is the last argument, then just ignore it.
*c
*      do 10 i = 1, ( ipxfargc() - 1 )
*c
*         call pxfgetarg( i, carg, last, ierr )
*         if ( ierr .ne. 0 ) then
*            return
*         endif
*c
*c        # string comparisons are folded to upper case.
*         if ( streq( '-m', carg) ) then
*c
*c           # found the option, now get the value.
*c
*            call pxfgetarg( (i+1), carg, last, ierr )
*            if ( ierr .ne. 0 ) then
*               return
*            endif
*c
*c           # put the value of last into the cfmt variable.
*            write( cfmt(3:3), '(i1)' ) last
*c
*c
*c           # in this version, the argument must be a nonnegative
*c           # integer;  character strings may be parsed in the future.
*c           # fmt should agree with lenc.
*c
*            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
*            if ( ierr .ne. 0 ) then
*               return
*            else
*c              # lcore was integer, check for validity.
*               if ( lcore .le. 0 ) then
*c                 # special meaning may be given to negative values
*c                 # later.  for now, treat this as an error.
*                  ierr = -2
*               endif
*               return
*            endif
*         endif
*10    continue
*c
*@elif defined  hp
*c
*c     # based on the "generic" unix version.
*c     # hp fortran does not add "_" suffix to external symbols
*c     # so they must be added manually here.
*c     #
*c     # usage:  % program_name [ -M lcore ] [ other ignored options]
*c     #
*c     # where lcore is a positive integer used for workspace allocation
*c     # in calling program.
*
*      integer   lenc
*      parameter(lenc=9)
*      character*(lenc) carg
*c
*      integer i, last
*      character*4 cfmt
*c
*      integer  iargc_, fstrlen
*      logical                 streq
*      external iargc_, fstrlen, streq
*c
*      ierr  = 0
*      lcore = lcored
*      cfmt  = '(i9)'
*c
*c     # loop over the arguments, looking for -M nnnnnn
*c     # if -M is the last argument, then just ignore it.
*c
*      do 10 i = 1, ( iargc_() - 1 )
*c
*c        # note: last=getarg() does not work on all machines.
*c        #       this is a lowest-common-denominator approach.
*         call getarg_( i, carg )
*c
*c        # string comparisons are folded to upper case.
*         if ( streq( '-m', carg) ) then
*c
*c           # found the option, now get the value.
*c
*            call getarg_( (i+1), carg )
*            last = fstrlen( carg )
*c
*c           # put the value of last into the cfmt variable.
*c           # added 24-apr-92. rs6000 does not allow fixed '(i9)'. -rls
*            write( cfmt(3:3), '(i1)' ) last
*c
*c
*c           # in this version, the argument must be a nonnegative
*c           # integer;  character strings may be parsed in the future.
*c           # fmt should agree with lenc.
*c
*            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
*            if ( ierr .ne. 0 ) then
*               return
*            else
*c              # lcore was integer, check for validity.
*               if ( lcore .le. 0 ) then
*c                 # special meaning may be given to negative values
*c                 # later.  for now, treat this as an error.
*                  ierr = -2
*               endif
*               return
*            endif
*         endif
*10    continue
*c
*
*@elif defined  unix
c
c     # "generic"  version.
c     # this isn't pretty, but it works almost everywhere.
c     #
c     # usage:  % program_name [ -M lcore ] [ other ignored options]
c     #
c     # where lcore is a positive integer used for workspace allocation
c     # in calling program.

      integer   lenc
      parameter(lenc=9)
      character*(lenc) carg
c
      integer i, last
      character*4 cfmt
c
      integer  iargc, fstrlen
      logical                 streq
      external iargc, fstrlen, streq
c
      ierr  = 0
      lcore = lcored
      cfmt  = '(i9)'
c
c     # loop over the arguments, looking for -M nnnnnn
c     # if -M is the last argument, then just ignore it.
c
      do 10 i = 1, ( iargc() - 1 )
c
c        # note: last=getarg() does not work on all machines.
c        #       this is a lowest-common-denominator approach.
         call getarg( i, carg )
c
c        # string comparisons are folded to upper case.
         if ( streq( '-m', carg) ) then
c
c           # found the option, now get the value.
c
            call getarg( (i+1), carg )
            last = fstrlen( carg )
c
c           # put the value of last into the cfmt variable.
c           # added 24-apr-92. rs6000 does not allow fixed '(i9)'. -rls
            write( cfmt(3:3), '(i1)' ) last
c
c           # in this version, the argument must be a nonnegative
c           # integer;  character strings may be parsed in the future.
c           # fmt should agree with lenc.
c
            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
            if ( ierr .ne. 0 ) then
               return
            else
c              # lcore was integer, check for validity.
               if ( lcore .le. 0 ) then
c                 # special meaning may be given to negative values
c                 # later.  for now, treat this as an error.
                  ierr = -2
               endif
               return
            endif
         endif
10    continue
c
*@elif defined  vax
*c
*c     # vax vms version.
*c     #
*c     # usage: $ program_name [ /m=lcore ]
*c     #
*c     # lcored is a fixed array size in the calling program,
*c     # and is the maximum allowed value.
*c
*      integer   lenc
*      parameter(lenc=9+3)
*      character*(lenc) carg
*      integer cmdlen
*c
*      logical  streq
*      external streq
*c
*      ierr  = 0
*      lcore = lcored
*c
*c     # get the foreign command line.
*c     # program must be defined as a foreign command with
*c     # the vms statement:
*c     # $ prog:==$device:[directory]prog.exe
*c
*      call lib$get_foreign( carg,, cmdlen, )
*c
*c     # string comparisons are folded to upper case.
*c
*      if ( cmdlen .eq. 0 ) then
*c        # no options on the command line.
*         return
*      elseif ( cmdlen .le. 3 ) then
*c        # /m= must be followed by an integer value.
*         ierr = -2
*         return
*      elseif ( streq( '/m=', carg(1:3) ) ) then
*c
*c        # in this version, the argument must be a nonnegative
*c        # integer;  character strings may be parsed in the future.
*c        # fmt should agree with lenc.
*c
*         read( carg(4:cmdlen), fmt='(i9)', iostat=ierr ) lcore
*         if ( ierr .ne. 0 ) then
*            return
*         else
*c           # lcore was integer, check for validity.
*            if ( (lcore .le. 0) .or. (lcore .gt. lcored) ) then
*c              # special meaning may be given to negative values
*c              # later.  for now, treat this as an error.
*               ierr = -3
*            endif
*            return
*         endif
*      else
*c        # unrecognized option.
*         ierr = -4
*         return
*      endif
*@else
*c
*c     # don't know how to do anything else, so punt.
*c
*      ierr  = 0
*      lcore = lcored
*c
*@endif
c
      return
      end


c
      subroutine getfvd( lcored, lcore, ierr )
c
c  get the lcore parameter from the command line
c  for workspace allocation.
c
c  input:
c  lcored = default value.  this is used if no replacement value can
c           be found on the command line.
c
c  output:
c  lcore = amount of workspace to be allocated in the calling program.
c  ierr  = return code. 0 for normal return, <>0 for error.
c
c  24-apr-92 variable format used in the internal read statement.
c            this is a workaround for an ibm rs6000 bug. -rls
c  13-mar-91 posix version added. -rls
c  30-nov-90 vms code added. -rls
c  29-nov-90 ierr returned. explicit error processing removed. -rls
c  20-mar-90 written by ron shepard.
c
       implicit logical(a-z)
      integer lcored, lcore, ierr
c
*@if defined(T3E64) || defined ( T3D)
*c
*      integer   lenc
*      parameter(lenc=9)
*      character*(lenc) carg
*c
*      integer i, last
*      character*4 cfmt
*c
*      integer  ipxfargc, fstrlen
*      logical                    streq
*      external ipxfargc, fstrlen, streq
*c
*      ierr  = 0
*      lcore = lcored
*      cfmt  = '(i9)'
*c
*c     # loop over the arguments, looking for -V nnnnnn
*c     # if -V is the last argument, then just ignore it.
*c
*      do 10 i = 1, (ipxfargc() - 1 )
*c
*         call pxfgetarg( i, carg, last, ierr )
*         if ( ierr .ne. 0 ) then
*            return
*         endif
*c
*c        # string comparisons are folded to upper case.
*         if ( streq( '-v', carg) ) then
*c
*c           # found the option, now get the value.
*c
*            call pxfgetarg( (i+1), carg, last, ierr )
*            if ( ierr .ne. 0 ) then
*               return
*            endif
*c
*c           # put the value of last into the cfmt variable.
*            write( cfmt(3:3), '(i1)' ) last
*c
*c           # in this version, the argument must be a nonnegative
*c           # integer;  character strings may be parsed in the future.
*c           # fmt should agree with lenc.
*c
*            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
*            if ( ierr .ne. 0 ) then
*               return
*            else
*c              # lcore was integer, check for validity.
*               if ( lcore .lt. 0 ) then
*c                 # special meaning may be given to negative values
*c                 # later.  for now, treat this as an error.
*                  ierr = -2
*               endif
*               return
*            endif
*         endif
*10    continue
*@elif defined  hp
*c
*c     # "generic" unix version.
*c     # this isn't pretty, but it works almost everywhere.
*c     #
*c     # usage:  % program_name [ -V lcore ] [ other ignored options]
*c     #
*c     # where lcore is a positive integer used for workspace allocation
*c     # in calling program.
*
*      integer   lenc
*      parameter(lenc=9)
*      character*(lenc) carg
*c
*      integer i, last
*      character*4 cfmt
*c
*      integer  iargc_, fstrlen
*      logical                 streq
*      external iargc_, fstrlen, streq
*c
*      ierr  = 0
*      lcore = lcored
*      cfmt  = '(i9)'
*c
*c     # loop over the arguments, looking for -V nnnnnn
*c     # if -V is the last argument, then just ignore it.
*c
*      do 10 i = 1, ( iargc_() - 1 )
*c
*c        # note: last=getarg() does not work on all machines.
*c        #       this is a lowest-common-denominator approach.
*         call getarg_( i, carg )
*c
*c        # string comparisons are folded to upper case.
*         if ( streq( '-v', carg) ) then
*c
*c           # found the option, now get the value.
*c
*            call getarg_( (i+1), carg )
*            last = fstrlen( carg )
*c
*c           # put the value of last into the cfmt variable.
*c           # added 24-apr-92. rs6000 does not allow fixed '(i9)'. -rls
*            write( cfmt(3:3), '(i1)' ) last
*c
*c           # in this version, the argument must be a nonnegative
*c           # integer;  character strings may be parsed in the future.
*c           # fmt should agree with lenc.
*c
*            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
*            if ( ierr .ne. 0 ) then
*               return
*            else
*c              # lcore was integer, check for validity.
*               if ( lcore .lt. 0 ) then
*c                 # special meaning may be given to negative values
*c                 # later.  for now, treat this as an error.
*                  ierr = -2
*               endif
*               return
*            endif
*         endif
*10    continue
*c
*@elif defined  (unix and (not unicos))
*c
*c     # "generic" unix version.
*c     # this isn't pretty, but it works almost everywhere.
*c     #
*c     # usage:  % program_name [ -V lcore ] [ other ignored options]
*c     #
*c     # where lcore is a positive integer used for workspace allocation
*c     # in calling program.
*
*      integer   lenc
*      parameter(lenc=9)
*      character*(lenc) carg
*c
*      integer i, last
*      character*4 cfmt
*c
*      integer  iargc, fstrlen
*      logical                 streq
*      external iargc, fstrlen, streq
*c
*      ierr  = 0
*      lcore = lcored
*      cfmt  = '(i9)'
*c
*c     # loop over the arguments, looking for -V nnnnnn
*c     # if -V is the last argument, then just ignore it.
*c
*      do 10 i = 1, ( iargc() - 1 )
*c
*c        # note: last=getarg() does not work on all machines.
*c        #       this is a lowest-common-denominator approach.
*         call getarg( i, carg )
*c
*c        # string comparisons are folded to upper case.
*         if ( streq( '-v', carg) ) then
*c
*c           # found the option, now get the value.
*c
*            call getarg( (i+1), carg )
*            last = fstrlen( carg )
*c
*c           # put the value of last into the cfmt variable.
*c           # added 24-apr-92. rs6000 does not allow fixed '(i9)'. -rls
*            write( cfmt(3:3), '(i1)' ) last
*c
*c           # in this version, the argument must be a nonnegative
*c           # integer;  character strings may be parsed in the future.
*c           # fmt should agree with lenc.
*c
*            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
*            if ( ierr .ne. 0 ) then
*               return
*            else
*c              # lcore was integer, check for validity.
*               if ( lcore .lt. 0 ) then
*c                 # special meaning may be given to negative values
*c                 # later.  for now, treat this as an error.
*                  ierr = -2
*               endif
*               return
*            endif
*         endif
*10    continue
*c
*@elif defined  (unix and unicos)
*c
*c     # "generic" unix version.
*c     # this isn't pretty, but it works almost everywhere.
*c     #
*c     # usage:  % program_name [ -V lcore ] [ other ignored options]
*c     #
*c     # where lcore is a positive integer used for workspace allocation
*c     # in calling program.
*
*      integer   lenc
*      parameter(lenc=9)
*      character*(lenc) carg
*c
*      integer i, last
*      character*4 cfmt
*c
*      integer  iargc, fstrlen
*      logical                 streq
*      external iargc, fstrlen, streq
*c
*      ierr  = 0
*      lcore = lcored
*      cfmt  = '(i9)'
*c
*c     # loop over the arguments, looking for -V nnnnnn
*c     # if -V is the last argument, then just ignore it.
*c
*      do 10 i = 1, ( iargc() - 1 )
*c
*c        # note: last=getarg() does not work on all machines.
*c        #       this is a lowest-common-denominator approach.
*         call getarg( i, carg )
*c
*c        # string comparisons are folded to upper case.
*         if ( streq( '-v', carg) ) then
*c
*c           # found the option, now get the value.
*c
*            call getarg( (i+1), carg )
*            last = fstrlen( carg )
*c
*c           # put the value of last into the cfmt variable.
*c           # added 24-apr-92. rs6000 does not allow fixed '(i9)'. -rls
*            write( cfmt(3:3), '(i1)' ) last
*c
*c           # in this version, the argument must be a nonnegative
*c           # integer;  character strings may be parsed in the future.
*c           # fmt should agree with lenc.
*c
*            read( carg(1:last), fmt=cfmt, iostat=ierr ) lcore
*            if ( ierr .ne. 0 ) then
*               return
*            else
*c              # lcore was integer, check for validity.
*               if ( lcore .lt. 0 ) then
*c                 # special meaning may be given to negative values
*c                 # later.  for now, treat this as an error.
*                  ierr = -2
*               endif
*               return
*            endif
*         endif
*10    continue
*c
*@elif defined  vax
*c
*c     # vax vms version.
*c     #
*c     # usage: $ program_name [ /v=lcore ]
*c     #
*c     # lcored is a fixed array size in the calling program,
*c     # and is the maximum allowed value.
*c
*      integer   lenc
*      parameter(lenc=9+3)
*      character*(lenc) carg
*      integer cmdlen
*c
*      logical  streq
*      external streq
*c
*      ierr  = 0
*      lcore = lcored
*c
*c     # get the foreign command line.
*c     # program must be defined as a foreign command with
*c     # the vms statement:
*c     # $ prog:==$device:[directory]prog.exe
*c
*      call lib$get_foreign( carg,, cmdlen, )
*c
*c     # string comparisons are folded to upper case.
*c
*      if ( cmdlen .eq. 0 ) then
*c        # no options on the command line.
*         return
*      elseif ( cmdlen .le. 3 ) then
*c        # /v= must be followed by an integer value.
*         ierr = -2
*         return
*      elseif ( streq( '/v=', carg(1:3) ) ) then
*c
*c        # in this version, the argument must be a nonnegative
*c        # integer;  character strings may be parsed in the future.
*c        # fmt should agree with lenc.
*c
*         read( carg(4:cmdlen), fmt='(i9)', iostat=ierr ) lcore
*         if ( ierr .ne. 0 ) then
*            return
*         else
*c           # lcore was integer, check for validity.
*            if ( (lcore .lt. 0) .or. (lcore .gt. lcored) ) then
*c              # special meaning may be given to negative values
*c              # later.  for now, treat this as an error.
*               ierr = -3
*            endif
*            return
*         endif
*      else
*c        # unrecognized option.
*         ierr = -4
*         return
*      endif
*@else
c
c     # don't know how to do anything else, so punt.
c
      ierr  = 0
      lcore = lcored
c
*@endif
c
      return
      end


      subroutine trnfln( nunits, fnames )
c
c  perform any machine-specific filename translations.
c
c  input:  nunits = number of filenames.
c          fnames(1:nunits) = character filenames to be translated.
c
c  output: fnames(1:nunits) = updated filenames.
c
c  03-sep-91 unicos 6.0 interface added. -rls
c  13-mar-91 posix code added. -rls
c  01-dec-90 written by ron shepard.
c
       implicit logical(a-z)
      integer nunits
      character*(*) fnames(nunits)
c
      integer i
c
*@if defined(T3E64) || defined ( T3D) || defined (posix)
*c     # posix version.
*      integer vlen, ierr
*      character*255 envirn, value
*c
*c     # bummer error types.
*      integer   wrnerr,  nfterr,  faterr
*      parameter(wrnerr=0,nfterr=1,faterr=2)
*c
*      do 10 i = 1, nunits
*         if ( fnames(i) .ne. ' ' ) then
*            envirn = fnames(i)
*c
*c           # convert to upper case.
*            call allcap( envirn )
*c
*c           # look for a logical name translation.
*c           write(6,*) 'trnfln disabled ...'
*c           call pxfgetenv( envirn, 0, value, vlen, ierr )
*            vlen=0
*c           if ( ierr .ne. 0 ) then
*c              call bummer('trnfln: from pxfgetenv(), ierr=',
*c    &          ierr, faterr )
*c           endif
*c
*c           # if found, replace the filename with the value.
*            if ( vlen .ne. 0 ) fnames(i) = value(1:vlen)
*         endif
*10    continue
*@elif defined  unicos
*c     # with unicos 6.0, character variables can be used.  however, the
*c     # getenv() interface is still nonstandard. -rls
*      character*255 envirn, value
*      integer iname
*c
*      integer  getenv
*      external getenv
*c
*      do 10 i = 1, nunits
*         if ( fnames(i) .ne. ' ' ) then
*            envirn = fnames(i)
*c
*c           # convert to upper case.
*            call allcap( envirn )
*c
*c           # look for a logical name translation.
*c           # note: a separate statement is used to avoid
*c           #       "value"-related side-effects. -rls
*            iname = getenv( envirn, value )
*c
*c           # if found, replace the filename with the value.
*            if ( iname .ne. 0 ) fnames(i) = value
*         endif
*10    continue
*@elif defined  fujitsu
*c     # unix version for fujitsu vp.
*c     # envirn must be null-terminated in getenv() call.
*c     # 09-apr-92 (Ross Nobes, Roger Edberg) -rls
*      character*255 envirn, value
*      integer ntlen
*c
*      integer  fstrlen
*      external fstrlen
*c
*      do 10 i = 1, nunits
*         if ( fnames(i) .ne. ' ' ) then
*            envirn = fnames(i)
*c
*c           # convert to upper case.
*            call allcap( envirn )
*c
*c           # null-terminate. truncate if necessary.
*            ntlen = min( fstrlen( envirn ) + 1, len(envirn) )
*            envirn(ntlen:ntlen) = char(0)
*c
*c           # look for a logical name translation.
*            call getenv( envirn, value )
*c
*c           # if found, replace the filename with the value.
*            if ( value .ne. ' ') fnames(i) = value
*         endif
*10    continue
*@elif defined  unix
c     # generic bsd  version.
      character*255 envirn, value
c
      do 10 i = 1, nunits
         if ( fnames(i) .ne. ' ' ) then
            envirn = fnames(i)
c
c           # convert to upper case.
            call allcap( envirn )
c
c           # look for a logical name translation.
            call getenv( envirn, value )
c
c           # if found, replace the filename with the value.
            if ( value .ne. ' ') fnames(i) = value
         endif
10    continue
*@elif defined  vax
*c     # vax vms version.
*c     # logical name translations are done automatically by
*c     # the fortran library, so the explicit translation is
*c     # not necessary.
*c     # furthermore, the version number associated with the file
*c     # depends on the "status=" value used in the open statement,
*c     # so it is not possible to determine at this time the fully
*c     # qualified filename.
*@else
*c     # default case...just return.
*@endif
c
      return
      end
c deck h2ozzz
c deck h2oini
c deck h2oset
c deck h2opsh
c deck h2opop
c deck h2omtr
c deck h2ofin
      subroutine h2ozzz( hih2o, ierr, nlist )
c
c  this set of entry points computes high-water marks
c  for workspace allocation.
c
c  the use of this set of entry points eliminates the need for
c  maintaining global variables which are updated in each subroutine
c  which performs memory allocation from a workspace array.
c
c  the correct calling sequence is:
c
c  call h2oini( 0 )             # initialize the stacks. [only once
c                               # for each main branch.]
c  do n times...
c     chunk = ...               # allocate workspace, and update the
c     call h2oset( chunk)       # high-water mark for each chunk value.
c  enddo
c  do n times...
c     chunk = ...               # allocate new workspace, and update the
c     ...                       # high-water mark.
c     call h2opsh( chunk )      # prepare for a subroutine call.
c     call sub(...)             # this routine computes high-water marks
c                               #   using this same scheme.
c     call h2opop               # update high-water marks again after
c                               #   accounting for sub() usage.
c  enddo
c  call h2omtr( hih2o )         # returns hih2o, the high-water mark.
c  call h2ofin( ierr )          # check for correct calling sequence.
c
c  note in particular that repeated calls to h2oset() should be used
c  to determine high-water marks for different workspace allocations
c  relative to the same initial location.
c
c  the h2opsh/h2opop sequence should be used when a lower level
c  subroutine is to perform additional workspace allocations
c  relative to a new location.  The total number of h2opsh() calls
c  should be the same as the total number h2opop calls.
c
c  14-dec-90 added to colib for general use. -rls
c  19-nov-90 substantially revised, actual memory allocation removed,
c            and incorporated into argos(). -rls
c  01-jun-83 (approx.) memory allocation version written by ron shepard.
c
       implicit logical(a-z)
      integer hih2o, ierr, nlist
c
      integer    nlevmx
      parameter( nlevmx=100 )
c
      integer nevent, level, nerror, print
      save    nevent, level, nerror, print
c
c     # high(*) is used to store old high-water marks.
c     # branch(*) is used to store active branch lengths.
c
      integer high(nlevmx), branch(nlevmx)
      save    high,         branch
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      data nevent / 0 /
      data level  / 1 /
      data nerror / 1 /
      data print  / 6 /
      data branch / nlevmx*0 /
      data high   / nlevmx*0 /
c
c     # h2ozzz() is just to set up the entry points and
c     # encapsulate the data.
      call bummer('h2ozzz() called', 0, wrnerr )
      nerror = nerror + 1
      return
c
c**********************************************************************
      entry h2oini( nlist )
c**********************************************************************
c
c     # initialize the stack to begin computing high-water marks.
c     #
c     # input: nlist = output unit for writing debugging info.
c     #              = 0 operate silently.  this is the usual mode
c     #                  of operation after the calling sequence is
c     #                  debugged.
c
      nevent        = nevent + 1
      level         = 1
      branch(level) = 0
      high(level)   = 0
      nerror        = 0
      print         = nlist
      if ( print .ne. 0 ) then
         write(print,6010) nevent, nlevmx
      endif
      return
c
c**********************************************************************
      entry h2oset( hih2o )
c**********************************************************************
c
c     # set a new high-water mark at the current stack level.
c     #
c     # input: hih2o = workspace being allocated at the current stack
c     #                level.
c     #
c
      nevent = nevent + 1
      if ( print .ne. 0 ) then
         write(print,6020) nevent, level, high(level),
     &    branch(level), hih2o
      endif
      high(level)    = max( high(level), hih2o )
      branch(level)  = hih2o
      return
c
c**********************************************************************
      entry h2opsh( hih2o )
c**********************************************************************
c
c     # set a new high-water mark at the current stack level and
c     # push to a new stack level.
c     #
c
      nevent = nevent + 1
      if ( print .ne. 0 ) then
         write(print,6030) nevent, level, high(level),
     &    branch(level), hih2o
      endif
      if ( level .ge. nlevmx ) then
         call bummer('h2opsh: nlevel exceeded, nevent=',
     &    nevent, wrnerr )
         nerror = nerror + 1
         return
      endif
c
      high(level)   = max( high(level), hih2o )
      branch(level) = hih2o
      level         = level + 1
      branch(level) = 0
      high(level)   = 0
      return
c
c**********************************************************************
      entry h2opop
c**********************************************************************
c
c     # pop back to the previous level, updating the cumulative
c     # high-water mark.
c     #
c
      nevent = nevent + 1
      if ( print .ne. 0 ) then
         write(print,6040) nevent, level,
     &    branch(level-1), branch(level), high(level-1)
      endif
      if ( level .eq. 1 ) then
         call bummer('h2opop: level=1, nevent=', nevent,  wrnerr )
         nerror = nerror + 1
         return
      endif
      branch(level-1) = branch(level-1) + high(level)
      high(level-1)   = max( high(level-1), branch(level-1) )
      branch(level)   = 0
      high(level)     = 0
      branch(level-1) = 0
      level           = level - 1
      return
c
c**********************************************************************
      entry h2omtr( hih2o )
c**********************************************************************
c
c     # read the meter at the current stack level.
c     #
c     # output: hih2o = high-water mark at the current level.
c
      nevent = nevent + 1
      if ( print .ne. 0 ) then
         write(print,6050) nevent, level, branch(level), high(level)
      endif
      hih2o  = high(level)
      return
c
c**********************************************************************
      entry h2ofin( ierr )
c**********************************************************************
c
c     # make sure that upon completion, the stack level is correct
c     # and that no errors were generated during the previous calls.
c     #
c     # output: ierr = 0 for normal return.
c     #              = nerror if nerror errors were detected since the
c     #                       last h2oini() call.
c     #
c     # this entry should only be called at the lowest level; otherwise
c     # the stack pointer will not be the expected value, and this is
c     # recorded as an error.
c
      nevent = nevent + 1
      if ( print .ne. 0 ) then
         write(print,6060) nevent, level, nerror
      endif
      if ( level .ne. 1 ) then
         nerror = nerror + 1
      endif
c
      ierr = nerror
c
      return
6010  format(' h2oini: nevent=',i6,
     & ' high-water mark stack initialized with nlevmx=',i4)
6020  format(' h2oset: nevent=',i6,' level=',i2,' high(level)=',i9,
     & ' branch(level)=',i9,' hih2o=',i9)
6030  format(' h2opsh: nevent=',i6,' level=',i2,' high(level)=',i9,
     & ' branch(level)=',i9,' hih2o=',i9)
6040  format(' h2opop: nevent=',i6,' level=',i2,' branch(level-1)=',i9,
     & ' branch(level)=',i9,' high(level-1)=',i9)
6050  format(' h2omtr: nevent=',i6,' level=',i2,' branch(level)=',i9,
     & ' high(level)=',i9)
6060  format(' h2ofin: nevent=',i6,' level=',i2,' nerror=',i2)
      end
c deck adazhd
      subroutine adazhd( iunit, fname, len, dispos,
     & acctyp, next, null, irec, reqnum, buffer, ierr)
c
c  the entry points in this routine are for manipulating asynchronous
c  direct access files.
c  the correct calling sequence is:
c
c  call adaopn(...next,null,...)
c  reqnum = null
c  do
c     if ( reqnum .ne. null ) call adaw8t( ...,reqnum,...)
c     # ...wait until buffer(*) is free, then define its contents...
c     call adawrt(...buffer,...,reqnum...)
c  enddo
c  call adaw8t( ...,reqnum,...) # wait until the record is defined.
c  call adaclw()                # close for writing.(see below)
c  call adaopr()                # open for random read/write.(see below)
c  call adard(...reqnum...)     # start the first read.
c  do
c     if ( reqnum .ne. null ) call adaw8t( ...,reqnum,...)
c     # ...wait until buffer(*) is ready. then use the contents.
c     call adard(...reqnum...)  # start the next read operation.
c  enddo
c  call adacls(...)
c
c  this i/o model allows queuing of multiple i/o operations, provided
c  the buffers are not modified after they are written until an adaw8t()
c  call, and they are not accessed after an adard() call until they are
c  defined by an adaw8t() call.  for maximum flexibility, adaw8t() uses
c  both the record number and the i/o request number.
c
c  if the file is opened for sequential writes, then all records should
c  be written before any are read.  if the file is open for random
c  writes, then reads and writes may be intersperced, provided each
c  individual record is written before it is read.
c
c  note that the i/o request number, reqnum, is not necessarily unique.
c  it is used to identify pending i/o requests.  it is usually a pointer
c  into one of only a few operating system structures associated with
c  the i/o device controller.  reqnum could also be associated with a
c  channel number in a message-passing parallel-cpu environment.
c  the value of reqnum should not be modified by the calling program.
c
c  note that record pointers are not assumed to follow fortran
c  conventions.  in particular, the calling program should not assume
c  that record numbers are incremented by 1 on each successive call.
c  returned record pointers should be stored by the calling program
c  after calls to adawrt() and used to identify the needed records
c  later when calling adard().
c
c  08-oct-90 (columbus day) written by ron shepard.
c
       implicit logical(a-z)
c     # disposition and access type parameters for async da files.
      integer    delete,     keep
      parameter( delete = 0, keep   = 1 )
      integer    seqwrt,     random
      parameter( seqwrt = 0, random = 1 )
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer iunit, len, dispos, acctyp, next, null, irec, reqnum, ierr
      character*(*) fname
      real*8 buffer(len)
c
      integer iparm
c
c     # synchronous colib i/o routines use character arguments.
      character*8 chrdis, chracc
c
c     # nullp  = reqnum for null requests to be used with
c     #          synchronous i/o routines.
c     # reqrd  = reqnum for outstanding read.
c     # reqwrt = reqnum for outstanding write.
c     # first  = first record to be written using colib i/o routines.
      integer    nullp,   reqrd,   reqwrt,   first
      parameter( nullp=0, reqrd=1, reqwrt=2, first=1 )
c
c     # adaqhd() is just to set up entry points.  entry points are
c     # used to encapsulate local data. -rls
      call bummer('adazhd: illegal call',0,faterr)
      return
c
c**********************************************************************
c deck adaopn
      entry adaopn(
     & iunit,  fname, len,  dispos,
     & acctyp, next,  null, ierr )
c
c     # open an asynchronous direct access file.
c
c     # input: iunit = fortran unit number.
c     #        fname = file name.
c     #        len =  record length.
c     #        dispos = intended file disposition type.
c     #               = 0 for 'delete'. file will be deleted with the
c     #                   adacls() call.
c     #               = 1 for 'keep'.  file may be either kept or
c     #                   deleted with the adacls() call.
c     #        acctyp = file access type.
c     #               = 0 for sequential writes followed by random
c     #                   reads. all records must be written before
c     #                   the first read operation.
c     #               = 1 for intersperced random writes and random
c     #                   reads.  an individual record must be written
c     #                   before it is read, but otherwise writes and
c     #                   reads are in any order and may access records
c     #                   in any order.
c     #
c     # output: next = record pointer of the next (i.e. first)
c     #                available record.
c     #         null = recnum identifier for this unit that signifies
c     #                no outstanding or unsatisfied i/o requests.
c     #         ierr = return status. 0 for normal return.
c
*@ifdef future
*@else
c
c     # portable implementation: use synchronous colib i/o routines
c     # to implement a portable version of these routines.
c
      if ( dispos .eq. delete ) then
         chrdis = 'scratch'
      elseif ( dispos .eq. keep ) then
         chrdis = 'keep'
      else
         call bummer('adaopn: illegal dispos=',dispos,faterr)
      endif
      if ( acctyp .eq. seqwrt ) then
         chracc = 'seqwrt'
      elseif ( acctyp .eq. random ) then
         chracc = 'random'
      else
         call bummer('adaopn: illegal acctyp=',acctyp,faterr)
      endif
c
c     # openda() does not return next or ierr.
c     # assume fortran conventions for record numbers.
      next = first
      ierr = 0
      call openda( iunit, fname, len, chrdis, chracc )
c
c     # set the null record request.
      null = nullp
*@endif
c
      return
c**********************************************************************
c deck adawrt
      entry adawrt( iunit, irec, buffer, len, reqnum, ierr )
c
c     # write an asynchronous direct access record.
c     #
c     # input: iunit = fortran unit number.
c     #        irec  = record pointer.
c     #        buffer(1:len) = buffer to be written.  the contents
c     #                        should not be modified until a adaw8t()
c     #                        call.
c     #        len = buffer length. also must equal record length.
c     #
c     # output: irec   = pointer to the next available record.
c     #         reqnum = i/o request number for this record.
c     #         ierr   = return status. 0 for normal return.
c
*@ifdef future
*c  put real asynchronous routines here later.
*@else
c
c     # portable implementaton: use colib synchronous i/o routines.
c
c     # writda() does not return reqnum or ierr.
      reqnum = reqwrt
      ierr   = 0
      call writda( iunit, irec, buffer, len )
c     # writda() does not update irec.
c     # assume fortran record number conventions.
      irec   = irec + 1
*@endif
c
      return
c**********************************************************************
c deck adaclw
      entry adaclw( iunit, ierr )
c
c     # close asynchronous iunit for writing. prepare for random reads.
c
c     # this call is necessary only if the file was originally
c     # opened as seqwrt.
c
*@ifdef future
*c     # put async close statement here if necessary.
*@else
c     # portable implementation: use colib synchronous i/o operations.
      chrdis = 'keep'
      call closda( iunit, chrdis, iparm )
      ierr = 0
*@endif
c
      return
c**********************************************************************
c deck adaopr
      entry adaopr( iunit, fname, len, dispos, ierr )
c
c     # open iunit for asynchronous reads.
c
c     # this call is necessary only if the file was originally
c     # opened as seqwrt.
c
*@ifdef future
*c     # put async open statement here if necessary.
*@else
c     # portable implementation: use synchronous colib routines.
      if ( dispos .eq. delete ) then
         chrdis = 'scratch'
      elseif ( dispos .eq. keep ) then
         chrdis = 'keep'
      else
         call bummer('adaopr: illegal dispos=',dispos,faterr)
      endif
      chracc = 'random'
      call openda( iunit, fname, len, chrdis, chracc )
      ierr = 0
*@endif
c
      return
c**********************************************************************
c deck adard
      entry adard( iunit, irec, buffer, len, reqnum, ierr )
c
c     # read a record from an asynchronous direct access file.
c     #
c     # input: iunit = fortran unit number.
c     #        irec  = record pointer.
c     #        len   = buffer length. also must equal record length.
c     #
c     # output: irec = pointer to the next available record.
c     #         buffer(1:len) = output buffer to be filled.  the
c     #                         contents are not defined until an
c     #                         adaw8t() call.
c     #         reqnum = i/o request number for this record.
c     #         ierr = return status. 0 for normal return.
c
*@ifdef future
*@else
c
c     # portable implementation: use colib synchrohous i/o routines.
c
c     # readda() does not return reqnum or ierr.
      reqnum = reqrd
      ierr = 0
      call readda(  iunit, irec, buffer, len )
*@endif
c
      return
c**********************************************************************
c deck adaw8t
      entry adaw8t( iunit, reqnum, ierr )
c
c     # wait until pending i/o request number reqnum has completed.
c     #
c     # input: iunit  = fortran unit number.
c     #        reqnum = i/o operation request number.
c     #
c     # output: reqnum = nullp.  reset to nullp on return.
c     #         ierr = return status. 0 for normal return.
c
*@ifdef future
*@else
c
c     # portable implementation: use colib synchronous i/o routines.
c
      reqnum = nullp
      ierr   = 0
*@endif
c
      return
c**********************************************************************
c deck adacls
      entry adacls( iunit, dispos, ierr )
c
c     # close an asynchronous direct access file.
c     #
c     # input: iunit = fortran unit number.
c     #       dispos = file disposition.
c     #              = 0 to delete the file.  this is valid for files
c     #                  that were opened either scrtch or keep in
c     #                  adaopn().
c     #              = 1 to keep the file.   this is valid only for
c     #                  files that were opened as keep in adaopn().
c     #
c     # output: ierr = return status. 0 for normal return.
c
*@ifdef future
*@else
c
c     # portable implementation: use colib synchronous i/o routines.
c
c     # colib i/o routines use character arguments.
      if ( dispos .eq. delete ) then
         chrdis = 'delete'
      elseif ( dispos .eq. keep ) then
         chrdis = 'keep'
      else
         call bummer('adacls: illegal dispos=',dispos,faterr)
      endif
c
c     # closda() does not return ierr.
      ierr = 0
      call closda( iunit, chrdis, iparm )
*@endif
c
      return
      end
c deck echoin
      subroutine echoin( nin, nlist, ierr )
c
c  read and echo the input.
c
c  the input is assumed to be on 80-character records.
c
c  input:  nin   = input unit (assumed to be correctly positioned).
c          nlist = output listing unit unit (correctly positioned).
c
c  output: ierr  = return code.
c                = 0 for normal return. input file is positioned after
c                    the eof mark.  output file is positioned after the
c                    last output record.
c                > 0 for an external iostat error.
c
c  08-oct-90 (columbus day) written by ron shepard.
c
       implicit logical(a-z)
      integer nin, nlist, lenl, ierr
c
      integer ixerr
      character*80 line
c
      integer  fstrlen
      external fstrlen
c
      ierr = 0
      write( unit = nlist, fmt = 6010 )
10    continue
      read( unit = nin, fmt = '(a)', iostat = ixerr ) line
      if ( ixerr .eq. 0 ) then
c        # normal line.
c        # trim trailing blanks, and write line using fortran
c        # carriage control.
         lenl = max( 1, fstrlen(line) )
         write( unit = nlist, fmt = '(1x,a)' ) line(1:lenl)
         goto 10
      elseif ( ixerr .lt. 0 ) then
c        # eof on nin.  normal return.
         write( unit = nlist, fmt = 6010 )
         ierr = 0
      elseif ( ixerr .gt. 0 ) then
c        # error return.
         ierr = ixerr
      endif
c
      return
6010  format(1x,72('-') )
      end
c deck plblkc
      subroutine plblkc( title, z, nr, labr, ifmt, nlist )
c
c  print a lower-triangular packed matrix with row labels.
c  this version prints eight columns at a time with three formats.
c  parameter ncol and formats 10 and 1-3 should be modified to print
c  a different number of columns or to use different formats.
c
c  input:
c  title   = optional title.  only printed if (title.ne.' ').
c  z(*)    = matrix to be printed.
c  nr      = row and column dimension.
c  labr(*) = character*8 row and column labels.
c  ifmt    = format type (1:f, 2:e, 3:g).
c  nlist   = output unit nubmer.
c
c  18-oct-90 plblk() modified. title,labr(*) change. -rls
c  format statement assignment version 5-jun-87 (rls).
c  ron shepard 17-may-84.
c
       implicit logical(a-z)
      integer    ncol
      parameter (ncol=8)
10    format(/8x,8(6x,a8,1x))
1     format(1x,a8,8f15.8)
2     format(1x,a8,1p,8e15.6)
3     format(1x,a8,1p,8g15.6)
20    format(/10x,a)
c
      integer nr, ifmt, nlist
      character*(*) title, labr(*)
      real*8 z(*)
c
      integer fmtz, jlast, jstrt, j, ij0, i, j2, jt
      real*8    zero
      parameter(zero=0d0)
c
      if(ifmt.le.1)then
         assign 1 to fmtz
      elseif(ifmt.eq.2)then
         assign 2 to fmtz
      else
         assign 3 to fmtz
      endif
c
      if ( title .ne. ' ' ) write(nlist,20) title
c
      jlast = 0
      do 400 jstrt = 1, nr, ncol
         jlast = min( nr, jlast+ncol )
c
         write(nlist,10) ( labr(j), j = jstrt, jlast )
c
         ij0=(jstrt*(jstrt-1))/2
         do 300 i = jstrt, nr
            j2 = min( i, jlast )
c
c           # print the row if a nonzero element is found.
c
            do 100 j = jstrt, j2
               if ( z(ij0+j) .ne. zero ) then
                  write(nlist,fmtz) labr(i), ( z(ij0+jt), jt=jstrt,j2)
                  go to 101
               endif
100         continue
101         continue
            ij0 = ij0 + i
300      continue
400   continue
c
      return
      end
c deck plsbkc
      subroutine plsbkc(
     & gtitle, z, nblk, btitle, nrow, labr, ifmt, nlist )
c
c  print a lower-triangular subblocked matrix with row labels.
c
c  input:
c  gtitle    = general character title.
c  z(*)      = subblocked matrix to be printed.
c  nblk      = number of subblocks in the matrix z(*).
c  btitle(*) = specific character title of each subblock.
c  nrow(*)   = number of rows in each block.
c  labr(*)   = character*8 row labels.
c  ifmt      = format type (1:f, 2:e, 3:g).
c  nlist     = output unit number.
c
c  18-oct-90 plblks() modified. btitle(*),labr(*) change. -rls
c  ron shepard 17-may-84.
c
       implicit logical(a-z)
      integer nblk, ifmt, nlist
      character*(*) gtitle, btitle(1:nblk), labr(*)
      integer nrow(nblk)
      real*8 z(*)
c
      integer ipt, zpt, i, nr
c
      ipt = 1
      zpt = 1
      do 100 i = 1, nblk
         write(nlist,6010) gtitle, btitle(i), i
         nr = nrow(i)
         if ( nr .gt. 0 ) then
            call plblkc(' ', z(zpt), nr, labr(ipt), ifmt, nlist )
         endif
         ipt = ipt + nr
         zpt = zpt + ( nr * ( nr + 1 ) ) / 2
100   continue
      return
6010  format(/10x,a,1x,a,' block',i4)
      end
c deck prblkc
      subroutine prblkc(title, z, ldz, nr, nc, labr, labc, ifmt, nlist)
c
c  print a sub-block of a rectangular matrix with row and col labels.
c
c  this version prints eight columns at a time with 1 of 3 formats.
c  parameter ncol and formats 10 and 1-3 should be modified to print
c  a different number of columns or to use different formats.
c
c  input:
c  title   = optional character title.   only printed if (title.ne.' ').
c  z(*)    = matrix to be printed.
c  ldz     = effective leading dimension in the calling program..
c  nr      = number of rows to print.
c  nc      = column dimension.
c  labr(*) = character*8 row labels.
c  labc(*) = character*8 column labels.
c  ifmt    = format type (1:f, 2:e, 3:g).
c  nlist   = output unit number.
c
c  18-oct-90 prblk() modified. title,labr(*),labc(*) change. -rls
c  format statement assignment version 5-jun-87 (rls).
c  ron shepard 17-may-84.
c
       implicit logical(a-z)
      integer    ncol
      parameter (ncol=8)
10    format(/8x,8(6x,a8,1x))
1     format(1x,a8,8f15.8)
2     format(1x,a8,1p,8e15.6)
3     format(1x,a8,1p,8g15.6)
20    format(/10x,a)
c
      integer ldz, nr, nc, ifmt, nlist
      character*(*) title, labr(*), labc(*)
      real*8 z(ldz,nc)
c
      integer fmtz, i, j, jt, jstrt, jlast
      real*8     zero
      parameter (zero=0d0)
c
      if(ifmt.le.1)then
         assign 1 to fmtz
      elseif(ifmt.eq.2)then
         assign 2 to fmtz
      else
         assign 3 to fmtz
      endif
c
      if ( title .ne. ' ' ) write(nlist,20) title
c
      jlast = 0
      do 400 jstrt = 1, nc, ncol
         jlast = min( nc, jlast+ncol )
c
         write(nlist,10) ( labc(j), j = jstrt, jlast )
c
         do 300 i = 1, nr
c
c           # print the row if a nonzero element is found.
c
            do 100 j = jstrt, jlast
               if ( z(i,j) .ne. zero ) then
                  write(nlist,fmtz) labr(i), ( z(i,jt), jt=jstrt,jlast)
                  go to 101
               endif
100         continue
101         continue
300      continue
400   continue
c
      return
      end
c deck prsbkc
      subroutine prsbkc(
     & gtitle, z, nblk, btitle, nrow, ncol, labr, labc, ifmt, nlist )
c
c  print a rectangular-packed, diagonal subblocked matrix with labels.
c
c  input:
c  gtitle    = general character title.
c  z(*)      = subblocked rectangular matrix to be printed.
c  nblk      = number of subblocks in the matrix z(*).
c  btitle(*) = specific character title of each subblock.
c  nrow(*)   = number of rows in each subblock.
c  ncol(*)   = number of columns in each subblock.
c  labr      = character*8 row label.
c  labc      = character*8 column label.
c  ifmt      = format type (1:f, 2:e, 3:g).
c  nlist     = output unit number.
c
c  18-oct-90 prblks() modified. btitle(*),labr(*) change. -rls
c  ron shepard 17-may-84.
c
       implicit logical(a-z)
      integer nblk, ifmt, nlist
      character*(*) gtitle, btitle(*), labr(*), labc(*)
      integer nrow(nblk), ncol(nblk)
      real*8 z(*)
c
      integer ipt, jpt, zpt, i, nr, nc, nrnc
c
      ipt = 1
      jpt = 1
      zpt = 1
      do 100 i = 1, nblk
         write(nlist,6010) gtitle, btitle(i), i
         nr   = nrow(i)
         nc   = ncol(i)
         nrnc = nr * nc
         if ( nrnc .gt. 0 ) then
            call prblkc(' ', z(zpt), nr, nr, nc,labr(ipt), labc(jpt),
     &       ifmt, nlist )
         endif
         ipt = ipt + nr
         jpt = jpt + nc
         zpt = zpt + nrnc
100   continue
      return
6010  format(/10x,a,1x,a,' block',i4)
      end
c deck prvbkc
      subroutine prvbkc(
     & title, z, v, ldz, nr, nc, labr, labc, labv, ifmt, nlist )
c
c  print a rectangular matrix and a corresponding vector with labels.
c
c  this version prints 8 columns at a time with one of three formats.
c  parameter ncol and the appropriate formats should be modified to
c  print a different number of columns or to use different formats.
c
c  input:
c  title    = optional title.  only printed if (title.ne.' ').
c  z(*,*)   = matrix to be printed.
c  v(*)     = vector to be printed. v(i) corresponds to z(*,i).
c  ldz      = effective leading dimension of z(*).
c  nr       = number of rows to print.
c  nc       = column dimension.
c  labr(*)  = character row labels.
c  labc(*)  = character column labels.
c  labv     = character vector label.
c  ifmt     = format type (1:f, 2:e, 3:g).
c  nlist    = output unit number.
c
c  18-oct-90 prvblk() modified. title,labr(*),labc(*) change. -rls
c  format statement assignment version 5-jun-87 (rls).
c  ron shepard 17-may-84.
c
       implicit logical(a-z)
      integer    ncol
      parameter (ncol=8)
10    format(/8x,8(6x,a8,1x))
1     format(1x,a8,8f15.8)
2     format(1x,a8,1p,8e15.6)
3     format(1x,a8,1p,8g15.6)
11    format(/1x,a8,8f15.8)
12    format(/1x,a8,1p,8e15.6)
13    format(/1x,a8,1p,8g15.6)
20    format(/10x,a)
c
      integer ldz, nr, nc, ifmt, nlist
      character*(*) title, labr(*), labc(*), labv
      real*8 z(ldz,nc), v(nc)
c
      integer fmtz, fmtv, jt, jlast, jstrt, j, i
      real*8     zero
      parameter (zero=0d0)
c
      if(ifmt.le.1)then
         assign 1 to fmtz
         assign 11 to fmtv
      elseif(ifmt.eq.2)then
         assign 2 to fmtz
         assign 12 to fmtv
      else
         assign 3 to fmtz
         assign 13 to fmtv
      endif
c
      if ( title .ne. ' ' ) write(nlist,20) title
c
      jlast = 0
      do 400 jstrt = 1, nc, ncol
         jlast = min( nc, jlast+ncol )
c
         write(nlist,10) ( labc(j), j = jstrt, jlast )
         write(nlist,fmtv) labv, ( v(j), j = jstrt, jlast)
         write(nlist,*)
c
         do 300 i = 1, nr
c
c           # print the row if a nonzero element is found.
c
            do 100 j = jstrt, jlast
               if ( z(i,j) .ne. zero ) then
                  write(nlist,fmtz) labr(i), ( z(i,jt), jt=jstrt,jlast)
                  go to 101
               endif
100         continue
101         continue
300      continue
400   continue
c
      return
      end
c deck rddbl
      subroutine rddbl( ndrt, lenbuf, iv, len )
c
c  read an integer vector from the drt file.
c
c  ndrt = drt unit number.
c  lenbuf = maximum logical record length.  vectors longer than
c           lenbuf are logically blocked.
c  iv(*) = output vector.
c  len = vector length parameter.
c      = >0 then read iv(1:len).
c      = <0 then skip over a vector of length (-len).
c           elements iv(1: min( lenbuf, (-len)) ) are referenced.
c
c  06-oct-90 len<0 added. -rls
c  written by ron shepard.
c
       implicit logical(a-z)
      integer ndrt, lenbuf, len
      integer iv(*)
c
      integer ifact, istart, nleft, ibuf, nread, ierr
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      if ( len .ge. 0 ) then
         ifact = 1
      else
         ifact = 0
      endif
      istart = 1
      nleft  = abs(len)
      do 10 ibuf = 1, nleft, lenbuf
         nread = min(nleft,lenbuf)
         call rdivf( ndrt, iv(istart), nread, ierr )
         if ( ierr .ne. 0 )
     &    call bummer( 'rddbl: from rdiv, ierr=', ierr, faterr)
         istart = istart + ifact * nread
         nleft = nleft - nread
10    continue
      return
      end
c deck rdivf
      subroutine rdivf( nunit, iv, len, ierr )
c
c  formatted read of an integer vector from unit nunit.
c
c  nunit = input unit number.
c  iv(*) = integer vector.
c  len   = vector length.  must be greater than zero.
c  ierr  = return status.  0 for normal return.
c
c  01-dec-90 (*) format added. -rls
c
       implicit logical(a-z)
      integer nunit, len, ierr
      integer iv(len)
c
      integer i1
      integer   lenc
      parameter(lenc=20)
      character*(lenc) cfmt
c
c     # initial record is: "cfmt*20             comments..."
c
      read(nunit,fmt='(a)',iostat=ierr) cfmt
      if ( ierr .ne. 0 ) then
          return
      endif
      i1 = 1
      call strskp( cfmt, i1, ' ' )
      if ( i1 .gt. lenc ) then
         ierr = -3
         return
      endif
c
c     # the following records contain iv(1:len) consistent with cmft.
c
      if ( cfmt(i1:) .eq. '(*)' ) then
         read( nunit,    *, iostat=ierr ) iv
      else
         read( nunit, cfmt, iostat=ierr ) iv
      endif
c
      return
      end
c deck rdmoc
c *** this routine will become obsolete when  ***
c *** the new mocoef file format is adopted.  ***
      subroutine rdmoc(
     & nlist, mocoef, nsym, nbpsy,
     & nmpsy, c,      ierr )
c
c  read the rectangular, symmetry-blocked matrix c(*,*).
c
c  input:
c  nlist    = listing file.
c  mocoef   = formatted mocoef file.
c  nsym     = number of symmetry blocks.
c  nbpsy(*) = number of basis functions in each block.
c  nmpsy(*) = number of orbitals in each block.
c             it is assumed that nbpsy(i).ge.nmpsy(i).
c             [these arrays should eventually be read from mocoef. -rls]
c
c  output:
c  c(*) = symmetry-blocked coefficients.
c  ierr =  0 for normal return. this includes some allowed eof.
c       = -1 for unexpected eof error.
c       = -2 for symmetry blocking inconsistency.
c       = -3 for fmt error.
c       =  iostat for the last read statement executed for errors.
c
c  08-oct-90 (columbus day) ifmt1 added, symmetry check added. -rls
c  19-apr-89 (*), (0), and (1) added. -rls
c  13-jun-80 written by ron shepard.
c
       implicit logical(a-z)
      integer nlist, mocoef, nsym, ierr
      integer nbpsy(*), nmpsy(*)
      real*8 c(*)
c
      integer numr, isym, ifmt1, i0, nmo, nbf, imo, i, ibfn
      real*8 coeff
      integer nsbm(8)
      character*40 cfmt
c
      real*8    zero,    one
      parameter(zero=0d0,one=1d0)
c
      ierr = 0
c
      numr = 0
      do 10 isym = 1,nsym
         if ( nmpsy(isym) .gt. nbpsy(isym) ) then
c           # symmetry error.
            ierr = -2
            return
         endif
         nsbm(isym) = numr
         numr = numr + nbpsy(isym) * nmpsy(isym)
10    continue
c
      call wzero( numr, c, 1 )
c
      read(mocoef,'(a)',iostat=ierr) cfmt
      if ( ierr .ne. 0 ) return
c
      write(nlist,6020) cfmt
c
      ifmt1 = 1
      call strskp( cfmt, ifmt1, ' ' )
      if ( ifmt1 .gt. len(cfmt) ) then
         ierr = -3
         return
      endif
c
      if ( cfmt(ifmt1:) .eq. '(*)' ) then
c
c        # list directed read.
c
         i0 = 0
         do 120 isym = 1, nsym
            nmo = nmpsy(isym)
            nbf = nbpsy(isym)
            do 110 imo = 1, nmo
               read(mocoef,*,iostat=ierr) (c(i0+i),i=1,nbf)
               if ( ierr .ne. 0 ) return
               i0 = i0 + nbf
110         continue
120      continue
c
      elseif ( cfmt(ifmt1:) .eq. '(0)'
     &    .or. cfmt(ifmt1:) .eq. '(1)' ) then
c
c        # initialize and read in individual elements.
c        # reduced labels are read, ending with isym.le.0.
c
         if ( cfmt(ifmt1:) .eq. '(1)' ) then
c           # initialize to a unit matrix.
            do 200 isym = 1, nsym
               if ( nmpsy(isym) .ne. 0 ) call wset
     &          (nmpsy(isym), one, c(nsbm(isym)+1), (nbpsy(isym)+1) )
200         continue
         endif
         numr = 0
210      continue
         isym = -1
         ibfn = -1
         imo  = -1
         coeff = zero
         read(mocoef,*,iostat=ierr) isym, ibfn, imo, coeff
         if ( ierr .gt. 0 ) then
c           # error return.
            return
         elseif ( (ierr .lt. 0) .or. (isym .le. 0) .or. (ibfn .le. 0)
     &       .or. (imo .le. 0) ) then
c           # eof normal return.
            ierr = 0
         else
            numr = numr + 1
            c( nsbm(isym) + (imo-1) * nbpsy(isym) + ibfn ) = coeff
            go to 210
         endif
c
      else
c
c        # formatted read.
c
         i0 = 0
         do 320 isym = 1, nsym
            nmo = nmpsy(isym)
            nbf = nbpsy(isym)
            do 310 imo = 1, nmo
               read( mocoef, cfmt, iostat=ierr ) (c(i0+i),i=1,nbf)
               if ( ierr .ne. 0 ) return
               i0 = i0 + nbf
310         continue
320      continue
c
      endif
c
      write(nlist,6030) numr
c
      return
6020  format(/' rdmoc: orbital coefficients will be read using',
     & ' the format:'/1x,a)
6030  format(' rdmoc:',i8,' coefficients were successfully read.')
      end
c deck srtiad
      subroutine srtiad( n, a, indx )
c
c         indexed absolute value sort in decreasing order
c
c  sort the elements a(*) by rearranging the index entries in indx(*)
c  into decreasing order of the absolute values of the vector entries.
c  this routine uses the o(n*log(n)) heapsort method.
c
c  input:
c  n  = vector length of indx(*).
c  a(*) = vector elements.
c  indx(1:n) = initial order of the elements of a(*).
c
c  output:
c  a(*) = unchanged.
c  indx(1:n) = initial entries are rearranged such that
c              abs(a(indx(i))).ge.abs(a(indx(i+1))) for i=1,(n-1).
c
c  written by ron shepard 23-july-87.
c  based on "numerical recipes, the art of scientific computing" by
c  w. h. press, b. p. plannery, s. a. teukolsky, and w. t. vetterling.
c
       implicit logical(a-z)
      integer n, indx(n)
      real*8 a(*)
c
      integer l, ir, indxt, i, j
      real*8 q
c
      l  = n / 2 + 1
      ir = n
10    continue
      if ( l .gt. 1 ) then
         l     = l - 1
         indxt = indx(l)
         q     = abs(a(indxt))
      else
         indxt    = indx(ir)
         q        = abs(a(indxt))
         indx(ir) = indx(1)
         ir       = ir - 1
         if ( ir .le. 1 ) then
            indx(1) = indxt
            return
         endif
      endif
      i = l
      j = l + l
20    if ( j .le. ir ) then
         if ( j .lt. ir ) then
            if ( abs(a(indx(j))) .gt. abs(a(indx(j+1))) ) j = j + 1
         endif
         if ( q .gt. abs(a(indx(j))) ) then
            indx(i) = indx(j)
            i       = j
            j       = j + j
         else
            j = ir + 1
         endif
         go to 20
      endif
      indx(i) = indxt
      go to 10
c
      end
c deck srtii
      subroutine srtii( n, a, indx )
c
c         indexed sort in increasing order
c
c  sort the elements a(*) by rearranging the index entries in indx(*)
c  into increasing order of the vector entries.
c  this routine uses the o(n*log(n)) heapsort method.
c
c  input:
c  n  = vector length of indx(*).
c  a(*) = vector elements.
c  indx(1:n) = initial order of the elements of a(*).
c
c  output:
c  a(*) = unchanged.
c  indx(1:n) = initial entries are rearranged such that
c              a(indx(i)).le.a(indx(i+1)) for i=1,(n-1).
c
c  written by ron shepard 23-july-87.
c  based on "numerical recipes, the art of scientific computing" by
c  w. h. press, b. p. plannery, s. a. teukolsky, and w. t. vetterling.
c
       implicit logical(a-z)
      integer n, indx(n)
      real*8 a(*)
c
      integer l, ir, indxt, i, j
      real*8 q
c
      l  = n / 2 + 1
      ir = n
10    continue
         if ( l .gt. 1 ) then
            l     = l - 1
            indxt = indx(l)
            q     = a(indxt)
         else
            indxt    = indx(ir)
            q        = a(indxt)
            indx(ir) = indx(1)
            ir       = ir - 1
            if ( ir .le. 1 ) then
               indx(1) = indxt
               return
            endif
         endif
         i = l
         j = l + l
20       if ( j .le. ir ) then
            if ( j .lt. ir ) then
               if ( a(indx(j)) .lt. a(indx(j+1)) ) j = j + 1
            endif
            if ( q .lt. a(indx(j)) ) then
               indx(i) = indx(j)
               i       = j
               j       = j + j
            else
               j = ir + 1
            endif
            go to 20
         endif
         indx(i) = indxt
      go to 10
c
      end
c deck iargc
*@ifdef fujitsu
*      integer function iargc()
*c
*c iargc.f - equiv of iargc() in Sun Fortran
*c roger edberg  14-nov-91
*c
*      integer*4     j
*      character*1   cargj
*c
*      j = 1
*1000  continue
*         call getarg( j, cargj )
*c        # cargj is returned left-justified, so ' ' implies
*c        # that j has overrun the command-line argument list.
*         if ( cargj .eq. ' ' ) goto 1100
*         j = j + 1
*      goto 1000
*1100  continue
*c
*      iargc = j - 1
*c
*      return
*      end
*@endif
c deck hostnm
c  Is this routine now redundant with $COLUMBUS/special//hostnm.c ?
c  If so, then we should consider removing it here. 28-apr-92 -rls
*@ifdef iris
*      integer function hostnm( name )
*c
*c  get name of current host ( IRIS )
*c
*c  11-sep-91 written by matthias schueler.
*c
*      character*(*) name
*c
*      intrinsic len
*c
*      integer  gethostname
*      external gethostname
*c
*      hostnm = gethostname( name, len(name) )
*c
*      return
*      end
*@endif
c deck skpx01
      subroutine skpx01( lenv, vector, icode, numv1 )
c
c  transform between the skip-vector form and the 0/1 form
c  of an index vector.  this is done in-place.
c  0/1 form:    1 0 0 1 1 0 0 0 1 0 1 ...
c  skip-vector: 0 2 1 0 0 3 2 1 0 1 0 ...
c
c  input:
c  lenv = vector length.
c  vector(*) = input vector.
c  icode = 0  convert from 0/1 form to skip-vector form.
c        = 1  convert from skip-vector form to 0/1 form.
c
c  output:
c  vector(*) = transformed vector.
c  numv1 = number of 1 entries in the 0/1 form.
c
c  05-jun-89 written by ron shepard.
c
       implicit logical(a-z)
       integer lenv, icode, numv1
      integer vector(lenv)
c
      integer i
      integer nx(0:1)
c
      integer    toskp,   to01
      parameter( toskp=0, to01=1 )
c
      if ( icode .eq. toskp ) then
c        # 0/1 to skip-vector form.
         numv1 = 0
         nx(0) = 0
         nx(1) = -1
         do 10 i = lenv, 1, -1
            numv1 = numv1 + vector(i)
            nx(0) = nx(vector(i)) + 1
            vector(i) = nx(0)
10       continue
      else
c        # skip-vector to 0/1 form.
         numv1 = 0
         do 20 i = 1, lenv
            if ( vector(i) .eq. 0 ) then
               vector(i) = 1
               numv1 = numv1 + 1
            else
               vector(i) = 0
            endif
20       continue
      endif
c
      return
      end
c deck indx01
      subroutine indx01( nlist, numi, len01, vec01, indxv, icode )
c
c  transform between the index representation of a vector and the
c  0/1 representation.
c  indxv(*) = 1,2,  4,    7
c  vec01(*) = 1,1,0,1,0,0,1
c
c  arguments:
c  numi = number of index-representation entries.
c  len01 = length of the 0/1 representation vector.
c  vec01(1:len01) = 0/1 representation.
c  indxv(1:numi) = index vector representation.
c  icode = 0  convert from vec01(*) to indxv(*),
c        = 1  convert from indxv(*) to vec01(*).
c
c  05-jun-89 written by ron shepard.
c
       implicit logical(a-z)
      integer nlist, len01, numi, icode
      integer vec01(len01),indxv(numi)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer    toindx,   to01
      parameter( toindx=0, to01=1 )
c
      integer i, icnt
c
      if ( icode .eq. toindx ) then
c        # vec01(*) to indxv(*).
         icnt = 0
         do 10 i = 1, len01
            if ( vec01(i) .eq. 1 ) then
               icnt = icnt + 1
               indxv(icnt) = i
            endif
10       continue
         if ( icnt .ne. numi )
     &    call bummer('indx01: (numi-icnt)=',(numi-icnt),faterr)
         write(nlist,'(/1x,a,i6,a)')'indx01:',
     &    numi,' indices saved in indxv(*)'
      else
c        # indxv(*) to vec01(*).
         do 20 i = 1, len01
            vec01(i) = 0
20       continue
         do 30 icnt = 1, numi
            vec01(indxv(icnt)) = 1
30       continue
         write(nlist,'(1x,a,i6,a)')'indx01:',
     &    numi,' elements set in vec01(*)'
      endif
      return
      end
c deck symcvt
      subroutine symcvt(
     & nsym,   ldso,   symorb, nmpsy,
     & nmskp,  norb,   ldout,  outorb,
     & ierr )
c
c  convert from symmetry-reduced orbitals to global orbitals.
c
c  input:
c  nsym = number of symmetry blocks.
c  ldso = leading dimension of the symorb(*,*) array.
c  symorb(1:ldso,1:norb) = symmetry orbitals to be converted.
c                 symorb(1,i) is the orbital symmetry; symorb(1,i)<=0
c                 is the end of data marker.  symorb(2,i) is the
c                 reduced orbital index within the symmetry block.
c                 if symorb(2,i)<0, then the orbital is counted
c                 from the end of the symmetry block instead of
c                 the beginning.
c  nmpsy(1:nsym) = number of orbitals in each symmetry block.
c  nmskp(1:nsym) = global-orbital offsets.
c  ldout = leading dimension of the outorb(*,*) array.  This is the
c          increment between consecutive elements.
c
c  output:
c  norb = number of orbitals that were converted.
c  outorb(1:norb) = global orbital indices.  errors are indicated by
c                   zero entries.
c  ierr = total number of errors generated; 0 for normal return.
c
c  14-dec-90 nsym, ierr ldso, ldout, added to the argument list. -rls
c  21-sep-87 written by ron shepard.
c
       implicit logical(a-z)
      integer nsym, ldso, norb, ldout, ierr
      integer symorb(ldso,*), nmpsy(nsym), nmskp(nsym), outorb(ldout,*)
c
      integer i, isym, imo, next
c
      ierr = 0
      i    = 0
      next = 1
10    if ( symorb(1,next) .gt. 0 ) then
         i    = next
         next = next + 1
         isym = symorb(1,i)
         imo  = symorb(2,i)
         if ( isym .le. nsym ) then
            if ( (imo .gt. 0) .and. (imo .le. nmpsy(isym)) ) then
               outorb(1,i) = nmskp(isym) + imo
            elseif ( (imo .lt. 0) .and. (-imo .le. nmpsy(isym)) ) then
               outorb(1,i) = nmskp(isym) + imo + nmpsy(isym) + 1
            else
               ierr = ierr + 1
               outorb(1,i) = 0
            endif
         else
            ierr = ierr + 1
            outorb(1,i) = 0
         endif
         goto 10
      endif
c
      norb = i
c
      return
      end
c deck flushx
      subroutine flushx( iunit )
c
c  flush any file buffers associated with fortran unit "iunit".
c
c  this routine is used primarily for flushing listing files
c  during iterative procedures.
c
c  06-may-92 written by ron shepard.
c
       implicit logical(a-z)
      integer iunit
c
*@if defined ( cray) || defined (T3E64) || defined ( T3D )
*c     # must include ierr to avoid aborts.
*      integer ierr
*      call flush( iunit, ierr )
*@elif defined( rs6000)  &&  (! defined( extname))
*c     # must access the fortran library routine, with the underscore,
*c     # not the C library function, without the underscore.
*      call flush_( iunit )
*@elif defined (  sun) || defined ( sgipower) || defined ( rs6000) || defined ( macabsoft)
*      call flush( iunit )
*@else
c     # no-op call.
      continue
*@endif
c
      return
      end
c deck skpend
c
c  usage:
c      open(iunit,file=tmpfile,form='formatted')
c      iform = 1   ! iform = 0 for 'unformatted'
c      call iskpend( iunit, iform, ierr )   !this overwrites tmpfile
c      if ( ierr .ne. 0 ) stop 'error in iskpend()'
c      ...
c      open(iunit2,file=filename,form='formatted')
c      iform = 1   ! iform = 0 for 'unformatted'
c      call skpend( iunit2, iform )
c      write(iunit2, ...) appended_record
c
c  only one iskpend() call is required for an arbitrary number of
c  skpend() calls and for an arbitrary number of units.
c
      subroutine skpend( iunit, iform )
c
c  skip to the end of file on unit=iunit.
c
c  input:
c  iunit = unit number.  may be either formatted or unformatted.
c          (this routine assumes the same behavior for both types
c          of files.)
c  iform =  0 for unformatted i/o.
c        <> 0 for formatted i/o.
c
c  upon return, iunit is positioned so that a subsequent
c  write(iunit) will succeed and will write after the last
c  record.
c
c  the fortran standard defines how this should be done, but
c  unfortunately some vendors do not follow this convention.
c  this routine tries to account for this nonstandard behavior.
c
c  version log:
c  25-oct-01 written by ron shepard.
c
      implicit logical(a-z)
       integer iunit, iform, ierr
c
c     # local:
      integer i
c
      character*(*) cfmt
      parameter( cfmt='(i1)' )
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # valid values for nback are either 0 or 1.
      integer nback
      save    nback
      data    nback/-1/
c
      if ( nback .lt. 0 ) then
c        # isetend() was not called or was unsuccessful.
         call bummer('setback: isetend() was uncessful',0,faterr)
      endif

      rewind iunit
c
      if ( iform .ne. 0 ) then
c        # formatted file
100      continue
         read(iunit,*,end=102)
         goto 100
      else
c        # unformatted file
101      continue
         read(iunit,end=102)
         goto 101
      endif
102   continue
c     # the file should now be positioned after the endfile record.
c     # backspace the appropriate number of times for this
c     # machine.
      do i = 1, nback
         backspace iunit
      enddo
c
      return
c
      entry iskpend( iunit, iform, ierr )
c
c     # initialize local variable nback in the skpend() routine.
c
c     # determine the the number of backspaces required for this
c     # machine.  Only one iskpend() call is required to
c     # initialize the variable nback, then an arbitrary number
c     # of skpend() calls can be made on an arbitrary number of
c     # files and unit numbers.
c
c     # input:
c     # iunit = unit number.
c     # iform =  0 for unformatted i/o.
c     #       <> 0 for formatted i/o.
c
c     # output:
c     # ierr =  0 for success
c     #      <> 0 for failure.
c     # the contents of the file connected to iunit are overwritten and
c     # iunit is left positioned at the beginning of the file.
c
      ierr = 0
      rewind iunit
      if ( iform .ne. 0 ) then
c
c        # formatted file
c
c        # write two formatted records.
         write(iunit,cfmt) 1
         write(iunit,cfmt) 2
         rewind iunit
200      continue
         read(iunit,*,end=201)
         goto 200
201      continue
c        # the file is positioned at the end, after the endfile
c        # mark, if it exists.
c        # backspace twice, then see which record is read next.
         backspace iunit
         backspace iunit
         read(iunit,cfmt) i
         if ( i .eq. 1 ) then
c           # this machine is not standard conforming
            nback = 0
         elseif ( i .eq. 2 ) then
c           # this machine is standard conforming.
            nback = 1
         else
c           # could not determine the correct value.
            nback = -1
            ierr  = 1
         endif
      else
c
c        # unformatted file
c
c        # write two unformatted records.
         write(iunit) 1
         write(iunit) 2
         rewind iunit
300      continue
         read(iunit,end=301)
         goto 300
301      continue
c        # the file is positioned at the end, after the endfile
c        # mark if it exists.
c        # backspace twice, then see which record is read next.
         backspace iunit
         backspace iunit
         read(iunit) i
         if ( i .eq. 1 ) then
c           # this machine is not standard conforming
            nback = 0
         elseif ( i .eq. 2 ) then
c           # this machine is standard conforming.
            nback = 1
         else
c           # could not determine the correct value.
            nback = -1
            ierr  = 1
         endif
      endif
c
c     # leave the file empty and positioned at the beginning so that
c     # write(iunit,...) will be a valid operation.
c
      rewind  iunit
      endfile iunit
      rewind  iunit
      return
      end
c deck mulpop
      subroutine mulpop(
     & an,     c,      scr1,    ityp,
     & mnl,    ms,     mtype,   ni,
     & nd_at,  scr2,   s,       nlist,
     & nbft,   nmpsy,  nsym)
c
c  10-dec-01 modified and included into colib by m.dallos
c  simple mulliken population analysis
c
c  ##  valiable meaning:
c     an     - occupation number
c     c      - NOs
c     scr1   - scratch array, dim=9*nd_at
c     ityp   - symmetry label
c     mnl    - pointer to the type of ao in an so (e.g. 1s,2p,3d)
c     ms     - pointer to the label for the atom(s) the so is on
c     mtype  - AO orbital label
c     ni     - number of internal orbitals
c     nd_at  - number of symmtery distinct itypes of atoms
c     scr2   - scratch array, dim=9*nd_at*nbpsy_mx
c     s      - sao - overlap intergrals
c     nlist  - listing file number
c     nbft   - total number of basis functions
c     nmpsy  - number of SAO per irrep
c     nsym   - number of irreducible representations
c
      implicit logical(a-z)
c  ##  parameter section
c
      real*8 a0
      integer lp1u
      parameter (a0=0.0d0, lp1u=9)
c
c  ##  integer section
c
      integer iim1,ip, isym, iis, isl, ii, ir, is, iq, isq, ipu,
     &        iu, isp, i, il, isr
      integer lapqis, lapi, lai, lap, lapq, lp1, lp1p, laqi, lapis,
     &        lapips, laq, lp1q
      integer ms(*), mnl(*)
      integer nd_at, ni(*), nmpsy(*), nbft, nlist, nsym
c
c  ## real*8
c
      real*8 an(*)
      real*8 c(*)
      real*8 numel
      real*8 prd
      real*8 scr1(nd_at,lp1u), scr2(nd_at,lp1u,*), s(*)
      real*8 tot(6)
c
c  ##  character section
c
      character*4 ityp(*)
      character*3 mtype(*)
c
c  ##  data section
c
      character*2 lbll(21)
      integer lfrnl(121)
      data lbll /' s',' p',' d',' f',' g',' h',' i',' k',' l',' m',
     &      ' n',' o',' q',' r',' t',' u',' v',' w',' x',' y',' z'/,
     &  lfrnl /1,2,1,3,2,4,1,3,5,2,4,6,1,3,5,7,2,4,6,8,1,3,5,7,9,
     & 2,4,6,8,10,1,3,5,7,9,11,2,4,6,8,10,12,1,3,5,7,9,11,13,
     & 2,4,6,8,10,12,14,1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16,
     & 1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18,
     & 1,3,5,7,9,11,13,15,17,19,2,4,6,8,10,12,14,16,18,20,
     & 1,3,5,7,9,11,13,15,17,19,21/
c
      iim1(i)=(i*(i-1))/2
c
c----------------------------------------------------------------------
c
      write (nlist,"(//10x,'Mulliken population analysis'//)")
      write(nlist,*)
     & ' NOTE: For HERMIT use spherical harmonics basis sets !!!'
      write(nlist,*)
      call wzero((lp1u*nd_at),scr1,1)
      lapq=0
      lap=0
      lai=0
      numel=0.0d+00
      lapi=0
      do 250 isym=1,nsym
      ipu=nmpsy(isym)
      if(ipu.eq.0) go to 250
      iu=ni(isym)
      if(iu.eq.0) then
        lap=lap+ipu
        lapq=lapq+iim1(ipu+1)
        go to 240
      endif
      call wzero((lp1u*nd_at)*iu,scr2,1)
      lapqis=lapq
      lapis=lap
      do 160 i=1,iu
        lapq=lapqis
        lap=lapis
        lai=lai+1
        lapips=lapi
        do 120 ip=1,ipu
          lap=lap+1
          isp=ms(lap)
          lp1p=lfrnl(mnl(lap))
          lapi=lapi+1
          laq=lapis
          laqi=lapips
          do 110 iq=1,ip
            lapq=lapq+1
            laq=laq+1
            isq=ms(laq)
            lp1q=lfrnl(mnl(laq))
            laqi=laqi+1
            prd=c(lapi)*s(lapq)*c(laqi)
            scr2(isp,lp1p,i)=scr2(isp,lp1p,i)+prd
            if(ip.ne.iq) scr2(isq,lp1q,i)=scr2(isq,lp1q,i)+prd
  110     continue
  120   continue
        do 150 is=1,nd_at
          do 140 lp1=1,lp1u
            scr2(is,lp1,i)=an(lai)*scr2(is,lp1,i)
            scr1(is,lp1)=scr1(is,lp1)+scr2(is,lp1,i)
  140     continue
  150   continue
  160 continue
      write (nlist,910) ityp(isym)(2:4)
      do 230 il=1,iu,6
        ir=min(il+5,iu)
        if(il.ne.1) write (nlist,*)
        write (nlist,920) (i, ityp(isym)(2:4), i=il,ir)
        do 220 is=1,nd_at
          do 210 lp1=1,lp1u
            do 200 i=il,ir
              if(scr2(is,lp1,i).ne.a0) then
                write (nlist,930) mtype(is),lbll(lp1),
     &            (scr2(is,lp1,ii),ii=il,ir)
                go to 210
              endif
  200       continue
  210     continue
  220   continue
  230 continue
  240 lai=lai+(ipu-iu)
      lapi=lapi+ipu*(ipu-iu)
  250 continue
      write (nlist,940)
      do 330 isl=1,nd_at,6
        isr = min(isl+5,nd_at)
        if(isl.ne.1) write (nlist,*)
        call wzero((isr-isl+1),tot,1)
        write (nlist,950) (mtype(is), is=isl,isr)
        do 320 lp1=1,lp1u
          do 310 is=isl,isr
            if(scr1(is,lp1).ne.a0) then
              write (nlist,960) lbll(lp1), (scr1(iis,lp1), iis=isl,isr)
              do 300 iis=isl,isr
                tot(iis-isl+1) = tot(iis-isl+1) + scr1(iis,lp1)
                numel = numel + scr1(iis,lp1)
  300         continue
              go to 320
            endif
  310     continue
  320   continue
        write(nlist,970) (tot(is-isl+1), is=isl,isr)
        write(nlist,*)
  330 continue
      write(nlist,'(/,'' Total number of electrons: '',f13.8/)')numel
c
      return
  910 format(/24x,a3,33h partial gross atomic populations)
  920 format(3x,8hao class,6(i8,a3))
  930 format(4x,a3,a2,4x,6f11.6)
  940 format(//24x,24hgross atomic populations)
  950 format(5x,2hao,3x,6(8x,a3))
  960 format(5x,a2,6x,6f11.6)
  970 format(4x,'total',4x,6f11.6)
      end
c
c
c****************************************************************
c
      subroutine mulpop_map(
     & bfnl,   imtype,   map,     mnl,
     & ms,     nmap,     nd_at,   mtype,
     & nbft,   nsym,     nlist)
c
c  extract mnl(*) and ms(*) from map(*,*).
c  ##  OUTPUT:
c       # nd_at = number of symmetry-distinct centers
c       # mtype(i) = character center-label for each distinct center
c       # nd_at = number of symmtery distinct itypes of atoms
c
c  10-dec-01 modified and included into colib by m.dallos
c  25-may-91 written by russ pitzer.
c
      implicit logical(a-z)
c  ## parameter & common block section
c
      integer mnltyp,mstyp
      parameter (mnltyp=4,mstyp=3)
c
c  ##  integer section
c
      integer i,imnl,ims, imtype(*), isym, idebug
      integer nmap, nd_at, nsym, nbft, nlist
      integer map(nbft,*),mnl(nbft),ms(nbft)
c
c  ##  character
c
      character*8 bfnl(*)
      character*3 mtype(*)
c
c------------------------------------------------------------------
c
      idebug=0 ! for debug set to 1
c
      imnl=0
      ims=0
      do i=1,nmap
        if (imtype(i).eq.mnltyp) then
          imnl=i
        else if (imtype(i).eq.mstyp) then
          ims=i
        endif
      enddo
c
c     # extract the bfntyp(*) array from map(*,*).
      if (imnl.eq.0) then
c       # default: set everything to an s-type function.
        call iset (nbft,1,mnl,1)
      else
        call icopy (nbft,map(1,imnl),1,mnl,1)
      endif
c
c     # extract the bfn-to-center array from map(*,*).
      if (ims.eq.0) then
c       # default: set everything to one center.
        call iset (nbft,1,ms,1)
      else
        call icopy (nbft,map(1,ims),1,ms,1)
      endif
c
c       # set nd_at = number of symmetry-distinct centers
c       # set mtype(i) = character center-label for each distinct center
        nd_at = 0
        do isym=1,nbft
          nd_at = max( nd_at, ms(isym) )
            mtype( ms(isym) ) = bfnl(isym)(4:5)//'_'
        enddo
c
c     # write out mtype and the map vectors ms(*) and mnl(*):
      if (idebug.eq.1) then
c
       write(nlist,*)'mtype:',(mtype(isym),isym=1,nd_at)
c
       write(nlist,"('input orbital labels, i:bfnlab(i)=')")
       write(nlist,"(6(i4,':',a8))") (i,bfnl(i),i=1,nbft)
c
       write(nlist,"('bfn_to_center map(*), i:map(i)')")
       write(nlist,"(10(i4,':',i3))") (i,ms(i),i=1,nbft)
c
       write(nlist,"('bfn_to_orbital_type map(*), i:map(i)')")
       write(nlist,"(10(i4,':',i3))") (i,mnl(i),i=1,nbft)
       write(nlist,*)' Number of symmetry-distinct centers:',nd_at
c
      endif ! end of: if (idebug.eq.1) then
c
      return
      end
c


