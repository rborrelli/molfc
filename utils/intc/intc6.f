      Subroutine pentahex(i,nlig1,nlig2,ilig1,ilig2,ipun,iptyp,
     1    ncg,isymb,scale)
cversion 5.0
      implicit real*8 (a-h,o-z)
c     build coordinates around a penta-  or hexavalent atom i
c     nlig1,2- number of ligands, ilig1,2(*) - list of ligands
c     there are two types (1,2), e.g. terminal and non-terminal
c     neigbors of i (this is not always significant)
      character*2 isymb
      parameter (mxneig=6,mxlin=4)
      parameter (ztol=1.d-4,collin=2.8)
      dimension ilig1(nlig1), ilig2(nlig2)
      dimension scale(*),isymb(*)
      dimension ilig(mxneig),it(mxneig),iz(mxneig),ieq(mxneig)
      dimension xz(mxneig),dist(mxneig),it1(mxneig)
      dimension klin(mxlin),jlin(mxlin)
c
c     in the following rearrangment of atoms, ilig will be the
c     list considered; in case of terminal atoms only, ilig=ilig1,
c     obviously;
c     for the moment, terminal and non-terminal ligands
c     will not be differentiated, and ilig1 and ilig2
c      will be concatenated in ilig; may be refined later
c
      nlig=nlig1+nlig2
      do 15 k=1,nlig1
        ilig(k)=ilig1(k)
 15   continue
      do 25 k=1,nlig2
        kk=nlig1+k
        ilig(kk)=ilig2(k)
 25   continue
c
      call cartes(.false.,i,xi,yi,zi,izi)
c  never used:      zaver=0.0d0
      daver=0.0d0
      do 10 k=1,nlig
        call cartes(.false.,ilig(k),xx,yy,zz,iz(k))
        xz(k)=iz(k)
c  never used:zaver=zaver+xz(k)
        dist(k)=sqrt((xx-xi)**2+(yy-yi)**2+(zz-zi)**2)
        daver=daver+dist(k)
10    continue
c  never used    zaver=zaver/nlig
      daver=daver/nlig
c     rearrange the ligands  acc. to their atomic numbers,
c     and acc. to their distance from center, resp., in descending
c     order
      do 20 k=1,nlig
        it(k)=ilig(k)
        it1(k)=ilig(k)
20    continue
      call sequence(nlig,xz,it)
      call sequence(nlig,dist,it1)
c
c     array 'it' is the rearranged list based on the nucl.charges,
c      it1 another rearranged list based on the distances from
c      the central atom
c      xz and dist are rearranged accordingly
c
c      determine the linear units around the central atom:
      alphmax=0.
      nlin=0
      do 30 k=2,nlig
      do 30 j=1,k
        call angle(it(k),i,it(j),alph)
        if(alph.gt.collin) then
          nlin=nlin+1
          klin(nlin)=it(k)
          jlin(nlin)=it(j)
c        these are the indices of the terminal atoms of a linear unit
          if(alph.gt.alphmax) then
            alphmax=alph
c  what was this?, it's  never used            nlmax=nlin
c         index defining the unit of maximum linearity
          endif
        endif
30    continue
c
c...   pentavalent case:
      if(nlig.eq.5) then
c     before the general case, somewhat arbitrary treatment of 1
c     and 2 hydrogen ligands, resp.:
c
        if((xz(5)-1.).lt.ztol.and.(xz(4)-1.).gt.1.) then
c       one hydrogen, treated as AX4H, independent of geometry
          my=it(5)
          call ax4y (ipun,iptyp,ncg,i,it(1),it(2),it(3),it(4),my,
     1    isymb,scale)
        else if((xz(5)-1.).lt.ztol.and.(xz(4)-1.).lt.ztol.and.
     1     (xz(3)-1.).gt.ztol) then
c       two hydrogens, treated as AX3H2, hydrogens axial, indep-
c       endent of geometry
          nax1=it(4)
          nax2=it(5)
          call ax3y2(ipun,iptyp,ncg,i,it(1),it(2),it(3),nax1,nax2,
     1    isymb,scale)
        else
c       in the general case,geometry checked in a simplified way:
c       decide between tetragonal pyramid and trigon. bipyr.,
c       on the basis of linear units:
          if(nlin.eq.0) then
c TEMP
c           print *, 'no linear unit in pentavalent'
c           tetragon.pyramid will be assumed, the unique atom selected
c           based on distances,either the largest or the smallest:
            if((it1(1)-daver).gt.(daver-it1(5))) then
             my=it1(1)
             call ax4y (ipun,iptyp,ncg,i,it1(2),it1(3),it1(4),it1(5),my,
     1       isymb,scale)
            else
             my=it1(5)
             call ax4y (ipun,iptyp,ncg,i,it1(1),it1(2),it1(3),it1(4),my,
     1      isymb,scale)
            endif
          else if(nlin.eq.1) then
c         axial atoms are those that formed the linear unit:
            nax1=klin(1)
            nax2=jlin(1)
c          now find the rest,these will be the equatorial atoms:
            neq=0
            do 40 k=1,nlig
             if(it(k).eq.nax1.or.it(k).eq.nax2) then
               goto 40
             else
               neq=neq+1
               ieq(neq)=it(k)
             endif
40          continue
            if(neq.ne.3) then
              print *, 'logical program error in pentahex'
              call bummer('logical error in pentahex',0,2)
              STOP
            else
              neq1=ieq(1)
              neq2=ieq(2)
              neq3=ieq(3)
              call ax3y2(ipun,iptyp,ncg,i,neq1,neq2,neq3,nax1,nax2,
     1        isymb,scale)
            endif
          else if(nlin.eq.2) then
c   find the fifth atom,which is not part of a linear unit,
c   this will be the top of the square pyramid:
            do 50 k=1,nlig
            do 52 kk=1,2
             if(it(k).eq.klin(kk).or.it(k).eq.jlin(kk)) then
               goto 50
             endif
 52         continue
             my=it(k)
             goto 55
  50        continue
  55        continue
            call ax4y(ipun,iptyp,ncg,i,klin(1),klin(2),jlin(1),jlin(2),
     1      my,  isymb,scale)
          else if(nlin.gt.2) then
            print *, 'more than 3 linear units in pentavalent???'
            print *, 'nlin=',nlin
          endif
        endif
c
      elseif(nlig.eq.6) then
c... hexavalent case:
c    TEMP
c       print *, 'hexavalent'
c    all cases considered simply as AX4Y2 tetrag. bipyramid;
c TEMP
c       print *, 'nlin and array it1', nlin, it1
        if(nlin.eq.0) then
c   in case of no linear units, selection of axial and equat.
c   atoms is arbitrary
c   the two longest distances will be taken axial
          my1=it1(1)
          my2=it1(2)
          call ax4y2(ipun,iptyp,ncg,i,it1(3),it1(4),it1(5),it1(6),
     1    my1,my2,isymb,scale)
        else
c    take the atom at longest distance, and its collinear
c    partner , as the axial unit:
          my1=it1(1)
          ifound=0
          do 60 k=1,nlin
           if(klin(k).eq.my1) then
             my2=jlin(k)
             ifound=1
             goto 65
           elseif(jlin(k).eq.my1) then
             my2=klin(k)
             ifound=1
             goto 65
           endif
  60      continue
  65      continue
c TEMP
c          print *, 'ifound', ifound
          if(ifound.eq.0) then
c     the longest bond has no collinear partner
c     take simply the existing linear unit as axial:
            my1=klin(1)
            my2=jlin(1)
          endif
c        now find the rest,these will be the equatorial atoms:
          neq=0
          do 45 k=1,nlig
           if(it1(k).eq.my1.or.it1(k).eq.my2) then
             goto 45
           else
             neq=neq+1
             ieq(neq)=it1(k)
           endif
45         continue
c TEMP
c          print *, 'neq=', neq
c          print *, 'ieq', ieq
c          print *, 'my1,my2', my1, my2
          call ax4y2(ipun,iptyp,ncg,i,ieq(1),ieq(2),ieq(3),ieq(4),
     1      my1,my2,isymb,scale)
        endif
      else
        print *, '??? number of ligands must be 5 or 6 in pentahex'
      endif
      return
      end
c end of pentahex
c
c
      Subroutine sequence (n,a,ib)
      implicit real*8 (a-h,o-z)
      dimension a(n),ib(n)
c......
      do 10 i=1,n
        amax=0.
        kmax=1
        ipartner=ib(1)
      do 20 k=i,n
        if(a(k).gt.amax) then
          kmax=k
          amax=a(k)
          ipartner=ib(k)
        endif
20    continue
      temp=a(i)
      itempb=ib(i)
      a(i)=amax
      ib(i)=ipartner
      a(kmax)=temp
      ib(kmax)=itempb
10    continue
      return
      end

      Subroutine ax3y2(ipun,iptyp,ncg,i,mx1,mx2,mx3,my1,my2,
     1  isymb,scale)
      implicit real*8 (a-h,o-z)
c
c     this generates the 7 deformational coordinates around
c     a pentavalent atom i, for the trigonal bipyramidal case
      character*2 isymb
      dimension scale(*),isymb(*)
      ncg=ncg+1
      write(ipun,100) ncg,1,mx1,mx2,mx3,i
      write(ipun,110) 1,mx2,mx3,mx1,i
      write(ipun,110) 1,mx3,mx1,mx2,i
        write(iptyp,500) ncg,isymb(i),i,scale(ncg)
 500     format(i4,'. ','trigbipyrmdA=',a2,i2,
     1     ' OUT   ','sc=',f5.4,2x)
c....
      ncg=ncg+1
      write(ipun,200)  ncg,2,mx2,mx3,i
      write(ipun,210) -1,mx1,mx2,i
      write(ipun,210) -1,mx1,mx3,i
        write(iptyp,510) ncg,isymb(i),i,scale(ncg)
 510     format(i4,'. ','trigbipyrmdA=',a2,i2,
     1     ' ADEFa ','sc=',f5.4,2x)
c..
      ncg=ncg+1
      write(ipun,200)  ncg,1,mx1,mx2,i
      write(ipun,210) -1,mx1,mx3,i
        write(iptyp,520) ncg,isymb(i),i,scale(ncg)
 520     format(i4,'. ','trigbipyrmdA=',a2,i2,
     1     ' ADEFb ','sc=',f5.4,2x)
c...
      ncg=ncg+1
      write(ipun,200)  ncg,2,mx1,my1,i
      write(ipun,210) -1,mx2,my1,i
      write(ipun,210) -1,mx3,my1,i
      write(ipun,210)  2,mx1,my2,i
      write(ipun,210) -1,mx2,my2,i
      write(ipun,210) -1,mx3,my2,i
        write(iptyp,530) ncg,isymb(i),i,scale(ncg)
 530     format(i4,'. ','trigbipyrmdA=',a2,i2,
     1     7H ROCKa','sc=',f5.4,2x)
c...
      ncg=ncg+1
      write(ipun,200)  ncg,2,mx1,my1,i
      write(ipun,210) -1,mx2,my1,i
      write(ipun,210) -1,mx3,my1,i
      write(ipun,210) -2,mx1,my2,i
      write(ipun,210) +1,mx2,my2,i
      write(ipun,210) +1,mx3,my2,i
        write(iptyp,540) ncg,isymb(i),i,scale(ncg)
 540     format(i4,'. ','trigbipyrmdA=',a2,i2,
     1     ' ROCKa"','sc=',f5.4,2x)
c..
      ncg=ncg+1
      write(ipun,200)  ncg,1,mx2,my1,i
      write(ipun,210) -1,mx3,my1,i
      write(ipun,210)  1,mx2,my2,i
      write(ipun,210) -1,mx3,my2,i
        write(iptyp,550) ncg,isymb(i),i,scale(ncg)
 550     format(i4,'. ','trigbipyrmdA=',a2,i2,
     1     7H ROCKb','sc=',f5.4,2x)
c..
      ncg=ncg+1
      write(ipun,200)  ncg,1,mx2,my1,i
      write(ipun,210) -1,mx3,my1,i
      write(ipun,210) -1,mx2,my2,i
      write(ipun,210) +1,mx3,my2,i
        write(iptyp,560) ncg,isymb(i),i,scale(ncg)
 560     format(i4,'. ','trigbipyrmdA=',a2,i2,
     1     ' ROCKb"','sc=',f5.4,2x)
c.....
100   format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
110   format(' ',9x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
200   format(i4,2x,'K',3x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
210   format(' ',9x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
      return
      end

      Subroutine ax4y(ipun,iptyp,ncg,i,mx1,mx2,mx3,mx4,my,
     1 isymb,scale)
      implicit real*8 (a-h,o-z)
c     this generates the coordinates around a pentavalent atom i,
c     for the tetragonal pyramidal case
      character*2 isymb
      parameter (ztol=1.d-4)
      dimension scale(*),isymb(*)
c
c     first: rearrange the ligands in cyclic order
c     by finding the atom opposite to mx1 based on the largest angle:
      call angle(mx1,i,mx2,bet12)
      call angle(mx1,i,mx3,bet13)
      call angle(mx1,i,mx4,bet14)
      bmax=bet13
      icase=13
c in the default case 3. is opposite to 1.
c and indexes will not be changed; else:
      if((bet12-bmax).gt.ztol) then
         bmax=bet12
         icase=12
      elseif ((bet14-bmax).gt.ztol)then
        bmax=bet14
        icase=14
      endif
      if(icase.eq.12) then
        mxtemp=mx3
        mx3=mx2
        mx2=mxtemp
      elseif (icase.eq.14) then
        mxtemp=mx3
        mx3=mx4
        mx4=mxtemp
      endif
c
c  out1, (symm)
      ncg=ncg+1
      write(ipun,100) ncg,1,mx1,mx2,mx3,i
      write(ipun,110) 1,mx2,mx3,mx4,i
      write(ipun,110) 1,mx3,mx4,mx1,i
      write(ipun,110) 1,mx4,mx1,mx2,i
        write(iptyp,500) ncg,isymb(i),i,scale(ncg)
 500    format(i4,'. ','tetragpyrmdA=',a2,i2,
     1       ' OUTa1 ','sc=',f5.4,2x)
c  out2
      ncg=ncg+1
      write(ipun,100) ncg,1,mx1,mx2,mx3,i
      write(ipun,110) -1,mx2,mx3,mx4,i
      write(ipun,110) 1,mx3,mx4,mx1,i
      write(ipun,110) -1,mx4,mx1,mx2,i
        write(iptyp,510) ncg,isymb(i),i,scale(ncg)
 510    format(i4,'. ','tetragpyrmdA=',a2,i2,
     1       ' OUTb1 ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1, mx1,mx2,i
      write(ipun,210) -1, mx2,mx3,i
      write(ipun,210) +1, mx3,mx4,i
      write(ipun,210) -1, mx4,mx1,i
        write(iptyp,520) ncg,isymb(i),i,scale(ncg)
 520    format(i4,'. ','tetragpyrmdA=',a2,i2,
     1       ' XAXb2 ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1, mx1,mx2,i
      write(ipun,210) -1, mx3,mx4,i
        write(iptyp,530) ncg,isymb(i),i,scale(ncg)
 530    format(i4,'. ','tetragpyrmdA=',a2,i2,
     1       ' XAXea ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200)  ncg,1,mx2,mx3,i
      write(ipun,210) -1, mx4,mx1,i
        write(iptyp,540) ncg,isymb(i),i,scale(ncg)
 540    format(i4,'. ','tetragpyrmdA=',a2,i2,
     1       ' XAXeb ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1,my,mx3,i
      write(ipun,210) 1,my,mx4,i
      write(ipun,210) -1,my,mx2,i
      write(ipun,210) -1,my,mx1,i
        write(iptyp,550) ncg,isymb(i),i,scale(ncg)
 550    format(i4,'. ','tetragpyrmdA=',a2,i2,
     1       ' YAXea ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1,my,mx2,i
      write(ipun,210) 1,my,mx3,i
      write(ipun,210) -1,my,mx1,i
      write(ipun,210) -1,my,mx4,i
        write(iptyp,560) ncg,isymb(i),i,scale(ncg)
 560    format(i4,'. ','tetragpyrmdA=',a2,i2,
     1       ' YAXeb ','sc=',f5.4,2x)
100   format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
110   format(' ',9x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
200   format(i4,2x,'K',3x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
210   format(' ',9x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
      return
      end
c
      Subroutine ax4y2(ipun,iptyp,ncg,i,mx1,mx2,mx3,mx4,my1,my2,
     1 isymb,scale)
      implicit real*8 (a-h,o-z)
c     this generates the coordinates around a hexavalent atom i,
c     for the tetragonal bipyramidal case
      character*2 isymb
      parameter (ztol=1.d-4)
      dimension scale(*),isymb(*)
c
c     first: rearrange the ligands in cyclic order
c     by finding the atom opposite to mx1 based on the largest angle:
      call angle(mx1,i,mx2,bet12)
      call angle(mx1,i,mx3,bet13)
      call angle(mx1,i,mx4,bet14)
      bmax=bet13
      icase=13
c in the default case 3. is opposite to 1.
c and indexes will not be changed; else:
      if((bet12-bmax).gt.ztol) then
         bmax=bet12
         icase=12
      elseif ((bet14-bmax).gt.ztol)then
        bmax=bet14
        icase=14
      endif
      if(icase.eq.12) then
        mxtemp=mx3
        mx3=mx2
        mx2=mxtemp
      elseif (icase.eq.14) then
        mxtemp=mx3
        mx3=mx4
        mx4=mxtemp
      endif
c
c  out1, (symm)
      ncg=ncg+1
      write(ipun,100) ncg,1,mx1,mx2,mx3,i
      write(ipun,110) 1,mx2,mx3,mx4,i
      write(ipun,110) 1,mx3,mx4,mx1,i
      write(ipun,110) 1,mx4,mx1,mx2,i
        write(iptyp,500) ncg,isymb(i),i,scale(ncg)
 500    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       ' OUTa1 ','sc=',f5.4,2x)
c  out2
      ncg=ncg+1
      write(ipun,100) ncg,1,mx1,mx2,mx3,i
      write(ipun,110) -1,mx2,mx3,mx4,i
      write(ipun,110) 1,mx3,mx4,mx1,i
      write(ipun,110) -1,mx4,mx1,mx2,i
        write(iptyp,510) ncg,isymb(i),i,scale(ncg)
 510    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       ' OUTb1 ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1, mx1,mx2,i
      write(ipun,210) -1, mx2,mx3,i
      write(ipun,210) +1, mx3,mx4,i
      write(ipun,210) -1, mx4,mx1,i
        write(iptyp,520) ncg,isymb(i),i,scale(ncg)
 520    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       ' XAXb2 ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1, mx1,mx2,i
      write(ipun,210) -1, mx3,mx4,i
        write(iptyp,530) ncg,isymb(i),i,scale(ncg)
 530    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       ' XAXea ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1, mx2,mx3,i
      write(ipun,210) -1, mx4,mx1,i
        write(iptyp,540) ncg,isymb(i),i,scale(ncg)
 540    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       ' XAXeb ','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1,my1,mx3,i
      write(ipun,210) 1,my1,mx4,i
      write(ipun,210) -1,my1,mx2,i
      write(ipun,210) -1,my1,mx1,i
      write(ipun,210) 1,my2,mx3,i
      write(ipun,210) 1,my2,mx4,i
      write(ipun,210) -1,my2,mx2,i
      write(ipun,210) -1,my2,mx1,i
        write(iptyp,550) ncg,isymb(i),i,scale(ncg)
 550    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       7H YAXea','sc=',f5.4,2x)
      ncg=ncg+1
      write(ipun,200) ncg,1,my1,mx3,i
      write(ipun,210) 1,my1,mx4,i
      write(ipun,210) -1,my1,mx2,i
      write(ipun,210) -1,my1,mx1,i
      write(ipun,210) -1,my2,mx3,i
      write(ipun,210) -1,my2,mx4,i
      write(ipun,210) +1,my2,mx2,i
      write(ipun,210) +1,my2,mx1,i
        write(iptyp,555) ncg,isymb(i),i,scale(ncg)
 555    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       7H YAXea",'sc=',f5.4,2x)
c
      ncg=ncg+1
      write(ipun,200) ncg,1,my1,mx2,i
      write(ipun,210) 1,my1,mx3,i
      write(ipun,210) -1,my1,mx1,i
      write(ipun,210) -1,my1,mx4,i
      write(ipun,210) 1,my2,mx2,i
      write(ipun,210) 1,my2,mx3,i
      write(ipun,210) -1,my2,mx1,i
      write(ipun,210) -1,my2,mx4,i
        write(iptyp,560) ncg,isymb(i),i,scale(ncg)
 560    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       7H YAXeb','sc=',f5.4,2x)
c
      ncg=ncg+1
      write(ipun,200) ncg,1,my1,mx2,i
      write(ipun,210) 1,my1,mx3,i
      write(ipun,210) -1,my1,mx1,i
      write(ipun,210) -1,my1,mx4,i
      write(ipun,210) -1,my2,mx2,i
      write(ipun,210) -1,my2,mx3,i
      write(ipun,210) +1,my2,mx1,i
      write(ipun,210) +1,my2,mx4,i
        write(iptyp,565) ncg,isymb(i),i,scale(ncg)
 565    format(i4,'. ','tetrbipyrmdA=',a2,i2,
     1       7H YAXeb",'sc=',f5.4,2x)
100   format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
110   format(' ',9x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
200   format(i4,2x,'K',3x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
210   format(' ',9x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
      return
      end
      Subroutine angle (ii,jj,kk,alpha)
      implicit real*8 (a-h,o-z)
c
c     this routine determines the angle between atoms ii-jj-kk
c
      dimension r1(3),r2(3),r3(3)
      call cartes (.false.,ii,r1(1),r1(2),r1(3),iz1)
      call cartes (.false.,jj,r2(1),r2(2),r2(3),iz1)
      call cartes (.false.,kk,r3(1),r3(2),r3(3),iz1)
      call vector (r1,r2,r1)
      call vector (r3,r2,r3)
      call nom (r1)
      call nom (r3)
      alpha=arcos(scalar(r1,r3))
      end
c
      FUNCTION ARCOS (X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (X.GE.1.0D0) GO TO 10
      IF (X.LE.-1.0D0) GO TO 20
      X1=SQRT(1.0D0-X**2)
      IF (ABS(X).LT.1.0D-11) GO TO 30
      S=ATAN(X1/X)
      IF (X.LT.0.0) S=S+3.1415926536D0
      ARCOS=S
      RETURN
   10 ARCOS=0.0D0
      RETURN
   20 ARCOS=3.1415926536D0
      RETURN
   30 ARCOS=1.5707963268D0
      RETURN
C
      END
      Subroutine cartes (in,i,xx,yy,zz,iz)
      implicit real*8 (a-h,o-z)
      logical in
      parameter (maxat=200)
      dimension x(3,maxat),ian(maxat)
      save x,ian
c     if in=true, x(1,i)=xx, x(2,i)=yy, x(3,i)=zz, ian(i)=iz
c     if in=false, xx=x(1,i),yy=x(2,i),zz=x(3,i),iz=ian(i)
      if (in) then
         x(1,i)=xx
         x(2,i)=yy
         x(3,i)=zz
         ian(i)=iz
      else
         xx=x(1,i)
         yy=x(2,i)
         zz=x(3,i)
         iz=ian(i)
      end if
      end
      Subroutine dihed (itm,i,ibd,noth,theta)
c     theta will be the cosine of the dihedral angle itm-i-ibd-noth
      implicit real*8 (a-h,o-z)
      dimension x1(3),x2(3),x3(3),x4(3),y1(3),y2(3)
      call cartes (.false.,itm, x1(1),x1(2),x1(3),iz)
      call cartes (.false.,i, x2(1),x2(2),x2(3),iz)
      call cartes (.false.,ibd, x3(1),x3(2),x3(3),iz)
      call cartes (.false.,noth,x4(1),x4(2),x4(3),iz)
      call vector (x1,x2,y1)
      call vector (x2,x3,y2)
      call normal  (y1,y2,x1)
      call vector (x2,x3,y1)
      call vector (x3,x4,y2)
      call normal (y1,y2,x4)
      theta = scalar (x1,x4)
      end
      Subroutine NOM (U)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3)
      X=1.0000000000/SQRT(SCALAR(U,U))
      DO 10 I=1,3
         U(I)=U(I)*X
   10 CONTINUE
      RETURN
C
      END
      Subroutine NORMAL (U,V,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3), V(3), W(3)
C
C     99999...  W WIRD EIN SENKRECHT AUF DIE EBENE(U,V) STEHENDER EINHEI
C      TOR
C
      W(1)=U(2)*V(3)-U(3)*V(2)
      W(2)=U(3)*V(1)-U(1)*V(3)
      W(3)=U(1)*V(2)-U(2)*V(1)
      CALL NOM (W)
      RETURN
C
      END
      FUNCTION S2 (X)
      IMPLICIT REAL*8 (A-H,O-Z)
      S2=SQRT(1.0-X**2)
      RETURN
C
      END
      FUNCTION SCALAR (U,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3), V(3)
      SCALAR=0.0
      DO 10 I=1,3
         SCALAR=SCALAR+U(I)*V(I)
   10 CONTINUE
      RETURN
C
      END
      Subroutine vector (r1,r2,r3)
      implicit real*8 (a-h,o-z)
c
c     r3 = r1 - r2 , all three 3-vectors
c
      dimension r1(3),r2(3),r3(3)
      do 100 i=1,3
         r3(i)=r1(i)-r2(i)
 100  continue
      end
