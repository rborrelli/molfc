      Subroutine subsets(natom,nbond,ibond,msub,indsub,nsub)
cversion 5.0
      implicit real*8(a-h,o-z)
       parameter(mxneig=6,maxat=200)
       dimension nbond(natom),ibond(mxneig,natom)
       dimension indsub(natom),nsub(natom)
       dimension mv(maxat)
c
c      output:
c      msub = number of subsystems found
c      indsub - indices of the atoms in the subsystems all collected
c      in one vector, the first nsub(1) defining subset 1., nsub(2)
c      subset 2., etc.
c
       do 10 j=1,natom
       mv(j)=0
       indsub(j)=j
       nsub(j)=0
  10   continue
       nmark=0
c   nmark will count the total number of markings
       isub=0
c this will count the number of subsystems
c
       do 90 jstart=1,natom
c   jstart is the start atom of the marking process
         if(mv(jstart).eq.1) goto 90
c   skip members of the previous subset
         isub=isub+1
         do 700 k=1,natom
           mv(k)=0
 700     continue
c  mv  for markers vector, will contain the indices of marked atoms
         mv(jstart)=1
         nmark=nmark+1
         indsub(nmark)=jstart
         nsub(isub)=nsub(isub)+1
c   dimension of the subset will be equal to the number of markers

 750     continue
         iadd=0
         do 900 k=1,natom
           if (mv(k).ne.1) then
              do 800 kk=1,nbond(k)
                 if (mv(ibond(kk,k)).eq.1.and.mv(k).ne.1) then
                     mv(k)=1
c
c                    if an atom has a marked neigbor, mark it
c  the .and. is to avoid that due to two marked neigbors an
c  atom will be marked twice - this would confuse the counting by nmark
c
                     iadd=1
                     nmark=nmark+1
                     indsub(nmark)=k
                     nsub(isub)=nsub(isub)+1
                  end if
 800          continue
           end if
 900     continue
        if (iadd.ne.0) go to 750
c     end of a subsystem
c
c     Now check the sum of subsystems:
        if(nmark.eq.natom) then
          goto 99
        endif
  90  continue
  99  continue
      msub=isub
c TEMP
c     print *, 'jstart=',jstart,'nmark=',nmark,'natom=',natom
c     print *, 'nsub'
c     print *, nsub
c     print *, 'indsub'
c     print *, indsub
C     print *, 'number of subsystems', msub
      return
      end
c
      Subroutine unique (nat,iat,n)
      implicit real*8 (a-h,o-z)
      dimension iat(*)
      dimension xx(3,6),iz(6),x(3),rr(6)
c      This routine tries to select the unique atom from the list
c      Usually, n is a central atom and iat(1..nat) is a set of
c      atoms bonded to it
c      The routine rearranges the atoms in iat so that the
c      unique center is now the first one
c      Strategy: First try if any of the atoms has a n atomic
c      number which is different from the rest. Choose the most
c      distinct one
c      Next, try the same with the bond lengths n...iat(k)
      call cartes(.false.,n,x(1),x(2),x(3),ii)
      aver=0.0d0
      do 100 i=1,nat
        call cartes(.false.,iat(i),xx(1,i),xx(2,i),xx(3,i),iz(i))
        aver=aver+iz(i)
       rr(i)=sqrt((xx(1,i)-x(1))**2+(xx(2,i)-x(2))**2+(xx(3,i)-x(3))**2)
 100  continue
      aver=aver/nat
      ximum=0.0d0
      iuniq=0
      do 200 i=1,nat
        if(abs(aver-iz(i)).gt.ximum) then
          ximum=abs(aver-iz(i))
          iuniq=i
        end if
 200  continue
        if (ximum.lt.1.0d-6.or.iuniq.eq.0) then
c       Try the bond lengths
        aver=0.0d0
        do 300 i=1,nat
          aver=aver+rr(i)
 300    continue
        aver=aver/nat
        ximum=0.0d0
        iuniq=0
        do 400 i=1,nat
          if(abs(aver-rr(i)).gt.ximum) then
            ximum=abs(aver-rr(i))
            iuniq=i
          end if
 400    continue
      end if
      if (iuniq.gt.1) then
        ii=iat(iuniq)
        iat(iuniq)=iat(1)
        iat(1)=ii
      end if
      end
      Subroutine tercoor (i,nterm,neigterm,nnon,neignon,ipun,
     1iptyp,iring,ncg,isymb,ibord,scale)
c     Internal coordinates for tertiary and quaternary atoms
      implicit real*8 (a-h,o-z)
      character*2 isymb
c      nterm - no of terminal neigbrs; nnon - nonterm.neigbrs
      parameter (planar=5.6d0)
      parameter(mxnr=40,mxdr=30,mxneig=6)
c     max no of rings, max dim of a ring, max no of neigbors of i
      dimension neigterm(*),neignon(*),iring(mxnr,*)
      dimension isymb(*), ibord(*), scale(*)
      dimension irnb(mxneig),inrnb(mxneig), ibb(mxneig)
c
c TEMP
      do 11 il=1,mxneig
       irnb(il)=0
       inrnb(il)=0
       ibb(il)=0
 11   continue
c TEMP end
      nbtot=nnon+nterm
      if(nbtot.gt.6) then
        print *, 'higher than hexavalent atom found in tercoor'
        print *, ' this is not yet implemented'
        call bummer('higher than hexavalent atom found',0,2)
        STOP
      endif
      call unique (nnon,neignon,i)
c     This selects one of the non-terminal neigbors as unique.
c     The unique atom will be in neignon(1)
      if (nnon.eq.3) then
         ib1=neignon(1)
         ib2=neignon(2)
         ib3=neignon(3)
         if (nterm.ne.0)  it1=neigterm(1)
      else
c   nnon is 4 here
         call unique (3,neignon(2),i)
         ib1=neignon(2)
         ib2=neignon(3)
         ib3=neignon(4)
         it1=neignon(1)
      end if
c     Determine if this is a planar center
      call angle (ib2,i,ib3,al1)
      call angle (ib1,i,ib3,al2)
      call angle (ib2,i,ib3,al3)
      sum=al1+al2+al3
c     This is the sum of the angles around the central atom
c     This is meaningful only for tertiary atoms
      isum=0
      do 50 j=1,mxnr
         isum=isum+iring(j,i)
 50   continue
c      isum gives the number of rings in which atom i  is present
      if (isum.eq.0) then
c      Non-ring central atom
c
        if(nbtot.gt.4) then
c FG: the following jumpout in penta- hexavalent case was added later:
c      (designed for non-ring system only)
        call pentahex(i,nterm,nnon,neigterm,neignon,ipun,
     1    iptyp,ncg,isymb,scale)
        RETURN
        endif
c
c   now the original, for max. tetravalent cases:
c        Test for planarity
C     write(*,*) 'angles',i,ib1,ib2,ib3,al1,al2,al3
c  FG
        if(nterm.eq.0.and.nnon.eq.3) then
c    XY3 case (where  the Y's are non-terminal atoms):
           if (sum.gt.planar) then
c  planar case
            ncg=ncg+1
            write (ipun,700) ncg,1,ib1,ib2,ib3,i
            write (ipun,800) 1,ib2,ib1,ib3,i
            write (ipun,800) 1,ib3,ib1,ib2,i
c            scale(ncg)=....
             write(iptyp,510) ncg,isymb(i),i,scale(ncg)
 510      format(i4,'. ','t-planXY3umbr',a2,i2,' SOUT  ','sc=',f5.4,2x)
 700     format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
 800     format(' ',9x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
c

           else
c          Non-planar case
            ncg=ncg+1
            write (ipun,100) ncg,1,ib2,ib3,i
            write (ipun,200) 1,ib1,ib3,i
            write (ipun,200) 1,ib1,ib2,i
c
c            scale(ncg)=....
            write(iptyp,520) ncg,isymb(i),i,scale(ncg)
 520      format(i4,'. ','tert XY3 , X=',a2,i2,' SDEF  ','sc=',f5.4,2x)
c
 100     format(i4,2x,'K',3x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
 200     format(' ',9x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
c     Symmetric def
           end if
c     end of XY3 case
c
        elseif(nterm.eq.1.or.nnon.eq.4) then
c
c  FG   XY3T (or XY3Z) case(Y's are non-term., T term. (or Z non-term.)
c   the two cases are treated equally;
c
           ncg=ncg+1
           write(ipun,100) ncg,1,ib2,ib3,i
           write(ipun,200) 1,ib1,ib3,i
           write(ipun,200) 1,ib1,ib2,i
           write(ipun,200) -1,ib1,it1,i
           write(ipun,200) -1,ib2,it1,i
           write(ipun,200) -1,ib3,it1,i
c
c            scale(ncg)=....
             if(nnon.ne.4) then
             write(iptyp,610) ncg,isymb(i),i,scale(ncg)
 610         format(i4,'. ','tert XY3T, X=',a2,i2,' SDEF  ',
     1       'sc=',f5.4,2x)
             else
             write(iptyp,612) ncg,isymb(i),i,scale(ncg)
 612         format(i4,'. ','quat XY3Z, X=',a2,i2,' SDEF  ',
     1       'sc=',f5.4,2x)
             endif
c
        else
          write(*,*) 'impossible number of atoms around a'
          write(*,*) 'tertiary or quaternary atom in tercoor'
          call bummer('tertiary or quaternary atom in tercoor',0,2)
          STOP
        endif
c  symm def above for XY3T or XY3Z, treated equally
c
c   the following as def is the same for XY3 or XY3T:
         ncg=ncg+1
         write (ipun,100) ncg,2,ib2,ib3,i
         write (ipun,200) -1,ib1,ib3,i
         write (ipun,200) -1,ib1,ib2,i
c
c           scale(ncg)=....
            if(nnon.ne.4) then
              if(nterm.eq.1) then
               write(iptyp,620) ncg,isymb(i),i,scale(ncg)
 620           format(i4,'. ','tert XY3T, X=',a2,i2,' ADEFa ',
     1         'sc=',f5.4,2x)
              else
               write(iptyp,621) ncg,isymb(i),i,scale(ncg)
 621           format(i4,'. ','tert XY3 , X=',a2,i2,' ADEFa ',
     1         'sc=',f5.4,2x)
              endif
            else
            write(iptyp,622) ncg,isymb(i),i,scale(ncg)
 622        format(i4,'. ','quat XY3Z, X=',a2,i2,' ADEFa ',
     1      'sc=',f5.4,2x)
            endif
c      As def (a)
c
         ncg=ncg+1
         write (ipun,100)  ncg,1,ib1,ib3,i
         write (ipun,200) -1,ib1,ib2,i
c           scale(ncg)=....
            if(nnon.ne.4) then
              if(nterm.eq.1) then
               write(iptyp,630) ncg,isymb(i),i,scale(ncg)
 630           format(i4,'. ','tert XY3T, X=',a2,i2,' ADEFb ',
     1         'sc=',f5.4,2x)
              else
               write(iptyp,631) ncg,isymb(i),i,scale(ncg)
 631           format(i4,'. ','tert XY3 , X=',a2,i2,' ADEFb ',
     1         'sc=',f5.4,2x)
              endif
            else
            write(iptyp,632) ncg,isymb(i),i,scale(ncg)
 632        format(i4,'. ','quat XY3Z, X=',a2,i2,' ADEFb ',
     1      'sc=',f5.4,2x)
            endif
c      As def (b)
c
        if (nterm.eq.1.or.nnon.eq.4) then
c          XY3T ---  add the motion of the terminal atom
          ncg=ncg+1
          write (ipun,100) ncg,2,ib1,it1,i
          write (ipun,200) -1,ib2,it1,i
          write (ipun,200) -1,ib3,it1,i
c
c           scale(ncg)=....
            if(nnon.ne.4) then
            write(iptyp,640) ncg,isymb(i),i,scale(ncg)
 640        format(i4,'. ','tert XY3T, X=',a2,i2,' ROCKa ',
     1      'sc=',f5.4,2x)
            else
            write(iptyp,642) ncg,isymb(i),i,scale(ncg)
 642        format(i4,'. ','quat XY3Z, X=',a2,i2,' ROCKa ',
     1      'sc=',f5.4,2x)
            endif
c       rocking (a)
c
          ncg=ncg+1
          write (ipun,100)  ncg,1,ib2,it1,i
          write (ipun,200) -1,ib3,it1,i
c
c           scale(ncg)=....
            if(nnon.ne.4) then
            write(iptyp,650) ncg,isymb(i),i,scale(ncg)
 650        format(i4,'. ','tert XY3T, X=',a2,i2,' ROCKb ',
     1      'sc=',f5.4,2x)
            else
            write(iptyp,652) ncg,isymb(i),i,scale(ncg)
 652        format(i4,'. ','quat XY3Z, X=',a2,i2,' ROCKb ',
     1      'sc=',f5.4,2x)
            endif
c       rocking (b)
        end if
      else
c       i is a ring atom
c   this is special in that the ring coordinates were already
c  generated in ring, thus only the atoms attached to
c     the ring will be handled here:
c
c         Locate the atoms in the same ring as i
        kk1=0
        do 260 kk=1,nnon
            kk1=kk1+1
            ibb(kk1)=neignon(kk)
 260    continue
        do 270 kk=1,nterm
           kk1=kk1+1
           ibb(kk1)=neigterm(kk)
 270    continue
        nnei=kk1
c  nnei is the total number of neigbors, including terminal ones
c
c FG  check if all neighbors are ring atoms; this means condensed
c  rings or spiro-systems where the relative motions of rings
c   are treated separately
c
        call ringmemb(mxneig,mxnr,nnei,ibb,iring,nrgnb,nonrgnb,
     1    irnb,inrnb)
c
c   (nrgnb number of ring neigbors stored in irnb)
        jspec=inrnb(1)
c
        if(nrgnb.eq.nnei) then
c
c     all neighbours are ring atoms,
c     condensed rings and spiro systems are treated separately
c
          RETURN
c... just 2 special cases will be handled here:
        else if(nrgnb.eq.3.and.nnei.eq.4) then
c one special case is bicyclo [n.m.0] compounds, with one terminal
c atom (H) attached to atom i, the annelation tertiary carbon:
          ncg=ncg+1
          write(ipun,100) ncg,2,jspec,irnb(1),i
          write(ipun,200) -1,jspec,irnb(2),i
          write(ipun,200) -1,jspec,irnb(3),i
           write(iptyp,712) ncg,isymb(jspec),jspec,scale(ncg)
 712    format(i4,'. ','2ringcornr-Y=',a2,i2,' ADEFa ','sc=',f5.4,2x)
            ncg=ncg+1
            write(ipun,100)  ncg,1,jspec,irnb(2),i
            write(ipun,200) -1,jspec,irnb(3),i
           write(iptyp,722) ncg,isymb(jspec),jspec,scale(ncg)
 722      format(i4,'. ','2ringcornr-Y=',a2,i2,' ADEFb ','sc=',f5.4,2x)
           RETURN
c......
        else if(nrgnb.eq.4.and.nterm.eq.1) then
c
c   a second special case: think of a silatrane-type
c   3-ring system  with a N..Si-c1c2c3-Y corner,where we want to
c   describe the motion of  Y:
c
c
c   (there is only one non-rg neigb -Y- this is the first element)
c   ringmemb has put N  - ngb of highest 'ring-multiplicity'
c   (member of the highest number of rings) into
c   irnb(1), thus the carbons start with irnb(2)
c
          ncg=ncg+1
          write(ipun,100) ncg,2,jspec,irnb(2),i
          write(ipun,200) -1,jspec,irnb(3),i
          write(ipun,200) -1,jspec,irnb(4),i
           write(iptyp,710) ncg,isymb(jspec),jspec,scale(ncg)
 710    format(i4,'. ','3ringcornr-Y=',a2,i2,' ADEFa ','sc=',f5.4,2x)
c
          ncg=ncg+1
          write(ipun,100)  ncg,1,jspec,irnb(3),i
          write(ipun,200) -1,jspec,irnb(4),i
           write(iptyp,720) ncg,isymb(jspec),jspec,scale(ncg)
 720      format(i4,'. ','3ringcornr-Y=',a2,i2,' ADEFb ','sc=',f5.4,2x)
          RETURN
        endif
c
        if(nbtot.gt.4) then
          print *, 'higher than tetravalent ring-atom found in'
          print *, 'tercoor, and this is not a silatrane case'
          print *, 'not yet implemented'
          call bummer('higher than tetravalent ring-atom found',0,2)
          STOP
        endif
c end FG.
c        Now the original treatment for ring-atom i:
c        The part below was apparently designed for isolated rings
c        The array ibb contains simply the nnei neighbors of i
c        Rearrange it so that the first two will be the ring connections
        kk1=0
        do 300 k=1,mxnr
           if (iring(k,i).eq.0) go to 300
           do 290 kk=1,nnei
              if (iring(k,ibb(kk)).ne.1) go to 290
c             The neighbor kk is in the same ring as i
              kk1=kk1+1
              itmp=ibb(kk1)
              ibb(kk1)=ibb(kk)
              ibb(kk)=itmp
 290       continue
 300    continue
        if (nnei.eq.3) then
           ncg=ncg+1
           write (ipun,100) ncg,1,ibb(1),ibb(3),i
           write (ipun,200) -1,ibb(2),ibb(3),i
             write(iptyp,810) ncg, isymb(ibb(3)),ibb(3),isymb(i),
     1       i,scale(ncg)
 810         format(i4,'. ','ring tXY',a2,i2,'-',a2,i2,' ROCK  ',
     1       'sc=',f5.4,2x)
c          Rocking
c
           ncg=ncg+1
           write (ipun,700) ncg,1,ibb(3),ibb(1),ibb(2),i
             write(iptyp,820) ncg, isymb(ibb(3)),ibb(3),isymb(i),
     1       i,scale(ncg)
 820         format(i4,'. ','ring tXY',a2,i2,'-',a2,i2,' WAGG  ',
     1       'sc=',f5.4,2x)
c          Out-of-plane wagging
c
        else
c FG      the ring atom i is quaterner, i.e. has two neigbors outside
c         the ring:
          ncg=ncg+1
          write(ipun,100) ncg,4,ibb(3),ibb(4),i
          write (ipun,200) -1,ibb(3),ibb(1),i
          write (ipun,200) -1,ibb(3),ibb(2),i
          write (ipun,200) -1,ibb(4),ibb(1),i
          write (ipun,200) -1,ibb(4),ibb(2),i
           write(iptyp,910) ncg, isymb(i),i,scale(ncg)
 910       format(i4,'. ','ring tXY2, X=',a2,i2,' SCIS  ',
     1     'sc=',f5.4,2x)
c      symm. def of the XY2 group
c
          ncg=ncg+1
          write (ipun,100)  ncg,1,ibb(3),ibb(1),i
          write (ipun,200)  1,ibb(3),ibb(2),i
          write (ipun,200) -1,ibb(4),ibb(1),i
          write (ipun,200) -1,ibb(4),ibb(2),i
           write(iptyp,920) ncg, isymb(i),i,scale(ncg)
 920       format(i4,'. ','ring tXY2, X=',a2,i2,' ROCK  ',
     1     'sc=',f5.4,2x)
c      XY2 rocking
          ncg=ncg+1
          write (ipun,100)  ncg,1,ibb(3),ibb(1),i
          write (ipun,200) -1,ibb(3),ibb(2),i
          write (ipun,200)  1,ibb(4),ibb(1),i
          write (ipun,200) -1,ibb(4),ibb(2),i
           write(iptyp,930) ncg, isymb(i),i,scale(ncg)
 930       format(i4,'. ','ring tXY2, X=',a2,i2,' WAGG  ',
     1     'sc=',f5.4,2x)
c      XY2 wagging
c
          ncg=ncg+1
          write (ipun,100)  ncg,1,ibb(3),ibb(1),i
          write (ipun,200) -1,ibb(3),ibb(2),i
          write (ipun,200) -1,ibb(4),ibb(1),i
          write (ipun,200)  1,ibb(4),ibb(2),i
           write(iptyp,940) ncg, isymb(i),i,scale(ncg)
 940       format(i4,'. ','ring tXY2, X=',a2,i2,' TWIST ',
     1     'sc=',f5.4,2x)
c      XY2 twisting
c
 900     format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
        end if
      end if
c END of TERCOOR
      end
c
      Subroutine ringmemb(mxneig,mxnr,nbvec,ibvec,iring,
     1nyes,nnon,iyes,inon)
c
c     given the indices of some atoms in the array ibvec,
c     establish for each atom the number of rings in which it
c     is a member -let's call it 'ring multiplicity'
c     put the ring atoms in iyes, with iyes(1) the atom of
c     highest ring multiplicity
c     specific use:
c     find those neigbors of an atom 'i' which are ring-atoms
c     (in polycyles, each may be in different rings)
c     atom i  is not explicitely defined, but this usage
c     tacitly assumes that the list ibvec refers to a given atom
c input:   nbvec = number of neigbors of atom i
c          ibvec - indices of the neigbors
c          iring(p,l) is 1 if atom l is part of ring p
c output:  nyes = no of ring neigbors; nnon = non-rg neigbors
c          iyes - indices of ring neigbors; inon - non-rg. ngb.
c             (iyes(1) has the highest ring-multiplicity)
c
      implicit real*8(a-h,o-z)
      dimension ibvec(nbvec),iring(mxnr,*),iyes(mxneig),inon(mxneig)
      nyes=0
      nnon=0
      mxsum=0
C   problem in this routine, if compilation with optimizer;
C TEMP try this, giving start value to jmax(but this never makes sense):
      jmax=ibvec(1)
      do 50 k=1,mxneig
        iyes(k)=0
        inon(k)=0
 50   continue
      do 100 k=1,nbvec
        j=ibvec(k)
        jsum=0
        do 200 ir=1,mxnr
          jsum=jsum+iring(ir,j)
200     continue
        if(jsum.gt.0) then
          nyes=nyes+1
          if(nyes.gt.mxneig) then
            print *, 'impossible number of neigbors in ringnb'
            call bummer('impossible number of neighbours',0,2)
            STOP
          endif
          iyes(nyes)=j
          if(jsum.gt.mxsum) then
            mxsum=jsum
            jmax=j
          endif
        else
c            j is not a ring-atom
          nnon=nnon+1
          inon(nnon)=j
        endif
100   continue
      if(nyes.gt.0) then
        if (iyes(1).ne.jmax) then
c          put the ring-atom with highest 'ring-multiplicity' into the
c          first element:
          do 300 k=2,nyes
            if(iyes(k).eq.jmax) then
              iyes(k)=iyes(1)
              iyes(1)=jmax
              goto 310
            endif
300     continue
        endif
      endif
310   continue
      return
      end
      Subroutine order (natom,nterm,iterm,nbond,ibond,ipun,iptyp,iring,
     1 ncent,icent,nprim,iprim,nsec,isec,ntert,itert,nquat,iquat,
     2 ncg,isymb,ibord,scale,dummy)
c     This routine calculates the number and indices of primary,
c       secondary , tertiary and quaternary atoms
c     Input: natom=number of atoms
c            nterm=number of terminal atoms, iterm(1..neterm)=
c            the indices of terminal atoms
c            nbond(1.. natom)= the number of bonds an atom participates
c            in; ibond(1..6,1.. natom)= their numbers
c            ipun=punch file
c            iring(i,j) is 1 if atom j is in ring i
c     Output: ncent= number of central atoms (in small molecules)
c             icent(1..ncent) = their serial numbers
c             nprim=number of primary atoms
c             iprim(1..nprim) = their indices, etc. for secondary
c             tertiary & quaternary atoms
      implicit real*8 (a-h,o-z)
      character*2 isymb
      logical dummy
      parameter (mxnr=40)
      parameter(mxneig=6)
      dimension nbond(natom),ibond(mxneig,natom)
      dimension iprim(*),isec(*),itert(*),iquat(*)
      dimension icent(*),iterm(*),iring(mxnr,*)
      dimension isymb(natom), ibord(*), scale(*)
      dimension neignon(mxneig),neigterm(mxneig)
      ncent=0
      nprim=0
      nsec=0
      ntert=0
      nquat=0
      do 500 i=1,natom
        do 100 j=1,nterm
           if (iterm(j).eq.i)  go to 500
c          exclude terminal atoms
 100    continue
        neigh=0
        nnterm=0
        nnon=0
c       count the non-terminal neighbors ONLY
        do 300 k=1,nbond(i)
c         scan the neighbors of i
          kk=ibond(k,i)
c         kk is a neighbor of i
          do  200 l=1,nterm
             if (iterm(l).eq.kk) then
                nnterm=nnterm+1
                neigterm(nnterm)=kk
                go to 300
             end if
c            leave out terminal neighbors
 200      continue
          nnon=nnon+1
          neignon(nnon)=kk
          neigh=neigh+1
 300    continue
        if (neigh.eq.0) then
           ncent=ncent+1
           icent(ncent)=i
           call centcoor (i,nnterm,neigterm,nnon,neignon,ipun,iptyp,
     1     natom,ncg,isymb,ibord,scale,dummy)
        else if (neigh.eq.1) then
           nprim=nprim+1
           iprim(nprim)=i
           call primcoor (i,nnterm,neigterm,nnon,neignon,ipun,iptyp,
     1       nbond,ibond,natom,ncg,isymb,ibord,scale,dummy)
        else if (neigh.eq.2) then
           nsec=nsec+1
           isec(nsec)=i
           call secoor (i,nnterm,neigterm,nnon,neignon,ipun,iptyp,iring,
     1ncg,natom,nbond,ibond,isymb,ibord,scale)
        else if (neigh.ge.3) then
           ntert=ntert+1
           itert(ntert)=i
           call tercoor (i,nnterm,neigterm,nnon,neignon,ipun,iptyp,
     1iring,ncg,isymb,ibord,scale)
        else
C     write(*,400)  i, (ibond(ll,i),ll=1,nbond(i))
 400       format (' atom higher than tetravalent found, atom=',
     1      i4, ' bonds=',6i3)
        end if
 500  continue
c TEMP
C      print *,'ncent,nprim,etc.', ncent,nprim,nsec,ntert
C     if (ncent.ne.0) write (*,550) (icent(k),k=1,ncent)
 550  format (' central atoms',/,3x,(5i4,3x,5i4,/))
C     if (nprim.ne.0) write (*,600) (iprim(k),k=1,nprim)
 600  format (' primary atoms',/,3x,(5i4,3x,5i4,/))
C     if (nsec.ne.0) write (*,700) (isec(k),k=1,nsec)
 700  format (' secondary atoms',/,3x,(5i4,3x,5i4,/))
C     if (ntert.ne.0) write (*,800) (itert(k),k=1,ntert)
 800  format (' tertiary or quaternary atoms',/,3x,(5i4,3x,5i4,/))
      end
c
      Subroutine hring (natom,nbond,ibond,nbondo,ibondo,ipun,iptyp,
     1nrings,iring,ncg,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,
     2aromat)
      implicit real*8 (a-h,o-z)
c   FG, Aug, 94:
c   hring is the  modification of my ring routine
c   by Richi(Richard Hargitai)
c   Richi died a few weeks ago at the age of 29,
c   one of the nicest guys I have ever known!!
      character*2 isymb
      character*4 aromat
      logical inextok
c     This routine locates rings, builds the internal coordinates
c       for them, an puts markers in iring to show that a given atom
c        is part of a ring
c      Input: natom=number of atoms
c             nbond(i): the number of bonds for atom i AFTER PRUNING
c             ibond(1..nbond(i),i): the atoms i is bound AFTER PRUNING
c             nbondo(i): the number of bonds for atom i originally
c             ibondo(1..nbondo(i),i): the atoms i is bound originally
c..FG
c             for subr fring, called by ringcoor, we need also:
c             nb = no of bonds in the whole molecule (from subr bonds)
c             ian(1..natom) - atomic numbers
c             li,lj,rr(1..nb) - indices of atoms in bond ij, and
c                bond length, all from subr bonds
c             fcb(1..nb) force constants for bonds,from subr fbonds
c...
c      Output: nrings= number of rings
c              iring(i,j) is 1 if atom j is in ring i
c
c FG  nbond and ibond are modified within this routine; I need them,
c  however, in subr polycycl: to avoid rewriting, they will be
c  stored in nbondrg and ibondrg, and the latter put back into
c  nbond and ibond at returning. These are then the original lists
c for the ring (ring only) atoms(terminal atoms etc pruned earlier
c  in these lists,of course)
c
      parameter(mxneig=6,maxat=200)
c        max no of neighbors, max no of atoms
      parameter(mxdr=30,mxnr=40)
c        max dim of a ring; max number of rings
      dimension  nbond(natom),ibond(mxneig,natom),nbondo(natom),
chr+
     1  ibondo(mxneig,natom), iring(mxnr,*),itemp(mxdr),jtemp(mxdr),
     2  mat1(maxat,maxat),mat2(maxat,maxat),ito(maxat),iback(maxat)
      integer pot1(maxat),pot2(maxat)
      logical log1(maxat),log2(maxat)
chr-
      dimension ian(natom),li(nb),lj(nb),rr(nb),fcb(nb)
      dimension isymb(natom)
      dimension fcq(*), scale(*)
chr+
      logical ffring,fanel
      dimension ipath(mxneig,2,maxat)
chr-
c FG to store the original ring-atom lists:
      dimension nbondrg(maxat),ibondrg(mxneig,maxat)
c FG number of connections for the current ring atoms in itemp
c will also be needed in polycycles below
      dimension nbonditp(mxdr)
      do 10 k=1,mxdr
      nbonditp(k)=0
 10   continue
      do 20 i=1,natom
        nbondrg(i)=nbond(i)
        do 20 k=1,mxneig
          ibondrg(k,i)=ibond(k,i)
chr+
          ipath(k,1,i)=0
          ipath(k,2,i)=0
chr-
 20   continue
      nrings=0
      do 50 i=1,natom
         do 50 k=1,mxnr
         iring(k,i)=0
 50   continue
      iat1=1
 80   continue
c     Return address for renewed search
c
c
c TEMPORARY
cD      print *, 'varying bondlist ibond in subr ring'
cD      do 250 j=1,natom
cD        write (*,3000) j,nbond(j),(ibond(k,j),k=1,nbond(j))
cD250   continue
cD3000  format(1x,2i3,3x,6i4)
c TEMPOR vege
c
chr+
cd     write(*,*) 'Number of bonds: ',nb
111   fanel=.false.
      ffring=.false.
      do i=1,natom
        if (nbond(i).gt.2) then
          fanel=.true.
          do k=1,nbond(i)
             j=1
             lo=i
             itemp(1)=lo
             l=ibond(k,i)
             do while (nbond(l).eq.2)
               j=j+1
               itemp(j)=l
               if (ibond(1,l).ne.lo) then
                 lo=l
                 l=ibond(1,l)
               else
                 lo=l
                 l=ibond(2,l)
               endif
             enddo
             ipath(k,1,i)=l
             ipath(k,2,i)=j
             if (l.eq.i) then
               ffring=.true.
C              write(*,133) j,(itemp(l),l=1,j)
 133           format('xx ',i2,' membered Ring found: ',(10i3))
c onmagaba visszatero hurok
               if (j.ge.3) then
                 nrings=nrings+1
                 do l=1,j
                   iring(nrings,itemp(l))=1
                 enddo
               else
                 write(*,*) 'Caution! Two membered ring were found '
               endif
               if (j.gt.3) then
                 do l=1,j
                   jtemp(l)=itemp(l)
                 enddo
                 call ringcoor (j,jtemp,nbondo,ibondo,ipun,iptyp,ncg,
     1natom,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
cd                write(*,*) 'ringcoor over'
               endif
               itemp(j+1)=i
               call hclear(j,itemp,nbond,ibond,natom)
               goto 111
             endif
          enddo
cd         write(*,*) 'Path ',i
cd         write(*,333) (ipath(k,1,i),k=1,mxneig)
cd         write(*,333) (ipath(k,2,i),k=1,mxneig)
cd 333     format(6i3)
          do k=1,nbond(i)-1
            do j=1,nbond(i)
              itemp(j)=-1
            enddo
            l=ipath(k,1,i)
            itemp(k)=ipath(k,2,i)
            do j=k+1,nbond(i)
              if (l.eq.ipath(j,1,i)) then
                 ffring=.true.
                 itemp(j)=ipath(j,2,i)
c ket anellacios pont kozotti hurok
              endif
            enddo
            if (ffring) then
              l=100
              m=0
              l2=100
              n=0
              do j=1,nbond(i)
                if ((itemp(j).lt.l).and.(itemp(j).ge.0)) then
                  l2=l
                  n=m
                  l=itemp(j)
                  m=j
                else
                  if ((itemp(j).lt.l2).and.(itemp(j).gt.0)) then
                    n=j
                    l2=itemp(j)
                  endif
                endif
              enddo
              itemp(1)=i
              itemp(l+l2+1)=i
              ii=i
              do j=2,l+1
                jj=ibond(m,ii)
                itemp(j)=jj
                if (ibond(1,jj).eq.ii) then
                  m=2
                  ii=jj
                else
                  m=1
                  ii=jj
                endif
              enddo
              ii=i
              m=n
              do j=l+l2,l+2,-1
                jj=ibond(m,ii)
                itemp(j)=jj
                if (ibond(1,jj).eq.ii) then
                  m=2
                  ii=jj
                else
                  m=1
                  ii=jj
                endif
              enddo
              j=l+l2
C              write(*,134) j,(itemp(l),l=1,j)
 134           format('xx ',i2,' membered rIng found: ',(10i3))
               if (j.ge.3) then
                 nrings=nrings+1
                 do l=1,j
                   iring(nrings,itemp(l))=1
                 enddo
               else
                 write(*,*) 'Caution! Two membered ring were found '
               endif
               if (j.gt.3) then
                 do l=1,j
                   jtemp(l)=itemp(l)
                 enddo
                 call ringcoor (j,jtemp,nbondo,ibondo,ipun,iptyp,ncg,
     1natom,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
cd              write(*,*) 'ringcoor over'
               endif
               j=j-l2
               do l=1,l2+1
                 itemp(l)=itemp(j+l)
               enddo
               call hclear(l2,itemp,nbond,ibond,natom)
               goto 111
            endif
          enddo
        endif
      enddo
      if (fanel) then
c        write(*,*) 'Ughh'
        i=0
        do j=1,natom
          if (nbond(j).ge.3) then
            i=i+1
            ito(i)=j
            iback(j)=i
          endif
        enddo
        nane=i
        do i=1,nane
          do j=1,nane
            mat1(i,j)=1000
            mat2(i,j)=1000
          enddo
cd         write(*,*) ito(i)
        enddo
        mxpath=0
        ista=0
        iend=0
        do i=1,nane
          j=ito(i)
          k=nbond(j)
          do l=1,k
            i1=iback(ipath(l,1,j))
            i2=ipath(l,2,j)
            if (i2.gt.mxpath) then
              ista=i
              iend=i1
              mxpath=i2
            endif
            mat1(i,i1)=i2
            mat2(i,i1)=i2
          enddo
        enddo
        mat1(ista,iend)=1000
        mat2(ista,iend)=1000
        mat1(iend,ista)=1000
        mat2(iend,ista)=1000
cd       write(*,*) ista,iend
cd       do i=1,nane
cd         write(*,*) (mat1(i,j),j=1,nane)
cd       enddo
        do i=1,nane
          log1(i)=.false.
          log2(i)=.false.
          pot1(i)=0
          pot2(i)=0
        enddo
        log1(ista)=.true.
        log2(iend)=.true.
        call hsea(nane,mat1,log1,pot1,iend,ffring)
        call hsea(nane,mat2,log2,pot2,ista,ffring)
cd       do i=1,nane
cd         write(*,*) log1(i),pot1(i),log2(i),pot2(i)
cd       enddo
        j=pot1(iend)
        do i=1,nane
          log1(i)=log1(i).and.log2(i)
          if (log1(i)) then
            if ((pot1(i)+pot2(i)).ne.j) then
              log1(i)=.false.
            endif
          endif
        enddo
cd       do i=1,nane
cd         write(*,*) ito(i),log1(i),pot1(i)
cd       enddo
        j=ito(iend)
        jj=ito(ista)
        itemp(1)=j
        jtemp(1)=j
        do i=1,nbond(j)
          if (ipath(i,1,j).eq.jj) then
            k=i
          endif
        enddo
        l=1
cd       write(*,*) jtemp(l)
        do while (j.ne.jj)
          j=ibond(k,j)
          if (itemp(l).eq.ibond(1,j)) then
            k=2
          else
            k=1
          endif
          l=l+1
          itemp(l)=j
          jtemp(l)=j
cd       write(*,*) jtemp(l)
        enddo
        l1=l-1
        if (.not.ffring) then
          call hclear(l1,itemp,nbond,ibond,natom)
          goto 111
        endif
        do while (ista.ne.iend)
          i=1
cd         write(*,*) ista,iend
cd         do i=1,nane
cd           write(*,*) (mat1(i,j),j=1,nane)
cd         enddo
          do while (.not.((log1(i)).and.(mat1(ista,i).eq.0).and.
     1(pot1(i).gt.pot1(ista))))
            i=i+1
          enddo
          j=ito(ista)
          jj=ito(i)
          ista=i
          do i=1,nbond(j)
            if (ipath(i,1,j).eq.jj) then
              k=i
            endif
          enddo
          do while (j.ne.jj)
            j=ibond(k,j)
            if (jtemp(l).eq.ibond(1,j)) then
              k=2
            else
              k=1
            endif
            l=l+1
            jtemp(l)=j
cd       write(*,*) jtemp(l)
          enddo
        enddo
        l=l-1
C       write(*,135) l,(jtemp(i),i=1,l)
 135    format('xx ',i2,' membered rinG found: ',(10i3))
        if (l.ge.3) then
          nrings=nrings+1
          do i=1,l
            iring(nrings,jtemp(i))=1
          enddo
        else
          write(*,*) 'Caution! Two membered ring were found '
        endif
        if (l.gt.3) then
          call ringcoor (l,jtemp,nbondo,ibondo,ipun,iptyp,ncg,
     1natom,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
cd         write(*,*) 'ringcoor over'
        endif
        call hclear(l1,itemp,nbond,ibond,natom)
        goto 111
      else
        ffring=.false.
        do i=1,natom
          if (nbond(i).eq.2) then
            ffring=.true.
            itemp(1)=i
            ii=ibond(1,i)
            jj=i
            j=1
            do while (ii.ne.i)
              j=j+1
              itemp(j)=ii
              if (ibond(1,ii).eq.jj) then
                jj=ii
                ii=ibond(2,jj)
              else
                jj=ii
                ii=ibond(1,jj)
              endif
            enddo
            itemp(j+1)=i
C           write(*,137) j,(itemp(l),l=1,j)
 137           format('xx ',i2,' membered riNg found: ',(10i3))
            if (j.ge.3) then
              nrings=nrings+1
              do l=1,j
                iring(nrings,itemp(l))=1
              enddo
            else
              write(*,*) 'Caution! Two membered ring were found '
            endif
            if (j.gt.3) then
              do l=1,j
                jtemp(l)=itemp(l)
              enddo
              call ringcoor (j,jtemp,nbondo,ibondo,ipun,iptyp,ncg,
     1natom,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
cd              write(*,*) 'ringcoor over'
            endif
            call hclear(j,itemp,nbond,ibond,natom)
          endif
        enddo
      endif
      do i=1,natom
        nbond(i)=nbondrg(i)
        do k=1,mxneig
          ibond(k,i)=ibondrg(k,i)
        enddo
      enddo
      return
      end
c
      Subroutine hclear(j,itemp,nbond,ibond,natom)
      parameter(mxneig=6,maxat=200)
c        max no of neighbors, max no of atoms
      parameter(mxdr=30,mxnr=40)
c        max dim of a ring; max number of rings
      dimension  nbond(natom),ibond(mxneig,natom),itemp(mxdr)

      ii=itemp(1)
      jj=itemp(2)
      k=nbond(ii)
      n=0
      do i=k,1,-1
        m=ibond(i,ii)
        if (m.eq.jj) then
          ibond(i,ii)=n
          goto 111
        else
          ibond(i,ii)=n
          n=m
        endif
      enddo
111   nbond(ii)=k-1
      ii=itemp(j+1)
      jj=itemp(j)
      k=nbond(ii)
      n=0
      do i=k,1,-1
        m=ibond(i,ii)
        if (m.eq.jj) then
          ibond(i,ii)=n
          goto 113
        else
          ibond(i,ii)=n
          n=m
        endif
      enddo
113   nbond(ii)=k-1
      do l=2,j
        jj=itemp(l)
        do i=1,nbond(jj)
          ibond(i,jj)=0
        enddo
        nbond(jj)=0
      enddo
      do i=1,natom
        if (nbond(i).eq.1) then
          nbond(i)=0
          m=ibond(1,i)
          ibond(1,i)=0
          n=i
          do while (nbond(m).eq.2)
            nbond(m)=0
            if (ibond(1,m).eq.n) then
              n=m
              m=ibond(2,n)
            else
              n=m
              m=ibond(1,n)
            endif
            ibond(1,n)=0
            ibond(2,n)=0
          enddo
          jj=n
          ii=m
          k=nbond(ii)
          n=0
          do l=k,1,-1
            m=ibond(l,ii)
            if (m.eq.jj) then
              ibond(l,ii)=n
              goto 119
            else
              ibond(l,ii)=n
              n=m
            endif
          enddo
119       nbond(ii)=k-1
        endif
      enddo
cd     write(*,*) 'hclear : end'
cD      print *, 'varying bondlist ibond in subr ring'
cD      do 250 j=1,natom
cD        write (*,3000) j,nbond(j),(ibond(k,j),k=1,nbond(j))
cD250   continue
cD3000  format(1x,2i3,3x,6i4)
      return
      end

      Subroutine hsea(nane,mat,log,pot,iend,ffring)
      parameter(mxneig=6,maxat=200)
c        max no of neighbors, max no of atoms
      parameter(mxdr=30,mxnr=40)
c        max dim of a ring; max number of rings
      dimension mat(maxat,maxat)
      integer pot(maxat)
      logical log(maxat),ffring
      ffring=.true.
      do while (.not.log(iend))
        imin=1000
        i1=0
        i2=0
        do i=1,nane
          if (log(i)) then
            do j=1,nane
              if (.not.log(j)) then
                if (mat(i,j).lt.imin) then
                  imin=mat(i,j)
                  i1=i
                  i2=j
                endif
              endif
            enddo
          endif
        enddo
        if (imin.gt.800) then
          ffring=.false.
          goto 555
        endif
        do i=1,nane
          if (log(i)) then
            do j=1,nane
              if (.not.log(j)) then
                mat(i,j)=mat(i,j)-imin
                mat(j,i)=mat(j,i)+imin
              endif
            enddo
          endif
        enddo
        do i=1,nane
          if (log(i)) then
            do j=1,nane
              if (.not.log(j)) then
                if (mat(i,j).eq.0) then
                  log(j)=.true.
                  pot(j)=pot(i)+imin
                endif
              endif
            enddo
          endif
        enddo
c       log(i2)=.true.
c       pot(i2)=pot(i1)+imin
      enddo
 555  return
      end
