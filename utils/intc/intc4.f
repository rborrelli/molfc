      Subroutine terminal (natom,nterm,iterm,nbond,ibond)
cversion 5.0
      parameter(mxneig=6)
      implicit real*8 (a-h,o-z)
C     this subroutine prepares the list of terminal atoms
c     input:  natom=number of atoms
c            nbond(j)= the number of bonds atom j is on
c            ibond(l,j),l=1,nbond(j)= the atoms j is bound to
c     output: nterm= number of terminal atoms, iterm (their numbers)
      dimension iterm(*),ibond(mxneig,natom),nbond(natom)
       nterm=0
       do 400 i=1,natom
          if (nbond(i).eq.1) then
             nterm=nterm+1
             iterm(nterm)=i
           end if
 400   continue
C      write (*,*) ' terminal atoms'
C      write (*,500) (iterm(j),j=1,nterm)
 500   format (5(2x,i4),3x,5(2x,i4))
       return
       end
C
      Subroutine prune (natom,nbond,ibond,nterm,iterm)
      parameter(mxneig=6)
      implicit real*8 (a-h,o-z)
c     This routine eliminates terminal atoms from the bond list
c     Input: natom is the number of atoms,nbond(i) is the number
c           of bonds atom i has, ibond (1..nbond(i),i) are the atoms
c           i is bound to
c           nterm is the number of terminal atoms
c           iterm is their serial numbers
c     Output: nbond and ibond: updated lists
c             They retain the original numbering but the atoms
c             which are eliminated have no neighbors
c
       dimension nbond(natom),ibond(mxneig,natom),iterm(*)
        do 500 i=1,natom
           if (nbond(i).eq.0) go to 500
           do 100 j=1,nterm
             if (iterm(j).eq.i) then
                 nbond(i)=0
                 go to 500
              end if
 100       continue
           k1=0
           do 300 k=1,nbond(i)
              kk=ibond(k,i)
              do 200 j=1,nterm
                 if (iterm(j).eq.kk) go to 300
 200          continue
              k1=k1+1
              ibond(k1,i)=kk
 300       continue
           nbond(i)=k1
 500   continue
C     write(*,*) 'bonds for the atoms after pruning'
       do 250 j=1,natom
C        write (*,400) j,nbond(j),(ibond(k,j),k=1,nbond(j))
 250   continue
 400  format(1x,2i3,3x,6i4)
      end
c
       Subroutine chainrin (natom,nbond,ibond,indsub,nsub,icycl)
c NEW: modified to handle SUBSTRUCTURES
       parameter(mxneig=6)
      implicit real*8 (a-h,o-z)
       dimension nbond(natom),ibond(mxneig,natom)
       dimension indsub(natom),nsub(natom)
       parameter (maxat=200)
       dimension mv(maxat),noring(maxat)
c..FG  when this subr is called with nbond1 and ibond1, these
c      bond lists contain only those atoms that have
c     survived the iterative pruning of terminal atoms.
c        (More precisely: nbond(i)=0 for all atoms which have
c         been cut off)
c     Surviving atoms can be: a) part of a ring
c       b) part of a chain closed by rings on all ends.
c     In this routine, type b atoms will be eliminated,
c      leaving only ring-atoms (part of possibly more
c     rings simultaneously)
c..   removal of a type b atom cuts the mol. into two parts
c     Now try if the removal of an atom divides the molecule into
c     two parts
c      logic: pick atom i; start with another atom mstart and mark it
c      [put 1 into the marking vector  - mv(mstart)=1]; scan the
c      neigbors of mstart and mark them too; continue to mark all
c      atoms that have already marked neigbors; if after this
c      cyclic process there remain unmarked atoms, atom i cuts
c      the molecule into two.
c....
c
c     in case of isolated substructures:
c     see routine subsets where the following was established:
c      indsub - indices of the atoms in the subsystems all collected
c      in one vector, the first nsub(1) defining subset 1., nsub(2)
c      subset 2., etc.
c     find the true index of the current atom:

      idim=nsub(icycl)
      nprev=0
      if(icycl.gt.1) then
        do 10 ic=1,icycl-1
        nprev=nprev+nsub(ic)
  10    continue
      endif

      do 1000 iat=1,idim
      i=indsub(nprev+iat)
c  the redefined actual index of the atom picked
        noring(i)=0
c
C..FG    the original below was:   if (nbond(i).ne.2) go to 1000
c       leave out obviously non-chain atoms
c  FG   the original thus left out the ring atom to which the chain
c      is attached (lets call it ring-corner atom)
c       thus avoiding that this atom  - which also cuts
c      the mol. in two  - would be qualified as non-ring
C..FG   consider, however:e.g. triphenyl-methan,
C..FG  an atom with 3 neighbors(the center) is also non-ring.
c      thus, all atoms will be included in the test now
c      and the ring-corner atom will be identified below by the fact
c      that it has at least two marked neigbors
c
        if (nbond(i).lt.2) go to 1000
c        this skips the atoms not part of the chain-ring system
c
c FG        j=0
c FG 600    j=j+1
c FG        if (j.gt.natom) go to 1000
c FG        if (nbond(j).eq.0.or.j.eq.i) go to 600
c FGc       j is a suitable starting atom
c
c new:   the start atom of the marking process will be a neigbor
c    of  i; for atoms with 3 or more neigbors the process
c     must be repeated starting with all neigbors
c    to find out if in any case i  gets two marked neigbors
c
        do 650 jk=1,nbond(i)
          mstart=ibond(jk,i)
          do 700 k=1,natom
            mv(k)=0
 700      continue
          mv(mstart)=1
 750      continue
          iadd=0
          do 900 kat=1,idim
           k=indsub(nprev+kat)
c   redefined value of k, the actual index of the running atom
            if (mv(k).ne.1) then
               do 800 kk=1,nbond(k)
                  if (mv(ibond(kk,k)).eq.1.and.k.ne.i) then
                      mv(k)=1
c                     if an atom has a marked neigbor, mark it
                      iadd=1
                   end if
 800           continue
            end if
 900      continue
        if (iadd.ne.0) go to 750
c
        if(nbond(i).lt.3) goto 660
c
c  FG   for atoms with 3 or more neigbors: decide in triphenyl-methane
c       whether i is the central (chain) atom or the ring-corner atom:
c      skip the test below if atom  i  has at least
c     two marked neigbors(this must be checked from all sides of i
c     that's why the marking process is now done repeatedly)
c     - if in any case two marked neigbors are found i  is a ring atom
c     even if it cuts the mol. in two
c     count the number of marked neigbors of i:
c
          markednb=0
          do 910 k=1,nbond(i)
            neighb=ibond(k,i)
            if(mv(neighb).eq.1) markednb=markednb+1
 910      continue
c
        if(markednb.ge.2) goto 1000
 650    continue
 660    continue
c
c     Now establish if the skeleton is still continous
c     if unmarked atoms are found [mv(k).eq.0]  atom i  cuts
c     the molecule in two
        do 950 kat=1,idim
           k=indsub(nprev+kat)
c   redefined value of k, the actual index of the running atom
           if (nbond(k).lt.2) go to 950
           if (mv(k).eq.0.and.k.ne.i) then
              noring(i)=1
c FG          go to 950
              goto 1000
           end if
 950    continue
 1000 continue
c     remove these atoms from the neighbor list
      do 1500 iat=1,idim
      i=indsub(nprev+iat)
        k1=0
        do 1400 k=1,nbond(i)
           neigh=ibond(k,i)
           if (noring(neigh).ne.1) then
              k1=k1+1
              ibond(k1,i)=neigh
            end if
 1400   continue
        nbond(i)=k1
 1500 continue
      do 1600 iat=1,idim
      i=indsub(nprev+iat)
         if (noring(i).eq.1) nbond(i)=0
 1600 continue
C     write(*,*) 'bonds for the atoms after chainrin'
C     write(*,*) 'i.e. after removing non-ring atoms'
       do 250 j=1,natom
C        write (*,300) j,nbond(j),(ibond(k,j),k=1,nbond(j))
 250   continue
 300  format(1x,2i3,3x,6i4)
      end
c
c
      Subroutine ringcoor (n,itemp,nbond,ibond,ipun,iptyp,ncg,
     1natom,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
      implicit real*8 (a-h,o-z)
      character*4 aromat
      character*2 isymb
c
c     This subroutine builds the ring coordinates
c     It tries to orient them in a sensible way
c     Input: n=ring size
c            itemp(1..n)= atoms constituting the ring
c            nbond(i)=the number of atoms i is bound to
c            ibond(1..nbond(i),i)=the atoms i is bound to
c            ipun=punch file
c
c..FG        added for use in subr fring:
c            natom - no of atoms, nb - no of bonds (in the whole molec)
c            arrays li,lj,rr coming from subr bonds
c            fcb - force constants from subr fbonds
c..
      parameter(mxneig=6)
c     mxneig=max no of neighbours for an atom
      dimension itemp(n),nbond(natom),ibond(mxneig,natom)
      dimension li(nb),lj(nb),rr(nb),fcb(nb)
      dimension fcq(*), scale(*)
      dimension isymb(natom)
c     Try first if there is a special atom
C FG  : reordering was wrong;
C     i is the relative index in the ring, while itemp(i) the abs.index
C     change i to itemp(i) every time when calling cartes
C
      zero=0.0d0
      aver=zero
      do 150 i=1,n
         call cartes (.false.,itemp(i),xx,yy,zz,iz)
         aver=aver+iz
 150  continue
      aver=aver/n
      atmax=zero
      iatmax=1
      do 210 i=1,n
         call cartes (.false.,itemp(i),xx,yy,zz,iz)
         del=abs(iz-aver)
         if (del.gt.atmax) then
            atmax=del
            iatmax=i
         end if
 210  continue
c   In case this is not successful, try to locate the atom with the
c   heaviest substituent(s)
      haver=zero
      do 260 i=1,n
        ib=nbond(itemp(i))
        do 240 j=1,ib
          call cartes (.false.,ibond(ib,itemp(i)),xx,yy,zz,iz)
          haver=haver+iz
 240    continue
 260  continue
      haver=haver/n
      havmax=0
      ihavmax=1
      do 280 i=1,n
        ib=nbond(itemp(i))
         havi=zero
        do 270 j=1,ib
          call cartes (.false.,ibond(ib,itemp(i)),xx,yy,zz,iz)
          havi=havi+iz
 270    continue
        del=abs(havi-haver)
        if (del.gt.havmax) then
           havmax=del
           ihavmax=i
         end if
 280  continue
      if (atmax.gt.0.0001) then
         i1=iatmax-1
         do 400 i=1,i1
            kk=itemp(1)
            do 300 j=1,n-1
              itemp(j)=itemp(j+1)
 300        continue
            itemp(n)=kk
 400     continue
      else if (havmax.gt.0.0001) then
         i1=ihavmax-1
         do 600 i=1,i1
            kk=itemp(1)
            do 500 j=1,n-1
              itemp(j)=itemp(j+1)
 500        continue
            itemp(n)=kk
 600     continue
      end if
c      We assume now that the ring begins with the highest substituent
C     write(*,700) (itemp(k),k=1,n)
 700   format (' reordered ring ',(10i4))
CHANGE.... the original program (subr) is substituted from here
c     by calling a subroutine ringdef, which works with general
c      formulae for any rings, rather than using individual
c      formulae for each ring
c....
      call ringdef(n,itemp,ipun,iptyp,ncg,nbond,ibond,natom,nb,ian,
     1  li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
      return
      end
c
      Subroutine ringdef(n,itemp,ipun,iptyp,ncg,nbond,ibond,natom,
     1nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
c..
c     generation of bend and torsional coordinates for larger rings
c     redundancies removed as commented below
c     n = dimension of the ring; itemp(n) - indices of the
c     ring atoms
c     ipun is the punch file, onto which the coordnts are written
c     ncg counts the no of coordts generated
c..
      implicit real*8 (a-h,o-z)
      character*4 aromat
      character*2 isymb
      parameter(mxdr=30,mxneig=6)
c     mxdr=max dim of a ring, mxneig=max no of neighbors for an atom
      dimension itemp(n)
c     the following needed for fring called in this routine
c     natom=no of atoms, nb=no of bonds, in the whole molec
c     arrays li,lj,rr coming from subr bonds
c     fcb force constants for bonds only
c     fcq contains the force constants for the entire set of coordts
      dimension nbond(natom),ibond(mxneig,natom)
      dimension li(nb),lj(nb),rr(nb),fcb(nb)
      dimension fcq(*), scale(*)
      dimension isymb(natom)
      dimension c(mxdr,mxdr), ctau(mxdr)
chr+
cd     write(*,*) 'Start of ringdef'
chr-
      one=1.0d0
      two=2.0d0
      delta=1.0d-5
      pi=acos(-one)
      pi2=two*pi
      n2=n/2
      mm=0
c          old:    do 100 m=0,n2
c.. the totally symm species (m=0) and the first 2-dim species (m=1)
c   are the redundancies, ALWAYS. Therefore, start with m=2
      do 100 m=2,n2
      mm=mm+1
      do 200 k=1,n
      c(mm,k)=cos(m*(k-1)*pi2/n)
 200  continue
c            old:    if(m.eq.0) goto 100
      if(m.eq.n2.and.mod(n,2).eq.0) goto 100
c     even membered ring
      mm=mm+1
      do 300 k=1,n
      c(mm,k)=sin(m*(k-1)*pi2/n)
 300  continue
 100  continue
      do 500 i=1,n-3
c..   be careful: if a coefficient is zero, bmat program interprets
c       it as 1
      kfirst=0
      do 550 k=1,n
      if (abs(c(i,k)).gt.delta) then
        kfirst=k
        goto 555
      endif
550   continue
555   k1=kfirst
      if(k1.eq.1) then
        km=n
      else
        km=k1-1
      endif
      if(k1.lt.n) then
        kp=k1+1
      else
        kp=k1+1-n
      endif
      ncg=ncg+1
      write(ipun,1000) ncg,c(i,kfirst),itemp(km),itemp(kp),itemp(kfirst)
      ityp=2
      if(i.eq.(n-3).and.mod(n,2).eq.0) ityp=1
c       2- and 1-dim species will be treated separately by fring
c       the only 1-dim sp. is the last one, in an even memebered
c       ring only
      call fring(n,'BEND',ityp,natom,nb,itemp,ian,li,lj,rr,fcb,
     1  fcq(ncg),scale(ncg),aromat)
      write(iptyp,1070) ncg, n,scale(ncg)
 1070 format(i4,'. ',i2,'-membered ring ',' BEND  ','sc=',f5.4,2x)
      do 500 k=kfirst+1,n
      km=k-1
      if(k.lt.n) then
        kp=k+1
      else
        kp=k+1-n
      endif
      if (abs(c(i,k)).gt.delta) then
        write(ipun,2000) c(i,k),itemp(km),itemp(kp),itemp(k)
      endif
  500 continue
c...   torsions follow
      if (mod(n,2).ne.0) then
c     odd-membered ring
c     torsions oriented such that first tau is opposite to alfa1(bend)
        do 600 i=1,n-3
c..   be careful: if a coefficient is zero, bmat program interprets
c       it as 1
          kfirst=0
          do 650 k=1,n
          if (abs(c(i,k)).gt.delta) then
            kfirst=k
            goto 655
          endif
650       continue
655       continue
          k1=kfirst-1+n2
           if(k1.gt.n) k1=k1-n
          k2=k1+1
           if(k2.gt.n) k2=k2-n
          k3=k1+2
           if(k3.gt.n) k3=k3-n
          k4=k1+3
           if(k4.gt.n) k4=k4-n
          ncg=ncg+1
          write(ipun,8000) ncg,c(i,kfirst),itemp(k1),itemp(k2),
     1      itemp(k3), itemp(k4)
c       there are only 2-dim species in odd-membered rings:
          call fring(n,'TORS',2,natom,nb,itemp,ian,li,lj,rr,fcb,
     1     fcq(ncg),scale(ncg),aromat)
          write(iptyp,1080) ncg, n,scale(ncg)
 1080     format(i4,'. ',i2,'-membered ring ',' TORS  ','sc=',f5.4,2x)
            do 600 k=kfirst+1,n
              k1=k1+1
               if(k1.gt.n) k1=k1-n
              k2=k2+1
               if(k2.gt.n) k2=k2-n
              k3=k3+1
               if(k3.gt.n) k3=k3-n
              k4=k4+1
               if(k4.gt.n) k4=k4-n
           if(abs(c(i,k)).gt.delta) then
            write(ipun,9000) c(i,k),itemp(k1),itemp(k2),itemp(k3),
     1      itemp(k4)
           endif
  600   continue
      else
c     in even-membered rings start with the antisymm. combintn of the
c     two taus opposite to alfa1
c     e.g. in a six-ring, start with tau4-tau3
c...
c     the logic will be different from that above for odd-rings
c     because each tau appears twice, its coefficients will be
c     collected and each tau punched only once; otherwise strange
c     renormalization  would also be necessary
c...
        do 700 i=1,n-3
          k2=n2
          do 710 k=1,n
c       k2 is the index of a given tau;e.g. tau1234=tau2
          k2=k2+1
          if(k2.gt.n) k2=k2-n
          kp=k+1
          if(kp.gt.n) kp=kp-n
          ctau(k2) = c(i,k)-c(i,kp)
  710     continue
c       ctau now contains, in this order, the coeffts of tau1, tau2,..
c
c..   be careful: if a coefficient is zero, bmat program interprets
c       it as 1, start with the first non-zero
c
          kfirst=0
          do 750 k=1,n
          if (abs(ctau(k)).gt.delta) then
            kfirst=k
            goto 755
          endif
750       continue
755       continue
          k1=kfirst-1
           if(k1.eq.0) k1=n
          k2=k1+1
           if(k2.gt.n) k2=k2-n
          k3=k1+2
           if(k3.gt.n) k3=k3-n
          k4=k1+3
           if(k4.gt.n) k4=k4-n
          ncg=ncg+1
          write(ipun,8000) ncg,ctau(k2),itemp(k1),itemp(k2),itemp(k3),
     1    itemp(k4)
          ityp=2
          if(i.eq.(n-3).and.mod(n,2).eq.0) ityp=1
          call fring(n,'TORS',ityp,natom,nb,itemp,ian,li,lj,rr,fcb,
     1     fcq(ncg),scale(ncg),aromat)
          write(iptyp,1080) ncg, n,scale(ncg)
        do 700 k=kfirst+1,n
          k1=k1+1
           if(k1.gt.n) k1=k1-n
          k2=k2+1
           if(k2.gt.n) k2=k2-n
          k3=k3+1
           if(k3.gt.n) k3=k3-n
          k4=k4+1
           if(k4.gt.n) k4=k4-n
         if(abs(ctau(k2)).gt.delta) then
          write(ipun,9000) ctau(k2),itemp(k1),itemp(k2),itemp(k3),itemp
     1    (k4)
         endif
  700   continue
      endif
 1000    format(i4,2x,'K',3x,f10.7,2x,'BEND',6x,3(i4,'.',5x))
 2000    format(' ',9x,f10.7,2x,'BEND',6x,3(i4,'.',5x))
 8000    format(i4,2x,'K',3x,f10.7,2x,'TORS',6x,4(i4,'.',5x))
 9000    format(' ',9x,f10.7,2x,'TORS',6x,4(i4,'.',5x))
chr+
cd     write(*,*) 'End of ringdef'
chr-
      return
      end
      Subroutine centcoor (i,nterm,neigterm,nnon,neignon,
     1ipun,iptyp,natom,ncg,isymb,ibord,scale,dummy)
      implicit real*8 (a-h,o-z)
      character*2 isymb
      logical dummy
      parameter (mxneig=6)
      dimension neigterm(*),neignon(*)
      dimension isymb(*), ibord(*), scale(*)
      dimension it(mxneig),iz(mxneig)
      parameter (planar=5.6d0,pi=3.141)
      parameter(scal1=0.8,scal2=0.8,scal3=0.8,scal4=0.8,scal5=0.8)
c..FG this routine treats a central atom (with index i) with terminal
c      neighbours only; nnon and neignon, referring to non-terminals
c..    are in fact redundant, not used
c
      if(nterm.gt.6) then
        print *, 'more than 6 ligands in centcoor not implemented'
        call bummer('more than 6 ligands in centcoor',0,2)
        STOP
c
      elseif(nterm.eq.5.or.nterm.eq.6) then
c      penta-hexavalent case added, these are
c      treated by an independnt routine
c
       call pentahex(i,nterm,nnon,neigterm,neignon,ipun,iptyp,
     1    ncg,isymb,scale)
      endif
c
c     now the cases with 4 or less ligands from the original:
c...  Try to locate a unique atom
      aver=0.0d0
      do 50 kk=1,nterm
        call cartes (.false.,neigterm(kk),xx,yy,zz,iz(kk))
        aver=aver+iz(kk)
 50   continue
      aver=aver/nterm
      devi=0.0d0
      maxi=1
      do 60 kk=1,nterm
        if (abs(iz(kk)-aver).gt.devi ) then
          devi=abs(iz(kk)-aver)
          maxi=kk
         end if
 60   continue

      if (nterm.eq.4.and.devi.lt.1.0d-9) then
c   XY4 type
 100     format(i4,2x,'K',3x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
 200     format(' ',9x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
         ncg=ncg+1
         write(ipun,100) ncg,2, neigterm(1),neigterm(2),i
         scale(ncg)=scal1
         write(iptyp,410) ncg,i,scale(ncg)
  410    format(i4,'. ','XY4, ctr=',i3,' Ea def  ','sc=',f5.4,2x)
         write(ipun,200) 2,neigterm(3),neigterm(4),i
         write(ipun,200) -1,neigterm(1),neigterm(3),i
         write(ipun,200) -1,neigterm(1),neigterm(4),i
         write(ipun,200) -1,neigterm(2),neigterm(3),i
         write(ipun,200) -1,neigterm(2),neigterm(4),i
c    E def. a
c
         ncg=ncg+1
         write(ipun,100)  ncg,1,neigterm(1),neigterm(3),i
         scale(ncg)=scal1
         write(iptyp,412) ncg,i,scale(ncg)
 412     format(i4,'. ','XY4, ctr=',i3,' Eb def  ','sc=',f5.4,2x)
         write(ipun,200) -1,neigterm(1),neigterm(4),i
         write(ipun,200) -1,neigterm(2),neigterm(3),i
         write(ipun,200)  1,neigterm(2),neigterm(4),i
c     E def b
c
         ncg=ncg+1
         write(ipun,100) ncg,1,neigterm(1),neigterm(2),i
         scale(ncg)=scal1
         write(iptyp,414) ncg,i,scale(ncg)
 414     format(i4,'. ','XY4, ctr=',i3,' F2a def ','sc=',f5.4,2x)
         write(ipun,200) -1,neigterm(3),neigterm(4),i
c    F2 1
c
         ncg=ncg+1
         write(ipun,100) ncg,1,neigterm(1),neigterm(3),i
         scale(ncg)=scal1
         write(iptyp,416) ncg,i,scale(ncg)
 416     format(i4,'. ','XY4, ctr=',i3,' F2b def ','sc=',f5.4,2x)
         write(ipun,200) -1,neigterm(2),neigterm(4),i
c    F2 2
c
         ncg=ncg+1
         write(ipun,100) ncg,1,neigterm(1),neigterm(4),i
         scale(ncg)=scal1
         write(iptyp,418) ncg,i,scale(ncg)
 418     format(i4,'. ','XY4, ctr=',i3,' F2b def ','sc=',f5.4,2x)
         write(ipun,200) -1,neigterm(3),neigterm(2),i
c    F2 3
c
       else if (nterm.eq.4.and. devi.ge.1.0d-9) then
c   XYZ3 type
       it(1)=neigterm(maxi)
       l1=1
       do 90 ll=1,4
         if (ll.eq.maxi) go to 90
         l1=l1+1
         it(l1)=neigterm(ll)
 90    continue
          ncg=ncg+1
          write(ipun,100) ncg,1,it(3),it(4),i
         scale(ncg)=scal2
         write(iptyp,510) ncg,i,scale(ncg)
 510     format(i4,'. ','XYZ3,ctr=',i3,' sym def ','sc=',f5.4,2x)
          write(ipun,200) 1,it(2),it(4),i
          write(ipun,200) 1,it(2),it(3),i
          write(ipun,200) -1,it(1),it(2),i
          write(ipun,200) -1,it(1),it(3),i
          write(ipun,200) -1,it(1),it(4),i
c  sym. def.
c
          ncg=ncg+1
          write(ipun,100) ncg,1,it(3),it(4),i
         scale(ncg)=scal2
         write(iptyp,512) ncg,i,scale(ncg)
 512     format(i4,'. ','XYZ3,ctr=',i3,' asdef,a ','sc=',f5.4,2x)
          write(ipun,200) -1,it(2),it(4),i
          write(ipun,200) -1,it(2),it(3),i
c  as. def a
c
          ncg=ncg+1
          write(ipun,100)  ncg,1,it(2),it(4),i
         scale(ncg)=scal2
         write(iptyp,514) ncg,i,scale(ncg)
 514     format(i4,'. ','XYZ3,ctr=',i3,' asdef,b ','sc=',f5.4,2x)
          write(ipun,200) -1,it(2),it(3),i
c  as def b
c
          ncg=ncg+1
          write(ipun,100)  ncg,2,it(1),it(2),i
         scale(ncg)=scal2
         write(iptyp,516) ncg,i,scale(ncg)
 516     format(i4,'. ','XYZ3,ctr=',i3,' rock ,a ','sc=',f5.4,2x)
          write(ipun,200) -1,it(1),it(3),i
          write(ipun,200) -1,it(1),it(4),i
c  rocking a
c
          ncg=ncg+1
          write(ipun,100)  ncg,1,it(1),it(3),i
         scale(ncg)=scal2
         write(iptyp,518) ncg,i,scale(ncg)
 518     format(i4,'. ','XYZ3,ctr=',i3,' rock ,b ','sc=',f5.4,2x)
          write(ipun,200) -1,it(1),it(4),i
c  rocking b
c
        else if (nterm.eq.3) then
        call angle (neigterm(1),i,neigterm(2),al1)
        call angle (neigterm(1),i,neigterm(3),al2)
        call angle (neigterm(2),i,neigterm(3),al3)
        sumang=al1+al2+al3
        it(1)=neigterm(maxi)
c      re-order the ligands - the first one is the unique one
        k1=1
        do 150 kk=1,3
          if (kk.eq.maxi) go to 150
          k1=k1+1
          it(k1)=neigterm(kk)
 150    continue
        if (sumang.lt.planar) then
c         sym. def.  non-planar case
          ncg=ncg+1
          write (ipun,100) ncg,1,it(1),it(2),i
          scale(ncg)=scal3
          write(iptyp,520) ncg,i,scale(ncg)
 520      format(i4,'. ','XY3 ,ctr=',i3,' sym def ','sc=',f5.4,2x)
          write(ipun,200)  1,it(1),it(3),i
          write(ipun,200)  1,it(2),it(3),i
        else
          ncg=ncg+1
          write (ipun,300) ncg,1,it(1),it(2),it(3),i
          scale(ncg)=scal4
          write(iptyp,522) ncg,i,scale(ncg)
 522      format(i4,'. ','XY3 ,ctr=',i3,' sym out ','sc=',f5.4,2x)
          write (ipun,310) 1,it(2),it(3),it(1),i
          write (ipun,310) 1,it(3),it(1),it(2),i
 300      format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
 310      format(' ',9x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
        end if
c
        ncg=ncg+1
        write (ipun,100) ncg,2,it(2),it(3),i
        scale(ncg)=scal3
        write(iptyp,524) ncg,i,scale(ncg)
 524    format(i4,'. ','XY3 ,ctr=',i3,' asdef,a ','sc=',f5.4,2x)
        write(ipun,200) -1,it(1),it(3),i
        write(ipun,200) -1,it(1),it(2),i
c  asym. def. a
c
        ncg=ncg+1
        write(ipun,100)  ncg,1,it(1),it(3),i
        scale(ncg)=scal3
        write(iptyp,526) ncg,i,scale(ncg)
 526    format(i4,'. ','XY3 ,ctr=',i3,' asdef,b ','sc=',f5.4,2x)
        write(ipun,200) -1,it(1),it(2),i
c  asym. def. b
c
      else if (nterm.eq.2) then
c FG      test for linearity
        call angle(neigterm(1),i,neigterm(2),alpha)
        if(abs(alpha).lt.pi) then
          ncg=ncg+1
          write (ipun,100) ncg,1,neigterm(1),neigterm(2),i
          scale(ncg)=scal5
          write(iptyp,530) ncg,i,scale(ncg)
 530      format(i4,'. ','XY2 ,ctr=',i3,' bending ','sc=',f5.4,2x)
        else
c    linear triatomic molecule
          if(dummy) then
            write(ipun,810)
 810        format('?note:existing dummy used to define LIN1,LIN2')
          else
            write(ipun,820)
 820        format('?CAUTION,nonexisting dummy introduced!')
          endif
          ifour=natom+1
          it1=neigterm(1)
          it2=neigterm(2)
          ncg=ncg+1
          write (ipun,700) ncg,1,it1,it2,i,ifour
           write(iptyp,630) ncg,isymb(it1),isymb(i),isymb(it2),
     1     scale(ncg)
 630       format(i4,'. ','lin triat  ',3a2,' LIN1  ','sc=',f5.4,2x)
          ncg=ncg+1
          write (ipun,800) ncg,1,it1,it2,i,ifour
           write(iptyp,632) ncg,isymb(it1),isymb(i),isymb(it2),
     1     scale(ncg)
 632       format(i4,'. ','lin triat  ',3a2,' LIN2  ','sc=',f5.4,2x)
 700     format(i4,2x,'K',3x,i4,'.',5x,'LIN1',6x,4(i4,'.',5x))
 800     format(i4,2x,'K',3x,i4,'.',5x,'LIN2',6x,4(i4,'.',5x))
        endif
      end if
      RETURN
c end of Centcoor
      end
c
      Subroutine primcoor (i,nterm,neigterm,nnon,neignon,ipun,iptyp,
     1  nbond,ibond,natom,ncg,isymb,ibord,scale,dummy)
      implicit real*8 (a-h,o-z)
      character*2 isymb
      parameter (obtuse=2.8d0)
      parameter(mxneig=6,atomdiff=0.2,zr=1.d-4)
      logical found, dummy
      dimension neigterm(nterm),neignon(nnon),iat(3)
      dimension nbond(natom),ibond(mxneig,natom)
      dimension isymb(natom), ibord(*), scale(*)
c
c     This routine generates internal coordinates for a primary
c     center
c     Parameters: i is the primary atom, nterm is the number
c                 of its terminal neighbors, neigterm(1..nterm)
c                 are the terminal neigbors
c                 nnon is the number of non-terminal neigbors (this
c                 must be 1 for a primary atom), and nbond(1) is
c                 the non-terminal neighbor
c                 nbond(i) is the number of bonds for atom i,
c                 ibond (k,i) is the k-th atom bound to i
c                 (k must be less than 7)
c                 natom is the total number of atoms
c
      ibd=neignon(1)
c     this is the connection to the frame of the molecule
c      This is about 160 degrees and is the linit for approximate
c       collinearity of two bonds
c  FG
      atdiff=atomdiff
      ipln=0
c  FG
      if(nterm.gt.5) then
        print *, 'stop in primcoor'
        print *, 'higher than hexavalent not implemented'
        call bummer('higher than hexavalent not implemented',0,2)
        STOP
c
      elseif(nterm.eq.4.or.nterm.eq.5) then
        call pentahex(i,nterm,nnon,neigterm,neignon,ipun,
     1    iptyp,ncg,isymb,scale)
c
      elseif (nterm.eq.3) then
c     Try later to adapt the coordinates to global symmetry
c FG  (see a special  test added below)
         it1=neigterm(1)
         it2=neigterm(2)
         it3=neigterm(3)
c     determine the unique terminal atom
         call cartes (.false.,it1,xx,yy,zz,iat(1))
c FG   preparing a special new test for unique atom further below,
c      check if one of the coordinates is zero
          if(abs(xx).lt.zr.or.abs(yy).lt.zr.or.abs(zz).lt.zr) ipln=1
         call cartes (.false.,it2,xx,yy,zz,iat(2))
          if(abs(xx).lt.zr.or.abs(yy).lt.zr.or.abs(zz).lt.zr) ipln=2
         call cartes (.false.,it3,xx,yy,zz,iat(3))
          if(abs(xx).lt.zr.or.abs(yy).lt.zr.or.abs(zz).lt.zr) ipln=3
         aver=(iat(1)+iat(2)+iat(3))/3.0d0
         nmax=0
         do  50 k=1,3
c  FG
         diff=(abs(aver-iat(k)))/aver
         if(diff.gt.atdiff) then
           atdiff=diff
           nmax=k
         end if
 50      continue
           if (nmax.eq.0) then
c....      there is no unique atom in the XY3 group
c  FG
c        the original scheme below  - based on maximum coplanarity -
c        is not always appropriate (e.g., good for one conformer
c        of toluene, not for the other;methylamine also wrong)
c        as a temporary solution, the atom which lies in one
c        of the cartesian planes is taken as unique  - works
c        of course only with right selection of the XYZ system
c        ipln from above specifies that Y which is in  a plane
c        the original test below will be skipped if there is such
c        an Y .... reconsider it later, may not be ideal
c
             if(ipln.gt.0) then
               nmax=ipln
               goto 82
             endif
c FG end, now comes the original test IF there was no Y in plane
c
c....      try to locate a unique atom by the connection X3Y-ZU(n)
c
           kbond=nbond(ibd)
c
c....      This is the number of bonds of ibd
c....      Calculate the dihedral angles. Select the one with the
c....      largest cosine abs. value  (this is most co-planar)
c
           dih=0.0D0
c never used:           linunit=0
           do 70 itm=1,3
c  remove if ok without  this    if(linunit.eq.1) goto 84
             ifour=0
             itt=neigterm(itm)
             do 80 noth=1,kbond
                no11=ibond(noth,ibd)
                if(no11.eq.i) go to 80
                no1=no11
c..
c..FG    take ch3cn as example; when calling dihed below,
c        itt is a hydrogen, i=C, ibd=C, no1=N; the latter 3 should
c        define a plane, i.e. cannot be collinear
c        check for collinearity
                pi=acos(-1.d0)
                call angle(i,ibd,no1,alfa)
                if(abs(alfa-pi).gt.0.0005) then
                  call dihed (itt,i,ibd,no1,theta)
c           theta is the cosine of the dihedral angle. Select the
c           largest value
c
                  d1=abs(theta)
                  if (d1.gt.dih) then
                     dih=d1
                     nmax=itm
                   end if
                endif
 80           continue
c
c  if nmax  stayed zero at this place, it means that no non-linear
c  i-ibd-no1 unit was found, ( considering all possible no1 neighbours
c  of ibd, although if this unit is linear, there can hardly be more
c  than one no1). In case of ch3ccn this means that none of the
c  hydrogens is a unique atom. Think of ch3-cc-phenyl, however:
c  (or: ch3-c=c-nh, h out of line)
c  check for atoms farther away which may distinguish one of the
c  ch3 hydrogens
c
             if(nmax.eq.0) then
c never used:            linunit=1
  84        call ncollin(natom,nbond,ibond,no1,ibd,i,ifour,found,it1,
     1        it2,it3)
                   if(ifour.gt.0) then
c     ifour takes now the role of no1, to define a plane with i and ibd
                     call dihed(itt,i,ibd,ifour,theta)
                   else
c  ..behind the xy3 group the whole mol is linear, no unique y atom
                     goto 82
                   endif
               endif
 70        continue
 82      continue
         endif
c....  nmax (1 to 3) is the unique atom
             if (nmax.eq.2) then
                it1=neigterm(2)
                it2=neigterm(1)
             else if (nmax.eq.3) then
                it1=neigterm(3)
                it3=neigterm(1)
             end if
c
         ncg=ncg+1
         write(ipun,100) ncg,1,it2,it3,i
         write (ipun,200) 1,it1,it3,i
         write (ipun,200) 1,it1,it2,i
         write (ipun,200) -1,ibd,it1,i
         write (ipun,200) -1,ibd,it2,i
         write (ipun,200) -1,ibd,it3,i
c        scale(ncg)=.....
         write(iptyp,520) ncg,isymb(i),i,scale(ncg)
 520     format(i4,'. ','prim XY3 , X=',a2,i2,' SDEF  ','sc=',f5.4,2x)
 100     format(i4,2x,'K',3x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
 200     format(' ',9x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
c      symmetric def coordinate
c
         ncg=ncg+1
         write(ipun,100) ncg,2,it2,it3,i
         write (ipun,200) -1,it1,it3,i
         write (ipun,200) -1,it1,it2,i
c          scale(ncg)=....
         write(iptyp,522) ncg,isymb(i),i,scale(ncg)
 522     format(i4,'. ','prim XY3 , X=',a2,i2,' ADEFa ','sc=',f5.4,2x)
c      asymm. def a'
c
         ncg=ncg+1
         write (ipun,100) ncg,1,it1,it3,i
c          scale(ncg)=....
         write (ipun,200) -1,it1,it2,i
         write(iptyp,524) ncg,isymb(i),i,scale(ncg)
 524     format(i4,'. ','prim XY3 , X=',a2,i2,' ADEFb ','sc=',f5.4,2x)
c      asymm. def a"
c
         ncg=ncg+1
         write (ipun,100)  ncg,2,ibd,it1,i
c          scale(ncg)=....
         write(iptyp,526) ncg,isymb(i),i,scale(ncg)
 526     format(i4,'. ','prim XY3 , X=',a2,i2,' ROCKa ','sc=',f5.4,2x)
         write (ipun,200) -1,ibd,it2,i
         write (ipun,200) -1,ibd,it3,i
c      rock a'
c
         ncg=ncg+1
         write (ipun,100)  ncg,1,ibd,it2,i
         write (ipun,200) -1,ibd,it3,i
c          scale(ncg)=....
         write(iptyp,528) ncg,isymb(i),i,scale(ncg)
 528     format(i4,'. ','prim XY3 , X=',a2,i2,' ROCKb ','sc=',f5.4,2x)
c      rock a"
      else if  (nterm.eq.2) then
         it1=neigterm(1)
         it2=neigterm(2)
c
c FG  in the original, there was no search for a unique atom here
c  and just it1 and it2 were used to define with i the first angle(alfa)
c  this is OK, e.g. in ethylene, but think of formamide:
c  the OCN angle should be used, rather than the HCO
c
         call cartes (.false.,it1,xx,yy,zz,iat(1))
         call cartes (.false.,it2,xx,yy,zz,iat(2))
         aver=(iat(1)+iat(2))/2.d0
         diff=abs((iat(1)-iat(2))/aver)
         if(diff.lt.atomdiff) then
           inb1=it1
           inb2=it2
           inb3=ibd
         else
           inb1=ibd
           if(iat(1).gt.iat(2)) then
             inb2=it1
             inb3=it2
           else
             inb2=it2
             inb3=it1
           endif
         endif
c
         ncg=ncg+1
         write (ipun,100) ncg,2, inb1,inb2,i
        write(iptyp,620) ncg,isymb(i),i,scale(ncg)
 620    format(i4,'. ','pXY2,2-1-1,X=',a2,i2,' SCIS  ','sc=',f5.4,2x)
         write (ipun,200) -1,inb3,inb1,i
         write (ipun,200) -1,inb3,inb2,i
c      symm. def
c
         ncg=ncg+1
         write (ipun,100) ncg,1,inb3,inb1,i
        write(iptyp,622) ncg,isymb(i),i,scale(ncg)
 622    format(i4,'. ','pXY2, 1-1, X=',a2,i2,' ROCK  ','sc=',f5.4,2x)
         write (ipun,200) -1,inb3,inb2,i
c      rock
c
         ncg=ncg+1
         write (ipun,300) ncg,1,inb3,inb1,inb2,i
        write(iptyp,624) ncg,isymb(i),i,scale(ncg)
 624    format(i4,'. ','pXY2, oop, X=',a2,i2,' WAGG  ','sc=',f5.4,2x)
 300     format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
c      wagging
       else if (nterm.eq.1) then
         it1=neigterm(1)
c
         call angle (ibd,i,it1,alpha)
         if (alpha.lt.obtuse)  then
           ncg=ncg+1
           write (ipun,100) ncg,1,ibd,it1,i
c        scale(ncg)=...
         write(iptyp,626) ncg,isymb(ibd),ibd,isymb(it1),it1,
     1   isymb(i),i,scale(ncg)
 626     format(i4,'. ','pXY, ',3(a2,i2),' BEND  ','sc=',f5.4,2x)
         else
c          try to find a non-collinear fourth atom bound to ibd
           call ncollin (natom,nbond,ibond,ibd,i,it1,ifour,
     1     found,-1,-1,-1)
c          This routine tries to locate an atom (if possible bound to
c            ibd) which is not collinear with the i-ibd bond
c             the result is ifour
           if(found) then
             goto 90
           else
             if(dummy) then
               write(ipun,810)
 810           format('?note:existing dummy used to define LIN1,LIN2')
             else
               write(ipun,820)
 820           format('?CAUTION,nonexisting dummy introduced!')
             endif
           endif
  90       continue
           ncg=ncg+1
           write (ipun,700) ncg,1,ibd,it1,i,ifour
           write(iptyp,630) ncg,isymb(ibd),ibd,isymb(it1),it1,
     1     isymb(i),i,scale(ncg)
 630       format(i4,'. ','pXY, ',3(a2,i2),' LIN1  ','sc=',f5.4,2x)
           ncg=ncg+1
           write (ipun,800) ncg,1,ibd,it1,i,ifour
           write(iptyp,632) ncg,isymb(ibd),ibd,isymb(it1),it1,
     1     isymb(i),i,scale(ncg)
 632       format(i4,'. ','pXY, ',3(a2,i2),' LIN2  ','sc=',f5.4,2x)
 700     format(i4,2x,'K',3x,i4,'.',5x,'LIN1',6x,4(i4,'.',5x))
 800     format(i4,2x,'K',3x,i4,'.',5x,'LIN2',6x,4(i4,'.',5x))
c      deformation
         end if
       else
          write (*,*) 'impossible (in primcoor)'
       end if
      return
      end
c
      Subroutine secoor (i,nterm,neigterm,nnon,neignon,ipun,iptyp,
     1iring,ncg,natom,nbond,ibond,isymb,ibord,scale)
      implicit real*8 (a-h,o-z)
      character*2 isymb
      logical found
      parameter (obtuse=2.8d0)
      parameter (mxnr=40,mxdr=30)
      parameter(mxneig=6)
      dimension neigterm(*),neignon(*),iring(mxnr,*)
      dimension nbond(natom),ibond(mxneig,natom)
      dimension isymb(natom), ibord(*), scale(*)
c     This routine generates internal coordinates for a secondary
c     center
c FG  if i is a ring-atom, only the terminal atoms need be handled here
c     (ring coordinates were handled by ring)
c     for a non-ring atom, the skeletal BEND (or LIN) is generated
c
      if(nterm.gt.4) then
        print *, 'higher than hexavalent case not allowed in secoor'
        call bummer('higher than hexavalent case',0,2)
        STOP
      endif
c
      ib1=neignon(1)
      ib2=neignon(2)
c     these are the connections to the frame of the molecule
      isum=0
      do 50 j=1,mxnr
         isum=isum+iring(j,i)
 50   continue
      if (isum.eq.0) then
c   non-ring system:
c
c FG: the following jumpout in penta- hexavalent case was added later:
c      (designed for non-ring system only)
        if(nterm.gt.2) then
        call pentahex(i,nterm,nnon,neigterm,neignon,ipun,
     1    iptyp,ncg,isymb,scale)
        RETURN
        endif
c
c FG end,  now the original treatment of max tetravalent cases:
c   the skeletal BEND (or LIN) created here:
c
         call angle (ib1,i,ib2,beta)
         if (beta.lt.obtuse) then
           ncg=ncg+1
           write (ipun,100) ncg,1,ib1,ib2,i
c          scale(ncg)=...
           write(iptyp,510) ncg,isymb(ib1),ib1,isymb(ib2),ib2,
     1     isymb(i),i,scale(ncg)
 510     format(i4,'. ','skel,',3(a2,i2),' BEND  ','sc=',f5.4,2x)
c        Skeletal deformation in non-ring system
          else
c
            call ncollin (natom,nbond,ibond,ib1,i,ib2,ifour,found,
     1        -1,-1,-1)
             if (.not.found) then
               call ncollin (natom,nbond,ibond,ib2,i,ib1,ifour,found,
     1          -1,-1,-1)
             endif
 700     format(i4,2x,'K',3x,i4,'.',5x,'LIN1',6x,4(i4,'.',5x))
 800     format(i4,2x,'K',3x,i4,'.',5x,'LIN2',6x,4(i4,'.',5x))
             ncg=ncg+1
             write (ipun,700) ncg,1,ib1,ib2,i,ifour
             scale(ncg)=0.7
             write(iptyp,520) ncg,isymb(ib1),ib1,isymb(ib2),ib2,
     1       isymb(i),i,scale(ncg)
 520     format(i4,'. ','skel,',3(a2,i2),' LIN1  ','sc=',f5.4,2x)
             ncg=ncg+1
             write (ipun,800) ncg,1,ib1,ib2,i,ifour
             scale(ncg)=0.7
             write(iptyp,530) ncg,isymb(ib1),ib1,isymb(ib2),ib2,
     1       isymb(i),i,scale(ncg)
 530     format(i4,'. ','skel,',3(a2,i2),' LIN2  ','sc=',f5.4,2x)
          end if
      end if
c FG  end of non-ring skel. bend
c   from here, the terminal atoms are handled; whether i  is
c   a ring-atom or not, is indifferent for the original, max.
c   tetravalent case:
c
      if(nterm.gt.2.and.isum.gt.1) then
        print *,'a ring atom of valency higher than 4 found in secoor'
        print *, 'this is not yet not handled by pentahex'
        call bummer('ring atom of valency higher than 4 ',0,2)
        STOP
c note:in case of isum=0,i.e. non-ring, there was a jumpout above
c   for the penta- hexaval. case
      endif
c
      if (nterm.eq.2) then
         it1=neigterm(1)
         it2=neigterm(2)
c
         ncg=ncg+1
         write(ipun,100) ncg,4,it1,it2,i
c         scale(ncg)=...
          if (isum.eq.0) then
           write(iptyp,610) ncg, isymb(i),i,scale(ncg)
 610      format(i4,'. ','secd XY2 , X=',a2,i2,' SCIS  ','sc=',f5.4,2x)
          else
           write(iptyp,612) ncg, isymb(i),i,scale(ncg)
 612      format(i4,'. ','ring sXY2, X=',a2,i2,' SCIS  ',
     1    'sc=',f5.4,2x)
          endif
         write (ipun,200) -1,it1,ib1,i
         write (ipun,200) -1,it1,ib2,i
         write (ipun,200) -1,it2,ib1,i
         write (ipun,200) -1,it2,ib2,i
 100     format(i4,2x,'K',3x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
 200     format(' ',9x,i4,'.',5x,'BEND',6x,3(i4,'.',5x))
c      symm. def of the XY2 group
c
         ncg=ncg+1
         write (ipun,100)  ncg,1,it1,ib1,i
         write (ipun,200)  1,it1,ib2,i
         write (ipun,200) -1,it2,ib1,i
         write (ipun,200) -1,it2,ib2,i
c      XY2 rocking
c         scale(ncg)=....
          if (isum.eq.0) then
           write(iptyp,620) ncg, isymb(i),i,scale(ncg)
 620      format(i4,'. ','secd XY2 , X=',a2,i2,' ROCK  ','sc=',f5.4,2x)
          else
           write(iptyp,622) ncg, isymb(i),i,scale(ncg)
 622      format(i4,'. ','ring sXY2, X=',a2,i2,' ROCK  ',
     1    'sc=',f5.4,2x)
          endif
         ncg=ncg+1
         write (ipun,100)  ncg,1,it1,ib1,i
         write (ipun,200) -1,it1,ib2,i
         write (ipun,200)  1,it2,ib1,i
         write (ipun,200) -1,it2,ib2,i
c      XY2 wagging
c         scale(ncg)=...
          if (isum.eq.0) then
           write(iptyp,630) ncg, isymb(i),i,scale(ncg)
 630      format(i4,'. ','secd XY2 , X=',a2,i2,' WAGG  ','sc=',f5.4,2x)
          else
           write(iptyp,632) ncg, isymb(i),i,scale(ncg)
 632      format(i4,'. ','ring sXY2, X=',a2,i2,' WAGG  ',
     1    'sc=',f5.4,2x)
          endif
c
         ncg=ncg+1
         write (ipun,100)  ncg,1,it1,ib1,i
         write (ipun,200) -1,it1,ib2,i
         write (ipun,200) -1,it2,ib1,i
         write (ipun,200)  1,it2,ib2,i
c      XY2 twisting
c         scale(ncg)=...
          if (isum.eq.0) then
           write(iptyp,640) ncg, isymb(i),i,scale(ncg)
 640      format(i4,'. ','secd XY2 , X=',a2,i2,' TWIST ','sc=',f5.4,2x)
          else
           write(iptyp,642) ncg, isymb(i),i,scale(ncg)
 642      format(i4,'. ','ring sXY2, X=',a2,i2,' TWIST ','sc=',f5.4,2x)
          endif
c
       else if (nterm.eq.1) then
         it1=neigterm(1)
c
         ncg=ncg+1
         write (ipun,100) ncg,1,ib1,it1,i
         write (ipun,200) -1,ib2,it1,i
c      XY rocking
c
c         scale(ncg)=...
          if (isum.eq.0) then
           write(iptyp,710) ncg, isymb(i),i,isymb(it1),it1,scale(ncg)
 710       format(i4,'. ','secd XY,',a2,i2,'-',a2,i2,' ROCK  ',
     1     'sc=',f5.4,2x)
          else
           write(iptyp,712) ncg, isymb(i),i,isymb(it1),it1,scale(ncg)
 712       format(i4,'. ','ring sXY',a2,i2,'-',a2,i2,' ROCK  ',
     1     'sc=',f5.4,2x)
          endif
c
         ncg=ncg+1
         write (ipun,300) ncg,1,it1,ib1,ib2,i
          scale(ncg)=0.7
          if (isum.eq.0) then
           write(iptyp,720) ncg, isymb(i),i,isymb(it1),it1,scale(ncg)
 720       format(i4,'. ','secd XY,',a2,i2,'-',a2,i2,' OUT   ',
     1     'sc=',f5.4,2x)
          else
           write(iptyp,722) ncg, isymb(i),i,isymb(it1),it1,scale(ncg)
 722       format(i4,'. ','ring sXY',a2,i2,'-',a2,i2,' OUT   ',
     1     'sc=',f5.4,2x)
          endif
 300     format(i4,2x,'K',3x,i4,'.',5x,'OUT ',6x,4(i4,'.',5x))
        end if
      end
c
      Subroutine ncollin (natom,nbond,ibond,ibd,i,it1,ifour,found,iy1,
     1  iy2,iy3)
      implicit real*8 (a-h,o-z)
c          This routine tries to locate an atom (if possible bound to
c            ibd) which is not collinear with the i-ibd bond and is not
c            identical with the terminal atom it1.  The result is ifour
c
       logical found
       parameter (obtuse=2.8d0)
       parameter(mxneig=6)
       dimension nbond(natom),ibond(mxneig,natom)
        ifour=-1
        pi=acos(-1.d0)
           found=.false.
           do 400 inn=1,nbond(ibd)
             ine=ibond(inn,ibd)
c  FG   omit neighbor i of ibd
             if(ine.eq.i) goto 400
c...end
             call angle (i,ibd,ine,beta)
C    FG here just for safety, important below for atoms of
c    general position     if (beta.lt.obtuse) then
             if (beta.lt.obtuse.and.beta.gt.(pi-obtuse)) then
               ifour=ine
               found=.true.
               go to 600
             end if
 400       continue
c         no BOUND non-collinear atom found; now try any atom
           do 500 inn=1,natom
              if (inn.eq.i.or.inn.eq.ibd.or.inn.eq.it1) go to 500
c  the following is relevant when this routine is called in the
c  first part of primcoor, trying to find a unique atom in an xy3
c  grouping: omit the y atoms themselves
c   (in other cases the routine is called by iy's = -1)
              if(inn.eq.iy1.or.inn.eq.iy2.or.inn.eq.iy3) goto 500
              call angle (inn,ibd,i,beta)
c
C FG  checking ANY atom, exclude parallel case  if (beta.lt.obtuse) then
              if (beta.lt.obtuse.and.beta.gt.(pi-obtuse)) then
c
                ifour=inn
                found=.true.
                go to 600
              end if
 500       continue
c FG if this routine was called by primcoor to find a unique atom
c  within an xy3 grouping, then only part of the molec is linear
c  think of ch3cn, e.g. This is indicated by non-zero iy1(or iy2,iy3)
           if(iy1.gt.0) goto 600
c...       linear molecule
          print *, 'linear molecule, needs a dummy atom'
          print *, 'dummy added automatically, but add it also to'
          print *, 'Texas input, with ZERO charge'
          ifour=natom+1
 600      continue
       end
c
      Subroutine torsions (natom,nbond,ibond,iring,ipun,iptyp,ncg,
     1 nb,ian,li,lj,rr,fcb,fcq,nbond1,ibond1,isymb,ibord,scale)
c     This routine generates the torsions, except the ring
c       torsions which are generated in ringcoor
c
c..FG  an estimate of the torsional force const is also made
c      and put into fcq(ncg)
c     arrays ian -atomic numbers-, and rr - bond lengths - are not
c     used yet, but may be useful for more sophisticated formula
c..
      implicit real*8 (a-h,o-z)
      character*2 isymb
      parameter(mxneig=6)
      parameter(mxnr=40,mxdr=30)
      parameter (obtuse=2.8d0)
c    This is about 160 degrees and it is the limit for accepting a
c     torsional angle
      logical ifirst
      dimension nbond(natom),ibond(mxneig,natom),iring(mxnr,*)
      dimension ian(natom),li(nb),lj(nb),rr(nb),fcb(nb)
      dimension fcq(*), scale(*)
      dimension nbond1(natom),ibond1(mxneig,natom)
      dimension isymb(natom)
c
c FG   new: in case of two collinear bonds, the central atom
c      will be skipped and a new bondlist prepared in ibond1
c TEMP
c      print *, 'at entering torsions, nbond, nbond1'
c      write(*,1111) (nbond(i),i=1,natom)
c      write(*,1111) (nbond1(i),i=1,natom)
c1111   format(20i3)
c end TEMP
       do 10 i=1,natom
        nbond1(i)=nbond(i)
        do 10  k=1, nbond(i)
        ibond1(k,i)=ibond(k,i)
   10  continue
       do 80 i=1,natom
c      if(nbond(i).le.1) goto 80
c  terminal atoms cannot be the looked-for central atom
c newer: this scheme was devised for, e.g., methyl-phenyl-acetylene
c  it is wrong, however, in something like AX4OH tetragon.pyramid
c  where the central atom is part of another linear unit(equatorial)
c  but must not be eliminated if we want torsion like XAOH
c  therefore: consider only divalent atoms in the following test:
       if(nbond(i).ne.2) goto 80
       do 70 j=2,nbond(i)
       do 70 k=1,j-1
c  the updated bondlist ibond1 is needed here
       ngb1=ibond1(j,i)
       ngb2=ibond1(k,i)
       call angle (ngb1,i,ngb2,alpha)
       if(alpha.gt.obtuse) then
c   i is the central atom of a collinear unit and will be replaced
c   in the modified bondlist
         nbond1(i)=0
         do 20 jk=1,nbond(ngb1)
           if(ibond1(jk,ngb1).eq.i) then
            ibond1(jk,ngb1)=ngb2
            goto 22
           endif
 20      continue
 22      continue
         do 30 jk=1,nbond(ngb2)
           if(ibond1(jk,ngb2).eq.i) then
            ibond1(jk,ngb2)=ngb1
            goto 32
           endif
 30       continue
 32       continue
       endif
 70     continue
 80    continue
c FG the above was added to treat collinear bonds.
      do 1000 i=1,natom
       if (nbond1(i).le.1) go to 1000
c      Leave out terminal atoms
c  leave also out the central atom of a linear unit, for which
c   nbond1 was put 0
       do 800 ii=1,nbond1(i)
          j=ibond1(ii,i)
          if (j.gt.i) go to 800
c TEMPOR
c      print *, 'i=,j=', i,j
c         This is to avoid double counting of the torsions
c         Check if i and j are in the same ring
          do 200 kk=1,mxnr
             if (iring(kk,i).eq.1.and.iring(kk,j).eq.1) go to 800
c         Leave out bonds which are in the same ring
 200      continue
c.... TEMP  Start of insert
c NEW:  torsions will be normalized by 1/n rather than 1/sqr(n)
c  to achieve this the coordinate will be multiplied by 1/sqr(n)
c  and the original standard normalization is kept in the BMAT
c  routine (used in several programs)
          compnt=0.
          renorm=-1.
c TEMP
c      print *, 'before 7000 in torsions, nbond, nbond1'
c      write(*,1111) (nbond(ipr),ipr=1,natom)
c      write(*,1111) (nbond1(ipr),ipr=1,natom)
c end TEMP
          do 7000 jneig=1,nbond1(j)
c           scan the neigbors of j
            k=ibond1(jneig,j)
            if (k.eq.i) go to 7000
c TEMPOR
c      print *, 'k=', k
            do 6000 ineig=1,nbond1(i)
               l=ibond1(ineig,i)
c  TEMP
c      print *, 'l=',l
               if (l.eq.j) go to 6000
            compnt=compnt+1.
c  TEMP
c      print *, 'l=',l
 6000        continue
 7000      continue
           if(compnt.gt.0.9) then
            renorm=1./dsqrt(compnt)
           endif
          ifirst=.true.
c
          if(nbond1(j).eq.0) goto 800
c
c.... this is to avoid the generation of a torsion in, e.g., CH3CN
c     after the middle C has been taken off from the list, N has
c     no neigbor
c
          do 700 jneig=1,nbond1(j)
c           scan the neigbors of j
            k=ibond1(jneig,j)
            if (k.eq.i) go to 700
c TEMPOR
c      print *, 'k=', k
c           Check for collinearity
c       with the above new treatment of collinearity this
c       check is  perhaps obsolete, but may be important in
c       some very exotic case (one connection non-linear,
c       but another around the same bond practically linear)
             call angle (i,j,k,alpha)
             if (alpha.gt.obtuse) go to 700
            do 600 ineig=1,nbond1(i)
               l=ibond1(ineig,i)
               if (l.eq.j) go to 600
c  TEMP
c      print *, 'l=',l
c           Check for collinearity
c           see note above at k
              call angle (j,i,l,beta)
              if (beta.gt.obtuse) go to 600
               if (ifirst) then
c I have removed the normalization coefficient.
                  ncg=ncg+1
                  write(ipun,410) ncg,l,i,j,k
           write(iptyp,450) ncg,isymb(i),i,isymb(j),j,scale(ncg)
 450     format(i4,'. ','torsion ',a2,i2,'-',a2,i2,
     1     ' TORS  ','sc=',f5.4,2x)
                  call ftors(i,j,natom,nb,ian,li,lj,rr,fcb,fcq(ncg))
c   with renormalized torsions,
c  a possibility would be:     fcq(ncg)=fcq(ncg)*compnt
                  ifirst=.false.
               else
                  write (ipun,500) l,i,j,k
c400     format(i4,2x,'K',1x,f9.7,3x,'1.',5x,'TORS',6x,4(i4,'.',5x))
 400     format(i4,2x,'K',3x,'1.',5x,'TORS',6x,4(i4,'.',5x))
 410     format(i4,2x,'K',3x,3x,'1.',5x,'TORS',6x,4(i4,'.',5x))
 500     format(' ',9x,3x,'1.',5x,'TORS',6x,4(i4,'.',5x))
               end if
 600        continue
 700      continue
 800    continue
 1000 continue
      end

