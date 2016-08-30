      Subroutine fring(nr,type,ityp,natom,nb,itemp,ian,li,lj,rr,fcb,
cversion 5.0
     1 forcon,sci,aromat)
c input: nr - dimension of the ring (no of atoms)
c        type - BEND or TORS
c        ityp = 1,2,3,.. serial no of coord within given type
c        natom - no of atoms in the whole molecule
c        itemp(1..nr) - indices of the ring-atoms
c        ian(1..natom) - atomic numbers (for the whole molec)
c        li,lj,rr(1..nb) - indices of atom i and j, defining bond ij,
c             and its length, all created in subr bonds
c        fcb(1..nb) estimated force constants for bond-stretchings,
c             coming from subr fbonds
c output: forcon -scalar, the force const estimated for the ring coord
c             inate in question
c
      implicit real*8 (a-h,o-z)
      character*4 aromat
      parameter(mxdr=30,nscal=20)
      parameter(drlim=0.04,rarom1=1.38,rarom2=1.42)
c       drlim for checking aromaticity
      dimension itemp(nr),ian(natom),li(nb),lj(nb),rr(nb),fcb(nb)
      dimension ftemp(mxdr),rtemp(mxdr)
c       temporary arrays of force constants and bond length in the ring
      dimension scaldat(nscal)
      character*4 type,bend,tors
      data fref,sctors/5.,1./
c     aug, 94: sctors multiplies the torsional forcon
c     here in rings they seem OK,thus 1.
c     but see subr ftors
      data scaldat /0.785,1.0,0.790,0.96,0.808,0.768,14*0.8/
c  above scale factors: cyclobutane bend,tors;c-hexane bend,tors;
c     benzene bend,tors;
chr+
cd     write(*,*) 'Start of fring ',type,' ',nr,' ',nb
chr-
      aromat='NOAR'
c     first thing to do: identify the bonds in the ring with bonds on
c      the list for the whole molec (li, lj) and assign the correspond.
c      force constants and bond lengths
      do 100 kr=1,nr
        if(kr.lt.nr) then
          kr1=kr+1
        else
          kr1=1
        endif
        do 200 kb=1,nb
          if(li(kb).eq.itemp(kr).and.lj(kb).eq.itemp(kr1)) then
            ftemp(kr)=fcb(kb)
            rtemp(kr)=rr(kb)
            goto 100
          else if(lj(kb).eq.itemp(kr).and.li(kb).eq.itemp(kr1)) then
            ftemp(kr)=fcb(kb)
            rtemp(kr)=rr(kb)
            goto 100
          endif
  200       continue
C            print *, 'bond not found in subr fring'
C            print *, 'kb=',kb,'kr=',kr,'nb=',nb,'nr=',nr
            call bummer('bond not found in subr fring',0,2)
            STOP
  100 continue
chr+
cd     write(*,*) 'first half of fring'
chr-
c    first simple scheme: adjustment of ring force constants
c    (to be done later below) based on stretching
c    force constants only; later, bond lengths may also
c    be used
c...
      sum=0.
      sumr=0.
      do 300 kr=1,nr
      sum=sum+ftemp(kr)
      sumr=sumr+rtemp(kr)
 300  continue
      favr=sum/nr
      adjust=favr/fref
      ravr=sumr/nr
c  above adjustment is a primitive formal test,try other forms
c   it may also be different for bendings and torsions
c...
c      in estimating force constants below, 4- and 5-membered
c      rings are treated individually, all others from 6-
c      upwards equally
c      there is,of course, a basic difference between bendings and
c      torsions. the 1-dim species (present in even-rings only)
c      is also treated separately from the 2-dim specii, but
c      this may be superfluous
c...
      bend='BEND'
      tors='TORS'
      if(nr.eq.4) then
c..
c   values taken from cyclobutane
c..
        if(type.eq.bend) then
          forcon=1.89
          sci=scaldat(1)
        else if(type.eq.tors) then
          forcon=0.11
          sci=scaldat(2)
        else
C      print *, 'no match for type in subr fring'
        endif
      else if(nr.eq.5) then
c..
c   values are just a guess
c..
        if(type.eq.bend) then
            forcon=1.5
          sci=scaldat(1)
        else if (type.eq.tors) then
            forcon=0.3
          sci=scaldat(2)
        else
C      print *, 'no match for type in fring'
        endif
      else if(nr.ge.6) then
c...
c     a simlified check for aromaticity:
        if(mod(nr,2).eq.0) then
          drmax=0.
          do 410 kr=1,nr
            rtp=rtemp(kr)
            if(rtp.lt.rarom1.or.rtp.gt.rarom2) goto 420
            dr=abs(rtemp(kr)-ravr)
            if(dr.gt.drmax) drmax =dr
  410     continue
          if (drmax.gt.drlim) goto 420
          aromat='AROM'
  420     continue
         endif
c..
c    for force consts all rings from 6 and larger treated equally
c    values based on cyclohexane; bends are similar in benzene too,
c    torsions should be OK after adjustment below
c..  aromaticity checked only for estimating scale factors,not for
c    the force constants.
c...
         if(type.eq.bend) then
           if(aromat.eq.'AROM') then
             sci=scaldat(5)
           else
             sci=scaldat(3)
           endif
           if(ityp.eq.2) then
             forcon=1.4
           else if (ityp.eq.1) then
             forcon=1.5
           endif
         else if(type.eq.tors) then
           if(aromat.eq.'AROM') then
             sci=scaldat(6)
           else
             sci=scaldat(4)
           endif
           if(ityp.eq.2) then
             forcon=0.25
           else if (ityp.eq.1) then
             forcon=0.33
           endif
         else
C         print *, 'no match for type in subr fring'
         endif
      endif
c
      if(type.eq.tors) then
      forcon=adjust*forcon
c     Aug ,94: a further possible adjustment
c     if torsional fconstns should be increased:
      forcon=sctors*forcon
      else
c..
c  based on cyclohexane and benzene, bends seem not to require
c  adjustment, while torsions are more sensitive
c..
      continue
      endif
chr+
cd     write(*,*) 'End of fring'
chr-
      return
      end
c...
      Subroutine ftors(ni,nj,natom,nb,ian,li,lj,rr,fcb,
     1 forcon)
      implicit real*8 (a-h,o-z)
c..
c     goal: estimate a force const for a torsion which is
c           not part of a ring  (for rings see subr fring)
c..
c input: ni,nj - indices of the atoms defining the central bond
c        natom - no of atoms in the whole molecule
c        nb - no of bonds in the molec, determined in subr bonds
c        ian(1..natom) - atomic numbers (for the whole molec)
c        li,lj,rr(1..nb) - indices of atom i and j, defining bond ij,
c             and its length, all created in subr bonds
c        fcb(1..nb) estimated force constants for bond-stretchings,
c             coming from subr fbonds
c output: forcon -scalar, the force const estimated for the torsion
c
      dimension ian(natom),li(nb),lj(nb),rr(nb),fcb(nb)
      data sctors/3./
c      Aug, 94: torsional force constants were too small
c      for geom optimzns; a physically good value at the equil.
c      may cause too large steps, be more conservative:
c      the above factor simply multiplies forcon at the end
c..
c     first thing to do:  find the bond of torsion  on
c      the list for the whole molec (li, lj) and assign the correspond.
c      force constant and bond length (the latter is actually not
c      used in the formula below, but may be useful later)
c..
       fbond=0.
       rbond=0.
        do 200 kb=1,nb
          if(li(kb).eq.ni.and.lj(kb).eq.nj) then
            fbond=fcb(kb)
            rbond=rr(kb)
            goto 200
          else if(lj(kb).eq.ni.and.li(kb).eq.nj)  then
            fbond=fcb(kb)
            rbond=rr(kb)
            goto 200
          endif
  200   continue
       if(rbond.lt.0.001) then
C          write(*,*) 'bond not found in subr ftors'
C      print *, 'kb=',kb,'nb=',nb
c
Caution,new in subr torsions: if an i-j-k unit is collinear, j is
c  skipped and i-k is used to define the central 'bond' for a torsion
c  this is then not found here as a bond
c  temporary solution: a fixed force const will be used in this
c  case; the value of 0.2 is reasonable in c=c=c, but too large if
c  one of the bonds is single
c
       forcon=0.2
C     write(*,*) 'CAUTION, it is assumed that a collinear unit'
C     write(*,*)  'x-y-z was found, and x--z used to define the'
C     write(*,*)  'torsion; check the force constant!, f=', forcon
c
       else
c....
c    first simple scheme: scaling based on the stretching force const
c     of the central bond
c    for the simple formula, torsional force const around single and
c    double bonds in ethylene, butadiene, hexatriene, glyoxal,acrolein
c    were considered, and it seems to work roughly also for formamide
c...
c         forcon=(fbond/5.)**2*0.033
c         the tors forconst for a single bond of fstre=5 is about .033
c.....
c....new: the above is OK near the equil., but is dangerously low
c  in general; the new 1/n normalization (see subr. torsions)
c      also requires higher value
c
        forcon=(fbond/5.)**2*0.1
       endif
      forcon=sctors*forcon
      return
      end
c
      Subroutine polycycl(natom,nrings,iring,ipun,iptyp,ncg,fcq,nbond1,
     1ibond1,isymb,ibord,scale)
c
c  generating internal coordinates in condensed and spiro-ring systems
c  to describe the relative motions of rings (within each ring,
c  coordinates were generated earlier by ringcoor, called by ring)
c  natom -no of atoms in the mol; nrings -no of rings; iring - matrix
c  containing the basic info on rings, see subr ring;
c  fcq -vector of force constants estimated; nbond1 and ibond1
c  are here the bond list for ring-atoms (only the ring atoms), see
c  subr chainrin and ring
c
c  parameters used also in other routines:
c  max no of rings, max no of neighbors for an atom:
      parameter(mxnr=40,mxneig=6)
c  internal param:
c     max no of atoms on a joint (common edge):
      parameter(maxjoint=8)
c     2 neigbours of an edge in a given ring:
      parameter(mx2=2)
c     max no of propellane type edges:
      parameter(mxprop=7)
      implicit real*8 (a-h,o-z)
      character*2 isymb
      dimension iring(mxnr,*), fcq(*), scale(*)
      dimension isymb(natom)
      dimension nbond1(natom), ibond1(mxneig,natom)
      dimension joint(maxjoint)
      dimension jointnbi(mx2),jointnbj(mx2),jointnbk(mx2)
c these will be the indices of the joint (common) atoms and their
c neigbors in rings ir and jr (and ring kr in a propellane type case)
      dimension iprop(mxprop,3)
      do 5 i=1,maxjoint
        joint(i)=0
5     continue
      do 8 i=1,mx2
        jointnbi(i)=0
        jointnbj(i)=0
        jointnbk(i)=0
8     continue
      nprop=0
      do 5000 ir=1,nrings-1
      do 5000 jr=ir+1,nrings
       if(nprop.gt.0) then
c         skip the ring pair ir,jr if they were treated previously
c         as part of a propellane type unit
        irprop=0
        jrprop=0
        do 120 kprop=1,nprop
         kp1=iprop(kprop,1)
         kp2=iprop(kprop,2)
         kp3=iprop(kprop,3)
         if(ir.eq.kp1.or.ir.eq.kp2.or.ir.eq.kp3) irprop=1
         if(jr.eq.kp1.or.jr.eq.kp2.or.jr.eq.kp3) jrprop=1
         if(irprop.eq.1.and.jrprop.eq.1) goto 5000
c          already handled propellane rings skipped
 120    continue
       endif
       ijt=0
       do 100 k=1,natom
         if(iring(ir,k).eq.1.and.iring(jr,k).eq.1) then
            ijt=ijt+1
            joint(ijt)=k
c              atom k is a common atom for rings ir and jr
c
C          print *, 'anellation atoms',k
          endif
 100   continue
       njoints=ijt
       ibiph=0
       if(njoints.gt.0) goto 140
c....
c  new:  treat biphenyl-type systems also here
c   check for biphenyl-type connections follows:
c
       do 110 k=1,natom
        if(iring(ir,k).eq.1) then
         do 112 l=1,natom
          if(iring(jr,l).eq.1) then
           do 115 lnb=1,nbond1(l)
            if(ibond1(lnb,l).eq.k) then
c   biphenyl type connection between k and l found:
             ibiph=k
             jbiph=l
             goto 130
            endif
115        continue
          endif
112      continue
        endif
110    continue
130    continue
       if(ibiph.ne.0) then
c
c newer: problem in phenanthrene, biphen type connctn found between
c     rings 1. and 3.; these cases must be neglected by checking
c     if ibiph and jbiph are members of a third ring other than
c     ir and jr; before that,
c     also neglect the case when ibiph are quaterner:
c     they are then spiro-centers  -hopefully spiro works then
c
       if(nbond1(ibiph).eq.4.and.nbond1(jbiph).eq.4) then
C        print *, 'observation in polycycl:'
C        print *, 'the atoms of a biphenyl type connection are'
C        print *, 'quaterner, i.e. they are spiro-centers'
C        print *, 'ibiph=', ibiph, 'jbiph=', jbiph
         goto 5000
       endif
       do 150 lr=1,nrings
        if(lr.eq.ir.or.lr.eq.jr) goto 150
        if(iring(lr,ibiph).eq.1.and.iring(lr,jbiph).eq.1) then
C        print *, 'observation in polycycl:'
C        print *, 'the atoms of a biphenyl type connection are'
C        print *, 'members of a third ring, lr; biphenyl connection'
C        print *, 'is neglected and polycondensed case assumed'
C        print *, 'ir=',ir,'jr=',jr,'ibiph=',ibiph,'jbiph=',jbiph
C        print *, 'the third ring: lr=', lr
        goto 5000
        endif
150    continue
c
C      print *, 'biphenyl type connection between :', ibiph,jbiph
        il=0
        jl=0
       do 135 ijn=1,3
Caution: both ibiph and jbiph MUST have 3 and only 3 ring neigbors here
        if(ibond1(ijn,ibiph).ne.jbiph) then
         il=il+1
         jointnbi(il)=ibond1(ijn,ibiph)
        endif
        if(ibond1(ijn,jbiph).ne.ibiph) then
         jl=jl+1
         jointnbj(jl)=ibond1(ijn,jbiph)
        endif
135    continue
       i1=jointnbi(1)
       i2=jointnbi(2)
       j1=jointnbj(1)
       j2=jointnbj(2)
c ROCK of the biphenyl-type rings relative to each other:
       ncg=ncg+1
       write(ipun,400) ncg,jbiph,i1,ibiph
c      fcq(ncg)=....
       write(ipun,420) jbiph,i2,ibiph
c      scale(ncg)=....
       write(iptyp,430) ncg,isymb(jbiph),jbiph,scale(ncg)
430     format(i4,'. ','biphenyl-type',a2,i2,' ROCK  ','sc=',f5.4,2x)
       ncg=ncg+1
       write(ipun,400) ncg,ibiph,j1,jbiph
c      fcq(ncg)=....
       write(ipun,420) ibiph,j2,jbiph
c      scale(ncg)=....
       write(iptyp,430) ncg,isymb(ibiph),ibiph,scale(ncg)
c   OUT follows:
       ncg=ncg+1
       write(ipun,500) ncg,jbiph,i1,i2,ibiph
500    format(i4,2x,'K',3x,' 1.',7x,'OUT ',6x,4(i4,'.',5x))
c      fcq(ncg)=....
c      scale(ncg)=....
       write(iptyp,510) ncg,isymb(jbiph),jbiph,scale(ncg)
510     format(i4,'. ','biphenyl-type',a2,i2,' OUT   ','sc=',f5.4,2x)
c   OUT follows:
       ncg=ncg+1
       write(ipun,500) ncg,ibiph,j1,j2,jbiph
c      fcq(ncg)=....
c      scale(ncg)=....
       write(iptyp,510) ncg,isymb(ibiph),ibiph,scale(ncg)
       endif
       goto 5000
140    continue
C      print *, 'no of joints,njoints=', njoints
       if(njoints.le.2) then
c              /n.m.0/ type bicyclo  systems okay this way,
c              /n.m.k/ with k>1 NOT YET
c
c           we need also the neighbors to define the coordinates  below:
         irnb=0
         jrnb=0
         krnb=0
         do 50 ijt=1,njoints
          jat=joint(ijt)
          do 60 k=1,nbond1(jat)
c           a neighbor of the joint atom jat:
           knb=ibond1(k,jat)
           do 70 kk=1,njoints
            if(knb.eq.joint(kk)) goto 60
 70        continue
c  leave out neighbors that are joint (annelation) atoms  themselves
c
c  find out if the neighbor knb is in ring ir or jr:
           if(iring(ir,knb).eq.1) then
            irnb=irnb+1
            jointnbi(irnb)=knb
           else if(iring(jr,knb).eq.1) then
            jrnb=jrnb+1
            jointnbj(jrnb)=knb
           else
C      print *, 'caution, in subr polycycl a ring neighbor'
C      print *, 'doesnot belong to either of the 2 annelated rings'
C      print *, 'more than two rings meet here?'
C      print *, 'propellane type is treated, but not yet safe'
c   first we assume propellane but this will finally be checked
c  further below
            krnb=krnb+1
            jointnbk(krnb)=knb
           endif
 60       continue
 50      continue
         jt1=joint(1)
         jt2=joint(2)
         i1=jointnbi(1)
         i2=jointnbi(2)
         j1=jointnbj(1)
         j2=jointnbj(2)
c  the following two will be used only for propellane case later below
         k1=jointnbk(1)
         k2=jointnbk(2)
c
         if(njoints.eq.1) then
c  spiro compound:
c
c  rock:
          ncg=ncg+1
          write(ipun,400) ncg,i1,j1,jt1
c   add force const later  fcq(ncg)=
          write(ipun,410) i2,j1,jt1
          write(ipun,420) i1,j2,jt1
          write(ipun,420) i2,j2,jt1
c            scale(ncg)=....
             write(iptyp,405) ncg,isymb(jt1),jt1,scale(ncg)
  405     format(i4,'. ','spiro cntr at',a2,i2,' ROCK  ','sc=',f5.4,2x)
c
c  wagg(or twist):
          ncg=ncg+1
          write(ipun,400) ncg,i1,j1,jt1
c    add force const later  fcq(ncg)=
          write(ipun,420) i2,j1,jt1
          write(ipun,410) i1,j2,jt1
          write(ipun,420) i2,j2,jt1
c            scale(ncg)=....
             write(iptyp,407) ncg,isymb(jt1),jt1,scale(ncg)
  407     format(i4,'. ','spiro cntr at',a2,i2,' WG/TW ','sc=',f5.4,2x)
c
c  twist(or wagg):
          ncg=ncg+1
          write(ipun,400) ncg,i1,j1,jt1
c    add force const later  fcq(ncg)=
          write(ipun,420) i2,j1,jt1
          write(ipun,420) i1,j2,jt1
          write(ipun,410) i2,j2,jt1
c            scale(ncg)=....
             write(iptyp,409) ncg,isymb(jt1),jt1,scale(ncg)
  409     format(i4,'. ','spiro cntr at',a2,i2,' TW/WG ','sc=',f5.4,2x)
400    format(i4,2x,'K',3x,'+1.',7x,'BEND',6x,3(i4,'.',5x))
410    format(' ',9x,'+1.',7x,'BEND',6x,3(i4,'.',5x))
420    format(' ',9x,'-1.',7x,'BEND',6x,3(i4,'.',5x))
c
         else
c njoints=2 here
c    check for propellane: are k1 and k2 (see above) in the same ring?
          krfound=0
          if(k1.gt.0.and.k2.gt.0) then
           do 40 kkr=1,nrings
            if(iring(kkr,k1).eq.1.and.iring(kkr,k2).eq.1) then
             krfound=1
             kr=kkr
c   kr is the index of the third ring(besides ir and jr)
             goto 45
            endif
  40       continue
          endif
  45      continue
          if(krfound.eq.0) then
c                 not a propellane type edge
            ncg=ncg+1
            write(ipun,200) ncg,i1,jt1,jt2,j2
c           add later: fcq(ncg)=..force const
            write(ipun,300) j1,jt1,jt2,i2
              write(iptyp,307) ncg,isymb(jt1),jt1,isymb(jt2),
     1        jt2,scale(ncg)
  307       format(i4,'. ','annelated',2(a2,i2),' BUTFLY','sc=',f5.4,2x)
 200   format(i4,2x,'K',3x,'1.',8x,'TORS',6x,4(i4,'.',5x))
 300   format(10x,'-1.',7x,'TORS',6x,4(i4,'.',5x))
          else
c                   propellane type edge, 3 rings join
           nprop=nprop+1
           iprop(nprop,1)=ir
           iprop(nprop,2)=jr
           iprop(nprop,3)=kr
c
c                   2 combinations of the 3 butterfly
c                   coordinates wil be constructed as in C3 symm
c
c   2,-1,-1 combination:
            ncg=ncg+1
            write(ipun,210) ncg,i1,jt1,jt2,j2
 210         format(i4,2x,'K',3x,'2.',8x,'TORS',6x,4(i4,'.',5x))
            write(ipun,310) j1,jt1,jt2,i2
 310         format(10x,'-2.',7x,'TORS',6x,4(i4,'.',5x))
            write(ipun,300) i1,jt1,jt2,k2
            write(ipun,320) k1,jt1,jt2,i2
            write(ipun,300) k1,jt1,jt2,j2
            write(ipun,320) j1,jt1,jt2,k2
              write(iptyp,317) ncg,isymb(jt1),jt1,isymb(jt2),
     1        jt2,scale(ncg)
  317         format(i4,'. ','3-rg,prop',2(a2,i2),' FAN   ',
     1       'sc=',f5.4,2x)
c
c   0,1,-1 combtn:
            ncg=ncg+1
            write(ipun,200) ncg,i1,jt1,jt2,k2
            write(ipun,300) k1,jt1,jt2,i2
            write(ipun,300) k1,jt1,jt2,j2
            write(ipun,320) j1,jt1,jt2,k2
 320        format(10x,'+1.',7x,'TORS',6x,4(i4,'.',5x))
              write(iptyp,317) ncg,isymb(jt1),jt1,isymb(jt2),
     1        jt2,scale(ncg)
          endif
         endif
        else
       print *, 'CAUTION:'
       print *, 'polycyclic rings with more than 2 joints'
       print *, 'like norbornane, not yet implemented'
        endif
 5000  continue
       return
       end
