      Subroutine INTC(x,ian,natom,inp,ipun,iptyp,ncg)
cversion 5.0
c.......... summer 94:
c    1) finding rings is safer now (due to Richi), which
c       does NOT mean that the coordinates are always OK
c       in macrorings
c    2a) hydrogen bonds see subr BONDH
c    2b) note that any additional bond, i.e. also H-bonds
c     can also be defined 'manually' in the input,
c     see BOND card here in INTC
c     Note for PP: a starting H-bond, which originally
c     may be non-linear, can go to linear; this
c     requires redefinition of int.coordints
c     (from BEND to LIN)
c     (rerunning INTC after a few steps)
c     It would be nice to make this AUTOMATIC
c    3) torsional fconstants increased, see subr FTORS
c     (ring torsions left unchanged, see FRING
c..............
C
C     PROGRAM FOR GENERATING 'NATURAL' INTERNAL COORDINATES
C
C         Authors: Peter Pulay and Geza Fogarasi
C
C     COPYRIGHT 1991 by the authors.
C
C     No distribution or adaptation of this program is allowed
C     without the written permission of the authors.
c
c        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     The basic reference to this work is:
c     Fogarasi, G., X. Zhou, Taylor, P.W., and Pulay, P.,
c     J.Am.Chem.Soc., vol. 114 (1992),  p. 8191.
c        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
c...  x(3,natom)  = Cartesians
c...  ian(1..natom)=atomic numbers
c...  natom= number of nuclei
c...  inp=input file; note: coordinates etc come from driver, but
c        INTC ITSELF READS ITS OWN OPTIONS,
c     the latter option cards can be put after the INTC card
c     even if intc is called from outside (by TX90)
c...  ipun=punch file for the coordinates generated
c..   iptyp=onto this file a description of coordinates
c     and a very crude estimate for scale factors will be written
c
      implicit real*8 (a-h,o-z)
      parameter (maxat=200,nmax=maxat/2,maxq=3*maxat)
      parameter(mxnr=40)
      parameter(mxneig=6)
      parameter(mxsymb=92)
      character*4 aromat
      character*4  offs,fdiag,covrad,keyw,hbond
      character*2 symb(mxsymb), isymb(maxat)
      logical dummy
      dimension x(3,natom),ian(natom)
      dimension li(maxat),lj(maxat),rr(maxat)
      dimension irad(maxat),radii(maxat)
      dimension nbond(maxat),ibond(mxneig,maxat)
      dimension nbond1(maxat),ibond1(mxneig,maxat)
      dimension iring(mxnr,maxat)
c
c..FG mxnr= max no of rings, mxneig= max no of neigbours for an atom
c...  x is the Cartesian,  ian the atomic numbers,
c    li(k) and lj(k) are the indices of the atompair forming bond k,
c      rr(k) is the bond length
c     radii - as a very special option, non-default atomic radii
c      can be read into this array from inp
c...  ibond(1,j) is the first atom j is bound to
c...  nbond(j) is the total number of atoms j is  bound to
c.... iring (i,j) is 1 if atom j is in ring i.
c
      dimension iterm(nmax),iprim(nmax),isec(nmax),itert(nmax),
     1  iquat(nmax),iterm1(nmax)
      dimension icent(nmax)

c...  iterm(1..nterm) lists the terminal atoms, iprim(1..nprim) the
c...  primary ones etc.
c...
c.FG  fcb and fb will contain estimated force constants and forces
c ..  for the bond stretching coordinates, ibord the bond order
c       (see subr fbonds)

      dimension fcb(maxat),fb(maxat),ibord(maxat)

c..   fcq will contain the force constants for the complete set of
c        coordinates, scale - scale factors for vibr.calcn.

      dimension fcq(maxq), scale(maxq)

c...  ncg will count the number of coordinates generated

      dimension indsub(maxat), nsub(maxat)
c...  needed for subroutine subsets that locates submolecules

      data symb /' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',
     1 'NA','MG','AL','SI',' P',' S','CL','AR',74*'XX'/
      data offs,fdiag,hbond/'NOOF','FDIA','NOHB'/
      data blownew,nrad /-1., 0 /
c
c     check dimension:
      if(natom.gt.maxat) then
       print *, 'the maximum number of atoms in subr INTC is set to:'
       print *, 'maxat=', maxat,  '(easy to increase)'
       print *, 'the number of atoms in the calling program is too high'
       print *, 'natom=', natom
       call bummer('natom.gt.maxat',natom,2)
      STOP
      endif
c
c  handling dummy atom for linear molecules:
c  Note: the dummy MUST have zero charge
c  CAUTION: in all routines below INTC, the number of atoms will
c  NOT include the dummy:
c
      dummy=.false.
      if(ian(natom).eq.0) then
        dummy=.true.
        natom=natom-1
      endif
c
      ncg=0
      nb=0
c
c... read options for INTC from input:
c
C FG  95 jan, next backspace and required by new fortran:
      !backspace(inp)
10    read(inp,370,end=190) keyw
370   format(a4)
      print *, 'KEYW', keyw
C FG 95 jan, also needed in connection with the above:
      if(keyw(1:2).eq.'N=') goto 10
c     this was still the last atom card
      if(keyw.eq.'B631'.or.keyw.eq.'B421') then
        offs=keyw
        goto 10
      else if(keyw.eq.'NOFD') then
        fdiag=keyw
        goto 10
      else if(keyw.eq.'BLOW'.or.keyw.eq.'RADI') then
        print *, 'reading RADI'
c   this is a special option to change covalent radii (in routine BONDS)
        !backspace (inp)
        if(keyw.eq.'BLOW') then
          read(inp,390) blownew
390     format(10x,f10.5)
        else
          read(inp,390) xnrad
          print *, xnrad
          nrad=xnrad+0.1
        endif
        goto 10
c  new: for specifying extra bonds, better than above:
c      simply define them in the input by indeces of atoms
c
      else if(keyw.eq.'BOND') then
        backspace (inp)
        read(inp,390) xnb
C      this start value of nb is the number of extra bonds
        nb=xnb+0.1
        do 20 inb=1,nb
          read(inp,395) xl,xj
          li(inb)=xl+0.1
          lj(inb)=xj+0.1
 20     continue
        goto 10
      else if(keyw.eq.'HALL'.or.keyw.eq.'HSUB') then
          hbond=keyw
        goto 10
      else
c   no more option cards found
        backspace (inp)
      endif
      if(nrad.gt.0) then
c            nrad number of atomic radii will be read,
c            these will replace default values in BONDS
      do 210 i=1,nrad
        print *, 'reading radii'
        read(inp,*) xind, radii(i)
        !read(inp,395) xind, radii(i)
395     format(2f10.5)
        irad(i)=xind+0.1
210   continue
      endif
190   continue

c...put chem symbols into isymb:
      do 160 i=1,natom
      iatno=ian(i)
      isymb(i)=symb(iatno)
 160  continue

c     defaults:
      do 110 i=1,maxq
      fcq(i)=1.
      scale(i)=0.8
 110  continue
      do  100 i=1,natom
         call cartes (.true.,i,x(1,i),x(2,i),x(3,i),ian(i))
 100  continue
c... cartes stores the Cartesians and the atomic numbers in an own
c... array if called with a .true. first argument
c... it recovers the same quantities if the first argument is false.
c... This is to avoid a COMMON and too many arguments
c
      call bonds (x,ian,natom,ipun,iptyp,maxat,nb,li,lj,rr,ncg,
     1fcb,fb,offs,symb,ibord,scale,blownew,nrad,irad,radii)
c
      if(hbond.eq.'HALL'.or.hbond.eq.'HSUB') then
c
c    new: in connection with H-bonds: establish bondlists
c       first WITHOUT H-bonds, put lists temorarily into nbond1, ibond1
c       after BONDH these routines will be called again,
c       and the complete lists put into nbond, ibond;
c       nbond1, ibond1 will be overwritten after that
c
C       print *, '  '
        call bondlist (nb,li,lj,natom,nbond1,ibond1)
        call terminal (natom,nterm,iterm,nbond1,ibond1)
c
c     check for isolated submolecules in the system:
c     the problem arouse in chainrin, in the connectivity tests there
c     and,at present, calling subsets has consequences for chainrin only
c     .... new: subsets is now used in connection with H-bonds (BONDH)
c
          call subsets(natom,nbond1,ibond1,msub,indsub,nsub)
c
          call BONDH (X,ian,isymb,natom,nb,li,lj,rr,nbond1,ibond1,
     1             keyw,indsub,nsub,msub,fcb,fb,ncg,ipun,iptyp,iout)
        endif
c t6272
c
c         end of treating H-bonds, now establish the complete lists
c         in nbond,ibond
c         nbond1 and ibond1 will be overwritten later below!
c
      call bondlist (nb,li,lj,natom,nbond,ibond)
      call terminal (natom,nterm,iterm,nbond,ibond)
c               subsets will be updated too, may be unnecessary:
      call subsets(natom,nbond,ibond,msub,indsub,nsub)
c
      if (nterm.eq.0) go to 130
      nterm1=nterm
      do 150 ii=1,nterm
        iterm1(ii)=iterm(ii)
 150  continue
      do 180 ii=1,natom
        nbond1(ii)=nbond(ii)
        do 170 jj=1,nbond(ii)
           ibond1(jj,ii)=ibond(jj,ii)
 170    continue
 180  continue
 200  call prune (natom,nbond1,ibond1,nterm1,iterm1)
      call terminal (natom,nterm1,iterm1,nbond1,ibond1)
      if (nterm1.gt.0) go to 200
c
c     after the above pruning, only ring atoms and chains between
c     rings survive in nbond1, ibond1
c
c NEW for substructures:

      do 135 istruct=1,msub
      call chainrin (natom,nbond1,ibond1,indsub,nsub,istruct)
 135  continue
c
 130  continue
c..
c FG after chainrin,chains between rings have been eliminated and
c the bondlists nbond1 and ibond1 refer to ring c atoms only.
c ( these lists will be changed in subr ring, but restored)
c..
c..FG original version did not work for a 'bare' ring, which
c   has no atoms attached to the ring, i.e. no terminal atoms
c...
      if(nterm.gt.0) then
chr+
c        call ring (natom,nbond1,ibond1,nbond,ibond,ipun,iptyp,nrings,
        call hring (natom,nbond1,ibond1,nbond,ibond,ipun,iptyp,nrings,
chr-
     1iring,ncg,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
        if(nrings.gt.1) call polycycl(natom,nrings,iring,ipun,iptyp,
     1ncg,fcq,nbond1,ibond1,isymb,ibord,scale)
        else
c.. in case of a bare ring, nbond1=nbond, ibond1=ibond
chr+
c        call ring (natom,nbond,ibond,nbond,ibond,ipun,iptyp,nrings,
        call hring (natom,nbond,ibond,nbond,ibond,ipun,iptyp,nrings,
chr-
     1iring,ncg,nb,ian,li,lj,rr,fcb,fcq,isymb,ibord,scale,aromat)
        if(nrings.gt.1) call polycycl(natom,nrings,iring,ipun,iptyp,
     1ncg,fcq,nbond,ibond,isymb,ibord,scale)
C    Bp
Caution: for 2 bare rings in biphenyltype connection
Calling torsion is necessary
C   lets keep also calling order, maybe unnecessary but ...
C       goto 400
       endif
c
      call order (natom,nterm,iterm,nbond,ibond,ipun,iptyp,iring,
     1 ncent,icent,nprim,iprim,nsec,isec,ntert,itert,nquat,iquat,
     2ncg,isymb,ibord,scale,dummy)
c
      call torsions (natom,nbond,ibond,iring,ipun,iptyp,ncg,
     1 nb,ian,li,lj,rr,fcb,fcq,nbond1,ibond1,isymb,ibord,scale)
c
C 400   continue
C     write(*,330) ncg
 330  format(/,1x,'number of coordinates generated=',I5)
       ncgok=3*natom-6
       if(ncg.ne.ncgok) then
        write(*,335) ncgok
 335    format('?CAUTION: number of coordinates is not 3N-6= ',i3)
cmd     write(ipun,335) ncgok
cmd     write(ipun,333) ncg
333     format('?actually generated',i5,'   coordinates')
       endif
C      print *, 'forces'
C     write(*,340) (fb(i),i=1,nb)
 340  format (8F10.6)
C      print *, 'force constants'
      do 140 i=1,nb
      fcq(i)=fcb(i)
 140  continue
C     write(*,340) (fcq(i),i=1,ncg)
      if(offs.eq.'B421'.or.offs.eq.'B631') then
        write(ipun,340) (fb(i),i=1,nb)
      endif
      if(fdiag.eq.'FDIA') then
      write(ipun,360) (fcq(i),i=1,ncg)
      endif
 360  format(8E10.2)
      end
c.........endof INTC
c
      Subroutine BONDS (X,IAN,natom,ipun,iptyp,IDIM,nb,LI,LJ,RR,ncg,
     1fcb,fb,offs,symb,ibord,scale,blownew,nrad,irad,radii)

      IMPLICIT REAL*8 (A-H,O-Z)
      character*4 offs
      character*2 symb(*), bsymb
      parameter(blow0=1.25)
c       tolerance for accepting atom distances as bonds
      DIMENSION X(3,IDIM), XX(3), IAN(natom)
      DIMENSION LI(*), LJ(*), RR(*)
      dimension fcb(*),fb(*), scale(*)
      dimension ibord(*)
      dimension irad(nrad),radii(nrad)
      DIMENSION IP(6), U(3), V(3), W(3), Z(3)
      DIMENSION RADIUS(19)
C
C.... IDIM IS THE CORRESPONDING DIMENSION STATEMENT IN THE CALLING ROUTI
C.... IT MAY ALSO READ AS X(IDIM),Y(IDIM),Z(IDIM) CONSECUTIVELY
C.... MAXIMUM 100  bonds... EASY TO CHANGE
C.... X CONTAINS THE CARTESIAN COORDINATES IN ANGSTROM
C.... natom=NUMBER OF ATOMS, IAN(I)=ATOMIC NUMBERS, ipun=punch FILE
c     symb - chemical symbols (for first 18 elements only)
c     blownew,nrad,irad,radii  control optional changes in the
c     covalent radii

      DATA RADIUS(1),RADIUS(2),RADIUS(3),RADIUS(4),RADIUS(5),RADIUS(6),R
     1ADIUS(7),RADIUS(8),RADIUS(9),RADIUS(10)/0.4,0.3,1.25,0.97,0.82,0.7
     25,0.69,0.66,0.62,0.4/
C     MODIFIED by Raffaele Borrelli
      DATA RADIUS(11),RADIUS(12),RADIUS(13),RADIUS(14),RADIUS(15),RADIUS
     1(16),RADIUS(17),RADIUS(18),RADIUS(19)/1.7,1.1,1.2,1.05,1.0,1.0,1.0
     2,1.0,1.3/
C
C     ORIGINAL DATA
c      DATA RADIUS(11),RADIUS(12),RADIUS(13),RADIUS(14),RADIUS(15),RADIUS
c     1(16),RADIUS(17),RADIUS(18),RADIUS(19)/1.7,1.6,1.2,1.05,1.0,1.0,1.0
c     2,1.0,1.3/
C
C.... Output arguments:
C.... nb is the number of bonds
C.... LI and LJ contain the indices of atoms connected by bonds
C.... RR the bond length
c     fcb - estimtd force constnts, scale - scale factors;
c     ibord - bond order ( just single,double..)

C     write(*,10)
   10 FORMAT (////,1X,'bonds no., atom i, atom j,    R=,  force const=,
     1force=',/)
c
c    FG
C     nb's start value is now input parameter
C      (necessary if BOND option in INTC reads extra bonds)
C             nb=0
c     Optional changes in covalent radii:
c     blowing up all radii uniformly:
      if(blownew.lt.0.0) then
        blow=blow0
      else
        blow=blownew
      endif
c
C     handle extra bonds for which indices were read in input
c      in case of option BOND
C       atom indices of extra bonds are in the first nb elements
C      of LI and LJ
      do 22 ib=1,nb
        i=li(ib)
        j=lj(ib)
        R=SQRT((X(1,i)-X(1,j))**2+(X(2,i)-X(2,j))**2+(X(3,i)-X(3,j))**2)
        II=IAN(i)
        JJ=IAN(j)
        rr(nb)=R
C
         ncg=ncg+1
         write (ipun,100) ncg,i,j
c
C    Caution: for the artificial extra bonds, routine fbonds
C      is replaced by explicit statements here:
C         call fbonds(II,JJ,R,fcb(nb),fb(nb),offs,symb,ibord(nb),
C     1   scale(ncg))
c
         fcb(ib)=2.
         fb(ib)=0.
         scale(ncg)=1.
         bsymb='..'
         write(iptyp,150) ncg, symb(ii),i,bsymb,symb(jj),j,scale(ncg)
C     write(*,110) ib,I,J,R,fcb(ib),fb(ib),scale(ib)
 22   continue
C
      DO 20 I=2,natom
      DO 20 J=1,I-1
         II=IAN(I)
         JJ=IAN(J)
         IF (II.GT.19) II=19
         IF (JJ.GT.19) JJ=19
c
c FG don't build bond with a dummy atom of zero charge
            if(ii.eq.0.or.jj.eq.0) goto 20

c    default values for radii:
         rad1=RADIUS(II)
         rad2=RADIUS(JJ)
c
c        changing individual radii:
         if(nrad.gt.0) then
           do 25 k=1,nrad
             if(i.eq.irad(k)) then
               rad1=radii(k)
               goto 25
             else if (j.eq.irad(k)) then
               rad2=radii(k)
             endif
25       continue
         endif

        R=SQRT((X(1,i)-X(1,j))**2+(X(2,i)-X(2,j))**2+(X(3,i)-X(3,j))**2)
         IF (R.GT.(blow*(rad1+rad2))) GO TO 20
         nb=nb+1
         LI(nb)=I
         LJ(nb)=J
         RR(nb)=R
         ncg=ncg+1
         write (99,*) 1.0,ncg,i,j
         write (ipun,100) ncg,1.0,i,j
C..FG
         call fbonds(II,JJ,R,fcb(nb),fb(nb),offs,symb,ibord(nb),
     1   scale(ncg))
         if (ibord(nb).eq.1) then
           bsymb=' -'
         else if(ibord(nb).eq.2) then
           bsymb=' ='
         else if(ibord(nb).eq.3) then
           bsymb='=-'
         else if(ibord(nb).eq.9) then
           bsymb='-.'
         endif
         write(iptyp,150) ncg, symb(ii),i,bsymb,symb(jj),j,scale(ncg)
150   format(i4,'.',4x,a2,i2,1x,a2,1x,a2,i2,3x,'STRE  ','sc=',f5.4,2x)
C     write(*,110) nb,I,J,R,fcb(nb),fb(nb),scale(nb)
110   format(1x,i5,3x,i5,3x,i5,2x,4f10.3)
 100     format(i4,2x,'K',3x,f10.7,2x'STRE',6x,i3,'.',6x,i3,'.')
   20 CONTINUE
c     endof BONDS
      END
c
      Subroutine BONDLIST (nb,li,lj,natom,nbond,ibond)
      parameter(mxneig=6)
      implicit real*8 (a-h,o-z)
      dimension li(*),lj(*),ibond(mxneig,natom),nbond(natom)
C     this subroutine inverts the bond list
c     input: nb=number of bonds, li(k),lj(k)= atoms on bond k,
c           natom=number of atoms,
c     output: nbond(j)= the number of bonds atom j is on
c            ibond(l,j),l=1,nbond(j)= the atoms j is bound to
c  FG ...nbond and ibond will contain the basic info on topology
c  FG    to be used in several subroutines
      do 100 j=1,natom
         nbond(j)=0
         do 100 i=1,6
           ibond(i,j)=0
 100  continue
      do 200 ib=1,nb
         j1=li(ib)
         nbj=nbond(j1)+1
         nbond(j1)=nbj
         ibond(nbj,j1)=lj(ib)
c FG...  write(*,*) 'ibond(nbj,j1)=ib', nbj,j1,ib
         j2=lj(ib)
         nbj=nbond(j2)+1
         nbond(j2)=nbj
         ibond(nbj,j2)=li(ib)
  200  continue
C     write(*,*) 'bonds for the atoms'
       do 250 j=1,natom
C        write (*,300) j,nbond(j),(ibond(k,j),k=1,nbond(j))
 250   continue
 300  format(1x,2i3,3x,6i4)
c     endof BONDLIST
      end
c
      Subroutine FBONDS(ni,nj,r,fcbi,fbi,offs,symb,ibordi,sci)
c
c     estimating force consts, offset force and scale factors for bonds
c...
c...  ni and nj are the atomic numbers, known in subr bonds
c...  where fbonds is called
c...  fcbi and fbi will be the force const and force,resp.,for bond i
c...  sci the scale factor
c     ibordi will be the bond order: 1 - unspecified, default;
c          1..3, single..triple,  9 - aromatic.
c...
      implicit real*8(a-h,o-z)
      parameter (maxtyp=40, maxtyp2=2*maxtyp)
      parameter(maxtypd=20,maxtypd2=2*maxtypd)
      character*2 symb(*),symb1,symb2
      character*4 type(maxtyp), type12, type21
      character*4 offs
      dimension  fconst(maxtyp), forces(maxtyp2),scalfac(maxtyp)
      dimension faromat(maxtypd2),fdouble(maxtypd2),ftriple(maxtypd2)
      logical flag
      save flag
      data flag /.true./
      data type /' C C',' C N',' C O',' C F',' N O',' N F',' O F',
     1 ' C H',' N H',' O H',30*'XXXX'/

c       estimated force constants for regular bonds: (some special
c       cases will be calculated separately farther below)
      data fconst /7.32, 7.8, 6.5, 7.3, 5., 5., 5.,
     1 5.5, 7., 8.3, 30*5./

c      forces for single,aromatic,double,triple bonds follow:
c      in the data array, first group for 4-21* basis, second
c      half for 6-31* basis
      data forces /-.05,+.02,+.02,0.,0.,0.,0.,
     1  +.06, +.02, -.02, 30*0., 0.0,+.10,+.14,0.,0.,0.,0.,
     2  +.04,+.07,+.09, 30*0./
      data faromat/+.10,+.10, 18*0., +.09,+.16, 18*0./
      data fdouble/+.23,+.26,+.20, 17*0.,+.19,+.28,+.32, 17*0./
      data ftriple/+.35,+.52, 18*0., +.34,+.47, 18*0./
c..    lower limits for cc single bond length, cc aromatic, etc.
      data ccsgl,ccarom,ccdbl, cnsgl,cnarom,cndbl, cosgl/
     1  1.42,1.38,1.28,  1.35,1.30,1.20,  1.30/
c  ..the scale factors refer to: CC single,arom,double,triple
c   CN single,arom,double,triple,  CO sgl,dbl,  CF
      data scalfac/0.920,0.911,0.86,0.86, 0.9,0.86,0.85,0.85,
     1 0.85,0.826, 0.75, 29*0.9/
c
c...defaults for unspecified cases:
      ibordi=1
      fbi=0.
      fcbi=5.
      sci=0.9
c      default in case of 421* basis:
        njump=0
        mjump=0
      if(offs.eq.'B421') then
         continue
      else if (offs.eq.'B631') then
        njump=maxtyp
        mjump=maxtypd
      else
        if(flag) then
C      print *, 'basis set unspecified for selecting offset forces'
        flag=.false.
        endif
      endif
      symb1=symb(ni)
      symb2=symb(nj)
      if(symb1.eq.'XX'.or.symb2.eq.'XX') then
C      print *, 'atomic number greater than handled'
C      print *, 'in subr fbonds   - force const for stretchings'
C      print *, 'defaults will be used'
        return
      endif
      type12=symb1//symb2
      type21=symb2//symb1
      itype=0
      do 100 i=1,maxtyp
        if(type12.eq.type(i).or.type21.eq.type(i)) then
          itype=i
          goto 200
        else
          if(type(i).eq.'XXXX') goto 200
        endif
  100 continue
  200 continue
      if (itype.eq.0) then
c      defaults for unspecified cases see above at beginning of subr
      continue
      else
c     single bonds assumed at this stage:
      fbi=forces(itype+njump)
      fcbi=fconst(itype)
      ibordi=1
      endif
c...  some special force constants will be handled separately, and
c     the value taken from data overwritten
c...  for the types of bonds see  data  array type
      if(itype.eq.1) then
c..cc
        sci=scalfac(1)
        if(r.le.ccsgl.and.r.gt.ccarom) then
          fbi=faromat(itype+mjump)
          ibordi=9
          sci=scalfac(2)
        endif
        if(r.le.ccarom.and.r.gt.ccdbl) then
          fbi=fdouble(itype+mjump)
          ibordi=2
          sci=scalfac(3)
        endif
        if(r.le.ccdbl) then
          fbi=ftriple(itype+mjump)
          ibordi=3
          sci=scalfac(4)
        endif
c
        fcbi=7.3227-34.1961*(r-1.4)+113.8478*(r-1.4)**2
      else if(itype.eq.2) then
c..cn
        sci=scalfac(5)
        if(r.le.cnsgl.and.r.gt.cnarom) then
          fbi=faromat(itype+mjump)
          ibordi=9
          sci=scalfac(6)
        endif
        if(r.le.cnarom.and.r.gt.cndbl) then
          fbi=fdouble(itype+mjump)
          ibordi=2
          sci=scalfac(7)
        endif
        if(r.le.cndbl) then
          fbi=ftriple(itype+mjump)
          ibordi=3
          sci=scalfac(8)
        endif
c
        fcbi = 7.8-33.*(r-1.35)+110.*(r-1.35)**2
      else if(itype.eq.3) then
c..co
        sci=scalfac(9)
        if(r.le.cosgl) then
          fbi=fdouble(itype+mjump)
          ibordi=2
          sci=scalfac(10)
        endif
c
        fcbi = 6.5-33.*(r-1.4)+110.*(r-1.4)**2
      else if(itype.eq.4) then
c..cf
        sci=scalfac(11)
        fcbi = 7.3 - 25.*(r-1.38)
      else
       continue
      endif
      if(offs.ne.'B421'.and.offs.ne.'B631') then
       fbi=0.
c  basis not specified, no offset force
      endif
      return
c     endof FBONDS
      end
c
c
      Subroutine BONDH (X,IAN,isymb,natom,nb,LI,LJ,RR,nbond1,ibond1,
     1keyw,indsub,nsub,msub,fcb,fb,ncg,ipun,iptyp,iout)

      IMPLICIT REAL*8 (A-H,O-Z)
      parameter(nh=2.8,oh=2.8,fh=2.6,dmax=5.0)
      parameter(mxneig=6)
      character*2 isymb(*), bsymb, type
      character*4 keyw
      DIMENSION X(3,natom), XX(3), IAN(natom)
      DIMENSION LI(*), LJ(*), RR(*)
      dimension nbond1(natom), ibond1(mxneig,natom)
      dimension indsub(*),nsub(*)
      dimension fcb(*), fb(*)
C
c     purpose: find H-bonds and add these to the bond lists
c     note: this is basically for ab initio calcns, because
c     it is assumed that in the starting geometry the H-bonds
c     as determined by the atom-pair distance are already present:
c     the possibility of a H-bond is determined by the H...X
c     distance  - see parameters nh,oh,fh above
c     in the moiety  Y-H...X  both X and H can only be N,O, or F
c
c     input arguments:
C.... X : cartesian coordinates in angstrom, ian(i):atomic numbers,
c     isymb(i) : chem. symbols (for first 18 elements only)
C.... natom=no of atoms,nb= no of bonds before H-bonds, will be changed
c     ipun,iptyp = files
C.... LI and LJ contain the indices of atoms connected by bonds
C.... RR the bond lengths,
c     nbond1,ibond1 = see BONDLISTS, before H-bonds,
c     keyw: if.keyw.eq.'HALL' all H-bonds will be formed,
c           if.keyw.eq.'HSUB' only between DIFFERENT SUBSETS
c     indsub,nsub,msub   refer to subsytems, see SUBSETS routine
c     fcb - estimtd force constnts, fb - offset forces, not really
c           important in this routine but needed
c     ncg = no of internal coordinates generated (increased here
c           by H-bond stretchings)
C
c     parameters: nh,oh,fh the maximum distances for an X...H bond
c                 dmax: a pair at a distance above this is discarded
c                 logically, dmax could be the largest of nh,oh,fh
c                 but a larger value is used, should the others change
c
C.... Output arguments:
C.... nb is the new, total number of bonds, incl. H-bonds
c     li,lj,rr,fcb updated with H-bonds
c
      type='XX'
      nb0=nb
      do i=1,natom
      do j=i+1,natom
        iani=ian(i)
        ianj=ian(j)
c
        if(iani.ne.1.and.ianj.ne.1)  then
          goto 110
        elseif(iani.eq.1.and.ianj.eq.1) then
          goto 110
        endif
c
        r=SQRT((X(1,i)-X(1,j))**2+(X(2,i)-X(2,j))**2+(X(3,i)-X(3,j))**2)
        if(r.gt.dmax) then
          goto 110
        endif
        if(keyw.eq.'HSUB') then
c          build H-bonds between atoms IN DIFFERENT SUBSETS ONLY
          iset=0
          jset=0
          incr=nsub(1)
          ki=1-incr
          kv=0
          do ks=1,msub
            ki=ki+incr
            kv=kv+incr
            do kss=ki,kv
              if(i.eq.indsub(kss)) then
                iset=ks
              elseif(j.eq.indsub(kss)) then
                jset=ks
              endif
            enddo
            if(iset.gt.0.and.jset.gt.0) then
c             subsets identified for both i and j
              goto 120
            endif
          enddo
 120      continue
          if(iset.eq.jset) then
            goto 110
          endif
        endif
c
        if(iani.eq.1) then
c          first, exclude hydrogen not attached to a donor (N,O,F)
          neigb=ibond1(1,i)
          id=ian(neigb)
          if(id.ne.7.and.id.ne.8.and.id.ne.9) then
            goto 110
          endif
c           OK, hydrogen i sits on a donor, now find its partner
          if(ianj.eq.7) then
            type='NH'
          elseif(ianj.eq.8) then
            type='OH'
          elseif(ianj.eq.9) then
            type='FH'
          else
            goto 110
          endif
        else
c              atom j is a hydrogen at this place,
c              repeat the procedures done above for i
c
          neigb=ibond1(1,j)
          id=ian(neigb)
          if(id.ne.7.and.id.ne.8.and.id.ne.9) then
            goto 110
          endif
          if(iani.eq.7) then
            type='NH'
          elseif(iani.eq.8) then
            type='OH'
          elseif(iani.eq.9) then
            type='FH'
          else
            goto 110
          endif
        endif
c      exclude existing bonds:
        do kb=1,nb0
          if((li(kb).eq.i.and.lj(kb).eq.j).or.
     1       (li(kb).eq.j.and.lj(kb).eq.i)) then
            goto 110
           endif
        enddo
c       R=SQRT..   was calctd at the beginning
c
        if((type.eq.'NH'.and.r.lt.nh).or.(type.eq.'OH'.and.r.lt.oh)
     1 .or.(type.eq.'FH'.and.r.lt.fh)) then
c
          nb=nb+1
          LI(nb)=I
          LJ(nb)=J
          RR(nb)=R
          ncg=ncg+1
          bsymb='..'
          fcb(nb)=2.0
          fb(nb)=0.0
          write (ipun,500) ncg,i,j
          write(iptyp,510) ncg, isymb(i),i,bsymb,isymb(j),j, '****'
          write(*,520) i,j
 520   format(/,'H-bond,  i=',i4,'  j=',i4)
        endif
 110   continue
       enddo
       enddo
 500   format(i4,2x,'K',13x,'STRE',6x,i3,'.',6x,i3,'.')
 510   format(i4,'.',4x,a2,i2,1x,a2,1x,a2,i2,3x,'STRE  ','sc=',a4)
c
c      endof BONDH
      end
