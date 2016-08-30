module intc

  ! Le function stretching(), bending() etc etc restituiscono una matrice B avente dimensioni B(1:N,1:3) dove
  ! N = numero di atomi coinvolti nella definizione della coordinata.

  use parameters
  use iomatrix
  use system_type
!  use sysop, only : get_bmat
  use vecmath

  implicit none

  private

  real(kind=dp), allocatable :: u(:,:), b(:,:), ub(:,:)
  integer, allocatable :: ibas(:,:), mult(:)
  integer line(12), ich, itors
  type(molecule_t), pointer :: molecule

  public :: bmat, nat_bmat, set_equilibrium_intc, set_equilibrium_natint

  contains

!==============================================================================

      subroutine bmat(molec)

      !+------------------------------------------------------------------------+
      ! Calculates transformation matrix B in the matrix equation               !
      !                            r = B * x                                    !
      ! where r and x are column vectors of internal and  cartesian coordinates !
      ! respectively.                                                           !
      !+------------------------------------------------------------------------+

      implicit none

      real(kind=dp) e1(3),e2(3),e3(3),en(3)
      real(kind=dp) bs(2,3), bb(3,3), bw(4,3), bl(3,3), bd(4,3), db(12,3)
      integer i, ii, j, nta, np, nq, la !,list(4)
      integer, pointer :: list(:)

      type(molecule_t), target, intent(in) :: molec
!     TODO: This must be changed since now I have the module ACM.
!
!     initialize array elements to zero
!
!     ibas    is an array in which is stored information concerning
!             the connectivity of the molecule
!     mult    is an array containing the multiplicity of each atom
!
      molecule => molec

      allocate (ibas(1:6,1:molecule%structure%numat))
      allocate (mult(1:molecule%structure%numat))

      ibas = 0
      mult = 0
!      write(fout,2000)
!
!     specify a "type" of coordinate involving the atoms in "list"
!
!     types available are :
!
!     s     stretch
!     b     bend
!     w     wag (out-of-plane bend)
!     l     linear bend
!     t     torsion
!     d     dihedral angle
!
INTC: do ii = 1, molecule%intcoord%nintc
!
!       branch according to the type of coordinate in order to
!       compute a row of the b matrix

        list => molecule%intcoord%coord(ii)%list

        select case(molecule%intcoord%coord(ii)%type)

          case("s","S","stre","STRE")
         
            !+----------------------------------------------------+ 
            !   stretching coordinate for bond p-q specified by :
            !   list(1) = p , list(2) = q     
            !+----------------------------------------------------+ 
            !write(iw,2100) ii,lab1(ii),lab2(ii),lmnt(np),np,lmnt(nq),nq,sf(ii)
            !forall (i=1:2) list(i) = molecule%intcoord%coord(ii)%list(i) 
            molecule%intcoord%coord(ii)%bmat = stretching(list(1),list(2))
!      
!         construct a topological map IBAS of the bonded atoms
!         and set the multiplicity MULT of each atom
!         IBAS is used in the function torsion().
            if (any(mult(list(1:2)) /= list(1:2))) then
            !if (any(mult(list(1:2)) == 0 )) then
              mult(list(1:2)) = mult(list(1:2)) + 1
              ibas(mult(list(1)),list(1)) = list(2)
              ibas(mult(list(2)),list(2)) = list(1)
            !end if
            end if

          case("b","B","bend","BEND")
         
            !+---------------------------------------------------
            ! bending coordinate for angle p-q-r specified by :
            !    list(1) = p , list(2) = q , list(3) = r
            ! (see : wilson,decius and cross, page 56 )
            !+---------------------------------------------------
         
            !write(iw,2200) ii,lab1(ii),lab2(ii),lmnt(np),np,lmnt(nq),nq,lmnt(nr),nr,sf(ii)
            !forall (i=1:3) list(i) = molecule%intcoord%coord(ii)%list(i) 
            molecule%intcoord%coord(ii)%bmat = bending(list(1),list(2),list(3)) 
         
          case("w","W","outp","OUTP")
         
            !+---------------------------------------------------------------
            ! out-of-plane bending coordinate for four coplanar atoms,
            ! apical atom q bonded to end atom p and anchor atoms r1 and r2,
            ! p is at the center.
            ! specified by: list(1)=p ,list(2)=r1 ,list(3)=r2 ,list(4)=q
            ! (see : mcintosh et al., can.j.chem. 56 (1978) 1289 )
            !+---------------------------------------------------------------

            !write(iw,2300) ii,lab1(ii),lab2(ii),lmnt(np),np,lmnt(nq),nq,lmnt(nr1), &
            !               nr1,lmnt(nq),nq,lmnt(nr2),nr2,sf(ii)
            !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
            !molecule%intcoord%coord(ii)%bmat = wag(list(1),list(2),list(3),list(4))
            !molecule%intcoord%coord(ii)%bmat = opb_wilson(list(4),list(3),list(2),list(1))
            molecule%intcoord%coord(ii)%bmat = opb_wilson(list(1),list(2),list(3),list(4))
         
          case("l","L","linb","LINB")
         
            !+----------------------------------------------------------------
            ! linear bending coordinate for collinear atoms p-q-r in plane of
            ! dummy atom d specified by :
            ! list(1) = p , list(2) = q , list(3) = r , list(4) = d
            !+----------------------------------------------------------------
            !write(iw,2400) ii,lab1(ii),lab2(ii),lmnt(np),np,lmnt(nq), & 
            !               nq,lmnt(nr),nr,lmnt(nd),nd,sf(ii)
            !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
            molecule%intcoord%coord(ii)%bmat = bendli(list(1),list(2),list(3),list(4)) 
         
         
          case("linc")
         
            !+----------------------------------------------------------------
            ! linear bending coordinate for collinear atoms p-q-r in plane of
            ! dummy atom d specified by :
            ! list(1) = p , list(2) = q , list(3) = r , list(4) = d
            !+----------------------------------------------------------------
            !write(iw,2400) ii,lab1(ii),lab2(ii),lmnt(np),np,lmnt(nq), & 
            !               nq,lmnt(nr),nr,lmnt(nd),nd,sf(ii)
            !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
            molecule%intcoord%coord(ii)%bmat = bendlic(list(1),list(2),list(3),list(4)) 
            print *, 'LINC ', molecule%intcoord%coord(ii)%bmat
          case("linp")
         
            !+----------------------------------------------------------------
            ! linear bending coordinate for collinear atoms p-q-r in plane of
            ! dummy atom d specified by :
            ! list(1) = p , list(2) = q , list(3) = r , list(4) = d
            !+----------------------------------------------------------------
            !write(iw,2400) ii,lab1(ii),lab2(ii),lmnt(np),np,lmnt(nq), & 
            !               nq,lmnt(nr),nr,lmnt(nd),nd,sf(ii)
            !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
            molecule%intcoord%coord(ii)%bmat = bendlip(list(1),list(2),list(3),list(4)) 
            print *, 'LINP ', molecule%intcoord%coord(ii)%bmat
 
          case("d","D")
         
            !+---------------------------------------------------
            ! dihedral angle p-q-r-s specified by:
            ! list(1) = p, list(2) = q, list(3) = r, list(4) = s
            ! (see: wilson, decius and cross, pp. 60-61.
            !+---------------------------------------------------

            !write(iw,2600) ii,lab1(ii),lab2(ii),lmnt(np),np,lmnt(nq),nq,lmnt(nr), &
            !               nr,lmnt(ns),ns,sf(ii)
            !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
            molecule%intcoord%coord(ii)%bmat = dihedral(list(1),list(2),list(3),list(4)) 
            !forall (i=1:4) 
            !  molecule%intcoord%coord(ii)%bmat(3*(list(i)-1)+1:3*list(i)) = bd(i,:)
            !end forall

          case("t","T")
			! itors can be accessed in the whole module	
		    itors = ii
            !+---------------------------------------------------------
            ! torsional internal coordinate about bond j-k
            ! specified by :  list(1) = j  ,  list(2) = k
            ! (see : williams et al., j. mol. struct., 55 (1979) 147 )
            !+---------------------------------------------------------

            !write(iw,2500) ii,lab1(ii),lab2(ii),lmnt(list(1)),list(1),lmnt(list(2)),list(2),sf(ii)
            !forall (i=1:2) list(i) = molecule%intcoord%coord(ii)%list(i) 
            db = torsion(list(1),list(2))
            ! NTA e' il numero di atomi coinvolti nella coordinata di torsione.
            nta = count(line > 0)
            !--------------------------------------------------------------!
            ! Per la torsione e' necessario prima capire quanti atomi      !
            ! sono coinvolti (max 12) quindi si puo' allocare l'array bmat !
            !--------------------------------------------------------------!
            deallocate(molecule%intcoord%coord(ii)%list)
            allocate(molecule%intcoord%coord(ii)%list(1:nta))
            molecule%intcoord%coord(ii)%list(1:nta) = line(1:nta)
            allocate(molecule%intcoord%coord(ii)%bmat(1:nta,1:3))
            forall (i=1:nta) 
              molecule%intcoord%coord(ii)%bmat(i,1:3) = db(i,1:3)
            end forall
            ! Normalizzazione : ICH e' il numero di catene di diedri
            ! coinvolti nella torsione.
            molecule%intcoord%coord(ii)%bmat = molecule%intcoord%coord(ii)%bmat / ich
            molecule%intcoord%coord(itors)%val = molecule%intcoord%coord(itors)%val / ich

        end select

      end do INTC
          
      !allocate(b(1:molecule%intcoord%nintc,1:3*molecule%structure%numat))
      !b = get_bmat(molecule)
      !write(fout,2004)
      !call layout(b,molecule%intcoord%nintc,3*molecule%structure%numat, &
      !            molecule%intcoord%nintc,3*molecule%structure%numat)
      !deallocate(b)
      deallocate(ibas,mult)

      include 'formats'
      return
      end subroutine bmat

!==============================================================================

      subroutine nat_bmat(molec)

      !+------------------------------------------------------------------------+
      ! Calculates transformation matrix N * B in the matrix equation           !
      !                            r = N * B * x  = R * x                       !
      ! where r and x are column vectors of internal and  cartesian coordinates !
      ! respectively.                                                           !
      !+------------------------------------------------------------------------+

      implicit none

      real(kind=dp) e1(3),e2(3),e3(3),en(3)
      real(kind=dp) bs(2,3), bb(3,3), bw(4,3), bl(3,3), bd(4,3), db(12,3)
      integer i, ii, j, ij, nta, np, nq, la, natc
      integer, pointer :: list(:)

      type(molecule_t), target, intent(in) :: molec
!     TODO: This must be changed since now I have the module ACM.
!
      molecule => molec

!      write(fout,2000)
!
!     specify a "type" of coordinate involving the atoms in "list"
!
!     types available are :
!
!     s     stretch
!     b     bend
!     w     wag (out-of-plane bend)
!     l     linear bend
!     d     dihedral angle
!
      natc = size(molecule%intcoord%ncoord)

INTC: do ii = 1, natc
!
!       branch according to the type of coordinate in order to
!       compute a row of the b matrix

INTC2:  do ij = 1, size(molecule%intcoord%ncoord(ii)%coord)

          list => molecule%intcoord%ncoord(ii)%coord(ij)%list
    !print *, 'TYPE ', molecule%intcoord%ncoord(ii)%coord(ij)%type
          select case(molecule%intcoord%ncoord(ii)%coord(ij)%type)
            case("s","S","stre","STRE")
              !+----------------------------------------------------+ 
              !   stretching coordinate for bond p-q specified by :
              !   list(1) = p , list(2) = q     
              !+----------------------------------------------------+ 
              !forall (i=1:2) list(i) = molecule%intcoord%coord(ii)%list(i) 
              molecule%intcoord%ncoord(ii)%coord(ij)%bmat = stretching(list(1),list(2))
            case("b","B","bend","BEND")
              !+---------------------------------------------------
              ! bending coordinate for angle p-q-r specified by :
              !    list(1) = p , list(2) = q , list(3) = r
              ! (see : wilson,decius and cross, page 56 )
              !+---------------------------------------------------
              !forall (i=1:3) list(i) = molecule%intcoord%coord(ii)%list(i) 
              molecule%intcoord%ncoord(ii)%coord(ij)%bmat = bending(list(1),list(2),list(3)) 
            case("w","W","outp","OUTP")
              !+---------------------------------------------------------------
              ! out-of-plane bending coordinate for four coplanar atoms,
              ! apical atom q bonded to end atom p and anchor atoms r1 and r2,
              ! specified by: list(1)=p ,list(2)=r1 ,list(3)=r2 ,list(4)=q
              ! (see : mcintosh et al., can.j.chem. 56 (1978) 1289 )
              !+---------------------------------------------------------------
              !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
              !molecule%intcoord%coord(ii)%bmat = wag(list(1),list(2),list(3),list(4))
              !molecule%intcoord%ncoord(ii)%coord(ij)%bmat = opb_wilson(list(4),list(3),list(2),list(1))
              molecule%intcoord%ncoord(ii)%coord(ij)%bmat = opb_wilson(list(1),list(2),list(3),list(4))
            case("l","L","linb","LINB")
              !+----------------------------------------------------------------
              ! linear bending coordinate for collinear atoms p-q-r in plane of
              ! dummy atom d specified by :
              ! list(1) = p , list(2) = q , list(3) = r , list(4) = d
              !+----------------------------------------------------------------
              !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
              molecule%intcoord%ncoord(ii)%coord(ij)%bmat = bendli(list(1),list(2),list(3),list(4)) 
              print *, 'lin ', molecule%intcoord%ncoord(ii)%coord(ij)%bmat
            case("linc")
              !+----------------------------------------------------------------
              ! linear bending coordinate for collinear atoms p-q-r in plane of
              ! dummy atom d specified by :
              ! list(1) = p , list(2) = q , list(3) = r , list(4) = d
              !+----------------------------------------------------------------
              !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
              molecule%intcoord%ncoord(ii)%coord(ij)%bmat = bendlic(list(1),list(2),list(3),list(4)) 
              ! print *, 'linc ', molecule%intcoord%ncoord(ii)%coord(ij)%bmat
            case("linp")
              !+----------------------------------------------------------------
              ! linear bending coordinate for collinear atoms p-q-r in plane of
              ! dummy atom d specified by :
              ! list(1) = p , list(2) = q , list(3) = r , list(4) = d
              !+----------------------------------------------------------------
              !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
              molecule%intcoord%ncoord(ii)%coord(ij)%bmat = bendlip(list(1),list(2),list(3),list(4)) 
              !print *, 'linp ', molecule%intcoord%ncoord(ii)%coord(ij)%bmat
            case("d","D","dih")
              !+---------------------------------------------------
              ! dihedral angle p-q-r-s specified by:
              ! list(1) = p, list(2) = q, list(3) = r, list(4) = s
              ! (see: wilson, decius and cross, pp. 60-61.)
              !+---------------------------------------------------
              !forall (i=1:4) list(i) = molecule%intcoord%coord(ii)%list(i) 
              molecule%intcoord%ncoord(ii)%coord(ij)%bmat = dihedral(list(1),list(2),list(3),list(4)) 

            end select

        end do INTC2
      end do INTC
          
      include 'formats'
      return
      end subroutine nat_bmat

!==============================================================================

      function stretching(np,nq) result(bs)

      real(kind=dp) bs(2,3)
      integer, intent(in) :: np, nq
      integer lp,lq
      real(kind=dp) e1(1:3), r1

      bs = zero

      call dc(np,nq,e1,r1)
      bs(1,1:3)= -e1(1:3)
      bs(2,1:3)=  e1(1:3)

      return
      end function stretching

!==============================================================================

      function bending(np,nq,nr) result(bb)

      real(kind=dp) bb(3,3)
      integer, intent(in) :: np, nq, nr
      real(kind=dp) e1(1:3), e2(1:3), qp, qr, br(1:3), bp(1:3), c, s

      bb = zero
      
      call dc(nq,np,e1,qp)
      call dc(nq,nr,e2,qr)
      c=dot_product(e1,e2)
      s = dsqrt( one - c*c )

      bp = ( c*e1 - e2 ) / (qp*s)
      br = ( c*e2 - e1 ) / (qr*s)

      bb(1,:) = bp
      bb(2,:) = -(bp + br)
      bb(3,:) = br

      return
      end function bending

!==============================================================================

      function bendli(np,nq,nr,nd) result(bl)

      real(kind=dp) bl(3,3)
      integer, intent(in) :: np, nq, nr, nd
      real(kind=dp) qp, qr, rd, s, r, e1(3), e2(3), e3(3), en(3), bp(3), br(3)

      bl = zero

      call dc(nq,np,e1,qp) ! calcolo qp, e1 non mi serve
      call dc(nq,nr,e2,qr)
      call dc(nr,nd,e3,rd)
      en = cross(e2,e3)
      en = en/norm(en)
      e1 = cross(en,e2)
      e1 = e1/norm(e1)

      bl(1,:) = -(e1+en)/qp
      bl(3,:) = -(e1+en)/qr
      bl(2,:) = -(bl(1,:)+bl(3,:))

      end function bendli


!==============================================================================

      function bendlic(np,nq,nr,nd) result(bl)

      real(kind=dp) bl(3,3)
      integer, intent(in) :: np, nq, nr, nd
      real(kind=dp) qp, qr, rd, s, r, e1(3), e2(3), e3(3), en(3), bp(3), br(3)

      bl = zero

      call dc(nq,np,e1,qp) ! calcolo qp, e1 non mi serve
      call dc(nq,nr,e2,qr)
      call dc(nr,nd,e3,rd)
      en = cross(e2,e3)
      en = en/norm(en)
      e1 = cross(en,e2)
      e1 = e1/norm(e1)

      bl(1,:) = -(e1)/qp
      bl(3,:) = -(e1)/qr
      bl(2,:) = -(bl(1,:)+bl(3,:))

      end function bendlic


!==============================================================================

      function bendlip(np,nq,nr,nd) result(bl)

      real(kind=dp) bl(3,3)
      integer, intent(in) :: np, nq, nr, nd
      real(kind=dp) qp, qr, rd, s, r, e1(3), e2(3), e3(3), en(3), bp(3), br(3)

      bl = zero

      call dc(nq,np,e1,qp) ! calcolo qp, e1 non mi serve
      call dc(nq,nr,e2,qr)
      call dc(nr,nd,e3,rd)
      en = cross(e2,e3)
      en = en/norm(en)
      e1 = cross(en,e2)
      e1 = e1/norm(e1)

      bl(1,:) = -(en)/qp
      bl(3,:) = -(en)/qr
      bl(2,:) = -(bl(1,:)+bl(3,:))

      end function bendlip

!==============================================================================

      function dihedral(np,nq,nr,ns) result(bd)

      integer, intent(in) :: np, nq, nr, ns
      real(kind=dp) bd(4,3)
      real(kind=dp) r1, r2, r3, ca, cb, sa, sb
      real(kind=dp) e1(3), e2(3), e3(3), en(3)

      bd = zero

      call dc(np,nq,e1,r1)
      call dc(nq,nr,e2,r2)
      call dc(nr,ns,e3,r3)
      en = cross(e1,e2)
      sa = norm(en)
      en = en/sa
      ca = -dot_product(e1,e2)
      e1 = cross(e3,e2)
      sb = norm(e1)
      e1 = e1 /sb
      cb = -dot_product(e2,e3)
!      if( sa*sa .lt. tm4 ) go to 999
!      if( sb*sb .lt. tm4 ) go to 999

      bd(1,:) = -en / (r1*sa)
      bd(2,:) = (en*(r2-r1*ca))/(r1*r2*sa) + (cb*e1)/(r2*sb)  
      bd(4,:) = -e1 / (r3*sb) 
      bd(3,:) = -(bd(1,:)+bd(2,:)+bd(4,:))

      end function dihedral

!==============================================================================

      function wag(np,nr1,nr2,nq) result(bw)

      integer, intent(in) :: np, nr1, nr2, nq
      real(kind=dp) bw(4,3)
      real(kind=dp) pq, r1, r2, sp, c1, c2, s1, s2
      real(kind=dp) en(3), e1(3), e2(3), e3(3)
      integer ia

      bw = zero

      call dc(nq,np,e3,pq)
      call dc(nq,nr1,e1,r1)
      call dc(nq,nr2,e2,r2)
      en = cross(e1,e2)
      sp = norm(en)
      en = en / sp
      c1 = dot_product(e3,e2)
      c2 = dot_product(e3,e1)
      s1 = sqrt(one - c1*c1)
      s2 = sqrt(one - c2*c2)

      bw(1,:) = en / pq
      bw(2,:) = en*s1 / (r1*sp)
      bw(3,:) = en*s2 / (r2*sp)
      bw(4,:) = -(bw(1,:)+bw(2,:)+bw(3,:))
       
      end function wag

!==============================================================================

      function opb_wilson(np,nr1,nr2,nq) result(bw)

      !                                      r1                2   
      !                                     /                 /
      ! Out of plane bending of p/1 :  p - q      o      1 - 4
      !                                     \                 \
      !                                      r2                3
      integer, intent(in) :: np, nr1, nr2, nq
      real(kind=dp) bw(4,3)
      real(kind=dp) pq, r1, r2, sp, c1, c2, s1, s2
      real(kind=dp) e1(3), e2(3), e3(3)
      real(kind=dp) e1x2(3), sTheta, cTheta, tTheta, cTs1
      integer ia

      !print *, 'OPB :', np, nr1, nr2, nq
      bw = zero

      ! According to Wilson. D. C. pag. 59
      ! e3 = e41; e1 = e42; e2 = e43
      call dc(nq,np,e3,pq)
      call dc(nq,nr1,e1,r1)
      call dc(nq,nr2,e2,r2)

      c1 = dot_product(e1,e2)
      s1 = sqrt(one - c1**2)

      e1x2 = cross(e1,e2)

      sTheta = dot_product(e1x2,e3)/s1
      cTheta = sqrt(one-sTheta**2)
      tTheta = sTheta/cTheta
      cTs1 = cTheta*s1 

      bw(1,:) = (e1x2/(cTs1) - tTheta*e3) / pq
      bw(2,:) = (cross(e2,e3)/(cTs1) - (tTheta/s1**2)*(e1-c1*e2)) / r1
      bw(3,:) = (cross(e3,e1)/(cTs1) - (tTheta/s1**2)*(e2-c1*e1)) / r2
      bw(4,:) = -(bw(1,:)+bw(2,:)+bw(3,:))
       
      end function opb_wilson


!==============================================================================

      function torsion(ia,ib) result(db)

      real(kind=dp) db(12,3)
      integer, intent(in) :: ia, ib
      real(kind=dp) dels1, dels2, dels3, dels4, sa, sb, r1, r2, r3, ca, cb
      real(kind=dp) e1(3), e2(3), e3(3), en(3)
      integer i, j, k, l, ita, itb, ina, la, lb, np, nq
      real(kind=dp) b1(3), b2(3), b3(3), b4(3), n1(3), n2(3), phi

      line = 0
      db = zero
!
!     ia,ib  are the 'central' atoms defining the bond about which
!     torsion occurs
!
      line(1)=ia
      line(2)=ib
      ich=0
      ina=3

!     select two 'terminal' atoms, np and nq
ACH:  do la=1,mult(ia)

BCH:    do lb=1,mult(ib)

          np=ibas(la,ia)
          nq=ibas(lb,ib)
!        
!         check that terminal atoms are not central atoms
!        
          if(np == ib .or. nq == ia)  cycle  ! questo controllo non mi convince...rivedere!
!        
!         calc. direction cosines,vector and scalar products of bond vectors
!        
          call dc(np,ia,e1,r1)
          call dc(ia,ib,e2,r2)
          call dc(ib,nq,e3,r3)
          en = cross(e1,e2)
          sa = norm(en)
          ca= -dot_product(e1,e2)
          e1 = cross(e3,e2)
          cb= -dot_product(e2,e3)
          sb = norm(e1)
!        
!         Check for the presence of a collinear subsection :
!         if so then reject this chain
!        
          if((sa*sa <= tm4) .or. (sb*sb <= tm4)) cycle 
!        
!         ich   is the number of four-atom chains not containing
!               collinear subsections
!        
          ich=ich+1
!        
!         Test whether each terminal atom has occurred in a previous chain :
!         if not then add the atom number to list
!        
          do 540 l=3,ina
          ita=l
          if(np.eq.line(l)) go to 541
  540     continue

          line(ina)=np
          ina=ina+1

  541     do 542  l=3,ina
          itb=l
          if(nq.eq.line(l)) go to 550
  542     continue

          line(ina)=nq
          ina=ina+1
!        
!         calculate components of s vector for each atom of chain
!         (see : wilson,decius and cross, page 60) and then
!         add contributions from this chain to b matrix elements obtained
!         from previous chains
!        
  550     do j=1,3
            dels1 =  -en(j) / (r1*sa)
            dels2 =  (en(j)*(r2 - r1*ca)) / (r1*r2*sa) + (cb*e1(j))/(r2*sb)
            dels4 =  -e1(j) / (r3*sb)
            dels3 = -(dels1+dels2+dels4)
            db(1,j) = db(1,j) + dels2
            db(2,j) = db(2,j) + dels3
            db(ita,j) = db(ita,j) + dels1
            db(itb,j) = db(itb,j) + dels4
          end do

		  !----------------------------------------
		  ! Find equilibrium value of the torsion
		  !----------------------------------------
          b1 = molecule%structure%atom(np)%coord(:) - molecule%structure%atom(ia)%coord(:)  
          b1 = b1/norm(b1)
          b2 = molecule%structure%atom(ib)%coord(:) - molecule%structure%atom(ia)%coord(:)
          b2 = b2/norm(b2)
          b3 = molecule%structure%atom(ib)%coord(:) - molecule%structure%atom(nq)%coord(:)
          b3 = b3/norm(b3)
          n1 = cross(b1,b2) ! normale al piano 1 2 3
          n1 = n1/norm(n1)
          n2 = cross(b2,b3) ! normale al piano 2 3 4
          n2 = n2/norm(n2)
          phi = atan2 (dot_product(b2,cross(n1,n2)), dot_product(n1,n2)) 
          ! ATAN2 restituisce phi in [-pi,pi], io lo rimetto in [0,2pi] 
          ! mi serve per definire il verso di percorrenza della rotazione: vedi itransf.f90
          if (phi < (-epsilon(one))) phi = phi + 2*PI
          molecule%intcoord%coord(itors)%val = molecule%intcoord%coord(itors)%val + phi
		  !----------------------------------------

        end do BCH

      end do ACH

      end function torsion

!==============================================================================

      subroutine dc(i,j,e,r)
!
!     calculates direction cosines and magnitude of bond vector i-j
!
!     on exit:
!     e: vector from atom i to atom j
!     r: norm of the vector, i.e. distance between atoms i and j
!
      integer, intent(in) :: i, j
      real(kind=dp), intent(out) :: e(3)
      real(kind=dp), intent(out) :: r

      e = zero
      e = molecule%structure%atom(j)%coord - molecule%structure%atom(i)%coord
      r = sqrt(dot_product(e,e))
      e = e/r

      return
      end subroutine dc

!==============================================================================

    subroutine set_equilibrium_intc(molecule)

    type(molecule_t), intent(in out) :: molecule
    real(kind=dp) intc(1:molecule%intcoord%nintc)
    real(kind=dp) b1(3), b2(3), b3(3), b4(3)
    real(kind=dp) n1(3), n2(3), phi, cosang
    integer i, ierr
    integer, pointer :: list(:)
 
    do i = 1, molecule%intcoord%nintc
      list => molecule%intcoord%coord(i)%list 
      select case(molecule%intcoord%coord(i)%type) 
        case("s","S")
          b1 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(1))%coord(:)  
          molecule%intcoord%coord(i)%val = sqrt(dot_product(b1,b1))
        case("b","B")
          b1 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(1))%coord(:)  
          b1 = b1/sqrt(dot_product(b1,b1)) 
          b2 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(3))%coord(:)
          b2 = b2/sqrt(dot_product(b2,b2)) 
          molecule%intcoord%coord(i)%val = acos(dot_product(b1,b2))
        case("d","D")
          !if (molecule%intcoord%coord(i)%assigned) cycle
          b1 = molecule%structure%atom(list(1))%coord(:) - molecule%structure%atom(list(2))%coord(:)  
          b1 = b1/norm(b1)
          b2 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(2))%coord(:)
          b2 = b2/norm(b2)
          b3 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(4))%coord(:)
          b3 = b3/norm(b3)
          n1 = cross(b1,b2) ! normale al piano 1 2 3
          n1 = n1/norm(n1)
          n2 = cross(b2,b3) ! normale al piano 2 3 4
          n2 = n2/norm(n2)
          phi = atan2 (dot_product(b2,cross(n1,n2)), dot_product(n1,n2)) 
          ! ATAN2 restituisce phi in [-pi,pi], io lo rimetto in [0,2pi] 
          ! mi serve per definire il verso di percorrenza della rotazione: vedi itransf.f90
          if (phi < (-epsilon(one))) phi = phi + 2*PI
          molecule%intcoord%coord(i)%val = phi
        case("l","L")
          molecule%intcoord%coord(i)%val = zero
        case("w","W")
          !if (molecule%intcoord%coord(i)%assigned) cycle
!          write(fout,*)"Out-of-plane coordinate value read from input"
          b1 = molecule%structure%atom(list(1))%coord(:) - molecule%structure%atom(list(4))%coord(:)  
          b1 = b1/norm(b1)
          b2 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(4))%coord(:)
          b2 = b2/norm(b2)
          b3 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(4))%coord(:)  
          b3 = b3/norm(b3)
          n2 = cross(b2,b3) ! normale al piano 2 3 4
          n2 = n2/norm(n2)  ! nu = sin(phi1), see WDC pag 59.
          molecule%intcoord%coord(i)%val = asin(dot_product(b1,n2))  ! Devo controllare se il segno va bene.
        case("t","T")
          write(fout,*)"================================================"
          write(fout,*)"Double check the value of the Torsion coordinate"
          write(fout,*)"================================================"
      end select
    end do
    
    return
    end subroutine set_equilibrium_intc
    
!==============================================================================

    subroutine set_equilibrium_natint(molecule)

    type(molecule_t), intent(in out) :: molecule
    real(kind=dp) intc(1:molecule%intcoord%nintc)
    real(kind=dp) b1(3), b2(3), b3(3), b4(3)
    real(kind=dp) n1(3), n2(3), phi, cosang, rcnorm
    integer i, ij,  ierr
    integer, pointer :: list(:)
 
    !print *, 'in equilibroium natint'
    !print *, 'nintc ', molecule%intcoord%nintc

    do i = 1, molecule%intcoord%nintc
      rcnorm = 1.0
      !if (molecule%intcoord%ncoord(i)%type == "d") then
      !  rcnorm = size(molecule%intcoord%ncoord(i)%coord)
      !else
      !  rcnorm = sqrt(sum(molecule%intcoord%ncoord(i)%coord(:)%c**2))
      !end if
      do ij = 1, size(molecule%intcoord%ncoord(i)%coord)
        list => molecule%intcoord%ncoord(i)%coord(ij)%list 
        select case(molecule%intcoord%ncoord(i)%coord(ij)%type) 
        case("s","S")
          b1 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(1))%coord(:)  
          molecule%intcoord%ncoord(i)%coord(ij)%val = sqrt(dot_product(b1,b1))
        case("b","B")
          b1 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(1))%coord(:)  
          b1 = b1/sqrt(dot_product(b1,b1)) 
          b2 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(3))%coord(:)
          b2 = b2/sqrt(dot_product(b2,b2)) 
          molecule%intcoord%ncoord(i)%coord(ij)%val = acos(dot_product(b1,b2))
        case("d","D")
          !if (molecule%intcoord%coord(i)%assigned) cycle
          b1 = molecule%structure%atom(list(1))%coord(:) - molecule%structure%atom(list(2))%coord(:)  
          b1 = b1/norm(b1)
          b2 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(2))%coord(:)
          b2 = b2/norm(b2)
          b3 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(4))%coord(:)
          b3 = b3/norm(b3)
          n1 = cross(b1,b2) ! normale al piano 1 2 3
          n1 = n1/norm(n1)
          n2 = cross(b2,b3) ! normale al piano 2 3 4
          n2 = n2/norm(n2)
          phi = atan2 (dot_product(b2,cross(n1,n2)), dot_product(n1,n2)) 
          ! ATAN2 restituisce phi in [-pi,pi], io lo rimetto in [0,2pi] 
          ! mi serve per definire il verso di percorrenza della rotazione: vedi itransf.f90
          if (phi < (-epsilon(one))) phi = phi + 2*PI
          molecule%intcoord%ncoord(i)%coord(ij)%val = phi
        case("l","L")
          molecule%intcoord%ncoord(i)%coord(ij)%val = zero
        case("w","W")
          !if (molecule%intcoord%coord(i)%assigned) cycle
!          write(fout,*)"Out-of-plane coordinate value read from input"
          b1 = molecule%structure%atom(list(1))%coord(:) - molecule%structure%atom(list(4))%coord(:)  
          b1 = b1/norm(b1)
          b2 = molecule%structure%atom(list(2))%coord(:) - molecule%structure%atom(list(4))%coord(:)
          b2 = b2/norm(b2)
          b3 = molecule%structure%atom(list(3))%coord(:) - molecule%structure%atom(list(4))%coord(:)  
          b3 = b3/norm(b3)
          n2 = cross(b2,b3) ! normale al piano 2 3 4
          n2 = n2/norm(n2)  ! nu = sin(phi1), see WDC pag 59.
          molecule%intcoord%ncoord(i)%coord(ij)%val = asin(dot_product(b1,n2))  ! Devo controllare se il segno va bene.
        end select
        molecule%intcoord%ncoord(i)%val = molecule%intcoord%ncoord(i)%val + &
        molecule%intcoord%ncoord(i)%coord(ij)%val*molecule%intcoord%ncoord(i)%coord(ij)%c
      end do
      molecule%intcoord%ncoord(i)%val = molecule%intcoord%ncoord(i)%val/rcnorm
    end do
    
    return
    end subroutine set_equilibrium_natint

end module intc
