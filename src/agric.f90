module AGRIC

  !-------------------------------------------------------+
  ! AUTOMATIC GENERATION OF REDUNDANT INTERNAL COORDINATE |
  !-------------------------------------------------------+

  !use acmatrix
  use system_type
  use transf_type

  implicit none

  private

  public:: gric

  contains

!====================================================================!

      subroutine gric(molecule, trns)
    
      ! This subroutine generats the set of redundant internal coordinates.

      type(molecule_t), intent(inout) :: molecule
      type(transf_t), intent(in) :: trns

      integer nstre, nbend, ndih, nwag, nintc, numat, nat1, nat2
      integer i, j, k, l, ib, ic, il, alloc_err, ncat
      integer, allocatable :: linki(:), linkj(:), link(:), atlist(:)
      integer, allocatable :: acm(:,:)
      integer, allocatable :: list(:)
    
      ! Determining number of internal coordinates to be used.
      ! We consider the stretching of bonded atoms, the bendings between
      ! bonds and the torsions about bonds.
    
      numat = molecule%structure%numat
      allocate(acm(1:numat,1:numat))
      acm = molecule%acm
    
      nstre = 0; nbend = 0;  ndih = 0 ; nwag = 0
      do i = 1, molecule%structure%numat
        nstre = nstre + sum(acm(i,i+1:numat))
        if (acm(i,i) >= 2) nbend = nbend + acm(i,i)*(acm(i,i)-1)/2
        if (trns%lopb) then
            if (acm(i,i) >= 3) nwag = nwag + acm(i,i)
        end if
        if (trns%ldih) then
            do j = i+1, molecule%structure%numat
              if (acm(i,j) == 1) then 
                ndih = ndih + (acm(i,i)-1)*(acm(j,j)-1)
                ! The number calculated above is the maximum number of dihedrals.
                ! Then we must subtract to this number the number of dihedrals
                ! corresonding to the same end atom (which happens in cyclic systems)
                !
                ! Check wheter atom i and j have a common bonded atom
                ncat = 0
                do k = 1, molecule%structure%numat
                  if ((k /= i) .and. (k /= j)) then
                    if ((acm(i,k) == 1) .and. (acm(j,k) == 1)) ncat = ncat + 1
                  end if
                end do
                ndih = ndih - ncat
              end if
            end do
        end if
      end do
    
      nintc = nstre + nbend + ndih + nwag
    
      molecule%intcoord%nintc = nintc
      allocate(molecule%intcoord%coord(1:nintc))

      !---------------------------------------------------------+
      ! IC: Internal coordinate counter.                        |
      ! At the end of the subroutine it must be equal to nintc. |
      !---------------------------------------------------------+
      ic = 0
    
      !-----------------+
      ! Set stretchings |
      !-----------------+
      do i = 1, molecule%structure%numat
        do j = i+1, molecule%structure%numat
          if (acm(i,j) == 1) then
            ic = ic + 1
            molecule%intcoord%coord(ic)%type = "s"
            allocate(molecule%intcoord%coord(ic)%list(1:2))
            allocate(molecule%intcoord%coord(ic)%bmat(1:2,1:3),stat=alloc_err)
            molecule%intcoord%coord(ic)%list = (/i, j/)
            molecule%intcoord%coord(ic)%bmat(1:2,1:3) = zero
          end if 
        end do
      end do 
    
      !--------------+
      ! Set bendings |
      !--------------+
      do i = 1, molecule%structure%numat
        if (acm(i,i) >= 2) then
          allocate(link(1:acm(i,i)))
    
          il = 0
          do ib = 1, molecule%structure%numat
            if (ib == i) cycle
            if (acm(i,ib) == 1) then
              il = il + 1
              link(il) = ib
            end if
          end do
    
          !print *, 'link ', link
          allocate(atlist(1:acm(i,i)))
          atlist = 0
          do
            call combi(atlist,acm(i,i),2)
            if (atlist(1) == 0) exit 
            ic = ic + 1
            molecule%intcoord%coord(ic)%type = "b"
            allocate(molecule%intcoord%coord(ic)%list(1:3))
            molecule%intcoord%coord(ic)%list = (/link(atlist(1)), i, link(atlist(2))/)
            allocate(molecule%intcoord%coord(ic)%bmat(1:3,1:3),stat=alloc_err)
            molecule%intcoord%coord(ic)%bmat(1:3,1:3) = zero
            !print *, 'b ', molecule%intcoord%coord(ic)%list
          end do
          deallocate(link)
          deallocate(atlist) ! Crea un invalid pointer alla prima chiamata, ma non capisco perche'!
        end if
      end do
        
      !---------------+
      ! Set dihedrals |
      !---------------+
      if (trns%ldih) then

          do i = 1, molecule%structure%numat
         
            do j = i+1, molecule%structure%numat
              ! If i and j are bond.
              if (acm(i,j) == 1) then
                ! If i and j have at least another bond, with a different atom.
                if ((acm(i,i) >= 2) .and. (acm(j,j) >= 2)) then
         
                  nat1 = (acm(i,i)-1)
                  nat2 = (acm(j,j)-1)
         
                  allocate(linki(1:nat1),linkj(1:nat2))
         
                  il = 0 ! i1 can at most be equal to nat1
                  do ib = 1, molecule%structure%numat
                    if (ib == i .or. ib == j) cycle ! Skip atoms i and j
                    if (acm(i,ib) == 1) then
                      il = il + 1
                      linki(il) = ib 
                    end if
                  end do
                 
                  il = 0
                  do ib = 1, molecule%structure%numat
                    if (ib == i .or. ib == j) cycle
                    if (acm(j,ib) == 1) then
                      il = il + 1
                      linkj(il) = ib
                    end if
                  end do
         
                  do k = 1, nat1
                    do l = 1, nat2
                      if (linki(k) /= linkj(l)) then ! the two end atoms must not be the same atom!
                        ic = ic + 1
                        molecule%intcoord%coord(ic)%type = "d"
                        allocate(molecule%intcoord%coord(ic)%list(1:4))
                        molecule%intcoord%coord(ic)%list = (/linki(k), i, j, linkj(l)/)
                        allocate(molecule%intcoord%coord(ic)%bmat(1:4,1:3),stat=alloc_err)
                        molecule%intcoord%coord(ic)%bmat(1:4,1:3) = zero
                        !print *, 'd ', molecule%intcoord%coord(ic)%list
                      end if
                    end do
                  end do
         
                  deallocate(linki,linkj)
         
                end if
              end if
         
            end do
          end do
      end if
    
      
      if (trns%lopb) then
          !---------------------------+
          ! Set out of plane bendings |
          !---------------------------+
          do i = 1, molecule%structure%numat
            ! Find the i-th atom with more then three bonds
            if (acm(i,i) >= 3) then
              allocate(link(1:acm(i,i)))
         
              il = 0
              do ib = 1, molecule%structure%numat
                if (ib == i) cycle
                if (acm(i,ib) == 1) then
                  il = il + 1
                  link(il) = ib
                end if
              end do
         
              allocate(atlist(1:acm(i,i)))
              atlist = 0
              ! We add the out pf plane for each of three atoms.
              do l = 1, 3
                if (l == 1) atlist = (/1, 2, 3/)
                if (l == 2) atlist = (/2, 1, 3/)
                if (l == 3) atlist = (/3, 2, 1/)
                ic = ic + 1
                molecule%intcoord%coord(ic)%type = "w"
                allocate(molecule%intcoord%coord(ic)%list(1:4))
                ! atom i is in list(4) and is the anchor atom. See the call to
                ! opb_wilson in the sunroutine BMAT.
                molecule%intcoord%coord(ic)%list = (/link(atlist(1)), link(atlist(2)), link(atlist(3)), i/)
                allocate(molecule%intcoord%coord(ic)%bmat(1:4,1:3),stat=alloc_err)
                molecule%intcoord%coord(ic)%bmat(1:4,1:3) = zero
              end do
              deallocate(link,atlist)
            end if
          end do
      end if

      return
      end subroutine gric

!======================================================================

      subroutine combi(ia,n,j)

      integer, intent(inout) :: ia(1:n)
      integer, intent(in) :: j, n
      integer i, i1

      if(n .le. 0 .or. j .le. 0) return

      if(j .gt. n) then
        write(fout,103) j,n
      else if(ia(1) .eq. 0) then
        do i = 1,j
          ia(i)=i
        end do
        if(j<n) ia(j+1)=0
      else
        do i1 = 1,n
          i=i1
          if (i < n) then
            if(ia(i+1) .ne. ia(i)+1) exit
          end if
          ia(i)=i
        end do
        ia(i)=ia(i)+1
        if(ia(j) .eq. n+1) ia(1)=0
      endif

      return

  103 format('j = ',i5,' > ',i5,' = n is not permitted')
      end subroutine combi


end module AGRIC

