      program driveINTC
cversion 5.0
c  this is just a driver program
c  for all notes see the main subroutine INTC
c
c  version log:
c  09-nov-01 filename added to unit 17 open statement. -rls
c
      implicit real*8 (a-h,o-z)
      parameter (maxat=200,maxelm=18)
      parameter(nkey=8)
c      maxat- maximum number of atoms;maxelm- first 18 elements handled
c       note: this is just a driver program,
c       basic dimensions are set in subr intc
      dimension xyz(3,maxat),ian(maxat),xian(maxat)
      character*100 line
      character*80 title
      character*2 name,elem(maxelm)
      character*4 keyw(nkey),code
      dimension x(maxat),y(maxat),z(maxat),name(maxat)
      data elem /'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     1 'NA','MG','AL','SI','P ','S ','CL','AR'/
      data keyw /'B631','B421','NOFD','BLOW','RADI',
     1 'BOND','HALL','HSUB'/
c      these are the keywords that can be read by INTC (not this
c      driver main pr.) as options;they are needed here only to
c      recognize the end of coordinates.
c FG   ipun will be the punch file for the coordinates generated
c      iptyp will contain a description of the coordinates generated
c matruc 2007/01/17:
c most f90 compilers do not allow to REWIND or BACKSPACE STDIN (unit=5)
c intcin is now read directly
      inp=5
      !open(unit=inp,file='intcin',status='unknown',form='formatted')
      ipun=7
c PSZ    open(unit=7,status='unknown')
      open(unit=ipun,file='intcfl',status='unknown')
      iptyp=17
      open(unit=iptyp,file='icoordtyp', status='unknown')
c
      natom=0
      read(inp,*) natom
      read(inp,*) ! read a blak line
      do i=1, natom
          read(inp,'(a100)',end=1000) line
          read(line,*) name(i),xian(i),x(i),y(i),z(i)
          ian(i)=xian(i)+0.1
      end do

      read(inp,'(a4)',end=1000) code
      do k=1,nkey
c  there may be option cards (for INTC) after the cartesians, these
c  indicate then the end of input for driver main
        if(code.eq.keyw(k)) then
          backspace (inp)
          exit
        endif
      end do 
      !backspace (inp)

      !do 70 i=1,natom
      !  write (ipun,400) name(i),ian(i),x(i),y(i),z(i)
!  70  !continue
! 400  format('N=',a2,6x,i2,'.',7x,3F10.3)

1000  continue
      do 80 i=1,natom
        xyz(1,i)=x(i)
        xyz(2,i)=y(i)
        xyz(3,i)=z(i)
  80  continue

      call intc(xyz,ian,natom,inp,ipun,iptyp,ncg)
      call bummer('normal termination',0,3)

      STOP 'end of intc'
      end
c
