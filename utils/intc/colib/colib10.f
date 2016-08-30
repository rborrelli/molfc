c deck getchr
      subroutine getchr (prompt,string,default)
c     # reads a character string string from the keyboard with prompt
c     # and default.
c     # if default.eq.' ', there is no default.  reprompt on <cr>.
      character*(*) prompt, default, string
      if (default.ne.' ') then
*@ifdef f90
         write (6,"(1x,a,' <',a,'> ')",advance='no') prompt, default
*@else
*         write (6,100) prompt, default
*@ifdef $format
*100      format (' ',a,' <',a,'> ',$)
*@else
*100      format (' ',a,' <',a,'>')
*@endif
*@endif
         read (5,'(a)') string
c
         if (string.eq.'&Basis') then
            read (5,'(a)') string
            read (5,'(a)') string
         elseif (string.eq.'&End') then
            read (5,'(a)') string
            if (string.eq.'&SO') then
               read (5,'(a)') string
            endif
         endif
         if (string.eq.'&So') then
            read (5,'(a)') string
         endif
c
         write(13,'(a)') string
         if (string.eq.' ') string=default
      else
10       continue
*@ifdef f90
         write (6,'(1x,a,1x)',advance='no') prompt
*@else
*         write (6,110) prompt
*@ifdef $format
*110      format (' ',a,' ',$)
*@else
*110      format (' ',a)
*@endif
*@endif
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') goto 10
      endif
      call strupp (string)
      return
      end
c deck getitle
      subroutine getitle (prompt,title)
c     # reads a character string from the keyboard with prompt,
c     # and defaults to the input value.
      character*(*) prompt, title
      character*130 line
c
      integer  fstrlen
      external fstrlen
c
      if (title.ne.' ') then
         ltitle=fstrlen(title)
*@ifdef f90
         write (6,100,advance='no') prompt, title(1:ltitle)
  100    format (1x,a,' <',a,'>',/,' --> ')
*@else
*         write (6,100) prompt,title(1:ltitle)
*@ifdef $format
*100      format (' ',a,' <',a,'> ',/,' --> ',$)
*@else
*  100    format (' ',a,' <',a,'>')
*@endif
*@endif
      else
*@ifdef f90
         write (6,"(1x,a,/,' --> ')",advance='no') prompt
*@else
*         write (6,110) prompt
*@ifdef $format
*110      format (' ',a,/,' --> ',$)
*@else
*  110    format (' ',a)
*@endif
*@endif
      endif
      read (5,'(a)') line
c
      if ((line.eq.'&Basis').or.(line.eq.'&End').or.
     & (line.eq.'&SO').or.(line.eq.'&Title')) then
         read (5,'(a)') line
      endif
c
      write(13,'(a)') line
      if (line.ne.' ') title=line
      return
      end
c deck getreal
      subroutine getreal ( prompt, anum, default)
c     # reads a real number anum from the keyboard with prompt
c     # and default
c     # default.eq.'$', 1000000 is returned on <cr>
c     # default.eq.' ', there is no default.  reprompt on <cr>.
      implicit real*8 (a-h,o-z)
      character*(*) prompt, default
      character*16 string
9     if (default.eq.' ') then
10       continue
*@ifdef f90
         write (6,110,advance='no') prompt
110      format(1x,a,1x)
*@else
*         write (6,110) prompt
*@ifdef $format
*110      format (' ',a,' ',$)
*@else
*110      format (' ',a)
*@endif
*@endif
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') goto 10
      elseif (default.eq.'$') then
         write (6,110) prompt
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') string='1000000'
      else
*@ifdef f90
         write (6,"(1x,a,' <',a,'> ')",advance='no') prompt, default
*@else
*         write (6,100) prompt, default
*@ifdef $format
*100      format (' ',a,' <',a,'> ',$)
*@else
*100      format (' ',a,' <',a,'>')
*@endif
*@endif
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') string=default
      endif
      read (string,'(bn,f16.0)',err=9) anum
      return
      end
c deck getint
      subroutine getint (prompt,inum,default)
c     # reads an integer inum from the keyboard with prompt and default
c     # if default.eq.'$', 1000000 is returned on <cr>
c     # if default.eq.' ', there is no default.  reprompt on <cr>.
      character*(*) prompt, default
      character*16 string
9     if (default.eq.' ') then
10       continue
*@ifdef f90
         write (6,110,advance='no') prompt
110      format(1x,a)
*@else
*         write (6,110) prompt
*@ifdef $format
*110      format (' ',a,' ',$)
*@else
*110      format (' ',a)
*@endif
*@endif
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') goto 10
      elseif (default.eq.'$') then
         write (6,110) prompt
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') string='1000000'
      else
*@ifdef f90
         write (6,"(1x,a,' <',a,'> ')",advance='no') prompt, default
*@else
*         write (6,100) prompt, default
*@ifdef $format
*100      format (' ',a,' <',a,'> ',$)
*@else
*100      format (' ',a,' <',a,'>')
*@endif
*@endif
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') string=default
      endif
      read (string,'(bn,i16)',err=9) inum
      return
      end
c deck getint2
      subroutine getint2 (prompt, inum, default)
c     # reads an integer inum from the keyboard with prompt and default
c     # if default.eq.'$', 1000000 is returned on <cr>
c     # if default.eq.' ', there is no default.  reprompt on <cr>.
      character*(*) prompt, default
      character*16 string
9     if (default.eq.' ') then
10       continue
*@ifdef f90
         write (6,110,advance='no') prompt
110      format(1x,a,1x)
*@else
*         write (6,110) prompt
*@ifdef $format
*110      format (' ',a,' ',$)
*@else
*110      format (' ',a)
*@endif
*@endif
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') goto 10
      elseif (default.eq.'$') then
         write (6,110) prompt
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') string='1000000'
      else
*@ifdef f90
         write (6,"(1x,a,' <',a,'> ')",advance='no') prompt, default
*@else
*         write (6,100) prompt, default
*@ifdef $format
*100      format (' ',a,' <',a,'> ',$)
*@else
*100      format (' ',a,' <',a,'>')
*@endif
*@endif
         read (5,'(a)') string
         write(13,'(a)') string
         if (string.eq.' ') string=default
      endif
      read (string,'(bn,i16)',err=9) inum
      return
      end
c deck getyn
      subroutine getyn (prompt, iyn)
c     # ask a yes or no question with prompt.  use old iyn
c     # ('no'.eq.0,'yes'.ne.0) as default, return result in iyn.
      character*(*) prompt
      character*3 yn
c
      logical  compabb
      external compabb
c
10    if (iyn.lt.0) then
         call getchr (prompt,yn,' ')
      elseif (iyn.gt.0) then
         call getchr (prompt,yn,'YES')
      else
         call getchr (prompt,yn,'NO')
      endif
      if (compabb(yn,'YES',1)) then
         iyn=1
      elseif (compabb(yn,'NO',1)) then
         iyn=0
      else
         goto 10
      endif
      return
      end
c deck compabb
      function compabb(answer,refrnc,minlen)
c     # test if answer equals refrnc to as many characters as are
c     # in answer.  answer must be at lest minlen characters long
      logical compabb
      character*(*) answer, refrnc
c
      integer  fstrlen
      external fstrlen
c
      i=fstrlen(answer)
      if (i.gt.0) then
         compabb=(i.ge.minlen).and.(answer(1:i).eq.refrnc(1:i))
      else
         compabb=(i.ge.minlen)
      endif
      return
      end

      SUBROUTINE UPRCASE (ANSWER)
C  converts ANSWER to upper case letters
      CHARACTER*(*) ANSWER
      IA=ICHAR('a')
      IZ=ICHAR('z')
      IDIF=ICHAR('A')-IA
      DO 10 I=1,LEN(ANSWER)
        J=ICHAR(ANSWER(I:I))
        IF ((J.GE.IA).AND.(J.LE.IZ)) ANSWER(I:I)=CHAR(J+IDIF)
   10 CONTINUE
      RETURN
      END



      subroutine headwr(unit,program,version)

      implicit logical(a-z) 
      integer unit
      character*(*) program,version
      character*8  scr1,scr2,scr3
      integer i,ilenx,ilen1,ilen2
       external ilenx
      ilen1=ilenx(program)
      ilen2=ilenx(version)
       do i=1,min(8,ilen1)
        scr1(i:i)=program(i:i)
       enddo
       do i=ilen1+1,8
        scr1(i:i)=' '
       enddo
       do i=1,min(8,ilen2)
        scr2(i:i)=version(i:i)
       enddo
       do i=ilen2+1,8
        scr2(i:i)=' '
       enddo

        scr3='5.9.a  '

      write(unit,100) scr1,scr2,scr3
 100  format(//,
     .  5x,'******************************************'/
     .  5x,'**    PROGRAM:             ',1x,a8,4x,'**'/
     .  5x,'**    PROGRAM VERSION:     ',1x,a8,4x,'**'/
     .  5x,'**    DISTRIBUTION VERSION:',1x,a8,4x,'**'/
     . ,5x,'******************************************'/)

      return
      end

      FUNCTION ILENX(IC)
C  returns the length of IC excluding trailing blanks
      CHARACTER*(*) IC
      DO 100 I=LEN(IC),1,-1
        IF (.NOT.(IC(I:I).EQ.' '.OR.IC(I:I).EQ.CHAR(0))) THEN
          ILENX=I
          RETURN
        ENDIF
  100 CONTINUE
      ILENX=0
      RETURN
      END

