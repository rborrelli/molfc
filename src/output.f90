module output

    use parameters
    use system_type
    use sysop

    implicit none

	private

	public :: print_headers, print_structure, print_normal_modes, job_headers, job_endings, &
	          print_system, print_job, print_input

    contains

        subroutine print_headers (state, iuf)
        
        type(state_t), intent(in) :: state
        integer iuf, J

        J = int((44-(len_trim(state%id)+5))/2)
        write(iuf,165) trim(state%id)

165   	format(//,19x,'+------------------------------------------+',/,&
     	          19x,'|            ELECTRONIC STATE              |',/,&
                  19x,'|',<J>x    '| ',a,' >'<J>x,    '|',/,&
     		      19x,'+------------------------------------------+',/)

        end subroutine print_headers

!==============================================================================!

        subroutine print_structure (state, iuf)

        type(state_t), intent(in) :: state
        integer iuf
        integer i, j, nm, nato, natof, natoi

        write(fout,30)
        natoi = 1
        do nm = 1, state%nmol
          write(fout,50)  trim(state%molecule%id)
          write(fout,35)
          do i = 1, state%molecule%structure%numat
            write (iuf,40) &
            state%molecule%structure%atom(i)%elem%sym, &
            state%molecule%structure%atom(i)%elem%AM,  &
            (state%molecule%structure%atom(i)%coord(j), j=1,3)
          end do
          write(iuf,'(/)')
        end do

30      format (/,2x,'----------------------',/, &
                  2x,'Equilibrium geometries',/,&
                  2x,'----------------------',/)
35      format(2x,'atom',5x,'mass',9x,'x',14x,'y',14x,'z',/)
40      format(2x,a,5x,f9.4,3(f12.5,2x))
50      format(2x,'---------',/,&
               2x,'Molecule:',2x,a,/,&
               2x '---------',/)

        return
        end subroutine print_structure 

!==============================================================================!

      	subroutine print_normal_modes ( state, iuf ) 

      	implicit none

        type(state_t), intent(in) :: state
        integer iuf ! IUF is currently not used: for future use...
      	real(kind=dp), allocatable :: aa(:,:), bb(:)
      	character (len=2), allocatable :: ca(:)
        character xx(3)
        integer ka, kc, kb, i, j, k, ia, la, lb, lc, ix, i4, i3, np
        integer nato, nvib, imol

        write(fout,165) trim(state%id)

      	xx(1) = 'x'; xx(2) = 'y'; xx(3) = 'z'
      	ix = 0; i3 = 0; i4 = 0
	
        write(fout,166) trim(state%molecule%id)

        nato = state%molecule%structure%numat
        nvib = state%molecule%nvib

        allocate(aa(1:3*nato,1:nvib),bb(1:nvib),ca(1:nato))
        ca = get_symbol(state)
        aa = get_vibrations(state)
        bb = get_frequencies(state)

      	ka=1
      	kc=6
   10 	kb=min0(kc,nvib)
      	write(fout,50) (i,i=ka,kb)
      	write(fout,60) (bb(i),i=ka,kb)
      	write(fout,70)
      	la=1
      	lc=40
   20 	lb=min0(lc,3*nato)
      	np=0
      	do 30 i=la,lb
		i3 = int ((i-1)/3 )
		ix = ix + 1
		if (ix > 3) ix = 1
	    if ( i4 == 0 ) then
              write(fout,80) ca(i3+1),xx(ix),(aa(i,j),j=ka,kb)
	      i4 = i4 + 1
	    else
              write(fout,80) ' ',xx(ix),(aa(i,j),j=ka,kb)
	      i4 = i4 + 1
	    end if
		if (i4 > 2) i4 = 0
          np=np+1
          if (np.lt.10) go to 30
          np=0
   30 	continue
      	if (lb == 3*nato) go to 40
      	la=lc+1
      	lc=lc+40
      	go to 20
   40 	if (kb == nvib) then
          deallocate(aa,bb,ca)
		  return
        end if
      	ka=kc+1
      	kc=kc+6
      	go to 10

   50 	format (//x,9h Mode    ,i5,9i12)
   60 	format (/,2x,'Freq.',2x,f8.2,10(4x,f8.2))
   70 	format (2h  )
   80 	format (2x,a3,a2,10f12.5)
   90 	format (1h1)
165   	format(//,2x,'--------------------------------------------------------------',/,&
     	          2x,'Mass-weighted normal coordinates of electronic state  | ',a '>',/,&
     		  2x,'--------------------------------------------------------------',/)
166   	format(//,2x,'---------',/,&
           	  2x,'Molecule:',2x,a,/,&
     		  2x,'---------',/)

        return
      	end subroutine print_normal_modes

!==============================================================================!

       subroutine job_headers

       external hostnm, getlog, getpid
       integer getpid
       character (len=80) getlog, hostnm

       integer date_time (8)
       character (len = 12) real_clock (3), month(12), &
                          actual_month

       data  month(1), month(2), month(3), month(4), month(5),  & 
             month(6), month(7), month(8), month(9), month(10), &
             month(11), month(12)/                             &
             'January', 'February', 'March', 'April', 'May',    &
             'June', 'July', 'August', 'September', 'October',  &
             'November', 'December'                           /      

       call date_and_time (real_clock (1), real_clock (2), &
                    real_clock (3), date_time)

       actual_month = month(date_time(2))

       write(fout,10)
       write(fout,12)date_time(3),trim(actual_month),date_time(1), &
                 date_time(5),date_time(6),date_time(7) 
       !write(fout,199)getpid(), getlog(), hostnm()
       write(fout,198)getpid(), hostnm()
!       write(fout,11)

10   format( 19X,'+--------------------------------------------------+',/,&
             19X,'|          Molecular Franck-Condon Factors         |',/,&
             19X,'|               Calculation  Package               |',/,&  
             19x,'|             Version 3.0       Jan 2015           |',/,&
             19X,'|                                                  |',/,&  
             19X,'|              Theoretical Chemistry               |',/,&
             19X,'|               University of Torino               |',/,&
             19X,'|                                                  |',/,&  
             19x,'|                                                  |',/,&
             19x,'|            Implemented and mantained by          |',/,&
             19X,'|                Raffaele Borrelli                 |',/,&
             19x,'|            raffaele.borrelli@unito.it            |',/,&
             19X,'+--------------------------------------------------+',//) 
11     format(2X,'-----------------------',/, &
              2X,'INPUT DATA SET COUNTERS',/, &
              2X,'-----------------------',/) 
12     format(/,22x,'==============================================',&
              /,22x,'Job started on:',2x,I2,1X,A9,1X,I4, &
              ' at ',I2,':',I2,':',I2,/, &
                22x,'==============================================',//)
198    format (//,2x,'=========================================',/,&
                  2x,'Process ID: ', i7, /, & 
                  2x,'Running on host: ', a50,/,&
                  2x,'=========================================',/)
199    format (//,2x,'=========================================',/,&
                  2x,'Process ID: ', i7, /, &
                  2x,'User: ', a80, /, &
                  2x,'Running on host: ', a50,/,&
                  2x,'=========================================',/)

       return
       end subroutine job_headers

!==============================================================================!

       subroutine job_endings

       integer, parameter :: nmonths = 12
       integer date_time (8)
       character (len = 12) real_clock (3), month(12),  actual_month
       
       data  month(1), month(2), month(3), month(4), month(5),  & 
             month(6), month(7), month(8), month(9), month(10), &
             month(11), month(12)/                             &
             'Jenuary', 'February', 'March', 'April', 'May',    &
             'June', 'July', 'August', 'September', 'October',  &
             'November', 'December'                           /      

       call date_and_time (real_clock (1), real_clock (2), &
                    real_clock (3), date_time)

       actual_month = month(date_time(2))

       write(fout,12)date_time(3),trim(actual_month),date_time(1), &
                 date_time(5),date_time(6),date_time(7) 

       include 'formats'
       return
       end subroutine job_endings

!==============================================================================!

       subroutine print_system(system)

       type(system_t), intent(in) :: system
       
       integer is, i, j
       
       ! Print system data
       write(fout,872) system%nstate, system%state(1)%nmol, system%state(1)%nvib
 
       ! if input is a model system then exit
       if (system%model) return
        
!      if (control%printlevel >= 2) then
        do i = 1, system%nstate
          !------------------------------------
          ! Scrive le informazioni sugli stati.
          !------------------------------------
          call print_headers(system%state(i), fout)
          !------------------------------------------------------------
          ! Scrive le strutture delle molecole presenti nei vari stati.
          !------------------------------------------------------------
          call print_structure(system%state(i), fout)
        end do
!      end if
    
!      if (control%printlevel >= 3) then
        do i = 1, system%nstate
          !--------------------------------------------------------------
          ! Scrive i modi normali delle molecole presenti nei vari stati.
          !--------------------------------------------------------------
          call print_normal_modes(system%state(i), fout)
        end do
!      end if


40     format(2x,'-----',/&
              2x,'State |',x,a,x,'>'/,&
              2x,'-----',/)
872    format (/,2x,'+-------------------+',/,&
                 2x,'| System Parameters |',/,& 
                 2x,'+-------------------+',//,&
                 2x,'Number of Electronic States :',i4,//,&
                 2x,'Number of Molecules in each state :',i4,//,&
                 2x,'Total Number of Vibrational Modes :',i4,/)

       return
       end subroutine print_system

!==============================================================================!

       subroutine print_job(fcjob)

       use fc_type
       use xmvar, only : system

       implicit none

       type(fc_t) :: fcjob
       integer i, j, nfctot, is, im, ig

       write(fout,10) !trim(adjustl(fcjob%bra)), trim(adjustl(fcjob%ket))

       nfctot = 1
       do ig = 1, fcjob%ngroup
           write(fout,15) ig
           ! --- Totale stati vibronici
           if (allocated(fcjob%group(ig)%active)) then
             do i = 1, size(fcjob%group(ig)%active)
                nfctot = nfctot * product(fcjob%group(ig)%active(i)%mode(:)%nq + 1)
                write(fout,20) fcjob%group(ig)%active(i)%state
                write(fout,30) fcjob%group(ig)%active(i)%molecule
                write(fout,8723)
                do j = 1, fcjob%group(ig)%active(i)%nact
                    is = get_state_from_id(system,fcjob%group(ig)%active(i)%state)
                    write(fout,8722) fcjob%group(ig)%active(i)%mode(j)%id,   &
                                     system%state(is)%molecule%normodes%vibration(fcjob%group(ig)%active(i)%mode(j)%id)%freq, &
                                     fcjob%group(ig)%active(i)%mode(j)%nq
                end do
                write(fout,'(/)')
             end do

           end if
       end do

       write(fout,875) nfctot

10     format (//2x,'+--------------------------------------+',/, &
                 2x,'| Franck-Condon Integrals Calculation  |',/, &
                 2x,'+--------------------------------------+',/)
15     format (/,2x,'+---------------+',/,&
                 2x,'| Group ', i2,6x,'|',/,&
                 2x,'+---------------+')
20     format (/2x,'==> Excited state:',x,a,/)
30     format (4x,'--> Excited Molecule:',2x,a,/)
870    format (/,2x,'+---------------------+',/,&
                 2x,'| Current Job Options |',/,&
                 2x,'+---------------------+',//,&
                 2x,'Total Number of Vibronic States :',i10,/)
8741   format (2x,'Number of Vibrational States in state |',x,a,x'>:',x,i10,/)
875    format (2x,'Num. Franck-Condon Integrals :',x,i10,/)
876    format (  2x,'Excited Vibrations in Electronic State |',a,'> :',i3,/)
8721   format (/,2x,'+--------------------+',/,&
                 2x,'| Excited Vibrations |',/,& 
                 2x,'+--------------------+-',/) 
8723   format (4x,'Vibration ID',4x,'Frequency',4x,'Quanta (States=Quanta+1)',/) 
8722   format (4x,i3,13x,f7.2,4x,i3)

       end subroutine print_job

!==============================================================================!

       subroutine print_input (fname)

       implicit none

       integer i, k, iu
       character(len=40), intent(in) :: fname
       character (len = 100) string

       open (finp,file=fname,status='unknown')

       do 
         read (*,10,iostat=k) string
         if (k /= 0) exit
       !  write(fout,20)string(1:len_trim(string))
         write (fout,10)string(1:len_trim(string))
       end do

! Reset the read position of the file.
       rewind iu

       close (iu)

10     format (A)
20     format (2x,A)

       return        
       end subroutine print_input

end module output
