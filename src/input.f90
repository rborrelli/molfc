module input

    use parameters
    use proc_type
    use job_type
    use fc_type
    use dos_type
    use transf_type
    use errors
    use iomodes
    use iogeom
    use iofiles
    use util
    use sysop, only : get_state_from_id
    use xmvar, only : system, proc, job
    use FoX_dom

  
    implicit none

    private

    public :: read_input

    logical, parameter :: input_debug = .true.
  
contains
  
    !===============================================================================

    subroutine read_input(fname)

        character(len=*), intent(in) :: fname
        integer :: status

        call read_system(fname)  ! this uses the DOM API

        call read_proc(fname)  ! this uses the DOM API

        call read_jobs(fname)  ! this uses the DOM API

        return
    end subroutine read_input

    !===============================================================================

    subroutine read_system(fname)

        character(len=*), intent(in) :: fname

        character(len=4) :: units, coord, lbl
        character(len=80) :: id, ftype, fpath, line, filename
        character(len=STRLEN) :: stateid, molid, value
        character(len=200) :: lline, chunk
    
        integer :: status, ierr, alloc_err
        integer :: numat, nvib, nstate, nmol, i, j, ios, n0
        integer :: iu, is, im, ratom, istate, imol, ifile
        integer :: isol, ieol, iam, anum, ist, lenstr, iz, xnodeid
        logical :: linear, readin, hpmodes, bool
        real(kind=dp) :: rval, zcpl, muc
        real(kind=dp) :: amass, x, y, z, fscale
        real(kind=dp), allocatable :: freq(:)

        type(Node), pointer :: xmlDoc
        type(Node), pointer :: xnode, np, nsNode, nstru, stateNode, cNode, muNode, dofNode
        type(NodeList), pointer :: xnodeList, nList, nsList, nmList, nmuList, ntmList, cList, dofList

        !    type(dictionary_t) :: attributes
        character(len=STRLEN) :: str, substr

        filename = adjustl(fname)
        xmlDoc => parsefile(filename(1:len_trim(filename)))
        call normalize(xmlDoc)

        xnodeList => getElementsByTagName(xmlDoc, "input")
        if (getLength(xnodeList) == 0) ierr = error(0,"Missing <input> tag.")
        if (getLength(xnodeList) > 1) ierr = error(0,"Too many <input> tags (max = 1).")
    
        ! xnode => <input>
        xnode => item(xnodeList,0)
        ! get childs of <input>, ie all the files xnodes
        nList => getElementsByTagName(xnode, "system")
        if (getLength(xnodeList) == 0) ierr = error(0,"Missing <system>.")
        if (getLength(xnodeList) > 1) ierr = error(0,"Too many <system> tags (max = 1).")

        xnode => item(nList,0)
        nsList => getElementsByTagName(xnode,"state")
        ! get List of <tm> xnodes
        ntmList => getElementsByTagName(xnode,"tm")
   
        nstate = getLength(nsList)
        system%nstate = nstate
        allocate(system%state(1:nstate)) !Alloca gli stati elettronici

        value = getAttribute(xnode,"type")
        if (.not.isempty(value)) then
            if (value == "model") system%model = .true.
        !      if (ios /= 0) ierr = error(0,"Cannot read system type")
        end if

        !-----------------------------------------------------------------------!
        ! If we use a model system and not real molecule then read the model,   !
        ! set up fictious states and molecules and exit                         !
        !-----------------------------------------------------------------------!
        if (system%model) then
            do is = 1, system%nstate
                xnode => item(nsList,is-1)

                value = getAttribute(xnode,"id")
                if (.not.isempty(value)) system%state(is)%id = value
	
                value = getAttribute(xnode,"energy")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%energy
                    if (ios /= 0) ierr = error(0,"Cannot read state energy")
                end if

                ! alloca una molecola "finta"
                system%state(is)%nmol = 1
                system%state(is)%molecule%id = "molecule" ! set molecule id to "molecule"
            
                ! Legge <freq>
                nmList => getElementsByTagName(xnode,"freq")
                xnode => item(nmList,0)

                im = 1 ! set molecule number
            
                value = getAttribute(xnode,"nvib")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) nvib
                    system%state(is)%molecule%nvib = nvib
                    system%state(is)%nvib = nvib
                    if (ios /= 0) ierr = error(0,"Cannot read number of vibrations in model.")
                end if

                ! allocate NVIB vibrations
                allocate(system%state(is)%molecule%normodes%vibration(1:nvib))

                ! Check if we have a file attribute
                value = ""
                if (hasAttribute(xnode,"file")) then
                    value = getAttribute(xnode,"file")
                end if
                if (.not.isempty(value)) then
                    !---------------------------------+
                    ! Read frequencies from file      |
                    !---------------------------------+
                    call openfile(fkm, trim(value))
                    !open(unit=fkm,file=trim(value),status="old")
                    read(fkm,*)(system%state(is)%molecule%normodes%vibration(i)%freq,i=1,nvib)
                    close(fkm)
                else
                    !---------------------------------+
                    ! Read frequencies from xml file  |
                    !---------------------------------+
                    allocate(freq(1:nvib))
                    call extractDataContent(xnode,freq)
                    system%state(is)%molecule%normodes%vibration(1:nvib)%freq = freq
                    deallocate(freq)
                end if
                ! allocate NVIB vibrations
                !allocate(system%state(is)%molecule%normodes%vibration(1:nvib))
                !xnode => getFirstChild(xnode)
                !str = getNodeValue(xnode)
                !call build_data_array(str,system%state(is)%molecule%normodes%vibration(1:nvib)%freq)
                ! controlla le freqeuenze del modello non siano troppo basse o negative (per errore di input...)
                if (any(system%state(is)%molecule%normodes%vibration(1:nvib)%freq <= 1.0e-2)) &
                ierr = error(0,"Vibrational frequencies too low!!")
            end do
            ! Read transition moment
            if (associated(ntmList)) then
                call read_tm(ntmList)
            end if
            return
        end if
    
        ifile = 0
        ! Nei cicli do non si devono modificare nsList e nmList !!
        ST: do is = 1, system%nstate

            xnode => item(nsList, is - 1)
            value = getAttribute(xnode,"id")
            if (.not.isempty(value)) system%state(is)%id = value

            value = getAttribute(xnode,"energy")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) system%state(is)%energy
                if (ios /= 0) ierr = error(0,"Cannot read state energy")
            end if

            value = getAttribute(xnode,"zcpl")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) system%state(is)%zcpl
                if (ios /= 0) ierr = error(0,"Cannot read state zcpl")
            end if

            ! Get number of molecule in state IS.
            nmList => getElementsByTagName(xnode,"molecule")
            nmol = getLength(nmList)
            system%state(is)%nmol = nmol
            if (nmol > 1) ierr = error(0,"Only one molecule allowed per calculation.")

            MOL:    do im = 1, nmol

                xnode => item(nmList, im - 1)
                !+----------------+
                ! Get molecule id |
                !+----------------+
                value = getAttribute(xnode,"id")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%molecule%id
                    if (ios /= 0) ierr = error(0,"Cannot read molecule id")
                end if

                !+------------------------+
                !  Parse xnode <structure> |
                !+------------------------+
                nList => GetElementsByTagName(xnode, "structure")
                if (getLength(nList) == 0) ierr = error(0,"Missing <structure>.")
                if (getLength(nList) > 1) ierr = error(0,"Too many <structure> tags (max = 1).")
                np => item(nList, 0)

                value = getAttribute(np,"nat")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) numat
                    system%state(is)%molecule%structure%numat = numat
                    if (ios /= 0) ierr = error(0,"Cannot read atom number")
                else
                    ierr = error(0,"Missing number of atoms (nat) in <structure>")
                end if

                !-----------------------------+
                ! Allocate molecule structure |
                !-----------------------------+
                allocate(system%state(is)%molecule%structure%atom(1:numat))
                value = getAttribute(np,"linear")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%molecule%structure%linear
                    if (ios /= 0) ierr = error(0,"Cannot read linear")
                end if

                value = getattribute(np,"units")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%molecule%structure%units
                    if (ios /= 0) ierr = error(0,"Cannot read geometry units")
                end if

                value = getattribute(np,"coord")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%molecule%structure%coord
                    if (ios /= 0) ierr = error(0,"Cannot read geometry units")
                end if

                nList => GetElementsByTagName(np,"file")
                np => item(nList, 0)

                if (associated(np)) then
                    system%state(is)%molecule%structure%fromfile = .true.
                    ifile = ifile + 1
                    system%state(is)%molecule%structure%file%iu = inpf + ifile
                    value = getAttribute(np,"path")
                    if (.not.isempty(value)) then
                        system%state(is)%molecule%structure%file%path = value
                      !read(value,*,iostat=ios) system%state(is)%molecule%structure%file%path
                      !if (ios /= 0) ierr = error(0,"Cannot read file path")
                    end if

                    value = getAttribute(np,"type")
                    if (.not.isempty(value)) then
                        read(value,*,iostat=ios) system%state(is)%molecule%structure%file%type
                        if (ios /= 0) ierr = error(0,"Cannot read file type")
                    end if

                    np => getFirstChild(np)
                    if (.not.associated(np)) exit
                    if (getNodeType(np) /= TEXT_NODE) exit
                    str = getNodeValue(np)
                    system%state(is)%molecule%structure%file%name = str
                end if

                !+--------------------------+
                ! Parse xnode <normal_modes> |
                !+--------------------------+
                nList => GetElementsByTagName(xnode, "normal_modes")
                if (getLength(nList) == 0) ierr = error(0,"Missing <normal_modes>.")
                if (getLength(nList) > 1) ierr = error(0,"Too many <normal_modes> tags (max = 1).")
                np => item(nList, 0)
                value = getAttribute(np,"massw")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%molecule%normodes%massw
                    if (ios /= 0) ierr = error(0,"Cannot read MASSW attribute: must be either T or F.")
                end if

                value = getAttribute(np,"units")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%molecule%normodes%units
                    if (ios /= 0) ierr = error(0,"Cannot read units attribute.")
                end if

                value = getAttribute(np,"sort")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) system%state(is)%molecule%normodes%sort
                    if (ios /= 0) ierr = error(0,"Cannot read SORT attribute: must be either T or F.")
                end if


                value = getAttribute(np,"scale")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) fscale
                    system%state(is)%molecule%normodes%fscale = fscale
                    if (ios /= 0) ierr = error(0,"Cannot read scale factors.")
                end if

                nList => GetElementsByTagName(np,"file")
                np => item(nList, 0)

                if (associated(np)) then
                    system%state(is)%molecule%normodes%fromfile = .true.
                    ifile = ifile + 1
                    system%state(is)%molecule%normodes%file%iu = inpf + ifile
                    value = getAttribute(np,"path")
                    if (.not.isempty(value)) then
                        system%state(is)%molecule%normodes%file%path = value
                      !print *, value(1:len_trim(value))
                      !print *, system%state(is)%molecule%normodes%file%path
                      !if (ios /= 0) ierr = error(0,"Cannot read file path")
                    end if

                    value = getAttribute(np,"type")
                    if (.not.isempty(value)) then
                        read(value,*,iostat=ios) system%state(is)%molecule%normodes%file%type
                        if (ios /= 0) ierr = error(0,"Cannot read file type")
                    end if
               
                    np => getFirstChild(np)
                    if (.not.associated(np)) exit
                    if (getNodeType(np) /= TEXT_NODE) exit
                    str = getNodeValue(np)
                    system%state(is)%molecule%normodes%file%name = str
                end if

                !---------------------------------------------------+
                ! Alloca 3N-6 (3N-5 per specie lineari) vibrazioni  |
                !---------------------------------------------------+
                nvib = 3*numat - 6
                if (system%state(is)%molecule%structure%linear .or.(numat == 2)) nvib = 3*numat - 5
                system%state(is)%molecule%nvib = nvib
                allocate(system%state(is)%molecule%normodes%vibration(1:nvib))
                do i = 1, nvib
                    allocate (system%state(is)%molecule%normodes%vibration(i)%atom(1:numat),stat=alloc_err)
                    if (alloc_err /= 0) ierr = error(0,"Cannot allocate atom in system%state.")
                end do

            end do MOL
    
        end do ST

        ! Read transition moments
        call read_tm(ntmList)

        !------------------------------------------------------------------------------+
        ! A questo punto legge geometria e modi normali dai file specificati in input. |
        ! oppure direttamente dal file xml.                                            |
        !------------------------------------------------------------------------------+
        do is = 1, system%nstate
            stateNode => item(nsList,is-1)
            nmList => getElementsByTagName(stateNode,"molecule")
            do im = 1, system%state(is)%nmol

                if (system%state(is)%molecule%structure%fromfile) then
                    iu = system%state(is)%molecule%structure%file%iu
                    fpath = system%state(is)%molecule%structure%file%path
                    filename = system%state(is)%molecule%structure%file%name
                    write(fout,89) adjustl(system%state(is)%id(1:len_trim(system%state(is)%id))), &
                    adjustl(system%state(is)%molecule%id(1:len_trim(system%state(is)%molecule%id))), &
                    adjustl(fpath(1:len_trim(fpath)))//&
                    adjustl(filename(1:len_trim(filename)))
                    call openfile(iu, filename, filedir=fpath)
                    ftype = trim(system%state(is)%molecule%structure%file%type)
                    select case(ftype)
                        case("gaussian","g94","g98","g03","g09")
                            call read_geometry_gau(iu,is,im,ratom)
                        case("turbo","tm","tmol")
                            call read_geometry_tm(iu,is,im,ratom)
                        case("xyz","cart")
                            call read_geometry_xyz(iu,is,im,ratom)
                        case("molden")
                            call read_geometry_molden(iu,is,im,ratom)
                        case("mopac")
                            call read_geometry_mopac(iu,is,im,ratom)
                        case("d2k","demon")
                            call read_geometry_demon(iu,is,im,ratom)
                        case default
                            call read_geometry_gau(iu,is,im,ratom)
                    end select
                else
                    !+----------------------------------------+
                    ! Read structure directly from xml file.
                    !+----------------------------------------+
                    np => item(nmList,im-1)
                    nList => GetElementsByTagName(np, "structure")
                    np => item(nList, 0)
                    nstru => getFirstChild(np)
                    if (.not.associated(nstru)) ierr = error(0,"Cannot read structure.")
                    str = getNodeValue(nstru)
                    lenstr = len(str)
                    ieol = 0
                    do iam = 1, numat
                        isol = ieol + 1
                        ieol = ieol + index(extract(str,isol,lenstr),char(10)) ! find first EOL in substr
                        chunk = extract(str,isol,ieol-1) ! extract single line
                        !+---------------------------+
                        ! I Format:
                        ! lbl, amass, x, y, z
                        !+---------------------------+
                        read(chunk(1:len_trim(chunk)),*,iostat=ist)lbl, amass, x, y, z
                        if (ist == 0) then
                            system%state(is)%molecule%structure%atom(iam)%elem%Sym = lbl
                            system%state(is)%molecule%structure%atom(iam)%elem%AM = amass
                            system%state(is)%molecule%structure%atom(iam)%elem%Z = elmnts(mapel(lbl))%Z
                            system%state(is)%molecule%structure%atom(iam)%coord(1:3) = (/x,y,z/)
                        else
                            !+---------------------------+
                            ! II Format:
                            ! anum, x, y, z
                            !+---------------------------+
                            read(chunk(1:len_trim(chunk)),*,iostat=ist)anum, x, y, z
                            ! --- Il compilatore IFC 9.0 ha problemi senza questo controllo.
                            if (anum == 0) ist = 1
                            if (ist == 0) then
                                system%state(is)%molecule%structure%atom(iam)%elem = elmnts(mapel(anum))
                                system%state(is)%molecule%structure%atom(iam)%coord(1:3) = (/x,y,z/)
                            else
                                !+---------------------------+
                                ! III Format:
                                ! lbl, x, y, z
                                !+---------------------------+
                                read(chunk(1:len_trim(chunk)),*,iostat=ist)lbl, x, y, z
                                if (ist == 0) then
                                    system%state(is)%molecule%structure%atom(iam)%elem = elmnts(mapel(lbl))
                                    system%state(is)%molecule%structure%atom(iam)%coord(1:3) = (/x,y,z/)
                                else
                                    ierr = error(0,"Cannot read atoms in structure.")
                                end if
                            end if
                        end if
                    end do
                end if

                close(iu)

                if (system%state(is)%molecule%normodes%fromfile) then
                    iu = system%state(is)%molecule%normodes%file%iu
                    fpath = system%state(is)%molecule%normodes%file%path
                    filename = system%state(is)%molecule%normodes%file%name
                    write(fout,88) adjustl(system%state(is)%id(1:len_trim(system%state(is)%id))), &
                    adjustl(system%state(is)%molecule%id(1:len_trim(system%state(is)%molecule%id))), &
                    adjustl(fpath(1:len_trim(fpath)))//&
                    adjustl(filename(1:len_trim(filename)))
                    call openfile(iu, filename, filedir=fpath)
                    ftype = trim(system%state(is)%molecule%normodes%file%type)
                    select case(ftype(1:len_trim(ftype)))
                        case("ghp","g94hp","g98hp","g03hp","g09hp")
                            system%state(is)%molecule%normodes%massw = .false. ! ?? superfluo perch?? ?? di default
                            call read_normal_modes_gau(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                        case("gaussian","g94","g98","g03","g09")
                            system%state(is)%molecule%normodes%massw = .false. ! ?? superfluo perch?? ?? di default
                            call read_normal_modes_gaunohp(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                        case("turbo","tm","tmol")
                            system%state(is)%molecule%normodes%massw = .false. ! ?? superfluo perch?? ?? di default
                            call read_normal_modes_tm(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                        case("molden")
                            system%state(is)%molecule%normodes%massw = .true. ! vale per quelle scritte da ACES 2
                            call read_normal_modes_molden(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                        case("mopac","mopac09")
                            system%state(is)%molecule%normodes%massw = .true. ! vale per quelle scritte da MOPAC09
                            call read_normal_modes_mopac(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                        case("gamess","firefly")
                            system%state(is)%molecule%normodes%massw = .false. ! vale per quelle scritte da Firefly8.0
                            call read_normal_modes_gamess(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                        case("d2k","demon")
                            system%state(is)%molecule%normodes%massw = .false. ! check this for deMon2K
                            call read_normal_modes_demon_nohp(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                        case default
                            system%state(is)%molecule%normodes%massw = .false. ! ?? superfluo perch?? ?? di default
                            call read_normal_modes_gau(iu,is,im,readin)
                            if (.not.readin) ierr = error(0,"Cannot read normal modes from input file "//filename)
                    end select
                end if
                close(iu)
            end do
        end do

        do is = 1, system%nstate
            do im = 1, system%state(is)%nmol
                system%state(is)%nvib = system%state(is)%nvib + system%state(is)%molecule%nvib
            end do
        end do


	    include 'formats'
	
        return
    end subroutine read_system
    
    !==============================================================================

    subroutine read_proc(fname)

        character(len=*), intent(in) :: fname

        type(Node), pointer :: xmlDoc
        type(Node), pointer :: xnode, np
        type(NodeList), pointer :: xnodeList, nList

        character(len=STRLEN) :: filename
        character(len=STRLEN) :: value
        integer :: i, i1, alloc_err, ierr

        filename = adjustl(fname)
        xmlDoc => parsefile(trim(adjustl(filename)))
        call normalize(xmlDoc)

        xnodeList => getElementsByTagName(xmlDoc, "proc")
        if (getLength(xnodeList) > 1) ierr = error(0,"Too many <proc> tags.")
        if (getLength(xnodeList) == 0) return

        !print *, 'reading proc'
        ! xnode => <proc>
        xnode => item(xnodeList,0)

        nList => getElementsByTagName(xnode,"reorder")

        if (getLength(nList) > 0) then
            allocate(proc%reorder(1:getLength(nList)),stat=alloc_err)
            if (alloc_err /= 0) ierr = error(0,"Cannot allocate atom sequence in <reorder>")

            ! read <reorder>
            do i = 0, getLength(nList) - 1
            !print *, 'reading reorder '
                np => item(nList,i)
                i1 = i + 1
                proc%reorder(i1)%data =  getAttribute(np,"data")
                if (isempty(proc%reorder(i1)%data)) proc%reorder(i)%data = "atoms"
     
                value = getAttribute(np,"reference")
                if (.not.isempty(value)) proc%reorder(i1)%reference =  getAttribute(np,"reference")

                proc%reorder(i1)%state =  getAttribute(np,"state")
                if (isempty(proc%reorder(i1)%state)) ierr = error(0,"Cannot read state in <reorder>")
                ! if state is "all" reorder all states (see sysop).
                !if (proc%reorder(i1)%state == "all") proc%reorder(i1)%all = .true.

                proc%reorder(i1)%molecule =  getAttribute(np,"molecule")
                if (isempty(proc%reorder(i1)%molecule)) ierr = error(0,"Cannot read molecule <reorder>")
     
                np => getFirstChild(np)
                if ((.not.associated(np))) then
     !                           print *, 'reo cycle'

                    ! if ((.not.associated(np)) .and. isempty(proc%reorder(i1)%state)) then
                    cycle
                else if (associated(np) .and. isempty(proc%reorder(i1)%state)) then
                    ierr = error(0,"Found order without a state specification. Check <reorder>.")
                    !else if ((.not.associated(np)) .and. (.not.isempty(proc%reorder(i1)%state))) then
                    !    ierr = error(0,"Found state without an order sequence. You must probably  &
                    !	                remove the attribute state. Check <reorder>.")
                else
                !print *, 'reo 1'
                    ! Read order sequence
                    proc%reorder(i1)%ord = getNodeValue(np)
                end if
            end do
        end if

        ! xnode => <frame>
        nList => getElementsByTagName(xnode,"frame")
        if (getLength(nList) > 0) then
          allocate(proc%frame(1:getLength(nList)),stat=alloc_err)
          if (alloc_err /= 0) ierr = error(0,"Cannot allocate frame sequence in <frame>")
        end if
        ! read <frame>
        do i = 0, getLength(nList) - 1
            np => item(nList,i)
            i1 = i + 1
            proc%frame(i1)%state =  getAttribute(np,"state")
            if (isempty(proc%frame(i1)%state)) ierr = error(0,"Cannot read state in <frame>")

            proc%frame(i1)%molecule =  getAttribute(np,"molecule")
            if (isempty(proc%frame(i1)%molecule)) ierr = error(0,"Cannot read molecule in <frame>")

            np => getFirstChild(np)
            if (.not.associated(np)) exit

            call extractDataContent(np,proc%frame(i1)%axis)
            !proc%frame(i1)%axis = getNodeValue(np)

        end do

        ! Read <subset>: this is used if we want to completely remove some vibrations from the calculation
        nList => getElementsByTagName(xnode,"subset")
        if (getLength(nList) > 0) then
        if (getLength(nList) == 2) then
            allocate(proc%subset(1:getLength(nList)),stat=alloc_err)
            if (alloc_err /= 0) ierr = error(0,"Cannot allocate subset sequence in <subset>")
            !print *, 'subset ', getLength(nList)
            do i = 0, getLength(nList) - 1
                np => item(nList,i)
                i1 = i + 1
                call get_subset(np,proc%subset(i1))
            end do
        else
            ierr = error(0,"Two <subset> directive must be specified. One for each electronic state.")
        end if
        end if

        return
    end subroutine read_proc

    !==============================================================================
    !
    subroutine read_jobs(fname)

        character(len=*), intent(in) :: fname

        type(Node), pointer :: xmlDoc
        type(Node), pointer :: xnode, np
        type(NodeList), pointer :: xnodeList, nList

        !    type(dictionary_t) :: attributes
        character(len=STRLEN) :: str

        character(len=80) :: filename
        character(len=STRLEN) :: stateid, molid
        integer :: njob, nact, is, im, ij, status, ierr
        integer :: ivib, nq, ival, i, alloc_err
        logical :: bool

        real(kind=dp) rval

        filename = adjustl(fname)
        xmlDoc => parsefile(filename(1:len_trim(filename)))
        call normalize(xmlDoc)

        xnodeList => getElementsByTagName(xmlDoc, "input")
        if (getLength(xnodeList) == 0) ierr = error(0,"Missing <input> tag.")
        if (getLength(xnodeList) > 1) ierr = error(0,"Too many <input> tags (max = 1).")
    
        ! xnode => <input>
        xnode => item(xnodeList,0)
        ! get childs of <input>, ie all the files xnodes
        nList => getElementsByTagName(xnode, "job")
        njob = getLength(nList)
        if (njob == 0) ierr = error(0,"One <job> must be specified")
        !if (njob == 0) ierr = error(0,"At least one <job> must be specified")
        !allocate(job(1:njob),stat=alloc_err) ! Allocate job array
        !if (alloc_err /= 0) ierr = error(0,"Cannot allocate JOBS.")

        !do i = 0, getLength(nList) - 1
        !    np => item(nList,i)
            !call get_job(np, job(i+1))
            np => item(nList,0)
            call get_job(np, job)
        !end do

        return
    end subroutine read_jobs

    !==============================================================================

    subroutine get_job(root_xnode, job)
    
        type(Node),  pointer, intent(in) ::  root_xnode
        type(job_t), intent(inout) :: job

        type(Node), pointer :: np, myNode
        type(NodeList), pointer ::  nList, nsList, nmList
        character(len=STRLEN) :: str
        character(len=STRLEN) :: value

        integer i, ierr, ios, ie, ifc, itrn
        integer nfcj, nevj, ntrnj

        !job%nmeth = 0
        job%nmeth = getLength(getElementsByTagName(root_xnode,"fc")) + &
        getLength(getElementsByTagName(root_xnode,"transform")) + &
        getLength(getElementsByTagName(root_xnode,"dos"))
        !getLength(getElementsByTagName(xnode,"fcwd"))
        !nList => getElementsByTagName(xnode,"fc")
        !job%nmeth = getLength(nList)
        nfcj = getLength(getElementsByTagName(root_xnode,"fc"))
        !nList => getElementsByTagName(xnode,"transform")
        !job%nmeth = job%nmeth +  getLength(nList)
        ntrnj = getLength(getElementsByTagName(root_xnode,"transform"))
        if (ntrnj > 1) ierr = error(0,"Only one <transform> allowed per job.")
        !job%nmeth = job%nmeth +  getLength(nList)
        !nevj =  getLength(nList)
    
        ! TODO: migliorare l'allocazione dei metodi nei singoli job.
        allocate(job%method(1:job%nmeth))
        if (nfcj > 0) allocate(job%fc(1:nfcj))
        allocate(job%trns)

        nList => getChildNodes(root_xnode)

        ie = 0
        ifc = 0
        itrn = 0
        do i = 0, getLength(nList) - 1
            np => item(nList,i)
            str = getNodeName(np)
            select case(str)
                case("fc")
                    ie = ie + 1
                    ifc = ifc + 1
                    job%method(ie) = str
                    call get_job_fc(np,job%fc(ifc))
                case("transform")
                    ie = ie + 1
                    itrn = itrn + 1
                    job%method(ie) = str
                    call get_job_transf(np,job%trns)
                case("dos")
               	    print *, str//' preliminary implementation'
            	    ie = ie + 1
            	    job%method(ie) = str
            	    call get_job_dos(np,job%dos)
            end select
        end do

        ! If non transformation are specified set an error. In future will switch to
        ! default settings.
        if (itrn == 0) ierr = error(0,"At least one <transform> tag must be used.")

        return
    end subroutine get_job

        !==============================================================================

    ! This subroutine allows to select only a few vibrations in the computation of the Duschinsky transformation

    SUBROUTINE get_subset(root_xnode, subset)

        type(Node), pointer, intent(in) :: root_xnode
        type(Node), pointer :: np, myNode
        type(NodeList), pointer :: npList
        type(subset_t), intent(out) :: subset

        integer :: i, i1, ierr
        character(len=3*STRLEN) :: value

            ! the molecule defining group IG
            subset%state = getAttribute(root_xnode,"state")
            value = getAttribute(root_xnode,"nvib")
            if (isempty(value)) then
                ! include all vibrations
                subset%nvib = system%state(1)%molecule%nvib
            else
                read(value,*) subset%nvib
                print *, 'SUBSET ',subset%nvib
            end if
            allocate(subset%incvib(1:subset%nvib))
            ! get included vibrations
            np => getFirstChild(root_xnode)
            if (associated(np)) then
                value = getNodeValue(np)
                ! get included vibrations
                subset%incvib(1:subset%nvib) = get_incvib(subset%nvib,value)
            else
                forall (i=1:subset%nvib) subset%incvib(i) = i
            end if

        return
    END SUBROUTINE get_subset

    !==============================================================================

	! This subroutine reads the  density of states (dos) job.

	SUBROUTINE get_job_dos(root_xnode, dosjob)

	type(Node), pointer :: root_xnode
	type(dos_t), intent(out) :: dosjob

    character(len=STRLEN) :: value
	integer :: ios, ierr

	dosjob%state = getAttribute(root_xnode,"state")
    if (isempty(dosjob%state)) dosjob%state = system%state(1)%id ! default value: state 1

    ! Controllare che la soluzione con "all" funzioni...
	dosjob%molecule = getAttribute(root_xnode,"molecule")
    if (isempty(dosjob%molecule)) dosjob%molecule = "all" ! default value all moecules included

	dosjob%file = getAttribute(root_xnode,"file")
	if (isempty(dosjob%file)) dosjob%file = "dos_states.dos"

	dosjob%method = getAttribute(root_xnode,"method")
    if (isempty(dosjob%method)) dosjob%method = "brsw" ! default: Beyer-Swinehart algorthm

	value = getAttribute(root_xnode,"emin")
    if (.not.isempty(value)) then
      read(value,*,iostat=ios) dosjob%emin
      if (ios /= 0) ierr = error(0,"Cannot read EMIN")
    end if

	value = getAttribute(root_xnode,"emax")
    if (.not.isempty(value)) then
      read(value,*,iostat=ios) dosjob%emax
      if (ios /= 0) ierr = error(0,"Cannot read EMAX")
    end if

	value = getAttribute(root_xnode,"egrain")
    if (.not.isempty(value)) then
      read(value,*,iostat=ios) dosjob%egrain
      if (ios /= 0) ierr = error(0,"Cannot read EGRAIN")
    end if

	return
	END SUBROUTINE get_job_dos

    !==============================================================================

    subroutine get_job_fc (root_xnode, fc_job)

        type(Node), pointer, intent(in) ::  root_xnode
        type(Node), pointer :: np, myNode, gxnode
        type(fc_t), intent(out) :: fc_job

        type(NodeList), pointer ::  nList, nsList, nmList, vibList, npList

        character(len=3*STRLEN) :: value, stateid, molid
        integer i, j, k, l, ivib, nq, is, im, ig, ios, ierr, alloc_err, ntvib
        integer ioff, imol
        integer, pointer :: exvib(:)
    
        ! Get some attributes of <fc>
        value =  getAttribute(root_xnode,"plevel")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%printlevel
            if (ios /= 0) ierr = error(0,"Cannot read plevel in <fc>")
        end if

        value =  getAttribute(root_xnode,"fcht")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%fcht
            if (ios /= 0) ierr = error(0,"Cannot read fcht in <fc>")
        end if

        value =  getAttribute(root_xnode,"class")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%class
            if (ios /= 0) ierr = error(0,"Cannot read class in <fc>")
        end if

        value =  getAttribute(root_xnode,"Temp")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%Temp
            if (ios /= 0) ierr = error(0,"Cannot read class in <fc>")
        end if

              !value =  getAttribute(root_xnode,"berk")
              !if (.not.isempty(value)) then
              !    read(value,*,iostat=ios) fc_job%berk
              !    if (ios /= 0) ierr = error(0,"Cannot read berk in <fc>")
              !end if

        value =  getAttribute(root_xnode,"pert")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%pert
            if (ios /= 0) ierr = error(0,"Cannot read pert")
        end if

        value = getAttribute(root_xnode,"order")
        if (.not.isempty(value)) then
            if (.not.fc_job%pert) ierr = error(0,'Please remove order attirbute or select a perturbative calculation.')
            read(value,*,iostat=ios) fc_job%order
            if (ios /= 0) ierr = error(0,"Cannot read perturbation order to be used.")
        end if

        ! E' inutile leggere questa variabile: viene automaticamente messa .true.
        ! se nell'input ?? presente <spectrum>
        !value =  getAttribute(xnode,"fcspec")
        !if (.not.isempty(value)) then
        !  read(value,*,iostat=ios) fc_job%fcspec
        !  if (ios /= 0) ierr = error(0,"Cannot read fcspec")
        !end if

        value =  getAttribute(root_xnode,"printfc")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%printfc
            if (ios /= 0) ierr = error(0,"Cannot read print")
        end if

        value =  getAttribute(root_xnode,"ftol")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%ftol
            if (ios /= 0) ierr = error(0,"Cannot read ftol")
        end if

        value =  getAttribute(root_xnode,"nclasses")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) fc_job%nclasses
            if (ios /= 0) ierr = error(0,"Cannot read nclasses")
        end if

        nsList => getElementsByTagName(root_xnode,"kubo")
        if (getLength(nsList) > 0) then
            fc_job%kubo%on = .true.
            myNode =>item(nsList,0)
            value = getAttribute(myNode,"file")
            if (.not.isempty(value)) fc_job%kubo%file = value
            value = getAttribute(myNode,"Temp")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%kubo%temp
                if (ios /= 0) ierr = error(0,"Cannot read Temp in <kubo>")
            end if
            value = getAttribute(myNode,"omr")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%kubo%omr
                if (ios /= 0) ierr = error(0,"Cannot read omr in <kubo>")
            end if
            value = getAttribute(myNode,"pow")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%kubo%pow
                if (ios /= 0) ierr = error(0,"Cannot read pow in <kubo>")
            end if
            value = getAttribute(myNode,"fwhm")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%kubo%fwhm
                if (ios /= 0) ierr = error(0,"Cannot read fwhm in <kubo>")
            end if
        end if

        ! Definisce il file dove salvare lo spettro e le caratteristiche dello
        ! spettro: FWHM, e forma (shape), che puo' essere lorentziana o gaussiana.
        nsList => getElementsByTagName(root_xnode,"spectrum")
        if (getLength(nsList) >= 1)  then
            fc_job%fcspec = .true.
            myNode =>item(nsList,0)
            value = getAttribute(myNode,"file")
            if (.not.isempty(value)) fc_job%spectrum%file = value
            value = getAttribute(myNode,"shape")
            if (.not.isempty(value)) fc_job%spectrum%shape = value
            value = getAttribute(myNode,"form")
            if (.not.isempty(value)) fc_job%spectrum%form= value
            value = getAttribute(myNode,"temp")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%spectrum%temp
                if (ios /= 0) ierr = error(0,"Cannot read temp in <spectrum>")
            end if
            value = getAttribute(myNode,"fwhm")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%spectrum%fwhm
                if (ios /= 0) ierr = error(0,"Cannot read fwhm in <spectrum>")
            end if
            value = getAttribute(myNode,"dE")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%spectrum%dE
                if (ios /= 0) ierr = error(0,"Cannot read dE in <spectrum>")
            end if
            value = getAttribute(myNode,"tol")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%spectrum%tol
                if (ios /= 0) ierr = error(0,"Cannot read tol in <spectrum>")
            end if
            value = getAttribute(myNode,"emin")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%spectrum%emin
                if (ios /= 0) ierr = error(0,"Cannot read emin in <spectrum>")
            end if
            value = getAttribute(myNode,"emax")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) fc_job%spectrum%emax
                if (ios /= 0) ierr = error(0,"Cannot read emax in <spectrum>")
            end if
        
            value = getAttribute(myNode,"type")
            if (.not.isempty(value)) fc_job%spectrum%type = value
        
        end if
    
        !------------------------------------------------------------------------------------------+
        ! Define <group> of vibrations which should be included in the FC calculation.             |
        ! It is possible to define several group and then take the tensor product of all the FCs   |
        ! as the overall FC integral.                                                              |
        !------------------------------------------------------------------------------------------+
        nsList => getElementsByTagName(root_xnode,"group")
        if (getLength(nsList) >= 1) then
            fc_job%ngroup = getLength(nsList)
        else
            ! il numero di gruppi deve essere almeno pari al numero di molecole
            ! (che ?? lo stesso per tutti gli stati elettronici).
            fc_job%ngroup = system%state(1)%nmol
        end if

        allocate(fc_job%group(1:fc_job%ngroup))

        ! Default setting: all vibrations included in each group, i.e. in each molecule.
        if (fc_job%ngroup == system%state(1)%nmol) then
            forall (ig=1:fc_job%ngroup) fc_job%group(ig)%nvib = system%state(1)%molecule%nvib
        end if

        ig = 0
NG:     do ! Start looping over <group>
		
            if (ig == fc_job%ngroup) exit
        
            ig = ig + 1
		        
            gxnode => item(nsList,ig-1) ! point to <group>
            if (associated(gxnode)) then
                value = getAttribute(gxnode,"id")
                read(value,*,iostat=ios) fc_job%group(ig)%id
                if (ios /= 0) then
                    ! automatically assign group id if missing
                    fc_job%group(ig)%id = ig
                end if
            !ierr = error(0,"Cannot read group id in <group>")
            end if
        
            ! Per ogni gruppo controlla se esiste una directive <include> e la analizza.

            ! <include> pu?? essere letto sia in <group> che nel nodo principale <fc>.
            ! serve per compatibilit?? con gli input vecchi. Si consiglia di mettere <include> sempre
            ! dentro <group>.
            if (associated(gxnode)) then
                npList => getElementsByTagName(gxnode,"include") ! cerca nel nodo <group>
            else
                npList => getElementsByTagName(root_xnode,"include") ! cerca nel nodo <fc>
            end if

            if (getLength(npList) == 1) then
                ! Ha trovato il nodo <include> e lo legge
                myNode => item(npList,0) ! point to <include> xnode
                ! the molecule defining group IG
                if (hasAttribute(myNode,"molecule")) then
                    fc_job%group(ig)%incvib%mol = getAttribute(myNode,"molecule")
                else
                    fc_job%group(ig)%incvib%mol = system%state(1)%molecule%id ! default molecule assignment
                end if
                value = getAttribute(myNode,"nvib")
                if (isempty(value)) then
                    ! include all vibrations
                    fc_job%group(ig)%nvib = system%state(1)%molecule%nvib
                else
                    read(value,*) fc_job%group(ig)%nvib
                end if
                allocate(fc_job%group(ig)%incvib%id(1:fc_job%group(ig)%nvib))
                ! get included vibrations
                np => getFirstChild(myNode)
                if (associated(np)) then
                    value = getNodeValue(np)
                    ! get included vibrations
                    fc_job%group(ig)%incvib%id(1:fc_job%group(ig)%nvib) = get_incvib(fc_job%group(ig)%nvib,value)
                else
                    forall (i=1:fc_job%group(ig)%nvib) fc_job%group(ig)%incvib%id(i) = i
                end if

            else if (getLength(npList) == 0) then
                ! Se il nodo <include> ?? assente allora si includono tutte le vibrazioni.
                fc_job%group(ig)%nvib = system%state(1)%molecule%nvib
                fc_job%group(ig)%incvib%mol = system%state(1)%molecule%id
                allocate(fc_job%group(ig)%incvib%id(1:fc_job%group(ig)%nvib))
                forall (i=1:fc_job%group(ig)%nvib) fc_job%group(ig)%incvib%id(i) = i
                !print *, 'All vibrations included'
            else if (getLength(npList) > 1) then
                ierr = error(0,"Only one <include> directive per <group> is allowed")
            end if

            ! Dopo aver letto <include> npList ora viene associata al tag <active>
            if (associated(gxnode)) then
                npList => getElementsByTagName(gxnode,"active") ! find <active> inside <group>
            else
                npList => getElementsByTagName(root_xnode,"active") ! find <active> inside <fc>
            end if
        
            allocate(fc_job%group(ig)%active(1:getLength(npList)),stat=alloc_err)
            if (alloc_err /= 0) ierr = error(0,"Cannot allocate memory for <active>.")

                    !----------------------------------------+
                    ! Cycle over <active> tag  in <group> ig |
                    !----------------------------------------+
                    ! If two or more <ative> tags with different molecules are found in the <fc> xnode then
                    ! the group number IG is changed accordingly.
                    ! This allow to use old (and simpler) input style without the definition of the <group> xnode.
                    ! TODO: Devo aggiungere un controllo su eventuali ripetizioni dei modi attivi in fase di input.
                    ! Il calcolo in tal caso va avanti senza accorgersi dell'errore.
AC:         do i = 0, getLength(npList) - 1

                myNode => item(npList,i)
                fc_job%group(ig)%active(i+1)%state = getAttribute(myNode,"state")

                if (isempty(fc_job%group(ig)%active(i+1)%state)) ierr = error(0,"Cannot read active state id")

                if (isempty(fc_job%group(ig)%incvib%mol)) then
                    fc_job%group(ig)%active(i+1)%molecule = getAttribute(myNode,"molecule")
                else
                    fc_job%group(ig)%active(i+1)%molecule = fc_job%group(ig)%incvib%mol
                end if
			
                if (isempty(fc_job%group(ig)%active(i+1)%molecule)) ierr = error(0,"Cannot read active molecule id")

                ! This option has not yet been implemented
                !value = getAttribute(myNode,"excitation")
                !if (.not.isempty(value)) then
                !    read(value,*) fc_job%group(ig)%active(i+1)%autex
                !    if (fc_job%group(ig)%active(i+1)%autex) then
                !        ! All vibrations excited
                !        allocate(fc_job%group(ig)%active(i+1)%mode(1:fc_job%group(ig)%nvib),stat=alloc_err)
                !        if (alloc_err /= 0) ierr = error(0,"Cannot allocate memory for excited vibration.")
                !    end if
                !end if

                ! List of active vibrations
                vibList => getElementsByTagName(myNode,"vibration")
                if (.not.allocated(fc_job%group(ig)%active(i+1)%mode)) then
                    ntvib = 0
                    do k = 0, getLength(vibList) - 1
                        myNode => item(vibList,k)
                        value = getAttribute(myNode,"id")
                        if (isempty(value)) then
                            np => getFirstChild(myNode)
                            if (associated(np)) then
                                value = getNodeValue(np)
                            end if
                        end if
                        if (.not.isempty(value)) then
                            exvib => get_numseq(value)
                            !if (ios /= 0) ierr = error(0,"Cannot read active vibration id sequence.")
                            ntvib = ntvib + size(exvib)
                        end if
                    end do
                    !allocate(fc_job%group(ig)%active(i+1)%mode(1:getLength(vibList)),stat=alloc_err)
                    allocate(fc_job%group(ig)%active(i+1)%mode(1:ntvib),stat=alloc_err)
                    if (alloc_err /= 0) ierr = error(0,"Cannot allocate memory for excited vibration.")
                end if

                fc_job%group(ig)%active(i+1)%nact = ntvib !getLength(vibList)

                is = get_state_from_id(system,fc_job%group(ig)%active(i+1)%state)
                !print *, 'state ', is,' molecule ',  im
            
                ! Le vibrazioni dello stato 1 vanno numerate da (NVIB + 1) a 2*NVIB dove
                ! NVIB e' il numero di vibrazioni della molecola IM.

                if (all((/FINAL_STATE, INITIAL_STATE/) /= is)) then
                    ierr = error(0,"Active state different from initial or finale states.")
                end if

                ioff = 0
VB:             do k = 0, getLength(vibList) - 1

                    myNode => item(vibList,k)

                    value = getAttribute(myNode,"id")
                    if (isempty(value)) then
                        np => getFirstChild(myNode)
                        if (associated(np)) then
                            value = getNodeValue(np)
                        end if
                    end if
                    if (.not.isempty(value)) then
                        exvib => get_numseq(value)
                      !read(value,*,iostat=ios) ivib
                      !if (ios /= 0) ierr = error(0,"Cannot read active vibration id")
                    else 
                        ierr = error(0,"Undefined active vibration id.")
                    end if

                    value = getAttribute(myNode,"nq")
                    nq = 1 ! default value of nq.
                    if (.not.isempty(value)) then
                        read(value,*,iostat=ios) nq
                        if (ios /= 0) ierr = error(0,"Cannot read vibration quanta (nq)")
                    end if

                    ! obsolete
                    !fc_job%active(i+1)%mode(k+1)%id = ivib + ioff

                    ! Forma nuova
                    fc_job%group(ig)%active(i+1)%mode(ioff+1:ioff+size(exvib))%id = exvib
                    fc_job%group(ig)%active(i+1)%mode(ioff+1:ioff+size(exvib))%nq = nq
                    fc_job%group(ig)%active(i+1)%mode(ioff+1:ioff+size(exvib))%freq = &
                    system%state(is)%molecule%normodes%vibration(exvib)%freq

                    ioff = ioff + size(exvib)
                ! Forma vecchia
                    !fc_job%group(ig)%active(i+1)%mode(k+1)%id = ivib
                    !fc_job%group(ig)%active(i+1)%mode(k+1)%nq = nq
                    !fc_job%group(ig)%active(i+1)%mode(k+1)%freq = system%state(is)%molecule%normodes%vibration(ivib)%freq

                ! obsolete
                    !system%state(is)%molecule%normodes%vibration(ivib)%nq = nq
                end do VB
                ! Check for duplicate entries.
                ! Controlla che non sia stato inserito due o pi?? volte lo stesso modo.
                do k = 1, getLength(vibList)
                    if (count(fc_job%group(ig)%active(i+1)%mode(:)%id == fc_job%group(ig)%active(i+1)%mode(k)%id) > 1) &
                    ierr = error(0,"Found duplicate entries in <active>!")
                end do
                ! Controlla che i modi attivi siano compatibili con quelli inclusi:
                do k = 1, getLength(vibList)
                    if (all(fc_job%group(ig)%incvib%id(:) /= fc_job%group(ig)%active(i+1)%mode(k)%id)) &
                    ierr = error(0,"Included modes and <active> modes do not match.")
                end do
            end do AC
        end do NG

        return
    end subroutine get_job_fc

    !==============================================================================

    subroutine get_job_transf (root_xnode, trns_job)

        type(Node), pointer, intent(in) ::  root_xnode
        type(Node), pointer ::  np, npc
        type(transf_t), intent(out) :: trns_job

        type(Node), pointer ::  cxnode
        type(NodeList), pointer ::  nList, nsList, ncList

        character(len=STRLEN) :: value, molid, alist, valuec
        integer ios, ierr, ic, intc, i, imol, istate, ilat, nintc, alloc_err, natc, ica, nintca
 
        ! This variables and parameters are used to assign internal coordinates (see
        ! below).
        character(len=4) :: ctype
        character(len=4), parameter :: ictype(8) = (/'s   ', 'b   ', 'l   ', 'w   ', 't   ', 'd   ', 'linc', 'linp'/)
        integer, parameter :: ntype(8) = (/2, 3, 4, 4, 2, 4, 4, 4/)
        real(kind=dp) :: coeff

        ! Get some attributes of <transformation>
        trns_job%molecule = getAttribute(root_xnode,"molecule")
        ! Set molecule to default molecule, i.e. the first of the first electronic state.
        if (isempty(trns_job%molecule)) trns_job%molecule = system%state(1)%molecule%id

        value = getAttribute(root_xnode,"tm_rotate")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) trns_job%tm_rotate
            if (ios /= 0) ierr = error(0,"Cannot read tm_rotate in <transformation>.")
        end if

        value = getAttribute(root_xnode,"type")
        if (.not.isempty(value)) then
            if (value == "model") then
                trns_job%model = .true.
                trns_job%axsw%on = .false. ! Do not perform axis switching
                trns_job%cartesian = .false. ! Do not perform cartesian transf.
                trns_job%internal = .false. ! Do not perform internal transf.
            end if
        end if

        !+------------------------------+
        ! Se usiamo un sistema modello  |
        !+------------------------------+
        if (trns_job%model) then
            ! Read normal modes rotations
            nsList => getElementsByTagName(root_xnode,"dusch")
            np => item(nsList,0)
            if (associated(np)) then
                trns_job%modr%file = getAttribute(np,"file")
                trns_job%modr%coord = getAttribute(np,"coord")
                np => getFirstChild(np)
                if (associated(np)) trns_job%modr%data = getNodeValue(np)
            end if
            ! Read normal modes shift
            nsList => getElementsByTagName(root_xnode,"shift")
            np => item(nsList,0)
            if (associated(np)) then
                trns_job%modk%file = getAttribute(np,"file")
                trns_job%modk%coord = getAttribute(np,"coord")
                trns_job%modk%type = getAttribute(np,"type")
                np => getFirstChild(np)
                if (associated(np)) trns_job%modk%data = getNodeValue(np)
            end if
            return
        end if
    
        trns_job%coord = getAttribute(root_xnode,"coord")
        ! Default transformation type is Cartesian.
        if (isempty(trns_job%coord)) trns_job%coord = "cartesian"

        if (trns_job%coord == "internal") then
            trns_job%cartesian= .false.
            trns_job%internal = .true.
        end if

        if (trns_job%coord == "natint") then
            trns_job%cartesian= .false.
            trns_job%internal = .false. ! optional
            trns_job%natint = .true.
        end if

        ! The next two options can be used if we want to modify internal coordinates after they have
        ! been generated.
        value = getAttribute(root_xnode,"icfile")
        if (.not.isempty(value)) then
            trns_job%setic = .true.
            trns_job%icfile = value
        end if

        value = getAttribute(root_xnode,"setic")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) trns_job%setic
            if (ios /= 0) ierr = error(0,"Cannot read SETIC in <transformation>.")
        end if

        if (trns_job%coord == "curv") then
            trns_job%cartesian= .false.
            trns_job%internal = .true.
            trns_job%nonlin   = .true.
        end if

        ! Read attribute for printing level (not yet used).
        value = getAttribute(root_xnode,"plevel")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) trns_job%printlevel
            if (ios /= 0) ierr = error(0,"Cannot read PLEVEL in <transformation>.")
        end if

        value = getAttribute(root_xnode,"debug")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) trns_job%debug
            if (ios /= 0) ierr = error(0,"Cannot read DEBUG in <transformation>.")
        end if

        value = getAttribute(root_xnode,"dusch")
        if (.not.isempty(value)) then
            read(value,*,iostat=ios) trns_job%dusch
            if (ios /= 0) ierr = error(0,"Cannot read DUSCH in <transformation>.")
        end if

        ! Possiamo specificare manualmente quale soluzione del problema di axis switching
        ! scegliamo.
        nList => getElementsByTagName(root_xnode,"axsw")
        if (getLength(nList) > 0) then

            trns_job%axsw%on = .true.

            np => item(nList,0)

            value = getAttribute(np,"on")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) trns_job%axsw%on
                if (ios /= 0) ierr = error(0,"Cannot read axsw method.")
            end if


            value = getAttribute(np,"solution")
            if (.not.isempty(value)) then
                read(value,*,iostat=ios) trns_job%axsw%solution
                if (ios /= 0) ierr = error(0,"Cannot set solution method of the axis switching procedure.")
            end if

        end if


        ! In questo nodo si seleziona l'utilizzo delle coordinate Out Of Plane Bengind
        ! e vengono definite manualmente le coordinate interne
        nList => getElementsByTagName(root_xnode,"internal_coordinates")

INTCIF: if (getLength(nList) > 0) then

            trns_job%cartesian= .false.
            trns_job%internal = .true.
            trns_job%intauto = .false.

INTC_DO:    do intc = 0, getLength(nList) - 1
           
                np => item(nList,0)

                value = getAttribute(np,"auto")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) trns_job%intauto
                    if (ios /= 0) ierr = error(0,"Cannot set automatic generation of redundant coordinates.")
                end if

                ! Out of plane bendings are not included by default. OPB must be set
                ! to true if we want to use them.
                value = getAttribute(np,"opb")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) trns_job%lopb
                    if (ios /= 0) ierr = error(0,"Out of plane bendings error. Check the OPB flag.")
                end if

                ! Read if we want to include dihedrals. Default is .true.
                value = getAttribute(np,"dih")
                if (.not.isempty(value)) then
                    read(value,*,iostat=ios) trns_job%ldih
                    if (ios /= 0) ierr = error(0,"Dihedral error. Check the DIH flag.")
                end if

                !+---------------------------------------------------------------------------------
                ! Here we start reading the internal coordinates if they have been manually defined
                !+---------------------------------------------------------------------------------
                ! Define list of internal coordinates tags
                nsList => getElementsByTagName(np,"coord")
                ! Number of internal coordinates
                nintc = getLength(nsList)

NINTC_IF:       if (nintc > 0) then
                  ! Allocate internal coordinates in each electronic state of the
                  ! system
                  do istate = 1, size(system%state)
                      allocate(system%state(istate)%molecule%intcoord%coord(1:nintc),stat=alloc_err)
                      if (alloc_err /= 0) ierr = error(71,error_string="coord() in intcoord")
                  end do
               
                    ! Read internal coordinates for the molecule specified by
                    ! trns_job%molecule
COORD_DO:         do ic = 1, getLength(nsList)
                    np => item(nsList,ic-1)
                    value = getAttribute(np,"type")
                    if (isempty(value)) ierr = error(0,error_string='Missing internal coordinate coordinate type')
                    select case (value)
                        case ("s","b","w","l","t","d","linp","linc")
                        continue
                    case default
                        ierr = error(0,"Wrong internal coordinate specification")
                    end select
                    
                    ! read the atom list
                    np => getFirstChild(np)
                    if (.not.associated(np)) exit
                    alist = getNodeValue(np)
               
                    do istate = 1, size(system%state)
                      ! this number (nintc) is used in other subroutine. I
                      ! could get rid of it by using: nintc = size(molecule(imol)%intcoord%coord)
                      system%state(istate)%molecule%intcoord%nintc = nintc
                      ! type e' un singolo carattere.
                      system%state(istate)%molecule%intcoord%coord(ic)%type = adjustl(value(1:1))
                      ctype = system%state(istate)%molecule%intcoord%coord(ic)%type
                      ctype = ToLower(ctype)
                      !-------------------------------------------------------!
                      ! In base al tipe di coordinata specificata alloca la   !
                      ! lista di atomi da cui e' compasta e l'array BMAT.     !
                      !-------------------------------------------------------!
                      do i = 1, size(ntype)
                          if(ictype(i) == ctype) ilat = ntype(i)
                      end do
                      !print *, ic
                      !print *, istate, imol, ilat
                      allocate(system%state(istate)%molecule%intcoord%coord(ic)%list(1:ilat),stat=alloc_err)
                      if (alloc_err /= 0) ierr = error(71,error_string="coord(intc).list() internal coordinates")
                      ! Read atom list for internal coord. definition.
                      read(alist,*) system%state(istate)%molecule%intcoord%coord(ic)%list(1:ilat)
                      !------------------------------------------------------------!
                      ! Se la coordinata e' una torsione bmat viene allocata nella !
                      ! subroutine intc.                                           !
                      !------------------------------------------------------------!
                      if (ctype /= 't') then
                            if (ctype == "linc" .or. ctype == "linp" .or. ctype == "l") then
                                allocate(system%state(istate)%molecule%intcoord%coord(ic)%bmat(1:3,1:3),stat=alloc_err)
                            else
                                allocate(system%state(istate)%molecule%intcoord%coord(ic)%bmat(1:ilat,1:3),stat=alloc_err)
                            end if
                            if (alloc_err /= 0) ierr = error(71,error_string="coord(intc).bmat(:,:)")
                            system%state(istate)%molecule%intcoord%coord(ic)%bmat(1:ilat,1:3) = zero
                      end if
                    end do
                  end do COORD_DO
                end if NINTC_IF

        end do INTC_DO
    end if INTCIF

        nList => getElementsByTagName(root_xnode,"natural_internal_coordinates")

NATINT_IF: if (getLength(nList) > 0) then

            print *, 'natint'
            trns_job%cartesian= .false.
            trns_job%internal = .false.
            trns_job%natint = .true.
            trns_job%intauto = .false.

                np => item(nList,0)
                !+-------------------------------------------------------------------------------------
                ! Here we start reading natural internal coordinates if they have been manually defined
                !+-------------------------------------------------------------------------------------
                ! Define list of internal coordinates tags
                nsList => getElementsByTagName(np,"ncoord")
                ! Number of natural internal coordinates
                natc = getLength(nsList)

print *, 'NATC ', natc
NATINTC_IF:     if (natc > 0) then
                    ! Allocate internal coordinates in each electronic state of the
                    ! system
                    do istate = 1, size(system%state)
                        allocate(system%state(istate)%molecule%intcoord%ncoord(1:natc),stat=alloc_err)
                        if (alloc_err /= 0) ierr = error(71,error_string="ncoord error")
                    end do
               

                    ! Read internal coordinates for the molecule specified by
                    ! trns_job%molecule
NCA_DO:             do ica = 1, getLength(nsList)
                      
                      !print *, 'NATINT ', ica
                      npc => item(nsList,ica-1)
                      ! Each <ncoord> is composed of several <coord>
                      ncList => getElementsByTagName(npc,"coord")
                      nintca = getLength(ncList)
                      ! Allocate ncoord()
                      do istate = 1, size(system%state)
                        allocate(system%state(istate)%molecule%intcoord%ncoord(ica)%coord(1:nintca),stat=alloc_err)
                        if (alloc_err /= 0) ierr = error(71,error_string="ncoord allocate error")
                      end do
 
                      value = getAttribute(npc,"type")
                      if (isempty(value)) ierr = error(0,error_string='Missing internal coordinate coordinate type')
                      select case (value)
                          case ("s","b","w","l","d","linc","linp")
                          continue
                      case default
                          ierr = error(0,"Wrong internal coordinate specification")
                      end select

                      system%state(1)%molecule%intcoord%ncoord(ica)%type = adjustl(value(1:1))
                      system%state(2)%molecule%intcoord%ncoord(ica)%type = adjustl(value(1:1))

NCB_DO:               do ic = 1, getLength(ncList)

                        np => item(ncList,ic-1)
                        value = getAttribute(np,"type")
                        if (isempty(value)) ierr = error(0,error_string='Missing internal coordinate type')
                        select case (value)
                            case ("s","b","w","l","d","linc","linp")
                            continue
                        case default
                            ierr = error(0,"Wrong internal coordinate specification")
                        end select
                       
                        valuec = getAttribute(np,"c")
                        coeff = one
                        if (.not.isempty(valuec)) then
                            read(valuec,*,iostat=ios) coeff
                            !print *, 'reading coeff : ', coeff
                        end if
 
                        ! Read atom list
                        np => getFirstChild(np)
                        if (.not.associated(np)) exit
                        alist = getNodeValue(np)
                       
                        do istate = 1, size(system%state)
                          ! this number (nintc) is used in other subroutines. I
                          ! could get rid of it by using: nintc = size(molecule(imol)%intcoord%coord)
                          system%state(istate)%molecule%intcoord%nintc = natc
                          ! type e' un singolo carattere.
                          system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%type = trim(adjustl(value(1:4)))
                          system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%c = coeff
                          ctype = system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%type
                          ctype = ToLower(ctype)
                          !-------------------------------------------------------!
                          ! In base al tipo di coordinata specificata alloca la   !
                          ! lista di atomi da cui e' composta  e l'array BMAT.    !
                          !-------------------------------------------------------!
                          do i = 1, size(ntype)
                              if(ictype(i) == ctype) ilat = ntype(i)
                          end do
                          !print *, ic
                          !print *, istate, imol, ilat
                          allocate(system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%list(1:ilat),stat=alloc_err)
                          if (alloc_err /= 0) ierr = error(71,error_string="ncoord(intc).list() internal coordinates")
                          ! Read atom list for internal coord. definition.
                          read(alist,*) system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%list(1:ilat)
                          ! Allocate BMAT
                          ! In the case of linear bending, even if we read in the position of 4 atoms the coordinate
                          ! only involves 3 atoms, thus we redefine ilat, and then we allocate the bmat for that coordinate.
                          !print *, 'ctype ', ctype
                          if (ctype == "linc" .or. ctype == "linp" .or. ctype == "l") ilat = 3
                          allocate(system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%bmat(1:ilat,1:3),stat=alloc_err)
                          if (alloc_err /= 0) ierr = error(71,error_string="coord(intc).bmat(:,:)")
                          system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%bmat(1:ilat,1:3) = zero
                        end do
                      end do NCB_DO
                   end do NCA_DO
              end if NATINTC_IF
    end if NATINT_IF

    return
end subroutine get_job_transf

!==============================================================================

function get_numseq(stringa) result(nsid)

    ! Exapand a numerical sequence of the type: 1-2,5,7,12-56
    ! in the whole sequence of numbers
    character(len=*) stringa
    character(len=1) sep, dash
    integer i, ia, iv, is, su, k, ins, inf, ierr
    integer, pointer :: nsid(:)

    !print *, stringa
    open(unit=tmpseq,status="scratch")

    sep = ','; dash = '-'

    is = 1; ia = 1; iv = 1; k = 0
    do while (iv /= 0)
        iv = index(stringa(ia:),sep)
        is = ia + iv
        if (ia == is) is = len_trim(stringa)
        su = index(stringa(ia:is),dash)
        if (su /= 0) then
            su = su + ia - 1
            read(stringa(ia:su-1),*) ins
            read(stringa(su+1:is),*) inf
        else
            read(stringa(ia:is),*) ins
            inf = ins
        end if
        if (inf < ins) ierr = error(0,"Cannot determine sequence.")
        do i = ins,inf 
            k = k + 1
            write(tmpseq,*) i
        end do
        ia = is
        ins = 0; inf = 0
        !print *, "B", si, is, string(si:)
    end  do

    rewind(tmpseq)
    allocate(nsid(1:k))
    read(tmpseq,*) (nsid(i),i=1,k)

    close(tmpseq) ! this also deletes the file

    do i = 1, size(nsid)
        if (count(nsid == nsid(i)) > 1) ierr = error(0,"Found duplicate entries in numerical sequence.")
    end do

end function get_numseq

!==============================================================================

function get_incvib(nvib,stringa) result(vibid)

    integer, intent(in) :: nvib
    character(len=*) stringa
    character(len=1) sep, dash
    integer i, ia, iv, is, su, k, ins, inf, ierr, icycle
    integer :: vibid(1:nvib)

    vibid = 0
    sep = ','; dash = '-'

    is = 1; iv = 1; k = 0
    icycle = 0
    do 
        if (len_trim(stringa) == 0) exit
        !print *, 'cycle ' , icycle
        iv = index(stringa,sep)
        if (iv == 0) then
            !print *,  'no coma '
            iv = len_trim(stringa)+1
        end if
        is = iv
        !print *, 'iv is', iv, is
        su = index(stringa(1:is-1),dash)
        if (su /= 0) then
            !print *, 'su ins', su, is, stringa(1:su-1)
            read(stringa(1:su-1),*) ins
            !print *, 'read ins', ins
            !print *, 'su inf', su, is
            read(stringa(su+1:is-1),*) inf
            !print *, '-', ins, inf
        else
            !print *, 'stringa ', stringa(1:is-1)
            read(stringa(1:is-1),*) ins
            inf = ins
            !print *, ',', ins, inf
        end if
        if (inf < ins) ierr = error(0,"Cannot determine included modes. Check <include> directive.")
        do i = ins,inf 
            k = k + 1
            if (k > size(vibid)) ierr = error(0,"Error in <include> directive. Too many vibrations.")
            vibid(k) = i
            !print *, 'VIB ', k, i
        end do
        if (is > 0) stringa = stringa(is+1:)
        !print *, 'str after ', stringa
        ins = 0; inf = 0
    end  do

    do i = 1, size(vibid)
        if (count(vibid == vibid(i)) > 1) then
            write(fout,*) 'Check entry: ', i
            write(fout,*), 'VIBID ', vibid
            ierr = error(0,"Found duplicate entries in <active> space and/or <subset>.")
        end if
    end do

end function get_incvib

!--------------------------------------------------------------------------!

subroutine read_tm(ntmList)
    !-------------------------!
    ! Read transition moments
    !-------------------------!
    type(NodeList), intent(in), pointer :: ntmList

    integer :: i, j, xnodeid, ierr
    real(kind=dp) :: muc
    character(len=1), parameter :: xyz(1:3) = (/'x', 'y', 'z'/)
    character(len=STRLEN) :: value

    type(Node), pointer :: xnode, cNode, muNode, dofNode
    type(NodeList), pointer :: nmuList, cList, dofList


    if (getLength(ntmList) == 1) then
        allocate (system%tm%dmudQ(1:3,1:system%state(1)%nvib))
        system%tm%dmudQ(1:3,1:system%state(1)%nvib) = zero
        xnode => item(ntmList,0)

        nmuList => getElementsByTagName(xnode,"mu0")
        do i = 1, getLength(nmuList)
            muNode => item(nmuList,i-1)
            do j = 1, 3
                cList => getElementsByTagName(muNode,xyz(j))
                if (getLength(cList) > 0) then
                    cNode => item(cList,0)
                    cNode => getFirstChild(cNode)
                    if (.not.associated(cNode)) ierr = error(0,"error in <tm><mu0><x|y|z>")
                    if (getNodeType(cNode) /= TEXT_NODE) ierr = error(0,"error in <tm><mu0><x|y|z>")
                    value = getNodeValue(cNode)
                    read(value,*) muc
                    system%tm%mu0(j) = muc
                end if
            end do
        end do

        nmuList => getElementsByTagName(xnode,"mu1")
        do i = 1, getLength(nmuList)
            muNode => item(nmuList,i-1)
            dofList => getElementsByTagName(muNode,"dof")
            dofNode => item(dofList,0)
            value = getAttribute(dofNode,"id")
            if (isempty(value)) ierr = error(0,"Error in tag <tm><dof>.")
            read(value,*) xnodeid
            do j = 1, 3
                cList => getElementsByTagName(muNode,xyz(j))
                if (getLength(cList) > 0) then
                    cNode => item(cList,0)
                    cNode => getFirstChild(cNode)
                    if (.not.associated(cNode)) ierr = error(0,"error in <tm><mu1><x|y|z>")
                    if (getNodeType(cNode) /= TEXT_NODE) ierr = error(0,"error in <tm><mu1><x|y|z>")
                    value = getNodeValue(cNode)
                    read(value,*) muc
                    !print *, xnodeid
                    !print *, muc
                    !print *, 'nvib ', system%state(1)%nvib
                    system%tm%dmudQ(j,xnodeid) = muc
                end if
            end do
        end do
        !print *, system%tm%dmudQ
    end if

end subroutine read_tm

end module input

