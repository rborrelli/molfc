MODULE natint

  SUBROUTINE read_natural_internal()



  END SUBROUTINE read_natural_internal

NATINTC_IF:     if (natc > 0) then
                    ! Allocate internal coordinates in each electronic state of the
                    ! system
                    do istate = 1, size(system%state)
                        allocate(system%state(istate)%molecule%intcoord%ncoord(1:natc),stat=alloc_err)
                        if (alloc_err /= 0) ierr = error(71,error_string="ncoord error")
                    end do
               

                    ! Read internal coordinates for the molecule specified by
                    ! trns_job%molecule
NCA_DO:          do ica = 1, getLength(nsList)
                      
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
                          case ("s","b","w","l","d")
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
                            case ("s","b","w","l","d")
                            continue
                        case default
                            ierr = error(0,"Wrong internal coordinate specification")
                        end select
                       
                        valuec = getAttribute(np,"c")
                        coeff = one
                        if (.not.isempty(valuec)) then
                            read(valuec,*,iostat=ios) coeff
                            !print *, 'reading coeff : ', coeff, value
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
                          system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%type = adjustl(value(1:1))
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
                          allocate(system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%bmat(1:ilat,1:3),stat=alloc_err)
                          if (alloc_err /= 0) ierr = error(71,error_string="coord(intc).bmat(:,:)")
                          system%state(istate)%molecule%intcoord%ncoord(ica)%coord(ic)%bmat(1:ilat,1:3) = zero
                        end do
                      end do NCB_DO
                   end do NCA_DO
              end if NATINTC_IF
    end if NATINT_IF

END MODULE natint

