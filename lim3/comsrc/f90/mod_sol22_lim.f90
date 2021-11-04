module mod_sol22_lim


  implicit none

  private


  public :: sol22


contains


  subroutine sol22
    use mod_params
    use mod_comt2
    use mod_comtor
    use allocate_arrays
    use mod_calcsol_interface
    use yreflection
    use mod_solparams
    use mod_allocate_sol22_storage
    use mod_sol22_input_lim
    use mod_sol22_input
    use mod_sol22_utils
    use mod_comt2
    use mod_comxyt

    integer :: ir, ix, iy, in
    character*100 :: sol22_fname

    integer :: ierr
    integer :: ncnt

    real*8 :: nb0,te0,ti0,ringlen,r8cizb,r8crmb
    real*8,allocatable :: nb(:),te(:),ti(:),vb(:),spts(:)
    real :: cs
    
    integer :: pz

    integer :: new_unit
    integer :: iqx, ixout
    integer,external :: ipos
    real :: xbnd1,xbnd2,pbnd1,pbnd2

    call allocate_sol22_storage

    !call allocate_array(nb,nys/2,'nb',ierr)
    !call allocate_array(te,nys/2,'te',ierr)
    !call allocate_array(ti,nys/2,'ti',ierr)
    !call allocate_array(vb,nys/2,'vb',ierr)
    !call allocate_array(spts,nys/2,'spts',ierr)

    call allocate_array(nb,mxspts,'nb',ierr)
    call allocate_array(te,mxspts,'te',ierr)
    call allocate_array(ti,mxspts,'ti',ierr)
    call allocate_array(vb,mxspts,'vb',ierr)
    call allocate_array(spts,mxspts,'spts',ierr)

    write(0,*) 'SOL22A: start',nys,mxspts
    
    if (nys/2.gt.mxspts) then
       call errmsg('MOD_SOL22_LIM:','MXSPTS IS SMALLER THAN NYS/2 - INCREASE MXSPTS IN MOD_SOLCOMMON.F90')
       write(0,*) 'MXSPTS=',mxspts,' NYS/2=',nys/2
       write(0,*) 'LIM exiting'
       stop 'MOD_SOL22_LIM: MXSPTS TOO SMALL'
    endif

    !
    ! Find IX first outboard index

    IXOUT = IPOS (-1.E-10, XS, NXS-1)

    ! background plasma and geometry information
    r8cizb = cizb
    r8crmb = crmb
    ringlen = 2.0 * cl

    ! initialize pzone to 1 for now
    pz = 1 

    do ir = 1,nsol22_opt

       ! Initialize unstructured inputs
       call sol22_initialize_unstructured_input

       ! pick a unit number - reassign stdin to the new unit number
       ! open the file with the SOL22 option specifications

       call find_free_unit_number(new_unit)
       call set_unit_numbers(in_stdin=new_unit)
       open(file=trim(sol22_filenames(ir)),unit=new_unit,form='formatted',status='old')

       ! read input options for solver 
       call readsol(ierr)

       close(new_unit)
       ! reset the unit numbers
       call reset_unit_numbers

       xbnd1 = sol22_regions(ir,1)
       xbnd2 = sol22_regions(ir,2)

       ! ignore pbnds for now
       ! create two zones to overlay SOL22
       ! pzone = 1 -> |P| =< CPC0 (width of probe)
       ! pzone = 2 -> |P| > CPC0
       ! CPC0 
       ! if maxpzones = 1 OR not 3D (MAXNPS =1) then only
       ! calculate for pzone = 1
       ! 
       ! if there are absorbing surfaces - use them
       !    - calculate plasma over [yabsorb1a,yabsorb2a]
       !    - NOTE: there will need to be duplication of
       !      all plasma arrays for 1:pzone areas
       !      e.g. ctembs(ix,iz,pzone)
       !           ctegs ... and forces by iz and pzone
       !
       ! if not ignore
       !    - calculate plasma over whole range and duplicate
       !    - plasma from limiter surface to limiter surface
       !    - SOL22 is NOT applied to pzone 2 since there are
       !      no surfaces there.
       !    - This means this is only applicable for non-3D cases 
       !          
       !
       pbnd1 = sol22_regions(ir,3)
       pbnd2 = sol22_regions(ir,4)
       ! SOL22 input specifies which poloidal zone to use the model 
       pz = sol22_regions(ir,5)
       
       do ix = 1,nxs

          ! calculate y locations where the plasma needs to be calculated and
          ! map them to 0 ... ymax (needs to be mapped back afterward)
          ! Also need starting nb, Te, Ti values extracted from existing plasma
          ! or input. Take limits from SOLEDGE code. 

          ! ctembs(ix,iy) ... iy = -nys,nys 
          ! Two sides 
          ! Y < 0 = nys ... nys/2+1 OR -1 ... -nys

          ! only replace specified regions with SOL22 plasma 

          ! take qedges into account for Y coordinates
          ! qedges(iqx,1) for Y< 0 limiter surface
          ! qedges(iqx,2) for Y> 0 limiter surface
          ! HOWEVER - values are +ve for both so edge for Y<0 is -qedges(iqx,1)
          !

          if (xouts(ix).ge.xbnd1.and.xouts(ix).le.xbnd2) then 

             iqx = iqxs(ix)
             ncnt = 0

             !write(0,*) 'SOL22:',ix,iqx,ixout,xouts(ix),xbnd1,xbnd2


             ! only overlay SOL22 outboard when no absorbing surfaces are present
             if (yabsorb_opt.eq.0.and.ix.le.ixout) then 

                ncnt = 0
                do iy = nys,nys/2+1,-1
                   ! values with spts < 0 will be assigned the value of spts = 0
                   if ((2.0*cl - qedges(iqx,1) - youts(iy)) .gt. 0.0) then 
                      ncnt = ncnt +1
                      spts(ncnt) = 2.0 * cl - qedges(iqx,1) - youts(iy)
                   endif
                end do

                te0 = qtembs(iqx,1)  ! Y < 0 
                ti0 = qtembsi(iqx,1) ! Y < 0 
                nb0 = qrnbs(iqx,1)   ! Y < 0 
                cs  = getcs_sol22(sngl(te0),sngl(ti0))  ! positive flow toward Y<0 surface
                
                ! calculate SOL

                call calcsol_interface(nb0,te0,ti0,ringlen,ncnt,spts,nb,te,ti,vb,r8crmb,r8cizb)

                ! assign values back to arrays

                ncnt = 0
                do iy = nys,nys/2+1,-1
                   if ((2.0*cl - qedges(iqx,1) - youts(iy)) .gt. 0.0) then 
                      ncnt = ncnt +1
                      crnbs(ix,iy,pz) = nb(ncnt)
                      ctembs(ix,iy,pz) = te(ncnt)
                      ctembsi(ix,iy,pz) = ti(ncnt)
                      velplasma(ix,iy,pz) = -vb(ncnt)
                   else
                      crnbs(ix,iy,pz) = nb0
                      ctembs(ix,iy,pz) = te0
                      ctembsi(ix,iy,pz) = ti0
                      velplasma(ix,iy,pz) = cs
                   endif
                end do


                ! Y > 0 is 1...nys/2 or -nys...-nys/2-1

                ncnt = 0
                do iy = 1,nys/2
                   if ((youts(iy)-qedges(iqx,2)) .gt. 0.0) then 
                      ncnt = ncnt + 1
                      spts(ncnt) = youts(iy) - qedges(iqx,2)
                   endif
                end do


                te0 = qtembs(iqx,2)  ! Y > 0 
                ti0 = qtembsi(iqx,2) ! Y > 0 
                nb0 = qrnbs(iqx,2)   ! Y > 0 
                cs  = -getcs_sol22(sngl(te0),sngl(ti0))  ! negative flow toward Y>0 surface

                ! calculate SOL

                call calcsol_interface(nb0,te0,ti0,ringlen,ncnt,spts,nb,te,ti,vb,r8crmb,r8cizb)

                ncnt = 0
                do iy = 1,nys/2
                   if ((youts(iy)-qedges(iqx,2)) .gt. 0.0) then 
                      ncnt = ncnt +1
                      crnbs(ix,iy,pz) = nb(ncnt)
                      ctembs(ix,iy,pz) = te(ncnt)
                      ctembsi(ix,iy,pz) = ti(ncnt)
                      velplasma(ix,iy,pz) = vb(ncnt)
                   else
                      crnbs(ix,iy,pz) = nb0
                      ctembs(ix,iy,pz) = te0
                      ctembsi(ix,iy,pz) = ti0
                      velplasma(ix,iy,pz) = cs
                   endif

                end do


                ! copy over modified ctembs etc to -nys:-1

                do iy = 1,nys
                   crnbs(ix,-nys+iy-1,pz) = crnbs(ix,iy,pz)
                   ctembs(ix,-nys+iy-1,pz) = ctembs(ix,iy,pz)
                   ctembsi(ix,-nys+iy-1,pz) =  ctembsi(ix,iy,pz)
                   velplasma(ix,-nys+iy-1,pz) =  velplasma(ix,iy,pz)
                end do

                ! copy velplasma from pz = 1 to pz = 2
                !do in = 2,maxpzone
                !   velplasma(:,:,in) = velplasma(:,:,1)
                !end do
                
                


             elseif (yabsorb_opt.ge.2) then






             endif


             pz = 1
             ! print out results
             do iy = -nys,nys

                write(6,'(a,2i8,20(1x,g12.5))') 'SOL22:',ix,iy,iqxs(ix),qrnbs(iqx,1),qtembs(iqx,1),qtembs(iqx,1),&
                     qrnbs(iqx,2),qtembs(iqx,2),qtembs(iqx,2),youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),velplasma(ix,iy,pz)
             end do


          endif

       end do


    end do

    ! deallocate
    if (allocated(nb)) deallocate(nb)
    if (allocated(te)) deallocate(te)
    if (allocated(ti)) deallocate(ti)
    if (allocated(vb)) deallocate(vb)
    if (allocated(spts)) deallocate(spts)

    call deallocate_sol22_storage

  end subroutine sol22


end module mod_sol22_lim
