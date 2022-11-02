module mod_readgrid

  implicit none



contains


  subroutine get_grid_parameters(gridunit,nr,ir_cut,nk,ik_cut1,ik_cut2)
    implicit none

    integer :: gridunit,nr,nk,ir_cut,ik_cut1,ik_cut2


    integer :: ios
    character*1000 :: buffer, line
    integer :: in,ik,ir
    real :: rv(4),zv(4),rc,zc,brat,psin

    ! read cells from the grid record the in, ik, ir for the first and last cells found
    !              1    4
    ! checking that 3,4 of last = 1,2 of current to find cuts
    !
    ! grid contains boundary cells
    !
    ! first ring without cuts is the separatrix
    !
    real :: rv_last(4), zv_last(4)
    integer:: in_start,ik_start,ir_start,in_end,ir_end,ik_end
    integer:: last_ik_cut1,last_ik_cut2 
    integer :: ir_current,np,ix
    logical :: first

    ios = 0

    ! read until line of '====='
    !gridunit = 10
    !call find_free_unit_number(gridunit)

    ! open (gridunit,file='sasv-190422-6-v002.div',status='old')

    rewind(gridunit)

    read(gridunit,'(a)') buffer

    do while (buffer(4:8).ne.'=====') 

       read(gridunit,'(a)',iostat=ios) buffer

       if (ios.ne.0) then
          write(0,*) 'Error reading grid file',ios
          stop 'Error reading grid file'
       endif
    end do

    ! at start of element listing

    ir_cut = 0
    ik_end = 0
    ir_end = 0
    in_end = 0
    last_ik_cut1 = 0
    last_ik_cut2 = 0
    ik_cut1 = 0
    ik_cut2 = 0

    first = .true.

    do while (ios.eq.0)

       call read_cell(gridunit,buffer,line,in,ik,ir,rv,zv,rc,zc,brat,psin,ios)

       if (ios.ne.0) exit

       !write(0,*) 'Cell:',in,ik,ir,ios,(rv(ix),zv(ix),ix=1,4)
       write(6,*) 'Cell:',in,ik,ir,ios,(rv(ix),zv(ix),ix=1,4)

       if (first) then
          in_start =  in
          ik_start =  ik
          ir_start =  ir
          ir_current = ir_start
          first = .false.
          rv_last = rv
          zv_last = zv
       else
          if (ir_current.ne.ir) then
             ! changing rings - don't check for cuts

             if (ik_cut1.eq.0.and.ik_cut2.eq.0.and.ir_cut.eq.0) then
                ! if the previous ring had no cuts but ir_cut has not been set then
                ! set ir_cut to the last ring before updating ir_current
                ir_cut = ir_current
             endif

             if (ir_cut.eq.0) then
                ! check to see if ik_cuts are consistent 
                if (last_ik_cut1.ne.0.and.last_ik_cut2.ne.0) then
                   if (last_ik_cut1.ne.ik_cut1.or.last_ik_cut2.ne.ik_cut2) then
                      write(0,*) 'ERROR: Cut point different on consecutive rings:',ir_current,last_ik_cut1,ik_cut1,last_ik_cut2,ik_cut2
                      write(6,*) 'ERROR: Cut point different on consecutive rings:',ir_current,last_ik_cut1,ik_cut1,last_ik_cut2,ik_cut2
                   endif
                endif

                if (ik_cut1.ne.0.and.ik_cut2.eq.0) then
                   write(0,*) 'ERROR: Only one cut point found on ring:',ir_current,ik_cut1,ik_cut2
                   write(6,*) 'ERROR: Only one cut point found on ring:',ir_current,ik_cut1,ik_cut2
                endif
             endif

             ir_current = ir
             rv_last = rv
             zv_last = zv
             if (ik_cut1.ne.0) then 
                last_ik_cut1 = ik_cut1
             endif
             if (ik_cut2.ne.0) then 
                last_ik_cut2 = ik_cut2
             endif
             ik_cut1 = 0
             ik_cut2 = 0

          else
             ! check to see if any vertex does not match - indicates a cut
             if (.not.((rv_last(3).eq.rv(2).and.zv_last(3).eq.zv(2)).and.&
                  ((rv_last(4).eq.rv(1).and.zv_last(4).eq.zv(1))))) then
                !write(0,*) 'Matches:',rv_last(3)-rv(2),zv_last(3)-zv(2),rv_last(4)-rv(1),zv_last(4)-zv(1),ik_cut1,ik_cut2,ik,ir
                if (ik_cut1.eq.0) then
                   ik_cut1 = ik
                elseif (ik_cut2.eq.0) then
                   ik_cut2 = ik
                else
                   write(0,*) 'ERROR: Extra cuts found on ring:',ir_current,ik_cut1,ik_cut2,ik
                   write(6,*) 'ERROR: Extra cuts found on ring:',ir_current,ik_cut1,ik_cut2,ik
                endif
             endif

             rv_last = rv
             zv_last = zv

          endif
       endif

       !if (ios.ne.0) exit

       ! save last indices
       in_end = in
       ik_end = ik
       ir_end = ir

    end do

    nr = ir_end-ir_start+1
    nk = ik_end-ik_start+1
    np = in_end-in_start+1

    ir_cut = ir_cut-ir_start
    ik_cut1 = last_ik_cut1-ik_start
    ik_cut2 = last_ik_cut2-ik_start+1  ! shift by 1 to point to cell after the cut

    rewind(gridunit)


  end subroutine get_grid_parameters



  subroutine read_cell(gridunit,buffer,line,in,ik,ir,rv,zv,rc,zc,brat,psin,ios)
    implicit none
    character*(*) buffer,line
    integer :: in,ik,ir,ios,gridunit
    real rv(4),zv(4)
    real rc,zc,psin
    real brat
    !
    !     This routine reads in one cell from a sonnet grid formatted file.
    !     This routine was added because Carre grid output format was changed
    !     and it is desirable that DIVIMP be able to continue reading these
    !     grids with at least some flexibility as to format.       
    !
    !     locals
    integer :: p1,p2
    !     
    ios = 0


    ! first line

    read(gridunit,'(a)',iostat=ios) line

    !write(0,*) 'line1:',trim(line),':'


    if (ios.ne.0) then
       ! set error code to 2 if error encountered reading first line of cell data
       ios = 2
       return
    endif

    p1 = index(line,'Element') +len_trim('Element') -1
    
    if (p1.lt.len_trim('Element')) then
       ! No more elements found - exit
       ios = 4
       return
    endif

    p2 = index(line(p1+1:),'=') + p1
    !      write(0,*) trim(line)
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >               line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) in

    p1 = index(line(p2+1:),'(') + p2 
    p2 = index(line(p1+1:),',') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >               line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) ik

    p1 = p2 
    p2 = index(line(p1+1:),')') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >               line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) ir

    p1 = index(line(p2+1:),'(') + p2
    p2 = index(line(p1+1:),',') + p1
    !     write(0,*) p1,':',p2,':',line(p1:p1),':',
    !    >               line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) rv(2)

    p1 = p2
    p2 = index(line(p1+1:),')') + p1 
    !     write(0,*) p1,':',p2,':',line(p1:p1),':',
    !    >               line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) zv(2)

    p1 = index(line(p2+1:),'(') +p2
    p2 = index(line(p1+1:),',') +p1
    !     write(0,*) p1,':',p2,':',line(p1:p1),':',
    !    >               line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) rv(3)

    p1 = p2
    p2 = index(line(p1+1:),')') + p1
    !     write(0,*) p1,':',p2,':',line(p1:p1),':',
    !    >               line(p2:p2),':',line(p1+1:p2-1),':'

    read(line(p1+1:p2-1),*) zv(3)
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !    >           line(p2:p2),':',line(p1+1:p2-1),':'

    !      write(0,*) 'cell:',in,ik,ir,rv(2),zv(2),rv(3),zv(3)


    ! second line
    read(gridunit,'(a)',iostat=ios) line

    !      write(0,*) 'line2:',trim(line),':'

    if (ios.ne.0) then
       ! set error code to 1 if error encountered reading first line of cell data
       ios = 3
       return
    endif

    p1 = index(line,'=') 
    p2 = index(line(p1+1:),'(') +p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) brat

    p1 = p2
    p2 = index(line(p1+1:),',') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) rc

    p1 = p2
    p2 = index(line(p1+1:),')') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) zc

    line = trim(line) // '    0.0000' ! put a 0.000 psi value at the end of the line for old meshes that did not have this 

    p1 = p2
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:),*) psin

    ! third line

    read(gridunit,'(a)',iostat=ios) line

    !      write(0,*) 'line3:',trim(line),':'

    if (ios.ne.0) then
       ! set error code to 1 if error encountered reading first line of cell data
       ios = 3
       return
    endif

    p1 = index(line,'(') 
    p2 = index(line(p1+1:),',') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'

    read(line(p1+1:p2-1),*) rv(1)

    p1 = p2
    p2 = index(line(p1+1:),')') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'


    read(line(p1+1:p2-1),*) zv(1)

    p1 = index(line(p2+1:),'(') + p2
    p2 = index(line(p1+1:),',') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'
    read(line(p1+1:p2-1),*) rv(4)

    p1 = p2
    p2 = index(line(p1+1:),')') + p1
    !      write(0,*) p1,':',p2,':',line(p1:p1),':',
    !     >           line(p2:p2),':',line(p1+1:p2-1),':'

    read(line(p1+1:p2-1),*) zv(4)

    ! read expected separator into buffer 

    read(gridunit,'(a)',iostat=ios) buffer

    !      write(0,*) 'line4:',trim(buffer),':'

    if (ios.ne.0) then
       ! set error code to 1 if error encountered reading first line of cell data
       ios = 3
       return
    endif


    !      read(line,9000) in,ik,ir,rvert(2),zvert(2),rvert(3),zvert(3)
    !     
    !      read(gridunit,'(a)',end=2000) line
    !      read(line,9001) brat,rcent,zcent
    !
    !      jdemod - read in the psi value for the cell if present
    !      
    !      line = trim(line)//'                         0.0   '
    !      read(line(87:),*) psin_cent
    !
    !      read(gridunit,'(a)',end=2000) line
    !      read(line,9002) rvert(1),zvert(1),rvert(4),zvert(4)
    !      read(gridunit,'(a)',end=2000) buffer


  end subroutine read_cell




end module mod_readgrid
