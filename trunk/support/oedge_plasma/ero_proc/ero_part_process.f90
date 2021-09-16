module ero_part_process
  use error_handling
  use utilities

  integer,parameter :: nvert=4
  real*8 :: vert(nvert,2)
  real*8 :: roffset,zoffset
  real*8 :: tor_extent

  real*8 :: xhat(2),yhat(2) ! vector coordinate transformation used both in the plasma and geometry code - calculated in plasma code
  real*8 :: xhatinv(2),yhatinv(2) ! vector coordinate transformation used both in the plasma and geometry code - calculated in plasma code

  integer, parameter :: line_length = 512
  character*6, parameter :: line_form='(a512)'
  character*(line_length) :: buffer,ero_part_data_fn,erodiv_part_data_fn

public read_transform_data,process_ero_part_data


contains


  subroutine read_transform_data
    implicit none
    character*100 :: coordinate_filename
    logical :: done
    integer :: inunit,i,ios

    coordinate_filename ='coordinate_transform_data.txt'

    ! this routine writes out the coordinate transformation data so that 
    ! the particle processing routing can run it in reverse to obtain the 
    ! DIVIMP coordinates from the ERO particle coordinates. 

    call find_free_unit_number(inunit)

    open(inunit,file=trim(coordinate_filename),form='formatted',iostat=ios)

    ! Need to read in
    ! 1) The vertices
    ! 2) The transform vectors
    ! 3) The inverse transform vectors 
    ! 4) The offsets

    done = .false.

    do while (.not.done)  

       read(inunit,line_form,iostat=ios) buffer

       !write(0,'(a,a,a,i5)') 'Buff:',trim(buffer),':',ios

       if (ios.ne.0) then 
          ! need to exit loop if eof reached - not an error since could be last block
          done = .true. 
          cycle
       endif

       if (buffer(1:8).eq.'#OFFSETS') then 
          ! read R,Z coordinate offset 
          read(buffer(9:),*) roffset,zoffset
          !write(0,'(a,10(1x,g12.5))') 'offset:',roffset,zoffset

       elseif (buffer(1:10).eq.'#RVERTICES') then 
          ! read R vertices
          read(buffer(11:),*) (vert(i,1),i=1,4)
          !write(0,'(a,10(1x,g12.5))') 'rvert:',(vert(i,1),i=1,4)

       elseif (buffer(1:10).eq.'#ZVERTICES') then 
          ! read Z vertices
          read(buffer(11:),*) (vert(i,2),i=1,4)
          !write(0,'(a,10(1x,g12.5))') 'zvert:',(vert(i,2),i=1,4)

       elseif (buffer(1:11).eq.'#TOR_EXTENT') then 
          ! read toroidal extent
          read(buffer(12:),*) tor_extent
          !write(0,'(a,10(1x,g12.5))') 'tor_extent:',tor_extent

       elseif (buffer(1:5).eq.'#XHAT') then 

          read(buffer(6:),*) (xhat(i),i=1,2)
          !write(0,'(a,10(1x,g12.5))') 'xhat:',(xhat(i),i=1,2)

       elseif (buffer(1:5).eq.'#YHAT') then 

          read(buffer(6:),*) (yhat(i),i=1,2)
          !write(0,'(a,10(1x,g12.5))') 'yhat:',(yhat(i),i=1,2)

       elseif (buffer(1:5).eq.'#XINV') then 

          read(buffer(6:),*) (xhatinv(i),i=1,2)
          !write(0,'(a,10(1x,g12.5))') 'xhatinv:',(xhatinv(i),i=1,2)

       elseif (buffer(1:5).eq.'#YINV') then 

          read(buffer(6:),*) (yhatinv(i),i=1,2)
          !write(0,'(a,10(1x,g12.5))') 'yhatinv:',(yhatinv(i),i=1,2)

       elseif (buffer(1:10).eq.'#EROPARTFN') then 

          read(buffer(11:),*) ero_part_data_fn
          !write(0,'(a,a,a,a)') 'ero part output:',trim(ero_part_data_fn),':'

       elseif (buffer(1:13).eq.'#ERODIVPARTFN') then 

          read(buffer(14:),*) erodiv_part_data_fn
          !write(0,'(a,a,a,a)') 'ero->div part input:',trim(erodiv_part_data_fn),':'

       endif

    end do

    close(inunit)


  end subroutine read_transform_data


  subroutine process_ero_part_data
    ! rewrite the ERO particle data transformed to DIVIMP coordinates
    implicit none

    ! Open the input and output files
    ! Transform coordinates (and units if required) 
    ! Write to output file

    integer :: partin,partout,ios,in
    logical :: done
    integer :: part_count
    real*8 :: total_atoms
    real*8 :: weight
    real*8 :: cumulative_weight

    real*8 :: rprime, zprime, vrprime, vzprime
    real*8 :: tprime

    real*8 :: x,y,z,vx,vy,vz,v0,energy,atoms
    integer :: charge,lost,pnum


    ! open files

    call find_free_unit_number(partin)

    open(partin,file=trim(ero_part_data_fn),status='old',form='formatted',iostat=ios)

    if (ios.ne.0) then 
       call errmsg('ERROR opening ERO particle input file: Case = ',trim(ero_part_data_fn))
       return
    endif

    call find_free_unit_number(partout)

    open(partout,file=trim(erodiv_part_data_fn),form='formatted',iostat=ios)

    if (ios.ne.0) then 
       call errmsg('ERROR opening DIVIMP-> particle transfer file: Case = ',trim(erodiv_part_data_fn))
       return
    endif


    part_count = 0
    total_atoms = 0.0
    done = .false.

    !
    ! First scan the file and count particles and total atoms
    !

    do while (.not.done)

       read(partin,line_form,iostat=ios) buffer

       if (ios.ne.0) then 
          ! need to exit loop if eof reached 
          done = .true. 
          cycle
       endif

       if (buffer(1:1).ne.'#'.and.buffer(1:1).ne.'%') then 

          part_count = part_count+1
          read(buffer,*) pnum,lost,charge,x,y,z,vx,vy,vz,v0,energy,atoms

          total_atoms = total_atoms+atoms

       endif
    end do

    ! Now process the file

    rewind(partin)

    ! convert particle source strength to particles/s/m-tor by dividing by the toroidal extent of the 
    ! ERO simulation space

    write(partout,'(a,(1x,i12))') '#NPARTICLES',part_count
    write(partout,'(a,2(1x,g18.8))') '#TOTAL_ATOMS', total_atoms,total_atoms/(2.0*tor_extent)
    write(partout,'(a,1x,g18.8)') '#FULL_TOR_EXTENT',2.0*tor_extent

    part_count = 0
    done =.false.

    do while (.not.done)

       read(partin,line_form,iostat=ios) buffer

       if (ios.ne.0) then 
          ! need to exit loop if eof reached 
          done = .true. 
          exit
       endif

       if (buffer(1:1).ne.'#'.and.buffer(1:1).ne.'%') then 

          part_count = part_count + 1

        ! tor,pol,radial, vx, vy,vz, energy, charge,atoms, surface_lost
        !    fprintf(erodiv_out,"%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %d %12.6e %d\n",
	!	      work_ptr->location[0]/1000.0,work_ptr->location[1]/1000.0,work_ptr->location[2]/1000.0,     /* jdecmt: Convert to m from mm */
 	!	      work_ptr->v[1]/100.0,work_ptr->v[2]/100.0,work_ptr->v[3]/100.0,work_ptr->v[0]/100.0,        /*         Convert from cm/s to m/s */
	!              work_ptr->energy,work_ptr->charge,work_ptr->former_atoms,work_ptr->reason_lost);


          read(buffer,*) pnum,lost,charge,x,y,z,vx,vy,vz,v0,energy,atoms

          ! convert units - distance mm to m
          x = x/1000.0
          y = y/1000.0
          z = z/1000.0

          ! convert units - velocity cm/s -> m/s
          vx = vx/100.0
          vy = vy/100.0
          vz = vz/100.0
          v0 = v0/100.0


          ! X->T  Y -> Z  Z-> R
          rprime = roffset - (z * xhatinv(1) + y * xhatinv(2)) 
          zprime = zoffset - (z * yhatinv(1) + y * yhatinv(2)) 

          vrprime = vz * xhat(1) + vy * xhat(2)
          vzprime = vz * yhat(1) + vy * yhat(2)
          
          weight = atoms/total_atoms

          cumulative_weight=cumulative_weight+weight

          ! convert individual particle atom counts to number/s/m-tor
          atoms = atoms/(2.0*tor_extent)

             ! T->X  Z -> Y  R-> Z
          write(partout,'(i8,i8,10(1x,g18.8),1x,g20.12)') part_count,charge,rprime,zprime,x, &
                  & vrprime,vzprime,vx,v0,energy,atoms,weight,cumulative_weight


       endif

    end do

    ! close files
    close(partin)
    close(partout)

  end subroutine process_ero_part_data


end module ero_part_process
