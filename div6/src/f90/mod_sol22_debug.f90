module sol22_debug
  use error_handling
  implicit none

  private


  !integer,public :: debug_sol22 = .false.
  !integer,public :: debug_sol22_ir, debug_sol22_ikopt



  logical,public :: debug_sol22_on = .false.
  integer :: debug_ir,debug_ikopt

  real*8, allocatable :: sol22_hr_data(:,:)

  integer :: init_size = 1000
  integer :: current_size
  integer :: data_counter= 0
  integer :: max_size = 100000
  integer :: n_data = 9


  public :: init_sol22_debug,check_init_record_data, save_s22_data,check_print_data




  contains

    subroutine init_sol22_debug(irdebug,ikoptdebug)
      implicit none
      integer :: irdebug, ikoptdebug

      debug_ir = irdebug
      debug_ikopt = ikoptdebug

      return

    end subroutine init_sol22_debug



    subroutine check_init_record_data(ir,ikopt)

      implicit none
      integer :: ir,ikopt

      if (debug_ir.eq.ir.and.debug_ikopt.eq.ikopt) then 
         ! Start S22 debug allocates the array and sets the debug flag
         call start_s22_debug
         data_counter =0
      endif

      return

    end subroutine check_init_record_data

    subroutine check_print_data(ir,ikopt)

      implicit none
      integer :: ir,ikopt

      !write(0,*) 'Checking print data:',ir,ikopt,debug_ir,debug_ikopt

      if (debug_ir.eq.ir.and.debug_ikopt.eq.ikopt) then 
         ! Print out the debugging data and turn off debugging
         ! 
         call print_s22_data(6)
      endif

      return

    end subroutine check_print_data



    subroutine start_s22_debug
      ! perform initial allocation of data storage array
      implicit none
      integer :: ierr

      debug_sol22_on = .true. 

      !write (0,*) 'Allocating sol22_hr_data:',init_size,n_data

      if (allocated(sol22_hr_data)) deallocate(sol22_hr_data)

      allocate(sol22_hr_data(init_size,n_data),stat=ierr)

      if (ierr.ne.0) call errmsg('START_S22_DEBUG: ALLOCATION ERROR:',ierr)

      current_size = init_size


    end subroutine start_s22_debug


    subroutine end_s22_debug
      ! clean up after S22 debug has been run

      implicit none
      if (allocated(sol22_hr_data)) deallocate(sol22_hr_data)

      data_counter= 0

    end subroutine end_s22_debug

    subroutine save_s22_data(h,s,ne,te,ti,vb,gamma,srci,srcf)
      implicit none
      real*8 :: h,s,ne,te,ti,vb,gamma,srci,srcf
      
      integer :: ierr

      data_counter = data_counter + 1

      if (data_counter.gt.max_size) then 
         ! array can not grow any further - exit routine
         return
      endif

      if (data_counter.gt.current_size) then 
         call grow_data_array
      endif

      sol22_hr_data(data_counter,1) = h
      sol22_hr_data(data_counter,2) = s
      sol22_hr_data(data_counter,3) = ne
      sol22_hr_data(data_counter,4) = te
      sol22_hr_data(data_counter,5) = ti
      sol22_hr_data(data_counter,6) = vb
      sol22_hr_data(data_counter,7) = gamma
      sol22_hr_data(data_counter,8) = srci
      sol22_hr_data(data_counter,9) = srcf

      return

    end subroutine save_s22_data

    subroutine grow_data_array
      implicit none

      integer :: next_size,ierr
      real*8,allocatable :: tmp_data(:,:)

      integer :: in,id
 
      ! If the array is as big as it can get - exit - should not be called in this case anyway
      if (current_size.eq.max_size) then 
         return
      endif

      if (.not.allocated(sol22_hr_data)) then
         ! real bug if this is called with sol22_hr_data not allocated
         return
      endif

      next_size = min(current_size*2, max_size)

      write(0,*) 'GROWING STORAGE:',next_size

      if (allocated(tmp_data)) deallocate(tmp_data)
      
      allocate(tmp_data(next_size,n_data),stat=ierr)

      if (ierr.ne.0) call errmsg('GROW_ARRAY: TMP_DATA ALLOCATION ERROR:',ierr)

      ! copy data so far into tmp_data

      do in=1,current_size
         do id = 1,n_data
            tmp_data(in,id) = sol22_hr_data(in,id)
         end do
      end do

      deallocate(sol22_hr_data)

      ! re-allocate sol22_hr_data with the new size

      allocate(sol22_hr_data(next_size,n_data),stat=ierr)

      if (ierr.ne.0) call errmsg('GROW_ARRAY: TMP_DATA ALLOCATION ERROR:',ierr)

      ! assign data to new array

      sol22_hr_data = tmp_data

      ! update current size of array
      current_size = next_size

      ! deallocate the temporary storage

      deallocate(tmp_data)


      return


    end subroutine grow_data_array

    subroutine print_s22_data(outunit)
      implicit none
      integer :: outunit

      integer :: in,id
      ! print out the data - usually to unit 6

      !write(0,*) 'Printing data to unit =',outunit

      write(outunit,'(a)') 'SOL22: HI RESOLUTION PLASMA PROFILE'

      write(outunit,'(15(1x,a14))') 'COUNT','STEP','S','NE','TE','TI','VB','GAMMA','INTIONSRC','INTSRCF','NE*VB'

      do in = 1,data_counter
         write(outunit,'(1x,i14,14(1x,g14.6))') in,(sol22_hr_data(in,id),id=1,n_data),sol22_hr_data(in,3)*sol22_hr_data(in,6)
      end do

      call end_s22_debug

      debug_sol22_on = .false.

      return

    end subroutine print_s22_data


  end module sol22_debug
  
