module allocate_arrays
  use error_handling
  implicit none


  ! note: the 1d array with dual bounds specified would have the same signature as the 2d array routine
  !       - to avoid a conflict in signatures the extra dimention is added after the description even though it 
  !         won't make obvious sense. It will disambiguate the routine signatures in the interface

  interface allocate_array

     module procedure allocate_i_1d_array, allocate_i_1db_array,&
          allocate_i_2d_array,allocate_i_2db_array,&
          allocate_i_3d_array,allocate_i_3db_array,&
          allocate_i_4db_array,allocate_i_5db_array,&

          allocate_l_1d_array, allocate_l_1db_array,&
          allocate_l_2d_array,allocate_l_2db_array,&
          allocate_l_3d_array,allocate_l_3db_array,&
          allocate_l_4db_array,allocate_l_5db_array,&

          allocate_r4_1d_array,allocate_r4_1db_array,&
          allocate_r4_2d_array,allocate_r4_2db_array,&
          allocate_r4_3d_array,allocate_r4_3db_array,&
          allocate_r4_4db_array,allocate_r4_5db_array,&

          allocate_r8_1d_array,allocate_r8_1db_array,&
          allocate_r8_2d_array,allocate_r8_2db_array,&
          allocate_r8_3d_array,allocate_r8_3db_array,&
          allocate_r8_4db_array,allocate_r8_5db_array,&

          allocate_c_2d_array

  end interface allocate_array


contains

  ! Integer

  subroutine allocate_i_1d_array(array,dim1,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,ierr
    integer,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_1d_array

  subroutine allocate_i_1db_array(array,dim1a,desc,dim1b,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,ierr
    integer,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_1db_array


  subroutine allocate_i_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    integer,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_2d_array

  subroutine allocate_i_2db_array(array,dim1a,dim1b,dim2a,dim2b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,ierr
    integer,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_2db_array

  subroutine allocate_i_3d_array(array,dim1,dim2,dim3,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,dim3,ierr
    integer,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2,dim3),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_3d_array

  subroutine allocate_i_3db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,ierr
    integer,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_3db_array

  subroutine allocate_i_4db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,ierr
    integer,allocatable :: array(:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_4db_array

  subroutine allocate_i_5db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,ierr
    integer,allocatable :: array(:,:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b,dim5a:dim5b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0
    endif

  end subroutine allocate_i_5db_array


  ! Logical

  subroutine allocate_l_1d_array(array,dim1,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,ierr
    logical,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_1d_array

  subroutine allocate_l_1db_array(array,dim1a,desc,dim1b,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,ierr
    logical,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_1db_array


  subroutine allocate_l_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    logical,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_2d_array

  subroutine allocate_l_2db_array(array,dim1a,dim1b,dim2a,dim2b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,ierr
    logical,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_2db_array

  subroutine allocate_l_3d_array(array,dim1,dim2,dim3,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,dim3,ierr
    logical,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2,dim3),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_3d_array

  subroutine allocate_l_3db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,ierr
    logical,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_3db_array

  subroutine allocate_l_4db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,ierr
    logical,allocatable :: array(:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_4db_array

  subroutine allocate_l_5db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,ierr
    logical,allocatable :: array(:,:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b,dim5a:dim5b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = .false.
    endif

  end subroutine allocate_l_5db_array


  ! Real*4

  subroutine allocate_r4_1d_array(array,dim1,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,ierr
    real,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_1d_array

  subroutine allocate_r4_1db_array(array,dim1a,desc,dim1b,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,ierr
    real,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_1db_array


  subroutine allocate_r4_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    real,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_2d_array

  subroutine allocate_r4_2db_array(array,dim1a,dim1b,dim2a,dim2b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,ierr
    real,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_2db_array


  subroutine allocate_r4_3d_array(array,dim1,dim2,dim3,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,dim3,ierr
    real*4,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2,dim3),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_3d_array

  subroutine allocate_r4_3db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,ierr
    real*4,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_3db_array

  subroutine allocate_r4_4db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,ierr
    real*4,allocatable :: array(:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_4db_array

  subroutine allocate_r4_5db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,ierr
    real*4,allocatable :: array(:,:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b,dim5a:dim5b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0
    endif

  end subroutine allocate_r4_5db_array



  ! Real*8 

  subroutine allocate_r8_1d_array(array,dim1,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,ierr
    real*8,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_1d_array

  subroutine allocate_r8_1db_array(array,dim1a,desc,dim1b,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,ierr
    real*8,allocatable :: array(:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_1db_array

  subroutine allocate_r8_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    real*8,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_2d_array

  subroutine allocate_r8_2db_array(array,dim1a,dim1b,dim2a,dim2b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,ierr
    real*8,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_2db_array

  subroutine allocate_r8_3d_array(array,dim1,dim2,dim3,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,dim3,ierr
    real*8,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2,dim3),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_3d_array

  subroutine allocate_r8_3db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,ierr
    real*8,allocatable :: array(:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_3db_array

  subroutine allocate_r8_4db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,ierr
    real*8,allocatable :: array(:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_4db_array

  subroutine allocate_r8_5db_array(array,dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1a,dim1b,dim2a,dim2b,dim3a,dim3b,dim4a,dim4b,dim5a,dim5b,ierr
    real*8,allocatable :: array(:,:,:,:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1a:dim1b,dim2a:dim2b,dim3a:dim3b,dim4a:dim4b,dim5a:dim5b),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = 0.0d0
    endif

  end subroutine allocate_r8_5db_array

  subroutine allocate_c_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    character*(*),allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
    else
       array = ''
    endif

  end subroutine allocate_c_2d_array


end module allocate_arrays
