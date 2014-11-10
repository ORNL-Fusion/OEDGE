module allocate_arrays
  use error_handling
  implicit none


  interface allocate_array

     module procedure allocate_i_1d_array, allocate_i_2d_array,allocate_i_2db_array,&
                      allocate_r4_1d_array,allocate_r4_2d_array,allocate_r4_2db_array,&
                      allocate_r8_1d_array,allocate_r8_2d_array,allocate_r8_2db_array,&
                      allocate_r8_3d_array,allocate_r8_3db_array

  end interface


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
    endif

end subroutine allocate_i_1d_array

  subroutine allocate_i_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    integer,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
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
    endif

end subroutine allocate_i_2db_array


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
    endif

end subroutine allocate_r4_1d_array

  subroutine allocate_r4_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    real,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
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
    endif

end subroutine allocate_r4_2db_array


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
    endif

end subroutine allocate_r8_1d_array

  subroutine allocate_r8_2d_array(array,dim1,dim2,desc,ierr)
    implicit none
    character*(*) :: desc
    integer :: dim1,dim2,ierr
    real*8,allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)
    allocate(array(dim1,dim2),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim(desc)//' IERR =',ierr)
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
    endif

end subroutine allocate_r8_3db_array




end module allocate_arrays
