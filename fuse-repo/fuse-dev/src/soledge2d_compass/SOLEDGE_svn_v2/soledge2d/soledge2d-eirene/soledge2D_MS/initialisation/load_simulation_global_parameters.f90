subroutine load_simulation_global_parameters()
  use all_variables, only : global_variables
  use hdf5
  implicit none
  !HDF5 variable
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dspace_id
  integer :: error
  character(50) :: filename
  integer(HSIZE_T) :: dim1d(1)
  !other variable
  integer*4 :: k,n
  logical :: dir_e
  call h5open_f(error)
  filename='Results/globals'
  dim1d(1)=1
  call h5fopen_f(trim(filename),H5F_ACC_RDWR_F,file_id,error) 
  call h5dopen_f(file_id,"tempus",dataset_id,error)
  call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,global_variables%tempus,dim1d,error)
  call h5dclose_f(dataset_id,error)
  call h5fclose_f(file_id,error)
  call h5close_f(error)
end subroutine load_simulation_global_parameters
