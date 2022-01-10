program run_nc_test
use nc_utils_generic

!
! Note: version 1 is not used so testing has been commented out and it is not included in the compilation
!
!use nc_utils_generic_v1

integer :: res


!write(0,*) 'TESTING VERSION 1:'
!res=test_nc_utils_v1()
!write(0,*) 'Test nc_utils_v1:',res


write(0,*) 'TESTING VERSION 2:'
res=test_nc_utils()
write(0,*) 'Test nc_utils_v2:',res


end program run_nc_test
