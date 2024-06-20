subroutine read_kepsilon_param()
  use all_variables, only : kepsilon_param
  implicit none
  open(unit=10,file='kepsilon_parameters',status='unknown')
  read(10,1) kepsilon_param%Cmu
  read(10,1) kepsilon_param%C1e
  read(10,1) kepsilon_param%C2e
  read(10,1) kepsilon_param%C3e
  read(10,2) kepsilon_param%sigma_k
  read(10,2) kepsilon_param%sigma_epsilon
  read(10,2) kepsilon_param%sigma_n
  read(10,2) kepsilon_param%sigma_v
  read(10,2) kepsilon_param%sigma_T
  read(10,3) kepsilon_param%C_interchange
  read(10,4) kepsilon_param%th_interchange
  read(10,7) kepsilon_param%gradT_weight
  read(10,5) kepsilon_param%C_kh
  read(10,6) kepsilon_param%mu_max
  read(10,6) kepsilon_param%mu_min
  read(10,9) kepsilon_param%deltaOmega
  read(10,8) kepsilon_param%tauV
1 format(6X,F5.2)
2 format(10X,F5.2)
3 format(16X,es9.2e2)
4 format(17X,F8.2)
5 format(7X,es9.2e2)
6 format(9X,F8.2)
7 format(15X,F8.2)
8 format(7X,es9.2e2)
9 format(13X,es9.2e2)
  close(10)
end subroutine read_kepsilon_param
