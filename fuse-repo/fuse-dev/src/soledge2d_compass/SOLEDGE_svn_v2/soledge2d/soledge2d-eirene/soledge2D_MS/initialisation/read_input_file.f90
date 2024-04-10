subroutine read_input_file
  use all_variables
  implicit none
  integer :: finput
  finput=100
  open(unit=finput,file='input.txt',status='unknown')
  call read_simulation_parameters_block(finput)
  call read_chemistry_block(finput)
  call read_boundary_conditions_block(finput)
  call read_transport_parameters_block(finput)
  call read_reference_parameters_block(finput)
  close(finput)
end subroutine read_input_file


subroutine read_simulation_parameters_block(finput)
  use all_variables, only : global_parameters, flags
  use Mreaders
  implicit none
  integer,intent(in) :: finput
  call skip_line(finput,10)
  read(finput,1) flags%is_SLAB
  call skip_line(finput,2)
  read(finput,2) flags%is_to_the_centre
  call skip_line(finput,2)
  read(finput,3) flags%is_Pen
  call skip_line(finput,2)
  read(finput,1) flags%solve_temperature
  call skip_line(finput,2)
  read(finput,5) flags%neutral_model
  call skip_line(finput,2)
  read(finput,10) flags%turbulence_model
  call skip_line(finput,2)
  read(finput,4) flags%restart
  call skip_line(finput,2)
  read(finput,6) global_parameters%N_iterations
  call skip_line(finput,2)
  read(finput,11) global_parameters%N_save
  call skip_line(finput,2)
  read(finput,7) global_parameters%CFL
  call skip_line(finput,2)  
1 format(9X,L1)
2 format(11X,L1)
3 format(8X,L1)
4 format(10X,L1)
5 format(16X,I1)
6 format(7X,I10)
7 format(6X,F5.2)
8 format(8X,I10)
9 format(14X,L1)
10 format(19X,I10)
11 format(8X,I10)
end subroutine read_simulation_parameters_block


subroutine read_chemistry_block(finput)
  use all_variables, only : global_parameters
  use Mreaders
  use Mlog_message
  implicit none
  integer,intent(in) :: finput
  integer*4 :: N_species
  character(128) :: String
  call skip_line(finput,5)
  read(finput,1) global_parameters%N_species
  N_species=global_parameters%N_species
  call allocate_species_properties(N_species)
  call skip_line(finput,2)
  read(finput,2) String
  call parse_line_elements(String,N_species)
  call check_element_list()
  call compute_ion_number()
  call affect_element_and_charge_to_species()
  call write_log_message(msg_element_list)
  call write_log_message(msg_species_list)
  call skip_line(finput,2)
1 format(11X,I2)
2 format(A128)
end subroutine read_chemistry_block


subroutine read_boundary_conditions_block(finput)
  use all_variables, only : global_parameters, boundary_conditions
  use Mreaders
  implicit none
  integer,intent(in) :: finput
  character(128) :: String
  integer*4 :: N_species
  integer*4,allocatable :: bufferI(:)
  real*8,allocatable :: bufferF(:)
  integer*4 :: n
  N_species=global_parameters%N_species
  allocate(bufferI(1:N_species))
  allocate(bufferF(1:N_species))
  call skip_line(finput,6)
  read(finput,1) String
  call parse_line_integer(String,N_species,bufferI)
  do n=1,N_species
     boundary_conditions%BCn_model(n)=bufferI(n)
  end do
  call skip_line(finput,4)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     boundary_conditions%BCn(n)=bufferF(n)
  end do
  call skip_line(finput,3)
  read(finput,1) String
  call parse_line_integer(String,N_species,bufferI)
  do n=1,N_species
     boundary_conditions%BCT_model(n)=bufferI(n)
  end do
  call skip_line(finput,4)
  read(finput,2) boundary_conditions%BCTe
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     boundary_conditions%BCTi(n)=bufferF(n)
  end do
  call skip_line(finput,2)
1 format(A128)  
2 format(8X,ES9.2E2)
  deallocate(bufferI,bufferF)
end subroutine read_boundary_conditions_block


subroutine read_transport_parameters_block(finput)
  use all_variables, only : global_parameters, transport_parameters,&
       ballooning_parameters, flags
  use Mreaders
  implicit none
  integer,intent(in) :: finput
  character(128) :: String
  integer*4 :: N_species
  real*8,allocatable :: bufferF(:)
  integer*4 :: n
  N_species=global_parameters%N_species
  allocate(bufferF(1:N_species))
  call skip_line(finput,5)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%Dn_p(n)=bufferF(n)
  end do
  call skip_line(finput,2)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%Dn_t(n)=bufferF(n)
  end do
  call skip_line(finput,2)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%nu_p(n)=bufferF(n)
  end do
  call skip_line(finput,2)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%nu_t(n)=bufferF(n)
  end do
  call skip_line(finput,2)
  read(finput,2) transport_parameters%chie_p
  call skip_line(finput,2)
  read(finput,3) transport_parameters%chie_t
  call skip_line(finput,2)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%chii_p(n)=bufferF(n)
  end do
  call skip_line(finput,2)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%chii_t(n)=bufferF(n)
  end do
  call skip_line(finput,2)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%v_pinch(n)=bufferF(n)
  end do
  call skip_line(finput,2)
  read(finput,4) ballooning_parameters%ballooning_model
  if(ballooning_parameters%ballooning_model.eq.3) then
     flags%radialFeedback=.true.
  else
     flags%radialFeedback=.false.
  end if
  call skip_line(finput,2)
  read(finput,5) ballooning_parameters%zbal
  call skip_line(finput,1)
  read(finput,6) ballooning_parameters%minmaxbal
  call skip_line(finput,1)
  read(finput,7) ballooning_parameters%sigmabal
  call skip_line(finput,8)
  read(finput,8) transport_parameters%Coulomb_log
  call skip_line(finput,2)
  read(finput,11) transport_parameters%delta_e
  call skip_line(finput,2)
  read(finput,11) transport_parameters%gamma_i
  call skip_line(finput,3)
  read(finput,9) transport_parameters%Flux_limiter(0)
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%Flux_limiter(n)=bufferF(n)
  end do
  read(finput,1) String
  call parse_line_float(String,N_species,bufferF)
  do n=1,N_species
     transport_parameters%Flux_limiter_nu(n)=bufferF(n)
  end do
  call skip_line(finput,2)
 read(finput,10) flags%non_zero_parallel_viscosity
  call skip_line(finput,2)
1 format(A128)  
2 format(13X,F5.2)
3 format(15X,F5.2)
4 format(13X,I1)
5 format(8X,F5.2)
6 format(9X,F5.2)
7 format(11X,F5.2)
8 format(12X,F6.1)
9 format(6X,F5.2)
10 format(10X,L1)
11 format(10X,F5.2)
  deallocate(bufferF)
end subroutine read_transport_parameters_block


subroutine read_reference_parameters_block(finput)
  use all_variables, only : reference_parameters
  use Mreaders
  implicit none
  integer,intent(in) :: finput
  call skip_line(finput,5)
  read(finput,2) reference_parameters%fields%n0
  call skip_line(finput,2)
  read(finput,1) reference_parameters%fields%T0eV
2 format(5X,ES9.2E2)
1 format(5X,F5.2)
end subroutine read_reference_parameters_block
