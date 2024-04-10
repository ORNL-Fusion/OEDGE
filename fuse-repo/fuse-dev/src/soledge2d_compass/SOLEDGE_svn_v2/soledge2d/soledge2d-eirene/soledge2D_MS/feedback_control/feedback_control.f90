module Mfeedback_control

  type :: TControl_data
     real*8 T_target
     real*8 density_target
     integer*4 i_target
     integer*4 j_target
     integer*4 k_target
     real*8 Gain
     real*8 tau_d
     real*8 tau_i
     real*8 max_puff
     real*8 min_puff
 end type TControl_data

  type :: Tcontrol_error_data
     real*8 integral_error
     real*8 error_old
  end type Tcontrol_error_data

  type(TControl_data) :: Control_data
  type(TControl_error_data) :: error_data

contains

  subroutine init_control(Control_dat,error_data)
    use all_variables, only : global_parameters, zones, reference_parameters, flags
    implicit none
    type(TControl_data),intent(inout) :: Control_dat
    type(Tcontrol_error_data),intent(inout) :: error_data
    integer*4 i,j,k
    open(unit=10,file='feedback_data.txt',status='unknown')
    read(10,*) Control_dat%density_target
    read(10,*) Control_dat%i_target
    i=Control_dat%i_target
    read(10,*) Control_dat%j_target
    j=Control_dat%j_target
    read(10,*) Control_dat%k_target
    k=Control_dat%k_target
    read(10,*) Control_dat%Gain
    read(10,*) Control_dat%tau_i
    read(10,*) Control_dat%tau_d
    read(10,*) Control_dat%max_puff
    read(10,*) Control_dat%min_puff
    close(10)
    if(flags%restart) then
       open(unit=10,file='feedback_puff_save',status='unknown')
       read(10,*) error_data%integral_error
       close(10)
    else
       error_data%integral_error=0.D0
    end if
    error_data%error_old=Control_dat%density_target-Zones(k)%species(1)%var(1)%density(i,j)&
         *reference_parameters%fields%n0
  end subroutine init_control

  subroutine puff_feedback(Control_dat,Puff,error_data)
    use all_variables, only : global_parameters, zones, reference_parameters, global_variables
    implicit none
    type(TControl_data),intent(inout) :: Control_dat
    real*8,intent(inout) :: Puff
    type(Tcontrol_error_data),intent(inout) :: error_data
    real*8 error
    integer*4 i,j,k
    i=Control_dat%i_target
    j=Control_dat%j_target
    k=Control_dat%k_target
    error=Control_dat%density_target-Zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0
    error_data%integral_error=error_data%integral_error+global_variables%dt*reference_parameters%fields%tau0&
         *error
    Puff=Control_dat%Gain*(error)+&
         1.D0/Control_dat%tau_i*error_data%integral_error!+&
    !         Control_dat%tau_d*(error-error_data%error_old)/(dt*dat%tau0))
    if(Puff.le.Control_dat%min_puff) then
       Puff=Control_dat%min_puff
    end if
    if(Puff.ge.Control_dat%max_puff) then
       Puff=Control_dat%max_puff
    end if
    if(isnan(Puff)) then
       Puff=1.D17
    end if
  end subroutine Puff_feedback

  subroutine write_feedback(tempus,Puff,Control_dat,ntps)
    use all_variables, only : global_parameters, zones, reference_parameters
    Implicit none
    type(TControl_data),intent(in) :: Control_dat
    real*8,intent(in) :: Puff
    real*8,intent(in) :: tempus
    integer*4,intent(in) :: ntps
    integer*4 i,j,k
    if((modulo(dble(ntps),dble(global_parameters%N_iterations)/1000.D0).eq.0).or.(global_parameters%N_iterations.lt.1000)) then
       i=Control_dat%i_target
       j=Control_dat%j_target
       k=Control_dat%k_target
       open(unit=10,file='puff_data',status='unknown',position='append')
       write(10,100) tempus, Puff, Control_dat%density_target,&
            Zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0
100    format(512es15.7)
    end if
  end subroutine write_feedback

  subroutine save_control(Control_dat,error_data)
    use all_variables, only : global_parameters, zones, reference_parameters, flags
    implicit none
    type(TControl_data),intent(inout) :: Control_dat
    type(Tcontrol_error_data),intent(inout) :: error_data
    open(unit=10,file='feedback_puff_save',status='unknown')
    write(10,110) error_data%integral_error
    close(10)
110 format(es15.7)
  end subroutine save_control

end module Mfeedback_control
