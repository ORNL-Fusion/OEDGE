module hc_kinetics

  !use mod_params - some routines use hc_init_lib_data which uses comhc which in turn uses the global parameters 
  use hc_velocity_type
  use hc_kinetics_options
  use error_handling
  !
  ! This module implements the changes in velocity vector of the HC fragment as it evolves
  ! Each reaction releases a certain amount of energy - a mass dependent fraction of that energy
  ! is imparted to the HC fragment in the center of mass frame.
  ! Treatment of CX reactions is somewhat different from the electron reactions. Since the energy 
  ! release is not isotropic in the center of mass frame. 
  !
  save

  integer,private,parameter :: veclen=3


  integer,private :: transition_count
  integer,private :: last_iprod = 0



contains



  subroutine update_hc_kinetics(hc_v,reaction_type,current_cell,current_ring,  & 
       & cur_hc_spec,reaction_index,followed_state,h_isotope_composition,  &
       & current_r,current_z,current_s, current_theta,current_cross,current_angle,  &
       & current_velocity_in_s,current_velocity_in_r,current_velocity_in_z, &
       & current_velocity, kin_energy_added,s_start,&
       & neutral_time_step,ion_time_step,sput_weight,iprod)

    use hc_init_lib_data
    use hc_get
    use bfield
    use velocity_dist
    use hc_diag_data

    implicit none

    type (hc_velocity_type1) :: hc_v
    integer,intent(in) :: reaction_type
    ! current_cell can be changed by getscross_approx for particles in the first/last cell of core rings
    integer :: current_cell,current_ring
    integer,intent(in) :: cur_hc_spec,reaction_index,followed_state
    Integer, Dimension (Number_H_Species), Intent (In) :: H_Isotope_Composition
    !Integer, Dimension (Number_H_Products), Intent (In) :: H_Isotope_Composition
    real :: current_r,current_z,current_s,current_theta,current_cross,current_angle
    real :: current_velocity_in_s,current_velocity_in_r,current_velocity_in_z
    real :: current_velocity
    !real :: current_energy
    real :: kin_energy_added
    real :: s_start
    real,intent(in) :: neutral_time_step,ion_time_step
    real,intent(in) :: sput_weight
    integer :: iprod

    !
    ! Local variables
    !
    real :: cur_hc_mass,total_mass
    real :: e_fraction,e_release,e_hc_gained
    integer :: e_type
    integer :: reaction_kind
    real :: b(3)  ! local magnetic field vector

    real,external :: atan2c
    real :: beta,psi,anglan,tanlan

    character :: hc_species_name*10

    ! Krieger IPP/07 - SUN compiler insists on 132 column limit
    !if (debug_kinetics) &
    !   & write(6,'(a,a4,4i4,25(1x,g9.2))') 'UPDATE_START:',hc_ident_species(cur_hc_spec),current_cell,current_ring,reaction_index,&
    !   &             reaction_type, hc_v%v(1),hc_v%v(2),hc_v%v(3),hc_v%vpara,hc_v%vperp,hc_v%vtot,mag_v(hc_v%v),current_velocity, &
    !   &             current_velocity_in_s,current_velocity_in_r/neutral_time_step,current_velocity_in_z/neutral_time_step, &
    !   &             current_r,current_z,current_s, current_theta*raddeg,current_cross,current_angle*raddeg,s_start,
    !   &             kin_energy_added,e_fraction,e_release,cur_hc_mass

    ! Krieger IPP/07 - SUN compiler insists on 132 column limit
    if (debug_kinetics) &
       & write(6,'(a,a4,4i4,a30,25g12.3)') 'UPDATE_START:',hc_ident_species(cur_hc_spec),current_cell,current_ring,reaction_index,&
       &             reaction_type, &
       &             hc_state_transform_table(reaction_index)%reaction_desc&
       &             (1:len_trim(hc_state_transform_table(reaction_index)%reaction_desc)),&
       &             hc_v%v(1),hc_v%v(2),hc_v%v(3),hc_v%vpara,hc_v%vperp,hc_v%vtot,current_velocity, &
       &             current_velocity_in_s,current_velocity_in_r/neutral_time_step,current_velocity_in_z/neutral_time_step, &
       &             current_r,current_z,current_s, current_theta,current_cross,current_angle*raddeg

    !
    ! Initialization
    !
    !
    ! Initialize transition count if this is a new particle - otherwise increment it
    !
    if (iprod.ne.last_iprod) then 
       transition_count = 1
    else
       transition_count = transition_count + 1
    endif

    !
    ! Define the kind of the reaction based on whether it is an electron or proton interaction
    !
    if (HC_State_Transform_Table(reaction_index) % Reaction_Type .eq.'e') then 
       reaction_kind = e_reaction
    elseif (HC_State_Transform_Table(reaction_index) % Reaction_Type .eq.'p') then 
       reaction_kind = p_reaction
    endif

    !
    ! assign name of current species
    !

    hc_species_name = hc_ident_species(cur_hc_spec)

    !
    ! UPDATE PARTICLE COORDINATES IF NECESSARY
    !
    ! Update the positional variables for ion to neutral and neutral to ion reactions. For ion to ion or neutral 
    ! to neutral - no change is required in the coordinates of the particles since they are still in the same
    ! coordinate system. (Ions use S,CROSS while neutrals use R,Z)
    !
    ! Neutral to ion
    !
    if (reaction_type.eq.ni_reaction) then

       ! SET Initial S and CROSS postion for particles.
       ! Note:  Check for control option is removed given that any HC can
       ! be launched with control option 0 so that should always be used.

       Call getscross_approx (Current_R, Current_Z, Current_S, Current_Cross, Current_Cell, Current_Ring)

       ! Record starting S-distance from nearest target.

       S_Start = MIN (Current_S, gksmaxs (Current_Ring) - Current_S)

       !
       !   Calculate theta is in ion_parallel_transport.f
       !
       call calculate_theta(current_cell,current_ring,current_s,current_theta)

       current_angle = 0.0

       !
       ! Ion to neutral
       !
    elseif (reaction_type.eq.in_reaction) then 

       !
       ! Hard code option 2 to getrz for now
       !
       call getrz(current_cell,current_ring,current_s,current_cross,current_r,current_z,2)

    endif

    !
    ! Update the hc velocity structure to reflect the current values in use for the particle motion
    ! Notes: For ions - the parallel velocity will have been evolving over time and will not be reflected in the vpara value
    !                 - this is relevant for ii and in reactions
    !        For neutrals - the only way the R and Z velocities can change is through a momentum transfer collision - code
    !                       is also implemented to update the hc_v structure for these interactions - at least for the advanced
    !                       model - however, to ensure the proper values these should be copied over to the 
    !                       R,Z velocities in case any additional code is added to change these values that does not get 
    !                       reflected in hc_v
    !

    hc_v%vpara = current_velocity_in_S  ! current_velocity_in_s does NOT implicitly contain the time step in the HC code
    hc_v%v(1) = current_velocity_in_R/ neutral_time_step  ! vr
    hc_v%v(2) = current_velocity_in_Z/ neutral_time_step  ! vz


    !
    ! Get current hc mass
    !

    cur_hc_mass = Find_HC_Mass (cur_hc_spec)

    ! Determine the energy released to the HC fragment. 
    ! - first check the type information regarding the hc_e value
    ! - hc_e_type = 0 -> no HC energy gain
    !             = 1 -> Energy listed is amount given to HC fragment
    !             = 2 -> Energy given is a fraction of the background plasma temperature (not sure this makes sense)
    !             = 3 -> Energy given is divided among all fragments inversely proportional to mass
    !             = 4 -> Energy given is divided only among charged fragments


    e_type     = hc_state_transform_table(reaction_index) % hc_e_type
    e_release  = HC_State_Transform_Table (Reaction_index) % HC_E

    !
    ! note: e_release and thus e_hc_gained are expected to be in units of eV
    !

    if (e_type.eq.0) then
       e_hc_gained = 0.0
    elseif (e_type.eq.1) then 
       e_hc_gained = e_release
    elseif (e_type.eq.2) then 
       ! this is copied from Adam's code for compatibility - doesn't make a lot of sense 
       if (hc_state_transform_table(reaction_index)%reaction_type.eq.'p') then
          e_hc_gained = e_release * gktebs (Current_Cell,Current_Ring)
       Else
          e_hc_gained = e_release * gktibs (Current_Cell,Current_Ring)
       End If
    elseif (e_type.eq.3.or.e_type.eq.4) then 
       !
       ! divide energy inversely proportional to mass
       !
       ! Sum up product particle masses - either all or just charged
       !
       ! Note: it will require substantial code modifications to record exactly which 
       !       H_isotopes were lost as a by-product of the state change in order to precisely
       !       calculate the mass of the CH fragment and the other reactants. On the assumption
       !       that mixed isotope molecules are unlikely to be used yet - the code here assumes
       !       the base molecular make up and background hydrogen istotope composition for the molecule. 
       !    
       ! It uses the precalculated values for state_mass and state_charge as well as the 
       ! reaction_product_mass and reaction_product_charge
       !
       if (e_type.eq.3) then 
          ! Energy divided among all reaction products
          e_fraction = (hc_state_transform_table(reaction_index)%reaction_product_mass - cur_hc_mass)/ &
               & hc_state_transform_table(reaction_index)%reaction_product_mass

       elseif (e_type.eq.4) then 
          !
          ! Energy divided among charged reaction products
          !

          if (hc_state_table(cur_hc_spec)%state_charge.eq.0) then 

             e_fraction = 0.0

          else

             e_fraction = (hc_state_transform_table(reaction_index)%charged_reaction_product_mass - cur_hc_mass)/ &
                  & hc_state_transform_table(reaction_index)%charged_reaction_product_mass
             ! jdemod - if the HC fragment is the only charged fragment then it should get all the energy - if there
             !          are more than one charged fragment then the energy is distributed inversely proportional
             !          to mass. 
             if (e_fraction.lt.0.0) e_fraction = 1.0


          endif

       endif

       e_hc_gained = e_fraction * e_release

    endif

    !
    ! Record energy gained
    !
    kin_energy_added = kin_energy_added + e_hc_gained

    !
    ! Obtain the local magnetic field vector (3D)
    !
    call get_bvector(current_cell,current_ring,b)

    ! UPDATE PARTICLE VELOCITY

    call calculate_hc_velocity(hc_v,cur_hc_mass,e_hc_gained,reaction_type,reaction_kind,b)


    !
    ! UPDATE CURRENT_ANGLE, CURRENT_VELOCITY VALUES
    !


    if (reaction_type.eq.ni_reaction.or.reaction_type.eq.ii_reaction) then 
       !
       ! Update ion related velocities
       !
       ! current_velocity_in_s does NOT contain the time step in the HC code - except in the transport routine ion_move
       !
       current_velocity_in_s = hc_v%vpara

       !
       ! For ions - set the current velocity to the total velocity when the ion was created
       !            This is not used in the code at present since the code uses current_velocity only for neutrals - in which
       !            case it is the velocity in the poloidal or R,Z plane. 
       !
       current_velocity = hc_v%vtot

    elseif (reaction_type.eq.in_reaction.or.reaction_type.eq.nn_reaction) then 
       !
       ! Update neutral related velocities
       !
       current_velocity_in_r = hc_v%v(1) * neutral_time_step
       current_velocity_in_z = hc_v%v(2) * neutral_time_step

       current_angle = atan2c(hc_v%v(2),hc_v%v(1))

       current_velocity = sqrt(hc_v%v(1)**2 + hc_v%v(2)**2)

    endif

    !
    ! Just assign kinetic energy for now
    !

    !current_energy = 0.5 * cur_hc_mass * amu * hc_v%vtot**2


    !
    ! Invoke velocity distribution debugging if the state is 'C' - carbon neutral
    !
    if (hc_species_name(1:len_trim(hc_species_name)).eq.'C') then 

       !
       ! Note - tanlan is meaningless - beta and psi have to be mapped from the velocity 
       !        components. 
       !
       ! Also - in neut - vin is the velocity in the R,Z plane - in the current case I will be
       !        passing in the total 3D velocity for recording instead.
       !
       ! The definition of beta and psi is not simple - anglan is the angle from the surface normal OR 
       ! from the positive Z axis - counterclockwise. 
       ! BETA is an angle measured from the normal (or positive R axis) while PSI is an angle around 
       ! the normal or R axis (in the absense of a reference tanlan defining the normal). 
       ! PSI could be measured relative to either the Z or T axes and could be defined in either direction.
       ! To be consistent - for now we will define PSI as counter clockwise from the positive Z axis.
       ! Given these definitions - the components of the velocity can be used to calculate a beta and psi
       ! value as well as an anglan value which lies in the R,Z plane relative to the positive R axis. 
       !
       beta = acos(hc_v%v(1)/hc_v%vtot)

       psi = atan2c(hc_v%v(3),hc_v%v(2))
       
       anglan = atan2c(sin(beta)*cos(psi),cos(beta))

       ! Krieger IPP/07 - SUN compiler insists on 132 column limit
       if (debug_kinetics) &
            & write(6,'(a,10(1x,g12.5))') 'DEBUGV:',hc_v%vtot,hc_v%v(1),hc_v%v(2),hc_v%v(3),beta*raddeg,psi*raddeg, &
            & anglan*raddeg,current_angle*raddeg

       tanlan = 0.0

       call record_vdist(hc_v%vtot,beta,psi,tanlan,anglan,sput_weight)

       call record_c_detailed_vdist(hc_v,sput_weight)


    endif


    call record_state_energy_diag_data(cur_hc_spec,hc_v,sput_weight,e_hc_gained,kin_energy_added,cur_hc_mass,transition_count)


    !
    ! Record the particle identifier
    !
    last_iprod = iprod


    ! Krieger IPP/07 - SUN compiler insists on 132 column limit
    if (debug_kinetics) &
       & write(6,'(a,a4,2i4,25(1x,g14.5))') 'UPDATE_E   :',hc_ident_species(cur_hc_spec),reaction_index,e_type,e_fraction, &
       &                                    e_release,e_hc_gained,cur_hc_mass

    ! Krieger IPP/07 - SUN compiler insists on 132 column limit
    !if (debug_kinetics) &
    !   & write(6,'(a,a4,4i4,25(1x,g9.2))') 'UPDATE_END  :',hc_ident_species(cur_hc_spec),current_cell,current_ring,reaction_index,&
    !   &             reaction_type, hc_v%v(1),hc_v%v(2),hc_v%v(3),hc_v%vpara,hc_v%vperp,hc_v%vtot,mag_v(hc_v%v),current_velocity, &
    !   &             current_velocity_in_s,current_velocity_in_r,current_velocity_in_z, &
    !   &             current_r,current_z,current_s, current_theta*raddeg,current_cross,current_angle*raddeg,s_start,
    !   &             kin_energy_added

    ! Krieger IPP/07 - SUN compiler insists on 132 column limit
    if (debug_kinetics) &
       & write(6,'(a,a4,4i4,a30,25g12.3)') 'UPDATE_END  :',hc_ident_species(cur_hc_spec),current_cell,current_ring,reaction_index,&
       &             reaction_type,   &
       &             hc_state_transform_table(reaction_index)%reaction_desc(1:len_trim(hc_state_transform_table&
       &             (reaction_index)%reaction_desc)),&
       &             hc_v%v(1),hc_v%v(2),hc_v%v(3),hc_v%vpara,hc_v%vperp,hc_v%vtot,current_velocity, &
       &             current_velocity_in_s,current_velocity_in_r/neutral_time_step,current_velocity_in_z/neutral_time_step, &
       &             current_r,current_z,current_s, current_theta,current_cross,current_angle*raddeg


  end subroutine update_hc_kinetics






  subroutine assign_hc_velocity(hc_v,                  &
       &                        launch_velocity,       &
       &                        velocity_multiplier,   &
       &                        current_angle,         & 
       &                        psi,                   &
       &                        hc_temperature,        &
       &                        current_velocity,      &
       &                        current_velocity_in_s, &
       &                        current_velocity_in_r, &
       &                        current_velocity_in_z, &
       &                        neutral_time_step,     &
       &                        ion_time_step,         &
       &                        charge,                &
       &                        cur_hc_mass)
    use mod_params
    implicit none
    !
    ! This routine sets the initial values of the velocity based on the velocity/angle flags selected
    ! for the particle launch and the type of launch - the routine is passed the angles and magnitude of the 
    ! assigned velocity - which are then used to come up with a 3-space representation for the initial 
    ! particle velocity. In the case where the launch is completely 2D - the T component of the velocity is
    ! set to zero. 
    !
    type (hc_velocity_type1) :: hc_v
    real,intent(in)    :: cur_hc_mass
    integer,intent(in) :: charge
    real,intent(in) :: launch_velocity,velocity_multiplier,current_angle,psi,neutral_time_step,hc_temperature,ion_time_step
    real,intent(out) :: current_velocity,current_velocity_in_s,current_velocity_in_r,current_velocity_in_z

    !
    ! Local variables 
    !
    real vr,vz,vt
    real temp

    !
    ! Initialize transition counter
    !
    transition_count = 0

    !
    ! Assign initial temperature values 
    !

    hc_v%t     = hc_temperature
    hc_v%tperp = hc_temperature
    hc_v%tpara = hc_temperature


    ! Assign the basic velocity quantities used in the original code


    ! Remember original velocity.  Current_Velocity will be updated as the particle moves.
    If (charge .eq. 0) Then

       !
       ! Assign initial temperature values 
       !

       hc_v%t     = hc_temperature
       hc_v%tperp = 0.0
       hc_v%tpara = 0.0


       Current_Velocity = Launch_Velocity

       vr =  Current_Velocity * COS (Current_Angle)
       vz =  Current_Velocity * SIN (Current_Angle)

       ! Calculate R,Z components of velocity, probability of ionization, etc.
       Current_Velocity_In_R = vr *  Neutral_Time_Step
       Current_Velocity_In_Z = vz *  Neutral_Time_Step

       !
       ! Ion velocity is not significant for neutrals - set to zero to initialize it
       !
       current_velocity_in_s = 0.0

       !
       ! Calculate the components of hc_v - the hc velocity structure
       !

       hc_v%v(1) = vr   ! vr
       hc_v%v(2) = vz   ! vz 

       temp = (current_velocity/velocity_multiplier)**2 - vr**2 -vz**2
       !
       ! velocity multiplier is related to the component of velocity in the T direction - if velocity multiplier is 1.0 then
       ! all of the particle velocity should be in the R,Z plane.
       !
       if (temp.lt.0.0.or.velocity_multiplier.eq.1.0) then 
          hc_v%v(3) = 0.0
       else
          hc_v%v(3) = sqrt(temp) * sign(1.0,psi-PI)  ! vt
       endif

       ! Krieger IPP/07 - SUN compiler insists on 132 column limit
       if (debug_kinetics) &
            & write(6,'(a,25(1x,g12.5))') 'ASSIGN V   :',vr,vz,hc_v%v(3),current_velocity,current_angle*raddeg, &
            & velocity_multiplier,psi,PI,sign(1.0,psi-PI), &
            & (current_velocity/velocity_multiplier)**2 - vr**2 -vz**2


       hc_v%vperp  = 0.0
       hc_v%vpara  = 0.0

       hc_v%vtot = sqrt(sum(hc_v%v**2))
       
    Else

       !
       ! Assign initial temperature values 
       !

       hc_v%t     = hc_temperature
       hc_v%tperp = hc_temperature
       hc_v%tpara = hc_temperature



       ! UPDATE: The fact that the launch_velocity for ions was scaled by the time_step was a bug introduced by carrying code over from DIVIMP 
       !         then changing the design. launch_velocity is now in m/s.
       !
       ! For ions - launch_velocity comes from HC_injection_velocity and the ion time step has already been applied
       ! NOTE! The ion time step should NOT be applied since current_velocity_in_s in the HC code is in m/s not dist/time-step - except
       !       in the actual time step code. 


       Current_Velocity_In_S = Launch_Velocity

       hc_v%vpara = current_velocity_in_s

       !
       ! Neutral velocities are not significant for ions - set to zero to initialize them
       !
       current_velocity_in_r = 0.0
       current_velocity_in_z = 0.0

       !
       ! Calculate the components of hc_v - the hc velocity structure
       !

       hc_v%v(1) = 0.0   ! vr
       hc_v%v(2) = 0.0   ! vz 
       hc_v%v(3) = 0.0   ! vt

       if (hc_vperp_opt.eq.0) then 
          hc_v%vperp  = 0.0
       elseif (hc_vperp_opt.eq.1.or.hc_vperp_opt.eq.3) then 
          ! hc_verp_opt 3 diffuses vperp in the same way as vpara - on particle creation set them to the same value
          hc_v%vperp = hc_v%vpara
       elseif (hc_vperp_opt.eq.2) then 
          call convert_temp_velocity(hc_v%vperp,hc_v%tperp,cur_hc_mass,conv_to_v)
       endif

       hc_v%vtot = sqrt(hc_v%vperp**2 + hc_v%vpara**2)

    End If

    if (debug_kinetics) &
! slmod begin - bug
! The '25g' causes problems for the Sun compiler, just removed 
! the 'g' since the format specifier looks invalid.  -SL 27.03.07
       & write(6,'(a,25(1x,g9.2))') 'ASSIGN      :',  &
!
!       & write(6,'(a,25g(1x,g9.2))') 'ASSIGN      :',  &
! slmod end
       &             hc_v%v(1),hc_v%v(2),hc_v%v(3),hc_v%vpara,hc_v%vperp,hc_v%vtot,mag_v(hc_v%v),current_velocity, &
       &             current_velocity_in_s,current_velocity_in_r,current_velocity_in_z, &
       &             hc_v%t,hc_v%tpara,hc_v%tperp,cur_hc_mass


  end subroutine assign_hc_velocity



  subroutine calculate_hc_velocity(v,m,energy,reaction_type,reaction_kind,b)
    use mod_params
    implicit none

    !
    ! energy is the energy released in the interaction
    ! v is the velocity of the HC fragment 
    ! m is the mass of the HC fragment 
    !
    !
    ! reaction_type - 1 - neutral -> neutral
    !                 2 - neutral -> ion
    !                 3 - ion -> neutral
    !                 4 - ion -> ion
    !
    ! reaction_kind - 1 - electron interaction
    !               - 2 - proton/HDT interaction
    !
    !
    ! NOTE: All these routine assume that the simulation time step is not 
    !       factored into the velocities for the calculations.Since this code
    !       is being layered on top of existing code - it will be necessary to 
    !       apply the appropriate scalings when the assignment to existing variables
    !       takes place. 
    !


    type (hc_velocity_type1) :: v
    real :: b(3)     ! magnetic field components at current location
    real :: energy   ! units are in eV
    real :: m        ! units are in amu
    integer :: reaction_type ! 1 to 4
    integer :: reaction_kind ! 1 = electron  2=H/D/T (eg CX)


    !
    ! local variables
    !
    real deltav

    !
    ! for ion source reactions - need to convert ion velocity to 3D velocities
    !

    if (debug_kinetics) &
         & write(6,'(a,10(1x,g12.5))') 'V1:', v%v(1),v%v(2),v%v(3),v%vpara,v%vperp,mag_v(v%v),sqrt(v%vpara**2+v%vperp**2)


    if (reaction_type.eq.in_reaction.or.reaction_type.eq.ii_reaction) then 

       !
       ! Need to adjust the value of vperp depending on the option selected
       ! option 0 - vperp does not change between reactions
       ! option 1 - equi-partition - vperp = vpara
       ! option 2 - ensemble heating/cooling - vperp is calculated from tperp 
       ! option 3 - vperp changes diffusively with the same diffusive step sizes as vpara but randomly assigned independently
       !            vperp is unaffected by other parallel forces - vperp does not need to be changed here
       !
       if (hc_vperp_opt.eq.1) then 
          v%vperp = v%vpara
       elseif (hc_vperp_opt.eq.2) then 
          call convert_temp_velocity(v%vperp,v%tperp,m,conv_to_v)
       endif

       call map_hcv(v,b,map_to_3D)  

    endif

    if (debug_kinetics) &
         & write(6,'(a,10(1x,g12.5))') 'V2:', v%v(1),v%v(2),v%v(3),v%vpara,v%vperp,mag_v(v%v),sqrt(v%vpara**2+v%vperp**2)

    !
    ! Update the 3-space velocity by the energy released in the reaction
    !

    if (energy.ge.0.0) then 
       deltav = sqrt(2.0 * emi * energy / m ) 
    else
       call errmsg('Calculate_HC_Velocity','Error: Negative energy value passed in!')
    endif


    call apply_deltav(v,deltav,reaction_kind)

    !
    ! Calculate total magnitude of the vector before mapping back to bfield (if req'd)
    !

    v%vtot = mag_v(v%v)


    if (debug_kinetics) &
         & write(6,'(a,10(1x,g12.5))') 'V3:', v%v(1),v%v(2),v%v(3),v%vpara,v%vperp,mag_v(v%v),sqrt(v%vpara**2+v%vperp**2)


    !
    ! for ion destination reactions - need to convert the 3D velocity into the parallel and perpendicular components
    !

    if (reaction_type.eq.ni_reaction.or.reaction_type.eq.ii_reaction) then 

       call map_hcv(v,b,map_to_bfield)

       !
       ! Need to adjust the value of vperp depending on the option selected
       ! option 0 - vperp does not change between reactions
       ! option 1 - equi-partition - vperp = vpara
       ! option 2 - ensemble heating/cooling - vperp is calculated from tperp 
       !
       if (hc_vperp_opt.eq.2) then 
          call convert_temp_velocity(v%vperp,v%tperp,m,conv_to_t)
       endif

       if ((sqrt(v%vpara**2+v%vperp**2)-v%vtot).gt.0.01) then
          write(6,'(a,4(1x,g12.5))') 'HC IONV ERR:',v%vtot,sqrt(v%vpara**2+v%vperp**2),v%vpara,v%vperp
       endif

    endif

    if (debug_kinetics) &
         & write(6,'(a,10(1x,g12.5))') 'V4:', v%v(1),v%v(2),v%v(3),v%vpara,v%vperp,mag_v(v%v),sqrt(v%vpara**2+v%vperp**2)

    if (debug_kinetics) &
       & write(6,'(a,2i4,25(1x,g14.5))') 'UPDATE_VEL  :',  reaction_type,reaction_kind,&
       &             deltav,energy,m,emi,b(1),b(2),b(3),mag_v(b)



  end subroutine calculate_hc_velocity

  subroutine reset_hc_velocity(hc_v,current_velocity,current_angle,toroidal_angle,hc_temperature)

    use mod_params
    implicit none

    ! Note:!: In this routine current_velocity is the total 3D velocity of the particle - it is subsequently adjusted by
    !         the cos(toroidal_angle) in the calling routine to convert it to the total velocity in the poloidal plane.
    !
    ! This routine updates the contents of hc_v when the values of current_velocity and current_angle
    ! change in the main code - e.g. at wall impacts. This is a stopgap solution until better wall interaction
    ! code is availble.

    type (hc_velocity_type1) hc_v

    real :: current_velocity,current_angle
    real :: toroidal_angle
    real :: hc_temperature

    !
    ! Reset transition count to zero (?)
    !
    transition_count = 0


    if (debug_kinetics) &
         & write(6,'(a,25(g14.5))') 'RESET HC_V1:',hc_v%v(1),hc_v%v(2),hc_v%v(3),current_velocity,&
         &            mag_v(hc_v%v),current_angle,toroidal_angle


    hc_v%v(1) = current_velocity * cos(toroidal_angle) * cos (current_angle)  ! vr
    hc_v%v(2) = current_velocity * cos(toroidal_angle) * sin (current_angle)  ! vz
    hc_v%v(3) = current_velocity * sin(toroidal_angle)                        ! vt
    hc_v%vtot = mag_v(hc_v%v)

    hc_v%vperp = 0.0
    hc_v%vpara = 0.0
    hc_v%t     = hc_temperature
    hc_v%tperp = hc_temperature
    hc_v%tpara = hc_temperature

    if (debug_kinetics) &
         & write(6,'(a,25(g14.5))') 'RESET HC_V2:',hc_v%v(1),hc_v%v(2),hc_v%v(3),current_velocity,&
         &       mag_v(hc_v%v),current_angle*raddeg,toroidal_angle*raddeg


  end subroutine reset_hc_velocity


  subroutine convert_temp_velocity(vperp,tperp,m,conv_opt)
    use mod_params
    implicit none
    real :: vperp,tperp
    real :: m
    integer :: conv_opt


    !
    ! There are several possibilities for converting - 
    ! 1/2 m v**2 =  1/2 kT, kT or 3/2 kT 
    !
    ! For now the code uses kT 
    ! scale_fact = 1.0 - for kT
    !
    ! scale_fact = 1.0/2.0 - for 1/2 kT
    ! scale_fact = 3.0/2.0 - for 3/2 kT
    !
    ! T (eV), m (amu), vperp (m/s)
    !

    real, parameter :: scale_fact = 1.0

    

    if (conv_opt.eq.conv_to_t) then 

       tperp = m * vperp**2 / (scale_fact * 2.0 * emi) 

    elseif (conv_opt.eq.conv_to_v) then 

       vperp = sqrt (scale_fact * 2.0 * emi * tperp / m) 

    endif


    if (debug_kinetics) &
       & write(6,'(a,i4,10(1x,g12.5))') 'CVRT:',conv_opt,vperp,tperp





  end subroutine convert_temp_velocity


  subroutine map_hcv(v,b,map_option)
    implicit none

    type(hc_velocity_type1) :: v
    real :: b(3)
    integer map_option

    !
    ! Local variables
    !
    real :: bperp(3),vtot2

    if (map_option.eq.map_to_3D) then 
    !----------------------------------------------------------
    !
    ! Map the contents of Vperp and Vpara to equivalent vr,vz,vt components 
    ! - at the present time this involves
    !   assuming a random direction for Vperp in the plane perpendicular to b
    !


       !
       ! Vpara component is along the direction of b - b is assumed to be a unit vector
       ! Vperp is in a random direction in the plane perpendicular to b
       !
       ! first - obtain a random bperp unit vector
       !
       ! 
       call get_random_perp_vector(b,bperp) 

       !
       !  Perform vector sum to obtain vr,vz and vt components 
       !
       ! R
       v%v(1) = v%vpara * b(1) + v%vperp * bperp(1)
       ! Z
       v%v(2) = v%vpara * b(2) + v%vperp * bperp(2)
       ! T
       v%v(3) = v%vpara * b(3) + v%vperp * bperp(3)


    elseif (map_option.eq.map_to_bfield) then 

       !----------------------------------------------------------
       !
       ! Map the contents of vr,vz,vt to equivalent values of vperp and vpara - note that 
       ! the direction data for vperp is discarded at the present time. Future implementations
       ! may retain this data. 
       !


       !
       ! Calculate dot product to get vpara
       !
       !   br * vr + bz * vz + bt * vt
       v%vpara = b(1) * v%v(1) + b(2) * v%v(2) + b(3) * v%v(3) 

       vtot2 = v%v(1)**2 + v%v(2)**2 + v%v(3)**2

       v%vperp = sqrt(vtot2 - v%vpara**2)


    endif


  end subroutine map_hcv



  subroutine apply_deltav(v,deltav,reaction_kind)
    implicit none


    !
    ! This routine only modifies the 3D mapped velocities not bfield projected
    ! velocities
    !
    ! deltaV is the velocity gained by the particle
    !
    ! reaction_kind determines the direction of the gained energy
    ! 0 - 3D isotropic along a random vector
    ! 1 - along the direction of v
    !
    ! Discussions at the present time seem to indicate that the electron 
    ! reactions will release energy 3D isotropically. On the other hand, it
    ! isn't clear what happens for the proton/CX reactions. It appears that the 
    ! energy may be absorbed by the internal states of the molecule and to first
    ! approximation there would be neither an energy gain or change in trajectory - 
    ! at least for resonant reactions. It is possible that the intermediate 
    ! molecular state for the non-resonant reactions may exist for a sufficiently 
    ! long time that the particle directions upon molecular breakup are affected. 
    ! However, one would expect these to still be subject to conservation laws. 
    !
    !

    type (hc_velocity_type1) :: v
    real :: deltav
    integer reaction_kind

    !
    ! local variables
    !
    real vadd(3)

    !
    ! At the present time do nothing for p_reactions - Detlev's advice was that energy released in all
    ! p-reactions would likely be absorbed by the energy states of the system rather than converted into
    ! kinetic energy of the products
    !

    if (debug_kinetics) &
       & write(6,'(a,10(1x,g12.5))') 'DV1:',v%v(1),v%v(2),v%v(3),mag_v(v%v),deltav


    if (reaction_kind.eq.e_reaction) then 
       !
       ! Obtain random unit vector 
       !

       call get_random_vector(vadd)

       !
       ! Calculate the new 3-space velocity vector
       !

       v%v(1) = v%v(1) + deltav * vadd(1)  ! R
       v%v(2) = v%v(2) + deltav * vadd(2)  ! Z
       v%v(3) = v%v(3) + deltav * vadd(3)  ! T

    endif


    if (debug_kinetics) &
       & write(6,'(a,10(1x,g12.5))') 'DV2:',v%v(1),v%v(2),v%v(3),mag_v(v%v),deltav
    if (debug_kinetics) &
       & write(6,'(a,10(1x,g12.5))') 'VADD:',vadd(1),vadd(2),vadd(3),mag_v(vadd),deltav
    if (debug_kinetics) &
       & write(6,'(a,10(1x,g12.5))') 'VADD:',deltav*vadd(1),deltav*vadd(2),deltav*vadd(3),mag_v(vadd),deltav


  end subroutine apply_deltav


  subroutine get_random_vector(v)
    use mod_params
    implicit none
    real :: v(veclen)
    real :: ran1, ran2
    real :: theta,phi

    ! Obtain a random normalized 3 vector
    ! Coordinates of the vector are R,Z,T (tokomak coordinates assuming cylindrical simulation geometry)
    ! Need two random angles - theta in (-PI/2, PI/2) and phi in (0,2PI)
    !
    ! In the coordinate system I am using here - theta is the angle from the RT plane to the vector and 
    ! phi is the angle to the vector projected in the RT plane to the positive R axis.

    call getran(ran1)
    call getran(ran2)

    theta = PI * ran1 - PI/2.0
    phi = 2.0*PI * ran2

    ! R
    v(1) = cos(theta) * cos (phi)
    ! Z
    v(2) = sin(theta)
    ! T 
    v(3) = cos(theta) * sin(phi) 

    ! Note that |v| = 1

  end subroutine get_random_vector



  subroutine get_random_perp_vector(b,bperp)
    use mod_params
    implicit none
    real b(veclen),bperp(veclen)

    !
    ! Calculate a random vector perpendicular to b by first choosing 
    ! a random vector then finding the normalized cross product vector 
    ! between the random vector and B - one needs to check that the 
    ! random vector is not parallel to B. 

    real :: v(veclen)
    logical :: vfound
    integer :: tries
    integer, parameter :: maxtries = 100
    real :: bperp_mag

    tries = 0
    vfound = .false.

    do while (.not.vfound)
       call get_random_vector(v)

       if (dot_product(b,v).eq.1.0) then 
          tries = tries + 1
          if (tries.ge.maxtries) then 
             call errmsg('ERROR in random_bperp_vector: Can not find appropriate random vector','B may be zero?')
             write(6,'(a,3(1x,g12.5),a,3(1x,g12.5))') 'ERROR in random_bperp_vector: B=',b(1),b(2),b(3),' V= ',v(1),v(2),v(3)
             v(1) = b(1) - 1.0
             v(2) = b(2) + 1.0
             v(3) = b(3) 
          else
             cycle
          endif
       endif

       vfound = .true.

    end do

    !
    ! Calculate cross product - this gives a random vector in the plane perpendicular to B
    !

    bperp(1) = b(2) * v(3) - b(3) * v(2) 
    bperp(2) = b(3) * v(1) - b(1) * v(3) 
    bperp(3) = b(1) * v(2) - b(2) * v(1) 

    !
    ! Ensure the vector is normalized
    !

    bperp_mag = sqrt(bperp(1)**2 + bperp(2)**2 + bperp(3)**2)

    bperp = bperp/bperp_mag


    if (debug_kinetics) &
       & write(6,'(a,10(1x,g12.5))') 'BP1:',dot_product(b,bperp),bperp(1),bperp(2),bperp(3),mag_v(bperp),b(1),b(2),b(3),mag_v(b)


  end subroutine get_random_perp_vector








end module hc_kinetics
