c       -*-Fortran-*-
c
c
      subroutine global_hc_init
      use comhc
      use hc_kinetics_options
      implicit none
c
      include 'hc_global_opts'
c


! ammod begin.
c -----------------------------------------------------------------------
c
c     Hydrocarbon Common Block Defaults
c
c -----------------------------------------------------------------------
c
c     Initialize ComHC common block arrays to zeros.
c
      Call RZero (hc_reflection_coefs, maxnds)
      Call RZero (hc_sticking_coefs, maxnds)
c
c -----------------------------------------------------------------------
c
c     Hydrocarbon Option Defaults
c
c -----------------------------------------------------------------------
c
c     TAG H15 (Integer)
c     hc_follow_option - Activates hydrocarbon following routine (Call HC) 
c                        instead of LAUNCH from initial point in neutbatch
c                        subroutine in neut.d6a.  Set to 'off' by default
c                        unless optional tag exists in input file.
c                        hc_follow_option = 0 - not activated
c                        hc_follow_option= 1 - activated
c
      hc_follow_option = 0
c
c     Also initialize the globa_hc_follow_option in comtor
c 
      global_hc_follow_option = 0

c
c -----------------------------------------------------------------------
c
c     TAG H16 (Integer)
c     hc_higher_hcs_option - Determines whether to follow evolution of hydro-
c                            carbons beyond Methane (CH4) - species 10.
c                            hc_higher_hcs_option = 0 - do not follow HCs larger than CH4
c                            hc_higher_hcs_option = 1 - follow HCs up to those available
c                                                 in data library (C3H8)
c
      hc_higher_hcs_option = 0
c
c -----------------------------------------------------------------------
c
c     TAG H17 (Integer)
c     hc_wbc_comp_option - Determines whether to attempt to match WBC output
c                          by following particles until a defined boundary is
c                          struck, then printing a WBC-like output chart for
c                          particles followed.
c                          hc_wbc_comp_option = 0 - do not implement comparison (default)
c                          hc_wbc_comp_option = 1 - use boundaries and print WBC table
c
      hc_wbc_comp_option = 0
c
c -----------------------------------------------------------------------
c
c     TAG H20 (Integer)
c     hc_sputtering_model - Model used for initial chemical sputtering species release.
c                      hc_sputtering_model = 0 - preset HC species used
c                      hc_sputtering_model = 1 - Mech, Haasz, Davis model
c                                                (nuc.mat., 1998).
c
      hc_sputtering_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H21 (Integer)
c     hc_sputtered_hc_species - Preset hydrocarbon species to be sputtered
c                               upon chemical sputtering event.
c                               (10 = Methane, CH4 is common default).
c
      hc_sputtered_hc_species = 10
c
c -----------------------------------------------------------------------
c
c     TAG H22 (Integer)
c     hc_evolution_model_primary - Primary dataset to use for reaction
c                                  rate profiles (1, E&L is the default).
c                      .
      hc_evolution_model_primary = 1
c
c     This set of reaction rate data will be read first by the data loading
c     routine for use in the code.
c
c -----------------------------------------------------------------------
c
c     TAG H23 (Integer)
c
c     hc_evolution_model_secondary - Secondary dataset to use for reaction
c                                    rate profiles (0, none is the default).
c                      .
      hc_evolution_model_secondary = 0
c
c     This set of reaction rate data will be read second, and used only to
c     compliment reaction rates missing from the primary set.
c
c -----------------------------------------------------------------------
c
c     TAG H24 (Integer)
c     hc_launch_location - Type of launch from various locations inside the
c                          vessel.  Options include:
c                         -1: Equal to CNEUTB 
c                          0: Distributed launch along target.
c                          1: Given at R,Z.
c                          2: Homogeneously along walls.
c                          3: Distributed along target.
c                          4: distributed along wall surfaces.
c                          5: 2D distributed launch from cell centres on the grid.
c                       ** 6: Launch from one specified wall index (see H60)
c                          The defualt is the value in the DIVIMP input file.
c
c     Note: This initialization here is USELESS since the values of 
c           CNEUTB and CNEUTC for the present simulation have NOT yet been read 
c           in the input file.      
c
c      hc_launch_location = CNEUTB
c
      hc_launch_location = -1
c
c -----------------------------------------------------------------------
c
c     TAG H25 (Integer)
c     hc_launch_angle_velocity - Type of launch angle and velocity option.
c                                Options are the same as the DIVIMP Vel/angle
c                                flag (CNEUTC) and the default is the CNEUTC
c                                value in the DIVIMP input file. Additionally:
c                                -1: Equal to CNEUTC
c
c      hc_launch_angle_velocity = CNEUTC
c
      hc_launch_angle_velocity = -1
c
c -----------------------------------------------------------------------
c
c     TAG H26 (Integer)
c     hc_launch_angle_velocity - Type of launch velocity option.  This option
c                                overtakes that of H25 if a value of 1 or 2
c                                (ie. a MB distribution) is used.  If H26
c                                is 0 (ie. a constant velocity of injection)
c                                then the option for H25 will be applied.
c                                Options include:
c                                0: Constant launch velocity.
c                                1: Single MB distribution with surface temp mean.
c                                2: Dual MB distribution with secondary described
c                                   in H27/H28.
c
      hc_launch_velocity_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H27 (Real)
c     hc_dual_mb_pri_vel_flux - Primary velocity flux. Range: 0.0 - 1.0.
c                               Note:  A value >0.6 is not recommended.
c                               for any surface temperature <800 deg K.
c                               See Vietzke, Nuc.Mat. 290-293 (2001) 158-161.
c
      hc_dual_mb_pri_vel_flux = 0.6
c
c -----------------------------------------------------------------------
c
c     TAG H28 (Real)
c     hc_dual_mb_sec_mean_temp - Mean temperature of the second, hotter MB
c                                distribution which is added to the primary
c                                associated with the surface temperature.
c                                Range: 0.0 - 2000.0 Deg. K
c                                Note:  A value >2,000 is not recommended
c                                for any surface temperature <800 deg K.
c                                See Vietzke, Nuc.Mat. 290-293 (2001) 158-161.
c
      hc_dual_mb_sec_mean_temp = 2000.0
c
c -----------------------------------------------------------------------
c
c     TAG H30 (Integer)
c     hc_neut_ion_velocity - Neutral->ion transition initial velocity 
c                            Options are the same as the DIVIMP Initial Ion
c                            velocity option flag (CNEUTG) and the default
c                            is the CNEUTG value in the DIVIMP input file.
c                            -1: Equal to CNEUTG
c
c      hc_neut_ion_velocity = CNEUTG
c
      hc_neut_ion_velocity = -1
c
c -----------------------------------------------------------------------
c
c     TAG H31 (Integer)
c     hc_ion_neut_angle - Ion->neutral transition angle option.
c                         Available options include:
c                         0 - isentropic distribution upon fragment
c                             neutralization
c                         1 - sine distribution biased parallel to
c                             field line at point of neutralization
c                         2 - direction of S at point of neutralization
c
      hc_ion_neut_angle = 0
c
c -----------------------------------------------------------------------
c
c     TAG H32 (Integer)
c     hc_ion_neut_velocity - Ion->neutral transition initial velocity
c                            option.  Available options include:
c                            0 - equal to ion at point of neutralization (default)
c                            1 - 
c
      hc_ion_neut_velocity = 0
c
c -----------------------------------------------------------------------
c
c     TAG H33 (Integer)
c     hc_lambda_calc - Use lambda=15 everywhere, or improved calculation
c                      of Sivukhin ("Coulomb collisions in a full ionized
c                      plasma", Reviews of Physics, Vol.4, 1966).
c                      hc_lambda_calc = 0 - use lambda = 15 everywhere (default)
c                      hc_lambda_calc = 1 - use calculated lambda for low T
c                                           L = 30.0-1/2ln(n*)+1.5ln(Ti)
c                                           Note use of Ti as ion-ion collisions
c                                           dominate FF,FiG,FPG forces (Stangby pg. 301).
c                                   
      hc_lambda_calc = 0
c
c -----------------------------------------------------------------------
c
c     TAG H34 (Integer)
c     hc_disable_transitions - Turn off ability for initial hydrocarbon 
c                              to evolve to lower hydrocarbons.
c                              hc_disable_transitions = 0 - off (default)
c                              hc_disable_transitions = 1 - on
c                                   
      hc_disable_transitions = 0
c
c -----------------------------------------------------------------------
c
c     TAG H35 (Integer)
c     hc_presheath_efield - Use improved model of electric field of Brooks
c                           ("Near-surface sputtered particle transport",
c                           Phys.Fluids B 2 (8), August 1990) 
c                           hc_presheath_efield = 0 - off (default)
c                           hc_presheath_efield = 1 - on
c                                   
      hc_presheath_efield = 0
c
c -----------------------------------------------------------------------
c
c     TAG H36 (Real)
c     hc_efield_drop_fraction - Fraction of potential drop in Debye region
c                               hc_efield_drop_fraction = 0.0-1.0 (0.25 typical)
c
      hc_efield_drop_fraction = 0.25
c
c -----------------------------------------------------------------------
c
c     TAG H37 (Integer)
c     hc_efield_cells - Number of cells from the divertor plate to apply
c                       the Brooks model for near-sheath E-field.
c                       hc_efield_drop_fraction = 0-5 (1 typical)
c
      hc_efield_cells = 1
c
c -----------------------------------------------------------------------
c
c     TAG H40 (Integer)
c     hc_neutral_reflection_option - Switch to turn on and off reflection of neutral
c                             impacting particles.
c                             0 - off (default)
c                             1 - on
c
      hc_neutral_reflection_option = 0
c
c -----------------------------------------------------------------------
c
c     TAG H41 (Integer)
c     hc_ion_reflection_option - Switch to turn on and off reflection of ionized
c                         impacting particles.
c                         0 - off (default)
c                         1 - on
c
      hc_ion_reflection_option = 0
c
c -----------------------------------------------------------------------
c
c     TAG H42 (Integer)
c     hc_reflection_coef_model - Method to assign or calculate reflection
c                                coefficient once particle has run into wall
c                                or divertor.
c                              hc_reflection_coef_model = 0 - preset sticking
c                                                         coefficient (default)
c                                                         (0.5 or half is common).
c                              hc_reflection_coef_model = 1 - Janev, et al. model
c                                                         (NIFS, 2001).
c                              hc_reflection_coef_model = 2 - Alman, Ruzic model
c                                                         (PSI, 2002).
c
      hc_reflection_coef_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H43 (Real)
c     hc_reflection_coef_preset - Preset reflection coefficient for both
c                                 neutral and ionized impacting hydrocarbons.
c                                 hc_reflection_coef_preset = 0.0-1.0 (0.5 default)
c                                 Note that for HC's, one cannot tell the
c                                 difference between reflection and sputtering.
c
      hc_reflection_coef_preset = 0.5
c
c -----------------------------------------------------------------------
c  
c     TAG H44 (Integer)
c     hc_reflection_species_model - Decide which species an impinging
c                           hydrocarbon becomes after reflecting
c                           off the vessel wall or divertor.
c                           hc_reflection_species_model = 0 - use pre-initalized 
c                                                         reflecting table set to 1
c                                                         increase in H of C species
c                           hc_reflection_species_model = 1 - use Alman, Ruzic model
c                                                         (PSI, 2002).
c
      hc_reflection_species_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H45 (Integer)
c     hc_reflection_energy_model - Calculate or assign energy of released
c                                  particle after reflection from the
c                                  wall or divertor.
c                                  hc_reflection_energy_model = 0 - preset.
c                                  hc_reflection_energy_model = 1 - impact energy.
c                                  hc_reflection_energy_model = 2 - substrate thermal.
c                                  hc_reflection_energy_model = 3 - Alman, Ruzic
c                                                                   (PSI, 2002).
c
      hc_reflection_energy_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H46 (Real)
c     hc_refl_energy_neutral_preset - Hydrocarbon molecular energy after a neutral
c                                     reflecting from vessel wall or divertor
c                                     (typically 1-10 eV for detached plasma).
c
      hc_refl_energy_neutral_preset = 2.0
c
c -----------------------------------------------------------------------
c
c     TAG H47 (Real)
c     hc_refl_energy_ion_preset - Hydrocarbon molecular energy after an ion
c                                 reflecting (as a neutral) from vessel wall
c                                 or divertor (typically higher for detached
c                                 plasma).
c
      hc_refl_energy_ion_preset = 30.0
c
c -----------------------------------------------------------------------
c
c     TAG H48 (Integer)
c     hc_reflection_angle_model - Model to decide at what angle to the 
c                                 perpendicular an ejected neutral 
c                                 particle is released at when not stuck
c                                 to the wall or divertor. 
c                                 
c                                 hc_reflection_angle_model =-1 - NRFOPT (default).
c                                 hc_reflection_angle_model = 0 - off.
c                                 hc_reflection_angle_model = 1 - specular.
c                                 hc_reflection_angle_model = 2 - isentropic.
c                                 hc_reflection_angle_model = 3 - normal.
c                                 hc_reflection_angle_model = 4 - use Alman
c                                                                 & Ruzic
c                                                                 (PSI, 2002).
c
c      hc_reflection_angle_model = NRFOPT
c
      hc_reflection_angle_model = -1
c   
c -----------------------------------------------------------------------
c
c     TAG H50 (Integer)
c     hc_sputtering_option - Switch to turn on and off sputtering of
c                            impacting hydrocarbon particles.
c                            0 - off
c                            1 - on (default)
c
      hc_sputtering_option = 1
c
c -----------------------------------------------------------------------
c
c     TAG H51 (Integer)
c     hc_sticking_coef_model - Method to assign or calculate sticking
c                              coefficient once particle has run into wall
c                              or divertor.
c                              hc_sticking_coef_model = 0 - preset sticking
c                                                      coefficient (default)
c                                                      (0.5 or half is common).
c                              hc_sticking_coef_model = 1 - Janev, et al. model
c                                                      (NIFS, 2001).
c                              hc_sticking_coef_model = 2 - Alman, Ruzic model
c                                                      (PSI, 2002).
c
      hc_sticking_coef_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H52 (Real)
c     hc_sticking_coef_preset - Preset default sticking coefficient used
c                               instead of experimental or theoretical model.
c				hc_sticking_coef_preset = -1.0 - CTRESH 
c                                           (DIVIMP sticking threshold CTRESH
c                                           beyond which the particle sputters)
c				hc_sticking_coef_preset = 0.0-1.0 (default 0.5)
c
      hc_sticking_coef_preset = 0.5
c
c -----------------------------------------------------------------------
c
c     TAG H53 (Integer)
c     hc_sputtering_species_model - Decide which species an impinging
c                           hydrocarbon becomes after self-sputtering
c                           off the vessel wall or divertor.
c                           hc_sputtering_species_model = 0 - use pre-initalized 
c                                                         self-sputtering table set to 1
c                                                         increase in H of C species
c                           hc_sputtering_species_model = 1 - use Alman, Ruzic model
c                                                         (PSI, 2002).
c
      hc_sputtering_species_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H54 (Integer)
c     hc_sputtering_energy_model - Calculate or assign energy of released
c                                  particle after sputtering from the
c                                  wall or divertor.
c                                  hc_sputtering_energy_model = 0 - preset.
c                                  hc_sputtering_energy_model = 1 - impact energy.
c                                  hc_sputtering_energy_model = 2 - substrate thermal.
c                                  hc_sputtering_energy_model = 3 - Alman, Ruzic
c                                                                   (PSI, 2002).
c
      hc_sputtering_energy_model = 0
c
c -----------------------------------------------------------------------
c
c     TAG H55 (Real)
c     hc_sput_energy_neutral_preset - Hydrocarbon molecular energy after a netural
c                                     sputtering from vessel wall or divertor
c                                     (typically 1-10 eV for detached plasma).
c
      hc_sput_energy_neutral_preset = 2.0
c
c -----------------------------------------------------------------------
c
c     TAG H56 (Real)
c     hc_sput_energy_ion_preset - Hydrocarbon molecular energy after an ion
c                                 sputtering from vessel wall or divertor
c                                 (typically 1-10 eV for detached plasma).
c
      hc_sput_energy_ion_preset = 2.0
c
c -----------------------------------------------------------------------
c
c     TAG H57 (Integer)
c     hc_sputtering_angle_model - Model to decide at what angle to the 
c                                 perpendicular an ejected neutral 
c                                 particle is released at when not stuck
c                                 to the wall or divertor.                              
c                                 hc_sputtering_angle_model =-1 - NRFOPT
c                                 hc_sputtering_angle_model = 0 - off.
c                                 hc_sputtering_angle_model = 1 - specular.
c                                 hc_sputtering_angle_model = 2 - isentropic.
c                                 hc_sputtering_angle_model = 3 - normal.
c                                 hc_sputtering_angle_model = 4 - SQRT(sin) (default).
c
      hc_sputtering_angle_model = 4
c
c -----------------------------------------------------------------------
c
c     TAG H60 (Integer)
c
c     Wall segment launch index for HC Launch option 6
c
c     The default value is 52 which is applicable for one case only - the
c     value of this quantity should always be specified in the input if 
c     hc_launch_location is set to 6. 
c
      hc_launch6_wall_index=52
      hc_launch6_wall_index_set=.false.
c
c -----------------------------------------------------------------------
c
c     TAG H61 (Integer)
c
c     HC reaction kinetics option 
c
c     The default value is 1 which will result in the additional 
c     kinetic energy for the HC reaction as specified in the HC data 
c     input being applied to the reaction products.
c     
c     NOTE: This option only affects Adam's original kinetics code. In the
c           new code excluding the reaction energy is not an option. 
c           hc_kintics_opt = 0 selects the original code.
c
      hc_reaction_kinetics = 1
c
c -----------------------------------------------------------------------
c
c     TAG H62 (Integer)
c
c     HC reaction kinetics model selection option 
c
c     This option chooses the model to be used for the HC reaction
c     kinetics. 
c        Option 0: Adam's old code 
c        Option 1: 3D velocity following code for ions and neutrals
c
c     Set the default to the new 3D kinetics code. 
c
      hc_kinetics_opt = 1
c
c -----------------------------------------------------------------------
c
c     TAG H63 (Integer)
c
c     HC perpendicular velocity option applying to hc_kinetics_opt=1
c
c     This option selects the option to be used when evaluating the
c     ion velocity perpendicular to the magnetic field during the 
c     transitions. 
c     Option 0: Vperp does not change between HC reactions
c     Option 1: Vperp is set equal to Vpara assuming collisional heating
c               of some sort comparable to vpara
c     Option 2: Vperp estimated based on separately evolving Tperp
c     Option 3: Perpendicular velocity experiences velocity diffusion 
c               independently of the parallel velocity.
c     Option 4: Vperp experiences exactly the same velocity diffusion 
c               as Vpara - but are otherwise uncoupled. 
c
      hc_vperp_opt = 0
c
c -----------------------------------------------------------------------
c
c     jdemod
c
c     TAG H64 (Real)
c     input_HC_H_mass        - Mass of H isotopes associated with the 
c                              HC molecules. Default value is set to 
c                              1.0 for hydrogen. 
c                            - This value is set to crmb if it is not
c                              given an explicit value in the input file
c
      input_HC_H_mass = -1.0
c
c -----------------------------------------------------------------------
c
c     Hydrocarbon output options.
c
c -----------------------------------------------------------------------
c
c     TAG H90 (Integer)
c     hc_coord_print_option - Prints r,z coordinates for hydrocarbon being
c                             followed for each timestep.
c                             hc_coord_print_option = 0 - not activated
c                             hc_coord_print_option = 1 - activated
c
      hc_coord_print_option = 0
c
c -----------------------------------------------------------------------
c
c     TAG H91 (Integer)
c     hc_evolve_print_option - Prints r,z coordinates and HC species to-
c                              from data for each hydrocarbon transition step.
c                              hc_evolve_print_option = 0 - not activated
c                              hc_evolve_print_option = 1 - activated
c
      hc_evolve_print_option = 0
c
c -----------------------------------------------------------------------
c
c
c     End of hydrocarbon-related intialization
c      
! ammod end.

      return
      end
c
c
c
      subroutine global_hc_assign_inputs
      use comhc
      implicit none
      include 'comtor'
c
c     For HC input values that can be set to a value requiring the 
c     assignment of the base DIVIMP input values - these values
c     need to be assigned after the input file is read in. 
c -----------------------------------------------------------------------
c
c     H24 -
c
      if(hc_launch_location==-1) hc_launch_location=CNEUTB
c
c -----------------------------------------------------------------------
c     H25 -
c
      if(hc_launch_angle_velocity==-1) hc_launch_angle_velocity=CNEUTC
c
c -----------------------------------------------------------------------
c     H30 -
c
      if(hc_neut_ion_velocity==-1) hc_neut_ion_velocity = CNEUTG
c
c -----------------------------------------------------------------------
c     H48 - neutral reflection angle option
c
c     NOTE: A reflection option of 0 is NOT valid in the HC code - however, this value is 
c           used to turn reflection off in the main code. SO - if NRFOPT is zero - assign
c           a default value of 2 (cosine reflection)
c
      if (hc_reflection_angle_model==-1) then 
         if (nrfopt.eq.0) then 
            hc_reflection_angle_model = 2
            write(0,*) 'H48 WARNING: HC neutral reflection model'//
     >                 ' set to 2 (cosine) since DIVIMP neutral'//
     >                 ' reflection is OFF'
            write(6,*) 'H48 WARNING: HC neutral reflection model'//
     >                 ' set to 2 (cosine) since DIVIMP neutral'//
     >                 ' reflection is OFF'
         else
            hc_reflection_angle_model = nrfopt
         endif
      endif

c
c -----------------------------------------------------------------------
c     H57 - ion sputtered fragment angle option
c
c     NOTE: A sputtering angle option of 0 is NOT valid in the HC code - however, this value is 
c           used to turn reflection off in the main code. SO - if NRFOPT is zero - assign
c           a default value of 2. (cosine reflection)
c
      if(hc_sputtering_angle_model==-1) then 
         if (nrfopt.eq.0) then 
            hc_sputtering_angle_model = 2
            write(0,*) 'H57 WARNING: HC sputtering angle model'//
     >                 ' set to 2 (cosine) since DIVIMP neutral'//
     >                 ' reflection is OFF'
            write(6,*) 'H57 WARNING: HC sputtering angle model'//
     >                 ' set to 2 (cosine) since DIVIMP neutral'//
     >                 ' reflection is OFF'
         else
            hc_sputtering_angle_model = nrfopt
         endif
      endif
c
c -----------------------------------------------------------------------
c     H60 - HC launch option 6 - wall segment index for launch
c
c     Check to see that if the launch option specified is option 6 - then
c     the wall index must also be set or issue an error message
c
      if (hc_launch_location.eq.6.and.
! slmod begin
     >    .not.hc_launch6_wall_index_set ) then 
!
!     >    hc_launch6_wall_index_set.eq..false.) then 
! slmod end
c          
            write(0,*) 'H24 AND H60: WARNING: HC launch option 6 '//
     >                 ' specified in H24 but wall segment index'//
     >                 ' for launch NOT specified by H60'
            write(6,*) 'H24 AND H60: WARNING: HC launch option 6 '//
     >                 ' specified in H24 but wall segment index'//
     >                 ' for launch NOT specified by H60'
c
      endif
c
c     jdemod 
c
c     H64 : input_HC_H_mass
c
c     If this value has not been assigned through unstructured input - set it
c     equal to the background plasma ion mass specified in the input file
c
      if (input_HC_H_mass == -1.0) then 
         input_HC_H_mass = crmb
      endif
c
c     jdemod - assign the value of cprint to hc_cprint for use inside the HC routines
c
      hc_cprint = cprint
c


      return
      end
c
c
c
      subroutine hc_electric_field_mod(ik,ir,iz,s,local_electric_field)
      use comhc
      Use HC_Utilities ! Sheath E-field calc by Brooks.
      implicit none
c
      include 'cgeom'      
      include 'comtor'
c
      integer ik,ir,iz,id
      real s
c
c     Local variables
c
      Real :: Local_Electric_Field, Ediv, Efact
      integer it

! ammod begin.
      ! Modify electric field strength if within applicable range of target.
      ! Check for applicability of sheath electric field modification.
      ! Note:  Default for HC_Presheath_Efield is 0 or off.
      ! Do only for SOL rings (IR .ge. IRSEP).
      if (HC_Presheath_Efield .eq. 1 .and. IR .ge. IRSEP) Then
         ! Check if within magnetic pre-sheath zone.
         ! Note:  Default for HC_EField_Cells is 1.
         If (MIN (ik, nks (ir) - ik) .le. HC_EField_Cells) Then
c
c           Find nearest target element - assume for now it is closest target
c           on current ring.
c
            if (ik.lt.nks(ir)/2) then 
               it = 2
            else
               it = 1
            endif
c
c           Define target index
c
            id = idds(ir,it) 
c
            ! Add Brooks modification to local electric field strength.
            ! Note:  Default for E-field drop fraction is 0.25.
c
c           jdemod - calculate the DIVIMP E-field in V/m - take out scaling factors
c
            efact =  QTIM * QTIM * EMI / CRMI

            Ediv = kes(ik,ir) / efact
c
            Local_Electric_Field = Sheath_E_Field (ik,ir,s,
     >        kteds(id), ktebs (ik,ir), ktibs (ik,ir),
     >        knbs (ik,ir), crmb, rizb,
     >        HC_EField_Drop_Fraction, bts (ik,ir), Ediv)
c
c           jdemod - *** BUG ***
c                    the sheath_e_field function can return kes if
c                    the Brooks Efield calculation is less than the 
c                    DIVIMP estimate. KES already contains the scaling 
c                    factors of QTIM - the following line which rescales
c                    the Efield is wrong when the routine returns the 
c                    DIVIMP Efield value 
c                  - fixed the bug by passing in a revised Ediv that could
c                    be directly compared to the calculated sheath value 
c
               Local_Electric_Field = Local_Electric_Field * efact            
c
         Else
            Local_Electric_Field = kes (ik,ir)
         End If
      Else
         ! Use standard DIVIMP electric field.
         Local_Electric_Field = kes (ik,ir)
      End If

! ammod end.

      return
      end
c
c
c
      subroutine global_hc_wbc_comp(iz,crmi,vel,temi,sputy)
      use comhc
      Use HC_WBC_Comp ! Records ion death statistics for WBC comparison.
      implicit none
      integer iz
      real crmi,vel,temi,sputy

! ammod begin.       
      ! WBC comparison addition for ion struck target.
      Call Record_WBC_Ion_Event (IZ,CRMI,VEL,TEMI,SPUTY)
! ammod end.	       

      return
      end
c
c     
c
      subroutine global_hc_check_WBC_Ion_Pos(IK,IR,S,CROSS,IZ,CRMI,
     >                      VEL,TEMI,SPUTY,NIMPS,NIMPS2,NATIZ,IMP,flag)
      use comhc
      use HC_WBC_Comp
      implicit none
      integer ik,ir,iz,nimps,nimps2,natiz,imp,flag
      real s,cross,crmi,vel,temi,sputy
c
      flag = 0
c
      flag= Check_WBC_Ion_Position (IK,IR,S,CROSS,IZ,CRMI,VEL,
     >    TEMI,SPUTY,NIMPS,NIMPS2,NATIZ,IMP)

      return
      end
c
c
c
      subroutine global_hc_check_WBC_Neut_Pos(R,Z,CRMI,VIN,TEMN,
     >                          SPUTY,
     >                          NPROD,LPROD,IPROD,CIST,flag)
      use comhc
      use HC_WBC_Comp
      implicit none
      integer nprod,lprod,iprod,flag
      real r,z,crmi,vin,temn,sputy,cist
c
      flag = 0
c
      flag= Check_WBC_Neut_Position (R,Z,CRMI,VIN,TEMN,SPUTY,
     >                             NPROD,LPROD,IPROD,CIST)
c
      return
      end
c
c
c
      subroutine global_hc_begin
      use comhc
      use hc_start
      implicit none

!     Initialize all data which needs to be set only once in the case.
      Call HC_Begin ()

      return
      end
c
c
c
      subroutine global_hc_end(NIMPS,NIMPS2,REXIT,RMAIN,RNEUT,NEUTIM,
     >                         NIZS)
      use comhc
      use hc_start
      implicit none
      integer nimps,nimps2,nizs
      real rexit,rmain,rneut,neutim

      ! Complete launch, follow and recording of all particles.
      ! Print hydrocarbon diagnostics.
      Call HC_End (NIMPS,NIMPS2,REXIT,RMAIN,RNEUT,NEUTIM,NIZS)

      return
      end
c
c
c      
      subroutine global_hc_launch(NN1,NN2,NI1,NI2,
     >     RSTRUKA,MTCSTRUKA,RMAINA,REXITA,
     >     RATIZA,RNEUTA,RWALLNA,MTCWALLNA,RCENTA,RTMAXA,
     >     SEED,NRAND,NEUTIM,RFAILA,STATUS,MATP,MATT,
     >     neuttype,cneutb,cneutc)
c
      use comhc
      use hc_batch 
c
      implicit none
c
      integer nn1,nn2,ni1,ni2,neuttype,matp,matt,
     >        cneutc,cneutb,status,nrand
      real     RSTRUKa,RMAINa,REXITa,RATIZa,RNEUTa,RWALLNa,
     >         RCENTa,RTMAXa,rfaila
      real     mtcstruka, mtcwallna,neutim  
      double precision seed 
c
      write(0,*) 'Calling HC_LAUNCH' 
c
           Call HC_Launch (NN1,NN2,NI1,NI2,
     >     RSTRUKA,MTCSTRUKA,RMAINA,REXITA,
     >     RATIZA,RNEUTA,RWALLNA,MTCWALLNA,RCENTA,RTMAXA,
     >     SEED,NRAND,NEUTIM,RFAILA,STATUS,MATP,MATT,
     >     neuttype,cneutb,cneutc)
c
      return
      end 

c
c
c
      subroutine global_hc_store_raw_data
      Use ComHC ! HC constants.
      Use HC_Init_DIV_Diag ! Included to re-set hc_state, hc_density, hc_output, Number_HC_Species.
      implicit none
c
! ammod begin.
c     Addition of hydrocarbon data handling.
c
      call iinout ('W HC_OPT',hc_follow_option,1)
      call iinout ('W HC_HIG',hc_higher_hcs_option,1)
      call iinout ('W HC_PRI',hc_evolution_model_primary,1)
      call iinout ('W HC_SEC',hc_evolution_model_secondary,1)

      call rinout ('W HC_STA', HC_State_List,
     >  Number_HC_Species)
      call dinout ('W HC_DEN', HC_Density,
     >  maxnks*maxnrs*(Number_HC_Species))
      !
      ! jdemod - changed last paramter from number_h_species to num_h_states since that 
      !          is what the array is declared
      !
      call dinout ('W H_DEN', H_Density,
     >  maxnks*maxnrs*(Num_H_States))
      call iinout ('W HC_OUT', HC_Output_List,
     >  Number_HC_Species+1)
      call rinout ('W HC_WLK', HC_Walks,
     >  Max_Number_Walks*2)
      !
      ! jdemod - remove the extra storage in hc_factor (0 and -1) indices
      !
      call rinout ('W HC_FACTA', HC_Factor_A,
     >  Number_HC_Species)
      call rinout ('W HC_FACTB', HC_Factor_B,
     >  Number_HC_Species)
c      call rinout ('W HC_FACTA', HC_Factor_A,
c     >  Number_HC_Species+2)
c      call rinout ('W HC_FACTB', HC_Factor_B,
c     >  Number_HC_Species+2)
      call rinout ('W HC_TIZS_CH', HC_TIZS_CH,
     >  maxnks*maxnrs*2)
c      call rinout ('W HC_TIZS_C2', HC_TIZS_C2,
c     >  maxnks*maxnrs*2)
c
c      CALL RINOUT ('W FYTOT ',FYTOT ,1)      
c
c
c     End addition of hydrocarbon data to BIN unit 8 file for read in OUT.
! ammod end.

      return 
      end
c
c
c
      subroutine update_cross(ik,ir,ikold,irold,kk,s,theta,cross,
     >                        adjust,dcross,ckkmin,smax,k,
     >                        nrand,imp,cist,debug)         
      implicit none 
      integer ik,ir,ikold,irold,kk
      real s,theta,adjust,dcross(4),cross,ckkmin,smax,k
      real*8 cist
      integer nrand, imp
      logical debug
c
c*************************************************************************
c
c     UPDATE_CROSS: This routine replaces the code that was originally
c                   in the DIV module to implement the cross-field
c                   transport. The base option uses the code that was 
c                   in DIV - the other options implement different 
c                   algorithms or methods for calculating the cross-field
c                   transport. The routine takes all of the particle 
c                   position related quantites from DIV as it's input. In
c                   addition some variables that are used to accumulate 
c                   statistics are also passed in order to keep the 
c                   existing code as intact as possible.          
c
c
c     David Elder, Jan 22, 1997
c
c*************************************************************************
c
      include    'params'
      include    'comtor'
      include    'cgeom'
      include    'crand' 
c
c     Local variables
c
      real kprob,theta1,pinchvel
      real kdrefin,kdrefout
      real kratin,kratout
      real sdperptmp,tmpran 
      real kperpstepin,kperpstepout
      integer jk,ikreftmp
      integer flag
c
      real tmpout,tmpin,crossfrac,oldcross,tmptheta,oldtheta
      real tmpcross,tmpcross2
      integer iktmp,irtmp,jktmp
C
      logical vr_assigned
      integer ierr
c

c
C
C-------- FIND NEAREST IK CORRESPONDING TO DISTANCE S ALONG CONTOUR IR
C-------- FIND THETA VALUE CORRESPONDING TO DISTANCE S ALONG CONTOUR IR
C-------- ADJUST CROSS FIELD TERM FOR NEW DISTANCES BETWEEN CONTOURS
C
          IKOLD = IK
          IROLD = IR
c
c         Find new cell
c
          IF (PDOPT.EQ.1) THEN
  620       IF (IK.LT.NKS(IR).AND.S.GT.KSB(IK,IR)) THEN
              IK = IK + 1
              GOTO 620
            ENDIF
  625       IF (IK.GT.1.AND.S.LT.KSB(IK-1,IR)) THEN
              IK = IK - 1
              GOTO 625
            ENDIF
          ELSEIF (PDOPT.eq.0) then 
  630       IF (IK.LT.NKS(IR).AND.S.GT.KSS(IK,IR)) THEN
              IK = IK + 1
              GOTO 630
            ENDIF
  635       IF (IK.GT.1.AND.S.LE.KSS(IK-1,IR)) THEN
              IK = IK - 1
              GOTO 635
            ENDIF
            IF (IK.GT.1.AND.S-KSS(IK-1,IR).LT.KSS(IK,IR)-S) IK = IK - 1   
          ENDIF
	!if (ik .ne. ikold) then
	!  write (0,*) "MOVED cfield!!!!",ik,ikold
	!end if
c
c         If using non-orthogonal transport - find theta value
c
c         NOTE: For the core rings the first and last cells are the same.
c               S=0 starts at the cell centre and S=SMAX is at the same
c               cell centre - thus - in the core - the particle should never
c               have an S-value with S < KSS(1,ir) or S > KSS(nks(ir),ir)
c
c

          IF (S.GT.KSS(IK,IR)) THEN
            IF (IK.LT.NKS(IR).or.ir.lt.irsep) THEN
              THETA = THETAG(IK,IR) + (S-KSS(IK,IR))/KFORDS(IK,IR)
     >                                *(THETAG(IK+1,IR)-THETAG(IK,IR))
            ELSE
              THETA = THETAG(IK,IR) + (S-KSS(IK,IR))/KFORDS(IK,IR)
     >                              *(THETAT(IDDS(IR,1))-THETAG(IK,IR))
            ENDIF
          ELSE
            IF (IK.GT.1.or.ir.lt.irsep) THEN
c
c             This was added to cover the case of S=0 occuring for an ion initially
c             injected in the first cell of a core ring.  
c
              if (s.eq.kss(ik,ir)) then
                 theta = thetag(ik,ir)
              else
                THETA = THETAG(IK,IR) - (KSS(IK,IR)-S)/KBACDS(IK,IR)
     >                                *(THETAG(IK,IR)-THETAG(IK-1,IR))
              endif  
            ELSE
              THETA = THETAG(IK,IR) - (KSS(IK,IR)-S)/KBACDS(IK,IR)
     >                              *(THETAG(IK,IR)-THETAT(IDDS(IR,2)))
            ENDIF
          ENDIF
c
c
c         Adjust cross-field term - if necessary
c
          call adjust_cross(cross,adjust,ik,ir,ikold,irold,debug) 
c
c         Record some statistics in the core. 
c
          if (ir.lt.irsep.and.adjust.ne.0.0) then 
             IF (ADJUST.LE.0.0) THEN
                DCROSS(1) = DCROSS(1) + 1.0
                DCROSS(2) = DCROSS(2) + ADJUST
             ELSE
                DCROSS(3) = DCROSS(3) + 1.0
                DCROSS(4) = DCROSS(4) + ADJUST
             ENDIF
          endif
C
          IF (DEBUG) WRITE (6,1001) 'D1A:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,
     >      'UPDATED S'

C
C-------- UPDATE CROSS FIELD DIFFUSION. 
c
c         Record position
c
          IKOLD = IK
          IROLD = IR
          oldcross = cross
          oldtheta = theta 
c
c         Perform Cross-field diffusive step. - various methods
c
c         Applying different methods depending on options. 
c
c
c        if (debug) write(6,*) '3d:',ik,ir,s,cross,theta,adjust
          

c
c         Calculate the probability of inward step and perform 
c         the appropriate cross-field step.
c
c---------------------------------------------------------------
c
c         Regular cross-field diffusion - cioptj = 0,1 ...
c
       if (cioptj.eq.0.or.cioptj.eq.1) then 
c
c
          if (cdiffopt.eq.0.or.
     >        (cdiffopt.eq.1.and.ir.ge.irsep)) then 
c
            kprob = kins(ik,ir) 
c
c
          elseif (cdiffopt.eq.1.and.ir.lt.irsep) then 
c
c
            if ((ir.eq.irsep-1.and.cross.le.0.0).or.ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI.or.tdistin(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI.or.tdistout(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.2) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI.or.tdistin(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI.or.tdistout(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.3) then 
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI.or.distin(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/distin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/distin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI.or.distout(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/distout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/distout(ik,ir)))
              endif
c
c              write(6,*) 'Kpout:',ik,ir,kprat2(ik,ir,2),distout(ik,ir)
c

c
            endif 
c
          endif

c
          call set_pinch_velocity(ik,ir,nrand,imp,cist,pinchvel,
     >                              vr_assigned,ierr)



c
c         Pinch Velocity
c
c          pinchvel = 0.0 
c
c          if (pinchopt.eq.0.0) then 
c
c             pinchvel = 0.0
c
c          elseif (pinchopt.eq.1.or.
c     >           (pinchopt.eq.2.and.
c     >           (ir.ge.irsep.and.ir.le.irwall))) then 
c
c             pinchvel = vpinch
c 
c          elseif (pinchopt.eq.3) then 
c
c             pinchvel = kpinchs(ik,ir)
c
c          endif
c
c         Update Cross
c
c
c          write (6,*) 'Kprob:',ik,ir,kprob,cdiffopt,cross
c
c
          KK = KK +1 

          if (pinchopt.eq.4.and.vr_assigned) then 
             CROSS = CROSS + PINCHVEL
          else
             CROSS = CROSS + SIGN (KPERPS(IK,IR), kprob-ranv(kk))
     >                      + PINCHVEL
          endif

c
c          CROSS = CROSS + SIGN (KPERPS(IK,IR), kprob-ranv(kk))
c     >                      + PINCHVEL
c
          tmpran = ranv(kk)
c
c---------------------------------------------------------------
c
c         Cioptj = 2 - reference line proportional Dperp
c
c         This option provides an alternative way of calculating 
c         "CROSS" - the cross-field displacement. 
c
c
       elseif (cioptj.eq.2) then
c
c         Calculate step direction probability
c
c
c         Calculate the ACTUAL size of the Dperp step that will
c         be taken.
c
            if (ir.lt.irsep) then 
c
               kdrefin = tdistin(ikrefcore,ir)
               kdrefout= tdistout(ikrefcore,ir)
               sdperptmp = sdperpref
c
            elseif (ir.ge.irsep.and.ir.le.irwall) then 
c
               kdrefin =  tdistin(ikrefsol,ir)
               kdrefout= tdistout(ikrefsol,ir)
               sdperptmp  = sdperpref
c
            elseif (ir.ge.irtrap.and.ir.le.nrs) then 
c
               kdrefin = tdistin(ikrefpp,ir)
               kdrefout= tdistout(ikrefpp,ir)
c
c              Keep the Dperp in the Private Plasma constant on
c              each set of knots matching the corresponding cell 
c              on the separatrix so that a pinch does not develop - 
c              This matches the procedure that appears to have been
c              used in EDGE2D.
c
               ikreftmp= ikouts(ik,nrs) 
c
               sdperptmp=sdperpref
     >               *(distin(ikreftmp,irsep)
     >               /(distin(ikrefsol,irsep)))
c
c               sdperptmp  = sdperppp
c
            endif
c
            if (kdrefout.eq.0.0) then 
               kratout = 1.0
            else   
               kratout = tdistout(ik,ir)/kdrefout
            endif
c
            if (kdrefin.eq.0.0) then 
               kratin = 1.0
            else   
               kratin = tdistin(ik,ir)/kdrefin
            endif
c
c           These are estimates necessary for approximating the step probability
c
c           The estimate process in the Private plasma keeps the Dperp constant
c           across a set of knots starting on the separatrix. This matches what
c           is done in EDGE2D and will avoid pinch effects but does not seem very 
c           physical.  
c   
            if (ir.ge.irtrap.and.ir.le.nrs) then 

               kperpstepin = sdperptmp 
               kperpstepout = sdperptmp 
               kratin = 1.0
               kratout = 1.0

            else  

               kperpstepin = sdperptmp * kratin
               kperpstepout = sdperptmp* kratout

            endif
c
c 
c
          if (cdiffopt.eq.0.or.
     >        (cdiffopt.eq.1.and.ir.ge.irsep)) then 
c
            kprob = kins(ik,ir) 
c
          elseif (cdiffopt.eq.1.and.ir.lt.irsep) then 
c
            if ((ir.eq.irsep-1.and.cross.le.0.0).or.ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperpstepin/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperpstepout/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistin(ik,ir)))
              endif
c
            endif 
c
         elseif (cdiffopt.eq.2) then 
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperpstepin/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperpstepout/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.3) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperpstepin/2.0)/distin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/distin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperpstepout/2.0)/distout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/distout(ik,ir)))
              endif
c
            endif 
c
          endif

c
          call set_pinch_velocity(ik,ir,nrand,imp,cist,pinchvel,
     >                              vr_assigned,ierr)

c
c         Pinch Velocity
c
c          pinchvel = 0.0 
c
c          if (pinchopt.eq.0.0) then 
c
c             pinchvel = 0.0
c
c          elseif (pinchopt.eq.1.or.
c     >           (pinchopt.eq.2.and.
c     >           (ir.ge.irsep.and.ir.le.irwall))) then 
c
c             pinchvel = vpinch
c 
c          elseif (pinchopt.eq.3) then 
c
c             pinchvel = kpinchs(ik,ir)
c
c          endif
c  
c         Map to cross on reference line
c
          if (cross.gt.0.0) then 
c
             tmpcross = cross / kratin
c
          elseif (cross.le.0.0) then 
c
             tmpcross = cross / kratout
c
          endif
c
c
c         Update TmpCross
c
          kk = kk +1
c
          if (pinchopt.eq.4.and.vr_assigned) then 
             tmpCROSS = tmpCROSS + PINCHVEL
          else
             tmpCROSS = tmpCROSS + SIGN (sdperptmp,kprob-ranv(kk))
     >                      + PINCHVEL
          endif
c
c          tmpCROSS = tmpCROSS + SIGN (sdperptmp,kprob-ranv(kk))
c     >                      + PINCHVEL
c
c
c         Map cross back using the reference line ratios again
c
          if (cross.gt.0.0) then 
c
             cross = tmpcross * kratin
c
          elseif (cross.le.0.0) then 
c
             cross = tmpcross * kratout
c
          endif
c
c
c
c     Spatially varying Dperp - more complicated 
c
      elseif (cioptj.eq.3.or.cioptj.eq.4) then 
c
c         Uses Kperps(ik,ir) with spatially
c         varying values.  
c
c
          if (cdiffopt.eq.0.or.
     >        (cdiffopt.eq.1.and.ir.ge.irsep)) then 
c
            kprob = kins(ik,ir) 
c
c
          elseif (cdiffopt.eq.1.and.ir.lt.irsep) then 
c
c
            if ((ir.eq.irsep-1.and.cross.le.0.0).or.ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir) ))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.2) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.3) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/distin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/distin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/distout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/distout(ik,ir)))
              endif
c
            endif 
c
          endif

c
          call set_pinch_velocity(ik,ir,nrand,imp,cist,pinchvel,
     >                              vr_assigned,ierr)


c
c         Pinch Velocity
c
c          pinchvel = 0.0 
c
c          if (pinchopt.eq.0.0) then 
c
c             pinchvel = 0.0
c
c          elseif (pinchopt.eq.1.or.
c     >           (pinchopt.eq.2.and.
c     >           (ir.ge.irsep.and.ir.le.irwall))) then 
c
c             pinchvel = vpinch
c 
c          elseif (pinchopt.eq.3) then 
c
c             pinchvel = kpinchs(ik,ir)
c
c          endif
c
c
c         This is where the process gets fancy - check to
c         see if the Dperp step will carry cross into the 
c         next cell - if it does - use a portion of the
c         step-size in the next cell to calculate the 
c         final position of the particle. I.E. Use the
c         differeing cross-field step sizes in the 
c         adjoining cells in proportion to amounts already
c         used.
c
c         e.g. for example if a dperp step carries a particle
c              across the boundary of cell A at a distance
c              that is 40% of a cell A step - then beyond the
c              boundary - the particle will be moved 60% of
c              a cell B step.  
c
c         The Dperp step is done first and THEN a pinch velocity
c         (if any) is added on. 
c          
          oldcross = cross
c
          kk = kk +1
          tmpCROSS = CROSS + SIGN (KPERPS(IK,IR), kprob-RANV(KK))
c
c
c         Check to see if has crossed a cell boundary. 
c
c         Check for only ONE boundary for now - if the particle
c         is crossing multiple rings with one cross-field step
c         then the time-step is too large. However, there is 
c         code that checks this condition and prints an error
c         message. 
c
c
c         First - find out exactly what cell the particle will
c         be in - so we can use the appropriate Dperp for the 
c         second portion of the step.
c          
c
          tmpcross2 = tmpcross
          tmptheta = theta
          iktmp = ik
          irtmp = ir
          jktmp = jk
c
c         Use the DO_CFSTEP routine to find the new cell into which the 
c         particle will be transported with it's current value of cross. 
c
          call do_cfstep(jktmp,iktmp,irtmp,irold,tmpcross2,adjust,
     >                   tmptheta,flag,debug)
c
          if (flag.gt.0) then 
c
c            Particle moved inward - use the Dperp of the cell 
c            it would have wound up in to calculate the final 
c            value of cross. Can double check this if one wants
c            to by calling do-cfstep again with the corrected
c            value of cross. 
c      
c            Note: cosali and cosalo have been set to 1.0 for all
c            cells on grids that have been assumed orthogonal. 
c
             crossfrac = (distin(ik,ir)- oldcross) / kperps(ik,ir)
c
             tmpcross =  distin(ik,ir)
     >                  + (1.0-crossfrac) * kperps(iktmp,irtmp)
c 
             if (flag.gt.1) then 
                write (6,*) 'ERROR: Cross-field step error -' 
                write (6,*) '       Crossed multiple rings = ',flag
                write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
             end if 
c
c            Info for debugging purposes ...
c             
             tmpout = distout(iktmp,irtmp)
             tmpin = distin(ik,ir)
c
c
c
c             write (6,'(a,5g13.6)') 'DperpIN:',tmpcross,tmpin,
c     >                 tmpcross2,
c     >                 tmpcross-tdistin(ik,ir),
c     >                 tmpout
c             write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
c
          elseif (flag.lt.0) then  
c      
c            Note: cosali and cosalo have been set to 1.0 for all
c            cells on grids that have been assumed orthogonal. 
c
             crossfrac = abs(-distout(ik,ir)-
     >                     oldcross) / kperps(ik,ir)
c
             tmpcross = -distout(ik,ir)
     >                  - (1.0-crossfrac) * kperps(iktmp,irtmp)
c 
             if (flag.lt.-1) then 
                write (6,*) 'ERROR: Cross-field step error -' 
                write (6,*) '       Crossed multiple rings = ',flag
                write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
             end if 
c
c            Info for debugging purposes ...
c             
             tmpout = distout(iktmp,irtmp)
             tmpin =  distin(ik,ir)
c
c             write (6,'(a,5g13.6)') 'DperpIN:',tmpcross,tmpin,
c     >                 tmpcross2,
c     >                 tmpcross-tdistin(ik,ir),
c     >                 tmpout
c             write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
c
          endif
c
          if (flag.eq.0.and.(ir.ne.irtmp.or.ik.ne.iktmp)) then 
             
             write (6,*) 'FLAG=0:IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp

          endif    
c
c         Finalize the updating of cross by adding the pinch velocity
c
          if (pinchopt.eq.4.and.vr_assigned) then 
             CROSS = cross + PINCHVEL
          else
             CROSS = tmpcross + PINCHVEL
          endif
c
c          CROSS = TMPCROSS + PINCHVEL
c
c
c         End of cioptj IF statement
c
       endif
c
c
C
          IF (DEBUG) WRITE (6,1000) 'D2A:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,kprob,
     >      kprob-tmpran,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'

c
c         Move the particle across rings if the CROSS value
c         has exceeded the distance to the next ring. Recalculate
c         THETA if using NON-orthogonal transport
C
          call do_cfstep(jk,ik,ir,irold,cross,adjust,theta,flag,debug)
c
c         Record some statistics in the core. 
c
          if (ir.lt.irsep.and.adjust.ne.0.0) then 
             IF (ADJUST.LE.0.0) THEN
                DCROSS(1) = DCROSS(1) + 1.0
                DCROSS(2) = DCROSS(2) + ADJUST
            ELSE
                DCROSS(3) = DCROSS(3) + 1.0
                DCROSS(4) = DCROSS(4) + ADJUST
             ENDIF
          endif
C
          IF (DEBUG) WRITE (6,1000) 'D3A:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,kprob,
     >      kprob-tmpran,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'
c
c         Re-calculate S if the particle has changed rings.  
c
          IF (IR.NE.IROLD) THEN
            K      = KKS(IR)
            CKKMIN = MIN (CKKMIN, K)
            SMAX   = KSMAXS(IR)
c 
c           ITER grid 
c
            if (cgridopt.eq.2) then
               if ( (((ir.ge.irsep.and.ir.le.irwall2).or.
     >               (ir.ge.irsep2.and.ir.le.irwall)).and.
     >               (irold.ge.irtrap)).or.
     >              (((irold.ge.irsep.and.irold.le.irwall2).or.
     >               (irold.ge.irsep2.and.irold.le.irwall)).and.
     >               (ir.ge.irtrap))) then
                  IF (S.GT.KSS(IKOLD,IROLD)) THEN
                    S = KSS(IK,IR) - (S-KSS(IKOLD,IROLD)) *
     >                         (KBACDS(IK,IR)/KFORDS(IKOLD,IROLD))
                  ELSE
                    S = KSS(IK,IR) + (KSS(IKOLD,IROLD)-S) *
     >                         (KFORDS(IK,IR)/KBACDS(IKOLD,IROLD))
                  ENDIF
               elseIF (S.GT.KSS(IKOLD,IROLD)) THEN
                    S = KSS(IK,IR) + (S-KSS(IKOLD,IROLD)) *
     >                         (KFORDS(IK,IR)/KFORDS(IKOLD,IROLD))
               ELSE
                  S = KSS(IK,IR) - (KSS(IKOLD,IROLD)-S) *
     >                         (KBACDS(IK,IR)/KBACDS(IKOLD,IROLD))
               ENDIF
c
c           Orthogonal Transport 
c
            ELSEIF (NORTHOPT.EQ.0.or.northopt.eq.2) THEN
c
              IF (S.GT.KSS(IKOLD,IROLD)) THEN
                  S = KSS(IK,IR) + (S-KSS(IKOLD,IROLD)) *
     >                         (KFORDS(IK,IR)/KFORDS(IKOLD,IROLD))
              ELSE
                 S = KSS(IK,IR) - (KSS(IKOLD,IROLD)-S) *
     >                         (KBACDS(IK,IR)/KBACDS(IKOLD,IROLD))
              ENDIF
c
c
c           Non-orthogonal Transport
c
c           This block generates a new value for S after non-orthogonal 
c           cross-field diffusion (that results in changing rings).
c
            ELSEIF (NORTHOPT.EQ.1.or.northopt.eq.3) THEN
c
c             Adjustment to Theta moved to DO_CFSTEP where
c             code decides new IK value of cell.
c
c             Special for the last knot on a core ring:
c
c              IF (IKOLD.EQ.NKS(IROLD).AND.IR.GE.IRSEP
c     >            .and.irold.lt.irsep) THEN
c                THETA = THETAG(IK,IR) - 
c     >                  (THETAG(IK,IR) - THETAG(IK-1,IR)            ) *
c     >                  (THETAG(IKOLD,IROLD) - THETA                ) /
c     >                  (THETAG(IKOLD,IROLD) - THETAG(IKOLD-1,IROLD))
c              ENDIF
c
              IF (IR.LT.IRSEP.AND.IK.EQ.1.AND.
     +            THETA.LT.THETAG(IK,IR)) THEN

                IK = NKS(IR)
c
c               Handle an error condition when IKOLD was also 1.
c               The spacing from nks(ir) to nks(ir) -1 on core rings
c               is the same as cell 1 to a mythical cell 0. since cell
c               1 and nks(ir) coincide. IKOLD should not be 1 for a 
c               non-core 
c
c
                if (ikold.eq.1) then 
                   THETA = THETAG(IK,IR) - 
     +                  (THETAG(IK,IR)       - THETAG(IK-1,IR)      ) *
     +                  (THETAG(IKOLD,IROLD) - THETA                ) /
     +          (THETAG(nks(irold),IROLD) - THETAG(nks(irold)-1,IROLD))
c
c               Regular case
c
                else 
                   THETA = THETAG(IK,IR) - 
     +                  (THETAG(IK,IR)       - THETAG(IK-1,IR)      ) *
     +                  (THETAG(IKOLD,IROLD) - THETA                ) /
     +                  (THETAG(IKOLD,IROLD) - THETAG(IKOLD-1,IROLD))
                endif               


          IF (DEBUG) WRITE (6,1000) 'D4:',IKold,IRold,S,K,
     >      THETA,thetag(ik,ir),thetag(ik-1,ir),
     >      thetag(ikold,irold),thetag(ikold-1,irold),
     >      kprob-tmpran,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'



              ENDIF


          IF (DEBUG) WRITE (6,1000) 'D5:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,kprob,
     >      kprob-tmpran,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'

c
c             Re-calculate THETA and S.
c
c
c             First Half of cell
c
              IF (THETA.LT.THETAG(IK,IR)) THEN
c
                IF (IK.EQ.1) THEN
                  IF (IR.LT.IRSEP) THEN
                    IK = NKS(IR)             
 
                    THETA = THETA + (THETAG(IK,IR) - THETAG(1,IR))

                    THETA1 = (THETAG(IK,IR) - THETA          ) /
     >                       (THETAG(IK,IR) - THETAG(IK-1,IR))
                    S = KSS(IK,IR) - KBACDS(IK,IR) * THETA1
                  ELSE
                    IF( THETA.LE.THETAT(IDDS(IR,2)) )THEN
c
c                     Particle has struck target cross-field 
c
c                     If Mirror target option is ON - place particle
c                     back in old ring with cross set to it's
c                     previous value.
c
c                     Otherwise do the usual. 
          IF (DEBUG) WRITE (6,'(a,2i6,10(1x,g20.12))') 'D5A IK=1:',
     >      IK,IR,S,K,
     >      THETA,SMAX,CROSS,thetat(idds(ir,2))


c
                      if (cmiropt.eq.0) then 

                         S = 0.0

                      elseif (cmiropt.eq.1) then                      
c
c                        Do not change S - reset CROSS and IR
c
                         cross = oldcross 
                         theta = oldtheta  
                         ir    = irold
                         ik    = ikold
c
                      endif
c
                    ELSE
                      S = KSS(IK,IR) * 
     >                    (THETA         - THETAT(IDDS(IR,2))) /
     >                    (THETAG(IK,IR) - THETAT(IDDS(IR,2)))
                    ENDIF
                  ENDIF
                ELSE
                  THETA1 = (THETAG(IK,IR) - THETA          ) /
     >                     (THETAG(IK,IR) - THETAG(IK-1,IR))
                  S = KSS(IK,IR) - KBACDS(IK,IR) * THETA1
                ENDIF
c
c             Particle in second half of cell. 
c
              ELSE
c
                IF( IK.EQ.NKS(IR) )THEN
                  IF (IR.LT.IRSEP) THEN
                    IK = 1
                    THETA = THETA - (THETAG(NKS(IR),IR) - THETAG(1,IR))
  
                    THETA1 = (THETA           - THETAG(IK,IR)) /
     >                       (THETAG(IK+1,IR) - THETAG(IK,IR))
                    S = KSS(IK,IR) + KFORDS(IK,IR) * THETA1
                  ELSE
                    IF (THETA.GE.THETAT(IDDS(IR,1))) THEN

c
c                     Particle has struck target cross-field 
c
c                     If Mirror target option is ON - place particle
c                     back in old ring with cross set to it's
c                     previous value.
c
          IF (DEBUG) WRITE (6,'(a,2i6,10(1x,g20.12))') 'D5A IK=NKS:',
     >      IK,IR,S,K,
     >      THETA,SMAX,CROSS,thetat(idds(ir,1))


c                     Otherwise do the usual. 
c
                      if (cmiropt.eq.0) then 

                         S = SMAX

                      elseif (cmiropt.eq.1) then                      
c
c                        Do not change S - reset CROSS and IR
c
                         cross = oldcross 
                         theta = oldtheta
                         ir    = irold
                         ik    = ikold 
c
                      endif
c
                    ELSE
                      S = KSS(IK,IR) + KFORDS(IK,IR) *
     >                    (THETAT(IDDS(IR,1)) - THETA        ) /
     >                    (THETAT(IDDS(IR,1)) - THETAG(IK,IR))
                    ENDIF
                  ENDIF
                ELSE
                  THETA1 = (THETA           - THETAG(IK,IR)) /
     >                     (THETAG(IK+1,IR) - THETAG(IK,IR))
                  S = KSS(IK,IR) + KFORDS(IK,IR) * THETA1
                ENDIF
c
              ENDIF
c    
            ENDIF


          IF (DEBUG) WRITE (6,1000) 'D6:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,kprob,
     >      kprob-tmpran,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'



          ENDIF

c        if (debug) write(6,*) '3d:',ik,ir,cross,s


c
c     End of routine 
c
      return
c
c     Format statements
c
 1000 format(a,2i4,1p,10(g11.4),1x,a) 
 1001 format(a,2i4,1p,6(g11.4),44x,1x,a) 
c
      end
c
c
c
      subroutine hc_update_line_profile(ik,ir,r,z,vr,vz,sputy)

      !                                 cion,rizb)

      implicit none
      integer ik,ir
      real vr,vz,sputy,r,z

      ! Note that cion and rizb are not included in the argument list. Since we are dealing with HC - the 
      ! only possible value of CION corresponds to carbon - which will be loaded into CION in the main 
      ! simulation - so these data are loaded from common blocks and passed to the update_line_profile routine. 
      ! In general however, the update_line_profile routine can be used for any species by changing the value of 
      ! CION. 

      include 'params'
      include 'comtor'
      include 'line_profile'
      
      
      !
      ! Return if the line profile option is not active
      !
      if (line_profile_opt.eq.0) return

      !
      ! Scaling of vr and vz is performed on the call to hc_update_line_profile
      !

      call update_line_profile(ik,ir,r,z,vr,vz,sputy,cion,rizb)


      return
      end

