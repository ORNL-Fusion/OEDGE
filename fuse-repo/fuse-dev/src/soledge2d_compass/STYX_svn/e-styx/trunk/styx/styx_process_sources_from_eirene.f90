subroutine styx_process_sources_from_eirene(coupling_mode)
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_ccona
  use eirmod_comusr
  use eirmod_comxs
  use styx2eirene
  implicit none  
  integer, intent(in) :: coupling_mode
  integer :: iatm,imol,iion
  integer :: irc,irrc,ipls
  real(dp), allocatable :: Sn(:,:)
  real(dp), allocatable :: Sm(:,:)
  real(dp), allocatable :: SE(:,:)


  !!!  prepare calculation of neutral species temperatures for output (only when diag requested ?)
  !!!  some part needed for preaveraged mode too ...

  atom_density = eov(0)%pdena
  E0_at        = eov(0)%edena
  Gamma0X_at   = eov(0)%vxdena
  Gamma0Y_at   = eov(0)%vydena
  Gamma0Z_at   = eov(0)%vzdena

  mol_density  = eov(0)%pdenm
  E0_mol       = eov(0)%edenm
  Gamma0X_mol  = eov(0)%vxdenm
  Gamma0Y_mol  = eov(0)%vydenm
  Gamma0Z_mol  = eov(0)%vzdenm
        
  tion_density = eov(0)%pdeni
  E0_tion      = eov(0)%edeni
  Gamma0X_tion = eov(0)%vxdeni
  Gamma0Y_tion = eov(0)%vydeni
  Gamma0Z_tion = eov(0)%vzdeni        

! EIRENE yields <E0>, VX, VY, VZ from which T0 and vOpar must be calculated

  call styx_recalculate_neutrals_velocity_fields

! calculate neutral species temperatures (units, eV)

  do iatm=1,natmi
    T0_at(iatm,:)=2._dp/3._dp*E0_at(iatm,:)/elcha & 
           -1._dp/3._dp*rmassa(iatm)*amuakg*(V0X_at(iatm,:)**2+V0Y_at(iatm,:)**2+V0Z_at(iatm,:)**2)/elcha
  enddo
  do imol=1,nmoli
    T0_mol(imol,:)=2._dp/3._dp*E0_mol(imol,:)/elcha &
           -1._dp/3._dp*rmassm(imol)*amuakg*(V0X_mol(imol,:)**2+V0Y_mol(imol,:)**2+V0Z_mol(imol,:)**2)/elcha
  enddo
  do iion=1,nioni
     T0_tion(iion,:)=2._dp/3._dp*E0_tion(iion,:)/elcha & 
           -1._dp/3._dp*rmassi(iion)*amuakg*(V0X_tion(iion,:)**2+V0Y_tion(iion,:)**2+V0Z_tion(iion,:)**2)/elcha
  enddo




  select case (coupling_mode)
    
    case(0)
    
      Sn_at(:,0)        = eov(0)%pael
      SE_at(:,0)        = eov(0)%eael
      Sn_mol(:,0)       = eov(0)%pmel
      SE_mol(:,0)       = eov(0)%emel
      Sn_tion(:,0)      = eov(0)%piel
      SE_tion(:,0)      = eov(0)%eiel
      do ipls=1,nplsi
        Sn_at(:,ipls)  = eov(0)%papl(ipls,:)
        Sn_mol(:,ipls) = eov(0)%pmpl(ipls,:)
        Sn_tion(:,ipls)= eov(0)%pipl(ipls,:)
        Sm_at(:,ipls)  = eov(0)%mapl(ipls,:)
        Sm_mol(:,ipls) = eov(0)%mmpl(ipls,:)
        Sm_tion(:,ipls)= eov(0)%mipl(ipls,:)
        SE_at(:,ipls)  = eov(0)%eaplr(ipls,:)
        SE_mol(:,ipls) = eov(0)%emplr(ipls,:)
        SE_tion(:,ipls)= eov(0)%eiplr(ipls,:)
      enddo

      allocate(Sn(Neir_cells,0:nplsi),Sm(Neir_cells,1:nplsi))
      allocate(SE(Neir_cells,0:nplsi))
 
! conversion from eV.s-1 to W
      EELRC1=EELRC1*elcha

      IRRC=0

       do ipls=1,nplsi
         do irc=1,NRCP(ipls)
           IRRC=IRRC+1

           Sn=0._dp
           Sm=0._dp
           SE=0._dp

           Sn(1:Neir_cells,0)=-1._dp*eiv_i(ipls)%dens(1:Neir_cells)*TABRC1(IRRC,1:Neir_cells)
           SE(1:Neir_cells,0)= eiv_i(ipls)%dens(1:Neir_cells)*EELRC1(IRRC,1:Neir_cells)

           Sn(1:Neir_cells,ipls)=Sn(1:Neir_cells,0)
           Sm(1:Neir_cells,ipls)=(RMASSP(ipls)*amuakg)*vpar_tri(1:Neir_cells,ipls)*Sn(1:Neir_cells,ipls)
           SE(1:Neir_cells,ipls)=-1._dp*(1.5_dp*eiv_i(ipls)%T(1:Neir_cells)+ & 
                                      edrift(ipls,1:Neir_cells))*elcha*Sn(1:Neir_cells,ipls)

           Sn_pls = Sn_pls + Sn
           Sm_pls = Sm_pls + Sm
           SE_pls = SE_pls + SE
         enddo
       enddo

       nrrci=0

       Sn_tot=Sn_at+Sn_mol+Sn_tion+Sn_pls

! remove mass from momentum sources

       do ipls=1,nplsi
         Sm_at(:,ipls)  =Sm_at(:,ipls)/(RMASSP(ipls)*amuakg)
         Sm_mol(:,ipls) =Sm_mol(:,ipls)/(RMASSP(ipls)*amuakg)
         Sm_tion(:,ipls)=Sm_tion(:,ipls)/(RMASSP(ipls)*amuakg)
         Sm_pls(:,ipls)=Sm_pls(:,ipls)/(RMASSP(ipls)*amuakg)
       
! correct signs (Sm from EIRENE has modified signs for plots...)

         Sm_at(:,ipls)  =Sm_at(:,ipls)*sign(1._dp,vpar_tri(:,ipls))
         Sm_mol(:,ipls) =Sm_mol(:,ipls)*sign(1._dp,vpar_tri(:,ipls)) 
         Sm_tion(:,ipls)=Sm_tion(:,ipls)*sign(1._dp,vpar_tri(:,ipls))
         Sm_tot=Sm_at+Sm_mol+Sm_tion+Sm_pls
       enddo
       
       SE_tot=SE_at+SE_mol+SE_tion+SE_pls

    case(1)

       call styx_preav_sources

  end select

  ! and now the radiated power [old option] RECOMBINATION MISSING !! Sn_at instead of Sn_at+Sn_pls

!  if (.not.allocated(Srad_h)) allocate(Srad_h(1:Neir_cells))
!  Srad_h=0._dp
! here all data from amjuel, ionization + recombination
!  Srad_h = SE_at(:,0) + 13.6_dp*elcha*(Sn_at(:,1)+Sn_pls(:,1)) 
!  call styx_total_sources(Srad_h,Srad_h_intg)

  do iatm=1,natmi
    call styx_total_sources(Srad_at(:,iatm),Srad_at_intg(iatm))
  enddo

  do ipls=0,nplsi
    call styx_total_sources(Sn_at(:,ipls),Sn_at_intg(ipls))
    call styx_total_sources(Sn_mol(:,ipls),Sn_mol_intg(ipls))
    call styx_total_sources(Sn_tion(:,ipls),Sn_tion_intg(ipls)) 
    call styx_total_sources(Sn_pls(:,ipls),Sn_pls_intg(ipls))
    call styx_total_sources(SE_at(:,ipls),SE_at_intg(ipls))
    call styx_total_sources(SE_mol(:,ipls),SE_mol_intg(ipls))
    call styx_total_sources(SE_tion(:,ipls),SE_tion_intg(ipls))
    call styx_total_sources(SE_pls(:,ipls),SE_pls_intg(ipls))
  enddo

! integrals of momentum sources, has to be done after sign correction
  do ipls=1,nplsi
    call styx_total_sources(Sm_at(:,ipls),Sm_at_intg(ipls))
    call styx_total_sources(Sm_mol(:,ipls),Sm_mol_intg(ipls))
    call styx_total_sources(Sm_tion(:,ipls),Sm_tion_intg(ipls))
    call styx_total_sources(Sm_pls(:,ipls),Sm_pls_intg(ipls))
  enddo





end subroutine styx_process_sources_from_eirene
