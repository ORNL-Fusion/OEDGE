module styx2eirene
  use eirmod_precision
  use eirmod_comxs

  logical :: is_3D
  character(12) :: fluid_code
  integer :: Ntor_cells,Neir_cells
  ! full torus = 360
  real(dp) :: ang_max 

  integer :: numiter,curiter
  logical :: all_tal

  type :: surfaces
  	integer :: itri
  	integer :: iside
  	integer :: iprop
  	real(dp) :: R1
  	real(dp) :: R2
  	real(dp) :: Z1
  	real(dp) :: Z2
  	integer :: v1
  	integer :: v2
  	real(dp) :: ds
  end type surfaces

  type(surfaces), allocatable :: surface(:),recsurf(:)
  real(dp), allocatable :: Rmin(:),Rmax(:),Zmin(:),Zmax(:)
  integer, allocatable :: smin(:),smax(:)
  integer, allocatable :: ksurf(:),krecsurf(:),ke2r(:),kr2e(:)
  integer, allocatable :: recsurfinv(:,:)
  real(dp), allocatable :: ssurf(:)
  real(dp) :: Rmaxg,Rming,Zmaxg,Zming

  real(dp), allocatable :: pflux_in(:,:,:)
  
  integer, allocatable :: map_species_index(:)
  

  type :: sheath_data
    real(dp) :: alphaB
    real(dp) :: uparX
    real(dp) :: uparY
    real(dp) :: uparZ
    real(dp) :: tau
    real(dp) :: ksi0
    real(dp) :: ksi
    real(dp) :: alpha_V
    real(dp) :: beta_V

    integer :: hit_V

    real(dp), allocatable :: E(:,:)
    real(dp), allocatable :: alpha(:,:)
    real(dp), allocatable :: beta(:,:)
  
  end type sheath_data


  type :: sheath1D_database
    integer  :: NalphaB
    integer  :: Ntau
    integer  :: Nksi
    real(dp), allocatable :: alphaB(:)
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: ksi(:)
    real(dp), allocatable :: E(:,:,:)
    real(dp), allocatable :: alpha(:,:,:)
    real(dp), allocatable :: beta(:,:,:)
  end type sheath1D_database

  type(sheath_data), allocatable :: sheath1D(:)
  type(sheath1D_database) :: sheath1D_av

  real(dp), allocatable :: Bnu2xyz(:,:,:),xyz2Bnu(:,:,:)

  integer :: n_pfc_types

  type :: pfc_model
    character(2) :: material
    integer :: iatm
    integer :: sputer_model
    real(dp) :: sputer_yield_phys
    real(dp) :: sputer_yield_chem
    logical :: sputer_on
    real(dp) :: T
    real(dp), allocatable :: R(:)
  end type pfc_model

  type(pfc_model), allocatable :: pfc_models(:)

  integer :: n_mat

  type :: material
    character(2) :: symbol
    integer :: iatm,ipfc
  end type material

  integer :: ispc_add

  type(material), allocatable :: materials(:)

  integer :: n_pumps

  type :: pump
    character(2) :: material 
    real(dp), allocatable :: R(:)
    real(dp) :: T
    real(dp) :: Surface ! in m2
    real(dp),allocatable :: species_Temperature(:) ! in K
    real(dp) :: pumping_Speed ! in m3/s
    integer,allocatable :: triangle_index(:)
    real(dp),allocatable :: triangle_weight(:)
    integer :: Ntriangles
    logical :: isSpeedSet
  end type pump

  type(pump), allocatable :: pumps(:)

  integer :: sheath_model

  type :: puff
   character(2) :: species
   real(dp) :: rate
   real(dp) :: T0
   real(dp) :: divergence
   integer :: itri
   integer :: iside
   integer :: itor
  end type puff

  type(puff), allocatable :: puffs(:)
  
  real(dp) :: total_puff 
  integer, allocatable :: singly_charged_ions(:)
 
  type :: eirene_in_vol
    real(dp), allocatable :: dens(:)
    real(dp), allocatable :: parmom(:)
    real(dp), allocatable :: vx(:)
    real(dp), allocatable :: vy(:)
    real(dp), allocatable :: vz(:)
    real(dp), allocatable :: T(:)
  end type eirene_in_vol
  

  type :: eirene_out_vol
    ! pre-averaged coupling data
    real(dp), allocatable :: pdena(:,:)
    real(dp), allocatable :: pdenm(:,:)
    real(dp), allocatable :: pdeni(:,:)
    real(dp), allocatable :: edena(:,:)
    real(dp), allocatable :: edenm(:,:)
    real(dp), allocatable :: edeni(:,:)
    real(dp), allocatable :: vxdena(:,:)
    real(dp), allocatable :: vxdenm(:,:)
    real(dp), allocatable :: vxdeni(:,:)
    real(dp), allocatable :: vydena(:,:)
    real(dp), allocatable :: vydenm(:,:)
    real(dp), allocatable :: vydeni(:,:)
    real(dp), allocatable :: vzdena(:,:)
    real(dp), allocatable :: vzdenm(:,:)
    real(dp), allocatable :: vzdeni(:,:)
          
    ! direct coupling data
    real(dp), allocatable :: pael(:)
    real(dp), allocatable :: pmel(:)
    real(dp), allocatable :: piel(:)
    real(dp), allocatable :: papl(:,:)
    real(dp), allocatable :: pmpl(:,:)
    real(dp), allocatable :: pipl(:,:)
    real(dp), allocatable :: eael(:)
    real(dp), allocatable :: emel(:)
    real(dp), allocatable :: eiel(:)
    real(dp), allocatable :: eapl(:)
    real(dp), allocatable :: empl(:)
    real(dp), allocatable :: eipl(:)
    real(dp), allocatable :: eaplr(:,:)
    real(dp), allocatable :: emplr(:,:)
    real(dp), allocatable :: eiplr(:,:)
    real(dp), allocatable :: mapl(:,:)
    real(dp), allocatable :: mmpl(:,:)
    real(dp), allocatable :: mipl(:,:)
 
  end type eirene_out_vol


  type :: eirene_out_surf
     real(dp), allocatable :: potat(:,:)
     real(dp), allocatable :: prfaat(:,:)
     real(dp), allocatable :: prfmat(:,:)
     real(dp), allocatable :: prfpat(:,:)
     real(dp), allocatable :: potml(:,:)
     real(dp), allocatable :: prfaml(:,:)
     real(dp), allocatable :: prfmml(:,:)
     real(dp), allocatable :: prfpml(:,:)
     real(dp), allocatable :: potio(:,:)
     real(dp), allocatable :: potpl(:,:)
     real(dp), allocatable :: eotat(:,:)
     real(dp), allocatable :: erfaat(:,:)
     real(dp), allocatable :: erfmat(:,:)
     real(dp), allocatable :: erfiat(:,:)
     real(dp), allocatable :: erfpat(:,:)
     real(dp), allocatable :: eotml(:,:)
     real(dp), allocatable :: erfaml(:,:)
     real(dp), allocatable :: erfmml(:,:)
     real(dp), allocatable :: erfiml(:,:)
     real(dp), allocatable :: erfpml(:,:)
     real(dp), allocatable :: eotio(:,:)
     real(dp), allocatable :: erfaio(:,:)
     real(dp), allocatable :: erfmio(:,:)
     real(dp), allocatable :: erfiio(:,:)
     real(dp), allocatable :: erfpio(:,:)
     real(dp), allocatable :: eotpl(:,:)
     ! sputtering data
     real(dp), allocatable :: sptat(:,:)
     real(dp), allocatable :: sptml(:,:)
     real(dp), allocatable :: sptio(:,:)
     real(dp), allocatable :: sptpl(:,:)
     real(dp), allocatable :: spttot(:)


     real(dp), allocatable :: spump(:,:)

  end type eirene_out_surf

  type(eirene_in_vol), allocatable :: eiv_i(:)
  type(eirene_in_vol) :: eiv_e
  type(eirene_out_vol), save, allocatable :: eov(:)
  type(eirene_out_surf), save, allocatable :: eos(:)

  real(dp) :: Temin_eirene

  integer, allocatable :: italv(:)

  real(dp), allocatable :: Delta_Omega(:,:)
  real(dp), allocatable :: Radflx(:)


  character(8), allocatable :: texts_styx(:)

  integer :: NX,NY,NRKNOT_styx,NTRI_styx,NQUADS
  integer :: Nsou

  integer :: NPUFFS,NSTRATA
  integer, allocatable :: Npart_eirene(:),seed_eirene(:)
  integer :: n_short_cycles,ns_refresh,nrefresh

  real(dp), allocatable :: scresc(:,:,:)

  integer, allocatable :: XMCP_styx(:)
  integer :: nprs_styx
  real(dp), allocatable :: time_para(:)
  real(dp) :: timpara
  real(dp) :: ptimax,ptimin
  integer :: rtimax,rtimin
  real(dp) :: ptmean

  integer :: stratum_plot,nhist_plot 
      
  integer, allocatable :: NVERT(:,:), NEIGH(:,:),NSIDE(:,:),IPROP(:,:)
  integer, allocatable :: NTRIVOIS(:), TRIVOIS(:,:)
 
  integer :: NTTRA_styx
  real(dp) :: ROA_styx
      
  real(dp), allocatable :: NL(:)
  real(dp), allocatable :: TeL(:)
  real(dp), allocatable :: TiL(:,:)
  real(dp), allocatable :: TeLm(:),TeLp(:)
  real(dp), allocatable :: TiLm(:,:),TiLp(:,:)
  
  real(dp), allocatable :: Gammapar(:)
  real(dp), allocatable :: vpar_styx(:),vpar_tri(:,:)
  real(dp), allocatable :: Gammapar_save(:)

  real(dp), allocatable :: BX_tri(:)
  real(dp), allocatable :: BY_tri(:)
  real(dp), allocatable :: BZ_tri(:)
  real(dp), allocatable :: BF_tri(:)

! Eflx_calc_pen = 0 : Te,Ti calculated by interpolation of extrapolated values (send fluxes)
!               = 1 : Te,Ti calculated such that energy conservation in ion channel ensured

  integer, parameter :: Eflx_calc_pen=1

  logical :: tweak_chemistry !  T: allow tweaks of chemistry
  logical :: hardwired       !  T: use fast atomic data calculation (equivalent to general one but faster for standard input files)
  logical :: interface_test_mode

  character(4), allocatable :: reac_switch(:),reac_switchE(:)
  real(dp), allocatable :: fakt(:),faktE(:)
  real(dp) :: fa1,fa2,fm1,fm2,fm3,fi1,fi2,fi3,fp1
  real(dp) :: addtls

  logical :: rad_mat

  ! this is the flux send to EIRENE, from routine send_fluxes
  real(dp), allocatable :: recflux_save(:),puff_save(:)
  ! this is the flux strength in EIRENE (part/s), strata 1
  real(dp) :: fluxt_eirene
  ! HUGO BUFFERAND
  real(dp) flux_neutral_out
  ! END HUGO BUFFERAND

  real(dp) :: RescaleE_save
 
  real(dp) :: Sn_tot_eirene,SE_tot_eirene,SEi_tot_eirene,SEi_in
  real(dp) :: n_tot_styx
  real(dp), allocatable :: vol_tri_eirene(:)   
  real(dp) :: voltot_eirene
      
  real(dp) :: v0 
  integer :: source_species

  real(dp) :: LX,LY,LZ
  real(dp), allocatable :: xtriang(:),ytriang(:)

  real(dp) :: dt_eirene
  integer :: ITN

  logical :: transfer_debug
  logical :: enter_short_cycle
  logical :: cx_checks
  logical :: direct_coupling
  logical :: timedep
  logical :: internal_energy
  logical :: short_cycle

  integer :: levgeo_styx,sc_level
  integer :: am_database

  integer :: ntime_styx
  integer :: nprnli_styx
  integer :: nptst_styx 
  integer :: ntmstp_styx
  integer :: nsnvi_styx
  integer :: icalleir

  real(dp) :: dtimv_styx,time0_styx

  integer :: nscycles

  integer :: n_call_eir
  integer :: seed_gap


  logical, allocatable :: acc_set(:)


  integer :: Nrecyc,Nrecomb

  integer :: NSURF_TAL,nsurf0
  integer :: ipropmin,ipropmax

  real(dp), allocatable :: wv2c(:,:)

  real(dp), allocatable :: atom_density(:,:)

  real(dp), allocatable :: Gamma0X_at(:,:),Gamma0Y_at(:,:),Gamma0Z_at(:,:)
  real(dp), allocatable :: V0X_at(:,:),V0Y_at(:,:),V0Z_at(:,:)
  real(dp), allocatable :: Gammapar0_at(:,:)
  real(dp), allocatable :: v0par_at(:,:)
  real(dp), allocatable :: E0_at(:,:),T0_at(:,:)

  real(dp), allocatable :: mol_density(:,:)

  real(dp), allocatable :: Gamma0X_mol(:,:),Gamma0Y_mol(:,:),Gamma0Z_mol(:,:)
  real(dp), allocatable :: V0X_mol(:,:),V0Y_mol(:,:),V0Z_mol(:,:)
  real(dp), allocatable :: Gammapar0_mol(:,:)
  real(dp), allocatable :: v0par_mol(:,:)
  real(dp), allocatable :: E0_mol(:,:),T0_mol(:,:)

  real(dp), allocatable :: tion_density(:,:)

  real(dp), allocatable :: Gamma0X_tion(:,:),Gamma0Y_tion(:,:),Gamma0Z_tion(:,:)
  real(dp), allocatable :: V0X_tion(:,:),V0Y_tion(:,:),V0Z_tion(:,:)
  real(dp), allocatable :: Gammapar0_tion(:,:)
  real(dp), allocatable :: v0par_tion(:,:)
  real(dp), allocatable :: E0_tion(:,:),T0_tion(:,:)

  real(dp), allocatable :: Sn_intg(:),Sm_intg(:),SE_intg(:)

  real(dp), allocatable :: Sn_at(:,:)
  real(dp), allocatable :: Sm_at(:,:)
  real(dp), allocatable :: SE_at(:,:)
  real(dp), allocatable :: Srad_at(:,:),Srad_at_intg(:)

  real(dp), allocatable :: Sn_at_intg(:),Sm_at_intg(:),SE_at_intg(:)
  
!  real(dp), allocatable :: SE_i_io(:)
!  real(dp), allocatable :: SE_i_cx(:)

!  real(dp) :: SE_i_io_intg,SE_i_cx_intg

  real(dp), allocatable :: Sn_mol(:,:)
  real(dp), allocatable :: Sm_mol(:,:)
  real(dp), allocatable :: SE_mol(:,:)

  real(dp), allocatable :: Sn_mol_intg(:),Sm_mol_intg(:),SE_mol_intg(:)

  real(dp), allocatable :: Sn_tion(:,:)
  real(dp), allocatable :: Sm_tion(:,:)
  real(dp), allocatable :: SE_tion(:,:)

  real(dp), allocatable :: Sn_tion_intg(:),Sm_tion_intg(:),SE_tion_intg(:)

  real(dp), allocatable :: Sn_pls(:,:)
  real(dp), allocatable :: Sm_pls(:,:)
  real(dp), allocatable :: SE_pls(:,:)

  real(dp), allocatable :: Sn_pls_intg(:),Sm_pls_intg(:),SE_pls_intg(:)

  real(dp), allocatable :: Spf(:)

  real(dp), allocatable :: Sn_tot(:,:)
  real(dp), allocatable :: Sm_tot(:,:)
  real(dp), allocatable :: SE_tot(:,:)

  real(dp) :: Nplasma_save,Eplasma_save

  integer :: icountTe,icountTi,icountDe


end module styx2eirene
