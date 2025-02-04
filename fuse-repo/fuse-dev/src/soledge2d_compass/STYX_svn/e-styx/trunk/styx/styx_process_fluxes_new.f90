subroutine styx_process_fluxes_new()
  use all_variables, only : global_parameters, interp_data2, reference_parameters, global_variables, globals
  use Mphysics
  use eirmod_precision
  use styx2eirene
  use eirmod_comprt, only : iunout
  use eirmod_csdvi
  implicit none

#ifdef S2D

  integer :: iion,isp
  character(2) :: ispc

  integer :: isurf,ipro,k,krec,i,istra
  integer :: IP1,IP2
  real(dp) :: total_flux_wall_ions,total_flux_wall_at,total_flux_core_at,total_Eflux_lost_wall_pl
  real(dp) :: total_flux_wall_mol,total_flux_core_mol,total_flux_core_tion,total_flux_wall_tion
  real(dp) :: total_flux_emitted_at,total_flux_emitted_mol,total_flux_emitted_tion
  real(dp) :: total_flux_emitted_atat,total_flux_emitted_molat, total_flux_emitted_molmol
  real(dp) :: total_Eflux_in,total_Eflux_emitted_at,total_Eflux_emitted_mol,total_Eflux_emitted_tion
  real(dp) :: total_Eflux_lost_cei,total_Eflux_lost_core_at,total_Eflux_lost_core_mol,total_Eflux_lost_core_tion
  real(dp) :: total_Eflux_lost_wall_at,total_Eflux_lost_wall_mol,total_Eflux_lost_wall_tion
  real(dp) :: total_Eflux_in_wall,total_Eflux_in_wall_at,total_Eflux_in_wall_mol,total_Eflux_in_wall_tion
  real(dp) :: total_Eflux_net_wall,total_Eflux_lost_wall
  real(dp) :: total_Eflux_styx,total_net_cx_power_flux,total_Eflux_styx_i

  real(dp) :: recyct,recycta,recyctm,pbr1,pbr2,pbr3
  real(dp) :: dneut_w,dneut_cei,Reff,Reffa
  real(dp), allocatable :: total_pflux_ions(:),total_pflux_at(:),total_pflux_mol(:)
  real(dp), allocatable :: total_eflux_ions(:),total_eflux_at(:),total_eflux_mol(:)

  real(dp), allocatable :: net_power_flux(:,:),net_cx_power_flux(:,:)
  real(dp), allocatable :: ref_power_flux(:,:),recomb_power_flux(:,:)
  real(dp), allocatable :: Vrad_power_flux(:,:,:),Vplasma_power_flux(:,:,:)
  real(dp), allocatable :: Vref_power_flux(:,:,:),Vrecomb_power_flux(:,:,:)
  real(dp), allocatable :: Vnet_cx_power_flux(:,:,:)
  real(dp), allocatable :: Vrad_at_power_flux(:,:,:) 
  real(dp), allocatable :: flux_in_eirene(:,:)
  real(dp) :: total_net_power_flux,total_ion_power_flux
  real(dp) :: warea,pi_,warea_tot
  real(dp) :: Tetri1,Tetri2,Titri1,Titri2,denstri,phi_i,phi_e,velparatri
  real(dp) :: phi_i_wall,peak_flux

  real(dp) :: dNplasma,dNplasma_dt,dEplasma,dEplasma_dt

  integer :: itri,iside,isou

  real(dp) :: ion_Eflux_wall,total_Eflux_se,ion_Eflux_se,Rsheath,REeff,Reat,Reio
  real(dp), allocatable :: flux_Estyx_e(:),flux_Estyx_i(:)
  real(dp) :: fsheath1,fsheath2,Mach

  real(dp), allocatable :: sigintv(:),sigints(:),cov(:),sig1(:),sig2(:)
  real(dp) :: sigrel_Sm,sigrel_SEe,sigrel_SEi

  real(dp) :: chk_sig1,chk_sig2,chk_sig3,chk_sig4,chk_sig5,chk_sig6,chk_sig7,chk_sig8,chk_sig9

  real(dp), parameter :: elcha = 1.6022e-19_dp

  !real(dp) :: Mach
  real(dp) :: Radflx_intg


  real(dp), allocatable :: Srad_Imp(:),Srad_at_tot(:)
  real(dp), allocatable :: Radflx2(:)
  integer*4 :: ntri, nspe, nspecies
  real(dp) :: molEnergy, ionEnergy
  integer*4 :: Zion, index0, Zspecies, iat, natom
  integer*4 :: itor,icell

  pi_=4._dp*atan(1._dp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   particle balance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! for the moment get the total fluxes out

  total_flux_wall_ions = eos(0)%potpl(1,2)+eos(0)%potpl(1,4)
  total_flux_wall_at   = eos(0)%potat(1,2)+eos(0)%potat(1,4)
  total_flux_wall_mol  = eos(0)%potml(1,2)+eos(0)%potml(1,4)
  total_flux_wall_tion = eos(0)%potio(1,2)+eos(0)%potio(1,4)
  total_flux_core_at   = eos(0)%potat(1,3)
  total_flux_core_mol  = eos(0)%potml(1,3)
  total_flux_core_tion = eos(0)%potio(1,3)

  if(abs(globals%flux_tot_out_ac(1)).ge.1.D-15) then
     write(*,*) total_flux_wall_ions
     write(iunout,*) 'total_flux_wall_ions/abs(globals%flux_tot_out_ac) = ',total_flux_wall_ions/globals%flux_tot_out_ac(1)
  else
     write(iunout,*) 'total_flux_wall_ions/abs(globals%flux_tot_out_ac) = infinity'
  end if
  if(abs(globals%source_ionz_tot(1)).ge.1.D-15) then
     write(iunout,*) 'Sn_intg/globals%source_ionz_tot = ',Sn_intg/globals%source_ionz_tot(1)
  else
     write(iunout,*) 'Sn_intg/globals%source_ionz_tot = infinity'
  end if

  ! fluxes emitted on the surfaces (bulk ions -> atoms, bulk ions -> molecules)

  total_flux_emitted_at = eos(0)%prfpat(1,2)+eos(0)%prfpat(1,4)
  total_flux_emitted_mol = eos(0)%prfpml(1,2)+eos(0)%prfpml(1,4)
  total_flux_emitted_tion= 0._dp !(to be defined) 

  ! fluxes emitted on the surfaces (atoms -> atoms)

  total_flux_emitted_atat = eos(0)%prfaat(1,2)+eos(0)%prfaat(1,4)

  ! fluxes emitted on the surfaces (atoms -> molecules)

  total_flux_emitted_molat = eos(0)%prfaml(1,2)+eos(0)%prfaml(1,4)

  ! fluxes emitted on the surfaces (molecules -> molecules)

  total_flux_emitted_molmol = eos(0)%prfmml(1,2)+eos(0)%prfmml(1,4)

  ! recycling coefficient for main ions

  if (total_flux_wall_ions > 0._dp) then
     recyct = (total_flux_emitted_at+2.d0*total_flux_emitted_mol)/total_flux_wall_ions
  else
     recyct = 0._dp
  endif
  ! recycling coefficient for atoms

  if (total_flux_wall_at > 0._dp) then
     recycta = (total_flux_emitted_atat+2.d0*total_flux_emitted_molat)/total_flux_wall_at
  else
     recycta = 0._dp
  endif
  ! recycling coefficient for molecules

  if (total_flux_wall_mol > 0._dp) then
     recyctm = total_flux_emitted_molmol/total_flux_wall_mol
  else
     recyctm = 0._dp
  endif

  ! quantity of neutrals absorbed on the wall (this is the pumped flux)

  dneut_w = (1.d0-recycta)*total_flux_wall_at+(1.d0-recyctm)*2.d0*total_flux_wall_mol

  ! quantity of neutrals absorbed on the cei

  dneut_cei = total_flux_core_at+2.d0*total_flux_core_mol

  ! effective recycling coefficient = fraction of the ion influx ionized in the domain

  if (recyct*total_flux_wall_ions > 0._dp) then
     Reff = recyct*(recyct*total_flux_wall_ions-dneut_w-dneut_cei)/(recyct*total_flux_wall_ions)
  else
     Reff = 0._dp
  endif

  dNplasma = globals%tot_N(1) -Nplasma_save
  dNplasma_dt = dNplasma/(float(nscycles+1)*global_variables%dt*reference_parameters%fields%tau0)

  ! read statistical noise for integral quantities (last call only)

  if (curiter == numiter) then

     !$     allocate(sigintv(100),sigints(59))
     !     allocate(cov(9),sig1(9),sig2(9))
     !
     !     sigintv=0
     !     sigints=0
     !
     !     do i=1,NSIGVI
     !     	sigintv(iih(i))=sgms(i)
     !     enddo
     !
     !     do i=1,NSIGSI
     !        sigints(iihw(i))=sgmws(i)
     !     enddo
     !
     !     do i=1,NSIGCI
     !        cov(i) =sgmcs(0,i)
     !        sig1(i)=sgmcs(1,i)
     !        sig2(i)=sgmcs(2,i)
     !     enddo
     !
     !     sig1(1:3)=sig1(1:3)*1e1_dp/elcha
     !     sig2(1:3)=sig2(1:3)*1e1_dp/elcha
     !     cov(1:3)=cov(1:3)*(1e1_dp/elcha)*(1e1_dp/elcha)
     !
     !    ! remove the mass in momentum source
     !
     !     sig1(1:3)=sig1(1:3)/(global_parameters%element_list(1)%mass*m_u)
     !     sig2(1:3)=sig2(1:3)/(global_parameters%element_list(1)%mass*m_u)
     !     cov(1:3)=cov(1:3)/((global_parameters%element_list(1)%mass*m_u)**2)
     !
     !     ! Energy sources, from W.cm-3 -> W.m-3
     !     sig1(4:9)=sig1(4:9)*1e6_dp
     !     sig2(4:9)=sig2(4:9)*1e6_dp
     !     cov(4:9)=cov(4:9)*(1e6_dp*1e6_dp)

     ! check results; below as in input file
     !	1    97     1    98
     !	1    97     1    99
     !	1    98     1    99
     !	1    33     1    39
     !	1    33     1    45
     !	1    39     1    45
     !	1    38     1    44
     !	1    38     1    50
     !	1    44     1    50

     !     chk_sig1 = (abs(SE_e_at_intg)/voltot_eirene)*sigintv(33)*1e-2_dp/sig1(4)
     !     chk_sig2 = (abs(SE_e_mol_intg)/voltot_eirene)*sigintv(39)*1e-2_dp/sig1(6)
     !     chk_sig3 = (abs(SE_e_tion_intg)/voltot_eirene)*sigintv(45)*1e-2_dp/sig2(5)

     !     chk_sig4 = (abs(SE_i_at_intg)/voltot_eirene)*sigintv(38)*1e-2_dp/sig1(7)
     !     chk_sig5 = (abs(SE_i_mol_intg)/voltot_eirene)*sigintv(44)*1e-2_dp/sig1(9)
     !     chk_sig6 = (abs(SE_i_tion_intg)/voltot_eirene)*sigintv(50)*1e-2_dp/sig2(8)

     !     chk_sig7 = (abs(Sm_at_intg)/voltot_eirene)*sigintv(97)*1e-2_dp/sig1(1)
     !     chk_sig8 = (abs(Sm_mol_intg)/voltot_eirene)*sigintv(98)*1e-2_dp/sig1(3)
     !     chk_sig9 = (abs(Sm_tion_intg)/voltot_eirene)*sigintv(99)*1e-2_dp/sig2(2)

     ! standard deviation = sqrt(variances) of the total sources
     !     sigrel_Sm  = sqrt(sig1(1)**2+sig2(2)**2+sig1(3)**2+2._dp*(cov(1)+cov(2)+cov(3)))
     !     sigrel_SEe = sqrt(sig1(4)**2+sig2(5)**2+sig1(6)**2+2._dp*(cov(4)+cov(5)+cov(6)))
     !     sigrel_SEi = sqrt(sig1(7)**2+sig2(8)**2+sig1(9)**2+2._dp*(cov(7)+cov(8)+cov(9)))

     ! relative noise level for total sources
     !     sigrel_Sm  = sigrel_Sm/abs(Sm_intg)*voltot_eirene*100._dp
     !     sigrel_SEe = sigrel_SEe/abs(SE_e_intg)*voltot_eirene*100._dp
     !     sigrel_SEi = sigrel_SEi/abs(SE_i_intg)*voltot_eirene*100._dp


  endif


  !  Reffa=Sn_intg/total_flux_wall_ions
  if(abs(globals%flux_tot_out_ac(1)).ge.1.D-15) then
     Reffa=globals%source_ionz_tot(1)/globals%flux_tot_out_ac(1)
  else
     Reffa=0.D0
  end if

  ! output

  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*                                 RECYCLING                                     *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'Total ion flux on the wall (part/s)                                = ',total_flux_wall_ions
  write(iunout,'(A70,ES12.4)') 'flux of atoms reemitted from recycling ions (at/s) =               = ',eos(0)%prfpat(1,2)
  write(iunout,'(A70,ES12.4)') 'mean energy of reemitted atoms (eV) =                              = '!,eos(0)%erfpat(1,2)/elcha/eos(0)%prfpat(1,2)
  write(iunout,'(A70,ES12.4)') 'flux of molecules reemitted from recycling ions (molecules/s)      = ',eos(0)%prfpml(1,2)
  write(iunout,'(A70,ES12.4)') 'mean energy of reemitted molecules (eV)                            = '!,eos(0)%erfpml(1,2)/elcha/eos(0)%prfpml(1,2)
  write(iunout,'(A70,ES12.4)') 'fraction of atoms in the recycling flux                            = '!,eos(0)%prfpat(1,2)/(eos(0)%prfpat(1,2)+eos(0)%prfpml(1,2))
  write(iunout,*)
  write(iunout,'(A70,ES12.4)') 'Recycling coefficient for main ions                                = ',recyct
  write(iunout,'(A70,ES12.4)') 'Effective Recycling coefficient for atoms                          = ',recycta
  write(iunout,'(A70,ES12.4)') 'Effective Recycling coefficient for molecules                      = ',recyctm
  write(iunout,'(A70,ES12.4)') 'Effective recycling coefficient (recalculated)                     = ',Reff
  write(iunout,'(A70,ES12.4)') 'Effective recycling coefficient (styx balance)                  = ',Reffa
  write(iunout,'(A70,ES12.4)') 'Expected flux amplification                                        = ',1._dp/(1._dp-Reffa)

  if (transfer_debug) then
     write(iunout,'(A82,ES12.4)') '****      WARNING= Exp. flux. Amp. inaccurate, transfer_debug=.true.       ****  '
  endif

  if(abs(globals%flux_tot_in_ac(1)).ge.1.D-15) then
     write(iunout,'(A70,ES12.4)') 'Current flux amplification                                         = ',abs(globals%flux_tot_out_ac(1)/globals%flux_tot_in_ac(1))
  else
     write(iunout,'(A70,ES12.4)') 'Current flux amplification                                         = infinity'
  end if
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*                                 PUMPING                                       *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'Pumped flux (recalculated) (H-D-T/s)                               = ',eos(0)%potat(1,4)+2._dp*eos(0)%potml(1,4)- & 
       (eos(0)%prfaat(1,4)+2._dp*eos(0)%prfaml(1,4)+2._dp*eos(0)%prfmat(1,4)+2._dp*eos(0)%prfmml(1,4))

  !write(iunout,*) ' Pumped flux in EIRENE wall, atoms (part/s) =', SPUMP_styx(1,2)
  !write(iunout,*) ' Pumped flux in EIRENE core, atoms (part/s) =', SPUMP_styx(1,3)
  write(iunout,'(A70,ES12.4)') 'Pumped flux in EIRENE pump, atoms (H-D-T/s)                        = ', eos(0)%spump(1,4)

  !write(iunout,*) ' Pumped flux in EIRENE wall, molecules (part/s) =', SPUMP_styx(2,2)
  !write(iunout,*) ' Pumped flux in EIRENE core, molecules (part/s) =', SPUMP_styx(2,3)
  write(iunout,'(A70,ES12.4)') 'Pumped flux in EIRENE pump, molecules (H-D-T/s)                    = ', 2._dp*eos(0)%spump(2,4)

  !write(iunout,*) ' Pumped flux in EIRENE wall, bulk ions (part/s) =', SPUMP_styx(4,2)
  !write(iunout,*) ' Pumped flux in EIRENE core, bulk ions (part/s) =', SPUMP_styx(4,3)
  write(iunout,'(A70,ES12.4)') 'Pumped flux in EIRENE pump, bulk ions (H+-D+-T+/s)                 = ', eos(0)%spump(4,4)

  write(iunout,'(A70,ES12.4)') 'Pump albedo (recalculated)                                         = '!,(eos(0)%prfaat(1,4)+2._dp*eos(0)%prfaml(1,4)+2._dp*eos(0)%prfmat(1,4)+2._dp*eos(0)%prfmml(1,4))/ &
  !(eos(0)%potat(1,4)+2._dp*eos(0)%potml(1,4))
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'


  write(iunout,*)
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*                           NEUTRAL PARTICLE FLUXES                             *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'Total cx flux (atoms) on the wall (part/s)                         = ',total_flux_wall_at
  write(iunout,'(A70,ES12.4)') 'Total cx flux (atoms) on the core interface (part/s)               = ',total_flux_core_at
  write(iunout,'(A70,ES12.4)') 'Total molecule flux on the wall (part/s)                           = ',total_flux_wall_mol
  write(iunout,'(A70,ES12.4)') 'Total molecule flux on the core interface (part/s)                 = ',total_flux_core_mol
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)


  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*                           SPUTTERED FLUXES                                    *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'Total flux sputtered by D+ ions         (part/s)                    = ',eos(0)%sptpl(1,2)
  write(iunout,'(A70,ES12.4)') 'Total flux sputtered by D atoms        (part/s)                     = ',eos(0)%sptat(1,2)
  write(iunout,'(A70,ES12.4)') 'Total flux sputtered by molecules  (part/s)                         = ',eos(0)%sptml(1,2)
  write(iunout,'(A70,ES12.4)') 'Total flux sputtered (sum over species) (part/s)                    = ',eos(0)%sptpl(1,2)


  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*                        GLOBAL PARTICLE BALANCE                                *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'Total ion/e source (part/s)                                        = ',globals%source_ionz_tot(1)
  if (curiter == numiter) then
     write(iunout,'(A70,ES12.4)') 'mean value of statistical noise, empirical estimation (%)          = ',sigintv(9)
  endif
  write(iunout,'(A70,ES12.4)') 'Total ion/e influx                                                 = ',abs(globals%flux_tot_in_ac(1))
  write(iunout,'(A70,ES12.4)') 'Total ion/e influx + ion/e source - Total ion/e outflow            = ',abs(globals%flux_tot_in_ac(1))+globals%source_ionz_tot(1)-abs(globals%flux_tot_out_ac(1))
  write(iunout,'(A70,ES12.4)') 'dNp/dt (part/s)                                                    = ',dNplasma_dt
  write(iunout,'(A70,ES12.4)') 'Total number of charged particles  Np                              = ',globals%tot_N(1)

  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)


  ! evaluate particle balance in EIRENE (should be satisfied at every time step !)
  ! pbr1 says: total numbers of D lost (ionsation, retention, pumps) is equal to what enters the system
  !            through recycling, recombination and gas puffs
  !            note: the total electron source Sn_intg includes electrons consummed by recombination

  ! the expressions below are not correct in time dependent mode: must include contribution of census
  ! i.e. particles stopped because time horizon reached, particles lauched from the census
  ! see in tmstep.f  : SPZ_FLX(ISTRA,ITYP,IATM) -> fluxes on census [this variable is not in a module]
  !                    FLUX(NSTRATA)            -> flux from census


  pbr1 = (sum(Sn_intg) + total_flux_core_at + 2.d0*total_flux_core_mol+ & 
       (1.d0-recycta)*total_flux_wall_at+(1.d0-recyctm)*2.d0*total_flux_wall_mol)/ & 
       (recyct*total_flux_wall_ions + total_puff )
  pbr2 = sum(Sn_intg) /(Reff*total_flux_wall_ions + total_puff )
  pbr3 = sum(Sn_intg) / (recyct*total_flux_wall_ions + total_puff  &
       -(total_flux_core_at + 2.d0*total_flux_core_mol+ & 
       (1.d0-recycta)*total_flux_wall_at+(1.d0-recyctm)*2.d0*total_flux_wall_mol))

  write(iunout,*) 'Eirene internal particle balance (1) (should be 1) = ',pbr1
  write(iunout,*) 'Eirene internal particle balance (2) (should be 1) = ',pbr2
  write(iunout,*) 'Eirene internal particle balance (3) (should be 1) = ',pbr3
  write(iunout,*)


  open(unit=555,file='eirene.particle_balance',position='append')
  write(555,'(4(ES14.7,1x))') recyct*100.d0,Reff*100.d0,1._dp/(1._dp-Reff)*abs(globals%flux_tot_in_ac(1)),globals%flux_tot_out_ac(1)
  close(555)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                    Momentum transfer					   !	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*            MOMENTUM TRANSFER BETWEEN SOLEDGE2D AND EIRENE                     *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'Total parallel momentum source (m/s2)                              = ',Sm_intg
  write(iunout,'(A70,ES12.4)') '  - atomic contribution (m/s2)                                     = ',Sm_at_intg
  write(iunout,'(A70,ES12.4)') '  - molecular contribution (m/s2)                                  = ',Sm_mol_intg
  write(iunout,'(A70,ES12.4)') '  - test ion contribution (m/s2)                                   = ',Sm_tion_intg
  !  if (curiter == numiter) then
  !     write(iunout,'(A70,ES12.4)') 'Mean value of statistical noise, empirical estimation (%)          = ',sigrel_Sm
  !     write(iunout,'(A70,ES12.4)') '  - for atomic contribution (%)                                    = ',sigintv(97)
  !     write(iunout,'(A70,ES12.4)') '  - for molecular contribution (%)                                 = ',sigintv(98)
  !     write(iunout,'(A70,ES12.4)') '  - for test ion contribution (%)                                  = ',sigintv(99)
  !  endif

  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                    Energy balance					   !	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! energy flux at sheath entrance

  total_Eflux_se=abs(globals%flux_totE_out_ac(1)+globals%flux_totE_out_ac(0))
  ion_Eflux_se=abs(globals%flux_totE_out_ac(1))

  ! energy flux in ion channel after sheath acceleration

  ion_Eflux_wall= eos(0)%eotpl(1,2)

  ! sheath amplification coefficient for ion energy flux

  if(abs(ion_Eflux_se).ge.1.D-15) then
     Rsheath=ion_Eflux_wall/ion_Eflux_se
  else
     Rsheath=1.D0
  end if

  ! energy in ( = energy from atoms/molecules recyling from bulk ions)

  total_Eflux_emitted_at   = eos(0)%erfpat(1,2)
  total_Eflux_emitted_mol  = eos(0)%erfpml(1,2)
  total_Eflux_emitted_tion = eos(0)%erfpio(1,2)

  total_Eflux_in = total_Eflux_emitted_at+total_Eflux_emitted_mol+total_Eflux_emitted_tion

  ! losses on the wall (energy influx - energy outflux)

  ! first on the cei (only outflux)

  total_Eflux_lost_core_at = eos(0)%eotat(1,3)
  total_Eflux_lost_core_mol = eos(0)%eotml(1,3)
  total_Eflux_lost_core_tion = eos(0)%eotio(1,3)

  total_Eflux_lost_cei = total_Eflux_lost_core_at + total_Eflux_lost_core_mol + total_Eflux_lost_core_tion

  ! then on the wall

  total_Eflux_lost_wall_pl = eos(0)%eotpl(1,2)
  total_Eflux_lost_wall_at = eos(0)%eotat(1,2)
  total_Eflux_lost_wall_mol = eos(0)%eotml(1,2)
  total_Eflux_lost_wall_tion = eos(0)%eotio(1,2)

  total_Eflux_lost_wall = total_Eflux_lost_wall_pl+total_Eflux_lost_wall_at + total_Eflux_lost_wall_mol + total_Eflux_lost_wall_tion

  total_Eflux_in_wall_at   = eos(0)%erfaat(1,2) + eos(0)%erfmat(1,2) + eos(0)%erfiat(1,2)
  total_Eflux_in_wall_mol  = eos(0)%erfaml(1,2) + eos(0)%erfmml(1,2) + eos(0)%erfiml(1,2)
  total_Eflux_in_wall_tion = eos(0)%erfaio(1,2) + eos(0)%erfmio(1,2) + eos(0)%erfiio(1,2)

  total_Eflux_in_wall = total_Eflux_in_wall_at + total_Eflux_in_wall_mol + total_Eflux_in_wall_tion

  total_Eflux_net_wall = total_Eflux_lost_wall - total_Eflux_in_wall

  ! effective energy recycling coefficient

  ! REeff=(total_Eflux_net_wall-total_Eflux_lost_cei)/total_Eflux_se
  if(abs(total_Eflux_se+eos(0)%eotat(1,2)+eos(0)%eotml(1,2)+eos(0)%eotio(1,2)).ge.1.D-15) then
     REeff=(total_Eflux_in_wall+total_Eflux_emitted_at)/(total_Eflux_se+eos(0)%eotat(1,2)+eos(0)%eotml(1,2)+eos(0)%eotio(1,2))
  else
     REeff=1.d0
  end if

  ! idem but for ions only
  if(abs(total_Eflux_se).ge.1.D-15) then
     REio = total_Eflux_emitted_at/total_Eflux_se
  else
     REio=1.D0
  end if
  ! cx only
  if(abs(total_Eflux_lost_wall_at).ge.1.d-15) then
     REat = total_Eflux_in_wall/total_Eflux_lost_wall_at
  else
     REat=1.d0
  end if


  ! try to asses energy conservation in EIRENE / to be reworked ...

  ! losses through interaction with electrons (radiation + ionization)
  ! affects only internal enegy when only atoms in the system !!!

  ! losses/gain through interaction with ions - SE_i_intg

  !ebr1 = (total_Eflux_in -SE_i_intg)/(abs(total_Eflux_lost_cei)+abs(total_Eflux_net_wall))

  ! total energy gain, algebraic (kinetic + internal)
  !Egain = total_Eflux_emitted_at+total_Eflux_emitted_mol+total_Eflux_emitted_tion + & 
  !     (13.6_dp*total_flux_emitted_at + 31.7_dp*total_flux_emitted_mol+16.3_dp*total_flux_emitted_tion)*elcha - &
  !     SE_i_intg !- (SE_e_intg-Srad_h_intg)

  ! total losses (at the wall + volumetric losses, by which internal energy is removed)
  ! only ionisation losses should be included, since radiation is only a transfer from electrons to photons
  !Eloss = abs(total_Eflux_lost_cei)+(total_flux_core_at*13.6_dp+total_flux_core_mol*31.7_dp+total_flux_core_tion*16.3_dp)*elcha + &
  !     abs(total_Eflux_net_wall)+((total_flux_wall_at-total_flux_emitted_at)*13.6_dp+ &
  !     (total_flux_wall_mol-total_flux_emitted_mol)*31.7_dp+ &
  !     (total_flux_wall_tion-total_flux_emitted_tion)*16.3_dp)*elcha + &
  !     Sn_intg*13.6_dp*elcha 

  !ebr1=Egain/Eloss
  !write(iunout,*) ' EIRENE internal energy balance (should be 1, but not yet correct) = ',ebr1


  ! calculate variation of energy between calls to styx_process_fluxes

  dEplasma = (globals%tot_E(0)+globals%tot_E(1))- Eplasma_save
  dEplasma_dt = dEplasma/(float(nscycles+1)*global_variables%dt*reference_parameters%fields%tau0)



!!! first check total fluxes on each kind of wall elements (recycling models, block 3A)

  if (.not.allocated(total_pflux_ions)) then
     allocate(total_pflux_ions(ipropmin:ipropmax))
     allocate(total_pflux_at(ipropmin:ipropmax))
     allocate(total_pflux_mol(ipropmin:ipropmax))
     allocate(total_eflux_ions(ipropmin:ipropmax))
     allocate(total_eflux_at(ipropmin:ipropmax))
     allocate(total_eflux_mol(ipropmin:ipropmax))
  endif

  total_pflux_ions=0._dp
  total_pflux_at=0._dp
  total_pflux_mol=0._dp
  total_eflux_ions=0._dp
  total_eflux_at=0._dp
  total_eflux_mol=0._dp

  do ipro=ipropmin,ipropmax
     do isurf=nsurf0,nsurf_tal
        if (surface(isurf)%iprop == ipro) then
           total_pflux_ions(ipro) = total_pflux_ions(ipro)+eos(0)%potpl(1,isurf)
           total_pflux_at(ipro) = total_pflux_at(ipro)+eos(0)%potat(1,isurf)
           total_pflux_mol(ipro) = total_pflux_mol(ipro)+eos(0)%potml(1,isurf)
           total_eflux_ions(ipro) = total_eflux_ions(ipro)+eos(0)%eotpl(1,isurf)
           total_eflux_at(ipro) = total_eflux_at(ipro)+eos(0)%eotat(1,isurf)
           total_eflux_mol(ipro) = total_eflux_mol(ipro)+eos(0)%eotml(1,isurf)
        endif
     enddo
  enddo


  ! calculate radiative power flux on wall (Watt) -- 2D only ---

  allocate(Radflx2(1:nsurf_tal))

  if (.not.is_3D) then

    Radflx=0._dp
    allocate(Srad_at_tot(1:Ntri_styx))
    Srad_at_tot=0._dp
    ! sum not running on possible extra species
    do isp=1,global_parameters%N_species
      Srad_at_tot(:) = Srad_at_tot(:) + Srad_at(:,isp)
    enddo

   ! this is now the total contribution from neutrals
    Radflx(1:nsurf_tal) = - matmul(Delta_Omega,Srad_at_tot*vol_tri_eirene)/(2._dp*pi_)

    allocate(Srad_Imp(1:NTRI_styx))
    Srad_Imp=0.D0
    do nspe=2,global_parameters%N_ions
      do ntri=1,NTRI_styx
        Srad_Imp(ntri)=Srad_Imp(ntri)+(interp_data2%knots_radiation(interp_data2%tri_knots(ntri,1),nspe)&
             +interp_data2%knots_radiation(interp_data2%tri_knots(ntri,2),nspe)&
             +interp_data2%knots_radiation(interp_data2%tri_knots(ntri,3),nspe))/3.D0
      end do
    end do
    Srad_Imp=Srad_Imp*reference_parameters%fields%n0*reference_parameters%fields%T0eV&
       *eV/reference_parameters%fields%tau0
    Radflx2(1:nsurf_tal) = - matmul(Delta_Omega,Srad_Imp*vol_tri_eirene)/(2._dp*pi_)

    ! total
    Radflx_intg = sum(Radflx(nsurf0:nsurf_tal)) +sum(Radflx2(1:nsurf_tal))

  else

    Radflx = 0.D0
    Radflx2 = 0.D0

  endif ! not is_3D

  ! calculate net power fluxes (Watt)

  if (.not.allocated(net_power_flux)) then
     allocate(net_power_flux(1:NSURF_TAL,1:Ntor_cells),net_cx_power_flux(1:NSURF_TAL,1:Ntor_cells))
     allocate(ref_power_flux(1:NSURF_TAL,1:Ntor_cells),recomb_power_flux(1:NSURF_TAL,1:Ntor_cells))
     allocate(Vrad_power_flux(1:NSURF_TAL,1:global_parameters%N_ions,1:Ntor_cells), &
		Vplasma_power_flux(1:NSURF_TAL,1:global_parameters%N_ions,1:Ntor_cells))
     allocate(Vref_power_flux(1:NSURF_TAL,1:global_parameters%N_species,1:Ntor_cells), &
		Vrecomb_power_flux(1:NSURF_TAL,1:global_parameters%N_ions,1:Ntor_cells))
     allocate(Vnet_cx_power_flux(1:NSURF_TAL,1:global_parameters%N_species,1:Ntor_cells), &
		Vrad_at_power_flux(1:NSURF_TAL,1:global_parameters%N_species,1:Ntor_cells))
  endif

  net_power_flux=0._dp
  net_cx_power_flux=0._dp

  ! incident power flux (ions+atoms+molecules) - (energy flux reemitted) + recombination energy for incident H ions 

  total_net_power_flux = 0._dp
  total_ion_power_flux = 0._dp
  total_Eflux_styx = 0._dp
  total_Eflux_styx_i = 0._dp
  total_net_cx_power_flux=0._dp

  do itor = 1, Ntor_cells
    do isurf=nsurf0,nsurf_tal
     ! exlude the cei from the calculation
       if (surface(isurf)%iprop /= 1) then
         itri=surface(isurf)%itri
         iside=surface(isurf)%iside
   
         icell = itri + (itor-1)*(Ntri_styx+1)

         net_cx_power_flux(isurf,itor) = eos(0)%eotat(1,isurf) + eos(0)%eotml(1,isurf) + eos(0)%eotio(1,isurf) &
                                ! energy reemitted by neutrals recycling from atoms, molecules, test ions
             - eos(0)%erfaat(1,isurf) - eos(0)%erfmat(1,isurf) - eos(0)%erfiat(1,isurf) &
             - eos(0)%erfaml(1,isurf) - eos(0)%erfmml(1,isurf) - eos(0)%erfiml(1,isurf) &
             - eos(0)%erfaio(1,isurf) - eos(0)%erfmio(1,isurf) - eos(0)%erfiio(1,isurf)

         ref_power_flux(isurf,itor) =  - eos(0)%erfpat(1,isurf)-eos(0)%erfpml(1,isurf)-eos(0)%erfpio(1,isurf)
         recomb_power_flux(isurf,itor) =  13.6_dp * eos(0)%potpl(1,isurf)*elcha &
             + 2.2_dp * (eos(0)%prfpml(1,isurf)+eos(0)%prfaml(1,isurf))*elcha 


        ! electron + ion power fluxes
         net_power_flux(isurf,itor)= Interp_Data2%tri_fluxE(itri,iside,0,itor)+Interp_Data2%tri_fluxE(itri,iside,1,itor) &
                                ! energy reemitted by neutrals recycling from bulk ions
             - eos(0)%erfpat(1,isurf)-eos(0)%erfpml(1,isurf)-eos(0)%erfpio(1,isurf) &             
                                ! recombination losses at wall (ions -> atoms)
             + 13.6_dp * eos(0)%potpl(1,isurf)*elcha &
                                ! recombination losses at wall (atoms -> molecules)
             + 2.2_dp * (eos(0)%prfpml(1,isurf)+eos(0)%prfaml(1,isurf))*elcha &
                                ! CX fluxes, including molecules and molecular ions
             + net_cx_power_flux(isurf,itor) &
                                ! radiative flux
             + Radflx(isurf) + Radflx2(isurf)


         total_Eflux_styx=total_Eflux_styx + & 
             Interp_Data2%tri_fluxE(itri,iside,0,itor)+Interp_Data2%tri_fluxE(itri,iside,1,itor)

         total_Eflux_styx_i = total_Eflux_styx_i + Interp_Data2%tri_fluxE(itri,iside,1,itor)

         total_net_power_flux = total_net_power_flux + net_power_flux(isurf,itor)
         total_net_cx_power_flux=total_net_cx_power_flux + net_cx_power_flux(isurf,itor)
         total_ion_power_flux = total_ion_power_flux + eos(0)%eotpl(1,isurf)
       endif
    enddo
  enddo ! itor

  ! total radiation added here for output
  !total_net_power_flux=total_net_power_flux + abs(Srad_h_intg)


  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*            ENERGY TRANSFER BETWEEN SOLEDGE2D AND EIRENE                       *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'total energy flux at sheath entrance from styx (MW)             = ',total_Eflux_se*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'idem, recalculated from InterpData2 (MW)                           = ',total_Eflux_styx*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'idem, in ion channel (MW)                                          = ',ion_Eflux_se*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'idem, recalculated from InterpData2 (MW)                           = ',total_Eflux_styx_i*1e-6_dp                       
  write(iunout,'(A70,ES12.4)') 'total power in ion channel on the wall, (reconstructed) (MW)       = ',total_ion_power_flux*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'total energy flux on wall in ion channel (EIRENE) (MW)             = ',ion_Eflux_wall*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'global sheath amplification coefficient for ion energy flux        = ',Rsheath
  write(iunout,'(A70,ES12.4)') 'Effective energy recycling coefficient RE                          = ',REeff
  write(iunout,'(A70,ES12.4)') '  -for ions only                                                   = ',REio
  write(iunout,'(A70,ES12.4)') '  -for atoms only                                                  = ',REat
  write(iunout,'(A70,ES12.4)') 'Kinetic Energy influx (neutrals) in EIRENE (MW)                    = ',total_Eflux_in*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'Total energy source for electrons (MW)                             = ',SE_intg(0)*1e-6_dp
  if (curiter == numiter) then
     write(iunout,'(A70,ES12.4)') '  - mean value of statistical noise, empirical estimation (%)      = ',(sigintv(33)+sigintv(39)+sigintv(45))/3._dp
  endif
  write(iunout,'(A70,ES12.4)') 'Total energy source for ions (MW)                                  = ',(sum(SE_intg)-SE_intg(0))*1e-6_dp
  !write(iunout,'(A70,ES12.4)') '  - contribution of ionization (MW)                                = ',SE_i_io_intg*1e-6_dp
  !write(iunout,'(A70,ES12.4)') '  - contribution of charge exchange (MW)                           = ',SE_i_cx_intg*1e-6_dp
  !if (curiter == numiter) then
  !   write(iunout,'(A70,ES12.4)') '  - mean value of statistical noise, empirical estimation (%)      = ',(sigintv(38)+sigintv(44)+sigintv(50))/3._dp
  !endif

  write(iunout,'(A70,ES12.4)') 'Total net cx power flux  (MW)                                      = ',total_net_cx_power_flux*1e-6_dp
  ! has to be integrated (here size Ntri)
  !write(iunout,'(A70,ES12.4)') 'total radiated power, volume integral (MW)                         = ',abs(Srad_at_tot+Srad_Imp)*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'total radiated power, received by the wall (MW)                    = ',abs(Radflx_intg)*1e-6_dp
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A82,ES12.4)') '*                            GLOBAL ENERGY BALANCE                              *'
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'

  write(iunout,'(A70,ES12.4)') 'Pin (MW)                                                           = ',abs(globals%flux_totE_in_ac(1)+globals%flux_totE_in_ac(0))*1e-6_dp

  write(iunout,'(A70,ES12.4)') 'Pnet = total net power flux on the wall                            = ',total_net_power_flux*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'Pcei = total power flux on the cei (MW)                            = ',total_Eflux_lost_cei*1e-6_dp

  write(iunout,'(A70,ES12.4)') 'Pout = Pnet+Pcei      (MW)                                         = ',(total_net_power_flux+total_Eflux_lost_cei)*1e-6_dp

  write(iunout,'(A70,ES12.4)') 'Pin-Pout (MW)                                                      = ',abs(globals%flux_totE_in_ac(1)+globals%flux_totE_in_ac(0))*1e-6_dp- &
       (total_net_power_flux+total_Eflux_lost_cei)*1e-6_dp    
  write(iunout,'(A70,ES12.4)') 'dEp/dt   (MW)                                                      = ',dEplasma_dt*elcha*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'Total energy content in plasma : Ep   (J)                          = ',(globals%tot_E(0)+globals%tot_E(1))*elcha  
  write(iunout,'(A70,ES12.4)') 'total recombination energy contribution (ions -> atoms) (MW)       = ',13.6_dp * eos(0)%potpl(1,2)*elcha*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'total recombination energy contribution (atoms-> mol. ) (MW)       = ',2.2_dp * (eos(0)%prfpml(1,2)+eos(0)%prfaml(1,2))*elcha*1e-6_dp


  ! chi=2.2 eV (Stangeby) for D+D->D2 is consistent with potential curves, which gives ~4.5 eV to be divided by 2 ...

!!! get the fluxes used as input by EIRENE along wall contour
!!! eos(0)%potpl is sampled according to this distriburtion, and should tend to it as N-part _> infinity

  allocate(flux_in_eirene(nsurf0:Nsurf_tal+1,1:Ntor_cells))

  if (.not. allocated(flux_Estyx_e)) then
     allocate(flux_Estyx_e(nsurf0:nsurf_tal),flux_Estyx_i(nsurf0:nsurf_tal))
  endif

  ! atom energy flux
  do iat=1,global_parameters%n_species
     Radflx2(1:nsurf_tal) = - matmul(Delta_Omega,Srad_at(:,iat)*vol_tri_eirene)/(2._dp*pi_)
     do isurf=nsurf0,nsurf_tal
        ! exlude the cei from the calculation
        if (surface(isurf)%iprop /= 1) then
           if(iat.eq.1) then !special treatment for H
              Vnet_cx_power_flux(isurf,iat,1) = eos(0)%eotat(iat,isurf) + eos(0)%eotml(1,isurf) + eos(0)%eotio(1,isurf) &
                   ! energy reemitted by neutrals recycling from atoms, molecules, test ions
                   - eos(0)%erfaat(iat,isurf) - eos(0)%erfmat(iat,isurf) - eos(0)%erfiat(iat,isurf) &
                   - eos(0)%erfaml(1,isurf) - eos(0)%erfmml(1,isurf) - eos(0)%erfiml(1,isurf) &
                   - eos(0)%erfaio(1,isurf) - eos(0)%erfmio(1,isurf) - eos(0)%erfiio(1,isurf)
              Vref_power_flux(isurf,iat,1) = - eos(0)%erfpat(iat,isurf)-eos(0)%erfpml(1,isurf)-eos(0)%erfpio(1,isurf) 
           else
              Vnet_cx_power_flux(isurf,iat,1) = eos(0)%eotat(iat,isurf) &
                   ! energy reemitted by neutrals recycling from atoms, molecules, test ions
                   - eos(0)%erfaat(iat,isurf) 
              Vref_power_flux(isurf,iat,1) = - eos(0)%erfpat(iat,isurf)
           end if
           Vrad_at_power_flux(isurf,iat,1)=Radflx2(isurf)
        end if
     end do
  end do


  ! ions energy flux
  do iion=1,global_parameters%n_ions
     do ntri=1,NTRI_styx
        Srad_Imp(ntri)=-(interp_data2%knots_radiation(interp_data2%tri_knots(ntri,1),iion)&
             +interp_data2%knots_radiation(interp_data2%tri_knots(ntri,2),iion)&
             +interp_data2%knots_radiation(interp_data2%tri_knots(ntri,3),iion))/3.D0
     end do
     Srad_Imp=Srad_Imp*reference_parameters%fields%n0*reference_parameters%fields%T0eV&
          *eV/reference_parameters%fields%tau0
     Radflx2(1:nsurf_tal) = - matmul(Delta_Omega,Srad_Imp*vol_tri_eirene)/(2._dp*pi_)
     do isurf=nsurf0,nsurf_tal
        ! exlude the cei from the calculation
        if (surface(isurf)%iprop /= 1) then
           itri=surface(isurf)%itri
           iside=surface(isurf)%iside

           nspecies=global_parameters%ions_list(iion,1)
           Zspecies=global_parameters%ions_list(iion,2)
           ionEnergy = 0.
           do Zion=1,Zspecies
              ionEnergy=ionEnergy+global_parameters%element_list(nspecies)%amdatas(Zion)%ionization_potential
           end do
           if (iion.eq.1) then ! Hydrogen ##### careful to be improved 
              molEnergy = 2.2
           else
              molEnergy = 0.
           end if

           Vplasma_power_flux(isurf,iion,1)=Interp_Data2%tri_fluxE(itri,iside,iion,1)
           Vrecomb_power_flux(isurf,iion,1)=ionEnergy * eos(0)%potpl(iion,isurf)*elcha &
                +molEnergy * (eos(0)%prfpml(iion,isurf)+eos(0)%prfaml(iion,isurf))*elcha
!           if(iion.eq.1) then !special treatment for hydrogen
!              Vref_power_flux(isurf,iion)=-eos(0)%erfpat(iion,isurf)-eos(0)%erfpml(iion,isurf)-eos(0)%erfpio(iion,isurf) 
!           else
!              Vref_power_flux(isurf,iion)=-eos(0)%erfpat(iion,isurf)
!           end if
           Vrad_power_flux(isurf,iion,1)=Radflx2(isurf)
        end if
     enddo
  enddo

!!! convert fluxes to fluxes densities on the fly (divide by surface area, m2)

  if (.not. allocated(flux_Estyx_e)) then
     allocate(flux_Estyx_e(nsurf0:nsurf_tal),flux_Estyx_i(nsurf0:nsurf_tal))
  endif
  warea_tot=0._dp

  do isurf=nsurf0,nsurf_tal
     if (surface(isurf)%iprop /= 1) then
        itri = surface(isurf)%itri
        iside = surface(isurf)%iside
        ! ds in m, R1/R2 in cm...
        warea = pi_*(surface(isurf)%R1+surface(isurf)%R2)*surface(isurf)%ds*1e-2_dp
        if (warea == 0._dp) then
           write(*,*) ' zero surface area calculated for wall surface element isurf = ',isurf
           call eirene_exit_own(1)
        endif
        flux_in_eirene(isurf,1) = pflux_in(itri,iside,1)/warea
        !particles
        eos(0)%potpl(:,isurf)=eos(0)%potpl(:,isurf)/warea
        eos(0)%potat(:,isurf)=eos(0)%potat(:,isurf)/warea
        eos(0)%potml(:,isurf)=eos(0)%potml(:,isurf)/warea
        !reflections particles
        eos(0)%prfpat(:,isurf)=eos(0)%prfpat(:,isurf)/warea
        eos(0)%prfpml(:,isurf)=eos(0)%prfpml(:,isurf)/warea
        eos(0)%prfaat(:,isurf)=eos(0)%prfaat(:,isurf)/warea
        eos(0)%prfaml(:,isurf)=eos(0)%prfaml(:,isurf)/warea
        !sputtering
        eos(0)%sptpl(:,isurf)=eos(0)%sptpl(:,isurf)/warea
        eos(0)%sptat(:,isurf)=eos(0)%sptat(:,isurf)/warea
        !energy
        eos(0)%eotpl(:,isurf)=eos(0)%eotpl(:,isurf)/warea
        eos(0)%eotat(:,isurf)=eos(0)%eotat(:,isurf)/warea
        eos(0)%eotml(:,isurf)=eos(0)%eotml(:,isurf)/warea
        flux_Estyx_e(isurf)=Interp_Data2%tri_fluxE(itri,iside,0,1)/warea
        flux_Estyx_i(isurf)=Interp_Data2%tri_fluxE(itri,iside,1,1)/warea
        net_power_flux(isurf,1)=net_power_flux(isurf,1)/warea
        net_cx_power_flux(isurf,1)=net_cx_power_flux(isurf,1)/warea
        ref_power_flux(isurf,1)=ref_power_flux(isurf,1)/warea
        recomb_power_flux(isurf,1)=recomb_power_flux(isurf,1)/warea
        Radflx(isurf)=Radflx(isurf)/warea
        Radflx2(isurf)=Radflx2(isurf)/warea
        warea_tot=warea_tot+warea
        
     endif
  enddo

  do isurf=nsurf0,nsurf_tal
     if (surface(isurf)%iprop /= 1) then
        itri = surface(isurf)%itri
        iside = surface(isurf)%iside
	warea = pi_*(surface(isurf)%R1+surface(isurf)%R2)*surface(isurf)%ds*1e-2_dp
        if (warea == 0._dp) then
           write(*,*) ' zero surface area calculated for wall surface element isurf = ',isurf
           call eirene_exit_own(1)
        endif
        do iion=1,global_parameters%n_ions
           Vplasma_power_flux(isurf,iion,1)=Vplasma_power_flux(isurf,iion,1)/warea
           Vrecomb_power_flux(isurf,iion,1)=Vrecomb_power_flux(isurf,iion,1)/warea
           Vrad_power_flux(isurf,iion,1)=Vrad_power_flux(isurf,iion,1)/warea
        end do
        do iat=1,global_parameters%n_species
           Vrad_at_power_flux(isurf,iat,1)=Vrad_at_power_flux(isurf,iat,1)/warea
           Vnet_cx_power_flux(isurf,iat,1)=Vnet_cx_power_flux(isurf,iat,1)/warea
           Vref_power_flux(isurf,iat,1)=Vref_power_flux(isurf,iat,1)/warea
        end do
     end if
  end do
  !separate loop for Radflx, already in right order
  ! do i=1,Nsou
  !       isurf=i+nsurf0
  !       k=ksurf(isurf)
  !       warea=pi_*(surface(k)%R1+surface(k)%R2)*surface(k)%ds*1e-2_dp
  ! 	Radflx(k+nsurf0)=Radflx(k+nsurf0)/warea
  ! enddo

  peak_flux = maxval(net_power_flux)

  write(iunout,'(A70,ES12.4)') 'Peak heat flux (MW/m2)                                             = ',peak_flux*1e-6_dp
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)



!!! output to assess parralelization efficiency

  if (nprs_styx > 1) then

     write(iunout,'(A82,ES12.4)') '*********************************************************************************'
     write(iunout,'(A82,ES12.4)') '*            	       EIRENE PARALLELISATION                                 *'
     write(iunout,'(A82,ES12.4)') '*********************************************************************************'

     write(iunout,'(A70,i6)') 'number of cores                                                    = ',nprs_styx
     write(iunout,'(A70)')    'number of histories calculated:                                      '
     do istra=1,nstrata
        write(iunout,'(A14,i2,A54,i8)') '    - strata #',istra,'                                                        ',xmcp_styx(istra)
     enddo
     write(iunout,'(A70)')    'Shortest computational time:                                         '
     write(iunout,'(A15,i5,A49,es12.4)') 'process rank = ',rtimin,' ,              time(s)                         = ',ptimin
     write(iunout,'(A70)')    'Longest computational time:                                         '
     write(iunout,'(A15,i5,A49,es12.4)') 'process rank = ',rtimax,' ,              time(s)                         = ',ptimax
     write(iunout,'(A70,es12.4)') 'Mean computationnal time  (s)                                      = ',ptmean

     write(iunout,'(A82,ES12.4)') '*********************************************************************************'
     write(iunout,*)
  endif


!!! get fluxes with number and sides of triangles


  !  open(unit = 555,file =  trim(adjustl(fluid_code))//'.particle_fluxes',status = 'replace')

  ! output in part/m2/s
  do isp=1,global_parameters%n_species
     write(ispc,'(i2)') isp
     open(unit=555,file = trim(adjustl(fluid_code))//'.atoms_fluxes_wall_'//trim(adjustl(ispc)),status = 'replace')
     do isou=1,Nsou
        isurf=isou+nsurf0
        k=ksurf(isurf)
        itri=surface(k)%itri
        iside=surface(k)%iside
        write(555,'(3i6,1x,8(es14.7,1x))') itri,iside,surface(k)%iprop, & 
             ssurf(isurf),eos(0)%potat(isp,k),eos(0)%potml(1,k),eos(0)%prfpat(isp,k),eos(0)%prfpml(1,k),eos(0)%prfaat(isp,k),eos(0)%prfaml(1,k)
     enddo
     close(555)
  enddo


  do iion=1,global_parameters%n_ions
     write(ispc,'(i2)') iion
     open(unit=555,file =  trim(adjustl(fluid_code))//'.ion_fluxes_wall_'//trim(adjustl(ispc)),status = 'replace')
     do isou=1,Nsou
        isurf=isou+nsurf0
        k=ksurf(isurf)
        itri=surface(k)%itri
        iside=surface(k)%iside
        natom=global_parameters%ions_list(iion,1)
        write(555,'(3i6,1x,6(es14.7,1x))') itri,iside,surface(k)%iprop, & 
             ssurf(isurf),eos(0)%potpl(iion,k),eos(0)%prfpat(natom,k) !,eos(0)%potml(isp,k)
     enddo
     close(555)
  enddo



!!!!! sanity check !!!!!!!

  !  open(unit = 555,file = 'styx2D.particle_fluxes',status = 'replace')

  ! output in part/m2/s
  !  open(unit=555,file = 'styx2D.particle_fluxes_wall_recsurf',status = 'replace')
  !  do isou=1,Nsou
  !     krec=krecsurf(isou)
  !     itri=recsurf(krec)%itri
  !     iside=recsurf(krec)%iside

  !     write(555,'(3i6,1x,6(es14.7,1x))') itri,iside,recsurf(krec)%iprop, & 
  !          ssurf(isurf),flux_in_eirene(kr2e(krec)),Interp_Data2%tri_fluxN(itri,iside,1),eos(0)%potpl(1,kr2e(krec)), &
  !          eos(0)%potat(1,kr2e(krec)),eos(0)%potml(1,kr2e(krec))
  !  enddo

  !  close(555)


  ! output in part/m2/s
  do iion=1,global_parameters%n_ions
     write(ispc,'(i2)') iion
     open(unit=555,file =  trim(adjustl(fluid_code))//'.fluxes_sputtered_by_ions_'//trim(adjustl(ispc)),status = 'replace')
     do isou=1,Nsou
        isurf=isou+nsurf0
        k=ksurf(isurf)
        itri=surface(k)%itri
        iside=surface(k)%iside
        write(555,'(3i6,1x,6(es14.7,1x))') itri,iside,surface(k)%iprop, & 
             ssurf(isurf),eos(0)%sptpl(iion,k)
     enddo
     close(555)
  enddo

  do isp=1,global_parameters%n_species
     write(ispc,'(i2)') isp
     open(unit=555,file =  trim(adjustl(fluid_code))//'.fluxes_sputtered_by_atoms_'//trim(adjustl(ispc)),status = 'replace')
     do isou=1,Nsou
        isurf=isou+nsurf0
        k=ksurf(isurf)
        itri=surface(k)%itri
        iside=surface(k)%iside
        write(555,'(3i6,1x,6(es14.7,1x))') itri,iside,surface(k)%iprop, &
             ssurf(isurf),eos(0)%sptat(isp,k)
     enddo
     close(555)
  enddo


  open(unit = 555,file =  trim(adjustl(fluid_code))//'.energy_fluxes_1',status = 'replace')
  open(unit = 556,file =  trim(adjustl(fluid_code))//'.energy_fluxes_details_1',status = 'replace')
  write(556,*) '%s fe fi fnet fcx frad fref frec'

  ! output in units of W
  do isou=1,Nsou
     isurf=isou+nsurf0
     k=ksurf(isurf)
     itri=surface(k)%itri
     iside=surface(k)%iside
     warea = pi_*(surface(k)%R1+surface(k)%R2)*surface(k)%ds*1e-2_dp
     write(555,'(3i6,1x,9(es14.7,1x))') itri,iside,surface(k)%iprop, & 
          ssurf(isurf),eos(0)%eotpl(1,k),flux_Estyx_e(k), & 
          flux_Estyx_i(k),net_power_flux(k,1),net_cx_power_flux(k,1),Radflx(k)+Radflx2(k),warea

     write(556,'(8(es14.7,1x))') ssurf(isurf),flux_Estyx_e(k), & 
          flux_Estyx_i(k),net_power_flux(k,1),net_cx_power_flux(k,1),Radflx(k),&
          ref_power_flux(k,1),recomb_power_flux(k,1)
  enddo

  close(555)
  close(556)


  !atom and molecules energy flux
  do isp=1,global_parameters%n_species
     write(ispc,'(i2)') isp
     open(unit=555,file =  trim(adjustl(fluid_code))//'.atoms_Efluxes_wall_'//trim(adjustl(ispc)),status = 'replace')
     do isou=1,Nsou
        isurf=isou+nsurf0
        k=ksurf(isurf)
        itri=surface(k)%itri
        iside=surface(k)%iside
        write(555,'(3i6,1x,6(es14.7,1x))') itri,iside,surface(k)%iprop, & 
             ssurf(isurf),Vnet_cx_power_flux(k,isp,1),Vrad_at_power_flux(k,isp,1),Vref_power_flux(k,isp,1)
     enddo
     close(555)
  enddo



  do iion=1,global_parameters%n_ions
     write(ispc,'(i2)') iion
     open(unit=555,file =  trim(adjustl(fluid_code))//'.ion_Efluxes_wall_'//trim(adjustl(ispc)),status = 'replace')
     do isou=1,Nsou
        isurf=isou+nsurf0
        k=ksurf(isurf)
        itri=surface(k)%itri
        iside=surface(k)%iside
        Write(555,'(3i6,1x,6(es14.7,1x))') itri,iside,surface(k)%iprop, ssurf(isurf),&
             Vplasma_power_flux(k,iion,1), & ! plasma flux
             Vrecomb_power_flux(k,iion,1), & ! recombination flux
!	     Vref_power_flux(k,iion), &    ! reflected flux
             Vrad_power_flux(k,iion,1)       ! radiated flux
     end do
     close(555)
  end do

  ! electrons energy flux
  write(ispc,'(i2)') 0
  open(unit=555,file =  trim(adjustl(fluid_code))//'.ion_Efluxes_wall_'//trim(adjustl(ispc)),status = 'replace')
  do isou=1,Nsou
     isurf=isou+nsurf0
     k=ksurf(isurf)
     itri=surface(k)%itri
     iside=surface(k)%iside
     Write(555,'(3i6,1x,6(es14.7,1x))') itri,iside,surface(k)%iprop, & 
          ssurf(isurf),flux_Estyx_e(k) ! plasma flux

  enddo
  close(555)



  ! output of plasma parameters on the wall surfaces
  ! temperatures in eV, and densities in m-3

  open(unit=555, file=  trim(adjustl(fluid_code))//'.plasma_param_surfaces',status= 'replace')

  phi_e=0._dp
  phi_i=0._dp
  phi_i_wall=0._dp

  do isou=1,Nsou
     isurf=isou+nsurf0
     k=ksurf(isurf)
     itri=surface(k)%itri
     iside=surface(k)%iside

     if (ISIDE == 1) then
        IP1=NVERT(1,ITRI)
        IP2=NVERT(2,ITRI)
     elseif (ISIDE == 2) then
        IP1=NVERT(2,ITRI)
        IP2=NVERT(3,ITRI)
     else
        IP1=NVERT(3,ITRI)
        IP2=NVERT(1,ITRI)
     endif

     denstri=0.5_dp*(Interp_data2%knots_density(IP1,0,1)+Interp_data2%knots_density(IP2,0,1))
     denstri=denstri*reference_parameters%fields%n0


     ! interpolated data
     Tetri1=0.5_dp*(Interp_data2%knots_temperature(IP1,0,1)+Interp_data2%knots_temperature(IP2,0,1))  			
     Titri1=0.5_dp*(Interp_data2%knots_temperature(IP1,1,1)+Interp_data2%knots_temperature(IP2,1,1))
     velparatri=0.5_dp*(Interp_data2%knots_velocity(IP1,1,1)+Interp_data2%knots_velocity(IP2,1,1))
     TeTri1=TeTri1*reference_parameters%fields%T0*kB/elcha
     TiTri1=TiTri1*reference_parameters%fields%T0*kB/elcha

     if (Interp_data2%tri_fluxN(ITRI,ISIDE,1,1)*Interp_data2%tri_fluxE(ITRI,ISIDE,0,1)*Interp_data2%tri_fluxE(ITRI,ISIDE,1,1) /= 0.d0) then
        Tetri2=Interp_data2%tri_fluxE(ITRI,ISIDE,0,1)/(4.5_dp*Interp_Data2%tri_fluxN(ITRI,ISIDE,1,1))
        Mach=1._dp   !0.5_dp*(Interpdata2_%knots_M(IP1)+Interpdata2_%knots_M(IP2))
        Titri2=(Interp_data2%tri_fluxE(ITRI,ISIDE,1,1)-0.5_dp*Mach*Interp_data2%tri_fluxN(ITRI,ISIDE,1,1)*Tetri2)/ &
             ((2.5_dp+0.5_dp*Mach)*Interp_Data2%tri_fluxN(ITRI,ISIDE,1,1))
     else
        Tetri2=-1._dp*elcha
        Titri2=-1._dp*elcha
     endif

     ! same threshold as in styx_vertices ... minT, not necessarily relevant
     Tetri2=max(Tetri2/elcha,3d-1)
     Titri2=max(Titri2/elcha,3d-1)

     ! reconstruct total fluxes at sheath entrance
     phi_e=phi_e+Interp_data2%tri_fluxE(ITRI,ISIDE,0,1)
     phi_i=phi_i+Interp_data2%tri_fluxE(ITRI,ISIDE,1,1)

     ! sheath acceleration factor (energy gain for ions : -fsheath*Te)
     ! as in EIRENE (assumes ion distribution is a drifting maxwellian, cs perpendicular to the wall)

     ! has to be calculated with the correct multifluid extention !!
     fsheath1=0.5_dp !*log(2.*pi_*9.1094e-31_dp/(global_parameters%element_list(1)%mass*m_u)*(1._dp+Titri1/Tetri1))
     fsheath2=0.5_dp !*log(2.*pi_*9.1094e-31_dp/(global_parameters%element_list(1)%mass*m_u)*(1._dp+Titri2/Tetri2))
     ! reconstruct total ion flux on the wall after sheath acceleration

     phi_i_wall=phi_i_wall+Interp_data2%tri_fluxE(ITRI,ISIDE,1,1)&
          +Interp_Data2%tri_fluxN(ITRI,ISIDE,1,1)*(-1._dp*fsheath2)*Tetri2*elcha


     write(555,'(3i6,1x,9(es14.7,1x))') itri,iside,surface(k)%iprop, & 
          ssurf(isurf),Tetri1,Titri1,Tetri2,Titri2,denstri,fsheath1,fsheath2,velparatri

  enddo

  close(555)

  ! check if what happens in the sheath is physical

  if (ion_Eflux_wall > total_Eflux_se) then
     write(*,*) '---------------------------------------------------------'
     write(*,*) ' Energy is created in the sheath !!!!'
     write(*,*) ' total flux at sheath entrance (MW) = ',total_Eflux_se
     write(*,*) ' ion flux on the wall, sampled (MW) = ',ion_Eflux_wall*1e-6_dp
     write(*,*) ' ion flux on the wall, fsheath (MW) = ',phi_i_wall*1e-6_dp
     write(*,*) '---------------------------------------------------------'
  endif

  write(iunout,'(A82,ES12.4)') '*********************************************************************************'  
  write(iunout,'(A82,ES12.4)') '*               CHECKS ON SHEATH IMPLEMENTATION                                 *' 
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,'(A70,ES12.4)') 'total electron flux from styx2D             (MW, se)            = ',phi_e*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'total ion flux flux from styx2D             (MW, se)            = ',phi_i*1e-6_dp
  write(iunout,'(A70,ES12.4)') 'total ion flux recalculated (MW, at wall)                          = ',phi_i_wall*1e-6_dp
  write(iunout,'(A82,ES12.4)') '*********************************************************************************'
  write(iunout,*)

  deallocate(flux_in_eirene)
  if (allocated(sigintv)) deallocate(sigintv,sigints)

  Nplasma_save=globals%tot_N(1)
  Eplasma_save=globals%tot_E(1)+globals%tot_E(0)

!!!!!!!!!!! checks on incidence angles !!!!!!!!!!!!!!!

  write(*,*) ' Reporting on incident angles ...  '

  ! calculating average angles (degrees)
  where (sheath1D(:)%hit_V /=0)
     sheath1D(:)%alpha_V=sheath1D(:)%alpha_V/float(sheath1D(:)%hit_V)*180._dp/pi_
     sheath1D(:)%beta_V=sheath1D(:)%beta_V/float(sheath1D(:)%hit_V)*180._dp/pi_
  elsewhere
     sheath1D(:)%alpha_V=0._dp
     sheath1D(:)%beta_V=0._dp
  end where

  open(unit=555, file=  trim(adjustl(fluid_code))//'.incidence_angles',status= 'replace')


  do isou=1,Nsou
     isurf=isou+nsurf0
     k=ksurf(isurf)
     krec = krecsurf(isou)
     write(555,'(6ES12.4,i6)') ssurf(k),sheath1D(krec)%alphaB,sheath1D(krec)%tau,sheath1D(krec)%ksi, &
          sheath1D(krec)%alpha_V,sheath1D(krec)%beta_V,sheath1D(krec)%hit_V
  enddo

  close(555)

#endif

end subroutine styx_process_fluxes_new
