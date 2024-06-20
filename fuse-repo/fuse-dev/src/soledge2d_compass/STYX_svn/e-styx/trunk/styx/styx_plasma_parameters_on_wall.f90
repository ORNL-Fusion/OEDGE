subroutine styx_plasma_parameters_on_wall(ip1,ip2,itri,iside,isurf,itor,ipls,isp,dens,Te,Ti,csx,csy,csz,pflux)
  use all_variables, only : global_parameters, interp_data2, reference_parameters
  use Mphysics
  use styx2eirene
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_ccona
  use eirmod_ctrig
  implicit none
  integer, intent(in) :: ip1,ip2,itri,iside,ipls,isurf,itor,isp
  real(dp), intent(out) :: dens,Te,Ti,csx,csy,csz,pflux
  real(dp) :: de,Tec,Tec2,Tic,Tic2,cstri,Mach
  logical :: neg_Te,neg_Ti,neg_de

  neg_Te=.false.
  neg_Ti=.false.
  neg_de=.false.

  if (Eflx_calc_pen == 0) then

   Tec=0.5_dp*(Interp_data2%knots_temperature(IP1,0,itor)+Interp_data2%knots_temperature(IP2,0,itor))
   Tic=0.5_dp*(Interp_data2%knots_temperature(IP1,ipls,itor)+Interp_data2%knots_temperature(IP2,ipls,itor))

   Tec=Tec*reference_parameters%fields%T0*kB
   Tic=Tic*reference_parameters%fields%T0*kB

   cstri=sqrt((Tec+Tic)/(global_parameters%element_list(isp)%mass*m_u))*100._dp

   Tec=Tec/elcha
   Tic=Tic/elcha

   if (Tec<0._dp) Tec=1e-2_dp
   if (Tic<0._dp) Tic=1e-2_dp

   ! use cs as velocity perpendicular to the wall (Bohm-Chodura)
   csx=cstri*PTRIX(ISIDE,ITRI)
   csy=cstri*PTRIY(ISIDE,ITRI)
   csz=0._dp

  else 
    if (Interp_Data2%tri_fluxN(ITRI,ISIDE,ipls,itor)*Interp_Data2%tri_fluxE(ITRI,ISIDE,0,itor)*Interp_Data2%tri_fluxE(ITRI,ISIDE,ipls,itor) /= 0._dp) then
    ! sheath transmission as in styx (Tetri first in J)
      Tec=Interp_data2%tri_fluxE(ITRI,ISIDE,0,itor)/(4.5_dp*Interp_Data2%tri_fluxN(ITRI,ISIDE,ipls,itor))       			
      Mach=1.d0
      Tic=(Interp_data2%tri_fluxE(ITRI,ISIDE,ipls,itor)-0.5_dp*Mach*Interp_data2%tri_fluxN(ITRI,ISIDE,ipls,itor)*Tec)/ &
                     ((2.5_dp+0.5_dp*Mach)*Interp_Data2%tri_fluxN(ITRI,ISIDE,ipls,itor))

      Tec2=0.5_dp*(Interp_data2%knots_temperature(IP1,0,itor)+Interp_data2%knots_temperature(IP2,0,itor))
      Tic2=0.5_dp*(Interp_data2%knots_temperature(IP1,ipls,itor)+Interp_data2%knots_temperature(IP2,ipls,itor))

      Tec2=Tec2*reference_parameters%fields%T0*kB
      Tic2=Tic2*reference_parameters%fields%T0*kB
      ! to eVs
      Tec=Tec/elcha
      Tic=Tic/elcha
      Tec2=Tec2/elcha
      Tic2=Tic2/elcha

      ! correction if Tetri and Titri is too high
      !  --> possible in region where tri_fluxn very very small and flu_energy just very small
      ! threshold for choice set to 5
      if (Tec.gt.5._dp*Tec2) then
        Tec=Tec2
      end if
      if(Tic.gt.5._dp*Tic2) then
        Tic=Tic2
      end if                    

      if (Tec<0._dp) then
        Tec=1e-2_dp
        neg_Te=.true.
      endif
      if (Tic<0._dp) then
        Tic=1e-2_dp
        neg_Ti=.true.
      endif

      cstri=sqrt((Tec+Tic)*elcha/(global_parameters%element_list(isp)%mass*m_u))*100._dp

      ! use cs as velocity perpendicular to the wall (Bohm-Chodura)
      csx=cstri*PTRIX(ISIDE,ITRI)*Mach
      csy=cstri*PTRIY(ISIDE,ITRI)*Mach
      csz=0._dp
    else
      Tec=1._dp
      Tic=1._dp
      cstri=sqrt((Tec+Tic)*elcha/(global_parameters%element_list(isp)%mass*m_u))*100._dp
      csx=cstri*PTRIX(ISIDE,ITRI)
      csy=cstri*PTRIY(ISIDE,ITRI)
      csz=0._dp
    endif

  endif

  de=0.5_dp*(Interp_data2%knots_density(IP1,ipls,itor)+Interp_data2%knots_density(IP2,ipls,itor))
  de=de*reference_parameters%fields%n0*1d-6

  if (de<0) then
    de=1.e10_dp
    neg_de=.true.
  endif

  ! conversion of the particle flux (part.s-1 -> Amp)
  pflux=abs(Interp_data2%tri_fluxN(ITRI,ISIDE,ipls,itor))*ELCHA

  ! avoid sending NaNs to EIRENE
  if (isnan(Tec)) then
    write(*,*) ' Tec isnan, isurf = ',isurf
    call eirene_exit_own(1)
  endif 
       
  if (isnan(Tic)) then
    write(*,*) ' Tic isnan, isurf = ',isurf
    call eirene_exit_own(1)
  endif 

  if (isnan(De)) then
    write(*,*) ' De isnan, isurf = ',isurf
    call eirene_exit_own(1)
  endif 

  if (isnan(csx)) then
    write(*,*) 'csx isnan, isurf = ',isurf
    write(*,*) ' Te = ',Tec
    write(*,*) ' Ti = ',Tic
    call eirene_exit_own(1)
  endif

  if (isnan(csy)) then
    write(*,*) ' csy isnan, isurf = ',isurf
    call eirene_exit_own(1)
  endif

  if (isnan(csz)) then
    write(*,*) ' csz isnan, isurf = ',isurf
    call eirene_exit_own(1)
  endif

  if (neg_de .and. icountDe ==0) then
    write(*,*) ' negative density defaulted to 1e10m-3 in plasma_parameters_on_wall ...'
    icountDe=1
  endif

  if (neg_Te .and. icountTe == 0) then
    write(*,*) ' negative Te defaulted to 0.01 eV in plasma_parameters_on_wall ...'
    icountTe=1
  endif     

  if (neg_Ti .and. icountTi == 0) then
    write(*,*) ' negative Ti defaulted to 0.01 eV in plasma_parameters_on_wall ...'
    icountTi=1
  endif

  ! finally, output

  dens=de
  Te=Tec
  Ti=Tic
 

end subroutine
