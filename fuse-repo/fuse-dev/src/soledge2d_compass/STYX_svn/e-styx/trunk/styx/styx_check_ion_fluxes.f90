subroutine styx_check_ion_fluxes
  use all_variables, only : interp_data2
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use eirmod_ccona
  use eirmod_cgrid
  use eirmod_comusr
  implicit none
  integer :: isurf,ipls,ip1,ip2,itri,iside,k,isou,kr,itor
  character(3) :: nplsic
  character(20) :: format
  real*8, allocatable :: pflux_s2d(:),pflux_eir(:)
  real*8 :: pi_,warea

! compares for each ion species the particle fluxes from styx
! to the flux sampled by EIRENE

  pi_=4._dp*atan(1._dp)

  write(nplsic,'(i3)') nplsi
  format='('//trim(adjustl(nplsic))//'e14.7)'

  allocate(pflux_s2d(nplsi),pflux_eir(nplsi))

  open(unit=666,file='fluxes_in.txt',status='replace')
  open(unit=667,file='fluxes_sm.txt',status='replace')

  do itor=1,Ntor_cells
    do isou=1,Nsou
      ! get vertex number for the segment making the wall
      kr=krecsurf(isou)
      IP1   = recsurf(kr)%v1
      IP2   = recsurf(kr)%v2
      ! corresponding triangle number and side
      ITRI  = recsurf(kr)%ITRI
      ISIDE = recsurf(kr)%ISIDE

      isurf=isou+nsurf0
      k=ksurf(isurf)

      warea = ang_max*pi_/180._dp*(surface(k)%R1+surface(k)%R2)*0.5_dp*surface(k)%ds*1e-2_dp
  
      do ipls=1,nplsi  
        ! fluxes in part/s
        pflux_s2d(ipls)=abs(Interp_data2%tri_fluxN(ITRI,ISIDE,ipls,itor))/warea
        pflux_eir(ipls)=eos(0)%potpl(ipls,k)/warea
      enddo
      write(666,format) (pflux_s2d(ipls),ipls=1,nplsi)
      write(667,format) (pflux_eir(ipls),ipls=1,nplsi)
    enddo     
  enddo
  
  close(666)
  close(667)
  deallocate(pflux_s2d,pflux_eir)

end subroutine styx_check_ion_fluxes
