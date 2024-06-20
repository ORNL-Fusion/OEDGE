subroutine styx_compute_pump_surface()
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_clgin
  use eirmod_ccona
  use styx2eirene

  implicit none

  integer*4 ipump, isurf, itri
  real*8 warea,pi_
  integer*4 ind

  pi_=4.D0*atan(1.d0)

  do ipump=1,n_pumps
     pumps(ipump)%Surface=0.d0
     pumps(ipump)%Ntriangles=0
     do isurf=nsurf0,nsurf_tal
        if (surface(isurf)%iprop == ipump+1) then
           ! ds in m, R1/R2 in cm...
           warea = pi_*(surface(isurf)%R1+surface(isurf)%R2)*surface(isurf)%ds*1e-2_dp
           pumps(ipump)%Surface=pumps(ipump)%Surface+warea
           pumps(ipump)%Ntriangles=pumps(ipump)%Ntriangles+1
        end if
     end do
     allocate(pumps(ipump)%triangle_index(1:pumps(ipump)%Ntriangles))
     allocate(pumps(ipump)%triangle_weight(1:pumps(ipump)%Ntriangles))
  end do

  do ipump=1,n_pumps
     ind=0
     do isurf=nsurf0,nsurf_tal
        if (surface(isurf)%iprop == ipump+1) then
           itri=surface(isurf)%itri
           ! ds in m, R1/R2 in cm...
           warea = pi_*(surface(isurf)%R1+surface(isurf)%R2)*surface(isurf)%ds*1e-2_dp
           ind=ind+1
           pumps(ipump)%triangle_index(ind)=itri
           pumps(ipump)%triangle_weight(ind)=warea/pumps(ipump)%Surface
        end if
     end do
  end do

  do ipump=1,n_pumps
     !number of eirene species (including molecules etc)
     allocate(pumps(ipump)%species_temperature(1:NSPZ))
     write(600+ipump,*) pumps(ipump)%surface
     write(600+ipump,*) pumps(ipump)%Ntriangles
  end do

end subroutine styx_compute_pump_surface
