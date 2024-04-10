subroutine interpolate_wall_fluxes()
#include "compile_opt.inc"
  use all_variables, only : global_parameters, interp_data2, zones, reference_parameters
  use Mphysics
  implicit none
  real*8,allocatable :: flux_tot_out(:)
  real*8,allocatable :: flux_totE_out(:)
  integer*4 :: i,j,k,ntri,n,nface,nion
  real*8 :: R0,rs0,n0,c0,T0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  n0=reference_parameters%fields%n0
  c0=reference_parameters%fields%c0
  T0=reference_parameters%fields%T0
  allocate(flux_tot_out(1:global_parameters%N_ions))
  allocate(flux_totE_out(0:global_parameters%N_ions))
  !reset
  interp_data2%tri_fluxn=0.d0
  flux_tot_out=0.d0
  interp_data2%tri_fluxE=0.d0
  flux_totE_out=0.d0

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,ntri,nface,i,j,k,nion)
  !$OMP DO
  do n=1,interp_data2%wall_data%n_triangles_on_wall
     ntri=interp_data2%wall_data%back_interp_on_wall(n,1)
     do nface=1,3
        if(interp_data2%type_face(ntri,nface).ne.0) then
           !north
           i=interp_data2%wall_data%s2d_to_use(n,1,1)
           j=interp_data2%wall_data%s2d_to_use(n,1,2)
           k=interp_data2%wall_data%s2d_to_use(n,1,3)
           if(k.ne.0) then
              do nion=1,global_parameters%N_ions
                 interp_data2%tri_fluxn(ntri,nface,nion,1)=interp_data2%tri_fluxn(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,1)*(Zones(k)%species(nion)%fluxes%fluxn(i,j,1)&
                      *zones(k)%metric_coefficients%ds_north_DD(i,j)*n0*c0*rs0**2/R0)
              end do
              do nion=0,global_parameters%N_ions
                 interp_data2%tri_fluxE(ntri,nface,nion,1)=interp_data2%tri_fluxE(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,1)*(Zones(k)%species(nion)%fluxes%fluxE(i,j,1)&
                      *zones(k)%metric_coefficients%ds_north_DD(i,j)*kb*n0*c0*T0*rs0**2/R0)
              end do
           end if
           !south
           i=interp_data2%wall_data%s2d_to_use(n,2,1)
           j=interp_data2%wall_data%s2d_to_use(n,2,2)
           k=interp_data2%wall_data%s2d_to_use(n,2,3)
           if(k.ne.0) then
              do nion=1,global_parameters%N_ions
                 interp_data2%tri_fluxn(ntri,nface,nion,1)=interp_data2%tri_fluxn(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,2)*(-Zones(k)%species(nion)%fluxes%fluxn(i,j,2)&
                      *zones(k)%metric_coefficients%ds_south_DD(i,j)*n0*c0*rs0**2/R0)
              end do
              do nion=0,global_parameters%N_ions
                 interp_data2%tri_fluxE(ntri,nface,nion,1)=interp_data2%tri_fluxE(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,2)*(-Zones(k)%species(nion)%fluxes%fluxE(i,j,2)&
                      *zones(k)%metric_coefficients%ds_south_DD(i,j)*kb*n0*c0*T0*rs0**2/R0)
              end do
           end if
           !east
           i=interp_data2%wall_data%s2d_to_use(n,3,1)
           j=interp_data2%wall_data%s2d_to_use(n,3,2)
           k=interp_data2%wall_data%s2d_to_use(n,3,3)
           if(k.ne.0) then
              do nion=1,global_parameters%N_ions
                 interp_data2%tri_fluxn(ntri,nface,nion,1)=interp_data2%tri_fluxn(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,3)*(Zones(k)%species(nion)%fluxes%fluxn(i,j,3)&
                      *zones(k)%metric_coefficients%ds_east_DD(i,j)*n0*c0*rs0**2/R0)
              end do
              do nion=0,global_parameters%N_ions
                 interp_data2%tri_fluxE(ntri,nface,nion,1)=interp_data2%tri_fluxE(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,3)*(Zones(k)%species(nion)%fluxes%fluxE(i,j,3)&
                      *zones(k)%metric_coefficients%ds_east_DD(i,j)*kb*n0*c0*T0*rs0**2/R0)
              end do
           end if
           !west
           i=interp_data2%wall_data%s2d_to_use(n,4,1)
           j=interp_data2%wall_data%s2d_to_use(n,4,2)
           k=interp_data2%wall_data%s2d_to_use(n,4,3)
           if(k.ne.0) then
              do nion=1,global_parameters%N_ions
                 interp_data2%tri_fluxn(ntri,nface,nion,1)=interp_data2%tri_fluxn(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,4)*(-Zones(k)%species(nion)%fluxes%fluxn(i,j,4)&
                      *zones(k)%metric_coefficients%ds_west_DD(i,j)*n0*c0*rs0**2/R0)
              end do
              do nion=0,global_parameters%N_ions
                 interp_data2%tri_fluxE(ntri,nface,nion,1)=interp_data2%tri_fluxE(ntri,nface,nion,1)+&
                      interp_data2%wall_data%weight_s2d_to_use(n,4)*(-Zones(k)%species(nion)%fluxes%fluxE(i,j,4)&
                      *zones(k)%metric_coefficients%ds_west_DD(i,j)*kb*n0*c0*T0*rs0**2/R0)
              end do
           end if
        end if
     end do

  end do
  !$OMP END DO
  !$OMP BARRIER

  !$OMP DO
  do n=1,interp_data2%n_triangles
     do nface=1,3
        if(interp_data2%type_face(n,nface).ne.0) then
           do nion=1,global_parameters%N_ions
              if (interp_data2%tri_fluxN(n,nface,nion,1) < 0) then
                 interp_data2%tri_fluxN(n,nface,nion,1)=0.d0
#if VERBOSE==1
                 write(*,*) 'negative particle outflux, ntri = ',n,' nface = ',nface
#endif
              endif
           end do
           do nion=0,global_parameters%N_ions
              if (interp_data2%tri_fluxE(n,nface,nion,1) < 0) then
                 interp_data2%tri_fluxE(n,nface,nion,1)=0.d0
                 !                write(*,*) 'negative electron energy outflux, ntri = ',n,' nface = ',nface
              endif
           end do
        endif
     enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER

  !$OMP DO REDUCTION(+:flux_tot_out,flux_totE_out)
  do n=1,interp_data2%n_triangles
     do nface=1,3
        do nion=1,global_parameters%N_ions
           flux_tot_out(nion)=flux_tot_out(nion)+interp_data2%tri_fluxn(n,nface,nion,1)
        end do
        do nion=0,global_parameters%N_ions
           flux_totE_out(nion)=flux_totE_out(nion)+interp_data2%tri_fluxE(n,nface,nion,1)
        end do
     end do
  end do
  !$OMP END DO

  !$OMP MASTER
  open(unit=680,file='fluxn',status='unknown')
  open(unit=690,file='fluxE',status='unknown')
  do n=1,interp_data2%wall_data%n_triangles_on_wall
     ntri=interp_data2%wall_data%back_interp_on_wall(n,1)
     do nface=1,3
        if(interp_data2%type_face(ntri,nface).ne.0) then 
           write(680,100) ntri, nface, (interp_data2%tri_fluxn(ntri,nface,nion,1),nion=0,global_parameters%N_ions)
           write(690,100) ntri, nface, (interp_data2%tri_fluxE(ntri,nface,nion,1),nion=0,global_parameters%N_ions)
100        format(I6,I6,512es15.7)
        end if
     end do
  end do
  close(680)
  close(690)
  !$OMP END MASTER

  !$OMP END PARALLEL
end subroutine interpolate_wall_fluxes
