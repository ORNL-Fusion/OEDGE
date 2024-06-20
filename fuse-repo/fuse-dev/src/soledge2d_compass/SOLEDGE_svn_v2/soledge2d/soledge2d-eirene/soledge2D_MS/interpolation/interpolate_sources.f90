subroutine interpolate_sources()
  use styx2eirene, only : vol_tri_eirene
  use all_variables, only : global_parameters, interp_data2, zones
  implicit none
  integer*4 :: k,i,j,Nx,Nz
  real*8 :: sumSn,sumSG,sumSE
  integer*4 :: itri,ntri,tri_num
  integer*4 :: nsp,ind_ion
  real*8 :: Vol
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,Nx,Nz,i,j,sumSn,sumSG,sumSE,ntri,Vol,itri,tri_num,ind_ion,nsp)
  !$OMP DO
  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           !ions
           do nsp=1,global_parameters%N_species
              ind_ion=global_parameters%ind_ion_0(nsp)
              sumSn=0.d0
              sumSG=0.d0
              sumSE=0.d0
              ntri = interp_data2%zones_data(k)%num_triangles(i,j)
              Vol=zones(k)%metric_coefficients%dvol_PU(i,j)
              do itri=1,ntri
                 tri_num = Interp_data2%zones_data(k)%triangles(i,j,itri)
                 sumSn=SumSn+Interp_data2%tri_Sn(tri_num,nsp,1)*vol_tri_eirene(tri_num)/vol
                 sumSG=SumSG+Interp_Data2%tri_SG(tri_num,nsp,1)*vol_tri_eirene(tri_num)/vol
                 sumSE=SumSE+Interp_Data2%tri_SE(tri_num,nsp,1)*vol_tri_eirene(tri_num)/vol
              enddo
              zones(k)%species(ind_ion)%sources%Sn_n(i,j)=sumSn
              zones(k)%species(ind_ion)%sources%Sn_G(i,j)=sumSG
              zones(k)%species(ind_ion)%sources%Sn_E(i,j)=sumSE
           end do
           !electrons
           sumSE=0.d0
           ntri = interp_data2%zones_data(k)%num_triangles(i,j)
           Vol=zones(k)%metric_coefficients%dvol_PU(i,j)
           do itri=1,ntri
              tri_num = Interp_data2%zones_data(k)%triangles(i,j,itri)
              sumSE=SumSE+Interp_Data2%tri_SE(tri_num,0,1)*vol_tri_eirene(tri_num)/vol
           enddo
           zones(k)%species(0)%sources%Sn_E(i,j)=sumSE
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine interpolate_sources
