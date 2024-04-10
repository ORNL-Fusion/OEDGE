subroutine report_source_in_pen_cells_to_closest_plasma_cells()
  use all_variables, only : global_parameters, zones, interp_data2
  integer*4 i,j,k,n
  integer*4 i_,j_,k_
  integer*4 Nx,Nz
  real*8 Vtot
  integer*4 :: nion,ind
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,Nx,Nz,i,j,Vtot,n,i_,j_,k_,nion,ind)
  !$OMP DO
  do k=1,global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           !Eirene for neutrals
           if(Zones(k)%masks%chi2(i,j).eq.1) then
              Vtot=0.d0
              do n=1,Zones(k)%masks%npts_around_penwall(i,j)
                 i_=Zones(k)%masks%pts_around_penwall(i,j,n,1)
                 j_=Zones(k)%masks%pts_around_penwall(i,j,n,2)
                 k_=Zones(k)%masks%pts_around_penwall(i,j,n,3)
                 Vtot=Vtot+zones(k_)%metric_coefficients%dvol_pu(i_,j_)
              end do
              ! reporting source in non-penalized cells
              do n=1,Zones(k)%masks%npts_around_penwall(i,j)
                 i_=Zones(k)%masks%pts_around_penwall(i,j,n,1)
                 j_=Zones(k)%masks%pts_around_penwall(i,j,n,2)
                 k_=Zones(k)%masks%pts_around_penwall(i,j,n,3)
                 do nion=1,global_parameters%n_species
                    ind=global_parameters%ind_ion_0(nion)
                    Zones(k_)%species(ind)%sources%Sn_n(i_,j_)= Zones(k_)%species(ind)%sources%Sn_n(i_,j_)&
                         +Zones(k)%species(ind)%sources%Sn_n(i,j)*zones(k)%metric_coefficients%dvol_pu(i,j)/Vtot
                    Zones(k_)%species(ind)%sources%Sn_G(i_,j_)= Zones(k_)%species(ind)%sources%Sn_G(i_,j_)&
                         +Zones(k)%species(ind)%sources%Sn_G(i,j)*zones(k)%metric_coefficients%dvol_pu(i,j)/Vtot
                    Zones(k_)%species(ind)%sources%Sn_E(i_,j_)= Zones(k_)%species(ind)%sources%Sn_E(i_,j_)&
                         +Zones(k)%species(ind)%sources%Sn_E(i,j)*zones(k)%metric_coefficients%dvol_pu(i,j)/Vtot
                 end do
                 Zones(k_)%species(0)%sources%Sn_E(i_,j_)= Zones(k_)%species(0)%sources%Sn_E(i_,j_)&
                      +Zones(k)%species(0)%sources%Sn_E(i,j)*zones(k)%metric_coefficients%dvol_pu(i,j)/Vtot
              end do
              do nion=1,global_parameters%n_species
                 ind=global_parameters%ind_ion_0(nion)
                 Zones(k)%species(ind)%sources%Sn_n(i,j)=0.d0 
                 Zones(k)%species(ind)%sources%Sn_G(i,j)=0.d0 
                 Zones(k)%species(ind)%sources%Sn_E(i,j)=0.d0 
              end do
              Zones(k)%species(0)%sources%Sn_E(i,j)=0.d0 
           end if
           !                end if
        end do !i
     end do !j
  end do !k
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine report_source_in_pen_cells_to_closest_plasma_cells
