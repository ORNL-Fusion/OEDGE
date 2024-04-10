subroutine compute_laplacian()
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables, penalisation_parameters
  use Mneutral_vars
  use Mlist
  use MDefinitions
  use Mphysics
  implicit none
  integer*4 i,j,k
  integer*4 Nx,Nz
  real*8 Vol,SE,SW,SN,SS
  real*8 SE_mod,SW_mod,SN_mod,SS_mod
  real*8,allocatable :: D(:,:)
  real*8,allocatable :: x(:,:),z(:,:)
  real*8,allocatable :: cpp(:,:),cpt(:,:),ctt(:,:),G(:,:)
  real*8 A,dt
  real*8 :: rs0, R0
  real*8 :: coef_east, coef_west
  real*8 :: DE,DW,DS,DN
  mat_neutral_nnz=0
  dt=global_variables%dt
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0

  do k=1,global_parameters%N_zones

     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz

     allocate(D(0:Nx+1,0:Nz+1))
     allocate(x(0:Nx+1,0:Nz+1))
     allocate(z(0:Nx+1,0:Nz+1))
     allocate(cpp(0:Nx+1,0:Nz+1))
     allocate(cpt(0:Nx+1,0:Nz+1))
     allocate(ctt(0:Nx+1,0:Nz+1))
     allocate(G(0:Nx+1,0:Nz+1))

     x=zones(k)%mesh%x
     z=zones(k)%mesh%z
     x(Nx+1,Nz+1)=zones(k)%mesh%cornerx(1,2,1)
     z(Nx+1,Nz+1)=zones(k)%mesh%cornerz(1,2,1)
     x(0,Nz+1)=zones(k)%mesh%cornerx(2,2,1)
     z(0,Nz+1)=zones(k)%mesh%cornerz(2,2,1)
     x(Nx+1,0)=zones(k)%mesh%cornerx(1,1,1)
     z(Nx+1,0)=zones(k)%mesh%cornerz(1,1,1)
     x(0,0)=zones(k)%mesh%cornerx(2,1,1)
     z(0,0)=zones(k)%mesh%cornerz(2,1,1)

     cpp=zones(k)%metric_coefficients%cpp
     cpt=zones(k)%metric_coefficients%cpt
     ctt=zones(k)%metric_coefficients%ctt
     G=zones(k)%metric_coefficients%G

     do i=0,Nx+1
        do j=0,Nz+1
           D(i,j)=zones(k)%neutrals%Dn(i,j)
        end do
     end do

     do i=1,Nx
        do j=1,Nz
           if(zones(k)%masks%chi2(i,j).eq.0) then
              Vol=zones(k)%metric_coefficients%dvol_dd(i,j)
              SN=zones(k)%metric_coefficients%ds_north_dd(i,j)
              SS=zones(k)%metric_coefficients%ds_south_dd(i,j)
              SE=zones(k)%metric_coefficients%ds_east_dd(i,j)
              SW=zones(k)%metric_coefficients%ds_west_dd(i,j)

              if(zones(k)%Neighbors(1).lt.-1) then
                 if(i.eq.zones(k)%mesh%Nx) then
                    SN_mod=0.d0
                 else
                    SN_mod=zones(k)%metric_coefficients%ds_north_dd(i,j)
                 end if
              else
                 SN_mod=zones(k)%metric_coefficients%ds_north_dd(i,j)
              end if
              if(zones(k)%Neighbors(2).lt.-1) then
                 if(i.eq.1) then
                    SS_mod=0.d0
                 else
                    SS_mod=zones(k)%metric_coefficients%ds_south_dd(i,j)
                 end if
              else
                 SS_mod=zones(k)%metric_coefficients%ds_south_dd(i,j)
              end if
              if(zones(k)%Neighbors(3).lt.0) then
                 if(j.eq.zones(k)%mesh%Nz) then
                    SE_mod=0.d0
                 else
                    SE_mod=zones(k)%metric_coefficients%ds_east_dd(i,j)
                 end if
              else
                 SE_mod=zones(k)%metric_coefficients%ds_east_dd(i,j)
              end if
              if(zones(k)%Neighbors(4).lt.0) then
                 if(j.eq.1) then
                    SW_mod=0.d0
                 else
                    SW_mod=zones(k)%metric_coefficients%ds_west_dd(i,j)
                 end if
              else
                 SW_mod=zones(k)%metric_coefficients%ds_west_dd(i,j)
              end if

              SN_mod=SN_mod*(1.D0-zones(k)%masks%chi2(i+1,j))
              SS_mod=SS_mod*(1.D0-zones(k)%masks%chi2(i-1,j))
              SE_mod=SE_mod*(1.D0-zones(k)%masks%chi2(i,j+1))
              SW_mod=SW_mod*(1.D0-zones(k)%masks%chi2(i,j-1))

              DN=(D(i+1,j)+D(i,j))*0.5d0
              DS=(D(i-1,j)+D(i,j))*0.5d0
              DE=(D(i,j+1)+D(i,j))*0.5d0
              DW=(D(i,j-1)+D(i,j))*0.5d0


              !point i,j
              if(.not.associated(neutral_mat)) then 
                 allocate(neutral_mat) ; nm_ptr => neutral_mat
              else
                 allocate(nm_ptr) ; ptr_prev%p=>nm_ptr
              end if
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i,j)
              nm_ptr%val=DE/Vol/(z(i,j+1)-z(i,j))&
                   *(ctt(i,j)+ctt(i,j+1))*0.5d0&
                   /sqrt((ctt(i,j)+ctt(i,j+1))*0.5d0)*SE_mod
              nm_ptr%val=nm_ptr%val&
                   +DW/Vol/(z(i,j)-z(i,j-1))&
                   *(ctt(i,j)+ctt(i,j-1))*0.5d0&
                   /sqrt((ctt(i,j)+ctt(i,j-1))*0.5d0)*SW_mod
              nm_ptr%val=nm_ptr%val&
                   +DN/Vol/(x(i+1,j)-x(i,j))&
                   *(cpp(i+1,j)+cpp(i,j))*0.5d0&
                   /sqrt((cpp(i+1,j)+cpp(i,j))*0.5d0)*SN_mod
              nm_ptr%val=nm_ptr%val&
                   +DS/Vol/(x(i,j)-x(i-1,j))&
                   *(cpp(i-1,j)+cpp(i,j))*0.5d0&
                   /sqrt((cpp(i-1,j)+cpp(i,j))*0.5d0)*SS_mod
              nm_ptr%val=nm_ptr%val+1.D0/dt
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i,j+1
              allocate(nm_ptr) ; ptr_prev%p=>nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i,j+1)
              nm_ptr%val=-DE/Vol/(z(i,j+1)-z(i,j))&
                   *(ctt(i,j)+ctt(i,j+1))*0.5d0&
                   /sqrt((ctt(i,j)+ctt(i,j+1))*0.5d0)*SE_mod
              nm_ptr%val=nm_ptr%val&
                   -DN/Vol/(z(i,j+1)-z(i,j-1))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((cpp(i+1,j)+cpp(i,j))*0.5d0)*SN_mod
              nm_ptr%val=nm_ptr%val&
                   +DS/Vol/(z(i,j+1)-z(i,j-1))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((cpp(i-1,j)+cpp(i,j))*0.5d0)*SS_mod
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i,j-1
              allocate(nm_ptr) ; ptr_prev%p=>nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i,j-1)
              nm_ptr%val=-DW/Vol/(z(i,j)-z(i,j-1))&
                   *(ctt(i,j)+ctt(i,j-1))*0.5d0&
                   /sqrt((ctt(i,j)+ctt(i,j-1))*0.5d0)*SW_mod
              nm_ptr%val=nm_ptr%val&
                   +DN/Vol/(z(i,j+1)-z(i,j-1))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((cpp(i+1,j)+cpp(i,j))*0.5d0)*SN_mod
              nm_ptr%val=nm_ptr%val&
                   -DS/Vol/(z(i,j+1)-z(i,j-1))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((cpp(i-1,j)+cpp(i,j))*0.5d0)*SS_mod
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i+1,j
              allocate(nm_ptr) ; ptr_prev%p=>nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i+1,j)
              nm_ptr%val=-DN/Vol/(x(i+1,j)-x(i,j))&
                   *(cpp(i+1,j)+cpp(i,j))*0.5d0&
                   /sqrt((cpp(i+1,j)+cpp(i,j))*0.5d0)*SN_mod
              nm_ptr%val=nm_ptr%val&
                   -DE/Vol/(x(i+1,j)-x(i-1,j))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((ctt(i,j)+ctt(i,j+1))*0.5d0)*SE_mod
              nm_ptr%val=nm_ptr%val&
                   +DW/Vol/(x(i+1,j)-x(i-1,j))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((ctt(i,j)+ctt(i,j-1))*0.5d0)*SW_mod
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i-1,j
              allocate(nm_ptr) ; ptr_prev%p=>nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i-1,j)
              nm_ptr%val=-DS/Vol/(x(i,j)-x(i-1,j))&
                   *(cpp(i-1,j)+cpp(i,j))*0.5d0&
                   /sqrt((cpp(i-1,j)+cpp(i,j))*0.5d0)*SS_mod
              nm_ptr%val=nm_ptr%val&
                   +DE/Vol/(x(i+1,j)-x(i-1,j))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((ctt(i,j)+ctt(i,j+1))*0.5d0)*SE_mod
              nm_ptr%val=nm_ptr%val&
                   -DW/Vol/(x(i+1,j)-x(i-1,j))*0.5d0&
                   *cpt(i,j)&
                   /sqrt((ctt(i,j)+ctt(i,j-1))*0.5d0)*SW_mod
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i+1,j+1
              allocate(nm_ptr) ; ptr_prev%p=>nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i+1,j+1)
              nm_ptr%val=-DE/Vol/(x(i+1,j+1)-x(i-1,j+1))*0.5d0&
                   *cpt(i,j+1)&
                   /sqrt((ctt(i,j)+ctt(i,j+1))*0.5d0)*SE_mod
              nm_ptr%val=nm_ptr%val&
                   -DN/Vol/(z(i+1,j+1)-z(i+1,j-1))*0.5d0&
                   *cpt(i+1,j)&
                   /sqrt((cpp(i+1,j)+cpp(i,j))*0.5d0)*SN_mod
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i+1,j-1
              allocate(nm_ptr) ; ptr_prev%p => nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i+1,j-1)
              nm_ptr%val=DW/Vol/(x(i+1,j-1)-x(i-1,j-1))*0.5d0&
                   *cpt(i,j-1)&
                   /sqrt((ctt(i,j)+ctt(i,j-1))*0.5d0)*SW_mod
              nm_ptr%val=nm_ptr%val&
                   +DN/Vol/(z(i+1,j+1)-z(i+1,j-1))*0.5d0&
                   *cpt(i+1,j)&
                   /sqrt((cpp(i+1,j)+cpp(i,j))*0.5d0)*SN_mod
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i-1,j-1
              allocate(nm_ptr) ; ptr_prev%p => nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i-1,j-1)
              nm_ptr%val=-DW/Vol/(x(i+1,j-1)-x(i-1,j-1))*0.5d0&
                   *cpt(i,j-1)&
                   /sqrt((ctt(i,j)+ctt(i,j-1))*0.5d0)*SW_mod
              nm_ptr%val=nm_ptr%val&
                   -DS/Vol/(z(i-1,j+1)-z(i-1,j-1))*0.5d0&
                   *cpt(i-1,j)&
                   /sqrt((cpp(i-1,j)+cpp(i,j))*0.5d0)*SS_mod
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

              !point i-1,j+1
              allocate(nm_ptr) ; ptr_prev%p => nm_ptr
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i-1,j+1)
              nm_ptr%val=DE/Vol/(x(i+1,j+1)-x(i-1,j+1))*0.5d0&
                   *cpt(i,j+1)&
                   /sqrt((ctt(i,j)+ctt(i,j+1))*0.5d0)*SE_mod
              nm_ptr%val=nm_ptr%val&
                   +DS/Vol/(z(i-1,j+1)-z(i-1,j-1))*0.5d0&
                   *cpt(i-1,j)&
                   /sqrt((cpp(i-1,j)+cpp(i,j))*0.5d0)*SS_mod
              Ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1

           else

              !point i,j
              if(.not.associated(neutral_mat)) then 
                 allocate(neutral_mat) ; nm_ptr => neutral_mat
              else
                 allocate(nm_ptr) ; ptr_prev%p=>nm_ptr
              end if
              nm_ptr%row=Zones(k)%mesh%index(i,j)
              nm_ptr%col=Zones(k)%mesh%index(i,j)
              nm_ptr%val=1.D0
              ptr_prev => nm_ptr
              mat_neutral_nnz=mat_neutral_nnz+1
              
           end if

        end do
     end do

     deallocate(D,x,z,cpp,cpt,ctt,G)
  end do
end subroutine compute_laplacian

  
