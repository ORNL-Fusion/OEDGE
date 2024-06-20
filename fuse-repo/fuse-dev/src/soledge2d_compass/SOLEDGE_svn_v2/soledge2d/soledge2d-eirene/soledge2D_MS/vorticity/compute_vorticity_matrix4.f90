subroutine compute_vorticity_matrix4(mode_BC,mode_BC_para)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables, penalisation_parameters, drift_flags
  use Mvorticity_vars
  use Mlist
  use MDefinitions
  use Mphysics
  implicit none

  integer*4,intent(in) :: mode_BC, mode_BC_para

  integer*4 i,j,k
  integer*4 Nx,Nz
  real*8 Vol,SE,SW,SN,SS
  real*8 SE_mod,SW_mod,SN_mod,SS_mod
  real*8,allocatable :: D(:,:),D2(:,:)
  real*8,allocatable :: x(:,:),z(:,:)
  real*8,allocatable :: x2(:,:),z2(:,:)
  real*8,allocatable :: cpp(:,:),cpt(:,:),ctt(:,:),G(:,:)
  real*8 delta
  real*8 A,eta_para0,dt,eta_perp0
  real*8 etaT
  real*8 :: rs0, R0
  real*8 :: coef_east, coef_west
  real*8,parameter :: eps=1.d-6
  real*8 :: Teps
  real*8 :: Lambda,beta,r1,r2

  !if(mode_BC_para.eq.LAMBDA_TE) then
  !   write(*,*) 'Parallel BC for potential (matrix) : Lambda*Te'
  !else
  !   write(*,*) 'Parallel BC for potential (matrix) : Bohm current'
  !end if

  Teps=global_variables%Teps
  mat_vort_nnz=0

  etaT=penalisation_parameters%eta

  eta_para0=8.e-4*reference_parameters%fields%T0eV**(-1.5d0)
  delta=(reference_parameters%fields%tau0*reference_parameters%geometry%rs0)/&
       (reference_parameters%fields%n0*m_u*eta_para0*(2.d0*pi*reference_parameters%geometry%R0))

  ! delta=(reference_parameters%fields%tau0*reference_parameters%geometry%rs0)/&
  !      (reference_parameters%fields%n0*m_u*eta_para0*(reference_parameters%geometry%rs0))

  dt=global_variables%dt_vort

  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0

  do k=1,global_parameters%N_zones

     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     allocate(D(0:Nx+1,0:Nz+1))
     allocate(D2(0:Nx+1,0:Nz+1))
     allocate(x(0:Nx+1,0:Nz+1))
     allocate(z(0:Nx+1,0:Nz+1))
     allocate(x2(0:Nx+1,0:Nz+1))
     allocate(z2(0:Nx+1,0:Nz+1))
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

     x2=zones(k)%mesh%x
     z2=zones(k)%mesh%z
     x2(Nx+1,Nz+1)=zones(k)%mesh%cornerx(1,2,2)
     z2(Nx+1,Nz+1)=zones(k)%mesh%cornerz(1,2,2)
     x2(0,Nz+1)=zones(k)%mesh%cornerx(2,2,2)
     z2(0,Nz+1)=zones(k)%mesh%cornerz(2,2,2)
     x2(Nx+1,0)=zones(k)%mesh%cornerx(1,1,2)
     z2(Nx+1,0)=zones(k)%mesh%cornerz(1,1,2)
     x2(0,0)=zones(k)%mesh%cornerx(2,1,2)
     z2(0,0)=zones(k)%mesh%cornerz(2,1,2)

     cpp=zones(k)%metric_coefficients%cpp
     cpt=zones(k)%metric_coefficients%cpt
     ctt=zones(k)%metric_coefficients%ctt
     G=zones(k)%metric_coefficients%G
     A=reference_parameters%geometry%A    

     eta_perp0=drift_flags%eta_perp0
     do i=0,Nx+1
        do j=0,Nz+1
           D2(i,j)=max(Zones(k)%species(0)%var(STEP_NEW)%temperature(i,j),Teps)**(1.5d0) 
           D(i,j)=Zones(k)%species(0)%var(STEP_NEW)%density(i,j)/(zones(k)%mesh%B(i,j)**2+eps)*global_parameters%element_list(1)%mass
           D(i,j)= D(i,j)-eta_perp0*dt
        end do
     end do


     do i=1,Nx
        do j=1,Nz
           if(zones(k)%masks%chi2(i,j).eq.0) then
              coef_east = &
                   ((D(i,j)*(ctt(i,j)-1.d0/A**2.d0*G(i,j)**2.d0)&
                   )&
                   +(D(i,j+1)*(ctt(i,j+1)-1.d0/A**2.d0*G(i,j+1)**2.d0)&
                   ))*0.5d0&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))
              coef_west = &
                   ((D(i,j)*(ctt(i,j)-1.d0/A**2.d0*G(i,j)**2.d0)&
                   )&
                   +(D(i,j-1)*(ctt(i,j-1)-1.d0/A**2.d0*G(i,j-1)**2.d0)&
                   ))*0.5d0&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))

              Vol=zones(k)%metric_coefficients%dvol_dd(i,j)
              SN=zones(k)%metric_coefficients%ds_north_dd(i,j)
              SS=zones(k)%metric_coefficients%ds_south_dd(i,j)
              SE=zones(k)%metric_coefficients%ds_east_dd(i,j)
              SW=zones(k)%metric_coefficients%ds_west_dd(i,j)

              if(mode_BC.eq.ZERO_FLUX) then
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
              else ! dirichlet
                 SN_mod=SN
                 SS_mod=SS
                 SE_mod=SE
                 SW_mod=SW
              end if

!!$           SN_mod=SN
!!$           SS_mod=SS
!!$           SE_mod=SE
!!$           SW_mod=SW

!!$           div=1/Vol*&
!!$                ((
!!$               *((Field(i,j)-Field(i+1,j))/(x(i,j)-x(i+1,j))&
!!$               *0.5d0*(cpp(i,j)+cpp(i+1,j))&
!!$               +0.5d0*((Field(i,j+1)-Field(i,j-1))/(z(i,j+1)-z(i,j-1))&
!!$               *cpt(i,j)&
!!$               +(Field(i+1,j+1)-Field(i+1,j-1))/(z(i+1,j+1)-z(i+1,j-1))&
!!$               *cpt(i+1,j))&
!!$               )/sqrt(0.5D0*(cpp(i,j)+cpp(i+1,j))))&
!!$               *SN*(0.5d0*(D(i+1,j)+D(i,j)))&
!!$               -(
!!$               *((Field(i,j)-Field(i-1,j))/(x(i,j)-x(i-1,j))&
!!$               *0.5d0*(cpp(i,j)+cpp(i-1,j))&
!!$               +0.5d0*((Field(i,j+1)-Field(i,j-1))/(z(i,j+1)-z(i,j-1))&
!!$               *cpt(i,j)&
!!$               +(Field(i-1,j+1)-Field(i-1,j-1))/(z(i-1,j+1)-z(i-1,j-1))&
!!$               *cpt(i-1,j))&
!!$               )/sqrt(0.5D0*(cpp(i,j)+cpp(i-1,j))))&
!!$               *SS&
!!$               +((
!!$               *(0.5d0*((Field(i+1,j)-Field(i-1,j))&
!!$               /(x(i+1,j)-x(i-1,j))*cpt(i,j)&
!!$               +(Field(i+1,j+1)-Field(i-1,j+1))&
!!$               /(x(i+1,j+1)-x(i-1,j+1))*cpt(i,j+1)))&
!!$               /sqrt(0.5d0*(ctt(i,j)+ctt(i,j+1))))&
!!$               +coef_east(i,j)*(Field(i,j+1)-Field(i,j))&
!!$                /(z(i,j+1)-z(i,j)))&
!!$               *SE&
!!$               -((
!!$               *(0.5d0*((Field(i+1,j)-Field(i-1,j))&
!!$               /(x(i+1,j)-x(i-1,j))*cpt(i,j)&
!!$               +(Field(i+1,j-1)-Field(i-1,j-1))&
!!$               /(x(i+1,j-1)-x(i-1,j-1))*cpt(i,j-1)))&
!!$               /sqrt(0.5d0*(ctt(i,j)+ctt(i,j-1))))&
!!$               +coef_west(i,j)*(Field(i,j)-Field(i,j-1))&
!!$                /(z(i,j)-z(i,j-1)))&
!!$               *SW)

!!$           div_para=1/vol*(&
!!$             (G(i,j)+G(i,j+1))*0.5D0&
!!$             *(Field(i,j+1)-Field(i,j))&
!!$             /(z(i,j+1)-z(i,j))&
!!$             *delta*(D2(i,j)+D2(i,j+1))*0.5D0*zone%metric_coefficients%sinepitch_east(i,j)*SE&
!!$             -(G(i,j)+G(i,j-1))*0.5D0&
!!$             *(Field(i,j)-Field(i,j-1))&
!!$             /(z(i,j)-z(i,j-1))&
!!$             *delta*(D2(i,j)+D2(i,j-1))*0.5D0*zone%metric_coefficients%sinepitch_west(i,j)*SW)



              !point i,j
              if(.not.associated(vorticity_mat)) then
                 allocate(vorticity_mat) ; vm_ptr => vorticity_mat
              else
                 allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              end if
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i,j)
              vm_ptr%val=1.d0/vol/(x(i,j)-x(i+1,j))&
                   *0.5d0*(cpp(i,j)+cpp(i+1,j))&
                   /(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i+1,j))))*SN_mod*(0.5d0*(D(i+1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   -1.d0/vol/(x(i,j)-x(i-1,j))&
                   *0.5d0*(cpp(i,j)+cpp(i-1,j))&
                   /(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i-1,j))))*SS_mod*(0.5d0*(D(i-1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   -1/vol*coef_east/(z2(i,j+1)-z2(i,j))*SE_mod
              vm_ptr%val=vm_ptr%val&
                   -1/vol*coef_west/(z2(i,j)-z2(i,j-1))*SW_mod
              vm_ptr%val=vm_ptr%val
              !j perp
              !j parallel
              vm_ptr%val=vm_ptr%val&
                   -1.d0/vol*(G(i,j)+G(i,j+1))*0.5d0/(z2(i,j+1)-z2(i,j))&
                   *delta*(D2(i,j)+D2(i,j+1))*0.5d0&
                   *zones(k)%metric_coefficients%sinepitch_east(i,j)*SE*dt
              vm_ptr%val=vm_ptr%val&
                   -1.d0/vol*(G(i,j)+G(i,j-1))*0.5d0/(z2(i,j)-z2(i,j-1))&
                   *delta*(D2(i,j)+D2(i,j-1))*0.5d0&
                   *zones(k)%metric_coefficients%sinepitch_west(i,j)*SW*dt

              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i+1,j
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i+1,j)
              vm_ptr%val=-1.d0/vol/(x(i,j)-x(i+1,j))&
                   *0.5d0*(cpp(i,j)+cpp(i+1,j))&
                   /(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i+1,j))))*SN_mod*(0.5d0*(D(i+1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   +1/vol*0.5d0/(x2(i+1,j)-x2(i-1,j))*cpt(i,j)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))*SE_mod*(0.5d0*(D(i,j+1)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   -1.d0/vol*0.5d0/(x2(i+1,j)-x2(i-1,j))*cpt(i,j)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))*SW_mod*(0.5d0*(D(i,j-1)+D(i,j)))
              vm_ptr%val=vm_ptr%val
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i-1,j
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i-1,j)
              vm_ptr%val=1.d0/vol/(x(i,j)-x(i-1,j))&
                   *0.5d0*(cpp(i,j)+cpp(i-1,j))&
                   /(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i-1,j))))*SS_mod*(0.5d0*(D(i-1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   -1/vol*0.5d0/(x2(i+1,j)-x2(i-1,j))*cpt(i,j)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))*SE_mod*(0.5d0*(D(i,j+1)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   +1.d0/vol*0.5d0/(x2(i+1,j)-x2(i-1,j))*cpt(i,j)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))*SW_mod*(0.5d0*(D(i,j-1)+D(i,j)))
              vm_ptr%val=vm_ptr%val
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i,j+1
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i,j+1)
              vm_ptr%val=1/vol*0.5d0/(z(i,j+1)-z(i,j-1))&
                   *cpt(i,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i+1,j))))*SN_mod*(0.5d0*(D(i+1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   -1/vol*0.5d0/(z(i,j+1)-z(i,j-1))&
                   *cpt(i,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i-1,j))))*SS_mod*(0.5d0*(D(i-1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   +1/vol*coef_east/(z2(i,j+1)-z2(i,j))*SE_mod
              vm_ptr%val=vm_ptr%val
              !j parallel
              vm_ptr%val=vm_ptr%val&
                   +1.d0/vol*(G(i,j)+G(i,j+1))*0.5d0/(z2(i,j+1)-z2(i,j))&
                   *delta*(D2(i,j)+D2(i,j+1))*0.5d0&
                   *zones(k)%metric_coefficients%sinepitch_east(i,j)*SE*dt
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i,j-1
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i,j-1)
              vm_ptr%val=-1/vol*0.5d0/(z(i,j+1)-z(i,j-1))&
                   *cpt(i,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i+1,j))))*SN_mod*(0.5d0*(D(i+1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   +1/vol*0.5d0/(z(i,j+1)-z(i,j-1))&
                   *cpt(i,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i-1,j))))*SS_mod*(0.5d0*(D(i-1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   +1/vol*coef_west/(z2(i,j)-z2(i,j-1))*SW_mod
              vm_ptr%val=vm_ptr%val
              !j parallel
              vm_ptr%val=vm_ptr%val&
                   +1.d0/vol*(G(i,j)+G(i,j-1))*0.5d0/(z2(i,j)-z2(i,j-1))&
                   *delta*(D2(i,j)+D2(i,j-1))*0.5d0&
                   *zones(k)%metric_coefficients%sinepitch_west(i,j)*SW*dt
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i+1,j+1
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i+1,j+1)
              vm_ptr%val=1/vol*0.5d0/(z(i+1,j+1)-z(i+1,j-1))&
                   *cpt(i+1,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i+1,j))))*SN_mod*(0.5d0*(D(i+1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   +1/vol*0.5d0/(x2(i+1,j+1)-x2(i-1,j+1))*cpt(i,j+1)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))*SE_mod*(0.5d0*(D(i,j+1)+D(i,j)))
              vm_ptr%val=vm_ptr%val
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i+1,j-1
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i+1,j-1)
              vm_ptr%val=-1/vol*0.5d0/(z(i+1,j+1)-z(i+1,j-1))&
                   *cpt(i+1,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i+1,j))))*SN_mod*(0.5d0*(D(i+1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   -1.d0/vol*0.5d0/(x2(i+1,j-1)-x2(i-1,j-1))*cpt(i,j-1)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))*SW_mod*(0.5d0*(D(i,j-1)+D(i,j)))
              vm_ptr%val=vm_ptr%val
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i-1,j+1
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i-1,j+1)
              vm_ptr%val=-1/vol*0.5d0/(z(i-1,j+1)-z(i-1,j-1))&
                   *cpt(i-1,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i-1,j))))*SS_mod*(0.5d0*(D(i-1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   -1/vol*0.5d0/(x2(i+1,j+1)-x2(i-1,j+1))*cpt(i,j+1)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j+1))))*SE_mod*(0.5d0*(D(i,j+1)+D(i,j)))
              vm_ptr%val=vm_ptr%val
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              !point i-1,j-1
              allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i-1,j-1)
              vm_ptr%val=1/vol*0.5d0/(z(i-1,j+1)-z(i-1,j-1))&
                   *cpt(i-1,j)/(0.5d0*(sqrt(cpp(i,j))+sqrt(cpp(i-1,j))))*SS_mod*(0.5d0*(D(i-1,j)+D(i,j)))
              vm_ptr%val=vm_ptr%val&
                   +1.d0/vol*0.5d0/(x2(i+1,j-1)-x2(i-1,j-1))*cpt(i,j-1)&
                   /(0.5d0*(sqrt(ctt(i,j))+sqrt(ctt(i,j-1))))*SW_mod*(0.5d0*(D(i,j-1)+D(i,j)))
              vm_ptr%val=vm_ptr%val
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1


           else

              !point i,j
              if(.not.associated(vorticity_mat)) then
                 allocate(vorticity_mat) ; vm_ptr => vorticity_mat
              else
                 allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
              end if
              vm_ptr%row=Zones(k)%mesh%index(i,j)
              vm_ptr%col=Zones(k)%mesh%index(i,j)
              vm_ptr%val=1.d0
              if(zones(k)%masks%chi1(i,j).eq.1) then
                 if(mode_BC_para.ne.LAMBDA_TE) then
                    Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+max(min(Zones(k)%species(1)%var(1)%temperature(i,j-1)&
                         /Zones(k)%species(0)%var(1)%temperature(i,j-1),10.D0),0.1D0)))
                    beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                         *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
!                    r1=-beta*zones(k)%species(0)%var(1)%Gamma(i,j-1)&
!                         *exp(Lambda-zones(k)%electric_fields(1)%phi(i,j-1)/zones(k)%species(0)%var(1)%temperature(i,j-1))&
!                         /zones(k)%species(0)%var(1)%temperature(i,j-1)
                    r1=-beta*zones(k)%species(0)%var(1)%Gamma(i,j-1)/zones(k)%species(0)%var(1)%temperature(i,j-1)
                    r2=-(zones(k)%metric_coefficients%G(i,j)+zones(k)%metric_coefficients%G(i,j-1))*0.5D0&
                         *(zones(k)%species(0)%var(1)%temperature(i,j-1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,j)**(1.5d0))*0.5D0&
                         /(zones(k)%mesh%z(i,j)-zones(k)%mesh%z(i,j-1))
                    vm_ptr%val=r1+r2
                 end if
              end if
              if(zones(k)%masks%chi3(i,j).eq.1) then
                 if(mode_BC_para.ne.LAMBDA_TE) then
                    Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+max(min(Zones(k)%species(1)%var(1)%temperature(i,j+1)&
                         /Zones(k)%species(0)%var(1)%temperature(i,j+1),10.D0),0.1D0)))
                    beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                         *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
!                    r1=-beta*zones(k)%species(0)%var(1)%Gamma(i,j+1)&
!                         *exp(Lambda-zones(k)%electric_fields(1)%phi(i,j+1)/zones(k)%species(0)%var(1)%temperature(i,j+1))&
!                         /zones(k)%species(0)%var(1)%temperature(i,j+1)
                    r1=-beta*zones(k)%species(0)%var(1)%Gamma(i,j+1)/zones(k)%species(0)%var(1)%temperature(i,j+1)
                    r2=-(zones(k)%metric_coefficients%G(i,j)+zones(k)%metric_coefficients%G(i,j+1))*0.5D0&
                         *(zones(k)%species(0)%var(1)%temperature(i,j+1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,j)**(1.5d0))*0.5D0&
                         /(zones(k)%mesh%z(i,j)-zones(k)%mesh%z(i,j+1))
                    vm_ptr%val=r1+r2
                 end if
              end if
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1

              if(zones(k)%masks%chi1(i,j).eq.1) then
                 if(mode_BC_para.ne.LAMBDA_TE) then
                    !point i,j-1
                    if(.not.associated(vorticity_mat)) then
                       allocate(vorticity_mat) ; vm_ptr => vorticity_mat
                    else
                       allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
                    end if
                    vm_ptr%row=Zones(k)%mesh%index(i,j)
                    vm_ptr%col=Zones(k)%mesh%index(i,j-1)
                    vm_ptr%val=(zones(k)%metric_coefficients%G(i,j)+zones(k)%metric_coefficients%G(i,j-1))*0.5D0&
                         *(zones(k)%species(0)%var(1)%temperature(i,j-1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,j)**(1.5d0))*0.5D0&
                         /(zones(k)%mesh%z(i,j)-zones(k)%mesh%z(i,j-1))
                    ptr_prev => vm_ptr
                    mat_vort_nnz=mat_vort_nnz+1
                 end if
              end if

              if(zones(k)%masks%chi3(i,j).eq.1) then
                 if(mode_BC_para.ne.LAMBDA_TE) then
                    !point i,j+1
                    if(.not.associated(vorticity_mat)) then
                       allocate(vorticity_mat) ; vm_ptr => vorticity_mat
                    else
                       allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
                    end if
                    vm_ptr%row=Zones(k)%mesh%index(i,j)
                    vm_ptr%col=Zones(k)%mesh%index(i,j+1)
                    vm_ptr%val=(zones(k)%metric_coefficients%G(i,j)+zones(k)%metric_coefficients%G(i,j+1))*0.5D0&
                         *(zones(k)%species(0)%var(1)%temperature(i,j+1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,j)**(1.5d0))*0.5D0&
                         /(zones(k)%mesh%z(i,j)-zones(k)%mesh%z(i,j+1))
                    ptr_prev => vm_ptr
                    mat_vort_nnz=mat_vort_nnz+1
                 end if
              end if
              
              if(zones(k)%masks%chi5(i,j).eq.1) then
                 !point i+1,j
                 if(.not.associated(vorticity_mat)) then
                    allocate(vorticity_mat) ; vm_ptr => vorticity_mat
                 else
                    allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
                 end if
                 vm_ptr%row=Zones(k)%mesh%index(i,j)
                 vm_ptr%col=Zones(k)%mesh%index(i+1,j)
                 vm_ptr%val=-1.d0
                 ptr_prev => vm_ptr
                 mat_vort_nnz=mat_vort_nnz+1                 
              end if

              if(zones(k)%masks%chi6(i,j).eq.1) then
                 !point i-1,j
                 if(.not.associated(vorticity_mat)) then
                    allocate(vorticity_mat) ; vm_ptr => vorticity_mat
                 else
                    allocate(vm_ptr) ; ptr_prev%p=>vm_ptr
                 end if
                 vm_ptr%row=Zones(k)%mesh%index(i,j)
                 vm_ptr%col=Zones(k)%mesh%index(i-1,j)
                 vm_ptr%val=-1.d0
                 ptr_prev => vm_ptr
                 mat_vort_nnz=mat_vort_nnz+1                 
              end if

           end if
        end do
     end do

     deallocate(D,D2,x,z,x2,z2,cpp,cpt,ctt,G)
  end do
end subroutine compute_vorticity_matrix4
