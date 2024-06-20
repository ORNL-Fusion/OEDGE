subroutine add_vorticity_matrix_BC4(mode_BC,mode_BC_para)
  use all_variables, only : global_parameters, global_variables, reference_parameters, zones
  use Mvorticity_vars
  use Mlist
  use Mphysics
  implicit none

  integer*4,intent(in) :: mode_BC
  integer*4,intent(in) :: mode_BC_para
  integer*4 i,j,k
  integer*4 Nx,Nz
  integer*4 North,South,East,West
  real*8 beta
  real*8 Lambda,r1,r2

  real*8 :: R0,rs0,Gamma
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0

  do k=1,global_parameters%N_Zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz

     !Equations for the neighboring points
     North=Zones(k)%Neighbors(1)
     South=Zones(k)%Neighbors(2)
     East=Zones(k)%Neighbors(3)
     West=Zones(k)%Neighbors(4)
     !North
     do j=1,Nz
        if(North.gt.0) then                  
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,j)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,j)
           vm_ptr%val=1.D0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,j)
           vm_ptr%col=Zones(North)%mesh%index(1,j)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else !North == -1 or North == -2
           if(mode_BC.eq.ZERO_FLUX) then
              !Delta phi = - (Delta p_i) /n or whatever
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(Nx+1,j)
              vm_ptr%col=Zones(k)%mesh%index(Nx+1,j)
              vm_ptr%val=1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(Nx+1,j)
              vm_ptr%col=Zones(k)%mesh%index(Nx,j)
              vm_ptr%val=-1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           else
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(Nx+1,j)
              vm_ptr%col=Zones(k)%mesh%index(Nx+1,j)
              vm_ptr%val=1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           end if
        end if
     end do
     !treating corner NE
     if(North.gt.0) then
        if(East.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(North)%mesh%index(1,Nz+1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(East)%mesh%index(Nx+1,1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(North)%mesh%index(1,Nz+1)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        end if
     else
        if(East.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(East)%mesh%index(Nx+1,1)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else !real corner
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,Nz)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(Nx,Nz+1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
        end if
     end if

     !treating corner NW
     if(North.gt.0) then
        if(West.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(North)%mesh%index(1,0)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(West)%mesh%index(Nx+1,Zones(West)%mesh%Nz)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(North)%mesh%index(1,0)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        end if
     else
        if(West.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(West)%mesh%index(Nx+1,Zones(West)%mesh%Nz)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else !real corner
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(k)%mesh%index(Nx+1,1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(Nx+1,0)
           vm_ptr%col=Zones(k)%mesh%index(Nx,0)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
        end if
     end if

     !South
     do j=1,Nz
        if(South.gt.0) then                  
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,j)
           vm_ptr%col=Zones(k)%mesh%index(0,j)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,j)
           vm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,j)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else ! South == -1 or South==-2
           if(mode_BC.eq.ZERO_FLUX) then
              !Flux=0
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(0,j)
              vm_ptr%col=Zones(k)%mesh%index(0,j)
              vm_ptr%val=1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(0,j)
              vm_ptr%col=Zones(k)%mesh%index(1,j)
              vm_ptr%val=-1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           else
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(0,j)
              vm_ptr%col=Zones(k)%mesh%index(0,j)
              vm_ptr%val=1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           end if
        end if
     end do
     !treating corner SE
     if(South.gt.0) then
        if(East.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,Nz+1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(East)%mesh%index(0,1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,Nz+1)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        end if
     else
        if(East.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(East)%mesh%index(0,1)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else !real corner
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(0,Nz)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(1,Nz+1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
        end if
     end if

     !treating corner SW
     if(South.gt.0) then
        if(West.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(k)%mesh%index(0,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,0)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(West)%mesh%index(0,Zones(West)%mesh%Nz)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(k)%mesh%index(0,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,0)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        end if
     else
        if(West.gt.0) then
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(k)%mesh%index(0,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(West)%mesh%index(0,Zones(West)%mesh%Nz)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else !real corner
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(k)%mesh%index(0,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(k)%mesh%index(0,1)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(0,0)
           vm_ptr%col=Zones(k)%mesh%index(1,0)
           vm_ptr%val=-0.5d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1             
        end if
     end if

     !East
     do i=1,Nx
        if(East.gt.0) then                  
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(i,Nz+1)
           vm_ptr%col=Zones(k)%mesh%index(i,Nz+1)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(i,Nz+1)
           vm_ptr%col=Zones(East)%mesh%index(i,1)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else ! Bohm boundary condition
           if(mode_BC_para.eq.LAMBDA_TE) then
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,Nz+1)
              vm_ptr%col=Zones(k)%mesh%index(i,Nz+1)
              vm_ptr%val=1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           else
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,Nz+1)
              vm_ptr%col=Zones(k)%mesh%index(i,Nz+1)
              Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+max(min(Zones(k)%species(1)%var(1)%temperature(i,Nz)&
                   /Zones(k)%species(0)%var(1)%temperature(i,Nz),10.D0),0.1D0)))
              beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                   *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
              Gamma=zones(k)%species(0)%var(1)%Gamma(i,Nz)+(zones(k)%species(0)%drifts%uEt(i,Nz)+zones(k)%species(0)%drifts%uBt(i,Nz))&
                   *(2.d0*pi*R0/rs0)*sqrt(zones(k)%metric_coefficients%ctt(i,Nz))/zones(k)%metric_coefficients%G(i,Nz)&
                   *zones(k)%species(0)%var(1)%density(i,Nz)
              r1=-beta*Gamma/zones(k)%species(0)%var(1)%temperature(i,Nz)
              r2=-(zones(k)%metric_coefficients%G(i,Nz+1)+zones(k)%metric_coefficients%G(i,Nz))*0.5D0&
                   *(zones(k)%species(0)%var(1)%temperature(i,Nz)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,Nz+1)**(1.5d0))*0.5D0&
                   /(zones(k)%mesh%z(i,Nz+1)-zones(k)%mesh%z(i,Nz))
              vm_ptr%val=r1+r2
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,Nz+1)
              vm_ptr%col=Zones(k)%mesh%index(i,Nz)
              vm_ptr%val=(zones(k)%metric_coefficients%G(i,Nz+1)+zones(k)%metric_coefficients%G(i,Nz))*0.5D0&
                   *(zones(k)%species(0)%var(1)%temperature(i,Nz)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,Nz+1)**(1.5d0))*0.5D0&
                   /(zones(k)%mesh%z(i,Nz+1)-zones(k)%mesh%z(i,Nz))
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           end if
        end if
     end do

     !West
     do i=1,Nx
        if(West.gt.0) then                  
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(i,0)
           vm_ptr%col=Zones(k)%mesh%index(i,0)
           vm_ptr%val=1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
           allocate(vm_ptr) ; ptr_prev%p => vm_ptr
           vm_ptr%row=Zones(k)%mesh%index(i,0)
           vm_ptr%col=Zones(West)%mesh%index(i,Zones(West)%mesh%Nz)
           vm_ptr%val=-1.d0
           ptr_prev => vm_ptr
           mat_vort_nnz=mat_vort_nnz+1
        else !Bohm boundary condition:
           if(mode_BC_para.eq.LAMBDA_TE) then
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,0)
              vm_ptr%col=Zones(k)%mesh%index(i,0)
              vm_ptr%val=1.D0
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           else
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,0)
              vm_ptr%col=Zones(k)%mesh%index(i,0)
              Lambda=-0.5d0*log(2.D0*pi*m_e/(global_parameters%element_list(1)%mass*m_u)*(1.D0+max(min(Zones(k)%species(1)%var(1)%temperature(i,1)&
                   /Zones(k)%species(0)%var(1)%temperature(i,1),10.D0),0.1D0)))
              beta=reference_parameters%fields%c0*reference_parameters%fields%n0*eV*8.e-4*reference_parameters%fields%T0eV**(-1.5d0)&
                   *2.D0*pi*reference_parameters%geometry%R0/reference_parameters%fields%phi0
              !                    r1=-beta*zones(k)%species(0)%var(1)%Gamma(i,j+1)&
              !                         *exp(Lambda-zones(k)%electric_fields(1)%phi(i,j+1)/zones(k)%species(0)%var(1)%temperature(i,j+1))&
              !                         /zones(k)%species(0)%var(1)%temperature(i,j+1)
              Gamma=zones(k)%species(0)%var(1)%Gamma(i,1)+(zones(k)%species(0)%drifts%uEt(i,1)+zones(k)%species(0)%drifts%uBt(i,1))&
                   *(2.d0*pi*R0/rs0)*sqrt(zones(k)%metric_coefficients%ctt(i,1))/zones(k)%metric_coefficients%G(i,1)&
                   *zones(k)%species(0)%var(1)%density(i,1)
              r1=-beta*Gamma/zones(k)%species(0)%var(1)%temperature(i,1)
              r2=-(zones(k)%metric_coefficients%G(i,0)+zones(k)%metric_coefficients%G(i,1))*0.5D0&
                   *(zones(k)%species(0)%var(1)%temperature(i,1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,0)**(1.5d0))*0.5D0&
                   /(zones(k)%mesh%z(i,0)-zones(k)%mesh%z(i,1))
              vm_ptr%val=r1+r2
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
              allocate(vm_ptr) ; ptr_prev%p => vm_ptr
              vm_ptr%row=Zones(k)%mesh%index(i,0)
              vm_ptr%col=Zones(k)%mesh%index(i,1)
              vm_ptr%val=(zones(k)%metric_coefficients%G(i,0)+zones(k)%metric_coefficients%G(i,1))*0.5D0&
                   *(zones(k)%species(0)%var(1)%temperature(i,1)**(1.5d0)+zones(k)%species(0)%var(1)%temperature(i,0)**(1.5d0))*0.5D0&
                   /(zones(k)%mesh%z(i,0)-zones(k)%mesh%z(i,1))
              ptr_prev => vm_ptr
              mat_vort_nnz=mat_vort_nnz+1
           end if
        end if
     end do

  end do
end subroutine add_vorticity_matrix_BC4
