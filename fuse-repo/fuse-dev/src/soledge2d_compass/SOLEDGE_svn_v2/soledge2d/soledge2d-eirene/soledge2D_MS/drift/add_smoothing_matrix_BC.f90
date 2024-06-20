subroutine add_smoothing_matrix_BC()
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables, penalisation_parameters
  use Msmoothing_vars
  use Mlist
  implicit none
  integer*4 i,j,k
  integer*4 Nx,Nz
  integer*4 North,South,East,West

  do k=1,global_parameters%N_zones

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
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,j)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,j)
           sm_ptr%val=1.D0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,j)
           sm_ptr%col=Zones(North)%mesh%index(1,j)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else
!!$           if(North.eq.-1) then !n=0
!!$              allocate(sm_ptr) ; ptr_prev%p => sm_ptr
!!$              sm_ptr%row=zones(k)%mesh%index(Nx+1,j)
!!$              sm_ptr%col=zones(k)%mesh%index(Nx+1,j)
!!$              sm_ptr%val=1.D0
!!$              ptr_prev => sm_ptr
!!$              mat_smoothing_nnz=mat_smoothing_nnz+1
!!$           else !grad = 0
              allocate(sm_ptr) ; ptr_prev%p => sm_ptr
              sm_ptr%row=zones(k)%mesh%index(Nx+1,j)
              sm_ptr%col=zones(k)%mesh%index(Nx+1,j)
              sm_ptr%val=1.d0
              ptr_prev => sm_ptr
              mat_smoothing_nnz=mat_smoothing_nnz+1
              allocate(sm_ptr) ; ptr_prev%p => sm_ptr
              sm_ptr%row=zones(k)%mesh%index(Nx+1,j)
              sm_ptr%col=zones(k)%mesh%index(Nx,j)
              sm_ptr%val=-1.D0
              ptr_prev => sm_ptr
              mat_smoothing_nnz=mat_smoothing_nnz+1
!           end if
        end if
     end do
     !treating corner NE
     if(North.gt.0) then
        if(East.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=Zones(North)%mesh%index(1,Nz+1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=Zones(East)%mesh%index(Nx+1,1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=Zones(North)%mesh%index(1,Nz+1)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        end if
     else
        if(East.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=Zones(East)%mesh%index(Nx+1,1)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else !real corner
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,Nz)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(Nx,Nz+1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
        end if
     end if

     !treating corner NW
     if(North.gt.0) then
        if(West.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=Zones(North)%mesh%index(1,0)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=Zones(West)%mesh%index(Nx+1,Zones(West)%mesh%Nz)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=Zones(North)%mesh%index(1,0)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        end if
     else
        if(West.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=Zones(West)%mesh%index(Nx+1,Zones(West)%mesh%Nz)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else !real corner
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=zones(k)%mesh%index(Nx+1,1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(Nx+1,0)
           sm_ptr%col=zones(k)%mesh%index(Nx,0)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
        end if
     end if

     !South
     do j=1,Nz
        if(South.gt.0) then                  
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,j)
           sm_ptr%col=zones(k)%mesh%index(0,j)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,j)
           sm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,j)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else
!!$           if(South.eq.-1) then !n=0
!!$              allocate(sm_ptr) ; ptr_prev%p => sm_ptr
!!$              sm_ptr%row=zones(k)%mesh%index(0,j)
!!$              sm_ptr%col=zones(k)%mesh%index(0,j)
!!$              sm_ptr%val=1.D0
!!$              ptr_prev => sm_ptr
!!$              mat_smoothing_nnz=mat_smoothing_nnz+1
!!$           else !grad=0
              allocate(sm_ptr) ; ptr_prev%p => sm_ptr
              sm_ptr%row=zones(k)%mesh%index(0,j)
              sm_ptr%col=zones(k)%mesh%index(0,j)
              sm_ptr%val=-1.D0
              ptr_prev => sm_ptr
              mat_smoothing_nnz=mat_smoothing_nnz+1
              allocate(sm_ptr) ; ptr_prev%p => sm_ptr
              sm_ptr%row=zones(k)%mesh%index(0,j)
              sm_ptr%col=zones(k)%mesh%index(1,j)
              sm_ptr%val=1.D0
              ptr_prev => sm_ptr
              mat_smoothing_nnz=mat_smoothing_nnz+1
!           end if
        end if
     end do
     !treating corner SE
     if(South.gt.0) then
        if(East.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,Nz+1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=Zones(East)%mesh%index(0,1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,Nz+1)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        end if
     else
        if(East.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=Zones(East)%mesh%index(0,1)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else !real corner
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(0,Nz)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(1,Nz+1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
        end if
     end if

     !treating corner SW
     if(South.gt.0) then
        if(West.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=zones(k)%mesh%index(0,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,0)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=Zones(West)%mesh%index(0,Zones(West)%mesh%Nz)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=zones(k)%mesh%index(0,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=Zones(South)%mesh%index(Zones(South)%mesh%Nx,0)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        end if
     else
        if(West.gt.0) then
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=zones(k)%mesh%index(0,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=Zones(West)%mesh%index(0,Zones(West)%mesh%Nz)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else !real corner
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=zones(k)%mesh%index(0,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=zones(k)%mesh%index(0,1)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(0,0)
           sm_ptr%col=zones(k)%mesh%index(1,0)
           sm_ptr%val=-0.5d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1             
        end if
     end if

     !East
     do i=1,Nx
        if(East.gt.0) then                  
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(i,Nz+1)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,Nz+1)
           sm_ptr%col=Zones(East)%mesh%index(i,1)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else !grad=0
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(i,Nz+1)
           sm_ptr%val=1.D0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,Nz+1)
           sm_ptr%col=zones(k)%mesh%index(i,Nz)
           sm_ptr%val=-1.D0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        end if
     end do

     !West
     do i=1,Nx
        if(West.gt.0) then                  
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,0)
           sm_ptr%col=zones(k)%mesh%index(i,0)
           sm_ptr%val=1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,0)
           sm_ptr%col=Zones(West)%mesh%index(i,Zones(West)%mesh%Nz)
           sm_ptr%val=-1.d0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        else !grad=0
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,0)
           sm_ptr%col=zones(k)%mesh%index(i,0)
           sm_ptr%val=-1.D0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
           allocate(sm_ptr) ; ptr_prev%p => sm_ptr
           sm_ptr%row=zones(k)%mesh%index(i,0)
           sm_ptr%col=zones(k)%mesh%index(i,1)
           sm_ptr%val=1.D0
           ptr_prev => sm_ptr
           mat_smoothing_nnz=mat_smoothing_nnz+1
        end if
     end do
  end do
end subroutine add_smoothing_matrix_BC
