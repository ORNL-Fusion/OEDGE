subroutine interpolate(temp,num)
  use all_variables, only : interp_data2, zones, global_parameters
  use Minterpolation_types
  implicit none
  integer*4,intent(in) :: num
  Type(interpolation_temp),intent(inout) :: temp
  integer*4 k,n
  integer*4 i,j
  integer*4 i1,j1,k1
  integer*4 i2,j2,k2
  integer*4 i3,j3,k3
  integer*4 nneigh
  real*8 d,inter
  real*8 sumX(num),sumd
  integer*4 ii,ij,ik
  real*8 X1,X2,X3
  integer*4 nst,nk
  real*8 :: dist(1:3), val(1:3)
  integer*4 pt
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nk,sumX,sumd,ik,ii,ij,i1,j1,k1,i2,j2,k2,i3,j3,k3,nneigh,d,inter,i,j,k,n,nst,X1,X2,X3,dist,val,pt)
  !$OMP DO
  do nk=1,interp_data2%N_knots
     if(interp_data2%knots_interp_points(nk)%pass.eq.5) then        
        ! centre pass
        sumX=0.D0
        sumd=0.D0
        do ik=1,global_parameters%N_zones
           if(Zones(ik)%Neighbors(1).eq.-5) then
              do ij=1,Zones(ik)%mesh%Nz
                 do nst=1,num
                    sumX(nst)=sumX(nst)+temp%zones(ik)%val(Zones(ik)%mesh%Nx,ij,nst)
                 end do
                 sumd=sumd+1.D0
              end do
           end if
        end do
        do nst=1,num
           temp%knots_val(nk,nst)=sumX(nst)/sumd
        end do
     end if
  end do
  !$OMP END DO
  !$OMP DO
  do nk=1,interp_data2%N_knots
     if(interp_data2%knots_interp_points(nk)%pass.eq.1) then
        ! first pass
        nneigh=interp_data2%knots_interp_points(nk)%n_soledge
        sumX=0.d0
        sumd=0.d0
        do n=1,nneigh !##CAREFUL## the case nneigh=1 is quite raw
           ii=interp_data2%knots_interp_points(nk)%sol(n,1)
           ij=interp_data2%knots_interp_points(nk)%sol(n,2)
           ik=interp_data2%knots_interp_points(nk)%sol(n,3)
           d=sqrt((interp_data2%knots_r(nk)-zones(ik)%mesh%Rgeom(ii,ij))**(2.d0)+&
                (interp_data2%knots_z(nk)-zones(ik)%mesh%Zgeom(ii,ij))**(2.d0))
           do nst=1,num
              sumX(nst)=sumX(nst)+1.d0/d*temp%zones(ik)%val(ii,ij,nst)
           end do
           sumd=sumd+1.d0/d
        end do
        do nst=1,num
           temp%knots_val(nk,nst)=sumX(nst)/sumd
        end do
     end if
  end do
  !$OMP END DO
  !$OMP DO
  do nk=1,interp_data2%N_knots
     if(interp_data2%knots_interp_points(nk)%pass.eq.2) then
        !second pass
        i1=interp_data2%knots_interp_points(nk)%sol(1,1)
        j1=interp_data2%knots_interp_points(nk)%sol(1,2)
        k1=interp_data2%knots_interp_points(nk)%sol(1,3)
        i2=interp_data2%knots_interp_points(nk)%sol(2,1)
        j2=interp_data2%knots_interp_points(nk)%sol(2,2)
        k2=interp_data2%knots_interp_points(nk)%sol(2,3)
        i3=interp_data2%knots_interp_points(nk)%sol(3,1)
        j3=interp_data2%knots_interp_points(nk)%sol(3,2)
        k3=interp_data2%knots_interp_points(nk)%sol(3,3)
        do nst=1,num
           X1=temp%zones(k1)%val(i1,j1,nst)
           X2=temp%zones(k2)%val(i2,j2,nst)
           X3=temp%zones(k3)%val(i3,j3,nst)
           call project(nk,X1,X2,X3,inter)
           temp%knots_val(nk,nst)=inter
        end do
     end if
  end do
  !$OMP END DO
  !$OMP DO
  do nk=1,interp_data2%N_knots
     if(interp_data2%knots_interp_points(nk)%pass.eq.3) then
        !third pass
        if(interp_data2%knots_interp_points(nk)%n_soledge.eq.2) then
           i1=interp_data2%knots_interp_points(nk)%sol(1,1)
           j1=interp_data2%knots_interp_points(nk)%sol(1,2)
           k1=interp_data2%knots_interp_points(nk)%sol(1,3)
           i2=interp_data2%knots_interp_points(nk)%sol(2,1)
           j2=interp_data2%knots_interp_points(nk)%sol(2,2)
           k2=interp_data2%knots_interp_points(nk)%sol(2,3)
           do nst=1,num
              X1=temp%zones(k1)%val(i1,j1,nst)
              X2=temp%zones(k2)%val(i2,j2,nst)
              X3=temp%knots_val(interp_data2%knots_interp_points(nk)%eir(1),nst)
              call project(nk,X1,X2,X3,inter)
              ! choose closest value
              dist(1)=abs(X1-inter)
              dist(2)=abs(X2-inter)
              dist(3)=abs(X3-inter)
              val(1)=X1
              val(2)=X2
              val(3)=X3
              pt=minloc(dist,1)
              inter=val(pt)
              temp%knots_val(nk,nst)=inter
           end do
        else
           i=interp_data2%knots_interp_points(nk)%sol(1,1)
           j=interp_data2%knots_interp_points(nk)%sol(1,2)
           k=interp_data2%knots_interp_points(nk)%sol(1,3)
           do nst=1,num
              X1=temp%zones(k)%val(i,j,nst)
              X2=temp%knots_val(interp_data2%knots_interp_points(nk)%eir(1),nst)
              X3=temp%knots_val(interp_data2%knots_interp_points(nk)%eir(2),nst)
              call project(nk,X1,X2,X3,inter)
              temp%knots_val(nk,nst)=inter
           end do
        end if
     end if
  end do
  !$OMP END DO
  !$OMP DO
  do nk=1,interp_data2%N_knots
     if(interp_data2%knots_interp_points(nk)%pass.eq.4) then
        !fourth pass
        if(interp_data2%knots_interp_points(nk)%n_soledge.eq.2) then
           i1=interp_data2%knots_interp_points(nk)%sol(1,1)
           j1=interp_data2%knots_interp_points(nk)%sol(1,2)
           k1=interp_data2%knots_interp_points(nk)%sol(1,3)
           i2=interp_data2%knots_interp_points(nk)%sol(2,1)
           j2=interp_data2%knots_interp_points(nk)%sol(2,2)
           k2=interp_data2%knots_interp_points(nk)%sol(2,3)
           do nst=1,num
              X1=temp%zones(k1)%val(i1,j1,nst)
              X2=temp%zones(k2)%val(i2,j2,nst)
              temp%knots_val(nk,nst)=(X1+X2)*0.5d0
           end do
        else
           i=interp_data2%knots_interp_points(nk)%sol(1,1)
           j=interp_data2%knots_interp_points(nk)%sol(1,2)
           k=interp_data2%knots_interp_points(nk)%sol(1,3)
           do nst=1,num
              X1=temp%zones(k)%val(i,j,nst)
              temp%knots_val(nk,nst)=X1
           end do
        end if
     end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine interpolate
