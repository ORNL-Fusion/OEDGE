subroutine correct_perp_velocity(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: n,Nx,Nz,i,j,np
  integer*4 :: i1,j1,k1
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           if(zone%masks%chi2(i,j).eq.0) then
              if(.not.zone%DriftsExtra%OK_points(i,j)) then
                 zone%species(n)%drifts%uEp(i,j)=0.d0
                 zone%species(n)%drifts%uEt(i,j)=0.d0
                 zone%species(n)%drifts%uBp(i,j)=0.d0
                 zone%species(n)%drifts%uBt(i,j)=0.d0
                 do np=1,zone%DriftsExtra%Nref_points(i,j)
                    i1=zone%DriftsExtra%ref_points(i,j,1,np)
                    j1=zone%DriftsExtra%ref_points(i,j,2,np)
                    k1=zone%DriftsExtra%ref_points(i,j,3,np)
                    zone%species(n)%drifts%uEp(i,j)=zone%species(n)%drifts%uEp(i,j)&
                         +zone%species(n)%drifts%uEp(i1,j1)&
                         *sqrt(zone%metric_coefficients%c_pp(i1,j1)/zone%metric_coefficients%c_pp(i,j))
                    zone%species(n)%drifts%uEt(i,j)=zone%species(n)%drifts%uEt(i,j)&
                         +zone%species(n)%drifts%uEt(i1,j1)&
                         *sqrt(zone%metric_coefficients%c_tt(i1,j1)/zone%metric_coefficients%c_tt(i,j))
                    !idem for uB
                    zone%species(n)%drifts%uBp(i,j)=zone%species(n)%drifts%uBp(i,j)&
                         +zone%species(n)%drifts%uBp(i1,j1)&
                         *sqrt(zone%metric_coefficients%c_pp(i1,j1)/zone%metric_coefficients%c_pp(i,j))
                    zone%species(n)%drifts%uBt(i,j)=zone%species(n)%drifts%uBt(i,j)&
                         +zone%species(n)%drifts%uBt(i1,j1)&
                         *sqrt(zone%metric_coefficients%c_tt(i1,j1)/zone%metric_coefficients%c_tt(i,j))

                 end do
                 zone%species(n)%drifts%uEp(i,j)=zone%species(n)%drifts%uEp(i,j)&
                      /real(zone%DriftsExtra%Nref_points(i,j))
                 zone%species(n)%drifts%uEt(i,j)=zone%species(n)%drifts%uEt(i,j)&
                      /real(zone%DriftsExtra%Nref_points(i,j))
                 zone%species(n)%drifts%uBp(i,j)=zone%species(n)%drifts%uBp(i,j)&
                      /real(zone%DriftsExtra%Nref_points(i,j))
                 zone%species(n)%drifts%uBt(i,j)=zone%species(n)%drifts%uBt(i,j)&
                      /real(zone%DriftsExtra%Nref_points(i,j))

              end if
           end if
        end do
     end do
  end do
end subroutine correct_perp_velocity
