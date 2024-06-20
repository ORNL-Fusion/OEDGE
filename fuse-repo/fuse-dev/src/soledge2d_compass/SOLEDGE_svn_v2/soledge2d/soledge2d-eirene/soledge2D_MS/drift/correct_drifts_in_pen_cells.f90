subroutine correct_drifts_in_pen_cells(zone)
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
           if(zone%masks%chi2(i,j).eq.1) then
              if(zone%masks%chi2(i,j+1).eq.0) then
                 zone%species(n)%drifts%uEt(i,j)=zone%species(n)%drifts%uEt(i,j+1)
                 zone%species(n)%drifts%uBt(i,j)=zone%species(n)%drifts%uBt(i,j+1)
              else
                 if(zone%masks%chi2(i,j-1).eq.0) then
                    zone%species(n)%drifts%uEt(i,j)=zone%species(n)%drifts%uEt(i,j-1)
                    zone%species(n)%drifts%uBt(i,j)=zone%species(n)%drifts%uBt(i,j-1)
                 else
                    zone%species(n)%drifts%uEt(i,j)=0.d0
                    zone%species(n)%drifts%uBt(i,j)=0.d0
                 end if
              end if
              if(zone%masks%chi2(i+1,j).eq.0) then
                 zone%species(n)%drifts%uEp(i,j)=zone%species(n)%drifts%uEp(i+1,j)
                 zone%species(n)%drifts%uBp(i,j)=zone%species(n)%drifts%uBp(i+1,j)
              else
                 if(zone%masks%chi2(i-1,j).eq.0) then
                    zone%species(n)%drifts%uEp(i,j)=zone%species(n)%drifts%uEp(i-1,j)
                    zone%species(n)%drifts%uBp(i,j)=zone%species(n)%drifts%uBp(i-1,j)
                 else
                    zone%species(n)%drifts%uEp(i,j)=0.d0
                    zone%species(n)%drifts%uBp(i,j)=0.d0
                 end if
              end if
           end if
        end do
     end do
  end do
end subroutine correct_drifts_in_pen_cells
   
