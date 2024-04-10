subroutine find_ref_point_for_NOK_points_step2(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: i,j,k,Nx,Nz,ii,jj,n
  integer*4 :: checksum
  ! points are OK if surrounded by plasma points
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  k=zone%number
  do i=1,Nx
     do j=1,Nz
        if((.not.zone%DriftsExtra%OK_points(i,j)).and.(zone%DriftsExtra%Nref_points(i,j).eq.0)) then
           do ii=-1,1
              do jj=-1,1
                 if(zone%DriftsExtra%OK_points(i+ii,j+jj)) then
                    zone%DriftsExtra%Nref_points(i,j)=zone%DriftsExtra%Nref_points(i,j)+1
                    zone%DriftsExtra%ref_points(i,j,1,zone%DriftsExtra%Nref_points(i,j))=i+ii
                    zone%DriftsExtra%ref_points(i,j,2,zone%DriftsExtra%Nref_points(i,j))=j+jj
                    zone%DriftsExtra%ref_points(i,j,3,zone%DriftsExtra%Nref_points(i,j))=zone%number
                 end if
              end do
           end do
        end if
     end do
  end do
end subroutine find_ref_point_for_NOK_points_step2
