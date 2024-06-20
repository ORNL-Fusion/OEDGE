subroutine detect_problem_with_ref_points(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: i,j,k,Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  k=zone%number
  do i=1,Nx
     do j=1,Nz
        if(.not.zone%DriftsExtra%OK_points(i,j)) then
           if((zone%masks%chi2(i,j).eq.0).and.(zone%driftsExtra%Nref_points(i,j).eq.0)) then
              write(*,*) 'problem to find ref point',i,j,k
              write(*,*) 'details',zone%DriftsExtra%OK_points(i,j),zone%masks%chi2(i,j),zone%driftsExtra%Nref_points(i,j)
              stop
           end if
        end if
     end do
  end do
end subroutine detect_problem_with_ref_points
