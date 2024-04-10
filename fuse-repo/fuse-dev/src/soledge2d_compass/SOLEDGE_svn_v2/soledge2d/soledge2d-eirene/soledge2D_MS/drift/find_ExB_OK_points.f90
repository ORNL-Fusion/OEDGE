subroutine find_ExB_OK_points(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: i,j,Nx,Nz
  integer*4 :: checksum
! points are OK if surrounded by plasma points
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(zone%DriftsExtra%OK_points(0:Nx+1,0:Nz+1))
  zone%DriftsExtra%OK_points=.false.
  do i=1,Nx
     do j=1,Nz
        checksum=zone%masks%chi2(i,j)+&
             zone%masks%chi2(i,j+1)+&
             zone%masks%chi2(i,j-1)+&
             zone%masks%chi2(i-1,j)+&
             zone%masks%chi2(i+1,j)+&
             zone%masks%chi2(i-1,j+1)+&
             zone%masks%chi2(i-1,j-1)+&
             zone%masks%chi2(i+1,j+1)+&
             zone%masks%chi2(i+1,j-1)
        if(checksum.eq.0) then
           zone%DriftsExtra%OK_points(i,j)=.true.
        end if
     end do
  end do
  if(zone%Neighbors(1).lt.0) then
     zone%DriftsExtra%OK_points(Nx,:)=.false.
  end if
  if(zone%Neighbors(2).lt.0) then
     zone%DriftsExtra%OK_points(1,:)=.false.
  end if
  if(zone%Neighbors(3).lt.0) then
     zone%DriftsExtra%OK_points(:,Nz)=.false.
  end if
  if(zone%Neighbors(4).lt.0) then
     zone%DriftsExtra%OK_points(:,1)=.false.
  end if
  !the extreme corner bidouille
  !NW
  if(zones(zone%neighbors(1))%neighbors(4).lt.0) then
     zone%DriftsExtra%OK_points(Nx,1)=.false.
  end if
  if(zones(zone%neighbors(4))%neighbors(1).lt.0) then
     zone%DriftsExtra%OK_points(Nx,1)=.false.
  end if
  !SW
  if(zones(zone%neighbors(2))%neighbors(4).lt.0) then
     zone%DriftsExtra%OK_points(1,1)=.false.
  end if
  if(zones(zone%neighbors(4))%neighbors(2).lt.0) then
     zone%DriftsExtra%OK_points(1,1)=.false.
  end if
  !NE
  if(zones(zone%neighbors(1))%neighbors(3).lt.0) then
     zone%DriftsExtra%OK_points(Nx,Nz)=.false.
  end if
  if(zones(zone%neighbors(3))%neighbors(1).lt.0) then
     zone%DriftsExtra%OK_points(Nx,Nz)=.false.
  end if
  !SE
  if(zones(zone%neighbors(2))%neighbors(3).lt.0) then
     zone%DriftsExtra%OK_points(1,Nz)=.false.
  end if
  if(zones(zone%neighbors(3))%neighbors(2).lt.0) then
     zone%DriftsExtra%OK_points(1,Nz)=.false.
  end if
end subroutine find_ExB_OK_points
