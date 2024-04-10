subroutine find_ref_point_for_NOK_points(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4 :: i,j,k,Nx,Nz
  integer*4 :: checksum
! points are OK if surrounded by plasma points
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  k=zone%number
  allocate(zone%DriftsExtra%ref_points(Nx,Nz,3,8))
  allocate(zone%DriftsExtra%Nref_points(Nx,Nz))
  zone%DriftsExtra%ref_points=0
  zone%DriftsExtra%Nref_points=0
  do i=1,Nx
     do j=1,Nz
        if(.not.zone%DriftsExtra%OK_points(i,j)) then
           ! look around for OK point
           ! best guess direction first
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i+1,j).eq.1)) then
              if(zone%DriftsExtra%OK_points(i-1,j)) then
                   zone%DriftsExtra%Nref_points(i,j)=1
                   zone%DriftsExtra%ref_points(i,j,1,1)=i-1
                   zone%DriftsExtra%ref_points(i,j,2,1)=j
                   zone%DriftsExtra%ref_points(i,j,3,1)=k
              end if
           end if
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i-1,j).eq.1)) then
              if(zone%DriftsExtra%OK_points(i+1,j)) then
                   zone%DriftsExtra%Nref_points(i,j)=1
                   zone%DriftsExtra%ref_points(i,j,1,1)=i+1
                   zone%DriftsExtra%ref_points(i,j,2,1)=j
                   zone%DriftsExtra%ref_points(i,j,3,1)=k
              end if
           end if
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i,j+1).eq.1)) then
              if(zone%DriftsExtra%OK_points(i,j-1)) then
                   zone%DriftsExtra%Nref_points(i,j)=1
                   zone%DriftsExtra%ref_points(i,j,1,1)=i
                   zone%DriftsExtra%ref_points(i,j,2,1)=j-1
                   zone%DriftsExtra%ref_points(i,j,3,1)=k
              end if
           end if
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i,j-1).eq.1)) then
              if(zone%DriftsExtra%OK_points(i,j+1)) then
                   zone%DriftsExtra%Nref_points(i,j)=1
                   zone%DriftsExtra%ref_points(i,j,1,1)=i
                   zone%DriftsExtra%ref_points(i,j,2,1)=j+1
                   zone%DriftsExtra%ref_points(i,j,3,1)=k
              end if
           end if
        end if
     end do
  end do

  if(zone%neighbors(1).lt.0) then
     do j=1,Nz
        if(zone%DriftsExtra%OK_points(Nx-1,j)) then
           zone%DriftsExtra%Nref_points(Nx,j)=1
           zone%DriftsExtra%ref_points(Nx,j,1,1)=Nx-1
           zone%DriftsExtra%ref_points(Nx,j,2,1)=j
           zone%DriftsExtra%ref_points(Nx,j,3,1)=k
        end if
     end do
  end if
  if(zone%neighbors(2).lt.0) then
     do j=1,Nz
        if(zone%DriftsExtra%OK_points(2,j)) then
           zone%DriftsExtra%Nref_points(1,j)=1
           zone%DriftsExtra%ref_points(1,j,1,1)=2
           zone%DriftsExtra%ref_points(1,j,2,1)=j
           zone%DriftsExtra%ref_points(1,j,3,1)=k
        end if
     end do
  end if
  if(zone%neighbors(3).lt.0) then
     do i=1,Nx
        if(zone%DriftsExtra%OK_points(i,Nz-1)) then
           zone%DriftsExtra%Nref_points(i,Nz)=1
           zone%DriftsExtra%ref_points(i,Nz,1,1)=i
           zone%DriftsExtra%ref_points(i,Nz,2,1)=Nz-1
           zone%DriftsExtra%ref_points(i,Nz,3,1)=k
        end if
     end do
  end if
  if(zone%neighbors(4).lt.0) then
     do i=1,Nx
        if(zone%DriftsExtra%OK_points(i,2)) then
           zone%DriftsExtra%Nref_points(i,1)=1
           zone%DriftsExtra%ref_points(i,1,1,1)=i
           zone%DriftsExtra%ref_points(i,1,2,1)=2
           zone%DriftsExtra%ref_points(i,1,3,1)=k
        end if
     end do
  end if

end subroutine find_ref_point_for_NOK_points
