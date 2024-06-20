subroutine find_pts_around_penwall()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: i,j,k,Nx,Nz
  integer*4 :: North,South,East,West
  do k=1,global_parameters%N_zones
     Nx=Zones(k)%mesh%Nx
     Nz=Zones(k)%mesh%Nz
     North=Zones(k)%Neighbors(1)
     South=Zones(k)%Neighbors(2)
     East=Zones(k)%Neighbors(3)
     West=Zones(k)%Neighbors(4)
     Zones(k)%masks%npts_around_penwall=0
     Zones(k)%masks%pts_around_penwall=0
     do i=1,Nx
        do j=1,Nz
           if(Zones(k)%masks%chi(i,j).eq.1) then
              !East point
              if(j.lt.Nz) then
                 if(Zones(k)%masks%chi(i,j+1).eq.0) then
                    Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=i
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=j+1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=k
                 end if
              else
                 if(East.gt.0) then
                    if(Zones(East)%masks%chi(i,1).eq.0) then
                       Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=i
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=1
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=East
                    end if
                 end if
              end if

              !West point
              if(j.gt.1) then
                 if(Zones(k)%masks%chi(i,j-1).eq.0) then
                    Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=i
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=j-1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=k
                 end if
              else
                 if(West.gt.0) then
                    if(Zones(West)%masks%chi(i,Zones(West)%mesh%Nz).eq.0) then
                       Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=i
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=Zones(West)%mesh%Nz
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=West
                    end if
                 end if
              end if

              !North point
              if(i.lt.Nx) then
                 if(Zones(k)%masks%chi(i+1,j).eq.0) then
                    Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=i+1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=j
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=k
                 end if
              else
                 if(North.gt.0) then
                    if(Zones(North)%masks%chi(1,j).eq.0) then
                       Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=1
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=j
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=North
                    end if
                 end if
              end if

              !South point
              if(i.gt.1) then
                 if(Zones(k)%masks%chi(i-1,j).eq.0) then
                    Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=i-1
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=j
                    Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=k
                 end if
              else
                 if(South.gt.0) then
                    if(Zones(South)%masks%chi(Zones(South)%mesh%Nx,j).eq.0) then
                       Zones(k)%masks%npts_around_penwall(i,j)=Zones(k)%masks%npts_around_penwall(i,j)+1
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),1)=Zones(South)%mesh%Nx
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),2)=j
                       Zones(k)%masks%pts_around_penwall(i,j,Zones(k)%masks%npts_around_penwall(i,j),3)=South
                    end if
                 end if
              end if
           end if !chi=1
        end do !j
     end do !i
  end do !k
end subroutine find_pts_around_penwall
