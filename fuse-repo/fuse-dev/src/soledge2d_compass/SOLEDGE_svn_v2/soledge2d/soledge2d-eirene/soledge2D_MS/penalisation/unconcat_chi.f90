subroutine unconcat_chi(n_mega,bchis,bNx,bNz)
  use all_variables, only : zones, megazones
  implicit none
  integer*4,intent(in) :: n_mega
  integer*4,intent(in) :: bNx
  integer*4,intent(in) :: bNz
  integer*4,intent(in) :: bchis(1:4,1:bNx,1:bNz)
  integer*4 :: ideb,size
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  ideb=0
  size=megazones(n_mega)%size
  do k=1,size
     Nx=zones(megazones(n_mega)%zone_number(k))%mesh%Nx
     Nz=zones(megazones(n_mega)%zone_number(k))%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           zones(megazones(n_mega)%zone_number(k))%masks%chi1(i,j)=bchis(1,i,ideb+j)
           zones(megazones(n_mega)%zone_number(k))%masks%chi2(i,j)=bchis(2,i,ideb+j)
           zones(megazones(n_mega)%zone_number(k))%masks%chi3(i,j)=bchis(3,i,ideb+j)
           zones(megazones(n_mega)%zone_number(k))%masks%chi4(i,j)=bchis(4,i,ideb+j)
        end do
     end do
     ideb=ideb+Nz
  end do
end subroutine unconcat_chi
