subroutine concat_chi(n_mega,bchi,bNx,bNz)
  use all_variables, only : zones, megazones
  implicit none
  integer*4,intent(in) :: n_mega
  integer*4,intent(in) :: bNx
  integer*4,intent(in) :: bNz
  integer*4,intent(inout) :: bchi(1:bNx,1:bNz)
  integer*4 :: ideb,size
  integer*4 :: i,j,k
  integer*4 :: zn
  size=megazones(n_mega)%size
  ideb=0
  do k=1,size
     zn=megazones(n_mega)%zone_number(k)
     do i=1,bNx
        do j=1,zones(zn)%mesh%Nz
           bchi(i,ideb+j)=zones(zn)%masks%chi(i,j)
        end do
     end do
     ideb=ideb+zones(zn)%mesh%Nz
  end do
end subroutine concat_chi
