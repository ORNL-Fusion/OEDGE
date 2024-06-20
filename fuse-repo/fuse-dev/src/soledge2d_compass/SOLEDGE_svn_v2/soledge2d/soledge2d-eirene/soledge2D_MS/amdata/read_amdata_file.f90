subroutine read_amdata_file(filename,polynom_buffer,Z,degree)
  use all_variables, only : global_parameters
  use Mphysics
  implicit none
  character(50), intent(in) :: filename
  integer*4,intent(in) :: Z
  integer*4,intent(in) :: degree
  type(TAm_polynom),intent(inout) :: polynom_buffer(1:Z) 
  integer :: i,j,iz,izmin,izmax,nl,nr,nf,k
  real*8 :: r
  character(2) :: ncolc,nrc
  character(11) :: format,format2
  real*8 :: Nemin,Nemax
  integer*4 :: nne,nte
  integer*4,parameter :: ncol=5
  open(unit=666,file=filename,status='old')

  do i=1,9
    read(666,*)
  enddo
  read(666,'(14x,i2,8x,i3)') izmin,izmax
  read(666,*)
  read(666,*) 
  read(666,*)
  read(666,'(10x,I3)') nne
  read(666,'(10x,e19.12)') Nemin
  read(666,'(10x,e19.12)') Nemax
  polynom_buffer(1:Z)%Ne_min=Nemin
  polynom_buffer(1:Z)%Ne_max=Nemax
  read(666,*)
 
  r=nne/ncol  
  nl=floor(r)

  do i=1,nl
    read(666,*)
  enddo

  nr=nne-nl*ncol

  if (nr>0) then
    read(666,*)
  endif

  ! now read temperature data

  do iz=izmin,izmax

    read(666,*)
    read(666,*)
    read(666,*)
    
    read(666,'(10x,i3)') nTe
    read(666,'(10x,e19.12)') polynom_buffer(iz)%Te_min
    read(666,'(10x,e19.12)') polynom_buffer(iz)%Te_max

    read(666,*)

    r=nTe/ncol  
    nl=floor(r)

    do i=1,nl
      read(666,*)
    enddo

    nr=nTe-nl*ncol

    if (nr>0) then
      read(666,*) 
    endif

    read(666,*)

    ! number of fitting coefficients

    nf=(degree+1)*degree/2
   
    r=nf/ncol  
    nl=floor(r)

    write(ncolc,'(I2)') ncol
    format='('//trim(adjustl(ncolc))//'(e19.12))'

    do i=1,nl
      read(666,format) (polynom_buffer(iz)%coefficients(j+(i-1)*ncol),j=1,ncol)
    enddo

    nr=nf-nl*ncol

    if (nr>0) then
      write(nrc,'(i1)') nr
      format2='('//trim(adjustl(nrc))//'(e19.12))'  
      read(666,format2) (polynom_buffer(iz)%coefficients(j+nl*ncol),j=1,nr)
    endif

    read(666,*)
    read(666,*)
    read(666,*)
    read(666,*)

  enddo
  
  close(666)  


end subroutine read_amdata_file
