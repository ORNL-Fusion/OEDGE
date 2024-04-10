subroutine read_feedback_data()
  use all_variables, only : zones
  use MradialFeedback
  use Mreaders
  implicit none
  integer*4,allocatable :: bufferI(:)
  integer*4 :: i,n,k
  integer*4 :: Nxtot,Nx
  character(128) :: String
  integer :: finput
  real*8 :: dum1, dum2, dum3
  real*8,allocatable :: xv(:)
  integer*4 :: Nxshift
  ! read profile location
  finput=10
  open(unit=finput,file='radialFeedback.txt',status='unknown')
  call skip_line(finput,9)
  read(finput,1) radialFeedbackData%Nzones
  allocate(radialFeedbackData%zone_profile(1:radialFeedbackData%Nzones))
  allocate(radialFeedbackData%set(1:radialFeedbackData%Nzones))
  allocate(bufferI(1:radialFeedbackData%Nzones))
  call skip_line(finput,2)
  read(finput,2) String
  call parse_line_integer(String,radialFeedbackData%Nzones,bufferI)
  do i=1,radialFeedbackData%Nzones
     radialFeedbackData%zone_profile(i)=bufferI(i)
  end do
  deallocate(bufferI)
  call skip_line(finput,2)
  read(finput,4) radialFeedbackData%Nz
  call skip_line(finput,6)
  read(finput,5) radialFeedbackData%Dmin
  call skip_line(finput,2)
  read(finput,5) radialFeedbackData%Dmax
  call skip_line(finput,2)
  read(finput,5) radialFeedbackData%keep
  call skip_line(finput,2)
  read(finput,5) radialFeedbackData%Gain
  call skip_line(finput,2)
  read(finput,6) radialFeedbackData%GainG
  close(finput)

  Nxtot=0
  do n=1,radialFeedbackData%Nzones
     Nx=zones(radialFeedbackData%zone_profile(n))%mesh%Nx
     Nxtot=Nxtot+Nx
     allocate(radialFeedbackData%set(n)%fluxN(1:Nx))
     allocate(radialFeedbackData%set(n)%fluxTe(1:Nx))
     allocate(radialFeedbackData%set(n)%fluxTi(1:Nx))
  end do
  allocate(radialFeedbackData%input_gradN(1:Nxtot))
  allocate(radialFeedbackData%input_gradTe(1:Nxtot))
  allocate(radialFeedbackData%input_gradTi(1:Nxtot))
  allocate(radialFeedbackData%input_N(1:Nxtot))
  allocate(radialFeedbackData%input_Te(1:Nxtot))
  allocate(radialFeedbackData%input_Ti(1:Nxtot))
  allocate(radialFeedbackData%input_D(1:Nxtot))
  allocate(radialFeedbackData%input_Chi(1:Nxtot))
  allocate(radialFeedbackData%x(1:Nxtot))
  allocate(radialFeedbackData%D(1:Nxtot))
  allocate(radialFeedbackData%chie(1:Nxtot))
  allocate(radialFeedbackData%chii(1:Nxtot))
  ! read gradient profiles
  open(unit=10,file='radialFeedbackProfiles.txt',status='unknown')
  do i=1,Nxtot
     read(10,*) radialFeedbackData%input_gradN(i), radialFeedbackData%input_gradTe(i)&
          , radialFeedbackData%input_gradTi(i), radialFeedbackData%input_N(i), &
          radialFeedbackData%input_Te(i), radialFeedbackData%input_Ti(i),&
          radialFeedbackData%input_D(i), radialFeedbackData%input_Chi(i)
  end do
  close(10)

  radialFeedbackData%Nxtot=Nxtot

  allocate(xv(1:radialFeedbackData%Nxtot))
  Nxshift=0
  do n=1,radialFeedbackData%Nzones
     k=radialFeedbackData%zone_profile(n)
     Nx=zones(k)%mesh%Nx
     do i=1,Nx
        xv(i+Nxshift)=zones(k)%mesh%x(i,1)
     end do
     Nxshift=Nxshift+Nx
  end do
  radialFeedbackData%x=xv

  do i=1,Nxtot
     write(120,100) radialFeedbackData%x(i),radialFeedbackData%input_n(i),radialFeedbackData%input_Te(i), radialFeedbackData%input_Ti(i)
  end do
100 format(512es15.7)

1 format(9X,I2)    
2 format(A128)
4 format(5X,I4)
5 format(7X,F6.2)
6 format(8X,F6.2)

end subroutine read_feedback_data
