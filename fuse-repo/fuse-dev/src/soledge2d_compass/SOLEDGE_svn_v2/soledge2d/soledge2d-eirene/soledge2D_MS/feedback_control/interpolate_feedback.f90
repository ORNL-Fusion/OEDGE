function interpolate_feedback(x,type) result(D)
  use MradialFeedback
  implicit none
  real*8, intent(in) :: x
  integer*4, intent(in) :: type
  real*8 :: D
  integer*4 :: index,index2,k
  real*8 :: weight
  integer*4 :: Nxtot
  real*8,allocatable :: Ds(:)
  Nxtot=radialFeedbackData%Nxtot
  allocate(Ds(1:Nxtot))
  if(x.lt.radialFeedbackData%x(1)) then
     index=1
     index2=1
     weight=1.
  else
     if(x.gt.radialFeedbackData%x(Nxtot)) then
        index=Nxtot
        index2=Nxtot
        weight=1.
     else
        k=2
        do while(x.gt.radialFeedbackData%x(k))
           k=k+1
        end do
        index2=k
        index=k-1
        weight=1.D0-(x-radialFeedbackData%x(index))/(radialFeedbackData%x(index2)-radialFeedbackData%x(index))
     end if
  end if
  
  select case(type)
     case (1) !D
        Ds=radialFeedbackData%D
!!$        Ds(1)=(radialFeedbackData%D(1)+radialFeedbackData%D(2))/2.D0
!!$        Ds(2:Nxtot-1)=(radialFeedbackData%D(1:Nxtot-2)+radialFeedbackData%D(2:Nxtot-1)+radialFeedbackData%D(3:Nxtot))/3.D0
!!$        Ds(Nxtot)=(radialFeedbackData%D(Nxtot-1)+radialFeedbackData%D(Nxtot))/2.D0
        D=weight*Ds(index)+(1.D0-weight)*Ds(index2)
     case (2) !chie
        Ds=radialFeedbackData%chie
!!$        Ds(1)=(radialFeedbackData%chie(1)+radialFeedbackData%chie(2))/2.D0
!!$        Ds(2:Nxtot-1)=(radialFeedbackData%chie(1:Nxtot-2)+radialFeedbackData%chie(2:Nxtot-1)+radialFeedbackData%chie(3:Nxtot))/3.D0
!!$        Ds(Nxtot)=(radialFeedbackData%chie(Nxtot-1)+radialFeedbackData%chie(Nxtot))/2.D0
        D=weight*Ds(index)+(1.D0-weight)*Ds(index2)
     case (3) !chii
        Ds=radialFeedbackData%chii
!!$        Ds(1)=(radialFeedbackData%chii(1)+radialFeedbackData%chii(2))/2.D0
!!$        Ds(2:Nxtot-1)=(radialFeedbackData%chii(1:Nxtot-2)+radialFeedbackData%chii(2:Nxtot-1)+radialFeedbackData%chii(3:Nxtot))/3.D0
!!$        Ds(Nxtot)=(radialFeedbackData%chii(Nxtot-1)+radialFeedbackData%chii(Nxtot))/2.D0
        D=weight*Ds(index)+(1.D0-weight)*Ds(index2)
     case (4) !density
        D=weight*radialFeedbackData%input_n(index)+(1.D0-weight)*radialFeedbackData%input_n(index2)
     case (5) !Te
        D=weight*radialFeedbackData%input_Te(index)+(1.D0-weight)*radialFeedbackData%input_Te(index2)
     case (6) !Ti
        D=weight*radialFeedbackData%input_Ti(index)+(1.D0-weight)*radialFeedbackData%input_Ti(index2)
  end select
end function interpolate_feedback
