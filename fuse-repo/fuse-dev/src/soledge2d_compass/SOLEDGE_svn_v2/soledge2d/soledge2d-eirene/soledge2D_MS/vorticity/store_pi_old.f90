subroutine store_pi_old(zone)
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  zone%electric_fields(1)%pi_old=zone%species(1)%var(1)%density&
       *zone%species(1)%var(1)%temperature
  zone%electric_fields(1)%cornersPi_old=zone%species(1)%corners%density&
       *zone%species(1)%corners%temperature
end subroutine store_pi_old
