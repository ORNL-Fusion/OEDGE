subroutine update_vorticity(zone)
  use MZone
  implicit none
  Type(TZone) :: zone
  zone%electric_fields(1)%phi=zone%electric_fields(2)%phi
  zone%electric_fields(1)%vorticity=zone%electric_fields(2)%vorticity
end subroutine update_vorticity
