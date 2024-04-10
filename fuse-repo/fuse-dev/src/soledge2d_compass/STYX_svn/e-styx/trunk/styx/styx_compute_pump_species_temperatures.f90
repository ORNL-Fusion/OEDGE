subroutine styx_compute_pump_species_temperature()
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_clgin
  use eirmod_ccona
  use styx2eirene
  use all_variables, only : global_parameters
  use Mphysics
  use eirmod_comusr
  use eirmod_ccona
  use eirmod_cpes

  implicit none

  integer :: isp,ispz,ipfc,ipump,itri
  integer :: iatm,imol,iion

  !default to 0.
  do ispz=1,nspz
     do ipump=1,n_pumps
        pumps(ipump)%species_temperature(ispz)=0.d0
     end do
  end do

  ! set recycling fluxs (in the future may evolve from call to call)
  do iatm=1,natmi
     ispz=iatm+nsph
     isp=map_species_index(ispz)
     if (isp >0) then
        do ipump=1,n_pumps
           pumps(ipump)%species_temperature(ispz)=0.d0
           do itri=1,pumps(ipump)%Ntriangles
              pumps(ipump)%species_temperature(ispz)=pumps(ipump)%species_temperature(ispz)&
                   +T0_at(isp,pumps(ipump)%triangle_index(itri))*pumps(ipump)%triangle_weight(itri)
           end do
        enddo
     else
        do ipump=1,n_pumps
           pumps(ipump)%species_temperature(ispz)=0.d0
        end do
     endif
  enddo

  do imol=1,nmoli
     ispz=imol+nspa
     isp=map_species_index(ispz)
     if (isp >0) then
        do ipump=1,n_pumps
           pumps(ipump)%species_temperature(ispz)=0.d0
           do itri=1,pumps(ipump)%Ntriangles
              pumps(ipump)%species_temperature(ispz)=pumps(ipump)%species_temperature(ispz)&
                   +T0_mol(isp,pumps(ipump)%triangle_index(itri))*pumps(ipump)%triangle_weight(itri)
           end do
        enddo
     else
        do ipump=1,n_pumps
           pumps(ipump)%species_temperature(ispz)=0.d0
        end do
     endif
  enddo

  do iion=1,nioni
     ispz=iion+nspam
     isp=map_species_index(ispz)
     if (isp >0) then
        do ipump=1,n_pumps
           pumps(ipump)%species_temperature(ispz)=0.d0
           do itri=1,pumps(ipump)%Ntriangles
              pumps(ipump)%species_temperature(ispz)=pumps(ipump)%species_temperature(ispz)&
                   +T0_tion(isp,pumps(ipump)%triangle_index(itri))*pumps(ipump)%triangle_weight(itri)
           end do
        enddo
     else
        do ipump=1,n_pumps
           pumps(ipump)%species_temperature(ispz)=0.d0
        end do
     endif
  enddo

  ! ev to kelvin
  do ipump=1,n_pumps
     do ispz=1,nspz
        pumps(ipump)%species_temperature(ispz)=pumps(ipump)%species_temperature(ispz)&
             *eV/kb
     end do
  enddo

end subroutine styx_compute_pump_species_temperature
