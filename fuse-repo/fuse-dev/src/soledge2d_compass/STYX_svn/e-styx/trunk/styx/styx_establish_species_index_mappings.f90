subroutine styx_establish_species_index_mappings
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use eirmod_comusr
  implicit none
  integer*4 :: iatm,imol,iion,ispz,ipls

! establish mapping between species number in EIRENE and species in styx
  allocate(map_species_index(nspz))

  do iatm=1,natmi
    ispz=iatm+nsph
    if (iatm <= global_parameters%n_species) then
      map_species_index(ispz)=iatm
    else
      ! this occurs when sputtered species is not followed in Soledge2D
      map_species_index(ispz)=0
    endif
  enddo

! for the moment molecules are attributed to species 1 (same wall recycling)
  do imol=1,nmoli
    ispz=imol+nspa
    map_species_index(ispz)=1
  enddo

  do iion=1,nioni  
    ispz=iion+nspam
    map_species_index(ispz)=1
  enddo

  do ipls=1,nplsi
    ispz=ipls+nspami
    if (ipls <= global_parameters%n_ions) then
      map_species_index(ispz)=global_parameters%ions_list(ipls,1)
    else
      ! this occurs when sputtered species is not followed in Soledge2D
      map_species_index(ispz)=0
    endif
  enddo

end subroutine styx_establish_species_index_mappings

