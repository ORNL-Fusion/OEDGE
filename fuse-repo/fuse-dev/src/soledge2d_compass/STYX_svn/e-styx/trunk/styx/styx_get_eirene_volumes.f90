subroutine styx_get_eirene_volumes
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comprt
  use eirmod_comusr
  use eirmod_comxs
  use eirmod_comprt
  use eirmod_clgin
  use eirmod_cgrid
  use eirmod_cpolyg
  use eirmod_cplot
  use eirmod_clogau
  use eirmod_ctetra
  use eirmod_ctrig
  use eirmod_cgeom
  use eirmod_ccona
  use eirmod_ctrcei
  use styx2eirene
  implicit none
  integer :: IR


  vol_tri_eirene=vol
! compute total volume
  
  do IR=1,Neir_cells
      	!vol_tri_eirene(IR)=vol(IR)
      	voltot_eirene = voltot_eirene + vol_tri_eirene(IR)
  enddo

#ifdef TK3X
! in styx, this is done in check_volume consistency; should be moved here.
  vol_tri_eirene = vol_tri_eirene/1d6
#endif

end subroutine
