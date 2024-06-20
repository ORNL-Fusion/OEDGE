subroutine styx_test_preav_sources(coupling_mode)
! purpose : compare sources calculated by EIRENE to those recalculated in pre-averaged mode
  use all_variables, only : global_parameters
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use eirmod_ccona
  use styx2eirene
  
  implicit none
  
  integer, intent(in) :: coupling_mode
  integer :: itri,ipls
  real(dp) :: dum



  dum=0._dp

  if (coupling_mode == 0) call styx_preav_sources


  open(unit=666,file='preav_sources_test_atoms',status='replace')

  do itri=1,Neir_cells
    write(666,'(4es12.4)') Sn_at(itri,0),dum,SE_at(itri,0)
  enddo

  do ipls=1,nplsi 
    do itri=1,Neir_cells
      write(666,'(4es12.4)') Sn_at(itri,ipls),Sm_at(itri,ipls),SE_at(itri,ipls)
    enddo
  enddo

  close(666)
  
  open(unit=666,file='preav_sources_test_molec',status='replace')

  do itri=1,Neir_cells
    write(666,'(4es12.4)') Sn_mol(itri,0),dum,SE_mol(itri,0)
  enddo


  do ipls=1,nplsi
    do itri=1,Neir_cells
      write(666,'(4es12.4)') Sn_mol(itri,ipls),Sm_mol(itri,ipls),SE_mol(itri,ipls)
    enddo
  enddo

  close(666)

  open(unit=666,file='preav_sources_test_tions',status='replace')
 
  do itri=1,Neir_cells
    write(666,'(4es12.4)') Sn_tion(itri,0),dum,SE_tion(itri,0)
  enddo
 
  do ipls=1,nplsi
    do itri=1,Neir_cells
      write(666,'(4es12.4)') Sn_tion(itri,ipls),Sm_tion(itri,ipls),SE_tion(itri,ipls)
    enddo
  enddo
 
  close(666)

  
  open(unit=666,file='direct_sources_test_atoms',status='replace')

  do itri=1,Neir_cells
    write(666,'(3es12.4)') eov(0)%pael(itri),dum,eov(0)%eael(itri)
  enddo

  do ipls=1,nplsi
    do itri=1,Neir_cells
      write(666,'(3es12.4)') eov(0)%papl(ipls,itri),eov(0)%mapl(ipls,itri)/&
                       (RMASSP(ipls)*amuakg)*sign(1._dp,vpar_tri(itri,ipls)), &
                           eov(0)%eaplr(ipls,itri)
    enddo
  enddo
 
  close(666)

  open(unit=666,file='direct_sources_test_molec',status='replace')

  do itri=1,Neir_cells
    write(666,'(3es12.4)') eov(0)%pmel(itri),dum,eov(0)%emel(itri)
  enddo

  do ipls=1,nplsi
    do itri=1,Neir_cells
      write(666,'(3es12.4)') eov(0)%pmpl(ipls,itri),eov(0)%mmpl(ipls,itri)/&
                       (RMASSP(ipls)*amuakg)*sign(1._dp,vpar_tri(itri,ipls)), &
                           eov(0)%emplr(ipls,itri)
    enddo
  enddo

  close(666)

  open(unit=666,file='direct_sources_test_tions',status='replace')
 
  do itri=1,Neir_cells
    write(666,'(3es12.4)') eov(0)%piel(itri),dum,eov(0)%eiel(itri)
  enddo

  do ipls=1,nplsi
    do itri=1,Neir_cells
      write(666,'(3es12.4)') eov(0)%pipl(ipls,itri),eov(0)%mipl(ipls,itri)/&
                       (RMASSP(ipls)*amuakg)*sign(1._dp,vpar_tri(itri,ipls)), &
                           eov(0)%eiplr(ipls,itri)
    enddo
  enddo

  close(666)


end subroutine styx_test_preav_sources
