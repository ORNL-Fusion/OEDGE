subroutine styx_rescale_neutral_fields
  use all_variables, only : global_parameters,interp_data2
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use eirmod_comprt, only : iunout
  use styx2eirene
  use Meirene_vars
  implicit none
  ! use nspami as dimension, module comusr
  real(dp) :: tot_flx(Nrecyc)
  integer :: itri,iside,ipuff,isurf,istra,itor
  integer :: iatm,imol,iion,imolt,iiont,ipls
  integer :: nstra0,iplsa,iplse

! how can "new" species appear in a given stratum (recycling or recomb for iatm)
! sputtering -> wall_mat
! CX He+D+ -> D + He+ and D is followed
! so necessary to rescale all species in all strata, or use an indicator array 
! in all cases, exept for sputtering, amounts to rescale with the variation of recycling flux

! total outlux (plasma) part/s, to be used to rescale the recycling source
  tot_flx=0._dp

  
! this is for recycling fluxes, summed over ionization stages for each species
  do istra=1,Nrecyc
    ! istra = isp species number
    iplsa=singly_charged_ions(istra)
    iplse=iplsa+global_parameters%element_list(istra)%Z-1
    do isurf=1,Nsou
      itri=recsurf(isurf)%itri
      iside=recsurf(isurf)%iside
      do itor=1,Ntor_cells
        do ipls=iplsa,iplse
          tot_flx(istra)=tot_flx(istra)+abs(Interp_data2%tri_fluxN(ITRI,ISIDE,ipls,itor))
        enddo
      enddo
    enddo
  enddo
 
! scresc(natmi+nmoli+niioni,1:5,1:NSTRATA)
! initialize to 1 and repair flux save

  scresc=1._dp
    
  do istra=1,Nrecyc
    if(recflux_save(istra).ne.0.D0) then
        scresc(1:nspami,1,istra) = tot_flx(istra)/recflux_save(istra)
    else
     write(*,*) ' zero saved flux, no rescaling done for recycling stratum # ',istra
    end if
  enddo

  write(*,*) '------------------------------------------------------------------------------'
  write(*,*) '  refreshing sources ...'
  do istra=1,Nrecyc
    write(*,*) '  flux rescaling (recycling) in short cycle : scresc(istra,1,istra) = ', scresc(istra,1,istra)
  enddo
   write(*,*) '------------------------------------------------------------------------------'

  nstra0=Nrecyc+Nrecomb

  ! this assumes feedback acting on strata 3 only, not general !
  if (eirene_vars%feedback) then
    do ipuff=1,Npuffs
      if (puff_save(ipuff).ne.0.d0) then
        scresc(1:nspami,1,nstra0+ipuff) = Interp_data2%Puff(ipuff)/puff_save(ipuff)       
      elseif (puff_save(ipuff) == 0._dp) then
        write(*,*) ' zeros saved flux for puff #',ipuff,', no rescaling done for puff stratum ... '
      endif
    enddo
  endif
   
! initialization (thanks Patrick)

  atom_density=0._dp
  E0_at=0._dp
  Gamma0X_at=0._dp
  Gamma0Y_at=0._dp
  Gamma0Z_at=0._dp
  mol_density=0._dp
  E0_mol=0._dp
  Gamma0X_mol=0._dp
  Gamma0Y_mol=0._dp
  Gamma0Z_mol=0._dp
  tion_density=0._dp
  E0_tion=0._dp
  Gamma0X_tion=0._dp
  Gamma0Y_tion=0._dp
  Gamma0Z_tion=0._dp

     
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(itri,iatm,imol,iion,imolt,iiont,ipuff,istra)
!$OMP DO
  do itri=1,Neir_cells
    do iatm=1,natmi
        do istra=1,nstra0        
          atom_density(iatm,itri)    = atom_density(iatm,itri)+scresc(iatm,1,istra)*eov(istra)%pdena(iatm,itri)
          E0_at(iatm,itri)           = E0_at(iatm,itri)       +scresc(iatm,2,istra)*eov(istra)%edena(iatm,itri)
          Gamma0X_at(iatm,itri)      = Gamma0X_at(iatm,itri)  +scresc(iatm,3,istra)*eov(istra)%vxdena(iatm,itri)
          Gamma0Y_at(iatm,itri)      = Gamma0Y_at(iatm,itri)  +scresc(iatm,4,istra)*eov(istra)%vydena(iatm,itri)
          Gamma0Z_at(iatm,itri)      = Gamma0Z_at(iatm,itri)  +scresc(iatm,5,istra)*eov(istra)%vzdena(iatm,itri)
        enddo
        do ipuff=nstra0+1,nstra0+Npuffs          
          atom_density(iatm,itri)= atom_density(iatm,itri)+scresc(iatm,1,ipuff)*eov(ipuff)%pdena(iatm,itri)
          E0_at(iatm,itri)       = E0_at(iatm,itri)       +scresc(iatm,2,ipuff)*eov(ipuff)%edena(iatm,itri)
          Gamma0X_at(iatm,itri)  = Gamma0X_at(iatm,itri)  +scresc(iatm,3,ipuff)*eov(ipuff)%vxdena(iatm,itri)
          Gamma0Y_at(iatm,itri)  = Gamma0Y_at(iatm,itri)  +scresc(iatm,4,ipuff)*eov(ipuff)%vxdena(iatm,itri)
          Gamma0Z_at(iatm,itri)  = Gamma0Z_at(iatm,itri)  +scresc(iatm,5,ipuff)*eov(ipuff)%vxdena(iatm,itri)
        enddo
      enddo
    
      do imol=1,nmoli
        imolt=imol+natmi
        do istra=1,nstra0
          mol_density(imol,itri)      = mol_density(imol,itri)+scresc(imolt,1,istra)*eov(istra)%pdenm(imol,itri)
          E0_mol(imol,itri)           = E0_mol(imol,itri)     +scresc(imolt,2,istra)*eov(istra)%edenm(imol,itri)
          Gamma0X_mol(imol,itri)      = Gamma0X_mol(imol,itri)+scresc(imolt,3,istra)*eov(istra)%vxdenm(imol,itri)
          Gamma0Y_mol(imol,itri)      = Gamma0Y_mol(imol,itri)+scresc(imolt,4,istra)*eov(istra)%vydenm(imol,itri)
          Gamma0Z_mol(imol,itri)      = Gamma0Z_mol(imol,itri)+scresc(imolt,5,istra)*eov(istra)%vzdenm(imol,itri)
        enddo
        do ipuff=nstra0+1,nstra0+Npuffs
          mol_density(imol,itri)  = mol_density(imol,itri)+scresc(imolt,1,ipuff)*eov(ipuff)%pdenm(imol,itri)
          E0_mol(imol,itri)       = E0_mol(imol,itri)     +scresc(imolt,2,ipuff)*eov(ipuff)%edenm(imol,itri)
          Gamma0X_mol(imol,itri)  = Gamma0X_mol(imol,itri)+scresc(imolt,3,ipuff)*eov(ipuff)%vxdenm(imol,itri)
          Gamma0Y_mol(imol,itri)  = Gamma0Y_mol(imol,itri)+scresc(imolt,4,ipuff)*eov(ipuff)%vxdenm(imol,itri)
          Gamma0Z_mol(imol,itri)  = Gamma0Z_mol(imol,itri)+scresc(imolt,5,ipuff)*eov(ipuff)%vxdenm(imol,itri)
        enddo
      enddo
    
      do iion=1,nioni
        iiont = iion + natmi +nmoli
        do istra=1,nstra0
          tion_density(iion,itri)     = tion_density(iion,itri)+scresc(iiont,1,istra)*eov(istra)%pdeni(iion,itri)
          E0_tion(iion,itri)          = E0_tion(iion,itri)     +scresc(iiont,2,istra)*eov(istra)%edeni(iion,itri) 
          Gamma0X_tion(iion,itri)     = Gamma0X_tion(iion,itri)+scresc(iiont,3,istra)*eov(istra)%vxdeni(iion,itri)
          Gamma0Y_tion(iion,itri)     = Gamma0Y_tion(iion,itri)+scresc(iiont,4,istra)*eov(istra)%vydeni(iion,itri)
          Gamma0Z_tion(iion,itri)     = Gamma0Z_tion(iion,itri)+scresc(iiont,5,istra)*eov(istra)%vzdeni(iion,itri)
        enddo
        do ipuff=1,Npuffs
          tion_density(iion,itri) =tion_density(iion,itri)  +scresc(iiont,1,ipuff)*eov(ipuff)%pdeni(iion,itri)
          E0_tion(iion,itri)      = E0_tion(iion,itri)      +scresc(iiont,2,ipuff)*eov(ipuff)%edeni(iion,itri)
          Gamma0X_tion(iion,itri) = Gamma0X_tion(iion,itri) +scresc(iiont,3,ipuff)*eov(ipuff)%vxdeni(iion,itri)
          Gamma0Y_tion(iion,itri) = Gamma0Y_tion(iion,itri) +scresc(iiont,4,ipuff)*eov(ipuff)%vxdeni(iion,itri)
          Gamma0Z_tion(iion,itri) = Gamma0Z_tion(iion,itri) +scresc(iiont,5,ipuff)*eov(ipuff)%vxdeni(iion,itri)
        enddo
      enddo
    
    enddo
!$OMP END DO

!$OMP END PARALLEL
    recflux_save=tot_flx
    puff_save=Interp_Data2%Puff

end subroutine


