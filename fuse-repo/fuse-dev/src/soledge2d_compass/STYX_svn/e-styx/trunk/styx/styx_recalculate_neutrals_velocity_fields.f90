subroutine styx_recalculate_neutrals_velocity_fields
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_ccona
  use eirmod_comusr
  use styx2eirene
  implicit none
  integer :: iatm,imol,iion

  V0X_at=0._dp
  V0Y_at=0._dp
  V0Z_at=0._dp

  V0X_mol=0._dp
  V0Y_mol=0._dp
  V0Z_mol=0._dp

  V0X_tion=0._dp
  V0Y_tion=0._dp
  V0Z_tion=0._dp

  do iatm=1,natmi
    where (atom_density(iatm,:) > 0._dp)
      E0_at(iatm,:)=E0_at(iatm,:)/atom_density(iatm,:)
      V0X_at(iatm,:)=Gamma0X_at(iatm,:)/(atom_density(iatm,:)*rmassa(iatm)*amuakg)
      V0Y_at(iatm,:)=Gamma0Y_at(iatm,:)/(atom_density(iatm,:)*rmassa(iatm)*amuakg)
      V0Z_at(iatm,:)=Gamma0Z_at(iatm,:)/(atom_density(iatm,:)*rmassa(iatm)*amuakg)
    elsewhere
      E0_at(iatm,:)=0._dp
      V0X_at(iatm,:)=0._dp
      V0Y_at(iatm,:)=0._dp
      V0Z_at(iatm,:)=0._dp
    end where
  enddo

  do imol=1,nmoli
    where (mol_density(imol,:) > 0._dp)
      E0_mol(imol,:)=E0_mol(imol,:)/mol_density(imol,:)
      V0X_mol(imol,:)=Gamma0X_mol(imol,:)/(mol_density(imol,:)*rmassm(imol)*amuakg)
      V0Y_mol(imol,:)=Gamma0Y_mol(imol,:)/(mol_density(imol,:)*rmassm(imol)*amuakg)
      V0Z_mol(imol,:)=Gamma0Z_mol(imol,:)/(mol_density(imol,:)*rmassm(imol)*amuakg)
    elsewhere
     E0_mol(imol,:)=0._dp
     V0X_mol(imol,:)=0._dp
     V0Y_mol(imol,:)=0._dp
     V0Z_mol(imol,:)=0._dp
    end where
  enddo

  do iion=1,nioni
    where (tion_density(iion,:) > 0._dp)
      E0_tion(iion,:)=E0_tion(iion,:)/tion_density(iion,:)
      V0X_tion(iion,:)=Gamma0X_tion(iion,:)/(tion_density(iion,:)*rmassi(iion)*amuakg)
      V0Y_tion(iion,:)=Gamma0Y_tion(iion,:)/(tion_density(iion,:)*rmassi(iion)*amuakg)
      V0Z_tion(iion,:)=Gamma0Z_tion(iion,:)/(tion_density(iion,:)*rmassi(iion)*amuakg)
    elsewhere
      E0_tion(iion,:)=0._dp
      V0X_tion(iion,:)=0._dp
      V0Y_tion(iion,:)=0._dp
      V0Z_tion(iion,:)=0._dp
    end where
  enddo

! now:
! atom_density in m-3
! E0 is in J (quadratic mean energy per neutral)
! V0X,Y,Z in m.s-1 (neutral average velocity)

  do iatm=1,natmi
     V0par_at(iatm,1:Neir_cells)=V0X_at(iatm,1:Neir_cells)*BX_tri(1:Neir_cells)+ &
                                     V0Y_at(iatm,1:Neir_cells)*BY_tri(1:Neir_cells)+ &
                                     V0Z_at(iatm,1:Neir_cells)*BZ_tri(1:Neir_cells)
  enddo

  do imol=1,nmoli
    V0par_mol(imol,1:Neir_cells)=V0X_mol(imol,1:Neir_cells)*BX_tri(1:Neir_cells)+ &
                                      V0Y_mol(imol,1:Neir_cells)*BY_tri(1:Neir_cells)+ &
                                      V0Z_mol(imol,1:Neir_cells)*BZ_tri(1:Neir_cells)
  enddo

  do iion=1,nioni
    V0par_tion(iion,1:Neir_cells)=V0X_tion(iion,1:Neir_cells)*BX_tri(1:Neir_cells)+ &
                                       V0Y_tion(iion,1:Neir_cells)*BY_tri(1:Neir_cells)+ &
                                       V0Z_tion(iion,1:Neir_cells)*BZ_tri(1:Neir_cells)
  enddo
end subroutine
