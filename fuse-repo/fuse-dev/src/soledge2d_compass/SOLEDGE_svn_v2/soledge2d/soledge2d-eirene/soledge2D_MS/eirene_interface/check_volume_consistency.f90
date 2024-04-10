subroutine check_volume_consistency()
  use all_variables, only : zones, global_parameters
  use Meirene_vars
  use styx2eirene, only : voltot_eirene,vol_tri_eirene, NTRI_styx
  use eirmod_comprt, only : iunout
  implicit none
  real*8 :: voltot, vol_recalc_tot,vol_recalc_nonpen,mdiff,mdiff_ed,md
  real*8 :: voltot_check,voltot_cut,vol
  real*8 :: pi
  real*8, parameter :: zer=0.d0
  Type(volume),allocatable :: vol_recalc(:),volumes(:)
  integer*4 :: ml(2),ldiff(2),ldiff_ed(2)
  integer*4 :: k,ix,iz,j,kdiff,kdiff_ed
  integer*4 :: kr,ixr,izr,i
  integer*4 :: Ntritot,itr,Ncut,icell
  integer*4 :: Nx,Nz
  integer*4 :: itri(Ntrimax)
  integer*4, allocatable :: ked(:),ixed(:),ized(:),itrip(:,:)
  character(3) :: Nintc
  character(6) :: format
  character(7) filen
  logical :: is_edge,dir_e

  pi=4.D0*atan(1.D0)
  allocate(vol_recalc(1:global_parameters%N_zones))
  allocate(volumes(1:global_parameters%N_zones))

  call compute_volumes(voltot,volumes,global_parameters%N_zones)

  ! units /cm3 -> /m3

  voltot_eirene = voltot_eirene/1d6
  vol_tri_eirene = vol_tri_eirene/1d6

  ! first total volume (voltot for soledge2D)

  write(iunout,*) 'total simulation volume (m3), soledge = ',voltot
  write(iunout,*) 'total simulation volume (m3), eirene = ',voltot_eirene
  write(iunout,*) 'relative difference (%) = ',abs(voltot-voltot_eirene)/voltot_eirene*100.d0

!!$  ! now per soledge2D cell (volumes)
!!$
!!$  ! first reconstruct volumes calculated by EIRENE for triangles, read soledge2D.triangles mapping
!!$
!!$  open(unit=333,file = 'triangles/soledge2D.triangles_mapping',form = 'formatted')
!!$
!!$  write(Nintc,'(i3)') Ntrimax+3
!!$  format = '('//trim(adjustl(Nintc))//'i6)'
!!$
!!$  Ntritot=0
!!$  vol_recalc_tot=0.d0
!!$  vol_recalc_nonpen =0.d0
!!$
!!$  do k=1,global_parameters%N_zones
!!$     Nx=Zones(k)%mesh%Nx
!!$     Nz=Zones(k)%mesh%Nz
!!$     allocate(vol_recalc(k)%cell(0:Nx+1,0:Nz+1))
!!$     vol_recalc(k)%cell=0.d0
!!$     do ix=1,Nx
!!$        do iz=1,Nz
!!$           read(333,format) kr,ixr,izr,(itri(j),j=1,Ntrimax)
!!$           if (kr /= k .or. ixr /= ix .or. izr /= iz) then
!!$              write(*,*) 'inconsistency while reading soledge2D.triangles_mapping'
!!$              write(*,*) ' ixr = ',ixr,' ix = ',ix
!!$              write(*,*) ' izr = ',izr,' iz = ',iz
!!$              write(*,*) ' kr = ',kr,' k = ',k
!!$              call eirene_exit_own(1)
!!$           endif
!!$           ! sum over triangles corresponding to the current cell
!!$           do itr=1,Ntrimax
!!$              if (itri(itr) /= 0) then
!!$                 vol_recalc(k)%cell(ix,iz)= vol_recalc(k)%cell(ix,iz) + vol_tri_eirene(itri(itr))
!!$                 Ntritot=Ntritot+1
!!$              endif
!!$           enddo
!!$           ! sum over all triangles
!!$           vol_recalc_tot = vol_recalc_tot + vol_recalc(k)%cell(ix,iz)
!!$           ! sum over all triangles that are in non-penalized cells
!!$           if (Zones(k)%masks%chi2(ix,iz) == 0) then
!!$              vol_recalc_nonpen = vol_recalc_nonpen + vol_recalc(k)%cell(ix,iz)
!!$           endif
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  close(333)
!!$
!!$  ! identify soledge2D cells that are cut (intercepted by wall), listed in soledge2D.backinterp
!!$
!!$  open(unit=333,file = 'triangles/soledge2D.triangles_mapping_edge',form = 'formatted')
!!$
!!$  read(333,*) Ncut
!!$
!!$  allocate(ked(Ncut),ixed(Ncut),ized(Ncut),itrip(Ncut,Ntrimax))
!!$  kdiff_ed=0
!!$  ixed=0
!!$  ized=0
!!$
!!$  do icell=1,Ncut
!!$     read(333,format) ked(icell),ixed(icell),ized(icell),(itrip(icell,j),j=1,Ntrimax)
!!$  enddo
!!$
!!$#if GFORTRAN==1
!!$  inquire(File='Volumes',exist=dir_e)
!!$#endif
!!$#if GFORTRAN==0
!!$  inquire(Directory='Volumes',exist=dir_e)
!!$#endif
!!$  if(.not.dir_e) then
!!$     call system("mkdir Volumes")
!!$  end if
!!$
!!$
!!$  open(unit=342,file='Volumes/check_volumes_edge',form = 'formatted')
!!$  do icell=1,Ncut
!!$     if(volumes(ked(icell))%cell(ixed(icell),ized(icell)).ne.0.D0) then
!!$        write(342,'(3i6,5(e14.7))') ked(icell),ixed(icell),ized(icell),volumes(ked(icell))%cell(ixed(icell),ized(icell)),&
!!$             vol_recalc(ked(icell))%cell(ixed(icell),ized(icell)),vol_recalc(ked(icell))%cell(ixed(icell),ized(icell))/volumes(ked(icell))%cell(ixed(icell),ized(icell))
!!$     else
!!$        write(342,'(3i6,5(e14.7))') ked(icell),ixed(icell),ized(icell),volumes(ked(icell))%cell(ixed(icell),ized(icell)),&
!!$             vol_recalc(ked(icell))%cell(ixed(icell),ized(icell)),vol_recalc(ked(icell))%cell(ixed(icell),ized(icell)),0.D0
!!$     end if
!!$  enddo
!!$  close(342)
!!$
!!$  ! now get the maximum differences (all cells)
!!$
!!$  mdiff = 0.d0
!!$  mdiff_ed=0.d0
!!$
!!$  open(unit=342,file='Volumes/volumes.log',form='formatted')
!!$
!!$  do k=1,global_parameters%N_Zones
!!$     Nx=Zones(k)%mesh%Nx
!!$     Nz=Zones(k)%mesh%Nz
!!$     md = maxval(abs(vol_recalc(k)%cell(1:Nx,1:Nz)-volumes(k)%cell(1:Nx,1:Nz))/vol_recalc(k)%cell(1:Nx,1:Nz),&
!!$          mask = vol_recalc(k)%cell(1:Nx,1:Nz) /= 0.d0 )
!!$     ml = maxloc(abs(vol_recalc(k)%cell(1:Nx,1:Nz)-volumes(k)%cell(1:Nx,1:Nz))/vol_recalc(k)%cell(1:Nx,1:Nz),&
!!$          mask = vol_recalc(k)%cell(1:Nx,1:Nz) /= 0.d0 )
!!$     is_edge=.false.
!!$
!!$     if (md > mdiff) then
!!$        ! identify wether the cell is on the boundary or not
!!$        do icell=1,Ncut
!!$           if (k == ked(Ncut) .and. ml(1) == ixed(icell) .and. & 
!!$                ml(2) == ized(icell)) then
!!$              is_edge=.true.
!!$           endif
!!$        enddo
!!$
!!$        if (is_edge) then
!!$           mdiff_ed = md
!!$           ldiff_ed = ml
!!$           kdiff_ed = k
!!$        else
!!$           if (Zones(k)%masks%chi2(ml(1),ml(2)) == 0) then
!!$              mdiff = md
!!$              ldiff = ml
!!$              kdiff = k
!!$           else
!!$              write(342,*) 'cell penalized but having a triangle in EIRENE : '
!!$              write(342,*) ' Zone = ',k,' ix = ',ml(1),' iz = ',ml(2)
!!$           endif
!!$        endif
!!$     endif
!!$  enddo
!!$
!!$  ! maximum relative difference (%)
!!$
!!$  mdiff = mdiff*100.d0
!!$  mdiff_ed = mdiff_ed*100.d0
!!$
!!$  ! now calculate the correction factor applied to sources (on solege mesh, edge cells only)
!!$
!!$  voltot_cut=0.d0
!!$
!!$  do icell=1,Ncut
!!$     k=ked(icell)
!!$     ix=ixed(icell)
!!$     iz=ized(icell)
!!$     if (volumes(k)%cell(ix,iz) == 0.d0) then
!!$        write(342,*) ' cell ix = ',ix,' iz = ',iz,' in zone k = ',k,' has zero volume'
!!$        if (Zones(k)%masks%chi2(ix,iz) == 1) then
!!$           write(342,*) ' ok,this cell is penalized...'
!!$           write(342,*) ' triangles involved = ',(itrip(icell,j),j=1,Ntrimax)
!!$
!!$        else
!!$           !  strange, stop the code
!!$           write(*,*) ' this cell is not penalized !'
!!$           call eirene_exit_own(1)
!!$        endif
!!$     endif
!!$  enddo
!!$
!!$  deallocate(ked,ixed,ized,itrip)
!!$
!!$  ! final check on total volume
!!$
!!$  voltot_check=0.d0
!!$  close(342)
!!$
!!$  open(unit=342,file='Volumes/volume.report',form='formatted')
!!$
!!$  write(342,*) '---------------------------------------------------------------------------------------------'
!!$  write(342,*) 'total number of triangles used to calculate the soledge2D cell volumes = ',Ntritot
!!$  write(342,*) 'total number of triangles in the simulation = ',NTRI_styx
!!$  write(342,*) 'total volume recalculated from EIRENE triangles (m3) = ',vol_recalc_tot
!!$  write(342,*) '(this should be equal to the volume from EIRENE)'
!!$  write(342,*) ' total volume recalculated from EIRENE triangles, excluding penalized cells (m3) = ',vol_recalc_nonpen
!!$  write(342,*) ' relative contribution of these cells (%) = ',(vol_recalc_tot-vol_recalc_nonpen)/vol_recalc_tot
!!$  write(342,*)
!!$  write(342,*) 'maximum relative difference for cell volumes appart from edge cells (%) = ', mdiff
!!$  write(342,*) ' zone = ',kdiff,' ix = ',ldiff(1),' iy = ',ldiff(2)
!!$  write(342,*) ' volume in EIRENE (cm3) = ',vol_recalc(kdiff)%cell(ldiff(1),ldiff(2))*1d6
!!$  write(342,*) ' volume in SOLEDGE (cm3) = ',volumes(kdiff)%cell(ldiff(1),ldiff(2))*1d6
!!$  write(342,*)
!!$  write(342,*) 'maximum relative difference for cell volumes, edge cells (%) = ',mdiff_ed
!!$  if(kdiff_ed /= 0) then
!!$     write(342,*) ' zone = ', kdiff_ed ,' ix = ',ldiff_ed(1),' iy = ',ldiff_ed(2)
!!$     write(342,*) ' volume in EIRENE (cm3) = ',vol_recalc(kdiff_ed)%cell(ldiff_ed(1),ldiff_ed(2))*1d6
!!$     write(342,*) ' volume in SOLEDGE (cm3) = ',volumes(kdiff_ed)%cell(ldiff_ed(1),ldiff_ed(2))*1d6
!!$  end if
!!$  write(342,*) '-----------------------------------------------------------------------------------------------'
!!$
!!$  close(342)
!!$
!!$  !write(iunout,*) ' ------------------after volume correction foer edge cells  ---------------------'
!!$
!!$  !write(iunout,*) ' Total volume in soledge after correction (m3) = ',voltot_check
!!$  !write(iunout,*) ' difference of volume introduced by edge cell correction (m3) = ',voltot_cut
!!$  !write(iunout,*) ' Relative difference with EIRENE (%) = ',abs(voltot_check-voltot_eirene)/voltot_eirene*100.d0
!!$  !write(iunout,*) ' this includes the difference in the cells at the edge/core interface, picket/fence'
!!$
!!$  ! now save the volumes (SOLEDGE and EIRENE)
!!$
!!$  do k=1,global_parameters%N_zones
!!$     filen='VS2D0'
!!$     write(filen(6:7),'(I2.2)') k
!!$     open(unit=10,file='Volumes/'//filen//'.txt',status='unknown')
!!$     do i=0,Zones(k)%mesh%Nx+1
!!$        write(10,120) (volumes(k)%cell(i,j),j=0,Zones(k)%mesh%Nz+1)
!!$     end do
!!$     close(10)
!!$     filen='VEir0'
!!$     write(filen(6:7),'(I2.2)') k
!!$     open(unit=10,file='Volumes/'//filen//'.txt',status='unknown')
!!$     do i=0,Zones(k)%mesh%Nx+1
!!$        write(10,120) (vol_recalc(k)%cell(i,j),j=0,Zones(k)%mesh%Nz+1)
!!$     end do
!!$     close(10)
!!$  enddo
!!$
!!$120 format( 512es15.7 )
!!$
!!$  do k=1,global_parameters%N_Zones
!!$     deallocate(volumes(k)%cell)
!!$     deallocate(vol_recalc(k)%cell)
!!$  enddo
!!$

end subroutine check_volume_consistency
