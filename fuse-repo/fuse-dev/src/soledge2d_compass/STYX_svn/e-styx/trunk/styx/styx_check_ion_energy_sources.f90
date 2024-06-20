subroutine styx_check_ion_energy_sources
! checks that the sum of ion energy sources over species is
! equal to the standard eirene input
! tallying made in user_dummy.f, routine uptusr 
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use styx2eirene
  implicit none
  integer :: ipls,itri
  real(dp), allocatable :: check_SEat(:)
  real(dp), allocatable :: check_SEml(:)
  real(dp), allocatable :: check_SEti(:)
  character(20) :: format
  character(3) :: nplsic

  allocate(check_SEat(1:Neir_cells))
  allocate(check_SEml(1:Neir_cells))
  allocate(check_SEti(1:Neir_cells))
  check_SEat=0._dp
  check_SEml=0._dp
  check_SEti=0._dp
  do ipls=1,nplsi
    check_SEat=check_SEat+eov(0)%eaplr(ipls,1:Neir_cells)
    check_SEml=check_SEml+eov(0)%emplr(ipls,1:Neir_cells)
    check_SEti=check_SEti+eov(0)%eiplr(ipls,1:Neir_cells)
  enddo

  ! output for sources resolved per ion

  write(nplsic,'(i3)') nplsi
  format='('//trim(adjustl(nplsic))//'(e14.7,1x))'

  open(unit=666,file='SE_at_resolved.txt',status='replace')
  do itri=1,Neir_cells
    write(666,format) (eov(0)%eaplr(ipls,itri),ipls=1,nplsi)
  enddo  
  close(666)

  open(unit=666,file='SE_ml_resolved.txt',status='replace')
  do itri=1,Neir_cells
    write(666,format) (eov(0)%emplr(ipls,itri),ipls=1,nplsi)
  enddo
  close(666)

  open(unit=666,file='SE_ti_resolved.txt',status='replace')
  do itri=1,Neir_cells
    write(666,format) (eov(0)%eiplr(ipls,itri),ipls=1,nplsi)
  enddo
  close(666)
  
  ! out of sums over species, comparison to EIRENE tallies

  
  open(unit=666,file='SE_at_sums.txt',status='replace')
  do itri=1,Neir_cells
    write(666,'(3(e14.7,1x))') eov(0)%eapl(itri),check_SEat(itri),& 
                          eov(0)%eapl(itri)-check_SEat(itri)
  enddo
  close(666)

  open(unit=666,file='SE_ml_sums.txt',status='replace')
  do itri=1,Neir_cells
    write(666,'(3(e14.7,1x))') eov(0)%empl(itri),check_SEml(itri),&
                          eov(0)%empl(itri)-check_SEml(itri)
  enddo
  close(666)

  open(unit=666,file='SE_ti_sums.txt',status='replace')
  do itri=1,Neir_cells
    write(666,'(3(e14.7,1x))') eov(0)%eipl(itri),check_SEti(itri),&
                          eov(0)%eipl(itri)-check_SEti(itri)
  enddo
  close(666)


  deallocate(check_SEat,check_SEml,check_SEti)





end subroutine styx_check_ion_energy_sources
