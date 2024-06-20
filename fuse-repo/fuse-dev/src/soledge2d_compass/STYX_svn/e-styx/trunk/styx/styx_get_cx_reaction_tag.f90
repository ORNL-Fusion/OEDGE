subroutine styx_get_cx_reaction_data(reac_tag,cx_database1,cx_database2,fit_coeffs)
  use styx2eirene, only : direct_coupling
  use Mreaders
  implicit none
  character(8), intent(in) :: reac_tag
  character(29), intent(out) :: cx_database1,cx_database2
  character(8) :: tag
  character(7) :: buff
  integer :: icount,ireac,nreac,unit_in
  real*8, dimension(9) :: fit_coeffs

  unit_in=667

  open(unit=unit_in,file='cx_database',status='old')
 
  call skip_line(unit_in,1)
  read(unit_in,'(A7,i3)') buff,nreac
  call skip_line(unit_in,2)  
  
  icount=0
  do ireac=1,nreac
    read(unit_in,'(A8)') tag
    if (trim(adjustl(reac_tag)) == trim(adjustl(tag))) then
      if (direct_coupling) then
        call skip_line(unit_in,2)
        read(unit_in,'(A29)') cx_database1
        read(unit_in,'(A29)') cx_database2
        if (ireac<nreac) call skip_line(unit_in,7)
        fit_coeffs=0.d0
      else
        call skip_line(unit_in,2)
        read(unit_in,'(A29)') cx_database1
        read(unit_in,'(A29)') cx_database2
        call skip_line(unit_in,2)
        read(unit_in,'(5(1x,e20.12))') fit_coeffs(1:5)  
        read(unit_in,'(4(1x,e20.12))') fit_coeffs(6:9)
        if (ireac<nreac) call skip_line(unit_in,3)
      endif
      icount=icount+1    
    else
      if (ireac<nreac) call skip_line(unit_in,11)      
    endif
  enddo

  if (icount==0) then
    write(*,*) 'cx reaction not found in cx_database for reac_tag = ',reac_tag
    stop
  endif

  if (icount>1) then
    write(*,*) 'Several hits in cx_database for reac_tag = ',reac_tag
  endif


  close(unit_in)

end subroutine
