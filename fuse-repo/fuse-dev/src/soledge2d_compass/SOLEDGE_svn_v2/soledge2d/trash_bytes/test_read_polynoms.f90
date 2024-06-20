
program test_read_polynoms
  use data_structures
  implicit none
  
  character(2) :: element_name
  character(30) :: file_recombination,file_ionization
  character(30) :: file_plt,file_prb
  integer Z 

  element_name='O '
  Z=8

  allocate(element_%recombination(Z))
  allocate(element_%ionization(Z))
  allocate(element_%plt(Z))
  allocate(element_%prb(Z))

  file_recombination = 'acd96_o_fits_dpoly8.dat'
  file_ionization    = 'scd96_o_fits_dpoly8.dat'
  file_plt           = 'plt96_o_fits_dpoly8.dat'
  file_prb           = 'prb96_o_fits_dpoly8.dat' 

  call read_amdata(element_name,Z,'recombination',file_recombination,element_%recombination)
  call read_amdata(element_name,Z,'ionization   ',file_ionization,element_%ionization)
  call read_amdata(element_name,Z,'plt          ',file_plt,element_%plt)
  call read_amdata(element_name,Z,'prb          ',file_prb,element_%prb)

! now test of data

  call write_rate_coeffs(element_%recombination,'recombination',Z)
  call write_rate_coeffs(element_%ionization,'ionization',Z)
  call write_rate_coeffs(element_%plt,'plt',Z)
  call write_rate_coeffs(element_%prb,'prb',Z)

end program

subroutine read_amdata(element_name,Z,data_type,file_name,fit_data_)
  use data_structures
  implicit none
  character(2), intent(in) :: element_name
  integer, intent(in) :: Z
  character(30), intent(in) :: file_name
  character(13), intent(in) :: data_type
  type(fit_data) :: fit_data_(Z)  !intent(out)
  integer :: i,j,Zr,iz,izmin,izmax,nl,nr,nf
  real*8 :: r
  character(2) :: ncolc,nrc
  character(11) :: format,format2

  open(unit=666,file=file_name,status='old')

  do i=1,8
    read(666,*)
  enddo

  read(666,'(4x,i2)') Zr

  if (Zr /= Z) then
    write(*,*) ' inconsistency between Z values !'
    stop 
  endif

  read(666,'(14x,i2,8x,i2)') izmin,izmax
  read(666,*)
  read(666,'(20x,i2)') fit_data_(1)%degree

 
  read(666,*)
  read(666,'(10x,i3)') fit_data_(1)%nne
  read(666,'(10x,e19.12)') fit_data_(1)%Nemin
  read(666,'(10x,e19.12)') fit_data_(1)%Nemax

  do i=1,Zr
    allocate(fit_data_(i)%Ne(fit_data_(1)%nne))
  enddo
  
  read(666,*)
 
  r=fit_data_(1)%nne/ncol  
  nl=floor(r)

  write(ncolc,'(i2)') ncol
  format='('//trim(adjustl(ncolc))//'(e19.12))'

  do i=1,nl
    read(666,format) (fit_data_(1)%Ne(j+(i-1)*ncol),j=1,ncol)
  enddo

  nr=fit_data_(1)%nne-nl*ncol

  if (nr>0) then
    write(nrc,'(i1)') nr
    format2='('//trim(adjustl(nrc))//'(e19.12))'
    read(666,format2) (fit_data_(1)%Ne(j+nl*ncol),j=1,nr)
  endif

  do i=2,Zr
    fit_data_(i)%nne    = fit_data_(1)%nne
    fit_data_(i)%Nemin  = fit_data_(1)%Nemin
    fit_data_(i)%Nemax  = fit_data_(i)%Nemax
    fit_data_(i)%Ne     = fit_data_(1)%Ne
    fit_data_(i)%degree = fit_data_(1)%degree
  enddo

  ! now read temperature data

  do iz=izmin,izmax

    read(666,*)
    read(666,*)
    read(666,*)
    
    read(666,'(10x,i3)') fit_data_(iz)%nTe
    read(666,'(10x,e19.12)') fit_data_(iz)%Temin
    read(666,'(10x,e19.12)') fit_data_(iz)%Temax

    read(666,*)

    allocate(fit_data_(iz)%Te(fit_data_(iz)%nTe))

    r=fit_data_(iz)%nTe/ncol  
    nl=floor(r)

    do i=1,nl
      read(666,format) (fit_data_(iz)%Te(j+(i-1)*ncol),j=1,ncol)
    enddo

    nr=fit_data_(iz)%nTe-nl*ncol

    if (nr>0) then
      write(nrc,'(i1)') nr
      format2='('//trim(adjustl(nrc))//'(e19.12))'
  
      read(666,format2) (fit_data_(iz)%Te(j+nl*ncol),j=1,nr)
    endif

    read(666,*)

    ! number of fitting coefficients

    nf=(fit_data_(iz)%degree+1)*fit_data_(iz)%degree/2
   
    allocate(fit_data_(iz)%coeffs(nf))

    r=nf/ncol  
    nl=floor(r)

    do i=1,nl
      read(666,format) (fit_data_(iz)%coeffs(j+(i-1)*ncol),j=1,ncol)
    enddo

    nr=nf-nl*ncol

    if (nr>0) then
      write(nrc,'(i1)') nr
      format2='('//trim(adjustl(nrc))//'(e19.12))'  
      read(666,format2) (fit_data_(iz)%coeffs(j+nl*ncol),j=1,nr)
    endif

    ! reading done, now goodness of fit

    read(666,*)
    read(666,'(33x,e8.2)') fit_data_(iz)%goodness
    read(666,*)
    read(666,*)


  enddo
  

  close(666)  


end subroutine read_amdata

subroutine write_rate_coeffs(fit_data_,data_type,Z)
  use data_structures
  implicit none
  integer,intent(in) :: Z
  type(fit_data), intent(in) :: fit_data_(Z)
  character(*), intent(in) :: data_type
  real*8, allocatable :: Q(:,:),pow_logNe(:,:),pow_logTe(:,:),logNe(:),logTe(:),Te_a(:),rate_a(:)
  integer :: Nx,Ny,ix,iy,i,j,iz,n,degree1,nTemax
  type(rate), allocatable :: rates(:)
  character(2) :: nTec,nTemaxc
  character(12) :: format
  character(21) :: file_name
  character(6) :: file_name2

  ! recalculate the rates on the ADAS (ne,Te) grid for checking purposes
  ! see this like a loop on Nx,Nz in a zone
  
 
  allocate(rates(Z))

  do iz=1,Z
    Nx=fit_data_(iz)%nne
    Ny=fit_data_(iz)%nTe

    allocate(logNe(Nx),logTe(Ny))
    allocate(Q(Nx,Ny),rates(iz)%values(Nx,Ny))
    
    degree1=fit_data_(iz)%degree-1

    allocate(pow_logNe(0:degree1,Nx))
    allocate(pow_logTe(0:degree1,Ny))
 
    ! precalculate the powers

    logNe=log10(fit_data_(iz)%Ne)
    logTe=log10(fit_data_(iz)%Te)

    pow_logne(0,:)=1.d0
    pow_logTe(0,:)=1.d0
    do i=1,degree1
      pow_logne(i,:)=pow_logne(i-1,:)*logNe
      pow_logTe(i,:)=pow_logTe(i-1,:)*logTe  
    enddo

    deallocate(logNe,logTe)

    Q=0.d0

    do ix=1,Nx
      do iy=1,Ny
        n=0
        do i=0,degree1
          do j=0,degree1-i
            n=n+1
            Q(ix,iy)=Q(ix,iy)+fit_data_(iz)%coeffs(n)*pow_logne(i,ix)*pow_logTe(j,iy)
          enddo
        enddo                
      enddo
    enddo
  
    rates(iz)%values(:,:)=exp(log(10.d0)*Q(:,:))

    deallocate(Q,pow_logNe,pow_logTe)

  enddo


  select case(trim(data_type))

    case('recombination')
      file_name = 'rate_coefficients_acd'
      file_name2 = 'Te_acd'
    case('ionization')
      file_name = 'rate_coefficients_scd'
      file_name2 = 'Te_scd'
    case('plt')
      file_name = 'rate_coefficients_plt'
      file_name2 = 'Te_plt'
    case('prb')
      file_name = 'rate_coefficients_prb'
      file_name2 = 'Te_prb'

  end select
  

! get the temperature grid out for each iz (in a format easy to read in matlab)

  nTemax=0.d0
  do iz=1,Z
    if (fit_data_(iz)%nTe> nTemax) nTemax=fit_data_(iz)%nTe
  enddo

  allocate(Te_a(nTemax))
 
  write(nTemaxc,'(i2)') nTemax
  format='('//nTemaxc//'(e14.7))'

  allocate(rate_a(nTemax))

  open(unit=666,file=file_name,status='replace')

  do iz=1,Z   
    Nx=fit_data_(iz)%nne
    Ny=fit_data_(iz)%nTe
    do ix=1,Nx
      rate_a=0.d0
      rate_a(1:Ny)=rates(iz)%values(ix,:)
      write(666,format) (rate_a(iy),iy=1,NTemax)
    enddo
  enddo

  close(666)


  open(unit=666,file=file_name2,status='replace')

  do iz=1,Z
    Te_a=0.d0
    Te_a(1:fit_data_(iz)%nTe)=fit_data_(iz)%Te
    write(666,format) (Te_a(j),j=1,NTemax)
  enddo

  close(666)

  deallocate(Te_a,rate_a)



end subroutine write_rate_coeffs
