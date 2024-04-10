subroutine styx_fill_sheath_data
  use all_variables, only : global_parameters
  use Mphysics
  use eirmod_precision
  use eirmod_parmmod
  use styx2eirene
  use eirmod_comusr
  use eirmod_ctrig, only : PTRIX,PTRIY
  use eirmod_ccona, only : eps60
  implicit none
  integer*4 :: isurf,itri,iside,i,j
  real*8 :: Bn,BsX,BsY,BsZ,Bs,sina
  real*8 :: uperpX,uperpY,uperpZ
  real*8 :: n_X,n_Y
  character(15) :: format,format2
  character(9) :: dums
  character(3) :: Nksis,Ntaus
  integer :: iu,im,il

!!!   the following was in init_surface_data   !!!!!!!!!!!!

    do isurf=1,Nsou
      itri=recsurf(isurf)%itri
      iside=recsurf(isurf)%iside
    
      ! normal towards the plasma
      n_X=-PTRIX(iside,itri)
      n_Y=-PTRIY(iside,itri)

      Bn=BXIN(itri)*n_X+BYIN(itri)*n_Y

      sina=Bn
      sina=min(sina,1._dp)
      sina=max(sina,-1._dp)

      BsX = BXIN(itri)-Bn*n_X
      BsY = BYIN(itri)-Bn*n_Y
      BsZ = BZIN(itri)
      Bs = sqrt(BsX*BsX+BsY*BsY+BsZ*BsZ) + eps60
      
      sheath1D(isurf)%alphaB = Asin(sina)*180._dp/pi
      sheath1D(isurf)%uparX   = BsX/Bs
      sheath1D(isurf)%uparY   = BsY/Bs
      sheath1D(isurf)%uparZ   = BsZ/Bs

      !write(*,*) ' upar**2 = ',sheath1D(isurf)%uparX**2+sheath1D(isurf)%uparY**2+sheath1D(isurf)%uparZ**2
      !write(*,*) ' upar.n = ',sheath1D(isurf)%uparX*PTRIX(iside,itri)+sheath1D(isurf)%uparY*PTRIY(iside,itri)

      ! calculation of ksi0 (SI units*sqrt(1e6) -> density in cm-3 to calculate ksi
      if (BFIN(itri) /= 0) then
        !##NS
        sheath1D(isurf)%ksi0 = sqrt(global_parameters%element_list(1)%mass*m_u&
             /(epsilon_0*BFIN(itri)*BFIN(itri)))*1e3_dp
      else
        write(*,*) ' Warning from init surface, Zero B field strength in triangle # itri = ',itri
      endif

      ! direction orthogonal to normal and upar (PTRIZ = 0 !)
 
      uperpX = -sheath1D(isurf)%uparZ*n_Y
      uperpY = sheath1D(isurf)%uparZ*n_X
      uperpZ = sheath1D(isurf)%uparX*n_Y-sheath1D(isurf)%uparY*n_X

      !write(*,*) ' uper**2 = ',uperpX**2+uperpY**2+uperpZ**2
      !write(*,*) ' uper.n = ',uperpX*PTRIX(iside,itri)+uperpY*PTRIY(iside,itri)
      !write(*,*) ' upar.uper = ',uperpX*sheath1D(isurf)%uparX+uperpY*sheath1D(isurf)%uparY+uperpZ*sheath1D(isurf)%uparZ


      ! construct change of coordinate frame matrix (Bnu2xyz)ij = (ui.uj')

      Bnu2xyz(1,1,isurf)=sheath1D(isurf)%uparX
      Bnu2xyz(1,2,isurf)=n_X
      Bnu2xyz(1,3,isurf)=uperpX

      Bnu2xyz(2,1,isurf)=sheath1D(isurf)%uparY
      Bnu2xyz(2,2,isurf)=n_Y
      Bnu2xyz(2,3,isurf)=uperpY

      Bnu2xyz(3,1,isurf)=sheath1D(isurf)%uparZ
      Bnu2xyz(3,2,isurf)=0._dp
      Bnu2xyz(3,3,isurf)=uperpZ

      xyz2Bnu(:,:,isurf)=transpose(Bnu2xyz(:,:,isurf))    

    enddo


    if (sheath_model == 1) then

    !! read the sheath database

      open(unit=555,file='Sheathdata/sheath1D_database',status='old')
      
      read(555,*)
      read(555,*)
      read(555,'(i3,1x,i3,1x,i3)') sheath1D_av%Nalphab,sheath1D_av%Ntau,sheath1D_av%Nksi
      
      ! grids
      allocate(sheath1D_av%alphaB(sheath1D_av%Nalphab))
      allocate(sheath1D_av%tau(sheath1D_av%Ntau))
      allocate(sheath1D_av%ksi(sheath1D_av%Nksi))
      ! data
      allocate(sheath1D_av%E(sheath1D_av%Ntau,sheath1D_av%NalphaB,sheath1D_av%Nksi))
      allocate(sheath1D_av%alpha(sheath1D_av%Ntau,sheath1D_av%NalphaB,sheath1D_av%Nksi))
      allocate(sheath1D_av%beta(sheath1D_av%Ntau,sheath1D_av%NalphaB,sheath1D_av%Nksi))

      write(Nksis,'(i3)') sheath1D_av%Nksi
      format='('//adjustl(trim(Nksis))//'(f10.2,1x))'
      format=trim(format)
 
      write(Ntaus,'(i3)') sheath1D_av%Ntau
      format2='('//adjustl(trim(Ntaus))//'(f10.2,1x))'
      format2=trim(format2)
  
      read(555,*)
      read(555,format2) sheath1D_av%tau(1:sheath1D_av%Ntau)

      read(555,*)
      read(555,format)  sheath1D_av%ksi(1:sheath1D_av%Nksi)
    
      do j=1,sheath1D_av%NalphaB
        read(555,'(A9,f10.2)') dums,sheath1D_av%alphaB(j)
        do i=1,sheath1D_av%Ntau
          read(555,format) sheath1D_av%E(i,j,1:sheath1D_av%Nksi)
        enddo
        do i=1,sheath1D_av%Ntau
          read(555,format) sheath1D_av%alpha(i,j,1:sheath1D_av%Nksi)
        enddo
        do i=1,sheath1D_av%Ntau
          read(555,format) sheath1D_av%beta(i,j,1:sheath1D_av%Nksi)
        enddo
      enddo

      close(555)

    
      if (any(sheath1D_av%ksi == 0._dp)) then
        write(*,*) ' ksi = 0 in sheath database, would crash when denormalizing energy. Exit ...'
        call eirene_exit_own(1)
      endif


    ! now determine alphaB for each wall elements, and store the corresponding data 
    ! (E,alpha,beta)=f(tau,ksi|alphaB(isurf))

      do i=1,Nsou
      ! look for sheath1D(i)%alphaB in alphaB grid, binary search
     
        il=0
        iu=sheath1D_av%NalphaB
     
        do while (iu-il > 1)
          im=0.5*(iu+il)
          if (abs(sheath1D(i)%alphaB) >= sheath1D_av%alphaB(im)) then
            il=im
          else
            iu=im
          endif
        enddo

        if (iu < 1) iu=1
        if (iu > sheath1D_av%NalphaB) iu=sheath1D_av%NalphaB

        allocate(sheath1D(i)%E(sheath1D_av%Ntau,sheath1D_av%Nksi))
        allocate(sheath1D(i)%alpha(sheath1D_av%Ntau,sheath1D_av%Nksi))
        allocate(sheath1D(i)%beta(sheath1D_av%Ntau,sheath1D_av%Nksi)) 

        ! for the moment takes the closest value (from above)
        sheath1D(i)%E(:,:)      =  sheath1D_av%E(:,iu,:)
        sheath1D(i)%alpha(:,:)  =  sheath1D_av%alpha(:,iu,:)
        ! determine sign of beta from sign of alphaB (drift change direction with sign of B)
        sheath1D(i)%beta(:,:)   =  sheath1D_av%beta(:,iu,:)*sign(1._dp,sheath1D(i)%alphaB)

      enddo

      ! put angles alpha and beta in radians

      do i=1,Nsou
        sheath1D(i)%alpha = sheath1D(i)%alpha*pi/180._dp
        sheath1D(i)%beta  = sheath1D(i)%beta*pi/180._dp
      enddo


     ! deallocate data not needed an_Ymore

      deallocate(sheath1D_av%E,sheath1D_av%alpha,sheath1D_av%beta)

  endif

end subroutine styx_fill_sheath_data

