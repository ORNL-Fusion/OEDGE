module interpolation_soleir

  use domain_types
  use equimag
  implicit none

  type :: InterpData1
     integer*4,allocatable :: Corners(:,:,:) ! table that is associated to a zone
                                             ! size: Nx * Nz * 4
                                             ! for each soledge2d grid point (cell), contains 
                                             ! eirene grid knot number for each soledge2d cell corner
  end type InterpData1

  type :: InterpData2
     integer*4 nknots
     integer*4 ntriangles
     integer*4,allocatable :: knotsNeighNumb(:) ! for each knot tells how many neighboring soledge2d grid points
     integer*4,allocatable :: Knots(:,:,:)      ! for each knot saves neighboring soledge2d grid point coord (i,j,Nzone)
     real*8,allocatable :: Knots_R(:) ! Radial coordinate of knots  
     real*8,allocatable :: Knots_Z(:) ! Vertical coordinate of knots
     real*8,allocatable :: Knots_N(:) ! Density on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_G(:) ! Gamma on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_Te(:) ! Te on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_Ti(:) ! Ti on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_Br(:) ! Br on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_Bz(:) ! Bz on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_Bphi(:) ! Bphi on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_B(:) ! B on knots (Soledge-->Eirene)
     real*8,allocatable :: Knots_Sn(:) ! Density source on knots (Eirene-->Soledge)
     real*8,allocatable :: Knots_SG(:) ! Momentum source on knots (Eirene-->Soledge)
     real*8,allocatable :: Knots_SUe(:) ! Electrons internal energy source on knots (Eirene-->Soledge)
     real*8,allocatable :: Knots_SUi(:) ! Ions internal energy source on knots (Eirene-->Soledge)
     integer*4,allocatable:: zone_info(:,:) ! soledge2d cell that contains triangle
     integer*4,allocatable:: flux_info(:,:) ! for each triangle side, tells which soledge2d cell side it is:
                            ! 1:N   2:S   3:E   4:W   0:None
  end type InterpData2

contains

  ! ######################################################################################### !
  ! ###  Subroutine that initiate the tables used to interpolate from soledge2d grid      ### !
  ! ###  to eirene triangles grid                                                         ### !
  ! ######################################################################################### !
  subroutine init_interp(Zones,NZones,InterpData1s,InterpData2_)
    implicit none
    integer*4 NZones
    Type(Zone) :: Zones(NZones)
    Type(InterpData1) :: InterpData1s(NZones)
    Type(Interpdata2) :: InterpData2_
    integer*4 V2i(2),V4i(4),V8i(8)
    integer*4 V2j(2),V4j(4),V8j(8)
    integer*4 V2k(2),V4k(4),V8k(8)
    real*8 V3(3)
    integer*4 V7(7)
    integer*4 Nx,Nz
    integer*4 k,j,i,nbr
    character(8) file
    character(12) file2
    integer*4 nknots,ntri
    ! Loading data concerning soledge to eirene interpolation    
    open(unit=20,file='triangles/vert2knots.txt',status='unknown')
    open(unit=21,file='triangles/soledge2D.npco_char',status='unknown')
    open(unit=22,file='triangles/soledge2D.zones',status='unknown')
    read(20,*) InterpData2_%nknots
    nknots=InterpData2_%nknots
    read(22,*) InterpData2_%ntriangles
    ntri=InterpData2_%ntriangles
    allocate(InterpData2_%knotsNeighNumb(1:nknots))
    allocate(InterpData2_%knots(1:nknots,1:8,1:3))
    allocate(InterpData2_%knots_R(1:nknots))
    allocate(InterpData2_%knots_Z(1:nknots))
    allocate(InterpData2_%knots_N(1:nknots))
    allocate(InterpData2_%knots_G(1:nknots))
    allocate(InterpData2_%knots_Te(1:nknots))
    allocate(InterpData2_%knots_Ti(1:nknots))
    allocate(InterpData2_%knots_Br(1:nknots))
    allocate(InterpData2_%knots_Bz(1:nknots))
    allocate(InterpData2_%knots_Bphi(1:nknots))
    allocate(InterpData2_%knots_B(1:nknots))
    allocate(InterpData2_%zone_info(1:ntri,1:3))
    allocate(InterpData2_%flux_info(1:ntri,1:3))
    do k=1,Nknots
       read(20,*) nbr
       read(20,*) InterpData2_%knotsNeighNumb(k)
       if(InterpData2_%knotsNeighNumb(k).eq.1) then
          read(20,*) InterpData2_%knots(k,1,3)
          read(20,*) InterpData2_%knots(k,1,1)
          read(20,*) InterpData2_%knots(k,j,2)
       else
          if(InterpData2_%knotsNeighNumb(k).eq.2) then
             read(20,*) (V2k(i),i=1,2)
             read(20,*) (V2i(i),i=1,2)
             read(20,*) (V2j(i),i=1,2)
             do j=1,2
                InterpData2_%knots(k,j,1)=V2i(j)
                InterpData2_%knots(k,j,2)=V2j(j)
                InterpData2_%knots(k,j,3)=V2k(j)
             end do
          else
             if(InterpData2_%knotsNeighNumb(k).eq.4) then
                read(20,*) (V4k(i),i=1,4)
                read(20,*) (V4i(i),i=1,4)
                read(20,*) (V4j(i),i=1,4)
                do j=1,4
                   InterpData2_%knots(k,j,1)=V4i(j)
                   InterpData2_%knots(k,j,2)=V4j(j)
                   InterpData2_%knots(k,j,3)=V4k(j)
                end do
             else
                read(20,*) (V8k(i),i=1,8)
                read(20,*) (V8i(i),i=1,8)
                read(20,*) (V8j(i),i=1,8)
                do j=1,8
                   InterpData2_%knots(k,j,1)=V8i(j)
                   InterpData2_%knots(k,j,2)=V8j(j)
                   InterpData2_%knots(k,j,3)=V8k(j)
                end do
             end if
          end if
       end if
       read(21,*) (V3(i),i=1,3)
       InterpData2_%knots_R(k)=V3(2)
       InterpData2_%knots_Z(k)=V3(3)
    end do
    do k=1,ntri
       read(22,*) (V7(i),i=1,7)
       InterpData2_%zone_info(k,1)=V7(2)
       InterpData2_%zone_info(k,2)=V7(3)
       InterpData2_%zone_info(k,3)=V7(4)
       InterpData2_%flux_info(k,1)=V7(5)
       InterpData2_%flux_info(k,2)=V7(6)
       InterpData2_%flux_info(k,3)=V7(7)
    end do
    close(20)
    close(21)
    close(22)
    ! Loading data concerning eirene to soledge interpolation    
    do k=1,NZones
       Nx=Zones(k)%Nx
       Nz=Zones(k)%Nz
       allocate(InterpData1s(k)%Corners(1:Nx,1:Nz,1:4))
       ! loading corners data
       file='triA_000'
       write(file(6:8),'(i3.3)') k
       open(unit=10,file='triangles/'//file//'.txt',status='unknown')
       file='triB_000'
       write(file(6:8),'(i3.3)') k
       open(unit=11,file='triangles/'//file//'.txt',status='unknown')
       file='triC_000'
       write(file(6:8),'(i3.3)') k
       open(unit=12,file='triangles/'//file//'.txt',status='unknown')
       file='triD_000'
       write(file(6:8),'(i3.3)') k
       open(unit=13,file='triangles/'//file//'.txt',status='unknown')
       do i=1,Nx
          read(10,*) (InterpData1s(k)%Corners(i,j,1),j=1,Nz)
          read(11,*) (InterpData1s(k)%Corners(i,j,2),j=1,Nz)
          read(12,*) (InterpData1s(k)%Corners(i,j,3),j=1,Nz)
          read(13,*) (InterpData1s(k)%Corners(i,j,4),j=1,Nz)
       end do
       close(10)
       close(11)
       close(12)
       close(13)
    end do
  end subroutine init_interp


  ! ######################################################################################### !
  ! ###  Subroutine that interpolates density, Gamma, Te and Ti fields solved on          ### !
  ! ###  soledge2d grid to the corners of the triangles constituing eirene grid.          ### !
  ! ###  This routine also gives the value of fluxes on triangles segments where BC       ### !
  ! ###  apply.                                                                           ### !
  ! ######################################################################################### !
  subroutine sol2eir(InterpData2_,NZones,Zones,MagFields)
    implicit none
    integer*4 NZones
    Type(Zone):: Zones(NZones)
    Type(InterpData2) :: InterpData2_
    Type(MagField) :: MagFields(Nzones)
    integer*4 k,n
    integer*4 nneigh
    real*8 d
    real*8 sumN,sumG,sumTe,sumTi,sumd
    real*8 sumBr,sumBz,sumBphi,sumB
    integer*4 ii,ij,ik
    do k=1,InterpData2_%Nknots
       nneigh=InterpData2_%knotsNeighNumb(k)
       sumN=0.d0
       sumG=0.d0
       sumTe=0.d0
       sumTi=0.d0
       sumd=0.d0
       sumBr=0.d0
       sumBz=0.d0
       sumBphi=0.d0
       sumB=0.d0
       do n=1,nneigh !##CAREFUL## the case nneigh=1 is quite raw
          ii=InterpData2_%knots(k,n,1)
          ij=InterpData2_%knots(k,n,2)
          ik=InterpData2_%knots(k,n,3)
          d=sqrt((InterpData2_%knots_r(k)-MagFields(ik)%R(ii,ij))**(2.d0)+&
               (InterpData2_%knots_z(k)-MagFields(ik)%z(ii,ij))**(2.d0))
          sumN=sumN+1.d0/d*Zones(ik)%plasma(2)%density(ii,ij)
          sumG=sumG+1.d0/d*Zones(ik)%plasma(2)%Gamma(ii,ij)
          sumTe=sumTe+1.d0/d*Zones(ik)%plasma(2)%Te(ii,ij)
          sumTi=sumTi+1.d0/d*Zones(ik)%plasma(2)%Ti(ii,ij)
          !could be done just once
          sumBr=sumBr+1.d0/d*MagFields(ik)%Br(ii,ij)
          sumBz=sumBz+1.d0/d*MagFields(ik)%Bz(ii,ij)
          sumBphi=sumBphi+1.d0/d*MagFields(ik)%Bphi(ii,ij)
          sumB=sumB+1.d0/d*MagFields(ik)%B(ii,ij)
          sumd=sumd+1.d0/d
       end do
       InterpData2_%knots_n(k)=sumN/sumd
       InterpData2_%knots_G(k)=sumG/sumd
       InterpData2_%knots_Te(k)=sumTe/sumd
       InterpData2_%knots_Ti(k)=sumTi/sumd
       InterpData2_%knots_Br(k)=sumBr/sumd
       InterpData2_%knots_Bz(k)=sumBz/sumd
       InterpData2_%knots_Bphi(k)=sumBphi/sumd
       InterpData2_%knots_B(k)=sumB/sumd
    end do
  end subroutine sol2eir


  ! ######################################################################################### !
  ! ###  Subroutine that interpolates the results obtained on eirene triangles corners    ### !
  ! ###  to soledge2d grid                                                                ### !
  ! ######################################################################################### !
  subroutine eir2sol(Zones,NZones,InterpData1s,MagFields,InterpData2_)
    implicit none
    integer*4 NZones
    Type(Zone) :: Zones(Nzones)
    Type(InterpData1) :: InterpData1s(NZones)
    Type(InterpData2) :: InterpData2_
    Type(MagField) :: MagFields(NZones)
    integer*4 i,j,k,n,Nx,Nz
    integer*4 nknot
    real*8 d
    real*8 sumd,sumSn,sumSG,sumSUe,sumSUi
    do k=1,NZones
       Nx=Zones(k)%Nx
       Nz=Zones(k)%Nz
       do i=1,Nx
          do j=1,Nz
             sumd=0.d0
             sumSn=0.d0
             sumSG=0.d0
             sumSUe=0.d0
             sumSUi=0.d0
             do n=1,4
                nknot=InterpData1s(k)%Corners(i,j,n)
                d=sqrt((MagFields(k)%R(i,j)-InterpData2_%knots_R(nknot))**2.d0+&
                     (MagFields(k)%Z(i,j)-InterpData2_%knots_z(nknot))**2.d0)
                sumSn=sumSn+1.d0/d*InterpData2_%knots_Sn(nknot)
                sumSG=sumSG+1.d0/d*InterpData2_%knots_SG(nknot)
                sumSUe=sumSn+1.d0/d*InterpData2_%knots_SUe(nknot)
                sumSUi=sumSn+1.d0/d*InterpData2_%knots_SUi(nknot)
                sumd=sumd+1.d0/d
             end do
             Zones(k)%plasma(2)%Sn_n(i,j)=sumSn/sumd
             Zones(k)%plasma(2)%Sn_G(i,j)=sumSG/sumd
             Zones(k)%plasma(2)%Sn_Ue(i,j)=sumSUe/sumd
             Zones(k)%plasma(2)%Sn_Ui(i,j)=sumSUi/sumd
          end do
       end do
    end do
  end subroutine eir2sol

  subroutine write_res_knots(InterpData2_)
    implicit none
    Type(Interpdata2) :: InterpData2_
    integer*4 V(4)
    integer*4 i,j
    open(unit=20,file='triangles/soledge2D.elemente',status='unknown')
    open(unit=21,file='Results/density.txt',status='unknown')
    open(unit=24,file='Results/Gamma.txt',status='unknown')
    open(unit=25,file='Results/Te.txt',status='unknown')
    open(unit=26,file='Results/Ti.txt',status='unknown')
    open(unit=22,file='Results/R.txt',status='unknown')
    open(unit=23,file='Results/Z.txt',status='unknown')
    do i=1,InterpData2_%ntriangles
       read(20,*) (V(j),j=1,4)
       write(21,100) InterpData2_%knots_n(V(2)), InterpData2_%knots_n(V(3)),&
            InterpData2_%knots_n(V(4))
       write(22,100) InterpData2_%knots_R(V(2)), InterpData2_%knots_R(V(3)),&
            InterpData2_%knots_R(V(4))
       write(23,100) InterpData2_%knots_Z(V(2)), InterpData2_%knots_Z(V(3)),&
            InterpData2_%knots_Z(V(4))
       write(24,100) InterpData2_%knots_G(V(2)), InterpData2_%knots_G(V(3)),&
            InterpData2_%knots_G(V(4))
       write(25,100) InterpData2_%knots_Te(V(2)), InterpData2_%knots_Te(V(3)),&
            InterpData2_%knots_Te(V(4))
       write(26,100) InterpData2_%knots_Ti(V(2)), InterpData2_%knots_Ti(V(3)),&
            InterpData2_%knots_Ti(V(4))
    end do
100 format(3es15.7)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
    close(26)
  end subroutine write_res_knots

end module interpolation_soleir
