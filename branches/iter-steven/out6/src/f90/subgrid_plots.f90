module subgrid_plots

  use error_handling
  use subgrid_options
  use subgrid



contains


  subroutine load_subgrid_array(subgrid_data,iselect,istate,itype,ylab,blab,ref,nizs,ierr)
    use hc_get
    implicit none

    integer :: iselect,istate,itype,nizs,ierr
    real :: subgrid_data(sg_rdim,sg_zdim)
    character*(*) :: ylab,blab,ref
    include 'adas_data_spec'

    !
    ! This routine loads the subgrid data selected as well as loading the 
    ! appropriate labels. To load the labels it uses the same routines
    ! as the load_divdata_array routine
    !
    ! 
    ! Supported ISELECT values for subgrids:
    !
    !              32 = Subgrid impurity density - STATE = IZ
    !              33 = Subgrid HC density - STATE = HC STATE INDEX
    !              34 = Subgrid impurity ADAS based emissions - additional data read 
    !              35 = Subgrid HC Emission

    real :: mfact
    integer :: cion
    real :: absfac
    !
    !
    !     ADAS variables
    !
    CHARACTER ADASID*80,graph3*80
    CHARACTER XFESYM*2
    character adasex*3
    integer   adasyr
    INTEGER ISELE,ISELR,ISELX,iseld,ircode
    integer line 
    REAL WLNGTH
    integer :: len,lenstr
    external lenstr

    !
    ! Create arrays for local storage
    !
    real,allocatable :: exc_state_density(:,:)
    real,allocatable :: rec_state_density(:,:)
    real,allocatable :: raxis(:),zaxis(:)
    integer :: iflag

    !
    ! Initialize 
    !
    ierr=0
    mfact = 1.0
    cion = gcion()
    absfac = gabsfac()

    !
    !     Set the YLAB value
    !
    call set_ylab(iselect,istate,itype,nizs,ylab) 
    !
    !     Set the BLAB value
    !
    call set_blab(iselect,istate,itype,nizs,blab) 


    if (iselect.eq.32.or.iselect.eq.33) then 

       !
       ! All options require loading the subgrid data
       !

       call load_subgrid_xsection(subgrid_data,iselect,istate) 

       !
       !        Scaling factor 
       !
       !WRITE (0,*) 'SCALE:',ABSFAC,MFACT

       IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC


       subgrid_data = subgrid_data * mfact

    elseif (iselect.eq.34) then 



       !
       !----------------------------------------------------------
       !
       !     ADAS based - Impurity spectral line
       !
       !----------------------------------------------------------
       !
       !
       !
       !        Need to read in ADAS data spec to calculate radiation
       !
       if (cadas_switch.eq.0) then  

          CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,ISELD,IERR)
          !
          !           Save the ADAS data read in into the common block for 
          !           possible re-use. Do not set the cadas_switch.
          !
          cadasid = adasid
          cadasyr = adasyr
          cadasex = adasex
          cisele  = isele
          ciselr  = iselr
          ciselx  = iselx
          ciseld  = iseld
          !
          !        Use ADAS data in common instead of reading from input 
          !
       elseif (cadas_switch.eq.1) then 
          !            
          adasid = cadasid
          adasyr = cadasyr
          adasex = cadasex
          isele  = cisele
          iselr  = ciselr
          iselx  = ciselx
          iseld  = ciseld

       endif

       if (ierr.ne.0) return 

       IF (ISTATE.GE.0.AND.ISTATE.LE.NIZS.AND.ISTATE.LT.CION)THEN

          !
          ! Allocate storage for emission data arrays and axes
          !
          allocate(exc_state_density(sg_rdim,sg_zdim),stat=iflag)
          allocate(rec_state_density(sg_rdim,sg_zdim),stat=iflag)
          allocate(raxis(sg_rdim),stat=iflag)
          allocate(zaxis(sg_zdim),stat=iflag)

          call calc_subgrid_axis(raxis,sg_rmin,sg_rmax,sg_rdim)
          call calc_subgrid_axis(zaxis,sg_zmin,sg_zmax,sg_zdim)

          call load_subgrid_xsection(exc_state_density,iselect,istate)

          if (istate.lt.nizs) then 
             call load_subgrid_xsection(rec_state_density,iselect,istate+1)
          else
             rec_state_density = 0.0
          endif

          call LDADAS_RZ(CION,ISTATE,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,Wlngth,IRCODE,&
               &subgrid_data,exc_state_density,rec_state_density,raxis,zaxis,sg_rdim,sg_zdim)

          IF (IRCODE.NE.0) THEN
             call errmsg('LOAD_SUBGRID_ARRAY:LDADAS ERROR',ircode)
             ierr = 1
             return   
          ENDIF

          REF = 'ADAS PLRP XX XXXXX ('

          WRITE(REF(11:12),'(I2)') ISTATE
          WRITE(REF(14:18),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(REF)
          IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
          LEN = LENSTR(REF)
          IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
          LEN = LENSTR(REF)
          IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // ') '
          LEN = LENSTR(REF)

       endif
       !
       !        Scale by MFACT
       !
       IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC

       !
       ! Scale subgrid_data
       !

       subgrid_data = subgrid_data * mfact

       !
       ! Deallocate storage for emission data arrays and axes
       !
       deallocate(exc_state_density)
       deallocate(rec_state_density)
       deallocate(raxis)
       deallocate(zaxis)

    elseif (iselect.eq.35) then 

       !
       !----------------------------------------------------------
       !
       !     CH Emission Calculation
       !
       !     Only supports ISTATE=4 (CH) at the moment 
       !
       !
       !     HC subgrid data array
       !     1 = C+
       !     2 = C
       !     3 = CH+
       !     4 = CH
       !   ====  CURRENTLY NOT STORED ON HC SUBGRIDS ====
       !     5 = CH2+
       !     6 = CH2
       !     7 = CH3+
       !     8 = CH3
       !     9 = CH4+
       !     10= CH4
       !
       !
       !----------------------------------------------------------
       !

       if (istate.ne.4) then 

          call errmsg('SUBGRID_PLOTS: INCORRECT VALUE OF ISTATE SPECIFIED FOR CH EMISSION',istate)
          ierr = 1
          return
       endif

       !
       ! Issue warning on use of this plot - removed warning when alternate emission calculation was added
       !
       !call errmsg('SUBGRID_PLOTS:','WARNING: LOAD_HC_DATA_ARRAY: OPT 35:'//&
       !     &' SOME VALUES HARD CODED FOR HC EMISSION CALCULATIONS')


       !
       ! Load CH state density and data axes
       !

       allocate(exc_state_density(sg_rdim,sg_zdim),stat=iflag)
       allocate(raxis(sg_rdim),stat=iflag)
       allocate(zaxis(sg_zdim),stat=iflag)

       call calc_subgrid_axis(raxis,sg_rmin,sg_rmax,sg_rdim)
       call calc_subgrid_axis(zaxis,sg_zmin,sg_zmax,sg_zdim)
       call load_subgrid_xsection(exc_state_density,iselect,istate)

       if (ierr.ne.0) return 




       ! jdemod
       !
       ! Change the CH emission calculation to use the 0-0 band emission data from Fantz
       !
       !
       !    http://www.eirene.de/eigen/Sample_ex/ax/ax.html
       !
       !    U. Fantz J. Nucl. Mat. 337-339 (2005) p1087 - 1091, and: U. Fantz priv. comm. 2005
       !
       !    Fit <ex> = 1.51e-14 * exp(-3.19/Te) / Te**0.29 * 1.e6   [cm3/s]
       !               remove 1e6 factor for m3/s
       !

       !
       ! NOTE: CD data is currently state 4 in the stored 
       !       HC data. 
       !


       call calc_ch_emission(subgrid_data,exc_state_density,raxis,zaxis,sg_rdim,sg_zdim)


       !
       !        Scale by MFACT
       !
       IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC

       !
       ! Scale subgrid_data
       !

       subgrid_data = subgrid_data * mfact

       !
       ! Deallocate storage for emission data arrays and axes
       !
       deallocate(exc_state_density)
       deallocate(raxis)
       deallocate(zaxis)

    endif


  end subroutine load_subgrid_array




  subroutine calc_ch_emission(subgrid_data,exc_state_density,raxis,zaxis,nr,nz)
    implicit none
    integer :: nr,nz
    real subgrid_data(nr,nz)
    real exc_state_density(nr,nz)
    real raxis(nr),zaxis(nz)




    integer :: iraxis,izaxis
    real :: ne,te,ti,vb,ef,nh,nh_mol
    real :: fCD,Gaunt,deltaE,XCD
    logical :: outofgrid

    subgrid_data = 0.0


    !
    ! jdemod
    !
    ! Change the CH emission calculation to use the 0-0 band emission data from Fantz
    !
    !
    !    http://www.eirene.de/eigen/Sample_ex/ax/ax.html
    !
    !    U. Fantz J. Nucl. Mat. 337-339 (2005) p1087 - 1091, and: U. Fantz priv. comm. 2005
    !
    !    Fit <ex> = 1.51e-14 * exp(-3.19/Te) / Te**0.29 * 1.e6   [cm3/s]
    !               remove 1e6 factor for m3/s
    !

    do iraxis = 1,nr
       do izaxis = 1,nz

          if (exc_state_density(iraxis,izaxis).le.0.0) cycle

          call get_plasma_rz(raxis(iraxis),zaxis(izaxis),ne,te,ti,vb,ef,nh,nh_mol)

          xcd =  1.51e-14 * exp(-3.19/te) / te**0.29 

          subgrid_data(iraxis,izaxis) = ne * exc_state_density(iraxis,izaxis) * xcd

          !write(6,'(a,2i6,10g12.5)') 'CD:',iraxis,izaxis,xcd,ne,te,ti,exc_state_density,subgrid_data(iraxis,izaxis)

       end do
    end do

    return
  end subroutine calc_ch_emission




end module subgrid_plots
