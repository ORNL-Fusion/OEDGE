module unstructured_input

  use debug_options
  use error_handling
  implicit none

  private

  integer,parameter,private :: maxlen=1024

  INTEGER :: tagnum=0,fp
  INTEGER,parameter ::    MAXTAG=1000
  INTEGER     ntaglist
  CHARACTER*3 taglist(MAXTAG)

  !
  !     this common block contains variables used in out and read in using the
  !     unstructured input methodology
  !
  !     -*-fortran-*-
  ! common /out_unstructured_input/ new_absfac,psi1_reg,psi2_reg
  ! save /out_unstructured_input/
  real,public :: new_absfac
  integer,public :: absfac_opt,e2dizs_offset  
  real ,public :: psi1_reg,psi2_reg
  integer,public :: absfac_iz, absfac_ir, absfac_ikstart, absfac_ikend
  integer,public :: scale_1d

  public :: allocate_unstructured_input,deallocate_unstructured_input,init_out_unstruc_input, ReadUnstructuredInput,scan_input_file
  real,public,allocatable :: transport_area(:,:)

contains



  subroutine scan_input_file
    use mod_params
    use mod_io_units
    use mod_reader
    use mod_io
    implicit none

    ! Update some parameters based on input file contents
    ! Also attempt to determine whether the input file is structured or tagged

    character*(maxlen) :: line
    character*3 :: tag
    character*2 :: tag_tn
    character*1 :: tag_mod

    integer :: ierr
    logical :: err

    ! jdemod
    ! - access the input file
    ! - read in any updates to the parameter values from their default values
    ! - determine the type of input file - tagged or untagged 
    ! parameter tags start with the # sign and are in series Z
    !
    
    ierr = 0
    err = .false.
    
    do while (.not.err)

       ! loop through input file

       call divrd(ierr,line)
      
       if (ierr.eq.0) then
          ! check for parameter tags - series Z
          
          WRITE(TAG,'(A3)') LINE(3:5) 
          tag_mod = line(2:2)
          tag_tn  = line(2:3)
          ! jdemod - add support for case insensitive tags
          call upcase(tag)
          call upcase(tag_tn)
          
          ! possible parameter over ride found - parameters use series Z tags
          if (tag_mod == '#') then 

             ! OUT reads most parameters from the raw file. However, maxdatx is used in OUT to define the maximum
             ! number of experimental data points in mod_expt_data
             if (tag(1:3).eq.'Z07') then 
                ! Z07 : MAXDATX - Max number of experimental data points allowed to be loaded
                CALL divrd(maxdatx,.TRUE.,500,.false.,0,'Max number of experimental data points allowed to be loaded',IERR)
                
             endif
          endif
       else
          err = .true.
       endif
          
    end do

    ! rewind input file
    call divrd

    
  end subroutine scan_input_file


  
  subroutine allocate_unstructured_input
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_out_unstruc','ALLOCATE')
    call allocate_array(transport_area,maxnks,maxnrs,'transport_area',ierr)

  end subroutine allocate_unstructured_input


  subroutine deallocate_unstructured_input
    implicit none

    call pr_trace('mod_out_unstruc','DEALLOCATE')
    if (allocated(transport_area)) deallocate(transport_area)

  end subroutine deallocate_unstructured_input



  subroutine init_out_unstruc_input
    use mod_params
    use mod_comtor
    implicit none
    ! 
    !     This routine assigns the default values to any
    !     OUT unstructured input values.The OUT unstructured
    !     input tags start with the letter "O" (oh)
    !
    !     jdemod - it also assigns default values to any other 
    !              unstructured input values shared by DIV and OUT
    !
    !
    !     Note: The netcdf output option A07 from DIVIMP is also supported
    !           since it would be useful to create a netcdf version of the 
    !           raw output from an OUT run
    !
    !     include 'params'
    !
    !     include 'out_unstruc'
    !     include 'comtor'
    !     
    !------------------------------------------------------
    !
    !     O01 - alternate absolute factor specification 
    !         - this option is deactivated by setting the 
    !           default value to zero.  
    !
    new_absfac=0.0
    !
    !
    !------------------------------------------------------
    !
    !     Core fueling code calculates integrated ionization
    !     profiles in the core ... these parameters allow the 
    !     PSIN inner bound of the integration regions to be set
    !     This is used in the pr_eirene_analysis routine
    !
    !     O02 - PSIN bound for calculating core ionization 
    !           profile 1 (psi1_reg)
    !     O03 - PSIN bound for calculating core ionization 
    !           profile 2 (psi2_reg)
    !

    psi1_reg = 0.9
    psi2_reg = 0.95


    ! -----------------------------------------------------------------------
    !
    !     O04:
    !     absfac_opt - options specifying how to calculate scaling factor for DIVIMP
    !                  results in OUT         
    !
    !     This option allows the absolute scaling factor for the DIVIMP
    !     run results to be specified in the OUT routine. It's default
    !     value is zero.
    !      

    absfac_opt = 0

    ! -----------------------------------------------------------------------
    !
    !     O05:
    !     e2dizs_offset - offset to match fluid code impurity charge state results data to 
    !                     DIVIMP results (needed sometimes since the fluid code results
    !                     sometimes contain multiple fluids). 
    !

    e2dizs_offset = 0


    ! -----------------------------------------------------------------------
    !
    !     O06->O09:
    !
    !     Indices to control the ABSFAC calculation
    !     Default values are initially zero. These will be reset to specific default
    !     values depending on the absfac_opt option when the code is executed
    !

    absfac_iz = 0       ! goes to maxizs
    absfac_ir = 0       ! goes to irsep-1
    absfac_ikstart = 0  ! goes to 1
    absfac_ikend = 0    ! goes to nks(absfac_ir) or nks(absfac_ir)-1 inside separatrix

    ! -----------------------------------------------------------------------
    !
    !  O10: Scale plots for 1D (convert from kareas to distance along the field line
    !       0 = off,
    !       1 = on (scale by distance along the field line)
    !       2 = on (scale by transport cell area instead of mesh cell area)
    scale_1d = 0


    !
    ! -----------------------------------------------------------------------
    !
    !     TAG A07: 
    !
    !
    !     Option to write a netcdf version of the raw data file 
    !     0 = off   1 = 0n .. default OFF
    !
    !
    netcdf_opt = 0
    !
    !
    !------------------------------------------------------
    !
    return
  end subroutine init_out_unstruc_input


  subroutine ReadUnstructuredInput(line2)

    use mod_params
    use mod_io
    use mod_io_units
    IMPLICIT none

    CHARACTER line2*(*),LINE*128,TAG*3,cdum1*1024
    REAL      R,vol,z1,version
    INTEGER   I,ir,ierr,i1,i2


    INTEGER    MAXTAG
    PARAMETER (MAXTAG=1000)
    COMMON /INPCHK/ ntaglist,taglist
    INTEGER     ntaglist,idum1
    REAL        rdum1
    CHARACTER*3 taglist(MAXTAG)

    INTEGER    itag
    LOGICAL    status,sol28_first
    CHARACTER  local_buffer*1024
    external upcase

    DATA sol28_first /.TRUE./
    SAVE

    WRITE(line,'(A128)') line2

    WRITE(TAG,'(A3)') LINE(3:5)

    call  upcase(tag)

    ierr = 0

    ! jdemod - fp is not used in OUT - set to ipinout instead of pinout from divimp version of mod_slcom
    ! revisit later if needed
    !fp = PINOUT
    fp = IPINOUT

    ntaglist = ntaglist + 1
    taglist(ntaglist) = tag


    !
    ! Series O
    !     - this tag refers to unstructured input related to OUT
    !     - OUT shares this code with DIVIMP but the initialization routines 
    !       are separate. This does cause some overhead storing variables not 
    !       used in the specific codes but for now this is only 3 reals in 
    !       OUT or about 12 bytes ... revisit later if it becomes an issue
    !
    if (tag(1:1).eq.'O') then 
       call read_out_unstructured_input(line,tag,fp)
    else
       call errmsg('MOD_OUT_UNSTRUC: ReadUnstructuredInput','INVALID TAG ='//trim(tag))
    endif

    !...  Check tag list:
    DO i1 = 1, ntaglist-1
       IF (tag.EQ.taglist(i1)) THEN
          call errmsg('ReadUnstructuredInput','Duplicate tag detected (Program Stopping):'//trim(tag))
          WRITE(0,*) 'DUPLICATE TAG: "'//tag//'"'
          STOP 'OUT STOP: DUPLICATE TAG'//trim(tag)
       ENDIF
    ENDDO


    return

  end subroutine ReadUnstructuredInput

  subroutine read_out_unstructured_input(line,tag,fp)
    use mod_params
    use mod_slcom
    use mod_io
    implicit none
    !
    !     READ "O" Series Unstructured input
    !
    !     The Oh is used for OUT related tagged input
    !
    INTEGER   fp
    CHARACTER line*(*),tag*3

    !
    ! -----------------------------------------------------------------------
    !
    !     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
    !
    ! -----------------------------------------------------------------------
    !

    if (tag(1:3).eq.'O01') then 

       !     new_absfac - use to change scaling of plots from OUT

       !     This option allows the absolute scaling factor for the DIVIMP
       !     run results to be specified in the OUT routine. It's default
       !     value is zero.


       CALL ReadR(line,new_absfac,0.0,HI,'Imposed ABSFAC in OUT')

    elseif (tag(1:3).eq.'O02') then 
       !        
       !     Core fueling code calculates integrated ionization
       !     profiles in the core ... these parameters allow the 
       !     PSIN inner bound of the integration regions to be set
       !     This is used in the pr_eirene_analysis routine
       !     

       !        O02 - PSIN bound for calculating core ionization 
       !              profile 1 (psi1_reg)

       CALL ReadR(line,psi1_reg,0.0,HI,'PSIN bound for core ionization profile 1')

    elseif (tag(1:3).eq.'O03') then 

       !        O03 - PSIN bound for calculating core ionization 
       !              profile 2 (psi2_reg)

       CALL ReadR(line,psi2_reg,0.0,HI,'PSIN bound for core ionization profile 2')

    elseif (tag(1:3).eq.'O04') then 

       !     absfac_opt - options specifying how to calculate scaling factor for DIVIMP
       !                  results in OUT         

       !     This option allows the absolute scaling factor for the DIVIMP
       !     run results to be specified in the OUT routine. It's default
       !     value is zero.


       CALL ReadI(line,absfac_opt,0,3,'ABSFAC calculation option in OUT')

    elseif (tag(1:3).eq.'O05') then 

       !     e2dizs_offset - offset to match fluid code impurity charge state results data to 
       !                     DIVIMP results (needed sometimes since the fluid code results
       !                     sometimes contain multiple fluids). 

       CALL ReadI(line,e2dizs_offset,0,100,'FC Impurity offset index')

    elseif (tag(1:3).eq.'O06') then 

       !     absfac_iz - 

       CALL ReadI(line,absfac_iz,0,maxizs+1,'FC ABSFAC IZ (no offset)')

    elseif (tag(1:3).eq.'O07') then 

       !     absfac_ir - 

       CALL ReadI(line,absfac_ir,1,maxnrs,'FC ABSFAC IR')

    elseif (tag(1:3).eq.'O08') then 

       !     absfac_ikstart - 

       CALL ReadI(line,absfac_ikstart,1,maxnks,'FC ABSFAC IK-START')

    elseif (tag(1:3).eq.'O09') then 

       !     absfac_ikend - 

       CALL ReadI(line,absfac_ikend,1,maxnks,'FC ABSFAC IK-END')

    elseif (tag(1:3).eq.'O10') then 

       !     scale_1d: Scale generalized results plots from load_divdata by using the
       !               cell length along the field line instead of the cell volume.
       !               This applies to code density results and others normalized
       !               by cell area.          

       CALL ReadI(line,scale_1d,0,2,'SCALE 1D')

    ELSE
       CALL ER('ReadUnstructuredInput','Unrecognized tag',*99)
    ENDIF

    return

99  WRITE(SLOUT,'(5X,3A)') 'LINE = "',line,'"'
    WRITE(SLOUT,'(5X,3A)') 'TAG  = "',tag ,'"'
    WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
    WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
    WRITE(0,*) '    DIVIMP HALTED'
    STOP
  END subroutine read_out_unstructured_input



end module unstructured_input
