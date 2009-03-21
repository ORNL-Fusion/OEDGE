c     -*-Fortran-*- 
c
c
c ======================================================================
c
c subroutine: ValidateUnstructuredInput
c
c
      SUBROUTINE ValidateUnstructuredInput
      IMPLICIT none
c 
      include 'params'
      include 'comtor'
c
c
c     No validation at present for LIM
c

      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE InitializeOUTUnstructuredInput
      IMPLICIT none
c
c     This routine sets the OUT related Unstructured inputs to their 
c     default values. 
c
c
      INCLUDE 'params'
      include 'comtor'


c
c -----------------------------------------------------------------------
c
c     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
c
c -----------------------------------------------------------------------
c
c
c -----------------------------------------------------------------------
c
c     TAG O01:
c
c     This option allows an absolute scaling factor for the LIM
c     run results to be specified in the OUT routine. It's default
c     value is zero.
c
      new_absfac = 0.0
c         
c -----------------------------------------------------------------------
c
c     TAG O99:
c
c     Net erosion plot scaling option
c     0 = normal (which is particles/m /particle entering the system)
c     1 = mm/hr   ABSFAC * 3600 / 1.22e26 for Beryllium!!
c     2 = Not yet implemented
c
      erosion_scaling_opt = 0
c      

      
      return
      end
c
c ======================================================================
c
      SUBROUTINE InitializeUnstructuredInput
      use iter_bm
      use variable_wall
      IMPLICIT none
c
c     This routine sets the Unstructured inputs to their 
c     default values. 
c
c
      INCLUDE 'params'
      include 'comtor'

c
c -----------------------------------------------------------------------
c
c     TAG D07:
c
c     Sputter data option - this option specifies which set of 
c     data will be used for calculating sputtering yields. The 
c     available options are:
c     1 - original LIM - Bohdansky
c     2 - Eckstein IPP9/82 (1993)
c     3 - Eckstein IPP9/82 + Adjustments from Garcia/Rosales-Roth 1996
c     4 - specified constant yield
c     5 - As 3 except a custom routine is used for W.  
c
c
      csputopt = 3

c
c -----------------------------------------------------------------------
c
c     TAG D08:
c
c     Chemical Sputter data option - this option specifies which set of 
c     data will be used for calculating chemical sputtering yields. The 
c     available options are:
c     1 - Garcia-Rosales/Roth 1994
c     2 - Garcia-Rosales/Roth 1996
c     3 - JET 1 - Garcia-Rosales Formula EPS94
c     4 - JET 2 - Popiesczyk EPS95
c     5 - JET 3 - Vietzke (from Phys.Processes.of.Interaction.Fus.Plas.with.Solids)
c     6 - JET 4 - Haasz - submitted to JNM Dec 1995
c     7 - JET 5 - Roth & Garcia-Rosales - submitted to JNM March 1996
c     8 - JET 6 - Haasz 1997 - Brian Mech's PhD thesis data
c     9 - Constant yield with value set to const_yield
c    10 - Modified Haasz 1997 - mass dependence made for H 
c                             - lower yield at low temperatures
c    11 - Modified Haasz 1997 - mass dependence
c                             - reduced yield at low plasma temps 
c                             - modified surface temp dependence
c
c     Set default to unmodified Haasz data 
c
      cchemopt = 8
c
c -----------------------------------------------------------------------
c
c     TAG D23:
c
c     Constant value for use with CPUTOPT option = 4 
c
      const_yield = 0.01
c     
c
c -----------------------------------------------------------------------
c
c     TAG D39 : Alternate Sputter data specifier - used to select one of 
c               several custom sputter datasets - usually based
c               on different impact angles
c
c     Set to normal incidence data as default
c 
      extra_sputter_angle = 0.0
c
c -----------------------------------------------------------------------
c
c     TAG D98 : Initial sputtered particle Y-coordinate
c               This quantity over rides the Y-coordinate generated
c               from the neutral following routines and replaces it with
c               the specified value. 
c
      init_y_coord = 0.0
c
c -----------------------------------------------------------------------
c
c     TAG D99 : Sputter impact energy option -
c               0 = LIM standard
c               1 = remove particle kinetic energy term 
c                   leaving sheath and temperature
c
c     Start counting down to avoid collision with DIVIMP tags
c
c     Set to standard as default
c
      impact_energy_opt = 0

c
c -----------------------------------------------------------------------
c
c     TAG L01: The "L" series is designated for LIM unstructured input
c
c              This option is only in effect for limiters with a
c              specified poloidal extent for 3D LIM. Basically, particles
c              crossing the L_con/2 point which are on a flux tube not
c              connecting to the limiter will have their P coordinate
c              adjusted so they will be placed on a flux tube which
c              connects to the limiter. 
c
c              Shear_short_circuit_opt 1: P = CPCO * ( 2*ran-1) = (-CPCO,+CPCO) 
c
      shear_short_circuit_opt=0
c
c -----------------------------------------------------------------------
c
c     TAG L02: LIM Wall shape option - allow CAW to vary as a function of Y
c              0 = off ... Wall = CAW
c              1..n = on      CAW  = function of Y (option specifies function)
c
c              Option 1 = linear wall from CAW at ywall_start to 
c                         CAW_MIN at YHALF (half way point between limiters)
c
      lim_wall_opt = 0
c
c     TAG L03: LIM Wall shape option - starting Y value for revised wall value
c
      ywall_start = 0.0

c
c     TAG L04: LIM Wall shape option - value of CAW reached at midpoint (YHALF)
c
      caw_min = HI
c
c -----------------------------------------------------------------------
c
c     TAG L05 to L9: Limiter shape parameters for EDGE option 11 - ITER
c     L05: rtor_setback - radial setback from LCFS at the toroidal half width of the BM
c     L06: rslot_setback - radial setback from LCFS at slot half width
c     L07: bm_tor_wid - toroidal half width of the BM (blanket module)
c     L08: slot_tor_wid - toroidal half witdth of the center slot of BM 
c     L09: lambda_design - design SOL decay length
c
      rtor_setback  = 0.07
      rslot_setback = 0.01
      bm_tor_wid    = 0.5954 
      slot_tor_wid  = 0.03
      lambda_design = 0.015
c
c -----------------------------------------------------------------------
c
c     TAG Q26:
c
c     Specification of a density multiplier (gradient) to be applied
c     to the outboard region. 
c
C     READ IN DENSITY GRADIENT INFORMATION, IF ANY
C     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
C     OF THE DENSITY
c
c     Turned off by default 
C
      nnbg = 0 
C
C
C
C -----------------------------------------------------------------------
C
C
C     TAG M01:
C
C     Set the initial velocity angle of neutrals as 0.0 unless otherwise
C     specified when CNEUTC=17
C
      CIANGN=0.0
C
c
c -----------------------------------------------------------------------
c
c
c     End of initialization 
c
      return
      end  
c
c ======================================================================
c
      SUBROUTINE ReadUnstructuredInput(line2)
      use iter_bm
      use variable_wall
      IMPLICIT none

      CHARACTER line2*(*),LINE*72,TAG*3,COMENT*72,cdum1*1024
      REAL      R,vol,z1
      INTEGER   I,ir,ierr,i1,i2
c
c jdemod - added variable to hold line read when calling RDG1 to get 
c          ADAS data.
c
      character line3*512

c
      INCLUDE 'params'
      include 'comtor'
c
c
c      COMMON /INPUTCHK/ inputflag
c      INTEGER           inputflag(100)
c
c      COMMON /MACHCOM/ machine
c      CHARACTER*64     machine
c
c      INTEGER    MAXTAG
c      PARAMETER (MAXTAG=1000)
c      COMMON /INPCHK/ ntaglist,taglist
c      INTEGER     ntaglist,idum1
c
c      CHARACTER*3 taglist(MAXTAG)
c
c     Function declaration for TAG T29 
c
c      real vtest,res,vr_pdf_int
c      external vr_pdf_int
c

      integer in
c
      WRITE(line,'(A72)') line2

      WRITE(TAG,'(A3)') LINE(3:5)

      ierr = 0

c
c -----------------------------------------------------------------------
c
c     TAG D07 : Physical Spuuter Data option
c
      IF (tag(1:3).EQ.'D07') THEN
c
c
c     Physical Sputter data option - this option specifies which set of 
c     data will be used for calculating physical sputtering yields. The 
c     available options are:
c     1 - original LIM - Bohdansky
c     2 - Eckstein IPP9/82 (1993)
c     3 - Eckstein IPP9/82 + Adjustments from Garcia/Rosales-Roth 1996
c     4 - specified constant yield
c     5 - As 3 except a custom routine is used for W.  
c
c
        CALL ReadI(line,csputopt,1,5,'Sputter Data option')
c
c
c -----------------------------------------------------------------------
c
c     TAG D08 : Chemical Spuuter Data option
c
      ELSEIF (tag(1:3).EQ.'D08') THEN
c
c
c     Chemical Sputter data option - this option specifies which set of 
c     data will be used for calculating chemical sputtering yields. The 
c     available options are:
c     1 - Garcia-Rosales/Roth 1994
c     2 - Garcia-Rosales/Roth 1996
c     3 - JET 1 - Garcia-Rosales Formula EPS94
c     4 - JET 2 - Popiesczyk EPS95
c     5 - JET 3 - Vietzke (from Phys.Processes.of.Interaction.Fus.Plas.with.Solids)
c     6 - JET 4 - Haasz - submitted to JNM Dec 1995
c     7 - JET 5 - Roth & Garcia-Rosales - submitted to JNM March 1996
c     8 - JET 6 - Haasz 1997 - Brian Mech's PhD thesis data
c     9 - Constant yield with value set to const_yield
c    10 - Modified Haasz 1997 - mass dependence made for H 
c                             - lower yield at low temperatures
c    11 - Modified Haasz 1997 - mass dependence
c                             - reduced yield at low plasma temps 
c                             - modified surface temp dependence
c
c
        CALL ReadI(line,cchemopt,1,11,'Chemical Sputter Data option')
c
c -----------------------------------------------------------------------
c
c     TAG D23 : Yield for Sputter option 4 - constant
c 
      elseif (tag(1:3).EQ.'D23') THEN
c
c     Constant value for use with CPUTOPT option = 4 
c
        CALL ReadR(line,const_yield,0.0,1.0,'Specified constant yield')

c
c -----------------------------------------------------------------------
c
c     TAG D39 : Alternate Sputter data specifier - used to select one of 
c               several custom sputter datasets - usually based
c               on different impact angles
c 
      elseif (tag(1:3).EQ.'D39') THEN
c
c     Secondary sputter data specifier
c
        CALL ReadR(line,extra_sputter_angle,-10.0,90.0,
     >             'Extra Sputter Angle Opt')
c
c
c -----------------------------------------------------------------------
c
c     TAG D98 : Initial sputtered particle Y-coordinate
c               This quantity over rides the Y-coordinate generated
c               from the neutral following routines and replaces it with
c               the specified value. 
c
      elseif (tag(1:3).EQ.'D98') THEN
c
c     Sputtered/launched particle initial Y coordinate
c
        CALL ReadR(line,init_y_coord,-HI,HI,
     >       'Specified Initial Y coordinate for all Launched Ions')
c
c -----------------------------------------------------------------------
c
c     TAG D99 : Sputter impact energy option -
c               0 = LIM standard
c               1 = remove particle kinetic energy term 
c                   leaving sheath and temperature
c
c     Start label counting down to avoid collision with DIVIMP tags
c     Set to standard as default
c
      elseif (tag(1:3).EQ.'D99') THEN
c
c     Impurity particle impact energy option
c
        CALL ReadI(line,impact_energy_opt,0,1,
     >                'Impurity impact Energy calculation option')
c
c
c
c -----------------------------------------------------------------------
c
c     TAG L01: The "L" series is designated for LIM unstructured input
c
c              This option is only in effect for limiters with a
c              specified poloidal extent for 3D LIM. Basically, particles
c              crossing the L_con/2 point which are on a flux tube not
c              connecting to the limiter will have their P coordinate
c              adjusted so they will be placed on a flux tube which
c              connects to the limiter. 
c
c              Shear_short_circuit_opt 0: OFF
c              Shear_short_circuit_opt 1: P = CPCO * ( 2*ran-1) = (-CPCO,+CPCO) 
c
      elseif (tag(1:3).EQ.'L01') THEN
        CALL ReadI(line,shear_short_circuit_opt,0,1,
     >                'Shear Short Circuit Option')
c
c -----------------------------------------------------------------------
c
c     TAG L02: LIM Wall shape option - allow CAW to vary as a function of Y
c              0 = off ... Wall = CAW
c              1..n = on      CAW  = function of Y (option specifies function)
c
c              Option 1 = linear wall from CAW at ywall_start to 
c                         CAW_MIN at YHALF (half way point between limiters)
c
      elseif (tag(1:3).EQ.'L02') THEN
        CALL ReadI(line,lim_wall_opt,0,1,'LIM wall option')
c
c     TAG L03: LIM Wall shape option - starting Y value for revised wall value
c
      elseif (tag(1:3).EQ.'L03') THEN
        CALL ReadR(line,ywall_start,0.0,HI,
     >               'Starting Y value for alternate wall')
c
c     TAG L04: LIM Wall shape option - value of CAW reached at midpoint (YHALF)
c
      elseif (tag(1:3).EQ.'L04') THEN
        CALL ReadI(line,caw_min,-HI,0.0,'Distance to wall at Yhalf')
c
c -----------------------------------------------------------------------
c
c     TAG L05 to L09: Limiter shape parameters for EDGE option 11 - ITER
c     L03: rtor_setback - radial setback from LCFS at the toroidal half width of the BM
c     L04: rslot_setback - radial setback from LCFS at slot half width
c     L05: bm_tor_wid - toroidal half width of the BM (blanket module)
c     L06: slot_tor_wid - toroidal half witdth of the center slot of BM 
c     L07: lambda_design - design SOL decay length
c
      elseif (tag(1:3).EQ.'L05') THEN
        CALL ReadR(line,rtor_setback,0.0,HI,
     >          'Radial setback at BM edge (M)')
      elseif (tag(1:3).EQ.'L06') THEN
        CALL ReadR(line,rslot_setback,0.0,HI,
     >                    'Radial setback at inner slot edge (M)')
      elseif (tag(1:3).EQ.'L07') THEN
        CALL ReadR(line,bm_tor_wid,0.0,HI,
     >                    'Toroidal Half width of BM (M)')
      elseif (tag(1:3).EQ.'L08') THEN
        CALL ReadR(line,slot_tor_wid,0.0,HI,
     >                    'Toroidal half width of slot (M)')
      elseif (tag(1:3).EQ.'L09') THEN
        CALL ReadR(line,lambda_design,0.0,HI,
     >                    'Design decay length (M)')
c
c -----------------------------------------------------------------------
c
c     TAG Q26:
c
c     Specification of a density multiplier (gradient) to be applied
c     to the outboard region. 
c
      elseif (tag(1:3).EQ.'Q26') THEN
c
c
C     READ IN DENSITY GRADIENT INFORMATION, IF ANY
C     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
C     OF THE DENSITY
C

         CALL RDRARN(MNBG,NNBG,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,            
     >                                     1,'SET OF Y,MNB VALUES',IERR)
C
C
C
C -----------------------------------------------------------------------
C
C     TAG M01:
C
C     Set initial velocity angle of neutrals as 0.0 unless otherwise
C     specified when CNEUTC=17
C
      ELSEIF (TAG(1:3).EQ.'M01') THEN
        CALL ReadR(line,CIANGN,-180.0,180.0,'Initial neutral velocity
     > angle')       
        ciangn = ciangn * degrad  
c
c
c -----------------------------------------------------------------------
c
c     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
c
c -----------------------------------------------------------------------
c
c     TAG O01: 
c
c     Specify an alternate absolute factor in LIM
c
      elseif (tag(1:3).eq.'O01') then 
c
c     This option allows an absolute scaling factor for the LIM
c     run results to be specified in the OUT routine. It's default
c     value is zero.
c
        CALL ReadR(line,new_absfac,0.0,HI,
     >                   'Imposed ABSFAC in OUT')
c
c -----------------------------------------------------------------------
c
c     TAG O99:
c
c     Net erosion plot scaling option
c     0 = normal (which is particles/m /particle entering the system)
c     1 = mm/hr   ABSFAC * 3600 / 1.22e26 for Beryllium!!
c     2 = Not yet implemented
c
      elseif (tag(1:3).eq.'O99') then 
c      
        CALL ReadI(line,erosion_scaling_opt,0,2,
     >                'Erosion Scaling Option')
c         
c
c -----------------------------------------------------------------------
c
c     TAG does not match available input - signal ERROR and EXIT
c

      ELSE
        CALL ER('ReadUnstructuredInput','Unrecognized tag',*99)
      endif
c

      RETURN
c
c There is an error:
c
99    WRITE(6,'(5X,3A)') 'LINE = "',line,'"'
      WRITE(6,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
      WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0,*) '    DIVIMP HALTED'
      STOP
      END


c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadIR(line,ival,rval,imin,imax,tag)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival,imin,imax
      REAL    rval

      INTEGER i
      REAL    r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i,r

      IF (i.LT.imin.OR.i.GT.imax)
     .  CALL ER('ReadI','Out of bounds: '//line,*99)

      ival = i
      rval = r

      WRITE(DBGUNIT,'(A)') line
      WRITE(DBGUNIT,'(5X,2A,I4,1P,E10.2)') tag,' = ',ival,rval

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadI(line,ival,imin,imax,tag)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival,imin,imax

      INTEGER i
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i

      IF (i.LT.imin.OR.i.GT.imax) then 

        write (0,*)  'READI:ERROR:',i,imin,imax 
        CALL ER('ReadI','Out of bounds: '//line,*99)

      endif

      ival = i

      WRITE(DBGUNIT,'(A)')        line
      WRITE(DBGUNIT,'(5X,2A,I4)') tag,' = ',ival

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadC(line,cval,tag)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*(*),tag*(*),cval*(*)
      INTEGER fp,ival,imin,imax

      INTEGER i
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,cval

      WRITE(DBGUNIT,'(A)')        line
      WRITE(DBGUNIT,'(5X,2A,A)') tag,' = ',cval

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,2A)')    'CVAL = ''',cval,''''
      STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE Read2I(line,ival1,ival2,imin,imax,tag)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival1,ival2,imin,imax

      INTEGER i1,i2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i1,i2

      IF (i1.LT.imin.OR.i1.GT.imax.OR.
     .    i2.LT.imin.OR.i2.GT.imax)
     .  CALL ER('Read2I','Out of bounds: '//line,*99)

      ival1 = i1
      ival2 = i2

      WRITE(DBGUNIT,'(A)')        line
      WRITE(DBGUNIT,'(5X,2A,I4)') tag,' = ',ival1
      WRITE(DBGUNIT,'(5X,2A,I4)') tag,' = ',ival2

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE ReadR(line,rval,rmin,rmax,tag)
      use error_handling
      IMPLICIT none

      CHARACTER line*72,tag*(*)
      REAL rval,rmin,rmax

      INCLUDE 'params'
      INCLUDE 'slcom'

      REAL r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r

      IF (r.LT.rmin.OR.r.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval = r

      WRITE(DBGUNIT,'(A)')        line
      WRITE(DBGUNIT,'(2A,G10.3)') tag,' = ',rval

      RETURN

 98   call errmsg('READR','Problem reading unstructured input')
      WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,3G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r,rval,rmin,rmax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE Read2R(line,rval1,rval2,rmin,rmax,tag)

      IMPLICIT none

      CHARACTER line*72,tag*(*)
      REAL rval1,rval2,rmin,rmax

      INCLUDE 'params'
      INCLUDE 'slcom'

      REAL r1,r2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r1,r2

      IF (r1.LT.rmin.OR.r1.GT.rmax.OR.
     .    r2.LT.rmin.OR.r2.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval1 = r1
      rval2 = r2

      WRITE(DBGUNIT,'(A)')        line
      WRITE(DBGUNIT,'(2A,2G10.3)') tag,' = ',rval1,rval2

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,6G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r1,r2,rval1,rval2,rmin,rmax
      STOP
      END







