module mod_soledge

  use mod_params

  implicit none
  private
  

  real*8 :: lensrc, lamsrc
  real*8 :: s0out, s0in,s0aout, s0ain,s0bout, s0bin


  integer :: maxn = 1000

  real*8,allocatable :: te(:), ti(:), ne(:),vb(:),ef(:),sd(:),teg(:),tig(:),ga(:)
  real, allocatable :: yd(:)

  real*8 :: soffset

  real*8 :: ymin, ymax



    !     CSOLLS    - DECAY COEFFICIENT FOR SOURCE
    !     CSOLPR    - INITIAL RADIATION COEFFICIENT
    !     CSOLLR    - DECAY COEFFICIENT FOR RADIATION
  
  real*8 :: ionint, ioninti, radint, radinti
  real*8 :: sionl, iionl
  real*8 :: sradl, iradl

  integer :: pionl, pradl
  
  integer,public :: colprobe3d
  
  integer:: fluxropt = 0
  integer:: srootopt = 0
  integer:: cpowopt  = 1

  real*8 :: plensrc, plamsrc,p0in,p0out

  ! jdemod - option to turn on/off use of SOL 12,13 etc
  integer,public :: soledge_opt
  
  integer,public :: cioptf_soledge


  integer,public :: csopt != 0
  integer,public :: cpopt != 0

  real,public :: csolls != 0.5
  real,public :: csollt != 0.5

  real,public :: csollr != 0.5  ! length of radiation source
  
  real,public :: csolfr != 0.0  ! strength of radiation as fraction of target power  (Popt 2 and 3)
  real,public :: csolpr != 0.0  ! strength of radiation source in absolute terms (popt 0 and 1) 
  
  real,public :: cfiz !=  0.0   ! fraction split between two ionization sources in sopt 4 and 5
  real*8 :: cfsrc = 1.0  ! leave for later - related to imposing over-ionization
  
  real*8 :: ck0  = 2e3
  real*8 :: ck0i = 58.9

  real,public :: sol13_pdist != 0.0
  real,public :: sol13_padd  != 0.0
 
  
  !
  ! Hard code some of the SOL option parameters for now. 
  !

  public :: soledge,init_soledge



contains

  subroutine init_soledge(y1,y2)
    implicit none
    real :: y1,y2

    if (y2.gt.y1) then
       ymax = y2
       ymin = y1
    else
       ymax = y1
       ymin = y2
    endif


    !
    ! Allocate storage for each field line of plasma data
    !
    allocate(te(maxn))
    allocate(ti(maxn))
    allocate(teg(maxn))
    allocate(tig(maxn))
    allocate(ne(maxn))
    allocate(vb(maxn))
    allocate(ef(maxn))
    allocate(sd(maxn))
    allocate(yd(maxn))
    allocate(ga(maxn))



  end subroutine init_soledge

  subroutine end_soledge
    implicit none

    deallocate(te)
    deallocate(ti)
    deallocate(teg)
    deallocate(tig)
    deallocate(ne)
    deallocate(vb)
    deallocate(ef)
    deallocate(sd)
    deallocate(yd)
    deallocate(ga)

  end subroutine end_soledge





  SUBROUTINE SETUPVAL(SOPT,POPT,FSRC,LNSRC,LMSRC,SMAX,&
       N0,V0,T0,N0I,V0I,T0I,RCF,RCFI,LSSIZ,LPSIZ,IX,&
       PAOUT,PAIN)

    use mod_params
    
    ! sazmod - need these for the variables I think its okay to add them?
    !use mod_comtor
    !use mod_comxyt
    !use mod_comt2
    !use yreflection
    
    
    IMPLICIT NONE
    INTEGER SOPT,POPT,IX
    DOUBLE PRECISION FSRC,LNSRC,LMSRC,SMAX
    DOUBLE PRECISION N0,V0,T0
    DOUBLE PRECISION N0I,V0I,T0I
    DOUBLE PRECISION RCF,RCFI
    DOUBLE PRECISION LSSIZ,LPSIZ
    DOUBLE PRECISION PAOUT,PAIN
	integer ipos
	
    !
    !
    !
    !     SETUPVAL: THIS ROUTINE LOADS THE COMMON BLOCK VALUES NEEDED TO CAL
    !               THE IONIZATION AND RADIATION SOURCE VALUES. THIS ALSO IN
    !               THE EFFECT OF FLOW RECIRCULATION - IF ANY.
    !
    !     SOPT      - SOURCE IONIZATION OPTION
    !     POPT      - SOURCE RADIATION OPTION
    !     N0        - DENSITY AT PLATE (FOR THIS RING)
    !     V0        - VELOCITY AT PLATE (FOR THIS RING)
    !     T0        - TEMPERATURE AT PLATE (FOR THIS RING)
    !     CSOLLS    - DECAY COEFFICIENT FOR SOURCE
    !     CSOLPR    - INITIAL RADIATION COEFFICIENT
    !     CSOLLR    - DECAY COEFFICIENT FOR RADIATION
    !
    !
    !     DAVID ELDER   OCT  1991
    !
    !
    !     OTHER DECLARATIONS
    !
    INTEGER I,J,K,NS,in
    DOUBLE PRECISION TVAL,S0,STEMP,DS
    double precision cisint

    !real*8, EXTERNAL :: CIS1,SRCION,srcrad

	! Why can't we do the vary_absorb adjustment to smax here? Because
	! at the end of the day we will have still a 100 length array, that
	! instead of going from say -10 to 10 will go from -5 to 5. In effect
	! we stretch out our shortened field line to still cover the entire
	! array, which since the array represents the actual length of the
	! plasma we are doing, this doesn't make sense. It's almost like
	! we say 5 = 10, which is not a good thing to do.
	
	! Delete the below when I get things working.
	
	!if (vary_absorb.eq.1) then
	  ! sazmod: want to find the correct smax. It would be smaller if
	  ! say you were in the part where the step in the wall shortens
	  ! the connection length.
	  !ix_step1 = ipos(xabsorb1a_step, xs, nxs-1)
	  !if (ix.le.ix_step1) then
	    ! you are in the step region with shorter conn length. This
	    ! assumes yabsorb2a > yabsorb1a (so 2a is the positive value).
	    !smax = yabsorb2a_step - yabsorb1a_step
	    !write(0,*) 'in step smax = ',smax
	  !else
	    ! you are in the normal region. Could save space and just use
	    ! the code below this, but I want to make it explicit what 
	    ! I'm doing here.
	    !smax = yabsorb1a - yabsorb2a
	    !write(0,*) 'out step smax = ',smax
	  !endif
	!else
      smax = ymax-ymin
	!endif
	
    soffset = ymin

    ! set up s axis
    ! Yaxis is s(in) - ymin
    do in = 1, maxn
    !write(0,*) 'checkpoint6 in = ',in,'/',maxn
       sd(in) = (in-1) * smax/(maxn-1)
       !write(0,*) 'SD:',in,sd(in)
    end do
	!write(0,*) 'checkpoint5'
    yd = sd + soffset

    !do in = 1, maxn
    !   write(0,*) 'SD:',in,sd(in),yd(in)
    !end do

    !write(0,*) 'Setup opts:', sopt,popt
    !
    !
    !     SET UP IONIZATION ARRAY FROM 0 TO SMAX FOR SPECIFIC RING
    !     USING THE GIVEN OPTION FLAG
    !
    !

    IF (SOPT.EQ.0) THEN
       LENSRC = LNSRC * SMAX
       LAMSRC = LMSRC * SMAX
       S0OUT = -FSRC*N0*V0/LENSRC
       S0IN  = -FSRC*N0I*V0I/LENSRC
    ELSEIF (SOPT.EQ.1) THEN
       LENSRC = LNSRC * SMAX
       LAMSRC = LMSRC * SMAX
       S0OUT = -FSRC*N0*V0/(LAMSRC* (1-EXP(-LENSRC/LAMSRC)))
       S0IN  = -FSRC*N0I*V0I/(LAMSRC* (1-EXP(-LENSRC/LAMSRC)))
    ELSEIF (SOPT.EQ.4) THEN
       LENSRC = SMAX* LNSRC
       LAMSRC = SMAX* LMSRC
       S0AOUT = -(1.0-CFIZ)*FSRC*N0*V0/LAMSRC - FSRC*CFIZ*N0*V0/LENSRC
       S0BOUT = -FSRC*CFIZ*N0*V0/LENSRC
       S0AIN  = -(1.0-CFIZ)*FSRC*N0I*V0I/LAMSRC - FSRC*CFIZ*N0I*V0I/LENSRC
       S0BIN = -FSRC*CFIZ*N0I*V0I/LENSRC
    ELSEIF (SOPT.EQ.5) THEN
       LENSRC = LNSRC * SMAX
       LAMSRC = LMSRC * SMAX
       S0AOUT = -FSRC*(1.0-CFIZ)*N0*V0/(LAMSRC* (1.0-EXP(-LENSRC/LAMSRC)))
       S0BOUT = -FSRC*CFIZ*N0*V0/LENSRC
       S0AIN  = -FSRC*(1.0-CFIZ)*N0I*V0I/(LAMSRC* (1.0-EXP(-LENSRC/LAMSRC)))
       S0BIN = -FSRC*CFIZ*N0I*V0I/LENSRC
    ENDIF

    LSSIZ = LENSRC

    IF (POPT.EQ.0) THEN
       PLENSRC = CSOLLR * SMAX
       PLAMSRC = 0.0
       P0OUT = CSOLPR
       P0IN  = CSOLPR
    ELSEIF (POPT.EQ.1) THEN
       PLENSRC = SMAX/2.0
       PLAMSRC = CSOLLR * SMAX
       P0OUT = CSOLPR
       P0IN  = CSOLPR
    ELSEIF (POPT.EQ.2) THEN
       PLENSRC = CSOLLR * SMAX
       PLAMSRC = 0.0
       P0OUT =  CSOLFR * PAOUT / PLENSRC
       P0IN = CSOLFR * PAIN / PLENSRC
    ELSEIF (POPT.EQ.3) THEN
       PLENSRC = SMAX/2.0
       PLAMSRC = CSOLLR * SMAX
       P0OUT = CSOLFR * PAOUT / (PLAMSRC*(1.0D0-DEXP(-PLENSRC/PLAMSRC)))
       P0IN  = CSOLFR * PAIN / (PLAMSRC*(1.0D0-DEXP(-PLENSRC/PLAMSRC)))
    ENDIF
    LPSIZ = PLENSRC


    !
    !     CALCULATE CROSS FIELD TERMS IF FLUX RECIRCULATION HAS BEEN
    !     SPECIFIED. DO NOT CALCULATE IT IF A NORMALIZED IONIZATION SOURCE
    !     HAS BEEN SPECIFIED.
    !
    !IF (FLUXROPT.EQ.1.AND.FSRC.NE.1.0.AND.SOPT.NE.3) THEN
    !   RCF = -(N0*V0 + CIS1(SMAX/2.0,SRCION,SOPT,0) )
    !   >        /(SMAX/2.0)
    !   RCFI = -(N0I*V0I + CIS1(SMAX/2.0,SRCION,SOPT,1))
    !   >        /(SMAX/2.0)
    !ELSE
       RCF = 0.0
       RCFI = 0.0
    !ENDIF
    !
    !     Calculate full source integrals
    !
    !write(0,*) 'checkpoint7.1'
    ionint  = cis1(lensrc,srcion,sopt,0)
    !write(0,*) 'checkpoint7.2'
    ioninti = cis1(lensrc,srcion,sopt,1)
    !write(0,*) 'checkpoint7.3'
    radint  = cis1(plensrc,srcrad,popt,0)
    !write(0,*) 'checkpoint7.4'
    radinti = cis1(plensrc,srcrad,popt,1)
    !
    !     Set-up the values needed for source function
    !     integration.
    !
    !     Ion source
    !
    sionl = 0.0d0
    pionl = 0
    iionl = 0.0d0
    !
    !     Rad source
    !
    sradl = 0.0d0
    pradl = 0
    iradl = 0.0d0
    !
    !
    !
    !      write(6,*) 'setup:end:',SMAX,SOPT,POPT,FSRC,LNSRC,LMSRC,
    !     >           N0,V0,T0,N0I,V0I,T0I,RCF,RCFI,LSSIZ,LPSIZ,IR,
    !     >           PAOUT,PAIN
    !      write(6,*) 'setup:end2:',p0out,p0in,plensrc,plamsrc
    !
    !
    RETURN
  END SUBROUTINE SETUPVAL




  DOUBLE PRECISION FUNCTION CIS1(S,FVAL,SOPT,PLATEOPT)
    IMPLICIT NONE
    INTEGER PLATEOPT,SOPT
    DOUBLE PRECISION S,SMAX,FVAL
    EXTERNAL FVAL
    !
    !     CIS1 - PERFORMS ONE DIMENSIONAL NUMERICAL INTEGRATION ON FUNCTION
    !            FVAL. THE RANGE OF INTEGRATION IS FROM 0 TO S.
    !
    !            THIS IS CALCULATED USING THE FORMULA:
    !
    !            SIGMA(I=1TON) SI * DELTAS
    !
    !            WHERE N IS THE NTH BIN AND CORRESPONDS TO POSITION S.
    !
    !            SOPT IS THE SOURCE OPTION THAT NEEDS TO BE
    !            PASSED TO THE ROUTINE FVAL - IN ORDER TO EVALUATE THE
    !            VALUE OF EITHER IONIZATION OR RADIATION CORRECTLY.
    !
    !
    !     NSVAL -  SPECIFIES THE NUMBER OF SEGMENTS WHICH 0 TO S WILL BE
    !     DIVIDED INTO FOR THE INTEGRATION. CAN BE MADE ARBITRARILY LARGE
    !     FOR INCREASED ACCURACY AT THE COST OF PERFORMANCE.
    !
    !     THE SOURCE FUNCTION IS KNOWN TO BE IDENTICALLY ZERO BEYOND THE
    !     POINT LSIZ.
    !
    INTEGER NSVAL
    PARAMETER (NSVAL=1000)
    DOUBLE PRECISION DS,REM,INT1,STMP,F1,F2,RESULT,tmpval
    double precision strt,sret
    INTEGER NS,START,STOP,STEP,I,IND
    INTEGER INTTYPE
    !
    !
    !
    !      write(6,*) 'cis1:',s,sopt,plateopt 
    IF (S.EQ.0.0) THEN
       CIS1 = 0.0D0
       RETURN
    ENDIF
    !
    !     This call tests if the integration range is greater
    !     than the source length ... if it is a pre-calculated
    !     value for the whole source is assigned and returned
    !     Otherwise it falls through to the whole integration routine.
    !
    IND = 1
    CIS1 = FVAL(S,SOPT,PLATEOPT,S,IND)
    IF (IND.EQ.-1) RETURN
    !
    !     This call retrieves the previous value of the integral
    !     if it is appropriate.
    !
    IND = 2
    tmpval = fval(s,sopt,plateopt,sret,ind)
    if (ind.eq.-2) then
       strt = sret
    else
       strt = 0.0d0
       tmpval = 0.0d0
    endif
    !
    !      write(6,*) 'before q:',strt,s,tmpval
    INTTYPE = 1
    !write(0,*) 'checkpoint8'
    CALL QSIMP(FVAL,strt,S,RESULT, SOPT,PLATEOPT,INTTYPE)
    CIS1 = RESULT + tmpval
    
    !
    !     Update the record in the source function with the
    !     value of the integral up to this point.
    !
    tmpval = result + tmpval
    ind = 3
    tmpval = FVAL(S,SOPT,PLATEOPT,tmpval,IND)
    !
    !      WRITE(6,*) 'COMPARISON:',RESULT,INT1
    !
    !      WRITE(6,*) 'LSIZ,STMP,DS:', LSIZ,STMP,DS,INT1
    !
    !
    RETURN
  END FUNCTION CIS1



  DOUBLE PRECISION FUNCTION CIS2(S,FVAL,SOPT,PLATEOPT)
    IMPLICIT NONE
    INTEGER SOPT,PLATEOPT
    DOUBLE PRECISION S,FVAL
    EXTERNAL FVAL
    !
    !
    !     CIS2 - PERFORM A DOUBLE INTEGRATION OVER THE GIVEN SOURCE
    !            FUNCTION.
    !
    !            SOPT IS THE SOURCE OPTION THAT NEEDS TO BE
    !            PASSED TO THE ROUTINE FVAL - IN ORDER TO EVALUATE THE
    !            VALUE OF EITHER IONIZATION OR RADIATION CORRECTLY.
    !
    !     NSVAL - SPECIFIES THE NUMBER OF PARTITIONS CONTRIBUTING TO THE
    !     INTEGRATION.
    !
    INTEGER NSVAL
    PARAMETER (NSVAL=1000)
    INTEGER I,NS,J,START,STOP,STEP
    DOUBLE PRECISION DS,INT2,INT1,STMP,DS2,F1,F2,RESULT
    INTEGER INTTYPE
    !
    !
    IF (S.EQ.0.0) THEN
       CIS2 = 0.0D0
       RETURN
    ENDIF
    !
    !
    !     ADD CODE TO USE QUADRATURE FOR SECOND INTEGRAL AS WELL...
    !
    !     FOR N = 2    I = INT (X-T) F(T) DT  FOR T ON 0 TO X
    !     SRCRAD RETURNS F(T) ... FUNCTION NEEDS TO BE (X-T)F(T)
    !
    !     NEED TO CREATE A FUNCTION WHICH RETURNS (X-T)F(T) ... OR ANOTHER
    !     OPTION TO SRCRAD - INTTYPE
    !
    INTTYPE = 2
    CALL QSIMP(FVAL,0.0D0,S,RESULT, SOPT,PLATEOPT,INTTYPE)
    CIS2 = RESULT
    !
    !      WRITE(6,*) 'COMPARISON:',RESULT
    !
    !
    !
    RETURN
  END FUNCTION CIS2



  SUBROUTINE QSIMP(FUNC,A,B,S,OPT,PLATE,ITYP)
    implicit none
    !      INCLUDE 'params'
    !      INCLUDE 'slcom'

    INTEGER JMAX,OPT,PLATE,ITYP
    DOUBLE PRECISION A,B,FUNC,S,EPS
    EXTERNAL FUNC
    PARAMETER (EPS=1.0D-5,JMAX=22)
    !
    !     QSIMP: THIS ROUTINE IS TAKEN FROM THE TEXT NUMERICAL RECIPES
    !     IN FORTRAN. IT IMPLEMENTS SIMPSON'S METHOD OF CALCULATING
    !     NUMERICAL QUADRATURE. IT CALLS THE ROUTINE TRAPZD TO REFINE
    !     THE VALUE OF THE NUMERICAL INTEGRAL.
    !
    !     RETURNS AS S THE INTEGRAL OF THE FUNCTION FUNC FROM A TO B.
    !     THE PARAMETER EPS CAN BE SET TO THE DESIRED FRACTIONAL ACCURACY
    !     AND JMAX SO THAT 2 TO THE POWER JMAX-1 IS THE MAXIMUM ALLOWED NUMB
    !     OF STEPS. INTEGRATION IS PERFORMED BY SIMPSONS RULE.
    !
    !     There are two types of integration performed by these
    !     routines and the choice is controlled by the ityp
    !     argument.
    !
    !     ITYP=1 : standard integration of f(x) from a to b
    !     ITYP=2 : integration of (b-x)f(x) from a to b
    !
    INTEGER J
    DOUBLE PRECISION OS,OST,ST
    OST = -1.0D30
    OS  = -1.0D30
    DO J = 1,JMAX
       CALL TRAPZD (FUNC,A,B,ST,J,OPT,PLATE,ITYP)
       S = (4.0*ST-OST)/3.0
       IF (ABS(S-OS).LE.(EPS*ABS(OS))) RETURN

       OS = S
       OST = ST
       !
       !        WRITE(6,*) 'QSIMP:',OST,OS,ST
       !
       !WRITE(0,*) 'checkpoint9 J = ',J,'/',JMAX
    end do
    !
    !     ERROR CONDITION - INCOMPLETE CONVERGENCE
    !
    WRITE(6,*) 'QSIMP: ERROR IN QUADRATURE',A,B,ST
    WRITE(6,*) 'QSIMP:',OST,OS,ST ,S
    write(6,*) 'ALSO :',A,B,J,OPT,PLATE,ITYP

	write(0,*) 'Error: QSIMP convergence or something like that'
    STOP

  END SUBROUTINE QSIMP



  SUBROUTINE TRAPZD (FUNC,A,B,S,N,OPT,PLATE,ITYP)
    implicit none
    INTEGER N,OPT,PLATE,ITYP
    DOUBLE PRECISION A,B,S,FUNC
    EXTERNAL FUNC
    !
    !     TRAPZD: THIS SUBROUTINE IS TAKEN FROM THE TEXT NUMERICAL RECIPES
    !     IN FORTRAN. IT IS PART OF A SET OF ROUTINES FOR CALCULATING
    !     NUMERICAL QUADRATURE OF A DISCRETE FUNCTION.
    !
    !     THIS ROUTINE COMPUTES THE NTH STAGE REFINEMENT OF AN EXTENDED
    !     TRAPEZOIDAL RULE. FUNC IS INPUT AS THE NAME OF THE FUNCTION TO
    !     BE INTEGRATED BETWEEN THE LIMITS OF A AND B, ALSO INPUT. WHEN
    !     CALLED WITH N=1, THE ROUTINE RETURNS THE CRUDEST ESTIMATE OF
    !     THE INTEGRAL. SUBSEQUENT CALLS WITH N=2,3... (IN THAT
    !     SEQUENTIAL ORDER) WILL IMPROVE THE ACCURACY OF S BY ADDING
    !     2**(N-2) ADDITIONAL INTERIOR POINTS. S SHOULD NOT BE MODIFIED
    !     BETWEEN SEQUENTIAL CALLS.
    !
    INTEGER IT,J,ind,dummy
    DOUBLE PRECISION DEL,SUM,TNM,X
    !
    !     Set ind = -1 so that the functions that return the
    !     value of the source function will not return the
    !     integrated value if S > Src length is passed. This is
    !     necessary for the 2D integration through CIS2.
    !
    !     This function is also performed by setting slim = -1.0.
    !

    ind = 0
    dummy = 0
    IF (N.EQ.1) THEN
       IF (ITYP.EQ.1) THEN
          S = 0.5*(B-A)*(FUNC(A,OPT,PLATE,-1.0d0,ind)+ FUNC(B,OPT,PLATE,-1.0d0,ind))
       ELSEIF (ITYP.EQ.2) THEN
          S = 0.5*(B-A)*((B-A)* FUNC(A,OPT,PLATE,-1.0d0,ind) )
       ENDIF
       !
       !         WRITE(6,*) 'TRAPZD:',ITYP,A,B,&
       !                     FUNC(A,OPT,PLATE,-1.0d0,ind),&
       !                     FUNC(B,OPT,PLATE,-1.0d0,ind)
       !
    ELSE
       IT = 2**(N-2)
       TNM = IT
       DEL = (B-A)/TNM
       X = A + 0.5*DEL
       SUM = 0.0
       DO J = 1,IT
          IF (ITYP.EQ.1) THEN
             SUM = SUM + FUNC(X,OPT,PLATE,-1.0d0,ind)
          ELSEIF (ITYP.EQ.2) THEN
             SUM = SUM + (B-X) * FUNC(X,OPT,PLATE,-1.0d0,ind)
          ENDIF
          X = X + DEL
       end do
       S = 0.5*(S+(B-A)*SUM/TNM)
       !
       !         WRITE (6,*) 'TRAPZ:',S,SUM
       !
    ENDIF
    RETURN
  END SUBROUTINE TRAPZD



  DOUBLE PRECISION FUNCTION SRCION(S,SOPT,PLATEOPT,SLIM,IND)
    IMPLICIT NONE
    DOUBLE PRECISION S,SLIM
    INTEGER SOPT,PLATEOPT,IND
    !
    !     THIS FUNCTION RETURNS THE VALUE OF THE IONIZATION SOURCE
    !     STRENGTH AT THE POSITION S FROM ONE OF TWO PLATES.
    !     THIS IS NEEDED TO ENHANCE THE ACCURACY OF THE INTEGRATION
    !     ROUTINES - WHICH ARE SOMEWHAT INACCURATE USING THE CURRENT
    !     METHOD - HOWEVER - THIS ALTERNATE COULD BE COMPUTATIONALLY
    !     EXPENSIVE AND SO THE ALTERNATE CODE WILL BE RETAINED.
    !
    !     DAVID ELDER    MAY 1, 1992
    !
    INTEGER IPOS,IN
    EXTERNAL IPOS
    REAL*8 S0,S0A,S0B
    REAL STMP
    !
    !     The source function needs to maintain some data
    !     so that the integrations can go more efficiently
    !     Unfortunately ... since the integration routines
    !     are generic, this data needs to be maintained in the
    !     source functions themselves. This information
    !     include the last integrated value. The inital
    !     values for each of these would be set in the
    !     setupval subroutine.
    !
    !     Functions:
    !        ind = 1 : check for whole source integration
    !        ind = 2 : check for S > last S and return partial
    !                  integral
    !        ind = 3 : update integrated data
    !
    !     This could be done utilizing the ENTRY statement to
    !     allow multiple entry points to the same subroutine.
    !     OR ... one could write separete routines for each
    !     source function to do each item ... unfortunately this
    !     would make using a generic intergration routine quite
    !     difficult.
    !

    IF (IND.EQ.1.AND.SLIM.Gt.LENSRC) THEN
       ind = -1
       if (plateopt.eq.0) then
          SRCION = IONINT
       elseif (plateopt.eq.1) then
          srcion = ioninti
       endif
       return
    endif

    if (ind.eq.2.and.slim.gt.sionl.and.plateopt.eq.pionl) then
       ind = -2
       srcion = iionl
       slim = sionl
       return
    endif

    if (ind.eq.3) then
       iionl = slim
       pionl = plateopt
       sionl = s
       srcion = iionl
       return
    endif

    IF (S.LE.LENSRC) THEN
       IF (SOPT.EQ.0.OR.SOPT.EQ.1) THEN
          IF (PLATEOPT.EQ.0 ) THEN
             S0 = S0OUT
          ELSE
             S0 = S0IN
          ENDIF
          IF (SOPT.EQ.0) THEN
             SRCION = S0
          ELSEIF(SOPT.EQ.1) THEN
             SRCION = S0 * EXP(-S/LAMSRC)
          ENDIF
       ELSEIF (SOPT.EQ.4.OR.SOPT.EQ.5) THEN
          IF (PLATEOPT.EQ.0 ) THEN
             S0A = S0AOUT
             S0B = S0BOUT
          ELSE
             S0A = S0AIN
             S0B = S0BIN
          ENDIF
          IF (SOPT.EQ.4) THEN
             IF (S.LE.LAMSRC) THEN
                SRCION = S0A
             ELSE
                SRCION = S0B
             ENDIF
          ELSEIF(SOPT.EQ.5) THEN
             SRCION = S0A * EXP(-S/LAMSRC) + S0B
          ENDIF
       ENDIF
    ELSE
       SRCION = 0.0D0
    ENDIF

    RETURN
  END FUNCTION SRCION


  DOUBLE PRECISION FUNCTION SRCRAD(S,POPT,PLATEOPT,SLIM,IND)
    IMPLICIT NONE
    DOUBLE PRECISION S,SLIM
    INTEGER POPT,PLATEOPT,IND
    !
    !     THIS FUNCTION RETURNS THE VALUE OF THE RADIATION SOURCE
    !     STRENGTH AT THE POSITION S FROM ONE OF TWO PLATES.
    !     THIS IS NEEDED TO ENHANCE THE ACCURACY OF THE INTEGRATION
    !     ROUTINES - WHICH ARE SOMEWHAT INACCURATE USING THE CURRENT
    !     METHOD - HOWEVER - THIS ALTERNATE COULD BE COMPUTATIONALLY
    !     EXPENSIVE AND SO THE ALTERNATE CODE WILL BE RETAINED.
    !
    !     DAVID ELDER    MAY 1, 1992
    !

    REAL*8 P0

    !     The source function needs to maintain some data
    !     so that the integrations can go more efficiently
    !     Unfortunately ... since the integration routines
    !     are generic, this data needs to be maintained in the
    !     source functions themselves. This information
    !     include the last integrated value. The inital
    !     values for each of these would be set in the
    !     setupval subroutine.
    !
    !     Functions:
    !        ind = 1 : check for whole source integration
    !        ind = 2 : check for S > last S and return partial
    !                  integral
    !        ind = 3 : update integrated data
    !
    !     This could be done utilizing the ENTRY statement to
    !     allow multiple entry points to the same subroutine.
    !     OR ... one could write separete routines for each
    !     source function to do each item ... unfortunately this
    !     would make using a generic intergration routine quite
    !     difficult.
    !

    if (ind.eq.1.and.slim.gt.plensrc) then
       ind = -1
       if (plateopt.eq.0) then
          srcrad = radint
       elseif (plateopt.eq.1) then
          srcrad = radinti
       endif
       return
    endif

    if (ind.eq.2.and.slim.gt.sradl.and.plateopt.eq.pradl) then
       ind = -2
       srcrad = iradl
       slim = sradl
       return
    endif

    if (ind.eq.3) then
       iradl = slim
       srcrad = iradl
       pradl = plateopt
       sradl = s
       return
    endif

    IF (S.LE.PLENSRC) THEN
       IF (PLATEOPT.EQ.0 ) THEN
          P0 = P0OUT
       ELSE
          P0 = P0IN
       ENDIF
       IF (POPT.EQ.0.OR.POPT.EQ.2) THEN
          SRCRAD = P0
       ELSEIF (POPT.EQ.1.OR.POPT.EQ.3) THEN
          SRCRAD = P0 * EXP(-S/PLAMSRC)
       ENDIF
    ELSE
       SRCRAD = 0.0D0
    ENDIF

    RETURN
  END FUNCTION SRCRAD




  SUBROUTINE SOLEDGE(ix1,ix2,qtim)
    use mod_params
    use mod_comt2
    use mod_comtor
    use mod_comxyt
    use yreflection
    IMPLICIT  NONE
    integer   ix1,ix2
    real :: qtim
    !      integer :: cioptf_soledge
    !
    !
    ! *****************************************************************
    ! *                                                               *
    ! * SOLEDGE:  THIS ROUTINE SETS THE TEMPERATURE AND DENSITY OF    *
    ! *           SOL PLASMA- IT ALSO ASSIGNS THE VELOCITY AND        *
    ! *           ELECTRIC FIELD IN THE SOL. THIS IS DONE USING A SET *
    ! *           OF SIMULTANEOUS EQUATIONS WHICH FORM A SOMEWHAT     *
    ! *           CONSISTENT MODEL OF THE PLASMA EDGE.                *
    ! *                                                               *
    ! * DAVID ELDER     OCT  1991                                     *
    ! *                                                               *
    ! *****************************************************************
    !
    ! strip down soledge to code for sol12 and sol13 only to start

    DOUBLE PRECISION ARG1,ARG2,ARG3,ARG4
    double precision n,v
    INTEGER IR,J,IK,NS,PLATEOPT,IRLIMIT
    INTEGER IKMID,ikstart,ikend,ikfirst,iklast
    DOUBLE PRECISION DS ,V0,V0I,PINF,PINFI
    REAL*8 ACT_PRESS
    DOUBLE PRECISION RCF,RCFI
    !DOUBLE PRECISION CIS1,CIS2
    !EXTERNAL CIS1,CIS2
    DOUBLE PRECISION DS1,DS2,DP1,DP2,DT1,DT2,NB1,NB2,DTI1,DTI2
    DOUBLE PRECISION S,SMAX,NBP,NBPI,TEBP,TEBPI,TIBP,TIBPI
    double precision spredi,spredo,sprev,predls,sinj
    DOUBLE PRECISION LPPA,LPPAI,LPPAEB,LPPAEI,LPPAIB,LPPAII
    DOUBLE PRECISION MASSI,TMP
    DOUBLE PRECISION DELTAS,pmax
    DOUBLE PRECISION FSRC,LMSRC,LNSRC
    DOUBLE PRECISION LSSIZ,LPSIZ,MFACT,MFACT2

    double precision sl1,sl2,slv,v1,v1i,t1,t1i,t2,t2i,n1,n1i,lprad,lpradi,news,ti1,ti1i
    double precision sl1i,sl2i,slvi
    double precision sl1a,sl1ai,sl1b,sl1bi  

    real*8 :: y_1b,te_1b,ti_1b,cs_1b,e_1b,y_1t,te_1t,ti_1t,cs_1t,e_1t,cl_1
    real*8 :: y_2b,te_2b,ti_2b,cs_2b,e_2b,y_2t,te_2t,ti_2t,cs_2t,e_2t,cl_2
    real*8 :: mult
    real*8 :: powse,powsi
    
    real*8 :: tgscal,dstep
    
    real*8 :: soli,solprn,e_scale,v_scale
    real*8 :: dy,dt
    real :: rizb
    integer :: ix,iy

    integer,external :: ipos
    
    !DOUBLE PRECISION SRCION,SRCRAD
    !EXTERNAL SRCION,SRCRAD

    CHARACTER*80 INFO
    !     REAL SOLVTEMP
    !     EXTERNAL SOLVTEMP
    !
    !     Variables for SOL 16+
    !
    double precision rootn,gamman,sact,strt,stmp
    integer ikn,in
    !     double precision soltelast,soltecur,soltilast,solticur
    !     double precision solnelast,solnecur,solvellast,solvelcur
    !     double precision solcorlast,solcorcur
    !     double precision solprn,solprhn,solpcxn,solphn,solpein
    !     double precision helpi,soli
    !
    !
    !       SET IK values for inner loops        
    !

    rizb = real(cizb)

    ikstart = 1
    ikmid = maxn/2
    ikend = maxn

    do ix = ix1,ix2
	   !write(0,*) 'ix,x = ', ix,xs(ix)
       !
       !     jdemod - do not read the starting target conditions from the grid
       !            - they should be loaded from the specified target conditions
       !              for each ring
       !
       ! should use qtembs(iqx)?
       !
       
       tebp = ctembs(ix,0)
       tibp = ctembsi(ix,0)
       nbp  = crnbs(ix,0) 
       V0  = - SQRT(0.5*EMI*(TEBP+TIBP)*(1+RIZB)/CRMB)

       tebpi = ctembs(ix,0)
       tibpi = ctembsi(ix,0)
       nbpi  = crnbs(ix,0) 
       V0i  = - SQRT(0.5*EMI*(TEBPi+TIBPi)*(1+RIZB)/CRMB)

       ! Comment needed to explain what is being calculated here.
       ! Target power flux ...
       IF (cioptf_soledge.eq.11.or.CIOPTF_SOLEDGE.EQ.12.OR.CIOPTF_SOLEDGE.EQ.14.or.cioptf_soledge.eq.16.or.cioptf_soledge.eq.18.or.cioptf_soledge.eq.19) THEN
          LPPA = (2.0*TIBP+5.0*TEBP)*1.602192E-19*NBP* DABS(V0)
          LPPAI = (2.0*TIBPI+5.0*TEBPI)*1.602192E-19*NBPI* DABS(V0I)
       ENDIF

       IF (CIOPTF_SOLEDGE.EQ.13.OR.CIOPTF_SOLEDGE.EQ.15.or.cioptf_soledge.eq.17.or.cioptf_soledge.eq.20) THEN
          LPPAEB = 5.0*TEBP*1.602192E-19*NBP*DABS(V0)
          LPPAEI = 5.0*TEBPI*1.602192E-19*NBPI*DABS(V0I)

          LPPAIB = 2.0*TIBP*1.602192E-19*NBP*DABS(V0)
          LPPAII = 2.0*TIBPI*1.602192E-19*NBPI*DABS(V0I)

          LPPA = LPPAEB
          LPPAI= LPPAEI
       ENDIF


       !
       !       Set up pressure value
       !
       ! Te=Ti
       IF (cioptf_soledge.eq.11.or.CIOPTF_SOLEDGE.EQ.12.OR.CIOPTF_SOLEDGE.EQ.14.OR.CIOPTF_SOLEDGE.EQ.16.or.cioptf_soledge.eq.18.or.cioptf_soledge.eq.19) THEN
          PINF  = 4.0*NBP*ECH*TEBP
          PINFI = 4.0*NBPI*ECH*TEBPI
       ELSEIF (CIOPTF_SOLEDGE.EQ.13.OR.CIOPTF_SOLEDGE.EQ.15.or.cioptf_soledge.eq.17.or.cioptf_soledge.eq.20) THEN
          PINF  = NBP*ECH*(2.0*TEBP+2.0*TIBP)
          PINFI = NBPI*ECH*(2.0*TEBPI+2.0*TIBPI)
       ENDIF
       !

       MASSI = CRMB * AMU
       !
       !     SET UP THE IONIZATION SOURCE ARRAY FIRST. THIS HAS MORE
       !     POINTS THAN THE BACKGROUND N,V,T,E ARRAYS IN ORDER TO MAKE
       !     THE NUMERICAL INTEGRATION FOR THE SOURCE TERMS MORE ACCURATE.
       !
       !       Execute the ionization source code only for
       !       options that require it. i.e. NOT Sol 21
       !

       !
       ! Setting parameters 
       !


       FSRC = CFSRC

       IF (CSOPT.EQ.0) THEN
          LNSRC = CSOLLS
          IF (CSOLLS.GT.0.5) LNSRC = 0.5
          LMSRC = 0.0
       ELSEIF (CSOPT.EQ.1) THEN
          LNSRC = 0.5
          LMSRC = CSOLLS
       ELSEIF (CSOPT.EQ.4.OR.CSOPT.EQ.5) THEN
          LNSRC = CSOLLT
          LMSRC = CSOLLS
       ENDIF



       !       Write out parameters for testing ...
       !
       !write(6,*) nbp,tebp,v0,nbpi,tebpi,v0i
       !
       !       WILL REPLACE WITH SPECIFIED VALUES IF FOUND.

       !        IF (FLUXROPT.EQ.1)  CALL FLUXRLOOK(IR,FSRC,LMSRC,LNSRC)

       !         WRITE(6,*) 'FLUX:',FSRC,LNSRC,LMSRC

       !       CALL THE SETUP ROUTINE : IT ESTABLISHES ALL THE VALUES
       !       REQUIRED FOR THE VARIOUS IONIZATION AND RADIATION OPTIONS
       !       THAT UNDERLIE THE CALCULATION OF THE SOL CHARACTERISTICS.
       !

       CALL SETUPVAL(CSOPT,CPOPT,FSRC,LNSRC,LMSRC,SMAX,&
            NBP,V0,TEBP,NBPI,V0I,TEBPI,RCF,RCFI,LSSIZ,LPSIZ,IX,&
            LPPA,LPPAI)

       !
       !       OUTER PLATES ....
       !

       !        if (ikopt.eq.1.or.ikopt.eq.3) then 

       sprev = 0.0
       !
       !do ik = ikstart, ikmid-1
       !  write(0,*) 'ix,ik,sd = ', ix, ik, sd(ik)
       !enddo
       ! sazmod
       ! If this ix is a ring that is in the step region, so has a shorter
       ! smax, we need to account for how the step will affect the solution
       ! via shortening smax. Need to identify a new maxn, which will be
       ! the n in the original s array where s > yabsorb1a_step (or 2a_step).
       ! Need to also do in the second loop a few hundred lines down.
       ! BUG: This only checks if you are in the step2 region (right side).
       !      Obviously would like to be able to have this work for the
       !      left side too. Would likely require indicating regions of
       !      the plasma.
       if (vary_absorb.eq.1) then
       
         ! Check if in step region.
         ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         if (ix.le.ix_step2) then
         
           ! If we are then we need to adjust the s values accordingly.
           ! Going in positive s direction.
           !do ik = ikstart, ikmid-1 
           do ik = ikstart, ikend 
             s = sd(ik)
             if (s.ge.yabsorb2a_step) then
             
               ! Set our new smax and leave.
               !write(0,*) 'L53: Old smax: ', smax
               !smax = s
               !write(0,*) 'L53: ix,x,smax_old,smax_new=',ix,xs(ix),smax,smax/2.0+s
               smax = smax/2.0 + s
               !write(0,*) 'L53: New smax: ', smax
               exit
               
             endif
           end do
         endif
       endif
       !write(0,*) 'L53: ix,x,smax=', ix, xs(ix), smax
       DO IK = ikstart, IKMID-1

          S  = sd(ik)
		  
          !     Calculate revised pressure!
          act_press = pinf
          if (sol13_pdist.gt.0.0) then 
             act_press = pinf * min((s/(sol13_pdist * smax)), dble(1.0)) * sol13_padd + pinf     
          else
             act_press = pinf * (1.0+sol13_padd)
          endif

          ! sazmod - My attempt at creating a constant Te background
		  !          to simulate a sheath-limited regime (opt 11). Just
		  !          going to copy opt. 12 for the most part.
		  

          if (cioptf_soledge.eq.11) then
		     te(ik) = tebp  ! This is really the only different line.
		     ti(ik) = te(ik)
		     
		     !write(0,*) 's/smax=',s,smax
		     
		     SOLI = CIS1(S,SRCION,CSOPT,0)
             GAMMAN = NBP*V0 + SOLI + RCF * S
             !
             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0) GAMMAN = 0.0
             !
             call calcnv(te(ik),ti(ik),gamman,act_press,n,v)
             !
             ne(ik) = n
             vb(ik) = v
             ga(ik) = gamman
             
             
             ! print s coordinate
			 !if (mod(ik, ikmid/8).eq.0) then
			 !  write(0,*) 'ik,s,vb= ',ik,s,v
			 !endif

          !IF (CIOPTF_SOLEDGE.EQ.12) THEN
          elseif(cioptf_soledge.eq.12) then
             !
             !           CALCULATE BACKGROUND TEMPERATURE
             !
             if (cpowopt.eq.0) then
                powse = lppa*s
             elseif (cpowopt.eq.1) then
                powse = lppa*s - 0.5 * lppa * s**2/(smax/2.0)
             endif

             te(ik) = (TEBP**3.5 + 7.0/(2.0*CK0)* (powse+CIS2(S,SRCRAD,CPOPT,0)))**(2.0/7.0)
             ti(ik) = te(ik)
             !
             !           CALCULATE DENSITY
             !
             SOLI = CIS1(S,SRCION,CSOPT,0)
             GAMMAN = NBP*V0 + SOLI + RCF * S
             !
             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0) GAMMAN = 0.0
             !
             call calcnv(te(ik),ti(ik),gamman,act_press,n,v)
             !
             ne(ik) = n
             vb(ik) = v
             ga(ik) = gamman
             !
             !            WRITE(6,*) 'NUMBERS:',IK,IR,KTEBS(IK,IR),KTIBS(IK,IR), KNBS(IK,IR),KVHS(IK,IR)
             !            WRITE(6,*) 'OTHERS1:',ROOTN,GAMMAN,S,act_press
             !            WRITE(6,*) 'OTHERS2:',(PINF/(2*ECH*KTEBS(IK,IR)))**2,
             !     >            -2.0*(MASSI*GAMMAN**2)/(ECH*KTEBS(IK,IR)),
             !     >             NBP*V0,SOLI,RCF
             !

             !write(0,'(a,i8,12(1x,g18.8))') 'F:',ik,smax,s,sd(ik),te(ik),ti(ik),ne(ik),vb(ik)

          ELSEIF (CIOPTF_SOLEDGE.EQ.13.OR.CIOPTF_SOLEDGE.EQ.15) THEN
             !
             !
             !           SET DIFFERENCE BETWEEN 13 AND 15
             !
             IF (CIOPTF_SOLEDGE.EQ.13) THEN
                MFACT = 7.0/2.0
                MFACT2 = 1.0
             ELSEIF (CIOPTF_SOLEDGE.EQ.15) THEN
                MFACT = 7.0/4.0
                MFACT2 = S / (SMAX/2.0)
             ENDIF
             !
             !           CALCULATE ELECTRON TEMPERATURE
             !
             TMP = CIS2(S,SRCRAD,CPOPT,0)
             

             if (cpowopt.eq.0) then
                powse = lppaeb*s
             elseif (cpowopt.eq.1) then
                powse = lppaeb*s - 0.5 * lppaeb * s**2/(smax/2.0)
             endif


             te(ik) = (TEBP**3.5 + (MFACT/CK0)*(powse*MFACT2+TMP))**(2.0/7.0)
             !
             !             WRITE (6,*) 'NUM1:',TMP,LPPAEB,TEBP,MFACT,
             !     >                CK0,MFACT2
             !
             !           CALCULATE ION TEMPERATURE
             !

             if (cpowopt.eq.0) then
                powsi = lppaib*s
             elseif (cpowopt.eq.1) then
                powsi = lppaib*s - 0.5 * lppaib * s**2/(smax/2.0)
             endif

             ti(ik) = (TIBP**3.5 + (MFACT/CK0I)*(powsi*MFACT2))**(2.0/7.0)
             !
             !           CALCULATE DENSITY
             !
             SOLI = CIS1(S,SRCION,CSOPT,0)
             GAMMAN = NBP*V0 + SOLI + RCF * S

             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0) GAMMAN = 0.0

             call calcnv(te(ik),ti(ik), gamman,act_press,n,v)

             ne(ik) = n
             vb(ik) = v
             ga(ik) = gamman
             !
             !            WRITE(6,*) 'NUMBERS:',IK,IR,KTEBS(IK,IR),KTIBS(IK,IR),
             !     >             KNBS(IK,IR),KVHS(IK,IR),RCF
             !            WRITE(6,*) 'OTHERS1:',ROOTN,GAMMAN,S
             !            WRITE(6,*) 'OTHERS2:',(PINF/(ECH*(KTEBS(IK,IR)
             !     >            +KTIBS(IK,IR)) ) ) **2,
             !     >            -4.0*(MASSI*GAMMAN**2)/(ECH*(KTEBS(IK,IR)
             !     >            +KTIBS(IK,IR))),
             !     >             NBP*V0,SOLI,RCF*S
             !
             !write(0,'(a,i8,12(1x,g18.8))') 'F:',ik,smax,s,sd(ik),te(ik),ti(ik),ne(ik),vb(ik)



          ELSEIF (CIOPTF_SOLEDGE.EQ.14) THEN

             SOLI = CIS1(S,SRCION,CSOPT,0)
             GAMMAN = NBP*V0 + SOLI + RCF * S

             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0) GAMMAN = 0.0
             !
             !           CALCULATE THE BACKGROUND TEMPERATURE
             !
             !            write(6,*) 'soli:',ir,ik,soli,gamman

             IF (IK.EQ.1.AND.S.LE.0.0) THEN
                te(ik) = TEBP
             ELSEIF (IK.EQ.1.AND.S.GT.0.0) THEN
                ARG1 = sd(ik) * 5.0 * GAMMAN * ECH
                solprn  = CIS1(S,SRCRAD,CPOPT,0)
                ARG2 = sd(ik) * (LPPA+solprn)
                ARG3 = TEBP * CK0
                ARG4 = -CK0
                !              write(6,*) 'args:',arg1,arg2,arg3,arg4,solprn,lppa
                te(ik) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
             ELSE
                DELTAS = sd(ik) - sd(ik-1)
                ARG1 = DELTAS * 5.0 * GAMMAN * ECH
                solprn = CIS1(S,SRCRAD,CPOPT,0)
                ARG2 = DELTAS * ( LPPA + solprn)
                ARG3 = te(ik-1) * CK0
                ARG4 = -CK0
                !              write(6,*) 'args:',arg1,arg2,arg3,arg4,solprn,lppa
                te(ik) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
             ENDIF

             ti(ik) = te(ik)
             !
             !           CALCULATE DENSITY
             !
             call calcnv(te(ik),ti(ik),gamman,act_press,n,v)
             !
             ne(ik) = n
             vb(ik) = v
             ga(ik) = gamman

          endif

       end do



       sprev = 0.0

	   ! Same as previous loop, just for this region. 
	   ! Maybe this loop isn't needed...
	   if (vary_absorb.eq.1) then
       
         ! Check if in step region.
         ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         if (ix.ge.ix_step2) then
         
           ! If we are then we need to adjust the s values accordingly.
           ! Going in positive s direction.
           do ik = ikend, ikmid, -1  
             s = sd(ik)
             if (s.gt.yabsorb2a_step) then
             
               ! Set our new smax and leave.
               !write(0,*) 'L53: Old smax:', smax
               !smax = s
               !write(0,*) 'L53: New smax:', smax
               exit
               
             endif
           end do
         endif
       endif

       DO IK = ikend, IKMID ,-1

          S  = SMAX - sd(ik)

          !write(0,*) 'SMAX:',smax,s,sd(ik)
          
          !
          !     Calculate revised PINF
          !
                    if (sol13_pdist.gt.0.0) then 
                       act_press = pinf *  min((s/(sol13_pdist*smax)), dble(1.0)) * sol13_padd +   pinf     
                    else
                       act_press = pinf * (1.0+sol13_padd)
                    endif
          !     
          !
          
          ! sazmod - Same as the above section, constant Te opt 11.
          if (cioptf_soledge.eq.11) then
             te(ik) = tebp
             ti(ik) = te(ik)
             
             SOLI = CIS1(S,SRCION,CSOPT,1)
             GAMMAN = NBPI*V0I + SOLI + RCFI * S
             !
             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)   GAMMAN = 0.0


             call calcnv(te(ik),ti(ik),gamman,act_press,n,v)

             ne(ik) = n
             vb(ik) = -v
             ga(ik) = gamman
          
          !IF (CIOPTF_SOLEDGE.EQ.12) THEN
          elseif(cioptf_soledge.eq.12) then
             !
             !           CALCULATE BACKGROUND TEMPERATURE
             !


             if (cpowopt.eq.0) then
                powse = lppai*s
             elseif (cpowopt.eq.1) then
                powse = lppai*s - 0.5 * lppai * s**2/(smax/2.0)
             endif


             te(ik) = (TEBPI**3.5 + 7.0/(2.0*CK0)*(powse+CIS2(S,SRCRAD,CPOPT,1)))**(2.0/7.0)

             ti(ik) = te(ik)
             !
             !           CALCULATE DENSITY
             !
             SOLI = CIS1(S,SRCION,CSOPT,1)
             GAMMAN = NBPI*V0I + SOLI + RCFI * S
             !
             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)   GAMMAN = 0.0


             call calcnv(te(ik),ti(ik),gamman,act_press,n,v)

             ne(ik) = n
             vb(ik) = -v
             ga(ik) = gamman

             !write(0,'(a,i8,12(1x,g18.8))') 'S:',ik,smax,s,sd(ik),te(ik),ti(ik),ne(ik),vb(ik)


          ELSEIF (CIOPTF_SOLEDGE.EQ.13.OR.CIOPTF_SOLEDGE.EQ.15) THEN



             !
             !           SET DIFFERENCE BETWEEN 13 AND 15
             !
             IF (CIOPTF_SOLEDGE.EQ.13) THEN
                MFACT = 7.0/2.0
                MFACT2 = 1.0
             ELSEIF (CIOPTF_SOLEDGE.EQ.15) THEN
                MFACT = 7.0/4.0
                MFACT2 = S / (SMAX/2.0)
             ENDIF
             !
             !           CALCULATE ELECTRON TEMPERATURE
             !
             if (cpowopt.eq.0) then
                powse = lppaei*s
             elseif (cpowopt.eq.1) then
                powse = lppaei*s - 0.5 * lppaei * s**2/(smax/2.0)
             endif
             
             TMP = CIS2(S,SRCRAD,CPOPT,1)
             !write(0,'(a,i10,10(1x,g18.8))') 'Ti:',ik,tebpi, lppaei,ck0,s,mfact,mfact2,tmp,CIS2(S,SRCRAD,CPOPT,1),TEBPI**3.5 + (MFACT/CK0)*(LPPAEI*S*MFACT2 +CIS2(S,SRCRAD,CPOPT,1))
             te(ik) = (TEBPI**3.5 + (MFACT/CK0)*(powse*MFACT2 +TMP))**(2.0/7.0)
             !
             !           CALCULATE ION TEMPERATURE
             !

             if (cpowopt.eq.0) then
                powsi = lppaii*s
             elseif (cpowopt.eq.1) then
                powsi = lppaii*s - 0.5 * lppaii * s**2/(smax/2.0)
             endif

             ti(ik) = (TIBPI**3.5 + (MFACT/CK0I)* (powsi*MFACT2))**(2.0/7.0)

             !            CALCULATE DENSITY

             SOLI = CIS1(S,SRCION,CSOPT,1)
             GAMMAN = NBPI*V0I + SOLI + RCFI * S

             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)    GAMMAN = 0.0

             call calcnv(te(ik),ti(ik),gamman,act_press,n,v)

             ne(ik) = n
             vb(ik) = -v
             ga(ik) = gamman

             !write(0,'(a,i8,12(1x,g18.8))') 'S:',ik,smax,s,sd(ik),te(ik),ti(ik),ne(ik),vb(ik)

          ELSEIF (CIOPTF_SOLEDGE.EQ.14) THEN

             SOLI = CIS1(S,SRCION,CSOPT,1)
             GAMMAN = NBPI*V0I + SOLI + RCFI * S

             IF (GAMMAN.GT.-LO.AND.FLUXROPT.EQ.0)      GAMMAN = 0.0
             !
             !           CALCULATE THE BACKGROUND TEMPERATURE
             !
             IF (IK.EQ.maxn.AND.S.LE.0.0) THEN
                te(ik) = TEBPI
             ELSEIF (IK.eq.maxn.AND.S.GT.0.0) THEN
                ARG1 = (SMAX-sd(ik)) * 5.0 * GAMMAN * ECH
                solprn = CIS1(S,SRCRAD,CPOPT,1)
                ARG2 = (SMAX-sd(ik)) * ( LPPAI + solprn)
                ARG3 = TEBPI * CK0
                ARG4 = -CK0
                te(ik) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
             ELSE
                DELTAS = sd(ik+1) - sd(ik)
                ARG1 = DELTAS * 5.0 * GAMMAN * ECH
                solprn = CIS1(S,SRCRAD,CPOPT,1)
                ARG2 = DELTAS * ( LPPAI + solprn)
                ARG3 = te(ik+1) * CK0
                ARG4 = -CK0
                te(ik) = SOLVTEMP(ARG1,ARG2,ARG3,ARG4)
             ENDIF


             ti(ik) = te(ik)
             !
             !           CALCULATE DENSITY
             !
             call calcnv(te(ik),ti(ik),gamman,act_press,n,v)

             ne(ik) = n
             vb(ik) = -v
             ga(ik) = gamman


          endif

       end do


       !       CALCULATE ELECTRIC FIELD and temperature gradients
       !
       !
       !       IN THE FOLLOWING EQUATIONS THE FACTOR E CANCELS WITH THE
       !       SAME FACTOR USED IN CONVERTING T IN EV TO KT.
       !
       DO IK = ikstart,ikend
          IF (IK.EQ.1) THEN
             DS1 =sd(ik+1) - sd(ik)
             DP1 =ne(ik+1)*te(ik+1)-ne(ik)*te(ik)
             DT1 =(te(ik+1)-te(ik))
             DTI1 =(ti(ik+1)-ti(ik))
             NB1 =0.5*(ne(ik+1)+ne(ik))
          ELSE
             DS1 =sd(ik) - sd(ik-1)
             DP1 =ne(ik)*te(ik)-ne(ik-1)*te(ik-1)
             DT1 =(te(ik)-te(ik-1))
             DTI1 =(ti(ik)-ti(ik-1))
             NB1 =0.5*(ne(ik)+ne(ik-1))
          ENDIF

          if (ik.eq.maxn) then
             DS2 =sd(ik) - sd(ik-1)
             DP2 =ne(ik)*te(ik)-ne(ik-1)*te(ik-1)
             DT2 =(te(ik)-te(ik-1))
             DTI2 =(ti(ik)-ti(ik-1))
             NB2 =0.5*(ne(ik)+ne(ik-1))
          ELSE
             DS2 =sd(ik+1) - sd(ik)
             DP2 =ne(ik+1)*te(ik+1)-ne(ik)*te(ik)
             DT2 =(te(ik+1)-te(ik))
             DTI2 =(ti(ik+1)-ti(ik))
             NB2 =0.5*(ne(ik+1)+ne(ik))
          ENDIF


          if (nb1.eq.0.0.or.nb2.eq.0.0.or.ds1.eq.0.0.or.ds2.eq.0.0)  then 
             ef(ik) = 0.0            
             teg(ik) = 0.0
             tig(ik) = 0.0
          elseif (ik.eq.1) then 
             ef(1) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
             teg(1) = dt1/ds1
             tig(1) = dti1/ds1

          elseif (ik.eq.maxn) then

             ef(maxn) = -(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2
             teg(maxn) = dt2/ds2
             tig(maxn) = dti2/ds1

          else
             ef(ik) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1) + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
             teg(ik) = 0.5*(dt1/ds1 + dt2/ds2)
             tig(ik) = 0.5*(dti1/ds1 + dti2/ds2)

          endif


       end do

       !write(6,'(a,i8)') 'Plasma on surface:',ix
       do ik=1,maxn
          write(6,'(i8,12(1x,g18.8))') ik,sd(ik),yd(ik),te(ik),ti(ik),ne(ik),vb(ik),ga(ik),ef(ik),teg(ik),tig(ik)
       end do


       !
       ! Pull out the boundary conditions for applying Vb and ef on field lines that connect
       ! to the probe. Using the simple model for E and Vb from SOL option 4. (linear decay
       ! to stagnation between surfaces
       !

       ! interpolate at youts(iy) to obtain
       ! fix sd axis

       y_1b = yd(1)
       te_1b = te(1)
       ti_1b = ti(1)
       cs_1b =9.79E+03 * SQRT(((te_1b+ti_1b)/2)* (1.0+REAL(CIZB))/CRMB)  

       y_1t = -qedges(iqxs(ix),1)
       in = ipos(y_1t,yd,maxn)
       te_1t = te(in)
       ti_1t = ti(in)
       cs_1t =9.79E+03 * SQRT(((te_1t+ti_1t)/2)* (1.0+REAL(CIZB))/CRMB)  

	   !write(0,*) 'y_1b,y_1t,yabsorb1a_step',y_1b,y_1t,yabsorb1a_step

	   cl_1 = (y_1t+y_1b)/2.0

       ! sazmod - The conn length will change if in step region.
       if (vary_absorb.eq.1) then
       
         ! Check if in step region.
         ix_step1 = ipos(xabsorb1a_step, xs, nxs-1)
         
         ! Adjust connection length.
         if (ix.le.ix_step1) then
           cl_1 = (y_1t + yabsorb1a_step) / 2.0  ! Y < 0
           !write(0,*) 'cl_1 adjusted to ',cl_1
         endif
       endif
       !write(0,*) 'cl_1 = ',cl_1
       
       ! max efield
       e_1b = te_1b/cl_1
       e_1t = te_1t/cl_1


       y_2b = qedges(iqxs(ix),2)
       in = ipos(y_2b,yd,maxn)
       te_2b = te(in)
       ti_2b = ti(in)
       cs_2b =9.79E+03 * SQRT(((te_2b+ti_2b)/2)* (1.0+REAL(CIZB))/CRMB)  

       y_2t = yd(maxn)
       te_2t = te(maxn)
       ti_2t = ti(maxn)
       cs_2t =9.79E+03 * SQRT(((te_2t+ti_2t)/2)* (1.0+REAL(CIZB))/CRMB) 
       
       cl_2 = (y_2t+y_2b)/2.0 
       
       !write(0,*)'y_2b,y_2t,yabsorb2a_step',y_2b,y_2t,yabsorb2a_step
       ! sazmod - The conn length will vary if in step region.
       if (vary_absorb.eq.1) then
       
         ! Check if in step region.
         ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         
         ! Adjust connection length.
         if (ix.le.ix_step2) then           
           cl_2 = (y_2b + yabsorb2a_step) / 2.0  ! Y > 0   
           !write(0,*) 'cl_2 adjusted to ',cl_2       
         endif
       endif
       !write(0,*) 'ix, x, cl_1, cl_2 = ',ix, xs(ix), cl_1, cl_2
       
       ! max efield
       e_2b = te_2b/cl_2
       e_2t = te_2t/cl_2

       TGSCAL = (1.6E-19)/(CRMI*1.673E-27) * QTIM *QTIM 

       do iy = -nys,nys


          in = ipos(youts(iy),yd,maxn)

          ! only overwrite if youts is inside of range
          if (in.gt.1.and.in.lt.maxn) then 


             if (.not.(youts(iy).ge.yd(in-1).and.youts(iy).le.yd(in))) then 
                write(0,*) 'Warning Y not in correct bin:', yd(in),youts(iy),yd(in-1)
                stop 'Y index error'
             endif

                
             dy = youts(iy)-yd(in-1)
             dt = yd(in) - yd(in-1)
             ctembs(ix,iy) = te(in-1) + dy/dt * (te(in)-te(in-1))
             ctembsi(ix,iy)= ti(in-1) + dy/dt * (ti(in)-ti(in-1))

             DSTEP = TGSCAL *  QS(IQXS(IX)) * QS(IQXS(IX))
             ctegs(ix,iy) = (teg(in-1) + dy/dt * (teg(in)-teg(in-1))) * dstep
             ctigs(ix,iy) = (tig(in-1) + dy/dt * (tig(in)-tig(in-1))) * dstep

             crnbs(ix,iy)= ne(in-1) + dy/dt * (ne(in)-ne(in-1))

             ! These weird scaing factors are required due to the way the
             ! force coefficients were originally calculated in LIM.
             ! CFEXZS and CFVHXS both have temperature scaling dependencies since the
             ! original LIM code only calculates the velocity and efield along the field
             ! line at the separatrix and scales everything outboard proportional to the
             ! temperature at the limiter. To use the same coefficients with this option
             ! the inverse temperature scaling is applied here so it cancels out.

             if (ctembs(ix,iy).le.0.0.or.ctembsi(ix,iy).le.0.0.or.crnbs(ix,iy).le.0.0) then
                write(0,'(a,2i8,10(g12.5))') 'Less than zero:',ix,iy, ctembs(ix,iy),ctembsi(ix,iy),crnbs(ix,iy)
                stop
             endif
             
             ! jdemod - the scaling factors for LIM original arrays really make no sense for more complex backgrounds
             ! the vel_efield_opt is intended to transition to using efield and velplasma exclusively and includes revised
             ! coefficients that only include timestep scaling and not scaling to temperature at separatrix (CFVHXS, CFVEXS)
             if (vel_efield_opt.eq.0) then 
                e_scale = ctbin/ctembs(ix,iy)
                if ((ctembs(ix,iy).gt.0.0).and.(ctembsi(ix,iy).gt.0.0).and.ctbin.ge.0.0.and.ctibin.ge.0.0) then 
                   v_scale = sqrt((ctbin+ctibin)/(ctembs(ix,iy)+ctembsi(ix,iy)))
                else
                   v_scale = 0.0
                   write(0,'(a,2i8,10(1x,g12.5))') 'V_scale error:',ix,iy,ctbin,ctibin,ctembs(ix,iy),ctembsi(ix,iy)
                endif
             elseif (vel_efield_opt.eq.1) then 
                e_scale = 1.0
                v_scale = 1.0
             endif
                
             ! sazmod
             ! Apply the overall scaling of the plasma velocity. Obviously
             ! nothing will change if vel_mod = 1.0. Note: mod_v_fact
             ! can still do it's own thing, and is applied after vel_mod.
             v_scale = v_scale * vel_mod
               
               
             ! sazmod 
             ! Calculate velplasma differently so the step in the absorbing walls
             ! is factored in and has an appropriate stagnation point. 
             if (vary_absorb.eq.1) then  
                             
                ! This is just the equation for a line with the two points
                ! (cl_1, -cs_1b)
                ! (cl_2,  cs_1b)
                !mult = cs_1b * (2*(youts(iy)-cl_2)/(cl_2-cl_1) + 1)
                !write(0,*) 'youts(iy) =', iy, youts(iy)
                !mult = cs_1b * (2*youts(iy)/(cl_1+cl_2) - 1)
                mult = cs_1b * ((youts(iy)-2*cl_1)/(cl_2-cl_1) - 1)
                
                ! Multiply by the v_scale factor because 3DLIM is weird.
                velplasma(ix,iy,1) = mult * v_scale
                
                ! Printout some numbers.
                !if (mod(iy, 10).eq.0) then
                !  write(0,'(a,2i6,10(1x,g12.5))')'ix,iy,youts,cs1b,mult,vp1=',ix,iy,youts(iy),cs_1b,mult,velplasma(ix,iy,1)
                !endif
                
             else
	             velplasma(ix,iy,1)= (vb(in-1) + dy/dt * (vb(in)-vb(in-1))) * v_scale
	             efield(ix,iy,1)= (ef(in-1) + dy/dt * (ef(in)-ef(in-1))) * e_scale
	     endif
             

             ! calculate second set of velplasma and efield for field lines that end on probe surfaces. 

             ! Y < 0
			 !write(0,*) 'youts, qedges(iqxs(ix),2)', youts(iy), qedges(iqxs(ix),2)
             if (youts(iy).lt.ymin) then
                velplasma(ix,iy,2) = 0.0
                efield(ix,iy,2) = 0.0
             elseif (youts(iy).gt.ymax) then 
                velplasma(ix,iy,2) = 0.0
                efield(ix,iy,2) = 0.0
             elseif (youts(iy).lt.-qedges(iqxs(ix),1)) then 

                ! set sign and scaling ... "-" towards -y surface
                !mult =  youts(iy)/cl_1 + 1.0
                
                ! sazmod - I think this should be -youts(iy)/cl_1, otherwise
                ! there isn't a stagnation point on the left side of the probe.
                mult = -youts(iy)/cl_1 + 1.0

                ! Y < 0

				! y < 0 and y < --5
                if (youts(iy).lt.-cl_1) then
				   !write(0,*) 'Section 1'
                   ! half nearest absorbing surface
                   velplasma(ix,iy,2) = cs_1b * mult * v_scale
                   efield(ix,iy,2) = e_1b * mult * e_scale
                   !write(0,*) 'Section 1 ', youts(iy), velplasma(ix,iy,2)
                   
                ! y < 0 and y > --5 (i.e. never!). This should be fixed but
                ! it doesn't seem to cause errors... yet.
                else
                   !write(0,*) 'Section 2'
                   ! on probe side of midpoint
                   velplasma(ix,iy,2) = cs_1t * mult * v_scale
                   efield(ix,iy,2) = e_1t * mult * e_scale 
                   !write(0,*) 'Section 2 ', youts(iy), velplasma(ix,iy,2) 
                endif


             elseif (youts(iy).gt.qedges(iqxs(ix),2)) then
                ! Y > 0

                ! set sign and scaling ... "-" towards -y surface
                mult =  youts(iy)/cl_2 - 1.0
                !write(0,*) 'Y>0 mult = ', mult

                ! Y < 0

				! y < 5 and y > 0
                if (youts(iy).lt.cl_2) then
                   !write(0,*) 'Section 3'
                   ! on probe side of midpoint
                   velplasma(ix,iy,2) = cs_2b * mult * v_scale
                   efield(ix,iy,2) = e_2b * mult * e_scale
                   !write(0,*) 'Section 3 ', youts(iy), velplasma(ix,iy,2)
                   
                   ! sazmod
                   if (vary_absorb.eq.1) then
                     
                     ! Check if in right step region.
					 ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         
					 ! Mult. velplasma by factor.
					 if (ix.le.ix_step2) then
					   velplasma(ix,iy,2) = velplasma(ix,iy,2) * mod_v_fact
					   !write(0,*) 'mod_v_fact applied'
					 endif                
                   endif
                   
                ! y > 5 and y > 0
                else
                   !write(0,*) 'Section 4'
                   ! on probe side of midpoint
                   velplasma(ix,iy,2) = cs_2t * mult * v_scale
                   efield(ix,iy,2) = e_2t * mult * e_scale 
                   !write(0,*) 'Section 4 ', youts(iy), velplasma(ix,iy,2) 
                   
                   ! sazmod
                   ! My own little ad-hoc factor just to mess with flows
                   ! in the right step region.
                   if (vary_absorb.eq.1) then
                     
                     ! Check if in right step region.
					 ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         
					 ! Mult. velplasma by factor.
					 if (ix.le.ix_step2) then
					   velplasma(ix,iy,2) = velplasma(ix,iy,2) * mod_v_fact
					   !write(0,*) 'mod_v_fact applied'
					 endif                
                   endif
                   
                endif
             endif
          endif
       end do


       !write(0,'(a,i8)') 'Plasma data for field line:', ix
       !DO IY = -nys,nys
       !   write(0,'(i8,12(1x,g18.8))') iy,youts(iy),ctembs(ix,iy),ctembsi(ix,iy),crnbs(ix,iy),ctegs(ix,iy),ctigs(ix,iy),velplasma(ix,iy,1),efield(ix,iy,1),velplasma(ix,iy,2),efield(ix,iy,2)
       !end do



       

    end do ! ix

    ! deallocate storage for local variables.
    call end_soledge

    !stop 'debug'
    
    RETURN

  END SUBROUTINE SOLEDGE





  REAL FUNCTION SOLVTEMP (ARG1,ARG2,ARG3,ARG4)
    IMPLICIT NONE
    DOUBLE PRECISION ARG1,ARG2,ARG3,ARG4
    !
    !     THIS IS A FRONT-END TO THE EQUATION SOLVING ROUTINES.
    !     IT'S SOLE PURPOSE IS TO SEGREGATE THE ROUTINES USED TO
    !     SOLVE THE EQUATION FOR THE TEMPERATURE. THE ARGUMENTS ARE
    !     THE COEFFICIENTS OF THE TEMPERATURE EQUATION WITH THE
    !     TEMPERATURE IN EV.
    !
    INTEGER MAXROOTS,MAXINTS
    REAL LOWTEMP,HIGHTEMP,XACC
    PARAMETER (MAXROOTS=7,MAXINTS=200,LOWTEMP=0.0,HIGHTEMP=200.0,XACC=1.0E-3)
    DOUBLE PRECISION FARG1,FARG2,FARG3,FARG4
    COMMON /FPARAMS/ FARG1,FARG2,FARG3,FARG4
    REAL XB1(MAXROOTS),XB2(MAXROOTS)
    INTEGER NB,N,I,J,K
    !REAL SOL14FX,RTBIS
    !EXTERNAL SOL14FX,RTBIS

    FARG1 = ARG1
    FARG2 = ARG2
    FARG3 = ARG3
    FARG4 = ARG4

    NB = MAXROOTS
    N = MAXINTS
    CALL ZBRAK(SOL14FX,LOWTEMP,HIGHTEMP,N,XB1,XB2,NB)

    IF  (NB.LE.0) THEN
       WRITE(6,*) 'SOLVTEMP: ERROR IN EQUATION SOLVER - NO T ROOT FOUND IN RANGE'
       STOP

    ELSEIF (NB.GT.1) THEN

       WRITE(6,*) 'SOLVTEMP: ERROR - MULTIPLE ROOTS FOUND',NB
       DO I = 1, NB
          WRITE(6,*) 'ROOT IN THE FOLLOWING BRACKETS:',XB1(I),XB2(I)
       end do
       !
       !       SOLVE FOR LARGEST AVAILABLE ROOT
       !
       SOLVTEMP = RTBIS(SOL14FX,XB1(NB),XB2(NB),XACC)

    ELSE
       !
       !       SOLVE FOR ROOT
       !
       SOLVTEMP = RTBIS(SOL14FX,XB1(NB),XB2(NB),XACC)
    ENDIF

    RETURN
  END FUNCTION SOLVTEMP



  REAL FUNCTION SOL14FX(T2)
    IMPLICIT NONE
    REAL T2
    !
    !     THIS CALCULATES THE FUNCTIONAL VALUE OF F(T2) WHICH FOR
    !     A ROOT WILL EQUAL ZERO.
    !
    DOUBLE PRECISION FARG1,FARG2,FARG3,FARG4
    COMMON /FPARAMS/ FARG1,FARG2,FARG3,FARG4

    SOL14FX = FARG1*T2 + FARG4*T2**3.5 + FARG3*T2**2.5 + FARG2

    RETURN
  END FUNCTION SOL14FX



  SUBROUTINE ZBRAK(FX,X1,X2,N,XB1,XB2,NB)
    implicit none
    !
    !     GIVEN A FUNCTION FX DEFINED ON THE INTERVAL FROM X1 TO X2, SUBDIVI
    !     THE INTERVAL INTO N EQUALLY SPACED SEGMENTS AND SEARCH FOR ZERO
    !     CROSSINGS OF THE FUNCTION. NB IS INPUT AS THE MAXIMUM NUMBER OF RO
    !     SOUGHT, AND IS RESET TO THE NUMBER OF BRACKETING PAIRS XB1, XB2 TH
    !     ARE FOUND.
    !
    !     FROM "NUMERICAL RECIPES" - CHAPTER 9
    !
    integer n,nb,i,nbb
    EXTERNAL FX
    REAL FX
    REAL XB1(NB),XB2(NB),x1,x2,dx,x,fp,fc

    NBB = NB
    NB = 0
    X = X1
    DX = (X2-X1)/N
    FP = FX(X)
    DO I = 1,N
       X = X+DX
       FC = FX(X)
       IF (FC*FP.LT.0.0) THEN
          NB = NB + 1
          XB1(NB) = X-DX
          XB2(NB) = X
       ENDIF
       FP = FC
       IF (NBB.EQ.NB) RETURN
    end do
    RETURN
  END SUBROUTINE ZBRAK



  REAL FUNCTION RTBIS(FUNC,X1,X2,XACC)
    implicit none
    !
    !     USING BISECTION, FIND THE ROOT OF A FUNCTION FUNC KNOWN
    !     TO LIE BETWEEN X1 AND X2. THE ROOT RETURNED AS RTBIS, WILL BE
    !     REFINED UNTIL ITS ACCURACY IS +/- XACC
    !
    !     FROM "NUMERICAL RECIPES - CHAPTER 9"
    !
    !     LIMIT OF 40 BISECTIONS
    !
    EXTERNAL FUNC
    REAL FUNC,x1,x2,xacc,dx,f,xmid,fmid
    integer jmax,j
    PARAMETER (JMAX=40)


    FMID = FUNC(X2)
    F = FUNC(X1)
    IF (F*FMID.GE.0.0) THEN
       WRITE(6,*) 'ROOT MUST BE BRACKETED FOR BISECTION'
       rtbis = X1
       RETURN
    ENDIF
    IF (F.LT.0.0) THEN
       RTBIS = X1
       DX = X2-X1
    ELSE
       RTBIS = X2
       DX = X1-X2
    ENDIF
    DO J = 1,JMAX
       DX = DX * 0.5
       XMID = RTBIS+DX
       FMID = FUNC(XMID)
       IF (FMID.LE.0.0) RTBIS = XMID
       IF (ABS(DX).LT.XACC.OR.FMID.EQ.0.0) RETURN
    end do
    WRITE (6,*) 'TOO MANY BISECTIONS'
    RETURN
  END FUNCTION RTBIS



  subroutine calcnv(te,ti,gamma,pinf,n,v)
    use mod_comtor
    implicit none
    double precision te,ti,gamma,pinf,n,v
    real :: rizb
    !
    !     The purpose of this routine is to extract the
    !     code that calculates the n,v values since it is
    !     virtually identical for the SOL options avaliable.
    !     This is due to the fact that n,v are calculated
    !     from the pressure balance equation and the
    !     flow equation after Te and Ti have been found.
    !     Since the pressure balance and flow are unaffected
    !     by additional terms in the equations for Te,Ti ...
    !     It seems likely that this code can stay as it
    !     is for some time.
    !
    !     Note: v is returned as -ve only and needs to be
    !     adjusted depending on which side it is on.
    !
    !
    double precision rootn,massi
    !
    MASSI = CRMB * AMU
    rizb = real(cizb)
    !
    !     CALCULATE DENSITY
    !
    ROOTN = (PINF/(ECH*(te+ti)))**2 - 4.0* (MASSI*GAMMA**2) / (ECH*(te+ti))
    !write(0,*) 'ROOTN = ',ROOTN


    IF (SROOTOPT.EQ.0.OR. (SROOTOPT.EQ.1.AND.ROOTN.GE.0.0)) THEN

       IF (ROOTN.LT.0.0) ROOTN = 0.0
       !
       !        Calculate density
       !
       n = PINF/(2.0*ECH*(te+ti)) + 0.5 * SQRT(ROOTN)
       !
       !                 CALCULATE VELOCITY
       !
       v = GAMMA/n
       
       !write(0,*) 'n = ',n
       !write(0,*) 'v = ',v

    ELSEIF (SROOTOPT.EQ.1) THEN
       !
       !        SET VELOCITY TO LOCAL SOUND SPEED.
       !
       v = -SQRT(0.5*EMI*(te+ti)*(1+RIZB)/CRMB)
       !
       !        SET DENSITY
       !
       n = ABS(GAMMA / v)
       !
    ENDIF

    !write(6,'(a,10(1x,g12.5))') 'calnv:',rootn,pinf,te,ti,gamma,n,v, (PINF/(ECH*(te+ti)))**2 ,4.0* (MASSI*GAMMA**2) / (ECH*(te+ti)) ,(PINF/(ECH*(te+ti)))**2 - 4.0* (MASSI*GAMMA**2) / (ECH*(te+ti))

    !      write(6,'(a,10(1x,g12.5))') 'SOLEDGE:NV:',n,v,rootn,te,ti,gamma,pinf

    return
  end subroutine calcnv


end module mod_soledge
