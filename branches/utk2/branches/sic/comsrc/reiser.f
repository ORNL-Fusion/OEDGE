c----------------------------------------------------------------------c
c----------------------------------------------------------------------c

       Subroutine VBGRAD     

c----------------------------------------------------------------------c
c
c   This subroutine calculates the background plasma velocity gradient
c   for each grid cell. It is CALLED in the TAU.D6A files.
c
c----------------------------------------------------------------------c
c     Counters:
c    
c     (IK,IR)= Grid coordinates where IK = Knot no. and IR = Ring no.
c     
c
c----------------------------------------------------------------------c 
c     Input:
c
c     KVHS(IK,IR)   = The background plasma velocity for each cell  [m]
c                     since KVHS(IK,IR)=KVHS(IK,IR)*QTIM in TAU.D6A
c     KBACDS(IK,IR) = The distance from cell(IK-1) to cell(IK)      [m]
c     KFORDS(IK,IR) = The distance from cell(IK) to cell(IK+1)      [m]
c     QTIM          = The quantum time step                         [s]
c     NRS           = The number of rings 
c     NKS(IR)       = The number of knots per ring    
c     IRSEP         = Ring enclosing the separatrix
c
c----------------------------------------------------------------------c
c     Output:
c
c     KVHGS(IK,IR)  = The background plasma velocity gradient [unitless]
c
c----------------------------------------------------------------------c

      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_reiser_com
      IMPLICIT NONE

c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c     include 'reiser_com'  
      
      INTEGER  IK,IR,ID
      REAL     TARGETVB
    
      IK    = 0
      IR    = 0

c----------------------------------------------------------------------c
c     Do SOL and private plasma
         
      DO 200 IR = IRSEP, NRS

         ID = IDDS(IR,2)

         TARGETVB = -SQRT((KTEDS(ID)+KTIDS(ID))*1.602E-19/
     >              (CRMB*1.67E-27))*QTIM
 
         KVHGS(1,IR) = ((KVHS(1,IR)-TARGETVB)/KSS(1,IR) +
     >                 (KVHS(2,IR)-KVHS(1,IR))/KFORDS(1,IR)) * 0.5

c         KVHGS(1,IR) = (KVHS(2,IR)-KVHS(1,IR))/KFORDS(1,IR)

        DO 100 IK = 2,NKS(IR)-1

           KVHGS(IK,IR) = ((KVHS(IK,IR)-KVHS(IK-1,IR))/KBACDS(IK,IR) +
     >      (KVHS(IK+1,IR)-KVHS(IK,IR))/KFORDS(IK,IR)) * 0.5

  100   CONTINUE

         ID = IDDS(IR,1)

         TARGETVB = SQRT((KTEDS(ID)+KTIDS(ID))*1.602E-19/
     >              (CRMB*1.67E-27))*QTIM

          KVHGS(NKS(IR),IR) = 
     >     ((KVHS(NKS(IR),IR)-KVHS(NKS(IR)-1,IR))/KBACDS(NKS(IR),IR) +
     >     (TARGETVB - KVHS(NKS(IR),IR))/(KSMAXS(IR)-KSS(NKS(IR),IR)))
     >     * 0.5

c          KVHGS(NKS(IR),IR) = 
c     >     (KVHS(NKS(IR),IR)-KVHS(NKS(IR)-1,IR))/KBACDS(NKS(IR),IR)

  200  CONTINUE

c----------------------------------------------------------------------c
c     Do main plasma

      DO IR = 1,IRSEP-1 

        KVHGS(1,IR) = ((KVHS(NKS(IR),IR)-KVHS(NKS(IR)-1,IR))
     >                 /KBACDS(NKS(IR),IR) +
     >      (KVHS(2,IR)-KVHS(1,IR))/KFORDS(1,IR)) * 0.5

        DO  IK = 2,NKS(IR)-1
          
           KVHGS(IK,IR) = ((KVHS(IK,IR)-KVHS(IK-1,IR))/KBACDS(IK,IR) +
     >      (KVHS(IK+1,IR)-KVHS(IK,IR))/KFORDS(IK,IR)) * 0.5

        END DO

        KVHGS(NKS(IR),IR) = KVHGS(1,IR)

      END DO
  
      RETURN
      END

c----------------------------------------------------------------------c
c----------------------------------------------------------------------cc

      Subroutine COEFF(NIZS)

c----------------------------------------------------------------------c
c
c     This subroutine calculates the Hermite coefficients required to 
c     evaluate the two dimensional drift and diffusive terms of the 
c     Drift-kinetic Model, see Subroutine Coulomb_Coll below.  These 
c     coefficients will be accessible as global variables.
c
c----------------------------------------------------------------------c
c     Counters:
c    
c     (IK,IR)= Grid coordinates where IK = Knot no. and IR = Ring no.
c
c      IZ    = Impurity charge state
c     
c
c----------------------------------------------------------------------c 
c     Input:
c
c     KNBS(IK,IR)  = Background plasma density for each cell    [m^-3]
c     KTIBS(IK,IR) = Background plasma temperature for each cell[eV]
c     KFIGS(IK,IR) = Background plasma temperature gradient     [m]
c     KVHGS(IK,IR) = Background plasma velocity gradient        [unitless]
c     QTIM         = Quantum time step                          [s]
c     NRS          = Number of rings 
c     NKS(IR)      = Number of knots per ring 
c     IRSEP        = Ring enclosing the separatrix   
c     NIZS         = Number of impurity ionization states
c     CRMI         = Impurity ion mass          [a.m.u.]
c     CRMB         = Background plasma ion mass [a.m.u.]
c     IZ           = Impurity charge state
c     RIZB         = Background plasma charge state
c     Coulomb_log  = 15 (usually)
c
c----------------------------------------------------------------------c
c     Output:
c
c                     coulomb_log*Z^2*Zb^2*e^4
c     LAMBDA1(IZ)   = ------------------------               [m^6/s^4]
c                       4*Pi*epsilon^2*mi^2
c
c     ALPHA(IK,IR)  = SQRT(mb/(2*Tb)       Inverse thermal velocity of 
c                                          background particles [s/m]
c
c     ALPHAI(IK,IR) = SQRT(mb/(2*Tb))/QTIM                      [1/m]
c
c                      6*Pi^2*epsilon^2*k,,*Tb*GradTb
c     CC1(IK,IR)    = --------------------------------       [unitless]
c                     SQRT(Pi)*e^4*coulomb_log*Zb^4*nb
c
c                     -6*Pi^2*epsilon^2*n,,*SQRT(mb*Tb^3)*Gradvb
c     CC2(IK,IR)    = ------------------------------------------     
c                              SQRT(2*Pi)*e^4*Zb^4*nb        [unitless]
c
c     QTIM2         = QTIM * QTIM  Quantum time step (squared)  [s^2]
c
c     MASSFACTOR    = 1.0 + mi/mb                            [unitless]
c
c     PI2SQRT       = 2.0/SQRT(Pi)
c
c     where coulomb_log = 15
c           GradTb      = background temperature gradient
c           Gradvb      = background velocity gradient
c           epsilon     = permittivity of free space 
c           Tb          = background temperature
c           Zb          = background charge state
c           Z           = impurity charge state
c           mb          = background mass
c           mi          = impurity mass
c           nb          = background density
c           e           = electric charge
c           k,,         = ion heat conductivity  1/0.5657
c           n,,         = electron viscosity     1/(1.2+0.8485*Zb^-1) 
c
c     EXP_COEFF(8)  An array containing the expansion coefficients as
c                   required in the Reiser module for small chi limit
c           
c----------------------------------------------------------------------c
c
c     Note: Input values given in the units above will have to be
c           multiplied by various factors (eg. EMI) which are not
c           explicitly shown in the output formulations.  This ensures
c           that the output values will be in the proper units as
c           required by the Reiser subroutine: Coulomb_Coll.
c
c----------------------------------------------------------------------c             
c
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_reiser_com
      IMPLICIT NONE

      INTEGER NIZS

c     include 'params'
c     include 'cgeom'
c     include 'comtor' 
c     include 'reiser_com'  

      INTEGER IK,IR,IZ
      REAL    FACT
 
      DATA IK,IR,IZ,FACT /0,0,0,0.0/

c----------------------------------------------------------------------c
c     Initialization of numerical quantities

      QTIM2 = QTIM * QTIM

      PI2SQRT = 2.0/SQRT(PI)

      MASSFACTOR = 1.0 + CRMI/CRMB
      
c----------------------------------------------------------------------c
c     Expansion coefficients for small chi limit
 
       EXP_COEFF(1) = 1.0
       EXP_COEFF(2) = 2.0/3.0
       EXP_COEFF(3) = 4.0/15.0
       EXP_COEFF(4) = 8.0/105.0
       EXP_COEFF(5) = 16.0/945.0
       EXP_COEFF(6) = 32.0/10395.0
       EXP_COEFF(7) = 64.0/135135.0
       EXP_COEFF(8) = 128.0/2027025.0
      
c----------------------------------------------------------------------c
 
      DO IZ = 1, NIZS

c       LAMBDA1(IZ) = 3.58616*REAL(IZ*IZ)*REAL(RIZB*RIZB)/(CRMI*CRMI)
  
        LAMBDA1(IZ) = 0.23995*Coulomb_log*REAL(IZ*IZ)*REAL(RIZB*RIZB)/
     >                (CRMI*CRMI)           

      END DO

c----------------------------------------------------------------------c
c     Do SOL and Private Plasma

      FACT = QTIM * QTIM * EMI / CRMI

      DO IR = IRSEP, NRS
   
        DO IK = 1, NKS(IR)

          ALPHA(IK,IR) = SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))
 
          ALPHAI(IK,IR) = SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))/QTIM
          
c          CC1(IK,IR) = 1.20246E16*KTIBS(IK,IR)*KFIGS(IK,IR)/
c     >               (REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)*FACT)

          CC1(IK,IR) = 1.8024009E17*KTIBS(IK,IR)*KFIGS(IK,IR)/
     >       (Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)*FACT)

c   Background velocity set to be positive in all instances was again
c   implemented for debugging purposes and is not essential to the code
c
c          KVHGS(IK,IR) = SQRT(KVHGS(IK,IR)*KVHGS(IK,IR))
         
c          CC2(IK,IR) = -4.9113E11*KVHGS(IK,IR)*SQRT(CRMB*KTIBS(IK,IR)*
c     >                  KTIBS(IK,IR)*KTIBS(IK,IR))/
c     >            (RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*KNBS(IK,IR)*QTIM) 

          CC2(IK,IR)=-7.361216E12*KVHGS(IK,IR)*SQRT(CRMB*KTIBS(IK,IR)*
     >             KTIBS(IK,IR)*KTIBS(IK,IR))/(Coulomb_log*
     >             RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*KNBS(IK,IR)*QTIM) 

        END DO

      END DO

c----------------------------------------------------------------------c
c     Do Main Plasma

      DO IR = 1, IRSEP - 1

       DO IK = 1, NKS(IR) - 1

          ALPHA(IK,IR) = SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))
 
          ALPHAI(IK,IR) = SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))/QTIM

c          CC1(IK,IR) = 1.20246E16*KTIBS(IK,IR)*KFIGS(IK,IR)/
c     >               (REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)*FACT)

          CC1(IK,IR) = 1.8024009E17*KTIBS(IK,IR)*KFIGS(IK,IR)/
     >       (Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)*FACT)

c  Same as above

c          KVHGS(IK,IR) = SQRT(KVHGS(IK,IR)*KVHGS(IK,IR))
         
c          CC2(IK,IR) = -4.9113E11*KVHGS(IK,IR)*SQRT(CRMB*KTIBS(IK,IR)*
c     >                  KTIBS(IK,IR)*KTIBS(IK,IR))/
c     >            (RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*KNBS(IK,IR)*QTIM)

          CC2(IK,IR)=-7.361216E12*KVHGS(IK,IR)*SQRT(CRMB*KTIBS(IK,IR)*
     >             KTIBS(IK,IR)*KTIBS(IK,IR))/(Coulomb_log*
     >             RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*KNBS(IK,IR)*QTIM)  

        END DO

        ALPHAI(NKS(IR),IR) = ALPHAI(1,IR)
        CC1(NKS(IR),IR) = CC1(1,IR)
        CC2(NKS(IR),IR) = CC2(1,IR)

      END DO

      RETURN
      END

c----------------------------------------------------------------------c

      Subroutine GradScaleLengthCheck

c----------------------------------------------------------------------c
c
c     This subroutine compares the temperature and parallel plasma
c     velocity gradient scale lengths of each grid cell to a minimum 
c     value calculated from the background plasma quantities assigned 
c     to those cells to determine whether the choice of background 
c     quantities is appropriate for Dirk Reiser's drift-kinetic model.  
c     An alert is set to the screen and to the .dat file when the 
c     criteria fail.
c
c----------------------------------------------------------------------c

      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_reiser_com
      IMPLICIT NONE

      INTEGER NIZS

c     include 'params'
c     include 'cgeom'
c     include 'comtor' 
c     include 'reiser_com'  

      INTEGER IK,IR
      REAL    FACT
      LOGICAL WARNING
 
      DATA IK,IR,FACT,WARNING /0,0,0.0,.FALSE./   

      FACT = QTIM * QTIM * EMI / CRMI

c----------------------------------------------------------------------c
c     Do Main Plasma

c      DO IR = 1, IRSEP - 1

c       DO IK = 1, NKS(IR) - 1

c      Ltb(IK,IR)=ABS(KTIBS(IK,IR)/(KFIGS(IK,IR)/FACT))
c      Lvb(IK,IR)=ABS(SQRT(2.0*KTIBS(IK,IR)*EMI/CRMB)
c     >           /(2.0*(KVHGS(IK,IR)/QTIM)))

c      Ltbmin(IK,IR)=ABS(2.0e27*2.0e27*1.602e-19*1.602e-19*
c     >    (1.0/0.5657)*KTIBS(IK,IR)*KTIBS(IK,IR)
c     >    /(15.0*REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)))
c      Lvbmin(IK,IR)=ABS(2.0e27*2.0e27*1.602e-19*1.602e-19*
c     >    (1.0/(1.2+(0.8485/REAL(RIZB))))
c     >    *KTIBS(IK,IR)*KTIBS(IK,IR)
c     >    /(15.0*REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)))
 
c      IF(Ltb(IK,IR).GE.Ltbmin(IK,IR).AND.Lvb(IK,IR).GE.
c     >   Lvbmin(IK,IR))THEN

c       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'>=',
c     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'>=',
c     >   'Lvbmin:',Lvbmin(IK,IR)

c      ELSEIF(Ltb(IK,IR).LT.Ltbmin(IK,IR).AND.Lvb(IK,IR).GE.
c     >   Lvbmin(IK,IR))THEN

c       WARNING = .TRUE.
   
c       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'<',
c     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'>=',
c     >   'Lvbmin:',Lvbmin(IK,IR)


c      ELSEIF(Ltb(IK,IR).GE.Ltbmin(IK,IR).AND.Lvb(IK,IR).LT.
c     >   Lvbmin(IK,IR))THEN

c       WARNING = .TRUE.
    
c       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'>=',
c     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'<',
c     >   'Lvbmin:',Lvbmin(IK,IR)

c      ELSE

c       WARNING = .TRUE.  

c       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'<',
c     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'<',
c     >   'Lvbmin:',Lvbmin(IK,IR)

c      ENDIF

c       END DO
 
c      END DO

c---------------------------------------------------------------------c
c     Do SOL and Private Plasma

c     Write to the .lim file.

      WRITE(6,*)' '
      WRITE(6,*)'COMPARISON OF THE GRADIENT SCALE LENGTHS,'
     >,' Ltb and Lvb, TO THE MINIMUM VALUE, Ltbmin and Lvbmin:'

      DO IR = IRSEP, NRS
   
        DO IK = 1, NKS(IR)

      Ltb(IK,IR)=ABS(KTIBS(IK,IR)/(KFIGS(IK,IR)/FACT))
      Lvb(IK,IR)=ABS(SQRT(2.0*KTIBS(IK,IR)*EMI/CRMB)
     >           /(2.0*(KVHGS(IK,IR)/QTIM)))

      Ltbmin(IK,IR)=ABS(2.0e27*1.602e-19*1.602e-19*2.0e27*
     >    (1.0/0.5657)*KTIBS(IK,IR)*KTIBS(IK,IR)
     >    /(15.0*REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)))
      Lvbmin(IK,IR)=ABS(2.0e27*1.602e-19*1.602e-19*2.0e27*
     >    (1.0/(1.2+(0.8485/REAL(RIZB))))
     >    *KTIBS(IK,IR)*KTIBS(IK,IR)
     >    /(15.0*REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR)))


c      IF(Ltb(IK,IR).GE.Ltbmin(IK,IR).AND.Lvb(IK,IR).GE.
c     >   Lvbmin(IK,IR))THEN
 
c       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'>=',
c     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'>=',
c     >   'Lvbmin:',Lvbmin(IK,IR)

      IF(Ltb(IK,IR).LT.Ltbmin(IK,IR).AND.Lvb(IK,IR).GE.
     >   Lvbmin(IK,IR))THEN

       WARNING = .TRUE.  

       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'<',
     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'>=',
     >   'Lvbmin:',Lvbmin(IK,IR)


      ELSEIF(Ltb(IK,IR).GE.Ltbmin(IK,IR).AND.Lvb(IK,IR).LT.
     >   Lvbmin(IK,IR))THEN

       WARNING = .TRUE.  
  
       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'>=',
     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'<',
     >   'Lvbmin:',Lvbmin(IK,IR)

      ELSEIF(Ltb(IK,IR).LT.Ltbmin(IK,IR).AND.Lvb(IK,IR).LT.
     >   Lvbmin(IK,IR))THEN

       WARNING = .TRUE.  

       WRITE(6,*)'RING:',IR,'KNOT:',IK,' Ltb:',Ltb(IK,IR),'<',
     >   'Ltbmin:',Ltbmin(IK,IR),' Lvb:',Lvb(IK,IR),'<',
     >   'Lvbmin:',Lvbmin(IK,IR)

      ENDIF

       END DO

      END DO

c----------------------------------------------------------------------c

      IF(WARNING.EQV..TRUE.)THEN

c     Write warning to the screen.

       WRITE(0,*)'WARNING: GRADIENT SCALE LENGTHS ARE LESS THAN',
     >' THE MINIMUM VALUE REQUIRED BY THE DRIFT-KINETIC MODEL!!!'
       CALL PRB

c     Write warning to the .dat file.

       WRITE(7,*)'WARNING: GRADIENT SCALE LENGTHS ARE LESS THAN',
     >' THE MINIMUM VALUE REQUIRED BY THE DRIFT-KINETIC MODEL!!!'
       CALL PRB
      ELSE

c     Write that the criteria has been met to the .lim and .dat files.

       WRITE(6,*)'ALL GRADIENT SCALE LENGTHS ARE EQUAL TO OR LARGER',
     >' THAN THE REQUIRED MINIMUM VALUE OF THE DRIFT-KINETIC MODEL.'
       CALL PRB
       WRITE(7,*)'ALL GRADIENT SCALE LENGTHS ARE EQUAL TO OR LARGER',
     >' THAN THE REQUIRED MINIMUM VALUE OF THE DRIFT-KINETIC MODEL.'
       CALL PRB
      ENDIF

      WRITE(6,*)' '

      RETURN
      END

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c

      Subroutine DECISION(CIOPTION,IKK,IRR)

c----------------------------------------------------------------------c
c
c     This subroutine switches the drift-kinetic model to the fluid
c     approximation model whenever the impurity ion travels into a grid
c     cell for which the gradient scale length is less than the minimum 
c     value as calculated in the subroutine: GradScaleLengthCheck.
c
c----------------------------------------------------------------------c 

      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_reiser_com
      IMPLICIT NONE
      
c     include 'params'
c     include 'cgeom'
c     include 'comtor' 
c     include 'reiser_com'

      INTEGER CIOPTION,IKK,IRR

      CIOPTION = CIOPTR
 
      IF((Ltb(IKK,IRR).LT.Ltbmin(IKK,IRR).OR.
     >   Lvb(IKK,IRR).LT.Lvbmin(IKK,IRR)).AND.COMBINE.EQ.1)
     >  CIOPTION = 0
       
      RETURN
      END

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c       

      Subroutine Coulomb_Coll(Kpara,Dparapara,Xa,Lambda,
     >                        IKK,IRR,K11,K12,K13,D11,D12,D13)

c      Subroutine Coulomb_Coll(Koeff,DiffEW,DiffEV,
c     >     Xa,Xb,C0,C1,C2,Lambda,Alpha,Mit,Mib,
c     >     XX,Phi00,Phi01,Phi10,Phi11,Phi20,Phi21,     
c     >     Psi00,Psi01,Psi02,Psi10,Psi11,Psi12,Psi20,Psi21,Psi22, 
c     >     A01,A11,A12,A21,A22,B01,B11,B12,B21,B22,
c     >     C01,C02,C11,C12,C21,C22,C23,D1Phi,D1Psi,D2Psi)

c     ------------------------------------------------------------     *
c                                                                      *
c     Original code contributed by Dirk Reiser Feb, 1999.              *
c                                                                      *
c     Modifications were made and sections were commented out to       *
c     reduce the CPU time since only the drift and diffusive terms     *
c     in the parallel-to-B direction are being considered.             *
c                                                                      *
c     This subroutine calculates the Trubnikov potentials and their    *
c     derivatives determining the two dimensional drift kinetic        *
c     Coulomb collision term.                                          *
c     The equation numbers below refer to                              *
c     "Improved Kinetic Test Particle Model for Inpurity Transport     *
c      in Tokamaks" by D.Reiser, D.Reiter & M.Z.Tokar                  *
c     Nucl.Fusion Vol.38, No.2 (1998)                                  *
c     ------------------------------------------------------------     *
c     Input:                                                           *
c                                                                      *
c     Xa       = Normalized velocity chi parallel to B-field           * eq. (36a)
c     Xb       = Normalized velocity chi perpendicular to B-field      * eq. (36a)
c     C0       = Multiplying factor (usually = 1.0)                    * 
c     C1       = Multiplying factor                                    * eq. (45a)
c     C2       = Multiplying factor                                    * eq. (45b)
c     Lambda   = Coulomb-Logarithm                                     * 
c     Alpha    = Inverse thermal velocity of background particles      * eq. (36b)
c     Mit      = Impurity mass [kg]                                    *
c     Mib      = Background ion mass [kg]                              *
c                                                                      *
c     ------------------------------------------------------------     *
c     Output:                                                          *
c                                                                      *
c     Koeff(1) = K_parallel                                            * eq. (32a)
c     Koeff(2) = K_perp                                                * eq. (32b)
c     Koeff(3) = D_parallel,parallel                                   * eq. (32c)
c     Koeff(4) = D_parallel,perp                                       * eq. (32d)
c     Koeff(5) = D_perp,parallel                                       * eq. (32e)
c     Koeff(6) = D_perp,perp                                           * eq. (32f)
c     DiffEW   = Eigenvalues of tensor                                 * 
c     DiffEV   = Eigenvectors of tensor D                              * 
c     XX       = Normalized velocity chi                               *
c     Phi00    = Potential function phi0                               * eq. (46a) 
c     Phi01    = First derivative of phi0                              * eq. (48a) 
c     Phi10    = Potential function phi1                               * eq. (46b) 
c     Phi11    = First derivative of phi1                              * eq. (48b) 
c     Phi20    = Potential function phi2                               * eq. (46c) 
c     Phi21    = First derivative of phi2                              * eq. (48c) 
c     Psi00    = Potential function psi0                               * eq. (46d) 
c     Psi01    = First derivative of psi0                              * eq. (49a)
c     Psi02    = Second derivative of psi0                             * eq. (50a)
c     Psi10    = Potential function psi1                               * eq. (46e) 
c     Psi11    = First derivative of psi1                              * eq. (49b)
c     Psi12    = Second derivative of psi1                             * eq. (50b)
c     Psi20    = Potential function psi2                               * eq. (46f) 
c     Psi21    = First derivative of psi2                              * eq. (49c)
c     Psi22    = Second derivative of psi2                             * eq. (50c)
c     A01      = phi01                -> matrix component (A0)_1       * eq. (53a)
c     A11      = phi11                -> matrix component (A1)_1       * eq. (54a)
c     A12      = phi10/chi            -> matrix component (A1)_2       * eq. (54b)
c     A21      = phi21/chi            -> matrix component (A2)_1       * eq. (55a) 
c     A22      = phi20/chi            -> matrix component (A2)_2       * eq. (55b)
c     B01      = psi01                -> matrix component (B0)_1       * eq. (56a)
c     B11      = psi11                -> matrix component (B1)_1       * eq. (57a)
c     B12      = psi10/chi            -> matrix component (B1)_2       * eq. (57b)
c     B21      = psi21                -> matrix component (B2)_1       * eq. (58a)
c     B22      = psi20/chi            -> matrix component (B2)_2       * eq. (58b)
c     C01      = psi02                -> matrix component (C0)_11      * eq. (59a)
c     C02      = psi01/chi            -> matrix component (C0)_22      * eq. (59c)
c     C11      = psi12                -> matrix component (C1)_11      * eq. (60a)
c     C12      = psi11/chi-psi10/chi  ->  (C1)_12 & (C1)_22            * eqs. (60b,60c)
c     C21      = psi22                -> matrix component (C2)_11      * eq. (61a)
c     C22      = psi21/chi            ->  (C2)_12 & (C2)_22            * eqs. (61b,61c)
c     C23      = psi20/chi^2          ->  (C2)_12 & (C2)_22            * eqs. (61b,61c) 
c     D1Phi(1) = Derivative of phi with respect to chi_parallel        * eq. (51a)
c     D1Phi(2) = Derivative of phi with respect to chi_perp            * eq. (51a)
c     D1Psi(1) = Derivative of psi with respect to chi_parallel        * eq. (51b)
c     D1Psi(2) = Derivative of psi with respect to chi_perp            * eq. (51b)
c     D2Psi(1,1) = Second derivative of psi with respect to            *
c                  chi_parallel                                        * eq. (51c)
c     D2Psi(1,2) = Second derivative of psi with respect to            *
c                  chi_parallel and chi_perp                           * eq. (51c)
c     D2Psi(2,1) = Second derivative of psi with respect to            *
c                  chi_perp and chi_parallel                           * eq. (51c)
c     D2Psi(2,2) = Second derivative of psi with respect to            *
c                  chi_perp                                            * eq. (51c)
c                                                                      *
c     ------------------------------------------------------------     *
c     Used Functions: Erf                                              *
c     ------------------------------------------------------------     *

      use mod_params
      use mod_comtor
      use mod_reiser_com
      IMPLICIT NONE
      
      INTEGER IKK,IRR,IZZ,NIZZS
      REAL    Lambda,Xa,VEL,VEL1
      REAL    Kpara,Dparapara
      REAL    K11,K12,K13,D11,D12,D13

c     include 'params'
c     include 'comtor' 
c     include 'reiser_com'  

      INTEGER i
      
      REAL AlphaL,Alpha2L,CosX
      REAL Erf0(0:5),Erf1(0:5),X(14)
      REAL XX,Phi01,Phi11,Phi21
      REAL Psi02,Psi12,Psi22
      REAL A01,A11,A21
      REAL C01,C11,C21
      REAL D1Phi(2),D2Psi(2,2)

      DOUBLE PRECISION SEED

      REAL Erf

c Variable declarations from original Code Source kept as Reference
c
c      INTEGER i
c      REAL Mit,Mib,Lambda,Alpha,MassFactor
c      REAL AlphaL,Alpha2L,Xa,Xb,C0,C1,C2
c      REAL CosX,Cos2X,SinX,Sin2X
c      REAL Erf0(0:5),Erf1(0:5),X(14),A(8)
c      REAL XX,Phi00,Phi01,Phi10,Phi11,Phi20,Phi21
c      REAL Psi00,Psi01,Psi02,Psi10,Psi11,Psi12,Psi20,Psi21,Psi22
c      REAL A01,A11,A12,A21,A22,B01,B11,B12,B21,B22
c      REAL C01,C02,C11,C12,C21,C22,C23
*      REAL Koeff(6),DiffEW(2),DiffEV(2,2),D1Phi(2),D1Psi(2),D2Psi(2,2)
*      REAL Erf,aa,bb,cc,dummy1,dummy2

      EXTERNAL Erf

c      REAL Pi
c      PARAMETER( Pi=3.1415926 )
c      REAL SqrtPi
c      PARAMETER( SqrtPi=1.77724539)

c     ------------------------------------------------------------     *
c     Expansion coefficients for small chi limit
c     ------------------------------------------------------------     *

c      A(1) = 1.0                             ! Made Global in Coeff.d6a
c      A(2) = 2.0/3.0
c      A(3) = 4.0/15.0
c      A(4) = 8.0/105.0
c      A(5) = 16.0/945.0
c      A(6) = 32.0/10395.0
c      A(7) = 64.0/135135.0
c      A(8) = 128.0/2027025.0

c     ------------------------------------------------------------     *
c     Initialize some numerical quantities and direction angles
c     ------------------------------------------------------------     *

c      MassFactor = 1.0 + Mit/Mib             ! Made Global in Coeff

      AlphaL	 = ALPHA(IKK,IRR)*Lambda
      Alpha2L	 = ALPHA(IKK,IRR)*AlphaL     

c      Xb = Max(1.0E-4,Xb)		      ! Division by Xb below

      XX = Sqrt(Xa*Xa)

      Erf0(0) = Erf(XX)			      ! Error-Function
c      Erf1(0) = 2.0/SqrtPi * Exp(-(XX*XX))   ! Derivative of Error-Function
      Erf1(0) = PI2SQRT * Exp(-(XX*XX))

      if (XX.ne.0.0) then
      	CosX  = Xa/XX
c        Cos2X = CosX*CosX
c	SinX  = Xb/XX
c        Sin2X = SinX*SinX
      else
        CosX  = 1.0
c        Cos2X = 1.0
c        SinX  = 0.0
c        Sin2X = 0.0
      end if

c     ------------------------------------------------------------     *

c     Only those terms pertaining to the motion parallel to the        *
c     B-field line have been retained. Terms corresponding to the      *
c     perpendicular motion have been commented out to increase         *
c     efficiency. To consider perpendicular motion in future           *
c     applications, simply 'uncomment' the relevant lines in this      *
c     subroutine and change the expansion coefficients A(8) to         *
c     EXP_COEFF(8).                                                    *

c     ------------------------------------------------------------     *
c     Computation for large chi
c     ------------------------------------------------------------     *

      if (XX.gt.1.0) then       

	X(1) = XX
	do i=2,5
	  X(i) = X(i-1)*XX
	end do
	do i=1,5
	  Erf0(i) = Erf0(0)/X(i)
	  Erf1(i) = Erf1(0)/X(i)
	end do

c	Phi00 = Erf0(1)
c	Phi01 = Erf1(1) - Erf0(2)
c	Phi10 = X(1)*Erf1(0)
c	Phi11 = Erf1(0) - 2.0*X(2)*Erf1(0)
c	Phi20 = 2.0*Erf1(0)/3.0 + Erf1(2) - Erf0(3)
c	Phi21 = -4.0*X(1)*Erf1(0)/3.0 - 2.0*Erf1(1) - 3.0*Erf1(3) +
c    >		 3.0*Erf0(4)
c	Psi00 = Erf1(0)/2.0 + Erf0(1)/2.0 + X(1)*Erf0(0)
c	Psi01 = Erf1(1)/2.0 + Erf0(0) - Erf0(2)/2.0
c	Psi02 = -Erf1(2) + Erf0(3)
c	Psi10 = Erf1(1)/2.0 - Erf0(2)/2.0
c	Psi11 = -Erf1(0) - Erf1(2) + Erf0(3)
c	Psi12 = 2.0*X(1)*Erf1(0) + 2.0*Erf1(1) + 3.0*Erf1(3) -
c    >		 3.0*Erf0(4)
c	Psi20 = Erf1(2)/2.0 + Erf0(1)/3.0 - Erf0(3)/2.0
c	Psi21 = -2.0*Erf1(1)/3.0 - 3.0*Erf1(3)/2.0 - Erf0(2)/3.0 +
c    >		 3.0*Erf0(4)/2.0
c	Psi22 = 4.0*Erf1(0)/3.0	+ 10.0*Erf1(2)/3.0 + 6.0*Erf1(4) +
c    >		 2.0*Erf0(3)/3.0 - 6.0*Erf0(5)

        A01  = Erf1(1) - Erf0(2)
        A11  = CC1(IKK,IRR)*(Erf1(0) - 2.0*X(2)*Erf1(0))
        A21  = CC2(IKK,IRR)*(-(4.0*X(1)*Erf1(0)/3.0) - 2.0*Erf1(1) -
     >          3.0*Erf1(3) + 3.0*Erf0(4))

c	A01  = C0*(Erf1(1) - Erf0(2))					! Phi01
c	A11  = C1*(Erf1(0) - 2.0*X(2)*Erf1(0))				! Phi11
c	A12  = C1*Erf1(0)						! Phi10/X
c	A21  = C2*(-(4.0*X(1)*Erf1(0)/3.0) - 2.0*Erf1(1) - 3.0*Erf1(3) +! Phi21
c     >		3.0*Erf0(4))
c	A22  = C2*(2.0*Erf1(1)/3.0 + Erf1(3) - Erf0(4))			! Phi20/X
c            
c	B01  = C0*(Erf1(1)/2.0 + Erf0(0) - Erf0(2)/2.0)			! Psi01
c	B11  = C1*(-Erf1(0) - Erf1(2) + Erf0(3))			! Psi11
c	B12  = C1*(Erf1(2)/2.0 - Erf0(3)/2.0)				! Psi10/X
c	B21  = C2*(-(2.0*Erf1(1)/3.0) - 3.0*Erf1(3)/2.0 - Erf0(2)/3.0 +	! Psi21
c     >		3.0*Erf0(4)/2.0)
c	B22  = C2*(Erf1(3)/2.0 + Erf0(2)/3.0 - Erf0(4)/2.0)		! Psi20/X

        C01  = -Erf1(2) + Erf0(3)
        C11  = CC1(IKK,IRR)*(2.0*X(1)*Erf1(0) + 2.0*Erf1(1) +
     >          3.0*Erf1(3) - 3.0*Erf0(4))
       	C21  = CC2(IKK,IRR)*(4.0*Erf1(0)/3.0 + 10.0*Erf1(2)/3.0 +
     >          6.0*Erf1(4) + 2.0*Erf0(3)/3.0 - 6.0*Erf0(5))

c	C01  = C0*(-Erf1(2) + Erf0(3))					! Psi02
c	C02  = B01/X(1)							! Psi01/X
c	C11  = C1*(2.0*X(1)*Erf1(0) + 2.0*Erf1(1) + 3.0*Erf1(3) -	! Psi12
c     >		3.0*Erf0(4))
c	C12  = B11/X(1) - B12/X(1)					! Psi11/X-Psi10/X2
c	C21  = C2*(4.0*Erf1(0)/3.0 + 10.0*Erf1(2)/3.0 + 6.0*Erf1(4) +	! Psi22
c     >		2.0*Erf0(3)/3.0 - 6.0*Erf0(5))
c	C22  = B21/X(1)							! Psi21/X
c	C23  = B22/X(1)							! Psi20/X2

      end if

*     ------------------------------------------------------------     *
*     Computation for small chi 
*     (the formulas in comments are for completeness and not necessary here)
*     ------------------------------------------------------------     *

      if (XX.le.1.0) then

	X(1) = XX
	do i=2,14
	  X(i) = X(i-1)*XX
	end do

	Phi01 = 0.0
	do i=0,6
	  Phi01 = Phi01 - EXP_COEFF(i+2)*X(2*i+1)
	end do
	Phi01 = Erf1(0)*Phi01

	Phi11 = Erf1(0)*(1.0-2.0*X(2))

	Phi21 = 0.0
	do i=0,5
	  Phi21 = Phi21 + 3.0*EXP_COEFF(i+3)*X(2*i+1)
	end do
        Phi21 = Erf1(0)*(-(4.0*X(1)/3.0)+Phi21)

c
c       Phi21 = Erf1(0)*(-(2.0*X(1)/3.0)+Phi21)
c
c	Phi00 = 0.0
c	do i=0,6
c	  Phi00 = Phi00 + A(i+2)*X(2*i+2)
c	end do
c	Phi00 = Erf1(0)*(1.0 + Phi00)
c	Phi01 = 0.0
c	do i=0,6
c	  Phi01 = Phi01 - A(i+2)*X(2*i+1)
c	end do
c	Phi01 = Erf1(0)*Phi01
c	Phi10 = Erf1(0)*X(1)
c	Phi11 = Erf1(0)*(1.0-2.0*X(2))
c	Phi20 = 0.0
c	do i=0,5
c	  Phi20 = Phi20 - A(i+3)*X(2*i+2)
c	end do
c	Phi20 = Erf1(0)*Phi20
c	Phi21 = 0.0
c	do i=0,5
c	  Phi21 = Phi21 + 3.0*A(i+3)*X(2*i+1)
c	end do
c	Phi21 = Erf1(0)*(-(4.0*X(1)/3.0)+Phi21)
c	Psi00 = 0.0
c	do i=0,6
c	  Psi00 = Psi00 + (A(i+2)/2.0+A(i+1))*X(2*i+2)
c	end do
c	Psi00 = Erf1(0)*(1.0 + Psi00)
c	Psi01 = 0.0
c	do i=0,6
c	  Psi01 = Psi01 + (A(i+1)-A(i+2)/2.0)*X(2*i+1)
c	end do
c	Psi01 = Erf1(0)*Psi01

	Psi02 = 0.0
	do i=0,5
	  Psi02 = Psi02 + EXP_COEFF(i+3)*X(2*i+2)
	end do
	Psi02 = Erf1(0)*(2.0/3.0+Psi02)

	Psi12 = 0.0
	do i=0,5
	  Psi12 = Psi12 - 3.0*EXP_COEFF(i+3)*X(2*i+1)
	end do
	Psi12 = Erf1(0)*(2.0*X(1) + Psi12)
        
        Psi22 = 0.0
	do i=0,4
	  Psi22 = Psi22 + (2.0*EXP_COEFF(i+3)/3.0 - 
     >            6.0*EXP_COEFF(i+4))*X(2*i+2)
	end do
	Psi22 = Erf1(0)*(8.0/45.0 + Psi22)

c	Psi02 = 0.0
c	do i=0,5
c	  Psi02 = Psi02 + A(i+3)*X(2*i+2)
c	end do
c	Psi02 = Erf1(0)*(2.0/3.0+Psi02)
c	Psi10 = 0.0
c	do i=0,6
c	  Psi10 = Psi10 - A(i+2)*X(2*i+1)/2.0
c	end do
c	Psi10 = Erf1(0)*Psi10
c	Psi11 = 0.0
c	do i=0,5
c	  Psi11 = Psi11 + A(i+3)*X(2*i+2)
c	end do
c	Psi11 = Erf1(0)*(-(1.0/3.0) + Psi11)
c	Psi12 = 0.0
c	do i=0,5
c	  Psi12 = Psi12 - 3.0*A(i+3)*X(2*i+1)
c	end do
c	Psi12 = Erf1(0)*(2.0*X(1) + Psi12)
c	Psi20 = 0.0
c	do i=0,5
c	  Psi20 = Psi20 + (A(i+2)/3.0 - A(i+3)/2.0)*X(2*i+2)
c	end do
c	Psi20 = Erf1(0)*Psi20
c	Psi21 = 0.0
c	do i=0,5
c	  Psi21 = Psi21 + (-(A(i+2)/3.0) + 3.0*A(i+3)/2.0)*X(2*i+1)
c	end do
c	Psi21 = Erf1(0)*Psi21
c	Psi22 = 0.0
c	do i=0,4
c	  Psi22 = Psi22 + (2.0*A(i+3)/3.0 - 6.0*A(i+4))*X(2*i+2)
c	end do
c	Psi22 = Erf1(0)*(8.0/45.0 + Psi22)

        A01  = Phi01
        A11  = CC1(IKK,IRR)*Phi11
        A21  = CC2(IKK,IRR)*Phi21

c       A01  = C0*Phi01
c	A11  = C1*Phi11
c	A12  = C1*Erf1(0)
c	A21  = C2*Phi21
c	A22  = 0.0
c	do i=0,5
c	  A22 = A22 - A(i+3)*X(2*i+1)
c	end do
c	A22  = C2*Erf1(0)*A22

c	B01  = C0*Psi01
c	B11  = C1*Psi11
c	B12  = 0.0
c	do i=0,5
c	  B12 = B12 - A(i+3)*X(2*i+2)/2.0
c	end do
c	B12  = C1*Erf1(0)*(-(1.0/3.0) + B12)
c	B21  = C2*Psi21
c	B22  = 0.0
c	do i=0,5
c	  B22 = B22 + (A(i+2)/3.0 - A(i+3)/2.0)*X(2*i+1)
c	end do
c	B22  = C2*Erf1(0)*B22

        C01  = Psi02
       	C11  = CC1(IKK,IRR)*Psi12
        C21  = CC2(IKK,IRR)*Psi22

c	C01  = C0*Psi02
c	C02  = 0.0
c	do i=0,5
c	  C02 = C02 + (A(i+2)-A(i+3)/2.0)*X(2*i+2)
c	end do
c	C02  = C0*Erf1(0)*(2.0/3.0 + C02)
c	C11  = C1*Psi12
c	C12  = 0.0
c	do i=0,5
c	  C12 = C12 + 3.0*A(i+3)*X(2*i+1)/2.0
c	end do
c	C12  = C1*Erf1(0)*C12
c	C21  = C2*Psi22
c	C22 = 0.0
c	do i=0,4
c	  C22 = C22 + (-(A(i+3)/3.0) + 3.0*A(i+4)/2.0)*X(2*i+2)
c	end do
c	C22 = C2*Erf1(0)*(8.0/45.0 + C22)
c	C23 = 0.0
c	do i=0,4
c	  C23 = C23 + (A(i+3)/3.0 - A(i+4)/2.0)*X(2*i+2)
c	end do
c	C23 = C2*Erf1(0)*(4.0/45.0 + C23)

      end if
 
c     ------------------------------------------------------------     *
c     Now sum up all parts of the drift and diffusion coefficients
c     ------------------------------------------------------------     *

c     dummy1 = 1.0-3.0*Cos2X 
c     dummy2 = CosX*SinX

c      D1Phi(1)   = CosX*A01 + Cos2X*A11 + Sin2X*A12 + 
c     >             dummy1*CosX*A21 - 6.0*dummy2*SinX*A22
c      D1Phi(2)   = SinX*A01 + dummy2*A11 - dummy2*A12 + 
c     >             dummy1*SinX*A21 - 6.0*dummy2*CosX*A22
c      D1Psi(1)   = CosX*B01 + Cos2X*B11 + Sin2X*B12 + 
c     >             dummy1*CosX*B21 - 6.0*dummy2*SinX*B22
c      D1Psi(2)   = SinX*B01 + dummy2*B11 - dummy2*B12 + 
c     >             dummy1*SinX*B21 - 6.0*dummy2*CosX*B22

c      D2Psi(1,1) = Cos2X*C01 + Sin2X*C02 
c     >		  + CosX*Cos2X*C11 + 3.0*Sin2X*CosX*C12 +
c     >		   dummy1*Cos2X*C21 - 12.0*Sin2X*Cos2X*C22 +
c     >		   dummy1*Sin2X*C22 + 12.0*Sin2X*Cos2X*C23 -
c     >		   6.0*(1.0-2.0*Cos2X)*Sin2X*C23
c      D2Psi(1,2) = dummy2*C01 - dummy2*C02 
c     >		  + SinX*Cos2X*C11 + dummy1*SinX*C12 +
c     >		   dummy1*dummy2*C21 - 6.0*dummy2*(1.0-2.0*Cos2X)*C22 -
c     >		   dummy1*dummy2*C22 + 12.0*(1.0-2.0*Cos2X)*dummy2*C23
c      D2Psi(2,1) = D2Psi(1,2)
c      D2Psi(2,2) = Sin2X*C01 + Cos2X*C02 
c     >		  + CosX*Sin2X*C11 + dummy1*CosX*C12 +
c     >		   dummy1*Sin2X*C21 + 12.0*Sin2X*Cos2X*C22 +
c     >		   dummy1*Cos2X*C22 - 12.0*Sin2X*Cos2X*C23 -
c     >		   6.0*(1.0-2.0*Cos2X)*Cos2X*C23
c
      
c      Koeff(1) = Alpha2L * MassFactor * D1Phi(1)
c      Koeff(2) = Alpha2L * MassFactor * D1Phi(2) 
c     >		 + 0.5 * Alpha2L * D1Psi(2)/(Xb*Xb)
c      Koeff(3) = AlphaL * D2Psi(1,1)
c      Koeff(4) = AlphaL * D2Psi(1,2)
c      Koeff(5) = AlphaL * D2Psi(2,1)
c      Koeff(6) = AlphaL * D2Psi(2,2)

       K11 = Alpha2L * MassFactor * CosX*A01 
       K12 = Alpha2L * MassFactor * A11
       K13 = Alpha2L * MassFactor * (-2.0*CosX*A21)
 
       D11 =  AlphaL * C01
       D12 =  AlphaL * CosX*C11
       D13 =  AlphaL * (-2.0*C21)

c      Aswitch allows for the specific selection of
c      individual collision terms for K|| and D|| ||
c      -Selection is made in the input data file.

       IF (Aswitch.eq.1) then
        if (sk11.eq.0) k11 = 0.0
	if (sk12.eq.0) k12 = 0.0
	if (sk13.eq.0) k13 = 0.0
	if (sd11.eq.0) d11 = 0.0
	if (sd12.eq.0) d12 = 0.0
	if (sd13.eq.0) d13 = 0.0
       ENDIF	

       Kpara = K11 + K12 + K13

       Dparapara = D11 + D12 + D13

c      A check is implemented at this point to ensure that 
c      a negative value of the diffusion coefficient is not
c      returned to the calling subroutine.  Thus one does
c      not have to contend with a negative square root.
c      
       IF (Dparapara.lt.0.0) THEN

c          write (6,*) 'REISER ERROR: Dparapara(',ikk,',',irr,') < 0',
c     >                 dparapara,d11,d12,d13,Xa  
c          write (0,*) 'REISER ERROR: Dparapara < 0',
c     >                ' Dparapara=',dparapara,' D11=',
c     >      d11,' D12=',d12,' D13=',d13,' IK=',ikk,' IR=',irr,' CHI=',Xa 

c      Since D11(CHIpara = 0) = [vz(thermal)/tau(parallel)]^2 and since
c      D11 is always positive, let it be the value of Dparapara when 
c      Dparapara would otherwise be negative.
   
c          Dparapara = 0.0

           Dparapara = D11

       ENDIF

c     ------------------------------------------------------------     *
c     Compute Eigenvalues and Eigenvectors of diffusion matrix
c     ------------------------------------------------------------     *

c      aa = Koeff(3)
c      bb = Koeff(4)
c      cc = Koeff(6)

c      DiffEW(1) = 0.5*(aa+cc) + Sqrt(0.25*(aa-cc)*(aa-cc)+bb*bb)
c      DiffEW(2) = 0.5*(aa+cc) - Sqrt(0.25*(aa-cc)*(aa-cc)+bb*bb)

c      if (bb.ne.0.0) then
c	dummy1 = bb*bb + (DiffEW(1)-cc)*(DiffEW(1)-cc)
c	dummy1 = Sqrt(dummy1)
c	dummy2 = bb*bb + (DiffEW(2)-cc)*(DiffEW(2)-cc)
c	dummy2 = Sqrt(dummy2)
c	DiffEV(1,1) = (DiffEW(1)-cc)/dummy1
c	DiffEV(1,2) = bb/dummy1
c	DiffEV(2,1) = (DiffEW(2)-cc)/dummy2
c	DiffEV(2,2) = bb/dummy2
c	return
c      end if

c      if (bb.eq.0.0) then
c	DiffEV(1,1) = 1.0
c	DiffEV(1,2) = 0.0
c	DiffEV(2,1) = 0.0
c	DiffEV(2,2) = 1.0
c	return
c      end if

      return
      end

c     ------------------------------------------------------------     *

      REAL Function Erf(X)

c     ------------------------------------------------------------     *
c     --Error-Function for 0 <= x <= infinity --------------------     *
c     --error less than 2.5E-5 for all values of x ---------------     *
c     --Abramowicz 7.1.25 ----------------------------------------     *
c     ------------------------------------------------------------     *

      IMPLICIT NONE
      REAL X,A1,A2,A3,T,T2,T3,P
      DATA A1/0.34802/, A2/-0.09587/
      DATA A3/0.74785/, P/0.47047/
      T   =  1.0/(1.0 + P*X)
      T2  =  T*T
      T3  =  T*T2
      Erf = 1.0 - (A1*T+A2*T2+A3*T3)*Exp(-(X*X))
      return
      end

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c

      Subroutine interp_quant(ik,ir,s,vals,targvals,value,grad,scalef)

c----------------------------------------------------------------------c
c
c     This subroutine interpolates the values of the plasma quantities 
c     between the cell centers of a given cell and its nearest neighbour.
c     The value of a plasma quantity is assumed to be averaged over the 
c     entire cell.  It is further assumed that this value is then exact 
c     at the cell's center and can therefore be used to determine the 
c     values between the cell centers.  The plasma quantities are passed
c     in from the subroutine: update_reiser_coeff.
c
c----------------------------------------------------------------------c  

      use mod_params
      use mod_cgeom
      implicit none
c     include 'params' 
c     include 'cgeom' 
      integer ik,ir
      real s,vals(maxnks,maxnrs),value,grad,scalef
      real targvals(maxnds)

c----------------------------------------------------------------------c
c      
c     Linearly interpolate the given DIVIMP background quantity
c         
      integer id
c     
      if (ik.eq.1.and.s.le.kss(ik,ir)) then 
c         
         id = idds(ir,2)
       
         if (ir.ge.irsep) then 
      
            if (kss(ik,ir).ne.0.0) then 
               grad = (scalef*vals(ik,ir)-targvals(id))/kss(ik,ir)
               value = targvals(id) + s * grad
            else
               value = targvals(id) 
               grad = 0.0
            endif
c     
         else
c     
            value = scalef*vals(ik,ir)
            grad = 0.0
c     
         endif
c     
      elseif (ik.eq.nks(ir).and.s.ge.kss(ik,ir)) then  
c     
         id = idds(ir,1)
      
         if (ir.ge.irsep) then   
      
            if (ksmaxs(ir).ne.kss(ik,ir)) then 
      
               grad  =  (targvals(id)-scalef*vals(ik,ir))/
     >                  (ksmaxs(ir)-kss(ik,ir))
               value = scalef*vals(ik,ir) + (s-kss(ik,ir)) * grad 
            else
               value = targvals(id) 
               grad = 0.0
            endif 
c     
        else
c     
           value = scalef*vals(ik,ir)
           grad  = 0.0
c     
        endif 
c     
      elseif (s.ge.kss(ik,ir)) then 
      
         grad = (scalef*vals(ik+1,ir)-scalef*vals(ik,ir))/
     >           (kss(ik+1,ir)-kss(ik,ir))
      
         value = scalef*vals(ik,ir) + (s-kss(ik,ir)) * grad
      
           
      elseif (s.le.kss(ik,ir)) then  
      
         grad =  (scalef*vals(ik,ir)-scalef*vals(ik-1,ir))/
     >            (kss(ik,ir)-kss(ik-1,ir))
      
         value = scalef*vals(ik-1,ir) + (s-kss(ik-1,ir)) * grad 
      
      endif
c    
      return
      end

c-----------------------------------------------------------------c
c-----------------------------------------------------------------c

      Subroutine update_reiser_coeff(ik,ir,iz,s,lambda2,vbqtim)

c-----------------------------------------------------------------c
c
c     This routine updates the Reiser coefficients at each time 
c     step for the specific position of the particle - it uses 
c     a linear interpolation between the grid centers along S 
c     to estimate interpolated values and gradients.  This is 
c     achieved by calling the subroutine: interp_quant. 
c
c-----------------------------------------------------------------c

      use mod_params
      use mod_reiser_com
      use mod_cgeom
      use mod_comtor
      implicit none
c     include 'params'
c     include 'reiser_com'  
c     include 'cgeom'
c     include 'comtor' 
      integer ik,ir,iz,id
      real s, lambda2, vbqtim,targetvb

c
      real tb,tslope,vb,vslope,nb,nslope,INTERCEPT,MIDPOINT

c-----------------------------------------------------------------c
c     
c      call interp_quant(ik,ir,s,ktibs,ktids,tb,tslope,1.0)
      call interp_quant(ik,ir,s,kvhs,kvds,vb,vslope,1.0/QTIM)
      call interp_quant(ik,ir,s,knbs,knds,nb,nslope,1.0)

c     Special case where one has the occurence of a discontinity at the 
c     midpoint when selecting SOL opt. 7 to create a temp. grad.
c     profile which rises linearly from the target to the midpoint and 
c     then descends linearly back to the other target. 

      IF (CIOPTF.EQ.7.and.linearpeak.EQ.1) THEN 

          MIDPOINT = KSMAXS(IR)/2

          IF (S.LE.MIDPOINT) THEN
           TSLOPE = KFIGS(2,IR)*CRMI/(QTIM2*EMI)
           INTERCEPT = CFIBT * CTIB0                       !outer target
          ELSE         
           TSLOPE = KFIGS(nks(ir)-1,IR)*CRMI/(QTIM2*EMI)
           INTERCEPT = (CFIBT * CTIB0)-(TSLOPE*KSMAXS(IR)) !inner target
          ENDIF

c          VSLOPE = KVHGS(2,IR)/QTIM
c          VINTERCEPT = CVHOUT

          TB = (Tslope * S) + INTERCEPT
c          VB = (vslope * S) + VINTERCEPT          
      
      ELSE 
         
          call interp_quant(ik,ir,s,ktibs,ktids,tb,tslope,1.0)

      ENDIF

c-----Evaluation of velocity at the targets based on target temperature

      if (ik.eq.1.and.s.lt.kss(1,ir).and.ir.ge.irsep)then 
         id = idds(ir,2)
         targetvb = -SQRT((KTEDS(id)+KTIDS(id))*1.602E-19/
     >              (CRMB*1.67E-27))
         vslope = ((KVHS(ik,ir)/QTIM)-targetvb)/kss(ik,ir)
         vb = targetvb + s * vslope
c         write(0,*)kteds(id),ktids(id),targetvb,vslope,vb

      elseif (ik.eq.nks(ir).and.s.ge.kss(ik,ir).and.ir.ge.irsep) then       
         id = idds(ir,1)
         targetvb = SQRT((KTEDS(id)+KTIDS(id))*1.602E-19/
     >              (CRMB*1.67E-27))
         vslope = ((KVHS(ik,ir)/QTIM)-targetvb)/kss(ik,ir)
         vb = targetvb + s * vslope

      endif 
     
      VBQTIM = VB*QTIM
      ALPHA(IK,IR) = SQRT(CRMB/(2*TB*EMI))
      ALPHAI(IK,IR) = SQRT(CRMB/(2*TB*EMI))/QTIM
      LAMBDA2 = LAMBDA1(IZ)*nb
     
c      CC2(IK,IR) = -4.9113E11*KVHGS(IK,IR)*SQRT(CRMB*TB*TB*TB)/
c     >         (RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*KNBS(IK,IR))
c
c      CC1(IK,IR) = 1.20246E16*TB*SLOPE/
c     >              (REAL(RIZB*RIZB*RIZB*RIZB)*KNBS(IK,IR))

c      CC2(IK,IR) = -4.9113E11*vslope*SQRT(CRMB*TB*TB*TB)/
c     >         (RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*nb)

c      CC1(IK,IR) = 1.20246E16*TB*tSLOPE/
c     >              (REAL(RIZB*RIZB*RIZB*RIZB)*nb)

      CC1(IK,IR) = 1.8024009E17*TB*tslope/
     >       (Coulomb_log*REAL(RIZB*RIZB*RIZB*RIZB)*nb)
      
      CC2(IK,IR) = -7.361216E12*vslope*SQRT(CRMB*TB*TB*TB)/
     >       (Coulomb_log*RIZB*RIZB*RIZB*(1.2*RIZB + 0.8485)*nb) 

      return
      end
    
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c

      Subroutine DATA3DI

c----------------------------------------------------------------------c
c
c   This subroutine prints the arrays of background values to the .lim 
c   file in the RESULTS directory in such a manner that the average 
c   value per grid cell for each ring is written in a column with each 
c   adjacent column corresponding to the next ring--done to allow for 
c   easier down-loading of the data onto disk for tranferance to Excel 
c   to produce 3D-Contour plots.  (Called in TAU.D6A)
c
c----------------------------------------------------------------------c

      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_reiser_com
      IMPLICIT NONE

c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c     include 'reiser_com'
      
      INTEGER  IK,IR
      REAL     VHFACT,EFACT

      IK    = 0
      IR    = 0

      VHFACT = QTIM
      EFACT  = QTIM * QTIM * EMI / CRMI

c     Write nb, Tb, Te, vb, E-field values only for the SOL region
c     --rings IRSEP to IRWALL for a Jet Grid where the number of knots
c     for all the rings in the SOL is conveniently NKS(IRSEP)=78.


      WRITE (6,*) ' '
      WRITE (6,*) ' nb in the SOL: '
      WRITE (6,*) ' -------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(KNBS(IK,IR),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Te in the SOL: '
      WRITE (6,*) ' -------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(KTEBS(IK,IR),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Tbi in the SOL: '
      WRITE (6,*) ' --------------- '
      DO IK = 1, NKS(IRSEP)
         WRITE (6,10)(KTIBS(IK,IR),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Tgrad in the SOL: '
      WRITE (6,*) ' ----------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(KFIGS(IK,IR),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' vb in the SOL: '
      WRITE (6,*) ' -------------- '
      DO IK = 1, NKS(IRSEP)
         WRITE (6,10)(KVHS(IK,IR)/VHFACT,IR=IRWALL,IRSEP,-1)
      END DO

      WRITE (6,*) ' '
      WRITE (6,*) ' vbgrad in the SOL: '
      WRITE (6,*) ' ------------------ '
      DO IK = 1, NKS(IRSEP)
         WRITE (6,10)(KVHGS(IK,IR)/VHFACT,IR=IRWALL,IRSEP,-1)
      END DO

      WRITE (6,*) ' '
      WRITE (6,*) ' E-field in the SOL: '
      WRITE (6,*) ' ------------------- '
      DO IK = 1, NKS(IRSEP)
         WRITE (6,10)(KES(IK,IR)/EFACT,IR=IRWALL,IRSEP,-1)
      END DO
      WRITE (6,*) ' '

      WRITE (6,*) ' '
      WRITE (6,*) ' Cell-center distances (s) : '
      WRITE (6,*) ' --------------------------- '
      DO IK = 1, NKS(IRSEP)
         WRITE (6,10)(KSS(IK,IR),IR=IRSEP,IRWALL,1)
      END DO
      WRITE (6,*) ' '

      WRITE (6,*) ' '
      WRITE (6,*) ' Cell-areas (m^2) : '
      WRITE (6,*) ' ------------------ '
      DO IK = 1, NKS(IRSEP)
         WRITE (6,10)(KAREAS(IK,IR),IR=IRSEP,IRWALL,1)
      END DO
      WRITE (6,*) ' '



   10 FORMAT (13E13.8) 

      RETURN         
      END
 
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c   

      Subroutine DATA3DII(ZZ)

c----------------------------------------------------------------------c
c
c   This subroutine performs the same function as DATA3DI above except
c   it writes the triple-indexed arrays of the impurity density, velocity
c   and forces acting on them to the .lim file. (Called in DIV.D6A)  
c
c----------------------------------------------------------------------c

      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_reiser_com
      use mod_dynam1
      IMPLICIT NONE

c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c     include 'reiser_com'
c     include 'dynam1'
      
      INTEGER  IK,IR,IZ,ZZ

      IK    = 0
      IR    = 0

c     Write values only for the SOL region--rings IRSEP to IRWALL
c     for Jet Grid where total knots for all rings is NKS(IRSEP)=78.
c     MAIN and TRAP regions are also included.
c
c      Write (0,*)'NKS(8)',NKS(8),'IRWALL',IRWALL,'IRSEP',IRSEP
c
      DO IZ = 1,ZZ 
c
c     SOL REGION
c
      WRITE (6,*) ' '
      WRITE (6,*) ' Ftotal in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(Fcell(IK,IR,IZ),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' FF in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(Ffi(IK,IR,IZ),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' FIG in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(Fthi(IK,IR,IZ),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Fvbg in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(Fvbg(IK,IR,IZ),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' DIFF in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(DIFF(IK,IR,IZ),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Imp Velocity in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ------------------------ '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(Velavg(IK,IR,IZ),IR=IRWALL,IRSEP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Imp Density in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ----------------------- '
      DO IK = 1, NKS(IRSEP) 
         WRITE (6,10)(DDLIMS(IK,IR,IZ),IR=IRWALL,IRSEP,-1)
      END DO

      WRITE (6,*) ' '
      WRITE (6,*) ' CHI Values in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------------- '
      DO IK = 1, NKS(IRSEP) 
         if (ktibs(ik,ir).ne.0.0) then 
       WRITE (6,10)(SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))*(Velavg(IK,IR,IZ)
     >              -KVHS(IK,IR)/QTIM),IR=IRWALL,IRSEP,-1)
         endif
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' CHI Values (when Vz=0) in the SOL (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------------------------- '
      DO IK = 1, NKS(IRSEP) 
         if (ktibs(ik,ir).ne.0.0) then 
       WRITE (6,10)(-SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))*
     >               KVHS(IK,IR)/QTIM,IR=IRWALL,IRSEP,-1)
         endif
      END DO
c
c     MAIN REGION
c
      WRITE (6,*) ' '
      WRITE (6,*) ' Ftotal in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(1) 
         WRITE (6,10)(Fcell(IK,IR,IZ),IR=IRSEP-1,1,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' FF in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(1) 
         WRITE (6,10)(Ffi(IK,IR,IZ),IR=IRSEP-1,1,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' FIG in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(1) 
         WRITE (6,10)(Fthi(IK,IR,IZ),IR=IRSEP-1,1,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Fvbg in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(1) 
         WRITE (6,10)(Fvbg(IK,IR,IZ),IR=IRSEP-1,1,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' DIFF in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(1) 
         WRITE (6,10)(DIFF(IK,IR,IZ),IR=IRSEP-1,1,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Imp Velocity in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ------------------------ '
      DO IK = 1, NKS(1) 
         WRITE (6,10)(Velavg(IK,IR,IZ),IR=IRSEP-1,1,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Imp Density in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ----------------------- '
      DO IK = 1, NKS(1) 
         WRITE (6,10)(DDLIMS(IK,IR,IZ),IR=IRSEP-1,1,-1)
      END DO

      WRITE (6,*) ' '
      WRITE (6,*) ' CHI Values in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------------- '
      DO IK = 1, NKS(1) 
         if (ktibs(ik,ir).ne.0.0) then 
       WRITE (6,10)(SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))*(Velavg(IK,IR,IZ)
     >              -KVHS(IK,IR)/QTIM),IR=IRSEP-1,1,-1)
         endif
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' CHI Values (when Vz=0) in MAIN (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------------------------- '
      DO IK = 1, NKS(1) 
         if (ktibs(ik,ir).ne.0.0) then 
       WRITE (6,10)(-SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))*
     >               KVHS(IK,IR)/QTIM,IR=IRSEP-1,1,-1)
         endif
      END DO
c
c     TRAP REGION
c
      WRITE (6,*) ' '
      WRITE (6,*) ' Ftotal in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRTRAP) 
         WRITE (6,10)(Fcell(IK,IR,IZ),IR=NRS,IRTRAP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' FF in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRTRAP) 
         WRITE (6,10)(Ffi(IK,IR,IZ),IR=NRS,IRTRAP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' FIG in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRTRAP) 
         WRITE (6,10)(Fthi(IK,IR,IZ),IR=NRS,IRTRAP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Fvbg in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRTRAP) 
         WRITE (6,10)(Fvbg(IK,IR,IZ),IR=NRS,IRTRAP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' DIFF in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------- '
      DO IK = 1, NKS(IRTRAP) 
         WRITE (6,10)(DIFF(IK,IR,IZ),IR=NRS,IRTRAP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Imp Velocity in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ------------------------ '
      DO IK = 1, NKS(IRTRAP) 
         WRITE (6,10)(Velavg(IK,IR,IZ),IR=NRS,IRTRAP,-1)
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' Imp Density in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ----------------------- '
      DO IK = 1, NKS(IRTRAP) 
         WRITE (6,10)(DDLIMS(IK,IR,IZ),IR=NRS,IRTRAP,-1)
      END DO

      WRITE (6,*) ' '
      WRITE (6,*) ' CHI Values in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------------- '
      DO IK = 1, NKS(IRTRAP) 
         if (ktibs(ik,ir).ne.0.0) then 
       WRITE (6,10)(SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))*(Velavg(IK,IR,IZ)
     >              -KVHS(IK,IR)/QTIM),IR=NRS,IRTRAP,-1)
         endif
      END DO 

      WRITE (6,*) ' '
      WRITE (6,*) ' CHI Values (when Vz=0) in TRAP (IZ=',IZ,'): '
      WRITE (6,*) ' ---------------------------------- '
      DO IK = 1, NKS(IRTRAP) 
         if (ktibs(ik,ir).ne.0.0) then 
       WRITE (6,10)(-SQRT(CRMB/(2*KTIBS(IK,IR)*EMI))*
     >               KVHS(IK,IR)/QTIM,IR=NRS,IRTRAP,-1)
         endif
      END DO
     
      WRITE (6,*) ' '

      END DO
       
   10 FORMAT (13E9.3)

      RETURN         
      END
 
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
   
      Subroutine VEL_DIST(Y,vel,vmax,vmin)

c----------------------------------------------------------------------c
c
c     This subroutine demonstrates the nature of the impurity ion 
c     velocity distribution ( for replication of a Maxwellian 
c     distritution, only the FF(K11) and FPG(D11) opt should be selected 
c     with a homogeneous background plasma for injection of particles 
c     on a sinlge ring.  The velocity values are binned into 201 bins. 
c     The array storing this info is then sent to the .lim file via the 
c     subroutine VEL_DATA. (Called in DIV.D6A)
c
c----------------------------------------------------------------------c
c

      use mod_params
      use mod_cgeom
      use mod_reiser_com
      use mod_dynam1
      IMPLICIT NONE

c     include 'params'
c     include 'cgeom'
c     include 'reiser_com'
c     include 'dynam1'

      REAL J,vel,vmax,vmin
      INTEGER Y(-100:100),K
      LOGICAL UNBIN

c   MAXWELLIAN DISTRIBUTION ANALYSIS 
          
      if (vel.gt.vmax) then
	  vmax=vel
      endif

      if(vel.lt.vmin) then
          vmin=vel
      end if
         
      UNBIN = .TRUE.
      J = 0.0
      K = 0

      DO WHILE (UNBIN.AND.K.LT.101)
        IF (vel.GE.-(6.0e-6+j).AND.vel.LE.(6.0e-6+j)) THEN
          IF(vel.lt.0.0)THEN
            Y(-K) = Y(-K)+1
          ELSE
            Y(K)  = Y(K)+1
          ENDIF    
            UNBIN = .FALSE.
        ELSE
          J=J+1.2e-5
          K=K+1
        ENDIF
      END DO

      RETURN
      END

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c

      Subroutine VEL_DATA(Y,Vmax,Vmin)

c----------------------------------------------------------------------c
c
c     This subroutine sends the resulting array of velocity values  
c     calculated in VEL_DATA above to the .lim file (Called in DIV.D6A)
c     when SWITCHV = 1 and additionally to the screen when SWITCHV = 2.
c     
c----------------------------------------------------------------------c

      use mod_params
      use mod_cgeom
      use mod_reiser_com
      use mod_dynam1
      IMPLICIT NONE

c     include 'params'
c     include 'cgeom'
c     include 'reiser_com'
c     include 'dynam1'

      REAL Vmax,Vmin
      INTEGER Y(-100:100),i,velcount

      WRITE(6,*)' '
      WRITE(6,*)'Velocity Distribution:'
      WRITE(6,*)'----------------------'
     
      DO i = -100, 100
       if (i.lt.101)then
         velcount = velcount + Y(i)
       endif
       write(6,7)i,Y(i)
       IF (SWITCHV.EQ.2) write(0,7)i,Y(i)
    7  format(1x,I4,1x,g9.3)
      END DO 
     
      WRITE(6,*)' Vmax: ',Vmax,' Vmin: ',Vmin,'Velcount: ',velcount
      WRITE(6,*)' '

      IF (SWITCHV.EQ.2) THEN
        WRITE(0,*)' Vmax: ',Vmax,' Vmin: ',Vmin,'Velcount: ',velcount
      ENDIF
      
      RETURN
      END  

c----------------------------------------------------------------------c  

