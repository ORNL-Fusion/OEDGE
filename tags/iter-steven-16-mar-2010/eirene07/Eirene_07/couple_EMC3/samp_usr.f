      SUBROUTINE SAMPLE_N0_W7(L_ADD_EIR,NELEM,X,Y,Z)                
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE 
      LOGICAL :: L_ADD_EIR
      INTEGER :: NELEM
      REAL*8  :: X,Y,Z

C L_ADD_EIR: .TRUE. particle starts from a additional surface defined
C                   in eirene
C     NELEM: addational surface number
C  X,Y,Z:    Birth position

      INTEGER,PARAMETER :: N_RUN=5
      INTEGER           :: NSRACT,JRUN,IERR
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
C--------------------------------------------------------------------
      L_ADD_EIR = NS_PLACE /= 0

      NSRACT = NSRACT_DET()
      IF    ( N0S_DEF < 0 ) THEN
        X   = XSOURP(NSRACT)                                 
        Y   = YSOURP(NSRACT)                               
        Z   = ZSOURP(NSRACT)                                
                                                                        
        NELEM = NADDSF(NSRACT)                                      
      ELSEIF( N0S_DEF == 0 )THEN 
        DO JRUN=1,N_RUN
          CALL LOAD_N0(NSRACT,NELEM,X,Y,Z)
          IF(NELEM /= 0) EXIT
          NSRACT = NSRACT_DET()
        ENDDO
        IF(NELEM == 0 ) THEN 
          WRITE(6,*)'No point source'           
          CALL STOP_ALL('Stopped in SAMPLE_N0_W7')
        ENDIF 
      ELSE                       
        WRITE(6,*)'N0S_DEF > 0 to be written'
        CALL STOP_ALL('Stopped in SAMPLE_N0_W7')
      ENDIF 

      IF(NS_PLACE.EQ.0) THEN
        ITRIA  = NELEM
        NSR_N0 =-NSSIDE
      ELSE
        ITRIA  = 0
        NSR_N0 = 0
      ENDIF 
      RETURN

      CONTAINS
C---------------------------------------------------------------------
      FUNCTION NSRACT_DET()
      IMPLICIT NONE

      INTEGER :: NSRACT_DET
      INTEGER :: JBEGIN,JEND,JTRY 
      REAL*8  :: WTMONT, RANF

      WTMONT = RANF()*GEWICH(NSSTOR)
      IF(WTMONT.LE.GEWICH(1)) THEN
         NSRACT_DET = 1
      ELSE 
        JBEGIN = 1
        JEND   = NSSTOR

        SEARCH_FOR_NSRACT : DO
           JTRY   = (JBEGIN+JEND)/2                                          
           IF(  WTMONT.GT.GEWICH(JBEGIN) .AND.
     +          WTMONT.LE.GEWICH(JTRY)         )THEN
              IF(JTRY.EQ.JBEGIN+1) THEN
                NSRACT_DET = JTRY
                EXIT SEARCH_FOR_NSRACT
              ELSE
                JEND = JTRY
              ENDIF
           ELSE
              IF(JTRY.EQ.JEND - 1) THEN
                NSRACT_DET = JEND
                EXIT SEARCH_FOR_NSRACT
              ELSE
                JBEGIN = JTRY
              ENDIF
           ENDIF
        ENDDO SEARCH_FOR_NSRACT  
      ENDIF
      END FUNCTION NSRACT_DET
C---------------------------------------------------------------------
      SUBROUTINE LOAD_N0(NSR,NELEM,X,Y,Z)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NSR
      INTEGER,INTENT(OUT):: NELEM
      REAL*8, INTENT(OUT):: X,Y,Z

      INTEGER :: IRUN
      INTEGER :: NZ0,JR0,JP0,JT0,ID0
      REAL*8  :: RJ0,PJ0, RANF
      
      NZ0 = IZO_N0(NSR)
      JR0 = IR0_N0(NSR)
      JP0 = IP0_N0(NSR)
      JT0 = IT0_N0(NSR)
      ID0 = ISD_N0(NSR)

      DO IRUN=1,N_RUN
       RJ0    = 2.*RANF() - 1.
       PJ0    = 2.*RANF() - 1.
       CALL LOAD_PARTICLE(NZ0,JR0,JP0,JT0,ID0,RJ0,PJ0,NELEM,X,Y,Z)
       IF( NELEM /=0) EXIT
      ENDDO
      END SUBROUTINE LOAD_N0
      END SUBROUTINE SAMPLE_N0_W7
C---------------------------------------------------------------------
      SUBROUTINE LOAD_PARTICLE(NZ0,JR0,JP0,JT0,ID0,RJ0,PJ0,NELEM,X,Y,Z)
      USE GEOMETRY_PL       
      USE SURFACE_PL        
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE
      INTEGER,INTENT(IN)  :: NZ0,JR0,JP0,JT0,ID0
      REAL*8 ,INTENT(IN)  :: RJ0,PJ0
      INTEGER,INTENT(OUT) :: NELEM
      REAL*8 ,INTENT(OUT) :: X,Y,Z

      INTEGER,PARAMETER   :: NT_S=5
      INTEGER :: IDT,ID,ISU,IST,K
      REAL*8  :: RJ,PJ,TJ,R1,R2,Z1,Z2,VL

      NELEM = 0            
      ITRIA = 0
      LOOP_POS_NEG_SEARCH : DO ID =-1,1,2
        IDT    =-ID0*ID             
        NZ0_N0 = NZ0
        JR0_N0 = JR0
        JP0_N0 = JP0
        RJ     = RJ0
        PJ     = PJ0
        TJ     = JT0
        JT0_N0 = TJ
        
        LOOP_K_SEARCH : DO K=0,NT_S
          IF( (JT0_N0 == 0                .AND. IDT < 0) .OR.   
     .        (JT0_N0 == ZON_TORO(NZ0_N0) .AND. IDT > 0) ) THEN
              ISU = JR0_N0+(JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*             ZON_RADI(NZ0_N0) + NTS_OFF(NZ0_N0)
              IST= IDSURT(ISU) 
          
              CALL T_SF_JUMP
     .        (IST,IDT,NZ0_N0,JR0_N0,JP0_N0,JT0_N0,RJ,PJ,TJ)

              CALL RZ_REAL_COORDINATES
     .        (NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,TJ,R1,Z1)
          ELSEIF(K == 0 ) THEN 
              CALL RZ_REAL_COORDINATES
     .        (NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,TJ,R1,Z1)
          ENDIF 

          JT0_N0 = TJ
          XP0 = R1*COSPHI(JT0_N0+PHI_PL_OS(NZ0_N0))
          YP0 = R1*SINPHI(JT0_N0+PHI_PL_OS(NZ0_N0))
          ZP0 = Z1

          TJ     = TJ + IDT
          JT0_N0 = TJ
              CALL RZ_REAL_COORDINATES
     .       (NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,TJ,R2,Z2)

          X   = R2*COSPHI(JT0_N0+PHI_PL_OS(NZ0_N0))
          Y   = R2*SINPHI(JT0_N0+PHI_PL_OS(NZ0_N0))
          Z   = Z2

          VL  = SQRT( (X-XP0)**2 + (Y-YP0)**2 + (Z-ZP0)**2 )
          VX  = (X-XP0)/VL
          VY  = (Y-YP0)/VL
          VZ  = (Z-ZP0)/VL
 
          CALL TLMAX_TRAVEL
          IF(TLMAX < VL) THEN 
            X = XP0 + TLMAX*VX
            Y = YP0 + TLMAX*VY
            Z = ZP0 + TLMAX*VZ
            IF(IDT == 1) JT0_N0 = JT0_N0-1
            NELEM = ITRIA
            EXIT LOOP_POS_NEG_SEARCH
          ENDIF 
          R1 = R2
          Z1 = Z2
        ENDDO LOOP_K_SEARCH
      ENDDO LOOP_POS_NEG_SEARCH

      RETURN 
      END SUBROUTINE LOAD_PARTICLE
