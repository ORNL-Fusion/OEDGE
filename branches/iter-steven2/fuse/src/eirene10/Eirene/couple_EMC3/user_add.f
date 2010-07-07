      SUBROUTINE VOLUSR(NRGEN,VOL)
      USE PHYSICAL_CELL 
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE              
      INTEGER,INTENT(IN) :: NRGEN 
      REAL*8, DIMENSION(*),INTENT(OUT) :: VOL

      IF( NRGEN < NCELL_N2+1) THEN
        WRITE(6,*)'NRGEN=',NRGEN,' < NCELL_N2+1=',NCELL_N2+1 
        CALL STOP_ALL('Stopped in VOLUSR')
      ENDIF  
     
      VOL(1:NC_PL   )        = VOLCEL(1:NC_PL)   
      VOL(NCELL_N1:NCELL_N2) = VOLCEL_N0(NCELL_N1:NCELL_N2)
      VOL(NCELL_N2+1:NRGEN)  = 1.

C Define additional tallies for coupling
      CALL INICOP
      RETURN
      END SUBROUTINE VOLUSR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SUBROUTINE PROUSR (PRO,INDEX,PRO0,PROS,P,Q,E,SEP,PROVAC,N)
      USE PHYSICAL_CELL
      USE PLASMA_PARAM     
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE

      REAL*8, DIMENSION(*) :: PRO
      REAL*8  :: PRO0,PROS,P,Q,E,SEP,PROVAC
      INTEGER :: INDEX,N,ITYPE,IC
      REAL*8, PARAMETER :: TE_MIN=0.1,TI_MIN=0.1,NE_MIN=1.E7
      
10    if (index.eq.0) then
          PRO(1:NC_PL) = MAX(TEMP0(1:NC_PL,0),TE_MIN) 
          DO ITYPE=1,NCTYPE_DEF                 
          DO IC=NCELL_AB(ITYPE-1)+1,NCELL_AB(ITYPE)
          PRO(IC) = MAX(PLASMA_ADD_CELL(2,ITYPE),TE_MIN)
          ENDDO
          ENDDO 
      elseif (index.eq.1) then
          PRO(1:NC_PL) = MAX(TEMP0(1:NC_PL,1),TI_MIN)
          DO ITYPE=1,NCTYPE_DEF                 
          DO IC=NCELL_AB(ITYPE-1)+1,NCELL_AB(ITYPE)
          PRO(IC) = MAX(PLASMA_ADD_CELL(3,ITYPE),TI_MIN) 
          ENDDO
          ENDDO 
      elseif (index.eq.2) then
          PRO(1:NC_PL) = MAX(DENS0(1:NC_PL,0),NE_MIN) 
          DO ITYPE=1,NCTYPE_DEF                 
          DO IC=NCELL_AB(ITYPE-1)+1,NCELL_AB(ITYPE)
          PRO(IC) = MAX(PLASMA_ADD_CELL(1,ITYPE),NE_MIN) 
          ENDDO
          ENDDO 
      endif
      RETURN
      END SUBROUTINE PROUSR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     NAME: VECUSR                                                    C
C FUNCTION: DEFINE LOCAL VECTORS                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SUBROUTINE VECUSR(MFLAG,VECTX,VECTY,VECTZ,IPLS)
      USE GEOMETRY_PL
      USE PHYSICAL_CELL
      USE PLASMA_PARAM     
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE
C   INPUT:
C   MFLAG: =1, B-LINE
C          =2, ION DRIFT VELOCITY
C          =3, VECTX=magnitude of ION DRIFT VELOCITY
C          =4, VECTX=VX*BX+VY*BY+VZ*BZ                         
C          >4, TO BE WRIITEN             
C  OUTPUT:
C  VECTX,-Y,Z: VECTOR REQUIRED                           
      INTEGER :: MFLAG,IPLS,IG,IC
      REAL*8  :: VECTX,VECTY,VECTZ,V
 
C 1. B-LINE 
      IG    = JR0_N0+(JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*       ZON_RADI(NZ0_N0) + MESH_P_OS(NZ0_N0)
      IF(MFLAG.EQ.1) THEN
        VECTX = B_VECT(1,IG)
        VECTY = B_VECT(2,IG)   
        VECTZ = B_VECT(3,IG)
C 2. ION DRIFT VELOCITY
      ELSEIF(MFLAG.EQ.2) THEN
        IC = IDCELL(IG)
 
        IF( IC <= NC_PL) THEN
          V  = ZMACH0(IC)*SOUND_S(IC)
        ELSE
          V  = 0.
        ENDIF 
 
        VECTX = B_VECT(1,IG)*V
        VECTY = B_VECT(2,IG)*V 
        VECTZ = B_VECT(3,IG)*V 
      ELSEIF(MFLAG.EQ.3) THEN
        IC = IDCELL(IG)
                       
        VECTX =ZMACH0(IC)*SOUND_S(IC)
C 3. VX*BX+VY*BY+VZ*BZ    
      ELSEIF(MFLAG.EQ.4) THEN
        IC = IDCELL(IG)
                               
        VECTX = B_VECT(1,IG)
        VECTY = B_VECT(2,IG)
        VECTZ = B_VECT(3,IG)
        VECTX = VX*VECTX+VY*VECTY+VZ*VECTZ
      ELSE
C 4. TO BE WRITTEN        
        WRITE(6,*)'VECTOR IS CALLED WITH AN UNKOWN MFLAG=',MFLAG
        CALL STOP_ALL('STOPPED IN VECUSR')
      ENDIF
      RETURN
      END SUBROUTINE VECUSR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      FUNCTION VDION(I)                
      USE PHYSICAL_CELL
      USE PLASMA_PARAM     

      IMPLICIT NONE                            
      REAL*8  :: VDION 
      INTEGER :: I
      
      VDION = 0.
      IF(I<=NC_PL) VDION = ZMACH0(I)*SOUND_S(I)
      RETURN
      END FUNCTION VDION
C------------------------------------------------------------------------
C This subroutine is called by PE=0
      SUBROUTINE OUT_FOR_EMC3(N0WR)
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE                  
 
      INTEGER,INTENT(IN) :: N0WR 
      INTEGER            :: IPLS,ITYP,ICA,ICB,I

      REAL*8,DIMENSION(:),ALLOCATABLE :: TEMP
       
      ALLOCATE ( TEMP(NCELL_N2) )
CCCCCC

      OUT_N0 = 0.

      IPLS = 1
C 1. Ionization source 
      CALL OUT_PROF(1,IPLS,1,NCELL_N2,TEMP) 
      OUT_N0(1:NCELL_N2,1) = TEMP(1:NCELL_N2)

C 2. Energy source due to neutral
C 2.1 for e  
      CALL OUT_PROF(2,IPLS,1,NCELL_N2,TEMP)   
      OUT_N0(1:NCELL_N2,2) = TEMP(1:NCELL_N2)

c 2.2 for i  
      CALL OUT_PROF(3,IPLS,1,NCELL_N2,TEMP)   
      OUT_N0(1:NCELL_N2,3) = TEMP(1:NCELL_N2)

C 3. Momentum loss due to neutral 
C CX rate 
C 1/S
c     CALL OUT_PROF(30,IPLS,1,NCELL_N2,TEMP)
C ATOM 
C G*CM/S*CM**-3
      CALL OUT_PROF(31,IPLS,1,NCELL_N2,TEMP)
      OUT_N0(1:NCELL_N2,4) = OUT_N0(1:NCELL_N2,4) + TEMP(1:NCELL_N2)
C MOLE 
C G*CM/S*CM**-3
      CALL OUT_PROF(32,IPLS,1,NCELL_N2,TEMP)
      OUT_N0(1:NCELL_N2,4) = OUT_N0(1:NCELL_N2,4) + TEMP(1:NCELL_N2)
C TION 
C G*CM/S*CM**-3
      CALL OUT_PROF(33,IPLS,1,NCELL_N2,TEMP)
      OUT_N0(1:NCELL_N2,4) = OUT_N0(1:NCELL_N2,4) + TEMP(1:NCELL_N2)

C 4. Atom density
      CALL OUT_PROF(-1,IPLS,1,NCELL_N2,TEMP)  
      OUT_N0(1:NCELL_N2,5) = TEMP(1:NCELL_N2)

C 5. Molecule density
      CALL OUT_PROF(-2,IPLS,1,NCELL_N2,TEMP)  
      OUT_N0(1:NCELL_N2,6) = TEMP(1:NCELL_N2)
 
      DEALLOCATE ( TEMP )
      RETURN 
      END
