C Name    : NEUTRAL_TRANSPORT
C Function:
C Contains:
C   Subroutines  : READ_SF_N0         UPDATA_SF_N0
C                  READ_CELL_NR_N0    CELL_DEF0           CELL_DEF1
C                  CELL_DEF2          UPDATA_CELL_NR_N0
C                  READ_SOURCE_N0     UPDATA_SOURCE_N0
C                  READ_ADD_SF_N0     UPDATA_ADD_SF_N0    LIM_ADD
C                  OPEN_NEUTRAL_TRANSPORT   COSE_NEUTRAL_TRANSPORT
C   Functions    : None
C Use     :
C   Subroutines  : None
C   Functions    : None
C   Common blocks: None
C   Modules      : None
C
       MODULE NEUTRAL_TRANSPORT
       IMPLICIT NONE
C---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
C 1. Geometry
C B_VECT: B-line unit vector
C SINPHI: sinus  (phi_plane)
C COSPHI: cosinus(phi_plane)
      REAL*8, ALLOCATABLE,DIMENSION(:,:), SAVE :: 
     R        B_VECT
      REAL*8, ALLOCATABLE,DIMENSION(  :), SAVE :: 
     R        SINPHI,      COSPHI
C---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
C 2. Neutral cell
C NCELL_N1 = NC_PL + 1
C NCELL_N2 = Maximal cell number for neutrals
C            NCELL_N2-NCELL_N1+1: Total number for neutrals
C NCTYPE_DEF: Type of defination
C NCELL_AB  : Cell number range
C VOLCEL_N0 : Volume of a cell
C PLASMA_ADD_CELL: Plasma data for n0-cells
      INTEGER,                         SAVE ::
     I        NCELL_N1,    NCELL_N2,    NCTYPE_DEF
      INTEGER,DIMENSION(:  ),ALLOCATABLE,SAVE ::
     I        NCELL_AB
      REAL*8, DIMENSION(:  ),ALLOCATABLE,SAVE ::
     R        VOLCEL_N0
      REAL*8, DIMENSION(:,:),ALLOCATABLE,SAVE ::
     R        PLASMA_ADD_CELL

      REAL*8, DIMENSION(:,:),ALLOCATABLE,SAVE ::
     R        OUT_N0                   
C---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
C 3. NON-DEFAULT, NON-TRANSPARENT (user facing)
      INTEGER,                           SAVE ::
     I        NONTRT_N0,  NONTPT_N0,  NONTTT_N0 
      INTEGER,DIMENSION(:  ),ALLOCATABLE,SAVE ::
     I  IZONR_N0,   IR0NR_N0      
     I, IP1NR_N0,   IP2NR_N0                     
     I, IT1NR_N0,   IT2NR_N0                     
     I,INDINR_N0    
c    poloidal                                  
     I, IZONP_N0,   IP0NP_N0      
     I, IR1NP_N0,   IR2NP_N0                     
     I, IT1NP_N0,   IT2NP_N0                     
     I,INDINP_N0 
c    toroidal                                  
     I, IZONT_N0,   IT0NT_N0      
     I, IR1NT_N0,   IR2NT_N0                     
     I, IP1NT_N0,   IP2NT_N0                     
     I,INDINT_N0
C---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
C 4. Additional surface         
      INTEGER,                           SAVE ::
     I        NPLATE,       NPSORT,      NTRIANG
      INTEGER,DIMENSION(:,:),ALLOCATABLE,SAVE ::
     I        NSTUCK
      INTEGER,DIMENSION(:  ),ALLOCATABLE,SAVE ::
     I        ISR_TYPE       
     I,       NP_NUM                                 
     I,       NP_GRP,       NPS_GRP    
     I,       NPN_GRP       

      REAL*8, DIMENSION(:,:),ALLOCATABLE,SAVE ::
     R        X_TRIA,       Y_TRIA,       Z_TRIA
      REAL*8, DIMENSION(:  ),ALLOCATABLE,SAVE ::
     R        A_PLAT,       A_TRIA
     R,       V_TRIA_X,     V_TRIA_Y,     V_TRIA_Z
     R,       X_SORT,       Y_SORT,       Z_SORT
     R,       D_MAXI,       D_MINI
     R,       X_SMIN,       X_SMAX
     R,       Y_SMIN,       Y_SMAX
     R,       Z_SMIN,       Z_SMAX
     R,       V_SORT_X,     V_SORT_Y,     V_SORT_Z
C---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
C 5. source distribution        
      INTEGER,                           SAVE ::
     I        NSSTOR,       N0S_DEF,      NS_PLACE,       NSSIDE
      INTEGER,DIMENSION(:  ),ALLOCATABLE,SAVE ::
     I        NADDSF,       IJK_N0,       IZO_N0,         ISD_N0        
     I,       IR0_N0,       IP0_N0,       IT0_N0         
      REAL*8, DIMENSION(:  ),ALLOCATABLE,SAVE ::
     R        GEWICH,       XSOURP,       YSOURP,         ZSOURP
C---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
C 6. transport geometry         
      INTEGER,                           SAVE :: 
     I        IDRJ=-1,      IDPJ=-1,     IDTJ=-1 
      REAL*8, DIMENSION(8  ),            SAVE ::
     R        XMY,          YMY,          ZMY
      INTEGER,DIMENSION(12),             SAVE :: 
     I         NEWSRF=(/ 5, 6, 8, 7, 1, 2, 4, 3,11,12, 9,10/)
     I,        I_JUMP=(/-1,-1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0/)
     I,        J_JUMP=(/ 0, 0,-1,-1, 0, 0, 1, 1, 0, 0, 0, 0/)
     I,        K_JUMP=(/ 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 1, 1/)
      INTEGER,DIMENSION(-1:1),           SAVE ::
     I         N_JUMP=(/0,0,1/) 
     I,        M_JUMP=(/1,0,0/)
     I,         INSRF=(/0,0,1/) 
      INTEGER,DIMENSION(12,3),           SAVE :: 
     I        IPUNKT=RESHAPE( 
     I           (/1,1,1,1,4,4,2,2,1,2,5,8,
     I             6,5,8,4,3,7,7,6,2,3,8,7,
     I             2,6,5,8,7,8,3,7,4,4,6,6/),(/12,3/) )
C---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
C 7. Neutral particle           
      INTEGER,DIMENSION(:  ),ALLOCATABLE      ::
     I        IEBENE,        IEB_PL,     IEB_YA,       MTRIA
     I,       NYBACK,        NYSTOR
      REAL*8 ,DIMENSION(:,:),ALLOCATABLE      ::
     R        TCR_PL
      REAL*8 ,DIMENSION(:  ),ALLOCATABLE      ::
     R        DISTAN,        SINFI,      TCROSS
      REAL*8 ,DIMENSION(3)                    ::
     R        x3eck,         y3eck,      z3eck
     R,       vpr1,          vpr2,       vpr3


      INTEGER,                           SAVE ::
     I        NSR_N0
     I,       ISRF=0,        MSURF,      IPLATE,       ITRIA
     I,       NZ0_N0,        JR0_N0,     JP0_N0,       JT0_N0
     I,       IC0_N0
      REAL*8,                            SAVE ::
     R        XP0,           YP0,        ZP0
     R,       VX,            VY,         VZ
     R,       XFIN,          YFIN,       ZFIN
     R,       TFIN,          TLMAX 

      CONTAINS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE OPEN_NEUTRAL_TRANSPORT()
      USE GEOMETRY_PL
      USE SURFACE_PL
      INTEGER :: I,J,K,IZ,IG,I1,I2,I3,I4,IS,ISURF
      REAL*8  :: R1,X1,Y1,Z1,R2,X2,Y2,Z2,DD

      ALLOCATE ( B_VECT(3,0:MESH_P_OS(NZONET)-1)
     .,          SINPHI(  0:PHI_PL_OS(NZONET)-1)
     .,          COSPHI(  0:PHI_PL_OS(NZONET)-1)
     .,          DISTAN(    NTRIANG            )
     .,          SINFI (    NTRIANG            )
     .,          TCROSS(    NTRIANG            )
     .,          TCR_PL(2,  NTRIANG            )
     .,          IEBENE(    NTRIANG            )
     .,          IEB_PL(    NTRIANG            )
     .,          IEB_YA(    NTRIANG            )
     .,          MTRIA (    NTRIANG            )
     .,          NYBACK(    NTRIANG            )
     .,          NYSTOR(    NTRIANG            )
     .,          OUT_N0(    NCELL_N2,   6      )
     .         )

      DO 10 IZ=0,NZONET-1
      DO 10 K =0,ZON_TORO(IZ)
        IG = K + PHI_PL_OS(IZ)
        SINPHI(IG) = SIN(PHI_PLANE(IG))
        COSPHI(IG) = COS(PHI_PLANE(IG))
10    CONTINUE

      DO 20 IZ=0,NZONET-1
      DO 20 K =0,ZON_TORO(IZ) - 1
      DO 20 J =0,ZON_POLO(IZ) - 1
      DO 20 I =0,ZON_RADI(IZ) - 1

        IG= I + (J+K*ZON_POLO(IZ))*ZON_RADI(IZ)+MESH_P_OS(IZ)

        I1= I + (J+K*SRF_POLO(IZ))*SRF_RADI(IZ)+GRID_P_OS(IZ)
        I2= I1+ SRF_RADI(IZ)
        I3= I2+ 1
        I4= I1+ 1

        R1 = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
        Z1 = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))
        X1 = R1*COSPHI(K+PHI_PL_OS(IZ))
        Y1 = R1*SINPHI(K+PHI_PL_OS(IZ))

        I1= I1+ SRF_POLO(IZ)*SRF_RADI(IZ)
        I2= I1+ SRF_RADI(IZ)
        I3= I2+ 1
        I4= I1+ 1

        R2 = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
        Z2 = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))
        X2 = R2*COSPHI(K+1+PHI_PL_OS(IZ))
        Y2 = R2*SINPHI(K+1+PHI_PL_OS(IZ))

        DD = 1./(SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2))

        B_VECT(1,IG) = (X2-X1)*DD
        B_VECT(2,IG) = (Y2-Y1)*DD
        B_VECT(3,IG) = (Z2-Z1)*DD
20    CONTINUE

C set up non-default and non-transparent surfaces for neutrals
      IDSURR = 0
      IDSURP = 0
      IDSURT = 0
C 1. non-default (same as for plasma)
      CALL UPDATA_NON_DEFAULT_SF()
C 2. non-transparent 
C user-defined 
      DO IS=1,NONTRT_N0
       IF(INDINR_N0(IS) >= 0) THEN 
         WRITE(6,*)'A radial non-transp. surface identified by IND>=0'
         CALL STOP_ALL('Stopped in OPEN_NEUTRAL') 
       ENDIF 
       IZ  = IZONR_N0(IS)
       I   = IR0NR_N0(IS)
       DO J= IP1NR_N0(IS),IP2NR_N0(IS)
       DO K= IT1NR_N0(IS),IT2NR_N0(IS)
       ISURF = I + (J+K*ZON_POLO(IZ))*SRF_RADI(IZ) + NRS_OFF(IZ)
       IDSURR(ISURF) = INDINR_N0(IS)
       ENDDO
       ENDDO
      ENDDO

      DO IS=1,NONTPT_N0
       IF(INDINP_N0(IS) >= 0) THEN 
         WRITE(6,*)'A polo. non-transp. surface identified by IND>=0'
         CALL STOP_ALL('Stopped in OPEN_NEUTRAL') 
       ENDIF 
       IZ  = IZONP_N0(IS)
       J   = IP0NP_N0(IS)
       DO I= IR1NP_N0(IS),IR2NP_N0(IS)
       DO K= IT1NP_N0(IS),IT2NP_N0(IS)
       ISURF = I + (J+K*SRF_POLO(IZ))*ZON_RADI(IZ) + NPS_OFF(IZ)
       IDSURP(ISURF) = INDINP_N0(IS)
       ENDDO
       ENDDO
      ENDDO

      DO IS=1,NONTTT_N0
       IF(INDINT_N0(IS) >= 0) THEN 
         WRITE(6,*)'A toro. non-transp. surface identified by IND>=0'
         CALL STOP_ALL('Stopped in OPEN_NEUTRAL') 
       ENDIF 
       IZ  = IZONT_N0(IS)
       K   = IT0NT_N0(IS)
       DO I= IR1NT_N0(IS),IR2NT_N0(IS)
       DO J= IP1NT_N0(IS),IP2NT_N0(IS)
       ISURF = I + (J+K*ZON_POLO(IZ))*ZON_RADI(IZ) + NTS_OFF(IZ)
       IDSURT(ISURF) = INDINT_N0(IS)
       ENDDO
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE OPEN_NEUTRAL_TRANSPORT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CLOSE_NEUTRAL_TRANSPORT()
      DEALLOCATE ( B_VECT
     .,            SINPHI
     .,            COSPHI
     .,            DISTAN
     .,            SINFI 
     .,            TCROSS
     .,            TCR_PL
     .,            IEBENE
     .,            IEB_PL
     .,            IEB_YA
     .,            MTRIA 
     .,            NYBACK
     .,            NYSTOR
     .,            OUT_N0
     .           )
      END SUBROUTINE CLOSE_NEUTRAL_TRANSPORT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READ_SF_N0(IUNIT)
      USE GEOMETRY_PL

      INTEGER,INTENT(IN) :: IUNIT
      INTEGER            :: I, IZ
      CHARACTER*72 :: READFI
      LOGICAL      :: OUT_OF_RANGE
C N_radial
      CALL SCRAPE(IUNIT,READFI,*99)
      READ(READFI,*) NONTRT_N0

      IF(NONTRT_N0.LT.0) NONTRT_N0 = 0
   
      IF( NONTRT_N0 > 0) ALLOCATE (
     I    IR0NR_N0(NONTRT_N0),IZONR_N0(NONTRT_N0),INDINR_N0(NONTRT_N0)
     I,   IP1NR_N0(NONTRT_N0),IP2NR_N0(NONTRT_N0)
     I,   IT1NR_N0(NONTRT_N0),IT2NR_N0(NONTRT_N0)
     &                            )

      WRITE(6,'(A10,I6)')'   radial:',NONTRT_N0
      DO I=1,NONTRT_N0
        CALL SCRAPE(IUNIT,READFI,*99)
        READ(READFI,*) IR0NR_N0(I),IZONR_N0(I),INDINR_N0(I)

        CALL SCRAPE(IUNIT,READFI,*99)
        READ(READFI,*) IP1NR_N0(I),IP2NR_N0(I),IT1NR_N0(I),IT2NR_N0(I)

        IZ = IZONR_N0(I)
        IF( OUT_OF_RANGE(IZ         ,0,NZONET-1         ) .OR.
     .      OUT_OF_RANGE(IR0NR_N0(I),0,SRF_RADI(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IP2NR_N0(I),0,ZON_POLO(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IP1NR_N0(I),0,IP2NR_N0(I)      ) .OR.
     .      OUT_OF_RANGE(IT2NR_N0(I),0,ZON_TORO(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IT1NR_N0(I),0,IT2NR_N0(I)      )      ) THEN
            WRITE(6,*)'RADIAL NON DEFAULT SF. OUT OF RANGE'
            CALL STOP_ALL('Stopped in READSF_N0:NEUTRAL_TRANSPORT')
         ENDIF
      ENDDO   
C N_poloidal
      CALL SCRAPE(IUNIT,READFI,*99)
      READ(READFI,*) NONTPT_N0

      IF(NONTPT_N0.LT.0) NONTPT_N0 = 0
   
      IF( NONTPT_N0 > 0) ALLOCATE (
     I    IP0NP_N0(NONTPT_N0),IZONP_N0(NONTPT_N0),INDINP_N0(NONTPT_N0)
     I,   IR1NP_N0(NONTPT_N0),IR2NP_N0(NONTPT_N0)
     I,   IT1NP_N0(NONTPT_N0),IT2NP_N0(NONTPT_N0)
     &                            )

      WRITE(6,'(A10,I6)')' poloidal:',NONTPT_N0
      DO I=1,NONTPT_N0
        CALL SCRAPE(IUNIT,READFI,*99)
        READ(READFI,*) IP0NP_N0(I),IZONP_N0(I),INDINP_N0(I)

        CALL SCRAPE(IUNIT,READFI,*99)
        READ(READFI,*) IR1NP_N0(I),IR2NP_N0(I),IT1NP_N0(I),IT2NP_N0(I)

        IZ = IZONP_N0(I)
        IF( OUT_OF_RANGE(IZ         ,0,NZONET-1         ) .OR.
     .      OUT_OF_RANGE(IP0NP_N0(I),0,SRF_POLO(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IR2NP_N0(I),0,ZON_RADI(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IR1NP_N0(I),0,IR2NP_N0(I)      ) .OR.
     .      OUT_OF_RANGE(IT2NP_N0(I),0,ZON_TORO(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IT1NP_N0(I),0,IT2NP_N0(I)      )      ) THEN
            WRITE(6,*)'TOROIDAL NON DEFAULT SF. OUT OF RANGE'
            CALL STOP_ALL('Stopped in READSF_N0:NEUTRAL_TRANSPORT')
         ENDIF
      ENDDO   
C N_radial
      CALL SCRAPE(IUNIT,READFI,*99)
      READ(READFI,*) NONTTT_N0

      IF(NONTTT_N0.LT.0) NONTTT_N0 = 0
   
      IF( NONTTT_N0 > 0) ALLOCATE (
     I    IT0NT_N0(NONTTT_N0),IZONT_N0(NONTTT_N0),INDINT_N0(NONTTT_N0)
     I,   IR1NT_N0(NONTTT_N0),IR2NT_N0(NONTTT_N0)
     I,   IP1NT_N0(NONTTT_N0),IP2NT_N0(NONTTT_N0)
     &                            )

      WRITE(6,'(A10,I6)')' toroidal:',NONTTT_N0
      DO I=1,NONTTT_N0
        CALL SCRAPE(IUNIT,READFI,*99)
        READ(READFI,*) IT0NT_N0(I),IZONT_N0(I),INDINT_N0(I)

        CALL SCRAPE(IUNIT,READFI,*99)
        READ(READFI,*) IR1NT_N0(I),IR2NT_N0(I),IP1NT_N0(I),IP2NT_N0(I)

        IZ = IZONT_N0(I)
        IF( OUT_OF_RANGE(IZ         ,0,NZONET-1         ) .OR.
     .      OUT_OF_RANGE(IT0NT_N0(I),0,SRF_TORO(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IR2NT_N0(I),0,ZON_RADI(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IR1NT_N0(I),0,IR2NT_N0(I)      ) .OR.
     .      OUT_OF_RANGE(IP2NT_N0(I),0,ZON_POLO(IZ)-1   ) .OR.
     .      OUT_OF_RANGE(IP1NT_N0(I),0,IP2NT_N0(I)      )      ) THEN
            WRITE(6,*)'TOROIDAL NON DEFAULT SF. OUT OF RANGE'
            CALL STOP_ALL('Stopped in READSF_N0:NEUTRAL_TRANSPORT')
         ENDIF
      ENDDO   
      RETURN
99    WRITE(6,*)' End_of_file'
      CALL STOP_ALL('Stopped in READSF_N0:NEUTRAL_TRANSPORT')
      END SUBROUTINE READ_SF_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UPDATA_SF_N0()
      END SUBROUTINE UPDATA_SF_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READ_CELL_NR_N0(IUNIT)
      USE GEOMETRY_PL
      USE PHYSICAL_CELL

      INTEGER,INTENT(IN) :: IUNIT
      CHARACTER*72 :: READFI,FLNAME
      REAL*8, PARAMETER :: TE_DEF=0.1,TI_DEF=0.1
     R,                   NE_DEF=1.E7,V_DEF=0.
      INTEGER :: ICA,ICB,IUM_D_C,NZ,NR1,NR2,NDR,IU
     I,          NP1,NP2,NDP,NT1,NT2,NDT,IND_CELL
     I,          IZ,IG,I,J,K,WR_CELL_G
C - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO 10 IZ=0,NZONET-1
      DO 10 K =0,ZON_TORO(IZ)-1
      DO 10 J =0,ZON_POLO(IZ)-1
      DO 10 I =0,ZON_RADI(IZ)-1
          IG =I+(J+K*ZON_POLO(IZ))*ZON_RADI(IZ)+MESH_P_OS(IZ)
          IF(IDCELL(IG)>0 .AND. IDCELL(IG)<=NC_PL) GOTO 10
          IDCELL(IG) = 0
10    CONTINUE

      NCELL_N1 = NC_PL + 1
      NCELL_N2 = NC_PL

      CALL SCRAPE(IUNIT,READFI,*99)
      READ(READFI,*) NCTYPE_DEF,WR_CELL_G

      IF(NCTYPE_DEF.LT.0) THEN
        CALL WRMESS('NCTYPE_DEF.LT.0, CELL CANNOT BE DEFINED')
        CALL STOP_ALL('Stopped in READ_CELL_NR_N0')
      ENDIF 

      ALLOCATE ( NCELL_AB       (0:NCTYPE_DEF)
     R,          PLASMA_ADD_CELL(4,NCTYPE_DEF) )

      NCELL_AB(0) = NC_PL

      CALL WRMESS(' Total types of additional cells for neutrals:')
      WRITE(6,*) '  NCTYPE_DEF=',NCTYPE_DEF
      CALL WRMESS(
     &'IND_CELL  CELL_1  CELL_2       Ne         Te        Ti       V')

      DO 20 IUM_D_C=1,NCTYPE_DEF
      CALL SCRAPE(IUNIT,READFI,*99)
      READ(READFI,*) IND_CELL
C
      SELECT CASE(IND_CELL) 
C CASE0 IND_CELL=0,geometric cell = physical cell
      CASE(0) 
         ICA = NCELL_N2 + 1
         CALL CELL_DEF0(ICA,ICB) 
         NCELL_N2 = ICB

C CASE1 IND_CELL=1,not for the coupled case
      CASE(1) 
         CALL WRMESS('IND_CELL.EQ.1 has already been defined')
         CALL STOP_ALL('Stopped in READ_CELL_NR_N0')

      CASE(2) 
C CASE2 IND_CELL=2,restrict a range and give the steps
         CALL SCRAPE(IUNIT,READFI,*99)
         READ(READFI,*) NZ,NR1,NR2,NDR,NP1,NP2,NDP,NT1,NT2,NDT

         ICA = NCELL_N2 + 1
         CALL CELL_DEF2(NZ,NR1,NR2,NDR,NP1,NP2,NDP,NT1,NT2,NDT
     +,                 ICA,ICB)
         NCELL_N2 = ICB

      CASE(3) 
C CASE3 IND_CELL=3,read cells from file FLNAME
         CALL SCRAPE(IUNIT,FLNAME,*99)

         IU = 99
         OPEN(IU,FILE=FLNAME)
           ICA = NCELL_N2 + 1
           CALL CELL_DEF1(IU,ICA,ICB)
           NCELL_N2 = ICB
         CLOSE(IU)
      CASE DEFAULT
         CALL WRMESS('3<IND_CELL<0,to be written')
         CALL STOP_ALL('Stopped in READ_CELL_NR_N0')
      END SELECT

      NCELL_AB(IUM_D_C) = ICB

      CALL SCRAPE(IUNIT,READFI,*99)
      READ(READFI,*)(PLASMA_ADD_CELL(I,IUM_D_C),I=1,4)

      WRITE(6,'(3I8,4e12.4)')
     & IND_CELL,NCELL_AB(IUM_D_C-1)+1,NCELL_AB(IUM_D_C)
     &,(PLASMA_ADD_CELL(I,IUM_D_C),I=1,4)

20    CONTINUE
     
      IF(     WR_CELL_G > 0 ) THEN
        OPEN (WR_CELL_G)   
        WRITE(WR_CELL_G,*) MESH_P_OS(NZONET),NC_PL,NCELL_N2
        WRITE(WR_CELL_G,'(12I6)')IDCELL
        CLOSE(WR_CELL_G) 
      ENDIF   
      RETURN 
99    WRITE(6,*)' End_of_file'
      CALL STOP_ALL('Stopped in READ_CELL_NR_N0')
      END SUBROUTINE READ_CELL_NR_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CELL_DEF0 *** CELL_DEF0  *** CELL_DEF0 *** CELL_DEF0 ***
C CELL DEFINITION WITH IND_CELL=0
      SUBROUTINE CELL_DEF0(ICA,ICB)
      USE GEOMETRY_PL
      USE PHYSICAL_CELL

      INTEGER,INTENT(IN) :: ICA
      INTEGER,INTENT(OUT):: ICB

      INTEGER ::  IG,IC,IZ,I,J,K

      IC = ICA - 1
      DO 10 IZ=0,NZONET-1
      DO 10 K =0,ZON_TORO(IZ)-1
      DO 10 J =0,ZON_POLO(IZ)-1
      DO 10 I =0,ZON_RADI(IZ)-1
         IG =I+(J+K*ZON_POLO(IZ))*ZON_RADI(IZ)+MESH_P_OS(IZ)
         IF(IDCELL(IG) == 0) THEN 
           IC = IC + 1
           IDCELL(IG) = IC
         ENDIF 
10    CONTINUE
      ICB = IC
      RETURN
      END SUBROUTINE CELL_DEF0
C CELL_DEF1 *** CELL_DEF1  *** CELL_DEF1 *** CELL_DEF1 ***
C CELL DEFINITION WITH IND_CELL=1
      SUBROUTINE CELL_DEF1(IU,ICA,ICB)

      INTEGER,INTENT(IN) :: IU,ICA
      INTEGER,INTENT(OUT):: ICB

      INTEGER ::  IG,IC,IZ,IRED,ICN,NC
     .,           NR1,NR2,NP1,NP2,NT1,NT2

      READ(IU,*) NC
      IC  = ICA 
      DO 15 ICN=1,NC
        READ(IU,*)IRED,IZ,NR1,NR2,NP1,NP2,NT1,NT2 
        CALL CELL_DEF2(IZ,  NR1,NR2,NR2-NR1 
     .,                     NR1,NR2,NR2-NR1  
     .,                     NR1,NR2,NR2-NR1  
     .,                IC , ICB               )
        IC = ICB + 1
15    CONTINUE 
      RETURN
      END SUBROUTINE CELL_DEF1
C CELL_DEF2 *** CELL_DEF2  *** CELL_DEF2 *** CELL_DEF2 ***
C CELL DEFINITION WITH IND_CELL=2
      SUBROUTINE CELL_DEF2(IZ,NR1,NR2,NDR
     .,                       NP1,NP2,NDP
     .,                       NT1,NT2,NDT
     .,                       ICA,ICB      )
      USE GEOMETRY_PL
      USE PHYSICAL_CELL

      INTEGER,INTENT(IN) :: IZ,NR1,NR2,NDR,NP1,NP2,NDP
     .,                     NT1,NT2,NDT,ICA
      INTEGER,INTENT(OUT):: ICB
      INTEGER            :: NI,NJ,NK,I,J,K,DEF,IG,IC
      LOGICAL            :: OUT_OF_RANGE

      IF( OUT_OF_RANGE(IZ       ,0,NZONET-1    ) .OR.
     .    OUT_OF_RANGE(NDR,1    ,NR2-NR1       ) .OR.
     .    OUT_OF_RANGE(NR1,0    ,NR2-1         ) .OR.
     .    OUT_OF_RANGE(NR2,NR1+1,SRF_RADI(IZ)-1) .OR.
     .    OUT_OF_RANGE(NDP,1    ,NP2-NP1       ) .OR.
     .    OUT_OF_RANGE(NP1,0    ,NP2-1         ) .OR.
     .    OUT_OF_RANGE(NP2,NP1+1,SRF_POLO(IZ)-1) .OR.
     .    OUT_OF_RANGE(NDT,1    ,NT2-NT1       ) .OR.
     .    OUT_OF_RANGE(NT1,0    ,NT2-1         ) .OR.
     .    OUT_OF_RANGE(NT2,NT1+1,SRF_TORO(IZ)-1) ) THEN 
         CALL WRMESS(' out of range')
         CALL STOP_ALL('Stopped in CELL_DEF2')
      ENDIF 

      IC  = ICA - 1
      DO 18 NI=NR1,NR2-1,NDR
      DO 18 NJ=NP1,NP2-1,NDP
      DO 18 NK=NT1,NT2-1,NDT
        IC = IC + 1
        DEF = 0
        DO 19 I=NI,NI+NDR-1
        DO 19 J=NJ,NJ+NDP-1
        DO 19 K=NK,NK+NDT-1
          IG =I+(J+K*ZON_POLO(IZ))*ZON_RADI(IZ)+MESH_P_OS(IZ)
          IF(IDCELL(IG) == 0) THEN 
            IDCELL(IG) = IC
            DEF = 1
          ENDIF 
19      CONTINUE
        IF(DEF == 0) IC = IC - 1
18    CONTINUE
      ICB = IC

      RETURN
      END SUBROUTINE CELL_DEF2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UPDATA_CELL_NR_N0()
      USE GEOMETRY_PL
      USE PHYSICAL_CELL

      INTEGER ::  IG,IZ,I,J,K,IC

      ALLOCATE ( VOLCEL_N0(NCELL_N1:NCELL_N2) )
      VOLCEL_N0 = 0.
      DO 10 IZ=0,NZONET-1
      DO 10 K =0,ZON_TORO(IZ)-1
      DO 10 J =0,ZON_POLO(IZ)-1
      DO 10 I =0,ZON_RADI(IZ)-1
         IG =I+(J+K*ZON_POLO(IZ))*ZON_RADI(IZ)+MESH_P_OS(IZ)
         IC = IDCELL(IG)

         IF(  IC >= NCELL_N1 .AND. IC <= NCELL_N2)
     .   VOLCEL_N0(IC) = VOLCEL_N0(IC) + VOL3D(IG)
10    CONTINUE
      RETURN
      END SUBROUTINE UPDATA_CELL_NR_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READ_SOURCE_N0(IUNIT)
      INTEGER, INTENT(IN) :: IUNIT
      CHARACTER*72 :: READFI
      INTEGER      :: I,J,IU
      REAL*8       :: SP_TOT
  
      CALL SCRAPE(IUNIT,READFI,*99)
      READ(READFI,*) N0S_DEF,NS_PLACE,NSSIDE
      IF    (N0S_DEF.GT.0) THEN
        CALL WRMESS('N0_Source defined with N0S_DEF>0 not included')
        GOTO 100 
      ELSEIF(N0S_DEF.EQ.0) THEN
        CALL WRMESS('Neutral Source from EMC3:')
        IU = 99
        OPEN(IU,FILE='PARTICLE_DEPO')
        READ(IU,*,END=99) NSSTOR


        ALLOCATE( IJK_N0(NSSTOR)
     .,           ISD_N0(NSSTOR)
     .,           IZO_N0(NSSTOR)
     .,           IR0_N0(NSSTOR)
     .,           IP0_N0(NSSTOR)
     .,           IT0_N0(NSSTOR)
     .,           GEWICH(NSSTOR)
     .           )

        GEWICH = 0.
        DO I=1,NSSTOR
        READ(IU,*)J,IJK_N0(I),ISD_N0(I),IZO_N0(I) 
     +,             IR0_N0(I),IP0_N0(I),IT0_N0(I),GEWICH(I)
        IF(I > 1) GEWICH(I) = GEWICH(I) + GEWICH(I-1)
        ENDDO
        CLOSE(IU)
        SP_TOT = GEWICH(NSSTOR)
        GEWICH = GEWICH/SP_TOT

        NS_PLACE = 0

      ELSE
        CALL SCRAPE(IUNIT,READFI,*99)

        IU = 99
        CALL WRMESS('Neutral Source from file:')
        CALL WRMESS(READFI)
        OPEN(IU,FILE=READFI)
        REWIND(IU)

        READ(IU,*,END=99) NSSTOR

        ALLOCATE( NADDSF(NSSTOR)
     .,           XSOURP(NSSTOR)
     .,           YSOURP(NSSTOR)
     .,           ZSOURP(NSSTOR)
     .,           GEWICH(NSSTOR)
     .           )

        READ(IU,*)
        GEWICH = 0.
        DO I=1,NSSTOR
        READ(IU,*)J,NADDSF(I),GEWICH(I)
     +,             XSOURP(I),YSOURP(I),ZSOURP(I)
        IF(I > 1) GEWICH(I) = GEWICH(I) + GEWICH(I-1)
        ENDDO
        CLOSE(IU)

        SP_TOT = GEWICH(NSSTOR)
        GEWICH = GEWICH/SP_TOT

        NS_PLACE = 0
      ENDIF
      IF(NS_PLACE.EQ.0) THEN
        CALL WRMESS('Source on add. sf. of geometric part')
        IF(NSSIDE.EQ.-1) THEN
         CALL WRMESS('Source on negative side of the srf.')
        ELSE
         CALL WRMESS('Source on positive side of the srf.')
        ENDIF
      ELSE
        CALL WRMESS('Source on add. sf. of atomic part')
      ENDIF

      RETURN
99    CALL WRMESS('End_of_file')
100   CALL STOP_ALL('Stopped in READ_SOURCE_N0')
      END SUBROUTINE READ_SOURCE_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UPDATA_SOURCE_N0()
      END SUBROUTINE UPDATA_SOURCE_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READ_ADD_SF_N0(IUNIT)
      INTEGER, INTENT(IN) :: IUNIT
      CHARACTER*72 :: READFI,FILNM
      INTEGER,PARAMETER :: IU1=98,IU2=99
      INTEGER           :: NTR,I,J,MP
C *** Additional surfaces
      CALL SCRAPE(IUNIT,READFI,*99)
      FILNM = READFI 
    
      NTRIANG = 0    
      OPEN(IU1,FILE=FILNM)
      CALL SCRAPE(IU1,READFI,*99)
      READ(READFI,*) NPLATE

      WRITE(6,*)'Number of the plates:',NPLATE
      IF(NPLATE <= 0) RETURN 

      ALLOCATE( 
     .          ISR_TYPE(NPLATE  )
     .,         NSTUCK  (NPLATE,2)
     .         ) 

      DO 10 I=1,NPLATE
       CALL SCRAPE(IU1,READFI,*99)
       READ(READFI,*) NTR,ISR_TYPE(I)

       NSTUCK(I,1) = NTRIANG + 1
       IF(NTR.LT.0) THEN
         CALL WRMESS('NTR.LT.0 to be written')
         GOTO 100
       ELSEIF(NTR.EQ.0) THEN
         CALL SCRAPE(IU1,READFI,*99)
         OPEN(IU2,FILE=READFI)
           CALL LIM_ADD(IU2,0,I,NTRIANG)
         CLOSE(IU2)
       ELSE
         DO J=1,NTR
           CALL SCRAPE(IU1,READFI,*99)
           CALL SCRAPE(IU1,READFI,*99)
         ENDDO 
         NTRIANG = NTRIANG + NTR
       ENDIF
       NSTUCK(I,2) = NTRIANG
10    CONTINUE      

      REWIND(IU1) 

      MP = 0 
      DO 
        READ(IU1,'(A)') READFI
        IF(READFI(1:1) /= '*') EXIT
      END DO
      READ(READFI,*)  ! NPLATE

      WRITE(6,*)'Plate  Triangles Type  <------->'
      
      ALLOCATE( 
     .          X_TRIA(3,NTRIANG)
     .,         Y_TRIA(3,NTRIANG)
     .,         Z_TRIA(3,NTRIANG)
     .         ) 

      DO 20 I=1,NPLATE
       WRITE(6,'(I4,6x,i5,4x,i3,2i6)')
     &        I,NTR,ISR_TYPE(I),NSTUCK(I,1),NSTUCK(I,2)
       DO 
        READ(IU1,'(A)') READFI
        IF(READFI(1:1) /= '*') EXIT
       END DO

       READ(READFI,*) NTR,J 

       IF(NTR.EQ.0) THEN
         DO 
          READ(IU1,'(A)') READFI
          IF(READFI(1:1) /= '*') EXIT
         END DO
         OPEN(IU2,FILE=READFI)
         REWIND(IU2)
           CALL LIM_ADD(IU2,1,I,MP)
         CLOSE(IU2)
       ELSE
         DO J=1,NTR
           MP = MP + 1
           DO 
            READ(IU1,'(A)') READFI
            IF(READFI(1:1) /= '*') EXIT
           END DO
           READ(READFI,*) X_TRIA(1,MP),Y_TRIA(1,MP),Z_TRIA(1,MP)
     &,                   X_TRIA(2,MP),Y_TRIA(2,MP),Z_TRIA(2,MP)
           DO 
            READ(IU1,'(A)') READFI
            IF(READFI(1:1) /= '*') EXIT
           END DO
           READ(READFI,*) X_TRIA(3,MP),Y_TRIA(3,MP),Z_TRIA(3,MP)
         ENDDO 
       ENDIF
20    CONTINUE      
      CLOSE(IU1)

      RETURN
99    WRITE(6,*)'End_of_file'           
100   CALL STOP_ALL('Stopped in READ_ADD_SF_N0')
      RETURN
      END SUBROUTINE READ_ADD_SF_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE UPDATA_ADD_SF_N0()
      INTEGER :: I,MP,IBEG,NPS1,NPS2,J,I3E,M,IPOSI
      REAL*8  :: V2,V,VVS,F,VX,VY,VZ,D,ANG,XP,YP,ZP,DIS_MAX
      INTEGER,DIMENSION(:),ALLOCATABLE :: IND
      REAL*8, PARAMETER                :: ANG_LIM=0.2,TOLER=1.E-12
      REAL*8, DIMENSION(3)             :: DIST

      IF(NTRIANG <= 0) RETURN 

      ALLOCATE ( 
     .          A_PLAT  (NPLATE )
     .,         A_TRIA  (NTRIANG)
     .,         V_TRIA_X(NTRIANG)
     .,         V_TRIA_Y(NTRIANG)
     .,         V_TRIA_Z(NTRIANG)
     .,         NP_NUM  (NTRIANG)
     .         )  

      DO I=1,NPLATE
        A_PLAT(I)   =  0.
        DO MP=NSTUCK(I,1),NSTUCK(I,2)
           V_TRIA_X(MP) =
     &     (Y_TRIA(2,MP)-Y_TRIA(1,MP))*(Z_TRIA(3,MP)-Z_TRIA(1,MP))
     &    -(Z_TRIA(2,MP)-Z_TRIA(1,MP))*(Y_TRIA(3,MP)-Y_TRIA(1,MP))

           V_TRIA_Y(MP) =
     &     (Z_TRIA(2,MP)-Z_TRIA(1,MP))*(X_TRIA(3,MP)-X_TRIA(1,MP))
     &    -(X_TRIA(2,MP)-X_TRIA(1,MP))*(Z_TRIA(3,MP)-Z_TRIA(1,MP))

           V_TRIA_Z(MP) =
     &     (X_TRIA(2,MP)-X_TRIA(1,MP))*(Y_TRIA(3,MP)-Y_TRIA(1,MP))
     &    -(Y_TRIA(2,MP)-Y_TRIA(1,MP))*(X_TRIA(3,MP)-X_TRIA(1,MP))

           V2 = V_TRIA_X(MP)**2 + V_TRIA_Y(MP)**2 + V_TRIA_Z(MP)**2
           IF(V2 == 0.) THEN
             CALL WRMESS('AN Add. SR NOT EXIST!')
             CALL STOP_ALL('STOPPED IN UPDATA_ADD_SF_N0')
           ENDIF
           V = SQRT(V2)

           A_TRIA(MP)   = 0.5*V
           A_PLAT(I)    = A_PLAT(I) + A_TRIA(MP)
           VVS = 1./V
           V_TRIA_X(MP) = V_TRIA_X(MP)*VVS
           V_TRIA_Y(MP) = V_TRIA_Y(MP)*VVS
           V_TRIA_Z(MP) = V_TRIA_Z(MP)*VVS

           NP_NUM(MP) = I
        ENDDO
      ENDDO
C Group 
      ALLOCATE ( 
     .          NP_GRP  (  NTRIANG)
     .,         NPS_GRP (0:NTRIANG)
     .,         NPN_GRP (  NTRIANG)
     .,         X_SORT  (  NTRIANG)
     .,         Y_SORT  (  NTRIANG)
     .,         Z_SORT  (  NTRIANG)
     .,         D_MAXI  (  NTRIANG)
     .,         D_MINI  (  NTRIANG)
     .,         X_SMIN  (  NTRIANG)
     .,         X_SMAX  (  NTRIANG)
     .,         Y_SMIN  (  NTRIANG)
     .,         Y_SMAX  (  NTRIANG)
     .,         Z_SMIN  (  NTRIANG)
     .,         Z_SMAX  (  NTRIANG)
     .,       V_SORT_X  (  NTRIANG)
     .,       V_SORT_Y  (  NTRIANG)
     .,       V_SORT_Z  (  NTRIANG)
     .         )  

      ALLOCATE ( IND(NTRIANG) )
      IND        = 0
      NPSORT     = 0
      NPS_GRP(0) = 0

      IBEG = 0
      SCAN_LOOP : DO
     
      IBEG_LOOP : DO
        IBEG = IBEG + 1
        IF    (IBEG > NTRIANG) THEN
          EXIT SCAN_LOOP
        ELSEIF(IND(IBEG) == 0) THEN
          EXIT IBEG_LOOP
        ENDIF 
      ENDDO IBEG_LOOP

      NPSORT = NPSORT + 1

      IND(IBEG)      = 1
      NP_GRP(NPSORT) = 1
      IPOSI          = NPS_GRP(NPSORT-1) + 1
      NPN_GRP(IPOSI) = IBEG

      X_SMIN(NPSORT)=MIN(X_TRIA(1,IBEG),X_TRIA(2,IBEG),X_TRIA(3,IBEG))
      X_SMAX(NPSORT)=MAX(X_TRIA(1,IBEG),X_TRIA(2,IBEG),X_TRIA(3,IBEG))
      Y_SMIN(NPSORT)=MIN(Y_TRIA(1,IBEG),Y_TRIA(2,IBEG),Y_TRIA(3,IBEG))
      Y_SMAX(NPSORT)=MAX(Y_TRIA(1,IBEG),Y_TRIA(2,IBEG),Y_TRIA(3,IBEG))
      Z_SMIN(NPSORT)=MIN(Z_TRIA(1,IBEG),Z_TRIA(2,IBEG),Z_TRIA(3,IBEG))
      Z_SMAX(NPSORT)=MAX(Z_TRIA(1,IBEG),Z_TRIA(2,IBEG),Z_TRIA(3,IBEG))

      X_SORT(NPSORT) = (X_TRIA(1,IBEG)+X_TRIA(2,IBEG)+X_TRIA(3,IBEG))/3.
      Y_SORT(NPSORT) = (Y_TRIA(1,IBEG)+Y_TRIA(2,IBEG)+Y_TRIA(3,IBEG))/3.
      Z_SORT(NPSORT) = (Z_TRIA(1,IBEG)+Z_TRIA(2,IBEG)+Z_TRIA(3,IBEG))/3.

      V_SORT_X(NPSORT) = V_TRIA_X(IBEG)
      V_SORT_Y(NPSORT) = V_TRIA_Y(IBEG)
      V_SORT_Z(NPSORT) = V_TRIA_Z(IBEG)

      DO 30 I=IBEG+1,NTRIANG
        IF(IND(I).NE.0) GOTO 30

        F = 1./DFLOAT(NP_GRP(NPSORT))
        VX = V_SORT_X(NPSORT)*F
        VY = V_SORT_Y(NPSORT)*F
        VZ = V_SORT_Z(NPSORT)*F

        V  = DSQRT(VX**2+VY**2+VZ**2)

        D  = (VX*V_TRIA_X(I) + VY*V_TRIA_Y(I) + VZ*V_TRIA_Z(I))/V

        D  = 0.9999999999*D

       ANG =ACOS(D)
       IF(ABS(ANG).LT.ANG_LIM) THEN
        XP = X_SORT(NPSORT)*F
        YP = Y_SORT(NPSORT)*F
        ZP = Z_SORT(NPSORT)*F

        DO J=1,3
         DIST(J)=ABS( (X_TRIA(J,I)-XP)*VX
     &               +(Y_TRIA(J,I)-YP)*VY
     &               +(Z_TRIA(J,I)-ZP)*VZ )
        ENDDO
         DIS_MAX=MAX(DIST(1),DIST(2),DIST(3))

         IF(DIS_MAX.LT.ANG_LIM*SQRT(A_TRIA(I))) THEN
           NP_GRP(NPSORT) = NP_GRP(NPSORT) + 1
           IPOSI          = NP_GRP(NPSORT) + NPS_GRP(NPSORT-1)
           NPN_GRP(IPOSI) = I
           IND(I)         = 1

           X_SMIN(NPSORT) =MIN(X_SMIN(NPSORT),
     &                     X_TRIA(1,I),X_TRIA(2,I),X_TRIA(3,I) )
           X_SMAX(NPSORT) =MAX(X_SMAX(NPSORT),
     &                     X_TRIA(1,I),X_TRIA(2,I),X_TRIA(3,I) )
           Y_SMIN(NPSORT) =MIN(Y_SMIN(NPSORT),
     &                     Y_TRIA(1,I),Y_TRIA(2,I),Y_TRIA(3,I) )
           Y_SMAX(NPSORT) =MAX(Y_SMAX(NPSORT),
     &                     Y_TRIA(1,I),Y_TRIA(2,I),Y_TRIA(3,I) )
           Z_SMIN(NPSORT) =MIN(Z_SMIN(NPSORT),
     &                     Z_TRIA(1,I),Z_TRIA(2,I),Z_TRIA(3,I) )
           Z_SMAX(NPSORT) =MAX(Z_SMAX(NPSORT),
     &                     Z_TRIA(1,I),Z_TRIA(2,I),Z_TRIA(3,I) )

           X_SORT(NPSORT) = X_SORT(NPSORT) +
     &                 (X_TRIA(1,I)+X_TRIA(2,I)+X_TRIA(3,I))/3.
           Y_SORT(NPSORT) = Y_SORT(NPSORT) +
     &                 (Y_TRIA(1,I)+Y_TRIA(2,I)+Y_TRIA(3,I))/3.
           Z_SORT(NPSORT) = Z_SORT(NPSORT) +
     &                 (Z_TRIA(1,I)+Z_TRIA(2,I)+Z_TRIA(3,I))/3.

           V_SORT_X(NPSORT) = V_SORT_X(NPSORT) + V_TRIA_X(I)
           V_SORT_Y(NPSORT) = V_SORT_Y(NPSORT) + V_TRIA_Y(I)
           V_SORT_Z(NPSORT) = V_SORT_Z(NPSORT) + V_TRIA_Z(I)

         ENDIF
       ENDIF
30    CONTINUE
      F = 1./FLOAT(NP_GRP(NPSORT))

      V_SORT_X(NPSORT) = V_SORT_X(NPSORT)*F
      V_SORT_Y(NPSORT) = V_SORT_Y(NPSORT)*F
      V_SORT_Z(NPSORT) = V_SORT_Z(NPSORT)*F

      V=DSQRT(
     &V_SORT_X(NPSORT)**2+V_SORT_Y(NPSORT)**2+V_SORT_Z(NPSORT)**2)
      V_SORT_X(NPSORT) = V_SORT_X(NPSORT)/V
      V_SORT_Y(NPSORT) = V_SORT_Y(NPSORT)/V
      V_SORT_Z(NPSORT) = V_SORT_Z(NPSORT)/V


      X_SORT(NPSORT)   = X_SORT(NPSORT)*F
      Y_SORT(NPSORT)   = Y_SORT(NPSORT)*F
      Z_SORT(NPSORT)   = Z_SORT(NPSORT)*F

      NPS_GRP(NPSORT) = NPS_GRP(NPSORT-1) + NP_GRP(NPSORT)

      ENDDO SCAN_LOOP

      DO I=1,NPSORT
        NPS1 = NPS_GRP(I-1) + 1
        NPS2 = NPS_GRP(I)

        D_MAXI(I) = TOLER
        D_MINI(I) =-TOLER

        DO J=NPS1,NPS2
         I3E = NPN_GRP(J)

         DO M=1,3
         DIST(M)= (X_TRIA(M,I3E)-X_SORT(I))*V_SORT_X(I)
     &           +(Y_TRIA(M,I3E)-Y_SORT(I))*V_SORT_Y(I)
     &           +(Z_TRIA(M,I3E)-Z_SORT(I))*V_SORT_Z(I)
         D_MAXI(I) = MAX(D_MAXI(I),DIST(M))
         D_MINI(I) = MIN(D_MINI(I),DIST(M))
         ENDDO

        ENDDO
      ENDDO
      DEALLOCATE ( IND )

      RETURN
      END SUBROUTINE UPDATA_ADD_SF_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LIM_ADD(IUNIT,IRUN,IPLA,TRI_ELE)
      USE GEOMETRY_PL
      USE PHYS_CONST

      INTEGER,INTENT(IN)    :: IUNIT,IRUN,IPLA !IPLA Plate Nr.
      INTEGER,INTENT(INOUT) :: TRI_ELE         ! Triangle nr.
C IRUN = 0: update TRI_ELE olny

      CHARACTER*72 ::  TITLE
 
      REAL*8, DIMENSION(:  ), ALLOCATABLE :: PHI
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: R, Z
      REAL*8  :: R0,Z0,SIN1,SIN2,COS1,COS2,PHI_PERIOD
     .,          PHI_A,PHI_B
      INTEGER :: N,NPOLO,MPER,I,J,IPE
C ---------------------------------------------------------------------
       READ(IUNIT,'(A72)')TITLE
       READ(IUNIT,*)N,NPOLO,MPER,R0,Z0   
         
       write(6,*) TRI_ELE
       IF(IRUN == 0) THEN
         WRITE(6,'(A72)')TITLE
         WRITE(6,'(3I6,3f6.2)') N,NPOLO,MPER,R0,Z0       
         IF( N<2 .OR. MPER<1 .OR. NPOLO<2 ) THEN
           WRITE(6,*)'N<2 .OR. MPER<1 .OR. NPOLO<2'
           call stop_all('Stopped in LIM_ADD')
         ENDIF 
         PHI_PERIOD=PI2/FLOAT(MPER)

         ALLOCATE (PHI(N))
         DO I=1,N
           READ(IUNIT,*) PHI(I)
           DO J=1,NPOLO;READ(IUNIT,*);ENDDO
         ENDDO
         PHI = PHI*PI/180.
C move the limiter into the computation domain if necessary
         DO I=1,N-1
           DO IPE=1,MPER
            PHI_A=PHI_PLANE(0                  )-PHI_PERIOD*(IPE-1)
            PHI_B=PHI_PLANE(PHI_PL_OS(NZONET)-1)-PHI_PERIOD*(IPE-1)
            IF( (PHI(I  )-PHI_A)*(PHI(I  )-PHI_B) < 0.
     ..OR.      (PHI(I+1)-PHI_A)*(PHI(I+1)-PHI_B) < 0.) THEN
               TRI_ELE = TRI_ELE+(NPOLO-1)*2
               EXIT
            ENDIF
           ENDDO
         ENDDO
         DEALLOCATE (PHI)
         RETURN 
       ENDIF 
 
      ALLOCATE ( PHI(2), R(2,NPOLO), Z(2,NPOLO) )
      PHI_PERIOD=PI2/FLOAT(MPER)

      READ(IUNIT,*) PHI(1)
      PHI(1) = PHI(1)*PI/180.
      DO J=1,NPOLO
        READ(IUNIT,*)R(1,J),Z(1,J)
        R(1,J) = R(1,J) + R0
        Z(1,J) = Z(1,J) + Z0
      ENDDO  

      DO 10 I = 2,N
        READ(IUNIT,*) PHI(2)
        PHI(2) = PHI(2)*PI/180.

        DO J=1,NPOLO
         READ(IUNIT,*)R(2,J),Z(2,J)
         R(2,J) = R(2,J) + R0
         Z(2,J) = Z(2,J) + Z0
        ENDDO  
        
C move the limiter into the computation domain if necessary
        DO IPE=1,MPER
           PHI_A=PHI_PLANE(0                  )-PHI_PERIOD*(IPE-1)
           PHI_B=PHI_PLANE(PHI_PL_OS(NZONET)-1)-PHI_PERIOD*(IPE-1)
           IF( (PHI(1)-PHI_A)*(PHI(1)-PHI_B) < 0.
     ..OR.     (PHI(2)-PHI_A)*(PHI(2)-PHI_B) < 0.) EXIT
        ENDDO

        IF(IPE <= MPER) THEN
          COS1=COS(PHI(1)+PHI_PERIOD*(IPE-1))
          COS2=COS(PHI(2)+PHI_PERIOD*(IPE-1))
          SIN1=SIN(PHI(1)+PHI_PERIOD*(IPE-1))
          SIN2=SIN(PHI(2)+PHI_PERIOD*(IPE-1))
          DO 20 J=1,NPOLO-1
           TRI_ELE = TRI_ELE + 1
           X_TRIA(1,TRI_ELE) = R(1,J  )*COS1
           Y_TRIA(1,TRI_ELE) = R(1,J  )*SIN1
           Z_TRIA(1,TRI_ELE) = Z(1,J  )
 
           X_TRIA(2,TRI_ELE) = R(2,J+1)*COS2
           Y_TRIA(2,TRI_ELE) = R(2,J+1)*SIN2
           Z_TRIA(2,TRI_ELE) = Z(2,J+1)
 
           X_TRIA(3,TRI_ELE) = R(2,J  )*COS2
           Y_TRIA(3,TRI_ELE) = R(2,J  )*SIN2
           Z_TRIA(3,TRI_ELE) = Z(2,J  )
 
           TRI_ELE = TRI_ELE + 1

           X_TRIA(1,TRI_ELE) = R(1,J  )*COS1
           Y_TRIA(1,TRI_ELE) = R(1,J  )*SIN1
           Z_TRIA(1,TRI_ELE) = Z(1,J  )
 
           X_TRIA(2,TRI_ELE) = R(1,J+1)*COS1
           Y_TRIA(2,TRI_ELE) = R(1,J+1)*SIN1
           Z_TRIA(2,TRI_ELE) = Z(1,J+1)
 
           X_TRIA(3,TRI_ELE) = R(2,J+1)*COS2
           Y_TRIA(3,TRI_ELE) = R(2,J+1)*SIN2
           Z_TRIA(3,TRI_ELE) = Z(2,J+1)
 
 20       CONTINUE 
        ENDIF 
 
        PHI(1) = PHI(2)
        DO J=1,NPOLO
         R(1,J) = R(2,J) 
         Z(1,J) = Z(2,J)
        ENDDO 
 10   CONTINUE
      DEALLOCATE ( PHI, R, Z )
      RETURN 
      END SUBROUTINE LIM_ADD

      END MODULE NEUTRAL_TRANSPORT
