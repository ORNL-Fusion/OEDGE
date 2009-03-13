C EIRENE07 COMPILATION
C ===== SOURCE: NEUTRAL.f
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
C ===== SOURCE: broadcast_n0.f
      SUBROUTINE BROADCAST_NEUTRAL 
      USE PARALLEL 
      USE GEOMETRY_PL
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT 

C1. surface
      ENTRY BROADCAST_SF_N0()

      CALL BROADCAST_INTE_S(NONTRT_N0)
      IF(NONTRT_N0 > 0) THEN
        IF(MYPE > 0) ALLOCATE(
     I    IR0NR_N0(NONTRT_N0),IZONR_N0(NONTRT_N0),INDINR_N0(NONTRT_N0)
     I,   IP1NR_N0(NONTRT_N0),IP2NR_N0(NONTRT_N0)
     I,   IT1NR_N0(NONTRT_N0),IT2NR_N0(NONTRT_N0)
     &                            )

          CALL BROADCAST_INTE( IR0NR_N0,NONTRT_N0)
          CALL BROADCAST_INTE( IZONR_N0,NONTRT_N0)
          CALL BROADCAST_INTE(INDINR_N0,NONTRT_N0)
          CALL BROADCAST_INTE( IP1NR_N0,NONTRT_N0)
          CALL BROADCAST_INTE( IP2NR_N0,NONTRT_N0)
          CALL BROADCAST_INTE( IT1NR_N0,NONTRT_N0)
          CALL BROADCAST_INTE( IT2NR_N0,NONTRT_N0)
      ENDIF 
   
      CALL BROADCAST_INTE_S(NONTPT_N0)
      IF(NONTPT_N0 > 0) THEN
        IF(MYPE > 0) ALLOCATE(
     I    IP0NP_N0(NONTPT_N0),IZONP_N0(NONTPT_N0),INDINP_N0(NONTPT_N0)
     I,   IR1NP_N0(NONTPT_N0),IR2NP_N0(NONTPT_N0)
     I,   IT1NP_N0(NONTPT_N0),IT2NP_N0(NONTPT_N0)
     &                            )

          CALL BROADCAST_INTE( IP0NP_N0,NONTPT_N0)
          CALL BROADCAST_INTE( IZONP_N0,NONTPT_N0)
          CALL BROADCAST_INTE(INDINP_N0,NONTPT_N0)
          CALL BROADCAST_INTE( IR1NP_N0,NONTPT_N0)
          CALL BROADCAST_INTE( IR2NP_N0,NONTPT_N0)
          CALL BROADCAST_INTE( IT1NP_N0,NONTPT_N0)
          CALL BROADCAST_INTE( IT2NP_N0,NONTPT_N0)
      ENDIF 

      CALL BROADCAST_INTE_S(NONTTT_N0)
      IF(NONTTT_N0 > 0) THEN
        IF(MYPE > 0) ALLOCATE(
     I    IT0NT_N0(NONTTT_N0),IZONT_N0(NONTTT_N0),INDINT_N0(NONTTT_N0)
     I,   IR1NT_N0(NONTTT_N0),IR2NT_N0(NONTTT_N0)
     I,   IP1NT_N0(NONTTT_N0),IP2NT_N0(NONTTT_N0)
     &                            )

          CALL BROADCAST_INTE( IT0NT_N0,NONTTT_N0)
          CALL BROADCAST_INTE( IZONT_N0,NONTTT_N0)
          CALL BROADCAST_INTE(INDINT_N0,NONTTT_N0)
          CALL BROADCAST_INTE( IR1NT_N0,NONTTT_N0)
          CALL BROADCAST_INTE( IR2NT_N0,NONTTT_N0)
          CALL BROADCAST_INTE( IP1NT_N0,NONTTT_N0)
          CALL BROADCAST_INTE( IP2NT_N0,NONTTT_N0)
      ENDIF 
      RETURN

C2. Cell nr.
      ENTRY BROADCAST_CELL_NR_N0()

      CALL BROADCAST_INTE_S(NCTYPE_DEF)
      CALL BROADCAST_INTE_S(NCELL_N1  )
      CALL BROADCAST_INTE_S(NCELL_N2  )
    
      IF( NCTYPE_DEF > 0) THEN
        IF(MYPE > 0) ALLOCATE(
     I           NCELL_AB       (0:NCTYPE_DEF)
     R,          PLASMA_ADD_CELL(4,NCTYPE_DEF) )

        CALL BROADCAST_INTE(NCELL_AB,NCTYPE_DEF+1)
        CALL BROADCAST_REAL(PLASMA_ADD_CELL,4*NCTYPE_DEF)
      ENDIF 
      CALL BROADCAST_INTE(IDCELL,MESH_P_OS(NZONET))
      RETURN

C3. Source distribution
      ENTRY BROADCAST_SOURCE_N0()

      CALL BROADCAST_INTE_S(N0S_DEF )
      CALL BROADCAST_INTE_S(NS_PLACE)
      CALL BROADCAST_INTE_S(NSSIDE  ) 
      CALL BROADCAST_INTE_S(NSSTOR  ) 

      IF    (N0S_DEF  > 0) THEN
      ELSEIF(N0S_DEF == 0) THEN
         IF(MYPE > 0) ALLOCATE (
     .            IJK_N0(NSSTOR)
     .,           ISD_N0(NSSTOR)
     .,           IZO_N0(NSSTOR)
     .,           IR0_N0(NSSTOR)
     .,           IP0_N0(NSSTOR)
     .,           IT0_N0(NSSTOR)
     .,           GEWICH(NSSTOR)
     .           )

         CALL BROADCAST_INTE(IJK_N0,NSSTOR)
         CALL BROADCAST_INTE(ISD_N0,NSSTOR)
         CALL BROADCAST_INTE(IZO_N0,NSSTOR)
         CALL BROADCAST_INTE(IR0_N0,NSSTOR)
         CALL BROADCAST_INTE(IP0_N0,NSSTOR)
         CALL BROADCAST_INTE(IT0_N0,NSSTOR)
         CALL BROADCAST_REAL(GEWICH,NSSTOR)
      ELSE
         IF(MYPE > 0) ALLOCATE (
     .            NADDSF(NSSTOR)
     .,           XSOURP(NSSTOR)
     .,           YSOURP(NSSTOR)
     .,           ZSOURP(NSSTOR)
     .,           GEWICH(NSSTOR)
     .           )

         CALL BROADCAST_INTE(NADDSF,NSSTOR)
         CALL BROADCAST_REAL(XSOURP,NSSTOR)
         CALL BROADCAST_REAL(YSOURP,NSSTOR)
         CALL BROADCAST_REAL(ZSOURP,NSSTOR)
         CALL BROADCAST_REAL(GEWICH,NSSTOR)
      ENDIF 
      RETURN 

C4. Additional surfaces
      ENTRY BROADCAST_ADD_SF_N0()

      CALL BROADCAST_INTE_S(NPLATE  )
      CALL BROADCAST_INTE_S(NTRIANG )

      IF(NPLATE  >  0) THEN
         IF(MYPE > 0) ALLOCATE (
     .          ISR_TYPE(NPLATE   )
     .,         NSTUCK  (NPLATE,2 )
     .,         X_TRIA  (3,NTRIANG)
     .,         Y_TRIA  (3,NTRIANG)
     .,         Z_TRIA  (3,NTRIANG)
     .         )
         
         CALL BROADCAST_INTE(ISR_TYPE, NPLATE )
         CALL BROADCAST_INTE(NSTUCK ,2*NPLATE )
         CALL BROADCAST_REAL(X_TRIA ,3*NTRIANG)
         CALL BROADCAST_REAL(Y_TRIA ,3*NTRIANG)
         CALL BROADCAST_REAL(Z_TRIA ,3*NTRIANG)
      ENDIF 
      RETURN 

      END SUBROUTINE BROADCAST_NEUTRAL

C ===== SOURCE: initial_n0.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: INITIAL_NEUTRAL                                          C
C Function: set up all data necessary for running Eirene             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INITIAL_NEUTRAL()
      USE PARALLEL 
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE 
      INTEGER,PARAMETER :: UNIT=4
      INTRINSIC SUM
C----------------------------------------------------------------------  
C 1. Geometry, mesh, grid
C NONTRT_N0 
C IZONR_N0 IR0NR_N0 IP1NR_N0 IP2NR_N0 IT1NR_N0 IT2NR_N0 INDINR_N0
C NONTPT_N0 
C IZONP_N0 IR0NP_N0 IP1NP_N0 IP2NP_N0 IT1NP_N0 IT2NP_N0 INDINP_N0
C NONTPT_N0 
C IZONT_N0 IR0NT_N0 IP1NT_N0 IP2NT_N0 IT1NT_N0 IT2NT_N0 INDINT_N0
      IF(MYPE == 0) CALL READ_SF_N0(UNIT)  
      IF (NPRS > 1) CALL BROADCAST_SF_N0() 
C IDSURR IDSURP IDSURT updated in OPEN_NEUTRAL
C     CALL UPDATA_SF_N0()

C 2. Cell for neutral
C NCELL_N1 NCELL_N2 NCTYPE_DEF NCELL_AB
C IDCELL 
C PLASMA_ADD_CELL
      IF(MYPE == 0) CALL READ_CELL_NR_N0(UNIT)
      IF (NPRS > 1) CALL BROADCAST_CELL_NR_N0() 
C VOLCEL_N0
      CALL UPDATA_CELL_NR_N0()
      IF(MYPE == 0)
     .WRITE(6,*)'  Total volume(+add. cells):'
     .,VOLCEL(0)+SUM(VOLCEL_N0),' cm**3 cm'

C3. Neutral Source distribution
C NSSTOR N0S_DEF NS_PLACE NSSIDE
C NADDSF IJK_N0 IZO_N0 ISD_N0 IR0_N0 IP0_N0 IT0_N0
C GEWICH XSOURP YSOURP ZSOURP
      IF(MYPE == 0) CALL READ_SOURCE_N0(UNIT)
      IF (NPRS > 1) CALL BROADCAST_SOURCE_N0()
      CALL UPDATA_SOURCE_N0()

C4. additional surfaces
C NPLATE   NTRIANG
C ISR_TYPE NSTUCK
C X_TRIA   Y_TRIA   Z_TRIA
      IF(MYPE == 0) CALL READ_ADD_SF_N0(UNIT)
      IF (NPRS > 1) CALL BROADCAST_ADD_SF_N0()
C NPSORT NP_NUM NP_GRP NPS_GRP NPN_GRP
C A_PLAT A_TRIA V_TRIA_X V_TRIA_Y V_TRIA_Z X_SORT Y_SORT Z_SORT
C D_MAXI D_MINI X_SMIN X_SMAX Y_SMIN Y_SMAX Z_SMIN Z_SMAX
C V_SORT_X V_SORT_Y V_SORT_Z
      CALL UPDATA_ADD_SF_N0()
      IF(MYPE == 0) THEN
        WRITE(6,*)' Groups:',NPSORT
        WRITE(6,*)' Elements:'
        WRITE(6,'(18I4)') NP_GRP(1:NPSORT)
      ENDIF 

      RETURN
      END 
C ===== SOURCE: neutral.f
      SUBROUTINE RUN_NEUTRAL()
      USE PARALLEL       
      USE PHYSICAL_CELL
      USE BOUNDARY_COND
      USE SOURCE_V_PL
      USE PLASMA_PARAM
      USE ION_SPECIES
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE
      REAL*8, PARAMETER :: RLF=0.5
      INTEGER           :: IPLS,ITYP,ICA,ICB,I
      REAL*8            :: FLUX,FLUX_M,SP
      LOGICAL           :: LRUN_FIRST

C  Call Eirene
      CALL EIRENE_EMC3_ENTRY(0.,.FALSE.)
      IF(MYPE == 0 ) THEN

C  Post_processing
      FLUX  = PFLUX_TOTAL(1)
      FLUX_M= FLUX/(1.6E-19*1.67E-24*ATMASS(1))

      OUT_N0(1:NCELL_N2,4  ) = OUT_N0(1:NCELL_N2,4  )*FLUX_M
      OUT_N0(1:NCELL_N2,5:6) = OUT_N0(1:NCELL_N2,5:6)*FLUX
      DO I=1,NC_PL
      IF( ZMACH0(I) > 0.) THEN 
        OUT_N0(I,4) = OUT_N0(I,4)/SOUND_S(I)
      ELSE 
        OUT_N0(I,4) =-OUT_N0(I,4)/SOUND_S(I)
      ENDIF 
      ENDDO 
C First run
      DO I=1,MIN(100,NC_PL)
        LRUN_FIRST = VSOUP0(I,1).EQ.0.
        IF(.NOT.LRUN_FIRST) EXIT
      ENDDO

      IF( LRUN_FIRST ) THEN
         VSOUP0(1:NC_PL,1  ) = OUT_N0(1:NC_PL,1)
         VSOUE0(1:NC_PL,0,0) = OUT_N0(1:NC_PL,2)      
         VSOUE0(1:NC_PL,0,1) = OUT_N0(1:NC_PL,3)      
         VSOUM0(1:NC_PL,1  ) = OUT_N0(1:NC_PL,4)
      ELSE
         VSOUP0(1:NC_PL,1  ) =    RLF *VSOUP0(1:NC_PL,1  )  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,1  )      
         VSOUE0(1:NC_PL,0,0) =    RLF *VSOUE0(1:NC_PL,0,0)  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,2  )      
         VSOUE0(1:NC_PL,0,1) =    RLF *VSOUE0(1:NC_PL,0,1)  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,3  )      
         VSOUM0(1:NC_PL,1  ) =    RLF *VSOUM0(1:NC_PL,1  )  
     &+                       (1.-RLF)*OUT_N0(1:NC_PL,4  )      
      ENDIF 

      SP_MAIN_NEW = 0.
      ICA= NCELL_AB(0)+1
      ICB= NCELL_AB(1)
      DO I=ICA,ICB
        SP_MAIN_NEW    = SP_MAIN_NEW+OUT_N0(I,1)*VOLCEL_N0(I)
      ENDDO
      IF(SP_MAIN_OLD.NE.0.) THEN
        SP_MAIN_NEW = RLF*SP_MAIN_OLD+(1.-RLF)*SP_MAIN_NEW
      ENDIF
      SP_MAIN_OLD   = SP_MAIN_NEW

       OPEN (40)
       WRITE(40,'(6E12.4)') VSOUP0(1:NC_PL,1)
       CLOSE(40)

       OPEN (42)
       WRITE(42,'(6E12.4)') VSOUE0(1:NC_PL,0,0)
       CLOSE(42)

       OPEN (43)
       WRITE(43,'(6E12.4)') VSOUE0(1:NC_PL,0,1)
       CLOSE(43)

       OPEN (47)
       WRITE(47,'(6E12.4)') VSOUM0(1:NC_PL,1)
       CLOSE(47)

       OPEN (99,FILE='DENSITY_A')
       WRITE(99,'(6E12.4)') OUT_N0(1:NCELL_N2,5)
       CLOSE(99)  

       OPEN (99,FILE='DENSITY_M')
       WRITE(99,'(6E12.4)') OUT_N0(1:NCELL_N2,6)
       CLOSE(99)  
          
       WRITE(6,*)'    Cell    S_P(Amp.)  S_E_E(Watt) S_E_I(Watt)   
     .    SM_P      SM_N'
          SP          =0.
          SE_SUM_N0_E =0.
          SE_SUM_N0_I =0.
          SM_SUMP_ALL =0.
          SM_SUMN_ALL =0.
          DO I=1,NC_PL
           SP          = SP          + OUT_N0(I,1)*VOLCEL(I)
           SE_SUM_N0_E = SE_SUM_N0_E + OUT_N0(I,2)*VOLCEL(I)
           SE_SUM_N0_I = SE_SUM_N0_I + OUT_N0(I,3)*VOLCEL(I)
           IF( ZMACH0(I) > 0.) THEN 
           SM_SUMP_ALL = SM_SUMP_ALL + OUT_N0(I,4)*VOLCEL(I)
           ELSE
           SM_SUMN_ALL = SM_SUMN_ALL + OUT_N0(I,4)*VOLCEL(I)
           ENDIF
          ENDDO
          WRITE(6,'(2I6,5E12.4)') 1,NC_PL,SP*FLUX,SE_SUM_N0_E*FLUX
     .,   SE_SUM_N0_I*FLUX,SM_SUMP_ALL,SM_SUMN_ALL

          DO ITYP=1,NCTYPE_DEF
          SP          =0.
          SE_SUM_N0_E =0.
          SE_SUM_N0_I =0.
          ICA= NCELL_AB(ITYP-1)+1
          ICB= NCELL_AB(ITYP)
          DO I=ICA,ICB
           SP          = SP          + OUT_N0(I,1)*VOLCEL_N0(I)
           SE_SUM_N0_E = SE_SUM_N0_E + OUT_N0(I,2)*VOLCEL_N0(I)
           SE_SUM_N0_I = SE_SUM_N0_I + OUT_N0(I,3)*VOLCEL_N0(I)
          ENDDO
          WRITE(6,'(2I6,5E12.4)') ICA,ICB,SP*FLUX,SE_SUM_N0_E*FLUX
     .,   SE_SUM_N0_I*FLUX,0.,0.                   
          ENDDO 
      ENDIF 
      IF(NPRS>1) CALL BROADCAST_SOURCE_V_PL()

      RETURN
      END SUBROUTINE RUN_NEUTRAL

C ===== SOURCE: plotdummy.f
        subroutine GRSTRT
        integer K
        character(len=*) :: name
        return
        entry  GREND
        return
        entry  CHCTRC
        return
        entry  PLT2D 
        return
        entry  RPSOUT 
        return
        entry  GRNXTB(K,name)
        return
        entry  PLTEIR 
        return
        entry  S15AEF
        return
        entry  GRMRKS   
        return
        entry  GRNWPN   
        return
        entry  GRJMPS   
        return
        entry  GRJMP    
        return
        entry  GRDRW    
        return
        entry  GRNXTF   
        return
        entry  GRSCLC   
        return
        entry  GRSPTS   
        return
        entry  GRDSH    
        return
        entry  GRSCLV   
        return
        entry  GRCHRC   
        return
        entry  GRTXT    
        return
        entry  GRLN     
        return
        entry  GRFILL   
        return
        entry  GRARRW   
        return
        entry  GRAXS    
        return
        entry  GRTXTC   
        return
        entry  GR3DIM   
        return
        entry  GR3EXT   
        return
        entry  GR3AXS   
        return
        entry  GR3ROT   
        return
        entry  GR3PLO   
        return
        entry  GR3NET   
        return
        entry  GR3NT1   
        return
        entry  KURVEF   
        return
        entry  GRBLD    
        return
        entry  S13AAF   
        return
        end
                   
   
C ===== SOURCE: samp_usr.f
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
C ===== SOURCE: trace_n0.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: TIMUSR                                                C
C FUNCTION: follow a particle untill it leaves the cell IRGEN     C
C     Y. FENG                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TIMUSR(IRGEN,X0,Y0,Z0,VELX,VELY,VELZ,NEW,
     .                  ICNXT,TIMET,ICOS,IER,NPANU,SURF_P)
      USE GEOMETRY_PL 
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE 

      INTEGER ::  IRGEN,NEW,ICNXT,ICOS,IER,NPANU
      REAL*8  ::  X0,Y0,Z0,VELX,VELY,VELZ,TIMET
      LOGICAL ::  SURF_P
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 1.   description of the variables                                  C
C Input:                                                             C
C    IRGEN: initial cell number                                      C
C X0,Y0,Z0: initial position                                         C
C VELX,Y,Z: velocity vector                                          C
C      NEW: = 0, new particle                                        C
C          <=>0, following particle, initially on a surface          C
C   SURF_P: .TRUE. a surface particle                                C
C    NPANU: number of the MC particle                                C
C Output:                                                            C
C    TIMET: time for the particle until it leaves this cell          C
C           (also the travelling distance upto now because VL=1)     C
C   ICNXT : next cell number                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER,SAVE :: IERR,IRNXT,IPNXT,ITNXT
      INTEGER      :: ISR1,ISTS,I,K
      REAL*8       :: XNORM,YNORM,ZNORM
      REAL*8       :: SCOS,VL
      LOGICAL,SAVE :: DIAGNO=.FALSE.
c     LOGICAL,SAVE :: DIAGNO=.TRUE. 
      LOGICAL      :: OUT_OF_RANGE
C--------------------------------------------------------------------
C begin
      IER = 0
      IERR = 0
c     DIAGNO = NPANU== 139
c     if(NPANU==2) stop

      IF( OUT_OF_RANGE(IC0_N0,1,NCELL_N2) ) THEN
        IF(DIAGNO) CALL WRMESS('IC0.LE.0 .OR. IC0.GT.NCELL_N2')
        IERR = 1
        TIMET=-1.
        WRITE(6,*)'*** ERROR IC0_N0 OUT OF RANGE'             
        WRITE(6,'(4I5)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0
        RETURN
      ENDIF
C New particle
      IF(NEW.EQ.0) THEN
        IRGEN = IC0_N0
        TIMET = 0.

C replace the initial coordinates and velocity
        XP0 = X0
        YP0 = Y0
        ZP0 = Z0

        VX = VELX
        VY = VELY
        VZ = VELZ

C the new particle located just inside a cell
        IF(.NOT.SURF_P) THEN
          ISRF = 0
          ITRIA=0
        ENDIF
C maximal travellin distance to a plate
        TLMAX = 1.D30
        IF(NTRIANG.GT.0) CALL TLMAX_TRAVEL

        XFIN = X0
        YFIN = Y0
        ZFIN = Z0

        IF(DIAGNO) THEN
         WRITE(6,*)'==>TIMUSR FOLLOWS A NEW PARTICLE'
         WRITE(6,*)'==>NPANU:',NPANU                   
         WRITE(6,*)'ISRF  NZ0   JR0  JP0  JT0   IRGEN'
         WRITE(6,'(5I5,I8)')ISRF,NZ0_N0,JR0_N0,JP0_N0,JT0_N0,IRGEN
         WRITE(6,*)'INITIAL COORDINATS FROM EIRENE P0 AND V0:'
         WRITE(6,'(1P,6E12.4)')X0,Y0,Z0,VELX,VELY,VELZ
         WRITE(6,*)'Maximal distance:',TLMAX
         IF(ITRIA.NE.0)  THEN 
           WRITE(6,*)'To the Plate:',IPLATE,' Element:',ITRIA
         ENDIF  
         WRITE(6,'(A53)')
     &   '  IZ0  JR0  JP0  JT0  M_SF IS0-->IS1    TIMET     IC0'
        ENDIF
        NEW = 1
      ELSE
C Old particle:
C the particle is still located in the last cell, on the surface ISRF.
C the particle entries the cell IRGEN
        CALL INTO_CELL()
        IF(IERR.NE.0) THEN
         IF(DIAGNO) CALL WRMESS('ERROR IN INTO_CELL')
         TIMET=-1.
         WRITE(6,*)'*** ERROR FROM INTO_CELL WITH IERR=',IERR
         WRITE(6,'(4I5)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0
         RETURN
        ENDIF
        IRGEN = IC0_N0
      ENDIF
c
c 4. follow particle
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c 5. determine the intersection point
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0: CELL ADDRESS
C        XFIN,YFIN,ZFIN   : INITIAL POSITION
C        VX,VY,VZ         : VELOCITY
C        ISRF             : INITILALLY LOCATED SURFACE NUMBER
C        JDRJ,JDPJ,JDTJ   : BOUNDARY SURFACE SEQUENCE
c  definition of the points in the corners of the cell:
       CALL CORNER_POINTS()

c initial surface number
       ISR1 = ISRF
       CALL INTERSECT()
COUTPUT: XFIN,YFIN,ZFIN   : NEW POSITION
C        ISRF             : SURFACE NUMBER OF INTERSECTION
C        TFIN             : TRAVELED DISTANCE
C        IERR >0          : ERROR
       IF(IERR.NE.0) THEN
        IF(DIAGNO) THEN 
          CALL WRMESS('ERROR IN INTERSECT')
          WRITE(6,*)'IERR=',IERR
          WRITE(6,'(7I5,E12.4,I6)')
     +    NZ0_N0,JR0_N0,JP0_N0,JT0_N0,MSURF,ISR1,ISRF,TIMET,IRGEN
        ENDIF
        TIMET=-1.
        WRITE(6,*)'*** ERROR FROM INTERSECT WITH IERR=',IERR
        WRITE(6,'(4I5)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0
        RETURN
       ENDIF

      TIMET  = TFIN + TIMET
C Intersect an additional surface
      IF(TIMET.GT.TLMAX) THEN
         XFIN = XP0 + VX*TLMAX
         YFIN = YP0 + VY*TLMAX
         ZFIN = ZP0 + VZ*TLMAX

         TIMET = TLMAX
         IF(DIAGNO) THEN
          WRITE(6,*)'The Particle reaches a limiter'
          WRITE(6,*)'Distance traveled till now:',TIMET
          WRITE(6,'(A15,3f10.4)')'Final position:',XFIN,YFIN,ZFIN
         ENDIF
         ICNXT=ISR_TYPE(IPLATE)
         ICOS=1
         ISRF = 0
         RETURN
      ENDIF

      IRNXT =-I_JUMP(ISRF)*IDRJ
      IPNXT =-J_JUMP(ISRF)*IDPJ
      ITNXT =-K_JUMP(ISRF)*IDTJ

      CALL CHECK_SF_N0()
      IF(DIAGNO) WRITE(6,'(7I5,E12.4,I6)')
     +NZ0_N0,JR0_N0,JP0_N0,JT0_N0,MSURF,ISR1,ISRF,TIMET,IRGEN

      IF    (MSURF == 0) THEN
C NORMAOL
        ICNXT = IRGEN + 0
      ELSEIF(MSURF <  0) THEN
C Nontransparent surface
        ICNXT = MSURF
        ICOS=1
        ITRIA = 0
        
        IF(DIAGNO) THEN
          IF    (IRNXT /= 0) THEN 
            WRITE(6,'(A31,I4)')
     &      'A NON-TRANSP. R-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JR0_N0+INSRF(IRNXT)
          ELSEIF(IPNXT /= 0) THEN
            WRITE(6,'(A31,I4)')
     &      'A NON-TRANSP. P-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JP0_N0+INSRF(IPNXT)
          ELSE
            WRITE(6,'(A31,I4)')
     &      'A NON-TRANSP. T-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JT0_N0+INSRF(ITNXT)
          ENDIF
        ENDIF

      ELSE
C Non-default surfaces 
        ICNXT = -1
        ICOS=1
        ITRIA = 0

        IF(DIAGNO) THEN
          IF    (IRNXT /= 0) THEN 
            WRITE(6,'(A31,I4)')
     &      'A NON-DEFAULT R-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JR0_N0+INSRF(IRNXT)
          ELSEIF(IPNXT /= 0) THEN
            WRITE(6,'(A31,I4)')
     &      'A NON-DEFAULT P-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JP0_N0+INSRF(IPNXT)
          ELSE
            WRITE(6,'(A31,I4)')
     &      'A NON-DEFAULT T-SURF. WITH IND=',MSURF
            WRITE(6,'(A11,I4)')'SURFACE NR:',JT0_N0+INSRF(ITNXT)
          ENDIF
        ENDIF
      ENDIF

      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      NAME: NORUSR                                                    C
C  FUNCTION: DETERMINE THE NORMAL VECTOR OF THE REFLECTING SURFACE     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ENTRY NORUSR(ISTS,X0,Y0,Z0,XNORM,YNORM,ZNORM,SCOS,
     .             VELX,VELY,VELZ,IRGEN)
cbp: here: ists, scos not used
C Additional surface ?
      IF(ITRIA /= 0) THEN
        XNORM = V_TRIA_X(ITRIA)
        YNORM = V_TRIA_Y(ITRIA)
        ZNORM = V_TRIA_Z(ITRIA)
C CALLED BY SAMUSR?
        IF(NSR_N0.NE.0) THEN
         XNORM= XNORM*NSR_N0
         YNORM= YNORM*NSR_N0
         ZNORM= ZNORM*NSR_N0

         NSR_N0 = 0
        ELSEIF(XNORM*VX+YNORM*VY+ZNORM*VZ.LT.0.) THEN
         XNORM=-XNORM
         YNORM=-YNORM
         ZNORM=-ZNORM
        ENDIF

        IF(DIAGNO) write(6,*)VELX,VELY,VELZ,XNORM,YNORM,ZNORM
      ELSEIF(MSURF < 0) THEN
C REFLECTING SURFACE ---------------------------------------------
c surface number:
         i=isrf
         do 2 k=1,3
            x3eck(k) = xmy(ipunkt(i,k))
            y3eck(k) = ymy(ipunkt(i,k))
            z3eck(k) = zmy(ipunkt(i,k))
 2       continue
c
C NORMAL VECTOR
         XNORM =(Y3ECK(3)-Y3ECK(1))*(Z3ECK(2)-Z3ECK(1)) -
     +          (Z3ECK(3)-Z3ECK(1))*(Y3ECK(2)-Y3ECK(1))
         YNORM =(Z3ECK(3)-Z3ECK(1))*(X3ECK(2)-X3ECK(1)) -
     +          (X3ECK(3)-X3ECK(1))*(Z3ECK(2)-Z3ECK(1))
         ZNORM =(X3ECK(3)-X3ECK(1))*(Y3ECK(2)-Y3ECK(1)) -
     +          (Y3ECK(3)-Y3ECK(1))*(X3ECK(2)-X3ECK(1))

         VL = DSQRT(XNORM**2+YNORM**2+ZNORM**2)
         XNORM=XNORM/VL
         YNORM=YNORM/VL
         ZNORM=ZNORM/VL

         IF(XNORM*VX+YNORM*VY+ZNORM*VZ.LT.0.) THEN
         XNORM=-XNORM
         YNORM=-YNORM
         ZNORM=-ZNORM
         ENDIF

         IF(DIAGNO) THEN
          WRITE(6,*)'===> NORUSR CALCULATES THE NORMAOL VECTOR'
          WRITE(6,*)'COS=',XNORM*VELX+YNORM*VELY+ZNORM*VELZ
         ENDIF

      ELSEIF( MSURF >  0) THEN
C Periodic surface, leading to a new particle
C Depending on MSURF:

         CALL INTO_CELL()
         IRGEN = IC0_N0
         IF(IERR.NE.0) THEN
            IF(DIAGNO) CALL WRMESS('ERROR IN INTO_CELL')
            TIMET=-1.
            RETURN
         ENDIF

         VELX = VX
         VELY = VY
         VELZ = VZ

         X0   = XFIN
         Y0   = YFIN
         Z0   = ZFIN
      ENDIF

      RETURN

      CONTAINS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CORNER_POINTS()

       IMPLICIT NONE

C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0: FINE CELL ADDRESS
C        JDRJ,JDPJ,JDTJ : CORNER POINTS SELECTION SEQUENCE
C
COUTPUT: XMY,YMY,ZMY(8): X- Y- Z-COORDINATES OF CORNER POINTS
C-----------------------------------------------------------------
C      1     I+N_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      2     I+N_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      3     I+M_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      4     I+M_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+N_JUMP(IDTJ)
C      5     I+N_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C      6     I+N_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C      7     I+M_JUMP(IDRJ)    J+M_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C      8     I+M_JUMP(IDRJ)    J+N_JUMP(IDPJ)   K+M_JUMP(IDTJ)
C------------------------------------------------------------------
C  where ID*J=-1,1, indicating the points sequence.
c  N_JUMP(L) = 0,1 and M_JUMP=1,0 when L=-1,1.
c  Standard case: IDRJ=IDPJ=IDTJ= -1.
c
c  The cell is surrounded by 12 triangles defined IPUNKT(I,J)
c  I: surface number (i=1,12)
c  J: point number   (J=1,3)
C
      INTEGER,DIMENSION(8) :: IG_P
      INTEGER :: I,I1,I2,J1,J2,K1,K2,ISTEPJ,ISTEPK,IG_BKGRD



      ISTEPJ=                  SRF_RADI(NZ0_N0)
      ISTEPK= SRF_POLO(NZ0_N0)*SRF_RADI(NZ0_N0)  
c two toroidal planes with M1 and M2
      I1 = JR0_N0 + N_JUMP(IDRJ)
      I2 = JR0_N0 + M_JUMP(IDRJ)

      J1 = JP0_N0 + N_JUMP(IDPJ)
      J2 = JP0_N0 + M_JUMP(IDPJ)

      K1 = JT0_N0 + N_JUMP(IDTJ)
      K2 = JT0_N0 + M_JUMP(IDTJ)

      IG_BKGRD= K1*ISTEPK + GRID_P_OS(NZ0_N0)
      IG_P(1) = I1 + J1*ISTEPJ + IG_BKGRD
      IG_P(2) = I1 + J2*ISTEPJ + IG_BKGRD
      IG_P(3) = I2 + J2*ISTEPJ + IG_BKGRD
      IG_P(4) = I2 + J1*ISTEPJ + IG_BKGRD

      IG_BKGRD= K2*ISTEPK + GRID_P_OS(NZ0_N0)
      IG_P(5) = I1 + J1*ISTEPJ + IG_BKGRD
      IG_P(6) = I1 + J2*ISTEPJ + IG_BKGRD
      IG_P(7) = I2 + J2*ISTEPJ + IG_BKGRD
      IG_P(8) = I2 + J1*ISTEPJ + IG_BKGRD

      DO I=1,4
         XMY(I) = RG(IG_P(I))*COSPHI(K1+PHI_PL_OS(NZ0_N0))
         YMY(I) = RG(IG_P(I))*SINPHI(K1+PHI_PL_OS(NZ0_N0))
         ZMY(I) = ZG(IG_P(I))
      ENDDO
      DO I=5,8
         XMY(I) = RG(IG_P(I))*COSPHI(K2+PHI_PL_OS(NZ0_N0))
         YMY(I) = RG(IG_P(I))*SINPHI(K2+PHI_PL_OS(NZ0_N0))
         ZMY(I) = ZG(IG_P(I))
      ENDDO
      RETURN
      END SUBROUTINE CORNER_POINTS

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: INTERSECT                                            C
C FUNCTION: DETERMINE THE INTERSECTION POINT OF FLIGHT OF PARTICEC
C           WITH A CELL BOUNDARY                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTERSECT()

      IMPLICIT NONE              
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0: CELL ADDRESS
C INPUT: XFIN,YFIN,ZFIN : INITIAL POSITION
C INPUT: VX,VY,VZ       : VELOCITY
C INPUT: ISRF           : INITILALLY LOCATED SURFACE NUMBER
C INPUT: JDRJ,JDPJ,JDTJ : BOUNDARY SURFACE SEQUENCE

COUTPUT: XFIN,YFIN,ZFIN : NEW POSITION
C        ISRF           : SURFACE NUMBER OF INTERSECTION
C        TFIN           : TRAVELED DISTANCE
C        IERR >0        : ERROR

C  parameters used in this program only
      REAL*8, PARAMETER :: TOLER=1.E-12
      REAL*8, DIMENSION(12) :: 
     R        xy,       xz,       yz 
      INTEGER,DIMENSION(12) ::
     I        J1,       J2,       J3
      INTEGER :: NMIN,N_CROSS,ICHECK,indmin,MINI,IN
      REAL*8  :: zeitmn,XCR,YCR,ZCR,CHECK_1,CHECK_2
C--------------------------------------------------------------------
      ierr = 0

c side surfaces:
      DO I=1,12
        J1(I)=ipunkt(i,1)
        J2(I)=ipunkt(i,2)
        J3(I)=ipunkt(i,3)
      ENDDO

      DO i=1,12

        xy(i) = (xmy(j2(i))-xmy(j1(i)))*(ymy(j3(i))-ymy(j1(i)))
     +-         (ymy(j2(i))-ymy(j1(i)))*(xmy(j3(i))-xmy(j1(i)))

        xz(i) = (zmy(j2(i))-zmy(j1(i)))*(xmy(j3(i))-xmy(j1(i)))
     +-         (xmy(j2(i))-xmy(j1(i)))*(zmy(j3(i))-zmy(j1(i)))

        yz(i) = (ymy(j2(i))-ymy(j1(i)))*(zmy(j3(i))-zmy(j1(i)))
     +-         (zmy(j2(i))-zmy(j1(i)))*(ymy(j3(i))-ymy(j1(i)))
      ENDDO
C
c DISTAN: Distance of the point (xfin,yfin,zfin) to the plane
c SINFI = SIN(FI): FI= the incidence angle
      DO i=1,12
       DISTAN(I)=- (xfin-xmy(j1(I)))*yz(I)
     &           - (yfin-ymy(j1(I)))*xz(I)
     &           - (zfin-zmy(j1(I)))*xy(I)
       SINFI(I) = vx*yz(I) + vy*xz(I) + vz*xy(I)
      ENDDO

c IF DABS(SINFI).LT.TOLER: either the plane does not exist or
c                    v parallel to the plane
      nmin=0

C     IF(DIAGNO)write(6,'(6E12.4)')SINFI(1:12),DISTAN(1:12)
C     IF(DIAGNO)write(6,'(3E12.4)')VX,VY,VZ 
      IF(ISRF/=0) DISTAN(ISRF) = 0.
      IF( ISRF >= 9) THEN
        IF    (ISRF == 9) THEN
          DISTAN(10) = 0.
        ELSEIF(ISRF ==10) THEN
          DISTAN( 9) = 0.
        ELSEIF(ISRF ==11) THEN
          DISTAN(12) = 0.
        ELSE
          DISTAN(11) = 0.
        ENDIF   
      ENDIF

      DO I=1,12
       IF(SINFI(I)*DISTAN(I) > 0.) THEN
        TCROSS(I) = DISTAN(I)/SINFI(I)
        nmin = nmin + 1
        iebene(nmin) = i
       ENDIF
      ENDDO
c=================================================================
c define the smallest time:
c-----------------------------------------------------------------
c check whether the point is "real" or "virtual":
C TOTAL NUMBER OF THE INTERSECTION POINTS: N_CROSS = nmin
      N_CROSS = nmin
      IF(N_CROSS.eq. 0) THEN
        IERR = 1
        RETURN
      ENDIF

      DO 100 ICHECK=1,N_CROSS

      indmin = iebene(1)
      zeitmn = tcross(indmin)
      mini   = 1

      do i=2,nmin
         in = iebene(i)
         if(tcross(in).lt.zeitmn) then
           zeitmn = tcross(in)
           indmin = in
           mini   = i
         endif
      enddo

      xcr = xfin + vx*zeitmn
      ycr = yfin + vy*zeitmn
      zcr = zfin + vz*zeitmn
c define a triangle:
      do i=1,3
      x3eck(i) = xmy(ipunkt(indmin,i)) - xcr
      y3eck(i) = ymy(ipunkt(indmin,i)) - ycr
      z3eck(i) = zmy(ipunkt(indmin,i)) - zcr
      enddo
c
      vpr1(1) = y3eck(1)*z3eck(2) - z3eck(1)*y3eck(2)
      vpr1(2) = z3eck(1)*x3eck(2) - x3eck(1)*z3eck(2)
      vpr1(3) = x3eck(1)*y3eck(2) - y3eck(1)*x3eck(2)
c
      vpr2(1) = y3eck(2)*z3eck(3) - y3eck(3)*z3eck(2)
      vpr2(2) = z3eck(2)*x3eck(3) - z3eck(3)*x3eck(2)
      vpr2(3) = x3eck(2)*y3eck(3) - x3eck(3)*y3eck(2)
c
c first check:
      CHECK_1=vpr1(1)*vpr2(1)+vpr1(2)*vpr2(2)+vpr1(3)*vpr2(3)
c     IF(DIAGNO) WRITE(6,'(A8,E12.4)')'CHECK_1=',CHECK_1
      IF(CHECK_1.GT.0.) THEN
c
        vpr3(1) = y3eck(3)*z3eck(1) - z3eck(3)*y3eck(1)
        vpr3(2) = z3eck(3)*x3eck(1) - x3eck(3)*z3eck(1)
        vpr3(3) = x3eck(3)*y3eck(1) - y3eck(3)*x3eck(1)
c
c last check:
        CHECK_2=vpr3(1)*vpr2(1)+vpr3(2)*vpr2(2)+vpr3(3)*vpr2(3)
c       IF(DIAGNO )WRITE(6,'(A8,E12.4)')'CHECK_2=',CHECK_2
        IF(CHECK_2.GT.0.) THEN

         XFIN = XCR
         YFIN = YCR
         ZFIN = ZCR 
         TFIN = ZEITMN
         ISRF = INDMIN
         RETURN
        ENDIF
      ENDIF

      nmin=nmin-1
      do i=mini,nmin
        iebene(I) = iebene(I+1)
      enddo

100   CONTINUE
      IERR = 2
      RETURN
      END SUBROUTINE INTERSECT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: CHECK_SF_N0                                          C
C FUNCTION: CHECK THE BOUNDARY SURFACE FOR NEUTRAL               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CHECK_SF_N0()
       USE SURFACE_PL

       IMPLICIT NONE             

       INTEGER :: ISURF
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0 : CELL ADDRESS
C        IRNXT,IPNXT,ITNXT: R,P,T JUMP INDEX

COUTPUT: MSURF          =0: NORMAL SURAFCE, PARTICLE CONTINUES
C                       <0: NON-TRANSPARENT SURFACE, |MSURF| SURFACE NUMBER
C-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C R1: Radial surface
      IF(IRNXT.NE.0) THEN
        ISURF= JR0_N0 + INSRF(IRNXT) 
     .+       (JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))*SRF_RADI(NZ0_N0) 
     .+        NRS_OFF(NZ0_N0)
        MSURF = IDSURR(ISURF)
C P1: Poloidal surface
      ELSEIF(IPNXT.NE.0) THEN
        ISURF= JR0_N0 + (JP0_N0+INSRF(IPNXT)
     .+        JT0_N0*SRF_POLO(NZ0_N0))*ZON_RADI(NZ0_N0) 
     .+        NPS_OFF(NZ0_N0)

        MSURF = IDSURP(ISURF)
C T1: Toroidal surface
      ELSEIF(ITNXT.NE.0) THEN
        ISURF= JR0_N0 + (JP0_N0+(JT0_N0+INSRF(ITNXT))
     .*        ZON_POLO(NZ0_N0))*ZON_RADI(NZ0_N0) 
     .+        NTS_OFF(NZ0_N0)

        MSURF = IDSURT(ISURF)
      ENDIF
      RETURN
      END SUBROUTINE CHECK_SF_N0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: INTO_CELL                                            C
C FUNCTION: DETERMINE THE CELL NUMBER INTO WHICH A PARTICLE      C
C           TRAVELS. IF NECESSARY, THE FLIGHT CAN BE CHANGED,    C
C           POSITION CAN BE SHIFTED ACCORDING TO THE SYMMETY OF  C
C           THE GIVEN GEOMETRY                                   C
C           THIS SUBROUTINE IS STANDARD FOR W7-AS AND -X CONF.   C
C           IN OTHER CASES, IT MUST BY REPLACED BY THE USER.     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE INTO_CELL()
       USE MAPPING_TOROIDAL

       IMPLICIT NONE               
C INPUT: JR0_N0,JP0_N0,JT0_N0,NZ0_N0  : CELL ADDRESS
C        XFIN,YFIN,ZFIN   : POSITION
C        VX,VY,VZ         : VELOCITY
C        MSURF            : Surface Type (see CHECK_SF_N0)
C        ISRF             : LOCATED SURFACE NUMBER
C        JDRJ,JDPJ,JDTJ   : BOUNDARY SURFACE SEQUENCE
C        IRNXT,IPNXT,ITNXT: R,P,T JUMP INDEX

COUTPUT: XFIN,YFIN,ZFIN   : NEW POSITION
C        ISRF             : NEW SURFACE NUMBER
C        JR0,JP0,JT0,NZ0  : NEW CELL ADDRESS
C        VX,VY,VZ         : NEW VELOCITY
C        JDRJ,JDPJ,JDTJ   : NEW BOUNDARY SURFACE SEQUENCE
C        IERR >0          : ERROR
C                           THE SYMMETRY OF THE GIVEN GEOMETRY
       INTEGER :: IG,IS,MAP_SF,PAI_SF,TSURF,IPHI
       REAL*8  :: R0,Z0,V_R,V_PHI,RJ,PJ,TJ
C-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      IERR = 0
C Case0: Normal surface
      IF    (MSURF.EQ.0) THEN
        JR0_N0 = JR0_N0 + IRNXT
        JP0_N0 = JP0_N0 + IPNXT
        JT0_N0 = JT0_N0 + ITNXT
        IG     = JR0_N0 + (JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*          ZON_RADI(NZ0_N0) + MESH_P_OS(NZ0_N0)
        IC0_N0 = IDCELL(IG)
        IF( OUT_OF_RANGE(IC0_N0,1,NCELL_N2) ) IERR = 1
        
        ISRF  = NEWSRF(ISRF)
        RETURN
      ENDIF 
         
      IF    ( IRNXT /=0 ) THEN
        CALL R_SF_JUMP(MSURF,IRNXT
     .,      NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )
        ISRF  = NEWSRF(ISRF)
      ELSEIF( IPNXT /=0 ) THEN
        CALL P_SF_JUMP(MSURF,IPNXT
     .,      NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )
        ISRF  = NEWSRF(ISRF)
      ELSEIF( ITNXT /=0 ) THEN
        TSURF = JT0_N0 + INSRF(ITNXT)

        IF    ( MSURF == 1 ) THEN 
C Periodic
          R0    = SQRT(XFIN**2+YFIN**2)
          IPHI  = TSURF  + PHI_PL_OS(NZ0_N0)
          V_R   = VX*COSPHI(IPHI) + VY*SINPHI(IPHI)
          V_PHI =-VX*SINPHI(IPHI) + VY*COSPHI(IPHI)

          IF(DIAGNO) THEN
           WRITE(6,*)'**>Periodic toroidal surface:'
           WRITE(6,*)'   IPHI  XFIN  YFIN  ZFIN  VX  VY  VZ'
           WRITE(6,'(A4,I4,6F9.3)')'old:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF


          CALL T_SF_JUMP(MSURF,ITNXT
     .,        NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )

          TSURF = JT0_N0 + INSRF(-ITNXT)
          IPHI  = TSURF  + PHI_PL_OS(NZ0_N0)

          XFIN  = R0*COSPHI(IPHI)
          YFIN  = R0*SINPHI(IPHI)

          VX   = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
          VY   = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)

          IF(DIAGNO) THEN
           WRITE(6,'(A4,I4,6F9.3)')'new:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF
          ISRF  = NEWSRF(ISRF)
        ELSEIF( MSURF == 2 ) THEN
C Up/down symmetric surface
          IPHI  = TSURF  + PHI_PL_OS(NZ0_N0)
          V_R   = VX*COSPHI(IPHI) + VY*SINPHI(IPHI)
          V_PHI =-VX*SINPHI(IPHI) + VY*COSPHI(IPHI)

          IF(DIAGNO) THEN
           WRITE(6,*)'**>Asymmetric toroidal surface:'
           WRITE(6,*)'   IPHI  XFIN  YFIN  ZFIN  VX  VY  VZ'
           WRITE(6,'(A4,I4,6F9.3)')'old:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF
          ZFIN  =-ZFIN
          V_PHI =-V_PHI

          VX    = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
          VY    = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)
          VZ    =-VZ

          IF(DIAGNO) THEN
           WRITE(6,'(A4,I4,6F9.3)')'new:'
     .,                  IPHI,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF

          CALL T_SF_JUMP(MSURF,ITNXT
     .,        NZ0_N0,JR0_N0,JP0_N0,JT0_N0,0.,0.,0. )

        ELSEIF( MSURF == 3 ) THEN 
C Mapping 
          MAP_SF = 0
          DO IS=1,TOTAL_MAP_SF_T
            IF(  NZ0_N0 == ZONE_NR_MAP_T(IS) .AND. 
     .           TSURF  == T_SF_NR_MAP_T(IS) ) THEN 
              MAP_SF = IS
              EXIT
            ENDIF
          ENDDO
          IF( MAP_SF == 0 ) THEN
            WRITE(6,*)'Mapping surface not found'
            CALL STOP_ALL('Stopped in INTO_CELL')
          ENDIF 

          IF(DIAGNO) THEN
           WRITE(6,*)'**>Mapping toroidal surface:'
           WRITE(6,*)'MAP NR.  XFIN  YFIN  ZFIN  VX  VY  VZ'
           WRITE(6,'(A4,I4,6F9.3)')'old:'
     .,                MAP_SF,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF

          PAI_SF = MAP_PAIR_SF_T(MAP_SF)
          JT0_N0 = T_SF_NR_MAP_T(PAI_SF)

          IPHI  = TSURF + PHI_PL_OS(NZ0_N0)
          V_R   = VX*COSPHI(IPHI) + VY*SINPHI(IPHI)
          V_PHI =-VX*SINPHI(IPHI) + VY*COSPHI(IPHI)

          R0 = SQRT(XFIN**2+YFIN**2)
          Z0 = ZFIN
          CALL PERFORM_MAPPING_TOROIDAL_N0
     .    (MAP_SF,NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,R0,Z0,IERR)

          IF(IERR > 0 .OR. IERR==-3) THEN  
            IERR = 2
            RETURN
          ELSEIF(IERR < 0) THEN
            RJ = -0.99999
            IF(IERR == -2) RJ=-RJ
            TJ = JT0_N0
            IPHI  = JT0_N0 + PHI_PL_OS(NZ0_N0)

            CALL RZ_REAL_COORDINATES
     .          (NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,TJ,R0,Z0)

            XFIN  = R0*COSPHI(IPHI)
            YFIN  = R0*SINPHI(IPHI)
            ZFIN  = Z0
            IERR  = 0
          ENDIF

C Mapping on itself
          IF( PAI_SF == MAP_SF ) THEN
            ZFIN  =-ZFIN
            V_PHI =-V_PHI

            VX    = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
            VY    = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)
            VZ    =-VZ
C Mapping on the last toroidal surfaces
          ELSEIF((JT0_N0 == 0                .AND. NZ0_N0 ==  0 ) .OR.
     .           (JT0_N0 == ZON_TORO(NZ0_N0) .AND. NZ0_N0 == NZONET-1)
     .          )THEN
            IPHI  = JT0_N0 + PHI_PL_OS(NZ0_N0)

            XFIN  = R0*COSPHI(IPHI)
            YFIN  = R0*SINPHI(IPHI)

            VX   = V_R*COSPHI(IPHI) - V_PHI*SINPHI(IPHI)
            VY   = V_R*SINPHI(IPHI) + V_PHI*COSPHI(IPHI)
            ISRF  = NEWSRF(ISRF)
          ELSE  
            ISRF  = NEWSRF(ISRF)
          ENDIF
          IF(JT0_N0 == ZON_TORO(NZ0_N0)) JT0_N0 = JT0_N0 - 1

          IF(DIAGNO) THEN
           WRITE(6,'(A4,I4,6F9.3)')'new:'
     .,                PAI_SF,XFIN,YFIN,ZFIN,VX,VY,VZ
          ENDIF
        ENDIF 

      ENDIF

      IG     = JR0_N0 + (JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*        ZON_RADI(NZ0_N0) + MESH_P_OS(NZ0_N0)
      IC0_N0 = IDCELL(IG)
      IF( OUT_OF_RANGE(IC0_N0,1,NCELL_N2) ) IERR = 1

      RETURN
      END SUBROUTINE INTO_CELL

      END SUBROUTINE TIMUSR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     NAME: TLMAX_TRAVEL                                         C
C FUNCTION: DETERMINE THE the travel distance of a particle to   C
C           a plate                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TLMAX_TRAVEL
      USE NEUTRAL_TRANSPORT 
      
      IMPLICIT NONE            
C INPUT: XP0,YP0,ZP0 : Particle location
C        VX,VY,VZ    : Velocity
      INTEGER :: IP,N_PL,IPM,MI,NMIN,N_CROSS,INDMIN,MINI,IN,ICHECK
     I,          I,J,L,LL,N,NY,N_YA,NY_B,L1,L2,IPOSI,ITR
      REAL*8  :: T_MIN,T_MI,T_MAX,ZEITMN,CHECK_1,CHECK_2,XCR,YCR,ZCR
     R,          T1,T2,DX1,DX2,DY1,DY2,DZ1,DZ2,XR1,XR2,YR1,YR2,ZR1,ZR2 

C Defaul value
      TLMAX = 1.D30

      IF(NTRIANG.LT.10) THEN
        CALL TLMAX_TRAVEL1
      ELSE
        CALL TLMAX_TRAVEL2
      ENDIF
      RETURN

      CONTAINS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TLMAX_TRAVEL1
      IMPLICIT NONE            
C INPUT: XP0,YP0,ZP0 : Particle location
C        VX,VY,VZ    : Velocity
C--------------------------------------------------------------------
C parameters used in this program only

C case1: if NTRIANG is mall, scan all the triangles
c DISTAN: Distance of the point (xfin,yfin,zfin) to the plane
c SINFI = SIN(FI): FI= the incidence angle
      DO I=1,NTRIANG
       DISTAN(I)=- (XP0-X_TRIA(1,I))*V_TRIA_X(I)
     &           - (YP0-Y_TRIA(1,I))*V_TRIA_Y(I)
     &           - (ZP0-Z_TRIA(1,I))*V_TRIA_Z(I)

       SINFI(I) =   VX*V_TRIA_X(I)
     &           +  VY*V_TRIA_Y(I)
     &           +  VZ*V_TRIA_Z(I)
      ENDDO

c IF SINFI=0.: either the plane does not exist or
c              v parallel to the plane
      nmin=0
      DO I=1,NTRIANG
       IF(SINFI(I).NE.0. .AND. I.NE.ITRIA) THEN
        TCROSS(I) = DISTAN(I)/SINFI(I)
        IF(TCROSS(I).GT.0.) THEN
         nmin = nmin + 1
         iebene(nmin) = i
        ENDIF
       ENDIF
      ENDDO
c=================================================================
c define the smallest time:
c-----------------------------------------------------------------
c check whether the point is "real" or "virtual":
C TOTAL NUMBER OF THE INTERSECTION POINTS: N_CROSS = nmin
      N_CROSS = nmin
      IPLATE= 0
      ITRIA = 0
      IF(N_CROSS == 0) RETURN

      CHECK_LOOP : DO ICHECK=1,N_CROSS

      indmin = iebene(1)
      zeitmn = tcross(indmin)
      mini   = 1

      do i=2,nmin
         in = iebene(i)
         if(tcross(in).lt.zeitmn) then
           zeitmn = tcross(in)
           indmin = in
           mini   = i
         endif
      enddo

      xcr = XP0 + vx*zeitmn
      ycr = YP0 + vy*zeitmn
      zcr = ZP0 + vz*zeitmn
c define a triangle:
      do i=1,3
      x3eck(i) = X_TRIA(I,indmin) - xcr
      y3eck(i) = Y_TRIA(I,indmin) - ycr
      z3eck(i) = Z_TRIA(I,indmin) - zcr
      enddo
c
      vpr1(1) = y3eck(1)*z3eck(2) - z3eck(1)*y3eck(2)
      vpr1(2) = z3eck(1)*x3eck(2) - x3eck(1)*z3eck(2)
      vpr1(3) = x3eck(1)*y3eck(2) - y3eck(1)*x3eck(2)
c
      vpr2(1) = y3eck(2)*z3eck(3) - y3eck(3)*z3eck(2)
      vpr2(2) = z3eck(2)*x3eck(3) - z3eck(3)*x3eck(2)
      vpr2(3) = x3eck(2)*y3eck(3) - x3eck(3)*y3eck(2)
c
c first check:
      CHECK_1=vpr1(1)*vpr2(1)+vpr1(2)*vpr2(2)+vpr1(3)*vpr2(3)
      IF(CHECK_1.GE.0.) THEN
c
        vpr3(1) = y3eck(3)*z3eck(1) - z3eck(3)*y3eck(1)
        vpr3(2) = z3eck(3)*x3eck(1) - x3eck(3)*z3eck(1)
        vpr3(3) = x3eck(3)*y3eck(1) - y3eck(3)*x3eck(1)
c
c last check:
        CHECK_2=vpr3(1)*vpr2(1)+vpr3(2)*vpr2(2)+vpr3(3)*vpr2(3)
        IF(CHECK_2.GE.0.) THEN

         TLMAX  = zeitmn
         ITRIA  = indmin
         IPLATE = NP_NUM(indmin)
         EXIT CHECK_LOOP
        ENDIF
      ENDIF

      nmin=nmin-1
      do i=mini,nmin
        iebene(I) = iebene(I+1)
      enddo

      ENDDO CHECK_LOOP 

      RETURN
      END SUBROUTINE TLMAX_TRAVEL1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TLMAX_TRAVEL2
      
      IMPLICIT NONE            
C INPUT: XP0,YP0,ZP0 : Particle location
C        VX,VY,VZ    : Velocity
C--------------------------------------------------------------------
C parameters used in this program only

      DO IP=1,NPSORT
       DISTAN(IP)=-(XP0-X_SORT(IP))*V_SORT_X(IP)
     &            -(YP0-Y_SORT(IP))*V_SORT_Y(IP)
     &            -(ZP0-Z_SORT(IP))*V_SORT_Z(IP)

       SINFI(IP) =  VX*V_SORT_X(IP)
     &           +  VY*V_SORT_Y(IP)
     &           +  VZ*V_SORT_Z(IP)
      ENDDO

c IF SINFI=0.: either the plane does not exist or
c              v parallel to the plane
      N_PL=0
      DO IP=1,NPSORT
       IF(SINFI(IP).NE.0.) THEN
        TCR_PL(1,IP)  = (DISTAN(IP)+D_MINI(IP))/SINFI(IP)
        TCR_PL(2,IP)  = (DISTAN(IP)+D_MAXI(IP))/SINFI(IP)
        IF(TCR_PL(1,IP) > 0. .OR. TCR_PL(2,IP) > 0.) THEN
         N_PL         = N_PL + 1
         IEB_PL(N_PL) = IP
        ENDIF
       ENDIF
      ENDDO
c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*c*
c     write(6,*)'N_PL=',N_PL
      IF(N_PL.EQ.0) THEN 
        IPLATE= 0
        ITRIA = 0
        RETURN 
      ENDIF
C smallst distance
250   IPM   = IEB_PL(1)
      T_MIN = MIN(TCR_PL(1,IPM),TCR_PL(2,IPM))
      MI    = 1

      DO I=2,N_PL
         IP = IEB_PL(I)
         T_MI = MIN(TCR_PL(1,IP),TCR_PL(2,IP))
         IF(T_MI.LT.T_MIN) THEN
           T_MIN = T_MI
           IPM   = IP
           MI    = I
         ENDIF
      ENDDO

c     write(6,*) MI,IPM,T_MIN

      N_YA  = 1
      IEB_YA(1) = MI

      T_MAX     = MAX(TCR_PL(1,IPM),TCR_PL(2,IPM))
      DO I=1,N_PL
        IF(I.NE.MI) THEN
         IP = IEB_PL(I)
         T_MI = MIN(TCR_PL(1,IP),TCR_PL(2,IP))

         IF(T_MI.LE.T_MAX) THEN
           N_YA = N_YA + 1
           IEB_YA(N_YA) = I
         ENDIF
        ENDIF
      ENDDO

      LL = 0
      DO 300 N=1,N_YA
       IPOSI = IEB_PL(IEB_YA(N))
       T1 = TCR_PL(1,IPOSI)
       T2 = TCR_PL(2,IPOSI)

       DX1 =  VX*T1
       DY1 =  VY*T1
       DZ1 =  VZ*T1
       DX2 =  VX*T2
       DY2 =  VY*T2
       DZ2 =  VZ*T2

       XR1 = XP0 + MAX(DX1,DX2)
       YR1 = YP0 + MAX(DY1,DY2)
       ZR1 = ZP0 + MAX(DZ1,DZ2)
       XR2 = XP0 + MIN(DX1,DX2)
       YR2 = YP0 + MIN(DY1,DY2)
       ZR2 = ZP0 + MIN(DZ1,DZ2)

       IF( (XR1.LT.X_SMIN(IPOSI) .OR. XR2.GT.X_SMAX(IPOSI))
     & .OR.(YR1.LT.Y_SMIN(IPOSI) .OR. YR2.GT.Y_SMAX(IPOSI))
     & .OR.(ZR1.LT.Z_SMIN(IPOSI) .OR. ZR2.GT.Z_SMAX(IPOSI))
     &   ) GOTO 300

       L1 = LL + 1
       DO I=NPS_GRP(IPOSI-1)+1,NPS_GRP(IPOSI)
         LL        = LL + 1
         MTRIA(LL) = NPN_GRP(I)
       ENDDO
       L2 = LL
       NYBACK(L1:L2) = N

       DO I=L1,L2
         DISTAN(I)=- (XP0-X_TRIA(1,MTRIA(I)))*V_TRIA_X(MTRIA(I))
     &             - (YP0-Y_TRIA(1,MTRIA(I)))*V_TRIA_Y(MTRIA(I))
     &             - (ZP0-Z_TRIA(1,MTRIA(I)))*V_TRIA_Z(MTRIA(I))

          SINFI(I)=  VX*V_TRIA_X(MTRIA(I))
     &            +  VY*V_TRIA_Y(MTRIA(I))
     &            +  VZ*V_TRIA_Z(MTRIA(I))
       ENDDO
300   CONTINUE

      NY_B=0

      IF(LL.EQ.0) GOTO 900

c IF SINFI=0.: either the plane does not exist or
c              v parallel to the plane

       NMIN = 0
       DO I=1,LL
        IF(SINFI(I).NE.0. .AND. MTRIA(I).NE.ITRIA) THEN
        TCROSS(I) = DISTAN(I)/SINFI(I)
        IF(TCROSS(I).GT.0.) THEN
         NMIN         = NMIN + 1
         IEBENE(NMIN) = I
        ENDIF
       ENDIF
      ENDDO
      IF(NMIN.EQ.0) GOTO 900

260   indmin = iebene(1)
      zeitmn = tcross(indmin)
      mini   = 1

      do i=2,nmin
         in = iebene(i)
         if(tcross(in).lt.zeitmn) then
           zeitmn = tcross(in)
           indmin = in
           mini   = i
         endif
      enddo
      ITR = MTRIA(indmin)

      xcr = XP0 + vx*zeitmn
      ycr = YP0 + vy*zeitmn
      zcr = ZP0 + vz*zeitmn
c define a triangle:
      do i=1,3
      x3eck(i) = X_TRIA(I,ITR) - xcr
      y3eck(i) = Y_TRIA(I,ITR) - ycr
      z3eck(i) = Z_TRIA(I,ITR) - zcr
      enddo
c
      vpr1(1) = y3eck(1)*z3eck(2) - z3eck(1)*y3eck(2)
      vpr1(2) = z3eck(1)*x3eck(2) - x3eck(1)*z3eck(2)
      vpr1(3) = x3eck(1)*y3eck(2) - y3eck(1)*x3eck(2)
c
      vpr2(1) = y3eck(2)*z3eck(3) - y3eck(3)*z3eck(2)
      vpr2(2) = z3eck(2)*x3eck(3) - z3eck(3)*x3eck(2)
      vpr2(3) = x3eck(2)*y3eck(3) - x3eck(3)*y3eck(2)
c
c first check:
      CHECK_1=vpr1(1)*vpr2(1)+vpr1(2)*vpr2(2)+vpr1(3)*vpr2(3)
      IF(CHECK_1.GE.0.) THEN
c
        vpr3(1) = y3eck(3)*z3eck(1) - z3eck(3)*y3eck(1)
        vpr3(2) = z3eck(3)*x3eck(1) - x3eck(3)*z3eck(1)
        vpr3(3) = x3eck(3)*y3eck(1) - y3eck(3)*x3eck(1)
c
c last check:
        CHECK_2=vpr3(1)*vpr2(1)+vpr3(2)*vpr2(2)+vpr3(3)*vpr2(3)
        IF(CHECK_2.GE.0.) THEN

         IF(zeitmn.GT.T_MAX) THEN
           NY_B = NY_B + 1
           NYSTOR(NY_B) = NYBACK(indmin)
         ELSE
           TLMAX  = zeitmn
           ITRIA  = ITR
           IPLATE = NP_NUM(ITRIA)
           return
         ENDIF
        ENDIF
      ENDIF
      nmin=nmin-1
      do i=mini,nmin
        iebene(I) = iebene(I+1)
      enddo
      IF(nmin.GT.0) GOTO 260

c900   IF(DIAGNO) THEN
c      WRITE(6,'(A12,I4,4X,16I3)')'Initial: NU=',N_PL,IEB_PL(1:N_PL)
c      WRITE(6,'(A12,I4,4X,16I3)')' Checks: NU=',N_YA,IEB_YA(1:N_YA)
c      WRITE(6,'(A12,I4,4X,16I3)')' virtau: NU=',NY_B,NYSTOR(1:NY_B)
c      ENDIF
900   DO 901 NY=1,N_YA
       DO L=1,NY_B
       IF(NY.EQ.NYSTOR(L)) GOTO 901
       ENDDO
       I=IEB_YA(NY)
       IEB_PL(I) = 0
901   CONTINUE

c      IF(DIAGNO)
c    &WRITE(6,'(A12,I4,4X,16I3)')'Correct: NU=',N_PL,IEB_PL(1:N_PL)

      I=1
902   IF(IEB_PL(I).EQ.0) THEN
         N_PL=N_PL-1
         IF(N_PL.LT.I) GOTO 903
         DO J=I,N_PL
          IEB_PL(J) = IEB_PL(J+1)
         ENDDO
      ELSE
         I = I+1
         IF(I.GT.N_PL) GOTO 903
      ENDIF
      GOTO 902

c903   IF(DIAGNO) THEN
c      WRITE(6,'(A12,I4,4X,16I3)')'  Final: NU=',N_PL,IEB_PL(1:N_PL)
c      WRITE(6,*)'zeitmn,T_MAX=',zeitmn,T_MAX
c      stop
c      ENDIF
903   IF(N_PL.EQ.0) THEN
        IPLATE= 0
        ITRIA = 0
        RETURN
      ELSE
        GOTO 250
      ENDIF
      END SUBROUTINE TLMAX_TRAVEL2

      END SUBROUTINE TLMAX_TRAVEL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      NAME: LEAUSR                                                    C
C  FUNCTION: DETERMINE THE CELL NUMBERS OF A GIVEN POINT(X0,Y0,Z0)     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION LEAUSR(X0,Y0,Z0)

      USE GEOMETRY_PL
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      IMPLICIT NONE            
      INTEGER ::  LEAUSR
     I,           J,IRGEN,NEW,ICNXT,NPANU,ICOS,IER,ISTS
     I,           ITR_STO,NTR_STO
     I,           IG,IG1,IG2,JP,JP1,JP2,JPD,L1,L2,ID,JD
      REAL*8  ::  X0,Y0,Z0,VELX,VELY,VELZ,TIMET
     R,           R_0,Z_0,PHI,PHI1,PHI2,RJ,PJ,TJ,DIST,DISTA
     R,           R,Z,XNORM,YNORM,ZNORM,SCOS
     R,           R_CENTER,Z_CENTER,D1,D2,D3,A,B,C1,C2
      LOGICAL ::  SURF_P

      LOGICAL:: DIAGNO=.FALSE.
      SAVE
C-----------------------------------------------------------------------
C  INPUT:
C     (X0,Y0,Z0 ):  A GIVEN POINT P0
C OUTPUT:
C          LEAUSR:  CELL NUMBER
C-----------------------------------------------------------------------

      LEAUSR = 0

      IF(DIAGNO) THEN
        WRITE(6,*)'====> LEAUSR DETERMINES CELL NUMBER'
        WRITE(6,'(A7,1p,3E12.4)')' INPUT:',X0,Y0,Z0
      ENDIF

      IF( N0S_DEF == 0) GOTO 200

      R_0  = SQRT(X0**2 + Y0**2)
      Z_0  = Z0

C  Phi-plane on which the particle is located
      IER = 1
      PHI = ATAN2(Y0,X0)
      LOOP_NZ0 : DO NZ0_N0=0,NZONET-1
      LOOP_JT0 : DO JT0_N0=0,ZON_TORO(NZ0_N0)-1
        PHI1 = PHI_PLANE(PHI_PL_OS(NZ0_N0)+JT0_N0  )
        PHI2 = PHI_PLANE(PHI_PL_OS(NZ0_N0)+JT0_N0+1)
        IF( PHI>= PHI1 .AND. PHI<= PHI2) THEN 
          IER = 0
          EXIT LOOP_NZ0
        ENDIF 
      ENDDO LOOP_JT0
      ENDDO LOOP_NZ0

      IF( IER/=0) THEN
        IF(DIAGNO)
     .  WRITE(6,*)'PHI,PHI_0,PHI_1=',PHI,PHI_PLANE(0)
     .,            PHI_PLANE(PHI_PL_OS(NZONET))
        RETURN
      ENDIF 

      JR0_N0 = (ZON_RADI(NZ0_N0)-1)/2
      
      IG1=  JT0_N0*SRF_POLO(NZ0_N0)*SRF_RADI(NZ0_N0)
     .+     GRID_P_OS(NZ0_N0)
      IG2= IG1 + SRF_RADI(NZ0_N0)*ZON_POLO(NZ0_N0)/2

      R_CENTER = 0.5*(RG(IG1)+RG(IG2))
      Z_CENTER = 0.

      A =  Z_0 - Z_CENTER
      B =-(R_0 - R_CENTER)

      JPD = (ZON_POLO(NZ0_N0)-1)/6
      JP1 = 0

      IG1 = JR0_N0 + (JP1+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+               GRID_P_OS(NZ0_N0)
      D1 = (RG(IG1)-R_CENTER)*A + (ZG(IG1)-Z_CENTER)*B
      C1 =-(RG(IG1)-R_CENTER)*B + (ZG(IG1)-Z_CENTER)*A

      LOOP_SREACH_JP1_2: DO L1=0,ZON_POLO(NZ0_N0)
        JP2 = JP1 + JPD
        IF(JP2 > ZON_POLO(NZ0_N0)) JP2 = ZON_POLO(NZ0_N0)
        
        IG2 = JR0_N0+(JP2+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+               GRID_P_OS(NZ0_N0)
        D2  = (RG(IG2)-R_CENTER)*A + (ZG(IG2)-Z_CENTER)*B
        C2  =-(RG(IG2)-R_CENTER)*B + (ZG(IG2)-Z_CENTER)*A

      
        IF( D1*D2.LE.0. .AND. (C1>0. .OR. C2>0.) ) THEN

         IG = IG1
         LOOP_SREACH_JP0_N0 : DO 
          IF(JP2-JP1<=1) THEN
           IF(-(RG(IG)-R_CENTER)*B+(ZG(IG)-Z_CENTER)*A > 0.) THEN
             JP0_N0 = JP1
             EXIT LOOP_SREACH_JP1_2 
           ELSE
             EXIT LOOP_SREACH_JP0_N0
           ENDIF
          ENDIF 

          JP= (JP1+JP2)/2
          IG= JR0_N0+(JP+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+                GRID_P_OS(NZ0_N0)
          D3= (RG(IG)-R_CENTER)*A + (ZG(IG)-Z_CENTER)*B
          IF    (D1*D3 <= 0.) THEN 
             JP2 = JP
             D2  = D3
          ELSEIF(D2*D3 <= 0.) THEN
             JP1 = JP
             D1  = D3
          ELSE  
             IER = 1
             EXIT LOOP_SREACH_JP1_2
          ENDIF
         ENDDO LOOP_SREACH_JP0_N0

        ELSEIF(  JP2==ZON_POLO(NZ0_N0) ) THEN  
             IER = 2
             EXIT LOOP_SREACH_JP1_2
        ENDIF 
        D1 = D2
        C1 = C2
        IG1=IG2
        JP1=JP2
      ENDDO LOOP_SREACH_JP1_2

      IF( IER/=0) RETURN 

C Now, we have JR0,JP0,JT0
C Find the grid point close to the (X0,Y0,Z0)
      IG= JR0_N0+(JP0_N0+JT0_N0*SRF_POLO(NZ0_N0))*SRF_RADI(NZ0_N0)
     .+           GRID_P_OS(NZ0_N0)

      DISTA = (RG(IG)-R_0)**2 + (ZG(IG)-Z_0)**2

      ID = 1
      JD = 1
      LOOP_SMALLEST_D: DO

       IF(JR0_N0 == 0 .OR. JR0_N0 == ZON_RADI(NZ0_N0)-1 
     ..OR.JP0_N0 == 0 .OR. JP0_N0 == ZON_POLO(NZ0_N0)-1) 
     . EXIT LOOP_SMALLEST_D

       IG1 = IG

       IG2 = IG + ID
       DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
       IF(DIST < DISTA) THEN
         DISTA = DIST
         JR0_N0= JR0_N0 + ID
         IG    = IG2
       ELSE
         IG2 = IG - ID
         DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
         IF(DIST < DISTA) THEN
           DISTA = DIST
           JR0_N0= JR0_N0 - ID
           IG    = IG2
           ID    =-ID
         ENDIF
       ENDIF

       IG2 = IG + SRF_RADI(NZ0_N0)*JD
       DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
       IF(DIST < DISTA) THEN
         DISTA = DIST
         JP0_N0= JP0_N0 + JD
         IG    = IG2
       ELSE
         IG2 = IG - SRF_RADI(NZ0_N0)*JD
         DIST= (RG(IG2)-R_0)**2 + (ZG(IG2)-Z_0)**2
         IF(DIST < DISTA) THEN
           DISTA = DIST
           JP0_N0= JP0_N0 - JD
           IG    = IG2
           JD    =-JD
         ENDIF
       ENDIF
       
       IF(IG == IG1) EXIT LOOP_SMALLEST_D
      ENDDO LOOP_SMALLEST_D
C*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
200   NTR_STO = NTRIANG
      ITR_STO = ITRIA
      NTRIANG = 0

      RJ = 0.
      PJ = 0.
      TJ = JT0_N0 + 0.00001               

      CALL RZ_REAL_COORDINATES
     .     (NZ0_N0,JR0_N0,JP0_N0,RJ,PJ,TJ,R,Z)
      IG1 = PHI_PL_OS(NZ0_N0) + JT0_N0
      IG2 = IG1 + 1
      PHI = PHI_PLANE(IG1)+(TJ-JT0_N0)*(PHI_PLANE(IG2)-PHI_PLANE(IG1))

      XP0 = R*COS(PHI)
      YP0 = R*SIN(PHI)
      ZP0 = Z
  
      DISTA= SQRT((X0-XP0)**2+(Y0-YP0)**2+(Z0-ZP0)**2)
      VELX = (X0-XP0)/DISTA
      VELY = (Y0-YP0)/DISTA
      VELZ = (Z0-ZP0)/DISTA

      SURF_P = .FALSE.
      IC0_N0 = IDCELL(JR0_N0+(JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))
     .*               ZON_RADI(NZ0_N0) + MESH_P_OS(NZ0_N0)   )
    
      NEW = 0

      IF(DIAGNO) THEN
        WRITE(6,*)'DISTANCE TO (X0,Y0,Z0):',DISTA,' CM'
        WRITE(6,*)'NZ0   JR0  JP0  JT0     DIST(CM)'
      ENDIF

      ICOS  = 1
      
      IF(DIAGNO) ICOS = 2
      TRACE_LOOP : DO 
            CALL TIMUSR(IRGEN,XP0,YP0,ZP0,VELX,VELY,VELZ,NEW,
     .                  ICNXT,TIMET,ICOS,IER,IRGEN,SURF_P)

                 IF(DIAGNO) WRITE(6,'(4I5,F12.6)')
     .           NZ0_N0,JR0_N0,JP0_N0,JT0_N0,TIMET

            IF(TIMET > DISTA) THEN
              LEAUSR = IC0_N0
              IF(DIAGNO) WRITE(6,*) '****LEAUSR=',IC0_N0
              ISRF  = 0
              EXIT TRACE_LOOP
            ENDIF

            IF(TIMET <= 0.) THEN
              WRITE(6,*)'IERR In CALLING TIMUSR'
              
              EXIT TRACE_LOOP
            ELSEIF(ICNXT == -1) THEN
              CALL NORUSR(ISTS,XP0,YP0,ZP0,XNORM,YNORM,ZNORM
     .,                   SCOS,VELX,VELY,VELZ,IRGEN)
              DISTA = DISTA - TIMET
              SURF_P=.TRUE.
              NEW   = 0
            ELSEIF(ICNXT <0 ) THEN
              WRITE(6,*)'A non-transparent surface'
              WRITE(6,*)'MSURF=',MSURF
              EXIT TRACE_LOOP
            ENDIF 
      END DO TRACE_LOOP
      NTRIANG = NTR_STO 
      ITRIA   = ITR_STO
      RETURN
      END FUNCTION LEAUSR
C ===== SOURCE: trace_test.f
      SUBROUTINE TEST_TRACE_N0()
      USE GEOMETRY_PL
      USE SURFACE_PL
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      LOGICAL :: SURF_P
      NZ0_N0 = 1
      JR0_N0 = 10
      JP0_N0 = 60
      JT0_N0 = 17
       
      RJ0 = 0.
      PJ0 = 0.
      TJ0 = JT0_N0 + 0.8

      IG = JR0_N0+(JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))*ZON_RADI(NZ0_N0)
     .+    MESH_P_OS(NZ0_N0)
      IC0_N0 = IDCELL(IG)
  
      CALL RZ_REAL_COORDINATES
     .          (NZ0_N0,JR0_N0,JP0_N0,RJ0,PJ0,TJ0,R0,Z0)
      K1 = JT0_N0 + PHI_PL_OS(NZ0_N0)
      K2 = K1  + 1
      X0 = R0*(COSPHI(K1)+(COSPHI(K2)-COSPHI(K1))*(TJ0-JT0_N0))
      Y0 = R0*(SINPHI(K1)+(SINPHI(K2)-SINPHI(K1))*(TJ0-JT0_N0))

      VELX =-1.
      VELY = 0.
      VELZ = 0.
      NEW  = 0
      SURF_P = .FALSE.
      NPANU  = 1
      
      ITR = 0
      OPEN(99,FILE='out')
      write(99,'(3f12.6)')x0,y0,z0
      DO 
            CALL TIMUSR(IRGEN,X0,Y0,Z0,VELX,VELY,VELZ,NEW,
     .                  ICNXT,TIMET,ICOS,IER,ITR,SURF_P)
            IF(TIMET <= 0.) THEN 
              WRITE(6,*)'TIMET=',TIMET
              EXIT
            ELSEIF(ICNXT == -1) THEN
              CALL NORUSR(ISTS,X0,Y0,Z0,XNORM,YNORM,ZNORM
     .,                   SCOS,VELX,VELY,VELZ,IRGEN)
              TIMET = 0.
              SURF_P=.TRUE.
              NEW   = 0
c             write(99,'(3f12.6)')xfin,yfin,zfin
c             ITR = ITR + 1
c             IF(ITR == 100) EXIT
            ELSEIF(ICNXT<0) THEN
              X0 = X0 + VELX*TIMET
              Y0 = Y0 + VELY*TIMET
              Z0 = Z0 + VELZ*TIMET
              CALL NORUSR(ISTS,X0,Y0,Z0,XNORM,YNORM,ZNORM
     .,                   SCOS,VELX,VELY,VELZ,IRGEN)
              TIMET = 0.
              SURF_P=.TRUE.
              NEW   = 0 
              A     = VELX*XNORM + VELY*YNORM + VELZ*ZNORM
              VELX  = VELX*(1.-XNORM**2)-(2.*A-VELX*XNORM)*XNORM
              VELY  = VELY*(1.-YNORM**2)-(2.*A-VELY*YNORM)*YNORM
              VELZ  = VELZ*(1.-ZNORM**2)-(2.*A-VELZ*ZNORM)*ZNORM

c             write(99,'(3f12.6)')xfin,yfin,zfin
              ITR = ITR + 1
              IF(ITR == 100000) EXIT
            ENDIF
      ENDDO 
      close(99)
      END SUBROUTINE TEST_TRACE_N0
      SUBROUTINE TEST_LEAUSR()
      USE GEOMETRY_PL
      USE SURFACE_PL
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      DO I=1,100000
      NZ0 = 2.*RANF()
      JR0 = RANF()*ZON_RADI(NZ0)
      JP0 = RANF()*ZON_POLO(NZ0)
      JT0 = RANF()*ZON_TORO(NZ0)

      RJ0 = 0.
      PJ0 = 0.
      TJ0 = JT0 + 0.00005

      IG = JR0+(JP0+JT0*ZON_POLO(NZ0))*ZON_RADI(NZ0)
     .+    MESH_P_OS(NZ0)
      IC0= IDCELL(IG)
 
      CALL RZ_REAL_COORDINATES
     .          (NZ0,JR0,JP0,RJ0,PJ0,TJ0,R0,Z0)
      K1 = JT0 + PHI_PL_OS(NZ0)
      K2 = K1  + 1
      X0 = R0*(COSPHI(K1)+(COSPHI(K2)-COSPHI(K1))*(TJ0-JT0))
      Y0 = R0*(SINPHI(K1)+(SINPHI(K2)-SINPHI(K1))*(TJ0-JT0))

      II = LEAUSR(X0,Y0,Z0,I)
      IF( NZ0_N0/=NZ0 .OR. JR0_N0/=JR0 .OR. JP0_N0/=JP0 .OR.
     .    JT0_N0/=JT0 .OR. IC0_N0/=IC0) THEN
      WRITE(6,*)'I=',I
      WRITE(6,'(5I6)')NZ0,JR0,JP0,JT0,IC0
      WRITE(6,'(5I6)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0,IC0_N0          
      ENDIF
     
      ENDDO 
      END
C ===== SOURCE: user.f
*DK USER
C
C   USER SUPPLIED SUBROUTINES
C
      SUBROUTINE ADDUSR
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
      USE PARMMOD
      USE CADGEO
      USE CCONA
      USE CGEOM
      USE CGRID
      USE CLGIN
      USE CINIT

      IMPLICIT NONE
C
C MODIFY GEOMETRY
C
C
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2, LGJUM3
C   SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   LGJUM3(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN IN ZELLE NCELL=J SITZT
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C             LGJUM3(J,I)=.FALSE. FUER ALLE I UND J
C
C  SETZE EINIGE VOLUMINA EXPLIZIT
C
C
      RETURN
      END
C
      SUBROUTINE GEOUSR
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CCONA
      USE CGEOM
      USE CGRID
      USE CLGIN
      USE CINIT

      IMPLICIT NONE
C
C MODIFY GEOMETRY
C
C
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2, LGJUM3
C   SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   LGJUM3(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN IN ZELLE NCELL=J SITZT
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C             LGJUM3(J,I)=.FALSE. FUER ALLE I UND J
C
C  SETZE EINIGE VOLUMINA EXPLIZIT
C
C
      RETURN
      END
C
C
      SUBROUTINE PLAUSR
      IMPLICIT NONE
C
      RETURN
      END
C
      SUBROUTINE PLTUSR(PLABLE,J)
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CPL3D
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: PLABLE
      INTEGER, INTENT(IN) :: J
      RETURN
      END
C
      SUBROUTINE UPTUSR(XSTOR2,XSTORV2,WV)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CSPEZ
      USE CUPD
      USE CGRID
      USE CZT1
      USE CLOGAU
      USE COMXS
      USE COMPRT
      USE CCONA
      USE CGEOM
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      REAL(DP) :: DIST, WTR, CVELBB
      INTEGER :: ICOU, IRD, IRDD
C
      IF (ITYP.EQ.1) THEN
        DO 51 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          CVELBB=SQRT(RMASSA(IATM))/CVELAA
          IRDD=NRCELL+NUPC(ICOU)*NR1ST+NBLCKA
          IRD=NCLTAL(IRDD)
          ADDV(IATM,IRD)=ADDV(IATM,IRD)+WTR*VELX*VEL*CVELBB
          ADDV(NATMI+IATM,IRD)=ADDV(NATMI+IATM,IRD)+WTR*VELY*VEL*
     .                         CVELBB
          ADDV(2*NATMI+IATM,IRD)=ADDV(2*NATMI+IATM,IRD)+WTR*VELZ*
     .                           VEL*CVELBB
51      CONTINUE
C
      ELSEIF (ITYP.EQ.2) THEN
c
      ELSEIF (ITYP.EQ.3) THEN
c
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE COMXS
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
C
C  COLLISION ESTIMATOR  FOR DENSITY
C
C     IF (ITYP.NE.3) RETURN
C
C     COLV(IION,NCELL)=COLV(IION,NCELL)+WEIGHT/SIGTOT
C
      RETURN
      END
C
      SUBROUTINE UPNUSR
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE COMPRT
      USE CGRID
      IMPLICIT NONE
      INTEGER :: IRD, IRDD
C
C  SNAPSHOT ESTIMATOR FOR DENSITY
C
      IF (ITYP.NE.1) RETURN
      IF (IATM.GT.NSNVI) RETURN
C
      IRDD=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2
      IRD=NCLTAL(IRDD)
      SNAPV(IATM,IRD)=SNAPV(IATM,IRD)+WEIGHT
C
      RETURN
      END
C
      SUBROUTINE UPSUSR(WT,IND)
C
C  SURFACE AVERAGED TALLIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WT
      INTEGER, INTENT(IN) :: IND
      RETURN
      END
C
C
      SUBROUTINE REFUSR
C  USER SUPPLIED SURFACE REFLECTION MODEL
C
C  INCIDENT PARTICLE WITH MASS NUMBER XMP AND NUCLEAR CHARGE NUMBER XCP
C
C  RETURN 1 : EIRENE DEFAULT SURFACE ANGULAR DISTRIBUTION
C  RETURN 2 : EIRENE DEFAULT THERMAL MOLECULE MODEL
C  RETURN 3 : EIRENE DEFAULT THERMAL ATOM MODEL
C  RETURN 4 : ABSORB PARTICLE
C
C
C  THIS ROUTINE: SEKI ET AL, NUCL.FUS. VOL 20 , NO 10, (1980)
C  FOR D,T AND HE ATOMS OR IONS INCIDENT ON COPPER AND ON IRON
C
C  MODREF=3 WALL (REFLECT TEST IONS, SPECULAR FOR FAST ATOMS)
C  MODREF=4 DEFLECTOR PLATE  (ABSORB TEST IONS, COSINE FOR FAST ATOMS)
C  EL: 1ST INDEX I1: SPECIES OF INCIDENT PARTICLE
C      2ND INDEX: SURFACE MATERIAL
C  I1=1: D,T ;  I1=2: HE
C  I2=1: SS ;  I2=2: COPPER
      USE PRECISION
      USE PARMMOD
      USE COMPRT
      USE CZT1
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XMW, XCW, XMPP, XCP, ZCOS, ZSIN, EXPI,
     .                        RPROB, E0TERM
      INTEGER, INTENT(IN) :: IGASF,IGAST
      REAL(DP), SAVE :: EL(2,2)=(/2610.,5510.,2990.,6290./)
      REAL(DP) :: RN, RE, ELI, XMP, ZCS, ZSN, RPR, ZEP1, E, EX,
     .            E0T 
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: I2, ITYPO, I1, IGT, IGF

      RN(E,ELI)=-0.237*LOG10(E/ELI)+0.19
      RE(E,ELI)=-0.22*LOG10(E/ELI)+0.06
C
      ENTRY RF0USR
      RETURN

      I2=MODREF-2
C  ABSORB TEST IONS AT DEFLECTOR PLATE
      IF (ITYP.EQ.3.AND.I2.EQ.2) RETURN 4
      ITYPO=ITYP
C
      I1=1
      IF (XMP.Ge.4.) I1=2
      write (6,*) 'refusr',ityp,i2,i1
      ELI=EL(I1,I2)
C
      ZCS=1.
      ZSN=0.
      RPR=RN(E0,ELI)
      ZEP1=RANF_EIRENE( )
C
C COSINE FOR CHARGED, SPECULAR FOR NEUTRALS
      EX=0.
      IF (ITYPO.LE.3) EX=100.
C
      IF (I1.EQ.1) THEN
C  INCIDENT D, T , T+ OR D+
        IGT=1
        IF (XMP.GT.2.) IGT=2
        IGT=IGF
        ITYP=1
C  FRANCK CONDON ATOM, 3 EV?
        IF (ZEP1.GT.RPR) THEN
          EX=0.
          E0T=3.
          RETURN 3
C  FAST ATOM
        ELSE
          IATM=1
          IF (XMP.GT.2.) IATM=2
          E0=RE(E0,ELI)/RPR*E0
          VEL=RSQDVA(IATM)*SQRT(E0)
          RETURN 1
        ENDIF
      ELSE
C  INCIDENT HE, HE+ OR HE++
        IGF=3
        IGT=IGF
        ITYP=1
C  THERMAL ATOM, 0.1 EV
        EX=0.
        E0T=0.1
        IF (I2.EQ.1) E0T=0.0612
        IF (ZEP1.GT.RPR) RETURN 3
C  FAST ATOM
        IATM=3
        E0=RE(E0,ELI)/RPR*E0
        VEL=RSQDVA(IATM)*SQRT(E0)
        RETURN 1
      ENDIF
C
      RETURN
C
      ENTRY RF1USR (XMW,XCW,XMPP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
C
      RETURN
      ENTRY SP0USR
      ENTRY SP1USR
      RETURN
      END
C
C
      SUBROUTINE RETUSR(SIG)
      USE PRECISION
      USE PARMMOD
      USE COMPRT
      USE COMUSR
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: SIG
      IF (E0.LE.0.5*TEIN(NCELL)) SIG=-1.

      RETURN
      END
C
C
      SUBROUTINE MODUSR(*)
      IMPLICIT NONE
C RETURN 1 Ohne Berechnung von neuen plasma
C Redefination of new sources
C FLUX(ISTRA) = Flux in AMP
C ISTRA     === 1
C NITER InTENT(IN)
C IITER 
C     CALL (FLUX(ISTRA),NLVOL,NLSRF,NLPNT,NLPLS,NLATM,NLMOL)
C     RETURN 1
      RETURN
      END
c
c
      SUBROUTINE SIGUSR(III,JJJ,ZDS,PEN,PSIG,TIMAX,ARGST)
      IMPLICIT NONE
      WRITE(6,*)'NOTHING HAS BEEN DONE IN SIGUSR!'
      RETURN
      END
c
c
      subroutine talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE PRECISION
      USE PARMMOD
      USE CGRID
      USE CGEOM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl,txtsp,txtun
      integer :: N0WR

C3. OUTPUT FOR EMC3
      N0WR = 1
      CALL OUT_FOR_EMC3(N0WR)

      ILAST = 1
      return 1
      end
c
C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGEOM
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      RETURN
      END
C
c
      SUBROUTINE TMSUSR(T)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: T
      RETURN
      END


      SUBROUTINE INIUSR
      USE PRECISION
      USE PARMUSR
      IMPLICIT NONE
      CHARACTER(72) :: ERRMES
 
      SAVE
C
C1. Updata 3D Geometry
c     IERR = 0
c     CALL UP_GEO_N0(N1ST,IERR,ERRMES)
c     IF(IERR.NE.0) THEN 
c       CALL WRMESS('YOU CANNOT RUN THE STELLARATOR OPTION DUE TO:')
c       CALL WRMESS(ERRMES)
c       call exit
c     ENDIF 
c     CALL WRMESS('*** Updata successful ***')

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* UPTCOP: stellarator option
C
C THIS SUBROUTINE CALCULATES 
C     1. CHARGE EXCHANGE RATES
C     2. MOMENTUM LOSS OR GAIN (CHARGE EXCHANGE AND IONIZATION BY e)        
C DUE ONLY TO NEUTRAL ATOM
C
C *1) CONTRIBUTION OF ION IMPACT IONIZATION IS NOT INCLUDED
C *2) CONTRIBUTION OF MOLECULE IS IN SUBP. COLLIDE              
C *3) NO CONTRIBUTION FROM TEST ION
C
C*** STELLARATOR OPTION
C    1. VECTORS(B-LINE, V_ion) ARE PROVIDED BY THE USER
C       1.1  V_ION = VDION (FUNCTION)
C       1.2  BX,BY,BZ PROVIDED BY SUB. VECUSR
C    2. VELOCITY OF THE NEUTRAL (VELX,-Y,-Z) IS REPLACED BY (VX,VY,VZ)
C       PROVIDED BY USER. THIS IS DUE TO THE PERIOD OF THE STELLA. 
C       GEOMETRY.                                                    
 
      SUBROUTINE UPTCOP(XSTOR2,WV)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CUPD
      USE COMXS
      USE CSPEZ
      USE CGRID
      USE CLOGAU
      USE CCONA
      USE CPOLYG
      USE CZT1
      USE CSDVI
      USE CINIT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                      XSTORV2(NSTORV,N2ND+N3RD), WV
      REAL(DP) :: P, WTRSIG, EION, V0_PARB, DIST, WTR
      INTEGER :: IAEL, IREL, IPL2, IAEI, IRDS, IBGK, IICX, IIEI, IIEL,
     .           IMEL, IPL1, I, IPL, IIO, IRD, IP, IR, IML, IAT, IFIRST,
     .           IRCX, IADD, ICOU, IACX, IRDD, IMCX, IMEI, IPL2, IPLSTI
      INTEGER, SAVE :: NMTSP
      REAL(DP), ALLOCATABLE, SAVE ::
     . CNDYNA(:), CNDYNM(:), CNDYNI(:), CNDYNP(:) 
      DATA IFIRST/0/
      SAVE
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        ALLOCATE (CNDYNA(NATM))
        ALLOCATE (CNDYNM(NMOL))
        ALLOCATE (CNDYNI(NION))
        ALLOCATE (CNDYNP(NPLS))
        DO 11 IAT=1,NATMI
11        CNDYNA(IAT)=AMUA*RMASSA(IAT)
        DO 12 IML=1,NMOLI
12        CNDYNM(IML)=AMUA*RMASSM(IML)
        DO 13 IIO=1,NIONI
13        CNDYNI(IIO)=AMUA*RMASSI(IIO)
        DO 14 IPL=1,NPLSI
14        CNDYNP(IPL)=AMUA*RMASSP(IPL)
C
        DO 2 IR=1,NR1STM
C         DO 2 IP=1,NP2NDM
CYHF
          IP_MAX = MAX(1,NP2NDM)
          DO 2 IP=1,IP_MAX
C END YHF 
            IRD=IR+(IP-1)*NR1P2

            IF ((INDPRO(4) == 8) .AND. (INDPRO(5) == 8)) THEN
              vion=vdion(ird)
              VSIG_PARB(1:NPLSI,ird)=CNDYNP(1:NPLSI)*vion
            ELSE
              IF (INDPRO(5) == 8) THEN
                CALL VECUSR (1,BX,BY,BZ,1)
              ELSE
                BX=BXIN(IRD)
                BY=BYIN(IRD)
                BZ=BZIN(IRD)
              END IF
              DO 3 IPL=1,NPLSI
                IF (INDPRO(4) == 8) THEN
                  CALL VECUSR (2,VX,VY,VZ,IPL)
                ELSE
                  VX=VXIN(IPL,IRD)
                  VY=VYIN(IPL,IRD)
                  VZ=VZIN(IPL,IRD)
                END IF
                VSIG_PARB(IPL,IRD)=CNDYNP(IPL)*(VX*BX+VY*BY+VZ*BZ)
3             CONTINUE
            END IF
2       CONTINUE

      NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI
C
      ENDIF
C
C  WV=WEIGHT/VEL
C
C  ATOMS
      IF (ITYP.EQ.1) THEN
        DO 20 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 20
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C  1,NPLSI:
C              PARTICLE CHARGE EXCHANGE RATE DUE TO IPLS: #/S
C              WITH ATOM SPECIES IATM=1,NATMI, PER ION
C  EACH RATE IS WEIGHTED WITH THE FACTOR (E0/EI-1), E0 BEING
C  THE NEUTRAL PARTCILE ENERGY, EI THE MEAN PLASMA ION ENERGY
C  THESE RATES ARE SCALED IN THE SHORT CYCLE WITH EI*NI
C
C
          IF (NCPVI.LT.NPLSI) GOTO 20
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 51
          DO 52 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 52
            IPLSTI= MPLSTI(IPLS)
            EION=1.5*TIIN(IPLSTI,IRD)+EDRIFT(IPLS,IRD)
            WTRSIG=WTR*SIGVCX(IRCX)/DIIN(IPLS,IRD)
CYHF 
C           COPV(IPLS,IRDD)=COPV(IPLS,IRDD)+WTRSIG*(E0/EION-1.)
CEND YHF 
            COPV(IPLS,IRDD)=COPV(IPLS,IRDD)+WTRSIG
            LMETSP(NMTSP+IPLS)=.TRUE.
52        CONTINUE
51        CONTINUE
C
C.........................................
C
C   MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C.........................................
C
C
C  CONTRIBUTIONS FROM ATOMS
C  NPLSI+1, 2*NPLSI:
C
          IF (NCPVI.LT.2*NPLSI) GOTO 20
C
          IADD=NPLSI
          IF (INDPRO(5) == 8) THEN
            CALL VECUSR (1,BX,BY,BZ,1)
          ELSE
            BX=BXIN(IRD)
            BY=BYIN(IRD)
            BZ=BZIN(IRD)
          END IF
          V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
          V0_PARB=V0_PARB*CNDYNA(IATM)
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 59
          DO 56 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 56
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 56
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHA
C
            VSIGCX_PARB=VSIG_PARB(IPLS,IRD)
C
            WTRSIG=WTR*SIGVCX(IRCX)
            WTRSIG=WTRSIG*SIGN(1.D0,VSIGCX_PARB)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*VSIGCX_PARB
            LMETSP(NMTSP+IPL1)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL2=IADD+N1STX(IRCX,2)
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*VSIGCX_PARB
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=IADD+N2NDX(IRCX,2)
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
56        CONTINUE
59        CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 61 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 62 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  VSIGEI_PARB=VSIG_PARB(IPL,IRD)
                  WTRSIG=WTR*SIGVEI(IRDS)*P
                  WTRSIG=WTRSIG*SIGN(1.D0,VSIGEI_PARB)
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
62            CONTINUE
            ENDIF
61        CONTINUE
C
C  ION IMPACT IONIZATION CONTRIBUTION: NOT INCLUDED
C
C
C  ELASTIC CONTRIBUTION FROM ATOMS
C
C
          IF (LGAEL(IATM,0,0).EQ.0) GOTO 80
C  DEFAULT TRACKLENGTH ESTIMATOR (BGK APPROXIMATION)
          DO 81  IAEL=1,NAELI(IATM)
            IREL=LGAEL(IATM,IAEL,0)
            IPLS=LGAEL(IATM,IAEL,1)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 81
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 81
C
            VSIGEL_PARB=VSIG_PARB(IPLS,IRD)
            WTRSIG=WTR*SIGVEL(IREL)
            WTRSIG=WTRSIG*SIGN(1.D0,VSIGEL_PARB)
C
            IPL1=IADD+IPLS
            IPL2=IPL1
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*VSIGEL_PARB
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
            LMETSP(NMTSP+IPL2)=.TRUE.
81        CONTINUE
80      CONTINUE
C
20      CONTINUE
C
C  MOLECULES
      ELSEIF (ITYP.EQ.2) THEN
C
        DO 200 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 200
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C             MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C
C
C
C  CONTRIBUTIONS FROM MOLECULES
C  2*NPLSI+1, 3*NPLSI:
C
          IF (NCPVI.LT.3*NPLSI) GOTO 200
C
          IADD=2*NPLSI
CYHF
C         V0_PARB=VEL*(VELX*BXIN(IRD)+VELY*BYIN(IRD)+VELZ*BZIN(IRD)) 
          IF (INDPRO(5) == 8) THEN
            CALL VECUSR (1,BX,BY,BZ,1)
          ELSE
            BX=BXIN(IRD)
            BY=BYIN(IRD)
            BZ=BZIN(IRD)
          END IF
          V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
CEND YHF 
          V0_PARB=V0_PARB*CNDYNM(IMOL)
C
          IF (LGMCX(IMOL,0,0).EQ.0) GOTO 590
          DO 560 IMCX=1,NMCXI(IMOL)
            IRCX=LGMCX(IMOL,IMCX,0)
            IPLS=LGMCX(IMOL,IMCX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 560
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 560
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHM
C
            VSIGCX_PARB=VSIG_PARB(IPLS,IRD)
C
            WTRSIG=WTR*SIGVCX(IRCX)
            WTRSIG=WTRSIG*SIGN(1.D0,VSIGCX_PARB)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*VSIGCX_PARB
            LMETSP(NMTSP+IPL1)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL2=IADD+N1STX(IRCX,2)
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*VSIGCX_PARB
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=IADD+N2NDX(IRCX,2)
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
560       CONTINUE
590       CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 610 IMEI=1,NMDSI(IMOL)
            IRDS=LGMEI(IMOL,IMEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 620 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  VSIGEI_PARB=VSIG_PARB(IPL,IRD)
                  WTRSIG=WTR*SIGVEI(IRDS)*P
                  WTRSIG=WTRSIG*SIGN(1.D0,VSIGEI_PARB)
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
620           CONTINUE
            ENDIF
610       CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM MOLECULES
C
C
          IF (LGMEL(IMOL,0,0).EQ.0) GOTO 800
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 810 IMEL=1,NMELI(IMOL)
            IREL=LGMEL(IMOL,IMEL,0)
            IPLS=LGMEL(IMOL,IMEL,1)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 810
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 810
C
            VSIGEL_PARB=VSIG_PARB(IPLS,IRD)
            WTRSIG=WTR*SIGVEL(IREL)
            WTRSIG=WTRSIG*SIGN(1.D0,VSIGEL_PARB)
C
            IPL1=IADD+IPLS
            IPL2=IPL1
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*VSIGEL_PARB
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
            LMETSP(NMTSP+IPL2)=.TRUE.
810       CONTINUE
800     CONTINUE
C
C
200     CONTINUE
C
C  TEST IONS
C
      ELSEIF (ITYP.EQ.3) THEN
C
        DO 2000 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRDD=NCLTAL(IRD)
C
          IF (LGVAC(IRD,0)) GOTO 2000
C
          XSTOR(:,:) = XSTOR2(:,:,ICOU)
          XSTORV(:) = XSTORV2(:,ICOU)
C
C             MOMENTUM EXCHANGE RATE: DYN/CM**3
C
C
C
C
C  CONTRIBUTIONS FROM TEST IONS
C  3*NPLSI+1, 4*NPLSI:
C
          IF (NCPVI.LT.4*NPLSI) GOTO 2000
C
          IADD=3*NPLSI
          IF (INDPRO(5) == 8) THEN
            CALL VECUSR (1,BX,BY,BZ,1)
          ELSE
            BX=BXIN(IRD)
            BY=BYIN(IRD)
            BZ=BZIN(IRD)
          END IF
          V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
          V0_PARB=V0_PARB*CNDYNI(IION)
C
          IF (LGICX(IION,0,0).EQ.0) GOTO 5900
          DO 5600 IICX=1,NICXI(IION)
            IRCX=LGICX(IION,IICX,0)
            IPLS=LGICX(IION,IICX,1)
            IF (LGVAC(IRD,IPLS)) GOTO 5600
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 5600
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) NOT AVAILABLE
C             FROM FUNCTION FPATHI
C
            VSIGCX_PARB=VSIG_PARB(IPLS,IRD)
C
            WTRSIG=WTR*SIGVCX(IRCX)
            WTRSIG=WTRSIG*SIGN(1.D0,VSIGCX_PARB)
C  PREVIOUS BULK ION IPLS, NOW LOST
            IPL1=IADD+IPLS
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*VSIGCX_PARB
            LMETSP(NMTSP+IPL1)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL2=IADD+N1STX(IRCX,2)
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*VSIGCX_PARB
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=IADD+N2NDX(IRCX,2)
              COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
              LMETSP(NMTSP+IPL2)=.TRUE.
            ENDIF
5600      CONTINUE
5900      CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 6100 IIEI=1,NIDSI(IION)
            IRDS=LGIEI(IION,IIEI)
            IF (PPLDS(IRDS,0).GT.0) THEN
              DO 6200 IPL=1,NPLSI
                P=PPLDS(IRDS,IPL)
                IF (P.GT.0) THEN
                  VSIGEI_PARB=VSIG_PARB(IPL,IRD)
                  WTRSIG=WTR*SIGVEI(IRDS)*P
                  WTRSIG=WTRSIG*SIGN(1.D0,VSIGEI_PARB)
C  NEW BULK ION IPL
                  IPL2=IADD+IPL
                  COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
                  LMETSP(NMTSP+IPL2)=.TRUE.
                ENDIF
6200          CONTINUE
            ENDIF
6100      CONTINUE
C
C
C  ELASTIC CONTRIBUTION FROM TEST IONS
C
          IF (LGIEL(IION,0,0).EQ.0) GOTO 8000
C  DEFAULT TRACKLENGTH ESTIMATOR
          DO 8100 IIEL=1,NIELI(IION)
            IREL=LGIEL(IION,IIEL,0)
            IPLS=LGIEL(IION,IIEL,1)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 8100
C  THIS TALLY IS A BGK TALLY. IT SHOULD NOT BE UPDATED HERE.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 8100
C
            VSIGEL_PARB=VSIG_PARB(IPLS,IRD)
            WTRSIG=WTR*SIGVEL(IREL)
            WTRSIG=WTRSIG*SIGN(1.D0,VSIGEL_PARB)
C
            IPL1=IADD+IPLS
            IPL2=IPL1
            COPV(IPL1,IRDD)=COPV(IPL1,IRDD)-WTRSIG*VSIGEL_PARB
            LMETSP(NMTSP+IPL1)=.TRUE.
            IPL2=IPL1
            COPV(IPL2,IRDD)=COPV(IPL2,IRDD)+WTRSIG*V0_PARB
            LMETSP(NMTSP+IPL2)=.TRUE.
8100      CONTINUE
8000    CONTINUE
C
C
2000    CONTINUE
C
C
      ENDIF
C
      RETURN
      END


      SUBROUTINE OUT_PROF(MFLAG,IPLAS,IC1,IC2,PROF) 
C MFLAG = 1: Ionisation source for IPLS 
C         2: Energy source for e 
C         3: Energy source for i 
C        30: CX Rate
C        31: Momentum source from atom
C        32: Momentum source from mole
C        33: Momentum source from TEST ION 
C
C        -1: Atom density    
C        -2: Molecule density 
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMSOU
      USE CADGEO
      USE CSPEZ
      USE CSPEI
      USE CSDVI
      USE CGEOM
      USE CTRCEI
      USE CTEXT
      USE COUTAU
      USE CGRID
      USE COMPRT
      USE CLOGAU
      USE CCONA

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: PROF(*)
      INTEGER, INTENT(IN) :: MFLAG, IPLAS, IC1, IC2
      INTEGER :: IAD
 
      IF    (MFLAG.EQ.1) THEN
        PROF(IC1:IC2)       
     .=                     PAPL(IPLAS,IC1:IC2)
     .+                     PMPL(IPLAS,IC1:IC2)
     .+                     PIPL(IPLAS,IC1:IC2)
      ELSEIF(MFLAG.EQ.2) THEN
        PROF(IC1:IC2)                  
     .=                     EAEL(IC1:IC2)
     .+                     EMEL(IC1:IC2)
     .+                     EIEL(IC1:IC2)
      ELSEIF(MFLAG.EQ.3) THEN
        PROF(IC1:IC2)                  
     .=                     EAPL(IC1:IC2)
     .+                     EMPL(IC1:IC2)
     .+                     EIPL(IC1:IC2)
      ELSEIF(MFLAG.EQ.30) THEN
C CX Rate
       IAD= 0*NPLSI + IPLAS
       PROF(IC1:IC2)
     .=                     COPV(IAD,IC1:IC2)
      ELSEIF(MFLAG.EQ.31) THEN
C CONTRIBUTION FROM ATOM
       IAD= 1*NPLSI + IPLAS
       PROF(IC1:IC2)
     .=                     COPV(IAD,IC1:IC2)
      ELSEIF(MFLAG.EQ.32) THEN
C CONTRIBUTION FROM MOLE 
       IAD= 2*NPLSI + IPLAS
       PROF(IC1:IC2)
     .=                     COPV(IAD,IC1:IC2)
      ELSEIF(MFLAG.EQ.33) THEN
C CONTRIBUTION FROM TEST ION 
       IAD= 3*NPLSI + IPLAS
       PROF(IC1:IC2)
     .=                     COPV(IAD,IC1:IC2)
 
      ELSEIF(MFLAG.EQ.-1) THEN
       PROF(IC1:IC2)
     .=                     PDENA(IPLAS,IC1:IC2)
      ELSEIF(MFLAG.EQ.-2) THEN
       PROF(IC1:IC2)
     .=                     PDENM(IPLAS,IC1:IC2)
      ENDIF
 
      RETURN 
      END
C*DK SAMUSR                                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SUBROUTINE SAMUSR(NI,X,Y,Z,AD1,AD2,AD3,AD4,AD5,AD6,               
     .                  IR,IP,IT,IA,IB,                            
     .                  TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   USER SUPPLIED SURFACE SAMPLING ROUTINE FOR EIRENE INPUT DS "W7AS" C 
C                                                                     C 
C   INPUT:        NI:  SOURCE NUMBER                                  C
C            AD1-AD6:  IRRELEVANT                                     C 
C                                                                     C 
C  OUTPUT:     X,Y,Z:  SOURCE POINT                                   C 
C           IR,IP,IT:  CELL NUMBER OF THE POINT SOURCE                C
C                                                                     C 
C  IA=0,  IB=1                                                        C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C---> COMMON BLOCKS                                                     
      USE PRECISION                     
      USE PARMMOD                      
      USE COMSOU                                                   
      USE CUPD                                
      USE CCONA                         
      USE CGRID                             
      USE CTRCEI                        
      USE CLOGAU                            
      USE CADGEO                               
      USE COMUSR
      IMPLICIT NONE                        
      REAL(DP), INTENT(IN) :: AD1, AD2, AD3, AD4, AD5, AD6, 
     .                        R1, R2, R3, R4, R5, R6
      REAL(DP), INTENT(OUT) :: X,Y,Z,TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,
     .                       WEISPZ
      INTEGER, INTENT(IN) :: NI,i1, i2
      INTEGER, INTENT(OUT) :: IR, IP, IT, IA, IB
      INTEGER :: NELEM
      LOGICAL L_ADD_EIR 
      SAVE
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      ENTRY SM0USR(I1,I2,R1,R2,R3,R4,R5,R6)
C Uebergabe des Flusses von EMC3
C     I2 = ISTRA
      IF( NLVOL(I2) ) THEN 
c       CALL UEBERSCHREIBUNG_VOL (I2, FLUX(I2) )
      ELSEIF( NLSRF(I2) ) THEN
c       CALL UEBERSCHREIBUNG_SRF (I2, FLUX(I2) )
      ENDIF 
      RETURN

      ENTRY SM1USR (NI,X,Y,Z,AD1,AD2,AD3,AD4,AD5,AD6,
     .              IR,IP,IT,IA,IB,
     .              TEWL,TIWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
C--------------------------------------------------------------------
      CALL SAMPLE_N0_W7(L_ADD_EIR,NELEM,X,Y,Z)                          
C Uebergabe von WEIGHT : Relative meaning

      IF(L_ADD_EIR) THEN 
       NASOR(NI,ISTRA) = 0                                               
       INSOR(NI,ISTRA) = NELEM                                           
       NRSOR(NI,ISTRA) = 0      
      ENDIF 
      RETURN
99    END
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE INICOP
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMSOU
      USE CSTEP
      USE COMPRT
      USE CINIT
      USE CGRID
      USE CCONA
      USE CGEOM
      USE CADGEO
      USE CSPEI
      USE CTRCEI
      USE CPLOT
      USE CLOGAU
      USE CTEXT
      USE CLGIN
      USE COUTAU
      USE CPOLYG
      USE CSDVI
      USE CGRPTL
      USE COMXS
      USE CCOUPL
      USE COMNNL
      USE CZT1
      IMPLICIT NONE
      INTEGER :: IPLS
     
      NCPVI=4*NPLSI
      IF (NCPVI.GT.NCPV) THEN
        WRITE (6,*) 'FROM INTERFACING SUBROUTINE INFCOP: '
        CALL MASPRM('NCPV',4,NCPV,'NCPVI',5,NCPVI,IERROR)
        STOP            
      ENDIF
      DO 70 IPLS=1,NPLSI
        ICPVE(IPLS)=1
        ICPRC(IPLS)=1
        TXTTAL(IPLS,NTALB)=
     .  'ENERGY WEIGHTED CX RATE OF ATOMS WITH IPLS                  '
        TXTSPC(IPLS,NTALB)=TEXTS(IPLS)
        TXTUNT(IPLS,NTALB)='AMPS                       '
C
        ICPVE(NPLSI+IPLS)=3
        ICPRC(NPLSI+IPLS)=1
        TXTTAL(NPLSI+IPLS,NTALB)=
     .  'PAR. MOM. SOURCE, FROM ATOMS, FOR IPLS             '
        TXTSPC(NPLSI+IPLS,NTALB)=TEXTS(IPLS)
        TXTUNT(NPLSI+IPLS,NTALB)='G*CM/S* AMP * CM**-3       '
C
        ICPVE(2*NPLSI+IPLS)=3
        ICPRC(2*NPLSI+IPLS)=2
        TXTTAL(2*NPLSI+IPLS,NTALB)=
     .  'PAR. MOM. SOURCE, FROM MOLECUELS, FOR IPLS         '
        TXTSPC(2*NPLSI+IPLS,NTALB)=TEXTS(IPLS)
        TXTUNT(2*NPLSI+IPLS,NTALB)='G*CM/S* AMP * CM**-3       '
C
        ICPVE(3*NPLSI+IPLS)=3
        ICPRC(3*NPLSI+IPLS)=1
        TXTTAL(3*NPLSI+IPLS,NTALB)=
     .  'PAR. MOM. SOURCE, FROM TEST IONS, FOR IPLS               '
        TXTSPC(3*NPLSI+IPLS,NTALB)=TEXTS(IPLS)
        TXTUNT(3*NPLSI+IPLS,NTALB)='G*CM/S* AMP * CM**-3       '
C
70    CONTINUE
      RETURN
      END

      SUBROUTINE BROAD_USR
      RETURN
      END
C ===== SOURCE: user_add.f
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
