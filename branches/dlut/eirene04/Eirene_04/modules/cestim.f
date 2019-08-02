      MODULE CESTIM

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CESTIM, DEALLOC_CESTIM, ASSOCIATE_CESTIM,
     P          INIT_CESTIM, EIRENE_SPECTRUM, SPECT_ARRAY, 
     P          ASSIGNMENT(=)

      PRIVATE :: SPEC_TO_SPEC
      TYPE EIRENE_SPECTRUM
        REAL(DP) :: SPCMIN, SPCMAX, SPCDEL, SPCDELI, SPCINT
        REAL(DP) :: SGMS, STVS, EES
        INTEGER :: NSPC, ISPCTYP, ISPCSRF, IPRTYP, IPRSP, IMETSP
        REAL(DP), DIMENSION(:), POINTER :: SPC, SDV, SGM
      END TYPE EIRENE_SPECTRUM

      TYPE SPECT_ARRAY
        TYPE(EIRENE_SPECTRUM), POINTER :: PSPC
      END TYPE SPECT_ARRAY

      INTERFACE ASSIGNMENT(=)  ! DEFINE ASSIGNMENT
        MODULE PROCEDURE SPEC_TO_SPEC
      END INTERFACE

      TYPE(SPECT_ARRAY), PUBLIC, ALLOCATABLE, SAVE :: ESTIML(:) 
      TYPE(SPECT_ARRAY), PUBLIC, ALLOCATABLE, SAVE :: SMESTL(:)

      INTEGER, PUBLIC, SAVE ::
     I NESTM1, NESTM2, NESTIM

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        ESTIMV(:,:), ESTIMS(:,:)

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: 
     R          CEMETERYV(:,:), CEMETERYS(:,:)

C  NESTM1, REAL, VOLUME AVERAGED TALLIES
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R PDENA(:,:), PDENM(:,:), PDENI(:,:), PDENPH(:,:),
     R EDENA(:,:), EDENM(:,:), EDENI(:,:), EDENPH(:,:),
     R PAEL(:),    PAAT(:,:),  PAML(:,:),  PAIO(:,:),  PAPHT(:,:), 
     R PAPL(:,:),
     R PMEL(:),    PMAT(:,:),  PMML(:,:),  PMIO(:,:),  PMPHT(:,:),
     R PMPL(:,:),
     R PIEL(:),    PIAT(:,:),  PIML(:,:),  PIIO(:,:),  PIPHT(:,:),
     R PIPL(:,:),
     R PPHEL(:),   PPHAT(:,:), PPHML(:,:), PPHIO(:,:), PPHPHT(:,:),
     R PPHPL(:,:),
     R EAEL(:),  EAAT(:),  EAML(:),  EAIO(:),  EAPHT(:),  EAPL(:), 
     R EMEL(:),  EMAT(:),  EMML(:),  EMIO(:),  EMPHT(:),  EMPL(:), 
     R EIEL(:),  EIAT(:),  EIML(:),  EIIO(:),  EIPHT(:),  EIPL(:),
     R EPHEL(:), EPHAT(:), EPHML(:), EPHIO(:), EPHPHT(:), EPHPL(:), 
     R ADDV(:,:),  COLV(:,:),  SNAPV(:,:),
     R COPV(:,:),  BGKV(:,:),  ALGV(:,:),
     R PGENA(:,:), PGENM(:,:), PGENI(:,:), PGENPH(:,:),
     R EGENA(:,:), EGENM(:,:), EGENI(:,:), EGENPH(:,:),
     R VGENA(:,:), VGENM(:,:), VGENI(:,:), VGENPH(:,:),
     R PPAT(:,:),  PPML(:,:),  PPIO(:,:),  PPPHT(:,:), PPPL(:,:),
     R EPAT(:),    EPML(:),    EPIO(:),    EPPHT(:),   EPPL(:)

C  NESTM2, REAL, SURFACE AVERAGED TALLIES
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R POTAT(:,:), 
     R PRFAAT(:,:), PRFMAT(:,:), PRFIAT(:,:), PRFPHAT(:,:),
     R PRFPAT(:,:),
C
     R POTML(:,:),
     R PRFAML(:,:), PRFMML(:,:), PRFIML(:,:), PRFPHML(:,:),
     R PRFPML(:,:),
C
     R POTIO(:,:),
     R PRFAIO(:,:), PRFMIO(:,:), PRFIIO(:,:), PRFPHIO(:,:),
     R PRFPIO(:,:),
C
     R POTPHT(:,:),
     R PRFAPHT(:,:), PRFMPHT(:,:), PRFIPHT(:,:), PRFPHPHT(:,:),
     R PRFPPHT(:,:), 
C
     R POTPL(:,:)
C
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R EOTAT(:,:),
     R ERFAAT(:,:), ERFMAT(:,:), ERFIAT(:,:), ERFPHAT(:,:),
     R ERFPAT(:,:),
C
     R EOTML(:,:),
     R ERFAML(:,:), ERFMML(:,:), ERFIML(:,:), ERFPHML(:,:),
     R ERFPML(:,:),
C
     R EOTIO(:,:),
     R ERFAIO(:,:), ERFMIO(:,:), ERFIIO(:,:), ERFPHIO(:,:),
     R ERFPIO(:,:),
C
     R EOTPHT(:,:),
     R ERFAPHT(:,:), ERFMPHT(:,:), ERFIPHT(:,:), ERFPHPHT(:,:),
     R ERFPPHT(:,:), 
C
     R EOTPL(:,:),
C
     R SPTAT(:,:), SPTML(:,:),
     R SPTIO(:,:), SPTPHT(:,:), SPTPL(:,:),
     R SPTTOT(:),
     R ADDS(:,:),  ALGS(:,:),
     R SPUMP(:,:)

C  FROM HERE: NO EQUIVALENCE
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NFIRST(:), NADDV(:),
     I IRESC1(:), IRESC2(:),
     I NFRSTW(:), NADDW(:)

      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE :: 
     L LIVTALV(:), LIVTALS(:)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LPDENA, LPDENM, LPDENI, LPDENPH,
     L LEDENA, LEDENM, LEDENI, LEDENPH,
     L LPAEL,  LPAAT,  LPAML,  LPAIO,   LPAPHT,  LPAPL,
     L LPMEL,  LPMAT,  LPMML,  LPMIO,   LPMPHT,  LPMPL,
     L LPIEL,  LPIAT,  LPIML,  LPIIO,   LPIPHT,  LPIPL,
     L LPPHEL, LPPHAT, LPPHML, LPPHIO,  LPPHPHT, LPPHPL,
     L LEAEL,  LEAAT,  LEAML,  LEAIO,   LEAPHT,  LEAPL, 
     L LEMEL,  LEMAT,  LEMML,  LEMIO,   LEMPHT,  LEMPL, 
     L LEIEL,  LEIAT,  LEIML,  LEIIO,   LEIPHT,  LEIPL,
     L LEPHEL, LEPHAT, LEPHML, LEPHIO,  LEPHPHT, LEPHPL, 
     L LADDV,  LCOLV,  LSNAPV,
     L LCOPV,  LBGKV,  LALGV,
     L LPGENA, LPGENM, LPGENI, LPGENPH,
     L LEGENA, LEGENM, LEGENI, LEGENPH,
     L LVGENA, LVGENM, LVGENI, LVGENPH,
     L LPPAT,  LPPML,  LPPIO,  LPPPHT,  LPPPL,
     L LEPAT,  LEPML,  LEPIO,  LEPPHT,  LEPPL

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LPOTAT, 
     L LPRFAAT, LPRFMAT, LPRFIAT, LPRFPHAT,
     L LPRFPAT,
C
     L LPOTML,
     L LPRFAML, LPRFMML, LPRFIML, LPRFPHML,
     L LPRFPML,
C
     L LPOTIO,
     L LPRFAIO, LPRFMIO, LPRFIIO, LPRFPHIO,
     L LPRFPIO,
C
     L LPOTPHT,
     L LPRFAPHT, LPRFMPHT, LPRFIPHT, LPRFPHPHT,
     L LPRFPPHT, 
C
     L LPOTPL
C
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LEOTAT,
     L LERFAAT, LERFMAT, LERFIAT, LERFPHAT, LERFPAT,
C
     L LEOTML,
     L LERFAML, LERFMML, LERFIML, LERFPHML, LERFPML,
C
     L LEOTIO,
     L LERFAIO, LERFMIO, LERFIIO, LERFPHIO, LERFPIO,
C
     L LEOTPHT,
     L LERFAPHT, LERFMPHT, LERFIPHT, LERFPHPHT, LERFPPHT, 
C
     L LEOTPL,
C
     L LSPTAT, LSPTML, LSPTIO, LSPTPHT, LSPTPL,
     L LSPTTOT,
     L LADDS,  LALGS,
     L LSPUMP
C

      CONTAINS


      SUBROUTINE ALLOC_CESTIM(ICAL)
      
      INTEGER, INTENT(IN) :: ICAL
      
      IF (ICAL == 1) THEN
        
        IF (ALLOCATED(LIVTALV)) RETURN

        ALLOCATE (LIVTALV(NTALV))
        ALLOCATE (LIVTALS(NTALS))

        ALLOCATE (NFIRST(NTALV)) 
        ALLOCATE (NADDV(NTALV))
        ALLOCATE (IRESC1(NTALV)) 
        ALLOCATE (IRESC2(NTALV))
        ALLOCATE (NFRSTW(NTALS)) 
        ALLOCATE (NADDW(NTALS))
        
        WRITE (55,'(A,T25,I15)')
     .      ' CESTIM ',4*(NTALV+NTALS) + (4*NTALV+2*NTALS)*4 

      ELSE IF (ICAL == 2) THEN
      
        IF (ALLOCATED(ESTIMV)) RETURN

        NESTM1=NVOLTL*NRTAL
        NESTM2=NSRFTL*NLMPGS
        NESTIM=NESTM1+NESTM2
        
        ALLOCATE (ESTIMV(NVOLTL,NRTAL))
        ALLOCATE (ESTIMS(NSRFTL,NLMPGS))
        
        ALLOCATE (CEMETERYV(0:0,NRTAL))
        ALLOCATE (CEMETERYS(0:0,NLMPGS))
        
        WRITE (55,'(A,T25,I15)')
     .      ' CESTIM ',(NESTIM+NRTAL+NLMPGS)*8 

      END IF

      CALL INIT_CESTIM(ICAL)

      RETURN
      END SUBROUTINE ALLOC_CESTIM


      SUBROUTINE ASSOCIATE_CESTIM

C  VOLUME AVERAGED TALLIES


      IF (LPDENA) THEN
        PDENA => ESTIMV(NADDV(1)+1:NADDV(2),:)
      ELSE
        PDENA => CEMETERYV(0:0,:)
      END IF
      IF (LPDENM) THEN
        PDENM => ESTIMV(NADDV(2)+1:NADDV(3),:)
      ELSE
        PDENM => CEMETERYV(0:0,:)
      END IF
      IF (LPDENI) THEN
        PDENI => ESTIMV(NADDV(3)+1:NADDV(4),:)
      ELSE
        PDENI => CEMETERYV(0:0,:)
      END IF
      IF (LPDENPH) THEN
        PDENPH => ESTIMV(NADDV(4)+1:NADDV(5),:)
      ELSE
        PDENPH => CEMETERYV(0:0,:)
      END IF

      IF (LEDENA) THEN
        EDENA => ESTIMV(NADDV(5)+1:NADDV(6),:)
      ELSE
        EDENA => CEMETERYV(0:0,:)
      END IF
      IF (LEDENM) THEN
        EDENM => ESTIMV(NADDV(6)+1:NADDV(7),:)
      ELSE
        EDENM => CEMETERYV(0:0,:)
      END IF
      IF (LEDENI) THEN
        EDENI => ESTIMV(NADDV(7)+1:NADDV(8),:)
      ELSE
        EDENI => CEMETERYV(0:0,:)
      END IF
      IF (LEDENPH) THEN
        EDENPH => ESTIMV(NADDV(8)+1:NADDV(9),:)
      ELSE
        EDENPH => CEMETERYV(0:0,:)
      END IF

      IF (LPAEL) THEN
        PAEL => ESTIMV(NADDV(10),:)
      ELSE
        PAEL => CEMETERYV(0,:)
      END IF
      IF (LPAAT) THEN
        PAAT => ESTIMV(NADDV(10)+1:NADDV(11),:)
      ELSE
        PAAT => CEMETERYV(0:0,:)
      END IF
      IF (LPAML) THEN
        PAML => ESTIMV(NADDV(11)+1:NADDV(12),:)
      ELSE
        PAML => CEMETERYV(0:0,:)
      END IF
      IF (LPAIO) THEN
        PAIO => ESTIMV(NADDV(12)+1:NADDV(13),:)
      ELSE
        PAIO => CEMETERYV(0:0,:)
      END IF
      IF (LPAPHT) THEN
        PAPHT => ESTIMV(NADDV(13)+1:NADDV(14),:)
      ELSE
        PAPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPAPL) THEN
        PAPL => ESTIMV(NADDV(14)+1:NADDV(15),:)
      ELSE
        PAPL => CEMETERYV(0:0,:)
      END IF

      IF (LPMEL) THEN
        PMEL => ESTIMV(NADDV(16),:)
      ELSE
        PMEL => CEMETERYV(0,:)
      END IF
      IF (LPMAT) THEN
        PMAT => ESTIMV(NADDV(16)+1:NADDV(17),:)
      ELSE
        PMAT => CEMETERYV(0:0,:)
      END IF
      IF (LPMML) THEN
        PMML => ESTIMV(NADDV(17)+1:NADDV(18),:)
      ELSE
        PMML => CEMETERYV(0:0,:)
      END IF
      IF (LPMIO) THEN
        PMIO => ESTIMV(NADDV(18)+1:NADDV(19),:)
      ELSE
        PMIO => CEMETERYV(0:0,:)
      END IF
      IF (LPMPHT) THEN
        PMPHT => ESTIMV(NADDV(19)+1:NADDV(20),:)
      ELSE
        PMPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPMPL) THEN
        PMPL => ESTIMV(NADDV(20)+1:NADDV(21),:)
      ELSE
        PMPL => CEMETERYV(0:0,:)
      END IF

      IF (LPIEL) THEN
        PIEL => ESTIMV(NADDV(22),:)
      ELSE
        PIEL => CEMETERYV(0,:)
      END IF
      IF (LPIAT) THEN
        PIAT => ESTIMV(NADDV(22)+1:NADDV(23),:)
      ELSE
        PIAT => CEMETERYV(0:0,:)
      END IF
      IF (LPIML) THEN
        PIML => ESTIMV(NADDV(23)+1:NADDV(24),:)
      ELSE
        PIML => CEMETERYV(0:0,:)
      END IF
      IF (LPIIO) THEN
        PIIO => ESTIMV(NADDV(24)+1:NADDV(25),:)
      ELSE
        PIIO => CEMETERYV(0:0,:)
      END IF
      IF (LPIPHT) THEN
        PIPHT => ESTIMV(NADDV(25)+1:NADDV(26),:)
      ELSE
        PIPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPIPL) THEN
        PIPL => ESTIMV(NADDV(26)+1:NADDV(27),:)
      ELSE
        PIPL => CEMETERYV(0:0,:)
      END IF

      IF (LPPHEL) THEN
        PPHEL => ESTIMV(NADDV(28),:)
      ELSE
        PPHEL => CEMETERYV(0,:)
      END IF
      IF (LPPHAT) THEN
        PPHAT => ESTIMV(NADDV(28)+1:NADDV(29),:)
      ELSE
        PPHAT => CEMETERYV(0:0,:)
      END IF
      IF (LPPHML) THEN
        PPHML => ESTIMV(NADDV(29)+1:NADDV(30),:)
      ELSE
        PPHML => CEMETERYV(0:0,:)
      END IF
      IF (LPPHIO) THEN
        PPHIO => ESTIMV(NADDV(30)+1:NADDV(31),:)
      ELSE
        PPHIO => CEMETERYV(0:0,:)
      END IF
      IF (LPPHPHT) THEN
        PPHPHT => ESTIMV(NADDV(31)+1:NADDV(32),:)
      ELSE
        PPHPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPPHPL) THEN
        PPHPL => ESTIMV(NADDV(32)+1:NADDV(33),:)
      ELSE
        PPHPL => CEMETERYV(0:0,:)
      END IF

      IF (LEAEL) THEN
        EAEL => ESTIMV(NADDV(34),:)
      ELSE
        EAEL => CEMETERYV(0,:)
      END IF
      IF (LEAAT) THEN
        EAAT => ESTIMV(NADDV(35),:)
      ELSE
        EAAT => CEMETERYV(0,:)
      END IF
      IF (LEAML) THEN
        EAML => ESTIMV(NADDV(36),:)
      ELSE
        EAML => CEMETERYV(0,:)
      END IF
      IF (LEAIO) THEN
        EAIO => ESTIMV(NADDV(37),:)
      ELSE
        EAIO => CEMETERYV(0,:)
      END IF
      IF (LEAPHT) THEN
        EAPHT => ESTIMV(NADDV(38),:)
      ELSE
        EAPHT => CEMETERYV(0,:)
      END IF
      IF (LEAPL) THEN
        EAPL => ESTIMV(NADDV(39),:)
      ELSE
        EAPL => CEMETERYV(0,:)
      END IF

      IF (LEMEL) THEN
        EMEL => ESTIMV(NADDV(40),:)
      ELSE
        EMEL => CEMETERYV(0,:)
      END IF
      IF (LEMAT) THEN
        EMAT => ESTIMV(NADDV(41),:)
      ELSE
        EMAT => CEMETERYV(0,:)
      END IF
      IF (LEMML) THEN
        EMML => ESTIMV(NADDV(42),:)
      ELSE
        EMML => CEMETERYV(0,:)
      END IF
      IF (LEMIO) THEN
        EMIO => ESTIMV(NADDV(43),:)
      ELSE
        EMIO => CEMETERYV(0,:)
      END IF
      IF (LEMPHT) THEN
        EMPHT => ESTIMV(NADDV(44),:)
      ELSE
        EMPHT => CEMETERYV(0,:)
      END IF
      IF (LEMPL) THEN
        EMPL => ESTIMV(NADDV(45),:)
      ELSE
        EMPL => CEMETERYV(0,:)
      END IF

      IF (LEIEL) THEN
        EIEL => ESTIMV(NADDV(46),:)
      ELSE
        EIEL => CEMETERYV(0,:)
      END IF
      IF (LEIAT) THEN
        EIAT => ESTIMV(NADDV(47),:)
      ELSE
        EIAT => CEMETERYV(0,:)
      END IF
      IF (LEIML) THEN
        EIML => ESTIMV(NADDV(48),:)
      ELSE
        EIML => CEMETERYV(0,:)
      END IF
      IF (LEIIO) THEN
        EIIO => ESTIMV(NADDV(49),:)
      ELSE
        EIIO => CEMETERYV(0,:)
      END IF
      IF (LEIPHT) THEN
        EIPHT => ESTIMV(NADDV(50),:)
      ELSE
        EIPHT => CEMETERYV(0,:)
      END IF
      IF (LEIPL) THEN
        EIPL => ESTIMV(NADDV(51),:)
      ELSE
        EIPL => CEMETERYV(0,:)
      END IF

      IF (LEPHEL) THEN
        EPHEL => ESTIMV(NADDV(52),:)
      ELSE
        EPHEL => CEMETERYV(0,:)
      END IF
      IF (LEPHAT) THEN
        EPHAT => ESTIMV(NADDV(53),:)
      ELSE
        EPHAT => CEMETERYV(0,:)
      END IF
      IF (LEPHML) THEN
        EPHML => ESTIMV(NADDV(54),:)
      ELSE
        EPHML => CEMETERYV(0,:)
      END IF
      IF (LEPHIO) THEN
        EPHIO => ESTIMV(NADDV(55),:)
      ELSE
        EPHIO => CEMETERYV(0,:)
      END IF
      IF (LEPHPHT) THEN
        EPHPHT => ESTIMV(NADDV(56),:)
      ELSE
        EPHPHT => CEMETERYV(0,:)
      END IF
      IF (LEPHPL) THEN
        EPHPL => ESTIMV(NADDV(57),:)
      ELSE
        EPHPL => CEMETERYV(0,:)
      END IF

      IF (LADDV) THEN
        ADDV => ESTIMV(NADDV(NTALA)+1:NADDV(NTALA+1),:)
      ELSE
        ADDV => CEMETERYV(0:0,:)
      END IF
      IF (LCOLV) THEN
        COLV => ESTIMV(NADDV(NTALC)+1:NADDV(NTALC+1),:)
      ELSE
        COLV => CEMETERYV(0:0,:)
      END IF
      IF (LSNAPV) THEN
        SNAPV => ESTIMV(NADDV(NTALT)+1:NADDV(NTALT+1),:)
      ELSE
        SNAPV => CEMETERYV(0:0,:)
      END IF
      IF (LCOPV) THEN
        COPV => ESTIMV(NADDV(NTALM)+1:NADDV(NTALM+1),:)
      ELSE
        COPV => CEMETERYV(0:0,:)
      END IF
      IF (LBGKV) THEN
        BGKV => ESTIMV(NADDV(NTALB)+1:NADDV(NTALB+1),:)
      ELSE
        BGKV => CEMETERYV(0:0,:)
      END IF
      IF (LALGV) THEN
        ALGV => ESTIMV(NADDV(NTALR)+1:NADDV(NTALR+1),:)
      ELSE
        ALGV => CEMETERYV(0:0,:)
      END IF

      IF (LPGENA) THEN
        PGENA => ESTIMV(NADDV(NTALV-21)+1:NADDV(NTALV-20),:)
      ELSE
        PGENA => CEMETERYV(0:0,:)
      END IF
      IF (LPGENM) THEN
        PGENM => ESTIMV(NADDV(NTALV-20)+1:NADDV(NTALV-19),:)
      ELSE
        PGENM => CEMETERYV(0:0,:)
      END IF
      IF (LPGENI) THEN
        PGENI => ESTIMV(NADDV(NTALV-19)+1:NADDV(NTALV-18),:)
      ELSE
        PGENI => CEMETERYV(0:0,:)
      END IF
      IF (LPGENPH) THEN
        PGENPH => ESTIMV(NADDV(NTALV-18)+1:NADDV(NTALV-17),:)
      ELSE
        PGENPH => CEMETERYV(0:0,:)
      END IF
      IF (LEGENA) THEN
        EGENA => ESTIMV(NADDV(NTALV-17)+1:NADDV(NTALV-16),:)
      ELSE
        EGENA => CEMETERYV(0:0,:)
      END IF
      IF (LEGENM) THEN
        EGENM => ESTIMV(NADDV(NTALV-16)+1:NADDV(NTALV-15),:)
      ELSE
        EGENM => CEMETERYV(0:0,:)
      END IF
      IF (LEGENI) THEN
        EGENI => ESTIMV(NADDV(NTALV-15)+1:NADDV(NTALV-14),:)
      ELSE
        EGENI => CEMETERYV(0:0,:)
      END IF
      IF (LEGENPH) THEN
        EGENPH => ESTIMV(NADDV(NTALV-14)+1:NADDV(NTALV-13),:)
      ELSE
        EGENPH => CEMETERYV(0:0,:)
      END IF
      IF (LVGENA) THEN
        VGENA => ESTIMV(NADDV(NTALV-13)+1:NADDV(NTALV-12),:)
      ELSE
        VGENA => CEMETERYV(0:0,:)
      END IF
      IF (LVGENM) THEN
        VGENM => ESTIMV(NADDV(NTALV-12)+1:NADDV(NTALV-11),:)
      ELSE
        VGENM => CEMETERYV(0:0,:)
      END IF
      IF (LVGENI) THEN
        VGENI => ESTIMV(NADDV(NTALV-11)+1:NADDV(NTALV-10),:)
      ELSE
        VGENI => CEMETERYV(0:0,:)
      END IF
      IF (LVGENPH) THEN
        VGENPH => ESTIMV(NADDV(NTALV-10)+1:NADDV(NTALV- 9),:)
      ELSE
        VGENPH => CEMETERYV(0:0,:)
      END IF

      IF (LPPAT) THEN
        PPAT => ESTIMV(NADDV(NTALV- 9)+1:NADDV(NTALV- 8),:)
      ELSE
        PPAT => CEMETERYV(0:0,:)
      END IF
      IF (LPPML) THEN
        PPML => ESTIMV(NADDV(NTALV- 8)+1:NADDV(NTALV- 7),:)
      ELSE
        PPML => CEMETERYV(0:0,:)
      END IF
      IF (LPPIO) THEN
        PPIO => ESTIMV(NADDV(NTALV- 7)+1:NADDV(NTALV- 6),:)
      ELSE
        PPIO => CEMETERYV(0:0,:)
      END IF
      IF (LPPPHT) THEN
        PPPHT => ESTIMV(NADDV(NTALV- 6)+1:NADDV(NTALV- 5),:)
      ELSE
        PPPHT => CEMETERYV(0:0,:)
      END IF
      IF (LPPPL) THEN
        PPPL => ESTIMV(NADDV(NTALV- 5)+1:NADDV(NTALV- 4),:)
      ELSE
        PPPL => CEMETERYV(0:0,:)
      END IF

      IF (LEPAT) THEN
        EPAT => ESTIMV(NADDV(NTALV- 3),:)
      ELSE
        EPAT => CEMETERYV(0,:)
      END IF
      IF (LEPML) THEN
        EPML => ESTIMV(NADDV(NTALV- 2),:)
      ELSE
        EPML => CEMETERYV(0,:)
      END IF
      IF (LEPIO) THEN
        EPIO => ESTIMV(NADDV(NTALV- 1),:)
      ELSE
        EPIO => CEMETERYV(0,:)
      END IF
      IF (LEPPHT) THEN
        EPPHT => ESTIMV(NADDV(NTALV  ),:)
      ELSE
        EPPHT => CEMETERYV(0,:)
      END IF
      IF (LEPPL) THEN
        EPPL => ESTIMV(NADDV(NTALV)+1 ,:)
      ELSE
        EPPL => CEMETERYV(0,:)
      END IF


C  SURFACE AVERAGED TALLIES


      IF (LPOTAT) THEN
        POTAT => ESTIMS(1:NADDW(2),:)
      ELSE
        POTAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAAT) THEN
        PRFAAT => ESTIMS(NADDW(2)+1:NADDW(3),:)
      ELSE
        PRFAAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMAT) THEN
        PRFMAT => ESTIMS(NADDW(3)+1:NADDW(4),:)
      ELSE
        PRFMAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIAT) THEN
        PRFIAT => ESTIMS(NADDW(4)+1:NADDW(5),:)
      ELSE
        PRFIAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHAT) THEN
        PRFPHAT => ESTIMS(NADDW(5)+1:NADDW(6),:)
      ELSE
        PRFPHAT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPAT) THEN
        PRFPAT => ESTIMS(NADDW(6)+1:NADDW(7),:)
      ELSE
        PRFPAT => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTML) THEN
        POTML => ESTIMS(NADDW(7)+1:NADDW(8),:)
      ELSE
        POTML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAML) THEN
        PRFAML => ESTIMS(NADDW(8)+1:NADDW(9),:)
      ELSE
        PRFAML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMML) THEN
        PRFMML => ESTIMS(NADDW(9)+1:NADDW(10),:)
      ELSE
        PRFMML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIML) THEN
        PRFIML => ESTIMS(NADDW(10)+1:NADDW(11),:)
      ELSE
        PRFIML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHML) THEN
        PRFPHML => ESTIMS(NADDW(11)+1:NADDW(12),:)
      ELSE
        PRFPHML => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPML) THEN
        PRFPML => ESTIMS(NADDW(12)+1:NADDW(13),:)
      ELSE
        PRFPML => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTIO) THEN
        POTIO => ESTIMS(NADDW(13)+1:NADDW(14),:)
      ELSE
        POTIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAIO) THEN
        PRFAIO => ESTIMS(NADDW(14)+1:NADDW(15),:)
      ELSE
        PRFAIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMIO) THEN
        PRFMIO => ESTIMS(NADDW(15)+1:NADDW(16),:)
      ELSE
        PRFMIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIIO) THEN
        PRFIIO => ESTIMS(NADDW(16)+1:NADDW(17),:)
      ELSE
        PRFIIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHIO) THEN
        PRFPHIO => ESTIMS(NADDW(17)+1:NADDW(18),:)
      ELSE
        PRFPHIO => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPIO) THEN
        PRFPIO => ESTIMS(NADDW(18)+1:NADDW(19),:)
      ELSE
        PRFPIO => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTPHT) THEN
        POTPHT => ESTIMS(NADDW(19)+1:NADDW(20),:)
      ELSE
        POTPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFAPHT) THEN
        PRFAPHT => ESTIMS(NADDW(20)+1:NADDW(21),:)
      ELSE
        PRFAPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFMPHT) THEN
        PRFMPHT => ESTIMS(NADDW(21)+1:NADDW(22),:)
      ELSE
        PRFMPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFIPHT) THEN
        PRFIPHT => ESTIMS(NADDW(22)+1:NADDW(23),:)
      ELSE
        PRFIPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPHPHT) THEN
        PRFPHPHT =>ESTIMS(NADDW(23)+1:NADDW(24),:)
      ELSE
        PRFPHPHT => CEMETERYS(0:0,:)
      END IF
      IF (LPRFPPHT) THEN
        PRFPPHT => ESTIMS(NADDW(24)+1:NADDW(25),:)
      ELSE
        PRFPPHT => CEMETERYS(0:0,:)
      END IF
C
      IF (LPOTPL) THEN
        POTPL => ESTIMS(NADDW(25)+1:NADDW(26),:)
      ELSE
        POTPL => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTAT) THEN
        EOTAT => ESTIMS(NADDW(26)+1:NADDW(27),:)
      ELSE
        EOTAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFAAT) THEN
        ERFAAT => ESTIMS(NADDW(27)+1:NADDW(28),:)
      ELSE
        ERFAAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFMAT) THEN
        ERFMAT => ESTIMS(NADDW(28)+1:NADDW(29),:)
      ELSE
        ERFMAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFIAT) THEN
        ERFIAT => ESTIMS(NADDW(29)+1:NADDW(30),:)
      ELSE
        ERFIAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHAT) THEN
        ERFPHAT => ESTIMS(NADDW(30)+1:NADDW(31),:)
      ELSE
        ERFPHAT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPAT) THEN
        ERFPAT => ESTIMS(NADDW(31)+1:NADDW(32),:)
      ELSE
        ERFPAT => CEMETERYS(0:0,:)
      END IF
C                  
      IF (LEOTML) THEN
        EOTML => ESTIMS(NADDW(32)+1:NADDW(33),:)
      ELSE
        EOTML => CEMETERYS(0:0,:)
      END IF
      IF (LERFAML) THEN
        ERFAML => ESTIMS(NADDW(33)+1:NADDW(34),:)
      ELSE
        ERFAML => CEMETERYS(0:0,:)
      END IF
      IF (LERFMML) THEN
        ERFMML => ESTIMS(NADDW(34)+1:NADDW(35),:)
      ELSE
        ERFMML => CEMETERYS(0:0,:)
      END IF
      IF (LERFIML) THEN
        ERFIML => ESTIMS(NADDW(35)+1:NADDW(36),:)
      ELSE
        ERFIML => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHML) THEN
        ERFPHML => ESTIMS(NADDW(36)+1:NADDW(37),:)
      ELSE
        ERFPHML => CEMETERYS(0:0,:)
      END IF
      IF (LERFPML) THEN
        ERFPML => ESTIMS(NADDW(37)+1:NADDW(38),:)
      ELSE
        ERFPML => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTIO) THEN
        EOTIO => ESTIMS(NADDW(38)+1:NADDW(39),:)
      ELSE
        EOTIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFAIO) THEN
        ERFAIO => ESTIMS(NADDW(39)+1:NADDW(40),:)
      ELSE
        ERFAIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFMIO) THEN
        ERFMIO => ESTIMS(NADDW(40)+1:NADDW(41),:)
      ELSE
        ERFMIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFIIO) THEN
        ERFIIO => ESTIMS(NADDW(41)+1:NADDW(42),:)
      ELSE
        ERFIIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHIO) THEN
        ERFPHIO => ESTIMS(NADDW(42)+1:NADDW(43),:)
      ELSE
        ERFPHIO => CEMETERYS(0:0,:)
      END IF
      IF (LERFPIO) THEN
        ERFPIO => ESTIMS(NADDW(43)+1:NADDW(44),:)
      ELSE
        ERFPIO => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTPHT) THEN
        EOTPHT => ESTIMS(NADDW(44)+1:NADDW(45),:)
      ELSE
        EOTPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFAPHT) THEN
        ERFAPHT => ESTIMS(NADDW(45)+1:NADDW(46),:)
      ELSE
        ERFAPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFMPHT) THEN
        ERFMPHT => ESTIMS(NADDW(46)+1:NADDW(47),:)
      ELSE
        ERFMPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFIPHT) THEN
        ERFIPHT => ESTIMS(NADDW(47)+1:NADDW(48),:)
      ELSE
        ERFIPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPHPHT) THEN
        ERFPHPHT =>ESTIMS(NADDW(48)+1:NADDW(49),:)
      ELSE
        ERFPHPHT => CEMETERYS(0:0,:)
      END IF
      IF (LERFPPHT) THEN
        ERFPPHT => ESTIMS(NADDW(49)+1:NADDW(50),:)
      ELSE
        ERFPPHT => CEMETERYS(0:0,:)
      END IF
C
      IF (LEOTPL) THEN
        EOTPL => ESTIMS(NADDW(50)+1:NADDW(51),:)
      ELSE
        EOTPL => CEMETERYS(0:0,:)
      END IF
C
      IF (LSPTAT) THEN
        SPTAT => ESTIMS(NADDW(51)+1:NADDW(52),:)
      ELSE
        SPTAT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTML) THEN
        SPTML => ESTIMS(NADDW(52)+1:NADDW(53),:)
      ELSE
        SPTML => CEMETERYS(0:0,:)
      END IF
      IF (LSPTIO) THEN
        SPTIO => ESTIMS(NADDW(53)+1:NADDW(54),:)
      ELSE
        SPTIO => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPHT) THEN
        SPTPHT => ESTIMS(NADDW(54)+1:NADDW(55),:)
      ELSE
        SPTPHT => CEMETERYS(0:0,:)
      END IF
      IF (LSPTPL) THEN
        SPTPL => ESTIMS(NADDW(55)+1:NADDW(56),:)
      ELSE
        SPTPL => CEMETERYS(0:0,:)
      END IF
      IF (LSPTTOT) THEN
        SPTTOT => ESTIMS(NADDW(56),:)
      ELSE
        SPTTOT => CEMETERYS(0,:)
      END IF
      IF (LADDS) THEN
        ADDS => ESTIMS(NADDW(57)+1:NADDW(58),:)
      ELSE
        ADDS => CEMETERYS(0:0,:)
      END IF
      IF (LALGS) THEN
        ALGS => ESTIMS(NADDW(58)+1:NADDW(59),:)
      ELSE
        ALGS => CEMETERYS(0:0,:)
      END IF
      IF (LSPUMP) THEN
        SPUMP => ESTIMS(NADDW(59)+1:,:)
      ELSE
        SPUMP => CEMETERYS(0:0,:)
      END IF

      RETURN
      END SUBROUTINE ASSOCIATE_CESTIM


      SUBROUTINE DEALLOC_CESTIM
C
      IF (ALLOCATED(ESTIMV)) THEN
         DEALLOCATE (ESTIMV)
         DEALLOCATE (ESTIMS)
         IF (NADSPC > 0) DEALLOCATE (ESTIML)
         DEALLOCATE (NFIRST)
         DEALLOCATE (NADDV)
         DEALLOCATE (IRESC1)
         DEALLOCATE (IRESC2)
         DEALLOCATE (NFRSTW)
         DEALLOCATE (NADDW)

         DEALLOCATE (CEMETERYV)
         DEALLOCATE (CEMETERYS)
      END IF

      IF (ALLOCATED(LIVTALV)) THEN
         DEALLOCATE (LIVTALV)
         DEALLOCATE (LIVTALS)
      END IF

      RETURN
      END SUBROUTINE DEALLOC_CESTIM


      SUBROUTINE INIT_CESTIM(ICAL)

      INTEGER, INTENT(IN) :: ICAL
C
      IF (ICAL == 1) THEN

        LIVTALV = .TRUE.
        LIVTALS = .TRUE.

! volume averaged tallies

        LPDENA  => LIVTALV(1) 
        LPDENM  => LIVTALV(2) 
        LPDENI  => LIVTALV(3) 
        LPDENPH => LIVTALV(4)
        LEDENA  => LIVTALV(5) 
        LEDENM  => LIVTALV(6) 
        LEDENI  => LIVTALV(7) 
        LEDENPH => LIVTALV(8)
        LPAEL   => LIVTALV(9)  
        LPAAT   => LIVTALV(10)  
        LPAML   => LIVTALV(11)  
        LPAIO   => LIVTALV(12)   
        LPAPHT  => LIVTALV(13)  
        LPAPL   => LIVTALV(14)
        LPMEL   => LIVTALV(15)  
        LPMAT   => LIVTALV(16)  
        LPMML   => LIVTALV(17)  
        LPMIO   => LIVTALV(18)   
        LPMPHT  => LIVTALV(19)  
        LPMPL   => LIVTALV(20)
        LPIEL   => LIVTALV(21)  
        LPIAT   => LIVTALV(22)  
        LPIML   => LIVTALV(23)  
        LPIIO   => LIVTALV(24)   
        LPIPHT  => LIVTALV(25)  
        LPIPL   => LIVTALV(26)
        LPPHEL  => LIVTALV(27) 
        LPPHAT  => LIVTALV(28) 
        LPPHML  => LIVTALV(29) 
        LPPHIO  => LIVTALV(30)  
        LPPHPHT => LIVTALV(31) 
        LPPHPL  => LIVTALV(32)
        LEAEL   => LIVTALV(33)  
        LEAAT   => LIVTALV(34)  
        LEAML   => LIVTALV(35)  
        LEAIO   => LIVTALV(36)   
        LEAPHT  => LIVTALV(37)  
        LEAPL   => LIVTALV(38) 
        LEMEL   => LIVTALV(39)  
        LEMAT   => LIVTALV(40)  
        LEMML   => LIVTALV(41)  
        LEMIO   => LIVTALV(42)   
        LEMPHT  => LIVTALV(43)  
        LEMPL   => LIVTALV(44) 
        LEIEL   => LIVTALV(45)  
        LEIAT   => LIVTALV(46)  
        LEIML   => LIVTALV(47)  
        LEIIO   => LIVTALV(48)   
        LEIPHT  => LIVTALV(49)  
        LEIPL   => LIVTALV(50)
        LEPHEL  => LIVTALV(51) 
        LEPHAT  => LIVTALV(52) 
        LEPHML  => LIVTALV(53) 
        LEPHIO  => LIVTALV(54)  
        LEPHPHT => LIVTALV(55) 
        LEPHPL  => LIVTALV(56) 
        LADDV   => LIVTALV(57)  
        LCOLV   => LIVTALV(58)  
        LSNAPV  => LIVTALV(59)
        LCOPV   => LIVTALV(60)  
        LBGKV   => LIVTALV(61)  
        LALGV   => LIVTALV(62)
        LPGENA  => LIVTALV(63) 
        LPGENM  => LIVTALV(64) 
        LPGENI  => LIVTALV(65) 
        LPGENPH => LIVTALV(66)
        LEGENA  => LIVTALV(67) 
        LEGENM  => LIVTALV(68) 
        LEGENI  => LIVTALV(69) 
        LEGENPH => LIVTALV(70)
        LVGENA  => LIVTALV(71) 
        LVGENM  => LIVTALV(72) 
        LVGENI  => LIVTALV(73) 
        LVGENPH => LIVTALV(74)
        LPPAT   => LIVTALV(75)  
        LPPML   => LIVTALV(76)  
        LPPIO   => LIVTALV(77)  
        LPPPHT  => LIVTALV(78)
        LPPPL   => LIVTALV(79)
        LEPAT   => LIVTALV(80)  
        LEPML   => LIVTALV(81)  
        LEPIO   => LIVTALV(82)  
        LEPPHT  => LIVTALV(83)
        LEPPL   => LIVTALV(84)
        
! surface averaged tallies

        LPOTAT    => LIVTALS(1) 
        LPRFAAT   => LIVTALS(2) 
        LPRFMAT   => LIVTALS(3) 
        LPRFIAT   => LIVTALS(4) 
        LPRFPHAT  => LIVTALS(5)
        LPRFPAT   => LIVTALS(6)
        LPOTML    => LIVTALS(7)
        LPRFAML   => LIVTALS(8) 
        LPRFMML   => LIVTALS(9) 
        LPRFIML   => LIVTALS(10) 
        LPRFPHML  => LIVTALS(11)
        LPRFPML   => LIVTALS(12)
        LPOTIO    => LIVTALS(13)
        LPRFAIO   => LIVTALS(14) 
        LPRFMIO   => LIVTALS(15) 
        LPRFIIO   => LIVTALS(16) 
        LPRFPHIO  => LIVTALS(17)
        LPRFPIO   => LIVTALS(18)
        LPOTPHT   => LIVTALS(19)
        LPRFAPHT  => LIVTALS(20) 
        LPRFMPHT  => LIVTALS(21) 
        LPRFIPHT  => LIVTALS(22) 
        LPRFPHPHT => LIVTALS(23)
        LPRFPPHT  => LIVTALS(24) 
        LPOTPL    => LIVTALS(25)
        LEOTAT    => LIVTALS(26)
        LERFAAT   => LIVTALS(27) 
        LERFMAT   => LIVTALS(28) 
        LERFIAT   => LIVTALS(29) 
        LERFPHAT  => LIVTALS(30) 
        LERFPAT   => LIVTALS(31)
        LEOTML    => LIVTALS(32)
        LERFAML   => LIVTALS(33) 
        LERFMML   => LIVTALS(34) 
        LERFIML   => LIVTALS(35) 
        LERFPHML  => LIVTALS(36) 
        LERFPML   => LIVTALS(37)
        LEOTIO    => LIVTALS(38)
        LERFAIO   => LIVTALS(39) 
        LERFMIO   => LIVTALS(40) 
        LERFIIO   => LIVTALS(41) 
        LERFPHIO  => LIVTALS(42) 
        LERFPIO   => LIVTALS(43)
        LEOTPHT   => LIVTALS(44)
        LERFAPHT  => LIVTALS(45) 
        LERFMPHT  => LIVTALS(46) 
        LERFIPHT  => LIVTALS(47) 
        LERFPHPHT => LIVTALS(48) 
        LERFPPHT  => LIVTALS(49) 
        LEOTPL    => LIVTALS(50)
        LSPTAT    => LIVTALS(51) 
        LSPTML    => LIVTALS(52) 
        LSPTIO    => LIVTALS(53) 
        LSPTPHT   => LIVTALS(54) 
        LSPTPL    => LIVTALS(55)
        LSPTTOT   => LIVTALS(56)
        LADDS     => LIVTALS(57)  
        LALGS     => LIVTALS(58)
        LSPUMP    => LIVTALS(59)

        NFIRST = 0
        NADDV  = 0
        IRESC1 = 0
        IRESC2 = 0
        NFRSTW = 0
        NADDW  = 0
              
      ELSE IF (ICAL == 2) THEN

        ESTIMV = 0._DP
        ESTIMS = 0._DP

        CEMETERYV = 0._DP
        CEMETERYS = 0._DP

      END IF

      RETURN
      END SUBROUTINE INIT_CESTIM

      
      SUBROUTINE SPEC_TO_SPEC (SPECA, SPECB)

      TYPE(EIRENE_SPECTRUM), INTENT(OUT) :: SPECA
      TYPE(EIRENE_SPECTRUM), INTENT(IN) :: SPECB
      
      SPECA%SPCMIN  = SPECB%SPCMIN
      SPECA%SPCMAX  = SPECB%SPCMAX
      SPECA%SPCDEL  = SPECB%SPCDEL
      SPECA%SPCDELI = SPECB%SPCDELI 
      SPECA%SPCINT  = SPECB%SPCINT
      SPECA%SGMS    = SPECB%SGMS
      SPECA%STVS    = SPECB%STVS
      SPECA%EES     = SPECB%EES 
      SPECA%NSPC    = SPECB%NSPC
      SPECA%ISPCTYP = SPECB%ISPCTYP
      SPECA%ISPCSRF = SPECB%ISPCSRF
      SPECA%IPRTYP  = SPECB%IPRTYP
      SPECA%IPRSP   = SPECB%IPRSP
      SPECA%IMETSP  = SPECB%IMETSP
      SPECA%SPC     = SPECB%SPC
      IF (ASSOCIATED(SPECA%SDV)) THEN
        SPECA%SDV     = SPECB%SDV
        SPECA%SGM     = SPECB%SGM
      END IF
      END SUBROUTINE SPEC_TO_SPEC

      END MODULE CESTIM
