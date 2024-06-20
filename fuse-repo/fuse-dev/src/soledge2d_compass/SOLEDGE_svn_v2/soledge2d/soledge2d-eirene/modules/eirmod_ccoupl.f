      MODULE EIRMOD_CCOUPL
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_CCOUPL, EIRENE_DEALLOC_CCOUPL, 
     P          EIRENE_INIT_CCOUPL
 
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCCPL(:)
 
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R FCTE(:), D(:),  FL(:), BMASS(:), XMCP_OLD(:),
     R CHGP,    CHGEE, CHGEI, CHGMOM
 
      REAL(DP), PUBLIC, SAVE ::
     R B2BREM,  B2RAD, B2QIE, B2VDP
 
      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         ICCPL1(:,:), ICCPL2(:)
 
      INTEGER, PUBLIC, POINTER, SAVE ::
     I  NDT(:,:),  NINCT(:,:), NIXY(:,:),
     I  NTIN(:,:), NTEN(:,:),  NIFLG(:,:),
     I  NPTC(:,:), NSPZI(:,:), NSPZE(:,:),
     I  NEMOD(:,:)
 
      INTEGER, PUBLIC, POINTER, SAVE ::
     I NTGPRT(:), IFLB(:), NAOTS(:), NAOTT(:),
     I NTARGI, NSTRI,  NFLA,   NCUTB,  NCUTL,  NDXA,   NDYA,
     I NCLMI,  NBLCKI, NPRNVI, NPRTVI, NPRDVI,
     I NMODEI, NFILNN, NCUTB_SAVE,
     I NAINB,  NAOTB
 
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NAINS(:), NAINT(:)
 
      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE :: LCCPL(:)
 
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L LBALAN, LSYMET, LPRSOU,
     L LNLPLG, LNLDRF, LTRCFL, LNLVOL(:)
 
      INTEGER, PUBLIC, SAVE ::
     I NCOUPL, MCOUPL1, MCOUPL2, LCOUPL
 
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_CCOUPL (ICAL)
 
      INTEGER, INTENT(IN) :: ICAL
 
      IF (ICAL == 1) THEN
 
        IF (ALLOCATED(RCCPL)) RETURN
 
        NCOUPL  = 4*NPLS+4+NSTRA
        MCOUPL1 = 10*NSTEP*NPTRGT
        MCOUPL2 = 1*NPLS+NSTEP+17+2*NLIMPS
        LCOUPL  = 6+NSTRA
 
        ALLOCATE (RCCPL(NCOUPL))
        ALLOCATE (ICCPL1(10*NSTEP,NPTRGT))
        ALLOCATE (ICCPL2(MCOUPL2))
        ALLOCATE (LCCPL(LCOUPL))
 
        WRITE (55,'(A,T25,I15)')
     .        ' CCOUPL ',NCOUPL*8 + (10*NSTEP*NPTRGT+MCOUPL2)*4
     .                   + LCOUPL*4
 
        FCTE     => RCCPL(1+0*NPLS : 1*NPLS)
        D        => RCCPL(1+1*NPLS : 2*NPLS)
        FL       => RCCPL(1+2*NPLS : 3*NPLS)
        BMASS    => RCCPL(1+3*NPLS : 4*NPLS)
        XMCP_OLD => RCCPL(1+4*NPLS : 4*NPLS+NSTRA)
        CHGP     => RCCPL(1+4*NPLS+NSTRA)
        CHGEE    => RCCPL(2+4*NPLS+NSTRA)
        CHGEI    => RCCPL(3+4*NPLS+NSTRA)
        CHGMOM   => RCCPL(4+4*NPLS+NSTRA)
 
        NDT    => ICCPL1(1+0*NSTEP : 1*NSTEP,:)
        NINCT  => ICCPL1(1+1*NSTEP : 2*NSTEP,:)
        NIXY   => ICCPL1(1+2*NSTEP : 3*NSTEP,:)
        NTIN   => ICCPL1(1+3*NSTEP : 4*NSTEP,:)
        NTEN   => ICCPL1(1+4*NSTEP : 5*NSTEP,:)
        NIFLG  => ICCPL1(1+5*NSTEP : 6*NSTEP,:)
        NPTC   => ICCPL1(1+6*NSTEP : 7*NSTEP,:)
        NSPZI  => ICCPL1(1+7*NSTEP : 8*NSTEP,:)
        NSPZE  => ICCPL1(1+8*NSTEP : 9*NSTEP,:)
        NEMOD  => ICCPL1(1+9*NSTEP :10*NSTEP,:)
 
        NTARGI     => ICCPL2( 1)
        NSTRI      => ICCPL2( 2)
        NFLA       => ICCPL2( 3)
        NCUTB      => ICCPL2( 4)
        NCUTL      => ICCPL2( 5)
        NDXA       => ICCPL2( 6)
        NDYA       => ICCPL2( 7)
        NCLMI      => ICCPL2( 8)
        NBLCKI     => ICCPL2( 9)
        NPRNVI     => ICCPL2(10)
        NPRTVI     => ICCPL2(11)
        NPRDVI     => ICCPL2(12)
        NMODEI     => ICCPL2(13)
        NFILNN     => ICCPL2(14)
        NCUTB_SAVE => ICCPL2(15)
        NAINB      => ICCPL2(16)
        NAOTB      => ICCPL2(17)
        NTGPRT     => ICCPL2(18 : 17+NSTEP)
        IFLB       => ICCPL2(18+NSTEP : 17+NSTEP+NPLS)
        NAOTS      => ICCPL2(18+NSTEP+NPLS :
     .                       17+NSTEP+NPLS+NLIMPS)
        NAOTT      => ICCPL2(18+NSTEP+NPLS+NLIMPS :
     .                       17+NSTEP+NPLS+2*NLIMPS)
 
        LBALAN => LCCPL(1)
        LSYMET => LCCPL(2)
        LPRSOU => LCCPL(3)
        LNLPLG => LCCPL(4)
        LNLDRF => LCCPL(5)
        LTRCFL => LCCPL(6)
        LNLVOL => LCCPL(7:6+NSTRA)
 
      ELSE IF (ICAL == 2) THEN
 
        IF (ALLOCATED(NAINS)) RETURN
 
        ALLOCATE (NAINS(NAIN))
        ALLOCATE (NAINT(NAIN))
 
      END IF
 
      CALL EIRENE_INIT_CCOUPL (ICAL)
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_CCOUPL
 
 
      SUBROUTINE EIRENE_DEALLOC_CCOUPL
 
      IF (.NOT.ALLOCATED(RCCPL)) RETURN
 
      DEALLOCATE (RCCPL)
      DEALLOCATE (ICCPL1)
      DEALLOCATE (ICCPL2)
      DEALLOCATE (LCCPL)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CCOUPL
 
 
      SUBROUTINE EIRENE_INIT_CCOUPL(ICAL)
 
      INTEGER, INTENT(IN) :: ICAL
 
      IF (ICAL == 1) THEN
 
        RCCPL  = 0._DP
        B2BREM = 0._DP
        B2RAD  = 0._DP
        B2QIE  = 0._DP
        B2VDP  = 0._DP
        ICCPL1 = 0
        ICCPL2 = 0
        LCCPL  = .FALSE.
 
      ELSE IF (ICAL == 2) THEN
 
        NAINS = 0
        NAINT = 0
 
      END IF
 
      RETURN
      END SUBROUTINE EIRENE_INIT_CCOUPL
 
      END MODULE EIRMOD_CCOUPL
