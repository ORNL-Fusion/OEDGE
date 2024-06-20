      SUBROUTINE EIRENE_CARBON_LINES (IST,IAD1,IAD2,IAD3,IAD4,IAD5,IADS)
c
C  SUBROUTINE FOR CARBON LINES EMISSIVITY.
C  CALLED FROM EIRENE, SECTION DIAGNO, SUBR. SIGAL
C  THE EMISSIVITY PROFILE (PHOTONS/S/CM**3) IS COMPUTED
C  AND WRITTEN ONTO TALLIES ADDV(IAD1,...),... FOR STRATUM NO. IST
C  PEC UNIT IS PHOTONS.CM**3/S, TO BE MULTIPLIED BY N_C+.Ne 
C  IAD1: LINE @ 6581A LINEAR IN C+ DENSITY (NO METASTABLE)
C  IAD2: 
C  IAD3: 
C  IAD4: 
C  IAD5: 
C  IADS: 
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CUPD
      USE EIRMOD_COMSIG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
 
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: IAD1, IAD2, IAD3, IAD4, IAD5, IADS, IST

      INTEGER :: IFIRST, NCELC, IERROR, IR, JFEXMN, JFEXMX, ntemp,ICON
      CHARACTER(8) :: FILNAM
      CHARACTER(4) :: H123
      CHARACTER(9) :: REAC
      CHARACTER(3) :: CRC
      CHARACTER(6) :: CISTRA

      TYPE(ADAS_DATA) , POINTER :: PEC1
 
      real(dp) :: TE, DE, DI, PEC_intp, RCMIN, RCMAX, FP(6)
C
      SAVE
C
      DATA IFIRST/0/

      interface
        function EIRENE_intp_adas (ad,p1,p2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function EIRENE_intp_adas
      end interface

C  INITIALIZE ATOMIC DATA ARRAYS
C
      IF (IFIRST.EQ.0) THEN

      	write(iunout,*) 'CARBON LINES CALLED, NOT BALMER ALPHA !'
        IFIRST=1
        IERROR=0
C
C  READ REDUCED POPULATION COEFFICIENT FOR HYDR. ATOMS FROM FILE AMJUEL
C  AND PUT THEM FROM CREAC(..,..,IR) ONTO DA,DPP,DM,DI, AND DN ARRAY
C
        IR=NREACI+1
        IF (IR.GT.NREAC) THEN
          WRITE (iunout,*) 'FROM SUBROUTINE EIRENE_Carbon_lines: '
          CALL EIRENE_MASPRM('NREAC',5,NREAC,'IR',2,IR,IERROR)
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      END IF


      IF (MAX(IAD1,IAD2,IAD3,IAD4,IAD5,IADS) > NADV) GOTO 999
      ADDV(IAD1,1:NRTAL) = 0.D0
      ADDV(IAD2,1:NRTAL) = 0.D0
      ADDV(IAD3,1:NRTAL) = 0.D0
      ADDV(IAD4,1:NRTAL) = 0.D0
      ADDV(IAD5,1:NRTAL) = 0.D0
 
      ADDV(IADS,1:NRTAL) = 0.D0

C  LINES @ 6581 A & 4268 A

      REAC='PEC96   '
      FILNAM='ADAS    '
      H123='H.4 '
      CRC='OT '
      FP = 0._dp
      RCMIN = -20._dp
      RCMAX =  20._dp
      JFEXMN = 0
      JFEXMX = 0

C the last number specifies the data block (i.e. the line) in the ADAS data file

      DO ICON=1,6     

c allows reading several times reaction NREACI+1
      	REACDAT(NREACI+1)%LRTC = .FALSE.


      	CALL EIRENE_SLREAC(NREACI+1,FILNAM,H123,REAC,CRC, 
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'c1 ',ICON)

      	PEC1 => reacdat(NREACI+1)%rtc%adas

        PEC1%fit = TRANSPOSE(PEC1%fit)

  	Ntemp=reacdat(NREACI+1)%rtc%adas%Ntemp
C
      		DO 1000 NCELL=1,NSBOX
C
C  LOCAL BACKGROUND DATA ARE IN CELL NCELL
C  LOCAL TEST PARTICLE DATA ARE IN CELL NCELC
C
        	NCELC=NCLTAL(NCELL)
C
        	IF (NSTGRD(NCELL) > 0) CYCLE
 
        	IF (LGVAC(NCELL,NPLS+1)) THEN
          		TE=TVAC
          		DE=DVAC
      			DI=DVAC
        	ELSE
          		TE=TEIN(NCELL)
          		DE=DEIN(NCELL)
c H+ density for CX contribution ?
      			DI=DIIN(1,NCELL)
        	ENDIF

C  INTERPOLATE TO FIND THE PEC VALUE IN THE CELL
        
        	if (TE > reacdat(NREACI+1)%rtc%adas%temp(Ntemp) .or.
     .              TE < reacdat(NREACI+1)%rtc%adas%temp(1)) then
        
        		PEC_intp=0._dp
        	else

        		PEC_intp = eirene_intp_adas(PEC1,TE,DE)
        	end if
        
        	if (PEC_intp < 0._dp) then
        		write(iunout,*) 'Negative value for PEC !'
                	call eirene_exit_own(1)
        	end if        

                select case (ICON)        	
			case(1)
     				ADDV(IAD1,NCELC)=ADDV(IAD1,NCELC)
     . 		+PDENI(3,NCELC)*DE*PEC_intp*vol(ncell)
      			case(2)
   				ADDV(IAD2,NCELC)=ADDV(IAD2,NCELC)
     .		+PDENI(4,NCELC)*DE*PEC_intp
      			case(3)
c CX refers to H + C++ -> H+ +C+(n) (see ADAS manual)
                        	ADDV(IAD3,NCELC)=ADDV(IAD3,NCELC)
     .		+PDENA(1,NCELC)*PDENI(4,NCELC)*PEC_intp*vol(ncell)
      			case(4)
     				ADDV(IAD4,NCELC)=ADDV(IAD4,NCELC)
     . 		+PDENI(3,NCELC)*DE*PEC_intp*vol(ncell)
      			case(5)
   				ADDV(IAD5,NCELC)=ADDV(IAD5,NCELC)
     .		+PDENI(4,NCELC)*DE*PEC_intp*vol(ncell)
      			case(6)
                        	ADDV(IADS,NCELC)=ADDV(IADS,NCELC)
     .		+PDENA(1,NCELC)*PDENI(4,NCELC)*PEC_intp*vol(ncell)
     		end select
1000  		CONTINUE
      END DO

c      ADDV(IADS,1:NSBOX_TAL)=ADDV(IAD1,1:NSBOX_TAL)+
c     .    ADDV(IAD2,1:NSBOX_TAL)+ADDV(IAD3,1:NSBOX_TAL)

      ADDV(IAD1,1:NSBOX_TAL)=ADDV(IAD1,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD2,1:NSBOX_TAL)=ADDV(IAD2,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD3,1:NSBOX_TAL)=ADDV(IAD3,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD4,1:NSBOX_TAL)=ADDV(IAD4,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD5,1:NSBOX_TAL)=ADDV(IAD5,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
 
      ADDV(IADS,1:NSBOX_TAL)=ADDV(IADS,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
 
      CALL EIRENE_LEER(2)
      CALL EIRENE_FTCRI(IST,CISTRA)
      IF (IST.GT.0) CALL EIRENE_MASBOX
     .   ('SUBR. BA_ALPHA CALLED, FOR STRATUM NO. '//CISTRA)
      IF (IST.EQ.0)
     .CALL EIRENE_MASBOX('SUBR. BA_ALPHA CALLED, FOR SUM OVER STRATA')
      CALL EIRENE_LEER(1)

      CALL EIRENE_LEER(2)
 
      CALL EIRENE_INTTAL
     .  (ADDV,VOLTAL,IAD1,NADV,NSBOX_TAL,ADDVI(IAD1,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL EIRENE_INTTAL
     .  (ADDV,VOLTAL,IAD2,NADV,NSBOX_TAL,ADDVI(IAD2,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL EIRENE_INTTAL
     .  (ADDV,VOLTAL,IAD3,NADV,NSBOX_TAL,ADDVI(IAD3,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL EIRENE_INTTAL
     .  (ADDV,VOLTAL,IAD4,NADV,NSBOX_TAL,ADDVI(IAD4,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
      CALL EIRENE_INTTAL
     .  (ADDV,VOLTAL,IAD5,NADV,NSBOX_TAL,ADDVI(IAD5,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
 
 
      CALL EIRENE_INTTAL
     .  (ADDV,VOLTAL,IADS,NADV,NSBOX_TAL,ADDVI(IADS,IST),
     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)

      RETURN

999   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. Carbon lines '
      WRITE (iunout,*) 'NO STORAGE AVAILBALE ON ADDITIONAL TALLY ADDV '
      WRITE (iunout,*) 'STORAGE REQUESTED FOR IADV= ',
     .             IAD1,IAD2,IAD3,IAD4,
     .             IAD5,IADS
      WRITE (iunout,*) 'CHECK INPUT BLOCK 10A '
      CALL EIRENE_EXIT_OWN(1)













      END
