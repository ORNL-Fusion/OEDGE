      SUBROUTINE EIRENE_READ_TRIANG (CASENAME)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CTRIG
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      CHARACTER*(*), INTENT(IN) :: CASENAME
      CHARACTER(100) :: FILENAME, ZEILE
      INTEGER :: LL, I, IND, IT, IS, IS1, IER, NRK, ISTS,
     .           ISTMIN, ISTMAX, IC, n2, n3, nb1, ns1, nb3, ns3,
     .           jt, js, js1, imin, imax, jmin, jmax, in1, in3
      real(dp) :: ar
 
      LL=LEN_TRIM(CASENAME)
 
      	FILENAME=CASENAME(1:LL) // '.npco_char'
      	OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')
 
      	ZEILE='*   '
      	DO WHILE (ZEILE(1:1) == '*')
         	READ (30,'(A100)') ZEILE
      	END DO
 
      	READ (ZEILE,*) NRK
 
      	IF (NRK /= NRKNOT) THEN
        	WRITE (iunout,*) ' NRKNOT IS WRONG IN EIRENE INPUT FILE'
        	WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        	CALL EIRENE_EXIT_OWN(1)
      	END IF
 
      	DO I=1,NRKNOT
        	READ(30,*) IND, XTRIAN(I), YTRIAN(I)
      	END DO
 
!pb      XTRIAN(1:NRKNOT) = XTRIAN(1:NRKNOT) * 100._DP
!pb      YTRIAN(1:NRKNOT) = YTRIAN(1:NRKNOT) * 100._DP
 
      	CLOSE (UNIT=30)
 
      	FILENAME=CASENAME(1:LL) // '.elemente'
      	OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')
 
      	ZEILE='*   '
      	DO WHILE (ZEILE(1:1) == '*')
         	READ (30,'(A100)') ZEILE
      	END DO
 
      	READ (ZEILE,*) NTRII
 
      	IF (NTRII > NR1ST) THEN
        	WRITE (iunout,*) ' NR1ST IS WRONG IN EIRENE INPUT FILE'
        	WRITE (iunout,*) ' CHECK FOR CORRECT NUMBER IN FILE ',FILENAME
        	CALL EIRENE_EXIT_OWN(1)
      	END IF
 
      	DO I=1,NTRII
        	READ (30,*) IND, NECKE(1,I), NECKE(2,I), NECKE(3,I)
     	END DO
 
      	CLOSE (UNIT=30)
 
      	FILENAME=CASENAME(1:LL) // '.neighbors'
      	OPEN (UNIT=30,FILE=FILENAME,ACCESS='SEQUENTIAL',FORM='FORMATTED')
 
      	ZEILE='*   '
      	DO WHILE (ZEILE(1:1) == '*')
         READ (30,'(A100)') ZEILE
      	END DO
 
      	DO I=1,NTRII
        	READ (30,*) IND, NCHBAR(1,I), NSEITE(1,I), INMTI(1,I),
     .                   NCHBAR(2,I), NSEITE(2,I), INMTI(2,I),
     .                   NCHBAR(3,I), NSEITE(3,I), INMTI(3,I),
     .                   IXTRI(I),    IYTRI(I)
!       IF (INMTI(1,I) /= 0) INMTI(1,I) = ABS(INMTI(1,I)) + NLIM
!       IF (INMTI(2,I) /= 0) INMTI(2,I) = ABS(INMTI(2,I)) + NLIM
!       IF (INMTI(3,I) /= 0) INMTI(3,I) = ABS(INMTI(3,I)) + NLIM
        	IF (INMTI(1,I) /= 0) INMTI(1,I) = INMTI(1,I) + NLIM
        	IF (INMTI(2,I) /= 0) INMTI(2,I) = INMTI(2,I) + NLIM
        	IF (INMTI(3,I) /= 0) INMTI(3,I) = INMTI(3,I) + NLIM

      	END DO
 
      	CLOSE (UNIT=30)

 
      IER = 0
      IF ((MAXVAL(NECKE(1:3,1:NTRII)) > NRKNOT) .OR.
     .    (MINVAL(NECKE(1:3,1:NTRII)) <= 0 )) THEN
        WRITE (iunout,*) ' WRONG COORDINATE NUMBER IS DEFINITION OF',
     .              ' TRIANGLES FOUND '
        IER = 2
      END IF
 
      DO IT=1,NTRII
        DO IS=1,3
          IF (NCHBAR(IS,IT).EQ.0.AND.INMTI(IS,IT).EQ.0) THEN
            WRITE (iunout,*) ' ERROR IN READ_TRIANG '
            WRITE (iunout,*) ' OPEN SIDE OF TRIANGLE ',IT,' SIDE ',IS
            IS1=IS+1
            IF (IS.EQ.3) IS1=1
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS,IT)),
     .                                    YTRIAN(NECKE(IS,IT))
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS1,IT)),
     .                                    YTRIAN(NECKE(IS1,IT))
            IER = 3
          ENDIF
        ENDDO
      ENDDO
 
! check and correct the orientation
 
      do it=1,ntrii
         AR=0.5*(XTRIAN(NECKE(2,IT))*(YTRIAN(NECKE(3,IT))
     >          -YTRIAN(NECKE(1,IT)))+XTRIAN(NECKE(3,IT))*
     >          (YTRIAN(NECKE(1,IT))
     >          -YTRIAN(NECKE(2,IT)))+XTRIAN(NECKE(1,IT))*
     >          (YTRIAN(NECKE(2,IT))-YTRIAN(NECKE(3,IT))))
         if (ar < 0._dp) then
           n2 = necke(2,it)
           n3 = necke(3,it)
           necke(2,it) = n3
           necke(3,it) = n2
 
           nb1 = nchbar(1,it)
           ns1 = nseite(1,it)
           in1 = inmti(1,it)
           nb3 = nchbar(3,it)
           ns3 = nseite(3,it)
           in3 = inmti(3,it)
           nchbar(1,it) = nb3
           nseite(1,it) = ns3
           inmti(1,it) = in3
           nchbar(3,it) = nb1
           nseite(3,it) = ns1
           inmti(3,it) = in1
 
           if (nb3 > 0) nseite(ns3,nb3) = 1
           if (nb1 > 0) nseite(ns1,nb1) = 3
         end if
      end do
 
      do it=1,ntrii
        do is=1,3
          if (nchbar(is,it) == 0) cycle
          is1 = is+1
          if (is1 > 3) is1 = 1
          imin = min(necke(is,it),necke(is1,it))
          imax = max(necke(is,it),necke(is1,it))
          jt = nchbar(is,it)
          js = nseite(is,it)
          js1 = js+1
          if (js1 > 3) js1 = 1
          jmin = min(necke(js,jt),necke(js1,jt))
          jmax = max(necke(js,jt),necke(js1,jt))
          if ((imin /= jmin) .or. (imax /= jmax)) then
            write (iunout,*) ' neighbor triangles do not match '
            write (iunout,*) ' it = ',it, ' is = ',is
            write (iunout,*) ' jt = ',jt, ' js = ',js
            write (iunout,*) ' imin = ',imin, ' imax = ',imax
            write (iunout,*) ' jmin = ',jmin, ' jmax = ',jmax
            ier = 4
          end if
        end do
      end do
 
      IF (IER /= 0) CALL EIRENE_EXIT_OWN(1)
 
      ISTMIN=MINVAL(INMTI(1:3,1:NTRII),MASK=(INMTI(1:3,1:NTRII)/=0))
      ISTMAX=MAXVAL(INMTI(1:3,1:NTRII),MASK=(INMTI(1:3,1:NTRII)/=0))
 
      IC=0
      DO ISTS=ISTMIN,ISTMAX
        DO IS=1,3
          DO IT=1,NTRII
            IF (INMTI(IS,IT) == ISTS) THEN
              IC=IC+1
              INSPAT(IS,IT)=IC
            END IF
          END DO
        END DO
      END DO
 
      RETURN
      END SUBROUTINE EIRENE_READ_TRIANG
