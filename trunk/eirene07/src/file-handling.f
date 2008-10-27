C EIRENE07 COMPILATION
C ===== SOURCE: rdtrim.f
C
C
      SUBROUTINE RDTRIM
C
C  THIS SUBROUTINE READS SELECTIVELY SOME
C  REFLECTION DATA PRODUCED BY MONTE CARLO CODES
C
      USE PRECISION
      USE PARMMOD
      USE CREFMOD
      USE COMPRT, ONLY: IUNOUT
      USE CREF
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PID180
      INTEGER :: I1, I2, I3, I4, I5, IUN, IFILE, I, IWWW
C
C
      INE=12
      INW=7
      INR=5
C
      IF (INE*INW*NFLR.GT.NH0 .OR.
     .    INE*INW*INR*NFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*NFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*NFLR.GT.NH3) THEN
        WRITE (iunout,*) 
     .    'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EXIT_OWN(1)
      ENDIF
C
      DO 7 IFILE=1,NFLR
        IUN=20
        OPEN (UNIT=IUN,FILE=REFFIL(IFILE),ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')
C
        READ (IUN,*)
        READ (IUN,*)
        READ (IUN,*)
        READ (IUN,*)
        DO 2 I1=1,INE
          DO 3 I2=1,INW
            READ (IUN,*)
            READ (IUN,*)
            READ (IUN,*) TC(IFILE),TM(IFILE),WC(IFILE),WM(IFILE),
     .                   enar(i1),wiar(i2),HFTR0(I1,I2,IFILE)
C  FIND NEAREST INTEGER FOR NUCLEAR MASS NUMBER
            IWWW=NINT(WM(IFILE))
            WM(IFILE)=IWWW
            READ (IUN,*)
            READ (IUN,*) (HFTR1(I1,I2,I3,IFILE),I3=1,INR)
            READ (IUN,*)
            DO 5 I3=1,INR
              READ (IUN,*) (HFTR2(I1,I2,I3,I4,IFILE),I4=1,INR)
5           CONTINUE
            READ (IUN,*)
            DO 6 I3=1,INR
            DO 6 I4=1,INR
              READ (IUN,*) (HFTR3(I1,I2,I3,I4,I5,IFILE),I5=1,INR)
6           CONTINUE
3         CONTINUE
2       CONTINUE
        CLOSE (UNIT=IUN)
7     CONTINUE
C
      INEM=INE-1
      DO 11 I=1,INEM
11      DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      PID180=ATAN(1.)/45.
      DO 12 I=1,INW
12      WIAR(I)=COS(WIAR(I)*PID180)
      INWM=INW-1
      DO 13 I=1,INWM
13      DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
      INRM=INR-1
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
      DO 15 I=1,INRM
15      DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
C
      RETURN
      END
C ===== SOURCE: read_adas.f
!pb  21.11.06: index error corrected in defintion of ap%dte

      subroutine read_adas (ir,reac,isw,iz1)

      use precision
      use parmmod
      use comxs
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, isw, iz1
      character(len=*), intent(in) :: reac
      integer :: nz, nde, nte, iza, ize, io, lc, ind, ian, ien,
     .           ide, ite, iz
      character(132) :: zeile
      type(adas_data), pointer :: ap

      read (29,*,iostat=io) nz, nde, nte, iza, ize

      if (io .ne. 0) then
        write (iunout,*) ' ERROR READING FILE FROM ADAS DATABASE '
        write (iunout,*) ' DIRECTORY IS ',reac
        call exit_own(1)
      end if

      if ((iz1 < iza) .or. (iz1 > ize)) then
        write (iunout,*) ' ERROR READING FILE FROM ADAS DATABASE '
        write (iunout,*) ' REQUESTED Z1 IS NOT AVAILABLE '
        write (iunout,*) ' Z1, ZA, ZE ',IZ1, IZA, IZE
        call exit_own(1)
      end if
      
      allocate (ap)
      allocate (ap%dens(nde))
      allocate (ap%temp(nte))
      allocate (ap%dde(nde))
      allocate (ap%dte(nte))
      allocate (ap%fit(nte,nde))

      ap%ndens = nde
      ap%ntemp = nte

      read (29,*)

      lc = len_trim(reac)
      if (reac(lc:lc) == 'r') then
        read (29,*)
        read (29,*)
      end if

! read densities
      read (29,*) (ap%dens(ide),ide=1,nde)

! read temperatures
      read (29,*) (ap%temp(ite),ite=1,nte)
      
! find appropriate Z1-block

      do

! read line between data blocks
        read (29,'(A132)') zeile
        if (zeile(1:5) == '-----') then
          ind = index(zeile,'Z1')
          if (ind == 0) cycle
          ian = ind + scan(zeile(ind+1:),' ')
          ien = ian + scan(zeile(ian+1:),'/') - 1
          read (zeile(ian:ien),*) iz
          if (iz == iz1) exit
        end if  
 
      end do

      do ite = 1, nte
        read (29,*) (ap%fit(ite,ide), ide = 1,nde)
      end do

      close (29)

! set up differenz arrays for densities and temperatures
      
      do ide=1,nde-1
        ap%dde(ide) = 1._dp / (ap%dens(ide+1) - ap%dens(ide))
      end do
      
      do ite=1,nte-1
        ap%dte(ite) = 1._dp / (ap%temp(ite+1) - ap%temp(ite))
      end do


      select case (isw)

      case (0)
        IF (REACDAT(IR)%LPOT) THEN
          WRITE (IUNOUT,*) ' POTENTIAL ALREADY SPECIFIED FOR REACTION',
     .                       IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lpot = .true.
        allocate (reacdat(ir)%pot)
        reacdat(ir)%pot%adas => ap
        reacdat(ir)%pot%ifit = 3

      case (1)
        IF (REACDAT(IR)%LCRS) THEN
          WRITE (IUNOUT,*) ' CROSS SECTION ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lcrs = .true.
        allocate (reacdat(ir)%crs)
        reacdat(ir)%crs%adas => ap
        reacdat(ir)%crs%ifit = 3

      case (2:4)
        IF (REACDAT(IR)%LRTC) THEN
          WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtc = .true.
        allocate (reacdat(ir)%rtc)
        reacdat(ir)%rtc%adas => ap
        reacdat(ir)%rtc%ifit = 3

      case (5:7)
        IF (REACDAT(IR)%LRTCMW) THEN
          WRITE (IUNOUT,*) ' MOMEMTUM WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtcmw = .true.
        allocate (reacdat(ir)%rtcmw)
        reacdat(ir)%rtcmw%adas => ap
        reacdat(ir)%rtcmw%ifit = 3

      case (8:10)
        IF (REACDAT(IR)%LRTCEW) THEN
          WRITE (IUNOUT,*) ' ENERGY WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtcew = .true.
        allocate (reacdat(ir)%rtcew)
        reacdat(ir)%rtcew%adas => ap
        reacdat(ir)%rtcew%ifit = 3

      case (11:12)
        IF (REACDAT(IR)%LOTH) THEN
          WRITE (IUNOUT,*) ' OTHER POLYNOMIAL FIT COEFFICIENTS',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF
        reacdat(ir)%loth = .true.
        allocate (reacdat(ir)%oth)
        reacdat(ir)%oth%adas => ap
        reacdat(ir)%oth%ifit = 3

      case default
        WRITE (IUNOUT,*) ' WRONG REACTION TYPE SPCIFIED '
        WRITE (IUNOUT,*) ' REACTION NO. ', IR
        WRITE (IUNOUT,*) ' REACTION TYPE H.', ISW
        CALL EXIT_OWN(1)
      end select
      
      return
      end subroutine read_adas
      
C ===== SOURCE: read_photdbk.f
      subroutine read_photdbk (ir, reac, isw)
!  16.2.05:  write statement taken out
!  2.11.05:  database handle introduced for file POLARI
      use precision
      use parmmod
      USE COMXS
      USE COMPRT
      USE CCONA
      USE CINIT
      USE PHOTON

      implicit none

      integer, intent(in) :: ir, isw
      CHARACTER(50), INTENT(IN) :: REAC

      real(dp) :: wl, aik, ei, ej, c2, c3, c4, c6_theo, b12, b21,
     .            c3_theo, c6_qs_mess, c6_ar_mess, c3_mess, c6
      real(dp) :: c6a(12), rdata(9,1)
      real(dp), save :: polari_fac(120)
      integer :: gi, gj
      integer :: ianf, iend, iblnk, lr, ic, i, j, iplsc3, iprftype,
     .           imess, ifremd, ii, i1, lel, nrjprt
      integer :: ik6, ipc6(12)
      integer, save :: ifrst=0, npolari
      character(1000) :: zeile
      character(20) :: elementname
      character(1) :: cha
      character(8) :: kenn, kennc3
      character(8) :: c6names(12)=(/'C6ARBORN','C6THEOR ','C6QS    ',
     .                              'C6AR    ','C6ZN    ',
     .                              'C6I     ','C6XE    ','C6DY    ',
     .                              'C6HO    ','C6TM    ','C6HG    ',
     .                              'C6TL    '/)
      character(2), save :: polari_elnam(120)
      type(line_data), pointer :: phline

      IF (REACDAT(IR)%LPHR) THEN
          WRITE (IUNOUT,*) ' PARAMETER FOR PHOTONIC REACTION ALREADY',
     .                     ' SPECIFIED FOR REACTION', IR
          WRITE (IUNOUT,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
          CALL EXIT_OWN(1)
        END IF

      if (ifrst == 0) then
         ifrst = 1
         call read_polari(polari_elnam,polari_fac,npolari)
      end if

      lr=len_trim(reac)

      read (29,*)
      read (29,*)
      read (29,*)

      do
        read (29,'(A1000)',end=990) zeile
        if (zeile(1:2) == '--') cycle
        call subcomma(zeile)

!  read Element
        ianf = 3
        iend = ianf + scan(zeile(ianf:),'|') - 1

        call delete_blanks(zeile(ianf:iend))
        iblnk = scan(zeile(ianf:iend),'|') - 1
        if (iblnk < 0) iblnk = iend-ianf+1
        if (iblnk == 0) then
           write (iunout,*) ' ERROR IN DATABASE PHOTON'
           write (iunout,*) ' NO ELEMENTNAME FOUND '
           stop
        end if

        elementname = repeat(' ',20)
        elementname(1:iblnk) = zeile(ianf:ianf+iblnk-1)

!  read Wellenlaenge
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) wl

!  read Zaehler
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        ic = ianf + verify(zeile(ianf:iend-1),' ') - 1
        cha = zeile(ic:ic)

        write (elementname(iblnk+1:),'(f10.4,a1)') wl,cha
        lel = len_trim(elementname)
        i1 = index(elementname,' ')
        do while (i1 < lel)
           elementname(i1:lel-1) = elementname(i1+1:lel)
           elementname(lel:lel) = ' '
           lel = lel - 1
           i1 = index(elementname,' ')
        end do

        if (elementname(1:iblnk+10) /= reac(1:iblnk+10)) cycle

!  skip Uebergang
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read aik
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) aik

!  skip fij
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read gj
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) gj

!  read gi
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) gi

!  skip ll
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1


!  skip lu
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read ej
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          ej = 0._dp
        else
          read (zeile(ianf:iend-1),*) ej
        end if

!  read ei
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          ei = 0._dp
        else
          read (zeile(ianf:iend-1),*) ei
        end if

!  read c2
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          c2 = 0._dp
        else
          read (zeile(ianf:iend-1),*) c2
        end if

!  read c3 (theo)
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          c3_theo = 0._dp
        else
          read (zeile(ianf:iend-1),*) c3_theo
        end if

!  read c4
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          c4 = 0._dp
        else
          read (zeile(ianf:iend-1),*) c4
        end if

!  read c6 (theo)
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          c6_theo = 0._dp
        else
          read (zeile(ianf:iend-1),*) c6_theo
        end if

!  read c6qs (mess)
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          c6_qs_mess = 0._dp
        else
          read (zeile(ianf:iend-1),*) c6_qs_mess
        end if

!  read c6 Ar (mess)
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          c6_ar_mess = 0._dp
        else
          read (zeile(ianf:iend-1),*) c6_ar_mess
        end if

!  read c3 (mess)
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          c3_mess = 0._dp
        else
          read (zeile(ianf:iend-1),*) c3_mess
        end if

        exit

      end do

      close (unit=29)

      read (iunin,'(12i6)') iprftype, iplsc3, imess, ifremd, nrjprt

      if (imess == 1) then
         c3 = c3_mess
         c6 = c6_ar_mess
      else
         c3 = c3_theo
         c6 = c6_theo
      end if

      if (ifremd > 12) then
         write (iunout,*) 
     .     ' too many fremddruckverbreiterungen specified'
         write (iunout,*) ' calculation abandonned '
      end if

!pb      if (nrjprt == 0) nrjprt = 1

      ipc6 = 0
      c6a = 0._dp
      do i=1,ifremd
        read (iunin,'(i6,1x,a2,3x,i6)') ii,kenn,ik6
        call uppercase(kenn)
        if (kenn == 'QS') then
          ipc6(i) = ik6
          c6a(i) = c6_qs_mess
        else
          do j=1,npolari
            if (kenn == polari_elnam(j)) then
              ipc6(i) = ik6
              c6a(i) = c6*polari_fac(j)
              exit
             end if
          end do
          if (j>npolari) then
            write (iunout,*) ' wrong code for foreign gas pressure',
     .                  ' broadening specified '
            write (iunout,*) ' code ',kenn,' not found in',
     .                  ' polarisation database '
            call exit_own(1)
        end if
        end if
      end do

      c2=0._dp

      reac_name(ir)(1:len_trim(reac)) = reac(1:len_trim(reac))

      allocate (phline)

      phline%aik = aik 
!     wl is in nm
      phline%e0 = hpcl / wl *1.E7_DP
      phline%g1 = gj
      phline%g2 = gi
!     Ej is in [1/cm]
      phline%e1 = ej * clight*hplanck*erg_to_ev
      phline%c2 = c2
      phline%c3 = c3
      phline%c4 = c4

      select case (isw)
        case (1)    ! absorbtion
           phline%ircart = 4
        case (2)    ! emission
           phline%ircart = 4
        case (3)    ! stimulated emission
        case (4)    ! ph_cabs ot (case(7) bei sven)
           phline%ircart = 7
        case default
           phline%ircart = 0
      end select

      phline%c6 = c6a(1)
      phline%c6a(1:12) = c6a(1:12)
      phline%iplsc6(1:12) = ipc6(1:12)
      phline%ifremd = count(phline%iplsc6(1:12)>0)

      phline%ignd = iplsc3
      phline%iprofiletype = iprftype

      reacdat(ir)%lphr = .true.

      allocate (reacdat(ir)%phr)
      allocate (reacdat(ir)%phr%line)
      reacdat(ir)%phr%ifit = -1
      reacdat(ir)%phr%line => phline

      call get_reaction(ir)
      b21=ph_b21()
      b12=b21*gi/gj
      phline%b21 = b21
      phline%b12 = b12
      phline%nrjprt = nrjprt

      phline%reacname = reac_name(ir)

      modclf(ir) = 100
      iftflg(ir,2) = 110
      
      rdata = 0._dp
      rdata(1,1) = aik

      call set_reaction_data (ir,isw,iftflg(ir,2),rdata,iunout,.false.)

      return

 990  continue
      write (iunout,*) 'REACTION ',reac,' NOT FOUND IN FILE PHOTON'
      call exit_own(1)


      contains


      subroutine subcomma (str)
      implicit none
      character(len=*), intent(in out) :: str
      integer :: i

      do
        i=scan(str,',')
        if (i == 0) exit
        str(i:i)='.'
      end do
      return
      end subroutine subcomma


      subroutine delete_blanks (str)
      implicit none
      character(len=*), intent (in out) :: str
      character, allocatable :: compact(:)
      integer :: i, ic, l

      l=len(str)
      allocate (compact(l))
      compact=' '

      i=1
      ic=1
      do while (i<=l)
        if (str(i:i) /= ' ') then
          compact(ic) = str(i:i)
          ic=ic+1
          i=i+1
        else
          i=i+1
        end if
      end do

      do i=1,l
        str(i:i) = compact(i)
      end do

      return
      end subroutine delete_blanks


      subroutine read_polari (name, factor, n)
      implicit none
      character(*), intent(out) :: name(:)
      real(dp), intent(out) :: factor(:)
      integer, intent(out) :: n
      character(200) :: zeile
      integer :: nmax, ifile

      nmax = min(size(name),size(factor))

      DO IFILE=1, NDBNAMES
        IF (INDEX(DBHANDLE(IFILE),'POLARI') /= 0) EXIT
      END DO

      IF (IFILE > NDBNAMES) THEN
        WRITE (IUNOUT,*) ' NO DATABASENAME FOR POLARI DEFINED '
        WRITE (IUNOUT,*) ' CALCULATION ABANDONNED '
        CALL EXIT_OWN(1)
      END IF

      open (unit=28,file=DBFNAME(IFILE),access='sequential',
     .      form='formatted')

      n=0
      read (28,*)
      read (28,*)

      do
        read (28,'(A200)',end=990) zeile
        if (zeile(1:2) == '--') cycle
        call subcomma(zeile)

        n = n + 1
        if (n > nmax) then
           write (iunout,*) ' too many lines in polarization file '
           write (iunout,*) ' calculation stopped '
           call exit_own(1)
        end if

        ianf = 3
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read element name
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        call delete_blanks(zeile(ianf:iend))
        iblnk = scan(zeile(ianf:iend),'|') - 1
        if (iblnk < 0) iblnk = iend-ianf+1
        if (iblnk == 0) then
           write (iunout,*) ' ERROR IN DATABASE POLARI'
           write (iunout,*) ' NO ELEMENTNAME FOUND '
           stop
        end if

        name(n) = repeat(' ',2)
        name(n)(1:iblnk) = zeile(ianf:ianf+iblnk-1)
        call uppercase(name(n))

!  skip polarization
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  skip estimated error
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read factor
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) factor(n)

      end do

 990  continue

      return
      end subroutine read_polari

      end subroutine read_photdbk

C ===== SOURCE: refdat.f
C
C
      SUBROUTINE REFDAT(TMM,TCC,WMM,WCC)
C
C  THIS SUBROUTINE READS REFLECTION DATA PRODUCED BY MONTE CARLO CODES
C    IFILE=1  H ON FE
C    IFILE=2  D ON FE
C    IFILE=3  H ON C
C    IFILE=4  D ON C
C    IFILE=5  HE ON FE
C    IFILE=6  HE ON C
C    IFILE=7  T ON FE
C    IFILE=8  T ON C
C    IFILE=9  D ON W
C    IFILE=10 HE ON W
C    IFILE=11 H ON W
C    IFILE=12 T ON W
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE CREFMOD
      USE CREF
      USE CSPEI
      USE CINIT

      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: TMM(*), TCC(*), WMM(*), WCC(*)
      REAL(DP) :: TML(12), TCL(12), WML(12), WCL(12),
     .          FELD(1092)
      REAL(DP) :: PID180
      INTEGER :: I1, I2, I3, I4, I5, NRECL, IUN, IFILE, I, J
C
      NFLR=NHD6
      IF (NHD6.GT.12) THEN
        WRITE (iunout,*) 'STORAGE ERROR. NHD6 MUST BE LESS OR EQUAL 12 '
        WRITE (iunout,*) 'NHD6= ',NHD6
        CALL EXIT_OWN(1)
      ENDIF
C
      NRECL=1092
      INE=12
      INEM=INE-1
      TML(1)=1.
      TCL(1)=1.
      WML(1)=56.
      WCL(1)=26.
      TML(2)=2.
      TCL(2)=1.
      WML(2)=56.
      WCL(2)=26.
      TML(3)=1.
      TCL(3)=1.
      WML(3)=12.
      WCL(3)=6.
      TML(4)=2.
      TCL(4)=1.
      WML(4)=12.
      WCL(4)=6.
      TML(5)=4.
      TCL(5)=2.
      WML(5)=56.
      WCL(5)=26.
      TML(6)=4.
      TCL(6)=2.
      WML(6)=12.
      WCL(6)=6.
      TML(7)=3.
      TCL(7)=1.
      WML(7)=56.
      WCL(7)=26.
      TML(8)=3.
      TCL(8)=1.
      WML(8)=12.
      WCL(8)=6.
      TML(9)=2.
      TCL(9)=1.
      WML(9)=184.
      WCL(9)=74.
      TML(10)=4.
      TCL(10)=2.
      WML(10)=184.
      WCL(10)=74.
      TML(11)=1.
      TCL(11)=1.
      WML(11)=184.
      WCL(11)=74.
      TML(12)=3.
      TCL(12)=1.
      WML(12)=184.
      WCL(12)=74.
      DO 10 J=1,NHD6
        TMM(J)=TML(J)
        TCC(J)=TCL(J)
        WMM(J)=WML(J)
        WCC(J)=WCL(J)
10    CONTINUE
      ENAR(1)=1.
      ENAR(2)=2.
      ENAR(3)=5.
      ENAR(4)=10.
      ENAR(5)=20.
      ENAR(6)=50.
      ENAR(7)=100.
      ENAR(8)=200.
      ENAR(9)=500.
      ENAR(10)=1000.
      ENAR(11)=2000.
      ENAR(12)=5000.
      DO 11 I=1,INEM
11      DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      INW=7
      INWM=INW-1
      PID180=ATAN(1.)/45.
      WIAR(1)=COS(0.D0)
      WIAR(2)=COS(30.*PID180)
      WIAR(3)=COS(45.*PID180)
      WIAR(4)=COS(60.*PID180)
      WIAR(5)=COS(70.*PID180)
      WIAR(6)=COS(80.*PID180)
      WIAR(7)=COS(85.*PID180)
      DO 12 I=1,INWM
12      DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
C
      INR=5
      INRM=INR-1
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
      DO 15 I=1,INRM
15      DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
C
      IF (INE*INW*NFLR.GT.NH0 .OR.
     .    INE*INW*INR*NFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*NFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*NFLR.GT.NH3) THEN
        WRITE (iunout,*) 
     .    'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EXIT_OWN(1)
      ENDIF
C
      DO IFILE=1, NDBNAMES
        IF (INDEX(DBHANDLE(IFILE),'TRIM') /= 0) EXIT
      END DO

      IF (IFILE > NDBNAMES) THEN
        WRITE (IUNOUT,*) ' NO DATABASENAME FOR TRIM DEFINED '
        WRITE (IUNOUT,*) ' CALCULATION ABANDONNED '
        CALL EXIT_OWN(1)
      END IF

      IUN=21
      OPEN (UNIT=IUN,FILE=DBFNAME(IFILE))
      REWIND IUN
C
660   FORMAT (1X,1P,10E12.4)
661   FORMAT (4E20.12)
      DO 1 IFILE=1,NFLR
        DO 2 I1=1,INE
          I=0
          READ (IUN,661) (FELD(J),J=1,NRECL)
          DO 3 I2=1,INW
            I=I+1
            HFTR0(I1,I2,IFILE)=FELD(I)
3         CONTINUE
          DO 4 I2=1,INW
            DO 4 I3=1,INR
              I=I+1
              HFTR1(I1,I2,I3,IFILE)=FELD(I)
4         CONTINUE
          DO 5 I2=1,INW
            DO 5 I3=1,INR
              DO 5 I4=1,INR
                I=I+1
                HFTR2(I1,I2,I3,I4,IFILE)=FELD(I)
5         CONTINUE
          DO 6 I2=1,INW
            DO 6 I3=1,INR
              DO 6 I4=1,INR
                DO 6 I5=1,INR
                  I=I+1
                  HFTR3(I1,I2,I3,I4,I5,IFILE)=FELD(I)
6         CONTINUE
2       CONTINUE
1     CONTINUE
C
      RETURN
      END
C ===== SOURCE: slreac.f
!  03.08.06:  data structure for reaction data redefined
C
C
      SUBROUTINE SLREAC (IR,FILNAM,H123,REAC,CRC,
     .                   RCMIN, RCMAX, FP, JFEXMN, JFEXMX, ELNAME, IZ1)
c
C  input
C    FILNAM: read a&m data from file filnam, e.g. AMJUEL, HYDHEL, METHAN, CONST
c    IR    : store data on eirene array CREAC(...,...,IR)
c    H123  : identifyer for data type in filnam, e.g. H.1, H.2, H.3, ...
c    REAC  : number of reaction in filnam, e.g. 2.2.5
c    CRC   : type of process, e.g. EI, CX, OT, etc
C  internal
C    ISW   <-- H123
C    IO    derived from ISW, initial value of 2nd index in CREAC
C  output
c    ISWR  : eirene flag for type of process  (1,2,...7)
c    CREAC : eirene storage array for a&m data CREAC(9,-1:9,IR)
c    MODCLF: see below
c    DELPOT: ionisation potential (for H.10 data),
c            currently handeled in input.f. not nice!
c
C    IFTFLG=IFTFLG(IR,IFLG)
C    IFLG  derived from ISW
C          0 for potential, 1 for cross section, 2 for rate-coeff,
C          3 for mom-weighted rate coeff. 4 for energy weighted rate coeff.
c    IFTFLG: eirene flag for type of fitting expression ("fit-flag=...")
C            DEFAULTS: =2 IFLG=0
C                         FOR POTENTIAL (GEN. MORSE)
C                         IFLG > 0:
C                         NOT IN USE
C                      =0 FOR ALL OTHERS (POLYNOM, DOUBLE POLYNOM)
C                      =3, IFLG=1:
C                          ionisation/exciation FORMULA (METHANE,...)
C                          IFLG > 1, IFLG=0:
C                          NOT IN USE
C                      =L10 (L=0,1): ONLY ONE CONSTANT RATE OR RATE-COEFF.
C                      =LMN L=0: rate coefficient. 
c                                multiply with density
C                      =LMN L=1: rate, not rate coefficient. 
c                                no multiplication of density
c
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM" INTO ARRAY "CREAC"
C    DEFINE PARAMETER MODCLF(IR) (5 DIGITS NMLKJ)
C    FIRST DEZIMAL  J           =1  POTENTIAL AVAILABLE
C                                   (ON CREAC(..,-1,IR))
C                   J           =0  ELSE
C    SECOND DEZIMAL K           =1  CROSS SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    THIRD  DEZIMAL L           =1  <SIGMA V> FOR ONE
C                                   PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  <SIGMA V> FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  <SIGMA V> FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   L           =0  ELSE
C    FOURTH DEZIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FIFTH  DEZIMAL N           =1  DELTA E FOR ONE PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  DELTA E FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  DELTA E FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   N           =0  ELSE
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT
      USE COMXS
      USE CINIT
      USE PHOTON

      IMPLICIT NONE

      INTEGER,      INTENT(IN) :: IR, IZ1
      INTEGER,      INTENT(IN OUT) :: JFEXMN, JFEXMX
      CHARACTER(8), INTENT(IN) :: FILNAM
      CHARACTER(4), INTENT(IN) :: H123
      CHARACTER(LEN=*), INTENT(IN) :: REAC, ELNAME
      CHARACTER(3), INTENT(IN) :: CRC
      REAL(DP), INTENT(IN OUT) :: RCMIN, RCMAX, FP(6) 
      CHARACTER(11) :: REACSTR
      REAL(DP) :: CONST
      REAL(DP) :: CREACD(9,9)
      INTEGER :: I, IND, J, K, IH, I0P1, I0, IC, IREAC, ISW, INDFF,
     .           IFLG, INC, IANF, IFILE, IL
      CHARACTER(80) :: ZEILE
      CHARACTER(6) :: AMJUEL, HYDHEL, H2VIBR, SPECTR
      CHARACTER(7) :: METHANE
      CHARACTER(2) :: CHR
      CHARACTER(3) :: CHRL, CHRR
      CHARACTER(200) :: DSN, DIR
      CHARACTER(1) :: CUT
      LOGICAL :: LCONST,LGEMIN,LGEMAX
C
      LGEMIN=.FALSE.
      LGEMAX=.FALSE.
      ISWR(IR)=0
      CONST=0.
      CHR='l0'
      I0=0
      CREACD = 0._DP
C
      AMJUEL='AMJUEL'
      HYDHEL='HYDHEL'
      METHANE='METHANE'
      H2VIBR='H2VIBR'
      SPECTR='SPECTR  '
C
      IF (INDEX(CRC,'EI').NE.0.OR.
     .    INDEX(CRC,'DS').NE.0) ISWR(IR)=1
      IF (INDEX(CRC,'CX').NE.0) ISWR(IR)=3
      IF (INDEX(CRC,'II').NE.0.OR.
     .    INDEX(CRC,'PI').NE.0) ISWR(IR)=4
      IF (INDEX(CRC,'EL').NE.0) ISWR(IR)=5
      IF (INDEX(CRC,'RC').NE.0) ISWR(IR)=6
      IF (INDEX(CRC,'OT').NE.0) ISWR(IR)=7
C
      IF (INDEX(FILNAM,'CONST').NE.0) THEN
        LCONST=.TRUE.
      ELSE
        DO IFILE=1,NDBNAMES
          IF (INDEX(FILNAM,DBHANDLE(IFILE)).NE.0) EXIT
        END DO
        IF (IFILE <= NDBNAMES) THEN
          LCONST=.FALSE.
          IF (INDEX(FILNAM,'ADAS') == 0) THEN
            OPEN (UNIT=29,FILE=DBFNAME(IFILE))
          ELSE
! FIND NAME OF ADAS-FILE TO BE READ 
            DIR = ' '
            IL = 0
            IF (VERIFY(DBFNAME(IFILE),' ') .NE. 0) THEN 
              IC = SCAN(DBFNAME(IFILE),'/\\')
              CUT = DBFNAME(IFILE)(IC:IC)
              DIR=TRIM(DBFNAME(IFILE)) // CUT // 
     .            ADJUSTL(TRIM(REAC)) // CUT 
              IL = INDEX(DIR,CUT,.TRUE.)
            END IF
            IF (IL == 0) THEN
              DSN = ADJUSTL(TRIM(REAC)) // '_' // TRIM(ELNAME) // '.dat'
            ELSE
              DSN = DIR(1:IL) // ADJUSTL(TRIM(REAC)) // '_' // 
     .              TRIM(ELNAME) // '.dat'
            END IF
            OPEN (UNIT=29,FILE=DSN)
          END IF
        ELSE
          WRITE (iunout,*) 
     .      ' NO SPECIFICATION FOR FILENAME IN REACTION CARD'
          WRITE (iunout,*) ' CHOOSE EITHER '
          WRITE (iunout,*) ' AMJUEL, METHAN, HYDHEL, H2VIBR, SPECTR '
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) 
     .      ' CONST FOR ENTERING REACTION COEFFICIENTS VIA '
          WRITE (iunout,*) ' EIRENE INPUT-FILE '
          CALL EXIT_OWN(1)
        END IF
      ENDIF
C
      IF (H123(4:4).EQ.' ') THEN
        READ (H123(3:3),'(I1)') ISW
      ELSE
        READ (H123(3:4),'(I2)') ISW
      ENDIF

C
      IF (INDEX(FILNAM,'PHOTON').NE.0) THEN
        CALL READ_PHOTDBK (IR,REAC,ISW)
        RETURN
      END IF

      REACSTR=REPEAT(' ',11)
      IANF=VERIFY(REAC,' ')
      IF (IANF > 0) THEN
        IREAC=INDEX(REAC(IANF:),' ')-1
        IF (IREAC.LT.0) IREAC=LEN(REAC(IANF:))
        REACSTR(2:IREAC+1)=REAC(IANF:IREAC+IANF-1)
C  ADD ONE MORE BLANK, IF POSSIBLE
        IREAC=IREAC+2
      ELSE
        IF (.NOT.LCONST) THEN
          WRITE (iunout,*) ' NO REACTION SPECIFIED IN REACTION CARD ',IR
          CALL EXIT_OWN (1)
        END IF
      END IF

      REAC_NAME(IR) = REACSTR(2:)  

C  H.0
      IF (ISW.EQ.0) THEN
        CHR='p0'
        CHRL='pl0'
        CHRR='pr0'
        I0=-1
        MODCLF(IR)=MODCLF(IR)+1
        IFLG=0
C  DEFAULT POTENTIAL: GENERALISED MORSE
        IFTFLG(IR,IFLG)=2
C  H.1
      ELSEIF (ISW.EQ.1) THEN
        CHR='a0'
        CHRL='al0'
        CHRR='ar0'
        I0=0
        MODCLF(IR)=MODCLF(IR)+10
        IFLG=1
C  DEFAULT CROSS SECTION: 8TH ORDER POLYNOM OF LN(SIGMA)
        IFTFLG(IR,IFLG)=0
C  H.2
      ELSEIF (ISW.EQ.2) THEN
        CHR='b0'
        CHRL='bl0'
        CHRR='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
C  DEFAULT RATE COEFFICIENT: 8TH ORDER POLYNOM OF LN(<SIGMA V>) FOR E0=0.
        IFTFLG(IR,IFLG)=0
C  H.3
      ELSEIF (ISW.EQ.3) THEN
        MODCLF(IR)=MODCLF(IR)+200
        I0=1
        IFLG=2
        IFTFLG(IR,IFLG)=0
C  H.4
      ELSEIF (ISW.EQ.4) THEN
        MODCLF(IR)=MODCLF(IR)+300
        I0=1
        IFLG=2
        IFTFLG(IR,IFLG)=0
C  H.5
      ELSEIF (ISW.EQ.5) THEN
        CHR='e0'
        CHRL='el0'
        CHRR='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.6
      ELSEIF (ISW.EQ.6) THEN
        MODCLF(IR)=MODCLF(IR)+2000
        I0=1
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.7
      ELSEIF (ISW.EQ.7) THEN
        MODCLF(IR)=MODCLF(IR)+3000
        I0=1
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.8
      ELSEIF (ISW.EQ.8) THEN
        CHR='h0'
        CHRL='hl0'
        CHRR='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10000
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.9
      ELSEIF (ISW.EQ.9) THEN
        MODCLF(IR)=MODCLF(IR)+20000
        I0=1
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.10
      ELSEIF (ISW.EQ.10) THEN
        MODCLF(IR)=MODCLF(IR)+30000
        I0=1
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.11
      ELSEIF (ISW.EQ.11) THEN
        CHR='k0'
        CHRL='kl0'
        CHRR='kr0'
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
C  H.12
      ELSEIF (ISW.EQ.12) THEN
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
      ENDIF

      IF (INDEX(FILNAM,'ADAS').NE.0) THEN
        CALL READ_ADAS (IR,REAC,ISW,IZ1)
        RETURN
      END IF
C
      IF (LCONST) THEN
        IND=INDEX(REACSTR,'FT')
        IF (IND /= 0) THEN
          READ (REACSTR(IND+2:),*) IFTFLG(IR,IFLG)
        END IF

        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
C
C  READ ONLY ONE FIT COEFFICIENT FROM INPUT FILE
          READ (IUNIN,6664) CREACD(1,1)
        ELSE
C
C  READ 9 FIT COEFFICIENTS FROM INPUT FILE
          READ (IUNIN,6664) (CREACD(IC,1),IC=1,9)

        END IF
        CALL SET_REACTION_DATA(IR,ISW,IFTFLG(IR,IFLG),CREACD,
     .                         IUNOUT,.FALSE.)
        RETURN
C
C  READ FROM DATA FILE
C
      ELSEIF (.NOT.LCONST) THEN
100     READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0) GOTO 100

1       READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,H123).EQ.0) GOTO 1
C
2       READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'H.').NE.0) GOTO 990
        IF (INDEX(ZEILE,'Reaction ').EQ.0.or.
     .      INDEX(ZEILE,REACSTR(1:ireac)).EQ.0) GOTO 2
      ENDIF
C
C  SINGLE PARAM. FIT, ISW=0,1,2,5,8,11
      IF (ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.
     .    ISW.EQ.11) THEN
        IF (.NOT.LCONST) THEN
3         READ (29,'(A80)',END=990) ZEILE
          INDFF=INDEX(ZEILE,'fit-flag')
          IF (INDEX(ZEILE,CHR)+INDFF.EQ.0) GOTO 3
          IF (INDFF > 0) THEN
            READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
            GOTO 3
          ENDIF
          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
            IND=INDEX(ZEILE,CHR(1:1))
            READ (ZEILE((IND+2):80),'(E20.12)') CREACD(1,1)
          ELSE
            DO 9 J=0,2
              IND=0
              DO 4 I=1,3
                IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
                READ (ZEILE((IND+2):80),'(E20.12)') CREACD(J*3+I,1)
4             CONTINUE
              READ (29,'(A80)',END=990) ZEILE
9           CONTINUE
          END IF
C
C  READ ASYMPTOTICS, IF AVAILABLE
C  I0P1=1 FOR CROSS SECTION
C  I0P1=2 FOR (WEIGHTED) RATE
          I0P1=I0+1
          IF (ISW.EQ.0) GOTO 12 ! NO ASYMPTOTICS FOR POTENTIALS
          IF (INDEX(ZEILE,CHRL).NE.0.AND.JFEXMN.EQ.0) THEN
            IND=0
            DO 5 I=1,3
              INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
              IF (INC.GT.0) THEN
                IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
                READ (ZEILE((IND+3):80),'(E20.12)') FP(I)
              ENDIF
5           CONTINUE
            LGEMIN=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
          IF (INDEX(ZEILE,CHRR).NE.0.AND.JFEXMX.EQ.0) THEN
            IND=0
            DO 7 I=4,6
              INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
              IF (INC.GT.0) THEN
                IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
                READ (ZEILE((IND+3):80),'(E20.12)') FP(I)
              ENDIF
7           CONTINUE
            LGEMAX=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
c
          if (lgemin.and.jfexmn.eq.0) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmin
            rcmin=log(rcmin)
            jfexmn=5
            READ (29,'(A80)',END=990) ZEILE
          endif
          if (lgemax.and.jfexmx.eq.0) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmax
            rcmax=log(rcmax)
            jfexmx=5
            READ (29,'(A80)',END=990) ZEILE
          endif
C
C  ANY OTHER ASYMPTOTICS INFO ON FILE?  SEARCH FOR Tmin, or Emin
          IF ((INDEX(ZEILE,'Tmin').NE.0.and.I0P1==2).or.
     .        (INDEX(ZEILE,'Emin').NE.0.and.I0P1==1)) then
            IND=INDEX(ZEILE,'n')
            READ (ZEILE((IND+2):80),'(E9.2)') rcmin
            rcmin=log(rcmin)
C  extrapolation from subr. CROSS
            if (I0P1.eq.1.and.iswr(ir).eq.1) jfexmn=1
            if (I0P1.eq.1.and.iswr(ir).eq.3) jfexmn=-1
            if (I0P1.eq.1.and.iswr(ir).eq.5) jfexmn=-1
C  extrapolation from subr. CDEF
C   ??      if (I0PT.eq.2) jfexmn=-1
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
12        CONTINUE
C       ELSEIF (LCONST) THEN
C  NOTHING TO BE DONE
        ENDIF
C
C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12
      ELSEIF (ISW.EQ.3.OR.ISW.EQ.4.OR.ISW.EQ.6.OR.ISW.EQ.7.OR.
     .        ISW.EQ.9.OR.ISW.EQ.10.OR.ISW.EQ.12) THEN
        DO 11 J=0,2
16        READ (29,'(A80)',END=990) ZEILE
          INDFF=INDEX(ZEILE,'fit-flag')
          IF (INDEX(ZEILE,'Index')+INDFF.EQ.0) GOTO 16
          IF (INDFF > 0) THEN
            READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
            GOTO 16
          ENDIF
          READ (29,'(1X)')
          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
            READ (29,*) IH,CREACD(1,1)
            EXIT
          ELSE
            DO 17 I=1,9
              READ (29,*) IH,(CREACD(I,K),K=J*3+1,J*3+3)
17          CONTINUE
          END IF
11      CONTINUE
C   NO ASYMPTOTICS AVAILABLE YET
C
      ENDIF

      CALL SET_REACTION_DATA(IR,ISW,IFTFLG(IR,IFLG),CREACD,IUNOUT,
     .                       .TRUE.,RCMIN,RCMAX,FP,JFEXMN,JFEXMX)
C
      CLOSE (UNIT=29)
C
      RETURN
C
990   WRITE (iunout,*) ' NO DATA FOUND FOR REACTION ',H123,' ',REAC,
     .            ' IN DATA SET ',FILNAM
      WRITE (iunout,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29)
      CALL EXIT_OWN(1)
991   WRITE (iunout,*) ' INVALID CONSTANT IN SLREAC. CONST= ',CONST
      WRITE (iunout,*) ' CHECK "REACTION CARDS" FOR REACTION NO. ',IR
      CLOSE (UNIT=29)
      CALL EXIT_OWN(1)
6664  FORMAT (6E12.4)
      END
C ===== SOURCE: wrgeom.f
C
      SUBROUTINE WRGEOM(TRCFLE)
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE COMPRT, ONLY: IUNOUT
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CLGIN
      USE CTRIG
      IMPLICIT NONE
      LOGICAL TRCFLE
C
      OPEN (UNIT=12,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 12
      WRITE (12) RCGM1,RCGM2,NPOINT,NSTGRD,NGHPLS,NGHPOL,NCLTAL
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 12: RCGM,ICGM '
      WRITE (12) RCGRID,ICGRID
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 12: RCGRID,ICGRID'
      WRITE (12) RCPLYG,RCPLY2,ICPLYG
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 12: RCPLYG,ICPLYG'
      IF (NRKNOT > 0) THEN
        WRITE (12) XTRIAN,YTRIAN,VTRIX,VTRIY,PTRIX,PTRIY,
     .             NECKE,NCHBAR,NSEITE,INMTI,NRKNOT,NTRII
        IF (TRCFLE) WRITE (iunout,*) 'WRITE 12: RCTRIG,ICTRIG'
      END IF
      WRITE (12)
     R RLWMN,RLWMX,EWALL,EWBIN,TRANSP,FSHEAT,
     R ZNML,ZNCL,
     R RECYCF,RECYCT,RECPRM,EXPPL,EXPEL,EXPIL,
     R RECYCS,RECYCC,SPTPRM,
     I ILSWCH,ILEQUI,ILTOR,ILSIDE,ILIIN,ILREF,
     I ILSPT,ILCOL,ILFIT,ILCELL,ILBOX,ILPLG,ISPUT,
     I NLIMII,NLIMIE,ISWICH,ILBLCK,ILACLL,JUMLIM,
     I NSTSI,INUMP,IRPTA,IRPTE,ISRF,ISRT,ISRS,ISRC,
     I INMP1I,INMP2I,INMP3I,
     I IGFIL,IGJUM0,IGJUM1,IGJUM2,IGJUM3
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 12: RCLGN,ICLGN,LCLGN'
      WRITE (12) RADGEO,IADGEO,NLIMI
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 12: RCADG,ICADG'
      CLOSE (UNIT=12)
      RETURN
C
      ENTRY RGEOM(TRCFLE)
      OPEN (UNIT=12,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 12
      READ (12) RCGM1,RCGM2,NPOINT,NSTGRD,NGHPLS,NGHPOL,NCLTAL
      IF (TRCFLE) WRITE (iunout,*) 'READ 12: RCGM,ICGM '
      READ (12) RCGRID,ICGRID
      IF (TRCFLE) WRITE (iunout,*) 'READ 12: RCGRID,ICGRID'
      READ (12) RCPLYG,RCPLY2,ICPLYG
      IF (TRCFLE) WRITE (iunout,*) 'READ 12: RCPLYG,ICPLYG'
      IF (NRKNOT > 0) THEN
        READ (12) XTRIAN,YTRIAN,VTRIX,VTRIY,PTRIX,PTRIY,
     .            NECKE,NCHBAR,NSEITE,INMTI,NRKNOT,NTRII
        IF (TRCFLE) WRITE (iunout,*) 'READ 12: RCTRIG,ICTRIG'
      END IF
      READ (12)
     R RLWMN,RLWMX,EWALL,EWBIN,TRANSP,FSHEAT,
     R ZNML,ZNCL,
     R RECYCF,RECYCT,RECPRM,EXPPL,EXPEL,EXPIL,
     R RECYCS,RECYCC,SPTPRM,
     I ILSWCH,ILEQUI,ILTOR,ILSIDE,ILIIN,ILREF,
     I ILSPT,ILCOL,ILFIT,ILCELL,ILBOX,ILPLG,ISPUT,
     I NLIMII,NLIMIE,ISWICH,ILBLCK,ILACLL,JUMLIM,
     I NSTSI,INUMP,IRPTA,IRPTE,ISRF,ISRT,ISRS,ISRC,
     I INMP1I,INMP2I,INMP3I,
     I IGFIL,IGJUM0,IGJUM1,IGJUM2,IGJUM3
      IF (TRCFLE) WRITE (iunout,*) 'READ 12: RCLGN,ICLGN,LCLGN'
      READ (12) RADGEO,IADGEO,NLIMI
      IF (TRCFLE) WRITE (iunout,*) 'READ 12: RCADG,ICADG'
      CLOSE (UNIT=12)
      RETURN
      END
C ===== SOURCE: wrmesh.f
c  nov 03:  use relative distances to find neighbor segment,
c           otherwise sometimes problems with non-closing polygons encountered.
      SUBROUTINE WRMESH
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE COMPRT, ONLY: IUNOUT
      USE CPLOT
      USE CPOLYG
      USE CGEOM
      USE CLGIN
      USE CGRPTL
      USE CTRIG
      USE CGRID
      IMPLICIT NONE


      INTEGER, PARAMETER :: MAXPOIN=1000
      REAL(DP) :: partcont(maxpoin,2,2), maxlen
      REAL(DP) :: XPE, YPE, HELP, XT, YT, PHI1, X1, X2, Y1, Y2, PHI2
      REAL(DP) :: DISTQI,DISTQJ1,DISTQJ2
      INTEGER  :: MAXCONT, ICONT, IPOIN, IWST, IWEN, IWL, IWP,
     .            IWAN, IMN, I, NCONT, J, IUHR, ISTORE, IP, IH, IFOUND,
     .            ICO, IPO, IN, IS, IS1, ITRI
      INTEGER  :: IDIAG(MAXPOIN),irip(maxpoin,2)
      REAL(SP) :: xmin,xmax,ymin,ymax,deltax,deltay,delta,xcm,ycm
      REAL(SP) :: XP,YP

C INITIALISIERUNG DER PLOTDATEN
      xmin = CH2X0-CH2MX
      ymin = CH2Y0-CH2MY
      xmax = CH2X0+CH2MX
      ymax = CH2Y0+CH2MY
      deltax = abs(xmax-xmin)
      deltay = abs(ymax-ymin)
      delta = max(deltax,deltay)
      xcm = 24. * deltax/delta
      ycm = 24. * deltay/delta

C ANZAHL DER KONTOUREN BESTIMMEN
C ILPLG WIRD IM INPUT-BLOCK 3 EINGELESEN
      CALL LEER(2)
      WRITE (iunout,*) 'SUBROUTINE WRMESH CALLED '
      CALL LEER(1)
      NCONT = 0
      DO I=1,NLIMI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO
      DO I=NLIM+1,NLIM+NSTSI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO

      WRITE (iunout,*) 'NUMBER OF CONTOURS FOR FEM MESH: ',NCONT
      CALL LEER(1)
      if (ncont == 0) return

      ALLOCATE (NCONPOINT(NCONT))
      ALLOCATE (XCONTOUR(MAXPOIN,NCONT))
      ALLOCATE (YCONTOUR(MAXPOIN,NCONT))
      NCONPOINT = 0
      ICO = 0

      call grnxtf
      call grsclc(3.,3.,3.+real(xcm,kind(1.e0)),3.+real(ycm,kind(1.e0)))
      call grsclv(real(xmin,kind(1.e0)),real(ymin,kind(1.e0)),
     .            real(xmax,kind(1.e0)),real(ymax,kind(1.e0)))

      DO ICONT = 1,NCONT
        IPOIN = 0
        MAXLEN = 0.
        irip=0
C AKTUELLE KONTOUR BESTIMMEN, STUECKE MIT ILPLG=ICONT GEHOEREN ZUR
C AKTUELLEN CONTOUR, ANFANGS UND ENDPUNKT DIESES STUECKES WERDEN AUF
C PARTCONT GESPEICHERT
        DO I=1,NLIMI
          IF (ABS(ILPLG(I)) .EQ. ICONT) THEN
            IUHR=ILPLG(I)
C 0 < RLB(I) < 2
C 2-PUNKT OPTION WIRD IM TIMEA0 AUF RLB=1 ZURUECKGEFUEHRT
            IF ((RLB(I) .GT. 0.) .AND. (RLB(I) .LT. 2.) .AND.
     >          (P3(1,I) .EQ. 1.D55 .OR. P3(2,I) .EQ. 1.D55
     >          .OR. P3(3,I) .EQ. 1.D55)) THEN
              IPOIN = IPOIN + 1
              IF (A3LM(I) .EQ. 0.) THEN
C               X,Y-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(2,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(2,I)
                idiag(ipoin)=i
              ELSEIF (A2LM(I) .EQ. 0.) THEN
C               X,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ELSEIF (A1LM(I) .EQ. 0.) THEN
C               Y,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(2,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(2,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ENDIF
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
            ELSE
C  ERROR
              WRITE(iunout,*) 'FALSCHE ANGABE FUER RLB, RLB = ',RLB(I),
     >                    ILPLG(I),I
            ENDIF
          ENDIF
        ENDDO

        IF (LEVGEO == 3) THEN
        DO I=1,NSTSI
          IF (ABS(ILPLG(NLIM+I)) .EQ. ICONT) THEN
            IUHR=ILPLG(NLIM+I)
            IF (INUMP(I,2) .NE. 0) THEN
C             POLOIDAL
              DO J=IRPTA(I,1),IRPTE(I,1)-1
                IF ((XPOL(J,INUMP(I,2)) .NE. XPOL(J+1,INUMP(I,2))) .OR.
     >              (YPOL(J,INUMP(I,2)) .NE. YPOL(J+1,INUMP(I,2)))) THEN
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,1,2) = YPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,2,1) = XPOL(J+1,INUMP(I,2))
                  PARTCONT(IPOIN,2,2) = YPOL(J+1,INUMP(I,2))
                idiag(ipoin)=-i
                irip(ipoin,1)=j
                irip(ipoin,2)=INUMP(I,2)
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSEIF (INUMP(I,1) .NE. 0) THEN
C             RADIAL
              DO J=IRPTA(I,2),IRPTE(I,2)-1
                IF ((XPOL(INUMP(I,1),J) .NE. XPOL(INUMP(I,1),J+1)) .OR.
     >              (YPOL(INUMP(I,1),J) .NE. YPOL(INUMP(I,1),J+1))) THEN
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,1,2) = YPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,2,1) = XPOL(INUMP(I,1),J+1)
                  PARTCONT(IPOIN,2,2) = YPOL(INUMP(I,1),J+1)
                  idiag(ipoin)=-i
                  irip(ipoin,1)=INUMP(I,1)
                  irip(ipoin,2)=j
                  maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSE
C  ERROR
              WRITE(iunout,*) 'CASE NOT FORESEEN: INUMP: ',
     >                     (INUMP(I,J),J=1,3)
            ENDIF
          ENDIF
        ENDDO
        ELSEIF (LEVGEO == 4) THEN
          DO ITRI = 1, NTRII
            DO IS = 1, 3
              IN=INMTI(IS,ITRI)
              IF (IN /= 0) THEN
                IF (ABS(ILPLG(IN)) == ICONT) THEN
                  IUHR=ILPLG(IN)
                  IS1 = IS+1
                  IF (IS1 > 3) IS1=1
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,1,2) = YTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,2,1) = XTRIAN(NECKE(IS1,ITRI))
                  PARTCONT(IPOIN,2,2) = YTRIAN(NECKE(IS1,ITRI))
                  idiag(ipoin)=IN
                  irip(ipoin,1)=itri
                  irip(ipoin,2)=is
                  maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)

                END IF
              END IF
            END DO
          END DO
        END IF

        IF (IPOIN.LE.0) THEN
          WRITE(iunout,*) 'CONTOUR ',ICONT,' NOT FOUND'
          GOTO 1000
        ENDIF

        call grnwpn(icont)
C STUECKE DER AKTUELLEN KONTOUR WERDEN SORTIERT
        DO I=1,IPOIN-1
          XPE = PARTCONT(I,2,1)
          YPE = PARTCONT(I,2,2)
          DISTQI=(PARTCONT(I,2,1)-PARTCONT(I,1,1))**2+
     .           (PARTCONT(I,2,2)-PARTCONT(I,1,2))**2
          IFOUND=0
          DO J=I+1,IPOIN
            DISTQJ1=(XPE-PARTCONT(J,1,1))**2+
     .              (YPE-PARTCONT(J,1,2))**2
            DISTQJ2=(XPE-PARTCONT(J,2,1))**2+
     .              (YPE-PARTCONT(J,2,2))**2
CDR         IF ((XPE .EQ. PARTCONT(J,1,1)) .AND.
CDR  >          (YPE .EQ. PARTCONT(J,1,2))) THEN
            IF (DISTQJ1/DISTQI.LE.1.D-10) THEN
              IFOUND=1
              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
            ELSEIF (DISTQJ2/DISTQI.LE.1.D-10) THEN
CDR         ELSEIF ((XPE .EQ. PARTCONT(J,2,1)) .AND.
CDR  >              (YPE .EQ. PARTCONT(J,2,2))) THEN
              IFOUND=1
              HELP = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
            ENDIF
          ENDDO
          IF (IFOUND.EQ.0) THEN
            WRITE (iunout,*) 'NO MATCHING POINT FOUND FOR CONTOUR ',
     >                        ICONT
            write(iunout,*) i,idiag(i),irip(i,1),irip(i,2),
     >                   partcont(i,1,1),partcont(i,1,2),
     >                   partcont(i,2,1),partcont(i,2,2)
            WRITE (iunout,*) 'USE NEXT POINT '
            IP=I+1
            write(iunout,*) iP,idiag(iP),irip(ip,1),irip(ip,2),
     >                    partcont(iP,1,1),partcont(iP,1,2),
     >                    partcont(iP,2,1),partcont(iP,2,2)
          ENDIF
        ENDDO
        IF ((PARTCONT(1,1,1) .NE. PARTCONT(IPOIN,2,1)) .OR.
     >      (PARTCONT(1,1,2) .NE. PARTCONT(IPOIN,2,2))) THEN
          WRITE(iunout,*) 'CONTOUR ',ICONT,' IS NOT CLOSED'
        ELSE
          WRITE(iunout,*) 'CLOSED CONTOUR ',ICONT
        ENDIF
        do i=1,ipoin
          write(iunout,*) i,idiag(i),irip(i,1),irip(i,2),
     >                 partcont(i,1,1),partcont(i,1,2),
     >                 partcont(i,2,1),partcont(i,2,2)
        enddo
        XP = PARTCONT(1,1,1)
        YP = PARTCONT(1,1,2)
        call grjmp(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        DO I=2,IPOIN
          XP = PARTCONT(I,1,1)
          YP = PARTCONT(I,1,2)
          call grdrw(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        ENDDO
        XP = PARTCONT(1,1,1)
        YP = PARTCONT(1,1,2)
        call grDRW(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))

C  BERECHNUNG VON DELTA ALS MITTLERE LAENGE DER TEILSTUECKE
C  DELTA IST MASS FUER DIE GROESSE DER DREIECKE
        IMN=0
        IF (ICONT .EQ. 1) THEN
          maxlen = maxlen / REAL(IPOIN,KIND(1.D0))
          WRITE(78,*) maxlen
          WRITE(78,*)
        endif

C BESTIMMUNG DES UHRZEIGERSINNS DER KONTOUR
C KONTOUR MUSS FUER DIE TRIANGULIERUNG FOLGENDERMASSEN AUSGEGEBEN
C WERDEN:
C  - IM UHRZEIGERSINN FUER INNERE BEGRENZUNGEN DES GEBIETES (POSITIV)
C  - GEGEN UHRZEIGERSINN FUER AEUSSERE BEGRENZUNGEN DES GEBIETES (NEGATIV)
        YMIN=PARTCONT(1,1,2)
        DO I=1,IPOIN
          IF (PARTCONT(I,2,2) .LT. YMIN) THEN
            YMIN = PARTCONT(I,2,2)
            IMN=I
          ENDIF
        ENDDO

        IF (IMN .EQ. 0) THEN
          XT = PARTCONT(1,1,1)
          YT = PARTCONT(1,1,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
          X1 = PARTCONT(IPOIN,1,1)
          Y1 = PARTCONT(IPOIN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
          X2 = PARTCONT(1,2,1)
          Y2 = PARTCONT(1,2,2)
        ELSE
C  SONDERFALL IMN=IPOIN ENTFAELLT, DA ERSTER PUNKT GLEICH LETZTER
C  PUNKT GILT
          XT = PARTCONT(IMN,2,1)
          YT = PARTCONT(IMN,2,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
          X1 = PARTCONT(IMN,1,1)
          Y1 = PARTCONT(IMN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
          X2 = PARTCONT(IMN+1,2,1)
          Y2 = PARTCONT(IMN+1,2,2)

        ENDIF

C  BESTIMME POLARWINKEL VON (X1,Y1) UND (X2,Y2) MIT (XT,YT) ALS URSPRUNG
        PHI1 = ATAN2 (Y1-YT,X1-XT)
        PHI2 = ATAN2 (Y2-YT,X2-XT)

        IF (PHI2 .GT. PHI1) THEN
C  ABSPEICHERUNG ERFOLGTE IM UHRZEIGERSINN
          ISTORE = 1
        ELSE
          ISTORE = -1
        ENDIF
C  IUHR=ILPLG > 0 ==> IM UHRZEIGERSINN AUSGEBEN
C  IUHR=ILPLG < 0 ==> ENTGEGEN DEM UHRZEIGERSINN AUSGEBEN
        IWAN=1
        IWEN=IPOIN
        IWST=1
        IWP=1
        IWL=2
        IF (ISTORE*IUHR .LT. 0.) THEN
          IWAN=IPOIN
          IWEN=1
          IWST=-1
          IWP=2
          IWL=1
        ENDIF

        WRITE(78,*) IPOIN+1
        IF (IUHR > 0) THEN
          ICO=ICO+1
          NCONPOINT(ICO)=IPOIN+1
          IPO=0
        END IF
        DO I=iwan,iwen,iwst
          WRITE(78,'(1P,2(2X,E21.14))')
     >          PARTCONT(I,IWP,1),PARTCONT(I,IWP,2)
          IF (IUHR > 0) THEN
            IPO=IPO+1
            XCONTOUR(IPO,ICO) = PARTCONT(I,IWP,1)
            YCONTOUR(IPO,ICO) = PARTCONT(I,IWP,2)
          END IF
        ENDDO
        WRITE(78,'(1P,2(2X,E21.14))') PARTCONT(IWEN,IWL,1),
     >                               PARTCONT(IWEN,IWL,2)
        IF (IUHR > 0) THEN
          IPO=IPO+1
          XCONTOUR(IPO,ICO) = PARTCONT(IWEN,IWL,1)
          YCONTOUR(IPO,ICO) = PARTCONT(IWEN,IWL,2)
        END IF
1000    CONTINUE
      ENDDO
      NCONTOUR=ICO
C  END OF NCONT LOOP
      call leer(1)
      write (iunout,*) 
     .  'input file fort.78 for FEM mesh generator written '
      call leer(2)
      call grnwpn(1)
      call grnxtf
      END
C ===== SOURCE: wrplam.f
c  sept. 05:  five more tallies added to step function, see also CSTEP.f
c  nov.  05:  add eltot and ve to step function data

C  write plasma (background) data, source distribution and atomic data
C  on unit 13.
C
c  at entry RPLAM:
C  read plasma (background) data, source distribution and atomic data
C  from unit 13.
C
C  trcfle:  confirm writing on printout on unit IUNOUT
C  IFLG  :

      SUBROUTINE WRPLAM(TRCFLE,IFLG)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CZT1
      USE COMSOU
      USE CSTEP
      USE COMXS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      LOGICAL TRCFLE
      REAL(DP), ALLOCATABLE :: RDUM(:)
      INTEGER, ALLOCATABLE :: IDUM(:)
      LOGICAL, ALLOCATABLE :: LDUM(:)
      INTEGER :: NRDUM, NIDUM, NLDUM
C
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      WRITE (13)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,BXPERP,BYPERP,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,DPHD,
     R           DION,DATM,DMOL,DPLS,DPHOT,
     R           TVAC,DVAC,VVAC,ALLOC,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,
C  LUSR, LOGICAL
     L           LGVAC,LGDFT
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module COMUSR.f '
      CALL WRITE_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMDTA,ICMDTA'
      CALL WRITE_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMAMF,ICMAMF'
      WRITE (13) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCZT1,RCZT2,ZT1,ZRG'
      WRITE (13) RCMSOU,SREC,EIO,EEL,
     .           ICMSOU,INGRDA,INGRDE,NSTRAI,
     .           LCMSOU,NLSYMP,NLSYMT
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMSOU,ICMSOU,LCMSOU'
      IF (ALLOCATED(FLSTEP))
     .  WRITE (13) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module CSTEP.f'
      CLOSE (UNIT=13)
      RETURN
C
      ENTRY RPLAM(TRCFLE,IFLG)
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      READ (13)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,BXPERP,BYPERP,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,DPHD,
     R           DION,DATM,DMOL,DPLS,DPHOT,
     R           TVAC,DVAC,VVAC,ALLOC,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,
C  LUSR, LOGICAL
     L           LGVAC,LGDFT
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: module COMUSR.f '
      CALL READ_CMDTA
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMDTA,ICMDTA'
      CALL READ_CMAMF
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMAMF,ICMAMF'
      READ (13) RCZT1,RCZT2,ZT1,ZRG
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCZT1,RCZT2,ZT1,ZRG'
      IF (IFLG == 0) THEN
        READ (13) RCMSOU,SREC,EIO,EEL,
     .            ICMSOU,INGRDA,INGRDE,NSTRAI,
     .            LCMSOU,NLSYMP,NLSYMT
      ELSE
        NRDUM = SIZE(RCMSOU) + SIZE(SREC) + SIZE(EIO) + SIZE(EEL)
        NIDUM = SIZE(ICMSOU) + SIZE(INGRDA) + SIZE(INGRDE) + 1
        NLDUM = SIZE(LCMSOU) + SIZE(NLSYMP) + SIZE(NLSYMT)
        ALLOCATE (RDUM(NRDUM))
        ALLOCATE (IDUM(NIDUM))
        ALLOCATE (lDUM(NLDUM))
        READ (13) RDUM, IDUM, LDUM
        DEALLOCATE (RDUM)
        DEALLOCATE (IDUM)
        DEALLOCATE (lDUM)
      END IF
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMSOU,ICMSOU,LCMSOU'
      IF (ALLOCATED(FLSTEP))
     .   READ (13) FLSTEP,ELSTEP,FLTOT,ELTOT,VF,VE,
     .             QUOT,ADD,QUOTI,ADDIV,
     .             TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .             FESTEP,FISTEP,SHSTEP,VPSTEP,MCSTEP,
     .             IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .             ISTUF,NSMAX,NSPSTI,NSPSTE
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: module CSTEP.f '
      CLOSE (UNIT=13)
      RETURN
      END
C ===== SOURCE: wrplam_xdr.f
C  sept. 05:  5 tallies added to step functions.  see also: CSTEP.f
C  nov.  05:  add eltot and ve to step function data.
C
C  same as  WRPLAM.f, RPLAM.f, but for XDR format
C
      SUBROUTINE WRPLAM_XDR(TRCFLE,IFLG)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CZT1
      USE COMSOU
      USE CSTEP
      USE COMXS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFLG
      LOGICAL TRCFLE, LHELP(1)
      REAL(DP), ALLOCATABLE :: RDUM(:)
      REAL(DP) :: RHELP(1)
      INTEGER, ALLOCATABLE :: IDUM(:)
      LOGICAL, ALLOCATABLE :: LDUM(:)
      INTEGER :: NRDUM, NIDUM, NLDUM,IHELP(1)
      INTEGER :: IUN, ISPZ
C
      CALL FXDROPN ('fort13','ENCODE',iun)

      CALL FXDRDBL (IUN,TEIN,NRAD)
      CALL FXDRDBL (IUN,TIIN,NPLSTI*NRAD)
      CALL FXDRDBL (IUN,DEIN,NRAD)
      CALL FXDRDBL (IUN,DIIN,NPLS*NRAD)
      CALL FXDRDBL (IUN,VXIN,NPLSV*NRAD)
      CALL FXDRDBL (IUN,VYIN,NPLSV*NRAD)
      CALL FXDRDBL (IUN,VZIN,NPLSV*NRAD)
      CALL FXDRDBL (IUN,BXIN,NRAD)
      CALL FXDRDBL (IUN,BYIN,NRAD)
      CALL FXDRDBL (IUN,BZIN,NRAD)
      CALL FXDRDBL (IUN,BFIN,NRAD)
      CALL FXDRDBL (IUN,ADIN,NAIN*NRAD)
      CALL FXDRDBL (IUN,EDRIFT,NPLS*NRAD)
      CALL FXDRDBL (IUN,VOL,NRAD)
      CALL FXDRDBL (IUN,WGHT,NSPZMC*NRAD)
      CALL FXDRDBL (IUN,BXPERP,NRAD)
      CALL FXDRDBL (IUN,BYPERP,NRAD)
      CALL FXDRDBL (IUN,FLXOUT,NLMPGS)
      CALL FXDRDBL (IUN,SAREA,NLMPGS)
      CALL FXDRDBL (IUN,TEINL,NRAD)
      CALL FXDRDBL (IUN,TIINL,NPLSTI*NRAD)
      CALL FXDRDBL (IUN,DEINL,NRAD)
      CALL FXDRDBL (IUN,DIINL,NPLS*NRAD)
      CALL FXDRDBL (IUN,BVIN,NPLSV*NRAD)
      CALL FXDRDBL (IUN,PARMOM,NPLS*NRAD)
      CALL FXDRDBL (IUN,RMASSI,NION)
      CALL FXDRDBL (IUN,RMASSA,NATM)
      CALL FXDRDBL (IUN,RMASSM,NMOL)
      CALL FXDRDBL (IUN,RMASSP,NPLS)
      CALL FXDRDBL (IUN,DIOD,NION)
      CALL FXDRDBL (IUN,DATD,NATM)
      CALL FXDRDBL (IUN,DMLD,NMOL)
      CALL FXDRDBL (IUN,DPLD,NPLS)
      CALL FXDRDBL (IUN,DPHD,NPHOT)
      CALL FXDRDBL (IUN,DION,NION)
      CALL FXDRDBL (IUN,DATM,NATM)
      CALL FXDRDBL (IUN,DMOL,NMOL)
      CALL FXDRDBL (IUN,DPLS,NPLS)
      CALL FXDRDBL (IUN,DPHOT,NPHOT)
      RHELP(1) = TVAC
      CALL FXDRDBL (IUN,RHELP,1)
      RHELP(1) = DVAC
      CALL FXDRDBL (IUN,RHELP,1)
      RHELP(1) = VVAC
      CALL FXDRDBL (IUN,RHELP,1)
      RHELP(1) = ALLOC
      CALL FXDRDBL (IUN,RHELP,1)

      DO ISPZ=1,NSPZ
        CALL FXDRCHR (IUN,TEXTS(ISPZ))
      END DO

      IHELP(1) = NSPH
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NPHOTI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NPHOTIM
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRINT (IUN,NFOLPH,NPHOT)
      CALL FXDRINT (IUN,NGENPH,NPHOT)
      IHELP(1) = NSPA
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NATMI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NATMIM
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRINT (IUN,NMASSA,NATM)
      CALL FXDRINT (IUN,NCHARA,NATM)
      CALL FXDRINT (IUN,NFOLA,NATM)
      CALL FXDRINT (IUN,NGENA,NATM)
      IHELP(1) = NSPAM
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NMOLI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NMOLIM
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRINT (IUN,NMASSM,NMOL)
      CALL FXDRINT (IUN,NCHARM,NMOL)
      CALL FXDRINT (IUN,NFOLM,NMOL)
      CALL FXDRINT (IUN,NGENM,NMOL)
      IHELP(1) = NSPAMI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NIONI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NIONIM
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRINT (IUN,NMASSI,NION)
      CALL FXDRINT (IUN,NCHARI,NION)
      CALL FXDRINT (IUN,NCHRGI,NION)
      CALL FXDRINT (IUN,NFOLI,NION)
      CALL FXDRINT (IUN,NGENI,NION)
      IHELP(1) = NSPTOT
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NPLSI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NPLSIM
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRINT (IUN,NMASSP,NPLS)
      CALL FXDRINT (IUN,NCHARP,NPLS)
      CALL FXDRINT (IUN,NCHRGP,NPLS)
      IHELP(1) = NBITS
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NSNVI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NCPVI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NADVI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NBGVI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NALVI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NCLVI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NADSI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NALSI
      CALL FXDRINT (IUN,IHELP,1)
      IHELP(1) = NAINI
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRINT (IUN,NPRT,NSPZ)
      CALL FXDRINT (IUN,ISPEZ,6*(NPHOT+2)*(NATM+2)*(NMOL+2)*(NION+2)*
     .              (NPLS+2))
      CALL FXDRINT (IUN,ISPEZI,6*NSPZ)
      CALL FXDRINT (IUN,MPLSTI,NPLS)
      CALL FXDRINT (IUN,MPLSV,NPLS)

      CALL FXDRLOG (IUN,LGVAC,(NPLS+2)*NRAD)
      CALL FXDRLOG (IUN,LGDFT,NRAD)

      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module COMUSR.f '

      CALL CMDTA_XDR (IUN)
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMDTA,ICMDTA'

      CALL CMAMF_XDR (IUN,0)
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMAMF,ICMAMF'

      CALL FXDRDBL (IUN,RCZT1,NZT1)
      CALL FXDRDBL (IUN,RCZT2,NZT2)
      CALL FXDRDBL (IUN,ZT1,NPLS*NRAD)
      CALL FXDRDBL (IUN,ZRG,NPLS*NRAD)
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCZT1,RCZT2,ZT1,ZRG'

      CALL FXDRDBL (IUN,RCMSOU,(12+11*NSRFS)*NSTRA)
      CALL FXDRDBL (IUN,SREC,(NPLS+1)*(NREC+1))
      CALL FXDRDBL (IUN,EIO,(NPLS+1)*(NREC+1))
      CALL FXDRDBL (IUN,EEL,(NPLS+1)*(NREC+1))
      CALL FXDRINT (IUN,ICMSOU,(14+9*NSRFS)*NSTRA)
      CALL FXDRINT (IUN,INGRDA,NSRFS*NSTRA*3)
      CALL FXDRINT (IUN,INGRDE,NSRFS*NSTRA*3)
      IHELP(1) = NSTRAI
      CALL FXDRINT (IUN,IHELP,1)
      CALL FXDRLOG (IUN,LCMSOU,13*NSTRA)
      CALL FXDRLOG (IUN,NLSYMP,NSTRA+1)
      CALL FXDRLOG (IUN,NLSYMT,NSTRA+1)
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: RCMSOU,ICMSOU,LCMSOU'

      IF (ALLOCATED(FLSTEP)) THEN
        CALL FXDRDBL (IUN,FLSTEP,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,ELSTEP,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,FLTOT,(NSPZ+1)*NSTEP)
        CALL FXDRDBL (IUN,ELTOT,(NSPZ+1)*NSTEP)
        CALL FXDRDBL (IUN,VF,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,VE,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,QUOT,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,ADD,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,QUOTI,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,ADDIV,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,TESTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,TISTEP,NPLSTI*NSTEP*NGITT)
        CALL FXDRDBL (IUN,RRSTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,VXSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,VYSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,VZSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,DISTEP,NPLS*NSTEP*NGITT)

c  next five tallies added in sept. 05, !dr
        CALL FXDRDBL (IUN,FESTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,FISTEP,NPLS*NSTEP*NGITT)
        CALL FXDRDBL (IUN,SHSTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,VPSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,MCSTEP,NPLSV*NSTEP*NGITT)

        CALL FXDRINT (IUN,IRSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IPSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,ITSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IASTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IBSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IGSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,ISTUF,NSTEP)
        CALL FXDRINT (IUN,NSMAX,NSTEP)
        CALL FXDRINT (IUN,NSPSTI,NSTEP)
        CALL FXDRINT (IUN,NSPSTE,NSTEP)
        IF (TRCFLE) WRITE (iunout,*) 'WRITE 13: module CSTEP.f'
      END IF

      CALL FXDRCLS (IUN)
      RETURN
C
C
C
      ENTRY RPLAM_XDR(TRCFLE,IFLG)
C
      CALL FXDROPN ('fort13','DECODE',iun)

      CALL FXDRDBL (IUN,TEIN,NRAD)
      CALL FXDRDBL (IUN,TIIN,NPLSTI*NRAD)
      CALL FXDRDBL (IUN,DEIN,NRAD)
      CALL FXDRDBL (IUN,DIIN,NPLS*NRAD)
      CALL FXDRDBL (IUN,VXIN,NPLSTI*NRAD)
      CALL FXDRDBL (IUN,VYIN,NPLSTI*NRAD)
      CALL FXDRDBL (IUN,VZIN,NPLSTI*NRAD)
      CALL FXDRDBL (IUN,BXIN,NRAD)
      CALL FXDRDBL (IUN,BYIN,NRAD)
      CALL FXDRDBL (IUN,BZIN,NRAD)
      CALL FXDRDBL (IUN,BFIN,NRAD)
      CALL FXDRDBL (IUN,ADIN,NAIN*NRAD)
      CALL FXDRDBL (IUN,EDRIFT,NPLS*NRAD)
      CALL FXDRDBL (IUN,VOL,NRAD)
      CALL FXDRDBL (IUN,WGHT,NSPZMC*NRAD)
      CALL FXDRDBL (IUN,BXPERP,NRAD)
      CALL FXDRDBL (IUN,BYPERP,NRAD)
      CALL FXDRDBL (IUN,FLXOUT,NLMPGS)
      CALL FXDRDBL (IUN,SAREA,NLMPGS)
      CALL FXDRDBL (IUN,TEINL,NRAD)
      CALL FXDRDBL (IUN,TIINL,NPLSTI*NRAD)
      CALL FXDRDBL (IUN,DEINL,NRAD)
      CALL FXDRDBL (IUN,DIINL,NPLS*NRAD)
      CALL FXDRDBL (IUN,BVIN,NPLSV*NRAD)
      CALL FXDRDBL (IUN,PARMOM,NPLS*NRAD)
      CALL FXDRDBL (IUN,RMASSI,NION)
      CALL FXDRDBL (IUN,RMASSA,NATM)
      CALL FXDRDBL (IUN,RMASSM,NMOL)
      CALL FXDRDBL (IUN,RMASSP,NPLS)
      CALL FXDRDBL (IUN,DIOD,NION)
      CALL FXDRDBL (IUN,DATD,NATM)
      CALL FXDRDBL (IUN,DMLD,NMOL)
      CALL FXDRDBL (IUN,DPLD,NPLS)
      CALL FXDRDBL (IUN,DPHD,NPHOT)
      CALL FXDRDBL (IUN,DION,NION)
      CALL FXDRDBL (IUN,DATM,NATM)
      CALL FXDRDBL (IUN,DMOL,NMOL)
      CALL FXDRDBL (IUN,DPLS,NPLS)
      CALL FXDRDBL (IUN,DPHOT,NPHOT)
      CALL FXDRDBL (IUN,RHELP,1)
      TVAC = RHELP(1)
      CALL FXDRDBL (IUN,RHELP,1)
      DVAC = RHELP(1)
      CALL FXDRDBL (IUN,RHELP,1)
      VVAC = RHELP(1)
      CALL FXDRDBL (IUN,RHELP,1)
      ALLOC = RHELP(1)

      DO ISPZ=1,NSPZ
        CALL FXDRCHR (IUN,TEXTS(ISPZ))
      END DO

      CALL FXDRINT (IUN,IHELP,1)
      NSPH = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NPHOTI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NPHOTIM = IHELP(1)
      CALL FXDRINT (IUN,NFOLPH,NPHOT)
      CALL FXDRINT (IUN,NGENPH,NPHOT)
      CALL FXDRINT (IUN,IHELP,1)
      NSPA = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NATMI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NATMIM = IHELP(1)
      CALL FXDRINT (IUN,NMASSA,NATM)
      CALL FXDRINT (IUN,NCHARA,NATM)
      CALL FXDRINT (IUN,NFOLA,NATM)
      CALL FXDRINT (IUN,NGENA,NATM)
      CALL FXDRINT (IUN,IHELP,1)
      NSPAM = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NMOLI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NMOLIM = IHELP(1)
      CALL FXDRINT (IUN,NMASSM,NMOL)
      CALL FXDRINT (IUN,NCHARM,NMOL)
      CALL FXDRINT (IUN,NFOLM,NMOL)
      CALL FXDRINT (IUN,NGENM,NMOL)
      CALL FXDRINT (IUN,IHELP,1)
      NSPAMI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NIONI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NIONIM = IHELP(1)
      CALL FXDRINT (IUN,NMASSI,NION)
      CALL FXDRINT (IUN,NCHARI,NION)
      CALL FXDRINT (IUN,NCHRGI,NION)
      CALL FXDRINT (IUN,NFOLI,NION)
      CALL FXDRINT (IUN,NGENI,NION)
      CALL FXDRINT (IUN,IHELP,1)
      NSPTOT = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NPLSI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NPLSIM = IHELP(1)
      CALL FXDRINT (IUN,NMASSP,NPLS)
      CALL FXDRINT (IUN,NCHARP,NPLS)
      CALL FXDRINT (IUN,NCHRGP,NPLS)
      CALL FXDRINT (IUN,IHELP,1)
      NBITS = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NSNVI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NCPVI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NADVI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NBGVI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NALVI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NCLVI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NADSI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NALSI = IHELP(1)
      CALL FXDRINT (IUN,IHELP,1)
      NAINI = IHELP(1)
      CALL FXDRINT (IUN,NPRT,NSPZ)
      CALL FXDRINT (IUN,ISPEZ,6*(NPHOT+2)*(NATM+2)*(NMOL+2)*(NION+2)*
     .              (NPLS+2))
      CALL FXDRINT (IUN,ISPEZI,6*NSPZ)
      CALL FXDRINT (IUN,MPLSTI,NPLS)
      CALL FXDRINT (IUN,MPLSV,NPLS)

      CALL FXDRLOG (IUN,LGVAC,(NPLS+2)*NRAD)
      CALL FXDRLOG (IUN,LGDFT,NRAD)

      IF (TRCFLE) WRITE (iunout,*) 'READ 13: module COMUSR.f '

      CALL CMDTA_XDR (IUN)
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMDTA,ICMDTA'

      CALL CMAMF_XDR (IUN,1)
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMAMF,ICMAMF'

      CALL FXDRDBL (IUN,RCZT1,NZT1)
      CALL FXDRDBL (IUN,RCZT2,NZT2)
      CALL FXDRDBL (IUN,ZT1,NPLS*NRAD)
      CALL FXDRDBL (IUN,ZRG,NPLS*NRAD)
      IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCZT1,RCZT2,ZT1,ZRG'

      IF (IFLG == 0) THEN
        CALL FXDRDBL (IUN,RCMSOU,(12+11*NSRFS)*NSTRA)
        CALL FXDRDBL (IUN,SREC,(NPLS+1)*(NREC+1))
        CALL FXDRDBL (IUN,EIO,(NPLS+1)*(NREC+1))
        CALL FXDRDBL (IUN,EEL,(NPLS+1)*(NREC+1))
        CALL FXDRINT (IUN,ICMSOU,(14+9*NSRFS)*NSTRA)
        CALL FXDRINT (IUN,INGRDA,NSRFS*NSTRA*3)
        CALL FXDRINT (IUN,INGRDE,NSRFS*NSTRA*3)
        CALL FXDRINT (IUN,IHELP,1)
        NSTRAI = IHELP(1)
        CALL FXDRLOG (IUN,LCMSOU,13*NSTRA)
        CALL FXDRLOG (IUN,NLSYMP,NSTRA+1)
        CALL FXDRLOG (IUN,NLSYMT,NSTRA+1)
        IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMSOU,ICMSOU,LCMSOU'
      ELSE
        NRDUM = MAX(SIZE(RCMSOU), SIZE(SREC), SIZE(EIO), SIZE(EEL))
        NIDUM = MAX(SIZE(ICMSOU), SIZE(INGRDA), SIZE(INGRDE))
        NLDUM = MAX(SIZE(LCMSOU), SIZE(NLSYMP), SIZE(NLSYMT))
        ALLOCATE (RDUM(NRDUM))
        ALLOCATE (IDUM(NIDUM))
        ALLOCATE (LDUM(NLDUM))
        LDUM = .FALSE.
        CALL FXDRDBL (IUN,RDUM,(12+11*NSRFS)*NSTRA)
        CALL FXDRDBL (IUN,RDUM,(NPLS+1)*(NREC+1))
        CALL FXDRDBL (IUN,RDUM,(NPLS+1)*(NREC+1))
        CALL FXDRDBL (IUN,RDUM,(NPLS+1)*(NREC+1))
        CALL FXDRINT (IUN,IDUM,(14+9*NSRFS)*NSTRA)
        CALL FXDRINT (IUN,IDUM,NSRFS*NSTRA*3)
        CALL FXDRINT (IUN,IDUM,NSRFS*NSTRA*3)
        CALL FXDRINT (IUN,IDUM,1)
        CALL FXDRLOG (IUN,LDUM,13*NSTRA)
        CALL FXDRLOG (IUN,LDUM,NSTRA+1)
        CALL FXDRLOG (IUN,LDUM,NSTRA+1)
        DEALLOCATE (RDUM)
        DEALLOCATE (IDUM)
        DEALLOCATE (LDUM)
      END IF

      IF (ALLOCATED(FLSTEP)) THEN
        CALL FXDRDBL (IUN,FLSTEP,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,ELSTEP,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,FLTOT,(NSPZ+1)*NSTEP)
        CALL FXDRDBL (IUN,ELTOT,(NSPZ+1)*NSTEP)
        CALL FXDRDBL (IUN,VF,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,VE,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,QUOT,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,ADD,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,QUOTI,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,ADDIV,(NSPZ+1)*NSTEP*NGITT)
        CALL FXDRDBL (IUN,TESTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,TISTEP,NPLSTI*NSTEP*NGITT)
        CALL FXDRDBL (IUN,RRSTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,VXSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,VYSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,VZSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,DISTEP,NPLS*NSTEP*NGITT)

c  next five tallies added in sept. 05, !dr
        CALL FXDRDBL (IUN,FESTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,FISTEP,NPLS*NSTEP*NGITT)
        CALL FXDRDBL (IUN,SHSTEP,NSTEP*NGITT)
        CALL FXDRDBL (IUN,VPSTEP,NPLSV*NSTEP*NGITT)
        CALL FXDRDBL (IUN,MCSTEP,NPLSV*NSTEP*NGITT)

        CALL FXDRINT (IUN,IRSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IPSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,ITSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IASTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IBSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,IGSTEP,NSTEP*NGITT)
        CALL FXDRINT (IUN,ISTUF,NSTEP)
        CALL FXDRINT (IUN,NSMAX,NSTEP)
        CALL FXDRINT (IUN,NSPSTI,NSTEP)
        CALL FXDRINT (IUN,NSPSTE,NSTEP)
        IF (TRCFLE) WRITE (iunout,*) 'READ 13: module CSTEP.f'
      END IF

      CALL FXDRCLS (IUN)
      RETURN
      END
C ===== SOURCE: wrrec.f
C
      SUBROUTINE WRREC
C
C  EVALUATE EIRENE RECOMMENDATIONS FOR A NEXT RUN OF THE SAME MODEL
C
      USE PRECISION
      USE PARMMOD
      USE CAI
      USE CCONA
      USE CTRCEI
      USE COMSOU
      USE COMPRT, ONLY: IUNOUT
      USE COUTAU

      IMPLICIT NONE

      REAL(DP) :: WSUM, FTOT, XNSUM
      INTEGER :: NREQ, ISTRA
      REAL(DP) :: WTOTT(NSTRA),WMEAN(NSTRA),WREC(NSTRA),XNEXP(NSTRA),
     .          CPUFAC(NSTRA)

      OPEN (UNIT=14,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 14
C
C  FIRSTLY: STRATIFIED SOURCE SAMPLING
C
      IF (NSTRAI.EQ.1) THEN
        RATIO(1)=1.
        GOTO 350
      ENDIF
C
      WSUM=0.
      NREQ=0
      FTOT=0.
      DO 100 ISTRA=1,NSTRAI
        IF (XMCP(ISTRA).LE.0.D0) GOTO 100
        WTOTT(ISTRA)=-WTOTP(0,ISTRA)+WTOTA(0,ISTRA)+
     .                WTOTM(0,ISTRA)+WTOTI(0,ISTRA)
        WTOTT(ISTRA)=WTOTT(ISTRA)/(FLXFAC(ISTRA)+EPS60)
        WMEAN(ISTRA)=WTOTT(ISTRA)/(XMCP(ISTRA)+EPS60)
        WSUM=WSUM+WTOTT(ISTRA)
        NREQ=NREQ+NPTS(ISTRA)
        FTOT=FTOT+FLUXT(ISTRA)
100   CONTINUE
C
C  PROPORTIONAL ALLOCATION: RECOMMENDED REL. WEIGHT PER STRATUM: WREC
C                           EXPECTED REL. NO OF PARTICLES NEEDED: XNEXP
C                           RECOMMENDED NO. OF PARTICLES: NRECOM
C  NOT THE NO. OF PARTICLES BUT THE SUM OF BIRTH WEIGHTS PER STRATUM
C  IS ALLOCATED  PROPORTIONAL TO THE RELATIVE STRATUM POPULATION
C
C
      XNSUM=0.
      DO 200 ISTRA=1,NSTRAI
        IF (XMCP(ISTRA).LE.0.D0) GOTO 200
        WREC(ISTRA)=FLUXT(ISTRA)/(FTOT+EPS60)
        XNEXP(ISTRA)=WREC(ISTRA)/(WMEAN(ISTRA)+EPS60)
C  ACCOUNT FOR DIFFERENT COMPUTING SPEED AT DIFFERENT STRATA
C  ASSUME THEREFORE THAT ALLOCATED CPU TIME WAS PROPORTIONAL NPTS(ISTRA)
        CPUFAC(ISTRA)=XMCP(ISTRA)/(DBLE(NPTS(ISTRA))+EPS60)
        XNEXP(ISTRA)=XNEXP(ISTRA)/(CPUFAC(ISTRA)+EPS60)
        XNSUM=XNSUM+XNEXP(ISTRA)
200   CONTINUE
      DO 300 ISTRA=1,NSTRAI
        RATIO(ISTRA)=0.
        IF (XMCP(ISTRA).LE.0.D0) GOTO 300
C  SCALE XNEXP TO CONSERVE TOTAL NUMBER OF REQUESTED TRACKS
        XNEXP(ISTRA)=XNEXP(ISTRA)*DBLE(NREQ)/(XNSUM+EPS60)
C  CONVERT XNEXP TO AN INTERGER
        NRECOM(ISTRA)=IDINT(REAL(XNEXP(ISTRA),KIND(1.D0)))
        IF (XNEXP(ISTRA)-DBLE(NRECOM(ISTRA)).GT.0.5)
     .      NRECOM(ISTRA)=NRECOM(ISTRA)+1
        RATIO(ISTRA)=WREC(ISTRA)/(WTOTT(ISTRA)/(WSUM+EPS60)+EPS60)
300   CONTINUE
C
350   CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 14: RATIO,NRECOM '
      WRITE (14) RATIO,NRECOM
C
      IF (.NOT.TRCREC.OR.NSTRAI.EQ.1) GOTO 1000
      CALL PAGE
      WRITE (iunout,*) 
     .  '=================================================='
      WRITE (iunout,*) 
     .  '= RECOMMENDED INPUT MODIFICATIONS FOR STRATIFIED ='
      WRITE (iunout,*) 
     .  '= SOURCE SAMPLING                                ='
      WRITE (iunout,*) 
     .  '=================================================='
      CALL LEER(3)
      WRITE (iunout,*) 'NPTS OLD = OLD INPUT VALUE FOR NPTS '
      WRITE (iunout,*) 'NPTS REC = EIRENE RECOMMENDATION FOR NEXT RUN '
      WRITE (iunout,*) 
     .  'RATIO    = RECOM. REL.WEIGHT / ACTUAL REL. WEIGHT'
      WRITE (iunout,*) 'VALUES OF RATIO CLOSE TO ONE MEAN '
      WRITE (iunout,*) '"ALMOST PROPORTIONAL ALLOCATION OF WEIGHTS"'
      CALL LEER(2)
      WRITE (iunout,*) 'STRAT. NO, NPTS OLD, NPTS REC, RATIO '
      DO 400 ISTRA=1,NSTRAI
        WRITE (iunout,'(1X,I2,8X,I6,4X,I6,5X,1P,E12.4)')
     .                ISTRA,NPTS(ISTRA),NRECOM(ISTRA),RATIO(ISTRA)
400   CONTINUE
      CALL LEER(2)
1000  CONTINUE
C
C  STRATIFIED SOURCE SAMPLING ACCESSMENT FINISHED
C
C  SECONDLY: WEIGHT WINDOWS
C
C  TO BE WRITTEN
C
      RETURN
C
      ENTRY RREC
C
      OPEN (UNIT=14,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 14
      IF (TRCFLE) WRITE (iunout,*) 'READ 14: RATIO,NRECOM '
      READ (14) RATIO,NRECOM
C
      RETURN
      END
C ===== SOURCE: wrsnap.f
!pb  26.10.06: close file after read or write
!pb  31.10.06:  definition of RPART, RPARTC, IPART, IPARTC changed
!               RPART(NPARTT,NPRNL) --> RPART(NPRNL,NPARTT)
C
C
      SUBROUTINE WRSNAP
C
C  SAVE SNAPSHOT POPULATION AT END OF TIMESTEP
C
      USE PRECISION
      USE PARMMOD
      USE CTRCEI
      USE COMPRT
      USE COMNNL
      USE COMSOU

      IMPLICIT NONE

      INTEGER :: I, J
C
      OPEN (UNIT=15,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15
C
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 15: IPRNL,FLUX,DTIMV '
      WRITE (15) IPRNL,FLUX(NSTRAI),DTIMV
      WRITE (15) ((RPARTC(J,I),J=1,NPARTT),I=1,IPRNL)
      WRITE (15)  (RPARTW(I)              ,I=0,IPRNL)
      WRITE (15) ((IPARTC(J,I),J=1,MPARTT),I=1,IPRNL)
      CLOSE (UNIT=15)
C
      RETURN
C
      ENTRY RSNAP
C
      OPEN (UNIT=15,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15
      READ (15) IPRNL,FLUX(NSTRAI),DTIMV
      IF (TRCFLE) WRITE (iunout,*) 'READ 15: IPRNL,FLUX,DTIMV '
      READ (15) ((RPARTC(J,I),J=1,NPARTT),I=1,IPRNL)
      READ (15)  (RPARTW(I)              ,I=0,IPRNL)
      READ (15) ((IPARTC(J,I),J=1,MPARTT),I=1,IPRNL)
      CLOSE (UNIT=15)
C
      RETURN
      END
