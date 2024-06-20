      subroutine EIRENE_read_photdbk (ir, reac, isw)
!  16.2.05:  write statement taken out
!  2.11.05:  database handle introduced for file POLARI
      use EIRMOD_precision
      use EIRMOD_parmmod
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_CINIT
      USE EIRMOD_PHOTON
 
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
      character(2) :: kenn(12)
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
          CALL EIRENE_EXIT_OWN(1)
        END IF
 
      if (ifrst == 0) then
         ifrst = 1
         call EIRENE_read_polari(polari_elnam,polari_fac,npolari)
      end if
 
      lr=len_trim(reac)
 
      read (29,*)
      read (29,*)
      read (29,*)
 
      do
        read (29,'(A1000)',end=990) zeile
        if (zeile(1:2) == '--') cycle
        call EIRENE_subcomma(zeile)
 
!  read Element
        ianf = 3
        iend = ianf + scan(zeile(ianf:),'|') - 1
 
        call EIRENE_delete_blanks(zeile(ianf:iend))
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
      kenn = '  '
      do i=1,ifremd
        read (iunin,'(i6,1x,a2,3x,i6)') ii,kenn(i),ik6
        call EIRENE_uppercase(kenn(i))
        if (kenn(i) == 'QS') then
          ipc6(i) = ik6
          c6a(i) = c6_qs_mess
        else
          do j=1,npolari
            if (kenn(i) == polari_elnam(j)) then
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
            call EIRENE_exit_own(1)
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
      phline%kenn(1:12) = kenn(1:12)
      phline%ifremd = count(phline%iplsc6(1:12)>0)
 
 
      phline%ignd = iplsc3
      phline%iprofiletype = iprftype
      phline%imess = imess
 
      reacdat(ir)%lphr = .true.
 
      allocate (reacdat(ir)%phr)
      allocate (reacdat(ir)%phr%line)
      reacdat(ir)%phr%ifit = -1
      reacdat(ir)%phr%line => phline
 
      call EIRENE_get_reaction(ir)
      b21=EIRENE_ph_b21()
      b12=b21*gi/gj
      phline%b21 = b21
      phline%b12 = b12
      phline%nrjprt = nrjprt
 
      phline%reacname = reac_name(ir)
 
      modclf(ir) = 100
      iftflg(ir,2) = 110
 
      rdata = 0._dp
      rdata(1,1) = aik
 
      call EIRENE_set_reaction_data
     .  (ir,isw,iftflg(ir,2),rdata,iunout,.false.)
 
      return
 
 990  continue
      write (iunout,*) 'REACTION ',reac,' NOT FOUND IN FILE PHOTON'
      call EIRENE_exit_own(1)
 
 
      contains
 
 
      subroutine EIRENE_subcomma (str)
      implicit none
      character(len=*), intent(in out) :: str
      integer :: i
 
      do
        i=scan(str,',')
        if (i == 0) exit
        str(i:i)='.'
      end do
      return
      end subroutine EIRENE_subcomma
 
 
      subroutine EIRENE_delete_blanks (str)
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
      end subroutine EIRENE_delete_blanks
 
 
      subroutine EIRENE_read_polari (name, factor, n)
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
        CALL EIRENE_EXIT_OWN(1)
      END IF
 
      open (unit=28,file=DBFNAME(IFILE),access='sequential',
     .      form='formatted')
 
      n=0
      read (28,*)
      read (28,*)
 
      do
        read (28,'(A200)',end=990) zeile
        if (zeile(1:2) == '--') cycle
        call EIRENE_subcomma(zeile)
 
        n = n + 1
        if (n > nmax) then
           write (iunout,*) ' too many lines in polarization file '
           write (iunout,*) ' calculation stopped '
           call EIRENE_exit_own(1)
        end if
 
        ianf = 3
        iend = ianf + scan(zeile(ianf:),'|') - 1
 
!  read element name
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        call EIRENE_delete_blanks(zeile(ianf:iend))
        iblnk = scan(zeile(ianf:iend),'|') - 1
        if (iblnk < 0) iblnk = iend-ianf+1
        if (iblnk == 0) then
           write (iunout,*) ' ERROR IN DATABASE POLARI'
           write (iunout,*) ' NO ELEMENTNAME FOUND '
           stop
        end if
 
        name(n) = repeat(' ',2)
        name(n)(1:iblnk) = zeile(ianf:ianf+iblnk-1)
        call EIRENE_uppercase(name(n))
 
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
      end subroutine EIRENE_read_polari
 
      end subroutine EIRENE_read_photdbk
 
