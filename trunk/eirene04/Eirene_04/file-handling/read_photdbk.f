      subroutine read_photdbk (ir, reac, isw)

      use precision
      use parmmod
      USE COMXS
      USE COMPRT
      USE CCONA
      USE PHOTON

      implicit none

      integer, intent(in) :: ir, isw
      CHARACTER(50), INTENT(IN) :: REAC

      real(dp) :: wl, aik, ei, ej, c2, c3, c4, c6_theo, b12, b21,
     .            c3_theo, c6_qs_mess, c6_ar_mess, c3_mess, c6
      real(dp) :: c6a(12)
      real(dp), save :: polari_fac(120)
      integer :: gi, gj
      integer :: ianf, iend, iblnk, lr, ic, i, j, iplsc3, iprftype,
     .           imess, ifremd, ii, i1, lel
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

     
!pb      open (unit=29,file='PHOTON',form='FORMATTED',access='SEQUENTIAL')
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
           write (6,*) ' ERROR IN DATABASE PHOTON'
           write (6,*) ' NO ELEMENTNAME FOUND '
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

!        write (elementname(iblnk+1:),'(f8.4,a1)') wl,cha
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

      read (iunin,'(12i6)') iprftype, iplsc3, imess, ifremd

      if (imess == 1) then
         c3 = c3_mess
         c6 = c6_ar_mess
      else
         c3 = c3_theo
         c6 = c6_theo
      end if

      if (ifremd > 12) then
         write (6,*) ' too many fremddruckverbreiterungen specified'
         write (6,*) ' calculation abandonned '
      end if

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
	 end if
      end do
      
      c2=0._dp
      creac(:,:,ir) = 0._dp
      creac(1,1,ir) = aik
!     wl is in nm
      creac(2,1,ir) = hpcl / wl *1.E7_DP
      creac(3,1,ir) = gj
      creac(4,1,ir) = gi
!     Ej is in [1/cm]
      creac(5,1,ir) = ej * cspeed*hplanck*erg_to_ev
      creac(6,1,ir) = c2
      creac(7,1,ir) = c3
      creac(8,1,ir) = c4
      
      select case (isw)
        case (1)    ! absorbtion
           creac(9,1,ir) = 4
        case (2)    ! emission
           creac(9,1,ir) = 4
        case (3)    ! stimulated emission
        case (4)    ! ph_cabs ot (case(7) bei sven)
           creac(9,1,ir) = 7
        case default
           creac(9,1,ir) = 0
      end select

      creac(1:9,2,ir) = c6a(1:9)
      creac(1:3,3,ir) = c6a(10:12)
      creac(4:9,3,ir) = ipc6(1:6)
      creac(1:6,4,ir) = ipc6(7:12)

      creac(7,4,ir) = iplsc3
      creac(8,4,ir) = iprftype

      call get_reaction(ir)
      b21=ph_b21()
      b12=b21*gi/gj
      creac(9,4,ir) = b21
      creac(1,5,ir) = b12

      modclf(ir) = 100
      iftflg(ir,2) = 110

      write (6,*) reac, wl, aik, gi, gj, ei, ej, c3, c4, c6a
      return

 990  continue
      write (6,*) 'REACTION ',reac,' NOT FOUND IN FILE PHOTON'
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
      integer :: nmax
      
      nmax = min(size(name),size(factor))

      open (unit=28,file='POLARI',access='sequential',form='formatted')

      n=0
      read (28,*)
      read (28,*)

      do 
        read (28,'(A200)',end=990) zeile
        if (zeile(1:2) == '--') cycle
        call subcomma(zeile)

        n = n + 1
        if (n > nmax) then
           write (6,*) ' too many lines in polarization file '
           write (6,*) ' calculation stopped '
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
           write (6,*) ' ERROR IN DATABASE POLARI'
           write (6,*) ' NO ELEMENTNAME FOUND '
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

