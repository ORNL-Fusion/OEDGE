      subroutine read_hydkin (ir,filename,h123,reac,crc,rmn,rmx,
     .                        e_el,e_k,lffl)

      use precision
      use parmmod
      use comxs
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: ir
      character(len=*), intent(in) :: reac, filename
      logical, intent(in) :: lffl
      character(4), intent(inout) :: h123
      character(3), intent(inout) :: crc
      real(dp) , intent(out) :: rmn, rmx, e_el, e_k
      character(132) :: zeile
      character(12) :: chr
      character(len=len(reac)+10) :: cpreac
      integer :: ianf, iend, ier, ll, io, ie, iflg
      type(hydkin_data), pointer :: hp

      open (unit=28,file=filename)

      zeile = repeat(' ',len(zeile))

! find number of temperatures
      do while (index(zeile,'Default energy mesh') == 0)
         read (28,'(A132)') zeile
      end do
      
      allocate (hp)
      hp%reacname = reac
      
      read (28,'(A132)') zeile
      ianf = index(zeile,'=')
      read (zeile(ianf+1:),*) hp%ntemps

      do while (index(zeile,'EeVDef') == 0)
        read (28,'(A132)') zeile
      end do
      
      allocate (hp%temps(hp%ntemps))
      allocate (hp%rates(hp%ntemps))
      allocate (hp%ratio(hp%ntemps))

! read energies
      do ie=1, hp%ntemps
        read (28,*) hp%temps(ie)
      end do

      rmn = hp%temps(1)
      rmx = hp%temps(hp%ntemps)

      e_el = 0._dp
      e_k = 0._dp

! find specified reaction

      cpreac = 'Eirname = '//adjustl(reac)
      ll = len_trim(cpreac)

      zeile = repeat(' ',len(zeile))
      do while (index(zeile,cpreac) == 0)
        read (28,'(A132)',iostat=io) zeile
      end do

      if (io .ne. 0) then
        write (iunout,*) ' ERROR READING REACTION FROM HYDKIN DATABASE '
        write (iunout,*) ' FILE IS ',filename
        write (iunout,*) ' REACTION IS ',reac
        call exit_own(1)
      end if

! reaction found 
! now read data
      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'E_el') == 0)
        read (28,'(A132)',iostat=io) zeile
      end do
      ianf = index(zeile,'=')+1
      read (zeile(ianf:),*) e_el

      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'E_K') == 0)
        read (28,'(A132)',iostat=io) zeile
      end do
      ianf = index(zeile,'=')+1
      read (zeile(ianf:),*) e_k

      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'RPrT') == 0)
        read (28,'(A132)',iostat=io) zeile
      end do
      
      ianf = index(zeile,'''')+1
      iend = ianf-1 + index(zeile(ianf:),'''') -1
      hp%rprt = zeile(ianf:iend)

      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'RName') == 0)
        read (28,'(A132)',iostat=io) zeile
      end do
      
      ianf = index(zeile,'''')+1
      iend = ianf-1 + index(zeile(ianf:),'''') -1
      hp%reac_string = zeile(ianf:iend)

      zeile = repeat(' ',len(zeile))
      do while (index(zeile,'RData') == 0)
        read (28,'(A132)',iostat=io) zeile
      end do

! read rates
      do ie=1, hp%ntemps
        read (28,*) hp%rates(ie)
      end do
    
      close (unit=28)

      do ie=1, hp%ntemps-1
        hp%ratio(ie) = (hp%rates(ie+1) - hp%rates(ie)) /
     .                 (hp%temps(ie+1) - hp%temps(ie))
      end do 

      if (reacdat(ir)%lrtc) then
        write (iunout,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                   ' FOR REACTION ',ir
        write (iunout,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
        deallocate (hp)
        call exit_own(1)
      end if 

      reacdat(ir)%lrtc = .true.
      allocate(reacdat(ir)%rtc)
      nullify (reacdat(ir)%rtc%adas)
      nullify (reacdat(ir)%rtc%line)
      nullify (reacdat(ir)%rtc%poly)
      reacdat(ir)%rtc%hyd => hp
      reacdat(ir)%rtc%ifit = 4

      if (lffl) then
        h123 = 'H.2 '
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
C  DEFAULT RATE COEFFICIENT: 8TH ORDER POLYNOM OF LN(<SIGMA V>) FOR E0=0.
        IFTFLG(IR,IFLG)=0
      
! find crc

        if (index(hp%rprt,'CX') > 0) then
           crc = 'CX '
           iswr(ir) = 3
        else if (index(hp%rprt,'_R') + index(hp%rprt,'R-DR') > 0) then
           crc = 'RC '
           iswr(ir) = 6
        else if (index(hp%rprt,'_DE') + index(hp%rprt,'_DI') + 
     .          index(hp%rprt,'_I') + index(hp%rprt,'I-DI') +
     .          index(hp%rprt,'_CAD') > 0) then
           crc = 'EI '
           iswr(ir) = 1
        else
           write (iunout,*) ' UNKNOWN REACTION TYPE ',HP%RPRT 
           write (iunout,*) ' USED IN REACTION ',reac
           write (iunout,*) ' PLEASE CHECK SPECIFICATION OF REACTIONS'
           call exit_own(1)
        end if
      endif        

      return
      end subroutine read_hydkin
