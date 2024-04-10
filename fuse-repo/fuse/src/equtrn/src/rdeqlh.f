      subroutine rdeqlh(lun,nr,nz,btf,rtf,*)
c
c  version : 16.09.2000 00:19
c
c=======================================================================
c*** Reads the header of the equilibrium file and returns dimensions
c*** and the value of the toroidal field, btf, measured at radius rtf -
c*** provided that this information is present in the equilibrium file.
c
      implicit none
      integer lun,nr,nz
      real btf,rtf
      integer i,j,k,l
      character ss*80, s(80), hh*80, h(80), utb
      equivalence (s,ss),(h,hh)
      logical ll
c=======================================================================
c
      utb=char(9)
      rtf=-1.
 10   read(lun,'(a)',end=90) ss
      if(ss.eq.' ') go to 10
      k=0
      hh=' '
      do 100 i=1,80
          if(s(i).eq.'(') then
              if(hh.eq.'r') return
              go to 90
          end if
          if(s(i).eq.':') go to 10
          if(s(i).eq.' ' .or. s(i).eq.utb) go to 100
          if(s(i).eq.'=') go to 20
          k=k+1
          h(k)=s(i)
 100  continue

      if(hh.eq.'r(1:jm);') return
      go to 90
c-----------------------------------------------------------------------
c*** Parse the string
c
 20   l=i+1
      if(hh.eq.'jm') then
          j=1
      else if(hh.eq.'km') then
          j=2
      else if(hh.eq.'btf') then
          j=3
      else if(hh.eq.'rtf') then
          j=4
      else
          go to 10
      end if
      k=0
      hh=' '
      ll=.false.
      do i=l,80 !{
        if(s(i).eq.';') go to 30
        if(s(i).eq.' ' .or. s(i).eq.utb) then !{
          if(ll) go to 30
        else !}{
          ll=.true.
          k=k+1
          h(k)=s(i)
        end if !}
      end do !}

 30   go to (31,32,33,34),j
c-----------------------------------------------------------------------
 31   read(hh,'(bn,i4)') nr
      go to 10
 32   read(hh,'(bn,i4)') nz
      go to 10
 33   read(hh,'(bn,e12.0)') btf
      go to 10
 34   read(hh,'(bn,e12.0)') rtf
      go to 10
c=======================================================================
c*** Error encountered
c
 90   write (6,*) 'Wrong format of the equilibrium file - sorry!'
      return 1
c=======================================================================
      end
