      subroutine mcrates(tmpe,tmpi,ni,za,zamax,zn,rion,rrec,rcxr,opt)
      use mod_params
      use mod_adpak_com
      implicit none
      real tmpe,tmpi,ni
      integer za,zamax,zn,opt
      real rion, rrec, rcxr
c
c     MCRATES: This routine has been modified and added to DIVIMP to
c              allow DIVIMP to utilize the ADPAK formatted (Braams)
c              atomic physics data.
c
c
c ... Inputs are:
c        electron temperature, tmpe;
c        **neutral hydrogen temperature, tmpi;
c          - replaced by ion temperature for now - divided by ion mass
c
c        background ion density, ni;
c
c        atomic charge state, za;
c        maximum atomic charge state, zamax;
c        nuclear charge state, zn;
c
c        opt - switch allowing for density indexing to be turned
c              on or off.
c
c
c ... Outputs are:
c        rate parameters (sigma*v) for ionization, recombination,
c        and charge-exchange recombination on neutral hydrogen
c
c     The tables used in this subroutine are generated with a code
c     supplied by Bas Braams.  The data file it produces is called
c     'b2frates' by default.  This gives rates that may depend on both
c     density and temperature.  Here, we use only the rates given
c     for the lowest density in the table.
c
c     Input electron temperature is given in [eV].
c     Input temperature for neutral hydrogen is given in [eV/AMU].
c     ** replaced by hydrogen ion temperature/AMU - estimates velocity of
c        approach.
c
c
c     Table temperatures are given in [eV].
c     Output and table rates are all given in [m**3/s].
c
c     The tables are read in separately and stored in common blocks for
c     ease of access. To conserve storage this methodology can be
c     changed to perform a file look up each time at the cost of
c     execution speed and disk I/O.
c
c
c     include 'params'
c     include 'adpak_com'
c
c
c      Use(Share)               # cutlo
c      Use(Physical_constants)  # ev
c      Use(Multicharge)         # rtlt,rtlsa,rtlra,rtlcx
c
      integer i1e,i1i,i,ii,i2
      real tmpenonz,tmpinonz,xte,xti,fxte,fxti,ninonz
      real dlogt,dlogn,xni,fxni
c
c     SET opt to zero for now
c
c
c     jdemod - opt is an intent(in) value - some routines call this with
c              a constant - this routine should not change the passed
c              in value 
c
c      opt = 0
c
      rion = 0.
      rrec = 0.
      rcxr = 0.
c
c     Avoid non-zero values of the arguments
c
      tmpenonz = max(tmpe, lo)
      tmpinonz = max(tmpi, lo)
      ninonz   = max(ni, lo)
c
c     Remove conversion for now - DIVIMP works directly in eV
c
c      xte = log(tmpenonz/ev)
c      xti = log(tmpinonz/ev)
c
      xte = log(tmpenonz)
      xti = log(tmpinonz)
      xni = log(ninonz)
c
      dlogt = rtlt(1) - rtlt(0)
      dlogn = rtln(1) - rtln(0)
c
c ... Find index i1 in temperature table such that
c                 rtlt(i1) .le. xt .lt. rtlt(i1+1)
c     or, equivalently,
c                  rtt(i1) .le. tmp .lt. rtt(i1+1).
c
      i1e = int( (xte-rtlt(0))/dlogt )
      i1i = int( (xti-rtlt(0))/dlogt )
c
c ... For temperatures below minimum table temperature, extrapolate
c     downwards from table entries 0 and 1.
c
      i1e = max(0, i1e)
      i1i = max(0, i1i)
c
c
c ... For temperatures above maximum table temperature, extrapolate
c     upwards from table entries rtnt-1 and rtnt.
c
      i1e = min(rtnt-1, i1e)
      i1i = min(rtnt-1, i1i)
c
c     Find index i2 in density table where
c                 rtln(i2) .le. xni .lt. rtln(i2+1)
c     or, equivalently,
c                  rtn(i2) .le. tmp .lt. rtn(i2+1).
c
      i2 = int( (xni-rtln(0))/dlogn )
c slmod begin - new
c ARRAY BOUNDS: rtln below has dimensions (0:MAXRTNT), but i2 assigned
c               to -1 in at least one instance
      i2 = max(0, i2)
c
c      i2 = max(-1, i2)
c slmod end
      i2 = min(rtnn,i2)
c
c ... Compute coefficient for linear interpolation.
c
c     Extrapolation is only performed in temperature space - for densities
c     outside the given range - the values for the high or low densities
c     are used as appropriate.
c
      fxte = (xte-rtlt(i1e))/(rtlt(i1e+1)-rtlt(i1e))
      fxti = (xti-rtlt(i1i))/(rtlt(i1i+1)-rtlt(i1i))
      fxni = (xni-rtln(i2))/(rtln(i2+1)-rtln(i2))
c
c ... For given za and zn, find the species index, ii, in the table.
c
      ii = -1
      do i=0,rtnsd-1
         if ((zn .eq. nint(rtzn(i))) .and. (za .eq. nint(rtza(i)))) then
            ii = i
            goto 100
         endif
      enddo
c

 100  if (ii .lt. 0) then
         write (*,*) '*** mcrates could not find za=',za,' zn=',zn
         write (*,*) '*** check mcfilenames array'
         stop
      endif
c
c     Compute rate parameters for transitions from table species ii.
c
      if (opt.eq.0) then
c
         if (za .lt. zamax) then
            rion = exp((1-fxte)*rtlsa(i1e,0,ii)+fxte*rtlsa(i1e+1,0,ii))
            if (za .eq. 0) return
         endif
c
         rrec = exp((1-fxte)*rtlra(i1e,0,ii)+fxte*rtlra(i1e+1,0,ii))
         rcxr = exp((1-fxti)*rtlcx(i1i,0,ii)+fxti*rtlcx(i1i+1,0,ii))
c
      elseif (opt.eq.1) then
c
         if (za .lt. zamax) then
            if (i2.eq.rtnn) then
c
               rion = exp((1-fxte)*rtlsa(i1e,rtnn,ii)
     >                       +fxte*rtlsa(i1e+1,rtnn,ii))
c
            elseif (i2.eq.-1) then
c
               rion = exp((1-fxte)*rtlsa(i1e,0,ii)
     >                       +fxte*rtlsa(i1e+1,0,ii))
c
            else
c
               rion = exp(
     >            (1.0-fxni) * ((1-fxte)*rtlsa(i1e,i2,ii)
     >                          +fxte*rtlsa(i1e+1,i2,ii))
     >            + fxni     * ((1-fxte)*rtlsa(i1e,i2+1,ii)
     >                          +fxte*rtlsa(i1e+1,i2+1,ii)))
c
            endif
c
            if (za .eq. 0) return
c
         endif
c
c        Recombination and CX
c
         if (i2.eq.rtnn) then
c
            rrec = exp((1-fxte)*rtlra(i1e,rtnn,ii)
     >                    +fxte*rtlra(i1e+1,rtnn,ii))
            rcxr = exp((1-fxte)*rtlcx(i1e,rtnn,ii)
     >                    +fxte*rtlcx(i1e+1,rtnn,ii))
c
         elseif (i2.eq.-1) then
c
            rrec = exp((1-fxte)*rtlra(i1e,0,ii)
     >                    +fxte*rtlra(i1e+1,0,ii))
            rcxr = exp((1-fxte)*rtlcx(i1e,0,ii)
     >                    +fxte*rtlcx(i1e+1,0,ii))
c
         else
c
            rrec = exp(
     >         (1.0-fxni) * ((1-fxte)*rtlra(i1e,i2,ii)
     >                       +fxte*rtlra(i1e+1,i2,ii))
     >         + fxni     * ((1-fxte)*rtlra(i1e,i2+1,ii)
     >                       +fxte*rtlra(i1e+1,i2+1,ii)))
c
            rcxr = exp(
     >         (1.0-fxni) * ((1-fxte)*rtlcx(i1e,i2,ii)
     >                       +fxte*rtlcx(i1e+1,i2,ii))
     >         + fxni     * ((1-fxte)*rtlcx(i1e,i2+1,ii)
     >                       +fxte*rtlcx(i1e+1,i2+1,ii)))
c
         endif
c
      endif
c
      return
      end
c
c
c
      subroutine readmc
      use mod_params
      use mod_adpak_com
      implicit none
c
c     include 'params'
c     include 'adpak_com'
c
c     READMC - this routine and READMC1 are responsible for reading the
c              ADPAK type data into the common blocks in the ADPAK common
c              block.
c
c
c     Use(Input)            # nzdf,mcfilename
c     Use(Multicharge)
c
c     local variables --
      integer i, ios, kstart, n, nget, rtnt_old, rtnn_old, rtnsd_old
      character idcod*8, idtyp*8, id1*32
c
c     For now - assign unit 73 to the ADPAK data file.
c
      parameter (nget=73)
c
      character*120 mcfilenames(5)
c
      integer nzdf
      parameter(nzdf=1)
c
c     procedures --
c      external freeus, gchange, kaboom
c
c----------------------------------------------------------------------c
c     Read rate data from 'un*formatted' files.
c           (file format from b2.5 code by B. Braams)
c----------------------------------------------------------------------c
c
c     Set input atomic data file name from selector
c
      mcfilenames(1)= mcfile
c
      rtnt=0
      rtnn=0
      rtnsd=0
c
      do i=1,nzdf
         rtnt_old  = rtnt
         rtnn_old  = rtnn
         rtnsd_old = rtnsd
c
c        call freeus(nget)
c
c        Open file
c
         open (nget, file=mcfilenames(i), form='formatted', iostat=ios,
     >        status='old')
c
c        check for errors
c
         if (ios .ne. 0) then
            write(*,*)
     >       '*** Input file mcfilename=',mcfilenames(i),' not found.'
     >       //' Index = ',i
            stop
c
         endif
c
c        read header --
c        un*formatted read for header data
c
         read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
         read (nget,'(1x,1a120)') labelrt(i)
c
c        read dimensions --
c        un*formatted read for integer data
c
         read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
         read (nget,*) rtnt,rtnn,rtns
c
         if (rtnt.gt.maxrtnt.or.rtnn.gt.maxrtnn.or.
     >      (rtns+rtnsd_old).gt.maxrtnsd) then
            write(6,*) 'SUBROUTINE READMC - Storage required'
     >                 //' for database exceeds allocated amount'
            write(6,*) 'Sizes:  rtnt =',rtnt,' maxrtnt =',maxrtnt
            write(6,*) '        rtnn =',rtnt,' maxrtnn =',maxrtnn
            write(6,*) '        rtns =',rtns,' maxrtnsd=',maxrtnsd
            write(6,*) '        rtnsd=',rtnsd,' maxrtnsd=',maxrtnsd
c
            write(7,*) 'SUBROUTINE READMC - Storage required'
     >                 //' for database exceeds allocated amount'
            write(7,*) 'Sizes:  rtnt =',rtnt,' maxrtnt =',maxrtnt
            write(7,*) '        rtnn =',rtnt,' maxrtnn =',maxrtnn
            write(7,*) '        rtns =',rtns,' maxrtnsd=',maxrtnsd
            write(7,*) '        rtnsd=',rtnsd,' maxrtnsd=',maxrtnsd
c
            write(0,*) 'SUBROUTINE READMC - Storage required'
     >                 //' for database exceeds allocated amount'
            write(0,*) 'Sizes:  rtnt =',rtnt,' maxrtnt =',maxrtnt
            write(0,*) '        rtnn =',rtnt,' maxrtnn =',maxrtnn
            write(0,*) '        rtns =',rtns,' maxrtnsd=',maxrtnsd
            write(0,*) '        rtnsd=',rtnsd,' maxrtnsd=',maxrtnsd
c
            stop
c
         endif
c
c        Test for compatibility of (rtnt,rtnn) from different tables:
c
         if ( (i .gt. 1) .and.
     >       ((rtnt .ne. rtnt_old) .or. (rtnn .ne. rtnn_old)) ) then
            write(*,*)
     >       '*** subroutine readmc: incompatible table dimensions in',
     >       mcfilenames(i),' and ',mcfilenames(i-1)
            stop
c
         endif
c
c
c        allocate storage --
c         call gchange("Multicharge",0)
c
c
         rtnsd=rtnsd+rtns
c
c        read abscissae and rates --
c
c        kstart = starting species index for this table
c
         kstart=rtnsd_old
c
         call readmc1 (nget, kstart)
c
         close (nget)

      enddo
c
      write (6,*) 'readmc:',rtnt,rtnn,rtnsd,kstart
      write (6,*) (rtza(i),i=0,rtnsd-1)
      write (6,*) (rtzn(i),i=0,rtnsd-1)
      write (6,*) (rtza2(i),i=0,rtnsd-1)
c
      return
      end

c-----------------------------------------------------------------------

      subroutine readmc1 (nget, kstart)
      use mod_params
      use mod_adpak_com
      implicit none
      integer nget, kstart
c
c     include 'params'
c     include 'adpak_com'
c
c     Use(Multicharge)
c
c     local variables --
c
      integer i, j, k, n
      character idcod*8, idtyp*8, id1*32
c
c     read abscissae --
c     un*formatted read for real data
c
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtza(k),k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtzn(k),k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtza2(k),k=kstart,kstart+rtns-1)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtt(i),i=0,rtnt)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtn(j),j=0,rtnn)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtlt(i),i=0,rtnt)

      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (rtln(j),j=0,rtnn)
c
c     read rate coefficients --
c     un*formatted read for real data
c
c
c     Ionization
c
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlsa(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)
c
c     Recombination
c
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlra(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)
c
c     Electron energy transfer
c
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlqa(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)
c
c     CX recombination
c
      read (nget,'(2a8,i12,4x,a32)') idcod, idtyp, n, id1
      read (nget,*) (((rtlcx(i,j,k),i=0,rtnt),j=0,rtnn),
     .                              k=kstart,kstart+rtns-1)
c
      return
      end
c
c
c
      subroutine inelinput(cion)
c
c     This subroutine reads INEL impurity rate information
c
      use mod_params
      use mod_inel
      use mod_adpak_com
      implicit none
      integer cion
c
c     include 'params'
c     include 'inel'
c     include 'adpak_com'
c
c ... Local variables:
c
      integer us, ios
c
c     Use unit number 73 since it is also used for the adpak data.
c
      parameter(us=73)
      integer i, j, k
c
      nz = cion
c
c     read in multi charge state data (e.g., file 'carbmc.dat')
c     File name is in variable mcfile
c
      open (us, file=mcfile, form='formatted', iostat=ios,
     .      status='old')
      if (ios .ne. 0) then
         write (*,*) '*** Input file inelmc=', mcfile, ' not found'
         stop
      endif
c
c     ionization rate
c
      do 200 i = 1,ntev
          read(us,1000) tevb(i),(rsi(i,k),k=0,nz-1)
c
c          do 201 k = 0,nz-1
c 201         rsi(i,k) = rsi(i,k)
c
 200  continue
c
c     recombination rate
c
      do 210 i = 1,ntev
         read(us,1000) tevb(i),(rre(i,k),k=1,nz)
c
c         do 211 k = 1,nz
c 211        rre(i,k) = rre(i,k)
c
 210  continue
c
c     radiative power rate
c
      do 220 i = 1,ntev
         read(us,1001) tevb(i),(rpwr(i,k),k=0,nz)
         do 221 k = 0,nz
            rpwr(i,k) = rpwr(i,k)*1.602e-19/1.5
 221     continue
 220  continue
c
c   CX recomb. rate
c
      do 230 i = 1,ntev
         read(us,1000) tevb(i),(rrcx(i,k),k=1,nz)
c
c         do 231 k = 1,nz
c 231        rrcx(i,k) = rrcx(i,k)
c
 230  continue
c
c ... Scale temperature scale for multi-charge rates to code units.
c
c     DIVIMP uses Te in eV - so leave as is
c
c      do i = 1, ntev
c         tevb(i)=tevb(i)*1.602e-19
c      enddo
c
      close (us)

 1000 format(7(1pe12.4))
 1001 format(8(1pe12.4))

      return
      end
c
c
c
      subroutine imprates(temp,kk,nzarg,rioniz,rrecomb,rcxrecom)
      use mod_params
      use mod_inel
      implicit none
c
      real temp
      integer kk, nzarg
      real rioniz, rrecomb, rcxrecom
c
c ... Given temperature "temp" and a charge state "kk", which is less
c     than or equal to the highest state "nzarg", interpolate from
c     tabulated rates for a particular impurity of ionization,
c     recombination, and charge-exchange recombination.
c     Note:  no scaling of temperatures is done here, so temp and tevb
c     must be provided in the same set of units.
c
c     include 'params'
c     include 'inel'
c
c     Local variables
c
      integer itemp
      real xltemn, dlogte
c
      rioniz   = 0.0
      rrecomb  = 0.0
      rcxrecom = 0.0
c
      xltemn = log10(tevb(1))
      dlogte = log10(tevb(2)) - xltemn
c
c ... Find index itemp into table such that
c        tevb(itemp) .le. temp .lt. tevb(itemp+1)
c
      itemp = int( ( log10( temp ) - xltemn ) / dlogte + 1. )
c
c ... For temperatures below minimum table temperature, extrapolate
c     downwards from table entries 1 and 2.

      itemp = max(1, itemp)
c
c ... For temperatures above maximum table temperature, extrapolate
c     upwards from table entries ntev-1 and ntev.

      itemp = min(ntev-1, itemp)
c
       if(kk .lt. nzarg)then
          rioniz = rsi(itemp,kk) + ( temp - tevb(itemp) )
     .   * ( rsi(itemp+1,kk) - rsi(itemp,kk) )
     .   / ( tevb(itemp+1) - tevb(itemp) )
          if(kk .eq. 0) return
       else
          rioniz = 0.0
       endif
c
          rrecomb = rre(itemp,kk) + ( temp - tevb(itemp) )
     .   * ( rre(itemp+1,kk) - rre(itemp,kk) )
     .   / ( tevb(itemp+1) - tevb(itemp) )
c
          rcxrecom = rrcx(itemp,kk) + ( temp - tevb(itemp) )
     .   * ( rrcx(itemp+1,kk) - rrcx(itemp,kk) )
     .   / ( tevb(itemp+1) - tevb(itemp) )

      return
      end


