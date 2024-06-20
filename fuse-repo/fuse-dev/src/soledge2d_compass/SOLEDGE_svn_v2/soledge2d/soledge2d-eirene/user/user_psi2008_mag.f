
*DK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ***********
C           *  TORE   *
C           ***********
C
      SUBROUTINE EIRENE_PROUSR (PRO,INDEX,P0,P1,P2,P3,P4,P5,PROVAC,N)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CCONA

      USE EIRMOD_CVARUSR_MAG

      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: PRO(*)
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      INTEGER, INTENT(IN) :: INDEX, N

      REAL(DP) :: ANGL, ANGLB
      INTEGER :: IR, IP, IT, j, j_mat
      integer :: ir0, ip0, it0
      integer :: ir1, ir2, ip1, ip2, it1, it2

c  first uniform profile with P0, then laser plume with P2

C
c  te profile
      if (index.eq.0) then

        do ir=1,nr1st
        	do ip=1,np2nd-1
        		do it=1,nt3rd-1
          			j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          			pro(j)=p0
        		enddo
        	enddo
        enddo
c      now loop on plasma plume
c      position of the source
        ir0 = 45
        ip0 = 151
        it0 = 32
c      size of the plasma plume
        ir1=ir0
        ir2=ir0+3
        ip1=ip0-1
        ip2=ip0+1
        it1=it0-1
        it2=it0+1
                 
        do ir=ir1,ir2
        	do ip=ip1,ip2
        		do it=it1,it2
          			j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          			pro(j)=p1
        		enddo
        	enddo
        enddo

c  ti profile, all species
      elseif (index.eq.1) then
        do ir=1,nr1st
        	do ip=1,np2nd-1
        		do it=1,nt3rd-1
          			j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          			pro(j)=p0
        		enddo
        	enddo
        enddo
c      now loop on plasma plume
c      position of the source
        ir0 = 45
        ip0 = 151
        it0 = 32
c      size of the plasma plume
        ir1=ir0
        ir2=ir0+3
        ip1=ip0-1
        ip2=ip0+1
        it1=it0-1
        it2=it0+1
                 
        do ir=ir1,ir2
        	do ip=ip1,ip2
        		do it=it1,it2
          			j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          			pro(j)=p1
        		enddo
        	enddo
        enddo
c  ni profile, all species
      elseif (index.eq.1+npls) then
        do ir=1,nr1st
        	do ip=1,np2nd-1
        		do it=1,nt3rd-1
          			j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          			pro(j)=p0
        		enddo
        	enddo
        enddo
c      now loop on plasma plume
c      position of the source
        ir0 = 45
        ip0 = 151
        it0 = 32
c      size of the plasma plume
        ir1=ir0
        ir2=ir0+3
        ip1=ip0-1
        ip2=ip0+1
        it1=it0-1
        it2=it0+1
                 
        do ir=ir1,ir2
        	do ip=ip1,ip2
        		do it=it1,it2
          			j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          			pro(j)=p1
        		enddo
        	enddo
        enddo
      elseif (index.eq.1+5*npls) then
c  B_x profile

c       first read the whole file and store it for subsequent calls

      call alloc_cvarusr(3)

      open(unit=324,file='B_field_centers_RZPhi_42920_t6',
     .   access='SEQUENTIAL',form='FORMATTED')

      icel=1
      do 
         read (324,*,END=1) icel, BTSx(icel), BTSy(icel), BTSz(icel),
     .   BTSf(icel)
         celnumB(icel)=icel
      end do

1     close(324)
     
      do ir=1,nr1st-1
        do ip=1,np2nd-1
        	do it=1,nt3rd-1
          		j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
c temporary correction: labeling problem in matlab routine 
c      			j_mat=IR+((IP-1)+(IT-1)*(NP2T3-1))*(NR1P2-1)
c                       pro(j)=BTSx(j_mat)
                        pro(j)=BTSx(j)
        		if (j.ne.celnumB(j)) then
      				write(6,*) 'inconsistency in cell numbers, Bx
     . , j= ',j
      				call eirene_exit_own(1)
      			end if  
                end do
        end do
      end do

      elseif (index.eq.2+5*npls) then
c  B_y profile
      do ir=1,nr1st-1
        do ip=1,np2nd-1
        	do it=1,nt3rd-1
          		j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
c      			j_mat=IR+((IP-1)+(IT-1)*(NP2T3-1))*(NR1P2-1)
c                        pro(j)=BTSy(j_mat)
      			pro(j)=BTSy(j)
        		if (j.ne.celnumB(j)) then
      				write(6,*) 'inconsistency in cell numbers, By'
      				call eirene_exit_own(1)
      			end if  
                end do
        end do
      end do

      elseif (index.eq.3+5*npls) then
c  B_z profile
      do ir=1,nr1st-1
        do ip=1,np2nd-1
        	do it=1,nt3rd-1
          		j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
c      			j_mat=IR+((IP-1)+(IT-1)*(NP2T3-1))*(NR1P2-1)
c                       pro(j)=BTSz(j_mat)
      			pro(j)=BTSz(j)
        		if (j.ne.celnumB(j)) then
      				write(6,*) 'inconsistency in cell numbers, Bz'
      				call eirene_exit_own(1)
      			end if  
                end do
        end do
      end do

      elseif (index.eq.4+5*npls) then
c  B_f profile (field strength)
      do ir=1,nr1st-1
        do ip=1,np2nd-1
        	do it=1,nt3rd-1
          		j=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
c      			j_mat=IR+((IP-1)+(IT-1)*(NP2T3-1))*(NR1P2-1)
c                       pro(j)=BTSf(j_mat)
      			pro(j)=BTSf(j)
        		if (j.ne.celnumB(j)) then
      				write(6,*) 'inconsistency in cell numbers, Bf'
      				call eirene_exit_own(1)
      			end if  
                end do
        end do
      end do

      call dealloc_cvarusr_mag      

      endif
      RETURN
      END
C
C
      SUBROUTINE EIRENE_GEOUSR

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRPTL
      USE EIRMOD_COMUSR

      IMPLICIT NONE

      real(dp), save :: emin, emax, ecut, angmin, angmax, dele1, dele2,
     .                  angmin3, angmax3, dela3, dela
      integer, save :: ifirst=0
      integer :: nener, nener1, nener2, j, jj, i1, i2, ic, ispz, nang
      integer :: nvolpl
c  make poloidal polygon no. 1 absorbing for d+ ions
      ispz=natmi+nmoli+2
      recyct(ispz,nlim+3)=0.
C
C   PREPARE DATA FOR LIMITER-SURFACES
C
C
C MODIFY GEOMETRY
C
C
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2 SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C
C
C
C  SETZE EINIGE VOLUMINA EXPLIZIT
C
C  SETZE PLOT GITTER FUER ADDITIONAL SURFACE TALLIES
C

C     yannick ... check what it means (removed because of run time checks)
C     if (nvolpl == 0) return
 
      if (ifirst.eq.0) then
        ifirst=1
        emin=0.
        emax=1.6d3
        nener=nr1st
c  gridpoint nener1 is at ecut
        ecut=200.
        nener1=21
c
        if (nener1.le.1) then
          ecut=emin
          nener1=0
        endif
        if (ecut.le.emin) then
          ecut=emin
          nener1=0
        endif
c
        nener2=nener-nener1
        if (ecut.gt.emin) dele1=(ecut-emin)/real(nener1-1,dp)
        dele2=(emax-ecut)/real(nener2-1,dp)
c
        nang=np2nd
        angmin=-180.
        angmax=180.
        dela=(angmax-angmin)/real(nang-1,dp)
c
        angmin3=-135.
        angmax3=-45.
        dela3=(angmax3-angmin3)/real(nang-1,dp)
c
c  set grid for plot no. 3
        do j=1,nener1
          xxp2d_usr(j,3) = emin + (j-1)/real(nener1-1,dp)*(ecut-emin)
C         write (6,*) 'xxp2d_usr ',j,xxp2d_usr(j,3)
        enddo
        do j=nener1+1,nener
          jj=j-nener1
          xxp2d_usr(j,3) = ecut + (jj)/real(nener2,dp)*(emax-ecut)
C         write (6,*) 'xxp2d_usr ',j,xxp2d_usr(j,3)
        enddo
C
        do i1=1,nr1st
          do i2=2,np2nd
            ic=i1+(i2-1)*nr1st
            xxp2d_usr(ic,3)=xxp2d_usr(i1,3)
          enddo
        enddo
C  set grid for plot no. 5
        do i1=1,nr1st
          do i2=1,np2nd
            ic=i1+(i2-1)*nr1st
            xxp2d_usr(ic,5)=xxp2d_usr(i1,3)
          enddo
        enddo
      endif
      RETURN
      END
C
      SUBROUTINE EIRENE_PLAUSR
      IMPLICIT NONE
      RETURN
      END

c      SUBROUTINE PLTUSR(PLABLE,J)
c      IMPLICIT NONE
c      LOGICAL, INTENT(IN) :: PLABLE
c      INTEGER, INTENT(IN) :: J
c      RETURN
c      END
C
C
      SUBROUTINE EIRENE_SAMUSR (NLSF,XX0,YY0,ZZ0,
     .              SORD1,SORD2,SORD3,SORD4,SORD5,SORD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,WEISPZ)
C
C  SAMPLE INITAL COORDIANTES X,Y,Z ON ADDITIONAL SURFACE NLLI
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_CGRID
      USE EIRMOD_COMUSR
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_CESTIM
      USE EIRMOD_CVARUSR_MAG
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: SORD1,SORD2,SORD3,SORD4,SORD5,SORD6
      REAL(DP), INTENT(OUT) :: XX0,YY0,ZZ0,TEWL,TIWL(*),DIWL(*),
     .                       VXWL(*),VYWL(*),VZWL(*),WEISPZ(*),
     .                       EFWL(*),SHWL
      INTEGER, INTENT(IN) :: NLSF,is1, is2
      INTEGER, INTENT(OUT) :: IRUSR, IPUSR, ITUSR, IAUSR, IBUSR
      REAL(DP) :: X, Y, T, B0, B1, B2, Z1, Z2, zep1, rad, pphi, wink,
     .            aout, ain, aslice, frac, tot_area
      real(dp), allocatable :: dummy(:)
      REAL(DP), EXTERNAL :: RANF_EIRENE
      integer, external :: eirene_learca, eirene_learc1
      INTEGER :: IER, i, irc, itc, icell, ic, il, im, iu, nnt

      entry eirene_sm0usr (is1,is2,sord1,sord2,sord3,sord4,sord5,sord6)

      open (unit=32,file='radius.txt',access='SEQUENTIAL',
     .      form='FORMATTED')

      nrciel = 0
      do 
         read (32,*,END=1) i, radciel(i) 
         radciel(i) = radciel(i) * 100._dp
      end do

 1    continue

      close (unit=32)
      nrciel = i

      allocate (flxciel(nrciel-1,nttra-1))
      allocate (tciel(nrciel-1,nttra-1))
      allocate (flxYchem(nrciel-1,nttra-1))

      open (unit=32,file='flx.txt',access='SEQUENTIAL',
     .      form='FORMATTED')

      do irc = 1, nrciel-1
         read (32,*,END=2) (flxciel(irc,itc),itc=1,nttra-1) 
      end do

 2    if (irc < nrciel-1) then
        write (iunout,*) ' premature end of data in file flx.txt'
        write (iunout,*) ' calculation abandonned '
        call eirene_exit_own(1)
      end if

      close (unit=32)

      open (unit=32,file='temp.txt',access='SEQUENTIAL',
     .      form='FORMATTED')

      do irc = 1, nrciel-1
         read (32,*,END=3) (tciel(irc,itc),itc=1,nttra-1) 
      end do

 3    if (irc < nrciel-1) then
        write (iunout,*) ' premature end of data in file temp.txt'
        write (iunout,*) ' calculation abandonned '
!yannick        call EIRENE_eirene_exit_own(1)
        call eirene_exit_own(1)
      end if

      close (unit=32)

c reads flux file for chemical erosion (includes CX contribution)

      open(unit=32,file='flxYchem.txt',access='SEQUENTIAL',
     .     form='FORMATTED')

      do irc=1, nrciel-1
          read (32,*,END=4) (flxYchem(irc,itc),itc=1,nttra-1)
      end do

 4    if  (irc < nrciel-1) then
        write (iunout,*) 'premature end of data in file flxYchem.txt '
        write (iunout,*) 'calculation abondonned'
        call eirene_exit_own(1)
      end if

      close(unit=32) 

! convert temperature distribution on CIEL from Kelvin to eV
! change sign of temperature for technical reasons
! Tciel > 0   ==>  Tciel = 3/2 k Twall => twall = 2/3*Tciel
! Tciel < 0   ==>  twall = -Tciel
      tciel = -tciel * evkel

      do irc = 1, nrciel-1
        drad(irc) = radciel(irc+1) - radciel(irc)
      end do

      allocate (delph(nttra-1))
      do itc = 1,nttra-1
        delph(itc) = zsurf(itc+1) - zsurf(itc)
      end do

      nbincll = nrciel * nttra
      
      if (nbincll > nrad) then
        write (iunout,*) ' NUMBER OF CELLS IN GEOMETRY TOO SMALL '
        write (iunout,*) ' NO BINNING OF RESULTS PERFORMED '
        write (iunout,*) ' NRAD =   ',NRAD
        write (iunout,*) ' NBINCLL =', NBINCLL
      end if

      allocate (area_ciel(nrciel-1,nttra-1))
      allocate (area1d(nbincll))

      do itc = 1,nttra-1

        frac = delph(itc) / pi2a

        do irc = 1, nrciel-1

          aout = radciel(irc+1)**2 * PIA * 1.E-4_DP
          ain = radciel(irc)**2 * PIA * 1.E-4_DP
          aslice = aout - ain       
! area in m**2
          area_ciel(irc,itc) = aslice * frac
          icell = irc + (itc-1)*nrciel
          area1d(icell) = area_ciel(irc,itc)

        end do
      end do

      allocate (dummy(nbincll))
      dummy = 1._dp
      call eirene_inttal (area1d,dummy,1,1,nbincll,
     .             tot_area,nrciel,nttra,1,1)
      deallocate (dummy)

      nflxcl = count(flxciel > 0._dp)

      allocate (vfunc(0:nflxcl))
      allocate (ifunc(nflxcl))

      vfunc = 0._dp
      ic = 0

      do itc = 1, nttra-1
         do irc = 1, nrciel-1
            if (flxciel(irc,itc) > 0._dp) then
               icell = irc + (itc-1)*nrciel
               ic = ic + 1
               ifunc(ic) = icell
               vfunc(ic) = vfunc(ic-1) + 
     .                     flxciel(irc,itc)*area_ciel(irc,itc)*elcha
            end if
         end do
      end do
      
      flux(is2) = vfunc(nflxcl)
      
      return

      entry eirene_SM1USR (NLSF,XX0,YY0,ZZ0,
     .              SORD1,SORD2,SORD3,SORD4,SORD5,SORD6,
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)

      zep1 = ranf_eirene()*vfunc(nflxcl)

      IL=0
      IU=nflxcl

c  binary search
      DO WHILE (IU-IL.gt.1)
        IM=(IU+IL)*0.5
        IF(ZEP1.GE.VFUNC(IM)) THEN
          IL=IM
        ELSE
          IU=IM
        ENDIF
      END DO
c
      ICELL=ifunc(IU)

      itc = icell / nrciel + 1
      irc = mod(icell,nrciel)

! cell on ciel limiter has been identified
! now find coordiantes in cell

      rad = radciel(irc) + ranf_eirene()*drad(irc)
      pphi = zsurf(itc) + ranf_eirene()*delph(itc)

!pb      xx0 = rad * cos(pphi) - rmtor
      xx0 = rad - rmtor
      yy0 = p1(2,3)
!pb      call fzrtra (xx0,zz0,pphi,nnt)
      zz0 = pphi * raddeg
      
      irusr = eirene_learc1(XX0,YY0,ZZ0,IPOLG,1,NR1STM,
     .    .FALSE.,.FALSE.,NPANU,'SM1USR      ')
      WINK=MOD(ATAN2(YY0,XX0)+PI2A-PSURF(1),PI2A)+PSURF(1)
      ipusr=eirene_learca(WINK,PSURF,1,NP2ND,1,'SM1USR')
      itusr = 1

      iausr = 0
      ibusr = 1
      NBLCKA=NSTRD*(ibusr-1)+iausr
      NCELL=irusr+((ipusr-1)+(itusr-1)*NP2T3)*NR1P2+NBLCKA

      TEWL = tein(ncell)
      TIWL(1:nplsti) = tiin(1:nplsti,ncell)
      DIWL(1:nplsi) = diin(1:nplsi,ncell)
      VXWL(1:nplsv) = vxin(1:nplsv,ncell)
      VYWL(1:nplsv) = vyin(1:nplsv,ncell)
      VZWL(1:nplsv) = vzin(1:nplsv,ncell)
      EFWL(1:nplsi) = 0._dp
      SHWL = 0._dp
      WEISPZ(1:nspz) = 1._dp

      RETURN
      END
C
C
C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
C
CC
C  WV=WEIGHT/VEL
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS

      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG

      RETURN
      END
C
C
      SUBROUTINE EIRENE_UPSUSR(WPR,IND)
C  USER SUPPLIED ESTIMATOR, SURFACE AVERAGED
C  (COLLISION - AND TRACKLENGTH ESTIMATORS ARE IDENTICAL, IF SURFACE
C   AVERAGED)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CSPEZ
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRPTL
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_CVARUSR_MAG
      USE EIRMOD_CLGIN

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: WPR
      INTEGER, INTENT(IN) :: IND

      integer :: itc, irc, il, im, iu, icell,i
      real(dp) :: rad, aout, ain, aslice, frac

! to be done
! setze twall der additional surface auf lokalen wert aus ortsauflösung


c  roof of CIEL is hit
      if (masurf.ne.3) return

      if (.not.allocated(delph)) then 
      	
c  do the intialization which was done in sm0usr
      	open (unit=32,file='radius.txt',access='SEQUENTIAL',
     .      form='FORMATTED')

      	nrciel = 0
      	do 
         read (32,*,END=1) i, radciel(i) 
         radciel(i) = radciel(i) * 100._dp
      	end do

 1    	continue

      	close (unit=32)
      	nrciel = i

      	do irc = 1, nrciel-1
        	drad(irc) = radciel(irc+1) - radciel(irc)
      	end do

      	allocate (delph(nttra-1))
      	do itc = 1,nttra-1
        	delph(itc) = zsurf(itc+1) - zsurf(itc)
      	end do

      	nbincll = nrciel * nttra

      	allocate (area_ciel(nrciel-1,nttra-1))
      	allocate (area1d(nbincll))

      	do itc = 1,nttra-1

        	frac = delph(itc) / pi2a

        		do irc = 1, nrciel-1

          			aout = radciel(irc+1)**2 * PIA * 1.E-4_DP
          			ain = radciel(irc)**2 * PIA * 1.E-4_DP
          			aslice = aout - ain       
! area in m**2
          			area_ciel(irc,itc) = aslice * frac
          			icell = irc + (itc-1)*nrciel
          			area1d(icell) = area_ciel(irc,itc)

        		end do
      	end do

c  this needs to be done after nbincll is defined
        

      end if

      if (.not.allocated(pfl_in_at)) call alloc_cvarusr(1)

      itc = iperid
      rad = x0 + rmtor

      IL=0
      IU=nrciel

c  binary search
      DO WHILE (IU-IL.gt.1)
        IM=(IU+IL)*0.5
        IF(rad.GE.radciel(IM)) THEN
          IL=IM
        ELSE
          IU=IM
        ENDIF
      END DO
      
      irc = iu
      icell = irc + (itc-1)*nrciel

!===================================================================

! UGLY!!!!!! One shouldn't do that !!!!!!!!!
! Only for Yannicks PSI paper!!!!!!!!

c set wall temperature
c      ewall(masurf) = tciel(irc,itc)
c set flux for chemical erosion (conversion to cm2 to compensate 
c reconversion in sputer.f)
c      FLXOUT(masurf) = flxYchem(irc,itc)*1.D-4
      
!===================================================================


c  ingoing particles only
      if (ind == 1) then

        if (ityp == 1) then
          pfl_in_at(0,icell) = pfl_in_at(0,icell) + wpr
          pfl_in_at(iatm,icell) = pfl_in_at(iatm,icell) + wpr
          enfl_in_at(0,icell) = enfl_in_at(0,icell) + wpr*e0
          enfl_in_at(iatm,icell) = enfl_in_at(iatm,icell) + wpr*e0
        elseif (ityp == 2) then
          pfl_in_ml(0,icell) = pfl_in_ml(0,icell) + wpr
          pfl_in_ml(imol,icell) = pfl_in_ml(imol,icell) + wpr
          enfl_in_ml(0,icell) = enfl_in_ml(0,icell) + wpr*e0
          enfl_in_ml(imol,icell) = enfl_in_ml(imol,icell) + wpr*e0
        elseif (ityp == 3) then
          pfl_in_io(0,icell) = pfl_in_io(0,icell) + wpr
          pfl_in_io(iion,icell) = pfl_in_io(iion,icell) + wpr
          enfl_in_io(0,icell) = enfl_in_io(0,icell) + wpr*e0
          enfl_in_io(iion,icell) = enfl_in_io(iion,icell) + wpr*e0
        elseif (ityp == 4) then
          pfl_in_pl(0,icell) = pfl_in_pl(0,icell) - wpr
          pfl_in_pl(ipls,icell) = pfl_in_pl(ipls,icell) - wpr
          enfl_in_pl(0,icell) = enfl_in_pl(0,icell) - wpr*e0
          enfl_in_pl(ipls,icell) = enfl_in_pl(ipls,icell) - wpr*e0
        end if

      else

        if (ityp == 1) then
          pfl_out_at(0,icell) = pfl_out_at(0,icell) + wpr
          pfl_out_at(iatm,icell) = pfl_out_at(iatm,icell) + wpr
          enfl_out_at(0,icell) = enfl_out_at(0,icell) + wpr*e0
          enfl_out_at(iatm,icell) = enfl_out_at(iatm,icell) + wpr*e0
        elseif (ityp == 2) then
          pfl_out_ml(0,icell) = pfl_out_ml(0,icell) + wpr
          pfl_out_ml(imol,icell) = pfl_out_ml(imol,icell) + wpr
          enfl_out_ml(0,icell) = enfl_out_ml(0,icell) + wpr*e0
          enfl_out_ml(imol,icell) = enfl_out_ml(imol,icell) + wpr*e0
        elseif (ityp == 3) then
          pfl_out_io(0,icell) = pfl_out_io(0,icell) + wpr
          pfl_out_io(iion,icell) = pfl_out_io(iion,icell) + wpr
          enfl_out_io(0,icell) = enfl_out_io(0,icell) + wpr*e0
          enfl_out_io(iion,icell) = enfl_out_io(iion,icell) + wpr*e0
        elseif (ityp == 4) then
          pfl_out_pl(0,icell) = pfl_out_pl(0,icell) + wpr
          pfl_out_pl(ipls,icell) = pfl_out_pl(ipls,icell) + wpr
          enfl_out_pl(0,icell) = enfl_out_pl(0,icell) + wpr*e0
          enfl_out_pl(ipls,icell) = enfl_out_pl(ipls,icell) + wpr*e0
        end if

      end if
 
      RETURN
      END
C
C
      SUBROUTINE EIRENE_UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
      RETURN
      END
C
C
      SUBROUTINE EIRENE_CRVUSR (ILIMI,TA,TB,IND,P,Q)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: P(*), Q(*)
      REAL(DP), INTENT(IN) :: T, TA, TB
      REAL(DP), INTENT(OUT) :: X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ
      INTEGER, INTENT(IN) :: ILIMI, IND
      ENTRY PNTUSR(T,X,Y,Z,IND,P,Q)
        X=5.+T
        Y=20.
        Z=2.5
      RETURN
      ENTRY DPTUSR(T,DX,DY,DZ,IND,P,Q)
        DX=1.
        DY=0.
        DZ=0.
      RETURN
      ENTRY DDPUSR(T,DDX,DDY,DDZ,IND,P,Q)
        DDX=0.
        DDY=0.
        DDZ=0.
      RETURN
      END

      SUBROUTINE EIRENE_sigusr
      return
      end

      SUBROUTINE EIRENE_upnusr
      return
      end
c
c
      SUBROUTINE EIRENE_talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL1,TXTSP1,TXTUN1,ILAST,*)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CGRID
      USE EIRMOD_CCONA
      USE EIRMOD_COUTAU
      USE EIRMOD_CVARUSR_MAG

      implicit none

      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(out) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl1,txtsp1,txtun1
      real(dp) :: add, dummy(nbincll)
      real(dp), allocatable :: aver(:)
      integer :: i, ip, ir, icell, ndim
      character(24) :: txat(0:natm), txml(0:nmol), txio(0:nion),
     .                 txpl(0:npls), txtun
      character(72) :: txttl

      ilast = 0

      if (.not.allocated(pfl_in_at)) goto 4711

! scale additional tallies sampled on top of CIEL

      dummy(1:nbincll) = flxfac(istra) / (area1d(1:nbincll)*elcha)

      do iatm = 0, natmi
        pfl_in_at(iatm,1:nbincll) = pfl_in_at(iatm,1:nbincll)*
     .                              dummy(1:nbincll)
        pfl_out_at(iatm,1:nbincll) = pfl_out_at(iatm,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_in_at(iatm,1:nbincll) = enfl_in_at(iatm,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_out_at(iatm,1:nbincll) = enfl_out_at(iatm,1:nbincll)*
     .                                dummy(1:nbincll)
      end do

      do imol = 0, nmoli
        pfl_in_ml(imol,1:nbincll) = pfl_in_ml(imol,1:nbincll)*
     .                              dummy(1:nbincll)
        pfl_out_ml(imol,1:nbincll) = pfl_out_ml(imol,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_in_ml(imol,1:nbincll) = enfl_in_ml(imol,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_out_ml(imol,1:nbincll) = enfl_out_ml(imol,1:nbincll)*
     .                                dummy(1:nbincll)
      end do

      do iion = 0, nioni
        pfl_in_io(iion,1:nbincll) = pfl_in_io(iion,1:nbincll)*
     .                              dummy(1:nbincll)
        pfl_out_io(iion,1:nbincll) = pfl_out_io(iion,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_in_io(iion,1:nbincll) = enfl_in_io(iion,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_out_io(iion,1:nbincll) = enfl_out_io(iion,1:nbincll)*
     .                                dummy(1:nbincll)
      end do

      do ipls = 0, nplsi
        pfl_in_pl(ipls,1:nbincll) = pfl_in_pl(ipls,1:nbincll)*
     .                              dummy(1:nbincll)
        pfl_out_pl(ipls,1:nbincll) = pfl_out_pl(ipls,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_in_pl(ipls,1:nbincll) = enfl_in_pl(ipls,1:nbincll)*
     .                               dummy(1:nbincll)
        enfl_out_pl(ipls,1:nbincll) = enfl_out_pl(ipls,1:nbincll)*
     .                                dummy(1:nbincll)
      end do

      ndim = max(natm, nmol, nion, npls)
      allocate (aver(0:ndim))

      txat(0) = 'total'
      txat(1:natmi) = texts(nsph+1:nspa)
      txml(0) = 'total'
      txml(1:nmoli) = texts(nspa+1:nspam)
      txio(0) = 'total'
      txio(1:nioni) = texts(nspam+1:nspami)
      txpl(0) = 'total'
      txpl(1:nplsi) = texts(nspami+1:nsptot)

      dummy = 1._dp

      call eirene_headng ('PARTICLE FLUXES ONTO CIEL',25)

      txttl='ATOMIC PARTICLE FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_in_at, dummy, 0, natm, natmi,
     .                tpfl_in_at, aver, txttl, txat, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='MOLECULAR PARTICLE FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_in_ml, dummy, 0, nmol, nmoli,
     .                tpfl_in_ml, aver, txttl, txml, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='TEST ION PARTICLE FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_in_io, dummy, 0, nion, nioni,
     .                tpfl_in_io, aver, txttl, txio, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='BULK ION PARTICLE FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_in_pl, dummy, 0, npls, nplsi,
     .                tpfl_in_pl, aver, txttl, txpl, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)


      call eirene_headng ('ENERGY FLUXES ONTO CIEL',23)

      txttl='ATOMIC ENERGY FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_in_at, dummy, 0, natm, natmi,
     .                tenfl_in_at, aver, txttl, txat, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='MOLECULAR ENERGY FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_in_ml, dummy, 0, nmol, nmoli,
     .                tenfl_in_ml, aver, txttl, txml, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='TEST ION ENERGY FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_in_io, dummy, 0, nion, nioni,
     .                tenfl_in_io, aver, txttl, txio, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='BULK ION ENERGY FLUXES ONTO TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_in_pl, dummy, 0, npls, nplsi,
     .                tenfl_in_pl, aver, txttl, txpl, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)
           

      call eirene_headng ('PARTICLE FLUXES FROM CIEL',25)

      txttl='ATOMIC PARTICLE FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_out_at, dummy, 0, natm, natmi,
     .                tpfl_out_at, aver, txttl, txat, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='MOLECULAR PARTICLE FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_out_ml, dummy, 0, nmol, nmoli,
     .                tpfl_out_ml, aver, txttl, txml, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='TEST ION PARTICLE FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_out_io, dummy, 0, nion, nioni,
     .                tpfl_out_io, aver, txttl, txio, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='BULK ION PARTICLE FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (pfl_out_pl, dummy, 0, npls, nplsi,
     .                tpfl_out_pl, aver, txttl, txpl, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)


      call eirene_headng ('ENERGY FLUXES FROM CIEL',23)

      txttl='ATOMIC ENERGY FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_out_at, dummy, 0, natm, natmi,
     .                tenfl_out_at, aver, txttl, txat, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='MOLECULAR ENERGY FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_out_ml, dummy, 0, nmol, nmoli,
     .                tenfl_out_ml, aver, txttl, txml, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='TEST ION ENERGY FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_out_io, dummy, 0, nion, nioni,
     .                tenfl_out_io, aver, txttl, txio, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)

      txttl='BULK ION ENERGY FLUXES FROM TOP OF LIMITER'
      txtun='#/s/m**2'
      call prt_array (enfl_out_pl, dummy, 0, npls, nplsi,
     .                tenfl_out_pl, aver, txttl, txpl, txtun,
     .                nrciel, 1, nttra, 1, nbincll, 2, 0)
           

! write raps outputs

c      call write_geo

c      call write_raps


 4711 continue

! cumulated source rate 
      do i=1,nrad
        vector(i) = 0
      enddo
      do ip = 1,np2nd
        add = 0
        do ir = 1,nr1stm
          icell=ir+(ip-1)*nr1st
          vector(icell)=add+(pael(icell)+pmel(icell)+piel(icell))*
     .                      vol(icell)
          add = vector(icell)
        enddo

      enddo
      taltot=vector(nr1stm+(np2nd-1)*nr1st)
      talav = 0._dp
      txttl1='cumulated source rate'
      txtsp1='electrons'
      txtun1='amp '
      return
      
      contains

      SUBROUTINE prt_array (ar, wei, lldim, ludim, nsp,
     .                      total, aver, txtl, txsp, txun,
     .                      nr, np, nt, nb, nc, ifl1, ifl2)

      implicit none
      integer, intent(in) :: lldim, ludim, nsp, nr, np, nt, nb,
     .                       nc, ifl1, ifl2
      real(dp), intent(inout) :: ar(lldim:ludim,nc)
      real(dp), intent(out) :: total(lldim:ludim), aver(lldim:ludim)
      real(dp), intent(in) :: wei(nc)
      character(len=*), intent(in) :: txtl,txsp(lldim:ludim),txun

      real(dp) :: vector(nc)
      integer :: isp

      aver = 0._dp

      do isp = lldim, nsp

        vector(1:nc) = ar(isp,1:nc)
        call eirene_inttal (vector,wei,1,1,nc,total(isp),nr,np,nt,1)     
        ar(isp,1:nc) = vector(1:nc)

        if (abs(total(isp)) > eps30) then 
          CALL EIRENE_PRTTAL(TXTL,TXSP(isp),TXUN,
     .                VECTOR,NR,NP,NT,1,NC,ifl1,ifl2)
          CALL eirene_leer(2)
          CALL eirene_masage
     .       ('TOTAL ("UNITS*CM**3), AND MEAN VALUE ("UNITS") ')
          CALL eirene_masr2 ('TOTAL, MEAN:    ',total(isp),aver(isp))       
          CALL eirene_leer(3)
        else
          CALL EIRENE_PRTTAL(TXTL,TXSP(isp),TXUN,
     .                VECTOR,NR,NP,NT,1,Nc,-1,0)
          CALL eirene_masage
     .           ('IDENTICAL ZERO, NOT PRINTED                  ')
          CALL eirene_leer(2)
        end if
      
      end do
      
      return
      end SUBROUTINE prt_array


      end SUBROUTINE EIRENE_talusr


      SUBROUTINE EIRENE_write_geo
      
      USE EIRMOD_precision
      USE EIRMOD_parmmod
      USE EIRMOD_cgrid
      USE EIRMOD_cvarusr_mag

      implicit none
      
      integer :: num, itc, irc, i1, i2, i3, i4
      real(dp) :: x, z
      
      open (unit=27,file='ciel.npco_char')

      WRITE(27,'(1X,A5,8X,A4,11X,A1,11X,A1,11X,A1,11X,A1)') '-1111',
     .      'NPCO','1','2','1','1'

      num = 0
      do itc = 1, nttra
        do irc = 1, nrciel
          num = num + 1
          x = radciel(irc)*cos(zsurf(itc))
          z = radciel(irc)*sin(zsurf(itc))
          WRITE(27,'(I6,1P,3ES12.4)') num, x, z         
        end do
      end do

      WRITE(27,'(1X,A5,8X,A3,12X,A1,11X,A1,11X,A1,11X,A1)') '-9999',
     .           'FIN','0','0','0','0'
     
      close (unit=27)

      
      open (unit=27,file='ciel.elemente')

      WRITE(27,'(A3,I3,A60,I6)') 'PSS',1,'BEISPIELDATEN',-3334
      WRITE(27,'(A8,I6,I9,I6)') 'QUAM4   ',1,(nrciel-1)*(nttra-1),4

      do itc = 1, nttra-1
        do irc = 1, nrciel-1
          i1 = irc + (itc-1)*nrciel
          i2 = irc+1 + (itc-1)*nrciel
          i3 = irc+1 + (itc+1-1)*nrciel
          i4 = irc + (itc+1-1)*nrciel
          WRITE(27,'(1X,A1,4I10)') '0',i1,i2,i3,i4
        end do
      end do
     
      close (unit=27)
      return
      
      end SUBROUTINE EIRENE_write_geo



      SUBROUTINE write_raps
      
      USE EIRMOD_precision
      USE EIRMOD_parmmod
      USE EIRMOD_cgrid
      USE EIRMOD_ccona
      USE EIRMOD_comusr
      USE EIRMOD_comprt
      USE EIRMOD_cvarusr_mag

      implicit none
      
      real(dp), allocatable :: flx_in(:,:), pl_flx(:)
      integer :: num, itc, irc, iraps, icell, i, ifirst,
     .           lz, l1, l2, irp
      character(6006) :: zout
      
      if (.not.allocated(cpfl_in_at)) call alloc_cvarusr(2)
      if (.not.allocated(flx_in)) then
        allocate (flx_in(1,nbincll))
        allocate (pl_flx(nbincll))
        flx_in = 0._dp
        pl_flx = 0._dp
      end if

      do itc = 1, nttra-1
        do irc = 1, nrciel-1
          icell = irc + (itc-1)*nrciel
          flx_in(1,icell) = flxciel(irc,itc)
        end do
      end do

      call to_corner (flx_in, pl_flx, vfunc(nflxcl), 1, 1, 1, nbincll)

      call to_corner (pfl_in_at, cpfl_in_at, tpfl_in_at, 
     .                0, natm, natmi, nbincll)
      call to_corner (pfl_in_ml, cpfl_in_ml, tpfl_in_ml, 
     .                0, nmol, nmoli, nbincll)
      call to_corner (pfl_in_io, cpfl_in_io, tpfl_in_io, 
     .                0, nion, nioni, nbincll)
      call to_corner (pfl_in_pl, cpfl_in_pl, tpfl_in_pl, 
     .                0, npls, nplsi, nbincll)

      call to_corner (enfl_in_at, cenfl_in_at, tenfl_in_at, 
     .                0, natm, natmi, nbincll)
      call to_corner (enfl_in_ml, cenfl_in_ml, tenfl_in_ml, 
     .                0, nmol, nmoli, nbincll)
      call to_corner (enfl_in_io, cenfl_in_io, tenfl_in_io, 
     .                0, nion, nioni, nbincll)
      call to_corner (enfl_in_pl, cenfl_in_pl, tenfl_in_pl, 
     .                0, npls, nplsi, nbincll)

      call to_corner (pfl_out_at, cpfl_out_at, tpfl_out_at, 
     .                0, natm, natmi, nbincll)
      call to_corner (pfl_out_ml, cpfl_out_ml, tpfl_out_ml, 
     .                0, nmol, nmoli, nbincll)
      call to_corner (pfl_out_io, cpfl_out_io, tpfl_out_io, 
     .                0, nion, nioni, nbincll)
      call to_corner (pfl_out_pl, cpfl_out_pl, tpfl_out_pl, 
     .                0, npls, nplsi, nbincll)

      call to_corner (enfl_out_at, cenfl_out_at, tenfl_out_at, 
     .                0, natm, natmi, nbincll)
      call to_corner (enfl_out_ml, cenfl_out_ml, tenfl_out_ml, 
     .                0, nmol, nmoli, nbincll)
      call to_corner (enfl_out_io, cenfl_out_io, tenfl_out_io, 
     .                0, nion, nioni, nbincll)
      call to_corner (enfl_out_pl, cenfl_out_pl, tenfl_out_pl, 
     .                0, npls, nplsi, nbincll)

      iraps = count(abs(tpfl_in_at(0:natmi)) > 0) +
     .        count(abs(tpfl_in_ml(0:nmoli)) > 0) +
     .        count(abs(tpfl_in_io(0:nioni)) > 0) +
     .        count(abs(tpfl_in_pl(0:nplsi)) > 0) +
     .        count(abs(tenfl_in_at(0:natmi)) > 0) +
     .        count(abs(tenfl_in_ml(0:nmoli)) > 0) +
     .        count(abs(tenfl_in_io(0:nioni)) > 0) +
     .        count(abs(tenfl_in_pl(0:nplsi)) > 0) +
     .        count(abs(tpfl_out_at(0:natmi)) > 0) +
     .        count(abs(tpfl_out_ml(0:nmoli)) > 0) +
     .        count(abs(tpfl_out_io(0:nioni)) > 0) +
     .        count(abs(tpfl_out_pl(0:nplsi)) > 0) +
     .        count(abs(tenfl_out_at(0:natmi)) > 0) +
     .        count(abs(tenfl_out_ml(0:nmoli)) > 0) +
     .        count(abs(tenfl_out_io(0:nioni)) > 0) +
     .        count(abs(tenfl_out_pl(0:nplsi)) > 0) + 1

      if (iraps > 500) then
         write (iunout,*) ' too many columns in raps output requested '
         write (iunout,*) ' only 500 columns written '
         iraps = 500
      end if

      open (unit=27,file='ciel.npst_char')
      open (unit=28,file='ciel.txt')

      WRITE(27,'(1X,A5,8X,A4,500(9X,I3))') '-1111',
     .        'NPST',1,IRAPS,1,(1,I=1,IRAPS)

      lz = iraps*12+6
      ifirst = 0

      num = 0

      tor:do itc = 1, nttra
        rad:do irc = 1, nrciel
          num = num + 1

          zout(1:lz) = repeat(' ',lz)
          l1 = 1
          l2 = 6
          write (zout(l1:l2),'(i6)') num
          l1 = l1 + 6
          l2 = l2 + 12
          if (ifirst == 0) then
             write (28,*) ' particle fluxes on CIEL (input) '
             write (28,*)
          end if
          write (zout(l1:l2),'(es12.4)') pl_flx(num)
          l1 = l1 + 12
          l2 = l2 + 12
          irp = 1

          do iatm=0,natmi
            if (abs(tpfl_in_at(iatm)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' atomic particle fluxes onto CIEL '
                write (28,*) ' iatm = ',iatm
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_in_at(iatm,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do imol=0,nmoli
            if (abs(tpfl_in_ml(imol)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' molecular particle fluxes onto CIEL '
                write (28,*) ' imol = ',imol
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_in_ml(imol,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do iion=0,nioni
            if (abs(tpfl_in_io(iion)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' test ion particle fluxes onto CIEL '
                write (28,*) ' iion = ',iion
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_in_io(iion,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do ipls=0,nplsi
            if (abs(tpfl_in_pl(ipls)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' bulk ion particle fluxes onto CIEL '
                write (28,*) ' ipls = ',ipls
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_in_pl(ipls,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do iatm=0,natmi
            if (abs(tpfl_out_at(iatm)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' atomic particle fluxes from CIEL '
                write (28,*) ' iatm = ',iatm
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_out_at(iatm,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do imol=0,nmoli
            if (abs(tpfl_out_ml(imol)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' molecular particle fluxes from CIEL '
                write (28,*) ' imol = ',imol
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_out_ml(imol,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do iion=0,nioni
            if (abs(tpfl_out_io(iion)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' test ion particle fluxes from CIEL '
                write (28,*) ' iion = ',iion
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_out_io(iion,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do ipls=0,nplsi
            if (abs(tpfl_out_pl(ipls)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' bulk ion particle fluxes from CIEL '
                write (28,*) ' ipls = ',ipls
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cpfl_out_pl(ipls,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do iatm=0,natmi
            if (abs(tenfl_in_at(iatm)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' atomic energy fluxes onto CIEL '
                write (28,*) ' iatm = ',iatm
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_in_at(iatm,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do imol=0,nmoli
            if (abs(tenfl_in_ml(imol)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' molecular energy fluxes onto CIEL '
                write (28,*) ' imol = ',imol
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_in_ml(imol,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do iion=0,nioni
            if (abs(tenfl_in_io(iion)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' test ion energy fluxes onto CIEL '
                write (28,*) ' iion = ',iion
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_in_io(iion,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do ipls=0,nplsi
            if (abs(tenfl_in_pl(ipls)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' bulk ion energy fluxes onto CIEL '
                write (28,*) ' ipls = ',ipls
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_in_pl(ipls,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do iatm=0,natmi
            if (abs(tenfl_out_at(iatm)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' atomic energy fluxes from CIEL '
                write (28,*) ' iatm = ',iatm
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_out_at(iatm,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do imol=0,nmoli
            if (abs(tenfl_out_ml(imol)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' molecular energy fluxes from CIEL '
                write (28,*) ' imol = ',imol
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_out_ml(imol,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do iion=0,nioni
            if (abs(tenfl_out_io(iion)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' test ion energy fluxes from CIEL '
                write (28,*) ' iion = ',iion
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_out_io(iion,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

          do ipls=0,nplsi
            if (abs(tenfl_out_pl(ipls)) > 0._dp) then
              if (ifirst == 0) then
                write (28,*) ' bulk ion energy fluxes from CIEL '
                write (28,*) ' ipls = ',ipls
                write (28,*)
              end if
              write (zout(l1:l2),'(es12.4)') cenfl_out_pl(ipls,num) 
              l1 = l1 + 12
              l2 = l2 + 12
              irp = irp + 1
              if (irp >= 500) goto 1000
            end if
          end do

 1000     continue
          WRITE(27,'(a6006)') zout 

          ifirst = 1
        end do rad
      end do tor

      WRITE(27,'(1X,A5,8X,A3,500(11X,I1))') '-9999',
     .           'FIN',(0,I=1,IRAPS)

      close (unit=27)
      close (unit=28)
 
      return
 
      contains

      SUBROUTINE to_corner (ar_in, ar_out, total, lldim, ludim, nsp, nc)

      implicit none
      integer, intent(in) :: lldim, ludim, nsp, nc
      real(dp), intent(in) :: ar_in(lldim:ludim,nc), total(lldim:ludim)
      real(dp), intent(out) :: ar_out(lldim:ludim,nc)
      real(dp), allocatable :: wgh(:), pl_flx(:)

      integer :: isp, itc, irc, i1, i2, i3, i4, icell
      real(dp) :: arinv

      allocate (wgh(nbincll))

      ar_out(lldim:ludim,1:nc) = 0._dp

      do isp = lldim, nsp

        if (abs(total(isp)) < eps30) cycle

        wgh = 0._dp

        do itc = 1, nttra-1
          do irc = 1, nrciel-1
            i1 = irc + (itc-1)*nrciel
            i2 = irc+1 + (itc-1)*nrciel
            i3 = irc+1 + (itc+1-1)*nrciel
            i4 = irc + (itc+1-1)*nrciel

            icell = irc + (itc-1)*nrciel
            arinv = 1._dp / area_ciel(irc,itc)

            wgh(i1) = wgh(i1) + arinv
            wgh(i2) = wgh(i2) + arinv
            wgh(i3) = wgh(i3) + arinv
            wgh(i4) = wgh(i4) + arinv

            ar_out(isp,i1) = ar_out(isp,i1) + ar_in(isp,icell)*arinv  
            ar_out(isp,i2) = ar_out(isp,i2) + ar_in(isp,icell)*arinv  
            ar_out(isp,i3) = ar_out(isp,i3) + ar_in(isp,icell)*arinv  
            ar_out(isp,i4) = ar_out(isp,i4) + ar_in(isp,icell)*arinv  
          end do
        end do

        ar_out(isp,1:nc) = ar_out(isp,1:nc) / wgh(1:nc)

      end do

      deallocate (wgh)
      
      return
      end SUBROUTINE to_corner

      end SUBROUTINE write_raps
c
c
c     SUBROUTINE EIRENE_diagno
c     return
c     end
c
c
      SUBROUTINE EIRENE_modusr
      return
      end
c
c
      SUBROUTINE EIRENE_retusr(sig)
      USE EIRMOD_PRECISION
      implicit none
      real(dp), intent(in) :: sig
      return
      end
C
C
      SUBROUTINE EIRENE_REFUSR
C
C   USER SUPPLIED REFLECTION MODEL
C
      USE EIRMOD_PRECISION
      use eirmod_clgin
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XMW,XCW,XMP,XCP,ZCOS,ZSIN,EXPI,RPROB,
     .                        E0TERM
      INTEGER, INTENT(IN) :: IGASF,IGAST

      real :: stick_C
      logical :: particle_mod

      ENTRY eirene_RF0USR
      ENTRY eirene_SPTUSR
      ENTRY eirene_SP0USR
 
! distinguish physical/self sputtering on C surfaces

      particle_mod=.true.

      if (particle_mod) then
  
      ISRS(3,3)=6
      ISRS(4,3)=6
      ISRS(5,3)=6
      ISRS(6,3)=6
      ISRS(10,3)=6
      ISRS(11,3)=6
      ISRS(13,3)=6
      ISRS(14,3)=6
      ISRS(15,3)=6

      ISRS(3,4)=6
      ISRS(4,4)=6
      ISRS(5,4)=6
      ISRS(6,4)=6
      ISRS(10,4)=6
      ISRS(11,4)=6
      ISRS(13,4)=6
      ISRS(14,4)=6
      ISRS(15,4)=6

      ISRS(3,5)=6
      ISRS(4,5)=6
      ISRS(5,5)=6
      ISRS(6,5)=6
      ISRS(10,5)=6
      ISRS(11,5)=6
      ISRS(13,5)=6
      ISRS(14,5)=6
      ISRS(15,5)=6

      ISRS(3,6)=6
      ISRS(4,6)=6
      ISRS(5,6)=6
      ISRS(6,6)=6
      ISRS(10,6)=6
      ISRS(11,6)=6
      ISRS(13,6)=6
      ISRS(14,6)=6
      ISRS(15,6)=6

! (fast) reflected particle species index -> C_refl for all C species

      ISRF(3,3)=5
      ISRF(4,3)=5
      ISRF(5,3)=5
      ISRF(6,3)=5
      ISRF(10,3)=5
      ISRF(11,3)=5
      ISRF(13,3)=5
      ISRF(14,3)=5
      ISRF(15,3)=5

      ISRF(3,4)=5
      ISRF(4,4)=5
      ISRF(5,4)=5
      ISRF(6,4)=5
      ISRF(10,4)=5
      ISRF(11,4)=5
      ISRF(13,4)=5
      ISRF(14,4)=5
      ISRF(15,4)=5

      ISRF(3,5)=5
      ISRF(4,5)=5
      ISRF(5,5)=5
      ISRF(6,5)=5
      ISRF(10,5)=5
      ISRF(11,5)=5
      ISRF(13,5)=5
      ISRF(14,5)=5
      ISRF(15,5)=5

      ISRF(3,6)=5
      ISRF(4,6)=5
      ISRF(5,6)=5
      ISRF(6,6)=5
      ISRF(10,6)=5
      ISRF(11,6)=5
      ISRF(13,6)=5
      ISRF(14,6)=5
      ISRF(15,6)=5

! (thermal) reflected particle species index -> C_refl for all C species

      ISRT(3,3)=5
      ISRT(4,3)=5
      ISRT(5,3)=5
      ISRT(6,3)=5
      ISRT(10,3)=5
      ISRT(11,3)=5
      ISRT(13,3)=5
      ISRT(14,3)=5
      ISRT(15,3)=5

      ISRT(3,4)=5
      ISRT(4,4)=5
      ISRT(5,4)=5
      ISRT(6,4)=5
      ISRT(10,4)=5
      ISRT(11,4)=5
      ISRT(13,4)=5
      ISRT(14,4)=5
      ISRT(15,4)=5

      ISRT(3,5)=5
      ISRT(4,5)=5
      ISRT(5,5)=5
      ISRT(6,5)=5
      ISRT(10,5)=5
      ISRT(11,5)=5
      ISRT(13,5)=5
      ISRT(14,5)=5
      ISRT(15,5)=5

      ISRT(3,6)=5
      ISRT(4,6)=5
      ISRT(5,6)=5
      ISRT(6,6)=5
      ISRT(10,6)=5
      ISRT(11,6)=5
      ISRT(13,6)=5
      ISRT(14,6)=5
      ISRT(15,6)=5

! sticking coefficient of C on C

      stick_C=0.5D0

      RECYCT(3,3)=1D0-stick_C
      RECYCT(4,3)=1D0-stick_C
      RECYCT(5,3)=1D0-stick_C
      RECYCT(6,3)=1D0-stick_C
      RECYCT(10,3)=1D0-stick_C
      RECYCT(11,3)=1D0-stick_C
      RECYCT(13,3)=1D0-stick_C
      RECYCT(14,3)=1D0-stick_C
      RECYCT(15,3)=1D0-stick_C

      RECYCT(3,4)=1D0-stick_C
      RECYCT(4,4)=1D0-stick_C
      RECYCT(5,4)=1D0-stick_C
      RECYCT(6,4)=1D0-stick_C
      RECYCT(10,4)=1D0-stick_C
      RECYCT(11,4)=1D0-stick_C
      RECYCT(13,4)=1D0-stick_C
      RECYCT(14,4)=1D0-stick_C
      RECYCT(15,4)=1D0-stick_C
 
      RECYCT(3,5)=1D0-stick_C
      RECYCT(4,5)=1D0-stick_C
      RECYCT(5,5)=1D0-stick_C
      RECYCT(6,5)=1D0-stick_C
      RECYCT(10,5)=1D0-stick_C
      RECYCT(11,5)=1D0-stick_C
      RECYCT(13,5)=1D0-stick_C
      RECYCT(14,5)=1D0-stick_C
      RECYCT(15,5)=1D0-stick_C

      RECYCT(3,6)=1D0-stick_C
      RECYCT(4,6)=1D0-stick_C
      RECYCT(5,6)=1D0-stick_C
      RECYCT(6,6)=1D0-stick_C
      RECYCT(10,6)=1D0-stick_C
      RECYCT(11,6)=1D0-stick_C
      RECYCT(13,6)=1D0-stick_C
      RECYCT(14,6)=1D0-stick_C
      RECYCT(15,6)=1D0-stick_C

      end if


      ENTRY eirene_SP1USR
      ENTRY eirene_RF1USR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,EXPI,
     .              RPROB,E0TERM,*,*,*,*)
      RETURN
      END
C
C
c      function leausr(a,b,c)
c      USE EIRMOD_PRECISION
c      IMPLICIT NONE
c      REAL(DP), INTENT(IN) :: A, B, C
c      INTEGER :: LEAUSR
c      leausr=1
c      RETURN
c      END


      SUBROUTINE EIRENE_TIMUSR(N,X,Y,Z,VX,VY,VZ,N1,N2,T,IC,IE,NP,NL)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: X,Y,Z,VX,VY,VZ,T,cx,cy,cz,sc
      INTEGER, INTENT(IN) :: N, N1, N2, IC, IE, NP, IS, NRCELL
      LOGICAL :: NL

      ENTRY eirene_NORUSR(is,x,y,z,cx,cy,cz,sc,VX,VY,VZ,NRCELL)

      RETURN
      END


      SUBROUTINE EIRENE_VOLUSR(N,A)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: A(*)
      INTEGER, INTENT(IN) :: N
      RETURN
      END


      SUBROUTINE EIRENE_VECUSR (I,VX,VY,VZ,IPLS)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(IN) :: VX,VY,VZ
      RETURN
      END


      SUBROUTINE EIRENE_INIusr
      IMPLICIT NONE
      RETURN
      END


c      SUBROUTINE TMSUSR (T0)
c      USE EIRMOD_PRECISION
c      IMPLICIT NONE
c      REAL(DP), INTENT(IN) :: T0
c      RETURN
c      END


c      function vdion(i)
c      USE EIRMOD_PRECISION
c      IMPLICIT NONE
c      INTEGER, INTENT(IN) :: I
c      REAL(DP) :: VDION
c      vdion=0.
c      RETURN
c      END
