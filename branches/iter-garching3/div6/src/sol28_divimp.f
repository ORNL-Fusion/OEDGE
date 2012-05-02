c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: divYield_Be_W
c
c      From an email from Sophie (16/11/2011), forwarded from Andreas Mutzke [mailto:Andreas.Mutzke@ipp.mpg.de] (16/11/2011)
c     
c       Be -> W 1 keV   e=  
c                2           1  ncp ncp_proj
c      SDTrimSP: VERSION 5.01    15.11.2011
c      Particles(NH*NR):        10000000 elements: Be   W      
c      CPT A-Z  A-MASS        DNS0         RHO    E_CUTOFF     E_DISPL     E_BULKB     E_SURFB       INEL0        ISBV    IINT
c       1   4   9.0122     0.12347     1.84770     2.00000    15.00000     0.00000     3.31000          3          1       2 
c       2  74 183.8400     0.06306    19.25000     2.00000    38.00000     0.00000     8.79000          3          1       2
c      !
c      !
c               20  number energy
c                1  number alpha
c                1  number mean depth
c                1  number refl.coeff  / energ.coeff
c                2  number sputt.coeff / energ.coeff
c            energy        alpha   mean depth    refl.coef  energ.coeff  sputt.coeff  energ.coeff  sputt.coeff  energ.coeff
c                                                    Be           Be           Be           Be           W            W    
c              1.00         0.00     0.753677     0.007354     0.000962     0.000000     0.000000     0.000000     0.000000
c              2.00         0.00     0.978780     0.133563     0.048357     0.000000     0.000000     0.000000     0.000000
c              3.00         0.00     1.248317     0.299761     0.131835     0.000000     0.000000     0.000000     0.000000
c              5.00         0.00     1.927753     0.544548     0.281648     0.000000     0.000000     0.000000     0.000000
c              7.00         0.00     2.684676     0.657730     0.370540     0.000000     0.000000     0.000000     0.000000
c             10.00         0.00     3.645109     0.709874     0.429249     0.000000     0.000000     0.000000     0.000000
c             20.00         0.00     5.587461     0.703014     0.450005     0.000000     0.000000     0.000000     0.000000
c             30.00         0.00     6.959517     0.676663     0.432105     0.000000     0.000000     0.000000     0.000000
c             50.00         0.00     9.133810     0.639375     0.399921     0.000000     0.000000     0.000000     0.000000
c             70.00         0.00    10.937870     0.615286     0.377816     0.000000     0.000000     0.000316     0.000004
c            100.00         0.00    13.249392     0.591254     0.355606     0.000000     0.000000     0.005947     0.000166
c            200.00         0.00    19.347540     0.551074     0.318757     0.000000     0.000000     0.045688     0.001636
c            300.00         0.00    24.276594     0.530139     0.300258     0.000000     0.000000     0.083110     0.002790
c            500.00         0.00    32.476271     0.505771     0.279909     0.000000     0.000000     0.135549     0.003844
c            700.00         0.00    39.537266     0.490152     0.267559     0.000000     0.000000     0.170108     0.004175
c           1000.00         0.00    48.928215     0.474044     0.254973     0.000000     0.000000     0.206025     0.004233
c           2000.00         0.00    75.322349     0.440764     0.230529     0.000000     0.000000     0.268974     0.003730
c           3000.00         0.00    98.145896     0.418819     0.215232     0.000000     0.000000     0.297168     0.003217
c           5000.00         0.00   139.313794     0.387903     0.194778     0.000000     0.000000     0.319288     0.002495
c           7000.00         0.00   177.686492     0.365954     0.181153     0.000000     0.000000     0.324028     0.002032
c
c
      REAL FUNCTION divYield_Be_W(e0)
      IMPLICIT none

      REAL, INTENT(IN) :: e0

      REAL yield(2,20),result
 
      DATA yield /
     .         1.00,  0.000000,
     .         2.00,  0.000000,
     .         3.00,  0.000000,
     .         5.00,  0.000000,
     .         7.00,  0.000000,
     .        10.00,  0.000000,
     .        20.00,  0.000000,
     .        30.00,  0.000000,
     .        50.00,  0.000000,
     .        70.00,  0.000316,
     .       100.00,  0.005947,
     .       200.00,  0.045688,
     .       300.00,  0.083110,
     .       500.00,  0.135549,
     .       700.00,  0.170108,
     .      1000.00,  0.206025,
     .      2000.00,  0.268974,
     .      3000.00,  0.297168,
     .      5000.00,  0.319288,
     .      7000.00,  0.324028  /

      divYield_Be_W = 0.0

      CALL Fitter(20,yield(1,1:20),yield(2,1:20),1,e0,result,'LINEAR')

      divYield_Be_W = result

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: divCompileSputteringYields
c
      SUBROUTINE divCompileSputteringYields
      USE mod_interface
      USE eckstein_2007_yield_data
      USE mod_divimp
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'cneut2'
c      INCLUDE 'pindata'
      INCLUDE 'slcom'


      REAL divYield_Be_W

      INTEGER       i,j,k,nseg,mion,fp,status,id,iz,ik,ir,matt,matp,
     .              ncore,in,nsputter,ierr,nymfs,nds2,n
      REAL          yield,rdum1,ratio,frac
      CHARACTER*512 fname,command,resdir
      CHARACTER     t1*3,t2*6

      nymfs = nymfs_global ! so lame, but need NYMFS and it's not in a common block
 
      nsputter = sputter_ndata

c
c     Set defaults for walls
c
      DO ID = 1, wallpts
        if (nymfs.gt.0.and.cymfs(1,1).eq.0) then
           if (wallpt(id,16).eq.1.or.wallpt(id,16).eq.4) then
              KMFPWS(ID) = CYMFS(1,3)
           else
              KMFPWS(ID) = CYMFS(1,6)
           endif
        else
           KMFPWS(ID) = 1.0
        endif
      end do
      do in=1,nymfs
         do id = NINT(cymfs(in,1)),NINT(cymfs(in,2))
            if (id.ge.1.and.id.le.wallpts) then
               if (wallpt(id,16).eq.1.or.wallpt(id,16).eq.4) then
                  KMFPWS(ID) = CYMFS(in,3)
               else
                  KMFPWS(ID) = CYMFS(in,6)
               endif
            endif
         enddo
      enddo

c
c     ------------------------------------------------------------------
c     LOAD PARTICLE FLUX DATA
c

      CALL GetEnv('RESDIRTOP',resdir)


      DO i = 1, nsputter

        IF (sputter_data(i)%data_type.EQ.4) THEN

          nds2 = 0
          DO id = 1, nds
            ir = irds(id)
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            nds2 = nds2 + 1
          ENDDO

          nseg = nds2
          mion = sputter_data(i)%atomic_number

          ALLOCATE(sputter_data(i)%ik   (nseg))
          ALLOCATE(sputter_data(i)%ir   (nseg))
          ALLOCATE(sputter_data(i)%id   (nseg))
          ALLOCATE(sputter_data(i)%r    (nseg))
          ALLOCATE(sputter_data(i)%z    (nseg))
          ALLOCATE(sputter_data(i)%dds  (nseg))
          ALLOCATE(sputter_data(i)%te   (nseg))
          ALLOCATE(sputter_data(i)%ti   (nseg))
          ALLOCATE(sputter_data(i)%flux (nseg,mion))
          ALLOCATE(sputter_data(i)%e0   (nseg,mion))
          ALLOCATE(sputter_data(i)%yield(nseg,mion))

          sputter_data(i)%type       = sputter_data(i)%data_type
          sputter_data(i)%nsegments  = nseg
          sputter_data(i)%charge_max = sputter_data(i)%atomic_number
	  sputter_data(i)%absfac     = 1.0
	  sputter_data(i)%flux       = 0.0

          k = 0
          DO j = 1, nds
            ik = ikds(j)
            ir = irds(j)
            IF (idring(ir).EQ.BOUNDARY) CYCLE

            k = k + 1
            sputter_data(i)%ik (k) = ik
	    sputter_data(i)%ir (k) = ir
	    sputter_data(i)%id (k) = wallindex(j)
	    sputter_data(i)%r  (k) = rp(j)
	    sputter_data(i)%z  (k) = zp(j)
	    sputter_data(i)%dds(k) = dds2(j)
	    sputter_data(i)%te (k) = kteds(j)
	    sputter_data(i)%ti (k) = ktids(j)

            iz = sputter_data(i)%charge
	    sputter_data(i)%flux(k,iz) = 
     .        knds(j) * ABS(kvds(j)) / kbfs(ik,ir) * costet(j) *   ! same as in NEUT.f
     .        (sputter_data(i)%fraction / 100.0)

            DO iz = 1, sputter_data(i)%atomic_number
    	      sputter_data(i)%e0(k,iz) = 
     .          3.0 * kteds(j) * REAL(iz) +   ! Missing contribution from ion velocity at sheath entrance...
     .          2.0 * ktids(j) 
            ENDDO
          ENDDO      

        ELSE

          fname = TRIM(resdir)//'/'//
     .            TRIM(sputter_data(i)%case_name)//'.'//
     .            TRIM(sputter_data(i)%extension)
c...      Copy file to run directory:
          command = 'cp -p '//TRIM(fname)//' .'
          WRITE(0,*) 'command: ',TRIM(fname)
          CALL CIssue(TRIM(command),status)
          IF (status.NE.0) 
     .      CALL ER('divCompileSputteringYields','Unable to copy '//
     .              'data file',*99)
	  
          fname = TRIM(sputter_data(i)%case_name)
          j = 1
          DO k = 1, LEN_TRIM(fname)
            IF (fname(k:k).EQ.'/') j = k
          ENDDO
          fname = fname(j+1:LEN_TRIM(fname))
          fname = TRIM(fname)//'.'//TRIM(sputter_data(i)%extension)
          WRITE(0,*) 'fname: '//TRIM(fname)
	  
c...      Unzip file if necessary:
          IF (fname(LEN_TRIM(fname)-2:LEN_TRIM(fname)).EQ.'zip'.OR.
     .        fname(LEN_TRIM(fname)-1:LEN_TRIM(fname)).EQ.'gz') 
     .      CALL UnzipFile(fname)
	  
          CALL find_free_unit_number(fp)
          OPEN(UNIT=fp,FILE=TRIM(fname),ACCESS='SEQUENTIAL',
     .         FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ierr,ERR=97)            
	  
          READ(fp) sputter_data(i)%version,
     .             sputter_data(i)%format, 
     .             sputter_data(i)%type
	  
          IF (sputter_data(i)%version.EQ.1.0) THEN
	  
            READ(fp) sputter_data(i)%absfac,         ! ABSFAC
     .               sputter_data(i)%atomic_number,  ! CION REAL()
     .               sputter_data(i)%atomic_mass,    ! CRMI
     .               sputter_data(i)%charge_min,     ! CIZSC
     .               sputter_data(i)%charge_max,     ! NIZS
     .               sputter_data(i)%nsegments,      ! NDS
     .               sputter_data(i)%ncore           ! number of radial core data points
	    
            SELECTCASE (sputter_data(i)%format)
c             ----------------------------------------------------------
              CASE (1)
                nseg = sputter_data(i)%nsegments
                mion = MIN(sputter_data(i)%atomic_number,
     .                     sputter_data(i)%charge_max    )
                ALLOCATE(sputter_data(i)%ik   (nseg))        ! this is stupid ... just make these a fixed size and avoid this pointless dynamic allocation
                ALLOCATE(sputter_data(i)%ir   (nseg))
                ALLOCATE(sputter_data(i)%id   (nseg))
                ALLOCATE(sputter_data(i)%r    (nseg))
                ALLOCATE(sputter_data(i)%z    (nseg))
                ALLOCATE(sputter_data(i)%dds  (nseg))
                ALLOCATE(sputter_data(i)%te   (nseg))
                ALLOCATE(sputter_data(i)%ti   (nseg))
                ALLOCATE(sputter_data(i)%flux (nseg,mion))
                ALLOCATE(sputter_data(i)%e0   (nseg,mion))
                ALLOCATE(sputter_data(i)%yield(nseg,mion))
                ncore = sputter_data(i)%ncore
                ALLOCATE(sputter_data(i)%core_rho (ncore))
                ALLOCATE(sputter_data(i)%core_psin(ncore))
                ALLOCATE(sputter_data(i)%core_ne  (ncore))
                ALLOCATE(sputter_data(i)%core_te  (ncore))
                ALLOCATE(sputter_data(i)%core_ti  (ncore))
                ALLOCATE(sputter_data(i)%core_percent_nfrac(ncore))
                ALLOCATE(sputter_data(i)%core_percent_efrac(ncore))
	  
                WRITE(0,*) 'ncore',ncore,nseg,mion
	  
                DO id = 1, nseg
                  READ(fp,END=10) 
     .              sputter_data(i)%ik (id),
     .              sputter_data(i)%ir (id),
     .              sputter_data(i)%id (id),
     .              sputter_data(i)%r  (id),
     .              sputter_data(i)%z  (id),
     .              sputter_data(i)%dds(id),
     .              sputter_data(i)%te (id),
     .              sputter_data(i)%ti (id)
                  DO iz = 1, mion
                    READ(fp) 
     .                sputter_data(i)%flux (id,iz),
     .                sputter_data(i)%e0   (id,iz),
     .                sputter_data(i)%yield(id,iz)
                  ENDDO
                ENDDO
c               Data on the impurity distribution in the core, for rescaling the 
c               wall flux later, if desired:
                ncore = sputter_data(i)%ncore
                IF (ncore.GT.0) THEN 
                  READ(fp) 
     .              sputter_data(i)%core_rho (1:ncore),
     .              sputter_data(i)%core_psin(1:ncore),
     .              sputter_data(i)%core_ne  (1:ncore),	  
     .              sputter_data(i)%core_te  (1:ncore), 
     .              sputter_data(i)%core_ti  (1:ncore), 
     .              sputter_data(i)%core_percent_nfrac(1:ncore),
     .              sputter_data(i)%core_percent_efrac(1:ncore)
                ENDIF
c             ----------------------------------------------------------
              CASEDEFAULT
                CALL ER('divCompileSputteringYields','Unrecognized '//
     .                  'file TYPE',*99)
c             ----------------------------------------------------------
            ENDSELECT
          ELSE
            CALL ER('divCompileSputteringYields','Unrecognized '//
     .              'file VERSION number',*99)

          ENDIF
	  
10        CONTINUE
          IF (id.NE.nseg+1) sputter_data(i)%nsegments = id - 1
          CLOSE(fp)

        ENDIF

      ENDDO
c
c     ------------------------------------------------------------------
c     DERIVE YIELDS
c
      DO i = 1, nsputter
        nseg = sputter_data(i)%nsegments
        mion = MIN(sputter_data(i)%atomic_number,
     .             sputter_data(i)%charge_max    )

        IF (sputter_data(i)%type.EQ.3) CYCLE   ! yield is already calculated

        sputter_data(i)%target_number = cion
 
c        matt = get_target_index(cion)
c        matp = get_plasma_index(    sputter_data(i)%atomic_number ,
c     .                          INT(sputter_data(i)%atomic_mass  ))

        CALL SYLD96(matt,matp,cneutd,-1,-1,cion,
     .              sputter_data(i)%atomic_number,
     .              sputter_data(i)%atomic_mass  ,rdum1)

        WRITE(0,*) 'yield 1',i,matt,matp,cizb,
     .              sputter_data(i)%atomic_number

        IF (matt.EQ.-1) 
     .    CALL ER('divCompileSputteringYields','Target material not '//
     .            'identified',*99)
        IF (matp.EQ.-1) 
     .    CALL ER('divCompileSputteringYields','Plasma material not '//
     .            'identified',*99)

c       Standard setup if not Be (cizb=4) incident on W (MATT=9) 
c         (MATP = -1 for Be, i.e. it's not in the standard setup at the moment)
        IF (sputter_data(i)%atomic_number.NE.4.OR.matt.NE.9) 
     .    CALL init_eckstein_2007(matt,matp)

        WRITE(0,*) 'yield 2',i,matt,matp

        DO j = 1, nseg
          DO iz = 1, mion
            IF (sputter_data(i)%atomic_number.EQ.4.AND.matt.EQ.9) THEN
              sputter_data(i)%yield(j,iz) = 
     .          divYield_Be_W(sputter_data(i)%e0(j,iz))
            ELSE
              sputter_data(i)%yield(j,iz) = 
     .          yield_2007(matp,matt,sputter_data(i)%e0(j,iz))
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c
c     ------------------------------------------------------------------
c     SET absfac FOR SPECIFIED IMPURITY FRACTIONS
c
      DO i = 1, nsputter

        IF (sputter_data(i)%type    .NE. 1  .OR.
     .      sputter_data(i)%absfac  .NE. 1.0.OR.
     .      sputter_data(i)%fraction.EQ.-1.0) CYCLE
        ncore = sputter_data(i)%ncore
        j = 0

c        DO j = 1, ncore
c          ratio = sputter_data(i)%core_rho(j    ) / 
c     .            sputter_data(i)%core_rho(ncore)
c          IF (ratio.LE.0.95.AND.ratio.GT.0.95) EXIT
c        ENDDO
c        IF (j.EQ.ncore+1) 
c     .    CALL ER('divCompileSputteringData','q95 not found',*99)  ! not really q95?

        frac = sputter_data(i)%fraction / 
     .         sputter_data(i)%core_percent_nfrac(ncore-3) 

        write(0,'(A,3I6,1P,4E10.2,0P)') 
     .    'process',i,j,ncore,frac,sputter_data(i)%fraction,
     .    sputter_data(i)%core_percent_nfrac(ncore-3),
     .    sputter_data(i)%fraction /
     .    sputter_data(i)%core_percent_nfrac(ncore-3) 

        sputter_data(i)%absfac = frac

      ENDDO
c
c     ------------------------------------------------------------------
c     ASSIGN wlprob
c
      WRITE(0 ,*) 'sputter data',nsputter
      WRITE(88,*) 'sputter data',nsputter
      DO i = 1, nsputter
        WRITE(88,'(F5.1,I4,1P,E10.2,0P,I4,F6.2,3I4,2I6)') 
     .    sputter_data(i)%version,
     .    sputter_data(i)%format   ,
     .    sputter_data(i)%absfac ,
     .    sputter_data(i)%atomic_number,
     .    sputter_data(i)%atomic_mass,
     .    sputter_data(i)%target_number,
     .    sputter_data(i)%charge_min,   
     .    sputter_data(i)%charge_max,   
     .    sputter_data(i)%nsegments,
     .    sputter_data(i)%type 
          nseg = sputter_data(i)%nsegments
          mion = MIN(sputter_data(i)%atomic_number,
     .               sputter_data(i)%charge_max    )
        DO id = 1, nseg
          WRITE(88,'(5X,F5.2,3I6,2F7.3,1P,E10.2,0P,2F8.2,
     .              20(2X,1P,E10.2,0P,F8.2,1P,E10.2,0P))')
     .      kmfpws(sputter_data(i)%id(id)),
     .      sputter_data(i)%ik (id),
     .      sputter_data(i)%ir (id),
     .      sputter_data(i)%id (id),
     .      sputter_data(i)%r  (id),
     .      sputter_data(i)%z  (id),
     .      sputter_data(i)%dds(id),
     .      sputter_data(i)%te (id),
     .      sputter_data(i)%ti (id),
     .      (sputter_data(i)%flux (id,iz),
     .       sputter_data(i)%e0   (id,iz),
     .       sputter_data(i)%yield(id,iz),iz=1,mion)
        ENDDO
        ncore = sputter_data(i)%ncore
        IF (ncore.GT.0) THEN 
          WRITE(0,*) 'doing the ncore dance, again'
          DO id = 1, ncore
            WRITE(88,'(2F10.3,1P,E10.2,0P,2F9.2,2F15.7)')
     .        sputter_data(i)%core_rho (id),
     .        sputter_data(i)%core_psin(id),
     .        sputter_data(i)%core_ne  (id),	  
     .        sputter_data(i)%core_te  (id), 
     .        sputter_data(i)%core_ti  (id), 
     .        sputter_data(i)%core_percent_nfrac(id),
     .        sputter_data(i)%core_percent_efrac(id)
          ENDDO
        ENDIF
      ENDDO


      wlprob  = 0.0
      nwlprob = wallpts
      DO id = 1, nwlprob
        wlprob(id,1) = id
        wlprob(id,2) = id
      ENDDO

      DO i = 1, nsputter
        nseg = sputter_data(i)%nsegments
        mion = MIN(sputter_data(i)%atomic_number,
     .             sputter_data(i)%charge_max    )

c       Calcualte the particle influx from each wall segment in [s-1 m-1 toroidally]:
        DO j = 1, nseg
          id = sputter_data(i)%id(j)
          DO iz = 1, mion
            wlprob(id,3) = wlprob(id,3) +
     .        sputter_data(i)%absfac      *    
     .        sputter_data(i)%flux (j,iz) *    
     .        sputter_data(i)%yield(j,iz) *
     .        sputter_data(i)%dds  (j   ) *
     .        kmfpws(id)
          ENDDO
        ENDDO

c       Add up all the segments to get ABSFAC value for all of the 
c       sputtering species processed so far in this loop (because
c       WLPROB is an integral over all sputtering species):
        sputter_data(i)%absfac_net = 0.0
        DO j = 1, nwlprob
          sputter_data(i)%absfac_net = sputter_data(i)%absfac_net + 
     .                                 wlprob(j,3)
        ENDDO

      ENDDO

c     Isolate the net influx for each individual species:
      DO i = nsputter, 2, -1
        sputter_data(i)%absfac_net = sputter_data(i  )%absfac_net - 
     .                               sputter_data(i-1)%absfac_net
      ENDDO

      DO i = 1, nsputter
        WRITE(0 ,*) 'nabsfac=',i,sputter_data(i)%absfac_net,nwlprob
        WRITE(88,*) 'nabsfac=',i,sputter_data(i)%absfac_net,nwlprob
      ENDDO

c     Add up the influxes for the individual sputtering species and 
c     assign the total influx over-ride value for DIVIMP:
      nabsfac = SUM(sputter_data(1:nsputter)%absfac_net)

c     Normalize the wall launch probabilities:
      wlprob(:,3) = wlprob(:,3) / nabsfac

      wallpt(:,13) = 0.0
c
c     ------------------------------------------------------------------     
c     Send data to IDL:

      CALL inOpenInterface('idl.divimp_sputter_data',ITF_WRITE)
      i = nsputter
      CALL inPutData(nabsfac,'ABSFAC_TOTAL','s-1')	  
      CALL inPutData(sputter_data(1:i)%absfac_net   ,'ABSFAC'  ,'s-1')	  
      CALL inPutData(sputter_data(1:i)%atomic_number,'BOMB_Z'  ,'NA')
      CALL inPutData(sputter_data(1:i)%atomic_mass  ,'BOMB_A'  ,'NA')	  
      CALL inPutData(sputter_data(1:i)%target_number,'TARGET_Z','NA') 
      CALL inPutData(sputter_data(1:i)%charge_min   ,'C_MIN'   ,'NA')   
      CALL inPutData(sputter_data(1:i)%charge_max   ,'C_MAX'   ,'NA')  
      CALL inPutData(sputter_data(1:i)%type         ,'TYPE'    ,'NA')   
c     Assume the data all comes from the same geometry and plasma:
      n = sputter_data(1)%nsegments
      CALL inPutData(kmfpws(sputter_data(1)%id(1:n)  ),'FACT','NA')
      CALL inPutData(wlprob(sputter_data(1)%id(1:n),3),'PROB','NA')
      CALL inPutData(sputter_data(1)%id (1:n),'ID'    ,'NA')	   
      CALL inPutData(sputter_data(1)%r  (1:n),'POS_R' ,'m' ) 
      CALL inPutData(sputter_data(1)%z  (1:n),'POS_Z' ,'m' )   
      CALL inPutData(sputter_data(1)%dds(1:n),'LENGTH','m' )  
c     Take IK,IR and Te,i for the first data set from standard DIVIMP 
c     sputtering (rather than EIRENE calculated yields), otherwise 
c     just take what's in the first data set, which is probably zero:
      DO i = 1, nsputter
        IF (sputter_data(i)%type.EQ.1) THEN
          write(0,*) 'debug: dumping Te,i from',i
          CALL inPutData(sputter_data(i)%ik(1:n),'IK','eV')   
          CALL inPutData(sputter_data(i)%ir(1:n),'IR','eV')	  
          CALL inPutData(sputter_data(i)%te(1:n),'TE','eV')   
          CALL inPutData(sputter_data(i)%ti(1:n),'TI','eV')	  
          EXIT
        ENDIF
      ENDDO
      IF (i.EQ.nsputter+1) THEN
        CALL inPutData(sputter_data(1)%ik(1:n),'IK','eV')   
        CALL inPutData(sputter_data(1)%ir(1:n),'IR','eV')	  
        CALL inPutData(sputter_data(1)%te(1:n),'TE','eV')   
        CALL inPutData(sputter_data(1)%ti(1:n),'TI','eV')	  
      ENDIF
      DO i = 1, nsputter
        WRITE(t1,'(A,I0.2)') '_',i
        mion = MIN(sputter_data(i)%atomic_number,
     .             sputter_data(i)%charge_max    )
        CALL inPutData(mion,'MAX_CHARGE'//t1,'NA')	  
        DO iz = 1, mion
          WRITE(t2,'(2A,I0.2)') t1,'_',iz
          CALL inPutData(sputter_data(i)%flux (1:n,iz),'FLUX'//t2,'s-1')
          CALL inPutData(sputter_data(i)%e0   (1:n,iz),'E0'//t2   ,'eV')	  
          CALL inPutData(sputter_data(i)%yield(1:n,iz),'YIELD'//t2,'NA')	  
        ENDDO
      ENDDO
      CALL inCloseInterface

      WRITE(0,*) 'nabsfac total=',nabsfac,nwlprob

c      stop 'here in the shit'

      RETURN
97    WRITE(0,*) 'OPEN error, IOSTAT=',ierr
      STOP 'shiiiit'
99    WRITE(0,*) '  FILE NAME = ',TRIM(fname)
      STOP
      END

c
c ======================================================================
c
      SUBROUTINE ImportOSMGrid
      USE mod_geometry
      USE mod_sol28_global
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'pindata'

      INTEGER GetObject

      INTEGER        status,n,ic,it,iobj,ik,ir,id,in
      REAL*8         rv(4),zv(4)
      CHARACTER*1024 fname,command

c...  Load reference OSM data file (binary form):
      fname = TRIM(opt%f_osm_load)
      n = LEN_TRIM(fname)
      command = 'cp '//TRIM(opt%f_osm_dir)//TRIM(fname)//' .'
      CALL CIssue(TRIM(command),status)
      CALL UnzipFile(fname)
      CALL LoadGrid (fname)

c...  Load the OSM geometry file:
      fname = TRIM(opt%f_osm_load)
      n = LEN_TRIM(fname)
      IF (fname(n-2:n).EQ.'zip') fname(n-6:n) = 'ogd.zip'
      IF (fname(n-1:n).EQ.'gz' ) fname(n-5:n) = 'ogd.gz'
      IF (fname(n-2:n).EQ.'osm') fname(n-2:n) = 'ogd'
c      WRITE(0,*) 'fname>'//TRIM(fname)//'<'
      command = 'cp '//TRIM(opt%f_osm_dir)//TRIM(fname)//' .'
      CALL CIssue(TRIM(command),status)
      CALL UnzipFile  (fname)
      CALL LoadObjects(fname,status)
      IF (status.LT.0) CALL ER('ImportOSMGrid','Halting DIVIMP',*99)

c...
      ALLOCATE(d_rvertp(5,MAXNKS*MAXNRS))
      ALLOCATE(d_zvertp(5,MAXNKS*MAXNRS))
      id = 0
      DO it = 1, ntube
        ir = it
        ik = 0
        DO ic = tube(it)%cell_index(1), tube(it)%cell_index(2) 
          id = id + 1
          ik = ik + 1          
          iobj = GetObject(ic,IND_CELL)
          DO in = 1, 4
            CALL GetVertex(iobj,in,rv(in),zv(in))
          ENDDO
          nvertp(id) = 4
          d_rvertp(1:4,id) = rv
          d_zvertp(1:4,id) = zv
          korpg(ik,ir) = id
          rs(ik,ir) = cell(ic)%cencar(1)
          zs(ik,ir) = cell(ic)%cencar(2)
          bratio(ik,ir) = field(ic)%bratio 
          kbfs  (ik,ir) = 1.0 / bratio(ik,ir)
        ENDDO
        IF (it.LT.grid%isep) THEN
          ik = ik + 1          
          id = id + 1
          nvertp(id) = 4
          d_rvertp(1:4,id) = d_rvertp(1:4,korpg(1,ir))
          d_zvertp(1:4,id) = d_zvertp(1:4,korpg(1,ir))
          korpg(ik,ir) = id
          rs(ik,ir) = rs(1,ir)
          zs(ik,ir) = zs(1,ir)
          bratio(ik,ir) = bratio(1,ir)
          kbfs  (ik,ir) = kbfs  (1,ir)
        ENDIF
        nks(ir) = ik
        psitarg(ir,:) = tube(it)%psin
      ENDDO
      rvertp = SNGL(d_rvertp)
      zvertp = SNGL(d_zvertp)

c...  Add virtual rings:
c      WRITE(0,*) 'GRID%IPFZ=',grid%ipfz

      irsep  = grid%isep
      irwall = grid%ipfz - 1
      irtrap = grid%ipfz
      nrs    = ntube

      ikto = grid%ikto
      ikti = grid%ikti

      npolyp  = id
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

      r0  = SNGL(grid%r0)
      z0  = SNGL(grid%z0)
      rxp = rvertp(4,korpg(ikto,irsep))
      zxp = zvertp(4,korpg(ikto,irsep))

c      WRITE(0,*) 'whoa=',r0,z0,rxp,zxp
c      stop 'dsdfsf'

      CALL OutputData(85,'Linear 1')
       
      IF (irwall.EQ.nrs) THEN 
        cgridopt = LINEAR_GRID
        DO ir = 2, nrs
          IF (nks(ir).NE.nks(1)) THEN 
            cgridopt = RIBBON_GRID
            EXIT
          ENDIF
        ENDDO
      ELSE
        cgridopt = 3
      ENDIF
c      WRITE(0,*) 'cgridopt=',cgridopt

      CALL InsertRing(1     ,BEFORE,PERMANENT)
      CALL InsertRing(irwall,AFTER ,PERMANENT)
      IF (cgridopt.EQ.LINEAR_GRID.OR.cgridopt.EQ.RIBBON_GRID) THEN 
      ELSE
        CALL InsertRing(irtrap,BEFORE,PERMANENT)
      ENDIF
      idring(1     ) = BOUNDARY
      idring(irtrap) = BOUNDARY
      idring(irwall) = BOUNDARY
      psitarg(1     ,:) = 0.0
      psitarg(irtrap,:) = 0.0
      psitarg(irwall,:) = 0.0

      CALL OutputData(86,'Linear 2')

      cutpt1   = ikto
      cutpt2   = ikti             ! These are semi-bogus for a connected double-null...?
      cutring  = irsep - 1
      maxkpts  = nks(irsep)
      maxrings = irwall

c...  Here to avoid the need to call TailorGrid:
c      IF (grdnmod.GT.0) CALL TailorGrid

      CALL SetupGrid
      CALL SequenceGrid
      CALL FindGridBreak
      CALL SetupGrid

c      CALL OutputData(85,'JUST FINISHED LOADING OSM GRID')
c      CALL DumpGrid  ('LOADING OSM GRID')

c...  Add virtual boundary cells, which will be stripped off later:
      IF (CTARGOPT.EQ.0.OR.CTARGOPT.EQ.1.OR.CTARGOPT.EQ.2.OR.
     .    CTARGOPT.EQ.3.OR.CTARGOPT.EQ.6) 
     .   CALL AddPoloidalBoundaryCells

c...  Clear arrays:
      CALL geoClean
      CALL osmClean

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE SetupSOL28
      USE mod_sol28
      USE mod_sol28_global
      USE mod_legacy
      USE mod_geometry
      IMPLICIT none

c      RETURN

      CALL MapRingstoTubes
      CALL AssignOSMWall

      CALL SaveGeometryData('osm_geometry.raw')

      CALL LoadLegacyData('osm_legacy.raw')

      CALL GenerateTubeGroups
      CALL DumpData_OSM('output.grid_tubes','Done analysing tubes')

      CALL GenerateTargetGroups
      CALL DumpData_OSM('output.grid_targets','Done analysing targets')

c      CALL SetTargetConditions(itube)

      IF (opt%osm_load.NE.0) CALL LoadReferenceSolution(1)

c...  Clear geometry arrays:
      CALL geoClean

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE DumpData(cdum1,cdum2)
      IMPLICIT none

      CHARACTER, INTENT(IN) :: cdum1*(*),cdum2*(*)

      CALL OutputData(85,cdum2)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE ExecuteSOL28(irstart,irend,ikopt,sloutput)
      USE mod_sol28
      USE mod_sol28_global
      USE mod_geometry
      USE mod_legacy
      IMPLICIT none

      INTEGER irstart,irend,ikopt
      LOGICAL sloutput

      INTEGER itube,itube1,itube2,status
      LOGICAL cont

c...  Need to reload geometry data in case Eirene has been called (hopefully
c     this won't be required in the future...):
 
      CALL LoadObjects('osm_geometry.raw',status)
      IF (status.NE.0) CALL ER('ExecuteSOL28','Unable to load '//
     .                         'geometry data',*99)

c...  Load up PIN data if available:
      IF (opt%pin_data) CALL MapNeutralstoTubes

c      CALL MapRingstoTubes
c      CALL LoadLegacyData('osm_legacy.raw')
c      CALL SetTargetConditions

      itube1 = 0
      itube2 = 0
      DO itube = 1, ntube
        IF (tube(itube)%ir.EQ.irstart) itube1 = itube
        IF (tube(itube)%ir.EQ.irend  ) itube2 = itube
      ENDDO
      IF (itube2.EQ.0.AND.irend.EQ.ntube+2) itube2 = ntube
      IF (itube1.EQ.0.OR.itube2.EQ.0) 
     .  CALL ER('ExecuteSOL28','Tube index error',*99)

c...  Call SOL28 plasma solver:
      CALL MainLoop(itube1,itube2,ikopt,sloutput)

c...  Generate output files:
c      CALL GenerateOutputFiles  ! this is now done at the end of bgplasma

c...  Fill DIVIMP arrays:
      CALL MapTubestoRings(irstart,irend)

c...  Clear geometry arrays:
      CALL geoClean

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE CloseSOL28
      USE mod_geometry
      USE mod_sol28_global
      IMPLICIT none

c      RETURN
      INTEGER status

      CALL LoadObjects('osm_geometry.raw',status)
      IF (status.NE.0) CALL ER('ExecuteSOL28','Unable to load '//
     .                         'geometry data',*99)

c...  Generate output files:
      CALL GenerateOutputFiles

c...  Save solution:
c      CALL SaveGrid('osm.raw')

c...  Clear memory:
      CALL osmClean
      CALL geoClean

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignNodeValues(itube,nnode,mnode,node)
      USE mod_sol28
      USE mod_sol28_global
      IMPLICIT none

      INTEGER itube,nnode,mnode       
      TYPE(type_node) :: node(*)

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL GetJsat

      INTEGER ik,ir,i1,id,id1,id7
      LOGICAL new
      REAL    te(0:6),ne(0:6),nf,s(0:6)

      STOP 'ROUTINE TAGGED FOR DELETION'
    

      ir = tube(itube)%ir

      new = .TRUE.

c     get rid of do loop, put job of finding where in the SOL28 input listing should be used into a subroutine

      CALL FindS28Parameters_V3(ir,te,ne,nf,s,new)



      IF (.TRUE.) THEN

        s(0) = 0.0
        s(6) = ksmaxs(ir)

        DO i1 = 0, 6
          DO ik = 1, nks(ir)
            IF (s(i1).GE.ksb(ik-1,ir).AND.s(i1).LE.ksb(ik,ir)) THEN
              node(i1+1)%icell = ik
              IF (i1.EQ.3) s(3) = kss(ik,ir)
            ENDIF
          ENDDO
        ENDDO

        nnode = 7
        mnode = 4

c...    Assign values to nodes:
        node(1:7)%s  = s (0:6)
        node(1:7)%ne = ne(0:6)
        node(1:7)%te = te(0:6)
c...    Assign other quantites:
        node(1:7)%jsat(1) = 0.0
        node(1:7)%pe      = 0.0
        node(1:7)%ni(1)   = 0.0
        node(1:7)%pi(1)   = 0.0
        node(1:7)%ti(1)   = 0.0
        node(1:7)%machno  = 0.0
        node(1:7)%epot    = 0.0
        node(1:7)%efield  = 0.0

        id1 = idds(ir,2)
        id7 = idds(ir,1)
        IF (node(1)%ne.EQ.0.0) node(1)%ne = knds (id1)
        IF (node(7)%ne.EQ.0.0) node(7)%ne = knds (id7)
        IF (node(1)%te.EQ.0.0) node(1)%te = kteds(id1)
        IF (node(7)%te.EQ.0.0) node(7)%te = kteds(id7)
        IF (node(1)%ti(1).EQ.0.0) node(1)%ti(1) = ktids(id1)
        IF (node(7)%ti(1).EQ.0.0) node(7)%ti(1) = ktids(id7)

c REAL FUNCTION GetJsat(te,ti,ne,v)
        id = id1
        node(1)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(1)%jsat(1) = -ABS(node(1)%jsat(1))
        id = id7
        node(7)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(7)%jsat(1) = -ABS(node(7)%jsat(1))
      ENDIF


      RETURN
99    CONTINUE
      WRITE(0,*) 'IK,IR=',ik,ir,i1,osms28(i1,1)
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE MapNeutralstoTubes
      USE mod_sol28_global
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ion,ncell1,ike,ir,cind1,cind2

c...  Copy PIN data:
      ion = 1
      ncell1 = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1

        cind1 = ncell1 + 1
        cind2 = ncell1 + ike 

        pin(cind1:cind2,ion)%ion = pinion(1:ike,ir)
        pin(cind1:cind2,ion)%rec = pinrec(1:ike,ir)
        pin(cind1:cind2,ion)%mom = pinmp (1:ike,ir)

        ncell1 = cind2
      ENDDO

      IF (ncell.NE.ncell1) 
     .  CALL ER('MapNeutraltoTubes','NCELL1.NE.NCELL',*99)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE AssignOSMWall
      USE mod_sol28_wall
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'walls_com'

      INTEGER iwall

      nwall = wallpts

      ALLOCATE(wall(nwall))

      DO iwall = 1, nwall
        wall(iwall)%class             =  1       
        wall(iwall)%index(WAL_GROUP ) = -1
        wall(iwall)%index(WAL_INDEX ) = -1    
        wall(iwall)%index(WAL_TUBE  ) = -1    
        wall(iwall)%index(WAL_TARGET) = NINT(wallpt(iwall,16))    
        wall(iwall)%material_tag      = 'no_set'
        wall(iwall)%material          = -1
        wall(iwall)%temperature       =      wallpt(iwall,19)
        wall(iwall)%v1(1)             = DBLE(wallpt(iwall,20))
        wall(iwall)%v1(2)             = DBLE(wallpt(iwall,21))
        wall(iwall)%v2(1)             = DBLE(wallpt(iwall,22))
        wall(iwall)%v2(2)             = DBLE(wallpt(iwall,23))
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE MapRingstoTubes
      USE mod_sol28_global
      USE mod_grid_divimp
      USE mod_geometry
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

      INTEGER ik,ike,ir,id,cind1,cind2,i1,i2,iside,ik1,ir1,ir2,ion,iobj,
     .        fp,idum1
      LOGICAL pfz_ring,load_bfield_data
      REAL    rdum1
      REAL, ALLOCATABLE:: bfield_data(:,:,:)
      REAL*8  a(3)
      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj

      REAL tube_3D_data(5,MAXNKS,MAXNRS),dangle

      load_bfield_data = .FALSE.

      IF (load_bfield_data) THEN
        ALLOCATE(bfield_data(MAXNKS,MAXNRS,5))
        fp = 99
        OPEN(UNIT=fp,FILE='objects.bfield',ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=95)      
        READ(fp,*)
        READ(fp,*)
        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          DO ik = 1, nks(ir)
            READ(fp,*) (rdum1,i1=1,4),ik1,ir1,rdum1,
     .                 bfield_data(ik,ir,1:5)
            IF (ik1.NE.ik.OR.ir1.NE.ir) 
     .        CALL ER('MapRingsToTubes','Index mismatch',*99)
          ENDDO
        ENDDO
        CLOSE(fp)
      ENDIF

c *CRUDE* Should be elsewhere... can IDRING take the place of RINGTYPE for EIRENE04? 
      ringtype = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) THEN
          ringtype(ir) = BOUNDARY
          CYCLE
        ENDIF
        IF (ir.LT.irsep                  ) ringtype(ir) = CORE
        IF (ir.GE.irsep .AND.ir.LT.irwall) ringtype(ir) = SOL1
        IF (ir.GT.irtrap.AND.ir.LE.nrs   ) ringtype(ir) = PFZ
      ENDDO
c...  Special for secondary PFZ in generalized grids:
      IF (grdnmod.NE.0) THEN
        DO ik = 1, nks(irwall)
c...      Scan over IRWALL and check all rings that point to IRWALL to
c         see if *all* of that ring is pointing at IRWALL, in which
c         case it must be a PFZ ring (could just be a discontinuous
c         SOL target otherwise):
          IF (irins(ikins(ik,irwall),irins(ik,irwall)).EQ.irwall.AND.
     .        ringtype(irins(ik,irwall)).NE.PFZ) THEN
            ir1 = irins(ik,irwall)  
            pfz_ring = .TRUE.
            DO ik1 = 1, nks(ir1)
              IF (irins(ik1,ir1).NE.irwall) pfz_ring = .FALSE.
            ENDDO
            IF (pfz_ring) THEN
              DO ir2 = ir1, nrs
                ringtype(ir2) = PFZ
                IF     (irouts(1,ir2).EQ.irouts(nks(ir2),ir2))THEN
                ELSEIF (irouts(1       ,ir2).EQ.irwall.OR.
     .                  irouts(nks(ir2),ir2).EQ.irwall) THEN
c...              Discontinuous target!
                  STOP 'THIS IS A PROBLEM'
                ELSE
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDIF

        ENDDO
      ENDIF

c...  Count number of cells:
      ntube = 0
      ncell = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ntube = ntube + 1
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1
        ncell = ncell + ike 
      ENDDO

c...  Declare global arrays:
      nfield    = ncell
      npin      = ncell
      nphoton   = 1
      ndrift    = 1
      nkinetic  = 1     
      nfluid    = ncell
      nimpurity = 1
      ALLOCATE(tube    (ntube ))
      ALLOCATE(tube2   (ntube ))
      tube2(:)%state        = 0
      tube2(:)%target_pe(1) = 0
      tube2(:)%target_pe(2) = 0
c      ALLOCATE(tube_state(ntube))
c      tube_state = 0
      ALLOCATE(cell    (ncell ))
      ALLOCATE(field   (nfield))
      ALLOCATE(pin     (npin     ,nion))
      ALLOCATE(photon  (nphoton  ,nion))
      ALLOCATE(drift   (ndrift   ,nion))
      ALLOCATE(kinetic (nkinetic ,nion))
      ALLOCATE(fluid   (nfluid   ,nion))
      ALLOCATE(impurity(nimpurity,nion)) 
c...  Reference plasma solution:
      ref_nion   = 1
      ref_ntube  = 1
      ref_nfluid = 1
      ALLOCATE(ref_tube(ref_ntube))
      ALLOCATE(ref_fluid(ref_nfluid,ref_nion))
c...  Copy DIVIMP grid:
      ion = 1
      ntube = 0
      ncell = 0

      tube_3D_data = 0.0
      dangle = 0.0
c      CALL CalcTubeDimensions(tube_3D_data,dangle)

      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1

        cind1 = ncell + 1
        cind2 = ncell + ike

        ntube = ntube + 1
        tube(ntube)%index  = ntube
        tube(ntube)%ir     = ir
        tube(ntube)%type   = ringtype(ir)
        tube(ntube)%ikti   = ikti2(ir)
        tube(ntube)%ikto   = ikto2(ir)
        tube(ntube)%n      = ike
        tube(ntube)%smax   = ksmaxs(ir)
        tube(ntube)%pmax   = kpmaxs(ir)
        tube(ntube)%dangle = dangle
        tube(ntube)%rho  = rho(ir,CELL1)
        tube(ntube)%psin = psitarg(ir,1)
        tube(ntube)%cell_index(1) = cind1  ! The "1" must be consistent with "LO" in mod_osm
        tube(ntube)%cell_index(2) = cind2  ! and same for the "2"...
        IF (ir.LT.irsep) THEN
          tube(ntube)%bratio = 0.0
          tube(ntube)%dds    = 0.0
          tube(ntube)%rp     = 0.0
          tube(ntube)%costet = 0.0
        ELSE
          id = idds(ir,2)
          tube(ntube)%bratio(1) = bratio(1,ir)
          tube(ntube)%dds   (1) = dds2(id)
          tube(ntube)%rp    (1) = rp(id)
          tube(ntube)%costet(1) = costet(id)
          tube(ntube)%metric(1) = thetat(id)
          id = idds(ir,1)
          tube(ntube)%bratio(2) = bratio(ike,ir)
          tube(ntube)%dds   (2) = dds2(id)
          tube(ntube)%rp    (2) = rp(id)
          tube(ntube)%costet(2) = costet(id)
          tube(ntube)%metric(2) = thetat(id)
        ENDIF
        
        ncell = cind2
        fluid(cind1:cind2,ion)%te = ktebs(1:ike,ir)
        fluid(cind1:cind2,ion)%ni = knbs (1:ike,ir)
        fluid(cind1:cind2,ion)%vi = kvhs (1:ike,ir)
        fluid(cind1:cind2,ion)%ti = ktibs(1:ike,ir)

        cell(cind1:cind2)%cencar(1) = rs(1:ike,ir)
        cell(cind1:cind2)%cencar(2) = zs(1:ike,ir)
        cell(cind1:cind2)%cencar(3) = 0.0
        cell(cind1:cind2)%vol       = kvols (1:ike,ir)


c        IF (cind1.LE.1025.AND.cind2.GE.1025) THEN
c          DO id = 1,ike
c            WRITE(0,*) 'VOLS:',id,ir,kvols(id,ir)
c          ENDDO
c          CALL OutputData(85,'Zero volume')
c          STOP 'sdkjlsdkfsd'
c        ENDIF

        cell(cind1:cind2)%s         = kss   (1:ike,ir)
        cell(cind1:cind2)%p         = kps   (1:ike,ir)
        cell(cind1:cind2)%sbnd(1)   = ksb   (0:ike-1,ir)
        cell(cind1:cind2)%sbnd(2)   = ksb   (1:ike  ,ir)
        cell(cind1:cind2)%pbnd(1)   = kpb   (0:ike-1,ir)
        cell(cind1:cind2)%pbnd(2)   = kpb   (1:ike  ,ir)
        cell(cind1:cind2)%metric    = thetag(1:ike,ir)

        cell(cind1:cind2)%sbnd_3D(1) = tube_3D_data(1,1:ike,ir)
        cell(cind1:cind2)%sbnd_3D(2) = tube_3D_data(2,1:ike,ir)
        cell(cind1:cind2)%s_3D       = tube_3D_data(3,1:ike,ir)
        cell(cind1:cind2)%area_3D    = tube_3D_data(4,1:ike,ir)
        cell(cind1:cind2)%volume_3D  = tube_3D_data(5,1:ike,ir)

        DO ik = 1, nks(ir)
c
c         jdemod - add if statement to test if dangle=0.0 
c                - divide by zero error otherwise
c
c          if (dangle.ne.0.0) then 
c            WRITE(6,'(A,2I6,4(2F12.6,2X))') 'CHECK 3D:',ik,ir,
c     .      cell(cind1+ik-1)%s      ,cell(cind1+ik-1)%s_3D      ,
c     .      cell(cind1+ik-1)%sbnd(1),cell(cind1+ik-1)%sbnd_3D(1),
c     .      cell(cind1+ik-1)%sbnd(2),cell(cind1+ik-1)%sbnd_3D(2),
c     .      cell(cind1+ik-1)%vol    ,cell(cind1+ik-1)%volume_3D *
c     .                               (360.0/dangle) 
c          else
c            WRITE(6,'(A,2I6,4(2F12.6,2X))') 'CHECK 3D:',ik,ir,
c     .      cell(cind1+ik-1)%s      ,cell(cind1+ik-1)%s_3D      ,
c     .      cell(cind1+ik-1)%sbnd(1),cell(cind1+ik-1)%sbnd_3D(1),
c     .      cell(cind1+ik-1)%sbnd(2),cell(cind1+ik-1)%sbnd_3D(2),
c     .      cell(cind1+ik-1)%vol    ,cell(cind1+ik-1)%volume_3D 
cc     .                               * (360.0/dangle) 
c          endif
c
c         end jdemod
c
        ENDDO

        field(cind1:cind2)%bratio = bratio(1:ike,ir)
        IF (load_bfield_data) THEN
          field(cind1:cind2)%b    = bfield_data(1:ike,ir,2)
          field(cind1:cind2)%br   = bfield_data(1:ike,ir,3)
          field(cind1:cind2)%bphi = bfield_data(1:ike,ir,4)
          field(cind1:cind2)%bz   = bfield_data(1:ike,ir,5)
        ELSE
          field(cind1:cind2)%b    = 0.0
          field(cind1:cind2)%br   = 0.0
          field(cind1:cind2)%bphi = 0.0
          field(cind1:cind2)%bz   = 0.0
        ENDIF

        DO ik = 1, ike
          cell(cind1+ik-1)%ik = ik
          cell(cind1+ik-1)%ir = ir
          cell(cind1+ik-1)%ds = ksb(ik,ir) - ksb(ik-1,ir)
        ENDDO
c        pin(cind1:cind2,ion)%ion = pinion(1:ike,ir)
c        pin(cind1:cind2,ion)%rec = pinrec(1:ike,ir)
c        pin(cind1:cind2,ion)%mom = pinmp (1:ike,ir)
      ENDDO
c...  Set grid quantities:
      grid%r0   = DBLE(r0)
      grid%z0   = DBLE(z0)
      grid%n    = ntube
      grid%isep = irsep  - 1
      grid%ipfz = irtrap - 2
      grid%ikto = ikto
      grid%ikti = ikti
c...  Geometry:
      ngrp = 1
      grp(1)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
      grp(1)%type   = GRP_QUADRANGLE
c...  
      nobj = 0
      nsrf = 0
      nvtx = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        ike = nks(ir)
        IF (ir.LT.irsep) ike = ike - 1

        DO ik = 1, ike
          id = korpg(ik,ir)         
          newobj%group         = ngrp
          newobj%index         = 0
          newobj%index(IND_IK) = ik
          newobj%index(IND_IR) = ir
          newobj%index(IND_IS) = 0
          newobj%index(IND_CELL   ) = nobj + 1
          newobj%index(IND_FLUID  ) = nobj + 1
          newobj%index(IND_KINETIC) = nobj + 1
          newobj%index(IND_NEUTRAL) = nobj + 1
          newobj%index(IND_FIELD  ) = nobj + 1
          newobj%segment(1) = 0
          newobj%phi        = 0.0
          newobj%nside      = 4
          DO iside = 1, 4
            newsrf%type = SPR_LINE_SEGMENT
            newsrf%obj  = nobj + 1
            newsrf%side = iside
            newsrf%nvtx = 2
            i1 = iside
            i2 = iside + 1
            IF (i2.EQ.5) i2 = 1
            IF (ALLOCATED(d_rvertp)) THEN
              a(1) = d_rvertp(i1,id)
              a(2) = d_zvertp(i1,id)
            ELSE
              a(1) = DBLE(rvertp(i1,id))
              a(2) = DBLE(zvertp(i1,id))
            ENDIF
            a(3) = 0.0D0
            newsrf%ivtx(1) = AddVertex(a) 
            IF (ALLOCATED(d_rvertp)) THEN
              a(1) = d_rvertp(i2,id)
              a(2) = d_zvertp(i2,id)
            ELSE
              a(1) = DBLE(rvertp(i2,id))
              a(2) = DBLE(zvertp(i2,id))
            ENDIF
            a(3) = 0.0D0
            newsrf%ivtx(2) = AddVertex(a) 
            newobj%iside(iside) = AddSurface(newsrf)
          ENDDO
          idum1 = AddObject(newobj)
        ENDDO
      ENDDO

c...  Build connection map:
      CALL BuildConnectionMap(1,nobj)
c...  Force idenfication of the radial boundaries in the OBJ connection
c     map, since generalised grids with poorly defined cross-field nearest
c     neighbours will come up with a 0 mapping, the same as for the grid
c     boundaries:
      DO iobj = 1, nobj
        ik =obj(iobj)%index(IND_IK)
        ir =obj(iobj)%index(IND_IR)
        IF (idring(irins (ik,ir)).EQ.BOUNDARY) obj(iobj)%omap(4) = -1
        IF (idring(irouts(ik,ir)).EQ.BOUNDARY) obj(iobj)%omap(2) = -1
        IF (ir.GE.irsep) THEN
          IF (ik.EQ.1      ) obj(iobj)%omap(1) = -2
          IF (ik.EQ.nks(ir)) obj(iobj)%omap(3) = -2
        ENDIF
      ENDDO


c...  Setup the old format SOL28 input data:
      osmns28 = osmnnode
      osms28 = 0.0
      DO i1 = 1, osmnnode
        osms28(i1,1) = osmnode(i1)%type 
        osms28(i1,5) = osmnode(i1)%rad_x
        osms28(i1,6) = osmnode(i1)%rad_y
      ENDDO

c...  Save legacy data:
      CALL SaveLegacyData('osm_legacy.raw')

c.... Clear arrays:
      IF (ALLOCATED(bfield_data)) DEALLOCATE(bfield_data)

      RETURN
 95   CALL ER('MapRingsToTubes','Error accessing B-field data file',*99)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE MapTubestoRings(irstart,irend)
      USE mod_sol28_global
      IMPLICIT none
 
      INTEGER irstart,irend,cind1,cind2,ike

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER ir,itube,ion,in

      ion = 1

      DO ir = irstart, irend
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        DO itube = 1, ntube
          IF (tube(itube)%ir.NE.ir) CYCLE

c...      Need to do volume weighted averaging...

          cind1 = tube(itube)%cell_index(1) 
          cind2 = tube(itube)%cell_index(2)

          IF (ir.GE.irsep) THEN
            in = idds(ir,2)
            knds (in) = tube(itube)%ne(1)      ! ne(1,ion) - 12/10/09 SL
            kvds (in) = tube(itube)%vi(1,ion) 
            kteds(in) = tube(itube)%te(1)
            ktids(in) = tube(itube)%ti(1,ion)

            in = idds(ir,1)
            knds (in) = tube(itube)%ne(2)      ! ni(1,ion) - 12/10/09 SL
            kvds (in) = tube(itube)%vi(2,ion)       
            kteds(in) = tube(itube)%te(2)
            ktids(in) = tube(itube)%ti(2,ion)

            ike = nks(ir)
          ELSE
            ike = nks(ir) - 1
          ENDIF

          knes (1:ike,ir) = fluid(cind1:cind2,ion)%ne
          knbs (1:ike,ir) = fluid(cind1:cind2,ion)%ne  ! *** USING ne NOT ni ***
          kvhs (1:ike,ir) = fluid(cind1:cind2,ion)%vi
          ktebs(1:ike,ir) = fluid(cind1:cind2,ion)%te
          ktibs(1:ike,ir) = fluid(cind1:cind2,ion)%ti

          osmion(1:ike,ir) = fluid(cind1:cind2,ion)%parion
          osmrec(1:ike,ir) = fluid(cind1:cind2,ion)%parrec
          osmcfp(1:ike,ir) = fluid(cind1:cind2,ion)%parano
          osmmp (1:ike,ir) = fluid(cind1:cind2,ion)%momsrc
          osmcfe(1:ike,ir) = fluid(cind1:cind2,ion)%eneano
          osmqe (1:ike,ir) = fluid(cind1:cind2,ion)%eneion

c...      Finish off core rings:
          IF (ir.LT.irsep) THEN
            knes (nks(ir),ir) = knes (1,ir)
            knbs (nks(ir),ir) = knbs (1,ir)
            kvhs (nks(ir),ir) = kvhs (1,ir) 
            ktebs(nks(ir),ir) = ktebs(1,ir)
            ktibs(nks(ir),ir) = ktibs(1,ir)
          ENDIF

        ENDDO
      ENDDO

c...  
c      nlpdato = 0
c      nlpdati = 0
c      DO ir = irsep, nrs
c        IF (idring(ir).EQ.BOUNDARY) CYCLE
c        nlpdato = nlpdato + 1
c        nlpdati = nlpdati + 1
c        in = idds(ir,2)
c        lpdato(nlpdato,1) = REAL(ir)
c        lpdato(nlpdato,2) = kteds(in)
c        lpdato(nlpdato,3) = ktids(in)
c        lpdato(nlpdato,4) = knds(in)
c        in = idds(ir,1)
c        lpdati(nlpdati,1) = REAL(ir)
c        lpdati(nlpdati,2) = kteds(in)
c        lpdati(nlpdati,3) = ktids(in)
c        lpdati(nlpdati,4) = knds(in)
c      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE SetupSOL28Options
      USE mod_sol28_global
      IMPLICIT none


      INCLUDE 'params'
      INCLUDE 'slcom'




      RETURN
 99   STOP
      END

c
c ======================================================================
c
      SUBROUTINE AssignSOL28Nodes_Old(itube,nnode,mnode,node)
      USE mod_sol28
      USE mod_sol28_global
      IMPLICIT none

      INTEGER itube,nnode,mnode       
      TYPE(type_node) :: node(*)

c      SUBROUTINE FindS28Parameters_V4(ir,te,ne,nf,s,new)

      INTEGER ir
      LOGICAL new
      REAL    te(0:6),ne(0:6),nf,s(0:6),isat(0:6)

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL    GetRelaxationFraction,GetJsat

      INTEGER i0,i1,i2,i3,ik,id,index,mode,ik1,ir1,id1,
     .        ikcell(3),ircell(3),id7
      LOGICAL tc,nc,density,tetarget,firsttime
      REAL    frac,t0,t1,n0,n1,A,B,C,tetmp1,tetmp2,coord,expon,
     .        psin0,psin1,prb1,tmp1,val,val0,val1,val2
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd
      DATA firsttime /.FALSE./

c...  
      ir = tube(itube)%ir
c      IF (output) WRITE(0,*) ' ASSIGN PARAMS: IR=',ir
     

      s (0) = 0.0
      ne(1) = 0.0
      te(1) = 0.0
      te(2) = 0.0
      te(3) = 0.0
      te(4) = 0.0
      te(5) = 0.0


      frac = GetRelaxationFraction()


      IF (new.AND.frac.NE.0.0) THEN
c        WRITE(0,*) 'FRAC:',frac,rel_step,rel_nstep
      ENDIF


      DO i1 = 2, osmns28
        i0 = i1 - 1
        IF ((osms28(i0,1).NE.osms28(i1,1)).OR.osms28(i1,2).EQ.0.0) CYCLE

c...    Do not apply data if IR is outside specified range:
        IF ((osms28(i1,11).NE.0.0.AND.REAL(ir).LT.osms28(i1,11) ).OR.
     .      (osms28(i1,12).NE.0.0.AND.REAL(ir).GT.osms28(i1,12)))
     .    CYCLE

c...    Check that rings from different grid regions are not in the same group
c       of rings:
        DO i2 = NINT(osms28(i1,11)), NINT(osms28(i1,12))-1
          IF (ringtype(i2).NE.ringtype(i2+1)) THEN
c            WRITE(0,*) 'RINGTYPE:',i1
c            WRITE(0,*) 'RINGTYPE:',i2,ringtype(i2)
c            WRITE(0,*) 'RINGTYPE:',i2+1,ringtype(i2+1)
c            WRITE(0,*) 'RINGTYPE:',osmns28
c            CALL ER('FindS28Parameters_V3','Mixed RINGTYPE',*99)
            IF (sloutput.AND.firsttime) THEN
              firsttime = .FALSE.
              WRITE(0,*)
              WRITE(0,*) '-----------------------------------------'
              WRITE(0,*) ' THAT FUNNY THING ABOUT MIXING REGIONS!? '
              WRITE(0,*) '-----------------------------------------'
              WRITE(0,*)
            ENDIF
          ENDIF      
        ENDDO


        IF (.FALSE..AND.osms28(i1,2).EQ.99.0) THEN
        ELSE

c MODE                          P1           P2
c   1 - power law               coordinate   index
c   2 - exponential v0-v2       coordinate   index
c   3 - exponential to infinity
c   4 - from probe data         coordinate   probe number
c
c
c   coord = 1 - linear on line segment
c         = 2 - linear on line segment, but from first to last ring intersection
c         = 3 - PSIn 
c         = 4 - RHO
c

          index = NINT(osms28(i1,1)) - 1
          mode  = NINT(osms28(i1,2))
          coord = osms28(i1,3)
          expon = osms28(i1,4)

c...  Decide if specified upstream data is density or pressure:
          density = .TRUE.   ! *** I DON'T LIKE THIS SOLUTION ***
          IF (index.EQ.3.AND.
     .        (ringtype(ir).EQ.SOL1.AND.s28nemode   .EQ.1.OR.
     .         ringtype(ir).EQ.PFZ .AND.s28nemodepfz.EQ.1)) 
     .      density = .FALSE.

          a1 = DBLE(osms28(i1-1,5))
          a2 = DBLE(osms28(i1-1,6))
          b1 = DBLE(osms28(i1  ,5))
          b2 = DBLE(osms28(i1  ,6))

          DO ik = 1, nks(ir)

            id = korpg(ik,ir)
            c1 = 0.5D0 * DBLE(rvertp(1,id) + rvertp(2,id))
            c2 = 0.5D0 * DBLE(zvertp(1,id) + zvertp(2,id))
            d1 = 0.5D0 * DBLE(rvertp(3,id) + rvertp(4,id))
            d2 = 0.5D0 * DBLE(zvertp(3,id) + zvertp(4,id))
       
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
 
            IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .          tcd.GE.0.0.AND.tcd.LT.1.0) THEN
c...          Intersecting between the line segment and the ring is found:

c...          These are here in case psin0 gets picked up when linking exponential
c             decay data to neighbouring ring:
              IF (ringtype(ir).EQ.PFZ) THEN
                psin0 = -HI
                psin1 =  HI
              ELSE
                psin0 = HI
                psin1 = HI
              ENDIF

              IF (.FALSE.) THEN
c...            Need a check to know if the ring is over-constrained:
              ELSEIF (index.GE.0.AND.index.LE.6) THEN
c              ELSEIF (index.GE.1.AND.index.LE.5) THEN

                IF (index.GE.1.AND.index.LE.5) THEN
                  s(index)=ksb(ik-1,ir)+
     .                     SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
                ENDIF

c...            Find data boundary values -- NEEDS WORK!:
                i2 = i0

                i3 = i1

c...            Flag...
                tetarget = .FALSE.
                IF ((index.EQ.1.OR.index.EQ.5).AND.osms28(i3,7).LT.0.0) 
     .            tetarget = .TRUE.

                IF     (mode.EQ.1.OR.mode.EQ.2.OR.mode.EQ.3) THEN
c...              Interpolation boundary values provided in the input file:
c                    WRITE(0,*) 'DATA1:',index,osms28(i2,7),osms28(i2,8)

                  t0 = osms28(i2,7) + frac*(osms28(i2,9) -osms28(i2,7))
                  n0 = osms28(i2,8) + frac*(osms28(i2,10)-osms28(i2,8))
                  t1 = osms28(i3,7) + frac*(osms28(i3,9) -osms28(i3,7))
                  n1 = osms28(i3,8) + frac*(osms28(i3,10)-osms28(i3,8))

                  IF ((index.EQ.3.OR.index.EQ.4.OR.index.EQ.5).AND.
     .                (osms28(i2,7).EQ.-99.0.OR.
     .                 osms28(i2,8).EQ.-99.0)) THEN
c...                Linking to another plasma region that is *already* calculated: 
                    CALL FindCell(i2,i3,ir,ikcell,ircell)
                    IF (ringtype(ir).EQ.PFZ) THEN
                      ik1 = ikouts(ikcell(2),ircell(2))
                      ir1 = irouts(ikcell(2),ircell(2))
c                      WRITE(0,*) 'FOUND YOU PFZ:',ik1,ir1,
c     .  density,tetarget,index
                    ELSE
                      ik1 = ikins(ikcell(2),ircell(2))
                      ir1 = irins(ikcell(2),ircell(2))
                    ENDIF
                    IF (osms28(i2,7).EQ.-99.0) t0 = ktebs(ik1,ir1)
                    IF (osms28(i2,8).EQ.-99.0) THEN
                      IF (density) THEN
                        n0 = knbs (ik1,ir1)
c                        WRITE(0,*) 'denisity POWER:',n0
                      ELSE
                        n0 = 2.0 * ktebs(ik1,ir1) * knbs(ik1,ir1)  ! Add Ti and M? 
                      ENDIF
                    ENDIF
c...                Shouldn't really be outer target (all this would go away if PSITARG was
c                   assigned properly):
                    IF (coord.EQ.3) psin0 = psitarg(ir1,1)
c                   IF (coord.EQ.4) rho0 =
c                    WRITE(0,*) 'DATA2:',t0,psin0,ik1,ir1,index
c                    WRITE(0,*) 'DATA2:',n0,psin0,ik1,ir1
c                    STOP 'rgsd'
                  ENDIF

c...              Make sure that t0,1 are positive:
                  IF (tetarget) THEN
                    t0 = ABS(t0)
                    t1 = ABS(t1)
                  ENDIF

                  IF (coord.EQ.1) THEN
c...                Linear along the line segment, nice and simple:
                    val0 = 0.0
                    val1 = 1.0
                    val = SNGL(tab)
c                    val = SNGL(tcd) ! BUG
                  ELSE
c...              
                    IF (NINT(osms28(i2,11)).NE.NINT(osms28(i3,11)).OR.
     .                  NINT(osms28(i2,12)).NE.NINT(osms28(i3,12)).OR.
     .                  NINT(osms28(i2,11)).GT.NINT(osms28(i2,12)).OR.
c...                    This is here because PSITARG is currently not assigned in the core... 
     .                  NINT(osms28(i2,11)).LT.irsep)
     .                CALL ER('FindS28Parameters_V3','Invalid '//
     .                        'IR range',*99)

c...                Need range of PSIn over the segment:
                    DO ir1 = NINT(osms28(i2,11)), NINT(osms28(i2,12)) 
                      DO ik1 = 1, nks(ir1)
                        id1 = korpg(ik1,ir1)
                        c1 = 0.5D0 * DBLE(rvertp(1,id1) + rvertp(2,id1))
                        c2 = 0.5D0 * DBLE(zvertp(1,id1) + zvertp(2,id1))
                        d1 = 0.5D0 * DBLE(rvertp(3,id1) + rvertp(4,id1))
                        d2 = 0.5D0 * DBLE(zvertp(3,id1) + zvertp(4,id1))
                        CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
                        IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .                      tcd.GE.0.0.AND.tcd.LT.1.0) THEN
c...                      Should be more careful which PSITARG is used and when...
c                          WRITE(0,*) 'MARK:',ir1,psin0,psin1
                          IF     (ringtype(ir).EQ.SOL1) THEN
                            IF (psin0.EQ.HI) psin0 = psitarg(ir1,1)
                            IF (psin0.NE.HI) psin1 = psitarg(ir1,1)
c                            WRITE(0,*) 'PSIN:',ir1,psin0,psin1
                          ELSEIF (ringtype(ir).EQ.PFZ) THEN
                            psin0 = MAX(psin0,psitarg(ir1,1))
                            psin1 = MIN(psin1,psitarg(ir1,1))
                          ELSE
                            CALL ER('FindS28Parameters_V3','Invalid '//
     .                              'RINGTYPE',*99)
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDDO


c                    WRITE(0,*) 'PSIn:',ir,psin0,psin1

                  
                    IF     (coord.EQ.2) THEN
c...                  Spatial along the IR-range of applicability:
                      STOP 'NOT READY 2'
                    ELSEIF (coord.EQ.3) THEN
c...                  PSIn:
                      val0 = 0.0
                      val1 = ABS(psin1         - psin0)
                      val  = ABS(psitarg(ir,1) - psin0)
                    ELSEIF (coord.EQ.4) THEN
c...                  RHO:
                      STOP 'NOT READY 4'
                    ELSE
                      WRITE(0,*) 'I0,I1=',i0,i1
                      CALL ER('S28params_v3','Invalid COORD A',*99)
                    ENDIF
                  ENDIF

                ELSEIF (mode.EQ.4) THEN
c...              Load probe data, dummy values here:
                  t0 = osms28(i2,7)
                  t1 = t0
                  n0 = osms28(i2,8)
                  n1 = n0

                ELSE
                  CALL ER('S28params_v3','Invalid MODE',*99)   
                ENDIF

c...            Check if quantities should be assigned:
                tc = .TRUE.
                nc = .TRUE.
                IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
                IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

c                IF (ir.GT.22.AND.ir.LT.irwall.AND.index.EQ.1.AND.
c     .              .NOT..FALSE..AND.(.TRUE..OR.index.EQ.4)) THEN
c                  WRITE(0,*) 'PSIn:',ir,psitarg(ir,1),psin0,psin1
c                  WRITE(0,*) 'i2,3:',i2,i3
c                  WRITE(0,*) 'val :',ir,val,val0,val1
c                  WRITE(0,*) 'Te  :',ir,t0,t1
c                  WRITE(0,*) 'ne  :',ir,n0,n1
c                ENDIF

                IF     (mode.EQ.1) THEN
c...              Power law between v1 and v2:
                  val2 = (val - val0) / (val1 - val0)
c                  IF (index.EQ.5) WRITE(0,*) '5:',val2,ir

                  IF (tc) te(index) = t0 + val2**expon * (t1 - t0)
                  IF (nc) ne(index) = n0 + val2**expon * (n1 - n0)

c                  IF (index.EQ.3.AND.ir.EQ.71) THEN
c                    WRITE(0,*) 'A:',val,n0,n1,density
c                  ENDIF
      

                ELSEIF (mode.EQ.2) THEN
c...              Exponential decay between v1 and v2:
                  C = expon  ! BUG!
c                  C = -expon
c                  val = 1.0

                  IF (tc) THEN
                    A = (t1 - t0) / (EXP(-val1 / C) - 1.0)
                    B = t0 - A
                    te(index) = A * EXP(-val / C) + B
                  ENDIF
                  IF (nc) THEN
                    A = (n1 - n0) / (EXP(-val1 / C) - 1.0)
                    B = n0 - A
                    ne(index) = A * EXP(-val / C) + B
                  ENDIF



                ELSEIF (mode.EQ.3) THEN
c...              Exponential decay to infinity:
                  C = expon
                  A = t0 - t1
                  B = t1 
                  IF (tc) te(index) = A * EXP(-val / C) + B
                  A = n0 - n1
                  B = n1
                  IF (nc) ne(index) = A * EXP(-val / C) + B

                ELSEIF (mode.EQ.4) THEN
c...              Load probe data from the .experiments file:
                  IF     (coord.EQ.3) THEN
                    prb1 = -1.0
                  ELSEIF (coord.EQ.4) THEN
                    prb1 = -3.0
                  ELSE
                      WRITE(0,*) 'I0,I1=',i0,i1,coord
                    CALL ER('S28params_v3','Invalid COORD B',*99)   
                  ENDIF

                  IF (tc) THEN
                    tmp1 = prb1
                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,7)),2,tmp1) 
                    te(index) = tmp1
                    WRITE(0,*) 'TE PROBE:',tmp1,ir,psitarg(ir,1)
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY'
                    ELSEIF (expon.EQ.2.0) THEN
c...                  
                      te(0) = -te(index)
                      te(6) = -te(index)
                    ELSEIF (expon.EQ.3.0) THEN
c...                  
                      te(0) = 97.0
                      te(6) = 97.0
                    ENDIF
                  ENDIF

                  IF (nc) THEN
                    tmp1 = prb1
                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,8)),1,tmp1) 
                    ne(index) = tmp1
                    WRITE(0,*) 'NE PROBE:',tmp1,ir,psitarg(ir,1)
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY'
                    ELSEIF (expon.EQ.2.0) THEN
c...                  
                      ne(0) = -ne(index)
                      ne(6) = -ne(index)
                    ELSEIF (expon.EQ.3.0) THEN
                    ENDIF
                  ENDIF

c                  WRITE(0,*) 'PROBIN:',prb1,te(index),ne(index),tc,nc

                ELSE
                  CALL ER('S28params_v3','Invalid MODE',*99)   
                ENDIF

              ELSE
                CALL ER('FindSOL28Parameters','Invalid parameter '//
     .                  'index',*99)
              ENDIF


              IF (tetarget) THEN
c                WRITE(0,*) 'TETARGET:',tetarget,te(index),i3,index
                te(index) = -te(index)
              ENDIF

              IF (.NOT.density.AND.index.EQ.3.AND.nc.AND.mode.NE.4) THEN
c...            Convert the pressure value to density:
                IF (te(3).EQ.0.0) THEN
                  CALL ER('FindSOL28Parameters','Te3 not assigned',*99)
                ELSE
c...              Assumes that the Mach no. is low (note: Ti.NE.Te is not a problem
c                 since that situation is currently resolved in the SOL28 routine):
c                  WRITE(0,*) 'CONVERTING PRESSURE TO DENSITY',ne(3)
c                  WRITE(0,*) 'OLD D:',ne(3)
                  ne(3) = ne(3) / (2.0 * te(3))
c                  WRITE(0,*) 'NEW D:',ne(3)
                ENDIF
              ENDIF

            ENDIF

          ENDDO

c          IF (.NOT.density.AND.index.EQ.3.AND.nc.AND.mode.NE.4) THEN
cc...        Convert the pressure value to density:
c            IF (te(3).EQ.0.0) THEN
c              CALL ER('FindSOL28Parameters','Te(3) not assigned',*99)
c            ELSE
cc...          Assumes that the Mach no. is low (note: Ti.NE.Te is not a problem
cc             since that situation is currently resolved in the SOL28 routine):
cc              WRITE(0,*) 'CONVERTING PRESSURE TO DENSITY',ne(3)
c              WRITE(0,*) 'OLD D:',ne(3)
c              ne(3) = ne(3) / (2.0 * te(3))
c              WRITE(0,*) 'NEW D:',ne(3)
c            ENDIF
c          ENDIF

c...    End of OSMS28(I1,1).EQ.98.0 block:      
        ENDIF

      ENDDO




c          WRITE(0,*) 'NE3=',ne(3)
      IF (.TRUE.) THEN
c...    Specify te(1) from target data.  The -ve is to trigger
c       te(0)=te(1):
        IF (te(1).EQ.-98.0) te(1) = -kteds(idds(ir,2))

c...    Specify te(5) from target data.  The -ve is to trigger
c       te(6)=te(5):
        IF (te(5).EQ.-98.0) te(5) = -kteds(idds(ir,1))
      ENDIF 


     

      IF (.NOT..TRUE.) THEN


        IF (te(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
          te(3) = tetmp1
c...      Don't necessarily want both targets done...:
          te(0) = 97.0
          te(6) = 97.0
        ENDIF
        IF (ne(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
          ne(3) = tetmp1
        ENDIF


        IF (te(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          te(3) = tetmp1
          te(0) = -te(3)
          te(6) = -te(3)
        ENDIF
        IF (ne(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          ne(3) = tetmp1
          ne(0) = -ne(3)
          ne(6) = -ne(3)
        ENDIF

c        WRITE(0,*) 'NE(3):',ne(3),te(3)

        IF (te(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,2,tetmp1) 
          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c         te(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
        ENDIF
        IF (ne(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,1,tetmp1) 
          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          WRITE(PINOUT,*) 'COMPROMISE:',ir,psitarg(ir,2),ne(3),te(3)
        ENDIF

        IF (te(3).LE.-1.0) CALL LoadProbeDataS28(ir,s28probe,2,te(3))
c *HACKISH*
        IF (ne(3).LE.-1.0.AND.ne(3).GT.-100.0) 
     .    CALL LoadProbeDataS28(ir,s28probe,1,ne(3))

        IF (ne(3).LT.-100.0) THEN
c...      Pressure specified, not density.  Calculate ne(3) from the pressure
c         and the specified te(3) value, assuming Ti=Te:
          ne(3) = -ne(3) / (2.0 * te(3))
        ENDIF

c...    Specify te(1) from target data.  The -ve is to trigger
c       te(0)=te(1):
        IF (te(1).EQ.-99.0) te(1) = -kteds(idds(ir,2))

c...    Specify te(5) from target data.  The -ve is to trigger
c       te(6)=te(5):
        IF (te(5).EQ.-99.0) te(5) = -kteds(idds(ir,1))

      ELSE
      ENDIF


      IF (.TRUE.) THEN

        s(0) = 0.0
        s(6) = ksmaxs(ir)

        DO i1 = 0, 6
          DO ik = 1, nks(ir)
            IF (s(i1).GE.ksb(ik-1,ir).AND.s(i1).LE.ksb(ik,ir)) THEN
              node(i1+1)%icell = ik
              IF (i1.EQ.3) s(3) = kss(ik,ir)
            ENDIF
          ENDDO
        ENDDO

        nnode = 7
        mnode = 4

c...    Assign values to nodes:
        node(1:7)%s  = s (0:6)
        node(1:7)%ne = ne(0:6)
        node(1:7)%te = te(0:6)
c...    Assign other quantites:
        node(1:7)%jsat(1) = 0.0
        node(1:7)%pe      = 0.0
        node(1:7)%ni(1)   = 0.0
        node(1:7)%pi(1)   = 0.0
        node(1:7)%ti(1)   = 0.0
        node(1:7)%machno  = 0.0
        node(1:7)%epot    = 0.0
        node(1:7)%efield  = 0.0

        id1 = idds(ir,2)
        id7 = idds(ir,1)
        IF (node(1)%ne.EQ.0.0) node(1)%ne = knds (id1)
        IF (node(7)%ne.EQ.0.0) node(7)%ne = knds (id7)
        IF (node(1)%te.EQ.0.0) node(1)%te = kteds(id1)
        IF (node(7)%te.EQ.0.0) node(7)%te = kteds(id7)
        IF (node(1)%ti(1).EQ.0.0) node(1)%ti(1) = ktids(id1)
        IF (node(7)%ti(1).EQ.0.0) node(7)%ti(1) = ktids(id7)

c REAL FUNCTION GetJsat(te,ti,ne,v)
        id = id1
        node(1)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(1)%jsat(1) = -ABS(node(1)%jsat(1))
        id = id7
        node(7)%jsat(1)= GetJsat(kteds(id),ktids(id),knds(id),kvds(id))
        node(7)%jsat(1) = -ABS(node(7)%jsat(1))


      ENDIF


      RETURN
99    CONTINUE
      WRITE(0,*) 'IK,IR=',ik,ir,i1,osms28(i1,1)
      WRITE(0,*) 'TAB,TCD=',tab,tcd
      STOP
      END
c
c ======================================================================
c
      SUBROUTINE SaveLegacyData(fname)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'

      CHARACTER*(*) fname
      INTEGER fp,i1,i2

      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)            
      WRITE(fp,ERR=98) 1.00 
      WRITE(fp,ERR=98) MAXNRS
      WRITE(fp,ERR=98) (grdntreg(i1),grdntseg(1:MAXNRS,i1),
     .                  (grdtseg(1:MAXNRS,i2,i1),i2=1,MAXNRS),i1=1,2)
      CLOSE (fp)
      
      RETURN
 98   CALL ER('SaveGrid','Problem saving legacy data file',*99)
 99   STOP
      END
c
c ======================================================================
c
c      SUBROUTINE SetTargetConditions_Legacy
c...  Dummy
c      END
c
c ======================================================================
c
      SUBROUTINE BuildLinearGrid_OLD
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'comtor'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'


      INTEGER id,ik,ir,i1,grid_option,nrings_inner,nrings_outer
      REAL*8  r,delr,L,r1,r2,z1,z2,frac1,frac2,
     .        vessel_radius,brat,frac,r_inner,r_outer,delta


      grid_option = 3  ! 7  ! 6

      brat = 1.0

      r0 = 0.0001D0  ! Need this tiny displacement to keep EIRENE04 from falling over 
c      r0 = 0.0000001D0  ! Need this tiny displacement to keep EIRENE04 from falling over 

      SELECTCASE (grid_option)
        CASE (1)  ! Full vessel, mirrored
          vessel_radius = 0.02D0
          L = 3.6D0                   ! Total length of mirrored plasma column (m)
          r = 0.015D0                 ! Plasma radius (m)
          z0 = L / 2.0D0              ! Height of the centre of the plasma (m)
          delr = (vessel_radius - r)  ! Distance from plasma to outer wall (m)
          maxrings = 10               ! Number of flux tubes (if changed, also need to change triangle grid in input file)
          nks(1:maxrings) = 100       ! Number of cells on each tube
        CASE (2)  ! Full vessel
          vessel_radius = 0.02D0
          L = 1.8D0
          r = 0.015D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 10      
          nks(1:maxrings) = 100         
        CASE (3) ! Target chamber
          brat = 0.05 ! 0.985 ! 0.5
  
          vessel_radius = 0.02D0 ! 0.05D0 ! 0.02D0
          L = 0.55D0  ! 0.56D0
          r = 0.015D0 ! 0.03D0  ! 0.015D0
          z0 = L / 2.0D0 ! 0.0 ! L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 10      
          nks(1:maxrings) = 20  ! 50  ! 175
        CASE (4) ! Target chamber: fancy
          vessel_radius = 0.16D0 
          L = 0.55D0
          r = 0.08D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (5) ! Target chamber: fancy #2, full vessel and small volume 
          vessel_radius = 0.16D0 
          L = 0.55D0
          r = 0.03D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (6) ! Target chamber: fancy #3, small volume 
          vessel_radius = 0.05D0 
          L = 0.55D0
          r = 0.03D0          
          z0 = L / 2.0D0              
          delr = (vessel_radius - r)  
          maxrings = 50      
          nks(1:maxrings) = 150
        CASE (7) ! Grid with 2 radial regions
          vessel_radius = 0.05D0 
          L = 0.55D0
          z0 = L / 2.0D0              
          nrings_inner = 30
          nrings_outer = 20
          maxrings = nrings_inner + nrings_outer
          r_inner = 0.01D0
          r_outer = 0.03D0     
          delr = (vessel_radius - r_outer)  
          nks(1:maxrings) = 150
      ENDSELECT

      id = 0

      DO ir = 1, maxrings

        IF (.TRUE.) THEN
          IF     (grid_option.EQ.4) THEN
            frac1 = (DBLE(ir-1) / DBLE(maxrings))**0.7
            frac2 = (DBLE(ir  ) / DBLE(maxrings))**0.7
            r1 = frac1 * r
            r2 = frac2 * r
          ELSEIF (grid_option.EQ.7) THEN
            IF (ir.LE.nrings_inner) THEN
              frac  = DBLE(ir-1) / DBLE(nrings_inner)
              delta = r_inner / DBLE(nrings_inner)
              r1 = frac * r_inner
              r2 = r1 + delta
            ELSE
              frac  = DBLE(ir-nrings_inner-1) / DBLE(nrings_outer)
              delta = (r_outer - r_inner) / DBLE(nrings_outer)
              r1 = r_inner + frac * (r_outer - r_inner)
              r2 = r1 + delta
            ENDIF
          ELSE
            frac  = DBLE(ir-1) / DBLE(maxrings)
            delta = r / DBLE(maxrings)
            r1 = frac * r 
            r2 = r1 + delta
          ENDIF
          IF (ir.EQ.1) r1 = r1 + r0
        ENDIF

C       Krieger IPP/07 - SUN compiler does not know SNGL, replaced by REAL  -strange since SNGL is used elsewhere... -SL
C       psitarg(ir,1) = ABS(0.5*(SNGL(r1+r2)))
C       psitarg(ir,2) = ABS(0.5*(SNGL(r1+r2)))
        psitarg(ir,1) = ABS(0.5*(REAL(r1+r2)))
        psitarg(ir,2) = ABS(0.5*(REAL(r1+r2)))
        idring(ir) = TARTOTAR

c        WRITE(0,*) 'IR:',ir,psitarg(ir,1)

c       nks(ir) = 100

        DO ik = 1, nks(ir)

          SELECTCASE (grid_option)
            CASE (1)  ! Full vessel, mirrored
              IF (.TRUE.) THEN
                frac = DBLE(ik-1) / DBLE(nks(ir)) 
                delta = L / DBLE(nks(ir)) 
                z1 = (0.5 - frac) * L
                z2 = z1 - delta
              ENDIF
            CASE (2)  ! Full vessel
              IF (.TRUE.) THEN
                frac = DBLE(ik-1) / DBLE(nks(ir)) 
                delta = L / DBLE(nks(ir)) 
                z1 = (1.0 - frac) * L
                z2 = z1 - delta       
              ENDIF
            CASE (3)  ! Target chamber
              IF (.TRUE.) THEN
                frac = DBLE(ik-1) / DBLE(nks(ir)) 
                delta = L / DBLE(nks(ir)) 
                z1 = (0.5 - frac) * L + z0
c                z1 = (1.0 - frac) * L 
                z2 = z1 - delta       
              ENDIF
            CASE (4:7) ! Target chamber: fancy
              frac1 = DBLE(ik-1) / DBLE(nks(ir)) 
              frac2 = DBLE(ik  ) / DBLE(nks(ir)) 
              frac1=SIGN(.5D0,frac1-.5D0)*(ABS(frac1-0.5)/0.5)**1.00+0.5
              frac2=SIGN(.5D0,frac2-.5D0)*(ABS(frac2-0.5)/0.5)**1.00+0.5
              z1 = (1.0 - frac1) * L
              z2 = (1.0 - frac2) * L     
          ENDSELECT

c          frac = ((ABS(0.5 * (z1 + z2) - z0) + 0.001) / L * 2.0)**0.05
c          IF (ir.EQ.2) WRITE(0,*) frac
c          bratio(ik,ir) = SNGL(brat * frac)
          bratio(ik,ir) = SNGL(brat)

          kbfs  (ik,ir) = 1.0 / brat
          bts   (ik,ir) = cbphi 

          id = id + 1

          korpg(ik,ir) = id

          nvertp(id) = 4

!          frac = 1.0D0 + 1.0D0 * DBLE(ik-1) / DBLE(nks(ir) - 1)
          frac = 1.0D0

          IF (ik.EQ.1) THEN
            rvertp(1,id) = SNGL(r1)
            rvertp(2,id) = SNGL(r2)
            zvertp(1,id) = SNGL(z1)
            zvertp(2,id) = SNGL(z1)
          ELSE
            rvertp(1,id) = rvertp(4,id-1)
            rvertp(2,id) = rvertp(3,id-1)
            zvertp(1,id) = zvertp(4,id-1)
            zvertp(2,id) = zvertp(3,id-1)
          ENDIF

          rvertp(3,id) = SNGL(r2 * frac)
          rvertp(4,id) = SNGL(r1 * frac)
          zvertp(3,id) = SNGL(z2)
          zvertp(4,id) = SNGL(z2)

          rs(ik,ir) = 0.0
          zs(ik,ir) = 0.0
          DO i1 = 1, nvertp(id)
            rs(ik,ir) = rs(ik,ir) + rvertp(1,id)
            zs(ik,ir) = zs(ik,ir) + zvertp(1,id)
          ENDDO
          rs(ik,ir) = rs(ik,ir) / REAL(nvertp(id))
          zs(ik,ir) = zs(ik,ir) / REAL(nvertp(id))

        ENDDO

      ENDDO

      npolyp  = id
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

      ikto = 2
      ikti = 3

      rves = 0
      rvesm = 0

      irsep  = 1
      irwall = maxrings 
      irtrap = irwall
      nrs    = irwall
      nbr    = 0

      WRITE(0,*) 'NVERT:',nvertp(5)

      CALL InsertRing(1         ,BEFORE,PERMANENT)
      CALL InsertRing(maxrings+1,AFTER ,PERMANENT)

      WRITE(0,*) 'NVERT:',nvertp(5)

c...  Necessary..? 
      cutring = 1
      cutpt1 = ikto
      cutpt2 = ikti

      idring(1) = -1
      idring(nrs) = -1

c...  Modify the grid based on entries in the GRDMOD array assigned 
c     from the input file:
c      IF (grdnmod.NE.0) CALL TailorGrid

      rmin = HI
      rmax = LO
      zmin = HI
      zmax = LO
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          rmin = MIN(rmin,rs(ik,ir))
          rmax = MAX(rmax,rs(ik,ir))
          zmin = MIN(zmin,zs(ik,ir))
          zmax = MAX(zmax,zs(ik,ir))
        ENDDO
      ENDDO

      rxp =  0.0 
      zxp =  0.25 * zmin + 0.75 * zmax
c      zxp = -0.25 * L

c...  Neutral wall
      IF (cneur.EQ.4) THEN

         SELECTCASE (grid_option)
           CASE (1)  ! Full vessel, mirrored
             nves = 20
             ir = irwall-1
             r1 = DBLE(rvertp(2,korpg(1      ,ir)))
             r2 = r1 + delr
             z1 = DBLE(zvertp(2,korpg(1      ,ir)))
             z2 = DBLE(zvertp(3,korpg(nks(ir),ir)))
  
             rves(1)  =  SNGL(r1)
             zves(1)  =  SNGL(z1)
             rves(2)  =  SNGL(r2)
             zves(2)  =  SNGL(z1)
  
             rves(3)  =  SNGL(r2)
             zves(3)  =  SNGL(L) / 2.0 - 0.56
             rves(4)  =  SNGL(r1) + 0.0101 
             zves(4)  =  SNGL(L) / 2.0 - 0.56
             rves(5)  =  SNGL(r1) + 0.0001
             zves(5)  =  SNGL(L) / 2.0 - 0.57
             rves(6)  =  SNGL(r2)
             zves(6)  =  SNGL(L) / 2.0 - 0.57
  
             rves(7)  =  SNGL(r2)
             zves(7)  =  0.03
             rves(8)  =  SNGL(r1) + 0.0001
             zves(8)  =  0.03
             rves(9)  =  SNGL(r1) + 0.0001
             zves(9)  =  0.02
             rves(10) =  SNGL(r2)
             zves(10) =  0.02
  
             rves(11) =  r2
             zves(11) = -0.02
             rves(12) =  r1 + 0.0001
             zves(12) = -0.02
             rves(13) =  r1 + 0.0001
             zves(13) = -0.03
             rves(14) =  r2
             zves(14) = -0.03
  
             rves(15) =  r2
             zves(15) = -L / 2.0 + 0.56
             rves(16) =  r1 + 0.0001
             zves(16) = -L / 2.0 + 0.56
             rves(17) =  r1 + 0.0101
             zves(17) = -L / 2.0 + 0.57
             rves(18) =  r2
             zves(18) = -L / 2.0 + 0.57
  
             rves(19) =  r2
             zves(19) =  z2
             rves(20) =  r1
             zves(20) =  z2
           CASE (2)  ! Full vessel
             nves = 12
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1)  =  r1
             zves(1)  =  z1
             rves(2)  =  r2
             zves(2)  =  z1
  
             rves(3)  =  r2
             zves(3)  =  L - 0.56
             rves(4)  =  r1 + 0.0101 
             zves(4)  =  L - 0.56
             rves(5)  =  r1 + 0.0001 
             zves(5)  =  L - 0.57
             rves(6)  =  r2
             zves(6)  =  L - 0.57
  
             rves(7)  =  r2
             zves(7)  =  0.03
             rves(8)  =  r1 + 0.0001
             zves(8)  =  0.03
             rves(9)  =  r1 + 0.0001
             zves(9)  =  0.02
             rves(10) =  r2
             zves(10) =  0.02
  
             rves(11) =  r2
             zves(11) =  z2
             rves(12) =  r1
             zves(12) =  z2
           CASE (3)  ! Target chamber
             nves = 7
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir)) - 0.0001 ! So that the clipping code is required / activated
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir)) 
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.55 * z1 + 0.45 * z2
             rves(4) =  r2
             zves(4) =  0.50 * z1 + 0.50 * z2
             rves(5) =  r2
             zves(5) =  0.45 * z1 + 0.55 * z2
  
             rves(6) =  r2
             zves(6) =  z2 
             rves(7) =  rvertp(3,korpg(nks(ir),ir)) - 0.0001 ! r1
             zves(7) =  z2
           CASE (4:5)  ! Target chamber: fancy
             nves = 10
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.11

             rves(4) =  0.21
             zves(4) =  0.11

             rves(5) =  0.21
             zves(5) =  0.06

             rves(6) =  0.212
             zves(6) =  0.06

             rves(7) =  0.212
             zves(7) =  0.05

             rves(8) =  0.21
             zves(8) =  0.05

             rves(9) =  0.21
             zves(9) =  z2
  
             rves(10) =  r1
             zves(10) =  z2
           CASE (6:7)  ! Target chamber: fancy #3, small volume
             nves = 10
             ir = irwall-1
             r1 = rvertp(2,korpg(1      ,ir))
             r2 = r1 + delr
             z1 = zvertp(2,korpg(1      ,ir))
             z2 = zvertp(3,korpg(nks(ir),ir))
  
             rves(1) =  r1
             zves(1) =  z1
             rves(2) =  r2
             zves(2) =  z1

             rves(3) =  r2
             zves(3) =  0.11

             rves(4) =  r2 + 0.01
             zves(4) =  0.11

             rves(5) =  r2 + 0.01
             zves(5) =  0.06

             rves(6) =  r2 + 0.015
             zves(6) =  0.06

             rves(7) =  r2 + 0.015
             zves(7) =  0.05

             rves(8) =  r2 + 0.01
             zves(8) =  0.05

             rves(9) =  r2 + 0.01
             zves(9) =  z2
  
             rves(10) =  r1
             zves(10) =  z2
        ENDSELECT
      ENDIF


      IF (.TRUE.) THEN
        nvesm = nves - 1
        DO i1 = 1, nves-1
          rvesm(i1,1) = rves(i1)
          zvesm(i1,1) = zves(i1)
          rvesm(i1,2) = rves(i1+1)
          zvesm(i1,2) = zves(i1+1)
        ENDDO
      ENDIF
 
c      CALL DumpGrid('BUILDING LINEAR GRID')

      IF (grdnmod.GT.0) CALL TailorGrid


      CALL OutputData(85,'Linear')


c...  Add virtual boundary cells, which will be stripped off later:
      IF (CTARGOPT.EQ.0.OR.CTARGOPT.EQ.1.OR.CTARGOPT.EQ.2.OR.
     .    CTARGOPT.EQ.3.OR.CTARGOPT.EQ.6) 
     .   CALL AddPoloidalBoundaryCells

c      STOP 'WHA-WHO!'

      RETURN
 99   STOP
      END
c
c
c
c
c
c
c
c
c ========================================================================
c
      SUBROUTINE FindKnot_SL(NUMZONE,izone,condition,
     .                       index1,index2)
      USE mod_grid
      IMPLICIT none
      INCLUDE 'params'  ! for SLOUTPUT
c
c     jdemod
c
c     Condition=1 finds two different knots which have the same 
c                 R,Z values for knot number 1. This condition
c                 should only occur at an Xpoint
c
c     Condition=2 finds knot 1=2 and knot 4=3 for test cell vs. other cells
c
c     Condition=3 finds knot 3=2 and knot 4=1 for test cell vs. other cells
c
c     Condition=4 finds knot 2=1 and knot 3=4 for test cell vs. other cells
c
c     Condition=5 finds knot 1=4 and knot 2=3 for test cell vs. other cells
c
c      TYPE type_grid_cell
c        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
c        REAL*8  :: rcen,zcen,bratio,rv(4),zv(4)
c      ENDTYPE type_grid_cell
c      TYPE type_cell_old
c        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
c        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
c      ENDTYPE type_cell_old

      INTEGER index1,index2,NUMZONE,izone(NUMZONE+1,NUMZONE),
     .        condition
c      INTEGER nknot,index1,index2,NUMZONE,izone(NUMZONE+1,NUMZONE),
c      TYPE(type_grid_cell) :: knot(0:nknot)      

      REAL*8, PARAMETER :: DTOL=1.0D-06

      INTEGER i1,i2,r1,z1

      i1 = index1

      index2 = -1

      DO z1 = knot(i1)%zzone-1, knot(i1)%zzone+1
        IF (z1.LT.1.OR.z1.GT.NUMZONE) CYCLE
        DO r1 = knot(i1)%rzone-1, knot(i1)%rzone+1
          IF (r1.LT.1.OR.r1.GT.NUMZONE) CYCLE  
c          DO i2 = 1, nknot
          DO i2 = izone(r1,z1), izone(r1+1,z1)-1
            IF (i1.EQ.i2) CYCLE
c...
            IF     (condition.EQ.1.AND.
     .              DABS(knot(i1)%rv(1)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(1)-knot(i2)%zv(1)).LT.DTOL) THEN
c              IF (sloutput) WRITE(0,*) 'XPOINT:',i1,i2
              index2 = i2
              RETURN
            ELSEIF (condition.EQ.2.AND.
     .              DABS(knot(i1)%rv(1)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(1)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(4)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(4)-knot(i2)%zv(3)).LT.DTOL) THEN
c              WRITE(0,*) 'SIDE14:',i1,i2
              index2 = i2
              RETURN
            ELSEIF (condition.EQ.3.AND.
     .              DABS(knot(i1)%rv(3)-knot(i2)%rv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(3)-knot(i2)%zv(2)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(4)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(4)-knot(i2)%zv(1)).LT.DTOL) THEN
c              WRITE(0,*) 'SIDE14:',i1,i2
              index2 = i2
              RETURN
            ELSEIF (condition.EQ.4.AND.
     .              DABS(knot(i1)%rv(2)-knot(i2)%rv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(2)-knot(i2)%zv(1)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(3)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(3)-knot(i2)%zv(4)).LT.DTOL) THEN
c              WRITE(0,*) 'SIDE14:',i1,i2
              index2 = i2
              RETURN
            ELSEIF (condition.EQ.5.AND.
     .              DABS(knot(i1)%rv(1)-knot(i2)%rv(4)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(1)-knot(i2)%zv(4)).LT.DTOL.AND.
     .              DABS(knot(i1)%rv(2)-knot(i2)%rv(3)).LT.DTOL.AND.
     .              DABS(knot(i1)%zv(2)-knot(i2)%zv(3)).LT.DTOL) THEN
c...          Matching sides 12 and 34:
c              WRITE(0,*) 'SIDE14:',i1,i2
              index2 = i2
              RETURN
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      RETURN
 99   STOP
      END

c
c ========================================================================
c
      SUBROUTINE MoveKnot_SL(knot1,knot2)
      USE mod_grid
      IMPLICIT none

c      TYPE type_cell_old
c        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
c        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
c      ENDTYPE type_cell_old
      TYPE(type_grid_cell) :: knot1,knot2      

      INTEGER i1

      knot2%index  = knot1%index
      knot2%ik     = knot1%ik
      knot2%ir     = knot1%ir
      knot2%rzone  = knot1%rzone
      knot2%zzone  = knot1%zzone
      knot2%xpt    = knot1%xpt
      knot2%rcen   = knot1%rcen
      knot2%zcen   = knot1%zcen
      knot2%bratio = knot1%bratio
      DO i1 = 1, 4
        knot2%rv(i1) = knot1%rv(i1)
        knot2%zv(i1) = knot1%zv(i1)
      ENDDO

      RETURN
 99   STOP
      END
c
c ========================================================================
c
      SUBROUTINE ReadGeneralisedGrid_SL(gridunit,ik,ir,
     .                                  rshift,zshift,indexiradj)
      USE mod_sol28_global
      USE mod_grid
      USE mod_grid_divimp
      IMPLICIT none

      INTEGER gridunit,ik,ir,indexiradj
      REAL    rshift,zshift

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'pindata'

c..TMP
      CHARACTER title*174,desc*1024,job*72,equil*60
      REAL      facta(-1:MAXIZS),factb(-1:MAXIZS)

      INTEGER, PARAMETER :: NUMZONE = 5
      REAL*8,  PARAMETER :: TOL = 1.0D-06

      INTEGER   i1,i2,z1,r1,kind,nxpt,ixpt(0:2),cxpt(0:2),i3,
c      INTEGER   nknot,i1,i2,z1,r1,kind,nxpt,ixpt(0:2),cxpt(0:2),i3,
     .          izone(NUMZONE+1,NUMZONE),newi1,icore(0:2),id,tmpnks,
     .          ikmax,irmax,ir1,istart,
     .          numpsi,ikpsi(MAXNRS),irpsi(MAXNRS)
      LOGICAL   cont,deleteknot,output,swap
      REAL      vrmin,vzmin,vrmax,vzmax,rspan,zspan,area,valpsi(MAXNRS)
      REAL*8    rvdp(4),zvdp(4),areadp
      CHARACTER buffer*1000

c      INTEGER, ALLOCATABLE :: imap(:,:)

c      TYPE type_cell_old
c        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
c        REAL    :: rcen,zcen,bratio,rv(4),zv(4)
c      ENDTYPE type_cell_old
c      TYPE(type_cell_old),ALLOCATABLE :: knot(:)

      output = .FALSE.

      ALLOCATE(knot(0:MAXNKS*MAXNRS))
      ALLOCATE(imap(MAXNKS,0:MAXNRS))

c...  Find the start of the cell/knot information in the grid file:
      WRITE(buffer,'(1000X)')
      DO WHILE (buffer(4:8).NE.'=====')
        READ(gridunit,'(A10)',END=98) buffer
      ENDDO

c...  Read the knot data:
      nknot = 0
      DO WHILE(nknot.EQ.0.OR.buffer(4:10).EQ.'Element')
        nknot = nknot + 1
        READ(gridunit,80,END=97) knot(nknot)%index,
     .                           knot(nknot)%ik   ,knot(nknot)%ir,
     .                           knot(nknot)%rv(2),knot(nknot)%zv(2),
     .                           knot(nknot)%rv(3),knot(nknot)%zv(3)
        READ(gridunit,81,END=97) knot(nknot)%bratio,
     .                           knot(nknot)%rcen ,knot(nknot)%zcen
        READ(gridunit,82,END=97) knot(nknot)%rv(1),knot(nknot)%zv(1),
     .                           knot(nknot)%rv(4),knot(nknot)%zv(4)
        knot(nknot)%nv = 4
c...    Dividing line:       
        READ(gridunit,*)
        READ(gridunit,'(A10)',END=20) buffer
        BACKSPACE(gridunit)
      ENDDO
 80   FORMAT(10X,I5,4X,I3,1x,i3,4x,e17.10e2,1x,e17.10e2,8x,e17.10e2,1x,
     .       E17.10E2)
 81   FORMAT(18x,e17.10e2,14x,e17.10e2,1x,e17.10e2)
 82   FORMAT(30x,e17.10e2,1x,e17.10e2,8x,e17.10e2,1x,e17.10e2)
c...  End of file continuation:
 20   CONTINUE

c...  Delete zero volume cells:

c...  Strip those boundary cells:
c     jdemod
      IF (output) then
         WRITE(0,*) 'CELL 1'
         WRITE(0,*) '2     :',knot(1)%rv(2),knot(1)%zv(2)
         WRITE(0,*) '3     :',knot(1)%rv(3),knot(1)%zv(3)
         WRITE(0,*) '1     :',knot(1)%rv(1),knot(1)%zv(1)
         WRITE(0,*) '4     :',knot(1)%rv(4),knot(1)%zv(4)
         WRITE(0,*) 'cen   :',knot(1)%rcen,knot(1)%zcen
         WRITE(0,*) 'bratio:',knot(1)%bratio

         WRITE(0,*) 'STRIPPING...'
         WRITE(6,*) 'STRIPPING...'
 

      endif
c
      ikmax = 0
      irmax = 0
      DO i1 = 1, nknot
        IF (knot(i1)%ik.GT.ikmax) ikmax = knot(i1)%ik
        IF (knot(i1)%ir.GT.irmax) irmax = knot(i1)%ir
      ENDDO
c...  (This is not the most efficient way of doing things -- certaintly it would
c      be better just to avoid storing the cells as the grid file is read in --    <---WHAT?
c      but it is general, which is the goal here, and makes minimal assumptions 
c      about the structure of the grid file):
c
c     jdemod - the following code loops through the grid and removes all cells for which the 
c              coordinates of vertex 3=2 and 4=1 - these are the parallel to the field line 
c              vertices - this will remove cells at the targets but not any boundary rings
c
      i1 = 1
      DO WHILE(i1.LE.nknot)
c        IF     (.TRUE..AND.  ! FOR SXD GRID
        IF     (opt%f_grid_strip.EQ.1.AND.
     .          (knot(i1)%ik.EQ.0.OR.knot(i1)%ik.EQ.ikmax.OR.     ! virtual cells on end of rings
     .           knot(i1)%ir.EQ.0.OR.knot(i1)%ir.EQ.irmax)) THEN  ! virtual rings
c...      Remove boundary knots:
          deleteknot = .TRUE.
        ELSEIF (ABS(knot(i1)%rv(3)-knot(i1)%rv(2)).LT.TOL.AND.
     .          ABS(knot(i1)%zv(3)-knot(i1)%zv(2)).LT.TOL.AND.
     .          ABS(knot(i1)%rv(4)-knot(i1)%rv(1)).LT.TOL.AND.
     .          ABS(knot(i1)%zv(4)-knot(i1)%zv(1)).LT.TOL) THEN
c...      Also get rid of zero volume cells, which can be present in UEDGE
c         double null grids.  The above condition is the best identifier
c         for these (for grids generated with UEDGE anyway):
          deleteknot = .TRUE.
        ELSE
c...      Cell to be kept, advance index:
          deleteknot = .FALSE.
          i1 = i1 + 1
        ENDIF
c...    Delete the knot:
        IF (deleteknot) THEN
          IF (i1.LT.nknot) THEN
            DO i2 = i1, nknot
              CALL MoveKnot_SL(knot(i2+1),knot(i2))
            ENDDO
          ENDIF
          nknot = nknot - 1
c          jdemod
c          IF (output) then
c             WRITE(0,*) 'GONE',i1
c             WRITE(6,*) 'GONE',i1
c          endif

        ENDIF
      ENDDO
c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'DONE'
         WRITE(6,*) 'DONE'
      endif


c...  Search the grid for remaining virtual cells (typically zero-volume):
c      DO i1 = 1, nknot
c        DO i2 = 1, knot(i1)%nv
c          rvdp(i2) = DBLE(knot(i1)%rv(i2))
c          zvdp(i2) = DBLE(knot(i1)%zv(i2))
c        ENDDO
c        areadp = 0.0
c        DO i2 = 1, knot(i1)%nv
c          i3 = i2 + 1
c          IF (i2.EQ.knot(i1)%nv) i3 = 1
c          areadp = areadp + (rvdp(i3) * zvdp(i2) -
c     .                       rvdp(i2) * zvdp(i3))
c        ENDDO
c        area = 0.5 * SNGL(DABS(areadp))
c
c         cont = .FALSE.
c          IF (DABS(rvdp(3)-rvdp(2)).LT.TOL.AND.
c     .        DABS(zvdp(3)-zvdp(2)).LT.TOL.AND.
c     .        DABS(rvdp(4)-rvdp(1)).LT.TOL.AND.
c     .        DABS(zvdp(4)-zvdp(1)).LT.TOL) cont = .TRUE.
c
c
c        IF (area.LT.1.0E-08) 
c     .    WRITE(0,*) 'I!,AREA:',knot(i1)%index,area,cont
c      ENDDO
c
c      STOP 'sdfsd'

c...  R,Z shifts:

c...  Assign knot sector (for improved efficiency in the search routines):
      vrmin =  HI
      vrmax = -HI
      vzmin =  HI
      vzmax = -HI
      DO i1 = 1, nknot      
        DO i2 = 1, knot(i1)%nv
          IF (knot(i1)%rv(i2).LT.vrmin) vrmin = knot(i1)%rv(i2)
          IF (knot(i1)%rv(i2).GT.vrmax) vrmax = knot(i1)%rv(i2)
          IF (knot(i1)%zv(i2).LT.vzmin) vzmin = knot(i1)%zv(i2)
          IF (knot(i1)%zv(i2).GT.vzmax) vzmax = knot(i1)%zv(i2)
        ENDDO
      ENDDO
      vrmin = vrmin - 0.001
      vrmax = vrmax + 0.001
      vzmin = vzmin - 0.001
      vzmax = vzmax + 0.001
c
c     jdemod
      IF (output) then
         WRITE(0,*) 'MIN,MAX:',vrmin,vrmax,vzmin,vzmax
         WRITE(6,*) 'MIN,MAX:',vrmin,vrmax,vzmin,vzmax
      endif
c      IF (output) WRITE(0,*) 'MIN,MAX:',vrmin,vrmax,vzmin,vzmax
c...  Assign knot to search zone:
      rspan = (vrmax - vrmin) / REAL(NUMZONE)
      zspan = (vzmax - vzmin) / REAL(NUMZONE)
      DO i1 = 1, nknot
        knot(i1)%rzone = INT( (knot(i1)%rcen - vrmin) / rspan ) + 1
        knot(i1)%zzone = INT( (knot(i1)%zcen - vzmin) / zspan ) + 1
c        WRITE(0,*) 'SPAN:',i1,knot(i1)%rzone,knot(i1)%zzone
      ENDDO

c...  Index cells by zone, work now for saved time later:
      kind = 1
      DO z1 = 1, NUMZONE
        DO r1 = 1, NUMZONE
          izone(r1,z1) = kind
          DO i1 = kind, nknot
            IF     (knot(i1)%rzone.EQ.r1.AND.knot(i1)%zzone.EQ.z1) THEN
              IF (i1.EQ.kind) THEN
c...            Do nothing:
              ELSE
c...            Swap knots:
                CALL MoveKnot_SL(knot(kind),knot(0)   )
                CALL MoveKnot_SL(knot(i1)  ,knot(kind))
                CALL MoveKnot_SL(knot(0)   ,knot(i1)  )
c                WRITE(0,*) 'ZONING SWAP:',kind,i1
              ENDIF
              kind = kind + 1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c...  
      DO i1 = 1, NUMZONE-1
        izone(NUMZONE+1,i1) = izone(1,i1+1)
      ENDDO
      izone(NUMZONE+1,NUMZONE) = nknot + 1

c      DO i1 = 1, nknot
c        WRITE(0,*) 'ZONED:',i1,knot(i1)%rzone,knot(i1)%zzone
c      ENDDO
c      DO z1 = 1, NUMZONE
c        DO r1 = 1, NUMZONE
c          WRITE(0,*) 'IZONE:',izone(r1,z1) 
c        ENDDO
c      ENDDO

c...  Search for x-point(s):
c
c     jdemod - it appears that the X-point finding algorithm used is the following:
c            - the only cells on the grid which can have an identical vertex - both index and value 
c              and not be the same cell must occur at the Xpoint - the shared vertex in this 
c              case IS the Xpoint. This is the condition tested for when FindKnot_SL is called with 
c              a 1. Search efficiency has been enhanced by using the zones set up above. 
c            - Zero volume cells would be an issue with this algorithm - code above this has 
c              apparently removed zero volume cells where vertices 3=2 and 1=4 - however, it would 
c              appear that boundary rings around the plasma and in the core/PFZ have not been 
c              removed.
c
      nxpt = 0
      DO i1 = 1, nknot
        IF (knot(i1)%xpt.NE.0) CYCLE

        IF (nxpt.EQ.2) THEN
c          WRITE(0,*)
c          WRITE(0,*) '--------------------------------------------'
c          WRITE(0,*) '-   MORE THAN 2 XPTS FOUND, IGNORING...    -'
c          WRITE(0,*) '--------------------------------------------'
c          WRITE(0,*)
          CALL WN('ReadGeneralisedGrid_SL','More than 2 x-points '//
     .            'found, ignoring...')
          EXIT
        ENDIF

        CALL FindKnot_SL(NUMZONE,izone,1,i1,i2)

        IF (i2.NE.-1) THEN
c
c     jdemod - the code appears to assume that the midplane is at 0.0 - this 
c              should probably be replaced with the zc value defining the 
c              center of the confined plasma. 
c
c...      Select the appropriate cell, whichever is closest to the midplane will
c         be the cell associated with the core (which we want to build first):
          IF    (knot(i1)%zcen.LT.0.0.AND.knot(i2)%zcen.LT.0.0) THEN
            nxpt = nxpt + 1
            IF (knot(i1)%zcen.GT.knot(i2)%zcen) THEN
              ixpt(nxpt) = i1
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ELSE
              ixpt(nxpt) = i2
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ENDIF
          ELSEIF(knot(i1)%zcen.GT.0.0.AND.knot(i2)%zcen.GT.0.0) THEN
            nxpt = nxpt + 1
c
c           jdemod
            IF (output) then 
               WRITE(0,*) '  >> ',i1,i2
               WRITE(0,*) '     ',knot(i1)%zcen,knot(i2)%zcen
               WRITE(6,*) '  >> ',i1,i2
               WRITE(6,*) '     ',knot(i1)%zcen,knot(i2)%zcen
            endif
c            IF (output) WRITE(0,*) '  >> ',i1,i2
c            IF (output) WRITE(0,*) '     ',knot(i1)%zcen,knot(i2)%zcen
            IF (knot(i1)%zcen.LT.knot(i2)%zcen) THEN
              ixpt(nxpt) = i1
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ELSE
              ixpt(nxpt) = i2
              knot(i1)%xpt = i2
              knot(i2)%xpt = i1
            ENDIF
          ELSE
            CALL ER('Readgeneralisedgrid','Unrecognized '//
     .              'x-point configuration',*99)
          ENDIF
        ENDIF
      ENDDO

c...  No x-points found, for some reason:
      IF (nxpt.EQ.0) 
     .  CALL ER('Readgeneralisedgrid','No x-points found',*99)

      IF (output) THEN
        DO i1 = 1, nxpt
          WRITE(0,*) 'XPTS:',i1,
     .               ixpt(i1),knot(ixpt(i1))%index,
     .               knot(ixpt(i1))%xpt,
     .               knot(knot(ixpt(i1))%xpt)%index
          WRITE(6,*) 'XPTS:',i1,
     .               ixpt(i1),knot(ixpt(i1))%index,
     .               knot(ixpt(i1))%xpt,
     .               knot(knot(ixpt(i1))%xpt)%index
        ENDDO
      ENDIF
c
c     jdemod - the Xpoint finding algorithm returns the two cells that share the 
c              Xpoint vertex as index 1. One of these cells should be the cell
c              in the SOL adjacent to the first cell on the core ring at the
c              Xpoint. The second of these cells is below the Xpoint adjacent to 
c              PFZ. By using the cell "closer to the midplane" it should choose
c              the cell adjacent to the confined plasma - the cell sharing the side
c              with vertices where 1 = 2 and 4 = 3 should be the first cell on the 
c              core ring. 
c
c              The code then walks inward from the Xpoint finding the cell on the
c              innermost core ring corresponding to the first cell on the ring at the 
c              Xpoint. The variable cxpt records the number of rings from the Xpoint
c              to the innermost ring. If this value is the same for two different
c              Xpoints then the grid is a connected double null. If not - the difference
c              in the two values should define the number of rings in the secondary 
c              plasma between the two Xpoints for the DDN plasma configuration. 
c
c...  Searching for the start of the core center ring:
      cxpt = 0
      icore = 0
      DO i3 = 1, nxpt
        newi1 = ixpt(i3)
        cont = .TRUE.
        DO WHILE (cont)
          i1 = newi1
          cxpt(i3) = cxpt(i3) + 1
          cont = .FALSE.
c
c     jdemod - when the search has moved all the way inward and 
c              can no longer find a cell with a matching side then
c              i2=-1 is returned and the code moves onto any other
c              Xpoints. 
c
          CALL FindKnot_SL(NUMZONE,izone,2,i1,i2)

          IF (i2.NE.-1) THEN
            cont = .TRUE.
            newi1 = i2 
            icore(i3) = i2
          ENDIF
        ENDDO
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'Cxpt:',ixpt(i3),cxpt(i3),icore(i3)
           WRITE(6,*) 'Cxpt:',ixpt(i3),cxpt(i3),icore(i3)
        endif
c        IF (output) WRITE(0,*) 'Cxpt:',ixpt(i3),cxpt(i3),icore(i3)
      ENDDO

      IF (nxpt.GT.1) THEN
c...    Check that the x-points are ordered properly, with the primary x-point
c       at index 1, and whether or not the double-null grid is connected:
        swap = .FALSE.
        IF     (nxpt.GT.1.AND.cxpt(1).EQ.cxpt(2)) THEN
c...      Connected:
          IF (knot(ixpt(1))%zcen.GT.0.0) swap = .TRUE.
          connected = .TRUE.
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'CONNECTED DN DETECTED'
             WRITE(6,*) 'CONNECTED DN DETECTED'
          endif
c          IF (output) WRITE(0,*) 'CONNECTED DN DETECTED'
        ELSEIF (nxpt.GT.1.AND.cxpt(1).GT.cxpt(2)) THEN
          swap = .TRUE.
          connected = .FALSE.
        ENDIF
        IF (swap) THEN
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'SWAPPING X-POINTS'
             WRITE(6,*) 'SWAPPING X-POINTS'
          endif
c          IF (output) WRITE(0,*) 'SWAPPING X-POINTS'
          ixpt (0) = ixpt (1)
          ixpt (1) = ixpt (2)
          ixpt (2) = ixpt (0)
          cxpt (0) = cxpt (1)
          cxpt (1) = cxpt (2)
          cxpt (2) = cxpt (0)
          icore(0) = icore(1)
          icore(1) = icore(2)
          icore(2) = icore(0)
        ENDIF
      ENDIF

c...  Make sure that x-point knot indeces are ordered properly, with the inner (lower radius)
c     of each pair listed in IXPT:
c      DO i1 = 1, nxpt
c        IF (knot(ixpt(i1))%rcen.GT.knot(knot(ixpt(i1))%xpt)%rcen) THEN
c          jdemod
c          IF (output) then 
c              WRITE(0,*) 'SWAPPING X-POINT PAIR',i1
c              WRITE(6,*) 'SWAPPING X-POINT PAIR',i1
c          endif       
c          IF (output) WRITE(0,*) 'SWAPPING X-POINT PAIR',i1
c          ixpt(i1) = knot(ixpt(i1))%xpt
c        ENDIF
c      ENDDO

c...  Location of the primary separatrix is known:
      irsep  = cxpt(1)
      irsep2 = irsep

c...  Build the grid:
c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING CORE RINGS'
         WRITE(6,*) 'PROCESSING CORE RINGS'
      endif
c      IF (output) WRITE(0,*) 'PROCESSING CORE RINGS'

c...  Start with the core rings:
      ik = 1
      ir = 1
      i1 = icore(1) 
      DO WHILE(ir.LT.irsep)
        imap(1,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move along the ring:
          CALL FindKnot_SL(NUMZONE,izone,3,i1,i2)
c          WRITE(0,*) '>>>',i1,i2,ir,imap(1,ir)  
          IF (i2.NE.-1) THEN
            IF (i2.NE.imap(1,ir)) THEN
              i1 = i2 
              ik = ik + 1
              imap(ik,ir) = i1 
              cont = .TRUE.
c
c             jdemod
              IF (output) then 
                 WRITE(0,*) 'CORE MAP:',ik,ir,i1
                 WRITE(6,*) 'CORE MAP:',ik,ir,i1
              endif
c              IF (output) WRITE(0,*) 'CORE MAP:',ik,ir,i1
            ENDIF
          ELSE
            CALL ER('ReadGeneralisedGrid','Bad IK step',*99)
          ENDIF
        ENDDO
        nks(ir) = ik
c...    Step outward, still in the core:        
        CALL FindKnot_SL(NUMZONE,izone,4,imap(1,ir),i2)
        IF (i2.NE.-1) THEN        
          i1 = i2
          ik = 1
          ir = ir + 1
        ELSE
          CALL ER('Readgeneralisedgrid','Bad IR step',*99)
        ENDIF
      ENDDO

c...  SOL rings: 
c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING SOL RINGS'
         WRITE(6,*) 'PROCESSING SOL RINGS'
      endif
c      IF (output) WRITE(0,*) 'PROCESSING SOL RINGS'
c
c     jdemod - Doesn't this assume an Xpoint down configuration? At least as far as the 
c              "high field side" reference goes? I think the code itself still works. 
c

c...  Step out of the core on the high field side:
      i1 = imap(1,irsep-1)
      CALL FindKnot_SL(NUMZONE,izone,4,i1,i2)
      IF (i2.NE.-1) THEN  
c...    Move down to the target:
        i1 = i2
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindKnot_SL(NUMZONE,izone,5,i1,i2)
          IF (i2.NE.-1) THEN 
            i1 = i2
            cont = .TRUE.
          ENDIF
        ENDDO
      ELSE
        CALL ER('Readgeneralisedgrid','Bad IR step to SOL',*99)
      ENDIF
c...  Target located, start mapping the SOL:
      ik = 1
      ir = irsep
      imap(ik,ir) = i1
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
c...    Move along the ring:
        CALL FindKnot_SL(NUMZONE,izone,3,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2 
          ik = ik + 1
          imap(ik,ir) = i1
          cont = .TRUE.
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'INNER SOL MAP:',ik,ir,i1
             WRITE(6,*) 'INNER SOL MAP:',ik,ir,i1
          endif
c          IF (output) WRITE(0,*) 'INNER SOL MAP:',ik,ir,i1
        ENDIF
c...    Step radially outward if ring is finished:
        IF (.NOT.cont) THEN
          nks(ir) = ik
          i1 = imap(1,ir)
          CALL FindKnot_SL(NUMZONE,izone,4,i1,i2)          
          IF (i2.NE.-1) THEN
            i1 = i2
            ik = 1
            ir = ir + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
               WRITE(0,*) 'INNER SOL MAP NEW RING:',ik,ir,i1
               WRITE(6,*) 'INNER SOL MAP NEW RING:',ik,ir,i1
            endif
c            IF (output) WRITE(0,*) 'INNER SOL MAP NEW RING:',ik,ir,i1
          ENDIF
        ENDIF
      ENDDO
      irwall = ir
      irtrap = ir + 1
      nrs = ir

      IF (nxpt.GT.1) THEN

        irsep2 = irsep + cxpt(2) - cxpt(1) - 1
c
c       jdemod
        IF (output) then
           WRITE(0,*) 'PROCESSING LOW FIELD SOL'
           WRITE(6,*) 'PROCESSING LOW FIELD SOL'
        endif
c        IF (output) WRITE(0,*) 'PROCESSING LOW FIELD SOL'

c...    Process the low field side looking for any rings that
c       were not processed when looking around the high field side (which
c       is usually the case for double-nulls):

c...    Register all knots that have been mapped to the grid:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            knot(imap(ik,ir))%map = 1
          ENDDO
        ENDDO
c
c     jdemod - why not use - i1=imap(nks(irsep-1),irsep-1) ?
c
c...    Start with the first cell on the outer-most core ring and move to 
c       the last cell on the same ring:
        i1 = imap(1,irsep-1)
        CALL FindKnot_SL(NUMZONE,izone,5,i1,i2)
        IF (i2.EQ.-1) 
     .    CALL ER('Readgeneralisedgrid','Core map problems',*99)
        i1 = i2
c...    Move into the SOL:
        CALL FindKnot_SL(NUMZONE,izone,4,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2
          IF (knot(i2)%map.EQ.1) THEN  
c...        Keep moving outward until a cell with no assigned mapping is
c           found:
            cont = .TRUE.
            DO WHILE(cont)
              cont = .FALSE.
              CALL FindKnot_SL(NUMZONE,izone,4,i1,i2)
              IF (i2.NE.-1) THEN 
                i1 = i2
                IF (knot(i1)%map.NE.0) cont = .TRUE.
              ELSE
c...            Either an error or a single-null grid is being tested:
                STOP 'SINGLE NULL GRID BEING TESTED OR ERROR?'
              ENDIF
            ENDDO
          ENDIF
        ELSE
          CALL ER('Readgeneralisedgrid','Bad IR step to SOL',*99)
        ENDIF

c...    An unmapped cell has been found, proceed to target:
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindKnot_SL(NUMZONE,izone,5,i1,i2)
          IF (i2.NE.-1) THEN 
            i1 = i2
            cont = .TRUE.
          ENDIF
        ENDDO

c...    Target located, start mapping the low field SOL:
        ik = 1
        ir = nrs + 1
        IF (connected) irsep2 = ir
        imap(ik,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move along the ring:
          CALL FindKnot_SL(NUMZONE,izone,3,i1,i2)
          IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
               WRITE(0,*) 'OUTER SOL MAP:',ik,ir,i1
               WRITE(6,*) 'OUTER SOL MAP:',ik,ir,i1
            endif
c            IF (output) WRITE(0,*) 'OUTER SOL MAP:',ik,ir,i1
          ENDIF
c...      Step radially outward if ring is finished:
          IF (.NOT.cont) THEN
c
c           jdemod
            IF (output) then 
               WRITE(0,*) 'STEPPING OUT:',ik,ir,i1
               WRITE(6,*) 'STEPPING OUT:',ik,ir,i1
            endif
c            IF (output) WRITE(0,*) 'STEPPING OUT:',ik,ir,i1
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindKnot_SL(NUMZONE,izone,4,i1,i2)          
            IF (i2.NE.-1) THEN
              i1 = i2
              ik = 1
              ir = ir + 1
              imap(ik,ir) = i1
              cont = .TRUE.
c
c             jdemod
              IF (output) then 
                 WRITE(0,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
                 WRITE(6,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
              endif
c              IF (output) WRITE(0,*) 'OUTER SOL MAP NEW RING:',ik,ir,i1
            ELSE
c...          Assume the outer boundary of the grid:
c
c             jdemod
              IF (output) then 
                 WRITE(0,*) 'ASSUMING OUTER GRID BOUNDARY'
                 WRITE(6,*) 'ASSUMING OUTER GRID BOUNDARY'
              endif
c              IF (output) WRITE(0,*) 'ASSUMING OUTER GRID BOUNDARY'
            ENDIF
          ENDIF
        ENDDO
        irwall = ir
        irtrap = ir + 1
        nrs = ir

c...    Register all knots that have been mapped to the grid:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            knot(imap(ik,ir))%map = 1
          ENDDO
        ENDDO



c       jdemod
        IF (output) then 
           WRITE(0,*) 'PROCESSING SECONDARY PFZ'
           WRITE(6,*) 'PROCESSING SECONDARY PFZ'
        endif
c        IF (output) WRITE(0,*) 'PROCESSING SECONDARY PFZ'
c...    Process the secondary x-point PFR, which is just considered part of the
c       SOL for generalized grids:
        IF (.FALSE..AND.connected) THEN  ! NEEDED FOR SXD GRID
          i1 = ixpt(2)  ! OK, not sure if this always needed to be hardcoded for connected grids,
                        ! but adding the option here - quasi-BUG!  16.06.09 - SL
        ELSE
          i1 = knot(ixpt(2))%xpt
        ENDIF
        IF (output) WRITE(6,*) 'I1 XPT:',i1

c...    Move into the PFR:
        CALL FindKnot_SL(NUMZONE,izone,2,i1,i2)
        IF (i2.EQ.-1) 
     .    CALL ER('Readgeneralisedgrid','PFR2 problems',*99)
        i1 = i2
        ik = 1
        ir = nrs + 1
        imap(ik,ir) = i1


c...    Proceed to target:
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'LOOKING FOR TARGET'
           WRITE(6,*) 'LOOKING FOR TARGET'
        endif
c        IF (output) WRITE(0,*) 'LOOKING FOR TARGET'
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
          CALL FindKnot_SL(NUMZONE,izone,5,i1,i2)
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'MOVING',i1,i2,istart
             WRITE(6,*) 'MOVING',i1,i2,istart
             IF (i2.GT.0) THEN
               WRITE(6,*) '  KNOT R:',SNGL(knot(i1)%rv(1:4))
               WRITE(6,*) '        :',SNGL(knot(i2)%rv(1:4))
               WRITE(6,*) '  KNOT Z:',SNGL(knot(i1)%zv(1:4))
               WRITE(6,*) '        :',SNGL(knot(i2)%zv(1:4))
            ELSE
               WRITE(6,*) '        : i2 = -1'
            ENDIF
          endif
c          IF (output) WRITE(0,*) 'MOVING',i1,i2,istart
          IF (i2.NE.-1) THEN 
            i1 = i2
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
          ENDIF
        ENDDO
 
c        STOP 'sdfsd'

c...    Target located, start mapping the secondary PFR:
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'TARGET LOCATED'
           WRITE(6,*) 'TARGET LOCATED'
        endif
c        IF (output) WRITE(0,*) 'TARGET LOCATED'
        ik = 1
        ir = nrs + 1
        imap(ik,ir) = i1
        cont = .TRUE.
        DO WHILE(cont)
          cont = .FALSE.
c...      Move along the ring:
          CALL FindKnot_SL(NUMZONE,izone,3,i1,i2)
c
c         jdemod
          IF (output) then 
             WRITE(0,*) 'MOVING',i1,i2
             WRITE(6,*) 'MOVING',i1,i2
          endif
c          IF (output) WRITE(0,*) 'MOVING',i1,i2
          IF (i2.NE.-1) THEN
            i1 = i2 
            ik = ik + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
               WRITE(0,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
               WRITE(6,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
            endif
c            IF (output) WRITE(0,*) '2ND PFR MAP:',ik,ir,i1,knot(i1)%zcen
          ENDIF
c...      Step radially outward if ring is finished:
          IF (.NOT.cont) THEN
            nks(ir) = ik
            i1 = imap(1,ir)
            CALL FindKnot_SL(NUMZONE,izone,2,i1,i2)          
            IF (i2.NE.-1) THEN
              i1 = i2
              ik = 1
              ir = ir + 1
              imap(ik,ir) = i1
              cont = .TRUE.
c
c             jdemod
              IF (output) then  
                WRITE(0,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
     >             knot(i1)%zcen
                WRITE(6,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
     >             knot(i1)%zcen
              endif
c              IF (output) 
c     .          WRITE(0,*) '2ND PFR MAP NEW RING:',ik,ir,i1,
c     .          knot(i1)%zcen
            ELSE
c...          Assume the outer boundary of the grid:
            ENDIF
          ENDIF
        ENDDO
        irwall = ir
        irtrap = ir + 1
        nrs = ir

      ENDIF ! Done processing double-null rings


      IF (.FALSE.) THEN   ! DEBUG
        nks(ir) = ik
        irwall = ir
        irtrap = ir + 1
        nrs = ir
        id = 0
        DO ir = 1, nrs
          DO ik = 1, nks(ir)        
            i1 = imap(ik,ir)
            rs(ik,ir) = knot(i1)%rcen
            zs(ik,ir) = knot(i1)%zcen
            bratio(ik,ir) = knot(i1)%bratio
            id = id + 1
            korpg(ik,ir) = id
            nvertp(id) = knot(i1)%nv
            DO i2 = 1, nvertp(id)
              rvertp(i2,id) = knot(i1)%rv(i2)
              zvertp(i2,id) = knot(i1)%zv(i2)
            ENDDO
          ENDDO
        ENDDO
        ikto = -1
        ikti = -1
        DO ik = 1, nks(irsep)
          IF (connected) THEN
            IF (imap(ik,irsep).EQ.ixpt(1)          ) ikto = ik    ! Not sure this will always work... 
            IF (imap(ik,irsep).EQ.knot(ixpt(2))%xpt) ikti = ik -1
c            IF (imap(ik,irsep).EQ.ixpt(1)) ikto = ik - 1  ! Not sure this will always work...
c            IF (imap(ik,irsep).EQ.ixpt(2)) ikti = ik     
          ENDIF
        ENDDO
        CALL SaveSolution
        CALL OutputData(86,'MAST!')
        title = '...'
        desc  = 'Call to STORE from DumpGrid'
        job   = 'Call to STORE from DumpGrid'
        equil = 'Call to STORE from DumpGrid'
        WRITE(0,*) 'CALLING STORE'
        CALL Store(title,desc,1,job,equil,facta,factb,1,1)
        WRITE(0,*) 'FUN WITH MAST GRIDS!'
        STOP
      ENDIF



c



c
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING PRIMARY PFZ'
         WRITE(6,*) 'PROCESSING PRIMARY PFZ'
      endif   
c      IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ'

c...  Process the primary x-point PFR:
c      IF (connected) THEN
c        i1 = ixpt(1)
c      ELSE
        i1 = knot(ixpt(1))%xpt
c      ENDIF
c...  Move into the PFR:
      CALL FindKnot_SL(NUMZONE,izone,2,i1,i2)
      IF (i2.EQ.-1) 
     .  CALL ER('Readgeneralisedgrid','PFR1 problems',*99)
c
c     jdemod
      IF (output) then
         WRITE(0,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
         WRITE(6,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
      endif


c      IF (output) 
c     .  WRITE(0,*) '  INTO PFZ',knot(i1)%index,knot(i2)%index
      i1 = i2
c...  Proceed to target:
      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
        CALL FindKnot_SL(NUMZONE,izone,5,i1,i2)
c        jdemod
c        IF (output) then 
c           WRITE(0,*) '  TO TARGET',i1,i2
c           WRITE(0,*) '  TO TARGET',knot(i1)%index
c           WRITE(0,*) '  TO TARGET',knot(i2)%index
c           WRITE(6,*) '  TO TARGET',i1,i2
c           WRITE(6,*) '  TO TARGET',knot(i1)%index
c           WRITE(6,*) '  TO TARGET',knot(i2)%index
c        endif
c        IF (output) WRITE(0,*) '  TO TARGET',i1,i2
c        IF (output) WRITE(0,*) '  TO TARGET',knot(i1)%index
c        IF (output) WRITE(0,*) '  TO TARGET',knot(i2)%index
c        STOP 'tet'
        IF (i2.NE.-1) THEN 
          i1 = i2
          cont = .TRUE.
        ENDIF
      ENDDO
c...  Target located, start mapping the primary PFR:
c     
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
         WRITE(6,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
      endif


c      IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ, TARGET LOCATED'
      ik = 1
      ir = nrs + 1
      imap(ik,ir) = i1

      cont = .TRUE.
      DO WHILE(cont)
        cont = .FALSE.
c...    Move along the ring:
        CALL FindKnot_SL(NUMZONE,izone,3,i1,i2)
        IF (i2.NE.-1) THEN
          i1 = i2 
          ik = ik + 1
          imap(ik,ir) = i1
          cont = .TRUE.
c          jdemod
c          IF (output) then 
c             WRITE(0,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
c             WRITE(6,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
c          endif
c          IF (output) WRITE(0,*) '1ST PFR MAP:',ik,ir,i1,knot(i1)%zcen
        ENDIF
c...    Step radially outward if ring is finished:
        IF (.NOT.cont) THEN
          nks(ir) = ik
          i1 = imap(1,ir)
          CALL FindKnot_SL(NUMZONE,izone,2,i1,i2)          
          IF (i2.NE.-1) THEN
            i1 = i2
            ik = 1
            ir = ir + 1
            imap(ik,ir) = i1
            cont = .TRUE.
c
c           jdemod
            IF (output) then 
               WRITE(0,*) '1ST PFR MAP NEW RING:',ik,ir,i1,knot(i1)%zcen
               WRITE(6,*) '1ST PFR MAP NEW RING:',ik,ir,i1,knot(i1)%zcen
            endif
c            IF (output) 
c     .        WRITE(0,*) '1ST PFR MAP NEW RING:',ik,ir,i1,knot(i1)%zcen
          ELSE
c...        Assume the outer boundary of the grid:
          ENDIF
c         jdemod
          IF (output) then 
             WRITE(0,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
             WRITE(6,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
          endif
c          IF (output) WRITE(0,*) 'PROCESSING PRIMARY PFZ, BUZZING...'
        ENDIF
      ENDDO
      nrs = ir

c...  Need to reorder the rings in the primary PFZ:
      DO i1 = 0, (nrs-irtrap+1)/2-1
c
c       jdemod
        IF (output) then 
           WRITE(0,*) 'I1???=',i1
           WRITE(6,*) 'I1???=',i1
        endif
c        IF (output) WRITE(0,*) 'I1???=',i1
        tmpnks = nks(irtrap+i1)
        DO ik = 1, nks(irtrap+i1)
          imap(ik,0) = imap(ik,irtrap+i1)
        ENDDO
        nks(irtrap+i1) = nks(nrs-i1)
        DO ik = 1, nks(nrs-i1)
          imap(ik,irtrap+i1) = imap(ik,nrs-i1)
        ENDDO
        nks(nrs-i1) = tmpnks
        DO ik = 1, tmpnks
          imap(ik,nrs-i1) = imap(ik,0)
        ENDDO
      ENDDO

c...  Put grid together:
c     jdemod
      IF (output) then 
         WRITE(0,*) 'PUTTING GRID TOGETHER'
         WRITE(6,*) 'PUTTING GRID TOGETHER'
      endif
c      IF (output) WRITE(0,*) 'PUTTING GRID TOGETHER'






c...TMP
      IF (.FALSE.) THEN
        WRITE(0,*) 'CELL IN QUESTION:',i1
        WRITE(6,*) 'CELL IN QUESTION:',i1
c        i1 = i2
c        ik = 1
c        ir = nrs + 1
c        imap(ik,ir) = i1
c        nks(ir) = 1
c        nks(ir) = ik
c        nrs = nrs + 2

c        ik = 1
c        imap(ik,ir) = i2
c        irwall = nrs
c        irtrap = nrs + 1
c        nks(irwall) = 1 
c        nks(irtrap) = 1 
c        nrs = ir 
c        WRITE(0,*) 'i2:',i2,knot(i2)%map
        id = 0
        DO ir = 1, nrs
          WRITE(0,*) 'NKS:',ir,nks(ir)
          DO ik = 1, nks(ir)        
            i1 = imap(ik,ir)
            rs(ik,ir) = knot(i1)%rcen
            zs(ik,ir) = knot(i1)%zcen
            bratio(ik,ir) = knot(i1)%bratio
            id = id + 1
            korpg(ik,ir) = id
            nvertp(id) = knot(i1)%nv
            DO i2 = 1, nvertp(id)
              rvertp(i2,id) = knot(i1)%rv(i2)
              zvertp(i2,id) = knot(i1)%zv(i2)
              IF (i1.EQ.588) THEN
                IF (i2.EQ.1) zvertp(i2,id) = 0.0
c                IF (i2.EQ.1.OR.i2.EQ.2) zvertp(i2,id) = 0.0
              ENDIF
c              IF (i1.EQ.2843) THEN
c                IF (i2.EQ.1) zvertp(i2,id) = 0.0
c                IF (i2.EQ.1.OR.i2.EQ.2) zvertp(i2,id) = 0.0
c              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CALL SaveSolution
        CALL OutputData(86,'MAST!')
        title = '...'
        desc  = 'Call to STORE from DumpGrid'
        job   = 'Call to STORE from DumpGrid'
        equil = 'Call to STORE from DumpGrid'
        WRITE(0,*) 'CALLING STORE'
        CALL Store(title,desc,1,job,equil,facta,factb,1,1)
        WRITE(0,*) 'FUN WITH MAST GRIDS!'
        STOP
      ENDIF























      CALL ALLOC_GRID(MAXNKS,MAXNRS)
      id = 0
      DO ir = 1, nrs
        DO ik = 1, nks(ir)        
          i1 = imap(ik,ir)
          rs(ik,ir) = knot(i1)%rcen
          zs(ik,ir) = knot(i1)%zcen
          bratio(ik,ir) = knot(i1)%bratio
          id = id + 1
          korpg(ik,ir) = id
          nvertp(id) = knot(i1)%nv
          DO i2 = 1, nvertp(id)
            rvertp(i2,id) = knot(i1)%rv(i2)
            zvertp(i2,id) = knot(i1)%zv(i2)
          ENDDO
c...      Store these in case B2 data from Rhozansky is being loaded:
          divimp_ik(ik,ir) = knot(i1)%ik 
          divimp_ir(ik,ir) = knot(i1)%ir
        ENDDO
      ENDDO


c...  Set NPOLYP:
      npolyp  = id
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

c...  Find IKTO,IKTI:
      ikto = -1
      ikti = -1
      DO ik = 1, nks(irsep)
        IF (connected) THEN
c          IF (imap(ik,irsep).EQ.ixpt(1)) ikto = ik       ! SXD GRID!
c          IF (imap(ik,irsep).EQ.ixpt(2)) ikti = ik - 1 
          IF (imap(ik,irsep).EQ.ixpt(1)          ) ikto = ik    ! Not sure this will always work...
          IF (imap(ik,irsep).EQ.knot(ixpt(2))%xpt) ikti = ik - 1
c          IF (imap(ik,irsep).EQ.ixpt(1)) ikto = ik - 1     ! Not sure this will always work...
c          IF (imap(ik,irsep).EQ.ixpt(2)) ikti = ik - 1
c          IF (imap(ik,irsep).EQ.knot(ixpt(2))%xpt) ikti = ik - 1  ! Not sure this will always work...
        ELSE
          IF (imap(ik,irsep).EQ.ixpt(1)          ) ikto = ik - 1
          IF (imap(ik,irsep).EQ.knot(ixpt(1))%xpt) ikti = ik 
        ENDIF
      ENDDO


      IF (ikto.EQ.-1.OR.ikti.EQ.-1)
     .  CALL ER('Readgeneralisedgrid','IKTI or IKTO not found',*99)

c     CUTPT1
c     CUTPT2
c     MAXKPTS
c     MAXRINGS
c     CUTRING
c
      IF (ikto.EQ.0) THEN
        nopriv = .TRUE.
        CALL InsertRing(1  ,BEFORE,PERMANENT)
        CALL InsertRing(nrs,AFTER ,PERMANENT)
      ELSE
        nopriv = .FALSE.
c...    Add virtual rings 1 (core boundary), IRWALL (SOL) and IRTRAP (PFZ):
        CALL InsertRing(1,BEFORE,PERMANENT)
        CALL InsertRing(nrs-irsep+2,AFTER,PERMANENT)
        CALL InsertRing(nrs-irsep+3,BEFORE,PERMANENT)
      ENDIF

      cutpt1 = ikto
      cutpt2 = ikti             ! These are semi-bogus for a connected double-null...?
      cutring = irsep - 1
      maxkpts = nks(irsep)
      maxrings = irwall
      indexiradj = 1

c...TMP
c      ik = 1
c      ir = nrs + 1
c      nks(ir) = 1 
c      imap(ik,ir) = i2
c      irwall = ir
c      irtrap = ir + 1
c      nrs = ir 
c      WRITE(0,*) 'i2:',i2,knot(i2)%map

      IF (.NOT..TRUE.) THEN
c        id = 0
c        DO ir = 1, nrs
c          DO ik = 1, nks(ir)        
c            i1 = imap(ik,ir)
c            rs(ik,ir) = knot(i1)%rcen
c            zs(ik,ir) = knot(i1)%zcen
c            bratio(ik,ir) = knot(i1)%bratio
c            id = id + 1
c            korpg(ik,ir) = id
c            nvertp(id) = knot(i1)%nv
c            DO i2 = 1, nvertp(id)
c              rvertp(i2,id) = knot(i1)%rv(i2)
c              zvertp(i2,id) = knot(i1)%zv(i2)
c            ENDDO
c          ENDDO
c        ENDDO
        CALL SaveSolution
        CALL OutputData(86,'MAST!')
        title = '...'
        desc  = 'Call to STORE from DumpGrid'
        job   = 'Call to STORE from DumpGrid'
        equil = 'Call to STORE from DumpGrid'
        WRITE(0,*) 'CALLING STORE'
        CALL Store(title,desc,1,job,equil,facta,factb,1,1)
        WRITE(0,*) 'FUN WITH MAST GRIDS!'
        STOP
      ENDIF

c...  Add virtual boundary cells, which will be stripped off later:
      IF (ctargopt.EQ.0.OR.ctargopt.EQ.1.OR.ctargopt.EQ.2.OR.
     .    ctargopt.EQ.3.OR.ctargopt.EQ.6) 
     .  CALL AddPoloidalBoundaryCells      

c...  Look for PSIn data for full double null grids (code mostly 
c     from tau.d6a):
      READ(gridunit,'(A)',END=25) buffer
      IF (buffer(1:16).EQ.'PSI-DOUBLE-NULLd') THEN   ! direct assignment to each ring, post mortem...
        READ(buffer(17:),*) numpsi
c...    The PSI values are to be loaded in TailorGrid...glorious hack!
        DO i1 = 1, numpsi
          READ(gridunit,*,END=97)
        ENDDO
      ELSEIF (buffer(1:15).EQ.'PSI-DOUBLE-NULL') THEN
        READ(buffer(16:),*) numpsi
c...    The PSI values are listed with one on each line
c       indexed by knot and ring index based on the SONNET 
c       grid coordinates:
        DO i1 = 1, numpsi
          READ(gridunit,*,END=97) ikpsi(i1),irpsi(i1),valpsi(i1)
        ENDDO
c...    Assign to grid rings:
        DO ir = 2, irwall-1
c...      (Need the "-1" because a virtual core ring has been added to the grid)
          ir1 = knot(imap(1,ir-1))%ir 
          DO i1 = 1, numpsi
            IF (irpsi(i1).EQ.ir1) EXIT
          ENDDO
          IF (i1.EQ.numpsi+1) 
     .      CALL ER('Readgeneralisedgrid','Problem with PSIn',*99)
c          WRITE(0,*) '--',ir,valpsi(i1)
          psitarg(ir,1) = valpsi(i1)
          psitarg(ir,2) = valpsi(i1)
          IF (ir.LT.irsep) THEN
            psitarg(ir-1+irtrap,1) = valpsi(i1)
            psitarg(ir-1+irtrap,2) = valpsi(i1)
          ENDIF
        ENDDO
        WRITE(0,*)
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) 'BOGUS PSITARG -- ALSO USING INNER TARGET DATA ONLY'
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*)
      ELSE
        BACKSPACE gridunit
      ENDIF
 25   CONTINUE

      DEALLOCATE(knot)
      DEALLOCATE(imap)

c      IF (nrs.EQ.60) THEN
c        WRITE(0,*)
c        WRITE(0,*) '--------------------------------------------------'
c        WRITE(0,*) 'HARDCODING IRSEP2 = 30 (FOR IR = 60) '
c        WRITE(0,*) '--------------------------------------------------'
c        WRITE(0,*)
c        irsep2 = 35
c      ENDIF

c...  For consistency with original SONNET code in tau.d6a:
      ir = maxrings
      ik = maxkpts

      RETURN
 97   CALL ER('Readgeneralisedgrid','Unexpected end-of-file',*99)
 98   CALL ER('Readgeneralisedgrid','Problem accessing grid file',*99)
 99   WRITE(0,*) 'IXPT:',ixpt(1),ixpt(2)
      STOP  
      END
c
c ========================================================================
c
      SUBROUTINE ReadGeneralisedGrid_OSM(gridunit,ik,ir,
     .                                   rshift,zshift,indexiradj)
      USE mod_sol28_global
      USE mod_grid
      USE mod_grid_divimp
      IMPLICIT none

      INTEGER, INTENT(IN)  :: gridunit
      INTEGER, INTENT(OUT) :: ik,ir,indexiradj
      REAL   , INTENT(OUT) :: rshift,zshift

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'
      INCLUDE 'pindata'

      INTEGER   id,i1,i2,ir1,
     .          numpsi,ikpsi(MAXNRS),irpsi(MAXNRS)
      CHARACTER buffer*1024
      REAL      valpsi(MAXNRS)

      rshift = 0.0
      zshift = 0.0

      CALL LoadGeneralisedGrid

      irsep      = grid_load%irsep
      irsep2     = grid_load%irsep2
      irwall     = grid_load%irwall
      irtrap     = grid_load%irtrap
      nrs        = grid_load%nrs
      nks(1:nrs) = grid_load%nks(1:nrs)
      psitarg(1:nrs,1) = grid_load%psin(1:nrs)
      psitarg(1:nrs,2) = psitarg(1:nrs,1)

      IF (irsep2.EQ.-1) irsep2 = irsep  ! *** is this OK? ***


c      write(0,*) 'irsep1 ',irsep
c      write(0,*) 'irsep2 ',irsep2
c      write(0,*) 'irwall1',irwall
c      write(0,*) 'irtrap1',irtrap
c      write(0,*) 'nrs1   ',nrs


      id = 0
      CALL ALLOC_GRID(MAXNKS,MAXNRS)
      DO ir = 1, nrs
        DO ik = 1, nks(ir)        
          i1 = imap(ik,ir)
          rs(ik,ir) = knot(i1)%rcen
          zs(ik,ir) = knot(i1)%zcen
          bratio(ik,ir) = knot(i1)%bratio
          id = id + 1
          korpg(ik,ir) = id
          nvertp(id) = knot(i1)%nv
          DO i2 = 1, nvertp(id)
            rvertp(i2,id) = knot(i1)%rv(i2)
            zvertp(i2,id) = knot(i1)%zv(i2)
          ENDDO
c...      Store these in case B2 data from Rhozansky is being loaded:
          divimp_ik(ik,ir) = knot(i1)%ik 
          divimp_ir(ik,ir) = knot(i1)%ir
        ENDDO
      ENDDO

c...  Set NPOLYP:
      npolyp  = id
      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin

c...  Find IKTO,IKTI:
      ikto = grid_load%ikto
      ikti = grid_load%ikti
      IF (ikto.EQ.-1.OR.ikti.EQ.-1)
     .  CALL ER('Readgeneralisedgrid','IKTI and/or IKTO not found',*99)

      CALL OutputData(85,'JUST FINISHED LOADING OSM GRID')


c...  Add virtual rings 1 (core boundary), IRWALL (SOL) and IRTRAP (PFZ):
c      CALL InsertRing(1          ,BEFORE,PERMANENT)
c      CALL InsertRing(nrs-irsep+2,AFTER ,PERMANENT)
c      CALL InsertRing(nrs-irsep+3,BEFORE,PERMANENT)

c     CUTPT1
c     CUTPT2
c     MAXKPTS
c     MAXRINGS
c     CUTRING
c
      IF (opt%f_grid_format.EQ.2) THEN

        IF (ikto.EQ.0) THEN
          STOP 'NOT READY YET, HERE!'
          nopriv = .TRUE.
          CALL InsertRing(1  ,BEFORE,PERMANENT)
          CALL InsertRing(nrs,AFTER ,PERMANENT)
        ELSE
          nopriv = .FALSE.
c...      Add virtual rings 1 (core boundary), IRWALL (SOL) and IRTRAP (PFZ):
          CALL InsertRing(1       ,BEFORE,PERMANENT)
          CALL InsertRing(irwall+0,AFTER ,PERMANENT)
          CALL InsertRing(irwall+1,BEFORE,PERMANENT)
        ENDIF

c        WRITE(0,*) '-------------'
c        WRITE(0,*) '--GRID HACK--'
c        WRITE(0,*) '-------------'
c        irwall = irwall + 1
c        irtrap = irwall + 1
      ELSE
c...    This is fine for SONNET grids, where the number of core and PFZ
c       rings are the same:
        IF (ikto.EQ.0) THEN
          nopriv = .TRUE.
          CALL InsertRing(1  ,BEFORE,PERMANENT)
          CALL InsertRing(nrs,AFTER ,PERMANENT)
        ELSE
          nopriv = .FALSE.
c...      Add virtual rings 1 (core boundary), IRWALL (SOL) and IRTRAP (PFZ):
          CALL InsertRing(1          ,BEFORE,PERMANENT)
          CALL InsertRing(nrs-irsep+2,AFTER ,PERMANENT)
          CALL InsertRing(nrs-irsep+3,BEFORE,PERMANENT)
        ENDIF
      ENDIF

      cutpt1     = ikto
      cutpt2     = ikti             ! These are semi-bogus for a connected double-null...?
      cutring    = irsep - 1
      maxkpts    = nks(irsep)
      maxrings   = irwall
      indexiradj = 1

c      write(0,*) 'cut,max',cutring,maxrings

      IF (.TRUE.) THEN
c        id = 0
c        DO ir = 1, nrs
c          DO ik = 1, nks(ir)        
c            i1 = imap(ik,ir)
c            rs(ik,ir) = knot(i1)%rcen
c            zs(ik,ir) = knot(i1)%zcen
c            bratio(ik,ir) = knot(i1)%bratio
c            id = id + 1
c            korpg(ik,ir) = id
c            nvertp(id) = knot(i1)%nv
c            DO i2 = 1, nvertp(id)
c              rvertp(i2,id) = knot(i1)%rv(i2)
c              zvertp(i2,id) = knot(i1)%zv(i2)
c            ENDDO
c          ENDDO
c        ENDDO
c        CALL SaveSolution
c        CALL OutputData(86,'MAST!')
c        title = '...'
c        desc  = 'Call to STORE from DumpGrid'
c        job   = 'Call to STORE from DumpGrid'
c        equil = 'Call to STORE from DumpGrid'
c        WRITE(0,*) 'CALLING STORE'
c        CALL Store(title,desc,1,job,equil,facta,factb,1,1)
c        WRITE(0,*) 'FUN WITH MAST GRIDS!'
c        STOP 'WHAT?'
      ENDIF

c...  Add virtual boundary cells, which will be stripped off later:
      IF (ctargopt.EQ.0.OR.ctargopt.EQ.1.OR.ctargopt.EQ.2.OR.
     .    ctargopt.EQ.3.OR.ctargopt.EQ.6) 
     .  CALL AddPoloidalBoundaryCells      

c...  Look for PSIn data for full double null grids (code mostly 
c     from tau.d6a):
      DO WHILE (.TRUE.)  
        READ(gridunit,'(A)',END=25) buffer
        IF     (buffer(1:16).EQ.'PSI-DOUBLE-NULLd') THEN   ! direct assignment to each ring, post mortem...
          READ(buffer(17:),*) numpsi
c...      The PSI values are to be loaded in TailorGrid...glorious hack!
          DO i1 = 1, numpsi
            READ(gridunit,*,END=97)
          ENDDO
        ELSEIF (buffer(1:15).EQ.'PSI-DOUBLE-NULL') THEN
          READ(buffer(16:),*) numpsi
c...      The PSI values are listed with one on each line
c         indexed by knot and ring index based on the SONNET 
c         grid coordinates:
          DO i1 = 1, numpsi
            READ(gridunit,*,END=97) ikpsi(i1),irpsi(i1),valpsi(i1)
          ENDDO
c...      Assign to grid rings:
          DO ir = 2, irwall-1
c...        (Need the "-1" because a virtual core ring has been added to the grid)
            ir1 = knot(imap(1,ir-1))%ir 
            DO i1 = 1, numpsi
              IF (irpsi(i1).EQ.ir1) EXIT
            ENDDO
            IF (i1.EQ.numpsi+1) 
     .        CALL ER('Readgeneralisedgrid','Problem with PSIn',*99)
            psitarg(ir,1) = valpsi(i1)
            psitarg(ir,2) = valpsi(i1)
            IF (ir.LT.irsep) THEN
              psitarg(ir-1+irtrap,1) = valpsi(i1)
              psitarg(ir-1+irtrap,2) = valpsi(i1)
            ENDIF
          ENDDO
          WRITE(0,*)
          WRITE(0,*) '------------------------------------------------'
          WRITE(0,*) 'BAD PSITARG -- AND, USING INNER TARGET DATA ONLY'
          WRITE(0,*) '------------------------------------------------'
          WRITE(0,*)
        ENDIF
      ENDDO
 25   REWIND(gridunit)  ! Reset the grid file to the beginning for the PSI:
                        ! and NEUTRAL WALL: data that's loaded in the the
                        ! RAUG subroutine after this routine has been called.   


c...  For consistency with original SONNET code in tau.d6a:
      ir = maxrings
      ik = maxkpts

      DEALLOCATE(knot)
      DEALLOCATE(imap)



      RETURN
 97   CALL ER('ReadGeneralisedGrid_OSM','Unexpected end of file',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: FindGridBreak
c
c Note - there are no virtual rings at this point -- should add them first?
c
c
      SUBROUTINE FindGridBreak
      USE mod_grid_divimp
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER i1

c... dicy...
      irbreak = MAXNRS
c...  ...
      DO i1 = 2, grdntseg(1,IKLO)
        IF (grdtseg(i1,1,IKLO).NE.grdtseg(i1-1,1,IKLO)+1) THEN
          irbreak = grdtseg(i1-1,1,IKLO) + 1
          EXIT
        ENDIF
      ENDDO
c...  Search though the target regions and select the first one that
c     does not end on a virtual ring (which is always IRWALL here):
      IF (grdtseg(grdntseg(1,IKLO),1,IKLO)+1.LT.irbreak.AND.  ! * NOT TESTED*
     .    grdntreg(IKLO).GT.2) THEN
        DO i1 = 2, grdntreg(IKLO) 
          IF (grdtseg(1,i1,IKLO).NE.irtrap) THEN 
            irbreak = grdtseg(1,i1,IKLO)
            EXIT
          ENDIF
        ENDDO
      ENDIF
c...  ...
      DO i1 = 2, grdntseg(1,IKHI)
        IF (grdtseg(i1,1,IKHI).NE.grdtseg(i1-1,1,IKHI)+1.AND. 
     .      grdtseg(i1-1,1,IKHI)+1.LT.irbreak) THEN
          irbreak = grdtseg(i1-1,1,IKHI) + 1
          EXIT
        ENDIF
      ENDDO
      IF (grdtseg(grdntseg(1,IKHI),1,IKHI)+1.LT.irbreak.AND.
     .    grdntreg(IKHI).GT.2) THEN
        DO i1 = 2, grdntreg(IKHI) 
c          WRITE(fp,*) '???',i1,grdtseg(1,i1,IKHI),irtrap
          IF (grdtseg(1,i1,IKHI).NE.irtrap) THEN 
            irbreak = grdtseg(1,i1,IKHI)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      IF (irbreak.EQ.MAXNRS) irbreak = 0

c...  Assign NBR:
      IF     (irbreak.EQ.0) THEN
        nbr = 0
      ELSEIF (irbreak.LT.irwall) THEN
        nbr = irwall - irbreak
      ELSE
        nbr = nrs - irbreak + 1
      ENDIF

      RETURN
 99   STOP
      END
