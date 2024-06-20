      subroutine EIRENE_setup_hydkin_reactions(hydkin_default, cadapt)
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_ccona
      use EIRMOD_comprt, only : iunout
      use EIRMOD_comxs
      use EIRMOD_comusr
      use EIRMOD_clgin
      use EIRMOD_cinit
      use EIRMOD_cplot
      use EIRMOD_photon
 
      implicit none
 
      TYPE PROPERTIES
        CHARACTER(8) :: NAME
        INTEGER :: TYP, MASSE, NCHARGE, NPRT, NCHRG, ISRF, ISRT,
     .             ID1, NRC, NFOL, NGEN, NHSTS, ID3, INDEX,
     .             MASS_DB
      END TYPE PROPERTIES
 
      character(len=*), intent(inout) :: hydkin_default, cadapt
      CHARACTER(15), ALLOCATABLE :: HYDSPEC(:)
      CHARACTER(1000) :: HLINE
      CHARACTER(50) :: FILENAME, REAC, RSTRING, LEFT, RIGHT
      CHARACTER(12) :: CHR
      CHARACTER(8) :: FILNAM, PCHR='H+      ', CLPARTS(6), CRPARTS(6),
     .                CHELP, oriname
      CHARACTER(4) :: H123
      CHARACTER(3) :: CRC
      CHARACTER(2) :: COR, CREP
      CHARACTER(10), PARAMETER :: DIGIT   = '0123456789'
      CHARACTER(6), PARAMETER  :: SPECIAL = '_^{}+-'
      CHARACTER(16) :: NONLETTER = SPECIAL//DIGIT
      INTEGER, ALLOCATABLE :: IEIGEN(:), INRC(:), IRC_PART(:)
      integer :: i2, ll, n_reac, n_spec, n_atoms, n_mol, n_ions, ifile,
     .           n_testions, n_bulkions, isp, ipos, iposm, mult, i,
     .           iel, itok, ier, irc, il, ind, ispz, iatm, imol, ig,
     .           iion, ipls, iphot, ndumm, ibulkno, nlparts, nrparts,
     .           ibulk, iscd1, iscd2, iscd3, iscde, iestm, ibgk, k,
     .           ip, ipi, numsec, ispl, iscd4, ir, nre, j, nsc, ic,
     .           iori, ityp, iln
      integer :: massec(4), multsec(4), itypar(0:4)
      real(dp) :: rmn, rmx, eelec, ebulk, escd1, escd2, freac, fldlm,
     .            escd3, e_el, e_k
      logical :: ladapt
      type(properties), allocatable :: species(:)
 
      TVAC=0.02
      DVAC=1.D2
 
      HYDKIN_DEFAULT=ADJUSTL(HYDKIN_DEFAULT)
      I2=INDEX(HYDKIN_DEFAULT,' ')
      LL=LEN_TRIM(HYDKIN_DEFAULT)
      FILENAME=HYDKIN_DEFAULT(1:LL) // '.reactions'
 
      OPEN (UNIT=27,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .      FORM='FORMATTED')
      READ (27,*)
      READ (27,*) CHR,n_reac
      READ (27,*) CHR,n_spec
      READ (27,*) CHR,n_atoms
      READ (27,*) CHR,n_ions
      READ (27,*) CHR,n_mol
 
      ALLOCATE (HYDSPEC(N_SPEC))
      ALLOCATE (IEIGEN(N_SPEC))
      ALLOCATE (INRC(N_SPEC))
      ALLOCATE (IRC_PART(N_SPEC))
 
      READ (27,*) HYDSPEC(1:N_SPEC)
      READ (27,'(A1000)') HLINE
      READ (HLINE(52:),*) IEIGEN(1:N_SPEC)
      READ (27,'(A1000)') HLINE
      READ (HLINE(52:),*) INRC(1:N_SPEC)
 
!      N_BULKIONS = COUNT((SCAN(HYDSPEC(1:N_SPEC),'+') > 0) .AND.
!     .                   (IEIGEN(1:N_SPEC) == 0))
!      N_TESTIONS = N_IONS - N_BULKIONS
 
 
      cor = '  '
      crep = '  '
      ladapt = len_trim(cadapt) > 0
      if (ladapt) then
        ic = index(cadapt,'-->')
        if (ic == 0) then
          ladapt = .false.
        else
          cor = cadapt(1:ic-1)
          crep = cadapt(ic+3:ic+4)
          call EIRENE_replace_string(pchr,cor,crep)
        end if
      end if
 
      allocate (species(n_spec))
 
!  find type of species
      do isp = 1, n_spec
 
        call EIRENE_remove_char (hydspec(isp),species(isp)%name,'_^')
 
        ll=len_trim(hydspec(isp))
 
        if (ladapt) then
          oriname = hydspec(isp)
          iori = EIRENE_find_element(oriname(1:2))
          call EIRENE_replace_string(hydspec(isp),cor,crep)
          call EIRENE_replace_string(species(isp)%name,cor,crep)
        end if
 
 
!  check periodic table of elements
        iel = EIRENE_find_element(species(isp)%name(1:2))
        if (.not.ladapt) iori = iel
 
        if (iel > 0) then
!  species is atom
          species(isp)%typ = 1
          species(isp)%masse   = pte(iel)%el_mass
          species(isp)%mass_db = pte(iori)%el_mass
          species(isp)%ncharge = pte(iel)%el_charge
          species(isp)%nprt    = 1
          species(isp)%nchrg   = 0
          species(isp)%isrf    = 0
          species(isp)%isrt    = 0
          species(isp)%id1     = 0
          species(isp)%nrc     = 0
          species(isp)%nfol    = 0
          species(isp)%ngen    = 0
          species(isp)%nhsts   = 0
          species(isp)%id3     = 0
 
        else
 
          ipos = scan(hydspec(isp),'+-')
          if (ipos > 0) then
!  species is ion
!  assume species is test ion
            species(isp)%typ = 3
!  find charge stage
            iposm = scan(hydspec(isp)(:ipos),'^{',.true.)
            mult = 1
            if (ipos-iposm > 1)
     .        read (hydspec(isp)(iposm+1:ipos-1),*) mult
            species(isp)%nchrg = mult
          else
!  species has to be molecule
            species(isp)%typ = 2
            species(isp)%nchrg = 0
          end if
 
!  now find the ingredients
          if (ladapt)
     .      call EIRENE_find_ingred (oriname, species(isp)%mass_db,
     .                        species(isp)%ncharge, species(isp)%nprt)
 
          call EIRENE_find_ingred (hydspec(isp), species(isp)%masse,
     .                      species(isp)%ncharge, species(isp)%nprt)
 
          if (.not.ladapt) species(isp)%mass_db = species(isp)%masse
 
          species(isp)%isrf    = 0
          species(isp)%isrt    = 0
          species(isp)%id1     = 0
          species(isp)%nrc     = 0
          species(isp)%nfol    = 0
          species(isp)%ngen    = 0
          species(isp)%nhsts   = 0
          species(isp)%id3     = 0
        end if
 
!  species is bulk ion
        if (ieigen(isp) == 0) species(isp)%typ = 4
 
!  check if species was already defined in eirene input and overrule type
!  if necessary
 
        ityp = 0
        itypar = (/ nsph, nspa, nspam, nspami, nsptot /)
        do ispz = 1, nsptot
          do while (ispz > itypar(ityp))
            ityp = ityp + 1
          end do
          if (texts(ispz) == species(isp)%name) exit
        end do
 
        if (ispz <= nsptot) species(isp)%typ = ityp
 
      end do
 
 
!  insert species into list of species specified in Eirene input
 
       write (iunout,'(a3,1x,a15,7a8)') 'no','name','short','ityp',
     .        'index','nmass','ncharg','nprt','nchrg'
 
      do isp = 1, n_spec
        select case (species(isp)%typ)
 
!  atom
        case (1)
          ind = natmi+1
          do iatm = 1,natmi
            ispz = nsph + iatm
            if (texts(ispz) == species(isp)%name) ind = iatm
          end do
 
          iatm = ind
          ispz = nsph + iatm
          if (iatm <= natmi) then
!  species found, check specifications
             if ((nmassa(iatm) == species(isp)%masse) .and.
     .           (nchara(iatm) == species(isp)%ncharge).and.
     .           (nprt(ispz) == species(isp)%nprt)) then
               species(isp)%index = iatm
             else
               write (iunout,*) ' SPECIFICATIONS DIFFER FOR ',
     .           texts(nsph+iatm),' SPECIFIED IN INPUT AND ',
     .           species(isp)%name,' FOUND IN ', HYDKIN_DEFAULT
!               write (iunout,*) ' NEW SPECIES ADDED '
!               iatm = natmi + 1
               write (iunout,*) ' PLEASE CORRECT AND RETRY '
               call EIRENE_exit_own(1)
             end if
          end if
          if (iatm > natmi) then
!  species not found, add to list
            nmassa(iatm) = species(isp)%masse
            nchara(iatm) = species(isp)%ncharge
            nrca(iatm) = species(isp)%nrc
            nfola(iatm) = species(isp)%nfol
            ngena(iatm) = species(isp)%ngen
            datd(iatm) = 1._dp
!  shift species
            do ispz = nsptot, nspa+1, -1
              texts(ispz+1) = texts(ispz)
              nprt(ispz+1) = nprt(ispz)
              isrf(ispz+1,1) = isrf(ispz,1)
              isrt(ispz+1,1) = isrt(ispz,1)
              nhsts(ispz+1) = nhsts(ispz)
            end do
            ispz = nsph + iatm
            texts(ispz) = species(isp)%name
            nprt(ispz) = species(isp)%nprt
            isrf(ispz,1) = species(isp)%isrf
            isrt(ispz,1) = species(isp)%isrt
            nhsts(ispz) = species(isp)%nhsts
 
            species(isp)%index = iatm
 
            natmi = natmi + 1
            nspa = nspa + 1
            nspam = nspam + 1
            nspami = nspami + 1
            nsptot = nsptot + 1
          end if
 
!  molecule
        case (2)
          ind = nmoli+1
          do imol = 1,nmoli
            ispz = nspa + imol
            if (texts(ispz) == species(isp)%name) ind = imol
          end do
 
          imol = ind
          ispz = nspa + imol
          if (imol <= nmoli) then
!  species found, check specifications
             if ((nmassm(imol) == species(isp)%masse) .and.
     .           (ncharm(imol) == species(isp)%ncharge) .and.
     .           (nprt(ispz) == species(isp)%nprt)) then
               species(isp)%index = imol
             else
               write (iunout,*) ' SPECIFICATIONS DIFFER FOR ',
     .           texts(ispz),' SPECIFIED IN INPUT AND ',
     .           species(isp)%name,' FOUND IN ', HYDKIN_DEFAULT
!               write (iunout,*) ' NEW SPECIES ADDED '
!               imol = nmoli + 1
               write (iunout,*) ' PLEASE CORRECT AND RETRY '
               call EIRENE_exit_own(1)
             end if
          end if
          if (imol > nmoli) then
!  species not found, add to list
            nmassm(imol) = species(isp)%masse
            ncharm(imol) = species(isp)%ncharge
            nrcm(imol) = species(isp)%nrc
            nfolm(imol) = species(isp)%nfol
            ngenm(imol) = species(isp)%ngen
            dmld(imol) = 1._dp
!  shift species
            do ispz = nsptot, nspam+1, -1
              texts(ispz+1) = texts(ispz)
              nprt(ispz+1) = nprt(ispz)
              isrf(ispz+1,1) = isrf(ispz,1)
              isrt(ispz+1,1) = isrt(ispz,1)
              nhsts(ispz+1) = nhsts(ispz)
            end do
            ispz = nspa + imol
            texts(ispz) = species(isp)%name
            nprt(ispz) = species(isp)%nprt
            isrf(ispz,1) = species(isp)%isrf
            isrt(ispz,1) = species(isp)%isrt
            nhsts(ispz) = species(isp)%nhsts
 
            species(isp)%index = imol
 
            nmoli = nmoli + 1
            nspam = nspam + 1
            nspami = nspami + 1
            nsptot = nsptot + 1
          end if
 
!  test ion
        case (3)
          ind = nioni+1
          do iion = 1,nioni
            ispz = nspam + iion
            if (texts(ispz) == species(isp)%name) ind = iion
          end do
 
          iion = ind
          ispz = nspam + iion
          if (iion <= nioni) then
!  species found, check specifications
             if ((nmassi(iion) == species(isp)%masse) .and.
     .           (nchari(iion) == species(isp)%ncharge) .and.
     .           (nprt(ispz) == species(isp)%nprt) .and.
     .           (nchrgi(iion) == species(isp)%nchrg)) then
               species(isp)%index = iion
             else
               write (iunout,*) ' SPECIFICATIONS DIFFER FOR ',
     .           texts(ispz),' SPECIFIED IN INPUT AND ',
     .           species(isp)%name,' FOUND IN ', HYDKIN_DEFAULT
!               write (iunout,*) ' NEW SPECIES ADDED '
!               iion = nioni + 1
               write (iunout,*) ' PLEASE CORRECT AND RETRY '
               call EIRENE_exit_own(1)
             end if
          end if
          if (iion > nioni) then
!  species not found, add to list
            nmassi(iion) = species(isp)%masse
            nchari(iion) = species(isp)%ncharge
            nchrgi(iion) = species(isp)%nchrg
            nrci(iion) = species(isp)%nrc
            nfoli(iion) = species(isp)%nfol
            ngeni(iion) = species(isp)%ngen
            diod(iion) = 1._dp
!  shift species
            do ispz = nsptot, nspami+1, -1
              texts(ispz+1) = texts(ispz)
              nprt(ispz+1) = nprt(ispz)
              isrf(ispz+1,1) = isrf(ispz,1)
              isrt(ispz+1,1) = isrt(ispz,1)
              nhsts(ispz+1) = nhsts(ispz)
            end do
            ispz = nspam + iion
            texts(ispz) = species(isp)%name
            nprt(ispz) = species(isp)%nprt
            isrf(ispz,1) = species(isp)%isrf
            isrt(ispz,1) = species(isp)%isrt
            nhsts(ispz) = species(isp)%nhsts
 
            species(isp)%index = iion
 
            nioni = nioni + 1
            nspami = nspami + 1
            nsptot = nsptot + 1
          end if
 
!  bulk ion
        case (4)
          ind = nplsi+1
          do ipls = 1,nplsi
            ispz = nspami + ipls
            if (texts(ispz) == species(isp)%name) ind = ipls
          end do
 
          ipls = ind
          ispz = nspami + ipls
          if (ipls <= nplsi) then
!  species found, check specifications
             if ((nmassp(ipls) == species(isp)%masse) .and.
     .           (ncharp(ipls) == species(isp)%ncharge) .and.
     .           (nprt(ispz) == species(isp)%nprt) .and.
     .           (nchrgp(ipls) == species(isp)%nchrg)) then
               species(isp)%index = ipls
             else
               write (iunout,*) ' SPECIFICATIONS DIFFER FOR ',
     .           texts(ispz),' SPECIFIED IN INPUT AND ',
     .           species(isp)%name,' FOUND IN ', HYDKIN_DEFAULT
!               write (iunout,*) ' NEW SPECIES ADDED '
!               ipls = nplsi + 1
               write (iunout,*) ' PLEASE CORRECT AND RETRY '
               call EIRENE_exit_own(1)
             end if
          end if
          if (ipls > nplsi) then
!  species not found, add to list
            nmassp(ipls) = species(isp)%masse
            ncharp(ipls) = species(isp)%ncharge
            nchrgp(ipls) = species(isp)%nchrg
            nrcp(ipls) = species(isp)%nrc
            dpld(ipls) = 1._dp
 
            texts(ispz) = species(isp)%name
            nprt(ispz) = species(isp)%nprt
            isrf(ispz,1) = species(isp)%isrf
            isrt(ispz,1) = species(isp)%isrt
            nhsts(ispz) = species(isp)%nhsts
 
            species(isp)%index = ipls
 
            nplsi = nplsi + 1
            nsptot = nsptot + 1
 
            cdenmodel(ipls) = 'CONSTANT'
 
            ALLOCATE (TDMPAR(IPLS)%TDM)
            TDMPAR(IPLS)%TDM%NRE=1
            ALLOCATE (TDMPAR(IPLS)%TDM%ISP(TDMPAR(IPLS)%TDM%NRE))
            ALLOCATE (TDMPAR(IPLS)%TDM%ITP(TDMPAR(IPLS)%TDM%NRE))
            ALLOCATE (TDMPAR(IPLS)%TDM%ISTR(TDMPAR(IPLS)%TDM%NRE))
            ALLOCATE (TDMPAR(IPLS)%TDM%FNAME(TDMPAR(IPLS)%TDM%NRE))
            ALLOCATE (TDMPAR(IPLS)%TDM%H2(TDMPAR(IPLS)%TDM%NRE))
            ALLOCATE (TDMPAR(IPLS)%TDM%REACTION(TDMPAR(IPLS)%TDM%NRE))
            ALLOCATE (TDMPAR(IPLS)%TDM%CR(TDMPAR(IPLS)%TDM%NRE))
 
            tdmpar(ipls)%tdm%tval = tvac
            tdmpar(ipls)%tdm%dval = dvac
            tdmpar(ipls)%tdm%vxval = 0._dp
            tdmpar(ipls)%tdm%vyval = 0._dp
            tdmpar(ipls)%tdm%vzval = 1._dp
          end if
          if (texts(ispz)(1:2) == pchr(1:2)) ibulkno = ipls
 
        case default
          write (iunout,*) ' SPECIES TYPE ', species(isp)%typ,
     .                     ' FOUND IN ',HYDKIN_DEFAULT
          write (iunout,*) ' SO FAR THIS TYPE IS NOT FORESEEN '
        end select
 
        write (iunout,'(i3,1x,a15,a8,6i8)') isp, hydspec(isp),
     .        species(isp)%name,
     .        species(isp)%typ, species(isp)%index, species(isp)%masse,
     .        species(isp)%ncharge, species(isp)%nprt,
     .        species(isp)%nchrg
      end do
 
      call EIRENE_leer(2)
      FILNAM = 'HYDRTC'
      DO IFILE=1,NDBNAMES
        IF (INDEX(FILNAM,DBHANDLE(IFILE)).NE.0) EXIT
      END DO
 
!pb      il = nreac_lines - n_reac
      il = irlines
      reac_loop: do irc=1, n_reac
        ir = nreaci+irc
        IRC_PART = 0
        READ (27,'(A1000)') HLINE
        READ (HLINE(52:),*) IRC_PART(1:N_SPEC)
        reac = repeat(' ',50)
        ll = len_trim(adjustl(hline(1:52)))
        reac = adjustl(hline(1:52))
        call EIRENE_read_hydkin
     .  (ir,dbfname(ifile),h123,reac,crc,rmn,rmx,
     .                    e_el,e_k,.true.)
 
!  set REACLINES for output of input block 4
        IL = IL + 1
        REACLINES(IL)%NO = ir
        REACLINES(IL)%FILE = FILNAM
        REACLINES(IL)%H_SELECT = H123
        REACLINES(IL)%REAC_STRING = REAC
        REACLINES(IL)%REACTYP = CRC
        REACLINES(IL)%MP = 0
        REACLINES(IL)%MT = 0
        REACLINES(IL)%DPP = 0._dp
        REACLINES(IL)%RMN = rmn
        REACLINES(IL)%RMX = rmx
        REACLINES(IL)%ELEMENT = ' '
        REACLINES(IL)%IZ = 0
        REACLINES(IL)%JFEXMN = 0
        REACLINES(IL)%JFEXMX = 0
        REACLINES(IL)%FP = 0._DP
 
        irlines = il
 
!  setup reaction specifictions for the reacting species
 
        rstring = reacdat(ir)%rtc%hyd%reac_string
        write (iunout,'(i6,1x,a)') ir, rstring
        if (ladapt) then
          call EIRENE_replace_string(rstring,cor,crep)
          write (iunout,'(i6,1x,a/1x)') ir, rstring
        end if
 
        ig = index(rstring,'=')
        left = rstring(1:ig-1)
        right = rstring(ig+1:)
 
        call EIRENE_split (left,'&',nlparts,clparts)
        call EIRENE_split (right,'&',nrparts,crparts)
 
        if (nlparts /= 2) then
          call EIRENE_leer(1)
          write (iunout,*) left(1:len_trim(left)),' = ',right
          write (iunout,*) ' IS NOT A VALID REACTION '
          write (iunout,*) ' THIS REACTION IS NOT USED EIRMOD_'
          cycle reac_loop
        end if
 
!  remove leading multiplication factors
        do i=1,nlparts
          ll = verify(clparts(i),digit)
          if (ll > 0) then
            chelp = repeat(' ',len(chelp))
            chelp = clparts(i)(ll:)
            clparts(i) = repeat(' ',len(clparts(i)))
            clparts(i) = chelp
          end if
        end do
 
        do i=1,nrparts
          ll = verify(crparts(i),digit)
          if (ll > 0) then
            chelp = repeat(' ',len(chelp))
            chelp = crparts(i)(ll:)
            crparts(i) = repeat(' ',len(crparts(i)))
            crparts(i) = chelp
          end if
        end do
 
!  check left hand side of reaction
        if (scan(clparts(1),'pe') > 0) then
          ip = 2
          ipi = 1
        else if (scan(clparts(2),'pe') > 0) then
          ip = 1
          ipi = 2
        else
          write (iunout,*) ' REACTION WITHOUT INVOLVEMENT OF'
          write (iunout,*) ' p OR e IS NOT FORESEEN '
          write (iunout,*) ' THIS REACTION IS NOT USED EIRMOD_'
          cycle reac_loop
        end if
 
!  find the species where this reaction should be added
        ispl = 0
        do i=1,n_spec
          if (species(i)%name == clparts(ip)) then
            ispl = i
            exit
          end if
        end do
 
        if (ispl == 0) then
          write (iunout,*) ' SPECIES NOT KNOWN IN REACTION'
          write (iunout,*) left(1:len_trim(left)), ' = ', right
          write (iunout,*) clparts(ip), ' NOT FOUND IN LIST '
          write (iunout,*) ' THIS REACTION IS NOT USED EIRMOD_'
          cycle reac_loop
        end if
 
        REACLINES(IL)%MT = species(ispl)%mass_db
        MASST(ir) = species(ispl)%mass_db
        DELPOT(ir) = 0._DP
 
!  find precollision bulk particle identifier
        if (index(clparts(ipi),'p') > 0 ) then
          ibulk = ibulkno*100 + 14
          REACLINES(IL)%MP = 1
          MASSP(ir) = 1
        elseif (index(clparts(ipi),'e') > 0 ) then
          ibulk = 115
          REACLINES(IL)%MP = 0
          MASSP(ir) = 0
        else
          ibulk = 0
        end if
 
!  find heavy secondary particle groups
        numsec = 0
        iscd1 = 0
        iscd2 = 0
        iscd3 = 0
        iscd4 = 0
        nsc = 0
        massec = 0
        multsec = 0
        do ip = 1, nrparts
          isp = 0
          do i=1,n_spec
            if (species(i)%name == crparts(ip)) then
              isp = i
              exit
            end if
          end do
          if (isp > 0) then
            numsec = numsec + 1
            if (numsec == 1) then
              iscd1 = species(isp)%typ + irc_part(isp)*10 +
     .                species(isp)%index*100
!              massec(1) = species(isp)%masse
              massec(1) = species(isp)%mass_db
              multsec(1) = irc_part(isp)
            else if (numsec == 2) then
              iscd2 = species(isp)%typ + irc_part(isp)*10 +
     .                species(isp)%index*100
!              massec(2) = species(isp)%masse
              massec(2) = species(isp)%mass_db
              multsec(2) = irc_part(isp)
            else if (numsec == 3) then
              iscd3 = species(isp)%typ + irc_part(isp)*10 +
     .                species(isp)%index*100
!              massec(3) = species(isp)%masse
              massec(3) = species(isp)%mass_db
              multsec(3) = irc_part(isp)
              nsc = 3
            else if (numsec == 4) then
              iscd4 = species(isp)%typ + irc_part(isp)*10 +
     .                species(isp)%index*100
!              massec(4) = species(isp)%masse
              massec(4) = species(isp)%mass_db
              multsec(4) = irc_part(isp)
              nsc = 4
            else
              call EIRENE_leer(1)
              write (iunout,*) ' TOO MANY SECONDARIES FOUND IN REACTION'
              write (iunout,*) left(1:len_trim(left)), ' = ', right
              write (iunout,*) ' THIS REACTION IS NOT USED EIRMOD_'
              cycle reac_loop
            end if
          end if
        end do
 
        reacdat(ir)%nosec = nsc
 
! check mass conservation
        if (massp(ir)+masst(ir) /= sum(massec*multsec)) then
          write (iunout,*)
     .      ' MASS CONSERVATION VIOLATED FOR REACTION ',IR
          write (iunout,*) left(1:len_trim(left)), ' = ', right
        end if
 
 
        iscde = 0
        iestm = 0
        ibgk = 0
 
!  differentiate between CX and PI
        if (crc == 'CX ') then
          if ((all(massec /= massp(ir))) .or.
     .        (all(massec /= masst(ir)))) then
!  PI reaction
            REACLINES(IL)%REACTYP = 'PI '
            ISWR(IR) = 4
          else
!  CX reaction
            iscde = 1000
          end if
        end if
 
 
        eelec = -e_el
        ebulk = 0._dp
        escd1 = e_k
        escd2 = 0._dp
        escd3 = 0._dp
        freac = 0._dp
        fldlm = 0._dp
 
        select case (species(ispl)%typ)
        case (1)
          iatm = species(ispl)%index
          nrca(iatm) = nrca(iatm) + 1
          ireaca(iatm,nrca(iatm)) = ir
          ibulka(iatm,nrca(iatm)) = ibulk
          iscd1a(iatm,nrca(iatm)) = iscd1
          iscd2a(iatm,nrca(iatm)) = iscd2
          iscd3a(iatm,nrca(iatm)) = iscd3
          iscd4a(iatm,nrca(iatm)) = iscd4
          iscdea(iatm,nrca(iatm)) = iscde
          iestma(iatm,nrca(iatm)) = iestm
          ibgka(iatm,nrca(iatm)) = ibgk
          eeleca(iatm,nrca(iatm)) = eelec
          ebulka(iatm,nrca(iatm)) = ebulk
          escd1a(iatm,nrca(iatm)) = escd1
          freaca(iatm,nrca(iatm)) = freac
          fldlma(iatm,nrca(iatm)) = fldlm
        case (2)
          imol = species(ispl)%index
          nrcm(imol) = nrcm(imol) + 1
          ireacm(imol,nrcm(imol)) = ir
          ibulkm(imol,nrcm(imol)) = ibulk
          iscd1m(imol,nrcm(imol)) = iscd1
          iscd2m(imol,nrcm(imol)) = iscd2
          iscd3m(imol,nrcm(imol)) = iscd3
          iscd4m(imol,nrcm(imol)) = iscd4
          iscdem(imol,nrcm(imol)) = iscde
          iestmm(imol,nrcm(imol)) = iestm
          ibgkm(imol,nrcm(imol)) = ibgk
          eelecm(imol,nrcm(imol)) = eelec
          ebulkm(imol,nrcm(imol)) = ebulk
          escd1m(imol,nrcm(imol)) = escd1
          freacm(imol,nrcm(imol)) = freac
          fldlmm(imol,nrcm(imol)) = fldlm
        case (3)
          iion = species(ispl)%index
          nrci(iion) = nrci(iion) + 1
          ireaci(iion,nrci(iion)) = ir
          ibulki(iion,nrci(iion)) = ibulk
          iscd1i(iion,nrci(iion)) = iscd1
          iscd2i(iion,nrci(iion)) = iscd2
          iscd3i(iion,nrci(iion)) = iscd3
          iscd4i(iion,nrci(iion)) = iscd4
          iscdei(iion,nrci(iion)) = iscde
          iestmi(iion,nrci(iion)) = iestm
          ibgki(iion,nrci(iion)) = ibgk
          eeleci(iion,nrci(iion)) = eelec
          ebulki(iion,nrci(iion)) = ebulk
          escd1i(iion,nrci(iion)) = escd1
          freaci(iion,nrci(iion)) = freac
          fldlmi(iion,nrci(iion)) = fldlm
        case (4)
          ipls = species(ispl)%index
          nrcp(ipls) = nrcp(ipls) + 1
          ireacp(ipls,nrcp(ipls)) = ir
          ibulkp(ipls,nrcp(ipls)) = ibulk
          iscd1p(ipls,nrcp(ipls)) = iscd1
          iscd2p(ipls,nrcp(ipls)) = iscd2
          iscd3p(ipls,nrcp(ipls)) = iscd3
          iscd4p(ipls,nrcp(ipls)) = iscd4
          iscdep(ipls,nrcp(ipls)) = iscde
          eelecp(ipls,nrcp(ipls)) = eelec
          ebulkp(ipls,nrcp(ipls)) = ebulk
          escd1p(ipls,nrcp(ipls)) = escd1
          freacp(ipls,nrcp(ipls)) = freac
          fldlmp(ipls,nrcp(ipls)) = fldlm
        case default
          write (iunout,*) ' SPECIES TYPE ', species(isp)%typ,
     .                     ' FOUND IN ',HYDKIN_DEFAULT
          write (iunout,*) ' SO FAR THIS TYPE IS NOT FORESEEN '
          cycle reac_loop
        end select
 
      end do reac_loop
 
      nreaci = nreaci + n_reac
      call EIRENE_leer(2)
 
      DEALLOCATE (HYDSPEC)
      DEALLOCATE (IEIGEN)
      DEALLOCATE (INRC)
      DEALLOCATE (IRC_PART)
 
! check test ions for recombinations
 
      do il = 1, irlines
        if (reaclines(il)%reactyp == 'RC ') then
          ir = reaclines(il)%no
! recombination found, check if used by test ions
          do iion=1,nioni
            do k=1,nrci(iion)
              if (ir == ireaci(iion,k)) then
                if (any(ireaci(iion,k+1:nrci(iion)) == ir) .or.
     .              any(ireaci(iion+1:nioni,:) == ir) .or.
     .              any(ireaca(1:natmi,:) == ir) .or.
     .              any(ireacm(1:nmoli,:) == ir) .or.
     .              any(ireacph(1:nphoti,:) == ir) .or.
     .              any(ireacp(1:nplsi,:) == ir)) then
! reaction ir is used elsewhere ==> copy reaction
                  write (iunout,*) ' problem in setup_hydkin_reactions'
                  write (iunout,*) ' option is still to be written '
                  call EIRENE_exit_own(1)
                else
! reaction ir is only used in conjunction with ion iion
! change reaction type to 'EI '
                  reaclines(il)%reactyp = 'EI '
                  iswr(ir) = 1
                  do iln = il+1, irlines
                    if (reaclines(iln)%no == ir) then
                      reaclines(iln)%reactyp = 'EI '
                    end if
                  end do
                end if
              end if
            end do
          end do
        end if
      end do
 
      CLOSE(27)
 
      NDUMM = 0
 
      open (unit=27,file='block4.'//HYDKIN_DEFAULT)
!pb      write (27,'(a)') '*** 4. '
      WRITE (27,*)
     .  '*      ATOMIC REACTION CARDS, NREACI DATA FIELDS'
      write (27,'(i6)') nreaci
!pb      do il = 1, nreac_lines
      do il = 1, irlines
        ll = max(9,len_trim(REACLINES(IL)%REAC_STRING))
        ir = REACLINES(IL)%NO
        write (27,'(I3,1X,A6,1X,A4,A,1X,A3,2I3,3E12.4)')
     .    REACLINES(IL)%NO, REACLINES(IL)%FILE, REACLINES(IL)%H_SELECT,
     .    REACLINES(IL)%REAC_STRING(1:LL), REACLINES(IL)%REACTYP,
     .    REACLINES(IL)%MP, REACLINES(IL)%MT,
     .    REACLINES(IL)%DPP, REACLINES(IL)%RMN,
     .    REACLINES(IL)%RMX
        if (verify(REACLINES(IL)%ELEMENT,' ') > 0)
     .    write (27,'(4X,A2,1X,I3)')
     .      REACLINES(IL)%ELEMENT, REACLINES(IL)%IZ
        if (index(REACLINES(IL)%FILE,'CONST') > 0) then
          write (27,'(6es12.4)') (REACLINES(IL)%CONST(i),
     .                            I=1,REACLINES(IL)%NCONST)
        end if
        if (index(REACLINES(IL)%FILE,'PHOTON') > 0) then
          write (27,'(12i6)')
     .       reacdat(ir)%phr%line%iprofiletype,
     .       reacdat(ir)%phr%line%ignd,
     .       reacdat(ir)%phr%line%imess, reacdat(ir)%phr%line%ifremd,
     .       reacdat(ir)%phr%line%nrjprt
          do i = 1, reacdat(ir)%phr%line%ifremd
            write(27,'(i6,1x,a2,3x,i6)') i,reacdat(ir)%phr%line%kenn(i),
     .        reacdat(ir)%phr%line%iplsc6(i)
          end do
        end if
        if (REACLINES(IL)%RMN > 0._DP)
     .    write (27,'(I6,6X,5E12.4)')
     .      REACLINES(IL)%JFEXMN,REACLINES(IL)%FP(1:3)
        if (REACLINES(IL)%RMX > 0._DP)
     .    write (27,'(I6,6X,5E12.4)')
     .      REACLINES(IL)%JFEXMX,REACLINES(IL)%FP(4:6)
      end do
 
      write (27,'(a)') '* 4A. atom species cards '
      write (27,'(i6)') natmi
      do iatm = 1, natmi
        ispz = nsph + iatm
        numsec = 0
        if (any(iscd3a(iatm,1:nrca(iatm)) /= 0)) numsec = 3
        if (any(iscd4a(iatm,1:nrca(iatm)) /= 0)) numsec = 4
        write (27,'(I2,1X,A8,12(I3),1X,A10,1X,I2)')
     .    IATM,TEXTS(ISPZ),NMASSA(IATM),NCHARA(IATM),
     .    NDUMM,NDUMM,
     .    ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .    NRCA(IATM),NFOLA(IATM),NGENA(IATM),
     .    NHSTS(ISPZ)
        do k=1,nrca(iatm)
          if (numsec < 3) then
            write (27,'(12i6)')
     .        ireaca(iatm,k),ibulka(iatm,k),iscd1a(iatm,k),
     .        iscd2a(iatm,k),iscdea(iatm,k),iestma(iatm,k),
     .        ibgka(iatm,k)
          else if (numsec == 3) then
            write (27,'(12i6)')
     .        ireaca(iatm,k),ibulka(iatm,k),iscd1a(iatm,k),
     .        iscd2a(iatm,k),iscd3a(iatm,k),iscdea(iatm,k),
     .        iestma(iatm,k),ibgka(iatm,k)
          else if (numsec == 4) then
            write (27,'(12i6)')
     .        ireaca(iatm,k),ibulka(iatm,k),iscd1a(iatm,k),
     .        iscd2a(iatm,k),iscd3a(iatm,k),iscd4a(iatm,k),
     .        iscdea(iatm,k),iestma(iatm,k),ibgka(iatm,k)
          end if
          write (27,'(6es12.4)')
     .      eeleca(iatm,k),ebulka(iatm,k),escd1a(iatm,k),
     .      escd2,freaca(iatm,k),fldlma(iatm,k)
        end do
      end do
 
      write (27,'(a)') '* 4B. molecule species cards '
      write (27,'(i6)') nmoli
      do imol = 1, nmoli
        ispz = nspa + imol
        numsec = 0
        if (any(iscd3m(imol,1:nrcm(imol)) /= 0)) numsec = 3
        if (any(iscd4m(imol,1:nrcm(imol)) /= 0)) numsec = 4
        write (27,'(I2,1X,A8,12(I3),1X,A10,1X,I2)')
     .    IMOL,TEXTS(ISPZ),NMASSM(IMOL),NCHARM(IMOL),
     .    NPRT(ISPZ),NDUMM,
     .    ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .    NRCM(IMOL),NFOLM(IMOL),NGENM(IMOL),
     .    NHSTS(ISPZ)
        do k=1,nrcm(imol)
          if (numsec < 3) then
            write (27,'(12i6)')
     .        ireacm(imol,k),ibulkm(imol,k),iscd1m(imol,k),
     .        iscd2m(imol,k),iscdem(imol,k),iestmm(imol,k),
     .        ibgkm(imol,k)
          else if (numsec == 3) then
            write (27,'(12i6)')
     .        ireacm(imol,k),ibulkm(imol,k),iscd1m(imol,k),
     .        iscd2m(imol,k),iscd3m(imol,k),iscdem(imol,k),
     .        iestmm(imol,k),ibgkm(imol,k)
          else if (numsec == 4) then
            write (27,'(12i6)')
     .        ireacm(imol,k),ibulkm(imol,k),iscd1m(imol,k),
     .        iscd2m(imol,k),iscd3m(imol,k),iscd4m(imol,k),
     .        iscdem(imol,k),iestmm(imol,k),ibgkm(imol,k)
          end if
          write (27,'(6es12.4)')
     .      eelecm(imol,k),ebulkm(imol,k),escd1m(imol,k),
     .      escd2,freacm(imol,k)
        end do
      end do
 
      write (27,'(a)') '* 4C. test ion species cards '
      write (27,'(i6)') nioni
      do iion = 1, nioni
        ispz = nspam + iion
        numsec = 0
        if (any(iscd3i(iion,1:nrci(iion)) /= 0)) numsec = 3
        if (any(iscd4i(iion,1:nrci(iion)) /= 0)) numsec = 4
        write (27,'(I2,1X,A8,12(I3),1X,A10,1X,I2)')
     .    IION,TEXTS(ISPZ),NMASSI(IION),NCHARI(IION),
     .    NPRT(ISPZ),NCHRGI(IION),
     .    ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .    NRCI(IION),NFOLI(IION),NGENI(IION),
     .    NHSTS(ISPZ)
        do k=1,nrci(iion)
          if (numsec < 3) then
            write (27,'(12i6)')
     .        ireaci(iion,k),ibulki(iion,k),iscd1i(iion,k),
     .        iscd2i(iion,k),iscdei(iion,k),iestmi(iion,k),
     .        ibgki(iion,k)
          else if (numsec == 3) then
            write (27,'(12i6)')
     .        ireaci(iion,k),ibulki(iion,k),iscd1i(iion,k),
     .        iscd2i(iion,k),iscd3i(iion,k),iscdei(iion,k),
     .        iestmi(iion,k),ibgki(iion,k)
          else if (numsec == 4) then
            write (27,'(12i6)')
     .        ireaci(iion,k),ibulki(iion,k),iscd1i(iion,k),
     .        iscd2i(iion,k),iscd3i(iion,k),iscd4i(iion,k),
     .        iscdei(iion,k),iestmi(iion,k),ibgki(iion,k)
          end if
          write (27,'(6es12.4)')
     .      eeleci(iion,k),ebulki(iion,k),escd1i(iion,k),
     .      escd2,freaci(iion,k)
        end do
      end do
 
      write (27,'(a)') '* 4D. photon species cards '
      write (27,'(i6)') nphoti
      do iphot = 1, nphoti
        ispz = iphot
        numsec = 0
        if (any(iscd3ph(iphot,1:nrcph(iphot)) /= 0)) numsec = 3
        if (any(iscd4ph(iphot,1:nrcph(iphot)) /= 0)) numsec = 4
        write (27,'(I2,1X,A8,12(I3),1X,A10,1X,I2)')
     .    IPHOT,TEXTS(ISPZ),NDUMM,NDUMM,
     .    NDUMM,NDUMM,
     .    ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .    NRCPH(IPHOT),NFOLPH(IPHOT),NGENPH(IPHOT),
     .    NHSTS(ISPZ)
        do k=1,nrcph(iphot)
          if (numsec < 3) then
            write (27,'(12i6)')
     .        ireacph(iphot,k),ibulkph(iphot,k),iscd1ph(iphot,k),
     .        iscd2ph(iphot,k),iscdeph(iphot,k),iestmph(iphot,k),
     .        ibgkph(iphot,k)
          else if (numsec == 3) then
            write (27,'(12i6)')
     .        ireacph(iphot,k),ibulkph(iphot,k),iscd1ph(iphot,k),
     .        iscd2ph(iphot,k),iscd3ph(iphot,k),iscdeph(iphot,k),
     .        iestmph(iphot,k),ibgkph(iphot,k)
          else if (numsec == 4) then
            write (27,'(12i6)')
     .        ireacph(iphot,k),ibulkph(iphot,k),iscd1ph(iphot,k),
     .        iscd2ph(iphot,k),iscd3ph(iphot,k),iscd4ph(iphot,k),
     .        iscdeph(iphot,k),iestmph(iphot,k),ibgkph(iphot,k)
          end if
          write (27,'(6es12.4)')
     .      eelecph(iphot,k),ebulkph(iphot,k),escd1ph(iphot,k),
     .      escd2,freacph(iphot,k),fldlmph(iphot,k)
        end do
      end do
 
      write (27,'(a)') '*** 5. bulk ion species cards '
      WRITE (27,'(a)') '*5A.   BULK ION SPECIES CARDS, NPLSI SPECIES '
      write (27,'(i6)') nplsi
      do ipls = 1, nplsi
        ispz = nspami + ipls
        numsec = 0
        if (any(iscd3p(ipls,1:nrcp(ipls)) /= 0)) numsec = 3
        if (any(iscd4p(ipls,1:nrcp(ipls)) /= 0)) numsec = 4
        nre = 0
        if (LEN_TRIM(CDENMODEL(IPLS)) > 0) then
          nre = TDMPAR(IPLS)%TDM%NRE
        end if
        write (27,'(I2,1X,A8,12(I3),1X,A10,1X,I2)')
     .    IPLS,TEXTS(ISPZ),NMASSP(IPLS),NCHARP(IPLS),
     .    NPRT(ISPZ),NCHRGP(IPLS),
     .    ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .    NRCP(IPLS),NDUMM,NDUMM,
     .    NHSTS(ISPZ),NDUMM,CDENMODEL(IPLS),NRE
        do k=1,nrcp(ipls)
          if (numsec < 3) then
            write (27,'(12i6)')
     .        ireacp(ipls,k),ibulkp(ipls,k),iscd1p(ipls,k),
     .        iscd2p(ipls,k),iscdep(ipls,k)
          else if (numsec == 3) then
            write (27,'(12i6)')
     .        ireacp(ipls,k),ibulkp(ipls,k),iscd1p(ipls,k),
     .        iscd2p(ipls,k),iscd3p(ipls,k),iscdep(ipls,k)
          else if (numsec == 4) then
            write (27,'(12i6)')
     .        ireacp(ipls,k),ibulkp(ipls,k),iscd1p(ipls,k),
     .        iscd2p(ipls,k),iscd3p(ipls,k),iscd4p(ipls,k),
     .        iscdep(ipls,k)
          end if
          write (27,'(6es12.4)')
     .      eelecp(ipls,k),ebulkp(ipls,k),escd1p(ipls,k),
     .      escd2,freacp(ipls,k)
        end do
 
        SELECT CASE (CDENMODEL(IPLS))
        CASE ('FORT.13   ')
          WRITE (27,'(12i6)') TDMPAR(IPLS)%TDM%ISP(1)
        CASE ('FORT.10   ')
          WRITE (27,'(12i6)') TDMPAR(IPLS)%TDM%ISP(1),
     .                        TDMPAR(IPLS)%TDM%ITP(1),
     .                        TDMPAR(IPLS)%TDM%ISTR(1)
        CASE ('CONSTANT  ')
          WRITE (27,'(6es12.4)') TDMPAR(IPLS)%TDM%TVAL,
     .                           TDMPAR(IPLS)%TDM%DVAL,
     .                           TDMPAR(IPLS)%TDM%VXVAL,
     .                           TDMPAR(IPLS)%TDM%VYVAL,
     .                           TDMPAR(IPLS)%TDM%VZVAL
        CASE ('MULTIPLY  ')
          WRITE (27,'(3I6,6x,3E12.4)')
     .           TDMPAR(IPLS)%TDM%ISP(1),
     .           TDMPAR(IPLS)%TDM%ITP(1),
     .           TDMPAR(IPLS)%TDM%ISTR(1),
     .           TDMPAR(IPLS)%TDM%DFACTOR,
     .           TDMPAR(IPLS)%TDM%TFACTOR,
     .           TDMPAR(IPLS)%TDM%VFACTOR
        CASE ('SAHA      ')
!PB   TO BE WRITTEN
        CASE ('BOLTZMANN ')
          WRITE (27,'(3I6,6x,2E12.4)')
     .           TDMPAR(IPLS)%TDM%ISP(1),
     .           TDMPAR(IPLS)%TDM%ITP(1),
     .           TDMPAR(IPLS)%TDM%ISTR(1),
     .           TDMPAR(IPLS)%TDM%G_BOLTZ,
     .           TDMPAR(IPLS)%TDM%DELTAE
        CASE ('CORONA    ')
          WRITE (27,'(3I6,1X,A6,1X,A4,A9,A3,E12.4)')
     .           TDMPAR(IPLS)%TDM%ISP(1),
     .           TDMPAR(IPLS)%TDM%ITP(1),
     .           TDMPAR(IPLS)%TDM%ISTR(1),
     .           TDMPAR(IPLS)%TDM%FNAME(1),
     .           TDMPAR(IPLS)%TDM%H2(1),
     .           TDMPAR(IPLS)%TDM%REACTION(1),
     .           TDMPAR(IPLS)%TDM%CR(1),
     .           TDMPAR(IPLS)%TDM%A_CORONA
        CASE ('COLRAD    ')
          DO I=1, TDMPAR(IPLS)%TDM%NRE
            WRITE (27,'(3I6,1X,A6,1X,A4,A9,A3)')
     .             TDMPAR(IPLS)%TDM%ISP(I),
     .             TDMPAR(IPLS)%TDM%ITP(I),
     .             TDMPAR(IPLS)%TDM%ISTR(I),
     .             TDMPAR(IPLS)%TDM%FNAME(I),
     .             TDMPAR(IPLS)%TDM%H2(I),
     .             TDMPAR(IPLS)%TDM%REACTION(I),
     .             TDMPAR(IPLS)%TDM%CR(I)
          END DO
        CASE DEFAULT
!PB  NOTHING TO BE DONE
        END SELECT
       end do
 
      WRITE (27,'(a)') '*5B.   PLASMA BACKGROUND DATA '
      WRITE (27,'(12i6)') (INDPRO(J),J=1,12)
      IF (INDPRO(1).LE.5.AND.NPLSI.GT.0)
     .  WRITE (27,'(6es12.4)') TE0,TE1,TE2,TE3,TE4,TE5
      IF (INDPRO(2).LE.5.AND.NPLSI.GT.0) THEN
        IF (NPLSTI == 1) THEN
          WRITE (27,'(6es12.4)') TI0(1),TI1(1),TI2(1),
     .                           TI3(1),TI4(1),TI5(1)
        ELSE
          WRITE (27,'(6es12.4)') (TI0(I),TI1(I),TI2(I),
     .                            TI3(I),TI4(I),TI5(I),I=1,NPLSI)
        END IF
      END IF
      IF (INDPRO(3).LE.5)
     .  WRITE (27,'(6es12.4)') (DI0(I),DI1(I),DI2(I),
     .                          DI3(I),DI4(I),DI5(I),I=1,NPLSI)
      IF (INDPRO(4).LE.5) THEN
        IF (NPLSV == 1) THEN
          WRITE (27,'(6es12.4)') VX0(1),VX1(1),VX2(1),
     .                           VX3(1),VX4(1),VX5(1)
          WRITE (27,'(6es12.4)') VY0(1),VY1(1),VY2(1),
     .                           VY3(1),VY4(1),VY5(1)
          WRITE (27,'(6es12.4)') VZ0(1),VZ1(1),VZ2(1),
     .                           VZ3(1),VZ4(1),VZ5(1)
        ELSE
          WRITE (27,'(6es12.4)') (VX0(I),VX1(I),VX2(I),
     .                            VX3(I),VX4(I),VX5(I),I=1,NPLSI)
          WRITE (27,'(6es12.4)') (VY0(I),VY1(I),VY2(I),
     .                            VY3(I),VY4(I),VY5(I),I=1,NPLSI)
          WRITE (27,'(6es12.4)') (VZ0(I),VZ1(I),VZ2(I),
     .                            VZ3(I),VZ4(I),VZ5(I),I=1,NPLSI)
        ENDIF
      ENDIF
      IF (INDPRO(5).LE.5)
     .  WRITE (27,'(6es12.4)') B0,B1,B2,B3,B4,B5
      IF (INDPRO(12).LE.5)
     .  WRITE (27,'(6es12.4)') VL0,VL1,VL2,VL3,VL4,VL5
 
      close (unit=27)
 
      return
 
      contains
 
 
      subroutine EIRENE_replace_string (inchar,rem,rep)
      character(len=*), intent(inout) :: inchar
      character(len=*), intent(in) :: rem, rep
      character(len=2*len(inchar)) :: outchar
      integer :: lin, lout, lrem, lrep, i, io, ii, linc
 
      linc = len(inchar)
 
      lout = len(outchar)
      outchar = repeat(' ',lout)
 
      lin=len_trim(inchar)
      lrem=len_trim(rem)
      lrep=len_trim(rep)
 
      io = 0
      ii = 0
 
      i = index(inchar(ii+1:lin),rem(1:lrem)) - 1
 
      do while (i >= 0)
 
        if (io+i > lout) exit
        outchar(io+1:io+i) = inchar(ii+1:ii+i)
        io = io + i
        if (io+lrep > lout) exit
        outchar(io+1:io+lrep) = rep(1:lrep)
        io = io + lrep
 
        ii = ii + i + lrem
 
        i = index(inchar(ii+1:lin),rem(1:lrem)) - 1
 
      end do
 
      if (i > = 0) then
        write (iunout,*) ' ERROR IN REPLACE_STRING '
        write (iunout,*) ' STRING IS TOO SHORT TO HOLD ALL REPLACEMENTS'
        write (iunout,*) ' INCHAR = ',inchar
        write (iunout,*) ' REMCHAR = ',rem
        write (iunout,*) ' REPCHAR = ',rep
        write (iunout,*) ' STRING SHORTENED TO ',outchar(1:linc)
        inchar(1:linc) = outchar(1:linc)
        return
      end if
 
      outchar(io+1:io+lin-ii) = inchar(ii+1:lin)
      io = io + lin-ii
 
      if (io > linc) then
        write (iunout,*) ' ERROR IN REPLACE_STRING '
        write (iunout,*) ' STRING IS TOO SHORT TO HOLD ALL REPLACEMENTS'
        write (iunout,*) ' INCHAR = ',inchar
        write (iunout,*) ' REMCHAR = ',rem
        write (iunout,*) ' REPCHAR = ',rep
        write (iunout,*) ' STRING SHORTENED TO ',outchar(1:linc)
      end if
 
      inchar = outchar(1:linc)
 
      return
      end subroutine EIRENE_replace_string
 
 
 
      subroutine EIRENE_find_ingred (name,mass,icharge,ipart)
 
      character(len=*), intent(in) :: name
      integer, intent(out) :: mass, icharge, ipart
 
      CHARACTER(26), PARAMETER :: LETTER  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      CHARACTER(26), PARAMETER :: SLETTER = 'abcdefghijklmnopqrstuvwxyz'
      CHARACTER(10), PARAMETER :: DIGIT   = '0123456789'
      CHARACTER(6), PARAMETER  :: SPECIAL = '_^{}+-'
      CHARACTER(16) :: NONLETTER = SPECIAL//DIGIT
      CHARACTER(12) :: CHR
      integer :: ianf, ipos, il, it, imult, il2, ll
 
      mass = 0
      icharge = 0
      ipart = 0
 
      ianf = 1
      ll = len(name)
 
      do
 
        ipos = scan(name(ianf:),letter)
        if (ipos == 0) exit
 
        il = ianf + ipos -1
        it = il + 1
 
        chr = '  '
        chr(1:1) = name(il:il)
 
        if (index(sletter,name(it:it)) > 0) then
          chr(2:2) = name(it:it)
          it = it + 1
        end if
 
        iel = EIRENE_find_element(chr(1:2))
 
        if (name(it:it) == '_') it = it + 1
 
        imult = 1
 
        il2 = 0
        if (index(digit,name(it:it)) > 0) then
          il2 = scan(name(it+1:),letter//special)
          if (il2 <= 0) il2 = ll-it+1
          read (name(it:it+il2-1),*) imult
        end if
 
        mass = mass + imult*pte(iel)%el_mass
        icharge = icharge + imult*pte(iel)%el_charge
        ipart = ipart + imult
 
        ianf = it+max(il2-1,0)
 
      end do
 
      return
      end subroutine EIRENE_find_ingred
 
 
      subroutine EIRENE_split(string, splchr, nparts, parts)
      character(len=*), intent(in) :: string, splchr
      character(len=*), intent(out) :: parts(6)
      integer, intent(out) :: nparts
      integer :: ianf, ipos, ll
 
      nparts = 0
      parts = repeat(' ',len(parts(1)))
 
      ianf = 1
      ll=len_trim(string)
 
      do
 
        ipos = scan(string(ianf:),splchr)
 
        if (ipos == 0) exit
 
        nparts = nparts + 1
 
        if (nparts > 5) then
          write (iunout,*) ' TOO MANY PARTS IN STRING '
          write (iunout,*) ' STRING = ''',string,''''
          write (iunout,*) ' SPLITTING ABANDONNED '
          nparts = nparts - 1
          exit
        end if
 
        parts(nparts) = adjustl(string(ianf:ianf+ipos-2))
        ianf = ianf + ipos
 
      end do
 
      nparts = nparts + 1
      parts(nparts) = adjustl(string(ianf:ll))
 
      return
      end subroutine EIRENE_split
 
      end subroutine EIRENE_setup_hydkin_reactions
