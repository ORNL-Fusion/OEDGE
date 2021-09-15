module lim_netcdf
  use error_handling
  use debug_options
  use nc_utils_generic


contains

  !
  !     Netcdf output
  !
  subroutine write_netcdf_output(TITLE,NIZS,NOUT,IERR,JOB,IMODE,PLAMS,PIZS,NLS,FACTA,FACTB,ITER,NITERS)   

    use mod_params
    use mod_cadas2
    use mod_dynam1
    use mod_dynam3
    use mod_comt2
    use mod_comnet
    use mod_cnoco
    use mod_comtor
    use mod_cadas
    use mod_comtau
    use mod_comxyt
    use mod_coords
    use yreflection


    !     
    !  *********************************************************************        
    !  *                                                                   *        
    !  *  DUMP:  STORE RESULTS OF LIM RUN IN UNFORMATTED FILE "NOUT".      *        
    !  *         THIS FILE SHOULD BE REWOUND BEFORE DUMP IS CALLED FOR THE *        
    !  *         FIRST TIME.                                               *        
    !  *                                                                   *        
    !  *********************************************************************        
    !                                                                               
    IMPLICIT  none

    CHARACTER TITLE*(*),JOB*(*)
    INTEGER   NIZS,IMODE,NLS                                                  
    REAL      PLAMS(MAXNLS),FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)                 
    INTEGER   NOUT,PIZS(MAXNLS),ITER,NITERS,JY                      

    INTEGER   IOS,JBLOCK,IL,II,IP,J,KBLOCK,IT,IO                              
    INTEGER   IQX,IX,IY,IZ,IYB,IYE,IZS,IZE,IQS,IQE                     
    !integer,parameter :: iblock = 1540
    integer :: ierr

    !
    ! jdemod - initialization resulting from complaints by INTEL compiler
    ! 
    ios = 0
    !



    !
    ! jdemod - open the netcdf file for output
    !
    ierr=open_nc_file('lim_netcdf_out.nc',NC_WRITE)
    if (ierr.ne.0) then 
       call errmsg('WRITE_NETCDF_OUTPUT: ERROR OPENING OUTPUT FILE:',ierr)
       return
    endif

    !                                                                               
    !-----------------------------------------------------------------------        
    !     WRITE COMMONS, LINE DATA, SHORT ARRAYS, PARAMETERS, ETC                   
    !-----------------------------------------------------------------------        
    !
    write(0,*) 'Writing to file: |      |'
    ierr = write_nc('VERSION',verson,'Code version')

    ierr = write_nc('TITLE',title,'Case Title')
    ierr = write_nc('JOB',job,'Job Description')

    !
    ! Write out LIM parameters - these were not previously in the LIM raw file
    ! 
    !  PARAMETER (MAXNXS=100,  MAXNYS=500, MAXIZS=74,   MAXQYS=5000,   &            
    !   MAXQXS=500,  MAXNLS=8,   MAXNTS=1,    MAXIMP=1000000,&               
    !   MAXY3D=500,  MAXNPS=31,  MAXOS =500,  VERSON='L3/04',&           
    !   MAXINS=100,   ISECT =128, MAXPUT=1000, BIG=.TRUE.,   &
    !   MAXLPD=20,   MAXT=10,    MAXLEN=100,  MAXADS=60,     &

	write(0,*) 'Writing to file: |*     |'
    ierr = write_nc('MAXNXS',maxnxs,'Maximum Number of X-bin boundaries')
    ierr = write_nc('MAXNYS',maxnys,'Maximum Number of Y-bin boundaries')
    ierr = write_nc('MAXIZS',maxizs,'Maximum Number of impurity charge states')
    ierr = write_nc('MAXQXS',maxqxs,'Maximum Number of X-bins for plasma background')
    ierr = write_nc('MAXQYS',maxqys,'Maximum Number of Y-bins for plasma background')

    ierr = write_nc('MAXNLS',maxnls,'Maximum number of specific emission line profiles')
    ierr = write_nc('MAXNTS',maxnts,'Maximum number of time steps')
    ierr = write_nc('MAXIMP',maximp,'Maximum Number of impurity particles to follow')

    ierr = write_nc('MAXY3D',maxy3d,'Maximum Number of Y-bins for 3D')
    ierr = write_nc('MAXNPS',maxnps,'Maximum Number of Poloidal bins (3D)')
    ierr = write_nc('MAXOS',maxos,'Maximum Number of limiter surface elements')
    ierr = write_nc('MAXINS',maxins,'Maximum Number of elements in input arrays')

    ierr = write_nc('ISECT',isect,'Chunking step for random number arrays')

    ierr = write_nc('MAXPUT',maxput,'Maximum Split/Roulette particles')
    ierr = write_nc('MAXLPD',maxlpd,'Max size of the LPDION array')
    ierr = write_nc('MAXT',maxt,'Max size of debug steps')
    ierr = write_nc('MAXLEN',maxlen,'Max particle track length for debugging')
    ierr = write_nc('MAXADS',maxads,'Max size of ADAS coefficient data array')


    !
    ! Actually used sizes
    !
    write(0,*) 'Writing to file: |**    |'
    ierr = write_nc('NY3D',NY3D  ,'Number of 3D Y-bins in use')
    ierr = write_nc('ITER',ITER  ,'LIM Iteration count')
    ierr = write_nc('NITERS',NITERS,'Number of LIM Iterations')
    ierr = write_nc('NXS',NXS   ,'Number of X-bins')
    ierr = write_nc('NYS',NYS   ,'Number of Y-bins')
    ierr = write_nc('NQXSO',NQXSO ,'Number of radial elements outboard')
    ierr = write_nc('NQXSI',NQXSI ,'Number of radial elements inboard')
    ierr = write_nc('NQYS',NQYS  ,'Number of Y positions for EY and VHY')
    ierr = write_nc('NTS',NTS   ,'Number of time bins')
    ierr = write_nc('NIZS',NIZS  ,'Number of impurity charge states')
    ierr = write_nc('NLS',NLS   ,'Number of particular spectroscopic lines')
    ierr = write_nc('IMODE',IMODE ,'LIM simulation mode')



    !
    !     Save the scaling factor for the case - the scaling factor is
    !     kept in DEFACT - however, absfac was added for compatibility 
    !     with DIVIMP code. ABSFAC was added to the comtor common include
    !     file. DEFACT is assigned to ABSFAC in LIM3.L3I
    !

    ierr = write_nc('ABSFAC',absfac     ,'Absolute scaling factor for case results')

    !                                                                               
    !---- SOME OF COMTOR COMMON, SOME OF COMXYT COMMON, MISCELLANEOUS               
    !                                                                                      

    ierr = write_nc('CA',CA     ,'Radial position of plamsa center','m')
    ierr = write_nc('CAW',CAW    ,'Radial position of wall','m')
    ierr = write_nc('CL',CL     ,'Connection length (1/2 distance between limiter surfaces)','m')
    ierr = write_nc('CRMB',CRMB   ,'Mass of background plasma ions','amu')
    ierr = write_nc('CIZB',CIZB   ,'Background ion charge')
    ierr = write_nc('CRMI',CRMI   ,'Impurity ion mass','amu')
    ierr = write_nc('CISEED',CISEED ,'Random number seed used for case')

    !
    ! LIM writes a number of quantites to the raw file that have very limited applicability
    ! - These will not be written to the netcdf file at the present time. The netcdf file
    !   will mostly contain results rather than a collection of input values to LIM
    ! - It is also important to note that LIM appears to be very inconsistent in terms
    !   of which input parameters have been written to the raw file .. so for the time
    !   being all of these will generally be left out. 
    !
    !ierr = write_nc('CTBIN',CTBIN  ,'Inboard plasma temperature base value','eV')
    !ierr = write_nc('CGTIN1',CGTIN1 ,'Inboard plasma temperature linear increase','eV/m')
    !ierr = write_nc('CLTIN1',CLTIN1 ,'Inboard plasma temperature exponential increase','m')
    !ierr = write_nc('CTBOUL',CTBOUL ,'Outboard plasma temperature base value Y<0','eV')
    !ierr = write_nc('CLTOUL',CLTOUL ,'Outboard plasma temperature exponential decay','m'
    !ierr = write_nc('CATIN',CATIN  ,'')

    !ierr = write_nc('CNBIN',CNBIN  ,'Inboard plasma density base value','m-3')
    !ierr = write_nc('CGNIN1',CGNIN1 ,'Inboard plasma density linear decay 1','m-4')
    !ierr = write_nc('CLNIN1',CLNIN1 ,'Inboard plasma density exponental decay 1','m')
    !ierr = write_nc('CNBOUL',CNBOUL ,'Outboard plasma density base Y<0','m-3')
    !ierr = write_nc('CLNOUL',CLNOUL ,'Outboard plasma density exponential decay Y<0','m')
    !ierr = write_nc('CANIN',CANIN  ,'')

    !ierr = write_nc('CTIBIN',CTIBIN ,'')
    !ierr = write_nc('CGTIIN1',CGTIIN1,'')
    !ierr = write_nc('CLTIIN1',CLTIIN1,'')
    !ierr = write_nc('CTIBOUL',CTIBOUL,'')
    !ierr = write_nc('CLTIOUL',CLTIOUL,'')
    !ierr = write_nc('CVIN',CVIN   ,'Inboard pinch parameter')
    !ierr = write_nc('CYFAR',CYFAR  ,'')
    !ierr = write_nc('CRDD',CRDD   ,'Inboard diffusion decay rate','m2/s')
    !ierr = write_nc('CXSC',CXSC   ,'X injection position','m')
    !ierr = write_nc('CYSC',CYSC   ,'Y injection position','m')
    !ierr = write_nc('CFIMP',CFIMP  ,'')
    !ierr = write_nc('CONO',CONO   ,'')
    !ierr = write_nc('CTEMSC',CTEMSC ,'')
    !ierr = write_nc('CTIMSC',CTIMSC ,'')
    !ierr = write_nc('CEBD',CEBD   ,'')
    !ierr = write_nc('CTGAS',CTGAS  ,'')
    !ierr = write_nc('CEMAXF',CEMAXF ,'')
    !ierr = write_nc('CENGSC',CENGSC ,'')
    !ierr = write_nc('CION',CION   ,'')
    !ierr = write_nc('CIZSC',CIZSC  ,'')
    !ierr = write_nc('CTRESH',CTRESH ,'')
    !ierr = write_nc('CIZEFF',CIZEFF ,'')
    !ierr = write_nc('CNHC',CNHC   ,'')
    !ierr = write_nc('CNHO',CNHO   ,'')
    !ierr = write_nc('CLAMHX',CLAMHX ,'')
    !ierr = write_nc('CLAMHY',CLAMHY ,'')
    !ierr = write_nc('CXNEAR',CXNEAR ,'')
    !ierr = write_nc('CYNEAR',CYNEAR ,'')
    !ierr = write_nc('CKO',CKO    ,'')
    !ierr = write_nc('CHALFL',CHALFL ,'')
    !ierr = write_nc('CKI',CKI    ,'')
    !ierr = write_nc('CCUT',CCUT   ,'')
    !ierr = write_nc('CPFIR',CPFIR  ,'')
    !ierr = write_nc('CPSUB',CPSUB  ,'')
    !ierr = write_nc('CPSC',CPSC   ,'')
    !ierr = write_nc('CSNORM',CSNORM ,'')
    !ierr = write_nc('CDPOL',CDPOL  ,'')
    !ierr = write_nc('CPLSMA',CPLSMA ,'')
    !ierr = write_nc('CEYIN',CEYIN  ,'')
    !ierr = write_nc('CVHYIN',CVHYIN ,'')
    !ierr = write_nc('CSTEPN',CSTEPN ,'')
    !ierr = write_nc('CIZSET',CIZSET ,'')
    !ierr = write_nc('CZENH',CZENH  ,'')
    !ierr = write_nc('CEYOUT',CEYOUT ,'')
    !ierr = write_nc('CVHOUT',CVHOUT ,'')
    !ierr = write_nc('CDIFOP',CDIFOP ,'')
    !ierr = write_nc('CTWOL',CTWOL  ,'')
    !ierr = write_nc('CYSTAG',CYSTAG ,'')
    !ierr = write_nc('CONI',CONI   ,'')
    !ierr = write_nc('XSCALO',XSCALO ,'')
    !ierr = write_nc('XSCALI',XSCALI ,'')
    !ierr = write_nc('YSCALE',YSCALE ,'')
    !ierr = write_nc('CANAL',CANAL  ,'')
    !ierr = write_nc('CTHETB',CTHETB ,'')
    !ierr = write_nc('CSINTB',CSINTB ,'')
    !ierr = write_nc('CLFACT',CLFACT ,'')
    !ierr = write_nc('CFBGFF',CFBGFF ,'')
    !ierr = write_nc('CLNIN2',CLNIN2 ,'')
    !ierr = write_nc('CLTIIN2',CLTIIN2,'')
    !ierr = write_nc('CLTIN2',CLTIN2 ,'')
    !ierr = write_nc('CVPOL',CVPOL  ,'')
    !ierr = write_nc('CIOPTJ',cioptj,'')
    !ierr = write_nc('CPCO',cpco  ,'')




    ierr = write_nc('XS',xs,['MAXNXS'],[maxnxs],'X bin boundaries','m')
    ierr = write_nc('YS',ys,['MAXNYS'],[maxnys],'Y bin boundaries','m')
    ierr = write_nc('PS',ps,['2MAXNPSP1'],[2*maxnps+1],'P bin boundaries','m')
    ierr = write_nc('XWIDS',xwids,['MAXNXS'],[maxnxs],'X bin widths','m')
    ierr = write_nc('YWIDS',ywids,['MAXNYS'],[maxnys],'Y bin widths','m')
    ierr = write_nc('PWIDS',pwids,['2MAXNPSP1'],[2*maxnps+1],'P bin widths','m')
    ierr = write_nc('PZONE',pzone,['2MAXNPSP1'],[2*maxnps+1],'Surface zone associated with P slice','m')
    ierr = write_nc('XOUTS',xouts,['MAXNXS'],[maxnxs],'X bin center','m')

    ierr = write_nc('YOUTS',youts,['2MAXNYSP1'],[2*maxnys+1],'Y bin center','m')
    ! Note: P bin centers are not currently calculated
    !ierr = write_nc('POUTS',pouts,['2MAXNPSP1'],[2*maxnps+1],'P bin center','m')

    !(DWELTS(IZ),IZ=0,NIZS)
    !(DWELFS(IT),IT=1,NTS)             

    ierr = write_nc('PIZS',pizs,['MAXNLS'],[maxnls],'Charge state for specific emission calculations')
    ierr = write_nc('PLAMS',plams,['MAXNLS'],[maxnls],'Wavelength for specific emission calculations')


    !
    !     Write out some 3D option information
    !
    !     CIOPTJ= 3D limiter extent option 
    !     CPCO  = 3D extent of limiter
    !
    !


    write(0,*) 'Writing to file: |***   |'
    ! subset of ddlims integrated over a smaller volume. 
    ierr = write_nc('SAVES',saves,['MAXNXS  ','MAXIZSP4'],[maxnxs,maxizs+4],'subset of ddlims integrated over a smaller volume','m-3 scaled')
    ierr = write_nc('DEPS',deps,['MAXNXS  ','MAXIZSP1','3       '],[maxnxs,maxizs+1,3],'Particle deposition')
    ierr = write_nc('NEROXS',neroxs,['MAXNXS','5     ','3     '],[maxnxs,5,3],'Net erosion along the X axis')
    ierr = write_nc('NEROYS',neroys,['MAXOS','6    '],[maxos,6],'Net erosion along Y?')
    ierr = write_nc('NERODS',nerods,['MAXOS','5    '],[maxos,5],'Net erosion along the surface')
    ierr = write_nc('NERODS3',nerods3,['MAXOS    ','2MAXNPSP1','6        '],[maxos,2*maxnps+1,6],'Net erosion along the 3D surface')
    ierr = write_nc('WALLS',walls,['2MAXNYSP1','MAXIZSP4 '],[2*maxnys+1,maxizs+4],'Deposition on walls')

    !tiz3(nxs,-ny3d:ny3d,-1:nizs,-maxnps:maxnps)
    ierr = write_nc('TIZ3',tiz3,    ['MAXNXS   ','2MAXY3DP1','MAXIZSP2 ','2MAXNPSP1'],[maxnxs,2*maxy3d+1,maxizs+2,2*maxnps+1],'Impurity ionization density results for 3D')
    
    !ddlim3(nxs,-ny3d:ny3d,-1:nizs,-maxnps:maxnps)
    ierr = write_nc('DDLIM3',ddlim3,['MAXNXS   ','2MAXY3DP1','MAXIZSP2 ','2MAXNPSP1'],[maxnxs,2*maxy3d+1,maxizs+2,2*maxnps+1],'Impurity density results for 3D')

    !plrps(nxs,-nys:nys,nls)
    ierr = write_nc('PLRPS',plrps,['MAXNXS   ','2MAXNYSP1','MAXNLS   '],[maxnxs,2*maxnys+1,maxnls],'Impurity particular line radiation profile emission')

    !plrp3(nxs,-ny3d:ny3d,nls,-maxnps:maxnps)
    ierr = write_nc('PLRP3',plrp3,['MAXNXS   ','2MAXY3DP1','MAXNLS   ','2MAXNPSP1'],[maxnxs,2*maxy3d+1,maxnls  ,2*maxnps+1],'Impurity particular line radiation profile emission 3D')
  
    write(0,*) 'Writing to file: |****  |'

    !oys(maxos)
    ierr = write_nc('OYS',oys,['MAXOS'],[maxos],'Y cell boundaries along surface','m')
    !ods(maxos)
    ierr = write_nc('ODS',ods,['MAXOS'],[maxos],'Distance cell boundaries along surface','m')

    !oyouts(maxos)
    ierr = write_nc('OYOUTS',oyouts,['MAXOS'],[maxos],'Y cell centers along surface','m')
    !odouts(maxos)
    ierr = write_nc('ODOUTS',odouts,['MAXOS'],[maxos],'Distance cell centers along surface','m')

    !oywids(maxos)
    ierr = write_nc('OYWIDS',oywids,['MAXOS'],[maxos],'Y cell widths along surface','m')
    !odwids(maxos)
    ierr = write_nc('ODWIDS',odwids,['MAXOS'],[maxos],'Distance cell widths along surface','m')

    !cdflux(maxos,3)
    ierr = write_nc('CDFLUX',cdflux,['MAXOS','3    '],[maxos,3],'Deuterium surface fluxes: parallel, cross-field and total')

    !ddlims(nxs,-nys:nys,-1:nizs)
    ierr = write_nc('DDLIMS',ddlims,['MAXNXS   ','2MAXNYSP1','MAXIZSP2 '],[maxnxs,2*maxnys+1,maxizs+2],'Impurity density results')

    !ddts(nxs,-nys:nys,1:nizs)
    ierr = write_nc('DDTS',ddts,['MAXNXS   ','2MAXNYSP1','MAXIZS   '],[maxnxs,2*maxnys+1,maxizs],'Impurity temperature','eV')

    !powls(nxs,-nys:nys,-1:nizs)
    ierr = write_nc('POWLS',powls,['MAXNXS   ','2MAXNYSP1','MAXIZSP2 '],[maxnxs,2*maxnys+1,maxizs+2],'Impurity radiated power')

    !lines(nxs,-nys:nys,-1:nizs)
    ierr = write_nc('LINES',lines,['MAXNXS   ','2MAXNYSP1','MAXIZSP2 '],[maxnxs,2*maxnys+1,maxizs+2],'Impurity line radiation')

    !tizs(nxs,-nys:nys,-1:nizs)
    ierr = write_nc('TIZS',tizs,['MAXNXS   ','2MAXNYSP1','MAXIZSP2 '],[maxnxs,2*maxnys+1,maxizs+2],'Impurity ionization density')

    !zeffs(nxs,-nys:nys,6)
    ierr = write_nc('ZEFFS',zeffs,['MAXNXS   ','2MAXNYSP1','6        '],[maxnxs,2*maxnys+1,6       ],'Impurity Z_effective')

    ! jdemod NC code doesn't have 5D arrays yet (nc_utils_generic.f90) ... since this array is not used much - leave it out for now
    ! August 21, 2018
    !lim5(nxs,-ny3d:ny3d,-1:nizs,-maxnps:maxnps,nts)
    !ierr = write_nc('LIM5',lim5,['MAXNXS   ','2MAXY3DP1','MAXIZSP2 ','2MAXNPSP1','MAXNTS   '],[maxnxs,2*maxy3d+1,maxizs+2,2*maxnps+1,maxnts],'Impurity density results for time dependent 3D')      

    !
    ! powl3 and line3 are not currently calculated and saved to the nc file. If these quantities are desired they can be added later or calculated
    ! as part of a post processor. The calculation of these can be found in the DMPOUT routine in iolim.f
    !
    !


    !sdtxs(nxs,1:nizs)
    ierr = write_nc('SDTXS',sdtxs,['MAXNXS','MAXIZS'],[maxnxs,maxizs],'Average impurity temperature along X','eV')


    !sdtys(-nys:nys,1:nizs)
    ierr = write_nc('SDTYS',sdtys,['2MAXNYSP1','MAXIZS   '],[2*maxnys+1,maxizs],'Average impurity temperature along Y','eV')



    !sdyxs(nxs,1:nizs)
    ierr = write_nc('SDYXS',sdyxs,['MAXNXS','MAXIZS'],[maxnxs,maxizs],'Average parallel diffusive step size along X')


    !sdyys(-nys:nys,1:nizs)
    ierr = write_nc('SDYYS',sdyys,['2MAXNYSP1','MAXIZS   '],[2*maxnys+1,maxizs],'Average parallel diffusive step size along Y')

    !sctxs(nxs,1:nizs)
    ierr = write_nc('SCTXS',sctxs,['MAXNXS','MAXIZS'],[maxnxs,maxizs],'Average spectroscopic weighted temperature along X','eV')
    !sctys(-nys:nys,1:nizs)
    ierr = write_nc('SCTYS',sctys,['2MAXNYSP1','MAXIZS   '],[2*maxnys+1,maxizs],'Average spectroscopic weighted temperature along Y','eV')


    !qxs(-nqxso:nqxsi)
    ierr = write_nc('QXS',qxs,['2MAXQXSP1'],[2*maxqxs+1],'X coordinates for detailed data')

    !qs(-nqxso:nqxsi)
    ierr = write_nc('QS',qxs,['2MAXQXSP1'],[2*maxqxs+1],'Timestep multiplier at X coordinates in QXS')

    !svybar(-nqxso:nqxsi)
    ierr = write_nc('SVYBAR',svybar,['2MAXQXSP1'],[2*maxqxs+1],'Average impurity velocity at X coordinates in QXS')


    !svyacc(-nqxso:nqxsi)
    ierr = write_nc('SVYACC',svyacc,['2MAXQXSP1'],[2*maxqxs+1],'Accumulated impurity weight at X coordinates in QXS')
    write(0,*) 'Writing to file: |***** |'

    !qedges(-nqxso:0,2)
    ierr = write_nc('QEDGES',qedges,['MAXQXSP1','2       '],[maxqxs+1,2],'Y coordinates of limiter edges along X - both sides')
    !qtans(-nqxso:0,2)
    ierr = write_nc('QTANS',qtans,['MAXQXSP1','2       '],[maxqxs+1,2],'Tangent angle of limiter edges along X - both sides')
    !qdists(-nqxso:0,2)
    ierr = write_nc('QDISTS',qdists,['MAXQXSP1','2       '],[maxqxs+1,2],'Total distance along limiter surface at each X - both sides')


    !qtembs(-nqxso:1,2)
    ierr = write_nc('QTEMBS',qtembs,['MAXQXSP2','2       '],[maxqxs+2,2],'Te along limiter surface along X - both sides')
    !qtembsi(-nqxso:1,2)
    ierr = write_nc('QTEMBSI',qtembsi,['MAXQXSP2','2       '],[maxqxs+2,2],'Ti along limiter surface along X - both sides')
    !qrnbs(-nqxso:1,2)
    ierr = write_nc('QRNBS',qrnbs,['MAXQXSP2','2       '],[maxqxs+2,2],'Plasma density along limiter surface along X - both sides')



    !facta(-1:nizs)

    ierr = write_nc('FACTA',facta,['MAXIZSP2'],[maxizs+2],'Charge state scaling factor')

    !factb(-1:nizs)
    ierr = write_nc('FACTB',factb,['MAXIZSP2'],[maxizs+2],'Charge state scaling factor with timestep')

    !TC
    ierr = write_nc('TC',tc,'Belt limiter geometry definition TC')
    !SC
    ierr = write_nc('SC',sc,'Belt limiter geometry definition SC')
    !TO
    ierr = write_nc('TO',to,'Belt limiter geometry definition TO')
    !SO
    ierr = write_nc('SO',so,'Belt limiter geometry definition SO')
    !TV
    ierr = write_nc('TV',tv,'Belt limiter geometry definition TV')
    !SV
    ierr = write_nc('SV',sv,'Belt limiter geometry definition SV')
    !GC
    ierr = write_nc('GC',gc,'Belt limiter geometry definition GC')
    !RP 
    ierr = write_nc('RP',rp,'Belt limiter geometry definition RP')


    !ctembs(nxs,-nys:nys)
    ierr = write_nc('CTEMBS',ctembs,['MAXNXS   ','2MAXNYSP1'],[maxnxs,2*maxnys+1],'Background Electron temperature','eV')

    !ctembsi(nxs,-nys:nys)
    ierr = write_nc('CTEMBSI',ctembsi,['MAXNXS   ','2MAXNYSP1'],[maxnxs,2*maxnys+1],'Background Ion temperature','eV')
    !crnbs(nxs,-nys:nys)
    ierr = write_nc('CRNBS',crnbs,['MAXNXS   ','2MAXNYSP1'],[maxnxs,2*maxnys+1],'Background Plasma density','m-3')



    ! write track data if available
    if (cstept.gt.0) then 
       !   cstept
       ierr = write_nc('CSTEPT',CSTEPT ,'Number of particle tracks recorded')
       !   ptracl(maxt)
       ierr = write_nc('PTRACL',ptracl,['MAXT'],[maxt],'Length of each recorded particle track')
       !   ptracs(maxlen,maxt,2)
       ierr = write_nc('PTRACS',ptracs,['MAXLEN','MAXT  ','2     '],[maxlen,maxt,2],'X,Y coordinates for each particle track')
	endif

	! sazmod - Write varying absorbing wall data if available.
	if (vary_absorb.eq.1) then
	  ierr = write_nc('yabsorb1a_step', yabsorb1a_step, 'Y location of left step in absorbing boundary', 'm')
	  ierr = write_nc('yabsorb2a_step', yabsorb2a_step, 'Y location of right step in absorbing boundary', 'm')
	  ierr = write_nc('xabsorb1a_step', xabsorb1a_step, 'X location of left step in absorbing boundary', 'm')
      ierr = write_nc('xabsorb2a_step', xabsorb2a_step, 'X location of right step in absorbing boundary', 'm')

    endif
    write(0,*) 'Writing to file: |******|'


    !
    ! jdemod - close the netcdf file
    !
    ierr=close_nc_file()
    if (ierr.ne.0) then 
       call errmsg('WRITE_NETCDF_OUTPUT: ERROR CLOSING OUTPUT FILE:',ierr)
       return
    endif

    return
  end subroutine write_netcdf_output




end module lim_netcdf
