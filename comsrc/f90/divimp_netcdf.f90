module divimp_netcdf
  use error_handling
  use debug_options
  use nc_utils_generic


contains

  !
  !     Netcdf output
  !
  subroutine write_netcdf_output(TITLE,desc,NIZS,JOB,EQUIL,ITER,NITERS)
  !subroutine write_netcdf_output(TITLE,desc,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)
    use mod_params

    !
    ! DIVIMP common blocks used as modules
    !

    use subgrid
    use mod_divimp
    use mod_comtor
    use mod_cneut
    use mod_commv
    use mod_cadas
    use mod_cadas2
    use mod_cneut2
    use mod_cgeom
    use mod_dynam1
    use mod_dynam3
    use mod_dynam4
    use mod_printopt
    use mod_pindata
    use mod_divxy
    use mod_grbound
    use mod_cedge2d
    use mod_transcoef
    use mod_cioniz
    use mod_promptdep
    use mod_reiser_com
    use mod_line_profile
    use mod_hc_global_opts
    use mod_driftvel
    use mod_diagvel
    use mod_slcom
    use mod_fp_data

    implicit none


    CHARACTER TITLE*(*),desc*(*),JOB*(*),EQUIL*(*)
    INTEGER   NIZS,ITER,NITERS
    !REAL      FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)
    real :: tmp_pinline(maxnks,maxnrs,6)
    integer :: ierr
    INTEGER IR,IZ,IT,IN
    INTEGER      i1,i2,i3,ik,i
    REAL         slver

    !
    !real :: test_array(5,5,5)
    !test_array = 1.0

    !
    ! Fortran note: static string array assignments require all elements to have the same length
    !

    call pr_trace('WRITE_NETCDF_OUTPUT','START')

    !
    ! jdemod - open the netcdf file for output
    !
    ierr=open_nc_file('divimp_netcdf_out.nc',NC_WRITE)
    if (ierr.ne.0) then 
       call errmsg('WRITE_NETCDF_OUTPUT: ERROR OPENING OUTPUT FILE:',ierr)
       return
    endif

    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER OPEN')

    !c slmod end
    !C
    !c
    !c     Version number 
    !c     
    !c     
    !c     Basic information - title, job description,
    !c                         equilibrium, shot number and time slice
    !c     
    !      write(8) verson
    ierr = write_nc('VERSION',verson,'Code version')


    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER VERSION')

    !c
    !c     Write out the global parameters used to write the file
    !c     - this is the first step in parameterizing the read statements in 
    !c       OUT and removing the dependency on having identical parameter 
    !c       values in both the DIVIMP and OUT compiles
    !c
    !      WRITE  (8) TITLE,desc,JOB,EQUIL,ISHOT,TSLICE
    
    !write(0,*) 'DESC:',len_trim(desc),':',trim(desc),':'

    ierr = write_nc('TITLE',title,'Case Title')
    ierr = write_nc('DESC',desc,'Case Description')
    ierr = write_nc('JOB',job,'Job Description')
    ierr = write_nc('EQUIL',equil,'Equil Description')
    ierr = write_nc('ISHOT',equil,'May be shot')
    ierr = write_nc('TSLICE',equil,'May be time slice')

    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER TITLE BLOCK')

    !c
    !c     These are fundamental array dimensions ... will be stored 
    !c     separately as dimension data but can also be stored as 
    !c     netcdf variables - these are compile time parameters in DIVIMP
    !c     at the moment
    !c
    !
    !write(8) MAXNKS,MAXNRS,MAXNDS,MAXNGS,MAXIZS,MAXINS,                
    !>  MAXIMP,ISECT,MAXNTS,MAXNOC,MAXNWS,MAXNXS,MAXNYS,MAXVMF,         
    !>  MAXTHE,MAXSN,MAXPTS,MAXPLRP,MSOLPT,MAXADS,      
    !>  MAXGXS,MAXGYS,MAXCH3,MAXIXS,MAXIYS,MAXSEG,maxplts,MAXNFLA,
    !>  maxpiniter,mbufle,mbufx,mves,maxe2dizs

    !      write(8) MAXNKS,MAXNRS,MAXNDS,MAXNGS,MAXIZS,MAXINS,MAXIMP,ISECT,MAXNTS,MAXNOC,MAXNWS,MAXNXS,MAXNYS,MAXVMF,MAXTHE,MAXSN,MAXPTS,MAXPLRP,MSOLPT,MAXADS,MAXGXS,MAXGYS,MAXCH3,MAXIXS,MAXIYS,MAXSEG,maxplts,MAXNFLA,maxpiniter,mbufle,mbufx,mves,maxe2dizs

    ierr = write_nc('MAXNKS',maxnks,'Maximum Knots on a ring')
    ierr = write_nc('MAXNRS',maxnrs,'Maximum Rings on a grid')
    ierr = write_nc('MAXNDS',maxnds,'Maximum Number of target elements')
    ierr = write_nc('MAXNGS',maxngs,'Maximum Number of plots')
    ierr = write_nc('MAXIZS',maxizs,'Maximum Number of impurity charge states')
    ierr = write_nc('MAXINS',maxngs,'Maximum Number of input values for some arrays')
    ierr = write_nc('MAXIMP',maximp,'Maximum Number of impurity particles to follow')

    ierr = write_nc('ISECT',isect,'Chunking step for random number arrays')
    ierr = write_nc('MAXNTS',maxnts,'Maximum number of time steps')
    ! not used - jdefix
    !ierr = write_nc('MAXNOC',maxnoc,'Does nothing apparently ')
    ierr = write_nc('MAXNWS',maxnws,'Maximum number of particle steps recorded for debugging')
    ierr = write_nc('MAXNXS',maxnxs,'Maximum xbins in XY grid - deprecated')
    ierr = write_nc('MAXNYS',maxnys,'Maximum ybins in XY grid - deprecated')
    ierr = write_nc('MAXVMF',maxvmf,'Maximum number of velocity multiplier factors')

    ierr = write_nc('MAXTHE',maxthe,'Maximum number of points in line of sight plots')

    ! jdemod - maxsn - not used
    !ierr = write_nc('MAXSN',maxsn,'Maximum Rings on a grid')

    ierr = write_nc('MAXPTS',maxpts,'Maximum Number of target elements')
    ierr = write_nc('MAXPLRP',maxplrp,'Maximum Number of particular line radiation profiles (OLD)')

    ierr = write_nc('MSOLPT',msolpt,'Number of points in high resolution SOL profiles')
    ierr = write_nc('MAXADS',maxads,'Array sizes for ADAS calls')

    ierr = write_nc('MAXGXS',maxgxs,'Maximum X resolution of plotting grid (OLD)')
    ierr = write_nc('MAXGYS',maxgys,'Maximum Y resolution of plotting grid (OLD)')
    ierr = write_nc('MAXCH3',maxch3,'Maximum Rings on a grid')
    ierr = write_nc('MAXIXS',maxixs,'Maximum X resolution of OUT XY grid (OLD)')
    ierr = write_nc('MAXIYS',maxiys,'Maximum Y resolution of OUT XY grid (OLD)')
    ierr = write_nc('MAXSEG',maxseg,'Maximum Number of wall segments for H flux arrays')
    ierr = write_nc('MAXPLTS',maxplts,'Maximum Number of plots in along ring plots in OUT')

    ierr = write_nc('MAXNFLA',maxnfla,'Maximum Number of fluids in fluid code solution')
    ierr = write_nc('MAXPINITER',maxpiniter,'Maximum Number of PIN iterations')
    ierr = write_nc('MBUFLE',mbufle,'Maximum size for baffle structure arrays (OLD)')
    ierr = write_nc('MBUFX',mbufx,'Maximum size of another baffle related array (OLD)')
    ierr = write_nc('MVES',mves,'Maximum Number of vessel wall vertices')
    ierr = write_nc('MAXE2DIZS',maxe2dizs,'Maximum Number of fluid code charge states')

    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 1')

    !c     
    !c     Simulation values
    !c     

    !c     
    !c     Geometry values 
    !c     
    !      write(8) ITER,NITERS,NIZS,NTS,CRMB,CRMI,CIZB,CION,IMODE,NITERSOL

    !c     
    !c     Simulation values
    !c     
    !      
    !      write(8) ITER,NITERS,NIZS,NTS,CRMB,CRMI,CIZB,CION,IMODE,NITERSOL
    !      write(8) ITER,NITERS,NIZS,NTS,CRMB,CRMI,CIZB,CION,IMODE,
    !     >         NITERSOL

    ierr = write_nc('ITER',iter,'Iteration count')
    ierr = write_nc('NITERS',niters,'Number of DIVIMP iterations')
    ierr = write_nc('NIZS',nizs,'Number of impurity charge states followed')
    ierr = write_nc('NTS',nts,'Number of time windows in time dependent simulations')
    ierr = write_nc('CRMB',crmb,'Mass of plasma species in amu','amu')
    ierr = write_nc('CRMI',crmi,'Mass of impurity species in amu','amu')
    ierr = write_nc('CION',cion,'Atomic number of impurity species')
    ierr = write_nc('IMODE',imode,'Operation mode - 1=time dependent, 2=steady state')
    ierr = write_nc('NITERSOL',nitersol,'Number of plasma solver iterations')


    !c     
    !c     Geometry values 
    !c     
    !
    !
    !      write(8) R0,Z0,RXP,ZXP,NRS,
    !     >         RMIN,RMAX,ZMIN,ZMAX,DR,DZ,NDS,NDSIN,NXS,NYS,
    !     >         IRCENT,IRSEP,IRWALL,IRTRAP,IKT,IKREF,IKTO,IKTI,
    !     >         irsep2,irwall2,nrs2,ndsin2,ndsin3,irtrap2,
    !     >         nves,nvesm,nvesp,inmid,oumid,refct,CIRHR,NPOLYP,
    !     >         cxnear,cynear
    !c
    !      write(8) R0,Z0,RXP,ZXP,NRS,RMIN,RMAX,ZMIN,ZMAX,DR,DZ,NDS,NDSIN,NXS,NYS,IRCENT,IRSEP,IRWALL,IRTRAP,IKT,IKREF,IKTO,IKTI,irsep2,irwall2,nrs2,ndsin2,ndsin3,irtrap2,nves,nvesm,nvesp,inmid,oumid,refct,CIRHR,NPOLYP,cxnear,cynear

    ierr = write_nc('R0',r0,'Plasma center R coordinate','m')
    ierr = write_nc('Z0',z0,'Plasma center Z coordinate','m')
    ierr = write_nc('RXP',rxp,'Xpoint R coordinate','m')
    ierr = write_nc('ZXP',zxp,'Xpoint Z coordinate','m')
    ierr = write_nc('NRS',nrs,'Number of rings on the grid')
    ierr = write_nc('RMIN',rmin,'Minimum R for grid','m')
    ierr = write_nc('RMAX',rmax,'Maximum R for grid','m')
    ierr = write_nc('ZMIN',zmin,'Minimum Z for grid','m')
    ierr = write_nc('ZMAX',zmax,'Maximum Z for grid','m')

    !ierr = write_nc('DR',dr,'Delta R bin size (not used)')
    !ierr = write_nc('DZ',dz,'Delta Z bin size (not used)')

    ierr = write_nc('NDS',nds,'Number of elements of target surfaces')
    ierr = write_nc('NDSIN',ndsin,'First index of inner target')

    !ierr = write_nc('NXS',nxs,'Number of X bins in XY grid (not used)')
    !ierr = write_nc('NYS',nys,'Number of Y bins in XY grid (not used)')

    ierr = write_nc('IRCENT',ircent,'Innermost ring of core')
    ierr = write_nc('IRSEP',irsep,'Index of first ring in main SOL')
    ierr = write_nc('IRWALL',irwall,'Index of outermost ring in main SOL')
    ierr = write_nc('IRTRAP',irtrap,'Index of outermost ring in PFZ')

    ierr = write_nc('IKT',ikt,'Knot index offset for mapping')
    ierr = write_nc('IKREF',ikref,'Midpoint cell on separatrix by cell index')
    ierr = write_nc('IKTO',ikto,'Outer leg knot index offset')
    ierr = write_nc('IKTI',ikti,'Inner leg knot index offset')
    ierr = write_nc('IRSEP2',irsep2,'Second separatrix in double null grid')
    ierr = write_nc('IRWALL2',irwall2,'Second wall ring in double null grid')
    ierr = write_nc('NRS2',nrs2,'Total number of rings in double null grid')
    ierr = write_nc('NDSIN2',ndsin2,'First index of second target in double null')
    ierr = write_nc('NDSIN3',ndsin3,'First index of third target in double null')
    ierr = write_nc('IRTRAP2',irtrap2,'Index of outermost ring in second PFZ')

    ierr = write_nc('NVES',nves,'Number of vessel elements')
    ierr = write_nc('NVESM',nvesm,'Number of main wall elements')
    ierr = write_nc('NVESP',nvesp,'Number of pfz wall elements')
    ierr = write_nc('INMID',inmid,'Inner midplane IK index')
    ierr = write_nc('OUMID',oumid,'Outer midplane IK index')

    ierr = write_nc('REFCT',refct,'Mass of impurity species in amu','amu')

    !ierr = write_nc('CIRHR',cirhr,'Ring number of high resolution background profiles')

    ierr = write_nc('NPOLYP',npolyp,'Number of grid polygons')

    !ierr = write_nc('CXNEAR',cxnear,'Number of plasma solver iterations')
    !ierr = write_nc('CYNEAR',cynear,'Number of plasma solver iterations')



    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 2')

    !c     
    !c     Scaling factors  
    !c     
    !      CALL IINOUT ('W NKS    ',nks ,maxnrs)
    ierr = write_nc('NKS',nks,['MAXNRS'],[maxnrs],'Number of knots on each ring')
    !c     
    !      write(8) QTIM,FSRATE,ABSFAC,absfac_neut,CBPHI,CK0
    ierr = write_nc('QTIM',qtim,'Time step for ions','s')
    ierr = write_nc('FSRATE',fsrate,'Time step for neutrals','s')
    ierr = write_nc('ABSFAC',absfac,'Absolute scaling factor (part/s/m-tor)','part/s/m-tor')
    ierr = write_nc('ABSFAC_NEUT',absfac_neut,'Absolute scaling factor based on neutrals (part/s/m-tor)','part/s/m-tor')
    ierr = write_nc('CBPHI',cbphi,'Toroidal On-axis B-field','T')
    ierr = write_nc('CK0',ck0,'Electron heat conduction coefficient')
    ierr = write_nc('CK0I',ck0i,'Ion heat conduction coefficient')
    !      CALL RINOUT ('W FACTA  ',facta ,maxizs+2)
    ierr = write_nc('FACTA',facta,['MAXIZSP2'],[maxizs+2],'Ion scaling factor')
    !      CALL RINOUT ('W FACTB  ',factb ,maxizs+2)
    ierr = write_nc('FACTB',factb,['MAXIZSP2'],[maxizs+2],'Electron scaling factor')
    !c slmod begin
    !      CALL RINOUT ('W DWELTS ',dwelts,maxizs+2)
    ierr = write_nc('DWELTS',dwelts,['MAXIZSP2'],[maxizs+2],'Dwell times')

    IF (IMODE.EQ.1) THEN
       !        CALL RINOUT ('W DWELFS ',dwelfs,maxnts)
       ierr = write_nc('DWELFS',dwelfs,['MAXNTS'],[maxnts],'Dwell time factors')
       !ELSE
       !        CALL RINOUT ('W DWELFS ',dwelfs,1     )
       !ierr = write_nc('DWELFS',dwelfs,['1'],[1],'')
       !c
       !c      CALL RINOUT ('W DWELFS ',dwelfs,maxnts)
       !c slmod end
    ENDIF
    !      CALL RINOUT ('W KALPHS ',kalphs,maxizs)
    ierr = write_nc('KALPHS',kalphs,['MAXIZS'],[maxizs],'Electron temperature gradient transport coefficients')
    !c     
    !c     DIVIMP Options
    !c     
    !      CALL RINOUT ('W KBETAS ',kbetas,maxizs)
    ierr = write_nc('KBETAS',kbetas,['MAXIZS'],[maxizs],'Ion temperature gradient transport coefficients')
    !c     
    !c     Result numbers
    !c     
    !      write(8) CIOPTF,cdatopt,cneur,cgridopt,cre2d,cre2dizs,xygrid
    ierr = write_nc('CIOPTF',cioptf,'SOL option specified in input file')
    ierr = write_nc('CDATOPT',cdatopt,'Atomic data source option')
    ierr = write_nc('CNEUR',cneur,'Neutral wall location option')
    ierr = write_nc('CGRIDOPT',cgridopt,'Grid option')
    ierr = write_nc('CRE2D',cre2d,'Fluid code plasma solution read option')
    ierr = write_nc('CRE2DIZS',cre2dizs,'Fluid code impurity read option')
    ierr = write_nc('XYGRID',xygrid,'XY grid option - not used')
    !c     
    !c     ADAS       
    !c     
    !      write(8) nleakcore,cleaksn
    ierr = write_nc('NLEAKCORE',nleakcore,'Particles crossing separatrix')
    ierr = write_nc('CLEAKSN',cleaksn,'Number of leakage boundaries in SOL')
    !c     
    !c     
    !c     The following are the piece-wise wall specifications
    !c     from the GA15 routines.
    !c     
    !      write(8) useridh,iyearh,useridz,iyearz
    ierr = write_nc('USERIDH',useridh,'ADAS H user id')
    ierr = write_nc('IYEARH',iyearh,'ADAS H year')
    ierr = write_nc('USERIDZ',useridz,'ADAS Impurity User ID')
    ierr = write_nc('IYEARZ',iyearz,'ADAS impurity year')


    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 3')

    !c     
    !c     Wall definitions
    !c     
    !      WRITE(8) ioncpts,ionwpts
    ierr = write_nc('IONCPTS',ioncpts,'Number of points in core ion boundary')
    ierr = write_nc('IONWPTS',ionwpts,'Number of points in core SOL boundary')
    !      CALL RINOUT ('W RIW    ',riw ,maxpts)
    ierr = write_nc('RIW',riw,['MAXPTS'],[maxpts],'R coordinates of Ion wall','m')
    !      CALL RINOUT ('W ZIW    ',ziw ,maxpts)
    ierr = write_nc('ZIW',ziw,['MAXPTS'],[maxpts],'Z coordinates of Ion wall','m')
    !      CALL RINOUT ('W RCW    ',rcw ,maxpts)
    ierr = write_nc('RCW',rcw,['MAXPTS'],[maxpts],'R coordinates of ion core boundary','m')
    !      CALL RINOUT ('W ZCW    ',zcw ,maxpts)
    ierr = write_nc('ZCW',zcw,['MAXPTS'],[maxpts],'Z coordinates of ion core boundary','m')
    !      CALL RINOUT ('W RW     ',rw  ,maxpts)
    ierr = write_nc('RW',rw,['MAXPTS'],[maxpts],'R coordinates of neutral wall','m')
    !      CALL RINOUT ('W ZW     ',zw  ,maxpts)
    ierr = write_nc('ZW',zw,['MAXPTS'],[maxpts],'Z coordinates of neutral wall','m')
    !      CALL RINOUT ('W RVES   ',rves,maxpts)
    ierr = write_nc('RVES',rves,['MAXPTS'],[maxpts],'R vessel coordinates','m')
    !c     
    !c     GA15 workspace 
    !c     
    !      CALL RINOUT ('W ZVES   ',zves,maxpts)
    ierr = write_nc('ZVES',zves,['MAXPTS'],[maxpts],'Z vessel coordinates','m')


    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 4')
    !
    ! These are working storage for GA15 ... should not be in raw file in the first place
    ! I would think. 
    !c
    !c               The following are the piece-wise wall specifications
    !c               from the GA15 routines.
    !c
    !c     >          ioncpts,ionwpts,(riw(in),ziw(in),rcw(in),zcw(in),
    !c     >          rw(in),zw(in),rves(in),zves(in),in=1,maxpts),
    !c     >          ((iwindw(it,in),icindw(it,in),it=1,2),in=1,maxpts),
    !c     >          (iwwork(in),icwork(in),in=1,4*maxpts),
    !c     >          (iwtdum(in),iwxdum(in),iwydum(in),in=1,maxpts),
    !c     >          (ictdum(in),icxdum(in),icydum(in),in=1,maxpts),
    !
    !      CALL IINOUT ('W IWINDW ',iwindw,2*maxpts)
    !ierr = write_nc('IWINDW',iwindw,['2','MAXPTS'],[2,maxpts],'')
    !      CALL RINOUT ('W IWWORK ',iwwork,4*maxpts)
    !ierr = write_nc('IWWORK',iwwork,['4','MAXPTS'],[4,maxpts],'')
    !      CALL RINOUT ('W IWTDUM ',iwtdum,maxpts)
    !ierr = write_nc('IWTDUM',iwtdum,['MAXPTS'],[maxpts],'')
    !      CALL RINOUT ('W IWXDUM ',iwxdum,maxpts)
    !ierr = write_nc('IWXDUM',iwxdum,['MAXPTS'],[maxpts],'')
    !c     
    !      CALL RINOUT ('W IWYDUM ',iwydum,maxpts)
    !ierr = write_nc('IWYDUM',iwydum,['MAXPTS'],[maxpts],'')
    !      CALL IINOUT ('W ICINDW ',icindw,2*maxpts)
    !ierr = write_nc('ICINDW',icindw,['2','MAXPTS'],[2,maxpts],'')
    !      CALL RINOUT ('W ICWORK ',icwork,4*maxpts)
    !ierr = write_nc('ICWORK',icwork,['4','MAXPTS'],[4,maxpts],'')
    !      CALL RINOUT ('W ICTDUM ',ictdum,maxpts)
    !ierr = write_nc('ICTDUM',ictdum,['MAXPTS'],[maxpts],'')
    !      CALL RINOUT ('W ICXDUM ',icxdum,maxpts)
    !ierr = write_nc('ICXDUM',icxdum,['MAXPTS'],[maxpts],'')
    !      CALL RINOUT ('W ICYDUM ',icydum,maxpts)
    !ierr = write_nc('ICYDUM',icydum,['MAXPTS'],[maxpts],'')
    !c     
    !c     
    !c
    !c
    !c               The following relate to the wall definition as
    !c               calculated in DIVIMP.
    !c
    !c     >          wlwall1,wlwall2,wltrap1,wltrap2,wallpts
    !c
    !C
    !c
    !      WRITE(8) wlwall1,wlwall2,wltrap1,wltrap2,wallpts
    ierr = write_nc('WLWALL1',wlwall1,'First main wall element index')
    ierr = write_nc('WLWALL2',wlwall2,'Last main wall element index')
    ierr = write_nc('WLTRAP1',wltrap1,'First pfz wall element index')
    ierr = write_nc('WLTRAP2',wltrap2,'Last pfz wall element index')
    ierr = write_nc('WALLPTS',wallpts,'Number of elements in wall definition')
    !c     >  TITLE,JOB,ITER,IMODE,refct
    !c
    !WRITE (6,9001) NXS,NYS,NRS,NDS,NIZS,NTS,MAXNXS,MAXNYS,MAXNRS,MAXNDS,MAXIZS,MAXNTS,TITLE,JOB,ITER,IMODE,refct,maxseg,nvesm,nvesp
    !
    !      CALL RINOUT ('W POWLS ',POWLS ,MAXNKS*MAXNRS*(MAXIZS+2))
    ierr = write_nc('POWLS',powls,['MAXNKS  ','MAXNRS  ','MAXIZSP2'],[maxnks,maxnrs,maxizs+2],'Total impurity radiated power','W/m3/particle')
    !      CALL RINOUT ('W LINES ',LINES ,MAXNKS*MAXNRS*(MAXIZS+2))
    ierr = write_nc('LINES',lines,['MAXNKS  ','MAXNRS  ','MAXIZSP2'],[maxnks,maxnrs,maxizs+2],'Total impurity line radiation','W/m3/particle')
    !      CALL RINOUT ('W HPOWLS',HPOWLS,MAXNKS*MAXNRS*2)
    ierr = write_nc('HPOWLS',hpowls,['MAXNKS','MAXNRS','2     '],[maxnks,maxnrs,2],'Total hydrogen radiated power','W/m3')
    !      CALL RINOUT ('W HLINES',HLINES,MAXNKS*MAXNRS*2)
    ierr = write_nc('HLINES',hlines,['MAXNKS','MAXNRS','2     '],[maxnks,maxnrs,2],'Total hydrogen line radiation','W/m3')
    !      CALL RINOUT ('W TIZS  ',TIZS  ,MAXNKS*MAXNRS*(MAXIZS+2))
    ierr = write_nc('TIZS',tizs,['MAXNKS  ','MAXNRS  ','MAXIZSP2'],[maxnks,maxnrs,maxizs+2],'Impurity ionization')
    !      CALL RINOUT ('W ZEFFS ',ZEFFS ,MAXNKS*MAXNRS*3)
    ierr = write_nc('ZEFFS',zeffs,['MAXNKS','MAXNRS','3     '],[maxnks,maxnrs,3],'Code calculated Zeffective')
    !      CALL RINOUT ('W WALLS ',WALLS ,MAXNKS*MAXNRS*(MAXIZS+2))
    ierr = write_nc('WALLS',walls,['MAXNKS  ','MAXNRS  ','MAXIZSP2'],[maxnks,maxnrs,maxizs+2],'Impurity wall deposition (grid exits)')
    !      CALL RINOUT ('W DEPS  ',DEPS  ,MAXNDS*MAXIZS)
    ierr = write_nc('DEPS',deps,['MAXNDS','MAXIZS'],[maxnds,maxizs],'Impurity deposition on targets')
    !      CALL RINOUT ('W NEROS ',NEROS ,MAXNDS*5)
    ierr = write_nc('NEROS',neros,['MAXNDS','5     '],[maxnds,5],'Impurity erosion/depositon/net erosion from targets')
    !      CALL RINOUT ('W PRDEPS',PROMPTDEPS,MAXNDS*9)
    ierr = write_nc('PRDEPS',promptdeps,['MAXNDS','9     '],[maxnds,9],'Impurity prompt deposition statistics')
    !      CALL RINOUT ('W WALLSN',WALLSN,MAXPTS+1)
    ierr = write_nc('WALLSN',wallsn,['MAXPTSP1'],[maxpts+1],'Neutral deposition on walls')
    !      CALL RINOUT ('W WALLSE',WALLSE,MAXPTS+1)
    ierr = write_nc('WALLSE',wallse,['MAXPTSP1'],[maxpts+1],'Total erosion (source) by wall index')
    !      CALL RINOUT ('W WALLSEI',WALLSE_I,MAXPTS+1)
    ierr = write_nc('WALLSEI',wallse_i,['MAXPTSP1'],[maxpts+1],'Total eroded particles that are ionized')
    !      CALL RINOUT ('W WALLSI',WALLSI,MAXPTS+1)
    ierr = write_nc('WALLSI',wallsi,['MAXPTSP1'],[maxpts+1],'Ion deposition on walls')
    !c
    !      CALL RINOUT ('W WALLSIL',WALLSIL,MAXPTS+1)
    ierr = write_nc('WALLSIL',wallsil,['MAXPTSP1'],[maxpts+1],'Ionized particles that leak out of the divertor')
    !      CALL RINOUT ('W WALLPT',WALLPT,MAXPTS*32)
    ierr = write_nc('WALLPT',wallpt,['MAXPTS','32    '],[maxpts,32],'Detailed wall element information')
    !      CALL RINOUT ('W RS    ',RS    ,MAXNKS*MAXNRS)
    ierr = write_nc('RS',rs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'R coordiante of cell centers','m')
    !      CALL RINOUT ('W ZS    ',ZS    ,MAXNKS*MAXNRS)
    ierr = write_nc('ZS',zs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Z coordinate of cell centers','m')
    !      CALL RINOUT ('W KSB   ',KSB   ,(MAXNKS+1)*MAXNRS)
    ierr = write_nc('KSB',ksb,['MAXNKSP1','MAXNRS  '],[maxnks+1,maxnrs],'S (parallel) coordinate boundaries for cells (0 to nks(ir)','m')
    !      CALL RINOUT ('W KPB   ',KPB   ,(MAXNKS+1)*MAXNRS)
    ierr = write_nc('KPB',kpb,['MAXNKSP1','MAXNRS  '],[maxnks+1,maxnrs],'Poloidal boundaries of each cell','m')



    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 5')

    !c
    !c     The storing of these arrays needed to be customized
    !c     because of a likely size mismatch between the
    !c     array in DIVIMP and that in OUT.
    !c
    !c      CALL IINOUT ('W IKXYS ',IKXYS ,MAXNXS*MAXNYS)
    !c      CALL IINOUT ('W IRXYS ',IRXYS ,MAXNXS*MAXNYS)
    !c      CALL IINOUT ('W IFXYS ',IFXYS ,MAXNXS*MAXNYS)
    !c
    !
    ! XY grid data not needed 
    !
    !CALL IINOUT2 ('W IKXYS ',IKXYS ,MAXNXS,MAXNYS,MAXNXS,MAXNYS)
    !CALL IINOUT2 ('W IRXYS ',IRXYS ,MAXNXS,MAXNYS,MAXNXS,MAXNYS)
    !c
    !CALL IINOUT2 ('W IFXYS ',IFXYS ,MAXNXS,MAXNYS,MAXNXS,MAXNYS)
    !      CALL IINOUT ('W IKDS  ',IKDS  ,MAXNDS)
    !
    ierr = write_nc('IKDS',ikds,['MAXNDS'],[maxnds],'IK grid index of target index ID')
    !C
    !      CALL IINOUT ('W IRDS  ',IRDS  ,MAXNDS)
    ierr = write_nc('IRDS',irds,['MAXNDS'],[maxnds],'IR grid index of target index ID')
    !      CALL IINOUT ('W IKINS ',IKINS ,MAXNKS*MAXNRS)
    ierr = write_nc('IKINS',ikins,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Connection map - IK index of next inward cell')
    !      CALL IINOUT ('W IKOUTS',IKOUTS,MAXNKS*MAXNRS)
    ierr = write_nc('IKOUTS',ikouts,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Connection map - IK index of next outward cell')
    !      CALL IINOUT ('W IRINS ',IRINS ,MAXNKS*MAXNRS)
    ierr = write_nc('IRINS',irins,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Connection map - IR index of next inward cell')
    !      CALL IINOUT ('W IROUTS',IROUTS,MAXNKS*MAXNRS)
    ierr = write_nc('IROUTS',irouts,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Connection map - IR index of next outward cell')


    !      CALL IINOUT ('W KORY  ',KORY  ,MAXNKS*MAXNRS)
    !ierr = write_nc('KORY',kory,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'NIMBUS neutral code index mapping (not used)')
    !      CALL IINOUT ('W KORPG ',KORPG ,MAXNKS*MAXNRS)

    ierr = write_nc('KORPG',korpg,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'IK,IR mapping to polygon index')
    !      CALL IINOUT ('W NVERTP',NVERTP,MAXNKS*MAXNRS)

    ierr = write_nc('NVERTP',nvertp,['MAXNKS*MAXNRS'],[maxnks*maxnrs],'Number of vertices of grid polygons')
    !      CALL RINOUT ('W RVERTP',RVERTP,5*MAXNKS*MAXNRS)
    ierr = write_nc('RVERTP',rvertp,['5            ','MAXNKS*MAXNRS'],[5,maxnks*maxnrs],'R corner coordinates of grid polygons')
    !C
    !      CALL RINOUT ('W ZVERTP',ZVERTP,5*MAXNKS*MAXNRS)
    ierr = write_nc('ZVERTP',zvertp,['5            ','MAXNKS*MAXNRS'],[5,maxnks*maxnrs],'Z corner coordinates of grid polygons')

    !c
    !CALL CHECK_DDLIM(nizs,3)
    !      CALL DINOUT ('W DDLIMS',DDLIMS,MAXNKS*MAXNRS*(MAXIZS+2))
    ierr = write_nc('DDLIMS',ddlims,['MAXNKS  ','MAXNRS  ','MAXIZSP2'],[maxnks,maxnrs,maxizs+2],'Impurity density by charge state','/m3/part')
    !      CALL DINOUT ('W DDTS  ',DDTS  ,MAXNKS*MAXNRS*(MAXIZS+2))
    ierr = write_nc('DDTS',ddts,['MAXNKS  ','MAXNRS  ','MAXIZSP2'],[maxnks,maxnrs,maxizs+2],'Impurity temperurature by charge state','eV')
    !      CALL RINOUT ('W ELIMS ',ELIMS ,MAXNKS*3*(MAXIZS+2))
    ierr = write_nc('ELIMS',elims,['MAXNKS  ','3       ','MAXIZSP2'],[maxnks,3,maxizs+2],'Locations of particles crossing into core')
    !c
    !      CALL RINOUT ('W WALKS ',WALKS ,MAXNWS*2)
    ierr = write_nc('WALKS',walks,['MAXNWS','2     '],[maxnws,2],'Particle tracks')
    !      call rinout ('W CHEM D',chemden,maxnks*maxnrs)
    ierr = write_nc('CHEM D',chemden,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Density of chemically sputtered species')
    !C
    !      call rinout ('W CHEMIZ',chemizs,maxnks*maxnrs)
    ierr = write_nc('CHEMIZ',chemizs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Ionization of chemically sputtered species')
    !      CALL RINOUT ('W KKS   ',KKS   ,MAXNRS)

    ierr = write_nc('KKS',kks,['MAXNRS'],[maxnrs],'Arbitrary "K" value for each ring')
    !      CALL RINOUT ('W KSS   ',KSS   ,MAXNKS*MAXNRS)
    ierr = write_nc('KSS',kss,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'S coordinate of cell centers along the field lines','m')
    !      CALL RINOUT ('W KPS   ',KPS   ,MAXNKS*MAXNRS)
    ierr = write_nc('KPS',kps,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Poloidal coordinate of cell centers along the field lines','m')
    !
    ! knorms - not used -jdefix
    !
    !      CALL RINOUT ('W KNORMS',KNORMS,MAXNKS*MAXNRS)
    !ierr = write_nc('KNORMS',knorms,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')
    !      CALL RINOUT ('W KPERPS',KPERPS,MAXNKS*MAXNRS)

    ierr = write_nc('KPERPS',kperps,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Cross field step size by cell','m')

    !  kcurvs - not used - jdefix
    !
    !      CALL RINOUT ('W KCURVS',KCURVS,MAXNKS*MAXNRS)
    !ierr = write_nc('KCURVS',kcurvs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')


    !      CALL RINOUT ('W KVOLS ',KVOLS ,MAXNKS*MAXNRS)
    ierr = write_nc('KVOLS',kvols,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Cell volume (2PI*R*KAREA)','m3')
    !      CALL RINOUT ('W KAREAS',KAREAS,MAXNKS*MAXNRS)
    ierr = write_nc('KAREAS',kareas,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Cell area','m2')
    !      CALL RINOUT ('W KTOTAS',KTOTAS,MAXNRS)
    ierr = write_nc('KTOTAS',ktotas,['MAXNRS'],[maxnrs],'Total area of each flux tube','m2')
    !      CALL RINOUT ('W KTOTVS',KTOTVS,MAXNRS)
    ierr = write_nc('KTOTVS',ktotvs,['MAXNRS'],[maxnrs],'Total volume of each flux tube','m3')
    !      CALL RINOUT ('W KVOL2 ',KVOL2 ,MAXNKS*MAXNRS)
    ierr = write_nc('KVOL2',kvol2,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Alternate calculation of KVOLS using polygon data (Currently KVOLS=KVOL2)','m2')
    !      CALL RINOUT ('W KAREA2',KAREA2,MAXNKS*MAXNRS)
    ierr = write_nc('KAREA2',karea2,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Alternate calculation of KAREA using polygon data (currently KAREA=KAREA2)','m3')
    !      CALL RINOUT ('W KTOTA2',KTOTA2,MAXNRS)
    ierr = write_nc('KTOTA2',ktota2,['MAXNRS'],[maxnrs],'Total area of each flux tube using polygon data','m2')
    !      CALL RINOUT ('W KTOTV2',KTOTV2,MAXNRS)
    ierr = write_nc('KTOTV2',ktotv2,['MAXNRS'],[maxnrs],'Total volume of each flux tube using polygon data','m3')
    !      CALL RINOUT ('W KBFS  ',KBFS  ,MAXNKS*MAXNRS)
    ierr = write_nc('KBFS',kbfs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Magnetic field ratio in each cell')

    !      CALL RINOUT ('W KINS  ',KINS  ,MAXNKS*MAXNRS)
    ierr = write_nc('KINS',kins,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Inward step probability')
    !      CALL RINOUT ('W KSMAXS',KSMAXS,MAXNRS)
    ierr = write_nc('KSMAXS',ksmaxs,['MAXNRS'],[maxnrs],'S max value for each ring','m')
    !      CALL RINOUT ('W KPMAXS',KPMAXS,MAXNRS)
    ierr = write_nc('KPMAXS',kpmaxs,['MAXNRS'],[maxnrs],'P max value for each ring','m')


    !
    ! Background plasma 
    !
    !      CALL RINOUT ('W KTEBS ',KTEBS ,MAXNKS*MAXNRS)
    ierr = write_nc('KTEBS',ktebs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Electron temperature for each cell','eV')
    !      CALL RINOUT ('W KTIBS ',KTIBS ,MAXNKS*MAXNRS)
    ierr = write_nc('KTIBS',ktibs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Ion temperature for each cell','eV')
    !      CALL RINOUT ('W KNBS  ',KNBS  ,MAXNKS*MAXNRS)
    ierr = write_nc('KNBS',knbs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Plasma density for each cell','/m3')
    !      CALL RINOUT ('W KVHS  ',KVHS  ,MAXNKS*MAXNRS)
    ierr = write_nc('KVHS',kvhs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Plasma flow velocity','m/s')
    !      CALL RINOUT ('W KES   ',KES   ,MAXNKS*MAXNRS)
    ierr = write_nc('KES',kes,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Electric field','V/m')



    !      CALL RINOUT ('W KFIZS ',KFIZS ,MAXNKS*MAXNRS*(MAXIZS+1))
    ierr = write_nc('KFIZS',kfizs,['MAXNKS  ','MAXNRS  ','MAXIZSP1'],[maxnks,maxnrs,maxizs+1],'Impurity ionization rate')
    !      CALL RINOUT ('W KINDS ',KINDS ,MAXNKS*MAXNRS)
    ierr = write_nc('KINDS',kinds,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Inward distance to next cell center','m')
    !c
    !c     jdemod - write cross field width of cells - v 46
    !c
    !      CALL RINOUT ('W KOUTDS',KOUTDS,MAXNKS*MAXNRS)
    ierr = write_nc('KOUTDS',koutds,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Outward distance to next cell center','m')


    !      CALL RINOUT ('W DISTIN',distin ,MAXNKS*MAXNRS)
    !
    ! Being written twice -jdefix
    !
    !ierr = write_nc('DISTIN',distin,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Inward distance to cell boundary','m')
    !c
    !      CALL RINOUT ('W DISTOU',distout,MAXNKS*MAXNRS)
    !ierr = write_nc('DISTOT',distout,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Outward distance to cell boundary','m')
    !
    !
    !      CALL RINOUT ('W KFORDS',KFORDS,MAXNKS*MAXNRS)
    ierr = write_nc('KFORDS',kfords,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Distance forward along field line to next cell center','m')
    !c
    !c     More geometry data
    !c
    !      CALL RINOUT ('W KBACDS',KBACDS,MAXNKS*MAXNRS)
    ierr = write_nc('KBACDS',kbacds,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Distance backward along field line to last cell center','m')
    !      CALL RINOUT ('W COSALI',COSALI,MAXNKS*MAXNRS)
    ierr = write_nc('COSALI',cosali,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Inward cosine from cell center to next inward ring')
    !      CALL RINOUT ('W COSALO',COSALO,MAXNKS*MAXNRS)
    ierr = write_nc('COSALO',cosalo,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Outward cosine from cell center to next outward ring')

    !      CALL RINOUT ('W DISTIN',DISTIN,MAXNKS*MAXNRS)
    ierr = write_nc('DISTIN',distin,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Inward distance to cell boundary','m')
    !      CALL RINOUT ('W DISTOU',DISTOUT,MAXNKS*MAXNRS)
    ierr = write_nc('DISTOU',distout,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Outward distance to cell boundary','m')


    !      CALL RINOUT ('W OKTEBS',OKTEBS,MAXNKS*MAXNRS)
    ierr = write_nc('OKTEBS',oktebs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Old KTEBS from last SOL iteration','eV')
    !      CALL RINOUT ('W OKTIBS',OKTIBS,MAXNKS*MAXNRS)
    ierr = write_nc('OKTIBS',oktibs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Old KTIBS from last SOL iteration','eV')
    !      CALL RINOUT ('W OKNBS ',OKNBS ,MAXNKS*MAXNRS)
    ierr = write_nc('OKNBS',oknbs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Old KNBS from last SOL iteration','eV')
    !      CALL RINOUT ('W OKVHS ',OKVHS ,MAXNKS*MAXNRS)
    ierr = write_nc('OKVHS',okvhs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Old KVHS from last SOL iteration','m/s')
    !      CALL RINOUT ('W OKES  ',OKES  ,MAXNKS*MAXNRS)
    ierr = write_nc('OKES',okes,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Old KES from last SOL iteration','V/m?')


    !      CALL RINOUT ('W KFEGS ',KFEGS ,MAXNKS*MAXNRS)
    ierr = write_nc('KFEGS',kfegs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Electron temperature gradient force component')
    !      CALL RINOUT ('W KFIGS ',KFIGS ,MAXNKS*MAXNRS)
    ierr = write_nc('KFIGS',kfigs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Ion tempterature gradient force component')

    !c
    !c     Average force arrays 
    !c
    !c
    !c     Coulomb logarithm for the Drift-Kinetic Model: CIOPTR
    !c
    !      Reiser model for transport
    !c
    !      WRITE(8) Coulomb_log
    ierr = write_nc('COULOMB_LOG',coulomb_log,'Value used for coulomb logarithm')
    !      CALL RINOUT ('W Fcell ',Fcell ,MAXNKS*MAXNRS*MAXIZS)
    ierr = write_nc('FCELL',fcell,['MAXNKS','MAXNRS','MAXIZS'],[maxnks,maxnrs,maxizs],'Average force in cell (reiser)')
    !      CALL RINOUT ('W Fthi  ',Fthi  ,MAXNKS*MAXNRS*MAXIZS)
    ierr = write_nc('FTHI',fthi,['MAXNKS','MAXNRS','MAXIZS'],[maxnks,maxnrs,maxizs],'Average ion thermal force (reiser)')
    !      CALL RINOUT ('W Ffi   ',Ffi   ,MAXNKS*MAXNRS*MAXIZS)
    ierr = write_nc('FFI',ffi,['MAXNKS','MAXNRS','MAXIZS'],[maxnks,maxnrs,maxizs],'Average ion friction force (reiser)')
    !      CALL RINOUT ('W Fvbg  ',Fvbg  ,MAXNKS*MAXNRS*MAXIZS)
    ierr = write_nc('FVBG',fvbg,['MAXNKS','MAXNRS','MAXIZS'],[maxnks,maxnrs,maxizs],'Average background flow velocity (reiser)')
    !      CALL RINOUT ('W DIFF  ',DIFF  ,MAXNKS*MAXNRS*MAXIZS)
    ierr = write_nc('DIFF',diff,['MAXNKS','MAXNRS','MAXIZS'],[maxnks,maxnrs,maxizs],'Average velocity diffusive step (reiser)')
    !      CALL RINOUT ('W VELavg',VELavg,MAXNKS*MAXNRS*MAXIZS)
    ierr = write_nc('VELAVG',velavg,['MAXNKS','MAXNRS','MAXIZS'],[maxnks,maxnrs,maxizs],'Average particle velocity (reiser)')

    !c
    !c     Background data at plates
    !c
    !      CALL RINOUT ('W RP    ',RP    ,MAXNDS)
    ierr = write_nc('RP',rp,['MAXNDS'],[maxnds],'R coordinate of target elements','m')
    !      CALL RINOUT ('W ZP    ',ZP    ,MAXNDS)
    ierr = write_nc('ZP',zp,['MAXNDS'],[maxnds],'Z coordinate of target elements','m')
    !      call iinout ('W IDDS  ',idds  ,maxnrs*2)
    ierr = write_nc('IDDS',idds,['MAXNRS','2     '],[maxnrs,2],'Target index by ring and target number')
    !      call rinout ('W PSITAR',psitarg,maxnrs*2)
    ierr = write_nc('PSITAR',psitarg,['MAXNRS','2     '],[maxnrs,2],'PSIn values of target elements')
    !      CALL RINOUT ('W KTEDS ',KTEDS ,MAXNDS)
    ierr = write_nc('KTEDS',kteds,['MAXNDS'],[maxnds],'Electron temperature at the targets','eV')
    !      CALL RINOUT ('W KTIDS ',KTIDS ,MAXNDS)
    ierr = write_nc('KTIDS',ktids,['MAXNDS'],[maxnds],'Ion temperature at the targets','eV')

    !
    ! not used in OUT - delete? - jdefix
    !
    !      CALL RINOUT ('W KTI3LS',KTI3LS,MAXNDS)
    !ierr = write_nc('KTI3LS',kti3ls,['MAXNDS'],[maxnds],'Ion temperature at 3X input decay length of source','eV')
    !      CALL RINOUT ('W KTINJ ',KTINJ ,MAXNDS)
    !ierr = write_nc('KTINJ',ktinj,['MAXNDS'],[maxnds],'Ion temperature at particle injection location','eV')

    !      CALL RINOUT ('W KNDS  ',KNDS  ,MAXNDS)
    ierr = write_nc('KNDS',knds,['MAXNDS'],[maxnds],'Density at target elements','/m3')
    !      CALL RINOUT ('W KFEDS ',KFEDS ,MAXNDS)
    ierr = write_nc('KFEDS',kfeds,['MAXNDS'],[maxnds],'Electron temperature gradient force at target')
    !      CALL RINOUT ('W KFIDS ',KFIDS ,MAXNDS)
    ierr = write_nc('KFIDS',kfids,['MAXNDS'],[maxnds],'Ion temperature gradient force at target')
    !      CALL RINOUT ('W KEDS  ',KEDS  ,MAXNDS)
    ierr = write_nc('KEDS',keds,['MAXNDS'],[maxnds],'Electric field at target elements','V/m')
    !c
    !      CALL RINOUT ('W KVDS  ',KVDS  ,MAXNDS)
    ierr = write_nc('KVDS',kvds,['MAXNDS'],[maxnds],'Plasma flow velocity at target','m/s')
    !c
    !      call rinout ('W HEATF ',targfluxdata,(maxnds+3)*4*4)      
    ierr = write_nc('HEATF',targfluxdata,['MAXNDSP3','4       ','4       '],[maxnds+3,4,4],'Heat flux at target')
    !c
    !      call rinout ('W KPREDB',kpredbar,maxnds*3*2)
    ierr = write_nc('KPREDB',kpredbar,['MAXNDS','3     ','2     '],[maxnds,3,2],'Near plate retention predictor values')
    !      CALL RINOUT ('W KFLUX ',KFLUX ,MAXNDS)
    ierr = write_nc('KFLUX',kflux,['MAXNDS'],[maxnds],'Target element ion flux','/m2/s')
    !      CALL RINOUT ('W KENER ',KENER ,MAXNDS)
    ierr = write_nc('KENER',kener,['MAXNDS'],[maxnds],'Target element average impact energy','eV')
    !      CALL RINOUT ('W KYIELD',KYIELD,MAXNDS)
    ierr = write_nc('KYIELD',kyield,['MAXNDS'],[maxnds],'Target element impurity sputtering yield')
    !      CALL RINOUT ('W KFY   ',KFY   ,MAXNDS)
    ierr = write_nc('KFY',kfy,['MAXNDS'],[maxnds],'Target element flux*yield')

    !      CALL RINOUT ('W KRMAX ',KRMAX ,MAXNDS)
    ierr = write_nc('KRMAX',krmax,['MAXNDS'],[maxnds],'Max random number for particle launch')
    !      CALL RINOUT ('W KCUM  ',KCUM  ,MAXNDS)
    ierr = write_nc('KCUM',kcum,['MAXNDS'],[maxnds],'Cumulative launch probability')

    !      CALL RINOUT ('W DDS   ',DDS   ,MAXNDS)
    ierr = write_nc('DDS',dds,['MAXNDS'],[maxnds],'Length of target elements (DDS=DDS2)','m')
    !      CALL RINOUT ('W DDS2  ',DDS2  ,MAXNDS)
    ierr = write_nc('DDS2',dds2,['MAXNDS'],[maxnds],'Alternate target length calculation (DDS=DDS2)','m')

    !      CALL RINOUT ('W THETAS',THETAS,MAXNDS)
    ierr = write_nc('THETAS',thetas,['MAXNDS'],[maxnds],'Angle of target element','radians')
    !      CALL RINOUT ('W THETA2',THETAS2,MAXNDS)
    ierr = write_nc('THETA2',thetas2,['MAXNDS'],[maxnds],'Alternate target angle calculation (THETAS=THETA2)','radians')

    !c
    !      CALL RINOUT ('W COSTET',COSTET,MAXNDS)
    ierr = write_nc('COSTET',costet,['MAXNDS'],[maxnds],'Target element angle cosine')
    !      call rinout ('W RHOG  ',rhog  ,maxnks*maxnrs)
    ierr = write_nc('RHOG',rhog,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Grid RHO values for JET grids')
    !      call rinout ('W THETAG',thetag,maxnks*maxnrs)
    ierr = write_nc('THETAG',thetag,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Cross field transport coordinate')
    !      call rinout ('W HRO   ',hro   ,maxnrs*maxnks)
    !ierr = write_nc('HRO',hro,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'No idea - only read in for JET grids - not used')
    !      call rinout ('W HTETA ',hteta ,maxnrs*maxnks)
    !ierr = write_nc('HTETA',hteta,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'No idea - only read in for JET grids - not used')
    !      call rinout ('W BTS   ',bts   ,maxnrs*maxnks)
    ierr = write_nc('BTS',bts,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Toroidal field on grid','T')
    !      call rinout ('W PSIFL ',psifl ,maxnrs*maxnks)
    ierr = write_nc('PSIFL',psifl,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PSI value for each cell (when available-depends on grid)')
    !      call iinout ('W TAGDV ',tagdv ,maxnrs*maxnks)
    ierr = write_nc('TAGDV',tagdv,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Cell orthogonality tag')
    !c
    !c     Chi Squared data - for studying convergence on sol iteration - needs initialization
    !c
    !      CALL DINOUT ('W CHISQ1',CHISQ1,25)
    ierr = write_nc('CHISQ1',chisq1,['MAXPINITER'],[maxpiniter],'KTEBS chisq for SOL iteration')
    !      CALL DINOUT ('W CHISQ2',CHISQ2,25)
    ierr = write_nc('CHISQ2',chisq2,['MAXPINITER'],[maxpiniter],'KTIBS chisq for SOL iteration')
    !      CALL DINOUT ('W CHISQ3',CHISQ3,25)
    ierr = write_nc('CHISQ3',chisq3,['MAXPINITER'],[maxpiniter],'KNBS chisq for SOL iteration')
    !      CALL DINOUT ('W CHISQ4',CHISQ4,25)
    ierr = write_nc('CHISQ4',chisq4,['MAXPINITER'],[maxpiniter],'KVHS chisq for SOL iteration')
    !      CALL DINOUT ('W CHISQ5',CHISQ5,25)
    ierr = write_nc('CHISQ5',chisq5,['MAXPINITER'],[maxpiniter],'KES chisq for SOL iteration')

    !
    ! Eirene/neutral hydrogen code quantities
    !
    !      CALL RINOUT ('W PINATO',PINATOM ,MAXNKS*MAXNRS)
    ierr = write_nc('PINATO',pinatom,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'H atom density','/m3')
    !      CALL RINOUT ('W PINION',PINION  ,MAXNKS*MAXNRS)
    ierr = write_nc('PINION',pinion,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'H ionization','/m3/s')
    !      CALL RINOUT ('W PINALP',PINALPHA,MAXNKS*MAXNRS)
    ierr = write_nc('PINALP',pinalpha,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Dalpha emission','ph/m3/s')
    !      CALL RINOUT ('W PINMOL',PINMOL  ,MAXNKS*MAXNRS)
    ierr = write_nc('PINMOL',pinmol,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'H2 molecule density','/m3')
    !      CALL RINOUT ('W PINZ0 ',PINZ0   ,MAXNKS*MAXNRS)
    ierr = write_nc('PINZ0',pinz0,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN neutral impurity density (if any)','/m3')
    !      CALL RINOUT ('W PININZ',PINIONZ ,MAXNKS*MAXNRS)
    ierr = write_nc('PININZ',pinionz,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN impurity ionization','/m3/s')
    !      CALL RINOUT ('W PINENA',PINENA  ,MAXNKS*MAXNRS)
    ierr = write_nc('PINENA',pinena,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN H atom temperature','eV')
    !      CALL RINOUT ('W PINENM',PINENM  ,MAXNKS*MAXNRS)
    ierr = write_nc('PINENM',pinenm,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN H2 molecule temperature','eV')
    !      CALL RINOUT ('W PINENZ',PINENZ  ,MAXNKS*MAXNRS)
    ierr = write_nc('PINENZ',pinenz,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN neutral impurity temperature','eV')
    !      CALL RINOUT ('W PINQI ',PINQI   ,MAXNKS*MAXNRS)
    ierr = write_nc('PINQI',pinqi,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN ion energy loss term')
    !      CALL RINOUT ('W PINQE ',PINQE   ,MAXNKS*MAXNRS)
    ierr = write_nc('PINQE',pinqe,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN electron energy loss term')
    !      CALL RINOUT ('W PINMP ',PINMP   ,MAXNKS*MAXNRS)
    ierr = write_nc('PINMP',pinmp,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN momentum loss term')
    !      CALL RINOUT ('W PINVDI',PINVDIST,3*14*MAXNKS*MAXNRS)
    ierr = write_nc('PINVDI',pinvdist,['3     ','14    ','MAXNKS','MAXNRS'],[3,14,maxnks,maxnrs],'PIN veocity distribution function')
    !      CALL RINOUT ('W PINREC',PINREC  ,MAXNKS*MAXNRS)
    ierr = write_nc('PINREC',pinrec,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'PIN calculated H+ recombination','/3/s')
    !      CALL RINOUT ('W PININF',PINIZ_INFO,MAXNRS*4)
    ierr = write_nc('PININF',piniz_info,['MAXNRS','4     '],[maxnrs,4],'PIN ionization and atom content on each ring - total and areal average')


    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 6')
    !
    !
    !
    !      call rinout ('W RVESM ',RVESM   ,2*MAXSEG)
    ierr = write_nc('RVESM',rvesm,['MAXSEG','2     '],[maxseg,2],'R coordinates of vessel wall segment end points','m')
    !      call rinout ('W ZVESM ',ZVESM   ,2*MAXSEG)
    ierr = write_nc('ZVESM',zvesm,['MAXSEG','2     '],[maxseg,2],'Z coordinates of vessel wall segment end points','m')
    !      call iinout ('W JVESM ',JVESM   ,MAXSEG)
    ierr = write_nc('JVESM',jvesm,['MAXSEG'],[maxseg],'Vessel wall segment type identification - wall/target/pfz etc')

    !      call rinout ('W FLUXHW',FLUXHW  ,MAXSEG)
    ierr = write_nc('FLUXHW',fluxhw,['MAXSEG'],[maxseg],'H atom+molecule flux to walls','/m2/s')
    !      call rinout ('W FLXHW2',FLXHW2  ,MAXSEG)
    ierr = write_nc('FLXHW2',flxhw2,['MAXSEG'],[maxseg],'H atom+ion flux to walls','/m2/s')
    !      call rinout ('W FLXHW3',FLXHW3  ,MAXSEG)
    ierr = write_nc('FLXHW3',flxhw3,['MAXSEG'],[maxseg],'Impurity wall sputtered flux','/m2/s')
    !      call rinout ('W FLXHW4',FLXHW4  ,MAXSEG)
    ierr = write_nc('FLXHW4',flxhw4,['MAXSEG'],[maxseg],'Impurity redeposited wall flux','/m2/s')
    !      call rinout ('W FLXHW5',FLXHW5  ,MAXSEG)
    ierr = write_nc('FLXHW5',flxhw5,['MAXSEG'],[maxseg],'Average energy of H atoms hitting walls','eV')
    !      call rinout ('W FLXHW6',FLXHW6  ,MAXSEG)
    ierr = write_nc('FLXHW6',flxhw6,['MAXSEG'],[maxseg],'Flux of H atoms to the wall','/m2/s')
    !      call rinout ('W FLXHW7',FLXHW7  ,MAXSEG)
    ierr = write_nc('FLXHW7',flxhw7,['MAXSEG'],[maxseg],'Average energy of H molecules hitting walls','eV')
    !      call rinout ('W FLXHW8',FLXHW8  ,MAXSEG)
    ierr = write_nc('FLXHW8',flxhw8,['MAXSEG'],[maxseg],'PIN calculated H ion fluxes to the walls','/m2/s')
    !
    !
    ! Leave out
    !
    !      CALL RINOUT ('W HWALKS',HWALKS  ,MAXNWS*2)
    !ierr = write_nc('HWALKS',hwalks,['MAXNWS','2     '],[maxnws,2],'H atom particle tracks','m')
    !      CALL rINOUT ('W SOLTE ',solte,maxnks*msolpt+msolpt+1)
    !ierr = write_nc('SOLTE',solte,['MAXNKS         ','MSOLPTPMSOLPTP1'],[maxnks,msolpt+msolpt+1],'High resolution Te background along one ring')
    !      CALL rINOUT ('W SOLTI ',solti,maxnks*msolpt+msolpt+1)
    !ierr = write_nc('SOLTI',solti,['MAXNKS         ','MSOLPTPMSOLPTP1'],[maxnks,msolpt+msolpt+1],'High resolution Ti background along one ring')
    !      CALL rINOUT ('W SOLNE ',solne,maxnks*msolpt+msolpt+1)
    !ierr = write_nc('SOLNE',solne,['MAXNKS         ','MSOLPTPMSOLPTP1'],[maxnks,msolpt+msolpt+1],'High resolution ne background along one ring')
    !      CALL rINOUT ('W SOLVEL',solvel,maxnks*msolpt+msolpt+1)
    !ierr = write_nc('SOLVEL',solvel,['MAXNKS         ','MSOLPTPMSOLPTP1'],[maxnks,msolpt+msolpt+1],'High resolution vb background along one ring')
    !      CALL rINOUT ('W SOLCOR',solcor,maxnks*msolpt+msolpt+1)
    !ierr = write_nc('SOLCOR',solcor,['MAXNKS         ','MSOLPTPMSOLPTP1'],[maxnks,msolpt+msolpt+1],'High resolution S value along one ring')
    !
    !C
    !c     Leakage data
    !c
    !      CALL RINOUT ('W CLEAKS',cleaks  ,MAXpts)
    ierr = write_nc('CLEAKS',cleaks,['MAXPTS'],[maxpts],'Parallel leakage boundaries for recording SOL transport','m')
    !      CALL RINOUT ('W CLEAKN',cleakn  ,MAXpts*(maxizs+1))
    ierr = write_nc('CLEAKN',cleakn,['MAXPTS  ','MAXIZSP1'],[maxpts,maxizs+1],'Number of particles passing leakage boundaries by charge state')
    !      call rinout ('W LEAKPS',cleakpos,maximp*2)
    ierr = write_nc('LEAKPS',cleakpos,['MAXIMP','2     '],[maximp,2],'Start R,Z locations of particles entering the confined plasma','m')
    !c
    !c     More arrays related to leakage results
    !c
    !      call rinout ('W ncore ',ncore,maxnks*maxnrs)
    ierr = write_nc('NCORE',ncore,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Start cells of particles entering core')
    !      call rinout ('W nedge ',nedge,maxnks*maxnrs)
    ierr = write_nc('NEDGE',nedge,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Start cells of particles entering edge')
    !      call rinout ('W ntrap ',ntrap,maxnks*maxnrs)
    ierr = write_nc('NTRAP',ntrap,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Start cells of particles entering PFZ')
    !      call rinout ('W ndivt ',ndivert,maxnks*maxnrs)
    ierr = write_nc('NDIVT',ndivert,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Start cells of particles entering divertor')
    !      call rinout ('W nmsol ',nmsol,maxnks*maxnrs)
    ierr = write_nc('NMSOL',nmsol,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Start cells of particles entering main SOL')

    !c
    !c     NOTE: Not all of wtsource is being read or written - if 
    !c           this is ever needed both write and read routines
    !c           need to be adjusted.
    !c
    !      call rinout ('W WTSOU ',wtsource,maxpts*maxnrs*4*5)
    ierr = write_nc('WTSOU',wtsource,['MAXPTS','MAXNRS','4     ','6     '],[maxpts,maxnrs,4,6],'Wall and target particle source data')
    !      call rinout ('W WTDEP ',wtdep,maxpts*(maxpts+1)*3)
    ierr = write_nc('WTDEP',wtdep,['MAXPTS  ','MAXPTSP1','3       '],[maxpts,maxpts+1,3],'Wall and target particle deposition data')
    !      call rinout ('W TSOUR ',targsrc,3*4)
    ierr = write_nc('TSOUR',targsrc,['3','4'],[3,4],'Integrated target source data')
    !      call rinout ('W TLEAK ',targleak,3*4)
    ierr = write_nc('TLEAK',targleak,['3','4'],[3,4],'Integrated leakage data over targets')
    !      call rinout ('W WSOUR ',wallsrc,5*3)
    ierr = write_nc('WSOUR',wallsrc,['5','3'],[5,3],'Integrated wall source data')
    !      call rinout ('W WLEAK ',wallleak,5*3)
    ierr = write_nc('WLEAK',wallleak,['5','3'],[5,3],'Integtated leakage data over walls')


    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 7')


    !c
    !c     Store any EDGE2D BG data that has been read in
    !c
    !c
    !      Saving externally loaded plasma solution

    !if (cre2d.eq.1.or.cre2d.eq.2.or.cre2d.eq.3.or.cre2d.eq.5) then
    !   !        call rinout ('W E2D N ',e2dnbs,maxnks*maxnrs)
    !   ierr = write_nc('E2D N',e2dnbs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma density')
    !   !        call rinout ('W E2D TE',e2dtebs,maxnks*maxnrs)
    !   ierr = write_nc('E2D TE',e2dtebs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma Te')
    !   !        call rinout ('W E2D TI',e2dtibs,maxnks*maxnrs)
    !   ierr = write_nc('E2D TI',e2dtibs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma Ti')
    !   !        call rinout ('W E2D VB',e2dvhs,maxnks*maxnrs)
    !   ierr = write_nc('E2D VB',e2dvhs,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma Vb')
    !   !        call rinout ('W E2D E ',e2des,maxnks*maxnrs)
    !   ierr = write_nc('E2D E',e2des,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma Efield')
    !   !        call rinout ('W E2D I ',e2dion,maxnks*maxnrs)
    !   ierr = write_nc('E2D I',e2dion,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma ionization source')
    !   !        call rinout ('W E2D A ',e2datom,maxnks*maxnrs)
    !   ierr = write_nc('E2D A',e2datom,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma neutral H density')
    !   !        call rinout ('W E2D TA',e2dtarg,maxnrs*8*2)
    !   ierr = write_nc('E2D TA',e2dtarg,['MAXNRS','8     ','2     '],[maxnrs,8,2],'External plasma target conditions')
    !   !        call rinout ('W E2D GP',e2dgpara,maxnks*maxnrs)
    !   ierr = write_nc('E2D GP',e2dgpara,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma parallel fluxes')
    !   !        call rinout ('W E2D GD',e2dgdown,maxnks*maxnrs)
    !   ierr = write_nc('E2D GD',e2dgdown,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma down fluxes')
    !   !        call rinout ('W E2D G ',e2dflux,(maxnks+1)*maxnrs)
    !   ierr = write_nc('E2D G',e2dflux,['MAXNKSP1','MAXNRS  '],[maxnks+1,maxnrs],'External plasma cell center flux')
    !   !c
    !   !        call rinout ('W E2D VE',e2dbvel,(maxnks+1)*maxnrs)
    !   ierr = write_nc('E2D VE',e2dbvel,['MAXNKSP1','MAXNRS  '],[maxnks+1,maxnrs],'External plasma flow velocity')
    !   !c
    !   !        call rinout ('W E2D Z0',e2dz0,maxnks*maxnrs)
    !   ierr = write_nc('E2D Z0',e2dz0,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma neutral impurity density')
    !   !        call rinout ('W E2D RC',e2dhrec,maxnks*maxnrs)
    !   ierr = write_nc('E2D RC',e2dhrec,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma H recombination')
    !   !        call rinout ('W E2D RC',e2drec,maxnks*maxnrs)
    !   ierr = write_nc('E2D RC',e2drec,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma impurity recombination')
    !   !c
    !   !        call rinout ('W E2D CX',e2dcxrec,maxnks*maxnrs)
    !   ierr = write_nc('E2D CX',e2dcxrec,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma CX recombination')
    !   !c
    !   if (cre2dizs.gt.0) then
    !      !          call rinout ('W E2D NZ ',e2dnzs,maxnks*maxnrs*(maxe2dizs+1))
    ! jdemod - save the fluid code impurity results to the netcdf file
          ierr = write_nc('E2D NZ',e2dnzs,['MAXNKS     ','MAXNRS     ','MAXE2DIZSP1'],[maxnks,maxnrs,maxe2dizs+1],'External plasma impurity density')
          ierr = write_nc('E2D VZ',e2dvzs,['MAXNKS     ','MAXNRS     ','MAXE2DIZSP1'],[maxnks,maxnrs,maxe2dizs+1],'External plasma impurity velocity')
    !      !          call rinout ('W E2D PW',e2dpowls,maxnks*maxnrs*(maxe2dizs+1))
    !      ierr = write_nc('E2D PW',e2dpowls,['MAXNKS     ','MAXNRS     ','MAXE2DIZSP1'],[maxnks,maxnrs,maxe2dizs+1],'External plasma radiated power')
    !      !          call rinout ('W E2D LI',e2dlines,maxnks*maxnrs*(maxe2dizs+1))
    !      ierr = write_nc('E2D LI',e2dlines,['MAXNKS     ','MAXNRS     ','MAXE2DIZSP1'],[maxnks,maxnrs,maxe2dizs+1],'External plasma line radiation')
    !
    !      !c
    !   endif
    !   !C
    !   !c     Store Data related to transport coefficient calculations
    !   !c
    !endif

    !
    ! Dperp/Xperp extractor values
    !
    !      call rinout ('W DPERP ',DPERP,maxnrs)
    ierr = write_nc('DPERP',dperp,['MAXNRS'],[maxnrs],'Extracted Dperp')
    !      call rinout ('W DPERPO',odperp,maxnrs)
    ierr = write_nc('DPERPO',odperp,['MAXNRS'],[maxnrs],'Extracted Dperp '//OUTER)
    !      call rinout ('W DPERPI',idperp,maxnrs)
    ierr = write_nc('DPERPI',idperp,['MAXNRS'],[maxnrs],'Extracted Dperp '//INNER)
    !      call rinout ('W XPERP ',xPERPt,maxnrs)
    ierr = write_nc('XPERP',xperpt,['MAXNRS'],[maxnrs],'Extracted Xperp_total')
    !      call rinout ('W XPERPO',oxperpt,maxnrs)
    ierr = write_nc('XPERPO',oxperpt,['MAXNRS'],[maxnrs],'Extracted Xperp_total '//OUTER)
    !      call rinout ('W XPERPI',ixperpt,maxnrs)
    ierr = write_nc('XPERPI',ixperpt,['MAXNRS'],[maxnrs],'Extracted Xperp_total '//INNER)
    !      call rinout ('W XPI   ',chiperpi,maxnrs)
    ierr = write_nc('XPI',chiperpi,['MAXNRS'],[maxnrs],'Extracted Xperp ion')
    !      call rinout ('W XPI  O',ochiperpi,maxnrs)
    ierr = write_nc('XPI  O',ochiperpi,['MAXNRS'],[maxnrs],'Extracted Xperp ion '//OUTER)
    !      call rinout ('W XPI  I',ichiperpi,maxnrs)
    ierr = write_nc('XPI  I',ichiperpi,['MAXNRS'],[maxnrs],'Extracted Xperp ion '//INNER)
    !      call rinout ('W XPE   ',chiperpe,maxnrs)
    ierr = write_nc('XPE',chiperpe,['MAXNRS'],[maxnrs],'Extracted Xperp electron')
    !      call rinout ('W XPE  O',ochiperpe,maxnrs)
    ierr = write_nc('XPE  O',ochiperpe,['MAXNRS'],[maxnrs],'Extracted Xperp electron '//OUTER)
    !      call rinout ('W XPE  I',ichiperpe,maxnrs)
    ierr = write_nc('XPE  I',ichiperpe,['MAXNRS'],[maxnrs],'Extracted Xperp electron '//INNER)


    !c
    !c     GEometric quantities originally associated with transport extractor
    !c
    !      call rinout ('W RC OUT',rcouter,maxnrs)
    ierr = write_nc('RC OUT',rcouter,['MAXNRS'],[maxnrs],'R coordinate outer midplane','m')
    !      call rinout ('W RC IN ',rcinner,maxnrs)
    ierr = write_nc('RC IN',rcinner,['MAXNRS'],[maxnrs],'R coordinate inner midplane','m')
    !      call rinout ('W ZC OUT',zcouter,maxnrs)
    ierr = write_nc('ZC OUT',zcouter,['MAXNRS'],[maxnrs],'Z coordinate outer midplane','m')
    !      call rinout ('W ZC IN ',zcinner,maxnrs)
    ierr = write_nc('ZC IN',zcinner,['MAXNRS'],[maxnrs],'Z coordinate inner midplane','m')
    !      call rinout ('W MIDIST',middist,maxnrs*2)
    ierr = write_nc('MIDIST',middist,['MAXNRS','2     '],[maxnrs,2],'Midplane distance to separatrix 1=inner 2=outer','m')



    !      call rinout ('W GRADNE',gradn,maxnks*maxnrs)
    ierr = write_nc('GRADNE',gradn,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Radial density gradient')
    !      call rinout ('W GRADTE',gradte,maxnks*maxnrs)
    ierr = write_nc('GRADTE',gradte,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Radial Te gradient')
    !      call rinout ('W GRADTI',gradti,maxnks*maxnrs)
    ierr = write_nc('GRADTI',gradti,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Radial Ti gradient')
    !      call rinout ('W E2DGNE',e2dgradn,maxnks*maxnrs)
    ierr = write_nc('E2DGNE',e2dgradn,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma - radial density gradient')
    !      call rinout ('W E2DGTE',e2dgradte,maxnks*maxnrs)
    ierr = write_nc('E2DGTE',e2dgradte,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma - radial Te gradient')
    !      call rinout ('W E2DGTI',e2dgradti,maxnks*maxnrs)
    ierr = write_nc('E2DGTI',e2dgradti,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'External plasma - radial Ti gradient')

    !c
    !c     IF Line profile data was calculated - store it 
    !c
    !c
    !      write(8) line_profile_opt
    ierr = write_nc('LINE_PROFILE_OPT',line_profile_opt,'Calculate doppler shifted line profile option')
    if (line_profile_opt.ne.0) then 
       !c slmod begin
       !c         Descriptor needs to be 8 characters long or generates
       !c         a runtime error in R8INOUT. -SL, 07/10/2011
       !          write(8) lp_wave,lp_instrument_width,lp_bin_width,lp_robs,lp_zobs,lp_theta,lp_dtheta
       ierr = write_nc('LP_WAVE',lp_wave,'Line profile wavelength')
       ierr = write_nc('LP_INSTRUMENT_WIDTH',lp_instrument_width,'Line profile instrument width')
       ierr = write_nc('LP_BIN_WIDTH',lp_bin_width,'Line profile - wavelength bin width')
       ierr = write_nc('LP_ROBS',lp_robs,'Line profile - R observation position','m')
       ierr = write_nc('LP_ZOBS',lp_zobs,'Line profile - Z observation position','m')
       ierr = write_nc('LP_THETA',lp_theta,'Line profile - Viewing angle','degrees')
       ierr = write_nc('LP_DTHETA',lp_dtheta,'Line profile - viewing cone angular width','degrees')
       !c
       !c          CALL R8INOUT ('W LP',line_profile,max_lp_bins*2+1)
       !c slmod end
       !          CALL R8INOUT('W LP    ',line_profile,max_lp_bins*2+1)
       ierr = write_nc('LP',line_profile,['2MAX_LP_BINSP1'],[2*max_lp_bins+1],'Doppler shifted line profile of requested line')
       !          CALL R8INOUT('W MOD_LP',modified_line_profile,max_lp_bins*2+1)
       ierr = write_nc('MOD_LP',modified_line_profile,['2MAX_LP_BINSP1'],[max_lp_bins*2+1],'Modified line profile')
    endif


    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER BLOCK 8')


    !c
    !c     Store the pressure - from SOL option 22
    !c
    !c     IPP/11 - first arg must be 8 chars string
    !      call rinout ('W KPRESS',kpress,maxnks*maxnrs*2)
    ierr = write_nc('KPRESS',kpress,['MAXNKS','MAXNRS','2     '],[maxnks,maxnrs,2],'Pressure along field lines from SOL22')
    !      call rinout ('W KPRAD ',kprad,maxnks*maxnrs)
    ierr = write_nc('KPRAD',kprad,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Prad term from SOL22 (or input)')

    !c
    !c     Write out the global HC activation option
    !c
    !c
    !c ammod begin
    !c
    !      write(8) global_hc_follow_option
    ierr = write_nc('GLOBAL_HC_FOLLOW_OPTION',global_hc_follow_option,'HC following option')
    !c
    !c Write out HC module related quantities 
    !c
    if (global_hc_follow_option.ne.0) then  
       !c
       call hc_store_raw_data_netcdf
       !c
       !c     Added to raw file primarily for use by HC - code
       !c 
    endif

    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER HC')



    !
    ! Add far periphery data to netcdf file
    !
    call fp_write_netcdf
    
    call pr_trace('WRITE_NETCDF_OUTPUT','AFTER FP')


    !      call rinout ('W BRATIO',BRATIO,maxnks*maxnrs)
    ierr = write_nc('BRATIO',bratio,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Magnetic field ratio')

    !c
    !c ammod end
    !c
    !c
    !c     Write out subgrid data - at the least the option value
    !c     will be added to the raw file. If the subgrid was in use - 
    !c     its data will be saved and any storage used deallocated.  
    !c
    !
    ! jdemod - Leave out for now
    !
    ! call save_subgrid(8)

    !c
    !c     jdemod - March 2016 - version 45
    !c
    !c     Write out potential and drift related results
    !c     
    !      call rinout ('W POT',osmpot2,maxnks*maxnrs)
    ierr = write_nc('POT',osmpot2,['MAXNKSP2','MAXNRS  '],[maxnks+2,maxnrs],'Electric potential by cell','V')
    !      call rinout ('W E_RAD',e_rad,maxnks*maxnrs)
    ierr = write_nc('E_RAD',e_rad,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Radial electric field','V/m')
    !      call rinout ('W E_POL',e_pol,maxnks*maxnrs)
    ierr = write_nc('E_POL',e_pol,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Poloidal electric field','V/m')
    !      call rinout ('W ExB_R',exb_rad_drft,maxnks*maxnrs)
    ierr = write_nc('EXB_R',exb_rad_drft,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Radial drift velocity','m/s')
    !      call rinout ('W ExB_P',exb_pol_drft,maxnks*maxnrs)
    ierr = write_nc('EXB_P',exb_pol_drft,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'Poloidal drift velocity','m/s')


    !c
    !c
    !c     Temporarily Add the following (?) 
    !c
    !      jdemod - not needed 
    !c
    !      call rinout ('W FLUXES',fluxes,maxnks*maxnrs*16)
    !ierr = write_nc('FLUXES',fluxes,['MAXNKS','MAXNRS','16'],[maxnks,maxnrs,16],'Detailed flux comparison data for debugging External code vs SOL22')
    !

    IF (IMODE.EQ.1) THEN
       !      CALL RINOUT ('W LIMS  ',LIMS  ,MAXNKS*MAXNRS*(MAXIZS+2)*MAXNTS)
       ierr = write_nc('LIMS',lims,['MAXNKS  ','MAXNRS  ','MAXIZSP2','MAXNTS  '],[maxnks,maxnrs,maxizs+2,maxnts],'')
       !c
       !c
       !c slmod begin - new
    ENDIF

    slver = 3.6

    !      WRITE(8) slver
    ierr = write_nc('SLVER',slver,'Steve Lisgo (SL) output version identifier')
    !
    ! jdemod - the following contains a lot of mostly EIRENE related detailed debugging output that 
    !          is not generally needed so it is commented out for now. 
    !
    !
    !!      WRITE(8) MAXASD,MAXNAS,eirnpgdat,((eirpgdat(i1,i2),i2=1,MAXASD),i1=1,MAXNAS)
    !ierr = write_nc('MAXASD',maxasd,'')
    !ierr = write_nc('MAXNAS',maxnas,'')
    !ierr = write_nc('EIRNPGDAT',eirnpgdat,'')
    !!      WRITE(8) asc_ncell,MAXASC
    !ierr = write_nc('ASC_NCELL',asc_ncell,'')
    !ierr = write_nc('MAXASC',maxasc,'')
    !!      CALL IINOUT('W CELL  ',asc_cell  ,MAXASC)
    !ierr = write_nc('CELL',asc_cell,['MAXASC'],[maxasc],'')
    !!      CALL IINOUT('W REGION',asc_region,MAXASC)
    !ierr = write_nc('REGION',asc_region,['MAXASC'],[maxasc],'')
    !!      CALL IINOUT('W LINK  ',asc_link  ,MAXASC*4)
    !ierr = write_nc('LINK',asc_link,['MAXASC','4'],[maxasc,4],'')
    !!      CALL IINOUT('W GRID  ',asc_grid  ,MAXASC*2)
    !ierr = write_nc('GRID',asc_grid,['MAXASC','2'],[maxasc,2],'')
    !!      CALL IINOUT('W NVP   ',asc_nvp   ,MAXASC)
    !ierr = write_nc('NVP',asc_nvp,['MAXASC'],[maxasc],'')
    !!      CALL RINOUT('W RVP   ',asc_rvp   ,MAXASC*8)
    !ierr = write_nc('RVP',asc_rvp,['MAXASC','8'],[maxasc,8],'')
    !!c...  slver 2.1:
    !!      CALL RINOUT('W ZVP   ',asc_zvp   ,MAXASC*8)
    !ierr = write_nc('ZVP',asc_zvp,['MAXASC','8'],[maxasc,8],'')
    !!c.... slver 2.2, switched to 40 from 20:
    !!      CALL IINOUT('W NVERTX',ascnvertex,MAXASC)
    !ierr = write_nc('NVERTX',ascnvertex,['MAXASC'],[maxasc],'')
    !!      CALL RINOUT('W VERTEX',ascvertex ,MAXASC*40)       
    !ierr = write_nc('VERTEX',ascvertex,['MAXASC','40'],[maxasc,40],'')

    !
    ! Write out components of EIRENE calculated emission - might be useful
    !
    !c...  6.34:
    !      CALL RINOUT ('W PINLN1',pinline(1,1,1,H_BALPHA),MAXNKS*MAXNRS*6)
    tmp_pinline(:,:,:) = pinline(:,:,:,H_BALPHA)
    ierr = write_nc('PINLN1',tmp_pinline,['MAXNKS','MAXNRS','6     '],[maxnks,maxnrs,6],'')
    !      CALL RINOUT ('W PINLN2',pinline(1,1,1,H_BBETA ),MAXNKS*MAXNRS*6)
    tmp_pinline(:,:,:) = pinline(:,:,:,H_BBETA)
    ierr = write_nc('PINLN2',tmp_pinline,['MAXNKS','MAXNRS','6     '],[maxnks,maxnrs,6],'')
    !      CALL RINOUT ('W PINLN3',pinline(1,1,1,H_BGAMMA),MAXNKS*MAXNRS*6)
    tmp_pinline(:,:,:) = pinline(:,:,:,H_BGAMMA)
    ierr = write_nc('PINLN3',tmp_pinline,['MAXNKS','MAXNRS','6     '],[maxnks,maxnrs,6],'')




    !ierr = write_nc('TEST_ARRAY',test_array,['5','5','5'],[5,5,5],'test array with same dimension names')




    !!      WRITE(8) pincode
    !ierr = write_nc('PINCODE',pincode,'')
    !!      CALL RINOUT('W PINMOI',pinmoi,MAXNKS*MAXNRS)
    !ierr = write_nc('PINMOI',pinmoi,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')
    !!      CALL RINOUT('W OSMCDE',osmcde,MAXNKS*MAXNRS)
    !ierr = write_nc('OSMCDE',osmcde,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')
    !!      CALL RINOUT('W OSMCDI',osmcdi,MAXNKS*MAXNRS)      
    !ierr = write_nc('OSMCDI',osmcdi,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')
    !!      CALL RINOUT('W OSMCVE',osmcve,MAXNKS*MAXNRS)
    !ierr = write_nc('OSMCVE',osmcve,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')
    !!c...  slver 1.6:
    !!      CALL RINOUT('W OSMCVI',osmcvi,MAXNKS*MAXNRS)      
    !ierr = write_nc('OSMCVI',osmcvi,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')
    !!      WRITE(8) tarshift(IKLO),tarshift(IKHI)
    !ierr = write_nc('TARSHIFT',tarshift,'')
    !!c...  slver 1.7:
    !!      WRITE(8) te_mult_o,ti_mult_o,n_mult_o,te_mult_i,ti_mult_i,n_mult_i
    !ierr = write_nc('TE_MULT_O',te_mult_o,'')
    !ierr = write_nc('TI_MULT_O',ti_mult_o,'')
    !ierr = write_nc('N_MULT_O',n_mult_o,'')
    !ierr = write_nc('TE_MULT_I',te_mult_i,'')
    !ierr = write_nc('TI_MULT_I',ti_mult_i,'')
    !ierr = write_nc('N_MULT_I',n_mult_i,'')


    !      WRITE(8) lpdatsw
    ierr = write_nc('LPDATSW',lpdatsw,'LP probe data interpretation switch 0=density 1=Jsat')
    !      CALL RINOUT('W SEPDIS',sepdist ,MAXNDS)
    ierr = write_nc('SEPDIS',sepdist,['MAXNDS'],[maxnds],'Distance along target to separatrix','m')
    !      CALL RINOUT('W SEPDI2',sepdist2,MAXNDS)
    ierr = write_nc('SEPDI2',sepdist2,['MAXNDS'],[maxnds],'Alternate calculation of distance to separatrix','m')

    !
    !c...  slver 1.8:
    !!      CALL IINOUT('W PRBNUM',prb_num ,NUMPRB) 
    !ierr = write_nc('PRBNUM',prb_num,['NUMPRB'],[numprb],'')
    !!      CALL RINOUT('W PRBRHO',prb_rho ,MAXPRB*NUMPRB)
    !ierr = write_nc('PRBRHO',prb_rho,['MAXPRB','NUMPRB'],[maxprb,numprb],'')
    !!      CALL RINOUT('W PRBTE ',prb_te  ,MAXPRB*NUMPRB)
    !ierr = write_nc('PRBTE',prb_te,['MAXPRB','NUMPRB'],[maxprb,numprb],'')
    !!      CALL RINOUT('W PRBTI ',prb_ti  ,MAXPRB*NUMPRB)
    !ierr = write_nc('PRBTI',prb_ti,['MAXPRB','NUMPRB'],[maxprb,numprb],'')
    !!      CALL RINOUT('W PRBNE ',prb_ne  ,MAXPRB*NUMPRB)
    !ierr = write_nc('PRBNE',prb_ne,['MAXPRB','NUMPRB'],[maxprb,numprb],'')
    !!      CALL RINOUT('W PRBR  ',prb_r   ,MAXPRB*NUMPRB)
    !ierr = write_nc('PRBR',prb_r,['MAXPRB','NUMPRB'],[maxprb,numprb],'')
    !!c...  slver 1.9:
    !!      CALL RINOUT('W PRBZ  ',prb_z   ,MAXPRB*NUMPRB)
    !ierr = write_nc('PRBZ',prb_z,['MAXPRB','NUMPRB'],[maxprb,numprb],'')
    !!      WRITE(8) eirnres
    !ierr = write_nc('EIRNRES',eirnres,'')
    !!c...  slver 2.4 (removed for 2.17):
    !!      CALL RINOUT('W EIRRES',eirres,6*7*MAXPINITER)
    !ierr = write_nc('EIRRES',eirres,['6','7','MAXPINITER'],[6,7,maxpiniter],'')
    !!c...  slver 2.5 & 3.4::
    !!      CALL RINOUT('W PINPLO',pinploss ,MAXNKS*MAXNRS*NMOMCHA)
    !ierr = write_nc('PINPLO',pinploss,['MAXNKS ','MAXNRS ','NMOMCHA'],[maxnks,maxnrs,nmomcha],'')
    !!      WRITE(8) eirnpuff,eirpmode
    !ierr = write_nc('EIRNPUFF',eirnpuff,'')
    !ierr = write_nc('EIRPMODE',eirpmode,'')
    !!c...  slver 2.7:
    !!      CALL RINOUT('W PUFF  ',eirpuff  ,MAXNAS*MAXPUFF)
    !ierr = write_nc('PUFF',eirpuff,['MAXNAS','MAXPUFF'],[maxnas,maxpuff],'')
    !!      CALL RINOUT('W BGK   ',pinbgk   ,MAXNKS*MAXNRS*MAXBGK)
    !ierr = write_nc('BGK',pinbgk,['MAXNKS','MAXNRS','MAXBGK'],[maxnks,maxnrs,maxbgk],'')
    !!c...  slver 3.0:
    !!      WRITE(8) MAXASD,MAXNAS3,eirnspdat,((eirspdat(i1,i2),i2=1,MAXASD),i1=1,MAXNAS3)
    !ierr = write_nc('MAXASD',maxasd,'')
    !ierr = write_nc('MAXNAS3',maxnas3,'')
    !ierr = write_nc('EIRNSPDAT',eirnspdat,'')
    !!      WRITE(8) asc_3dmode
    !ierr = write_nc('ASC_3DMODE',asc_3dmode,'')
    !!      CALL IINOUT('W NVP3D ',asc_nvp3D ,MAXASC3D)
    !ierr = write_nc('NVP3D',asc_nvp3d,['MAXASC3D'],[maxasc3d],'')
    !!      CALL IINOUT('W LIMK3D',asc_link3D,MAXASC3D*6)
    !ierr = write_nc('LIMK3D',asc_link3d,['MAXASC3D','6       '],[maxasc3d,6],'')
    !!      CALL RINOUT('W ZMIN3D',asc_zmin3D,MAXASC3D)
    !ierr = write_nc('ZMIN3D',asc_zmin3d,['MAXASC3D'],[maxasc3d],'')
    !!      CALL RINOUT('W ZMAX3D',asc_zmax3D,MAXASC3D)
    !ierr = write_nc('ZMAX3D',asc_zmax3d,['MAXASC3D'],[maxasc3d],'')
    !!      CALL RINOUT('W XVP3D ',asc_xvp3D ,MAXASC3D*6*8)
    !ierr = write_nc('XVP3D',asc_xvp3d,['MAXASC3D','6       ','8       '],[maxasc3d,6,8],'')
    !!      CALL RINOUT('W YVP3D ',asc_yvp3D ,MAXASC3D*6*8)
    !ierr = write_nc('YVP3D',asc_yvp3d,['MAXASC3D','6       ','8       '],[maxasc3d,6,8],'')
    !!c...  slver 3.1:
    !!      CALL RINOUT('W ZVP3D ',asc_zvp3D ,MAXASC3D*6*8)
    !ierr = write_nc('ZVP3D',asc_zvp3d,['MAXASC3D','6       ','8       '],[maxasc3d,6,8],'')
    !!c...  slver 3.2:
    !!c     Resetting the 3D data reading in OUT.
    !!c...  slver 3.3:
    !!      WRITE(8) ascncut,MAXASC3D,MAXASCDAT
    !ierr = write_nc('ASCNCUT',ascncut,'')
    !ierr = write_nc('MAXASC3D',maxasc3d,'')
    !ierr = write_nc('MAXASCDAT',maxascdat,'')
    !!      WRITE(8) eirnsdtor,(eirsdtor(i1),eirsdvol(i1),i1=1,eirnsdtor)
    !ierr = write_nc('EIRNSDTOR',eirnsdtor,'')

    !
    ! Leave out of file for now
    !
    !DO i1 = 2, eirnsdtor
    !   i2 = (i1-1)*MAXBGK+1
    !   ! writes slices of the array as separate entries
    !   !        CALL RINOUT('W BGKTOR',pinbgk(1,1,i2),MAXNKS*MAXNRS*MAXBGK)
    !
    !   ierr = write_nc('BGKTOR',pinbgk(1,['1'],[1],'')
    !   !c...  6.13:
    !ENDDO


    !!      WRITE(8) eirniontime,MAXIONTIME,MAXBIN
    !ierr = write_nc('EIRNIONTIME',eirniontime,'')
    !ierr = write_nc('MAXIONTIME',maxiontime,'')
    !ierr = write_nc('MAXBIN',maxbin,'')
    !!c...  6.14:
    !!      CALL RINOUT('W IONTIM',eiriontime,MAXIONTIME*(20+MAXBIN*3))
    !ierr = write_nc('IONTIM',eiriontime,['MAXIONTIME','20PMAXBIN','3'],[maxiontime,20+maxbin,3],'')
    !!      WRITE(8) MAXASD2,MAXASS,asc_ncell,ascncut
    !ierr = write_nc('MAXASD2',maxasd2,'')
    !ierr = write_nc('MAXASS',maxass,'')
    !ierr = write_nc('ASC_NCELL',asc_ncell,'')
    !ierr = write_nc('ASCNCUT',ascncut,'')
    !!c...  6.15:
    !!      WRITE(8)(asc_vol(i1),(ascdata(i1,i2),i2=1,5),((pinasd(i1,i2,i3,1),pinasd(i1,i2,i3,2),i2=1,MAXASD2),i3=1,MAXASS ),i1=1,asc_ncell*ascncut+1+eirnpgdat),999999
    !ierr = write_nc('ASC_VOL',asc_vol,'')
    !!      WRITE(8) eirntally,MAXNKS,MAXNRS,MAXTALLY
    !ierr = write_nc('EIRNTALLY',eirntally,'')
    !ierr = write_nc('MAXNKS',maxnks,'')
    !ierr = write_nc('MAXNRS',maxnrs,'')
    !ierr = write_nc('MAXTALLY',maxtally,'')
    !!c...  6.16:
    !!      WRITE(8)((eirtally(i1,i2),i2=1,4),((pinalgv(ik,ir,i1),ik=1,nks(ir)),ir=1,nrs),i1=1,eirntally),999999
    !ierr = write_nc('EIRTALLY',eirtally,'')
    !!c...  6.22:
    !!      WRITE(8) eirzaa,eirtorfrac,eirsrcmul
    !ierr = write_nc('EIRZAA',eirzaa,'')
    !ierr = write_nc('EIRTORFRAC',eirtorfrac,'')
    !ierr = write_nc('EIRSRCMUL',eirsrcmul,'')
    !!      WRITE(8) osmns28,8,osm_nfnt,3,MAXFNT
    !ierr = write_nc('OSMNS28',osmns28,'')
    !ierr = write_nc('8',8,'')
    !ierr = write_nc('OSM_NFNT',osm_nfnt,'')
    !ierr = write_nc('3',3,'')
    !ierr = write_nc('MAXFNT',maxfnt,'')
    !!c...  6.25:
    !!
    !      jdemod - writing osm28 data - not included in netcdf output for now
    !
    !      WRITE(8) ((osms28(i1,i2),i2=1,8),i1=1,osmns28),((osm_ionfnt(i1,i2),i2=1,3),i1=1,osm_nfnt)
    !  ierr = write_nc('',,'')
    !      WRITE(8) eirnstrdat,eirnstrai,eirnatmi,eirnmoli,MAXSTRDAT,MAXSTRATA,100
    !ierr = write_nc('EIRNSTRDAT',eirnstrdat,'')
    !ierr = write_nc('EIRNSTRAI',eirnstrai,'')
    !ierr = write_nc('EIRNATMI',eirnatmi,'')
    !ierr = write_nc('EIRNMOLI',eirnmoli,'')
    !ierr = write_nc('MAXSTRDAT',maxstrdat,'')
    !ierr = write_nc('MAXSTRATA',maxstrata,'')
    !ierr = write_nc('100',100,'')
    !c...  6.26:
    !
    !      jdemod - writing eirene strata data - not included in netcdf output for now
    !
    !      WRITE(8) (((eirstrdat(i1,i2,i3),i3=1,100      ),i2=1,MAXSTRATA),eirstrlis(i1)        ,i1=1,MAXSTRDAT),999999
    !      ierr = write_nc('',,'')
    !
    !      CALL IINOUT('W IKBOUN',ikbound,MAXNRS*2)
    !ierr = write_nc('IKBOUN',ikbound,['MAXNRS','2'],[maxnrs,2],'')
    !c...  6.30:
    !c     This is temporary.  I want the ability to store ionisation data
    !c     for a greater number of strata, so I made room in the PINDATA array
    !c     by not storing the atom and molecule density information for each
    !c     stratum - SL, Sep 19, 2002:
    !      CALL IINOUT('W IKBOUN',ikfluid,MAXNRS*2)
    !ierr = write_nc('IKBOUN',ikfluid,['MAXNRS','2'],[maxnrs,2],'')
    !DO i1 = H_ION1, H_ION1+11
    !   !        CALL IINOUT('W PINDAT',pindata(1,1,i1),MAXNKS*MAXNRS)
    !   ierr = write_nc('PINDAT',pindata(1,['1'],[1],'')
    !   !c...  6.28:
    !ENDDO
    !DO i1 = 1, 3
    !   !        CALL IINOUT('W PINICP',pinioncomp(1,1,i1),MAXNKS*MAXNRS)
    !   ierr = write_nc('PINICP',pinioncomp(1,['1'],[1],'')
    !   !c...  6.29:
    !ENDDO
    !c...  6.33:
    !      WRITE(8) eirntorseg
    !ierr = write_nc('EIRNTORSEG',eirntorseg,'')


    !
    ! jdemod - these values define some X,Y boxes and ciopte is the ion injection option
    !        - not generally needed in output
    !
    !      WRITE(8) ciopte,cxsc,cysc,cxsca,cysca,cxscb,cyscb
    !ierr = write_nc('CIOPTE',ciopte,'')
    !ierr = write_nc('CXSC',cxsc,'')
    !ierr = write_nc('CYSC',cysc,'')
    !ierr = write_nc('CXSCA',cxsca,'')
    !ierr = write_nc('CYSCA',cysca,'')
    !ierr = write_nc('CXSCB',cxscb,'')
    !ierr = write_nc('CYSCB',cyscb,'')

    !c...  slver = 3.5: *TEMP*
    !!      CALL RINOUT('W EIRPH1',eirpho1,MAXNKS*MAXNRS)
    !ierr = write_nc('EIRPH1',eirpho1,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')
    !!      CALL RINOUT('W EIRPH2',eirpho2,MAXNKS*MAXNRS)
    !ierr = write_nc('EIRPH2',eirpho2,['MAXNKS','MAXNRS'],[maxnks,maxnrs],'')

    !c...  6.41:

    !
    ! jdemod - debug flags - logicals - netcdf does not support storing logicals natively ... would need conversion to 
    !                                   integer or byte or other compatible data type ... so leave out for now - not really needed
    !
    !      WRITE(8) debugv,cstepv
    !ierr = write_nc('DEBUGV',debugv,'')
    !ierr = write_nc('CSTEPV',cstepv,'')
    !
    ! 
    !
    !IF (debugv) then 
    !   !CALL RINOUT ('W SDVS',sdvs,MAXNKS*MAXNRS*(MAXIZS+2))
    !
    ! jdemod - save ddvs to netcdf since it is now calculated all the time
    ierr = write_nc('DDVS',ddvs,['MAXNKS  ','MAXNRS  ','MAXIZSP2'],[maxnks,maxnrs,MAXIZS+2],'Impurity ion average velocity','m/s')
    !
    !endif
    !c...  slver 3.6:      
    !
    ! jdemod - not sure what this is - not included in netcdf for now
    !
    !
    !IF (ALLOCATED(wall_flx)) THEN 
    !   !        WRITE(8) 1
    !   ierr = write_nc('1',1,'')
    !   !        WRITE(8) wall_n,1.0  ! this 1.0 is a version number
    !   ierr = write_nc('WALL_N',wall_n,'')
    !   ierr = write_nc('1.0  ! THIS 1.0 IS A VERSION NUMBER',1.0  ! this 1.0 is a version number,'')
    !   !        WRITE(8) MAXNBLK,MAXNATM,MAXNMOL,MAXNSRC,MAXNLAUNCH
    !   ierr = write_nc('MAXNBLK',maxnblk,'')
    !   ierr = write_nc('MAXNATM',maxnatm,'')
    !   ierr = write_nc('MAXNMOL',maxnmol,'')
    !   ierr = write_nc('MAXNSRC',maxnsrc,'')
    !   ierr = write_nc('MAXNLAUNCH',maxnlaunch,'')
    !   !        WRITE(8) wall_flx
    !   ierr = write_nc('WALL_FLX',wall_flx,'')
    !ELSE
    !   !        WRITE(8) 0
    !   ierr = write_nc('0',0,'')
    !ENDIF
    !c...  6.14 (end of file flag):
    !      WRITE(8) 123456789
    !ierr = write_nc('123456789',123456789,'')
    
    ! sazmod - Write out some of the blobby transport stuff.
    ierr = write_nc('fblob',fblob,'Blobby transport: Blob frequency')
    ierr = write_nc('div_vr_fact',div_vr_fact,'Blobby transport: Divertor region transport multiplier')
    !ierr = write_nc('in_blob_switch',in_blob_switch,'Blobby transport: Turn off parallel transport in blob')
    ierr = write_nc('pinch_correlation_time',pinch_correlation_time,'Blobby transport: Blob correlation time')
    ierr = write_nc('pinch_pdf',pinch_pdf_data,['MAXPTS', '2     '], [maxpts, 2], 'Blobby transport: Blob vr distribution', 'm/s')
    

    
    call pr_trace('WRITE_NETCDF_OUTPUT','BEFORE CLOSE')


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


  subroutine hc_store_raw_data_netcdf

    Use ComHC ! HC constants.
    Use HC_Init_DIV_Diag ! Included to re-set hc_state, hc_density, hc_output, Number_HC_Species.
    implicit none

    integer :: ierr
    !c
    ! ammod begin.
    !c     Addition of hydrocarbon data handling.
    !c
    !      call iinout ('W HC_OPT',hc_follow_option,1)
    ierr = write_nc('HC_OPT',hc_follow_option,'HC follow option')
    !      call iinout ('W HC_HIG',hc_higher_hcs_option,1)
    ierr = write_nc('HC_HIG',hc_higher_hcs_option,'HC higher hydrocarbon option')
    !      call iinout ('W HC_PRI',hc_evolution_model_primary,1)
    ierr = write_nc('HC_PRI',hc_evolution_model_primary,'HC evolution primary model option')
    !      call iinout ('W HC_SEC',hc_evolution_model_secondary,1)
    ierr = write_nc('HC_SEC',hc_evolution_model_secondary,'HC evolution secondary model option')

    !
    ! jdemod - hc_state_list is a character array - not currently in my netcdf interface
    !
    !      call rinout ('W HC_STA', HC_State_List,Number_HC_Species)
    !ierr = write_nc('HC_STA',hc_state_list,['NUMBER_HC_SPECIES'],[number_hc_species],'HC state list')

    !      call dinout ('W HC_DEN', HC_Density,maxnks*maxnrs*(Number_HC_Species))
    ierr = write_nc('HC_DEN',hc_density,['MAXNKS           ','MAXNRS           ','NUMBER_HC_SPECIES'],[maxnks,maxnrs,number_hc_species],'HC density')
    !      call dinout ('W H_DEN', H_Density,maxnks*maxnrs*(Number_H_Species))
    ierr = write_nc('H_DEN',h_density,['MAXNKS      ','MAXNRS      ','NUM_H_STATES'],[maxnks,maxnrs,num_h_states],'H Density (from HC)')
    !      call iinout ('W HC_OUT', HC_Output_List,Number_HC_Species+1)
    ierr = write_nc('HC_OUT',hc_output_list,['NUMBER_HC_SPECIESP1'],[number_hc_species+1],'HC species output list')
    !      call rinout ('W HC_WLK', HC_Walks,Max_Number_Walks*2)
    ierr = write_nc('HC_WLK',hc_walks,['MAX_NUMBER_WALKS','2               '],[max_number_walks,2],'HC particle tracks')
    !
    ! jdemod - remove the extra storage in hc_factor (0 and -1) indices
    !
    !      call rinout ('W HC_FACTA', HC_Factor_A,Number_HC_Species)
    ierr = write_nc('HC_FACTA',hc_factor_a,['NUMBER_HC_SPECIES'],[number_hc_species],'HC scale factor A')
    !      call rinout ('W HC_FACTB', HC_Factor_B,Number_HC_Species)
    ierr = write_nc('HC_FACTB',hc_factor_b,['NUMBER_HC_SPECIES'],[number_hc_species],'HC scale factor B')
    !      call rinout ('W HC_TIZS_CH', HC_TIZS_CH,maxnks*maxnrs*2)
    ierr = write_nc('HC_TIZS_CH',hc_tizs_ch,['MAXNKS','MAXNRS','2     '],[maxnks,maxnrs,2],'CH species ionization by cell')

    !
    !c      call rinout ('W HC_TIZS_C2', HC_TIZS_C2,
    !c     >  maxnks*maxnrs*2)
    !c
    !c      CALL RINOUT ('W FYTOT ',FYTOT ,1)      
    !c
    !c
    !c     End addition of hydrocarbon data to BIN unit 8 file for read in OUT.
    ! ammod end.

    return 
  end subroutine hc_store_raw_data_netcdf


end module divimp_netcdf
