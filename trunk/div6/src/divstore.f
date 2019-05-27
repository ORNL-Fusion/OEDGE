c     -*Fortran*-
C
      SUBROUTINE STORE (TITLE,desc,NIZS,JOB,EQUIL,
     >                  FACTA,FACTB,ITER,NITERS)
      use debug_options
      use subgrid
c slmod begin
      use mod_divimp
c slmod end
      use mod_fp_data
      use divimp_netcdf
      use mod_params
      use mod_comtor
      use mod_cadas
      use mod_cneut2
      use mod_cgeom
      use mod_dynam1
      use mod_dynam3
      use mod_dynam4
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
      IMPLICIT  NONE
C     INCLUDE   "PARAMS"
c     include    'params'
      CHARACTER TITLE*(*),desc*(*),JOB*(*),EQUIL*(*)
c
      INTEGER   NIZS,ITER,NITERS
      REAL      FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)
C
C  *********************************************************************
C  *                                                                   *
C  *  STORE: WRITE RESULTS OF DIVIMP RUN ON UNFORMATTED UNIT 8         *
C  *         UNIT 8 SHOULD BE REWOUND BEFORE STORE IS CALLED FOR THE   *
C  *         FIRST TIME.  STORE MAY BE CALLED SEVERAL TIMES, ONCE FOR  *
C  *         EACH ITERATION.                                           *
C  *                                   CHRIS FARRELL    MARCH 1989     *
C  *********************************************************************
C
C     INCLUDE   "COMTOR"
c     include    'comtor'
c     include    'cadas'
C     INCLUDE   "CNEUT2"
c     include    'cneut2'
C     INCLUDE   "CGEOM"
c     include    'cgeom'
C     INCLUDE   "DYNAM1"
c     include    'dynam1'
C     INCLUDE   "DYNAM3"
c     include    'dynam3'
C     INCLUDE   "DYNAM4"
c     include    'dynam4'
C     INCLUDE   "PINDATA"
c     include    'pindata'
c
c     include    'divxy'
c
c     include    'grbound'
c
c     include    'cedge2d'
c
c     include    'transcoef'
c
c     include    'cioniz'
c     include    'promptdep'
c     include    'reiser_com'
c     include    'line_profile' 
c     include    'hc_global_opts'
c     include    'driftvel'

      ! Add periphery variables
      !include 'fperiph_com'

C
      INTEGER IR,IZ,IT,IN
c slmod begin 
c     INCLUDE 'diagvel'
c     INCLUDE 'slcom'



      INTEGER      i1,i2,i3,ik,i
      REAL         slver

c
      call pr_trace('STORE','START')

c
c slmod end
C
c
c     Version number 
c     
      write(8) verson
c     
c     Basic information - title, job description,
c                         equilibrium, shot number and time slice
c     
c      write(0,*) 'writing:',len(title),len(desc),len(job),len(equil),
c     >                ishot,tslice

      WRITE  (8) TITLE,desc,JOB,EQUIL,ISHOT,TSLICE
c
c     Write out the global parameters used to write the file
c     - this is the first step in parameterizing the read statements in 
c       OUT and removing the dependency on having identical parameter 
c       values in both the DIVIMP and OUT compiles
c
      write(8) MAXNKS,MAXNRS,MAXNDS,MAXNGS,MAXIZS,MAXINS,                
     >  MAXIMP,ISECT,MAXNTS,MAXNOC,MAXNWS,MAXNXS,MAXNYS,MAXVMF,         
     >  MAXTHE,MAXSN,MAXPTS,MAXPLRP,MSOLPT,MAXADS,      
     >  MAXGXS,MAXGYS,MAXCH3,MAXIXS,MAXIYS,MAXSEG,maxplts,MAXNFLA,
     >  maxpiniter,mbufle,mbufx,mves,maxe2dizs

c     
c     Simulation values
c     
      
      write(8) ITER,NITERS,NIZS,NTS,CRMB,CRMI,CIZB,CION,IMODE,
     >         NITERSOL
c     
c     Geometry values 
c     
      write(8) R0,Z0,RXP,ZXP,NRS,
     >         RMIN,RMAX,ZMIN,ZMAX,DR,DZ,NDS,NDSIN,NXS,NYS,
     >         IRCENT,IRSEP,IRWALL,IRTRAP,IKT,IKREF,IKTO,IKTI,
     >         irsep2,irwall2,nrs2,ndsin2,ndsin3,irtrap2,
     >         nves,nvesm,nvesp,inmid,oumid,refct,CIRHR,NPOLYP,
     >         cxnear,cynear
c
      CALL IINOUT ('W NKS    ',nks ,maxnrs)
c     
c     Scaling factors  
c     
      write(8) QTIM,FSRATE,ABSFAC,absfac_neut,CBPHI,CK0,CK0I
c     
      CALL RINOUT ('W FACTA  ',facta ,maxizs+2)
      CALL RINOUT ('W FACTB  ',factb ,maxizs+2)
      CALL RINOUT ('W DWELTS ',dwelts,maxizs+2)
c slmod begin
      IF (IMODE.EQ.1) THEN
        CALL RINOUT ('W DWELFS ',dwelfs,maxnts)
      ELSE
        CALL RINOUT ('W DWELFS ',dwelfs,1     )
      ENDIF
c
c      CALL RINOUT ('W DWELFS ',dwelfs,maxnts)
c slmod end
      CALL RINOUT ('W KALPHS ',kalphs,maxizs)
      CALL RINOUT ('W KBETAS ',kbetas,maxizs)
c     
c     DIVIMP Options
c     
      write(8) CIOPTF,cdatopt,cneur,cgridopt,cre2d,cre2dizs,
     >         xygrid
c     
c     Result numbers
c     
      write(8) nleakcore,cleaksn
c     
c     ADAS       
c     
      write(8) useridh,iyearh,useridz,iyearz
c     
c     
c     The following are the piece-wise wall specifications
c     from the GA15 routines.
c     
      WRITE(8) ioncpts,ionwpts
c     
c     Wall definitions
c     
      CALL RINOUT ('W RIW    ',riw ,maxpts)
      CALL RINOUT ('W ZIW    ',ziw ,maxpts)
      CALL RINOUT ('W RCW    ',rcw ,maxpts)
      CALL RINOUT ('W ZCW    ',zcw ,maxpts)
      CALL RINOUT ('W RW     ',rw  ,maxpts)
      CALL RINOUT ('W ZW     ',zw  ,maxpts)
      CALL RINOUT ('W RVES   ',rves,maxpts)
      CALL RINOUT ('W ZVES   ',zves,maxpts)
c     
c     GA15 workspace 
c     
      CALL IINOUT ('W IWINDW ',iwindw,2*maxpts)
      CALL RINOUT ('W IWWORK ',iwwork,4*maxpts)
      CALL RINOUT ('W IWTDUM ',iwtdum,maxpts)
      CALL RINOUT ('W IWXDUM ',iwxdum,maxpts)
      CALL RINOUT ('W IWYDUM ',iwydum,maxpts)
c     
      CALL IINOUT ('W ICINDW ',icindw,2*maxpts)
      CALL RINOUT ('W ICWORK ',icwork,4*maxpts)
      CALL RINOUT ('W ICTDUM ',ictdum,maxpts)
      CALL RINOUT ('W ICXDUM ',icxdum,maxpts)
      CALL RINOUT ('W ICYDUM ',icydum,maxpts)
c     
c     
c     The following relate to the wall definition as 
c     calculated in DIVIMP.
c     
      WRITE(8) wlwall1,wlwall2,wltrap1,wltrap2,wallpts
c
c
c      WRITE (8) VERSON,ITER,NITERS,NIZS,TITLE,JOB,NTS,R0,Z0,RXP,ZXP,
c     >          NRS,(NKS(IR),IR=1,MAXNRS),RMIN,RMAX,ZMIN,ZMAX,DR,DZ,NDS,
c     >          (FACTA(IZ),FACTB(IZ),DWELTS(IZ),IZ=-1,MAXIZS),
c     >          (DWELFS(IT),IT=1,MAXNTS),QTIM,FSRATE,NXS,NYS,ABSFAC,
c     >          CRMB,CRMI,CIZB,CION,IMODE,IRCENT,IRSEP,IRWALL,IRTRAP,
c     >          IKT,NDSIN,(KALPHS(IZ),KBETAS(IZ),IZ=1,MAXIZS),IKREF,
c     >          CK0,CIRHR,IKTO,IKTI,NPOLYP,ISHOT,TSLICE,EQUIL,
c     >          NITERSOL,CIOPTF,CBPHI,cdatopt,cneur,nves,nleakcore,
c     >          irsep2,irwall2,nrs2,ndsin2,ndsin3,cgridopt,
c     >          irtrap2,xygrid,cleaksn,cxnear,cynear,refct,
c     >          nvesm,nvesp,inmid,oumid,useridh,iyearh,useridz,iyearz,
c     >          cre2d,cre2dizs,
c
c               The following are the piece-wise wall specifications
c               from the GA15 routines.
c
c     >          ioncpts,ionwpts,(riw(in),ziw(in),rcw(in),zcw(in),
c     >          rw(in),zw(in),rves(in),zves(in),in=1,maxpts),
c     >          ((iwindw(it,in),icindw(it,in),it=1,2),in=1,maxpts),
c     >          (iwwork(in),icwork(in),in=1,4*maxpts),
c     >          (iwtdum(in),iwxdum(in),iwydum(in),in=1,maxpts),
c     >          (ictdum(in),icxdum(in),icydum(in),in=1,maxpts),
c
c               The following relate to the wall definition as
c               calculated in DIVIMP.
c
c     >          wlwall1,wlwall2,wltrap1,wltrap2,wallpts
c
C
c
      call pr_trace('STORE','PART 1')


      WRITE (6,9001) NXS,NYS,NRS,NDS,NIZS,NTS,
     >  MAXNXS,MAXNYS,MAXNRS,MAXNDS,MAXIZS,MAXNTS,
     >  TITLE,JOB,ITER,IMODE,refct,maxseg,nvesm,nvesp
c     >  TITLE,JOB,ITER,IMODE,refct
c
      CALL RINOUT ('W POWLS ',POWLS ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('W LINES ',LINES ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('W HPOWLS',HPOWLS,MAXNKS*MAXNRS*2)
      CALL RINOUT ('W HLINES',HLINES,MAXNKS*MAXNRS*2)
      CALL RINOUT ('W TIZS  ',TIZS  ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('W ZEFFS ',ZEFFS ,MAXNKS*MAXNRS*3)
      CALL RINOUT ('W WALLS ',WALLS ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('W DEPS  ',DEPS  ,MAXNDS*MAXIZS)
      CALL RINOUT ('W NEROS ',NEROS ,MAXNDS*5)
      CALL RINOUT ('W PRDEPS',PROMPTDEPS,MAXNDS*9)
c
      CALL RINOUT ('W WALLSN',WALLSN,MAXPTS+1)
      CALL RINOUT ('W WALLSE',WALLSE,MAXPTS+1)
      CALL RINOUT ('W WALLSEI',WALLSE_I,MAXPTS+1)
      CALL RINOUT ('W WALLSI',WALLSI,MAXPTS+1)
      CALL RINOUT ('W WALLSIL',WALLSIL,MAXPTS+1)
c
      CALL RINOUT ('W WALLPT',WALLPT,MAXPTS*32)
      CALL RINOUT ('W RS    ',RS    ,MAXNKS*MAXNRS)
      CALL RINOUT ('W ZS    ',ZS    ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KSB   ',KSB   ,(MAXNKS+1)*MAXNRS)
      CALL RINOUT ('W KPB   ',KPB   ,(MAXNKS+1)*MAXNRS)
c
c     The storing of these arrays needed to be customized
c     because of a likely size mismatch between the
c     array in DIVIMP and that in OUT.
c
c      CALL IINOUT ('W IKXYS ',IKXYS ,MAXNXS*MAXNYS)
c      CALL IINOUT ('W IRXYS ',IRXYS ,MAXNXS*MAXNYS)
c      CALL IINOUT ('W IFXYS ',IFXYS ,MAXNXS*MAXNYS)
c
      CALL IINOUT2 ('W IKXYS ',IKXYS ,MAXNXS,MAXNYS,MAXNXS,MAXNYS)
      CALL IINOUT2 ('W IRXYS ',IRXYS ,MAXNXS,MAXNYS,MAXNXS,MAXNYS)
      CALL IINOUT2 ('W IFXYS ',IFXYS ,MAXNXS,MAXNYS,MAXNXS,MAXNYS)
c
      CALL IINOUT ('W IKDS  ',IKDS  ,MAXNDS)
      CALL IINOUT ('W IRDS  ',IRDS  ,MAXNDS)
C
      CALL IINOUT ('W IKINS ',IKINS ,MAXNKS*MAXNRS)
      CALL IINOUT ('W IKOUTS',IKOUTS,MAXNKS*MAXNRS)
      CALL IINOUT ('W IRINS ',IRINS ,MAXNKS*MAXNRS)
      CALL IINOUT ('W IROUTS',IROUTS,MAXNKS*MAXNRS)
C
      CALL IINOUT ('W KORY  ',KORY  ,MAXNKS*MAXNRS)
      CALL IINOUT ('W KORPG ',KORPG ,MAXNKS*MAXNRS)
      CALL IINOUT ('W NVERTP',NVERTP,MAXNKS*MAXNRS)
      CALL RINOUT ('W RVERTP',RVERTP,5*MAXNKS*MAXNRS)
      CALL RINOUT ('W ZVERTP',ZVERTP,5*MAXNKS*MAXNRS)
C
      CALL CHECK_DDLIM(nizs,3)
c
      CALL DINOUT ('W DDLIMS',DDLIMS,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL DINOUT ('W DDTS  ',DDTS  ,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RINOUT ('W ELIMS ',ELIMS ,MAXNKS*3*(MAXIZS+2))
      CALL RINOUT ('W WALKS ',WALKS ,MAXNWS*2)
c
      call rinout ('W CHEM D',chemden,maxnks*maxnrs)
      call rinout ('W CHEMIZ',chemizs,maxnks*maxnrs)
C
      CALL RINOUT ('W KKS   ',KKS   ,MAXNRS)
      CALL RINOUT ('W KSS   ',KSS   ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KPS   ',KPS   ,MAXNKS*MAXNRS)

      ! jdemod - not used at all - must be development leftover
      !CALL RINOUT ('W KNORMS',KNORMS,MAXNKS*MAXNRS)

      CALL RINOUT ('W KPERPS',KPERPS,MAXNKS*MAXNRS)

      ! jdemod - not used at all - must be development leftover
      !CALL RINOUT ('W KCURVS',KCURVS,MAXNKS*MAXNRS)

      CALL RINOUT ('W KVOLS ',KVOLS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KAREAS',KAREAS,MAXNKS*MAXNRS)
      CALL RINOUT ('W KTOTAS',KTOTAS,MAXNRS)
      CALL RINOUT ('W KTOTVS',KTOTVS,MAXNRS)
      CALL RINOUT ('W KVOL2 ',KVOL2 ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KAREA2',KAREA2,MAXNKS*MAXNRS)
      CALL RINOUT ('W KTOTA2',KTOTA2,MAXNRS)
      CALL RINOUT ('W KTOTV2',KTOTV2,MAXNRS)
      CALL RINOUT ('W KBFS  ',KBFS  ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KINS  ',KINS  ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KSMAXS',KSMAXS,MAXNRS)
      CALL RINOUT ('W KPMAXS',KPMAXS,MAXNRS)
      CALL RINOUT ('W KTEBS ',KTEBS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KTIBS ',KTIBS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KNBS  ',KNBS  ,MAXNKS*MAXNRS)
c
      CALL RINOUT ('W KFIZS ',KFIZS ,MAXNKS*MAXNRS*(MAXIZS+1))
c
c      CALL RINOUT ('W KFSSMO',KFSSMOD,MAXNKS*MAXNRS)
c
      CALL RINOUT ('W KINDS ',KINDS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KOUTDS',KOUTDS,MAXNKS*MAXNRS)
c
c     jdemod - write cross field width of cells - v 46
c             removed on version 47 since already 
c             written below
c
c      CALL RINOUT ('W DISTIN',distin ,MAXNKS*MAXNRS)
c      CALL RINOUT ('W DISTOU',distout,MAXNKS*MAXNRS)
c
      CALL RINOUT ('W KFORDS',KFORDS,MAXNKS*MAXNRS)
      CALL RINOUT ('W KBACDS',KBACDS,MAXNKS*MAXNRS)
c
c     More geometry data
c
      CALL RINOUT ('W COSALI',COSALI,MAXNKS*MAXNRS)
      CALL RINOUT ('W COSALO',COSALO,MAXNKS*MAXNRS)
      CALL RINOUT ('W DISTIN',DISTIN,MAXNKS*MAXNRS)
      CALL RINOUT ('W DISTOU',DISTOUT,MAXNKS*MAXNRS)
c
      CALL RINOUT ('W OKTEBS',OKTEBS,MAXNKS*MAXNRS)
      CALL RINOUT ('W OKTIBS',OKTIBS,MAXNKS*MAXNRS)
      CALL RINOUT ('W OKNBS ',OKNBS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W OKVHS ',OKVHS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W OKES  ',OKES  ,MAXNKS*MAXNRS)
C
      CALL RINOUT ('W KFEGS ',KFEGS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KFIGS ',KFIGS ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KES   ',KES   ,MAXNKS*MAXNRS)
      CALL RINOUT ('W KVHS  ',KVHS  ,MAXNKS*MAXNRS)
c
      call pr_trace('STORE','PART 2')
c
c     Average force arrays 
c
c
c     Coulomb logarithm for the Drift-Kinetic Model: CIOPTR
c
      WRITE(8) Coulomb_log
c
      CALL RINOUT ('W Fcell ',Fcell ,MAXNKS*MAXNRS*MAXIZS)
      CALL RINOUT ('W Fthi  ',Fthi  ,MAXNKS*MAXNRS*MAXIZS)
      CALL RINOUT ('W Ffi   ',Ffi   ,MAXNKS*MAXNRS*MAXIZS)
      CALL RINOUT ('W Fvbg  ',Fvbg  ,MAXNKS*MAXNRS*MAXIZS)
      CALL RINOUT ('W DIFF  ',DIFF  ,MAXNKS*MAXNRS*MAXIZS)
      CALL RINOUT ('W VELavg',VELavg,MAXNKS*MAXNRS*MAXIZS)
c
c
c     Background data at plates
c
      CALL RINOUT ('W RP    ',RP    ,MAXNDS)
      CALL RINOUT ('W ZP    ',ZP    ,MAXNDS)
      call iinout ('W IDDS  ',idds  ,maxnrs*2)
      call rinout ('W PSITAR',psitarg,maxnrs*2)
      CALL RINOUT ('W KTEDS ',KTEDS ,MAXNDS)
      CALL RINOUT ('W KTIDS ',KTIDS ,MAXNDS)

      ! jdemod - these are not used in out at all 
      !          removed for now
      !CALL RINOUT ('W KTI3LS',KTI3LS,MAXNDS)
      !CALL RINOUT ('W KTINJ ',KTINJ ,MAXNDS)
c

      CALL RINOUT ('W KNDS  ',KNDS  ,MAXNDS)
      CALL RINOUT ('W KFEDS ',KFEDS ,MAXNDS)
      CALL RINOUT ('W KFIDS ',KFIDS ,MAXNDS)
      CALL RINOUT ('W KEDS  ',KEDS  ,MAXNDS)
      CALL RINOUT ('W KVDS  ',KVDS  ,MAXNDS)
c
      call rinout ('W HEATF ',targfluxdata,(maxnds+3)*4*4)      
c
      call rinout ('W KPREDB',kpredbar,maxnds*3*2)
c
      CALL RINOUT ('W KFLUX ',KFLUX ,MAXNDS)
      CALL RINOUT ('W KENER ',KENER ,MAXNDS)
      CALL RINOUT ('W KYIELD',KYIELD,MAXNDS)
      CALL RINOUT ('W KFY   ',KFY   ,MAXNDS)
      CALL RINOUT ('W KRMAX ',KRMAX ,MAXNDS)
      CALL RINOUT ('W KCUM  ',KCUM  ,MAXNDS)
      CALL RINOUT ('W DDS   ',DDS   ,MAXNDS)
      CALL RINOUT ('W THETAS',THETAS,MAXNDS)
      CALL RINOUT ('W DDS2  ',DDS2  ,MAXNDS)
      CALL RINOUT ('W THETA2',THETAS2,MAXNDS)
      CALL RINOUT ('W COSTET',COSTET,MAXNDS)
c
      call rinout ('W RHOG  ',rhog  ,maxnrs*maxnks)
      call rinout ('W THETAG',thetag,maxnrs*maxnks)
      call rinout ('W HRO   ',hro   ,maxnrs*maxnks)
      call rinout ('W HTETA ',hteta ,maxnrs*maxnks)
      call rinout ('W BTS   ',bts   ,maxnrs*maxnks)
      call rinout ('W PSIFL ',psifl ,maxnrs*maxnks)
      call iinout ('W TAGDV ',tagdv ,maxnrs*maxnks)
c
c     Chi Squared data
c
      CALL DINOUT ('W CHISQ1',CHISQ1,maxpiniter)
      CALL DINOUT ('W CHISQ2',CHISQ2,maxpiniter)
      CALL DINOUT ('W CHISQ3',CHISQ3,maxpiniter)
      CALL DINOUT ('W CHISQ4',CHISQ4,maxpiniter)
      CALL DINOUT ('W CHISQ5',CHISQ5,maxpiniter)
c
c      CALL DINOUT ('W CHISQ1',CHISQ1,25)
c      CALL DINOUT ('W CHISQ2',CHISQ2,25)
c      CALL DINOUT ('W CHISQ3',CHISQ3,25)
c      CALL DINOUT ('W CHISQ4',CHISQ4,25)
c      CALL DINOUT ('W CHISQ5',CHISQ5,25)
C
      CALL RINOUT ('W PINATO',PINATOM ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINION',PINION  ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINALP',PINALPHA,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINMOL',PINMOL  ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINZ0 ',PINZ0   ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PININZ',PINIONZ ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINENA',PINENA  ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINENM',PINENM  ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINENZ',PINENZ  ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINQI ',PINQI   ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINQE ',PINQE   ,MAXNKS*MAXNRS)
      CALL RINOUT ('W PINMP ',PINMP   ,MAXNKS*MAXNRS)
c
      CALL RINOUT ('W PINVDI',PINVDIST,3*14*MAXNKS*MAXNRS)
      CALL RINOUT ('W PINREC',PINREC  ,MAXNKS*MAXNRS)
c
      CALL RINOUT ('W PININF',PINIZ_INFO,MAXNRS*4)
c
      call pr_trace('STORE','PART 3')
C
      call rinout ('W RVESM ',RVESM   ,2*MAXSEG)
      call rinout ('W ZVESM ',ZVESM   ,2*MAXSEG)
      call iinout ('W JVESM ',JVESM   ,MAXSEG)
      call rinout ('W FLUXHW',FLUXHW  ,MAXSEG)
      call rinout ('W FLXHW2',FLXHW2  ,MAXSEG)
      call rinout ('W FLXHW3',FLXHW3  ,MAXSEG)
      call rinout ('W FLXHW4',FLXHW4  ,MAXSEG)
      call rinout ('W FLXHW5',FLXHW5  ,MAXSEG)
      call rinout ('W FLXHW6',FLXHW6  ,MAXSEG)
      call rinout ('W FLXHW7',FLXHW7  ,MAXSEG)
      call rinout ('W FLXHW8',FLXHW8  ,MAXSEG)
      CALL RINOUT ('W HWALKS',HWALKS  ,MAXNWS*2)
C
      CALL rINOUT ('W SOLTE ',solte,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('W SOLTI ',solti,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('W SOLNE ',solne,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('W SOLVEL',solvel,maxnks*msolpt+msolpt+1)
      CALL rINOUT ('W SOLCOR',solcor,maxnks*msolpt+msolpt+1)
C
c     Leakage data
c
      CALL RINOUT ('W CLEAKS',cleaks  ,MAXpts)
      CALL RINOUT ('W CLEAKN',cleakn  ,MAXpts*(maxizs+1))
      call rinout ('W LEAKPS',cleakpos,maximp*2)
c
c     More arrays related to leakage results
c
      call rinout ('W ncore ',ncore,maxnks*maxnrs)
      call rinout ('W nedge ',nedge,maxnks*maxnrs)
      call rinout ('W ntrap ',ntrap,maxnks*maxnrs)
      call rinout ('W ndivt ',ndivert,maxnks*maxnrs)
      call rinout ('W nmsol ',nmsol,maxnks*maxnrs)
c
c     NOTE: Not all of wtsource is being read or written - if 
c           this is ever needed both write and read routines
c           need to be adjusted.
c
c     jdemod - changed wtsou to *6 from *5 in 6.47
c
      call rinout ('W WTSOU ',wtsource,maxpts*maxnrs*4*6)
      call rinout ('W WTDEP ',wtdep,maxpts*(maxpts+1)*3)
      call rinout ('W TSOUR ',targsrc,3*4)
      call rinout ('W TLEAK ',targleak,3*4)
      call rinout ('W WSOUR ',wallsrc,5*3)
      call rinout ('W WLEAK ',wallleak,5*3)
c
c     Store any EDGE2D BG data that has been read in
c
      if (cre2d.eq.1.or.cre2d.eq.2.or.cre2d.eq.3.or.cre2d.eq.5) then
c
        call rinout ('W E2D N ',e2dnbs,maxnks*maxnrs)
        call rinout ('W E2D TE',e2dtebs,maxnks*maxnrs)
        call rinout ('W E2D TI',e2dtibs,maxnks*maxnrs)
        call rinout ('W E2D VB',e2dvhs,maxnks*maxnrs)
        call rinout ('W E2D E ',e2des,maxnks*maxnrs)
        call rinout ('W E2D I ',e2dion,maxnks*maxnrs)
        call rinout ('W E2D A ',e2datom,maxnks*maxnrs)
        call rinout ('W E2D TA',e2dtarg,maxnrs*8*2)
        call rinout ('W E2D GP',e2dgpara,maxnks*maxnrs)
        call rinout ('W E2D GD',e2dgdown,maxnks*maxnrs)
        call rinout ('W E2D G ',e2dflux,(maxnks+1)*maxnrs)
        call rinout ('W E2D VE',e2dbvel,(maxnks+1)*maxnrs)
c
        call rinout ('W E2D Z0',e2dz0,maxnks*maxnrs)
c
        call rinout ('W E2D RC',e2dhrec,maxnks*maxnrs)
        call rinout ('W E2D RC',e2drec,maxnks*maxnrs)
        call rinout ('W E2D CX',e2dcxrec,maxnks*maxnrs)
c
      if (.true.) then
         write(6,*) 'E2DNZS2:',nrs,nfla-1,cre2dizs
         do ir = 1,nrs
            do ik = 1,nks(ir)
               write(6,'(2i8,30(1x,g12.5))') ik,ir,
     >               (e2dnzs(ik,ir,iz),iz=1,nfla-1)
            end do
            write(6,*) '----------------'
         end do
      endif
c
      if (cre2dizs.gt.0) then
c
          call rinout ('W E2D NZ ',e2dnzs,maxnks*maxnrs
     >                 *(maxe2dizs+1)) 
c
c         jdemod - add e2dvzs to output
c
          call rinout ('W E2D VZ ',e2dvzs,maxnks*maxnrs
     >                 *(maxe2dizs+1))
          call rinout ('W E2D PW',e2dpowls,maxnks*maxnrs
     >                   *(maxe2dizs+1))
          call rinout ('W E2D LI',e2dlines,maxnks*maxnrs
     >                   *(maxe2dizs+1))

        endif
c
      endif
c
      call pr_trace('STORE','PART 4')
C
c     Store Data related to transport coefficient calculations
c
      call rinout ('W DPERP ',DPERP,maxnrs)
      call rinout ('W DPERPO',odperp,maxnrs)
      call rinout ('W DPERPI',idperp,maxnrs)
      call rinout ('W XPERP ',xPERPt,maxnrs)
      call rinout ('W XPERPO',oxperpt,maxnrs)
      call rinout ('W XPERPI',ixperpt,maxnrs)
      call rinout ('W XPI   ',chiperpi,maxnrs)
      call rinout ('W XPI  O',ochiperpi,maxnrs)
      call rinout ('W XPI  I',ichiperpi,maxnrs)
      call rinout ('W XPE   ',chiperpe,maxnrs)
      call rinout ('W XPE  O',ochiperpe,maxnrs)
      call rinout ('W XPE  I',ichiperpe,maxnrs)
c
c     GEometric quantities originally associated with transport extractor
c
      call rinout ('W RC OUT',rcouter,maxnrs)
      call rinout ('W RC IN ',rcinner,maxnrs)
      call rinout ('W ZC OUT',zcouter,maxnrs)
      call rinout ('W ZC IN ',zcinner,maxnrs)
      call rinout ('W MIDIST',middist,maxnrs*2)
c
      call rinout ('W GRADNE',gradn,maxnks*maxnrs)
      call rinout ('W GRADTE',gradte,maxnks*maxnrs)
      call rinout ('W GRADTI',gradti,maxnks*maxnrs)
      call rinout ('W E2DGNE',e2dgradn,maxnks*maxnrs)
      call rinout ('W E2DGTE',e2dgradte,maxnks*maxnrs)
      call rinout ('W E2DGTI',e2dgradti,maxnks*maxnrs)
c
c     IF Line profile data was calculated - store it 
c
      write(8) line_profile_opt
c
      if (line_profile_opt.ne.0) then 
          write(8) lp_wave,lp_instrument_width,
     >             lp_bin_width,lp_robs,lp_zobs,lp_theta,lp_dtheta
c slmod begin
c         Descriptor needs to be 8 characters long or generates
c         a runtime error in R8INOUT. -SL, 07/10/2011
          CALL R8INOUT('W LP    ',line_profile,max_lp_bins*2+1)
c
c          CALL R8INOUT ('W LP',line_profile,max_lp_bins*2+1)
c slmod end
          CALL R8INOUT('W MOD_LP',modified_line_profile,max_lp_bins*2+1)
      endif     
c
c     Store the pressure - from SOL option 22
c
      call rinout ('W KPRESS',kpress,maxnks*maxnrs*2)
c     IPP/11 - first arg must be 8 chars string
      call rinout ('W KPRAD ',kprad,maxnks*maxnrs)
c
c     Write out the global HC activation option
c
      write(8) global_hc_follow_option
c
c ammod begin
c
      call pr_trace('STORE','PART 5')
c
      if (global_hc_follow_option.ne.0) then  
c
c Write out HC module related quantities 
c
         call global_hc_store_raw_data
c
      endif
c
c     Added to raw file primarily for use by HC - code
c 
      call rinout ('W BRATIO',BRATIO,maxnks*maxnrs)

c
c ammod end
c
c
c     Write out subgrid data - at the least the option value
c     will be added to the raw file. If the subgrid was in use - 
c     its data will be saved and any storage used deallocated.  
c
      call save_subgrid(8)

c
c     jdemod - March 2016 - version 45
c
c     Write out potential and drift related results
c     
      call rinout ('W POT',osmpot2,(maxnks+2)*maxnrs)
      call rinout ('W E_RAD',e_rad,maxnks*maxnrs)
      call rinout ('W E_POL',e_pol,maxnks*maxnrs)
      call rinout ('W ExB_R',exb_rad_drft,maxnks*maxnrs)
      call rinout ('W ExB_P',exb_pol_drft,maxnks*maxnrs)
c
c     jdemod - version 48 
c
c     Add writing of far periphery related quantities
c     Only written if the option is active
c
      call fp_write_raw(8)
c
c     Temporarily Add the following (?) 
c
      call rinout ('W FLUXES',fluxes,maxnks*maxnrs*16)
c
      IF (IMODE.EQ.1) THEN
      CALL RINOUT ('W LIMS  ',LIMS  ,MAXNKS*MAXNRS*(MAXIZS+2)*MAXNTS)
      ENDIF


      call pr_trace('STORE','PART 6')
c
c
c slmod begin - new
      slver = 3.6

      WRITE(8) slver
      WRITE(8) MAXASD,MAXNAS,
     .         eirnpgdat,((eirpgdat(i1,i2),i2=1,MAXASD),i1=1,MAXNAS)
      WRITE(8) asc_ncell,MAXASC
      CALL IINOUT('W CELL  ',asc_cell  ,MAXASC)
      CALL IINOUT('W REGION',asc_region,MAXASC)
      CALL IINOUT('W LINK  ',asc_link  ,MAXASC*4)
      CALL IINOUT('W GRID  ',asc_grid  ,MAXASC*2)
      CALL IINOUT('W NVP   ',asc_nvp   ,MAXASC)
      CALL RINOUT('W RVP   ',asc_rvp   ,MAXASC*8)
      CALL RINOUT('W ZVP   ',asc_zvp   ,MAXASC*8)
c...  slver 2.1:
      CALL IINOUT('W NVERTX',ascnvertex,MAXASC)
c.... slver 2.2, switched to 40 from 20:
      CALL RINOUT('W VERTEX',ascvertex ,MAXASC*40)       

      CALL RINOUT ('W PINLN1',pinline(1,1,1,H_BALPHA),MAXNKS*MAXNRS*6)
c...  6.34:
      CALL RINOUT ('W PINLN2',pinline(1,1,1,H_BBETA ),MAXNKS*MAXNRS*6)
      CALL RINOUT ('W PINLN3',pinline(1,1,1,H_BGAMMA),MAXNKS*MAXNRS*6)
      WRITE(8) pincode
      CALL RINOUT('W PINMOI',pinmoi,MAXNKS*MAXNRS)
      CALL RINOUT('W OSMCDE',osmcde,MAXNKS*MAXNRS)
      CALL RINOUT('W OSMCDI',osmcdi,MAXNKS*MAXNRS)      
      CALL RINOUT('W OSMCVE',osmcve,MAXNKS*MAXNRS)
      CALL RINOUT('W OSMCVI',osmcvi,MAXNKS*MAXNRS)      
c...  slver 1.6:
      WRITE(8) tarshift(IKLO),tarshift(IKHI)
      WRITE(8) te_mult_o,ti_mult_o,n_mult_o,
     .         te_mult_i,ti_mult_i,n_mult_i
c...  slver 1.7:
      WRITE(8) lpdatsw
      CALL RINOUT('W SEPDIS',sepdist ,MAXNDS)
      CALL RINOUT('W SEPDI2',sepdist2,MAXNDS)
c...  slver 1.8:
      CALL IINOUT('W PRBNUM',prb_num ,NUMPRB) 
      CALL RINOUT('W PRBRHO',prb_rho ,MAXPRB*NUMPRB)
      CALL RINOUT('W PRBTE ',prb_te  ,MAXPRB*NUMPRB)
      CALL RINOUT('W PRBTI ',prb_ti  ,MAXPRB*NUMPRB)
      CALL RINOUT('W PRBNE ',prb_ne  ,MAXPRB*NUMPRB)
      CALL RINOUT('W PRBR  ',prb_r   ,MAXPRB*NUMPRB)
      CALL RINOUT('W PRBZ  ',prb_z   ,MAXPRB*NUMPRB)
c...  slver 1.9:
      WRITE(8) eirnres
      CALL RINOUT('W EIRRES',eirres,6*7*MAXPINITER)
c...  slver 2.4 (removed for 2.17):
      CALL RINOUT('W PINPLO',pinploss ,MAXNKS*MAXNRS*NMOMCHA)
c...  slver 2.5 & 3.4::
      WRITE(8) eirnpuff,eirpmode
      CALL RINOUT('W PUFF  ',eirpuff  ,MAXNAS*MAXPUFF)
c...  slver 2.7:
      CALL RINOUT('W BGK   ',pinbgk   ,MAXNKS*MAXNRS*MAXBGK)
      WRITE(8) MAXASD,MAXNAS3,
     .         eirnspdat,((eirspdat(i1,i2),i2=1,MAXASD),i1=1,MAXNAS3)
c...  slver 3.0:
      WRITE(8) asc_3dmode
      CALL IINOUT('W NVP3D ',asc_nvp3D ,MAXASC3D)
      CALL IINOUT('W LIMK3D',asc_link3D,MAXASC3D*6)
      CALL RINOUT('W ZMIN3D',asc_zmin3D,MAXASC3D)
      CALL RINOUT('W ZMAX3D',asc_zmax3D,MAXASC3D)
      CALL RINOUT('W XVP3D ',asc_xvp3D ,MAXASC3D*6*8)
      CALL RINOUT('W YVP3D ',asc_yvp3D ,MAXASC3D*6*8)
      CALL RINOUT('W ZVP3D ',asc_zvp3D ,MAXASC3D*6*8)
c...  slver 3.1:
      WRITE(8) ascncut,MAXASC3D,MAXASCDAT
c...  slver 3.2:
c     Resetting the 3D data reading in OUT.
c...  slver 3.3:
      WRITE(8) eirnsdtor,(eirsdtor(i1),eirsdvol(i1),i1=1,eirnsdtor)
      DO i1 = 2, eirnsdtor
        i2 = (i1-1)*MAXBGK+1
        CALL RINOUT('W BGKTOR',pinbgk(1,1,i2),MAXNKS*MAXNRS*MAXBGK)
      ENDDO
c...  6.13:
      WRITE(8) eirniontime,MAXIONTIME,MAXBIN
      CALL RINOUT('W IONTIM',eiriontime,MAXIONTIME*(20+MAXBIN*3))
c...  6.14:
      WRITE(8) MAXASD2,MAXASS,asc_ncell,ascncut
      WRITE(8) 
     .  (asc_vol(i1),
     .   (ascdata(i1,i2),i2=1,5),
     .   ((pinasd(i1,i2,i3,1),pinasd(i1,i2,i3,2),i2=1,MAXASD2),
     .                                           i3=1,MAXASS ), 
     .  i1=1,asc_ncell*ascncut+1+eirnpgdat),999999
c...  6.15:
      WRITE(8) eirntally,MAXNKS,MAXNRS,MAXTALLY
      WRITE(8) 
     .  ((eirtally(i1,i2),i2=1,4),
     .   ((pinalgv(ik,ir,i1),ik=1,nks(ir)),ir=1,nrs),
     .   i1=1,eirntally),999999
c...  6.16:
      WRITE(8) eirzaa,eirtorfrac,eirsrcmul
c...  6.22:
      WRITE(8) osmns28,8,osm_nfnt,3,MAXFNT
      WRITE(8) ((osms28(i1,i2),i2=1,8),i1=1,osmns28),
     .         ((osm_ionfnt(i1,i2),i2=1,3),i1=1,osm_nfnt)
c...  6.25:
      WRITE(8) eirnstrdat,eirnstrai,eirnatmi,eirnmoli,
     .         MAXSTRDAT,MAXSTRATA,100
      WRITE(8) (((eirstrdat(i1,i2,i3),i3=1,100      ),
     .                                i2=1,MAXSTRATA),
     .          eirstrlis(i1)        ,i1=1,MAXSTRDAT),999999
c...  6.26:
      CALL IINOUT('W IKBOUN',ikbound,MAXNRS*2)
      CALL IINOUT('W IKBOUN',ikfluid,MAXNRS*2)
c...  6.30:
c     This is temporary.  I want the ability to store ionisation data
c     for a greater number of strata, so I made room in the PINDATA array
c     by not storing the atom and molecule density information for each
c     stratum - SL, Sep 19, 2002:
      DO i1 = H_ION1, H_ION1+11
        CALL IINOUT('W PINDAT',pindata(1,1,i1),MAXNKS*MAXNRS)
      ENDDO
c...  6.28:
      DO i1 = 1, 3
        CALL IINOUT('W PINICP',pinioncomp(1,1,i1),MAXNKS*MAXNRS)
      ENDDO
c...  6.29:
      WRITE(8) eirntorseg
c...  6.33:
      WRITE(8) ciopte,cxsc,cysc,cxsca,cysca,cxscb,cyscb

c...  slver = 3.5: *TEMP*
      CALL RINOUT('W EIRPH1',eirpho1,MAXNKS*MAXNRS)
      CALL RINOUT('W EIRPH2',eirpho2,MAXNKS*MAXNRS)
 
c...  6.41:
      WRITE(8) debugv,cstepv
c
c     jdemod - save this all the time
c      
c      IF (debugv) CALL RINOUT ('W SDVS',sdvs,MAXNKS*MAXNRS*(MAXIZS+2))
      CALL DINOUT ('W DDVS',ddvs,MAXNKS*MAXNRS*(MAXIZS+2))

c...  slver 3.6:      
      IF (ALLOCATED(wall_flx)) THEN 
        WRITE(8) 1
        WRITE(8) wall_n,1.0  ! this 1.0 is a version number
        WRITE(8) MAXNBLK,MAXNATM,MAXNMOL,MAXNSRC,MAXNLAUNCH
        WRITE(8) wall_flx
      ELSE
        WRITE(8) 0
      ENDIF

c...  6.14 (end of file flag):
      WRITE(8) 123456789


      call pr_trace('STORE','AFTER RAW FILE WRITTEN')


      if (netcdf_opt.eq.1) then 
         call write_netcdf_output(TITLE,desc,NIZS,JOB,EQUIL,
     >                  FACTA,FACTB,ITER,NITERS)
      endif

      call pr_trace('STORE','AFTER NETCDF WRITE')


c slmod end
c
c     LEAVE RETURN AT END OF ROUTINE  
c
      RETURN
c
 9001 FORMAT(//1X,'GET:     NXS    NYS    NRS    NDS   NIZS    NTS',
     >  /6X,6I7,/1X,'(MAX)',6I7,/6X,A,/6X,A,
     >        /1X,'        ITER   MODE   REFCT  MAXSEG NVESM  NVESP',
     >  /6x,6I7,/)
c 9001 FORMAT(//1X,'STORE:   NXS    NYS    NRS    NDS   NIZS    NTS',
c     >  /6X,6I7,/1X,'(MAX)',6I7,/6X,A,/6X,A,/1X,'ITER',I4,/1X,
c     >  'MODE',I4,'REFCT',I4,/)
      END
C
C
C
      subroutine check_ddlim(nizs,num)
      use mod_params
      use mod_dynam1
      use mod_cgeom
      implicit none
c
c     CHECK_DDLIM: Verifies that all elements of DDLIM contain 
c                  values that are greater than or equal to zero.
c                  A negative value would indicate a coding error. 
c
c
      integer nizs,num
c     include 'params'
c     include 'dynam1'
c     include 'cgeom'
c
      integer ik,ir,iz
c
      do iz = 0,nizs
         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (ddlims(ik,ir,iz).lt.0.0) then
                  write (6,'(a,4i4,g16.8)') 'DDLIM error:',
     >                     num,ik,ir,iz,ddlims(ik,ir,iz)
               endif
            end do
         end do
      end do

      return
      end
