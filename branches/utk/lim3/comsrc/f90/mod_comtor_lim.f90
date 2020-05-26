module mod_comtor

  use mod_params

  implicit none


  private

!c -*-Fortran-*-
!C                                                                               
!c     jdemod - cprint moved to global_options
!c
!      COMMON /COMTOR/ CA,    CAW,   CL,    CRMB,  CTBIN, CGTIN1,CLTIN1,         
!     >                CTBOUL,CLTOUL,CNBIN, CGNIN1,CLNIN1,CNBOUL,CLNOUL,         
!     >                CVIN,  CRDXO, CRDD,  CRMI,  CXSC,  CFIMP, CQS,            
!     >                CYSC,  CTEMSC,CTIMSC,CEBD,  CTGAS, CEMAXF,CENGSC,         
!     >                CIZB,  CION,  CIZSC, CATIN, CANIN, CTRESH,CIZEFF,         
!     >                CNHC,  CNHO,  CLAMHX,CLAMHY,CXNEAR,CYNEAR,CORECT,         
!     >                CTBINS,CNBINS,CCUT,  CPFIR, CPSUB, CPSC,  CSNORM,         
!     >    CDPOL, CPLSMA,CEYIN,CVHYIN,CSTEPN,CIZSET,CZENH,debugt,ptracl,          
!     >    CEYOUT,CVHOUT,CISEED,CDIFOP,CTWOL, CYSTAG,CBOMBF,PTRACS,         
!     >    CBOMBZ,CIRF,  CSEF,  CGTIN2,CGNIN2,CLTIN2,CLNIN2,CSTEPT,         
!     >    CEIN2, CSOLEF,CVCX,CYMFS, CTHETB,CSINTB,CVPOUT,CVXMIN,         
!     >    CXSPLS,CNSPL, CTBOUG,CLTOUG,CNBOUG,CONI  ,CPRUL,CVXMAX,          
!     >    CLNOUG,CLARMR,CKI   ,CTBI,  CHALFL,DEBUGN,CKO, SVYBAR,SVYACC,        
!     >    CONO,CYMFLG,CFTCUT,CNBA  ,CGAMMA,CANAL ,CAP, VPV0,VPALPH,         
!     >    CWL,CSTEPL,DEBUGL,CLFACT,CYFAR ,CRDXI ,CTSUB, CVPCUT,DPBETA,        
!     >    CQPL,CQSL, CIOPTJ, CPCO,CTIBIN,CLTIIN1,CGTIIN1,CATIIN,CDPERP,
!     >    DPALPH,
!     >    CTIBOUL,CLTIOUL,CTIBOUG,CLTIOUG,CGTIIN2,CLTIIN2,CVPOPT,
!     >    CTIBINS,CPROB,CTICHG,CLPD,LPDION,LPDCUM,CVSA,RLC,CDPSTP,
!     >    CBRK,RLEDGE7,CALPHE,CBETAI,QMULTP,QMULTS,CIANGN,
!     >    YCFADD,CSPUMAX,CFBGFF,C3HALFL,CSVYMIN,CMAXGENS,CDCALC,CVPOL
!C
!      REAL            CVPOUT,CVXMIN,CVXMAX
!      REAL            SVYBAR(-MAXQXS:MAXQXS),SVYACC(-MAXQXS:MAXQXS)
!      REAL            CA,CAW,CL,CRMB,CTBIN,CGTIN1,CLTIN1,CTBOUL,CLTOUL          
!      REAL            CNBIN,CGNIN1,CLNIN1,CNBOUL,CLNOUL,CFIMP,CRDD,CDPOL        
!      REAL            CVIN,CRMI,CXSC,CYSC,CTEMSC,CTIMSC,CEBD,CTGAS              
!      REAL            CEMAXF,CENGSC,CATIN,CANIN,CTRESH,CQS(MAXINS,2)            
!      REAL            CNHC,CNHO,CLAMHX,CLAMHY,CXNEAR,CYNEAR,CCUT,CSTEPN         
!      REAL            CPFIR,CPSUB,CPSC,CSNORM,CPLSMA,CEYIN,CVHYIN,CZENH         
!      REAL            CEYOUT,CVHOUT,CTWOL,CYSTAG,CIRF,CSEF,CSOLEF,CVCX          
!      REAL            CGTIN2,CIANGN,CGNIN2,CLTIN2,CLNIN2,CEIN2,CTHETB,
!     >                CSINTB           
!      REAL            CTBINS(MAXINS,2),CNBINS(MAXINS,2),CYMFS(MAXINS,3)         
!      REAL            CXSPLS(0:MAXINS+1),CTBOUG,CLTOUG,CNBOUG,CLNOUG            
!      REAL            CLARMR,CKI,CKO,CTBI,CPRUL,CHALFL,CONI,CONO,CFTCUT         
!      REAL            CNBA,CGAMMA,CANAL,CAP,CWL,CSTEPL,CLFACT,CYFAR             
!      REAL            CRDXO(3),CRDXI(3),CTSUB,CQPL,CQSL,CPCO,CPROB
!      REAL            CTIBIN,CLTIIN1,CGTIIN1,CATIIN,CTIBOUL,CLTIOUL
!      REAL            CTIBOUG,CLTIOUG,CGTIIN2,CLTIIN2,CTIBINS(MAXINS,2)
!      REAL            LPDION(MAXLPD,2),LPDCUM(MAXLPD)
!      REAL            CVSA(MAXINS,3),RLEDGE7
!      REAL            CALPHE(MAXIZS),CBETAI(MAXIZS)
!      REAL            YCFADD,CSPUMAX,CFBGFF,C3HALFL,CSVYMIN
!      REAL            QMULTP,QMULTS  
!      REAL            CVPOL,RLC, CDPSTP
!      REAL            DPALPH,DPBETA,VPV0,VPALPH,CVPCUT 
!      REAL            PTRACS(maxlen,maxt,2)
!c
!      integer         cioptj
!      INTEGER         CMAXGENS,CDCALC,CDPERP,CVPOPT
!      INTEGER         CIZB,CION,CIZSC,CIZEFF,CIZSET,CISEED,CDIFOP,CORECT        
!      INTEGER         CBOMBF,CBOMBZ,CNSPL,CYMFLG,CLPD,CBRK,CSTEPT
!      INTEGER         PTRACL(MAXT) 
!      LOGICAL         DEBUGN,DEBUGL,debugt,CTICHG                        
!c
!c     This file contains declarations for the unstructured 
!c     input values. If it becomes too unwieldy this file will 
!c     be split into separate common blocks for different 
!c     unstructured input values. 
!c
!      common /sputdat/ csputopt,extra_sputter_angle,
!     >                 cchemopt, const_yield, impact_energy_opt,
!     >                 init_y_coord,cselfs,
!     >                 ss_nymfs,ss_cymfs
!      integer csputopt,cchemopt,impact_energy_opt,cselfs
!      integer ss_nymfs
!      real    extra_sputter_angle,init_y_coord
!      real    const_yield        
!      real    ss_cymfs(maxins,3) 
!c
!c     Gradient multipliers
!c
!      common /gradmult/ NTIG,NTEG,nnbg,TMIG,TMEG,mnbg
!      integer NTIG,NTEG,nnbg
!      REAL    TMIG(MAXINS,2),TMEG(MAXINS,2),MNBG(maxins,2)
!c
!c     Scaling factor
!c     
!      common /scaling_factor/ absfac
!c
!c     In DIVIMP absfac is a real - should LIM be the same? Since all use of absfac should derive 
!c     from this common block include - the declaration should remain consistent within LIM - however
!c     using real*8 here runs the risk that imported code might have an issue. 
!c
!      real*8 absfac
!c      real absfac
!c
!c     Common block for generic LIM related unstructured input
!c      
!      common /lim_unstruc/ shear_short_circuit_opt,calc_3d_power,
!     >                     extfluxopt,nextfluxdata,extfluxdata,
!     >                     vpflow_3d
!c
!      integer shear_short_circuit_opt,calc_3d_power,extfluxopt,
!     >        nextfluxdata
!      real extfluxdata(maxins,3),vpflow_3d
!c
!c     Common block for OUT related unstructured input
!c
!      common /out_unstruc/ new_absfac,erosion_scaling_opt
!c
!      real*8 new_absfac
!c
!c      real new_absfac
!c
!      integer erosion_scaling_opt
!


  
!      COMMON /COMTOR/ CA,    CAW,   CL,    CRMB,  CTBIN, CGTIN1,CLTIN1,         
!     >                CTBOUL,CLTOUL,CNBIN, CGNIN1,CLNIN1,CNBOUL,CLNOUL,         
!     >                CVIN,  CRDXO, CRDD,  CRMI,  CXSC,  CFIMP, CQS,            
!     >                CYSC,  CTEMSC,CTIMSC,CEBD,  CTGAS, CEMAXF,CENGSC,         
!     >                CIZB,  CION,  CIZSC, CATIN, CANIN, CTRESH,CIZEFF,         
!     >                CNHC,  CNHO,  CLAMHX,CLAMHY,CXNEAR,CYNEAR,CORECT,         
!     >                CTBINS,CNBINS,CCUT,  CPFIR, CPSUB, CPSC,  CSNORM,         
!     >    CDPOL, CPLSMA,CEYIN,CVHYIN,CSTEPN,CIZSET,CZENH,debugt,ptracl,          
!     >    CEYOUT,CVHOUT,CISEED,CDIFOP,CTWOL, CYSTAG,CBOMBF,PTRACS,         
!     >    CBOMBZ,CIRF,  CSEF,  CGTIN2,CGNIN2,CLTIN2,CLNIN2,CSTEPT,         
!     >    CEIN2, CSOLEF,CVCX,CYMFS, CTHETB,CSINTB,CVPOUT,CVXMIN,         
!     >    CXSPLS,CNSPL, CTBOUG,CLTOUG,CNBOUG,CONI  ,CPRUL,CVXMAX,          
!     >    CLNOUG,CLARMR,CKI   ,CTBI,  CHALFL,DEBUGN,CKO, SVYBAR,SVYACC,        
!     >    CONO,CYMFLG,CFTCUT,CNBA  ,CGAMMA,CANAL ,CAP, VPV0,VPALPH,         
!     >    CWL,CSTEPL,DEBUGL,CLFACT,CYFAR ,CRDXI ,CTSUB, CVPCUT,DPBETA,        
!     >    CQPL,CQSL, CIOPTJ, CPCO,CTIBIN,CLTIIN1,CGTIIN1,CATIIN,CDPERP,
!     >    DPALPH,
!     >    CTIBOUL,CLTIOUL,CTIBOUG,CLTIOUG,CGTIIN2,CLTIIN2,CVPOPT,
!     >    CTIBINS,CPROB,CTICHG,CLPD,LPDION,LPDCUM,CVSA,RLC,CDPSTP,
!     >    CBRK,RLEDGE7,CALPHE,CBETAI,QMULTP,QMULTS,CIANGN,
!     >    YCFADD,CSPUMAX,CFBGFF,C3HALFL,CSVYMIN,CMAXGENS,CDCALC,CVPOL
!C
      REAL,public ::  CVPOUT,CVXMIN,CVXMAX,&
                  SVYBAR(-MAXQXS:MAXQXS),SVYACC(-MAXQXS:MAXQXS),&
                  CA,CAW,CL,CRMB,CTBIN,CGTIN1,CLTIN1,CTBOUL,CLTOUL,&          
                  CNBIN,CGNIN1,CLNIN1,CNBOUL,CLNOUL,CFIMP,CRDD,CDPOL,&        
                  CVIN,CRMI,CXSC,CYSC,CTEMSC,CTIMSC,CEBD,CTGAS,&              
                  CEMAXF,CENGSC,CATIN,CANIN,CTRESH,CQS(MAXINS,2),&            
                  CNHC,CNHO,CLAMHX,CLAMHY,CXNEAR,CYNEAR,CCUT,CSTEPN,&         
                  CPFIR,CPSUB,CPSC,CSNORM,CPLSMA,CEYIN,CVHYIN,CZENH,&         
                  CEYOUT,CVHOUT,CTWOL,CYSTAG,CIRF,CSEF,CSOLEF,CVCX,&          
                  CGTIN2,CIANGN,CGNIN2,CLTIN2,CLNIN2,CEIN2,CTHETB,CSINTB,&           
                  CTBINS(MAXINS,2),CNBINS(MAXINS,2),CYMFS(MAXINS,3),&         
                  CXSPLS(0:MAXINS+1),CTBOUG,CLTOUG,CNBOUG,CLNOUG,&            
                  CLARMR,CKI,CKO,CTBI,CPRUL,CHALFL,CONI,CONO,CFTCUT,&         
                  CNBA,CGAMMA,CANAL,CAP,CWL,CSTEPL,CLFACT,CYFAR,&             
                  CRDXO(3),CRDXI(3),CTSUB,CQPL,CQSL,CPCO,CPROB,&
                  CTIBIN,CLTIIN1,CGTIIN1,CATIIN,CTIBOUL,CLTIOUL,&
                  CTIBOUG,CLTIOUG,CGTIIN2,CLTIIN2,CTIBINS(MAXINS,2),&
                  LPDION(MAXLPD,2),LPDCUM(MAXLPD),&
                  CVSA(MAXINS,3),RLEDGE7,&
                  CALPHE(MAXIZS),CBETAI(MAXIZS),&
                  YCFADD,CSPUMAX,CFBGFF,C3HALFL,CSVYMIN,&
                  QMULTP,QMULTS,&  
                  CVPOL,RLC, CDPSTP,&
                  DPALPH,DPBETA,VPV0,VPALPH,CVPCUT,& 
                  PTRACS(maxlen,maxt,2)


      real,public :: te_prof_shift,ti_prof_shift,ne_prof_shift
      real,public :: te_prof_mult,ti_prof_mult,ne_prof_mult

      
      integer,public:: cioptj,&
               CMAXGENS,CDCALC,CDPERP,CVPOPT,&
               CIZB,CION,CIZSC,CIZEFF,CIZSET,CISEED,CDIFOP,CORECT,&        
               CBOMBF,CBOMBZ,CNSPL,CYMFLG,CLPD,CBRK,CSTEPT,&
               PTRACL(MAXT) 

      LOGICAL,public::  DEBUGN,DEBUGL,debugt,CTICHG                        
!c
!c     This file contains declarations for the unstructured 
!c     input values. If it becomes too unwieldy this file will 
!c     be split into separate common blocks for different 
!c     unstructured input values. 
!c
!      common /sputdat/ csputopt,extra_sputter_angle,
!     >                 cchemopt, const_yield, impact_energy_opt,
!     >                 init_y_coord,cselfs,
!     >                 ss_nymfs,ss_cymfs

      integer,public:: csputopt,cchemopt,impact_energy_opt,cselfs,ss_nymfs
      real,public::    extra_sputter_angle,init_y_coord,const_yield,ss_cymfs(maxins,3) 
!c
!c     Gradient multipliers
!c
!      common /gradmult/ NTIG,NTEG,nnbg,TMIG,TMEG,mnbg
      integer,public ::  NTIG,NTEG,nnbg
      REAL,public::    TMIG(MAXINS,2),TMEG(MAXINS,2),MNBG(maxins,2)
!c
!c     Scaling factor
!c     
!      common /scaling_factor/ absfac
!c
!c     In DIVIMP absfac is a real - should LIM be the same? Since all use of absfac should derive 
!c     from this common block include - the declaration should remain consistent within LIM - however
!c     using real*8 here runs the risk that imported code might have an issue. 
!c
      real*8,public ::  absfac
!c      real absfac
!c
!c     Common block for generic LIM related unstructured input
!c      
!      common /lim_unstruc/ shear_short_circuit_opt,calc_3d_power,
!     >                     extfluxopt,nextfluxdata,extfluxdata,
!     >                     vpflow_3d
!c
      integer,public:: shear_short_circuit_opt,calc_3d_power,extfluxopt,nextfluxdata
      real,public:: extfluxdata(maxins,3),vpflow_3d
!c
!c     Common block for OUT related unstructured input
!c
!      common /out_unstruc/ new_absfac,erosion_scaling_opt
!c
      real*8,public:: new_absfac
!c
!c      real new_absfac
!c
      integer,public:: erosion_scaling_opt

!c    sazmod - switch to vary absorbing boundary that affects plasma solution
!c             and the accompanying values.
	  !integer, public:: vary_absorb, ix_step1, ix_step2
	  !real,    public:: yabsorb1a_step, yabsorb2a_step, xabsorb1a_step, xabsorb2a_step
	  
      ! Add in regions for varying radial Dperp values.
	  integer, public:: dperp_reg_switch
	  real,    public:: dperp_reg1, dperp_reg2, dperp_reg3, dperp_reg4
	  
	  ! If we don't care about the .raw file save some time and skip it.
	  integer, public:: skip_raw
	  
	  ! Modify the velplasma in the step region (right half only right now).
	  real,    public:: mod_v_fact
	  
	  ! Option to chose from an exponential distribution in the Y direction (3D only).
	  integer, public:: choose_exp
	  real,    public:: choose_exp_lambda, choose_exp_fact
	  
	  ! Overall scaling factor to apply to the background plasma velocity.
	  real,    public:: vel_mod



  public :: allocate_mod_comtor, deallocate_mod_comtor


contains

  subroutine allocate_mod_comtor
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(PTES,maxnxs,'PTES',ierr)

  end subroutine allocate_mod_comtor


  subroutine deallocate_mod_comtor
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_comtor



end module mod_comtor
