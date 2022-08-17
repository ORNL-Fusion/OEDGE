!     -*-Fortran-*-

!    sazmod - The code is very messy with tons of commented out code,
!      much of which is probably long-obsolete anyways. Some possibly
!      relevant commented code has been left in, but some of the older
!      commented out code has been removed. Other miscellenous 
!      aesthetic notes include: 
!      - As I run through code I make it lowercase since all caps seems 
!         to be an antiquidated and ugly convention. 
!      - Generally try to align the spacing, i.e. if and endif start in 
!         the same column.
!      - Comments swapped from "c" to "!" so f90 syntax highlighting 
!         works.
!      - Spaces added to let code breathe some, according to modern 
!         Fortran style conventions.
!      - Replacing old style do-continue loops with do-enddo loops.
!      1/18/21.

      subroutine lim3 (imode,nizs,nimps,impadd,impcf,qtim,cpulim,prinps,         
     >  xwidm,ywidm,fsrate,iontim,neutim,seed,igeom,ntbs,ntibs,nnbs,
     >  nymfs,ncvs,facta,factb,iter,defact,nrand,title)                     

!      use iter_bm
      use mod_params
      use debug_options
      use eckstein_2007_yield_data
      use variable_wall
      use yreflection
      use mod_dynam1
      use mod_dynam3
      use mod_comt2
      use mod_comnet
      use mod_cneut
      use mod_cnoco
      use mod_comtor
      use mod_cadas
      use mod_commv
      use mod_comtau
      use mod_comxyt
      use mod_coords
      use mod_zommv
      use mod_save
      use mod_crand
      use mod_printr
      use mod_global_options
      use mod_slcom
      use mod_soledge
      use mod_lim3_local
      use mod_diagvel
      
      implicit none                                                                                                            
                                                                         
!  *********************************************************************        
!  *                                                                   *        
!  *   LIM3: Main controlling routine                                  *        
!  *   ------------------------------                                  *
!  *   Shawn's update to this block of comment (12/14/21)              *
!  *                                                                   *
!  *   LIM, or now known as 3DLIM, has changed a lot since the code    *
!  *   was originally written, so there may be quite a few areas with  *
!  *   an outdated description of the way the code works. Such is the  *
!  *   the reality when working with a code that was written in the    *
!  *   late 80's (this comment written in late 2021, so the code is    *
!  *   over 30 years old at the time of this writing!). The code has   *
!  *   been revived numerous time for specific purposes, those         *
!  *   purposes guiding the development of the code. Examples include  *
!  *   making the the code 3D with a poloidal direction, additional    *
!  *   upgrades to look at the ITER blanket module, incorporation of a *
!  *   collector probe, and most recently allowing the Y-absorbing     *
!  *   bounds to have structure in the radial and poloidal dimensions. *
!  *   Thus, not all the options will work together, some options not  *
!  *   working at all anymore.                                         *
!  *   I expect to be the last developer of this code, as full-scale   *
!  *   kinetic codes are finally arriving on the scene, but until they *
!  *   are truly production-ready, 3DLIM remains an industry standard  *
!  *   of far-SOL impurity transport.                                  *
!  *   In homage to those who gave their lives to this code,           *
!  *   developers of LIM/3DLIM include:                                *
!  *    - Chris Farrell (Hunterskil)                                   *
!  *    - Peter Stangeby                                               *
!  *    - David Elder                                                  *
!  *    - Steve Lisgo                                                  *
!  *    - Shawn Zamperini                                              *
!  *                                                                   *        
!  *   This  routine follows the diffusion with time                   *        
!  *     of a set of injected impurity ions and returns                *        
!  *     either the state of the ion cloud at a set of time            *        
!  *     points (impulse mode), or the steady state distribution       *        
!  *    (steady state mode) or both.                                   *        
!  *   Each ion possesses (X,Y,P) coordinates, a Y direction           *        
!  *     velocity, a temperature and an ionisation state, all of which *        
!  *     change with time.                                             *        
!  *   The ions are followed until they are absorbed, ionise           *        
!  *     beyond a given state or until a cutoff time is reached.       *        
!  *   Note that the tau factors, the background temperature           *        
!  *     and density and the limiter edge are found for a set          *        
!  *     of X positions before the iterations begin.  during the       *        
!  *     iterations, the values nearest the current (X,Y,P) position   *        
!  *     are used.                                                     *        
!  *   Similarly the outboard electric field and drift velocity        *        
!  *     are calculated for a set of Y positions along with a          *        
!  *     scaling factor for each X position outboard so that           *        
!  *     during the iteration a value can be calculated for any        *        
!  *     position by taking the product of the applicable X and Y      *        
!  *     position values.                                              *        
!  *   This requires making the assumption that the electric field     *        
!  *     values are proportional to temb/l and that the background     *        
!  *     velocities are proportional to sqrt(temb).                    *        
!  *   Times are scaled by 1/qtim. this means that time values         *        
!  *     can be stored in integers. input time values are rounded to   *        
!  *     the nearest integer (ie the nearest timestep).                *        
!  *   Y direction velocity values are scaled by qtim.                 *        
!  *     This again saves an inner loop multiplication as the change   *        
!  *     in Y at each time step becomes y = y + vy.                    *        
!  *   LIM monitors various aspects of the diffusion and               *        
!  *     prints out a set of diagnostics at the end of the run.        *        
!  *                                                                   *        
!  *  Arguments :-                                                     *        
!  *  imode  : Set to 1 for impulse mode, 2 for steady state mode      *        
!  *               and 0 for both                                      *        
!  *  nizs   : Maximum ionization state to be followed                 *        
!  *  nimps  : Number of impurity ions to be used in monte carlo       *        
!  *  impadd : Number of additional neutrals to be launched to         *
!  *            simulate cross-field fluxes                            *
!  *  impcf  : Number of cross-field sputtered primary impurities      *
!  *            the number is calculated in neut if cfbgff.gt.0.0      *
!  *            otherwise it is zero (initialized at start of runlm3)  *
!  *  qtim   : Size of quantum timestep to be used following ions      *        
!  *  cpulim : Maximum amount of time to be used by routine (secs)     *        
!  *            a value of 0 indicates infinite time                   *        
!  *  xwidm  : Minimum X bin width, calculated from bin boundaries     *        
!  *  ywidm  : Minimum Y bin width, calculated from bin boundaries     *        
!  *  fsrate : Size of timestep to be used following neutrals          *        
!  *  iontim : CPU time spent following ions accumulator               *        
!  *  neutim : CPU time spent following neutrals accumulator           *        
!  *  seed   : Current random number seed value for surand (d.p)       *        
!  *  igeom  : Radial geometry flag:  0 slab,  1 cylinder              *        
!  *                                                                   *        
!  *  defact : Factor for converting ddlims to actual physical density *        
!  *                                                                   *        
!  *                        Chris Farrell (Hunterskil)  March 1988     *        
!  *                                                                   *        
!  *********************************************************************

!     *** Old local variable declarations used to be here ***                                                            
!     jdemod - NOTE: Most local variables here were moved to the module 
!       mod_lim_local to faciliate the conversion to dynamic storage 
!       allocation.

      integer :: imode, nizs, nimps, igeom, ntbs, ntibs, nnbs, nymfs 
      integer :: iter, nrand, ncvs, impadd, impcf, pz, pz1, pz2, ip2
      integer :: outunit
      real :: qtim, cpulim, xwidm, ywidm, fsrate, iontim, neutim  
      real :: velplasma_val, efield_val, cx_start, cdf_sum
      real :: sum_divimp_prob, max_divimp_s, tmp_width, max_width
      real :: facta(-1:maxizs), factb(-1:maxizs)
      character :: prinps(-maxnps-1:maxnps)*7, title*80, what(51)*10
      character :: fate(11)*16, string*21  
      double precision :: seed, defact 
      logical, external :: res
      real, external :: za02as, yield
      integer, external :: ipos, jpos

      ! Local temporary time variables
      real*8 :: time_frac,time_start,time_end,time_win,dtime

      data  fate  /'REACHED X=AW',        'HIT Y=0 FROM Y>0',                   
     >             'REACHED Y=2L',        'HIT Y=0 FROM Y<0',                   
     >             'REACHED Y=-2L',       'REACHED TIME CUT',                   
     >             'SPLITTING ION',       'ROULETTE DISCARD',                   
     >             'CHARGE CHECK',        'X-ABSORPTION',
     >             'Y-ABSORPTION'/                                              
                                                                               
      data  what  /'  PRIMARY ',   ' SECONDARY',   ' TERTIARY ',                
     >             'QUATERNARY',   '  QUINARY ',   '  SIXTH   ',                
     >             '  SEVENTH ',   '  EIGHTH  ',   '  NINTH   ',                
     >             '  TENTH   ',   ' ELEVENTH ',   ' TWELFTH  ',                
     >             'THIRTEENTH',   'FOURTEENTH',   ' FIFTEENTH',                
     >             ' SIXTEENTH',   ' SEVENTEEN',   'EIGHTEENTH',
     >             'NINETEENTH',   ' TWENTIETH',   'TWENTY-1ST',
     >             'TWENTY-2ND',   'TWENTY-3RD',   'TWENTY-4TH',
     >             'TWENTY-5TH',   'TWENTY-6TH',   'TWENTY-7TH',
     >             'TWENTY-8TH',   'TWENTY-9TH',   ' THIRTIETH',
     >             'THIRTY-1ST',   'THIRTY-2ND',   'THIRTY-3RD',
     >             'THIRTY-4TH',   'THIRTY-5TH',   'THIRTY-6TH',
     >             'THIRTY-7TH',   'THIRTY-8TH',   'THIRTY-9TH',
     >             ' FORTIETH ',   ' FORTY-1ST',   ' FORTY-2ND',
     >             ' FORTY-3RD',   ' FORTY-4TH',   ' FORTY-5TH',
     >             ' FORTY-6TH',   ' FORTY-7TH',   ' FORTY-8TH',
     >             ' FORTY-9TH',   ' FIFTIETH ',   
     >             '    ALL   '/                                                
                                                                              
!-----------------------------------------------------------------------        
!                   Initialisation                                              
!-----------------------------------------------------------------------        

      write(0,*) 'Begin LIM3'

      if (optdp.eq.1) then
        write(0,*) 'Warning! Hard code adjustment to bin location',
     >             ' for DIVIMP ion profile.'
      endif 

      ! Comment needed.
      do ii = 1, nbin
        bsbin(ii) = bsbin(ii) + 0.5
        ysbin(ii) = ysbin(ii) + 0.5
      enddo

      ! Comment needed to describe each variable.
      avgtrac = 0.0
      tgloss  = 0.0
      wloss   = 0.0
      lloss   = 0.0
      izloss  = 0.0
      tsloss  = 0.0
      aloss   = 0.0
      mark    = 0.005
      target  =-5.0

      ! Comment needed.
      do iy=-nys-1,nys+1
        injbint(iy) = 0
      enddo

      ! Comment needed.
      tavxpos = 0.0
      if (cpulim.le.0.0) cpulim = 1.0e7                                         
      dwol = dble (ctwol)                                                       
      if (cdpstp.ne.0.0) then
         dpprob = 2.0 * qtim / cdpstp /cdpstp
      else
         dpprob = 0.0
      endif
                                                                              
      ! Plasma elongation - set-up delps array and cono, coni constants.                                                                               
      coni = (cki - 1.0) / chalfl                                                 
      cono = (cko - 1.0) / chalfl                                                 
      do iy = 1, nys                                                        
        yy = mod (youts(iy), cl)                                                
        if (yy.gt.chalfl) yy = cl - yy                                           
        do ix = 1, nxs                                                       
          if (xs(ix).gt.0.0) then                                               
            delps(ix,iy) = 1.0 / (yy * coni + 1.0)                                
          else                                                                  
            delps(ix,iy) = 1.0 / (yy * cono + 1.0)                                
          endif                                                                 
        end do                                                               
      end do                                                                  
   
      ! jdemod - Moved these calculations to before TAU is called so 
      ! that inboard flow and efield defaults are available in TAU. 
      do iqx = 1-nqxso, nqxsi                                               
        polods(iqx) = sqrt (2.0 * cdpol * qtim * qs(iqx))                       
        svpols(iqx) = qtim * qs(iqx) * cvpol
        svhins(iqx) = qtim * qs(iqx) * cvhyin                                   
        do iz = 1, nizs                                                      
          seyins(iqx,iz) = real (iz) * qtim * qs(iqx) * qtim * qs(iqx) *        
     >      (1.602192e-19 / 1.672614e-27) * ceyin / crmi           
        end do                                                                
      end do                                                                 
    
      ! Set up factors in common comtau. Further comment needed on the 
      ! role of each variable.                                                                                                                        
      if (iter.eq.1) call tauin1 (qtim, nizs, icut, fsrate, igeom, 
     >  ntbs, ntibs, nnbs)                           
                                                                               
      ! Set up vfluid, the fluid velocity.                                                                                                                       
      if(abs(cvhys(1)).eq.0.0) then                                             
         vfluid = 1.56e4 * sqrt (ctbin / crmb)                                    
      else                                                                      
         vfluid = abs (cvhys(1))                                                
      endif                                                                     
      write (6,*) 'LIM3: VFLUID=', vfluid                                       
                                                                              
      ! Set up fvycol: Post collision velocity factor    *** not used             
      ! so no need for an array of iqx vals !          
      ! vy0   : Initial velocity for non-neut cases                        
      ! vy02  : Second velocity for injection 4
      ! polods: Poloidal diffusion factor                                  
      ! svpols: Scaled poloidal drift velocity 
      ! svhins: Scaled inboard plasma flow velocity                        
      ! seyins: Scaled inboard electric field                              
      !  (deltat.deltat.zi.e/mp/mi is scale factor)  
      !               
      ! Note: 1.56e4 is associated with a velocity calculated from
      !     v = sqrt ( 8 k T / (PI m) ) 
      !
      ! 1.38e4 is associated with a velocity calculated from 
      !     v = sqrt ( 2 kT / m ) 
      fvycol = qtim * 1.56e+04 / sqrt(crmi)                                     
      if (ciopte.eq.1.or.ciopte.eq.3.or.ciopte.eq.6.or.ciopte.eq.8
     >     .or.ciopte.eq.9.or.ciopte.eq.13) then        
        vy0 = 1.56e4 * sqrt (ctemsc / crmi)

      ! jdemod - all the injection options are pretty messed up but using
      ! a "-" sign on the velocity which changes sign due to porm doesn't
      ! make much sense. ciopte=12 changed to vy0=0
      !ELSEIF  (CIOPTE.EQ.12) THEN
      !  VY0 = -1.56E4 * SQRT (CTEMSC/CRMI)

      elseif (ciopte.eq.4) then

        ! ctemsc - Usually used for ion temperatures
        ! cengsc - Usually used for neutral energies
        ! For injection option 4 the energy of the ions is specified
        ! by the value entered in cengsc
        ! David Elder, 1990, Feb 15 
        vy0  = 1.38e4 * sqrt(cengsc / crmi) 
        vy02 = 1.38e4 * sqrt(cein2 / crmi)
      else                                                          
        vy0 = 0.0                                                               
      endif                                                                     

      ! slmod begin - now defunkt i think - april 28, 97
      if (slopt.eq.1) then
        do alpha = 0.0, ca, ca / real(nqxsi) 
          iqx = min(int(alpha * xscali) + 1, nqxsi)                          
          write(0,*) 'Before:', iqx, alpha, ca, svhins(iqx)
          if (alpha.lt.catin.or.alpha.gt.0.02) svhins(iqx) = 0.0
          write(0,*) 'After :', iqx, alpha, ca, svhins(iqx)
          write(0,*) ' '
        enddo
      endif
      ! slmod end
                                                                               
      ! Check cplsma, point outside which we have special 
      ! characteristics. If not required, artificially set so that it 
      ! never arises.                                                                                               
      if (cioptb.eq.0.and.cioptc.eq.0.and.cioptd.eq.0) then
        cplsma = 2.0 * caw         
      endif
      jx = ipos (cplsma * 0.999, xs, nxs-1)                                       
                                                                               
      ! Set up monitoring variables                                                                                                              
      rstmin = ctimsc / qtim                                                    
      rstmax_win= ctimsc_win /qtim

      if (nizs.gt.0) call monini (rstmin, cstmax, nizs, ctbin)                  
                                                                               
      ! Zero arrays
      call rzero (tag2,   maximp)
      call dzero (ddlims, maxnxs*(2*maxnys+1)*(maxizs+2))                       
      call dzero (ddts,   maxnxs*(2*maxnys+1)*maxizs)                           
      call dzero (ddys,   maxnxs*(2*maxnys+1)*maxizs)                           
      call rzero (deps,   maxnxs*(maxizs+1)*3)
      call rzero (neroxs, maxnxs*5*3)                                           
      call rzero (neroys, maxos*6)                                              
      call rzero (nerods, maxos*5)                                              
      nerods3 = 0.0
      !call rzero (nerods3, maxos*5*(2*maxnps+1))                                              
      call rzero (walls,  (2*maxnys+1)*(maxizs+4))                              
      call rzero (tizs,   maxnxs*(2*maxnys+1)*(maxizs+2))                       
      call rzero (zeffs,  maxnxs*(2*maxnys+1)*6)                                
      call dzero (ddlim3, maxnxs*(2*maxy3d+1)*(maxizs+2)*(2*maxnps+1))          
      call rzero (tiz3,   maxnxs*(2*maxy3d+1)*(maxizs+2)*(2*maxnps+1))          
      call rzero (sdtxs,  maxnxs*maxizs)                                        
      call rzero (sdtys,  (2*maxnys+1)*maxizs)                                  
      call rzero (sdtzs,  maxizs)                                               
      call rzero (sdyxs,  maxnxs*maxizs)                                        
      call rzero (sdyys,  (2*maxnys+1)*maxizs)                                  
      call rzero (sdyzs,  maxizs)                                               
      call dzero (douts,  maxizs*10)                                               
      call rzero (rions,  maxizs)                                               
      if (imode.ne.2) then                                                     
        call rzero (lim5, maxnxs * (2 * maxy3d + 1) * (maxizs + 2) * 
     >    (2 * maxnps + 1) * maxnts)                    
      endif                           
      call rzero (svybar,2*maxqxs+1)
      call rzero (svyacc,2*maxqxs+1)  
                                                                                                   
      ! Print selected parameters set up on first iteration only                                                                                                      
      if (iter.eq.1) then                                                       
        if (imode.eq.0) then                                                
          call prc ('OPERATION MODE 0  COMBINED IMPULSE & STEADY STATE')         
        elseif (imode.eq.1) then                                                
          call prc ('OPERATION MODE 1  IMPULSE')                                 
        elseif (imode.eq.2) then                                                
          call prc ('OPERATION MODE 2  STEADY STATE')                            
        endif                                                                   
                                                                            
        call pri ('  NO OF X BINS                     ', nxs)                   
        call pri ('  NO OF Y BINS                     ', 2 * nys)                 
        call pri ('  NO OF Y BINS FOR 3D RESULTS      ', 2 * ny3d)                
        call pri ('  NO OF P BINS                     ', 2 * maxnps + 1)            
        call pri ('  NO OF X PTS FOR OUTBOARD FACTORS ', nqxso)                 
        call pri ('  NO OF X PTS FOR INBOARD FACTORS  ', nqxsi)                 
        call pri ('  NO OF Y POINTS FOR FACTORS       ', nqys)                  
        call pri ('  NO OF T POINTS FOR TIME RESULTS  ', nts)                   
        call pri ('  MAXIMUM IONIZATION STATE         ', nizs)                  
        call pri ('  NO OF IMPURITY IONS TO FOLLOW    ', nimps)                 
        if (impadd.gt.0) then 
          call pri ('  NO OF EXTRA IMPURITY NEUTALS     ', impadd)
          call prr ('     LAUNCHED ON Y =           +/- ', ycfadd)
          if (cexneut.eq.0) then
            call prc ('     WITH A UNIFORM DISTRIBUTION') 
          elseif (cexneut.eq.1) then 
            call prc ('     WITH A NORMAL DISTRIBUTION') 
          endif
        endif  
        call prr ('  SIZE OF NEUT QUANTUM TIMESTEP (S)', fsrate)                
        call prr ('  BASE VALUE FOR LIM TIMESTEP   (S)', qtim)                  
        write (7,'(1X,''  RANDOM NUMBER SEED'',10X,I15)') ciseed                
                                                                             
        if (igeom.eq.0) then                                                
          call prc ('  RADIAL GEOMETRY OPTION              SLAB')               
        elseif (igeom.eq.1) then                                                
          call prc ('  RADIAL GEOMETRY OPTION              CYLINDER')           
        endif                                                                   
                                                                               
        if (cftcut.gt.caw)                                                      
     >    call prr ('  STOP CROSS FIELD TRANSPORT AT X =', cftcut)              
                                                                               
        if     (cprint.eq.0) then                                               
          call prc ('  PRINT OPTION                        REDUCED')            
        elseif (cprint.eq.1.or.cprint.eq.9) then                                               
          call prc ('  PRINT OPTION                        FULL')               
        endif                                                                   
                                                                               
        if (abs(cthetb-90.0).lt.1.0e-3) then                                    
          call prr ('  LIMITER GEOMETRY POLOIDAL: THETAB', cthetb)              
        else                                                                    
          call prr ('  LIMITER GEOMETRY TOROIDAL: THETAB', cthetb)              
        endif                                                                   
                                                                               
        if (cxspls(1).lt.ca) then                                                   
          call prc ('  SPLITTING & ROULETTING IN OPERATION')                    
          split = .true. 
        else
          split = .false.
        endif  
                                                                              
        call prr ('  ELONGATION PARAMETER OUTBOARD  KO', cko)                   
        call prr ('  ELONGATION PARAMETER INBOARD   KI', cki)                   
                                                                               
        if (canal.lt.ca)                                                        
     >    call prr ('  ANALYTIC EXTENSION INBOARD OF X =', canal)               
                                                                               
        if (cprint.eq.1.or.cprint.eq.9) then                                                   
          call prr ('  SMALLEST X BIN SIZE  (M)         ', xwidm)               
          call prr ('  DELTA X INBOARD FOR FACTORS  (M) ', 1.0 / xscali)          
          call prr ('  DELTA X OUTBOARD FOR FACTORS  (M)', 1.0 / xscalo)          
          call prr ('  SMALLEST Y BIN SIZE  (M)         ', ywidm)               
          call prr ('  DELTA Y ALONG X=0 FOR FACTORS (M)', 1.0 / yscale)          
          call prr ('  STOP TIME FOR ITERATION  (S)     ', timmax)              
          call prr ('  ALLOCATED CPU TIME  (S)          ', cpulim)              
          call prr ('  "NEAR LIMITER" MEANS  0.0 <= X <=', cxnear)              
          call prr ('                AND YP- BETWEEN +/-',cynear/clfact)        
          call prb                                                              
          call prc ('BOUNDARIES OF X BINS')                                     
          write (7,'(8(1X,F7.4))') caw,(xs(ix),ix=1,nxs)                        
          call prb                                                              
          call prc ('BOUNDARIES OF Y BINS (MIRRORED FOR -Y REGION)')            
          write (7,'(8F8.4)') 0.0,(ys(iy),iy=1,nys)                             
          call prb                                                              
          call prc ('BOUNDARIES OF P BINS')                                     
          write (7,'(8(1X,A7))') (prinps(ip),ip=-maxnps-1,maxnps)               
        endif                                                                   
                                                                               
        if (imode.ne.2) then                                                    
          call prb                                                              
          call prc ('DWELL TIME FACTORS AND SPECIFIC OUTPUT TIMES (S)')         
          write (7,9020) (dwelfs(it), it=1, nts)                                  
!          do 30 iz = 0, nizs   
          do iz = 0, nizs                                                 
!            if (iz.eq.0.and.cneuta.ne.0) goto 30 
            if (iz.eq.0.and.cneuta.ne.0) cycle                             
            write (7,9021) iz, (dwelts(iz) * dwelfs(it), it=1, nts)                  
            write (6,9022) iz, (ctimes(it, iz), it=1, nts)                          
!   30     continue
          end do                                                              
        endif
        
        if (vary_2d_bound.eq.1) then
          call prc ('3D ABSORBING BOUNDRIES IN USE')
        endif                                                                   
                                                                               
        rizb = real (cizb)                                                      

        ! Place call to initialize the sputtering yield data so that 
        ! it is available for ion as well as neut launch cases.
        ! Load yield common block with appropriate data.
        if (csputopt.eq.1.or.csputopt.eq.7.or.csputopt.eq.8) then
!          call syield(matlim,mat1,mat2,cneutd,cbombf,cbombz,cion,
!     >            cizb,crmb,ctsub)   
          call syield (matlim,mat1,cneutd,0,cbombf,cbombz,1.0,cion,cizb,
     >      crmb,cebd,csputopt)
        else if (csputopt.eq.2) then
          call syld93 (matlim,mat1,cneutd,0,cbombf,cbombz,1.0,cion,cizb,
     >      crmb,cebd)
        else if (csputopt.eq.3.or.csputopt.eq.4.or.csputopt.eq.5.or.
     >    csputopt.eq.6) then
          call syld96 (matlim,mat1,cneutd,0,cbombf,cbombz,1.0,cion,cizb,
     >      crmb,cebd)
          call init_eckstein_2007(matlim,mat1)
        endif

        ! Set up mat2 if cneutd.eq.2 - this is code from the original 
        ! lim version of syield.
        call syield_set_mat2(mat2,cneutd,cbombf,cbombz)

        ! Write out materials in use.
        write(0,*) 'Materials'
        write(0,'((1x,a15,1x,i3),(3x,a15,1x,f6.2))') 
     >    'Data Option:', csputopt, 'Incident Angle:', 
     >    extra_sputter_angle
        write(0,'((1x,a15,1x,i3),(3x,a15,1x,i6))') 
     >    'Plasma Mat1:', mat1, 'Plasma Mat2:', mat2
        write(0,'((1x,a15,1x,i3),(3x,a15,1x,f6.3))') 
     >    'Limiter Mat:', matlim, 'Binding Energy:', cebd

        ! Comment needed.
        call test_phys_yld(matlim,mat1)

        ! Comment needed.
        cqpl =  122. * exp (-9048. / ctsub)                                       
        cqsl = 1014. * exp (-9048. / ctsub)                                       
        call prdata (nizs, xscalo, xscali, nnbs, ntbs, ntibs, nymfs)
                                                                               
        ! Convert csnorm from degree into radians (prdata prints it out) 
        ! Comment needed: what is csnorm?                                                                                      
        csnorm = csnorm / raddeg                                                
        if (cprint.eq.1.or.cprint.eq.9) then                                                   
          call prb                                                              
          call prc ('SIMPLE FACTORS     ')                                      
          call prr ('  FACTOR FOR POST COLLISION VELOCITY  ', fvycol)           
          call prr ('  FACTOR FOR POLOIDAL DIFFUSION       ',                   
     >                      polods(0) / sqrt(qs(0)))                              
          call prr ('  FACTOR FOR INBOARD PLASMA FLOW VEL. ',                   
     >                      svhins(0) / qs(0))                                    
          call prr ('  FACTOR FOR INBOARD E, IZ STATE 1    ',                   
     >                      seyins(0,1) / (qs(0) * qs(0)))                          
        endif                                                                   
                                                                               
        call taupr1 (qtim, nizs)
      
      ! End of first iteration to-do list.                                                 
      endif                                                                     
                                                                               
      !-----------------------------------------------------------------    
      ! Launch primary neutrals en masse                                           
      !-----------------------------------------------------------------       
                                                                             
      ! tneut : Total number of neutrals launched                                 
      ! tatiz : Total no of ions created                                          
      ! twall : Total no of ions reaching wall                                    
      ! twalln: Total no of neutrals reaching wall                                
      ! tdep  : Total no of ions deposited on limiters                            
      ! ttmax : Total no of neutrals existing at tmax                             
      ! tcent : Total no of neutrals reaching centre                              
      ! tbyond: Total no of ions ionised beyond limit                             
      ! tbelow: Total no of ions recombining to neutrals                          
      ! tcut  : Total no of ions existing at tmax                                 
      ! tstruk: Total no of neutrals striking limiter                             
      ! tfail : Total no of failed neutral launches                                                                                                        
      tneut  = 0.0                                                              
      tres   = 0.0                                                              
      tatiz  = 0.0                                                              
      twall  = 0.0                                                              
      twalln = 0.0                                                              
      tdep   = 0.0                                                              
      ttmax  = 0.0                                                              
      tcent  = 0.0                                                              
      tbyond = 0.0                                                              
      tbelow = 0.0                                                              
      tcut   = 0.0                                                              
      tstruk = 0.0                                                              
      tfail  = 0.0                                                              
                                                                               
      ! tsplit, trulet, nsplit, nrulet: For table of splitting data              
      ! mput: Greatest size of splitting bank.                                                                                                                   
      iput = 0                                                                  
      mput = 0                                                                  
      do is = 1, maxins                                                      
        tsplit(is) = 0.0                                                        
        trulet(is) = 0.0                                                        
        nsplit(is) = 0                                                          
        nrulet(is) = 0                                                          
      end do                                                                 
                                                                                                                                                  
      ! Fit interpolating curve to set of yield modifier values (YMF)                   
      ! Calculate interpolated/extrapolated ymf at each outboard x posn.          
      ! Different yield modifiers for each side of y = 0                          
      ! Flag determines whether to apply to primaries, secondaries, both          
      
      ! Set primary (cymfps) and secondary (cymfss) YMF to 1.0 to start. 
      do j = 1, 2                                                            
        do iqx = 1-nqxso, 0                                                  
          cymfps(iqx,j) = 1.0                                                   
          cymfss(iqx,j) = 1.0  
        end do                                                 
      end do                                                                 
                                                                              
      ! Calculating YMFs
      ! cymfs(in,1) - X coordinate 
      ! cymfs(in,2) - Y < 0 side data
      ! cymfs(in,3) - Y > 0 side data
      if (cymflg.ne.-2) then                                                    
        call fitter (nymfs, cymfs(1,1), cymfs(1,2), nqxso, qxs(1-nqxso), 
     >    cymfps(1-nqxso,1), 'LINEAR')             
        call fitter (nymfs, cymfs(1,1), cymfs(1,3), nqxso, qxs(1-nqxso), 
     >    cymfps(1-nqxso,2), 'LINEAR')             
      endif                                                                     
                                                                            
      if (cymflg.ne.-1) then                                                    
        call fitter (nymfs, cymfs(1,1), cymfs(1,2), nqxso, qxs(1-nqxso), 
     >    cymfss(1-nqxso,1), 'LINEAR')             
        call fitter (nymfs, cymfs(1,1), cymfs(1,3), nqxso, qxs(1-nqxso), 
     >    cymfss(1-nqxso,2), 'LINEAR')             
      endif                                                                     

      ! jdemod - if specific self-sputtering yields are specified they 
      ! override the cymflg flag specification. 
      if (ss_nymfs.gt.0) then
         write(6,'(a)') 'Specified self-sputtering'//
     >     ' yield modifiers will be used:'
         do in = 1, nymfs
           write(6,'(a,i6,10(1x,g12.5))') 'CYMFS:', in, cymfs(in,1),
     >       cymfs(in,2),cymfs(in,3)
         end do

        ! Apply yield modifiers to cymfss
        call fitter (ss_nymfs, ss_cymfs(1,1), ss_cymfs(1,2), nqxso, 
     >    qxs(1-nqxso), cymfss(1-nqxso,1), 'LINEAR')             
        call fitter (ss_nymfs, ss_cymfs(1,1), ss_cymfs(1,3), nqxso, 
     >    qxs(1-nqxso), cymfss(1-nqxso,2), 'LINEAR')             
      endif

      ! Write out YMFs and self-sputtering YMFs.
      do in = 1,nymfs
         write(6,'(a,i6,10(1x,g12.5))') 'CYMFS:',in,cymfs(in,1),
     >     cymfs(in,2),cymfs(in,3)
      end do
      do in = 1,ss_nymfs
         write(6,'(a,i6,10(1x,g12.5))') 'SS_CYMFS:',in,ss_cymfs(in,1),
     >     ss_cymfs(in,2),ss_cymfs(in,3)
      end do

!      do iqx=1-nqxso,0
!         write(6,'(a,i8,10(1x,g12.5))') 
!     >        'YMFS:',nymfs,qxs(iqx),CYMFPS(iqx,1),cymfps(iqx,2),
!     >            cymfss(iqx,1),cymfss(iqx,2)
!      end do

      ! Call neut to simulate neutral particle production.
      ! rneut1: Number of primary neutrals launched                               
      ! rneut : Number of neutrals launched in current batch                      
      ! ratiz : Number of ions created in current batch                           
      ! rstruk: Number of neutrals striking limiter                               
      ! status: 1,2,..cmaxgens for primary,2nd,3rd launches etc, 31 
      !   for total. 
      if (cneuta.eq.0) then                                                     
        status = 1                                                              
        write (0,9012) '***  LAUNCHING ',what(status),' NEUTRALS  ***'          
        write (6,9012) '***  LAUNCHING ',what(status),' NEUTRALS  ***'          
        write (7,9012) '***  LAUNCHING ',what(status),' NEUTRALS  ***'          

        call neut (natiz, fsrate, rres, icut, matlim, mat1, mat2, nimps,
     >    impadd, impcf, qtim, gtot1, gytot1, rstruk, ratiz, rneut, 
     >    rwalln, rcent, rtmax, seed, nrand, neutim, rfail, nymfs, ncvs,
     >    status)         
        write(0,'(a20,i9)')'After NEUT: natiz',natiz
        if (natiz.eq.0) goto 806                                                
      else                                                                      
        status = 1                                                             
        natiz  = nimps                                                          
        ratiz  = real (nimps)                                                   
        rneut  = 0.0                                                            
        rwalln = 0.0                                                            
        rcent  = 0.0                                                            
        rtmax  = 0.0                                                            
        rstruk = 0.0                                                            
        rfail  = 0.0                                                            
        rres   = 0.0                                                            
      endif                                                                     

      rneut1 = rneut                                                            
      rres1  = rres                                                             

      ! Set up minimum velocity by multiplying by qtim. This  minimum is 
      ! constant over space (i.e. not multiplied by the qs(iqx) factor).
      svymin = csvymin * qtim  
      
!-----------------------------------------------------------------------        
!    Follow ions to absorption / eventual fate ...                              
!    This continuation point is taken after secondary neutrals are              
!    launched from those ions that struck the limiters.                         
!-----------------------------------------------------------------------                                                                                   
  200 continue                                                                  

      ! jdemod - moved inside the main loop so that self sputtering is 
      ! fractioned off properly as well. 

      ! Comment needed (this is an slmod).
      ioncnt  = 0
      if (cneuta.eq.0) then 
         ionpnt  = 0.1 * real(natiz)
      else
         ionpnt  = 0.1 * real(nimps)
      endif
                                                                              
      ! Set up tau parallel, heating, stopping: may depend on ctemsc                
      ! value calculated in launch/neut. Print sample values.                                                                                                  
      if (nizs.gt.0 .and. iter.eq.1) then
          call tauin2 (qtim,nizs)                    

        ! The below is just to print out the forces. They aren't applied
        ! to the impurity here.      
        if (cprint.eq.9) then
          ciz = nizs
          
          write(6,'(a,10(1x,g12.5))') 'Force balance:', calphe(ciz),
     >      cbetai(ciz)
         
          if (vary_2d_bound.eq.0) then
            write(6,'(a6,3(2x,a4),a6,40a13)') 'IX','IY','IQX','IQY',
     >        'IQYTMP','XOUT','YOUT','FEG','FIG','FF','FE','FVH','FF2',
     >        'FE2','fvh2','FTOT1','FTOT2','TEGS','TIGS','CFSS',
     >        'CFVHXS','VP1','VP2','FFB','FEB','CVHYS','CEYS','TE','TI',
     >        'NE','VELB','CVHYS2'
     
          ! Arrays that use vary_2d_bounds have additional dimension for 
          ! ip. Add column for ip index.
          else
            write(6,'(a6,4(2x,a4),a6,40a13)') 'IP','IX','IY','IQX',
     >        'IQY','IQYTMP','XOUT','YOUT','FEG','FIG','FF','FE','FVH',
     >        'FF2','FE2','fvh2','FTOT1','FTOT2','TEGS','TIGS','CFSS',
     >        'CFVHXS','VP1','VP2','FFB','FEB','CVHYS','CEYS','TE','TI',
     >        'NE','VELB','CVHYS2'
          endif
          
          do ix = 1,nxs
            write(6,*) 'Static forces:',ix
            do iy = -nys,nys
              iqx = iqxs(ix) 
              if (y.lt.0.0) then 
                iqy_tmp = max(min(int((youts(iy) + ctwol) * yscale) + 1,
     >            nqys), 1)
              else
                iqy_tmp = max(min(int(youts(iy) * yscale) + 1, nqys), 1)
              endif

              y = youts(iy)
              if (iqx.le.ixout) then 
                if (y.gt.0.0) then 
                  iqy = int ((y-qedges(iqx,2)) * cyscls(iqx)) + 1                    
                else
                  iqy = int((-y-qedges(iqx,1)) * cyscls(iqx)) + 1                     
                endif
              else
                iqy  = iqy_tmp
              endif

              if (vary_2d_bound.eq.0) then
             
                ! Forces in the first pzone region (no collector probe).
                pz1 = 1
                feg = calphe(ciz) * ctegs(ix,iy)
                fig = cbetai(ciz) * ctigs(ix,iy)
                ff  = (cfss(ix,iy,ciz) * (cfvhxs(ix,iy)
     >            * velplasma(ix,iy,pz1) - 0.0))
                fe  = (cfexzs(ix,iy,ciz) * efield(ix,iy,pz1))
                fvh = cfvhxs(ix,iy) * velplasma(ix,iy,pz1)

                ! Forces in the second pzone region (collector probe).
                if (maxpzone.gt.1) then 
                  pz2 = 2
                  ff2 = (cfss(ix,iy,ciz) * (cfvhxs(ix,iy)
     >              * velplasma(ix,iy,pz2) - 0.0))
                  fe2 = (cfexzs(ix,iy,ciz) * efield(ix,iy,pz2))
                  fvh2 = cfvhxs(ix,iy) * velplasma(ix,iy,pz2)
                else
                  pz2  = 1
                  ff2  = 0.0
                  fe2  = 0.0
                  fvh2 = 0.0
                endif
                
                ! Bit of a wall of variables, but print them out.
                write(6,'(5i8,40(1x,g12.5))') ix,iy,iqx,iqy,iqy_tmp,
     >            xouts(ix),youts(iy),feg,fig,ff,fe,fvh,ff2,fe2,fvh2,
     >            feg+fig+ff+fe,feg+fig+ff2+fe2,ctegs(ix,iy),
     >            ctigs(ix,iy),cfss(ix,iy,ciz),cfvhxs(ix,iy),
     >            velplasma(ix,iy,pz1),velplasma(ix,iy,pz2),
     >            (cfss(ix,iy,ciz)*(cfvhxs(ix,iy)*cvhys(iqy_tmp)+0.0)),
     >            (cfexzs(ix,iy,ciz)*ceys(iqy_tmp)),cvhys(iqy_tmp),
     >            ceys(iqy_tmp),ctembs(ix,iy),ctembsi(ix,iy),
     >            crnbs(ix,iy),(cfvhxs(ix,iy)*cvhys(iqy_tmp)+0.0),
     >            cvhys(iqy)
     
              ! If using a varying bound then we just copy/paste the 
              ! above and swap out with the respective 3D/4D arrays.
              ! Calculate values at the poloidal middle.
              else
                ip = npbins / 2
                
                ! Forces in the first pzone region (no collector probe).
                pz1 = 1
                feg = calphe(ciz) * ctegs_3d(ip,ix,iy)
                fig = cbetai(ciz) * ctigs_3d(ip,ix,iy)
                ff  = (cfss_4d(ip,ix,iy,ciz) * (cfvhxs_3d(ip,ix,iy)
     >            * velplasma_4d(ip,ix,iy,pz1) - 0.0))
                fe  = (cfexzs_4d(ip,ix,iy,ciz) * 
     >            efield_4d(ip,ix,iy,pz1))
                fvh = cfvhxs_3d(ip,ix,iy) * velplasma_4d(ip,ix,iy,pz1)

                ! Forces in the second pzone region (collector probe).
                if (maxpzone.gt.1) then 
                  pz2 = 2
                  ff2 = (cfss_4d(ip,ix,iy,ciz) * (cfvhxs_3d(ip,ix,iy)
     >              * velplasma_4d(ip,ix,iy,pz2) - 0.0))
                  fe2 = (cfexzs_4d(ip,ix,iy,ciz) * 
     >              efield_4d(ip,ix,iy,pz2))
                  fvh2 = cfvhxs_3d(ip,ix,iy) * 
     >              velplasma_4d(ip,ix,iy,pz2)
                else
                  pz2  = 1
                  ff2  = 0.0
                  fe2  = 0.0
                  fvh2 = 0.0
                endif
                
                ! Bit of a wall of variables, but print them out.
                write(6,'(6i8,40(1x,g12.5))') ip,ix,iy,iqx,iqy,iqy_tmp,
     >            xouts(ix),youts(iy),feg,fig,ff,fe,fvh,ff2,fe2,fvh2,
     >            feg+fig+ff+fe,feg+fig+ff2+fe2,ctegs_3d(ip,ix,iy),
     >            ctigs_3d(ip,ix,iy),cfss_4d(ip,ix,iy,ciz),
     >            cfvhxs_3d(ip,ix,iy),velplasma_4d(ip,ix,iy,pz1),
     >            velplasma_4d(ip,ix,iy,pz2),(cfss_4d(ip,ix,iy,ciz)*
     >            (cfvhxs_3d(ip,ix,iy)*cvhys(iqy_tmp)+0.0)),
     >            (cfexzs_4d(ip,ix,iy,ciz)*ceys(iqy_tmp)),
     >            cvhys(iqy_tmp),ceys(iqy_tmp),ctembs_3d(ip,ix,iy),
     >            ctembsi_3d(ip,ix,iy),crnbs_3d(ip,ix,iy),
     >            (cfvhxs_3d(ip,ix,iy)*cvhys(iqy_tmp)+0.0),cvhys(iqy)              
                  
              endif
             end do
           end do 
         endif
       endif 

      ! Additional printing.
      if (nizs.gt.0) call taupr2 (qtim,nizs)                                    
                                                                               
      ! Initialize (average?) values for location info.                                                                                                                
      nprod  = 0                                                                
      avxpos = 0.0                                                              
      avapos = 0.0                                                              
      avypos = 0.0                                                              
      avppos = 0.0                                                              
      do j = 1, 2                                                           
        rwall (j) = 0.0                                                         
        rdep  (j) = 0.0                                                         
        yldtot(j) = 0.0                                                         
        ythtot(j) = 0.0                                                         
        yldmax(j) = 0.0                                                         
      end do      
  
      ! rdifft: Time to first diffusion                                                 
      randep = 0.0
      rdifft = 0.0      
                          
      ! Print out of which generation of ions we are about to follow.                                   
      if (status.le.50) then 
         write (0,9012) '***  FOLLOWING ',what(status),'   IONS    ***'      
         write (6,9012) '***  FOLLOWING ',what(status),'   IONS    ***'      
         write (7,9012) '***  FOLLOWING ',what(status),'   IONS    ***'      
      else 
         write (0,9013) '***  FOLLOWING ',status,' GENERATION IONS ***'  
         write (6,9013) '***  FOLLOWING ',status,' GENERATION IONS ***'  
         write (7,9013) '***  FOLLOWING ',status,' GENERATION IONS ***'  
      endif
                                                                              
      ! Switch on debug if required.                                                                                                                        
      debugl = .false.                                                          
      if (cstepl.gt.0.0) then                                                   
        write (6,9004) nint(cstepl),qtim                                        
        debugl = .true.                                                         
      endif                                                                     
      debugt = .false.
      if (cstept.gt.0) then
        write (6,9006) cstept                                       
        debugt = .true.                                                         
      endif    
      write (6,*) 'cstept:' , cstept, debugt, chalfl
      cistot = 0.0                                                              
      cismax = 0.0                                                              

      if (debugl) write(78,'(A15)') 'SEYINS'
      if (debugl) write(78,'(G15.4)') seyins(1,1)
      write(78,*) ' '

      if (debugl) write(78,'(A4,A7,13A12)')
     >   'IMP','STEP','Y','SVY','SVYMOD',
     >   'RGAUSS','VPARA','VPARAT','DY1','DY2','QS','DTEMI','QUANT',
     >   'CFSS','YFACT'

      if (debugl)
     >     write(77,'(a,a6,4a8,20a13)') 'Forces:','ix','iy',
     >             'ip','iqx','iqy',
     >     'cist','alpha','y','p','svy','quant','ff',
     >     'fe','feg','fig','ftot','fvh','svg',
     >     'qs','yfact','svymod','spara','delta_y1','delta_y2',
     >     'vpara','vparaqt'

      ! sazmod
      ! Save a little bit of computation time by calculating this 
      ! constant for the exponential 3D injection option.
      ! Update 1/18/22: This option has pretty much been generalized by
      ! the divimp_probs inputs, so is now redundant.
      if (choose_exp.eq.1) then
        choose_exp_fact = choose_exp_lambda * (exp(y0l / 
     >      choose_exp_lambda) - exp(y0s / choose_exp_lambda))
      endif
      
      ! CDF = Cumulative Probability Density Function
      ! PDF = Probability Density Function
      ! Set up CDF for Y injection probabilities between Y0S and Y0L.
      !write(0,*) 'ndivimp_probs = ',ndivimp_probs
      if (ndivimp_probs.gt.0) then
      
        ! Before doing anything, scale the S values (the distance from 
		! target values) to the distance between Y0S and Y0L. Ideally
		! they line up exactly and the scaling = 1, but realistically
		! there will always be a little scaling. The smaller the better.
		max_divimp_s = maxval(divimp_probs(1:maxnys, 1:1))
        do in=1, maxnys
          divimp_probs(in, 1) = divimp_probs(in, 1) * (y0l - y0s) 
     >      / max_divimp_s
        end do
        
        ! Now need to weigh the (unnormalized) probabilities by the 
        ! width of each S bin. The following assumes S starts after 
        ! zero, in accordance with DIVIMP output. Therefore the first
        ! bin width is just the first S value (S1 - 0.0 = S1).   
        do in=1, maxnys
          if (in.eq.1) then
            tmp_width = divimp_probs(in, 1)
          else
            tmp_width = divimp_probs(in, 1) - divimp_probs(in-1, 1)
          endif
          divimp_probs(in, 2) = divimp_probs(in, 2) * tmp_width 
        end do
        
        ! Sum of probabilities to normalize PDF in the next step.
        sum_divimp_prob = sum(divimp_probs(1:maxnys, 2:2))
        
        ! Here we go through one entry at a time and create the CDF 
        ! from the PDF, where normalizing the PDF is done on the fly.
        call rzero(yinj_cdf, maxnys)
        cdf_sum = 0.0
        do in=1, maxnys
          yinj_cdf(in) = cdf_sum + divimp_probs(in, 2) / sum_divimp_prob
          cdf_sum = yinj_cdf(in)
        end do          
      endif
!      open(unit=44, file=
!     >   "/fusion/projects/ird/3dlim/zamperinis/results/saz_debug.txt")
                    
      ! Comment needed.                                                             
      implim = 4                                                                
      statim = za02as (1)                                                       
      porm   = -1.0                                                             
      kk     = 1000 * isect                                                     
      kklim  = kk - 10                                                          

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
      !                                                                +         
      !                     For each ion do                            +         
      !															       +
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do 800  imp = 1, natiz                                                                      

        ! Print update every 10% of particles
        !write(0,*) 'Particle: ',imp,'/',natiz
        if ((natiz / 10).gt.0) then 
          if (mod(imp, natiz / 10).eq.0) then 
            perc = int((imp * 10)/(natiz / 10))
              write(0,'(a,i3,a,i8)') 
     >          'Following Ions: ',perc,' % complete. Particle # =',imp
          endif
        endif

        ! Comment needed to explain the role of these variables.
        ! svybit: scaled velocity, to be multiplied by qtim * qs(iqx)             
        !   to form svy a few lines later when iqx is known
        ! tstepl:
        ! porm:
        ! oldy: 
        ! old_y_position:
        svybit = 0.0
        tstepl = cstepl                                                         
        porm   = -1.0 * porm                                                    
        oldy   = 0.0                                                            
        old_y_position = 0.0
        oldalp = 0.0                                                            

        ! jdemod - initialize particle reflection
        call init_part_reflection

        ! Initialize tracking variables - if debugt
        if (debugt) then

           ! Turn track debugging off if it has done the required number 
           ! of particles. 
           if (imp.gt.cstept) then 
              debugt = .false.
           else
              bigtrac = .false.
              traclen = 1
              avgtrac = 0.0
           endif
           ! write(6,*) 'imp:', imp,bigtrac,traclen,cstept,debugt 
        endif 
                                                                               
        ! Set initial characteristics of ion                                      
        ! cx     x position                                                    
        ! y      y position                                                    
        ! p      p position                                                    
        ! svybit scaled velocity, to be multiplied by qtim * qs(iqx)             
        !   to form svy a few lines later when iqx is known ...           
        ! ciz    ionization level                                              
        ! cist   scaled time - now a real                                      
        ! itime  number of next measurement time                               
        ! ix     x bin ion in                                                  
        ! iy     y bin ion in                                                  
        ! ip     p bin ion in                                                  
        ! iqx    array offset for reading factors                              
        ! alpha  x position translated to x axis using elongation data         
        !                                                                        
        ! use data from neut to inject particles (+/- already accounted 
        ! for)           
                                                                              
        ctemi = ctemsc                                                          
        ciz   = cizsc                                                           
        maxciz = cizsc                                                           
        if (cneuta.eq.0) then                                                   
                                                            

          ! This option was added to simulate particle production which 
          ! does not recycle locally - the entire target production is 
          ! introduced into the SOL at a specified y coordinate.
          ! init_y_coord is an optional input with a default value 
          ! of 0.0.
          if (init_y_coord.ne.0.0) then 
             y = sign(init_y_coord, yatizs(imp))
          else
             y = yatizs(imp)                                                   
          endif

          ! Extract starting X, Y (above), P, scaled velocity and 
          ! particle weight values for the ions as calculated in neut.
          ! sputy: This is a very Monte Carlo-esque variable. It is the
          !   particle's weight, which can be changed in a few ways. If
          !   a particle is injected as an ion option, circumventing
          !   neut and all the sputtering code, then it is started with
          !   a weight of sputy = 1. The particle represents 1 particle.
          !   But if this particle were to cause self-sputtering, then
          !   the sputtered particle would be sputy = sputy * yield.
          !   The particle thus represents less than a single particle.
          !   sputy is also changed when rouletting and when the 
          !   particle is split. See further discussion where ddlim3
          !   is scored for an example.
          cx     = xatizs(imp) 
          p      = patizs(imp)                                                   
          svybit = vins(imp)                                                     
          sputy  = sputys(imp)                                                   

		  ! Comment needed (this is an slmod).
          ioncnt     = ioncnt + 1.0      
          tloss(imp) = 0.0
          if (ioncnt.gt.ionpnt) then
            write(0,'(a,i5,a)') ' ',int(ionpnt/real(natiz)*100.0),' %'
            ionpnt = ionpnt + 0.1*real(natiz)
          endif
      
        ! Use injection coordinates specified in datafile. Launch 
        ! alternately on +y,-y regions. Set initial ionisation state.                                                                                        
        elseif (cneuta.eq.1) then                                               
           cx  = cxsc                                                            

          ! svybit is specified first : for injection 4 if a random 
          ! number is greater than the probability that the velocity is 
          ! vy0 then svybit is set to vy02. direction alternates +/-. 
          if (ciopte.eq.1.or.ciopte.eq.13) then 
            call surand (seed,1,ran)
            svybit= vy0 * sign(1.0,ran-0.5)
          else
            svybit= vy0 * porm                                                    
          endif   

          if (ciopte.eq.0.or.ciopte.eq.1) then
            y = cysc
          elseif (ciopte.eq.5.or.ciopte.eq.6.or.
     >        ciopte.eq.12.or.ciopte.eq.13) then                   
            y = cysc * porm                                                     
          elseif (ciopte.eq.4) then
            call surand (seed,1,ran)
            y = cysc*2.0*(ran-0.5)
            call surand (seed,1,ran) 
            if (ran.gt.cprob)  svybit = vy02 * porm 
            nrand = nrand +2    
          elseif (ciopte.eq.9) then 

            ! Overwrite the assigned x injection value of cxsc,
            ! assign the y-injection value
            y = cysc * porm
            call surand(seed,1,ran)
            cx = (x0l - x0s) * ran + x0s 
            nrand = nrand + 1
            
          ! Gaussian ion injection option
          elseif (ciopte.eq.10) then

250         call surand(seed,1,ran)
            nrand = nrand + 1
            y = 6.0e-2 * ran * porm
          
            call surand(seed,1,ran)
            nrand = nrand + 1
            if (ran.gt.exp(-((y / 0.02)**2))) goto 250

            ioncnt = ioncnt + 1.0      
            if (ioncnt.gt.ionpnt) then
              write(0,'(a,i5,a)') ' ',int(ionpnt/real(nimps)*100.0),' %'
              ionpnt = ionpnt + 0.1*real(nimps)
            endif

          ! 3D Gaussian ion launch
          elseif (ciopte.eq.11) then
            write(0,*)'Error! CIOPTE = 11 not implemented.'
            stop

          else                                                                  
            call surand (seed, 1, ran)                                          
            y = ctwol * ran * porm                                              
            nrand = nrand + 1                                                   
            yy = mod (abs(y), cl)                                               
            if (yy.gt.chalfl) yy = cl - yy                                      
            if (cx.ge.0.0) then                                                 
              cx = cx * (yy*coni + 1.0)                                         
            else                                                                
              cx = cx * (yy*cono + 1.0)                                         
            endif                                                               
          endif                                                                 
          p     = cpsc                                                          
          sputy = 1.0                                                           
                                                                              
        ! Simulate neut using a rectangular injection region                                                                                                   
        elseif (cneuta.eq.2) then                                               
          call surand (seed, 1, ran)                                            
          cx  = (x0s + ran * (x0l - x0s))                                         
          call surand (seed, 1, ran)                                            
          y   = (y0s + ran * (y0l - y0s))                                 
          p   = cpsc                                                            
          nrand = nrand + 2                                                     
          svybit= vy0 * porm                                                    
          sputy = 1.0      
                                                               
        elseif (cneuta.eq.3) then                                               

         ! Allow for injection over a 3D volume
          call surand (seed, 1, ran)                                            
          cx  = (x0s + ran * (x0l - x0s))
          call surand (seed, 1, ran)                                            
          p   = (p0s + ran * (p0l - p0s))                                  
          
          ! Choose uniformly between Y0S and Y0L. 
          call surand (seed, 1, ran)     
          if (choose_exp.eq.0) then                                                
            y = (y0s + ran * (y0l - y0s))   
          else          

            ! 1/18/22 Update: Generalized by the divimp_prob option. Use 
            ! that instead.
            ! Choose from exponential. Equation below is from choosing
            ! directly from the pdf exp(y/lambda). I.e., normalize this
            ! pdf, then find the cdf, then set it equal to random number
            ! and solve for y.
            y = choose_exp_lambda * log(ran * choose_exp_fact / 
     >          choose_exp_lambda + exp(y0s / choose_exp_lambda))
            !write(69,*) y
          endif
          
          if (ndivimp_probs.gt.0) then
          
			! Probably a more efficient way to do this, but I'm a 
			! spoiled python programmer and don't know a better way.
            do in=1, maxnys
            
              ! Once the random number is less than the CDF, take
              ! the last encountered S value (which is the distance from 
              ! one of the targets) as the injection location. Offset it
              ! by y0s to line it up with the simulation domain.
              if (ran.lt.yinj_cdf(in)) then
                
                ! To prevent choosing at the discrete S values provided 
                ! in divimp_probs, we will uniformally choose between 
                ! the chosen value and the one before it. 
                call surand (seed, 1, ran)
                if (in.eq.0) then
                  y = y0s + ran * divimp_probs(in, 1)
                else
                  y = y0s + divimp_probs(in-1, 1) + ran * 
     <              (divimp_probs(in, 1) - divimp_probs(in-1, 1))
                end if
                nrand = nrand + 1
                !write(44,*) ran,yinj_cdf(in),y0s,divimp_probs(in,1),y 
                exit
              endif
            end do
          endif                                       
          
          nrand = nrand + 3
     
          ! Set initial velocity to range of -vel to +vel assigned 
          ! randomly.
          vy0 = 1.56e4 * sqrt(ctemsc / crmi)
          call surand (seed, 1, ran)
          svybit = vy0 * (2.0 * ran - 1.0)
          !svybit= vy0 * porm                                                    
          sputy = 1.0                                                           

        endif                                                                                                                                       

        absy = abs (y)                                                          
        yy = mod (absy, cl)                                                     
        if (yy.gt.chalfl) yy = cl - yy     
        
        ! alpha is the x location after accounting for elongation. If
        ! no elongation then alpha = cx.                                     
        if (cx.ge.0.0) then                                                     
          alpha = cx / (yy * coni + 1.0)                                          
        else                                                                    
          alpha = cx / (yy * cono + 1.0)                                          
        endif                                                                   

        !write(0,'(a,i8,10(1x,g18.8))') 'Inject:',imp,cx,y,p,vy0
        !write(0,'(a,5(1x,g18.10))') 'ALPHA:',ALPHA,CX,YY,CONI,CONO

        dsputy = dble (sputy)                                                   
                                                                               
        ! Spreading: Spread out ions by a partial iteration                            
        ! of frqtim * qtim (note may be negative)                                                                                                            
        if (cneuta.eq.1 .and.(ciopte.eq.5.or.
     >    ciopte.eq.6.or.ciopte.eq.7.or.ciopte.eq.8)) then

          frqtim = (imp - 0.5) / ratiz - 0.5                                    
          if (alpha.ge.0.0) then                                                
            iqx = min (int(alpha * xscali) + 1, nqxsi)                              
          else                                                                  
            iqx = max (int(alpha * xscalo), 1 - nqxso)                              
          endif                                                                 
          call surand (seed, 1, ran)                                            

          ! Decide which y-region the particle is in and then use  
          ! the index to access the x-diff data for the appropriate region
          ! D. Elder Nov 23 1990
          if (y.le.0.0) then                                                    
            if (y.gt.-chalfl) then 
              j = 1
            elseif (y.lt.-c3halfl) then 
              j = 2
            else 
              j = 3
            endif
          else      
            if (y.gt.c3halfl) then 
              j = 1
            elseif (y.lt.chalfl) then 
              j = 2
            else 
              j = 3
            endif
          endif

          cx = cx + frqtim * (sign(cxbfs(iqx, j), cxcfs(iqx, j) - ran) 
     >      + cxafs(iqx, j))  

          nrand  = nrand + 1                                                     
          y      = y + frqtim * svybit * qtim * qs(iqx)                         
          absy   = abs (y)                                                      
          yy = mod (absy, cl)                                                   

          if (yy.gt.chalfl) yy = cl - yy                                        
          if (cx.ge.0.0) then                                                   
            alpha = cx / (yy * coni + 1.0)                                        
          else                                                                  
            alpha = cx / (yy * cono + 1.0)                                        
          endif                                                                 
        endif                                                                   

        ! jdemod - At this point the intial CX, Y, P particle 
        ! coordinates are set. Verify valid P if reflection option is on 
        ! which places bounds on the P value.
        if (preflect_opt.eq.1.and.abs(P).gt.preflect_bound) then 
           write(0,'(a,5(1x,g12.5))') 'WARNING: Particle P coordinate'//
     >                         ' outside of P bound:', P, preflect_bound
           P = sign(preflect_bound, p)
        endif

        ! Get starting indices of the particle's location.
        ix = ipos (alpha, xs, nxs-1)                                         
        iy = ipos (absy,  ys, nys-1)                                         

        ! We have indexed iy using the absolute value of Y (presumably
        ! this must've been quicker back in the day since it meant
        ! searching a smaller array). Since Y is symmetrical, if Y is
        ! negative then just use the negative index.
        if (y.lt.0.0) iy = -iy                                                  
        jy = iabs (iy)                                                      
        ip = ipos (p, ps, 2 * maxnps) - maxnps - 1    
        
        ! What's this? Well, the arrays associated with vary_2d_bound
        ! go from 1 to 2*maxnps+1, sorta by accident but also because 
        ! negative indices are strange by modern conventions (they won't
        ! always correspond to a negative poloidal coordinate either 
        ! since you could in theory have asymmetric poloidal bounds). 
        ! So instead, we are using a different ip to correctly index the 
        ! newer arrays. But note: ip or ip2 refer to the same location, 
        ! it's just a means of correctly indexing the newer arrays.
        ip2 = ipos(p, ps, npbins)
                              
        ! jdemod      
        ! RTIME: Starting time of the particle relative to t=0
        ! CIST: Elapsed time for the specific particle and starts at 0.0
        ! DTIME and DIST: Double precision versions of the same 
        !   variables
        ! Both increment by the timestep. RTIME is used to determine the 
        ! time bin while CTIME is used for particle lifetime statistics. 
        !
        ! Note: Currently time spent as neutrals in the NEUT routine is 
        ! not included in total particle elapsed time. RTIME is ION 
        ! elapsed time from t=0.             
        if (rstmax_win.eq.0.0) then
          rtime = rstmin
        else
          !call surand (seed, 1, ran)                                            
          !nrand = nrand + 1
          !rtime = (rstmax_win - rstmin) * ran + rstmin
          ! Distribute particles evenly over time
          time_frac = (dble(imp - 1) / dble(nimps - 1))
          time_end = dble(rstmax_win)
          time_start = dble(rstmin)
          time_win = time_end - time_start
          dtime = time_win * time_frac + time_start
          rtime = sngl(dtime)
        endif
        it = ipos (rtime, ctimes(1,ciz), nts) 
        !it = ipos (rstmin, ctimes(1, ciz), nts)                              

!        if (1000*(imp/1000).eq.imp) then
!           write(6,'(a,2(1x,i10),5(1x,g12.5),2(1x,i10))')
!     >            'DEBUG:',imp,it,rtime,time_start,time_end,
!     >                time_frac, time_win, imp-1,nimps-1
!        endif
        
        ! Comment needed. What is this for?                      
        is = 1         
                                                              
        ! slmod begin - divimp ion profile
        if (optdp.eq.1) then
          injbint(iy) = injbint(iy) + 1
          !injbinp(ip) = injbinp(ip) + 1
          !ciz = ciz + 1
          !write(0 ,*) 'injection ',y,iy,exp(-(y/0.02)**2),porm
          write(60,*) 'Injection ',y,iy,exp(-(y/0.02)**2),porm
          write(79,*) cx,y,p,ix,iy,ip
        endif
                                                                              
        ! Set initial iqx, iqy values. Reminder alpha = cx if no 
        ! elongation has been specified.                         
        ! iqx: Index of the X location in regards to the more detailed 
        !   set of x coordinates that describe a location along the limiter.
        ! iqy: Similar, but just a more detailed set of y coordinates.                                                                   
        if (alpha.ge.0.0) then                                                  
          iqx = min (int(alpha * xscali) + 1, nqxsi)                                
        else                                                                    
          iqx = max (int(alpha * xscalo), 1 - nqxso)                                
        endif                                                                   
        iqy = 0             
        
        ! Can now finish off svy by multiplying by timestep factor.                                                    
        svy = svybit * qtim                                           
                                                                              
        ! Record injection position.                                                                                                                            
        avxpos = avxpos + cx * sputy                                            
        avapos = avapos + alpha * sputy                                         
        avypos = avypos + absy * sputy                                          
        avppos = avppos + abs(p) * sputy                                        
                                                                               
        ! Check if ion has started above max ionization state.                                                                                                     
        if ((ciz.gt.cion).or.(ciz.gt.nizs)) then                         
          tbyond = tbyond + sputy                                               
          ifate = 9                                                             
          goto 790                                                              
        endif                                                                   

!        write(0,'(a,3i8,10(1x,g18.8))') 'Ion start:',
!     >      ix,iy,ip,cx,y,p,alpha,svybit
                                                                        
        ! If set Ti=Tb for state ciz applies, better do it                                      
        if (vary_2d_bound.eq.0) then                                                              
          if (ciz.eq.cizset) ctemi = max(ctemi, ctembs(ix, iy)) 
        else
          if (ciz.eq.cizset) ctemi = max(ctemi, ctembs_3d(ip2, ix, iy))
        endif                   
                                                                               
        ! Calculate point at which diffusion will be first applied.               
        ! Depends on diffusion option as follows :-                               
        ! 0) Immediate diffusion                                                  
        ! 1) After a time based randomly on initial temperature,                  
        !    -taupara.log$/2, where $ in (0,1)                                   
        ! 2) After time taupara, taking into account changes in taupara           
        !    as ion heats up                                                                                                                                     
        rconst = 1.e20                                                          
        if (cdifop.eq.0) then                                                   
          if (cioptb.ne.1) rconst = 0.0                                         
        elseif (cdifop.eq.1) then            
          if (vary_2d_bound.eq.0) then                                   
            if (cfps(ix,iy,ciz).gt.0.0) then                                      
              nrand = nrand + 1                                                   
              call surand (seed, 1, ran)                                       
              rconst = -ctemi * log (ran) / cfps(ix,iy,ciz) * qs(iqx) 
            endif 
          else
            if (cfps_4d(ip2,ix,iy,ciz).gt.0.0) then                                      
              nrand = nrand + 1                                                   
              call surand (seed, 1, ran)     
              rconst = -ctemi * log (ran) / cfps_4d(ip2,ix,iy,ciz) * 
     >          qs(iqx) 
            endif            
          endif                                                                 
        endif                                                                   
        spara  = 0.0                                                            
        diffus = .false.                                                        
                                                                              
        ! Set initial (x,y) coordinates & Ti in double precision.                 
        ! dy1 accumulates non-diffusive changes, dy2 diffusion changes.                                                                                          
        ! jdemod - changed variables
        ! dy1   = 0.0d0                                                           
        ! dy2   = dble (y)                                                        

        delta_y1 = 0.0d0
        delta_y2 = 0.0d0
        y_position = dble(y)

        ! Comment needed. What is quant for? Also, boo! Bad variable 
        ! name!
        quant = 0.0                                                             
        dtemi = dble (ctemi)                                                    
                          
        ! Record ion as being injected.                                                     
        call monup (8,sputy)                                                    
                                                                                                                                                             
        ! Iterate up to next event in variable steps of qs(iqx)*qtim            
        ! dependent on current x position.  based on multiples of               
        ! "standard" iteration time qtim itself.                                
        ! dist records the iteration number, cist same but single prec.         
                                                                        
        ! jdemod - add option to specify the injection time between 0.0 
        !   and rstmin
        ! cist is particle elapsed time from 0.0 so it always starts 
        !   at 0.0
        ! rtime is the particle time relative to t=0 for the simulation        
        ! rtime is initialized with particle initialization

        cist   = 0.0
        dist   = dble (cist)                                                  
        qfact  = qs(iqx)                                                      
        dqfact = dble (qfact)                                                 

        if (debugl) then                                                      
          write (6,9005)                                                      
          write (6,9003) imp,cist,iqx,iqy,ix,iy,cx,alpha,y,p,svy,ctemi,
     >      spara,sputy,ip,it,is,'ION APPEARED'          
        endif     
          
        ! Before particle tracking starts.                                                                 
                                                                              
          ! Particle tracking loop begin.                                                                  
  500     continue                                                              

          ! Five random numbers are always used for each iteration                
          ! in loop 500. An occasional extra one is required for                 
          ! testing recombination, for rouletting or for determining              
          ! whether to apply a deltay diffusion step. Hence a maximum of         
          ! eight randoms are used in each iteration. The random numbers         
          ! vector thus has to be regenerated whenever we are within 8            
          ! numbers from the end of the vector (ie. kklim).                     
          ! When we jump out of loop 500 a proportion of this vector of           
          ! randoms will be wasted. To prevent this, we only generate            
          ! enough here to replace those used in the last pass through            
          ! loop 500  (ie the value of kk).   
          if (kk.gt.kklim) then                                                 
            call surand (seed, kk, ranv)                                        
            nrand = nrand + kk                                                  
            kk = 0                                                              
          endif                                                                 

          ! Record particle track if track debugging is turned on.
          if (debugt) then 
             tptrac(traclen,1) = cx
             tptrac(traclen,2) = y 
             traclen = traclen + 1
             if (traclen.gt.maxlen) then
                bigtrac = .true.
                traclen = 1
             endif   
             
             ! slmod begin - not sure what this does
             if (traclen.gt.maxlen) avgtrac = 0.0

             if (traclen.gt.2) then
               avgtrac = avgtrac + tptrac(traclen-1,2) - 
     +                             tptrac(traclen-2,2)
             endif

             ! jdemod - the print out of the particle tracks can use a 
             ! lot of space 300Mb+ for a 6 particle debug for example - 
             ! 6 particles from each generation are printed. 
             ! To allow collection of track information in the code 
             ! without the overhead in the LIM file... I have commented 
             ! this out for now... it could be added back with a specific 
             ! print option if that would help.  

!             write(6,'(a,i6,i8,10(1x,g12.5))') 
!     >              'trac:',imp,traclen-1,tptrac(traclen-1,1),
!     >              tptrac(traclen-1,2),
!     >             (tptrac(traclen-1,2)-tptrac(traclen-2,2)),
!     >             AVGTRAC/(TRACLEN-2)
!     >                  tptrac(traclen-1,2),bigtrac
          endif
                                                                              
          ! Check for changes in characteristic times data. Will occur            
          ! when we use a special plasma. Each time an ion enters                
          ! this region new coefficients are calculated based on its              
          ! temperature at that time, providing temp is 10% different             
          ! from the previous time this region was entered.                                                                                     
          if (ix.le.jx) then                                                    
            if (cioptb.ge.2 .or. cioptc.eq.2 .or. cioptd.eq.3) then   
              if (vary_2d_bound.eq.0) then          
                temold = ctolds(ix,ciz)                                           
                if (ctemi.gt.1.1*temold .or. ctemi.lt.0.9*temold) then            
                  call taufix (ix,temold,ctemi)                                   
                  ctolds(ix,ciz) = ctemi                                          
!                 if (debugl) write (6,9003) imp,cist,iqx,iqy,ix,iy,              
!    >              cx,alpha,y,p,svy,ctemi,spara,sputy,ip,it,is,                 
!    >              'fix',0,0,temold                                              
                endif 
              else
                temold = ctolds_3d(ip2,ix,ciz)                                           
                if (ctemi.gt.1.1*temold .or. ctemi.lt.0.9*temold) then       
                  call taufix_3d(ip2,ix,temold,ctemi)                                   
                  ctolds_3d(ip2,ix,ciz) = ctemi 
                endif 
              endif                                                            
            endif                                                               
          endif                                                                 
                                                                               
          ! In most cases, calculate parallel diffusion coefficient           
          ! and move on.  But to start with, each ion must exist for          
          ! "rconst" iterations before parallel diffusion is 
          ! applied. The value of "rconst" depends on which "first 
          ! diffusion" option was chosen 0,1 or 2. Once the ion has 
          ! existed long enough, set "diffus" flag true for subsequent 
          ! iterations.         
          ! Note: cccfps = sqrt(4.88e8/(cfps*crmi))* qtim*qs * ...            
          ! Fixed 14/7/88: If rconst < 1  (ie taupara < deltat), then         
          ! diffusion should be switched on straight away.                                          
          if (vary_2d_bound.eq.0) then                                       
            if (diffus) then                                                  
              spara = ctemi * cccfps(ix,iy,ciz)                               
            else                                                              
              if (cdifop.eq.2) then                                           
                if (cfps(ix,iy,ciz).gt.0.0) then                              
                  rconst = 2.0 * ctemi / cfps(ix,iy,ciz) * qfact              
                else                                                          
                  rconst = 1.e20                                              
                endif                                                         
              endif                                                           
              if (cist.ge.rconst .or. rconst.lt.1.0) then                     
                rdifft = rdifft + cist * qtim * sputy                         
                diffus = .true.                                               
                spara  = ctemi * cccfps(ix,iy,ciz)                            
              else                                                            
                spara  = 0.0                                                  
              endif                                                           
            endif                                                             
                                                                               
            if (ix.le.jx) then                                                
              if (cioptb.eq.3) then                                           
                kk = kk + 1                                                   
                if (ranv(kk).gt.cfps(ix,iy,ciz) / (2.0*ctemi)) 
     >            spara = 0.0        
              elseif (cioptb.eq.4) then
                kk = kk + 1
                if (ranv(kk).gt.(cfps(ix,iy,ciz) / (2.0*ctemi))) then
                  spara = 0.0
                else
                  spara = spara * sqrt(ctemi / ctemsc) 
                  dtemi = dtemi + (dble(ctembsi(ix,iy)) - dtemi) 
     >              * dmin1( dble(dtemi / ctembsi(ix,iy)), 0.5d0)
                  ctemi = sngl(dtemi)
                endif
              endif                                                           
            endif   
              
          ! vary_2d_bound routine.
          else
            if (diffus) then                                                  
              spara = ctemi * cccfps_4d(ip2,ix,iy,ciz)                               
            else                                                              
              if (cdifop.eq.2) then                                           
                if (cfps_4d(ip2,ix,iy,ciz).gt.0.0) then                              
                  rconst = 2.0 * ctemi / cfps_4d(ip2,ix,iy,ciz) 
     >              * qfact              
                else                                                          
                  rconst = 1.e20                                              
                endif                                                         
              endif                                                           
              if (cist.ge.rconst .or. rconst.lt.1.0) then                     
                rdifft = rdifft + cist * qtim * sputy                         
                diffus = .true.                                               
                spara  = ctemi * cccfps_4d(ip2,ix,iy,ciz)                            
              else                                                            
                spara  = 0.0                                                  
              endif                                                           
            endif                                                             
                                                                               
            if (ix.le.jx) then                                                
              if (cioptb.eq.3) then                                           
                kk = kk + 1                                                   
                if (ranv(kk).gt.cfps_4d(ip2,ix,iy,ciz) / (2.0 * ctemi)) 
     >            spara = 0.0        
              elseif (cioptb.eq.4) then
                kk = kk + 1
                if (ranv(kk).gt.(cfps_4d(ip2,ix,iy,ciz)/
     >            (2.0*ctemi))) then 
                  spara = 0.0
                else
                  spara = spara * sqrt(ctemi / ctemsc) 
                  dtemi = dtemi + (dble(ctembsi_3d(ip2,ix,iy))
     >              -dtemi) * dmin1(dble(dtemi
     >              / ctembsi_3d(ip2,ix,iy)),0.5d0)
                  ctemi = sngl(dtemi)
                endif
              endif                                                           
            endif 
          endif                                                          

          if (cioptb.eq.13) then

            ! Velocity diffusion:
            spara  = 0.0
            if (vary_2d_bound.eq.0) then
              vparat = cccfps(ix,iy,ciz)
            else
              vparat = cccfps_4d(ip2,ix,iy,ciz)
            endif

 7702       nrand = nrand + 1

            call surand(seed,1,ran1)
            if (ran1.eq.0.0) goto 7702
            nrand = nrand + 1

            call surand(seed,1,ran2)
            rgauss = sqrt(-2.0* log(ran1))*cos(2.0*pi*ran2)

            vpara = vparat * rgauss

            ! timestep is part of cccfps now
            !svy = svy + vpara * qtim
            svy = svy + vpara 
              
          endif
                                                                              
          ! Update y position of ion.                                                        
          ! jdemod - set Y_position to Y in case code has adjusted 
          ! the Y value outside of the particle movement loop. The 
          ! single precision Y variable holds the definitive version 
          ! of the particle Y position. The update only is performed 
          ! in double precision. 
          y_position = dble(y)
          old_y_position = y_position

          absy = abs(y)
          oldy = y                                                         
          
          ! We will let cynear = 0.0 default to setting yfact = 1. This 
          ! is because if we have absorbing boundaries on, which we
          ! typically do, we almost certainly have already incorporated
          ! the fact that boundaries are field-aligned distances. 
          ! Therefore we do not want to decrease the strength of the
          ! parallel transport, it would be double jeopardy. 
          if (cynear.eq.0) then
            yfact = 1.0
          elseif (absy.le.cynear .or. absy.ge.cyfar) then                       
            yfact = csintb                                                
          else                                                              
            yfact = csintb * clfact                                         
          endif                                                            

          ! Ensure particle stays at or above minimum velocity.
          if (svymin.eq.0.0) then
            svymod = svy
          else
            svymod = sign(max(abs(svy),svymin),svy)
            if (debugl) 
     >        write(6,*) 'SVYMOD,SVY,SVYMIN:',svymod,svy,svymin
          endif

          ! Accumulate velocity information.
          svybar(iqx) = svybar(iqx) + abs(svymod * qs(iqx)) * sputy
          svyacc(iqx) = svyacc(iqx) + sputy

          ! jdemod
          ! Comment needed.
          delta_y1 = dble ((svymod + 0.5 * quant) * qs(iqx) * yfact)             
          !dy1 = dy1 + dble ((svymod + 0.5 * quant)*qs(iqx)*yfact)             

          if (spara.gt.0.0) then                                            
            kk  = kk + 1                                                    
            spara = spara * yfact                                           

            ! jdemod
            delta_y2 = dble (sign(spara, ranv(kk) - 0.5))                    
            !dy2 = dy2 + dble (sign (spara,ranv(kk)-0.5))                    
          endif                                                             

          ! Updating Y coordinate - 
          ! Note: DY2 contains the initial Y coordinate PLUS all 
          !         spatial diffusive steps
          !       DY1 contains all forces and velocity diffusive 
          !         steps   
          !write(0,*) 'Y_position = ', Y_position
          !write(0,*) 'delta_y1   = ', delta_y1
          !write(0,*) 'delta_y2   = ', delta_y2
          y_position = y_position + delta_y1 + delta_y2
          y          = sngl (y_position)                                            

          if (debugl) then
            write(77,'(a,5i8,30(1x,g12.5))') 'Forces:',ix,iy,ip,iqx,iqy,
     >            cist,alpha,y,p,svy,quant,ff,fe,feg,fig,ff+fe+fig+feg,
     >            fvh,svg,
     >            qs(iqx),yfact,svymod,spara,delta_y1,delta_y2,vpara,
     >            vpara*qtim
          endif 
!          write(6,'(a,5i8,30(1x,g12.5))') 'Forces:',ix,iy,ip,iqx,iqy,
!     >            cist,alpha,y,p,svy,quant,ff,fe,feg,fig,ff+fe+fig+feg,
!     >            fvh,svg,
!     >            qs(iqx),yfact,svymod,spara,delta_y1,delta_y2,vpara,
!     >            vpara*qtim  
     
          ! jdemod - Check for Y absorption
          if (yabsorb_opt.ne.0) then 

            ! Choose correct absorption subroutine.
            if (vary_2d_bound.eq.1) then
              call check_y_absorption_2d(ip2, ix, cx, y, oldy, 
     >          sputy, ciz, ierr)
            else
              call check_y_absorption(cx,y,oldy,sputy,ciz,ierr)
            endif

            if (ierr.eq.1) then 
                 
              ! Particle absorbed - exit tracking loop - y absorption
              ifate = 11
              goto 790
            endif 
          endif
             
          ! jdemod
          ! Y-boundary is checked in the inboard/outboard code 
          ! because the constraints are different for the 
          ! different regions. Y boundary checking is not 
          ! desired outboard where a limiter surface is present. 
          ! However - we can check for reflections here. 
          if (yreflection_opt.ne.0) then 
            if (abs(y).gt.ctwol) then 
              write(6,*) 'Y > CTWOL'
              write (string,'(1x,f10.6,f10.5)') oldalp,oldy                       
              write (6,9003) imp,cist,iqx,iqy,ix,iy,                              
     >          cx,alpha,y,p,svy,ctemi,spara,sputy,ip,it,is,string               
            endif

            call check_reflection(cx,y,oldy,svy,sputy,2,debugl,ierr)

            if (ierr.eq.1) then
                 
              ! Write some debugging info.
              write (string,'(1x,f10.6,f10.5)') oldalp,oldy                       
              write (6,9003) imp,cist,iqx,iqy,ix,iy,                              
     >          cx,alpha,y,p,svy,ctemi,spara,sputy,ip,it,is,string               
            endif
          endif

          absy  = abs (y)                                                   
          yy    = mod (absy, cl)                                            
          if (yy.gt.chalfl) yy = cl - yy                                    

          if (debugl) write(78,'(i4,f7.1,13g12.5)') 
     >      imp,cist,y,svy,
     >      svymod,rgauss,vpara*qtim,vparat,delta_y1,delta_y2,
     >      qs(iqx),dtemi,quant,cfss(ix,iy,ciz),yfact
                                                               
          ! Update X position of ion, allowing for elongation                 
          ! Note yycon and alpha are updated in inboard/out loop below                                                                            
          oldalp = alpha                                                    
          if (cx.ge.0.0) then                                               
            yycon = yy*coni + 1.0                                           
          else                                                              
            yycon = yy*cono + 1.0                                           
          endif                                                             
          kk = kk + 1                                                       

          ! Decide which Y-region the particle is in and then use  
          ! the index to access the X-diff data for the 
          ! appropriate region
          !
          ! D.Elder Nov 23 1990
          if (y.le.0.0) then                                          
            if (y.gt.-chalfl) then 
              j = 1
              elseif (y.lt.-c3halfl) then 
                j = 2
              else 
                j = 3
            endif
          else      
            if (y.gt.c3halfl) then 
              j = 1
            elseif (y.lt.chalfl) then 
              j = 2
            else 
              j = 3
            endif
          endif

          cx_start = cx
          if (cioptn.eq.0) then 
            cx = alpha * yycon + cxafs(iqx,j) +                             
     >        sign (cxbfs(iqx,j),cxcfs(iqx,j)-ranv(kk))         
          elseif (cioptn.eq.1) then 
            cx = alpha * yycon + cxafs(iqx,j)                              
            kk = kk + 1
            if (ranv(kk).lt.cxdps(iqx,j)) then 
              kk = kk +1
              cx = cx + sign(cdpstp,cxcfs(iqx,j)-ranv(kk))
            endif
          endif                
     
          ! Add check for X absorption here
          if (xabsorb_opt.ne.0) then 
            call check_x_absorption(cx,y,sputy,ciz,ierr)
            if (ierr.eq.1) then 
            
              ! Particle absorbed - exit tracking loop - X absorption
              ifate = 10
              goto 790 
            endif 
          endif

          ! Add check for X reflection       
          if (xreflection_opt.ne.0) then
            call check_x_reflection(cx,cx_start)
          endif
                                                                              
          ! Do not need the following two lines for the quick standard        
          ! LIM version with no poloidal diffusion. Then kk will only        
          ! be incremented 4 times per loop count, but this is allowed        
          ! for in subsequent calls to surand, which will only replace        
          ! the kk actual random numbers used.                                                                                                              
          if (big) then                                                     
            kk = kk + 1                                                     
            p  = p  + sign (polods(iqx),ranv(kk)-0.5) + svpols(iqx)    

            ! jdemod - check for p reflection if the option is set
            !   this is done in the routine                   
            call check_p_reflection(p)

            absp = abs(p) 
          endif                                                             
                                                                               
          ! Iterate ctemi for temperature change                                                                                                             
          if ((cioptb.ne.4).or.(cioptb.eq.4.and.ix.gt.jx)) then
            if (vary_2d_bound.eq.0) then
              dtemi = dtemi + (dble(ctembsi(ix,iy))-dtemi) *                
     >          dble(cfts(ix,iy,ciz))                          
              ctemi = sngl(dtemi)   
            else
              dtemi = dtemi + (dble(ctembsi_3d(ip2,ix,iy))-dtemi) *                
     >          dble(cfts_4d(ip2,ix,iy,ciz))                          
              ctemi = sngl(dtemi) 
            endif                                         
          endif

          if (alpha.lt.cftcut) then                                         
            cx    = cftcut                                                     
            alpha = 0.0                                                   
            iqx   = ipos (cx, xs, nxs-1)
          endif

          ! Check if ion has penetrated into the DIVIMP grid.
          if (optdp.eq.1) then
            if (tag2(imp).eq.0.0.and.alpha.gt.mark) then

              ! The ion has not entered the DIVIMP grid region yet:
              tag2(imp) = 1.0
              do ii = 1, nbin                 
                if (y.lt.bsbin(ii)) then 
                  nsbin(ii,ciz) = nsbin(ii,ciz) + 1 
                  izbin(ii)     = izbin(ii)     + ciz
                  exit
                endif
              enddo
              tloss(imp) = 0.0 
            elseif (tag2(imp).eq.1.0) then
              if (alpha.lt.mark) then

               ! The ion has left the DIVIMP grid region.  Record the
               ! time since leaving the grid to get an estimate for 
               ! tauFP to input into DIVIMP:
                tloss(imp) = tloss(imp) + qtim
              else            

                ! The ion has left the DIVIMP grid region and re-entered:
                tloss(imp) = 0.0
              endif
            endif
          endif

          ! Add code to check whether the ion has crossed Y=+/-L - if 
          ! the "shear short circuit" option is active and the current 
          ! poloidal position of the particle does not coincide with the 
          ! limiter - i.e. |P| > CPCO then reset P = (-CPCO,CPCO) 
          ! randomly distributed. 
          ! Use oldy and y to determine if the particle cross the chalfl 
          ! boundaries.
          if (shear_short_circuit_opt.eq.1) then 
            if (absp.gt.cpco.and.(
     >        (y.lt.-chalfl.and.oldy.gt.-chalfl).or.
     >        (y.gt.-chalfl.and.oldy.lt.-chalfl).or.
     >        (y.lt.chalfl.and.oldy.gt.chalfl).or.
     >        (y.gt.chalfl.and.oldy.lt.chalfl))) then
     
              kk = kk+1 
              p = cpco * (2.0*ranv(kk) -1.0)
!              write(6,'(a,i10,10g12.5)') 'Shear1:',imp,
!     >          cx,y,chalfl,oldy,p,cpco,absp
              absp = abs(p)
            endif
          endif
                                                                               
          !-------------------------------------------------------------       
          ! ION INBOARD                                                             
          !-------------------------------------------------------------       
          !write(0,*)'x,y,p,pzone = ',cx, y, p, pzones(ip)                                              
          if (cx.ge.0.0) then                                               
            yycon = yy*coni + 1.0                                           
            alpha = cx / yycon                                              
                                                                               
            ! Reflect off X=a if required; Set iqx pointer                    
            ! Ensure iqx pointer never exceeds array bounds                   
            ! (in case of rounding errors etc.)                                                                                                               
            if (alpha.ge.ca) then                                           
              cx    = 2.0 * ca * yycon - cx                                 
              alpha = cx / yycon                                            
              if (cflrxa) then                                              
                cicrxa = cicrxa + sputy                                     
                cisrxa = cisrxa + cist * sputy                              
                citrxa = citrxa + ctemi * sputy                             
                if (cist.lt.cifrxa) cifrxa = cist                           
                  cflrxa = .false.                                            
                endif                                                         
              endif                                                           

              iqx = min (int(alpha*xscali)+1, nqxsi)                          
                
              !if (iqx.lt.0) then
              !   write(0,*) 'IQX < 0:',ca,cx,yycon,alpha,xscali,iqx
              !endif
                      
              ! Boundary condition y>=2L or y<=-2L                              
              ! If qtim is large its possible that adding 2L still              
              ! leaves the particle outside the region of interest:             
              ! check for this and add another 2l if necessary ...              

              tmp_y = y

              ! sazmod - To be honest, I don't think this boundary
              ! checking code nor the one for the outboard side are 
              ! needed since there is already a boundary checking
              ! part repeated before this, but I'll leave this
              ! here just as  note just in case I'm wrong.
              ! If using a 2D boundary then use the newer ip2.
              if (vary_2d_bound.eq.1) then
                call check_y_boundary(ip2,ix,cx,y,oldy,absy,svy,alpha,
     >            ctwol,sputy,ciz,debugl,ierr,vary_2d_bound)
              else
                call check_y_boundary(ip,ix,cx,y,oldy,absy,svy,alpha,
     >            ctwol,sputy,ciz,debugl,ierr,vary_2d_bound)
              endif
               
              if (ierr.eq.1) then 
               
                ! Write some debugging info
                write (string,'(1x,f10.6,f10.5)') oldalp,oldy                       
                write (6,9003) imp,cist,iqx,iqy,ix,iy,                              
     >            cx,alpha,y,p,svy,ctemi,spara,sputy,ip,it,is,string
                   
              elseif (ierr.eq.2) then
               
                ! Particle Y-absorbed - exit tracking loop
                ifate = 11
                goto 790
              endif

              ! If crossed 2L 
              if (y.ne.tmp_y) then 
                if (cfly2l) then                                              
                  cicy2l = cicy2l + sputy                                     
                  if (cist.lt.cify2l) cify2l = cist                           
                  cfly2l = .false.                                            
                endif                                                         
              endif

              ! sazmod - There was a huge block of comments and 
              ! commented out code here. I made the decision to remove
              ! it all since it presumably does not apply to the 
              ! modern usage of 3DLIM and took up a lot of space. Check 
              ! repository if it actually is relevant. 12/16/21.

              ! sazmod - We're inboard, so shouldn't pzone be set to 1 
              ! (i.e. we're not in CP region)?
              pz = 1
              !pz = pzones(ip)

              ! Force balance with simple collector probe model or no 
              ! collector probe 
              if (colprobe3d.eq.0) then 
                svg = 0.0
                fvel = svy
                if (vary_2d_bound.eq.1) then
                  quant = -seyins(iqx,ciz) - cfss_4d(ip2,ix,iy,ciz) * 
     >              (svy - svhins(iqx)) 
                else
                  quant = -seyins(iqx,ciz) - cfss(ix,iy,ciz) * 
     >              (svy - svhins(iqx))   
                endif          

              elseif (colprobe3d.eq.1) then
            
                ! The 3D collector probe plasma conditions inboard are 
                ! not typical core plasma conditions and so efields and 
                ! gradients are present the fixed efield option for core 
                ! is ignored and inboard flow is added to any local 
                ! plasma velocity. Note: No differences between Y>0, Y<0 
                ! ... need to be careful when spatially varying inboard 
                ! plasmas are used.
              
                if (vary_2d_bound.eq.1) then
                  ff = cfss_4d(ip2,ix,iy,ciz) * (cfvhxs_3d(ip2,ix,iy)
     >              * velplasma_4d(ip2,ix,iy,pz) - svy)
                  fe = cfexzs_4d(ip2,ix,iy,ciz) * 
     >              efield_4d(ip2,ix,iy,pz)
                  fvh = cfvhxs_3d(ip2,ix,iy) * 
     >              velplasma_4d(ip2,ix,iy,pz)
                  fvel = svy

                  ! Record temperature gradient forces
                  feg = calphe(ciz) * ctegs_3d(ip2,ix,iy)
                  fig = cbetai(ciz) * ctigs_3d(ip2,ix,iy)
                  svg = feg+fig
                
                  ! Normal routine.  
                else
                  ff = cfss(ix,iy,ciz) * (cfvhxs(ix,iy)
     >              * velplasma(ix,iy,pz) - svy)
                  fe = cfexzs(ix,iy,ciz) * efield(ix,iy,pz)
                  fvh = cfvhxs(ix,iy) * velplasma(ix,iy,pz)
                  fvel = svy

                  ! Record temperature gradient forces 
                  feg = calphe(ciz) * ctegs(ix,iy)
                  fig = cbetai(ciz) * ctigs(ix,iy)
                  svg = feg+fig
                endif
               
              ! quant is the displacement due to the acceleration caused 
              ! by forces. d = vt + 1/2 a * t^2 essentially.
              quant = ff + fe + feg + fig

            endif
                                                        
            svy = svy + quant                                               
                                                                               
          !-------------------------------------------------------------       
          ! Ion outboard                                                            
          !-------------------------------------------------------------                                                                                       
          else                                                                

            ! jdemod - if poloidal extent limiters are in use in the SOL 
            ! then need to check the +/-2L boundaries which is not 
            ! normally needed in the SOL
            if (big.and.cioptj.eq.1.and.absp.gt.cpco) then 
                  
              ! If using 2D bound then use the newer ip2.
              if (vary_2d_bound.eq.1) then
                call check_y_boundary(ip2,ix,cx,y,oldy,absy,svy,alpha,
     >            ctwol,sputy,ciz,debugl,ierr,vary_2d_bound)
              else
                call check_y_boundary(ip,ix,cx,y,oldy,absy,svy,alpha,
     >            ctwol,sputy,ciz,debugl,ierr,vary_2d_bound)
              endif
                                
              if (ierr.eq.1) then 
                ! write some debugging info
                write (string,'(1x,f10.6,f10.5)') oldalp,oldy                       
                write (6,9003) imp,cist,iqx,iqy,ix,iy,cx,alpha,y,p,
     >            svy,ctemi,spara,sputy,ip,it,is,string 
                   
              elseif (ierr.eq.2) then
                ! Particle Y-absorbed - exit tracking loop
                ifate = 11
                goto 790
              endif
            endif

            yycon = yy * cono + 1.0                                             
            alpha = cx / yycon                                                
                                                                              
            ! Update iqx value.                                                                                                                               
            iqx  = max(int(alpha * xscalo), 1 - nqxso)                           

            ! Calculate an IQY_TMP value to access CAW_QYS which gives 
            ! the wall distance at a specific value of QYS - this 
            ! IQY_TMP has a different meaning than the IQY calculated 
            ! below (which is relative to the limiter faces). 
            if (y.lt.0.0) then 
              iqy_tmp = max(min(int((y + ctwol) * yscale) + 1, nqys), 1)
            else
              iqy_tmp = max(min(int(y * yscale) + 1, nqys), 1)
            endif
                                                                              
            ! Note 151,270: Stopping cross field transport.                    
            ! Originally, when x reached -4 lambda it was brought back          
            ! to 0. For flexibility, a cutoff is now specified in              
            ! the input data - can be switched off by setting to -99.0          
            if (alpha.lt.cftcut) then                                         
              cx    = 0.0                                                     
              alpha = 0.0                                                     
              iqx   = 0                                                       
C                                                                               
C------------ ION HAS REACHED WALL AND IS ABSORBED                              
C------------ SET TEMPERATURE TO SMALLEST POSSIBLE VALUE                        
C------------ SCORE PARTICLE IN "WALLS" ARRAY (REDIM'D TO -NYS:NYS)             
C                                                                               
c              ELSEIF (ALPHA.LE.CAW) THEN                                        
c
              ELSEIF (ALPHA.LE.CAW_QYS(IQY_TMP)) THEN                                        
                CIAB   = 0                                                      
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(1,IY)                                           
                IF (Y.GT.0.0) THEN                                              
                  CALL MONUP (6,SPUTY)                                          
                  RWALL(2) = RWALL(2) + SPUTY                                   
                ELSE                                                            
                  CALL MONUP (10,SPUTY)                                         
                  RWALL(1) = RWALL(1) + SPUTY                                   
                ENDIF                                                           
                WALLS(IY,CIZ) = WALLS(IY,CIZ) + SPUTY                           
                IFATE = 1                                                       
c slmod begin
c
c Some statistics for the DIVIMP grid source stuff:
c
                IF (optdp.EQ.1) THEN
                  IF (TLOSS(IMP).NE.0.0.AND.IFATE.EQ.1) THEN
                    ALOSS = ALOSS + TLOSS(IMP)
                    WLOSS = WLOSS + 1.0
                  ENDIF
                ENDIF
c slmod end
                GOTO 790                                                        
              ENDIF                                                             

C
C     IF POLOIDAL EXTENT LIMITS ARE IN EFFECT ONLY CHECK FOR COLLISION
C     IF THE P COORDINATE IS INSIDE +/- CPCO.  
C     NOTE: A MAJOR ASSUMPTION IN ALMOST ALL OF THIS CODE IS THAT STEP
C           SIZES ARE SMALL. THUS COLLISIONS OCCUR AT THE APPROXIMATE 
C           POSITIONS OF THE PARTICLES WHEN THEY ARE FOUND TO HAVE HIT 
C           THE LIMITER. FIRST ORDER IS THEN TO RELEASE NEWLY SPUTTERED
C           PARTICLES FROM THE LAST X,Y,P VALUES. SOME REFINING OF THE 
C           X,Y VALUES IS DONE IN HIT AND EDGINT. FOR NOW I WILL LEAVE
C           THE P VALUES UNREFINED.

          ! sazmod - Don't check for limiter collision when in 1DLIM mode.
          IF (((.NOT.BIG).OR.(.NOT.((CIOPTJ.EQ.1).AND.
     >          (ABSP.GT.CPCO)))).and.lim1d.eq.0) THEN

C                                                                               
C------------ ION HAS HIT LIMITER AT Y=0 FROM Y>0 REGION                        
C                                                                               

              EDGE1 = QEDGES(IQX,1)                                             
              EDGE2 = QEDGES(IQX,2)                                             
              IF (OLDY.GT.0.0 .AND. Y.LE.EDGE2) THEN                            
                CIAB   = 1                                                      
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (6,SPUTY)                                            
                IFATE = 2                                                       
                GOTO 780                                                        
C                                                                               
C------------ ION HAS HIT LIMITER AT Y=2L                                       
C                                                                               
              ELSEIF (Y .GE. CTWOL-EDGE1) THEN                                  
                CIAB   = 2                                                      
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (10,SPUTY)                                           
                IFATE = 3                                                       
                GOTO 780                                                        
C                                                                               
C------------ ION HAS HIT LIMITER AT Y=0 FROM Y<0 REGION                        
C                                                                               
              ELSEIF (OLDY.LT.0.0 .AND. Y.GE.-EDGE1) THEN                       
                CIAB   = -1                                                     
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (10,SPUTY)                                           
                IFATE = 4                                                       
                GOTO 780                                                        
C                                                                               
C------------ ION HAS HIT LIMITER AT Y=-2L                                      
C                                                                               
              ELSEIF (Y .LE. EDGE2-CTWOL) THEN                                  
                CIAB   = -2                                                     
                CALL HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)        
                CVABS  = SVY / QFACT                                            
                CTBIQX = CTEMBS(IX,IY)                                          
                CALL MONUP (6,SPUTY)                                            
                IFATE = 5                                                       
                GOTO 780                                                        
              ENDIF                                                             

C           MATCHING ENDIF FOR POLOIDAL EXTENT TEST
            ENDIF  

c slmod tmp
c              WRITE(0,*) 'Error! Ion outboard: ',IMP,CX,Y,ABSP
c              STOP
c slmod end

C                                                                               
C-------------- SET IQY VALUE AND                                               
C-------------- UPDATE ION VELOCITY, AFFECTED BY                                
C-------------- ELECTRIC FIELD AND DRIFT VELOCITY                               
C-------------- ION HAS COMPLETED OUTBOARD STEP   (MONUP1)                      
C                                                                               
c
c           jdemod - WARNING - there is an inconsistency in the 
c                    definition of IQY in the code. IQY is the
c                    index into the underlying QYS and related
c                    arrays. However, it appears that in the code
c                    below IQY was at some point redefined to 
c                    be a number of points between the limiter
c                    surfaces along the field lines - thus the
c                    use of EDGE1, EDGE2 and CYSCLS in the
c                    calculation of IQY here. However, the related
c                    variables that are indexed by IQY were NOT
c                    changed to reflect this usage. 
c                    e.g. CEYS(IQY), CVHYS(IQY) and QYS(IQY) are
c                         all calculated without taking the 
c                         location of the limiter edges into effect.
c                    This means that the IQY value calculated here
c                    does not match properly with these arrays - in
c                    particular QYS. 
c
c                    On the other hand, the 
c                    procedure here will map the IQY range between
c                    limiter surfaces 1:1 onto the CEYS and CVHQYS
c                    arrays thus allowing for a variable limiter
c                    shape still mapping the first data point at IQY=1
c                    to the limiter surface.
c                    The dependencies in CEYS and CVHYS are typically
c                    linear in QYS to CL and so would not change 
c                    much if they were calculated properly for the 
c                    actual limiter surface location. 
c
c                    The biggest concern is QYS - it can not be
c                    properly indexed by this IQY.                     
c
c
c            SVG = CALPHE(CIZ) * CTEGS(IX,IY) +
c     >            CBETAI(CIZ) * CTIGS(IX,IY) 
c
c
c           jdemod = - record temperature gradient forces
c
            if (vary_2d_bound.eq.0) then
              feg = calphe(ciz) * ctegs(ix,iy)
              fig = cbetai(ciz) * ctigs(ix,iy)
c            SVG = CALPHE(CIZ) * CTEGS(IX,IY) +
c     >            CBETAI(CIZ) * CTIGS(IX,IY)              
            else
              feg = calphe(ciz) * ctegs_3d(ip2,ix,iy)
              fig = cbetai(ciz) * ctigs_3d(ip2,ix,iy)
            endif
            svg = feg + fig

c           jdemod
c
c           Add frictional coupling to parallel flow beyond the limiter
c           extent if one is specified. (vpflow_3d (L28) - default is 0.0)
c           Only in SOL.
c
c
c           Outboard parallel force balance
c            

            ! jdemod - this is the original LIM calculation for IQY 
            ! gives the incorrect index for Y< 0
            ! CVHYS is defined for Y=0 to Y = 2L 
            ! Y = -2L to 0 should map onto the range [0,2L] 1:1
            ! so Y = -2L -> Y= 0
            ! However the IQY calculation below maps Y= -2L to an index for Y = +2L
            ! which inverts the velocity array for Y<0
            ! This bug was compensated for in the transport equation by using the
            ! opposite sign for the frictions force in Y<0. For simple or 
            ! symmetric velocity profiles this isn't an issue but for 
            ! more complex or assymmetric profiles it is a problem - so I am fixing
            ! the indexing and changing the sign of the frictional force in Y<0
            
            if (y.gt.0.0) then 
              IQY   = INT ((Y-EDGE2) * CYSCLS(IQX)) + 1                    
            else
              IQY   = INT((CTWOL+Y-EDGE1) * CYSCLS(IQX)) + 1                     
              !IQY   = INT((-Y-EDGE1) * CYSCLS(IQX)) + 1                     
            endif


            ! set pz = 1 for now  (1 = not cp region, 2 = cp region)
            !pz = pzones(ip)
            pz = 1

            ! sazmod - Moving this into the colprobe3d.eq.0 block.
            ! 12/16/21 - Undoing this, shouldn't matter though.
            ! 1/18/22 - Moving back into block. Leads to segmentation
            !   fault if left here sometimes.
            !if (vel_efield_opt.eq.0) then
            !   efield_val = CEYS(IQY)
            !   velplasma_val = CVHYS(IQY)
            !elseif (vel_efield_opt.eq.1) then 
            !   efield_val = efield(ix,iy,pz)
            !   velplasma_val = velplasma(ix,iy,pz)
            !endif

            ! Force balance with no collector probe. 
            if (colprobe3d.eq.0) then 
            
            ! sazmod - Moved here.
            if (vel_efield_opt.eq.0) then
            
			   ! sazmod - I get a segfault here, I think the above code
			   ! for iqy isn't quite right. Instead of dealing with it,
			   ! I'm just going to jump around it for SOL option 0.
			   if (cioptf.eq.0) then
			     efield_val = 0.0
			     velplasma_val = 0.0
			   else
                 efield_val = CEYS(IQY)
                 velplasma_val = CVHYS(IQY)
               endif
            elseif (vel_efield_opt.eq.1) then 
               efield_val = efield(ix,iy,pz)
               velplasma_val = velplasma(ix,iy,pz)
            endif
            
            IF (Y.GT.0.0) THEN                                              
              !IQY   = INT ((Y-EDGE2)  * CYSCLS(IQX)) + 1                    
              IF ((BIG).AND.(CIOPTJ.EQ.1).AND.(ABSP.GT.CPCO)) THEN
                ! jdemod - assign forces
                fe = 0.0
                ff = -CFSS(IX,IY,CIZ)*(SVY-vpflow_3d)
                fvel = svy
                fvh = vpflow_3d
c                QUANT = -CFSS(IX,IY,CIZ)*(SVY-vpflow_3d)
                quant = ff
             ELSE 
                !     jdemod - assign forces
                ff   = (CFSS(IX,IY,CIZ)*
     >                  (CFVHXS(IX,IY)*velplasma_val-SVY))
                fe   = (CFEXZS(IX,IY,CIZ) * efield_val)
                fvh  = CFVHXS(IX,IY)*velplasma_val
                fvel = svy
                
c                QUANT = (CFEXZS(IX,IY,CIZ) * CEYS(IQY)) + SVG +               
c     >           (CFSS(IX,IY,CIZ)*(CFVHXS(IX,IY)*CVHYS(IQY)-SVY))  
                quant = fe + svg + ff

             ENDIF 
            ELSE                                                            
              !IQY   = INT((-Y-EDGE1) * CYSCLS(IQX)) + 1                     
              IF ((BIG).AND.(CIOPTJ.EQ.1).AND.(ABSP.GT.CPCO)) THEN
                ! jdemod - assign forces
                fe = 0.0
                ff = -CFSS(IX,IY,CIZ)*(SVY-vpflow_3d)
                fvel = svy
                fvh = vpflow_3d
c                QUANT = -CFSS(IX,IY,CIZ)*(SVY-vpflow_3d)
                quant = ff
             ELSE
                ! jdemod - assign forces
                ff   = (CFSS(IX,IY,CIZ)*
     >                  (CFVHXS(IX,IY)*velplasma_val-SVY))
                fe   = (CFEXZS(IX,IY,CIZ) * efield_val)
                fvh  = CFVHXS(IX,IY)*velplasma_val
                fvel = svy
c                QUANT =-(CFEXZS(IX,IY,CIZ) * CEYS(IQY)) + SVG -              
c     >           (CFSS(IX,IY,CIZ)*(CFVHXS(IX,IY)*CVHYS(IQY)+SVY))      
!                jdemod - note sign change on ff to account for summation in quant
                quant = fe + svg + ff 
             ENDIF
            ENDIF                                                           
            
          ! Force balance with collector probe plasma.
          elseif (colprobe3d.eq.1) then 

            ! Determine if on a flux tube connected to probe
            ! since this affects the Efield and friction forces.
            pz = pzones(ip)
             
               ! use forces for areas not connected to a probe

!                QUANT = (CFEXZS(IX,IY,CIZ) * CEYS(IQY)) + SVG +               
!     >           (CFSS(IX,IY,CIZ)*(CFVHXS(IX,IY)*CVHYS(IQY)-SVY))  

            if (vary_2d_bound.eq.1) then
            
              ! Assign forces.
              ff   = (cfss_4d(ip2,ix,iy,ciz) * (cfvhxs_3d(ip2,ix,iy)
     >          * velplasma_4d(ip2,ix,iy,pz) - svy))
              fe = (cfexzs_4d(ip2,ix,iy,ciz) * efield_4d(ip2,ix,iy,pz))
              fvh = cfvhxs_3d(ip2,ix,iy) * velplasma_4d(ip2,ix,iy,pz)
              fvel = svy
            else
           
              ! Assign forces.
              ff   = (cfss(ix,iy,ciz) * (cfvhxs(ix,iy)
     >          * velplasma(ix,iy,pz) - svy))
              fe   = (cfexzs(ix,iy,ciz) * efield(ix,iy,pz))
              fvh  = cfvhxs(ix,iy) * velplasma(ix,iy,pz)
              fvel = svy
            endif
            quant = ff + fe + svg
                
          endif

             
            SVY = SVY + QUANT                                               

            DOUTS(CIZ,10) = DOUTS(CIZ,10) + DSPUTY * DQFACT                       

            ! jdemod - record some force statistics
            DOUTS(CIZ,1) = DOUTS(CIZ,1) + DSPUTY
            
            ! The numbers here are not really too reasonable with an
            ! actual 2D boundary, but can at least be used in a 
            ! debugging sense when using a "varying" boundary of all
            ! the same value.
            if (vary_2d_bound.eq.1) then
              douts(ciz,2) = douts(ciz,2) + dsputy * dtemi / 
     >          cfps_4d(ip2,ix,iy,ciz)
              douts(ciz,3) = douts(ciz,3) + dsputy / 
     >          cfss_4d(ip2,ix,iy,ciz)
     
            else
              DOUTS(CIZ,2) = DOUTS(CIZ,2)+DSPUTY * DTEMI/CFPS(IX,IY,CIZ)
              DOUTS(CIZ,3) = DOUTS(CIZ,3)+DSPUTY / CFSS(IX,IY,CIZ)
            endif
            DOUTS(CIZ,4) = DOUTS(CIZ,4) + DSPUTY * abs(FF)
            DOUTS(CIZ,5) = DOUTS(CIZ,5) + DSPUTY * abs(FE)
            DOUTS(CIZ,6) = DOUTS(CIZ,6) + DSPUTY * abs(FEG)
            DOUTS(CIZ,7) = DOUTS(CIZ,7) + DSPUTY * abs(FIG)
            DOUTS(CIZ,8) = DOUTS(CIZ,8) + DSPUTY * abs(FVEL)
            DOUTS(CIZ,9) = DOUTS(CIZ,9) + DSPUTY * abs(FVH)

            IF (ALPHA.LT.CRXMIN) CRXMIN = ALPHA                             

          ENDIF                                                             


C     
C-----------------------------------------------------------------------        
C    BOTH ROUTES CONTINUE HERE.     CHECK FOR COLLISION                         
C    THIS SECTION ACCOUNTS FOR 6% OF CPU TIME AND COULD BE COMMENTED OUT        
C-----------------------------------------------------------------------        
C                                                                               
              KK = KK + 1     
              if (vary_2d_bound.eq.0) then                                                  
                IF ((CTEMI*RANV(KK)) .LE. CFPS(IX,IY,CIZ)) THEN                   
                  CICCOL = CICCOL + SPUTY * QFACT  
                endif
              else
                 IF ((CTEMI*RANV(KK)) .LE. CFPS_4d(ip2,IX,IY,CIZ)) THEN                   
                  CICCOL = CICCOL + SPUTY * QFACT   
                 endif                      
              ENDIF                                                             
C                                                                               
C-----------------------------------------------------------------------        
C    FIND USER X,Y,P BINS ION LIES WITHIN                                       
C    DO NOT NEED LINES  "430 CONTINUE" THROUGH TO "GOTO 440" FOR                
C    THE QUICKER STANDARD LIM VERSION WITH NO POLOIDAL DIFFUSION.               
C-----------------------------------------------------------------------        
C                                                                               
c
c             jdemod - change structure of these statements to address
c                      intel fortran compiler issue
c
              IF (.NOT.BIG) GOTO 450                                            
  430         CONTINUE                                                          
                IF ((IP.LE.-MAXNPS))GOTO 440                 
                IF ((PS(IP-1).LT.P))GOTO 440                 
                IP = IP - 1                                                     
                GOTO 430                                                        
  440         CONTINUE                                                          
                IF ((IP.GE.MAXNPS)) GOTO 450                 
                IF ((PS(IP).GE.P)) GOTO 450                 
                IP = IP + 1                                                     
                GOTO 440                                                        

  450         CONTINUE                                                          
                IF ((JY.LE.1)) GOTO 460                 
                IF ((YS(JY-1).LT.ABSY)) GOTO 460                 
                JY = JY - 1                                                     
                GOTO 450                                                        
  460         CONTINUE                                                          
                IF ((JY.GE.NYS)) GOTO 470                 
                IF ((YS(JY).GE.ABSY)) GOTO 470                 
                JY = JY + 1                                                     
                GOTO 460                                                        

  470         CONTINUE                                                          
                IF ((IX.LE.1)) GOTO 480                 
                IF ((XS(IX-1).LT.ALPHA)) GOTO 480                 
                IX = IX - 1                                                     
                GOTO 470                                                        
  480         CONTINUE                                                          
                IF ((IX.GE.NXS)) GOTO 490                 
                IF ((XS(IX).GE.ALPHA)) GOTO 490                 
                IX = IX + 1                                                     
                GOTO 480                                                        
  490         CONTINUE                                                          
              IY = JY    
                                                                     
              ! jdemod - I'm not sure how the IY=-IY can be correct
              ! since this will mirror the particle
              ! location in terms of Y. 
              ! sazmod - A ways up we defined jy as abs(iy). Since we 
              ! reassigned iy = jy, we need to make it negative again if 
              ! necessary.
              if (y.lt.0.0) iy = -iy                                            
C                                                                               
C-----------------------------------------------------------------------        
C    SCORE PARTICLE IN DDLIMS "NUMBER DENSITY" ARRAY AND IN THE                 
C    CLOUD TEMPERATURES ARRAY DDTS (IN DOUBLE PRECISION)                        
C    ALSO STORE LIM5 COMPONENT, IE TIME DEPENDENT DISTRIBUTIONS, IF             
C    WE HAVE REACHED A GIVEN TIMEPOINT.                                         
C    SOME SECTIONS BELOW ARE ONLY NEEDED FOR 3D AND TIME DEPENDENT              
C    CASES - THEY ARE NOT REQUIRED FOR THE STANDARD VERSION.                    
C    SCORE IN Y STEPSIZES ARRAY THE PARALLEL DIFFUSION COEFFICIENT.             
C-----------------------------------------------------------------------        
C                                                                  

             ! Further comment helping explain sputy and the ddlim and
             ! ddlim3 array. These number density arrays at each 
             ! location contain the sum of the "number of particles"   
             ! at each location during each time step. Some further 
             ! scaling later on I think effectively normalizes this
             ! to some kind of unit of particles/m3 (m2 for ddlim),
             ! maybe some other units included. But this is taken care
             ! of in sputtering cases by absfac, which is calculated
             ! for you.    
             DDLIMS(IX,IY,CIZ) = DDLIMS(IX,IY,CIZ) + DSPUTY * DQFACT           

c            Update velocity diagnostics            
             call update_diagvel(ix,iy,ciz,dsputy*dqfact,dble(svy))

             ! sazmod - changing to le.maxnps. Not sure why ip couldn't
             ! equal the maxnps.
             if (jy.le.ny3d.and.abs(ip).le.maxnps) then
               ddlim3(ix,iy,ciz,ip) = ddlim3(ix,iy,ciz,ip) + 
     >           dsputy * dqfact
            else
              write(0,*) 'warning: ip > maxnps or jy > ny3d :',
     >          ip,maxnps,jy,ny3d
             endif
c
c              IF (JY.LE.NY3D) DDLIM3(IX,IY,CIZ,IP) =   
c     >           DDLIM3(IX,IY,CIZ,IP) + DSPUTY * DQFACT
c slmod end
c
c slmod begin
c                IF (DEBUGL) WRITE(79,*)
c     +            CX,ALPHA,Y,P,IX,IY,IP,CIZ,DDLIM3(IX,IY,CIZ,IP)
c slmod end
c              if (debugl) then
c                 write(6,'(a,8i8,3(1x,g12.5))')
c     >                'LIM5:',jy,ny3d,ix,iy,ciz,ip,it,cdwelt_sum,
c     >                    cist,ctimes(it,ciz),
c     >                        lim5(ix,iy,iz,ip,it)
c              endif
c             
c             jdemod - time relative to t=0 for simulation is used for assigning
c                      the time bins - only update for imode not equal to 2 - i.e. time dependent
c     
              IF (IMODE.ne.2.and.RTIME.GE.CTIMES(IT,CIZ)) THEN                                  
c              IF (CIST.GE.CTIMES(IT,CIZ)) THEN                                  
c
c     IF (DEBUGL) WRITE (6,9003) IMP,CIST+QFACT,IQX,IQY,IX,IY,                
c     >    CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,'UPDATE LIM5'
c
c
c     jdemod - remove cdwelt_sum option functionality because it isn't
c              physically meaningful.                  
c                if (cdwelt_sum.eq.0) then 

                    IF (JY.LE.NY3D)                                                 
     >                LIM5(IX,IY,CIZ,IP,IT) = LIM5(IX,IY,CIZ,IP,IT)
     >                  + SPUTY 
c                endif
                ! jdemod - the update to the time bin has to be AFTER the
                ! particle has been recorded!!
                IT = IT + 1                                                     

             ENDIF                                                             

              ! jdemod - move this outside the test for time bin for cdwelt_sum option 1
              !        - otherwise  data is recorded in this array only once
              !        - does this need to be double precision?
c              if (cdwelt_sum.eq.1) then 
c                 IF (JY.LE.NY3D)                                                 
c     >             LIM5(IX,IY,CIZ,IP,IT) = LIM5(IX,IY,CIZ,IP,IT)
c     >                  + SPUTY * QFACT        
c              endif
          
              DDTS(IX,IY,CIZ) =DDTS(IX,IY,CIZ)+DSPUTY*   DTEMI   *DQFACT        
              DDYS(IX,IY,CIZ) =DDYS(IX,IY,CIZ)+DSPUTY*DBLE(SPARA)*DQFACT        
C                                                                               
C-----------------------------------------------------------------------        
C    TEST IF IONISATION OR RECOMBINATION OCCURING                               
C    THE SECOND RANDOM NUMBER IS NOT USED IN EVERY ITERATION                    
C    THROUGH LOOP 500, HENCE EXTRA CALL USED FOR SURAND.                        
C-----------------------------------------------------------------------        
C                                                                               
              KK = KK + 1      
              if (vary_2d_bound.eq.0) then                                                 
              IF (RANV(KK).LE.CPCHS(IX,IY,CIZ)
     >            .and.ranv(kk).gt.0.0) THEN                            
                KK = KK + 1                                                     
c                write(6,*) 'CHS:',ix,iy,ciz,kk,
c     >               ranv(kk-1),cpchs(ix,iy,ciz),
c     >               ranv(kk),cprcs(ix,iy,ciz)
                IF (RANV(KK).LE.CPRCS(IX,IY,CIZ)
     >              .and.ranv(kk).gt.0.0) THEN                          
C                                                                               
C---------------- EVENT IS A RECOMBINATION.  UPDATE MONITOR VARS                
C---------------- POINTERS ETC.  CHECK FOR C+ --> C EVENT WHICH                 
C---------------- MEANS WE GO ON TO THE NEXT ION.                               
C                                                                               
                  CICRCS(CIZ) = CICRCS(CIZ) + SPUTY                             
                  IF (CIST .LT. CIFRCS(CIZ)) CIFRCS(CIZ) = CIST                 
                  IF (CIST .GT. CILRCS(CIZ)) CILRCS(CIZ) = CIST                 
                  CISRCS(CIZ) = CISRCS(CIZ) + CIST * SPUTY                      
                  CIZ  = CIZ - 1                                                
c
c                 jdemod - RTIME is particle time since t=0 for simulation
c                  
                  IF (BIG) IT = IPOS (RTIME, CTIMES(1,CIZ), NTS)                 
c                  IF (BIG) IT = IPOS (CIST, CTIMES(1,CIZ), NTS)                 
c slmark
                  IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,            
     >              CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,              
     >              'RECOMBINED:',CIZ                                           
                  IF (CIZ.LT.1) THEN                                            
                    TBELOW = TBELOW + SPUTY                                     
                    IFATE = 9                                                   
                    GOTO 790                                                    
                  ENDIF                                                         
                ELSE                                                            
C                                                                               
C---------------- EVENT IS AN IONISATION.    UPDATE MONITOR VARS                
C---------------- POINTERS ETC.  CHECK FOR GOING ABOVE MAXIMUM                  
C---------------- IONISATION STATE - GO ON TO THE NEXT ION.                     
C---------------- RECORD IONISATION POSITION IN "TIZS" ARRAY                    
C---------------- IF SET TI>=TB APPLIES, DO THAT AS WELL.                       
C                                                                               
                  CICIZS(CIZ) = CICIZS(CIZ) + SPUTY                             
                  IF (CIST .LT. CIFIZS(CIZ)) CIFIZS(CIZ) = CIST                 
                  IF (CIST .GT. CILIZS(CIZ)) CILIZS(CIZ) = CIST                 
                  CISIZS(CIZ) = CISIZS(CIZ) + CIST * SPUTY                      
                  TIZS(IX,IY,CIZ) = TIZS(IX,IY,CIZ) + SPUTY                     
                  IF (JY.LE.NY3D)                                               
     >              TIZ3(IX,IY,CIZ,IP) = TIZ3(IX,IY,CIZ,IP) + SPUTY             
                  CIZ  = CIZ + 1                                                
c
c                 jdemod - RTIME is particle time since t=0 for simulation
c                  
                  IF (BIG) IT = IPOS (RTIME, CTIMES(1,CIZ), NTS)                 
c                  IF (BIG) IT = IPOS (CIST, CTIMES(1,CIZ), NTS)                 
                  IF (CIZ.EQ.CIZSET) CTEMI = MAX (CTEMI,CTEMBS(IX,IY))          
                  IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,            
     >              CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,             
     >              'IONISED TO:',CIZ                                           
                  IF (CIZ .GT. NIZS) THEN                                       
                    TBYOND = TBYOND + SPUTY                                     
                    IFATE = 9                                                   
c slmod begin
c
c Check to see if the ion ionised beyond the ionisation state limit before
c it entered the DIVIMP grid.  If this happens a lot, then the validity of the
c ion source for input into DIVIMP would be questionable:
c
                    IF (optdp.EQ.1) THEN
                      IF (TAG2(IMP).NE.1.0) TSLOSS = TSLOSS + 1
                      IZLOSS = IZLOSS + 1.0
                    ENDIF
c slmod end
                    GOTO 790                                                    
                  ENDIF                                                         
                  MAXCIZ = MAX (MAXCIZ, CIZ)                                    
                ENDIF                                                           
              ENDIF      
              
              ! vary_2d_bound routine. Just copy/pasted from above with
              ! respective 3D/4D arrays swapped in. Fixed indentation.
              else
                if (ranv(kk).le.cpchs_4d(ip2,ix,iy,ciz)
     >            .and.ranv(kk).gt.0.0) then                            
                  kk = kk + 1                                                     
                  if (ranv(kk).le.cprcs_4d(ip2,ix,iy,ciz)
     >              .and.ranv(kk).gt.0.0) then                          
                  
                    ! Recombination.                                                        
                    cicrcs(ciz) = cicrcs(ciz) + sputy                             
                    if (cist .lt. cifrcs(ciz)) cifrcs(ciz) = cist                 
                    if (cist .gt. cilrcs(ciz)) cilrcs(ciz) = cist                 
                    cisrcs(ciz) = cisrcs(ciz) + cist * sputy                      
                    ciz = ciz - 1                                                
                  
                    if (big) it = ipos (rtime, ctimes(1,ciz), nts)                 

                    if (debugl) write (6,9003) imp,cist,iqx,iqy,ix,iy,            
     >                cx,alpha,y,p,svy,ctemi,spara,sputy,ip,it,is,              
     >                'recombined:',ciz                                           
                    if (ciz.lt.1) then                                            
                      tbelow = tbelow + sputy                                     
                      ifate = 9                                                   
                      goto 790                                                    
                    endif                                                         
                  else                                                            
                  
                    ! Ionization.                                                  
                    cicizs(ciz) = cicizs(ciz) + sputy                             
                    if (cist .lt. cifizs(ciz)) cifizs(ciz) = cist                 
                    if (cist .gt. cilizs(ciz)) cilizs(ciz) = cist                 
                    cisizs(ciz) = cisizs(ciz) + cist * sputy                      
                    tizs(ix,iy,ciz) = tizs(ix,iy,ciz) + sputy
                    
                    ! Similar to DDLIM3, this array is already 4D and
                    ! ready for use in the vary_2d_bound routines, but
                    ! the original ip is used to index it.                        
                    if (jy.le.ny3d)                          
     >                tiz3(ix,iy,ciz,ip) = tiz3(ix,iy,ciz,ip) + sputy             
                    ciz  = ciz + 1                                                

                    if (big) it = ipos (rtime, ctimes(1,ciz), nts)                 
                    if (ciz.eq.cizset) then
                      ctemi = max(ctemi,ctembs_3d(ip2,ix,iy))    
                    endif      
                    if (debugl) write (6,9003) imp,cist,iqx,iqy,ix,iy,            
     >                cx,alpha,y,p,svy,ctemi,spara,sputy,ip,it,is,             
     >                'ionised to:',ciz                                           
                    if (ciz .gt. nizs) then                                       
                      tbyond = tbyond + sputy                                     
                      ifate = 9                                                   

                      if (optdp.eq.1) then
                        if (tag2(imp).ne.1.0) tsloss = tsloss + 1
                        izloss = izloss + 1.0
                      endif
                      goto 790                                                    
                    endif                                                         
                    maxciz = max (maxciz, ciz)                                    
                  endif                                                           
                endif  
              endif                                                       
C                                                                               
C-----------------------------------------------------------------------        
C             SPLITTING PLANE CROSSED - SAVE ION DETAILS, LEAP TO 790           
C-----------------------------------------------------------------------        
C                                                                               
              IF (split.and.ALPHA.GT.CXSPLS(IS)) THEN                                     
                IF (IPUT.EQ.MAXPUT) THEN                                        
                  WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,CX,ALPHA,Y,             
     >              P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,'SPLIT FAILED'            
                ELSE                                                            
                  TSPLIT(IS) = TSPLIT(IS) + SPUTY                               
                  NSPLIT(IS) = NSPLIT(IS) + 1                                   
                  IS = IS + 1                                                   
                  IPUT = IPUT + 1                                               
                  MPUT = MAX (MPUT, IPUT)                                       
                  IGET(IPUT) = 1                                                
                  R(1,IPUT) = CX                                                
                  R(2,IPUT) = Y                                                 
                  R(3,IPUT) = P                                                 
                  R(4,IPUT) = QUANT                                             
                  R(5,IPUT) = SVY                                               
                  R(6,IPUT) = CTEMI                                             
                  R(7,IPUT) = CIST                                              
                  R(8,IPUT) = SPUTY                                             
                  R(9,IPUT) = TSTEPL                                            
                  R(10,IPUT)= RTIME   ! Add RTIME to splitting/rouletting data recorded
                  I(1,IPUT) = IQX                                               
                  I(2,IPUT) = IQY                                               
                  I(3,IPUT) = IX                                                
                  I(4,IPUT) = IY                                                
                  I(5,IPUT) = IP                                                
                  I(6,IPUT) = IT                                                
                  I(7,IPUT) = CIZ                                               
                  I(8,IPUT) = MAXCIZ                                            
                  I(9,IPUT) = IS                                                
                  L(1,IPUT) = DIFFUS                                            
                  L(2,IPUT) = CFLRXA                                            
                  L(3,IPUT) = CFLY2L                                            
                  IFATE = 7                                                     
                  GOTO 790                                                      
                ENDIF                                                           
C                                                                               
C-----------------------------------------------------------------------        
C             ROULETTING PLANE CROSSED - DISCARD OR KEEP WITH NEW WEIGHT        
C-----------------------------------------------------------------------        
C                                                                               
              ELSEIF (ALPHA.LT.CXSPLS(IS-1)) THEN                               
                IS = IS - 1                                                     
                KK = KK + 1                                                     
                IF (RANV(KK).GE.CPRUL) THEN                                     
                  IFATE = 8                                                     
                  GOTO 790                                                      
                ELSE                                                            
                  SPUTY  = SPUTY / CPRUL                                        
                  DSPUTY = DBLE (SPUTY)                                         
                  TRULET(IS) = TRULET(IS) + SPUTY                               
                  NRULET(IS) = NRULET(IS) + 1                                   
                  IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,            
     >              CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,               
     >              'ROULETTE KEPT'                                             
                ENDIF                                                           
              ENDIF                                                             
c slmod begin
              IF (optdp.EQ.1.AND.Y.LT.TARGET) THEN
                IF (TLOSS(IMP).NE.0.0) THEN
                  ALOSS  = ALOSS  + TLOSS(IMP)
                  TGLOSS = TGLOSS + 1.0
                ENDIF
c
c Cheating here - not really the wall, may throw off some statistics:
c
                IFATE = 1 
                GOTO 790
              ENDIF
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C       LOOP BACK TO 500 IF ION HAS NOT YET REACHED CUTOFF TIME.                
C       OTHERWISE RECORD THIS FACT WITH MONUP AND STOP ITERATING.               
C       DEBUG PRINTOUT EVERY CSTEPL'TH ITERATION  (EG EVERY 100)                
C       CSTEPL ENTERED AS 0 FOR NO DEBUG OPTION, >0 FOR DEBUG ON.               
C-----------------------------------------------------------------------        
C                                                                               
        IF (DEBUGL) THEN                                                        
          IF (CIST.GE.TSTEPL) THEN                                              
  495       TSTEPL = TSTEPL + CSTEPL                                            
            IF (TSTEPL.LE.CIST) GOTO 495                                        
            WRITE (STRING,'(1X,F10.6,F10.5)') OLDALP,OLDY                       
c slmod
c            WRITE (6,'(I4,F8.3,3I4,5G14.5  )') 
c     +        IMP,CIST,IX,IY,IP,Y,CTEMI,SPARA,CFTS(IX,IY,CIZ),QTIM
            WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,                              
     >        CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,STRING               
c slmod end
          ENDIF                                                                 
        ENDIF  
        
        ! sazmod - 1DLIM usage: Check if ion is killed off due to
        ! parallel sink action.
        if ((lim1d.eq.1).and.(ntausink.gt.0)) then
          call surand(seed,1,ran)
          nrand = nrand + 1
          if (ran.le.qtim/tausink(ix)) then
            !write(0,*) 'tausink: absorbed!'
            goto 790
          endif
        endif
                                                                         
C                                                                               
        DIST   = DIST + DQFACT                                                  
        CIST   = SNGL (DIST)                                                    
c
c       jdemod - update time simce t=0      
c     
        DTIME  = DTIME + DQFACT
        RTIME  = SNGL (DTIME)
c
        QFACT  = QS(IQX)                                                        
        DQFACT = DBLE (QFACT)                                                   
        !IF (CIST.LT.CSTMAX) GOTO 500 
        if (cist.lt.cstmax) then
          goto 500
        else
          write(0,*) 'Warning! Particle cutoff time reached. Consider'//
     >      ' increasing max iteration time.' 
        endif                                           
C                                                                               
        CALL MONUP (7,SPUTY)                                                    
        TCUT = TCUT + SPUTY                                                     
        IFATE = 6                                                               
        GOTO 790                                                                
C                                                                               
C-----------------------------------------------------------------------        
C   ION DEPOSITED ON LIMITER:  STORE DETAILS READY FOR RE-LAUNCH OF             
C   NEUTRAL FRAGMENTS WHEN ALL IONS HAVE BEEN EXHAUSTED.                        
C   EXTRACT EXACT X BIN PARTICLE LIES WITHIN AT THE END OF THIS CURRENT         
C   TIMESTEP  (PRESENT VALUE OF IX INDICATES BIN AT THE START OF THE            
C   CURRENT TIMESTEP, WHICH MIGHT BE FOR X > 0)                                 
C                                                                               
C   BEWARE !!!!! SELF-SPUTTERING COULDN'T BE USED WITH SPLITTING, SINCE         
C   WE WERE OVERWRITING THE ARRAY SPUTYS.  THIS IS ONLY OK SO LONG AS           
C   NPROD <= IMP AT ALL TIMES.  WITH SPLITTING, WE MIGHT HAVE TO RECORD         
C   SAY 3 ABSORPTIONS OF SUB-IONS FOR SOME SPLIT IONS, WHICH WOULD MESS         
C   EVERYTHING UP.  HENCE USE OF SNEWS TEMPORARY ARRAY BELOW.                   
C                                                                               
C   NOTE:  USE OF CX HERE RATHER THAN ALPHA IS DELIBERATE !                     
C   JAN89: USE YMF FACTOR WITH SELF-SPUTTERED NEUTRALS IF FLAG = 1              
C-----------------------------------------------------------------------        
C                                                                               
  780   CONTINUE                                                                
c slmod begin
c
c A little bit of statistics:
c
        IF (optdp.EQ.1) THEN
          IF (TLOSS(IMP).NE.0.0.AND.IFATE.GE.2.AND.IFATE.LE.5) THEN
            ALOSS = ALOSS + TLOSS(IMP)
            LLOSS = LLOSS + 1.0
          ENDIF
        ENDIF
c slmod end
        RMACH = ABS (CVABS/QTIM)                                                
c
c       Impact energy option
c
        if (impact_energy_opt.eq.0) then 
           ENERGY = 3.0 * REAL(CIZ) * CTBIQX +                                     
     >       5.22E-9 * CRMI * RMACH * RMACH + 2.0 * CTEMI                          
        elseif (impact_energy_opt.eq.1) then 
           ENERGY = 3.0 * REAL(CIZ) * CTBIQX + 2.0 * CTEMI                          
        endif
c
        NEROYS(IOY,1) = NEROYS(IOY,1) + SPUTY                                   
        NERODS(IOD,1) = NERODS(IOD,1) + SPUTY                                   
        NERODS3(IOD,IP,1) = NERODS3(IOD,IP,1) + SPUTY                                   
c        
c      jdemod
c
c        write(6,'(a,2i8,12(1x,g12.5))') 'DEP:',iod,ip,sputy,
c     >       nerods3(iod,ip,1),alpha,y,p

C                                                                               
        KK = KK + 1                                                             
c
c        write(0,'(a,2i8,10(1x,g18.8))') 'Dep:',imp,ciab,cx,y,p
c
        IF (CIAB.EQ.-1.OR.CIAB.EQ.2) THEN                                       
          DEPS(IX,CIZ,1) = DEPS(IX,CIZ,1) + SPUTY                               
          NEROXS(IX,1,1) = NEROXS(IX,1,1) + SPUTY                               
          RDEP(1)   = RDEP(1) + SPUTY                                     
          IF (CNEUTD.EQ.8)
     >       ENERGY = REAL(CIZ) * CVS(IQX,1)
     >              + 5.22E-9 * CRMI * RMACH * RMACH
     >              + 2.0 * CTEMI
          RYIELD    = YIELD (6, MATLIM, ENERGY,
     >                       ctembs(ix,iy),ctembsi(ix,iy))
     >                      *QMULTS*CYMFSS(IQX,1)           
          RESPUT    = RES (6,MATLIM,RYIELD,.FALSE.,CNEUTD,RANV(KK),
     >                     QMULTS)
          SPUNEW    = SPUTY * RYIELD                                            
C
C         LIMIT THE MAXIMUM VALUE OF A SPUTTERED FRAGMENT TO THE 
C         INPUT VALUE CSPUMAX. THIS CAN BE USED TO CONTROL 
C         SITUATIONS THAT MIGHT DEVELOP INTO A LOCAL RUNAWAY.
C         THE EFFECT CAN BE NEGATED BY MAKING THE VALUE OF 
C         CSPUMAX APPROPRIATELY LARGE. IN MOST CASES A VALUE OF 
C         2.0 WOULD BE ADEQUATE, HOWEVER, A DEFAULT VALUE OF 100.0
C         IS USUALLY ASSIGNED.
C
          SPUNEW = MIN(SPUNEW,CSPUMAX)
C
          YLDTOT(1) = YLDTOT(1) + SPUNEW                                        
          YLDMAX(1) = MAX (YLDMAX(1), SPUNEW)  
        ELSE                                                                    
          DEPS(IX,CIZ,2) = DEPS(IX,CIZ,2) + SPUTY                               
          NEROXS(IX,1,2) = NEROXS(IX,1,2) + SPUTY                               
          RDEP(2)   = RDEP(2) + SPUTY                                           
          IF (CNEUTD.EQ.8)
     >       ENERGY = REAL(CIZ) * CVS(IQX,2)
     >              + 5.22E-9 * CRMI * RMACH * RMACH
     >              + 2.0 * CTEMI
          RYIELD    = YIELD (6, MATLIM, ENERGY,
     >                       ctembs(ix,iy),ctembsi(ix,iy))
     >                      *QMULTS*CYMFSS(IQX,2)           
          RESPUT    = RES (6,MATLIM,RYIELD,.FALSE.,CNEUTD,RANV(KK),
     >                     QMULTS)             
          SPUNEW    = SPUTY * RYIELD                                            
C
C         CSPUMAX - AS ABOVE 
C
          SPUNEW = MIN(SPUNEW,CSPUMAX)
C
          YLDTOT(2) = YLDTOT(2) + SPUNEW                                        
          YLDMAX(2) = MAX (YLDMAX(2), SPUNEW)                                   
        ENDIF                                                                   
        RANDEP = RANDEP + RANV(KK) * SPUTY                                      
C                                                                               
        IF (SPUNEW.GT.CTRESH) THEN                                              
          IF (NPROD.GE.MAXIMP) THEN                                             
            WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,CX,ALPHA,Y,P,
     >        SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,'SELFSPUT FAILED'                  
          ELSE                                                                  
            IF (CIAB.EQ.-1.OR.CIAB.EQ.2) THEN                                   
              YTHTOT(1) = YTHTOT(1) + SPUNEW                                    
            ELSE                                                                
              YTHTOT(2) = YTHTOT(2) + SPUNEW                                    
            ENDIF                                                               
            NPROD = NPROD + 1                                                   
            SNEWS(NPROD) = SPUNEW                                               
            IF     (RESPUT) THEN                                                
              RMAXS(NPROD) =-1.0                                                
            ELSEIF (CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5.OR.               
     >              CNEUTD.EQ.4) THEN                                           
              EMAX = CEMAXF * ENERGY                                            
              RMAXS(NPROD) = 1.0 / ((1.0+CEBD/EMAX) * (1.0+CEBD/EMAX))          
            ELSE                                                                
              RMAXS(NPROD) = 1.0                                                
            ENDIF                                                               
            XPRODS(NPROD) = XM                                                  
            YPRODS(NPROD) = YM                                                  
            PPRODS(NPROD) = P                                                   
          ENDIF                                                                 
        ENDIF                                                                   
C      (GOTO 790)  NOT ACTUALLY NEEDED!                                         
C                                                                               
C-----------------------------------------------------------------------        
C       CURRENT ION OR SUB-ION FINISHED WITH.  UPDATE NUMBER OF IONS            
C       REACHING STATE, ETC.  CARE NEEDED WITH ROULETTED SPUTY VALUES...        
C       SEE IF ANY FURTHER SPLIT IONS EXIST - IF SO RETRIEVE THE DETAILS        
C       OF THE SPLITTING POSITION AND LAUNCH ANOTHER SUB-ION.                   
C       THE RESETTING OF ION DETAILS AT THE POINT OF SPLITTING IS               
C       SOMEWHAT OVER-ENGINEERED BELOW (EG ABSY DOESN'T REALLY NEED             
C       RESETTING), BUT THERE IS NO HARM IN THIS.                               
C-----------------------------------------------------------------------        
C                                                                               
  790   CONTINUE                                                                
        IF (IFATE.NE.7 .AND. IFATE.NE.8) THEN                                   
          CISTOT = CISTOT + CIST * SPUTY                                        
          CISMAX = MAX (CISMAX, CIST)                                           
          DO 792 IZ = CIZSC, MAXCIZ                                             
            RIONS(IZ) = RIONS(IZ) + SPUTY                                       
  792     CONTINUE                                                              
        ENDIF                                                                   


        !
        ! Update particle reflection statistics when this ion is finished. 
        !

        call update_part_refl_stats(sputy)


C                                                                               
        IF (DEBUGL) WRITE (6,9003) IMP,CIST+QFACT,IQX,IQY,IX,IY,                
     >    CX,ALPHA,Y,P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,FATE(IFATE)            
c
c       Record particle track to permanent array.)  
c
        if (debugt) then 
c
c          Deal with the last position of the particle
c
           tptrac(traclen,1) = cx
           tptrac(traclen,2) = y 
           traclen = traclen + 1
           if (traclen.gt.maxlen) then
              bigtrac = .true.
              traclen = 1
           endif   
c 
c          Move the data to the permanent array
c 
c          Map the positions to the middle limiter - if necessary
c
c     jdemod - this comment doesn't make sense since it shifts
c     the particle track based on the last position of the
c     particle ... potentially moving all the rest of the particle
c     track away from the middle limiter          
c
c     It might be reasonable if mapping a particle track that is near
c     one end or the other but doesn't work for ones that cross the
c     CL boundaries
c     
           if (bigtrac) then 
c
              if (traclen.eq.1) then
                 if (tptrac(maxlen,2).gt.cl) then
                    do 2040 in = 1,maxlen
                       tptrac(in,2) = tptrac(in,2) - ctwol     
 2040               continue
                 elseif (tptrac(maxlen,2).lt.-cl) then
                    do 2050 in = 1,maxlen
                       tptrac(in,2) = tptrac(in,2) + ctwol     
 2050               continue
                 endif
              elseif (tptrac(traclen-1,2).gt.cl) then
                 do 2060 in = 1,maxlen
                    tptrac(in,2) = tptrac(in,2) - ctwol     
 2060            continue
              elseif (tptrac(traclen-1,2).lt.-cl) then
                 do 2070 in = 1,maxlen
                    tptrac(in,2) = tptrac(in,2) + ctwol     
 2070            continue
              endif
c
              ptracl(imp) = maxlen
c
              do 2000 in = traclen, maxlen
                 ptracs(in-traclen+1,imp,1) = tptrac(in,1)
                 ptracs(in-traclen+1,imp,2) = tptrac(in,2)
 2000         continue
              do 2010 in = 1,traclen-1          
                 ptracs(maxlen-traclen+in+1,imp,1) = tptrac(in,1)
                 ptracs(maxlen-traclen+in+1,imp,2) = tptrac(in,2)
 2010         continue
c
           else
c
              ptracl(imp) = traclen-1
c
              if (tptrac(traclen-1,2).gt.cl) then
                 do 2080 in = 1,traclen -1
                    tptrac(in,2) = tptrac(in,2) - ctwol     
 2080            continue
              elseif (tptrac(traclen-1,2).lt.-cl) then
                 do 2090 in = 1,traclen -1
                    tptrac(in,2) = tptrac(in,2) + ctwol     
 2090            continue
              endif
c
              do 2020 in = 1,traclen-1          
                 ptracs(in,imp,1) = tptrac(in,1)
                 ptracs(in,imp,2) = tptrac(in,2)
 2020         continue
           endif
           write(6,*) 'imp2:',imp,ptracl(imp),bigtrac
     >                   ,traclen
        endif   
C                                                                               
  794   CONTINUE                                                                
        IF (IGET(IPUT).GT.0 .AND. IGET(IPUT).LE.CNSPL) THEN                     
          IF (IGET(IPUT).GT.1) THEN                                             
            CX     = R(1,IPUT)                                                  
            Y      = R(2,IPUT)                                                  
            P      = R(3,IPUT)                                                  
            QUANT  = R(4,IPUT)                                                  
            SVY    = R(5,IPUT)                                                  
            CTEMI  = R(6,IPUT)                                                  
            CIST   = R(7,IPUT)                                                  
            SPUTY  = R(8,IPUT)                                                  
            TSTEPL = R(9,IPUT)                                                  
c
c           jdemod - add time since t=0 to split/roulette data
c
            RTIME  = R(10,IPUT)
            IQX    = I(1,IPUT)                                                  
            IQY    = I(2,IPUT)                                                  
            IX     = I(3,IPUT)                                                  
            IY     = I(4,IPUT)                                                  
            IP     = I(5,IPUT)                                                  
            IT     = I(6,IPUT)                                                  
            CIZ    = I(7,IPUT)                                                  
            MAXCIZ = I(8,IPUT)                                                  
            IS     = I(9,IPUT)                                                  
            DIFFUS = L(1,IPUT)                                                  
            CFLRXA = L(2,IPUT)                                                  
            CFLY2L = L(3,IPUT)                                                  
            Delta_Y1    = 0.0D0                                                      
            Delta_Y2    = 0.0d0                                                   
            Y_position = dble(y)

            ABSY   = ABS (Y)                                                    
            YY     = MOD (ABSY, CL)                                             
            IF (YY.GT.CHALFL) YY = CL - YY                                      
            IF (CX.GT.0.0) THEN                                                 
              ALPHA  = CX / (YY*CONI + 1.0)                                     
            ELSE                                                                
              ALPHA  = CX / (YY*CONO + 1.0)                                     
            ENDIF                                                               
            DTEMI  = DBLE (CTEMI)                                               
            DIST   = DBLE (CIST)                                                
            ! jdemod - update DTIME for split/roulette
            DTIME  = DBLE (RTIME)
            JY     = IABS (IY)                                                  
          ENDIF                                                                 
          SPUTY  = SPUTY / REAL(CNSPL)                                          
          DSPUTY = DBLE (SPUTY)                                                 
          IF (DEBUGL) WRITE (6,9003) IMP,CIST,IQX,IQY,IX,IY,CX,ALPHA,Y,         
     >      P,SVY,CTEMI,SPARA,SPUTY,IP,IT,IS,'SUB-ION:',IGET(IPUT),IPUT       
          IGET(IPUT) = IGET(IPUT) + 1                                           
          GOTO 500                                                              
        ENDIF                                                                   
C                                                                               
        IF (IPUT.GT.0) THEN                                                     
          IPUT = IPUT - 1                                                       
          GOTO 794                                                              
        ENDIF                                                                   
C                                                                               
C-----------------------------------------------------------------------        
C       SEE IF TEST OF CPU TIME USED IS DUE                                     
C       TRAP CASE WHERE TIMUSD=0 OCCURS.                                        
C-----------------------------------------------------------------------        
C                                                                               
        IF (IMP.GE.IMPLIM) THEN                                                 
         TIMUSD = ZA02AS(1) - STATIM                                            
         IF (TIMUSD.GT.0.0) THEN                                                
           PARTIM = (CPULIM-NEUTIM-IONTIM) / TIMUSD                             
         ELSE                                                                   
           PARTIM = 10.0                                                        
         ENDIF                                                                  
         IF (PARTIM.GE.1.05) THEN                                               
           IMPLIM = INT (REAL(IMP) * MIN (4.0,0.25+0.75*PARTIM))                
         ELSE                                                                   
C                                                                               
C--------- HAVE RUN OUT OF CPU TIME, STOP ITERATION                             
C--------- THERE ARE SO MANY COUNTERS IT IS VIRTUALLY IMPOSSIBLE TO             
C--------- WIND UP THE ROUTINE CLEANLY.  JUST WORK OUT HOW MANY IONS            
C--------- HAVE BEEN FOLLOWED IN THE TIME ALLOTTED AND CALL IT A DAY.           
C                                                                               
           IMPLIM = IMP                                                         
           DO 791 J = IMP+1, NATIZ                                              
             RATIZ = RATIZ - SPUTYS(J)                                          
  791      CONTINUE                                                             
           NATIZ = IMP                                                          
           WRITE (6,'('' ERROR:  CPU TIME LIMIT REACHED'')')                    
           WRITE (6,'('' NUMBER OF IONS REDUCED TO'',I15)') NINT(RATIZ)          
           WRITE (0,'('' ERROR:  CPU TIME LIMIT REACHED'')')                    
           WRITE (0,'('' NUMBER OF IONS REDUCED TO'',I15)') NINT(RATIZ)          
           CALL PRB                                                             
           CALL PRC ('ERROR:  CPU TIME LIMIT REACHED')                          
           CALL PRI ('NUMBER OF IMPURITY IONS REDUCED TO ',NINT(RATIZ))         
           CALL PRC ('INCREASE CPU TIME OR INCREASE QUANTUM TIMESTEP')          
           CALL PRC ('RESULTS WILL BE PRINTED BUT THEY SHOULD BE TREATED        
     > WITH CAUTION ...')                                                       
           GOTO 805                                                             
         ENDIF                                                                  
        ENDIF                                                                   

        
C                                                                               
c     800 is the end of the main particle following loop
c
  800 CONTINUE                                                                  
C                                                                               
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
C                                                                               
  805 CONTINUE                                                                  
C                                                                               
C---- ITERATION OF ALL IONS COMPLETE.  LAUNCH SECONDARY, TERTIARY, ETC          
C---- NEUTRALS AND LEAP BACK TO FOLLOW NEXT LOT OF IONS CREATED.                
C---- YTHTOT AND YLDMAX ARE PRETTY MEANINGLESS IF SPLITTING/ROULETTING          
C---- ARE IN USE SO THERE IS NO POINT IN PRINTING THEM.                         
C---- DON'T BOTHER TO RELAUNCH IF "TOTAL YIELD" IS LESS THAN (SAY) 10.0.        
C                                                                               
      CALL PRB                                                                  
      CALL PRR('AVERAGE X POSITION OF ION INJECTION      ',AVXPOS/RATIZ)        
      CALL PRR('AVERAGE ALPHA POSITION OF ION INJECTION  ',AVAPOS/RATIZ)        
      CALL PRR('AVERAGE ABS(Y) POS OF ION INJECTION      ',AVYPOS/RATIZ)        
      CALL PRR('AVERAGE ABS(P) POS OF ION INJECTION      ',AVPPOS/RATIZ)        
      CALL PRR('AVERAGE TIME TO FIRST DIFFUSION          ',RDIFFT/RATIZ)        
      CALL PRI2('NUMBER OF IONS PLATING OUT ON WALLS',                          
     >                                 NINT(RWALL (1)), NINT(RWALL (2)))        
      CALL PRI2('NUMBER OF IONS PLATING ON LIMITERS ',                          
     >                                 NINT(RDEP  (1)), NINT(RDEP  (2)))        
      IF     (CYMFLG.EQ.-1) THEN                                                
        CALL PRI2('TOTAL SELF-SPUTTERING YIELD        ',                        
     >                                 NINT(YLDTOT(1)), NINT(YLDTOT(2)))        
      ELSE                                                                      
        CALL PRI2('TOTAL SELF-SPUT YIELD  (WITH YMF)  ',                        
     >                                 NINT(YLDTOT(1)), NINT(YLDTOT(2)))        
      ENDIF                                                                     
      CALL PRI2('TOTAL YIELDS ABOVE THRESHOLD YIELD ',                          
     >                                 NINT(YTHTOT(1)), NINT(YTHTOT(2)))        
      CALL PRR2('MAXIMUM YIELD ANY NEUTRAL FRAGMENT ',                          
     >                                              YLDMAX(1),YLDMAX(2))        
      WRITE (6,'(1X,A,I9 )') 'AVERAGE NUMBER OF ITERATIONS PER ION   ',         
     >                                               NINT(CISTOT/RATIZ)         
      WRITE (6,'(1X,A,I13)') 'MAXIMUM NUMBER OF ITERATIONS       ',             
     >                                                     NINT(CISMAX)
      IF ((RDEP(1).NE.0).OR.(RDEP(2).NE.0)) THEN   
        WRITE (6,*) ' AV. RANDOM USED FOR RES ',RANDEP/(RDEP(1)+RDEP(2))
      ENDIF
      
      TATIZ  = TATIZ  + RATIZ                                                   
      TAVXPOS = TAVXPOS + AVXPOS 
      TDEP   = TDEP   + RDEP(1)  + RDEP(2)                                      
      TWALL  = TWALL  + RWALL(1) + RWALL(2)                                     
      TNEUT  = TNEUT  + RNEUT                                                   
      TRES   = TRES   + RRES                                                    
      TWALLN = TWALLN + RWALLN                                                  
      TCENT  = TCENT  + RCENT                                                   
      TTMAX  = TTMAX  + RTMAX                                                   
      TSTRUK = TSTRUK + RSTRUK                                                  
      TFAIL  = TFAIL  + RFAIL                                                   
      IONTIM = IONTIM + ZA02AS (1) - STATIM                                     
C                                                                               
      IF ((CNEUTD.EQ.3.OR.CNEUTD.EQ.4.OR.CNEUTD.GE.5).AND.(NPROD.GT.0)         
     >    .AND.((YTHTOT(1)+YTHTOT(2)).GT.2.0).AND.(STATUS.LT.CMAXGENS)) 
     >  THEN          
        STATUS = STATUS + 1                                                     
        CNEUTA = 0                                                              
        CNEUTB = 2                                                              
        IF (CNEUTD.EQ.4) CNEUTC = 1                                             
        IF (CNEUTE.EQ.1) CNEUTE = 0                                             
        IF (STATUS.LE.50) THEN 
           WRITE(6,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'      
           WRITE(7,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'      
        ELSE
           WRITE (6,9013) '***  LAUNCHING ',STATUS,
     >                    ' GENERATION NEUTRALS  ***'      
           WRITE (7,9013) '***  LAUNCHING ',STATUS,
     >                    ' GENERATION NEUTRALS  ***'      
        ENDIF
        DO 8700 IMP = 1, NPROD                                                  
          SPUTYS(IMP) = SNEWS(IMP)                                              
 8700   CONTINUE                                                                
        CALL LAUNCH (FSRATE,1,NPROD,1,NATIZ,RSTRUK,RRES,                        
     >               RATIZ,RNEUT,RWALLN,RCENT,RTMAX,SEED,NRAND,                 
     >               NEUTIM,RFAIL,STATUS,6,MATLIM,QMULTS)                    
        IF (NATIZ.GT.0) GOTO 200                                                
      ENDIF                                                                     
C                                                                               
      SSEF = 0.0                                                                
      YEFF = 0.0                                                                
      IF (RNEUT1.GT.0.0) SSEF = (TNEUT-RNEUT1) / RNEUT1                         
      IF (GTOT1.GT.0.0)  YEFF = (1.0 + SSEF) * GYTOT1 / GTOT1                   
C                                                                               
  806 CONTINUE                                                                  
      CALL PRB                                                                  
      CALL PRC ('***  S U M M A R Y   D E T A I L S  ***')                      
      CALL PRB                                                                  
      IF (CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7) THEN            
      CALL PRI('NO OF (HIGH) PRIMARY NEUTRALS LAUNCHED   ',                     
     >                                  NINT(RNEUT1-RRES1))                     
      CALL PRI('NO OF (LOW)  PRIMARY NEUTRALS LAUNCHED   ',NINT(RRES1))         
      CALL PRI('NO OF (HIGH) SELF-SPUTTERED NEUTRALS     ',                     
     >                                  NINT(TNEUT-RNEUT1-(TRES-RRES1)))        
      CALL PRI('NO OF (LOW)  SELF-SPUTTERED NEUTRALS     ',                     
     >                                  NINT(TRES-RRES1))                       
      ELSE                                                                      
      CALL PRI('NUMBER OF PRIMARY NEUTRALS LAUNCHED      ',NINT(RNEUT1))        
      CALL PRI('NUMBER OF SELF-SPUTTERED NEUTRALS        ',                     
     >                                  NINT(TNEUT-RNEUT1))                     
      ENDIF                                                                     
      CALL PRI('TOTAL NUMBER OF NEUTRALS LAUNCHED        ',NINT(TNEUT))         
      CALL PRI('TOTAL NO OF NEUTRALS PLATING ON WALLS    ',NINT(TWALLN))        
      CALL PRI('TOTAL NO OF NEUTRALS REACHING CENTRE     ',NINT(TCENT))         
      CALL PRI('TOTAL NO OF NEUTRALS EXISTING AT TMAX    ',NINT(TTMAX))         
      CALL PRI('TOTAL NO OF NEUTRALS STRIKING LIMITER    ',NINT(TSTRUK))        
      CALL PRI('TOTAL NO OF FAILED NEUTRAL LAUNCHES      ',NINT(TFAIL))         
      CALL PRB                                                                  
      CALL PRI('TOTAL NUMBER OF IONS CREATED             ',NINT(TATIZ))         
      CALL PRI('TOTAL NO OF IONS PLATING ON WALLS        ',NINT(TWALL))         
      CALL PRI('TOTAL NO OF IONS PLATING ON LIMITERS     ',NINT(TDEP))          
      CALL PRI('TOTAL NO OF IONS IONISED BEYOND LIMIT    ',NINT(TBYOND))        
      CALL PRI('TOTAL IONS RECOMBINED TO FORM NEUTRALS   ',NINT(TBELOW))        
      CALL PRI('TOTAL NO OF IONS EXISTING AT TMAX        ',NINT(TCUT))          
      CALL PRR('TOTAL AVERAGE X ION INJECTION POSITION   ',
     >          TAVXPOS/TATIZ)
      CALL PRR('SELF-SPUTTERING ENHANCEMENT FACTOR       ',SSEF)                
      CALL PRR('EFFECTIVE YIELD                          ',YEFF)                
C                                                                               
      IF (NSPLIT(1).GT.0) THEN                                                  
        CALL PRB                                                                
        CALL PRB                                                                
        CALL PRC ('*** SPLITTING & ROULETTING ACTIVITY ***')                    
        CALL PRB                                                                
        CALL PRI ('NO OF SUB-IONS TO BE CREATED BY SPLITS   ', CNSPL)           
        CALL PRR ('PROBABILITY OF RETENTION ON ROULETTING   ', CPRUL)           
        CALL PRI ('GREATEST SIZE OF SPLITTING BANK          ', MPUT)            
        WRITE (7,9010)                                                          
        DO 8800 IS = 1, MAXINS                                                  
          IF (NSPLIT(IS).GT.0)                                                  
     >      WRITE (7,9011) IS,CXSPLS(IS),                                       
     >        NINT(TSPLIT(IS)),NSPLIT(IS),NINT(TRULET(IS)),NRULET(IS)           
 8800   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C     SCALE ARRAYS DDLIMS, TIZS, ETC   ...                                      
C-----------------------------------------------------------------------        
C                                                                               
      do 3990 iqx = 1-nqxso,nqxsi
        if (svyacc(iqx).gt.0.0) svybar(iqx) = svybar(iqx)/svyacc(iqx)
3990  continue 
C
C---- SET FACTORS ARRAY FOR NORMALISING GRAPHS TO 1 PARTICLE                    
C                                                                               
C
C     SET THE DEFAULT VALUE FOR THESE SCALING CONSTANTS TO 0.0
C     THEN IF THE NUMBER OF NEUTRALS LAUNCHED IS ZERO - AS OCCURS
C     FOR SOME ION INJECTION CASES - ALLOW THE SCALING FACTOR TO 
C     BE SET TO THE INVERSE OF THE NUMBER OF IONS INJECTED.
C
      FACTA(-1) = 0.0                                                           
      IF (RNEUT1.GT.0.0) THEN 
        FACTA(-1) = 1.0 / RNEUT1                               
      ELSEIF (TATIZ.GT.0.0) THEN
        FACTA(-1) = 1.0 / TATIZ
      ENDIF
      FACTB(-1) = FACTA(-1) * FSRATE                                            
      FACTA(0) = 0.0                                                            
      IF (TNEUT.GT.0.0) THEN 
         FACTA(0) = 1.0 / TNEUT                                  
      ELSEIF (TATIZ.GT.0.0) THEN 
         FACTA(0) = 1.0 / TATIZ
      ENDIF  
      FACTB(0) = FACTA(0) * FSRATE                                              
C                                                                               
      IF (NIZS.GT.0) THEN                                                       
        DO 4000 IZ = 1, NIZS                                                    
          IF (TATIZ.GT.0.0) THEN                                                
            FACTA(IZ) = 1.0 / TATIZ                                             
          ELSE                                                                  
            FACTA(IZ) = 0.0                                                     
          ENDIF                                                                 
          FACTB(IZ) = FACTA(IZ) * QTIM                                          
 4000   CONTINUE                                                                
      ENDIF                                                                     
      WRITE (6,*) '1:'
C                                                                               
C-----------------------------------------------------------------------        
C     SECTION FOR SYMMETRIC CONTRIBUTIONS.  NOTES 75,219                        
C-----------------------------------------------------------------------        
C                                                                               
      ! jdemod - symmetric contributions for ddvs and ddvs2 are
      ! done in the routine finish_diagvel below. 
      DO 4130 IZ = -1, NIZS                                                     
       DO 4120 IX = 1, NXS                                                      
        DO 4100 IY = 1, NYS                                                     
          IF (IZ.GT.0) THEN                                                     
            DEMP( IY,1) = DDTS(IX, IY,IZ) + DDTS(IX,-NYS-1+IY,IZ)               
            DEMP(-IY,1) = DDTS(IX,-IY,IZ) + DDTS(IX, NYS+1-IY,IZ)               
            DEMP( IY,2) = DDYS(IX, IY,IZ) + DDYS(IX,-NYS-1+IY,IZ)               
            DEMP(-IY,2) = DDYS(IX,-IY,IZ) + DDYS(IX, NYS+1-IY,IZ)               
            ! jdemod - symmetric contributions for ddvs and ddvs2 are
            ! done in the routine finish_diagvel below. 
            !if (debugv) then
            !  DEMP( IY,5) = DDVS(IX, IY,IZ) + DDVS(IX,-NYS-1+IY,IZ)               
            !  DEMP(-IY,5) = DDVS(IX,-IY,IZ) + DDVS(IX, NYS+1-IY,IZ)                           
            !endif
          ENDIF                                                                 
          DEMP( IY,3) = DBLE (TIZS(IX, IY,IZ) + TIZS(IX,-NYS-1+IY,IZ))          
          DEMP(-IY,3) = DBLE (TIZS(IX,-IY,IZ) + TIZS(IX, NYS+1-IY,IZ))          
          DEMP( IY,4) = DDLIMS(IX, IY,IZ) + DDLIMS(IX,-NYS-1+IY,IZ)             
          DEMP(-IY,4) = DDLIMS(IX,-IY,IZ) + DDLIMS(IX, NYS+1-IY,IZ)             
 4100   CONTINUE                                                                
        DO 4110 IY = -NYS, NYS                                                  
          IF (IZ.GT.0) THEN                                                     
            DDTS(IX,IY,IZ) = DEMP(IY,1)                                         
            DDYS(IX,IY,IZ) = DEMP(IY,2)                                         
            !if (debugv) then 
            !   DDVS(IX,IY,IZ) = DEMP(IY,5)                                         
            !endif
          ENDIF                                                                 
          TIZS(IX,IY,IZ) = SNGL(DEMP(IY,3))                                     
          DDLIMS(IX,IY,IZ) = DEMP(IY,4)                                         
 4110   CONTINUE                                                                
 4120  CONTINUE                                                                 
 4130 CONTINUE                                                                  

c
c     jdemod - in original LIM the 3D arrays could be smaller than the 2D arrays so they only recorded
c     part of the 3D space. This means that particles outside the region were ignored and the
c     symmetric contributions were excluded. The entire symmetric aspect will be reomoved in a code
c     rewrite to update the 3D features but until then - if ny3d = nys then the symmetric contributions in the 3D      
c     arrays ddlim3 and lim5 will be processed - otherwise half the statistics are being left out.
c
      if (ny3d.eq.nys) then 
         do iz = -1,nizs
            do ix = 1,nxs
               do iy = 1,ny3d
                 do ip = -maxnps,maxnps
                    ddlim3(ix,iy,iz,ip) = ddlim3(ix,iy,iz,ip) + 
     >                              ddlim3(ix,iy-ny3d-1,iz,ip)  
                    ddlim3(ix,iy-ny3d-1,iz,ip) = ddlim3(ix,iy,iz,ip)

                    do it = 1,nts
                       lim5(ix,iy,iz,ip,it) = lim5(ix,iy,iz,ip,it)+
     >                                   lim5(ix,iy-ny3d-1,iz,ip,it)
                       lim5(ix,iy-ny3d-1,iz,ip,it)=lim5(ix,iy,iz,ip,it)
                    end do

                end do
              end do
            end do
          end do
      endif
      
c
c     jdemod - symmetric contributions in the walls array
c            - nys goes to -2L and nys to +2L 
c            - the ranges -2L,0 maps directly onto 0,2L - they are the same
c            - so data is combined and overlaid
c            - I have no idea why this design decision was originally made
c
      do iz = -1,nizs
         do iy = 1,nys
            walls(iy,iz) = walls(iy,iz) + walls(-nys-1+iy,iz)
            walls(-nys-1+iy,iz) = walls(iy,iz)
         end do 
      end do

C                                                                               
C====================== DDTS AND DDYS ARRAYS ===========================        
C                                                                               
C---- CALCULATE STEADY STATE CLOUD TEMPERATURES AND AVERAGE DELTAY STEPS        
C---- THE VALUES MUST BE STORED IN D.P. TO GET ENOUGH ACCURACY,                 
C---- AND THE SUMMATION IS ALSO DONE IN D.P..                                   
C                                                                               
      DO 4290 IZ = 1, NIZS                                                      
C                                                                               
        DO 4210 IX = 1, NXS                                                     
          DSUM1 = 0.0D0                                                         
          DSUM2 = 0.0D0                                                         
          DSUM3 = 0.0D0                                                         
          DO 4200 IY = -NYS, NYS                                                
            DSUM1 = DSUM1 + DDTS  (IX,IY,IZ)                                    
            DSUM2 = DSUM2 + DDYS  (IX,IY,IZ)                                    
            DSUM3 = DSUM3 + DDLIMS(IX,IY,IZ)                                    
 4200     CONTINUE                                                              
          IF (DSUM3.GT.0.0D0) THEN                                              
            SDTXS(IX,IZ) = DSUM1 / DSUM3                                        
            SDYXS(IX,IZ) = DSUM2 / DSUM3                                        
          ENDIF                                                                 
 4210   CONTINUE                                                                
C                                                                               
        DO 4230 IY = -NYS, NYS                                                  
          DSUM1 = 0.0D0                                                         
          DSUM2 = 0.0D0                                                         
          DSUM3 = 0.0D0                                                         
          DO 4220 IX = 1, NXS                                                   
            DSUM1 = DSUM1 + DDTS  (IX,IY,IZ)                                    
            DSUM2 = DSUM2 + DDYS  (IX,IY,IZ)                                    
            DSUM3 = DSUM3 + DDLIMS(IX,IY,IZ)                                    
 4220     CONTINUE                                                              
          IF (DSUM3.GT.0.0D0) THEN                                              
            SDTYS(IY,IZ) = DSUM1 / DSUM3                                        
            SDYYS(IY,IZ) = DSUM2 / DSUM3                                        
          ENDIF                                                                 
 4230   CONTINUE                                                                
C
C      THE CENTRAL VALUE OF THE AVERAGE TEMPERATURE Y ARRAY IS 
C      NOT SIGNIFICANT BECAUSE THE CENTRAL BIN IS NOT USED. BUT 
C      IT DOES SHOW UP IN PLOTTING. THEREFORE, SET IT TO THE AVERAGE
C      OF THE +1 AND -1  VALUES.
C
C      D. ELDER APRIL 6, 1990
C
        SDTYS(0,IZ) = (SDTYS(1,IZ)+SDTYS(-1,IZ))/2.0
C                                                                               
        DSUM1 = 0.0D0                                                           
        DSUM2 = 0.0D0                                                           
        DSUM3 = 0.0D0                                                           
        DSUM4 = 0.0D0                                                           
        DO 4250 IY = -NYS, NYS                                                  
          DO 4240 IX = 1, NXS                                                   
            DSUM1 = DSUM1 + DDTS  (IX,IY,IZ)                                    
            DSUM3 = DSUM3 + DDLIMS(IX,IY,IZ)                                    
            IF (XS(IX).LE.0.0) THEN                                             
              DSUM2 = DSUM2 + DDYS  (IX,IY,IZ)                                  
              DSUM4 = DSUM4 + DDLIMS(IX,IY,IZ)                                  
            ENDIF                                                               
 4240     CONTINUE                                                              
 4250   CONTINUE                                                                
        IF (DSUM3.GT.0.0D0) SDTZS(IZ) = DSUM1 / DSUM3                           
        IF (DSUM4.GT.0.0D0) SDYZS(IZ) = DSUM2 / DSUM4                           
C                                                                               
        DO 4270 IY = -NYS, NYS                                                  
          DO 4260 IX = 1, NXS                                                   
            IF (DDLIMS(IX,IY,IZ).GT.0.0D0) THEN                                 
              DDTS(IX,IY,IZ) = DDTS(IX,IY,IZ) / DDLIMS(IX,IY,IZ)                
              DDYS(IX,IY,IZ) = DDYS(IX,IY,IZ) / DDLIMS(IX,IY,IZ)                
            ENDIF                                                               
 4260     CONTINUE                                                              
 4270   CONTINUE                                                                
 4290 CONTINUE                                                                  
c      WRITE(6,*) '3:'
C                                                                               
C
      ! jdemod - normalize the velocity diagnostic data before ddlims is converted to density
      ! jdemod - symmetric contributions for ddvs and ddvs2 are
      !          done in the routine finish_diagvel  
      call finish_diagvel(qtim,nizs)
C
C
C                                                                               
C======================== TIZS AND TIZ3 ARRAYS =========================        
C                                                                               
      DO 4430 IZ = -1, NIZS                                                     
       DO 4420 IX = 1, NXS                                                      
        DO 4410 IY = 1, NYS                                                     
          FACT = FACTA(IZ) /                                                    
     >           (XWIDS(IX) * XCYLS(IX) * YWIDS(IY) * DELPS(IX,IY))             
          TIZS(IX, IY,IZ) = FACT * TIZS(IX, IY,IZ)                              
          TIZS(IX,-IY,IZ) = FACT * TIZS(IX,-IY,IZ)                              
          IF (IY.LE.NY3D) THEN                                                  
            DO 4400 IP = -MAXNPS, MAXNPS                                        
              TIZ3(IX, IY,IZ,IP) = FACT/pwids(ip) * TIZ3(IX, IY,IZ,IP)                    
              TIZ3(IX,-IY,IZ,IP) = FACT/pwids(ip) * TIZ3(IX,-IY,IZ,IP)                    
 4400       CONTINUE                                                            
          ENDIF                                                                 
 4410   CONTINUE                                                                
 4420  CONTINUE                                                                 
 4430 CONTINUE                                                                  
      WRITE(6,*) '4:'
C                                                                               
C======================== LIM5 ARRAY ===================================        
C                                                                               
      IF (IMODE.NE.2) THEN                                                      

         !if (ANY(lim5.ne.0.0)) then
         !   write(0,*) 'LIM5A - non-zero elements found'
         !endif

       DO 4540 IZ = -1, NIZS                                                    
        DO 4530 IX = 1, NXS                                                     
         DO 4520 IY = 1, NY3D                                                   
!     jdemod - remove cdwelt_sum option functionality because it isn't
!              physically meaningful.                  
!           if (cdwelt_sum.eq.0) then 
              FACT = FACTA(IZ) /                                                    
     >           (XWIDS(IX) * XCYLS(IX) * YWIDS(IY) * DELPS(IX,IY))             
!           elseif (cdwelt_sum.eq.1) then 
!               FACT = FACTB(IZ) /                                                    
!     >           (XWIDS(IX) * XCYLS(IX) * YWIDS(IY) * DELPS(IX,IY))             
!           endif
           
           
           
           DO 4510 IT = 1, NTS                                                   
             DO 4500 IP = -MAXNPS, MAXNPS                                         

c               if (lim5(ix,iy,iz,ip,it).ne.0.0.or.
c     >              lim5(ix,-iy,iz,ip,it).ne.0.0) then
c
c                   write(6,*) 'LIM5:'
c
c                
c                   write(6,'(a,6i8,3(1x,g12.5),1x,2l5)')
c     >                 'LIM5:',ix,iy,iy-ny3d-1,iz,ip,it,
c     >                   fact,lim5(ix,iy,iz,ip,it),
c     >                  lim5(ix,iy-ny3d-1,iz,ip,it),
c     >               lim5(ix,iy,iz,ip,it).ne.0.0,
c     >               lim5(ix,iy-ny3d-1,iz,ip,it).ne.0.0
c
c                   write(6,'(a,6i8,3(1x,g12.5),1x,2l5)')
c     >                  'LIM5:',ix,-iy,ny3d-iy+1,iz,ip,it,
c     >               fact,lim5(ix,-iy,iz,ip,it),
c     >               lim5(ix,ny3d-iy+1,iz,ip,it),
c     >               lim5(ix,-iy,iz,ip,it).ne.0.0,
c     >               lim5(ix,ny3d-iy+1,iz,ip,it).ne.0.0
c
c               endif
                
                LIM5(IX, IY,IZ,IP,IT)=FACT * LIM5(IX, IY,IZ,IP,IT)                
                LIM5(IX,-IY,IZ,IP,IT)=FACT * LIM5(IX,-IY,IZ,IP,IT)                

              !     jdemod - only scale by poloidal bin width when poloidal transport
              !          is active
              if (cdpol.ne.0.0) then 
                 LIM5(IX, IY,IZ,IP,IT)=LIM5(IX, IY,IZ,IP,IT)/pwids(ip)
                 LIM5(IX,-IY,IZ,IP,IT)=LIM5(IX,-IY,IZ,IP,IT)/pwids(ip)                
              endif              
 4500      CONTINUE                                                             
 4510     CONTINUE                                                              
 4520    CONTINUE                                                               
 4530   CONTINUE                                                                
 4540  CONTINUE                                                                 

         !if (ANY(lim5.ne.0.0)) then
         !   write(0,*) 'LIM5B - non-zero elements found'
         !endif

      ! print time dependence

      !if (imode.ne.2) then 
      ! output file
      !   outunit = 6
      !open(outunit,file='density.nt',form='formatted')

      ! only outputing max charge state for now
      !do it = 1,nts
      !   if (ANY(lim5(:,:,:,:,it).ne.0.0)) then
      !      write(0,*) 'LIM5 ',it,' non-zero elements found'
      !   endif
      !do iz = nizs,nizs
      !   write(outunit,'(a)') ' '
      !      write(outunit,'(a,i8,2(a,g12.5))') ' DENSITY IZ= ',iz,
      !>            ' TIME= ',ctimes(it,iz) * qtim
      !   write(outunit,'(1000(1x,g12.5))') 0.0,0.0,(ywids(iy),iy=1,nys)
      !   write(outunit,'(1000(1x,g12.5))') 0.0,0.0,(youts(iy),iy=1,nys)
      !   do ix = 1,nxs
      !      write(outunit,'(1000(1x,g12.5))') xwids(ix),xouts(ix),
      !>                  (lim5(ix,iy,iz,0,it),iy=1,nys)
      !   end do
      !end do
      !end do
      !endif

      ENDIF                                                                     

      
      WRITE(6,*) '5:'
C                                                                               
c slmod begin - total volume
      DO IY = -NYS, NYS
        DO IX = 1, NXS
          IF (IY.NE.0) TOTVOL = TOTVOL + YWIDS(ABS(IY)) * XWIDS(IX)
        ENDDO
      ENDDO
c slmod end
c
C===================== DDLIMS AND DDLIM3 ARRAYS ========================        
C                                                                               
C---- DEAL WITH NEUTRALS FIRST:  SCALE BY BIN WIDTHS, NEUT TIMESTEP             
C---- AND NO. OF NEUTRALS LAUNCHED;                                             
C---- THEN DEAL WITH IONS:  SCALE BY BIN WIDTHS,                                
C---- LIM3 TIMESTEP AND NO. OF IONS FOLLOWED.                                   
C                                                                               
      DO 4630 IZ = -1, NIZS                                                     
       DO 4620 IX = 1, NXS                                                      
        DO 4610 IY = 1, NYS                                                     
          DACT = DBLE (FACTB(IZ) /                                              
     >           (XWIDS(IX) * XCYLS(IX) * YWIDS(IY) * DELPS(IX,IY)))            
          DDLIMS(IX, IY,IZ) = DACT * DDLIMS(IX, IY,IZ)                          
          DDLIMS(IX,-IY,IZ) = DACT * DDLIMS(IX,-IY,IZ)                          
          IF (IY.LE.NY3D) THEN                                                  
            DO 4600 IP = -MAXNPS, MAXNPS                                        
              DDLIM3(IX, IY,IZ,IP)=DACT/pwids(ip) * DDLIM3(IX, IY,IZ,IP)                
              DDLIM3(IX,-IY,IZ,IP)=DACT/pwids(ip) * DDLIM3(IX,-IY,IZ,IP)                
c slmod tmp
              IF (DEBUGL) WRITE(79,'(4i8,10(1x,g12.5))') 
     +          IX,IY,IZ,IP,DDLIM3(IX,IY,IZ,IP),DACT,XWIDS(IX),
     +          YWIDS(IY),PWIDS(IP),XCYLS(IX),DELPS(IX,IY)
c              write(6,'(4i8,10(1x,g12.5))') 
c     +          ix,iy,iz,ip,ddlim3(ix,iy,iz,ip),dact,xwids(ix),
c     +          ywids(iy),pwids(ip),xcyls(ix),delps(ix,iy)
c slmod end
 4600       CONTINUE                                                            
          ENDIF                                                                 
 4610   CONTINUE                                                                
 4620  CONTINUE                                                                 
 4630 CONTINUE                                                                  
      WRITE(6,*) '6:'
c slmod begin - density integration over all space
      IF (BIG)THEN
        DO IZ = 1, NIZS              
        
          TOTDEN(IZ) = 0.0

          DO IX = 1, NXS
            DO IY = 1, NYS
              DO IP = -MAXNPS, MAXNPS
                TOTDEN(IZ) = TOTDEN(IZ) + 
     +                DDLIM3(IX, IY,IZ,IP) * XWIDS(IX) * YWIDS(IY)
     +                     * PWIDS(IP)           
                TOTDEN(IZ) = TOTDEN(IZ) + 
     +            DDLIM3(IX,-IY,IZ,IP) * XWIDS(IX) * YWIDS(IY)
     +                     * PWIDS(IP)
             ENDDO
            ENDDO
          ENDDO
    
          TOTDEN(IZ) = TOTDEN(IZ) / TOTVOL

        ENDDO
      ENDIF
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C      ANALYTIC EXTENSION NOTE 282 ... APPLIES TO DDLIMS ARRAY ONLY             
C      (BUT THIS THEN AFFECTS POWLS, LINES AND PLRPS ARRAYS BELOW).             
C-----------------------------------------------------------------------        
C                                                                               
       IC = IPOS (CANAL, XS, NXS-1)                                             
       IF (IC.LT.NXS) THEN                                                      
C                                                                               
         CALL RZERO (SAVES, MAXNXS*(MAXIZS+4))                                  
         DO 1890 IZ = -1, NIZS                                                  
           DO 1880 IX = 1, NXS                                                  
             DSUM1 = 0.0D0                                                      
             DO 1870 IY = 1, NYS/2                                              
               DSUM1 = DSUM1 + DBLE(YWIDS(IY)) *                                
     >                 (DDLIMS(IX,IY,IZ) + DDLIMS(IX,-IY,IZ))                   
 1870        CONTINUE                                                           
             SAVES(IX,IZ) = SNGL(DSUM1)                                         
             IF (IZ.GT.0) SAVES(IX,NIZS+1)=SAVES(IX,NIZS+1)+SNGL(DSUM1)         
 1880      CONTINUE                                                             
 1890    CONTINUE                                                               
C                                                                               
         DSUM1 = 0.0D0                                                          
         DO 1910 IZ = 1, NIZS                                                   
           DO 1900 IY = 1, NYS                                                  
             DSUM1 = DSUM1 + DBLE(YWIDS(IY)) *                                  
     >               (DDLIMS(IC,IY,IZ) + DDLIMS(IC,-IY,IZ))                     
 1900      CONTINUE                                                             
 1910    CONTINUE                                                               
         DSUM1 = DSUM1 / (2.0D0 * DWOL) * DBLE (EXP (CVIN *                     
     >         ((CA-XOUTS(IC)) * (CA-XOUTS(IC)) / (CA*CA)) - 1.0))              
         WRITE (6,'('' LIM3: ANLY EXT IC,DSUM1'',I5,G12.4)') IC,DSUM1           
C                                                                               
         DO 1940 IX = IC+1, NXS                                                 
           DO 1930 IY = -NYS, NYS                                               
             DO 1920 IZ = 1, NIZS-1                                             
               DDLIMS(IX,IY,IZ) = 0.0D0                                         
 1920        CONTINUE                                                           
             DDLIMS(IX,IY,NIZS) = DSUM1 * DBLE (EXP (CVIN *                     
     >         (1.0 - (CA-XOUTS(IX)) * (CA-XOUTS(IX)) / (CA*CA))))              
 1930      CONTINUE                                                             
 1940    CONTINUE                                                               
       ENDIF                                                                    
       WRITE(6,*) '7:'
C                                                                               
C-----------------------------------------------------------------------        
C     SCALE DEPOSITION AND NET EROSION QUANTITIES                               
C-----------------------------------------------------------------------        
C                                                                               
C---- SCALE DEPS ARRAY BY DIVIDING THROUGH BY BIN-WIDTHS AND THE TOTAL          
C---- NUMBER OF ABSORBED PARTICLES AT Y=0,2L,-2L                                
C---- SET DEPS(,,3) TO TOTAL FOR BOTH SIDES OF Y=0                              
C                                                                               
C---- NOTE : THE FOLLOWING IS SCALED BY DIVIDING BY THE TOTAL 
C---- NUMBER OF NEUTRALS, HOWEVER FOR ION INJECTION CASES THE 
C---- TOTAL NUMBER OF NEUTRALS IS ZERO AND THE FACTA(0) 
C---- VALUE IS SET TO ZERO. IN ORDER TO KEEP THE DEPOSITION DATA
C---- AVAILABLE FOR PLOTTING, TEST IF TNEUT=0 AND THEN SET THE 
C---- MULTIPLYING FACTOR TO 1.0 FOR THESE CASES. THIS COULD BE
C---- CHANGED TO THE TOTAL NUMBER OF IONS INJECTED.
C
      FACTDEPS = FACTA(0)
      IF (FACTDEPS.EQ.0.0) FACTDEPS=1.0
      IF (NIZS.GT.0) THEN                                                       
        DO 881 IZ = 1, NIZS                                                     
          DO 880 IX = 1, NXS                                                    
            DEPS(IX,IZ,1) = DEPS(IX,IZ,1) / XWIDS(IX) * FACTDEPS                
            DEPS(IX,IZ,2) = DEPS(IX,IZ,2) / XWIDS(IX) * FACTDEPS                
            DEPS(IX,IZ,3) = DEPS(IX,IZ,1) + DEPS(IX,IZ,2)                       
            DEPS(IX,NIZS+1,1) = DEPS(IX,NIZS+1,1) + DEPS(IX,IZ,1)
            DEPS(IX,NIZS+1,2) = DEPS(IX,NIZS+1,2) + DEPS(IX,IZ,2)
            DEPS(IX,NIZS+1,3) = DEPS(IX,NIZS+1,3) + DEPS(IX,IZ,3)
  880     CONTINUE                                                              
  881   CONTINUE                                                                
      ENDIF                                                                     
      WRITE(6,*) '8:'
C                                                                               
C---- SCALE REMOVAL RATES AS STORED IN NEROXS (,2) AND (,3) LOCATIONS           
C---- BY BIN WIDTHS AND TO "1 ATOM LAUNCHED"                                    
C---- CALCULATE TOTAL DEPOSITION, NET EROSION AND NENNL (NOTE 289)              
C---- CALCULATE SUM OVER BOTH SIDES OF Y=0 AND STORE THIS ALSO                  
C                                                                               
      FACT = 0.0                                                                
      IF (TDEP.GT.0.0) FACT = TNEUT / TDEP                                      
      DO 883 J = 1, 2                                                           
       DO 883 IX = 1, NXS                                                       
        NEROXS(IX,1,J) =-NEROXS(IX,1,J) / XWIDS(IX) * FACTA(0)                  
c       jdemod - change normalization of primary removal to TNEUT instead of RNEUT1
c        NEROXS(IX,2,J) = NEROXS(IX,2,J) / XWIDS(IX) * FACTA(-1)                 
        NEROXS(IX,2,J) = NEROXS(IX,2,J) / XWIDS(IX) * FACTA(0)                 
        NEROXS(IX,3,J) = NEROXS(IX,3,J) / XWIDS(IX) * FACTA(0)                  
        NEROXS(IX,4,J) = NEROXS(IX,1,J) + NEROXS(IX,3,J)                        
        NEROXS(IX,5,J) = FACT * NEROXS(IX,1,J) + NEROXS(IX,3,J)                 
  883 CONTINUE                                                                  
      DO 884 IX = 1, NXS                                                        
       DO 884 II = 1, 5                                                         
        NEROXS(IX,II,3) = NEROXS(IX,II,1) + NEROXS(IX,II,2)                     
  884 CONTINUE                                                                  
      WRITE(6,*) '9:'
                                                                               
      ! Scale removal rates in neroys, nerods arrays ...                                                                                                          

      write(6,'(a)') 'TOTAL DEPOSITION:'
      write(6,'(101(1x,g12.5))') 
     >  0.0, 0.0, (ps(ip) - pwids(ip) * 0.5, ip = -maxnps, maxnps)
      do io = 1,maxos
        write(6,'(101(1x,g12.5))') odwids(io), odouts(io), 
     >    (-nerods3(io,ip,1), ip = -maxnps, maxnps)
      end do   

      write(6,*) 'NER:', facta(0), facta(-1)
      do 885 io = 1, maxos                                                      
        if (oywids(io).gt.0.0) then
          neroys(io,1) =-neroys(io,1) / oywids(io) * facta(0)     
                        
          ! jdemod - change normalization of primary removal to tneut 
          ! instead of rneut1
          !neroys(io,2) = neroys(io,2) / oywids(io) * facta(-1)                    
          neroys(io,2) = neroys(io,2) / oywids(io) * facta(0)                    
          neroys(io,3) = neroys(io,3) / oywids(io) * facta(0)                     
        endif
        neroys(io,4) = neroys(io,1) + neroys(io,3)                              
        neroys(io,5) = fact * neroys(io,1) + neroys(io,3)  
                            
        if (odwids(io).gt.0.0) then
          nerods(io,1) =-nerods(io,1) / odwids(io) * facta(0)       
                      
          ! jdemod - change normalization of primary removal to tneut 
          ! instead of rneut1
          !nerods(io,2) = nerods(io,2) / odwids(io) * facta(-1)                    
          nerods(io,2) = nerods(io,2) / odwids(io) * facta(0)                    
          nerods(io,3) = nerods(io,3) / odwids(io) * facta(0)                     

          ! Need to scale by the 3D bin width as well taking into 
          ! account any limiter poloidal extent.
          do ip = -maxnps,maxnps

            ! jdemod - Implicitly assume that the limiter poloidal 
            ! extents have been chosen to coincide with pbin boundaries. 
            ! Doesn't make much sense otherwise and this avoids issues 
            ! with the new p bin options.
            local_pwid = pwids(ip)
            nerods3(io,ip,1) =-nerods3(io,ip,1) / odwids(io) / 
     >        local_pwid  * facta(0)                     
           
            ! jdemod - Change normalization of primary removal to tneut 
            ! instead of rneut1.
!            nerods3(io,ip,2) = nerods3(io,ip,2) / odwids(io) / 
!     >        local_pwid * facta(-1)                    
            nerods3(io,ip,2) = nerods3(io,ip,2) / odwids(io) / 
     >        local_pwid * facta(0)                    
            nerods3(io,ip,3) = nerods3(io,ip,3) / odwids(io) / 
     >        local_pwid * facta(0)                     
          end do
        endif
        nerods(io,4) = nerods(io,1) + nerods(io,3)                              
        nerods(io,5) = fact * nerods(io,1) + nerods(io,3)                       
        do ip = -maxnps,maxnps
          nerods3(io,ip,4) = nerods3(io,ip,1) + nerods3(io,ip,3)                              
          nerods3(io,ip,5) = fact * nerods3(io,ip,1) + nerods3(io,ip,3)                       
        end do
  885 continue                                                                  

C                                                                               
C---- SCALE WALL DEPOSITION ARRAY TO "1 ATOM LAUNCHED"                          
C---- SET SECONDARY NEUTRALS LINE AND TOTALS LINE.                              
C                                                                               
      DO 889 IY = -NYS, NYS                                                     
        IF (IY.EQ.0) GOTO 889                                                   
        WALLS(IY,-2) = WALLS(IY,0) - WALLS(IY,-1)                               
        DO 888 IZ = -2, NIZS                                                    
          WALLS(IY,IZ) = WALLS(IY,IZ) / YWIDS(IABS(IY)) * FACTA(0)              
          IF (IZ.GT.0) WALLS(IY,NIZS+1)= WALLS(IY,NIZS+1) + WALLS(IY,IZ)        
  888   CONTINUE                                                                
  889 CONTINUE                                                                  
      WRITE(6,*) '11:'
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE RADIATIVE POWER LOSS AND LINE RADIATION                         
C-----------------------------------------------------------------------        
C                                                                               
C---- ZERO ARRAYS.                                                              
C                                                                               
      CALL RZERO (POWLS,  MAXNXS*(2*MAXNYS+1)*(MAXIZS+2))                       
      CALL RZERO (LINES,  MAXNXS*(2*MAXNYS+1)*(MAXIZS+2))                       
C     CALL RZERO (POWL3,  MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
C     CALL RZERO (LINE3,  MAXNXS*(2*MAXY3D+1)*(MAXIZS+2)*(2*MAXNPS+1))          
C                                                                               
C-----------------------------------------------------------------------        
C   FOR 
C-----------------------------------------------------------------------        
C                                                                               
      if (cdatopt.eq.0) then 

      DO 990 IY = -NYS, NYS                                                     
        IF (IY.EQ.0) GOTO 990                                                   
        JY = IABS (IY)                                                          
C                                                                               
C---- CALCULATE PLASMA TEMPERATURE AND ELECTRON DENSITY AT MID POINTS           
C---- OF EACH X BIN.  NOTE THIS INVOLVES A CONVERSION FROM THE REGULAR          
C---- SPACED QXS MESH TO THE USER SUPPLIED XS MESH; AND A CONVERSION            
C---- FROM M**3 TO CM**3 FOR NOCORONA.                                          
C---- THE IQX --> IX INDICES ARE TAKEN FROM COMMON /COMXYT/                     
C                                                                               

      DO 900 IX = 1, NXS                                                        
        PTES(IX) = CTEMBS(IX,IY)                                                
        PNES(IX) = CRNBS(IX,IY) * 1.E-6 * REAL (CIZB)                           
  900 CONTINUE                                                                  
C                                                                               
C                                                                               
C------ CALCULATE BIN VOLUMES CM**3 (ASSUME 1 METRE IN THIRD DIMENSION)         
C------ TRANSFER IMPURITY ION DENSITIES TO PNZS ARRAY FOR NOCORONA              
C------ CONVERTING TO CM**-3                                                    
C                                                                               
        DO 920 IX = 1, NXS                                                      
          PDVOLS(IX) = 1.0E6 * XWIDS(IX) * YWIDS(JY)                            
          DO 910 IZ = 0, NIZS                                                   
            PNZS(IZ+1,1,IX) = 1.0E-6 * SNGL(DDLIMS(IX,IY,IZ))                   
  910     CONTINUE                                                              
  920   CONTINUE                                                                
C                                                                               
C------ CALL ROUTINE FROM NOCORONA PACKAGE TO CALCULATE RADIATIVE               
C------ POWER LOSS (W CM**3) AND LINE RADIATION LOSS (W CM**3)                  
C------ THESE ARE CONVERTED TO UNITS (W M**3)                                   
C------ VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY              
C------ GOES OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.         
C                                                                               
        CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,NXS)                          
C                                                                               
C------ COPY INTO "POWLS,LINES" ARRAY:  POWER LOSS; LINE RADIATION LOSS         
C------ WE KEEP TRACK OF GRAND TOTALS FOR POWER LOSS, POWER LOSS IN SOL,        
C------ LINE RADIATION LOSS & LINE RADIATION IN SOL IN DTOTS ARRAY BELOW        
C                                                                               
        DO 940 IX = 1, NXS                                                      
          DO 930 IZ = 0, NIZS                                                   
            POWLS(IX,IY,IZ) = MAX (0.0, PRADIS(1,IZ+1,1,IX)*1.0E6)              
            LINES(IX,IY,IZ) = MAX (0.0, PRADIS(3,IZ+1,1,IX)*1.0E6)              
c           write(6,'(a,3i6,8(1x,g12.5))') 'POWLS:',ix,iy,iz,
c    >         pradis(1,iz+1,1,ix)*1e6,ddlims(ix,iy,iz),
c    >         ptes(ix),pnes(ix),pnzs(iz+1,1,ix),pdvols(ix)

  930     CONTINUE                                                              
  940   CONTINUE                                                                
C                                                                               
C------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1) LOCATIONS             
C                                                                               
        DO 960 IX = 1, NXS                                                      
          IF (DDLIMS(IX,IY,0).LE.0.0D0) THEN                                    
            POWLS(IX,IY,-1) = 0.0                                               
            LINES(IX,IY,-1) = 0.0                                               
          ELSE                                                                  
            POWLS(IX,IY,-1) = POWLS(IX,IY,0) / SNGL(DDLIMS(IX,IY,0)) *          
     >                                         SNGL(DDLIMS(IX,IY,-1))           
            LINES(IX,IY,-1) = LINES(IX,IY,0) / SNGL(DDLIMS(IX,IY,0)) *          
     >                                         SNGL(DDLIMS(IX,IY,-1))           
          ENDIF                                                                 
  960   CONTINUE                                                                
C                                                                               
C------ EXTRA SECTION FOR 3D ARRAYS POWL3 AND LINE3 ...                         
C                                                                               
C
C       MEMORY RESTRICTIONS ON THE CRAY MAKE IT IMPOSSIBLE TO RUN 3D
C       CASES AND MAINTAIN ALL OF THE 3D DATA ARRAYS
C       THE ARRAYS LINE3, POWL3 ARE ALL JUST CALCULATED 
C       AND THEN WRITTEN TO DISK FOR LATER ANALYSIS
C       TO SAVE STORAGE AT THE COST OF ADDITIONAL/INEFFICIENT 
C       COMPUTATION THE CALCULATION OF THESE QUANTITIES 
C       HAS BEEN MOVED TO THE DMPOUT ROUTINE WHERE THE TIZ3 AND LIM5 
C       ARRAY STORAGE ARE REUSED TO ALLOW THE CALCULATIONS TO PROCEED  
C
C       DAVID ELDER , JAN 29 , 1990 
C
C
C       IF (JY.LE.NY3D) THEN                                                    
C         DO 985 IP = -MAXNPS, MAXNPS                                           
C           DO 970 IX = 1, NXS                                                  
C             DO 970 IZ = 0, NIZS                                               
C               PNZS(IZ+1,1,IX) = 1.0E-6 * SNGL(DDLIM3(IX,IY,IZ,IP))            
C 970       CONTINUE                                                            
C           CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,NXS)                      
C           DO 975 IX = 1, NXS                                                  
C             DO 975 IZ = 0, NIZS                                               
C               POWL3(IX,IY,IZ,IP)=MAX(0.0, PRADIS(1,IZ+1,1,IX)*1.0E6)          
C               LINE3(IX,IY,IZ,IP)=MAX(0.0, PRADIS(3,IZ+1,1,IX)*1.0E6)          
C 975       CONTINUE                                                            
C           DO 980 IX = 1, NXS                                                  
C             IF (DDLIM3(IX,IY,0,IP).LE.0.0D0) THEN                             
C               POWL3(IX,IY,-1,IP) = 0.0                                        
C               LINE3(IX,IY,-1,IP) = 0.0                                        
C             ELSE                                                              
C               POWL3(IX,IY,-1,IP) = POWL3(IX,IY,0,IP) /                        
C    >            SNGL(DDLIM3(IX,IY,0,IP)) * SNGL(DDLIM3(IX,IY,-1,IP))          
C               LINE3(IX,IY,-1,IP) = LINE3(IX,IY,0,IP) /                        
C    >            SNGL(DDLIM3(IX,IY,0,IP)) * SNGL(DDLIM3(IX,IY,-1,IP))          
C             ENDIF                                                             
C 980       CONTINUE                                                            
C 985     CONTINUE                                                              
C       ENDIF                                                                   


  990 CONTINUE                                                                  


c
c     USE ADAS DATA
c

      elseif (cdatopt.eq.1) then 



C
      DO 1190 IY = -NYS,NYS
C
C---- LOAD POWER DATA ONE RING AT A TIME.
C
        IF (IY.EQ.0) GOTO 1190                                                   
        JY = IABS (IY)                                                          
C                                                                               


        DO 1100 IX = 1, NXS   
          PTESA(IX) = CTEMBS(IX,IY)
          PNESA(IX) = CRNBS(IX,IY) * RIZB
          PNBS(IX) =  CRNBS(IX,IY)
c
c         Set hydrogen density to zero for now - not available in LIM
c          PNHS(IX) = pinaton(ik,ir)
          pnhs(ix) = 0.0

 1100   CONTINUE
        DO 1120 IX = 1, NXS
          DO 1110 IZ = 0, NIZS
            PNZSA(IX,IZ) = SNGL(DDLIMS(IX,IY,IZ))
 1110     CONTINUE
 1120   CONTINUE
C
C------ GET POWER LOSS FROM ADAS DATA FILES. LOAD TOTAL LINE RADIATION
C------ INTO LINES AND ADD RECOMBINATION AND BREMSSTRAHLUNG POWER TO
C------ GET TOTAL RADIATIVE LOSSES
C
        write(year,'(i2.2)') iyearz
        call xxuid(useridz)
c        YEAR = '89'
c        YEARDF = '89'
        ICLASS = 5
        MIZS = MIN(CION-1,NIZS)
        DO 1130 IZ = 0, MIZS
          CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NXS,PTESA,PNESA,
     +                PCOEF(1,IZ+1))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ+1))
          DO 1135 IX = 1, NXS
            LINES(IX,IY,IZ) = PCOEF(IX,IZ+1)*PNESA(IX)*PNZSA(IX,IZ)
            POWLS(IX,IY,IZ) = LINES(IX,IY,IZ)
c            write (6,'(a,3i5,3g16.8)') 'Debug DIV:',ir,ik,iz,
c     >              pcoef(ik,iz+1),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      DIV:',
c     >            lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
 1135     CONTINUE
 1130   CONTINUE
        ICLASS = 4
        MIZS = MIN(CION,NIZS)
        DO 1140 IZ = 1, MIZS
          CALL ADASRD(YEAR,CION,IZ,ICLASS,NXS,PTESA,PNESA,
     +                PCOEF(1,IZ))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ))
          DO 1145 IX = 1, NXS
            POWLS(IX,IY,IZ) = POWLS(IX,IY,IZ)
     +                        + PCOEF(IX,IZ)*PNESA(IX)*PNZSA(IX,IZ)
 1145     CONTINUE
 1140   CONTINUE
C
C------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1)
C
        DO 1160 IX = 1, NXS
          IF (DDLIMS(IX,IY,0).LE.0.0) THEN
            POWLS(IX,IY,-1) = 0.0
            LINES(IX,IY,-1) = 0.0
          ELSE
            POWLS(IX,IY,-1) = POWLS(IX,IY,0) *
     +                        SNGL(DDLIMS(IX,IY,-1) / DDLIMS(IX,IY,0))
            LINES(IX,IY,-1) = LINES(IX,IY,0) *
     +                        SNGL(DDLIMS(IX,IY,-1) / DDLIMS(IX,IY,0))
c
c            write (6,'(a,3i5,3g16.8)') 'Debug POW:',ir,ik,iz,
c     >              pcoef(ik,iz),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      POW:',
c     >           lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
          ENDIF
 1160   CONTINUE


 1190 CONTINUE


      endif


      call pr_trace('LIM3','Before ABSFAC Calculation')
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE Z EFFECTIVE ETC OVER "NEAR" REGION... NOTE 107                  
c      WRITE (datunit,'(1X,A,I7,4X,I7)') NAME,I1,I2
c
C-----------------------------------------------------------------------        
C                                                                               
C---- THREE QUANTITIES PER BIN:  NIE,   ZB.NBT(X),   ZEFF                       
C---- THE AVERAGE VALUES OVER THE "NEAR" REGION ARE SUMMED IN THE DTOTS         
C---- ARRAY BELOW, JUST OVER THE REGION 0:CXNEAR, -CYNEAR:CYNEAR.               
C---- SEE ALSO NOTES 139, 221, 225, 288, 293                                    
C                                                                               
c
c     jdemod
c
c     Change Default scaling factor to 1.0 so that data is scaled to 
c     1 particle/s if a better scaling factor is not available     
c
c      DEFACT = 0.0D0                                                            

      DEFACT = 1.0D0                                                            
c
c      jdemod - removed TATIZ/TNEUT scaling of the absolute factor
c
c      IF (TNEUT.GT.0.0) DEFACT = DBLE (2.0*GTOT1*YEFF*CSEF*TATIZ/TNEUT)         
      IF (TNEUT.GT.0.0) DEFACT = DBLE (2.0*GTOT1*YEFF*CSEF)         
c
c     Assign DEFACT to ABSFAC which is in the common include file comtor 
c     This is for compatibility with some DIVIMP code.
c
      ABSFAC = DEFACT
c
c     Calculation of scaling factor:
c
      call prb
      call prc(' CALCULATION OF "ABSOLUTE" FACTOR:')
      call prc(' FORMULA USED: ABSFAC ='//
     >         ' 2.0*GTOT1*YEFF*CSEF')
c      call prc(' FORMULA USED: ABSFAC ='//
c     >         ' 2.0*GTOT1*YEFF*CSEF*TATIZ/TNEUT')
      call prr(' ABSFAC = ',real(absfac))
      call prr(' GTOT1  = ',gtot1)
      call prr(' YEFF   = ',yeff)
      call prr(' CSEF   = ',csef)
      call prc(' For Refererence: ') 
      call pri(' TATIZ  = ',nint(tatiz))
      call pri(' TNEUT  = ',nint(tneut))
c
C                                                                               
      DO 1020 IX = 1, NXS                                                       
        DO 1010 IY = -NYS, NYS                                                  
          JY = IABS(IY)                                                         
C                                                                               
C-------- SUM CLOUD DENSITIES OVER IONISATION STATES                            
C-------- CALCULATE TOTAL IMPURITY INFLUX, MULTIPLY BY SPUTTERING               
C-------- ENHANCEMENT FACTOR...                                                 
C                                                                               
          DSUM1 = 0.0D0                                                         
          DSUM2 = 0.0D0                                                         
          DO 1000 IZ = 1, NIZS                                                  
            DIZ   = DFLOAT (IZ)                                                 
            DSUM1 = DSUM1 +  DIZ  * DDLIMS(IX,IY,IZ)                            
            DSUM2 = DSUM2 +DIZ*DIZ* DDLIMS(IX,IY,IZ)                            
 1000     CONTINUE                                                              
          DSUM1 = DSUM1 * DEFACT                                                
          DSUM2 = DSUM2 * DEFACT                                                
C                                                                               
C-------- CALCULATE AND STORE THE 3 ZEFFS RELATED QUANTITIES.                   
C-------- PRIOR TO NOTE 139, ZEFF WAS CALCULATED FROM THE FORMULA :-            
C-------- ZEFFS(IX,IY,3) = (RIZB * MAX (0.0,ZEFFS(IX,IY,2)) + SUM2) /           
C--------                  (RIZB * CRNBS(IX,IY))                                
C                                                                               
          ZEFFS(IX,IY,1) = SNGL (DSUM1)                                         
          ZEFFS(IX,IY,2) = RIZB * CRNBS(IX,IY) - ZEFFS(IX,IY,1)                 
          IF (ZEFFS(IX,IY,2).GT.0.0) THEN                                       
            ZEFFS(IX,IY,3) = (RIZB * ZEFFS(IX,IY,2) + SNGL(DSUM2)) /            
     >                       (RIZB * CRNBS(IX,IY))                              
          ELSE                                                                  
            ZEFFS(IX,IY,3) = SNGL (DSUM2/DSUM1)                                 
          ENDIF                                                                 
C                                                                               
C-------- NOTE 297.  DUPLICATE SET OF ZEFFS RESULTS BASED ON LT NOT LP+         
C-------- NOTE 299.  JUST USE SIN(THETAB) FACTOR, REMOVE LP+/LP- FACTOR         
C                                                                               
          ZEFFS(IX,IY,4) = SNGL (DSUM1) * CSINTB                                
          ZEFFS(IX,IY,5) = RIZB * CRNBS(IX,IY) - ZEFFS(IX,IY,4)                 
          IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                       
            ZEFFS(IX,IY,6) = (RIZB * ZEFFS(IX,IY,5) + SNGL(DSUM2) *             
     >        CSINTB) / (RIZB * CRNBS(IX,IY))                                   
          ELSE                                                                  
            ZEFFS(IX,IY,6) = SNGL (DSUM2/DSUM1)                                 
          ENDIF                                                                 
C         WRITE (6,'('' LIM3: IX,IY,NB'',2I5,G11.4,'' ==>'',3G11.4)')           
C    >      IX,IY,CRNBS(IX,IY),(ZEFFS(IX,IY,J), J=1,3)                          
 1010   CONTINUE                                                                
 1020 CONTINUE                                                                  
      WRITE(6,*) '13:'
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE TOTALS OVER ENTIRE PLASMA, ETC ....                             
C-----------------------------------------------------------------------        
C                                                                               
      CALL DZERO (DTOTS, 20)                                                    
      DO 4050 IY = -NYS, NYS                                                    
        IF (IY.EQ.0) GOTO 4050                                                  
        JY = IABS (IY)                                                          
        DO 4040 IX = 1, NXS                                                     
          DACT = XWIDS(IX) * XCYLS(IX) * YWIDS(JY)                              
          DTOTS(1) = DTOTS(1) + DACT                                            
          DTOTS(2) = DTOTS(2) + DACT * DBLE(CRNBS(IX,IY))                       
          DO 4030 IZ = 0, NIZS                                                  
            DTOTS(3) = DTOTS(3) + DACT * DDLIMS(IX,IY,IZ)                       
            DTOTS(4) = DTOTS(4) + DACT * DBLE(POWLS(IX,IY,IZ))                  
            DTOTS(5) = DTOTS(5) + DACT * DBLE(LINES(IX,IY,IZ))                  
            IF (XS(IX).LE.0.0) THEN                                             
              DTOTS(6) = DTOTS(6) + DACT                                        
              DTOTS(7) = DTOTS(7) + DACT * DBLE(POWLS(IX,IY,IZ))                
              DTOTS(8) = DTOTS(8) + DACT * DBLE(LINES(IX,IY,IZ))                
            ENDIF                                                               
 4030     CONTINUE                                                              
          IF (XS(IX).GT.0.0.AND.XS(IX).LE.CXNEAR.AND.YS(JY).LE.CYNEAR)          
     >    THEN                                                                  
            DTOTS(9) = DTOTS(9) + DACT                                          
            DTOTS(10)= DTOTS(10)+ DACT * DBLE(ZEFFS(IX,IY,4))                   
            DTOTS(11)= DTOTS(11)+ DACT * DBLE(ZEFFS(IX,IY,5))                   
            DTOTS(12)= DTOTS(12)+ DACT * DBLE(ZEFFS(IX,IY,6))                   
          ENDIF                                                                 
 4040   CONTINUE                                                                
 4050 CONTINUE                                                                  
      DTOTS(10) = DTOTS(10) / DTOTS(9)                                          
      DTOTS(11) = DTOTS(11) / DTOTS(9)                                          
      DTOTS(12) = DTOTS(12) / DTOTS(9)                                          
      DTOTS(13) = DTOTS(2) / DTOTS(1)                                           
      DTOTS(14) = DTOTS(4) / (DTOTS(13) * DTOTS(3))                             
      DO 4090 J = 1, 20                                                         
        STOTS(J) = SNGL(DTOTS(J))                                               
 4090 CONTINUE                                                                  
      WRITE(6,*) '14:'
C                                                                               
C-----------------------------------------------------------------------        
C                     PRINT CLOSING MESSAGE                                     
C-----------------------------------------------------------------------        
C                                                                               


      IF (NIZS.GT.0)                                                            
     >     CALL MONPRI (QTIM,FACTA(1),VFLUID,NIZS,SDTZS,SDYZS,                  
     >           STOTS,DOUTS,RIONS,CTBIN,CRMI)                                  

c
c     Print yreflection statistics if the option is active
c

      call pr_yref_stats
c
c     Print velocity diagnostic data
c
      call print_diagvel(qtim,nizs)
      

      CALL PRB                                                                  
      CALL PRI ('NUMBER OF NEUTRALS FOLLOWED   ',NINT(TNEUT))                   
      CALL PRI ('NUMBER OF IONS FOLLOWED       ',NINT(TATIZ))                   
      WRITE (6,'('' NUMBER OF NEUTRALS FOLLOWED   '',G11.4)') TNEUT             
      WRITE (6,'('' NUMBER OF IONS FOLLOWED       '',G11.4)') TATIZ             
C                                                                               
C---- FORMATS ...                                                               
C                                                                               
 9002 FORMAT(1X,I5,F9.1,12X,I4,4X,F10.6,10X,F10.5)                              
c slmod
 9003 FORMAT(1X,I5,1x,f10.1,4(1x,I6),2F12.6,2F12.5,1P,G11.3,0P,F10.4,            
     >  1P,G10.3,0P,F7.4,3(1x,I4),1X,A,:,I3,I4,F8.2)                                 
c
c 9003 FORMAT(1X,I5,F9.1,I6,I6,I4,I4,2F10.6,2F10.5,1P,G11.3,0P,F8.2,            
c     >  1P,G10.3,0P,F7.4,3I3,1X,A,:,I3,I4,F8.2)                                 
c slmod end
 9004 FORMAT(//1X,'LIM DEBUG: DIAGNOSTICS TO BE PRINTED EVERY',I6,              
     >     ' TIMESTEPS  (DELTA T =',G10.3,' SECONDS).',//)
      
 9005 FORMAT(1X,'--ION-----TIME----IQX----IQY---IX---IY------X-----',
     >  '----ALPHA',         
     >  '----------Y-----------P-------DRIFT-VEL----TEMP--PARA-DIFF',
     >  '--FRACT--IP---IT---IS--',            
     >     12('-'))
      
 9006 FORMAT(//1X,'LIM DEBUG: TRACK DIAGNOSTICS TO BE RECORDED FOR FIRST
     >',I6,' PARTICLES')              
 9010 FORMAT(/5X,'                                        IONS SURVIV',         
     >                                          'ING     ',                     
     >       /5X,'SPLITTING    IONS CROSSING INWARDS  AFTER CROSSING ',         
     >                                          'OUTWARDS',                     
     >       /5X,'  PLANE       WEIGHTED  UNWEIGHTED    WEIGHTED  UNW',         
     >                                          'EIGHTED ')                     
 9011 FORMAT(1X,I3,'  X =',F6.3,'M',I9,I11,I13,I11)                             
 9012 FORMAT(/1X,A,A,A)                                                         
 9013 FORMAT(/1X,A,I7,A)
 9020 FORMAT(' IZ',7F8.3,/,(3X,7F8.3))                                          
 9021 FORMAT(I3,1P,1X,7G8.1,/,(4X,7G8.1))                                       
 9022 FORMAT(I3,7F9.1,/,(3X,7F9.1))                                             
c
c slmod begin
      WRITE(0,*) 'Done  LIM3'

      TITL2 = TITLE
      MIZS  = NIZS
      MIMPS = NIMPS
      MATIZ = NATIZ

      CALL OutputDiag
      CALL GetProfiles
c slmod end
      RETURN                                                                    
      END                                                                       
c     
c     
c     
      subroutine check_y_boundary(ip,ix,cx,y,oldy,absy,svy,alpha,ctwol,
     >                            sputy,ciz,debugl,ierr,vary_2d_bound)
      use error_handling
      use yreflection
      
      ! This routine checks to see if the particle has reached the 
      ! Y-bounds of the modeling space and then adjusts the Y coordinate 
      ! of the particle appropriately. In addition, if the Y-axis mirror 
      ! option is in use this code checks for reflections from the 
      ! mirrors at the specified Y values. 

      implicit none
      real :: cx,y,oldy,ctwol,absy,svy,alpha,sputy,tmp_oldy
      logical :: debugl
      integer :: ierr,ciz, ip, ix, vary_2d_bound

      tmp_oldy = oldy

      ierr = 0

      ! Check for crossing Y-absorbing surface before the y-coordinate 
      ! are updated.
      if (yabsorb_opt.ne.0) then 
         
        ! Call correct routine based off using a fully customizable 2D 
        ! wall or not.
        if (vary_2d_bound.eq.1) then
          call check_y_absorption_2d(ip, ix, cx, y, oldy, sputy, ciz, 
     >      ierr)
        else
          call check_y_absorption(cx,y,oldy,sputy,ciz,ierr)
        endif

        if (ierr.eq.1) then 
          ! Particle absorbed - exit tracking loop - y absorption
          ierr = 2 
          return
        endif 
      endif

      IF (Y.LE.-CTWOL) THEN                                           

 401     continue
         Y   = y + 2.0 * ctwol
         tmp_oldy = tmp_oldy + 2.0 * ctwol
         IF (Y.LE.-CTWOL) GOTO 401                                     
    
c        jdemod 
c     
c        Need to make sure that a particle 
c        does not enter a reflected region
c        inside the confined plasma.
c     
         if (yreflection_opt.ne.0) then 
            if (check_reflected_region(y)) then 
               write(6,'(a,5(1x,g18.10))') 
     >              'REFLECTION ERROR INBOARD < CTWOL:',alpha,y,oldy

               call check_reflection(cx,y,tmp_oldy,svy,sputy,
     >                               2,debugl,ierr)

               if (y.lt.-ctwol) then 
                  y = y+2.0*ctwol
               elseif (y.gt.ctwol) then 
                  y = y-2.0*ctwol
               endif
c     
               if (check_reflected_region(y)) then 
                  CALL errmsg('LIM3: ION INBOARD:',
     >                 'ION HAS ENTERED MIRROR BOUNDED REGION')
               endif

            endif  
         endif

      ELSEIF (Y.GE.CTWOL) THEN                                        

 402     continue
         Y   = Y - 2.0 * ctwol
         tmp_oldy = tmp_oldy - 2.0 * ctwol
         IF (Y.GE.CTWOL) GOTO 402                                      

         if (yreflection_opt.ne.0) then 
            if (check_reflected_region(y)) then 
               write(6,'(a,5(1x,g18.10))') 
     >              'REFLECTION ERROR INBOARD > CTWOL:',alpha,y,oldy

               call check_reflection(cx,y,tmp_oldy,svy,sputy,
     >                               2,debugl,ierr)

c     
               if (y.lt.-ctwol) then 
                  y = y+2.0*ctwol
               elseif (y.gt.ctwol) then 
                  y = y-2.0*ctwol
               endif
c     
               if (check_reflected_region(y)) then 
                  CALL errmsg('LIM3: ION INBOARD:',
     >                 'ION HAS ENTERED MIRROR BOUNDED REGION')
               endif
               
            endif
         endif

      ENDIF                                                           


      ABSY = ABS (Y)                                                

      return
      end
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE HIT (OLDALP,ALPHA,OLDY,Y,CIAB,IQX,IX,IOY,IOD,XM,YM)            
      use mod_params
      use error_handling
      use mod_comt2
      use mod_comnet
      use mod_comtor
      use mod_comxyt
      IMPLICIT none                                                    
      REAL    OLDALP,ALPHA,OLDY,Y,XM,YM                                         
      INTEGER CIAB,IQX,IX,IOY,IOD                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  HIT:  WE KNOW ROUGHLY WHERE AN ION HAS HIT THE LIMITER.  THIS    *        
C  *  ROUTINE BACKTRACKS ALONG THE FINAL TRAJECTORY TO DETERMINE MORE  *        
C  *  PRECISELY WHERE THE HIT OCCURED.  IT RETURNS THE INDICES IQX,    *        
C  *  IX, IOY AND IOD INDICATING THE POSITION WHERE THE ION STRUCK     *        
C  *  THE SURFACE.  THE ROUTINE RETURNS (XM,YM) THE IMPACT POINT, BUT  *        
C  *  SHOULD HAVE VERY LITTLE EFFECT EXCEPT WHERE THE CROSS FIELD      *        
C  *  DIFFUSION STEPS ARE LARGE  (IE, BIG TIMESTEPS).  WRITTEN IN      *        
C  *  RESPONSE TO NOTE 292.                                            *        
C  *                                                                   *        
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
c      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
c      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
c      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
c      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
c      INCLUDE 'comnet'                                                          
C     INCLUDE (COMNET)                                                          
      REAL    XT,XB,YT,YB,EDGE1,EDGE2,DIST1,DIST2                               
      INTEGER K,IPOS                                                            
      LOGICAL THERE                                                             
c
      real :: xt_org,xb_org
C                                                                               
c     jdemod - 
c     Collision must be on limiter - therefore set XT to 0.0 if OLDALP > 0
c      

      if (oldalp.gt.0.0) then 
         XT = 0.0
      else
         XT = oldalp
      endif
c
      XB = ALPHA                                                                
c
      xt_org = xt
      xb_org = xb
c
      YT = OLDY                                                                 
      YB = Y                                                                    
      K  = 1                                                                    
c
      IF (DEBUGL) THEN                                                          
C       WRITE (6,9001) OLDALP,OLDY                                              
C       WRITE (6,9001) ALPHA,Y,IQX,QEDGES(IQX,1),QEDGES(IQX,2),K,.TRUE.         
      ENDIF                                                                     
C                                                                               
  100 CONTINUE                                                                  

      IF (XT.GT.0.0.AND.XB.GT.0.0) THEN 
         write(error_message_data,
     >               '(a,10(a,g18.10))')
     >  'HIT CALLED WHEN OLDALP and ALPHA both greater than 0.0:',
     >               ' OLDALP =',oldalp,
     >               ' ALPHA =',alpha,
     >               ' XT    =',xt,
     >               ' XB    =',xb,
     >               ' OLDY =',oldy,
     >               ' Y =',y
         CALL errmsg('HIT:',error_message_data)


         xm = 0.0 - 1.0e-10
         xt = xm
         xb = xm
         YM = 0.5 * (YT + YB)                                                      
      else
         XM = 0.5 * (XT + XB)                                                      
         YM = 0.5 * (YT + YB)                                                      
      endif

      IQX = MAX (INT (XM*XSCALO), 1-NQXSO)                                      
C                                                                               
c
c     jdemod - changed (iqx.gt.0) to (iqx.ge.0) so that 0 values will not give a 
c              match - this seems to cause the search algorithm to walk off the 
c              end of the limiter and give an intersection with X>0
c     Change back to gt 0 when max xt = 0.0 imposed
c     
c

      IF (IQX.Gt.0) THEN                                                        
        XT    = XM                                                              
        YT    = YM                                                              
        THERE =.FALSE.                                                          
C                                                                               
      ELSE                                                                      
        CALL EDGINT (XM,IQX,1,EDGE1,DIST1)                                      
        CALL EDGINT (XM,IQX,2,EDGE2,DIST2)                                      
        IF ((CIAB.EQ. 1 .AND. YM.LE.EDGE2)        .OR.                          
     >      (CIAB.EQ. 2 .AND. YM.GE.CTWOL-EDGE1)  .OR.                          
     >      (CIAB.EQ.-1 .AND. YM.GE.-EDGE1)       .OR.                          
     >      (CIAB.EQ.-2 .AND. YM.LE.EDGE2-CTWOL)) THEN                          
          XB    = XM                                                            
          YB    = YM                                                            
          THERE =.TRUE.                                                         
        ELSE                                                                    
          XT    = XM                                                            
          YT    = YM                                                            
          THERE =.FALSE.                                                        
        ENDIF                                                                   
      ENDIF                                                                     

c      write(0,'(a,2i8,20(1x,g12.5))') 'HIT1:',ciab,iqx,
c     >                    oldalp,alpha,oldy,y,
c     >                    xm,ym,xt,yt,edge1,dist1,edge2,dist2

C                                                                               
      K = K + 1                                                                 
      IF (K.LE.100) THEN                                                        
C       IF (DEBUGL) WRITE (6,9001) XM,YM,IQX,EDGE1,EDGE2,K,THERE                
        IF ((.NOT.THERE) .OR. (K.LE.5)) GOTO 100                                
      ENDIF                                                                     
C                                                                               
      if (xm.gt.0.0) then 
c
c        correct error by causing impact at limiter tip
c
         write(error_message_data,'(a,5(a,g18.10),a)')
     >        'Calculated XM greater than 0.0:',
     >        ' XM = ',xm,' XB =',xb,' XT = ',xt,
     >        ' XB_ORG =',xb_org,' XT_ORG = ',xt_org,
     >        ' XM=MIN(XT_ORG,XB_ORG) Assigned'

         CALL errmsg('HIT:',error_message_data)
         !xm = 0.0-1.0e-10
         xm = min(xb_org,xt_org)
         IQX = MAX (INT (XM*XSCALO), 1-NQXSO)                                      
         CALL EDGINT (XM,IQX,1,EDGE1,DIST1)                                      
         CALL EDGINT (XM,IQX,2,EDGE2,DIST2)                                      
      endif
c
      IX = IPOS (XM, XS, NXS-1)                                                 
      IF (CIAB.EQ.-1 .OR. CIAB.EQ.2) THEN                                       
        YM  =-EDGE1 - 1.E-10                                                    
        IOY = IPOS (-EDGE1, OYS, MAXOS-1)                                       
        IOD = IPOS (-DIST1, ODS, MAXOS-1)                                       
      ELSE                                                                      
        YM  = EDGE2 + 1.E-10                                                    
        IOY = IPOS ( EDGE2, OYS, MAXOS-1)                                       
        IOD = IPOS ( DIST2, ODS, MAXOS-1)                                       
      ENDIF                                                                     


c      WRITE (0,'(a,3i8,l5,10(1x,g12.5))') 'HIT2:',IQX,IOY,IOD,there,
c     >          XM,YM,EDGE1,EDGE2,dist1,dist2

      RETURN                                                                    
 9001 FORMAT(1X,'HIT: XM',F10.6,' YM',F10.5,:,                                  
     >  ' IQX',I6,' EDGES',2F10.6,I5,L2)                                        
      END                                                                       
