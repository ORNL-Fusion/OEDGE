      MODULE PHOTON
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CESTIM
      USE CGEOM
      USE CGRID
      USE CLOGAU
      USE COMPRT
      USE COMSOU
      USE COMUSR
      USE COMXS
      USE COUTAU
      USE CRAND
      USE CSDVI
      USE CSPEI
      USE CSPEZ
      USE CTEXT
      USE CTRCEI
      USE CUPD
      USE CZT1

      IMPLICIT NONE

      PRIVATE

csw public funcs/subs
       PUBLIC :: ph_init, ph_energy, ph_getcoeff,
     .    ph_xsectph, ph_alloc_xsecta, ph_xsectp,
     .    ph_update_surf_spectrum, ph_post_energy, ph_alloc_xsectph, 
     .    ph_integrate, ph_lorvdw, ph_vdwprof,
     .    ph_b21, ph_philips_tnprofile
!      PUBLIC :: ph_mod_photon, ph_prepare_global,
!     .    ph_prepare_stratum, 
!     .    ph_xsecta, ph_update_spectrum,
!     .    ph_showdata, ph_finalize, 
!     .    ph_samvol_angle,ph_sm0usr,ph_sm1usr,
!     .    ph_post_stratum, ph_sam_gauss
csw public vars
      integer, public, save  :: phv_nrota, phv_nrotph, phv_muldens
      integer, public, save, allocatable :: 
     .   PHV_LGAOT(:,:,:),PHV_LGPHOT(:,:,:),
     .   PHV_NAOTI(:),PHV_NPHOTI(:),PHV_LGPRC(:,:),PHV_IESTOTPH(:,:,:),
     .   PHV_IESTOTAT(:,:,:),
     .   PHV_N1STOTPH(:,:,:),PHV_N2NDOTPH(:,:,:),
     .   PHV_N1STOTAT(:,:,:),PHV_N2NDOTAT(:,:,:),
     .   phv_xistra(:)
csw constants
      real(dp), public, save :: CLIGHT,EV2HZ,HPLNK,STEFBCON,hpcl
      real(dp), public, save :: hwvdw

csw external
      integer, external :: idez, learc1, learc2
      real(dp), external :: ranf_eirene

csw some help types
      type my_2d_realpointer
      real(dp), dimension(:,:), pointer :: p
      end type my_2d_realpointer

csw internal vars
      real(DP),SAVE,allocatable, dimension (:) ::
     .    eemax, eemin, eeint, splinemax,
     .    neweemax,neweemin

      real(DP),SAVE,allocatable,TARGET,dimension(:,:)::splinex,splinea,
     .    splineb,splinec,lorvdwnorm1,lorvdwnorm2,
     .    faddeeva2norm

      real(DP), POINTER, dimension(:) :: spa,spb,spc,spx

c input variables read from 'photon.dat'
      integer, save :: numlines,numvolspecs,numsrfspecs,
     .    iterstat,spline_ivs,interval_ivs
      integer, save, allocatable,dimension(:) ::
     .    volspecnum, volspectal,
     .    srfspecnum, srfspecsurf, iterarg,
     .    line2volspec,line2srfspec,volspec2line,srfspec2line,
     .    srfspeclne,volspeclne,
     .    phphilbulk2phil_h,phphilnivbulk,phphilnivflag,      
     .    phphilbulk2phil_l,phphilline2phil,phphilphil2bulk_h,
     .    phphilphil2bulk_l,phphilphil2line
      real(DP),save, allocatable,dimension(:) ::volspecmin,
     .    volspecmax,srfspecmin,srfspecmax
      character*20, save, allocatable, dimension(:) :: volspecfile,
     .    srfspecfile
      character*20, save, target :: iterstr(5)
      character*10, save, allocatable, dimension(:) :: phphilnivname

c spectrum tally pointer
      type(my_2d_realpointer), save, allocatable :: ptr_volspec(:),
     .                                              ptr_srfspec(:)

c dummy-args, should be removed with ph_setarg/ph_clrarg
      integer, save :: arg_ipls,arg_iaei,arg_iprt,arg_flag
      real(dp),save :: arg_e0,arg_diwl,arg_tiwl
      integer, public, save :: actual_linenum
      logical, save :: flag_adapt_intervals,flag_make_splines,
     .                 flag_planck,hermann_flag

c philips-stuff:
      real(dp), public, allocatable, save, dimension(:) :: philips_e0,
     .                                             philips_e1,
     .                                             philips_e2,
     .                                             philips_g1,
     .                                             philips_g2,
     .                                             philips_aik,
     .                                             philips_cs,
     .                                             philips_cr,
     .                                             philips_cw,
     .                                             philips_ss,
     .                                             philips_sr,
     .                                             philips_c3,
     .                                             philips_c4,
     .                                             philips_c6,
     .                                             philips_l0
      
      character*20, allocatable, save ::           philips_name(:)

      integer, save :: philips_numlines,phv_philips_numlines,
     .                 phv_philips_numnivs, phv_philips_ignd,
     .                 phv_philips_iion
      real(dp), save :: phv_philips_pressure,phv_philips_ion_threshold

c samusr-stuff
      integer, save :: phv_pnum, phv_ecnum, phv_pnummax, phv_directs,
     .    phv_coords, phv_enum, phv_cnum, phv_dnum
      real(dp), save, allocatable :: xc(:), yc(:), zc(:), vxc(:),
     .    vyc(:), vzc(:), ec(:), ecdelta(:)

c bgk-stuff
      integer, save :: bgk_num
      integer, save, allocatable :: bgk_iphot(:), bgk_iatm(:),
     .                              bgk_iprocph(:),bgk_iprocat(:)
      logical, save, allocatable :: bgk_flagr(:), bgk_flageat(:), 
     .    bgk_flageph(:)
      real(dp), save,allocatable :: bgk_tallies(:,:,:,:)

      CONTAINS

c I/O & MISC-ROUTINES
      SUBROUTINE PH_INIT(ICAL)

      real(DP) :: x,dx
      integer :: ios,i,ia,ip,j,ih,il,jh,jl,nr,ierr,jj,ipl,iat,iplh,
     .           iflag,ipll,pil,ignd,iion,nh,nl,iplh_pb, ipll_pb
      logical :: ex
      real(dp),allocatable,dimension(:,:) :: dummytgt
      integer, allocatable :: idummy(:)
      character(72) :: cline
      integer, intent(in) :: ical

      select case(ical)
csw ICAL == 0
      CASE(0)
c constants
      CLIGHT = 2.99792458D10
      EV2HZ  = 2.41798834D14
      HPLNK  = 1./EV2HZ
      STEFBCON=8.*PIA**5 / (15. * CLIGHT**3 * HPLNK**3)
c h*c [eV*cm]
      HPCL   = CLIGHT*HPLNK
c misc
      flag_adapt_intervals=.false.
      flag_make_splines=.false.

!pb      hermann_flag = .true.
      hermann_flag = .false.

c
c init_reacexp
c there are different types of processes:
c mode = -1: undefined
c      =  0: no bulk precoll partner (i.e. sp.emission)
c      =  1: 1 precoll partner, 1 secondary (ie. absorption)
c      =  2: 1 precoll partner, 2 secondaries (ie. stim.emission)
c
c updf = 0: flag for update: nothing happens
c      = 1:                  2 times wtrsig in update (stim.em.)
c      = 2:                  
      return

csw ICAL == 1
      case (1)

c allocate
      allocate(eemax(nrad), stat=ierr)
      allocate(eemin(nrad), stat=ierr)
      allocate(neweemax(nrad), stat=ierr)
      allocate(neweemin(nrad), stat=ierr)
      allocate(eeint(nrad), stat=ierr)
      spline_ivs=0
      interval_ivs=0 

c read photon inputfile 'photon.dat'
      INQUIRE (file='photon.dat', exist=ex)
      if(ex) then
         open(unit=77,file='photon.dat')

c lines-block
         read(77,FMT='(I6)',IOSTAT=ios) numlines
         
         
         read(77,FMT='(A72)',IOSTAT=ios) cline

c volume spectrum-block
         read(77,FMT='(I6)',IOSTAT=ios) numvolspecs
c         if(numvolspecs .ne. numlines) then
c            write(6,*)'PHOTON MODULE (PH_INIT): numvolspecs != numlines'
c            call exit(1)
c         endif
         allocate(volspecmin(0:numvolspecs))
         allocate(volspecmax(0:numvolspecs))
         allocate(volspecnum(0:numvolspecs))
         allocate(volspeclne(0:numvolspecs))
         allocate(volspectal(0:numvolspecs))
         allocate(volspecfile(0:numvolspecs))
         allocate(line2volspec(numlines))
         allocate(volspec2line(numvolspecs))         
         line2volspec=0
         volspec2line=0
         volspecnum=0

         do i=1,numvolspecs

            read(77,FMT='(2E15.4,3I6,1X,A20)',IOSTAT=ios) 
     .                                             volspecmin(i),
     .                                             volspecmax(i),
     .                                             volspecnum(i),
     .                                             volspeclne(i),
     .                                             volspectal(i),
     .                                             volspecfile(i)
c postprocessing in PH_INIT(2)
         enddo

         read(77,FMT='(A72)',IOSTAT=ios) cline

c surface spectrum-block
         read(77,FMT='(I6)',IOSTAT=ios) numsrfspecs
c         if(numsrfspecs .ne. numlines) then
c            write(6,*)'PHOTON MODULE (PH_INIT): numsrfspecs != numlines'
c            call exit(1)
c         endif
         allocate(srfspecmin(0:numsrfspecs))
         allocate(srfspecmax(0:numsrfspecs))
         allocate(srfspecnum(0:numsrfspecs))
         allocate(srfspeclne(0:numsrfspecs))
         allocate(srfspecsurf(0:numsrfspecs))
         allocate(srfspecfile(0:numsrfspecs))
         allocate(line2srfspec(numlines))
         allocate(srfspec2line(numsrfspecs))
         line2srfspec=0
         srfspec2line=0
         srfspecnum=0

         do i=1,numsrfspecs
            read(77,FMT='(2E15.4,3I6,1X,A20)',IOSTAT=ios) 
     .                                             srfspecmin(i),
     .                                             srfspecmax(i),
     .                                             srfspecnum(i),
     .                                             srfspeclne(i),
     .                                             srfspecsurf(i),
     .                                             srfspecfile(i)
         enddo

         read(77,FMT='(A72)',IOSTAT=ios) cline

c iteration-block
         read(77,FMT='(I6)',IOSTAT=ios)  iterstat
         allocate(iterarg(5))
         read(77,FMT='(5I6)',IOSTAT=ios) iterarg(1:5)
         do i=1,5
            read(77,FMT='(A72)',IOSTAT=ios) cline
            iterstr(i) = cline(1:20)
         enddo

         read(77,FMT='(5L1)', IOSTAT=ios)
     .       flag_adapt_intervals,
     .       flag_make_splines

         read(77,FMT='(A72)',IOSTAT=ios) cline

c bgk-block
         read(77,fmt='(I6)',iostat=ios) bgk_num
         if(bgk_num > 0) then
            allocate(bgk_iphot(bgk_num))
            allocate(bgk_iatm(bgk_num))
            allocate(bgk_iprocph(bgk_num))
            allocate(bgk_iprocat(bgk_num))
            allocate(bgk_flagr(bgk_num))
            allocate(bgk_flageat(bgk_num))
            allocate(bgk_flageph(bgk_num))
            bgk_iphot=0
            bgk_iatm=0
            bgk_iprocph=0
            bgk_iprocat=0
            bgk_flagr=.false.
            bgk_flageat=.false.
            bgk_flageph=.false.
            do i=1,bgk_num
               read(77,fmt='(4I6,1X,5L1)',iostat=ios)
     .              bgk_iphot(i), bgk_iprocph(i), bgk_iatm(i), 
     .              bgk_iprocat(i), bgk_flagr(i), bgk_flageat(i),
     .              bgk_flageph(i)
            enddo
         endif

         read(77,FMT='(A72)',IOSTAT=ios) cline

c philips block
         read(77,FMT='(I6)', IOSTAT=ios) phv_philips_numnivs
         if(phv_philips_numnivs > 0) then
            allocate(phphilnivbulk(phv_philips_numnivs))         
            allocate(phphilnivflag(phv_philips_numnivs))         
            allocate(phphilnivname(phv_philips_numnivs))         
            phphilnivbulk=0
            phphilnivflag=0
            do i=1,phv_philips_numnivs
               cline=' '
               read(77,FMT='(2I6,1X,A10)',IOSTAT=ios) ipl,ip,cline
               phphilnivbulk(i)=ipl
               phphilnivflag(i)=ip
               phphilnivname(i)=cline(1:10)
            enddo
         endif
         read(77,FMT='(I6)', IOSTAT=ios) phv_philips_numlines
         if(phv_philips_numlines>0) then
            allocate(phphilphil2line(phv_philips_numlines))
            allocate(phphilbulk2phil_h(nplsi))
            allocate(phphilbulk2phil_l(nplsi))
            allocate(phphilphil2bulk_h(phv_philips_numlines))
            allocate(phphilphil2bulk_l(phv_philips_numlines))
            allocate(phphilline2phil(numlines))
            phphilphil2line=0
            phphilbulk2phil_h=0
            phphilbulk2phil_l=0
            phphilphil2bulk_h=0
            phphilphil2bulk_l=0
            phphilline2phil=0
            do i=1,phv_philips_numlines
               read(77,FMT='(3I6)',IOSTAT=ios) il,iplh, ipll
               if(il <= 0 .or. il > numlines) then
                  write(6,*)'PHOTON MODULE (PH_INIT):'
                  write(6,*)'   philips-block: il=',il
                  call exit(1)
               endif
               phphilphil2line(i)=il
               phphilbulk2phil_h(iplh)=i
               phphilbulk2phil_l(ipll)=i

               phphilphil2bulk_h(i)=iplh
               phphilphil2bulk_l(i)=ipll

               phphilline2phil(il)=i
            enddo
         endif
         read(77,FMT='(2E15.4)', IOSTAT=ios) phv_philips_pressure,
     .                                       phv_philips_ion_threshold


c that's all
         close(77)
         write(6,*)'PHOTON MODULE (PH_INIT): photon.dat read.'
      else
!pb         write(6,*)'PHOTON MODULE (PH_INIT): cant open photon.dat, exit'
!pb         call exit(1)
         phv_philips_numlines=0   
      endif
c photon.dat read

c read PHILIPS lines.dat ?
      if(phv_philips_numlines > 0) then
         call ph_read_philips
         if(philips_numlines /= phv_philips_numlines) then
            write(6,*) 'PHOTON MODULE (PH_INIT):'
            write(6,*) '   philips_numlines /= phv_philips_numlines'
            call exit(1)
         endif
      endif
c set ground/ion-state
      ignd=0
      iion=0
      do i=1,phv_philips_numnivs
         if(phphilnivflag(i) ==  0) ignd=phphilnivbulk(i)
         if(phphilnivflag(i) == -1) iion=phphilnivbulk(i)
      enddo
      phv_philips_ignd=ignd
      phv_philips_iion=iion

c init some vars
      phv_nrota=0
      phv_nrotph=0



      return

      case(2)
c ICAL=2 
c allocate  post-processing
         allocate(ptr_volspec(0:numvolspecs))
         do i=1,numvolspecs
            if(volspeclne(i) /= 0) then
!pb avoid array bound violation
               if (volspeclne(i) > numlines) then
                  write (6,*) ' problem in defintion of volume spectra '
                  write (6,*) ' number of line ',volspeclne(i)
                  write (6,*) ' maximum number of lines available ',
     .                          numlines
                  call exit(1)
               else
                  line2volspec(volspeclne(i))=i           
                  volspec2line(i)=volspeclne(i)
                  allocate(ptr_volspec(i)%p(nrad,-1:volspecnum(i)+1),
     .                     STAT=ierr)
               end if
            else
               if(.not.associated(ptr_volspec(0)%p)) then
                  allocate(ptr_volspec(0)%p(nrad,
     .                     -1:volspecnum(i)+1),
     .                     STAT=ierr)
                  volspecmin(0) = volspecmin(i)
                  volspecmax(0) = volspecmax(i)
                  volspecnum(0) = volspecnum(i)
                  volspeclne(0) = volspeclne(i)
                  volspectal(0) = volspectal(i)
                  volspecfile(0) = volspecfile(i)
               else
                  write(6,*) 'PHOTON MODULE (PH_INIT,2):'
                  write(6,*) '    global volspec already allocated'
                  write(6,*) '    check photon.dat!!'
                  call exit(1)
               endif
            endif
         enddo

c srfspecs post-processing
         allocate(ptr_srfspec(0:numsrfspecs))
         do i=1,numsrfspecs
            if(srfspecsurf(i) > nlimps) then
               write(6,*) 'PHOTON MODULE (PH_INIT(2)):'
               write(6,*) '    srfspecsurf(',i,') > nlimps'
               call exit(1)
            endif
            if(srfspeclne(i) /= 0) then
               if (srfspeclne(i) > numlines) then
                  write (6,*) ' problem in defintion of surface spectra'
                  write (6,*) ' number of line ',srfspeclne(i)
                  write (6,*) ' maximum number of lines available ', 
     .                          numlines
                  call exit(1) 
               else 
                  line2srfspec(srfspeclne(i))=i           
                  srfspec2line(i)=srfspeclne(i)
                  allocate(ptr_srfspec(i)%p(0:nlimps,
     .                     -1:srfspecnum(i)+1),STAT=ierr)
                  ptr_srfspec(i)%p(0:nlimps,-1:srfspecnum(i)+1)=0._dp
               end if
            else
               if(.not.associated(ptr_srfspec(0)%p)) then
                  allocate(ptr_srfspec(0)%p(0:nlimps,
     .                     -1:srfspecnum(i)+1),STAT=ierr)
                  ptr_srfspec(0)%p(0:nlimps,-1:srfspecnum(i)+1)=0._dp
                  srfspecmin(0) = srfspecmin(i)
                  srfspecmax(0) = srfspecmax(i)
                  srfspecnum(0) = srfspecnum(i)
                  srfspeclne(0) = srfspeclne(i)
                  srfspecsurf(0)= srfspecsurf(i)
                  srfspecfile(0)= srfspecfile(i)
               else
                  write(6,*) 'PHOTON MODULE (PH_INIT,2):'
                  write(6,*) '    global srfspec already allocated'
                  write(6,*) '    check photon.dat!!'
                  call exit(1)
               endif
            endif
         enddo
c
         allocate(phv_xistra(0:numlines))
         phv_xistra=0

c
         if(bgk_num > 0) then
            allocate(bgk_tallies(bgk_num,0:1,3,0:nrad))
c x,0,1,nrad : iphot, flagr
c x,0,2,     : iphot, flageat
c x,0,3,     : iphot, flageph
c x,1,1,     : iatm , flagr
c x,1,2,     : iatm , flageat
c x,1,3,     : iatm , flageph
         endif
c
         return

      case default
         write(6,*) 'PHOTON MODULE (PH_INIT): ICAL=',ical,'N/A'
         call exit(1)
      end select
      RETURN
      END SUBROUTINE PH_INIT

c CROSS-SECTIONS

      SUBROUTINE PH_GETCOEFF(kkin,isp,ity,icell,iipl,idsc,fac,res) 
      IMPLICIT NONE
      integer, intent(in) :: kkin,isp,ity,icell,iipl, idsc
      real(dp), intent(out) :: fac,res
      integer :: iid,ii,nrc,n1,n2,n,i,iflag,iwarn,ivs,kk,n4,ipl,pil,
     .           iil,ignd,iptype, ipl2
      real(dp)::gam,e1,e2,e00,v,dv,val,w,valm,dnd,xx,yy,pnue,pnue0,
     .          hw,shift,dvdw,de,hnorm,g1,g2,d1,d2,cvel

      if (kkin /= idreac) call get_reaction(kkin)

      kk=kkin

      iid = reaction%ircart
      fac=0._dp
      hwvdw=0._dp

      select case(ity)
      case(0)
         ipl = reaction%ignd
      case(1)
         do nrc=1,nrca(isp)
            if(kk == ireaca(isp,nrc)) exit
         enddo
         ipl = idez(ibulka(iatm,nrc),3,3)
      case default
         ipl=0
      end select

      select case(iid)
      case(1)
c     P.2 TEST OT:
         write(6,*)'PHOTON MODULE (PH_GETCOEFF): IID=1 not in use'
         res=0.         

      case(2,3)
c     P.2 AT_ABS OT, P.2 AT_STIM OT
         if(.not.allocated(splinemax)) return

         if(splinemax(icell) <= 0.) then
            res=0.
            return
         endif

         spx => SPLINEX(icell,:)
         spa => SPLINEA(icell,:)
         spb => SPLINEB(icell,:)
         spc => SPLINEC(icell,:)

         e00 = reaction%e0
         e1 = reaction%e1
         e2 = e1 + e00

         select case(reaction%iprofiletype)
         case(0,1)
         case(2,3)
            call ph_homprof(icell, gam)            
         case(4,5)            
            call ph_vdwprof (icell,hw,shift,dvdw,.true.)
         end select

         n = volspecnum(ivs)
         n4= n*4
         v = eemin(icell)
         dv= (eemax(icell)-eemin(icell))/n4
         res=0.
         iwarn=0
         valm=0.
         do i=1,n4

            select case(reaction%iprofiletype)
            case(0)
               w=v
               fac=1.
            case(1)
               w=v*(1. - vel/clight)
               fac=1.
            case(2)
               w=v
               fac=ph_lorentz(w-e00, gam)
            case(3)
               w=v*(1. - vel/clight)
               fac=ph_lorentz(w-e00,gam)
            case(4)
               w=v
               iil = phphilphil2line(pil)
               if(iil > 0) then
c               if(iil /= il) then
c                  write(6,*) 'PHOTON MODULE (PH_GETCOEFF):'
c                  write(6,*) '   iil/=il'
c                  write(6,*) '   (2) il,iil,pil,iid=',il,iil,pil,iid
c                  call exit(1)
c               endif
                  de = w-e00
                  fac = ph_lorvdw(de, hw, shift, dvdw,icell,0)
               else
                  fac=1.
               endif
            case(5)
!pb   cvel?
               w=v - e00*cvel/clight
               iil = phphilphil2line(pil)
               if(iil > 0) then
c               if(iil /= il) then
c                  write(6,*) 'PHOTON MODULE (PH_GETCOEFF):'
c                  write(6,*) '   iil/=il'
c                  write(6,*) '   (2) il,iil,pil,iid=',il,iil,pil,iid
c                  call exit(1)
c               endif
                  de = w-e00
                  fac = ph_lorvdw(de, hw, shift, dvdw,icell,0)
               else
                  fac=1.
               endif
            end select

            val = ph_xquaval(n,v,iflag,icell)
            if(iflag /= 0) then
               select case(iflag)
               case(3)
                  iwarn=iwarn+1
                  val=0.
               case default
                  write(6,*)'PHOTON MODULE (PH_GETCOEFF), iflag=',iflag
                  call exit(1)
               end select
            endif

c           avoid neg. spline values (implies an error!)
            if(val>0.) then
               res=res+val*dv*fac
               valm=max(valm,val)
            endif
            v=v+dv
         enddo
         select case(iid)
         case(2)
            res=res*reaction%B12
         case(3)
            res=res*reaction%B21
         end select
         res=res/(4._dp*PIA)
c result should have units 1/s --> no multiplication of photon density needed,
c already included in background splines
         phv_muldens=0

      case(4,5)
c     P.2 PH_ABS OT, P.2 PH_STIM OT
         if(lgvac(icell,ipl)) then
            res=0.
            return
         endif

         e00=reaction%e0
!        write (40,*) ' e1, e2, e00 ', e1, e2, e00
         pnue0 = e00*EV2HZ
         pnue  = e0 *EV2HZ

         
         iptype = reaction%iprofiletype
         select case(iptype)
         case(0)
            fac=1.
         case(1)
            call ph_dopplerprof(iipl,icell,dnd)
            xx=e0-e00
            fac = PH_GAUSS(xx,dnd)
         case(2)
            call ph_homprof(icell,gam)
            fac = PH_LORENTZ(e0-e00, gam)
         case(3)
            call ph_voigtprof(iipl,icell,dnd,gam)
            xx=(e0-e00)/dnd
            yy=gam/dnd
            fac = DBLE(PH_FADDEEVA(xx,yy,dnd,icell,iipl))
!            fac = real(PH_FADDEEVA(xx,yy,dnd,icell,iipl),dp)
         case(4)
            call  ph_vdwprof(icell,hw,shift,dvdw,.true.)
            hwvdw = -hw
            de = e0-e00
            fac = ph_lorvdw(de, hw, shift, dvdw, icell,0)
         case(5)
            write (6,*) ' ph_getcoeff for photons,'
            write (6,*) ' iptype == 5 not implemented yet'
            call exit(1)
         end select   

         select case(iid)
         case(4)
            res=e00*reaction%b12
         case(5)
            res=e00*reaction%b21
         end select
         res=res*fac*clight/(4._dp*PIA)
c result should have units cm**3/s --> multiply backgr. density!
         phv_muldens=1

      case(6)
c     P.1 AT_AIK OT
         if(ipl > 0) then
            write(6,*) 'PHOTON MODULE (PH_GETCOEFF):'
            write(6,*) '       IID=6, ipl > 0'
         endif
         res=reaction%aik
c result should have units 1/s, no multiply of any density
         phv_muldens=0

      case(7)
c     P.2 PH_CABS OT
         if(lgvac(icell,ipl)) then
            res=0.
            return
         endif

         g1 = reaction%g1
         g2 = reaction%g2
         e00 = reaction%e0
         e1 = reaction%e1
         e2 = e00 + e1
         pnue0 = e00*EV2HZ
         pnue  = e0 *EV2HZ

         pil=phphilbulk2phil_l(iipl)

         select case(reaction%iprofiletype)
         case(0)
            fac=1.
         case(1)
            call ph_dopplerprof(iipl,icell,dnd)
            xx=e0-e00
            fac = PH_GAUSS(xx,dnd)
         case(2)
            call ph_homprof(icell,gam)
            fac = PH_LORENTZ(e0-e00, gam)
         case(3)
            call ph_voigtprof(iipl,icell,dnd,gam)
            xx=(e0-e00)/dnd
            yy=gam/dnd
            fac = DBLE(PH_FADDEEVA(xx,yy,dnd,icell,iipl))
         case(4)
            ignd=reaction%ignd
            call  ph_vdwprof(icell, hw,shift,dvdw,.true.)
            de = e0-e00
            fac = ph_lorvdw(de, hw, shift, dvdw, icell,0)
         end select   

         d1=DIIN(iipl,icell)
         if ((phv_n1stotph(isp,idsc,1) == 4) .and. 
     .       (phv_n1stotph(isp,idsc,3) == 1)) then    
            ipl2=phv_n1stotph(isp,idsc,2)
            d2=DIIN(ipl2,icell)
         else
            d2=0._dp
            write (6,*) ' WARNING !! '
            write (6,*) ' WRONG TYPE OF SECONDARY PARTICLE IN',
     .                  ' PH_GETCOEFF'
            write (6,*) ' RATE CHANGED ACCORDINGLY '
         END IF
         
         res=e00*reaction%B12*(1._dp- g1*d2/g2/d1)
         res=res*fac*clight/(4._dp*PIA)

c yes, multiply density
         phv_muldens=1
c         res=0.
      case default
         write(6,*)'PHOTON MODULE (PH_GETCOEFF): iid=',iid,'not in use'
         res=0.
      end select
      
      return
      END SUBROUTINE PH_GETCOEFF

      REAL(DP) FUNCTION PH_B12() result(res)
      IMPLICIT NONE
c calculates B12 Einstein coefficient in units: cm**2
      integer :: g1,g2,n1,n2

      res=PH_B21()

      g1=reaction%g1
      g2=reaction%g2

      res=res*g2/g1
      return      
      END FUNCTION PH_B12

      REAL(DP) FUNCTION PH_B21() result(res)
      IMPLICIT NONE
c calculates B21 Einstein coefficient in units: cm**2
      
      real(DP) :: e00,pnue0

      e00=reaction%e0         
      
      if(e00 > 0.0) THEN
         pnue0 = E00*EV2HZ
         res=reaction%aik
c         res=res * CLIGHT**2/(2.0*pnue0**3)
         res=res * (hplnk*clight)**2 / (2.*e00**3)
c convert [cm^2 / (eV*s) ] --> [cm^2]
         res=res*hplnk
      else
         write(6,*) 'PHOTON MODULE (PH_B21): e00 <= 0.'
         res=0.
      endif
      return
      END FUNCTION PH_B21

c SAMPLING
      REAL(dp) FUNCTION PH_ENERGY(il,icell,kk,ipl2,nldoppl) result(res)
      IMPLICIT NONE
cdr  sample frequency (here: energy) from emission profile
cdr  iprofiletype =0  only line centre (delta),  stationary atoms
cdr  iprofiletype =1  ditto plus doppler in calling program
cdr  iprofiletype =2  only lorentz (delta),  stationary atoms
cdr  iprofiletype =3  ditto plus doppler in calling program
cdr  iprofiletype =4  lorentz+vdw ,  stationary atoms
cdr  iprofiletype =5  ditto plus doppler in calling program (new, july 04)
cdr
cdr
cdr
cdr
cdr  july 04: nldoppl introduced
cdr             now ph_energy only samples in the rest frame of the
cdr             emitting atom.
cdr             it returns nldoppl=true, if doppler shift is to be
cdr                                      added in calling program
cdr             
      integer, intent(in) :: il,icell,kk,ipl2
      logical, intent(out) :: nldoppl 
      real(dp) :: e00,pnue0,dnd,gam,hw,shift,dvdw,zep1,l0
      integer :: pil,iil,dum

!  idreac: kk from last call to get_reaction
      if (idreac /= kk) call get_reaction(kk)
      e00=reaction%e0
c     pnue0=e00*ev2hz

      select case(reaction%iprofiletype)
      case(0)
c  delta, no doppler
         res=e00
         nldoppl=.false.
      case(1)
c  delta plus doppler
cdr      call ph_dopplerprof(ipl2,icell,dnd)
cdr      zep1 = ph_sam_gauss(dnd/sqrt(2.))
cdr      res=e00 + zep1
         res=e00
         nldoppl=.true.
      case(2)
c  lorentz, no doppler
         call ph_homprof(icell, gam)
         res = ph_sam_lorentz(gam,e00)
         nldoppl=.false.
      case(3)
c  lorentz plus doppler
cdr      call ph_voigtprof(ipl2,icell,dnd,gam)
cdr      zep1 = ph_sam_voigt(gam,dnd/sqrt(2.))
cdr      res=e00 + zep1
         call ph_homprof(icell, gam)
         res = ph_sam_lorentz(gam,e00)
         nldoppl=.true.
      case(4)
c  lorentz plus vanderWaals, no doppler
         call ph_vdwprof(icell,hw,shift,dvdw,.false.)
         l0 = hpcl/e00
         zep1 = ph_sam_lorvdw2(hw,shift,dvdw,l0)
         res = hpcl/zep1
         nldoppl=.false.
      case(5)
c  lorentz plus vanderWaals plus doppler
         call ph_vdwprof(icell,hw,shift,dvdw,.false.)
         l0 = hpcl/e00
         zep1 = ph_sam_lorvdw2(hw,shift,dvdw,l0)
         res = hpcl/zep1
         nldoppl=.true.

cdr   case(6) :zeeman, no doppler
cdr   case(7) :zeeman plus doppler

      end select
      weight = 1.

      return
      END FUNCTION PH_ENERGY
c
c
      REAL(dp) FUNCTION PH_SAMVOL_ANGLE(mode) result(res)
      IMPLICIT NONE

      integer, intent(in) :: mode
      real(dp) :: vx,vy,vz,velq
      integer :: dum

      select case(mode)
      case(1)
c        mode  1: isotropic
         IF (INIV3.EQ.0) CALL FISOTR
         vx=FI1(INIV3)*clight
         vy=FI2(INIV3)*clight
         vz=FI3(INIV3)*clight
         INIV3=INIV3-1

         VEL=CLIGHT
         VELQ=VEL*VEL
         velx=vx/vel
         vely=vy/vel
         velz=vz/vel
         res=vel        
      case(2)
c        mode  2: no spec. contr. (f2=0,f1=1)
         IF (INIV4.LE.0) CALL FCOSIN
         VX=FC1(INIV4)
         VY=FC2(INIV4)
         VZ=FC3(INIV4)
         INIV4=INIV4-1
         CALL ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
         VEL=CLIGHT
         res=vel
      case(-1)
c        mode -1: debugging purpose
         vel=clight
         velq = vel*vel
         velx=0.
         vely=0.
         velz=1.
         res=vel
      case(-2)
c        directions from vxc...arrays
         vel=clight
         velq=vel*vel
         
         velx=vxc(phv_dnum)
         vely=vyc(phv_dnum)
         velz=vzc(phv_dnum)
         res=vel
      case default
         res=0.
      end select
      return
      END FUNCTION PH_SAMVOL_ANGLE

c POST-COLLISION
      SUBROUTINE PH_POST_ENERGY(icell,kk,iflg,il,
     .                          iold,itypold,vxo,vyo,vzo,vlo,e0o,
     .                          itypnew)
c sample post collision energy.
c
c incident particle: (iold,itypold,vxo,....e0o)
c already decided: new test particle has type itypnew
c
c called from colatm (itypold=1), or called from colphot (itypold=0)
c
c enerphot
c iflg = 0: sp.emission, no bulk precoll partner
c        1: absorption
c        2: stim.em
      IMPLICIT NONE
      integer, intent(in) :: icell,kk,iflg,iold,il,itypold,itypnew
      real(dp),intent(in) :: vxo,vyo,vzo,vlo,e0o
      integer :: n1,n2,nrc,ipln,itypn,imax,ii,iwarn,n,ir,ivs,
     .    ityp0,ityp1,ityp2,ipl0,ipl1,ipl2,iflag,pil,iipl,iil,ignd,
     .    ipl0v
      real(dp) :: vx,vy,vz,vxn,vyn,vzn,cvrss1,velq,e1,e2,e00,gam,
     .    smax,swma,swmi,zep1,zep2,zzep1,cangl,val,een,
     .    hw,shift,dvdw,ee,cvel,zs,zc,vxoo,vyoo,vzoo,vloo,l0,e00s,
     .    velx_b, vely_b, velz_b, velparm, vel_b

      if (idreac /= kk) call get_reaction(kk)

      ir=kk
!pb      ivs=line2volspec(il)

      select case(itypold)
      case(0)
         do nrc=1,nrcph(iold)
            if(kk == ireacph(iold,nrc)) exit            
         enddo
         IPL0 =IDEZ(IBULKPH(IOLD,nrc),3,3)
         IPL1 =IDEZ(ISCD1PH(IOLD,nrc),3,3)
         IPL2 =IDEZ(ISCD2PH(IOLD,nrc),3,3)
         ITYP0=IDEZ(IBULKPH(IOLD,nrc),1,3)
         ITYP1=IDEZ(ISCD1PH(IOLD,nrc),1,3)
         ITYP2=IDEZ(ISCD2PH(IOLD,nrc),1,3)
      case(1)
         do nrc=1,nrca(iold)
            if(kk == ireaca(iold,nrc)) exit
         enddo
         IPL0 =IDEZ(IBULKA(IOLD,nrc),3,3)
         IPL1 =IDEZ(ISCD1A(IOLD,nrc),3,3)
         IPL2 =IDEZ(ISCD2A(IOLD,nrc),3,3)
         ITYP0=IDEZ(IBULKA(IOLD,nrc),1,3)
         ITYP1=IDEZ(ISCD1A(IOLD,nrc),1,3)
         ITYP2=IDEZ(ISCD2A(IOLD,nrc),1,3)
      case default
         nrc=-1
      end select

      ipl0v = mplsv(ipl0)
c
c
      select case(iflg)
      case(0)
c     spont. emission, no bulk particle as coll.partner
c
c  this part: itypold=1
c
         if (itypold==1) then
            if(ipl0 /= 0) then
               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
               write(6,*) ' sp.emission process, ipl0 /= 0 ?!?'
               call exit(1)
            endif
            if(ipl1 == 0 .or. ipl2 == 0) then
               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
               write(6,*) '       at leasttwo secondaries in sp.emis.'
               write(6,*) '       check input file'
               write(6,*) '       n=2 -> n=1 + ph'
               call exit(1)
            endif
         elseif (itypold==0) then
c     to be written
         endif

         if(ityp1 == itypnew) then
            itypn = ityp1
            ipln  = ipl1
         elseif(ityp2 == itypnew) then
            itypn = ityp2
            ipln  = ipl2
         else
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   ityp of secondary wrong'
            write(6,*) '   (sp.emission)'
            call exit(1)
         endif
         
         secs: select case(itypn)
         case(0)
c        photon followed:
c
c        new direction, isotropically
            vel=PH_POST_ANGLE(1)

            IF (ITYPOLD == 0) THEN
c  sample velocity from ipl0 bulk, shifted maxw
c  this case has been added (cdr, aug.23, 2004)
c  "complete re-distribution", 
c  for doppler here first sample a background velocity of bulk atom 
c  (similar to subr. locate)
               IF(INIV2.EQ.0) CALL FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1
               
               VELPARM=ZRG(ipl0,icell)
               
               vx=VELPARM*VXN+VXIN(ipl0v,icell)
               vy=VELPARM*VYN+VYIN(ipl0v,icell)
               vz=VELPARM*VZN+VZIN(ipl0v,icell)
               velq=vx*vx+vy*vy+vz*vz
               vel_b=dsqrt(velq)
               
               VELX_B=VX/vel_b
               VELY_B=VY/vel_b
               VELZ_B=VZ/vel_b
               cangl = (velx*velx_b + vely*vely_b + velz*velz_b)
               cvel  = vel_b*cangl
            ELSE
c  this case for incident atoms, ie, old velocity for doppler already known
c  already foreseen in Nov03 version.
               cangl = (velx*vxo + vely*vyo + velz*vzo)
               cvel  = vlo*cangl
            END IF

c        sample new freq, doppler shift due n=2 velocity
            e00=reaction%e0
            e1=reaction%e1
            e2=e00+e1
            
            labelp1: select case(reaction%iprofiletype)
            case(0)
               ee = e00
            case(1)
               ee = e00*(1.-cvel/clight)
            case(2)
               call PH_HOMPROF(icell,gam)
               zep1=ph_sam_lorentz(gam,e00)
c               ee = e00 + zep1
               ee = zep1
            case(3)
               call PH_HOMPROF(icell,gam)
               een = e00*(1.-cvel/clight)
               zep1=ph_sam_lorentz(gam,een)
c               ee = een + zep1
               ee = zep1
            case(4)
               ee = e00
               ignd=reaction%ignd

c constants in wavelengths 
               call ph_vdwprof(icell, hw,shift,dvdw,.true.)
               l0 = hpcl/ee
               zep1=ph_sam_lorvdw2(hw,shift,dvdw,l0)
c                ee = ee + zep1
c                ee = hpcl/(l0+zep1)
               ee = hpcl/zep1
            case(5)
               ee = e00
               ignd=reaction%ignd

c constants in wavelengths 
               call ph_vdwprof(icell, hw,shift,dvdw,.true.)
               l0 = hpcl/ee
               zep1=ph_sam_lorvdw2(hw,shift,dvdw,l0)
c                ee = ee + zep1
c                ee = hpcl/(l0+zep1)
               ee = hpcl/zep1
               ee = ee*(1.-cvel/clight)
            case(6)
               call ph_strkprof(icell,hw,shift)
               call PH_HOMPROF(icell,gam)
               gam=gam+hw

               e00s = e00+shift
               een = e00s
               zep1=ph_sam_lorentz(gam,een)
c               ee = ee + zep1
               ee = zep1
            case(7)
               call ph_strkprof(icell,hw,shift)
               call PH_HOMPROF(icell,gam)
               gam=gam+hw

               e00s = e00+shift
               een = e00s*(1.-cvel/clight)
               zep1=ph_sam_lorentz(gam,een)
c               ee = ee + zep1
               ee = zep1
            end select labelp1            
            e0=ee
            return
            
         case(1)
c     n=1 followed
c     use old values from n=2 testpart.
            velx=vxo
            vely=vyo
            velz=vzo
            vel=vlo
            e0=e0o
            return
         case(4)
c     bulk particle (photon or atom)
            vel=0.
            e0=0.
            return
         case default
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   default case?!, check input!'
            call exit(1)
         end select secs
c
c  spontane emission , for incident atom, finished 
c
      case(1)
c     aborption process
         if(ipl1 /= 0 .and. ipl2 /= 0) then
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '       only 2 secondaries allowed in abs.'
            write(6,*) '       process, check input file'
            write(6,*) '       ph + n=1 -> n=2'
            call exit(1)
         endif
         if(ityp1 == itypnew) then
            itypn = ityp1
            ipln  = ipl1
         elseif(ityp2 == itypnew) then
            itypn = ityp2
            ipln  = ipl2
         else
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   ityp of secondary wrong'
            write(6,*) '   (absorption)'
            call exit(1)
         endif

         abcs: select case(itypold)
         case(0)
c        photon on n=1 background, n=2 followed
            abcs2: select case(itypn)
            case(1)
c           followed n=2 is atom
c           sample new energy from ipl0 bulk, shifted maxw
               IF(INIV2.EQ.0) CALL FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1
               
               cvrss1=cvrssa(ipln)
               VEL=ZRG(ipl0,icell)

c               vx=VEL*VXN+VXIN(ipl0v,icell)
c               vy=VEL*VYN+VYIN(ipl0v,icell)
c               vz=VEL*VZN+VZIN(ipl0v,icell)
               vx=VEL*VXN
               vy=VEL*VYN
               vz=VEL*VZN
               velq=vx*vx+vy*vy+vz*vz
               vel=dsqrt(velq)
               
               VELX=VX/vel
               VELY=VY/vel
               VELZ=VZ/vel
               E0=CVRSS1*VELQ
               return
            case(4)
c               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c               write(6,*) '       secondary in abs.process is bulk'
c               write(6,*) '       somethings wrong in input file'
c               write(6,*) '       ph + n=1(bulk) --> n=2(bulk)'
c               call exit(1)
               e0=0.
               vel=0.
               return
            end select abcs2

         case(1)
c        n=1 on photon background, n=2 followed
            abcs3: select case(itypn)
            case(1)
c           followed n=2 is atom
c           use energy from n=1 testparticle
               cvrss1=cvrssa(ipln)
               velx=vxo
               vely=vyo
               velz=vzo
               vel =vlo
               velq=vel*vel
               e0=e0o
               return
            case(4)
c           followed n=2 is bulk
c               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c               write(6,*) '       secondary in abs.process is bulk'
c               write(6,*) '       somethings wrong in input file'
c               write(6,*) '       n=1 + ph(bulk) --> n=2(bulk)'
c               call exit(1)
               e0=0.
               vel=0.
               return
            end select abcs3
         case default
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   default case?!, check input!'
            call exit(1)
         end select abcs

      case(2)
c     stimulated process
         if(ipl1 == 0 .or. ipl2 == 0) then
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '       two secondaries are needed in stim.em.'
            write(6,*) '       check input file'
            write(6,*) '       (ph + n=2) -> n=1 + 2ph'
            call exit(1)
         endif
         if(ityp1 == itypnew) then
            itypn = ityp1
            ipln  = ipl1
         elseif(ityp2 == itypnew) then
            itypn = ityp2
            ipln  = ipl2
         else
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   ityp of secondary wrong'
            write(6,*) '   (stim.emission)'
            call exit(1)
         endif

         stcs: select case(itypold)
         case(0) 
c        photon on n=2 background
            stcs2: select case(itypn)

            case(0)
c           two photons followed
c           use old energy from photon(ipl0)
            e0=e0o
            velx=vxo
            vely=vyo
            velz=vzo
            vel =vlo
            velq=vlo*vlo
            return

            case(1)
c           n=1 followed
c           sample new energy from ipl0 bulk, shifted maxwellian
               IF(INIV2.EQ.0) CALL FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1
               
               cvrss1=cvrssa(ipln)
               VEL=ZRG(ipl0,icell)

c               vx=VEL*VXN+VXIN(ipl0v,icell)
c               vy=VEL*VYN+VYIN(ipl0v,icell)
c               vz=VEL*VZN+VZIN(ipl0v,icell)
               vx=VEL*VXN
               vy=VEL*VYN
               vz=VEL*VZN
               velq=vx*vx+vy*vy+vz*vz
               vel=dsqrt(velq)
               
               VELX=VX/vel
               VELY=VY/vel
               VELZ=VZ/vel
               E0=CVRSS1*VELQ
               return
c     bulk particle
            case(4)
               e0=0.
               vel=0.
               return
            case default
               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
               write(6,*) '   default case?!, check input!'
               call exit(1)
            end select stcs2

         case(1)
c        n=2 on photon background
            stcs3: select case(itypn)

            case(0)
c           two photons followed               
c           sample from photon background, doppler shift from iold
               e00=reaction%e0
               e1=reaction%e1
               e2=e00+e1

c           new direction, isotropically
               vloo=vel
               vxoo=velx
               vyoo=vely
               vzoo=velz
               vel=PH_POST_ANGLE(1)
               cangl = (velx*vxoo + vely*vyoo + velz*vzoo)
               cvel  = vloo*cangl

               labelp3: select case(reaction%iprofiletype)
               case(1)
                  ee = e00*(1.-cvel/clight)
               case(2)
                  call PH_HOMPROF(icell,gam)
                  ee = e00
               case(3)
                  call PH_HOMPROF(icell,gam)
                  ee = e00*(1.-cvel/clight)
               case(4)
                  ignd=reaction%ignd
                  ee = e00

c constants in wavelengths
                  call ph_vdwprof(icell,hw,shift,dvdw,.true.)
               case(5)
                  call ph_strkprof(icell,hw,shift)
                  call PH_HOMPROF(icell,gam)
                  gam=gam+hw 

                  e00s = e00+shift
                  ee = e00s*(1.-cvel/clight)
               case(0)
                  ee = e00
               end select labelp3

               smax=splinemax(icell)
c debug:
c               smax = ph_planck(tiin(10,icell),e00)
               if(smax <= 0.) then
                  write(6,*) 'PHOTON MODULE (PH_POST_ENERGY): WARNING'
                  write(6,*) '   splinemax for icell',icell,'<= 0.'
                  write(6,*) '   use of linecenter, isotropic'
                  e0=e00
                  return
               else
                  SPX => SPLINEX(icell,:)
                  SPA => SPLINEA(icell,:)
                  SPB => SPLINEB(icell,:)
                  SPC => SPLINEC(icell,:)
                  swma = eemax(icell)
                  swmi = eemin(icell)
                  n=volspecnum(ivs)
                  imax=10000
                  ii=0
                  iwarn=0
                  l0=0.
                  DO WHILE(ii<imax)

                     labelp4: select case(reaction%iprofiletype)
                     case(1)
                        zep1=0.
                     case(2)
                        zep1=ph_sam_lorentz(gam,ee)
                     case(3)
                        zep1=ph_sam_lorentz(gam,ee)
                     case(4)
                        iil=phphilphil2line(pil)
                        if(iil > 0) then
c                        if(iil /= il) then
c                           write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c                           write(6,*) '   iil/=il'
c                           write(6,*) '   (4) il,iil,pil=',il,iil,pil
c                           call exit(1)
c                        endif
                           l0 = hpcl/ee
                           zep1=ph_sam_lorvdw2(hw,shift,dvdw,l0)
                           
                        else
                           zep1=0.
                        endif
                     case(5)
                        zep1=ph_sam_lorentz(gam,ee)
                     case(0)
                        zep1=0.
                     end select labelp4

                     if(abs(l0) > 0.) then
c                        zzep1= hpcl/(l0+zep1)
                        zzep1 = hpcl/zep1
                     else
c                        ZZEP1=ee + zep1
                        zzep1 = zep1
                     endif

                     val = PH_XQUAVAL(n,zzep1,iflag,icell)
c debug:
c                     iflag=0
c                     val = ph_planck(tiin(10,icell),zzep1)

                     if(iflag.ne.0) THEN
                        ilferr: select case(iflag)
                        case(3)
                           iwarn=iwarn+1
                        case default    
                           write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
                           write(6,*) '   photon sampling, iflag=',iflag
                           call exit(1)
                        end select ilferr
                     endif                

                     if(val < 0.) val = 0.

                     ZEP2=ranf_eirene()*smax                 
                     if(zep2 <= val) then
                        e0=zzep1
                        exit
                     endif
                     ii=ii+1

                  ENDDO
                  IF(ii>=imax) THEN
                     write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
                     write(6,*) '   MAXIMUM NUM OF REJECTION-SAMPLINGS'
                     write(6,*) '   REACHED, imax,ii=',imax,ii
                     write(6,*) '   using e00 as energy'
                     e0 = e00
                  ENDIF
               endif               
               return

            case(1)
c           n=1 followed           
c           use energy from n=2 testparticle
               cvrss1=cvrssa(ipln)
               velx=vxo
               vely=vyo
               velz=vzo
               vel =vlo
               velq=vel*vel
               e0=e0o
               return
            case(4)
c bulk particle (photon or atom)
               vel=0.
               e0=0.
               return
            case default
               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
               write(6,*) '   default case?!, check input!'
               call exit(1)
            end select stcs3
         case default
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   default case?!, check input!'
            call exit(1)
         end select stcs
      end select
      
      write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
      write(6,*) '   end of subroutine reached'
      write(6,*) '   there has to be an error, check input!'
      call exit(1)
      END SUBROUTINE PH_POST_ENERGY


c POST-COLLISION
      SUBROUTINE PH_POST_ENERGY_OLD(icell,kk,iflg,il,
     .                          iold,itypold,vxo,vyo,vzo,vlo,e0o)
c enerphot
c iflg = 0: sp.emission, no bulk precoll partner
c        1: absorption
c        2: stim.em
      IMPLICIT NONE
      integer, intent(in) :: icell,kk,iflg,iold,il,itypold
      real(dp),intent(in) :: vxo,vyo,vzo,vlo,e0o
      integer :: n1,n2,nrc,ipln,itypn,imax,ii,iwarn,n,ir,ivs,
     .    ityp0,ityp1,ityp2,ipl0,ipl1,ipl2,iflag,pil,iipl,iil,ignd
      real(dp) :: vx,vy,vz,vxn,vyn,vzn,cvrss1,velq,e1,e2,e00,gam,
     .    pnue0,pnue,smax,swma,swmi,zep1,zep2,zzep1,cangl,val,
     .    hw,shift,dvdw,ee,pnuehz

      if (idreac /= kk) call get_reaction(kk)

      ir=kk
!pb      ivs=line2volspec(il)

      select case(itypold)
      case(0)
         do nrc=1,nrcph(iold)
            if(kk == ireacph(iold,nrc)) exit            
         enddo
         IPL0 =IDEZ(IBULKPH(IOLD,nrc),3,3)
         IPL1 =IDEZ(ISCD1PH(IOLD,nrc),3,3)
         IPL2 =IDEZ(ISCD2PH(IOLD,nrc),3,3)
         ITYP0=IDEZ(IBULKPH(IOLD,nrc),1,3)
         ITYP1=IDEZ(ISCD1PH(IOLD,nrc),1,3)
         ITYP2=IDEZ(ISCD2PH(IOLD,nrc),1,3)
      case(1)
         do nrc=1,nrca(iold)
            if(kk == ireaca(iold,nrc)) exit
         enddo
         IPL0 =IDEZ(IBULKA(IOLD,nrc),3,3)
         IPL1 =IDEZ(ISCD1A(IOLD,nrc),3,3)
         IPL2 =IDEZ(ISCD2A(IOLD,nrc),3,3)
         ITYP0=IDEZ(IBULKA(IOLD,nrc),1,3)
         ITYP1=IDEZ(ISCD1A(IOLD,nrc),1,3)
         ITYP2=IDEZ(ISCD2A(IOLD,nrc),1,3)
      case default
         nrc=-1
      end select
c
c
      select case(iflg)
      case(0)
c     spont. emission, no bulk particle as coll.partner
         if(ipl0 /= 0) then
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) ' sp.emission process, ipl0 /= 0 ?!?'
            call exit(1)
         endif
         if(ityp1 == 4 .and. ityp2 == 4) then
c            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c            write(6,*) '       both secondaries are bulk in sp.emis.'
c            write(6,*) '       check input file'
c            write(6,*) '       n=2 -> n=1(bulk)+ ph(bulk)'
c            call exit(1)
            e0=0.
            vel=0.
            return
         endif
         if(ityp1 /= 4) then
            itypn = ityp1
            ipln  = ipl1
         else
            itypn = ityp2
            ipln  = ipl2
         endif
         secs: select case(itypn)
         case(0)
c        photon followed:
c
c        new direction, isotropically
            vel=PH_POST_ANGLE(1)
            cangl = (velx*vxo + vely*vyo + velz*vzo)/(vlo*vel)

c        sample new freq, doppler shift due n=2 velocity
            e00=reaction%e0
            e1=reaction%e1
            e2=e00+e1
            pnue0 = e00*ev2hz               
            
            labelp1: select case(reaction%iprofiletype)
            case(1)
               pnue  = pnue0*(1.-vlo*cangl/clight)
            case(2)
               call PH_HOMPROF(icell,gam)
               pnue  = pnue0
            case(3)
               call PH_HOMPROF(icell,gam)
               pnue  = pnue0*(1.-vlo*cangl/clight)
            case(4)
               ignd=reaction%ignd
               call ph_vdwprof(icell,hw,shift,dvdw,.true.)
               pnue  = pnue0
            case(0)
               pnue  = pnue0
            end select labelp1

            swmi=eemin(icell)
            swma=eemax(icell)

            imax=1000
            ii=0
	    pnuehz=pnue/ev2hz
            DO WHILE(ii<imax)

               labelp2: select case(reaction%iprofiletype)
               case(1)
                  ee=pnuehz
               case(2)
                  ee=ph_sam_lorentz(gam,pnuehz)
               case(3)
                  ee=ph_sam_lorentz(gam,pnuehz)
               case(4)
!pb                  iil=phphilphil2line(pil)
!pb                  if(iil > 0) then
                     ee=ph_sam_lorvdw2(hw,shift,dvdw,pnuehz)
!pb                  else
!pb                     ee=pnuehz
!pb                  endif
               case(0)
                  ee=pnuehz
               end select labelp2

               E0=ee
               ii=ii+1
            enddo
            if(ii.ge.imax) then
               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
               write(6,*) '    too much samplings in spnt. decay'
               call exit(1)
            endif
            return
            
         case(1)
c        n=1 followed
c        use old values from n=2 testpart.
            velx=vxo
            vely=vyo
            velz=vzo
            vel=vlo
            e0=e0o
            return
         case default
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   default case?!, check input!'
            call exit(1)
         end select secs

      case(1)
c     aborption process
         abcs: select case(itypold)
         case(0)
c        photon on n=1 background, n=2 followed
            abcs2: select case(ityp1)
            case(1)
c           followed n=2 is atom
c           sample new energy from ipl0 bulk
               cvrss1=cvrssa(ipl1)
               IF(INIV2.EQ.0) CALL FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1
               
               VEL=SQRT(ZT1(ipl0,icell))
c               vx=VEL*VXN+VXIN(ipl0v,icell)
c               vy=VEL*VYN+VYIN(ipl0v,icell)
c               vz=VEL*VZN+VZIN(ipl0v,icell)
               vx=VEL*VXN
               vy=VEL*VYN
               vz=VEL*VZN
               velq=vx*vx+vy*vy+vz*vz
               vel=sqrt(velq)
               
               VELX=VX/vel
               VELY=VY/vel
               VELZ=VZ/vel
               E0=CVRSS1*VELQ
               return
            case(4)
c               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c               write(6,*) '       secondary in abs.process is bulk'
c               write(6,*) '       somethings wrong in input file'
c               write(6,*) '       ph + n=1(bulk) --> n=2(bulk)'
c               call exit(1)
               e0=0.
               vel=0.
               return
            end select abcs2

         case(1)
c        n=1 on photon background, n=2 followed
            abcs3: select case(ityp1)
            case(1)
c           followed n=2 is atom
c           use energy from n=1 testparticle
               cvrss1=cvrssa(ipl1)
               velx=vxo
               vely=vyo
               velz=vzo
               vel =vlo
               velq=vel*vel
               e0=cvrss1*velq
               return
            case(4)
c           followed n=2 is bulk
c               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c               write(6,*) '       secondary in abs.process is bulk'
c               write(6,*) '       somethings wrong in input file'
c               write(6,*) '       n=1 + ph(bulk) --> n=2(bulk)'
c               call exit(1)
               e0=0.
               vel=0.
               return
            end select abcs3
         case default
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   default case?!, check input!'
            call exit(1)
         end select abcs

      case(2)
c     stimulated process
         if(ityp1 == 4 .and. ityp2 == 4) then
c            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c            write(6,*) '       both secondaries are bulk in stim.em.'
c            write(6,*) '       check input file'
c            write(6,*) '       (ph + n=2) -> n=1(bulk)+ 2ph(bulk)'
c            call exit(1)
            e0=0.
            vel=0.
            return
         endif
         if(ityp1 /= 4) then
            itypn = ityp1
            ipln  = ipl1
         else
            itypn = ityp2
            ipln  = ipl2
         endif

         stcs: select case(itypold)
         case(0) 
c        photon on n=2 background
            stcs2: select case(itypn)

            case(0)
c           two photons followed
c           use old energy from photon(ipl0)
            e0=e0o
            velx=vxo
            vely=vyo
            velz=vzo
            vel =vlo
            velq=vlo*vlo
            return

            case(1)
c           n=1 followed
c           sample new energy from ipl0 bulk, shifted maxwellian
               cvrss1=cvrssa(ipl0)
               IF(INIV2.EQ.0) CALL FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1
               
               VEL=SQRT(ZT1(ipl0,icell))
c               vx=VEL*VXN+VXIN(ipl0v,icell)
c               vy=VEL*VYN+VYIN(ipl0v,icell)
c               vz=VEL*VZN+VZIN(ipl0v,icell)
               vx=VEL*VXN
               vy=VEL*VYN
               vz=VEL*VZN
               velq=vx*vx+vy*vy+vz*vz
               vel=sqrt(velq)
               
               VELX=VX/vel
               VELY=VY/vel
               VELZ=VZ/vel
               E0=CVRSS1*VELQ
               return

            case default
               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
               write(6,*) '   default case?!, check input!'
               call exit(1)
            end select stcs2

         case(1)
c        n=2 on photon background
            stcs3: select case(itypn)

            case(0)
c           two photons followed               
c           sample from photon background, doppler shift from iold
               e00=reaction%e0
               e1=reaction%e1
               e2=e00+e1
               pnue0 = e00*ev2hz

               labelp3: select case(reaction%iprofiletype)
               case(1)
                  pnue  = pnue0*(1.-vel/clight)
               case(2)
                  call PH_HOMPROF(icell,gam)
                  pnue  = pnue0
               case(3)
                  call PH_HOMPROF(icell,gam)
                  pnue  = pnue0*(1.-vel/clight)
               case(4)
                  ignd=reaction%ignd
                  call ph_vdwprof(icell,hw,shift,dvdw,.true.)
                  pnue  = pnue0
               case(0)
                  pnue  = pnue0
               end select labelp3

               imax=2000
               smax=splinemax(icell)
               if(smax <= 0.) then
                  write(6,*) 'PHOTON MODULE (PH_POST_ENERGY): WARNING'
                  write(6,*) '   splinemax for icell',icell,'<= 0.'
                  write(6,*) '   use of linecenter, isotropic'
                  e0=e00
                  vel=PH_POST_ANGLE(1)
               else
                  SPX => SPLINEX(icell,:)
                  SPA => SPLINEA(icell,:)
                  SPB => SPLINEB(icell,:)
                  SPC => SPLINEC(icell,:)
                  swma = eemax(icell)
                  swmi = eemin(icell)
                  n=volspecnum(ivs)
                  ii=0
                  iwarn=0
		  pnuehz=pnue/ev2hz
                  DO WHILE(ii<imax)

                     labelp4: select case(reaction%iprofiletype)
                     case(1)
                        zzep1=pnuehz
                     case(2)
                        zzep1=ph_sam_lorentz(gam,pnuehz)
                     case(3)
                        zzep1=ph_sam_lorentz(gam,pnuehz)
                     case(4)
                        iil=phphilphil2line(pil)
                        if(iil > 0) then
                           zzep1=ph_sam_lorvdw2(hw,shift,dvdw,pnuehz)
                        else
                           zzep1=pnuehz
                        endif
                     case(0)
                        zzep1=pnuehz
                     end select labelp4

                     val = PH_XQUAVAL(n,zzep1,iflag,icell)
                     if(iflag.ne.0) THEN
                        ilferr: select case(iflag)
                        case(3)
                           iwarn=iwarn+1
                        case default    
                           write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
                           write(6,*) '   photon sampling, iflag=',iflag
                           call exit(1)
                        end select ilferr
                     endif                

                     if(val < 0.) val = 0.

                     ZEP2=ranf_eirene()*smax                 
                     if(zep2 <= val) then
                        e0=zzep1
                        if(e0 < swmi.and.e0 > swma) cycle
                     endif
                     ii=ii+1

                  ENDDO
                  IF(ii>=imax) THEN
                     write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
                     write(6,*) 'MAXIMUM NUMBER OF REJECTION-SAMPLINGS'
                     write(6,*) 'REACHED, imax,ii=',imax,ii
                     call exit(1)
                  ENDIF
csw sample direction isotropically, set velx,vely...
                  vel=PH_POST_ANGLE(1)
               endif               
               return

            case(1)
c           n=1 followed           
c           use energy from n=2 testparticle
               cvrss1=cvrssa(iold)
               velx=vxo
               vely=vyo
               velz=vzo
               vel =vlo
               velq=vel*vel
               e0=cvrss1*velq
               return

            case default
               write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
               write(6,*) '   default case?!, check input!'
               call exit(1)
            end select stcs3
         case default
            write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
            write(6,*) '   default case?!, check input!'
            call exit(1)
         end select stcs
      end select
      
      write(6,*) 'PHOTON MODULE (PH_POST_ENERGY):'
      write(6,*) '   end of subroutine reached'
      write(6,*) '   there has to be an error, check input!'
      call exit(1)
      END SUBROUTINE PH_POST_ENERGY_OLD


      real(dp) FUNCTION PH_POST_ANGLE(mode) result(res)
      IMPLICIT NONE
      integer, intent(in) :: mode
c     not realy done, use volume routine
      res = PH_SAMVOL_ANGLE(mode)
      return
      END FUNCTION PH_POST_ANGLE

c PROFILE PARAMETERS                                                 

      SUBROUTINE PH_VDWPROF(icell,hw,shift,dvdw,lscale) 
      implicit none
      integer, intent(in) :: icell
      real(dp),intent(out) :: hw,shift,dvdw
      logical, intent(in) :: lscale
      integer :: pil,ipl,n1,n2, ip, iplti
      real(dp) :: nipl,ne,e00,e00s,dlst,dlqs,dlvw,te,tipl,dlres,l0,mipl
      real(dp) :: dlvws, dlqss
      hw=0.
      shift=0.
      dvdw=0.

      ipl=reaction%ignd
      e00=reaction%e0
      
c densities [m^-3], temperatures [eV]
      if(ipl < 0) then
         ipl=-ipl
         nipl=arg_DIWL*1.e6
         tipl=arg_tiwl
      else
         nipl=DIIN(ipl,icell)*1.e6
         tipl=tiin(mplsti(ipl),icell)
      endif
      ne  =DEIN(icell)*1.e6
      te  =tein(icell)

c mass ipl
      mipl=rmassp(ipl)

c broadening constants/shift [nm]-->[m]
      l0=hpcl/reaction%e0*1.E-2_dp

!      if(HERMANN_FLAG) then
C HERMANN:
!         shift=philips_sr(pil)*nipl + philips_ss(pil)*ne
!         hw   =philips_cr(pil)*nipl + philips_cs(pil)*ne
!         dvdw =philips_cw(pil)*nipl**2
!      else
C MARGARITA:
         dlres= l0*l0*PIA/2./clight*100.*reaction%c3*nipl

         dlst=(sqrt(8./PIA)*sqrt(te/pmasse)*cvel2a/100.)**(1./3.)
         dlst=dlst*l0*l0/4./PIA/clight*100.*11.37*
     .        reaction%c4**(2./3.)
         dlst=dlst*ne

         dlvws = 0._dp
         dlqss = 0._dp

         do ip=1,10
           ipl=reaction%iplsc6(ip)
           if (ipl > 0) then
             iplti=mplsti(ipl)
             nipl=DIIN(ipl,icell)*1.e6_dp
             tipl=tiin(iplti,icell)
             mipl=rmassp(ipl)
             dlvw=(sqrt(8._dp/PIA)*sqrt(tipl/(mipl/2._dp))*
     .            cvel2a/100._dp)**(3._dp/5._dp)
             dlvw=dlvw*l0*l0/4._dp/PIA/CLIGHT*100._dp*8.08_dp
     .            *reaction%c6a(ip)**(2._dp/5._dp)
             dlvw=dlvw*nipl

             dlvws = dlvws + dlvw

             dlqs = l0*l0/2._dp/PIA/clight*100._dp*
     .              (4._dp/3._dp*PIA)**2._dp*reaction%c6a(ip)
             dlqs = dlqs * nipl*nipl

             dlqss = dlqss + dlqs
           end if
         end do

         shift=-1.73*dlst -0.7*dlvws
         hw   =dlres + dlst + dlvws
         dvdw =dlqss
!      endif

c unit conversion: wavelength [m] --> [cm]
      shift=shift*100.
      hw=hw*100.
      dvdw=dvdw*100.

c transform wavelength [cm] --> energy [eV]
c dE / dlambda = -hc*lambda^-2 = -e^2/(hc)
      if (lscale) then
        shift=-shift*e00*e00/hpcl
c     e00s=e00-shift
        e00s=e00
        hw   =-hw*e00s*e00s/hpcl
        dvdw =-dvdw*e00s*e00s/hpcl
      end if

c      shift=0.

      return
      end subroutine ph_vdwprof

      subroutine ph_voigtprof(ipl,icell,dnd,gam)
      implicit none
      integer, intent(in) :: ipl,icell
      real(dp), intent(out) :: dnd,gam

      call ph_dopplerprof(ipl,icell,dnd)
      call ph_homprof(icell,gam)
      return
      end subroutine ph_voigtprof

      subroutine ph_dopplerprof(iipl,icell,dnd)
      implicit none
      integer, intent(in) :: iipl,icell
      real(dp), intent(out) :: dnd
      real(dp) :: e00,t,e1,e2
      integer :: n1,n2,ipl

      ipl=iipl
      e00=reaction%e0
      e1=reaction%e1
      e2=e00+e1
      if(ipl < 0) then
         ipl=-ipl
         t=arg_tiwl
      else
         t=tiin(mplsti(ipl),icell)     
      endif
      dnd=e00*sqrt(t)*rsqdvp(ipl)/clight
      return
      end subroutine ph_dopplerprof

      subroutine PH_HOMPROF(icell, gam)
      IMPLICIT NONE
      integer, intent(in) :: icell
      real(dp), intent(out) :: gam
      real(dp) :: gam1,gam2

      gam=0.
      
c natural
c      if(Lnatl(il)) then
         gam=gam+reaction%aik
c      endif

c resonant pressure broadening
!pb      if(LRESO(il)) then
c         gam=gam+ cr*diin(ipls,icell)
!pb      endif

c stark broadening
!pb      if(LSTRK(il)) then
c         gam=gam+ cs*diin(electrons,icell)
!pb      endif

c zeeman
c     if(LZEEMAN(il)) then
c     endif
      gam=gam/(2.*PIA)/ev2hz
      return
      END subroutine PH_HOMPROF

      subroutine PH_STRKPROF(icell, hw,shift)
      IMPLICIT NONE
      integer, intent(in) :: icell
      real(dp), intent(out) :: hw,shift
      real(dp) :: z,de,te,nn,nn2,nn3,nn4,eh,ec,sqeh,part1,part2,part3,
     .            lc,ve,fac,w
      real(dp) :: rd,rw,hw0

c
c stark-effect: for now, only hydrogen treated (linear stark)
c
      
      de = dein(icell)
      te = tein(icell)
c     upper state
      nn  = 2._dp
c compton length h^bar/m_e/c [cm]
      lc = 3.8616e-11_dp

csw
csw from Sobolev's Book, eq. 7.3.35
csw
c     debye radius, mean electron velocity, weisskopf radius
      rd = 7.43e2_dp*dsqrt(te/de)
      ve = dsqrt(8._dp/PIA)*4.19e7_dp*dsqrt(te)
      rw = dsqrt(2._dp/3._dp)*nn*nn/ve*lc*clight
      
      hw0= 32._dp/3._dp*de/ve*lc*lc*clight*clight
     .                 * (dlog(rd/rw)+.215_dp)*(1._dp+nn**4)
      hw0=hw0*hplnk/(2._dp*PIA)

csw
csw stark effect (Griem)      
csw
      eh = 13.606_dp
      z  = 1._dp
      nn2 = nn*nn
      nn3 = nn*nn*nn
      nn4 = nn2*nn2

c electron plasma frequency [eV]
      ec =8.98e3_dp*dsqrt(de)*hplnk
      
      sqeh=dsqrt(eh)

      part1 = te**1.5/((z-1._dp)*sqeh*ec + te*nn2*ec/z/sqeh)
      part1 = dlog(part1) *(3._dp*nn4-9._dp*nn2)

      part2 = nn4 + nn4*2._dp*te/eh/(1._dp+2._dp*te/eh+z*z/nn4)

      part3 = nn3/((z-1._dp)*z*z+te*nn2*z/eh)*(te/eh)**1.5+1.4_dp
      part3 = dlog(part3) *(nn4/3._dp + 17._dp*nn2/3._dp)
         
c mean velocity of electron sqrt(kT/me) [cm/s]
      ve = 4.19e7_dp*dsqrt(te)

      fac = dsqrt(2._dp*PIA)/z/z*lc*lc*clight*clight
      fac = fac / ve * de
         
      w  = hplnk * fac * (part1 + part2 + part3)/(2._dp*PIA)
c brd.const, shift [eV]
c      hw = w

      hw = hw0
      shift = hw*dsqrt(3._dp)/2._dp

      return
      end subroutine ph_strkprof

      REAL(dp) FUNCTION PH_LORENTZ(xx,gm) result(res)
      IMPLICIT NONE
      real(dp), intent(in) :: xx,gm
      res = (abs(gm)/(2.*PIA))/(xx**2 + (gm/2.)**2)
      return
      END FUNCTION PH_LORENTZ

      REAL(dp) FUNCTION PH_GAUSS(xx,dnd) result(res)
      IMPLICIT NONE
      real(dp), intent(in) :: xx,dnd
      res=exp(-(xx/dnd)**2)/(dnd*sqrt(PIA))
      return
      END FUNCTION PH_GAUSS

      complex(dp) function ph_faddeeva(x,y,dnd,icell,ipl) result(cres)
      implicit none
      real(dp), intent(in) :: x,y,dnd
      integer, intent(in) :: icell,ipl

      cres=ph_humlik(x,y)/(dnd*sqrt(pia))
      return
      end function ph_faddeeva

      complex(dp) function ph_faddeeva2(x,y,icell,ipl) result(cres)
      implicit none
      real(dp), intent(in) :: x,y
      integer, intent(in) :: icell,ipl
      real(dp) :: hnorm,u,v
      logical :: flag

      if(ipl == 0) then
         hnorm=1.
      else
         hnorm=faddeeva2norm(icell,ipl)
         if(hnorm <= 0.) hnorm=1.
         hnorm=1.
      endif

      !cres=ph_humlik(x,y)/hnorm

      call ph_wofz(x,y,u,v,flag)
      if(flag) then
         write(6,*) 'PHOTON MODULE (PH_FADDEEVA2):'
         write(6,*) '   flag=',flag
         call exit(1)
      endif
      cres=cmplx(u,v)/hnorm
      return
      end function ph_faddeeva2


      real(dp) function ph_lorvdw(de,hw,shift,dvdw,
     .                            icell,ipl) result(res)
      implicit none
      real(dp), intent(in) :: de,hw,shift,dvdw
      real(dp) :: a,b,ade,c,rlor,rvdw,lnorm,ahw,aa
      complex(dp) :: z1,z2,zz,sz1,fad,z,zz1,zz2,fad1,fad2,sz2
      integer, intent(in) :: icell,ipl

      if(ipl == 0) then
         lnorm=1.
      else
         !lnorm=lorvdwnorm1(icell,ipl)
         lnorm=lorvdwnorm2(icell,ipl)
         if(lnorm <= 0.) lnorm=1.
         lnorm=1.
      endif

      ahw = dabs(hw)

      a=(de-shift)/hw
      aa  = 1._dp/(1._dp+a*a)
      res = aa/(ahw*PIA)

      b = PIA*dvdw/4./hw
      c = dsqrt(abs(dvdw)) / (2._dp*PIA*ahw*dsqrt(ahw))

      z1= cmplx(-a, -1._dp) * aa
      zz = cdsqrt(b*z1) *(0._dp,1._dp)
      fad = ph_faddeeva2(dble(zz),aimag(zz),icell,0)

      sz1 = z1*cdsqrt(z1)
      z=sz1*fad
      rvdw = AIMAG(z)*c*PIA

c      z1= cmplx(-a, -1)/(1+a*a)
c      z2= cmplx(-a, +1)/(1+a*a)
c      zz1 = cdsqrt(-b*z1)
c      zz2 = cdsqrt(-b*z2)
c      sz1 = z1**1.5
c      sz2 = z2**1.5
c      fad1 = ph_faddeeva2(dble(zz1),aimag(zz1),icell,0)
c      fad2 = ph_faddeeva2(dble(zz2),aimag(zz2),icell,0)
c      z = sz1*fad1 - sz2*fad2
c      rvdw = dble(c*PIA*z/2.*(0.,1.))

      res=res+rvdw

      res=res/lnorm
      return
      end function ph_lorvdw                                                       w

      REAL(dp) FUNCTION PH_SAM_GAUSS(sig) result(res)
      IMPLICIT NONE
      real(dp), intent(in) :: sig
      real(dp) :: v1,v2,s,ar,f1
      do
         v1=2.*RANF_EIRENE()-1.
         v2=2.*RANF_EIRENE()-1.
         s=v1*v1+v2*v2
         if(s < 1.) exit            
      enddo
      ar=log(s)
      f1=v1*sqrt(-(ar+ar)/s)*sig
c     f2=v2*sqrt(-(ar+ar)/s)*sig
      res=f1
      END FUNCTION PH_SAM_GAUSS


      REAL(dp) FUNCTION PH_SAM_VOIGT(alph,sig) result(res)
      IMPLICIT NONE
      real(dp),intent(in) :: alph,sig
      real(dp) :: xx,f1,v1,v2,s,f2,ar
c lorentz
      xx = PIHA*(RANF_EIRENE()*2.-1.)
      f1=alph*tan(xx)
c gauss
      do 
         v1=2.*RANF_EIRENE()-1.
         v2=2.*RANF_EIRENE()-1.
         s=v1*v1+v2*v2
         if(s <= 1.) exit
      enddo
      ar=log(s)
      f2=v1*sqrt(-(ar+ar)/s)*sig
c gauss*lorentz:
      res=f1+f2
      return
      END FUNCTION PH_SAM_VOIGT
c
c
c
      SUBROUTINE PH_WOFZ (XI, YI, U, V, FLAG)
C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
C
C  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
C  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
C  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
C  MEANS SQRT(-1).
C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
C  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
C  OF THE FUNCTION.
C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
C
C
C  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY AR/POE DEFINED WILL 
C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
C
C
C  PARAMETER LIST
C     XI     = REAL      PART OF Z
C     YI     = IMAGINARY PART OF Z
C     U      = REAL      PART OF W(Z)
C     V      = IMAGINARY PART OF W(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
C  PUT TO 0 UPON UNDERFLOW;
C
C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
C  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
C
*
*
*
*
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      LOGICAL :: A, B
      LOGICAL, intent(inout) :: FLAG
      real(dp),PARAMETER   :: FACTOR   = 1.12837916709551257388D0,
     *                        RMAXREAL = 1.34078079299426D+154,
     *                        RMAXEXP  = 709.089565710919D0,
     *                        RMAXGONI = 3.53711887601422D+15
c     *                        RMAXREAL = 0.5D+154,
c     *                        RMAXEXP  = 708.503061461606D0,
c     *                        RMAXGONI = 3.53711887601422D+15
      real(dp), intent(in) :: xi,yi
      real(dp), intent(out):: u,v
      real(dp) :: xabs,yabs,x,y,qrho,xabsq,xquad,yquad,xsum,ysum,
     .            xaux,u1,v1,daux,u2,v2,h,rx,ry,sx,sy,tx,ty,c,qlambda,
     .            w1,h2
      integer :: n,j,i,kapn,nu,np1
*
      FLAG = .FALSE.
*
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
*
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
      QRHO = X**2 + Y**2
*
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
*
      A     = QRHO.LT.0.085264D0
*
      IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
 10     CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
*
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
*
      ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
*
        B = (H.GT.0.0)
*
        IF (B) QLAMBDA = H2**KAPN
*
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
*
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
 11     CONTINUE
*
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
*
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
*
      END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
      IF (YI.LT.0.0) THEN
*
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
*
C
C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
*
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
*
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
*
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END SUBROUTINE PH_WOFZ

      REAL(dp) FUNCTION PH_SAM_LORENTZ(alph,totshift) result(res)
      IMPLICIT NONE
      real(dp), intent(in) :: alph, totshift
      real(dp) :: xx
      do
        xx=PIHA*(RANF_EIRENE()*2._DP-1._DP)
        res=alph*tan(xx)+totshift
        if (res > 0._dp) exit
      end do
      return
      END FUNCTION PH_SAM_LORENTZ

      REAL(dp) FUNCTION PH_SAM_LORVDW2(hw,shift,dvdw,l0) 
     .         result(res)
      IMPLICIT NONE
      real(dp), intent(in) :: hw,shift,dvdw,l0
      real(dp) :: xx,sig,z1,z2,mx,beta
      integer :: i,imax

c lorentz  
      z1 = ph_sam_lorentz(hw,shift+l0)
      z2=0.

c van der waals only for red wing of line
c      if(z1 < 0.) then

         beta=pia*0.25_dp*dvdw
         sig=1._dp/sqrt(2._dp*beta)
         do
           xx=ph_sam_gauss(sig)
           if (abs(xx) > eps60) exit
         end do 
         z2=1._dp / (xx*xx)

c      endif

      res=z1+z2
      return
      end function ph_sam_lorvdw2
c
c
c
      complex(dp) FUNCTION PH_HUMLIK(x,y) result(cres)
      IMPLICIT NONE
c arguments
      real(dp), intent(in) :: x,y
c To calculate the FADDEEVA function with relative error less than 10^(-R).
c R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the user
c subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
      REAL(dp) :: K,L
      real(dp), PARAMETER :: R0 = 146.7, R1 = 14.67 ! for R=4, region boundaries
            
c Constants
      real(dp), PARAMETER :: RRTPI = 0.56418958     ! 1/sqrt(pi)
      real(dp), PARAMETER :: Y0 = 1.5, Y0PY0 = Y0+Y0, Y0Q = Y0*Y0 ! for cpf12 algor.
      REAL(dp), save :: C(0:5), S(0:5), T(0:5)
c SAVE preserves values of C, S and T (static) arrays between procedure calls

      DATA C / 1.0117281,     -0.75197147,        0.012557727,
     .     0.010022008,   -0.00024206814,     0.00000050084806 /
      DATA S / 1.393237,       0.23115241,       -0.15535147,
     .     0.0062183662,   0.000091908299,   -0.00000062752596 /
      DATA T / 0.31424038,     0.94778839,        1.5976826,
     .     2.2795071,      3.0206370,         3.8897249 /
      
c Local variables
      INTEGER :: I, J              ! Loop variables
      INTEGER :: RG1, RG2, RG3     ! y polynomial flags
      REAL(dp) :: ABX, XQ, YQ, YRRTPI ! |x|, x^2, y^2, y/SQRT(pi)
      REAL(dp) :: XLIM0, XLIM1, XLIM2, XLIM3, XLIM4 ! |x| on region boundaries
      REAL(dp) :: A0, D0, D2, E0, E2, E4, H0, H2, H4, H6 ! W4 temporary variables
      REAL(dp) :: P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
      real(dp) :: b1,f1,f3,f5,q1,q3,q5,q7
      REAL(dp) :: XP(0:5), XM(0:5), YP(0:5), YM(0:5) ! CPF12 temporary values
      REAL(dp) :: MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
      REAL(dp) :: D, YF, YPY0, YPY0Q  
      
c**** Start of executable code *****************************************
      
      RG1 = 1                   ! Set flags
      RG2 = 1
      RG3 = 1
      YQ  = Y*Y                 ! y^2
      YRRTPI = Y*RRTPI          ! y/SQRT(pi)
      
c Region boundaries when both K and L are required or when R<>4 
      XLIM0 = R0 - Y 
      XLIM1 = R1 - Y
      XLIM3 = 3.097*Y - 0.45
      
      XLIM2 = 6.8 - Y
      XLIM4 = 18.1*Y + 1.65
      IF ( Y .LE. 0.000001 ) THEN ! When y<10^-6
         XLIM1 = XLIM0          ! avoid W4 algorithm
         XLIM2 = XLIM0
      ENDIF
c.....
      ABX = ABS ( X )           ! |x|
      XQ  = ABX*ABX             ! x^2
      IF ( ABX .GT. XLIM0 ) THEN ! Region 0 algorithm
         K = YRRTPI / (XQ + YQ)
         L = k*x/y        
      ELSEIF ( ABX .GT. XLIM1 ) THEN ! Humlicek W4 Region 1
         IF ( RG1 .NE. 0 ) THEN ! First point in Region 1
            RG1 = 0
            A0 = YQ + 0.5       ! Region 1 y-dependents
            D0 = A0*A0
            D2 = YQ + YQ - 1.0
            b1 = yq - 0.5
         ENDIF
         D = RRTPI / (D0 + XQ*(D2 + XQ))
         K = D*Y * (A0 + XQ)
         L = d*x * (b1 + xq)         
      ELSEIF ( ABX .GT. XLIM2 ) THEN ! Humlicek W4 Region 2 
         IF ( RG2 .NE. 0 ) THEN ! First point in Region 2
            RG2 = 0
            H0 =  0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ))) ! Region 2 y-dependents
            H2 = -4.5    + YQ*(9.0 + YQ*( 6.0 + YQ* 4.0))
            H4 = 10.5    - YQ*(6.0 - YQ*  6.0)
            H6 = -6.0    + YQ* 4.0
            E0 =  1.875  + YQ*(8.25 + YQ*(5.5 + YQ))
            E2 =  5.25   + YQ*(1.0  + YQ* 3.0)
            E4 =  0.75*H6
            f1 = -1.875 + yq*(5.25 + yq*(4.5 + yq))
            f3 =  8.25  - yq*(1.0  - yq* 3.0)
            f5 = -5.5   + yq* 3.0
         ENDIF
         D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
         K = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)))
         L = d*x   *(f1 + xq*(f3 + xq*(f5 + xq)))         
      ELSEIF ( ABX .LT. XLIM3 ) THEN ! Humlicek W4 Region 3
         IF ( RG3 .NE. 0 ) THEN ! First point in Region 3
            RG3 = 0
            Z0 = 272.1014     + Y*(1280.829 + Y*(2802.870 + Y*(3764.966 ! Region 3 y-dependents
     &         + Y*(3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989
     &         + Y*(88.26741 + Y*(13.39880 + Y)))))))))
            Z2 = 211.678      + Y*(902.3066 + Y*(1758.336 + Y*(2037.310
     &           + Y*(1549.675 + Y*(793.4273 + Y*(266.2987
     &           + Y*(53.59518 + Y*5.0)))))))
            Z4 = 78.86585     + Y*(308.1852 + Y*(497.3014 + Y*(479.2576
     &           + Y*(269.2916 + Y*(80.39278 + Y*10.0)))))
            Z6 = 22.03523     + Y*(55.02933 + Y*(92.75679 + Y*(53.59518
     &           + Y*10.0)))
            Z8 = 1.496460     + Y*(13.39880 + Y*5.0)
            P0 = 153.5168     + Y*(549.3954 + Y*(919.4955 + Y*(946.8970
     &          + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941
     &          + Y*(4.264678 + Y*0.3183291))))))))
            P2 = -34.16955    + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730
     &           + Y*(139.4665 + Y*(56.81652 + Y*(12.79458
     &           + Y*1.2733163))))))
            P4 = 2.584042     + Y*(10.46332 + Y*(24.01655 + Y*(29.81482
     &           + Y*(12.79568 + Y*1.9099744))))
            P6= -0.07272979  + Y*(0.9377051+ Y*(4.266322 + Y*1.273316))
            P8 = 0.0005480304 + Y*0.3183291
            q1 = 173.2355  + y*(508.2585 + y*(685.8378 + y*(557.5178
     .                     + y*(301.3208 + y*(111.0528 + y*(27.62940
     .                     + y*(4.264130 + y*0.3183291)))))))
            q3 = 18.97431  + y*(100.7375 + y*(160.4013 + y*(130.8905
     .                     + y*(55.88650 + y*(12.79239+y*1.273316)))))
            q5 = 7.985877  + y*(19.83766 + y*(28.88480 + y*(12.79239
     .                     + y*1.909974)))
            q7 = 0.6276985 + y*(4.264130 + y*1.273316)            
         ENDIF
         D =1.7724538 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))))
         K = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
         L = d*x*(q1+xq*(q3+xq*(q5+xq*(q7 + xq*0.3183291))))         
      ELSE                      ! Humlicek CPF12 algorithm
         YPY0 = Y + Y0
         YPY0Q = YPY0*YPY0
         K = 0.0
         L = 0.0
         DO J = 0, 5
            D = X - T(J)
            MQ(J) = D*D
            MF(J) = 1.0 / (MQ(J) + YPY0Q)
            XM(J) = MF(J)*D
            YM(J) = MF(J)*YPY0
            D = X + T(J)
            PQ(J) = D*D
            PF(J) = 1.0 / (PQ(J) + YPY0Q)
            XP(J) = PF(J)*D
            YP(J) = PF(J)*YPY0
            L=L+c(j)*(xm(j)+xp(j)) + s(j)*(ym(j)-yp(j))
         ENDDO
         
         IF ( ABX .LE. XLIM4 ) THEN ! Humlicek CPF12 Region I
            DO J = 0, 5
               K = K + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
            ENDDO            
         ELSE                   ! Humlicek CPF12 Region II
            YF   = Y + Y0PY0
            DO J = 0, 5
               K = K
     &     + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)) / (MQ(J)+Y0Q)
     &     + (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF*XP(J)) / (PQ(J)+Y0Q)
            ENDDO
            K = Y*K + EXP ( -XQ )
         ENDIF
      ENDIF
*.....
      cres = CMPLX(K,L)
      RETURN
      END FUNCTION PH_HUMLIK

      REAL(dp) FUNCTION PH_XQUAVAL(n,v,iflag,ic) result(res)
c ic not needed
      IMPLICIT NONE
      integer, intent(in) :: n,ic
      integer, intent(out):: iflag
      real(dp), intent(inout) :: v
      real(dp) :: dx
      integer :: l,k,i
      data i/1/

      IFLAG=0
      IF(N.LT.2) THEN
         IFLAG=1
         RETURN
      ENDIF

c     inline PH_INTONE
      IF(I.GE.N) I=1
      IF(V.LT.SPX(1).OR.V.GT.SPX(N)) THEN
         IFLAG=3
         RETURN
      ENDIF

      IF(V.LT.SPX(I)) GOTO 10
      IF(V.LE.SPX(I+1)) GOTO 40
      L=N
      GOTO 30
 10   L=I
      I=1
 20   K=(I+L)/2
      IF(V.LT.SPX(K)) THEN
         L=K
      ELSE
         I=K
      ENDIF
 30   IF(L.GT.I+1) GOTO 20
c     end inline PH_INTONE

 40   IF(IFLAG.NE.0) RETURN
      DX=V-SPX(I)
      res=SPA(I)+DX*(SPB(I)+SPC(I)*DX)
      RETURN
      END FUNCTION PH_XQUAVAL

c EIRENE utilities:

      SUBROUTINE PH_ALLOC_XSECTPH(nnrot)
      IMPLICIT NONE
      integer, intent(in) :: nnrot
      integer :: iphot,nrc,kk
c allocate
      if(nnrot > 0) then
         allocate(PHV_LGPHOT(0:nphot,0:nnrot,0:5))
         phv_nrotph=nnrot
         phv_lgphot=0
      else
         allocate(PHV_LGPHOT(0:nphot,0:0,0:5))
         phv_nrotph=0
         phv_lgphot=0
      endif
      allocate(phv_nphoti(nphot))
      phv_nphoti=0
      if(.not.allocated(phv_iestotph)) then
         allocate(phv_iestotph(0:nphot,nnrot,3))
         phv_iestotph=0
      endif
      if(.not.allocated(phv_n1stotph)) then
         allocate(phv_n1stotph(0:nphot,nnrot,3))
         allocate(phv_n2ndotph(0:nphot,nnrot,3))
         phv_n1stotph=0
         phv_n2ndotph=0
      endif
      return
      END SUBROUTINE PH_ALLOC_XSECTPH

      SUBROUTINE PH_XSECTPH(ipht,nrc,idsc)
      IMPLICIT NONE
      integer, intent(in) :: ipht,nrc,idsc,ipl
      integer :: kk,ipl0,ipl1,ipl2,ityp0,ityp1,ityp2,il,n0,n1,n2,
     .    nh,nl,iid,ifnd,mode,updf, j, nseot4, ierr, ipl0ti
      real(dp) :: factkk, ebulk

      kk=ireacph(ipht,nrc)
      factkk=freacph(ipht,nrc)
      if(factkk == 0.) factkk=1.

      IPL0 =IDEZ(IBULKPH(ipht,nrc),3,3)
      IPL1 =IDEZ(ISCD1PH(ipht,nrc),3,3)
      IPL2 =IDEZ(ISCD2PH(ipht,nrc),3,3)
      ITYP0=IDEZ(IBULKPH(ipht,nrc),1,3)
      ITYP1=IDEZ(ISCD1PH(ipht,nrc),1,3)
      ITYP2=IDEZ(ISCD2PH(ipht,nrc),1,3)

c collect niveau information data, N0, N1, N2
      
      updf=1
      mode=0
      ifnd=999
      
      PHV_LGPHOT(ipht,idsc,0)=idsc
      PHV_LGPHOT(ipht,idsc,1)=ipl0
      PHV_LGPHOT(ipht,idsc,2)=ifnd
      PHV_LGPHOT(ipht,idsc,3)=kk
      PHV_LGPHOT(ipht,idsc,4)=updf
      PHV_LGPHOT(ipht,idsc,5)=mode

CDR  1ST SECONDARY
      PHV_N1STOTph(ipht,idsc,1) = ityp1
      PHV_N1STOTph(ipht,idsc,2) = ipl1
      PHV_N1STOTph(ipht,idsc,3) = 0
      IF (ityp1 < 4) 
     .  PHV_N1STOTph(ipht,idsc,3) = IDEZ(ISCD1PH(ipht,nrc),2,3)
CDR  2ND SECONDARY
      PHV_N2NDOTph(ipht,idsc,1) = ityp2
      PHV_N2NDOTph(ipht,idsc,2) = ipl2
      PHV_N2NDOTph(ipht,idsc,3) = 0
      IF (ityp2 < 4) 
     .  PHV_N2NDOTph(ipht,idsc,3) = IDEZ(ISCD2PH(ipht,nrc),2,3)
      

cdr  CROSS SECTION: HERE: BEAM-BEAM, NO DOPPLER FROM THERMAL MOTION
      MODCOL(7,1,IPHT,IPL0)=KK
cdr  COLLISION MODEL 4: BEAM-BEAM
      MODCOL(7,2,IPHT,IPL0)=4
c
c     DEFCX(IRCX)=LOG(CVELI2*PMASS)
c     EEFCX(IRCX)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4. BULK ION ENERGY LOSS RATE
C
cdr   NSEOT4=IDEZ(ISCDE,4,5)
cdr bulk energy loss not ready
      nseot4=0
      ebulk=0.
cdr
      IF (NSEOT4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK ION = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        write (6,*) ' in ph_xsectph, nseot4=0 '
        IF (EBULK.LE.0.D0) THEN
          write (6,*) ' in ph_xsectph, nseot4=0, ebulk <= 0 '
          IF (NSTORDR >= NRAD) THEN
            write (6,*) ' in ph_xsectph, nstordr>nrad '
            IPL0TI=MPLSTI(IPL0)
            DO J=1,NSBOX
              EPLOT3(Idsc,J,1)=1.5*TIIN(IPL0TI,J)+EDRIFT(IPL0,J)
            ENDDO
            NELROT(Idsc) = -3
          ELSE
            write (6,*) ' in ph_xsectph, nstordr<nrad '
            NELROT(Idsc) = -3
          END IF
        ELSE
          write (6,*) ' in ph_xsectph, nseot4=0, ebulk > 0 '
          IF (NSTORDR >= NRAD) THEN
            write (6,*) ' in ph_xsectph, nstordr>nrad '
            DO 251 J=1,NSBOX
              EPLOT3(Idsc,J,1)=EBULK+EDRIFT(IPL0,J)
251         CONTINUE
            NELROT(Idsc) = -2
          ELSE
            NELROT(Idsc) = -2
            EPLOT3(Idsc,1,1)=EBULK
            write (6,*) ' in ph_xsectph, nstordr<nrad '
          END IF
        ENDIF
        MODCOL(7,4,IPHT,IPL0)=3
        write (6,*) ' in ph_xsectph, Modcol(7,4,1,1) ',
     .                MODCOL(7,4,IPHT,IPL0)
      ELSEIF (NSEOT4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        write (6,*) ' in ph_xsectph, nseot4=1 '
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            IPL0TI=MPLSTI(IPL0)
            DO 252 J=1,NSBOX
              EPLOT3(Idsc,J,1)=1.5*TIIN(IPL0TI,J)+EDRIFT(IPL0,J)
  252       CONTINUE
            NELROT(Idsc) = -3
          ELSE
            NELROT(Idsc) = -3
          END IF
        ELSE
          WRITE (6,*) 'WARNING FROM SUBR. xsectph '
          WRITE (6,*) 'MODIFIED TREATMENT OF photon collision '
          WRITE (6,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (6,*) 'RATHER THEN WITH T = TIIN '
          CALL LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLOT3(Idsc,J,1)=EBULK+EDRIFT(IPL0,J)
2511        CONTINUE
            NELROT(Idsc) = -2
          ELSE
            NELROT(Idsc) = -2
            EPLOT3(Idsc,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(7,4,IPHT,IPL0)=1
C     ELSEIF (NSECX4.EQ.2) THEN
C  use i-integral expressions. to be written
c     ELSEIF (NSECX4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
      ELSE
        IERR=5
        GOTO 996
      ENDIF

C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION

      PHV_IESTOTph(ipht,idsc,1) = IDEZ(IESTMPH(IPHT,idsc),1,3)
      PHV_IESTOTph(ipht,idsc,2) = IDEZ(IESTMPH(IPHT,idsc),2,3)
      PHV_IESTOTph(ipht,idsc,3) = IDEZ(IESTMPH(IPHT,idsc),3,3)
C
cdr  not ready
c
c     ITYP1=N1STX(IRCX,1)
c     ITYP2=N2NDX(IRCX,1)
c     IF (IESTCX(IRCX,1).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       CALL LEER(1)
c       WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
c       WRITE (6,*) 'IRCX = ',IRCX
c       WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       IESTCX(IRCX,1)=0
c     ENDIF
c     IF (IESTCX(IRCX,2).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       CALL LEER(1)
c       WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
c       WRITE (6,*) 'IRCX = ',IRCX
c       WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       IESTCX(IRCX,2)=0
c     ENDIF
c     IF (IESTCX(IRCX,3).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       CALL LEER(1)
c       WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR EN.-BALANCE '
c       WRITE (6,*) 'IRCX = ',IRCX
c       WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       IESTCX(IRCX,3)=0
c     ENDIF
      RETURN
C
      ENTRY XSTPH_2(Idsc,IPL)
C
c     CALL LEER(1)
c     WRITE (6,*) 'Photon REACTION NO. Idsc= ',Idsc
c     CALL LEER(1)
c     WRITE (6,*) 'Collision WITH BULK IONS IPLS:'
c     WRITE (6,*) '1ST AND 2ND NEXT GEN. SPECIES I2ND1, I2ND2:'
c     ITYP1=N1STX(IRCX,1)
c     ITYP2=N2NDX(IRCX,1)
c     ISPZ1=N1STX(IRCX,2)
c     ISPZ2=N2NDX(IRCX,2)
c     IF (ITYP1.EQ.1) TEXTS1=TEXTS(NSPH+ISPZ1)
c     IF (ITYP1.EQ.2) TEXTS1=TEXTS(NSPA+ISPZ1)
c     IF (ITYP1.EQ.3) TEXTS1=TEXTS(NSPAM+ISPZ1)
c     IF (ITYP1.EQ.4) TEXTS1=TEXTS(NSPAMI+ISPZ1)
c     IF (ITYP2.EQ.1) TEXTS2=TEXTS(NSPH+ISPZ2)
c     IF (ITYP2.EQ.2) TEXTS2=TEXTS(NSPA+ISPZ2)
c     IF (ITYP2.EQ.3) TEXTS2=TEXTS(NSPAM+ISPZ2)
c     IF (ITYP2.EQ.4) TEXTS2=TEXTS(NSPAMI+ISPZ2)
c     WRITE (6,*) 'IPLS= ',TEXTS(NSPAMI+IPL),'I2ND1= ',TEXTS1,
c    .                    'I2ND2= ',TEXTS2
c     CALL LEER(1)
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph. EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR OT '
      CALL EXIT(1)
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
!pb   WRITE (6,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
      CALL EXIT(1)
993   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph: EXIT CALLED  '
      WRITE (6,*) 'EBULK_ION .LE.0, BUT MONOENERGETIC DISTRIBUTION?'
      WRITE (6,*) 'CHECK ENERGY FLAG ISCDEA'
!pb   WRITE (6,*) 'KK,IATM,IPLS,ISCDEA ',KK,IATM,IPLS,ISCDEA
      CALL EXIT(1)
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSectph: EXIT CALLED  '
      WRITE (6,*) 'NO CROSS SECTION AVAILABLE FOR NON DEFAULT OT'
      WRITE (6,*) 'KK,IPHT,IPL0 ',KK,IPHT,IPL0
      WRITE (6,*) 'EITHER PROVIDE CROSS SECTION OR USE DIFFERENT '
      WRITE (6,*) 'POST COLLISION SAMPLING FLAG ISCDEA'
      CALL EXIT(1)
      return
      END SUBROUTINE PH_XSECTPH

      SUBROUTINE PH_ALLOC_XSECTA(nnrot)
      IMPLICIT NONE
      integer, intent(in) :: nnrot
      if(allocated(phv_lgaot)) deallocate(phv_lgaot)
      if(allocated(phv_naoti)) deallocate(phv_naoti)
      if(allocated(phv_iestotat)) deallocate(phv_iestotat)
      if(allocated(phv_n1stotat)) deallocate(phv_n1stotat)
      if(allocated(phv_n2ndotat)) deallocate(phv_n2ndotat)
c allocate
      IF(NNROT > 0) then
         allocate(PHV_LGAOT(0:natm,0:nnrot,0:5))
         PHV_NROTA=nnrot
         phv_lgaot=0
      else
         allocate(PHV_LGAOT(0:natm,0:0,0:5))
         PHV_NROTA=0
         phv_lgaot=0         
      endif
      allocate(PHV_NAOTI(natm))
      phv_naoti=0

      allocate(phv_iestotat(0:natm,nnrot,3))
      phv_iestotat=0

      allocate(phv_n1stotat(0:natm,nnrot,3))
      allocate(phv_n2ndotat(0:natm,nnrot,3))
      phv_n1stotat=0
      phv_n2ndotat=0
      END SUBROUTINE PH_ALLOC_XSECTA

      SUBROUTINE PH_XSECTP(ipls,nrc,idsc,irrc)
      IMPLICIT NONE
      integer, intent(in) :: ipls,nrc,idsc,irrc
      integer :: kk,itp,ispz,ipl0,ipl1,ipl2,ityp0,ityp1,ityp2,
     .    np,n0,n1,n2,il,nh,nl,ifnd,j,iadd
      real(dp) :: factkk,aik

      kk = IREACP(ipls,nrc)
      if (kk /= idreac) call get_reaction(kk)
      FACTKK=FREACP(IPLS,NRC)      
      IF (FACTKK.EQ.0.D0) FACTKK=1.
      aik=reaction%aik

      IF (irrc.GT.NREC) then
         write(6,*) 'PHOTON MODULE (PH_XSECTP): NREC too small'
         call exit(1)
      endif
      if(.not.allocated(phv_lgprc)) then
         allocate(phv_lgprc(0:npls,0:nrec))
      endif
      
      IPL0 =IDEZ(IBULKP(ipls,nrc),3,3)
      IPL1 =IDEZ(ISCD1P(ipls,nrc),3,3)
      IPL2 =IDEZ(ISCD2P(ipls,nrc),3,3)
      ITYP0=IDEZ(IBULKP(ipls,nrc),1,3)
      ITYP1=IDEZ(ISCD1P(ipls,nrc),1,3)
      ITYP2=IDEZ(ISCD2P(ipls,nrc),1,3)

      LGPRC(IPLS,IDSC)=IRRC

      ITP=IDEZ(ISCD1P(IPLS,NRC),1,3)
      ISPZ=IDEZ(ISCD1P(IPLS,NRC),3,3)

      facrea(kk) = factkk
      NREARC(irrc) = kk
      do j=1,nrad
         tabrc1(irrc,j)=aik*factkk
      enddo
      modcol(6,2,nspami+ipls,1)=1

      select case(itp)
      case(0)
         NPHPRC(IRRC)=ISPZ
      case(1)
         NATPRC(IRRC)=ISPZ
      case default
         write(6,*) 'PHOTON MODULE (PH_XSECTP):'
         write(6,*) '   itp=',itp,' not allowed'
         call exit(1)
      end select
      return
      END SUBROUTINE PH_XSECTP
csw
      SUBROUTINE PH_UPDATE_SURF_SPECTRUM(wt,ind)
      IMPLICIT NONE
      real(dp), intent(in) :: wt
      integer, intent(in) :: ind
      integer :: il,iss,numd,surf,j,i
      real(dp) :: smin,smax,dele,delta,e
      real(dp), pointer :: dummyptr(:,:)

      if(ityp /= 0) return
      select case(ind)
      case(2)
c emitted
         return
      case(1)
c incident
      case default
         write(6,*) 'PHOTON MODULE (PH_UPDATE_SURF_SPECTRUM):'
         write(6,*) '    ind /= 1,2   exit.'
         call exit(1)
      end select

      il=actual_linenum
      if(il ==0) return

      iss=line2srfspec(il)
      if(iss == 0) return
      numd=srfspecnum(iss)
      if(numd == 0) return

      smin=srfspecmin(iss)
      smax=srfspecmax(iss)
      surf=srfspecsurf(iss)

      delta=smax-smin
      dele=delta/numd
      dummyptr => ptr_srfspec(iss)%p
      if    (e0 <= smin) then
         j=0
      elseif(e0 >= smax) then
         j=numd+1
      else
         if(.not.allocated(ecdelta)) then
            j=1+INT(dble(numd)*(e0-smin)/delta)
         else
            j=-1
            do i=phv_ecnum,1,-1
               e=ec(i)
               if(e0 >= ec(i) .and. e0 < ec(i)+ecdelta(i)) then
                  j=i
                  dele=ecdelta(i)
                  exit
               endif
            enddo
         endif
      endif         

      if(j>= 0) dummyptr(surf,j) = dummyptr(surf,j)+ wt*CLIGHT*e0/dele

      return
      END SUBROUTINE PH_UPDATE_SURF_SPECTRUM
c
      subroutine ph_integrate(istr, scale, fpht,fatm,fmol,fion)
!pb      use parmmod
      implicit none
      real(dp), intent(in) :: scale,fpht,fatm,fmol,fion
      integer, intent(in) :: istr
      if(bgk_num > 0) then
         call ph_integrate_bgk(istr,scale,fpht,fatm)
      endif
!pb as far as I know spectra are done by calc_spectrum
!pb      call ph_integrate_spectra(istr,scale,fpht)
      return
      end subroutine ph_integrate
c
      SUBROUTINE PH_INTEGRATE_BGK(istr,scale,fpht,fatm)
      implicit none
      real(dp), intent(in) :: scale, fpht,fatm
      integer, intent(in) :: istr
      real(dp) :: yp1,yp2,yp3, ya1,ya2,ya3
      integer :: nc,i,j,ibgk

c bgk-evaluation after all strata have been finished
      if(istr /= 0 .or. bgk_num == 0) return

      do ibgk=1,bgk_num
         yp1=0.
         yp2=0.
         yp3=0.
         ya1=0.
         ya2=0.
         ya3=0.
c scale # per cell
         bgk_tallies(ibgk,:,:,:)=bgk_tallies(ibgk,:,:,:)*scale
c integrate
         do nc =1,nrad
            yp1=yp1+bgk_tallies(ibgk,0,1,nc)
            yp2=yp2+bgk_tallies(ibgk,0,2,nc)
            yp3=yp3+bgk_tallies(ibgk,0,3,nc)
            ya1=ya1+bgk_tallies(ibgk,1,1,nc)
            ya2=ya2+bgk_tallies(ibgk,1,2,nc)
            ya3=ya3+bgk_tallies(ibgk,1,3,nc)
         enddo         
         bgk_tallies(ibgk,0,1,0) = yp1
         bgk_tallies(ibgk,0,2,0) = yp2
         bgk_tallies(ibgk,0,3,0) = yp3
         bgk_tallies(ibgk,1,1,0) = ya1
         bgk_tallies(ibgk,1,2,0) = ya2
         bgk_tallies(ibgk,1,3,0) = ya3
      enddo

      return
      end subroutine PH_INTEGRATE_BGK
c
      SUBROUTINE PH_INTEGRATE_SPECTRA(istr,scale,fpht)
      implicit none
      integer, intent(in) :: istr
      real(dp), intent(in) :: scale,fpht
      integer :: i,num,nc,j,surf,il,ivs
      real(dp),allocatable :: dummy(:)
      real(dp) :: yint,smax,smin,dele,fact,delta,delee
      real(dp),pointer :: dummyptr(:,:)

      fact=scale*fpht/clight
      il=actual_linenum
      if(il == 0) return

C VOLUMES, all
!pb      if(istr == 0) then
      if((istr == 0) .and. (numvolspecs > 0)) then
c global spectrum
         num=volspecnum(0)
         smax=volspecmax(0)
         smin=volspecmin(0)
         delta=smax-smin
         dele=delta/num
         allocate(dummy(-1:num+1))        
         dummyptr => ptr_volspec(0)%p
         if(.not.associated(dummyptr)) then
            write(6,*) 'PHOTON MODULE (PH_INTEGRATE_SPECTRA):'
            write(6,*) '    pointer not associated for iv=',0
            call exit(1)
         endif

         do nc =1,nrad
c scale # per cell
            dummy(0:num+1) = dummyptr(nc,0:num+1)
c integrate here!!!, include JUNK tallies
            yint=0.
            do j=0,num+1
               yint=yint+dummy(j)*dele
            enddo
            dummy(-1) = yint

            dummyptr(nc,-1:num+1)=dummy(-1:num+1)*fact
         enddo         
         if(allocated(dummy)) deallocate(dummy)
         return
      endif

      if (istr == 0) return

C VOLUMES
      do i=1,numvolspecs
         if(volspeclne(i) /= il) cycle
         num=volspecnum(i)
c         smax=volspecmax(i)
c         smin=volspecmin(i)
c         dele=(smax-smin)/num
         allocate(dummy(-1:num+1))
         dummyptr => ptr_volspec(i)%p
         if(.not.associated(dummyptr)) then
            write(6,*) 'PHOTON MODULE (PH_INTEGRATE_SPECTRA):'
            write(6,*) '    pointer not associated for ivs=',i
            call exit(1)
         endif
         
         do nc =1,nrad
            delta=eemax(nc)-eemin(nc)
            dele=delta/num
c scale # per cell
            dummy(0:num+1) = dummyptr(nc,0:num+1)
c integrate here!!!, include JUNK tallies
            yint=0.
            do j=0,num+1
               yint=yint+dummy(j)*dele
            enddo
            dummy(-1)=yint

            dummyptr(nc,-1:num+1)=dummy(-1:num+1)*fact
         enddo         
         if(allocated(dummy)) deallocate(dummy)
      enddo

C SURFACES
      do i=1,numsrfspecs
         if(srfspeclne(i) /= il) cycle
         num=srfspecnum(i)
         smax=srfspecmax(i)
         smin=srfspecmin(i)
         delta=smax-smin
         dele=delta/num
         surf=srfspecsurf(i)
         allocate(dummy(-1:num+1))
         dummyptr => ptr_srfspec(i)%p
         if(.not.associated(dummyptr)) then
            write(6,*) 'PHOTON MODULE (PH_INTEGRATE_SPECTRA):'
            write(6,*) '    pointer not associated for iss=',i
            call exit(1)
         endif
         
c scale # per cell
         dummy(0:num+1) = dummyptr(surf,0:num+1)
c integrate here!!!, include JUNK tallies
         yint=0.
         do j=0,num+1
            if(j>0 .and. j<= num .and. allocated(ecdelta)) then
               delee = ecdelta(j)
            else
               delee = dele
            endif

            yint=yint+dummy(j)*delee
         enddo
         dummy(-1)=yint

         dummyptr(surf,-1:num+1)=dummy(-1:num+1)*fact

         if(allocated(dummy)) deallocate(dummy)
      enddo
      end subroutine PH_INTEGRATE_SPECTRA
C
C
C read lines data from lines.dat file, provided from
c PHILIPS (lighting purposes)
C
      subroutine ph_read_philips
      implicit none
      integer :: numl,ios,i
      real(dp) :: g1,g2
      logical :: ex

      inquire (file='lines.dat',exist=ex)
      if(.not.ex) then
         write(6,*)'PHOTON MODULE (PH_READ_PHILIPS):'
         write(6,*)'   cant open lines.dat, exit.'
         call exit(1)
      else
         open(unit=77,file='lines.dat')
         read(77,FMT='(I6)', IOSTAT=ios) numl
         philips_numlines=numl
         allocate(philips_e0(numl))
         allocate(philips_e1(numl))
         allocate(philips_g1(numl))
         allocate(philips_e2(numl))
         allocate(philips_g2(numl))
         allocate(philips_aik(numl))
         allocate(philips_cs(numl))
         allocate(philips_cr(numl))
         allocate(philips_cw(numl))
         allocate(philips_ss(numl))
         allocate(philips_sr(numl))

         allocate(philips_l0(numl))
         allocate(philips_c3(numl))
         allocate(philips_c4(numl))
         allocate(philips_c6(numl))
         allocate(philips_name(numl))

         do i =1,numl
            if(HERMANN_FLAG) then
c HERMANN:
               read(77,FMT='(2E10.3,2E4.1,1E6.3,4E10.3)',IOSTAT=ios) 
     .           philips_l0(i),
     .           philips_e1(i),
     .           philips_g1(i),
     .           philips_g2(i),
     .           philips_aik(i),
     .           philips_cs(i),
     .           philips_cr(i),
     .           philips_cw(i),
     .           philips_ss(i)
            else
C MARGARITA:
               read(77,fmt='(8E12.4,A20)',iostat=ios)
     .           philips_l0(i),
     .           philips_aik(i),
     .           philips_g1(i),
     .           philips_g2(i),
     .           philips_e1(i),
     .           philips_c3(i),
     .           philips_c4(i),
     .           philips_c6(i),
     .           philips_name(i)
            endif
         enddo


c PHILIPS STUFF
c
c unit conversions:
         if(HERMANN_FLAG) then
C HERMANN:
            do i=1,numl
c  wavelength [m] --> Energy [eV]
               philips_e0(i) = hplnk*clight/(philips_l0(i)*100.)
c  energy levels [J] --> [eV]
               philips_e1(i) = philips_e1(i) / 1.6021892e-19
               philips_e2(i) = philips_e1(i) + philips_e0(i)
c  osc.strength   --> Aik Rate [1/s]
               g1=philips_g1(i)
               g2=philips_g2(i)
               philips_aik(i)= philips_aik(i)* 4.3392e7 *g1/g2
               philips_aik(i)= philips_aik(i)* philips_e0(i)**2
c  broadening constants Cs,Cr [m^4]
               philips_cs(i) = philips_cs(i)
               philips_cr(i) = philips_cr(i)
c  broadening constants Cw [m^7]
               philips_cw(i) = philips_cw(i)
c  line shift Ss [m^4]
               philips_ss(i) = philips_ss(i)
               philips_sr(i) = 0.
c  [m] --> [nm]
               philips_l0(i) = philips_l0(i)*1.e9
            enddo

            if (numl == 9) then
c remove inconsistency of upper levels for hermann's data
            philips_e2(4) = philips_e2(6)
            philips_e2(5) = philips_e2(6)

            philips_e0(4) = philips_e2(4)-philips_e1(4)
            philips_l0(4) = hplnk*clight*1e9/philips_e0(4)/100.
            philips_e0(5) = philips_e2(5)-philips_e1(5)
            philips_l0(5) = hplnk*clight*1e9/philips_e0(5)/100.

!pb special settings for checking 185nm, 254nm and 546nm lines
            philips_e2(9) = 6.7018_dp
            philips_e1(7) = philips_e2(9)
            philips_e1(8) = philips_e2(9)

            philips_e0(7) = philips_e2(7)-philips_e1(7)
            philips_l0(7) = hplnk*clight*1e9/philips_e0(7)/100.
            philips_e0(8) = philips_e2(8)-philips_e1(8)
            philips_l0(8) = hplnk*clight*1e9/philips_e0(8)/100.
            
            philips_e2(1) = 4.8812_dp
            philips_e0(1) = philips_e2(1)-philips_e1(1)
            philips_l0(1) = hplnk*clight*1e9/philips_e0(1)/100.

            end if

         else
C MARGARITA:
            do i=1,numl
c  wavelength l0 [nm] --> energy [eV]
               philips_e0(i) = hplnk*clight/(philips_l0(i)*1e-7)
c  upper energy level
               philips_e2(i) = philips_e1(i)+philips_e0(i)

!               write (40,*) ' i, philips_e0, philips_e1, philips_e2 ',
!     .           i, philips_e0(i), philips_e1(i), philips_e2(i)
            enddo
         endif

         close(77)
      endif
      return
      end subroutine ph_read_philips

      subroutine ph_free_philips
      implicit none
      if(allocated(philips_e0))  deallocate(philips_e0)
      if(allocated(philips_e1))  deallocate(philips_e1)
      if(allocated(philips_e2))  deallocate(philips_e2)
      if(allocated(philips_g1))  deallocate(philips_g1)
      if(allocated(philips_g2))  deallocate(philips_g2)
      if(allocated(philips_aik)) deallocate(philips_aik)
      if(allocated(philips_cs))  deallocate(philips_cs)
      if(allocated(philips_cr))  deallocate(philips_cr)
      if(allocated(philips_cw))  deallocate(philips_cw)
      if(allocated(philips_ss))  deallocate(philips_ss)
      if(allocated(philips_sr))  deallocate(philips_sr)
      if(allocated(philips_c3))  deallocate(philips_c3)
      if(allocated(philips_c4))  deallocate(philips_c4)
      if(allocated(philips_c6))  deallocate(philips_c6)
      if(allocated(philips_name))  deallocate(philips_name)
      if(allocated(philips_l0))  deallocate(philips_l0)
      return
      end subroutine ph_free_philips
c
c plasma profile    
c

      subroutine ph_philips_temp_dens(ipl)
      implicit none
      integer, intent(in) :: ipl
      real(dp) :: r,t,p,n
      integer :: nxm,nym,nzm,ir,ip,it,irad,iplti
      integer, save :: iipl=0
      
      if(iipl /= 0) return

c convert [bar] --> [Pa]
      p=phv_philips_pressure*1.e5
      if(p <= 0.) then
         write(6,*) 'PHOTON MODULE (PH_PHILIPS_TEMP_DENS):'
         write(6,*) '   PHV_PHILIPS_PRESSURE <= 0., exit.'
         call exit(1)
      endif

      iipl=ipl
      iplti=mplsti(ipl)
      NXM=MAX(1,NR1STM)
      NYM=MAX(1,NP2NDM)
      NZM=MAX(1,NT3RDM)
      do it=1,nzm
         do ip=1,nym
            do ir=1,nxm
               r=rhozne(ir)
               call ph_philips_tnprofile(r,p,t,n)
               irad = IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2
               TIIN(iplti,irad) = t
               DIIN(ipl,irad) = n
            enddo
         enddo
      enddo
      return
      end subroutine ph_philips_temp_dens

      subroutine ph_philips_tnprofile(rin, p, t, n)
      implicit none      
      real(dp), intent(in) :: rin,p
      real(dp), intent(out) :: t,n
      real(dp) :: r,res
c convert cm --> mm
      r=rin*10.
c res = temperature in [K}
      res=     7489.
      res=res- 10398. * r**2
      res=res+ 11106. * r**4
      res=res- 6126.  * r**6
      res=res+ 1790.  * r**8
      res=res- 260.6  * r**10
      res=res+ 14.24  * r**12
c n = density in [m**-3]..., p in [Pa]
      n = 7.12e22 * p/res
c ...convert to [cm**-3]
      n = n / 1e6

c convert temperature [K] --> [eV]
      t=res/1.1604e4

      return
      end subroutine ph_philips_tnprofile
      

      subroutine standard_van_der_waals
      use precision
      use ccona
      
      implicit none
      real(dp) :: aanf, aend, a, da, u, v, xint2, erg2,
     .            erg2old, erg3, erg3old, xint3, phimax, 
     .            phimin, dphi, phi, aanf_log, aend_log, da_log,
     .            xpos_phimax, xpos_log
      real(dp) :: dl12, dlqs, b, c, p1, p2, sqpii, p3
      real(dp), allocatable :: fp1(:), fp2(:), xip1(:), xip2(:), aar(:),
     .                         fp3(:), xip3(:)
      real(dp), allocatable :: phileft(:), phiright(:), xlleft(:), 
     .                         xlright(:), xintleft(:), xintright(:) 
      integer :: na, nphi, na2
      integer :: j, maxpos, ilook, iphi, 
     .           naleft, naright, namid, ipart, jj
      complex(dp) :: z, ii, zz, z2, z3h
      logical :: flag
      

!pb      call setcon 
      sqpii = 1._dp / sqrt(pia) 

      xpos_phimax = 1.243846358386715
      xpos_log = log(xpos_phimax)

      na = 10000
      na2 = na/2
      aend = 20000
     
      da_log = (log(aend) - xpos_log) / real(na2)

      allocate (aar(0:na))
!pb      allocate (fp1(0:na))
!pb      allocate (fp2(0:na))
      allocate (fp3(0:na))
!pb      allocate (xip1(0:na))
!pb      allocate (xip2(0:na))
      allocate (xip3(0:na))

      aar = 0._dp
!pb      fp1 = 0._dp
!pb      fp2 = 0._dp
      fp3 = 0._dp
!pb      xip1 = 0._dp
!pb      xip2 = 0._dp
      xip3 = 0._dp

      ii=cmplx(0._dp,1._dp)

      aar(na2) = xpos_phimax
      do j = 1, na2
        aar(na2+j) = exp(xpos_log + j*da_log)
      end do

      do j=1,na2
        aar(na2-j) = xpos_phimax - (aar(na2+j)-aar(na2))
      end do

      aanf = aar(0)

!pb      xint2 = 0._dp
      xint3 = 0._dp

      a = aanf
      z = -cmplx(a,1._dp) / (1._dp + a*a)
      z3h = z * cdsqrt(z)
      zz = cdsqrt(z)*ii
      
      call ph_wofz(real(zz,dp), aimag(zz), u, v, flag)
      z2 = cmplx(u,v) * z3h
      erg2old = aimag(z2)

      erg3old = -sqpii*aimag(z) + erg2old

!pb      fp1(0) = 1._dp / (1._dp + a*a)
!pb      fp2(0) = erg2old
      fp3(0) = erg3old
!pb      xip1(0) = 0._dp
!pb      xip2(0) = 0._dp
      xip3(0) = 0._dp
      
      j = 0

!pb      write (6,'(a)') ' j, a, (u v), z3h, erg1, xint1 '
      do j = 1, na

         a = aar(j) 
         da = aar(j) - aar(j-1)
         
         z = -cmplx(a,1._dp)  / (1._dp + a*a)
         z3h = z * cdsqrt(z)
         zz = cdsqrt(z)*ii
      
         call ph_wofz(real(zz,dp), aimag(zz), u, v, flag)
         z2 = cmplx(u,v) * z3h
         erg2 = aimag(z2)

         erg3 = -sqpii*aimag(z) + erg2
         
!pb         xint2 = xint2 + 0.5_dp*(erg2+erg2old)*da
         xint3 = xint3 + 0.5_dp*(erg3+erg3old)*da

!pb         fp1(j) = 1._dp / (1._dp + a*a)
!pb         fp2(j) = erg2
         fp3(j) = erg3
!pb         xip1(j) = xip1(j-1) + 0.5_dp*(fp1(j)+fp1(j-1))*da
!pb         xip2(j) = xint2
         xip3(j) = xint3

         erg2old = erg2
         erg3old = erg3

!pb         write (16,'(i6,4es12.4)') j, aar(j), fp1(j), fp2(j), fp3(j)

      end do

      dl12 = 1._dp
      dlqs = 4._dp / pia

      b = pia * dlqs / (4._dp * dl12)    ! b=1
      c = sqrt(dlqs) / (2._dp * pia * dl12**1.5_dp)

!pb      p1 = xip1(na) / (pia * dl12)

!pb      p2 = c * pia / b**2.5_dp * xip2(na)

      p3 = 1._dp /(pia*dl12) * sqrt(pia) / b**2 * xip3(na)

!pb      write (6,*) p1, p2, p1+p2, p3
      write (6,*) p3

!pb      maxpos = maxloc(fp3,dim=1)-1
!pb      write (6,*) ' position of maximum ',maxpos
!pb      write (6,*) ' a,f(a) = ',aar(maxpos), fp3(maxpos)

      nphi = 1001
      allocate (phileft(nphi))
      allocate (phiright(nphi))
      allocate (xlleft(nphi))
      allocate (xlright(nphi))
      allocate (xintleft(nphi))
      allocate (xintright(nphi))
      phileft = 0._dp
      phiright = 0._dp
      xlleft = 0._dp
      xlright = 0._dp
      xintleft = 0._dp
      xintright = 0._dp
      
!pb      phimax = log(fp3(maxpos))
      phimax = log(fp3(na2))
      phimin = log(fp3(0))
      dphi = (phimax - phimin) / real(nphi-1)
      
! linke Seite
      
      phileft(1) = phimax
!pb      xlleft(1) = aar(maxpos)
      xlleft(1) = aar(na2)
!pb      ilook = maxpos
      ilook = na2

      do iphi = 2,nphi

        phi = phimax - (iphi-1)*dphi
        do while ((ilook > 1) .and. (log(fp3(ilook-1)) > phi))
          ilook = ilook - 1
        end do

        phileft(iphi) = phi
        xlleft(iphi) = aar(ilook) - (log(fp3(ilook))-phi) / 
     .                 (log(fp3(ilook))-log(fp3(ilook-1))) *
     .                 (aar(ilook)-aar(ilook-1)) 
        xintleft(iphi) = xip3(ilook) - (log(fp3(ilook))-phi) / 
     .                 (log(fp3(ilook))-log(fp3(ilook-1))) *
     .                 (xip3(ilook)-xip3(ilook-1)) 
      end do
      
! rechte Seite
      
      phiright(1) = phimax
!pb      xlright(1) = aar(maxpos)
      xlright(1) = aar(na2)
!pb      ilook = maxpos
      ilook = na2

      do iphi = 2,nphi

        phi = phimax - (iphi-1)*dphi
        do while ((ilook > na) .and. (log(fp3(ilook+1)) > phi))
          ilook = ilook + 1
        end do

        phiright(iphi) = phi
        xlright(iphi) = aar(ilook) + (log(fp3(ilook))-phi) / 
     .                 (log(fp3(ilook+1))-log(fp3(ilook))) *
     .                 (aar(ilook+1)-aar(ilook)) 
        xintright(iphi) = xip3(ilook) + (log(fp3(ilook))-phi) / 
     .                 (log(fp3(ilook+1))-log(fp3(ilook))) *
     .                 (xip3(ilook+1)-xip3(ilook)) 
      end do

      deallocate (aar)
!pb      deallocate (fp1)
!pb      deallocate (fp2)
      deallocate (fp3)
!pb      deallocate (xip1)
!pb      deallocate (xip2)
      deallocate (xip3)

      return

      end subroutine standard_van_der_waals


      END MODULE PHOTON
