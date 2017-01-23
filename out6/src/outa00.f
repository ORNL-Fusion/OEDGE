c     -*-Fortran-*-
c
c
      subroutine loadm_axis(mouts,mwids,ir,ip,axistype,offset)
      implicit none
      include 'params'
      include 'cgeom'
c
      integer ir,ip,axistype,offset
      real mouts(maxdatx,maxplts,maxngs)
      real mwids(maxdatx,maxplts,maxngs)
c
c     This routine calculates the X-axis for m-series plots by
c     proceeding along the ring and copying the appropriate data
c     into the mouts and mwids arrays.
c
c     AxisType = 1 for S
c              = 2 for P
c     Offset   = 0 - do NOT load target values
c     Offset   = 1 - load target values at both ends of the ring
c
c     Local variables
c
      integer in,ik
c
c     Set in equal to zero for the default case - things are only
c     adjusted at the moment for offset values equal to 1.
c
      in = 0
c
      if (offset.eq.1) then
c
         in = 1
c
         mouts(1,ip,1) = 0.0
c
         if (axistype.eq.1) then
            mWIDS(1,ip,1) = kss(1,ir)
            mouts(nks(ir)+2,ip,1) = ksmaxs(ir)
            mwids(nks(ir)+2,ip,1) = ksmaxs(ir) - kss(nks(ir),ir)
         else
            mWIDS(1,ip,1) = kps(1,ir)
            mouts(nks(ir)+2,ip,1) = kpmaxs(ir)
            mwids(nks(ir)+2,ip,1) = kpmaxs(ir) - kps(nks(ir),ir)
         endif
c
      endif
c
      DO IK = 1, NKS(IR)

         IF (axistype.eq.1) THEN
            mOUTS(IK+in,ip,1) = KSS(IK,IR)
            mWIDS(IK+in,ip,1) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
         ELSE
            mOUTS(IK+in,ip,1) = KPS(IK,IR)
            mWIDS(IK+in,ip,1) = 0.0
            IF (IK.GT.1) mWIDS(IK+in,ip,1)=
     >                              0.5*(KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) mWIDS(IK+in,ip,1) =
     >                              mWIDS(IK+in,ip,1) +
     >                              0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
         ENDIF
c
      enddo
C
      return
      end
c
c
c
      real function osm_interpolate(osmval1,osmval2,osmaxis,osmvals,
     >                         exp_coord,err_default,ratio_opt)
      implicit none
      integer osmvals,ratio_opt
      real osmval1(osmvals),osmval2(osmvals),osmaxis(osmvals)
      real exp_coord
      real err_default
c
c     OSM_INTERPOLATE:
c
c     This routine takes two arrays of data and an ordinate axis
c     for these data as input. It will try to find the values
c     in the arrays which are adjacent to the input coordinate. It will
c     then interpolate between these values to obtain the value
c     required in the calling routine. If the value of exp_coord lies
c     beyond the end of the array then the routine will return the value of
c     the closest value entry depending on the value of ratio_opt.
c
c     If ratio_opt = 1 - the code returns the value in the osmval1 array
c                        closest to exp_coord and will interpolate between
c                        values for points within the data.
c
c     If ratio_opt = 2 - the code returns the interpolated value of
c                        the quantity osmval1/osmval2.
c
c     If ratio_opt = 3 - this retuens the interpolated value of osmval1 for
c                        data inside the range and the ratio of the
c                        values osmval1/osmval2 * err_default of the
c                        closest data for points beyond the range.
c
c
c     Local Variables
c
      real fact,val1,val2
      integer in,ipos,loc
      external ipos
c
c     Set default return value - it is always overwritten
c
      osm_interpolate = 0.0
c
c     Data point less than available axis
c
      if (exp_coord.lt.osmaxis(1))then
c
c        Extract value
c
         if (ratio_opt.eq.1) then
c
            osm_interpolate = osmval1(1)
c
c
c        Extract Ratio
c
         elseif (ratio_opt.eq.2) then
c
c           Verify that denominator is non-zero
c
            if (osmval2(1).ne.0.0) then

               osm_interpolate = osmval1(1)/osmval2(1)

            else

               write (6,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(1),osmval2(1)
c
               write (0,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(1),osmval2(1)
c
               osm_interpolate = err_default
c
            endif
c
         elseif (ratio_opt.eq.3) then
c
c           Verify that denominator is non-zero
c
            if (osmval2(1).ne.0.0) then

               osm_interpolate = osmval1(1)/osmval2(1) * err_default

            else

               write (6,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(1),osmval2(1)
c
               write (0,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(1),osmval2(1)
c
               osm_interpolate = err_default
c
            endif
c
         endif
c
c     Data point greater than available axis
c
      elseif (exp_coord.gt.osmaxis(osmvals)) then
c
c        Extract value
c
         if (ratio_opt.eq.1) then
c
            osm_interpolate = osmval1(osmvals)
c
c
c        Extract Ratio
c
         elseif (ratio_opt.eq.2) then
c
c           Verify that denominator is non-zero
c
            if (osmval2(osmvals).ne.0.0) then

               osm_interpolate = osmval1(osmvals)/osmval2(osmvals)

            else

               write (6,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(osmvals),osmval2(osmvals)

               write (0,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(osmvals),osmval2(osmvals)
c
               osm_interpolate = err_default
c
            endif
c
         elseif (ratio_opt.eq.3) then
c
c           Verify that denominator is non-zero
c
            if (osmval2(osmvals).ne.0.0) then

               osm_interpolate = osmval1(osmvals)/osmval2(osmvals)
     >                           * err_default
c
            else

               write (6,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(osmvals),osmval2(osmvals)

               write (0,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(osmvals),osmval2(osmvals)
c
               osm_interpolate = err_default
c
            endif
c
         endif



c
c     Data point within available axis
c
      else
c
c        Find cell first above given axis coordinate
c
         loc = ipos(exp_coord,osmaxis,osmvals)
c
         fact = (exp_coord-osmaxis(loc-1))/(osmaxis(loc)-osmaxis(loc-1))
c
c        Extract value
c
         if (ratio_opt.eq.1.or.ratio_opt.eq.3) then
c
            osm_interpolate = osmval1(loc-1)
     >                        + fact * (osmval1(loc)-osmval1(loc-1))
c
c        Extract Ratio
c
         elseif (ratio_opt.eq.2) then
c
            val1 = osmval1(loc-1)
     >                        + fact * (osmval1(loc)-osmval1(loc-1))
            val2 = osmval2(loc-1)
     >                        + fact * (osmval2(loc)-osmval2(loc-1))
c
            if (val2.ne.0.0) then

               osm_interpolate = val1/val2

            else

               write (6,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(osmvals),osmval2(osmvals)

               write (0,*) 'ERROR: Denominator is 0.0'//
     >                     ' in OSM_INTERPOLATE',
     >                      osmval1(osmvals),osmval2(osmvals)
c
               osm_interpolate = err_default
c
            endif
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
      subroutine load_rzdata(iseld,ndata,rzdata,maxnpts,max_incols,
     >                       axis_offset_r,axis_offset_z,datatitle)
      implicit none
c
      include 'params' 
c
      integer iseld,ndata,maxnpts,max_incols
      real rzdata(maxnpts,max_incols)
      real axis_offset_r,axis_offset_z
      character*(*) datatitle
c
c     load_rzdata :  This routine loads the results of a
c                    thompson scattering analysis (or other R,Z
c                    based diagnostic) assuming a
c                    set of (R,Z, Value) entries in a table.
c                    Where Value may be the density or
c                    temperature or another quantity.
c
c
c     Local variables
c
      integer dataunit,maxcols
      parameter (dataunit=13,maxcols=8)
c
c      integer maxdatx,dataunit,maxcols
c      parameter (maxdatx=1000,dataunit=13,maxcols=8)
c
      integer in,ik,axis_type
      integer num_expt,ncols
      real expt_axis(maxdatx)
      real expt_data(maxdatx,maxcols)
      real cs,osm_interpolate
      real jsato,jsati,teav
      external osm_interpolate
c
c     Load experimental data - ncols is expected to be 3 or more
c     3 = R + Z + value (+ value) ...
c
      write (6,*) 'Loading:',iseld
c
      call load_expt_data(dataunit,iseld,expt_axis,axis_type,
     >                    expt_data,maxcols,maxdatx,
     >                    num_expt,ncols,datatitle)
c
      write (6,*) 'Loaded:',num_expt,ncols,datatitle
c
      if ((ncols+2).gt.max_incols) then
c
c        Limit number of columns of data to maximum allowed in arrays
c        and issue error message - note - ncols MUST be increased by
c        one because the experimental data loading routine stores the
c        R coordinate as an axis coordinate and not a data value -
c        this is corrected by this loading routine.
c
         ncols = max_incols -2
c
         write(6,*) 'ERROR: LOADING RZ-DATA - NUMBER OF COLUMNS'//
     >              ' EXCEEDS AVAILABLE STORAGE',ncols+1,max_incols
      endif
c
c     Transfer data to rz data arrays
c
      if (num_expt.gt.maxnpts) then

         write(6,*) 'ERROR: LOADING RZ-DATA - NUMBER OF DATA'//
     >              ' ENTRIES EXCEEDS AVAILABLE STORAGE',
     >                      num_expt,maxnpts
         num_expt = maxnpts
c
      endif
c
      ndata = num_expt
c
      do in = 1, num_expt

         rzdata(in,1) = expt_axis(in) + axis_offset_r
         rzdata(in,2) = expt_data(in,1) + axis_offset_z

         do ik = 2,ncols

            rzdata(in,2+ik) = expt_data(in,ik)

         end do

      end do
c
c     Calculate an axis estimate and store in element 3
c
      rzdata(1,3) = 0.0
c
      do in = 2,ndata
c
         rzdata(in,3) = rzdata(in-1,3) +
     >              sqrt((rzdata(in,2)-rzdata(in-1,2))**2+
     >                   (rzdata(in,1)-rzdata(in-1,1))**2)
c
      end do
c
      write (6,*)  'RZDATA:',ndata,ncols
c
      do in = 1,ndata
         write(6,'(5(1x,g12.5))') (rzdata(in,ik),ik=1,ncols+2)
      end do
c
      return
      end
c
c
c 
      subroutine writedata
      implicit none
c
      include 'params'
      include 'outcom'
c
      include 'cgeom'
      include 'comtor'
c      include 'cneut2'
      include 'dynam2'
      include 'dynam3'
c      include 'dynam4'
      include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
      include 'printopt' 
c
c     Local variables
c
c
      integer prnizs
      real    bgcontent,totbgcontent
      real    impcontent(0:maxizs+1),totimpcontent(0:maxizs+1)
      REAL    ZSUM(max(MAXPLRP,maxizs))
      real tmpsum,tmpsum2
      real tote,toti,totn
      real r,z
      integer ik,ir,iz
      integer id,in 

c
c     Temporary variables
c
      real temppr1(maxnks,maxnrs),temppr2(maxnks,maxnrs)


c     Print out time dependent data for charge state 1
c
c
c      do it = 1,nts
c         CALL PRRMATDIV(lims(1,1,1,it),MAXNKS,nks(irsep),
c     >                  NRS,6,'Z=1 TIME DEP')
c      end do
c
C
C     WRITE OUT X,Y GRID DIAGRAM
C
C      DO 654 IY = 1,NYS
C        WRITE (6,'(500I1)') (IFXYS(IX,IY),IX=1,NXS/3)
C 654  CONTINUE
C      DO 655 IY = 1,NYS
C        WRITE (6,'(500I1)') (IFXYS(IX,IY),IX=NXS/3+1,(2*NXS)/3)
C 655  CONTINUE
C      DO 656 IY = 1,NYS
C        WRITE (6,'(500I1)') (IFXYS(IX,IY),IX=((2*NXS)/3)+1,NXS)
C 656  CONTINUE
C
c
c
c     Print out EDGE2D vs. DIVIMP Arrays
c
c
      if (cre2d.eq.1.and.cre2dizs.gt.-1) then
c
         write (6,*) 'EDGE2D - States Read - ' , cre2dizs
c
         do iz = 0,cre2dizs

            write (6,*) 'Ionization STATE IZ = ',iz

            do ir = 1,nrs
               do ik = 1,nks(ir)
                  temppr1(ik,ir) = sdlims(ik,ir,iz)*absfac
                  temppr2(ik,ir) = e2dnzs(ik,ir,iz)
               end do
            end do

         CALL PRRMATDIV(temppr1,MAXNKS,nks(irsep),NRS,6,'DIVIMP-Zn')
         CALL PRRMATDIV(temppr2,MAXNKS,nks(irsep),NRS,6,'EDGE2D-Zn')


            do ir = 1,nrs
               do ik = 1,nks(ir)
                  temppr1(ik,ir) = powls(ik,ir,iz)*absfac
                  temppr2(ik,ir) = e2dpowls(ik,ir,iz)
               end do
            end do

         CALL PRRMATDIV(temppr1,MAXNKS,nks(irsep),NRS,6,'DIVIMP-POW')
         CALL PRRMATDIV(temppr2,MAXNKS,nks(irsep),NRS,6,'EDGE2D-POW')

            do ir = 1,nrs
               do ik = 1,nks(ir)
                  temppr1(ik,ir)=powls(ik,ir,iz)*absfac*kareas(ik,ir)
                  temppr2(ik,ir)=e2dpowls(ik,ir,iz)*kareas(ik,ir)
               end do
            end do

         CALL PRRMATDIV(temppr1,MAXNKS,nks(irsep),NRS,6,
     >                                         'DIVIMP-POW-ABS')
         CALL PRRMATDIV(temppr2,MAXNKS,nks(irsep),NRS,6,
     >                                         'EDGE2D-POW-ABS')

         end do


         write (6,*) 'Sum over all ionization states'//
     >                   ' (except neutrals)'

c
c        Sum over density
c
         call rzero (temppr1,maxnks*maxnrs)
         call rzero (temppr2,maxnks*maxnrs)

         do iz = 1,cre2dizs
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  temppr1(ik,ir) = temppr1(ik,ir) +
     >                         sdlims(ik,ir,iz)*absfac*kareas(ik,ir)
                  temppr2(ik,ir) = temppr2(ik,ir) +
     >                         e2dnzs(ik,ir,iz)*kareas(ik,ir)
               end do
            end do
         end do
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (kareas(ik,ir).gt.0.0) then
                  temppr1(ik,ir) = temppr1(ik,ir) / kareas(ik,ir)
                  temppr2(ik,ir) = temppr2(ik,ir) / kareas(ik,ir)
               end if
            end do
         end do
c
         CALL PRRMATDIV(temppr1,MAXNKS,nks(irsep),NRS,6,'DIVIMP-ZT')
         CALL PRRMATDIV(temppr2,MAXNKS,nks(irsep),NRS,6,'EDGE2D-ZT')

c
c        Sum over POWER
c

         call rzero (temppr1,maxnks*maxnrs)
         call rzero (temppr2,maxnks*maxnrs)

         do iz = 1,cre2dizs
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  temppr1(ik,ir) = temppr1(ik,ir) +
     >                         powls(ik,ir,iz)*absfac*kareas(ik,ir)
                  temppr2(ik,ir) = temppr2(ik,ir) +
     >                         e2dpowls(ik,ir,iz)*kareas(ik,ir)
               end do
            end do
         end do
c
c        Absolute Power
c
         CALL PRRMATDIV(temppr1,MAXNKS,nks(irsep),NRS,6,
     >                               'DIVIMP-POW-ABS-T')
         CALL PRRMATDIV(temppr2,MAXNKS,nks(irsep),NRS,6,
     >                               'EDGE2D-POW-ABS-T')
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (kareas(ik,ir).gt.0.0) then
                  temppr1(ik,ir) = temppr1(ik,ir) / kareas(ik,ir)
                  temppr2(ik,ir) = temppr2(ik,ir) / kareas(ik,ir)
               end if
            end do
         end do
c
c        Power density
c
         CALL PRRMATDIV(temppr1,MAXNKS,nks(irsep),NRS,6,
     >                                  'DIVIMP-POW-DEN-T')
         CALL PRRMATDIV(temppr2,MAXNKS,nks(irsep),NRS,6,
     >                                  'EDGE2D-POW-DEN-T')
c

c
c        KAREAS
c
         CALL PRRMATDIV(kareas,MAXNKS,nks(irsep),NRS,6,
     >                                  'KAREAS')


      endif
c
c     Calculate the core content
c
      core_content=0.0
      core_area=0.0
c
      do ir = 1,irsep-1
         do ik = 1,nks(ir)-1
            do iz = 1,nizs
               core_content = core_content + sdlims(ik,ir,iz)
     >                                *kareas(ik,ir)
            end do
            core_area = core_area+ kareas(ik,ir)
         end do
      end do
c
      do ir = irsep,nrs
         do ik = 1,nks(ir)
            z = zs(ik,ir)
c
            do iz = 1,nizs

               edge_content = edge_content +   sdlims(ik,ir,iz)
     >                                *kareas(ik,ir)
c
c               write (6,*) ik,ir,iz,edge_content
c
            end do
            edge_area = edge_area+ kareas(ik,ir)
c

            if (ir.ge.irtrap.and.ir.le.nrs) then

               do iz = 1,nizs
                  pp_content = pp_content +   sdlims(ik,ir,iz)
     >                             *kareas(ik,ir)
               end do
               pp_area = pp_area+ kareas(ik,ir)

            elseif (z.le.zxp) then

               do iz = 1,nizs
                  div_content = div_content + sdlims(ik,ir,iz)
     >                             *kareas(ik,ir)
               end do
               div_area = div_area+ kareas(ik,ir)

            else

               do iz = 1,nizs
                  main_content = main_content + sdlims(ik,ir,iz)
     >                             *kareas(ik,ir)
           end do
               main_area = main_area+ kareas(ik,ir)
c
            endif
         end do
      end do
c
      write (6,*) 'Regions and Impurity Content:'
      write (6,*) 'Region       Area           Content'
      write (6,'(a8,2g16.8)') 'Core :', core_area,core_content
      write (6,'(a8,2g16.8)') 'Edge :', edge_area,edge_content
      write (6,'(a8,2g16.8)') 'Trap :', pp_area,pp_content
      write (6,'(a8,2g16.8)') 'Div  :', div_area,div_content
      write (6,'(a8,2g16.8)') 'Main :', main_area,main_content
      write (6,'(a8,2g16.8)') 'Total:', core_area+edge_area,
     >                       core_content+edge_content
c
c     Calculate content and C/D ratios near the inner target
c
      totbgcontent=0.0
      call rzero(totimpcontent,maxizs+2)
c
      do ir = irsep,irwall-1
c
         bgcontent=0.0
         call rzero(impcontent,maxizs+2)
c
c slmod begin
c...     Ha!  A ring with 8 or fewer cells is problems:
         do ik = MAX(1,nks(ir)-8),nks(ir)
c
c         do ik = nks(ir)-8,nks(ir)
c slmod begin
c
            bgcontent = bgcontent + knbs(ik,ir)*kareas(ik,ir)
c
            do iz = 0,nizs
c
               impcontent(iz) = impcontent(iz) + sdlims(ik,ir,iz)*
     >                          kareas(ik,ir) * absfac
c
            end do
c
         end do
c
c        Sum over charge states
c
         do iz = 0,nizs
c
               impcontent(nizs+1) = impcontent(nizs+1)
     >                            + impcontent(iz)
c
         end do
c
c        Print out results for ring
c
         write (6,*) 'CONTENT:', ir, ' HYDROGEN = ',bgcontent

         do iz = 0,nizs

            write(6,*) 'C ',iz,' - ',impcontent(iz),' RATIO = ',
     >                 impcontent(iz)/bgcontent

         end do
c
         write(6,*) 'C TOT - ',impcontent(nizs+1),' RATIO = ',
     >                 impcontent(nizs+1)/bgcontent
c
c        Add to grand totals
c
         totbgcontent = totbgcontent + bgcontent
c
         do iz = 0,nizs+1
            totimpcontent(iz) = totimpcontent(iz) + impcontent(iz)
         end do
c
      end do
c
c     Print out grand totals
c
      write (6,*) 'TOTAL CONTENT:  HYDROGEN = ',totbgcontent

      do iz = 0,nizs

         write(6,*) 'C ',iz,' - ',totimpcontent(iz),' RATIO = ',
     >              totimpcontent(iz)/totbgcontent

      end do

      write(6,*) 'C TOT - ',totimpcontent(nizs+1),' RATIO = ',
     >              totimpcontent(nizs+1)/totbgcontent
c
c 
c      Print out location specific impurity density profiles
c
      call pr_imp_density_profiles
c
c     print out impurity deposition
c     
      call print_deposition(cgridopt)

c
c     Print out background plasma and polygon information
c
      write (6,
     >  '(''# CORE PLASMA RINGS   : IR = '',i4,'' TO '',i4)') 
     >                      1,irsep-1 
      write (6,
     >  '(''# MAIN SOL RINGS      : IR = '',i4,'' TO '',i4)') 
     >                     irsep,irwall-1
      write (6,
     >  '(''# PRIVATE PLASMA RINGS: IR = '',i4,'' TO '',i4)') 
     >                     irtrap+1,nrs
      write (6,'(''# NOTE: RING IR = '','//
     >'i4,'' IS THE PP RING NEXT TO SEPARATRIX IR = '',i4)')
     > nrs,irsep
      write(6,'(''# RINGS WITHOUT POLYGONS ARE NOT LISTED'')')
      write(6,'(''# UNITS:'')')
      write(6,'(''# RC,ZC ...        M'')')
      write(6,'(''# TE,TI            eV'')')
      write(6,'(''# NB               M^-3'')')
      write(6,'(''# VB               MS^-1'')')
      write(6,'(''# E                V'')')
c
      write(6,'(''# NUMBER OF RINGS:'',i4)') nrs
      write(6,'(''# NUMBER OF CELLS ON EACH RING'//
     >             ' (INCLUDING RINGS W/O POLYGONS):'')')
      write(6,'(100i5)') (nks(ir),ir=1,nrs)
c
c     Target plasma conditions
c
      write (6,*)  
      write (6,'(''# TARGET PLASMA CONDITIONS and LOCATIONS'')')
      write (6,*)  
      write (6,
     >     '(''# ID'',2x,''IR'',2x,''IK'',5x,''RP'',11x,''ZP'','//
     >         '11x,''TE'',11x,''TI'',11x,''NE'',11x,''VB'','//
     >         '12x,''E'')')
c
      do id = 1,nds
c
         ik = ikds(id)
         ir = irds(id)
c
         in = korpg(ik,ir)
c            
         if (in.gt.0) then
c
               write(6,'(3i4,1p,16(1x,g12.5))') 
     >           id,ir,ik,
     >           rp(id),zp(id),kteds(id),
     >           ktids(id),knds(id),kvds(id),
     >           keds(id)

         endif  
c
      end do      

      write (6,*) 
      write (6,'(''# BACKGROUND PLASMA and POLYGON INFORMATION'')')
      write (6,*) 
c
      write (6,'(''# IR'',2x,''IK'',5x,''RC'',11x,''ZC'','//
     >         '11x,''TE'',11x,''TI'',11x,''NE'',11x,''VB'','//
     >         '12x,''E'',11x,''RV1'',10x,''ZV1'',10x,''RV2'','//
     >         '10x,''ZV2'',10x,''RV3'',10x,''ZV3'',10x,''RV4'','//
     >         '10x,''ZV4'')')
c
      do ir = 1,nrs
c
         do ik = 1,nks(ir)
c
            in = korpg(ik,ir)
            if (in.gt.0) then 
c
               write(6,'(2i4,1p,16(1x,g12.5))') 
     >           ir,ik,
     >           rs(ik,ir),zs(ik,ir),ktebs(ik,ir),
     >           ktibs(ik,ir),knbs(ik,ir),kvhs(ik,ir)/qtim,
     >           kes(ik,ir)/(QTIM*QTIM*EMI/CRMI),
     >          (rvertp(id,in),zvertp(id,in),id = 1,nvertp(in))
c
            endif 

         end do

      end do 

c
c     Assign value for number of printable charge states
c
      prnizs = min(6,nizs)


C
C---- PRINT TABLE OF ABSOLUTE POWER VALUES
C
      WRITE(6,'(//1X,''TABLE OF POWER VALUES:'')')
      WRITE (6,'(/1X,''ABSOLUTE POWER FACTOR ='',1P,E9.2,//)') ABSFAC
      DO 80 IR = 1, NRS
        WRITE (6,9031) (ZLABS(IZ)(5:11),IZ=-1,PRNIZS),
     >                 'SUM RAD','SUM DEN'
        WRITE (6,9032)
        DO 80 IK = 1, NKS(IR)
c
          tmpsum = 0.0
          tmpsum2 = 0.0
          do iz = 1,nizs
             tmpsum = tmpsum + powls(ik,ir,iz)
             tmpsum2 = tmpsum2 + sdlims(ik,ir,iz)
          end do
c
c 9033 FORMAT(1X,2I3,2F7.3,1P,E8.2,10E9.2)
c
          WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),KAREAS(IK,IR),
     >        (ABSFAC*POWLS(IK,IR,IZ),IZ=-1,prnizs),absfac*tmpsum,
     >        absfac*tmpsum2
   80 CONTINUE
c
c     Analyse the Radiation source. Determine mean ne and niz
c     from the radiating volume. (Assumed to be highest 2/3 of
c     radiating cells.
c
      call radproc(nizs,job,pradclev)
c
C
C---- PRINT TABLE OF DENSITY VALUES
C
      WRITE(6,'(//1X,''TABLE OF DENSITY VALUES:'')')
      DO 81 IR = 1, NRS
        CALL RZERO (ZSUM, NIZS)
        WRITE (6,9031) (ZLABS(IZ)(5:11),IZ=-1,PRNIZS),'SUM DEN'
        WRITE (6,9032)
        DO 82 IK = 1, NKS(IR)
c
          tmpsum = 0.0
          do iz = 1,nizs
             tmpsum = tmpsum + sdlims(ik,ir,iz)
          end do
c
          WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),KAREAS(IK,IR),
     >            (absfac*SDLIMS(IK,IR,IZ),IZ=-1,PRNIZS),absfac*tmpsum
          DO 83 IZ = 1 , NIZS
             ZSUM(IZ) = ZSUM(IZ) + ABSFAC*SDLIMS(IK,IR,IZ)*KAREAS(IK,IR)
   83     CONTINUE
   82   CONTINUE
        WRITE (6,9034) ( ZSUM(IZ),IZ=1,PRNIZS)
   81 CONTINUE
C
C---- PRINT TABLE OF IONIZATION VALUES
C
      WRITE(6,'(//1X,''TABLE OF IONIZATION VALUES:'')')
      DO 84 IR = 1, NRS
        CALL RZERO (ZSUM, NIZS)
        WRITE (6,9031) (ZLABS(IZ)(5:11),IZ=-1,PRNIZS),
     >                 'PIN HIZ','E2D HIZ'
        WRITE (6,9032)
        DO 85 IK = 1, NKS(IR)
          WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),KAREAS(IK,IR),
     >           (absfac*TIZS(IK,IR,IZ),IZ=-1,PRNIZS),pinion(ik,ir),
     >           e2dion(ik,ir)
          DO 86 IZ = 1 , NIZS
             ZSUM(IZ) =ZSUM(IZ)+absfac*TIZS(IK,IR,IZ)*KAREAS(IK,IR)
   86     CONTINUE
   85   CONTINUE
        WRITE (6,9034) ( ZSUM(IZ),IZ=1,PRNIZS)
   84 CONTINUE
C
C---- PRINT TABLE OF PLRPS VALUES
C
      WRITE(6,'(A)') 'TABLE OF PLRP VALUES:'
      DO 91 IR = 1, NRS
        CALL RZERO (ZSUM, PLRPCNT)
        WRITE (6,9031) (PLABS(IZ)(5:11),IZ=-1,PLRPCNT)
        WRITE (6,9032)
        DO 92 IK = 1, NKS(IR)
          WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),KAREAS(IK,IR),
     >                          (ABSFAC*PLRPS(IK,IR,IZ),IZ=-1,PLRPCNT)
          DO 93 IZ = 1 , PLRPCNT
             ZSUM(IZ) = ZSUM(IZ) + ABSFAC*PLRPS(IK,IR,IZ)*KAREAS(IK,IR)
   93     CONTINUE
   92   CONTINUE
        WRITE (6,9034) ( ZSUM(IZ),IZ=1,PLRPCNT)
   91 CONTINUE
c
c     Print wall erosion and deposition data
c
      write(6,*) 'Deposition and Erosion'

      tote = 0.0
      toti = 0.0
      totn = 0.0

      do in = 1,wallpts
c
c slmod begin - new
c ARRAY BOUNDS: WALLPT(,18) is 0 for non-target wall segments, which is out
c               of bounds for dds2.
        IF (wallpt(in,18).EQ.0.0) CYCLE
c slmod end
         write (6,'(a,i5,f7.2,3(1x,g12.5),4(1x,f10.3))')
     >       'Erosion:',in,wallpt(in,18),wallpt(in,1),
     >       wallpt(in,2),dds2(INT(wallpt(in,18))),
     >       wallsn(in),wallsi(in),wallse(in),
     >       wallsn(in)+wallsi(in)-wallse(in)
         tote = tote + wallse(in)
         toti = toti + wallsi(in)
         totn = totn + wallsn(in)


      end do
      write (6,'(a,4(1x,f10.3))')
     >       'Erosion:',
     >       wallsn(maxpts+1),wallsi(maxpts+1),wallse(maxpts+1),
     >       wallsn(maxpts+1)+wallsi(maxpts+1)-wallse(maxpts+1)
      write (6,'(a,4(1x,f10.3))')
     >       'Erosion:',
     >         totn,toti,tote,totn+toti-tote

c
c     Print out wall erosion and deposition data. 
c
      call pr_calc_walldep
c
c     Print eirene data analysis
c
      call pr_eirene_analysis
c
c     Print exb analysis
c
      call pr_exb_analysis
c


      return 

c
c     Format statements 
c



 9031 FORMAT(/1X,' IK IR    R      Z     AREA',10(2X,A7))
 9032 FORMAT(1X,131('-'))
 9033 FORMAT(1X,2I3,2F7.3,1P,E8.2,10E9.2)
 9034 FORMAT(29X , 1P , 10E9.2 )


      end
c
c
c

      subroutine init_plot(iref,graph,iopt)
      implicit none
      character*(*) graph
      integer iopt,iref
c
      include 'params'
      include 'outcom'
      include 'cgeom' 
      include 'comtor'      
c
c     Local Variables
c
      real aspect


c
c     write (6,*) 'Plot:',iref,iopt
c
      WRITE (iplot,9001) GRAPH(4:58),IOPT,GRAPH(1:3)
      NPLOTS = NPLOTS + 1
C
C     INITIALIZE SMOOTHING VARIABLES - EACH ITERATION
C
      ITEC = 1
      ISMOTH = 99
c
c     Initialize descriptive strings in case a plot has set them
c
      NVIEW  = ' '
      PLANE  = ' '
      ANLY   = ' '
      TABLE  = 'SYMBOL TABLE'
      SMOOTH = ' '
c
c     Set plotting limits
c
C
c slmod begin
      IF (2*(IREF/2).EQ.IREF.OR.IREF.EQ. 91.OR.IREF.EQ.93.OR.
     .                          IREF.EQ.981.OR.IREF.EQ.983.OR.
     .                          IREF.EQ.985.OR.IREF.EQ.987.OR.
     .                          IREF.EQ.989) THEN
c
c      IF (2*(IREF/2).EQ.IREF.OR.IREF.EQ.91.OR.IREF.EQ.93) THEN
c slmod end
        XPOINT= ' NEAR X PT'
        if (cgridopt.eq.2) then
           XXMIN = RXP - 1.0
           XXMAX = RXP + 1.5
           IF (ZXP-ZMIN.LT.1.25) THEN
             YYMIN = ZMIN
             YYMAX = YYMIN + 2.5
           ELSE IF (ZMAX-ZXP.LT.1.25) THEN
             YYMAX = ZMAX
             YYMIN = YYMAX - 2.5
           ELSE
             YYMIN = ZXP - 0.5
             YYMAX = ZXP + 0.5
           ENDIF
        elseif (cgridopt.eq.RIBBON_GRID) then 
           ! jdemod - add some options for plotting ribbon grids
           XXMIN = RMIN-0.05
           XXMAX = RMAX+0.05
           YYMIN = ZMIN-0.1
           YYMAX = ZMAX+0.1
        else
           XXMIN = RXP - xnear
           XXMAX = RXP + xnear
           IF (ZXP-ZMINp.LT.ynear) THEN
             YYMIN = ZMINp
             YYMAX = YYMIN + 2.0 * ynear
           ELSE IF (ZMAXp-ZXP.LT.ynear) THEN
             YYMAX = ZMAXp
             YYMIN = YYMAX - 2.0 * ynear
           ELSE
             YYMIN = ZXP - ynear
             YYMAX = ZXP + ynear
           ENDIF

        endif
CX      XXMIN = 2.45
CX      XXMAX = 2.65
CX      YYMIN = 1.85
CX      YYMAX = 2.05
c
c       Add code to read plotting range for plot type from the
c       input line.
c
c slmod begin - new
        IF (graph(4:4).EQ.'&'.OR.graph(4:4).EQ.'Z'.OR.
     .      graph(4:4).EQ.'z'.OR.zmode.EQ.1) THEN

          IF (zmode.EQ.0) THEN
            IF (iref.EQ. 12.OR.iref.EQ.602.OR.iref.EQ.606.OR.
     .          iref.EQ.604.OR.iref.EQ.618.OR.
     .          iref.EQ.964.OR.iref.EQ.972.OR.iref.EQ.974) THEN
              READ(graph(20:24),*) xcen
              READ(graph(26:30),*) ycen
              READ(graph(32:35),*) xnear2
              READ(graph(37:40),*) ynear2
            ELSE
              READ(graph(24:28),*) xcen
              READ(graph(30:34),*) ycen
              READ(graph(36:39),*) xnear2
              READ(graph(41:44),*) ynear2
            ENDIF
          ENDIF

          xxmin = xcen - xnear2
          xxmax = xcen + xnear2
          yymin = ycen - ynear2
          yymax = ycen + ynear2
        ELSEIF (iref.EQ.974) THEN
          XPOINT= ' '
          ASPECT= 0.5 * (ZMAXp-ZMINp)
          XXMIN = 0.5 * (RMAX+RMIN) - ASPECT
          XXMAX = 0.5 * (RMAX+RMIN) + ASPECT
          YYMIN = ZMINp
          YYMAX = ZMAXp
        ENDIF
c slmod end
c
c        IXMIN = MAX (1  , INT((XXMIN-RMIN)/DR)+1)
c        IXMAX = MIN (NXS, INT((XXMAX-RMIN)/DR)+1)
c        IYMIN = MAX (1  , INT((YYMIN-ZMIN)/DZ)+1)
c        IYMAX = MIN (NYS, INT((YYMAX-ZMIN)/DZ)+1)
        
      elseif (cgridopt.eq.RIBBON_GRID) then 
        ! jdemod - add some options for plotting ribbon grids
        XXMIN = RMIN-0.05
        XXMAX = RMAX+0.05
        YYMIN = ZMIN-0.1
        YYMAX = ZMAX+0.1
      ELSE
        XPOINT= ' '
c       IPP/10 - added option to define coordinates in full plots
        if (zmode.eq.1) then
          xxmin = xcen - xnear2
          xxmax = xcen + xnear2
          yymin = ycen - ynear2
          yymax = ycen + ynear2
        else
          ASPECT= 0.5 * (ZMAXp-ZMINp)
          XXMIN = 0.5 * (RMAX+RMIN) - ASPECT
          XXMAX = 0.5 * (RMAX+RMIN) + ASPECT
          YYMIN = ZMINp
          YYMAX = ZMAXp
        endif
c        IXMIN = 1
c        IXMAX = NXS
c        IYMIN = 1
c        IYMAX = NYS
      ENDIF
C     WRITE (6,'(A,4I5)') ' IXMIN ETC ...',IXMIN,IXMAX,IYMIN,IYMAX


c slmod begin 
c...  This is here to initialize some postscript something when
c     saving individual plots from Ghostview.  Without this, the
c     text strings on the first plot are rotated:
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL FULL
      CALL LINCOL(1)	
      CALL CTRMAG(10)
      CALL PLOTST(0.0,0.0,'[i]')
c slmod end


      return 
c
c     Format statements
c


 9001 FORMAT(/1X,'OUT:    ',A55,'  OPT',I3,' REF ',A3,/1X,79('-'))

      end
c
c
c
      subroutine outinit
      use divimp_netcdf
      implicit none

c
      include 'params'
      include 'outcom'
c
c     Other common blocks
c

      include 'cgeom'
      include 'comtor'

c      include 'cneut2'
c      include 'dynam2'
c      include 'dynam3'
      include 'dynam4'
c      include 'pindata'
c      include 'cadas'
      include 'grbound'
      include 'outxy'
c      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
      include 'printopt' 
      include 'plot_switches'
c
      include 'out_unstruc' 
c
c     Local Variables
c
c
      integer ix,iy
      integer ik,ir,iz,it
      integer in,ii
      real r,z

c
c     Local Variables
c
      integer ierr
      character*80 graph
      ! jdemod - size needs to be coordinated with 
      !          max contents of raw file
      character*1024 desc

c
c     Initialization
c
c     Initialize unit number for data unit output (.dag or .daga)
c
      datunit = 7
c
c     Initialize the first_contour plot switch to false
c
      first_contour = .false.  
c
      PSWITCH = .FALSE.
      ATYPE = 0
      NUMSMOOTH = 3
      outgrid = .false.
      zadj = -1.780
c
c     Set ignors all non-zero
c   
      do ii = 1,maxngs 
         ignors(ii) = 1
      end do
c
c     Initialize the signal output
c
      call setup_signal_output

c
c     Calculate format string for printing
c
c      write(prform,'''(1X,2I3,2F7.3,1P,E8.2,'',i3,''E9.2)''') 2+maxizs
c

      call rzero(pltmins,maxplts)
      call rzero(pltmaxs,maxplts)
C
      CALL GPSTOP (100)
      CALL PAPER  (1)
c
c     System dependent printer initialization
c
      call printerinit
c
      CALL XUFLOW (0)
C
      REWIND (8)

   10 CONTINUE

c
c     Load case data from RAW data file
c

      CALL GET (TITLE,desc,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)
C
c
C-----------------------------------------------------------------------
c
c     SET INNER/OUTER Variables for print purposes
c
C-----------------------------------------------------------------------
c
c     Set up the meaning of the words INNER and OUTER in the context
c     of this specific case. An X-point up JET case has
c     INNER = 'INNER' and OUTER = 'OUTER' while and X-point down
c     case has this orientation reversed - this will be easier than
c     adding "IF" statements throughout repeating the test to define
c     INNER and OUTER.
c
c     X-point down - rings start numbering from INNER target clockwise
c
      if (zxp.le.z0) then
         INNER = 'OUTER'
         OUTER = 'INNER'
         xpoint_up = .false.
         inner_targid=2
         outer_targid=1 
c
c     X-point UP - clockwise from OUTER target - original default
c
      else
         INNER = 'INNER'
         OUTER = 'OUTER'
         xpoint_up = .true.
         inner_targid=1
         outer_targid=2 
      endif
c
c   
C
      REWIND (5)
      ierr = 0
c
c     Initialize any unstructured input values
c
      call init_out_unstruc_input
c
c     READ the input data for the case
c

      CALL RDC (GRAPH ,'TITLE LINE',IERR)
      call rdi (rxygrid,.TRUE.,-1,.true., 99,'XYGRID READ/CALC ',ierr)
      call rdi (cgrprint,.true.,0,.true.,1  ,'NUMERICS FOR PLOT',ierr)
      call rdi (plrpopt,.true.,0,.true.,1,      'PLRP data source',ierr)
      call rdr (alphae,.false.,0.0,.false.,0.0,'Exp Alpha value',ierr)
      CALL RDI (NAVS  ,.TRUE., 0 ,.TRUE.,100,'NUMBER OF WEIGHTS',IERR)
      CALL RDR (CZD ,.FALSE.,0.0,.FALSE.,0.0,'DIVERTOR LIMIT ZD',IERR)
      CALL RDI (ICNTR,.TRUE.,0,.TRUE.,   2,  'CONTOUR PLOT TYPE',IERR)
      call rdi (cntropt,.true.,0,.true.,  4,'OPTION TO CALC CONT',ierr)
      call rdi (global_cngs, .true.,1,.false.,  0,'NUM. OF CONTOURS ',
     >                                  ierr)
      CALL RDRAR(conts,nconts,maxpts,0.0,HI,.TRUE.,'Contours',IERR)
      call rdi (clsup,.true.,0,.true.,  1 ,  'CLOSEUP PLOTS N/Y',ierr)
      call rdr (magstep,.true.,0.0,.true.,1.0,'MAG. FOR CLOSEUP',ierr)
      CALL RDR (xnear,.true.,0.0,.FALSE.,0.0,'X PLOTTING LIMIT',IERR)
      CALL RDR (ynear,.true.,0.0,.FALSE.,0.0,'Y PLOTTING LIMIT',IERR)
      call rdr (scalef,.false.,0.0,.false.,0.0,'Scaling Factor',IERR)
      call rdr (zadj,.false.,0.0,.false.,0.0,'Z-adj for R-plots',ierr)
      call rdi (iseldef,.true.,0,.false.,0,'DEFAULT EXPT DATASET',ierr)
c
c     Input error exit condition
c
      IF (IERR.NE.0) then 
         write(0,*) 'OUTINIT: ERROR READING OUT INPUT DATA:'//
     >                      ' PROGRAM EXITING'
         write(6,*) 'OUTINIT: ERROR READING OUT INPUT DATA:'//
     >                      ' PROGRAM EXITING'
         stop 1
      endif
c
c     Set absfac to value read in from unstructured input - if specified
c
      if (new_absfac.gt.0.0) then 
c
c        jdemod - write out a warning message
c
         write(0,*) 'OUTINIT: NEW ABSOLUTE SCALING FACTOR SPECIFIED =',
     >                      new_absfac
         write(6,*) 'OUTINIT: NEW ABSOLUTE SCALING FACTOR SPECIFIED =',
     >                      new_absfac
c
         absfac=new_absfac
c
      endif  

c
c     Depending on the value of SCALEF - assign a final value
c
c     IF SCALEF=1.0  - do nothing - this is the default value.
c     IF SCALEF=-1.0 - set scale factor to ABSFAC
c
c     IF SCALEF>0    - leave scaling factor at this value
c     IF SCALEF<0    - set scaling factor to ABSFAC/SCALEF
c
c     Adjust absfac - if it is equal to 0.0 - set it to 1.0 for
c     plotting purposes.
c
      if (absfac.eq.0.0) absfac = 1.0
c
      if (scalef.eq.-1.0) then
         scalef=absfac
      elseif (scalef.lt.0.0) then
         scalef= absfac/abs(scalef)
      endif
c
      write (6,*) 'SCALE FACTOR = ',scalef
      write (6,*) 'ABSFAC       = ',absfac
      write (iplot,*) 'SCALE FACTOR = ',scalef
      write (iplot,*) 'ABSFAC       = ',absfac


c
c     jdemod - if the netcdf output flag has been set at the beginning of the plot
c              file then write the netcdf file since all the required inputs 
c              should be available here. 
c
C
C-----------------------------------------------------------------------
C     GENERATE NETCDF VERSION OF RAW OUTPUT
C-----------------------------------------------------------------------
C
      if (netcdf_opt.eq.1) then 
          call write_netcdf_output(TITLE,desc,NIZS,JOB,EQUIL,
     >                             FACTA,FACTB,ITER,NITERS)
         
      endif


c
c     Set up colors - from dark to light - one for each contour - let
c     colour 1 be black - use the HSI system from GHOST
c
      if (cntropt.eq.0.or.cntropt.eq.2) then
         n_cols = global_cngs +1
      elseif (cntropt.eq.1) then
         n_cols = 10
      elseif (cntropt.eq.3.or.cntropt.eq.4) then
         n_cols = nconts +1
      endif
c
c     The following routine sets up and initializes the number of colours
c     specified.
c
      if (icntr.eq.0.or.icntr.eq.1) then
         col_opt = 1
      elseif (icntr.eq.2) then
         col_opt = 2
         icntr = 1
      endif
c
C     IPP/01 Krieger - bug in color selection fixed
c
      if (cntropt.eq.3.or.cntropt.eq.4) then
        col_opt = 3
      endif
c
      call setup_col(n_cols,col_opt)
c
c     Set magnification related values
c
      mgst = magstep
      mgnd = 1.0 - magstep
C
C
C     ADDITIONAL INITIALIZATION - FOR BACKGROUND PLOTS
C
      RIZB = REAL (CIZB)
C
      IF (ITER.EQ.1) THEN
        WRITE (6,9020) VERSON,TITLE(1:58),JOB(36:72),equil(1:54)
      ELSE
        WRITE (6,9021) TITLE(61:80)
      ENDIF
c
c     Temporary write statement
c
c      write (6,*) 'Table of Hydrogenic characteristics:'
c      write (6,*) 'IR  IK   Te (eV)   Ti (eV)    ne   nD+  nD0  nD2'
c
c      do ir = irsep,irwall-1
c         do ik = 1,nks(ir)
c            write (6,'(2i5,6(1x,g12.5))') ir,ik,ktebs(ik,ir),
c     >         ktibs(ik,ir),knbs(ik,ir),knbs(ik,ir),pinatom(ik,ir),
c     >         pinmol(ik,ir)
c         end do
c      end do
c
c
c      Generate the PLRP's based on the selected input (plrpopt)
c      The call to the PLRP module has been moved from DIVIMP and
c      the contents of SDLIMS have been substituted for DDLIMS.
c
c      Since the PLRP routines worked on the NON-REFLECTED DIVIMP
c      grid - this processing has been placed into OUT before the
c      reflection occurs.
C
C-----------------------------------------------------------------------
C   CALCULATE PARTICULAR NE RADIATION
C-----------------------------------------------------------------------
C   added rizb to argument list; Krieger, IPP 95
C
      CALL PLRP (NIZS,PLAMS,PIND,PIZS,PLRPCNT,CION,rizb,plrpopt)
C
C     FIND HOW MANY NEUTRAL PLRP LINES THERE ARE:
C
      PSHIFT = 1
      PNCNT = 0
      DO 5 IZ = -1,PLRPCNT
        IF (PIZS(IZ).EQ.0) PNCNT = PNCNT +1
 5    CONTINUE

c
C-----------------------------------------------------------------------
C   CALCULATE MEAN FREE PATHS AND CHARACTERISTIC SCALE LENGTHS
C-----------------------------------------------------------------------
c
      call calc_mfp(lgradti,lgradte,lmfpii,lmfpee)
c
C-----------------------------------------------------------------------
c     REFLECT GEOMETRY - IF REQUIRED
C-----------------------------------------------------------------------
C
C     IF REFCT=1 THEN GEOMETRY NEEDS TO BE REFLECTED IN R=0 PLANE
C
      IF (REFCT.EQ.1) THEN
        CALL REFLECT
        WRITE(6,'('' EQUILIBRIUM GEOMETRY HAS BEEN REFLECTED'')')
      ENDIF
c
C
C-----------------------------------------------------------------------
c     Calculte krb, kzb and ksb for DIVIMP runs
C-----------------------------------------------------------------------
c
c     NOTE: Calculation of KSB and KPB are no longer required since 
c           these are now included in the raw data file - at least for
c           rings with polygons. 
c
c
      if (cgridopt.eq.0.or.cgridopt.eq.3) then
c
c        Only possible for grids with polygon information
c
         DO IR = 1, NRS
c
c           Calculate beginning of ring.
c
            in = korpg(1,ir)
c
c           If a polygon exists for the cell.
c
c           The code assumes that all cells will exist for a given
c           ring and that there won't be gaps in rings without defined
c           polygons. If the first cell doesn't have a polygon it just
c           uses the KSS values found previously for KSS2.
c
c
            if (in.ne.0) then
c
               in = korpg(1,ir)
c
               krb(0,ir) = (rvertp(1,in) + rvertp(2,in)) /2.0
               kzb(0,ir) = (zvertp(1,in) + zvertp(2,in)) /2.0
c               ksb(0,ir) = 0.0
c               kpb(0,ir) = 0.0
c
               do ik = 1,nks(ir)
c
                  in = korpg(ik,ir)
c
                  krb(ik,ir) = (rvertp(3,in) + rvertp(4,in)) /2.0
                  kzb(ik,ir) = (zvertp(3,in) + zvertp(4,in)) /2.0
c
c                  ksb(ik,ir) = sqrt((krb(ik,ir)-krb(ik-1,ir))**2+
c     >                              (kzb(ik,ir)-kzb(ik-1,ir))**2)
c     >                         * kbfs(ik,ir) + ksb(ik-1,ir)
c
c                  kpb(ik,ir) = sqrt((krb(ik,ir)-krb(ik-1,ir))**2+
c     >                              (kzb(ik,ir)-kzb(ik-1,ir))**2)
c     >                              + kpb(ik-1,ir)
c
c                  write(6,'(a,3i4,1p,6(1x,g12.5))') 'KSS:',ik,ir,in,
c     >                 ksb(ik,ir),kss(ik,ir),ksb(ik-1,ir)
c
c
               end do

c
c           For rings without polygons defined - KSS2 -> KSS and the
c           KSB values -> (KSS(IK+1,ir)-KSS(ik,IR))/2.0
c
            elseif (in.eq.0) then
c
               ksb(0,ir) = 0.0
               kpb(0,ir) = 0.0
c
               if (ir.ge.irsep) then
                  krb(0,ir) = rp(idds(ir,2))
                  kzb(0,ir) = zp(idds(ir,2))
               else
                  krb(0,ir) = (rs(1,ir) + rs(nks(ir)-1,ir))/2.0
                  kzb(0,ir) = (zs(1,ir) + zs(nks(ir)-1,ir))/2.0
               endif
c
c              Loop through ring
c
               do ik = 1,nks(ir)
c
c                 Set cell boundaries for the ring without polygons.
c
                  if (ik.eq.nks(ir)) then
                     ksb(ik,ir) = ksmaxs(ir)
                     kpb(ik,ir) = kpmaxs(ir)
                     if (ir.ge.irsep) then
                        krb(ik,ir) = rp(idds(ir,1))
                        kzb(ik,ir) = zp(idds(ir,1))
                     else
                        krb(ik,ir) = krb(0,ir)
                        kzb(ik,ir) = kzb(0,ir)
                     endif
                  else
                     ksb(ik,ir) = (kss(ik+1,ir)+kss(ik,ir))/2.0
                     kpb(ik,ir) = (kps(ik+1,ir)+kps(ik,ir))/2.0
                     krb(ik,ir) = (rs(ik+1,ir) +rs(ik,ir))/2.0
                     kzb(ik,ir) = (zs(ik+1,ir) +zs(ik,ir))/2.0
                  endif
c
               end do

            endif
c
         end do
c
      endif
C
C-----------------------------------------------------------------------
C     RE-CALCULATE THE CTIMES ARRAY FOR TIME DEPENDENT RUNS
c     Leave this in terms of seconds for plotting purposes!
C-----------------------------------------------------------------------
c
      call rzero(ctimes,(maxnts+2)*(maxizs+2))
      CSTMAX = 0
      IF (IMODE.NE.2) THEN
        DO IT = 1, NTS
          CTIMES(IT,-1) = DWELTS(-1)*DWELFS(IT)
          CTIMES(IT, 0) = DWELTS( 0)*DWELFS(IT)
          CSTMAX = MAX (CTIMES(IT,-1),CTIMES(IT,0))
          IF (NIZS.GT.0) THEN
            DO IZ = 1, NIZS
              CTIMES(IT,IZ) = DWELTS(IZ)*DWELFS(IT)
              CSTMAX = MAX (CSTMAX, CTIMES(IT,IZ))
            end do
          ENDIF
        end do
c
c       Set a zero'th time point for calculating time-bin widths easily
c
        do iz = -1,nizs
c
           ctimes(0,iz) = 0.0
c
        end do
c
      ENDIF
C
C---- TAKE MAX OF THIS VALUE WITH CTIMMAX SECONDS, IF STEADY STATE REQUI
C
      IF (IMODE.NE.1) CSTMAX = MAX (CSTMAX, 10.0)
C
C---- INTRODUCE AN EXTRA TIMEPOINT WHICH WILL NEVER BE REACHED
C---- THIS TIMEPOINT WILL APPLY EQUALLY TO NON-IMPULSE MODE CASES
C
      CTIMES(NTS+1,-1) = 2.0 * CSTMAX
      CTIMES(NTS+1, 0) = 2.0 * CSTMAX
      IF (NIZS.GT.0) THEN
        DO  IZ = 1, NIZS
          CTIMES(NTS+1,IZ) = 2.0 * CSTMAX
        end do
      ENDIF
c
c     Print out times array
c
      write (6,*) 'Imode = ', imode
      write (6,*) 'CTIMES:'
      do it = 0,nts+1
         write (6,'(1p,10(1x,g12.4))') (ctimes(it,iz),iz=0,nizs)
      end do

c
C-----------------------------------------------------------------------
c
c     Read in input related to calculating the Dperp values in OUT
c
c      call rdi (outdp,.true.,0,.true.,1,'CALC D/XPERP IN OUT',ierr)
c
c      if (outdp.eq.1) then
c
c     Dperp extractor options
c
c      CALL RDI(dpmethod,.TRUE.,0  ,.true.,2 ,'DPERP EXT METHOD',IERR)
c      call rdi(dpsuml,.TRUE.,0  ,.true.,1 ,'DPERP EXT SUM LIMIT',IERR)
c      call rdi(dpouter,.TRUE.,0  ,.true.,1 ,'DP EXT OUTER RING',IERR)
c      call rdi(dpconv,.TRUE.,0  ,.true.,1 ,'DP EXT CONVECT LOSS',IERR)
c      call rdi(dpfluxopt,.TRUE.,0  ,.true.,1 ,'CELL CENTRE FLUX',IERR)
c      call rdi(dpavopt,.TRUE.,0  ,.false.,1 ,'AVERAGE DPERP OPT',IERR)
c      call rdi(dprcopt,.TRUE.,0  ,.true.,2 ,'MAJOR RADIUS CORR',IERR)
c      call rdi(dpsmooth,.true.,0,.false.,0 ,'GRADIENT SMOOTHING',IERR)
c      call rdi(dpnav,.true.,-1,.true.,2 ,   'GRADIENT CALC METH',IERR)
c      call rdi(dparea,.TRUE.,0  ,.true.,1 ,'CELL BOUND AREA',IERR)
c      call rdi(dpploss,.TRUE.,0  ,.true.,1 ,'POWER LOSS TERMS',IERR)
c
c     Call the dperp extractor routine.
c
c      call oskin
c
c      endif
c
c     Read and/or generate XY GRIDS if required.
c
c     If XY grids were off in DIVIMP - OUT will need to either
c     calculate or load them. Or for JET grids give it the option of
c     using the method in DIVIMP of calculating the position in the
c     array ... this could be more computationally intensive. Allow this
c     if the grid option for OUT is specified as -1 .. leave the OUTGRID
c     variable set to .false. so that the gridpos routine thinks that
c     it is dealing with DIVIMP ... this can only be done where the grid
c     and their vertices are well-defined. (i.e. available)
c
      if (xygrid.eq.0) then
C
C-----------------------------------------------------------------------
C     CALCULATE RECTANGULAR GRID EQUIVALENTS TO IK, IR INDICES
C     CALCULATE SET OF FLAGS IFXYS INDICATING PTS EXTERNAL TO SYSTEM.
C-----------------------------------------------------------------------
C
c     Read grid from file.
c
      IF (rxygrid.EQ.99) THEN
        WRITE (6,'(/1X,A)') 'RECTANGULAR GRID DATA READ FROM FILE :-'
        REWIND (13)
        READ (13,*) NXS,NYS,DR,DZ,RMIN,RMAX,ZMIN,ZMAX
        write(6,*)  NXS,NYS,DR,DZ,RMIN,RMAX,ZMIN,ZMAX,MAXGXS,MAXGYS
        IF (NXS.GT.MAXGXS) WRITE (6,*) ' ERROR! NXS =',NXS
        IF (NYS.GT.MAXGYS) WRITE (6,*) ' ERROR! NYS =',NYS
        READ (13,*) ((IKXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        READ (13,*) ((IRXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        READ (13,*) ((IFXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        outgrid = .true.
C
      ELSEif (rxygrid.ge.0) then
c
c       First - go through and calculate all points inside the
c       outer wall - then find all points inside the central plasma
c       region - using the Harwell area routine.
c
        WRITE (6,'(/1X,A)') 'RECTANGULAR GRID CALCULATED AND STORED :-'
        NXS = MAXGXS
        NYS = MAXGYS
        DR = (RMAX-RMIN) / REAL(NXS-4)
        DZ = (ZMAX-ZMIN) / REAL(NYS-4)
        RMIN = RMIN - 2.0 * DR
        RMAX = RMAX + 2.0 * DR
        ZMIN = ZMIN - 2.0 * DZ
        ZMAX = ZMAX + 2.0 * DZ
C
        DO 20 IX = 1, NXS
          R = (RMAX-RMIN) * REAL(IX)/REAL(NXS) + RMIN - 0.5 * DR
          DO 20 IY = 1, NYS
            Z = (ZMAX-ZMIN) * REAL(IY)/REAL(NYS) + ZMIN - 0.5 * DZ
            call gridpos(ik,ir,r,z,.false.,griderr)
            if (griderr) then
               ifxys(ix,iy) = 0
               ikxys(ix,iy) = 0
               irxys(ix,iy) = 0
            else
               ifxys(ix,iy) = 1
               ikxys(ix,iy) = ik
               irxys(ix,iy) = ir
            endif
c
 20     continue
C
C
        outgrid = .true.
c
c
        REWIND (13)
        WRITE (13,*) NXS,NYS,DR,DZ,RMIN,RMAX,ZMIN,ZMAX
        WRITE (13,'(1X,23I3)') ((IKXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        WRITE (13,'(1X,23I3)') ((IRXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        WRITE (13,'(1X,23I3)') ((IFXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        WRITE (6,*) NXS,NYS,DR,DZ,RMIN,RMAX,ZMIN,ZMAX
c
c        WRITE (6,'(2I8,'':'',3I8)') (( ix,iy,IFXYS(Ix,Iy),
c     >            IKXYS(IX,IY),
c     >            IRXYS(Ix,Iy),IX=1,NXS),IY=1,NYS)
c
      elseif (rxygrid.eq.-1) then
c
c       Set up the base values required even if the rectangular
c       grids are not calculated - since the GHOST plotting routines
c       are based on an underlying X,Y grid ... the size of this
c       grid and the range of R,Z values are required - even if the
c       locations of the X,Y points are not pre-mapped to the
c       IK,IR bins.
c
        WRITE (6,'(/1X,A)') 'RECTANGULAR GRID NOT CALCULATED :-'
        NXS = MAXGXS
        NYS = MAXGYS
        DR = (RMAX-RMIN) / REAL(NXS-4)
        DZ = (ZMAX-ZMIN) / REAL(NYS-4)
        RMIN = RMIN - 2.0 * DR
        RMAX = RMAX + 2.0 * DR
        ZMIN = ZMIN - 2.0 * DZ
        ZMAX = ZMAX + 2.0 * DZ
        outgrid = .false.
c
      ENDIF
c
      endif



c
c     Calculate the vessel wall for neutrals used in the case
c     based on the contents of the wallpts array.
c
       
      call calc_neutralwall

c
      return 
c
c     Format statements
c


 9020 FORMAT(/1X,62('*'),/1X,'*',60X,'*',
     >  /1X,'*',18X,'RUN OF OUT VERSION ',A5,18X,'*',
     >  /1X,'*',18X,24('-'),18X,'*',/1X,'*',60X,'*',
     >  /1X,'* ',A58,' *',/1X,'*',60X,'*',/1X,'* ',A53,5X,' *',
     >  /1X,'*',60X,'*',/1X,'* ',A54,4X,' *',
     >  /1X,'*',60X,'*',/1X,62('*'),/)
 9021 FORMAT(/1X,62('*'),/1X,'*',60X,'*',
     >  /1X,'*',18X,A20,22X,'*',/1X,'*',60X,'*',/1X,62('*'),/)


      end
c
c
c
      subroutine global_plotsetup 
c slmod begin
      use mod_out985
      use mod_out985_variables
c slmod end
      implicit none

c
      include 'params'
      include 'outcom'
c
c     Other common blocks
c
      include 'cgeom'
      include 'comtor'
c      include 'cneut2'
c      include 'dynam2'
c      include 'dynam3'
c      include 'dynam4'
      include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
c      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
c      include 'printopt' 
c
c     Initialize some global graph parameters
c
      include 'comgra'
c
c     Local Variables
c
c
      integer i,j
      integer ix,iy
      integer id,ig,jd
      integer ik,ir,iz

c
      real np,ns,nt 

c
c     Local Variables
c

      integer startid,endid,stepid,switchid

c
c     Initialize the line base line thickness for plots
c     This is in the common block 'comgra'
c
      thickness = 1
c


      ITEC = 1
      DO 15 IG = 0, NAVS
        AVS(IG) = REAL (2*NAVS-IG)
   15 CONTINUE
      WRITE (6,'('' SMOOTHING BY WEIGHTED AVERAGES WITH FUNCTION :-'')')
      WRITE (6,'((15X,7(F7.3,:,'' :'')))') (AVS(IABS(J)),J=-NAVS,NAVS)
C
      NP = 0.0
      NT = 0.0
      IF (FACTA(-1).GT.0.0) NP = 1.0 / FACTA(-1)
      IF (FACTA(0) .GT.0.0) NT = 1.0 / FACTA(0)
      NS = NT - NP
      IF (NS.GT.0.001) THEN
        FP = NP / NS
        FT = NT / NS
      ELSE
        FP = 0.0
        FT = 0.0
      ENDIF
      WRITE (6,'(/,'' OUT: NP,NS,NT,FP,FT='',5F10.3,/)') NP,NS,NT,FP,FT
c
c     Calculates and writes out some numbers - not apparently used elsewhere since they 
c     are locals
c
      call global_hc_plot_setup 
C
      NVIEW  = ' '
      PLANE  = ' '
      ANLY   = ' '
      TABLE  = 'SYMBOL TABLE'
      SMOOTH = ' '
c
      ZLABS(-2)= '  S SECNEUT'
      ZLABS(-1)= 'P   PRINEUT'
      ZLABS(0) = ' T  TOTNEUT'
      WRITE (PLABS(-2),'(A7,I4)') '  S SEC',INT(PLAMS(-1))
      IF (PLRPCNT.GT.0) THEN
        DO 35 IZ = -1,PLRPCNT
          IF (PIZS(IZ).EQ.-1) THEN
            WRITE (PLABS(IZ),'(A7,I4)') 'P   PRI',INT(PLAMS(IZ))
          ELSEIF (PIZS(IZ).EQ.0) THEN
            WRITE (PLABS(IZ),'(A7,I4)') ' T  TOT',INT(PLAMS(IZ))
          ELSEIF  (3*(PIZS(IZ)/3).EQ.PIZS(IZ)) THEN
            WRITE(PLABS(IZ), '(I4   ,I2,I5)') PIZS(IZ),
     >            PIZS(IZ),INT(PLAMS(IZ))
          ELSEIF (2*(PIZS(IZ)/2).EQ.PIZS(IZ)) THEN
            WRITE(PLABS(IZ), '(I3,1X,I2,I5)') PIZS(IZ),
     >            PIZS(IZ),INT(PLAMS(IZ))
          ELSE
            WRITE(PLABS(IZ), '(I2,2X,I2,I5)') PIZS(IZ),
     >            PIZS(IZ),INT(PLAMS(IZ))
          ENDIF
          WRITE(6,*) 'PLABS:',IZ,'(',PLABS(IZ),')',PLAMS(IZ),PIZS(IZ)
 35     CONTINUE
      ENDIF
C
      IF (NIZS.GT.0) THEN
        DO 40 IZ = 1, NIZS
          IF     (3*(IZ/3).EQ.IZ) THEN
            WRITE(ZLABS(IZ), '(I4   ,A5,I2)') IZ,'IONIZ',IZ
          ELSEIF (2*(IZ/2).EQ.IZ) THEN
            WRITE(ZLABS(IZ), '(I3,1X,A5,I2)') IZ,'IONIZ',IZ
          ELSE
            WRITE(ZLABS(IZ), '(I2,2X,A5,I2)') IZ,'IONIZ',IZ
          ENDIF
   40   CONTINUE
      ENDIF
      ZLABS(NIZS+1) = 'A   ALL IZS'
C
      HLABS(0) = ' T  TOTNEUT'
      HLABS(1) = '  1 IONIZ 1'
      HLABS(2) = 'A   ALL IZS'
C
      BLABS    = ' B  BREMS. '
C
      DO 50 IX = 1, NXS
        XOUTS(IX) = RMIN + DR * (REAL(IX)-0.5)
        XWIDS(IX) = DR
   50 CONTINUE
      DO 60 IY = 1, NYS
        YOUTS(IY) = ZMIN + DZ * (REAL(IY)-0.5)
        YWIDS(IY) = DZ
   60 CONTINUE
      IF (NDS+2.GT.MAXNDS) WRITE (6,*) ' ERROR! MAXNDS NOT >= ',NDS+2
      JD = 0
      WRITE (6,'('' NDSIN,NDS,NRS,NKS()='',(1X,30I3))')
     >              NDSIN,NDS,NRS,(NKS(IR),IR=1,NRS)
      WRITE (6,'('' IRSEP,IRWALL,IRTRAP,IKT,IKTO,IKTI,IKREF='',7I3)')
     >              IRSEP,IRWALL,IRTRAP,IKT,IKTO,IKTI,IKREF
c
c     Set IKT values correctly for various geometries ... this
c     may be temporary until the meaning of the IKT values
c     is well-defined for each geometry.
c
c      if (cgridopt.eq.0.or.cgridopt.eq.3) then
c         ikto = ikt
c         ikti = ikt
c      endif
c
c
c slmod begin
      if (nds.eq.0) then
c...    In case NDS is not set:
        startid = 1
        endid = 0
        stepid = 1
      elseif (rp(1).gt.rp(nds)) then
c
c      if (rp(1).gt.rp(nds)) then
c slmod end
         startid = nds
         endid   = 1
         stepid  = -1
         switchid = ndsin +1
      else
         startid = 1
         endid   = nds
         stepid = 1
         switchid = ndsin
      endif
c
      jd = 0
c
      write (6,*) 'Calc douts:',switchid,ndsin,startid,endid
c
      DO 70 ID = startid, endid, stepid
        JD = JD + 1
        DOUTS(JD) = rp(id)
        DWIDS(JD) = DDS(ID)
c        write (6,*) id,jd,douts(jd),rp(id)
        IF (ID.EQ.switchid) JD = JD + 2
   70 CONTINUE
c
c     Insert points to space between targets
c
c     If there are ever cases where there are
c     different number of target points on each target
c     then the following code will need to be changed.
c
c slmod begin
      if (grdnmod.ne.0) 
     .  WRITE(0,*) 'NOTE: UPDATE CODE FOR GENERALIZED '//
     .             'GRIDS'
      if (ndsin.gt.0) then
        DOUTS(NDSIN+1) = DOUTS(NDSIN)   + 0.001
        DWIDS(NDSIN+1) = 0.001
        DOUTS(NDSIN+2) = DOUTS(NDSIN+3) - 0.001
        DWIDS(NDSIN+2) = 0.001
      endif
c
c      DOUTS(NDSIN+1) = DOUTS(NDSIN)   + 0.001
c      DWIDS(NDSIN+1) = 0.001
c      DOUTS(NDSIN+2) = DOUTS(NDSIN+3) - 0.001
c      DWIDS(NDSIN+2) = 0.001
c slmod end
c
c     Print DOUTS
c
      write(6,*) 'Douts:',nds,ndsin
c
      do jd = 1,nds + 2
         write (6,*) jd,douts(jd)
      end do
c
c     Mirror the first two targets in the Z-axis for
c     ITER grids - just so all four targets can be plotted.
c
      if (cgridopt.eq.2) then
         DO 74 ID = 1, NDSIN2+2
           DOUTS(ID) = -DOUTS(ID)
   74    CONTINUE
      endif
c
c
c     Calculate plot ranges.
c
      zminp = zmin
      zmaxp = zmax
      if (nvesm.ne.0) then
        do i = 1, nvesm+nvesp
          zminp = amin1(zminp,zvesm(i,1),zvesm(i,2))
          zmaxp = amax1(zmaxp,zvesm(i,1),zvesm(i,2))
        enddo
      elseif (nves.ne.0) then
        do i = 1, nves
          zminp = min(zminp,zves(i))
          zmaxp = max(zmaxp,zves(i))
        enddo
      endif
c
      write (6,*) 'nvesm,nves:',nvesm,nves,zminp,zmin,zmaxp,zmax
c
c     Print out table of impurity content between given limits on
c     the SOL rings.
c
c      call calcnt(nizs)



      NPLOTS = 0
c slmod begin
      nobj      = 0
      stepopt   = 0
      nsteplist = 0
      mode      = 0
c slmod end





c 
      return 
      end 
c
c
c
      subroutine plotloopinit(iopt,ierr)
      implicit none
      integer iopt,ierr

c
      include 'params'
      include 'outcom'
      include 'comtor' 
      include 'adas_data_spec' 
c
c     Local Variables
c
      integer ir,ii


c slmod begin
c
c     Initialization
c
      pltfact = 0.0
      loadstep   = -1
      qt         = qtim
      zmode      = 0
      iopt       = 0
      iopt_ghost = 0
      grm_opt    = 0
      slopt      = 0
      slopt2     = 0
      slopt3     = 0
      slopt4     = 0
      slopt5     = 0
      losopt     = 0
      opt_xscale = 0
      opt_yscale = 0
      rel_step   = 0

      map1x = 0.10
      map2x = 0.90
c slmod begin
c...  True space
      map1y = 0.10
      map2y = 0.90
c      map1y = 0.11
c      map2y = 0.89
c slmod end
c
c     Initialize ADAS data switch to 0
c
      cadas_switch = 0
c
c     Initialize number of experimental datasets for this plot to zero
c
      expt_nsets = 0
c
      DO ir = 1, MAXNGS
        grm_shade(IKLO,ir) = 0.0
        grm_shade(IKHI,ir) = 0.0
      ENDDO
      CALL RZero(grm_cell(0,1),(MAXNKS+1)*MAXNGS)

c...  Shoddy business : should not have 30 hard coded...
      DO ii = 1, 30
        char (ii) = ' '
        ylab2(ii) = ' '
      ENDDO

      DO ii = 1, MAXNGS
        plottype(ii) = 1
        plotnorm(ii) = ii
        CALL ISet(plottype2(1,ii),8,1)
      ENDDO

c...  For colour overheads:
c     CALL RGB
c     ncols = ncols + 1
c     CALL ColSet(1.0,0.0,0.0,ncols+1)
c     CALL ColSet(0.0,1.0,0.0,ncols+2)
c     CALL ColSet(0.0,0.0,1.0,ncols+3)
c     CALL ColSet(0.5,0.7,0.2,ncols+4)
c     CALL ColSet(0.5,0.0,0.5,ncols+5)
c     CALL ColSet(0.0,0.5,0.5,ncols+6)
c     ncols = ncols - 1
c slmod end

      ierr = 0

      return 
      end
c
c
c
      subroutine load_additionalplotdata(iref,graph,iopt,ierr)
      implicit none
      integer iref,iopt,ierr
      character*(*) graph
c
      include 'params'     
      include 'outcom'
c
c     NOTE: This routine MUST contain some call to an RDG routine that 
c           will allow it to catch optional input data defined by the 
c           "#" in the input files. The RDG routines look for this and
c           will store the data correctly in OUTCOM. At present RDG7 is
c           performing this function - though this may change if 
c           RDG7 is phased out. 
c

c
c     Local variables
c
      character*80 label,graph6
c


c...  Look for zoom data:
      READ(5,'(A80)',END=150) graph6
      IF (graph6(8:11).EQ.'Zoom'.OR.graph6(8:11).EQ.'ZOOM'.OR.
     .    graph6(8:11).EQ.'zoom') THEN
        READ(graph6,*) label,xcen,ycen,xnear2,ynear2
        zmode = 1
        slopt = 1
      ELSE
        BACKSPACE 5
      ENDIF
150   CONTINUE
c slmod end

c
c kkmod     
c
c     For countour plots read min/max scale to be used if present
c     Krieger IPP/97
c
c     Try to read lines keyed with 001 for all plots - should
c     only be present and apply to contour plots. This is expected
c     to be the next line after the contour plot line. If present
c     for other plots the values will be loaded but not 
c     utilized. (Set number of contour lines)
c
      cngs = global_cngs      
c
c     Ignore any errors encountered in RDG7 since it could be 
c     an EOF condition and the input is optional. 
c
      call rdg7(graph,minscale,maxscale,localcngs,ierr)
      ierr = 0
c
      if (localcngs.ne.0) then 
         cngs = localcngs
      endif 
c
c
c kkmod


c
C
C  READ ADAS SELECTOR INPUT, IF REQUIRED
C
      IF (IREF.EQ.125 .OR. IREF.EQ.126 .OR.
     >    IREF.EQ.127 .OR. IREF.EQ.128 .OR.
     >    IREF.EQ.135 .OR. IREF.EQ.136 .OR.
     >    IREF.EQ.137 .OR. IREF.EQ.138 .OR.
     >    iref.eq.174 .or.
     >    iref.eq.176 .or. iref.eq.177 .or.
     >    IREF.EQ.212 .OR. IREF.EQ.214 .OR.
c slmod begin - new
     >    IREF.EQ.216 .OR. IREF.EQ.218 .OR. IREF.EQ.219 .OR.
     >    IREF.EQ.950 .OR. IREF.EQ.988 .OR.
c
c     >    IREF.EQ.216 .OR. IREF.EQ.218 .OR.
c slmod end
     >    IREF.EQ.246 .OR. IREF.EQ.266 .OR.
     >    IREF.EQ.276) THEN
        CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >             ISELE,ISELR,ISELX,ISELD,IERR)
        IF (IERR.NE.0) THEN
          WRITE(6,*) 'ERROR READING ADAS DETAILS, IREF = ',IREF
          IERR = 1
          return
        ENDIF
      ENDIF


c
c
c     Read multiple lines of ADAS selector input with ion type and
c     charge state for calculating ratio contours and LOS plots.
c     The ion type must correspond to either the background or
c     the impurity for the case in question.
c
      if (iref.eq.123.or.iref.eq.124.or.
     >    iref.eq.217) then
c
c        Read in two sets of ADAS data information and allow
c        for different species so that the ratios of carbon to
c        hydrogen lines may be plotted as well.
c
         call rdg5(graph5,adasid,adasyr,adasex,
     >             isele,iselr,iselx,iseld,iz_state,z_atom,ierr)
         call rdg5(graph5,adasid2,adasyr2,adasex2,
     >             isele2,iselr2,iselx2,iseld2,iz_state2,z_atom2,ierr)
c
      endif

c
c     Read additional plot data for certain HC plots
c
      if (iref.eq.513.or.iref.eq.514) then 

         call global_hc_read_additional_data(iref,graph,iopt,ierr)

      endif
c
      return 
      end

c
c
c
      subroutine calc_neutralwall
      implicit none
      include 'params'
      include 'comtor'
      include 'grbound'
c
c     Use the contents of wallpts to construct the neutral wall boundary
c     used for this simulation.
c
      integer in,kind
c 
      do in = 1,wallpts
c       
         rnw(in) = wallpt(in,20)
         znw(in) = wallpt(in,21)
c
      end do
c
c     Connect wall back to the beginning 
c     
      rnw(wallpts+1) = rnw(1) 
      znw(wallpts+1) = znw(1) 
c
c     Assign total number of points in the wall 
c
      neutwpts = wallpts+1 
c
c     Calculate boundary finding workspace  
c
      kind = 1
      CALL GA15A(neutwpts,KIND,nwWORK,4*MAXPTS,nwINDW,MAXPTS,
     >             RNW,ZNW,nwTDUM,nwXDUM,nwYDUM,6)
c
      return
      end
c
c
c 
      subroutine init_out_unstruc_input
      implicit none
c 
c     This routine assigns the default values to any
c     OUT unstructured input values.The OUT unstructured
c     input tags start with the letter "O" (oh)
c
c     jdemod - it also assigns default values to any other 
c              unstructured input values shared by DIV and OUT
c
c
c     Note: The netcdf output option A07 from DIVIMP is also supported
c           since it would be useful to create a netcdf version of the 
c           raw output from an OUT run
c
      include 'params'
c
      include 'out_unstruc'
      include 'comtor'
c     
c------------------------------------------------------
c
c     O01 - alternate absolute factor specification 
c         - this option is deactivated by setting the 
c           default value to zero.  
c
      new_absfac=0.0
c
c
c------------------------------------------------------
c
c     Core fueling code calculates integrated ionization
c     profiles in the core ... these parameters allow the 
c     PSIN inner bound of the integration regions to be set
c     This is used in the pr_eirene_analysis routine
c
c     O02 - PSIN bound for calculating core ionization 
c           profile 1 (psi1_reg)
c     O03 - PSIN bound for calculating core ionization 
c           profile 2 (psi2_reg)
c

      psi1_reg = 0.9
      psi2_reg = 0.95

c
c -----------------------------------------------------------------------
c
c     TAG A07: 
c
c
c     Option to write a netcdf version of the raw data file 
c     0 = off   1 = 0n .. default OFF
c
c
      netcdf_opt = 0
c
c
c------------------------------------------------------
c
      return
      end
c
c
c
      subroutine pr_eirene_analysis
      use error_handling
      implicit none
      include 'params'
      include 'outcom'
c
      include 'cgeom'
      include 'comtor'
c      include 'cneut2'
      include 'dynam2'
      include 'dynam3'
c      include 'dynam4'
      include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
      include 'printopt' 
c
c     OUT unstructured input
c
      include 'out_unstruc'
c
c     Local variables
c
c
      real*8 :: totsrc,totcore,totcore_poloidal(maxnks,3),
     >           totcore_radial(maxnrs,5)
      real*8 :: totcore_poloidal_vol(maxnks,3)
      real*8 :: totcore_radial_vol(maxnrs,5)
      real*8 :: totcore_poloidal_area(maxnks),
     >          totcore_radial_area(maxnrs)
      real*8 :: omp_reg,xpt_reg
      !real*8 :: psi1_reg,psi2_reg
      integer :: in

      integer :: nr
      real*8 :: plen, rc

c      integer prnizs
c      real    bgcontent,totbgcontent
c      real    impcontent(0:maxizs+1),totimpcontent(0:maxizs+1)
c      REAL    ZSUM(max(MAXPLRP,maxizs))
      real tmpsum,tmpsum2
c      real tote,toti,totn
c      real r,z
      integer ik,ir,iz
c      integer id,in 
c      
c     Calculate density profiles - first is outer midplane
c
      integer :: profcnt
      real*8 :: zprof,deltaz
      real*8 :: h_profiles(maxnrs,5)
      real*8 :: tmp_profiles(maxnrs,6)

      integer :: tu,ierr

c
c     Print out table of Eirene calcualted Hydrogenic values
c
      call find_free_unit_number(tu)
      open(tu,file='eirene_iz_analysis.dat',iostat=ierr,
     >     form='formatted')

      if (ierr.ne.0) then 
         call errmsg('PR_EIRENE_ANALYSIS: Error opening file',ierr)
         return
      endif


c
c     Core and pedestal analysis 
c
      write(tu,*)
      write(tu,'(a)') ' EIRENE Calculated Core Fueling'
      write(tu,*)
      write(tu,'(a,i6)')' SEPARATRIX RING                      = ',irsep 
      write(tu,'(a,4(1x,g18.8))') ' R0,Z0 :',r0,z0
      write(tu,'(a,4(1x,g18.8))') ' RXP,ZZP :',rxp,zxp

      
      ! integrate along the ring and print out the total ionization and scaled ionization

      ! First - sum total ionization source on grid and total in confined plasma 

      ! Use volumes for this calculation


      totsrc = 0.0
      totcore = 0.0
      totcore_poloidal = 0.0
      totcore_radial = 0.0
      totcore_poloidal_vol = 0.0
      totcore_radial_vol = 0.0
      totcore_poloidal_area =0.0
      totcore_radial_area = 0.0

      omp_reg = 0.05
      xpt_reg = 0.05
      !psi1_reg = 0.9
      !psi2_reg = 0.95

      write(tu,'(a,1x,g12.5)') ' XPT_REGION=Z<Z0_and_ABS(R-RXP)<= ',
     >                      xpt_reg
      write(tu,'(a,1x,g12.5)') ' OMP_REGION=R>R0_and_ABS(Z-Z0)<= ',
     >                      omp_reg
      write(tu,'(a,1x,g12.5)') ' NEAR_SEP_REGION1:_PSI_<=1.0'//
     >                         '_and_PSI_>',psi1_reg
      write(tu,'(a,1x,g12.5)') ' NEAR_SEP_REGION2:_PSI_<=1.0'//
     >                         '_and_PSI_>',psi2_reg


      do ir = 1,nrs
         do ik = 1,nks(ir)
            totsrc = totsrc + pinion(ik,ir) * kvol2(ik,ir)
            if (ir.lt.irsep.and.ik.ne.nks(ir)) then 
               totcore = totcore + pinion(ik,ir) * kvol2(ik,ir)
            endif
         end do 
      end do


      do ir = 2,irsep-1
         do ik = 1,nks(ir)-1
            totcore_poloidal(ik,1) = totcore_poloidal(ik,1) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
            totcore_radial(ir,1) = totcore_radial(ir,1) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
            totcore_poloidal_vol(ik,1) = totcore_poloidal_vol(ik,1) + 
     >                  kvol2(ik,ir)
            totcore_radial_vol(ir,1) = totcore_radial_vol(ir,1) + 
     >                  kvol2(ik,ir)

            ! look for poloidal subset regions - PSI > 0.9
            if (psitarg(ir,1).gt.0.9) then
                totcore_poloidal(ik,2) = totcore_poloidal(ik,2) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
                totcore_poloidal_vol(ik,2)=totcore_poloidal_vol(ik,2)+ 
     >                  kvol2(ik,ir)
            endif

            ! look for poloidal subset regions - PSI > 0.95
            if (psitarg(ir,1).gt.0.95) then
                totcore_poloidal(ik,3) = totcore_poloidal(ik,3) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
                totcore_poloidal_vol(ik,3)=totcore_poloidal_vol(ik,3)+ 
     >                  kvol2(ik,ir)
            endif


            ! Look for radial subset regions - near outer midplane

            ! Near outer midplane for radial
            if (abs(zs(ik,ir)-z0).lt.omp_reg.and.rs(ik,ir).gt.r0) then 
               totcore_radial(ir,2) = totcore_radial(ir,2) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
               totcore_radial_vol(ir,2) = totcore_radial_vol(ir,2) + 
     >                  kvol2(ik,ir)
            endif

            ! Near Xpoint for radial
            if (abs(rs(ik,ir)-rxp).lt.xpt_reg.and.
     >          abs(zs(ik,ir)-zxp).lt.abs(z0-zxp)
     >          ) then 
               totcore_radial(ir,3) = totcore_radial(ir,3) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
               totcore_radial_vol(ir,3) = totcore_radial_vol(ir,3) + 
     >                  kvol2(ik,ir)
            endif


            if (ik.lt.nks(ir)/2) then 
            ! Outer (JET Xpt up) for radial

               totcore_radial(ir,4) = totcore_radial(ir,4) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
               totcore_radial_vol(ir,4) = totcore_radial_vol(ir,4) + 
     >                  kvol2(ik,ir)

            else
            ! INNER (JET Xpt up) for radial

               totcore_radial(ir,5) = totcore_radial(ir,5) + 
     >                 pinion(ik,ir) * kvol2(ik,ir)
               totcore_radial_vol(ir,5) = totcore_radial_vol(ir,5) + 
     >                  kvol2(ik,ir)

            endif


            ! find the cell and take the length of the outer side
            ! which is between vertices 2 and 3

            nr = korpg(ik,ir)

            if (nr.ne.0) then 
c                 
               plen = sqrt((rvertp(2,nr)-rvertp(3,nr))**2+
     >                      (zvertp(2,nr)-zvertp(3,nr))**2)
c
               rc = (rvertp(2,nr)+rvertp(3,nr))/2.0
c
               ! calculate poloidal area at edge of confined plasma
               if (ir.eq.irsep-1) then 
                  totcore_poloidal_area(ik) = 2.0 * PI * rc * plen

               endif

               ! Calculate surface area of confined plasma at each radial location
               totcore_radial_area(ir) = totcore_radial_area(ir) 
     >                + 2.0 * PI * rc * plen
           endif


         end do
      end do

      ! Convert all the sources to density

      do ir = 2,irsep-1
         do in = 1,5
            if (totcore_radial_vol(ir,in).gt.0.0) then 
               totcore_radial(ir,in) = 
     >            totcore_radial(ir,in)/totcore_radial_vol(ir,in)
            else 
               totcore_radial(ir,in) = 0.0
            endif
         end do
      enddo

      do ik = 1,nks(irsep-1)
         do in = 1,3
            if (totcore_poloidal_vol(ik,in).gt.0.0) then 
               totcore_poloidal(ik,in) = 
     >            totcore_poloidal(ik,in)/totcore_poloidal_vol(ik,in)
            else 
               totcore_poloidal(ik,in) = 0.0
            endif
         end do
      enddo

      write(tu,'(a,g18.8)') ' Total Ionization     :',totsrc
      write(tu,'(a,g18.8)') ' Total Core_Ionization:',totcore
      if (totsrc.ne.0.) then 
         write(tu,'(a,g18.8)') ' Fraction Core_IZ  :',totcore/totsrc
      endif

      write(tu,*)
      write(tu,'(a)') ' NOTE: All integrated ionization'//
     >                ' profile data below is calculated in terms'//
     >                ' of density (ionization in region/'//
     >                'volume of region)'

      write(tu,*)

      write(tu,'(a)') ' Core poloidal'//
     >               ' ionization distributions'
      
      write(tu,*)

      write(tu,'(3x,a,3x,2(8x,a,8x),2(9x,a,9x),6x,a,5x,5x,a,5x
     >       2(4x,a,4x,4x,a,3x),2x,a,2x)') 
     >       'IK','KSS','KPS','R','Z','Total_IZ','Total_vol',
     >       'NearSEP1_IZ','NearSEP1_vol',
     >       'NearSEP2_IZ','NearSEP2_vol','Separatrix_area'
                    
      ir = irsep-1
      do ik = 1,nks(irsep-1)-1
         write(tu,'(i8,12(1x,g18.8))') ik,kss(ik,ir),
     >         kps(ik,ir),rs(ik,ir),zs(ik,ir),
     >         totcore_poloidal(ik,1),totcore_poloidal_vol(ik,1),
     >         totcore_poloidal(ik,2),totcore_poloidal_vol(ik,2),
     >         totcore_poloidal(ik,3),totcore_poloidal_vol(ik,3),
     >         totcore_poloidal_area(ik)
      end do
         
      write(tu,*)

      write(tu,'(a)') ' Core radial'//
     >               ' ionization distributions'
      write(tu,*)

      write(tu,'(3x,a,3x,5x,a,5x,6x,a,5x,5x,a,5x,
     >           2(7x,a,6x,6x,a,6x),
     >           2(6x,a,5x,5x,a,5x),
     >           4x,a,3x)')  
     >           'IR','PSIN','Total_IZ','Total-vol',
     >           'OMP_IZ','OMP-vol',
     >           'XPT_IZ','XPT-vol',
     >           OUTER//'_IZ ',OUTER//'-vol ',
     >           INNER//'_IZ ',INNER//'-vol ',
     >           'Surface_area'


      do ir = 2,irsep-1
         write(tu,'(i8,1x,f13.6,12(1x,g18.8))') ir,
     >         psitarg(ir,1),
     >         totcore_radial(ir,1),totcore_radial_vol(ir,1),
     >         totcore_radial(ir,2),totcore_radial_vol(ir,2),
     >         totcore_radial(ir,3),totcore_radial_vol(ir,3),
     >         totcore_radial(ir,4),totcore_radial_vol(ir,4),
     >         totcore_radial(ir,5),totcore_radial_vol(ir,5),
     >         totcore_radial_area(ir)
      end do



      ! Calculate outer midplane profiles
      ! 
      ! 

      deltaz = 0.05
      zprof = z0
      h_profiles = 0.0
      tmp_profiles = 0.0


      do ir = 1,nrs
         do ik = 1,nks(ir)
            if (rs(ik,ir).gt.r0.and.
     >         (zs(ik,ir).ge.zprof-deltaz).and.
     >         (zs(ik,ir).le.zprof+deltaz)) then

               if (pinatom(ik,ir).gt.0.0) then 
                  tmp_profiles(ir,1) = tmp_profiles(ir,1)+
     >                kvol2(ik,ir)
                  tmp_profiles(ir,2) = tmp_profiles(ir,2)+
     >                pinatom(ik,ir)*kvol2(ik,ir)
               endif

               if (pinmol(ik,ir).gt.0.0) then 
                  tmp_profiles(ir,3) = tmp_profiles(ir,3)+
     >                kvol2(ik,ir)
                  tmp_profiles(ir,4) = tmp_profiles(ir,4)+
     >                pinmol(ik,ir)*kvol2(ik,ir)
               endif


               if (pinion(ik,ir).gt.0.0) then 
                  tmp_profiles(ir,5) = tmp_profiles(ir,5)+
     >                kvol2(ik,ir)
                  tmp_profiles(ir,6) = tmp_profiles(ir,6)+
     >                pinion(ik,ir)*kvol2(ik,ir)
               endif

            endif
         end do
      end do

      ! Calculate average density
      profcnt = 0

      do ir = 1,nrs
         if (tmp_profiles(ir,1).gt.0.0) then 
            profcnt = profcnt + 1
            h_profiles(profcnt,1)=ir
            h_profiles(profcnt,2)=psitarg(ir,1)
            
            if (tmp_profiles(ir,1).gt.0.0) then 
               h_profiles(profcnt,3)=tmp_profiles(ir,2)
     >                              /tmp_profiles(ir,1)
            endif

            if (tmp_profiles(ir,3).gt.0.0) then 
               h_profiles(profcnt,4)=tmp_profiles(ir,4)
     >                              /tmp_profiles(ir,3)
            endif

            if (tmp_profiles(ir,5).gt.0.0) then 
               h_profiles(profcnt,5)=tmp_profiles(ir,6)
     >                              /tmp_profiles(ir,5)            
            endif
         endif
      end do

      ! write out the outer midplane density profiles

      write(tu,*)
      write(tu,'(a,f12.5,a,f12.5,a,f12.5)') 
     >                ' OUTER MID-PLANE DENSITY PROFILES: R > ',
     >                 r0,' : Z = ',z0,' +/- ',deltaz
      write(tu,'(3x,a,2x,9x,a,8x,8x,a,7x,5x,a,5x,5x,a,4x,7x,a,6x)') 
     >                'CNT','IR','PSIN',
     >                'D_Density','D2_Density',
     >                'Dioniz'

      do in = 1,profcnt
         write(tu,'(i8,12(1x,g18.8))') in,
     >         h_profiles(in,1),h_profiles(in,2),
     >         h_profiles(in,3),h_profiles(in,4),
     >         h_profiles(in,5)
      end do
c
c
c
      write(tu,*)
      write(tu,'(a)') ' EIRENE Calculated Hydrogenic Quantities'
      write(tu,*)
      write(tu,'(a,i6)')' SEPARATRIX RING                      = ',irsep 
      write(tu,'(a,i6)') ' PFZ RING ADJACENT TO SEPARATRIX RING = ',nrs 

      do ir = 1,nrs

         write(tu,*) 
c         write(tu,'(a,2x,a,i6)') ' START TARGET:'//OUTER,' RING = ',ir
         write(tu,'(2(3x,a,2x),6x,a,6x,4(9x,a,9x),
     >   2(8x,a,8x)  ,6x,a,6x,4(7x,a,6x)  )') 
     >   'IK','IR','R','Z','S','P',
     >   'DELTA-P','PSI','VOL',
     >   'PINATOM','PINENA','PINMOL','PINENM','PINION'

         tmpsum2 = 0.0

         do ik = 1,nks(ir)

            write(tu,'(2(1x,i5),1p,15(1x,g18.8))') ik,ir,
     >           rs(ik,ir),zs(ik,ir),kss(ik,ir),kps(ik,ir),
     >           kpb(ik,ir)-kpb(ik-1,ir),
     >           psifl(ik,ir),kvol2(ik,ir),
     >           pinatom(ik,ir),pinena(ik,ir),
     >           pinmol(ik,ir),pinenm(ik,ir),pinion(ik,ir)
c
c            tmpsum = 
c     >           sqrt((krb(ik,ir)-krb(ik-1,ir))**2+
c     >                (kzb(ik,ir)-kzb(ik-1,ir))**2)
c
c            tmpsum2= tmpsum2 + tmpsum
c
c
c            write(tu,'(2(1x,i5),1p,14(1x,g12.5))') ik,ir,
c     >           krb(ik,ir),kzb(ik,ir),
c     >           krb(ik-1,ir),kzb(ik-1,ir),rs(ik,ir),zs(ik,ir),
c     >           kpb(ik,ir)-kpb(ik-1,ir),tmpsum,
c     >           kpb(ik,ir), tmpsum2,ksb(ik,ir),ksb(ik-1,ir)
c

         end do  
c         write(tu,'(a)') 'END TARGET:'//INNER

      end do

      close(tu)

      ! end of pr_eirene_analysis
      return 
      end
c
c
c
      subroutine pr_exb_analysis
      use error_handling
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'driftvel'

      integer :: ik,ir,id,in,in_sep
      integer :: tu   ! temp unit number
      integer :: ierr
      character*1024 :: headings
      real :: efact,fact
      real :: exb_pol_tmp
c
c     jdemod - add local arrays for calculating fluxes
c

      real :: frad(maxnks,maxnrs),fpol(maxnks,maxnrs) 
      real :: fpol_int_sol(maxnks),fpol_int_pfz(maxnks)
      real :: frad_int_id(maxnrs), frad_int_od(maxnrs)
      integer :: ikref_id,ikref_od,ikref_pfz_id,ikref_pfz_od
c     
c     Print out table of ExB related values
c
      call find_free_unit_number(tu)
      open(tu,file='exb_source_terms.dat',iostat=ierr)

      if (ierr.ne.0) then 
         call errmsg('PR_ExB_ANALYSIS: Error opening file',ierr)
         return
      endif

      write(tu,*) 'ExB Analysis'

      headings = ' IK  IR   R   Z   S  S_Lower  S_Upper'//
     >           ' Cos_out cos_in  Dist_in  Dist_out  Tot_width'//
     >           '  Pol_length '//
     >           '  Btot  Btot/Bpol   FACT    Ne    Te     KES'//
     >           '  PHI    E-radial   E-poloidal   V-exb-radial'//
     >           '   V-exb-poloidal  V-exb-poloidal-Spara '//
     >           '  F-radial  F-poloidal '

c
c     Corrective scaling factor for KES 
c
      efact =  QTIM * QTIM * EMI / CRMI
c
      do ir = irsep, nrs
         write(tu,'(a)') trim(headings)
         do ik = 1,nks(ir)



            if (ik.eq.1) then
               id = idds(ir,2)
               ! print first target values - use ik = 0 
               write(tu,'(1x,2i8,30(1x,g18.8))')
     >          ik,ir,krb(ik-1,ir),kzb(ik-1,ir),0.0,
     >          ksb(ik-1,ir),ksb(ik,ir),
     >          0.0,0.0,
     >          0.0,0.0,0.0,
     >          0.0,          
     >          bts(ik,ir),kbfs(ik,ir),sqrt(kbfs(ik,ir)**2-1.0),
     >          knds(id),kteds(id),keds(id),
     >          osmpot2(ik-1,ir),0.0,0.0,0.0,0.0,0.0,
     >          0.0,0.0

            endif
            
            if ((kbfs(ik,ir)**2-1.0).lt.0.0) then 
                exb_pol_tmp = 0.0
            else
               fact = sqrt(kbfs(ik,ir)**2-1.0)
               if (fact.ne.0.0) then 
                  exb_pol_tmp = exb_pol_drft(ik,ir)/fact/qtim
               else
                   exb_pol_tmp = 0.0
               endif
            endif



            frad(ik,ir) =   knbs(ik,ir)*exb_rad_drft(ik,ir)/qtim*
     >           (kpb(ik,ir)-kpb(ik-1,ir)) * 2.0 * PI * rs(ik,ir)

            fpol(ik,ir) =   knbs(ik,ir) * exb_pol_tmp *
     >           (distin(ik,ir)+distout(ik,ir)) * 2.0 * PI * rs(ik,ir)


            write(tu,'(1x,2i8,30(1x,g18.8))')
     >          ik,ir,rs(ik,ir),zs(ik,ir),kss(ik,ir),
     >          ksb(ik-1,ir),ksb(ik,ir),
     >          cosalo(ik,ir),cosali(ik,ir),
     >          distin(ik,ir),distout(ik,ir),
     >          (distin(ik,ir)+distout(ik,ir)),
     >          (kpb(ik,ir)-kpb(ik-1,ir)),
     >          bts(ik,ir),kbfs(ik,ir),sqrt(kbfs(ik,ir)**2-1.0),
     >          knbs(ik,ir),ktebs(ik,ir),kes(ik,ir)/efact,
     >          osmpot2(ik,ir),e_rad(ik,ir),e_pol(ik,ir),
     >          exb_rad_drft(ik,ir)/qtim,exb_pol_tmp,
     >          exb_pol_drft(ik,ir)/qtim,
     >          frad(ik,ir), fpol(ik,ir)
            


            if (ik.eq.nks(ir)) then
               id = idds(ir,1)
               ! print target values - use nks(ir)+1
               write(tu,'(1x,2i8,30(1x,g18.8))')
     >          ik,ir,krb(ik,ir),kzb(ik,ir),ksmaxs(ir),
     >          ksb(ik-1,ir),ksb(ik,ir),
     >          0.0,0.0,
     >          0.0,0.0,0.0,
     >          0.0,          
     >          bts(ik,ir),kbfs(ik,ir),sqrt(kbfs(ik,ir)**2-1.0),
     >          knds(id),kteds(id),keds(id),
     >          osmpot2(ik+1,ir),0.0,0.0,0.0,0.0,0.0,
     >          0.0,0.0
            endif

         end do

      end do
c
c

c
c     Calculate integrated flux profiles
c
      fpol_int_sol = 0.0
      fpol_int_pfz = 0.0
      
      do ik = 1,nks(irsep)
         do ir = irsep,irwall-1
            fpol_int_sol(ik) = fpol_int_sol(ik)+ fpol(ik,ir)
         end do
      end do

      do ik = 1,nks(nrs)
         do ir = irtrap+1,nrs
            fpol_int_pfz(ik) = fpol_int_pfz(ik)+ fpol(ik,ir)
         end do
      end do

      frad_int_id = 0.0
      frad_int_od = 0.0
      
      ikref_id = ikouts(1,irsep-1) -1
      ikref_od = ikouts(nks(irsep-1)-1,irsep-1) + 1

      ikref_pfz_id = ikins(ikref_id,irsep)
      ikref_pfz_od = ikins(ikref_od,irsep)

c      write(0,*) 'EXB IKREFS:',ikref_id,ikref_od,
c     >             ikref_pfz_id,ikref_pfz_od

      do ir = irsep,irwall-1
         do ik = 1,ikref_id
            frad_int_id(ir) = frad_int_id(ir)+frad(ik,ir)
         end do 

         do ik = ikref_od,nks(ir)
            frad_int_od(ir) = frad_int_od(ir)+frad(ik,ir)
         end do 
      end do

      do ir = irtrap+1,nrs
         do ik = 1,ikref_pfz_id
            frad_int_id(ir) = frad_int_id(ir)+frad(ik,ir)
         end do 

         do ik = ikref_pfz_od,nks(ir)
            frad_int_od(ir) = frad_int_od(ir)+frad(ik,ir)
         end do 
      end do


c
c     Write out flux profiles
c
      write(tu,*) 
      write(tu,*) ' EXB_IK_REFERENCES:'
      write(tu,*) ' INNER_DIVERTOR_MAIN_SOL:',1,ikref_id
      write(tu,*) ' OUTER_DIVERTOR_MAIN_SOL:',ikref_od,nks(irsep)
      write(tu,*) ' INNER_DIVERTOR_PFZ:',1,ikref_pfz_id
      write(tu,*) ' OUTER_DIVERTOR_PFZ:',ikref_pfz_od,nks(nrs)
      write(tu,*) 
      write(tu,*) ' Integrated flux profiles'
      write(tu,*) 
      write(tu,*) ' Poloidal flux (main SOL)'
      write(tu,*) 
      write(tu,'(a)') ' IK  S(m)  P(m)  Rsep(m)  Zsep(m)'//
     >         ' Fpol(part/s)'
      ir = irsep
      do ik = 1,nks(ir)
         write(tu,'(1x,i8,20(1x,g18.8))') ik,kss(ik,ir),
     >        kps(ik,ir),rs(ik,ir),zs(ik,ir),fpol_int_sol(ik)
      end do

      write(tu,*) 
      write(tu,*) ' Poloidal flux (PFZ)'
      write(tu,*) 
      write(tu,'(a)') ' IK  S(m)  P(m)  Rsep(m)  Zsep(m)'//
     >         ' Fpol(part/s)'
      ir = nrs
      do ik = 1,nks(ir)
         write(tu,'(1x,i8,20(1x,g18.8))') ik,kss(ik,ir),
     >        kps(ik,ir),rs(ik,ir),zs(ik,ir),fpol_int_pfz(ik)
      end do

      write(tu,*) 
      write(tu,*) ' Radial flux (ID and OD including PFZ)'
      write(tu,*) 
      write(tu,'(a)') ' IN  IR Rinner(m) Zinner(m) Frad_inner(part/s)'//
     >         '    Router(m) Zouter(m) Frad_outer(part/s)'


      in = 1
      do ir = irwall-1,irsep,-1
         write(tu,'(1x,2i8,20(1x,g18.8))') in,ir,
     >         rp(idds(ir,2)),zp(idds(ir,2)),frad_int_id(ir),
     >         rp(idds(ir,1)),zp(idds(ir,1)),frad_int_od(ir)
         in = in+1
      end do

      in_sep = in-1

      do ir = nrs,irtrap+1,-1
         write(tu,'(1x,2i8,20(1x,g18.8))') in,ir,
     >         rp(idds(ir,2)),zp(idds(ir,2)),frad_int_id(ir),
     >         rp(idds(ir,1)),zp(idds(ir,1)),frad_int_od(ir)
         in = in+1
      end do

      write(tu,*) 
      write(tu,*) ' SEPARATRIX_INDEX:',in_sep


c
c     Close the file and free the unit number
c      
      close(tu)
      
      return
      end

      subroutine pr_imp_density_profiles
      implicit none
      include 'params'
      include 'outcom'
c     
      include 'cgeom'
      include 'comtor'
c     
      include 'dynam2'
      include 'dynam3'
c     
      include 'printopt' 


!     print two profiles - near outer midplane defined as Z = Z0 +/- 10cm with R>R0
!     - top of machine defined as +/- 10cm from the R,Z defined by the 
!     mid point of the separatrix ring.

      integer ir,in,ik,iz,ounit
      real omp_prof(maxnrs,3)
      real top_prof(maxnrs,3)
      real r_mid,z_mid,r_top,z_top
      real range
      character*1024 :: filename

      write(0,*) 'Printing impurity density profiles:'

      range = 0.1

      omp_prof = 0.0
      top_prof = 0.0

      r_mid = r0
      z_mid = z0

      r_top = rs(nks(irsep)/2,irsep)
      z_top = zs(nks(irsep)/2,irsep)
      
      
      do ir = 1, irwall
         do ik = 1,nks(ir)

!     omp
            if (rs(ik,ir).ge.r_mid.and.
     >           abs(zs(ik,ir)-z_mid).le.range) then 

               omp_prof(ir,1) = omp_prof(ir,1)+ kareas(ik,ir)

               do iz = 1,nizs

                  omp_prof(ir,2) = omp_prof(ir,2)+ 
     >                 kareas(ik,ir)* sdlims(ik,ir,iz)

               end do 

            endif

!     "top"
            if (((xpoint_up.and.zs(ik,ir).lt.z_mid).or.
     >           (.not.xpoint_up.and.zs(ik,ir).gt.z_mid)).and.
     >           abs(rs(ik,ir)-r_top).le.range) then 

            top_prof(ir,1) = top_prof(ir,1)+ kareas(ik,ir)

            do iz = 1,nizs

               top_prof(ir,2) = top_prof(ir,2)+ 
     >              kareas(ik,ir)* sdlims(ik,ir,iz)
            end do

         endif
         
      end do
      end do
      
!     calculate density

      do ir = 1,irwall

         if (omp_prof(ir,1).gt.0.0) then 
            omp_prof(ir,3) = omp_prof(ir,2)/omp_prof(ir,1)
         endif 


         if (top_prof(ir,1).gt.0.0) then 
            top_prof(ir,3) = top_prof(ir,2)/top_prof(ir,1)
         endif 

      end do 


!
!     convert to storing in a separate file
!
      filename='upstream_impurity_profiles.dat'
      call find_free_unit_number(ounit)

      open(ounit,file=trim(filename),form='formatted')

! write out the tabulated data

      write(ounit,'(a)')  ' Tabulated Impurity Density Profiles:'
      write(ounit,'(a)')  ' OMP'
      write(ounit,'(a)')  ' IR        PSIN      AREA'//
     >     '     CONTENT    DENSITY   DENSITY-ABS'

      do ir = 1,irwall

         if (omp_prof(ir,1).gt.0.0) then 
            write(ounit,'(1x,i8,5(1x,g18.8))') ir,psitarg(ir,2),
     >           (omp_prof(ir,in),in=1,3),omp_prof(ir,3)*absfac
         endif 
      end do



      write(ounit,'(a)')  ' TOP'
      write(ounit,'(a)')  ' IR        PSIN      AREA'//
     >     '     CONTENT    DENSITY   DENSITY-ABS'

      do ir = 1,irwall

         if (top_prof(ir,1).gt.0.0) then 
            write(ounit,'(1x,i8,5(1x,g18.8))') ir,psitarg(ir,2),
     >           (top_prof(ir,in),in=1,3),top_prof(ir,3)*absfac
         endif 
      end do


      close(ounit)

      end                       ! pr_imp_density_profiles


