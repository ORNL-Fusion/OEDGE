c     -*-Fortran-*-
c
c
      subroutine osmprobe(lvals,louts,osmvals,osmplots,
     >                 r1p,z1p,r2p,z2p,int_type,crmb,qtim)
      use mod_params
      use mod_cgeom
      implicit none
c
c     include 'params'
c     include 'cgeom'

c
      integer osmvals,osmplots,int_type,i1
      real lvals(maxseg,maxngs),louts(maxseg)
      real r1p,z1p,r2p,z2p,crmb,qtim
c
c     OSMPROBE: This routine extract the OSM results from the background
c               plasma along a specific line intended to model the location
c               of a probe - this line is defined from P(R1,Z1) to P(R2,Z2).
c
c               The INT_TYPE parameter will allow for future enhancement for
c               calculation of cuts through any part of the grid - for
c               now - intersections are only calculated for the main SOL
c               and only for the first intersection between the line
c               segment and the rings composing the grid.
c
c               **NOTE**: INT_TYPE now defines the type of axis to be used
c                         for the rcp/osm plot
c                         1 = mid-plane distance - Jsat - 0.35 multiplier
c                         2 = R-coordinate of intersections - Ne
c                         3 = Z-coordinate of intersections - Ne
c                         4 = mid-plane distance - Ne
c                         5 = mid-plane distance - Jsat - 0.5 multiplier
c                         6 = PSIN - Ne  
c
c
c               The probe values returned are found by interpolating
c               between the nearest grid points. This code assumes that
c               a properly defined polygonal grid is in use - it can
c               be modified to work for a cell center only grid in necessary
c               but since these are not in common use at this time it
c               is not a useful enhancement.
c
c               The extracted values are recorded relative to the outer
c               outer mid-plane coordinates.
C
C               osmvals - output index in lvals
C                     1 - Ne
c                     2 - Te
C                     3 - Ti
c                     4 - pressure
c                     5 - Vb
c
c
c     Local Variables
c
      integer in,ip,ik,ir,icnt
      real sint,psin,pint
      real tmp,rsect,zsect
c
c     Loop through the main SOL rings starting at the separatrix and the
c     outer target - find the first intersection with the specified LOS
c     and extract/interpolate the values of Ne, Te, Ti and pressure.
c
c     Counter for number of intersections found
c
      icnt = 0
c
c     jdemod - allow for collection of data from core as well - might help with plots
c
      do ir = 1,irwall-1
c
c      do ir = irsep,irwall-1
c
         do ik = 1,nks(ir)
c
            call find_intsect(ik,ir,r1p,z1p,r2p,z2p,rsect,zsect,
     >                        sint,pint,psin)
c
            if (sint.gt.0.0) then
c
c              Increment counter
c
               icnt = icnt + 1
c
c
c
               write(6,'(a,3i4,1p,8(1x,g12.5))') 'OSM1:',ik,ir,icnt,
     >               rsect,zsect,sint,kss(ik,ir),kss2(ik,ir)
              
c
c              Assign axis value depending on option
c
               if (int_type.eq.1.or.int_type.eq.4.or.int_type.eq.5) then

                  louts(icnt) = middist(ir,outer_targid)

               elseif (int_type.eq.2) then

                  louts(icnt) = rsect

               elseif (int_type.eq.3) then

                  louts(icnt) = zsect

               elseif (int_type.eq.6) then 

                  louts(icnt) = psitarg(ir,outer_targid)

               endif
c
c              Intersection greater than cell center
c
               if (sint.ge.kss(ik,ir)) then
c
c                 Grid point at target
c
                  if (ik.eq.nks(ir)) then
c
c                    Ne
c
                     lvals(icnt,1) = knbs(ik,ir)     +
     >                 (knds(idds(ir,1))-knbs(ik,ir) ) *
     >                 (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir))
c
c                    Te
c

                     lvals(icnt,2) = ktebs(ik,ir)     +
     >                 (kteds(idds(ir,1))-ktebs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir))
c
c                    Ti
c

                     lvals(icnt,3) = ktibs(ik,ir)     +
     >                 (ktids(idds(ir,1))-ktibs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir))
c
c                    Vb
c
                     lvals(icnt,5) = (kvhs(ik,ir)     +
     >                 (kvds(idds(ir,1))-kvhs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(ksmaxs(ir)-kss(ik,ir)))
     >                 /qtim
c
c                    Pressure
c
                     lvals(icnt,4) = lvals(icnt,1) * (
     >                         (lvals(icnt,2)+lvals(icnt,3)) *ech +
     >                         crmb * amu * lvals(icnt,5)**2)
c
c                 Grid points away from target
c
                  else
c
c                    Ne
c
                     lvals(icnt,1) = knbs(ik,ir)     +
     >                 (knbs(ik+1,ir)-knbs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
c
c                    Te
c

                     lvals(icnt,2) = ktebs(ik,ir)     +
     >                 (ktebs(ik+1,ir)-ktebs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
c
c                    Ti
c

                     lvals(icnt,3) = ktibs(ik,ir)     +
     >                 (ktibs(ik+1,ir)-ktibs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
c
c                    Vb
c
                     lvals(icnt,5) = (kvhs(ik,ir)     +
     >                 (kvhs(ik+1,ir)-kvhs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir)))
     >                 /qtim
c
c                    Pressure
c
                     lvals(icnt,4) = lvals(icnt,1) * (
     >                         (lvals(icnt,2)+lvals(icnt,3)) *ech +
     >                         crmb * amu * lvals(icnt,5)**2)
c
                  endif
c
c              Intersection less than cell center
c

               else
c
c                 Grid points at the target
c
                  if (ik.eq.1) then

c
c                    Ne
c
                     lvals(icnt,1) = knbs(ik,ir)     +
     >                 (knds(idds(ir,2))-knbs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir))
c
c                    Te
c

                     lvals(icnt,2) = ktebs(ik,ir)     +
     >                 (kteds(idds(ir,2))-ktebs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir))
c
c                    Ti
c

                     lvals(icnt,3) = ktibs(ik,ir)     +
     >                 (ktids(idds(ir,2))-ktibs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir))
c
c                    Vb
c
                     lvals(icnt,5) = (kvhs(ik,ir)     +
     >                 (kvds(idds(ir,2))-kvhs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)))
     >                 /qtim
c
c                    Pressure
c
                     lvals(icnt,4) = lvals(icnt,1) * (
     >                         (lvals(icnt,2)+lvals(icnt,3)) *ech +
     >                         crmb * amu * lvals(icnt,5)**2)
c
c                 Grid points away from the target
c
                  else

c
c                    Ne
c
                     lvals(icnt,1) = knbs(ik,ir)     +
     >                 (knbs(ik-1,ir)-knbs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
c
c                    Te
c

                     lvals(icnt,2) = ktebs(ik,ir)     +
     >                 (ktebs(ik-1,ir)-ktebs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
c
c                    Ti
c

                     lvals(icnt,3) = ktibs(ik,ir)     +
     >                 (ktibs(ik-1,ir)-ktibs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
c
c                    Vb
c
                     lvals(icnt,5) = (kvhs(ik,ir)     +
     >                 (kvhs(ik-1,ir)-kvhs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir)))
     >                 /qtim
c
c                    Pressure
c
                     lvals(icnt,4) = lvals(icnt,1) * (
     >                         (lvals(icnt,2)+lvals(icnt,3)) *ech +
     >                         crmb * amu * lvals(icnt,5)**2)
c
                  endif

               endif
c slmod begin - temp
c          WRITE(6,'(5X,A,3I4,2F8.3,1P,E10.2,0P,2F8.2,1P,3E10.2,0P)') 
c     .     '--> ',
c     .     ik,ir,icnt,louts(icnt),zsect,(lvals(icnt,i1),i1=1,5),sint
c slmod end
               goto 100
c
            endif

         end do
c
c
c        Proceed with next ring
c
100      continue

c
      end do
c
c     Wrap up processing
c
c     Set number of points on each plot
c
      osmvals = icnt
c
c     Set number of sets of data available
c
      osmplots= 5
c
c     Verify order of scaling on axis and reorder if necessary
c
      if (louts(1).gt.louts(osmvals) ) then
c
c        Reorder from lowest to highest
c
         do in = 1,osmvals
c
            tmp = louts(osmvals-in+1)
            louts(osmvals-in+1) = louts(in)
            louts(in) = tmp
c
            do ip = 1,osmplots
c
               tmp = lvals(osmvals-in+1,ip)
               lvals(osmvals-in+1,ip) = lvals(in,ip)
               lvals(in,ip) = tmp
c
            end do
c
         end do
c
      end if
c
c
c     Print results
c
      write (6,*) 'OSM probe results:',osmvals
c
      do in = 1,osmvals
c
         write(6,'(a,i4,6(1x,g15.6))') 'OSM:',in,louts(in),
     >        lvals(in,1),lvals(in,2),lvals(in,3),lvals(in,4),
     >        lvals(in,5)
c
      end do
c
      return
      end
c
c
c
      subroutine rcpprobe(tvals,touts,rcpvals,rcpplots,
     >                 exp_ds,exp_offset,exp_dataopt,exp_vcalcopt,
     >                 exp_tcalcopt,exp_param,
     >                 lvals,louts,osmvals,osmplots,
     >                 int_type,rizb,crmb,datatitle,r1p,z1p,r2p,z2p)
      use mod_params
      implicit none
c     include 'params'
      integer rcpvals,rcpplots,exp_ds,exp_dstype,exp_vcalcopt
      integer exp_tcalcopt,exp_dataopt,osmvals,osmplots,int_type
      real    exp_param,exp_offset,rizb,crmb
      real tvals(maxthe,maxngs),touts(maxthe)
      real lvals(maxseg,maxngs),louts(maxseg)
      real r1p,z1p,r2p,z2p
      character*(*) datatitle
c slmod begin
      REAL CalcPressure,GetCs
c slmod end
c
c     RCPPROBE: This routine reads in the experimental Te and Jsats or Ne
c               reported for a specific shot at a specific time from
c               the experimental data file associated with the case.
c               The data can be in one of several forms -
c
c                      Rmid   Jsat          Te
c                      Rmid   Jsat1  Jsat2  Te
c                      Rmid   Jsatu1 Jsatd1 Te1 ... 2 ... 3  ... 4
c                      R      Z      Ne     Te
c                      PSIN          Ne     Te
c
c               The second case reports both Jsats for a differential
c               probe which allows an estimate of the Mach number to
c               be extracted.
c
c               The fourth set is more typical of a Thomson result and is
c               thus treated somewhat differently. It can only be plotted as
c               function of R or Z - in addition some of the other options
c               may not apply.
c
c               The EXP_CALCOPT input specifies the method that will
c               be used to interpret the Jsat and TE measurements when
c               converting them to both Ne and pressure.
c
c               Some of the assumptions covered by this parameter are:
c               1) Value used for Vpara in the pressure calculation
c                  - V=0
c                  - Vosm
c                  - Vmach (from mach probe interpretation)
c               2) Value of Ti used in density calculation
c                  - TiRCP = TiOSM
c                  - TiRCP = EXP_PARAM * TeRCP
c                  - TiRCP = TiOSM/TeOSM * TeRCP
c
c               The code interprets the probe data and finds values
c               for Ne, Te, Ti and Pressure - which will subsequently be
c               plotted against the equivalent OSM values.
c
c               In addition, the R coordinate of the RCP data will
c               be shifted by the value of EXP_OFFSET which allows for
c               shifting of the data to overlay the OSM result.
c
c               NOTE: the results of the OSM extraction are rquired
c                     for some of the RCP calculation models.
c
c               If a calculation model is chosen that uses data from the
c               OSM and the data is not available for the position
c               specified on the RCP axis then, depending on the options
c               used - data for the nearest defined OSM ring may be used.
c
C               rcpvals - output index in tvals - same order as osm data
C                     1 - Ne
c                     2 - Te
C                     3 - Ti
c                     4 - pressure
c                     5 - Vb
c
c
c     Local Variables
c
      integer dataunit,maxcols
      parameter (dataunit=13,maxcols=8)
c
c      integer maxdatx,dataunit,maxcols
c      parameter (maxdatx=5000,dataunit=13,maxcols=8)
c
      integer in,ik,ip,axis_type
      integer num_expt,ncols,ndata
      real expt_axis(maxdatx)
      real expt_data(maxdatx,maxcols)
      real cs,osm_interpolate
      real jsato,jsati,teav
      external osm_interpolate
c
c     Two types of experimental data may be loaded - these are specified by
c     the int_type parameter which also specifies the axis type.
c
c      For INT_TYPE = 1 - the axis is calculated for midplane R-coordinates (rho)
c                       - probe data in terms of Jsat ...
c                   = 2 - probe R-coordinate
c                       - probe data in terms of Ne, Te
c                   = 3 - probe Z-coordinate
c                       - probe data in terms of Ne, Te
c                   = 4 - the axis is calculated for midplane R-coordinates (rho)
c                       - probe data in terms of Ne, Te
c                   = 5 - the axis is calculated for midplane R-coordinates (rho)
c                       - probe data in terms of Jsat ...
c                   = 6 - the axis is stored in the data file in terms of PSIN
c                       - probe data is in the form Ne, Te
c                       
c
c      The offset is applied to the scale as appropriate.
c
c      For INT_TYPE = 2 or 3 the experimental data is expected to be
c      in the form R, Z, Ne, Te and thus needs to be treated/loaded in a
c      somewhat different fashion. NCOLS will be 3 for data of this
c      type - but the first column of data is actually the Z coordinate.
c
c
c     Load experimental data - ncols is expected to be either 2, 3, 6 or 8
c     2 = Jsat + Te values
c     2 = Ne   + Te Values
c     3 = Jsat1 + Jsat2 + Te values
c slmod begin
c     3 = Te + ne + Jsat1/Jsat2 (for exp_vcalcopt = 4)
c slmod end
c     6 = Jsat1 + Te1 + Jsat2 + Te2 + Jsat3 + Te3
c     8 = Jsat1 + Te1 + Jsat2 + Te2 + Jsat3 + Te3 + Jsat4 + Te4
c
      write (6,*) 'Loading:',exp_ds
c
      call load_expt_data(dataunit,exp_ds,expt_axis,axis_type,
     >                    expt_data,maxcols,maxdatx,
     >                    num_expt,ncols,datatitle)
c
      write (6,*) 'Loaded:',num_expt,ncols,datatitle
c
c     Only convert if the axis data is not needed in terms of PSIN. 
c
c     If the axis_type of the experimental data is in terms of PSIN - then 
c     convert this to Z along the probe axis by interpolating and extrapolating
c     using the PSI values for each ring on the grid along the probe line
c     as a reference. 
c
c     AXIS_TYPE = 6 = PSIN
c
      if (axis_type.eq.6.and.int_type.ne.6) then 
c
c        Include exp_offset in the shifted axis
c
         call convert_axis(expt_axis,num_expt,exp_offset,
     >                     r1p,z1p,r2p,z2p)
c
      endif
c
c     Apply the specified axis offset as the first action
c     depending on type of axis in use.
c
c     R-midplane or PSIN
c
      if (int_type.eq.1.or.int_type.eq.4.or.int_type.eq.5.or.
     >    int_type.eq.6) then
c
         do in = 1,num_expt
            expt_axis(in) = expt_axis(in) + exp_offset
         end do
c
c     R-probe - remove Z-coordinate, shift data
c
      elseif (int_type.eq.2) then
c
c        Reduce number of columns in data
c
         ncols = ncols-1
c
         do in = 1,num_expt
c
c           Apply offset to R-data
c
            expt_axis(in) = expt_axis(in) + exp_offset
c
c           Remove Z-data
c
            do ip = 1,ncols

               expt_data(in,ip) = expt_data(in,ip+1)

            end do
c
         end do
c
c
c     Z-probe - remove R-coordinate, shift data
c
      elseif (int_type.eq.3.and.axis_type.ne.6) then
c
c        Reduce number of columns in data
c
         ncols = ncols-1
c
         do in = 1,num_expt
c
c           Apply offset to Z-data
c
            expt_axis(in) = expt_data(in,1) + exp_offset
c
c           Remove Z-data
c
            do ip = 1,ncols

               expt_data(in,ip) = expt_data(in,ip+1)

            end do
c
         end do
c
      endif
c
c
c     Depending on the number of columns of data and the
c     RCP data processing option that has been specified -
c     Adjust the results accordingly.
c
c
c     If the data has been loaded in a raw format it must be processed
c     first.
c
c      write (6,*) 'RCP EXPT DATA:',num_expt,ncols
c      do ik = 1,num_expt
c
c          write(6,'(a,i4,8(1x,g12.5))') 'EXPT:',ik,
c     >                  (expt_data(ik,in),in = 1,ncols)
c
c
c      end do
c

      if (ncols.eq.8) then
c
         if (exp_dataopt.eq.1.or.
     >      (exp_dataopt.eq.2.and.(8*num_expt).ge.maxthe)) then
c
c           Average both outer and inner probe results to
c           get one value each for Jsato, Jsati.
c
c           Average all four values of Te to get a Te value.
c
            ncols = 3
c
            do ik = 1,num_expt
c
               jsato = (expt_data(ik,1) + expt_data(ik,3))/2.0
               jsati = (expt_data(ik,5) + expt_data(ik,7))/2.0
               teav  = (expt_data(ik,2) + expt_data(ik,4) +
     >                  expt_data(ik,6) + expt_data(ik,8)) / 4.0
c
               expt_data(ik,1) = jsato
               expt_data(ik,2) = jsati
               expt_data(ik,3) = teav
c

c
            end do
c
         elseif (exp_dataopt.eq.2) then
c
c           Try to generate a plot demonstrating the possible scatter
c           in the data - not recommended for analysis purposes.
c
c           This will result in 8 times the original data.
c
            ncols = 3
c
            do ik = 1,num_expt
c
c              Map to 8 sets of data using all the combinations
c              of jsato to jsati and Te.
c
               expt_data(ik+7*num_expt,1) = expt_data(ik,3)
               expt_data(ik+7*num_expt,2) = expt_data(ik,7)
               expt_data(ik+7*num_expt,3) = expt_data(ik,8)
               expt_axis(ik+7*num_expt)   = expt_axis(ik)
c
               expt_data(ik+6*num_expt,1) = expt_data(ik,3)
               expt_data(ik+6*num_expt,2) = expt_data(ik,7)
               expt_data(ik+6*num_expt,3) = expt_data(ik,4)
               expt_axis(ik+6*num_expt)   = expt_axis(ik)
c
               expt_data(ik+5*num_expt,1) = expt_data(ik,3)
               expt_data(ik+5*num_expt,2) = expt_data(ik,5)
               expt_data(ik+5*num_expt,3) = expt_data(ik,6)
               expt_axis(ik+5*num_expt)   = expt_axis(ik)
c
               expt_data(ik+4*num_expt,1) = expt_data(ik,3)
               expt_data(ik+4*num_expt,2) = expt_data(ik,5)
               expt_data(ik+4*num_expt,3) = expt_data(ik,4)
               expt_axis(ik+4*num_expt)   = expt_axis(ik)
c
               expt_data(ik+3*num_expt,1) = expt_data(ik,1)
               expt_data(ik+3*num_expt,2) = expt_data(ik,7)
               expt_data(ik+3*num_expt,3) = expt_data(ik,8)
               expt_axis(ik+3*num_expt)   = expt_axis(ik)
c
               expt_data(ik+2*num_expt,1) = expt_data(ik,1)
               expt_data(ik+2*num_expt,2) = expt_data(ik,7)
               expt_data(ik+2*num_expt,3) = expt_data(ik,2)
               expt_axis(ik+2*num_expt)   = expt_axis(ik)
c
               expt_data(ik+1*num_expt,1) = expt_data(ik,1)
               expt_data(ik+1*num_expt,2) = expt_data(ik,5)
               expt_data(ik+1*num_expt,3) = expt_data(ik,6)
               expt_axis(ik+1*num_expt)   = expt_axis(ik)
c
               expt_data(ik+0*num_expt,1) = expt_data(ik,1)
               expt_data(ik+0*num_expt,3) = expt_data(ik,2)
               expt_data(ik+0*num_expt,2) = expt_data(ik,5)
c
c              expt_axis(ik+0*num_expt)   = expt_axis(ik)
c
            end do
c
c           Reset value of num_expt so it includes all data.
c
            num_expt = 8 * num_expt
c
         endif
c
         write (6,*) 'RCP PROC DATA:',num_expt,ncols
c
         do ik = 1,num_expt

            write(6,'(a,i4,9(1x,g11.4))') 'EXPT:',ik,expt_axis(ik),
     >               (expt_data(ik,in),in = 1,ncols)

         end do
c
      elseif (ncols.eq.6) then
c
         if (exp_dataopt.eq.1.or.exp_dataopt.eq.2) then
c
c           Average all probe results to
c           get one value each for Jsat and Te.
c
            ncols = 2
c
            do ik = 1,num_expt
c
               jsato = (expt_data(ik,1) + expt_data(ik,3) +
     >                  expt_data(ik,5))/3.0
               teav  = (expt_data(ik,2) + expt_data(ik,4) +
     >                  expt_data(ik,6))/3.0
c
               expt_data(ik,1) = jsato
               expt_data(ik,2) = teav
c
            end do
c
         endif
c
      endif
c
      ndata = ncols
c
      write (6,*) 'NDATA:',ndata
c
c     Assign the number of points in the RCP data
c
      rcpvals = num_expt
      rcpplots= 5
c
c     Assign the axes for ALL of the RCP data
c
      do ik = 1,rcpvals
            touts(ik) = expt_axis(ik)
      end do
c
c     If V calculation is specified as MACH probe and the differential data
c     is not present - reset the calculation option to use OSM Vpara.
c
      if (ndata.eq.2.and.exp_vcalcopt.eq.3) exp_vcalcopt = 2
c
c
c     Use the tvals array as working space for calculating the various
c     quantities.
c
c     Assign the Electron temperature from the probe data directly
c     to the results arrays.
c
c     Next - calculate the <Jsat> value and store it temporarily in the
c     density index of the output array.
c
c
c     Te
c     <Jsat>
c     Ti
c     Vb
c     Pressure
c
      do in = 1,rcpvals
c
c
c        Te
c
c
         tvals(in,2) = expt_data(in,ndata)
c
c
c        Jsat or Ne
c
c
         if (ndata.eq.2) then
c
            tvals(in,1) = expt_data(in,1)
c
         elseif (ndata.eq.3) then
c
            tvals(in,1) = (expt_data(in,1)+expt_data(in,2))/2.0
c
         endif
c
c
c        Ti
c
c
c        Fixed ratio
c
         if (exp_tcalcopt.eq.1) then
c
            tvals(in,3) = tvals(in,2) * exp_param
c
c        Use Nearest OSM Ti/Te ratio
c
         elseif (exp_tcalcopt.eq.2) then
c
            tvals(in,3) = tvals(in,2) *
     >                  osm_interpolate(lvals(1,3),lvals(1,2),
     >                    louts,osmvals,touts(in),exp_param,2)
c
c        TiRCP = OSM Ti -> OSM Ti/Te outside range
c
         elseif (exp_tcalcopt.eq.3) then
c
            tvals(in,3) = osm_interpolate(lvals(1,3),lvals(1,2),
     >                    louts,osmvals,touts(in),tvals(in,2),3)
c
         endif
c
c
c        Convert Jsat to Ne if necessary
c
c
         if (int_type.eq.1) then

            cs = SQRT(EMI*(tvals(in,2)+tvals(in,3))/CRMB)
c
c           0.35 comes from interpretation of reciprocating type probes.
c
            tvals(in,1) = tvals(in,1) / (0.35 * ech * cs)
c
         elseif (int_type.eq.5) then

            cs = SQRT(EMI*(tvals(in,2)+tvals(in,3))/CRMB)
c
c           0.5 comes from an alternate interpretation of reciprocating type probes.
c
            tvals(in,1) = tvals(in,1) / (0.5 * ech * cs)
c
         endif
c
c
c        Vb
c
c
c        Vpara = zero
c
         if (exp_vcalcopt.eq.1) then
c
            tvals(in,5) = 0.0
c
c        Vpara = Vosm
c
         elseif (exp_vcalcopt.eq.2) then
c
            tvals(in,5) = osm_interpolate(lvals(1,5),lvals(1,5),
     >                     louts,osmvals,touts(in),0.0,1)
c
c        Vpara = Vmachprobe
c
         elseif (exp_vcalcopt.eq.3) then
c
            tvals(in,5) = cs*0.4*log(expt_data(in,1)/expt_data(in,2))
c
         endif
c
c
c        Pressure
c
c
         tvals(in,4) = tvals(in,1) * (
     >                         (tvals(in,2)+tvals(in,3)) *ech +
     >                         crmb * amu * tvals(in,5)**2)
c
c slmod begin
         IF (exp_vcalcopt.EQ.4) THEN
c...       Ne:
           tvals(in,1) = expt_data(in,1)
c...       Te,Ti:
           tvals(in,2) = expt_data(in,2)
           tvals(in,3) = expt_data(in,2)
c...       Vb:
           cs = GetCs(tvals(in,2),tvals(in,3))
           tvals(in,5) = cs * 0.43 * LOG(expt_data(in,3))
c...       p:
           tvals(in,4) = CalcPressure(tvals(in,1),tvals(in,2),
     .                                tvals(in,3),tvals(in,5)) * ECH
         ENDIF
c slmod end
      end do
c
c     Print results
c
      write (6,*) 'RCP probe results:',rcpvals
c
      do in = 1,rcpvals
c
         write(6,'(a,i4,6(1x,g15.6))') 'RCP:',in,touts(in),
     >        tvals(in,1),tvals(in,2),tvals(in,3),tvals(in,4),
     >        tvals(in,5)
c
      end do
c

      return
      end
c
c
c
      subroutine find_intsect(ik,ir,r1p,z1p,r2p,z2p,rsect,zsect,
     >                        sint,pint,psin)
      use mod_params
      use mod_cgeom
      implicit none
      integer ik,ir
      real r1p,z1p,r2p,z2p,sint
      real rsect,zsect,psin,pint
c
c     include 'params'
c     include 'cgeom'
c
c     FIND_INTSECT: This routine checks to see if the line segment
c                   defined by the along ring axis of the given ik,ir
c                   cell intersects with the line segment defined
c                   by the points (r1p,z1p),(r2p,z2p) within the cell.
c
c                   If this occurs it then calculates the S-value
c                   of the intersection.
c
c     Local Variables
c
      logical check_intsect
      external check_intsect
c
c     Check to see if the line segments intersect and return an
c     intersection if one is available.
c
c     Initialize rsect and zsect to large negative values
c
      rsect = -HI
      zsect = -HI
c
      sint = 0.0
      psin = 0.0
c
      if (check_intsect(r1p,z1p,r2p,z2p,krb(ik-1,ir),kzb(ik-1,ir),
     >                  krb(ik,ir),kzb(ik,ir),rsect,zsect)) then
c
         sint = ((rsect-krb(ik-1,ir))**2+(zsect-kzb(ik-1,ir))**2) /
     >   ((krb(ik,ir)-krb(ik-1,ir))**2+(kzb(ik,ir)-kzb(ik-1,ir))**2)
     >     * (ksb(ik,ir)-ksb(ik-1,ir)) + ksb(ik-1,ir)
c
         pint = ((rsect-krb(ik-1,ir))**2+(zsect-kzb(ik-1,ir))**2) /
     >   ((krb(ik,ir)-krb(ik-1,ir))**2+(kzb(ik,ir)-kzb(ik-1,ir))**2)
     >     * (kpb(ik,ir)-kpb(ik-1,ir)) + kpb(ik-1,ir)
c
         if (ik.gt.nks(ir)/2) then 
            psin = psitarg(ir,1)
         else
            psin = psitarg(ir,2) 
         endif 
c
c         write(6,'(a,2i4,1p,9(1x,g12.5))') 'SINT:',ik,ir,sint,
c     >         krb(ik,ir),rsect,krb(ik-1,ir),kzb(ik,ir),zsect,
c     >         kzb(ik-1,ir),
c     >         ksb(ik,ir),kss(ik,ir),ksb(ik-1,ir),kss(ik-1,ir),
c     >         kss(ik,ir),
c     >         kss(ik+1,ir)
c
c
      endif
c
      return
c
      end
c
c
c
      logical function check_intsect(r1p,z1p,r2p,z2p,
     >                               r3p,z3p,r4p,z4p,
     >                               rsect,zsect)
      implicit none
      real r1p,z1p,r2p,z2p
      real r3p,z3p,r4p,z4p
      real rsect,zsect
c
c     CHECK_INTSECT: This function takes two sets of points
c                    (P1,P2) and (P3,P4) which define two
c                    line segments - it then calculates their
c                    intersection point - if any - the result
c                    is set to .true. if the intersection point
c                    lies within the ends of both line segments.
c
c     Local Variables:
c
      logical linxi,linxw
      double precision mi,bi,mw,bw,rsect2
c
      check_intsect = .false.
      rsect = 0.0
      zsect = 0.0
c
c     Calculate slope and intersection of first line segment
c     Set LINXI for vertical lines
C
      IF (R1P.NE.R2P) THEN
         MI = DBLE(Z1P-Z2P)/DBLE(R1P-R2P)
         BI = DBLE(Z1P) - DBLE(R1P) * DBLE(MI)
         LINXI = .FALSE.
      ELSE
         MI = 0.0
         BI = DBLE(R1P)
         LINXI = .TRUE.
      ENDIF
c
c     Calculate slope and intersection of second line segment
c     Set LINXW for vertical lines
C
      IF (R3P.NE.R4P) THEN
         MW = DBLE(Z3P-Z4P)/DBLE(R3P-R4P)
         BW = DBLE(Z3P) - DBLE(R3P) * DBLE(MW)
         LINXW = .FALSE.
      ELSE
         MW = 0.0
         BW = R3P
         LINXW = .TRUE.
      ENDIF
c
c      write(6,'(a,4(1x,g12.5))') 'CheckA:',r1p,z1p,r2p,z2p
c      write(6,'(a,4(1x,g12.5))') 'CheckB:',r3p,z3p,r4p,z4p
c      write(6,'(a,2(2g12.5,l5))') 'Check1:',mi,bi,linxi,mw,bw,linxw
c
c     Calculate intersections
c
c
c     If line segment 2 is not a vertical line
c
      IF (.not.linxw) THEN
c
c        If first line segment is a vertical line
c
         IF (LINXI) THEN
c
            rsect = bi
            zsect = MW * BI + BW
c
            IF ((ABS(zsect-Z4P).LE.ABS(Z3P-Z4P)) .AND.
     >          (ABS(zsect-Z3P).LE.ABS(Z3P-Z4P)) .AND.
     >          (ABS(zsect-Z1P).LE.ABS(Z1P-Z2P)) .AND.
     >          (ABS(zsect-Z2P).LE.ABS(Z1P-Z2P)) .AND.
     >          (ABS(rsect-r4P).LE.ABS(r3P-r4P)) .AND.
     >          (ABS(rsect-r3P).LE.ABS(r3P-r4P)) .AND.
     >          (ABS(rsect-r1P).LE.ABS(r1P-r2P)) .AND.
     >          (ABS(rsect-r2P).LE.ABS(r1P-r2P))) THEN
c
               check_intSECT = .TRUE.
            ELSE
               check_intSECT = .FALSE.
            ENDIF
c
c        If first line segment is not a vertical line
c
         ELSEIF (.NOT.LINXI) THEN
c
            IF (MW.EQ.MI) THEN
c
               rsect = 0.0
               zsect = 0.0
c
               check_intSECT = .FALSE.
c
            ELSE
c
c...HAD TO INCREASE PRECISION FOR HORIZONTAL PROBE - SL
               Rsect2 = (BI - BW) / (MW - MI)
               Rsect  = SNGL(Rsect2)
               Zsect = MW * Rsect2 + BW
c
               IF ((ABS(zsect-Z4P).LE.ABS(Z3P-Z4P)) .AND.
     >             (ABS(zsect-Z3P).LE.ABS(Z3P-Z4P)) .AND.
     >             (ABS(zsect-Z1P).LE.ABS(Z1P-Z2P)) .AND.
     >             (ABS(zsect-Z2P).LE.ABS(Z1P-Z2P)) .AND.
     >             (ABS(rsect-r4P).LE.ABS(r3P-r4P)) .AND.
     >             (ABS(rsect-r3P).LE.ABS(r3P-r4P)) .AND.
     >             (ABS(rsect-r1P).LE.ABS(r1P-r2P)) .AND.
     >             (ABS(rsect-r2P).LE.ABS(r1P-r2P))) THEN
c
                  check_intSECT = .TRUE.
               ELSE
                  check_intSECT = .FALSE.
               ENDIF
c
            ENDIF
         ENDIF
c
c     Second line segment is vertical
c
      ELSEif (linxw) then
c
c        First segment is also vertical
c
         IF (LINXI) THEN
c
c           If the two segments represent the same vertical line -
c           return one of the points as an intersection and set
c           check_intsect to false since it is not a proper intersection.
c
            IF (BI.EQ.BW) THEN
c
               Rsect = r3p
               Zsect = z3p

            ELSE
c
               rsect = 0.0
               zsect = 0.0

            ENDIF
c
            check_intsect = .false.
c
c        First segment is not vertical
c
         ELSEIF (.NOT.LINXI) THEN
c
            rsect = bw
            zsect = MI * BW + BI
c
            IF ((ABS(zsect-Z4P).LE.ABS(Z3P-Z4P)) .AND.
     >          (ABS(zsect-Z3P).LE.ABS(Z3P-Z4P)) .AND.
     >          (ABS(zsect-Z1P).LE.ABS(Z1P-Z2P)) .AND.
     >          (ABS(zsect-Z2P).LE.ABS(Z1P-Z2P)) .AND.
     >          (ABS(rsect-r4P).LE.ABS(r3P-r4P)) .AND.
     >          (ABS(rsect-r3P).LE.ABS(r3P-r4P)) .AND.
     >          (ABS(rsect-r1P).LE.ABS(r1P-r2P)) .AND.
     >          (ABS(rsect-r2P).LE.ABS(r1P-r2P))) THEN
c
               check_intSECT = .TRUE.
            ELSE
               check_intSECT = .FALSE.
            ENDIF
         ENDIF
      ENDIF
c
c
c      write(6,*) 'Check2:',rsect,zsect,check_intsect
c
      return
      end
c
c
c
      subroutine convert_axis(expt_axis,num_expt,exp_offset,
     >                        r1p,z1p,r2p,z2p)
      use mod_params
      use mod_cgeom
      implicit none
c
c     include 'params'
c     include 'cgeom'
c
      integer num_expt
      real expt_axis(num_expt),exp_offset
      real r1p,z1p,r2p,z2p 

c
c     Local variables
c 
      real psiaxis(maxnrs),psimin,psimax
      real zvals(maxnrs)
      integer icnt,ncnt,in,ic,ir,ik
      real psin,sint,rsect,zsect,pint
      real axisp,axisz
c
      integer ipos
      external ipos
c
c     Find the Z,Psi axis coordinates - within main SOL region
c
c     Loop through the main SOL rings starting at the separatrix and the
c     outer target - find the first intersection with the specified LOS
c     and extract the Z and PSI intersection points.
c
c     Counter for number of intersections found
c
      icnt = 0
      ic   = 1 
c
      do ir = irsep,irwall-1
c
         do ik = 1,nks(ir)
c
            call find_intsect(ik,ir,r1p,z1p,r2p,z2p,rsect,zsect,
     >                        sint,pint,psin)
c
            if (sint.gt.0.0) then
c
c              Increment counter
c
               icnt = icnt + 1
c
               write(6,'(a,3i4,1p,8(1x,g12.5))') 'CVRT_AXIS1:',icnt,
     >               ir,ik,rsect,zsect,sint,psin
c
               psiaxis(icnt) = psin
               zvals(icnt) = zsect 
c
               goto 100
c
            end if

         end do

100      continue      

      end do 
c
      ncnt = icnt
c
      psimin = psiaxis(1)  
      psimax = psiaxis(ncnt)
c 
      do in = 1,num_expt
c
         axisp = expt_axis(in)
c
         if (axisp.le.psimin) then 
c
            ic = 0
            axisz = zvals(1) 
     >            - (zvals(2)-zvals(1))/(psiaxis(2)-psiaxis(1))
     >                               * (psimin-axisp)
c
         elseif (axisp.ge.psimax) then   

            ic = ncnt+1
            axisz = zvals(ncnt) 
     >            + (zvals(ncnt)
     >            -zvals(ncnt-1))/(psiaxis(ncnt)-psiaxis(ncnt-1))
     >                  * (axisp-psimax)
c
         else
c
            ic = ipos(axisp,psiaxis,ncnt)

            axisz = zvals(ic) 
     >            - (zvals(ic)-zvals(ic-1))/(psiaxis(ic)-psiaxis(ic-1)) 
     >                 * (psiaxis(ic)-axisp)

         endif  

         expt_axis(in) = axisz  + exp_offset
c
         if (ic.gt.1.and.ic.le.ncnt) then 

            write(6,'(a,2i5,8(1x,g12.5))') 'CVRT_AXIS2A:',in,ic,
     >                    axisp,axisz,expt_axis(in),
     >                    zvals(ic),zvals(ic-1),exp_offset

         else

            write(6,'(a,2i5,8(1x,g12.5))') 'CVRT_AXIS2B:',in,ic,
     >                    axisp,axisz,expt_axis(in),
     >                    psimin,psimax,exp_offset

         endif
c
      end do


      return 
      end
c
c
c
   


