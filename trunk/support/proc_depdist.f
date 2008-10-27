      program proc_depdist
      implicit none
c     
c     PROC_DEPDIST
c     
c     READS in the wall distribution and totals the depsotion in
c     different regions as a fraction of the total deposition
c     
c     This routine reads in deposition data files and produces an 
c     Excel compatible plot file. It reads in a list of the files and 
c     comments from unit 5 and writes the output to unit 6. 
c     
c     The specification file has the following format:
c     
c     line 1:   <number of data sets to process>
c     line 2..line N: "filename" "comment1" "comment2"
c     
c     Options: 
c     
c     The code takes the options -i -o or -w which specify 
c     Inner, outer or wall deposition. The default for no flag is 
c     wall. 
c     
c     The plot axis is specified by the -a option
c     -a 1  = R coordinate 
c     2  = Z coordinate
c     3  = distance from strike point for targets
c     distance along wall for wall
c     4  = Index number
c     5  = Wall length 
c     
c     -h = help - lists available options
c     
c     
c     
c     
c     
c     
c     
c     The code assumes that the geometry for all cases is the same
c     unless the -g option is set indicating that the geometry information
c     should be included for each data set. 
c     
      integer datunit,maxsets
      parameter(datunit=50,maxsets=40)

c
      real HI
      parameter(HI=1.0e37)
c
      integer axis_opt,include_geometry,output_region,sepnum
      integer dat_opt
      real scalef
      real axis_shift
      character*10 axis
      logical data_allocated,axis_loaded,norm_data
      character*100 filenames(maxsets),com1(maxsets),com2(maxsets)

      real,allocatable :: outaxis(:,:),outdata(:,:),dep_summary(:,:)
c     
      integer len,lenstr
      external lenstr
c     
      integer flag 
      integer nregion
c     
      integer ndata,nsets,naxis,is,in,id,ierr
      integer get_data_elements
      external get_data_elements
      character*100 formd,formt
c
c     Start and end coordinates for output
c
      real outstart,outend,icount,minaxis
      integer minaxis_in,startin,stepin
c     
c     Option processing
c     
      integer numopt,maxopts
      parameter(maxopts=100) 
      character*100 opts(maxopts)
c     
c     Initialization
c     
c     Axis is set to the distance along surface as the default
c     
c     Include geometry for all data sets is off by default   
c     
c     The default output region is set to the wall 
c     
c     The default axis option is now the distance along the
c     wall from the inner midplane
c     
      axis_opt=5
      axis_shift = 0.0
      include_geometry=0 
      output_region = 3
      data_allocated=.false.
      norm_data=.false.
      sepnum=31
      scalef=1.0
      dat_opt=0
c
c     Set output limits to the maximum range
c
      outstart = -HI
      outend = HI
c     
c     Load option values
c     
      call load_options(numopt,opts,maxopts)
c     
c     Parse the options and assign the values
c     
c     Specify desired options and values
c     
      do in = 1, numopt

c
c         write(0,'(a,a,a)') ':',opts(in)(1:3),':'
c
c     -a: Axis Option (1 to 5)
c     
c     Set axis option
c     
         if (opts(in)(1:2).eq.'-a'.and.opts(in)(1:3).ne.'-as') then
c     
            if (in.lt.numopt) then 
               read(opts(in+1),*,err=200) axis_opt
            else
               axis_opt=3
            endif

            write(0,'(a,i6)') 'Using Axis Option:',axis_opt
c     
c     -i: Select inner target region
c     
         elseif (opts(in)(1:2).eq.'-i') then
c     
            output_region=2
c     
c     -o: Select outer target region
c     
         elseif (opts(in)(1:2).eq.'-o'.and.
     >           opts(in)(1:3).ne.'-os'.and.
     >           opts(in)(1:3).ne.'-oe') then
c     
            output_region=1
c
c     -os: Load output start
c
         elseif (opts(in)(1:3).eq.'-os') then 
            
c            write(0,*) 'OUTSTART:',trim(opts(in+1))
            read(opts(in+1),*,err=200) outstart

c
c     -oe: Load output end
c
         elseif (opts(in)(1:3).eq.'-oe') then

c            write(0,*) 'OUTEND:',trim(opts(in+1))
            read(opts(in+1),*,err=200) outend

c     
c     -w: Select wall region
c
         elseif (opts(in)(1:2).eq.'-w') then
c     
            output_region = 3 
c     
c     -g: include geometry option - includes scale axis for each different data set
c     
         elseif (opts(in)(1:2).eq.'-g') then
c
            include_geometry = 1
c     
c     -d: Output data option
c     
         elseif (opts(in)(1:2).eq.'-d') then
c     
            if (in.lt.numopt) then 
               read(opts(in+1),*,err=200) dat_opt
               output_region=3 
            else
               dat_opt=0
            endif
c     
c     Set separatrix number
c     
c     -ns: Set the number of the separatrix ring
c     
         elseif (opts(in)(1:3).eq.'-ns') then
c
            if (in.lt.numopt) then 
               read(opts(in+1),*,err=200) sepnum
            else
               sepnum=31
            endif
c     
c     -s: Set scaling factor
c     
         elseif (opts(in)(1:2).eq.'-s') then
c     
            if (in.lt.numopt) then 
               read(opts(in+1),*,err=200) scalef
            else
               scalef=1.0
            endif
c     
            write(0,*) 'Scaling Factor Applied = ',scalef

c     
c     -as: Read in axis shift - moves the axis by this constant offset
c     
         elseif (opts(in)(1:3).eq.'-as') then
c     
            if (in.lt.numopt) then 
               read(opts(in+1),*,err=200) axis_shift
            else
               axis_shift = 0 
            endif
c     
            write(0,*) 'Axis Shift Applied = ',axis_shift
c     
c     Print normalized data sets
c     
         elseif (opts(in)(1:2).eq.'-n') then
c     
c     Normalized deposition data will be printed.
c     
            norm_data=.true.
c     
c     Help text - exit afterward
c     
         elseif (opts(in)(1:2).eq.'-h') then
c     
            write(0,'(a)') 'PROC_DEPOSITION:'
            write(0,'(a)') ' Reads in DEPOSITION data'//
     >           ' formatted data output and'//
     >           ' converts it to a compatible EXCEL format'
            write(0,'(a)') ' Options:'
            write(0,'(a)') ' -h   this help text:'
            write(0,'(a)') ' -o   outer target'
            write(0,'(a)') ' -i   inner target'
            write(0,'(a)') ' -w   wall'
            write(0,'(a)') ' -a   <num> axis option'//
     >           ' - 1 to 5 (r,z,d,index,dist along wall)'
            write(0,'(a)') ' -as  <num>  Axis Shift'
            write(0,'(a)') ' -g   include axis for each dataset'
            write(0,'(a)') ' -ns  <num> Separatrix index'
            write(0,'(a)') ' -s   <num> Scale factor'
            write(0,'(a)') ' -n   Normalize deposition to 1'
            write(0,'(a)') ' -d   Data option - 0=/m2-tor 1=raw'
            
c     
            return
c     
         endif
c     
      end do
c     
      write(0,*) 'Finished options'
c     
c     Read in specification file - get number of datasets
c     
      call get_plotlist(filenames,com1,com2,maxsets,nsets)
c     
c     Find number of data elements and allocate storage 
c     
      if (.not.data_allocated) then 
c     
         ndata=get_data_elements(filenames(1),output_region,datunit)
c     
         if (ndata.eq.0) then 
            write(0,*) 'Exiting ...'
            return
         endif
c     
c     Allocate the arrays to store the output data
c     
      write(0,*) 'Before Allocation: Data Points =',ndata,
     >           ' Data Sets =',nsets
c     
         if (include_geometry.eq.0) then 
            allocate(outaxis(ndata,1),stat=flag)
            allocate(outdata(ndata,nsets),stat=flag)
            naxis = 1
         elseif (include_geometry.eq.1) then
            allocate(outaxis(ndata,nsets),stat=flag)
            allocate(outdata(ndata,nsets),stat=flag)
            naxis = nsets
         endif
c     
c     write(0,*) 'After Allocation:'
c     
         if (flag.ne.0) then 
            write(0,*) 'ERROR Allocating storage:', flag
            write(0,*) 'Exiting ...'
            return
         endif
c     
         data_allocated=.true.
c     
      endif
c     
c     Read in and load the data files - read the first one then continue with the rest
c     
      axis_loaded=.false.
c     
      do in=1,nsets
c     
         len=lenstr(filenames(in))
         write(0,*) 'Loading data set : ',in,':',filenames(in)(1:len)
c     
         call load_data(in,outaxis,outdata,ndata,nsets,naxis,
     >        filenames(in),axis_loaded,include_geometry,
     >        output_region,axis_opt,axis_shift,sepnum,datunit,dat_opt,
     >        norm_data,ierr)
c     
         if (ierr.ne.0) then
            write(0,*) 'ERROR found loading data. IERR=',ierr
            write(0,*) 'Exiting ...'
            return
         endif
c     
      end do
c     
c     Write formatted output file
c     
c     Each column with a heading and data 
c     
      if (dat_opt.eq.0) then 
c     

         formt='(a10,40(3x,a18))' 
c     
         select case (axis_opt)
         case(1)
            axis ='R'
            formd='(f10.5,40(3x,g18.8))'
         case(2)
            axis ='Z'
            formd='(f10.5,40(3x,g18.8))'
         case(3)
            axis ='DIST (M)'
            formd='(f10.5,40(3x,g18.8))'
         case(4)
            axis ='INDEX'
            formd='(i10,40(3x,g18.8))'
         case(5)
            axis ='DIST(M)'
            formd='(f10.5,40(3x,g18.8))'
         end select   
c     
         if (include_geometry.eq.0) then 

            write(6,formt) 'AXIS',(com1(in),in=1,nsets)
            write(6,formt)  axis,(com2(in),in=1,nsets)

         elseif (include_geometry.eq.1) then 

            write(6,formt) (('AXIS',com1(in)),in=1,nsets)
            write(6,formt) ((axis,com2(in)),in=1,nsets)

         endif 
c     
c        The following assumes axis data is ordered - either ascending or descending
c        Depending on the axis type selected the code will find the lowest axis value and
c        go from there - either ascending or descending as required. 
c     
c        Find values for wall distance axis
c
         if (axis_opt.eq.5.and.output_region.eq.3) then 
c
c           Find axis minimum
c
            minaxis = HI

            do in = 1,ndata
               if (outaxis(in,1).lt.minaxis) then 
                  minaxis_in = in
                  minaxis = outaxis(in,1)
               endif
            end do
c
            startin = minaxis_in
c
c           Wall distances are counter-clockwise 
c     
            stepin = -1
c
         else
            startin= 1
            stepin = 1
         endif
c        
         in = startin
         icount = 1

         do while (icount.le.ndata)
c           
            if (include_geometry.eq.0) then 
c     
               if (outaxis(in,1).ge.outstart.and.
     >             outaxis(in,1).le.outend) then 
c
                  write(6,formd) outaxis(in,1),
     >                 (outdata(in,is)*scalef,is=1,nsets) 
c
               endif
c     
            elseif (include_geometry.eq.1) then
c     
c              Only filter on first axis set of data if 
c              multiple data axes are provided
c
               if (outaxis(in,1).ge.outstart.and.
     >             outaxis(in,1).le.outend) then 
c
                  write(6,formd) ((outaxis(in,is),
     >              outdata(in,is)*scalef),is=1,nsets) 
c
               endif
c
            endif
c
            icount = icount + 1
            in = in + stepin
            if (in.gt.ndata) in = 1
            if (in.lt.1) in = ndata
c
         end do

c
c         do in = 1,ndata 
c     
c            if (include_geometry.eq.0) then 
c     
c               if (outaxis(in,1).ge.outstart.and.
c     >             outaxis(in,1).le.outend) then 
c
c                  write(6,formd) outaxis(in,1),
c     >                 (outdata(in,is)*scalef,is=1,nsets) 
c
c               endif
c     
c            elseif (include_geometry.eq.1) then
c     
c              Only filter on first axis set of data if 
c              multiple data axes are provided
c
c               if (outaxis(in,1).ge.outstart.and.
c     >             outaxis(in,1).le.outend) then 
c
c                  write(6,formd) ((outaxis(in,is),
c     >              outdata(in,is)*scalef),is=1,nsets) 
c
c               endif
c     
c            endif
c
c         end do
c     
c     Print out region deposition summary
c     
      elseif(dat_opt.eq.1) then  
c     
c     Allocate storage for summary data
c     
c     Wall section regions are hard coded for the 
c     13C grid at the moment.
c
c     Check to see if grid matches the 13C grid data - exit if not
c
         if (ndata.ne.154) then
            write(0,'(a,2i6)') 'ERROR GENERATING REGION SUMMARY:'//
     >              ' WALL DOES NOT MATCH HARD CODED VALUES:', ndata,154
            return
         endif
c     
         nregion = 6
c     
         allocate(dep_summary(nregion+1,nsets),stat=flag)

         call rzero(dep_summary,nregion+1*nsets)
c     
         do is = 1,nsets
c     
            do id = 1,ndata
c     
c     Inner wall 
c     
               if (id.ge.1.and.id.le.6) then  
c     
                  dep_summary(1,is) = dep_summary(1,is) +
     >                 outdata(id,is) 

c     
c     Top Wall 
c     
               elseif (id.ge.7.and.id.le.58) then  
                  
                  dep_summary(2,is) = dep_summary(2,is) +
     >                 outdata(id,is) 


c     
c     Outer Wall
c     
               elseif (id.ge.59.and.id.le.73) then  

                  dep_summary(3,is) = dep_summary(3,is) +
     >                 outdata(id,is) 

c     
c     Outer Target
c     
               elseif (id.ge.74.and.id.le.113) then  

                  dep_summary(4,is) = dep_summary(4,is) +
     >                 outdata(id,is) 

c     
c     PFZ Wall
c     
               elseif (id.ge.114.and.id.le.114) then  

                  dep_summary(5,is) = dep_summary(5,is) +
     >                 outdata(id,is) 


c     
c     Inner Target
c     
               elseif (id.ge.115.and.id.le.154) then  

                  dep_summary(6,is) = dep_summary(6,is) +
     >                 outdata(id,is) 


               endif 

               dep_summary(nregion+1,is) = dep_summary(nregion+1,is) +
     >              outdata(id,is) 

            enddo
            
         enddo 
c     
c     Normalize
c     
         do is = 1,nsets
c     
            do id =1,nregion
c     
               if (dep_summary(nregion+1,is).gt.0.0) then   
                  dep_summary(id,is) = dep_summary(id,is)
     >                 /dep_summary(nregion+1,is)
               else 
                  dep_summary(id,is) = 0.0
               endif

            enddo

         enddo
c     
c     Print out the deposition summary data
c     

         formt='(24x,40(2x,a10))' 
         formd='(2(a11,1x),40(2x,f8.3,2x))'
c     
         write(6,formt) 'INNER WALL',
     >        ' TOP WALL ',
     >        'OUTER WALL',
     >        'OUTER TARG',
     >        ' PFZ WALL ',
     >        'INNER TARG'

         do is = 1,nsets
            
            write(6,formd) com1(is),com2(is),
     >           (dep_summary(in,is),in=1,6)

         end do

c     
         deallocate(dep_summary)
c     
      endif


      if (data_allocated) then  
         deallocate(outaxis)
         deallocate(outdata)
      endif

c     
c     
c     
      return
 200  write(0,*) 'ERROR in Input Option Processing'
      write(0,*) 'Exiting ...'
      return

      end
c     
c     
c     
      integer function get_data_elements(filename,output_region,datunit)
      implicit none
      character*(*) filename
      integer output_region,datunit
c     
c     Open file name - read in the number of elements specification line
c     return the selected value. 
c     
      integer len, lenstr,ierr
      logical data_found
      external lenstr
      character*256 line
      integer nr1,nr2,nr3
c     
      data_found=.false.
c     
      len=lenstr(filename)
c     
      open(datunit,file=filename(1:len),form='formatted',status='old',
     >     iostat=ierr)
c     
      if (ierr.eq.0) then 
c     
c     Read the number of elements line
c     
         do while (.not.data_found.and.ierr.eq.0)
c     
            read(datunit,'(a)',iostat=ierr) line
c     
            if (line(1:39).eq.
     >           'DEPOSITION DATA: NUMBER OF DATA POINTS:') then
               read (line(40:),*) nr1,nr2,nr3
               data_found = .true.          
            endif

         end do
c     
         close(datunit)
c     
         if (data_found) then 
            select case (output_region)
            case(1)
               get_data_elements=nr1
            case(2)
               get_data_elements=nr2
            case(3)
               get_data_elements=nr3
            case default
               get_data_elements=0
            end select
         else
            len=lenstr(filename)
            write(0,*) 'ERROR: Number of data elements not'//
     >           ' found in ',filename(1:len)
            get_data_elements=0
         endif
c     
      else
         len=lenstr(filename)
         write(0,*) 'ERROR Opening :',filename(1:len)
         get_data_elements=0
      endif
c     
      return
      end
c     
c     
c     
      subroutine get_plotlist(filenames,com1,com2,maxsets,nsets)
      implicit none
      integer maxsets,nsets
      character*100 filenames(maxsets),com1(maxsets),com2(maxsets)
c
      character*512 line
c     
c     Read the number of data sets and specifications from UNIT 5
c     
      integer ierr,in
      integer len1,len2,len3,lenstr
      external lenstr
c     
 10   read(5,'(a512)') line
      if (line(1:1).eq.'$') goto 10
c
      read(line,*,iostat=ierr) nsets
c     
      write(0,*) 'Reading Plot Input:',nsets
c     
      if (ierr.eq.0) then  
c     
         do in= 1,nsets
            read(5,'(a512)') line
            if (line(1:1).eq.'$') cycle
c
            read(line,*,iostat=ierr) filenames(in),com1(in),com2(in)
c     
            if (ierr.ne.0) then 
               write(0,*) 'ERROR reading specifications: SET =',in
c     
               len1 = lenstr(filenames(in))
               len2 = lenstr(com1(in))
               len3 = lenstr(com2(in))
c     
               write(0,*) 'Input Data:',filenames(in)(1:len1),' | ',
     >              com1(in)(1:len2),' | ',            
     >              com2(in)(1:len3)
               write(0,*) 'NSETS reset to zero from ',nsets 
               nsets = 0
               return
            endif
         end do 
c     
      else
c     
         write(0,*) 'ERROR reading specifications'
         nsets=0
         return
      endif
c     
      return
      end

      subroutine load_data(iset,outaxis,outdata,ndata,nsets,naxis,
     >     filename,axis_loaded,include_geometry,
     >     output_region,axis_opt,axis_shift,sepnum,
     >     datunit,dat_opt,norm_data,ierr)
      implicit none
      integer iset,include_geometry,ierr,ndata,nsets,naxis
      integer output_region,datunit,axis_opt,sepnum,dat_opt
      logical axis_loaded,norm_data
      real axis_shift
      real outaxis(ndata,naxis)
      real outdata(ndata,nsets)
      character*(*) filename
c     
c     Local variables
c     
      integer in,len,lenstr
      external lenstr
      logical data_section_found
      character*512 line
c     
      integer id,wid
      real r,z,dds,sepdist,walldist,wallsi,wallsn,den_m2_tor
c     
      real mult,max_data_value
      real maxwalldist,minwalldist
c
      real HI
      parameter(HI=1.0e37)
c
c     
c     Initial values
c     
      max_data_value= -HI
      maxwalldist= -HI
      minwalldist= HI 
c     
c     Read in data 
c     
c     
      data_section_found=.false.

      len=lenstr(filename)

      open(datunit,file=filename(1:len),form='formatted',
     >     status='old',iostat=ierr)
c     
c     If no errors and open was successful 
c     
      if (ierr.eq.0) then
c     
c     Need to scan file to the correct header and then 
c     read in the data - it is assumed that it will always 
c     have the same number of elements. On the first time 
c     through the axis needs to be assigned for include_geometry
c     option 0 - when option 1 is supported a different axis will
c     be allowed for each plot element.
c     

         do while (.not.data_section_found.and.ierr.eq.0)
c     
            read(datunit,'(a)',iostat=ierr) line
c     
            if (ierr.eq.0) then 
c     
               if (
     >              (line(1:32) .eq. 
     >              'DEPOSITION DATA for OUTER target'
     >              .and.output_region.eq.1)
     >              .or.            
     >              (line(1:32) .eq. 
     >              'DEPOSITION DATA for INNER target'
     >              .and.output_region.eq.2)
     >              .or.
     >              (line(1:31) .eq. 
     >              'DEPOSITION DATA for Entire Wall'
     >              .and.output_region.eq.3)) then
c     
                  data_section_found=.true.
c     
               endif
c     
c     
c     Error reading from file
c     
            else
c     
               len=lenstr(filename)
               write(0,*) 'Error Reading Data: ',filename(1:len)
               close(datunit)
               return
c     
            endif
c     
         end do

c     
c     Read data section
c     
         if (ierr.eq.0.and.data_section_found) then 
c     
c     Read 2 junk lines in all formats
c     
            read(datunit,'(a)',iostat=ierr) line
            read(datunit,'(a)',iostat=ierr) line
c     
            do in = 1,ndata
c     
               if (output_region.eq.1.or.output_region.eq.2) then                                
                  read(datunit,200)
     >                 id, wid,
     >                 r,z,dds,sepdist,walldist,
     >                 wallsi,wallsn,den_m2_tor
               elseif (output_region.eq.3) then 

                  read(datunit,300)
     >                 id, wid,
     >                 r,z,sepdist,walldist,
     >                 wallsi,wallsn,den_m2_tor
               endif
c
               write(8,200) 
     >                 id, wid,
     >                 r,z,dds,sepdist,walldist,
     >                 wallsi,wallsn,den_m2_tor
c

c     
c     Depending on axis options - store data for now
c     
c     Note: AXIS LOADED only applies to include geometry 0 
c     
               if (.not.axis_loaded.and.include_geometry.eq.0) then 
c     
c     Store axis data - further
c     processing is required for option 3
c     once all data has been loaded 
c     
                  select case (axis_opt)
                  case (1)
                     outaxis(in,1) = r  
                  case (2)
                     outaxis(in,1) = z
                  case (3)
                     outaxis(in,1) = sepdist
                  case (4)
                     outaxis(in,1) = in
                  case (5)
                     outaxis(in,1) = walldist

                     if (output_region.eq.3) then 
                        maxwalldist = max(maxwalldist,
     >                                 walldist + sepdist/2.0)
                        minwalldist = min(minwalldist,
     >                                 walldist-sepdist/2.0)
                     endif

                  end select
c     
               elseif (include_geometry.eq.1) then 
c     
c     Store axis data - further
c     processing is required for option 3
c     once all data has been loaded 
c     
                  select case (axis_opt)
                  case (1)
                     outaxis(in,iset) = r  
                  case (2)
                     outaxis(in,iset) = z
                  case (3)
                     outaxis(in,iset) = sepdist
                  case (4)
                     outaxis(in,iset) = in
                  case (5)
                     outaxis(in,iset) = walldist

                     if (output_region.eq.3) then 
                        maxwalldist = max(maxwalldist,
     >                                 walldist + sepdist/2.0)
                        minwalldist = min(minwalldist,
     >                                 walldist-sepdist/2.0)
                     endif

                  end select
c     
               endif
c     
c     Save data
c     
               if (dat_opt.eq.0) then 

                  outdata(in,iset) = den_m2_tor

               elseif (dat_opt.eq.1) then 

                  outdata(in,iset) = wallsi+wallsn

               endif
c     
               write(8,'(i6,f12.6,g18.6)') in,outaxis(in,1),
     >                                     outdata(in,iset)

               max_data_value = max(max_data_value,outdata(in,iset))
c     
            end do   
c     
c     Finish processing on axes if required then mark as complete
c     
            if (.not.axis_loaded.and.include_geometry.eq.0) then 
c     
c     Deal with case 3 axis type
c     
               if (axis_opt.eq.3) then 
c     
c     
c     Wall - calculate along wall distances
c     
                  if (output_region.eq.3) then 
c     
                     sepdist = 0.0
c     
                     do in = 1,ndata
c     
c     Take half of segment length 
c     
                        dds = outaxis(in,1)/2.0
c     
c     Add this length before and after the center location
c     
                        sepdist = sepdist+dds
c     
                        outaxis(in,1) = sepdist 
c     
                        sepdist = sepdist+dds
c     
                     end do

                  endif

               endif
c     
c     Apply axis_shift if one is specified
c     
               if (axis_shift.ne.0.0) then 
                  do in = 1,ndata
c     
c                   Deal with wall distance axis for full wall data
c 
                     if (axis_opt.eq.5.and.output_region.eq.3) then 
c
                        if ((outaxis(in,1)+axis_shift).lt.
     >                       minwalldist) then 

                           outaxis(in,1) = outaxis(in,1) + axis_shift 
     >                                     + (maxwalldist-minwalldist)

                        elseif((outaxis(in,1)+axis_shift).gt.
     >                          maxwalldist) then 

                           outaxis(in,1) = outaxis(in,1) + axis_shift 
     >                                     - (maxwalldist-minwalldist)

                        else
                           
                           outaxis(in,1) = outaxis(in,1) + axis_shift

                        endif
                     else
                        outaxis(in,1) = outaxis(in,1) + axis_shift
                     endif

                  end do
               endif

               axis_loaded = .true.
c     
            elseif(include_geometry.eq.1) then
c     
c     Deal with case 3 axis type
c     
               if (axis_opt.eq.3) then 
c     
c     Wall - calculate along wall distances
c     
                  if (output_region.eq.3) then 
c     
                     sepdist = 0.0
c     
                     do in = 1,ndata
c     
c     Take half of segment length 
c     
                        dds = outaxis(in,iset)/2.0
c     
c     Add this length before and after the center location
c     
                        sepdist = sepdist+dds
c     
                        outaxis(in,iset) = sepdist 
c     
                        sepdist = sepdist+dds
c     
                     end do

                  endif

               endif
c     
c     
c     Apply axis_shift if one is specified
c     
               if (axis_shift.ne.0) then 
c     
c                   Deal with wall distance axis for full wall data
c 
                     if (axis_opt.eq.5.and.output_region.eq.3) then 
c
                        if ((outaxis(in,iset)+axis_shift).lt.
     >                       minwalldist) then 

                           outaxis(in,iset)=outaxis(in,iset)+axis_shift 
     >                                     + (maxwalldist-minwalldist)

                        elseif((outaxis(in,iset)+axis_shift).gt.
     >                          maxwalldist) then 

                           outaxis(in,iset)=outaxis(in,iset)+axis_shift 
     >                                     - (maxwalldist-minwalldist)

                        else
                           
                           outaxis(in,iset)=outaxis(in,iset)+axis_shift

                        endif
                     else
                        outaxis(in,iset)=outaxis(in,iset)+axis_shift
                     endif

               endif

            endif 
c     
c     Normalize the data set if requested 
c     
            if (norm_data.and.(max_data_value.ne.0.0)) then 
c     
               do in = 1,ndata
c     
                  outdata(in,iset) = outdata(in,iset)/max_data_value           
c     
               end do
c     
            endif
c     
         else
c     
            len=lenstr(filename)
            write(0,*) 'Error Reading Data: ',filename(1:len)
            close(datunit)
            return
c     
         endif

c     
c     Error opening file 
c     
      else
         len=lenstr(filename)
         write(0,*) 'Error Opening: ',filename(1:len)
         return
      endif
c     
c     Close data file
c     
      close(datunit)


      write(8,*) trim(filename)
      do in = 1,ndata
         write(8,'(f12.6,g18.6)') outaxis(in,1),outdata(in,iset)
      end do



                                ! Labels
                                !       'ID','IW','R','Z','DDS','SEPDIST','WALLDIST',
                                !       'ION-DEP','NEUT-DEP','DEN/M2-TOR'
 100  format(5x,a,5x,a,10x,a,10x,a,8x,a,4x,a,3x,a,8x,a,11x,a,10x,a)
                                ! Target Data
 200  format(2(1x,i6),5(1x,f10.6),3(1x,g18.8))
                                ! Wall data
 300  format(2(1x,i6),3(1x,f10.6),11x,(1x,f10.6),3(1x,g18.8))

c     
      return
      end

