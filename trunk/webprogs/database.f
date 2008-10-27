c      module logrecorddata
c
c         integer maxrefs
c         parameter(maxrefs=10)
c         type logrecord 
c            character*100 casename
c            character*100 seriesname
c            character*20  submitted  
c            character*20  date
c            character*20  time
c            character*512 desc
c            integer       nrefs
c            character*100 refs(maxrefs)
c            integer nextrecord 
c            integer lastrecord  
c         end type logrecord
c
c      end module logrecorddata   
c
c
c
      module logdatabase 
c
         use logrecorddata 
c
         integer maxcases
         parameter(maxcases=1000)
c
c         common /database/ logdata,nrecords,firstrecord,endrecord
c
         type (logrecord) logdata(maxcases)
         integer nrecords,firstrecord,endrecord
c
      end module logdatabase 
c
c     Case log database maintenance routines
c
      subroutine writedb(dbunit,tempdb)
      use logdatabase
      implicit  none
      integer dbunit
      character*(*) tempdb
c
c      include 'logdatabase'
c
c     Loops through the data base writing all of the records to 
c     the temporary file.
c      
      integer lenf,lenstr,ios,nextindex,lastindex,curindex,len,in
      external lenstr,nextindex
c 
      lenf = lenstr(tempdb)
c
      open(dbunit,FILE=tempdb(1:lenf),IOSTAT=ios)
c
      if (ios.ne.0) then 
c
         write(0,*) 'ERROR: Could not write database. Program exiting.'
         stop 
c 
      endif
c
c     Write information to database
c 
      write(dbunit,'(''CASELOG: '',i7)') nrecords
c
c     Write out non-empty database  
c
      if (nrecords.gt.0) then 
c
c        Set lastindex to beginning of file - curindex to first element
c 
         lastindex= 0
         curindex = nextindex(0)
c 
         do while (curindex.gt.0) 
c
c           Write out series name
c
            if (lastindex.eq.0) then  

               len = lenstr(logdata(curindex)%seriesname)  
               write(dbunit,'(''SERIES:'',a)')
     >               logdata(curindex)%seriesname(1:len)

            elseif(logdata(curindex)%seriesname.ne.
     >             logdata(lastindex)%seriesname) then 

               len = lenstr(logdata(curindex)%seriesname)  
               write(dbunit,'(''SERIES:'',a)')
     >               logdata(curindex)%seriesname(1:len)

            endif
c
c           Write out specific case information
c
c           CASE name
c
            len = lenstr(logdata(curindex)%casename)  
            write(dbunit,'(''CASE:'',a)')
     >               logdata(curindex)%casename(1:len)
c
c           Submitted by
c
            len = lenstr(logdata(curindex)%submitted)  
            write(dbunit,'(''SUBMITTED:'',a)')
     >               logdata(curindex)%submitted(1:len)
c
c           Date submitted
c
            len = lenstr(logdata(curindex)%date)  
            write(dbunit,'(''DATE:'',a)')
     >               logdata(curindex)%date(1:len)
c
c           Time submitted
c
            len = lenstr(logdata(curindex)%time)  
            write(dbunit,'(''TIME:'',a)')
     >               logdata(curindex)%time(1:len)
c
c           Case Description 
c
            len = lenstr(logdata(curindex)%desc)  
            write(dbunit,'(''DESC:'',a)')
     >               logdata(curindex)%desc(1:len)
c
c           Number of hypertext references available if any
c
            write(dbunit,'(''NREFS:'',i6)')
     >               logdata(curindex)%nrefs
c
c           String containing names of files to link to
c
            do in = 1, logdata(curindex)%nrefs 
               len = lenstr(logdata(curindex)%refs(in))  
               write(dbunit,'(''REFS:'',a)')
     >               logdata(curindex)%refs(in)(1:len)
            end do 
c
c           Go onto next record
c
            curindex = nextindex(curindex)
c
         end do
c
      endif
c
      close(dbunit) 
c
      return
      end
c
c
c
      subroutine addcase(newentry)
      use logdatabase
      implicit none 
c
c      include 'logdatabase'
c
      type (logrecord) newentry 
c
c     Add a case record to the in-memory database 
c
      integer newindex,getindex,in
      external getindex
      logical new
c
c     Does case exist already - if so just update refs
c      
c      write(0,*) 'Adding case to database:' 
c
      newindex = getindex(newentry%seriesname,newentry%casename,new)
c
c      write(0,*) 'At INDEX:',newindex 
c
c
c     Assign all values to the database 
c
      if (new) then
c
         logdata(newindex)%casename = newentry%casename
         logdata(newindex)%seriesname = newentry%seriesname
         logdata(newindex)%submitted = newentry%submitted
         logdata(newindex)%date = newentry%date
         logdata(newindex)%time = newentry%time
         logdata(newindex)%desc = newentry%desc
c
         logdata(newindex)%nrefs = newentry%nrefs
c
c
         do in = 1,newentry%nrefs 

            logdata(newindex)%refs(in) = newentry%refs(in)
c
         end do
c 
c     Just load the references in case they have changed or new plots have 
c     been added. 
c
      else 

         logdata(newindex)%nrefs = newentry%nrefs
c
         do in = 1,newentry%nrefs 

            logdata(newindex)%refs(in) = newentry%refs(in)
c
         end do
c 
      endif
c
      return
c
      end   
c
      integer function nextindex(curindex)
      use logdatabase
      implicit none
      integer curindex 
c
c      include 'logdatabase'
c 
      if (curindex.eq.0) then 
         nextindex = firstrecord
      elseif(curindex.eq.-1) then
         nextindex = endrecord
      elseif (curindex.ge.1.and.curindex.le.nrecords) then 
         nextindex = logdata(curindex)%nextrecord
      else
         nextindex = 0
      endif  
c
      return
      end 
c
c
c
      integer function getnrefs(curindex)
      use logdatabase
      implicit none
      integer curindex 
c
c     Returns the nrefs field of the current database index
c
      if (curindex.ge.1.and.curindex.le.nrecords) then 

         getnrefs = logdata(curindex)%nrefs

      else

         getnrefs = 0

      endif 
c
      return
      end
c
c
c
      integer function getindex(seriesname,basename,new)
      use logdatabase
      implicit none
      character*(*) seriesname,basename
      logical new
c
c      include 'logdatabase'
c
c     Loop through the database looking for a match to the case name
c     or the appropriate place to insert the new item. 
c            
      integer curindex,nextindex,len1,len2,lenstr
      intrinsic llt 
      external nextindex,lenstr
      logical found,llt
c
c     Start at beginning of database  
c      
      new = .false. 
      found = .false. 
      getindex = 0  
      curindex = nextindex(0)
c
c     Empty database 
c
      if (curindex.eq.0) then 
c
         firstrecord = 1
         nrecords = 1
         endrecord = 1
c
         getindex = 1
c
         logdata(getindex)%lastrecord = 0
         logdata(getindex)%nextrecord = -1
         new = .true.
         found = .true.
c
c     Some database entries - search for match or appropriate insert 
c     position
c
      else
c
         do while (.not.found) 
c
c           Insert record at end
c
            if (curindex.eq.-1) then 
c
               new = .true.
               found = .true.
c
c              Add new record
c 
               nrecords = nrecords +1 
               getindex = nrecords
c
c              Update links
c                
               logdata(getindex)%lastrecord = endrecord
               logdata(getindex)%nextrecord = -1
c
               logdata(endrecord)%nextrecord = getindex
c
               endrecord = getindex
c
            else

               len1 = lenstr(basename) 
               len2 = lenstr(logdata(curindex)%casename)
c
c
c              Do lexical name search to keep order
c  
               if (basename(1:len1).eq.
     >            logdata(curindex)%casename(1:len2)) then
c
                  new = .false.
                  found = .true.
                  getindex = curindex
c
c              Insert new record before next record
c
               elseif (llt(basename(1:len1),
     >            logdata(curindex)%casename(1:len2))) then 
c
c                 Get new record index - add at end
c
                  new = .true.
                  found = .true.
c 
                  nrecords = nrecords +1 
                  getindex = nrecords
c
c                 Check if first record
c                
                  if (curindex.eq.firstrecord) then
c
                     logdata(getindex)%lastrecord = 0
                     logdata(getindex)%nextrecord = firstrecord
c
                     logdata(firstrecord)%lastrecord = getindex
c
                     firstrecord = getindex
c
                  else
c
c                    Update current record links
c
                     logdata(getindex)%lastrecord = 
     >                       logdata(curindex)%lastrecord
                     logdata(getindex)%nextrecord = curindex
c
c                    Update adjacent record links
c
                     logdata(logdata(getindex)%lastrecord)%nextrecord=
     >                     getindex
                     logdata(logdata(getindex)%nextrecord)%lastrecord= 
     >                     getindex
c
                  endif 
c
               else
c
c              Get next index  
c
                  curindex = nextindex(curindex)

               endif                  

            endif

         end do

      endif 
c
      return
      end
c
c
c
      subroutine initdb
      use logdatabase
      implicit none
c
c     Initialize the record counting related quantities in the database 
c
c      include 'logdatabase'
c
      nrecords = 0
      firstrecord = 0
      endrecord = 0 
c
      return
      end 
c      
c
c
      subroutine loadlogdata(logunit)
      use logdatabase
      implicit none
      integer logunit
c
c      include 'logdatabase' 
c
c     Local variables
c
      integer ios,len,lenstr,nlogs,in
      external lenstr 
      character*512 line
      character*100 series  
c
c     Init
c
      ios = 0
c
c     Empty database with zero records will be trapped as an error 
c     Execution continues
c
      do while(nrecords.eq.0.and.ios.eq.0)

         read(logunit,'(a512)',end=100,err=100,IOSTAT=ios) line
c
         if (line(1:8).eq.'CASELOG:') then 
            len = lenstr(line)
            read(line(9:len),*) nrecords
         endif
c
      end do
c
c     Initialize record count
c
      nlogs = 0
c
c     Read in the database records
c
      do while(ios.eq.0) 
c
         read(logunit,'(a512)',IOSTAT=ios) line
c
         if (ios.eq.0.and.line(1:1).ne.'$'.and.
     >       line(1:1).ne.'!'.and.
     >       line(1:1).ne.'#') then 
c
c           Line length 
c
            len = lenstr(line) 
c
c           Extract the pieces of the individual records
c
c           Series name
c
            if (line(1:7).eq.'SERIES:') then 
c
               series = line(8:len)
c
c           Case name  
c
            elseif(line(1:5).eq.'CASE:') then
               nlogs = nlogs +1
               logdata(nlogs)%casename=line(6:len)
               logdata(nlogs)%seriesname=series
               logdata(nlogs)%lastrecord = nlogs -1 
               logdata(nlogs)%nextrecord = nlogs +1
c
c              Initialize rest of record  
c
               logdata(nlogs)%submitted= ' '
               logdata(nlogs)%date= ' '
               logdata(nlogs)%time= ' '
               logdata(nlogs)%desc= ' '
               logdata(nlogs)%refs= ' '
               logdata(nlogs)%nrefs= 0
               do in = 1,maxrefs 
                  logdata(nlogs)%refs(in)= ' '
               end do
c
c           Submitted by who  
c
            elseif(line(1:10).eq.'SUBMITTED:') then 
               logdata(nlogs)%submitted=line(11:len)
c
c           Time submitted
c
            elseif(line(1:5).eq.'TIME:') then 
               logdata(nlogs)%time=line(6:len)
c
c           Date submitted
c
            elseif(line(1:5).eq.'DATE:') then
               logdata(nlogs)%date=line(6:len)
c
c           Case description - extracted from .dat/.html file
c
            elseif(line(1:5).eq.'DESC:') then 
               logdata(nlogs)%desc=line(6:len)
c
c           Hypertext references for the html case log
c
            elseif(line(1:5).eq.'NREFS:') then 
               logdata(nlogs)%nrefs= 0
c
c            elseif(line(1:5).eq.'REFS:') then 
c               logdata(nlogs)%refs=' '
c
            endif
c
         endif
c
      end do     
c
c     Fix up link of last record if non-zero - set tail to -1
c
      if (nlogs.gt.0) then  

         logdata(nlogs)%nextrecord = -1

      endif 
c
      if (nlogs.ne.nrecords) then 
         write(0,*) 'LOADLOGDATA: Mismatch in'//
     >              ' number of database records',nrecords,nlogs
      endif
c
c     record number of items read.
c
      nrecords = nlogs 
c
c     Set first and last database pointers 
c
      firstrecord = 1
      endrecord = nrecords
c
      return 
c
c     Error condition 
c

 100  nrecords = 0
      write(0,*) 'LOG FILE not complete',ios
      return    
      end 
c
c
c
      subroutine getrecord(index,item)
      use logdatabase
      implicit none
      integer index  
c
c      include 'logdatabase'
c
      type (logrecord) item
c
c     This routine loads the database record at the given index
c     into the item that has been passed in.
c 
      item = logdata(index)
c
      return
      end 
