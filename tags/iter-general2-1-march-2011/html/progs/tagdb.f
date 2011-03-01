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
      write(dbunit,'(''TAGLOG: '',i7)') nrecords
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
c
c           Write out specific TAG information
c
c           TAG name
c
            len = lenstr(logdata(curindex)%tag)  
            write(dbunit,'(''CASE:'',a)')
     >               logdata(curindex)%tag(1:len)
c
c           HTML Reference
c
            len = lenstr(logdata(curindex)%htmlref)  
            write(dbunit,'(''HTML:'',a)')
     >               logdata(curindex)%htmlref(1:len)
c
c           Bookmark
c
            len = lenstr(logdata(curindex)%bookmark)  
            write(dbunit,'(''BOOK:'',a)')
     >               logdata(curindex)%bookmark(1:len)
c
c           Input Description
c
            len = lenstr(logdata(curindex)%desc)  
            write(dbunit,'(''DESC:'',a)')
     >               logdata(curindex)%desc(1:len)
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
      subroutine addrecord(newentry)
      use logdatabase
      implicit none 
c
c      include 'logdatabase'
c
      type (logrecord) newentry 
c
c     Add a case record to the in-memory database 
c
      integer newindex,getindex,in,len1,len2,len3,lenstr
      external getindex,lenstr
      logical new
c
c     Does case exist already - if so just update refs
c      
c      write(0,*) 'Adding case to database:',newentry%tag(1:3)
c
      newindex = getindex(newentry%tag,new)
c
c      write(0,*) 'At INDEX:',newindex 
c
c
c     Assign all values to the database 
c
      if (new) then
c
         logdata(newindex)%tag = newentry%tag
         logdata(newindex)%htmlref = newentry%htmlref
         logdata(newindex)%bookmark = newentry%bookmark
         logdata(newindex)%desc = newentry%desc
c 
c     Issue error message about duplicate tags
c
      else 
c
         write(0,*) 'Duplicate TAG found:'
c
         len1 = lenstr(logdata(newindex)%tag)
         len2 = lenstr(logdata(newindex)%htmlref)
         len3 = lenstr(logdata(newindex)%bookmark)
c
         write(0,'(a,a,a,a)') 'First TAG :',
     >              logdata(newindex)%tag(1:len1),
     >              logdata(newindex)%htmlref(1:len2),
     >              logdata(newindex)%bookmark(1:len3)
c
         len1 = lenstr(newentry%tag)
         len2 = lenstr(newentry%htmlref)
         len3 = lenstr(newentry%bookmark)
c
         write(0,'(a,a,a,a)') 'Second TAG:',
     >              newentry%tag(1:len1),
     >              newentry%htmlref(1:len2),
     >              newentry%bookmark(1:len3)
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
      integer function getindex(tag,new)
      use logdatabase
      implicit none
      character*(*) tag
      logical new
c
c      include 'logdatabase'
c
c     Loop through the database looking for a match to the tag
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

               len1 = lenstr(tag) 
               len2 = lenstr(logdata(curindex)%tag)
c
c
c              Do lexical name search to keep order
c  
               if (tag(1:len1).eq.
     >            logdata(curindex)%tag(1:len2)) then
c
                  new = .false.
                  found = .true.
                  getindex = curindex
c
c              Insert new record before next record
c
               elseif (llt(tag(1:len1),
     >            logdata(curindex)%tag(1:len2))) then 
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
         if (line(1:7).eq.'TAGLOG:') then 
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
            if(line(1:4).eq.'TAG:') then
c
               nlogs = nlogs +1
               logdata(nlogs)%tag=line(5:len)
               logdata(nlogs)%lastrecord = nlogs -1 
               logdata(nlogs)%nextrecord = nlogs +1
c
c              Initialize rest of record  
c
               logdata(nlogs)%htmlref= ' '
               logdata(nlogs)%bookmark= ' '
               logdata(nlogs)%desc= ' '
c
c           Submitted by who  
c
            elseif(line(1:5).eq.'HTML:') then 
               logdata(nlogs)%htmlref=line(6:len)
c
c           Time submitted
c
            elseif(line(1:5).eq.'BOOK:') then 
               logdata(nlogs)%bookmark=line(6:len)
c
c           Date submitted
c
            elseif(line(1:5).eq.'DESC:') then
               logdata(nlogs)%desc=line(6:len)
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
