      module logrecorddata 

         integer maxrefs
         parameter(maxrefs=20)
         type logrecord 
            character*100 casename
            character*100 seriesname
            character*20  submitted  
            character*20  date
            character*20  time
            character*512 desc
            integer       nrefs
            character*100 refs(maxrefs)
            integer nextrecord 
            integer lastrecord  
         end type logrecord

      end module logrecorddata   
c
c
c
      program buildindex
      use logrecorddata 
      implicit none
      integer numitems
      parameter(numitems=10)
c
c     Locations ... and file names
c
      character*100 resdir,fname,logname,dbname,defresdir
      character*100 caseback,caselock
      character*10 htmlstr,pdfstr,datstr
      character*100 tempfile,tempindex,tempdir,tempplot,templog,
     >              tempdb
      character line*256
      character*100 a(numitems) 
c
      character*200 indextitle,defindextitle
c
c     Table entry data 
c
c      include 'logrecord' 
c
      type (logrecord) newentry
c
      character casename*100,submitted*20,date*20,
     >          time*20,description*512
      integer nplots 
c
c     Indices and unit numbers
c
      integer rc,ios,nfound
      integer ierr,in
      integer len,lenstr
      external lenstr,iargc,openlock,removelock
      integer nargs,iargc,openlock,removelock
      character*(256) arg  
c
      integer logunit,tmpunit,inunit,outunit,pltunit,dbunit
c
      data tmpunit /8/
      data inunit  /9/
      data outunit /10/
      data pltunit /11/
      data logunit /12/
      data dbunit  /13/
c
c     Debugging
c
      logical debug  
c
c     Initilaization 
c
      debug =  .false.
      defresdir = '/nfs/divimp/web/results' 
      defindextitle= 'DIVIMP Results Index'
c
c     Fixed file names
c
      fname     = 'result_index.html' 
      logname   = 'caselog.html'
      dbname    = 'caselog.db'
      caseback  = 'caselog.db.back'
      caselock  = 'caselog.lock' 
c
c     Temporary file names and location
c
      tempdir   = '/tmp/'
c
      tempfile  = 'buildindex.tmp'
      tempindex = 'buildindex.index'
      tempplot  = 'buildindex.plts' 
      templog   = 'caselog.tmp' 
      tempdb    = 'casedb.tmp' 
c
c     File search extensions  
c
      htmlstr   = '*.htm?'
      datstr    = '*.dat'
c
c     Read in any command line arguments
c
      nargs = iargc()
c
      if (nargs.eq.1) then
c
c        Directory is first argument    
c
         call getarg(1,arg)
c
         len = lenstr(arg)
         resdir = arg(1:len)
c
      else
c
c        Read in any environment variables if present.
c
         call getenv('BUILDINDEXDIR',resdir) 
c
         call getenv('BUILDINDEXTITLE',indextitle) 
c
      endif 
c
c     Fix up resdir to make sure it ends with a / 
c
      len = lenstr(resdir)
c
      if (len.le.1) then 
         resdir = defresdir
         len = lenstr(resdir)
      endif
c
      if (resdir(len:len).ne.'/') then 
         resdir(len+1:len+1) = '/'
      endif 
c
c     Set indextitle to default if not set in environment 
c
      len = lenstr(indextitle)
c
      if (len.le.1) then 
         indextitle = defindextitle
      endif
c
c     Combine parts to create absolute path names  
c
      call makefilename(resdir,fname,fname)
      call makefilename(resdir,dbname,dbname)
      call makefilename(resdir,logname,logname)
      call makefilename(resdir,caselock,caselock)
      call makefilename(resdir,caseback,caseback)
c
c     Create a lock file in resdir 
c
      rc = openlock(caselock,dbname,caseback)
c
c     Exit code with error message if lock failed 
c
      if (rc.ne.0) then
c
         write(0,*) 'Build Update failed - LOCK not available'
         stop 'Failed: LOCK not available'
c
      endif
c
c     Open and load the database file
c
      call openlog(dbunit,dbname)
c
c     Generate temporary file names - and create the files
c
      call makenames(tempdir,tempfile,tempindex,tempplot,templog,
     >               tempdb)
c
c     Generate file list 
c
      call makelist(resdir,htmlstr,datstr,tempfile)
c
c     Open output file from previous command 
c        
      open(tmpunit,FILE=tempfile,ERR=100,IOSTAT=ios)
c
      do while (ios.eq.0) 
c
c        Loop through list extracting required information from HTML files
c
         read(tmpunit,'(a256)',IOSTAT=ios) line
c
         if (debug) then 
            len = lenstr(line) 
            write(0,*) 'DEBUG LINE:',line(1:len) 
         endif
c
         if (ios.eq.0) then 
c
            call parseline(line,a,numitems,nfound)
c
            if (debug) then 
               len = lenstr(a(nfound)) 
               write(0,*) 'DEBUG A:',a(nfound)(1:len) 
            endif
c
c           Extract and write information to temporary index file getdata
c
c            call getdata(inunit,pltunit,resdir,a,nfound,casename,
c     >                   plotnames,maxplots,nplots,tempplot,
c     >                   submitted,date,time,description,ierr)
c
            call getdata(inunit,pltunit,resdir,a,nfound,
     >                   tempplot,newentry,debug,ierr)
c
            if (ierr.eq.0) then  
c
c              Update database 
c
               call addcase (newentry)
c
            endif
c  
         endif

      end do
c
      close(tmpunit) 
c
c     Write out database file 
c
      call writedb(dbunit,tempdb)
c
c     Write out HTML Caselogs 
c     
      call writelog(logunit,templog,0,'DIVIMP POSTED CASE LOG')
      call writelog(outunit,tempindex,1,indextitle)
c
c     Clean up temp files and copy index file to results directory  
c
      call savefiles(fname,dbname,logname,
     >               tempindex,tempfile,tempplot,tempdb,templog) 
c
      rc = removelock(caselock)
c
      stop 'Finished Builing Index' 

 100  continue
      
      write (0,*) 'File OPEN Error:',ios,tempfile
      stop 'Error on OPEN' 
      end
c
c
c
      subroutine writelog(logunit,templog,opt,title)
      use logrecorddata
      implicit none
      integer logunit,opt
      character*(*) templog,title
c
c     Write the posted case database into a tabulated HTML file.
c 
c     OPT = 0 - write all database entries to html table
c     OPT = 1 - write only those entries with references
c
c     Local declarations
c
c      include 'logrecord'
      type (logrecord) item 
c
      integer curindex,nextindex,getnrefs,len,lenstr
      external nextindex,getnrefs,lenstr
c
c     Open the file and add HTML header
c    
      len = lenstr(title)
c
      call openhtml(logunit,templog,title(1:len))
c 
c     Start table 
c
      call tablestart(logunit,title(1:len))
c
c     Loop through database getting and printing each record   
c
      curindex = nextindex(0)
c
c     Only fill out table for non-empty database
c
      if (curindex.gt.0) then 
c
         do while(curindex.ne.-1) 
c
            call getrecord(curindex,item)
c
c            write (6,*) 'LOG:',curindex, getnrefs(curindex) 
c
            if (opt.eq.0.or.
     >         (opt.eq.1.and.getnrefs(curindex).gt.0)) then 

               call tableentry(logunit,item)
c
            endif
c
            curindex=nextindex(curindex)

         end do 
c
      endif
c
c     End table 
c
      call tableend(logunit)
c
c     Close the file and finish off HTML 
c
      call closehtml(logunit)
c
      return
      end
c
c
c
      subroutine makefilename(basename,fname,resname)
      implicit none
      character*(*) basename,fname,resname
c
c     Combine basedir and fname to get compound file name - place result in 
c     third argument
c
      integer lenb,lenf,lenstr
      external lenstr 
c
      lenb = lenstr(basename)
      lenf = lenstr(fname)
c
      resname = basename(1:lenb)//fname(1:lenf)
c
      return
      end 
c
c
c
      integer function openlock(caselock,casedb,caseback)
      implicit none
      character*(*) caselock,casedb,caseback
c
c     For now - just check if the file exists - if it does then 
c     return an error - if it doesn't then create it. This may be enhanced
c     later to add the directory to the file and only fail if the present
c     directory is listed.
c
c
c     Locals
c
      integer len1,len2,len3,lenstr,rc
      external lenstr 
      logical found
c
      len1 = lenstr(caselock)
c
      inquire(FILE=caselock(1:len1),EXIST=found)
c
      if (found) then
c
         openlock =1 
c
      else
c
         openlock = 0
c
         call cissue(' touch '//caselock(1:len1),rc) 
         call cissue('chmod g+rw '//caselock(1:len1),rc) 
c
         openlock = rc
c        
c        If openlock is successful - copy database to backup before 
c        code can make changes  
c
         if (rc.eq.0) then 
c
            len1 = lenstr(casedb)
            len2 = lenstr(caseback)
c
            call cissue('cp '//casedb(1:len1)//
     >                  ' '//caseback(1:len2),rc) 
            call cissue('chmod g+rw '//caseback(1:len2),rc) 
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
      integer function removelock(caselock)
      implicit none 
      character*(*) caselock 
c
c     Remove the lock file as the last action
c
      integer len1,lenstr,rc
      external lenstr
c
      len1 = lenstr(caselock)
c
      call cissue('rm '//caselock(1:len1),rc) 
c
      removelock = rc
c
      return
      end
c
c
c
      subroutine openlog(dbunit,caselog)
      implicit none
c
c     Load log to in-memory database  
c
      integer nrecords,dbunit
      character*(*) caselog
c
c     Local variables
c
      integer len,ios,lenstr
      external lenstr 
c
c     Initialize record count 
c
      call initdb
c
c     Open file
c
      len = lenstr(caselog)
      open(dbunit,FILE=caselog,ERR=100,IOSTAT=ios)
c
      if (ios.eq.0) then
         call loadlogdata(dbunit)
      endif
c
      close(dbunit)
c  
      return

 100  write(0,*) 'LOG FILE not found',ios
      return    
      end 
c
c
c
      subroutine savefiles(fname,dbname,logname,
     >               tempindex,tempfile,tempplot,tempdb,templog) 
      implicit none
      character*(*) fname,dbname,logname,
     >              tempindex,tempfile,tempplot,tempdb,templog
c
c
c     Save the output files and perform cleanup by deleting 
c     temporary files
c
c     Local variables
c             
      integer len2,len3,len4,lenstr,rc
      external lenstr
c
c     Record string lengths
c
c
c     Set permissions and move to results diectory
c
c     HTML index file 
c
      len2 = lenstr(fname)
      len3 = lenstr(tempindex)
      call cissue('chmod g+rw '//tempindex(1:len3),rc) 
      call cissue('mv '//tempindex(1:len3)//' '
     >                 //fname(1:len2),rc)
c
c     CASELOG HTML file
c
      len2 = lenstr(logname)
      len3 = lenstr(templog)
      call cissue('chmod g+rw '//templog(1:len3),rc) 
      call cissue('mv '//templog(1:len3)//' '
     >                 //logname(1:len2),rc)
c
c     CASELOG database file 
c
      len2 = lenstr(dbname)
      len3 = lenstr(tempdb)
      call cissue('chmod g+rw '//tempdb(1:len3),rc) 
      call cissue('mv '//tempdb(1:len3)//' '
     >                 //dbname(1:len2),rc)
c
c     Delete other temporary files  
c
      len4 = lenstr(tempfile)
      call cissue(' rm '//tempfile(1:len4),rc)
      len4 = lenstr(tempplot)
      call cissue(' rm '//tempplot(1:len4),rc)
c
      return
      end  
c
c
c
      subroutine makenames(tempdir,name1,name2,name3,name4,name5)
      implicit none
      character*(*) tempdir,name1,name2,name3,name4,name5
c
c     This routine generates several temporary file names 
c     that are reasonably unique. The tempdir is 
c     combined with a base name - if any - supplied in 
c     name1 and name2. If these files do not exist already
c     then they are created as a lock mechanism - if they do exist
c     already then variations with random numbers are created 
c     until a unique one is found. 
c
      character*(512) command,randstr,tmpname1,tmpname2,tmpname3,
     >                tmpname4,tmpname5
      integer rc,len,len1,len2,len3,len4,len5,lenr
      integer tlen1,tlen2,tlen3,tlen4,tlen5,lenstr
      real x
      external lenstr
      logical found1,found2,found3,found4,found5,files_created
c
      call makefilename(tempdir,name1,name1)
      call makefilename(tempdir,name2,name2)
      call makefilename(tempdir,name3,name3)
      call makefilename(tempdir,name4,name4)
      call makefilename(tempdir,name5,name5)
c  
      do while (.not.files_created)
c
         call random_number(x)
         write(randstr,'(f12.8)') x
c  
c        Take last 8 digits of random number only
c
         lenr = lenstr(randstr)
         randstr='.'//randstr(lenr-8+1:lenr)
c
c        Create temporary file names 
c
         call makefilename(name1,randstr,tmpname1)
         call makefilename(name2,randstr,tmpname2)
         call makefilename(name3,randstr,tmpname3)
         call makefilename(name4,randstr,tmpname4)
         call makefilename(name5,randstr,tmpname5)
c         
         tlen1 = lenstr(tmpname1)  
         tlen2 = lenstr(tmpname2)  
         tlen3 = lenstr(tmpname3)  
         tlen4 = lenstr(tmpname4)  
         tlen5 = lenstr(tmpname5)  
c
c         write(0,*) 'TMP1:',tmpname1(1:tlen1)
c         write(0,*) 'TMP2:',tmpname2(1:tlen2)
c         write(0,*) 'TMP3:',tmpname3(1:tlen3)
c         write(0,*) 'TMP4:',tmpname4(1:tlen4)
c         write(0,*) 'TMP5:',tmpname5(1:tlen5)
c
         inquire(FILE=tmpname1(1:tlen1),EXIST=found1) 
         inquire(FILE=tmpname2(1:tlen2),EXIST=found2) 
         inquire(FILE=tmpname3(1:tlen3),EXIST=found3) 
         inquire(FILE=tmpname4(1:tlen4),EXIST=found4) 
         inquire(FILE=tmpname5(1:tlen5),EXIST=found5) 
c
         if ((.not.found1).and.(.not.found2).and.
     >       (.not.found3).and.(.not.found4).and.(.not.found5)) then
            call cissue(' touch '//tmpname5(1:tlen5),rc) 
            call cissue(' touch '//tmpname4(1:tlen4),rc) 
            call cissue(' touch '//tmpname3(1:tlen3),rc) 
            call cissue(' touch '//tmpname2(1:tlen2),rc) 
            call cissue(' touch '//tmpname1(1:tlen1),rc) 
            name1= tmpname1(1:tlen1)  
            name2= tmpname2(1:tlen2)  
            name3= tmpname3(1:tlen3)  
            name4= tmpname4(1:tlen4)  
            name5= tmpname5(1:tlen5)  
            files_created=.true.
         endif  
      end do
c
      return
      end 
c
c
c
      subroutine makelist(searchdir,files,files2,outputfile)
      implicit none
      character*(*) searchdir,files,outputfile,files2
c
c     Generates a list of the files specified by the file
c     search specification files in the directory searchdir and
c     places the results in the outputfile. 
c
c
      integer len,len1,len2,len3,len4,lenstr
      external lenstr 
      character*(512) command 
      integer rc
      character*50 err_redirect
c
c     Initilize
c     
      err_redirect = ' 2> /dev/null'
c     
c     Compile list of files to include in index 
c
      len1 = lenstr(searchdir)
      len2 = lenstr(files)
      len4 = lenstr(files2)    
      len3 = lenstr(outputfile) 
c
      call cissue('cd '//searchdir(1:len1)
     >        //' ; ls -l '//files(1:len2)
     >        //' > '//outputfile(1:len3)
     >        //err_redirect,rc)
c
      if (files2.ne.' ') then 

         call cissue('cd '//searchdir(1:len1)
     >        //' ; ls -l '//files2(1:len4)
     >        //' >> '//outputfile(1:len3)
     >        // err_redirect,rc)

      endif 
c
      return
      end
c
c
c
      subroutine getdata(inunit,pltunit,resdir,a,nfound,
     >                   tempplot,newentry,debug,ierr)
      use logrecorddata
      implicit none
c
      type (logrecord) newentry
c
c     This routine extracts the data required from the html file.  
c
      logical debug
      integer nfound,inunit,ierr,nplots,pltunit
      character*(*) resdir,a(nfound),tempplot
c
      character*100  casename,plotnames(maxrefs),seriesname,
     >               basename
      character*20   submitted,date,time 
      character*512  description
c
c     Local variables
c
      integer len,len1,len2,lenstr,extstr,findchars,pos1,in
      external lenstr,extstr,findchars
      character*256 line,templine
c
c     Line Parsing variables
c
      integer maxitems,items
      parameter(maxitems = 10)
      character*100  b(maxitems)
c
      logical found
c
c     Initialize 
c
c     Init record 
c
      newentry%seriesname = ' ' 
      newentry%casename = ' ' 
      newentry%submitted = ' ' 
      newentry%date = ' ' 
      newentry%time = ' ' 
      newentry%desc = ' ' 
      newentry%nrefs= 0 
      newentry%lastrecord = 0 
      newentry%nextrecord= 0 
c
c     Init local values
c
      casename = ' '
      submitted= ' '
      date     = ' '
      time     = ' '
      description = ' '  
c
c     Init some counters
c
      nplots = 0
      ierr = 0
c
c     Assign the case name as well as owner - extract the plot name and 
c     check if it exists.
c
      len = lenstr(a(nfound))
c
c     Check for common non-case files in the html directory 
c
      if (a(nfound)(1:len).eq.'index.html') then 
         ierr = 1
         return
      elseif (a(nfound)(1:len).eq.'index1.html') then 
         ierr = 1
         return
      elseif (a(nfound)(1:len).eq.'caselog.html') then 
         ierr = 1
         return
      elseif (a(nfound)(1:len).eq.'result_index.html') then 
         ierr = 1
         return
      endif
c
c     Set case name
c
      casename = a(nfound)(1:len)
c  
c     Determine the plot name and see if it exists
c
c      write(0,*) 'Casename:',casename
c
      pos1 = findchars(casename,'.html',-1,-1)
c
      if (pos1.eq.-1) then   
         pos1 = findchars(casename,'.dat',-1,-1)
      endif 
c
      if (pos1.gt.1) then 
         basename = casename(1:pos1-1) 
      else
         basename = ' '
      endif 
c
      if (debug) then 
         write(0,*) 'Basename:',basename,':'
      endif 
c
c     No case name found - all spaces 
c
      if (lenstr(basename).eq.1) then  
c
         write(0,*) 'ERROR: No name found:',pos1
         ierr = 1 
         return 
c
      endif
c
      call findplots(pltunit,resdir,basename(1:pos1-1),tempplot,
     >               plotnames,maxrefs,nplots)
c
c
c     Pull out the owner of the file for the submitted field 
c 
c     Name should be third item 
c
      len = lenstr(a(3))
      submitted = a(3)(1:len)
c
c      write(0,*) 'Submitted:',submitted(1:len)
c
c     Open file and extract remaining data
c
      len = lenstr(resdir)
c
c     The last entry in a is the file name 
c
      len1= lenstr(casename)
c
      open(inunit,FILE=resdir(1:len)//casename(1:len1),
     >     STATUS='OLD',IOSTAT=ierr) 
c     
      if (debug.and.ierr.eq.1) then 
         write(0,*) 'FAILED HTML OPEN:',
     >              resdir(1:len)//casename(1:len1)
      endif
c
      if (ierr.eq.0) then
c     
c        Read in the html file and pull out the case description,
c        date and time from the header.
c     
c     
c        Loop through file to header information (if found)
c     
         found = .false.  
c     
         do while(ierr.eq.0.and.(.not.found))
c     
            read(inunit,'(a256)',IOSTAT=ierr) line
c     
            len2 = extstr(line,len1)
c     
            if (len2.ne.0) then 
               if (line(len1:len1+10).eq.'***********')
     >              found=.true.
            endif 
c     
         end do
c     
c        Case header found
c     
         if (found) then 
c     
c           Continue to read and extract description, date and time
c     
c           Description is line 5 until next blank line 
c     
            do in = 1,4
               read(inunit,'(a256)',IOSTAT=ierr) line
            end do
c
            if (ierr.eq.0) then 
c     
c              Read in description 
c     
               description = ' ' 
c     
               do while (len2.ne.0.and.ierr.eq.0) 
c     
                  read(inunit,'(a256)',IOSTAT=ierr) line
c     
c                  write (0,*) len_trim(line)
c                  
                  call extracttext(line,templine,len2,ierr) 
c     
                  if (len2.ne.0) then 
c     
                     len2 = extstr(templine,len1)
                     len  = lenstr(description)
                     description = description(1:len)//
     >                             templine(len1:len2)   
      
                  endif
c     
               end do                 
c     
               len = lenstr(description)
c
c               write(0,*) 'DESC:',description(1:len)
c     
            endif 
c     
            if (ierr.eq.0) then 
c     
c              Read and discard more lines until line contains the 
c              word DIVIMP in capital letters then pull out the
c              date and time from that line
c
c              Needs updating for LIM result files
c     
c              Initialize
c
               date=' '
               time=' '
c
               read(inunit,'(a256)',IOSTAT=ierr) line
c
               do while (findchars(line,'DIVIMP',1,1).eq.-1.and.
     >                   ierr.eq.0) 
c
                  read(inunit,'(a256)',IOSTAT=ierr) line
c
               end do 
c     
               if (ierr.eq.0) then 
c     
c                  read(inunit,'(a256)',IOSTAT=ierr) line
c     
                  len2 = extstr(line,len1)
c      
                  if (len2.gt.3.and.ierr.eq.0) then 
c     
c                    Remove *'s
c     
                     templine = line(len1+1:len2-1)  
c     
                     call parseline(templine(len1:len2),b,
     >                                 maxitems,items)              
c     
                     if (items.ge.2) then 
c     
                        len = lenstr(b(1))
                        date = b(1)(1:len)
c
c                        write(0,*) 'DATE:',date(1:len)
c
                        len = lenstr(b(2))
                        time = b(2)(1:len)
c
c                        write(0,*) 'TIME:',time(1:len)
c     
      
                     endif
                  endif
               endif 
c
c              Since this is not a fatal error reset ierr to 0
c
               ierr =0
 
            endif
         else
            len= lenstr(resdir)
            len1= lenstr(casename)
            write(0,*) 'WARNING: ',
     >                resdir(1:len)//casename(1:len1),
     >               ' IS NOT A VALID DIVIMP HTML FILE'
         endif
c     
c        Close file 
c     
         close(inunit)
c     
      endif
c
c     Find series name 
c
c
      pos1 = findchars(basename,'-',-1,-1)
c
      if (pos1.le.1) then 
         seriesname = basename
      else
         seriesname = basename(1:pos1-1)
      endif
c
c
c     Save data to entry 
c
      newentry%seriesname = seriesname
      newentry%casename = basename
      newentry%submitted = submitted
      newentry%date = date
      newentry%time = time
      newentry%desc = description
      newentry%nrefs= nplots + 1
c
      newentry%refs(1) = casename 
c
      do in = 1,nplots
         newentry%refs(in+1) = plotnames(in)
      end do 
c
      newentry%lastrecord = 0 
      newentry%nextrecord= 0 
c
      return
      end 
c
c
c
      subroutine findplots(pltunit,resdir,basename,tempplot,
     >                     plotnames,maxplots,nplots)
      implicit none
c
      integer pltunit,maxplots,nplots
      character*(*) resdir,basename,tempplot,plotnames(maxplots)
c
c     This routine searches the target directory for plots that are 
c     a variation of the basename.
c
      integer len,len1,len2,len3,lenstr,rc,ios,ierr
      external lenstr
c
c     Line parsing variable 
c
      integer maxitems,items
      parameter(maxitems = 10)
      character*100  b(maxitems)
      character*256 line
      character*50 err_redirect
c
c     Initialization
c
      err_redirect = ' 2> /dev/null'
      nplots = 0 
      ierr = 0
      ios = 0
c
c     Issue commands to fill the tempplot file with requisite file names
c      
      len1 = lenstr(resdir)
      len2 = lenstr(basename)
      len3 = lenstr(tempplot)  
c
c     .notes 
c
      call cissue('cd '//resdir(1:len1)//' ; '
     >     //'ls -l '//basename(1:len2)//'.notes > '
     >     //tempplot(1:len3)
     >     //err_redirect,rc)
c
c     append .pdf 
c
      call cissue('cd '//resdir(1:len1)//' ; '
     >     //'ls -l '//basename(1:len2)//'.pdf >> '
     >     //tempplot(1:len3)
     >     //err_redirect,rc)
c
c     append .*.pdf
c
      call cissue('cd '//resdir(1:len1)//' ; '
     >     //'ls -l '//basename(1:len2)//'.*.pdf >> '
     >     //tempplot(1:len3)
     >     //err_redirect,rc)
c
c     append -*.pdf
c
      call cissue('cd '//resdir(1:len1)//' ; '
     >     //'ls -l '//basename(1:len2)//'-*.pdf >> '
     >     //tempplot(1:len3)
     >     //err_redirect,rc)
c
c     append .jpg
c
      call cissue('cd '//resdir(1:len1)//' ; '
     >     //'ls -l '//basename(1:len2)//'.jpg >> '
     >     //tempplot(1:len3)
     >     //err_redirect,rc)
c
c     append .*.jpg
c
      call cissue('cd '//resdir(1:len1)//' ; '
     >     //'ls -l '//basename(1:len2)//'.*.jpg >> '
     >     //tempplot(1:len3)
     >     //err_redirect,rc)
c
c     append -*.jpg
c
      call cissue('cd '//resdir(1:len1)//' ; '
     >     //'ls -l '//basename(1:len2)//'-*.jpg >> '
     >     //tempplot(1:len3)
     >     //err_redirect,rc)
c
c     Loop through tempplot and add plotnames to list - 10 maximum
c
      open(pltunit,FILE=tempplot,IOSTAT=ios) 
c
      if (ios.eq.0) then 
c
         do while(ierr.eq.0)               
c            
            read(pltunit,'(a256)',IOSTAT=ierr) line
c
            if (ierr.eq.0) then 
c
c              Parseline
c
               len = lenstr(line)
               call parseline(line(1:len),b,maxitems,items)              
c
               if (items.gt.0) then 
c
                  nplots  = nplots+1 
c   
                  len = lenstr(b(items))
c   
                  plotnames(nplots) = b(items)(1:len)
c
                  if (nplots.eq.maxplots) ierr =2 
c
               endif   
c
            endif
c      
         end do
c
      end if  
c
      close (pltunit)
c
      return 
      end 
c
c
c
      subroutine extracttext(line,templine,length,ierr) 
      implicit none
      integer length,ierr  
      character*(*) line, templine
c
      integer len2,len1,len,lenstr,extstr
c
      len2 = extstr(line,len1)
c
      if (len2.lt.3) then 
c
         length = 0 
         templine = ' '
c
         return
c             
      endif  
c
      templine = line(len1+1:len2-1)
c  
      len2 = extstr(templine,len1)
c
      if (len2.eq.0) then
c
         length = 0 
         templine = ' '
         return
c
      endif
c
      templine = templine(len1:len2)
c
      length = lenstr(templine)
c
      return
      end  
c
c
c
      subroutine parseline (line,a,maxn,nfound)
      implicit none
      integer maxn, nfound
      character*(*) line
      character*(*) a(maxn)
c
c     PARSELINE: This routine parses the input line into its 
c                component strings.
c
c     Local variables 
c        
      integer len1,len2,in,ic,lenstr,extstr 
      integer pos1,pos2,findchars
      external lenstr,extstr,findchars
      character*512 templine 
c
      len1 = lenstr(line) 
c
      nfound = 0     
      pos1   = 0
c      
      templine = line
c
      len2 = extstr(templine,len1)
c
c     Check for blank string  
c
      if (len2.eq.0) then 
         do in = 1,maxn
            a(in) = ' '
         end do
         nfound = 0
         return
      endif     
c
c     Loop through string assigning space separated substrings -
c     eliminate extra spaces - start at beginning and work to the end. 
c
      do while (len2.ne.0.and.nfound.lt.maxn.and.pos1.ne.-1) 

         templine = templine(len1:len2)
         len1 = lenstr(templine) 
c
c        Find next space separator
c
         pos1 = findchars(templine(1:len1),' ',1,1) 
c
         if (pos1.gt.0) then 
c
c           Record substring
c 
            nfound = nfound +1 
c
            a(nfound) = templine(1:min(pos1-1,len(a(nfound))))
c
c           Remove string from start of templine and continue parsing 
c
            templine = templine(pos1+1:len2)
            len2 = extstr(templine,len1)
c
         elseif (pos1.eq.-1.and.len1.gt.0) then
c
c           Record last substring
c 
            nfound = nfound +1 
c
            a(nfound) = templine(1:len1)
c
c           Allow for exit by leaving pos1 as -1
c
         endif 
c
      end do
c
c
c
      return
      end

C
C  *********************************************************************
C  *  FINDCHARS: Search for characters                                 *
C  *********************************************************************
C
      integer function findchars(source,target,istart,idir)
      implicit none
      character*(*) source,target
      integer istart,idir
c
c     Routine searches for a substring defined by target inside the  
c     string defined by source. It starts at position start and proceeds
c     either forward or backward depending on the value of dir. 
c
c     start = start position for searching
c           =  1 = start at beginning
c           = -1 = start at end
c
c     dir  >0 = search forward
c     dir  <0 = search backward
c     dir  =0 = look for exact match only
c
c     Routine returns the starting position of match or -1 if no match
c
      integer start,dir 
      integer targlen,in,sourcelen,lenstr
      external lenstr
c
      start = istart
      dir = idir
c
      targlen = lenstr(target)
      sourcelen = lenstr(source)
c
      findchars = -1
c
      if (targlen.eq.0) return
c
c     Get length of source string for start = -1
c
      if (start.eq.-1) start = sourcelen 
c
      if (dir.eq.0) then 
c
         if (target(1:targlen).eq.source(1:sourcelen)) findchars = 1
c   
c     Search forward
c
      elseif (dir.gt.0) then  
c
         do in = start,sourcelen-targlen+1
c
            if (target(1:targlen).eq.source(in:in+targlen-1)) then
               findchars = in 
               return
            endif
c                           
         end do
c
c     Search backward
c
      elseif (dir.lt.0) then 
c
         do in = start,targlen,-1
c
            if (target(1:targlen).eq.source(in-targlen+1:in)) then
               findchars = in -targlen+1
               return
            endif
c                           
         end do
c
      endif
c
      return
      end 
C
C
C
      INTEGER FUNCTION LENSTR (ASTR)
      implicit none
      CHARACTER*(*) ASTR
C
C  *********************************************************************
C  *                                                                   *
C  *  LENSTR: RETURNS EFFECTIVE LENGTH OF STRING ASTR IGNORING         *
C  *          ANY TRAILING BLANKS.                                     *
C  *                                                                   *
C  *********************************************************************
C
      integer i
c
      DO 10 I = LEN(ASTR),1,-1
         IF (ASTR(I:I) .NE. ' ') THEN
            LENSTR = I
            RETURN
         ENDIF
   10 CONTINUE
      LENSTR = 1
      RETURN
      END
C
C
C
      INTEGER FUNCTION EXTSTR (ASTR,start)
      implicit none
      CHARACTER*(*) ASTR
      integer start,i
c
      extstr = 0
C
C  *********************************************************************
C  *                                                                   *
C  *  EXTSTR: RETURNS EFFECTIVE LENGTH OF STRING ASTR IGNORING         *
C  *          ANY TRAILING BLANKS AND THE EFFECTIVE STARTING POINT BY  *
c  *          IGNORING LEADING BLANKS.                                 *
C  *                                                                   *
C  *********************************************************************
C
      DO 10 I = LEN(ASTR),1,-1
         IF (ASTR(I:I) .NE. ' ') THEN
            extSTR = I
            goto 20
         ENDIF
   10 CONTINUE

 20   if (extstr.gt.0) then
         do i = 1, extstr
            if (astr(I:I).ne.' ') then
               start = i
               return
            endif
         end do
      else
         extSTR = 0
         start  = 0
      endif
c
      RETURN
      END
C
C
C
      SUBROUTINE CISSUE(command,code)
      implicit none
      CHARACTER command*(*)
      INTEGER code

      INTEGER System

      INTEGER len,lenstr
      external lenstr,system
      CHARACTER exeline*256

      len = lenstr(command)
c
      exeline = command(1:len)
c
c     Add NULL termination required for 'C' strings
c
      exeline(len+1:len+2) = '\0'
c
c      WRITE(0,*) 'CISSUE: "'//exeline(1:len+1)//'"'
c
      code = SYSTEM(exeline(1:len+1))

      RETURN
      END



C
C  *********************************************************************
C  *  OPENHTML:  ADDS HTML HEADER INFORMATION                          *
C  *********************************************************************
C
      SUBROUTINE openhtml(htmlunit,filename,title)
      implicit none
c
      integer htmlunit,len,lenstr,len1
      character*(*) filename,title
c
      len  = lenstr(filename)
      open(htmlunit,FILE=filename(1:len))  
c            
      rewind(htmlunit)
c  
      write (htmlunit,*) '<HTML>'
      write (htmlunit,*) '<HEAD>'
c
      len1 = lenstr(title)
      write (htmlunit,*) '  <TITLE>'//title(1:len1)//'</TITLE>'
c
      write (htmlunit,*) '   <META NAME="DIVIMP DAT File" CONTENT="">'
      write (htmlunit,*) '  <STYLE TYPE="text/css">'
      write (htmlunit,*) '    BODY     { color:black; '
      write (htmlunit,*) '               background:white }'
      write (htmlunit,*) '    A:link   { color:blue; '
      write (htmlunit,*) '               text-decoration: none;'
      write (htmlunit,*) '               font-weight:normal}'
      write (htmlunit,*) '    A:active { color:green; '
      write (htmlunit,*) '               text-decoration: none;'
      write (htmlunit,*) '               font-weight:bold}'
      write (htmlunit,*) '    A:visited{ color:red; '
      write (htmlunit,*) '               text-decoration: none;'
      write (htmlunit,*) '               font-weight:normal}'
      write (htmlunit,*) 
      write (htmlunit,*) '    H1       { font-weight:bold} '
      write (htmlunit,*) 
      write (htmlunit,*) '    H2       { font-weight:bold} '
      write (htmlunit,*) 
      write (htmlunit,*) '    H3       { font-weight:bold; '
      write (htmlunit,*) '               color:crimson} ' 
      write (htmlunit,*) '</STYLE>'  
      write (htmlunit,*) '</HEAD>'
      write (htmlunit,*) '<BODY>'
c
c
      RETURN
      END

C
C  *********************************************************************
C  *  CLOSEHTML:  ADDS HTML END OF DOCUMENT TAGS                       *
C  *********************************************************************
C
      SUBROUTINE closehtml(htmlunit)
      implicit none
      integer htmlunit
c
c      write(htmlunit,*) '</PRE>'
c 
      write(htmlunit,*) '</BODY>'
      write(htmlunit,*) '</html>'
c
      close(htmlunit) 
c
      RETURN
      END



      SUBROUTINE tablestart(htmlunit,tabletitle)
      implicit none
      integer htmlunit
      character*(*) tabletitle
      integer len,lenstr
      external lenstr 
c
      len = lenstr(tabletitle) 
c
      write (htmlunit,*) '<TABLE BORDER=2>'
      write (htmlunit,*) '<CAPTION>'
      write (htmlunit,*) '<P>&nbsp</P>'
      write (htmlunit,*) ' <H1>'//tabletitle(1:len)//'</H1>'
      write (htmlunit,*) '</CAPTION>'
      write (htmlunit,*) '<TBODY>'
      write (htmlunit,*) '<TR>'
      write (htmlunit,*) '   <TH>'
      write (htmlunit,*) '      Case Name'
      write (htmlunit,*) '   </TH>'
      write (htmlunit,*) '   <TH>'
      write (htmlunit,*) '      Submitted'
      write (htmlunit,*) '   </TH>'
      write (htmlunit,*) '   <TH>'
      write (htmlunit,*) '      Date'
      write (htmlunit,*) '   </TH>'
      write (htmlunit,*) '   <TH COLSPAN=3>'
      write (htmlunit,*) '      Description'
      write (htmlunit,*) '   </TH>'
      write (htmlunit,*) '</TR>'
c
      RETURN
      END
c
c
c
      SUBROUTINE tableend(htmlunit)
      implicit none
c
      integer htmlunit
c
      write (htmlunit,*) '</TABLE>'
c
      end
c
c
c
      SUBROUTINE tableentry(htmlunit,item)
      use logrecorddata
      implicit none
c
      integer htmlunit
c      include 'logrecord' 
      type (logrecord) item
c
      integer len,lenstr,in
c
      write (htmlunit,*) '<TR>'
      write (htmlunit,*) '   <TD>'
c
c      write(0,*) 'Adding Table entry:'
c
c
c     Case title and Hypertext References
c
      write (htmlunit,*) '   <UL>'
c
      len = lenstr(item%casename)
c
      write (htmlunit,*) '<LI><B>'//
     >                     item%casename(1:len)//'</B>'
c
      if (item%nrefs.gt.0) then 
c
          do in = 1,item%nrefs

             len = lenstr(item%refs(in))

             if (len.gt.1) then 
                   write (htmlunit,*)  
     >             '<LI><A HREF="'//item%refs(in)(1:len)//'">'//
     >                     item%refs(in)(1:len)//'</A>'
             endif
          end do
c
      endif 
c
      write (htmlunit,*) '   </UL>'
c
      write (htmlunit,*) '   </TD>'
      write (htmlunit,*) '   <TD>'
      len = lenstr(item%submitted) 
      write (htmlunit,*) item%submitted(1:len)
      write (htmlunit,*) '   </TD>'
      write (htmlunit,*) '   <TD>'
      write (htmlunit,*) '   <UL>'
      len = lenstr(item%date)    
      write (htmlunit,*) '<LI>'//item%date(1:len)
      len = lenstr(item%time)    
      write (htmlunit,*) '<LI>'//item%time(1:len)
      write (htmlunit,*) '   </UL>'
      write (htmlunit,*) '   </TD>'
      write (htmlunit,*) '   <TD COLSPAN=3>'
      len = lenstr(item%desc)
      write (htmlunit,*) item%desc(1:len) 
      write (htmlunit,*) '   </TD>'
      write (htmlunit,*) '</TR>'
c
      RETURN
      END










