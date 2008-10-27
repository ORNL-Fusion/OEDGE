      module logrecorddata 

         integer maxrefs
         type logrecord 
            character*3   tag 
            character*256 htmlref
            character*256 bookmark  
            character*512 desc
            integer nextrecord 
            integer lastrecord  
         end type logrecord

      end module logrecorddata   
c
c
c
      program buildtagdb
      use logrecorddata 
      implicit none
      integer maxfiles
c     
      parameter(maxfiles=3)
c
c     Locations ... and file names
c
      character*256 refnames(maxfiles),dbname,template
c
c     Table entry data 
c
c      include 'logrecord' 
c
      type (logrecord) newentry
c
c     Indices and unit numbers
c
      integer ios,ifile
      integer len,lenstr
      external lenstr,iargc
      integer nargs,iargc
      character*(256) arg  
c
      logical do_xref
c
      integer tmpunit,inunit,outunit,dbunit
c
      data inunit  /9/
      data outunit /10/
c
      data dbunit  /13/
c
c     Debugging
c
      logical debug  
c
c     Initilaization 
c
      debug =  .false.
      do_xref = .false. 
c
      write(0,'(a)') 'Building TAG database' 
c
      call initdb  
c
c     Fixed file names
c
      refnames(1) = 'divref.html'
      refnames(2) = 'divref-sol23.html'
      refnames(3) = 'divref-opt.html'
      dbname='tagdb.db' 
c
c     Template is initially empty.
c
      template=' '
c
c     Temporary file names and location
c
c      tempdir   = '/tmp/'
c
c     Read in any command line arguments
c
      nargs = iargc()
c
      if (nargs.eq.1) then
c
c        Template is first argument    
c
         call getarg(1,arg)
c
         len = lenstr(arg)
         template = arg(1:len)
c
         do_xref = .true.
c
      endif 
c
c     Loop through files list to extract tags and generate the 
c     tag database.
c
      ifile = 1
c
      do while (ios.eq.0.and.ifile.le.maxfiles) 
c
c
c        Loop through list of html files extracting tag references
c        and adding them to the database.
c
         open(inunit,FILE=refnames(ifile),
     >        STATUS='OLD',IOSTAT=ios) 
c
         len = lenstr(refnames(ifile))
         write(0,'(a,a)') 'Searching file: ',refnames(ifile)(1:len)
c
         if (ios.eq.0) then  
            call gettags(inunit,refnames(ifile))
         endif  
c
         close(inunit) 
c
         ifile = ifile +1 
c
      end do
c
c     Database of collected tags should be complete at this point.
c
c     Write the database to a file.
c
      call writedb(dbunit,dbname)
c
c     If a template file has been specified then apply the collected 
c     database information to the template and create an HTML 
c     cross-reference file based on the database contents.
c
      if (do_xref) then 
         call create_xref(template,inunit,outunit)
      endif 
c
c     End message
c
      stop 'Finished Builing Tag Database' 
c
      end
c
c    
c
      subroutine create_xref(template,inunit,outunit)
      implicit none
      character*(*) template
      integer inunit,outunit
c
c     This routine will read in a template file and 
c     create a new file with a .html extension which 
c     will contain hypertext references for each TAG in 
c     the file. 
c
c     The name of the output file will be the same as the 
c     input except a .html extension will be added. 
c      
      integer len,lenstr,ios
      external lenstr
      character*512 line
c
c     Open the input file 
c
      len = lenstr(template)
c
      open(inunit,FILE=template(1:len),
     >        STATUS='OLD',IOSTAT=ios) 
c
c     If the open on the input file is successful then continue
c
      if (ios.eq.0) then
c
         write(0,'(a,a)') 'Creating HTML Cross-reference from: ',
     >                  template(1:len)

c
c        Open the html output file 
c
         call openhtml(outunit,template(1:len)//'.html',
     >                 'CROSS-REFERENCED DIVIMP INPUT FILE')

         do while (ios.eq.0) 
c
            read(inunit,'(a512)',IOSTAT=ios) line
c
            call refhtml(outunit,line) 
c
         end do 

         call closehtml(outunit) 

      endif 
c
      return
      end      

C
C  *********************************************************************
C  *  REFHTML:  ADDS HTML HYPERTEXT REFERENCE TO LINE OF FILE          *
C  *********************************************************************
C
      SUBROUTINE refhtml(htmlunit,line)
      use logrecorddata
c
      implicit none
      character*(*) line 
c
      integer htmlunit
c
      character*3 tag 
      logical new,found
      integer in,len,lenh,lenb,lenstr,getindex
      external lenstr,getindex
c
      type (logrecord) item
c
      call extracttag(line,tag,found)
c
      if (found) then 
c
         in =  getindex(tag,new)
c
         if (new) then 
c
            len = lenstr(line)
            write(htmlunit,'(a)') line(1:len)
c
         else
c
            call getrecord(in,item)
c
            len = lenstr(line) 
            lenh = lenstr(item%htmlref)
            lenb = lenstr(item%bookmark)
c
            write(htmlunit,'(a,a,a)') 
     >       '<A HREF="'//
     >        item%htmlref(1:lenh)//item%bookmark(1:lenb)//'">',
     >        line(1:len), 
     >        '</A>'
c
         endif
c
      else
c
         len = lenstr(line)
         write(htmlunit,'(a)') line(1:len)
c
      endif 
c
      RETURN
      END
c
c
c
      subroutine extracttag(line,tag,found)
      implicit none
      character*(*) line,tag
      logical found
c
c     This routine looks for a TAG in the appropriate place
c     in the text string.
c 
      character*3 tmptag
      integer lenstr,len,len1,len2,extstr
      external lenstr,extstr
c
c     Initialize to false
c
      found = .false.
c
      len = lenstr(line)
      tag = ' '
c
      if (len.lt.5) return
c
c     Does the line start with a single quote "'" or a "c"?
c     Is the second character on the line a "+" or a "*"?
c     Are the next three characters after the "+" or "*" not spaces?
c
      if (line(1:1).eq.''''.or.line(1:1).eq.'c'.or.
     >    line(1:1).eq.'C') then
c
         if (line(2:2).eq.'+'.or.line(2:2).eq.'*') then 
c
            tmptag = line(3:5)
c
            len2 = extstr(tmptag,len1)
c
            if (len1.eq.1.and.len2.eq.3) then 
c 
               found = .true.
               tag = tmptag
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
      subroutine gettags(inunit,htmlref1)
      use logrecorddata
c
      implicit none
c
c     This routine extracts the tags required from the html file.  
c
      logical debug
      integer inunit
      character*(*) htmlref1
c
c
c     Local variables
c
      type (logrecord) newentry
c
      integer len,len1,len2,lenstr,extstr,findchars,len3
      integer in1,in2,ios,pos1,pos2,lenh,extstr
      external lenstr,extstr,findchars,extstr
      character*512 line,templine
      character*512 tag
      logical found,done
c
c     Line Parsing variables
c
c     Initialize 
c
      found =.false.
      newentry%tag = ' '
      newentry%htmlref=' '
      newentry%bookmark=' '
      newentry%desc=' '
c
c     Rewind the input file just in case
c
      rewind(inunit)
c
c     Scan for the <DIV ID=toc> tag 
c
      ios = 0
c
      do while (ios.eq.0.and..not.found)
c
         read(inunit,'(a512)',END=100,ERR=200,IOSTAT=ios) line
c
         if (findchars(line,'<DIV ID=toc>',1,1).gt.0) found = .true.
c
      end do   
c
c     Start looking for tags 
c
c     Tags are preceded by the specific TAG <FONT CLASS=c1> 
c
c     TAGS are expected with ONLY one on each line and the Hypertext reference
c     or bookmark is the immediately preceding tag enclosed
c     in " (double quotes)
c
      done = .false.
c
      do while (ios.eq.0.and..not.done)  
c
         read(inunit,'(a512)',END=100,ERR=200,IOSTAT=ios) line
c
         len = lenstr(line) 
c
c        Found mark up for tag
c
         if (findchars(line,'<FONT CLASS=c1>',1,1).gt.0) then 
c
c           Find start position of TAG 
c
            in1 = findchars(line,'<FONT CLASS=c1>',1,1)
c
c           Find END position of TAG string
c
            in2 = findchars(line,'</',in1+15,1)
c
c           Assign TAG String           
c
            tag = line(in1+15:in2-1)
c
c           Assign name of html reference document
c
            lenh = lenstr(htmlref1)
            newentry%htmlref = htmlref1(1:lenh)
c
c           Find bookmark - bookmarks are enclosed in " " - 
c           immediately preceding the ID TAG.
c
            pos2 = findchars(line,'"',in1-1,-1)
c
            pos1 = findchars(line,'"',pos2-1,-1)
c
            newentry%bookmark = line(pos1+1:pos2-1)
c
c           Extract Description from TOC as well - assuming that the 
c           entry is all on one line in the html file. Remove the end
c           Anchor tag at the end - otherwise - it should ust be the 
c           entire line after the TAG section.
c
            len1 = lenstr(line) 
c
            templine = line(in2+7:len1-4) 
c
            len3 = extstr(templine,len2)
c
            newentry%desc = templine(len2:len3)
c
c           Interpret the TAG string and add TAGS to the database.  
c
            call addtags(tag,newentry)
c
         endif 
c
c        Check for end condition 
c
         if (findchars(line,'</DIV>',1,1).gt.0) done=.true.
c
      end do          
c
c     Normal exit
c
      return

 100  len = lenstr(htmlref1)
      close(inunit)
      write(0,'(a,a)') 'END of File reached on :',htmlref1(1:len)
      write(0,'(a,i5)') 'ERROR Number Reported = ',ios
      return

 200  len = lenstr(htmlref1)
      close(inunit)
      write(0,'(a,a)') 'ERROR Encountered on File :',htmlref1(1:len)
      write(0,'(a,i5)') 'ERROR Number Reported    = ',ios
      return
c
c
      return 
      end
c
c
c
      subroutine addtags(tag,newentry)
      use logrecorddata
c
      implicit none
      character*(*) tag 
      type (logrecord) newentry
c
c     This routine interprets the TAG string and adds any 
c     tags found there to the database - most TAG strings
c     represent one identifier but some represent more. 
c
c
      integer in,in1,in2
      integer lent
      integer findchars,lenstr
      external findchars,lenstr 
      logical done
c
      lent = lenstr(tag)
c
      done = .false.
      in1  = 1 
c
      do while (.not.done)
c
         in = findchars(tag,',',in1,1)
c
c        TAG still contains a comma - process and continue
c
         if (in.gt.0) then 
c
            call proctag(tag(in1:in-1),newentry)
c
            in1 = in+1 
c
c        Adding last TAG - set done flag
c
         else
c
            call proctag(tag(in1:lent),newentry)

            done=.true.
c
         endif
c
      end do
c
      return
      end
c
c
c
      subroutine proctag(tag,newentry)
      use logrecorddata
c
      implicit none
      character*(*) tag
      type (logrecord) newentry
c
c     This routine interprets the TAG string and adds any 
c     tags found there to the database - most TAG strings
c     represent one identifier but some represent more. 
c
c
      integer in,in1,in2
      integer lent
      integer findchars,lenstr
      external findchars,lenstr 
      integer tagstart,tagend,tagcnt
      character*3 tmptag,tag1,tag2
      character*1 taglet 
      logical done
c
      lent = lenstr(tag)
c
      done = .false.
      in1  = 1 
c
      in = findchars(tag,'..',in1,1)
c
c     Add multiple consecutive tags 
c
      if (in.gt.0) then
c
         tag1=tag(1:in-1)
         tag2=tag(in+2:lent)
c
c        Read tag series letter 
c
         taglet=tag(1:1)
c
c        Read Tag number indices.
c
c        NOTE: This is assuming that the second two digits of the
c              TAG are always NUMERIC - this may not always be true.
c
         read(tag1(2:3),'(i2)') tagstart
         read(tag2(2:3),'(i2)') tagend
c
c        Loop through adding all tags to database
c     
         do tagcnt = tagstart,tagend
c
            if (tagcnt.lt.10) then
c
               write(tmptag,'(a1,''0'',i1)') taglet,tagcnt
c
            else 
c
               write(tmptag,'(a1,i2)') taglet,tagcnt
c
            endif
c
            newentry%tag = tmptag
c
            call addrecord(newentry) 
c
         end do
c
c     Add single TAG
c 
      else
c
         newentry%tag=tag(1:lent)
c 
         call addrecord(newentry) 
c
      endif
c
      return
      end
c
c
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
      write (htmlunit,*) '<PRE>'
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
      write(htmlunit,*) '</PRE>'
c 
      write(htmlunit,*) '</BODY>'
      write(htmlunit,*) '</html>'
c
      close(htmlunit) 
c
      RETURN
      END

