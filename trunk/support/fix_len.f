      program fix_len
      implicit none
c
      integer copyfiles
      parameter(copyfiles=1)
c
      
      integer ios,cnt,flen,llen,eof,fnum,fixed,total_fixed
      character*100 filename
      character*512 line
      integer lenstr,findchars
      external lenstr,findchars
      character*200 exeline
      integer elen,system,retcode
      external system

      open(10,file="linelength_output.txt")

      open(8,file="filelist")

      ios=0
      fnum = 11 
      total_fixed = 0

      do while (ios.eq.0)    

         read(8,'(a100)',iostat=ios) filename

         if (ios.eq.0) then 

            flen=lenstr(filename)

            if (flen.le.1) cycle

            open(9,file=filename(1:flen))
            open(fnum,file=filename(1:flen)//'.mod')

            eof=0
            cnt = 0

            do while (eof.eq.0)

               read(9,'(a512)',iostat=eof) line

               if (eof.eq.0) then

                  call rewriteline(line,fnum,fixed) 
                  cnt = cnt+fixed

               endif
            end do 
 
            close(9)
            close(fnum)

            write(6,'(a,a,a,i6,a)') 
     >              'FINISHED File:', filename(1:max(flen,30)),
     >              '   Fixed ',
     >               cnt,' lines '
            total_fixed = total_fixed+ cnt
c
c           If lines were fixed - copy the new file in place of the original
c
            if (cnt.gt.0.and.copyfiles.eq.1) then 
               exeline = 'cp '//filename(1:flen)//'.mod '
     >                        //filename(1:flen)
               elen=lenstr(exeline)
               retcode=system(exeline(1:elen))
            endif
c
             
         end if

      end do 

      write(6,'(a,i6,a)') 
     >              'FINISHED:', total_fixed,
     >             ' long lines were repaired.'




      close(8)
      close(10)

      return
      end
      

      subroutine rewriteline(line,fnum,fixed)
      implicit none
      integer fnum,fixed
      character*(*) line
      integer extstr,lenstr,findchars
      external extstr,lenstr,findchars
c
      integer llen,start,charend,endpos,comloc,endpos2,llen2
c
      character*100 endchars
c
      logical is_comment
c
      endchars = '+-*/() ,:;'
c
      is_comment=.false.
      fixed= 0
      charend=lenstr(endchars)
c
      llen=extstr(line,start)
      llen2=lenstr(line)
c
      if (llen.ne.llen2) then 
         write (0,*) 'LENGTH ERROR: lenstr,extstr',llen,llen2
      endif
c
      comloc = findchars(line(1:llen),'!',-1,-1)

c
      if (line(start:start).eq.'!') is_comment=.true. 

c
c     Break and reformat long lines in the source code
c     Exclude likely comment lines
c 
 
      if (llen.gt.132.and.(comloc.lt.1.or.comloc.gt.132).and.
     >      .not.is_comment) then 
c
c        First - remove inline spaces to see if this is 
c        sufficient
c
         call remove_spaces(line,start,llen,comloc,fixed)
c
c        If the line has been fixed - write it out and exit 
c
         if (fixed.eq.1) then 
            write(fnum,'(a)') line(1:llen)
         else
c
            fixed = 1 
c         
            endpos = scan(line(1:131),endchars(1:charend),.true.)
c
c           write out first section
c
            write(fnum,'(a,a1)') line(1:endpos),'&'
c
c           Assume that no lines are greater than 260 
c           Also do not indent the new line since 
c           if the previous cut was in a character string
c           then any additional indenting will affect 
c           formatting - need a more intelligent code than 
c           this one for that
c
c           Need to rewrite the remaining portion of the line again
c
            if ((llen-endpos).gt.130) then

               endpos2 = scan(line(1:endpos+130),
     >                           endchars(1:charend),.true.)
c
c               write(6,'(a,4i6)') 'Endpos2:',endpos,endpos2,
c     >                      llen
c               write(6,'(a)') line(1:endpos)
c               write(6,'(a)') line(endpos+1:endpos2)
c               write(6,'(a)') line(endpos2+1:llen)
c               write(6,'(a)') line(endpos2+1:llen+10)
c
c
c
               write(fnum,'(a1,a,a1)') '&',line(endpos+1:endpos2),
     >                                 '&'
               write(fnum,'(a1,a)') '&',line(endpos2+1:llen)
c
            else

               write(fnum,'(a1,a)') '&',line(endpos+1:llen)
 
            endif
c
         endif
c
c     Lines less than 132 characters are fine
c
      else
 
         write(fnum,'(a)') line(1:llen)

      endif
c
c
c1
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
c
c
C
C  *********************************************************************
C  *  FINDCHARS: Search for characters                                 *
C  *********************************************************************
C
      integer function findchars(source,target,startpos,dir)
      implicit none
      character*(*) source,target
      integer startpos,dir
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
      integer targlen,start,in
c
c     Initialization 
c
      start = startpos
      targlen = len(target)
c
      findchars = -1
c
      if (targlen.eq.0) return
c
c     Get length of source string for start = -1
c
      if (start.eq.-1) start = len(source)  
c
      if (dir.eq.0) then 
c
         if (target.eq.source) findchars = 1
c   
c     Search forward
c
      elseif (dir.gt.0) then  
c
         do in = start,len(source)-targlen+1
c
            if (target.eq.source(in:in+targlen-1)) then
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
            if (target.eq.source(in-targlen+1:in)) then
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

 20   if (extstr.gt.1) then
         do i = 1, extstr
            if (astr(I:I).ne.' ') then
               start = i
               return
            endif
         end do
      else
         extSTR = 1
         start  = 1
      endif
c
      RETURN
      END
c
c
c
      subroutine remove_spaces(line,start,llen,comloc,fixed)
      implicit none
      character*(*) line
      integer start,llen,comloc,fixed
c
      character*512 newline
c
      integer pa,pb,in,np,lenstr,newlen
      integer p1,p2,p3,p4,p5,p6,endp
      external lenstr
      logical pa_alpha,pb_alpha
c
c     Remove spaces from the line to see if that will make it short enough
c     try not to remove those separating alphabetic characters
c
      p1 = ichar('a')
      p2 = ichar('z')
      p3 = ichar('A')
      p4 = ichar('Z')
      p5 = ichar('0')
      p6 = ichar('9') 
c
      newline=' '
c
      fixed = 0
c
      np = start-1
c
      if (comloc.gt.1) then 
         endp = comloc-1
      else
         endp = llen
      endif
c
      do in=start,endp
c
         if (line(in:in).ne.' ') then 
            np=np+1
            newline(np:np) = line(in:in)
c
c        Handle spaces
c
         else
c
            pa=ichar(line(max(in-1,start):max(in-1,start)))
            pb=ichar(line(min(in+1,llen):min(in+1,llen)))
c
c           Remove the space only if one of the two 
c           characters surrounding it is non-alphabetic and non-numeric
c
            pa_alpha= (pa.ge.p1.and.pa.le.p2).or.
     >                (pa.ge.p3.and.pa.le.p4).or.
     >                (pa.ge.p5.and.pa.le.p6)

            pb_alpha= (pb.ge.p1.and.pb.le.p2).or.
     >                (pb.ge.p3.and.pb.le.p4).or.
     >                (pb.ge.p5.and.pb.le.p6)
c
c           Copy space if both alphabetic
c 
            if (pa_alpha.and.pb_alpha) then 
               np = np+1
               newline(np:np) = line(in:in)
            elseif(pa_alpha.and.line(in+1:in+1).eq.' ') then 
               np = np+1
               newline(np:np) = line(in:in)
            endif
c
         endif
c
      end do
c
c
c     If the problem isn't fixed - leave the spaces in for readability - don't change the inputs
c
      if (np.le.132) then
         fixed = 1
         if (comloc.gt.1) comloc = np+1
c
c        If endp is less than llen - copy the rest - should be a comment
c
         if (endp.lt.llen) then 
c
            do in = endp+1,llen
               np = np+1
               newline(np:np)=line(in:in)
            end do
         endif
c
         line=newline
         newlen= lenstr(newline)
         llen=newlen
c
      endif

      return
      end
