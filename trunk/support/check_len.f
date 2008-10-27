      program check_len
      implicit none
      
      integer ios,cnt,flen,llen,eof,lncnt,oklines,comloc,dlines
      character*100 filename
      character*512 line
      integer lenstr,findchars,extstr
      external lenstr,findchars,extstr
c
      logical is_comment
      integer start
c
      open(10,file="linelength_output.txt")

      open(8,file="filelist")

      ios=0
      lncnt = 0
      dlines = 0
      oklines = 0

      do while (ios.eq.0)    

         read(8,'(a100)',iostat=ios) filename

         if (ios.eq.0) then 

            flen=lenstr(filename)

            if (flen.le.1) cycle

            open(9,file=filename(1:flen))

            cnt=0
            eof=0

            do while (eof.eq.0)

               read(9,'(a512)',iostat=eof) line

               cnt=cnt+1

               if (eof.eq.0) then

                  llen = extstr(line,start)
            
                  if (line(start:start).eq.'!') then
                     is_comment=.true.
                  else
                     is_comment=.false.
                  endif
 
                  if (llen.gt.132) then 
                     lncnt=lncnt+1

                     if (llen.gt.260) then 
                        dlines=dlines+1
                     endif

                     comloc = findchars(line(1:llen),'!',-1,-1)
                     if (comloc.gt.0.and.comloc.lt.132.or.
     >                   is_comment) then 
                        oklines=oklines+1
                     else
                        
                        write(10,'(a,'':'',3i6,'':'',a,'':'')') 
     >                    filename(1:max(flen,20)),cnt,llen,comloc,
     >                    line(1:min(llen,131))
                     endif

                  endif
               endif
            end do 
 
            close(9)

         end if

      end do 

      write(10,'(a,i6)') 'Total Long Lines found:',lncnt
      write(10,'(a,i6)') 'Total double length Lines found:',dlines
      write(10,'(a,i6)') 'Estimate of acceptable lines:',oklines
      write(10,'(a,i6)') 'Estimate of actual lines too long:',
     >                    lncnt-oklines

 
      write(0,'(a,i6)') 'Total Long Lines found:',lncnt
      write(0,'(a,i6)') 'Total double length Lines found:',dlines
      write(0,'(a,i6)') 'Estimate of acceptable lines:',oklines
      write(0,'(a,i6)') 'Estimate of actual lines too long:',
     >                    lncnt-oklines



      close(8)
      close(10)

      

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
