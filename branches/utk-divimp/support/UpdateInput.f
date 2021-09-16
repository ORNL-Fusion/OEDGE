      program updateinput
      implicit none
c
c     This is a quick hack program to update input files to 
c     quickly create sets for parameter scans.  
c
c     Scan through files looking for strings and replace them 
c     with alternate values - echo the changed lines to the 
c     display. 
c
      character*512 line,val
      character*1 series
      integer in,len_v
      integer change_count,len_l,lenstr,contains_string
      external lenstr,contains_string
      real temp
      logical tchanged,dchanged
c
      series = 'F'
      temp   = 2.0
      
c
      
      tchanged = .false.
      dchanged = .false.
      change_count = 0
c
 10   read(5,'(a512)',end=100,err=100) line
c
c     Search for title and change it 
c
      if (.not.tchanged) then 
          val = 'D-105516-'
          len_v = lenstr(val)
          in = contains_string(line,val)
          if (in.gt.0) then 
 
            tchanged = .true.
            write(line(in+9:in+9),'(a1)') series

            len_l= lenstr(line)
            write(0,'(a,a)') 'CHANGED:',line(1:len_l)
         endif
c
      endif
c
c     Search for description and change it
c
      if (.not.dchanged) then 
          val = 'Te=0.8eV'
          len_v = lenstr(val)
          in = contains_string(line,val)
          if (in.gt.0) then 
 
            dchanged = .true.
            write(line(in+3:in+5),'(f3.1)') temp

            len_l= lenstr(line)
            write(0,'(a,a)') 'CHANGED:',line(1:len_l)
         endif
c
      endif


c
c     Search for base temperature data and change it
c
      if (change_count.lt.30) then 

         val = '0.8          0.8'
         len_v = lenstr(val)
         in = contains_string(line,val)
         if (in.gt.0) then 

            change_count= change_count+1
            write(line(in:in+len_v-1),'(f3.1,10x,f3.1)') temp,temp 
            len_l= lenstr(line)
            write(0,'(a,a)') 'CHANGED:',line(1:len_l)

         endif

      endif

      len_l = lenstr(line)
      write(6,'(a)') line(1:len_l)

      goto 10

 100  Continue

      write(0,'(a)') 'Finished Update'


      return
      end



      integer function contains_string(target,val)
      implicit none 
      character*(*) target,val
      integer len_t, len_v, lenstr,in
      external lenstr
c
      len_t = lenstr(target)
      len_v = lenstr(val)
c
      contains_string = -1
c
      do in = 1, len_t-len_v+1
c
         if (val(1:len_v).eq.target(in:in+len_v-1)) then 
            contains_string = in 
            exit
         endif
c
      enddo
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
