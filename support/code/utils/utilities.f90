module utilities



contains



  subroutine find_free_unit_number(unit)
    implicit none
    integer :: unit
    !
    !     FIND_FREE_UNIT_NUMBER:
    !
    !     This routine scans through unit numbers looking for one that
    !     is not currently in use. This number is returned. This code
    !     is based on the assumption that any unit numbers returned will
    !     be used before this routine is called again asking for another 
    !     number - otherwise it will likely return the previous value.
    !
    integer :: test_unit
    logical :: unit_open

    test_unit = 10
    unit_open = .true.

    ! Check for unit number assignment.  
    Do While (Unit_open.or.test_unit.gt.98)
       test_unit=test_unit + 1
       Inquire (Unit = test_unit, Opened = Unit_open)
    End Do

    if (unit_open) then 
       write(0,'(a)') 'ERROR: No unit numbers available:',test_unit
       stop
    else

       unit = test_unit
    endif

    return
  end subroutine find_free_unit_number


  function upcase(string) result(upper)
    character(len=*), intent(in) :: string
    character(len=len(string)) :: upper
    integer :: j

    do j = 1,len(string)
       if(string(j:j) >= "a" .and. string(j:j) <= "z") then
          upper(j:j) = achar(iachar(string(j:j)) - 32)
       else
          upper(j:j) = string(j:j)
       end if
    end do

  end function upcase


      REAL FUNCTION ATAN2C (ARGZ,ARGR)
      implicit none
      REAL ARGZ,ARGR
!     INCLUDE "PARAMS"
!      include 'params'
    real,parameter :: PI=3.141592654
!
!     THIS ACTS AS AN ERROR-CHECKING FRONT-END TO THE ATAN2
!     IMPLICIT FUNCTION. IT RETURNS APPROPRIATE ANGLE VALUES FOR
!     EITHER OF THE ARGUMENTS EQUAL TO ZERO AND RETURNS A
!     ZERO VALUE IF BOTH ARGUMENTS ARE EQUAL TO ZERO. SINCE
!     THE TANGENT IS UNDEFINED IN THIS CASE.
!
!     D. ELDER  SEPTEMBER 1992
!
      IF (ARGZ.EQ.0.0) THEN
         IF (ARGR.GT.0.0) THEN
            ATAN2C = 0.0
         ELSEIF (ARGR.LT.0.0) THEN
            ATAN2C = PI
         ELSE
            ATAN2C = 0.0
         ENDIF
      ELSEIF (ARGR.EQ.0.0) THEN
         IF (ARGZ.GT.0.0) THEN
            ATAN2C = PI /2.0
         ELSEIF (ARGZ.LT.0.0) THEN
            ATAN2C = - PI /2.0
         ELSE
            ATAN2C = 0.0
         ENDIF
      ELSE
         ATAN2C = ATAN2(ARGZ,ARGR)
      ENDIF
      RETURN
    END FUNCTION ATAN2C


    subroutine intsect2dp(ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect)
      implicit none   
      real*8 ra,za,rb,zb,r1,z1,r2,z2,rint,zint
      integer sect

      !     Calculates intersection of two lines and sets flag if the
      !     intersection is between the end points of both line segments.

      !     Instead of just true and false for the intersection - this 
      !     code has been generalized to 4 return values.

      !     0 - intersection does not occur between end-points or no intersection

      !     1 - intersection between the lines occurs between both sets of 
      !         specified end points 

      !     2 - intersection between the lines occurs between the first
      !         set of end-points but not the second

      !     3 - intersection between the lines occurs between the second
      !         set of end-points but not the first

      real*8 ma,m1,ba,b1
      real*8 rdsta,rdstb
      logical verta,vert1 
      real*8 eps
      parameter (eps=1.0d-8)

      integer warnings
      data warnings/0/

      logical :: debug = .false.


      verta = .false.
      vert1 = .false.

      rint = 0.0
      zint = 0.0

      sect = 0

      !     Define slopes of lines

      if (ra.eq.rb) then 
         verta = .true.
         ma = 0.0  
         ba = ra
      else
         ma = (zb-za)/(rb-ra)
         ba = za - ma * ra
      endif

      if (r1.eq.r2) then 
         vert1 = .true.
         m1 = 0.0  
         b1 = r1
      else
         m1 = (z2-z1)/(r2-r1)
         b1 = z1 - m1 * r1
      endif

      !     Debug:

      if (debug) then 
         write(6,'(a,8g20.12)') 'IS2DP:',ra,za,rb,zb,r1,z1,r2,z2
         write(6,'(a,8g20.12)') '      ',ma,ba,m1,b1
         write(6,'(a,6l6)')     '      ',verta,vert1,ba.eq.b1,ma.eq.m1,&
              &                     abs(ba-b1).lt.eps,abs(ma-m1).lt.eps
      endif


      !     Find intersection 

      if (verta.and.vert1) then 

         !        Line segments may overlap
         !        Error - set sect to false and issue error message

         !        Do nothing for parallel case

         if (abs(ba-b1).lt.eps) then 

            warnings = warnings + 1

            write(6,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//&
                 &                     ' COLINEAR-VERTICAL',warnings
            write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//&
                 &                     ' COLINEAR-VERTICAL',warnings
            write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,&
                 &                                  r2,z2,ma,m1,ba,b1

            !           Check for intersection region and return the mid-point of 
            !           the overlap as the intersection point. 

            call find_parallel_intersection(ra,za,rb,zb,r1,z1,r2,z2,&
                 &                                      ma,ba,rint,zint,0)


         endif

      elseif (verta) then 

         !        Line A is vertical - line 1 is not -> ra = rb = rint

         rint = ra
         zint = m1 * rint + b1         

      elseif (vert1) then  

         !        Line 1 is vertical - line A is not -> r1 = r2 = rint

         rint = r1
         zint = ma * rint + ba        

      else

         !        Neiither line vertical - check for parallel or numerically parallel lines
         !  
         if (abs(ma-m1).lt.eps) then 
            !         if (ma.eq.m1) then 

            !           Check for the same line 

            if (abs(ba-b1).lt.eps) then 
               !            if (ba.eq.b1) then 

               warnings = warnings+1

               if (ma.eq.0.0) then 
                  write(6,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//&
                       ' COLINEAR-HORIZONTAL',warnings
                  write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//&
                       ' COLINEAR-HORIZONTAL',warnings
                  write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,&
                       r2,z2,ma,m1,ba,b1
               else
                  write(6,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//&
                       ' COLINEAR-PARALLEL',warnings
                  write(0,'(a,i6)') 'INTSECT2DP:WARNING:LINE SEGMENTS'//&
                       ' COLINEAR-PARALLEL',warnings
                  write(6,'(a,12(1x,g18.10))') 'DATA:',ra,za,rb,zb,r1,z1,&
                       r2,z2,ma,m1,ba,b1
               endif

               !              Check for intersection region and return the mid-point of 
               !              the overlap as the intersection point. 

               call find_parallel_intersection(ra,za,rb,zb,r1,z1,r2,z2,ma,ba,rint,zint,1)

            endif

            !        Calculate intersection 

         else

            rint = (b1-ba)/(ma-m1)
            zint = ma*rint + ba

         endif

      endif

      !     Now determine if rint,zint lies inside both of the rectangles
      !     defined by the two pairs of end points. 

      if (&
           (dabs(dabs(rint-ra)+dabs(rint-rb)-dabs(ra-rb)).le.eps).and.&
           (dabs(dabs(zint-za)+dabs(zint-zb)-dabs(za-zb)).le.eps).and.&  
           (dabs(dabs(rint-r1)+dabs(rint-r2)-dabs(r1-r2)).le.eps).and.&
           (dabs(dabs(zint-z1)+dabs(zint-z2)-dabs(z1-z2)).le.eps)) then 
         sect = 1
      elseif (&
           (dabs(dabs(rint-ra)+dabs(rint-rb)-dabs(ra-rb)).le.eps).and.&
           (dabs(dabs(zint-za)+dabs(zint-zb)-dabs(za-zb)).le.eps)) then 
         sect = 2
      elseif (&
           (dabs(dabs(rint-r1)+dabs(rint-r2)-dabs(r1-r2)).le.eps).and.&
           (dabs(dabs(zint-z1)+dabs(zint-z2)-dabs(z1-z2)).le.eps)) then 
         sect = 3 
      endif

      !     Check for the case of a numerical interection where the point 
      !     is exceptionally close to the wall. 

      if (sect.eq.3) then 
         rdsta = sqrt((ra-rint)**2+(za-zint)**2)
         rdstb = sqrt((rb-rint)**2+(zb-zint)**2)

         !        If the point is close enough to wall then count it as the intersection

         if (rdsta.lt.eps.or.rdstb.lt.eps) then 
            sect=1 
         endif

         !        Print notification if this code triggers

         if (debug) then
            write(6,'(a,2g20.12,2l6)') '      ',rdsta,rdstb,&
                 &                                 rdsta.lt.eps,rdstb.lt.eps
         endif


      endif



      !      write(6,'(a,2g18.10,i8)') '      ',rint,zint,sect


      !      write (6,'(a,6l4,1p,10(g14.7))') 'DEBUG I2A:',verta,vert1,
      !     >  ((abs(rint-ra)+abs(rint-rb)-abs(ra-rb)).lt.eps),
      !     >  ((abs(zint-za)+abs(zint-zb)-abs(za-zb)).lt.eps),
      !     >  ((abs(rint-r1)+abs(rint-r2)-abs(r1-r2)).lt.eps),
      !     >  ((abs(zint-z1)+abs(zint-z2)-abs(z1-z2)).lt.eps),
      !     >  (abs(rint-ra)+abs(rint-rb)-abs(ra-rb)),
      !     >  (abs(zint-za)+abs(zint-zb)-abs(za-zb)),
      !     >  (abs(rint-r1)+abs(rint-r2)-abs(r1-r2)),
      !     >  (abs(zint-z1)+abs(zint-z2)-abs(z1-z2))

      !      write (6,'(a,1p,10(g14.7))') 'DEBUG I2B:',ra,rint,rb,za,zint,zb
      !      write (6,'(a,1p,10(g14.7))') 'DEBUG I2C:',r1,rint,r2,z1,zint,z2
      !      write (6,'(a,1p,10(g14.7))') 'DEBUG I2C:',ma,ba,m1,b1



      return 
    end subroutine intsect2dp

    subroutine find_parallel_intersection(ra,za,rb,zb,r1,z1,r2,z2,&
         ma,ba,rint,zint,&
         vert)
      implicit none
      real*8 ra,za,rb,zb,r1,z1,r2,z2,ma,ba,rint,zint
      integer vert,sect


      !     The lines are known to be co-linear.
      !     vert = 0 - vertical lines
      !     vert = 1 - not vertical lines

      !     The code returns the center of the overlap region of the 2 line
      !     segments if it exists

      real*8 zstart,zend,rstart,rend

      !     vertical lines

      if (vert.eq.0) then 

         !        Organize by Z-coordinate

         zstart=max(min(za,zb),min(z1,z2)) 
         zend  =min(max(za,zb),max(z1,z2))

         !        No overlap

         if (zend.lt.zstart) then 
            !            sect = 0
            rint = 0.0
            zint = 0.0

            !        Get center of overlap region

         else
            !            sect = 1
            zint = (zstart+zend)/2.0
            rint = ra
         endif
         !   
         !     Base analysis on R - and use line equation for Z. 

      else

         rstart=max(min(ra,rb),min(r1,r2)) 
         rend  =min(max(ra,rb),max(r1,r2))


         !        No overlap

         if (rend.lt.rstart) then 
            !            sect = 0
            rint = 0.0
            zint = 0.0

            !        Get center of overlap region

         else

            !            sect = 1

            rint = (rstart+rend)/2.0
            zint =  rint * ma + ba

         endif

      endif

      !      write(6,'(a,10g18.10,i6)') 'WARNING: INTSECT2:'//
      !     >                    ' INTERSECTION REGION IS CO-LINEAR:',
      !     >                      ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect

      return
    end subroutine find_parallel_intersection


!
!  *********************************************************************
!  *                                                                   *
!  *  IPOS   : FINDS NEAREST HIGHER VALUE IN RS ARRAY TO GIVEN R.      *
!  *           RS ARRAY MUST BE IN ASCENDING ORDER                     *
!  *                                                                   *
!  *  CHRIS FAR RELL    FEBRUARY 1989
!  *                                                                   *
!  *********************************************************************
!
      INTEGER FUNCTION IPOS (R, RS, NRS)
      implicit none
      INTEGER NRS,ILOW,IMID
      REAL    R,RS(NRS)
      integer :: in

!
!     NRS = 0 is an error condition - however, it appears that LIM
!     sometimes does this when calculating time points in cases where
!     the case being run is not time dependent so NTS=0. In any case, 
!     IPOS should return some value in error cases - so IPOS will be 
!     set to 1 initially. A fix has been added to LIM setting NTS to 1. 
!
      if (nrs.eq.0) then 
         ipos = 1 
         WRITE (6,'(a,i6,3(1x,g12.5))') ' IPOS ERROR:'//&
                 ' NUMBER OF ELEMENTS IS ZERO',&
                       nrs,r,rs(1),rs(nrs)
         return
      elseif (RS(1).GT.RS(NRS)) then 
         WRITE (6,'(a,i6,3(1x,g12.5))') ' IPOS ERROR: DESCENDING ORDER',&
                       nrs,r,rs(1),rs(nrs)
      endif
!
      ILOW = 0
      IPOS = NRS
      IF (NRS.EQ.1) RETURN
100   CONTINUE
      IMID = (IPOS + ILOW) / 2
      IF (R.GT.RS(IMID)) THEN
        ILOW = IMID
      ELSE
        IPOS = IMID
      ENDIF
      IF (IPOS-ILOW.GT.1) GOTO 100
!
      RETURN
    END FUNCTION IPOS




end module utilities
