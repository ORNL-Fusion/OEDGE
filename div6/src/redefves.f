      subroutine redefine_vessel(nves,
     >                           rves,zves,nbufle,
     >                           rbufle,zbufle,
     >                           nbufmx,nbufx,rbufx,zbufx,
     >                           node_origin,wallredef)
      use mod_params
      use mod_local_baffles
      implicit none
c     include 'params'
c     include 'local_baffles'
c
      integer nves,nbufle,nbufmx,nbufx(mbufx)
      integer node_origin(mves,2),wallredef
c
c     DIVIMP     - real*4 
c     NIMBUS/PIN - real*8  
c
      real*4 rves(mves),zves(mves),
     >       rbufle(mbufle),zbufle(mbufle),
     >       rbufx(mbufx,mbufle),zbufx(mbufx,mbufle)
c
c
c     REDEFINE_VESSEL: This routine loops through the baffles 
c     finding any points where they intersect the vessel wall
c     or each other. It then follows along all intersectiions and 
c     modifies the vessel wall to match these changes. Any baffles
c     that are not included by this process are left alone though
c     they may be moved to fill any space vacated by baffles that
c     have been incorporated into the vessel wall. DIVIMP is not, 
c     at present equipped to deal with free-standing baffles 
c     within the vessel. If this becomes necessary, DIVIMP will be 
c     appropriately modified. 
c
c     David Elder,   Feb  12, 1999.
c
c     Copy BAFFLE data to local common variables for processing   
c
c
c
c     Other local variables 
c
      logical checkpoint,point
      integer curindex,curbaffle,newindex,newbaffle,baffle_count
      integer in,id,ib,nvesn,ibstart
      integer baffleused(mbufx+2)
      real*8  rvesn(mves+mbufx*mbufle),zvesn(mves+mbufx*mbufle)
c
c     Check for existence of baffles - EXIT if none
c
      wallredef = 0
c
      if (nbufle.le.0.and.nbufmx.le.0) return
c 
c     Set starting point for baffle/wall counting 
c
      ibstart = 0    
      point = .false.
      checkpoint = .true.
c
c     Copy all baffles and wall into one array for easier processing
c
c     Wall
c
      if (nves.gt.0) then 
         ibstart = ibstart + 1
         nbufxl(ibstart) = nves
         do in = 1,nves
            rbufxl(ibstart,in) = rves(in)  
            zbufxl(ibstart,in) = zves(in)  
         end do
      endif
c
c     Primary baffle
c
      if (nbufle.gt.0) then 
         ibstart = ibstart + 1
         nbufxl(ibstart) = nbufle
         do in = 1,nbufle
            rbufxl(ibstart,in) = rbufle(in)  
            zbufxl(ibstart,in) = zbufle(in)  
         end do
      endif 
c
c     Copy rest of the baffles
c
      if (nbufmx.gt.0) then
         do ib = 1,nbufmx
            nbufxl(ib+ibstart) = nbufx(ib)
            do in = 1,nbufx(ib)
               rbufxl(ib+ibstart,in) = rbufx(ib,in)
               zbufxl(ib+ibstart,in) = zbufx(ib,in)
            end do
         end do
      endif 
c
c     Set total number of local baffles 
c 
      nbufmxl = nbufmx + ibstart
c
c
c     Assume for now that a given point will only lie on one
c     other baffle or wall and not intersect with two 
c     structures at the same point.  
c
c     Assume that all of the baffle structures are listed 
c     counter-clockwise. This is true for the ones that are 
c     in the g46529.vers1 grid file. 
c
c     Algorithm - loop around first baffle (main vessel wall) - 
c               - check point - is it on a baffle?
c               - check segment - does a baffle point lie on it
c               - follow baffle - check for wall or other baffle on 
c                 segment
c               - check next node - is it on wall or other baffle?
c               - if on another structure - follow it - otherwise 
c                 continue with current structure.
c
c
c     Start on wall which is indexed as the first baffle - stay on
c     that surface until hitting an intersecting surface - then follow
c     that surface - always incrementing the index for each structure. 
c        
      nvesn = 0
      curindex  = 1 
      curbaffle = 1 
c
c     Record baffles that have been incorporated into the wall
c     So that extras can be left as baffles.
c
      do in = 1,nbufmxl
         baffleused(in) = 0
      end do
c
c     Add the initial point to the wall
c
      nvesn = nvesn +1
      rvesn(nvesn) = rbufxl(curbaffle,curindex) 
      zvesn(nvesn) = zbufxl(curbaffle,curindex)
c         
      node_origin(nvesn,1) = curbaffle
      node_origin(nvesn,2) = curindex
c
c     Use GOTO structure since F77 does not support while constructs and 
c     some users still operate in that environment.
c
 100  continue

      call find_next(curindex,curbaffle,newindex,newbaffle,
     >               checkpoint,point)
c
      if (newbaffle.eq.curbaffle) checkpoint = .true. 
c
      if (rbufxl(newbaffle,newindex).ne.rvesn(nvesn).or.
     >    zbufxl(newbaffle,newindex).ne.zvesn(nvesn)) then 

         nvesn = nvesn +1
         rvesn(nvesn) = rbufxl(newbaffle,newindex) 
         zvesn(nvesn) = zbufxl(newbaffle,newindex)
         
         if (point) then 
c
c           If a point was found to lie on a segment of another 
c           surface then the segment is associated with the 
c           start point of the segment on which the intersection
c           occurs. This requires changing the values for 
c           the originating nodes of the previous point as
c           well as setting the values related to the new point.
c
            node_origin(nvesn-1,1) = newbaffle
            if (newindex.eq.1) then   
               node_origin(nvesn-1,2) = nbufxl(newbaffle)
            else
               node_origin(nvesn-1,2) = newindex -1
            endif 
            node_origin(nvesn,1) = newbaffle
            node_origin(nvesn,2) = newindex
c
         else        
c
c           Segment is associated with it's first listed node
c
            node_origin(nvesn,1) = newbaffle
            node_origin(nvesn,2) = newindex
         endif

         baffleused(newbaffle) = 1

      else
 
         checkpoint=.false.
      
      endif 

c
c      write(6,'(a,7i4,l4,4(1x,g12.5))') 'REDEF:',curbaffle,curindex,
c     >      newbaffle,newindex,nvesn,
c     >      node_origin(nvesn,1),node_origin(nvesn,2),checkpoint,
c     >   rvesn(nvesn),zvesn(nvesn)
c
      curindex = newindex
      curbaffle= newbaffle
c
c
c     Exit when the last point of the initial baffle is reached or
c     if nvesn reaches the maximum.
c
      if ((curindex.eq.nves.and.curbaffle.eq.1).or.
     >      (nvesn.ge.mves)) goto 200

      goto 100

c
c     Finish off processing of the new wall definition 
c

200   continue 

c
c     Check to see if nvesn is small enough  
c
      if (nvesn.ge.mves) then 
c
c        Issue error message and return without making any changes
c
         write(0,*) 'ERROR REDEFINING WALL: NVESN > MVES'//
     >              ' - NO CHANGES MADE',nvesn,mves
         write(6,*) 'ERROR REDEFINING WALL: NVESN > MVES'//
     >              ' - NO CHANGES MADE',nvesn,mves
         write(7,*) 'ERROR REDEFINING WALL: NVESN > MVES'//
     >              ' - NO CHANGES MADE',nvesn,mves
         return
c
      endif 
c
c     Set flad indicating wall redefined ...
c
      wallredef = 1  
c
c     Copy new wall into the original wall definition
c
      nves = nvesn
c
      do in = 1,nvesn
         rves(in) = rvesn(in)
         zves(in) = zvesn(in)
      end do
c
c     Zero counting variable - to keep track of remaining baffles
c
      baffle_count = 0
c
c     Do NOT copy in any remaining unused baffles 
c     Just find how many were not incorporated into the vessesl wall 
c
      do ib = 1,nbufmxl
c
         if (baffleused(ib).eq.0) then 
c
            baffle_count = baffle_count + 1
c
         end if
c 
      end do
c
c     Write out new wall for comparison  
c
      write(6,*) 'Remaining Baffle Count:',baffle_count
      write(6,*) 'Redefined Wall:',nves
c
      do in = 1,nves
         write (6,'(2(1x,g14.7))') rves(in),zves(in) 
      end do
c
c     Exit routine -  
c
      return
      end
c
c
c
      subroutine find_next(curindex,curbaffle,newindex,newbaffle,
     >                     checkpoint,point)
      use mod_params
      use mod_local_baffles
      implicit none
      logical checkpoint,point
      integer curindex,curbaffle,newindex,newbaffle
c     include 'params'
c     include 'local_baffles'
c
c     FIND_NEXT: This routine determines the next index or node
c                to be included in the wall by moving from baffle 
c                to baffle.
c
c                There are two conditions checked for:
c                1) Does the current node lie on a surface that
c                   is part of another baffle? 
c                   If so the newindex and newbaffle are set to the
c                   next greater point on the intersecting baffle.
c
c                2) Does the segment to the next node on the current 
c                   baffle contain a node belonging to another baffle?
c                   If yes - set the next baffle and index to the
c                   intersecting point - need to check for multiple 
c                   nodes on the one segment and choose the one
c                   closest to the current node.
c
c                3) If no intersections are found increment the
c                   current index and do not change baffles. 
c
c
c     Check for point interscetions
c
      point = .false.
c
      if (checkpoint) then
         call check_point(curindex,curbaffle,newindex,newbaffle,
     >                    point)
c
         if (curbaffle.ne.newbaffle) return
c
      endif
c
c     Check for segment intersections
c
      call check_segment(curindex,curbaffle,newindex,newbaffle)
c
      if (curbaffle.ne.newbaffle) then 
         checkpoint = .false.
         return 
      endif  
c
c     No interections found - return next element in current baffle
c
c     Explicitly set new indices      
c     
      if (curindex.eq.nbufxl(curbaffle)) then 
         newindex = 1
      else
         newindex = curindex +1  
      endif
c
      newbaffle= curbaffle
c
      return
      end
c
c
c
      subroutine check_point(curindex,curbaffle,newindex,newbaffle,
     >                       point)
      use mod_params
      use mod_local_baffles
      implicit none
      integer curindex,curbaffle,newindex,newbaffle
      logical point
c
c     CHECK_POINT: This routine checks to see if the node defined by
c                  curindex and curbaffle lies on any other segment
c                  of the remaining baffles. It then returns the 
c                  indices of the "higher numbered" end of the segment
c                  where the intersection occurs. 
c
      integer in,ib,nextin
      logical result
      real*8  angle
c
c     include 'params'
c     include 'local_baffles'
c
c     Loop through looking for intersections
c
c     Set default return value 
c
      point = .false.
c
      newbaffle = curbaffle
      newindex  = curindex 
c
      do ib = 1,nbufmxl
c
         if (ib.ne.curbaffle) then  
c 
            do in = 1,nbufxl(ib)
c
               if (in.eq.nbufxl(ib)) then 
                  nextin = 1
               else
                  nextin = in +1 
               endif
c
c              Check for intersection     
c
               call baffle_sect(rbufxl(curbaffle,curindex),
     >                          zbufxl(curbaffle,curindex),
     >                          rbufxl(ib,in),zbufxl(ib,in),
     >                          rbufxl(ib,nextin),zbufxl(ib,nextin),
     >                          angle,result,0)
c
               if (result) then 
c
c                 If the test point is exactly equal to the second
c                 endpoint - then pass back the next segment and 
c                 not the current one. 
c
                  if (rbufxl(curbaffle,curindex).eq.rbufxl(ib,nextin)
     >                  .and.
     >                zbufxl(curbaffle,curindex).eq.zbufxl(ib,nextin))
     >                   then

                     if (nextin.eq.nbufxl(ib)) then 
                        nextin = 1
                     else
                        nextin = nextin +1 
                     endif

c
                  endif 
c
c                  write(6,'(a,5i4,4(1x,g12.5))') 'checkpoint:',
c     >                       curbaffle,curindex,
c     >                       ib,in,nextin
c
c
                  newindex = nextin
                  newbaffle= ib
                  point = .true.
c
                  return 
c
               endif 
c 
            end do

         end if

      end do 
c
      return
      end      
c
c
c
      subroutine check_segment(curindex,curbaffle,newindex,newbaffle)
      use mod_params
      use mod_local_baffles
      implicit none
      integer curindex,curbaffle,newindex,newbaffle
c
c     CHECK_SEGMENT: This routine checks to see if any nodes lie 
c                    on the segment defined by the starting points of
c                    curindex and curbaffle.
c                    If there are multiple nodes on that segment
c                    then the closest to the current node is 
c                    chosen. 
c
c
      integer in,ib,nextin
      logical result
c
      real*8 dist,newdist,angle
      real*4 dz,dr,atan2c
      external atan2c
c
c     include 'params'
c     include 'local_baffles'
c
c     Loop through looking for intersections
c
c     Set default values
c
      newdist = 1.0d38 
c
      newbaffle = curbaffle
      newindex  = curindex 
c
c     Set next index for other end of segment 
c
      if (curindex.eq.nbufxl(curbaffle)) then 
         nextin = 1
      else
         nextin = curindex+1
      endif 
c
c     Calculate angle (slope) of line segment 
c
      dr = rbufxl(curbaffle,nextin)-rbufxl(curbaffle,curindex)
      dz = zbufxl(curbaffle,nextin)-zbufxl(curbaffle,curindex)
      angle = atan2c(dz,dr) 
c
      do ib = 1,nbufmxl
c
         if (ib.ne.curbaffle) then  
c 
            do in = 1,nbufxl(ib)
c
c              Check for intersection     
c
               call baffle_sect(rbufxl(ib,in),zbufxl(ib,in),
     >                          rbufxl(curbaffle,curindex),
     >                          zbufxl(curbaffle,curindex),
     >                          rbufxl(curbaffle,nextin),
     >                          zbufxl(curbaffle,nextin),
     >                          angle,result,1)
c
               if (result) then 
c
                 dist = 
     >             (rbufxl(ib,in)-rbufxl(curbaffle,curindex))**2
     >            +(zbufxl(ib,in)-zbufxl(curbaffle,curindex))**2
c
c                  write(6,'(a,4i4,4(1x,g12.5))') 'checkseg:',
c     >                       curbaffle,curindex,
c     >                       ib,in,dist,newdist
c

                  if (dist.lt.newdist) then                   

                      newindex = in
                      newbaffle= ib
                      newdist  = dist

                  endif 
c
               endif 
c 
            end do

         end if

      end do 
c
      return
      end      
c
c
c
      subroutine baffle_sect(rp,zp,rstart,zstart,rend,zend,
     >                       angle,result,opt)
      use mod_params
      implicit none
c     include 'params'
      integer opt
      logical result
      real*8 rp,zp,rstart,zstart,rend,zend,angle
c
c     BAFFLE_SECT: The purpose of this routine is to determine if the
c                  point rp,zp lies on the line segment defined by
c                  [rstart,zstart] to [rend,zend]  
c
c     Opt 0 : calculate both angles            
c     Opt 1 : Use angle argument as the segment angle
c
c     Check R-coordinate in-bounds
c
      real*8 seg_angle,test_angle 
      real*4 dz,dr,atan2c,dist,dist_to_point
      external atan2c,dist_to_point
c
      result = .false. 
c
      if ((abs(rend-rp)+abs(rp-rstart))
     >      .le.(abs(rend-rstart)+0.002)) then
c
c        Check Z-coordinate in-bounds
c
         if ((abs(zend-zp)+abs(zp-zstart))
     >               .le.(abs(zend-zstart)+0.002)) then
c
c           Compare angles from start point to test point and
c           segment end point - see if they are the same  
c
c
c           Calculate Segment angle
c
            if (opt.eq.0) then
c
c              Calculate segment angle 
c
               dr = rend-rstart 
               dz = zend-zstart
               seg_angle = atan2c(dz,dr) 
c
            elseif (opt.eq.1) then  
c
               seg_angle = angle
c
            endif
c
c           How close is the point to the line?  
c
            dist = dist_to_point(real(rp),real(zp),
     >                  real(rstart),real(zstart),
     >                  real(seg_angle-PI/2.0))
c
c           If point is a millimeter from the line then assume that 
c           intersection occurs ...
c
            if (dist.lt.0.001) then 
 
               result = .true.

c               write(6,'(a,8(1x,g12.5))') 'Intersect:',
c     >                     rp,zp,rstart,zstart,zstart,zend,dist,
c     >                     seg_angle-PI/2.0   
c


            endif 
c
         endif  
c
      endif
c
      return
c
      end
