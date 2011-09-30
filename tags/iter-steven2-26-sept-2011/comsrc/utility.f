c     -*-Fortran-*-
c
      subroutine gridpos(ik,ir,r,z,newinj,outofgrid)
      implicit none
      integer ik,ir
      real    r,z
      logical newinj,outofgrid
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'grbound'
c
c
c     GRIDPOS: This subroutine determines which ik,ir bin
c              the position R,Z is in. It utilizes the
c              checking function INCELL which uses
c              vector cross-products to determine if
c              the test-point is in the given cell.
c
c              Alternatively - for grids for which the cell
c              vertices are not available - it will use the
c              IKXYS,IRXYS values from the older
c              implementation.
c
c              The implicit assumption is that the cell
c              vertices are ordered clockwise.
c
c              If IK and IR are non-zero ... this routine
c              will first check the IK,IR cell and then
c              all immediately neighbouring cells before
c              it starts scanning the entire grid for
c              the bin containing the R,Z position.
c
c              IF an IK,IR bin is not found an error
c              message is issued and the values IR = IRSEP
c              and IK = 1 are returned. The argument
c              OUTOFGRID is set to TRUE.
c              These are the default return values for a
c              target launch. In the case of a wall launch
c              the IK,IR indices returned will be for the
c              centre of the closest outermost ring segment.
c              If newinj = .false. for a wall injection
c              particle that does not fall on the grid then
c              the values of ir and ik passed to the will be
c              routine will be returned if they are valid.
c
c              Note: Under error conditions - when the
c              particle is outside the grid - this routine
c              will return the nearest IN grid bin centre.
c              This means that rings 2, irwall-1 and irtrap+1
c              become the real rings because the rest are
c              additional non-polygons that are usually
c              used by the fluid codes to impose boundary
c              conditions and do not represent valid grid
c              elements. The 1,nks(ir) elements are Ok because
c              the virtual points have usually already been
c              stripped.
c
c              Code which uses this routine MUST be able
c              to deal with outofgrid condition being set. This
c              routine attempts to make reasonable guesses
c              for out of bounds values but these can not
c              be treated as reasonable under all circumstances.
c
c              In addition - in an attempt to make the code
c              more efficient for out of grid cases - this routine
c              has been modified to use the GA15 routines for
c              determining if a point is inside a closed n-sided
c              polygon. There are two used - one defining the outer
c              edge of the grid and the other defining the edge
c              of the core plasma. If the particle is not a new
c              injection then the value of the lastin variable
c              is used to determine the region search strategy.
c              If it is a new wall injection - outside is checked
c              first - then a grid location is searched for.
c
c
c              David Elder, Dec 8 1993
c
c
      external incell
      logical  incell
      real     rsq,zsq,dsq,best,result
      integer  ix,iy,in,jk,jr,ikorg,irorg
      integer  iktmp,iktmp1,irtmp ,iklim,ikdiff
c
c     Initialization
c
      outofgrid = .false.
      irorg = ir
      ikorg = ik
c
c
c     Make initial assumption that newinj particles are IN
c     the grid somewhere.
c
      if (newinj) lastin = 0
c
c
c      write(6,*) 'gridpos:',ik,ir,r,z,newinj,outofgrid,
c     >           cgridopt,xygrid,outgrid,lastin
c
      if (xygrid.eq.0.and.(.not.outgrid)) then
c
c       Check if the particle is in the current cell
c
c        WRITE(50,*) 'GRIDPOS:',ik,ir,lastin,newinj
        if (ik.ne.0.and.ir.ne.0.and.lastin.eq.0) then
           if (incell(ik,ir,r,z)) return
c
c          If it is not in the current cell - check the
c          adjacent cells.
c
           irtmp = ir
           iktmp = ik
c
c          Check same ring - forward and back
c
c           WRITE(50,*) 'GRIDPOS: FOR + BAK'
           if (iktmp.ne.1) then
              ik = iktmp -1
              if (incell(ik,ir,r,z)) return
           endif
           if (iktmp.ne.nks(ir)) then
              ik = iktmp +1
              if (incell(ik,ir,r,z)) return
           endif
c
c          Check next inner ring - same ik - forward and back
c
c           WRITE(50,*) 'GRIDPOS: SAME IK, FOR + BAK'
           ir = irins(iktmp,irtmp)
           ik = ikins(iktmp,irtmp)
           iktmp1 = ik
           if (incell(ik,ir,r,z)) return
           if (iktmp1.ne.1) then
              ik = iktmp1 -1
              if (incell(ik,ir,r,z)) return
           endif
           if (iktmp1.ne.nks(ir)) then
              ik = iktmp1 +1
              if (incell(ik,ir,r,z)) return
           endif
c
c          Check next outer ring - same ik - forward and back
c
c           WRITE(50,*) 'GRIDPOS: NEXT OUTER RING, FOR + BAK'
           ir = irouts(iktmp,irtmp)
           ik = ikouts(iktmp,irtmp)
           iktmp1 = ik
           if (incell(ik,ir,r,z)) return
           if (iktmp1.ne.1) then
              ik = iktmp1 -1
              if (incell(ik,ir,r,z)) return
           endif
           if (iktmp1.ne.nks(ir)) then
              ik = iktmp1 +1
              if (incell(ik,ir,r,z)) return
           endif
c
c          Reset ik,ir to original values
c
           ik = iktmp
           ir = irtmp
         endif
c
c        Perform a search over the entire grid to find cell.
c        This should be optimized by checking the more
c        likely locations first - e.g. near walls and
c        plates. So it will search plates, irwall and irtrap
c        then search the inner plasma regions.
c        In fact - depending on the launch option employed
c        different search options may be chosen. One of
c        which could include checking the original IK,IR of the last
c        particle checked - first.
c
c
c        New particle launches
c
         if (newinj) then
c           WRITE(50,*) 'GRIDPOS: NEWINJ = .TRUE., ENTIRE GRID'
c
c           CNEUTB = 0 check near plates first
c
            if (cneutb.eq.0.and.cneuta.eq.0) then
               lastin = 0
               do 10 ir = irsep,nrs
                  if (float(nks(ir)/2) .eq.float(nks(ir))/2.0) then
                     iklim = nks(ir)/2 - 1
                  else
                     iklim = nks(ir)/2
                  endif
                  do 20 ikdiff = 0,iklim
                     ik = ikdiff + 1
                     if (incell(ik,ir,r,z)) return
                     ik = nks(ir) - ikdiff
                     if (incell(ik,ir,r,z)) return
 20               continue
 10            continue
c
c              Check core plasma
c
               do 30 ir = 1,irsep-1
c
c                 First cell is repeat of last - but incell handles this - relevant to 
c                 S values
c
c                  do 30 ik = 1,nks(ir) -1
c
                  do 30 ik = 1,nks(ir)
c
                     if (incell(ik,ir,r,z)) return
30             continue
c
c              Particle was NOT found in grid -
c              Take default action. For Plate cases
c              this means issuing an error message
c              and determining IK,IR.
c
c
c              Check if ouside plasma region.
c
c               write(6,*) 'before ga15b:',ionwpts,iwindw(2,1)
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c               write(6,*) 'gp:ga15b1:',result,r,z
c
c              If < 0 particle is outside last ring.
c
               if (result.lt.0.0) then
                  lastin = -1
                  outofgrid = .true.
                  ik = 1
                  ir = irsep
c
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                 ' PARTICLE OFF GRID (1):',ikorg,irorg,r,z
c
c              Check if in core
c
               else
                  CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >              icINDW,MAXPTS,RCW,ZCW,icTDUM,icXDUM,icYDUM,6)
c
                  if (result.gt.0.0) then
c
c                    Particle in core
c
                     outofgrid = .true.
                     lastin =1
                     ik = 1
                     ir = 2
c                     write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                 ' PARTICLE OFF GRID (2):',ikorg,irorg,r,z
                  else
c
c                    Particle is in a small region
c                    between polygon edges and the
c                    ionwall defined above - because
c                    they use different methods - find
c                    the closest ik,ir indices and
c                    return those - hopefully this
c                    will not be required often.
c
                     lastin = 0
                     outofgrid = .false.
c
                     call findwall(ik,ir,r,z)
c
c                     write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                 ' PARTICLE OFF GRID (3):',ikorg,irorg,ik,ir,r,z
c
                  endif
               endif
               return
            elseif (cneutb.eq.2.and.cneuta.eq.0) then
c
c              New injection for wall launch - check perimeter
c              cells first. 
c
               lastin = 0
c
c              Check if ouside plasma region.
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c               write(6,*) 'gp:ga15b2:',ionwpts,result,r,z
c
c              If < 0 particle is outside last ring.
c
               if (result.lt.0.0) then
                  lastin = -1

c                 Particle not found in grid - since it is
c                 a new injection for a wall launch ... find
c                 the nearest wall bin to the injection
c                 position and set outofgrid.
c
                  outofgrid = .true.
c
                  call findwall(ik,ir,r,z)
c
c                  write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >             ' PARTICLE OFF GRID (4):',ikorg,irorg,ik,ir,r,z
c
                  return
               endif
c
c              Check rest of plasma for particle
c
               ir = irwall-1
               do 40 ik = 1,nks(ir)
                  if (incell(ik,ir,r,z)) return
 40            continue

c
c              jdemod - check for existence of PFZ before checking trap rings
c
               if (.not.nopriv) then
                  ir = irtrap+1
                  do 50 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) return
 50               continue
               endif

c
c              Check rest of grid
c
               do 60 ir = irwall -2,irsep,-1
                  do 60 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) return
 60            continue
c
c              jdemod - check for existence of PFZ before checking trap rings
c
               if (.not.nopriv) then 
                  do 70 ir = irtrap+2 ,nrs
                     do 70 ik = 1,nks(ir)
                        if (incell(ik,ir,r,z)) return
 70               continue
               endif
c
c              Check core
c
               do 80 ir = 1,irsep-1
c                  do 80 ik = 1,nks(ir)-1
                  do 80 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) return
 80            continue
c
c              Particle not found in grid - since it is
c              a new injection for a wall launch ... find
c              if it is in the core OR find nearest
c              boundary bin.
c
c
c              Check if actually in core
c
               CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >              icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c
               if (result.gt.0.0) then
c
c                 Particle in core
c
                  outofgrid = .true.
                  lastin =1
                  ik = 1
                  ir = 2
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >              ' PARTICLE OFF GRID (5):',ikorg,irorg,r,z
               else
c
c                 Particle is in a small region
c                 between polygon edges and the
c                 ionwall defined above - because
c                 they use different methods - find
c                 the closest ik,ir indices and
c                 return those - hopefully this
c                 will not be required often.
c
                  lastin = 0
                  outofgrid = .false.
c
                  call findwall(ik,ir,r,z)
c
                  write(6,'(a,4i6,2(1x,g12.5))') 'GRIDPOS:'//
     >                 ' PARTICLE OFF GRID (6):',
     >                  ikorg,irorg,ik,ir,r,z
c

               endif
               return
c
            else
c
c            elseif ((cneutb.eq.1.and.cneuta.eq.0).or.cneuta.eq.1) then
c
c              This was originally the intended strategy for the
c              conditions in the commented out elseif - but the
c              strategy seems to be a reasonable one for the
c              default case.
c
c              Check the last original start position first
c              just in case all the particles are originating
c              in the same bin.
c
               lastin = 0
c
               if ((irstold.ge.1.and.irstold.le.nrs).and.
c slmod begin
c...Problem for Intel compiler (non-starndard FORTRAN):
     >             (ikstold.ge.1.and.
     >              ikstold.le.nks(MAX(1,irstold)))) then
c
c     >             (ikstold.ge.1.and.ikstold.le.nks(irstold))) then
c slmod end
                  ik = ikstold
                  ir = irstold
                  if (incell(ik,ir,r,z)) return
               endif
c
c              Scan the full grid and set ikstold and irstold
c
c              Scan the SOL first - since injection is more likely
c              in this region.
c
               do 100 ir = irsep,nrs
                  do 100 ik = 1,nks(ir)
                     if (incell(ik,ir,r,z)) then
                        ikstold = ik
                        irstold = ir
                        return
                     endif
100            continue
c
c              Scan the core
c
               do 110 ir = 1,irsep-1
c
c                 First and last cell of core are the same except for S coordinate issue -
c                 Incell now checks which part of first or last cell to determine particle
c                 prescence.   
c
c                 do 110 ik = 1,nks(ir)-1
c
                  do 110 ik = 1,nks(ir)
c
                     if (incell(ik,ir,r,z)) then
                        ikstold = ik
                        irstold = ir
                        return
                     endif
110            continue
c
c              Error condition - particle may not be in grid
c
c
c              Check if ouside plasma region.
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c               write(6,*) 'gp:ga15b3:',ionwpts,result,r,z
c
c              If result < 0 particle is outside last ring.
c
               if (result.lt.0.0) then
                  lastin = -1
                  outofgrid = .true.
                  ir = irsep
                  ik = 1
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (7):',ikorg,irorg,r,z
c
c              Otherwise may be in core
c
               else

c
c                 Check if actually in core
c
                  CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >               icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c
                  if (result.gt.0.0) then
c
c                    Particle in core
c
                     outofgrid = .true.
                     lastin =1
                     ik = 1
                     ir = 2
c                     write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (8):',ikorg,irorg,r,z
                  else
c
c                    Particle is in a small region
c                    between polygon edges and the
c                    ionwall defined above - because
c                    they use different methods - find
c                    the closest ik,ir indices and
c                    return those - hopefully this
c                    will not be required often.
c
                     lastin = 0
                     outofgrid = .false.
c
                     call findwall(ik,ir,r,z)
c                     write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >                  ' PARTICLE OFF GRID (9):',
c     >                     ikorg,irorg,ik,ir,r,z
c
                  endif
               endif
c
               return
c
            endif
         elseif (.not.newinj) then
c           WRITE(50,*) 'GRIDPOS: NEWINJ = .FALSE., ENTIRE GRID',lastin
c
c           Not a new injection - but it has fallen through to
c           this point - This means either it was not near the
c           last recorded bin or the last recorded bin was not set.
c           i.e. IK = IR = 0 at the start of the routine.
c
c           Since this routine does not scan for ions stepping along
c           the field lines - only initial injections - and neutrals
c           are not expected to take large arbitrary steps :-) - There
c           is no reason to start at any specific place. So scan the
c           SOL and CORE grid in that order and set any appropriate
c           error conditions.
c
c           Check the lastin region - if the particle is still
c           outside the grid region the one can return quickly.
c
c           One does not expect the particle to move from the
c           core to outside the plasma in one step ... so if
c           the particle has changed regions - check for it in the
c           grid.
c
            if (lastin.eq.-1) then
c
c              Check if particle is still outside.
c              return the values passed to the subroutine.
c              Otherwise scan the grid region.
c
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c               write(6,*) 'gp:ga15b4:',ionwpts,result,r,z
c
c              If < 0 particle is outside last ring.
c
               ik = ikorg
               ir = irorg
c               WRITE(50,*) 'GRIDPOS: POLYGON CHECK',result
c               WRITE(50,*) 'GRIDPOS:              ',ionwpts,maxpts
c               WRITE(50,*) 'GRIDPOS:              ',r,z
               if (result.lt.0.0) then
                  outofgrid = .true.
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (10):',ikorg,irorg,r,z
                  return
               endif
c
            elseif (lastin.eq.1) then
               CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >              icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c
c              If > 0 particle is inside core ring.
c
c               write(6,*) 'gp:ga15b5:',ioncpts,result,r,z
c
               ik = ikorg
               ir = irorg
               if (result.gt.0.0) then
                  outofgrid = .true.
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >             ' PARTICLE OFF GRID (11):',ikorg,irorg,r,z
                  return
               endif
            endif
c
c
c           Assume particle is now in grid.
c
            lastin = 0
c
c           Scan the SOL first - since injection is more likely
c           in this region.
c
c            WRITE(50,*) 'GRIDPOS: SCAN SOL+PFZ'
            do 120 ir = irsep,nrs
               do 120 ik = 1,nks(ir)
                  if (incell(ik,ir,r,z)) return
120         continue
c
c           Scan the core
c
c            WRITE(50,*) 'GRIDPOS: SCAN CORE'
            do 130 ir = 1,irsep-1
c               do 130 ik = 1,nks(ir)-1
               do 130 ik = 1,nks(ir)
                  if (incell(ik,ir,r,z)) return
130         continue
c
c
c           Error condition - particle may not be in grid
c
c
c           Check if ouside plasma region.
c
c            WRITE(50,*) 'GRIDPOS: CHECK IF OUTSIDE PLASMA REGION'
            CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
c            write(6,*) 'gp:ga15b6:',ionwpts,result,r,z
c
c           If result < 0 particle is outside last ring.
c
c            WRITE(50,*) 'GRIDPOS: CHECKING POLYGON AGAIN',result
c            WRITE(50,*) 'GRIDPOS:                       ',r,z
            if (result.lt.0.0) then
               lastin = -1
               outofgrid = .true.


               if (cneutb.eq.2.and.cneuta.eq.0) then
c
c                 If a wall launch that is not in the grid
c                 and not a new injection then return the
c                 ik and ir values that were passed to
c                 the routine - if they are valid.
c
                  if ((irorg.ge.1.and.irorg.le.nrs).and.
     >               (ikorg.ge.1.and.ikorg.le.nks(irorg))) then

                     ir = irorg
                     ik = ikorg
                  else
                     ir = irsep
                     ik = 1
                  endif
               else
c
c                 Default failure values
c
                  ir = irsep
                  ik = 1
               endif

c               write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (12):',
c     >                       ikorg,irorg,r,z


c
c           Otherwise may be in core
c
            else
c
c              Check if actually in core
c
               CALL GA15B(R,Z,RESULT,IONcPTS,1,icWORK,4*MAXPTS,
     >               icINDW,MAXPTS,RcW,ZcW,icTDUM,icXDUM,icYDUM,6)
c               WRITE(50,*) 'GRIDPOS: IN CORE?',result
c
c               write(6,*) 'gp:ga15b8:',ioncpts,result,r,z
c
               if (result.gt.0.0) then
c
c                 Particle in core
c
                  outofgrid = .true.
                  lastin =1
                  ik = 1
                  ir = 2
c
c                  write(6,'(a,2i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (13):',
c     >                ikorg,irorg,r,z
c
               else
c
c                 Particle is in a small region
c                 between polygon edges and the
c                 ionwall defined above - because
c                 they may use different methods - find
c                 the closest ik,ir indices and
c                 return those - hopefully this
c                 will not be required often.
c
c                  WRITE(50,*) 'GRIDPOS: DISPIRATION SET IN...'
                  lastin = 0
                  outofgrid = .false.
c
                  call findwall(ik,ir,r,z)
c
c                  write(6,'(a,4i6,2(1x,g12.5)') 'GRIDPOS:'//
c     >               ' PARTICLE OFF GRID (14):',
c     >                ikorg,irorg,ik,ir,r,z

c
               endif
            endif
c
            return
c
         endif
      elseif (xygrid.eq.1.or.outgrid) then
        IX   = MAX (1, MIN (NXS, INT((R-RMIN)/DR)+1))
        IY   = MAX (1, MIN (NYS, INT((Z-ZMIN)/DZ)+1))
        call gridcoords(ix,iy,ik,ir,in)
c
c        IK   = IKXYS(IX,IY)
c        IR   = IRXYS(IX,IY)
c
c       Set outofgrid to true if particle is outside of the
c       defined grid region - even if one is using the
c       X,Y grid.
c
        if (in.ne.1) outofgrid = .true.
c
      endif

      return
      end
c
c
c

      subroutine findwall(ik,ir,r,z)
      implicit none
      include 'params'
      include 'cgeom'
      integer ik,ir
      real r,z
c
c     Check the targets, the last real wall ring
c     the last real trap ring and the last
c     core ring and return the bin indices of
c     the centre nearest to r,z
c
      real rsq,zsq,dsq,best
      integer jr,jk,id
c
      RSQ = R**2
      ZSQ = Z**2
      best = hi
c
c     Check wall rings
c
      do 200 jr = irwall-1,irtrap+1,
     >              ((irtrap+1)-(irwall-1))
         DO 200 JK = 1, NKS(JR)
            DSQ = RS(JK,JR)**2-2.0*RS(JK,JR)*R+RSQ +
     >            ZS(JK,JR)**2 - 2.0*ZS(JK,JR)*Z + ZSQ
            IF (DSQ.LT.BEST) THEN
               BEST = DSQ
               IK   = JK
               IR   = JR
            ENDIF
 200  CONTINUE
c
c     Check target plates too.
c
      do 210 id = 1,nds
         jr = irds(id)
         if (jr.eq.irwall.or.jr.eq.irtrap) goto 210
         jk = ikds(id)
         DSQ = RS(JK,JR)**2-2.0*RS(JK,JR)*R+RSQ +
     >         ZS(JK,JR)**2 - 2.0*ZS(JK,JR)*Z + ZSQ
         IF (DSQ.LT.BEST) THEN
            BEST = DSQ
            IK   = JK
            IR   = JR
         ENDIF
 210  CONTINUE
c
c     Check core ring
c
      jr = 2
      do 220 jk = 1,nks(jr)
         DSQ = RS(JK,JR)**2-2.0*RS(JK,JR)*R+RSQ +
     >         ZS(JK,JR)**2 - 2.0*ZS(JK,JR)*Z + ZSQ
         IF (DSQ.LT.BEST) THEN
            BEST = DSQ
            IK   = JK
            IR   = JR
         ENDIF
 220  CONTINUE
c
c
      return
      end


c
c
c
      logical function incell(ik,ir,r,z)
      implicit none
      integer ik,ir
      real r,z
      include 'params'
      include 'cgeom'
c
c     INCELL: This function returns a simple YES/NO decision
c             about whether the point R,Z is in the cell designated
c             by IK,IR with a set of vertices defined in an ordered
c             clockwise fashion. It takes the cross product
c             between the vector from the vertex to the test point
c             and the vector from the vertex to the next clockwise
c             vertex of the polygon. The cross-product must be
c             the same sign for all vertices - if the
c             point is outside the polygon it will fail this test
c             for at least one vertex. (i.e. the cross-product will
c             be less than zero.) (Suggested solution courtesy
c             of Ian Youle :-) )
c
c             David Elder, Dec 8, 1993
c
c             Note: the objectives of the solution method were
c             simplicity and reasonable computational cost.
c             This solution avoids the need for square roots
c             or trigonometric calculations.
c
c             Note: in the confined plasma the first and last cells
c                   are identical. However, S=0 and S=SMAX are at the 
c                   center of this cell. This causes some inconsistencies
c                   when calculating particle positions. In order
c                   to address this - this routine will consder a paricle
c                   in the first cell on a core ring when it is in the
c                   second half of the cell and in the last cell of a core
c                   ring when it is in the first half of the cell. 
c
      integer k,v,nextv,i,nv
      real vxr,vxz,vwr,vwz,cp,lastcp
c
      logical inpoly,res
      external inpoly
      real rc(4),zc(4)
c
      lastcp = 0.0
c
      incell = .false.
      k = korpg(ik,ir)
      if (k.eq.0) return
      nv = nvertp(k)
      if (nv.eq.0) return
      do 10 v = 1,nv
         if (v.eq.nv) then
            nextv = 1
         else
            nextv = v+1
         endif
c
c        Want the vector cross-product Rx X Rw
c
c         vxr = r - rvertp(v,k)
c         vxz = z - zvertp(v,k)
c         vwr = rvertp(nextv,k) - rvertp(v,k)
c         vwz = zvertp(nextv,k) - zvertp(v,k)
c
c         cp = vxr*vwz - vxz*vwr
c
c         if (cp.lt.0.0)  return
c
c          if (   (
c     >     ( (r-rvertp(v,k)) *
c     >       (zvertp(nextv,k)-zvertp(v,k)) )
c     >    -( (z-zvertp(v,k)) *
c     >       (rvertp(nextv,k)-rvertp(v,k)) )
c     >           )
c     >         .lt.0.0) return
c
          cp =    (
     >     ( (r-rvertp(v,k)) *
     >       (zvertp(nextv,k)-zvertp(v,k)) )
     >    -( (z-zvertp(v,k)) *
     >       (rvertp(nextv,k)-rvertp(v,k)) )
     >           )
c
c         There is a problem for points that should 
c         lie on the boundary of the cell - i.e. that 
c         are calculated based on the polygon corners and 
c         which are mathematically on the polygon surface. 
c         Numerically, these points can have a cross product
c         which is close to zero but can vary to either side. 
c         In order to consider these points in the cell - the 
c         cross products are set to zero for values less than
c         a specified limit. In this case the limit is set to 1.0e-7 
c
c         This value was determined by examining the range of cross 
c         product values generated when sampling 50,000 points 
c         calculated on a polygon with a scale size of 1.0m. 
c         The maximum error cross product in this case was 6e-8.
c
c         D. Elder, Dec 13, 2006
c
c         Upon consideration - it might be best to not allow these points
c         to be considered inside the cell since if they are detected they
c         can be moved slightly to an appropriate location. 
c
c          if (abs(cp).lt.1.0e-7) cp = 0.0 
c
          if (v.eq.1.and.cp.ne.0.0) lastcp = cp
c
          if ((lastcp * cp).lt.0.0) return
c
          if (cp.ne.0.0) lastcp = cp
c
10    continue
c
c     Particle has been found in cell
c
      incell = .true.
c
c     Check particles in the first or last cell of core rings
c     for more accurate assessement.
c
      if (ir.lt.irsep.and.
     >   (ik.eq.1.or.ik.eq.nks(ir))) then 
c
c        For a particle found to be in the first or last cell of 
c        a core ring - need to decide if it is in the first
c        half or second half and revise incell result accordingly. 
c
c        Only the second half of the first cell or the first half
c        of the last cell should return true
c
c
c        NOTE: Keep in mind that the first and last cells of core
c              rings are supposed to be identical - if this changes
c              then this code needs to be modified. 
c         
c        Check to see if particle is in the first half of the 
c        cell - set vertices for call to inpoly.-  
c         
c        k was set to cell geometry index at the beginning of this
c        routine.
c
         rc(1) = rvertp(1,k)
         rc(2) = rvertp(2,k)
         rc(3) = (rvertp(2,k) + rvertp(3,k)) /2.0
         rc(4) = (rvertp(1,k) + rvertp(4,k)) /2.0
c
         zc(1) = zvertp(1,k)
         zc(2) = zvertp(2,k)
         zc(3) = (zvertp(2,k) + zvertp(3,k)) /2.0
         zc(4) = (zvertp(1,k) + zvertp(4,k)) /2.0
c
c        Check to see of point is in first half of cell
c
         res = inpoly(r,z,4,rc,zc)
c
c        Change value of incell to false for either of the 
c        invalid cases
c        First half and in first cell
c        Last half and in last cell. 
c
c
         if ((ik.eq.1.and.res).or.
     >       (ik.eq.nks(ir).and.(.not.res))) then 
            incell = .false.
         endif

      endif


      return
      end




c
c
c
      logical function incell_debug(ik,ir,r,z)
      implicit none
      integer ik,ir
      real r,z
      include 'params'
      include 'cgeom'
c
c     INCELL: This function returns a simple YES/NO decision
c             about whether the point R,Z is in the cell designated
c             by IK,IR with a set of vertices defined in an ordered
c             clockwise fashion. It takes the cross product
c             between the vector from the vertex to the test point
c             and the vector from the vertex to the next clockwise
c             vertex of the polygon. The cross-product must be
c             the same sign for all vertices - if the
c             point is outside the polygon it will fail this test
c             for at least one vertex. (i.e. the cross-product will
c             be less than zero.) (Suggested solution courtesy
c             of Ian Youle :-) )
c
c             David Elder, Dec 8, 1993
c
c             Note: the objectives of the solution method were
c             simplicity and reasonable computational cost.
c             This solution avoids the need for square roots
c             or trigonometric calculations.
c
c             Note: in the confined plasma the first and last cells
c                   are identical. However, S=0 and S=SMAX are at the 
c                   center of this cell. This causes some inconsistencies
c                   when calculating particle positions. In order
c                   to address this - this routine will consder a paricle
c                   in the first cell on a core ring when it is in the
c                   second half of the cell and in the last cell of a core
c                   ring when it is in the first half of the cell. 
c
      integer k,v,nextv,i,nv
      real vxr,vxz,vwr,vwz,cp,lastcp
c
      logical inpoly,res
      external inpoly
      real rc(4),zc(4)
c
      lastcp = 0.0
c
      incell_debug = .false.
      k = korpg(ik,ir)
      if (k.eq.0) return
      nv = nvertp(k)
      if (nv.eq.0) return
      do 10 v = 1,nv
         if (v.eq.nv) then
            nextv = 1
         else
            nextv = v+1
         endif
c
c        Want the vector cross-product Rx X Rw
c
c         vxr = r - rvertp(v,k)
c         vxz = z - zvertp(v,k)
c         vwr = rvertp(nextv,k) - rvertp(v,k)
c         vwz = zvertp(nextv,k) - zvertp(v,k)
c
c         cp = vxr*vwz - vxz*vwr
c
c         if (cp.lt.0.0)  return
c
c          if (   (
c     >     ( (r-rvertp(v,k)) *
c     >       (zvertp(nextv,k)-zvertp(v,k)) )
c     >    -( (z-zvertp(v,k)) *
c     >       (rvertp(nextv,k)-rvertp(v,k)) )
c     >           )
c     >         .lt.0.0) return
c
          cp =    (
     >     ( (r-rvertp(v,k)) *
     >       (zvertp(nextv,k)-zvertp(v,k)) )
     >    -( (z-zvertp(v,k)) *
     >       (rvertp(nextv,k)-rvertp(v,k)) )
     >           )
c
c
c         There is a problem for points that should 
c         lie on the boundary of the cell - i.e. that 
c         are calculated based on the polygon corners and 
c         which are mathematically on the polygon surface. 
c         Numerically, these points can have a cross product
c         which is close to zero but can vary to either side. 
c         In order to consider these points in the cell - the 
c         cross products are set to zero for values less than
c         a specified limit. In this case the limit is set to 1.0e-7 
c
c         This value was determined by examining the range of cross 
c         product values generated when sampling 50,000 points 
c         calculated on a polygon with a scale size of 1.0m. 
c         The maximum error cross product in this case was 6e-8.
c
c         D. Elder, Dec 13, 2006
c
c         Upon consideration - it might be best to not allow these points
c         to be considered inside the cell since if they are detected they
c         can be moved slightly to an appropriate location. 
c
c          if (abs(cp).lt.1.0e-7) cp = 0.0 
c
          write(6,'(a,5i6,10(1x,g15.7))') 'INCELL_DEBUG:',ik,ir,
     >             korpg(ik,ir),v,k,
     >             cp,lastcp,lastcp*cp,r,z,rvertp(v,k),
     >             zvertp(v,k)
c
c         For cell scale sizes on the order of 10m - the maximum 
c
c          if (abs(cp).lt.1.0e-5)) cp = 0.0
c
          if (v.eq.1.and.cp.ne.0.0) lastcp = cp
c
          if ((lastcp * cp).lt.0.0) return
c
          if (cp.ne.0.0) lastcp = cp
c
10    continue
c
c     Particle has been found in cell
c
      incell_debug = .true.
c
c     Check particles in the first or last cell of core rings
c     for more accurate assessement.
c
      if (ir.lt.irsep.and.
     >   (ik.eq.1.or.ik.eq.nks(ir))) then 
c
c        For a particle found to be in the first or last cell of 
c        a core ring - need to decide if it is in the first
c        half or second half and revise incell result accordingly. 
c
c        Only the second half of the first cell or the first half
c        of the last cell should return true
c
c
c        NOTE: Keep in mind that the first and last cells of core
c              rings are supposed to be identical - if this changes
c              then this code needs to be modified. 
c         
c        Check to see if particle is in the first half of the 
c        cell - set vertices for call to inpoly.-  
c         
c        k was set to cell geometry index at the beginning of this
c        routine.
c
         rc(1) = rvertp(1,k)
         rc(2) = rvertp(2,k)
         rc(3) = (rvertp(2,k) + rvertp(3,k)) /2.0
         rc(4) = (rvertp(1,k) + rvertp(4,k)) /2.0
c
         zc(1) = zvertp(1,k)
         zc(2) = zvertp(2,k)
         zc(3) = (zvertp(2,k) + zvertp(3,k)) /2.0
         zc(4) = (zvertp(1,k) + zvertp(4,k)) /2.0
c
c        Check to see of point is in first half of cell
c
         res = inpoly(r,z,4,rc,zc)
c
c        Change value of incell to false for either of the 
c        invalid cases
c        First half and in first cell
c        Last half and in last cell. 
c
c
         if ((ik.eq.1.and.res).or.
     >       (ik.eq.nks(ir).and.(.not.res))) then 
            incell_debug = .false.
         endif

      endif


      return
      end


c
c
c
      subroutine position_on_target(r,z,cross,id)
      implicit none
      integer id 
      real r,z,cross
      include 'params'
      include 'cgeom'
      include 'comtor'  
c
c     POSITION_ON_TARGET: This routine returns an estimate of 
c                         the actual R,Z location where a 
c                         particle strikes the target - by 
c                         using the CROSS coordinate of the 
c                         particle to move the strike point from
c                         the cell center to a point that 
c                         is proportionally located along 
c                         the target with respect to the 
c                         value of CROSS.
c            
c     Local variables
c     
      integer iw,ik,ir,in
      real dist_frac,direction_factor
c     IPP/08 Krieger - additional variables for iteration of
c     out-of-cell particles
      real rstep,zstep,rorig,zorig
c
      logical test_result,incell
      external incell
c
      integer loop_cnt
      integer :: max_loop_cnt=1000
      real :: step_frac = 1.0/10000.0
c
      iw = wallindex(id)  
      ik = ikds(id)
      ir = irds(id)  
c
c     CROSS is positive when moving "INWARD" and negative when
c     moving "OUTWARD". However, all target and wall elements are
c     ordered clockwise. This means that on the "INNER" or first
c     target the positive CROSS values correspond to moving from 
c     the Rcenter and Zcenter to Rend and Zend of the wall element
c     while on the other target the same positive CROSS value 
c     corresponds to the element between Rcenter,Zcenter and
c     Rstart,Zstart for the target element. 
c
c     The code below calculates the fraction of the way to the next cell
c     that the particle has reached. It then adds a direction_factor to
c     the cross fraction so that positive dist_frac values will always
c     refer to clockwise fractions of the target element and negative
c     dist_frac values to counter-clockwise. 
c
c     The coordinates of Rstart,Zstart Rcenter, Zcenter and Rend,Zend are 
c     then used to determine an actual R,Z location that represents the 
c     approximate location of impact of the particle - this R,Z is stored
c     in the R,Z arguments passed to the routine. If the R,Z of the 
c     particle was passed in then these would be replaced by the target 
c     impact values calculated here.   
c
c
c      write(6,'(a,4i5,8(1x,g12.5))') 'TARGPOS:',ik,ir,iw,id,r,z
c
      if (cross.ge.0.0) then 
c
         if (distin(ik,ir).ne.0.0) then  
            dist_frac = min(cross/distin(ik,ir),1.0)
         else
            dist_frac = 0.0
         endif
c
      elseif (cross.lt.0.0) then 
c
         if (distin(ik,ir).ne.0.0) then  
            dist_frac = max(cross/distout(ik,ir),-1.0)
         else
            dist_frac = 0.0
         endif
c
      endif
c
      if (id.le.ndsin) then 
         direction_factor = 1.0
      else
         direction_factor = -1.0
      endif   
c
      dist_frac = dist_frac * direction_factor
c
c     Calculate R,Z values 
c
      if (dist_frac.ge.0.0) then 
c  
         r = wallpt(iw,1) + dist_frac*(wallpt(iw,22)-wallpt(iw,1))
         z = wallpt(iw,2) + dist_frac*(wallpt(iw,23)-wallpt(iw,2))
c
      else   
c
         r = wallpt(iw,1) - dist_frac*(wallpt(iw,20)-wallpt(iw,1))
         z = wallpt(iw,2) - dist_frac*(wallpt(iw,21)-wallpt(iw,2))
c
      endif 
c
c
c     Test to see if the point is actually inside the cell. 
c
      test_result = incell(ik,ir,r,z)
c
c     If the point is not inside the cell - move it along the line
c     to the cell center by small incremental amounts until it is
c     found to be inside the cell.
c
      loop_cnt =0

c
c     This implementation results in smaller consecutive step sizes but with 
c     a typical step of 1/10000 of the distance to the cell center this should 
c     not be a big effect. 
c
      if (.not.test_result) then 


c        IPP/08 Krieger - the previous iteration scheme was numerically
c        instable due to precision limits (large absolute value of r
c        combined with radially narrow cells lead to diverging error)
c        Solution: replace recursion scheme by explicit iteration
         rorig = r
         zorig = z
         rstep = (rs(ik,ir)-r) * step_frac
         zstep = (zs(ik,ir)-z) * step_frac

         do while (.not.test_result.and.loop_cnt.lt.max_loop_cnt) 

            loop_cnt = loop_cnt+1 
          
            r = rorig + rstep*loop_cnt
            z = zorig + zstep*loop_cnt
         
            test_result = incell(ik,ir,r,z)

         end do
c
c         do while (.not.test_result.and.loop_cnt.lt.max_loop_cnt) 
c
c            loop_cnt = loop_cnt+1 
c          
c            r = r + (rs(ik,ir)-r) * step_frac
c            z = z + (zs(ik,ir)-z) * step_frac
c         
c            test_result = incell(ik,ir,r,z)
c
c         end do
c
c        Point inside cell was not found
c     
         if (.not.test_result) then 
            write(6,'(a,i10,2i6,4g18.8)') 
     >       'ERROR:POSITION_ON_TARGET:'//
     >       ' POINT NOT FOUND IN CELL AFTER MAXIMUM ITERATIONS:',
     >        loop_cnt,ik,ir,r,z,rs(ik,ir),zs(ik,ir)
            write(0,'(a,i10,2i6,4g18.8)') 
     >       'ERROR:POSITION_ON_TARGET:'//
     >       ' POINT NOT FOUND IN CELL AFTER MAXIMUM ITERATIONS:',
     >        loop_cnt,ik,ir,r,z,rs(ik,ir),zs(ik,ir)
         else
            write(6,'(a,i10,2i6,4g18.8)') 
     >       'POSITION_ON_TARGET:'//
     >       ' POINT FOUND IN CELL AFTER N ITERATIONS:',
     >        loop_cnt,ik,ir,r,z,rs(ik,ir),zs(ik,ir)

         endif
      endif

c
c     Exit
c
      return
      end
c
c
c

c
c ======================================================================
c
c subroutine: GetRZ
c
c Given IK,IR,S and CROSS, GETRZ returns an approximate R,Z co-ordinate
c for the impurity... currently the R,Z position is returned assuming
c that CROSS is zero.
c
c ======================================================================
c
      SUBROUTINE GETRZ(IK,IR,S,CROSS,R,Z,opt)
      IMPLICIT none 

C     INCLUDE "PARAMS"
      INCLUDE 'params'
C     INCLUDE "CGEOM"
      INCLUDE 'cgeom'
c     INCLUDE "COMTOR"
      INCLUDE 'comtor'
c
c     jdemod
c
c     Opt is a sub-option to getrz - at the present time it is usually
c     set to rzopt from the input file
c
      integer ik,ir,opt
      real r,z,s,cross,dist
c
c     Local Variables:
c
      REAL ATAN2C,CENLEN

      INTEGER ID
c      REAL    HOLDR,HOLDZ
      REAL    THETA1,THETA2,DS,SR,SZ,DCROSS,CROSSR,CROSSZ,DVR,DVZ
      REAL    RCEN,ZCEN
      REAL    DELTAD,DELTAS
c
c     jdemod
c
      real rside,zside,rtmp,ztmp 
      integer in
c
c     jdemod
c
      LOGICAL message

      DATA message /.FALSE./

      SAVE

c
c     Return cell center values for cells which are listed with zero volume or those 
c     in the boundary cells of the grid - no matter what RZ option has been chosen
c
c Temporary fix? - this was not working for the JET CMOD grid!
c
      IF (KAREAS(IK,IR).EQ.0.0.OR.
     +    IR.EQ.IRWALL.OR.IR.EQ.IRTRAP.OR.IR.EQ.1) THEN
         R = RS(IK,IR)
         Z = ZS(IK,IR)
         RETURN
      ENDIF
c
c     Record values of present position
c 
c      HOLDR=R
c      HOLDZ=Z
c
c     Get cell index
c
      ID = KORPG(IK,IR)
c
c     jdemod - need better estimates in first and last cell on ring - added rzopt=2 for my changes
c
c     Also - code has been written which calls getrz without checking what rz option is in effect - as 
c     a result - I am adding support for all the rzoptions directly to getrz - hopefully this won't break
c     any code since most cases should already have the rz option correctly set in the input file.              
c
c     Option 0 - return cell centers
c
      if (opt.eq.0) then 
         R = RS(IK,IR)
         Z = ZS(IK,IR)
         RETURN
      endif



c
c slnote - Need an option setting here?
c jdemod - adding option settings ... opt 1 is the standard as seen here - opt 2 refines estimates 
c          in the first and last cells of the grids
c
c
      if (opt.eq.1) then 

         IF (S.LT.KSS(IK,IR)) THEN

           IF (IK.EQ.1) THEN
c
c Should be able to use target S values here:
c
             IF (IR.GE.IRSEP) THEN
               R = RS(IK,IR)
               Z = ZS(IK,IR)
               RETURN
             ELSE
c
c slnote - Not sure that this is the best place for the correction
c   for the first knot of the core rings.
c
c Needs some fixing...
c

c             CALL WARN('GetRZ','CORE! Last cell DS calculation',IK,IR,0)
              IF (.NOT.message) THEN
                WRITE(0,*) 'WARNING: LAST CELL DS CALCULATION NEEDS TO'
                WRITE(0,*) '         BE VERIFIED IN CORE (GETRZ)'
                message = .TRUE.
              ENDIF


               DELTAD = CENLEN(NKS(IR),IR,NKS(IR)-1,IR)
               DELTAS = KBACDS(NKS(IR),IR)
             ENDIF

           ELSE

             IF (PDOPT.EQ.1) THEN

               RCEN = 0.5 * (RVERTP(1,ID) + RVERTP(2,ID))
               ZCEN = 0.5 * (ZVERTP(1,ID) + ZVERTP(2,ID))

               DELTAD = SQRT((RCEN - RS(IK,IR))**2 +
     +                       (ZCEN - ZS(IK,IR))**2)
               DELTAS = KSS(IK,IR) - KSB(IK-1,IR)

             ELSE
               DELTAD = CENLEN(IK,IR,IK-1,IR)
               DELTAS = KBACDS(IK,IR)
             ENDIF

           ENDIF

           DS = (S - KSS(IK,IR)) * DELTAD / DELTAS


C
C       S greater than mid-point S in cell 
C
        ELSE

          IF (IK.EQ.NKS(IR)) THEN

            IF (IR.GE.IRSEP) THEN
              R = RS(IK,IR)
              Z = ZS(IK,IR)
              RETURN
            ELSE
c
c Needs some fixin...
c
c            CALL WARN('GetRZ','CORE! Last cell DS calculation',IK,IR,0)
              IF (.NOT.message) THEN
                WRITE(0,*) 'WARNING: LAST CELL DS CALCULATION NEEDS TO'
                WRITE(0,*) '         BE VERIFIED IN CORE (GETRZ)'
                message = .TRUE.
              ENDIF

              DS = (S - KSS(IK,IR)) * CENLEN(1,IR,2,IR) / KFORDS(1,IR)
            ENDIF

          ELSE

            IF (PDOPT.EQ.1) THEN

              RCEN = 0.5 * (RVERTP(3,ID) + RVERTP(4,ID))
              ZCEN = 0.5 * (ZVERTP(3,ID) + ZVERTP(4,ID))

              DS = (S - KSS(IK,IR)) *
     +           SQRT((RCEN - RS(IK,IR))**2 + (ZCEN - ZS(IK,IR))**2) /
     +           (KSB(IK,IR) - KSS(IK,IR))

             ELSE

               DS = (S - KSS(IK,IR)) *
     +              CENLEN(IK,IR,IK+1,IR) / KFORDS(IK,IR)

             ENDIF

           ENDIF

        ENDIF


        DCROSS = 0.0
  
        DVZ = 0.5 * (ZVERTP(3,ID) + ZVERTP(4,ID))-
     +        0.5 * (ZVERTP(1,ID) + ZVERTP(2,ID))

        DVR = 0.5 * (RVERTP(3,ID) + RVERTP(4,ID))-
     +        0.5 * (RVERTP(1,ID) + RVERTP(2,ID))
c
        THETA1 = ATAN2C(ABS(DVZ),ABS(DVR))
        THETA2 = PI / 2.0 - THETA1

        SR = DS * COS(THETA1)
        SZ = DS * SIN(THETA1)

        CROSSR = DCROSS * COS(THETA2)
        CROSSZ = DCROSS * SIN(THETA2)

        IF     (DVR.GE.0.0.AND.DVZ.GE.0.0) THEN
           CROSSZ = -CROSSZ
        ELSEIF (DVR.GE.0.0.AND.DVZ.LT.0.0) THEN
           SZ     = -SZ
           CROSSR = -CROSSR
           CROSSZ = -CROSSZ
        ELSEIF (DVR.LT.0.0.AND.DVZ.GE.0.0) THEN
           SR     = -SR
        ELSEIF (DVR.LT.0.0.AND.DVZ.LT.0.0) THEN
           SR     = -SR
           SZ     = -SZ
           CROSSR = -CROSSR
        ENDIF

        R = RS(IK,IR) + SR + CROSSR
        Z = ZS(IK,IR) + SZ + CROSSZ
c
c     jdemod
c
c     Added additional RZ option based on CELL 
c     Method is exact for a rectangular grid 
c
c     This option uses polygon information - without polygon information it defaults to option 0 
c
      elseif (opt.eq.2) then  
c
c        Find polygon index information
c
         in = korpg(ik,ir)
c
c        No polygon information, zero volume cell, or boundary cell   
c
         if (in.eq.0) then 
               r = rs(ik,ir)
               z = zs(ik,ir) 
               return
         ELSEif (nvertp(in).eq.0.or.kareas(ik,ir).eq.0.0.or.
     >          ir.eq.1.or.ir.eq.irwall.or.ir.eq.irtrap) then 
               r = rs(ik,ir)
               z = zs(ik,ir) 
               return
         endif    
c         
c
c        Need to be careful in identifying the cell polygon corners
c
c        Also - this code garantees that the R,Z value returned actually lies within
c               the specified cell. 
c
c        Brief outline:
c        - find fractional S along the cell centerline
c        - Calculate this centerline R,Z location - RC,ZC 
c        - find matching fractional distances along the parallel to field line sides
c          of the polygon - RS1,ZS1  RS2,ZS2
c        - These points all lie along the same line
c        - use the fractional cross distance to find an R,Z location along this line. 
c
c        Notes: This is not an accurate representation - it is a first order approximation
c               However, it should mimic the actual cells on the grid ... 
c
c        Find fractional S and CROSS distances from the edges of the cell
c        Use the values of distin and distout for cross and the ksb values for
c        the cell for S. 
c
         ds     =   (s - ksb(ik-1,ir)) / (ksb(ik,ir)-ksb(ik-1,ir))
c
c        The parallel to field line sides are: 
c        - OUTWARD23 and INWARD 41
c        - the outward cell face is between polygon vertices 2,3 
c        - the inward cell  face is between polygon vertices 4,1
c
c        - midpoint of the ends of the cells are on sides UP34 and DOWN12 
c
c        Calculate centerline point corresponding to S location
c
         rtmp = ds *  (krb(ik,ir)-krb(ik-1,ir)) + krb(ik-1,ir)
         ztmp = ds *  (kzb(ik,ir)-kzb(ik-1,ir)) + kzb(ik-1,ir)
c
c        Calculate location along line based on cross field location
c
c        INWARD motion direction - towards core and PFZ
c
         if (cross.ge.0.0) then 
c
            dcross = abs(cross/distin(ik,ir))
c
            rside = ds * (rvertp(4,in)-rvertp(1,in)) + rvertp(1,in)
            zside = ds * (zvertp(4,in)-zvertp(1,in)) + zvertp(1,in)
c
            r = dcross * (rside-rtmp) + rtmp
            z = dcross * (zside-ztmp) + ztmp

c
c        OUTWARD motion direction - towards IRWALL
c
         elseif (cross.lt.0.0) then  
c
c           Calculate cross fraction
c
            dcross = abs(cross/distout(ik,ir))
c
c           Calculate the side R,Z locations
c
            rside = ds * (rvertp(3,in)-rvertp(2,in)) + rvertp(2,in)
            zside = ds * (zvertp(3,in)-zvertp(2,in)) + zvertp(2,in)
c
c           Calculate R,Z
c 
            r = dcross * (rside-rtmp) + rtmp
            z = dcross * (zside-ztmp) + ztmp
c
         endif  
c
c         write(6,'(a,3i6,10g12.5)') 'GETRZ2:',ik,ir,in,r,z,s,cross,
c     >            ds,dcross,rside,zside
c 
c
c
c     This option more closely matches the "interpretation" of ion transport in DIVIMP where S and CROSS
c     are essentially orthogonal coordinates at the cell center. However, transport is not conservative
c     when mapping to R,Z and calculated positions outside the cell or even off the grid are possible. 
c
      elseif (opt.eq.3) then  

c
c        Find polygon index information
c
         in = korpg(ik,ir)
c
c        No polygon information, zero volume cell, or boundary cell   
c     
         if (in.eq.0) then 
            r = rs(ik,ir)
            z = zs(ik,ir) 
            return
         elseif (nvertp(in).eq.0.or.kareas(ik,ir).eq.0.0.or.
     >       ir.eq.1.or.ir.eq.irwall.or.ir.eq.irtrap) then 
            r = rs(ik,ir)
            z = zs(ik,ir) 
            return
         endif
c
c        Calculate the distance along the center line to position (from center)
c        - done in steps
c
         ds     =   s - kss(ik,ir)
c
         if (ds.lt.0.0) then 
            ds = ds / (kss(ik,ir)-ksb(ik-1,ir)) 
c
            dist = sqrt( (rs(ik,ir)-krb(ik-1,ir) )**2
     >                 + (zs(ik,ir)-kzb(ik-1,ir) )**2)             
c
            ds = ds * dist 
         elseif (ds.ge.0.0) then 
            ds = ds / (ksb(ik,ir)-kss(ik,ir))
            dist = sqrt((rs(ik,ir)-krb(ik,ir))**2
     >                 +(zs(ik,ir)-kzb(ik,ir))**2)             
            ds = ds * dist 
         endif
c
c        Assign the cross field distance
c
         DCROSS = cross
c
c        Calculate the angle corresponding to the line along the cell center
c
  
         DVZ = 0.5 * (ZVERTP(3,IN) + ZVERTP(4,IN))-
     +         0.5 * (ZVERTP(1,IN) + ZVERTP(2,IN))

         DVR = 0.5 * (RVERTP(3,IN) + RVERTP(4,IN))-
     +         0.5 * (RVERTP(1,IN) + RVERTP(2,IN))
c
c        Remove ABS so that the actual angle is returned - not just a positive one
c
         THETA1 = ATAN2C(DVZ,DVR)
c
c         THETA1 = ATAN2C(ABS(DVZ),ABS(DVR))
c
c        NOTE: changed calculation of THETA2 to be PI/2 clockwise from THETA1 - so
c              it corresponds to the inward direction.  
c
c         THETA2 = PI / 2.0 - THETA1
c
         THETA2 =  THETA1 - PI/2.0

         SR = DS * COS(THETA1)
         SZ = DS * SIN(THETA1)

         CROSSR = DCROSS * COS(THETA2)
         CROSSZ = DCROSS * SIN(THETA2)
c
c        These checks shouldn't be needed if the angles are set up
c        correctly. 
c
c         IF     (DVR.GE.0.0.AND.DVZ.GE.0.0) THEN
c           CROSSZ = -CROSSZ
c         ELSEIF (DVR.GE.0.0.AND.DVZ.LT.0.0) THEN
c           SZ     = -SZ
c           CROSSR = -CROSSR
c           CROSSZ = -CROSSZ
c         ELSEIF (DVR.LT.0.0.AND.DVZ.GE.0.0) THEN
c           SR     = -SR
c         ELSEIF (DVR.LT.0.0.AND.DVZ.LT.0.0) THEN
c           SR     = -SR
c           SZ     = -SZ
c           CROSSR = -CROSSR
c         ENDIF
c
         R = RS(IK,IR) + SR + CROSSR
         Z = ZS(IK,IR) + SZ + CROSSZ

      endif
c
c     jdemod
c

      RETURN
      END
c
c
c END SUBROUTINE: GETRZ ================================================
c
c
c
c ======================================================================
c FUNCTION: CENLEN
c ======================================================================
c      ENDIF

      REAL FUNCTION CENLEN(IK1,IR1,IK2,IR2)

C     INCLUDE "PARAMS"
      INCLUDE 'params'
C     INCLUDE "CGEOM"
      INCLUDE 'cgeom'

      INTEGER IK1,IK2,IR1,IR2
      REAL    RCEN,ZCEN,DELTAL
c
c If the cells are not adjacent and on the same ring then the distance
c returned is simply the R,Z displacement between cell centers:
c
      IF (ABS(IK1-IK2).GT.1.OR.IR1.NE.IR2) THEN
         CENLEN = SQRT((RS(IK1,IR1) - RS(IK2,IR2))**2 +
     +                 (ZS(IK1,IR1) - ZS(IK2,IR2))**2)
      ELSE
c
c slnote - Add an option reference
c
        IF (IK1.LT.IK2) THEN
          RCEN = 0.5 * (RVERTP(3,KORPG(IK1,IR1)) +
     +                  RVERTP(4,KORPG(IK1,IR1)))
          ZCEN = 0.5 * (ZVERTP(3,KORPG(IK1,IR1)) +
     +                  ZVERTP(4,KORPG(IK1,IR1)))
        ELSE
          RCEN = 0.5 * (RVERTP(3,KORPG(IK2,IR2)) +
     +                  RVERTP(4,KORPG(IK2,IR2)))
          ZCEN = 0.5 * (ZVERTP(3,KORPG(IK2,IR2)) +
     +                  ZVERTP(4,KORPG(IK2,IR2)))
        ENDIF

        CENLEN=SQRT((RS(IK1,IR1)-RCEN)**2 + (ZS(IK1,IR1)-ZCEN)**2) +
     +         SQRT((RS(IK2,IR2)-RCEN)**2 + (ZS(IK2,IR2)-ZCEN)**2)
      ENDIF

      RETURN
      END
c
c FUNCTION END: CENLEN ===============================================
c

      subroutine getscross_approx(r,z,s,cross,ik,ir)
      use error_handling
      implicit none
      integer ik,ir
      real r,z,s,cross
      include 'params'
      include 'cgeom'
c
c     GETSCROSS_APPROX: This routine calculates a value
c                       for S, CROSS for a given R,Z position
c                       starting in a given cell. The first order
c                       assumption used in DIVIMP maps the particle
c                       to the cell center resulting in an S value 
c                       equal to KSS(IK,IR) and a CROSS value of 0.0
c  
c                       This routine performs an estimation of the 
c                       position of a given R,Z location along the
c                       two major axes of the cell. For a rectangular
c                       cell this would map correctly to the appropriate 
c                       values of S and CROSS. This is not true for the
c                       non-orthogonal grids presently in use. However, 
c                       it should approximate the location somewhat 
c                       better than the default cell center assumption.
c 
c                       This routine calls the routine 
c                       position_in_poly which takes any 4-sided polygon
c                       and recursively quarters the section containing 
c                       the given R,Z - thus refining its estimate of 
c                       the location of R,Z along the cell axes at each 
c                       iteration. This routine returns a fraction in the
c                       range +/- 1.0 for both S and CROSS - representing
c                       the fraction of the distance IN/OUT or
c                       FORWARD/BACK within the cell towards the cell 
c                       boundaries. This is then used to assign a value for
c                       S and CROSS to be used for the initial ion 
c                       position. 
c        
c                       This is not a perfect solution - simply the 
c                       next level of approximation in trying to 
c                       map between S,CROSS and R,Z which are 
c                       fundamentatlly incompatible on a grid with 
c                       finite and defined cells. 
c
c
c                       NOTE: a special procedure is required for 
c                       particles in the first and last cells of the 
c                       core rings - these cells are identical - each
c                       only maps a half cell - thus
c                       S=0.0 for kss(1,ir) in core and ksb(0,ir) also = 0.0
c                       S=SMAX for kss(nks(ir),ir) in core and ksb(nks(ir),ir) also = SMAX    
c
c
c
c     Local Variables
c
      integer nv,maxnv,iv,maxiter,iter,in,ierr
      parameter(maxnv=4,maxiter=15)
      real rv(maxnv),zv(maxnv)
      real s_frac,cross_frac,base_frac 
c
      logical incell
      external incell
c
      ierr=0

c
c     Set default return values 
c
      s = kss(ik,ir) 
      cross = 0.0

c
c     Dobule check that r,z are actually in this cell - they should be if this routine is called
c     
      if (.not.incell(ik,ir,r,z)) then 
         call errmsg('Problem in Position in poly',
     >                'R,Z does not lie in assigned cell')
         return
      endif

c 
      in = korpg(ik,ir)
c
c     Check for existence of 4 - sided polygonal cell 
c
      if (in.eq.0) return
      if (nvertp(in).ne.4) return
c
      nv = nvertp(in)
c
      do iv = 1,nv
         rv(iv) = rvertp(iv,in)
         zv(iv) = zvertp(iv,in)
      end do
c
c     Initialization for call
c
      iter = 1
      s_frac = 0.0
      cross_frac = 0.0
      base_frac = 1.0
c      
c     Call routine to find the fraction along the cell axes of 
c     given r,z position.  
c
      call position_in_poly(r,z,nv,rv,zv,
     >             s_frac,cross_frac,base_frac,iter,maxiter,ierr)
c
c
c     Calculate approximate S and CROSS values.
c
c     S Coordinate
c
c     Deal with first and last cells of core ring under
c     certain circumstances for the S coordinate
c
      if (ir.lt.irsep.and.
     >   ((ik.eq.1.and.s_frac.lt.0.0).or.
     >     (ik.eq.nks(ir).and.s_frac.gt.0.0))) then 

c
c        First half of first cell - change to last cell
c
         if (ik.eq.1.and.s_frac.lt.0.0) then 
            ik = nks(ir)
            s = kss(ik,ir) + s_frac * (kss(ik,ir) - ksb(ik-1,ir)) 
c
c        Last half of last cell - change to first cell
c
         elseif (ik.eq.nks(ir).and.s_frac.gt.0.0) then
            ik = 1 
            s = kss(ik,ir) + s_frac * (ksb(ik,ir) - kss(ik,ir))        
         endif
c
      else
c
c     Deal with most cases
c

         if (s_frac.ge.0.0) then 
            s = kss(ik,ir) + s_frac * (ksb(ik,ir) - kss(ik,ir)) 
         else 
            s = kss(ik,ir) + s_frac * (kss(ik,ir) - ksb(ik-1,ir)) 
         endif
c
      endif

c
c     CROSS coordinate
c
      if (cross_frac.ge.0.0) then
         cross = cross_frac * distin(ik,ir)
      else
         cross = cross_frac * distout(ik,ir)
      endif



c
      if (ierr.ne.0) then 
         write(6,'(a,2i5,l4,10(1x,g12.5))') 'GETSC:',ik,ir,
     >             incell(ik,ir,r,z),r,z,s,
     >             kss(ik,ir),distin(ik,ir),cross,
     >             distout(ik,ir),s_frac,cross_frac
      endif
c
      return 
      end

      subroutine test_pinp(r,z)
      implicit none
      real r,z
c
c     Local Variables
c
      integer nv,maxnv,iv,maxiter,iter,in,ierr
      parameter(maxnv=4,maxiter=10)
      real rv(maxnv),zv(maxnv)
      real s_frac,cross_frac,base_frac 
c
      ierr = 0
c
c     Test position_in_poly routine
c
      nv = 4 
      do iv = 1,nv
         if (iv.eq.1.or.iv.eq.2) then  
            rv(iv) = 0.0
         else
            rv(iv) = 1.0
         endif 
c
         if (iv.eq.1.or.iv.eq.4) then  
            zv(iv) = 0.0
         else
            zv(iv) = 1.0
         endif 
c
      end do
c
      write(6,'(a,4(1x,g16.8))') 'TEST:R,Z:',r,z 
      write(6,'(a,4(1x,g16.8))') 'RSTART:',(rv(in),in=1,4)
      write(6,'(a,4(1x,g16.8))') 'ZSTART:',(zv(in),in=1,4)
c
      iter = 1
      s_frac = 0.0
      cross_frac = 0.0
      base_frac = 1.0
c         
      call position_in_poly(r,z,nv,rv,zv,
     >             s_frac,cross_frac,base_frac,iter,maxiter,ierr)   
c
      return
      end  
c
c
c
      recursive subroutine position_in_poly(r,z,nvert,rvert,zvert,
     >              s_frac,cross_frac,base_frac,iter,maxiter,ierr)
      implicit none
      integer nvert,iter,maxiter,ierr

      real r,z,rvert(nvert),zvert(nvert)

      real s_frac,cross_frac,base_frac 
c
c
c
c     POSITION_IN_POLY: This routine is invoked recursively - the 
c                       intention is to determine roughly at what fraction
c                       of the way across the two axes of the cell the
c                       input point R,Z lies. It does this by recursively 
c                       dividing the cell into quarters and iterating 
c                       a specified number of times. The point R,Z is always
c                       within the polygon that is passed on to the 
c                       routine when it invokes itself. The s_frac and 
c                       cross_frac values record the cumulative relative
c                       position on the "s" and "cross" axes - (though 
c                       "cross" is not an accurate description of the 
c                       second axis for non-orthogonal cells.) This
c                       routine can be iterated as long as desired to 
c                       obtain any desired level of accuracy.
c
c
c    Local variables 
c
c
c
      integer nv,nvmax,iv,ivf,ivnext,ivlast,in
      parameter(nvmax=4)
      real rv(nvmax),zv(nvmax),rs(nvmax),zs(nvmax),rcp,zcp
      logical found,inpoly
      external inpoly  
c
c     Code is designed to work for 4-sided polygons only.
c  
      if (nvert.ne.4) return
c
c
c     s_frac and cross_frac both start at 0.0
c 
c     base_frac starts initially at 1.0 when passed in
c       
      base_frac = base_frac/2.0
c
c     Split the polygon into 4 pieces and determine which part of the
c     cell the point R,Z lies in - adjust 
c
      nv = nvert
c
      do iv = 1,nv
c
         ivnext = iv+1
         if (iv.eq.nv) ivnext = 1 
c
         rs(iv) = (rvert(iv)+rvert(ivnext))/2.0
         zs(iv) = (zvert(iv)+zvert(ivnext))/2.0
c
      end do
c
      rcp = (rs(1) + rs(3))/2.0
      zcp = (zs(1) + zs(3))/2.0
c
      ivf= 0
      iv = 1
      found = .false. 
c
      do while(iv.le.4.and.(.not.found)) 
         ivlast = iv-1
         if (iv.eq.1) ivlast = 4 
c
c        Determine corners of polygon to check 
c
         rv(mod(iv-1,4)+1)   = rvert(iv) 
         rv(mod(1+iv-1,4)+1) = rs(iv)
         rv(mod(2+iv-1,4)+1) = rcp
         rv(mod(3+iv-1,4)+1) = rs(ivlast)
c
         zv(mod(iv-1,4)+1)   = zvert(iv) 
         zv(mod(1+iv-1,4)+1) = zs(iv)
         zv(mod(2+iv-1,4)+1) = zcp
         zv(mod(3+iv-1,4)+1) = zs(ivlast)
c
         found = inpoly(r,z,nv,rv,zv) 
c
         if (found) ivf = iv
c
         iv = iv+1
c
      end do
c
c
c     Check to make sure location found
c
      if (ivf.ne.0) then 
c
c        The definition of the +/- S and CROSS axes are set here
c        to match the conventions used in DIVIMP relative to the sides of the
c        cell. 
c
c        Side 41 (between polygon corners 1 and 4) is considered INWARD and 
c                corresponds to a positive CROSS displacement.
c        Side 23 (between polygon corners 2 and 3) is considered OUTWARD and 
c                corresponds to a negative CROSS displacement.
c        Side 34 (between polygon corners 3 and 4) is considered UP (larger
c                value of S along the field line) and 
c                corresponds to a positive S displacement (relative to the 
c                cell center)
c        Side 12 (between polygon corners 1 and 2) is considered DOWN (smaller
c                value of S along the field line) and 
c                corresponds to a negative S displacement (relative to the 
c                cell center)
c
         if (ivf.eq.1.or.ivf.eq.2) then 
            s_frac = s_frac - base_frac
         else
            s_frac = s_frac + base_frac
         endif
c
         if (ivf.eq.1.or.ivf.eq.4) then 
            cross_frac = cross_frac + base_frac
         else
            cross_frac = cross_frac - base_frac
         endif
c
         if (iter.eq.maxiter) then 
            return
c
         else
c
c           Increment iteration
c
            iter = iter + 1
c         
            call position_in_poly(r,z,nv,rv,zv,
     >             s_frac,cross_frac,base_frac,iter,maxiter,ierr)   
c
         endif
c
      else
c
c        Error    
c 
         ierr = 1
c
         write(6,'(a,8(1x,g12.5))') 
     >             'ERROR in "position_in_poly": point not'//
     >             ' found in cell: ITER=',iter
         write(6,'(a,8(1x,g12.5))') 'Last poly:',
     >                          (rv(iv),zv(iv),iv=1,4)
c     >                          ((rv(iv),zv(iv)),iv=1,4)  ! gfortran
         write(6,'(a,8(1x,g12.5))') 'Point R,Z,S,C:',r,z,
     >                              s_frac,cross_frac
c
      endif
c
      return 
      end        
c
c
c
      logical function inpoly(r,z,nv,rvert,zvert)
      implicit none
      integer ik,ir
      integer nv,maxvert   
c      parameter (maxvert=8)
      real r,z,rvert(nv),zvert(nv)
c
c     INPOLY: This function returns a simple YES/NO decision
c             about whether the point R,Z is in the cell designated
c             by a set of vertices defined in an ordered fashion.
c             It takes the cross product
c             between the vector from the vertex to the test point
c             and the vector from the vertex to the next
c             vertex of the polygon. The cross-product must be
c             the same sign for all vertices - if the
c             point is outside the polygon it will fail this test
c             for at least one vertex. (i.e. the cross-product will
c             change sign) (Suggested solution courtesy
c             of Ian Youle :-) )
c
c             David Elder, Dec 8, 1993
c
c             Note: the objectives of the solution method were
c             simplicity and reasonable computational cost.
c             This solution avoids the need for square roots
c             or trigonometric calculations.
c
      integer v,nextv
      real cp,lastcp
c
      lastcp = 0.0 
c
      inpoly = .false.

      if (nv.eq.0) return  
c
c     Loop through vertices
c
      do 10 v = 1,nv
c
         if (v.eq.nv) then
            nextv = 1
         else
            nextv = v+1
         endif
c
c        Want the vector cross-product Rx X Rw
c
c         vxr = r - rvert(v)
c         vxz = z - zvert(v)
c         vwr = rvert(nextv) - rvert(v)
c         vwz = zvert(nextv) - zvert(v)
c
c         cp = vxr*vwz - vxz*vwr
c
          cp =    (
     >     ( (r-rvert(v)) *
     >       (zvert(nextv)-zvert(v)) )
     >    -( (z-zvert(v)) *
     >       (rvert(nextv)-rvert(v)) )
     >           )
c
c         There is a problem for points that should 
c         lie on the boundary of the cell - i.e. that 
c         are calculated based on the polygon corners and 
c         which are mathematically on the polygon surface. 
c         Numerically, these points can have a cross product
c         which is close to zero but can vary to either side. 
c         In order to consider these points in the cell - the 
c         cross products are set to zero for values less than
c         a specified limit. In this case the limit is set to 1.0e-7 
c
c         This value was determined by examining the range of cross 
c         product values generated when sampling 50,000 points 
c         calculated on a polygon with a scale size of 1.0m. 
c         The maximum error cross product in this case was 6e-8.
c
c         D. Elder, Dec 13, 2006
c
          if (abs(cp).lt.1.0e-7) cp = 0.0 
c
          if (v.eq.1.and.cp.ne.0.0) lastcp = cp
c
c         Look for change in sign of cp  
c
          if ((lastcp * cp).lt.0.0) return 
c
          if (cp.ne.0.0) lastcp = cp  
c
10    continue
c
      inpoly = .true.
      return
      end
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
C  *********************************************************************
C  *  OPENHTML:  ADDS HTML HEADER INFORMATION                          *
C  *********************************************************************
C
      SUBROUTINE openhtml(casename)
      implicit none
      character*(*) casename
      include 'params'  
c
      integer start,end,extstr,findchars,lenstr
      external extstr,findchars,lenstr
c            
      rewind(htmlunit)
c  
      start = findchars(casename,'/',-1,-1)
      end = lenstr(casename)      
c
      if (start.eq.-1) start = 1
c
      write (htmlunit,*) '<HTML>'
      write (htmlunit,*) '<HEAD>'
      write (htmlunit,*) '  <TITLE>'//casename(start+1:end)//'</TITLE>'
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
      write(htmlunit,*)  '<PRE>'
c 
c
      RETURN
      END

C
C  *********************************************************************
C  *  CLOSEHTML:  ADDS HTML END OF DOCUMENT TAGS                       *
C  *********************************************************************
C
      SUBROUTINE closehtml
      implicit none
      include 'params' 
c
      write(htmlunit,*) '</PRE>'
c 
      write(htmlunit,*) '</BODY>'
      write(htmlunit,*) '</html>'
c
      RETURN
      END
C
C  *********************************************************************
C  *  PRCHTML:  PRINTS A CHARACTER STRING WITH SOME HTML ATTRIBUTES    *
C  *********************************************************************
C
      SUBROUTINE PRCHTML(STRING,REF,HREF,ATTRIB)
      implicit none
      include 'params' 
      CHARACTER STRING*(*),ref*(*),href*(*),attrib*(*)
c
c     This routine adds some simple HTML attributes to the string as 
c     well as any reference or link. 
c
      character*512 line,outstring
      integer stringlen,findchars,lenstr
      external lenstr,findchars
c
      stringlen = lenstr(string)
c
      line = string(1:stringlen)
c
      stringlen = lenstr(line)  
c
      if (findchars(attrib,'B',1,1).ne.-1) then
c         
         line = '<B>'//line(1:stringlen)//'</B>' 
         stringlen = lenstr(line)
c
      endif               
c
      if (findchars(attrib,'U',1,1).ne.-1) then
c         
         line = '<U>'//line(1:stringlen)//'</U>' 
         stringlen = lenstr(line)
c
      endif               
c
      if (findchars(attrib,'I',1,1).ne.-1) then
c         
         line = '<I>'//line(1:stringlen)//'</I>' 
         stringlen = lenstr(line)
c
      endif               
c
c     Add references and anchors
c
      if (ref.ne.'0'.and.href.ne.'0') then 
c
         line = '<A HREF ="#'//href//'" NAME = "'//ref//
     >       '">'//line(1:stringlen)//'</A>'
         stringlen = lenstr(line)
c
         
c
      elseif (ref.ne.'0') then
c
         line = '<A NAME = "'//ref//
     >       '">'//line(1:stringlen)//'</A>'
         stringlen = lenstr(line)
c
      elseif (href.ne.'0') then
c
         line = '<A HREF ="#'//href//
     >       '">'//line(1:stringlen)//'</A>'
         stringlen = lenstr(line)
c
      endif
c
c     Add paragraph formatting so each entry is on a different line
c     
      stringlen = lenstr(line)
      line  = '<P>'//line(1:stringlen)//'</P>'
      stringlen = lenstr(line)
c
      WRITE (tmpunit,'(a)') line(1:stringlen)
c
      stringlen = lenstr(string)
      outstring = '(*) '//string(1:stringlen)  
      stringlen = lenstr(outstring)
c
      WRITE (datunit,'(a)') outstring(1:stringlen)
c
      RETURN
      END
C
      SUBROUTINE TEST
      implicit none
c
C     IPP/01 Krieger - SUN compiler chokes on data initialization
C     in a variable declaration. Added appropiate data statements
c      REAL AS(8) /1.,2.,3.,4.,5.,6.,7.,8./
c      REAL BS(8) /8.,7.,6.,5.,4.,3.,2.,1./
c      REAL DS(7) /-0.1,2.0,2.5,6.5,7.5,8.0,10.0/
      REAL AS(8)
      REAL BS(8)
      REAL DS(7)
c
      integer ipos,jpos,n,i,ip,jp
      EXTERNAL IPOS,JPOS
c
      data AS /1.,2.,3.,4.,5.,6.,7.,8./
      data BS /8.,7.,6.,5.,4.,3.,2.,1./
      data DS /-0.1,2.0,2.5,6.5,7.5,8.0,10.0/
c
      WRITE (6,'(/1X,''TESTING IPOS AND JPOS FUNCTIONS ...'')')
      DO 200 N = 7, 8
        DO 100 I = 1, 7
          IP = IPOS(DS(I),AS,N)
          WRITE (6,9001) 'IPOS',DS(I),'AS',N,IP,AS(IP)
  100   CONTINUE
        DO 110 I = 1, 7
          JP = JPOS(DS(I),BS(8-N+1),N)
          WRITE (6,9001) 'JPOS',DS(I),'BS',N,JP,BS(8-N+JP)
  110   CONTINUE
  200 CONTINUE
 9001 FORMAT(1X,'TEST ... ',A4,' (',F4.1,',',A2,',',I2,') =',I2,
     >  '  (VALUE',F4.1,')')
      RETURN
      END
C
C
C
      SUBROUTINE RDVMF( NAME , IERR )
      implicit none
C
C-----------------------------------------------------------------------
C
C PURPOSE : TO READ IN VMF DATA BLOCKS.
C
C INPUT   : C*      NAME     = NAME FOR ERROR MESSAGES.
C
C COMMON  : I*4     CNVMF    = NUMBER OF VMF BLOCKS.
C           I*4     CIRNG0() = START VMF AT THIS RING.
C           I*4     CIRNG1() = STOP  VMF AT THIS RING.
C           I*4     CJ0()    = FIRST CJ0() POINTS ON A RING.
C           I*4     CJ1()    = LAST  CJ1() POINTS ON A RING.
C           R*4     CVMF0()  = VMF VALUE FOR FIRST CJ0() POINTS.
C           R*4     CVMF1()  = VMF VALUE FOR POINTS BETWEEN REGIONS.
C           R*4     CVMF2()  = VMF VALUE FOR LAST  CJ1() POINTS.
C
C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4866
C           JET/TESSELLA SUPPORT SERVICES PLC
C
C DATE    : 26/10/90
C
C-----------------------------------------------------------------------
C
C     INCLUDE  "PARAMS"
      include 'params'
C     INCLUDE  "DYNAM5"
      include 'dynam5'
C     INCLUDE  "READER"
      include 'reader'
C
      CHARACTER NAME*(*)
      CHARACTER MESAGE*72 , COMENT*72 , HEAD*22
C
      INTEGER   IERR , I
C
C-----------------------------------------------------------------------
C
      CALL RDC( COMENT , NAME , IERR )
C
      CALL RDI( CNVMF , .TRUE. , 0 , .TRUE. , MAXVMF , NAME , IERR )
C
      IF( CNVMF.EQ.0 ) THEN
          READ(5,'(A72)',ERR=9999,END=9999) BUFFER
          READ(5,'(A72)',ERR=9999,END=9999) BUFFER
          READ(5,'(A72)',ERR=9999,END=9999) BUFFER
          RETURN
      END IF
C
C-----------------------------------------------------------------------
C
      MESAGE = 'END OF FILE ON UNIT 5'
C
      DO 100 I = 1 , CNVMF
C
         READ(5,'(A72)',ERR=9999,END=9999) BUFFER
         WRITE(9,'(1X,A72,1X,A6)') BUFFER , 'RDVMF'
         READ(BUFFER,*,ERR=9999,END=9999) HEAD , CIRNG0(I) , CIRNG1(I)
C
         READ(5,'(A72)',ERR=9999,END=9999) BUFFER
         WRITE(9,'(1X,A72,1X,A6)') BUFFER , 'RDVMF'
         READ(BUFFER,*,ERR=9999,END=9999) HEAD , CJ0(I)    , CJ1(I)
         IF( CJ0(I).LT.0 ) CJ0(I) = 0
         IF( CJ1(I).LT.0 ) CJ1(I) = 0
C
         READ(5,'(A72)',ERR=9999,END=9999) BUFFER
         WRITE(9,'(1X,A72,1X,A6)') BUFFER , 'RDVMF'
         READ(BUFFER,*,ERR=9999,END=9999) HEAD , CVMF0(I) , CVMF1(I)
     >                                         , CVMF2(I)
C
         IF( CVMF0(I).LE.0.0E+00 ) CVMF0(I) = 1.0E+00
         IF( CVMF1(I).LE.0.0E+00 ) CVMF0(I) = 1.0E+00
         IF( CVMF2(I).LE.0.0E+00 ) CVMF0(I) = 1.0E+00
C
  100 CONTINUE
C
      RETURN
C
C-----------------------------------------------------------------------
C
 9999 IERR=1
      WRITE(7,'(1X,2A,3(/1X,A))')
     > 'RDVMF: ERROR READING ',NAME,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
C
      END
C
C
C
      SUBROUTINE PRVMF
      implicit none
C
C-----------------------------------------------------------------------
C
C PURPOSE : TO PRINT CVMF(POINT,RING)
C
C COMMON  : I*4     NRS      = NUMBER OF RINGS
C           I*4     NKS(I)   = NUMBER OF PINTS IN RING # I
C           R*4     CVMF(,)  = VMF FOR A POINT ON A RING
C
C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4866
C           JET/TESSELLA SUPPORT SERVICES PLC
C
C DATE    : 26/10/90
C
C-----------------------------------------------------------------------
C
C     INCLUDE   "PARAMS"
      include 'params'
C     INCLUDE   "CGEOM"
      include 'cgeom'
C
      CHARACTER BUFFER*40
C
      INTEGER   IR , IK , ILIMIT
C
C-----------------------------------------------------------------------
C
      DATA ILIMIT/8/
C
C-----------------------------------------------------------------------
C
      CALL PRB
      CALL PRC('VMF FOR EACH POINT ON A RING :-')
C
      DO 100 IR    = 1 , NRS
C
         IF( NKS(IR).LE.ILIMIT ) THEN
             WRITE(datunit,1000) IR , ( CVMF(IK,IR) , IK=1,NKS(IR) )
         ELSE
             WRITE(datunit,1000) IR , ( CVMF(IK,IR) , IK=1,ILIMIT  )
             WRITE(datunit,1010)    ( CVMF(IK,IR) , IK=ILIMIT+1,NKS(IR))
         END IF
  100 CONTINUE
C
C-----------------------------------------------------------------------
C
 1000 FORMAT( 1X , I4 , ' : ' , 8F7.3 )
 1010 FORMAT( 8X              , 8F7.3 )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
CNEXT
      SUBROUTINE FORMCII(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(6)

      POS(1) = INDEX(BUFFER,'''')
      POS(2) = POS(1)+INDEX(BUFFER(POS(1)+1:),'''')
      POS(3) = POS(2)+1
10    IF (BUFFER(POS(3):POS(3)) .EQ. ' ') THEN
         POS(3) = POS(3) + 1
         GOTO 10
      ENDIF
      POS(4) = POS(3)+1
20    IF (BUFFER(POS(4):POS(4)) .NE. ' ') THEN
         POS(4) = POS(4) + 1
         GOTO 20
      ENDIF
      POS(5) = POS(4)+1
30    IF (BUFFER(POS(5):POS(5)) .EQ. ' ') THEN
         POS(5) = POS(5) + 1
         GOTO 30
      ENDIF
      POS(6) = POS(5)+1
40    IF (BUFFER(POS(6):POS(6)) .NE. ' ') THEN
         POS(6) = POS(6) + 1
         GOTO 40
      ENDIF
      FORM = '(  X,A  ,1X,  X,I  ,  X,I  )'
      WRITE(FORM(2:3),'(I2)') POS(1)
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)-1
      WRITE(FORM(13:14),'(I2)') POS(3)-POS(2)-1
      WRITE(FORM(18:19),'(I2)') POS(4)-POS(3)
      WRITE(FORM(21:22),'(I2)') POS(5)-POS(4)
      WRITE(FORM(26:27),'(I2)') POS(6)-POS(5)
      END
CNEXT
      SUBROUTINE FORMCI(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(4)

      POS(1) = INDEX(BUFFER,'''')
      POS(2) = POS(1)+INDEX(BUFFER(POS(1)+1:),'''')
      POS(3) = POS(2)+1
10    IF (BUFFER(POS(3):POS(3)) .EQ. ' ') THEN
         POS(3) = POS(3) + 1
         GOTO 10
      ENDIF
      POS(4) = POS(3)+1
20    IF (BUFFER(POS(4):POS(4)) .NE. ' ') THEN
         POS(4) = POS(4) + 1
         GOTO 20
      ENDIF
      FORM = '(  X,A  ,1X,  X,I  )'
      WRITE(FORM(2:3),'(I2)') POS(1)
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)-1
      WRITE(FORM(13:14),'(I2)') POS(3)-POS(2)-1
      WRITE(FORM(18:19),'(I2)') POS(4)-POS(3)
      END
CNEXT
      SUBROUTINE FORMCR(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(4)

      POS(1) = INDEX(BUFFER,'''')
      POS(2) = POS(1)+INDEX(BUFFER(POS(1)+1:),'''')
      POS(3) = POS(2)+1
10    IF (BUFFER(POS(3):POS(3)) .EQ. ' ') THEN
         POS(3) = POS(3) + 1
         GOTO 10
      ENDIF
      POS(4) = POS(3)+1
20    IF (BUFFER(POS(4):POS(4)) .NE. ' ') THEN
         POS(4) = POS(4) + 1
         GOTO 20
      ENDIF
      FORM = '(  X,A  ,1X,  X,F  .0)'
      WRITE(FORM(2:3),'(I2)') POS(1)
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)-1
      WRITE(FORM(13:14),'(I2)') POS(3)-POS(2)-1
      WRITE(FORM(18:19),'(I2)') POS(4)-POS(3)
      END
CNEXT
      SUBROUTINE FORMCRRR(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(8)

      POS(1) = INDEX(BUFFER,'''')
      POS(2) = POS(1)+INDEX(BUFFER(POS(1)+1:),'''')
      POS(3) = POS(2)+1
10    IF (BUFFER(POS(3):POS(3)) .EQ. ' ') THEN
         POS(3) = POS(3) + 1
         GOTO 10
      ENDIF
      POS(4) = POS(3)+1
20    IF (BUFFER(POS(4):POS(4)) .NE. ' ') THEN
         POS(4) = POS(4) + 1
         GOTO 20
      ENDIF
      POS(5) = POS(4)+1
30    IF (BUFFER(POS(5):POS(5)) .EQ. ' ') THEN
         POS(5) = POS(5) + 1
         GOTO 30
      ENDIF
      POS(6) = POS(5)+1
40    IF (BUFFER(POS(6):POS(6)) .NE. ' ') THEN
         POS(6) = POS(6) + 1
         GOTO 40
      ENDIF
      POS(7) = POS(6)+1
50    IF (BUFFER(POS(7):POS(7)) .EQ. ' ') THEN
         POS(7) = POS(7) + 1
         GOTO 50
      ENDIF
      POS(8) = POS(7)+1
60    IF (BUFFER(POS(8):POS(8)) .NE. ' ') THEN
         POS(8) = POS(8) + 1
         GOTO 60
      ENDIF
      FORM = '(  X,A  ,1X,  X,F  .0,  X,F  .0,  X,F  .0)'
      WRITE(FORM(2:3),'(I2)') POS(1)
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)-1
      WRITE(FORM(13:14),'(I2)') POS(3)-POS(2)-1
      WRITE(FORM(18:19),'(I2)') POS(4)-POS(3)
      WRITE(FORM(23:24),'(I2)') POS(5)-POS(4)
      WRITE(FORM(28:29),'(I2)') POS(6)-POS(5)
      WRITE(FORM(33:34),'(I2)') POS(7)-POS(6)
      WRITE(FORM(38:39),'(I2)') POS(8)-POS(7)
      END
CNEXT
      SUBROUTINE FORMCC(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(4)

      POS(1)=INDEX(BUFFER,'''')
      POS(2)=POS(1)+INDEX(BUFFER(POS(1)+1:),'''')
      POS(3)=POS(2)+INDEX(BUFFER(POS(2)+1:),'''')
      POS(4)=POS(3)+INDEX(BUFFER(POS(3)+1:),'''')
      FORM = '(  X,A  ,  X,A  )'
      WRITE(FORM(2:3),'(I2)') POS(1)
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)-1
      WRITE(FORM(10:11),'(I2)') POS(3)-POS(2)+1
      WRITE(FORM(15:16),'(I2)') POS(4)-POS(3)-1
      END
CNEXT
      SUBROUTINE FORMI(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(2)

      POS(1) = 1
10    IF (BUFFER(POS(1):POS(1)) .EQ. ' ') THEN
         POS(1) = POS(1) + 1
         GOTO 10
      ENDIF
      POS(2) = POS(1)+1
20    IF (BUFFER(POS(2):POS(2)) .NE. ' ') THEN
         POS(2) = POS(2) + 1
         GOTO 20
      ENDIF
      FORM = '(  X,I  )'
      WRITE(FORM(2:3),'(I2)') POS(1)-1
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)
      END
CNEXT
      SUBROUTINE FORMMR(BUFFER,FORM,ANZ)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(22),ANZ,I

      POS(1) = 1
      FORM(1:1)='('
      DO 5 I=1,ANZ+1
         IF (I .NE. 1) POS(I*2-1) = POS(I*2-2)+1
10       IF (BUFFER(POS(I*2-1):POS(I*2-1)) .EQ. ' ') THEN
            POS(I*2-1) = POS(I*2-1) + 1
            GOTO 10
         ENDIF
         POS(2*I) = POS(I*2-1)+1
20       IF (BUFFER(POS(2*I):POS(2*I)) .NE. ' ') THEN
            POS(2*I) = POS(2*I) + 1
            GOTO 20
         ENDIF
         IF (I .NE. ANZ+1) THEN
            WRITE(FORM(10*(I-1)+2:10*I+1),'(A10)') '  X,F  .0,'
         ELSE
            WRITE(FORM(10*(I-1)+2:10*I+1),'(A10)') '  X,F  .0)'
         ENDIF
         IF (I .NE. 1) THEN
           WRITE(FORM(10*(I-1)+2:10*(I-1)+3),'(I2)')
     >           POS(2*I-1)-POS(2*I-2)
         ELSE
           WRITE(FORM(2:3),'(I2)') POS(1)-1
         ENDIF
         WRITE(FORM(10*(I-1)+7:10*(I-1)+8),'(I2)') POS(2*I)-POS(2*I-1)
5     CONTINUE
      END
CNEXT
      SUBROUTINE FORM2I4R(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(12)

      POS(1) = 1
10    IF (BUFFER(POS(1):POS(1)) .EQ. ' ') THEN
         POS(1) = POS(1) + 1
         GOTO 10
      ENDIF
      POS(2) = POS(1)+1
20    IF (BUFFER(POS(2):POS(2)) .NE. ' ') THEN
         POS(2) = POS(2) + 1
         GOTO 20
      ENDIF
      POS(3) = POS(2)+1
30    IF (BUFFER(POS(3):POS(3)) .EQ. ' ') THEN
         POS(3) = POS(3) + 1
         GOTO 30
      ENDIF
      POS(4) = POS(3)+1
40    IF (BUFFER(POS(4):POS(4)) .NE. ' ') THEN
         POS(4) = POS(4) + 1
         GOTO 40
      ENDIF
      POS(5) = POS(4)+1
50    IF (BUFFER(POS(5):POS(5)) .EQ. ' ') THEN
         POS(5) = POS(5) + 1
         GOTO 50
      ENDIF
      POS(6) = POS(5)+1
60    IF (BUFFER(POS(6):POS(6)) .NE. ' ') THEN
         POS(6) = POS(6) + 1
         GOTO 60
      ENDIF
      POS(7) = POS(6)+1
70    IF (BUFFER(POS(7):POS(7)) .EQ. ' ') THEN
         POS(7) = POS(7) + 1
         GOTO 70
      ENDIF
      POS(8) = POS(7)+1
80    IF (BUFFER(POS(8):POS(8)) .NE. ' ') THEN
         POS(8) = POS(8) + 1
         GOTO 80
      ENDIF
      POS(9) = POS(8)+1
90    IF (BUFFER(POS(9):POS(9)) .EQ. ' ') THEN
         POS(9) = POS(9) + 1
         GOTO 90
      ENDIF
      POS(10) = POS(9)+1
100   IF (BUFFER(POS(10):POS(10)) .NE. ' ') THEN
         POS(10) = POS(10) + 1
         GOTO 100
      ENDIF
      POS(11) = POS(10)+1
110   IF (BUFFER(POS(11):POS(11)) .EQ. ' ') THEN
         POS(11) = POS(11) + 1
         GOTO 110
      ENDIF
      POS(12) = POS(11)+1
120   IF (BUFFER(POS(12):POS(12)) .NE. ' ') THEN
         POS(12) = POS(12) + 1
         GOTO 120
      ENDIF
      FORM ='(  X,I  ,  X,I  ,  X,F  .0,  X,F  .0,  X,F  .0,  X,F  .0)'
      WRITE(FORM(2:3),'(I2)') POS(1)-1
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)
      WRITE(FORM(10:11),'(I2)') POS(3)-POS(2)
      WRITE(FORM(15:16),'(I2)') POS(4)-POS(3)
      WRITE(FORM(18:19),'(I2)') POS(5)-POS(4)
      WRITE(FORM(23:24),'(I2)') POS(6)-POS(5)
      WRITE(FORM(28:29),'(I2)') POS(7)-POS(6)
      WRITE(FORM(33:34),'(I2)') POS(8)-POS(7)
      WRITE(FORM(38:39),'(I2)') POS(9)-POS(8)
      WRITE(FORM(43:44),'(I2)') POS(10)-POS(9)
      WRITE(FORM(48:49),'(I2)') POS(11)-POS(10)
      WRITE(FORM(53:54),'(I2)') POS(12)-POS(11)
      END
CNEXT
      SUBROUTINE FORE2I4R(BUFFER,FORM)
      implicit none
      CHARACTER BUFFER*(*),FORM*(*)
      INTEGER POS(12)

      POS(1) = 1
10    IF (BUFFER(POS(1):POS(1)) .EQ. ' ') THEN
         POS(1) = POS(1) + 1
         GOTO 10
      ENDIF
      POS(2) = POS(1)+1
20    IF (BUFFER(POS(2):POS(2)) .NE. ' ') THEN
         POS(2) = POS(2) + 1
         GOTO 20
      ENDIF
      POS(3) = POS(2)+1
30    IF (BUFFER(POS(3):POS(3)) .EQ. ' ') THEN
         POS(3) = POS(3) + 1
         GOTO 30
      ENDIF
      POS(4) = POS(3)+1
40    IF (BUFFER(POS(4):POS(4)) .NE. ' ') THEN
         POS(4) = POS(4) + 1
         GOTO 40
      ENDIF
      POS(5) = POS(4)+1
50    IF (BUFFER(POS(5):POS(5)) .EQ. ' ') THEN
         POS(5) = POS(5) + 1
         GOTO 50
      ENDIF
      POS(6) = POS(5)+1
60    IF (BUFFER(POS(6):POS(6)) .NE. ' ') THEN
         POS(6) = POS(6) + 1
         GOTO 60
      ENDIF
      POS(7) = POS(6)+1
70    IF (BUFFER(POS(7):POS(7)) .EQ. ' ') THEN
         POS(7) = POS(7) + 1
         GOTO 70
      ENDIF
      POS(8) = POS(7)+1
80    IF (BUFFER(POS(8):POS(8)) .NE. ' ') THEN
         POS(8) = POS(8) + 1
         GOTO 80
      ENDIF
      POS(9) = POS(8)+1
90    IF (BUFFER(POS(9):POS(9)) .EQ. ' ') THEN
         POS(9) = POS(9) + 1
         GOTO 90
      ENDIF
      POS(10) = POS(9)+1
100   IF (BUFFER(POS(10):POS(10)) .NE. ' ') THEN
         POS(10) = POS(10) + 1
         GOTO 100
      ENDIF
      POS(11) = POS(10)+1
110   IF (BUFFER(POS(11):POS(11)) .EQ. ' ') THEN
         POS(11) = POS(11) + 1
         GOTO 110
      ENDIF
      POS(12) = POS(11)+1
120   IF (BUFFER(POS(12):POS(12)) .NE. ' ') THEN
         POS(12) = POS(12) + 1
         GOTO 120
      ENDIF
      FORM ='(  X,I  ,  X,I  ,  X,E  .4,  X,E  .4,  X,E  .4,  X,E  .4)'
      WRITE(FORM(2:3),'(I2)') POS(1)-1
      WRITE(FORM(7:8),'(I2)') POS(2)-POS(1)
      WRITE(FORM(10:11),'(I2)') POS(3)-POS(2)
      WRITE(FORM(15:16),'(I2)') POS(4)-POS(3)
      WRITE(FORM(18:19),'(I2)') POS(5)-POS(4)
      WRITE(FORM(23:24),'(I2)') POS(6)-POS(5)
      WRITE(FORM(28:29),'(I2)') POS(7)-POS(6)
      WRITE(FORM(33:34),'(I2)') POS(8)-POS(7)
      WRITE(FORM(38:39),'(I2)') POS(9)-POS(8)
      WRITE(FORM(43:44),'(I2)') POS(10)-POS(9)
      WRITE(FORM(48:49),'(I2)') POS(11)-POS(10)
      WRITE(FORM(53:54),'(I2)') POS(12)-POS(11)
      END
CNEXT
c
c
c
c     interpolation routines added by Krieger, IPP/95
c
      SUBROUTINE linint(xa,ya,n,x,y)
      implicit none
      INTEGER n
      REAL x,y,xa(n),ya(n)
      INTEGER i
      call locate2(xa,n,x,i)
      if (i.eq.0) then
      y=ya(1)
      else if (i.eq.n) then
      y=ya(n)
      else
      y=(ya(i+1)-ya(i))/(xa(i+1)-xa(i))*(x-xa(i))+ya(i)
      endif
      return
      END
c
c
c
      SUBROUTINE linin2(x1a,x2a,ya,m,n,x1,x2,y)
      implicit none
      INTEGER m,n,NMAX,MMAX
      REAL x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=100,MMAX=100)
      INTEGER j,k
      REAL ymtmp(MMAX),yntmp(NMAX)
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
11      continue
        call linint(x2a,yntmp,n,x2,ymtmp(j))
12    continue
      call linint(x1a,ymtmp,m,x1,y)
      return
      END
c
c
c
      SUBROUTINE locate2(xx,n,x,j)
      implicit none
      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
        goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END
c
c
c
      subroutine sortrad(radord,count,top,nizs)
      implicit none
      include 'params'
      include 'cgeom'
      real radord(maxnks*maxnrs,0:maxizs+1,5)
      integer count (0:maxizs+1)
      integer top (0:maxizs+1)
      integer nizs
c
c     SORTRAD: Sorts each IZ slice of the array radord based
c              on the first element.
c
c
      integer ik,ir,iz,comp,in,lastcomp
      real cur
c
      do iz = 0,nizs+1
c
c        Sort each segment of the array.
c
c         if (count(iz).eq.0) cycle
c
         if (count(iz).ne.0) then
c
         top(iz) = 1
c
         do in = 2,count(iz)
c
            cur = radord(in,iz,1)
c
            if (cur.ge.radord(top(iz),iz,1)) then
               radord(in,iz,2) = top(iz)
               radord(top(iz),iz,3) = in
               top(iz) = in
            else
c
               lastcomp = top(iz)
               comp = radord(top(iz),iz,2)
c
100            continue
c
               if (comp.eq.-1) then

                  radord(lastcomp,iz,2) = in
                  radord(in,iz,3) = lastcomp
                  goto 200

               elseif (cur.ge.radord(comp,iz,1)) then
                  radord(in,iz,2) = comp
                  radord(in,iz,3) = radord(comp,iz,3)
                  radord(int(radord(comp,iz,3)),iz,2) = in
                  radord(comp,iz,3) = in
                  goto 200
               endif
c
               lastcomp = comp
               comp = radord(comp,iz,2)
c
               goto 100
c
200            continue
c
            endif
c
         end do

c
c        Endif for count(iz).ne.0
c
         endif

      end do
c
c
c
c      do iz = 0,nizs+1
c
c         write (6,*) 'Ionization state (RAD): ',iz,count(iz),top(iz)
c
c         if (count(iz).gt.0) then
c
c            comp = top(iz)
c
c            ik = radord(num,iz,4)
c            ir = radord(num,iz,5)
c
c            do while(comp.gt.0)
c               write (6,'(3i4,3g13.5)') ik,ir,comp,radord(comp,iz,1)
c     >                           ,knbs(ik,ir),ktebs(ik,ir)
c               comp = radord(comp,iz,2)
c            end do
c
c         end if
c
c      end do
c

      return
      end
c
c
c
      SUBROUTINE INTCALC(R,Z,ROLD,ZOLD,RP,ZP,
     >                TBAC,TFOR,LBAC,
     >                LFOR,RNEW,ZNEW,TNEW,TNORMAL,SECT,nrfopt)
      IMPLICIT NONE
      integer nrfopt
      REAL R,Z,ROLD,ZOLD,RP,ZP,TFOR,TBAC,LFOR,LBAC,RNEW,ZNEW,TNEW
      real tnormal
      LOGICAL SECT
C
C     THIS ROUTINE TESTS WHETHER OR NOT A PARTICLE TRAJECTORY
C     DEFINED BY ITS LAST TWO POSITIONS CROSSES THE SPECIFIED
C     WALL SEGMENT. IF IT DOES THE POINT OF INTERSECTION AND THE
C     REFLECTION ANGLE ARE RETURNED IN THE VARIABLES RNEW, ZNEW
C     AND TNEW. AND SECT IS SET TO TRUE. IF AN INTERSECTION
C     IS NOT FOUND THE VALUES IN THE ABOVE VARIABLES ARE UNDEFINED
C     AND SECT IS SET TO FALSE.
c
c     The new position and reflection angle are based on the specified
c     value of NRFOPT (neutral reflection option).
C
C     DAVID ELDER , SEPT 23 , 1992
C
c     Convert to double precision for internal processing of
c     intersection calculations.
c
      double precision
     >   Rdp,Zdp,ROLDdp,ZOLDdp,RPdp,ZPdp,
     >   TFORdp,TBACdp,LFORdp,LBACdp,RNEWdp,ZNEWdp,TNEWdp,
     >   tnormaldp

      rnewdp = 0.0
      znewdp = 0.0
      tnewdp = 0.0
      tnormaldp = 0.0

      Rdp    = r
      Zdp    = z
      ROLDdp = rold
      ZOLDdp = zold
      RPdp   = rp
      ZPdp   = zp
      TFORdp = tfor
      TBACdp = tbac
      LFORdp = lfor
      LBACdp = lbac

      call intcalcdp(Rdp,Zdp,ROLDdp,ZOLDdp,RPdp,ZPdp,
     >                TBACdp,TFORdp,LBACdp,
     >                LFORdp,RNEWdp,ZNEWdp,TNEWdp,
     >                tnormaldp,SECT,nrfopt)

      RNEW   = rnewdp
      ZNEW   = znewdp
      TNEW   = tnewdp
      tnormal= tnormaldp

      return
      end
c
c
c
      SUBROUTINE INTCALCDP(R,Z,ROLD,ZOLD,RP,ZP,
     >                TBAC,TFOR,LBAC,
     >                LFOR,RNEW,ZNEW,TNEW,Tnormal,SECT,nrfopt)
      IMPLICIT NONE
      integer nrfopt
      REAL*8 R,Z,ROLD,ZOLD,RP,ZP,TFOR,TBAC,LFOR,LBAC,RNEW,
     >       ZNEW,TNEW,tnormal
      REAL ATAN2C
      EXTERNAL ATAN2C
      LOGICAL SECT
C     INCLUDE "PARAMS"
      include    'params'
C
C     THIS ROUTINE TESTS WHETHER OR NOT A PARTICLE TRAJECTORY
C     DEFINED BY ITS LAST TWO POSITIONS CROSSES THE SPECIFIED
C     WALL SEGMENT. IF IT DOES THE POINT OF INTERSECTION AND THE
C     REFLECTION ANGLE ARE RETURNED IN THE VARIABLES RNEW, ZNEW
C     AND TNEW. AND SECT IS SET TO TRUE. IF AN INTERSECTION
C     IS NOT FOUND THE VALUES IN THE ABOVE VARIABLES ARE UNDEFINED
C     AND SECT IS SET TO FALSE.
C
C     DAVID ELDER , SEPT 23 , 1992
C
      real*8 MI,BI,MW,BW ,TNORM,TIMP,r2,z2
      real   dz,dr
      LOGICAL LINX
      integer secti
c
c     IPP/08 Krieger - apparently there are conditions where "tnorm"
c     is not computed. Yet, tnormal is always set to tnorm at the end
c     of the routine -> leads to possible run time error. Fixed by
c     initializing tnorm to 0
      tnorm = 0.0
c
C
c      WRITE (6,*) 'INTCALC:',R,Z,ROLD,ZOLD
C
      IF (R.NE.ROLD) THEN
         MI = (Z-ZOLD)/(R-ROLD)
         BI = Z - R * MI
         LINX = .FALSE.
      ELSE
         MI = 0.0
         BI = R
         LINX = .TRUE.
      ENDIF
c
      dz = z-zold
      dr = r-rold
C
      IF (LFOR.GT.0.0) THEN
c
c         CALL INTSECTDP(TFOR,LFOR,RP,ZP,MI,BI,LINX,SECT,RNEW,ZNEW,
c     >                R,Z,ROLD,ZOLD)
c
         r2 = rp + lfor * cos(tfor)
         z2 = zp + lfor * sin(tfor)
c
         CALL INTSECT2DP(r,z,rold,zold,rp,zp,r2,z2,rnew,znew,secti)
c
         IF (SECTI.eq.1) THEN
            sect = .true. 
            TNORM = TFOR - PI/2.0
            TIMP  = ATAN2C(dz,dr)
            CALL REFANGDP(TNORM,TIMP,TNEW,nrfopt)
         else 
            sect = .false. 
         ENDIF
      ENDIF
      IF ((.NOT.SECT).AND.LBAC.GT.0.0) THEN
c
c         CALL INTSECTDP(TBAC,LBAC,RP,ZP,MI,BI,LINX,SECT,RNEW,ZNEW,
c     >                R,Z,ROLD,ZOLD)
c
         r2 = rp + lbac * cos(tbac)
         z2 = zp + lbac * sin(tbac)
c
         CALL INTSECT2DP(r,z,rold,zold,rp,zp,r2,z2,rnew,znew,secti)
c
         IF (SECTI.eq.1) THEN
            sect=.true.
            TNORM = TBAC + PI/2.0
            TIMP  = ATAN2C(dz,dr)
            CALL REFANGDP(TNORM,TIMP,TNEW,nrfopt)
         else
            sect = .false.
         ENDIF
      ENDIF
C
c      WRITE(6,'(a,11f8.3)')  'INTCALC:',R,Z,ROLD,ZOLD,
c     >              TFOR*180.0/PI,
c     >              TBAC*180.0/PI,rnew,znew,tnew*180.0/PI,
c     >              timp*180.0/PI,tnorm*180.0/PI
C

c
c     Assign value to tnormal  
c
      tnormal = tnorm 
c
      RETURN
      END
c
c
c
      SUBROUTINE INTCALC2(R,Z,ROLD,ZOLD,RS,ZS,
     >                RE,ZE,RNEW,ZNEW,TNEW,TNORMAL,SECT,nrfopt)
      IMPLICIT NONE
      integer nrfopt
      REAL R,Z,ROLD,ZOLD,RP,ZP,RS,ZS,RE,ZE,RNEW,ZNEW,TNEW
      real tnormal
      LOGICAL SECT
C
C     THIS ROUTINE TESTS WHETHER OR NOT A PARTICLE TRAJECTORY
C     DEFINED BY ITS LAST TWO POSITIONS CROSSES THE SPECIFIED
C     WALL SEGMENT. IF IT DOES THE POINT OF INTERSECTION AND THE
C     REFLECTION ANGLE ARE RETURNED IN THE VARIABLES RNEW, ZNEW
C     AND TNEW. AND SECT IS SET TO TRUE. IF AN INTERSECTION
C     IS NOT FOUND THE VALUES IN THE ABOVE VARIABLES ARE UNDEFINED
C     AND SECT IS SET TO FALSE.
c
c     The new position and reflection angle are based on the specified
c     value of NRFOPT (neutral reflection option).
C
C     DAVID ELDER , SEPT 23 , 1992
C
c     Convert to double precision for internal processing of
c     intersection calculations.
c
      double precision
     >   Rdp,Zdp,ROLDdp,ZOLDdp,
     >   rsdp,zsdp,redp,zedp,
     >   RNEWdp,ZNEWdp,TNEWdp,
     >   tnormaldp

      rnewdp = 0.0
      znewdp = 0.0
      tnewdp = 0.0
      tnormaldp = 0.0

      Rdp    = r
      Zdp    = z
      ROLDdp = rold
      ZOLDdp = zold
      Rsdp   = rs
      Zsdp   = zs
      Redp   = re
      Zedp   = ze

      call intcalc2dp(Rdp,Zdp,ROLDdp,ZOLDdp,
     >                rsdp,zsdp,redp,zedp,
     >                RNEWdp,ZNEWdp,TNEWdp,
     >                tnormaldp,SECT,nrfopt)

      RNEW   = rnewdp
      ZNEW   = znewdp
      TNEW   = tnewdp
      tnormal= tnormaldp

      return
      end
c
c
c
      SUBROUTINE INTCALC2DP(R,Z,ROLD,ZOLD,
     >                rs,zs,re,ze,
     >                RNEW,ZNEW,TNEW,Tnorm,SECT,nrfopt)
      IMPLICIT NONE
      integer nrfopt
      REAL*8 R,Z,ROLD,ZOLD,rs,zs,re,ze,RNEW,
     >       ZNEW,TNEW,tnorm
      REAL ATAN2C
      EXTERNAL ATAN2C
      LOGICAL SECT
C     INCLUDE "PARAMS"
      include    'params'
C
C     THIS ROUTINE TESTS WHETHER OR NOT A PARTICLE TRAJECTORY
C     DEFINED BY ITS LAST TWO POSITIONS CROSSES THE SPECIFIED
C     WALL SEGMENT. IF IT DOES THE POINT OF INTERSECTION AND THE
C     REFLECTION ANGLE ARE RETURNED IN THE VARIABLES RNEW, ZNEW
C     AND TNEW. AND SECT IS SET TO TRUE. IF AN INTERSECTION
C     IS NOT FOUND THE VALUES IN THE ABOVE VARIABLES ARE UNDEFINED
C     AND SECT IS SET TO FALSE.
C
C     DAVID ELDER , SEPT 23 , 1992
C
      real*8 TIMP
      real   dz,dr
      integer secti
c
c     Surface normal angle
c
      dz=ze-zs
      dr=re-rs
c
      tnorm = atan2c(dz,dr) - PI/2.0      
c
c     Angle of particle trajectory 
c
      dz=z-zold
      dr=r-rold
c
      timp = atan2c(dz,dr)
C
      CALL INTSECT2DP(r,z,rold,zold,rs,zs,re,ze,rnew,znew,secti)
c
      IF (SECTI.eq.1) THEN
         sect = .true. 
         CALL REFANGDP(TNORM,TIMP,TNEW,nrfopt)
      else 
         sect = .false. 
      ENDIF
C
c      WRITE(6,'(a,11f8.3)')  'INTCALC:',R,Z,ROLD,ZOLD,
c     >              TFOR*180.0/PI,
c     >              TBAC*180.0/PI,rnew,znew,tnew*180.0/PI,
c     >              timp*180.0/PI,tnorm*180.0/PI
C

      RETURN
      END
C
C
C
      SUBROUTINE INTSECTDP(THETA,DIST,RP,ZP,MI,BI,LINX,SECT,RNEW,ZNEW,
     >                   R,Z,ROLD,ZOLD)
      IMPLICIT NONE
      real*8 THETA,DIST,RP,ZP,MI,BI,RNEW,ZNEW
      real*8 r,z,rold,zold
      LOGICAL LINX,SECT
C     INCLUDE "PARAMS"
      include    'params'
c
c     EPS is required since the boundary finding routine GA15B will sometimes
c     return a result of exactly zero when the particle may be very close to
c     the boundary. 
c
      real*8 eps
      parameter (eps=1.0d-6)

C
C     THIS ROUTINE TAKES AS INPUT THE START POINT LENGTH AND
C     ANGLE OF A LINE SEGMENT AND THE EQUATION OF A LINE. IT THEN
C     CHECKS TO SEE IF THE TWO INTERSECT BETWEEN THE END-POINTS
C     OF THE LINE SEGMENT - IT THEN RETURNS A FLAG INDICATING THE
C     INTERSECTION AND THE CO-ORDINATES OF THE INTERSECTION POINT.
C     (IF IT EXISTS)
C
C     DAVID ELDER, SEPT 23, 1992
C
      REAL*8 MW,BW,YSTART,YEND,YTEST,xstart,xend,xtest
c
      IF (ABS(THETA).NE.(PI/2.0)) THEN
         MW = TAN(THETA)
         BW = ZP - RP * MW
c
         YSTART = ZP
         YEND = ZP + DIST * SIN(THETA)
         XSTART = RP
         XEND = RP + DIST * COS(THETA)
c
         IF (LINX) THEN
            YTEST = MW * BI + BW
            XTEST = BI
c
c              Requires both R and Z tests in order to deal with 
c              horizontal wall segments with MW=0.0
c
            IF ((ABS(YTEST-YEND).LE.ABS(YSTART-YEND)+eps) .AND.
     >          (ABS(YTEST-YSTART).LE.ABS(YSTART-YEND)+eps).and.
     >          (ABS(XTEST-XEND).LE.ABS(XSTART-XEND)+eps) .AND.
     >          (ABS(XTEST-XSTART).LE.ABS(XSTART-XEND)+eps)) THEN
c             
               SECT = .TRUE.
               RNEW = BI
               ZNEW = YTEST
            ELSE
               SECT = .FALSE.
            ENDIF
         ELSEIF (.NOT.LINX) THEN
            IF (MW.EQ.MI) THEN
               SECT = .FALSE.
            ELSE

               XTEST = (BI - BW) / (MW - MI)
               YTEST = MW * RNEW + BW
c
c               RNEW = (BI - BW) / (MW - MI)
c               ZNEW = MW * RNEW + BW
c
c              Requires both R and Z tests in order to deal with 
c              horizontal wall segments with MW=0.0
c
               IF ((ABS(YTEST-YEND).LE.ABS(YSTART-YEND)+eps) .AND.
     >            (ABS(YTEST-YSTART).LE.ABS(YSTART-YEND)+eps).and.
     >            (ABS(XTEST-XEND).LE.ABS(XSTART-XEND)+eps) .AND.
     >            (ABS(XTEST-XSTART).LE.ABS(XSTART-XEND)+eps)) THEN

                  SECT = .TRUE.
                  RNEW = XTEST
                  ZNEW = YTEST
               ELSE
                   SECT = .FALSE.
               ENDIF
            ENDIF
         ENDIF
      ELSE
c
         BW = RP
         YSTART = ZP
         YEND = ZP + SIGN(DIST,THETA)
c
         IF (LINX) THEN
            IF (BI.EQ.BW) THEN
               RNEW = RP
               ZNEW = ZP
            ELSE
               SECT = .FALSE.
            ENDIF
         ELSEIF (.NOT.LINX) THEN
            YTEST = MI * BW + BI
c
            IF ((ABS(YTEST-YEND).LE.ABS(YSTART-YEND)+eps) .AND.
     >           (ABS(YTEST-YSTART).LE.ABS(YSTART-YEND)+eps)) THEN
c
               SECT = .TRUE.
               RNEW = BW
               ZNEW = YTEST
            ELSE
               SECT = .FALSE.
            ENDIF
         ENDIF
      ENDIF
c
c     Check to see if the intersection point found also lies
c     between R,Z and ROLD,ZOLD - if not - change sect to false.
c

      if (sect.and.
     >    (
     >   ((r.eq.rold).and.
     >   (abs(z-znew).gt.abs(z-zold)+eps.or.
     >    abs(zold-znew).gt.abs(z-zold)+eps))
     >        .or.
     >   ((z.eq.zold).and.
     >   (abs(r-rnew).gt.abs(r-rold).or.
     >    abs(rold-rnew).gt.abs(r-rold)+eps))
     >        .or.
     >   ((r.ne.rold.and.z.ne.zold).and.
     >   ((abs(z-znew).gt.abs(z-zold)+eps.or.
     >     abs(zold-znew).gt.abs(z-zold)+eps)
     >                        .or.
     >   (abs(r-rnew).gt.abs(r-rold)+eps.or.
     >    abs(rold-rnew).gt.abs(r-rold)+eps)))
     >     )) then
c
         write (6,'(a,3l4,12(1x,g16.10))')
     >           'FIX INTSECT ERROR:',sect,linx,(abs(theta).eq.PI/2.0),
     >                  r,z,rnew,znew,
     >                  rold,zold,zp,ystart,yend,theta,PI/2.0

         write (6,'(a,10(1x,g16.10))')
     >           'FIX INTSECT DIFFS:',
     >        dabs(r-rnew),dabs(rold-rnew),dabs(r-rold),
     >        dabs(z-znew),dabs(zold-znew),dabs(z-zold),
     >        dabs(znew-ystart),dabs(znew-yend),dabs(ystart-yend)
c
         sect = .false.
c
      endif

         write (6,'(a,3l4,4x,1p,12(1x,g12.6))')
     >           'INTSECT DBG1:',sect,linx,(abs(theta).eq.PI/2.0),
     >                  r,z,rnew,znew,
     >                  rold,zold,xtest,ytest,theta,PI/2.0

         write (6,'(a,4l5,1p,12(1x,g12.6))')
     >           'INTSECT DBG2:',
     > dabs(r-rnew).le.dabs(r-rold).and.dabs(rold-rnew).le.dabs(r-rold),
     > dabs(z-znew).le.dabs(z-zold).and.dabs(zold-znew).le.dabs(z-zold),
     > dabs(xtest-xstart).le.dabs(xstart-xend).and.
     >                         dabs(xtest-xend).le.dabs(xstart-xend),
     > dabs(ytest-ystart).le.dabs(ystart-yend).and.
     >                         dabs(ytest-yend).le.dabs(ystart-yend),
     >        dabs(r-rnew),dabs(rold-rnew),dabs(r-rold),
     >        dabs(z-znew),dabs(zold-znew),dabs(z-zold),
     >        dabs(ytest-ystart),dabs(ytest-yend),dabs(ystart-yend),
     >        dabs(xtest-xstart),dabs(xtest-xend),dabs(xstart-xend)
C
c      WRITE (6,*) 'INTSECT:',RNEW,ZNEW,SECT,LINX
C
      RETURN
      END
C
C
C
      SUBROUTINE REFANGDP(TNORM,TIMP,TREF,nrfopt)
      IMPLICIT NONE
      integer nrfopt
      REAL*8 TNORM,TIMP,TREF
C     INCLUDE "PARAMS"
      include    'params'
C
C     THIS ROUTINE CALCULATES THE ANGLE OF REFLECTION FROM THE
C     SURFACE. THE THREE ANGLES ARE ALL RELATIVE TO THE POSITIVE
C     R-AXIS.
C              TNORM - SURFACE NORMAL
C              TIMP  - ANGLE OF IMPACT
C              TREF  - CALCULATED REFLECTION ANGLE
C
c
       real   ran1,ran2
       real*8 seed
       real :: n(2),ri(2),rr(2)
       real, external :: atan2c
c slmod begin
       logical first_warning
       data    first_warning /.true./
c slmod end
c
c       real*8 tnorm_tmp,timp_tmp,dthe,tref_tmp
c
c
c       jdemod - Set the maximum reflection angle to 89 degrees
c                from normal for the random angle cases  
c
       if (nrfopt.eq.0.or.nrfopt.eq.1) then
c slmod begin        
         if (.true.) then 
           if (first_warning) then
             call wn('REFANGDP','Using old TREF calculation')
             first_warning = .false.
           endif
           TREF = 2.0 * TNORM - SIGN((PI-ABS(TIMP)),-TIMP)
         else
c
c          jdemod - formula only works under specific circumstances
c                 - replace with general vector formulation
c
           n(1) = cos(tnorm)
           n(2) = sin(tnorm)
           ri(1) = -cos(timp)
           ri(2) = -sin(timp)
          
           rr = ri - 2.0 * n * dot_product(ri,n)
c
c          Calculate reflection angle from reflection vector
c
           tref = atan2c(rr(2),rr(1))
         endif
c
cc          TREF = 2.0 * TNORM - SIGN((PI-ABS(TIMP)),-TIMP)
cc          TREF = TNORM + (TNORM - SIGN((PI-ABS(TIMP)),-TIMP))
cc
cc          jdemod - formula only works under specific circumstances
cc                 - replace with general vector formulation
c
c           n(1) = cos(tnorm)
c           n(2) = sin(tnorm)
c           ri(1) = -cos(timp)
c           ri(2) = -sin(timp)
c
c           rr = ri - 2.0 * n * dot_product(ri,n)
cc
cc          Calculate reflection angle from reflection vector
cc
c           tref = atan2c(rr(2),rr(1))
c slmod end
       elseif (nrfopt.eq.2) then
c
          CALL SURAND2 (SEED, 1, RAN1)
          CALL SURAND2 (SEED, 1, RAN2)
c
          TREF = tnorm + SIGN (min(ACOS (ran1),DEGRAD*89.0), ran2-0.5)
c
c ammod - added options 3 and 4 to NRFOPT
c
       elseif (nrfopt.eq.3) then
c
	  CALL SURAND2 (SEED, 1, RAN1)
          CALL SURAND2 (SEED, 1, RAN2)
c
          TREF = tnorm + SIGN (min(ASIN (ran1),DEGRAD*89.0), ran2-0.5)
c
       elseif (nrfopt.eq.4) then
c
	  CALL SURAND2 (SEED, 1, RAN1)
          CALL SURAND2 (SEED, 1, RAN2)
c
          TREF = tnorm + SIGN (min(ASIN (SQRT(ran1)),degrad*89.0),
     >                         ran2-0.5)
c
c ammod
c
       endif
c

c
c      TREF = 2.0 * TNORM - SIGN((PI-ABS(TIMP)),-TIMP)
c
c      Adjust TNORM
c
c       tnorm_tmp = tnorm
c       if (tnorm_tmp.lt.0.0) tnorm_tmp = tnorm_tmp + 2.0 * PI
c       if (tnorm_tmp.ge.2.0*PI) tnorm_tmp = tnorm_tmp - 2.0 * PI
c
c      Adjust TIMP
c
c       timp_tmp = timp
c       timp_tmp = timp_tmp + PI
c       if (timp_tmp.lt.0.0) timp_tmp = timp_tmp + 2.0 * PI
c       if (timp_tmp.ge.2.0*PI) timp_tmp = timp_tmp - 2.0 * PI
c
c       if (timp_tmp.gt.(tnorm_tmp+PI/2.0)) timp_tmp = timp_tmp - 2.0*PI
c       if (timp_tmp.lt.(tnorm_tmp-PI/2.0)) timp_tmp = timp_tmp + 2.0*PI
c
c       dthe = tnorm_tmp-timp_tmp
c
c       TREF = TNORM_tmp + dthe
C
c
c        WRITE(6,'(a,i6,8(1x,g13.6))') 'REFANG:TREF:',nrfopt,TREF*raddeg,
c     >                 tnorm*raddeg,timp*raddeg
c
c
c     >             tnorm_tmp*raddeg,timp_tmp*raddeg,dthe*raddeg,
c     >             tref_tmp*raddeg
c      if (tref.ne.tref_tmp)
c     >        write(6,'(a,2(1x,g13.6))') 'TREF:',
c     >        tref*raddeg,tref_tmp*raddeg
c
C
      RETURN
      END

C
C
C
c      SUBROUTINE INTSECT(THETA,DIST,RP,ZP,MI,BI,LINX,SECT,RNEW,ZNEW,
c     >                   R,Z,ROLD,ZOLD)
c      IMPLICIT NONE
c      real THETA,DIST,RP,ZP,MI,BI,RNEW,ZNEW
c      real r,z,rold,zold
c      LOGICAL LINX,SECT
C     INCLUDE "PARAMS"
c      include    'params'
C
C     THIS ROUTINE TAKES AS INPUT THE START POINT LENGTH AND
C     ANGLE OF A LINE SEGMENT AND THE EQUATION OF A LINE. IT THEN
C     CHECKS TO SEE IF THE TWO INTERSECT BETWEEN THE END-POINTS
C     OF THE LINE SEGMENT - IT THEN RETURNS A FLAG INDICATING THE
C     INTERSECTION AND THE CO-ORDINATES OF THE INTERSECTION POINT.
C     (IF IT EXISTS)
C
C     DAVID ELDER, SEPT 23, 1992
C
c      REAL MW,BW,YSTART,YEND,YTEST
c
c      IF (ABS(THETA).NE.(PI/2.0)) THEN
c         MW = TAN(THETA)
c         BW = ZP - RP * MW
c         YSTART = ZP
c         YEND = ZP + DIST * SIN(THETA)
c         IF (LINX) THEN
c            YTEST = MW * BI + BW
c            IF ((ABS(YTEST-YEND).LE.ABS(YSTART-YEND)) .AND.
c     >          (ABS(YTEST-YSTART).LE.ABS(YSTART-YEND))) THEN
c               SECT = .TRUE.
c               RNEW = BI
c               ZNEW = YTEST
c            ELSE
c               SECT = .FALSE.
c            ENDIF
c         ELSEIF (.NOT.LINX) THEN
c            IF (MW.EQ.MI) THEN
c               SECT = .FALSE.
c            ELSe
c
c               RNEW = (BI - BW) / (MW - MI)
c               ZNEW = MW * RNEW + BW
c
c               IF ((ABS(ZNEW-YEND).LE.ABS(YSTART-YEND)) .AND.
c     >            (ABS(ZNEW-YSTART).LE.ABS(YSTART-YEND))) THEN
c                   SECT = .TRUE.
c               ELSE
c                   SECT = .FALSE.
c               ENDIF
c            ENDIF
c         ENDIF
c      ELSE
c         BW = RP
c         YSTART = ZP
c         YEND = ZP + SIGN(DIST,THETA)
c         IF (LINX) THEN
c            IF (BI.EQ.BW) THEN
c               RNEW = RP
c               ZNEW = ZP
c            ELSE
c               SECT = .FALSE.
c            ENDIF
c         ELSEIF (.NOT.LINX) THEN
c            YTEST = MI * BW + BI
c            IF ((ABS(YTEST-YEND).LE.ABS(YSTART-YEND)) .AND.
c     >           (ABS(YTEST-YSTART).LE.ABS(YSTART-YEND))) THEN
c               SECT = .TRUE.
c               RNEW = BW
c               ZNEW = YTEST
c            ELSE
c               SECT = .FALSE.
c            ENDIF
c         ENDIF
c      ENDIF
c
c     Check to see if the intersection point found also lies
c     between R,Z and ROLD,ZOLD - if not - change sect to false.
c
c
c      if (sect.and.
c     >    (
c     >   ((r.eq.rold).and.
c     >   (abs(z-znew).gt.abs(z-zold).or.abs(zold-znew).gt.abs(z-zold)))
c     >        .or.
c     >   ((z.eq.zold).and.
c     >   (abs(r-rnew).gt.abs(r-rold).or.abs(rold-rnew).gt.abs(r-rold)))
c     >        .or.
c     >   ((r.ne.rold.and.z.ne.zold).and.
c     >   ((abs(z-znew).gt.abs(z-zold).or.abs(zold-znew).gt.abs(z-zold))
c     >                        .or.
c     >   (abs(r-rnew).gt.abs(r-rold).or.abs(rold-rnew).gt.abs(r-rold))))
c     >     )) then
c
c         sect = .false.
c
c         write (6,'(a,l4,6(1x,g18.10))')
c     >           'FIX INTSECT ERROR:',sect,r,z,rnew,znew,
c     >                  rold,zold
c
c      endif
c
C
C      WRITE (6,*) 'INTSECT:',RNEW,ZNEW,SECT,LINX
C
c      RETURN
c      END
C
C
C
c      SUBROUTINE REFANG(TNORM,TIMP,TREF)
c      IMPLICIT NONE
c      REAL TNORM,TIMP,TREF
C     INCLUDE "PARAMS"
c      include    'params'
C
C     THIS ROUTINE CALCULATES THE ANGLE OF REFLECTION FROM THE
C     SURFACE. THE THREE ANGLES ARE ALL RELATIVE TO THE POSITIVE
C     R-AXIS.
C              TNORM - SURFACE NORMAL
C              TIMP  - ANGLE OF IMPACT
C              TREF  - CALCULATED REFLECTION ANGLE
C
c      TREF = 2.0 * TNORM - SIGN((PI-ABS(TIMP)),-TIMP)
c
c       real tnorm_tmp,timp_tmp,dthe,tref_tmp
c
c       TREF = 2.0 * TNORM - SIGN((PI-ABS(TIMP)),-TIMP)
c
c      Adjust TNORM
c
c       tnorm_tmp = tnorm
c       if (tnorm_tmp.lt.0.0) tnorm_tmp = tnorm_tmp + 2.0 * PI
c       if (tnorm_tmp.ge.2.0*PI) tnorm_tmp = tnorm_tmp - 2.0 * PI
c
c      Adjust TIMP
c
c       timp_tmp = timp
c       timp_tmp = timp_tmp + PI
c       if (timp_tmp.lt.0.0) timp_tmp = timp_tmp + 2.0 * PI
c       if (timp_tmp.ge.2.0*PI) timp_tmp = timp_tmp - 2.0 * PI
c
c       if (timp_tmp.gt.(tnorm_tmp+PI/2.0)) timp_tmp = timp_tmp - 2.0*PI
c       if (timp_tmp.lt.(tnorm_tmp-PI/2.0)) timp_tmp = timp_tmp + 2.0*PI
c
c       dthe = tnorm_tmp-timp_tmp
c
c       TREF = TNORM_tmp + dthe
C
c      WRITE(6,'(a,8(1x,g13.6))') 'REFANG:TREF:',TREF*raddeg,
c     >                 tnorm*raddeg,timp*raddeg,
c     >             tnorm_tmp*raddeg,timp_tmp*raddeg,dthe*raddeg,
c     >             tref_tmp*raddeg
c      if (tref.ne.tref_tmp)
c     >        write(6,'(a,2(1x,g13.6))') 'TREF:',
c     >        tref*raddeg,tref_tmp*raddeg
c
C
c      RETURN
c      END
c
c
c
c
      subroutine calc_jhfactors(jhtots)
      implicit none
      include 'params'
      include 'pindata'
      include 'cgeom'
c
      real jhtots(10,3)
c
c     CALC_JHFACTORS: This routine calculates a summary of the
c                     Johnson-Hinov factor for various regions of the 
c                     grid.
c
c                     The Johnson-Hinov factor is the ratio of dalpha photons
c                     emitted to ionizations.  
c
c     Local variables
c        
      integer ik,ir,in
c
      real n_photons(12),n_ionizations(12)
c
c      real    dalpha_energy,dalpha_lambda,h,c
c      data    dalpha_lambda /6561.9e-10/   
c      data    h /6.626e-34/
c      data    c /2.998e8/ 
c
c     Calculate energy of Dalpha photon 
c      
c      dalpha_energy = (h*c)/dalpha_lambda 
c
      call rzero(jhtots,10*3)
      call rzero(n_photons,12)
      call rzero(n_ionizations,12)
c
      do ir = 1,nrs
c        
         do ik = 1,nks(ir)
c
c           Core
c
            if (ir.lt.irsep) then 
c
               n_photons(1)     = n_photons(1) 
     >                            + pinalpha(ik,ir)
     >                              * kareas(ik,ir)  
               n_ionizations(1) = n_ionizations(1) 
     >                            + pinion(ik,ir) * kareas(ik,ir)  
c
c           Main SOL
c
            elseif (ir.le.irwall) then 

c
c              First half of ring - inner leg on DIIID grids (xpoint down)
c
               if (ik.lt.nks(ir)/2) then 
c
c                 Divertor Region
c 
                  if (ik.lt.ikto) then 
c
                     n_photons(2)     = n_photons(2) 
     >                            + pinalpha(ik,ir)
     >                                * kareas(ik,ir)  
                     n_ionizations(2) = n_ionizations(2) 
     >                            + pinion(ik,ir) * kareas(ik,ir)  
c
c                 Main Vessel region
c
                  else  
c
                     n_photons(3)     = n_photons(3) 
     >                            + pinalpha(ik,ir)
     >                              * kareas(ik,ir)  
                     n_ionizations(3) = n_ionizations(3) 
     >                            + pinion(ik,ir) * kareas(ik,ir)  
                  endif
c
c              Second half of ring
c
               else
c
c                 Divertor Region
c 
                  if (ik.gt.ikti) then 
c
                     n_photons(5)     = n_photons(5) 
     >                            + pinalpha(ik,ir)
     >                              * kareas(ik,ir)  
                     n_ionizations(5) = n_ionizations(5) 
     >                            + pinion(ik,ir) * kareas(ik,ir)  
c
c                 Main Vessel region
c
                  else  
c
                     n_photons(6)     = n_photons(6) 
     >                            + pinalpha(ik,ir)
     >                              * kareas(ik,ir)  
                     n_ionizations(6) = n_ionizations(6) 
     >                            + pinion(ik,ir) * kareas(ik,ir)  
                  endif

               endif
c
c           PFZ
c
            else
c
               n_photons(9)     = n_photons(9) 
     >                            + pinalpha(ik,ir)
     >                              * kareas(ik,ir)  
               n_ionizations(9) = n_ionizations(9) 
     >                            + pinion(ik,ir) * kareas(ik,ir)  
c
            endif
c
         end do
c
      end do 
c
c     Summarize certain regions        
c
c     First half of SOL  
c      
      n_photons(4) = n_photons(2)+n_photons(3)
      n_ionizations(4) = n_ionizations(2)+n_ionizations(3)
c
c     Second half of SOL
c         
      n_photons(7) = n_photons(5)+n_photons(6)
      n_ionizations(7) = n_ionizations(5)+n_ionizations(6)
c
c     ALL of SOL
c
      n_photons(8) = n_photons(4)+n_photons(7)
      n_ionizations(8) = n_ionizations(4)+n_ionizations(7)
c
c     ALL of plasma 
c
      n_photons(10) = n_photons(1)+n_photons(8)+n_photons(9)
      n_ionizations(10) = n_ionizations(1)+n_ionizations(8)
     >                    +n_ionizations(9)
c
c     Assign values to jhtots
c
      do in = 1,10
c
         jhtots(in,1) = n_photons(in) 
         jhtots(in,2) = n_ionizations(in) 
c
         if (n_ionizations(in).ne.0.0) then 
            jhtots(in,3) = n_ionizations(in)/n_photons(in)
         else
            jhtots(in,3) = 0.0
         endif 
c
      end do
c
      return
      end 
c
c
c
      subroutine calc_jh2d(jh2d)
      implicit none
      include 'params'
      include 'pindata'
      include 'cgeom'
c
      real jh2d(maxnks,maxnrs)
c
c     CALC_JHFACTORS: This routine calculates the Johnson-Hinov factor for
c                     each cell on the grid.
c
c                     The Johnson-Hinov factor is the ratio of dalpha photons
c                     emitted to ionizations.  
c
c     Local variables
c        
      integer ik,ir
c
c      real    dalpha_energy,dalpha_lambda,h,c
c      data    dalpha_lambda /6561.9e-10/   
c      data    h /6.626e-34/
c      data    c /2.998e8/ 
c
c     Calculate energy of Dalpha photon 
c      
c      dalpha_energy = (h*c)/dalpha_lambda 
c
      do ir = 1,nrs
c        
         do ik = 1,nks(ir)
c
c           Calculate 2D array
c
            if (pinion(ik,ir).ne.0.0) then 
               jh2d(ik,ir) = pinion(ik,ir)/pinalpha(ik,ir)
            else
               jh2d(ik,ir) = 0.0
            endif
c
         end do 
c
      end do  
c
      return
      end
c
c
c
      subroutine pr_calc_walldep
      implicit none
c
c     Adds contributions from wall element wallin - contained in the 
c     specified region to the totals in walldep.  
c          
      include 'params'
      include 'dynam3'
      include 'comtor' 
      include 'printopt' 
c
      integer in,ik,start_region
      real walldep(4,5)
      real walldep_i(4,5)
c
      call rzero(walldep,4*5) 
      call rzero(walldep_i,4*5) 
c
c     Need total source from each segment as well as total ion source from
c     each segment.  
c
      do ik = 1,wallpts
c
c        Starting on main wall 
c
         if (ik.ge.wlwall1.and.ik.le.wlwall2) then 
            start_region = 1  
c	 
c        Starting on PFZ wall
c	 
         elseif (ik.ge.wltrap1.and.ik.le.wltrap2) then 
            start_region = 2
c	 
c        Starting on target 1 
c	 
         elseif (ik.ge.(wlwall2+1).and.ik.le.(wltrap1-1)) then 
            start_region = 3
c	 
c        Starting on target 2
c	 
         elseif (ik.ge.(wltrap2+1).and.ik.le.wallpts) then 
            start_region = 4
         endif
c
c        Total erosion (ion+neutral) 
c
         walldep(start_region,5) = walldep(start_region,5) + wallse(ik)
         walldep_i(start_region,5) = walldep_i(start_region,5) +
     >                                       wallse_i(ik)

c
c        Loop over wall summing up end regions 
c
         do in = 1,wallpts 
c
c           Ending on main wall 
c
            if (in.ge.wlwall1.and.in.le.wlwall2) then 
               walldep(start_region,1) = walldep(start_region,1) 
     >                                 + wtdep(ik,in,3) 
               walldep_i(start_region,1) = walldep_i(start_region,1) 
     >                                 + wtdep(ik,in,1) 
c
c           Ending on PFZ wall
c
            elseif (in.ge.wltrap1.and.in.le.wltrap2) then 
               walldep(start_region,2) = walldep(start_region,2) 
     >                                 + wtdep(ik,in,3) 
               walldep_i(start_region,2) = walldep_i(start_region,2) 
     >                                 + wtdep(ik,in,1) 
c
c           Ending on target 1 
c
            elseif (in.ge.(wlwall2+1).and.in.le.(wltrap1-1)) then 
               walldep(start_region,3) = walldep(start_region,3) 
     >                                 + wtdep(ik,in,3) 
               walldep_i(start_region,3) = walldep_i(start_region,3) 
     >                                 + wtdep(ik,in,1) 
c
c           Ending on target 2
c
            elseif (in.ge.(wltrap2+1).and.in.le.wallpts) then 
               walldep(start_region,4) = walldep(start_region,4) 
     >                                 + wtdep(ik,in,3) 
               walldep_i(start_region,4) = walldep_i(start_region,4) 
     >                                 + wtdep(ik,in,1) 
            endif
c
         end do
      end do 

c
c     Print out raw deposition and erosion numbers 
c

      call prb
      call prb
      call prc(' Summary of Total Erosion/Deposition:') 
      call prb
c
      call prc('                      Deposition Region ')
      call prc(' Source           Main      PFZ     '
     >                  //INNER//'    '//OUTER//'   Total')
      call prc(' Region           Wall      Wall   '//
     >                  'Target   Target   Erosion') 
c
      write(coment,'(a12,2x,5(1x,f8.1))') 'Main Wall',
     >    (walldep(1,in),in=1,4),walldep(1,5)
      call prc(coment)
      write(coment,'(a12,2x,5(1x,f8.1))') 'PFZ Wall',
     >    (walldep(2,in),in=1,4),walldep(2,5)
      call prc(coment)
      write(coment,'(a12,2x,5(1x,f8.1))') INNER//' Target',
     >    (walldep(3,in),in=1,4),walldep(3,5)
      call prc(coment)
      write(coment,'(a12,2x,5(1x,f8.1))') OUTER//' Target',
     >    (walldep(4,in),in=1,4),walldep(4,5)
      call prc(coment)
c
      write(coment,'(a12,2x,5(1x,f8.1))') 'Total',
     > (walldep(1,1)+walldep(2,1)+walldep(3,1)+walldep(4,1)), 
     > (walldep(1,2)+walldep(2,2)+walldep(3,2)+walldep(4,2)), 
     > (walldep(1,3)+walldep(2,3)+walldep(3,3)+walldep(4,3)), 
     > (walldep(1,4)+walldep(2,4)+walldep(3,4)+walldep(4,4)), 
     > (walldep(1,5)+walldep(2,5)+walldep(3,5)+walldep(4,5))
     

      call prc(coment)
        

      call prb
      call prb
      call prc(' Summary of Ionized Erosion/Deposition:') 
c
      call prb
c
      call prc('                      Deposition Region ')
      call prc(' Source           Main      PFZ     '
     >                  //INNER//'    '//OUTER//'   Total')
      call prc(' Region           Wall      Wall   '//
     >                  'Target   Target   Erosion') 
c
      write(coment,'(a12,2x,5(1x,f8.1))') 'Main Wall',
     >    (walldep_i(1,in),in=1,4),walldep_i(1,5)
      call prc(coment)
      write(coment,'(a12,2x,5(1x,f8.1))') 'PFZ Wall',
     >    (walldep_i(2,in),in=1,4),walldep_i(2,5)
      call prc(coment)
      write(coment,'(a12,2x,5(1x,f8.1))') INNER//' Target',
     >    (walldep_i(3,in),in=1,4),walldep_i(3,5)
      call prc(coment)
      write(coment,'(a12,2x,5(1x,f8.1))') OUTER//' Target',
     >    (walldep_i(4,in),in=1,4),walldep_i(4,5)
      call prc(coment)
c
      write(coment,'(a12,2x,5(1x,f8.1))') 'Total',
     > (walldep_i(1,1)+walldep_i(2,1)
     > +walldep_i(3,1)+walldep_i(4,1)), 
     > (walldep_i(1,2)+walldep_i(2,2)
     > +walldep_i(3,2)+walldep_i(4,2)), 
     > (walldep_i(1,3)+walldep_i(2,3)
     > +walldep_i(3,3)+walldep_i(4,3)), 
     > (walldep_i(1,4)+walldep_i(2,4)
     > +walldep_i(3,4)+walldep_i(4,4)), 
     > (walldep_i(1,5)+walldep_i(2,5)
     > +walldep_i(3,5)+walldep_i(4,5))
     
      call prc(coment)


c
c     Normalize each line of walldep into probabilities
c
      do ik = 1,4
c
         do in = 1,4
c
            if (walldep(ik,5).gt.0.0) then 
               walldep(ik,in) = walldep(ik,in) 
     >                     / walldep(ik,5)
            else 
               walldep(ik,in) = 0.0
            endif

            if (walldep_i(ik,5).gt.0.0) then 
               walldep_i(ik,in) = walldep_i(ik,in) 
     >                     / walldep_i(ik,5)
            else 
               walldep_i(ik,in) = 0.0
            endif
c
         end do 
c
      end do
c
c     Print out table of results - all deposition
c
      call prb
      call prb
      call prc(' Summary of Total Deposition Probabilities:') 
      call prb
c
      call prc('                      Deposition Region ')
      call prc(' Source          Main     PFZ    '
     >                  //INNER//'   '//OUTER)
      call prc(' Region          Wall     Wall  '//
     >                  'Target  Target') 
c
      write(coment,'(a12,2x,4(1x,f7.4))') 'Main Wall',
     >                 (walldep(1,in),in=1,4)
      call prc(coment)
      write(coment,'(a12,2x,4(1x,f7.4))') 'PFZ Wall',
     >                 (walldep(2,in),in=1,4)
      call prc(coment)
      write(coment,'(a12,2x,4(1x,f7.4))') INNER//' Target',
     >                 (walldep(3,in),in=1,4)
      call prc(coment)
      write(coment,'(a12,2x,4(1x,f7.4))') OUTER//' Target',
     >                 (walldep(4,in),in=1,4)
      call prc(coment)

c
c     Print out table of results - ion deposition
c
      call prb
      call prb
      call prc(' Summary of Ion Deposition Probabilities:') 
      call prb
c
      call prc('                      Deposition Region ')
      call prc(' Source          Main     PFZ    '
     >                  //INNER//'   '//OUTER)
      call prc(' Region          Wall     Wall  '//
     >                  'Target  Target') 
c
      write(coment,'(a12,2x,4(1x,f7.4))') 'Main Wall',
     >                 (walldep_i(1,in),in=1,4)
      call prc(coment)
      write(coment,'(a12,2x,4(1x,f7.4))') 'PFZ Wall',
     >                 (walldep_i(2,in),in=1,4)
      call prc(coment)
      write(coment,'(a12,2x,4(1x,f7.4))') INNER//' Target',
     >                 (walldep_i(3,in),in=1,4)
      call prc(coment)
      write(coment,'(a12,2x,4(1x,f7.4))') OUTER//' Target',
     >                 (walldep_i(4,in),in=1,4)
      call prc(coment)

c
      return 
      end   












c
c
c
      subroutine  Read_AdditionalPlotData(buffer)
      implicit none
      character*(*) buffer
      include 'params'
      include 'outcom'
c
c     Read_AdditionalPlotData:
c
c     This routine reads in addtional optional plot data that
c     is defined in the input file by a "#" followed by a two
c     character designation. Each option is described below - an
c     unexpected option will generate a message but not an error. 
c
      integer tmp_nsets,in
c      integer lenstr
c      external lenstr
      character*80 graph
c
c
c     #01 - read in contour plot scaling data
c
      if (buffer(2:4).eq.'#01') then   

         READ (BUFFER,*,ERR=9999,END=9999) graph,
     >                      minscale,maxscale,localcngs
c
c     #02 - read in contour plot zoom data
c
      elseif (buffer(2:4).eq.'#02') then   

         READ(buffer,*,err=9999,end=9999) graph,
     >                        xcen,ycen,xnear2,ynear2

c
c     #03 - reads in a list of experimental datasets to be included
c           on the plot if possible. 
c

      elseif (buffer(2:4).eq.'#03') then   

         READ(buffer,*,err=9999,end=9999) graph,
     >                        tmp_nsets
c
c        Limit number of datasets to max specified in code
c 
         expt_nsets= min(tmp_nsets,max_expt_datasets)
c
         if (expt_nsets.gt.0) then 

            READ(buffer,*,err=9999,end=9999) graph,
     >          tmp_nsets,(expt_datasets(in),in=1,expt_nsets)
       
         else
            expt_nsets = 0
         endif
c
c     Issue message for unrecognized option 
c
      else 
c

         len = lenstr(buffer) 
         write(6,'(a,a)') 'Unrecognized Optional Plot Input Line:',
     >                    buffer(1:len)
         write(0,'(a,a)') 'Unrecognized Optional Plot Input Line:',
     >                    buffer(1:len)


      endif

      return 
c
c     Trap I/O errors and issue some message
c

 9999 len = lenstr(buffer)

      write(6,'(a,a)') 'ERROR reading additional plot data:',
     >                 buffer(1:len)
      write(0,'(a,a)') 'ERROR reading additional plot data:',
     >                 buffer(1:len)
     
      return 
      end 


C
C
C
      SUBROUTINE RDG1 (GRAPH,ADASID,ADASYR,ADASEX,
     >                 ISELE,ISELR,ISELX,ISELD,IERR)
      implicit none
      INTEGER   ISELE,ISELR,ISELX,ISELD,IERR,ADASYR
      CHARACTER GRAPH*(*), ADASID*(*),ADASEX*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RDG1 : READ IN SELECTOR SWITCHES FOR ADAS PLRP CALCULATIONS      *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
C
      IERR = 0
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RDG1'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     Feature Only useful in OUT
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
c
c      write(0,'(a,8i5)')
c     >  'RDG1:',len(adasid),len(adasex),adasyr,isele,iselr,iselx
C
      MESAGE = 'EXPECTING 2 CHAR, 1 INT, 1 CHAR  AND 4 INTEGERS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,ADASID,ADASYR,ADASEX,
     >                                  ISELE,ISELR,ISELX,ISELD
c
c      write(0,'(a,8i5)')
c     >  'RDG1:',len(adasid),len(adasex),adasyr,isele,iselr,iselx
c
c      write(0,'(3a)')
c     >  'RDG1:',buffer,':'
c      write(0,'(3a)')
c     >  'RDG1:',graph,':'
c      write(0,'(3a)')
c     >  'RDG1:',adasid,':'
c      write(0,'(3a)')
c     >  'RDG1:',adasex,':'
c

      RETURN
C
 9998 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
C
 9999 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RDG1: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END

C
C
C
      SUBROUTINE RD_lp_los (GRAPH,lp_robs,lp_zobs,lp_theta,lp_dtheta,
     >                      lp_instrument_width,lp_bin_width,ierr)
      implicit none
      INTEGER   IERR
      real lp_robs,lp_zobs,lp_theta,lp_dtheta,lp_instrument_width,
     >     lp_bin_width
      CHARACTER GRAPH*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  RD_LP_LOS : LOS DEFINITION FOR LINE PROFILE CALCULATION          *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE   "READER"
      include 'reader'
      CHARACTER MESAGE*72
C
      IERR = 0
      MESAGE = 'END OF FILE ON UNIT 5'
  100 IF (IBUF.EQ.0) READ (5,'(a512)',ERR=9998,END=9998) BUFFER
      WRITE (9,'(1X,A72,1X,A6)') BUFFER,'RD_LP'
      IF (BUFFER(1:1).EQ.'$') GOTO 100
c
c     jdemod - Added so that global plot modifiers could be read from
c              anywhere. 
c
      IF (BUFFER(2:2).EQ.'#') THEN
        CALL Read_AdditionalPlotData(BUFFER)
        GOTO 100
      ENDIF
C
      MESAGE = 'EXPECTING 1 CHAR, 6 REALS'
      READ (BUFFER,*,ERR=9999,END=9999) GRAPH,lp_robs,lp_zobs,lp_theta,
     >                                  lp_dtheta,lp_instrument_width,
     >                                  lp_bin_width
c
      RETURN
C
 9998 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))')
     >  'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
C
 9999 IERR = 1
      WRITE (6,'(1X,A,4(/1X,A))')
     >  'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      WRITE (7,'(1X,A,4(/1X,A))')
     >  'RD_LP: ERROR READING ',GRAPH,MESAGE,'LAST LINE READ :-',BUFFER
      RETURN
      END
C
C
C
      SUBROUTINE ADASRD(YEAR,IZ0,IZ1,ICLASS,NPTS,TE,NE,COEF)     
c      SUBROUTINE ADASRD(YEAR,YEARDF,IZ0,IZ1,ICLASS,NPTS,TE,NE,COEF)
C
C  READ THE REQUESTED RATE COEFFICIENT FROM THE ADAS MASTER ELEMENT
C  FILES:
C        ICLASS = 1: RECOMBINATION RATE COEFFICIENT
C                 2: IONISATION RATE COEFFICIENT
C                 3: CHARGE EXCHANGE RECOMBINATION COEFFICIENT
C                 4: POWER COEF. FOR RECOMBINATION AND BREMSSTRAHLU
C                 5: POWER COEFFICIENT FOR LINE RADIATION
C                 6: POWER COEFFICIENT FOR CHARGE EXCHANGE
C                 added by Krieger, IPP 5/95
C                 7: PHOTON EMISSIVITY FOR DIAGNOSTIC LINES
C  THIS ROUTINE USES THE STANDARD ADAS EXTRACTION ROUTINE D2DATA AND
C  REALLY ONLY PROVIDES A 'CLEAN' INTERFACE, TAKING CARE OF CHANGES
C  IN UNITS AND IN PRECISION OF VARIABLES.  IF THE REQUESTED DATA
C  DOESN'T EXIST (IFAIL=1 RETURNED FROM D2DATA) THE PROGRAM IS STOPPED.
C
      IMPLICIT NONE
C     INCLUDE   "PARAMS"
      include    'params'
C     INCLUDE   "CADAS2"
      include    'cadas2'
C
      CHARACTER*2 YEAR
      character*80 class
      INTEGER IZ0, IZ1, ICLASS, NPTS,len,lenstr
      external lenstr 
      REAL TE(NPTS), NE(NPTS), COEF(NPTS)
c
c      logical lintrp(maxpts)
C
      INTEGER I, J
C
C
c     additional diagnostic output; Krieger IPP/97
c
      class = ' '
      if (iclass.eq.1) then
        class = 'RECOMBINATION RATE COEFFICIENT'
      else if (iclass.eq.2) then
        class = 'IONISATION RATE COEFFICIENT'
      else if (iclass.eq.3) then
        class = 'CHARGE EXCHANGE RECOMBINATION COEFFICIENT'
      else if (iclass.eq.4) then
        class = 'POWER COEF. FOR RECOMBINATION AND BREMSSTRAHLUNG'
      else if (iclass.eq.5) then
        class = 'POWER COEFFICIENT FOR LINE RADIATION'
      else if (iclass.eq.6) then
        class = 'POWER COEFFICIENT FOR CHARGE EXCHANGE'
      else if (iclass.eq.7) then
        class = 'PHOTON EMISSIVITY FOR DIAGNOSTIC LINES'
      endif
c
      IEVCUT = 0
C
      DO J = 1, NPTS
        DTEV(J) = DBLE(ALOG10(TE(J)))
        DDENS(J) = DBLE(ALOG10(NE(J)*1.0E-6))
      ENDDO

      CALL D2DATA(YEAR, YEAR, TITLF, IFAIL,
     >            IZ0, IZ1, ICLASS, NPTS, IEVCUT,
     >            MAXADS, ITMAXD, IDMAXD, IZMAXD,
     >            DTEV, DDENS,
     >            DTEVD, DDENSD, DRCOFD, ZDATA,
     >            DRCOFI
c     >            , LINTRP
     >            )
c
      IF (IFAIL.EQ.1) THEN
        len = lenstr(class) 
        WRITE(6,1000) IZ0, IZ1,YEAR
        write(6,1001) iclass, class(1:len)
        WRITE(7,1000) IZ0, IZ1,YEAR
        write(7,1001) iclass, class(1:len)
        WRITE(0,1000) IZ0, IZ1,YEAR
        write(0,1001) iclass, class(1:len)
        STOP
      ENDIF
C
C  EXTRAPOLATED VALUES ARE RETURNED AS ZERO!
C
      DO J = 1, NPTS
        IF (DRCOFI(J).NE.0.0) THEN
          COEF(J) = 10.**SNGL(DRCOFI(J)) * 1.0E-6
        ELSE
          COEF(J) = 0.0
c
          len = lenstr(class)
          write (6,'(a,4i4,3(1x,g12.5))') 
     >           'Extrapolated coefficient:'//class(1:len),
     >            j,iz0,iz1,npts,te(j),ne(j)
c
        ENDIF

      ENDDO
C
 1000 FORMAT(' ERROR READING REQUESTED ATOMIC DATA!',/,
     >       ' MASTER ELEMENT FILE FOR NUCLEAR CHARGE ',I2,
     >       ' AND ION CHARGE ',I2,
     >       ' WAS NOT FOUND IN YEAR ',A2)
 1001 format(' CLASS ',i2,1x,a)
C
      RETURN
      END
c
c
c
      subroutine find_wall_intersection(ra,za,rb,zb,rint,zint,
     >                        reflection_angle,intersect_normal,
     >                        reflection_option,
     >                        intersect_index,intersect_result,
     >                        intersect_logical)
      implicit none
c
      real ra,za,rb,zb,rint,zint,reflection_angle,intersect_normal
      integer reflection_option,intersect_result,intersect_index
      logical intersect_logical
c
      include 'params'
      include 'walls_com'
c
c     Centralize the processing of the calculation of wall
c     intersections. This code is only called when a the 
c     particle trajectory is calculated to be outside the 
c     wall definition. This code then loops through the 
c     wall segments looking for the one struck - it tries
c     the segment closest to the last known particle position
c     first. If an explicit intersection is not found the 
c     code looks at the distance to the closest intersection - if
c     the distance is small it returns that value.  
c
      REAL TDUM(MAXPTS),XDUM(MAXPTS),YDUM(MAXPTS),WORK(4*MAXPTS)
      INTEGER INDWORK(2,MAXPTS)
      real resulta,resultb
c
      real best,dsq,dr,dz,atan2c
      external atan2c
      integer id,start_index,min_index,ind,sect_index
c
      real*8 rad,zad,rbd,zbd,rintd,zintd,ref_angd,min_dist
      real*8 ravd,zavd,new_dist
      real*8 min_dist_cutoff,min_rintd,min_zintd
      real*8 theta_normal,theta_impact,tnew
c
c     This parameter specifies the maximum distance that an 
c     intersection can be from the given location of the particle. 
c     This needs to be revised upward from 1e-6m since there are
c     too many false rejections under these conditions. 
c
      parameter(min_dist_cutoff=1.0d-4)
c
c     
c
      integer loop_cnt
c     IPP/08 Krieger - two variables for fixed iteration scheme below
      real rtest,ztest,rstep,zstep
c
c     Sometimes, for numerical reasons, the intersection found may be calculated
c     to be outside the boundary - if this is the case then code shifts the point
c     slightly in small steps until the calculation shows it within the boundary. 
c
      real,parameter :: step_dist = 0.000001
      integer,parameter :: max_loop_cnt = 100
c
c     Copy input to double precision 
c
      rad=ra 
      zad=za 
      rbd=rb 
      zbd=zb 
      rintd=rint 
      zintd=zint 
      ref_angd=reflection_angle
      min_dist=HI
      min_index=0
      sect_index=0
c
c     For initial debugging - verify that the point is outside
c     the defined wall. 
c
      CALL GA15A(PCNT,1,WORK,4*MAXPTS,INDWORK,MAXPTS,
     >             RW,ZW,TDUM,XDUM,YDUM,6)
      CALL GA15B(Ra,Za,RESULTa,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
      CALL GA15B(Rb,Zb,RESULTb,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)

      if (resulta*resultb.gt.0.0) then 
c
         write(0,'(a,4g18.10)') 'FIND_WALL_INTERSECTION:'//
     >                     ' ERROR: NO WALL INTERSECTION'//
     >                     ' POSSIBLE: RESULT CHECKS: ',
     >                     resulta,resultb
         write(0,'(a,4g18.10)') 'FIND_WALL_INTERSECTION:'//
     >                     ' NEAREST POINT ON WALL IS RETURNED'
c
         write(6,'(a,4g18.10)') 'FIND_WALL_INTERSECTION:'//
     >                     ' ERROR: NO WALL INTERSECTION'//
     >                     ' POSSIBLE: RESULT CHECKS (A): ',
     >                     resulta,ra,za
         write(6,'(a,4g18.10)') 'FIND_WALL_INTERSECTION:'//
     >                     ' ERROR: NO WALL INTERSECTION'//
     >                     ' POSSIBLE: RESULT CHECKS (B): ',
     >                      resultb,rb,zb
         write(6,'(a,4g18.10)') 'FIND_WALL_INTERSECTION:'//
     >                     ' NEAREST POINT ON WALL IS RETURNED'
c
c
c        Since both points are either inside or outside the boundary - 
c        no intersection on the line segment is possible. In fact, this 
c        code should not have been called once this condition has been 
c        reached since it should really be called as the particle crosses
c        the boundary - for this case - use the extrapolated closest wall
c        intersection to the RA, ZA point.
c

      endif        
c
      if (resulta.gt.0.0) then 
         write(0,'(a)') 'FIND_WALL_INTERSECTION: ERROR:'//
     >     ' WALL INTERSECTION CODE CALLED FOR POINT INSIDE WALL' 
         write(6,'(a)') 'FIND_WALL_INTERSECTION: ERROR:'//
     >     ' WALL INTERSECTION CODE CALLED FOR POINT INSIDE WALL' 
      endif     
c
c     Loop through and find the wall element closest to (ra,za)
c
      BEST = HI
      DSQ  = HI
      IND = 1
      DO 650 ID = 1,WALLPTS
         DSQ = (WALLPT(ID,1)-RA) ** 2 + (WALLPT(ID,2)-ZA) ** 2
         IF (DSQ.LT.BEST) THEN
           BEST = DSQ
           start_INDex   = ID
         ENDIF
650   CONTINUE
C
c         WRITE(6,'(a,2i5,10(1x,g14.8))') 
c     >              'DSQ1:',ind,WALLPTS,DSQ,Ra,Za,Rb,Zb
c
c         WRITE(6,'(a,i5,10(1x,g14.8))') 
c     >    'DSQ2:',ind,WALLPT(ind,20),wallpt(ind,21),wallpt(ind,1),
c     >                wallpt(ind,2),wallpt(ind,22),wallpt(ind,23)
C
c
c      Check nearest wall element first
c      Intersect=1 - intersection point found that lies on both lines
c                2 - intersection point found that lies on trajectory
c                3 - intersection point found that lies on wall
c
      call intsect2dp(rad,zad,rbd,zbd,
     >        dble(wallpt(start_index,20)),dble(wallpt(start_index,21)),
     >        dble(wallpt(start_index,22)),dble(wallpt(start_index,23)),
     >        rintd,zintd,intersect_result)
c
c     A value of 1 indicates that the appropriate intersection point has been found
c
      if (intersect_result.eq.3) then  
         min_dist=(rad-rintd)**2+(zad-zintd)**2
         min_index = start_index
         min_rintd=rintd
         min_zintd=zintd
      elseif(intersect_result.eq.1) then 
         sect_index=start_index  
      endif
c
      if (intersect_result.ne.1) then     
c
c        Loop around the rest of the wall looking for an intersection point
c
         id = 0
c
         do while (intersect_result.ne.1.and.id.le.wallpts/2)
c
            id = id +1
c
            ind = start_index+id
            if (ind.gt.wallpts) ind = ind-wallpts
c
            call intsect2dp(rad,zad,rbd,zbd,
     >          dble(wallpt(ind,20)),dble(wallpt(ind,21)),
     >          dble(wallpt(ind,22)),dble(wallpt(ind,23)),
     >          rintd,zintd,intersect_result)
c
c           Exit if intersection is 1
c

            if (intersect_result.eq.1) then 
               sect_index=ind
               exit
               !cycle
            endif
c 
c           Record distance if intersection is 3
c
            if (intersect_result.eq.3) then 
               new_dist=(rad-rintd)**2+(zad-zintd)**2
c
c              check if distance is less than current minimum and record
c              it if it is.
c
               if (new_dist.lt.min_dist) then 
                  min_dist=new_dist 
                  min_index=ind
                  min_rintd=rintd
                  min_zintd=zintd
               endif
            endif
c
c           Check index in the other direction
c
            ind = start_index - id
            if (ind.lt.1) ind = wallpts+ind 
c           
            call intsect2dp(rad,zad,rbd,zbd,
     >          dble(wallpt(ind,20)),dble(wallpt(ind,21)),
     >          dble(wallpt(ind,22)),dble(wallpt(ind,23)),
     >          rintd,zintd,intersect_result)
c
c           Exit if intersection is 1
c
            if (intersect_result.eq.1) then 
               sect_index=ind
               exit
               !cycle
            endif
c 
c           Record distance if intersection is 3
c
            if (intersect_result.eq.3) then 
               new_dist=(rad-rintd)**2+(zad-zintd)**2              
c
c              check if distance is less than current minimum and record
c              it if it is.
c
               if (new_dist.lt.min_dist) then 
                  min_dist=new_dist 
                  min_index=ind
                  min_rintd=rintd
                  min_zintd=zintd
               endif
            endif
      
         end do
      endif
c
c     Code has not found a proper intersection
c
      if (intersect_result.ne.1) then     
c
c        Check for a close intersection  
c        Min_index is non-zero if at least some intersections have been found
c
         if (min_index.ne.0) then
c
            min_dist=sqrt(min_dist)
c
c           If the distance is less than the minimum or if there was an error in input 
c           and the test line segment does not cross the wall - then use the intersection
c           closest to the test point. 
c
            if ((min_dist.lt.min_dist_cutoff).or.
     >          (resulta*resultb.gt.0.0)) then 
               intersect_result = 1
               sect_index = min_index  
               rintd=min_rintd
               zintd=min_zintd
            endif
c
         endif 
c
      endif
      
c
c     If intersection is 1 at this point then an acceptable intersection has been found
c     - otherwise an error condition has been encountered in which no intersection
c       point could be found
c
c
c     For intersect=1 - find the reflection angle
c
      if (intersect_result.eq.1) then 
c
c        Angle of particle trajectory 
c
         dz=za-zb
         dr=ra-rb
c
         theta_impact = atan2c(dz,dr)
c
c        Surface normal angle
c 
         dz=wallpt(sect_index,23)-wallpt(sect_index,21)
         dr=wallpt(sect_index,22)-wallpt(sect_index,20)
c
         theta_normal = atan2c(dz,dr) - PI/2.0      
c
         CALL REFANGDP(Theta_NORMal,Theta_IMPact,TNEW,reflection_option)
c
         write(6,'(a,2i10,6g18.10)') 'FIND_WALL_INTERSECTION: REFANG:',
     >        reflection_option,sect_index,
     >        theta_normal*raddeg,theta_impact*raddeg,tnew*raddeg
c
c        Copy outputs to input variables
c
c        Assign reflection angle and intersection point - and intersection wall index
c
         reflection_angle = tnew
         rint=rintd
         zint=zintd
         intersect_index = sect_index
         intersect_logical = .true.
         intersect_normal = theta_normal
c
c     ERROR condition
c
      else
c
c        Return the center point of the wall segment closest to the initial
c        position of the particle trajectory.
c
         rint=wallpt(start_index,1)
         zint=wallpt(start_index,2)
c
c        Set to normal to this wall segment
c
         reflection_angle = wallpt(start_index,9) - PI/2.0
c
c        Set intersection logical to false and intersection_index to start_index
c
         intersect_index = start_index
         intersect_logical = .false.
         intersect_normal = reflection_angle
c
c         write(0,'(a)') 'FIND_WALL_INTERSECTION: ERROR:'//
c     >                   ' NO INTERSECTION POINT FOUND'
c
         write(6,'(a)') 'FIND_WALL_INTERSECTION: ERROR:'//
     >                   ' NO INTERSECTION POINT FOUND'
        write(6,'(a,7g18.10,i6,4g18.10)') 'DATA:',ra,za,rb,zb,rint,zint,
     >                 reflection_angle*raddeg,start_index,
     >          wallpt(start_index,20),wallpt(start_index,21),
     >          wallpt(start_index,22),wallpt(start_index,23)

        write(6,'(a,2i6,5g18.10)')
     >       'FIND_WALL_INTERSECTION: INTSECT 3:',min_index,
     >        intersect_result,min_dist,
     >        min_rintd,min_zintd
c
      endif 
c
c
c     Check to see if the intersection point calculated is inside or outside the wall. 
c
      CALL GA15B(Rint,Zint,RESULTa,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c
      if (resulta.lt.0.0) then 
c         write(0,'(a,2g18.9)') 
c     >     'FIND_WALL_INTERSECTION: INTERSECTION IS OUTSIDE WALL:',
c     >      rint,zint
         write(6,'(a,2g18.9)') 
     >     'FIND_WALL_INTERSECTION: INTERSECTION IS OUTSIDE WALL:',
     >      rint,zint
c
c        Take corrective action to ensure point is inside wall
c
         loop_cnt = 0.0
c        IPP/08 Krieger - make this numerically more stable by
c        replacing recursive iteration with explicit iteration
         
         rstep = step_dist * cos(reflection_angle)
         zstep = step_dist * sin(reflection_angle)

         do while (resulta.lt.0.0.and.loop_cnt.lt.max_loop_cnt) 

            loop_cnt = loop_cnt + 1.0

            rtest = rint + loop_cnt * rstep
            ztest = zint + loop_cnt * zstep

            CALL GA15B(Rtest,Ztest,RESULTa,PCNT,1,WORK,4*MAXPTS,
     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)

         end do
c
c         rtest = rint
c         ztest = zint
c
c         do while (resulta.lt.0.0.and.loop_cnt.lt.max_loop_cnt) 
c
c            loop_cnt = loop_cnt + 1.0
c
c            rtest = rtest + step_dist * cos(reflection_angle)
c            ztest = ztest + step_dist * sin(reflection_angle)
c
c            CALL GA15B(Rtest,Ztest,RESULTa,PCNT,1,WORK,4*MAXPTS,
c     >             INDWORK,MAXPTS,RW,ZW,TDUM,XDUM,YDUM,6)
c
c         end do
c
c
c        Revised intersection point found within search parameters
c
         if (resulta.gt.0.0) then 
            rint = rtest
            zint = ztest
            
            write(6,'(a,3g18.10,i10)')
     >          'FIND_WALL_INTERSECTION: REVISED INTERSECTION FOUND  :',          
     >          rint,zint,reflection_angle*raddeg,loop_cnt

         else
            write(6,'(a,3g18.10,i10)')
     >          'ERROR: FIND_WALL_INTERSECTION: '//
     >          'REVISED INTERSECTION NOT  FOUND  :',          
     >          rint,zint,reflection_angle*raddeg,loop_cnt

            




         end if




      endif



      return
      end
c
c
c
      subroutine dist_to_segment(x1,y1,x2,y2,xp,yp,xi,yi,dist,sect,opt)
      implicit none
      real x1,y1,x2,y2,xp,yp,xi,yi,dist
      integer,intent(in) :: opt
      logical sect
c      
c     DIST_TO_SEGMENT: This routine returns the distanced from a point to line 
c                      segment along with the normal intersection point. The
c                      variable sect is set to true if this point lies
c                      between the ends of the line segment
c
c     Opt is a switch used to control the function of the routine - it controls the results obtained
c     when the normal point does not lie on the line segment. 
c
c     Opt = 0 : normal point is returned whether or not it lies on the line segment. 
c     Opt = 1 : normal point is returned as the nearest end point if point is not on line segment
c               distance is the distance to the endpoint
c     Opt = 2 : normal point is the center of the line segment 
c               dist is the distance to the center of the center of the line segment 
c     
      real dx,dy,u
c
      sect = .false.
c
      dx = x2-x1
      dy = y2-y1
c
c
c     Line segment is actually a point - error condition - return point and distance to point 
c
      if (dx.eq.0.0.and.dy.eq.0.0) then
         xi = x1
         yi = y1

         sect = .false.
c
c     Solve for point of tangency/normal intersection to line segment from xp,yp
c
      else
c
         u = ((xp-x1) * dx + (yp-y1) * dy) / (dx**2 + dy**2)
c
         if (u.ge.0.0.and.u.le.1.0) then 
            sect = .true.
            xi = x1 + u * dx
            yi = y1 + u * dy
c
         else
            sect = .false.
            if (opt.eq.0) then 
               ! Return normal point outside line segment
               xi = x1 + u * dx
               yi = y1 + u * dy
            elseif (opt.eq.1) then 
               ! return end points of line segment 
               if (u.lt.0.0) then 
                  xi = x1
                  y1 = y1
               else
                  xi = x2
                  y1 = y2
               endif
            elseif(opt.eq.2) then 
               ! return center point of line segment
               xi = (x2+x1)/2.0
               yi = (y2+y1)/2.0
            endif
         endif
c
      endif
c
c     Calculate distance
c
      dist = sqrt((xp-xi)**2+(yp-yi)**2)
c
      return
      end
c
c
c     
      subroutine find_nearest_point_on_wall(rsect,zsect,id_out,is_out)
      implicit none
c
      real rsect,zsect
      integer id_out,is_out
c
      include 'params'
      include 'walls_com'
c
c     FIND_NEAREST_POINT_ON_WALL: This routine scans through the 
c         segments forming the wall to find a point on the wall 
c         segment that is perpendicularly closest to the given test point.
c         The closest point must lie on the line segment forming the section
c         of wall to which it is closest. 
c
c         The test point data is then replaced with the calculated wall point.
c
c     In addition, the code also returns the index of the walls segment where the 
c     point is found and whether the point is in the first half or second half of 
c     that wall segment.
c
c     Note - this code is based on the premise that the test point is generally quite
c     close to the wall to start with and all that is required is a refinement of the 
c     position. If this is not the case, it is possible for this code to give some
c     odd results. The code will try to intelligently choose between the closest wall
c     segment and the nearest perpendicular walls segment if there is a large discrepancy in
c     the distances involved. 
c
c
      integer id,min_dist_id,min_perp_dist_id
      real min_dist,min_perp_dist
      real r_md,z_md,r_mpd,z_mpd
      real rtest,ztest,dist
      real d1,d2,r1,z1,r2,z2
      logical sect
c
c     Loop over wall
c
      min_dist_id = 0
      min_dist = HI
      min_perp_dist_id = 0
      min_perp_dist = HI
c
c     Loop over wall      
c
      do id = 1,wallpts
c
c        Coordinates of endpoints of wall segment 
c
         r1 = wallpt(id,20)
         z1 = wallpt(id,21)
         r2 = wallpt(id,22)
         z2 = wallpt(id,23)
c
         call dist_to_segment(r1,z1,r2,z2,rsect,zsect,rtest,ztest,dist,
     >                        sect,2)
c
c
         if (sect) then 
            if (dist.lt.min_perp_dist) then 
               min_perp_dist = dist
               min_perp_dist_id = id
               r_mpd = rtest
               z_mpd = ztest
            endif
         else
            if (dist.lt.min_dist) then 
               min_dist = dist
               min_dist_id = id
               r_md = rtest
               z_md = ztest
            endif
         endif
c
      end do
c
c     Perpendicular wall segment was found 
c
      if (min_perp_dist_id.ne.0) then 
c
c        Perpendicular distance is much larger than the 
c        center of the closest wall segment on the wall - use
c        closest point instead.
c
         if (min_perp_dist.gt.10.0* min_dist) then 
c
c           Since this is the center point of the nearest wall segment - IS can be set to 0 or 1
c
            rsect = r_md
            zsect = z_md
            id_out = min_dist_id
            is_out = 0 
         else
            rsect = r_mpd
            zsect = z_mpd
            id_out = min_perp_dist_id
c
c           Determine is_out
c            
            d1 = (rsect-wallpt(id,20))**2 + (zsect-wallpt(id,21))**2 
            d2 = (rsect-wallpt(id,22))**2 + (zsect-wallpt(id,23))**2 
c
c           First half of wall segment is 0 - second half is 1 - corresponds to the values in ISPRODS 
c
            if (d1.lt.d2) then 
               is_out = 0
            else
               is_out = 1
            endif

         endif
c
c     No perpendicular wall segment found - use closest
c
      else
c
c        Since this is the center point of the nearest wall segment - IS can be set to 0 or 1
c
         rsect = r_md
         zsect = z_md
         id_out = min_dist_id
         is_out = 0 
c         
      endif
c
      return
      end
c
c
c
      subroutine find_free_unit_number(unit)
      implicit none
      integer unit
c
c     FIND_FREE_UNIT_NUMBER:
c
c     This routine scans through unit numbers looking for one that
c     is not currently in use. This number is returned. This code
c     is based on the assumption that any unit numbers returned will
c     be used before this routine is called again asking for another 
c     number - otherwise it will likely return the previous value.
c
      integer test_unit
      logical unit_open

      test_unit = 10
      unit_open = .false.

      ! Check for unit number assignment.  
      Do While (Unit_open)
         test_unit=test_unit + 1
         Inquire (Unit = test_unit, Opened = Unit_open)
      End Do

      unit = test_unit

      return
      end

c
c
c
      subroutine calc_wall_length_coordinate(opt)
      implicit none
      integer opt
c
c     Use the data in the wallpt array to calculate the 
c     distance along the walls from the Inside mid-plane
c     counter clockwise. Since the simulation walls
c     do not match the actual geometry exactly - this 
c     coordinate will only be approximate and a 
c     shift may need to be applied to match any 
c     experimental system. This shift could be
c     calculated by determining the corresponding 
c     coordinate values at one common point and matching
c     the values there. 
c
c     OPT is available to allow for different calculation schemes later
c
c
      include 'params'
      include 'walls_com'
c
      integer nw,cnt,cin,in,startin
      real minr,totdist,rminw
      real r0,z0,dr,dz

c
c     Find the element of wall stradling the inside mid-plane. 
c      
      minr = hi
c
      do in = 1,wallpts
         if ((wallpt(in,21)*wallpt(in,23)).le.0.0) then 
            rminw = min(wallpt(in,20),wallpt(in,22))
            if (rminw.lt.minr) then 
               startin = in
               minr = rminw
            endif
         endif
      enddo
c
c     Now have the starting index - need to go counter clockwise
c     around the wall from this location recording the 
c     coordinate for the center of the wall element. 
c     Note: the center of the first wall element could either be
c           at the beginning or the end of the wall. 
c
      cnt = 0

      z0 = 0.0
      in = startin
      dr = wallpt(in,20)-wallpt(in,22)
      dz = wallpt(in,21)-wallpt(in,23)

      r0 = wallpt(in,22) + (z0-wallpt(in,23))/dz * dr
c
c     since working counter clockwise - need to move around the wall
c     backwards
c     Set initial distance to edge of first wall element 
c
      totdist = sqrt((wallpt(in,20)-r0)**2 
     >             + (wallpt(in,21)-z0)**2)
c
c     Set wall element counter
c
      cnt = 1

      do
         cin = startin - cnt
         if (cin.lt.1) cin = cin + wallpts
c
c        Check if finished wall and complete last wall element 
c
         if (cin.eq.startin) then 
            if (wallpt(cin,2).lt.0.0) then 
               wallpt(cin,32) = sqrt((wallpt(cin,1)-r0)**2 
     >                             + (wallpt(cin,2)-z0)**2)
            else
               wallpt(cin,32)= totdist +
     >                     sqrt((wallpt(cin,22)-wallpt(cin,1))**2 
     >                         +(wallpt(cin,23)-wallpt(cin,2))**2)
            endif
            exit
         endif
c
c        Continue if not the last element 
c

         totdist = totdist+ sqrt((wallpt(cin,22)-wallpt(cin,1))**2 
     >           + (wallpt(cin,23)-wallpt(cin,2))**2)

         wallpt(cin,32) = totdist
         
         totdist = totdist+ sqrt((wallpt(cin,20)-wallpt(cin,1))**2 
     >           + (wallpt(cin,21)-wallpt(cin,2))**2)


         cnt = cnt + 1

      end do

      write(6,*) 'DISTANCE ALONG WALL',startin,wallpts
      do in = 1,wallpts
         write(6,'(3i6,f12.5)') in,int(wallpt(in,16)),
     >                  int(wallpt(in,18)),wallpt(in,32)
      end do


      return 
      end
c
c
c
      subroutine get_plasma_rz(r,z,ne,te,ti,vb,ef,nh,nh_mol)
      implicit none
      real r,z,ne,te,ti,vb,ef,nh,nh_mol

      include 'params'
      include 'cgeom'
      include 'pindata'

      !
      ! GET_PLASMA_RZ: This routine returns the plasma data found on the 
      !                grid at location R,Z. If this location is off the 
      !                grid the first default action is to return a density
      !                of 1.0e10 and a temperature of 1.0eV 
      !                (essentially no plasma) - this could be modified in
      !                future to return some estimate of local off grid
      !                conditions. 
      !
      
      !
      ! Local variables
      !
      integer :: lastik = 0
      integer :: lastir = 0
      save lastik,lastir
      integer :: ik,ir
      logical :: outofgrid,newinj

      ik = lastik
      ir = lastir
      newinj = .true.
      outofgrid = .false.
      
      call gridpos(ik,ir,r,z,newinj,outofgrid)

      if (outofgrid) then
         ne = 1.0e14
         te = 1.0
         ti = 1.0
         vb = 0.0
         ef = 0.0
         nh = 0.0
         nh_mol = 0.0
         lastik = 0
         lastir = 0
      else
         ne = knbs(ik,ir)
         te = ktebs(ik,ir)
         ti = ktibs(ik,ir)
         vb = kvhs(ik,ir)
         ef = kes(ik,ir)
         nh = pinatom(ik,ir)
         nh_mol = pinmol(ik,ir)
         lastik = ik
         lastir = ir
      endif

      return
      end
