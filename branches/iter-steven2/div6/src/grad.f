c
c ======================================================================
c
c function: PolySideLen
c
c Return the length of a specified side of polygon IK,IR.
c
c ======================================================================
c
      real function PolySideLen(ik,ir,side,rc)
     
      implicit none
c
c     Side should be treated as a CONSTANT and it's value left 
c     unchanged. RC=return code - normally 0
c
      integer ik,ir,side,rc
c 
c     Global Variables:
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     PolySideLen - returns the length of the specified side of the 
c                  ik,ir polygon - if the a polygon or the specified
c                  side does not exist the function returns a value
c                  of 1.0.
c
c     side is defined to be one of four values that define the side
c         of interest in the specific geometry. 
c
c     side = INWARD41 - side facing core or into PP region.
c                       Polygon corners 1,4
c          = OUTWARD23- side facing outer wall or towards separatrix
c                       from inside the PP region. Polygon corners 2,3
c          = DOWN34   - side facing along field lines toward IK=NKS(IR)
c                       target. Polygon corners 4,3
c          = UP12     - side facing along the field lines toward the
c                       IK=1 target.Polygon corners 1,2
c
c     Local Variables:
c
      integer in,nside
      real    deltar, deltaz
c
      rc = 0
      in = korpg(ik,ir)
c
c     Check for polygon 
c 
      if (in.le.0) then
         write (6,*) 'PolySideLen: ERROR : No polygon'//
     >              ' exists for (ik,ir) ',ik,ir
         rc = 1 
         polysidelen = kpsiz(ik,ir)
         return
      endif 
c
      nside = nvertp(in)
c
c     Check number of sides
c 
      if (nside.lt.side) then
         write (6,*) 'PolySideLen: ERROR : Specified side'//
     >             ' does not exist for (ik,ir) ',ik,ir,side,nside
         polysidelen = 1.0
         return
      endif 
c
c     Calculate the length of the polygon side.
c      
      if (side.eq.INWARD41) then 

         deltar = rvertp(4,in)-rvertp(1,in)
         deltaz = zvertp(4,in)-zvertp(1,in)
 
      else

         deltar = rvertp(side,in)-rvertp(side+1,in)
         deltaz = zvertp(side,in)-zvertp(side+1,in)
         
      endif
c
      polysidelen = sqrt(deltar**2+deltaz**2)
c
      return
      end
c
c ======================================================================
c
c function: CellWidth
c
c Return an estimate of the width of cell IK,IR.
c
c ======================================================================
c
      REAL FUNCTION CellWidth(ik,ir,s,mode)

      IMPLICIT none
c
c Input:
c      
      INTEGER ik,ir,mode
      REAL    s

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL    r,z,r1,z1,r2,z2,c1,c2,wid(3)
      REAL    z3,z4,r3,r4
      REAL    widtop,widbot,cbot1,ctop1,ctop2,cbot2
      REAL    deltar,deltaz
      REAL    rside1,zside1,rside2,zside2
      INTEGER id
c
c Hard coded width option (need to restore this when adding to
c Dave's version, since slcom is not included there):
c
c      INTEGER cwcopt

      cwcopt = 3

      id = korpg(ik,ir)
c
c 
c
      IF (mode.NE.SIDE14.AND.mode.NE.side23.AND.mode.NE.TOTAL)
     .  STOP 'ERROR (CellWidth): Invalid MODE value'
c
c
c
      IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap) THEN 
        WRITE(0,*) 'ERROR (CellWidth): Virtual ring (IR = ',ir,')'
        STOP 
      ENDIF
c
c 
c
      IF (nvertp(id).NE.4) THEN
        write (0,*) 'ERROR (CellWidth): Invalid cell',ik,ir,id,
     >   nvertp(id)
        STOP 'ERROR (CellWidth): Invalid cell'
      ENDIF
c
c
c
      IF (s.EQ.CENTER) THEN
        r = rs(ik,ir)
        z = zs(ik,ir)
      ELSE
        IF (cwcopt.NE.3) THEN
          WRITE(0,*) 'ERROR (CellWidth): Non-central S values are ',
     .               'invalid (CWCOPT = ',cwcopt,')'
          STOP
        ENDIF
c change ID if changing IK
c assign IK to a non-refenced variable if altering its value
        STOP 'ERROR (CellWidth): Invalid S value'
      ENDIF
c
c
c
      IF (cwcopt.EQ.1) THEN
c
c Approximate the cell width as the sum of the distances between the
c cell side mid-points and the cell center:
c
        r1 = 0.5 * (rvertp(1,id) + rvertp(4,id))
        z1 = 0.5 * (zvertp(1,id) + zvertp(4,id))

        r2 = 0.5 * (rvertp(2,id) + rvertp(3,id))
        z2 = 0.5 * (zvertp(2,id) + zvertp(3,id))

        wid(SIDE14) = SQRT((r - r1)**2 + (z - z1)**2)
        wid(SIDE23) = SQRT((r - r2)**2 + (z - z2)**2)

c Debug:
c        WRITE(SLOUT,'(10X,A,2I3,F12.6)')
c     +    'WidDat:',IK,IR,CellWidth

      ELSEIF (cwcopt.EQ.2) THEN
c
c More 'advanced' method - DO NOT USE:
c
        rside1 = 0.5 * (rvertp(4,id) + rvertp(1,id))
        zside1 = 0.5 * (zvertp(4,id) + zvertp(1,id))

        deltar = rvertp(4,id) - rside1
        deltaz = zvertp(4,id) - zside1

        ctop1 = ((rvertp(3,id) - rside1) * deltar +
     .           (zvertp(3,id) - zside1) * deltaz) /
     .          (deltar**2 + deltaz**2)

        cbot1 = ((rvertp(2,id) - rside1) * deltar +
     .           (zvertp(2,id) - zside1) * deltaz) /
     .          (deltar**2 + deltaz**2)

        rside2 = 0.5 * (rvertp(3,id) + rvertp(2,id))
        zside2 = 0.5 * (zvertp(3,id) + zvertp(2,id))

        deltar = rvertp(3,id) - rside2
        deltaz = zvertp(3,id) - zside2
	
        ctop2 = ((rvertp(4,id) - rside2) * deltar +
     .           (zvertp(4,id) - zside2) * deltaz) /
     .          (deltar**2 + deltaz**2)
	
        cbot2 = ((rvertp(1,id) - rside2) * deltar +
     .           (zvertp(1,id) - zside2) * deltaz) /
     .          (deltar**2 + deltaz**2)
c
c Decide whether to use vertex 3 or 4 when estimating the width
c of the 'top' of the cell:
c
        IF (abs(ctop1).LT.abs(ctop2)) THEN
          r1 = rvertp(3,id)
          z1 = zvertp(3,id)
          r2 = ctop1 * (rvertp(4,id) - rside1) + rside1
          z2 = ctop1 * (zvertp(4,id) - zside1) + zside1
        ELSE
          r1 = rvertp(4,id)
          z1 = zvertp(4,id)
          r2 = ctop2 * (rvertp(3,id) - rside2) + rside2
          z2 = ctop2 * (zvertp(3,id) - zside2) + zside2
        ENDIF
c
c Decide whether to use vertex 1 or 2 when estimating the width
c of the 'bottom' of the cell:
c
        IF (abs(cbot1).LT.abs(cbot2)) THEN   
          r3 = rvertp(2,id)
          z3 = zvertp(2,id)
          r4 = cbot1 * (rvertp(4,id) - rside1) + rside1
          z4 = cbot1 * (zvertp(4,id) - zside1) + zside1
        ELSE
          r3 = rvertp(1,id)
          z3 = zvertp(1,id)
          r4 = cbot2 * (rvertp(3,id) - rside2) + rside2
          z4 = cbot2 * (zvertp(3,id) - zside2) + zside2
        ENDIF

        wid(SIDE14) = 0.25 * (SQRT((r2 - r1)**2.0 + (z2 - z1)**2.0) +
     +                        SQRT((r4 - r3)**2.0 + (z4 - z3)**2.0))
        wid(SIDE23) = wid(SIDE14)
c Debug:
c        WRITE(SLOUT,'(10X,A,2I3,4F10.4)')
c     +    'WidDat:',IK,IR,ctop1,ctop2,cbot1,cbot2
c        WRITE(SLOUT,'(14X,9F12.6)')
c     +    R1,Z1,R2,Z2,R3,Z3,R4,Z4,CellWidth

      ELSEIF (cwcopt.EQ.3) THEN
c
c Perpendicular distance from cell center to cell sides:
c
c
c Side 1-4:
c
        deltar = rvertp(4,id) - rvertp(1,id)
        deltaz = zvertp(4,id) - zvertp(1,id)

        c1 = ((r - rvertp(1,id)) * deltar +
     .        (z - zvertp(1,id)) * deltaz) /
     .       (deltar**2 + deltaz**2)

        r1 = rvertp(1,id) + c1 * deltar
        z1 = zvertp(1,id) + c1 * deltaz
c
c Side 2-3:
c
        deltar = rvertp(3,id) - rvertp(2,id)
        deltaz = zvertp(3,id) - zvertp(2,id)

        c2 = ((r - rvertp(2,id)) * deltar +
     .        (z - zvertp(2,id)) * deltaz) /
     .       (deltar**2 + deltaz**2)

        r2 = rvertp(2,id) + c2 * deltar
        z2 = zvertp(2,id) + c2 * deltaz
c
c Calculate total distance:
c
        wid(SIDE14) = SQRT((r - r1)**2 + (z - z1)**2)
        wid(SIDE23) = SQRT((r - r2)**2 + (z - z2)**2)

c Debug:
c        WRITE(SLOUT,'(10X,A,2I3,4F10.4)')
c     +    'WidDat:',ik,ir,c1,c2,CellWidth

      ELSEIF (cwcopt.EQ.4) THEN
        wid(SIDE14) = 0.5
        wid(SIDE23) = 0.5
      ELSE
        STOP 'ERROR (CellWidth): Invalid option'
      ENDIF

      wid(TOTAL) = wid(SIDE14) + wid(SIDE23)

      CellWidth = wid(mode)

      RETURN
      END
c 
c = function end: CellWidth ============================================
c 
c
c ======================================================================
c
c function: GetTheta
c
c Returns the THETA value corresponding to a given S value for the 
c ring IR.
c
c Aug 14, 97 - The current method of handling the x-point 
c region in the core will be inaccurate if the method of determining 
c the distance between points on neighbouring rings differs from
c from the current implimentation, which uses an estimate of 
c the cell widths (sorry about the run-on sentence).
c
c ======================================================================
c
      REAL FUNCTION GetTheta(ir,s)
      implicit none 
      REAL    s
      INTEGER ir,ik

      INCLUDE 'params'                                                  
      INCLUDE 'cgeom'                                                  
      INCLUDE 'comtor'
c
c Check for valid IR value (core and virtual rings not allowed):
c 
      IF (ir.LT.irsep.OR.ir.EQ.irtrap.OR.ir.EQ.irwall) THEN
        WRITE(0,*) 'ERROR (GetTheta): Invalid ring (IR = ',ir,')'
        STOP
      ENDIF
c
c Find ik:
c      
      DO ik = 1, nks(ir)
        IF (s.LE.kss(ik,ir)) GOTO 10
      ENDDO
10    CONTINUE
c
c Interpolate linearly to find GetTheta:
c
      IF (ik.EQ.1) THEN
        GetTheta = thetat(idds(ir,2)) +
     .               (s             - ksb   (0,ir)      ) * 
     .               (thetag(ik,ir) - thetat(idds(ir,2))) / 
     .               (kss   (ik,ir) - ksb   (0,ir)      )
      ELSEIF (ik.EQ.nks(ir)+1) THEN
        ik = nks(ir)
        GetTheta = thetag(ik,ir) +
     .               (s                  - kss   (ik,ir)) *
     .               (thetat(idds(ir,1)) - thetag(ik,ir)) *
     .               (ksb   (ik,ir)      - kss   (ik,ir))
      ELSE
        GetTheta = thetag(ik-1,ir) +
     .               (s             - kss   (ik-1,ir)) *
     .               (thetag(ik,ir) - thetag(ik-1,ir)) /
     .               (kss   (ik,ir) - kss   (ik-1,ir))
      ENDIF
c
c Check that GetTheta is not out of bounds:
c
      IF (GetTheta.LT.thetat(idds(ir,2)).OR.
     .    GetTheta.GT.thetat(idds(ir,1))) THEN
        WRITE(0,*) 'ERROR (GetTheta): Calculated THETA value is out ',
     .             'of bounds (IR = ',ir,'  THETA = ',GetTheta,')'   
        GetTheta = 1.0
c        STOP
      ENDIF

      RETURN
      END
c
c = function end: GetTheta =============================================
c
c
c ======================================================================
c
c function: GetS
c
c Returns the S value corresponding to a given THETA value for the 
c ring IR.
c
c some liberties around the x-point
c check that ksb is defined for core rings
c
c ======================================================================
c
      REAL FUNCTION GetS(ir,thetav)
      implicit none
      REAL    thetav
      INTEGER ir,ik

      INCLUDE 'params'                                                  
      INCLUDE 'cgeom'                                                  
      INCLUDE 'comtor'
c
c Check for valid IR value (virtual rings not allowed):
c 
      IF (ir.EQ.1.OR.ir.EQ.irtrap.OR.ir.EQ.irwall) THEN
        WRITE(0,*) 'ERROR (GetS): Virtual ring (IR = ',ir,')'
        STOP
      ENDIF
c
c Check if THETA is out of bounds:
c
      IF (ir.GE.irsep) THEN
        IF (thetav.LT.thetat(idds(ir,2))) THEN
          GetS = ksb(0,ir)
          RETURN
        ELSEIF (thetav.GT.thetat(idds(ir,1))) THEN
          GetS = ksb(nks(ir),ir)
          RETURN
        ENDIF
      ELSE
        IF (thetav.LE.thetag(1,ir)) THEN
          GetS = kss(1,ir)
          RETURN
        ELSEIF (thetav.GE.thetag(nks(ir),ir)) THEN
          GetS = kss(nks(ir),ir)
          RETURN
        ENDIF
      ENDIF
c
c Find IK:
c      
      DO ik = 1, nks(ir)
        IF (thetav.LE.thetag(ik,ir)) GOTO 10
      ENDDO
10    CONTINUE
c
c Interpolate linearly to find S:
c
      IF (ik.EQ.1) THEN
        GetS = ksb(0,ir) +
     .               (thetav        - thetat(idds(ir,2))) * 
     .               (kss   (ik,ir) - ksb   (0,ir)      ) / 
     .               (thetag(ik,ir) - thetat(idds(ir,2)))
      ELSEIF (ik.EQ.nks(ir)+1) THEN
        ik = nks(ir)
        GetS = kss(ik,ir) +
     .               (thetav             - thetag(ik,ir)) *
     .               (ksb   (nks(ir),ir) - kss   (ik,ir)) *
     .               (thetat(idds(ir,1)) - thetag(ik,ir))
      ELSE
        GetS = kss(ik-1,ir) +
     .               (thetav        - thetag(ik-1,ir)) *
     .               (kss   (ik,ir) - kss   (ik-1,ir)) /
     .               (thetag(ik,ir) - thetag(ik-1,ir))
      ENDIF

c Debug:
c	write(0,*) 's stuff: ',ir,ksb(ik-1,ir),GetS,ksb(ik,ir)

c
c Check that GetS is not out of bounds:
c
      IF (GetS.LT.ksb(0,ir).OR.GetS.GT.ksb(nks(ir),ir)) THEN
        WRITE(0,*) 'ERROR (GetS): Calculated S value is out ',
     .             'of bounds (IR = ',ir,'  S = ',GetS,')'   
        STOP
      ENDIF

      RETURN
      END
c
c = function end: GetS =================================================
c
c
c ======================================================================
c
c subroutine: GetPos
c
c Determine P, the poloidal co-ordinate, and IK, the cell index,
c for ring IR given THETA.
c
c ======================================================================
c
      SUBROUTINE GetPos(ir,thetav,p,ik)
      implicit none 
c
c Input:
c
      REAL    thetav
      INTEGER ir
c
c Output:
c
      REAL    p
      INTEGER ik

      INCLUDE 'params'                                                  
      INCLUDE 'cgeom'                                                  
      INCLUDE 'comtor'
c
c Check for valid IR value:
c 
      IF (ir.EQ.1.OR.ir.EQ.irtrap.OR.ir.EQ.irwall) THEN
        WRITE(0,*) 'ERROR (GetPos): Cannot calculate poloidal co-',
     .             'ordinate for virtual rings (IR = ',ir,')'
        STOP
      ENDIF
c
c Check that THETAV is valid for the ring, otherwise return target location:
c
      IF (ir.GE.irsep) THEN
c
c SOL and PFZ:
c
        IF (thetav.LT.thetat(idds(ir,2))) THEN
          ik = 1
          p  = kpb(0,ir)
          RETURN
        ELSEIF (thetav.GT.thetat(idds(ir,1))) THEN
          ik = nks(ir)
          p  = kpb(ik,ir)
          RETURN
        ENDIF
      ELSE
c
c Core:
c
        IF (thetav.LE.thetag(1,ir)) THEN 
          ik = 1
          p  = kps(ik,ir)
          RETURN
        ELSEIF (thetav.GT.thetag(nks(ir)-1,ir)) THEN
c Should this be nks(ir) instead of nks(ir)-1...?
          ik = nks(ir)-1
          p  = kps(ik,ir)
          RETURN
        ENDIF
      ENDIF
c
c Find IK:
c      
      DO ik = 1, nks(ir)
        IF (thetav.LE.thetag(ik,ir)) GOTO 10
      ENDDO
10    CONTINUE
c
c Interpolate linearly to find P:
c
      IF (ik.EQ.1) THEN
        p = kpb(0,ir) + (thetav        - thetat(idds(ir,2))) * 
     .                  (kps   (ik,ir) - kpb   (0,ir)      ) / 
     .                  (thetag(ik,ir) - thetat(idds(ir,2)))
      ELSEIF (ik.EQ.nks(ir)+1) THEN
        ik = nks(ir)
        p  = kps(ik,ir) + (thetav             - thetag(ik,ir)) * 
     .                    (kpb   (ik,ir)      - kps   (ik,ir)) / 
     .                    (thetat(idds(ir,1)) - thetag(ik,ir))
      ELSE
        p = kps(ik-1,ir) + (thetav        - thetag(ik-1,ir)) * 
     .                     (kps   (ik,ir) - kps   (ik-1,ir)) / 
     .                     (thetag(ik,ir) - thetag(ik-1,ir)) 
c
c Adjust IK based on P if necessary:
c
        IF (p.LT.kpb(ik-1,ir)) ik = ik - 1
      ENDIF
c
c Check that P and IK are not out of bounds:
c
      IF (p.LT.kpb(0,ir).OR.p.GT.kpb(nks(ir),ir)) THEN
        WRITE(0,*) 'ERROR (GetPos): Calculated P value is out ',
     .             'of bounds (IR = ',ir,'  P = ',p,')'   
        STOP
      ENDIF

      IF (ik.LT.1.OR.ik.GT.nks(ir)) THEN
        WRITE(0,*) 'ERROR (GetPos): IK value is out ',
     .             'of bounds (IR = ',ir,'  IK = ',ik,')'   
        STOP
      ENDIF

      RETURN
      END
c
c = subroutine end: GetPos =============================================
c
c
c ======================================================================
c
c subroutine: GetDensity
c
c
c ======================================================================
c
c      RECURSIVE SUBROUTINE GetDensity(ir,thetav,density,ik,
c     >                                quant,targval)
c
      SUBROUTINE GetDensity(ir,thetav,density,ik,
     >                                quant,targval,denoptin)
      implicit none 
c
c Commons:
c 
      INCLUDE 'params'                                                  
      INCLUDE 'cgeom'                                                  
      INCLUDE 'comtor'
c
c Input:
c
      REAL    thetav,quant(maxnks,maxnrs),targval(maxnds)
      INTEGER ir,denoptin
c
c Output:
c
      REAL    density
      INTEGER ik
      REAL    CellWidth

      INTEGER denopt,idum1,idum2,idum3
      REAL    den1,den2,den3,den4,p,p1,p2,theta1

      theta1 = thetav
c Option currently hard coded: not anymore
      denopt = denoptin

c
c Check for valid IR value (virtual rings not allowed):
c 
      IF (ir.EQ.1.OR.ir.EQ.irtrap.OR.ir.EQ.irwall) THEN
        WRITE(0,*) 'ERROR (GetDensity): Virtual ring (IR = ',ir,')'
        STOP
      ENDIF
c
c Make sure that THETA1 is not out of bounds:
c

c
c Find IK and P:
c

      IF (denopt.EQ.1) THEN
c
c Assign value at cell center:
c
        CALL GetPos(ir,theta1,p,ik)
        density = quant(ik,ir)

      ELSEIF (denopt.EQ.2) THEN
c
c Linear interpolation between cell centers:
c
        CALL GetPos(ir,theta1,p,ik)
c
        IF (p.LT.kps(ik,ir).AND.ik.EQ.1) THEN
          p1   = kps (ik,ir)
          p2   = kpb (0 ,ir) 
          den1 = quant(ik,ir)
          if (ir.ge.irsep) then 
             den2 = targval(idds(ir,2))
          else 
             den2 = (quant(nks(ir)-1,ir)+quant(ik,ir))/2.0
          endif
        ELSEIF (p.GT.kps(ik,ir).AND.ik.EQ.nks(ir)) THEN
          p1   = kpb (ik,ir) 
          p2   = kps (ik,ir)         
          if (ir.ge.irsep) then 
             den1 = targval(idds(ir,1))
          else
             den1 = (quant(2,ir)+quant(ik,ir))/2.0
          endif
          den2 = quant(ik,ir)
        ELSE
          IF (p.GE.kps(ik,ir)) THEN
            p1   = kps (ik+1,ir)
            p2   = kps (ik  ,ir)
            den1 = quant(ik+1,ir)
            den2 = quant(ik  ,ir)
          ELSE
            p1   = kps (ik  ,ir)
            p2   = kps (ik-1,ir)
            den1 = quant(ik  ,ir)
            den2 = quant(ik-1,ir)
          ENDIF
        ENDIF

        density = den2 + (p - p2) / (p1 - p2) * (den1 - den2)      
c
c Error checking:
c
      elseif (denopt.eq.3) then 

         density = quant(ik,ir)

      ELSE
        WRITE(0,*) 'ERROR: Illegal option in GetDensity'
      ENDIF

      RETURN
      END
c
c = subroutine end: GetDensity =========================================
c
c
c ======================================================================
c
c function: Quant2Grad
c
c Should I be working with P or S?  Does it really matter?
c
c
c ======================================================================
c
c
C     IPP/01 - Krieger: FUJI f90 compiler is very strict and requires
C     standard implementation of recursive functions with RESULT keyword
c     NOTE: I think this has to be a compiler bug - if there is no other
c           workaround we will change the code - however, it is possible
c           that the problem is related to some other code aspect - other
c           than recursion. Could also try this by specifying the function
c           name in the result statement which should be the default 
c           behaviour in any case.
c
c
      RECURSIVE REAL FUNCTION Quant2Grad(ir,s,quant,targval) ! RESULT(res)
     .               RESULT (res_Quant2Grad)
c     .               RESULT (Quant2Grad)

      IMPLICIT none
c     
c Commons:
c
      INCLUDE 'params'                                                  
      INCLUDE 'cgeom'                                                  
      INCLUDE 'comtor'      
      INCLUDE 'slcom'      
c 
c Input:
c
      INTEGER ir
      REAL    s,quant(maxnks,maxnrs),targval(maxnds)
c
c Locals:
c
      REAL GetTheta,GetS,CellWidth
c slmod begin
      INTEGER ik1
c slmod end
      INTEGER ik,iki,iko,iri,iro
      REAL    thetav,den,deni,deno,gradi,grado,tgrado,tgradi,tgrad
      REAL    widi,wido,wid,s1
c
c Output:
c
c      REAL    res_Quant2Grad
c
c
c Check for valid IR value:
c 
      IF (ir.LT.irsep.OR.ir.EQ.irtrap.OR.ir.EQ.irwall) THEN
        WRITE(6,*) 'ERROR (Quant2Grad): Cannot calculate second order ',
     .             'gradient for virtual rings (IR = ',ir,')'
c        quant2grad = 1.0
c
        res_quant2grad = 1.0
c
        return
      ELSEIF (nrs-irtrap.LT.2.OR.irwall-irsep.LT.2) THEN
c        quant2grad = 1.0
c
        res_quant2grad = 1.0
c
        return
      ENDIF
c
c Since calculating the second derivative for points along a ring
c requires information about the nighbouring rings, finding the
c derivative for rings next to virtual rings is a problem.  The
c approximation used below takes the gradient from the neighbouring
c ring (which is then one ring removed from the virtual ring):
c
c slmod begin
c
c...  Find IK value for the cell for which the gradient is calculated:
      IF (grdnmod.NE.0) THEN
        DO ik1 = 1, nks(ir)
          IF (ksb(ik1-1,ir).LE.s.AND.ksb(ik1,ir).GT.s) EXIT
        ENDDO
        IF (ik1.EQ.nks(ir)+1) 
     .    CALL ER('Quant2Grad','Cannot find IK index',*99)
      ENDIF
c slmod end
      IF (ir.EQ.2.OR.ir.EQ.irtrap+1) THEN
c... BUG - OCT 2, 2000
        s1 = (s - ksb(0,ir)) * (ksb(nks(ir+1),ir+1) - ksb(0,ir+1)) /
     .                         (ksb(nks(ir)  ,ir  ) - ksb(0,ir  ))
c        s = (s - ksb(0,ir)) * (ksb(nks(ir+1),ir+1) - ksb(0,ir+1)) /
c     .                         (ksb(nks(ir)  ,ir  ) - ksb(0,ir  ))
c        res = Quant2Grad(ir+1,s1,quant,targval)
c        Quant2Grad = Quant2Grad(ir+1,s1,quant,targval)  ! gfortran
c
        res_Quant2Grad = Quant2Grad(ir+1,s1,quant,targval)
c
        RETURN
      ELSEIF (ir.EQ.irwall-1) THEN
        s1 = (s - ksb(0,ir)) * (ksb(nks(ir-1),ir-1) - ksb(0,ir-1)) /
     .                         (ksb(nks(ir)  ,ir  ) - ksb(0,ir  ))
c        s = (s - ksb(0,ir)) * (ksb(nks(ir-1),ir-1) - ksb(0,ir-1)) /
c     .                         (ksb(nks(ir)  ,ir  ) - ksb(0,ir  ))
c        res = Quant2Grad(ir-1,s1,quant,targval)
c        Quant2Grad = Quant2Grad(ir-1,s1,quant,targval) !gfortran
c
        res_Quant2Grad = Quant2Grad(ir+1,s1,quant,targval)
c
        RETURN
c slmod begin
c     IPP/08 Krieger - Usual problem: ik1 exists only if grdnmod.ne.0
c     INTEL compiler runtime system always evaluates complete set of
c     conditions -> requires nested if clause
*     ELSEIF (grdnmod.NE.0.AND.
*    .        ir.LT.irwall.AND.irouts(ik1,ir).EQ.irwall) THEN
      ELSEIF (grdnmod.NE.0) THEN
        IF (ir.LT.irwall.AND.irouts(ik1,ir).EQ.irwall) THEN
          s1 = (s - ksb(0,ir)) * (ksb(nks(ir-1),ir-1) - ksb(0,ir-1)) /
     .                           (ksb(nks(ir)  ,ir  ) - ksb(0,ir  ))
c         s  = (s - ksb(0,ir)) * (ksb(nks(ir-1),ir-1) - ksb(0,ir-1)) /
c     .                          (ksb(nks(ir)  ,ir  ) - ksb(0,ir  ))
c          res = Quant2Grad(ir-1,s1,quant,targval)
c          Quant2Grad = Quant2Grad(ir-1,s1,quant,targval)  ! gfortran
c
        res_Quant2Grad = Quant2Grad(ir-1,s1,quant,targval)
c
          RETURN
        ENDIF
c slmod end
      ENDIF
c
c Find the THETAV value corresponding to given S:
c
      thetav = GetTheta(ir,s)
c
c Find density on ring IR, the IK index and the IR values 
c for the neighbouring cells:
c
      CALL GetDensity(ir,thetav,den,ik,quant,targval,2)

      iri = irins (ik,ir)
      iro = irouts(ik,ir)
c
c Find the density on neighbouring rings:
c
      IF (ir.EQ.irsep.AND.ik.GT.ikto) THEN
        IF (ik.GE.ikti) THEN
c
c         The THETA co-ordinate needs to be adjusted when referencing
c         the PFZ from the high index leg of the SOL:
c
          CALL GetDensity(iri,thetav-dthetg,deni,iki,quant,targval,2)
        ELSE
c
c         If IK,IR is on the separatrix ring and between the cut points,
c         then the poorly understood characteristics of cross-field 
c         transport between the core and the SOL become a problem.  So,
c         approximate the second derivative in this region by taking
c         the derivative on the neighbouring SOL ring:
c

c          res = Quant2Grad(ir+1,GetS(ir+1,thetav),quant,targval)
c          Quant2Grad = Quant2Grad(ir+1,GetS(ir+1,thetav),quant,targval) ! gfortran
c
          res_Quant2Grad = 
     >                Quant2Grad(ir+1,GetS(ir+1,thetav),quant,targval)
c
          RETURN
        ENDIF
      ELSE  
        CALL GetDensity(iri,thetav,deni,iki,quant,targval,2)
      ENDIF

      IF (ir.EQ.nrs.AND.ik.GT.ikto) THEN
c
c       The THETA co-ordinate needs to be adjusted when referencing
c       the high index leg of the SOL from the PFZ:
c
        CALL GetDensity(iro,thetav+dthetg,deno,iko,quant,targval,2)
      ELSE
        CALL GetDensity(iro,thetav,deno,iko,quant,targval,2)
      ENDIF    
c
c Determine cell widths:
c      
      wid  = CellWidth(ik ,ir ,CENTER,TOTAL)
      widi = CellWidth(iki,iri,CENTER,TOTAL)
      wido = CellWidth(iko,iro,CENTER,TOTAL)

c
c Approximate second derivative:
c
      gradi = (deni - den) / (0.5 * (widi + wid))
      grado = (deno - den) / (0.5 * (wido + wid))

c      res = (gradi + grado) / wid
c      Quant2Grad = (gradi + grado) / wid  ! gfortran
c
      res_Quant2Grad = (gradi + grado) / wid
c
c      Quant2Grad = (((deni - den) / (0.5 * (wid + widi))) +
c     .            ((deno - den) / (0.5 * (wid + wido)))) / wid

      tgradi = (quant(ikins(ik,ir),irins(ik,ir)) - quant(ik,ir)) /
     .       (0.5 * (CellWidth(ikins(ik,ir),irins(ik,ir),CENTER,TOTAL) +
     .               CellWidth(ik          ,ir          ,CENTER,TOTAL)))

      tgrado = (quant(ikouts(ik,ir),irouts(ik,ir)) - quant(ik,ir)) /
     .     (0.5 * (CellWidth(ikouts(ik,ir),irouts(ik,ir),CENTER,TOTAL) +
     .             CellWidth(ik           ,ir           ,CENTER,TOTAL)))

      tgrad = (tgradi + tgrado) / wid

c
c
c
      if (cprint.eq.8.or.cprint.eq.9) then     
  
         WRITE(SLOUT,'(2I3,A,E12.6,A,E12.6,A,F8.3,A)') 
c     .     ik,ir,' >> GRAD: ',Quant2Grad,' <<     TGRAD: ',tgrad,
c     .     '   DIFF = ',(Quant2Grad - tgrad) / Quant2Grad * 100.0,' %'
c
     .     ik,ir,' >> GRAD: ',res_Quant2Grad,' <<     TGRAD: ',tgrad,
     .     '   DIFF = ',
     .     (res_Quant2Grad - tgrad) / res_Quant2Grad * 100.0,' %'
c
      endif
c
c
c      IF (((Quant2Grad - tgrad) / Quant2Grad * 100.0).GT.0.01) THEN
c
      IF (((res_Quant2Grad-tgrad)/res_Quant2Grad*100.0).GT.0.01)THEN
c
        if (cprint.eq.8.or.cprint.eq.9) then   

        WRITE(SLOUT,'(3x,2I3,3X,2I3,2X,2I3,2X,2I3,3F10.6)')
     .    iki,iri,ik,ir,iko,iro,ikins(ik,ir),ikouts(ik,ir),widi,wid,wido
      
        WRITE(SLOUT,'(3X,A,1F12.6)') 'KSS:   ',
     .    ksb(ik-1,ir)
        WRITE(SLOUT,'(10X,2F12.6)')
     .    kss(ik,ir),s
        WRITE(SLOUT,'(10X,1F12.6)')
     .    ksb(ik,ir)

        WRITE(SLOUT,'(3X,A,3F12.6)') 'THETA: ',
     .    thetag(iki-1,iri),thetag(ik-1,ir),thetag(iko-1,iro)
        WRITE(SLOUT,'(10X,4F12.6)')
     .    thetag(iki,iri)  ,thetag(ik,ir)  ,thetag(iko,iro),thetav
        WRITE(SLOUT,'(10X,3F12.6)')
     .    thetag(iki+1,iri),thetag(ik+1,ir),thetag(iko+1,iro)

        WRITE(SLOUT,'(3X,A,3E12.6)') 'KBNS:  ',
     .    quant(iki-1,iri),quant(ik-1,ir),quant(iko-1,iro)
        WRITE(SLOUT,'(10X,6E12.6)')
     .    quant(iki,iri)  ,quant(ik,ir)  ,quant(iko,iro), deni,den,deno
        WRITE(SLOUT,'(10X,3E12.6)')
     .    quant(iki+1,iri),quant(ik+1,ir),quant(iko+1,iro)

        endif

      ENDIF


      RETURN
c slmod begin
99    STOP 
c slmod end
      END
c
c = function end: Quant2Grad ===========================================
c

c ======================================================================
c
c subroutine: QuantGrad
c
c ======================================================================
c
      subroutine QuantGrad(ir,s,quant,targval,gradi,grado,denopt)

      IMPLICIT none
c     
c Commons:
c
      INCLUDE 'params'                                                  
      INCLUDE 'cgeom'                                                  
      INCLUDE 'comtor'      
      INCLUDE 'slcom'      
c 
c Input:
c
      INTEGER ir,denopt
      REAL    s,quant(maxnks,maxnrs),targval(maxnds)
      real    gradi, grado
c
c Locals:
c
      REAL GetTheta,GetS,CellWidth

      INTEGER ik,iki,iko,iri,iro
      REAL    thetav,den,deni,deno,tgrado,tgradi,tgrad
      REAL    widi,wido,wid


c
c Check for valid IR value:
c 
      IF (ir.eq.1.OR.ir.EQ.irtrap.OR.ir.EQ.irwall) THEN
        WRITE(6,*) 'ERROR (QuantGrad): Cannot calculate ',
     .             'gradient for virtual rings (IR = ',ir,')'
        gradi  = 0.0
        grado  = 0.0
        return
      ENDIF
c
c Find the THETAV value corresponding to given S:
c
      thetav = GetTheta(ir,s)
c
c Find density and cell width on ring IR, the IK index and the IR values 
c for the neighbouring cells:
c
c
      CALL GetDensity(ir,thetav,den,ik,quant,targval,2)
c
      wid  = CellWidth(ik ,ir ,CENTER,TOTAL)
c
      iri = irins (ik,ir)
      iro = irouts(ik,ir)
c
c     Set initial values - may be changed in GetDensity
c
      iki = ikins(ik,ir)
      iko = ikouts(ik,ir)
c
c
c Since calculating the derivative for points along a ring
c requires information about the nighbouring rings, finding the
c derivative for rings next to virtual rings is a problem.
c The derivative on the side facing the virtual ring is set to zero.
c
c
c Find the density on and cell width neighbouring rings:
c
c     Inward
c
      IF (ir.ne.2.and.ir.ne.irtrap+1) THEN
c
        IF (ir.EQ.irsep.AND.ik.Ge.ikti) THEN
c	
c         The THETA co-ordinate needs to be adjusted when referencing
c         the PFZ from the high index leg of the SOL:
c	
          CALL GetDensity(iri,thetav-dthetg,deni,iki,quant,targval,
     >                    denopt)
c	
        ELSE  
          CALL GetDensity(iri,thetav,deni,iki,quant,targval,denopt)
        ENDIF
c
        widi = CellWidth(iki,iri,CENTER,TOTAL)
c
      endif 
c
c     Outward
c
      IF (ir.ne.irwall-1) THEN
c
        IF (ir.EQ.nrs.AND.ik.GT.ikto) THEN
c	
c         The THETA co-ordinate needs to be adjusted when referencing
c         the high index leg of the SOL from the PFZ:
c	
          CALL GetDensity(iro,thetav+dthetg,deno,iko,quant,targval,
     >                    denopt)
	
        ELSE
          CALL GetDensity(iro,thetav,deno,iko,quant,targval,denopt)
        ENDIF    
c
        wido = CellWidth(iko,iro,CENTER,TOTAL)
c
      endif
c
c Approximate derivatives:
c
      IF (ir.eq.2.or.ir.eq.irtrap+1) THEN
        gradi = 0.0
      else 
        gradi = (deni - den) / (0.5 * (widi + wid))
      endif
c
      write(6,'(a,4i4,5(1x,g12.5))') 'GRADI:',ik,ir,iki,iri,gradi,
     >             deni,den,widi,wid
c
      IF (ir.eq.irwall-1) THEN
        grado = 0.0
      else
        grado = (den - deno) / (0.5 * (wido + wid))
      endif
c
      write(6,'(a,4i4,5(1x,g12.5))') 'GRADO:',ik,ir,iko,iro,grado,
     >             deno,den,wido,wid
c
      RETURN
      END
c
c = function end: QuantGrad ===========================================
c



