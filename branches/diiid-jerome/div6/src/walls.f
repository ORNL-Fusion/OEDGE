c     -*-Fortran-*-
c
      SUBROUTINE DOTARG
      use debug_options
      IMPLICIT NONE
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "CGEOM"
      include 'cgeom'
C     INCLUDE "CNEUT"
c      include 'cneut'
C     INCLUDE "COMTOR"
      include 'comtor'
C     INCLUDE "CIONIZ"
      include 'cioniz'
C     INCLUDE "READER"
      include 'reader'
C     INCLUDE "DYNAM5"
      include 'dynam5'
      include 'printopt'
c slmod begin
      include 'slcom'

c...  For COSTET2 calcualtion (temporary):
      INTEGER CalcPoint
      INTEGER ii
      LOGICAL discrepancy
      REAL sadj,shyp,costet2(MAXNRS*2)
      REAL*8 a1,a2,b1,b2,c1,c2,t(4)
c slmod end
      INTEGER IK,IR,K,NP,L,J,I,NR,NC,NXW,IEXTRA,JK,JR,MIZS,IZ,IERR,ID
      INTEGER IX,IY,IKIN,IKOUT,IRIN,IROUT,MKS,NP1,ICOUNT,INEXT,KNEXT
      INTEGER IND,IW,KP,in
      REAL    DD1,DD2,DELTAZ,DELTAR,THIN,THOUT
      real    xcos
      INTEGER RES
      INTEGER FINDNUM
      REAL ATAN2C
      EXTERNAL FINDNUM,ATAN2C
c
      integer ikshift,vertex1,vertex2    
c
      real isepdist,osepdist,rspo,zspo,rspi,zspi
C
C-----------------------------------------------------------------------
C     CALCULATE SET OF POINTS ALONG TARGET PLATES
C-----------------------------------------------------------------------
C
C
C     IN THE ARRAY IDDS - THE INNER PLATE IS 1 AND THE OUTER IS 2
c     The inner plate is assumed to be for IK values = nks(ir) and
c     the outer plate is for IK values = 1.
c
      if (cprint.eq.3.or.cprint.eq.9) then 
         WRITE(6,*) 'DO TARGETS:'
      endif

      call pr_trace('WALLS','DOTARG START')

c
c     Target options:
c     0 - targets at at the last set of real points on each ring. The
c         virtual points are deleted.
c     4 - targets at the last points on each ring. (These are the
c         virtual points.)
c
c slmod begin
      IF (nbr.GT.0.OR.cgridopt.EQ.LINEAR_GRID.or.
     >    cgridopt.eq.RIBBON_GRID) THEN
c...    Generalized grid code:
        CALL BuildTargets
        CALL OutputData(85,'IN DOTARG')
c
c      if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then     
c slmod end
      elseif (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
c
c     JET or ASDEX grids
c
        IF (CTARGOPT.EQ.0.or.ctargopt.eq.4) THEN
           ID = 0
           DO 460 IR = IRWALL, IRSEP, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RP(ID) = RS(NKS(IR),IR)
             ZP(ID) = ZS(NKS(IR),IR)
  460      CONTINUE
           DO 462 IR = NRS, IRTRAP, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RP(ID) = RS(NKS(IR),IR)
             ZP(ID) = ZS(NKS(IR),IR)
  462      CONTINUE
           NDSIN = ID
           DO 464 IR = IRTRAP, NRS
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RP(ID) = RS(1,IR)
             ZP(ID) = ZS(1,IR)
  464      CONTINUE
           DO 466 IR = IRSEP, IRWALL
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RP(ID) = RS(1,IR)
             ZP(ID) = ZS(1,IR)
  466      CONTINUE
           NDS = ID
C
C     Target options:
c     1 - target midway between last real and virtual point. Target
c         points calculated in tau module. Virtual points deleted.
c     2 - target points specified in input data and connect to last
c         real point. Virtual points deleted.
c     3 - target points are loaded from the loadgeo routine. They
c         connect to the last real points and vitual points are deleted.
c     5 - target points specified in the input data. They connect to
c         the virtual points.
c     6 - the polygon edge along the target is used to define the target.
c         The target coordinates have been defined in PLATCO to be
c         the midpoint of the the polygon side. The virtual point has
c         been removed.
C
         ELSEIF (ctargopt.eq.1.or.CTARGOPT.EQ.2.or.
     >           ctargopt.eq.3.or.ctargopt.eq.5.or.
     >           ctargopt.eq.6) THEN
           ID = 0
           DO 480 IR = IRWALL, IRSEP, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,2)
             ZP(ID) = PLATCO(RES,3)
  480      CONTINUE
         
         


           DO 482 IR = NRS, IRTRAP, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,2)
             ZP(ID) = PLATCO(RES,3)
  482      CONTINUE
         
c        
c          Add target element for FRC grid
c        
           if (sonnet_grid_sub_type.eq.1) then 
             IR = IRSEP -1 
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,2)
             ZP(ID) = PLATCO(RES,3)
           endif
        
           NDSIN = ID
        
        
           DO 484 IR = IRTRAP, NRS
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,4)
             ZP(ID) = PLATCO(RES,5)
  484      CONTINUE
         
         
c        
c          Add target element for FRC grid
c        
           if (sonnet_grid_sub_type.eq.1) then 
             IR = IRSEP -1 
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,4)
             ZP(ID) = PLATCO(RES,5)
           endif
c        
         
           DO 486 IR = IRSEP, IRWALL
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,4)
             ZP(ID) = PLATCO(RES,5)
  486      CONTINUE
           NDS = ID
           ndsin2 = nds
           ndsin3 = ndsin
         ENDIF
c      
c      
        ELSEIF (cgridopt.eq.2) then
c    
c       ITER grid
c    
        IF (CTARGOPT.EQ.0.or.ctargopt.eq.4) THEN
           ID = 0
c    
c          Lower right target
c    
           DO 1460 IR = IRWALL, IRSEP2, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RP(ID) = RS(NKS(IR),IR)
             ZP(ID) = ZS(NKS(IR),IR)
 1460      CONTINUE
           DO 1462 IR = NRS2, IRTRAP2, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RP(ID) = RS(NKS(IR),IR)
             ZP(ID) = ZS(NKS(IR),IR)
 1462      CONTINUE
c        
           NDSIN = ID
c        
c          Lower left target
c        
           DO 1464 IR = IRTRAP2, NRS2
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RP(ID) = RS(1,IR)
             ZP(ID) = ZS(1,IR)
 1464      CONTINUE
           DO 1466 IR = IRSEP, IRWALL2
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RP(ID) = RS(1,IR)
             ZP(ID) = ZS(1,IR)
 1466      CONTINUE
c         
           ndsin2 = id
c         
c          Upper left target
c
           DO 1470 IR = IRWALL2, IRSEP, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RP(ID) = RS(NKS(IR),IR)
             ZP(ID) = ZS(NKS(IR),IR)
 1470      CONTINUE
           DO 1472 IR = NRS, IRTRAP, -1
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RP(ID) = RS(1,IR)
             ZP(ID) = ZS(1,IR)
 1472      CONTINUE
c        
           NDSIN3 = ID
c        
c            Upper right target
c       
           DO 1474 IR = IRTRAP, NRS
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = nks(ir)
             IRDS(ID) = IR
             RP(ID) = RS(nks(ir),IR)
             ZP(ID) = ZS(nks(ir),IR)
 1474      CONTINUE
           DO 1476 IR = IRSEP2, IRWALL
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RP(ID) = RS(1,IR)
             ZP(ID) = ZS(1,IR)
 1476      CONTINUE
c       
           NDS = ID
C
C     Target options:
c     1 - target midway between last real and virtual point. Target
c         points calculated in tau module. Virtual points deleted.
c     2 - target points specified in input data and connect to last
c         real point. Virtual points deleted.
c     3 - target points are loaded from the loadgeo routine. They
c         connect to the last real points and vitual points are deleted.
c     5 - target points specified in the input data. They connect to
c         the virtual points.
C
         ELSEIF (ctargopt.eq.1.or.CTARGOPT.EQ.2.or.
     >           ctargopt.eq.3.or.ctargopt.eq.5.or.ctargopt.eq.6) THEN
           ID = 0
c        
c          Lower right target
c        
           DO 1480 IR = IRWALL, IRSEP2, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,2)
             ZP(ID) = PLATCO(RES,3)
 1480      CONTINUE
           DO 1482 IR = NRS2, IRTRAP2, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,2)
             ZP(ID) = PLATCO(RES,3)
 1482      CONTINUE
c        
           NDSIN = ID
c        
c          Lower left target
c        
           DO 1484 IR = IRTRAP2, NRS2
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,4)
             ZP(ID) = PLATCO(RES,5)
 1484      CONTINUE
           DO 1486 IR = IRSEP, IRWALL2
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,4)
             ZP(ID) = PLATCO(RES,5)
 1486      CONTINUE
c        
           NDSin2 = ID
c        
c          Upper left target
c        
           DO 1490 IR = IRWALL2, IRSEP, -1
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = NKS(IR)
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,2)
             ZP(ID) = PLATCO(RES,3)
 1490      CONTINUE
           DO 1492 IR = NRS, IRTRAP, -1
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,4)
             ZP(ID) = PLATCO(RES,5)
 1492      CONTINUE
c
           NDSIN3 = ID
c       
c          Upper right target
c       
           DO 1494 IR = IRTRAP, NRS
             ID = ID + 1
             IDDS(IR,1) = ID
             IKDS(ID) = nks(ir)
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,2)
             ZP(ID) = PLATCO(RES,3)
 1494      CONTINUE
           DO 1496 IR = IRSEP2, IRWALL
             ID = ID + 1
             IDDS(IR,2) = ID
             IKDS(ID) = 1
             IRDS(ID) = IR
             RES = FINDNUM(PLATCO,NPLAT,IR)
             IF (RES.EQ.0) THEN
                WRITE(6,*) 'ERROR ERROR ERROR: NO PLATE'//
     >              ' COORDINATES SPECIFIED FOR RING:', IR
                STOP
             ENDIF
             RP(ID) = PLATCO(RES,4)
             ZP(ID) = PLATCO(RES,5)
 1496      CONTINUE
c
           NDS = ID
c
        ENDIF
c
c     End for grid options
c
      ENDIF
c
C
C     PLATE POINTS WHICH HAVE A LENGTH OF ZERO BELONG TO THE VIRTUAL
C     RING WHEN THIS REGION IS BEING EXCLUDED - I.E. WHEN WALLS ARE
C     SET TO HALF WAY BETWEEN THE OUTER TWO RINGS.
C
c     added for AUG (control of target options)
c
      if (cprint.eq.3.or.cprint.eq.9)
     >  write (6,'(1x,a)')
     >               'Rtarget  Ztarget  Index   Remaining stuff??'
c
      DO 468 ID = 1, NDS
c
c       Code for target option 0 to 5
c
        if (ctargopt.eq.0.or.ctargopt.eq.1.or.ctargopt.eq.2.or.
     >      ctargopt.eq.3.or.ctargopt.eq.4.or.ctargopt.eq.5) then
c
        IF (ID.EQ.1.OR.ID.EQ.NDSIN+1.or.
     >      id.eq.ndsin2+1.or.id.eq.ndsin3+1) THEN
           DELTAR = RP(ID+1)-RP(ID)
           DELTAZ = ZP(ID+1)-ZP(ID)
           IF (DELTAR.NE.0.0.AND.DELTAZ.NE.0.0)
     >        THIN  = ATAN2C (DELTAZ,DELTAR) - 0.5 * PI
           THOUT = THIN
           IF (CIONR.EQ.0.or.cionr.eq.2) then
              DD1 = 0.0
           ELSE
              DD1 = SQRT(DELTAR**2+DELTAZ**2)
           ENDIF
           DD2 = 0.0
        ELSEIF (ID.EQ.NDS.OR.ID.EQ.NDSIN.or.
     >          id.eq.ndsin2.or.id.eq.ndsin3) THEN
           DELTAR = RP(ID-1)-RP(ID)
           DELTAZ = ZP(ID-1)-ZP(ID)
           IF (DELTAR.NE.0.0.AND.DELTAZ.NE.0.0)
     >        THOUT  = ATAN2C (DELTAZ,DELTAR) + 0.5 * PI
           THIN = THOUT
           DD1 = 0.0
           IF (CIONR.EQ.0.or.cionr.eq.2) then
              DD2 = 0.0
           ELSE
              DD2 = SQRT(DELTAR**2+DELTAZ**2)
           ENDIF
        ELSE
           DELTAR = RP(ID+1)-RP(ID)
           DELTAZ = ZP(ID+1)-ZP(ID)
           IF (DELTAR.NE.0.0.AND.DELTAZ.NE.0.0)
     >        THIN  = ATAN2C (DELTAZ,DELTAR) - 0.5 * PI
           DD1 = SQRT(DELTAR**2+DELTAZ**2)
           DELTAR = RP(ID-1)-RP(ID)
           DELTAZ = ZP(ID-1)-ZP(ID)
           IF (DELTAR.NE.0.0.AND.DELTAZ.NE.0.0)
     >        THOUT  = ATAN2C (DELTAZ,DELTAR) + 0.5 * PI
           DD2 = SQRT(DELTAR**2+DELTAZ**2)
        ENDIF
c
        IF (THIN.GT.PI) THEN
           THIN = THIN - 2.0 * PI
        ELSEIF (THIN.LT.-PI) THEN
           THIN = THIN + 2.0 * PI
        ENDIF
        IF (THOUT.GT.PI) THEN
           THOUT = THOUT - 2.0 * PI
        ELSEIF (THOUT.LT.-PI) THEN
           THOUT = THOUT + 2.0 * PI
        ENDIF
        THETAS(ID) = 0.5 * (THIN + THOUT)
        DDS(ID)    = 0.5 * (DD1 + DD2)
c
        if (cprint.eq.3.or.cprint.eq.9)
     >      WRITE (6,9043) RP(ID),ZP(ID),id,ikds(id),irds(id),
     >         RADDEG*THIN,RADDEG*THOUT,RADDEG*THETAS(ID),DDS(ID)
c
c       For target option 6 - which uses the polygon edges.
c
        elseif (ctargopt.eq.6) then
c
c       For target option 6 - which uses the grid polygons to define the
c       target - then the code may not be valid - check for a non-zero
c       value of NVERTP - if it is 0.0 - then the segment length
c       must also be zero.
c
c
c       Calculate the inclination angle of each target
c       segment and the length of each target segment and
c       save these in the arrays THETAS(ID) and DDS(ID)
c       respectively.
c
        kp = korpg(ikds(id),irds(id))
        ir = irds(id)
c
c       Modified to check for existence of a cell - if a cell exists calculate
c       values - for all boundary cells the indices in korpg should have been
c       zeroed out - CHANGED since some transport code seems to need 
c       some defined quantities for the boundary rings.  
c 
c       Boundary rings are ir = 1, ir = irwall and ir = irtrap - check to see if 
c       irds(id) matches any of these
c
c        if (kp.ne.0.and.nvertp(kp).gt.0) then
c
        if (id.ge.1.and.id.le.ndsin) then
c
c          IPP/08 Krieger - ensure index of nvertp is not zero
           if (kp.ne.0.and.nvertp(max(1,kp)).gt.0.and.
     >         ir.ne.1.and.
     >         ir.ne.irwall.and.
     >         ir.ne.irtrap) then
              DELTAR = RVERTP(4,KP) - RVERTP(3,KP)
              DELTAZ = ZVERTP(4,KP) - ZVERTP(3,KP)
              DDS(ID) = SQRT(DELTAR**2 + DELTAZ**2)
              THETAS(ID) = ATAN2C(DELTAZ,DELTAR) - 0.5 * PI
              IF (THETAS(ID).LT.-PI) THETAS(ID) = THETAS(ID) + 2.0 * PI
           elseif (cionr.eq.0.or.cionr.eq.2) then
              dds(id) = 0.0
              thetas(id) = 0.0
           elseif (cionr.eq.1.and.(id.eq.1.or.id.eq.ndsin)) then
              if (id.eq.1) then
                 DELTAR = RP(ID+1)-RP(ID)
                 DELTAZ = ZP(ID+1)-ZP(ID)
              else
                 DELTAR = RP(ID)-RP(ID-1)
                 DELTAZ = ZP(ID)-ZP(ID-1)
              endif
              thetas(id) =  ATAN2C (DELTAZ,DELTAR) - 0.5 * PI
              dds(id) = SQRT(DELTAR**2+DELTAZ**2)
           else
              dds(id) = 0.0
              thetas(id) = 0.0
           endif
c
        elseif (id.ge.ndsin+1.and.id.le.nds) then
c
c          IPP/08 Krieger - ensure index of nvertp is not zero
           if (kp.ne.0.and.nvertp(max(1,kp)).gt.0.and.
     >         ir.ne.1.and.
     >         ir.ne.irwall.and.
     >         ir.ne.irtrap) then
              DELTAR = RVERTP(2,KP) - RVERTP(1,KP)
              DELTAZ = ZVERTP(2,KP) - ZVERTP(1,KP)
              DDS(ID) = SQRT(DELTAR**2 + DELTAZ**2)
              THETAS(ID) = ATAN2C(DELTAZ,DELTAR) - 0.5 * PI
              IF (THETAS(ID).LT.-PI) THETAS(ID) = THETAS(ID) + 2.0 * PI
           elseif (cionr.eq.0.or.cionr.eq.2) then
              dds(id) = 0.0
              thetas(id) = 0.0
           elseif (cionr.eq.1.and.(id.eq.ndsin+1.or.id.eq.nds)) then
              if (id.eq.ndsin+1) then
                 DELTAR = RP(ID+1)-RP(ID)
                 DELTAZ = ZP(ID+1)-ZP(ID)
              else
                 DELTAR = RP(ID)-RP(ID-1)
                 DELTAZ = ZP(ID)-ZP(ID-1)
              endif
              thetas(id) =  ATAN2C (DELTAZ,DELTAR) - 0.5 * PI
              dds(id) = SQRT(DELTAR**2+DELTAZ**2)
           else
              dds(id) = 0.0
              thetas(id) = 0.0
           endif
c
        else
           dds(id) = 0.0
           thetas(id) = 0.0
        endif
c
        if (cprint.eq.3.or.cprint.eq.9)
     >     WRITE (6,9043) RP(ID),ZP(ID),id,ikds(id),irds(id),
     >        RADDEG*THETAS(ID),DDS(ID)
c
        endif
c
  468 CONTINUE
c
c     changed for AUG (better readability)
c
 9043 FORMAT(1X,2(F7.4,2X),3I4,2X,3(F8.2,2X),F8.5)
c
c 9043 FORMAT(1X,2F11.6,1x,3I5,2X,3F8.3,F11.6)
C
C  CALCULATE TARGET SEGMENT LENGTHS AND ANGLES FROM THE PLASMA POLYGON
C  VERTEX POSITIONS
c
c  Applicable for JET and Sonnet grids - calculate target angles.
c
C
      if (cgridopt.eq.0.or.cgridopt.eq.3.or.
     >     cgridopt.eq.LINEAR_GRID.or.cgridopt.eq.RIBBON_GRID) then
c
C
C  CALCULATE 1.TARGET SEGMENT LENGTHS AND ANGLES FROM THE PLASMA
C              POLYGON VERTEX POSITIONS : DDS2(ID) , THETAS2(ID)
C            2.COSINES FOR DEVIATIONS FROM ORTHOGONALITY OF ANGLE
C              BETWEEN FIELD LINES AND TARGET PLATES - FROM ANGLES
C              OF FIRST AND LAST FIELD LINE SEGMENTS AND EITHER
C              ANGLES (THETAS2) OR ABOVE ANGLES (THETAS) : COSTET(ID)
C
C  NOTE: There are no plasma polygons defined on virtual rings
C
         if (cprint.eq.3.or.cprint.eq.9) 
     >         WRITE(6,9046)
c
c NOTE: Execute this entire block of code for ID = 1,NDS - don't recalculate the
c       ID index since it makse no sense to do so - the targets have been defined 
c       above.    
c
c
        do id = 1,nds
c       
           ir = irds(id)
           ik = ikds(id) 
           if (ik.gt.nks(ir)/2) then 
              ikshift = -1
              vertex1 = 4 
              vertex2 = 3 
           else
              ikshift = 1 
              vertex1 = 2
              vertex2 = 1 
           endif
c
           IF( CTARGOPT.EQ.4 .OR. CTARGOPT.EQ.5 )THEN
             KP = KORPG(IK+ikshift,IR)
           ELSE
             KP = KORPG(IK,IR)
           ENDIF

c          IPP/08 Krieger - ensure index of nvertp is not zero
           if (kp.ne.0.and.nvertp(max(1,kp)).gt.0.and.
     >         ir.ge.irsep.and.
     >         ir.ne.irwall.and.
     >         ir.ne.irtrap.and.
     >         ir.ne.1) then 
           

              DELTAR = RVERTP(vertex1,KP) - RVERTP(vertex2,KP)
              DELTAZ = ZVERTP(vertex1,KP) - ZVERTP(vertex2,KP)
              DDS2(ID) = SQRT(DELTAR**2 + DELTAZ**2)
              THETAS2(ID) = ATAN2C(DELTAZ,DELTAR) - 0.5 * PI
              IF (THETAS2(ID).LT.-PI) 
     >            THETAS2(ID) = THETAS2(ID) + 2.0 * PI
C	    
              IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4)THEN
                 DELTAR = RS(IK,IR) - RS(IK+ikshift,IR)
                 DELTAZ = ZS(IK,IR) - ZS(IK+ikshift,IR)
              ELSE
                 DELTAR = RP(ID) - RS(IK,IR)
                 DELTAZ = ZP(ID) - ZS(IK,IR)
              ENDIF
              XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c	   
              IF( CTARGOPT.EQ.1 )THEN
                 COSTET(ID) = COS( XCOS - THETAS2(ID) )
              ELSE
                 COSTET(ID) = COS( XCOS - THETAS(ID) )
              ENDIF
C           
           else
c
              DDS2(ID)    = DDS(ID)
              THETAS2(ID) = THETAS(ID)
              IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
                 DELTAR = RS(IK,IR) - RS(IK+ikshift,IR)
                 DELTAZ = ZS(IK,IR) - ZS(IK+ikshift,IR)
              ELSE
                 DELTAR = RP(ID) - RS(IK,IR)
                 DELTAZ = ZP(ID) - ZS(IK,IR)
              ENDIF
              XCOS = ATAN2C(DELTAZ,DELTAR) - PI
              IF( CTARGOPT.EQ.1 )THEN
                 COSTET(ID) = COS( XCOS - THETAS2(ID) )
              ELSE
                 COSTET(ID) = COS( XCOS - THETAS(ID) )
              ENDIF

           endif
c
           if (cprint.eq.3.or.cprint.eq.9) 
     >        WRITE(6,9045) ID, RP(ID), ZP(ID), RADDEG*THETAS(ID),
     >           RADDEG*THETAS2(ID), DDS(ID), DDS2(ID), COSTET(ID)

        end do

c      endif




c
c--------- Code replace by above - will be deleted later     
c
c
c      ID = 0 
c
c     First target in two sections  
c
c      DO IR = IRWALL,IRSEP,-1
c
c        IK = NKS(IR)
c
c        ID = ID + 1
c
c        IF( CTARGOPT.EQ.4 .OR. CTARGOPT.EQ.5 )THEN
c          KP = KORPG(IK-1,IR)
c        ELSE
c          KP = KORPG(IK,IR)
c        ENDIF
c
c        if (kp.ne.0.and.nvertp(kp).gt.0) then 
c           
c           DELTAR = RVERTP(4,KP) - RVERTP(3,KP)
c           DELTAZ = ZVERTP(4,KP) - ZVERTP(3,KP)
c           DDS2(ID) = SQRT(DELTAR**2 + DELTAZ**2)
c           THETAS2(ID) = ATAN2C(DELTAZ,DELTAR) - 0.5 * PI
c           IF (THETAS2(ID).LT.-PI) THETAS2(ID) = THETAS2(ID) + 2.0 * PI
C	   
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4)THEN
c             DELTAR = RS(IK,IR) - RS(IK-1,IR)
c             DELTAZ = ZS(IK,IR) - ZS(IK-1,IR)
c           ELSE
c             DELTAR = RP(ID) - RS(IK,IR)
c             DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c	   
c           IF( CTARGOPT.EQ.1 )THEN
c             COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c             COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
C           
c        else
c
c           DDS2(ID)    = DDS(ID)
c           THETAS2(ID) = THETAS(ID)
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c              DELTAR = RS(IK,IR) - RS(IK-1,IR)
c              DELTAZ = ZS(IK,IR) - ZS(IK-1,IR)
c           ELSE
c              DELTAR = RP(ID) - RS(IK,IR)
c              DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c           IF( CTARGOPT.EQ.1 )THEN
c              COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c              COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
c
c        endif
c
c
c        WRITE(6,9045) ID, RP(ID), ZP(ID), RADDEG*THETAS(ID),
c     >     RADDEG*THETAS2(ID), DDS(ID), DDS2(ID), COSTET(ID)
c      ENDDO
c
c
c      DO IR = NRS,IRTRAP,-1
c
c        IK = NKS(IR)
c
c        ID = ID + 1
c
c        IF( CTARGOPT.EQ.4 .OR. CTARGOPT.EQ.5)THEN
c          KP = KORPG(IK-1,IR)
c        ELSE
c          KP = KORPG(IK,IR)
c        ENDIF
c
c        if (kp.ne.0.and.nvertp(kp).gt.0) then 
c           
c           
c           DELTAR = RVERTP(4,KP) - RVERTP(3,KP)
c           DELTAZ = ZVERTP(4,KP) - ZVERTP(3,KP)
c           DDS2(ID) = SQRT(DELTAR**2 + DELTAZ**2)
c           THETAS2(ID) = ATAN2C(DELTAZ,DELTAR) - 0.5 * PI
c           IF (THETAS2(ID).LT.-PI) THETAS2(ID) = THETAS2(ID) + 2.0 * PI
C	   
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c             DELTAR = RS(IK,IR) - RS(IK-1,IR)
c             DELTAZ = ZS(IK,IR) - ZS(IK-1,IR)
c           ELSE
c             DELTAR = RP(ID) - RS(IK,IR)
c             DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c	   
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c           IF( CTARGOPT.EQ.1 )THEN
c             COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c             COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
c           
c        else
c
c           DDS2(ID)    = DDS(ID)
c           THETAS2(ID) = THETAS(ID)
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c             DELTAR = RS(IK,IR) - RS(IK-1,IR)
c             DELTAZ = ZS(IK,IR) - ZS(IK-1,IR)  
c           ELSE
c             DELTAR = RP(ID) - RS(IK,IR)
c             DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c           IF( CTARGOPT.EQ.1 )THEN
c             COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c             COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
c
c        endif  
C	   
c        WRITE(6,9045) ID, RP(ID), ZP(ID), RADDEG*THETAS(ID),
c     >     RADDEG*THETAS2(ID), DDS(ID), DDS2(ID), COSTET(ID)
c      ENDDO
c
c
c     Next target in two sections
c
c
c      DO IR = IRTRAP,NRS
c
c        IK = 1
c
c        ID = ID + 1
c        IF ( CTARGOPT.EQ.4 .OR. CTARGOPT.EQ.5 )THEN
c          KP = KORPG(IK+1,IR)
c        ELSE
c          KP = KORPG(IK,IR)
c        ENDIF
c
c        if (kp.ne.0.and.nvertp(kp).gt.0) then 
c           
c           DELTAR = RVERTP(2,KP) - RVERTP(1,KP)
c           DELTAZ = ZVERTP(2,KP) - ZVERTP(1,KP)
c           DDS2(ID) = SQRT(DELTAR**2 + DELTAZ**2)
c           THETAS2(ID) = ATAN2C(DELTAZ,DELTAR) - 0.5 * PI
c           IF (THETAS2(ID).LT.-PI) THETAS2(ID) = THETAS2(ID) + 2.0 * PI
C          
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c             DELTAR = RS(IK,IR) - RS(IK+1,IR)
c             DELTAZ = ZS(IK,IR) - ZS(IK+1,IR)
c           ELSEIF (CTARGOPT.EQ.1.OR.CTARGOPT.EQ.2) THEN
c           ELSE
c             DELTAR = RP(ID) - RS(IK,IR)
c             DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c           IF( CTARGOPT.EQ.1 )THEN
c             COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c             COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
c           
c        else
c
c           DDS2(ID)    = DDS(ID)
c           THETAS2(ID) = THETAS(ID)
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c             DELTAR = RS(IK,IR) - RS(IK+1,IR)
c             DELTAZ = ZS(IK,IR) - ZS(IK+1,IR)
c           ELSE
c             DELTAR = RP(ID) - RS(IK,IR)
c             DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c           IF( CTARGOPT.EQ.1 )THEN
c             COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c             COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
c      
c        endif 
C        
c        WRITE(6,9045) ID, RP(ID), ZP(ID), RADDEG*THETAS(ID),
c     >     RADDEG*THETAS2(ID), DDS(ID), DDS2(ID), COSTET(ID)
c      ENDDO
c
c      DO IR = IRSEP,IRWALL
c
c        IK = 1
c
c        ID = ID + 1
c       IF ( CTARGOPT.EQ.4 .OR. CTARGOPT.EQ.5 )THEN
c          KP = KORPG(IK+1,IR)
c        ELSE
c          KP = KORPG(IK,IR)
c        ENDIF
c
c        if (kp.ne.0.and.nvertp(kp).gt.0) then 
c            
c           DELTAR = RVERTP(2,KP) - RVERTP(1,KP)
c           DELTAZ = ZVERTP(2,KP) - ZVERTP(1,KP)
c           DDS2(ID) = SQRT(DELTAR**2 + DELTAZ**2)
c           THETAS2(ID) = ATAN2C(DELTAZ,DELTAR) - 0.5 * PI
c           IF (THETAS2(ID).LT.-PI) THETAS2(ID) = THETAS2(ID) + 2.0 * PI
C	   
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c             DELTAR = RS(IK,IR) - RS(IK+1,IR)
c             DELTAZ = ZS(IK,IR) - ZS(IK+1,IR)
c           ELSE
c             DELTAR = RP(ID) - RS(IK,IR)
c             DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c	   
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c           IF( CTARGOPT.EQ.1 )THEN
c             COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c             COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
c
c        else
c   
c           DDS2(ID)    = DDS(ID)
c           THETAS2(ID) = THETAS(ID)
c           IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c             DELTAR = RS(IK,IR) - RS(IK+1,IR)
c             DELTAZ = ZS(IK,IR) - ZS(IK+1,IR)
c           ELSE
c             DELTAR = RP(ID) - RS(IK,IR)
c             DELTAZ = ZP(ID) - ZS(IK,IR)
c           ENDIF
c           XCOS = ATAN2C(DELTAZ,DELTAR) - PI
c           IF( CTARGOPT.EQ.1 )THEN
c             COSTET(ID) = COS( XCOS - THETAS2(ID) )
c           ELSE
c             COSTET(ID) = COS( XCOS - THETAS(ID) )
c           ENDIF
c      
c        endif
c
C
c        WRITE(6,9045) ID, RP(ID), ZP(ID), RADDEG*THETAS(ID),
c     >     RADDEG*THETAS2(ID), DDS(ID), DDS2(ID), COSTET(ID)
c      ENDDO
c

c
C
C  FOR NON-ORTHOGONAL GEOMETRIES RECORD THETA VALUES AT TARGET PLATES
C  FOR CTARGOPT = 0 OR 4 (ALREADY DONE FOR CTARGOPT = 1 AND 6 IN SUBROUTINE
C  TAU)
C
C  NOTE: These values cannot be calculated for CTARGOPT = 2, 3, or 5
C        as the underlying orthogonal grid must be known to enable THETA
C        to be found by internal extrapolation.  For these options
C        THETAT must be specified externally along with the target plate
C        coordinates.
C
c      IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
c        DO ID = 1 , NDSIN
c          THETAT(ID) = THETAG(1,IRDS(ID))
c        ENDDO
c        DO ID = NDSIN+1 , NDS
c          THETAT(ID) = THETAG(NKS(IRDS(ID)),IRDS(ID))
c        ENDDO
c      ENDIF
c
c     The above code appears to have a bug in that the ik,ir indices
c     which correspond to id = 1,ndsin are nks(ir),ir and not 1,ir as
c     were being assigned above. Changed to code below.
c
c slmod begin
      IF (nbr.GT.0.OR.cgridopt.EQ.LINEAR_GRID.or.
     >    cgridopt.eq.RIBBON_GRID) THEN
c...    More reliable method of calculating COSTET for generalized
c       grids:

        costet2 = 0.0

        DO ir = irsep, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          ik = 1
          id = korpg(ik,ir)
          in = idds(ir,2)
          a1 = rp(in)
          a2 = zp(in)
          b1 = rs(ik,ir)
          b2 = zs(ik,ir)
          DO ii = 1, 2
            c1 = rvertp(ii,id)
            c2 = zvertp(ii,id)
            IF (CalcPoint(a1,a2,b1,b2,c1,c2,t(ii)).LT.0) 
     .        CALL ER('DOTARG','Low IK intersection not found',*99)
          ENDDO
          IF (t(1).GT.t(2)) THEN
            c1 = rvertp(1,id)
            c2 = zvertp(1,id)
          ELSE
            t(1) = t(2)
            c1 = rvertp(2,id)
            c2 = zvertp(2,id)
          ENDIF
          b1 = a1 + t(1) * (b1 - a1)
          b2 = a2 + t(1) * (b2 - a2)
          shyp = SQRT((a1 - c1)**2 + (a2 - c2)**2)
          sadj = SQRT((b1 - c1)**2 + (b2 - c2)**2)
          costet2(in) = sadj / shyp

          WRITE(SLOUT,'(5X,A,5I5,2F10.4,4F12.6)')
     .      'II,IN,ID,IK,IR = ',
     .      ii,in,id,ik,ir,costet(in),costet2(in),(t(ii),ii=1,4)

          ik = nks(ir)
          id = korpg(ik,ir)
          in = idds(ir,1)
          a1 = rp(in)
          a2 = zp(in)
          b1 = rs(ik,ir)
          b2 = zs(ik,ir)
          DO ii = 3, 4
            c1 = rvertp(ii,id)
            c2 = zvertp(ii,id)
            IF (CalcPoint(a1,a2,b1,b2,c1,c2,t(ii)).LT.0) 
     .        CALL ER('DOTARG','High IK intersection not found',*99)
          ENDDO
          IF (t(3).GT.t(4)) THEN
            c1 = rvertp(3,id)
            c2 = zvertp(3,id)
          ELSE
            t(3) = t(4)
            c1 = rvertp(4,id)
            c2 = zvertp(4,id)
          ENDIF
          b1 = a1 + t(3) * (b1 - a1)
          b2 = a2 + t(3) * (b2 - a2)
          shyp = SQRT((a1 - c1)**2 + (a2 - c2)**2)
          sadj = SQRT((b1 - c1)**2 + (b2 - c2)**2)
          costet2(in) = sadj / shyp

          WRITE(SLOUT,'(5X,A,5I5,2F10.4,4F12.6)')
     .      'II,IN,ID,IK,IR = ',
     .      ii,in,id,ik,ir,costet(in),costet2(in),(t(ii),ii=1,4)
        ENDDO        


        discrepancy = .FALSE.
        DO in = 1, nds
          IF (idring(irds(in)).EQ.BOUNDARY) CYCLE
          IF (ABS((costet(in)-costet2(in))/costet2(in)).GT.0.01) THEN
            discrepancy = .TRUE.
            WRITE(SLOUT,*) 'WARNING: DISCREPANCY BETWEEN COSTET '//
     .                     'AND COSTET2 FOR TARGET IN=',in
          ENDIF
          write(50,'(A,3I5,2F12.5)') 'costet ',in,ikds(in),irds(in),
     .                               costet(in),costet2(in)

c...      Overwrite COSTET with COSTET2, and also THETAS for now:
          costet(in) = costet2(in)
          thetas(in) = thetas2(in)
        ENDDO

        IF (sloutput) THEN
          WRITE(0,*) 'THETAS and COSTET set from THETAS2 and COSTET2'
          IF (discrepancy) 
     .      CALL WN('DoTarg','Discrepancy found between COSTET and '//
     .              'COSTET2')
        ENDIF
      ENDIF
c slmod end
      if (northopt.eq.1.or.northopt.eq.3) then
c
        IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4 )THEN
          DO ID = 1 , NDSIN
            THETAT(ID) = THETAG(nks(irds(id)),IRDS(ID))
          ENDDO
          DO ID = NDSIN+1 , NDS
            THETAT(ID) = THETAG(1,IRDS(ID))
          ENDDO
        ENDIF
c
c
c       Karl found the following code necessary for some reason - however
c       - not sure why - include but comment out for now. 
c      
c
c       to prevent discontinuities in the flux to the target
c       caused by inaccuracies in the calculation of the incl. angle,
c       we need to average over jumps
c       this is an ugly kludge so we warn the user by printing
c       message directly to screen output
c       Krieger IPP/97
c
c        IF( CTARGOPT.EQ.6)THEN
c          DO ID = 3 , NDSIN-2
c            if (abs((thetas(id)-thetas(id-1))/thetas(id)).gt.0.1) then
c              thetas(id) = (thetas(id-1)+thetas(id+1))/2.
c              thetas2(id) = (thetas2(id-1)+thetas2(id+1))/2.
c              costet(id) = (costet(id-1)+costet(id+1))/2.
c              write(0,*) 'Warning: target angle corrected for id=',id
c              write(6,*) 'Warning: angle corrected; new values:'
c              write(6,*) 'id, thetas, thetas2, costet: ',id,
c     >                   raddeg*thetas(id),raddeg*thetas2(id),
c     >                   costet(id)
c            endif
c          ENDDO
c          DO ID = NDSIN+3 , NDS-2
c            if (abs((thetas(id)-thetas(id-1))/thetas(id)).gt.0.1) then
c              thetas(id) = (thetas(id-1)+thetas(id+1))/2.
c              thetas2(id) = (thetas2(id-1)+thetas2(id+1))/2.
c              costet(id) = (costet(id-1)+costet(id+1))/2.
c              write(0,*) 'Warning: target angle corrected for id=',id
c              write(6,*) 'Warning: angle corrected; new values:'
c              write(6,*) 'id, thetas, thetas2, costet: ',id,
c     >                   raddeg*thetas(id),raddeg*thetas2(id),
c     >                   costet(id)
c            endif
c          ENDDO
c        ENDIF
c
c
c
      endif
c nonorth
c
c       Endif for northopt
c
c
c     ASDEX and ITER grids
c
      elseif (cgridopt.eq.1.or.cgridopt.eq.2) then
         do id = 1,nds
            dds2(id) = dds(id)
            thetas2(id) = thetas(id)
         enddo
      endif
c
c     Calculate the distance of each target point from the respective
c     separatrix target point.
c
      write (6,*) 'Calculating SEPDIST:'
c
      do ir = irsep,nrs
c
c        Do INNER targets  (Index 1)
c
         sepdist(idds(ir,1)) = sqrt(
     >                   (rp(idds(ir,1))-rp(idds(irsep,1)))**2
     >                +  (zp(idds(ir,1))-zp(idds(irsep,1)))**2)
c
c        Do OUTER targets  (Index 2)
c
         sepdist(idds(ir,2)) = sqrt(
     >                    (rp(idds(ir,2))-rp(idds(irsep,2)))**2
     >                 +  (zp(idds(ir,2))-zp(idds(irsep,2)))**2)
c
c        Print out Sepdist
c
         if (cprint.eq.3.or.cprint.eq.9) then
c
            write (6,'(''IN :'',3i4,g16.8,2f16.8)')
     >               ir,ikds(idds(ir,1)),idds(ir,1),
     >               sepdist(idds(ir,1)),rp(idds(ir,1)),zp(idds(ir,1))
            write (6,'(''OUT:'',3i4,g16.8,2f16.8)')
     >               ir,ikds(idds(ir,2)),idds(ir,2),
     >               sepdist(idds(ir,2)),rp(idds(ir,2)),zp(idds(ir,2))
         endif
      end do
c
c     Calculate distances ALONG the targets
c
      Write(6,*) 'Calculating Distance along target:'
c
      isepdist = 0.0
      osepdist = 0.0
c
      do ir = irsep,irwall
c
        isepdist = isepdist + dds(idds(ir,1)) /2.0
        osepdist = osepdist + dds(idds(ir,2)) /2.0

        sepdist2(idds(ir,1)) = isepdist
        sepdist2(idds(ir,2)) = osepdist

        isepdist = isepdist + dds(idds(ir,1)) /2.0
        osepdist = osepdist + dds(idds(ir,2)) /2.0


        if (cprint.eq.3.or.cprint.eq.9) then
c
            write (6,'(a,3i4,g16.8,2f16.8)') INNER//':',
     >               ir,ikds(idds(ir,1)),idds(ir,1),
     >               sepdist2(idds(ir,1))
            write (6,'(a,3i4,g16.8,2f16.8)') OUTER//':',
     >               ir,ikds(idds(ir,2)),idds(ir,2),
     >               sepdist2(idds(ir,2))
c
        endif
c
      end do
c
      isepdist = 0.0
      osepdist = 0.0
c
      do ir = nrs,irtrap,-1
c
        isepdist = isepdist - dds(idds(ir,1)) /2.0
        osepdist = osepdist - dds(idds(ir,2)) /2.0

        sepdist2(idds(ir,1)) = isepdist
        sepdist2(idds(ir,2)) = osepdist

        isepdist = isepdist - dds(idds(ir,1)) /2.0
        osepdist = osepdist - dds(idds(ir,2)) /2.0


        if (cprint.eq.3.or.cprint.eq.9) then
c
            write (6,'(a,3i4,g16.8,2f16.8)') INNER//':',
     >               ir,ikds(idds(ir,1)),idds(ir,1),
     >               sepdist2(idds(ir,1))
            write (6,'(a,3i4,g16.8,2f16.8)') OUTER//':',
     >               ir,ikds(idds(ir,2)),idds(ir,2),
     >               sepdist2(idds(ir,2))
c
        endif
c
      end do
c
c     Calculate R and Z distances from the R,Z location of
c     the strike points.
c 
c     Ik = 1 target 
c
      kp = korpg(1,irsep)
      rspo = rvertp(1,kp) 
      zspo = zvertp(1,kp)
      kp = korpg(nks(irsep),irsep)
      rspi = rvertp(4,kp) 
      zspi = zvertp(4,kp)
c
c     Loop across targets - "Inner"
c 
c     and Write out results - R-RSEP and Z-ZSEP values
c
      if (cprint.eq.3.or.cprint.eq.9) then 
         write(6,'(a)') 'R-Rsep and Z-Zsep for target elements'  
c 
         write(6,'(a,2(1x,g12.5))') INNER//' Strike Point:',
     >                           rspi,zspi  
      endif 
c
      do id = 1,ndsin
c
         rspdist(id) = rp(id) - rspi
         zspdist(id) = zp(id) - zspi
c
         if (cprint.eq.3.or.cprint.eq.9) then 
            write(6,'(a,3i5,2(1x,g12.5))') 'ELE:',id,
     >            ikds(id),irds(id),rspdist(id),zspdist(id)  
         endif
c
      end do 
c
c     "Outer"
c
      if (cprint.eq.3.or.cprint.eq.9) then 
         write(6,'(a,2(1x,g12.5))') OUTER//' Strike Point:',
     >                           rspo,zspo  
      endif
c
      do id = ndsin+1,nds
c
         rspdist(id) = rp(id) - rspo
         zspdist(id) = zp(id) - zspo
c
         if (cprint.eq.3.or.cprint.eq.9) then  
            write(6,'(a,3i5,2(1x,g12.5))') 'ELE:',id,
     >            ikds(id),irds(id),rspdist(id),zspdist(id)  
         endif   
c
      end do 
c
c
      



c
c
 9044 FORMAT(1X,I8,2F10.5,2X,2F8.2,2F10.7)
 9045 FORMAT(1X,I8,2F10.5,2X,2F8.2,2F10.7,F8.5)
 9046 FORMAT(7X,'ID',4x,'RP',7x,'ZP',6X,
     >       'THETAS',2x,'THETAS2',
     >       4X,'DDS',6X,'DDS2',5x,'COSTET')
C
      return
c slmod begin
c...  For COSTET2 calculation:
99    CONTINUE
      WRITE(EROUT,'(5X,A,5I5)')
     .  'II,IN,ID,IK,IR = ',ii,in,id,ik,ir
      WRITE(EROUT,'(5X,A,4E15.7)') 'T:',(t(ii),ii=1,4)
      DO in = 1, nds
        WRITE(EROUT,'(5X,A,I5,2F12.6)')
     .    'IN,COSTET,COSTET2 = ',in,costet(in),costet2(in)
      ENDDO
      CALL OutputData(85,'Error A in DOTARG')
      STOP
c slmod end
      END
C
C
C
      SUBROUTINE DOWALL
      implicit none
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "CGEOM"
      include 'cgeom'
C     INCLUDE "CNEUT"
c      include 'cneut'
C     INCLUDE "COMTOR"
      include 'comtor'
C     INCLUDE "CIONIZ"
      include 'cioniz'
C     INCLUDE "READER"
      include 'reader'
C     INCLUDE "DYNAM5"
      include 'dynam5'
      INTEGER IK,IR,K,NP,L,J,I,NR,NC,NXW,IEXTRA,JK,JR,MIZS,IZ,IERR,ID
      INTEGER IX,IY,IKIN,IKOUT,IRIN,IROUT,MKS,NP1,ICOUNT,INEXT,KNEXT
      INTEGER IND,IW,in
      INTEGER STARTID,ENDID,STARTIK,ENDIK,stepik,IKOUT2,IROUT2
      REAL    R1,Z1,R2,Z2,R3,Z3,RWT,ZWT,SRCPRB,R4,Z4
      REAL    DD1,DD2
      INTEGER RES
      INTEGER FINDNUM
      REAL    ATAN2C
      EXTERNAL FINDNUM,ATAN2C
c
      integer ikn,irn
      real dist,ri,zi
C
C     WALLPT (IND,1) = R
C     WALLPT (IND,2) = Z
C     WALLPT (IND,3) = WEIGHT FACTOR FOR ANTI-CLOCKWISE
C     WALLPT (IND,4) = WEIGHT FACTOR FOR CLOCKWISE
C     WALLPT (IND,5) = LENGTH OF 1/2 SEGMENT ANTI-CLOCKWISE
C     WALLPT (IND,6) = LENGTH OF 1/2 SEGMENT CLOCKWISE
C     WALLPT (IND,7) = TOTAL LENGTH OF LAUNCH SEGMENT
C     WALLPT (IND,8) = ANGLE FOR ANTI-CLOCKWISE LAUNCH
C     WALLPT (IND,9) = ANGLE FOR CLOCKWISE LAUNCH
C     WALLPT (IND,10) = NET PROBABILITY ANTI-CLOCKWISE
C     WALLPT (IND,11) = NET PROBABILITY CLOCKWISE
C     WALLPT (IND,12) = NET PROBABILITY FOR ENTIRE SEGMENT
C     WALLPT (IND,13) = FINAL PROBABILITY FOR SEGMENT
c
c     wallpt (ind,16) = TYPE OF WALL SEGMENT
c                       1 = Outer Target (JET) - inner for Xpt down
c                       4 = Inner Target (JET) - outer      "
c                       7 = Main Wall
c                       8 = Private Plasma Wall
c
c                       9 = Baffle Segment
c
c                       These are similar to the quantity in the JVESM
c                       array associated with the NIMBUS wall
c                       specification. The difference is that the
c                       Main Wall is split into Inner and Outer Divertor
c                       Wall as well as the Main (SOL) Wall - this
c                       is not done here.
c
c     WALLPT (ind,17) = INDEX into the NIMBUS flux data returned
c                       for each wall segment - ONLY if the NIMBUS
c                       wall option has been specified. NOTE: if
c                       the NIMBUS wall has been specified - it is
c                       still combined with the DIVIMP target polygon
c                       corners because rounding errors may result in
c                       small discrepancies between the coordinates.
c
c     WALLPT (IND,18) = Index of corresponding target segment if the wall
c                       segment is also a target segment.
c
c     WALLPT (IND,19) = Temperature of wall segment in Kelvin (K)
c
c     WALLPT (IND,20) = RSTART
c     WALLPT (IND,21) = ZSTART
c     WALLPT (IND,22) = REND
c     WALLPT (IND,23) = ZEND
c
c     wallpt (ind,24) = Used for additional indexing information - used
c                       as IK knot number for wall and trap wall option 7
c
c     wallpt (ind,25) = Value of reflection coefficient - if reflection
c                       for this segment is turned off the value here
c                       will be zero. If a positive value is specified
c                       then regular reflection occurs. If it is negative
c                       then a PTR (prompt thermal re-emission) type
c                       reflection is used. The value for this is
c                       set with the individual YMF's and is read from
c                       the CYMFS array.
c
c     wallpt (ind,26) = IK value of nearest plasma cell to wall segment
c     wallpt (ind,27) = IR value of nearest plasma cell to wall segment
c     wallpt (ind,28) = Minimum distance to outermost ring
c     wallpt (ind,29) = Plasma Te at wall segment - Temporary storage for RI
c     wallpt (ind,30) = Plasma Ti at wall segment - Temporary storage for ZI
c     wallpt (ind,31) = Plasma density at wall segment
c
c slmod begin
c     NOTE: ANY CHANGES MADE HERE MUST ALSO BE DONE IN BUILDNEUTRALWALL, THE
c     ROUTINE THAT ASSEMBLES THE VESSEL WALL DESCRIPTION FOR GENERALIZED
c     GRID GEOMETRIES. -SL, 24.2.2004
c slmod end
C
C     RW (PCNT) = R COORDINATES OF WALL BOUNDARY
C     ZW (PCNT) = Z COORDINATES OF WALL BOUNDARY
C
C     RW AND ZW CONTAIN THE "TURNING POINTS" OF THE WALL DEFINITION
C     FOR EXAMPLE, WHEN THE WALL IS SPECIFIED BY AN INPUT SET OF
C     DATA OR BY USING THE PRE-LOADED GEOMETRY DATA ... THE
C     INFORMATION SUPPLIED IS THE COORDINATES OF THE END POINTS
C     OF THE WALL SEGMENTS. THE CODE KEEPS TRACK OF THE MIDDLE OF
C     THE WALL SEGMENTS ... SINCE IT IS THE MID-POINTS OF THE BINS
C     THAT ARE USED FOR PARTICLE TRACKING. HOWEVER, IF THE MIDPOINTS
C     OF THE INPUT DATA ARE USED FOR THE WALL DEFINITION THEN SOME
C     SPACE WILL BE LEFT OUT ... THAT SHOULD BE INSIDE THE WALL.
C
C     THE RW,ZW VALUES DO NOT NECESSARILY MATCH THE WALLPT(*,1 AND 2)
C     VALUES SINCE THERE ARE SOME CORNERS WHICH MUST BE DEFINED TO
C     DEFINE THE "WALL EDGE" BUT WHICH ARE NOT REQUIRED TO DEFINE
C     THE PARTICLE LAUNCH.
C
C     THE FUNCTION ATAN2C - IS AN ERROR CHECKING FRONT END TO
C     THE ATAN2 FUNCTION. FOR EITHER DELTAZ OR DELTAR = 0 IT RETURNS
C     THE APPROPRIATE VALUE OF THE ANGLE - E.G. PI,PI/2,0,-PI/2
C     IN THE CASE WHERE BOTH ARGUMENTS ARE ZERO - IT CURRENTLY
C     RETURNS A ZERO VALUE.
C
C-----------------------------------------------------------------------
C     CALCULATE THE WALLS AND THEIR NORMALS ...
C-----------------------------------------------------------------------
C
c
c     Initialization
c
      call rzero(wallpt,maxpts*31)
      call izero(wallindex,maxnds)
c
c     After the targets have been calculated but before doing the walls-
c     if Neutral Wall option 4 (Vessel coordinates taken from grid file)
c     CNEUR=4 or CTRAP=4 it is necessary to map the RVES,ZVES arrays
c     onto the WALLCO and WALLCO2 arrays respectively.
c
      if (cneur.eq.4.or.ctrap.eq.4) then
         call dovessel
      endif
c
c
c     If CNEUR = 2 or 3 or 4 OR CTRAP is 3 or 4
c     and the target option (CTARGOPT) is 6
c     Add the corner points of the end polygons to the neutral wall
c     and privarte plasma wall specifications.
c
      if ((cneur.eq.2.or.cneur.eq.3.or.cneur.eq.4)
     >    .and.ctargopt.eq.6) then
         call fixwallco
      endif
c
      if ((ctrap.eq.3.or.ctrap.eq.4)
     >    .and.ctargopt.eq.6) then
         call fixwallco2
      endif
c
c
c     Check for repeated points in the wallco arrays and remove them
c
      if (nwall.gt.0) then
         call check_wallco(wallco,nwall)
      endif
c
      if (nwall2.gt.0) then
         call check_wallco(wallco2,nwall2)
      endif
c
c     Print out walls if option is set
c
      if (cprint.eq.3.or.cprint.eq.9) then
c
         write (6,*) 'WallcoA:',nwall
         do in = 1,nwall
            write (6,'(2(1x,g18.10),i4)')
     >        wallco(in,1),wallco(in,2),in
         end do
c
         write (6,*) 'Wallco2A:',nwall2
         do in = 1,nwall2
            write (6,'(2(1x,g18.10),i4)')
     >        wallco2(in,1),wallco2(in,2),in
         end do
c
      endif
c
c     Start processing the wall:
c
      WRITE(6,*) 'DO WALL:'
      IND = 1
      PCNT = 1
C
C     SET UP INITIAL VALUES AND LOOP BOUNDARIES DEPENDING ON OPTIONS
C
      IF (CNEUR.EQ.0) THEN
         if (ctargopt.eq.6) then
            in = korpg(ikds(nds-1),irds(nds-1))
            R1 = RVERTP(2,in)
            Z1 = ZVERTP(2,in)
         else
            R1 = 0.5 * (RP(NDS)+RP(NDS-1))
            Z1 = 0.5 * (ZP(NDS)+ZP(NDS-1))
         endif
         IR = IRWALL-1
      ELSEIF (CNEUR.EQ.1) THEN
         R1 = RP(NDS)
         Z1 = ZP(NDS)
         IR = IRWALL
      ELSEIF (CNEUR.EQ.2.OR.CNEUR.EQ.3.or.cneur.eq.4) THEN
         R1 = WALLCO(1,1)
         Z1 = WALLCO(1,2)
      elseif (cneur.eq.7) then
         in = korpg(ikds(nds-1),irds(nds-1))
         R1 = RVERTP(2,in)
         Z1 = ZVERTP(2,in)
         ir = irwall -1
      ENDIF
C
C
      IF (CNEUR.EQ.2.OR.CNEUR.EQ.3.or.cneur.eq.4) THEN
         STARTIK = 1
         ENDIK = NWALL-1
      ELSE
         STARTIK = 1
         ENDIK = NKS(IR)
      ENDIF
C
C
      RW(PCNT) = R1
      ZW(PCNT) = Z1
      PCNT = PCNT +1
      WLWALL1 = 1
C
C       THIS POINT IS A TURNING POINT FOR THE WALL BUT IS NOT A LAUNCH
C       POINT. THE PREVIOUS TARGET POINT (NDS-1) WILL HAVE THE PLATE
C       PORTION OF THIS SEGMENT AND THE FIRST WALL POINT WILL HAVE THE
C       WALL SEGMENT. IF TARGET OPTION 0 HAS BEEN SELECTED THEN R1,Z1
C       WILL BE THE SAME AS R2,Z2. WHICH WILL RESULT IN A ZERO LENGTH
C       ANTI-CLOCKWISE WALL LAUNCH SEGMENT FOR THE FIRST POINT.
C
      DO 2035 IK = STARTIK,ENDIK
        IF (CNEUR.EQ.0) THEN
           IKOUT = IKOUTS(IK,IR)
           IROUT = IROUTS(IK,IR)
           R2 = RS(IK,IR)+0.50*(RS(IKOUT,IROUT)-RS(IK,IR))
           Z2 = ZS(IK,IR)+0.50*(ZS(IKOUT,IROUT)-ZS(IK,IR))
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKOUTS(IK+1,IR)
              IROUT2 = IROUTS(IK+1,IR)
              R4 = RS(IK+1,IR)+0.50*(RS(IKOUT2,IROUT2)-RS(IK+1,IR))
              Z4 = ZS(IK+1,IR)+0.50*(ZS(IKOUT2,IROUT2)-ZS(IK+1,IR))
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              R3 = 0.5 * (RP(1) + RP(2))
              Z3 = 0.5 * (ZP(1) + ZP(2))
           ENDIF
           RWT = R2
           ZWT = Z2
        ELSEIF (CNEUR.EQ.1) THEN
           R2 = RS(IK,IR)
           Z2 = ZS(IK,IR)
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKOUTS(IK+1,IR)
              IROUT2 = IROUTS(IK+1,IR)
              R4 = RS(IK+1,IR)
              Z4 = ZS(IK+1,IR)
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              R3 = RP(1)
              Z3 = ZP(1)
           ENDIF
           RWT = R2
           ZWT = Z2
         ELSEIF (CNEUR.EQ.2.OR.CNEUR.EQ.3.or.cneur.eq.4) THEN
           R3 = WALLCO(IK+1,1)
           Z3 = WALLCO(IK+1,2)
           R2 = (R1+R3)/2.0
           Z2 = (Z1+Z3)/2.0
           RWT = R3
           ZWT = Z3
         elseif (cneur.eq.7) then

           in = korpg(ik,ir)
           R3 = RVERTP(3,in)
           Z3 = ZVERTP(3,in)

           R2 = (R1+R3)/2.0
           Z2 = (Z1+Z3)/2.0
           RWT = R3
           ZWT = Z3
c
c          Set index for now - to be ik'th wall element
c
           wallpt(ind,24) = ik
c
c          Record nearest plasma cell 
c
           wallpt(ind,26) = ik 
           wallpt(ind,27) = ir
c
c          Calculate minimum distance to nearest plasma element 
c 

c
         ENDIF
C
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
         call find_nearest_boundary_cell(ikn,irn,r2,z2,dist,ri,zi)  
c
         if (cneur.ne.7) then 
c
c          Record nearest plasma cell 
c
           wallpt(ind,26) = ikn
           wallpt(ind,27) = irn
c
         endif
c
c        Save distance and intersection point - intersection is temporary
c
         wallpt(ind,28) = dist
         wallpt(ind,29) = ri
         wallpt(ind,30) = zi
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         wallpt(ind,16)= 7
         wallpt(ind,18)= 0
         IND = IND+1
         RW(PCNT) = RWT
         ZW(PCNT) = ZWT
         PCNT = PCNT +1
         R1 = R3
         Z1 = Z3
2035  CONTINUE
      WLWALL2 = IND-1
C
C     INNER PLATE
C
C
C     SET UP INITIAL VALUES AND LOOP BOUNDARIES DEPENDING ON OPTIONS
C
      IF (CNEUR.EQ.0) THEN
         if (ctargopt.eq.6) then
            in = korpg(ikds(2),irds(2))
            R1 = RVERTP(3,in)
            Z1 = ZVERTP(3,in)
         else
            R1 = 0.5 * (RP(1)+RP(2))
            Z1 = 0.5 * (ZP(1)+ZP(2))
         endif
         STARTID = 2
      ELSEIF (CNEUR.EQ.1) THEN
         R1 = RP(1)
         Z1 = ZP(1)
         STARTID = 1
      ELSEIF (CNEUR.EQ.2.OR.CNEUR.EQ.3.or.cneur.eq.4) THEN
         R1 = WALLCO(NWALL,1)
         Z1 = WALLCO(NWALL,2)
         STARTID = 2
      elseif (cneur.eq.7) then
         in = korpg(ikds(2),irds(2))
         R1 = RVERTP(3,in)
         Z1 = ZVERTP(3,in)
         STARTID = 2
      ENDIF


C
C
      IF (CTRAP.EQ.0.OR.CTRAP.EQ.2.or.ctrap.eq.3.or.ctrap.eq.4.or.
     >    ctrap.eq.7) THEN
         ENDID = NDSIN-1
      ELSEIF (CTRAP.EQ.1.or.ctrap.eq.8) THEN
         ENDID = NDSIN
      ENDIF
C
C
      RW(PCNT) = R1
      ZW(PCNT) = Z1
      PCNT = PCNT +1
C
      DO 2040 ID = STARTID,ENDID
         R2 = RP(ID)
         Z2 = ZP(ID)
         IF (ID.EQ.NDSIN) THEN
            R4 = R2
            Z4 = Z2
         ELSE
            R4 = RP(ID+1)
            Z4 = ZP(ID+1)
         ENDIF
c
         in = korpg(ikds(id),irds(id))
         if (ctargopt.eq.6.and.in.ne.0.and.nvertp(in).gt.0) then
c
c         if (ctargopt.eq.6.and.id.ne.1.and.id.ne.ndsin) then
c            in = korpg(ikds(id),irds(id))
c
            R3 = rvertp(4,in)
            Z3 = zvertp(4,in)
         else
            R3 = (R2 + R4) /2.0
            Z3 = (Z2 + Z4) /2.0
         endif
c
         wallindex(id) = ind
c
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
         wallpt(ind,26) = ikds(id)
         wallpt(ind,27) = irds(id)
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
c
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         wallpt(ind,16)= 4
         wallpt(ind,18)= id
         IND = IND+1
c
         RW(PCNT) = R2
         ZW(PCNT) = Z2
         PCNT = PCNT +1
         RW(PCNT) = R3
         ZW(PCNT) = Z3
         PCNT = PCNT +1
c
         R1 = R3
         Z1 = Z3
 2040 CONTINUE

c
c      RW(PCNT) = R3
c      ZW(PCNT) = Z3
c      PCNT = PCNT +1
C
c
c     The code for the trap wall section is only executed
c     when CTRAP.ne.8 - TRAP wall option 8 deletes the PFZ
c     wall section and is used for LIMITER grids.
c
C     TRAP
C
C
C     SET UP INITIAL VALUES AND LOOP BOUNDARIES DEPENDING ON OPTIONS
C
      WLTRAP1 = IND
      IF (CTRAP.EQ.0.or.ctrap.eq.2) THEN
         if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
            IR = IRTRAP+1
         elseif (cgridopt.eq.2) then
            ir = irtrap2+1
         endif
c
         if (ctargopt.eq.6) then
            in = korpg(ikds(ndsin-1),irds(ndsin-1))
            R1 = RVERTP(4,in)
            Z1 = ZVERTP(4,in)
         else
            R1 = 0.5 * (RP(NDSIN)+RP(NDSIN-1))
            Z1 = 0.5 * (ZP(NDSIN)+ZP(NDSIN-1))
         endif
c
      ELSEIF (CTRAP.EQ.1) THEN
         if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
            IR = IRTRAP
         elseif (cgridopt.eq.2) then
            ir = irtrap2
         endif
         R1 = RP(NDSIN)
         Z1 = ZP(NDSIN)
c
c      ELSEIF (CTRAP.EQ.2) THEN
c         R1 = 0.5 * (RP(NDSIN)+RP(NDSIN-1))
c         Z1 = 0.5 * (ZP(NDSIN)+ZP(NDSIN-1))
c
c         write(6,*) 'Test Trap Option 3, r1=',r1,'  z1=',z1
c
      elseif (ctrap.eq.3.or.ctrap.eq.4) then
         if (ctargopt.eq.6) then
            R1 = wallco2(1,1)
            Z1 = wallco2(1,2)
         else
            R1 = 0.5 * (RP(NDSIN)+RP(NDSIN-1))
            Z1 = 0.5 * (ZP(NDSIN)+ZP(NDSIN-1))
         endif
      elseif (ctrap.eq.7) then
c
         in = korpg(ikds(ndsin-1),irds(ndsin-1))
         R1 = RVERTP(4,in)
         Z1 = ZVERTP(4,in)
c
         IR = IRTRAP+1
c
      elseif (ctrap.eq.8) then
c
         in = korpg(ikds(ndsin),irds(ndsin))
         R1 = RVERTP(4,in)
         Z1 = ZVERTP(4,in)
c
      ENDIF
C
      RW(PCNT) = R1
      ZW(PCNT) = Z1
      PCNT = PCNT +1
C
      IF (CTRAP.EQ.2) THEN
         STARTIK = 1
         ENDIK = 1
         stepik = 1
      elseif (ctrap.eq.3.or.ctrap.eq.4) then
         if (ctargopt.eq.6) then
            startik= 2
            endik= nwall2
         else
            startik= 1
            endik= nwall2+1
         endif
         stepik = 1
c
      elseif (ctrap.eq.8) then 
c
c        Set loop values so that the TRAP loop does NOT execute
c 
         startik = 1
         endik   = 0
         stepik  = 1
c
      ELSE
         STARTIK = NKS(IR)
         ENDIK = 1
         stepik = -1
      ENDIF
C
C
      DO 2045 IK = STARTIK,ENDIK,STEPIK
        IF (CTRAP.EQ.0) THEN
           IKOUT = IKINS(IK,IR)
           IROUT = IRINS(IK,IR)
           R2 = RS(IK,IR)+0.50*(RS(IKOUT,IROUT)-RS(IK,IR))
           Z2 = ZS(IK,IR)+0.50*(ZS(IKOUT,IROUT)-ZS(IK,IR))
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKINS(IK-1,IR)
              IROUT2 = IRINS(IK-1,IR)
              R4 = RS(IK-1,IR)+0.50*(RS(IKOUT2,IROUT2)-RS(IK-1,IR))
              Z4 = ZS(IK-1,IR)+0.50*(ZS(IKOUT2,IROUT2)-ZS(IK-1,IR))
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              if (ctargopt.eq.6) then
                 in = korpg(ikds(ndsin+2),irds(ndsin+2))
                 R3 = RVERTP(1,in)
                 Z3 = ZVERTP(1,in)
              else
                 R3 = 0.5* (RP(NDSIN+1)+RP(NDSIN+2))
                 Z3 = 0.5* (ZP(NDSIN+1)+ZP(NDSIN+2))
              endif
           ENDIF
           RWT = R2
           ZWT = Z2
        ELSEIF (CTRAP.EQ.1) THEN
           R2 = RS(IK,IR)
           Z2 = ZS(IK,IR)
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKOUTS(IK-1,IR)
              IROUT2 = IROUTS(IK-1,IR)
              R4 = RS(IK-1,IR)
              Z4 = ZS(IK-1,IR)
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              R3 = RP(NDSIN+1)
              Z3 = ZP(NDSIN+1)
           ENDIF
           RWT = R2
           ZWT = Z2
        ELSEIF (CTRAP.EQ.2) THEN
c
           if (ctargopt.eq.6) then
              in = korpg(ikds(ndsin+2),irds(ndsin+2))
              R3 = RVERTP(1,in)
              Z3 = ZVERTP(1,in)
           else
              R3 = 0.5* (RP(NDSIN+1)+RP(NDSIN+2))
              Z3 = 0.5* (ZP(NDSIN+1)+ZP(NDSIN+2))
           endif
c
           R2 = (R1+R3)/2.0
           Z2 = (Z1+Z3)/2.0
           RWT = R3
           ZWT = Z3
        elseif (ctrap.eq.3.or.ctrap.eq.4) then
            if (ik.ne.nwall2+1) then
              r3=wallco2(ik,1)
              z3=wallco2(ik,2)
           else
              r3 = 0.5* (rp(ndsin+1)+rp(ndsin+2))
              z3 = 0.5* (zp(ndsin+1)+zp(ndsin+2))
           endif
c
c           write(6,*) 'Test trap option 3, r3=',r3,'  z3=',z3
c
           r2 = (r1+r3)/2.0
           z2 = (z1+z3)/2.0
           rwt = r3
           zwt = z3
        elseif (ctrap.eq.7) then
c
           in = korpg(ik,ir)
           R3 = RVERTP(1,in)
           Z3 = ZVERTP(1,in)
c
           R2 = (R1+R3)/2.0
           Z2 = (Z1+Z3)/2.0
           RWT = R3
           ZWT = Z3
c
c          Set index for now - to be ik'th wall element
c
           wallpt(ind,24) = ik
c
c          Record nearest plasma cell 
c
           wallpt(ind,26) = ik 
           wallpt(ind,27) = ir
c
        ENDIF
c
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
c
         call find_nearest_boundary_cell(ikn,irn,r2,z2,dist,ri,zi)  
c
         if (ctrap.ne.7) then 
c
c          Record nearest plasma cell 
c
           wallpt(ind,26) = ikn
           wallpt(ind,27) = irn
c
         endif
c
c        Save distance and intersection point - intersection is temporary
c
         wallpt(ind,28) = dist
         wallpt(ind,29) = ri
         wallpt(ind,30) = zi
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         wallpt(ind,16)= 8
         wallpt(ind,18)= 0
         IND = IND+1
c
         RW(PCNT) = RWT
         ZW(PCNT) = ZWT
         PCNT = PCNT +1
c
         R1 = R3
         Z1 = Z3
2045  CONTINUE

      WLTRAP2 = IND - 1
C
C
C         RW(PCNT) = R3
C         ZW(PCNT) = Z3
C         PCNT = PCNT +1
C
C        OUTER PLATE
C
C
C     SET UP INITIAL VALUES AND LOOP BOUNDARIES DEPENDING ON OPTIONS
C
      IF (CTRAP.EQ.0.OR.CTRAP.EQ.2.or.ctrap.eq.3.or.ctrap.eq.4) THEN
         if (ctargopt.eq.6) then
            in = korpg(ikds(ndsin+2),irds(ndsin+2))
            r1 = rvertp(1,in)
            z1 = zvertp(1,in)
         else
            R1 = 0.5* (RP(NDSIN+1)+RP(NDSIN+2))
            Z1 = 0.5* (ZP(NDSIN+1)+ZP(NDSIN+2))
         endif
         STARTID = NDSIN+2
      ELSEIF (CTRAP.EQ.1) THEN
         R1 = RP(NDSIN+1)
         Z1 = ZP(NDSIN+1)
         STARTID = NDSIN+1
      elseif (ctrap.eq.7) then
         in = korpg(ikds(ndsin+2),irds(ndsin+2))
         r1 = rvertp(1,in)
         z1 = zvertp(1,in)
         STARTID = NDSIN+2
      elseif (ctrap.eq.8) then
         in = korpg(ikds(ndsin+1),irds(ndsin+1))
         r1 = rvertp(1,in)
         z1 = zvertp(1,in)
         STARTID = NDSIN+1
      ENDIF
C
C
      IF (CNEUR.EQ.0.OR.CNEUR.EQ.2.OR.CNEUR.EQ.3
     >    .or.cneur.eq.4.or.cneur.eq.7) THEN
         if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
            ENDID = NDS-1
         elseif (cgridopt.eq.2) then
            endid = ndsin2-1
         endif
      ELSEIF (CNEUR.EQ.1) THEN
         if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
            ENDID = NDS
         elseif (cgridopt.eq.2) then
            endid = ndsin2
         endif
      ENDIF
C
C
      RW(PCNT) = R1
      ZW(PCNT) = Z1
      PCNT = PCNT +1
C
C
      DO 2050 ID = STARTID,ENDID
         R2 = RP(ID)
         Z2 = ZP(ID)
         IF (ID.EQ.NDS.or.(cgridopt.eq.2.and.id.eq.ndsin2))THEN
            R4 = R2
            Z4 = Z2
         ELSE
            R4 = RP(ID+1)
            Z4 = ZP(ID+1)
         ENDIF
c
         IF ((ID.EQ.ENDID).AND.(CNEUR.EQ.2.OR.CNEUR.EQ.3
     >                          .or.cneur.eq.4)) THEN
            if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
               R3 = WALLCO(1,1)
               Z3 = WALLCO(1,2)
            elseif (cgridopt.eq.2) then
               R3 = WALLCO2(1,1)
               Z3 = WALLCO2(1,2)
            endif
         ELSEIF (id.eq.endid.and.cneur.eq.7) then
            in = korpg(ikds(id),irds(id))
            r3 = rvertp(2,in)
            z3 = zvertp(2,in)
         elseif (id.eq.endid.and.cneur.eq.1) then 
            R3 = (R2 + R4) /2.0
            Z3 = (Z2 + Z4) /2.0
         else  
c
            if (ctargopt.eq.6) then
               in = korpg(ikds(id),irds(id))
c
               r3 = rvertp(2,in)
               z3 = zvertp(2,in)
            else
               R3 = (R2 + R4) /2.0
               Z3 = (Z2 + Z4) /2.0
            endif
         ENDIF
C
         wallindex(id) = ind
c
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
         wallpt(ind,26) = ikds(id)
         wallpt(ind,27) = irds(id)
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
c
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         wallpt(ind,16)= 1
         wallpt(ind,18)= id
         IND = IND+1
c
         RW(PCNT) = R2
         ZW(PCNT) = Z2
         PCNT = PCNT +1
c
         RW(PCNT) = R3
         ZW(PCNT) = Z3
         PCNT = PCNT +1
c
         R1 = R3
         Z1 = Z3
 2050 CONTINUE
c
c     For an ITER geometry - one now has to add the second half
c     of the wall to the wall definition. This will require an
c     additional input of a second half wall specification ... if the
c     wall is being specified in the input file.
c
c     This code has NOT been updated for CTARGOPT = 6
c
      if (cgridopt.eq.2) then
c
c
c     ITER grids
c
c
      IF (CNEUR.EQ.0) THEN
         R1 = 0.5 * (RP(NDS)+RP(NDS-1))
         Z1 = 0.5 * (ZP(NDS)+ZP(NDS-1))
         IR = IRWALL2-1
      ELSEIF (CNEUR.EQ.1) THEN
         R1 = RP(NDSIN2)
         Z1 = ZP(NDSIN2)
         IR = IRWALL2
      ELSEIF (CNEUR.EQ.2.OR.CNEUR.EQ.3) THEN
         R1 = WALLCO2(1,1)
         Z1 = WALLCO2(1,2)
      ENDIF
C
C
      IF (CNEUR.EQ.2.OR.CNEUR.EQ.3) THEN
         STARTIK = 1
         ENDIK = NWALL2-1
      ELSE
         STARTIK = 1
         ENDIK = NKS(IR)
      ENDIF
C
C
      RW(PCNT) = R1
      ZW(PCNT) = Z1
      PCNT = PCNT +1
      WLWALL3 = IND-1
C
C       THIS POINT IS A TURNING POINT FOR THE WALL BUT IS NOT A LAUNCH
C       POINT. THE PREVIOUS TARGET POINT (NDS-1) WILL HAVE THE PLATE
C       PORTION OF THIS SEGMENT AND THE FIRST WALL POINT WILL HAVE THE
C       WALL SEGMENT. IF TARGET OPTION 0 HAS BEEN SELECTED THEN R1,Z1
C       WILL BE THE SAME AS R2,Z2. WHICH WILL RESULT IN A ZERO LENGTH
C       ANTI-CLOCKWISE WALL LAUNCH SEGMENT FOR THE FIRST POINT.
C
      DO 5035 IK = STARTIK,ENDIK
        IF (CNEUR.EQ.0) THEN
           IKOUT = IKOUTS(IK,IR)
           IROUT = IROUTS(IK,IR)
           R2 = RS(IK,IR)+0.50*(RS(IKOUT,IROUT)-RS(IK,IR))
           Z2 = ZS(IK,IR)+0.50*(ZS(IKOUT,IROUT)-ZS(IK,IR))
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKOUTS(IK+1,IR)
              IROUT2 = IROUTS(IK+1,IR)
              R4 = RS(IK+1,IR)+0.50*(RS(IKOUT2,IROUT2)-RS(IK+1,IR))
              Z4 = ZS(IK+1,IR)+0.50*(ZS(IKOUT2,IROUT2)-ZS(IK+1,IR))
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              R3 = 0.5 * (RP(1) + RP(2))
              Z3 = 0.5 * (ZP(1) + ZP(2))
           ENDIF
           RWT = R2
           ZWT = Z2
        ELSEIF (CNEUR.EQ.1) THEN
           R2 = RS(IK,IR)
           Z2 = ZS(IK,IR)
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKOUTS(IK+1,IR)
              IROUT2 = IROUTS(IK+1,IR)
              R4 = RS(IK+1,IR)
              Z4 = ZS(IK+1,IR)
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              R3 = RP(1)
              Z3 = ZP(1)
           ENDIF
           RWT = R2
           ZWT = Z2
         ELSEIF (CNEUR.EQ.2.OR.CNEUR.EQ.3) THEN
           R3 = WALLCO2(IK+1,1)
           Z3 = WALLCO2(IK+1,2)
           R2 = (R1+R3)/2.0
           Z2 = (Z1+Z3)/2.0
           RWT = R3
           ZWT = Z3
         ENDIF
C
C
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         RW(PCNT) = RWT
         ZW(PCNT) = ZWT
c
c        Default wall launch weighting factor
c
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
c
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         IND = IND+1
         PCNT = PCNT +1
         R1 = R3
         Z1 = Z3
5035  CONTINUE
      WLWALL4 = IND-1
C
C     Upper Left plate
C
C
C     SET UP INITIAL VALUES AND LOOP BOUNDARIES DEPENDING ON OPTIONS
C
      IF (CNEUR.EQ.0) THEN
         R1 = 0.5 * (RP(ndsin2+1)+RP(ndsin2+2))
         Z1 = 0.5 * (ZP(ndsin2+1)+ZP(ndsin2+2))
         STARTID = ndsin2+2
      ELSEIF (CNEUR.EQ.1) THEN
         R1 = RP(ndsin2+1)
         Z1 = ZP(ndsin2+1)
         STARTID = ndsin2+1
      ELSEIF (CNEUR.EQ.2.OR.CNEUR.EQ.3) THEN
         R1 = WALLCO2(NWALL2,1)
         Z1 = WALLCO2(NWALL2,2)
         STARTID = ndsin2+2
      ENDIF
C
C
      IF (CTRAP.EQ.0.OR.CTRAP.EQ.2) THEN
         ENDID = NDSIN3-1
      ELSEIF (CTRAP.EQ.1) THEN
         ENDID = NDSIN3
      ENDIF
C
C
      RW(PCNT) = R1
      ZW(PCNT) = Z1
      PCNT = PCNT +1
C
C
      DO 5040 ID = STARTID,ENDID
         R2 = RP(ID)
         Z2 = ZP(ID)
         IF (ID.EQ.NDSIN3) THEN
            R4 = R2
            Z4 = Z2
         ELSE
            R4 = RP(ID+1)
            Z4 = ZP(ID+1)
         ENDIF
c
         R3 = (R2 + R4) /2.0
         Z3 = (Z2 + Z4) /2.0
c
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         RW(PCNT) = R2
         ZW(PCNT) = Z2
c
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
c
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         IND = IND+1
         PCNT = PCNT +1
         R1 = R3
         Z1 = Z3
 5040 CONTINUE
      RW(PCNT) = R3
      ZW(PCNT) = Z3
      PCNT = PCNT +1
C
C     TRAP part 2
C
C
C     SET UP INITIAL VALUES AND LOOP BOUNDARIES DEPENDING ON OPTIONS
C
      WLTRAP3 = IND
      IF (CTRAP.EQ.0) THEN
         IR = IRTRAP+1
         R1 = 0.5 * (RP(NDSIN3)+RP(NDSIN3-1))
         Z1 = 0.5 * (ZP(NDSIN3)+ZP(NDSIN3-1))
      ELSEIF (CTRAP.EQ.1) THEN
         IR = IRTRAP
         R1 = RP(NDSIN3)
         Z1 = ZP(NDSIN3)
      ELSEIF (CTRAP.EQ.2) THEN
         R1 = 0.5 * (RP(NDSIN3)+RP(NDSIN3-1))
         Z1 = 0.5 * (ZP(NDSIN3)+ZP(NDSIN3-1))
      ENDIF
C
C
      IF (CTRAP.EQ.2) THEN
         STARTIK = 1
         ENDIK = 1
      ELSE
         STARTIK = 1
         ENDIK = nks(ir)
      ENDIF
C
C
      DO 5045 IK = STARTIK,ENDIK
        IF (CTRAP.EQ.0) THEN
           IKOUT = IKINS(IK,IR)
           IROUT = IRINS(IK,IR)
           R2 = RS(IK,IR)+0.50*(RS(IKOUT,IROUT)-RS(IK,IR))
           Z2 = ZS(IK,IR)+0.50*(ZS(IKOUT,IROUT)-ZS(IK,IR))
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKINS(IK-1,IR)
              IROUT2 = IRINS(IK-1,IR)
              R4 = RS(IK-1,IR)+0.50*(RS(IKOUT2,IROUT2)-RS(IK-1,IR))
              Z4 = ZS(IK-1,IR)+0.50*(ZS(IKOUT2,IROUT2)-ZS(IK-1,IR))
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              R3 = 0.5* (RP(NDSIN3+1)+RP(NDSIN3+2))
              Z3 = 0.5* (ZP(NDSIN3+1)+ZP(NDSIN3+2))
           ENDIF
           RWT = R2
           ZWT = Z2
        ELSEIF (CTRAP.EQ.1) THEN
           R2 = RS(IK,IR)
           Z2 = ZS(IK,IR)
           IF (IK.NE.ENDIK) THEN
              IKOUT2 = IKOUTS(IK-1,IR)
              IROUT2 = IROUTS(IK-1,IR)
              R4 = RS(IK-1,IR)
              Z4 = ZS(IK-1,IR)
              R3 = (R2 + R4) /2.0
              Z3 = (Z2 + Z4) /2.0
           ELSE
              R3 = RP(NDSIN3+1)
              Z3 = ZP(NDSIN3+1)
           ENDIF
           RWT = R2
           ZWT = Z2
         ELSEIF (CTRAP.EQ.2) THEN
           R3 = 0.5* (RP(NDSIN3+1)+RP(NDSIN3+2))
           Z3 = 0.5* (ZP(NDSIN3+1)+ZP(NDSIN3+2))
           R2 = (R1+R3)/2.0
           Z2 = (Z1+Z3)/2.0
           RWT = R3
           ZWT = Z3
         ENDIF
c
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         RW(PCNT) = RWT
         ZW(PCNT) = ZWT
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         IND = IND+1
         PCNT = PCNT +1
         R1 = R3
         Z1 = Z3
5045     CONTINUE
         WLTRAP4 = IND - 1
C
C
C         RW(PCNT) = R3
C         ZW(PCNT) = Z3
C         PCNT = PCNT +1
C
C     Upper right target
C
C
C     SET UP INITIAL VALUES AND LOOP BOUNDARIES DEPENDING ON OPTIONS
C
      IF (CTRAP.EQ.0.OR.CTRAP.EQ.2) THEN
         R1 = 0.5* (RP(NDSIN3+1)+RP(NDSIN3+2))
         Z1 = 0.5* (ZP(NDSIN3+1)+ZP(NDSIN3+2))
         STARTID = NDSIN3+2
      ELSEIF (CTRAP.EQ.1) THEN
         R1 = RP(NDSIN3+1)
         Z1 = ZP(NDSIN3+1)
         STARTID = NDSIN3+1
      ENDIF
C
C
      IF (CNEUR.EQ.0.OR.CNEUR.EQ.2.OR.CNEUR.EQ.3) THEN
         ENDID = NDS-1
      ELSEIF (CNEUR.EQ.1) THEN
         ENDID = NDS
      ENDIF
C
C
      RW(PCNT) = R1
      ZW(PCNT) = Z1
      PCNT = PCNT +1
C
C
      DO 5050 ID = STARTID,ENDID
         R2 = RP(ID)
         Z2 = ZP(ID)
         IF (ID.EQ.NDS) THEN
            R4 = R2
            Z4 = Z2
         ELSE
            R4 = RP(ID+1)
            Z4 = ZP(ID+1)
         ENDIF
         IF ((ID.EQ.ENDID).AND.(CNEUR.EQ.2.OR.CNEUR.EQ.3)) THEN
            R3 = WALLCO(1,1)
            Z3 = WALLCO(1,2)
         ELSE
            R3 = (R2 + R4) /2.0
            Z3 = (Z2 + Z4) /2.0
         ENDIF
C
C
         WALLPT(IND,1) =  R2
         WALLPT(IND,2) =  Z2
c
         WALLPT(IND,20) = R1
         WALLPT(IND,21) = Z1
         WALLPT(IND,22) = R3
         WALLPT(IND,23) = Z3
c
         RW(PCNT) = R2
         ZW(PCNT) = Z2
c
         WALLPT(IND,3) = 1.0
         WALLPT(IND,4) = 1.0
c
         WALLPT(IND,5) = SQRT((R1-R2)**2+(Z1-Z2)**2)
         WALLPT(IND,6) = SQRT((R2-R3)**2+(Z2-Z3)**2)
         WALLPT(IND,7) = WALLPT(IND,5) + WALLPT(IND,6)
         WALLPT(IND,8) = ATAN2C(Z1-Z2,R1-R2)
         WALLPT(IND,9) = ATAN2C(Z3-Z2,R3-R2)
         IND = IND+1
         PCNT = PCNT +1
         R1 = R3
         Z1 = Z3
 5050 CONTINUE
      RW(PCNT) = R3
      ZW(PCNT) = Z3
      PCNT = PCNT +1
c
c     End additions for ITER grids
c
      endif
c
c     Finish off processing
c

      WALLPTS = IND-1
      RW(PCNT) = RW(1)
      ZW(PCNT) = ZW(1)
c
      call check_wall(rw,zw,pcnt)
C
      if (cprint.eq.3.or.cprint.eq.9) then
c
      WRITE (6,*) 'WLWALL:',WLWALL1,WLWALL2,WLTRAP1,WLTRAP2,WALLPTS
c
c     added for AUG (control of wall options)
c
      write (6,'(1x,a)') '  Rwall      Zwall     Index  '//
     >         '(Wall segment end points and duplicate turning points)'
      do 5051 id=1,pcnt
        write (6,'(1x,2(f9.5,2x),2x,i4)') rw(id),zw(id),id
 5051 continue
      write (6,'(a)') 'Actual wall: (wallpt data)' 
      write (6,'(1x,a)') '  Rwall      Zwall      Index'//
     >                   '    Plasma-ik  Plasma-ir'//
     >                   ' (Wall segment center points used in code)'
      do 5052 id=1,wallpts
        write (6,'(1x,2(f9.5,2x),2x,3(2x,i5),4(1x,g18.10))') 
     >                    wallpt(id,1),wallpt(id,2),id,
     >                    int(wallpt(id,26)),int(wallpt(id,27)),
     >                    wallpt(id,20),wallpt(id,21),wallpt(id,22),
     >                    wallpt(id,23)
 5052 continue

c
c     Print wallindex array
c
      do id = 1,nds
         write (6,'(a,2i5)') 'Wall Index:',id,wallindex(id)
      end do


c
c     end for print option 3
c
      endif
c
c     Calculate the distance along the wall and store it in wallpt(in,32)
c
c     Opt - currently does nothing - will allow for different calculation
c                                    methods when required. 
c
      call calc_wall_length_coordinate(1)
c
C
C     THIS SHOULD COMPLETE THE WALL DEFINITION
C
c
      RETURN
      END
C
C
C
      SUBROUTINE CALCWP
c slmod begin
      USE mod_divimp
c slmod end
      implicit none
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "COMTOR"
      include 'comtor'
      INTEGER IK,IND
      REAL    SRCPRB
C
C-----------------------------------------------------------------------
C
C     CALCWP: THIS SUBROUTINE CALCULATES THE WALL LAUNCH PROBABILITIES
C     FOR EACH WALL SEGMENT. THE VALUES IN THE ARRAY WLPROB IN
C     CONJUNCTION WITH OPTION SELECTOR WLPABS ARE USED TO ASSIGN
C     THE PROBABILITY OF LAUNCH AT RESPECTIVE WALL SEGMENTS.
C
C     THIS IS SEPARATE FROM THE CODE THAT DEALS WITH THE REST OF THE
C     WALL BECAUSE IT MAY BE NECESSARY, DEPENDING ON OPTIONS, TO WAIT
C     AFTER PIN PROGRAM EXECUTION, IN ORDER TO OBTAIN THE WALL
C     WEIGHTING FUNCTION
C
C-----------------------------------------------------------------------
C
C     CALCULATE CUMULATIVE PROBABILITIES AND APPLY ANY INPUT
C     WEIGHT FUNCTIONS.
C
C-----------------------------------------------------------------------
C
C     THE DEFAULT PROBABILITY IS THE SIZE OF THE WALL SEGMENTS. THIS IS
C     RESET LATER IF OTHER WEIGHT OPTIONS ARE SPECIFIED.
C

      DO 2100 IK = 1,WALLPTS
         WALLPT(IK,10) = WALLPT(IK,3)* WALLPT(IK,5)
         WALLPT(IK,11) = WALLPT(IK,4)* WALLPT(IK,6)
         WALLPT(IK,12) = WALLPT(IK,10)+WALLPT(IK,11)
         WALLPT(IK,13) = WALLPT(IK,12)
 2100 CONTINUE
C
C-----------------------------------------------------------------------
C
C     CALCULATE THE CUMULATIVE AND MODIFIED WALL LAUNCH PROBABILITIES
C     FOR EACH WALL SEGMENT. BASED ON THE LENGTH OF THE WALL SEGMENT.
C     THE PLATE PORTIONS OF THE WALL CAN HAVE A DEFAULT PROBABILITY OF
C     0.0 FOR WALL LAUNCH NEUTRALS - THIS CAN BE CHANGED IN THE
C     CODE ABOVE. THE BASE PROBABILITY OF LAUNCH FROM THE WALL
C     SEGMENT REPRESENTED BY A GIVEN WALL POINT IS IN THE
C     VARIABLE WALLPT(I,3:4) WHERE I IS THE POINT OF INTEREST. IT
C     DEFAULTS TO 1.0 FOR THE REST OF THE OUTER WALL AND TRAP WALL.
C
C-----------------------------------------------------------------------
C
C     THE ARRAY WLPROB CONTAINS THE MODIFICATIONS TO THE BASE
C     PROBABILITY IN THE FORM  (INDEX1)  (INDEX2)  (PROB MULTIPLIER)
C     THESE ARE APPLIED DIRECTLY TO THE WALLPT ARRAY.
C
C     THIS PROCESSING ONLY NEEDS TO BE DONE IF A WALL LAUNCH HAS
C     BEEN SPECIFIED.
C
C     IF WLPABS HAS BEEN SPECIFED AS 1 THEN THE INPUT WALL LAUNCH
C     PROBABILITIES ARE TAKEN AS ABSOLUTE VALUES AND OVERRIDE
C     THE WALL SEGMENT LENGTH MULTIPLICATION.
C
c slmod begin
C-----------------------------------------------------------------------
C     CALCULATE SPUTTERING YIELDS BASED ON PARTICLE FLUXES CALCULATED
C     IN OTHER DIVIMP RUNS
      IF (SPUTTER_NDATA.GT.0) CALL divCompileSputteringYields
c slmod end
C-----------------------------------------------------------------------
C
C     THE PROCESSING FOR THE WLPABS OPTIONS 2 AND 3 IS DONE
C     IN THE READPIN ROUTINE.
C     IF WLPABS HAS BEEN SPECIFIED AS 2 - THE QUANTITIES IN
C     WLPROB MULTIPLY THE WALL-SPUTTERING DATA RETURNED BY
C     PIN.
C     IF WLPABS HAS BEEN SPECIFIED AS 3 - THE VALUES FROM PIN
C     ARE USED AS DIRECTLY FOR WALL SEGMENT WEIGHTS.
C
C     DAVID ELDER   JAN 23  1992, SEPT 3 1992, JAN 18 1993
C
      IF (CNEUTB.EQ.2.OR.CNEUTH.EQ.2.or.cneutb.eq.4.or.cneuth.eq.4) THEN
        CALL RZERO (FWLPROB,MAXPTS)
        IF (NWLPROB.GT.0) THEN
          DO 2200 IND = 1, NWLPROB
            DO 2200 IK = INT(WLPROB(IND,1)),INT(WLPROB(IND,2))
              IF (WLPABS.EQ.0) THEN
                 WALLPT(IK,3) = WALLPT(IK,3) * WLPROB(IND,3)
                 WALLPT(IK,4) = WALLPT(IK,4) * WLPROB(IND,3)
                 WALLPT(IK,10) = WALLPT(IK,3)* WALLPT(IK,5)
                 WALLPT(IK,11) = WALLPT(IK,4)* WALLPT(IK,6)
                 WALLPT(IK,12) = WALLPT(IK,10)+WALLPT(IK,11)
                 WALLPT(IK,13) = WALLPT(IK,12)
              ELSEIF (WLPABS.EQ.1.OR.WLPABS.EQ.2.OR.WLPABS.EQ.3) THEN
                 WALLPT(IK,13) = WLPROB(IND,3)
              ENDIF
 2200     CONTINUE
        ENDIF
C
C-----------------------------------------------------------------------
C
C       GENERATE CUMULATIVE PROBABILITY FUNCTION AND CROSS REFERENCE
C       INDEX
C
C-----------------------------------------------------------------------
C
        IND = 1
        TOTWL = 0.0
        DO 2210 IK = 1,WALLPTS
          IF (WALLPT(IK,13).GT.0.0) THEN
            SRCPRB = WALLPT(IK,13)
            TOTWL = TOTWL + SRCPRB
            IF (IND.EQ.1) THEN
              FWLPROB(IND) = SRCPRB
            ELSE
              FWLPROB(IND) = FWLPROB(IND-1) + SRCPRB
            ENDIF
            WLIND(IND) = IK
            IND = IND +1
          ENDIF
2210    CONTINUE
        NWLIND = IND -1

        DO 2220 IND = 1,NWLIND
c
c           WRITE(6,*) 'WALL PROB:',IND,FWLPROB(IND),TOTWL,
c     >                  FWLPROB(IND)/TOTWL
c
           FWLPROB(IND) = FWLPROB(IND) / TOTWL
2220    CONTINUE

      ENDIF
C
C-----------------------------------------------------------------------
C
      if (cprint.eq.3.or.cprint.eq.9) then
         write (6,*) 'WALL DEFINITION:'
         write (6,*) 'IND, WLIND, FWLPROB:',nwlind
         do ik = 1,nwlind
            WRITE(6,*) 'WLIND:',ik,WLIND(IK),fwlprob(ik)
         end do
         WRITE(6,*) 'WALLPTS:', WALLPTS
         WRITE(6,*) 'R,Z,BP1,BP2,L1,L2,LT,A1,A2:'
         WRITE(6,*) 'PART1:'
         DO 2500 IND = 1,WALLPTS
            WRITE(6,'(a,i5,9g13.5)') 'P1:',ind,(WALLPT(IND,IK),IK=1,9)
            WRITE(6,'(a,5x,9g13.5)') 'P2:',(WALLPT(IND,IK),IK=10,18)
            WRITE(6,'(a,5x,9g13.5)') 'P3:',(WALLPT(IND,IK),IK=19,25)
 2500    CONTINUE
         WRITE(6,*) 'PART2:',PCNT
         DO 2520 IND = 1,PCNT
            WRITE(6,*) RW(IND),ZW(IND)
 2520    CONTINUE
      endif
C
C-----------------------------------------------------------------------
C
      RETURN
      END
C
C
C
      SUBROUTINE IONWALL
      IMPLICIT NONE
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "COMTOR"
      include 'comtor'
C     INCLUDE "CGEOM"
      include 'cgeom'
c
      include 'grbound'
      include 'divxy'
C
C     THIS ROUTINE CALCULATES THE UNDERLYING XY GRID
C     AND ITS CHARCTERISTICS. THIS INCLUDES DEFINING THE
C     POINTS ON THE IFXYS ARRAY WHICH LIE INSIDE THE
C     SPACE WHERE IONS CAN TRAVEL.
c
c     As a separate aspect of the above it also calculates
c     the effective boundary for ions which is used to
c     determine the location of individual X,Y points.
c     Although the X,Y points are rarely used the outer
c     wall for ions is still useful.
C
C     DAVID ELDER , SEPT 24, 1992
C
      REAL RESULT ,BEST,DSQ,R,Z
      INTEGER KIND
      INTEGER IK,IR,IND,ID,IKIN,IRIN,IKOUT,IROUT,IX,IY,JR,JK
c
      integer     IONWI1,IONWI2,IONTI1,IONTI2,
     >            ionwi3,ionwi4,ionti3,ionti4
c
      integer endid,in
C
      REAL RICHTABX,RICHTABY,RICHTCDX,RICHTCDY,MUE
c
c      real rcw(maxpts),zcw(maxpts)
c      integer ioncpts
C
c     IPP/08 Krieger - initialized some variables, which lead to
c     runtime error (use without being initialized in write statment)
c
      ionwi1 = 0
      ionwi2 = 0
      ionwi3 = 0
      ionwi4 = 0
      ionti1 = 0
      ionti2 = 0
      ionti3 = 0
      ionti4 = 0
c
      WRITE(6,*) 'DO ION WALL:'
      IF (CIONR.EQ.1) THEN
        IR = IRWALL
        IND = 1
        IONWI1 = 1
        DO 2000 IK = 1,NKS(IR)
           RIW(IND) = RS(IK,IR)
           ZIW(IND) = ZS(IK,IR)
           IND = IND+1
 2000   CONTINUE
        IONWI2 = IND -1
        DO 2005 ID = 1,NDSIN
           RIW(IND) = RP(ID)
           ZIW(IND) = ZP(ID)
           IND = IND +1
2005    CONTINUE
        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
           IR = IRTRAP
        elseif (cgridopt.eq.2) then
           ir = irtrap2
        endif
        IONTI1 = IND
        DO 2015 IK = NKS(IR),1,-1
           RIW(IND) = RS(IK,IR)
           ZIW(IND) = ZS(IK,IR)
           IND = IND+1
2015    CONTINUE
        IONTI2 = IND -2
        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
           endid = nds
        elseif (cgridopt.eq.2) then
           endid = ndsin2
        endif
        DO 2020 ID = NDSIN+1,endid
           RIW(IND) = RP(ID)
           ZIW(IND) = ZP(ID)
           IND = IND +1
2020    CONTINUE
c
c       Additional wall for ITER grids
c
        if (cgridopt.eq.2) then
           IR = IRWALL2
           IONWI3 = IND-1
           DO 3000 IK = 1,NKS(IR)
              RIW(IND) = RS(IK,IR)
              ZIW(IND) = ZS(IK,IR)
              IND = IND+1
3000       CONTINUE
           IONWI4 = IND -1
           DO 3005 ID = ndsin2+1,NDSIN3
              RIW(IND) = RP(ID)
              ZIW(IND) = ZP(ID)
              IND = IND +1
3005       CONTINUE
           IR = IRTRAP
           IONTI3 = IND
           DO 3015 IK = 1,nks(ir)
              RIW(IND) = RS(IK,IR)
              ZIW(IND) = ZS(IK,IR)
             IND = IND+1
3015       CONTINUE
           IONTI4 = IND -2
           DO 3020 ID = NDSIN3+1,nds
              RIW(IND) = RP(ID)
              ZIW(IND) = ZP(ID)
              IND = IND +1
3020       CONTINUE
c
c          End of extra for ITER grids
c
        endif
c
c
        RIW(IND) = RIW(1)
        ZIW(IND) = ZIW(1)
        IONWPTS = IND
      ELSEIF (CIONR.EQ.0) THEN
        IR = IRWALL-1
        IND = 1
        IONWI1 = 1
        DO 2035 IK = 1,NKS(IR)
c           IKOUT = IKOUTS(IK,IR)
           IKOUT = IKOUTG(IK,IR)
           IROUT = IROUTS(IK,IR)
           RIW(IND) = 0.50*(RS(IKOUT,IROUT)+RS(IK,IR))
           ZIW(IND) = 0.50*(ZS(IKOUT,IROUT)+ZS(IK,IR))
           IND = IND+1
2035    CONTINUE
        IONWI2 = IND-2
        RIW(IND) = 0.5 * (RP(1) + RP(2))
        ZIW(IND) = 0.5 * (ZP(1) + ZP(2))
        IND = IND +1
        DO 2040 ID =  2,NDSIN-1
           RIW(IND) = RP(ID)
           ZIW(IND) = ZP(ID)
           IND = IND +1
2040    CONTINUE
        RIW(IND) = 0.5 * (RP(NDSIN-1) + RP(NDSIN))
        ZIW(IND) = 0.5 * (ZP(NDSIN-1) + ZP(NDSIN))
        IND = IND +1
        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
           IR = IRTRAP+1
        elseif (cgridopt.eq.2) then
           ir = irtrap2+1
        endif
        IONTI1 = IND
        DO 2050 IK = NKS(IR),1,-1
c           IKOUT = IKINS(IK,IR)
           IKOUT = IKING(IK,IR)
           IROUT = IRINS(IK,IR)
           RIW(IND) = 0.50*(RS(IKOUT,IROUT)+RS(IK,IR))
           ZIW(IND) = 0.50*(ZS(IKOUT,IROUT)+ZS(IK,IR))
           IND = IND+1
2050    CONTINUE
        IONTI2 = IND-2
        RIW(IND) = 0.5 * (RP(NDSIN+1) + RP(NDSIN+2))
        ZIW(IND) = 0.5 * (ZP(NDSIN+1) + ZP(NDSIN+2))
        IND = IND +1
c
        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
           endid = nds-1
        elseif (cgridopt.eq.2) then
           endid = ndsin2-1
        endif
c
        DO 2055 ID = NDSIN+2,endid
           RIW(IND) = RP(ID)
           ZIW(IND) = ZP(ID)
           IND = IND +1
2055    CONTINUE
        RIW(IND) = 0.5 * (RP(endid) + RP(endid+1))
        ZIW(IND) = 0.5 * (ZP(endid) + ZP(endid+1))
        IND = IND +1
c
c       Additional Wall Points for ITER grids
c
        if (cgridopt.eq.2) then
c
           IR = IRWALL2-1
           IONWI3 = IND -1
           DO 3035 IK = 1,NKS(IR)
              IKOUT = IKOUTS(IK,IR)
              IROUT = IROUTS(IK,IR)
              RIW(IND) = RS(IK,IR)+0.50*(RS(IKOUT,IROUT)-RS(IK,IR))
              ZIW(IND) = ZS(IK,IR)+0.50*(ZS(IKOUT,IROUT)-ZS(IK,IR))
             IND = IND+1
3035       CONTINUE
           IONWI4 = IND-2
           RIW(IND) = 0.5 * (RP(ndsin2+1) + RP(ndsin2+2))
           ZIW(IND) = 0.5 * (ZP(ndsin2+1) + ZP(ndsin2+2))
           IND = IND +1
           DO 3040 ID =  ndsin2+2,NDSIN3-1
              RIW(IND) = RP(ID)
              ZIW(IND) = ZP(ID)
              IND = IND +1
3040       CONTINUE
           RIW(IND) = 0.5 * (RP(NDSIN3-1) + RP(NDSIN3))
           ZIW(IND) = 0.5 * (ZP(NDSIN3-1) + ZP(NDSIN3))
           IND = IND +1
           IR = IRTRAP+1
           IONTI3 = IND
           DO 3050 IK = 1,nks(ir)
              IKOUT = IKINS(IK,IR)
              IROUT = IRINS(IK,IR)
              RIW(IND) = RS(IK,IR)+0.50*(RS(IKOUT,IROUT)-RS(IK,IR))
              ZIW(IND) = ZS(IK,IR)+0.50*(ZS(IKOUT,IROUT)-ZS(IK,IR))
              IND = IND+1
3050       CONTINUE
           IONTI4 = IND-2
           RIW(IND) = 0.5 * (RP(NDSIN3+1) + RP(NDSIN3+2))
           ZIW(IND) = 0.5 * (ZP(NDSIN3+1) + ZP(NDSIN3+2))
           IND = IND +1
c
           DO 3055 ID = NDSIN3+2,nds-1
              RIW(IND) = RP(ID)
              ZIW(IND) = ZP(ID)
              IND = IND +1
3055       CONTINUE
           RIW(IND) = 0.5 * (RP(nds-1) + RP(nds))
           ZIW(IND) = 0.5 * (ZP(nds-1) + ZP(nds))
           IND = IND +1
c
c       End of extra for ITER grids
c
        endif
c
        RIW(IND) = RIW(1)
        ZIW(IND) = ZIW(1)
        IONWPTS = IND
c
c     Ion Wall Option 2 - polygon edges.
c
      elseif (cionr.eq.2) then
c
        IR = IRWALL-1
        IND = 1
        IONWI1 = 1
        DO IK = 1,NKS(IR)
           in = korpg(ik,ir)
           RIW(IND) = rvertp(2,in)
           ZIW(IND) = zvertp(2,in)
           IND = IND+1
        end do
        RIW(IND) = rvertp(3,in)
        ZIW(IND) = zvertp(3,in)
        IND = IND +1
        IONWI2 = IND-1
c
        DO ID =  2,NDSIN-1
           RIW(IND) = RP(ID)
           ZIW(IND) = ZP(ID)
           IND = IND +1
        end do
c slmod begin
        if (.not.nopriv) then
          if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
             IR = IRTRAP+1
          elseif (cgridopt.eq.2) then
             ir = irtrap2+1
          endif
          in = korpg(nks(ir),ir)
          RIW(IND) = rvertp(4,in)
          ZIW(IND) = zvertp(4,in)
          IND = IND +1
          IONTI1 = IND-1
          DO IK = NKS(IR),1,-1
             in = korpg(ik,ir)
             RIW(IND) = rvertp(1,in)
             ZIW(IND) = zvertp(1,in)
             IND = IND+1
          end do
        endif
c
c        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
c           IR = IRTRAP+1
c        elseif (cgridopt.eq.2) then
c           ir = irtrap2+1
c        endif
c        in = korpg(nks(ir),ir)
c        RIW(IND) = rvertp(4,in)
c        ZIW(IND) = zvertp(4,in)
c        IND = IND +1
c        IONTI1 = IND-1
c        DO IK = NKS(IR),1,-1
c           in = korpg(ik,ir)
c           RIW(IND) = rvertp(1,in)
c           ZIW(IND) = zvertp(1,in)
c           IND = IND+1
c        end do
c slmod end
        IONTI2 = IND-1
c
c slmod begin
        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3.or.
     .      cgridopt.eq.LINEAR_GRID) then
c
c        if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
c slmod end
           endid = nds-1
        elseif (cgridopt.eq.2) then
           endid = ndsin2-1
        endif
c
        DO ID = NDSIN+2,endid
           RIW(IND) = RP(ID)
           ZIW(IND) = ZP(ID)
           IND = IND +1
        end do
c
        RIW(IND) = RIW(1)
        ZIW(IND) = ZIW(1)
        IONWPTS = IND
      ENDIF

c
      if (cprint.eq.3.or.cprint.eq.9) then
c
         WRITE(6,*) 'IONWPTS:',IONWPTS,IONWI1,IONWI2,IONTI1,IONTI2,
     >            ionwi3,ionwi4,ionti3,ionti4
c
c        changed for AUG (control of wall options)
c
         write (6,'(1x,a)') 'Rwall    Zwall    (Ion wall points)'
         DO 2100 IND = 1,IONWPTS
            write (6,'(1x,2(f7.4,2x),i4)') riw(ind),ziw(ind),ind
 2100    CONTINUE
c
      endif
c
      KIND = 1
      WRITE(6,*) 'BEFORE GA15:',IONWPTS,KIND,MAXPTS

      CALL GA15A(IONWPTS,KIND,iwWORK,4*MAXPTS,iwINDW,MAXPTS,
     >             RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
c
      WRITE(6,90) iwINDW(2,1)
 90   FORMAT(' RETURN FROM GA15A, TAU = ',I10)
c
c     Set up core plasma bounds
c
      ioncpts = nks(1)+1
c
c     Code assumes same number of elements in first and second
c     rings of core plasma region.
c
      do 455 ik = 1,nks(1)
         if (cionr.eq.2) then
            in = korpg(ik,2)
            rcw(ik) = rvertp(1,in)
            zcw(ik) = zvertp(1,in)
         else
            rcw(ik) = (rs(ik,1)+rs(ik,2)) /2.0
            zcw(ik) = (zs(ik,1)+zs(ik,2)) /2.0
         endif
 455  continue
c
      rcw(ioncpts) = rcw(1)
      zcw(ioncpts) = zcw(1)
c
      write (6,'(1x,a)') 'Rwall    Zwall    (Ion wall points)'
      DO IND = 1,IONCPTS
        write (6,'(1x,2(f7.4,2x),i4)') rcw(ind),zcw(ind),ind
      end do
c
      kind = 1
c
      CALL GA15A(IONCPTS,KIND,icWORK,4*MAXPTS,icINDW,MAXPTS,
     >             RCW,ZCW,icTDUM,icXDUM,icYDUM,6)
c
c     The quantities IFXYS, IKXYS and IRXYS are now only
c     calculated for grid geometries that do not support
c     or supply the entire set of cell vertex information to DIVIMP.
c     This would be original ASDEX and ITER grids. Asdex UPGRADE
c     and JET grids do not require this information - at least
c     for calculating the location of the particle in DIVIMP.
c
c     Note: This is being left in for the time-being so that
c     the plotting routines can remain unaffected.
c
      if (xygrid.eq.1) then
C
C-----------------------------------------------------------------------
C     CALCULATE RECTANGULAR GRID EQUIVALENTS TO IK, IR INDICES
C     CALCULATE SET OF FLAGS IFXYS INDICATING PTS EXTERNAL TO SYSTEM.
C-----------------------------------------------------------------------
C
      IF (CRECT.EQ.99) THEN
        WRITE (6,'(/1X,A)') 'RECTANGULAR GRID DATA READ FROM FILE :-'
        REWIND (13)
        READ (13,*) NXS,NYS,DR,DZ,RMIN,RMAX,ZMIN,ZMAX
        IF (NXS.GT.MAXNXS) WRITE (6,*) ' ERROR! NXS =',NXS
        IF (NYS.GT.MAXNYS) WRITE (6,*) ' ERROR! NYS =',NYS
        READ (13,*) ((IKXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        READ (13,*) ((IRXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        READ (13,*) ((IFXYS(IX,IY),IX=1,NXS),IY=1,NYS)
C
      ELSE
c
c       First - go through and calculate all points inside the
c       outer wall - then find all points inside the central plasma
c       region - using the Harwell area routine.
c
        WRITE (6,'(/1X,A)') 'RECTANGULAR GRID CALCULATED AND STORED :-'
        NXS = MAXNXS
        NYS = MAXNYS
        DR = (RMAX-RMIN) / REAL(NXS-4)
        DZ = (ZMAX-ZMIN) / REAL(NYS-4)
        RMIN = RMIN - 2.0 * DR
        RMAX = RMAX + 2.0 * DR
        ZMIN = ZMIN - 2.0 * DZ
        ZMAX = ZMAX + 2.0 * DZ
C
        DO 450 IX = 1, NXS
          R = (RMAX-RMIN) * REAL(IX)/REAL(NXS) + RMIN - 0.5 * DR
          DO 450 IY = 1, NYS
            Z = (ZMAX-ZMIN) * REAL(IY)/REAL(NYS) + ZMIN - 0.5 * DZ
            BEST = HI
            DO 440 JR = 1, NRS
              DO 440 JK = 1, NKS(JR)
                IF (JK.EQ.NKS(JR).AND.JR.LT.IRSEP) GOTO 440
                DSQ = RS(JK,JR)**2 - 2.0*RS(JK,JR)*R + R**2 +
     >                ZS(JK,JR)**2 - 2.0*ZS(JK,JR)*Z + Z**2
                IF (DSQ.LT.BEST) THEN
                  BEST = DSQ
                  IK   = JK
                  IR   = JR
                ENDIF
  440       CONTINUE
C
C
            IKXYS(IX,IY) = IK
            IRXYS(IX,IY) = IR
            IFXYS(IX,IY) = 1
C
C         CALCULATE THE INSIDE/OUTSIDE PLASMA GRID - THIS IS COMPUTATION
C         INTENSIVE SO RUN IT ONCE/GEOMETRY AND USE THE SAVED FILE.
C         (UNLESS CPU USAGE IS NOT IMPORTANT - IT TAKES ABOUT 5 MIN ON
C         AN IBM RS6000 WORKSTATION)
C
c            IF (IR.EQ.1) THEN
c              IFXYS(IX,IY) = 0
c            ELSE
               CALL GA15B(R,Z,RESULT,IONWPTS,1,iwWORK,4*MAXPTS,
     >              iwINDW,MAXPTS,RIW,ZIW,iwTDUM,iwXDUM,iwYDUM,6)
              IF (RESULT.LT.0.0) IFXYS(IX,IY) = 0
c
c             This code was taken from the original ASDEX grid
c             routine. In marginal cases some core plasma
c             points were being mapped to the second ring and not
c             first ... resulting in anomalous points in the core
c             region being associated with the second ring. The
c             RICHT??? variables are a method of dealing with this.
c
c              IF (IK .LT. NKS(IR)) THEN
c                 RICHTABX = RS(IK+1,1) - RS(IK,1)
c                 RICHTABY = ZS(IK+1,1) - ZS(IK,1)
c                 RICHTCDX = RS(IK,IR) - R
c                 RICHTCDY = ZS(IK,IR) - Z
c                 IF (RICHTCDY*RICHTABX-RICHTCDX*RICHTABY .NE. 0.) THEN
c                    MUE = (ZS(IK,1)*RICHTABX + (R-RS(IK,1))*RICHTABY
c     >               - Z*RICHTABX)/ (RICHTCDY*RICHTABX -
c     >                RICHTCDX*RICHTABY)
c                    IF ((MUE .GT. 0) .AND. (MUE .LT. 1)) then
c                       IFXYS(IX,IY) = 0
c                       write(6,*) 'mue:',ix,iy,mue,ik,ir
c                    endif
c                 endif
c              ELSEIF (IK .GT. 1) THEN
c                 RICHTABX = RS(IK-1,1) - RS(IK,1)
c                 RICHTABY = ZS(IK-1,1) - ZS(IK,1)
c                 RICHTCDX = RS(IK,IR) - R
c                 RICHTCDY = ZS(IK,IR) - Z
c                 IF (RICHTCDY*RICHTABX-RICHTCDX*RICHTABY .NE. 0.) THEN
c                    MUE = (ZS(IK,1)*RICHTABX + (R-RS(IK,1))*RICHTABY
c     >               - Z*RICHTABX)/ (RICHTCDY*RICHTABX -
c     >                RICHTCDX*RICHTABY)
c                    IF ((MUE .GT. 0) .AND. (MUE .LT. 1)) then
c                       IFXYS(IX,IY) = 0
c                       write(6,*) 'mue:',ix,iy,mue,ik,ir
c                    endif
c                ENDIF
c             endif
c          ENDIF

  450   CONTINUE
C
        CALL WALLEDGE (IONWPTS,RIW,ZIW)
C
C       Exclude core plasma points from X,Y grid
c
        DO 460 IX = 1, NXS
          R = (RMAX-RMIN) * REAL(IX)/REAL(NXS) + RMIN - 0.5 * DR
          DO 460 IY = 1, NYS
            Z = (ZMAX-ZMIN) * REAL(IY)/REAL(NYS) + ZMIN - 0.5 * DZ
C
C
            IF (IFXYS(IX,IY).eq.0) goto 460

            CALL GA15B(R,Z,RESULT,IONCPTS,1,icWORK,4*MAXPTS,
     >            icINDW,MAXPTS,RCW,ZCW,icTDUM,icXDUM,icYDUM,6)
            IF (RESULT.gt.0.0) IFXYS(IX,IY) = 0
c
  460   CONTINUE

C
        REWIND (13)
        WRITE (13,*) NXS,NYS,DR,DZ,RMIN,RMAX,ZMIN,ZMAX
        WRITE (13,'(1X,23I3)') ((IKXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        WRITE (13,'(1X,23I3)') ((IRXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        WRITE (13,'(1X,23I3)') ((IFXYS(IX,IY),IX=1,NXS),IY=1,NYS)
        WRITE (6,*) NXS,NYS,DR,DZ,RMIN,RMAX,ZMIN,ZMAX
c        WRITE (6,'(2I8,'':'',3I8)') (( ix,iy,IFXYS(Ix,Iy),
c     >            IKXYS(IX,IY),
c     >            IRXYS(Ix,Iy),IX=1,NXS),IY=1,NYS)
      ENDIF

c
c     Elseif corresponding to xygrid above.
c
      elseif (xygrid.eq.0) then
c
c        Set values for dr and dz
c        And adjust R,Z min/max values.
c
c
         DR = (RMAX-RMIN) / REAL(maxgxs-4)
         DZ = (ZMAX-ZMIN) / REAL(maxgys-4)
         RMIN = RMIN - 2.0 * DR
         RMAX = RMAX + 2.0 * DR
         ZMIN = ZMIN - 2.0 * DZ
         ZMAX = ZMAX + 2.0 * DZ

c
c     Endif corresponding to cgridopt above.
c
      endif

C
c      WRITE (6,*) 'GEO1'
c      DO 654 IY = 1,NYS
c        WRITE (6,'(120I1)') (IFXYS(IX,IY),IX=1,NXS/3)
c 654  CONTINUE
c      WRITE (6,*) 'GEO2'
c      DO 655 IY = 1,NYS
c        WRITE (6,'(120I1)') (IFXYS(IX,IY),IX=NXS/3+1,(2*NXS)/3)
c  655  CONTINUE
c      WRITE (6,*) 'GEO3'
c      DO 656 IY = 1,NYS
c        WRITE (6,'(120I1)') (IFXYS(IX,IY),IX=((2*NXS)/3)+1,NXS)
c 656  CONTINUE
c      WRITE (6,*) 'FINISHED'
C
      RETURN
      END
C
C
      SUBROUTINE WALLEDGE(IONWPTS,RIW,ZIW)
      implicit none
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "CGEOM"
      include 'cgeom'
      include 'divxy'
      INTEGER IONWPTS
      REAL    RIW(IONWPTS),ZIW(IONWPTS)
C
C     WALLEDGE: THIS ROUTINE REFINES THE EDGE OF THE
C               REGION DEFINING THE PLASMA - SO THAT
C               MOST IX,IY BINS WITH A PROTRUSION INTO
C               THE PLASMA REGION ARE MARKED AS BEING PART
C               OF THE PLASMA IN THE ARRAY IFXYS.
C
C               THIS IS DONE IN A SIMPLE-MINDED WAY BY
C               FOLLOWING ALONG EACH WALL SEGMENT USING
C               A SMALL STEP SIZE AND MARKING EVERY IX,IY BIN
C               IN WHICH THESE POINTS FALL AS BEING INSIDE THE
C               PLASMA.
C
C     DAVID ELDER , FEB 13 , 1992
C
C     NOTE: THIS IS NOT AN "ELEGANT" SOLUTION BUT IT SHOULD
C           DEAL WITH THE PROBLEMS CURRENTLY OBSERVED IN
C           NEUTRAL LAUNCHES FROM THE PLATES - WHICH HAVE
C           A SIGNIFICANT COMPONENT HITTING THE PLATES AFTER
C           ONLY A FEW TIME STEPS.
C
      INTEGER IX,IY,ID,IN
      REAL NPT
      PARAMETER (NPT=1000.0)
      REAL    RSTART,ZSTART,REND,ZEND,RSTEP,ZSTEP,R,Z
C
      DO 10 ID =1 , IONWPTS
        RSTART = RIW(ID)
        ZSTART = ZIW(ID)
        IF (ID.EQ.IONWPTS) THEN
          REND = RIW(1)
          ZEND = ZIW(1)
        ELSE
          REND = RIW(ID+1)
          ZEND = ZIW(ID+1)
        ENDIF
        RSTEP = (REND-RSTART)/NPT
        ZSTEP = (ZEND-ZSTART)/NPT
        DO 20 IN = 0,INT(NPT)
          R = RSTART + IN * RSTEP
          Z = ZSTART + IN * ZSTEP
          IX    = MAX (1, MIN (NXS, INT((R-RMIN)/DR)+1))
          IY    = MAX (1, MIN (NYS, INT((Z-ZMIN)/DZ)+1))
          IFXYS(IX,IY) = 1
20      CONTINUE
10    CONTINUE
      RETURN
      END
C
C
C
      INTEGER FUNCTION FINDNUM(ARS,NARS,TNUM)
      implicit none
C     INCLUDE "PARAMS"
      include 'params'
      REAL ARS(MAXNRS,5)
      INTEGER NARS,TNUM
C
C     THIS ROUTINE RETURNS THE INDEX INTO THE ARRAY ARS POINTING
C     TO THE LINE WHOSE CONTENTS RELATE TO THE QUANTITY TNUM.
C
C     THE FIRST ENTRY IN ARS IS THE NUMBER TNUM IN A REAL REPRESENTATION
C
      INTEGER I
C
      FINDNUM = 0
C
C      WRITE(6,*) 'FINDNUM:', NARS,TNUM
C
      DO 100 I = 1,NARS
         IF (TNUM.EQ.INT(ARS(I,1))) THEN
            FINDNUM = I
            GOTO 200
         ENDIF
 100  CONTINUE
 200  CONTINUE
      RETURN
      END
c
c
c
      subroutine gridcoords (ix,iy,ik,ir,in)
      integer ix,iy,ik,ir,in
      include 'params'
      include 'divxy'
c
c     Return the ik,ir cooridnates for the ix,iy bin.
c     The separate routine is necessary because div and
c     out now differ in their definition of the following
c     arrays.
c
      ik = ikxys(ix,iy)
      ir = irxys(ix,iy)
      in = ifxys(ix,iy)
      return
      end
c
c
c
      subroutine dovessel
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c***********************************************************************
c
c     DOVESSEL:
c
c     This subroutine takes the series of R.Z points that have been
c     read in from the GRID2D geometry file - which specifiy the
c     VESSEL wall location and map these onto the two arrays
c     WALLCO and WALLCO2 that contain these points for the outer
c     wall and trap wall regions respectively. The TARGET segments
c     themselves are used in the DOWALL routine to define the
c     vessel wall in the target region. This routine splits the
c     GRID2D vessel wall into OUTER WALL and PRIVATE PLASMA WALL
c     segments and discards vessel wall coordinates that lie
c     on or behind the targets.
c
c     David Elder,       June 27, 1995
c
c***********************************************************************
c
c     Local variables
c
      integer i,j,k,ik,ir,in,id,icnt,ind1,ind2
      real    solcorners(2,2),trapcorners(2,2)
      logical  intarg
      integer  nearwall
      external nearwall,intarg
c
c     The basic algorithm is to find the Vessel Wall coordinate
c     that is closest to each of the target end-points without
c     being on or behind the target itself.
c
c     The method used here is to find the closest point to the end-point
c     first - then check whether it is on the side of this point towards
c     the rest of the target - or if it is on the side of this point away
c     from the target. If it is away from the target - it is taken to
c     be the point that will connect to the target - if not then the
c     nearby points on the vessel wall - stepping up and down in indices -
c     are checked until one away from the target is found.
c
c
c     First - calculate the actual corners of the target that are
c     relevent - depending on target options.
c
      if (ctargopt.eq.6) then
c
c        This uses the end polygon corners of the next to last
c        ring - since polygons do not exist for the virtual rings.
c
         in = korpg(ikds(2),irds(2))
         solcorners(1,1) = rvertp(3,in)
         solcorners(1,2) = zvertp(3,in)
         in = korpg(ikds(nds-1),irds(nds-1))
         solcorners(2,1) = rvertp(2,in)
         solcorners(2,2) = zvertp(2,in)
         in = korpg(ikds(ndsin-1),irds(ndsin-1))
         trapcorners(1,1) = rvertp(4,in)
         trapcorners(1,2) = zvertp(4,in)
         in = korpg(ikds(ndsin+2),irds(ndsin+2))
         trapcorners(2,1) = rvertp(1,in)
         trapcorners(2,2) = zvertp(1,in)
      else
         solcorners(1,1) = rp(2)
         solcorners(1,2) = zp(2)
         solcorners(2,1) = rp(nds-1)
         solcorners(2,2) = zp(nds-1)

         IF (CTRAP.EQ.0.OR.CTRAP.EQ.2.or.ctrap.eq.3) THEN
            trapcorners(1,1) = rp(ndsin-1)
            trapcorners(1,2) = zp(ndsin-1)
            trapcorners(2,1) = rp(ndsin+2)
            trapcorners(2,2) = zp(ndsin+2)
         ELSEIF (CTRAP.EQ.1) THEN
            trapcorners(1,1) = rp(ndsin)
            trapcorners(1,2) = zp(ndsin)
            trapcorners(2,1) = rp(ndsin+1)
            trapcorners(2,2) = zp(ndsin+1)
         ENDIF
      endif
c
c     Set up NEUTRAL OUTER wall if required.
c
      if (cneur.eq.4) then
c
c        The OUTER wall numbers clockwise from the
c        OUTER target - i.e. IK = 1, and ID = NDSIN+1,NDS
c        So start at the second solcorner.
c
c        Find indices of VESSEL coordinates closest to
c        ends of targets.
c
c slmod begin
c *NOTE*
c 
c I think this is where the trouble was with Adams 116208.03750 grid, where
c the neutral wall was not being processed properly.  The call to INTARG in
c the NEARWALL routine was not correctly identifying some wall points
c as behind the outer target. -SL, April 30, 2004
c
c slmod end
         ind1 = nearwall(solcorners(2,1),solcorners(2,2),
     >                 trapcorners(2,1),trapcorners(2,2))
         ind2 = nearwall(solcorners(1,1),solcorners(1,2),
     >                 trapcorners(1,1),trapcorners(1,2))
c
c        Calculate number of elements in the array
c
         if (ind1.lt.ind2) then
            nwall = ind2 - ind1 + 1
         elseif (ind1.eq.ind2) then
            write (6,*) 'Error in DOVESSEL: ind1 = ind2',ind1,ind2
            write (6,*) 'CNEUR reset to 0 since wall data is invalid'
            nwall = 0
            cneur = 0
         else
            nwall = (nves - 1) - ind1 + 1 + ind2
         endif
c
         write (6,*) 'Vessel:',nves,nwall,ind1,ind2
c
c        Now - loop through from ind1 to ind2 assigning the
c        values to the wallco array.
c
         icnt = 1
         in = ind1
c
c        Loop through assigning values.
c
 100     if (icnt .gt.nwall) goto 200
         wallco(icnt,1) = rves(in)
         wallco(icnt,2) = zves(in)
         in = in + 1
         icnt = icnt + 1
c
c        Ignore first point since it is a repeat of the last 
c
         if (in.gt.nves) in = 2
         goto 100
c
 200     continue
c
      endif

c
c     Set up NEUTRAL PRIVATE PLASMA wall if required.
c
      if (ctrap.eq.4) then
c
c        The PRIVATE PLASMA wall numbers clockwise from the
c        INNER target - i.e. IK = 1, and ID = 1,NDSIN
c        So start at the first trapcorner.
c
c        FInd indices of VESSEL coordinates closest to
c        ends of targets.
c
         ind1 = nearwall(trapcorners(1,1),trapcorners(1,2),
     >                 solcorners(1,1),solcorners(1,2))
         ind2 = nearwall(trapcorners(2,1),trapcorners(2,2),
     >                 solcorners(2,1),solcorners(2,2))
c
c        Verify that the ind1 and ind2 values found are
c        not inside the other target - if there are NO
c        points between the ends of the PFZ targets then
c        the trap wall option is changed to 2.   
c
         if (intarg(ind1,trapcorners(2,1),trapcorners(2,2),
     >                 solcorners(2,1),solcorners(2,2)).or.
     >       intarg(ind2,trapcorners(2,1),trapcorners(2,2),
     >                 solcorners(2,1),solcorners(2,2))) then

c
c           Error - one or both of the PFZ endpoints lie inside the
c           adjacent target.  
c
            write (6,'(a,2i6)') 'Error in DOVESSEL - PFZ:'//
     >                  ' ind1 or ind2 is inside other target',ind1,ind2
            write (6,'(a)') 'CTRAP reset to 2 since wall data'//
     >                      ' is invalid'
            nwall2 = 0
            ctrap = 2
c
         else
c
c           Now - loop through from ind1 to ind2 assigning the
c           values to the wallco2 array.
c
            icnt = 1
            in = ind1
c
c           Calculate number of elements in the array
c
c           ONE point IS valid for the PFZ wall
c
            if (ind1.le.ind2) then
               nwall2 = ind2 - ind1 + 1
c            elseif (ind1.eq.ind2) then
c               write (6,*) 'Error in DOVESSEL: ind1 = ind2',ind1,ind2
c               write (6,*) 'CTRAP reset to 2 since wall'//
c     >                     ' data is invalid'
c            nwall2 = 0
c            ctrap = 2
            else
               nwall2 = (nves - 1) - ind1 + 1 + ind2
            endif
c
c           Loop through assigning values.
c
 300        if (icnt .gt.nwall2) goto 400
            wallco2(icnt,1) = rves(in)
            wallco2(icnt,2) = zves(in)
            in = in + 1
            icnt = icnt + 1
c
c           Ignore first point since it is a repeat of the last 
c
            if (in.gt.nves) in = 2
            goto 300
c
 400        continue
c
         endif
c
      endif
c
      if (cprint.eq.3.or.cprint.eq.9) then
c
         write (6,*) 'Wallco:',nwall
         do in = 1,nwall
            write (6,'(i4,2(1x,g12.6))') in,
     >        wallco(in,1),wallco(in,2)
         end do
c
         write (6,*) 'Wallco2:',nwall2
         do in = 1,nwall2
            write (6,'(i4,2(1x,g12.6))') in,
     >        wallco2(in,1),wallco2(in,2)
         end do
c
         write (6,*) 'Vessel:',nves
         do in = 1,nves
            write (6,'(i4,2(1x,g12.6))') in,
     >        rves(in),zves(in)
         end do
c
      endif
c
      return
      end
c
c
c
      integer function nearwall (rp1,zp1,rp2,zp2)
      implicit none
      real rp1,zp1,rp2,zp2
      include 'params'
      include 'cgeom'
c
c     Find the index of the element of the nearest
c     RVES,ZVES point that is outside the target.
c     The point of interest is RP1,ZP1 and the
c     line of the target is defined by joining
c     (rp2,zp2) and (rp1,zp1)
c
      real mindr2,dr2
      integer i,itmp,mini,alt
      logical intarg
      external intarg
c
      mindr2 = HI
      mini  = 0
c
      do i = 1,nves-1
c
c     I = 1 and I = NVES are the same point since the
c     vessel wall is a closed figure.
c
         dr2 = (rp1-rves(i))**2 + (zp1-zves(i))**2
         if (dr2.lt.mindr2) then
            mindr2 = dr2
            mini = i
         endif
      end do
c
c     Is the point towards or away from the target?
c
      alt = 0
      itmp = mini
 100  if (intarg(mini,rp1,zp1,rp2,zp2)) then
         alt = alt + 1
         if (alt.gt.nves/2) then
            write (6,*)  'ERROR in FINDVESSEL:'
     >                   //' Problem determining end-point'
            stop
         endif
         itmp = mini+alt
         if (itmp.gt.nves-1) then
            itmp = itmp-(nves-1)
         endif
         if (.not.intarg(itmp,rp1,zp1,rp2,zp2)) goto 200
         itmp = mini-alt
         if (itmp.lt.1) then
            itmp = itmp + (nves-1)
         endif
         if (.not.intarg(itmp,rp1,zp1,rp2,zp2)) goto 200
         goto 100
      endif
c
c     Unless an error has occurred - want to return the
c     value of itmp.
c
 200  nearwall = itmp
      return
      end
c
c
c
      logical function intarg(itmp,rp1,zp1,rp2,zp2)
      implicit none
      integer itmp
      real rp1,zp1,rp2,zp2
      include 'params'
      include 'cgeom'
c
c     This function tests to see if the point is inside
c     the target.
c
c     Use the COSINE law to determine angle between lines
c     Could also use equations of lines or direction cosines.
c
      real a2,b2,c2,calph
c
      a2 = (rves(itmp)-rp1)**2 + (zves(itmp)-zp1)**2
      b2 = (rp2-rp1)**2 + (zp2-zp1)**2
      c2 = (rves(itmp)-rp2)**2 + (zves(itmp)-zp2)**2
      if (a2.ne.0.0.and.b2.ne.0.0) then
         calph = (c2-a2-b2)/(-2.0 * sqrt(a2) * sqrt(b2))
         if (calph.gt.0.0) then
            intarg = .true.
         else
            intarg = .false.
         endif
      elseif (a2.eq.0.0) then
         intarg = .false.
      elseif (b2.eq.0.0) then
         intarg = .true.
         write (6,*) 'ERROR in INTARG (rp1,zp1) = (rp2,zp2)'
         stop
      endif
c
      return
      end
c
c
c
      subroutine fixwallco
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     FIXWALLCO:
c
c     This routine fixes the specified wall coordinates to add
c     the corners of the target polygons to the wall specification
c     if CTARGOPT 6 has been specified.
c
      integer i,in
c
c        Check to see if the target option is 6 - if it is
c        add the corners of the target polygons to the wall
c        definition so that the target and wall code interaction
c        in the DOWALL routine works correctly.
c
         if (ctargopt.eq.6.and.nwall.gt.0) then

            do i = nwall,1,-1
               wallco(i+1,1) = wallco(i,1)
               wallco(i+1,2) = wallco(i,2)
            end do
c
            nwall = nwall + 2
c
c           The OUTER target corner becomes the first point -
c           i.e. the cell where ID = NDS -1 and IK = 1
c
c           This is RVERTP(2,in) and ZVERTP(2,in)
c
            in = korpg(ikds(nds-1),irds(nds-1))
            wallco(1,1) = rvertp(2,in)
            wallco(1,2) = zvertp(2,in)
c
c            write (6,*) 'p1'
c            write (6,'(2(2x,g18.10))')
c     >           (rvertp(i,in),zvertp(i,in),i=1,4)
c
c           Add the other target corner - for ID = 2 and
c           IK = NKS(IR)
c
c           This is RVERTP(3,in) and ZVERTP(3,in)
c
            in = korpg(ikds(2),irds(2))
            wallco(nwall,1) = rvertp(3,in)
            wallco(nwall,2) = zvertp(3,in)
c
c            write (6,*) 'p2'
c            write (6,'(2(2x,g18.10))')
c     >           (rvertp(i,in),zvertp(i,in),i=1,4)
c
         endif

         if (cprint.eq.3.or.cprint.eq.9) then 
            write (6,*) 'FIXWALLCO:NWALL:',nwall  
         endif  
c
      return
      end
c
c
c
      subroutine fixwallco2
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     FIXWALLCO2:
c
c     This routine fixes the specified wall coordinates to add
c     the corners of the private plasma target polygons to the
c     trap wall specification if CTARGOPT 6 has been specified.
c
      integer i,in
c
c        Check to see if the target option is 6 - if it is
c        add the corners of the target polygons to the wall
c        definition so that the target and wall code interaction
c        in the DOWALL routine works correctly.
c
         if (ctargopt.eq.6.and.nwall2.gt.0) then
            do i = nwall2,1,-1
               wallco2(i+1,1) = wallco2(i,1)
               wallco2(i+1,2) = wallco2(i,2)
            end do
c
            nwall2 = nwall2 + 2
c
c           The INNER target trap corner becomes the first point -
c           i.e. the cell where ID = NDSIN -1 and IK = NKS(IR)
c
c           This is RVERTP(4,in) and ZVERTP(4,in)
c
            in = korpg(ikds(ndsin-1),irds(ndsin-1))
            wallco2(1,1) = rvertp(4,in)
            wallco2(1,2) = zvertp(4,in)
c
c           Add the other target corner - for ID = NDSIN+2  and
c           IK = 1
c
c           This is RVERTP(1,in) and ZVERTP(1,in)
c
            in = korpg(ikds(ndsin+2),irds(ndsin+2))
            wallco2(nwall2,1) = rvertp(1,in)
            wallco2(nwall2,2) = zvertp(1,in)
         endif
c
         if (cprint.eq.3.or.cprint.eq.9) then 
            write (6,*) 'FIXWALLCO2:NWALL:',nwall2  
         endif  
c
      return
      end
c
c
c
      subroutine dotemp
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c slmod begin
      include 'slcom'
c slmod end
c
c     DOTEMP: This routine takes the input temperatures for the wall
c             and targets and any segment over-rides  and assigns
c             them to the appropriate wall temperature arrays.
c
c             The default values are loaded first before any
c             specific wall segments over-rides are put in place.
c
c             David Elder, March 13, 1997
c
      integer in,id,id1,id2
c
c     Set default on targets.
c
      do in = 1,nds
         tempds(in) = ctargt
      end do
c slmod begin
      if (grdnmod.ne.0) then
c...    Generalized geometry:

c...    Everywhere on wall:
        do in = 1, wallpts
           wallpt(in,19) = cwallt
        end do
c...    PFZ wall, which is still well defined:
        do in = wltrap1,wltrap2
           wallpt(in,19) = cwallt
        end do
c...    Target segments:
        do in = 1, wallpts
           if (wallpt(in,18).NE.0.0) wallpt(in,19) = ctargt
        end do
      else
c
c       Set default on main wall
c
        do in = wlwall1,wlwall2
           wallpt(in,19) = cwallt
        end do
c
c       Set default on PP Wall
c
        do in = wltrap1,wltrap2
           wallpt(in,19) = cwallp
        end do
c
c       Set default for target segments of wall
c
        do in = wlwall2+1,wltrap1-1
           wallpt(in,19) = ctargt
        end do
c
        do in = wltrap2+1,wallpts
           wallpt(in,19) = ctargt
        end do
      endif
c
cc
cc     Set default on main wall
cc
c      do in = wlwall1,wlwall2
c         wallpt(in,19) = cwallt
c      end do
cc
cc     Set default on PP Wall
cc
c c     do in = wltrap1,wltrap2
c         wallpt(in,19) = cwallp
c      end do
cc
cc     Set default for target segments of wall
cc
c      do in = wlwall2+1,wltrap1-1
c         wallpt(in,19) = ctargt
c      end do
cc
c      do in = wltrap2+1,wallpts
c         wallpt(in,19) = ctargt
c      end do
c slmod end
c
c     Now - impose any over-rides in the walltemp array.
c
      if (nwltemp.gt.0) then
c
         do in = 1,nwltemp
c
c           Get indices
c
            id1 = walltemp(in,1)
            id2 = walltemp(in,2)
c
            do id = id1,id2
c
c              Check if index is in valid range
c              - then apply temperature over-rides
c
               if (id.ge.1.and.id.le.wallpts) then

                   wallpt(id,19) = walltemp(in,3)
c
c                  If there is a corresponding target element
c                  apply over-ride to target too.
c
                   if (wallpt(id,18).gt.0) then
                      tempds(int(wallpt(id,18))) = walltemp(in,3)
                   endif

               endif
c
            end do
c
         end do
c
      endif
c
c
c
      return
      end
c
c
c
      subroutine check_wallco(wallco,nwall)
      implicit none
      include 'params'
c
      integer nwall
      real wallco(maxpts,2)
c
c     Local variable
c
      integer i,j,new_nwall
      real new_wallco(maxpts,2)
c
      new_nwall = 0
c
      do i = 1, nwall-1
c
c        Check for repeated segments - only add unique ones
c
         if (.not.((wallco(i,1).eq.wallco(i+1,1)).and.
     >       (wallco(i,2).eq.wallco(i+1,2)))) then

             new_nwall = new_nwall +1

             new_wallco(new_nwall,1) = wallco(i,1)
             new_wallco(new_nwall,2) = wallco(i,2)

         else

           write(6,'(a,2i6,4(1x,g16.9))')
     >         'CHECK_WALLCO: DUPLICATE SEGMENT',nwall,new_nwall,
     >         i,wallco(i,1),wallco(i+1,1),wallco(i,2),wallco(i+1,2)


         endif
c
       end do
c
c      Add last point always
c
       new_nwall = new_nwall +1

       new_wallco(new_nwall,1) = wallco(nwall,1)
       new_wallco(new_nwall,2) = wallco(nwall,2)
c
c      Assign revised wall coordinates if necessary
c
       if (new_nwall.ne.nwall) then
c
          nwall = new_nwall
c
          do i = 1,nwall
c
             wallco(i,1) = new_wallco(i,1)
             wallco(i,2) = new_wallco(i,2)
c
          end do
c
       endif
c
       return
       end
c
c
c
      subroutine check_wall(rw,zw,pcnt)
      implicit none
c
      include 'params'
c
      integer pcnt,maxarrpts
      real rw(maxpts),zw(maxpts)
c
c     Local variables
c
      integer i,j,new_nwall
      real new_rw(maxpts),new_zw(maxpts)
c
      new_nwall = 0
c
      do i = 1, pcnt-1
c
c        Check for repeated points - only add unique ones
c
         if (.not. ((rw(i).eq.rw(i+1)).and.
     >       (zw(i).eq.zw(i+1) ) ) ) then

             new_nwall = new_nwall +1

             new_rw(new_nwall) = rw(i)
             new_zw(new_nwall) = zw(i)

         endif
c
       end do
c
c      Add last point always   
c
       new_nwall = new_nwall +1
       new_rw(new_nwall) = rw(pcnt)
       new_zw(new_nwall) = zw(pcnt)
c
c      Assign revised wall coordinates if necessary
c
       if (new_nwall.ne.pcnt) then
c
          pcnt = new_nwall
c
          do i = 1,pcnt
c
             rw(i) = new_rw(i)
             zw(i) = new_zw(i)
c
          end do
c
       endif
c
       return
       end
c
c
c 
       subroutine find_nearest_boundary_cell(ikn,irn,r,z,dist,ri,zi)
       implicit none
       integer ikn,irn
       real r,z,dist,ri,zi  
       include 'params'
       include 'cgeom'
c slmod begin
       include 'slcom'
c slmod end
c
c      This routine finds the ik,ir indices of the plasma 
c      cell on the boundary of the plasma which is closest to 
c      the r,z location. It then calculates the minimum distance
c      to the segment - either the perpendicular distance if the 
c      intersection lies on the axis of the cell or the distance to the
c      end points of the cell center axis.
c
       real   eps
       real   r1,z1,r2,z2  
       parameter(eps=1.0e-6) 
       real*8 best,dsq
       integer ik,ir,in,ik1,ir1
       real dt,d1,d2
c
       BEST = HI
       DSQ  = HI
c
c      Loop over IRWALL -1 first
c
c slmod begin
       if (grdnmod.NE.0) then
c...     Generalized geometry:
         ir1 = irwall 
c
         do ik1 = 1,nks(ir1) 
            ik = ikins(ik1,ir1)
            ir = irins(ik1,ir1)
c
            DSQ = (rs(ik,ir)-R) ** 2 + (zs(ik,ir)-Z) ** 2
            IF (DSQ.LT.BEST) THEN
              BEST = DSQ
              ikn = ik
              irn = ir  
            ENDIF
c
         end do
       else
         ir = irwall -1
c
         do ik = 1,nks(ir) 
c
            DSQ = (rs(ik,ir)-R) ** 2 + (zs(ik,ir)-Z) ** 2
            IF (DSQ.LT.BEST) THEN
              BEST = DSQ
              ikn = ik
              irn = ir  
            ENDIF
c
         end do
       endif
c
c       ir = irwall -1
c
c       do ik = 1,nks(ir) 
c
c          DSQ = (rs(ik,ir)-R) ** 2 + (zs(ik,ir)-Z) ** 2
c          IF (DSQ.LT.BEST) THEN
c            BEST = DSQ
c            ikn = ik
c            irn = ir  
c          ENDIF
c
c       end do
c slmod end
c
c      Loop over irtrap + 1 next
c
       ir = irtrap+1
c
c      Account for limiter case with no PFZ
c
       if (irtrap+1.lt.nrs) then 

          do ik = 1,nks(ir) 
c
             DSQ = (rs(ik,ir)-R) ** 2 + (zs(ik,ir)-Z) ** 2
             IF (DSQ.LT.BEST) THEN
               BEST = DSQ
               ikn = ik
               irn = ir  
             ENDIF
c
          end do
c
       endif
c
c      Loop over first target skipping virtual rings
c
       do ir = irsep,nrs
c
          ik = 1
c
          if (ir.eq.irwall.or.ir.eq.irtrap) cycle 
c
          DSQ = (rs(ik,ir)-R) ** 2 + (zs(ik,ir)-Z) ** 2
          IF (DSQ.LT.BEST) THEN
            BEST = DSQ
            ikn = ik
            irn = ir  
          ENDIF
c
       end do
c
c      Loop over second target skipping virtual rings
c
       do ir = irsep,nrs
c
          ik = nks(ir)
c
          if (ir.eq.irwall.or.ir.eq.irtrap) cycle 
c
          DSQ = (rs(ik,ir)-R) ** 2 + (zs(ik,ir)-Z) ** 2
          IF (DSQ.LT.BEST) THEN
            BEST = DSQ
            ikn = ik
            irn = ir  
          ENDIF
c
       end do
c
c      The results in ikn irn should be the plasma boundary cell closest 
c      to the test point
c
c      Calculate distance and intersection point
c
c begin afmod
c Geier IPP/02 ik, ir -> ikn, ikr in next line -bug? -SL
c
c       jdemod - incorporated the change to use ikn,irn since the other is clearly a bug
c
c       IF (ippchange) THEN
c
         in = korpg(ikn,irn)
c
c       ELSE
c         in = korpg(ik,ir)
c       ENDIF
c
c       in = korpg(ik,ir)
c afmod end
c
       if (in.ne.0) then 
c
c         Calculate distance from cell parallel axis
c
          r1 = (rvertp(1,in) + rvertp(2,in)) /2.0
          r2 = (rvertp(3,in) + rvertp(4,in)) /2.0
          z1 = (zvertp(1,in) + zvertp(2,in)) /2.0
          z2 = (zvertp(3,in) + zvertp(4,in)) /2.0
c
c         Call pointdist to find distance and intersection
c
c afmod begin
c          IF (ippchange) THEN
c Geier IPP/02 corrected order of variables in call of pintdist -bug
c
c           jdemod - fix bug in calling parameters
c
            call pointdist(r1,r2,z1,z2,r,z,dist,ri,zi)
c
c          ELSE
c            call pointdist(r1,z1,r2,z2,r,z,dist,ri,zi)
c          ENDIF  
c
c          call pointdist(r1,z1,r2,z2,r,z,dist,ri,zi)
c afmod end
c
          dt = (r1-r2)**2 + (z1-z2)**2
          d1 = (ri-r1)**2 + (zi-z1)**2
          d2 = (ri-r2)**2 + (zi-z2)**2
c
c         Intersection is not between end points
c
          if ((d1+d2).gt.(dt+eps)) then 
c
             if (d1.lt.d2) then 
c
                dist = sqrt((r-r1)**2+(z-z1)**2)
c
                ri = r1
                zi = z1 
c
             else
c
                dist = sqrt((r-r2)**2+(z-z2)**2)
                ri = r2 
                zi = z2
c
             endif
c
          endif  

c
c      No polygon - just use distance to cell center
c
       else

         dist = sqrt((r-rs(ik,ir))**2+(z-zs(ik,ir))**2)
         ri = rs(ik,ir)
         zi = zs(ik,ir)
c  
       endif
c
       return
       end












