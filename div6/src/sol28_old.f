c     -*Fortran*-
c
c ======================================================================
c
      SUBROUTINE FindCell(ind0,ind1,irgive,ikcell,ircell)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      IMPLICIT none

      INTEGER ind0,ind1,irgive,ikcell(3),ircell(3)

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'slcom'


      INTEGER ik,ir,id,i1,ik1(MAXNRS)
      REAL    dist,dist1(MAXNRS),clcell(3)
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd

      ik1 = 0
      dist1 = 0.0

      DO ir = NINT(osms28(ind1,11)), NINT(osms28(ind1,12))
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        DO i1 = ind0+1, ind1

          a1 = DBLE(osms28(i1-1,5))
          a2 = DBLE(osms28(i1-1,6))
          b1 = DBLE(osms28(i1  ,5))
          b2 = DBLE(osms28(i1  ,6))
          dist = REAL(DSQRT((a1 - b1)**2 + (a2 - b2)**2))

          DO ik = 1, nks(ir)

            id = korpg(ik,ir)
            c1 = 0.5D0 * DBLE(rvertp(1,id) + rvertp(2,id))
            c2 = 0.5D0 * DBLE(zvertp(1,id) + zvertp(2,id))
            d1 = 0.5D0 * DBLE(rvertp(3,id) + rvertp(4,id))
            d2 = 0.5D0 * DBLE(zvertp(3,id) + zvertp(4,id))

            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
            IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .          tcd.GE.0.0.AND.tcd.LT.1.0) THEN
              ik1(ir) = ik
              dist1(ir) = dist1(ir) + REAL(tab) * dist
            ENDIF

          ENDDO

          IF (ik1(ir).EQ.0) dist1(ir) = dist1(ir) + dist

        ENDDO

      ENDDO

c...  Sort intesections:
      ikcell = 0
      ircell = 0
      clcell = 0.0
      DO ir = NINT(osms28(ind1,11)), NINT(osms28(ind1,12))
c        WRITE(0,*) 'PICKENS:',ir,ik1(ir),dist1(ir)
        IF (idring(ir).EQ.BOUNDARY.OR.ik1(ir).EQ.0) CYCLE
        IF     (ir.EQ.irgive) THEN
          ikcell(1) = ik1(ir)
          ircell(1) = ir
        ENDIF
        IF (clcell(2).EQ.0.0.OR.dist1(ir).LT.clcell(2)) THEN
          ikcell(2) = ik1(ir)
          ircell(2) = ir
          clcell(2) = dist1(ir)
        ENDIF
        IF (clcell(3).EQ.0.0.OR.dist1(ir).GT.clcell(3)) THEN
          ikcell(3) = ik1(ir)
          ircell(3) = ir
          clcell(3) = dist1(ir)
        ENDIF
      ENDDO

c      WRITE(0,*) 'IKDATA:',ikcell(1),ikcell(2),ikcell(3)
c      WRITE(0,*) '      :',ircell(1),ircell(2),ircell(3)
c      WRITE(0,*) '      :',ind0,ind1
         
      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE FindS28Parameters_V3(ir,te,ne,nf,s,new)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      IMPLICIT none
       
      INTEGER ir
      LOGICAL new
      REAL    te(0:6),ne(0:6),nf,s(0:6),isat(0:6)

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'slcom'

      REAL    GetRelaxationFraction

      INTEGER i0,i1,i2,i3,i4,ifit,ik,id,index,mode,ik1,ir1,id1,
     .        ikcell(3),ircell(3)
      LOGICAL tc,nc,density,tetarget
      REAL    frac,t0,t1,n0,n1,A,B,C,tetmp1,tetmp2,coord,expon,
     .        psin0,psin1,psin2,
     .        prb1,tmp1,val,val0,val1,val2,p(0:5),v0,v1,v2
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd

      s (0) = 0.0
      ne(1) = 0.0
      te(1) = 0.0
      te(2) = 0.0
      te(3) = 0.0
      te(4) = 0.0
      te(5) = 0.0


      frac = GetRelaxationFraction()


      IF (new.AND.frac.NE.0.0) THEN
c        WRITE(0,*) 'FRAC:',frac,rel_step,rel_nstep
      ENDIF


      DO i1 = 2, osmns28
        i0 = i1 - 1
        IF ((osms28(i0,1).NE.osms28(i1,1)).OR.osms28(i1,1).EQ.0.0.OR.
     .      osms28(i1,2).EQ.0.0) CYCLE

c...    Do not apply data if IR is outside specified range:
        IF ((osms28(i1,11).NE.0.0.AND.REAL(ir).LT.osms28(i1,11) ).OR.
     .      (osms28(i1,12).NE.0.0.AND.REAL(ir).GT.osms28(i1,12)))
     .    CYCLE

c...    Check that rings from different grid regions are not in the same group
c       of rings:
        DO i2 = NINT(osms28(i1,11)), NINT(osms28(i1,12))-1
          IF (ringtype(i2).NE.ringtype(i2+1)) THEN
c            WRITE(0,*) 'RINGTYPE:',i1
c            WRITE(0,*) 'RINGTYPE:',i2,ringtype(i2)
c            WRITE(0,*) 'RINGTYPE:',i2+1,ringtype(i2+1)
c            WRITE(0,*) 'RINGTYPE:',osmns28
c            CALL ER('FindS28Parameters_V3','Mixed RINGTYPE',*99)
            IF (sloutput.AND.ringtype(i2).NE.CORE) THEN
              WRITE(0,*)
              WRITE(0,*) '-----------------------------------------'
              WRITE(0,*) ' THAT FUNNY THING ABOUT MIXING REGIONS!? '
              WRITE(0,*) '-----------------------------------------'
              WRITE(0,*)
            ENDIF
          ENDIF      
        ENDDO


        IF (.FALSE..AND.osms28(i1,2).EQ.99.0) THEN
        ELSE

c MODE                          P1           P2
c   1 - power law               coordinate   index
c   2 - exponential v0-v2       coordinate   index
c   3 - exponential to infinity
c   4 - from probe data         coordinate   probe number
c   5 - something strange...
c
c   coord = 1 - linear on line segment
c         = 2 - linear on line segment, but from first to last ring intersection
c         = 3 - PSIn 
c         = 4 - RHO
c

          index = NINT(osms28(i1,1)) - 1
          mode  = NINT(osms28(i1,2))
          coord = osms28(i1,3)
          expon = osms28(i1,4)

c...  Decide if specified upstream data is density or pressure:
          density = .TRUE.   ! *** I DON'T LIKE THIS SOLUTION ***
          IF (index.EQ.3.AND.
     .        (ringtype(ir).EQ.SOL1.AND.s28nemode   .EQ.1.OR.
     .         ringtype(ir).EQ.PFZ .AND.s28nemodepfz.EQ.1)) 
     .      density = .FALSE.

          a1 = DBLE(osms28(i1-1,5))
          a2 = DBLE(osms28(i1-1,6))
          b1 = DBLE(osms28(i1  ,5))
          b2 = DBLE(osms28(i1  ,6))

          DO ik = 1, nks(ir)

            id = korpg(ik,ir)
            c1 = 0.5D0 * DBLE(rvertp(1,id) + rvertp(2,id))
            c2 = 0.5D0 * DBLE(zvertp(1,id) + zvertp(2,id))
            d1 = 0.5D0 * DBLE(rvertp(3,id) + rvertp(4,id))
            d2 = 0.5D0 * DBLE(zvertp(3,id) + zvertp(4,id))
       
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
 
            IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .          tcd.GE.0.0.AND.tcd.LT.1.0) THEN
c...          Intersecting between the line segment and the ring is found:

c...          These are here in case psin0 gets picked up when linking exponential
c             decay data to neighbouring ring:
              IF (ringtype(ir).EQ.PFZ) THEN
                psin0 = -HI
                psin1 =  HI
              ELSE
                psin0 = HI
                psin1 = HI
              ENDIF

              IF (.FALSE.) THEN
c...            Need a check to know if the ring is over-constrained:
              ELSEIF (index.GE.0.AND.index.LE.6) THEN
c              ELSEIF (index.GE.1.AND.index.LE.5) THEN

                IF (index.GE.1.AND.index.LE.5) THEN
                  s(index)=ksb(ik-1,ir)+
     .                     SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
                ENDIF

c...            Find data boundary values -- NEEDS WORK!:
                i2 = i0

                i3 = i1

c...            Flag...
                tetarget = .FALSE.
                IF ((index.EQ.1.OR.index.EQ.5).AND.osms28(i3,7).LT.0.0) 
     .            tetarget = .TRUE.

                IF     (mode.EQ.1.OR.mode.EQ.2.OR.mode.EQ.3.OR.
     .                  mode.EQ.5) THEN
c...              Interpolation boundary values provided in the input file:
c                    WRITE(0,*) 'DATA1:',index,osms28(i2,7),osms28(i2,8)

                  t0 = osms28(i2,7) + frac*(osms28(i2,9) -osms28(i2,7))
                  n0 = osms28(i2,8) + frac*(osms28(i2,10)-osms28(i2,8))
                  t1 = osms28(i3,7) + frac*(osms28(i3,9) -osms28(i3,7))
                  n1 = osms28(i3,8) + frac*(osms28(i3,10)-osms28(i3,8))

                  IF ((index.EQ.3.OR.index.EQ.4.OR.index.EQ.5).AND.
     .                (osms28(i2,7).EQ.-99.0.OR.
     .                 osms28(i2,8).EQ.-99.0)) THEN
c...                Linking to another plasma region that is *already* calculated: 
                    CALL FindCell(i2,i3,ir,ikcell,ircell)
                    IF (ringtype(ir).EQ.PFZ) THEN
                      ik1 = ikouts(ikcell(2),ircell(2))
                      ir1 = irouts(ikcell(2),ircell(2))
c                      WRITE(0,*) 'FOUND YOU PFZ:',ik1,ir1,
c     .                                            density,tetarget,index
                    ELSE
                      ik1 = ikins(ikcell(2),ircell(2))
                      ir1 = irins(ikcell(2),ircell(2))
                    ENDIF
                    IF (osms28(i2,7).EQ.-99.0) t0 = ktebs(ik1,ir1)
                    IF (osms28(i2,8).EQ.-99.0) THEN
                      IF (density) THEN
                        n0 = knbs (ik1,ir1)
c                        WRITE(0,*) 'denisity POWER:',n0
                      ELSE
                        n0 = 2.0 * ktebs(ik1,ir1) * knbs(ik1,ir1)  ! Add Ti and M? 
                      ENDIF
                    ENDIF
c...                Shouldn't really be outer target (all this would go away if PSITARG was
c                   assigned properly):
                    IF (coord.EQ.3) psin0 = psitarg(ir1,1)
c                   IF (coord.EQ.4) rho0 =
c                    WRITE(0,*) 'DATA2:',t0,psin0,ik1,ir1,index
c                    WRITE(0,*) 'DATA2:',n0,psin0,ik1,ir1
c                    STOP 'rgsd'
                  ENDIF

c...              Make sure that t0,1 are positive:
                  IF (tetarget) THEN
                    t0 = ABS(t0)
                    t1 = ABS(t1)
                  ENDIF

                  IF (coord.EQ.1) THEN
c...                Linear along the line segment, nice and simple:
                    val0 = 0.0
                    val1 = 1.0
                    val = SNGL(tab)
c                    val = SNGL(tcd) ! BUG
                  ELSEIF (coord.EQ.5) THEN
c...                Just give me PSIn:
                    psin0 = psitarg(ir,1)
                  ELSE
c...              
                    IF (NINT(osms28(i2,11)).NE.NINT(osms28(i3,11)).OR.
     .                  NINT(osms28(i2,12)).NE.NINT(osms28(i3,12)).OR.
     .                  NINT(osms28(i2,11)).GT.NINT(osms28(i2,12)).OR.
c...                    This is here because PSITARG is currently not assigned in the core... 
     .                  NINT(osms28(i2,11)).LT.irsep)
     .                CALL ER('FindS28Parameters_V3','Invalid '//
     .                        'IR range',*99)

c...                Need range of PSIn over the segment:
                    DO ir1 = NINT(osms28(i2,11)), NINT(osms28(i2,12)) 
                      DO ik1 = 1, nks(ir1)
                        id1 = korpg(ik1,ir1)
                        c1 = 0.5D0 * DBLE(rvertp(1,id1) + rvertp(2,id1))
                        c2 = 0.5D0 * DBLE(zvertp(1,id1) + zvertp(2,id1))
                        d1 = 0.5D0 * DBLE(rvertp(3,id1) + rvertp(4,id1))
                        d2 = 0.5D0 * DBLE(zvertp(3,id1) + zvertp(4,id1))
                        CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
                        IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .                      tcd.GE.0.0.AND.tcd.LT.1.0) THEN
c...                      Should be more careful which PSITARG is used and when...
c                          WRITE(0,*) 'MARK:',ir1,psin0,psin1
                          IF     (ringtype(ir).EQ.SOL1) THEN
                            IF (psin0.EQ.HI) psin0 = psitarg(ir1,1)
                            IF (psin0.NE.HI) psin1 = psitarg(ir1,1)
c                            WRITE(0,*) 'PSIN:',ir1,psin0,psin1
                          ELSEIF (ringtype(ir).EQ.PFZ) THEN
                            psin0 = MAX(psin0,psitarg(ir1,1))
                            psin1 = MIN(psin1,psitarg(ir1,1))
                          ELSE
                            CALL ER('FindS28Parameters_V3','Invalid '//
     .                              'RINGTYPE',*99)
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDDO


c                    WRITE(0,*) 'PSIn:',ir,psin0,psin1

                  
                    IF     (coord.EQ.2) THEN
c...                  Spatial along the IR-range of applicability:
                      STOP 'NOT READY 2'
                    ELSEIF (coord.EQ.3) THEN
c...                  PSIn:
                      val0 = 0.0
                      val1 = ABS(psin1         - psin0)
                      val  = ABS(psitarg(ir,1) - psin0)
                    ELSEIF (coord.EQ.4) THEN
c...                  RHO:
                      STOP 'NOT READY 4'
                    ELSE
                      WRITE(0,*) 'I0,I1=',i0,i1
                      CALL ER('S28params_v3','Invalid COORD A',*99)
                    ENDIF
                  ENDIF

                ELSEIF (mode.EQ.4) THEN
c...              Load probe data, dummy values here:
                  t0 = osms28(i2,7)
                  t1 = t0
                  n0 = osms28(i2,8)
                  n1 = n0

                ELSE
                  CALL ER('S28params_v3','Invalid MODE',*99)   
                ENDIF

c...            Check if quantities should be assigned:
                tc = .TRUE.
                nc = .TRUE.
                IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
                IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

c                IF (ir.GT.22.AND.ir.LT.irwall.AND.index.EQ.1.AND.
c     .              .NOT..FALSE..AND.(.TRUE..OR.index.EQ.4)) THEN
c                  WRITE(0,*) 'PSIn:',ir,psitarg(ir,1),psin0,psin1
c                  WRITE(0,*) 'i2,3:',i2,i3
c                  WRITE(0,*) 'val :',ir,val,val0,val1
c                  WRITE(0,*) 'Te  :',ir,t0,t1
c                  WRITE(0,*) 'ne  :',ir,n0,n1
c                ENDIF

                IF     (mode.EQ.1) THEN
c...              Power law between v1 and v2:
                  val2 = (val - val0) / (val1 - val0)
c                  IF (index.EQ.5) WRITE(0,*) '5:',val2,ir

                  IF (tc) te(index) = t0 + val2**expon * (t1 - t0)
                  IF (nc) ne(index) = n0 + val2**expon * (n1 - n0)

c                  IF (index.EQ.3.AND.ir.EQ.71) THEN
c                    WRITE(0,*) 'A:',val,n0,n1,density
c                  ENDIF
      

                ELSEIF (mode.EQ.2) THEN
c...              Exponential decay between v1 and v2:
                  C = expon  ! BUG!
c                  C = -expon
c                  val = 1.0

                  IF (tc) THEN
                    A = (t1 - t0) / (EXP(-val1 / C) - 1.0)
                    B = t0 - A
                    te(index) = A * EXP(-val / C) + B
                  ENDIF
                  IF (nc) THEN
                    A = (n1 - n0) / (EXP(-val1 / C) - 1.0)
                    B = n0 - A
                    ne(index) = A * EXP(-val / C) + B
                  ENDIF



                ELSEIF (mode.EQ.3) THEN
c...              Exponential decay to infinity:
                  C = expon
                  A = t0 - t1
                  B = t1 
                  IF (tc) te(index) = A * EXP(-val / C) + B
                  A = n0 - n1
                  B = n1
                  IF (nc) ne(index) = A * EXP(-val / C) + B

c                  WRITE(0,*) 'DATA T:',ir,te(index)
c                  IF (index.EQ.3.AND.ir.EQ.71) THEN
c                    WRITE(0,*) 'A:',a,b,c
c                    WRITE(0,*) 'A:',val,A*EXP(-val / C)
c                    WRITE(0,*) 'A:',ne(index)
c                  ENDIF


                ELSEIF (mode.EQ.4) THEN
c...              Load probe data from the .experiments file:
                  IF     (coord.EQ.3) THEN
                    prb1 = -1.0
                  ELSEIF (coord.EQ.4) THEN
                    prb1 = -3.0
                  ELSE
                      WRITE(0,*) 'I0,I1=',i0,i1,coord
                    CALL ER('S28params_v3','Invalid COORD B',*99)   
                  ENDIF

                  IF (tc) THEN
                    tmp1 = prb1
                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,7)),2,tmp1) 
                    te(index) = tmp1
                    WRITE(0,*) 'TE PROBE:',tmp1,ir,psitarg(ir,1)
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY'
                    ELSEIF (expon.EQ.2.0) THEN
c...                  
                      te(0) = -te(index)
                      te(6) = -te(index)
                    ELSEIF (expon.EQ.3.0) THEN
c...                  
                      te(0) = 97.0
                      te(6) = 97.0
                    ENDIF
                  ENDIF

                  IF (nc) THEN
                    tmp1 = prb1
                    CALL LoadProbeDataS28(ir,NINT(osms28(i2,8)),1,tmp1) 
                    ne(index) = tmp1
                    WRITE(0,*) 'NE PROBE:',tmp1,ir,psitarg(ir,1)
                    IF     (expon.EQ.1.0) THEN
                      STOP 'OPTION NOT READY'
                    ELSEIF (expon.EQ.2.0) THEN
c...                  
                      ne(0) = -ne(index)
                      ne(6) = -ne(index)
                    ELSEIF (expon.EQ.3.0) THEN
                    ENDIF
                  ENDIF

c                  WRITE(0,*) 'PROBIN:',prb1,te(index),ne(index),tc,nc

                ELSEIF (mode.EQ.5) THEN
c...              Interpolations from polynomial and exponential fitting parameters 
c                 that are listed in the input file:

                  DO i4 = 1, NINT(osms28(i1,4))

                    ifit = i1 + i4

c                    WRITE(0,*) '  CHECKING!',
c     .                psin1,osms28(ifit,3),osms28(ifit,4)

                    IF     (osms28(ifit,2).EQ.-1.0) THEN 
                      psin1 = psin0 - osms28(ifit,3)

                    ELSEIF (psin1.GE.osms28(ifit,3).AND.
     .                      psin1.LE.osms28(ifit,4)) THEN

c                      WRITE(0,*) '  FOUND!',i4,ifit

c...                  Only need this for the exponential fit coefficients
c                     returned in IDL (ts.pro at the moment):
                      IF (osms28(ifit,2).EQ.2.0) THEN
                        psin2 = psin1 - osms28(ifit,3)
                      ELSE
                        psin2 = psin1
                      ENDIF
 
c                      WRITE(0,*) 'PSIN:',psin1,psin2

                      SELECTCASE (NINT(osms28(ifit,2)))
                        CASE ( 1)  ! Polynomial
                          val = osms28(ifit,7 )            +
     .                          osms28(ifit,8 ) * psin2    +
     .                          osms28(ifit,9 ) * psin2**2 +
     .                          osms28(ifit,10) * psin2**3 +
     .                          osms28(ifit,11) * psin2**4 +
     .                          osms28(ifit,12) * psin2**5

                     WRITE(0,*) 'COEF:',osms28(ifit,7:12)
                        CASE ( 2)  ! Exponential
                          val = osms28(ifit,7) * 
     .                          EXP(osms28(ifit,8)*psin2) + 
     .                          osms28(ifit,9)
                        CASE ( 3)  ! TANH 
                          p(0:4) = osms28(ifit,7:11)

                          v0 = (p(0) - psin2) / (2.0 * p(1))               ! from /home/mastts/lib/edgefunctionats.pro
                          v1 = ((1. + p(3) * v0) * EXP(v0) - EXP(-v0)) /   ! from /home/mastts/bck/mtanh.pro (right one?) 
     .                                            (EXP(v0) + EXP(-v0))
                          val = (p(2) - p(4)) / 2.0 * (v1 + 1.0) + p(4)

                        CASEDEFAULT
                          CALL ER('Find...','Unknown fit type',*99)
                      ENDSELECT

                      SELECTCASE (NINT(osms28(ifit,6)))
                        CASE (1)  
                          IF (osms28(ifit,2).EQ.3.0) val = val * 1.0E+19  ! Special for TANH fit from ts.pro
                          ne(index) = val
                        CASE (4)  
                          te(index) = val 

        WRITE(SLOUT,'(A,2I6,2F10.3,I6,1P,E10.2,0P,F10.2)')
     .   'FITS:',
     .  ir,NINT(osms28(ifit,2)),psin1,psin2,index,ne(index),te(index)

                        CASEDEFAULT
                          CALL ER('Find...','Unknown data type',*99)
                      ENDSELECT
             
                       

                    ENDIF

                  ENDDO

c                  WRITE(0,*) 'MODE=5:',ir,psin1,i1
c                  WRITE(0,*) '      :',index,te(index),ne(index)
 
                ELSE
                  CALL ER('S28params_v3','Invalid MODE',*99)   
                ENDIF

              ELSE
                CALL ER('FindSOL28Parameters','Invalid parameter '//
     .                  'index',*99)
              ENDIF


              IF (tetarget) THEN
c                WRITE(0,*) 'TETARGET:',tetarget,te(index),i3,index
                te(index) = -te(index)
              ENDIF

              IF (.NOT.density.AND.index.EQ.3.AND.nc.AND.mode.NE.4) THEN
c...            Convert the pressure value to density:
                IF (te(3).EQ.0.0) THEN
                  CALL ER('FindSOL28Parameters','Te3 not assigned',*99)
                ELSE
c...              Assumes that the Mach no. is low (note: Ti.NE.Te is not a problem
c                 since that situation is currently resolved in the SOL28 routine):
c                  WRITE(0,*) 'CONVERTING PRESSURE TO DENSITY',ne(3)
c                  WRITE(0,*) 'OLD D:',ne(3)
                  ne(3) = ne(3) / (2.0 * te(3))
c                  WRITE(0,*) 'NEW D:',ne(3)
                ENDIF
              ENDIF

            ENDIF

          ENDDO

c          IF (.NOT.density.AND.index.EQ.3.AND.nc.AND.mode.NE.4) THEN
cc...        Convert the pressure value to density:
c            IF (te(3).EQ.0.0) THEN
c              CALL ER('FindSOL28Parameters','Te(3) not assigned',*99)
c            ELSE
cc...          Assumes that the Mach no. is low (note: Ti.NE.Te is not a problem
cc             since that situation is currently resolved in the SOL28 routine):
cc              WRITE(0,*) 'CONVERTING PRESSURE TO DENSITY',ne(3)
c              WRITE(0,*) 'OLD D:',ne(3)
c              ne(3) = ne(3) / (2.0 * te(3))
c              WRITE(0,*) 'NEW D:',ne(3)
c            ENDIF
c          ENDIF

c...    End of OSMS28(I1,1).EQ.98.0 block:      
        ENDIF

      ENDDO




c          WRITE(0,*) 'NE3=',ne(3)
      IF (.TRUE.) THEN
c...    Specify te(1) from target data.  The -ve is to trigger
c       te(0)=te(1):
        IF (te(1).EQ.-98.0) te(1) = -kteds(idds(ir,2))

c...    Specify te(5) from target data.  The -ve is to trigger
c       te(6)=te(5):
        IF (te(5).EQ.-98.0) te(5) = -kteds(idds(ir,1))
      ENDIF 


     

      IF (.NOT..TRUE.) THEN


        IF (te(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
          te(3) = tetmp1
c...      Don't necessarily want both targets done...:
          te(0) = 97.0
          te(6) = 97.0
        ENDIF
        IF (ne(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
          ne(3) = tetmp1
        ENDIF


        IF (te(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          te(3) = tetmp1
          te(0) = -te(3)
          te(6) = -te(3)
        ENDIF
        IF (ne(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          ne(3) = tetmp1
          ne(0) = -ne(3)
          ne(6) = -ne(3)
        ENDIF

c        WRITE(0,*) 'NE(3):',ne(3),te(3)

        IF (te(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,2,tetmp1) 
          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c         te(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
        ENDIF
        IF (ne(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,1,tetmp1) 
          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          WRITE(PINOUT,*) 'COMPROMISE:',ir,psitarg(ir,2),ne(3),te(3)
        ENDIF

        IF (te(3).LE.-1.0) CALL LoadProbeDataS28(ir,s28probe,2,te(3))
c *HACKISH*
        IF (ne(3).LE.-1.0.AND.ne(3).GT.-100.0) 
     .    CALL LoadProbeDataS28(ir,s28probe,1,ne(3))

        IF (ne(3).LT.-100.0) THEN
c...      Pressure specified, not density.  Calculate ne(3) from the pressure
c         and the specified te(3) value, assuming Ti=Te:
          ne(3) = -ne(3) / (2.0 * te(3))
        ENDIF

c...    Specify te(1) from target data.  The -ve is to trigger
c       te(0)=te(1):
        IF (te(1).EQ.-99.0) te(1) = -kteds(idds(ir,2))

c...    Specify te(5) from target data.  The -ve is to trigger
c       te(6)=te(5):
        IF (te(5).EQ.-99.0) te(5) = -kteds(idds(ir,1))

      ELSE
      ENDIF


      RETURN
99    CONTINUE
      WRITE(0,*) 'IK,IR=',ik,ir,i1,osms28(i1,1)
      WRITE(0,*) 'TAB,TCD=',tab,tcd
      WRITE(0,*) 'MODE= ',mode
      WRITE(0,*) 'COORD=',coord
      STOP
      END
c
c
c ======================================================================
c
c sburoutine: CalcRadialDrift
c
c 
      SUBROUTINE CalcRadialDrift(ir1)
      USE mod_sol28_old
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_slcom
      use mod_slout
      IMPLICIT none

      INTEGER ir1

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'slcom'
c     INCLUDE 'slout'

      INTEGER ik,ir,ik2,ir2,irs,ire
      REAL    t1,t2,frac,tgrad,area,vperp,cfpflx,tmpflx1,tmpflx2



c...  Get the new cross-field flux:
      IF (ir1.EQ.-2) THEN
c *HACK* MOVE?

        DO ir = irsep, irwall-1

c          IF ((ringtype(ir).EQ.SOL1.AND.s28cfpdrft   .EQ.-3).OR.
c     .        (ringtype(ir).EQ.PFZ .AND.s28cfppfzdrft.EQ.-3)) THEN
          IF ((ringtype(ir).EQ.SOL1.AND.
     .         (s28cfpdrft   .EQ.2.OR.s28cfpdrft   .EQ.6)).OR.
     .        (ringtype(ir).EQ.PFZ .AND.
     .         (s28cfppfzdrft.EQ.2.OR.s28cfppfzdrft.EQ.6))) THEN

            DO ik = 1, nks(ir)

              ir2 = 0
c              tmpflx1 = tmpflx(ik,ir)
              tmpflx1 = osmcfpflx(ik,ir,2)
              tmpflx2 = 0.0

              IF (ik.LE.nks(ir)/2) THEN
                IF (ir.EQ.irsep) THEN
                ELSE
                  ir2 = irins(ik,ir)       ! need to do something with tgrad to know 
                  DO ik2 = 1, nks(ir2)-1   ! which way the cf wind is blowing... 
                    IF (thetag(ik,ir).GT.thetag(ik2  ,ir2).AND.  ! ignoring region near targets... 
     .                  thetag(ik,ir).LE.thetag(ik2+1,ir2)) THEN ! need a way of flagging errors in THETAG...
                      frac = (thetag(ik   ,ir ) - thetag(ik2,ir2)) /
     .                       (thetag(ik2+1,ir2) - thetag(ik2,ir2))                    
c                      tmpflx2 = (1.0 - frac) * tmpflx(ik2  ,ir2) +
c     .                                 frac  * tmpflx(ik2+1,ir2)
                      tmpflx2 = (1.0 - frac) * osmcfpflx(ik2  ,ir2,2) +
     .                                 frac  * osmcfpflx(ik2+1,ir2,2)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF

              ELSE

                IF (ir.EQ.irwall-1) THEN
                ELSE
                  ir2 = irouts(ik,ir) ! need to do something with tgrad to know which way the cf wind is blowing...
                  DO ik2 = 1, nks(ir2)-1 
                    IF (thetag(ik,ir).GT.thetag(ik2  ,ir2).AND.  ! ignoring region near targets... 
     .                  thetag(ik,ir).LE.thetag(ik2+1,ir2)) THEN ! need a way of flagging errors in THETAG...
                      frac = (thetag(ik   ,ir ) - thetag(ik2,ir2)) /
     .                       (thetag(ik2+1,ir2) - thetag(ik2,ir2))                    
c                      tmpflx2 = (1.0 - frac) * tmpflx(ik2  ,ir2) +
c     .                                 frac  * tmpflx(ik2+1,ir2)
                      tmpflx2 = (1.0 - frac) * osmcfpflx(ik2  ,ir2,2) +
     .                                 frac  * osmcfpflx(ik2+1,ir2,2)
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF  

              ENDIF

              osmcfpflx(ik,ir,1) = tmpflx1 - tmpflx2

c              IF (ir.Eq.irwall-1) 
c     .           WRITE(0,'(A,4I6,1P,3E10.2)') 
c     .           'TMPFLX:',ik,ir,ik2,ir2,tmpflx1,tmpflx2,
c     .           osmcfpflx(ik,ir,1)

            ENDDO

          ELSE
c...
            DO ik = 1, nks(ir)
              osmcfpflx(ik,ir,1) = osmcfpflx(ik,ir,2)
            ENDDO

          ENDIF
        ENDDO

c * MESSY*
        RETURN

      ELSEIF (ir1.EQ.-1) THEN
        irs = irsep
        ire = irwall-1
      ELSEIF (ir1.GE.irsep.AND.ir1.LT.irwall) THEN
        irs = ir1
        ire = ir1
      ELSE
        CALL ER('CalcRadialDrift','Invalid parameter',*99)
      ENDIF

      DO ir = irs, ire

        DO ik = nks(ir), 2, -1


c        DO ik = nks(ir), nks(ir)/2, -1

C        IF (zs(ik,ir).GT.zxp) EXIT

          IF (.TRUE.) THEN
c...        Temperature gradient calculated along center of cell:          

            IF     (ik.EQ.nks(ir)) THEN
              t1 = kteds(idds(ir,1))
            ELSEIF (ik.EQ.1) THEN
              t1 = kteds(idds(ir,2))
            ELSE
              frac = (ksb(ik  ,ir) - kss(ik,ir)) / 
     .               (kss(ik+1,ir) - kss(ik,ir))
              t1 = ktebs(ik,ir) * (1.0 - frac) + ktebs(ik+1,ir) * frac
            ENDIF

            frac = (ksb(ik-1,ir) - kss(ik-1,ir)) / 
     .             (kss(ik  ,ir) - kss(ik-1,ir))
            t2 = ktebs(ik-1,ir) * (1.0 - frac) + ktebs(ik,ir) * frac

            tgrad = (t1 - t2) / (kpb(ik,ir) - kpb(ik-1,ir))
 
            vperp = 0.71 * tgrad / 5.0 

c...        Approximation (slight):
            area = 2.0 * PI * rs(ik,ir) * (kpb(ik,ir) - kpb(ik-1,ir)) * 
     .             eirtorfrac

            cfpflx = -ABS(vperp * knbs(ik,ir) * area / kvols(ik,ir))

            IF ((ringtype(ir).EQ.SOL1.AND.
     .           (s28cfpdrft   .EQ.-2.OR.s28cfpdrft   .EQ.-3.OR.
     .            s28cfpdrft   .EQ.-2)).OR.
     .          (ringtype(ir).EQ.PFZ .AND.
     .           (s28cfppfzdrft.EQ.-2.OR.s28cfppfzdrft.EQ.-3.OR.
     .            s28cfppfzdrft.EQ.-2))) THEN

              osmcfpflx(ik,ir,1) = 0.0
              osmcfpflx(ik,ir,2) = 0.0
              tmpflx(ik,ir) = cfpflx

            ELSE
              osmcfpflx(ik,ir,2) = (1.0-rel_frac) * osmcfpflx(ik,ir,2) + 
     .                                  rel_frac  * cfpflx
            ENDIF

            IF (ir.EQ.irsep.AND.ik.EQ.nks(irsep))  
     .        WRITE(0,*) 'FLUX A:',osmcfpflx(nks(irsep),irsep,2)

c            WRITE(0,'(A,2I6,1P,4E10.2,0P,4X,3F10.2)') 
c     .        'FLUX:',ik,ir,vperp,tgrad,area,cfpflx,t1,t2,
c     .        ktebs(ik,ir)
          ELSE
            CALL ER('CalcRadialDrift','Invalid option',*99)
          ENDIF

        ENDDO

      ENDDO


      RETURN
99    STOP
      END
c
c ======================================================================
c
      SUBROUTINE CalcRecLineEmissions(ir)
      use mod_params
      use mod_cgeom
      use mod_slcom
      IMPLICIT none

      INTEGER ir

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'

      REAL GetEAD

      INTEGER ik
      REAL    pop,
     .        fac32,fac42,fac52,
     .        fac21,fac31,fac41,fac51

c...  Einstein coefficients for Balmer lines (from HALPHA.F 
c     in EIRENE):
      FAC32=4.410E+07
      FAC42=8.419E+06
      FAC52=2.530E+06

      FAC21=4.699E+08
      FAC31=5.575E+07
      FAC41=1.278E+07
      FAC51=4.125E+06

c...  Dgamma:
      DO ik = 1, nks(ir)
        pop = FAC52 * GetEAD(ktebs(ik,ir),knbs(ik,ir),20,'H.12')
        
        pinline(ik,ir,6,H_BGAMMA) = pop * knbs(ik,ir) 
c...    Unit conversion to W m-3:
c     .     (6.63E-34 * 3.0E+08) / (4340.0 * 1.0E-10)

c        WRITE(0,*) 'IK,IR,EM',ik,ir,pinline(ik,ir,6,H_BGAMMA),pop/FAC52
      ENDDO


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadProbeDataS28(ir,index,mode,val)
      use mod_params
      use mod_cgeom
      use mod_slcom
      IMPLICIT none

      INTEGER ir,index,mode
      REAL    val

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'

      INTEGER    DATAUNIT   ,MAXTDAT     ,MAXCOLS
      PARAMETER (DATAUNIT=13,MAXTDAT=1000,MAXCOLS=10)

      INTEGER   etype,ndata,ncol,iblock,i1,i2,lastindex
      REAL      xdata(MAXTDAT),edata(MAXTDAT,MAXCOLS),xval
      CHARACTER datatitle*128

      DATA lastindex /-1/

      SAVE

      IF (index.NE.lastindex) THEN
        iblock = index

        CALL Load_Expt_Data(DATAUNIT,iblock,xdata,etype,edata,
     .                      MAXCOLS,MAXTDAT,ndata,ncol,datatitle)

c        WRITE(0,*) '============================================='
c        WRITE(0,*) 'TITLE  = ',datatitle(1:LEN_TRIM(datatitle))
c        WRITE(0,*) 'INDEX  = ',iblock
c        WRITE(0,*) 'TYPE   = ',etype
c        WRITE(0,*) 'NDATA  = ',ndata
c        WRITE(0,*) 'NCOL   = ',ncol

c        DO i1 = 1, ndata
c          WRITE(0,*) i1,xdata(i1),(edata(i1,i2),i2=1,ncol)
c        ENDDO

        lastindex = index
      ENDIF

      IF     (val.EQ.-1.0) THEN
        xval = psitarg(ir,1)
      ELSEIF (val.EQ.-2.0) THEN
c        xval = rho(ir,CELL1) 
c *HARDCODED*
        xval = rho(ir,CELL1) - 0.002
      ELSEIF (val.EQ.-3.0) THEN
        xval = rho(ir,CELL1) 
      ELSE
        CALL ER('LoadProbeDataS28','Unrecognized mapping',*99)
      ENDIF

      IF     (xval.LT.xdata(1)) THEN
        WRITE(0,*) 'WARNING:  INTERPOLATION FAILED, X-DATA BEYOND RANGE'
        val = edata(1,mode)
      ELSEIF (xval.GT.xdata(ndata)) THEN
        WRITE(0,*) 'WARNING:  INTERPOLATION FAILED, X-DATA BEYOND RANGE'
        val = edata(ndata,mode)
      ELSE
        CALL Fitter(ndata,xdata,edata(1,mode),1,xval,val,'LINEAR')
      ENDIF

c      WRITE(0,*) ir,xval,val

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: NormalizeInitSrc
c
      SUBROUTINE NormalizeInitSrc(region,ir,source,flux,tag)
      use mod_params
      use mod_cgeom
      use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'

      INTEGER region,ir,ik
      REAL    source(MAXNKS,MAXNRS)
      REAL*8  flux,integral
      CHARACTER tag*(*)

      REAL    frac

      CALL DB('Normalizing PIN source estimates')

       

      IF     (region.EQ.IKLO) THEN
        CALL CalcIntegral2(source(1,ir),ikbound(ir,IKLO),
     .                     osm_sympt(ir),ir,integral,6)
c...GOG
c        IF (tag.EQ.'mp3') THEN
c          DO ik = 1, osm_sympt(ir)
c            WRITE(SLOUT,*) '-->',source(ik,ir)
c          ENDDO
c        ENDIF
      ELSEIF (region.EQ.IKHI) THEN
        CALL CalcIntegral2(source(1,ir),osm_sympt(ir)+1,
     .                     ikbound(ir,IKHI),ir,integral,7)
      ELSEIF (region.EQ.3) THEN
        CALL CalcIntegral2(source(1,ir),ikbound(ir,IKLO),
     .                     ikbound(ir,IKHI),ir,integral,2)
      ELSEIF (region.EQ.4) THEN
        CALL CalcIntegral2(source(1,ir),1,nks(ir),
     .                     ir,integral,2)
      ELSE
        CALL ER('NormalizeInitSrc','Invalid region',*99)
      ENDIF


      IF (integral.NE.0.0) THEN
c THIS IS SCARY!
        IF (tag.EQ.'mp3'.OR.tag.EQ.'mp4'.OR.tag.EQ.'mp5'.OR.
     .      tag.EQ.'cf5'.OR.tag.EQ.'cf6'.OR.tag.EQ.'mp6') THEN
          frac = flux / integral
        ELSE
          frac = DABS(flux / integral)
        ENDIF

c        WRITE(SLOUT,'(A,I4,1P,2E10.2,0P,F12.4,F10.2,A)')
c     .    'NORM SRC  IR FLUX INT FRAC %DIFF =',
c     .    ir,flux,integral,frac,
c     .    (DABS(flux) - DABS(integral))/DABS(flux)*100.0,'% '//tag

        IF (region.EQ.IKLO) THEN
          DO ik = ikbound(ir,IKLO), osm_sympt(ir)
            source(ik,ir) = source(ik,ir) * frac
          ENDDO
        ELSEIF (region.EQ.IKHI) THEN
          DO ik = osm_sympt(ir)+1, ikbound(ir,IKHI)
            source(ik,ir) = source(ik,ir) * frac
          ENDDO
        ELSEIF (region.EQ.3) THEN
          DO ik = ikbound(ir,IKLO), ikbound(ir,IKHI)
            source(ik,ir) = source(ik,ir) * frac
          ENDDO
        ELSEIF (region.EQ.4) THEN
          DO ik = 1, nks(ir)
            source(ik,ir) = source(ik,ir) * frac
          ENDDO
        ENDIF


        IF     (region.EQ.IKLO) THEN
          CALL CalcIntegral2(source(1,ir),ikbound(ir,IKLO),
     .                       osm_sympt(ir),ir,integral,6)
        ELSEIF (region.EQ.IKHI) THEN
          CALL CalcIntegral2(source(1,ir),osm_sympt(ir)+1,
     .                       ikbound(ir,IKHI),ir,integral,7)
        ELSEIF (region.EQ.3) THEN
          CALL CalcIntegral2(source(1,ir),ikbound(ir,IKLO),
     .                       ikbound(ir,IKHI),ir,integral,2)
        ELSEIF (region.EQ.4) THEN
          CALL CalcIntegral2(source(1,ir),1,nks(ir),
     .                       ir,integral,2)
        ENDIF


        IF (DABS(DABS(integral)-DABS(flux))/DABS(flux).GT.1.0E-6.AND.
     .      .NOT.(DABS(integral).EQ.0.0D0.AND.DABS(flux).EQ.0.0)) THEN
          WRITE(PINOUT,*) 'PROBLEM! ',flux,integral,tag
c          WRITE(0     ,*) 'PROBLEM! ',flux,integral,tag
        ENDIF

      ELSE
        WRITE(PINOUT,*) 'Integral is zero'
      ENDIF


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: Normalize
c
      SUBROUTINE Normalize(iks,ike,ir,source,flux,mode,tag)
      use mod_params
      use mod_cgeom
      use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'

      INTEGER iks,ike,ir,ik,mode
      REAL    source(MAXNKS,MAXNRS),frac
      REAL*8  flux,integral
      CHARACTER tag*(*)

      CALL DB('Start of Normalize')


      CALL CalcIntegral2(source(1,ir),iks,ike,ir,integral,4)
       
      IF (integral.NE.0.0) THEN

        IF     (mode.EQ.0) THEN
          frac = flux / integral
        ELSEIF (mode.EQ.1) THEN
          frac = DABS(flux / integral)
        ELSE
          CALL ER('Normalize','Invalid MODE specification',*99)
        ENDIF

        DO ik = iks, ike
          source(ik,ir) = source(ik,ir) * frac
c          WRITE(0,*) 'NORMAL:',ik,ir,frac,source(ik,ir)
        ENDDO


c...    Check:
        CALL CalcIntegral2(source(1,ir),iks,ike,ir,integral,4)
        IF (DABS(DABS(integral)-DABS(flux))/DABS(flux).GT.1.0E-6) THEN
          WRITE(PINOUT,*) 'NORMALIZE PROBLEM: ',flux,integral,tag
c          WRITE(0     ,*) 'NORMALIZE PROBLEM: ',flux,integral,tag
        ENDIF

      ELSE
        WRITE(PINOUT,*) 'Integral is zero'
      ENDIF


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: AdjustDensity
c
      SUBROUTINE AdjustDensity(n)
      use mod_params
      use mod_slcom
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'slcom'

      REAL n

      IF (ABS(muldensity).GE.1.0E+12) THEN
        n = n + muldensity
      ELSEIF (muldensity.NE.1.0) THEN
        n = n * muldensity
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
      SUBROUTINE FindS28Parameters(ir,te,ne,nf,s,new)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      IMPLICIT none
       
      INTEGER ir
      LOGICAL new
      REAL    te(0:6),ne(0:6),nf,s(0:6),isat(0:6)

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'slcom'


      INTEGER i0,i1,i2,ik,id
      LOGICAL tc,nc
      REAL    frac,t0,t1,n0,n1,A,B,C,tetmp1,tetmp2
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd

      s (0) = 0.0
      ne(1) = 0.0
      te(1) = 0.0
      te(2) = 0.0
      te(3) = 0.0
      te(4) = 0.0
      te(5) = 0.0

      frac = 0.0
      IF (new) THEN
        IF (rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN
          frac = MAX(0.0,REAL(rel_step - 1) / REAL(rel_nstep - 1))
          frac = rel_bound1 + frac * (rel_bound2 - rel_bound1) 
        ENDIF
      ELSE
        frac = MAX(0.0,REAL(rel_step - 1) / REAL(rel_nstep - 1))
        frac = rel_bound1 + frac * (rel_bound2 - rel_bound1) 
      ENDIF  
    
      IF (new) THEN
c        WRITE(0,*) 'FRAC:',frac
      ENDIF

      DO i1 = 2, osmns28
        i0 = i1 - 1
        IF (osms28(i0,1).NE.osms28(i1,1).AND.osms28(i1,2).NE.99.0) CYCLE

c...    Do not apply data if IR is outside specified range:
        IF ((osms28(i1,9) .NE.0.0.AND.REAL(ir).LT.osms28(i1,9) ).OR.
     .      (osms28(i1,10).NE.0.0.AND.REAL(ir).GT.osms28(i1,10)))
     .    CYCLE

        IF (osms28(i1,2).EQ.99.0) THEN

          IF (ir.EQ.NINT(osms28(i1,3))) THEN
            i2 = NINT(osms28(i1,1)) - 1
            IF (osms28(i1,5).NE.0.0) te(i2) = osms28(i1,5)
            IF (osms28(i1,6).NE.0.0) ne(i2) = osms28(i1,6)

           WRITE(0,*) 'OVERRIDE:',ir,i2,te(i2),ne(i2)
c           STOP 'sdfsdfsd'

          ENDIF

        ELSE

          a1 = DBLE(osms28(i1-1,3))
          a2 = DBLE(osms28(i1-1,4))
          b1 = DBLE(osms28(i1  ,3))
          b2 = DBLE(osms28(i1  ,4))
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            c1 = 0.5 * DBLE(rvertp(1,id) + rvertp(2,id))
            c2 = 0.5 * DBLE(zvertp(1,id) + zvertp(2,id))
            d1 = 0.5 * DBLE(rvertp(3,id) + rvertp(4,id))
            d2 = 0.5 * DBLE(zvertp(3,id) + zvertp(4,id))
       
            CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
 
            IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .          tcd.GE.0.0.AND.tcd.LT.1.0) THEN

              IF     (osms28(i1,1).EQ.1.0) THEN
                IF (te(1).EQ.0.0) THEN
                  s(0) =ksb(ik-1,ir)+SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
                  n0 = osms28(i0,6) + frac*(osms28(i0,8)-osms28(i0,6))
                  n1 = osms28(i1,6) + frac*(osms28(i1,8)-osms28(i1,6))
                  nf = n0 + SNGL(tab)**osms28(i1,2) * (n1 - n0)
                ELSE
                  CALL ER('FindSOL28Parameters','Type 1 parameter '//
     .                    'over-specified',*99)
                ENDIF
              ELSEIF (osms28(i1,1).EQ.2.0) THEN
                IF (te(1).EQ.0.0) THEN

                  s(1) =ksb(ik-1,ir)+SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
                  t0 = osms28(i0,5) + frac*(osms28(i0,7)-osms28(i0,5))
                  t1 = osms28(i1,5) + frac*(osms28(i1,7)-osms28(i1,5))
                  n0 = osms28(i0,6) + frac*(osms28(i0,8)-osms28(i0,6))
                  n1 = osms28(i1,6) + frac*(osms28(i1,8)-osms28(i1,6))
c...              Check if quantities should be assigned:
                  tc = .TRUE.
                  nc = .TRUE.
                  IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
                  IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

                  IF (osms28(i1,2).LT.0.0) THEN
c...                Exponential decay:
                    C = -osms28(i1,2)
                    A = (t1 - t0) / (EXP(-1.0 / C) - 1.0)
                    B = t0 - A
                    IF (tc) te(1) = A * EXP(-SNGL(tab) / C) + B

                    A = (n1 - n0) / (EXP(-1.0 / C) - 1.0)
                    B = n0 - A
                    IF (nc) ne(1) = A * EXP(-SNGL(tab) / C) + B
                  ELSE
                    IF (tc) te(1) = t0 + SNGL(tab)**osms28(i1,2)*(t1-t0)
                    IF (nc) ne(1) = n0 + SNGL(tab)**osms28(i1,2)*(n1-n0)
                  ENDIF

c                  te(1) = t0 + SNGL(tab)**osms28(i1,2) * (t1 - t0)
c                  ne(1) = n0 + SNGL(tab)**osms28(i1,2) * (n1 - n0)


c                  WRITE(0,*) 'FRAC:',tab,osms28(i1,2),
c     .              SNGL(tab)**osms28(i1,2)
c                   STOP 'sdfsd'
                ELSE
                  CALL ER('FindSOL28Parameters','Type 2 parameter '//
     .                    'over-specified',*99)
                ENDIF
              ELSEIF (osms28(i1,1).EQ.3.0) THEN
                IF (te(2).EQ.0.0) THEN
                  s(2) =ksb(ik-1,ir)+SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
 
                  t0 = osms28(i0,5) + frac*(osms28(i0,7)-osms28(i0,5))
                  t1 = osms28(i1,5) + frac*(osms28(i1,7)-osms28(i1,5))

c...              Check if quantities should be assigned:
                  tc = .TRUE.
c                  nc = .TRUE.
                  IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
c                  IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

                  IF (osms28(i1,2).LT.0.0) THEN
c...                Exponential decay:
                    C = -osms28(i1,2)
                    A = (t1 - t0) / (EXP(-1.0 / C) - 1.0)
                    B = t0 - A
                    IF (tc) te(2) = A * EXP(-SNGL(tab) / C) + B
c                    A = (n1 - n0) / (EXP(-1.0 / C) - 1.0)
c                    B = n0 - A
c                    IF (nc) ne(2) = A * EXP(-SNGL(tab) / C) + B
                  ELSE
                    IF (tc) te(2) = t0 + SNGL(tab)**osms28(i1,2)*(t1-t0)
c                    IF (nc) ne(2) = n0 + SNGL(tab)**osms28(i1,2) * (n1-n0)
                  ENDIF

c                  te(2) = t0 + SNGL(tab)**osms28(i1,2) * (t1 - t0)

                ELSE
                  CALL ER('FindSOL28Parameters','Type 3 parameter '//
     .                    'over-specified',*99)
                ENDIF
              ELSEIF (osms28(i1,1).EQ.4.0) THEN
                IF (.TRUE..OR.te(3).EQ.0.0) THEN
                  s(3) =ksb(ik-1,ir)+SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))
                  t0 = osms28(i0,5) + frac*(osms28(i0,7)-osms28(i0,5))
                  t1 = osms28(i1,5) + frac*(osms28(i1,7)-osms28(i1,5))
                  n0 = osms28(i0,6) + frac*(osms28(i0,8)-osms28(i0,6))
                  n1 = osms28(i1,6) + frac*(osms28(i1,8)-osms28(i1,6))
c...              Check if quantities should be assigned:
                  tc = .TRUE.
                  nc = .TRUE.
                  IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
                  IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

                  IF (osms28(i1,2).LT.0.0) THEN
c...                Exponential decay:
                    C = -osms28(i1,2)
                    A = (t1 - t0) / (EXP(-1.0 / C) - 1.0)
                    B = t0 - A
                    IF (tc) te(3) = A * EXP(-SNGL(tab) / C) + B

                    A = (n1 - n0) / (EXP(-1.0 / C) - 1.0)
                    B = n0 - A
                    IF (nc) ne(3) = A * EXP(-SNGL(tab) / C) + B
                  ELSE
                    IF (tc) te(3) = t0 + SNGL(tab)**osms28(i1,2)*(t1-t0)
                    IF (nc) ne(3) = n0 + SNGL(tab)**osms28(i1,2)*(n1-n0)
                  ENDIF
                ELSE
                  CALL ER('FindSOL28Parameters','Type 4 parameter '//
     .                    'over-specified',*99)
                ENDIF

              ELSEIF (osms28(i1,1).EQ.5.0) THEN
                IF (s(4).EQ.0.0) THEN
                  s(4) =ksb(ik-1,ir)+SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))

                  n0 = osms28(i0,6) + frac*(osms28(i0,8)-osms28(i0,6))
                  n1 = osms28(i1,6) + frac*(osms28(i1,8)-osms28(i1,6))
                  t0 = osms28(i0,5) + frac*(osms28(i0,7)-osms28(i0,5))
                  t1 = osms28(i1,5) + frac*(osms28(i1,7)-osms28(i1,5))

c...              Check if quantities should be assigned:
                  tc = .TRUE.
                  nc = .TRUE.
                  IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
                  IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

                  IF (osms28(i1,2).LT.0.0) THEN
c...                Exponential decay:
                    C = -osms28(i1,2)
                    A = (t1 - t0) / (EXP(-1.0 / C) - 1.0)
                    B = t0 - A
                    IF (tc) te(4) = A * EXP(-SNGL(tab) / C) + B

                  IF (ir.EQ.68) THEN
                    WRITE(0,*) 'A:',a,b,c
                    WRITE(0,*) 'A:',SNGL(tab),A*EXP(-SNGL(tab) / C)
                    WRITE(0,*) 'A:',te(4)
                  ENDIF

                    A = (n1 - n0) / (EXP(-1.0 / C) - 1.0)
                    B = n0 - A
                    IF (nc) ne(4) = A * EXP(-SNGL(tab) / C) + B
                  ELSE
                    IF (tc) te(4) = t0 + SNGL(tab)**osms28(i1,2)*(t1-t0)
                    IF (nc) ne(4) = n0 + SNGL(tab)**osms28(i1,2)*(n1-n0)
                  ENDIF

c                  ne(4) = n0 + SNGL(tab)**osms28(i1,2) * (n1 - n0)
c                  te(4) = t0 + SNGL(tab)**osms28(i1,2) * (t1 - t0)
                ELSE
                  CALL ER('FindSOL28Parameters','Type 2 parameter '//
     .                    'over-specified',*99)
                ENDIF

              ELSEIF (osms28(i1,1).EQ.6.0) THEN

                IF (s(5).EQ.0.0) THEN

                  s(5) =ksb(ik-1,ir)+SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))

                  n0 = osms28(i0,6) + frac*(osms28(i0,8)-osms28(i0,6))
                  n1 = osms28(i1,6) + frac*(osms28(i1,8)-osms28(i1,6))
                  t0 = osms28(i0,5) + frac*(osms28(i0,7)-osms28(i0,5))
                  t1 = osms28(i1,5) + frac*(osms28(i1,7)-osms28(i1,5))

c...              Check if quantities should be assigned:
                  tc = .TRUE.
                  nc = .TRUE.
                  IF (t0.EQ.0.0.OR.t1.EQ.0.0) tc = .FALSE.
                  IF (n0.EQ.0.0.OR.n1.EQ.0.0) nc = .FALSE.

                  IF (osms28(i1,2).LT.0.0) THEN
c...                Exponential decay:
                    C = -osms28(i1,2)
                    A = (t1 - t0) / (EXP(-1.0 / C) - 1.0)
                    B = t0 - A
                    IF (tc) te(5) = A * EXP(-SNGL(tab) / C) + B

                    A = (n1 - n0) / (EXP(-1.0 / C) - 1.0)
                    B = n0 - A
                    IF (nc) ne(5) = A * EXP(-SNGL(tab) / C) + B
                  ELSE
                    IF (tc) te(5) = t0 + SNGL(tab)**osms28(i1,2)*(t1-t0)
                    IF (nc) ne(5) = n0 + SNGL(tab)**osms28(i1,2)*(n1-n0)
                  ENDIF

                ELSE
                  CALL ER('FindSOL28Parameters','Type 6 parameter '//
     .                    'over-specified',*99)
                ENDIF

              ELSEIF (osms28(i1,1).EQ.10.0) THEN
              ELSEIF (osms28(i1,1).EQ.11.0) THEN
              ELSEIF (osms28(i1,1).EQ.12.0) THEN
              ELSE
                CALL ER('FindSOL28Parameters','Invalid parameter '//
     .                  'type',*99)
              ENDIF
            ENDIF

          ENDDO

c...    End of OSMS28(I1,1).EQ.98.0 block:      
        ENDIF

      ENDDO


      IF (new) THEN


        IF (te(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
          te(3) = tetmp1
c...      Don't necessarily want both targets done...:
          te(0) = 97.0
          te(6) = 97.0
        ENDIF
        IF (ne(3).EQ.-97.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
          ne(3) = tetmp1
        ENDIF


        IF (te(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,2,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          te(3) = tetmp1
          te(0) = -te(3)
          te(6) = -te(3)
        ENDIF
        IF (ne(3).EQ.-98.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,s28probe,1,tetmp1) 
c          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          ne(3) = tetmp1
          ne(0) = -ne(3)
          ne(6) = -ne(3)
        ENDIF

c        WRITE(0,*) 'NE(3):',ne(3),te(3)

        IF (te(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,2,tetmp1) 
          CALL LoadProbeDataS28(ir,77,2,tetmp2) 
c         te(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          te(3) = 0.50 * tetmp1 + 0.50 * tetmp2
        ENDIF
        IF (ne(3).EQ.-99.0) THEN
c *HARDCODED*
          tetmp1 = -1.0
          tetmp2 = -1.0
          CALL LoadProbeDataS28(ir,70,1,tetmp1) 
          CALL LoadProbeDataS28(ir,77,1,tetmp2) 
c          ne(3) = 1.00 * tetmp1 + 0.00 * tetmp2
          ne(3) = 0.50 * tetmp1 + 0.50 * tetmp2
          WRITE(PINOUT,*) 'COMPROMISE:',ir,psitarg(ir,2),ne(3),te(3)
        ENDIF

        IF (te(3).LE.-1.0) CALL LoadProbeDataS28(ir,s28probe,2,te(3))
c *HACKISH*
        IF (ne(3).LE.-1.0.AND.ne(3).GT.-100.0) 
     .    CALL LoadProbeDataS28(ir,s28probe,1,ne(3))

        IF (ne(3).LT.-100.0) THEN
c...      Pressure specified, not density.  Calculate ne(3) from the pressure
c         and the specified te(3) value, assuming Ti=Te:
          ne(3) = -ne(3) / (2.0 * te(3))
        ENDIF

c...    Specify te(1) from target data.  The -ve is to trigger
c       te(0)=te(1):
        IF (te(1).EQ.-99.0) te(1) = -kteds(idds(ir,2))

c...    Specify te(5) from target data.  The -ve is to trigger
c       te(6)=te(5):
        IF (te(5).EQ.-99.0) te(5) = -kteds(idds(ir,1))

      ELSE
      ENDIF


      RETURN
99    CONTINUE
      WRITE(0,*) 'IK,IR=',ik,ir,i1,osms28(i1,1)
      WRITE(0,*) 'TAB,TCD=',tab,tcd
      STOP
      END
c
c
c
c
c
c
c
      SUBROUTINE InitializeRingSOL28(ir,count)
      USE mod_sol28_old
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_slcom
      use mod_slcom_sol28
      use mod_sol22_support
      use mod_sl_oldplasma
      IMPLICIT none

      INTEGER ir,count

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'

!      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
!     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
!      REAL
!     .     oldknbs  (MAXNKS,MAXNRS),oldktebs (MAXNKS,MAXNRS),
!     .     oldktibs (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
!     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
!     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)

      COMMON /TIMESTATS/ timeint
      REAL               timeint(100)


      REAL Clock2,GetCs

      INTEGER ik,i1,ir1,pincount
      REAL    lastpinrec(MAXNKS)
      SAVE
c
c     Initialization:
c

c...  Identify cells that correspond to s-parameters:
      IF (new) THEN
        ik0 = 1
        ik1 = 1
        ik2 = 1
        ik3 = 0
        ik4 = nks(ir)
        ik5 = nks(ir)
        ik6 = nks(ir)
        s(0) = 0.0
        IF (s(4).EQ.0.0) s(4) = ksmaxs(ir)
        IF (s(5).EQ.0.0) s(5) = ksmaxs(ir)
        s(6) = ksmaxs(ir)
      ELSE
        ik4 = nks(ir)
        ik5 = nks(ir)
        s(5) = ksmaxs(ir)
      ENDIF

      s12 = 0.5 * (s(1) + s(2))
      s25 = 0.5 * (s(2) + s(3))
      s34 = 0.5 * (s(3) + s(4))
      s45 = 0.5 * (s(4) + s(5))

      DO ik = 1, nks(ir)
        IF (new) THEN
          IF (ksb(ik-1,ir).LT.s(1).AND.ksb(ik,ir).GT.s(1)) ik1 = ik
          IF (ksb(ik-1,ir).LT.s(2).AND.ksb(ik,ir).GT.s(2)) ik2 = ik
          IF (ksb(ik-1,ir).LT.s(3).AND.ksb(ik,ir).GT.s(3)) ik3 = ik
          IF (ksb(ik-1,ir).LT.s(4).AND.ksb(ik,ir).GT.s(4)) ik4 = ik
          IF (ksb(ik-1,ir).LT.s(5).AND.ksb(ik,ir).GT.s(5)) ik5 = ik

          IF (ksb(ik-1,ir).LT.s12 .AND.ksb(ik,ir).GT.s12 ) ik12= ik
          IF (ksb(ik-1,ir).LT.s45 .AND.ksb(ik,ir).GT.s45 ) ik45= ik
        ELSE
        ENDIF
      ENDDO

      IF (.NOT.loadsources) THEN

        CALL RZero(pinmp (1,ir),MAXNKS)
        CALL RZero(osmpmk(1,ir),MAXNKS)

        IF (count.GE.1) THEN

          pincount = pincount + 1

          IF (imaginaryL2) THEN
            v(1) = lastv1 
            DO ik = 1, nks(ir)
              knbs(ik,ir) = lastknbs(ik)
            ENDDO
          ENDIF

c...      Should load OSMREC here and leave PINREC to whatever EIRENE returns, for comparison later... 

c          CALL RZero(pinrec(1,ir),MAXNKS)
c          CALL CalcInitRecom(IKLO,ir)
c          CALL CalcInitRecom(IKHI,ir)

          IF (.FALSE..AND.pincount.GT.1.AND.
     .        ((ringtype(ir).EQ.SOL1.AND.s28superdet   .EQ.2).OR.
     .         (ringtype(ir).EQ.PFZ .AND.s28superdetpfz.EQ.2).OR.
     .         (ringtype(ir).EQ.SOL1.AND.s28superdet   .EQ.4).OR.
     .         (ringtype(ir).EQ.PFZ .AND.s28superdetpfz.EQ.4))) THEN
c...        Need to relax recombination source:
            lastpinrec = 0.0
            DO ik = 1, nks(ir)
              lastpinrec(ik) = pinrec(ik,ir)
            ENDDO

c...        Should load OSMREC here and leave PINREC to whatever EIRENE returns, for comparison later... 

            CALL RZero(pinrec(1,ir),MAXNKS)
            CALL CalcInitRecom(IKLO,ir)
            CALL CalcInitRecom(IKHI,ir)

            DO ik = 1, nks(ir)
              pinrec(ik,ir) = 0.1 * pinrec(ik,ir) + 0.9 * lastpinrec(ik)
            ENDDO

            WRITE(0,*) 'RELAXING PINREC'

          ELSEIF ((ringtype(ir).EQ.SOL1.AND.s28rec   .EQ.2).OR.
     .            (ringtype(ir).EQ.PFZ .AND.s28recpfz.EQ.2)) THEN
            WRITE(0,*) 
            WRITE(0,*) 'NOT RECALCULATING PINREC'
            WRITE(0,*)
          ELSE
            CALL RZero(pinrec(1,ir),MAXNKS)
            CALL CalcInitRecom(IKLO,ir)
            CALL CalcInitRecom(IKHI,ir)
c            WRITE(0,*) 'RELAXING PINREC SOLID'
c            WRITE(0,*) 'PINREC: CALCULATING!'
          ENDIF
          
        ELSEIF (count.EQ.0) THEN
c...      Initialization:
          machno = 1.0
          lastv1 = 0.0
          v(0) = 0.0
          v(1) = 0.0
          teval = 1.6
          pincount = 0

          CALL RZero(kvhs(1,ir),MAXNKS)  ! Needed to add when SOLOPT=6,8 sets background velocity

          IF     (ringtype(ir).EQ.SOL1.AND.s28cfpdrft   .EQ.-2.OR.
     .            ringtype(ir).EQ.PFZ .AND.s28cfppfzdrft.EQ.-2) THEN
            IF (.NOT.ALLOCATED(tmpflx)) ALLOCATE(tmpflx(MAXNKS,MAXNRS))
            CALL CalcRadialDrift(ir)
            CALL CalcRadialDrift(-2)
          ELSEIF (ir.EQ.irsep.AND.
     .            ringtype(ir).EQ.SOL1.AND.s28cfpdrft   .EQ.-3.OR.
     .            ringtype(ir).EQ.PFZ .AND.s28cfppfzdrft.EQ.-3) THEN
            DO ir1 = irsep, irwall-1
              DO ik = 1, nks(ir1)
                knbs (ik,ir1) = oldknbs (ik,ir1)
                ktebs(ik,ir1) = oldktebs(ik,ir1)
                ktibs(ik,ir1) = oldktibs(ik,ir1)
                kvhs (ik,ir1) = oldkvhs (ik,ir1)
              ENDDO
            ENDDO
            IF (.NOT.ALLOCATED(tmpflx)) ALLOCATE(tmpflx(MAXNKS,MAXNRS))
            CALL CalcRadialDrift(-1)
            CALL CalcRadialDrift(-2)
          ENDIF

          machcountL2 = 0
          machadjustL2 = -0.1

          IF (.NOT.new) scalerec = .FALSE.
          holdv2 = kvds(idds(ir,1))
          holdn2 = knds(idds(ir,1))
c          te(0) = kteds(idds(ir,2))

          timeint(6) = Clock2()

          IF (s28fit.EQ.1) fitline = .TRUE.
          calcmom1 = .TRUE.

          CALL RZero(addmp,MAXNKS)

          IF (new) THEN
            p6 = HI
            id0 = idds(ir,2)
            id6 = idds(ir,1)
            isat(0)=knds(id0)*ECH*GetCs(kteds(id0),ktids(id0))
            isat(6)=knds(id6)*ECH*GetCs(kteds(id6),ktids(id6))
          ELSE
            p5 = HI
          ENDIF

        ENDIF

      ELSE
c...    Loading sources from EIRENE.  Bit of a strange BUSINESS ALWAYS REPLACING PINREC!

        IF (count.EQ.0) THEN
c...      Initialization:
          machno = 1.0

          DO ik = 1, nks(ir)
            knbs (ik,ir) = oldknbs (ik,ir)
            ktebs(ik,ir) = oldktebs(ik,ir)
            ktibs(ik,ir) = oldktibs(ik,ir)
            kvhs (ik,ir) = oldkvhs (ik,ir)
          ENDDO


          IF     (ringtype(ir).EQ.SOL1.AND.s28cfpdrft   .EQ.-2.OR.
     .            ringtype(ir).EQ.PFZ .AND.s28cfppfzdrft.EQ.-2) THEN
            IF (.NOT.ALLOCATED(tmpflx)) ALLOCATE(tmpflx(MAXNKS,MAXNRS))
            CALL CalcRadialDrift(ir)
            CALL CalcRadialDrift(-2)
          ELSEIF (ir.EQ.irsep.AND. ! * DANGEROUS *
     .            ringtype(ir).EQ.SOL1.AND.s28cfpdrft   .EQ.-3.OR.
     .            ringtype(ir).EQ.PFZ .AND.s28cfppfzdrft.EQ.-3) THEN
            DO ir1 = irsep, irwall-1
              DO ik = 1, nks(ir1)
                knbs (ik,ir1) = oldknbs (ik,ir1)
                ktebs(ik,ir1) = oldktebs(ik,ir1)
                ktibs(ik,ir1) = oldktibs(ik,ir1)
                kvhs (ik,ir1) = oldkvhs (ik,ir1)
              ENDDO
            ENDDO
            IF (.NOT.ALLOCATED(tmpflx)) ALLOCATE(tmpflx(MAXNKS,MAXNRS))
            CALL CalcRadialDrift(-1)
            CALL CalcRadialDrift(-2)
          ENDIF


          lastv1 = kvhs(ik1,ir)
          v(0) = 0.0
          v(1) = kvhs(ik1,ir) 
          teval = ktebs(ik1,ir) + 0.3

          machcountL2 = 0
          machadjustL2 = -0.1

          IF (.NOT.new) scalerec = .FALSE.

          holdv2 = kvds(idds(ir,1))
          holdn2 = knds(idds(ir,1))

c          te(0) = kteds(idds(ir,2))

          timeint(6) = Clock2()

          IF (s28fit.EQ.1) fitline = .TRUE.
          calcmom1 = .TRUE.

          CALL RZero(addmp,MAXNKS)

          IF (new) THEN
            p6 = HI
            id0 = idds(ir,2)
            id6 = idds(ir,1)
            isat(0)=knds(id0)*ECH*GetCs(kteds(id0),ktids(id0))
            isat(6)=knds(id6)*ECH*GetCs(kteds(id6),ktids(id6))
          ELSE
            p5 = HI
          ENDIF

c...      Set PINION from existing OSMION values, which were likely generated
c         from a previous case with S28MOM.LT.1:


c REDUNDANT
          IF (loadsources) THEN
c          IF (.NOT.loadsources) THEN

            IF     (ringtype(ir).EQ.SOL1) THEN


              IF (s28ion.GT.0.AND..NOT.s28ionset) THEN
                DO ik = 1, nks(ir)
                  pinion(ik,ir) = osmion(ik,ir)
                  DO i1 = 1, 3               
                    pindata(ik,ir,H_ION1+i1-1) = osmion(ik,ir) / 3.0
                  ENDDO
c                  WRITE(0,*) 'DATA:',ik,ir,pinion(ik,ir)
                ENDDO
              ENDIF
c              WRITE(0,*) 'S28ION:',s28ion,s28ionset
c              STOP 'sdfsdfsd'
              IF (s28mom.GT.0.AND..NOT.s28momset) THEN
                DO ik = 1, nks(ir)
                  pinmp(ik,ir) = osmmp(ik,ir)
                ENDDO
              ENDIF

            ELSEIF (ringtype(ir).EQ.PFZ) THEN
              IF (s28ionpfz.GT.0.AND..NOT.s28ionsetpfz) THEN
                DO ik = 1, nks(ir)
                  pinion(ik,ir) = osmion(ik,ir)
                  DO i1 = 1, 3               
                    pindata(ik,ir,H_ION1+i1-1) = osmion(ik,ir) / 3.0
                  ENDDO
                ENDDO
              ENDIF  
              IF (s28mompfz.GT.0.AND..NOT.s28momsetpfz) THEN
                DO ik = 1, nks(ir)
                  pinmp(ik,ir) = osmmp(ik,ir)
                ENDDO
              ENDIF  

            ENDIF          

          ENDIF
        ENDIF

        IF (imaginaryL2) THEN
          v(1) = lastv1 
          DO ik = 1, nks(ir)
            knbs(ik,ir) = lastknbs(ik)
          ENDDO
        ENDIF

* DO THIS BECAUSE DON'T WANT INTERFERENCE FROM MULTIPLIERS?  NOT OVERLY REALISTIC OR CONSISTENT  


c... NEED A CHECK OF S28REC OPTION **********

c        IF (.FALSE.) THEN


        IF (ringtype(ir).EQ.SOL1.AND.s28superdet   .EQ.2.OR.
     .      ringtype(ir).EQ.PFZ .AND.s28superdetpfz.EQ.2) THEN
c...      Need to relax recombination source:
          lastpinrec = 0.0
          DO ik = 1, nks(ir)
            lastpinrec(ik) = pinrec(ik,ir)
          ENDDO

c...      Should load OSMREC here and leave PINREC to whatever EIRENE returns, for comparison later... 

          CALL RZero(pinrec(1,ir),MAXNKS)
          CALL CalcInitRecom(IKLO,ir)
          CALL CalcInitRecom(IKHI,ir)

          DO ik = 1, nks(ir)
            pinrec(ik,ir) = 0.1 * pinrec(ik,ir) + 0.9 * lastpinrec(ik)
          ENDDO

        ELSEIF ((ringtype(ir).EQ.SOL1.AND.s28rec   .EQ.2).OR.
     .          (ringtype(ir).EQ.PFZ .AND.s28recpfz.EQ.2)) THEN
          WRITE(0,*) 
          WRITE(0,*) 'NOT RECALCULATING PINREC'
          WRITE(0,*)

        ELSE
          CALL RZero(pinrec(1,ir),MAXNKS)
          CALL CalcInitRecom(IKLO,ir)
          CALL CalcInitRecom(IKHI,ir)
c          WRITE(0,*) 'PINREC: CALCULATING!'
        ENDIF

      ENDIF



      RETURN
 99   STOP
      END
c
c
c
c
      SUBROUTINE ParticleConservation(ir,count,machadjust1,machadjust2)

      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_cedge2d
      use mod_slcom
      use mod_slcom_sol28
      IMPLICIT none

      INTEGER ir,count,s28ion_,s28cfp_,s28cfpnt_,s28superdet_,s28rec_,
     .        s28cfpdrft_
      REAL    machadjust1,machadjust2


c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'cedge2d'
c... NEED TO GET THIS OUT OF HERE OR EVENTUALLY TROUBLE WITH _'s
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'

      REAL GetEAD,GetCs

      INTEGER ik,ike,iks
      REAL    sval,frac,rdum1,rdum2,rdum3,ion03,ion36,cfp03,cfp36,flxx,
     .        par03,par36,tmpcfp(MAXNKS,MAXNRS),ionbound,tt,nt,sv,av,
     .        mfp,rec56
      REAL*8  integral



c...  Process sources:
      CALL RZero(osmion(1,ir),MAXNKS)
      CALL RZero(osmrec(1,ir),MAXNKS)
c...  Cross-field particle source near target specified, not from continuity:
      CALL RZero(osmcfp(1,ir),MAXNKS)


c      DO ik = 1, nks(ir)
c        osmrec(ik,ir) = pinrec(ik,ir) * mulrecl(ir)
c      ENDDO
c *MOVE*
      s28rec_ = 99
      IF (ringtype(ir).EQ.SOL1) s28rec_ = s28rec
      IF (ringtype(ir).EQ.PFZ ) s28rec_ = s28recpfz
c
c     S28REC -1 -
c             0 - 
c             1 - load PINREC which is recomputed locally with each iteration
c             2 - load PINREC using the value returned from EIRENE
c

c *MOVE*
      s28superdet_ = 99
      IF (ringtype(ir).EQ.SOL1) s28superdet_ = s28superdet
      IF (ringtype(ir).EQ.PFZ ) s28superdet_ = s28superdetpfz
c
c
c

      IF     (s28rec_.EQ.20) THEN
c...    Assing volume recombination source from 2D fluid code:
        osmrec(1:nks(ir),ir) = e2drec(1:nks(ir),ir)

      ELSEIF (loadsources.AND.s28rec_.GT.0) THEN

        DO ik = ik0, ik6
          osmrec(ik,ir) = pinrec(ik,ir)
        ENDDO

      ELSEIF (s28rec_.EQ.1) THEN

        IF (s28superdet_.EQ.1) THEN
          DO ik = ik0, ik3
            osmrec(ik,ir) = pinrec(ik,ir) * mulrecl(ir) * machadjust1
          ENDDO
          DO ik = ik3+1, ik6
            osmrec(ik,ir) = pinrec(ik,ir) * mulrecl(ir) * machadjust2
          ENDDO
        ELSE
          DO ik = ik0, ik3
            osmrec(ik,ir) = pinrec(ik,ir) * mulrecl(ir) 
          ENDDO
          DO ik = ik3+1, ik6
            osmrec(ik,ir) = pinrec(ik,ir) * mulrecl(ir) 
          ENDDO
        ENDIF

      ELSEIF (s28rec_.EQ.0) THEN

        DO ik = 1, nks(ir)
          osmrec(ik,ir) = 0.0
        ENDDO

      ELSE
        CALL ER('ParticleConservation','Unrecognized S28REC option',
     .          *99)
      ENDIF




c *MOVE*
      s28cfpdrft_ = 99
      IF (ringtype(ir).EQ.SOL1) s28cfpdrft_ = s28cfpdrft
      IF (ringtype(ir).EQ.PFZ ) s28cfpdrft_ = s28cfppfzdrft
c
c     S28CFPDRFT  0 - none
c                 1 - ExB
c                 2 - 
c

      IF (loadsources.AND.s28cfpdrft_.NE.0) THEN

        IF (ir.EQ.irsep)  
     .    WRITE(0,*) 'FLUX B:',osmcfpflx(nks(irsep),irsep,2)



        IF (ir.EQ.irsep)  
     .    WRITE(0,*) 'FLUX C:',osmcfpflx(nks(irsep),irsep,2)

c ... NOT SURE THAT THE FINAL OSMCFP SHOULD BE ZERO FOR NON ZERO ALREADY ASSIGNED 
c     OSMCFP CELLS 

        IF     (s28cfpdrft_.EQ.1 .OR.s28cfpdrft_.EQ.-1.OR.
     .                               s28cfpdrft_.EQ.-2.OR.
     .                               s28cfpdrft_.EQ.-3) THEN

          DO ik = ik3, ik6
            osmcfp(ik,ir) = osmcfp(ik,ir) + osmcfpflx(ik,ir,2)
          ENDDO

        ELSEIF (s28cfpdrft_.EQ.2) THEN
 
          DO ik = ik3, ik6
            osmcfp(ik,ir) = osmcfp(ik,ir) + osmcfpflx(ik,ir,1) 
c            WRITE(0,*) 'OSM:',ik,ir,osmcfpflx(ik,ir,1)
          ENDDO

        ELSEIF (s28cfpdrft_.EQ.5) THEN

          DO ik = ik0, ik6
            osmcfp(ik,ir) = osmcfp(ik,ir) + osmcfpflx(ik,ir,2)
          ENDDO

        ELSEIF (s28cfpdrft_.EQ.6) THEN
 
          DO ik = ik0, ik6
            osmcfp(ik,ir) = osmcfp(ik,ir) + osmcfpflx(ik,ir,1) 
          ENDDO

        ELSE
          CALL ER('ParticleConservation','Unrecognized S28CFPDRFT '//
     .            'option',*99)
        ENDIF

      ENDIF


c *MOVE*
      s28cfpnt_ = 99
      IF (ringtype(ir).EQ.SOL1) s28cfpnt_ = s28cfpnt
      IF (ringtype(ir).EQ.PFZ ) s28cfpnt_ = s28cfppfznt
c
c     S28CFPNT  0 - uniform along the entire ring
c               1 - parabolic (?) along the entire ring
c               2 - 
c
      s28superdet_ = 99
      IF (ringtype(ir).EQ.SOL1) s28superdet_ = s28superdet
      IF (ringtype(ir).EQ.PFZ ) s28superdet_ = s28superdetpfz

      tmpcfp = 0.0

      IF (count.GT.0.AND.s28cfpnt_.GT.0) THEN
c      IF (s28superdet_.EQ.3.AND.count.GT.0) THEN
c...  
        IF (detach1) THEN

          IF     (s28cfpnt_.EQ.2.OR.s28cfpnt_.EQ.5) THEN
            DO ik = ik0, (ik1+ik2)/2
              tmpcfp(ik,ir) = 1.0
            ENDDO
          ELSEIF (s28cfpnt_.EQ.3) THEN
            DO ik = ik0, ik2 
              tmpcfp(ik,ir) = MAX(1.0,pinrec(ik,ir))
            ENDDO
          ELSEIF (s28cfpnt_.EQ.4) THEN
            DO ik = ik12, ik2
c            DO ik = (ik1+ik2)/2, ik2
              tmpcfp(ik,ir) = 1.0
            ENDDO
          ELSEIF (s28cfpnt_.EQ.6) THEN
            DO ik = ik0, ik2 
              tmpcfp(ik,ir) = MAX(1.0,knbs(ik,ir))
            ENDDO
          ELSEIF (s28cfpnt_.EQ.7) THEN
            DO ik = ik0,ik1
              tmpcfp(ik,ir) = 1.0
            ENDDO
          ELSE
            CALL ER('ParticleConservation','Invalid S28CFPNT value',*99)
          ENDIF

c...WORK IN PROGRESS
          IF (s28superdet_.EQ.4) THEN

            IF (iflexopt(5).EQ.20.AND.ir.LT.20) THEN

              frac = 1.0 / machadjust1

c              IF (ir.GE.14.AND.ir.LE.17) THEN
c                flxx = MIN(-4.0E+5/ECH,flx0)            ! Bit of a lame fix, but without this the target flux
c              ELSE
                flxx = MIN(-1.0+5/ECH,flx0)            ! Bit of a lame fix, but without this the target flux
c              ENDIF

              IF (.FALSE..AND.machadjust1.EQ.0.9) THEN
c...            RESET!
                CALL CalcIntegral2(osmrec(1,ir),ik0,ik2,ir,integral,4)  ! * HACK *
                rec06 = SNGL(integral)
                CALL CalcIntegral2(osmcfp(1,ir),ik0,ik2,ir,integral,4)
                cfp06 = SNGL(integral)
                flxx = (flx0 - rec06) * 0.5
c                flxx = flx0 - rec06 + cfp06
                WRITE(0,*) '......    :',machadjust1,flx0,-rec06,cfp06
              ENDIF

              WRITE(0,*) 'CF_SCALING:',frac,flx0,flxx ! may not be large enought to control the flow... 
              WRITE(0,*) '......    :',frac*flxx,machadjust1
            ELSE
              frac = 2.0 * (1.0 - machadjust1)
              flxx = MIN(-1.0+5/ECH,flx0)       
            ENDIF
            CALL Normalize(ik0,ik2,ir,tmpcfp,DBLE(frac*flxx),1,'CF9')
          ELSEIF (ir.LT.irwall) THEN
            CALL Normalize(ik0,ik2,ir,tmpcfp,DBLE(0.70*flx0),1,'CF9')
          ELSE
            CALL Normalize(ik0,ik2,ir,tmpcfp,DBLE(1.00*flx0),1,'CF9')
          ENDIF

        ENDIF


        IF (detach2) THEN

          IF     (s28cfpnt_.EQ.2) THEN
            DO ik = ik4+5,ik6
              tmpcfp(ik,ir) = 1.0
            ENDDO
          ELSEIF (s28cfpnt_.EQ.3) THEN
            DO ik = ik4+5,ik6
              tmpcfp(ik,ir) = MAX(1.0,pinrec(ik,ir))
            ENDDO
          ELSEIF (s28cfpnt_.EQ.4) THEN
            DO ik = ik4, ik45
c            DO ik = ik4, (ik4+ik5)/2
              tmpcfp(ik,ir) = 1.0
            ENDDO
          ELSEIF (s28cfpnt_.EQ.6) THEN
            DO ik = ik4, ik6
              tmpcfp(ik,ir) = MAX(1.0,knbs(ik,ir))
            ENDDO
          ELSEIF (s28cfpnt_.EQ.7) THEN
            DO ik = ik5,ik6
              tmpcfp(ik,ir) = 1.0
            ENDDO
          ELSE
            CALL ER('ParticleConservation','Invalid S28CFPNT value',*99)
          ENDIF

c...WORK IN PROGRESS
          IF (s28superdet_.EQ.4) THEN
            frac = 2.0 * (1.0 - machadjust2)
            WRITE(0,*) 'CF_SCALING:',frac
            CALL Normalize(ik4,ik6,ir,tmpcfp,DBLE(-frac*flx6),1,'CF9')

          ELSEIF (ir.LT.irwall) THEN
c            CALL Normalize(ik4,ik6,ir,tmpcfp,DBLE(-0.70*flx6),1,'CF9')
            CALL CalcIntegral3(osmrec,ik5+1,ik6,ir,rec56,4)
            CALL Normalize(ik5+1,ik6,ir,tmpcfp,DBLE(rec56),1,'RC1')
          ELSE
            CALL Normalize(ik4,ik6,ir,tmpcfp,DBLE(-1.00*flx6),1,'CF9')
          ENDIF

        ENDIF

        DO ik = 1, nks(ir)
          osmcfp(ik,ir) = osmcfp(ik,ir) + tmpcfp(ik,ir)
        ENDDO

      ENDIF



      IF (new) THEN
        flx0 = -isat(0) / ECH
        flx6 =  isat(6) / ECH
c INTEGRALS SET UP RIGHT?
        CALL CalcIntegral3(osmrec,ik0  ,ik3,ir,rec03,4)
        CALL CalcIntegral3(osmrec,ik3+1,ik6,ir,rec36,4)

c *MOVE*
        s28ion_ = 99
        IF (ringtype(ir).EQ.SOL1) s28ion_ = s28ion
        IF (ringtype(ir).EQ.PFZ ) s28ion_ = s28ionpfz
c
c       S28ION 
c              -4 - all ionisation between ik1-2 and ik4-5
c              -3 - all ionisation near the target
c              -2 -
c              -1 -
c               0 - standard
c               1 - as 0 for 1st iteration, then from PIN
c               2 - as 0 for 1st iteration, then from PIN with source bounded by S28IONBOUND
c               3 - as 2 but only on the inside
c               4 - as 2 but only on the outside
c
        IF     (s28ion_.EQ.20) THEN
c...      Assing ionisation source from 2D fluid code:
          osmion(1:nks(ir),ir) = e2dion(1:nks(ir),ir) 

        ELSEIF (loadsources.AND.s28ion_.GT.0) THEN

          IF     (s28ion_.EQ.1) THEN 
c...        
            DO ik = 1, nks(ir)
              osmion(ik,ir) = pinion(ik,ir)
            ENDDO

          ELSEIF (s28ion_.EQ.2.OR.s28ion_.EQ.3.OR.s28ion_.EQ.4) THEN
c...        

            DO ik = 1, nks(ir)
c              osmion(ik,ir) = 0.5*pinion(ik,ir)
              osmion(ik,ir) = pinion(ik,ir)

c              IF (ik.LT.10) WRITE(0,*) 'OSMION1:',ik,ir,osmion(ik,ir),
c     .          pinion(ik,ir)
            ENDDO

            CALL CalcIntegral3(osmcfp,ik0  ,ik3,ir,cfp03,4)
            CALL CalcIntegral3(osmcfp,ik3+1,ik6,ir,cfp36,4)

            CALL CalcIntegral3(osmion,ik0  ,ik3,ir,ion03,4)
            CALL CalcIntegral3(osmion,ik3+1,ik6,ir,ion36,4)

            WRITE(PINOUT,*) 'IONCHECK A:',ir,ion03,ion36

            IF (s28ion_.EQ.4) THEN
              ionbound = 1.0
            ELSE
              ionbound = s28ionbound
              IF (iflexopt(6).EQ.30) THEN
                WRITE(0,*)
                WRITE(0,*) '*** BOGUS S28IONBOUND ***'
                WRITE(0,*)
                IF (ir.GE.20) s28ionbound = 5.0
                IF (ir.GE.18.AND.ir.LE.19) s28ionbound = 2.0
                IF (ir.GE.16.AND.ir.LE.17) s28ionbound = 1.5
              ENDIF

            ENDIF

            par03 = -flx0 + rec03 - cfp03
c...        Low IK index:
            IF     (ion03.GT.par03*ionbound) THEN
              s28ionfrac(IKLO,ir) = (par03 * ionbound) / ion03
              ion03 = par03 * ionbound
              WRITE(PINOUT,*) 'IONCHECK: INNER TOO HIGH'
            ELSEIF (ion03.LT.par03/ionbound) THEN
              s28ionfrac(IKLO,ir) = (par03 / ionbound) / ion03
              ion03 = par03 / ionbound
              WRITE(PINOUT,*) 'IONCHECK: INNER TOO LOW'
            ENDIF

            IF (s28ion_.EQ.3) THEN
              ionbound = 1.0
            ELSE
              ionbound = s28ionbound
            ENDIF

c...        High IK index:
            par36 = flx6 + rec36 - cfp36
            IF     (ion36.GT.par36*ionbound) THEN
              s28ionfrac(IKHI,ir) = (par36 * ionbound)/ion36
              ion36 = par36 * ionbound
              WRITE(PINOUT,*) 'IONCHECK: OUTER TOO HIGH'
            ELSEIF (ion36.LT.par36/ionbound) THEN
              s28ionfrac(IKHI,ir) = (par36 / ionbound) / ion36
              ion36 = par36 / ionbound
              WRITE(PINOUT,*) 'IONCHECK: OUTER TOO LOW'
            ENDIF
           

            CALL Normalize(ik0  ,ik3,ir,osmion,DBLE(ion03),0,'IO7')
            CALL Normalize(ik3+1,ik6,ir,osmion,DBLE(ion36),0,'IO8')

            WRITE(PINOUT,*) 'IONCHECK B:',ir,ion03,ion36

            CALL CalcIntegral3(osmion,ik0  ,ik3,ir,ion03,4)
            CALL CalcIntegral3(osmion,ik3+1,ik6,ir,ion36,4)
 
            WRITE(PINOUT,*) '        :',ir,ion03,ion36
            WRITE(PINOUT,*) '        :',-flx0,rec03,-flx0+rec03
            WRITE(PINOUT,*) '        :', flx6,rec36, flx6+rec36

c            IF (ir.EQ.16) STOP 'sdfsdfsd'

          ELSE
            STOP 'SOL28: S28ION OPTION NOT SPECIFIED'
          ENDIF

        ELSE

          IF (detach1) THEN
            IF (count.GT.0.AND.s28ion_.EQ.-1) THEN
c...          Find region where Te < 5 eV and apply uniform momentum loss:
              iks = 1
              DO iks = 1, ik3
                IF (ktebs(iks+1,ir).GT.5.0) EXIT 
              ENDDO
              IF (iks.EQ.ik3+1) STOP 'PROBLEM C'
              ike = ikto
              IF (s28ion_.EQ.-1) THEN
c...            Triangular source:
                DO ik = iks, ike
                  osmion(ik,ir) = (kss(ik ,ir) - kss(ike,ir)) / 
     .                            (kss(iks,ir) - kss(ike,ir))
                ENDDO             
              ENDIF

            ELSEIF (s28ion_.EQ.-3) THEN

              IF (stopopt2.EQ.31) THEN
                DO ik = ik0, ik1        ! *** HACK ***
                   osmion(ik,ir) = 1.0 
                ENDDO
              ELSE
                DO ik = ik0, ik2  
                   osmion(ik,ir) = 1.0 
                ENDDO
              ENDIF
            ELSEIF (s28ion_.EQ.-4) THEN

              DO ik = ik1, ik2
                 osmion(ik,ir) = 1.0 
              ENDDO

            ELSE
c              DO ik = ik2-1, ik2+1
              DO ik = (ik1+ik2)/2-1,ik2+4
                 osmion(ik,ir) = 1.0 
              ENDDO
            ENDIF
            CALL Normalize(ik0,ik3,ir,osmion,DBLE(-flx0+rec03),0,'IO6')

          ELSE
c...        Ionisation profile:

c..MAKE INDEPENDENT OF GRID RESOLUTION, MAKE TRIANGE, ADD OPTIONS
c            STOP 'FIX IONISTAION PROFILE FOR ATTACHED TARGETS'


            IF (.NOT..TRUE.) THEN
c...          Ionisation rate coefficient: 
              tt = te(0)
              nt = -(flx0 * ECH) / ECH / GetCs(tt,tt)
              sv = GetEAD(tt,nt,1,'H.4 ') * 1.0E-06
              av = SQRT(2.0 * tt * ECH / crmb / AMU)
              mfp = av / (sv * nt) 
              DO ik = ik0, ik3
                osmion(ik,ir) = EXP( -kps(ik,ir) / mfp)
              ENDDO
            ELSE
              ike = ik0+5
c              ike = ik0+14
              DO ik = ik0, ike
                 osmion(ik,ir) = 1.0 - (kss(ik ,ir) - kss(ik0,ir)) /
     .                                 (kss(ike,ir) - kss(ik0,ir))
              ENDDO
            ENDIF

c...        Normalize the ionisation source to give particle balance on 
c           half-ring (as defined by IK3):
c            CALL CalcIntegral2(osmion(1,ir),ik0,ik3,ir,integral,4)
c            DO ik = ik0, ik3
c              osmion(ik,ir)=osmion(ik,ir)*(-flx0+rec03)/SNGL(integral)
c            ENDDO
            CALL Normalize(ik0,ik3,ir,osmion,DBLE(-flx0+rec03),0,'IO6')

          ENDIF

          IF (detach2) THEN
            IF (count.GT.0.AND.s28ion_.EQ.-1) THEN
c...          Find region where Te < 5 eV and apply uniform momentum loss:
              ike = nks(ir)
              DO ike = nks(ir), ik3+1, -1
                IF (ktebs(ike-1,ir).GT.5.0) EXIT 
              ENDDO
              IF (ike.EQ.ik3) STOP 'PROBLEM D'
              iks = ikti
              IF (s28ion_.EQ.-1) THEN
c...            Triangular source:
                DO ik = iks, ike
                  osmion(ik,ir) = (kss(ik ,ir) - kss(iks,ir)) / 
     .                            (kss(ike,ir) - kss(iks,ir))
                ENDDO             
              ENDIF

            ELSEIF (s28ion_.EQ.-2) THEN

              ike = ik6-10
              DO ik = ik6, ike, -1
                 osmion(ik,ir) = 1.0 - (kss(ik6,ir) - kss(ik ,ir)) /
     .                                 (kss(ik6,ir) - kss(ike,ir))
              ENDDO

            ELSEIF (s28ion_.EQ.-3) THEN

              DO ik = ik5, ik6
                 osmion(ik,ir) = 1.0 
              ENDDO

            ELSEIF (s28ion_.EQ.-4) THEN

              DO ik = ik4,ik5
                 osmion(ik,ir) = 1.0 
              ENDDO

            ELSE
c              DO ik = ik4-5, ik4
c...DIII-D default:
              DO ik = ik4-5, (ik4+ik5)/2+1
                 osmion(ik,ir) = 1.0 
              ENDDO
            ENDIF      
            CALL Normalize(ik3+1,ik6,ir,osmion,DBLE(flx6+rec36),0,'IO6')

          ELSE
c...        Ionisation profile:

c..MAKE INDEPENDENT OF GRID RESOLUTION, MAKE TRIANGE, ADD OPTIONS
c            STOP 'FIX IONISTAION PROFILE FOR ATTACHED TARGETS'


            IF (.NOT..TRUE.) THEN
c...          Ionisation rate coefficient: 
              tt = te(6)
              nt = (flx6 * ECH) / ECH / GetCs(tt,tt)
              sv = GetEAD(tt,nt,1,'H.4 ') * 1.0E-06
              av = SQRT(2.0 * tt * ECH / crmb / AMU)
              mfp = av / (sv * nt) 
              DO ik = ik6, ik3+1, -1
                osmion(ik,ir) = EXP( -(kpmaxs(ir) - kps(ik,ir)) / mfp)
              ENDDO
            ELSE
              ike = ik6-5
c              ike = ik6-20
              DO ik = ik6, ike, -1
                 osmion(ik,ir) = 1.0 - (kss(ik6,ir) - kss(ik ,ir)) /
     .                               (kss(ik6,ir) - kss(ike,ir))
              ENDDO
            ENDIF

c SHOULD IT REALLY BE IK3+1 -- THESE INTEGRALS SETUP RIGHT? TO INCLUDE HALF RINGS AT THE START?
c            CALL CalcIntegral2(osmion(1,ir),ik3+1,ik6,ir,integral,4)
c            DO ik = ik3+1, ik6
c              osmion(ik,ir)=osmion(ik,ir)*(flx6+rec35)/SNGL(integral)
c            ENDDO
            CALL Normalize(ik3+1,ik6,ir,osmion,DBLE(flx6+rec35),0,'IO6')
          ENDIF
c...      Copy ionsation profile into standard PIN array:
c *NOW DONE BELOW*
c          DO ik = 1, nks(ir)
c            pinion(ik,ir) = osmion(ik,ir)
c          ENDDO
        ENDIF



c...    Satisfy particle conservation:


        CALL CalcIntegral2(osmion(1,ir),ik0,ik6,ir,integral,4)
        ion06 = SNGL(integral)
        CALL CalcIntegral2(osmrec(1,ir),ik0,ik6,ir,integral,4)
        rec06 = SNGL(integral)
        CALL CalcIntegral2(osmcfp(1,ir),ik0,ik6,ir,integral,4)
        cfp06 = SNGL(integral)
        cfp06 = -(ion06 + flx0 - flx6 - rec06 + cfp06)


c *MOVE*
        s28cfp_ = 99
        IF (ringtype(ir).EQ.SOL1) s28cfp_ = s28cfp
        IF (ringtype(ir).EQ.PFZ ) s28cfp_ = s28cfppfz
c
c       S28CFP  0 - uniform along the entire ring
c               1 - parabolic (?) along the entire ring
c               2 - uniform but no cf source in the detached region 
c               3 - ballooning ...
c

c...    Default:
        tmpcfp = 0.0
        DO ik = ik0, ik6
       
          IF     (s28cfp_.EQ.0) THEN
c...        Linear for SOL:
            IF (osmcfp(ik,ir).EQ.0.0) tmpcfp(ik,ir)=1.0*SIGN(1.0,cfp06)

          ELSEIF (s28cfp_.EQ.1) THEN
c...        Parabolic for PFZ:
            deltas = 0.5 * (s(6) - s(0))
            sval = (kss(ik,ir) - 0.5 * (s(0) + s(6))) / deltas
            IF (osmcfp(ik,ir).EQ.0.0) 
     .        tmpcfp(ik,ir)=MAX(0.0,-sval**2+1.0)

          ELSEIF (s28cfp_.EQ.2) THEN

            IF (osmcfp(ik,ir).EQ.0.0.AND.
     .          ((ik.LE.ik3.AND.(.NOT.detach1.OR.ik.GT.ik2)).OR.
     .           (ik.GT.ik3.AND.(.NOT.detach2.OR.ik.LT.ik4))))
     .        tmpcfp(ik,ir)=1.0*SIGN(1.0,cfp06)

          ELSEIF (s28cfp_.EQ.3) THEN

            IF (osmcfp(ik,ir).EQ.0.0.AND.
     .          ((ik.LE.ik3.AND.(.NOT.detach1.OR.ik.GT.ik2)).OR.
     .           (ik.GT.ik3.AND.(.NOT.detach2.OR.ik.LT.ik4))))
     .        tmpcfp(ik,ir)=(rs(ik,ir)**2)*SIGN(1.0,cfp06)

          ELSEIF (s28cfp_.EQ.4) THEN
c...        Linear for SOL around core:
            IF (thetag(ik,ir).GE.thetag(ikto+1,irsep).AND.
     .          thetag(ik,ir).LE.thetag(ikti-1,irsep))
     .        tmpcfp(ik,ir)=1.0*SIGN(1.0,cfp06)

          ELSEIF (s28cfp_.EQ.20) THEN
c...        Just the outside:
            IF (thetag(ik,ir).GE.94.6                .AND.
     .          thetag(ik,ir).LE.thetag(ikti-1,irsep))
     .        tmpcfp(ik,ir)=1.0*SIGN(1.0,cfp06)

          ELSEIF (s28cfp_.EQ.21) THEN
c...        Just the outside:
            IF (thetag(ik,ir).GE.thetag(ikto+1,irsep).AND.
     .          thetag(ik,ir).LE.94.6                )
     .        tmpcfp(ik,ir)=1.0*SIGN(1.0,cfp06)

          ELSE
            CALL ER('ParticleConservation','Invalid S28CFP value',*99)
          ENDIF
        ENDDO

        CALL Normalize(ik0,ik6,ir,tmpcfp,DBLE(cfp06),1,'CF6')
c        CALL Normalize(ik0,ik6,ir,osmcfp,DBLE(cfp06),1,'CF6')
        DO ik = 1, nks(ir)
          osmcfp(ik,ir) = osmcfp(ik,ir) + tmpcfp(ik,ir)
        ENDDO


c        CALL Normalize(ik0,ik6,ir,osmcfp,DBLE(cfp06),0,'CF6')

c        WRITE(0,*) 'SOURCES:',-(ion06 + flx0 - flx6 - rec06)
c        WRITE(0,*) 'SOURCES:',ion06,rec06,cfp06,ion06-rec06+cfp06
c        DO ik = 1, nks(ir)
c          WRITE(0,*) 'OSMION',ik,ir,osmion(ik,ir)
c        ENDDO
c        STOP 'sdfsd'

c *TEMP*
c        CALL RZero(osmcfp(1,ir),MAXNKS)


        IF (.FALSE..AND.detach1) THEN
          CALL CalcIntegral3(osmcfp,1,ik1,ir,cfp01,2)
          flx1 = ion01 + flx0 - rec01 + cfp01
        ENDIF

c        WRITE(0,*) 'FLX0,6',flx0,flx6
        CALL CalcIntegral2(osmion(1,ir),ik0,ik6,ir,integral,4)
        rdum1 = SNGL(integral)
c        WRITE(0,*) 'INT0',rdum1
        CALL CalcIntegral2(osmrec(1,ir),ik0,ik6,ir,integral,4)
        rdum2 = SNGL(integral)
c        WRITE(0,*) 'INT6',rdum2
        CALL CalcIntegral2(osmcfp(1,ir),ik0,ik6,ir,integral,4)
        rdum3 = SNGL(integral)
c        WRITE(0,*) 'CFP6',rdum3
c        IF (sloutput) 
c     .    WRITE(0,*) 'PARTICLE:',rdum1-rdum2+rdum3+flx0-flx6

c        DO ik = 1, nks(ir)
c          WRITE(0,*) 'SOURCES:',ik,ir,osmion(ik,ir),osmcfp(ik,ir)
c        ENDDO


c            DO ik = 1, nks(ir)
c              IF (ik.LT.10) WRITE(0,*) 'OSMION2:',ik,ir,osmion(ik,ir),
c     .             pinion(ik,ir)
c            ENDDO

c
c * RETURN *
c



        RETURN
      ENDIF


      RETURN
 99   WRITE(0,*) 'S28CFP_=',s28cfp_
      WRITE(0,*) 'IR,RINGTYPE=',ir,ringtype(ir)
      STOP
      END
c
c
c
c
c
      SUBROUTINE CheckLineMatch(ik,ir,linedat,mode,line,tarden)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_slcom
      use mod_slcom_sol28
      IMPLICIT none

      INTEGER ik,ir,mode,line
      REAL    linedat,tarden


c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'


      REAL GetEAD

      INTEGER FAC32,FAC42,FAC52
      LOGICAL cont
      REAL    pop,emission,diff,adjust,lastsign

      tarden = 0.0
      adjust = 0.1
      lastsign = 0.0

      IF (linedat.LT.1.0E+20) THEN
        WRITE(0,*) 'LINEDAT TOO LOW',ik,ir  
        RETURN
      ENDIF

c...  Einstein coefficients for Balmer lines (from HALPHA.F in EIRENE):
      FAC32=4.410E+07
      FAC42=8.419E+06
      FAC52=2.530E+06

      tarden = knbs(ik,ir)

      cont = .TRUE.
      DO WHILE (cont)

        IF (line.EQ.H_BGAMMA) THEN
          pop = GetEAD(ktebs(ik,ir),tarden,20,'H.12')        
          emission = FAC52 * pop * tarden
        ELSE
          CALL ER('CheckLineMatch','Unrecognized line',*99)
        ENDIF

c...    % difference appropriate?
        diff = (emission - linedat) / linedat 

        WRITE(0,'(A,2I6,1P,3E10.2,0P,F10.4)') 
     .    'LINEDAT:',ik,ir,linedat,emission,tarden,diff

        IF (ABS(diff).GT.0.05) THEN

          IF (lastsign.NE.0.0.AND.SIGN(1.0,diff).NE.lastsign) THEN
            adjust = 0.3 * adjust 
          ENDIF
          lastsign = SIGN(1.0,diff)
          IF (diff.LT.0.0) THEN
            tarden = tarden * (1.0 + adjust)
          ELSE
            tarden = tarden * (1.0 - adjust)
          ENDIF
          lastsign = SIGN(1.0,diff)          

        ELSE
          cont = .FALSE.
        ENDIF

      ENDDO

c      IF (tarden.EQ.knbs(ik,ir)) tarden = 0.0

c      WRITE(0,'(A,2I6,1P,2E10.2,0P,F10.2)') 
c     .  'LINEDAT:',ik,ir,linedat,emission,tarden

      RETURN
 99   STOP
      END
c
c
c
c
c
c
c
c
      SUBROUTINE MomentumConservation(ir,count)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_cedge2d
      use mod_slcom
      use mod_slcom_sol28
      IMPLICIT none

      INTEGER ir,count

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'cedge2d'
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'

      REAL CalcPressure,GetCs

      INTEGER ik,iks,ike,s28mom_,s28cfm_,irlast,s28cfpnt_
      LOGICAL status,warn1
      REAL    tmpmp(MAXNKS,MAXNRS),frac,sval,cs,
     .        a1,b1,
     .        x(0:2),tmpmp15,tmpmp01,A,B,B2,bestA1,bestB1,bestB2,
     .        bestA2,val3,minval3,pr,tmppr,mp(0:2),A2,lasttmppr,
     .        mom01,mom15,mom12,mom45,mom06,alpha,mom03,mom36,
     .        mom35,mom56,mom13,
     .        lastx1,lastx2,lastmp1,lastA1,lastA2,lastB1,lastB2,v1,v5
      DATA irlast,warn1 /-1, .FALSE./
      SAVE


      IF (new) THEN
        lastp0 = p0
        lastp6 = p6

c...    *NOT* using velocity at peak Te location:
        p3 = CalcPressure(n(3),te(3),ti(3),0.0)*ECH
c        WRITE(0,*) 'n(:',n(3),te(3),ti(3)
c        WRITE(0,*) 'p3:',p3,n(3),te(3),ti(3)
c        STOP 'sdfsd'

        CALL RZero(osmmp(1,ir),MAXNKS)

c...
        DO ik = 1, nks(ir)
          osmmp(ik,ir) = osmmp(ik,ir) + addmp(ik)
        ENDDO  


c *MOVE*
c... Eventually pass S28MOM_, PINMP, etc. so that there are not potential mixups...
        IF (ringtype(ir).EQ.SOL1) s28mom_ = s28mom
        IF (ringtype(ir).EQ.PFZ ) s28mom_ = s28mompfz
c
c       S28MOM 
c              -3 - uniform to S28MOMTE with triangular component near the target
c              -2 - exponential decay from target
c              -1 - attached-none, detached-triangular between target and Te=S28MOMTE eV point
c               0 - attached-none, detached-uniform between target and Te = S28MOMTE eV point
c               1 - as  0 for 1st iteration, then from PIN
c               2 - as -1 for 1st iteration, then from PIN
c               3 - as -2 for 1st iteration, then from PIN
c               4 - as -3 for 1st iteration, then from PIN
c               5 - as -3 for 1st iteration, then from PIN with uniform source for target-upstream difference
c              
c
c...    PIN contribution:
        IF (s28mom_.EQ.20) THEN
          osmmp(1:nks(ir),ir) = e2dpepara(1:nks(ir),ir) 

        ELSEIF (loadsources.AND.s28mom_.GT.0) THEN
c...        
          IF     (s28mom_.EQ.1.OR.s28mom_.EQ.2.OR.
     .            s28mom_.EQ.5) THEN
c ***OUTSIDE ONLY***
            DO ik = ik3+1, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + pinmp(ik,ir)
            ENDDO

          ELSEIF (s28mom_.EQ.6) THEN
c ***INSIDE ONLY***
            DO ik = ik0, ik3
              osmmp(ik,ir) = osmmp(ik,ir) + pinmp(ik,ir)
            ENDDO
          ELSEIF (s28mom_.EQ.7) THEN
c...        Everywhere:
            DO ik = ik0, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + pinmp(ik,ir)
c              osmmp(ik,ir) = osmmp(ik,ir) + 2.0 * pinmp(ik,ir)
            ENDDO
          ELSEIF (s28mom_.EQ.8) THEN
c...        Everywhere and vol. rec. momentum loss included:
            DO ik = ik0, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + pinmp(ik,ir) -
     .                       crmb*AMU * kvhs(ik,ir) * osmrec(ik,ir)
c              osmmp(ik,ir) = osmmp(ik,ir) + 2.0 * pinmp(ik,ir)
            ENDDO
          ELSEIF (s28mom_.EQ.0) THEN
          ELSE
            STOP 'SOL28: S28MOM OPTION NOT SPECIFIED'
          ENDIF

        ENDIF

c *MOVE*
c... Eventually pass S28MOM_, PINMP, etc. so that there are not potential mixups...
        IF (ringtype(ir).EQ.SOL1) s28cfpnt_ = s28cfpnt
        IF (ringtype(ir).EQ.PFZ ) s28cfpnt_ = s28cfppfznt
c...    Another momentum source specified by outside considerations, so it
c       has to be included here, before the momentum balance sources are applied:
        IF (s28cfpnt_.EQ.5) THEN

          STOP 'CRAP CFPNT OPTION'

c  This is a very limited approach -- need to determine if the plasma is really
c  flowing towards the target, proper mass dependence, ...
          tmpmp = 0.0

          DO ik = ik0, (ik1+ik2)/2
            cs= GetCs(ktebs(ik,ir),ktibs(ik,ir))
            tmpmp(ik,ir) = 2.0 * AMU * kvhs(ik,ir)*osmrec(ik,ir)
          ENDDO

c...  
          DO ik = ik0, ik6
            osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
          ENDDO  
        ENDIF


c      IF (.TRUE.) THEN
c...    Specify momentum source:

        IF (detach1.AND.
     .      (s28mom_.LE.0.OR.s28mom_.EQ.5)) THEN

          IF (.TRUE..OR.p0.LT.p3) THEN

            tmpmp = 0.0

            IF (count.GT.0) THEN
c...          Find region where Te < 10 eV and apply uniform momentum loss:

              iks = 1
              DO ike = iks, nks(ir)-1
                IF (ktebs(ike+1,ir).GT.s28momTe) EXIT
c                IF (ktebs(ike+1,ir).GT.10.0) EXIT
              ENDDO
c              IF (ir.GT.irtrap) ike = ik3

              IF (s28mom_.NE.-3.AND.s28mom_.NE.-4.AND.
     .            ike.EQ.nks(ir)) THEN

                IF (ringtype(ir).EQ.PFZ) THEN
c...              Switch to triangular source (for now):
                  IF (s28mom_.NE.-1.AND.s28mom_.NE.-2.AND.
     .                s28mom_.NE.5 ) s28mom_ = -1

                  ike = ik3
                ELSE
                  STOP 'PROBLEMO A'
                ENDIF

              ELSEIF (s28mom_.EQ.-1.AND.ringtype(ir).EQ.PFZ.AND.
     .                ir.EQ.43) THEN
c...TMP: for compatability with PSI cases (d-hds-001*.d6i)
                  ike = ik3
              ENDIF
c              WRITE(0,*) '>>>> MOM IKE=',ike,kss(iks,ir),count,s28mom_


              IF     (s28mom_.EQ.-3.OR.

     .                (ringtype(ir).EQ.SOL1.AND.s28mom_.EQ.5).OR.

     .                ((ir.EQ.70.OR.ir.EQ.71).AND.
     .                 ringtype(ir).EQ.PFZ.AND.
     .                 (s28mom_.EQ.-1.OR.s28mom_.EQ.-2.OR.
     .                  s28mom_.EQ.5 ))) THEN

c...            Uniform with triangluar source near target:



                IF (ringtype(ir).EQ.SOL1) THEN
                  iks = 1
                  DO ike = iks, nks(ir)-1
                    IF (kss(ike,ir).GT.0.15*ksmaxs(ir)) EXIT
                  ENDDO
                ENDIF

                DO ik = iks, ike
                  tmpmp(ik,ir) = 1.0
                ENDDO             

                alpha = 5.0
c                IF (ir.EQ.70) alpha = 1.0
                IF (ir.EQ.70) alpha = 40.0
                IF (ir.EQ.71) alpha = 20.0

                DO ik = ik0, ik1
                  tmpmp(ik,ir) = tmpmp(ik,ir) + 
     .                           alpha * (kss(ik ,ir) - kss(ik1,ir)) / 
     .                                   (kss(ik0,ir) - kss(ik1,ir))
                ENDDO             

              ELSEIF (s28mom_.EQ.-4) THEN

c....           Uniform betwen ik4 and ik5:
                iks = ik1
                ike = ik2
                DO ik = iks, ike
                  tmpmp(ik,ir) = 1.0
                ENDDO


              ELSEIF (s28mom_.EQ.-2.OR.
     .                s28mom_.EQ.5 ) THEN
c...            Exponential decay:
                IF (.FALSE.) THEN
c...CMOD:
                  IF (ir.GT.irtrap) alpha = 0.07
                  IF (ir.EQ.38) alpha = 2.0
                  IF (ir.EQ.37) alpha = 2.0
                  IF (ir.EQ.36) alpha = 2.0
                  IF (ir.EQ.35) alpha = 0.9
                  IF (ir.EQ.34) alpha = 0.7
                  IF (ir.EQ.33) alpha = 0.6
                  IF (ir.EQ.32) alpha = 0.5
                  IF (ir.GE.27.AND.ir.LE.31) alpha = 0.07

                  IF (ir.LT.irwall) alpha = 1.0  
                  IF (ir.EQ.15) alpha = 0.5   ! 1.0
                  IF (ir.EQ.16) alpha = 1.0
                  IF (ir.EQ.17) alpha = 1.0
                  IF (ir.EQ.18) alpha = 1.0
                  IF (ir.EQ.19) alpha = 0.5
                  IF (ir.EQ.20) alpha = 0.100 !0.075
                  IF (ir.EQ.21) alpha = 0.01
                  IF (ir.EQ.22) alpha = 0.001
                  IF (ir.EQ.23) alpha = 0.01
                  IF (ir.EQ.24) alpha = 1.0
                ELSE
c...DIIID:
                  alpha = 0.1
c                  alpha = 0.5
c                  alpha = 2.0
                  IF (ir.EQ.70) alpha = 0.2
                  IF (ir.EQ.71) alpha = 0.1
                ENDIF

                A = 1.0 / (1.0 - EXP(-kss(ike,ir)/alpha))
                B = 1.0 - A
                DO ik = iks, ike
                  sval = kss(ik,ir)
                  tmpmp(ik,ir) = A * EXP(-sval/alpha) + B
c                  WRITE(0,*) '--->',tmpmp(ik,ir)
                ENDDO




              ELSEIF (s28mom_.EQ.-1.OR.s28mom_.EQ.2) THEN
c...            Triangular source:
                DO ik = iks, ike
                  tmpmp(ik,ir) = (kss(ik ,ir) - kss(ike,ir)) / 
     .                           (kss(iks,ir) - kss(ike,ir))
                ENDDO             
              ELSE
                DO ik = iks, ike
                  tmpmp(ik,ir) = 1.0
                ENDDO             
              ENDIF
            ELSE
              DO ik = ik0, ik2
                tmpmp(ik,ir) = 1.0
              ENDDO
            ENDIF
c *NEED TO ADD MOMENTUM LOSS HERE*

            IF (n(1).NE.0.0) THEN

              IF (ike.NE.0.AND.ike.LT.ik1) THEN
                IF (.NOT.warn1) THEN
                  WRITE(0,*) '***** MP PROBLEM A? *****'
                ENDIF
              ENDIF

              CALL CalcIntegral3(osmmp,ik0,ik1-1,ir,mom01,2)
c              CALL CalcIntegral3(osmmp,ik1,ik3  ,ir,mom13,2)  ! dev

c... set v ... problems if ike is less than ik1
              p1 = CalcPressure(n(1),te(1),ti(1),0.0)*ECH
c              p1 = CalcPressure(n(1),te(1),ti(1),ABS(kvhs(ik1,ir)))*ECH
c              v1 = MIN(ABS(kvhs(ik1,ir)),
c     .                 GetCs(ktebs(ik1,ir),ktibs(ik1,ir)))
c              p1 = CalcPressure(n(1),te(1),ti(1),v1)*ECH


              IF (s28mom_.EQ.-4) THEN

c... What a mess, needs to be overhauled:

                CALL CalcIntegral3(osmmp,ik1,ik3,ir,mom13,2)
                CALL Normalize(ik1,ik3,ir,tmpmp,DBLE(p3-p1-mom13),0,'8')

                DO ik = ik0, ik1-1
                  tmpmp(ik,ir) = 1.0
                ENDDO  

              ELSE

c                tmpmp = 0.0

                DO ik = ik1, ik2
                  tmpmp(ik,ir) = 0.0
                ENDDO  
                DO ik = ik2+1, ik3
                  tmpmp(ik,ir) = 1.0
                ENDDO  

              ENDIF



             CALL Normalize(ik0,ik1-1,ir,tmpmp,DBLE(p1-p0-mom01),0,'FA')
c             CALL Normalize(ik1,ik3  ,ir,tmpmp,DBLE(p3-p1-mom13),0,'FB')

              DO ik = ik0, ik6
                osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
              ENDDO  


c              WRITE(0,*) 'MOM:',n(1),te(1),ti(1)              
c              WRITE(0,*) 'MOM:',n(3),te(3),ti(3)      
        
c              WRITE(0,*) 'MOM:',ABS(kvhs(ik1,ir))/GetCs(te(1),ti(1))
c              WRITE(0,*) 'MOM:',p0/ECH,p1/ECH,p3/ECH
c              WRITE(0,*) 'MOM:',mom01/ECH,mom13/ECH
c              CALL CalcIntegral3(osmmp,ik0  ,ik1,ir,mom01,4)
c              CALL CalcIntegral3(osmmp,ik1+1,ik3,ir,mom13,4)
c              WRITE(0,*) 'MOM:',mom01/ECH,mom13/ECH
c              WRITE(0,*) 'MOM:',p3-p1,p1-p0
c              STOP 'sdfsddsf'

            ELSE

              CALL Normalize(ik0,ik6,ir,tmpmp,DBLE(p3-p0),0,'CF7')
              DO ik = ik0, ik6
                osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
              ENDDO  

c                DO ik = iks, ike
c                  WRITE(0,*) 'sdsdaa',ik,ir,tmpmp(ik,ir)
c                ENDDO
       

            ENDIF

c            WRITE(0,*) 'MOM->',p3,p0
          ENDIF
        ENDIF


        IF (detach2.AND.
     .      (s28mom_.LE.0.OR.s28mom_.EQ.5)) THEN

          mom36 = 0.0

          IF (.TRUE..OR.p6.LT.p3) THEN

            tmpmp = 0.0
c... 
            IF (count.GT.0) THEN

c...          Find region where Te < 10 eV and apply uniform momentum loss:
              ike = nks(ir)
              DO iks = ike, 2, -1
                IF (ktebs(iks-1,ir).GT.s28momTe) EXIT
c                IF (ktebs(iks-1,ir).GT.10.0) EXIT
              ENDDO
c              IF (ir.GT.irtrap) iks = ik3 + 1
              IF ((s28mom_.NE.-3.AND.s28mom_.NE.-4).AND.iks.EQ.1) THEN
c              IF (iks.EQ.0) THEN
                IF (ringtype(ir).EQ.PFZ) THEN
c...              Switch to triangular source (for now):
                  IF (s28mom_.NE.-1.AND.s28mom_.NE.-2.AND.
     .                s28mom_.NE.5 ) s28mom_ = -1
                  iks = ik3 + 1
                ELSE
                  STOP 'PROBLEMO B'
                ENDIF
              ELSEIF (s28mom_.EQ.-1.AND.ringtype(ir).EQ.PFZ.AND.
     .                ir.EQ.43) THEN
c...TMP: for compatability with PSI cases (d-hds-001*.d6i)
                  iks = ik3 + 1
              ENDIF
c              WRITE(0,*) '>>>> MOM IKS=',iks,kss(iks,ir),count,s28mom_
 
              IF     (s28mom_.EQ.5.AND.
     .                (ringtype(ir).EQ.SOL1.OR.
     .                 ir.EQ.70.OR.ir.EQ.71)) THEN

c... size of source needed...
                CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,2)

                IF (ringtype(ir).EQ.SOL1) THEN
                  ike = ik4
                  DO iks = ike, ik3+1, -1
                    IF (kss(iks,ir).LT.0.90*ksmaxs(ir)) EXIT
                  ENDDO                 
                ELSE
                  iks = ik3+1
                  ike = ik4
                ENDIF
                DO ik = iks, ike
                  tmpmp(ik,ir) = 1.0
                ENDDO

c                alpha = 5.0
c                IF (ir.EQ.70) alpha = 40.0
c                IF (ir.EQ.71) alpha = 20.0
c                DO ik = ik5, ik6
c                  tmpmp(ik,ir) = tmpmp(ik,ir) + 
c     .                           alpha * (kss(ik ,ir) - kss(ik5,ir)) / 
c     .                                   (kss(ik6,ir) - kss(ik5,ir))
c                ENDDO             

              ELSEIF (s28mom_.EQ.-4) THEN

c....           Uniform betwen ik4 and ik5:
                iks = ik4
                ike = ik5
                DO ik = iks, ike
                  tmpmp(ik,ir) = 1.0
                ENDDO

              ELSEIF (s28mom_.EQ.-3.OR.
     .                ((ir.EQ.70.OR.ir.EQ.71).AND.
     .                 (s28mom_.EQ.-1.OR.s28mom_.EQ.-2))) THEN
c...            Uniform with triangluar source near target:

                IF (ir.LT.irwall) THEN
                  ike = nks(ir)
                  DO iks = ike, 2, -1
                    IF (kss(iks,ir).LT.0.90*ksmaxs(ir)) EXIT
                  ENDDO                 
                ENDIF
                DO ik = iks, ike
                  tmpmp(ik,ir) = 1.0
                ENDDO

                alpha = 3.0
                IF (ir.GE.25.AND.ir.LT.irwall) alpha = 10.0
                IF (ir.EQ.70) alpha = 15.0

c...            Triangular source:
                DO ik = ik5, ik6
                  tmpmp(ik,ir) = tmpmp(ik,ir) + 
     .                         alpha * (kss(ik ,ir) - kss(ik5,ir)) / 
     .                                 (kss(ik6,ir) - kss(ik5,ir))
                ENDDO             


              ELSEIF (s28mom_.EQ.-2.OR.
     .                s28mom_.EQ.5 ) THEN
c...            Exponential decay:

                IF (.FALSE.) THEN
c...CMOD:
                  alpha = 5.0  
                ELSE
c...DIIID:
                  alpha = 0.1
c                  alpha = 0.5
                  IF (ir.EQ.66) alpha = 0.05
                  IF (ir.EQ.70) alpha = 0.20

                ENDIF
                
c *LATER*
c                IF (ir.EQ.30) ike = nks(ir) - 15

                A = 1.0 / (1.0 - EXP(-(ksb(ike,ir)-
     .                                 kss(iks,ir))/alpha))
                B = 1.0 - A
                DO ik = iks, ike
                  sval = ksb(ike,ir) - kss(ik,ir)
                  tmpmp(ik,ir) = A * EXP(-sval/alpha) + B
c                  WRITE(0,*) '--->',ik,sval,tmpmp(ik,ir)
                ENDDO

c                STOP 'sdfsd'                  

              ELSEIF (s28mom_.EQ.-2.OR.
     .                s28mom_.EQ.-1.OR.s28mom_.EQ.2) THEN
c...            Triangular source:
                DO ik = iks, ike
                  tmpmp(ik,ir) = (kss(ik ,ir) - kss(iks,ir)) / 
     .                           (kss(ike,ir) - kss(iks,ir))
                ENDDO             

              ELSE
                DO ik = iks, ike
                  tmpmp(ik,ir) = 1.0
                ENDDO             
              ENDIF

            ELSE
              DO ik = ik4, ik6
                tmpmp(ik,ir) = 1.0
              ENDDO
            ENDIF

c *NEED TO ADD MOMENTUM LOSS HERE*

c            WRITE(0,*) '---------------------------------------------'
c            WRITE(0,*) 'N5:',n(5)



            IF (n(5).NE.0.0) THEN

c              WRITE(0,*) 'MOM:',ik5,ike,iks
              IF (iks.NE.0.AND.iks.GT.ik5) THEN
c                tmpmp(ik3+1:MAXNKS,ir) = 0.0
                IF (.NOT.warn1) THEN
                  WRITE(0,*) '***** MP PROBLEM B? *****'
                  warn1 = .TRUE.
                ENDIF
c                STOP 'MP PROBLEM B?'
              ENDIF

c              CALL CalcIntegral3(osmmp,ik3+1,ik5,ir,mom35,2)
              CALL CalcIntegral3(osmmp,ik5+1,ik6,ir,mom56,2)
c... set v ... problems if iks is greater than ik5
c              p5 = CalcPressure(n(5),te(5),ti(5),0.0)*ECH
              v5 = MIN(ABS(kvhs(ik5,ir)),
     .                 GetCs(ktebs(ik5,ir),ktibs(ik5,ir)))
              p5 = CalcPressure(n(5),te(5),ti(5),v5)*ECH

              IF (sloutput) THEN
                WRITE(0,*) 'MOM P5:',p5/ECH,n(5),te(5),ti(5)
                WRITE(0,*) 'MOM P5:',kvhs(ik5,ir),mom56
              ENDIF

c              tmpmp(ik0:ik6,ir) = 0.0
              IF (s28mom_.EQ.-4) THEN

c... What a mess, needs to be overhauled:

                CALL CalcIntegral3(osmmp,ik3,ik5,ir,mom35,2)
                CALL Normalize(ik4,ik5,ir,tmpmp,DBLE(p5-p3-mom35),0,'8')

                DO ik = ik5+1, ik6
                  tmpmp(ik,ir) = 1.0
                ENDDO  

              ELSE

c                tmpmp = 0.0
                DO ik = ik3, ik4-1
                  tmpmp(ik,ir) = 1.0
                ENDDO  
                DO ik = ik4, ik5
                  tmpmp(ik,ir) = 0.0
                ENDDO  

              ENDIF

c             CALL Normalize(ik3+1,ik5,ir,tmpmp,DBLE(p5-p3-mom35),0,'F1')
              CALL Normalize(ik5+1,ik6,ir,tmpmp,DBLE(p6-p5-mom56),0,'F')

              DO ik = ik0, ik6
                osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
              ENDDO  

c              WRITE(0,*) 'MOM:',n(5),te(5),ti(5)              
c              WRITE(0,*) 'MOM:',p3/ECH,p5/ECH,p6/ECH
c              WRITE(0,*) 'MOM:',mom35/ECH,mom56/ECH
              CALL CalcIntegral3(osmmp,ik3+1,ik5,ir,mom35,4)
              CALL CalcIntegral3(osmmp,ik5  ,ik6,ir,mom56,4)
c              WRITE(0,*) 'MOM1:',mom56/ECH
c              WRITE(0,*) 'MOM2 P3:',p5/ECH,(p5-p3)/ECH,(p6-p5)/ECH
c              WRITE(0,*) 'MOM2 P3:',ABS(kvhs(ik5,ir))/
c     .                              GetCs(te(5),ti(5))
c              STOP 'sdfsddsf'

            ELSE

              IF (p3.EQ.0.0) STOP 'NO UPSTREAM PRESSURE B'

c...          Pressure loss based on upstream value:
              CALL Normalize(ik0,ik6,ir,tmpmp,DBLE(p6-p3-mom36),0,'CF8')
              DO ik = ik0, ik6
                osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
              ENDDO  

            ENDIF

          ENDIF

        ENDIF






c...    Pressure conservation:

c *** REALLY WANT THIS BETWEEN IK3 AND THE TARGETS!

c *MOVE*
        s28cfm_ = 99
        IF (ringtype(ir).EQ.SOL1) s28cfm_ = s28cfm
        IF (ringtype(ir).EQ.PFZ ) s28cfm_ = s28cfmpfz
c
c       S28CFM -1 - 
c               0 - uniform source applied between targets
c               1 - 
c               2 - 
c
        IF (.FALSE..AND..TRUE..AND.
     .      ringtype(ir).EQ.SOL1.AND.
     .      (.NOT.detach1.AND..NOT.detach2)) THEN

          IF (n(3).NE.0.0) THEN

            CALL CalcIntegral3(osmmp,ik0  ,ik3,ir,mom03,2)
            mom03 = p3 - p0 - mom03
            tmpmp = 0.0
            DO ik = ik0, ik3  
              tmpmp(ik,ir) = 1.0
            ENDDO
            CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom03),'mp6')
            DO ik = 1, nks(ir)
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

            CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,2)
            mom36 = p6 - p3 - mom36
            tmpmp = 0.0
            DO ik = ik3+1, ik6
              tmpmp(ik,ir) = 1.0
            ENDDO
            CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom36),'mp6')
            DO ik = 1, nks(ir)
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ELSE

            CALL CalcIntegral3(osmmp,ik0,ik6,ir,mom06,2)
            mom06 = p6 - p0 - mom06
            tmpmp = 0.0
            DO ik = 1, nks(ir)
              tmpmp(ik,ir) = 1.0
            ENDDO
            CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom06),'mp6')
            DO ik = 1, nks(ir)
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
c              WRITE(0,*) 'MOM:',ik,ir,osmmp(ik,ir),mom06
            ENDDO  

          ENDIF 

        ELSEIF (.FALSE..AND.ringtype(ir).EQ.PFZ) THEN

          IF (.NOT.detach1.AND..NOT.detach2) THEN

            IF (n(3).NE.0.0) THEN

              CALL CalcIntegral3(osmmp,ik0  ,ik3,ir,mom03,2)
              mom03 = p3 - p0 - mom03
              tmpmp = 0.0
              DO ik = ik0, ik3  
                tmpmp(ik,ir) = 1.0
              ENDDO
              CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom03),'mp6')
              DO ik = 1, nks(ir)
                osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
              ENDDO  

              CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,2)
              mom36 = p6 - p3 - mom36
              tmpmp = 0.0
              DO ik = ik3+1, ik6
                tmpmp(ik,ir) = 1.0
              ENDDO
              CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom36),'mp6')
              DO ik = 1, nks(ir)
                osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
              ENDDO  

            ELSE

              CALL CalcIntegral3(osmmp,ik0,ik6,ir,mom06,2)
              mom06 = p6 - p0 - mom06
              tmpmp = 0.0
              DO ik = 1, nks(ir)
                tmpmp(ik,ir) = 1.0
              ENDDO
              CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom06),'mp6')
              DO ik = 1, nks(ir)
                osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
              ENDDO  
            ENDIF           

          ELSEIF (.NOT.detach2) THEN

            CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,2)
            mom36 = p6 - p3 - mom36
            tmpmp = 0.0
            DO ik = ik3+1, ik6
              tmpmp(ik,ir) = 1.0
            ENDDO
            CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom36),'mp6')
            DO ik = 1, nks(ir)
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ENDIF

        ELSEIF (ringtype(ir).EQ.PFZ.OR.ringtype(ir).EQ.SOL1) THEN

          iks = ik0
          ike = ik6
          IF (detach1) iks = ik2 + 1
          IF (detach2) ike = ik4 - 1

          IF     (s28mom_.EQ.20) THEN
c...        Do nothing:
          ELSEIF (p3.EQ.0.0) THEN

            CALL CalcIntegral3(osmmp,ik0,ik6,ir,mom06,2)
            mom06 = p6 - p0 - mom06
            tmpmp = 0.0
            DO ik = iks, ike
              tmpmp(ik,ir) = 1.0
            ENDDO  
            CALL Normalize(iks,ike,ir,tmpmp,DBLE(mom06),0,'FA')
            DO ik = ik0, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ELSE

            CALL CalcIntegral3(osmmp,ik0,ik3,ir,mom03,2)
            mom03 = p3 - p0 - mom03
            tmpmp = 0.0
            DO ik = iks, ik3
              tmpmp(ik,ir) = 1.0
            ENDDO  
            CALL Normalize(ik0,ik3,ir,tmpmp,DBLE(mom03),0,'FA')  
            DO ik = ik0, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

            CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,2)
            mom36 = p6 - p3 - mom36
            tmpmp = 0.0
            DO ik = ik3+1, ike
              tmpmp(ik,ir) = 1.0
            ENDDO  
            CALL Normalize(ik3+1,ik6,ir,tmpmp,DBLE(mom36),0,'FA')
            DO ik = ik0, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ENDIF

        ELSEIF (.FALSE..AND.
     .          ringtype(ir).EQ.PFZ.OR.ringtype(ir).EQ.SOL1) THEN

          IF     (.NOT.detach1) THEN

            CALL CalcIntegral3(osmmp,ik0,ik3,ir,mom03,2)
            mom03 = p3 - p0 - mom03
            tmpmp = 0.0
            DO ik = ik0, ik3
              tmpmp(ik,ir) = 1.0
            ENDDO  
            CALL Normalize(ik0,ik3,ir,tmpmp,DBLE(mom03),0,'FA')
            DO ik = ik0, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ELSEIF (detach1) THEN

            CALL CalcIntegral3(osmmp,ik0,ik3,ir,mom03,2)
            mom03 = p3 - p0 - mom03
            tmpmp = 0.0
            DO ik = ik2+1, ik3
              tmpmp(ik,ir) = 1.0
            ENDDO  
            CALL Normalize(ik0,ik3,ir,tmpmp,DBLE(mom03),0,'FA')
            DO ik = ik0, ik6
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ENDIF

          IF     (.NOT.detach2) THEN

            CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,2)
            mom36 = p6 - p3 - mom36
            tmpmp = 0.0
            DO ik = ik3+1, ik6
              tmpmp(ik,ir) = 1.0
            ENDDO
            CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom36),'mp6')
            DO ik = 1, nks(ir)
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ELSEIF (detach2) THEN

            CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,2)
            mom36 = p6 - p3 - mom36
            tmpmp = 0.0
            DO ik = ik3+1, ik4-1
              tmpmp(ik,ir) = 1.0
            ENDDO
            CALL NormalizeInitSrc(4,ir,tmpmp,DBLE(mom36),'mp6')
            DO ik = 1, nks(ir)
              osmmp(ik,ir) = osmmp(ik,ir) + tmpmp(ik,ir) 
            ENDDO  

          ENDIF

        ELSE
          WRITE(0,*) 'MOMENTUM: WHAT TO DO?',ir,ringtype(ir),SOL1,PFZ
          STOP 'klkjhd'
        ENDIF



        CALL CalcIntegral3(osmmp,ik0  ,ik3,ir,mom13,4)
        CALL CalcIntegral3(osmmp,ik3+1,ik6,ir,mom36,4)

        IF (sloutput.AND.ir.NE.irlast) THEN
c          WRITE(0,*) 'MOMENTUM:',(p0-p3)-mom13,p0,p3,mom13
c          WRITE(0,*) 'MOMENTUM:',(p6-p3)-mom36,p3,p6,mom36
          irlast = ir
        ENDIF

c
c * RETURN *
c
        RETURN
      ENDIF


      RETURN
 99   STOP
      END
c
c
c
c
      SUBROUTINE MachProfile(ir,teprofile)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      use mod_slcom_sol28
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'

      INTEGER ir
      REAL    teprofile(MAXNKS)

      REAL    GetMach

      INTEGER ik,iks1,ike1,iks2,ike2
      REAL    machpro(MAXNKS),M0,Mx,x0,xx,frac,
     .        ion1x,rec1x,cfp1x,flx1x,mom1x
      REAL*8  alpha,beta,M,mass,gamma

      REAL
     .        a1,b1,c1,r1

      teprofile = 0.0
      machpro = 0.0


c...  Find range of KTEBS values to modify:
      iks1 = 0
      ike1 = 0
      iks2 = 0
      ike2 = 0
      DO ik = 1, nks(ir)
        IF (.NOT.detach1.AND.ik.LE.ik3) THEN
          IF (ktebs(ik,ir).EQ.te(0).AND.iks1.EQ.0) iks1 = ik
          IF (ktebs(ik,ir).EQ.te(0).AND.iks1.NE.0) ike1 = ik
        ENDIF
        IF (.NOT.detach2.AND.ik.GT.ik3) THEN
          IF (ktebs(ik,ir).EQ.te(6).AND.iks2.EQ.0) iks2 = ik
          IF (ktebs(ik,ir).EQ.te(6).AND.iks2.NE.0) ike2 = ik
        ENDIF
      ENDDO

c...  Assign Mach profile:
      
c...  Low IK index:
      IF (.NOT.detach1.AND.iks1.NE.0) THEN
        M0 = 1.0
        x0 = 0.0
        Mx = GetMach(kvhs(ike1+1,ir),ktebs(ike1+1,ir),ktibs(ike1+1,ir))
        xx = kss(ike1+1,ir)
        DO ik = iks1, ike1
          frac = (kss(ik,ir) - x0) / (xx - x0)
          machpro(ik) = M0 + frac * (Mx - M0)
        ENDDO
      ENDIF
c...  High IK index:
      IF (.NOT.detach2.AND.iks2.NE.0) THEN
        M0 = GetMach(kvhs(iks2-1,ir),ktebs(iks2-1,ir),ktibs(iks2-1,ir))
        x0 = kss(iks2-1,ir)
        Mx = 1.0
        xx = ksmaxs(ir)
        DO ik = iks2, ike2
          frac = (kss(ik,ir)- x0) / (xx - x0)
          machpro(ik) = M0 + frac * (Mx - M0)
        ENDDO
      ENDIF

       
c      DO ik = 1, nks(ir)
c        WRITE(0,*) 'MACHPRO:',ik,machpro(ik)
c      ENDDO
c      STOP 'sgsdsdg'

c...  Low IK index target:
      DO ik = 1, nks(ir)
        IF (machpro(ik).EQ.0.0) THEN
          teprofile(ik) = ktebs(ik,ir)
        ELSE

          CALL CalcIntegral3(osmion,1,ik,ir,ion1x,8)
          CALL CalcIntegral3(osmrec,1,ik,ir,rec1x,8)
          CALL CalcIntegral3(osmcfp,1,ik,ir,cfp1x,8)
          alpha = DBLE(ion1x + flx0 - rec1x + cfp1x)**2 * DBLE(ECH)
          flx1x = ion1x + flx0 - rec1x + cfp1x

          CALL CalcIntegral3(osmmp,1,ik,ir,mom1x,8)
          beta = DBLE(p0 + mom1x)**2 

          gamma = DBLE(machpro(ik) + 1.0 / machpro(ik))**2

          mass = DBLE(crmb * AMU) 
     

          teprofile(ik) = SNGL(beta / (2.0D0*alpha*gamma*mass))


c          WRITE(0,'(A,I6,F10.2,1P,4D10.2,0P)') 
c     .      'TEPRO:',ik,teprofile(ik),alpha,beta,gamma,mass


c...Test:
              a1 = (2*teprofile(ik)) * ECH
              b1 = -(p0 + mom1x)
              c1 = AMU * crmb * flx1x**2
              r1 = b1**2 - 4.0 * a1 * c1
              knbs(ik,ir) = (-b1 + SQRT(r1)) / (2.0 * a1)

              kvhs(ik,ir) = flx1x / knbs(ik,ir)

              mach1 = GetMach(kvhs(ik,ir),teprofile(ik),teprofile(ik))

              IF (ik.EQ.1) THEN
                WRITE(0,'(A,I6,2E10.2,1P,2E10.2,0P)') 
     .            'DATA:',ik,teprofile(ik),teprofile(ik),flx1x,mom1x
              ENDIF
             
c          WRITE(0,'(A,I6,2F10.4,1P,E10.2,0P)')
c     .      'TEST :',ik,machpro(ik),mach1,flx1x

                 

        ENDIF
      ENDDO

c      STOP 'sdfsd'

      RETURN
 99   STOP
      END

c
c ======================================================================
c
      SUBROUTINE CalcTeProfile(ir,ikflat1,ikflat2,machmode,count)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      use mod_slcom_sol28
      IMPLICIT none

      INTEGER ir,ikflat1,ikflat2,machmode,count

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'  
  
      INTEGER ik,iks,ike,ikm,ikb,temode,icount,ikadjust,ikmin,ikmax
      LOGICAL cont,bottom,final,output
      REAL    x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6,y2p,yA,yB,xB,
     .        aA,aB,bA,bB,aC,bC,deltat,frac,L,Psol,k,f,
     .        sval,texpt2,
     .        teprofile(MAXNKS),te1,te2,tres,tequ,cs,holdknbs(MAXNKS),
     .        dPsol,qconv,qcond,qcf,tgrad,qei,qrad,qpine,mPsol,
     .        flx1x,ion1x,rec1x,cfp1x,nrms,rms,fPsol,qtarg,
     .        mPsoladjust,fPsoladjust,smin,smax,smid

        x0 = s(0)
        x1 = s(1)
        x2 = s(2)
        x3 = s(3)
        x4 = s(4)
        x5 = s(5)
        x6 = s(6)

c...    Calculate temperature profile:
        y0 = te(0)
        y1 = te(1)
        y2 = te(2)  
        y3 = te(3)   
        y4 = te(4)   
        y5 = te(5)   
        y6 = te(6)   

        bA = 2.0
        aA =  (y2 - y1) / (x2 - x1)**bA

        y2p = -bB * aB * (x2 - x3)**(bB-1.0)


c...Conduction calculation of Te (used for PSI): 
c        temode = 2
c        IF (ringtype(ir).EQ.SOL1.AND.detach2) temode = 5

        temode = 3
        IF (ringtype(ir).EQ.PFZ             ) temode = 2

      IF (.NOT..TRUE..AND.irsep.EQ.10.AND.                ! **** HACK **** OFF?
     .    (ir.GE.10.AND.ir.LE.17)) THEN

c        WRITE(0,*) 'CHEAT!',techeat

        
        y6 =  techeat

      ENDIF



        DO ik = 1, nks(ir)

          te1 = 0.0
          te2 = 0.0


          IF (.FALSE.) THEN
          ELSE
c
c OLDIES!
c

            IF     (kss(ik,ir).LT.s(1)-0.0001) THEN

c              frac = (kss(ik,ir) - x5) / (x6 - x5)       
              frac = (kss(ik,ir) - x0) / (x1 - x0)
              IF (stopopt2.EQ.31) THEN               ! **** HACK ****
                frac = frac**10
c              ELSE
c                frac = 1.0
              ENDIF
c              ktebs(ik,ir) = y5 + frac * (y6 - y5)
              ktebs(ik,ir) = y0 + frac * (y1 - y0)

              ktibs(ik,ir) = ktebs(ik,ir)
c              WRITE(0,'(A,2I6,5F10.2,1P,2E10.2,0P)')
c     .          '1:',ik,ir,ktebs(ik,ir),te1,te2,y0,y1,sval,osmion(ik,ir)

c              ktebs(ik,ir) = y1 + (y0 - y1) *
c     .                       ((s(1) - kss(ik,ir))/kss(ik1,ir))**2.0
c              ktebs(ik,ir) = y1
            ELSEIF (kss(ik,ir).LE.s(2)) THEN

c              frac = (kss(ik,ir) - x1) / (x2 - x1)
c              ktebs(ik,ir) = y1 + frac * (y2 - y1)

              bC = 3.0
              aC = (y2 - y1) / ABS(x2 - x1)**bC
              ktebs(ik,ir) = aC * ABS(kss(ik,ir) - x1)**bC + y1

c              ktebs(ik,ir) =  aA * (kss(ik,ir) - x1)**bA + y1

              ktibs(ik,ir) = ktebs(ik,ir)
c              WRITE(0,'(A,2I6,5F10.2,1P,2E10.2,0P)')
c     .          '2:',ik,ir,ktebs(ik,ir),te1,te2,y1,y2,sval,osmion(ik,ir)

            ELSEIF (kss(ik,ir).LE.s(3)) THEN

c...          Attached plasma:
              IF     (machmode.EQ.3) THEN
                ktebs(ik,ir) = teprofile(ik)                  
c...
              ELSEIF (machmode.EQ.2.AND.ik.LE.ikflat1) THEN
                ktebs(ik,ir) = te(0)
              ELSE
                IF     (temode.EQ.1) THEN
                  bB = 2.0
                  aB = -(y2 - y3) / (x3 - x2)**bB
                  sval = kss(ik,ir) - x3
                  ktebs(ik,ir) = -aB * sval**bB + y3
                ELSEIF (temode.EQ.2.OR.temode.EQ.5) THEN
                  IF (detach1) THEN
                    bB = s28b34det
                  ELSE
                    bB = s28b34
                  ENDIF
                  aB = (y3 - y2) / (x3 - x2)**bB
                  sval = kss(ik,ir) - x2
                  ktebs(ik,ir) =  aB * sval**bB + y2

                  ktibs(ik,ir) = ktebs(ik,ir)
                ELSEIF (temode.EQ.3) THEN

                  bB = 4.0
                  aB = -(y2 - y3) / (x3 - x2)**bB
                  sval = kss(ik,ir) - x3
                  te1 = -aB * sval**bB + y3

                  IF (detach1) THEN
                    bB = s28b34det
                  ELSE
                    bB = s28b34
                  ENDIF

                  aB = (y3 - y2) / (x3 - x2)**bB
                  sval = kss(ik,ir) - x2
                  te2 =  aB * sval**bB + y2

                  frac = sval / (x3 - x2)
c                  frac = frac ** 0.90
                  frac = frac ** 0.75

                  ktebs(ik,ir) = frac * te1 + (1.0 - frac) * te2

                  frac = frac ** 2.0
                  ktibs(ik,ir) = frac * ktebs(ik,ir) * ti(3) / te(3) +
     .                           (1.0 - frac) * ktebs(ik,ir) 

                ELSEIF (temode.EQ.4) THEN

                  IF (detach1) THEN
                    bB = s28b34det
                  ELSE
                    bB = s28b34
                  ENDIF
                  aB = (y3 - y2) / (x3 - x2)**bB
                  sval = kss(ik,ir) - x2
                  ktebs(ik,ir) =  aB * sval**bB + y2

                  ktibs(ik,ir) = ktebs(ik,ir)

                ENDIF
              ENDIF

c              WRITE(0,'(A,2I6,1P,5E10.2,0P)')
c     .          '3:',ik,ir,ktebs(ik,ir),aB,bB,y3,kss(ik,ir)-x3

            ELSEIF (kss(ik,ir).LE.s(4)) THEN

c...          Attached:
              IF     (machmode.EQ.3) THEN

                ktebs(ik,ir) = teprofile(ik)                  

              ELSEIF (machmode.EQ.2.AND.ik.GE.ikflat2) THEN
                ktebs(ik,ir) = te(6)
              ELSE
                IF     (temode.EQ.1) THEN
                  bC = 2.0
                  aC = -(y4 - y3) / (x4 - x3)**bC
                  sval = kss(ik,ir) - x3
                  ktebs(ik,ir) = -aC * sval**bC + y3
                  ktibs(ik,ir) = ktebs(ik,ir)
                ELSEIF (temode.EQ.2) THEN

                  IF (detach2) THEN
                    bC = s28b34det
                  ELSE
                    bC = s28b34
                  ENDIF
                  aC = (y3 - y4) / (x4 - x3)**bC
                  sval = x4 - kss(ik,ir) 
                  ktebs(ik,ir) =  aC * sval**bC + y4
                  ktibs(ik,ir) = ktebs(ik,ir)

                ELSEIF (temode.EQ.3) THEN
c                  bC = 3.0
                  bC = 4.0
                  aC = -(y4 - y3) / (x4 - x3)**bC
                  sval = kss(ik,ir) - x3
                  te1 = -aC * sval**bC + y3

                  IF (detach2) THEN
                    bC = s28b34det
                  ELSE
                    bC = s28b34
                  ENDIF

                  aC = (y3 - y4) / (x4 - x3)**bC
                  sval = x4 - kss(ik,ir) 
                  te2 =  aC * sval**bC + y4

                  frac = sval / (x4 - x3)
c                  frac = frac ** 2.00
                  frac = frac ** 0.75
c                  frac = 1.0

                  ktebs(ik,ir) = frac * te1 + (1.0 - frac) * te2

                  frac = frac ** 2.0
                  ktibs(ik,ir) = frac * ktebs(ik,ir) * ti(3) / te(3) +
     .                           (1.0 - frac) * ktebs(ik,ir) 
c                  ktibs(ik,ir) = ktebs(ik,ir)

c                  ktebs(ik,ir) = te2
c                  ktibs(ik,ir) = te1
 
c                  ktebs(ik,ir) = MAX(te1,te2)
c                  ktibs(ik,ir) = ktebs(ik,ir)

                ELSEIF (temode.EQ.4) THEN

                  IF (detach2) THEN
                    bC = s28b34det
                  ELSE
                    bC = s28b34
                  ENDIF

                  bC = 2.0 / 7.0

                  aC = (y3**(1.0/bC) - y4**(1.0/bC)) / (x3 - x4)
                  sval = kss(ik,ir) - x4
                  ktebs(ik,ir) =  (y4**(1.0/bC) + aC * sval)**bC

c                  L = x4 - x3
c                  aC = (y3**(1.0/bC) - y4**(1.0/bC)) / L**2.0
c                  sval = L - (x4 - kss(ik,ir))
c                  ktebs(ik,ir) =  (y4**(1.0/bC) + 
c     .                             aC * (L**2.0 - sval**2.0))**bC

c                  aC = (y3**(1.0/bC) - y4**(1.0/bC)) / (x3 - x4)**2.0
c                  sval = kss(ik,ir) - x4
c                  ktebs(ik,ir) =  (y4**(1.0/bC) + aC * sval**2.0)**bC


c                  aC = (y3 - y4) / (x4 - x3)**bC
c                  sval = x4 - kss(ik,ir)
c                  ktebs(ik,ir) =  y4 + aC * sval**bC

                  ktibs(ik,ir) = ktebs(ik,ir)

                ELSEIF (temode.EQ.5) THEN

c...              Need a better expression for location of x-point:
                  IF (ik.LT.ikti) THEN
   
                     k = 2000.0

c...                 Doubt this is good:
c                     L  = kss(ikti,ir) - kss(ik3,ir)
                     L  = kss(ik4,ir) - kss(ik3,ir)
c                     L  = ksmaxs(ir) - kss(ik3,ir)

                     Psol = 2.0 / 7.0 * te(3)**(7.0/2.0) * k / L 
c                     Psol = 4.0 / 7.0 * te(3)**(7.0/2.0) * k / L 

                     sval = kss(ik,ir) - kss(ik3,ir)

                     te1 = te(3)**(7.0/2.0) - 
     .                     (7.0/2.0) * Psol * sval**2.0 / (L * k) 
c     .                     (7.0/4.0) * Psol * sval**2.0 / (L * k) 
                     te1 = te1**(2.0/7.0)

                     ktebs(ik,ir) = te1

                     texpt2 = te1



c                     WRITE(0,*) 'TE:',ik,ikti,te1
c                     WRITE(0,*) 'PSOL:',Psol
c                     STOP

                     ktibs(ik,ir) = ktebs(ik,ir)

                  ELSE

                     L    = kss(ik4,ir) - kss(ikti-1,ir)
                     sval = kss(ik4,ir) - kss(ik,ir)

                  
c                     deltat = texpt2**(7.0/2.0) - te(4)**(7.0/2.0)
c                     f = ((4.0/7.0) * deltat * k / Psol / L) - 1.0

                     f = 0.10

                     te1 = te(4)**(7.0/2.0) + 
     .                     (7.0/2.0) / k * Psol * sval *
     .                     (sval / (2.0 * L) * (1.0 - f) + f)           

                     te2 = (7.0/2.0) / k * Psol * sval *
     .                     (sval / (2.0 * L) * (1.0 - f) + f)           



c                     L    = kss(ik4,ir) - kss(ikti-1,ir)
c                     sval = kss(ik ,ir) - kss(ikti-1,ir)

c                     f = (texpt2**(7.0/2.0) - te(4)**(7.0/2.0)) /
c     .                    (Psol * 7.0 * L / ( 2.0 * k * 2.0))
c                     te1 = texpt2**(7.0/2.0) - 
c     .                     (7.0/2.0) * Psol * f * (1.0 / k) *
c     .                     sval * (1.0 - sval / (2 * L))
c                     te1 = te1**(2.0/7.0)

c                     te1 = texpt2**(7.0/2.0) - 
c     .                     (7.0/4.0) * Psol * sval**2.0 / (L * k) 
c                     te1 = te1**(2.0/7.0)

                     te1 = te1**(2.0/7.0)
                     te2 = te2**(2.0/7.0)
                     ktebs(ik,ir) = te1
                     ktibs(ik,ir) = ktebs(ik,ir)

c                    WRITE(0,*) 'TE2:',ik,f,te1,te2
c                    WRITE(0,*) 'TE2:',ik,f,texpt2


c                    IF (detach2) THEN
c                      bC = s28b34det
c                    ELSE
c                      bC = s28b34
c                    ENDIF
c                    aC = (y3 - y4) / (x4 - x3)**bC
c                    sval = x4 - kss(ik,ir) 
c                    ktebs(ik,ir) =  aC * sval**bC + y4
c                    ktibs(ik,ir) = ktebs(ik,ir)                  
                  ENDIF


                ENDIF
              ENDIF

c              WRITE(0,'(A,2I6,6F10.2,1P,2E10.2,0P)')
c     .          '4:',ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     .          te1,te2,y3,y4,sval,osmion(ik,ir)
            ELSEIF (kss(ik,ir).LE.s(5)) THEN


c              frac = (kss(ik,ir) - x4) / (x5 - x4)
c              ktebs(ik,ir) = y4 + frac * (y5 - y4)
c              ktibs(ik,ir) = ktebs(ik,ir)

              bC = 3.0
              aC = (y4 - y5) / ABS(x4 - x5)**bC
              ktebs(ik,ir) = aC * ABS(kss(ik,ir) - x5)**bC + y5
              ktibs(ik,ir) = ktebs(ik,ir)

c              WRITE(0,'(A,2I6,1P,5E10.2,0P)')
c     .          '5:',ik,ir,ktebs(ik,ir),aB,bB,y3,kss(ik,ir)-x3

            ELSEIF (kss(ik,ir).LE.s(6)) THEN

              frac = (kss(ik,ir) - x5) / (x6 - x5)       
              IF (stopopt2.EQ.30) THEN               ! **** HACK ****
                frac = frac**0.1
              ELSE
                frac = 1.0
              ENDIF
              ktebs(ik,ir) = y5 + frac * (y6 - y5)

              ktibs(ik,ir) = ktebs(ik,ir)
c              WRITE(0,'(A,2I6,1P,5E10.2,0P)')
c     .          '6:',ik,ir,ktebs(ik,ir),aB,bB,y3,kss(ik,ir)-x3
            ELSE
              STOP 'ERROR'
            ENDIF
c           ktibs(ik,ir) = ktebs(ik,ir)

c *TEMP*
            IF (ik.GT.ik3) THEN
c              ktebs(ik,ir) = te(6)
c              ktibs(ik,ir) = te(6)
            ELSE
c              ktebs(ik,ir) = te(0)
c              ktibs(ik,ir) = te(0)
            ENDIF
 

          ENDIF

c          WRIT(E0,*) 'TE:',ik,ik3,ir,ktebs(ik,ir)

        ENDDO



c * NEW *

c...  Don't do this until count .GT. 0..?


      IF (ir.GE.irsep.AND.ir.LE.irsep+2.AND.
     .    iflexopt(4).EQ.41.AND..NOT.detach2) THEN

        icount = 0
        final = .FALSE.
        cont = .FALSE.
        output = .FALSE.
        IF (count.GT.0) cont = .TRUE.
        mPsol = 1.0
        fPsol = 0.50
        mPsoladjust = 0.005
        fPsoladjust = 0.10

        ikadjust = 0

        qtarg = (5.0 / 2.0 * te(6) * ECH) * (isat(6) / ECH)

        qtarg = (5.0 / 2.0 * (te(6) + ti(6)) * ECH + 
     .           0.5 * AMU * crmb * kvhs(1,ir)**2.0) * ! no mi convection for electron cannel?  
     .          (isat(6) / ECH)                        ! really need to solve Te and Ti together 

        DO WHILE (cont)
          cont = .FALSE.          

          nrms = 0.0
          rms = 0.0

          IF (detach2) THEN
            te1 = te(3)
            te2 = te(4)
            iks = ik3
            ikm = ikti + ikadjust - 1
            ike = ik4
            smin = kss(iks,ir)
            smax = kss(ike,ir)
            smid = kss(ikm,ir)
            L  = smax - smin

          ELSE
            te1 = te(3)
            te2 = te(6)
            iks = ik3
            ike = ik6+1
            DO ikm = ik3, ik6-1
              IF (thetag(ikm  ,ir).LT.thetag(ikti,irsep).AND.
     .            thetag(ikm+1,ir).GE.thetag(ikti,irsep)) EXIT
            ENDDO
            IF (ikm.EQ.ik6) THEN
              WRITE(0,*) 'TE: IKM NOT FOUND, CONDUCTION OFF'
              ike = iks-1
            ELSE
              ikm = ikm + ikadjust
              IF (output) WRITE(0,*) 'IKM:',ikm,ikti
            ENDIF
c            ikm = ikti + ikadjust - 1
            smin = kss(iks,ir)
            smax = ksb(ik6,ir)
            smid = kss(ikm,ir)
            L  = smax - smin
          ENDIF

c... Need flux limiters?  Worth the trouble in this instance?
          k = 2000.0
c... Factor 0.5 needed here, half to inner target, half to outer? How can I estimate asymmetries? 
          Psol = 4.0/7.0 * (te1**(7.0/2.0) - te2**(7.0/2.0)) * k/L 
c          Psol = 4.0/7.0 * (te1**(7.0/2.0) - te2**(7.0/2.0)) * k/L * 2.0
          Psol = Psol * mPsol

c...      Validate upstream power (move out of loop?):

          IF (output) THEN
            WRITE(0,*) 'fPsol:',fPsol,qtarg,Psol,mPsol
          ENDIF

          bottom = .FALSE.
          
          DO ik = iks, ike
c...      Need a better expression for location of x-point:
c...      Trouble for double null grids, or rings without a clear association with the x-point

            qconv = 0.0
            qei   = 0.0
            qpine = 0.0
            qrad  = 0.0
            qcf   = 0.0

            IF     (ik.LE.ikm) THEN
c...          Power entering SOL uniformly 
              frac = (kss(ik ,ir) - smin) / (smid - smin)
              dPsol = frac * Psol

              deltas = kss(ik,ir) - kss(ik-1,ir)

            ELSEIF (ik.LE.ik6) THEN
c...          Power leaving SOL uniformly 
              frac  = (kss(ik,ir) - smid) / (smax - smid)
              dPsol = Psol * (1.0 - frac * fPsol)

              deltas = kss(ik,ir) - kss(ik-1,ir)

c...          On the fly updates (can do better with the flux since it is known):
              CALL CalcIntegral3(osmion,1,ik,ir,ion1x,8)
              CALL CalcIntegral3(osmrec,1,ik,ir,rec1x,8)
              CALL CalcIntegral3(osmcfp,1,ik,ir,cfp1x,8)
              flx1x = ion1x + flx0 - rec1x + cfp1x
              qconv = (5.0 / 2.0 * ktebs(ik,ir) * ECH) * flx1x  ! Use te1 instead of ktebs? 
c              qconv = (5.0 / 2.0 * (ktebs(ik,ir) + ktibs(ik,ir)) * ECH + 
c     .                 0.5 * AMU * crmb * kvhs(ik,ir)**2.0) * ! no mi convection for electron cannel?  
c     .                flx1x                                   ! really need to solve Te and Ti together 

            ELSEIF (ik.EQ.ik6+1) THEN
c...          Last half of target cell:

              frac = 1.0
              dPsol = Psol * (1.0 - fPsol)
              flx1x = isat(6) / ECH
              qconv = (5.0 / 2.0 * te(6) * ECH) * flx1x

              deltas = ksb(ik6,ir) - kss(ik6,ir)

            ELSE
              STOP 'ERROR IN TE PROFILE'
            ENDIF

c...        Update Te:
            qcond = dPsol - qconv 
            tgrad = -qcond / (k * te1**(5.0/2.0))
            deltat = deltas * tgrad
            te1 = te1 + deltat

            IF (te1.LT.te2) THEN
              bottom = .TRUE.
              te1 = te2
            ENDIF

            IF (final) THEN
c...          Assign new Te profile:
              IF (bottom) THEN
                STOP 'NOT POSSIBLE!'
              ENDIF              
              IF (ik.LT.ik6+1) THEN
                ktebs(ik,ir) = te1
                ktibs(ik,ir) = te1
              ENDIF
            ENDIF

            IF (output)
     .        WRITE(0,'(A,3I6,2F8.4,1P,3E10.2,0P,4F10.4)')
     .          'Te :',ik,iks,ikm,frac,dPsol/Psol,
     .          dPsol,qconv,qcond,
     .          tgrad,deltat,ktebs(ik,ir),te1

            nrms = nrms + 1.0
            rms = rms + ABS(te1 - ktebs(ik,ir))**2
          ENDDO


c...      Evaluate quality of the Te profile:
c            ikadjust = ikadjust + 1
          IF     (fPsoladjust.LT.0.0.AND.     bottom) THEN
            fPsoladjust = -0.33 * fPsoladjust
          ELSEIF (fPsoladjust.GT.0.0.AND..NOT.bottom) THEN
            fPsoladjust = -0.33 * fPsoladjust
          ENDIF

          IF (output) THEN
            WRITE(0,*) '     :',te1,te2,ABS(te1-te2)
          ENDIF

          IF     (final) THEN
            IF (output) WRITE(0,*) 'LEAVING!'

c            STOP 'DONE!'
          ELSEIF (.NOT.bottom.AND.ABS(te1-te2).LT.0.1) THEN
            IF (output) WRITE(0,*) 'DONE!',icount
            final = .TRUE.
            cont = .TRUE.
c            output = .TRUE.

          ELSE
            fPsol = fPsol + fPsoladjust

            cont = .TRUE.

            IF (output) WRITE(0,*) 'mPsoladjust:',fPsoladjust,bottom

          ENDIF

          icount = icount + 1
          IF (icount.EQ.50) THEN
            WRITE(0,*) 'GOTTA BE A PROBLEM'
            STOP 'HALTING CODE'
          ENDIF



        ENDDO

      ENDIF




      IF (.NOT..TRUE..AND.irsep.EQ.10.AND.                ! **** HACK ****  OFF?
     .    (ir.GE.10.AND.ir.LE.17)) THEN

        ikmax = 0
        ikmin = MAXNKS
        DO ik = nks(ir), 1, -1
          IF (zs(ik,ir).GT.zxp) EXIT

          IF (zs(ik,ir).LT.-1.330) THEN
            ikmin = MIN(ik,ikmin)
            ikmax = MAX(ik,ikmax)
          ENDIF
        ENDDO

c        WRITE(0,*) 'ktebs cheater:',ikmin,ikmax,ir
c        te(6) = tehold
c        ti(6) = te(6)
        DO ik = ikmin, ikmax
          frac = (kss(ik   ,ir) - kss(ikmin-1,ir)) /
     .           (kss(ikmax,ir) - kss(ikmin-1,ir)) 

          frac = 1.0

          ktebs(ik,ir) = (1.0-frac) * ktebs(ikmin-1,ir) + frac * te(6)
c          WRITE(0,*) 'ktebs cheater:',ik,ir,ktebs(ik,ir),te(6)
        ENDDO

      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE CalculatePlasma(ir,count,s28te_,s28ti_,
     .   machsuper ,imaginary ,
     .   machsuper1,imaginary1,
     .   machsuper2,imaginary2)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_cedge2d
      use mod_slcom
      use mod_slcom_sol28
      IMPLICIT none

      INTEGER ir,count,s28te_,s28ti_,s28superdet_
      LOGICAL machsuper,imaginary,machsuper1,imaginary1,
     .                            machsuper2,imaginary2

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'cedge2d'
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'

      REAL GetMach,GetCs

      INTEGER ik,iksuper1,iksuper2,machmode,ikflat1,ikflat2
      LOGICAL cont,status,check,localimag2,gosuper2
      REAL   ! x0,x1,x2,x3,x4,x5,x6,y0,y1,y2,y3,y4,y5,y6,y2p,yA,yB,xB,
     .       deltat,frac, ! aA,aB,bA,bB,aC,bC,deltat,frac,L,Psol,k,f,
     .       a1,b1,c1,r1,sval, ! a1,b1,c1,r1,sB,sP,tB,tP,sval,texpt2,
     .       flx1x,ion1x,rec1x,cfp1x,mom1x,
     .       teprofile(MAXNKS),tres,tequ,cs,holdknbs(MAXNKS)
       REAL mom56 ! TEMP

      SAVE

      IF (new) THEN
c        CALL RZero(ktebs(1,ir),MAXNKS)
        CALL RZero(ktibs(1,ir),MAXNKS)
c        CALL RZero(knbs(1,ir) ,MAXNKS)
c        CALL RZero(kvhs(1,ir) ,MAXNKS)
        CALL RZero(teprofile  ,MAXNKS)
 



        machmode = 1

        ikflat1 = 0
        ikflat2 = MAXNKS + 1
        status = .TRUE.
        DO WHILE (status) 
          status = .FALSE.


          IF (s28te_.EQ.20.OR.s28te_.EQ.21) THEN
            ktebs(1:nks(ir),ir) = e2dtebs(1:nks(ir),ir) 
          ELSE
            CALL CalcTeProfile(ir,ikflat1,ikflat2,machmode,count)
          ENDIF



          IF     (s28ti_.EQ.20) THEN
            ktibs(1:nks(ir),ir) = e2dtibs(1:nks(ir),ir) 

          ELSEIF (s28ti_.EQ.0) THEN
c...
            DO ik = 1, nks(ir)
              ktibs(ik,ir) = ktebs(ik,ir)
            ENDDO

          ELSEIF (s28ti_.EQ.1) THEN
c...
            DO ik = 1, nks(ir)
              ktibs(ik,ir) = ktebs(ik,ir) * s28tiratio
            ENDDO

          ELSEIF (s28ti_.EQ.2) THEN
c...
            DO ik = ik0, ik6
              ktibs(ik,ir) = ktebs(ik,ir)      
            ENDDO 
           

            DO ik = ik2+1, ik3
              sval = kss(ik,ir) - s(2)

              frac = sval / (s(3) - s(2))
              frac = frac ** 0.90
              frac = frac ** 2.0
              ktibs(ik,ir) = frac * ktebs(ik,ir) * ti(3) / te(3) +
     .                       (1.0 - frac) * ktebs(ik,ir) 

c              WRITE(0,'(A,2I6,3F10.4)') 
c     .          'TI:',ik3,ir,0.0,ktebs(ik,ir),ktibs(ik,ir)


            ENDDO

            DO ik = ik3+1, ik4-1
              sval = s(4) - kss(ik,ir) 

              frac = sval / (s(4) - s(3))
              frac = frac ** 0.75
              frac = frac ** 2.0

              ktibs(ik,ir) = frac * ktebs(ik,ir) * ti(3) / te(3) +
     .                       (1.0 - frac) * ktebs(ik,ir) 

c              WRITE(0,'(A,2I6,3F10.4)') 
c     .          'TI:',ik3,ir,0.0,ktebs(ik,ir),ktibs(ik,ir)
    

            ENDDO


c            DO ik = ik0, ik6
c              WRITE(0,'(A,2I6,3F10.4)') 
c     .          'TI:',ik3,ir,0.0,ktebs(ik,ir),ktibs(ik,ir)
c            ENDDO 

            IF (s28ti_.EQ.3) THEN
              DO ik = 1, nks(ir)
                ktibs(ik,ir) = MAX(ktibs(ik,ir),2.0*ktebs(ik,ir))
              ENDDO 
            ENDIF

          ELSEIF (s28ti_.EQ.3) THEN
c...
            DO ik = ik0, ik6
              ktibs(ik,ir) = ktebs(ik,ir)      
            ENDDO 
           

            DO ik = ik0, ik3
              sval = kss(ik,ir) - s(2)

              frac = sval / (s(3) - s(2))
              frac = frac ** 0.4
              ktibs(ik,ir) = frac * ktebs(ik,ir) * ti(3) / te(3) +
     .                       (1.0 - frac) * ti(0) 

c              WRITE(0,'(A,2I6,3F10.4)') 
c     .          'TI:',ik3,ir,0.0,ktebs(ik,ir),ktibs(ik,ir)


            ENDDO

            DO ik = ik3, ik6
              sval = s(4) - kss(ik,ir) 

              frac = sval / (s(4) - s(3))
              frac = frac ** 0.4

              ktibs(ik,ir) = frac * ktebs(ik,ir) * ti(3) / te(3) +
     .                       (1.0 - frac) * ti(6)

c              WRITE(0,'(A,2I6,3F10.4)') 
c     .          'TI:',ik3,ir,0.0,ktebs(ik,ir),ktibs(ik,ir)
    

            ENDDO




          ELSEIF (.FALSE..AND.s28ti_.EQ.3.AND.count.GT.1) THEN

c...Crude model:
c     -doesn't account for poloidal distribution of ions out of the core
c     -doesn't account for thermalization as particles move cross-field
c

c I DON'T LIKE THE SHARP DROP BETWEEN IK3 and IK3+1...THINK ABOUT HOW ELSE TO DO THIS
            ktibs(ik3,ir) = ti(3)

c              WRITE(0,'(A,2I6,3F10.4)') 
c     .          'TI:',ik3,ir,0.0,ktebs(ik3,ir),ktibs(ik3,ir)


            DO ik = ik3+1, ik6

              cs = GetCs(ktebs(ik,ir),ktebs(ik,ir))
              deltas = 0.5*(ksb(ik,ir) - ksb(ik-1,ir))

              tres = deltas / cs

              tequ = ktebs(ik,ir)*1.5 / (1.6E-15 * holdknbs(ik) * 18.0)

              frac = MIN(1.0,tres / tequ)
 
              deltat = ktibs(ik+1,ir) - ktebs(ik+1,ir)

              ktibs(ik,ir) = ktebs(ik,ir) + (1.0 - frac) * deltat

c...          Ti can end up increasing slighling in the above scheme, so
c             make sure it doesn't happen:
              ktibs(ik,ir) = MAX(ktebs(ik,ir),
     .                           MIN(ktibs(ik,ir),0.995*ktibs(ik-1,ir)))
            
c              WRITE(0,'(A,2I6,6F10.4)') 
c     .          'TI:',ik,ir,frac,ktebs(ik,ir),ktibs(ik,ir),
c     .          (1.0 - frac) * deltat,deltat,ktibs(ik-1,ir)

c              WRITE(0,'(A,2I6,1P,3E10.2,3X,3E10.2,0P)') 
c     .          'TI:',ik,ir,cs,deltas,tres,
c     .          ktebs(ik,ir),holdknbs(ik),tequ 
            ENDDO
          ELSE
            CALL ER('CalculatePlasma','Unrecognized S28TI option',*99)
          ENDIF



          IF (.FALSE..AND.detach1) THEN
            STOP 'DETACH1 CODE UNTESTED'
c...        Smoothing of the temperature profile at the slope discontinuity:
            deltas = 0.5 * MIN(s(2)-s(1),s(3)-s(2))
            deltat = 0.10 * (ktebs(ik3,ir) - ktebs(ik2,ir))
            DO ik = 1, nks(ir)
              IF (kss(ik,ir).GE.kss(ik2,ir)-deltas.AND.
     .            kss(ik,ir).LE.kss(ik2,ir)+deltas) THEN
                frac = 1.0 - MIN(1.0,ABS(kss(ik,ir)-kss(ik2,ir))/deltas)
                frac = frac**2
                ktebs(ik,ir) = ktebs(ik,ir) + frac * deltat
                ktibs(ik,ir) = ktebs(ik,ir)
              ENDIF
            ENDDO

            IF (s(4).NE.ksmaxs(ir)) THEN
              deltas = 0.5 * MIN((s(5)-s(4)),s(4)-s(3))
              deltat = 0.10 * (ktebs(ik3,ir) - kteds(idds(ir,1)))
              DO ik = 1, nks(ir)
                IF (kss(ik,ir).GE.s(4)-deltas.AND.
     .              kss(ik,ir).LE.s(4)+deltas) THEN

                  frac = 1.0 - MIN(1.0,ABS(kss(ik,ir)-s(4))/deltas)
                  frac = frac**2

                  ktebs(ik,ir) = ktebs(ik,ir) + frac * deltat
                  ktibs(ik,ir) = ktebs(ik,ir)
                ENDIF
              ENDDO
            ENDIF
          ENDIF

          cont = .TRUE.
          check = .TRUE.
          iksuper1 = 0
          iksuper2 = MAXNKS+1
          gosuper2 = .FALSE.
          imaginary1 = .FALSE.
          imaginary2 = .FALSE.

c              CALL CalcIntegral3(osmmp,ik5  ,ik6,ir,mom56,4)
c              WRITE(0,*) 'MOM2:',mom56/ECH

          DO WHILE (cont)
            cont = .FALSE.
c...        Density and velocity:
c            WRITE(0,*) 'OSMMP:',osmmp(141,ir)
            DO ik = 1, nks(ir)
              CALL CalcIntegral3(osmion,1,ik,ir,ion1x,8)
              CALL CalcIntegral3(osmrec,1,ik,ir,rec1x,8)
              CALL CalcIntegral3(osmcfp,1,ik,ir,cfp1x,8)
              flx1x = ion1x + flx0 - rec1x + cfp1x

c... BUG??? -SL, NOV 18, 2004
              CALL CalcIntegral3(osmmp,1,ik,ir,mom1x,8)

c              CALL CalcIntegral3(osmmp,1,ik,ir,mom1x,2)

c              IF (ik.EQ.1) THEN
c                WRITE(0,'(A,I6,2E10.2,1P,2E10.2,0P)') 
c     .            'DATA:',ik,ktebs(ik,ir),ktibs(ik,ir),flx1x,mom1x
c              ENDIF

              localimag2 = .FALSE.
          
              a1 = (ktebs(ik,ir) + ktibs(ik,ir)) * ECH
              b1 = -(p0 + mom1x) 

c...WORK IN PROGRESS
c              IF (ir.EQ.18.AND.ik.GT.128.AND.ik.LE.ik4+5) b1 = b1 * 2.0

              c1 = AMU * crmb * flx1x**2
              r1 = b1**2 - 4.0 * a1 * c1

              IF     (s28te_.EQ.20) THEN
                knbs(ik,ir) = e2dnbs(ik,ir)

              ELSEIF (stopopt2.EQ.30.AND.irsep.EQ.10.AND.                ! **** HACK ****  
     .            (ir.GE.10.AND.ir.LE.17).AND.
     .            ik.GE.ik5) THEN
  
                IF (n(5).NE.0.0) THEN
                  frac = (kss(ik ,ir) - kss(ik5,ir)) /
     .                   (kss(ik6,ir) - kss(ik5,ir))
                  knbs(ik,ir) = (1.0 - frac) * n(5) + frac * n(6)
                ELSE
                  knbs(ik,ir) = n(6)
                ENDIF
c                knbs(ik,ir) = n(6)

              ELSEIF (r1.GE.0.0) THEN
c...            Real solution:
                IF ((ik.LE.ik3.AND.ik.LT.iksuper1).OR.
     .              (ik.GT.ik3.AND.ik.GT.iksuper2.AND.gosuper2)) THEN
                  knbs(ik,ir) = (-b1 - SQRT(r1)) / (2.0 * a1)
                ELSE
                  knbs(ik,ir) = (-b1 + SQRT(r1)) / (2.0 * a1)
                ENDIF
              ELSE

c *MOVE*
                s28superdet_ = 99
                IF (ringtype(ir).EQ.SOL1) s28superdet_ = s28superdet
                IF (ringtype(ir).EQ.PFZ ) s28superdet_ = s28superdetpfz
c
c               S28SUPERDET  0 - none
c                            1 - scale recombination source
c                            2 - allow target/plasma to go supersonic
c                            3 - 
c                            4 - increase cf flow into ring near target to reduce upstream flow 
c

c...            Solution of quadratic is imaginary:
                IF (ik.LE.ik3) THEN

                  IF (detach1.AND.
     .                (s28superdet_.EQ.0.OR.s28superdet_.EQ.3)) THEN
                  ELSEIF (.NOT.detach1.OR.
     .                    (s28superdet_.EQ.1.AND.ik.GE.ik1).OR.
     .                    (s28superdet_.EQ.2.AND.ik.LE.ik2+10).OR.
c...Only applies near the target:
     .                    (s28superdet_.EQ.4.AND.ik.GE.ik1.AND.
     .                     kss(ik,ir).LT.0.35*ksmaxs(ir).AND.
     .                     flx1x.LE.0.0)
     .                   ) THEN
                    machsuper1 = .TRUE.
                    imaginary1 = .TRUE.
                  ENDIF

                  IF (.NOT.detach1) THEN
                    IF     (machmode.EQ.3) THEN

                    ELSEIF (.NOT.detach1.AND.machmode.EQ.2.AND.
     .                      check.AND.ik.GT.ikflat1) THEN
c *DEFUNCT*
                      STOP 'THIS CODE IS DEFUNCT'
c                      x2 = kss(ik,ir)
                      ikflat1 = ik
                      check = .FALSE.
                    ENDIF
                  ENDIF

                ELSE

                  IF (detach2.AND.
     .                (s28superdet_.EQ.0.OR.s28superdet_.EQ.3)) THEN
c                  ELSEIF ((.FALSE..AND..NOT.detach2).OR.
                  ELSEIF (.NOT.detach2.OR.
     .                    (s28superdet_.EQ.1.AND.ik.LE.ik5).OR.
     .                    (s28superdet_.EQ.2.AND.ik.GE.ik4-10).OR.
c...Only applies near the target:
c     .                    (s28superdet_.EQ.4.AND.ik.LE.ik6.AND.
     .                    (s28superdet_.EQ.4.AND.ik.LE.ik5.AND.
     .                     kss(ik,ir).GT.0.65*ksmaxs(ir).AND.
     .                     flx1x.GE.0.0)
     .                   ) THEN
                    machsuper2 = .TRUE.
                    imaginary2 = .TRUE.
                    localimag2 = .TRUE.
                  ENDIF

                  IF (.NOT.detach2) THEN
                    IF     (machmode.EQ.3) THEN
                    ELSEIF (.NOT.detach2.AND.machmode.EQ.2) THEN
c *DEFUNCT*
                      STOP 'THIS CODE IS ALSO DEFUNCT'
c                      x4 = kss(ik,ir)
                      ikflat2 = ik
                    ENDIF
                  ENDIF
                ENDIF

c...            An approximation to KNBS:
                knbs(ik,ir) = -b1 / (2.0 * a1)

                IF (.FALSE..AND.s28output.AND.ik.LT.ik3) 
     .            WRITE(0,'(A,5I6,1P,4E10.2,0P,L2)') 
     .              'IMAGINARY:',ik,ik4,ik5,ik6,ir,b1,a1,c1,flx1x,
     .              localimag2

                IF (.FALSE..AND.s28output.AND.ik.GT.ik3-3) 
     .            WRITE(PINOUT,'(2F10.2,1P,E10.2,0P,5I6,A)') 
     .              ktebs(ik,ir)+ktibs(ik,ir),p0+mom1x,flx1x,
     .              ik,ir,ik4,ik5,ik6,'IMAGINARY'


              ENDIF

c...          Calculate the plasma velocity:            
              IF     (s28te_.EQ.20) THEN
                kvhs(ik,ir) = e2dvhs(ik,ir)
              ELSE
                kvhs(ik,ir) = flx1x / knbs(ik,ir)
              ENDIF

              mach1 = GetMach(kvhs(ik,ir),ktebs(ik,ir),ktibs(ik,ir))
              mach2 = GetMach(kvhs(ik,ir),ktebs(ik,ir),ktibs(ik,ir))

              IF     (machmode.EQ.3) THEN

              ELSEIF (machmode.EQ.2) THEN

                IF (.NOT.detach1.AND.imaginary1.OR.
     .              .NOT.detach2.AND.imaginary2) status = .TRUE.

              ELSE

c...            Decide when to apply supersonic density value and then recalculate
c               plasma densities:
                IF     ((.NOT.detach1.OR.s28superdet_.EQ.2).AND.
     .                  ik.LE.ik3) THEN
                  IF (machsuper1.AND.mach1.GT.0.98.AND.
     .                iksuper1.EQ.0) THEN
                    iksuper1 = ik
                    cont = .TRUE.
                  ENDIF     
                ELSEIF ((.NOT.detach2.OR.s28superdet_.EQ.2).AND.
     .                  ik.GT.ik3) THEN
                  IF (machsuper2.AND.mach2.GT.0.98.AND.
c                  IF (.FALSE..AND.machsuper2.AND.mach2.GT.0.98.AND.
     .                iksuper2.EQ.MAXNKS+1) THEN
                    iksuper2 = ik
                    cont = .TRUE.
                  ENDIF
                ENDIF     
              ENDIF

              IF (s28output.AND.ik.LE.ik3)
c              IF (machsuper1.AND.ik.LE.15) 
     .             WRITE(0,'(A,3I6,L2,2F10.2,1P,2E10.2,0P,
     .                    F10.4,1P,4E14.6,0P)') 
     .             'N:',ik,ikflat1,
     .             iksuper1,imaginary1,
     .             ktebs(ik,ir),ktibs(ik,ir),knbs(ik,ir),kvhs(ik,ir),
     .             mach1,
     .             osmion(ik,ir),osmrec(ik,ir),osmcfp(ik,ir),
     .             osmmp(ik,ir)


c              IF (s28output.AND.ik.GE.ik3)
c              IF (machsuper2.AND.ik.GE.50) 
                IF (.FALSE..AND.ik.GT.ik4)
     .           WRITE(0,'(A,3I6,L2,2F10.2,1P,2E10.2,0P,
     .                     F10.4,1P,4E14.6,0P)') 
     .             'N:',ik,ikflat2,
     .             iksuper2,imaginary2,
     .             ktebs(ik,ir),ktibs(ik,ir),knbs(ik,ir),kvhs(ik,ir),
     .             mach2,
     .             osmion(ik,ir),osmrec(ik,ir),osmcfp(ik,ir),
     .             osmmp(ik,ir)



            ENDDO
c            WRITE(0,*) 'OSMMP:',osmmp(141,ir)

            gosuper2 = .TRUE.

c...        End of DO WHILE density,velocity loop:
          ENDDO

          IF (.NOT.status.AND.machmode.EQ.2) THEN 
            STOP 'OLD CODE'
            CALL MachProfile(ir,teprofile)
            status = .TRUE.
            machmode = 3

          ENDIF



c...      End of DO WHILE temperature loop:
        ENDDO

        IF (s28output) THEN
          CALL CalcIntegral3(osmion,1,nks(ir),ir,ion1x,8)
          CALL CalcIntegral3(osmrec,1,nks(ir),ir,rec1x,8)
          CALL CalcIntegral3(osmcfp,1,nks(ir),ir,cfp1x,8)
          WRITE(0,*) 'CHECK:',-flx0,flx6,-flx0+flx6   
          WRITE(0,*) 'CHECK:',-v(0)*n(0),v(6)*n(6),
     .                        -v(0)*n(0)+v(6)*n(6)
          WRITE(0,*) 'CHECK:',ion1x,rec1x,cfp1x,ion1x-rec1x+cfp1x   
  
          CALL CalcIntegral3(osmion,1,nks(ir),ir,ion1x,4)
          CALL CalcIntegral3(osmrec,1,nks(ir),ir,rec1x,4)
          CALL CalcIntegral3(osmcfp,1,nks(ir),ir,cfp1x,4)
          WRITE(0,*) 'CHECK:',ion1x,rec1x,cfp1x,ion1x-rec1x+cfp1x   
        ENDIF

        IF (machmode.EQ.2.OR.machmode.EQ.3) THEN
c...      Make sure target Mach number iterations are not done:
          machsuper1 = .FALSE.
          machsuper2 = .FALSE.
c          STOP 'sdfsdg'
        ENDIF


c...    Store KNBS values:
        DO ik = 1, nks(ir)
          holdknbs(ik) = knbs(ik,ir)
        ENDDO

c TMP          
c        machsuper2 = .FALSE.
c        imaginary2 = .FALSE.

c
c * RETURN *
c
        RETURN
      ENDIF



      RETURN
 99   STOP
      END
c
c
c
c
c
      SUBROUTINE SOL28(irs,ire,ikopt)
      USE mod_sol28_old
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_cedge2d
      use mod_slcom
      use mod_slcom_sol28
      use mod_sl_oldplasma
      IMPLICIT none

      INTEGER irs,ire,ikopt

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'cedge2d'
c     INCLUDE 'slcom'
c     INCLUDE 'slcom_sol28'

      REAL Clock2,GetMach


!      COMMON /OLDPLASMA/ oldknbs ,oldktebs ,oldktibs ,oldkvhs ,
!     .                   oldknbs2,oldktebs2,oldktibs2,oldkvhs2
!      REAL
!     .     oldknbs  (MAXNKS,MAXNRS),oldktebs (MAXNKS,MAXNRS),
!     .     oldktibs (MAXNKS,MAXNRS),oldkvhs  (MAXNKS,MAXNRS),
!     .     oldktebs2(MAXNKS,MAXNRS),oldktibs2(MAXNKS,MAXNRS),
!     .     oldknbs2 (MAXNKS,MAXNRS),oldkvhs2 (MAXNKS,MAXNRS)

      COMMON /TIMESTATS/ timeint
      REAL               timeint(100)


      REAL CalcPressure

      INTEGER index,type
      REAL    pressure

      INTEGER ik,ir,ir1,nd,i1,machcount,count,fp,lastik,s28te_,s28ti_,
     .        s28superdet_
      LOGICAL machsuper,iterate,imaginary,status,stoprun,
     .        imaginaryL,tempset,s28done(MAXNRS)
 
c...NEW
      INTEGER machcount1,
     .        machcount2
      LOGICAL machsuper1,imaginary1,lockrec1,
     .        machsuper2,imaginary2,lockrec2
      REAL    machadjust1,
     .        machadjust2

      REAL, ALLOCATABLE :: linedat(:,:)

      REAL    tarden,lastadjust1,lastadjust2,tmprec(MAXNKS),mi,
c     .        x1,x2,x3,x5,
     .        deltat,sval,lastn1,deltap,
     .        machval,frac,
     .        machadjust,nf,
     .        tempte,mulden,multe,tmpflx1,tmpflx2
      REAL*8  integral
      CHARACTER*10 note

      REAL csmax,getcs,lastcsmax,csmax1,csmax2

      INTEGER tmpcount,tmpcount2

      LOGICAL firstcall,firstcall2
      DATA    firstcall,firstcall2 /.TRUE.,.FALSE./

      SAVE

c      WRITE(0,*) 'SOL28 RINGS:',irs,ire

      IF (s28mode.GE.4.0) THEN 
c        IF (sloutput) WRITE(0,*) 'WHAT-WHO! CALLING SOL28_V4'
        CALL ExecuteSOL28(MAX(2,irs),ire,ikopt,sloutput)
        RETURN
      ENDIF


      s28output = .NOT..TRUE.
      new       = .FALSE.

      IF (s28mode.EQ.1.0) THEN
         IF (sloutput) WRITE(0,*) 'CALLING SOL28 CLASSIC'
         CALL SOL28_OLD(irs,ire,ikopt)
         RETURN
      ELSEIF (s28mode.EQ.1.1) THEN
         IF (sloutput) WRITE(0,*) 'TURNING ON SOL28 DGAMMA FITTING'      
        s28fit = 1
      ELSEIF (s28mode.GE.2.0) THEN
         IF (sloutput) WRITE(0,*) 'ACTIVATING NEW SOL28 ROUTINES'
         new = .TRUE.
      ENDIF

      IF (s28mode.LT.2.0.AND.eiropacity.EQ.5) THEN
        WRITE(0,*)
        WRITE(0,*) '---------------------------------------'
        WRITE(0,*) 'NOTE FULL LYMAN TRAPPING WITH SOL28_OLD'
        WRITE(0,*) '---------------------------------------'
        WRITE(0,*) 
c        WRITE(0,*) 'ERROR: INITRECOM NOT WORKING PROPERLY WITH'//
c     .             'EIROPACITY.EQ.5'
c        STOP 'HALT'
      ENDIF

c...THERE IS A HUGE PROBLEM ITERATING THIS ROUTINE
c      -RECOMBINATION RATES ARE NOT SCALED BY FACTOR 3 TO START, AS
c       THEY ARE IN EIRENE
c      -EIRENE SOURCES MAKE NO SENSE AS YET! 
c      -THE DGAMMA FIT MAY CRAP OUT IF IT IS APPLIED ONLY AFTER
c       THE SOLUTION IS FOUND, SINCE ON SUBSEQUENT ITERATIONS
c       THE SOLUTION MAY NOT BE FINDABLE WITHOUT THE SCALING ALREADY IN PLACE ... PERHAPS
c       SOME RELAXATION WILL HELP THIS
c       


      IF (.NOT.new.AND.firstcall2) THEN
c...    Reload old solution and go home:
        WRITE(0,*)
        WRITE(0,*) '*** BACKLOADING SOL28 SOLUTION ***'
        WRITE(0,*)
        WRITE(PINOUT,*)
        WRITE(PINOUT,*) '*** BACKLOADING SOL28 SOLUTION ***',rel_nstep
        WRITE(PINOUT,*)
        DO ir = irtrap+1, nrs
          DO ik = 1, nks(ir)
            ktebs(ik,ir) = oldktebs(ik,ir)
            ktibs(ik,ir) = oldktibs(ik,ir)
            knbs(ik,ir) = oldknbs(ik,ir)
            kvhs(ik,ir) = oldkvhs(ik,ir)
          ENDDO
        ENDDO
        RETURN
      ENDIF

      firstcall2 = .TRUE.


      CALL SetBounds

      tmpcount = 0
      tmpcount2 = 0


      mulden = 0.0
      multe = 0.0
c      mulden = 0.2
c      multe = 0.10
 
      frac = MAX(0.0,REAL(rel_step - 1) / REAL(rel_nstep - 1))
      frac = rel_bound1 + frac * (rel_bound2 - rel_bound1) 
c      WRITE(0,*) 'SOL28 FRAC:',frac
      frac = 0.0

      WRITE(PINOUT,*) 'SOL28: ITER',rel_iter,rel_niter,frac

      IF (relreset.GE.1.AND.rel_iter.EQ.rel_niter) firstcall = .TRUE.
 
      IF (firstcall) THEN
        IF (osm_store.GT.-1) THEN
          loadsources = .TRUE.
        ELSE
          loadsources = .FALSE.
        ENDIF
c        WRITE(0,*) 'SETTING MULRECL TO 1.0'
        CALL RSet(mulrecl,MAXNRS,1.0)
      ELSE
        loadsources = .TRUE.
        CALL RSet(mulrecl,MAXNRS,1.0)
      ENDIF

      IF (.NOT.new) THEN
        WRITE(0,*)
        WRITE(0,*) '*************************'
        WRITE(0,*) '* LOADSOURCES = .FALSE. *'
        WRITE(0,*) '*************************'
        WRITE(0,*)
        loadsources = .FALSE.
      ENDIF



      tempset = .FALSE.
      scalerec = .FALSE.
c      IF (new) scalerec = .TRUE.

      IF (sloutput) WRITE(0,*) 'SOL28:',loadsources,scalerec
      WRITE(PINOUT,*) 'SOL28:',loadsources,s28ion,s28mom
c
c     =================================
c                 MAIN LOOP
c     =================================
c
      fitline = .FALSE.


      IF (s28fit.EQ.1) THEN
c NEEDS TO BE FANCIER
        CALL LoadEIRENEAtomicData
        ALLOCATE(linedat(MAXNKS,MAXNRS))
        fp = 98
        OPEN(UNIT=fp,FILE='match.dat',ACCESS='SEQUENTIAL',STATUS='OLD')
       READ(fp,*)
        DO WHILE (.TRUE.)
          READ(fp,*,END=10,ERR=10) ik,ir,linedat(ik,ir)
c          WRITE(0,*) ik,ir,linedat(ik,ir)
        ENDDO
 10     CLOSE(fp)
      ENDIF

c *CRUDE* Should be elsewhere... can IDRING take the place of RINGTYPE for EIRENE04? 
      ringtype = 0
      DO ir = 1, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        IF (ir.LT.irsep) ringtype(ir) = CORE
        IF (ir.GE.irsep.AND.ir.LT.irwall) ringtype(ir) = SOL1
        IF (ir.GT.irtrap.AND.ir.LE.nrs) ringtype(ir) = PFZ
      ENDDO
c...  Special for secondary PFZ in generalized grids:
      IF (grdnmod.NE.0) THEN
        DO ik = 1, nks(irwall)
          IF (irins(ikins(ik,irwall),irins(ik,irwall)).EQ.irwall.AND.
     .        ringtype(irins(ik,irwall)).NE.PFZ) THEN
c...        Here we go, trying to use the connection map to identify 
c           the rather poorly defined secondary PFZ.  
            DO ir1 = irins(ik,irwall), nrs
c              WRITE(0,*) 'RING1:',ir1,irouts(1,ir1),irouts(nks(ir1),ir1)
              IF     (irouts(1,ir1).EQ.irouts(nks(ir1),ir1))THEN
                ringtype(ir1) = PFZ
              ELSEIF (irouts(1       ,ir1).EQ.irwall.OR.
     .                irouts(nks(ir1),ir1).EQ.irwall) THEN
c...            Discontinuous target!
                STOP 'THIS IS A PROBLEM'
              ELSE
                ringtype(ir1) = PFZ
                EXIT
              ENDIF
c              WRITE(0,*) '     :',ir1,irouts(1,ir1),irouts(nks(ir1),ir1)
            ENDDO
c            WRITE(0,*) 'PROBLEM FOUND!',irins(ik,irwall)
          ENDIF
        ENDDO
      ENDIF

c      DO ir1 = 2, nrs
c        WRITE(0,*) 'RINT:',ir1,ringtype(ir1)
c      ENDDO

      count = 0
      ir = irs
      lockrec1 = .FALSE.
      lockrec2 = .FALSE.
      tmprec = 0.0
      machsuper1 = .FALSE.
      machsuper2 = .FALSE.
      lastadjust1 = 0.9
      lastadjust2 = 0.9
      machadjust1 = 1.0
      machadjust2 = 1.0
      DO WHILE (ir.GE.irs.AND.ir.LE.ire)

        s28done(ir) = .TRUE.

        IF (idring(ir).EQ.BOUNDARY) THEN
          ir = ir + 1
          CYCLE
        ENDIF

c...    SPECIAL:
        IF (nrs.EQ.45) THEN
          IF (ir.GE.20.AND.ir.LE.31) THEN
            IF (sloutput) THEN
              WRITE(0,*)
              WRITE(0,*) '------------------------------------------'
              WRITE(0,*) 'SUPER CHEAT ON RINGTYPE FOR MAST UPPER PFZ'
              WRITE(0,*) '------------------------------------------'
              WRITE(0,*)
            ENDIF
            ringtype(ir) = PFZ
          ENDIF
        ENDIF
 

c...    Initialization:
        s = 0.0

c...    Dig out peak parameters:
        IF (s28mode.GE.3.0) THEN
          te = 0.0
          n = 0.0
          CALL FindS28Parameters_V3(ir,te,n,nf,s,new)
        ELSE
          CALL FindS28Parameters   (ir,te,n,nf,s,new)
        ENDIF
c        WRITE(0,*) 'TE(0):',te(0),te(5),te(6),n(1)

c...    CORE!
        IF (ringtype(ir).EQ.CORE) THEN
          WRITE(0,*) 'SOL28_SOL: CORE RING',ir,te(3),n(3)
          ti(3) = te(3)
          DO ik = 1, nks(ir)
            ktebs(ik,ir) = te(3)
            ktibs(ik,ir) = ti(3)
            knbs (ik,ir) = n (3)
            kvhs (ik,ir) = 0.0
          ENDDO
          ir = ir + 1
          CYCLE
        ENDIF

c...    Simple global scaling of density peak:
        CALL AdjustDensity(n(1))
c        WRITE(0,*) 'TE(0):',te(0),te(5),te(6),n(1)


c...    SOLPS:



c...    Check for over-rides of the SOL28 parameters:
        DO i1 = 1, osmns28over
          IF (NINT(osms28over(i1,1)).NE.ir) CYCLE
          index = NINT(osms28over(i1,2))
          type  = NINT(osms28over(i1,3))
          IF     (type.EQ.1 ) THEN
c...        P1 - 
c           P2 -
c           P3 -
            IF (osms28over(i1,4).NE.99.0) te(index) = osms28over(i1,4)
            IF (osms28over(i1,5).NE.99.0) n (index) = osms28over(i1,5)

          ELSEIF (type.EQ.2 ) THEN
c...        P1 - 
c           P2 -
c           P3 -
            pressure = 0.0
            IF (osms28over(i1,4).NE.99.0) te(index) = osms28over(i1,4)
            IF (osms28over(i1,5).NE.99.0) pressure  = osms28over(i1,5)
            n(index) = pressure / (2.0 * ABS(te(index))) / ECH

          ELSE
            CALL ER('SOL28','Unknown over-ride TYPE',*99)
          ENDIF 
        ENDDO


c...
        CALL InitializeRingSOL28(ir,count)
c        WRITE(0,*) 'TE(0):',te(0),te(6),te(6),n(1)

        IF (new) THEN
          detach1 = .TRUE.
          detach2 = .TRUE.
          IF (s(1).EQ.s(0)) detach1 = .FALSE.
          IF (s(5).EQ.s(6)) detach2 = .FALSE.
        ENDIF


        IF (new) THEN

c            WRITE(0,*) 'EH?',te(0),detach1

          IF (detach1.AND.te(1).LT.0.0) THEN
c            WRITE(0,*) 'WORKING:',te(1)
            te(0) = -te(1)
            te(1) = -te(1)
          ELSE

            IF     (te(0).EQ.-98.0) THEN
              te(0) = te(3)
            ELSEIF (te(0).EQ.97.0.OR.te(0).EQ.-97.0) THEN
              IF (sloutput) THEN
                WRITE(0,*) '*******************************************'
                WRITE(0,*) '*WARNING: NOT CONSISTENT WITH S28TI OPTION*'
                WRITE(0,*) '*******************************************'
              ENDIF
c              STOP 'SOL28: CHECK THAT UPSTREAM PRESSURE CORRECT'    
              p3 = n(3) * 2.0 * te(3) * ECH
              mi = crmb * AMU
              te(0) = (p3 / isat(0))**2 * ECH / (8.0 * mi)   
c              WRITE(0,*) 'TE=',p3,mi,isat(0)
c              WRITE(0,*) 'TE=',te(0),te(3)
            ELSEIF (te(0).LT.0.0) THEN
              te(0) = -te(0)
            ELSE
              te(0) = kteds(idds(ir,2))
            ENDIF
          ENDIF

          IF (detach2.AND.te(5).LT.0.0) THEN
            te(6) = -te(5)
            te(5) = -te(5)           
          ELSE
            IF     (te(6).EQ.-98.0) THEN
              te(6) = te(3)
            ELSEIF (te(6).EQ.97.0.OR.te(6).EQ.-97.0) THEN
              IF (sloutput) THEN
                WRITE(0,*) '*******************************************'
                WRITE(0,*) '*WARNING: NOT CONSISTENT WITH S28TI OPTION*'
                WRITE(0,*) '*******************************************'
              ENDIF
c              STOP 'SOL28: CHECK THAT UPSTREAM PRESSURE CORRECT'    
              p3 = n(3) * 2.0 * te(3) * ECH
              mi = crmb * AMU
              te(6) = (p3 / isat(6))**2 * ECH / (8.0 * mi)   
c              WRITE(0,*) 'TE=',p3/ECH,mi,isat(6)
c              WRITE(0,*) 'TE=',te(6),te(3)    
            ELSEIF (te(6).LT.0.0) THEN
              te(6) = -te(6)
            ELSE
              te(6) = kteds(idds(ir,1))
            ENDIF
          ENDIF

          IF (detach1) THEN
          ELSE
            te(1) = te(0)
            te(2) = te(0)
          ENDIF
          IF (detach2) THEN
          ELSE
            te(4) = te(6)
            te(5) = te(6)
          ENDIF


          IF (.NOT..TRUE..AND.irsep.EQ.10.AND.                ! **** HACK **** OFF?
     .        (ir.GE.10.AND.ir.LE.17)) THEN
 
            techeat = te(6)
            te(6) = 0.25
            ti(6) = te(6)

            WRITE(0,*) 'TECHEAT ASSIGNED:',techeat

          ENDIF



c...      *HACK* Need a better way to trigger this!
          IF (s28ti.EQ.0) THEN
c            WRITE(0,*) 'NE:',ir,n(0),n(6)
c            WRITE(0,*) 'NE:',ir,n(3),te(3)
c            WRITE(0,*) 'ISATi,o:',ir,isat(0),isat(6)
            IF (n(0).EQ.-98.0) 
     .        isat(0)=0.5 * n(3)*GetCs(te(3),te(3))*ECH    ! LEFT OFF!? IS THE RP DENSITY FACT 2 TOO LOW?
            IF (n(6).EQ.-98.0) 
     .        isat(6)=0.5 * n(3)*GetCs(te(3),te(3))*ECH    ! CALCULATON WRONG?
c            IF (isat(0).LT.0.0) isat(0)=0.5*n(3)*GetCs(te(3),te(3))*ECH 
c            IF (isat(6).LT.0.0) isat(6)=0.5*n(3)*GetCs(te(3),te(3))*ECH
          ELSE
c...        Mess...
            ti(3) = te(3) * s28tiratio
            IF (n(0).EQ.-98.0) 
     .        isat(0)=0.5 * n(3)*GetCs(te(3),ti(3))*ECH 
            IF (n(6).EQ.-98.0) 
     .        isat(6)=0.5 * n(3)*GetCs(te(3),ti(3))*ECH
c            IF (isat(0).LT.0.0) isat(0)=0.5*n(3)*GetCs(te(3),ti(3))*ECH 
c            IF (isat(6).LT.0.0) isat(6)=0.5*n(3)*GetCs(te(3),ti(3))*ECH
          ENDIF
c            WRITE(0,*) 'ISATi,o:',ir,isat(0),isat(6)





c          WRITE(0,'(A,7F10.4)') 'S= ',(s(i1),i1=0,6)
c          WRITE(0,'(A,7F10.4)') 'te=',(te(i1),i1=0,6)
c          WRITE(0,'(A,1P,7E10.2)') 'n=',(n(i1),i1=0,6)
c          WRITE(0,'(A,7E10.2)') 'isat=',(isat(i1),i1=0,6)
c          WRITE(0,*) 'DETACH:',detach1,detach2
c          WRITE(0,'(A,7I10)')   'IKs',ik0,ik1,ik2,ik3,ik4,ik5,ik6

c          STOP 'sdfgsd'
        ELSE
        ENDIF

c *MOVE*
        IF (ringtype(ir).EQ.SOL1) s28te_ = s28te
        IF (ringtype(ir).EQ.PFZ ) s28te_ = s28tepfz

        IF (ringtype(ir).EQ.SOL1) s28ti_ = s28ti
        IF (ringtype(ir).EQ.PFZ ) s28ti_ = s28tipfz


        IF (s28te_.EQ.20.OR.s28te_.EQ.21) THEN
          te(0) = e2dtebs(1      ,ir) 
          te(6) = e2dtebs(nks(ir),ir) 
        ENDIF

        IF     (s28ti_.EQ.20) THEN
          ti(0) = e2dtibs(1      ,ir) 
          ti(6) = e2dtibs(nks(ir),ir) 
        ELSEIF (s28ti_.EQ.1) THEN
          IF (ir.GT.irtrap.AND.s28tiratio.NE.1.0) 
     .      CALL ER('SOL28_OLD','Option not ready',*99)
          ti(0) = te(0) * s28tiratio
          ti(1) = te(1) * s28tiratio
          ti(2) = te(2) * s28tiratio
          ti(3) = te(3) * s28tiratio
          ti(4) = te(4) * s28tiratio
          ti(5) = te(5) * s28tiratio
          ti(6) = te(6) * s28tiratio
c...      THIS ASSUMES THAT THE N(3) DATA ORIGINATED WITH THE ASSUMPTION THAT TE=TI:
          n(3) = n(3) * SQRT((2.0 * te(3)) / (te(3) + ti(3)))
        ELSEIF (s28ti_.EQ.2) THEN
          IF (ir.GT.irtrap.AND.s28tiratio.NE.1.0) 
     .      CALL ER('SOL28_OLD','Option not ready',*99)
          ti(0) = te(0)
          ti(1) = te(1)
          ti(2) = te(2)
          ti(3) = te(3) * s28tiratio
          ti(4) = te(4)
          ti(5) = te(5)
          ti(6) = te(6)
c...      THIS ASSUMES THAT THE N(3) DATA ORIGINATED WITH THE ASSUMPTION THAT TE=TI:
          n(3) = n(3) * SQRT((2.0 * te(3)) / (te(3) + ti(3)))
        ELSEIF (s28ti_.EQ.3) THEN
          ti(0) = ktids(idds(ir,2))
          ti(1) = ktids(idds(ir,2))
          ti(2) = ktids(idds(ir,2))
          ti(3) = te(3) 
          ti(4) = ktids(idds(ir,1))
          ti(5) = ktids(idds(ir,1))
          ti(6) = ktids(idds(ir,1))
        ELSE
c          te(3) = te(3) * 0.50
c          n(3) = n(3) * SQRT(1.0 / 0.50)

          ti(0) = te(0)
          ti(1) = te(1)
          ti(2) = te(2)
          ti(3) = te(3)
          ti(4) = te(4)
          ti(5) = te(5)
          ti(6) = te(6)
        ENDIF

        IF (s28te_.EQ.20.OR.s28te_.EQ.21) THEN
          isat(0) = e2dnbs(1      ,ir) * ECH * ABS(e2dvhs(1      ,ir))
          isat(6) = e2dnbs(nks(ir),ir) * ECH * ABS(e2dvhs(nks(ir),ir))

c          isat(0) = e2dnbs(1      ,ir) * ECH * GetCs(te(0),ti(0))
c          isat(6) = e2dnbs(nks(ir),ir) * ECH * GetCs(te(6),ti(6))          

          WRITE(0,*) 'WHAT',e2dnbs(1,ir),e2dnbs(nks(ir),ir)
          WRITE(0,*) 'WHAT',e2dion(1,ir),e2dion(nks(ir),ir)
        ENDIF

c...
        CALL ParticleConservation(ir,count,machadjust1,machadjust2)



        IF (new) THEN
          IF (te(0).EQ.98.0) THEN
            te(1) = te(0)
            te(2) = te(0)
          ENDIF
          IF (te(6).EQ.98.0) THEN
            te(1) = te(0)
            te(2) = te(0)
          ENDIF
        ELSE
        ENDIF

        IF (new) THEN
        ELSE
        ENDIF

c...NECESSARY?  JUST DO THIS LATER FOR THE WHOLE RING?



        IF (new) THEN
          v(0) = -GetCs(te(0),ti(0))
          v(6) =  GetCs(te(6),ti(6))
c...      Ill-defined flx0 here?
          n(0) = flx0 / v(0)
          n(6) = flx6 / v(6)
        ELSE
        ENDIF

        IF (new.AND.count.EQ.0.AND.sloutput) THEN
          WRITE(0,*) 'DETACH  :',ir,detach1,detach2
          WRITE(0,'(A,7(F11.4:))'   ) ' s =',(s (i1),i1=0,6)
          WRITE(0,'(A,7(I11:  ))'   ) ' ik=',ik0,ik1,ik2,ik3,ik4,ik5,ik6
          WRITE(0,'(A,7(F11.4:))'   ) ' te=',(te(i1),i1=0,6)
          WRITE(0,'(A,1P,7(E11.2:))') ' ne=',(n (i1),i1=0,6)
          WRITE(0,'(A,7(F11.4:))'   ) ' v =',(v (i1),i1=0,6)
        ENDIF

c          WRITE(0,*) 'DETACH :',ir,detach1,detach2
c          WRITE(0,'(A,7(F11.4:))'   ) ' s =',(s (i1),i1=0,6)
c          WRITE(0,'(A,7(F11.4:))'   ) ' te=',(te(i1),i1=0,6)
c          WRITE(0,'(A,1P,7(E11.2:))') ' ne=',(n (i1),i1=0,6)
c          WRITE(0,'(A,7(F11.4:))'   ) ' v =',(v (i1),i1=0,6)
 
c        STOP 'HELLO!'



        IF (new) THEN
c p0 and p6 claculation used to be here ... perhaps it had to be outside the mach number iteration?
        ELSE
c...BUG! CHECK IT OUT!  THE 2's SHOULD NOT BE THERE IF V IS ASSIGNED!
c CHANGE
          p0 = CalcPressure(n(0),te(0),ti(0),v(0))*ECH
c          p0 = CalcPressure(n(0),2.0*te(0),2.0*ti(0),v(0))*ECH
          p1 = CalcPressure(n(1),    te(1),    ti(1),v(1))*ECH
          IF (s(4).NE.ksmaxs(ir)) THEN
c...        For p4:
            IF (n(4).LT.0.0) THEN
              p5 = CalcPressure(holdn2,te(5),ti(5),holdv2)*ECH
              p4 = -n(4) * p5
            ELSE
              CALL CalcIntegral3(osmion,1,ik4,ir,ion04,2)
              CALL CalcIntegral3(osmrec,1,ik4,ir,rec04,2)
              CALL CalcIntegral3(osmcfp,1,ik4,ir,cfp04,2)
              flx4 = ion04 + flx0 - rec04 + cfp04
              v(4) = flx4 / n(4)
              p4 = CalcPressure(n(4),te(4),ti(4),v(4))*ECH
            ENDIF
          ENDIF
        ENDIF


        IF (new) THEN

c *MOVE*
          s28superdet_ = 99
          IF (ringtype(ir).EQ.SOL1) s28superdet_ = s28superdet
          IF (ringtype(ir).EQ.PFZ ) s28superdet_ = s28superdetpfz

          IF (s28superdet_.EQ.1.OR.s28superdet_.EQ.3.OR.
     .        s28superdet_.EQ.4) THEN
            scalerec = .TRUE.
          ELSE
            scalerec = .FALSE.
          ENDIF

          IF (.NOT.detach1.OR.s28superdet_.EQ.2) THEN
            machsuper1 = .FALSE.
            machadjust1 = 1.01
            machcount1 = 0
          ENDIF
          IF (.NOT.detach2.OR.s28superdet_.EQ.2) THEN
            machsuper2 = .FALSE.
            machadjust2 = 1.01
            machcount2 = 0
          ENDIF
        ELSE
          machsuper  = .FALSE.
          machadjust = 1.01
          machcount  = 0
        ENDIF

        iterate    = .TRUE.


c        WRITE(0,*) 'NEW:',ir,tmpcount

        DO WHILE (iterate)


          IF (new) THEN
c...        Done in Calc...
          ELSE
            imaginaryL  = .FALSE.
            imaginaryL2 = .FALSE.
            imaginary = .FALSE.
          ENDIF



          IF (new) THEN
c...        This okay here?  Perhaps destabilizing?
            p0 = CalcPressure(n(0),te(0),ti(0),v(0))*ECH
            p6 = CalcPressure(n(6),te(6),ti(6),v(6))*ECH
c            WRITE(0,*) 'P0,6:',p0,p6
          ENDIF


c
c
c
c
c
          CALL MomentumConservation(ir,count)



          CALL CalculatePlasma(ir,count,s28te_,s28ti_,
     .           machsuper ,imaginary ,
     .           machsuper1,imaginary1,
     .           machsuper2,imaginary2)


           IF (.NOT..TRUE.) THEN
             IF (count.EQ.0) WRITE(0,*) '**** SUPER CHEAT ***',ir,n(6)
             DO ik = ik5, ik6
               knbs(ik,ir) = n(6)
             ENDDO
           ENDIF


              mach1 = ABS(v(0)) / GetCs(te(0),ti(0))
              mach2 = ABS(v(6)) / GetCs(te(6),ti(6))


c          IF (new) THEN
c          IF (new.AND.count.EQ.10) THEN
          IF (.FALSE..AND.new) THEN
            knds(id0) = n(0)
            knds(id6) = n(6)
            kvds(id0) = v(0)
            kvds(id6) = v(6)
            DO ir = 1, nrs
              DO ik = 1, nks(ir)
                pinion(ik,ir) = osmion(ik,ir)
              ENDDO
            ENDDO              
            CALL DumpGrid('NEW STUFF IS FUN!')
          ENDIF

c          WRITE(0,*) 'WHAT THE:',imaginary,machsuper,imaginaryL2
c          STOP 'sdfgsd'

c...      Turn on target mach number iteration as necessary:
          IF (new) THEN

            iterate = .FALSE.

            IF ((.NOT.detach1.OR.s28superdet_.EQ.2).AND.
     .          machsuper1.AND.machcount1.LT.5) THEN
c     .          machsuper1.AND.machcount1.LE.6) THEN

              iterate = .TRUE.

c              WRITE(0,*) 'MACH INNER:',machcount1

              IF (imaginary1) THEN

                IF (machadjust1.LT.1.0) 
     .            machadjust1 = 1.0 + 0.5 * (1.0 - machadjust1)

              ELSEIF (lastp0.LT.p0.AND.machadjust1.LT.1.0) THEN
                WRITE(0,*) 'FUNK - something is wrong'
                machadjust1 = 1.0 + 0.5 * (1.0 - machadjust1)
                p6 = HI

              ELSEIF (ABS(1.0-machadjust1).LT.1.0E-5) THEN
c...            Seems limiting machcount to 6 is not good enough when the target
c               is only slightly supersonic:
                iterate = .FALSE.

              ELSEIF (.TRUE.) THEN

                IF (machadjust1.GT.1.0) THEN
                   machadjust1 = 1.0 - 0.5 * (machadjust1 - 1.0)
                   machcount1  = machcount1 + 1
                ELSEIF (machadjust1.LT.1.0.AND.machno.LT.0.5) THEN
c...              Turn this silly business off and set target conditions to
c                 give Mach number unity since something else has
c                 made the outer target no longer supersonic:
                  STOP 'DFDC2'
                ENDIF
              ENDIF
c...          Adjust target conditions to account for the change in Mach number
c             at the target:
              n(0) = n(0) / machadjust1
              v(0) = v(0) * machadjust1
            ENDIF

            IF ((.NOT.detach2.OR.s28superdet_.EQ.2).AND.
     .          machsuper2.AND.machcount2.LT.5) THEN
c     .          machsuper2.AND.machcount2.LE.6) THEN

              iterate = .TRUE.


              IF (imaginary2) THEN

                IF (machadjust2.LT.1.0) 
     .            machadjust2 = 1.0 + 0.5 * (1.0 - machadjust2)
                
              ELSEIF (lastp6.LT.p6.AND.machadjust2.LT.1.0) THEN
                WRITE(0,*) 'FUNK - something is wrong'
                machadjust2 = 1.0 + 0.5 * (1.0 - machadjust2)
                p6 = HI
                STOP 'WAKADO!'

              ELSEIF (ABS(1.0-machadjust2).LT.1.0E-5) THEN
c...            Seems limiting machcount to 6 is not good enough when the target
c               is only slightly supersonic:
                iterate = .FALSE.

              ELSEIF (.TRUE.) THEN
                IF (machadjust2.GT.1.0) THEN
                   machadjust2 = 1.0 - 0.5 * (machadjust2 - 1.0)
                   machcount2  = machcount2 + 1
                ELSEIF (machadjust2.LT.1.0.AND.machno.LT.0.5) THEN
c...              Turn this silly business off and set target conditions to
c                 give Mach number unity since something else has
c                 made the outer target no longer supersonic:
                  machcount2 = 0
                  machsuper2 = .FALSE.
                  machadjust2 = 1.1
                  WRITE(0,*) 'CHECK:',knds(idds(ir,1)),holdn2
                  WRITE(0,*) '      ',kvds(idds(ir,1)),holdv2
 
                  CALL SaveSolution
                  STOP 'DFDC'

c...tmp:
c                  knds(idds(ir,2)) = holdn2 * machadjust
c                  kvds(idds(ir,2)) = holdv2 / machadjust
                   n(6) = n(6) * machadjust2
                   v(6) = v(6) / machadjust2
                ENDIF
              ENDIF
c...          Adjust target conditions to account for the change in Mach number
c             at the target:
              n(6) = n(6) / machadjust2
              v(6) = v(6) * machadjust2

c              WRITE(0,*) 'MACH OUTER:',
c     .                   ir,imaginary2,machadjust2,machcount2




c              knds(id0) = n(0)
c              knds(id6) = n(6)
c              kvds(id0) = v(0)
c              kvds(id6) = v(6)
c              DO ir = 1, nrs
c                DO ik = 1, nks(ir)
c                  pinion(ik,ir) = osmion(ik,ir)
c                  pinrec(ik,ir) = osmrec(ik,ir)
c                ENDDO
c              ENDDO              
c              CALL DumpGrid('FUSSING ABOUT BAD')

            ENDIF

            IF (s28output) THEN
              WRITE(0,'(A,I6,2L2,2X,2F12.6,1P,3E16.6,0P,F10.2,I6)') 
     .          'MACH:',ir,machsuper1,imaginary1,machadjust1,
     .          lastadjust1,
     .          n(0),v(0),p0,mach1,machcount1
              WRITE(0,'(A,6X,2L2,2X,2F12.6,1P,3E16.6,0P,F10.2,I6)') 
     .          '    :',machsuper2,imaginary2,machadjust2,lastadjust2,
     .          n(6),v(6),p6,mach2,machcount2
            ENDIF


          ELSE
c         *OLD CODE DELETED*
            STOP 'TRYING TO EXECUTE OLD CODE FOR SOME REASON'
          ENDIF

c...CHANGE
c          csmax = 0.0
c          DO ik = ik3, ik5
c            csmax = MAX(csmax,
c     .              ABS(kvhs(ik,ir))/GetCs(ktebs(ik,ir),ktibs(ik,ir))) 
c          ENDDO
c          csmax = MAX(csmax,
c     .            ABS(kvds(idds(ir,1)))/GetCs(kteds(idds(ir,1)),
c     .                                        ktids(idds(ir,1))))

          IF (Clock2()-timeint(6).GT.45) THEN
            WRITE(0     ,*) 'TIME: ',NINT(Clock2() - timeint(6))
            WRITE(PINOUT,*) 'TIME: ',NINT(Clock2() - timeint(6))
            CALL SaveSolution
            CALL DumpGrid('SOL28: Out of time')
c            CALL ER('SOL28','Out of time',*99)
          ENDIF



c...      Solve plasma:
        ENDDO


c            frac = 2.0 * (1.0 - machadjust1)
c        IF (new.AND.sloutput) THEN
        IF (new.AND.sloutput) THEN
c         WRITE(0,*)
c         WRITE(0,*)
c         WRITE(0,*)
c         WRITE(0,*)
          WRITE(0,'(A,I6,2L2,2X,2F12.6,1P,3E10.2,0P,F10.2,I6)') 
     .      'MACH:',ir,machsuper1,imaginary1,machadjust1,lastadjust1,
     .      n(0),v(0),p0,mach1,machcount1
          WRITE(0,'(A,6X,2L2,2X,2F12.6,1P,3E10.2,0P,F10.2,I6)') 
     .      '    :',machsuper2,imaginary2,machadjust2,lastadjust2,
     .      n(6),v(6),p6,mach2,machcount2
        ENDIF

        csmax1 = 99.0
        csmax2 = 99.0
c        machval = 0.95
        machval = 0.98
        IF (scalerec) THEN

c                  WRITE(0,*) 'MACH1: A',machsuper1

          IF (detach1.AND.machsuper1) THEN

c                  WRITE(0,*) 'MACH1: A',machadjust1,imaginary1 ! * LEFT OFF * 

            IF     (machadjust1.LT.0.01) THEN
c...          Recombination turned off but problem not resolved so give up:
              csmax1 = machval
            ELSEIF (imaginary1) THEN
              csmax1 = 100.0
            ELSE
              csmax1 = 0.0
              DO ik = ik1, ik3
c...            Only check flow toward the inner target:
c...  UGLY: SEE BELOW
c                  WRITE(0,*) 'MACH1:',ik,ir,kss(ik,ir),0.35*ksmaxs(ir)

                IF (s28superdet_.EQ.4.AND.kss(ik,ir).GT.0.35*ksmaxs(ir)) 
     .            CYCLE

                IF (iflexopt(5).EQ.20.AND.ir.LT.20) THEN
                  IF (kvhs(ik,ir).LT.0.0) THEN
                    mach1 = GetMach(kvhs(ik,ir),ktebs(ik,ir),
     .                              ktibs(ik,ir))
                    csmax1 = MAX(csmax1,mach1)
c                    WRITE(0,*) 'MACH2:',ik,ir,mach1,csmax1
c                  ELSE
c                    EXIT
                  ENDIF
                ELSE
                  IF (kvhs(ik,ir).LT.0.0) THEN
                    mach1 = GetMach(kvhs(ik,ir),ktebs(ik,ir),
     .                              ktibs(ik,ir))
                    csmax1 = MAX(csmax1,mach1)
c                    WRITE(0,*) 'MACH2:',ik,ir,mach1,csmax1
                  ELSE
                    EXIT
                  ENDIF
                ENDIF
              ENDDO
              IF (csmax1.EQ.0.0) THEN
                WRITE(0     ,*) 'SOL28: NO VALID CSMAX VALUE FOUND'
                WRITE(PINOUT,*) 'SOL28: NO VALID CSMAX VALUE FOUND'
                CALL DUMPGRID('CSMAX PROBLEMS')
                STOP 'PROGRAM HALTED 1'
              ENDIF
            ENDIF
          ENDIF

          IF (detach2.AND.machsuper2) THEN
            IF     (machadjust2.LT.0.01) THEN
c...          Recombination turned off but problem not resolved so give up:
              csmax2 = machval
            ELSEIF (imaginary2) THEN
              csmax2 = 100.0
            ELSE
              csmax2 = 0.0
              DO ik = ik5, ik3+1, -1
c...            Only check flow toward the outer target:
c...  UGLY: THIS CHECK IS NECESSARY SINCE MASSIVE FLOW REVERSAL CAN PRODUCE LARGE
c           MACH NUMBER FAR UPSTREAM.  BETTER TO CONTROL THE IONISATION SOURCES
c           THEN PUT THE ARBITRARY 0.65 RESTRICTION HERE (AND IN CALCUALTEPLASMA)
c           BUT THIS IS BEST FOR DEVELOPMENT
                IF (s28superdet_.EQ.4.AND.kss(ik,ir).LT.0.65*ksmaxs(ir)) 
     .            CYCLE
                IF (kvhs(ik,ir).GT.0.0) THEN
                  mach2 = GetMach(kvhs(ik,ir),ktebs(ik,ir),ktibs(ik,ir))
                  csmax2 = MAX(csmax2,mach2)
                ELSE
                  EXIT
                ENDIF
              ENDDO
              IF (csmax2.EQ.0.0) THEN
                WRITE(0     ,*) 'SOL28: NO VALID CSMAX VALUE FOUND'
                WRITE(PINOUT,*) 'SOL28: NO VALID CSMAX VALUE FOUND'
                CALL DumpGrid('PROGRAM HALTED 2')
                STOP 'PROGRAM HALTED 2' 
              ENDIF
            ENDIF
          ENDIF

        ENDIF

c...      Determine if the solution should continue to iterate:
c
c
c
c        WRITE(0,*) 'CSMAX:',csmax1,csmax2,count

        IF     (scalerec.AND.detach1.AND.machsuper1.AND.
     .          .NOT.imaginary1.AND.machadjust1.GT.1.0) THEN

c          WRITE(0,*)
          WRITE(0,*) '>>>> TURNING OFF MACHSUPER1'
c          WRITE(0,*)         

          machadjust1 = 1.0
          machsuper1 = .FALSE.

        ELSEIF (scalerec.AND.detach2.AND.machsuper2.AND.
     .          .NOT.imaginary2.AND.machadjust2.GT.1.0) THEN

c          WRITE(0,*)
          WRITE(0,*) '>>>> TURNING OFF MACHSUPER2'
c          WRITE(0,*)         

          machadjust2 = 1.0
          machsuper2 = .FALSE.


        ELSEIF (scalerec.AND.
     .    ((detach1.AND.machsuper1.AND.ABS(machval-csmax1).GT.0.01).OR.
     .     (detach2.AND.machsuper2.AND.ABS(machval-csmax2).GT.0.02)))
c     .     (detach2.AND.machsuper2.AND.ABS(machval-csmax2).GT.0.01)))
     .    THEN

          IF (detach1.AND.machsuper1.AND.ABS(machval-csmax1).GT.0.01)
     .      THEN

            IF (iflexopt(5).EQ.20.AND.ir.LT.20) THEN
              WRITE(0,*) 'DEBUG:',csmax1,machval,lastadjust1

              IF ((csmax1.LT.machval.AND.lastadjust1.LE.1.0).OR.
     .            (csmax1.GT.machval.AND.lastadjust1.GE.1.0)) THEN
c...        
                lastadjust1 = 1.0 + 0.67 * (1.0 - lastadjust1) 
                machcount1 = 0
              ENDIF
              machcount1 = machcount1 + 1
c...          Speed things up:
              IF (MOD(machcount1,5).EQ.0) lastadjust1 = lastadjust1**2
              machadjust1 = machadjust1 * lastadjust1

            ELSE
              IF ((csmax1.LT.machval.AND.lastadjust1.LE.1.0).OR.
     .            (csmax1.GT.machval.AND.lastadjust1.GE.1.0)) THEN
c...        
                lastadjust1 = 1.0 + 0.67 * (1.0 - lastadjust1) 
                machcount1 = 0
              ENDIF
              machcount1 = machcount1 + 1
c...          Speed things up:
              IF (MOD(machcount1,5).EQ.0) lastadjust1 = lastadjust1**2
              machadjust1 = machadjust1 * lastadjust1
            ENDIF


          ENDIF

          IF (detach2.AND.machsuper2.AND.ABS(machval-csmax2).GT.0.02)
c          IF (detach2.AND.machsuper2.AND.ABS(machval-csmax2).GT.0.01)
     .      THEN
            IF     (lastadjust2.EQ.1.0) THEN
c... 
              knds(id0) = n(0)
              knds(id6) = n(6)
              kvds(id0) = v(0)
              kvds(id6) = v(6)
              DO ir = 1, nrs
                DO ik = 1, nks(ir)
                  pinion(ik,ir) = osmion(ik,ir)
                  pinrec(ik,ir) = osmrec(ik,ir)
                ENDDO
              ENDDO              
              CALL DumpGrid('LASTADJUST2=1.0 -- BAD')
c              STOP 'LASTADJUST2=1.0 -- BAD'
            ELSEIF ((csmax2.LT.machval.AND.lastadjust2.LT.1.0).OR.
     .              (csmax2.GT.machval.AND.lastadjust2.GT.1.0)) THEN
c...        
              lastadjust2 = 1.0 + 0.67 * (1.0 - lastadjust2) 
              machcount2 = 1
            ENDIF
            machcount2 = machcount2 + 1
c...        Speed things up:
            IF (MOD(machcount2,10).EQ.0) lastadjust2 = lastadjust2**2
c            IF (MOD(machcount2,5).EQ.0) lastadjust2 = lastadjust2**2
            machadjust2 = machadjust2 * lastadjust2
          ENDIF

          count = 1
          tmpcount = tmpcount + 1

c          IF (tmpcount.EQ.2) THEN
          IF (tmpcount.EQ.300) THEN
            knds(id0) = n(0)
            knds(id6) = n(6)
            kvds(id0) = v(0)
            kvds(id6) = v(6)
            DO ir = 1, nrs
              DO ik = 1, nks(ir)
                pinion(ik,ir) = osmion(ik,ir)
                pinrec(ik,ir) = osmrec(ik,ir)
              ENDDO
            ENDDO              
            CALL DumpGrid('SCALING RECOMBINATION SOURCE')
          ENDIF

          IF (sloutput) THEN
            WRITE(0,*) '>>>> SCALE:',machadjust1,lastadjust1,
     .                 machcount1,tmpcount
            WRITE(0,*) '          :',machadjust2,lastadjust2,
     .                 machcount2
            WRITE(0,*)
          ENDIF

c          IF (tmpcount.EQ.90.OR.tmpcount.EQ.91) CALL SaveSolution


c        ELSEIF (new.AND.count.LT.1) THEN
        ELSEIF (new.AND.
     .          (count.LT.2.OR.
     .           (count.LT.5.AND.
     .           ((ringtype(ir).EQ.SOL1.AND.s28rec        .GT.0.AND.
     .                                      s28superdet   .EQ.0).OR.
     .            (ringtype(ir).EQ.SOL1.AND.s28rec        .GT.0.AND.
     .                                      s28superdet   .EQ.2).OR.
     .            (ringtype(ir).EQ.PFZ .AND.s28recpfz     .GT.0.AND.
     .                                      s28superdetpfz.EQ.0)).AND.
     .          (detach1.OR.detach2)).OR.
c...             Relax cf drift:
     .           (count.LT.10.AND.(s28cfpdrft   .EQ.-2.OR.
     .                             s28cfppfzdrft.EQ.-2.OR.
     .                             s28cfpdrft   .EQ.-3.OR.
     .                             s28cfppfzdrft.EQ.-3))))  THEN
c...      Needed for relaxation of recombination source:
          IF (sloutput) THEN
c            WRITE(0,*)
            WRITE(0,*) '>>>> COUNT:',count
c            WRITE(0,*)
          ENDIF
c            WRITE(0,*) '>>>> COUNT:',count

          IF (.FALSE..AND.count.EQ.1) THEN
            knds(id0) = n(0)
            knds(id6) = n(6)
            kvds(id0) = v(0)
            kvds(id6) = v(6)
            DO ir = 1, nrs
              DO ik = 1, nks(ir)
                pinion(ik,ir) = osmion(ik,ir)
                pinrec(ik,ir) = osmrec(ik,ir)
              ENDDO
            ENDDO              
            CALL DumpGrid('FUSSING ABOUT')
          ENDIF

          count = count + 1

          IF (s28cfpdrft.EQ.-2.OR.s28cfppfzdrft.EQ.-2) THEN
            DO ik = 1, nks(ir)
              osmcfpflx(ik,ir,2) = (1.0-rel_frac) * osmcfpflx(ik,ir,2) + 
     .                                  rel_frac  * tmpflx   (ik,ir)
              osmcfpflx(ik,ir,1) = osmcfpflx(ik,ir,2)
            ENDDO
          ENDIF 

          IF (s28cfpdrft.EQ.-3.OR.s28cfppfzdrft.EQ.-3) THEN

            DO ik = 1, nks(ir)
              tmpflx1 = tmpflx(ik,ir)

              IF (ik.GT.nks(ir)/2) THEN

                ir1 = irouts(ik,ir)
                tmpflx2 = 0.0
                DO ik1 = 1, nks(ir1)-1
                  IF (thetag(ik,ir).GT.thetag(ik1  ,ir1).AND.
     .                thetag(ik,ir).LE.thetag(ik1+1,ir1)) THEN

                    frac = (thetag(ik   ,ir ) - thetag(ik1,ir1)) /
     .                     (thetag(ik1+1,ir1) - thetag(ik1,ir1))                    

                    tmpflx2 = (1.0-frac) * tmpflx(ik1  ,ir1) +
     .                             frac  * tmpflx(ik1+1,ir1)

                    WRITE(0,'(A,4I6,1P,2E10.2)') 
     .                 'TMPFLX:',ik,ir,ik1,ir1,tmpflx1,tmpflx(ik1  ,ir1)
                  ENDIF
                ENDDO




              ELSE
                tmpflx2 = 0.0
              ENDIF


              osmcfpflx(ik,ir,2) = (1.0-rel_frac) * osmcfpflx(ik,ir,2) + 
     .                                  rel_frac  * (tmpflx1 - tmpflx2)
            ENDDO

          ENDIF 

c          IF (MOD(count,10).EQ.0) CALL SaveSolution

c          STOP 'HERE IN COUNT'
c *TEMP*
c          CALL SaveSolution
c          IF (count.EQ.10) CALL DumpGrid('drifts...')

        ELSEIF (.NOT.new.AND.count.LT.10) THEN
c...      Give the solution some time to relax -- necessary?:  LIKELY JUST FOR RECOMBINATION? CHEKC FOR NEW?
c          WRITE(0,*) 'COUNTING',count
c...WHY IS THIS HERE?
          IF (machadjust.LT.0.001) machadjustL2 = machadjustL2 * 10.0
          count = count + 1

        ELSEIF (addtemp.NE.1.0.AND.addtemp.NE.99.0.AND.
     .          .NOT.tempset.AND.osm_store.EQ.-1) THEN

         STOP 'NOT IN USE YET, WATCH OUT FOR ADDMP'

c        ELSEIF (addtemp.NE.1.0.AND.addtemp.NE.99.0.AND..NOT.tempset.AND.
c     .          .NOT.(loadsources.OR.osm_store.NE.-1.OR.
c     .                (relreset.GE.1.AND.rel_count.GT.1).OR.
c     .                (relreset.EQ.0.AND.rel_count.GT.0))) THEN

c...      Add ADDTEMP to the density peak temperature and solve again:
          WRITE(0,'(A,I6,1P,3E10.2,0P,L2)') 
     .      'SOL28 VALUE SET:',ir,csmax,machval,teval,scalerec
          count = 0
          tempte = te(1) * addtemp
          tempset = .TRUE.
c...      Turn off 
          scalerec = .FALSE.    

          IF (loadsources) STOP 'NO LOADSOURCES=.TRUE.!'

        ELSEIF (.FALSE..AND.fitline) THEN
c...      Fit the solution to a specified line radiation
c         profile:

c...      Base solution has been found, now turn off calculation of
c         density peak momentum loss:
c...THIS SEEMS TO SUCK
c          calcmom1 = .FALSE.

          tmpcount2 = tmpcount2 + 1

c...      Make sure this is off -- may need another way of making sure that
c         the target temperature is not floating in the future -- it must be
c         fixed by this point:
          scalerec = .FALSE.    
     
          WRITE(0,*) 'PRESSURE:',te(0),scalerec
          IF (tmpcount2.EQ.1) THEN
c...        Initialize            
            lastik = 0
          ELSEIF (tmpcount2.EQ.20) THEN
            kteds(idds(ir,2)) = te(0)
            ktids(idds(ir,2)) = te(0) * timul
            CALL DUMPGRID('DOING ADDMP')
          ENDIF

          DO ik = ik1+1, ik2-1
c...        Process cells between IK1 and IK2:
            CALL CheckLineMatch(ik,ir,linedat(ik,ir),1,H_BGAMMA,tarden)
            IF (tarden.NE.knbs(ik,ir)) THEN
c...          Modify density to improve agreement with line radiation:
              IF (ik.NE.lastik) THEN
c...            Initialize momentum source adjustment for this cell:
                
              ENDIF
              deltas = kss(ik,ir) - ksb(ik-1,ir)
              deltap = CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                              ktibs(ik,ir),kvhs (ik,ir))*ECH *
     .                 ((tarden / knbs(ik,ir)) - 1.0)

              WRITE(0,*) 'PRESSURE:',ik,ir,deltap,deltas
c              DO i1 = 1, ik1
c                WRITE(0,*) 'MOM:',i1,osmmp(i1,ir)
c              ENDDO
c...          Likely want to save this, so better improve declaration:
              addmp(ik) = deltap / deltas


              lastik = ik
c              WRITE(0,*) 'TARDEN .NE. 0',ik,ir
              EXIT
            ENDIF
          ENDDO

          IF (ik.EQ.ik2) THEN
c...        Sufficient agreement for all cells in range:
            fitline = .FALSE.
          ENDIF


          IF (tmpcount2.EQ.10) THEN
            fitline = .FALSE.
            tmpcount2 = 0
          ENDIF
c          WRITE(0,*) 'FITLINE:',fitline,tmpcount2

        ELSE

          IF (sloutput) THEN
            WRITE(0     ,*) 'SOL28 VALUE:',ir,csmax,machval,te(1)
            WRITE(PINOUT,*) 'SOL28 VALUE:',ir,csmax,machval,te(1)
          ENDIF

          IF (detach1) s28recfrac(IKLO,ir) = machadjust1
          IF (detach2) s28recfrac(IKHI,ir) = machadjust2


          IF ((ringtype(ir).EQ.SOL1.AND.s28rec   .GT.0).OR.
     .        (ringtype(ir).EQ.PFZ .AND.s28recpfz.GT.0)) THEN
            DO ik = 1, nks(ir)
c...          Always over-ride PINREC with OSMREC:
              pinrec(ik,ir) = osmrec(ik,ir) 
c...          Modifying KVHS:
c              kvhs(ik,ir) = -GetCs(ktibs(ik,ir),ktebs(ik,ir)) 
            ENDDO
          ENDIF

c          WRITE(0,*) '--->',ir,ringtype(ir),SOL1,s28mom

          IF (.NOT.loadsources) THEN
 
            IF ((ringtype(ir).EQ.SOL1.AND.s28ion   .GT.0).OR.
     .          (ringtype(ir).EQ.PFZ .AND.s28ionpfz.GT.0)) THEN
c...
              DO ik = 1, nks(ir)
                pinion(ik,ir) = osmion(ik,ir) 
                DO i1 = 1, 3               
                  pindata(ik,ir,H_ION1+i1-1) = osmion(ik,ir) / 3.0
                ENDDO
              ENDDO
            ENDIF

            IF ((ringtype(ir).EQ.SOL1.AND.s28mom   .GT.0).OR.
     .          (ringtype(ir).EQ.PFZ .AND.s28mompfz.GT.0)) THEN
c...
              DO ik = 1, nks(ir)
                pinmp(ik,ir) = osmmp(ik,ir) 
              ENDDO
            ENDIF

          ENDIF

          IF (eiropacity.EQ.-3.AND.
     .        ((ringtype(ir).EQ.SOL1.AND.s28rec   .EQ.0).OR.
     .         (ringtype(ir).EQ.PFZ .AND.s28recpfz.EQ.0))) THEN

c...        Make volume recombination very small in EIRENE (not zero, since that
c           will cause problems):
            Stop 'crap 2'

            DO ik = ik0, ik6           
              mulrec = 1.0E-05
            ENDDO
            DO ik = 1, nks(ir)
              pinrec(ik,ir) = 0.0
            ENDDO

          ELSEIF (eiropacity.EQ.-3.AND.
     .            (ringtype(ir).EQ.SOL1.AND.s28superdet   .EQ.1).OR.
     .            (ringtype(ir).EQ.PFZ .AND.s28superdetpfz.EQ.1)) THEN
c...        Set MULREC according to MACHADJUST1,2 for s28machsuper.EQ.2:
            CALL RSet(mulrec(1,ir),MAXNKS,1.0)
            IF (detach1) THEN
              DO ik = ik0, ik3
                mulrec(ik,ir) = machadjust1
              ENDDO
            ENDIF
            IF (detach2) THEN
              DO ik = ik3+1, ik6
                mulrec(ik,ir) = machadjust2
              ENDDO
            ENDIF
          ENDIF

          DO ik = 1, nks(ir)
            oldknbs (ik,ir) = knbs (ik,ir)
            oldktebs(ik,ir) = ktebs(ik,ir)
            oldktibs(ik,ir) = ktibs(ik,ir)
            oldkvhs (ik,ir) = kvhs (ik,ir)
          ENDDO

          count = 0
          tmpcount = 0

          tempset = .FALSE.

          kteds(idds(ir,2)) = te(0)
          ktids(idds(ir,2)) = ti(0)
          kteds(idds(ir,1)) = te(6)
          ktids(idds(ir,1)) = ti(6)

          IF (new) THEN
            knds(id0) = n(0)
            knds(id6) = n(6)
            kvds(id0) = v(0)
            kvds(id6) = v(6)
          ENDIF

          lockrec1 = .FALSE.
          lockrec2 = .FALSE.
          tmprec = 0.0
          machsuper1 = .FALSE.
          machsuper2 = .FALSE.
          lastadjust1 = 0.9
          lastadjust2 = 0.9
          machadjust1 = 1.0
          machadjust2 = 1.0



c...      If PIN is not being called, then calculated Dgamma excitation and 
c         recombination components:
          IF (cpinopt.EQ.0) THEN
            CALL CalcRecLineEmissions(ir)
          ENDIF



          ir = ir + 1


        ENDIF


      ENDDO
c
c     =================================
c              END OF MAIN LOOP
c     =================================
c



      IF (imaginaryL) THEN
        WRITE(0,*)
        WRITE(0,*) '*** LOW IMAGINARY DETECTED ***'
        WRITE(0,*)
      ENDIF

      WRITE(PINOUT,*)    
      WRITE(PINOUT,*) 'END OF SOL28:'
      CALL AnalyseSolution(PINOUT)

      IF (s28fit.EQ.1) THEN
        DEALLOCATE(linedat)
      ENDIF

c... 
c      IF (firstcall) THEN
c        IF (s28ion   .GT.0.AND..NOT.s28ionset   ) s28ionset    = .TRUE.
c        IF (s28ionpfz.GT.0.AND..NOT.s28ionsetpfz) s28ionsetpfz = .TRUE.
c        IF (s28mom   .GT.0.AND..NOT.s28momset   ) s28momset    = .TRUE.
c        IF (s28mompfz.GT.0.AND..NOT.s28momsetpfz) s28momsetpfz = .TRUE.
c      ENDIF


      IF (s28sym.EQ.1) THEN
c      IF (.TRUE.) THEN
c...    Symmeterize NON-SOL28 rings:
        IF (sloutput) THEN
          WRITE(0,*)
          WRITE(0,*) '***************************'
          WRITE(0,*) '* SOL28 SYMMETRY APPLIED! *'
          WRITE(0,*) '***************************'
          WRITE(0,*) 
        ENDIF
        DO ir = irsep, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          IF (osm_model(IKLO,ir).NE.28) THEN
            DO ik1 = 1, nks(ir) 
              IF (kss(ik1,ir).LT.0.5*ksmaxs(ir)) THEN
                DO ik2 = nks(ir), 1, -1
                  IF ((ksmaxs(ir)-kss(ik2,ir)).GT.kss(ik1,ir)) EXIT
                ENDDO
                ktebs(ik1,ir) =  ktebs(ik2,ir)
                ktibs(ik1,ir) =  ktibs(ik2,ir)
                knbs (ik1,ir) =  knbs (ik2,ir)
                kvhs (ik1,ir) = -kvhs (ik2,ir)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF




      IF (.FALSE..AND.loadsources) THEN
        knds(id0) = n(0)
        knds(id6) = n(6)
        kvds(id0) = v(0)
        kvds(id6) = v(6)
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            pinion(ik,ir) = osmion(ik,ir)
          ENDDO
        ENDDO              
        CALL DumpGrid('AFTER 1 CALL TO EIRENE')
      ENDIF


      IF (firstcall) THEN
c....   Make sure that FIRSTCALL is not set to .FALSE. if
c       there are still more entries in the ... array:
        firstcall = .FALSE.
        DO ir = irsep, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          IF (osm_model(IKLO,ir).EQ.28.AND..NOT.s28done(ir)) 
     .      firstcall = .TRUE.
        ENDDO

        IF (.NOT.firstcall) s28done = .FALSE.

      ENDIF

      IF (.FALSE..AND.
     .    .NOT.firstcall.AND.new.AND.tagpinatom.AND.nitersol.EQ.1) THEN
c      IF (new.AND.tagpinatom) THEN

        stoprun = .TRUE.
        DO ir = irsep, nrs
          IF (s28done(ir).OR.idring(ir).EQ.BOUNDARY) CYCLE

          IF (osm_model(IKLO,ir).EQ.28) stoprun = .FALSE.
        ENDDO

        IF (stoprun) CALL DumpGrid('CHECKING PIN DATA SOLUTION')
      ENDIF

c...  Clear arrays:
      IF (ALLOCATED(tmpflx)) DEALLOCATE(tmpflx)


      WRITE(0,*) 'SOL28 RINGS:',irs,ire


      RETURN
99    STOP
      END
c
c
c
c
c
c
c
c
c
c
c
