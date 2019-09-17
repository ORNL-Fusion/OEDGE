c
c ======================================================================
c taken from tau.d6a
c
      subroutine calc_wallprad_SL(nizs)
      implicit none
c
      integer nizs 
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cedge2d'
      include 'cadas'
      include 'cnoco'                  
      include 'dynam3'
      include 'pindata'
      include 'printopt'      
      include 'slcom' 
c     
c     CALC_WALLPRAD:
c
c     This routine calculates the total radiated power onto each wall segment.
c     It also calculates the total radiated H and IMP power. 
c     In order to calculate the amount on each wall segment it has to 
c     allocate the amount from each cell that would strike each wall elememt.
c     This routine uses the data in HPOWLS and POWLS to make these calculations.
C     The POWLS data is multiplied by the ABSFAC quantity to get absolute scalings
c     unless ABSFAC is zero. 
c
      integer in,ir,ik,iz,it
      real*8 imp_prad,h_prad,tot_prad
      real*8 rcent,zcent
c
      integer n_int,max_int,backward,forward,ivert
      parameter(max_int=10,backward=-1,forward=1)
      real*8 r_int(max_int),z_int(max_int)
      integer n_intersections(maxpts+1),next_vert
      external next_vert 
c
      real*8 datan2c
      real*8 wall_angle(maxpts+1),net_angle(maxpts),result_angle
      external datan2c
c
      integer nstart,nend,ntest,segtype
      real*8 rtest,ztest
c     slmod begin
      integer iatom,iz_max,iz_sft,atno,iz_max_ultimate
      integer kfail
      real*4  pmass
      character adas_id*80,adas_yr_i*2,adas_yr_h*2
c     slmod end
c
      
c     Initialization 
c
      REAL wallprad2(maxpts+7,3)

      wallprad2 = 0.0

      write(0,*) 'nizs',nizs,absfac

c...  Initialize nc calculation, where 20.0=Ne:
      powls  = 0.0
      lines  = 0.0
      hpowls = 0.0
      hlines = 0.0

      iz_max_ultimate = -1
      
      DO iatom = 2, 2

        IF     (iatom.EQ.1) THEN 
          pmass  = 4.0
          iz_max = 2
          iz_sft = 0
          atno = 2
          adas_yr_i = '89'

          adas_id = '*'
        ELSEIF (iatom.EQ.2) THEN 
          pmass  = 20.0
          iz_max = 10
          iz_sft = 2
          atno = 10
          adas_yr_i = '89'
        ELSE
          STOP 'no be here dude' 
        ENDIF

        iz_max_ultimate = MAX(iz_max,iz_max_ultimate)
        
        adas_id = '*'
        adas_yr_h = '89'          
        
c       From div6/src/div.f:
        pnzs = 0.0
        
        DO IR = 1, nrs ! irsep, irsep ! 1, NRS
          IF (idring(ir).EQ.BOUNDARY) CYCLE
        
          IF (.TRUE.) THEN

            PTES  (1:nks(ir)) = KTEBS (1:nks(ir),IR)
            PNES  (1:nks(ir)) = KNBS  (1:nks(ir),IR) * 1.E-6 * RIZB
            PDVOLS(1:nks(ir)) = KAREAS(1:nks(ir),IR) * 1.E+6
            
            if (ir.eq.irsep) then
              write(0,*) ' > ',iatom,1+iz_sft,iz_max+iz_sft
            endif
    
            DO IZ = 1, iz_max
               PNZS(IZ+1,1,1:nks(ir)) = e2dnzs(1:nks(ir),ir,iz+iz_sft) * 
     .                                  1.E-6          
               IF (ir.Eq.irsep) THEN
                 write(0,*) '> ',50,ir,iz+iz_sft,e2dnzs(50,ir,iz+iz_sft)
               ENDIF
            end do

c    ------ CALCULATE POWER LOSS (W.CM**-3) AND LINE RAD (W.CM**-3), RECONVERT TO W.M**-3
c           VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY
c           GOES OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.

            CALL AATOM   (pmass,1,9,KFAIL)
            call ncrdlong(nks(ir))

            DO IK = 1, NKS(IR)
              DO IZ = 1, iz_max
                powls(IK,IR,IZ) = powls(IK,IR,IZ) +
     .                 MAX (0.0, PRADIS(1,IZ+1,1,IK)*1.0E6)
                lines(IK,IR,IZ) = lines(IK,IR,IZ) +
     .                 MAX (0.0, PRADIS(3,IZ+1,1,IK)*1.0E6)
              end do
            end do
            
          ELSE

            pnzsa = 0.0
             
            DO IK = 1, NKS(IR)
              PTESA(IK) = KTEBS(IK,IR)
              PNESA(IK) = KNBS(IK,IR) * RIZB
              PNBS(IK) = KNBS(IK,IR)
              PNHS(IK) = PINATOM(IK,IR)
            ENDDO
            DO IK = 1, NKS(IR)
               DO IZ = 1, iz_max
c              DO 1110 IZ = 0, NIZS
                 PNZSA(IK,IZ) = e2dnzs(ik,ir,iz+iz_sft) ! * 1.E-6
c                PNZSA(IK,IZ) = SNGL(DDLIMS(IK,IR,IZ))             
              ENDDO
            ENDDO
C
C------ GET POWER LOSS FROM ADAS DATA FILES. LOAD TOTAL LINE RADIATION
C------ INTO LINES AND ADD RECOMBINATION AND BREMSSTRAHLUNG POWER TO
C------ GET TOTAL RADIATIVE LOSSES
C
            call xxuid(adas_id)
            ICLASS = 5
            DO IZ = 1, iz_max
              CALL ADASRD(adas_yr_i,atno,IZ+1,ICLASS,NKS(IR),PTESA,
     +                    PNESA,PCOEF(1,IZ+1))
              DO IK = 1, NKS(IR)
                LINES(IK,IR,IZ)=PCOEF(IK,IZ+1)*PNESA(IK)*PNZSA(IK,IZ)
                powls(IK,IR,IZ)=LINES(IK,IR,IZ)
              ENDDO
            ENDDO

            ICLASS = 4

            DO IZ = 1, iz_max
              CALL ADASRD(adas_yr_i,atno,IZ,ICLASS,NKS(IR),PTESA,
     +                    PNESA,PCOEF(1,IZ))
              DO IK = 1, NKS(IR)
                POWLS(IK,IR,IZ) = POWLS(IK,IR,IZ) +
     +                            PCOEF(IK,IZ)*PNESA(IK)*PNZSA(IK,IZ)
              ENDDO
            ENDDO
c
c------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1)
c
c            DO 1160 IK = 1, NKS(IR)
c              IF (DDLIMS(IK,IR,0).LE.0.0) THEN
c                POWLS(IK,IR,-1) = 0.0
c                LINES(IK,IR,-1) = 0.0
c              ELSE
c                POWLS(IK,IR,-1) = POWLS(IK,IR,0) *
c         +                        SNGL(DDLIMS(IK,IR,-1) / DDLIMS(IK,IR,0))
c                LINES(IK,IR,-1) = LINES(IK,IR,0) *
c         +                        SNGL(DDLIMS(IK,IR,-1) / DDLIMS(IK,IR,0))
cc       
cc                write (6,'(a,3i5,3g16.8)') 'Debug POW:',ir,ik,iz,
cc         >              pcoef(ik,iz),pnesa(ik),pnzsa(ik,iz)
cc                write (6,'(a,15x,3g16.8)') '      POW:',
cc         >           lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
cc       
c              ENDIF
c 1160       CONTINUE
       
          ENDIF
c
c...      Hydrogen:
          call xxuid(adas_id)
          ICLASS = 5
          CALL ADASRD(adas_yr_h,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
          DO IK = 1, NKS(IR)
            HLINES(IK,IR,0) = PCOEF(IK,1)*PNESA(IK)*PNHS(IK)
            HPOWLS(IK,IR,0) = HLINES(IK,IR,0)
          ENDDO
          ICLASS = 4
          CALL ADASRD(adas_yr_h,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
          DO IK = 1, NKS(IR)
            HPOWLS(IK,IR,1) = PCOEF(IK,1)*PNESA(IK)*PNBS(IK)
          ENDDO
          
        ENDDO  ! ir

      ENDDO  ! iatom


      write(0,*) ' scanning wall' 

c     
c     For every cell on the grid - this code must loop through every element of 
c     the wall.  
c                
      tot_prad = 0.0
      
      do ir = 1,nrs
        do ik = 1,nks(ir)

           imp_prad = 0.0
           do iz = 1, iz_max_ultimate
             imp_prad = imp_prad + powls(ik,ir,iz)
           end do
c
         !  if (absfac.gt.0.0) imp_prad = imp_prad * absfac

           h_prad = 0.0
           do iz = 0,1
             h_prad = h_prad + hpowls(ik,ir,iz)
           end do

c...       Convert to Watts:           
           imp_prad = imp_prad * kareas(ik,ir) * 2.0 * PI * rs(ik,ir) 
           h_prad   = h_prad   * kareas(ik,ir) * 2.0 * PI * rs(ik,ir) 
           tot_prad = h_prad + imp_prad

           if (tot_prad.eq.0.0) cycle

           wallprad2(maxpts+7,1) = wallprad2(maxpts+7,1) + imp_prad
           wallprad2(maxpts+7,2) = wallprad2(maxpts+7,2) + h_prad
           wallprad2(maxpts+7,3) = wallprad2(maxpts+7,3) + tot_prad

           rcent = rs(ik,ir)
           zcent = zs(ik,ir)  

           do in = 1,wallpts
c
c             Calculate initial intersections and whether LOS to each
c             vertex is clear. Last point and first point of the wall
c             should be identical
c
              n_intersections(in) = 0
c
              rtest = wallpt(in,20)
              ztest = wallpt(in,21)
c
              call calc_wall_intersections(ntest,max_int,r_int,z_int,
     >                     rcent,zcent,rtest,ztest,.true.)
c           
              n_intersections(in) = ntest 
c
c             Calculate wall_angle array
c
c             Get angles to the start and end of the wall segment
c 
              wall_angle(in) = datan2c(ztest-zcent,rtest-rcent)
c
c             Convert angles to 0 to 2 PI range  from -PI to PI
c
              if (wall_angle(in).lt.0.0)
     >            wall_angle(in)=wall_angle(in)+2.0*PI
c
           end do
c
c          Assign values to complete the wall
c
           n_intersections(wallpts+1) = n_intersections(1)
           wall_angle(wallpts+1) = wall_angle(1)
c
c          Calculate the net angle for each segment of the wall
c
           do in = 1,wallpts
c
              net_angle(in) = wall_angle(in+1)- wall_angle(in)
c
              if (net_angle(in).gt.PI) 
     >            net_angle(in) = net_angle(in) - 2.0 * PI
              if (net_angle(in).lt.-PI) 
     >            net_angle(in) = net_angle(in) + 2.0 * PI 
c
           end do 
c  
c          Basic LOS rules:
c          
c          1) Wall is listed in a basically clockwise order so assuming that all angles are positive in the range 0.0 to 2 PI then the angle to the leading point of a wall segment will always be less than the angle to the second point for a forward oriented wall segment. 
c          2) All wall elements that have some exposure to radiation from the cell being examined will be forward going assuming all radiation is from within the closed vessel. 
c          3) If the LOS from the cell to the ends of the wall segment cross the the wall at any point - then the view to that wall segment is assumed to be blocked.
c          4) If the LOS from the cell to the ends of the wall segment does not cross the wall then the LOS is considered unblocked and the entire segment is visible to the source. 
c          5) If either end of a backward going wall segment is blocked then none of that line segment can be seen from the cell. 
c          6) If the wall segment is forward going but one vertex is blocked then the wall element is partially obscured and the proportion of the segment exposed to the source cell is approximately calculated.  
c          7) Only clockwise oriented segments can receive radiation (net_angle<0)
c          8) The total exposure of a partially blocked segment is limited by the visible vertex and either the next or last visible vertex along the wall. 

           do in = 1,wallpts

              if (net_angle(in).lt.0.0) then 
c
c                If net_angle is less than zero then the wall element is forward going or clockwise relative to the observation point.
c            
c                Check number of intersections for the wall element vertices  
c
c                Totals only need updating for unobstructed and partially obstructed views.
c            
c                Wall segment is unobstructed:
                 if (n_intersections(in  ).eq.0.and.
     >               n_intersections(in+1).eq.0) then 
c
                    result_angle = -net_angle(in) 
c
c                Clockwise wall segment with one blocked vertex.
c                Need to calculate the required angle
c                Need angle of last unobstructed vertex
c
                 elseif (n_intersections(in  ).gt.0.and.
     >                   n_intersections(in+1).eq.0) then 
c     	    
                    ivert = next_vert(n_intersections,
     >                                wallpts+1,in,backward)
c     	    
                    result_angle = wall_angle(in+1)- wall_angle(ivert)
c     	    
                    if (result_angle.gt.PI) 
     >                  result_angle = result_angle - 2.0 * PI
                    if (result_angle.lt.-PI) 
     >                  result_angle = result_angle + 2.0 * PI 

                    result_angle = -result_angle

c                Need angle of next unobstructed vertex

                 elseif (n_intersections(in  ).eq.0.and.
     >                   n_intersections(in+1).gt.0) then
      	  
                    ivert = next_vert(n_intersections,
     >                                wallpts+1,in+1,forward)
c     	  
                    result_angle = wall_angle(ivert)- wall_angle(in)
c     	  
                    if (result_angle.gt.PI) 
     >                  result_angle = result_angle - 2.0 * PI
                    if (result_angle.lt.-PI) 
     >                  result_angle = result_angle + 2.0 * PI 

                    result_angle = -result_angle

c                LOS obstructed

                 elseif (n_intersections(in  ).gt.0.and.
     >                   n_intersections(in+1).gt.0) then

                     result_angle = 0.0 

                 endif
c
c                Calculate contributions to this wall segment
c
                 wallprad2(in,1) = wallprad2(in,1) 
     >                       + imp_prad * result_angle/(2.0*PI)
                 wallprad2(in,2) = wallprad2(in,2) 
     >                       + h_prad * result_angle/(2.0*PI)
                 wallprad2(in,3) = wallprad2(in,3) 
     >                       + tot_prad * result_angle/(2.0*PI)
              end if

           end do
        end do
      end do  
c
c     Print outs of results
c
      write(0,*) ' radiation output'
      
      do in = 1,wallpts
c
c       The following grand totals are calculated on a complete torus basis factoring in the major radius. 
c       
c       Check the tag for the wall element to see which type of segment it is:
c
c       1 = Target 1 (Outer target for X-point up grids)
c       4 = Target 2 (Inner target for X-point up grids)
c       7 = (and others 2,3) Main Vessel wall
c       8 = PFZ wall
c       9,10 = baffle segments - added to PFZ only applies for some JET grids 

        segtype = wallpt(in,16)
c
c       First Target (Outer for X-point up - inner for down)
c  
        if (segtype.eq.1) then          
          wallprad2(maxpts+1,1) = wallprad2(maxpts+1,1)+wallprad2(in,1)
          wallprad2(maxpts+1,2) = wallprad2(maxpts+1,2)+wallprad2(in,2)
          wallprad2(maxpts+1,3) = wallprad2(maxpts+1,3)+wallprad2(in,3)
c
c       Second Target (Inner for X-point up - Outer for down)
c  
        elseif (segtype.eq.4) then          
          wallprad2(maxpts+2,1) = wallprad2(maxpts+2,1)+wallprad2(in,1)
          wallprad2(maxpts+2,2) = wallprad2(maxpts+2,2)+wallprad2(in,2)
          wallprad2(maxpts+2,3) = wallprad2(maxpts+2,3)+wallprad2(in,3)
c
c       Main Vessel wall elements
c
        elseif (segtype.eq.7.or.segtype.eq.2.or.segtype.eq.3) then
          wallprad2(maxpts+3,1) = wallprad2(maxpts+3,1)+wallprad2(in,1)
          wallprad2(maxpts+3,2) = wallprad2(maxpts+3,2)+wallprad2(in,2)
          wallprad2(maxpts+3,3) = wallprad2(maxpts+3,3)+wallprad2(in,3)
c
c       Private Flux Zone wall elements   
c
        elseif (segtype.eq.8.or.segtype.eq.9.or.segtype.eq.10) then
          wallprad2(maxpts+4,1) = wallprad2(maxpts+4,1)+wallprad2(in,1)
          wallprad2(maxpts+4,2) = wallprad2(maxpts+4,2)+wallprad2(in,2)
          wallprad2(maxpts+4,3) = wallprad2(maxpts+4,3)+wallprad2(in,3)
        else
          wallprad2(maxpts+5,1) = wallprad2(maxpts+5,1)+wallprad2(in,1)
          wallprad2(maxpts+5,2) = wallprad2(maxpts+5,2)+wallprad2(in,2)
          wallprad2(maxpts+5,3) = wallprad2(maxpts+5,3)+wallprad2(in,3)
        endif

      end do 
c
c     Sum up grand totals
c
      do in = 1,5
         do it = 1,3
            wallprad2(maxpts+6,it) = wallprad2(maxpts+6 ,it) 
     >                             + wallprad2(maxpts+in,it)
         end do
      end do
c
c     Convert all of the wall element data to W/m2
      do in = 1,wallpts
         do it = 1,3
            if (wallpt(in,7).gt.0.0.and.wallpt(in,1).gt.0.0) then
               wallprad2(in,it) = wallprad2(in,it)
     >                           /(wallpt(in,7)*2.0*PI*wallpt(in,1))
            endif
         end do
      end do
c
c
c     Print summary to unit 6 
      write(0,*) 
      write(0,'(a)') 'WALL Radiation Flux Summary [W/m2]:'
      write(0,'(2(a4,2X),4A12)') 'ID','TYPE','IMPURITY','HYDROGEN',
     .                           'TOTAL'
c     Wall elements
      do in = 1,wallpts
         write(0,'(2(i4,2x),4(1x,g12.5))') in,int(wallpt(in,16)),
     >             (wallprad2(in,it),it=1,3)
      end do

      write(0,*) 
      write(0,*) 'Regional Radiation Totals in [MW]'
      write(0,*) 
      write(0,80) 'TOT  '//inner,(wallprad2(maxpts+1,it)/1.E6,it=1,3)
      write(0,80) 'TOT  '//outer,(wallprad2(maxpts+2,it)/1.E6,it=1,3)
      write(0,80) 'TOT   MAIN'  ,(wallprad2(maxpts+3,it)/1.E6,it=1,3)
      write(0,80) 'TOT    PFZ'  ,(wallprad2(maxpts+4,it)/1.E6,it=1,3)
      write(0,80) 'TOT   MISC'  ,(wallprad2(maxpts+5,it)/1.E6,it=1,3)
      write(0,*) 
      write(0,80) 'TOTAL  SEG'  ,(wallprad2(maxpts+6,it)/1.E6,it=1,3)
      write(0,80) 'TOTAL  SRC'  ,(wallprad2(maxpts+7,it)/1.E6,it=1,3)
      write(0,*) 

 80   FORMAT(a10,4(7x,f6.2))
c 80   FORMAT(a10,1P,4(1x,e12.5),0P)      

      return
      end
c     
c ======================================================================
c
c subroutine: AnalyseSolution
c
c ... assumes SetBounds has been called
c
      SUBROUTINE AnalyseSolution(fp)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INCLUDE 'solparams'
      INCLUDE 'solcommon'
      INCLUDE 'solswitch'

      INTEGER    IN_ION    ,IN_REC    ,IN_CFP,
     .           IN_PQE    ,IN_PEI    ,IN_CFE,
     .           IN_PQI    ,IN_RAD    ,IN_CFI
      PARAMETER (IN_ION = 1,IN_REC = 2,IN_CFP = 3,
     .           IN_PQE = 4,IN_PEI = 5,IN_CFE = 6,
     .           IN_PQI = 7,IN_RAD = 8,IN_CFI = 9)

      INTEGER GetModel,SymmetryPoint
      LOGICAL ChkInt
      REAL    CalcPressure,GetJsat,GetFlux,GetCs,ATan2C,GetHeatFlux,
     .        GetGamma

      INTEGER ik,ik1,ik2,ir,i1,i2,model1(MAXNRS),model2(MAXNRS),
     .        id,maxik1,maxik2,fp,ikm,model,ring,target
      REAL    parflx,eleflx,ionflx,p,p1,p2,jsat,area,
     .        partot,eletot,iontot,maxs,rdum(20),mach,beta,alpha,
     .        radpow(MAXNKS,MAXNRS),fact,cost,flux(2),deltar,deltaz,
     .        fluxmax(2),isat,qpara,gamma,sum_target,sum_total
      REAL*8  intg(20),ddum1,ddum2,a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,
     .        f1,f2,tab,tcd
      CHARACTER*7 irtag(0:MAXNRS)
      CHARACTER*1 note1,note2

      COMMON /OPTTEMP/ osm_matcht,forcet1
      INTEGER          osm_matcht,forcet1

      DOUBLE PRECISION EPS10
      PARAMETER       (EPS10=1.0D-10)

c..temp
      REAL    delta,r1,area2

c      IF (cioptg.NE.91) THEN
c        CALL WN('AnalyseSolution','Avoiding')
c        RETURN
c      ENDIF

c      WRITE(0,*) 'ANALYSING THE SOLUTION'

      CALL DB('Here in AnalyseSolution')

c
c jdemod - initialization
c
      qt = qtim 
c
c jdemod
c
c      fp = PINOUT

      DO ir = 1, MAXNRS
        irtag(ir) = '      '
      ENDDO

      irtag(irsep)   = 'IRSEP  '
      irtag(irbreak) = 'IRBREAK'
      irtag(irwall)  = 'IRWALL '
      irtag(irtrap)  = 'IRTRAP '
      irtag(nrs)     = 'NRS    '


      CALL HD(fp,'Solution analysis','SOLANAL-TITLE',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Solution analysis:'

      DO i1 = 1, 20
        intg(i1) = 0.0
      ENDDO

      DO ir = irsep, nrs
        model1(ir) = osm_model(IKLO,ir)
        model2(ir) = osm_model(IKHI,ir)
      ENDDO
c
c
c     Target parameters:
c
c
      CALL HD(fp,'Target conditions','SOLANAL-TARGETS',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Target conditions:'
      DO i1 = 2, 1, -1
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,A5,2A9,1X,2A12,2A10,1X,A10,A6,A12)')
     .    'ir','sol','psin','rho','jsat','ne','Te','Ti','v','M','p'
c        WRITE(fp,'(1X,A4,A5,2A9,1X,2A10,2A12,1X,A10,A6,A12)')
c     .    'ir','sol','psin','rho','Te','Ti','Ne','Jsat','v','M','p'
c        DO ir = irsep, nrs
        DO ring = 1, nrs-irsep+1
          IF (ring.LT.nrs-irtrap+1) THEN
            ir = irtrap + ring
          ELSE
            ir = ring - (nrs - irtrap + 1) + irsep 
          ENDIF
          IF (idring(ir).EQ.-1) CYCLE
          IF (i1.EQ.1) THEN
            model = model2(ir)
          ELSE
            model = model1(ir)
          ENDIF
          id   = idds(ir,i1)
          p    = CalcPressure(knds(id),kteds(id),ktids(id),kvds(id))
          jsat = GetJsat     (kteds(id),ktids(id),knds(id),kvds(id))
          mach = kvds(id) / GetCs(kteds(id),ktids(id))
          WRITE(fp,'(1X,I4,I5,2F9.5,1X,1P,2E12.4,0P,2F10.4,1X,1P,'//
     .             'E10.2,0P,F6.2,1P,E12.4,0P,1X,A,2F10.5,2X)')
     .      ir,model,psitarg(ir,1),rho(ir,CELL1),
     .      jsat,knds(id),kteds(id),ktids(id),kvds(id),mach,p,
     .      irtag(ir),rp(id),zp(id)
c          WRITE(fp,'(1X,I4,I5,2F9.5,1X,2F10.4,1P,2E12.4,1X,E10.2,0P'//
c     .             ' ,F6.2,1P,E12.4,0P,1X,A,2F10.5)')
c     .      ir,model,psitarg(ir,1),rho(ir,CELL1),
c     .      kteds(id),ktids(id),knds(id),jsat,kvds(id),mach,p,
c     .      irtag(ir),rp(id),zp(id)
        ENDDO
      ENDDO
c
c
c
      CALL HD(fp,'Target heat flux','SOLANAL-TARGETS-HEAT',5,67)
      sum_total = 0.0
      DO i1 = 2, 1, -1
        IF (i1.EQ.2) target = IKLO
        IF (i1.EQ.1) target = IKHI
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,A5,2A9,2X,1A10,2A7,2X,A7,2A10)')
     .    'ir','sol','psin','rho','jsat','Te','Ti','gamma',
     .    'P Flux','H Flux'
        sum_target = 0.0
        DO ring = 1, nrs-irsep+1
          IF (ring.LT.nrs-irtrap+1) THEN
            ir = irtrap + ring
          ELSE
            ir = ring - (nrs - irtrap + 1) + irsep 
          ENDIF
          IF (idring(ir).EQ.-1) CYCLE
          IF (i1.EQ.1) THEN
            model = model2(ir)
          ELSE
            model = model1(ir)
          ENDIF
          id = idds(ir,i1)
          jsat  = ABS(GetJsat(kteds(id),ktids(id),knds(id),kvds(id)))
          isat  = ABS(GetFlux(target,ir))
          gamma = GetGamma   (target,ir)
          qpara = GetHeatFlux(target,ir)
          WRITE(fp,'(1X,I4,I5,2F9.5,2X,1P,E10.2,0P,2F7.2,2X,'//
     .             'F7.2,1P,2E10.2,0P,F7.2,2X,F10.5,2X,3F6.2,2X,A)')
     .      ir,model,psitarg(ir,1),rho(ir,CELL1),
     .      jsat,kteds(id),ktids(id),gamma,isat,qpara,
     .      qpara/(2.0*PI*rp(id)*dds2(id)*1.0E+6),dds2(id),
     .      dds2(id)             /(rho(ir,OUT23)-rho(ir,IN14)),
     .      (dds2(id)*costet(id))/(rho(ir,OUT23)-rho(ir,IN14)),
     .      1.0/costet(id),
     .      irtag(ir)


          sum_target = sum_target + qpara
        ENDDO
        sum_total = sum_total + sum_target
        WRITE(fp,*) 
        WRITE(fp,*) 'TARGET TOTAL = ',sum_target/1.0E+6,' MW'
      ENDDO
      WRITE(fp,*) 
      WRITE(fp,*) 'TOTAL = ',sum_total/1.0E+6,' MW'
c
c     High index "symmetry" point parameters:
c
c
      CALL HD(fp,'High index symmetry point','SOLANAL-HISYM',5,67)
      DO i1 = 2, 1, -1
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,A5,1X,A10,2A10,A12,1X,2A12,1X,A10)')
     .    'ir','sol','s','Te','Ti','ne','p (with v)','p (no v)','Mach'
        DO ir = irsep, nrs
          IF (idring(ir).EQ.-1) CYCLE

          IF (i1.EQ.1) THEN
            model = model2(ir)
          ELSE
            model = model1(ir)
          ENDIF

          ik   = ikmids(ir) + 1
          p1   = CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                        ktibs(ik,ir),kvhs (ik,ir)/qt)
          p2   = CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                        ktibs(ik,ir),0.0         )
          mach = kvhs(ik,ir) / qt / 
     .           GetCs(ktebs(ik,ir),ktibs(ik,ir))

          WRITE(fp,'(1X,I4,I5,1X,F10.4,2F10.4,1P,E12.4,1X,2E12.4,0P,
     .               1X,F10.4,1X,A)')
     .      ir,model,kss(ik,ir),ktebs(ik,ir),ktibs(ik,ir),knbs(ik,ir),
     .      p1,p2,mach,irtag(ir)
        ENDDO
      ENDDO

c
c     FSP pressure error (by excluding the plasma velocity from the
c     calculation):
c...note: Presently dip_v is set to 0.0, so this doesn't show much.
c
      IF (cmodopt.EQ.1.AND.osm_matcht.NE.0) THEN
        WRITE(fp,*)
        WRITE(fp,'(1X,A4,2A12,2X,A12)') 'ir','n2T','mnv2','mnv2 / n2T'
        DO ir = irsep, irwall-1
          IF (idring(ir).EQ.-1) CYCLE

          p1 = CalcPressure(prp_ne(ir,osm_probe),prp_te(ir,osm_probe),
     .                      prp_ti(ir,osm_probe),0.0                 )
          p2 = CalcPressure(prp_ne(ir,osm_probe),0.0                 ,
     .                      0.0                 ,dip_v (ir,osm_probe))

          WRITE(fp,'(1X,I4,1P,2E12.4,0P,1X,F12.2,A)')
     .      ir,p1,p2,p2/p1*100.0,'% '//irtag(ir)
        ENDDO
      ENDIF
c
c
c
c
c
      CALL HD(fp,'Density peaks','SOLANAL-DENPEAK',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Density peaks:'
      WRITE(fp,'(2X,A2,2(1X,2(A4,A6,A10),2X,2A6,2X))')
     .  'ir','ik1','s1','ne1','ikm','sm','nem','te1','te1+1',
     .       'ik2','s2','ne2','ikm','sm','nem','te2','te2-1'
      DO ir = irsep, irwall-1
        IF (idring(ir).EQ.-1.OR.ikbound(ir,IKLO).EQ.0.OR.
     .                          ikbound(ir,IKHI).EQ.0) CYCLE

        ik1 = 1
        ik2 = ikmids(ir)
        maxik1 = 1
        DO ik = ik1, ik2
          IF (knbs(ik,ir).GT.knbs(maxik1,ir)) maxik1 = ik
        ENDDO

        ik1 = ik2 + 1
        ik2 = nks(ir)
        maxik2 = nks(ir)
        DO ik = ik1, ik2
          IF (knbs(ik,ir).GT.knbs(maxik2,ir)) maxik2 = ik
        ENDDO

        ik1 = ikbound(ir,IKLO)
        ik2 = ikbound(ir,IKHI)

        note1 = ' '
        note2 = ' '
        IF (ik1.NE.maxik1) note1 = '*'
        IF (ik2.NE.maxik2) note2 = '*'

        maxs = ksmaxs(ir)

        WRITE(fp,'(2X,I2,2(1X,2(I4,F6.3,1P,E10.2,0P),'//
     .           ' A,1X,2F6.2,2X))')
     .    ir,   ik1,      kss(   ik1,ir) /maxs,knbs(   ik1,ir),
     .       maxik1,      kss(maxik1,ir) /maxs,knbs(maxik1,ir),note1,
     .       ktebs(ik1,ir),ktebs(ik1+1,ir),
     .          ik2,(maxs-kss(   ik2,ir))/maxs,knbs(   ik2,ir),
     .       maxik2,(maxs-kss(maxik2,ir))/maxs,knbs(maxik2,ir),note2,
     .       ktebs(ik2,ir),ktebs(ik2-1,ir)
      ENDDO
c
c
c     Conservation:
c
c
c     Tally sources and fluxes:
c
      CALL HD(fp,'Continuity','SOLANAL-CONTINUITY',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Continuity:'
      WRITE(fp,'(1X,2A4,A5,4A12,2X,A15)')
     .  'ik1','ik2','ir','targ','ion','rec','cfp','net'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1.OR.
     .      .NOT.(model1(ir).EQ.22.OR.model1(ir).EQ.24).OR.
     .      .NOT.(model2(ir).EQ.22.OR.model2(ir).EQ.24)) CYCLE

        ik1 = ikbound(ir,IKLO)
        ik2 = ikbound(ir,IKHI)

        CALL CalcIntegral3(pinion,ik1,ik2,ir,rdum(1),2)
        CALL CalcIntegral3(pinrec,ik1,ik2,ir,rdum(2),2)
        CALL CalcIntegral3(osmcfp,ik1,ik2,ir,rdum(3),2)

        parflx = 0.0
        IF (ik1.EQ.1) THEN
          parflx = parflx + knds(idds(ir,2)) * kvds(idds(ir,2))
        ELSE
          parflx = parflx + knbs(ik1,ir)     * kvhs(ik1,ir)
        ENDIF
        IF (ik2.EQ.nks(ir)) THEN
          parflx = parflx - knds(idds(ir,1)) * kvds(idds(ir,1))
        ELSE
          parflx = parflx - knbs(ik2,ir)     * kvhs(ik2,ir)
        ENDIF

        partot = parflx + rdum(1) - rdum(2) + rdum(3)

        WRITE(fp,'(1X,2I4,I5,1P,4E12.3,2X,E15.5,0P)')
     .    ik1,ik2,ir,parflx,rdum(1),-rdum(2),rdum(3),partot/parflx
      ENDDO
c
c     Source accounting:
c
c
c
c
      CALL RZero(rdum,20)

c     Target flux in the SOL:
      DO ir = irsep, irwall-1
        rdum(1) = rdum(1) - GetFlux(IKLO,ir)
        rdum(2) = rdum(2) + GetFlux(IKHI,ir)
      ENDDO

c     Target flux in the PFZ:
      DO ir = irtrap+1, nrs
        rdum(3) = rdum(3) - GetFlux(IKLO,ir) + GetFlux(IKHI,ir)
      ENDDO

c     Target flux to the inner wall:
      IF (nbr.GT.0) THEN
        DO ir = irbreak, irwall-1
          rdum(5) = rdum(5) - GetFlux(IKLO,ir)
        ENDDO
      ENDIF
c
c     Recombination and ionisation in the SOL:
      CALL VolInteg(pinion,IKLO,irsep,irwall-1,rdum(6))
      CALL VolInteg(pinrec,IKLO,irsep,irwall-1,rdum(7))
      CALL VolInteg(pinion,IKHI,irsep,irwall-1,rdum(8))
      CALL VolInteg(pinrec,IKHI,irsep,irwall-1,rdum(9))
c
c     Recombination and ionisation in the PFZ:
      CALL VolInteg(pinion,3,irtrap+1,nrs,rdum(10))
      CALL VolInteg(pinrec,3,irtrap+1,nrs,rdum(11))
c
c     Recombination and ionisation in the core:
      CALL VolInteg(pinion,3,2,irsep-1,rdum(12))
      CALL VolInteg(pinrec,3,2,irsep-1,rdum(13))

c...Total ionisation:
      rdum(17) = rdum( 6) + rdum( 8) + rdum(10) + rdum(12)
c...Total recombination:
      rdum(14) = rdum( 7) + rdum( 9) + rdum(11) + rdum(13)
      rdum(15) = rdum( 1) + rdum( 2) + rdum( 3)
      rdum(16) = rdum(15) + rdum(14)



      ir    = 2
      area  = 0.0
      area2 = 0.0      
      DO ik = 1, nks(ir)-1
        id = korpg(ik,ir)
        r1 = 0.5 * (rvertp(1,id) + rvertp(4,id))
c
C       IPP/01 - Krieger: incredibly enough, the FUJI f90 compiler
C       chokes over v**2.0 if v is negative :-(
C       fixed by v**2.0 -> v**2
c
        delta = SQRT((rvertp(1,id) - rvertp(4,id))**2 +
     .               (zvertp(1,id) - zvertp(4,id))**2)      
c
c        delta = SQRT((rvertp(1,id) - rvertp(4,id))**2.0 +
c     .               (zvertp(1,id) - zvertp(4,id))**2.0)      
c
        area  = area  + 2.0 * PI * r1  * delta
        area2 = area2 + 2.0 * PI * rxp * delta
      ENDDO 
      WRITE(fp,*)
      WRITE(fp,*) 'AREA OF CORE RING (m2) = ',area
      WRITE(fp,*) 'AREA OF CORE RING (m2) = ',area2
      
c...balance
      CALL HD(fp,'Full plasma analysis','SOLANAL-FULLPLASMA',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Full plasma analysis:'
      WRITE(fp, 9) 'GLOBAL','LOW SOL','HIGH SOL','CORE','PFZ','IW'
9     FORMAT(3X,30X,6A11)
      WRITE(fp,10) 'Continuity      ion/(flx+rec)',
     .  rdum(17)/rdum(16),rdum( 6)/(rdum( 1)+rdum( 7)),
     .                    rdum( 8)/(rdum( 2)+rdum( 9))
      WRITE(fp,14) 'Ionisation           ion/tot ',
     .  rdum( 6)/(rdum(17)+EPS10),rdum( 8)/(rdum(17)+EPS10),
     .  rdum(12)/(rdum(17)+EPS10),rdum(10)/(rdum(17)+EPS10)
      WRITE(fp,10) 'Recombination   rec/(flx+rec)',
     .  rdum(14)/(rdum(16)+EPS10),rdum( 7)/(rdum(16)+EPS10),
     .  rdum( 9)/(rdum(16)+EPS10),rdum(13)/(rdum(16)+EPS10),
     .  rdum(11)/(rdum(16)+EPS10)
      WRITE(fp,16) 'Target flux     flx/(flx+rec)',
     .  rdum(15)/(rdum(16)+EPS10),rdum( 1)/(rdum(16)+EPS10),
     .  rdum( 2)/(rdum(16)+EPS10),rdum( 3)/(rdum(16)+EPS10),
     .  rdum( 5)/(rdum(16)+EPS10)
      WRITE(fp,10) 'Source ratio         rec/flx ',
     .  rdum(14)/(rdum(15)+EPS10)
10    FORMAT(3X,A30,5(F11.3:))
14    FORMAT(3X,A30,11X,5(F11.3:))
16    FORMAT(3X,A30,3F11.3,11X,2(F11.3:))

      WRITE(fp,*)
      WRITE(fp,11) 'Ionisation                   ',
     .  rdum(17),rdum( 6),rdum( 8),rdum(12),rdum(10)
      WRITE(fp,11) 'Recombination                ',
     .  rdum(14),rdum( 7),rdum( 9),rdum(13),rdum(11)
      WRITE(fp,15) 'Target flux                  ',
     .  rdum(15),rdum( 1),rdum( 2),rdum( 3),rdum( 5)
11    FORMAT(3X,A30,1P,5(E11.3:))
15    FORMAT(3X,A30,1P,3E11.3,11X,2E11.3)
c
c
c
c
c

      CALL RZero(rdum,20)
c
c     Target flux:
      DO ir = irsep, irwall-1
        IF (model1(ir).EQ.22) rdum(1) = rdum(1) - GetFlux(IKLO,ir)
        IF (model2(ir).EQ.22) rdum(2) = rdum(2) + GetFlux(IKHI,ir)
      ENDDO
c
c     Recombination and ionisation in the SOL:
      DO ir = irsep, irwall-1
        ikm = ikmids(ir)
        IF (model1(ir).EQ.22) THEN
          DO ik = 1, ikm
            rdum(3) = rdum(3) + pinion(ik,ir) * kvols(ik,ir)
            rdum(4) = rdum(4) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
        IF (model2(ir).EQ.22) THEN
          DO ik = ikm+1, nks(ir)
            rdum(5) = rdum(5) + pinion(ik,ir) * kvols(ik,ir)
            rdum(6) = rdum(6) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
      ENDDO

      rdum(7) = rdum(1) + rdum(4)
      rdum(8) = rdum(2) + rdum(6)

      CALL HD(fp,'Attached SOL plasma analysis','SOLANAL-ATTPLASMA',
     .        5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Attached SOL plasma analysis:'
      WRITE(fp, 9) 'LOW SOL','HIGH SOL'
      WRITE(fp,10) 'Continuity      ion/(flx+rec)',
     .  rdum(3)/(rdum(7)+EPS10),rdum(5)/(rdum(8)+EPS10)
      WRITE(fp,10) 'Source ratio         rec/flx ',
     .  rdum(4)/(rdum(1)+EPS10),rdum(6)/(rdum(2)+EPS10)
      WRITE(fp,11) 'Ionistaion                   ',rdum(3),rdum(5)
      WRITE(fp,11) 'Recombination                ',rdum(4),rdum(6)
      WRITE(fp,11) 'Target flux                  ',rdum(1),rdum(2)

c...temp
      IF (fp.EQ.PINOUT)
     .  WRITE(79,'(A,1X,3I4,1X,I4,1P,8E12.4,0P)') '''SRCANAATT 1.00''',
     .    rel_step,rel_iter,rel_count,ir,(rdum(i1),i1=1,8)
c
c
c
c
c
      CALL RZero(rdum,20)
c
c     Target flux:
      DO ir = irsep, irwall-1
        IF (model1(ir).EQ.24) rdum(1) = rdum(1) - GetFlux(IKLO,ir)
        IF (model2(ir).EQ.24) rdum(2) = rdum(2) + GetFlux(IKHI,ir)
      ENDDO
c
c     Recombination and ionisation in the SOL:
      DO ir = irsep, irwall-1
        ikm = ikmids(ir)
        IF (model1(ir).EQ.24) THEN
          DO ik = 1, ikm
            rdum(3) = rdum(3) + pinion(ik,ir) * kvols(ik,ir)
            rdum(4) = rdum(4) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
        IF (model2(ir).EQ.24) THEN
          DO ik = ikm+1, nks(ir)
            rdum(5) = rdum(5) + pinion(ik,ir) * kvols(ik,ir)
            rdum(6) = rdum(6) + pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDIF
      ENDDO

      rdum(7) = rdum(1) + rdum(4)
      rdum(8) = rdum(2) + rdum(6)

      CALL HD(fp,'Detached (prescription) SOL plasma analysis',
     .           'SOLANAL-DETPLASMA',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Detached (prescription) SOL plasma analysis:'
      WRITE(fp, 9) 'LOW SOL','HIGH SOL'
      WRITE(fp,10) 'Continuity      ion/(flx+rec)',
     .  rdum(3)/(rdum(7)+EPS10),rdum(5)/(rdum(8)+EPS10)
      WRITE(fp,10) 'Source ratio         rec/flx ',
     .  rdum(4)/(rdum(1)+EPS10),rdum(6)/(rdum(2)+EPS10)
      WRITE(fp,11) 'Ionistaion                   ',rdum(3),rdum(5)
      WRITE(fp,11) 'Recombination                ',rdum(4),rdum(6)
      WRITE(fp,11) 'Target flux                  ',rdum(1),rdum(2)

c...temp
      IF (fp.EQ.PINOUT)
     .  WRITE(79,'(A,1X,3I4,1X,I4,1P,8E12.4,0P)') '''SRCANADET 1.00''',
     .    rel_step,rel_iter,rel_count,ir,(rdum(i1),i1=1,8)
c
c
c
c
c
c
      IF (switch(SWION  ).EQ.1.0.OR.
     .    switch(SWION  ).EQ.2.0) THEN


        CALL HD(fp,'Low index ionisation distribution by stratum',
     .             'SOLANAL-LIION',5,67)
c        WRITE(fp,*)
c        WRITE(fp,*) 'Low index ionisation distribution by stratum:'
        WRITE(fp,'(3X,A4,4A10)') 'ir','1','2','3','4'

        DO ir = irsep, nrs
          IF (model1(ir).NE.22.AND.model1(ir).NE.24) CYCLE

          ik1 = 1
          ikm = ikmids (ir)

          CALL CalcIntegral3(pindata(1,1,H_ION1),ik1,ikm,ir,rdum(1),2)
          CALL CalcIntegral3(pindata(1,1,H_ION2),ik1,ikm,ir,rdum(2),2)
          CALL CalcIntegral3(pindata(1,1,H_ION3),ik1,ikm,ir,rdum(3),2)
          CALL CalcIntegral3(pindata(1,1,H_ION4),ik1,ikm,ir,rdum(4),2)
          CALL CalcIntegral3(pinion             ,ik1,ikm,ir,rdum(5),2)

c...temp
          CALL CalcIntegral3(pinion,ikbound(ir,IKLO),ikm,ir,rdum(6),2)
          CALL CalcIntegral3(pinrec,ikbound(ir,IKLO),ikm,ir,rdum(7),2)
          rdum(8) = knbs(ikbound(ir,IKLO),ir) *
     .              kvhs(ikbound(ir,IKLO),ir)
          rdum(9) = rdum(6) / (rdum(7) - rdum(8))

          WRITE(fp,'(3X,I4,4F10.4,2X,1P,5E10.2,0P,F10.2)')
     .      ir,rdum(1)/(rdum(5)+EPS10),rdum(2)/(rdum(5)+EPS10),
     .         rdum(3)/(rdum(5)+EPS10),rdum(4)/(rdum(5)+EPS10),
     .        (rdum(i1),i1=1,5),rdum(9)
        ENDDO

        CALL HD(fp,'High index ionisation distribution by stratum',
     .             'SOLANAL-HIION',5,67)
c        WRITE(fp,*)
c        WRITE(fp,*) 'High index ionisation distribution by stratum:'
        WRITE(fp,'(3X,A4,4A10)') 'ir','1','2','3','4'
        DO ir = irsep, nrs
          IF (model1(ir).NE.22.AND.model1(ir).NE.24) CYCLE

          ik1 = nks   (ir)
          ikm = ikmids(ir) + 1

          CALL CalcIntegral3(pindata(1,1,H_ION1),ikm,ik1,ir,rdum(1),2)
          CALL CalcIntegral3(pindata(1,1,H_ION2),ikm,ik1,ir,rdum(2),2)
          CALL CalcIntegral3(pindata(1,1,H_ION3),ikm,ik1,ir,rdum(3),2)
          CALL CalcIntegral3(pindata(1,1,H_ION4),ikm,ik1,ir,rdum(4),2)
          CALL CalcIntegral3(pinion             ,ikm,ik1,ir,rdum(5),2)

c...temp
          CALL CalcIntegral3(pinion,ikm,ikbound(ir,IKHI),ir,rdum(6),2)
          CALL CalcIntegral3(pinrec,ikm,ikbound(ir,IKHI),ir,rdum(7),2)
          rdum(8) = knbs(ikbound(ir,IKHI),ir) *
     .              kvhs(ikbound(ir,IKHI),ir)
          rdum(9) = rdum(6) / (rdum(7) + rdum(8))

          WRITE(fp,'(3X,I4,4F10.4,2X,1P,5E10.2,0P,F10.2)')
     .      ir,rdum(1)/(rdum(5)+EPS10),rdum(2)/(rdum(5)+EPS10),
     .         rdum(3)/(rdum(5)+EPS10),rdum(4)/(rdum(5)+EPS10),
     .        (rdum(i1),i1=1,5),rdum(9)
        ENDDO
      ENDIF
c
c
c
c
c
      CALL HD(fp,'Electron power','SOLANAL-ELECPOWER',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Electron power:'
      WRITE(fp,'(1X,2A4,A5,4A12,2X,A15)')
     .  'ik1','ik2','ir','targ','Qe','-Pei','cfe','net'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        IF ((model1(ir).EQ.22.OR.model1(ir).EQ.24).AND.
     .      (model2(ir).EQ.22.OR.model2(ir).EQ.24)) THEN

          ik1 = ikbound(ir,IKLO)
          ik2 = ikbound(ir,IKHI)

c...prad!
          IF (switch(SWPHELP).GE.2.0) THEN
            CALL CalcIntegral2(pinqe(1,ir),ik1,ik2,ir,ddum1,2)
            CALL CalcIntegral2(osmqe(1,ir),ik1,ik2,ir,ddum2,2)
            intg(IN_PQE) = ddum1 + ddum2
          ENDIF
          IF (switch(SWPEI  ).EQ.4.0.OR.
     .        switch(SWPEI  ).EQ.5.0)
     .      CALL CalcIntegral2(osmpei(1,ir),ik1,ik2,ir,intg(IN_PEI),2)
          IF (switch(SWPOW  ).EQ.12.0.OR.
     .        switch(SWPOW  ).EQ.13.0)
     .       CALL CalcIntegral2(osmcfe(1,ir),ik1,ik2,ir,intg(IN_CFE),2)

          eleflx = elecptarg(ir,3)
          eletot = eleflx + intg(IN_PQE) - intg(IN_PEI) + intg(IN_CFE)

          WRITE(fp,'(1X,2I4,I5,1P,4E12.3,2X,E15.5,0P)')
     .      ik1,ik2,ir,eleflx,intg(IN_PQE),
     .     -intg(IN_PEI),intg(IN_CFE),eletot/(eleflx+EPS10)
        ENDIF
      ENDDO

      CALL HD(fp,'Ion power','SOLANAL-IONPOWER',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Ion power:'
      WRITE(fp,'(1X,2A4,A5,4A12,2X,A15)')
     .  'ik1','ik2','ir','targ','Qi','Pei','cfi','net'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        IF ((model1(ir).EQ.22.OR.model1(ir).EQ.24).AND.
     .      (model2(ir).EQ.22.OR.model2(ir).EQ.24)) THEN

          ik1 = ikbound(ir,IKLO)
          ik2 = ikbound(ir,IKHI)


          IF (switch(SWPCX  ).GE.2.0)
     .      CALL CalcIntegral2(pinqi(1,ir),ik1,ik2,ir,intg(IN_PQI),2)
          IF (switch(SWPEI  ).EQ.4.0.OR.
     .        switch(SWPEI  ).EQ.5.0)
     .      CALL CalcIntegral2(osmpei(1,ir),ik1,ik2,ir,intg(IN_PEI),2)
          IF (switch(SWPOW  ).EQ.12.0.OR.
     .        switch(SWPOW  ).EQ.13.0)
     .      CALL CalcIntegral2(osmcfi(1,ir),ik1,ik2,ir,intg(IN_CFI),2)


          ionflx = ionptarg(ir,3)
          iontot = ionflx + intg(IN_PQI) + intg(IN_PEI) + intg(IN_CFI)

          WRITE(fp,'(1X,2I4,I5,1P,4E12.3,2X,E15.5,0P)')
     .      ik1,ik2,ir,ionflx,intg(IN_PQI),
     .      intg(IN_PEI),intg(IN_CFI),iontot/(ionflx+EPS10)
        ENDIF
      ENDDO


      CALL CalcRadiatedPower(radpow,1)

c REGION:
c
c 1 - IKLO
c 2 - IKHI
c 3 - ENTIRE RING
c 4 - INNER HALF-RING, BELOW X-POINT
c 5 - INNER HALF-RING, ABOVE X-POINT
c 6 - OUTER HALF-RING, BELOW X-POINT
c 7 - OUTER HALF-RING, ABOVE X-POINT
c

      CALL HD(fp,'EIRENE P_rad_h VOLUME INTEGRATED POWER',
     .           'SOLANA-HRADPOW-HD',5,67)

      CALL RZero(rdum,20)

c...  SOL:
      CALL VolInteg(radpow,4,irsep,irwall-1,rdum(1))
      CALL VolInteg(radpow,5,irsep,irwall-1,rdum(2))
      CALL VolInteg(radpow,6,irsep,irwall-1,rdum(3))
      CALL VolInteg(radpow,7,irsep,irwall-1,rdum(4))

c...  PFZ:
      CALL VolInteg(radpow,IKLO,irtrap+1,nrs,rdum(5))
      CALL VolInteg(radpow,IKHI,irtrap+1,nrs,rdum(6))

c...  Core:
      CALL VolInteg(radpow,IKLO,2,irsep-1,rdum(7))
      CALL VolInteg(radpow,IKHI,2,irsep-1,rdum(8))

      CALL VolInteg(radpow,3,2,nrs,rdum(9))

c...  Convert to mega-Watts:
      DO i1 = 1, 9
        rdum(i1) = rdum(i1) * 1.0E-06
      ENDDO

      fact = 2.0 * PI * rxp

      WRITE(fp,*) 
      WRITE(fp,'(4X,A,F8.3,A)') 'TOROIDAL CIR.= ',fact,' m'

      WRITE(fp,*) 
      WRITE(fp,20) 'REGION','INTEGRAL/m','INTEGRAL'
      WRITE(fp,20) ' ','(MW/m)','(MW)'

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index SOL below X-point',rdum(1),rdum(1)*fact
      WRITE(fp,21) 'Low  index SOL above X-point',rdum(2),rdum(2)*fact
      WRITE(fp,21) 'High index SOL below X-point',rdum(3),rdum(3)*fact
      WRITE(fp,21) 'High index SOL above X-point',rdum(4),rdum(4)*fact

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index PFZ',rdum(5),rdum(5)*fact
      WRITE(fp,21) 'High index PFZ',rdum(6),rdum(6)*fact

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index core',rdum(7),rdum(7)*fact
      WRITE(fp,21) 'High index core',rdum(8),rdum(8)*fact

      WRITE(fp,*) 
      WRITE(fp,21) 'Low  index SOL', rdum(1)+rdum(2),
     .                              (rdum(1)+rdum(2))*fact
      WRITE(fp,21) 'High index SOL', rdum(3)+rdum(4),
     .                              (rdum(3)+rdum(4))*fact
      WRITE(fp,21) 'SOL', rdum(1)+rdum(2)+rdum(3)+rdum(4),
     .                   (rdum(1)+rdum(2)+rdum(3)+rdum(4))*fact
      WRITE(fp,21) 'PFZ', rdum(5)+rdum(6),
     .                   (rdum(5)+rdum(6))*fact
      WRITE(fp,21) 'Core', rdum(7)+rdum(8),
     .                    (rdum(7)+rdum(8))*fact
      WRITE(fp,21) 'Entire plasma',rdum(9),rdum(9)*fact


20    FORMAT(4X,A29,2A12  )
21    FORMAT(4X,A29,2F12.3)








c
c     Some momentum loss analysis:
c
      CALL HD(fp,'Momentum loss (SOL22,24 only)','SOLANAL-MOMLOSS',5,67)
c      WRITE(fp,*)
c      WRITE(fp,*) 'Momentum loss (SOL22,24 only):'
      WRITE(fp,'(1X,2A4,A5,3A12,2X,3A12)')
     .  'ik1','ik2','ir','total rec','total mloss','model mloss',
     .  'total/rec','model/rec','model/total'

      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        ik1 = ikbound(ir,IKLO)
        ik2 = ikmids(ir)

        IF (model1(ir).EQ.22.OR.model1(ir).EQ.24) THEN
          CALL CalcIntegral3(pinrec,1  ,ik2,ir,rdum(1),2)
          CALL CalcIntegral3(osmmp ,1  ,ik2,ir,rdum(2),2)
          CALL CalcIntegral3(osmmp ,ik1,ik2,ir,rdum(3),2)

          WRITE(fp,'(1X,2I4,I5,1P,3E12.4,2X,3E12.4,0P)')
     .      ik1,ik2,ir,rdum(1),rdum(2),rdum(3),
     .      rdum(2)/rdum(1),rdum(3)/rdum(1),rdum(3)/rdum(2)
        ENDIF
      ENDDO

      WRITE(fp,*)
      DO ir = irsep, nrs
        IF (idring(ir).EQ.-1) CYCLE

        ik1 = ikmids(ir) + 1
        ik2 = ikbound(ir,IKHI)

        IF (model2(ir).EQ.22.OR.model2(ir).EQ.24) THEN
          CALL CalcIntegral3(pinrec,1  ,ik2,ir,rdum(1),2)
          CALL CalcIntegral3(osmmp ,1  ,ik2,ir,rdum(2),2)
          CALL CalcIntegral3(osmmp ,ik1,ik2,ir,rdum(3),2)

          WRITE(fp,'(1X,2I4,I5,1P,3E12.4,2X,3E12.4,0P)')
     .      ik1,ik2,ir,rdum(1),rdum(2),rdum(3),
     .      rdum(2)/rdum(1),rdum(3)/rdum(1),rdum(3)/rdum(2)
        ENDIF
      ENDDO


c      STOP 'Temp in AS'




c
c
c     Density peak flow and recombination:
c
c
      DO ir = 1, nrs
        IF (idring(ir).EQ.-1.OR.model1(ir).NE.24) CYCLE




      ENDDO
      DO ir = 1, nrs
        IF (idring(ir).EQ.-1.OR.model2(ir).NE.24) CYCLE



      ENDDO


c
c     Symmetry point continuity:
c

      IF (.FALSE.) THEN

c...    Figure out the maximum ion flux into the divertor.  This is 
c       specific to grid SL2 for 990429019 at the moment:

c...    Assign surface across the divertor throat:

c       Ring 19 is the outermost ring for the inner SOL that is still in
c       the divertor:
        ir = 19
        ik = 1
        a1 = DBLE(rvertp(2,korpg(ik,ir)))
        a2 = DBLE(zvertp(2,korpg(ik,ir)))

c       Ring 19 is the outermost ring for the outer SOL that is still in
c       the divertor:
        ir = 19
        ik = nks(ir)
        b1 = DBLE(rvertp(3,korpg(ik,ir)))
        b2 = DBLE(zvertp(3,korpg(ik,ir)))

        flux   (IKLO) = 0.0
        flux   (IKHI) = 0.0
        fluxmax(IKLO) = 0.0
        fluxmax(IKHI) = 0.0
        DO ir = irsep, irwall-1
          DO ik = 1, osm_sympt(ir)
            id = korpg(ik,ir)
            c1 = DBLE(rvertp(1,id))
            c2 = DBLE(zvertp(1,id))
            d1 = DBLE(rvertp(4,id))
            d2 = DBLE(zvertp(4,id))
            IF (ChkInt(a1,a2,b1,b2,c1,c2,d1,d2)) THEN
              deltar = rvertp(4,id) - rvertp(1,id)
              deltaz = zvertp(4,id) - zvertp(1,id)
	      alpha  = ATAN2C(deltar,deltaz)
              deltar = SNGL(b1 - a1)
              deltaz = SNGL(b2 - a2)
              beta   = ATAN2C(deltar,deltaz) - alpha
              cost   = COS(PI / 2.0 - beta)
c...          Find segment of divertor boundary in cell IK,IR:
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              e1 = c1 + tcd * (d1 - c1)
              e2 = c2 + tcd * (d2 - c2)
              c1 = DBLE(rvertp(2,id))
              c2 = DBLE(zvertp(2,id))
              d1 = DBLE(rvertp(3,id))
              d2 = DBLE(zvertp(3,id))
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              f1 = c1 + tcd * (d1 - c1)
              f2 = c2 + tcd * (d2 - c2)
              area = SQRT(SNGL((f1 - e1)**2 + (f2 - e2)**2)) * 
     .               2.0 * PI * rs(ik,ir) * eirtorfrac
              flux(IKLO) = flux(IKLO) + knbs(ik,ir) * kvhs(ik,ir) * 
     .                                  cost * bratio(ik,ir) * area
              fluxmax(IKLO) = fluxmax(IKLO) + knbs(ik,ir) * 
     .           GetCs(ktibs(ik,ir),ktebs(ik,ir)) * cost * 
     .           bratio(ik,ir) * area
              WRITE(0,'(A,2I6,2F12.6,1P,2E10.2,0P,F7.3)') 
     .          '-->',ik,ir,cost,area,flux(IKLO),fluxmax(IKLO),
     .          kvhs(ik,ir)/GetCs(ktibs(ik,ir),ktebs(ik,ir))
            ENDIF
          ENDDO
c...      High index contribution:
          DO ik = osm_sympt(ir)+1, nks(ir)
            id = korpg(ik,ir)
            c1 = DBLE(rvertp(1,id))
            c2 = DBLE(zvertp(1,id))
            d1 = DBLE(rvertp(4,id))
            d2 = DBLE(zvertp(4,id))
            IF (ChkInt(a1,a2,b1,b2,c1,c2,d1,d2)) THEN
              deltar = rvertp(1,id) - rvertp(4,id)
              deltaz = zvertp(1,id) - zvertp(4,id)
	      alpha  = ATAN2C(deltar,deltaz)
              deltar = SNGL(b1 - a1)
              deltaz = SNGL(b2 - a2)
              beta   = ATAN2C(deltar,deltaz) - alpha
              cost   = COS(PI / 2.0 - beta)
c...          Find segment of divertor boundary in cell IK,IR:
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              e1 = c1 + tcd * (d1 - c1)
              e2 = c2 + tcd * (d2 - c2)
              c1 = DBLE(rvertp(2,id))
              c2 = DBLE(zvertp(2,id))
              d1 = DBLE(rvertp(3,id))
              d2 = DBLE(zvertp(3,id))
              CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)
              f1 = c1 + tcd * (d1 - c1)
              f2 = c2 + tcd * (d2 - c2)
              area = SQRT(SNGL((f1 - e1)**2 + (f2 - e2)**2)) * 
     .               2.0 * PI * rs(ik,ir) * eirtorfrac
              flux(IKHI) = flux(IKHI) + knbs(ik,ir) * kvhs(ik,ir) * 
     .                                  cost * bratio(ik,ir) * area
              fluxmax(IKHI) = fluxmax(IKHI) + knbs(ik,ir) * 
     .           GetCs(ktibs(ik,ir),ktebs(ik,ir)) * cost * 
     .           bratio(ik,ir) * area
              WRITE(0,'(A,2I6,2F12.6,1P,2E10.2,0P,F7.3)') 
     .          '-->',ik,ir,cost,area,flux(IKHI),fluxmax(IKHI),
     .          kvhs(ik,ir)/GetCs(ktibs(ik,ir),ktebs(ik,ir))
            ENDIF
          ENDDO



        ENDDO

      ENDIF


c     Target flux by target region:
      WRITE(fp,*)
      WRITE(fp,*) 'Target flux by region:'

      rdum = 0.0
 
      DO i1 = 1, grdntreg(IKLO)
        rdum(1) = 0.0
        DO i2 = 1, grdntseg(i1,IKLO)
          ir = grdtseg(i2,i1,IKLO)
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          rdum(1) = rdum(1) - GetFlux(IKLO,ir)      
        ENDDO
        rdum(2) = rdum(2) + rdum(1)
        WRITE(fp,'(A,I4,1P,E11.3,0P,F11.1)') 
     .    '  IKLO:',i1,rdum(1),rdum(1)*ECH
      ENDDO

      DO i1 = 1, grdntreg(IKHI)
        rdum(1) = 0.0
        DO i2 = 1, grdntseg(i1,IKHI)
          ir = grdtseg(i2,i1,IKHI)
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          rdum(1) = rdum(1) + GetFlux(IKHI,ir)      
        ENDDO
        rdum(2) = rdum(2) + rdum(1)
        WRITE(fp,'(A,I4,1P,E11.3,0P,F11.1)') 
     .    '  IKHI:',i1,rdum(1),rdum(1)*ECH
      ENDDO
      WRITE(fp,'(A,4X,1P,E11.3,0P,F11.1)') 
     .  ' TOTAL:',rdum(2),rdum(2)*ECH

      CALL DB('Done')


      RETURN
99    STOP 'AnalyseSolution'
      END


