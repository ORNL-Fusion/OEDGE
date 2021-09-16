c ======================================================================
c
c block data program unit SLData
c
c ======================================================================
c
c      jdemod - intialization moved to mod_slcom_lim.f90
c
c     
c      BLOCK DATA SLData
c      use mod_params
c      use mod_slcom
c      
c      INCLUDE 'params'
c      INCLUDE 'slcom'
c     
c      DATA  YSBIN / -1.0E+06, -3.7696, -1.7560, -0.5338, 
c     +                0.3259,  1.1181,  1.9309,  2.7438,
c     +                3.2317, 1.0E+06 /
c      DATA  BSBIN /  -3.9985, -1.9847, -0.7625,  0.0971,
c     +                0.8893,  1.7021,  2.5149,  3.0029,
c     +                3.1231, 1.0E+06 /
c
c      DATA  YSBIN / -1.0E+06, -3.2317, -2.7438, -1.9309, 
c     +               -1.1181, -0.3259,  0.5338,  1.7560,
c     +                3.7696, 1.0E+06 /
c      DATA  BSBIN /  -3.1231, -3.0029, -2.5149, -1.7021,
c     +               -0.8893, -0.0971,  0.7625,  1.9847,
c     +                3.9985, 1.0E+06 /
c
c      DATA  NBIN /10/
c           
c      END BLOCK DATA
c
c ======================================================================
c
c subroutine OutputDiagnostics
c
c ======================================================================
c
      SUBROUTINE OutputDiag
      use mod_params
      use mod_dynam1
      use mod_dynam3
      use mod_comt2
      use mod_comtor
      use mod_comtau
      use mod_slcom
      implicit none

c      INCLUDE   'params'                                                         
c      INCLUDE   'dynam1'                                                        
c      INCLUDE   'dynam3'                                                        
c      INCLUDE   'comtau'                                                        
c      INCLUDE   'comtor'                                                        
c      INCLUDE   'comt2'                                                        
c      INCLUDE   'slcom'

      INTEGER II,SUM
c
c     jdemod - defined some variables
c
      integer iz,ix
      real count

      WRITE(63,*) ' '
      WRITE(63,'(2A)') 'CASE: ',TITL2
      WRITE(63,*) ' '

      WRITE(63,'(A)') ' '
      WRITE(63,'(A)') 'Inboard transport parameters:'
      WRITE(63,'(A)') ' '

      WRITE(63,'(A,G15.5,A)') 'Impurity inj temp   ',CTEMSC  ,' eV'
      WRITE(63,'(A,G15.5,A)') 'Electric field      ',CEYIN   ,' V/m'
      WRITE(63,'(A,G15.5,A)') 'Radial Dperp        ',CRDXI(1),' m2/s'
      WRITE(63,'(A,G15.5,A)') 'Parallel drift      ',CVHYIN  ,' m/s'
      WRITE(63,'(A,G15.5,A)') 'Poloidal Dperp      ',CDPOL   ,' m2/s'
      WRITE(63,'(A,G15.5,A)') 'Poloidal drift      ',CVPOL   ,' m/s'

      WRITE(63,*) ' '      
      WRITE(63,'(A,I6)') 'Maximum ionisation state ',MIZS

      WRITE(63,*) ' '      
c
c     jdemod - subscripts undefined
c
c      WRITE(63,'(A,3I4,A,G10.5,A,F10.3,A)') 'Ionisation rate peak (',
c     +  IXIRP,IYIRP,IPIRP,'): ',CRNBS(IXIRP,IYIRP),' m^-3, ',
c     +  CTEMBS(IXIRP,IYIRP),' eV'

      WRITE(63,*) ' '
      WRITE(63,*) 'Density integrated over all space:'
      WRITE(63,*) ' '          
    
      WRITE(63,'(3X,A,G12.3)') 'Total volume =',TOTVOL
      WRITE(63,*) ' ' 
      DO IZ = 1, MIZS
        WRITE(63,'(3X,I6,G12.3)') IZ,TOTDEN(IZ)
      ENDDO


      IF (optdp.EQ.1) THEN
      
        WRITE(63,*) ' '
        WRITE(63,'(A)') 'DIVIMP source distribution:'
        WRITE(63,*) ' '
        
        DO II = NBIN, 1, -1
          SUM = NSBIN(II,1) + NSBIN(II,2) + NSBIN(II,3) + 
     +          NSBIN(II,4) + NSBIN(II,5)

          WRITE(63,'(3x,2I6,G12.4,6I6)') 
     +      II, 33-II, BSBIN(II), NSBIN(II,1), 
     +      NSBIN(II,2), NSBIN(II,3), NSBIN(II,4), NSBIN(II,5),
     +      SUM
        ENDDO

        WRITE(63,*) ' '
        WRITE(63,*) '   Integrated density over all space:'
        WRITE(63,*) ' '          
      
        WRITE(63,'(3X,A,G12.3)') 'Total volume =',TOTVOL
        WRITE(63,*) ' '
        DO IZ = 1, MIZS
          WRITE(63,'(3X,I6,G12.3)') IZ,TOTDEN(IZ)
        ENDDO
        
        WRITE(63,*) ' '

        COUNT = 0
        DO IX = 1, MIMPS
          IF (TAG2(IX).EQ.1.0) COUNT = COUNT + 1.0
        ENDDO
        
        WRITE(63,'(3X,A,I7)') 
     +    'ST 1: Number of particles                              ',
     +    MIMPS     
        WRITE(63,'(3X,A,I7)') 
     +    'ST 2: Number of ions                                   ',
     +    MATIZ     
        WRITE(63,'(3X,A,G12.5)') 
     +    'ST 3: Radial location of DIVIMP grid                   ',
     +    MARK
        WRITE(63,'(3X,A,I7)') 
     +    'ST 4: Number of ions entering DIVIMP grid              ',
     +    INT(COUNT)
        WRITE(63,'(3X,A,F7.4)') 
     +    'ST 5: Fraction of particle source entering DIVIMP grid ',
     +    COUNT/REAL(MIMPS)
        WRITE(63,'(3X,A,I7,A,I7,A)') 
     +    'ST10: Number of ions ionised beyond limit (too soon)   ',
     +    INT(IZLOSS),' (',INT(TSLOSS),')'
        WRITE(63,'(3X,A,I7)') 
     +    'ST 6: Number of ions striking limiter                  ',
     +    INT(LLOSS)
        WRITE(63,'(3X,A,I7)') 
     +    'ST 8: Number of ions striking wall                     ',
     +    INT(WLOSS)
        WRITE(63,'(3X,A,I7)') 
     +    'ST 9: Number of ions striking target                   ',
     +    INT(TGLOSS)
        WRITE(63,'(3X,A,G12.4)') 
     +    'ST 7: Average loss time to limiter (for DIVIMP grid)   ',
     +    ALOSS/(LLOSS + WLOSS + TGLOSS)
      
      ENDIF
        
      RETURN
      END
c
c ======================================================================
c
c subroutine GetProfiles
c
c ======================================================================
c
      SUBROUTINE GetProfiles
c
c     jdemod - add implicit none
c
      use mod_params
      use mod_dynam1
      use mod_dynam3
      use mod_comt2
      use mod_cneut
      use mod_cnoco
      use mod_comtor
      use mod_comtau
      use mod_comxyt
      use mod_coords
      use mod_slcom
      implicit none

c      INCLUDE   'params'                                                         
c      INCLUDE   'dynam1'
c      INCLUDE   'dynam3'
c      INCLUDE   'comtor'                                                        
c      INCLUDE   'comtau'                                                        
c      INCLUDE   'comt2'                                                         
c      INCLUDE   'coords'                                                        
c      INCLUDE   'comxyt'                                                        
c      INCLUDE   'cneut'                                                         
c      INCLUDE   'cnoco'                                                         
c      INCLUDE   'slcom'


      REAL       PEAK1,PEAK2,PEAK3,PEAK4,PEAK5,PEAK6,PEAK7,PEAK8,PEAK9

      REAL       PEAK


      REAL       TMPVOL,NUMSUM,VOLSUM,YPROL,YPROR,YLENL,YLENR
      INTEGER    IXIRP,IYIRP,IPIRP,FLAG
      REAL       YPRO(-MAXNYS:MAXNYS,0:MAXIZS),
     +           PPRO(-MAXNPS:MAXNPS,0:MAXIZS),
     +           EPRO (MAXNXS,-MAXNPS:MAXNPS),
     +           EAREA(MAXNXS,-MAXNPS:MAXNPS),
     +           EVOL (MAXNXS,-MAXNPS:MAXNPS)
      REAL       NPRO(-MAXNYS:MAXNYS,-MAXNPS:MAXNPS)
      REAL       IPRO(MAXNXS,0:6),Cs
      REAL       COUNT

      INTEGER    MAXIX,MAXIP
      REAL       MAXE
c
c     jdemod - defined missing variables
c
      integer :: iy,iz,ip,ix,i,iqx

c
c Initialize arrays:
c
      DO IY = -NYS,NYS 
        DO IZ = 0, MIZS
          YPRO(IY,IZ)    = 0.0
        ENDDO
        NPRO(IY,1)  = 0.0

c
c       jdemod - fixed array bounds addressing error when iy=0 
c              - also check for ywids(iy) = 0
c
        if (iy.ne.0) then 
           if (ywids(abs(iy)).ne.0.0) then 
              INJBINT(IY) = INJBINT(IY) / YWIDS(ABS(IY))
           else
              INJBINT(IY) = 0.0
           endif
        else
           if (Ywids(1).ne.0.0) then 
              INJBINT(IY) = INJBINT(IY) / YWIDS(1)
           else
              INJBINT(IY) = 0.0
           endif
        endif
      ENDDO    
      
      DO IP = -MAXNPS,MAXNPS
        DO IZ = 0, MIZS
          PPRO(IP,IZ) = 0.0
        ENDDO
      ENDDO
c
c Integrate profiles:
c
c
c Toroidal (parallel to field line) profile:
c
      DO IZ = 0,MIZS

        PEAK = 0

        DO IY = -NYS, NYS
          IF (IY.NE.0) THEN
            NUMSUM = 0.0
            VOLSUM = 0.0
      
            DO IX = 1, NXS
              DO IP = -MAXNPS, MAXNPS 
      
                 ! jdemod - YWIDS only defined for IY>0
                 !          The width of IY=1 is from 0.0 to YS(1) - thus any array
                 !          ranging from -NYS to +NYS should contain no data in the IY=0 element
                 if (iy.ne.0) then 
                  TMPVOL = XWIDS(IX) * YWIDS(ABS(IY))
                 else
                  TMPVOL = 0.0
                 endif

               NUMSUM = NUMSUM + DDLIM3(IX,IY,IZ,IP) * TMPVOL
               VOLSUM = VOLSUM + TMPVOL
      
              ENDDO 
            ENDDO

            YPRO(IY,IZ) = NUMSUM / VOLSUM
          ENDIF
        ENDDO   
      ENDDO
c
c Find integral within viewing are for either side of source:
c
      WRITE(63,*) ' '
      WRITE(63,*) 'Total integral for +ve and -ve y (iz,-ve,+ve):'
      WRITE(63,*) ' '

      WRITE( 7,*) ' '
      WRITE( 7,*) 'PARALLEL PLUME ASYMMETRY'
      WRITE( 7,*) ' '

      DO IZ = 0,MIZS
        YPROL = 0.0
        YPROR = 0.0
        YLENL = 0.0
        YLENR = 0.0
        DO IY = -0.5*NYS,0.5*NYS
          ! jdemod - added test to avoid iy = 0 
          if (iy.ne.0) then 
            IF (YS(ABS(IY)).LE.0.16.AND.IY.LT.0) THEN
              YPROL = YPROL + (YPRO(IY,IZ)+YPRO(IY+NYS+1,IZ)) * 
     +                        YWIDS(ABS(IY))
              YLENL = YLENL + YWIDS(ABS(IY))
            ELSEIF (YS(ABS(IY)).LE.0.16.AND.IY.GT.0) THEN
              YPROR = YPROR + (YPRO(IY,IZ)+YPRO(IY-NYS-1,IZ)) * 
     +                        YWIDS(ABS(IY))
              YLENR = YLENR + YWIDS(ABS(IY))
            ENDIF
          endif

        ENDDO

        YPROL = YPROL / YLENL
        YPROR = YPROR / YLENR

        WRITE(63,'(A,I6,4G12.4)') 'tag1: ',
     +    IZ,YPROL,YPROR,YLENL,YLENR

        WRITE( 7,'(A,I6,4G12.4)') 'marker 1  ',
     +    IZ,YPROL,YPROR,YLENL,YLENR

      ENDDO
c
c Radial ionisation rate profile:
c      
      DO IX = 1, NXS
          NUMSUM = 0.0
          VOLSUM = 0.0

          DO IY = -NYS, NYS
            IF (IY.NE.0) THEN
              DO IP = -MAXNPS, MAXNPS 
                TMPVOL = XWIDS(IX) * YWIDS(ABS(IY))
                NUMSUM = NUMSUM + TIZ3(IX,IY,0,IP) * TMPVOL
                VOLSUM = VOLSUM + TMPVOL
              ENDDO 
            ENDIF
          ENDDO
            
          IPRO(IX,0) = NUMSUM / VOLSUM

          NUMSUM = 0.0
          VOLSUM = 0.0
      
          DO IY = -NYS, NYS
            IF (IY.NE.0) THEN
              DO IP = -MAXNPS, MAXNPS 
                TMPVOL = XWIDS(IX) * YWIDS(ABS(IY))
                NUMSUM = NUMSUM + TIZ3(IX,IY,2,IP) * TMPVOL
                VOLSUM = VOLSUM + TMPVOL
              ENDDO 
            ENDIF
          ENDDO
      
          IPRO(IX,1) = NUMSUM / VOLSUM
      
      ENDDO   
c
c Poloidal density profile:
c
      DO IZ = 0, MIZS
        DO IP = -MAXNPS, MAXNPS 
          NUMSUM = 0.0
          VOLSUM = 0.0
      
          DO IX = 1, NXS
            DO IY = -NYS, NYS
              IF (IY.NE.0) THEN
                TMPVOL = XWIDS(IX) * YWIDS(ABS(IY))
                NUMSUM = NUMSUM + DDLIM3(IX,IY,IZ,IP) * TMPVOL
                VOLSUM = VOLSUM + TMPVOL
              ENDIF
            ENDDO 
          ENDDO
      
          PPRO(IP,IZ) = NUMSUM / VOLSUM
        ENDDO   
      ENDDO
c
c Electron production rate integration:
c
      ! jdemod - these lines do not work for ix = maxnxs - change ix loop to maxnxs-1 from maxnxs
      DO IX = 1, MAXNXS-1
        DO IP = -MAXNPS+1, MAXNPS-1

c
c         jdemod - fix bug in indexing of epro - IY -> IP
c
c          EPRO(IX,IY) = 0.0
c
          EPRO(IX,IP) = 0.0
c
c Have to make sure this is correct if last charge state is the highest:
c      
          DO IZ = 0, MIZS

            NUMSUM = 0.0
            VOLSUM = 0.0

            DO IY = -NYS, NYS
              IF (IY.NE.0) THEN
                TMPVOL = XWIDS(IX) * YWIDS(ABS(IY))
                NUMSUM = NUMSUM + TIZ3(IX,IY,IZ,IP) * TMPVOL
                VOLSUM = VOLSUM + TMPVOL
              ENDIF
            ENDDO 

            EPRO(IX,IP) = EPRO(IX,IP) + NUMSUM / VOLSUM / 
     +                                  CFIZS(IX,1,IZ)
            EVOL(IX,IP) = VOLSUM

c          IF (VOLSUM.EQ.0.0.OR.CFIZS(IX,IY,IZ).EQ.0.0) THEN
c            WRITE(63,*) 'trouble: ',IX,IY,IZ,VOLSUM,CFIZS(IX,IY,IZ)
c            WRITE( 0,*) 'trouble: ',IX,IY,IZ,VOLSUM,CFIZS(IX,IY,IZ)
c          ENDIF

c          write(63,*) 'epro ',ix,iy,iz,numsum,volsum,
c     +      cfizs(ix,iy,iz),epro(ix,ip)

          ENDDO


          ! jdemod - these lines do not work for ix = maxnxs
          IF (IP.LE.0) THEN
            EAREA(IX,IP) = (XS(IX+1)-XS(IX)) * (PS(-IP)-PS(-IP-1))
          ELSEIF (IP.GT.0) THEN
            EAREA(IX,IP) = (XS(IX+1)-XS(IX)) * (PS( IP)-PS( IP-1))
          ENDIF
        ENDDO   
      ENDDO
c
c Find the peak ionisation rate:
c
        IXIRP = 0
        IYIRP = 0
        IPIRP = 0

        PEAK9 = 0.0
 
        DO IX = 1, NXS
          DO IY = -0.5*NYS, 0.5*NYS
            DO IP = -MAXNPS+1, MAXNPS-1
              IF (DDLIM3(IX,IY,0,IP)/CFIZS(IX,1,0).GT.PEAK9) THEN
                IXIRP = IX
                IYIRP = IY
                IPIRP = IP

                PEAK9 = DDLIM3(IX,IY,0,IP)/CFIZS(IX,1,0)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
c
c Find peak densities:
c
        PEAK1 = 0.0
        PEAK2 = 0.0 
        PEAK3 = 0.0
        PEAK4 = 0.0
        PEAK5 = 0.0
        PEAK6 = 0.0
        PEAK7 = 0.0
        PEAK8 = 0.0

        DO IY = -0.5*NYS, 0.5*NYS 
          IF (IY.NE.0) THEN
c            IF (PEAK1.LT.YPRO(IY))    PEAK1 = YPRO(IY)
            IF (PEAK2.LT.NPRO(IY,1))  PEAK2 = NPRO(IY,1)
            IF (PEAK3.LT.INJBINT(IY)) PEAK3 = INJBINT(IY)
            IF (PEAK4.LT.NPRO(IY,-1)) PEAK4 = NPRO(IY,-1)
          ENDIF
        ENDDO
       
        DO IX = 1, MAXNXS
          IF (PEAK6.LT.IPRO(IX,0)) PEAK6 = IPRO(IX,0)
          IF (PEAK7.LT.IPRO(IX,1)) PEAK7 = IPRO(IX,1)
        ENDDO

        DO IP = -MAXNPS, MAXNPS
c          IF (PEAK5.LT.PPRO(IP)) PEAK5 = PPRO(IP) 
        ENDDO
c
c Make sure there are no division by zero errors:
c
        IF (PEAK1.EQ.0.0) PEAK1 = 1.0
        IF (PEAK2.EQ.0.0) PEAK2 = 1.0
        IF (PEAK3.EQ.0.0) PEAK3 = 1.0
        IF (PEAK4.EQ.0.0) PEAK4 = 1.0
        IF (PEAK5.EQ.0.0) PEAK5 = 1.0
        IF (PEAK6.EQ.0.0) PEAK6 = 1.0
        IF (PEAK7.EQ.0.0) PEAK7 = 1.0
c
c X data:
c
        WRITE (63,*) ' '
        WRITE (63,'(A5,A8,2A10,4A8,2X,A10,A10)') 
     +    'IX','x','<ov>',
     +    'nb','Tb','Ti','Peak I','Peak III','Freq. I->II','Drift'

        CIZ = 1
        DO IX = 1, NXS

          IQX = MIN (INT(XS(IX)*XSCALI)+1, NQXSI)                          
          WRITE (63,'(I5,F8.5,2E10.3,2F8.3,2F8.3,E10.3,$)')   
     +      IX,XS(IX),1.0/(CFIZS(IX,1,CIZ)*CRNBS(IX,1)), CRNBS(IX,1), 
     +      CTEMBS(IX,1),CTEMBSI(IX,1),
     +      IPRO(IX,0)/PEAK6,IPRO(IX,1)/PEAK7,
     +      1.0/CFIZS(IX,1,0)

          IF (XS(IX).LT.0.0) THEN
            WRITE(63,'(G10.2)') CVHYS(IX)
          ELSE  
            WRITE(63,'(G10.2)') SVHINS(IX)
          ENDIF
        ENDDO
c
c Output Y integrated profile:
c
        WRITE(63,*) ' '      
        WRITE(63,'(A)') '[BEGIN TOROIDAL PROFILE]'      
        WRITE(63,'(A5,A9,$)') 'IY','Y'
        DO I = 0, MIZS 
          WRITE(63,'(I9,$)') I
        ENDDO
        WRITE(63,*) ' '

        DO IY= -0.5*NYS, 0.5*NYS
          IF (IY.LT.0) THEN

             ! jdemod - code appears to be wanting to write youts but is recalculating it and 
             !          going outside the bounds of YS
             if (iy.eq.-1) then
              WRITE(63,'(I5,F9.4,$)') 
     +          IY,-0.5*(YS(-IY))
            else
              WRITE(63,'(I5,F9.4,$)') 
     +          IY,-0.5*(YS(-IY)+YS(-IY-1))
            endif

            DO IZ = 0, MIZS
              WRITE(63,'(E10.3,$)') YPRO(IY,IZ)+YPRO(IY+NYS+1,IZ)
            ENDDO
            WRITE(63,*) ' ' 
           
c     +      (YPRO(IY,0)+YPRO(IY+NYS+1,0))/PEAK1,
c     +      (YPRO(IY,1)+YPRO(IY+NYS+1,1))/PEAK1,
c     +      (YPRO(IY,3)+YPRO(IY+NYS+1,2))/PEAK1
          ELSEIF (IY.GT.0) THEN

             ! jdemod - code appears to be wanting to write youts but is recalculating it and 
             !          going outside the bounds of YS

             if (iy.eq.1) then 
                WRITE(63,'(I5,F9.4,$)') 
     +           IY,0.5*(YS(IY))
             else
                WRITE(63,'(I5,F9.4,$)') 
     +            IY,0.5*(YS(IY)+YS(IY-1))
             endif

            DO IZ = 0, MIZS
              WRITE(63,'(E10.3,$)') YPRO(IY,IZ)+YPRO(IY-NYS-1,IZ)
            ENDDO
            WRITE(63,*) ' ' 

c     +      (YPRO(IY,0)+YPRO(IY-NYS-1,0))/PEAK1,
c     +      (YPRO(IY,1)+YPRO(IY+NYS+1,1))/PEAK1,
c     +      (YPRO(IY,2)+YPRO(IY+NYS+1,2))/PEAK1
          ENDIF
        ENDDO
        WRITE(63,'(A)') '[END TOROIDAL PROFILE]'      
c
c Output P integrated profile:
c   
        WRITE(63,*) ' '      
        WRITE(63,'(A10,$)') 'P'
        DO I = 0, MIZS 
          WRITE(63,'(I9,$)') I
        ENDDO
        WRITE(63,*) ' '
  
        DO IP= -MAXNPS+1, MAXNPS-1
          IF (IP.LE.0) THEN
            WRITE(63,'(G10.3,$)') 
     +        -0.5*(PS(-IP)+PS(-IP-1))

            DO IZ = 0, MIZS
              WRITE(63,'(E10.3,$)') PPRO(IP,IZ)/PEAK5
            ENDDO
            WRITE(63,*) ' ' 

c     +        PPRO(IP,0)/PEAK5,
c     +        PPRO(IP,1)/PEAK5,
c     +        PPRO(IP,2)/PEAK5
          ELSEIF (IP.GT.0) THEN
            WRITE(63,'(G10.3,$)') 
     +        0.5*(PS(IP)+PS(IP-1))

            DO IZ = 0, MIZS
              WRITE(63,'(E10.3,$)') PPRO(IP,IZ)/PEAK5
            ENDDO
            WRITE(63,*) ' ' 

c     +        PPRO(IP,0)/PEAK5,
c     +        PPRO(IP,1)/PEAK5,
c     +        PPRO(IP,2)/PEAK5
          ENDIF
        ENDDO
c
c Output trace analysis:
c   
        WRITE(63,*) ' '
        WRITE(63,'(A,A14)') 'Trace analysis:     ',' '
        WRITE(63,*) ' '

        MAXIX = 0
        MAXIP = 0
        MAXE  = 1.0E+31

        DO IX = 1, NXS
          FLAG = 0
          DO IP = -MAXNPS+1, MAXNPS-1
            IF (EPRO(IX,IP).GT.0.0) THEN

              WRITE(63,'(2I5,F8.4,$)') IX,IP,XS(IX)

              IF (IP.LE.0) THEN
                WRITE(63,'(F7.4,$)') -0.5*(PS(-IP)+PS(-IP-1))
              ELSE
                WRITE(63,'(F7.4,$)')  0.5*(PS( IP)+PS( IP-1))
              ENDIF
c 
c Only good if there are no parallel gradients:
c
              Cs = SQRT((CTEMBS(IX,1) + CTEMBSI(IX,1)) / CRMB * 
     +                 (1.6E-19/1.67E-27))

              IF (((CRNBS(IX,1) * Cs) / EPRO(IX,IP)).LT.MAXE) THEN
                MAXE  = (CRNBS(IX,1) * Cs) / EPRO(IX,IP)
                MAXIX = IX
                MAXIP = IP
              ENDIF

              WRITE(63,'(4E12.3)') 
     +          EPRO (IX,IP),
     +          CRNBS(IX,1) * Cs,
     +          EPRO(IX,IP) / (CRNBS(IX,1) * Cs),
     +          (CRNBS(IX,1) * Cs) / EPRO(IX,IP)

              FLAG = 1
            ENDIF
          ENDDO
          IF (FLAG.EQ.1) WRITE(63,*) ' '
        ENDDO

        WRITE(63,*) ' '
        WRITE(63,'(A)') 'Maximum injection rate: '
        WRITE(63,'(7X,E10.3,A,2I4,A)') 
     +    MAXE,' particles per second at (IX,IP) = (',
     +    MAXIX,MAXIP,')'
        WRITE(63,*) ' '

        WRITE(7,*) ' '
        WRITE(7,'(A)') 'TRACE IMPURITY ANALYSIS'
        WRITE(7,*) ' '
        WRITE(7,'(A)') 'Maximum injection rate: '
        WRITE(7,'(7X,E10.3,A,2I4,A)') 
     +    MAXE,' particles per second at (IX,IP) = (',
     +    MAXIX,MAXIP,')'
        WRITE(7,*) ' '
c
c Output trace analysis - second time around:
c   
        WRITE(63,*) ' '
        WRITE(63,'(A,A14)') 'Trace analysis II:     ',' '
        WRITE(63,*) ' '

        MAXIX = 0
        MAXIP = 0
        MAXE  = 1.0E+31

        WRITE(63,'(2A5,9A10)')
     +    'IX','IP','EPRO','VOL','NUMi','n','Cs','AREA',
     +    'NUMe','i/e','e/i'

        DO IX = 1, NXS
          FLAG = 0
          DO IP = -MAXNPS+1, MAXNPS-1
            IF (EPRO(IX,IP).GT.0.0) THEN

              WRITE(63,'(2I5,$)') IX,IP
c 
c Only good if there are no parallel gradients:
c
              Cs = SQRT((CTEMBS(IX,1) + CTEMBSI(IX,1)) / CRMB * 
     +                 (1.6E-19/1.67E-27))

              IF (((CRNBS(IX,1)*Cs*EAREA(IX,IP))/
     +             (EPRO(IX,IP)*EVOL(IX,IP))).LT.MAXE) THEN
                MAXE  = (CRNBS(IX,1)*Cs*EAREA(IX,IP))/
     +                  (EPRO(IX,IP)*EVOL(IX,IP))
                MAXIX = IX
                MAXIP = IP
              ENDIF

              WRITE(63,'(9E10.3)') 
     +          EPRO (IX,IP),EVOL(IX,IP),EPRO(IX,IP)*EVOL(IX,IP),
     +          CRNBS(IX,1),Cs,EAREA(IX,IP),
     +          CRNBS(IX,1)*Cs*EAREA(IX,IP),
     +          (EPRO(IX,IP)*EVOL(IX,IP))/(CRNBS(IX,1)*Cs*EAREA(IX,IP)),
     +          (CRNBS(IX,1)*Cs*EAREA(IX,IP))/(EPRO(IX,IP)*EVOL(IX,IP))

              FLAG = 1
            ENDIF
          ENDDO
          IF (FLAG.EQ.1) WRITE(63,*) ' '
        ENDDO

        WRITE(63,*) ' '
        WRITE(63,'(A)') 'Maximum injection rate: '
        WRITE(63,'(7X,E10.3,A,2I4,A)') 
     +    MAXE,' particles per second at (IX,IP) = (',
     +    MAXIX,MAXIP,')'
        WRITE(63,*) ' '

        WRITE(7,'(7X,E10.3,A,2I4,A)') 
     +    MAXE,' particles per second at (IX,IP) = (',
     +    MAXIX,MAXIP,') (number calculation)'
        WRITE(7,*) ' '

c
c Ion injection distribution:
c
c        WRITE (63,*) ' '
c        WRITE (63,'(A10,2A13)') 'Y     ','INJ','FORMULA'
c        DO IY= -0.5*NYS, 0.5*NYS
c          IF (IY.LT.0) THEN
c            WRITE(63,'(G10.3,2G13.4)') 
c     +        -0.5*(YS(-IY)+YS(-IY-1)),
c     +        REAL(INJBINT(IY))/PEAK3,
c     +        EXP(-((0.5*(YS(-IY)+YS(-IY-1))/0.02)**2))
c          ELSEIF (IY.GT.0) THEN
c            WRITE(63,'(G10.3,2G13.4)') 
c     +        0.5*(YS(IY)+YS(IY-1)),
c     +        REAL(INJBINT(IY))/PEAK3,
c     +        EXP(-((0.5*(YS(IY)+YS(IY-1))/0.02)**2))
c           ENDIF
c        ENDDO
c
c      WRITE(63,*) ' '
c      DO IY = -MAXNYS, MAXNYS
c        WRITE(63,*) IY, YWIDS(IY)
c      ENDDO
c
c      WRITE(63,*) ' '
c      DO IX = 1, MAXNXS
c        WRITE(63,*) IX, XWIDS(IX)
c      ENDDO

       RETURN
       END
c
c ======================================================================
c
