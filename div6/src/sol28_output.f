c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE ListTargetData(fp,title)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none      

      INTEGER  , INTENT(IN) :: fp
      CHARACTER, INTENT(IN) :: title*(*)

      INTEGER     itarget,itube,ion
      CHARACTER*2 target_tag(2) 

      ion = 1

      target_tag(LO) = 'LO'
      target_tag(HI) = 'HI'

      WRITE(fp,*)
      WRITE(fp,'(A)') 'TARGET DATA:'
      WRITE(fp,'(A)') '  '//TRIM(title)
      DO itarget = LO, HI
        WRITE(fp,*)
        WRITE(fp,'(A6,A16,3A10,A6,4A10,A8,2X,A)') 
     .    'TUBE','jsat','ne','ni','vi','M','pe','pi','Te','Ti','Gamma',
     .     target_tag(itarget)
        DO itube = 1, ntube
          IF (tube(itube)%type.EQ.GRD_CORE) CYCLE
          WRITE(fp,'(I6,1P,E16.6,3E10.2,0P,F6.2,
     .               1P,2E10.2,0P,2F10.6,F8.2)')
     .      itube,
     .      tube(itube)%jsat  (itarget,ion),
     .      tube(itube)%ne    (itarget),
     .      tube(itube)%ni    (itarget,ion),
     .      tube(itube)%vi    (itarget,ion),
     .      tube(itube)%machno(itarget),
     .      tube(itube)%pe    (itarget),
     .      tube(itube)%pi    (itarget,ion),
     .      tube(itube)%te    (itarget),
     .      tube(itube)%ti    (itarget,ion),
     .      tube(itube)%gamma (itarget,ion)
        ENDDO
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c ======================================================================
c
      SUBROUTINE GenerateOutputFiles
      IMPLICIT none

      CALL DumpData_OSM('output.end','Simulation complete')

c...  Save solution:
      CALL SaveGrid('osm.raw')

      CALL User_GenerateOutputFiles

      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE DumpData_OSM(fname,title)
      USE mod_sol28_params
      USE mod_sol28_global
      IMPLICIT none

      CHARACTER*(*) fname,title 

      INTEGER fp,it,ic,ion


      CALL User_DumpData(fp,fname,title)


      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=98)            

      
      WRITE(fp,'(A)') 'TITLE: '//TRIM(title)

      WRITE(fp,*)
      WRITE(fp,10) 'NTUBE     = ',nobj     ,'NGROUP   = ',ngrp
      WRITE(fp,10) 'NSRF      = ',nsrf     ,'NVTX     = ',nvtx
      WRITE(fp,*)
      WRITE(fp,10) 'NTUBE     = ',ntube    ,'NCELL    = ',ncell
      WRITE(fp,10) 'NION      = ',nion     ,'NFIELD   = ',nfield
      WRITE(fp,10) 'NFLUID    = ',nfluid   ,'NKINETIC = ',nkinetic
      WRITE(fp,10) 'NPIN      = ',npin     ,'NPHOTON  = ',nphoton
      WRITE(fp,10) 'NIMPURITY = ',nimpurity,'NDRIFT   = ',ndrift

      WRITE(fp,*)
      WRITE(fp,'(A)') 'GRID DATA:'
      WRITE(fp,*)

      WRITE(fp,*)
      WRITE(fp,'(A)') 'TUBE DATA:'
      WRITE(fp,*)
      WRITE(fp,'(A6,2A9,7A10)') 'Index','Cell_LO','Cell_HI','psi_n',
     .                          'costet_LO','costet_HI',
     .                          'rp_LO'    ,'rp_HI'    ,
     .                          'dds_LO'   ,'dds_HI'   

      DO it = 1, ntube
        WRITE(fp,'(I6,2I9,9F10.4)')
     .    it,tube(it)%cell_index(LO),tube(it)%cell_index(HI),
     .    tube(it)%psin,
     .    tube(it)%costet(1:2),
     .    tube(it)%rp    (1:2),
     .    tube(it)%dds   (1:2)
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A)') 'TARGET DATA:'
      WRITE(fp,*)

      WRITE(fp,*)
      WRITE(fp,'(A)') 'GEOMETRY DATA:'
      DO it = 1, ntube
        DO ion = 1, nion
          WRITE(fp,*)
          WRITE(fp,10) '   TUBE =',it
          WRITE(fp,10) '   ION  =',ion
          WRITE(fp,'(A8,7A11)') 
     .      'CELL ','s','sbnd1','sbnd2','p','q','R','Z'
          WRITE(fp,'(A8,7A11)') 
     .      '     ','(m)',' ',' ','(m)',' ','(m)','(m)'
          WRITE(fp,'(8X,1P,5E11.2,0P,2F11.6)')
     .       0.0,
     .      -1.0,
     .      -1.0,
     .       0.0,
     .      -1.0,
     .      -1.0,
     .      -1.0
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
            WRITE(fp,'(I8,1P,5E11.2,0P,2F11.6)') ic,
     .        cell(ic)%s,
     .        cell(ic)%sbnd(1:2),
     .        cell(ic)%p,
     .        -1.0,
     .        cell(ic)%cencar(1:2)
          ENDDO
          WRITE(fp,'(8X,1P,5E11.2,0P,2F11.6)')
     .      tube(it)%smax,
     .      -1.0,
     .      -1.0,
     .      tube(it)%pmax,
     .      -1.0,
     .      -1.0,
     .      -1.0
        ENDDO
      ENDDO


      WRITE(fp,*)
      WRITE(fp,'(A)') 'PLASMA DATA:'
      DO it = 1, ntube
        DO ion = 1, nion
          WRITE(fp,*)
          WRITE(fp,10) '   TUBE =',it
          WRITE(fp,10) '   ION  =',ion
          WRITE(fp,'(A8,3A11,2A10)') 
     .      'CELL ','ne','ni','vi','Te','Ti'
          WRITE(fp,'(A8,3A11,2A10)') 
     .      '     ','(m-3)','(m-3)','(m s-1)','(eV)','(eV)'
          WRITE(fp,'(8X,1P,3E11.2,0P,2F10.2)') 
     .      tube(it)%ne(LO),
     .      tube(it)%ni(LO,ion),
     .      tube(it)%vi(LO,ion),
     .      tube(it)%te(LO),
     .      tube(it)%ti(LO,ion)
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
            WRITE(fp,'(I8,1P,3E11.2,0P,2F10.2)') ic,
     .        fluid(ic,ion)%ne,
     .        fluid(ic,ion)%ni,
     .        fluid(ic,ion)%vi,
     .        fluid(ic,ion)%te,
     .        fluid(ic,ion)%ti
          ENDDO
          WRITE(fp,'(8X,1P,3E11.2,0P,2F10.2)') 
     .      tube(it)%ne(HI),
     .      tube(it)%ni(HI,ion),
     .      tube(it)%vi(HI,ion),
     .      tube(it)%te(HI),
     .      tube(it)%ti(HI,ion)
        ENDDO
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A)') 'NEUTRAL DATA:'
      DO it = 1, ntube
        DO ion = 1, nion
          WRITE(fp,*)
          WRITE(fp,10) '   TUBE =',it
          WRITE(fp,10) '   ION  =',ion
          WRITE(fp,'(A8,2(A11,6X),4A11)') 
     .      'CELL ','n_atm','n_mol','ion','rec','Dalpha','Dgamma'
          WRITE(fp,'(A8,2(A11,6X),4A11)') 
     .      '     ','(m-3)','(m-3)','(m-3 s-1)','(m-3 s-1)','(m-3 s-1)',
     .      '(m-3 s-1)'
          DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
            WRITE(fp,'(I8,2(1P,E11.2,0P,F6.2),1P,4E11.2,F6.2)') 
     .        ic,
     .        pin(ic,ion)%n_atm,
     .        pin(ic,ion)%n_atm / (fluid(ic,ion)%ne + EPS10),
     .        pin(ic,ion)%n_mol,
     .        pin(ic,ion)%n_mol / (fluid(ic,ion)%ne + EPS10),
     .        pin(ic,ion)%ion,
     .        pin(ic,ion)%rec,
     .        pin(ic,ion)%dalpha(1),
     .        pin(ic,ion)%dgamma(1),
     .        pin(ic,ion)%dgamma(1) / (pin(ic,ion)%dalpha(1) + EPS10)
          ENDDO
        ENDDO
      ENDDO

      WRITE(fp,*)
      WRITE(fp,'(A)') 'GEOMETRY DATA (FLUID GRID):'
      WRITE(fp,*)

      WRITE(fp,*)
      WRITE(fp,'(A)') 'FIELD DATA:'
      DO it = 1, ntube
        WRITE(fp,*)
        WRITE(fp,10) '   TUBE =',it
        WRITE(fp,'(A8,5A11)') 
     .    'CELL ','Brat','Btot','Br','Bphi','Bz'
        WRITE(fp,'(A8,5A11)') 
     .    '     ','(Tesla)','(Tesla)','(Tesla)','(Tesla)','(Tesla)'
        DO ic = tube(it)%cell_index(LO), tube(it)%cell_index(HI)        
          WRITE(fp,'(I8,5F11.5)') ic,
     .      field(ic)%bratio,
     .      field(ic)%b,
     .      field(ic)%br,
     .      field(ic)%bphi,
     .      field(ic)%bz
        ENDDO
      ENDDO

      CLOSE(fp)

 10   FORMAT(10(A,I8,4X))


      RETURN
 98   CALL ER('DumpSolution','Trouble accessing file',*99)
 99   WRITE(0,*) '  FILE NAME= "'//TRIM(fname)//'"'
      STOP
      END




