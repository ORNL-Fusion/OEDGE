
      SUBROUTINE OUTUSR
      USE precision
      USE parmmod
      USE comusr
      USE ccona
      USE cestim
      USE ctrcei
      USE comxs
      USE cgeom
      USE cinit
      USE ctrig
      USE csdvi
      USE clogau
      USE comsou
      IMPLICIT none

c      call prousr (clst,1+6*npls,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,
c     .             0._dp,ntrii)
c      ICELLRD=CLST
c      call prousr (clst,2+6*npls,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,
c     .             0._dp,ntrii)
c      IRINGRD=CLST

c      DO IR=1,NTRII
c         IF (IYTRI(IR) == -1) THEN
c            ITRI = IXTRI(IR)
c            NOTRI(ITRI) = NOTRI(ITRI) + 1
c            ICELLST(ITRI) = ICELLRD(IR)
c            IRINGST(ITRI) = IRINGRD(IR)
c         ENDIF
c      ENDDO

c      MTRI=COUNT(NOTRI>0)
c      DO ITRI=1,MTRI
c      ENDDO

c  write output file for OSM code, from last iteration

c      WRITE (96,'(I6)') MTRI
c      WRITE (96,'(A6,9A12)') 'NO.','VOL ',
c     .     'TABDS1','TABDS2', 'CELL INDEX','RING INDEX' 
c      DO IR=1,MTRI
c         WRITE (96,'(I6,3ES12.4,2I12)') IR, VOLOUT(IR),
c     .         TDNE(IR), TDNE2(IR), ICELLST(IR), IRINGST(IR)
c      END DO

c  output for OSM from last iteration: done

      WRITE(0,*) 'HERE IN OUTUSR'

      RETURN
 99   STOP
      END
