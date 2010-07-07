C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(110) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J

      ndxa=0
      ndya=0

      read (30,*)
      read (30,*)
      read (30,*)
      read (30,*)

1     continue
      read (30,'(a110)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a110)') zeile

      read (30,'(a110)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
!pb      DO 1015 IY=1,NDYA
!pb        DO 1014 IX=1,NDXA
!pb          X1(IX)=br(ix,iy,1)
!pb          Y1(IX)=bz(ix,iy,1)
!pb          X2(IX)=br(ix,iy,2)
!pb          Y2(IX)=bz(ix,iy,2)
!pb          X3(IX)=br(ix,iy,4)
!pb          Y3(IX)=bz(ix,iy,4)
!pb          X4(IX)=br(ix,iy,3)
!pb          Y4(IX)=bz(ix,iy,3)
!pb1014    CONTINUE
!pb        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
!pb     .                NR1ST,IY)
!pb1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IPART=1
      NPOINT(1,IPART)=1
      IY=1
      DO IX=1,NDXA
        DX=BR(IX+1,IY,1)-BR(IX,IY,2)
        DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
        IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
          NPOINT(2,IPART)=IX+1+IPART-1
          IPART=IPART+1
          NPOINT(1,IPART)=IX+1+IPART-1
        ENDIF
      ENDDO
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1

      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=XPOL(IY,IX)
          Y1(IX)=YPOL(IY,IX)
          X2(IX)=XPOL(IY,IX+1)
          Y2(IX)=YPOL(IY,IX+1)
          X3(IX)=XPOL(IY+1,IX)
          Y3(IX)=YPOL(IY+1,IX)
          X4(IX)=XPOL(IY+1,IX+1)
          Y4(IX)=YPOL(IY+1,IX+1)
1014    CONTINUE
C
        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
     .                NR1ST,IY)
1015  CONTINUE
c
C     do j=1,ndya+1
C       write (iunout,*)
C       write (iunout,*) 'in geomd polygon ',j
C       write (iunout,'(1p,6e12.4)') 
C    .        (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD//
      END
