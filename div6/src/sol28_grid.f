c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE ProcessGrid
      USE mod_geometry
      USE mod_legacy
      IMPLICIT none

      INTEGER status 

      CALL LoadObjects('osm_geometry.raw',status)  ! Change to LoadGeometryObjects...
      IF (status.EQ.-1) 
     .  CALL ER('SetupGrid','Geometry data not found',*99)

      CALL LoadGrid('osm.raw')

      CALL LoadLegacyData('osm_legacy.raw')

      RETURN
 99   STOP
      END
