C
C
      SUBROUTINE XSECTPH_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMUSR
      USE PHOTON
      IMPLICIT NONE
csw
csw  PHOTON COLLISION, OT - type
csw
      integer :: iphot,nrc,kk,nnrot

      nnrot=0
      do iphot=1,nphoti
         if(nrcph(iphot) > 0) then
            do nrc=1,nrcph(iphot)
               kk=ireacph(iphot,nrc)
               if(iswr(kk) == 7) then
                  nNROT=nNROT+1
               endif
            enddo
         endif
      enddo
      call PH_ALLOC_XSECTPH(nnrot)

      RETURN

      END SUBROUTINE XSECTPH_PARAM
