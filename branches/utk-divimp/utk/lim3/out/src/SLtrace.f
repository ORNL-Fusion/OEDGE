c
c ======================================================================
c
c subroutine: XLOGSCALE
c
      SUBROUTINE XLogScale
      use mod_params
      use mod_comgra
      use mod_slout
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slout'
c      include 'comgra'

c      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,NPLOTS,ISPOT
c      REAL            CXMIN,CXMAX,CYMIN,CYMAX
c      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT

      REAL decade,decval,t,val,x1

      cxmin = LOG10(cxmin)
      cxmax = LOG10(cxmax)

      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)

      DO decade = REAL(INT(cxmin)-1), REAL(INT(cxmax)+2)
        decval = 10.0**decade

        DO val = decval, decval*10.0, decval
          t = (LOG10(val) - cxmin) / (cxmax - cxmin)

          IF (t.LT.0.0.OR.t.GT.1.0) CYCLE

          x1 = 0.1 + t * (0.9 - 0.1)

          IF (val.EQ.decval) THEN
            CALL POSITN(x1,map1y)
            CALL JOIN  (x1,map1y+0.01)
            CALL POSITN(x1,map2y)
            CALL JOIN  (x1,map2y-0.01)

            CALL PCSCEN(x1,map1y-0.025,'10')
            CALL POSITN(x1+0.006,map1y-0.025+0.009)
            CALL TYPENI(INT(decade))
          ELSE
            CALL POSITN(x1,map1y)
            CALL JOIN  (x1,map1y+0.005)
            CALL POSITN(x1,map2y)
            CALL JOIN  (x1,map2y-0.005)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
c
c ======================================================================
c
c subroutine: YLOGSCALE
c
      SUBROUTINE YLogScale
      use mod_params
      use mod_comgra
      use mod_slout
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slout'
c      include 'comgra'

c      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,NPLOTS,ISPOT
c      REAL            CXMIN,CXMAX,CYMIN,CYMAX
c      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT

      REAL decade,decval,t,val,y1

      cymin = LOG10(cymin+LO)
      cymax = LOG10(cymax+LO)

      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)

      DO decade =  REAL(INT(cymin)-1), REAL(INT(cymax)+2)

        decval = 10.0**decade + LO

        DO val = decval, decval*10.0, decval

          t = (LOG10(val) - cymin) / (cymax - cymin)

          IF (t.LT.0.0.OR.t.GT.1.0) CYCLE

          y1 = 0.11 + t * (0.89 - 0.11)

          IF (val.EQ.decval) THEN
            CALL POSITN(map1x      ,y1)
            CALL JOIN  (map1x+0.010,y1)
            CALL POSITN(map2x      ,y1)
            CALL JOIN  (map2x-0.010,y1)

            CALL PCSCEN(map1x-0.037,y1,'10')
            CALL POSITN(map1x-0.037+0.008,y1+0.009)
            CALL TYPENI(INT(decade))
          ELSE
            CALL POSITN(map1x      ,y1)
            CALL JOIN  (map1x+0.005,y1)

            CALL POSITN(map2x      ,y1)
            CALL JOIN  (map2x-0.005,y1)
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END







