      SUBROUTINE TEST_TRACE_N0()
      USE GEOMETRY_PL
      USE SURFACE_PL
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      LOGICAL :: SURF_P
      NZ0_N0 = 1
      JR0_N0 = 10
      JP0_N0 = 60
      JT0_N0 = 17
       
      RJ0 = 0.
      PJ0 = 0.
      TJ0 = JT0_N0 + 0.8

      IG = JR0_N0+(JP0_N0+JT0_N0*ZON_POLO(NZ0_N0))*ZON_RADI(NZ0_N0)
     .+    MESH_P_OS(NZ0_N0)
      IC0_N0 = IDCELL(IG)
  
      CALL RZ_REAL_COORDINATES
     .          (NZ0_N0,JR0_N0,JP0_N0,RJ0,PJ0,TJ0,R0,Z0)
      K1 = JT0_N0 + PHI_PL_OS(NZ0_N0)
      K2 = K1  + 1
      X0 = R0*(COSPHI(K1)+(COSPHI(K2)-COSPHI(K1))*(TJ0-JT0_N0))
      Y0 = R0*(SINPHI(K1)+(SINPHI(K2)-SINPHI(K1))*(TJ0-JT0_N0))

      VELX =-1.
      VELY = 0.
      VELZ = 0.
      NEW  = 0
      SURF_P = .FALSE.
      NPANU  = 1
      
      ITR = 0
      OPEN(99,FILE='out')
      write(99,'(3f12.6)')x0,y0,z0
      DO 
            CALL TIMUSR(IRGEN,X0,Y0,Z0,VELX,VELY,VELZ,NEW,
     .                  ICNXT,TIMET,ICOS,IER,ITR,SURF_P)
            IF(TIMET <= 0.) THEN 
              WRITE(6,*)'TIMET=',TIMET
              EXIT
            ELSEIF(ICNXT == -1) THEN
              CALL NORUSR(ISTS,X0,Y0,Z0,XNORM,YNORM,ZNORM
     .,                   SCOS,VELX,VELY,VELZ,IRGEN)
              TIMET = 0.
              SURF_P=.TRUE.
              NEW   = 0
c             write(99,'(3f12.6)')xfin,yfin,zfin
c             ITR = ITR + 1
c             IF(ITR == 100) EXIT
            ELSEIF(ICNXT<0) THEN
              X0 = X0 + VELX*TIMET
              Y0 = Y0 + VELY*TIMET
              Z0 = Z0 + VELZ*TIMET
              CALL NORUSR(ISTS,X0,Y0,Z0,XNORM,YNORM,ZNORM
     .,                   SCOS,VELX,VELY,VELZ,IRGEN)
              TIMET = 0.
              SURF_P=.TRUE.
              NEW   = 0 
              A     = VELX*XNORM + VELY*YNORM + VELZ*ZNORM
              VELX  = VELX*(1.-XNORM**2)-(2.*A-VELX*XNORM)*XNORM
              VELY  = VELY*(1.-YNORM**2)-(2.*A-VELY*YNORM)*YNORM
              VELZ  = VELZ*(1.-ZNORM**2)-(2.*A-VELZ*ZNORM)*ZNORM

c             write(99,'(3f12.6)')xfin,yfin,zfin
              ITR = ITR + 1
              IF(ITR == 100000) EXIT
            ENDIF
      ENDDO 
      close(99)
      END SUBROUTINE TEST_TRACE_N0
      SUBROUTINE TEST_LEAUSR()
      USE GEOMETRY_PL
      USE SURFACE_PL
      USE PHYSICAL_CELL
      USE NEUTRAL_TRANSPORT

      DO I=1,100000
      NZ0 = 2.*RANF()
      JR0 = RANF()*ZON_RADI(NZ0)
      JP0 = RANF()*ZON_POLO(NZ0)
      JT0 = RANF()*ZON_TORO(NZ0)

      RJ0 = 0.
      PJ0 = 0.
      TJ0 = JT0 + 0.00005

      IG = JR0+(JP0+JT0*ZON_POLO(NZ0))*ZON_RADI(NZ0)
     .+    MESH_P_OS(NZ0)
      IC0= IDCELL(IG)
 
      CALL RZ_REAL_COORDINATES
     .          (NZ0,JR0,JP0,RJ0,PJ0,TJ0,R0,Z0)
      K1 = JT0 + PHI_PL_OS(NZ0)
      K2 = K1  + 1
      X0 = R0*(COSPHI(K1)+(COSPHI(K2)-COSPHI(K1))*(TJ0-JT0))
      Y0 = R0*(SINPHI(K1)+(SINPHI(K2)-SINPHI(K1))*(TJ0-JT0))

      II = LEAUSR(X0,Y0,Z0,I)
      IF( NZ0_N0/=NZ0 .OR. JR0_N0/=JR0 .OR. JP0_N0/=JP0 .OR.
     .    JT0_N0/=JT0 .OR. IC0_N0/=IC0) THEN
      WRITE(6,*)'I=',I
      WRITE(6,'(5I6)')NZ0,JR0,JP0,JT0,IC0
      WRITE(6,'(5I6)')NZ0_N0,JR0_N0,JP0_N0,JT0_N0,IC0_N0          
      ENDIF
     
      ENDDO 
      END
