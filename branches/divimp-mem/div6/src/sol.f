      SUBROUTINE SOL(irstart,irend,ikopt)
      use mod_comtor
      IMPLICIT  NONE
C
      integer irstart,irend,ikopt
c
C     INCLUDE   "PARAMS"
      include 'params'
C     INCLUDE   "CGEOM"
      include 'cgeom'
C     INCLUDE   "COMTOR"
c      include 'comtor'
c
C
C  *********************************************************************
C  *                                                                   *
C  *  SOL:   THIS ROUTINE CALCULATES E(S) AND VH(S) IN THE SCRAPE OFF  *
C  *         LAYER.                                                    *
C  *         USES A TEMPORARY ARRAY SS INSTEAD OF KSS TO ENSURE S IS   *
C  *         NEVER QUITE 0 OR SMAX IN CALCULATIONS OF E, VB.           *
C  *                                                                   *
C  *  CHRIS FARRELL  (HUNTERSKIL)  JUNE 1989                           *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IK,IR,ikstart,ikend
      REAL    TEB,TIB,Q,CS,TEMP,M,MSQR,EFACT,VHFACT,SHALF,SS(MAXNKS)
      REAL    SMAX,S,SRAT,K,QT,MT,CFE
      real    s1,s2,vb1,vb2,ds1
C
C
C     Exclude the code here for SOL opt 12. The velocities and
C     Electric fields for this case are set from PLASMA in the
C     SOLEDGE subroutine. However the values still need to be
C     adjusted by the scaling factors thus the if block ends
C     just before the scaling code at the end.
C
      IF (CIOPTF.LT.12) then
C
C
      DO 9900 IR = IRSTART, IREND
c  
       call set_ikvals(ir,ikstart,ikend,ikopt)   
c
       K     = KKS(IR)
       SHALF = 0.5 * KSMAXS(IR)
       SMAX  = KSMAXS(IR)

c
c      Leave this IK reference unchanged because it is a 
c      geometric reference and may be used for IK values that
c      are not strictly within the ik range in use.
c
       DO 10 IK = 1, nks(ir)
         SS(IK) = KSS(IK,IR)
   10  CONTINUE

       IF (SS(1).LE.0.0) SS(1) = SS(1) + 0.01 * (SS(2)-SS(1))

       IF (SS(NKS(IR)).GE.SMAX)
     >   SS(NKS(IR)) = SS(NKS(IR)) - 0.01 * (SS(NKS(IR))-SS(NKS(IR)-1))

C
c
c      Setup code for sol option 6 and 7:
c
       if (cioptf.eq.6.or.cioptf.eq.7) then
          s1 = cvbl1 * smax
          s2 = cvbl2 * smax
          vb1 = cvbm1 * cvhout
          vb2 = cvbm2 * cvhout
          ds1 = s2 - s1
       endif
c
c      Set up the loop index so that the one loop can be set to work 
c      for OUTER, INNER or ENTIRE ring.
c
       DO 9900 IK = ikstart, ikend

        S    = SS(IK)
        IF (S.GT.SHALF) S = SMAX - S
        SRAT = SS(IK) / SHALF - 1.0
        TEB  = KTEBS(IK,IR)
        TIB  = KTIBS(IK,IR)
        CS   = 9.79E+03 * SQRT (0.5*(TEB+TIB)*(1.0+RIZB)/CRMB)
C
C  *********************************************************************
C  *  SOL0:    E(S) AND VH(S) ARE SET TO 0 FOR ALL Y                   *
C  *********************************************************************
C
        IF (CIOPTF.EQ.0) THEN
          KVHS(IK,IR) = 0.0
          KES (IK,IR) = 0.0
C
C  *********************************************************************
C  *  SOL1: E(S) AND VH(S) ARE CALCULATED FROM THE FOLLOWING :-        *
C  *       FOR S < SMAX                                                *
C  *   E(S)  = -(TB/SMAX).M.(1+M.M)/(1-M.M)                            *
C  *   VH(S) = -CS.M                                                   *
C  *   CS    = 9.79E+03.SQRT(TB*(1+RIZB)/CRMB)                         *
C  *   M     = TEMP-SQRT(TEMP.TEMP-1)                                  *
C  *   TEMP  = 1/(1-S/SMAX)                                            *
C  *        AND FOR S > SMAX                                           *
C  *   E(S)  = -E(2SMAX-S)                                             *
C  *   VH(S) = -VH(2SMAX-S)                                            *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.1) THEN
          IF (SRAT.EQ.0.0) THEN
            KVHS(IK,IR) = 0.0
            KES(IK,IR)  = 0.0
          ELSE
            TEMP    = 1.0 / ABS(SRAT)
            M       = SQRT(TEMP * TEMP - 1.0) - TEMP
            MSQR    = M * M
            KVHS(IK,IR) = SIGN((CS * M), SRAT)
            KES(IK,IR)  = SIGN((TEB/SHALF*M*(1.0+MSQR))/(1.0-MSQR),SRAT)
          ENDIF
C
C  *********************************************************************
C  * SOL1A: E(S) AND VH(S) ARE CALCULATED FROM THE FOLLOWING :-        *
C  *       FOR S < SMAX                                                *
C  *   E(S)  = (TB/SMAX).M.(1+M.M)/(1-M.M)                             *
C  *   VH(S) = CS.M                                                    *
C  *   CS    = 9.79E+03.SQRT(TB*(1+RIZB)/CRMB)                         *
C  *   M     = 1/Q + SQRT (1/Q**2-1)                                   *
C  *   Q     = (FS/FL.(S/SMAX-FL)-1)/(1+FS) FOR  S < FL.SMAX           *
C  *           (S/SMAX-1)/(1-FL)/(1+FS)     FOR  FL.SMAX < S < SMAX    *
C  *        AND FOR S > SMAX                                           *
C  *   E(S)  = -E(2SMAX-S)                                             *
C  *   VH(S) = -VH(2SMAX-S)         (SEE NOTE 344)                     *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.-1) THEN
          IF (S.EQ.SHALF) THEN
            KVHS(IK,IR) = 0.0
            KES(IK,IR)  = 0.0
          ELSE
            IF (S.LT.CFL*SMAX) THEN
              Q = (CFS/(2.0*CFL) * (S/SHALF-2.0*CFL) - 1.0) / (1.0+CFS)
            ELSE
              Q = (S/SHALF-1.0) / (1.0+CFS) / (1.0-2.0*CFL)
            ENDIF
            M    = 1.0/Q + SQRT (1.0/(Q*Q) - 1.0)
            MSQR = M * M
            KVHS(IK,IR)  = SIGN (CS*M, SRAT)
            QT = -1.0 / (1.0+CFS)
            MT = (1.0/QT) + SQRT(1.0/(QT*QT) -1.0)
            IF (S.LT.CFL*SMAX) THEN
              CFE = (2.0/(2.0*CFL)) * (0.5+MT/(1.0+MT*MT))
            ELSE
              CFE = (2.0/(2.0*CFL-1.0)) * (MT/(1.0+MT*MT))
            END IF
            IF (MSQR.EQ.1.0) THEN
              KES(IK,IR) = 0.0
            ELSE
              KES(IK,IR) =
     >          SIGN (CFE*TEB/SHALF*M*(1.0+MSQR)/(1.0-MSQR),SRAT)
            ENDIF
          ENDIF
C
C  *********************************************************************
C  *  SOL2: E(S) AND VH(S) ARE CALCULATED FROM THE FOLLOWING FORMULAE  *
C  *       FOR S < SMAX                                                *
C  *  E(S)  = -TB/2L                                                   *
C  *  VH(S) = CS.(1-S/SMAX)                                            *
C  *  CS    = 9.79E+03.SQRT(TB*(1+RIZB)/CRMB)                          *
C  *       AND FOR S > SMAX                                            *
C  *  E(S)  = -E(2SMAX-S)                                              *
C  *  VH(S) = -VH(2SMAX-S)                                             *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.2) THEN
          IF (SRAT.EQ.0.0) THEN
            KVHS(IK,IR) = 0.0
            KES(IK,IR)  = 0.0
          ELSE
            KVHS(IK,IR) = CS * SRAT
            KES(IK,IR)  = SIGN (-TEB/SMAX, SRAT)
          ENDIF
C
C  *********************************************************************
C  *  SOL3:  E(S) AND VH(S) ARE CALCULATED FROM THE FOLLOWING FORMULAE *
C  *        FOR S < SMAX                                               *
C  *   E(S)  = -(TB/(PI.SMAX)).(1-S/SMAX)                              *
C  *   VH(S) = (4/3).CS.(1-S/SMAX)                                     *
C  *   CS    = 9.79E+03.SQRT(TB*(1+RIZB)/CRMB)                         *
C  *        AND FOR S > SMAX                                           *
C  *   E(S)  = -E(2SMAX-S)                                             *
C  *   VH(S) = -VH(2SMAX-S)                                            *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.3) THEN
          IF (SRAT.EQ.0.0) THEN
            KVHS(IK,IR) = 0.0
            KES(IK,IR)  = 0.0
          ELSE
            KVHS(IK,IR) = 1.333333 * CS * SRAT
            KES(IK,IR)  = TEB / (PI*SHALF) * SRAT
          ENDIF
C
C  *********************************************************************
C  *   SOL4: E(S) AND VH(S) ARE CALCULATED FROM THE FOLLOWING FORMULAE *
C  *        FOR S < SMAX                                               *
C  *   E(S)  = -(TB/(2.SMAX)).(1-S/SMAX)                               *
C  *   VH(S) = CS.(1-S/SMAX)                                           *
C  *   CS    = 9.79E+03.SQRT(TB*(1+RIZB)/CRMB)                         *
C  *        AND FOR S > SMAX                                           *
C  *   E(S)  = -E(2SMAX-S)                                             *
C  *   VH(S) = -VH(2SMAX-S)                                            *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.4) THEN
          IF (SRAT.EQ.0.0) THEN
            KVHS(IK,IR) = 0.0
            KES (IK,IR) = 0.0
          ELSE
            KVHS(IK,IR) = CS * SRAT
            KES (IK,IR) = TEB / SMAX * SRAT
          ENDIF
C
C  *********************************************************************
C  *  SOL5:E(S) AND VH(S) ARE CONSTANT, BEARING IN MIND ANY CHANGES OF *
C  *      SIGN IN THE TWO  REGIONS ...                                 *
C  *                                                                   *
C  *                             |   +VH(S)    |   -VH(S)   |          *
C  *                             |    +E(S)    |    -E(S)   |          *
C  *                             |-------------|------------|----> S   *
C  *                             0            SMAX/2        SMAX       *
C  *                                                                   *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.5) THEN
          IF (SRAT.LT.0.0) THEN
            KVHS(IK,IR) = CVHOUT
            KES(IK,IR)  = CEOUT
          ELSE
            KVHS(IK,IR) =-CVHOUT
            KES(IK,IR)  =-CEOUT
          ENDIF
c
        elseif (cioptf.eq.6) then
c
c          SOL option 6:
c
c          The electric field is specifed as a constant
c          as in SOL opt 5. The velocity is a set of fixed values
c          starting at CVHOUT and then set to
c          a value of CVBM1 * CVHOUT at CVBL1 * SMAX and
c          then to CVBM2 * CVHOUT at CVBL2*SMAX until
c          SMAX/2. from L2.
c
           if (s.lt.0.0) then
              kvhs(ik,ir) = cvhout
           elseif (s.le.s1) then
              kvhs(ik,ir) = cvhout
           elseif (s.le.s2) then
              kvhs(ik,ir) = vb1
           else
              kvhs(ik,ir) = vb2
           endif
c
c          Reverse the sign of Vb for the opposite
c          half of the SOL. The value as calculated
c          above is correct for the first half of the
c          SOL.   ( Also E-field)
c
           if (s.gt.s1.and.ofield.eq.2) then
              kes(ik,ir) = 0.0
           else
              kes(ik,ir) = ceout
           endif
c
           if (srat.gt.0.0) then
              kvhs(ik,ir) = -kvhs(ik,ir)
              kes(ik,ir)  = -kes(ik,ir)
           endif
c
c           write (6,*) 'sol6:', ik,ir,kvhs(ik,ir),srat,s,
c     >         s1,s2,vb1,vb2,kes(ik,ir),smax,shalf
c
        elseif (cioptf.eq.7) then
c
c          SOL option 7:
c
c          The electric field is specifed as a constant
c          as in SOL opt 5. The velocity is a set of linear
c          decays starting at CVHOUT and falling to
c          a value of CVBM1 * CVHOUT at CVBL1 * SMAX and
c          then to CVBM2 * CVHOUT at CVBL2*SMAX and then
c          falling to zero at SMAX/2. from L2.
c
           if (s.lt.0.0) then
              kvhs(ik,ir) = cvhout
           elseif (s.le.s1) then
              kvhs(ik,ir) = cvhout + (vb1-cvhout) * s/s1
           elseif (s.le.s2) then
              kvhs(ik,ir) = vb1 + (vb2-vb1) * (s-s1) / ds1
           else
              kvhs(ik,ir) = vb2 * (1.0 - (s-s2)/(shalf-s2))
           endif
c
c          Reverse the sign of Vb for the opposite
c          half of the SOL. The value as calculated
c          above is correct for the first half of the
c          SOL.   ( Also E-field)
c
           if (s.gt.s1.and.ofield.eq.2) then
              kes(ik,ir) = 0.0
           else
              kes(ik,ir) = ceout
           endif
c
           if (srat.gt.0.0) then
              kvhs(ik,ir) = -kvhs(ik,ir)
              kes(ik,ir)  = -kes(ik,ir)
           endif
c
c           write (6,*) 'sol6:', ik,ir,kvhs(ik,ir),srat,s,
c     >         s1,s2,vb1,vb2,kes(ik,ir),smax,shalf
c
c
C
C  *********************************************************************
C  *  SOL9: E(S) AND VH(S) ARE CALCULATED FROM THE FOLLOWING FORMULAE *
C  *       FOR S < SMAX                                                *
C  *  E(S)  = -TB/2L                                                   *
C  *  VH(S) = -CS                                                      *
C  *  CS    = 9.79E+03.SQRT(TB*(1+RIZB)/CRMB)                          *
C  *       AND FOR S > SMAX                                            *
C  *  E(S)  = -E(2SMAX-S)                                              *
C  *  VH(S) = -VH(2SMAX-S)                                             *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.9) THEN
          IF (SRAT.EQ.0.0) THEN
            KVHS(IK,IR) = 0.0
            KES (IK,IR) = 0.0
          ELSE
            KVHS(IK,IR) = SIGN (CS, SRAT)
            KES(IK,IR)  = SIGN (-TEB/SMAX, SRAT)
          ENDIF
C
C  *********************************************************************
C  *  SOL10:   SEE NOTE 353                                            *
C  *********************************************************************
C
        ELSEIF (CIOPTF.EQ.10) THEN
          IF (K.GE.CKIN.AND.K.LE.CKOUT.AND.
     >        S.GE.CFRMIN*SMAX.AND.S.LE.CFRMAX*SMAX) THEN
            KVHS(IK,IR) =
     >        SIGN (CFRM*CS*(CFRMAX-S/SMAX)/(CFRMAX-CFRMIN),-SRAT)
            KES(IK,IR)  = 0.0
          ELSEIF (S.EQ.SHALF) THEN
            KVHS(IK,IR) = 0.0
            KES(IK,IR)  = 0.0
          ELSE
            IF (S.LT.CFL*SMAX) THEN
              Q = (CFS/CFL * (S/SMAX-CFL) - 1.0) / (1.0+CFS)
            ELSE
              Q = (S/SMAX-1.0) / (1.0+CFS) / (1.0-CFL)
            ENDIF
            M       = 1.0/Q + SQRT (1.0/(Q*Q) - 1.0)
            MSQR    = M * M
            KVHS(IK,IR) = SIGN (CS * M, SRAT)
            IF (MSQR.EQ.1.0) THEN
              KES(IK,IR)= 0.0
            ELSE
              KES(IK,IR)= SIGN (TEB/SHALF*M*(1.0+MSQR)/(1.0-MSQR), SRAT)
            ENDIF
          ENDIF
        ENDIF
 9900 CONTINUE
c
c     Set target values to the value in the first cells for electric 
c     field and leave velocity at target equal to sound speed from 
c     initialization routines - this is needed for correct target fluxes.
c
      do ir = irstart,irend
c
c         kvds(idds(ir,2)) = kvhs(1,ir) 
c         kvds(idds(ir,1)) = kvhs(nks(ir),ir) 
c
         keds(idds(ir,2)) = kes(1,ir) 
         keds(idds(ir,1)) = kes(nks(ir),ir) 
c
      end do
c
C     End of SOL 12 block if
C
c
      ENDIF
c
      RETURN
      END

