 
 
      subroutine EIRENE_B_PROJI(B,EB,V,V_PARALLEL,V_PERP,PHI)
C
C  INVERS OF B_PROJ: GIVEN B(3), V_PARALLEL, V_PERP AND PHI
C  CALCULATE V(3)
C  ENTRY B_PROJI : FIRST NORMALIZE B TO EB.
C  ENTRY B_PROJIN: EB IS KNOWN TO BE THE UNIT VECTOR OF B
C  INPUT:  B(3), V_PARALLEL, V_PERP, PHI
C          SIGN OF V_PARALLEL IS WITH RESPECT TO B(3)
C  OUTPUT: V(3)
C
      USE EIRMOD_PRECISION
      implicit none
      REAL(DP), INTENT(IN) :: B(3), V_PARALLEL, V_PERP, PHI
      REAL(DP), INTENT(INOUT) :: EB(3),V(3)
      REAL(DP) :: V_P(3), E1(3), E2(3), E3(3)
      REAL(DP) :: VQ, BNI, VNI, B12I, B12
 
      BNI=1./SQRT(SUM(B*B)+1.D-30)
      EB=B*BNI
C
      ENTRY EIRENE_B_PROJIN(B,EB,V,V_PARALLEL,V_PERP,PHI)
C
      E1=EB
C
      B12=B(1)**2+B(2)**2
      IF (B12.GT.1.D-30) THEN
        B12I=1./SQRT(B12)
        E2(1)=-B(2)*B12I
        E2(2)= B(1)*B12I
        E2(3)= 0.
      ELSE
        E2(1)=1.
        E2(2)=0.
        E2(3)=0.
      ENDIF
C  E3 = E1 X E2
      E3(1)= E1(2)*E2(3)-E1(3)*E2(2)
      E3(2)=-E1(1)*E2(3)+E1(3)*E2(1)
      E3(3)= E1(1)*E2(2)-E1(2)*E2(1)
C
      V_P(1)=V_PARALLEL
      V_P(2)=COS(PHI)*V_PERP
      V_P(3)=SIN(PHI)*V_PERP
C  V = [E1,E2,E3] * V_PER
      V(1)=V_P(1)*E1(1)+V_P(2)*E2(1)+V_P(3)*E3(1)
      V(2)=V_P(1)*E1(2)+V_P(2)*E2(2)+V_P(3)*E3(2)
      V(3)=V_P(1)*E1(3)+V_P(2)*E2(3)+V_P(3)*E3(3)
      VNI=1./SQRT(SUM(V**2)+1.D-30)
      V = V * VNI
      return
      end
