c     -*-Fortran-*-
c
      DOUBLE PRECISION FUNCTION CONTIN(WAVE, NE, TEV, NZS, 
     >                                 MAXIZS, ZIMP, NH, NZ)
      IMPLICIT NONE
C
      INTEGER   NZS, MAXIZS, ZIMP(NZS)
      REAL*8    WAVE, NE, TEV, NH, NZ(0:MAXIZS,NZS)
C
C  PHYSICAL CONSTANTS
C
      REAL*8 A0, ALPHA, C, IH, H, K, E, PI
      PARAMETER (A0 = 5.2918D-9, ALPHA = 7.2974D-3)
      PARAMETER (C = 2.9979D10, IH = 13.606)
      PARAMETER (H = 6.6261D-27, K = 1.3807D-16)
      PARAMETER (E = 1.6022D-12, PI = 3.14159)
C
C  VARIABLES FOR CALCULATION EMISSIVITIES
C
      INTEGER    N, NMIN, NMAX
      INTEGER    IE, ISUM, NFULL, IZ, EREM, J
      REAL*8     Q1, LAMBDA, TE, NU, U, U2, GAM2
      REAL*8     Z, FACT
      REAL*8     EN, EPSBF, GBF, BOTH, GAV
c     
      real*8     tmpcontin, tmpboth, tmpzcontin,tmpzboth
c
      tmpcontin = 0.0
      tmpboth = 0.0
      tmpzcontin = 0.0
      tmpzboth = 0.0
C
C  CONSTANTS FOR THE CALCULATIONS
C
      Q1 = 32.0*ALPHA**4.0*C/(3.0*DSQRT(3.0D0)*PI*A0)
      LAMBDA = WAVE*1.0D-8
      TE = TEV*11604.0
      NU = C/LAMBDA
      U = H*NU/(K*TE)
      CONTIN = 0.0
C
C  START WITH HYDROGEN
C
      Z = 1.0
      NMIN = IDINT(DSQRT(Z*Z*E*IH/(H*NU))) + 1
C  USE DISCRETE CALCULATION FOR FIRST TEN LEVELS
      NMAX = NMIN + 9
C
C  FREE-BOUND CONTINUA - LOW LEVELS
C
      DO N = NMIN, NMAX
        EN = DFLOAT(N)
        U2 = EN*EN*H*NU/(Z*Z*E*IH) - 1.0
C
        EPSBF = NE*NH*2.0*Q1*(PI*A0*A0*IH/TEV)**1.5*GBF(EN,U2)*Z**4.0
     >          /(EN**3.0)*DEXP(Z*Z*IH/(EN*EN*TEV))*DEXP(-U)/WAVE
C
        CONTIN = CONTIN + EPSBF
      ENDDO
c
      tmpcontin = contin

C
C  HIGH LEVEL FREE-BOUND CONTINUA + FREE-FREE
C
      GAM2 = Z*Z*IH/TEV
      BOTH  = NE*NH*Z*Z*Q1*(PI*A0*A0*IH/TEV)**1.5*(TEV/IH)
     >        *GAV(U,GAM2,NMAX+1)*DEXP(-U)/WAVE
C
c
      tmpboth = both
c

      CONTIN = CONTIN + BOTH
C
C  ADD IMPURITY CONTRIBUTIONS
C
      DO J = 1, NZS
        DO IZ = 1, ZIMP(J)
          Z = DFLOAT(IZ)
          NMIN = IDINT(DSQRT(Z*Z*E*IH/(H*NU))) + 1
C  SKIP FULL LEVELS
          IE = ZIMP(J) - IZ
          ISUM = 0
          NFULL = 0
          EREM = IE
          DO N = 1, 5
            ISUM = ISUM + 2*N*N
            IF (ISUM.LT.IE) THEN
              NFULL = N
              EREM = IE - ISUM
            ENDIF
          ENDDO
          FACT = DFLOAT(2*(NFULL+1)*(NFULL+1)-EREM)
     >           /DFLOAT(2*(NFULL+1)*(NFULL+1))
          IF (NMIN.LT.NFULL+1) NMIN = NFULL + 1
C  USE DISCRETE CALCULATION FOR FIRST TEN LEVELS
          NMAX = NMIN + 9
C
C  FREE-BOUND CONTINUA - LOW LEVELS
C
          DO N = NMIN, NMAX
            EN = DFLOAT(N)
            U2 = EN*EN*H*NU/(Z*Z*E*IH) - 1.0
C
            EPSBF = NE*NZ(IZ,J)*2.0*Q1*(PI*A0*A0*IH/TEV)**1.5*GBF(EN,U2)
     >              *Z**4.0/(EN**3.0)*DEXP(Z*Z*IH/(EN*EN*TEV))
     >              *DEXP(-U)/WAVE
C
C  USE PARTIAL OCCUPANCY FACTOR, IF NECESSARY
C
            IF (N.EQ.NFULL+1) EPSBF = EPSBF*FACT
C
            CONTIN = CONTIN + EPSBF
c
            tmpzcontin = tmpzcontin + epsbf
c
          ENDDO
C
C  HIGH LEVEL FREE-BOUND CONTINUA + FREE-FREE
C
          GAM2 = Z*Z*IH/TEV
          BOTH  = NE*NZ(IZ,J)*Z*Z*Q1*(PI*A0*A0*IH/TEV)**1.5*(TEV/IH)
     >            *GAV(U,GAM2,NMAX+1)*DEXP(-U)/WAVE
C
          tmpzboth = tmpzboth + both 
c
          CONTIN = CONTIN + BOTH
        ENDDO
      ENDDO
c
c
      write (6,'(a,2i4,20(1x,g10.3))') 'CONTIN:',nzs,zimp(nzs),
     >             tev,ne,nh,
     >             (nz(iz,1),iz=1,zimp(1)),
     >             tmpcontin,tmpboth,tmpzcontin,tmpzboth,contin
c      
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION CONTFF(WAVE, NE, TEV, NZS,
     >                                 MAXIZS, ZIMP, NH, NZ)
      IMPLICIT NONE
C
      INTEGER   NZS, MAXIZS, ZIMP(NZS)
      REAL*8    WAVE, NE, TEV, NH, NZ(0:MAXIZS,NZS)
C
C  PHYSICAL CONSTANTS
C
      REAL*8 A0, ALPHA, C, IH, H, K, E, PI
      PARAMETER (A0 = 5.2918D-9, ALPHA = 7.2974D-3)
      PARAMETER (C = 2.9979D10, IH = 13.606)
      PARAMETER (H = 6.6261D-27, K = 1.3807D-16)
      PARAMETER (E = 1.6022D-12, PI = 3.14159)
C
C  VARIABLES FOR CALCULATION EMISSIVITIES
C
      INTEGER    IZ, J
      REAL*8     Q1, LAMBDA, TE, NU, U, GAM2
      REAL*8     Z
      REAL*8     BOTH, GIIIAV
C
C  CONSTANTS FOR THE CALCULATIONS
C
      Q1 = 32.0*ALPHA**4.0*C/(3.0*DSQRT(3.0D0)*PI*A0)
      LAMBDA = WAVE*1.0D-8
      TE = TEV*11604.0
      NU = C/LAMBDA
      U = H*NU/(K*TE)
      CONTFF = 0.0
C
C  START WITH HYDROGEN
C
      Z = 1.0
      GAM2 = Z*Z*IH/TEV
      BOTH  = NE*NH*Z*Z*Q1*(PI*A0*A0*IH/TEV)**1.5*(TEV/IH)
     >        *GIIIAV(U,GAM2)*DEXP(-U)/WAVE
C
      CONTFF = CONTFF + BOTH
C
C  ADD IMPURITY CONTRIBUTIONS
C
      DO J = 1, NZS
        DO IZ = 1, ZIMP(J)
          Z = DFLOAT(IZ)
          GAM2 = Z*Z*IH/TEV
          BOTH  = NE*NZ(IZ,J)*Z*Z*Q1*(PI*A0*A0*IH/TEV)**1.5*(TEV/IH)
     >            *GIIIAV(U,GAM2)*DEXP(-U)/WAVE
C
          CONTFF = CONTFF + BOTH
        ENDDO
      ENDDO
C
      RETURN
      END
C
C
C
       DOUBLE PRECISION FUNCTION GBF(EN,U)
       IMPLICIT REAL*8(A-H,O-Z)                                         
C  WHERE U = N*N/K*K = N*N*E/Z*Z                                        
       real*8 en,u,x,t1,t2
       X=(EN*(U+1.0))**0.6666667                                        
       T1=0.1728*(U-1.0)/X                                              
       T2=1.3888889*T1*T1+1.333333*(0.0496*(U*U+1.333333*U              
     1 +1.0))/(X*X)                                                     
       GBF=1.0/(1.0-1.333333*T1+T2)**0.75                               
       RETURN                                                           
      END                                                               
C									
C									
C									
      DOUBLE PRECISION FUNCTION GIIAV(U,U2,N0)				
      IMPLICIT NONE							
C									
C  AVERAGED FREE-BOUND GAUNT FACTOR FOR SUMMED CONTRIBUTIONS		
C  OF THE HIGH LYING STATES - SUMMERS AND HOOPER EQU. 19		
C  ******  L.D. HORTON, JET    24 JUN 1997  ***************		
C  INPUT								
C      U=HV/IH  WHERE HV IS PHOTON ENERGY				
C               AND IH IS THE RYDBERG ENERGY				
C      U2=K*TE/IH  WHERE TE IS THE ELECTRON TEMPERATURE			
C                  AND K IS THE BOLTZMANN CONSTANT			
C      N0       IS THE FIRST LEVEL TO BE COUNTED IN THIS		
C               QUASI-CONTINUUM						
C  OUTPUT								
C      GIIAV=AVERAGED FREE-BOUND GAUNT FACTOR.				
C									
      INTEGER N0							
      REAL*8  U, U2							
C									
      INTEGER NSIMP							
      PARAMETER (NSIMP = 400)						
C									
      INTEGER I								
      REAL*8  X, XMIN, XMAX, EN, H, U3, SUM, EXPX			
      REAL*8  GBF, GLIM							
C									
      XMIN = -1.0D0/(N0*N0*U2)						
      XMAX = 0.0							
      H = (XMAX-XMIN)/(2.0*NSIMP)					
      EN = N0								
      U3 = EN*EN*U - 1.0						
      SUM = GBF(EN,U3)*DEXP(-XMIN)					
      DO I = 1,2*NSIMP-1,2						
        X = XMIN + H*I							
        EN = 1.0/DSQRT(-X*U2)						
        U3 = EN*EN*U - 1.0						
C       WRITE(6,*) 'I,X,U3,EN:', I, X, U3, EN				
        SUM = SUM + 4.0*GBF(EN,U3)*DEXP(-X)				
      ENDDO								
      DO I = 2, 2*NSIMP-2,2						
        X = XMIN + H*I							
        EN = 1.0/DSQRT(-X*U2)						
        U3 = EN*EN*U - 1.0						
C       WRITE(6,*) 'I,X,U3,EN:', I, X, U3, EN				
        SUM = SUM + 2.0*GBF(EN,U3)*DEXP(-X)				
      ENDDO								
      GLIM = 1.0/(1.0-0.2304*U**0.3333333+0.1076*U**0.6666667)**0.75	
      SUM = SUM + GLIM							
      GIIAV = H/3.0*SUM							
      RETURN								
      END								
C									
C									
C									
       DOUBLE PRECISION FUNCTION GAV(U,GAM2,N0)
       IMPLICIT REAL*8(A-H,O-Z)						
C  CALCULATES TOTAL GAUNT FACTOR FOR FREE-FREE AND QUASI-CONTINUOUS	
C  FREE-BOUND TRANSITIONS						
C  ******  L.D. HORTON, JET    24 JUN 1997  ***************		
C  INPUT								
C      U=HV/KT  WHERE HV IS PHOTON ENERGY				
C               AND KT IS ELECTRON TEMPERATURE (ENERGY UNITS)		
C      GAM2=Z*Z*IH/KT  WHERE Z IS TARGET ION CHARGE			
C                      AND IH IS THE RYDBERG ENERGY			
C      N0       IS THE FIRST BOUND LEVEL TO BE INCLUDED IN THE		
C               INTEGRAL						
C  OUTPUT								
C      GAV=MAXWELLIAN AVERAGED GAUNT FACTOR.				
       integer i 
       real*8 u,gam2,n0,u0,rgam2,sum,x,e1,e2,en,u1,gbf,giii

       REAL*8 XA(8),WA(8)						
       DATA XA/0.17027963D0,0.90370178D0,2.25108663D0,4.26670017D0,	
     &7.04590540D0,10.75851601D0,15.74067864D0,22.86313174D0/		
       DATA WA/3.69188589D-1,4.18786781D-1,1.75794987D-1,		
     &3.33434923D-2,2.79453624D-3,9.07650877D-5,8.48574672D-7,		
     &1.04800117D-9/							
       U0 = GAM2/(N0*N0)						
C      WRITE(6,*) ' U0 ', U0						
       RGAM2=1.0D0/GAM2							
       SUM=0.0D0							
       DO 20 I=1,8							
       X=XA(I)								
C      WRITE(6,*) ' X, X-U0, X-U0+U ', X, X-U0, X-U0+U			
       E1=RGAM2*(X-U0)							
       E2=RGAM2*(X-U0+U)						
C      WRITE(6,*) E1,E2							
       IF (E1.LT.0.0) THEN						
         EN = 1.0/DSQRT(-E1)						
         U1 = -E2/E1							
C        WRITE(6,*) EN,U1,GBF(EN,U1)					
         SUM = SUM + WA(I)*GBF(EN,U1)					
       ELSE								
C        WRITE(6,*) GIII(1,1,E1,E2)					
         SUM=SUM+WA(I)*GIII(1,1,E1,E2)					
       ENDIF								
   20  CONTINUE								
       GAV=DEXP(U0)*SUM							
       RETURN								
      END								
C									
C									
C									
       REAL*8 FUNCTION GIIIAV(U,GAM2)			
       IMPLICIT REAL*8(A-H,O-Z)						
C  CALCULATES MAXWELLIAN AVERAGED FREE-FREE GAUNT FACTORS		
C  ******  H.P. SUMMERS, JET    22 NOV 1984  ***************		
C  INPUT								
C      U=HV/KT  WHERE HV IS PHOTON ENERGY				
C               AND KT IS ELECTRON TEMPERATURE (ENERGY UNITS)		
C      GAM2=Z*Z*IH/KT  WHERE Z IS TARGET ION CHARGE			
C               AND IH IS THE RYDBERG ENERGY				
C  OUTPUT								
C      GIIIAV=MAXWELLIAN AVERAGED FREE-FREE GAUNT FACTOR.		
       integer i
       real*8 gam2,rgam2,sum,u,x,e1,e2,giii

       REAL*8 XA(8),WA(8)						
       DATA XA/0.17027963D0,0.90370178D0,2.25108663D0,4.26670017D0,	
     &7.04590540D0,10.75851601D0,15.74067864D0,22.86313174D0/		
       DATA WA/3.69188589D-1,4.18786781D-1,1.75794987D-1,		
     &3.33434923D-2,2.79453624D-3,9.07650877D-5,8.48574672D-7,		
     &1.04800117D-9/							
       RGAM2=1.0D0/GAM2							
       SUM=0.0D0							
       DO 20 I=1,8							
       X=XA(I)								
       E1=RGAM2*X							
       E2=RGAM2*(X+U)							
   20  SUM=SUM+WA(I)*GIII(1,1,E1,E2)					
       GIIIAV=SUM							
       RETURN								
      END								
C									
C									
C     IPP/01 - Krieger: fixed function declaration
c     NOTE: Old syntax was acceptable F90 for some compilers       
C									
      REAL*8 FUNCTION GIII(JZ,L,E1,E2)                                  
      IMPLICIT REAL*8(A-H,O-Z)                                          
C CALCULATES GIII GIVEN IN EQUATIONS (11) AND (15) OF A. BURGESS,       
C J. PHYS. B7,?,1974. SET L=1.SET JZ ZERO FOR THE ZERO CHARGE (NEUTRAL  
C ATOM) CASE. FOR PARTIAL SUMS SET L TO LOWER LIMIT OF SUMMATION(MUST BE
C GREATER THAN ZERO). SET E1=(KAPPA1)**2 FOR NON ZERO CHARGE, =(K1)**2  
C FOR ZERO CHARGE. SET E2=(KAPPA2)**2 FOR NON ZERO CHARGE,              
C =(K2)**2 FOR ZERO CHARGE.                                             
      integer jzl,l1,l2,jz,l
      real*8 e1,e2,f1,fdip,f2,dl,fdip0,el
      L1=L                                                              
      L2=L-1                                                            
      IF(JZ)1,2,1                                                       
1     F1=FDIP(E1,L1,E2,L2)                                              
      F2=FDIP(E1,L2,E2,L1)                                              
      EL=L                                                              
      GIII=1.102658D0*((1.0+EL*EL*E1)*F1*F1-(1.0+EL*EL*E2)*F2*F2)       
      GIII=GIII/(EL*(E1-E2))						
      RETURN                                                            
2     F1=FDIP0(E1,L1,E2,L2,1.0D-12)                                     
      F2=FDIP0(E1,L2,E2,L1,1.0D-12)                                     
      EL=L                                                              
      GIII=1.102658D0*EL*(E1*F1*F1-E2*F2*F2)/(E1-E2)                    
      RETURN                                                            
      END                                                               
C									
C									
C									
      REAL*8 FUNCTION FDIP(E1,L1,E2,L2)                                 
       IMPLICIT REAL*8(A-H,O-Z)                                         
C ALAN BURGESS DEPT. OF APPLIED MATHS. AND THEORETICAL PHYSICS,CAMBRIDGE
C CALCULATES THE FUNCTION I(KAPPA1,L1,KAPPA2,L2,1) DEFINED IN PHIL.     
C TRANS. ROY. SOC. A226,255,1970, WHERE E1=KAPPA1**2 AND E2=KAPPA2**2.  
C IT IS SUITABLE FOR USE IN EQUATIONS (8),(9),(10) OR (11) OF           
C J. PHYS. B. 7,L364,1974.                                              
      integer l1,l2
      real*8 e1,e2,emin,emax,t,fdip1,fdip2  
      IF(E1+E2-1.0D-40) 11,11,12                                        
   11 FDIP=0.0D0                                                        
      RETURN                                                            
   12 IF(E1-E2) 1,1,2                                                   
    1 EMIN=E1                                                           
      EMAX=E2                                                           
      GO TO 3                                                           
    2 EMIN=E2                                                           
      EMAX=E1                                                           
    3 T=EMIN/EMAX                                                       
      IF(T-0.02944D0) 4,4,5                                             
    4 FDIP=FDIP1(E1,L1,E2,L2)                                           
      GO TO 9                                                           
    5 IF(T-0.16667D0) 7,6,6                                             
    6 FDIP=FDIP2(E1,L1,E2,L2)                                           
      GO TO 9                                                           
    7 FDIP=FDIP1(E1,L1,E2,L2)                                           
       IF(FDIP*FDIP-1.0D-40) 6,6,8                                      
    8 RETURN                                                            
    9  IF(FDIP*FDIP-1.0D-40) 10,10,8                                    
   10 WRITE(6,100)                                                      
      RETURN                                                            
  100 FORMAT(15H   FDIP FAILURE)                                        
      END                                                               
C									
C									
C									
      REAL*8 FUNCTION FDIP0(E1,L1,E2,L2,EPS)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
C ALAN BURGESS,DEPT OF APPLIED MATHS. AND THEORETICAL PHYSICS,CAMBRIDGE 
C CALCULATES THE FUNCTION I0(K1,L1,K2,L2,1) DEFINED IN PHIL. TRANS.     
C ROY. SOC. A266,255,1970, WHERE E1=K1*K1, E2=K2*K2, AND THE RELATIVE   
C ACCURACY IS APPROXIMATELY EPS.                                        
C IT IS SUITABLE FOR USE IN EQUATIONS (13) ETC. OF J.PHYS.B. 7,L364,1974
      integer l1,l2,l,i,i0
      real*8 e1,e2,eps,el,e,p,p1,t,h0,h,ti,x,s,a,b,c,d,t1,f,f21,c1,
     >       ai,aii
 

c
      IF(L1-L2)1,2,4                                                    
1     L=L1                                                              
      GO TO 5                                                           
2     WRITE(6,100)L1                                                    
      FDIP0=0.0D0                                                       
3     RETURN                                                            
4     L=L2                                                              
5     EL=L                                                              
      FDIP0=0.5D0/(EL+1.0D0)                                            
      IF(E1-E2)6,3,7                                                    
6     E=E1/E2                                                           
      P=L1-L                                                            
      GO TO 8                                                           
7     E=E2/E1                                                           
      P=L2-L                                                            
8     FDIP0=FDIP0*E**((EL+P+0.5D0)*0.5D0)                               
C TO OBTAIN THE FUNCTION E1 OF M.J. SEATON, PROC. PHYS. SOC. A68,457,   
C 1955, REMOVE THE 'C' ON THE NEXT LINE.                                
C     FDIP0=1.0D0                                                       
      IF(E-0.5D0)21,20,20                                               
20    P1=P-0.5D0                                                        
      T=P1*(EL+1.0D0)*(E-1.0D0)                                         
      I0=L+1                                                            
      H0=0.0D0                                                          
      DO 9 I=1,I0                                                       
      TI=I                                                              
      H0=H0+1.0D0/TI                                                    
9     CONTINUE                                                          
      X=1.0D0-E                                                         
      H=1.0D0-(P+P+H0+DLOG(0.25D0*X))                                   
      S=1.0D0+T*H                                                       
      A=EL+1.0D0                                                        
      B=P1                                                              
      C=1.0D0                                                           
      D=0.0D0                                                           
10    A=A+1.0D0                                                         
      B=B+1.0D0                                                         
      C=C+1.0D0                                                         
      D=D+1.0D0                                                         
      T=T*A*B*X/(C*D)                                                   
      H=H+P1/(D*B)+EL/(C*A)                                             
      T1=T*H                                                            
      S=S+T1                                                            
      IF(DABS(T1)-EPS*DABS(S))13,11,11                                  
11    IF(C-300.0D0)10,12,12                                             
12    WRITE(6,101)                                                      
13    FDIP0=FDIP0*S                                                     
      RETURN                                                            
21    A=EL+1.0D0                                                        
      B=P-0.5D0                                                         
      C=EL+P+1.5D0                                                      
      F=F21(A,B,C,E,EPS)                                                
      L=L+1                                                             
      EL=L                                                              
      IF(P-0.5D0)23,23,24                                               
23    C1=EL+EL+1.0D0                                                    
      GO TO 25                                                          
24    C1=1.0D0                                                          
25    DO 22 I=1,L                                                       
      AI=I                                                              
      AII=AI+AI                                                         
      C1=C1*AI*AI*4.0D0/(AII*(AII+1.0D0))                               
22    CONTINUE                                                          
      FDIP0=FDIP0*F*C1                                                  
      RETURN                                                            
100   FORMAT(' FAILED IN FDIP0, L1=L2=',I5)                             
101   FORMAT(' FAILED TO CONVERGE IN FDIP0')                            
      END                                                               
C									
C									
C									
       REAL*8 FUNCTION FDIP1(E1,L1,E2,L2)                               
       IMPLICIT REAL*8(A-H,O-Z)                                         
       integer l1,l2,l,lp
       real*8 e1,e2,a1,a2,elp,b1,b2,fmon1
       IF(L1-L2)1,2,3                                                   
    1  L=L1                                                             
       A1=E1                                                            
       A2=E2                                                            
       GO TO 4                                                          
    2  FDIP1=0.0D0                                                      
       RETURN                                                           
    3  L=L2                                                             
       A1=E2                                                            
       A2=E1                                                            
    4  LP=L+1                                                           
       ELP=LP                                                           
       B1=DSQRT(1.0D0+ELP*ELP*A2)*FMON1(E1,E2,L)                        
       B2=DSQRT(1.0D0+ELP*ELP*A1)*FMON1(E1,E2,LP)                       
       IF(B1*B2-1.0D-40)5,5,6                                           
    5  FDIP1=0.0D0                                                      
       RETURN                                                           
    6  FDIP1=(B1-B2)/ELP                                                
       RETURN                                                           
      END                                                               
C									
C									
C									
       REAL*8 FUNCTION FDIP2(E1,L1,E2,L2)                               
       IMPLICIT REAL*8(A-H,O-Z)                                         
       integer l1,l2,l,j1,j2,iw0
       real*8 e1,e2,wmax,eta1,eta2,w1,pi,a,b,c,c1,c2,t1,u0,u1,
     >        v0,v1,w0,x0,y2,y0,y1,z0,z1,t,q0,q1,x,t0,p,argam,p0,p1

       WMAX=200.0D0                                                     
       ETA1=1.0D0/DSQRT(E1)                                             
       ETA2=1.0D0/DSQRT(E2)                                             
       W1=ETA2-ETA1                                                     
       PI=3.141592653589793D0                                           
       A=DABS(W1)                                                       
       B=PI*A                                                           
       IF(B-0.01D0)1,1,2                                                
    1  C=3.0D0/(3.0D0-B*(3.0D0-B*(2.0D0-B)))                            
       C=DSQRT(C)                                                       
       GO TO 5                                                          
    2  IF(B-14.0D0)4,3,3                                                
    3  C=DSQRT(B+B)                                                     
       GO TO 5                                                          
    4  B=B+B                                                            
       C1=1.0D0-DEXP(-B)                                                
       C=DSQRT(B/C1)                                                    
    5  C=0.5D0*C/DSQRT(ETA1*ETA2)                                       
       C2=ETA1+ETA2                                                     
       C1=4.0D0*ETA1*ETA2/(C2*C2)                                       
       L=L1                                                             
       IF(L2-L1)6,6,7                                                   
    6  L=L2                                                             
       T1=ETA1                                                          
       ETA1=ETA2                                                        
       ETA2=T1                                                          
       W1=-W1                                                           
    7  C=C*C1**(L+1)                                                    
       U0=L+1                                                           
       U1=ETA1                                                          
       V0=U0                                                            
       V1=-ETA2                                                         
       W0=1.0D0                                                         
       X0=W1/(C2*C2)                                                    
       Y2=-ETA2-ETA2                                                    
       Y0=-U0*W1+Y2                                                     
       Y1=ETA2*W1                                                       
       T1=X0/(1.0D0+W1*W1)                                              
       Z0=U0*T1                                                         
       Z1=U1*T1                                                         
       T=Z0-Z1*W1                                                       
       Z1=Z0*W1+Z1                                                      
       Z0=T                                                             
       Q0=-1.0D0+Z0*Y0-Z1*Y1                                            
       Q1=Z0*Y1+Z1*Y0                                                   
       X=W1*X0                                                          
    8  U0=U0+1.0D0                                                      
       V0=V0+1.0D0                                                      
       W0=W0+1.0D0                                                      
       IF(W0-WMAX)21,21,20                                              
   20  FDIP2=0.0D0                                                      
       RETURN                                                           
   21  CONTINUE                                                         
       Y0=Y0+Y2                                                         
       T=Z0*U0-Z1*U1                                                    
       Z1=Z0*U1+Z1*U0                                                   
       Z0=T                                                             
       T=Z0*V0-Z1*V1                                                    
       Z1=Z0*V1+Z1*V0                                                   
       Z0=T                                                             
       T=Z0*W0-Z1*W1                                                    
       Z1=Z0*W1+Z1*W0                                                   
       Z0=T                                                             
       X0=X/(W0*(W0*W0+W1*W1))                                          
       Z0=Z0*X0                                                         
       Z1=Z1*X0                                                         
       T0=Z0*Y0-Z1*Y1                                                   
       T1=Z0*Y1+Z1*Y0                                                   
       Q0=Q0+T0                                                         
       Q1=Q1+T1                                                         
       T1=T0*T0+T1*T1                                                   
       T0=Q0*Q0+Q1*Q1                                                   
       IF(T0-1.0D24*T1)8,8,9                                            
    9  J1=0                                                             
       J2=L+1                                                           
       P=ARGAM(J1,W1)+ARGAM(L,ETA1)-ARGAM(J2,ETA2)                      
       IW0=W0                                                           
       IF(A-1.0D-40)11,11,10                                            
   10  P=P+W1*DLOG(C2/A)                                                
   11  P0=DCOS(P)                                                       
       P1=DSIN(P)                                                       
       T=P0*Q0-P1*Q1                                                    
       Q1=P0*Q1+P1*Q0                                                   
       Q0=T                                                             
       FDIP2=C*Q1                                                       
       RETURN                                                           
      END                                                               
C									
C									
C									
       REAL*8 FUNCTION FMON1(E1,E2,L)                                   
       IMPLICIT REAL*8(A-H,O-Z)                                         
       integer  l,iv
       real*8 e1,e2, vmax,x1,x2,x3,x4,x5,x6,x7,pi,eta,g,a1,a2,mg,
     >        ma1,ma2,m,em,t,emn,b,s0,s1,u,v,w,t0,t1,u0,u1,s,emm

       IF(E1+E2-1.0D-40)28,28,29                                        
   28  FMON1=1.0D50                                                     
       RETURN                                                           
   29  CONTINUE                                                         
       VMAX=200.0D0                                                     
       X1=DSQRT(E1)                                                     
       X2=DSQRT(E2)                                                     
       X3=X1+X2                                                         
       X4=X3*X3                                                         
       X5=X1*X2                                                         
       X6=X2-X1                                                         
       X7=4.0D0/X4                                                      
       PI=3.141592653589793D0                                           
       IF(E1-E2)1,1,2                                                   
    1  ETA=1.0D0/X2                                                     
       GO TO 3                                                          
    2  ETA=1.0D0/X1                                                     
    3  G=0.5D0*PI*DEXP(-PI*ETA)                                         
       A1=1.0D0                                                         
       A2=1.0D0                                                         
       MG=0                                                             
       MA1=0                                                            
       MA2=0                                                            
       M=-1                                                             
    4  M=M+1                                                            
       EM=M                                                             
       T=EM+EM+1.0D0                                                    
       G=G*X7/(T*(T+1.0D0))                                             
       EMM=EM*EM                                                        
       A1=A1*(1.0D0+EMM*E1)                                             
       A2=A2*(1.0D0+EMM*E2)                                             
   30  IF(G-0.015625D0) 31,32,32                                        
   31  G=64.0D0*G                                                       
       MG=MG-1                                                          
       GO TO 30                                                         
   32  IF(G-64.0D0) 34,34,33                                            
   33  G=0.015625D0*G                                                   
       MG=MG+1                                                          
       GO TO 32                                                         
   34  IF(A1-64.0D0) 36,36,35                                           
   35  A1=0.015625D0*A1                                                 
       MA1=MA1+1                                                        
       GO TO 34                                                         
   36  IF(A2-64.0D0) 38,38,37                                           
   37  A2=0.015625D0*A2                                                 
       MA2=MA2+1                                                        
       GO TO 36                                                         
   38  CONTINUE                                                         
       IF(M-L)4,5,5                                                     
    5  G=G*(T+1.0D0)                                                    
       IF(X1-300.0D0)7,6,6                                              
    6  B=PI/X1                                                          
       A1=1.5D0*A1/(B*(3.0D0-B*(3.0D0-B*(2.0D0-B))))                    
       GO TO 9                                                          
    7  IF(X1-0.2D0)9,9,8                                                
    8  B=-PI/X1                                                         
       A1=A1/(1.0D0-DEXP(B+B))                                          
    9  IF(X2-300.0D0)11,10,10                                           
   10  B=PI/X2                                                          
       A2=1.5D0*A2/(B*(3.0D0-B*(3.0D0-B*(2.0D0-B))))                    
       GO TO 13                                                         
   11  IF(X2-0.2)13,13,12                                               
   12  B=-PI/X2                                                         
       A2=A2/(1.0D0-DEXP(B+B))                                          
   13  G=G*DSQRT(A1*A2)*(8.0D0)**(MG+MG+MA1+MA2)                        
       S0=1.0D0                                                         
       S1=0.0D0                                                         
       U=L                                                              
       V=0.0D0                                                          
       W=U+U+1.0D0                                                      
       T0=1.0D0                                                         
       T1=0.0D0                                                         
   14  U=U+1.0D0                                                        
       V=V+1.0D0                                                        
       W=W+1.0D0                                                        
       IF(V-VMAX)21,21,20                                               
   20  FMON1=0.0D0                                                      
       RETURN                                                           
   21  CONTINUE                                                         
       U0=U*U*X5+1.0D0                                                  
       U1=U*X6                                                          
       T=T0*U0-T1*U1                                                    
       T1=T0*U1+T1*U0                                                   
       T0=T                                                             
       T=X7/(V*W)                                                       
       T0=T*T0                                                          
       T1=T*T1                                                          
       S0=S0+T0                                                         
       S1=S1+T1                                                         
       S=S0*S0+S1*S1                                                    
       T=T0*T0+T1*T1                                                    
CLDH   IF(S-1.0D24*T)14,15,15                                           
       IF(1.0D-24*S-T)14,15,15                                          
   15  FMON1=G*DSQRT(S)                                                 
       IV=V                                                             
       RETURN                                                           
      END                                                               
C									
C									
C									
      REAL*8 FUNCTION F21(A,B,C,D,EPS)                                  
      IMPLICIT REAL*8(A-H,O-Z)                                          
      integer i
      real*8 a,b,c,d,eps,t,dd,sum,tn1,ai,tn2,at,as 

      T=(A*B*D)/C                                                       
      DD=1.0D0/(1.0D0-D)                                                
      SUM=1.0D0+T                                                       
      TN1=0.0D0                                                         
      I=1                                                               
3     AI=I                                                              
      T=T*(A+AI)*(B+AI)*D/((C+AI)*(1.0D0+AI))                           
      TN2=T*DD                                                          
      F21=SUM+TN2                                                       
      SUM=SUM+T                                                         
      AT=DABS(T+TN2-TN1)                                                
      AS=DABS(F21)*EPS                                                  
      IF(AS-AT)1,2,2                                                    
1     TN1=TN2                                                           
      I=I+1                                                             
      IF(I-300)3,3,4                                                    
4     WRITE(6,100)                                                      
100   FORMAT(' FAILED TO CONVERGE IN F21')                              
2     RETURN                                                            
      END                                                               
C									
C									
C									
       REAL*8 FUNCTION ARGAM(L,A)                                       
C CALCULATES ARGGAMMA(L+1+I*A)                                          
C WHERE L IS AN INTEGER NOT LESS THAN ZERO                              
       IMPLICIT REAL*8(A-H,O-Z)                                         
       integer l,j0,j1,j
       real*8  a,b,c,d,z,d1,d2,d0,u,d3
       B=DABS(A)                                                        
       B=250.0*B**0.25-A*A                                              
       J0=L+1                                                           
       C=J0                                                             
       D=C*C                                                            
       Z=0.0                                                            
       IF(D-B)1,6,6                                                     
    1  B=DSQRT (B)                                                      
       J1=B                                                             
       DO 5 J=J0,J1                                                     
       D=J                                                              
       D=A/D                                                            
       D1=DABS(D)                                                       
       IF(D1-0.1)2,3,3                                                  
    2  D1=D*D                                                           
       D2=-35.0*D1+45.0                                                 
       D2=-D1*D2+63.0                                                   
       D2=-D1*D2+105.0                                                  
       D1=D-D*D1*D2/315.0                                               
       GO TO 4                                                          
    3  D1=DATAN (D)                                                     
    4  Z=Z+D1                                                           
    5  CONTINUE                                                         
       J0=J1+1                                                          
    6  D=J0                                                             
       D0=D*D                                                           
       U=A*A                                                            
       D1=1.0/(D0+U)                                                    
       D2=D1*D1                                                         
       D3=10.0*D0*D0-20.0*D0*U+2.0*U*U                                  
       D3=D3*D2-21.0*D0+7.0*U                                           
       D3=D3*D2+210.0                                                   
       D1=A*D3*D1/2520.0                                                
        ARGAM=-Z+0.5*A*DLOG(D0+U)+(D-0.5)*DATAN(A/D)-A-D1               
       RETURN                                                           
      END                                                               

