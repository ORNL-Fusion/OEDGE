*     ------------------------------------------------------------     *

      Subroutine Coulomb_Coll2(K11,K12,K12a,K12b,K13,KBetaf,D11,
     >D12,D13,Coulomb_log,Frictionf,CHIpara,Xper,VBG,TBG,TB,NB,
     >MB,MI,ZI,ZB)

*     ------------------------------------------------------------     *
*     This subroutine calculates the Trubnikov potentials and their    *
*     derivatives determining the two dimensional drift kinetic        *
*     Coulomb collision term.                                          *
*     The equation numbers below refer to                              *
*     "Improved Kinetic Test Particle Model for Inpurity Transport     *
*      in Tokamaks" by D.Reiser, D.Reiter & M.Z.Tokar                  *
*     Nucl.Fusion Vol.38, No.2 (1998)                                  *
*     ------------------------------------------------------------     *
*     Input:                                                           
*
*     CHIpara = Normalized velocity chi parallel to B-field
*     MI      = Impurity mass [kg]                                     
*     MB      = Background ion mass [kg]
*     VBG     = Background ion velocity [m/s]     
*     TBG     = Background temperature gradient [eV/M] 
*     TB      = Background temperature [eV]
*     NB      = Background density [m^-3]
*     ZI      = Impurity charge [C]
*     ZB      = Background ion charge [C]
* 
*     ------------------------------------------------------------     *
*     Some Calculated Parameters:
*
*     Xper     = Normalized velocity chi perpendicular to B-field
*     XX       = SQRT(CHIpara^2 + CHIper^2)
*     C0       = Multiplying factor (usually = 1.0)                     
*     C1       = Multiplying factor                                    
*     C2       = Multiplying factor
*     kappa    = Ion heat conductivity = 1.0/0.5657          
*     Lambda1  = Coulomb-Logarithm  (usually = 15)          
*     Alpha1   = Inverse thermal velocity of background particles
*     Erf0     = Error function
*     Erf1     = First derivative of the error function      
*                                                                      
*     ------------------------------------------------------------     *
*     Output:
*
*     K11   = Frictional Acceleration Parallel to B-field:
*               
*                K11 = (1+MI/MB)*Lambda*Alpha^2*(CHIpara/CHI)*
*
*                      [2*Erf1/(CHI*SQRT(PI))-Erf0/CHI^2)]
*
*     Frictionf =  Frictional Acceleration Parallel to B-field using 
*                  the fluid approximation model:
*
*           Frictionf = -CHIpara*SQRT(2)*(1.0+MB/MI)*NB*lambda1*ZI^2*ZB^2
*                       -------------------------------------------------
*                                        1.5e9*MI*TB             
*
*     K12   = Thermal Acceleration Parallel to B-field:
*
*                K12 = BETAkinetic*TBG/MI
*
*           BETAkinetic = (1+MB/MI)*(3/2)*(ZI^2/ZB^2)*kappa*
*                           
*                              (1-2*CHIpara^2)*Erf1
*
*     KBetaf  = Acceleration Parallel to B-field where here BETAkinetic
*                is replaced with BETAfluid (BETAf):
* 
*           BETAfluid = -3[1-u-5*SQRT(2)*ZI^2*(1.1*u^(5/2)-0.35*u^(3/2)]
*                       --------------------------------------------------
*                                       2.6-2*u+5.4*u^2
*
*             where u = MI/(MI+MB)
*
*     K13   = Viscous Acceleration Parallel to B-field:
* 
*                K13 = (1+MI/MB)*Lambda*Alpha^2*(-2)*C2*(CHIpara/CHI)*
*
*                    [{-(4*CHI/3)-(2/CHI)-(3/CHI^3)}*Erf1+(3/CHI^4)*Erf0]
*
*     D11   = Parallel Diffusive term, dependent on NB
*     D12   = Parallel Diffusive term, dependent on TBG
*     D13   = Parallel Diffusive term, dependent on VBG
*
      IMPLICIT NONE

      INTEGER MI,MB,ZB,ZI,Xper
      REAL    VBG,TBG,TB,NB,K11,K12,K13,KBETAf,CHIpara
      REAL    D11,D12,D13,Frictionf,Coulomb_log
      REAL    K12a,K12b

      include 'params'
c      include 'cgeom'
c      include 'comtor' 

      INTEGER i
      
      REAL MII,MBB,ZBB,ZII,erxx
      REAL LAMBDA,ALPHA1,AlphaL,Alpha2L,C1,C2,Massfactor1,XX,Xb,u1
      REAL Erf0(0:5),Erf1(0:5),X(14)
      REAL Phi01,Phi11,Phi21,CosX,Cos2X,SinX,SIN2X,A01,A11,A12,A21,A22
      REAL Betaf
      REAL Psi01,Psi11,Psi21,C01,C11,C21
      REAL A11a,A11b

      real exp_coeff(8)

      REAL Erff

      EXTERNAL Erff

      D11 = 0.0
      D12 = 0.0
      D13 = 0.0
      Psi01 = 0.0
      Psi11 = 0.0
      Psi21 = 0.0
      C01   = 0.0
      C11   = 0.0
      C21   = 0.0

      MII     = REAL(MI)
      MBB     = REAL(MB)
      ZII     = REAL(ZI)
      ZBB     = REAL(ZB)
      ALPHA1  = SQRT(MBB/(2*TB*EMI))
c
c  Using the coulomb logarithm value of 12.5 from Reiser's thesis
c
c      LAMBDA  = 3.58616*ZII*ZII*NB*ZBB*ZBB*(12.5/15.0)/(MII*MII)
c      C1      = 1.20246e16*TB*TBG/(ZBB*ZBB*ZBB*ZBB*NB*12.5/15.0)
c      C2      = -4.9113E11*VBG*SQRT(MBB*TB*TB*TB)/
c     >          (ZBB*ZBB*ZBB*(1.2*ZBB + 0.8485)*NB*12.5/15.0)
c
c  Using the coulomb logarithm value of 15 as done normally in DIVIMP
c 
c      LAMBDA  = 3.58616*ZII*ZII*NB*ZBB*ZBB/(MII*MII)
c      C1      = 1.20246e16*TB*TBG/(ZBB*ZBB*ZBB*ZBB*NB)
c      C2      = -4.9113E11*VBG*SQRT(MBB*TB*TB*TB)/
c     >          (ZBB*ZBB*ZBB*(1.2*ZBB + 0.8485)*NB)

      LAMBDA  = 0.23995*Coulomb_log*ZII*ZII*NB*ZBB*ZBB/(MII*MII)
      C1      = 1.8024009e17*TB*TBG/(ZBB*ZBB*ZBB*ZBB*NB*Coulomb_log)
      C2      = -7.361216e12*VBG*SQRT(MBB*TB*TB*TB)/
     >          (ZBB*ZBB*ZBB*(1.2*ZBB + 0.8485)*NB*Coulomb_log)

*     ------------------------------------------------------------     *
*     Initialize some numerical quantities and direction angles
*     ------------------------------------------------------------     *

      MassFactor1 = 1.0 + MII/MBB
      u1          = MII/(MII+MBB)

      AlphaL	 = ALPHA1*Lambda
      Alpha2L	 = ALPHA1*AlphaL
 
*      CHIper = Xb
       if(Xper.eq.1)then
         Xb = CHIpara
       else
         Xb = 0.0
       endif
*      Xb     = Max(1.0E-4,Xb)		        ! Division by Xb below
    
*      XX = Sqrt(CHIpara*CHIpara+CHIper*CHIper) ! For FF+Fth comparison
*      XX = Sqrt(2*CHIpara*CHIpara)             ! When CHIpara = CHIper
 
      XX = Sqrt(CHIpara*CHIpara + Xb*Xb)

      Erf0(0) = Erff(XX)	              ! Error-Function
      Erf1(0) =(2.0/SQRT(Pi)) * Exp(-(XX*XX)) ! Derivative of Error-Function

      if (XX.ne.0.0) then
      	CosX  = CHIpara/XX
        Cos2X = CosX*CosX
      	SinX  = Xb/XX
        Sin2X = SinX*SinX
      else
        CosX  = 1.0
        Cos2X = 1.0
        SinX  = 0.0
        Sin2X = 0.0
      end if

*     ------------------------------------------------------------     *
*     Expansion coefficients for small chi limit
*     ------------------------------------------------------------     *
 
       EXP_COEFF(1) = 1.0
       EXP_COEFF(2) = 2.0/3.0
       EXP_COEFF(3) = 4.0/15.0
       EXP_COEFF(4) = 8.0/105.0
       EXP_COEFF(5) = 16.0/945.0
       EXP_COEFF(6) = 32.0/10395.0
       EXP_COEFF(7) = 64.0/135135.0
       EXP_COEFF(8) = 128.0/2027025.0

*     ------------------------------------------------------------     *
*     Computation for large chi
*     ------------------------------------------------------------     *

      if (XX.gt.1.0) then

	X(1) = XX
	do i=2,5
	  X(i) = X(i-1)*XX
	end do
	do i=1,5
	  Erf0(i) = Erf0(0)/X(i)
	  Erf1(i) = Erf1(0)/X(i)
	end do

        A01  = Erf1(1) - Erf0(2)
        A11  = C1*(Erf1(0) - 2.0*X(2)*Erf1(0))
        A12  = C1*Erf1(0)
	A21  = C2*(-(4.0*X(1)*Erf1(0)/3.0) - 2.0*Erf1(1) - 3.0*Erf1(3) 
     >	       + 3.0*Erf0(4))
	A22  = C2*(2.0*Erf1(1)/3.0 + Erf1(3) - Erf0(4))	

*       A11  = C1*((2*Erf1(0)) - 2.0*X(2)*Erf1(0)) ! Error in Reiser's Graph 	
     
        C01  = -Erf1(2) + Erf0(3)
	C11  = C1*(2.0*X(1)*Erf1(0) + 2.0*Erf1(1) + 3.0*Erf1(3) -
     >		3.0*Erf0(4))
	C21  = C2*(4.0*Erf1(0)/3.0 + 10.0*Erf1(2)/3.0 + 6.0*Erf1(4) +
     >		2.0*Erf0(3)/3.0 - 6.0*Erf0(5))

*     BREAKING APART THE COMPONENTS

        A11a = C1*Erf1(0)
        A11b = -C1*2.0*X(2)*Erf1(0)
       
      end if

*     ------------------------------------------------------------     *
*     Computation for small chi 
*     (the formulas in comments are for completeness and not necessary here)
*     ------------------------------------------------------------     *

      if (XX.le.1.0) then

	X(1) = XX
	do i=2,14
	  X(i) = X(i-1)*XX
	end do

	Phi01 = 0.0
	do i=0,6
	  Phi01 = Phi01 - EXP_COEFF(i+2)*X(2*i+1)
	end do
	Phi01 = Erf1(0)*Phi01

	Phi11 = Erf1(0)*(1.0-2.0*X(2))

	Phi21 = 0.0
	do i=0,5
	  Phi21 = Phi21 + 3.0*EXP_COEFF(i+3)*X(2*i+1)
	end do
        Phi21 = Erf1(0)*(-(4.0*X(1)/3.0)+Phi21)

	Psi01 = 0.0
	do i=0,5
	  Psi01 = Psi01 + EXP_COEFF(i+3)*X(2*i+2)
	end do
	Psi01 = Erf1(0)*(2.0/3.0+Psi01)

	Psi11 = 0.0
	do i=0,5
	  Psi11 = Psi11 - 3.0*EXP_COEFF(i+3)*X(2*i+1)
	end do
	Psi11 = Erf1(0)*(2.0*X(1) + Psi11)
        
        Psi21 = 0.0
	do i=0,4
	  Psi21 = Psi21 + (2.0*EXP_COEFF(i+3)/3.0 - 
     >            6.0*EXP_COEFF(i+4))*X(2*i+2)
	end do
	Psi21 = Erf1(0)*(8.0/45.0 + Psi21)


        A01  = Phi01
        A11  = C1*Phi11
        A12  = C1*Erf1(0)
	A21  = C2*Phi21
	A22  = 0.0
	do i=0,5
	  A22 = A22 - EXP_COEFF(i+3)*X(2*i+1)
	end do
	A22  = C2*Erf1(0)*A22

*        A11  = Phi11                     ! Error in Reiser's Graph
*        A11  = (Erf1(0) - 2.0*X(2)*Erf1(0))
*        A11  = C1*((2*Erf1(0)) - 2.0*X(2)*Erf1(0))

        C01  = Psi01
        C11  = C1*Psi11
        C21  = C2*Psi21

*     BREAKING APART THE COMPONENTS

        A11a = C1*Erf1(0)
        A11b = -C1*2.0*X(2)*Erf1(0) 

      end if

*     ------------------------------------------------------------     *
*     Now sum up all parts of the drift and diffusion coefficients
*     ------------------------------------------------------------     *


      ERXX=ERFF(xx)
      IF (CHIpara.lt.0.0)then 
       ERXX = -erff(xx)
      ENDIF 

      K11 = Alpha2L * MassFactor1 * CosX*A01 

      Frictionf = -CHIpara*SQRT(2.0)*(1.0+(MBB/MII))*NB*15.0*
     >            ZII*ZII*ZBB*ZBB/(1.5e9*MII*TB)

      K12 = Alpha2L * MassFactor1 * ((Cos2X * A11)+(Sin2X * A12))

      K12a = Alpha2L * MassFactor1 * Cos2X * A11a
      K12b = Alpha2L * MassFactor1 * Cos2X * A11b

      K13 = Alpha2L * MassFactor1 * ((1.0-3.0*Cos2X)*CosX*A21 -
     >                                6.0*CosX*SinX*CosX*A22)

      BETAf = -3.0*(1.0-u1-(5.0*REAL(ZI*ZI)*SQRT(2.0))*
     >             (1.1*(u1**(2.5))-0.35*(u1**(1.5))))/
     >             (2.6 - (2*u1) + (5.4*u1*u1))

      KBetaf = (BETAf*TBG*EMI/MII)

      D11 = AlphaL*Cos2X*C01
      D12 = AlphaL*CosX*Cos2X*C11
      D13 = AlphaL*(1.0-3.0*Cos2x)*Cos2X*C21     
 
      RETURN
      END

*     ------------------------------------------------------------     *

      REAL Function Erff(X)

*     ------------------------------------------------------------     *
*     --Error-Function for 0 <= x <= infinity --------------------     *
*     --error less than 2.5E-5 for all values of x ---------------     *
*     --Abramowicz 7.1.25 ----------------------------------------     *
*     ------------------------------------------------------------     *

*      IMPLICIT NONE
*      REAL X,A1,A2,A3,T,T2,T3,P
*      DATA A1/0.34802/, A2/-0.09587/
*      DATA A3/0.74785/, P/0.47047/
*      T    =  1.0/(1.0 + P*X)
*      T2   =  T*T
*      T3   =  T*T2
*      Erff = 1.0 - (A1*T+A2*T2+A3*T3)*Exp(-(X*X))

      IMPLICIT NONE
      REAL X,A1,A2,A3,A4,A5,T,T2,T3,T4,T5,P
      DATA A1/0.254829592/, A2/-0.284496736/
      DATA A3/1.421413741/, A4/-1.453152027/
      DATA A5/1.061405429/, P/0.3275911/
      T    = 1.0/(1.0 + P*X)      
      T2   = T*T
      T3   = T*T2
      T4   = T*T3
      T5   = T*T4
      Erff = 1.0 - (A1*T+A2*T2+A3*T3+A4*T4+A5*T5)*Exp(-(X*X))

      return
      end

*     ------------------------------------------------------------     *

*      Program Root
*      real a,b,tol,g,bis,answer
*      real arr(100),rt(10),i
*      integer nvalue,pvalue,j,k
*      external g
*      j=0
*      k=0
*      nvalue = 0
*      pvalue = 0
*      write(0,*)'Enter the interval and the required tolerance'
*      read*,a,b,tol
*      answer=bis(a,b,tol,g)
*      print*,'Root is',answer,' tolerance',tol
*      print*,cos(a),cos(b),a,b,tol
*      do i=-2.0,5.0,0.1
*         k=k+1
*         if (nvalue.eq.-1.and.pvalue.eq.1) then
*           j=j+1
*           print*,'j',j,arr(k-2),arr(k-1)
*           rt(j) = bis(i-0.2,i-0.1,tol,g)
*           pvalue = 0
*           nvalue = 0
*         endif
*         arr(k) = g(i)
*         print*,'i',i,arr(k),nvalue,pvalue
*         if (arr(k).gt.0.0) then
*            pvalue = 1
*         else
*            nvalue = -1
*         endif
*      end do
*      write(0,*)rt
*      stop
*      end

      Real Function BISCT(a,b,vbg,tbg,tb,nb,mb,mi,zi,zb,xper,force)

      IMPLICIT NONE
      real a,b,tolerance,mid,Ya,Yb,Ymid
      real vbg,tbg,tb,nb,kaa,kab,kaba,kabb,kac,kbetaf
      real Daa,Dab,Dac,Frictionf,Coulomb_log
      integer counter,dcounter,mb,mi,zi,zb,xper,force
      counter=0
      dcounter=0
      tolerance=0.0001
      mid=(a+b)/2
     
    5 if(abs(b-a).ge.tolerance)then
       CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >    Daa,Dab,Dac,Coulomb_log,Frictionf,a,Xper,
     >                    VBG,TBG,TB,NB,MB,MI,ZI,ZB)

       if(force.eq.1) then
        Ya=Kaa*REAL(MI)*1.67e-27
        CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >    Daa,Dab,Dac,Coulomb_log,Frictionf,a,Xper,
     >                     VBG,TBG,TB,NB,MB,MI,ZI,ZB)
        Ymid = Kaa*REAL(MI)*1.67e-27

       elseif(force.eq.2) then
        Ya=Kab*REAL(MI)*1.67e-27
        CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >    Daa,Dab,Dac,Coulomb_log,Frictionf,a,Xper,
     >                     VBG,TBG,TB,NB,MB,MI,ZI,ZB)
        Ymid = Kab*REAL(MI)*1.67e-27

       elseif(force.eq.3) then
        Ya=Kac*REAL(MI)*1.67e-27
        CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >    Daa,Dab,Dac,Coulomb_log,Frictionf,a,Xper,
     >                     VBG,TBG,TB,NB,MB,MI,ZI,ZB)
        Ymid = Kac*REAL(MI)*1.67e-27

       elseif(force.eq.4) then
        Ya=(Kaa+Kab)*REAL(MI)*1.67e-27
*        Ya = cos(a)
        CALL Coulomb_Coll2(Kaa,Kab,Kaba,Kabb,Kac,KBetaf,
     >    Daa,Dab,Dac,Coulomb_log,Frictionf,a,Xper,
     >                     VBG,TBG,TB,NB,MB,MI,ZI,ZB)
        Ymid = (Kaa+Kab)*REAL(MI)*1.67e-27
*        Yb = cos(b)
*        Ymid = cos(mid)
*     WRITE(0,*)'a',a,'b',b,'mid',mid,'Ya',ya,'Yb',yb,'Ymid',ymid
       endif

        if (Ya*Ymid.lt.0.0)then
            b=mid
            counter=counter+1
            mid=(a+b)/2
            go to 5
        else if (Ya*Ymid.gt.0.0)then
            a=mid
            counter=counter+1
            mid=(a+b)/2
            go to 5
        else
            dcounter=dcounter+1
            if(Ya.ne.0.0.and.Ymid.ne.0.0)then
               bisct=b
*
*           jdemod - yb is not defined in any of this code and 
*                    does not appear to be used - so I have 
*                    commented out the following test.  
*
*            elseif(Yb.ne.0.0.and.Ymid.ne.0.0)then
*               bisct=a
            else
               dcounter=dcounter+1
               bisct=mid
*               write(0,*)Ya,Yb
            endif
        endif
      endif
*      WRITE(0,*)counter,dcounter
      if (dcounter.eq.0)then
         Bisct = mid
      endif
      return
      end
     
*      real function g(x)
*      real x
*      g=cos(x)
*      return
*      end  
 
*     ------------------------------------------------------------     *


      Subroutine Potentials(Xa,Xb,XX,Phi00,Phi01,Phi10,Phi11,Phi20,
     >   Phi21,Psi00,Psi01,Psi02,Psi10,Psi11,Psi12,Psi20,Psi21,Psi22)

*     ------------------------------------------------------------     *
*     This subroutine calculates the Trubnikov potentials and their    *
*     derivatives determining the two dimensional drift kinetic        *
*     Coulomb collision term.                                          *
*     The equation numbers below refer to                              *
*     "Improved Kinetic Test Particle Model for Inpurity Transport     *
*      in Tokamaks" by D.Reiser, D.Reiter & M.Z.Tokar                  *
*     Nucl.Fusion Vol.38, No.2 (1998)                                  *
*     ------------------------------------------------------------     *
*     Input:                                                           *
*                                                                      *
*     Xa       = Normalized velocity chi parallel to B-field           * eq. (36a)
*     Xb       = Normalized velocity chi perpendicular to B-field      * eq. (36a)
*     C0       = Multiplying factor (usually = 1.0)                    * 
*     C1       = Multiplying factor                                    * eq. (45a)
*     C2       = Multiplying factor                                    * eq. (45b)
*     Lambda   = Coulomb-Logarithm                                     * 
*     Alpha    = Inverse thermal velocity of background particles      * eq. (36b)
*     Mit      = Impurity mass [kg]                                    *
*     Mib      = Background ion mass [kg]                              *
*                                                                      *
*     ------------------------------------------------------------     *
*     Output:                                                          *
*                                                                      *
*     Koeff(1) = K_parallel                                            * eq. (32a)
*     Koeff(2) = K_perp                                                * eq. (32b)
*     Koeff(3) = D_parallel,parallel                                   * eq. (32c)
*     Koeff(4) = D_parallel,perp                                       * eq. (32d)
*     Koeff(5) = D_perp,parallel                                       * eq. (32e)
*     Koeff(6) = D_perp,perp                                           * eq. (32f)
*     DiffEW   = Eigenvalues of tensor                                 * 
*     DiffEV   = Eigenvectors of tensor D                              * 
*     XX       = Normalized velocity chi                               *
*     Phi00    = Potential function phi0                               * eq. (46a) 
*     Phi01    = First derivative of phi0                              * eq. (48a) 
*     Phi10    = Potential function phi1                               * eq. (46b) 
*     Phi11    = First derivative of phi1                              * eq. (48b) 
*     Phi20    = Potential function phi2                               * eq. (46c) 
*     Phi21    = First derivative of phi2                              * eq. (48c) 
*     Psi00    = Potential function psi0                               * eq. (46d) 
*     Psi01    = First derivative of psi0                              * eq. (49a)
*     Psi02    = Second derivative of psi0                             * eq. (50a)
*     Psi10    = Potential function psi1                               * eq. (46e) 
*     Psi11    = First derivative of psi1                              * eq. (49b)
*     Psi12    = Second derivative of psi1                             * eq. (50b)
*     Psi20    = Potential function psi2                               * eq. (46f) 
*     Psi21    = First derivative of psi2                              * eq. (49c)
*     Psi22    = Second derivative of psi2                             * eq. (50c)
*     A01      = phi01                -> matrix component (A0)_1       * eq. (53a)
*     A11      = phi11                -> matrix component (A1)_1       * eq. (54a)
*     A12      = phi10/chi            -> matrix component (A1)_2       * eq. (54b)
*     A21      = phi21/chi            -> matrix component (A2)_1       * eq. (55a) 
*     A22      = phi20/chi            -> matrix component (A2)_2       * eq. (55b)
*     B01      = psi01                -> matrix component (B0)_1       * eq. (56a)
*     B11      = psi11                -> matrix component (B1)_1       * eq. (57a)
*     B12      = psi10/chi            -> matrix component (B1)_2       * eq. (57b)
*     B21      = psi21                -> matrix component (B2)_1       * eq. (58a)
*     B22      = psi20/chi            -> matrix component (B2)_2       * eq. (58b)
*     C01      = psi02                -> matrix component (C0)_11      * eq. (59a)
*     C02      = psi01/chi            -> matrix component (C0)_22      * eq. (59c)
*     C11      = psi12                -> matrix component (C1)_11      * eq. (60a)
*     C12      = psi11/chi-psi10/chi  ->  (C1)_12 & (C1)_22            * eqs. (60b,60c)
*     C21      = psi22                -> matrix component (C2)_11      * eq. (61a)
*     C22      = psi21/chi            ->  (C2)_12 & (C2)_22            * eqs. (61b,61c)
*     C23      = psi20/chi^2          ->  (C2)_12 & (C2)_22            * eqs. (61b,61c) 
*     D1Phi(1) = Derivative of phi with respect to chi_parallel        * eq. (51a)
*     D1Phi(2) = Derivative of phi with respect to chi_perp            * eq. (51a)
*     D1Psi(1) = Derivative of psi with respect to chi_parallel        * eq. (51b)
*     D1Psi(2) = Derivative of psi with respect to chi_perp            * eq. (51b)
*     D2Psi(1,1) = Second derivative of psi with respect to            *
*                  chi_parallel                                        * eq. (51c)
*     D2Psi(1,2) = Second derivative of psi with respect to            *
*                  chi_parallel and chi_perp                           * eq. (51c)
*     D2Psi(2,1) = Second derivative of psi with respect to            *
*                  chi_perp and chi_parallel                           * eq. (51c)
*     D2Psi(2,2) = Second derivative of psi with respect to            *
*                  chi_perp                                            * eq. (51c)
*                                                                      *
*     ------------------------------------------------------------     *
*     Used Functions: Erff                                              *
*     ------------------------------------------------------------     *

      IMPLICIT NONE

      INTEGER i
      REAL Xa,Xb
      REAL CosX,Cos2X,SinX,Sin2X
      REAL Erf0(0:5),Erf1(0:5),X(14),A(8)
      REAL XX,Phi00,Phi01,Phi10,Phi11,Phi20,Phi21
      REAL Psi00,Psi01,Psi02,Psi10,Psi11,Psi12,Psi20,Psi21,Psi22
      REAL Erff

      EXTERNAL Erff

      REAL Pi
      PARAMETER( Pi=3.1415926 )
      REAL SqrtPi
      PARAMETER( SqrtPi=1.77724539)

*     ------------------------------------------------------------     *
*     Expansion coefficients for small chi limit
*     ------------------------------------------------------------     *

      A(1) = 1.0
      A(2) = 2.0/3.0
      A(3) = 4.0/15.0
      A(4) = 8.0/105.0
      A(5) = 16.0/945.0
      A(6) = 32.0/10395.0
      A(7) = 64.0/135135.0
      A(8) = 128.0/2027025.0

*     ------------------------------------------------------------     *
*     Initialize some numerical quantities and direction angles
*     ------------------------------------------------------------     *

      Xb = Max(1.0E-4,Xb)		      ! Division by Xb below

      XX = Sqrt(Xa*Xa+Xb*Xb)
      Erf0(0) = Erff(XX)	              ! Error-Function
      Erf1(0) = 2.0/SqrtPi * Exp(-(XX*XX))    ! Derivative of Error-Function
      if (XX.ne.0.0) then
	CosX  = Xa/XX
        Cos2X = CosX*CosX
	SinX  = Xb/XX
        Sin2X = SinX*SinX
      else
        CosX  = 1.0
        Cos2X = 1.0
        SinX  = 0.0
        Sin2X = 0.0
      end if

*     ------------------------------------------------------------     *
*     Computation for large chi
*     ------------------------------------------------------------     *

      if (XX.gt.1.0) then       

	X(1) = XX
	do i=2,5
	  X(i) = X(i-1)*XX
	end do
	do i=1,5
	  Erf0(i) = Erf0(0)/X(i)
	  Erf1(i) = Erf1(0)/X(i)
	end do

	Phi00 = Erf0(1)
	Phi01 = Erf1(1) - Erf0(2)
	Phi10 = X(1)*Erf1(0)
	Phi11 = Erf1(0) - 2.0*X(2)*Erf1(0)
	Phi20 = 2.0*Erf1(0)/3.0 + Erf1(2) - Erf0(3)
	Phi21 = -4.0*X(1)*Erf1(0)/3.0 - 2.0*Erf1(1) - 3.0*Erf1(3) +
     >		 3.0*Erf0(4)
	Psi00 = Erf1(0)/2.0 + Erf0(1)/2.0 + X(1)*Erf0(0)
	Psi01 = Erf1(1)/2.0 + Erf0(0) - Erf0(2)/2.0
	Psi02 = -Erf1(2) + Erf0(3)
	Psi10 = Erf1(1)/2.0 - Erf0(2)/2.0
	Psi11 = -Erf1(0) - Erf1(2) + Erf0(3)
	Psi12 = 2.0*X(1)*Erf1(0) + 2.0*Erf1(1) + 3.0*Erf1(3) -
     >		 3.0*Erf0(4)
	Psi20 = Erf1(2)/2.0 + Erf0(1)/3.0 - Erf0(3)/2.0
	Psi21 = -2.0*Erf1(1)/3.0 - 3.0*Erf1(3)/2.0 - Erf0(2)/3.0 +
     >		 3.0*Erf0(4)/2.0
	Psi22 = 4.0*Erf1(0)/3.0	+ 10.0*Erf1(2)/3.0 + 6.0*Erf1(4) +
     >		 2.0*Erf0(3)/3.0 - 6.0*Erf0(5)

      end if

*     ------------------------------------------------------------     *
*     Computation for small chi 
*     (the formulas in comments are for completeness and not necessary here)
*     ------------------------------------------------------------     *


      if (XX.le.1.0) then

	X(1) = XX
	do i=2,14
	  X(i) = X(i-1)*XX
	end do

	Phi00 = 0.0
	do i=0,6
	  Phi00 = Phi00 + A(i+2)*X(2*i+2)
	end do
	Phi00 = Erf1(0)*(1.0 + Phi00)
	Phi01 = 0.0
	do i=0,6
	  Phi01 = Phi01 - A(i+2)*X(2*i+1)
	end do
	Phi01 = Erf1(0)*Phi01
	Phi10 = Erf1(0)*X(1)
	Phi11 = Erf1(0)*(1.0-2.0*X(2))
	Phi20 = 0.0
	do i=0,5
	  Phi20 = Phi20 - A(i+3)*X(2*i+2)
	end do
	Phi20 = Erf1(0)*Phi20
	Phi21 = 0.0
	do i=0,5
	  Phi21 = Phi21 + 3.0*A(i+3)*X(2*i+1)
	end do
	Phi21 = Erf1(0)*(-(4.0*X(1)/3.0)+Phi21)
	Psi00 = 0.0
	do i=0,6
	  Psi00 = Psi00 + (A(i+2)/2.0+A(i+1))*X(2*i+2)
	end do
	Psi00 = Erf1(0)*(1.0 + Psi00)
	Psi01 = 0.0
	do i=0,6
	  Psi01 = Psi01 + (A(i+1)-A(i+2)/2.0)*X(2*i+1)
	end do
	Psi01 = Erf1(0)*Psi01
	Psi02 = 0.0
	do i=0,5
	  Psi02 = Psi02 + A(i+3)*X(2*i+2)
	end do
	Psi02 = Erf1(0)*(2.0/3.0+Psi02)
	Psi10 = 0.0
	do i=0,6
	  Psi10 = Psi10 - A(i+2)*X(2*i+1)/2.0
	end do
	Psi10 = Erf1(0)*Psi10
	Psi11 = 0.0
	do i=0,5
	  Psi11 = Psi11 + A(i+3)*X(2*i+2)
	end do
	Psi11 = Erf1(0)*(-(1.0/3.0) + Psi11)
	Psi12 = 0.0
	do i=0,5
	  Psi12 = Psi12 - 3.0*A(i+3)*X(2*i+1)
	end do
	Psi12 = Erf1(0)*(2.0*X(1) + Psi12)
	Psi20 = 0.0
	do i=0,5
	  Psi20 = Psi20 + (A(i+2)/3.0 - A(i+3)/2.0)*X(2*i+2)
	end do
	Psi20 = Erf1(0)*Psi20
	Psi21 = 0.0
	do i=0,5
	  Psi21 = Psi21 + (-(A(i+2)/3.0) + 3.0*A(i+3)/2.0)*X(2*i+1)
	end do
	Psi21 = Erf1(0)*Psi21
	Psi22 = 0.0
	do i=0,4
	  Psi22 = Psi22 + (2.0*A(i+3)/3.0 - 6.0*A(i+4))*X(2*i+2)
	end do
	Psi22 = Erf1(0)*(8.0/45.0 + Psi22)

      end if

      RETURN
      END

*     ------------------------------------------------------------     *

      Subroutine COEFF_C1(C1plus,C1minus,Xa,Xb,X,
     >                    Psi01,Psi02,Psi10,Psi11,Psi12)

      IMPLICIT NONE

      REAL Xa,Xb,X,Psi01,Psi02,Psi10,Psi11,Psi12
      REAL aa,bb,cc,C1plus,C1minus

      aa = Psi02*Psi01/X

      bb = (Xa/X)*(Psi12*Psi01/X + Psi02*(Psi11/X - 
     >     Psi10/(X*X)))

      cc = Xa*Xa/(X*X)*Psi12*(Psi11/X - Psi10/(X*X)) - 
     >     Xb*Xb/(X*X)*((Psi11/X - Psi10/(X*X))**2)

      C1plus = (-bb + SQRT(bb*bb - 4.0*aa*cc))/(2.0*cc)
      C1minus = (-bb - SQRT(bb*bb - 4.0*aa*cc))/(2.0*cc)

      RETURN
      END

*     ------------------------------------------------------------     *      
      
      Subroutine COEFF_C2(C2plus,C2minus,Xa,Xb,X,
     >                    Psi01,Psi02,Psi20,Psi21,Psi22)

      IMPLICIT NONE

      REAL Xa,Xb,X,Psi01,Psi02,Psi20,Psi21,Psi22
      REAL aa,bb,cc,C2plus,C2minus

      aa = Psi02*Psi01/X

      bb = (1.0 - 3.0*Xa*Xa/(X*X))*Psi02*Psi21/X -
     >     6.0*(1.0 - 2.0*Xa*Xa/(X*X))*Psi02*Psi20/(X*X) +
     >     (1.0 - 3.0*Xa*Xa/(X*X))*Psi22*Psi01/X

      cc = ((1.0 - 3.0*Xa*Xa/(X*X))**2)*Psi22*Psi21/X -
     >     6.0*(1.0 - 3.0*Xa*Xa/(X*X))*(1.0 - 2.0*Xa*Xa/(X*X))*
     >     Psi22*Psi20/(X*X) - 36.0*Xa*Xa*Xb*Xb/(X*X*X*X)*
     >     (Psi21/X - Psi20/(X*X))*(Psi21/X - Psi20/(X*X))

      C2plus = (-bb + SQRT(bb*bb - 4.0*aa*cc))/(2.0*cc)
      C2minus = (-bb - SQRT(bb*bb - 4.0*aa*cc))/(2.0*cc)

      RETURN
      END

*     ------------------------------------------------------------      *
