      program test_thom
      implicit none
c
c     Program to test and write out Thompson distributions 
c
      integer maxpts,maxsets,mnorm,maxsamp
      parameter(maxpts=300, mnorm=6, maxsets=2*mnorm,maxsamp=100000)
c
      real*8 thom_data(maxpts+1,maxsets),eimp,ebd,maxe,sampe,deltae,maxv
      real*8 deltav,sampv
      real*8 get_thompson_dist,find_thompson_velocity
      external get_thompson_dist,find_thompson_velocity
      integer in,thom_opt,iv,ic
      integer seed(3)
      real*8 thom_max(mnorm)
      real*8 vf,ef,vout
      external vf,ef
      common /massf/ crmi
      real*8 crmi
c
      seed(1) = 1
      seed(2) = 1
      seed(3) = 1 
      call random_seed(put=seed)
c
      crmi = 12.0
      eimp = 10000.0d0
      ebd  = 7.62d0
      maxe = 150.0d0 
      maxv = vf(maxe)
      deltae = maxe/maxpts
      deltav = maxv/maxpts 
c
c     Sample data using E option - look at 0.0 to 5eV range
c
      do in  = 1,maxpts
c
         sampv = deltav * in - deltav/2.0
         sampe = ef(sampv)
c
         deltae = sampv * deltav 
c
         thom_data(in,1) = get_thompson_dist(eimp,ebd,sampe,0)
         thom_data(in,2) = get_thompson_dist(eimp,ebd,sampe,1)
         thom_data(in,5) = thom_data(in,1) * deltae
         thom_data(in,6) = thom_data(in,2) * deltav
c
      end do
c
c     Random sampling
c 
      do in = 1,maxsamp
c
         if (1000*int(in/1000.0).eq.in) then
            write(0,*) 'GENERATING:',in 
         endif
c
         vout = find_thompson_velocity(eimp,ebd,0) 

         iv = min(int(vout/deltav)+1,maxpts+1)

         thom_data(iv,3) = thom_data(iv,3) + 1.0

         vout = find_thompson_velocity(eimp,ebd,1) 
         
         iv = min(int(vout/deltav)+1,maxpts+1)

         thom_data(iv,4) = thom_data(iv,4) + 1.0
c
      end do 
c
c     Create normalized functions
c
      do ic = 1,mnorm
         thom_max(ic) = 0.0 
      end do 
c
      do in = 1,maxpts
c
         do ic = 1,mnorm
c
            thom_max(ic) = 
     >         max(thom_max(ic),thom_data(in,ic))
         end do
      end do
c
      do in = 1,maxpts
c
         do ic = 1,mnorm

            thom_data(in,ic+mnorm) = 
     >         thom_data(in,ic)/thom_max(ic)

         end do
c
      end do 
c
c     Write out results
c
      do in = 1,maxpts

         sampv = deltav * in - deltav/2.0
         sampe = ef(sampv)

         write(6,'(i6,20(1x,g16.5))') in,
     >         sampv,sampe,(thom_data(in,ic),ic=1,maxsets)

      end do 

      stop
      end
c 
c
c 
      real*8 function vf(e)
      implicit none
      real*8 e
      real*8 crmi 
      common /massf/ crmi
      vf = 1.38d4 * sqrt(e/crmi)
      return
      end
c
c
c
      real*8 function ef(v)
      implicit none
      real*8 v
      real*8 crmi 
      common /massf/ crmi
      ef =   crmi *  (v/1.38d4)**2 
      return
      end
c
c
c
      real function get_thompson_dist(eimpi,ebdi,einput,
     >                                     thom_opt)
      implicit none 
      integer thom_opt
      real*8 eimpi,ebdi,einput
c
c      include 'params'
c      include 'cgeom'
c      include 'comtor'
c      include 'pindata'
c
      common /thom_ye_params/  eimp,gamma,ebd
      real*8 eimp,gamma,ebd 

      common /thom_yv_params/  vimp,vgamma,vbd
      real*8 vimp,vgamma,vbd 

c
c     -Need a random number
c     -Need to calculate the integral over the distribution for this 
c      segment or have it stored in a pre-calculated array. 
c     -Option used for physical sputtering only.      
c     -What about the difference between atom and ion flux physical 
c      sputtering energies? These are part of Eimpact
c     -What is the value of gamma? gamma = 4 (mC*mD) / (mC+mD)**2
c     
c      intye_targ(maxnds)
c      intye_wall(maxpts)
c
c
c
c     Local variables  
c
      integer in,ik,ir,targ,ierr,iter_cnt
      real*8 emin,emax,e_result
      real*8 vmin,vmax,vtmp,v_result,vinput
      common /massf/ crmi  
      real*8   crmb,crmi 
c
      integer ringno
      real*8 thom_ye,thom_yv
      external thom_ye,thom_yv
c
      eimp = eimpi
      ebd = ebdi
c 
      crmb =  2.0
c     
      vimp = 1.38e4 * sqrt(eimp/crmi)    
      vinput = 1.38e4 * sqrt(einput/crmi)    
c
c     Calculate Gamma and vgamma
c
      GAMMA  = 4.0 * CRMB * CRMI / ((CRMB+CRMI) * (CRMB+CRMI))
      vgamma = gamma
c
c     Calculate Ebd and vbd
c
      vbd = 1.38e4 * sqrt(ebd/crmi)    
c
c     Calculate emin and emax limits for integration
c     Also vmin and vmax limits 
c     - note vmax requires a 
c       square root and instead of getting an error when negative
c       it would be imaginary - thus the sqrt is taken to get vmax
c       only after the quantity is verified. 
c
c
      emin = 0.0
      vmin = 0.0  
      emax = gamma*(1.0d0-gamma)*eimp-ebd
      vtmp = gamma*(1.0d0-gamma)*vimp**2-vbd**2
c
c     Check validity of Emax 
c
      if (thom_opt.eq.0.and.emax.lt.0.0) then 
c
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : GET_THOMPSON_VELOCITY : EMAX<0',
     >                  eimp,ebd,gamma,emax
         write(0,'(a,4(1x,g12.5))')
     >                  'ERROR : GET_THOMPSON_VELOCITY : EMAX<0',
     >                  eimp,ebd,gamma,emax
         stop 
      endif
c
c     Check validity of vmax
c
      if (thom_opt.eq.1.and.vtmp.lt.0.0) then 
c
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : GET_THOMPSON_VELOCITY :'//
     >                  ' VMAX IMAGINARY',
     >                  vimp,vbd,gamma,vtmp
         write(0,'(a,4(1x,g12.5))')
     >                  'ERROR : GET_THOMPSON_VELOCITY :'//
     >                  ' VMAX IMAGINARY',
     >                  vimp,vbd,gamma,vtmp
         stop 
c
      endif 
c
      vmax = sqrt(vtmp)
c
c      write(0,'(a,6(1x,g12.5))') 'FT:',emin,emax,ebd,vmin,vmax,vbd
c
c
c     Energy desired is last_e - calculate the velocity
c
      if (thom_opt.eq.0) then 
         get_thompson_dist = thom_ye(einput)
      elseif (thom_opt.eq.1) then   
         get_thompson_dist = thom_yv(vinput)
      endif
c
      return
      end


      real function find_thompson_velocity(eimpi,ebdi,
     >                                     thom_opt)
      implicit none 
      integer thom_opt
      real*8 eimpi,ebdi
c
c      include 'params'
c      include 'cgeom'
c      include 'comtor'
c      include 'pindata'
c
      common /thom_ye_params/  eimp,gamma,ebd
      real*8 eimp,gamma,ebd 

      common /thom_yv_params/  vimp,vgamma,vbd
      real*8 vimp,vgamma,vbd 

c
c     -Need a random number
c     -Need to calculate the integral over the distribution for this 
c      segment or have it stored in a pre-calculated array. 
c     -Option used for physical sputtering only.      
c     -What about the difference between atom and ion flux physical 
c      sputtering energies? These are part of Eimpact
c     -What is the value of gamma? gamma = 4 (mC*mD) / (mC+mD)**2
c     
c      intye_targ(maxnds)
c      intye_wall(maxpts)
c
c
c
c     Local variables  
c
      real*8 seed
      integer nrand

      integer in,ik,ir,targ,ierr,iter_cnt
      real*8 emin,emax,e_result
      real*8 vmin,vmax,vtmp,v_result
      common /massf/ crmi  
      real*8   crmb,crmi 
c
      integer ringno
      real*8 thom_ye,thom_yv
      external thom_ye,thom_yv
c
      eimp = eimpi
      ebd = ebdi
c
      crmb =  2.0
c     
      vimp = 1.38e4 * sqrt(eimp/crmi)    
c
c     Calculate Gamma and vgamma
c
      GAMMA  = 4.0 * CRMB * CRMI / ((CRMB+CRMI) * (CRMB+CRMI))
      vgamma = gamma
c
c     Calculate Ebd and vbd
c
      vbd = 1.38e4 * sqrt(ebd/crmi)    
c
c     Calculate emin and emax limits for integration
c     Also vmin and vmax limits 
c     - note vmax requires a 
c       square root and instead of getting an error when negative
c       it would be imaginary - thus the sqrt is taken to get vmax
c       only after the quantity is verified. 
c
c
      emin = 0.0
      vmin = 0.0  
      emax = gamma*(1.0d0-gamma)*eimp-ebd
      vtmp = gamma*(1.0d0-gamma)*vimp**2-vbd**2
c
c     Check validity of Emax 
c
      if (thom_opt.eq.0.and.emax.lt.0.0) then 
c
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : FIND_THOMPSON_VELOCITY : EMAX<0',
     >                  eimp,ebd,gamma,emax
         write(0,'(a,4(1x,g12.5))')
     >                  'ERROR : FIND_THOMPSON_VELOCITY : EMAX<0',
     >                  eimp,ebd,gamma,emax
         stop 
      endif
c
c     Check validity of vmax
c
      if (thom_opt.eq.1.and.vtmp.lt.0.0) then 
c
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : FIND_THOMPSON_VELOCITY :'//
     >                  ' VMAX IMAGINARY',
     >                  vimp,vbd,gamma,vtmp
         write(0,'(a,4(1x,g12.5))')
     >                  'ERROR : FIND_THOMPSON_VELOCITY :'//
     >                  ' VMAX IMAGINARY',
     >                  vimp,vbd,gamma,vtmp
         stop 
c
      endif 
c
      vmax = sqrt(vtmp)
c
c      write(0,'(a,6(1x,g12.5))') 'FT:',emin,emax,ebd,vmin,vmax,vbd
c
c     The next section of code is common two the various 
c     options available for the Thompson disribution. 
c
      if (thom_opt.eq.0) then  
         call evaluate_thom(thom_ye,emin,emax,e_result,
     >                      seed,nrand,ierr)

      elseif (thom_opt.eq.1) then 
         call evaluate_thom(thom_yv,vmin,vmax,v_result,
     >                      seed,nrand,ierr)
      endif
c
c     Energy desired is last_e - calculate the velocity
c
      if (thom_opt.eq.0) then 
         find_thompson_velocity = 1.38E4 * SQRT (e_result/CRMI) 
      elseif (thom_opt.eq.1) then   
         find_thompson_velocity = v_result
      endif
c
      return
      end
c
c
c
      subroutine evaluate_thom(thom_func,range_min,range_max,
     >                         result_val,seed,nrand,ierr)
      implicit none
      integer nrand,ierr
      real*8 thom_func,range_min,range_max,result_val,seed
      external thom_func
c
c     EVALUATE_THOM: 
c
c     This routine takes a range minimum and maxium as well as a function.
c     It then selects a random number in the range [0,1] and finds the 
c     result_value that corresponds to the value of the range for which 
c     the integration of the function over the interval is equal to that 
c     fraction of the integration for the full range. 
c 
c     Local variables
c
      real*8 norm_val      
      real ran
      real*8 ran_frac,res_frac,last_r,rtest,rmin,rmax
      real*8 cumulative_result,current_result
      integer iter_cnt
c
      real*8 eps
      parameter(eps=0.00001d0)
c
      call gen_qsimp(thom_func,range_min,range_max,norm_val,0)
c
c     Check validity of integration value
c
      if (norm_val.le.0.0) then 
         write(0,'(a,4(1x,g12.5))')
     >                  'ERROR : EVALUATE THOM : NORM <= 0',
     >                  range_min,range_max,norm_val
         write(6,'(a,4(1x,g12.5))')
     >                  'ERROR : EVALUATE THOM : NORM <= 0',
     >                  range_min,range_max,norm_val
         stop 
c
      endif 
c
c     Get random number for fraction of integration
c
      call random_number(ran)
      ran_frac=ran
c
c     Perform binary search on integral until fraction is obtained
c     within error limit - start with rtest = 1/2 (range_max - range_min)
c
      rmin = range_min
      rmax = range_max 
      rtest = 0.5d0 * (rmax-rmin)
c
      res_frac = -1.0
c
c     Perform search loop - include maximum iteration test - set at 10,000 for now
c
      iter_cnt = 0
      cumulative_result  = 0.0d0
c
      do while ((abs(res_frac-ran_frac).gt.eps).and.(iter_cnt.lt.10000))
c
         iter_cnt = iter_cnt + 1
c
         call gen_qsimp(thom_func,rmin,rtest,current_result,0)
c
         last_r   = rtest 
c
         res_frac = (cumulative_result+current_result) / norm_val 
c
c         write(6,'(a,i4,8(1x,g18.7))') 'QS:',iter_cnt,current_result,
c     >                 cumulative_result,
c     >                 res_frac,ran_frac,rmin,rtest,rmax
c         write(0,'(a,i4,8(1x,g18.7))') 'QS:',iter_cnt,current_result,
c     >                 cumulative_result,
c     >                 res_frac,ran_frac,rmin,rtest,rmax
c
c
c        Integral is less than the random fraction - need to increase the Energy
c        - add current integral to stored value
c        - set rmin to rtest
c        - set rtest to (rmin+rmax)/2
c        - integrate from rmin to rmax
c
         if (res_frac.lt.ran_frac) then 
c
            cumulative_result = cumulative_result + current_result 
            rmin = rtest
            rtest = 0.5d0*(rmin+rmax)
c
c        Integral is greater than the random fraction
c        - set rmax to rtest
c        - set rtest to (rmin + rmax)/2
c        - integrate from rmin to rmax
c
         else
c
            rmax = rtest
            rtest = 0.5d0*(rmin+rmax)
c
         endif
c
      end do
c
c     The value of the range for which the fractional condition is 
c     satisfied is stored in last_r
c
      result_val = last_r
c
c
c     return
      end
c
C
C
      SUBROUTINE GEN_QSIMP(FUNC,A,B,S,OPT)
      implicit none
      INTEGER JMAX,OPT,PLATE,ITYP
      DOUBLE PRECISION A,B,FUNC,S,EPS
      EXTERNAL FUNC
      PARAMETER (EPS=1.0D-7,JMAX=30)
C
C     GEN_QSIMP: THIS ROUTINE IS TAKEN FROM THE TEXT NUMERICAL RECIPES
C     IN FORTRAN. IT IMPLEMENTS SIMPSON'S METHOD OF CALCULATING
C     NUMERICAL QUADRATURE. IT CALLS THE ROUTINE TRAPZD TO REFINE
C     THE VALUE OF THE NUMERICAL INTEGRAL.
C
C     RETURNS AS S THE INTEGRAL OF THE FUNCTION FUNC FROM A TO B.
C     THE PARAMETER EPS CAN BE SET TO THE DESIRED FRACTIONAL ACCURACY
C     AND JMAX SO THAT 2 TO THE POWER JMAX-1 IS THE MAXIMUM ALLOWED NUMB
C     OF STEPS. INTEGRATION IS PERFORMED BY SIMPSONS RULE.
C
C     This is a generic version based on the code in soledge.d6a. The
c     main difference is that any specific function related parameters
c     are passed to the function using a common block so that the 
c     code itself can work for any function of the specified type. 
C
      INTEGER J
      DOUBLE PRECISION OS,OST,ST
c
c     Check for equal integration bounds 
c
      if (a.eq.b) then 
         s = 0.0
         return
      endif
c
      OST = -1.0D30
      OS  = -1.0D30
      DO 10 J = 1,JMAX
        CALL GEN_TRAPZD (FUNC,A,B,ST,J,OPT)

c
c       As a numerical check when dealing with small numbers
c       also check for s = os since if the numbers are small 
c       eps * abs(os) could underflow.  
c
        S = (4.0*ST-OST)/3.0
c        write(0,*) 'QSIMP:',s,os,s-os,eps*dabs(os)

        IF (ABS(S-OS).LE.(EPS*ABS(OS))) RETURN
C
C
        OS = S
        OST = ST
C
C        WRITE(6,*) 'QSIMP:',OST,OS,ST
C
10    CONTINUE
C
C     ERROR CONDITION - INCOMPLETE CONVERGENCE
C
      WRITE(6,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP: ERROR IN QUADRATURE',A,B,ST
      WRITE(6,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP:',OST,OS,ST ,S
      write(6,'(a,5(1x,g12.5))') 
     >        '    ALSO :',A,B,J,OPT
c
      WRITE(0,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP: ERROR IN QUADRATURE',A,B,ST
      WRITE(0,*) 
     >       'GEN_QSIMP: ERROR IN QUADRATURE',A,B,ST
c
      WRITE(0,'(a,5(1x,g12.5))') 
     >       'GEN_QSIMP:',OST,st,os ,S
      WRITE(0,*)
     >       'GEN_QSIMP:',OST,st,os ,S,s-os
c
      write(0,'(a,5(1x,g12.5))') 
     >        '    ALSO :',A,B,J,OPT
      STOP
C      RETURN
      END
C
C
C
      SUBROUTINE GEN_TRAPZD (FUNC,A,B,S,N,OPT)
      implicit none
      INTEGER N,OPT
      DOUBLE PRECISION A,B,S,FUNC
      EXTERNAL FUNC
C
C     TRAPZD: THIS SUBROUTINE IS TAKEN FROM THE TEXT NUMERICAL RECIPES
C     IN FORTRAN. IT IS PART OF A SET OF ROUTINES FOR CALCULATING
C     NUMERICAL QUADRATURE OF A DISCRETE FUNCTION.
C
C     THIS ROUTINE COMPUTES THE NTH STAGE REFINEMENT OF AN EXTENDED
C     TRAPEZOIDAL RULE. FUNC IS INPUT AS THE NAME OF THE FUNCTION TO
C     BE INTEGRATED BETWEEN THE LIMITS OF A AND B, ALSO INPUT. WHEN
C     CALLED WITH N=1, THE ROUTINE RETURNS THE CRUDEST ESTIMATE OF
C     THE INTEGRAL. SUBSEQUENT CALLS WITH N=2,3... (IN THAT
C     SEQUENTIAL ORDER) WILL IMPROVE THE ACCURACY OF S BY ADDING
C     2**(N-2) ADDITIONAL INTERIOR POINTS. S SHOULD NOT BE MODIFIED
C     BETWEEN SEQUENTIAL CALLS.
C
      INTEGER IT,J
      DOUBLE PRECISION DEL,SUM,TNM,X,test
c
      IF (N.EQ.1) THEN
c
         S = 0.5*(B-A)*(FUNC(A)+FUNC(B))
C
      ELSE
         IT = 2**(N-2)
         TNM = IT
         DEL = (B-A)/TNM
         X = A + 0.5*DEL
         SUM = 0.0
         DO 10 J = 1,IT
            SUM = SUM + FUNC(X)
            X = X + DEL
 10      CONTINUE
         S = 0.5*(S+(B-A)*SUM/TNM)
C
C         WRITE (6,*) 'TRAPZ:',S,SUM
C
      ENDIF
      RETURN
      END
c
c
c
      double precision function thom_ye(e)
      implicit none
      double precision e
      common /thom_ye_params/  eimp,gamma,ebd
      real*8 eimp,gamma,ebd 
c
c     THOM_YE: 
c
c     This function returns the value of the Thomson 
c     energy yield function for a specific energy.
c
c
      thom_ye = e / (e+ebd)**3 * 
     >       (1.0d0 - sqrt((e+ebd)/(gamma*(1.0d0-gamma)*eimp)))
c
c      write(7,'(a,5(1x,g12.4))') 'THOM YE:',e,thom_ye,eimp,
c     >                  ebd,gamma
c
      return
      end 
c
c
c
      double precision function thom_yv(v)
      implicit none
      double precision v
      common /thom_yv_params/  vimp,vgamma,vbd
      real*8 vimp,vgamma,vbd 
c
c     THOM_YV: 
c
c     This function returns the value of the Thomson 
c     energy yield function for a specific velocity.
c
c
      thom_yv = v**3 / (v**2+vbd**2)**3 * 
     >       (1.0d0 - sqrt((v**2+vbd**2)/
     >            (vgamma*(1.0d0-vgamma)*vimp**2)))
c
c      write(7,'(a,5(1x,g12.4))') 'THOM YV:',v,thom_yv,vimp,
c     >                  vbd,vgamma
c
      return
      end 
c
