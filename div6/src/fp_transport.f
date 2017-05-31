c     -*-Fortran-*-
c
      subroutine fp_transport(imp,ik,ir,iz,istate,s,theta,
     >                        cross,vel,temi,
     >                        particle_mass,nrand,
     >                        cist,cistfp,cstmax,ctemav,rsect,zsect,
     >                        sputy,rc)
      use error_handling
      use debug_options
      use mod_fp_transport
      implicit none
      integer ik,ir,iz,istate,rc,imp,nrand
      real s,cross,vel,temi,theta
      real particle_mass
      real cstmax,ctemav
      real*8 cist,cistfp
      real rsect,zsect
      real sputy
c
c
      include 'params'
      include 'cgeom'
      include 'driftvel'
      include 'fperiph_com'
c
c     Every outermost ring of the grid that does not have a wall as its boundary - forms a 
c     far periphery region where this transport code can be used. 
c
c     The plasma conditions in the periphery are defined by the plasma data stored in separate
c     arrays. 
c
c
c     This subroutine implements parallel particle transport in the 
c     peripheral plasma region.The particle is always associated with the 
c     corresponding cell on the outermost real ring of the grid. These rings
c     are IRWALL-1 and IRTRAP+1. There are several options for the forces acting on
c     the particle:
c          1) All forces calculated as if the particle was in the outermost ring
c          2) No forces - just perpendicular transport
c          3) Simple velocity diffusion. 
c      
c
c     fp_plasma already exists and is assigned values - need to see how to couple this to the 
c     transport routines. 
c
c
c
c     When the particle returns to the grid it returns to the edge of the first real
c     cell in rings IRWALL-1 and IRTRAP+1 in the cell corresponding to its current 
c     location. 
c
c     Wall impact - there are a number of options for calculating this - 
c        1) Calculate an R,Z position for the particle - check to see if this is outside
c           the wall - if it outside - calculate a wall intersection point for the 
c           exit location of the particle. Problems: calculating R,Z is based on the 
c           associated cell geometry - when a cell is projected to the wall in a simple
c           way - the coverage on the wall is not complete - however, this would do a good
c           job when dealing with shadowed regions and complex wall geometries where 
c           only parts of the wall are exposed to the plasma. 
c        2) Precalculate Wdist(S) - the distance to the wall as a function of S along the
c           outermost ring - need a good algorithm for Wd(S) and also need a way to associate
c           the Wd(S) value to a specific wall element. 
c        3) 
c  
c
c        Status codes set to match those in FPERIPH
c
c        Is_Followed   : particle_status=0   (Particle is still being followed) 
c        Return to grid: particle_status=1   (Cross < 0 - assuming Cross is set to 0 when entering FP) 
c        Maximum time  : particle_status=2   (cistfp+cist > cstmax)
c        Wall Impact   : particle_status=3   (Cross > Wall_dist(S))
c        Target impact : particle_status=4   (S<0, S>Smax)
c
c
      real    fp_smax
      integer fp_ir,fp_reg
      integer fp_ran_used
      real fp_cross
      integer :: rc3_count = 0
      integer :: rc4_count = 0
c
c      write(6,'(a,i6,10(1x,g12.5)') 'FP:START:',imp,cist
c

c      call pr_trace('FP_TRANSPORT:','START FOLLOWING FP PARTICLE')
c
c     Assign periphery region involved
c
c     jdemod - this is a bit of a kludge to attempt to determine which 
c              periphery the particle exited. IF the exit condition
c              gives ir=the local boundary ring within 2 then it will
c              work
c
c
      if (ir.le.irwall.and.ir.ge.irwall-2) then
         fp_reg = fp_main
      elseif (ir.ge.irtrap.and.ir.le.irtrap+2) then 
         fp_reg = fp_pfz
      elseif (ir.le.irwall2.and.ir.ge.irwall2-2) then
         fp_reg = fp_main2
      elseif (ir.ge.irtrap2.and.ir.le.irtrap2+2) then 
         fp_reg = fp_pfz2
      else
         fp_reg = fp_main
         call errmsg('FP_TRANSPORT: Particle entering the'//
     >           ' FP not associated with a boundary ring region:',ir)
      endif

c
c     Code is assuming that the IK value passed in is compatible with the FP ring 
c     selected in fp_ir
c
      fp_ir = fp_rings(fp_reg)
      fp_smax = ksmaxs(fp_ir)
c
c     jdemod - make use of cross information to determine initial cross position of particle since setting it to zero 
c              will make back diffusion to the grid more likely. 
c
c     Set fp_cross to zero - all particles start at the edge of the fp and have a 50/50 chance more or less
c     of immediately diffusing back into the plasma - the former cross value isn't relavent at this point though 
c     a refinement would be to calculate how far into the fp the particle starts based on its cross coordinate and 
c     the cell characteristics.
c
c      if (cross.ne.0.0) then 
c      write(0,'(a,5i6,10(1x,g12.5))') 'FP Entry:', ik,ir,irwall,iz,
c     >              fp_reg,fp_cross_tmp,
c     >              cross, s, theta
c      write(6,'(a,5i6,10(1x,g12.5))') 'FP Entry:', ik,ir,irwall,iz,
c     >              fp_reg,fp_cross_tmp,
c     >              cross, s, theta
c      endif

c
c     Set the fp_cross value to the absolute value of fp_cross_tmp ... for transport into the PFZ
c     FP the value of cross > 0.0 while for the main SOL outermost ring cross < 0.0. In the FP, 
c     cross is positive everywhere. A value < 0 indicates a cross field step back onto the grid. 
c     fp_cross_tmp is set in ion_crossfield_transport.f
c
      fp_cross = abs(fp_cross_tmp)
c
c     Set the average neutral temperature for use in some of the
c     collision options
c
      if (ctemav.gt.0.0) then
         call fp_set_ctemav(ctemav)
      endif
c
c     Initialize the periphery code for the current particle
c     
c     Added the start and end of the drift region (if any) to the call
c
      call fp_init_particle(s,fp_cross,vel,temi,ik,fp_ir,
     >                     iz,istate,sputy,fp_reg,particle_mass,fp_smax,
     >                     fp_flow_velocity(fp_reg),
     >                     sdrft_start(fp_ir),sdrft_end(fp_ir))
c     >                      fp_sdrft_start(fp_reg),fp_sdrft_end(fp_reg))

c
      call fp_follow_particle(rc,imp,cist,cistfp,cstmax,
     >                        rsect,zsect,fp_ran_used)

c
c     Add random numbers used to total
c
      nrand = nrand + fp_ran_used
c
c     Assign fp results back to particle
c
      call fp_get_particle_data(s,cross,vel,temi,ik,ir,iz,istate,sputy)
c
c
c     Deal with particle status and accounting - set return code
c      
c     0) Continue following ...
c     1) Return to grid- calling code handles the particle - keep in mind that S-position and ik index
c                        may have changed as a result of far periphery transport. 
c     2) Out of time 
c     3) Struck Wall   - code returns the approximate wall element where particle struck the wall
c     4) Struck target - this is treated the same as the far periphery loss time
c                      - only difference is appropriate accounting of section of target
c                        stuck - calling code needs to be modified. 
c                      - this code returns wall element struck
c
c     When change state is enabled: (needs to be implemented)
c     - deal with neutral transitions ... pass back information? - allow for FP reionization?
c
c
c
c     Calling code needs to deal with return code and data appropriately. 
c
c      if (rc.eq.1) then 
c         write(6,'(a,4i6,6g12.5)') 'FP RTP :',rc,
c     >               ik,ir,iz,istate,s,cross,fp_flow_velocity(fp_reg)
c      endif
c
      if (rc.eq.3) then 
         rc3_count = rc3_count + 1
         write(6,'(a,4i6,6g12.5)') 'FP WALL:',rc3_count,
     >               ik,ir,iz,istate,s,cross,fp_flow_velocity(fp_reg)
      elseif (rc.eq.4) then 
         rc4_count = rc4_count + 1
         write(6,'(a,4i6,6g12.5)') 'FP TARG:',rc4_count,
     >               ik,ir,iz,istate,s,cross,fp_flow_velocity(fp_reg)
      endif

c      if (rc.eq.1) then 
c      write(0,'(a,5i6,10(1x,g12.5))') 'FP Exit :', ik,ir,irwall,iz,
c     >              fp_reg,fp_cross_tmp,
c     >              cross, s, theta
c      write(6,'(a,5i6,10(1x,g12.5))') 'FP Exit :', ik,ir,irwall,iz,
c     >              fp_reg,fp_cross_tmp,
c     >              cross, s, theta
c      endif
c
c     Note - the calling code expects the returned particle to be associated with the virtual boundary ring. 
c     Thus the IR of the particle needs to be mapped to the relavant ring. 
c
      if (fp_reg.eq.fp_main) then 
         ir = irwall
      elseif (fp_reg.eq.fp_pfz) then
         ir = irtrap
      elseif (fp_reg.eq.fp_main2) then 
         ir = irwall2
      elseif (fp_reg.eq.fp_pfz2) then
         ir = irtrap2
      endif
c
      return 
      end


c
c
c


      real function fp_delta_s_dperpz(ik,ir,nrand)
      implicit none
      integer ik,ir,nrand
      include 'params'
      include 'cgeom'
      include 'dperpz'
                                !
                                !     This routine returns a deltaS displacement that would result
                                !     from a cross-field step occurring in the Z or P (paramagnetic direction). 
                                !
      real fact,ran1
      real*8 seed

      fp_delta_s_dperpz = 0.0

      if (dperpz_opt.ne.0) then
                                !
                                !        Need to calculate Btor/Bpol from KBFS which is Btot/Bpol
                                !        Btot = sqrt(Bpol**2 + Btor**2)
                                !        KBFS**2 = (Bpol**2 + Btor**2) / Bpol**2 = 1 + (Btor/Bpol)**2
                                !        Btor/Bpol = sqrt(kbfs**2 -1)
                                !
         fact = sqrt(kbfs(ik,ir)**2-1.0)

         nrand = nrand + 1
         call getran(ran1)

         fp_delta_s_dperpz = sign(base_dperpz_step * fact,ran1-0.5)

      endif

      return
      end

C
C
C
      INTEGER FUNCTION FPERIPH(CIST,cistfp,FPXMAX,PROBLT,
     >                         CSTMAX,NRAND,DPERP,SEED,XSTART)
      IMPLICIT NONE
      real*8 cist,cistfp 
      REAL FPXMAX,PROBLT,CSTMAX,DPERP,XSTART
      DOUBLE PRECISION SEED
      INTEGER NRAND
c
C     FPERIPH:
C
C     THIS FUNCTION FOLLOWS THE ION IN THE FAR PERIPHERY
C     UNTIL IT IS ELIMINATED THROUGH ONE OF 4 ROUTES.
C
C     1) DIFFUSES BACK INTO THE PLASMA
C     2) EXCEEDS ION TIME-LIMIT
C     3) HITS THE WALL AT EDGE OF PERIPHERY
C     4) IS ASSUMED TO HIT TARGET PLATE IN FP REGION
C
      REAL X,RAN
      FPERIPH = 0
      X = XSTART
      cistfp = 0.0 
 
100   CONTINUE
C
      NRAND = NRAND + 1
      CALL SURAND2 (SEED, 1, RAN)
      IF (X.LT.0.0) THEN
         FPERIPH = 1
         GOTO 200
      ELSEIF ((CIST+cistfp).GT.CSTMAX) THEN
         FPERIPH = 2
         GOTO 200
      ELSEIF (X.GT.FPXMAX) THEN
         FPERIPH = 3
         GOTO 200
      ELSEIF (RAN.LE.PROBLT) THEN
         FPERIPH = 4
         GOTO 200
      ENDIF
C
      NRAND = NRAND + 1
      CALL SURAND2 (SEED, 1, RAN)
C
      X = X + SIGN(DPERP,RAN-0.5)
      cistfp = cistfp + 1.0
C
      GOTO 100
C
C     EXIT THROUGH COMMON POINT.
C
 200  CONTINUE
      RETURN
      END












