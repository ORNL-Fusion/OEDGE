module mtc

use error_handling
use hc_kinetics

integer,private :: local_mtcopt
real, private :: local_kelighi,local_kelighg

logical,private :: debug_mtc=.false.

save 

contains


  subroutine load_mtc_options(mtcopt,kelighi,kelighg)

    implicit none

    integer :: mtcopt
    real :: kelighi, kelighg ! Model coefficients  (m3/s)

    local_mtcopt = mtcopt
    local_kelighi = kelighi
    local_kelighg = kelighg

    if (debug_mtc) &
         & write(6,'(a,i5,10g12.5)') 'MTC OPTIONS LOADED:',local_mtcopt,&
         & local_kelighi, local_kelighg


  end subroutine load_mtc_options


  subroutine execute_mtc(mtcdat,mtccnt,mtccist,cist,mtcinf,xvelf,yvelf,& 
       &  sputy,vin,temn,cneutvel,fsrate,nrand,crmi,ti,v)


    use hc_get  ! for access to gktibs to get local ion temperature - could be added to argument list
    implicit none

    integer,intent(in) :: mtcdat
    !integer :: mtccnt,nrand,cneutvel,ik,ir
    integer :: mtccnt,nrand,cneutvel
    real*8 :: mtccist,cist
    real :: mtcinf(7,3)  ! Note mtcinf and hc_mtcinf must be kept to the same dimensions as used here
    real :: xvelf,yvelf,vin
    real :: xvelftmp, yvelftmp
    real :: sputy,temn
    real :: fsrate
    real :: crmi
    real :: ti
    real,dimension(3),optional :: v

    ! Local variables 
    real :: tmpvelf,vintmp
    real :: tmpvpol
    real :: vperp(3)
    !integer :: in
    real,external :: getranf


    !           Momentum Transfer Collision event
    !          
    !           Record event statistics

    if (debug_mtc.and.present(v)) then 
          write(6,'(a,3i5,10(1x,g12.5))') 'MTC2:',local_mtcopt,xvelf,yvelf,vin,v(1),v(2),v(3)
    elseif (debug_mtc) then 
          write(6,'(a,3i5,10(1x,g12.5))') 'MTC2:',local_mtcopt,xvelf,yvelf,vin
    endif


    mtccnt = mtccnt + 1 

    mtcinf(1,mtcdat) = mtcinf(1,mtcdat) + sputy

    if (mtccnt.eq.1) then 
       mtcinf(2,mtcdat) = mtcinf(2,mtcdat) + sputy * (cist-mtccist)
       mtcinf(4,mtcdat) = mtcinf(4,mtcdat) + sputy * (cist-mtccist)*abs(vin)*fsrate
    endif

    mtcinf(3,mtcdat) = mtcinf(3,mtcdat) + sputy * (cist-mtccist)
    mtcinf(5,mtcdat) = mtcinf(5,mtcdat) + sputy * (cist-mtccist)*abs(vin)*fsrate
    mtcinf(6,mtcdat) = mtcinf(6,mtcdat) + sputy * temn 
    mtcinf(7,mtcdat) = mtcinf(7,mtcdat) + sputy * abs(vin)

    !           Update MTC time 

    mtccist = cist

    ! Option 1 : Simple 90 degree turn in the 2D plane

    if (local_mtcopt.eq.1.or.(.not.present(v))) then 


       !           Calculate event (MTC)
       !          
       !           Choose randomly between +90 and -90 change in velocity
       !           - even probability.

       nrand = nrand + 1

       ! The random number counter kk is ONLY a counter in this routine - it can not be used as an index to a random number array
       ! since there is no supporting code to replenish the random number array or to check for out of bounds issues.
       !
       !if (granv(kk).ge.0.5) then 
       !
       if (getranf().ge.0.5) then 

          !              counter-clockwise 90 degrees

          tmpvelf = xvelf
          xvelf   = -yvelf
          yvelf   = tmpvelf

          if (present(v)) then
             tmpvelf = v(1)  ! tmp=Vr
             v(1) = -v(2)    ! Vr=-Vz
             v(2) = tmpvelf  ! Vz=tmp
             ! v(3) is not changed -> Vt
          endif

       else  

          !              clockwise 90 degrees

          tmpvelf = xvelf
          xvelf   = yvelf
          yvelf   = -tmpvelf

          if (present(v)) then
             tmpvelf = v(1)  ! tmp=Vr
             v(1) = v(2)     ! Vr=Vz
             v(2) = -tmpvelf ! Vz=-tmp
             ! v(3) is not changed -> Vt
          endif

       endif


    elseif (local_mtcopt.eq.2) then 

       ! Note VIN should be consistent with v in that |v(r,z)| = VIN.
       ! VIN is the total velocity in the R,Z plane. 
       ! This code has the side effect of changing VIN since the component in the 
       ! R,Z plane will change when a 90 degree change of path is executed in 3D. 

       ! This option requires a 3-space vector velocity v be specified for
       ! the particle - it then will perform a 90 degree collision along 
       ! a random vector in the plane perpendicular to v - this will then
       ! be mapped to the R,Z motion of the particle by assigning the R and 
       ! Z components of the vector to xvelf and yvelf with appropriate
       ! scaling. 

       ! Returns a random unit vector in the plane perpendicular to v

       call get_random_perp_vector(v,vperp)

       !
       ! magnitude of V is NOT the same as VIN since VIN is only the component of the
       ! velocity in the R,Z plane - really - VIN should equal sqrt(v(1)**2+v(2)**2) 
       ! and that is what should be checked.
       !

       tmpvelf =  sqrt(sum(v**2)) ! total magnitude of V

       tmpvpol =  sqrt(v(1)**2+v(2)**2) ! magnitude of V in the poloidal plane


       if (debug_mtc) &
            & write(6,'(a,15f18.10)') 'MTC:VPERP1:',tmpvelf,vin,tmpvpol,abs(tmpvpol-vin),&
            & v(1),v(2),v(3),vperp(1),vperp(2),vperp(3),&
            &                        dot_product(v,vperp)


       !
       ! Issue an error message if velocities differ by greater than 1m/s
       !
       if (abs(tmpvpol-vin)  .gt. 1.0) then 
          write(6,'(a,10f18.8)') 'ERROR: MTC VPERP RESULT:',tmpvpol,tmpvelf,vin,abs(tmpvpol-vin)
          call errmsg('EXECUTE_MTC:','POSSIBLE ERROR: VELOCITY VALUES DIFFER BY:',abs(tmpvpol-vin))
       endif

       !
       ! use total 3D velocity to determine new components
       !
       xvelf = vperp(1) * fsrate * tmpvelf ! R component
       yvelf = vperp(2) * fsrate * tmpvelf ! Z component
       !xvelf = vperp(1) * fsrate * vin ! R component
       !yvelf = vperp(2) * fsrate * vin ! Z component

       v(1) = tmpvelf * vperp(1)  ! new Vr
       v(2) = tmpvelf * vperp(2)  ! new Vz
       v(3) = tmpvelf * vperp(3)  ! new Vt

       ! Assign new value to vin
       vin =  sqrt(v(1)**2+v(2)**2) 

       if (debug_mtc) &
            & write(6,'(a,15f18.10)') 'MTC:VPERP2:',tmpvelf,vin,tmpvpol,abs(tmpvpol-vin),&
            & v(1),v(2),v(3),vperp(1),vperp(2),vperp(3),&
            &                        dot_product(v,vperp)


    endif

    !           If Neutral Velocity Option 2 is in effect - change
    !           the neutral velocity in magnitude based on local 
    !           conditions.

    if (cneutvel.eq.2) then 

       !              Record old values

       xvelftmp = xvelf
       yvelftmp = yvelf
       vintmp   = vin

       !              Assign new speed

       !vin = 1.38e4 * sqrt(gktibs(ik,ir)/crmi)
       vin = 1.38e4 * sqrt(ti/crmi)

       !              Assign new temperature

       TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)

       !              Reset velocity components

       xvelf = xvelftmp * vin/vintmp              
       yvelf = yvelftmp * vin/vintmp              

       !
       ! Also adjust v if required
       ! 
       if (present(v)) then 
          v(1) = v(1) * vin/vintmp
          v(2) = v(2) * vin/vintmp
          v(3) = v(3) * vin/vintmp
       endif

    endif


    if (debug_mtc.and.present(v)) then 
          write(6,'(a,3i5,10(1x,g12.5))') 'MTC2:',local_mtcopt,xvelf,yvelf,vin,v(1),v(2),v(3)
    elseif (debug_mtc) then 
          write(6,'(a,3i5,10(1x,g12.5))') 'MTC2:',local_mtcopt,xvelf,yvelf,vin
    endif

  end subroutine execute_mtc



  real function mtc_prob(background_mass,impurity_mass,ik,ir,&
                          &neutral_timestep)
    
    use hc_get
    implicit none
    
    !real :: mtc_prob   ! Probability of significant momentum transfer collision - i.e. 90 degrees
    integer :: ik,ir   ! cell indices
    real :: background_mass,impurity_mass  ! Impurity mass and background plasma mass
    !real :: kelighi, kelighg ! Model coefficients  (m3/s)
    real :: neutral_timestep     ! neutral time step


    ! NOTE: Using a linear form is based on the assumption that the momentum transfer collision 
    !       probability is small. 
    !       This is based on the assumption that the probability density for a momentum transfer 
    !       collision has an exponential form:
    !       f(t) = lam exp (-lam * t) 
    !       P(dt) = integral (0 to dt)  f(t) =   1.0 - exp (-lam dt) 
    !       lam = event rate in s-1
    !       lam = 16 * crmb / (3.0 * (crmi + crmb)) * (kelighi * ni + kelighg * nh) 
    !
    !       If we expand in taylor series and assume lam * dt is small - this reduces to the 
    !       expression used here. 
    !
    !       This is the same probability theory used in the DIVIMP change of state code - I am 
    !       documenting it here for reference. 
    !
    !       Also note!: For small probabilities these are additive
    !
    !       P(a+b) = 1.0 - exp (-(a+b) dt) 
    !       Expanding in taylor series assuming a+b is small yields 
    !       P(a+b) = a dt + b dt = P(a) + P(b) - thus the code uses this form when looking at combined
    !       change of state probabilities. 
    !
    !       kelighi and kelighg are model coefficients in units of m3/s
    !       knbs is plasma ion density in m-3
    !       nh is the neutral background density in m-3

    if (local_mtcopt.eq.0) then 
       mtc_prob = 0.0
    else

       mtc_prob = 16.0 * background_mass / (3.0 * (background_mass + impurity_mass)) *  &
              (local_kelighi * gknbs(ik,ir) + local_kelighg * gnh(ik,ir)) * neutral_timestep
    endif
       
  end function mtc_prob



end module mtc
