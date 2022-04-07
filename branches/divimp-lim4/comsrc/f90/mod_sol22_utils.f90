module mod_sol22_utils

  use mod_sol22_sources

  implicit none





contains



      SUBROUTINE SOL22Headers
! ======================================================================

! subroutine: SOL22Headers


      !use mod_params
      !use mod_slcom
      use mod_solparams
      use mod_solcommon
!     INCLUDE 'params'
!     INCLUDE 'slcom'
!     INCLUDE 'solparams'
!     INCLUDE 'solcommon'
      IMPLICIT none
      IF (sol22_osm_mode.LE.1) RETURN
      IF (miter.EQ.1) THEN
        WRITE(75,*)
        WRITE(75,'(A,A3,A8,1X,2A7,A10,1X,A10,A5,2(1X,A10))')'`','in','s','Ti','Te','ne','Vb','M','Ga','P'
        IF (sol22_outmode.GE.3) THEN
          WRITE(71,*)
          WRITE(71,'(A,A3,4(1X,A10,10X))')'`','in','pais','paes','peis','srcf'
        ENDIF
        IF (forcet.EQ.0.OR.forcet.EQ.2.OR.forcet.EQ.3) THEN
          WRITE(72,*)
          WRITE(72,'(A,A3,1X,A7,3A10,10X,A9,2X,3A10)')'`','in','Ti','Pcf','Pcx','Pei','Pu/Pt','Conv','Cond','Total'
        ENDIF
        WRITE(73,*)
!        WRITE(73,*)
!        WRITE(73,'(A,A3,1X,A7,3A10,A9,2X,3A10)')
!     .    '`','in','Te','Pcf','PHi','Pei','Pu/Pt',
!     .    'Conv','Cond','Total'
        WRITE(73,'(A,A3,1X,A7,4A10,A9,2X,3A10)')'`','in','Te','Pcf','PHi','Pei','Prad','Pu/Pt','Conv','Cond','Total'
      ENDIF
      RETURN
99    STOP



    END SUBROUTINE SOL22Headers




    SUBROUTINE SOL22Output(loopstart,spts,npts,conde,condi,conve,convi,pcxv,peiv,phelpiv,pradv,te,ti,ne,vb,ga,act_press,pmloss,note)
! ======================================================================

! subroutine: SOL22Output

      !use mod_params
      !use mod_comtor
      !use mod_slcom
      use mod_solparams
      use mod_solcommon
!     INCLUDE 'params'
!     INCLUDE 'comtor'
!     INCLUDE 'slcom'
!     INCLUDE 'solparams'
!     INCLUDE 'solcommon'
      IMPLICIT none


      ! jdemod - this common block does not exist in any of the other SOL22 source code modules
      !          in addition - mxspts is no longer a constant due to the shift to dynamic allocation
      !          and arrays with variable size are not allowed in common blocks.
      !          These quantities are also not assigned a value in this routine so I have commented
      !          out dp4 and dp6 and will remove references to them in this routine.
      !
      !COMMON /OUTPUTJUNK/ dp4        ,dp6
      !REAL                dp4(MXSPTS),dp6(MXSPTS)
      !REAL     GetCs_sol22

      !REAL*8   cond,conv,paes,pais,pmomloss,gperpf
      !EXTERNAL cond,conv,paes,pais,pmomloss,gperpf

      INTEGER i,loopstart,npts
      !REAL*8  srcf,powi,powe,mach,te(MXSPTS),spts (MXSPTS),pmloss(MXSPTS),exp_press(MXSPTS),ne(MXSPTS),&
      REAL*8  powi,powe,mach,te(MXSPTS),spts (MXSPTS),pmloss(MXSPTS),exp_press(MXSPTS),ne(MXSPTS),&
           prad (MXSPTS),pcxv  (MXSPTS),act_press(MXSPTS),ga(MXSPTS),peiv (MXSPTS),pradv (MXSPTS),&
           phelpiv  (MXSPTS),ti(MXSPTS),condi(MXSPTS),conde (MXSPTS),vb(MXSPTS),convi(MXSPTS),conve (MXSPTS)
      CHARACTER*2 note(MXSPTS)
      ! jdemod - dumpai1 is printed below, declared in slcom but never assigned any value so I am just replacing
      ! it with a local variable assigned a value of zero - I am commenting them out in slcom as well
      real*8 :: dumpai1,dumpae1,dumpei1
      dumpai1 = 0.0
      dumpae1 = 0.0
      dumpei1 = 0.0
      
      !     Pcx > 0 => cooling
!     PHi > 0 => cooling
!     Pei > 0 => electron cooling and ion heating
      IF (sol22_osm_mode.LE.1) RETURN
      DO i = loopstart, npts
        powi = (pai - pais(spts(i))) - pcxv   (i) + peiv(i)
        powe = (pae - paes(spts(i))) - phelpiv(i) - peiv(i)
        mach = DBLE(GetCs_sol22(SNGL(te(i)),SNGL(ti(i))))
!     .      srcf(spts(i)),gperpf(spts(i)),
        WRITE(70,'(1X,I3,F8.3,1X,2F7.2,1P,E10.2,1X,E10.2,0P,F5.2,1P,2(1X,E10.2),1X,2E10.2,0P,F10.4,A)')i,spts(i),ti(i),&
             te(i),ne(i),vb(i),DABS(vb(i)/mach),ga(i),act_press(i)/ECONV,ionsrc(i),gperpf(spts(i)),pmloss(i),note(i)
        IF (forcet.EQ.0.OR.forcet.EQ.2.OR.forcet.EQ.3) THEN
           WRITE(72,'(1X,I3,1X,F7.2,1P,3E10.2,0P,10X,F9.3,2X,1P,3E10.2,0P,A)')i,ti(i),-(pais(spts(i))-pai),-pcxv (i), &
                peiv(i),powi/pai,convi(i),condi(i),convi(i)+condi(i),note(i)
!...Prad!
        ENDIF                            
!        WRITE(73,'(1X,I3,1X,F7.2,1P,3E10.2,0P,F9.3,
!     .             2X,1P,3E10.2,0P,1X,2F8.2,A)')
!     .      i,te(i),-(paes(spts(i))-pae),-phelpiv(i),-peiv(i),powe/pae,
!     .      conve(i),conde(i),conve(i)+conde(i),dp4(i),dp6(i),note(i)
!
        ! jdemod - replace dp4 and dp6 with zeroes in teh following output - since they are undefined
        !
        !WRITE(73,'(1X,I3,1X,F7.2,1P,4E10.2,0P,F9.3,'//' 2X,1P,3E10.2,0P,1X,2F8.2,A)')i,te(i),-(paes(spts(i))-pae),&
        !     -phelpiv(i),-peiv(i),-pradv(i),powe/pae,conve(i),conde(i),conve(i)+conde(i),dp4(i),dp6(i),note(i)
        WRITE(73,'(1X,I3,1X,F7.2,1P,4E10.2,0P,F9.3,'//' 2X,1P,3E10.2,0P,1X,2F8.2,A)')i,te(i),-(paes(spts(i))-pae),&
             -phelpiv(i),-peiv(i),-pradv(i),powe/pae,conve(i),conde(i),conve(i)+conde(i),0.0,0.0,note(i)
        IF (sol22_outmode.GE.3)WRITE(71,'(1X,I3,1P,4(1X,2E10.2),3X,E10.2,A)')i,pais(spts(i)),dumpai1,paes(spts(i)),&
             dumpae1,peiv(i),dumpei1,srcf(spts(i)),-1.0,pmomloss(spts(i),1,vb(i),te(i),ti(i)),note(i)
      ENDDO
      RETURN
99    STOP

! ======================================================================


! ======================================================================

    END SUBROUTINE SOL22Output


  REAL FUNCTION GetCs_sol22(te,ti)
    use mod_solparams
    use mod_solcommon
    !use mod_params
    !use mod_comtor
    IMPLICIT none
    !     INCLUDE 'params'
    !     INCLUDE 'comtor'
    REAL, INTENT(IN) :: te,ti

    !     Te,i in eV
    !     a    in amu
    !     result is m s-1

    GetCs_sol22 = 9.78817E+03 * SQRT(0.5 * (1.0 + zb) * (te + ti) / mb)

    !write(0,*) 'GETCS:',getcs,te,ti,mb,zb

    RETURN
  END FUNCTION GetCs_sol22

    
end module mod_sol22_utils
