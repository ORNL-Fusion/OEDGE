module mod_assignpp


  implicit none







contains

  SUBROUTINE ThomPP(irlim1,irlim2,ikopt,targ_con_opt)
  ! ======================================================================

  ! subroutine: ThomPP

    use mod_params
    use mod_cgeom
    use mod_pindata
    use mod_slcom

    !     irlim1 - start ring
    !     irlim2 - end ring
    !     ikopt  - apply to whole or partial ring
    !     targ_con_opt   =3 - do not apply conditions to target
    !                    =4 - apply average values to target


    !     INCLUDE 'params'
    !     INCLUDE 'cgeom'
    !     INCLUDE 'pindata'
    !     INCLUDE 'slcom'

    !      INCLUDE 'solparams'
    !      INCLUDE 'solswitch'
    !      INCLUDE 'solcommon'

    IMPLICIT none
    INTEGER irlim1,irlim2,ikopt,targ_con_opt
    integer irstart,irend,ikstart,ikend
    REAL GetCs
    INTEGER    MAXTDAT     ,MAXCOLS
    PARAMETER (MAXTDAT=1000,MAXCOLS=10)
    ! slmod begin - new
    INTEGER thnnraw

    !      REAL    dummy(MAXNKS,MAXNRS),thnraw(MAXTDAT,MAXCOLS)
    ! slmod end
    REAL    thnraw(MAXTDAT,MAXCOLS)
    INTEGER ik,ir,ir1,ir2,in,tag(MAXNRS),num
    LOGICAL status,output,loaded
    ! slmod begin - new
    REAL    te(2,MAXNRS),ne(2,MAXNRS),d1,d2,f
    DATA loaded /.FALSE./
    ! slmod end
    SAVE
    output = .TRUE.
    ! slmod begin - new
    status = .FALSE.
    ! slmod end
    loaded = .false.
    CALL IZero(tag,MAXNRS)
    CALL RZero(te ,MAXNRS*2)

    CALL RZero(ne ,MAXNRS*2)
    ! slmod begin - new
    if (.not.loaded) then

       !         CALL LoadThomsonData(dummy,thnnraw,thnraw,MAXTDAT,MAXCOLS,2,2)
       ! slmod end
       CALL LoadThomsonData(thnnraw,thnraw,MAXTDAT,MAXCOLS,NIL,NIL,2)
       loaded = .true.
       IF (output) THEN
          DO in = 1, thnnraw
             WRITE(PINOUT,'(A,I6,4F10.4,1P,E10.2,0P,F10.4)')'THOMSON : ', in,(thnraw(in,ir),ir=1,6)
          ENDDO

          !...     Find average temperature and density for rings with Thomson data:

       ENDIF
       DO ir = 1, nrs
          num = 0
          DO in = 1, thnnraw
             IF (INT(thnraw(in,4)).EQ.ir) THEN
                num = num + 1
                f = 1.0 / REAL(num)
                te(IKLO,ir) =  f * thnraw(in,6)+ (1.0 - f) * te(IKLO,ir)
                ne(IKLO,ir) =  f * thnraw(in,5)+ (1.0 - f) * ne(IKLO,ir)
                te(IKHI,ir) = te(IKLO,ir)
                ne(IKHI,ir) = ne(IKLO,ir)
                if (ir.ge.irtrap.and.ir.le.nrs) status = .TRUE.
                IF (output)WRITE(PINOUT,'(A,2I6,F6.3,2F10.4,1P,2E10.2,0P)')' A : ',ir,num,f,thnraw(in,6),te(IKLO,ir),&
                     thnraw(in,5),ne(IKLO,ir)
             ENDIF
          ENDDO
          IF (num.GT.0) tag(ir) = 1
          IF (output)WRITE(PINOUT,'(A,3I6,F10.4,1P,E10.2,0P)') ' B : ',ir,tag(ir),num,te(IKLO,ir),ne(IKLO,ir)

          !...     Make sure there is some Thomson data in the PP:
       ENDDO
       !...     Assign values to rings without Thomson data:
       IF (ir.ge.irtrap.and.ir.le.nrs.and.(.NOT.status)) CALL ER('ThomPP','No data in PP',*99)
       DO ir = 1, nrs
          IF (ir.ge.irtrap.and.ir.le.nrs) then
             irstart = irtrap
             irend = nrs
          elseif (ir.ge.1.and.ir.le.irwall) then
             irstart = 1
             irend = irwall

          endif
          IF (tag(ir).EQ.0) THEN
             ! slmod begin - new
             !...fix
             ir1 = ir
             DO WHILE (ir1.GT.irstart.AND.tag(ir1).EQ.0)
                ir1 = ir1 - 1
             ENDDO
             ir2 = ir
             DO WHILE (ir2.LT.irend+1.AND.tag(ir2).EQ.0)
                ir2 = ir2 + 1

                !             DO WHILE (ir.GT.irstart.AND.tag(ir1).EQ.0)
                !               ir1 = ir1 - 1
                !             ENDDO
                !             ir2 = ir
                !             DO WHILE (ir.LT.irend+1.AND.tag(ir2).EQ.0)
                !               ir2 = ir2 + 1
                !             ENDDO
                ! slmod end

             ENDDO
             if (ir1.eq.irstart.and.ir2.eq.irend+1) then
                te(IKLO,ir) = 1.0
                ne(IKLO,ir) = 1.0e10
                te(IKHI,ir) = te(IKLO,ir)
                ne(IKHI,ir) = ne(IKLO,ir)
             elseif (ir1.EQ.irstart) THEN
                te(IKLO,ir) = te(IKLO,ir2)
                ne(IKLO,ir) = ne(IKLO,ir2)
                te(IKHI,ir) = te(IKLO,ir)
                ne(IKHI,ir) = ne(IKLO,ir)
             ELSEIF (ir2.EQ.irend+1) THEN
                te(IKLO,ir) = te(IKLO,ir1)
                ne(IKLO,ir) = ne(IKLO,ir1)
                te(IKHI,ir) = te(IKLO,ir)
                ne(IKHI,ir) = ne(IKLO,ir)
             ELSE
                d1 = REAL(ir  - ir1)
                d2 = REAL(ir2 - ir )
                te(IKLO,ir) = (d2*te(IKLO,ir1)+ d1*te(IKLO,ir2)) / (d1+d2)
                ne(IKLO,ir) = (d2*ne(IKLO,ir1)+ d1*ne(IKLO,ir2)) / (d1+d2)
                te(IKHI,ir) = te(IKLO,ir)
                ne(IKHI,ir) = ne(IKLO,ir)

             ENDIF
             IF (output)WRITE(PINOUT,'(A,3I6,F10.4,1P,E10.2,0P)') ' C : ',ir,ir1,ir2,te(IKLO,ir),ne(IKLO,ir)

          ENDIF

       ENDDO


       !...  Assign plasma to rings with Thomson data:

    endif

    !        IF (ir.LT.irtrap.OR.ir.GT.nrs) CYCLE

    DO ir = irlim1, irlim2

       call set_ikvals(ir,ikstart,ikend,ikopt)
       DO ik = ikstart, ikend
          ktebs(ik,ir) = te(IKLO,ir)
          ktibs(ik,ir) = te(IKLO,ir)
          knbs (ik,ir) = ne(IKLO,ir)
          kvhs (ik,ir) = 0.0
          !...    Assign target values using Thomson data:
       ENDDO

       IF (targ_con_opt.eq.4) THEN
          if (ikopt.eq.2.or.ikopt.eq.3) then
             kteds(idds(ir,1)) = ktebs(nks(ir),ir)
             ktids(idds(ir,1)) = ktibs(nks(ir),ir)
             knds (idds(ir,1)) = knbs (nks(ir),ir)
             kvds (idds(ir,1)) =GetCs(kteds(idds(ir,1)),ktids(idds(ir,1)))
          endif
          if (ikopt.eq.1.or.ikopt.eq.3) then
             kteds(idds(ir,2)) = ktebs(1      ,ir)
             ktids(idds(ir,2)) = ktibs(1      ,ir)
             knds (idds(ir,2)) = knbs (1      ,ir)
             kvds (idds(ir,2)) =-GetCs(kteds(idds(ir,2)),ktids(idds(ir,2)))
          endif

       ENDIF
    ENDDO
    RETURN
99  STOP



  END SUBROUTINE ThomPP

  SUBROUTINE AssignPP(irlim1,irlim2,ikopt,targ_con_opt)
    ! subroutine: AssignPP

    ! Assign a uniform plasma to the private flux zone from data listed in
    ! DIVIMP input file.
    use mod_params
    use mod_cgeom
    use mod_comtor
    use mod_pindata
    use mod_slcom

    !     irlim1 - start ring
    !     irlim2 - end ring
    !     ikopt  - apply to whole or partial ring
    !     targ_con_opt   =5 - do not apply plasma values to target
    !                    =6 - apply plasma values to target

    !     INCLUDE 'params'
    !     INCLUDE 'cgeom'
    !     INCLUDE 'comtor'
    !     INCLUDE 'pindata'
    !     INCLUDE 'slcom'
    IMPLICIT none
    INTEGER SymmetryPoint
    REAL    GetCs
    INTEGER irlim1,irlim2,ikopt,targ_con_opt
    INTEGER ik,ikm,ir,i1
    LOGICAL output
    REAL    te(2,MAXNRS),ti(2,MAXNRS),ne(2,MAXNRS),vb(2,MAXNRS)
    output = .TRUE.
    IF (osmnppv.EQ.0) CALL ER('AssignPP','No plasma data from input'//'file found',*99)
    CALL RZero(te,MAXNRS*2)
    CALL RZero(ti,MAXNRS*2)
    CALL RZero(ne,MAXNRS*2)
    !...  Search input file listing for ring plasma data:
    CALL RZero(vb,MAXNRS*2)
    DO ir = irlim1, irlim2
       DO i1 = 1, osmnppv
          IF (INT(osmppv(i1,1)).EQ.ir) THEN
             IF (osmppv(i1,2).EQ.1.0.OR.osmppv(i1,2).EQ.3.0) THEN
                te(IKLO,ir) =  osmppv(i1,3)
                ti(IKLO,ir) =  osmppv(i1,4)
                ne(IKLO,ir) =  osmppv(i1,5)
                vb(IKLO,ir) = -osmppv(i1,6)
             ENDIF
             IF (osmppv(i1,2).EQ.2.0.OR.osmppv(i1,2).EQ.3.0) THEN
                te(IKHI,ir) =  osmppv(i1,3)
                ti(IKHI,ir) =  osmppv(i1,4)
                ne(IKHI,ir) =  osmppv(i1,5)
                vb(IKHI,ir) =  osmppv(i1,6)
             ENDIF
          ENDIF
       ENDDO
       !...  Assign plasma to rings with Thomson data:
    ENDDO
    DO ir = irlim1, irlim2
       ikm = SymmetryPoint(ir)
       IF (ikopt.EQ.1.OR.ikopt.EQ.3) THEN
          IF (te(IKLO,ir).EQ.0.0)CALL ER('AssignPP','Required low index plasma data not '//'found',*99)
          DO ik = 1, ikm
             ktebs(ik,ir) = te(IKLO,ir)
             ktibs(ik,ir) = te(IKLO,ir)
             knbs (ik,ir) = ne(IKLO,ir)
             kvhs (ik,ir) = vb(IKLO,ir)
          ENDDO
       ENDIF
       IF (ikopt.EQ.2.OR.ikopt.EQ.3) THEN
          IF (te(IKHI,ir).EQ.0.0)CALL ER('AssignPP','Required high index plasma data not '//'found',*99)
          DO ik = ikm+1, nks(ir)
             ktebs(ik,ir) = te(IKHI,ir)
             ktibs(ik,ir) = te(IKHI,ir)
             knbs (ik,ir) = ne(IKHI,ir)
             kvhs (ik,ir) = vb(IKHI,ir)
          ENDDO
          !...    Assign target values using Thomson data:
       ENDIF
       IF (targ_con_opt.EQ.6) THEN
          IF (ikopt.EQ.2.OR.ikopt.EQ.3) THEN
             kteds(idds(ir,1)) = ktebs(nks(ir),ir) * te_mult_i
             ktids(idds(ir,1)) = ktibs(nks(ir),ir) * ti_mult_i
             !...          Sometimes Isat or Ne data in the input file uses a multiplier
             !             so that the numbers listed in the target data arrays do not
             !             have to include the exponents (e.g. 1.0E+20 can be entered as
             !             1.0 if a 1.0E+20 multiplier is being used).  However,
             !             such vile conduct will break this code, which is assigning
             !             the target data from the bulk plasma prescription but also applies
             !             the target multiplier:
             IF (n_mult_i.GT.10) THEN
                CALL ER('AssignPP','Suspicious high index target data '//'multiplier detected',*99)
             ELSE
                knds(idds(ir,1)) = knbs(nks(ir),ir) * n_mult_i
             ENDIF
             kvds (idds(ir,1)) =GetCs(kteds(idds(ir,1)),ktids(idds(ir,1)))
          ENDIF
          IF (ikopt.EQ.1.OR.ikopt.EQ.3) THEN
             kteds(idds(ir,2)) = ktebs(1,ir) * te_mult_o
             ktids(idds(ir,2)) = ktibs(1,ir) * ti_mult_o
             !...          See the above comment for the high index target:
             IF (n_mult_o.GT.10) THEN
                CALL ER('AssignPP','Suspicious low index target data '//'multiplier detected',*99)
             ELSE
                knds(idds(ir,2)) = knbs(1,ir) * n_mult_o
             ENDIF
             kvds (idds(ir,2)) =-GetCs(kteds(idds(ir,2)),ktids(idds(ir,2)))
          ENDIF
       ENDIF
    ENDDO
    RETURN
99  STOP



  END SUBROUTINE AssignPP









end module mod_assignpp
