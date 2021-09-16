module mod_sol22_divimp



  implicit none




contains





  subroutine print_edge2d_flux_analysis(rconst,gtarg,areasum,ir,gperpa,oldknbs,grad_oldknbs,srcinteg,recinteg,ike2d_start)
    use mod_params
    use mod_solparams
    use mod_solswitch
    use mod_comtor
    use mod_cgeom
    use mod_pindata
    use mod_cedge2d
    !     include 'params'
    !     include 'solparams'
    !     include 'solswitch'
    implicit none
    integer ir,ike2d_start
    real    gperpa(maxnks,maxnrs),oldknbs(maxnks,maxnrs)
    real    grad_oldknbs(maxnks,maxnrs)
    real*8 rconst(mxspts),areasum(mxspts)
    real*8 gtarg(mxspts,3)
    real srcinteg (maxnks+3,maxnrs)

    !     include 'comtor'
    !     include 'cgeom'
    !     include 'pindata'
    !     include 'cedge2d'

    !     PRINT_EDGE2D_FLUX_ANALYSIS:


    !     This routine prints a number of different analyses of the
    !     EDGE2D fluxes in the context of SOL option 22.



    !     Local variables

    real recinteg (maxnks+3,maxnrs)
    integer ik,in,id,startik,endik

    !     New variables for summaries

    real rmean1,rmean2
    real intgrad,intgradpos,intgradneg,inte2dsrc,intdist
    real inte2dhrec,intgrad2,ds,dcell,mulfact
    integer ikin,irin,ikout,irout
    real    influx,outflux,netflux,intnflux,flux1,flux2,flux3
    real    flux4
    real fluxst,fluxend,ioniz,rec,gperpd,gperpn,gdiv,sider
    real intgperpd,intgnet,intgdiv
    real initflux,dp,dp1,dp2,brat1,brat2,startflux

    real flux_const,endflux
    if (switch(swe2d).eq.0.0) then
       startik = 1
       endik   = nks(ir)
    elseif (switch(swe2d).gt.0.0) then
       startik = ike2d_start
       endik   = nks(ir) - ike2d_start + 1

    endif
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    !     ------------------------------------------------------------------

    intnflux = 0.0

    write(6,'(20a)') 'FLUX:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ',&
         '  |GPERP   ','  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ',&
         '    GPERPD  ','    E2DION  ','    RECSRC   ','   NETFLUX   ','   SUM2'

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)


       if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          endif

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
       flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1

       flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,15(1x,g11.4))') 'F:',ir,ik,kss(ik,ir),ds,gtarg(ir,2),inte2dsrc,&
            0.15*intgrad2,-inte2dhrec,intnflux,flux1,flux2,flux1/flux2,0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,&
            e2dhrec(ik,ir)*ds,netflux*ds,flux3,flux3/flux2
       if (ik.eq.endik) then
          if (intgrad2.gt.0.0) then 
             mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
          else
             mulfact = 0.0
          endif

       endif


       if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          endif
       endif

    end do

    ik = nks(ir)
    flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
    flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)

    flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec

    !     ------------------------------------------------------------------

    write (6,'(a,5g12.5)') 'NUM:',gtarg(ir,2),0.15*intgrad2,flux1,flux2,mulfact

    write(6,'(a,2i4,1p,8(1x,g13.5))') 'F:',ir,ik,gtarg(ir,2),inte2dsrc,0.15*intgrad2,-inte2dhrec,&
         gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec,e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),&
         (gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec)/(e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),gtarg(ir,1)

    write (6,*)
    write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact

    write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact

    write(6,'(20a)') 'FLUX:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ',&
         '  |GPERP   ','  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ',&
         '    GPERPD  ','    E2DION  ','    RECSRC   ','   NETFLUX   ','   SUM2'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)


       if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1 * mulfact

          intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          endif

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif

       write(6,'(a,2i4,2f9.4,1p,15(1x,g11.4))') 'FLUX:',ir,ik,kss(ik,ir),ds,gtarg(ir,2),inte2dsrc,0.15*intgrad2,&
            -inte2dhrec,intnflux,gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec,e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1,&
            (gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec)/(e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1),&
            0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,e2dhrec(ik,ir)*ds,netflux*ds,&
            gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec


       if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2 * mulfact

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          endif
       endif

    end do

    !     ------------------------------------------------------------------

    !           Analysis of the actual EDGE2D boundary fluxes


    write (6,*)
    write (6,500)

500 format('BND:',2x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',&
         2x,1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,3x,&
         'D-RAT',4x,2x,'E2DV-RAT')
    intgperpd = 0.0
    intgnet   = 0.0
    intgdiv   = 0.0

    do ik = startik,endik
       fluxst  = e2dflux(ik,ir)

       fluxend = e2dflux(ik+1,ir)

       ioniz = e2dion(ik,ir)

       rec   = e2dhrec(ik,ir)

       ds    = ksb(ik,ir) - ksb(ik-1,ir)

       gperpn= -(fluxst - fluxend+ e2dion(ik,ir)*ds - e2dhrec(ik,ir)*ds)

       gperpd= gperpn/ds
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)

       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       gdiv =  0.15*grad_oldknbs(ik,ir)
       intgperpd = intgperpd + ds * gperpd
       intgnet   = intgnet   + ds * netflux

       intgdiv   = intgdiv   + ds * gdiv

       write (6,'(a,2i4,1p,15(1x,g11.4))') 'BND:',ik,ir,fluxst,ioniz*ds,-rec*ds,gperpn,fluxend,gperpd,&
            netflux,gdiv,intgperpd,intgnet,intgdiv,gdiv/gperpd,netflux/gperpd

       !     ------------------------------------------------------------------

    end do
    write (6,*)
    write (6,*) 'DOWN FLUX ANALYSIS:  1/2 CELL'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0

    write(6,'(20a)') 'DFH:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ','  |GPERP   ',&
         '  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ','    E2DION  ',&
         '    RECSRC   ','   NETFLUX   ','   SUM2'
    startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)


       if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          endif

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'DFH:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,0.15*intgrad2,&
            -inte2dhrec,intnflux,flux1,flux2,flux1/flux2,0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,&
            e2dhrec(ik,ir)*ds,netflux*ds,flux3,flux3/flux2
       if (ik.eq.endik) then
          mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0

       endif


       if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          endif
       endif

    end do

    ik = nks(ir)
    flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
    flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
    flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec

    flux4 =  e2dgpara(ik+1,ir)

    !     ------------------------------------------------------------------

    write (6,'(a,5g12.5)') 'NUM:',startflux,0.15*intgrad2,flux1,flux2,mulfact

    write(6,'(a,2i4,1p,8(1x,g13.5))') 'F:',ir,ik,startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,&
         startflux+inte2dsrc+0.15*intgrad2-inte2dhrec,e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),&
         (startflux+inte2dsrc+0.15*intgrad2-inte2dhrec)/(e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),flux4

    write (6,*)
    write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact

    write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact

    write(6,'(20a)') 'MDFH:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ',&
         '  |GPERP   ','  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',&
         '    E2DION  ','    RECSRC   ','   NETFLUX   ','   SUM2'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0

    startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)


       if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1 * mulfact

          intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          endif

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MDFH:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,&
            0.15*intgrad2,-inte2dhrec,intnflux,flux1,flux2,flux1/flux2,0.15 * grad_oldknbs(ik,ir)*ds,&
            e2dion(ik,ir)*ds,e2dhrec(ik,ir)*ds,netflux*ds,flux3,flux3/flux2


       if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2 * mulfact

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          endif
       endif

    end do

    !     ------------------------------------------------------------------
    !     ------------------------------------------------------------------

    write (6,*)
    write (6,*)
    write (6,*) 'DOWN FLUX ANALYSIS:  FULL CELL'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0

    write(6,'(20a)') 'DFF:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ',&
         '  |GPERP   ','  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ',&
         '    GPERPD  ','    E2DION  ','    RECSRC   ','   NETFLUX   ','   SUM2'
    startflux = e2dgpara(1,ir)

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)

       intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'DFF:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,0.15*intgrad2,&
            -inte2dhrec,intnflux,flux1,flux2,flux1/flux2,0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,e2dhrec(ik,ir)*ds,&
            netflux*ds,flux3,flux3/flux2
       if (ik.eq.endik) then
          mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0

       endif

       intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       endif
       if (ik.eq.endik) then
          mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0

       endif

    end do

    ik = nks(ir)
    flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
    flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
    flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec

    flux4 = e2dgpara(ik+1,ir)

    !     ------------------------------------------------------------------


    write (6,'(a,5g12.5)') 'NUM:',startflux,0.15*intgrad2,flux1,flux2,mulfact

    write(6,'(a,2i4,1p,8(1x,g13.5))') 'DFF:',ir,ik,startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,&
         startflux+inte2dsrc+0.15*intgrad2-inte2dhrec,e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),&
         (startflux+inte2dsrc+0.15*intgrad2-inte2dhrec)/(e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),flux4

    write (6,*)
    write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact

    write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact

    write(6,'(20a)') 'MDFF:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ','  |GPERP   ',&
         '  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ','    E2DION  ',&
         '    RECSRC   ','   NETFLUX   ','   SUM2'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0

    startflux = e2dgpara(1,ir)

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)

       intgrad = intgrad+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1 * mulfact

       intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0


       endif
       flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MDFF:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,&
            0.15*intgrad2,-inte2dhrec,intnflux,flux1,flux2,flux1/flux2,0.15 * grad_oldknbs(ik,ir)*ds,&
            e2dion(ik,ir)*ds,e2dhrec(ik,ir)*ds,netflux*ds,flux3,flux3/flux2

       intgrad = intgrad+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2 * mulfact

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       endif

    end do

    !     ------------------------------------------------------------------
    !     ------------------------------------------------------------------

    write (6,*)
    write (6,*)
    write (6,*) 'FLUX ANALYSIS:  1/2 CELL - REGULAR DIST.'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0

    write(6,'(20a)') 'GFH:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ','  |GPERP   ',&
         '  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ','    E2DION  ',&
         '    RECSRC   ','   NETFLUX   ','   SUM2'
    startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0

    endflux   = (e2dgpara(nks(ir),ir)+ e2dgpara(nks(ir)+1,ir))/2.0

    flux_const= (startflux - endflux + srcinteg(nks(ir)+1,ir)- recinteg(nks(ir)+1,ir))/areasum(ir)

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)


       if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ flux_const*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          endif

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'GFH:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,intgrad,&
            -inte2dhrec,intnflux,flux1,flux2,flux1/flux2,flux_const*ds,e2dion(ik,ir)*ds,e2dhrec(ik,ir)*ds,&
            netflux*ds,flux3,flux3/flux2
       if (ik.eq.endik) then
          mulfact = ((flux2-flux1)/intgrad) + 1.0

       endif


       if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ flux_const*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          endif
       endif

    end do

    ik = nks(ir)
    flux1 =  gtarg(ir,2)+inte2dsrc+intgrad-inte2dhrec
    flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
    flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec

    flux4 =  e2dgpara(ik+1,ir)

    write (6,'(a,5g12.5)') 'NUM:',startflux,0.15*intgrad2,flux1,flux2,mulfact

    !     ------------------------------------------------------------------


    write (6,*)
    write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact

    write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact

    write(6,'(20a)') 'MGFH:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ','  |GPERP   ',&
         '  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ','    E2DION  '&
         ,'    RECSRC   ','   NETFLUX   ','   SUM2'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0
    startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0

    endflux   = (e2dgpara(nks(ir),ir)+ e2dgpara(nks(ir)+1,ir))/2.0

    flux_const= (startflux -endflux + srcinteg(nks(ir)+1,ir)- recinteg(nks(ir)+1,ir))/areasum(ir)

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)


       if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ flux_const*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1 * mulfact

          intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

          endif

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MGFH:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,intgrad,&
            -inte2dhrec,intnflux,flux1,flux2,flux1/flux2,flux_const*ds,e2dion(ik,ir)*ds,e2dhrec(ik,ir)*ds,&
            netflux*ds,flux3,flux3/flux2


       if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.&
            switch(swe2d).eq.9))) then

          intgrad = intgrad+ flux_const*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2 * mulfact

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

          endif
       endif

    end do

    !     ------------------------------------------------------------------
    !     ------------------------------------------------------------------

    write (6,*)
    write (6,*)
    write (6,*) 'FLUX ANALYSIS:  FULL CELL - REGULAR DIST.'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0

    write(6,'(20a)') 'GFF:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ','  |GPERP   ',&
         '  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ','    E2DION  ',&
         '    RECSRC   ','   NETFLUX   ','   SUM2'
    startflux = e2dgpara(1,ir)

    endflux   = e2dgpara(nks(ir)+1,ir)
    flux_const= (startflux  -endflux + srcinteg(nks(ir)+1,ir)- recinteg(nks(ir)+1,ir))/areasum(ir)

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)

       intgrad = intgrad+ flux_const*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'GFF:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,intgrad,-inte2dhrec,&
            intnflux,flux1,flux2,flux1/flux2,flux_const*ds,e2dion(ik,ir)*ds,e2dhrec(ik,ir)*ds,netflux*ds,flux3,flux3/flux2
       if (ik.eq.endik) then
          mulfact = ((flux2-flux1)/(intgrad)) + 1.0

       endif

       intgrad = intgrad+ flux_const*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       endif
       if (ik.eq.endik) then
          mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0

       endif

    end do

    ik = nks(ir)
    flux1 =  gtarg(ir,2)+inte2dsrc+intgrad-inte2dhrec
    flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
    flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec

    flux4 = e2dgpara(ik+1,ir)

    !     ------------------------------------------------------------------

    write (6,'(a,5g12.5)') 'NUM:',startflux,0.15*intgrad2,flux1,flux2,mulfact

    write(6,'(a,2i4,1p,8(1x,g13.5))') 'GFF:',ir,ik,startflux,inte2dsrc,intgrad,-inte2dhrec,&
         startflux+inte2dsrc+intgrad-inte2dhrec,e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),&
         (startflux+inte2dsrc+intgrad-inte2dhrec)/(e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),flux4

    write (6,*)
    write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact

    write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact

    write(6,'(20a)') 'MGFF:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ','  |GPERP   ',&
         '  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ','    E2DION  ',&
         '    RECSRC   ','   NETFLUX   ','   SUM2'
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0
    startflux = e2dgpara(1,ir)

    endflux   = e2dgpara(nks(ir)+1,ir)

    flux_const= (startflux -endflux + srcinteg(nks(ir)+1,ir)- recinteg(nks(ir)+1,ir))/areasum(ir)

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       ds = ksb(ik,ir)-ksb(ik-1,ir)

       intgrad = intgrad+ flux_const*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1 * mulfact

       intnflux = intnflux+ netflux*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       intdist = intdist+ (kss2(ik,ir)-ksb(ik-1,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))* rmean1

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0


       endif
       flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
       !               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
       flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0

       flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec

       write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MGFF:',ir,ik,kss(ik,ir),ds,startflux,inte2dsrc,&
            intgrad,-inte2dhrec,intnflux,flux1,flux2,flux1/flux2,flux_const*ds,e2dion(ik,ir)*ds,&
            e2dhrec(ik,ir)*ds,netflux*ds,flux3,flux3/flux2

       intgrad = intgrad+ flux_const*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2 * mulfact

       inte2dsrc = inte2dsrc+ e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intnflux = intnflux+ netflux*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       intdist = intdist+ (ksb(ik,ir)-kss2(ik,ir))
       if (gperpa(ik,ir).le.0.0) then
          intgradneg = intgradneg+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2
       else
          intgradpos = intgradpos+ gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))* rmean2

       endif

    end do

    !     ------------------------------------------------------------------
    !     ------------------------------------------------------------------


    !           Analysis of the actual EDGE2D boundary fluxes


    write (6,*)
    write (6,*) 'DOWN FLUX ANALYSIS:'
    write (6,510)

510 format('D-G:',2x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,3x,'D-RAT',4x,2x,'E2DV-RAT')
    intgperpd = 0.0
    intgnet   = 0.0
    intgdiv   = 0.0

    do ik = startik,endik
       fluxst  = e2dgpara(ik,ir)

       fluxend = e2dgpara(ik+1,ir)

       ioniz = e2dion(ik,ir)

       rec   = e2dhrec(ik,ir)

       ds    = ksb(ik,ir) - ksb(ik-1,ir)

       gperpn= -(fluxst - fluxend+ e2dion(ik,ir)*ds - e2dhrec(ik,ir)*ds)

       gperpd= gperpn/ds
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)

       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       gdiv =  0.15*grad_oldknbs(ik,ir)
       intgperpd = intgperpd + ds * gperpd
       intgnet   = intgnet   + ds * netflux

       intgdiv   = intgdiv   + ds * gdiv

       write (6,'(a,2i4,1p,15(1x,g11.4))') 'D-G:',ik,ir,fluxst,ioniz*ds,-rec*ds,gperpn,fluxend,gperpd,netflux,gdiv,intgperpd,intgnet,intgdiv,gdiv/gperpd,netflux/gperpd



       !           Analysis of the actual EDGE2D boundary fluxes


       !     ------------------------------------------------------------------

    end do
    write (6,*) 'DOWN FLUX ACTUAL PARTICLE ANALYSIS: KAREAS'
    write (6,520)

520 format('DPGA:',1x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,3x,'D-RAT',4x,2x,'E2DV-RAT')
    intgperpd = 0.0
    intgnet   = 0.0
    intgdiv   = 0.0

    do ik = startik,endik

       in = korpg(ik,ir)
       sider = (rvertp(1,in) + rvertp(2,in)) /2.0

       fluxst  = e2dgdown(ik,ir)/ (2.0* PI * sider)
       sider = (rvertp(3,in) + rvertp(4,in)) /2.0

       fluxend = e2dgdown(ik+1,ir)/ (2.0 * PI * sider)

       ioniz = e2dion(ik,ir) * kareas(ik,ir)

       rec   = e2dhrec(ik,ir) * kareas(ik,ir)

       ds    = ksb(ik,ir) - ksb(ik-1,ir)

       gperpn= -(fluxst - fluxend+ ioniz - rec)

       gperpd= gperpn/kareas(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)

       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       gdiv =  0.15*grad_oldknbs(ik,ir)
       intgperpd = intgperpd + ds * gperpd
       intgnet   = intgnet   + ds * netflux

       intgdiv   = intgdiv   + ds * gdiv

       write (6,'(a,2i4,1p,15(1x,g11.4))') 'D-G:',ik,ir,fluxst,ioniz,-rec,gperpn,fluxend,gperpd,netflux,gdiv,intgperpd,intgnet,intgdiv,gdiv/gperpd,netflux/gperpd



       !           Analysis of the actual EDGE2D boundary fluxes

       !     ------------------------------------------------------------------

    end do
    write (6,*) 'DOWN FLUX ACTUAL PARTICLE ANALYSIS: E2DAREAS'
    write (6,530)

530 format('DGPE:',1x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,1x,&
         'G-REQ-DEN',2x,2x,'G-VPERP',3x,2x,'G-DIV-E',3x,3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,3x,&
         'D-RAT',4x,2x,'E2DV-RAT',2x,3x,'GDIV-A')
    intgperpd = 0.0
    intgnet   = 0.0
    intgdiv   = 0.0

    do ik = startik,endik

       in = korpg(ik,ir)
       sider = (rvertp(1,in) + rvertp(2,in)) /2.0

       fluxst  = e2dgdown(ik,ir)/ (2.0* PI * sider)
       sider = (rvertp(3,in) + rvertp(4,in)) /2.0

       fluxend = e2dgdown(ik+1,ir)/ (2.0 * PI * sider)

       ioniz = e2dion(ik,ir) * e2dareas(ik,ir)

       rec   = e2dhrec(ik,ir) * e2dareas(ik,ir)

       ds    = ksb(ik,ir) - ksb(ik-1,ir)

       gperpn= -(fluxst - fluxend+ ioniz - rec)

       gperpd= gperpn/e2dareas(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)

       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       gdiv =  0.15*grad_oldknbs(ik,ir)*kareas(ik,ir)/e2dareas(ik,ir)
       intgperpd = intgperpd + ds * gperpd
       intgnet   = intgnet   + ds * netflux

       intgdiv   = intgdiv   + ds * gdiv

       write (6,'(a,2i4,1p,15(1x,g11.4))') 'D-G:',ik,ir,fluxst,ioniz,-rec,gperpn,fluxend,gperpd,netflux,gdiv,intgperpd,&
            intgnet,intgdiv,gdiv/gperpd,netflux/gperpd,0.15*grad_oldknbs(ik,ir)


       !     ------------------------------------------------------------------
       !     ------------------------------------------------------------------

    end do
    intgrad = 0.0
    intgrad2 = 0.0
    intgradpos = 0.0
    intgradneg = 0.0
    inte2dsrc = 0.0
    inte2dhrec = 0.0
    intdist = 0.0

    intnflux = 0.0
    write(6,*) 'Poloidal plane analysis:'

    write(6,'(20a)') 'FLUX:','  IR ',' IK ','      S    ','   DS   ','   GTARG     ',' |E2DION   ','  |GPERP   ',&
         '  |RECSRC     ','|NETFLUX ','     SUM    ','   E2D FLUX    ' , 'RATIO   ','    GPERPD  ','    E2DION  ',&
         '    RECSRC   ','   NETFLUX   ','   SUM2'

    initflux = gtarg(ir,2) / kbfst(ir,2)

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell
       dp = kpb(ik,ir)-kpb(ik-1,ir)
       dp1 = kps2(ik,ir) - kpb(ik-1,ir)

       dp2 = kpb(ik,ir) - kps(ik,ir)


       if (.not.(ik.eq.startik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*dp1* rmean1

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*dp1* rmean1

          intnflux = intnflux+ netflux*dp1* rmean1

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*dp1* rmean1

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*dp1* rmean1

          intdist = intdist+ dp1
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*dp1* rmean1
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*dp1* rmean1

          endif

       endif
       if (switch(swmajr).eq.4.0) then
          rmean1 = rs(ik,ir)
       else
          rmean1 = 1.0

       endif
       flux1 =  initflux+inte2dsrc+0.15*intgrad2-inte2dhrec
       flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1/kbfs(ik,ir)

       flux3 =  initflux+inte2dsrc+intnflux-inte2dhrec

       !               if (ik.eq.endik) then
       !                 mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
       !               endif

       write(6,'(a,2i4,2f9.4,1p,15(1x,g11.4))') 'F:',ir,ik,kps(ik,ir),dp,initflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,&
            intnflux,flux1,flux2,flux1/flux2,0.15 * grad_oldknbs(ik,ir)*dp,e2dion(ik,ir)*dp,e2dhrec(ik,ir)*dp,&
            netflux*dp,flux3,flux3/flux2


       if (.not.(ik.eq.endik.and.(switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.switch(swe2d).eq.3.or.&
            switch(swe2d).eq.4.or.switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then

          intgrad = intgrad+ gperpa(ik,ir)*dp2* rmean2

          intgrad2 = intgrad2+ grad_oldknbs(ik,ir)*dp2* rmean2

          inte2dsrc = inte2dsrc+ e2dion(ik,ir)*dp2* rmean2

          inte2dhrec = inte2dhrec+ e2dhrec(ik,ir)*dp2* rmean2

          intnflux = intnflux+ netflux*dp2* rmean2

          intdist = intdist+ dp2
          if (gperpa(ik,ir).le.0.0) then
             intgradneg = intgradneg+ gperpa(ik,ir)*dp2* rmean2
          else
             intgradpos = intgradpos+ gperpa(ik,ir)*dp2* rmean2

          endif
       endif


       !           Analysis of the actual EDGE2D boundary fluxes - poloidal plane


       !     ------------------------------------------------------------------

    end do
    write (6,*) 'Poloidally calculated:'
    write (6,501)

501 format('BND:',2x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,1x,&
         'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,3x,'D-RAT',4x,2x,&
         'E2DV-RAT')
    intgperpd = 0.0
    intgnet   = 0.0

    intgdiv   = 0.0

    do ik = startik,endik
       if (ik.eq.startik) then
          ! slmod begin - bug fix, 04.03.2010
          !... Correct correction, swapping KBFST for KBFS when IR is used to index the array?
          brat1 = kbfst(ir,2)
          brat2 = (kbfs(ik  ,ir)+kbfs(ik+1,ir))/2.0
       elseif (ik.eq.endik) then
          brat1 = (kbfs(ik-1,ir)+kbfs(ik  ,ir))/2.0
          brat2 = kbfst(ir,1)
       else
          brat1 = (kbfs(ik-1,ir)+kbfs(ik  ,ir))/2.0
          brat2 = (kbfs(ik  ,ir)+kbfs(ik+1,ir))/2.0

          !                  brat2 = (kbfst(ik,ir)+kbfs(ik+1,ir))/2.0
          !               elseif (ik.eq.endik) then
          !                  brat1 = (kbfst(ik-1,ir)+kbfs(ik,ir))/2.0
          !                  brat2 = kbfst(ir,1)
          !               else
          !                  brat1 = (kbfst(ik-1,ir)+kbfs(ik,ir))/2.0
          !                  brat2 = (kbfst(ik,ir)+kbfs(ik+1,ir))/2.0
          !               endif
          ! slmod end
       endif
       fluxst  = e2dflux(ik,ir) /brat1

       fluxend = e2dflux(ik+1,ir) /brat2

       ioniz = e2dion(ik,ir)

       rec   = e2dhrec(ik,ir)

       dp    = kpb(ik,ir) - kpb(ik-1,ir)

       gperpn= -(fluxst - fluxend+ e2dion(ik,ir)*dp - e2dhrec(ik,ir)*dp)

       gperpd= gperpn/dp
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)

       dcell = distin(ik,ir) + distout(ik,ir)

       influx = e2dnbs(ik,ir) * e2dvro(ik,ir)

       outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       netflux = (influx - outflux)/dcell

       gdiv =  0.15*grad_oldknbs(ik,ir)
       intgperpd = intgperpd + dp * gperpd
       intgnet   = intgnet   + dp * netflux

       intgdiv   = intgdiv   + dp * gdiv

       write (6,'(a,2i4,1p,15(1x,g11.4))') 'BND:',ik,ir,fluxst,ioniz*dp,-rec*dp,gperpn,fluxend,gperpd,netflux,gdiv,&
            intgperpd,intgnet,intgdiv,gdiv/gperpd,netflux/gperpd


       !     ------------------------------------------------------------------

    end do

    !           Analyse the EDGE2D perpendicular fluxes and compare them
    !           to the calculated values.

    write (6,*)

    write(6,*) '  FLUX:','  IR ',' IK ','      S     ','      DS     ','     FLUX1   ','     FLUX2   ','      FLUX3 ',&
         '  (F1-F2)/dcell',' (F2-F3)/dcell','  DIVFLUXD   ', '   E2D NBS  ','    VRO'

    do ik = startik,endik
       if (switch(swmajr).eq.4.0) then
          rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
          rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0

       endif
       ikin = ikins(ik,ir)

       irin = irins(ik,ir)
       ikout = ikouts(ik,ir)

       irout = irouts(ik,ir)
       ds = ksb(ik,ir)-ksb(ik-1,ir)

       dcell = distin(ik,ir) + distout(ik,ir)
       flux1 = e2dnbs(ikin,irin) * e2dvro(ikin,irin)
       flux2 = e2dnbs(ik,ir) * e2dvro(ik,ir)

       flux3 = e2dnbs(ikout,irout) * e2dvro(ikout,irout)

       write(6,'(a,2i4,1p,15(1x,g12.4))') 'E2DFLUX:',ir,ik,kss(ik,ir),ds,flux1,flux2,flux3,(flux1-flux2)/dcell,&
            (flux2-flux3)/dcell,0.15*grad_oldknbs(ik,ir),e2dnbs(ik,ir),e2dvro(ik,ir)

       !     ------------------------------------------------------------------


    end do
    return
  end subroutine print_edge2d_flux_analysis



  subroutine ioniz_comp
    use mod_params
    use mod_comtor
    use mod_cgeom
    use mod_pindata
    use mod_cedge2d
    use mod_printopt

    !     include 'params'
    !     include 'comtor'
    !     include 'cgeom'
    !     include 'pindata'
    !     include 'cedge2d'
    !     include 'printopt'

    !     IONIZ_COMP: This routine prints a comparison of the
    !                 ionization data read in from EDGE2D and
    !                 the ionization calculated by NIMBUS.


    implicit none

    integer ir,ik,in
    real ioniz(maxnrs,3,21)

    real sumiz(9,21),sum1,sum2
    call rzero(ioniz,maxnrs*3*21)
    call rzero(sumiz,9*21)
    sum1 = 0.0

    !     Sum up over rings

    sum2 = 0.0
    write(6,*)
    write(6,*)  'RAW DATA:'
    write(6,*)

    write(6,*)  ' IK  IR   PINION  E2DION   PIN-AREA  E2D-AREA'

    DO IR=1,nrs

       do ik = 1,nks(ir)

          write (6,'(2(i4,1x),5g12.4)') ik,ir,pinion(ik,ir),e2dion(ik,ir),karea2(ik,ir),e2dareas(ik,ir)

          if (.not.(ir.lt.irsep.and.ik.eq.nks(ir))) then
             if (ik.le.ikmids(ir)) then
                ioniz(ir,1,1) = ioniz(ir,1,1) +karea2(ik,ir) * pinion(ik,ir)
                ioniz(ir,1,2) = ioniz(ir,1,2) +karea2(ik,ir) * pinion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,1,3) = ioniz(ir,1,3) +karea2(ik,ir) * pinion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,1,4) = ioniz(ir,1,4) +karea2(ik,ir) * e2dion(ik,ir)
                ioniz(ir,1,5) = ioniz(ir,1,5) +karea2(ik,ir) * e2dion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,1,6) = ioniz(ir,1,6) +karea2(ik,ir) * e2dion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,1,7) = ioniz(ir,1,7) +karea2(ik,ir) * pinrec(ik,ir)
                ioniz(ir,1,8) = ioniz(ir,1,8) +karea2(ik,ir) * pinrec(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,1,9) = ioniz(ir,1,9) +karea2(ik,ir) * pinrec(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,1,10) = ioniz(ir,1,10) +e2dareas(ik,ir) * pinion(ik,ir)
                ioniz(ir,1,11) = ioniz(ir,1,11) +e2dareas(ik,ir) * pinion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,1,12) = ioniz(ir,1,12) +e2dareas(ik,ir) * pinion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,1,13) = ioniz(ir,1,13) +e2dareas(ik,ir) * e2dion(ik,ir)
                ioniz(ir,1,14) = ioniz(ir,1,14) +e2dareas(ik,ir) * e2dion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,1,15) = ioniz(ir,1,15) +e2dareas(ik,ir) * e2dion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,1,16) = ioniz(ir,1,16) +e2dareas(ik,ir) * pinrec(ik,ir)
                ioniz(ir,1,17) = ioniz(ir,1,17) +e2dareas(ik,ir) * pinrec(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,1,18) = ioniz(ir,1,18) +e2dareas(ik,ir) * pinrec(ik,ir)*rs(ik,ir)/r0
                if (hcorr(ik,ir).ne.0.0) then 
                   ioniz(ir,1,19) = ioniz(ir,1,19) +karea2(ik,ir) * pinion(ik,ir) / hcorr(ik,ir)
                   ioniz(ir,1,20) = ioniz(ir,1,20) +karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)*r0/rs(ik,ir)
                   ioniz(ir,1,21) = ioniz(ir,1,21) +karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)*rs(ik,ir)/r0
                else
                   ioniz(ir,1,19) = 0.0
                   ioniz(ir,1,20) = 0.0
                   ioniz(ir,1,21) = 0.0


                endif

             else
                ioniz(ir,2,1) = ioniz(ir,2,1) +karea2(ik,ir) * pinion(ik,ir)
                ioniz(ir,2,2) = ioniz(ir,2,2) +karea2(ik,ir) * pinion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,2,3) = ioniz(ir,2,3) +karea2(ik,ir) * pinion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,2,4) = ioniz(ir,2,4) +karea2(ik,ir) * e2dion(ik,ir)
                ioniz(ir,2,5) = ioniz(ir,2,5) +karea2(ik,ir) * e2dion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,2,6) = ioniz(ir,2,6) +karea2(ik,ir) * e2dion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,2,7) = ioniz(ir,2,7) +karea2(ik,ir) * pinrec(ik,ir)
                ioniz(ir,2,8) = ioniz(ir,2,8) +karea2(ik,ir) * pinrec(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,2,9) = ioniz(ir,2,9) +karea2(ik,ir) * pinrec(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,2,10) = ioniz(ir,2,10) +e2dareas(ik,ir) * pinion(ik,ir)
                ioniz(ir,2,11) = ioniz(ir,2,11) +e2dareas(ik,ir) * pinion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,2,12) = ioniz(ir,2,12) +e2dareas(ik,ir) * pinion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,2,13) = ioniz(ir,2,13) +e2dareas(ik,ir) * e2dion(ik,ir)
                ioniz(ir,2,14) = ioniz(ir,2,14) +e2dareas(ik,ir) * e2dion(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,2,15) = ioniz(ir,2,15) +e2dareas(ik,ir) * e2dion(ik,ir)*rs(ik,ir)/r0
                ioniz(ir,2,16) = ioniz(ir,2,16) +e2dareas(ik,ir) * pinrec(ik,ir)
                ioniz(ir,2,17) = ioniz(ir,2,17) +e2dareas(ik,ir) * pinrec(ik,ir)*r0/rs(ik,ir)

                ioniz(ir,2,18) = ioniz(ir,2,18) +e2dareas(ik,ir) * pinrec(ik,ir)*rs(ik,ir)/r0
                if (hcorr(ik,ir).ne.0.0) then 
                   ioniz(ir,2,19) = ioniz(ir,2,19) +karea2(ik,ir) * pinion(ik,ir) / hcorr(ik,ir)
                   ioniz(ir,2,20) = ioniz(ir,2,20) +karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)*r0/rs(ik,ir)
                   ioniz(ir,2,21) = ioniz(ir,2,21) +karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)*rs(ik,ir)/r0
                else
                   ioniz(ir,2,19) = 0.0
                   ioniz(ir,2,20) = 0.0
                   ioniz(ir,2,21) = 0.0

                endif

             endif

          endif
       end do

       !     Summations

    end do

    do ir = 1,nrs
       do in = 1,21
          ioniz(ir,3,in) = ioniz(ir,1,in) + ioniz(ir,2,in)

          !        Core totals

       end do
       if (ir.lt.irsep) then
          if (ir.ge.ircent) then

             do in = 1,21

                sumiz(1,in) = sumiz(1,in) + ioniz(ir,3,in)
             end do
          else

             do in = 1,21

                sumiz(9,in) = sumiz(9,in) + ioniz(ir,3,in)
             end do

             !        Main SOL totals

          endif

       elseif (ir.ge.irsep.and.ir.lt.irwall) then

          do in = 1,21
             sumiz(2,in) = sumiz(2,in) + ioniz(ir,3,in)
             sumiz(3,in) = sumiz(3,in) + ioniz(ir,1,in)

             sumiz(4,in) = sumiz(4,in) + ioniz(ir,2,in)

             !        PP totals

          end do

       elseif (ir.gt.irtrap) then

          do in = 1,21
             sumiz(5,in) = sumiz(5,in) + ioniz(ir,3,in)
             sumiz(6,in) = sumiz(6,in) + ioniz(ir,1,in)

             sumiz(7,in) = sumiz(7,in) + ioniz(ir,2,in)

          end do
       endif

       !     Grand totals

    end do
    do in = 1,21
       sumiz(8,in) = sumiz(1,in) + sumiz(2,in) + sumiz(5,in)

       !     Print out summaries for ALL

    end do
    write (6,*)
    write (6,*) 'Ionization SUMMARY - PIN Areas:'
    write (6,*)

    write (6,100)
    write (6,110) 'PIN     :',(sumiz(in,1),in=1,8),1.0

    if (sumiz(8,1).ne.0.0) then 
       write (6,110) 'PIN R0/R:',(sumiz(in,2),in=1,8),sumiz(8,2)/sumiz(8,1)
       write (6,110) 'PIN R/R0:',(sumiz(in,3),in=1,8),sumiz(8,3)/sumiz(8,1)
    endif
    if (sumiz(8,1).ne.0.0.and.sumiz(8,2).ne.0.0.and.sumiz(8,3).ne.0.0) then
       write (6,110) 'E2D     :',(sumiz(in,4),in=1,8),sumiz(8,4)/sumiz(8,1),sumiz(8,4)/sumiz(8,2),sumiz(8,4)/sumiz(8,3)
       write (6,110) 'E2D R0/R:',(sumiz(in,5),in=1,8),sumiz(8,5)/sumiz(8,1),sumiz(8,5)/sumiz(8,2),sumiz(8,5)/sumiz(8,3)
       write (6,110) 'E2D R/R0:',(sumiz(in,6),in=1,8),sumiz(8,6)/sumiz(8,1),sumiz(8,6)/sumiz(8,2),sumiz(8,6)/sumiz(8,3)
    endif
    write (6,110) 'REC     :',(sumiz(in,7),in=1,8),1.0
    if (sumiz(8,7).ne.0.0) then 
       write (6,110) 'REC R0/R:',(sumiz(in,8),in=1,8),sumiz(8,8)/sumiz(8,7)
       write (6,110) 'REC R/R0:',(sumiz(in,9),in=1,8),sumiz(8,9)/sumiz(8,7)
    endif

    write(6,*)
    write (6,*)
    write (6,*) 'Ionization SUMMARY - E2D Areas - Estimated:'
    write (6,*)

    write (6,100)
    write (6,110) 'PIN     :',(sumiz(in,10),in=1,8),1.0
    if (sumiz(8,10).ne.0.0) then 
       write (6,110) 'PIN R0/R:',(sumiz(in,11),in=1,8),sumiz(8,11)/sumiz(8,10)
       write (6,110) 'PIN R/R0:',(sumiz(in,12),in=1,8),sumiz(8,12)/sumiz(8,10)
    endif
    if (sumiz(8,10).ne.0.0.and.sumiz(8,11).ne.0.0.and.sumiz(8,12).ne.0.0) then
       write (6,110) 'E2D     :',(sumiz(in,13),in=1,8),sumiz(8,13)/sumiz(8,10),sumiz(8,13)/sumiz(8,11),sumiz(8,13)/sumiz(8,12)
       write (6,110) 'E2D R0/R:',(sumiz(in,14),in=1,8),sumiz(8,14)/sumiz(8,10),sumiz(8,14)/sumiz(8,11),sumiz(8,14)/sumiz(8,12)
       write (6,110) 'E2D R/R0:',(sumiz(in,15),in=1,8),sumiz(8,15)/sumiz(8,10),sumiz(8,15)/sumiz(8,11),sumiz(8,15)/sumiz(8,12)
    endif
    write (6,110) 'REC     :',(sumiz(in,16),in=1,8),1.0
    if (sumiz(8,16).ne.0.0) then 
       write (6,110) 'REC R0/R:',(sumiz(in,17),in=1,8),sumiz(8,17)/sumiz(8,16)
       write (6,110) 'REC R/R0:',(sumiz(in,18),in=1,8),sumiz(8,18)/sumiz(8,16)
    endif
    write(6,*)
    write (6,*) 'Ionization SUMMARY - PIN Areas / HCORR :'
    write (6,*)
    write (6,100)
    write (6,110) 'PIN     :',(sumiz(in,19),in=1,8),1.0
    if (sumiz(8,19).ne.0.0) then 
       write (6,110) 'PIN R0/R:',(sumiz(in,20),in=1,8),sumiz(8,20)/sumiz(8,19)
       write (6,110) 'PIN R/R0:',(sumiz(in,21),in=1,8),sumiz(8,21)/sumiz(8,19)
    endif
    if (sumiz(8,19).ne.0.0.and.sumiz(8,20).ne.0.0.and.sumiz(8,21).ne.0.0) then
       write (6,110) 'E2D     :',(sumiz(in,13),in=1,8),sumiz(8,13)/sumiz(8,19),sumiz(8,13)/sumiz(8,20),sumiz(8,13)/sumiz(8,21)
       write (6,110) 'E2D R0/R:',(sumiz(in,14),in=1,8),sumiz(8,14)/sumiz(8,19),sumiz(8,14)/sumiz(8,20),sumiz(8,14)/sumiz(8,21)
       write (6,110) 'E2D R/R0:',(sumiz(in,15),in=1,8),sumiz(8,15)/sumiz(8,19),sumiz(8,15)/sumiz(8,20),sumiz(8,15)/sumiz(8,21)


    endif
    write (6,*)
    write (6,'(a,1x,3g12.4)') 'PIN I-Core Ioniz (D-A):',sumiz(9,1),sumiz(9,2),sumiz(9,3)
    write (6,'(a,1x,3g12.4)') 'PIN I-Core Recom (D-A):',sumiz(9,7),sumiz(9,8),sumiz(9,9)
    write (6,'(a,1x,3g12.4)') 'PIN I-Core Ioniz (D-E):',sumiz(9,10),sumiz(9,11),sumiz(9,12)
    write (6,'(a,1x,3g12.4)') 'PIN I-Core Recom (D-E):',sumiz(9,16),sumiz(9,17),sumiz(9,18)

    write (6,'(a,1x,3g12.4)') 'E2D Inner core Ioniz:',sumiz(9,4)
    write (6,*)
    write (6,*) 'Grand Total Ionizations:'
    write (6,*)
    write (6,*) 'Polygon Areas:'
    write (6,*)
    write (6,110) 'PIN     :',sumiz(8,1) + sumiz(9,1)
    write (6,110) 'PIN R0/R:',sumiz(8,2) + sumiz(9,2)
    write (6,110) 'PIN R/R0:',sumiz(8,3) + sumiz(9,3)
    write (6,110) 'E2D     :',sumiz(8,4) + sumiz(9,4)
    write (6,110) 'E2D R0/R:',sumiz(8,5) + sumiz(9,5)
    write (6,110) 'E2D R/R0:',sumiz(8,6) + sumiz(9,6)
    write (6,*)
    write (6,*) 'E2D Areas:'
    write (6,*)
    write (6,110) 'PIN     :',sumiz(8,10) + sumiz(9,10)
    write (6,110) 'PIN R0/R:',sumiz(8,11) + sumiz(9,11)
    write (6,110) 'PIN R/R0:',sumiz(8,12) + sumiz(9,12)
    write (6,110) 'E2D     :',sumiz(8,13) + sumiz(9,13)
    write (6,110) 'E2D R0/R:',sumiz(8,14) + sumiz(9,14)
    write (6,110) 'E2D R/R0:',sumiz(8,15) + sumiz(9,15)

    write (6,*)
100 format('TYPE',5x,3x,'O-Core',3x,2x,'SOL Total',1x,2x,'SOL Outer',1x,2x,'SOL Inner',1x,2x,'PP Total',2x,2x,'PP Outer',&
         2x,2x,'PP Inner',2x,1x,'Grand Total',6x,'Ratios')
110 format(a8,1x,8g12.4,3f9.5)
120 format('TYPE',5x,4x,a5,3x,4x,a5,3x,4x,'Total',6x,'Ratio')
    write(6,*)
    write(6,*) 'Ionization Summary by ring'

    write(6,*)

    write(6,120) outer,inner
    sum1 = 1.0

    sum2 = 0.0

    do ir = 1,nrs
       write(6,*)
       write(6,*) 'RING = ',ir
       write(6,*)
       write (6,*)
       write (6,*) '     PIN/DIV - Areas'
       write (6,*)
       write (6,110) 'PIN     :',(ioniz(ir,in,1),in=1,3),1.0
       if (ioniz(ir,3,1).ne.0.0) then 
          write (6,110) 'PIN R0/R:',(ioniz(ir,in,2),in=1,3),ioniz(ir,3,2)/ioniz(ir,3,1)
          write (6,110) 'PIN R/R0:',(ioniz(ir,in,3),in=1,3),ioniz(ir,3,3)/ioniz(ir,3,1)
       endif
       if (ioniz(ir,3,1).ne.0.0.and.ioniz(ir,3,2).ne.0.0.and.ioniz(ir,3,3).ne.0.0) then
          write (6,110) 'E2D     :',(ioniz(ir,in,4),in=1,3),ioniz(ir,3,4)/ioniz(ir,3,1),ioniz(ir,3,4)/ioniz(ir,3,2),&
               ioniz(ir,3,4)/ioniz(ir,3,3)
          write (6,110) 'E2D R0/R:',(ioniz(ir,in,5),in=1,3),ioniz(ir,3,5)/ioniz(ir,3,1),ioniz(ir,3,5)/ioniz(ir,3,2),&
               ioniz(ir,3,5)/ioniz(ir,3,3)
          write (6,110) 'E2D R/R0:',(ioniz(ir,in,6),in=1,3),ioniz(ir,3,6)/ioniz(ir,3,1),ioniz(ir,3,6)/ioniz(ir,3,2),&
               ioniz(ir,3,6)/ioniz(ir,3,3)
       endif
       write (6,110) 'REC     :',(ioniz(ir,in,7),in=1,3),1.0
       if (ioniz(ir,3,7).ne.0.0) then 
          write (6,110) 'REC R0/R:',(ioniz(ir,in,8),in=1,3),ioniz(ir,3,8)/ioniz(ir,3,7)
          write (6,110) 'REC R/R0:',(ioniz(ir,in,9),in=1,3),ioniz(ir,3,9)/ioniz(ir,3,7)
       endif
       write (6,*)
       write (6,*) '     E2D - Areas         IRCENT=',ircent
       write (6,*)
       if (ir.ge.ircent.and.ir.ne.irtrap.and.ir.ne.irwall) then
          sum1 = sum1 + ioniz(ir,3,10)
          sum2 = sum2 + ioniz(ir,3,13)
       endif
       write (6,110) 'PIN     :',(ioniz(ir,in,10),in=1,3),1.0
       if (ioniz(ir,3,10).ne.0.0) then 
          write (6,110) 'PIN R0/R:',(ioniz(ir,in,11),in=1,3),ioniz(ir,3,11)/ioniz(ir,3,10)
          write (6,110) 'PIN R/R0:',(ioniz(ir,in,12),in=1,3),ioniz(ir,3,12)/ioniz(ir,3,10)
       endif
       if (ioniz(ir,3,10).ne.0.0.and.ioniz(ir,3,11).ne.0.0.and.ioniz(ir,3,12).ne.0.0) then
          write (6,110) 'E2D     :',(ioniz(ir,in,13),in=1,3),ioniz(ir,3,13)/ioniz(ir,3,10),ioniz(ir,3,13)/ioniz(ir,3,11),&
               ioniz(ir,3,13)/ioniz(ir,3,12)
          write (6,110) 'E2D R0/R:',(ioniz(ir,in,14),in=1,3),ioniz(ir,3,14)/ioniz(ir,3,10),ioniz(ir,3,14)/ioniz(ir,3,11),&
               ioniz(ir,3,14)/ioniz(ir,3,12)
          write (6,110) 'E2D R/R0:',(ioniz(ir,in,15),in=1,3),ioniz(ir,3,15)/ioniz(ir,3,10),ioniz(ir,3,15)/ioniz(ir,3,11),&
               ioniz(ir,3,15)/ioniz(ir,3,12)
       endif
       write (6,110) 'REC     :',(ioniz(ir,in,16),in=1,3),1.0
       if (ioniz(ir,3,16).ne.0.0) then 
          write (6,110) 'REC R0/R:',(ioniz(ir,in,17),in=1,3),ioniz(ir,3,17)/ioniz(ir,3,16)
          write (6,110) 'REC R/R0:',(ioniz(ir,in,18),in=1,3),ioniz(ir,3,18)/ioniz(ir,3,16)
       endif
       write(6,*)
       if (sum1.ne.0.0) then 
          write (6,110) 'RUNTOT:',sum1,sum2,sum2/sum1

       endif

    end do
    write (6,*)
    write (6,*) 'EDGE 2D without boundary rings:',sumiz(8,4)-ioniz(1,3,4)-ioniz(irwall,3,4)-ioniz(irtrap,3,4)
    write (6,*)
    write (6,*) 'Total Ratios - PIN/DIV Areas:'
    write (6,*)
    if (sumiz(8,1).ne.0.0) then 
       write (6,'(a,g12.4)')'E2D  /DIVIMP = ', sumiz(8,4)/sumiz(8,1)
       write (6,'(a,g12.4)')'E2DR1/DIVIMP = ', sumiz(8,5)/sumiz(8,1)
       write (6,'(a,g12.4)')'E2DR2/DIVIMP = ', sumiz(8,6)/sumiz(8,1)
    endif
    write(6,*)
    write (6,*)
    write (6,*) 'Total Ratios - E2D Areas:'
    write (6,*)
    if (sumiz(8,10).ne.0.0) then 
       write (6,'(a,g12.4)')'E2D  /DIVIMP = ', sumiz(8,13)/sumiz(8,10)
       write (6,'(a,g12.4)')'E2DR1/DIVIMP = ', sumiz(8,14)/sumiz(8,10)
       write (6,'(a,g12.4)')'E2DR2/DIVIMP = ', sumiz(8,15)/sumiz(8,10)
    endif

    write(6,*)


    write (6,*)
    return



  end subroutine ioniz_comp



  subroutine detached_plasma(spts,npts,errcode,serr,te,ti,ne,vb,ir,swtarg,te0,ti0,n0,v0,act_press,mb)
    use mod_params
    use mod_cgeom
    use mod_comtor
    use mod_solparams

    !     include 'params'
    !     include 'cgeom'
    !     include 'comtor'

    implicit none
    integer errcode,npts,ir

    !     include 'solparams'

    real    swtarg
    real*8 serr,te0,ti0,n0,v0,mb
    real*8 spts(npts)

    !     Local variables

    real*8 te(npts),ti(npts),ne(npts),vb(npts),act_press(npts)
    real*8 smax,pmax,news,lppa,lprad
    real*8 t1,ti1,t2,n1,v1,s,sl1,sl2,slv,na

    real*8 sl1a,sl1b

    !     Extra parameters - have no effect unless specified in
    !     for S21            additional per ring data sets.

    !     l1a ... allow for more detail in the description of the
    !              density evolution in region A - three linear
    !              fitted regions.

    double precision l1r,l1ri,l2r,l2ri,lvr,lvri,ter,teri,tir,tiri,nr,nri,vbm,vbmi,qr,qri,n_exp,n_expi

    double precision l1a,l1ai,l1b,l1bi,nr1a,nr1ai,nr1b,nr1bi,ter1a,ter1ai,tir1a,tir1ai,ter1b,ter1bi,tir1b,tir1bi

    integer ik
    serr = 0.0

    !     Assign lengths

    errcode = 0
    smax = ksmaxs(ir)

    !     Set up scaling lengths

    pmax = kpmaxs(ir)

    !     Note : OUTER target values are loaded by default - only overwritten
    !            if INNER target is being calculated.

    call load_s21params(ir,l1r,l1ri,l2r,l2ri,lvr,lvri,ter,teri,tir,tiri,nr,nri,vbm,vbmi,qr,qri,n_exp,n_expi,l1a,nr1a,&
         l1ai,nr1ai,l1b,nr1b,l1bi,nr1bi,ter1a,ter1ai,tir1a,tir1ai,ter1b,ter1bi,tir1b,tir1bi)

    !        Outer target - ik = 1

    !         l1r = l1r
    !         l2r = l2r
    !         lvr = lvr

    !         tir = tir
    !         ter = ter

    !        Pressure matching not allowed for SOL 21 option included in
    !        SOL 22.

    !         nr  = abs(nr)

    if (swtarg.eq.1.0) then


       !         qr  = qr

       na  = n_exp

       !        Inner target - ik = nks(ir)

    elseif (swtarg.eq.2.0) then
       l1a = l1ai
       l1b = l1bi
       l1r = l1ri
       l2r = l2ri
       lvr = lvri
       tir = tiri

       ter = teri
       ter1a = ter1ai 
       ter1b = ter1bi 
       tir1a = tir1ai 


       !        Pressure matching not allowed for SOL 21 option included in
       !        SOL 22.

       tir1b = tir1bi 
       nr  = abs(nri)
       na  = n_expi
       nr1a = abs(nr1ai)
       nr1b = abs(nr1bi)
       qr  = qri

       !     Length reference switch for SOL 21 - options

       !        0 = relative to SMAX
       !        1 = relative to PMAX
       !        2 = absolute SMAX units
       !        3 = absolute PMAX units

    endif

    if (s21refsw.eq.0) then
       sl1a = l1a * smax
       sl1b = l1b * smax
       sl1  = l1r * smax
       sl2  = l2r * smax

       slv  = lvr * smax

    elseif (s21refsw.eq.1) then
       sl1a = l1a * pmax
       sl1b = l1b * pmax
       sl1  = l1r * pmax
       sl2  = l2r * pmax

       slv  = lvr * pmax

    elseif (s21refsw.eq.2.or.s21refsw.eq.3) then
       sl1a = l1a 
       sl1b = l1b 
       sl1  = l1r
       sl2  = l2r

       slv  = lvr

       !     Convert P-coordinates to S-coordinates

    endif

    if (s21refsw.eq.1.or.s21refsw.eq.3) then
       call cnvrtptos(sl1a,news,ir)

       sl1a = news
       call cnvrtptos(sl1b,news,ir)

       sl1b = news
       call cnvrtptos(sl1,news,ir)

       sl1 = news
       call cnvrtptos(sl2,news,ir)

       sl2 = news
       call cnvrtptos(slv,news,ir)

       slv = news

    endif

    lppa = (5.0+2.0*ti0/te0+15.0/te0)*n0*abs(v0)*te0*1.602192e-19

    lprad = qrat * lppa / (sl2-sl1)

    t1    = ter * te0

    ti1   = tir * ti0

    n1    = nr * n0

    t2    = (t1**3.5+7.0/(2.0*ck0)*(lppa*(sl2-sl1)+0.5*(sl2-sl1)**2*lprad))**(2.0/7.0)

    !      write (6,*) 'SOL22/21o:',sl1,sl2,slv,smax,lppa,lprad
    !      write (6,*) 'SOL22/21b:',t1,ti1,n1,t2,v1

    !     Calculate the background for the given S-values

    v1    = (n0 * v0) / n1

    do ik = 1,npts

       !        Test S-value and then calculate Te,Ti,N and V

       s = spts(ik)

       if (s.lt.sl1a) then

          te(ik) = te0 + (te0*ter1a-te0) * (s/sl1a)

          ti(ik) = ti0 + (ti0*tir1a-ti0) * (s/sl1a)

          ne(ik)  = n0 + (n0*nr1a-n0)  * (s/sl1a)**na

          vb(ik)  = (n0 * v0)/ ne(ik)

       elseif (s.lt.sl1b) then

          te(ik) = te0*ter1a + (te0*ter1b-te0*ter1a)* (s-sl1a)/(sl1b-sl1a)

          ti(ik) = ti0*tir1a + (ti0*tir1b-ti0*tir1a)* (s-sl1a)/(sl1b-sl1a)

          ne(ik)  = n0*nr1a + (n0*nr1b-n0*nr1a)* ((s-sl1a)/(sl1b-sl1a))**na

          vb(ik)  = (n0 * v0)/ ne(ik)

       elseif (s.le.sl1) then

          te(ik) = te0*ter1b + (t1-te0*ter1b)* (s-sl1b)/(sl1-sl1b)

          !            ne(ik)  = n0 + (n1-n0)  * (s/sl1)**na

          ti(ik) = ti0*tir1b + (ti1-ti0*tir1b)* (s-sl1b)/(sl1-sl1b)

          ne(ik)  = n0*nr1b + (n1-n0*nr1b)* ((s-sl1b)/(sl1-sl1b))**na

          vb(ik)  = (n0 * v0)/ ne(ik)

       elseif (s.le.sl2) then

          te(ik) = (t1**3.5+7.0/(2.0*ck0)*(lppa*(s-sl1)+0.5*(s-sl1)**2*lprad))**(2.0/7.0)

          ti(ik) = te(ik)

          ne(ik)  = (n1*t1) / te(ik)
          if (s.le.slv) then
             vb(ik)= v1 * (slv-s)/(slv-sl1)
          else
             vb(ik)= 0.0

          endif

       else

          te(ik) = (t2**3.5+7.0/(2.0*ck0)*((1.0+qr)*lppa*(s-sl2)))**(2.0/7.0)

          ti(ik) = te(ik)

          ne(ik)  = (n1*t1) / te(ik)
          if (s.le.slv) then
             vb(ik)= v1 * (slv-s)/(slv-sl1)
          else
             vb(ik)= 0.0

          endif

          !        Calculate pressure

       endif

       act_press(ik) = ne(ik) * econv * (te(ik) + ti(ik))+ mb * mconv * ne(ik) * vb(ik)**2

       !     End of Detached Plasma calculation

    end do
    return



  end subroutine detached_plasma







end module mod_sol22_divimp
