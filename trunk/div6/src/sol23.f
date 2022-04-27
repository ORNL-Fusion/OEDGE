      subroutine sol23_interface (irlim1,irlim2)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_sol23_input
      use mod_pindata
      use mod_pin_cfd
      use mod_sol23_com
      implicit none
      integer irlim1, irlim2
c
c     SOL23_INTERFACE:
c
c     This subroutine provides the interface to the CFD_OSM solver
c     developed by Wojciech Fundamenski. It loads source terms
c     and boundary conditions into a set of one-dimensional arrays
c     which are then passed to the solver which solves the fluid
c     equations. This series is called iteratively - on the first
c     iteration the solver uses some assumptions to generate a seed
c     plasm solution.
c
c     David Elder             March 21, 1997
c
c
c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c     include 'sol23_input'
c     include 'pindata'
c     include 'pin_cfd'
c     include 'sol23_com'
c
      real ringpow(maxnrs)
      real powfact
c
      integer ir,ik,midnks,opt,pin_count,pin_temp,Te_flag
c
      data pin_count /0/ 
      data pin_temp /0/ 
c
c     Integration functions
c
      real*8 src_int, src_int_ab, Q_CX
      external src_int, src_int_ab, Q_CX
      real*8 srcsum1, srcsum2
c
c     Variables for calculating the electric field
c
      real ds1,ds2,nb1,nb2,dt1,dt2,dp1,dp2,start_time,
     >     end_time,elapsed_time,fract2_last,fract3_last,
     >     zhi_1, zhi_2, zhi_3, zhi_4, zhi_n, temp_1, temp_2,
     >     temp_3, pinqe_mult, fract2_mp, temp_4, temp_5,
     >     temp_6, temp_7, temp_8, temp_9, power_ave,
     >     zhi_temp(maxnrs), calcdist, calcwidth2,
     >     Mach_a, Mach_b, Mach_c, Mach_d, f_rec_0, f_rec_L
      integer j, dummy_1, err_dummy, power_sector, i_a, i_b,
     >         i_c, i_d
      real     za02as
      external za02as
c
      err_dummy  = 0
      start_time = za02as(1)
c
c     Set DIVIMP call - 0 for DIVIMP call - 1 for stand-alone
c
c      write (6,*) 'SOL23 Start:',irlim1,irlim2
c
c     Set powfact temporarily here - may be parameterized later
c
      powfact = 1.0
c
c     INTOPT -> Integration option - controls the calculation
c     method used for calculating integrals
c
      intopt           = sol23_intopt
      s23_adaptgrid    = sol23_adaptgrid
      s23_maxtol       = sol23_maxtol
      s23_rmstol       = sol23_rmstol
      s23_perp         = sol23_perp
      s23_relax        = sol23_relax
      s23_symm         = 0
      s23_2T           = 1
c
      s23_par_ptipte   = sol23_par_ptipte
      s23_par_adaptnr  = sol23_par_adaptnr
      s23_par_debugflag= sol23_par_debugflag
      s23_par_debugnr  = sol23_par_debugnr
      s23_par_refresh  = sol23_par_refresh
      s23_par_artvisc2 = sol23_par_artvisc2
      s23_par_artvisc4 = sol23_par_artvisc4
      s23_par_dtf      = sol23_par_dtf
      s23_par_dtg      = sol23_par_dtg
      s23_par_grid0    = sol23_par_grid0
      s23_par_gridexp  = sol23_par_gridexp
      s23_par_itermax  = sol23_par_itermax
      s23_par_ga1      = sol23_par_ga1
      s23_par_ga2      = sol23_par_ga2
      s23_par_ga3      = sol23_par_ga3
      s23_par_ga4      = sol23_par_ga4
      s23_par_updtqpit = sol23_par_updtqpit
      s23_par_updtqp2  = sol23_par_updtqp2
      s23_par_updtqp3  = sol23_par_updtqp3
      s23_par_updtqpte = sol23_par_updtqpte
      s23_par_garelax  = sol23_par_garelax
      s23_par_gaiter   = sol23_par_gaiter
      s23_par_limitte  = sol23_par_limitte
      s23_par_celldte  = sol23_par_celldte
      s23_par_updtbcte = sol23_par_updtbcte
      s23_par_updtdel0 = sol23_par_updtdel0
      s23_par_updtdel1 = sol23_par_updtdel1
      s23_par_updtdelm = sol23_par_updtdelm
      s23_par_qbrelax  = sol23_par_qbrelax
      s23_par_gridg    = sol23_par_gridg
      s23_par_grid_dx0 = sol23_par_grid_dx0
      s23_par_gae      = sol23_par_gae
      s23_par_gai      = sol23_par_gai
      s23_par_tectrl   = sol23_par_tectrl
      s23_par_drflag   = sol23_par_drflag
      s23_par_dsflag   = sol23_par_dsflag
      s23_par_limrel1  = sol23_par_limrel1
      s23_par_limrel2  = sol23_par_limrel2
      s23_par_limrel3  = sol23_par_limrel3
      s23_par_limrel4  = sol23_par_limrel4
      s23_par_tulimit  = sol23_par_tulimit
      s23_par_g0relax  = sol23_par_g0relax
      s23_par_p0relax  = sol23_par_p0relax
      s23_par_dpdrflag = sol23_par_dpdrflag
      s23_par_dpdrtemin= sol23_par_dpdrtemin
      s23_par_dpdrstep = sol23_par_dpdrstep
      s23_par_nulimit  = sol23_par_nulimit
      s23_par_nuflag   = sol23_par_nuflag
      s23_par_nulimit  = sol23_par_nulimit
      s23_par_vnmult   = sol23_par_vnmult
      s23_par_pinmom   = sol23_par_pinmom
      s23_par_emolec   = sol23_par_emolec
      s23_par_rec_heat = sol23_par_rec_heat
      s23_par_pinqimult= sol23_par_pinqimult
      s23_par_pinqiflag= sol23_par_pinqiflag
      s23_par_prring0  = sol23_par_prring0
      s23_par_prring1  = sol23_par_prring1
      s23_par_qperp34  = sol23_par_qperp34
      s23_par_qeiflag  = sol23_par_qeiflag
      s23_par_chie     = sol23_par_chie
      s23_par_joule    = sol23_par_joule
      s23_par_fluxexp  = sol23_par_fluxexp
      s23_par_qrec     = sol23_par_qrec
      s23_par_fzrad    = sol23_par_fzrad
      s23_par_dvmrel   = sol23_par_dvmrel
c
c      write(71,'(1G12.6)')
c     >  s23_par_ptipte,
c     >  s23_par_adaptnr , s23_par_debugflag, s23_par_debugnr,
c     >  s23_par_refresh , s23_par_artvisc2 , s23_par_artvisc4,
c     >  s23_par_dtf     , s23_par_dtg      , s23_par_grid0  ,
c     >  s23_par_gridexp , s23_par_itermax  , s23_par_ga1    ,
c     >  s23_par_ga2     , s23_par_ga3      , s23_par_ga4    ,
c     >  s23_par_updtqpit, s23_par_updtqp2  , s23_par_updtqp3,
c     >  s23_par_updtqpte, s23_par_garelax  , s23_par_gaiter ,
c     >  s23_par_limitte , s23_par_celldte  , s23_par_updtbcte,
c     >  s23_par_updtdel0, s23_par_updtdel1 , s23_par_updtdelm,
c     >  s23_par_qbrelax , s23_par_gridg    , s23_par_grid_dx0,
c     >  s23_par_gae     , s23_par_gai      , s23_par_tectrl ,
c     >  s23_par_drflag  , s23_par_dsflag   , s23_par_limrel1,
c     >  s23_par_limrel2 , s23_par_limrel3  , s23_par_limrel4,
c     >  s23_par_tulimit , s23_par_g0relax  , s23_par_p0relax,
c     >  s23_par_dpdrflag, s23_par_dpdrtemin, s23_par_dpdrstep,
c     >  s23_par_nuflag  , s23_par_nulimit  , s23_par_vnmult ,
c     >  s23_par_pinmom  , s23_par_emolec   , s23_par_rec_heat,
c     >  s23_par_pinqimult,s23_par_pinqiflag, s23_par_prring0,
c     >  s23_par_prring1 , s23_par_qperp34  , s23_par_qeiflag,
c     >  s23_par_chie    , s23_par_joule    , s23_par_fluxexp ,
c     >  s23_par_qrec    , s23_par_fzrad    , s23_par_dvmrel,
c     >  intopt          , s23_adaptgrid    , sol23_seed,
c     >  sol23_izlen     , sol23_izlam      , sol23_izoffset,
c     >  sol23_momlen    , s23_relax        , s23_maxtol,
c     >  s23_rmstol      , s23_perp
c
      mb = crmb
      cb = cizb
      gae_const  = s23_par_gae
      gai_const  = s23_par_gai
      s23_zhi0   = 1.0 / (1.0 + s23_par_fzrad)
c
c      do ir = irlim1, irlim2
c         write(71,'(A10,5G12.6)') 'TARGN', kteds(idds(ir,1)),
c     >         ktids(idds(ir,1)), knds(idds(ir,1)),
c     >         kvds(idds(ir,1)), cmachno(ir,1)
c       enddo
c       write(71,*)
c
c       do ir = irlim1, irlim2
c         write(71,'(A10,5G12.6)') 'TARG0', kteds(idds(ir,2)),
c     >         ktids(idds(ir,2)), knds(idds(ir,2)),
c     >         kvds(idds(ir,2)), cmachno(ir,2)
c      enddo
c
c     IF PIN has not been run then we must use the seedplasma in the
c     solver. Eventually seedplasma will have a specified value in the
c     input file - for now just set it to 1.
c
c     Use seedplasma to determine if it is necessary to load sources
c
c
      if (piniter) then
         pin_count  = pin_count + 1
         pin_temp   = pin_temp + 1
         seedplasma = 0
      else
         pin_count  = 0
         pin_temp   = 0
         seedplasma = sol23_seed
c
         do ir = irlim1, irlim2
            finished(ir)    = 0
            call calctargfluxes(ir)
            gamma_Lang(ir,1) = gtarg(1)
            gamma_Lang(ir,2) = gtarg(2)
            gamma_Lang(ir,3) = gtarg(1) + gtarg(2)
            Te_Lang(ir,1) = kteds(idds(ir,2))
            Te_Lang(ir,2) = kteds(idds(ir,1))
            Te_ctrl(ir,1) = Te_Lang(ir,1)
            Te_ctrl(ir,2) = Te_Lang(ir,2)
            mom_Lang(ir,1)   = momtarg(1)
            mom_Lang(ir,2)   = momtarg(2)
            mom_Lang(ir,3)   = momtarg(2) + momtarg(1)
            poweri_Lang(ir,1) = ionptarg(1)
            poweri_Lang(ir,2) = ionptarg(2)
            poweri_Lang(ir,3) = (ionptarg(1) + ionptarg(2))
            powere_Lang(ir,1) = elecptarg(1)
            powere_Lang(ir,2) = elecptarg(2)
            powere_Lang(ir,3) = (elecptarg(1) + elecptarg(2))
            power_Lang(ir,1) = poweri_Lang(ir,1) + powere_Lang(ir,1)
            power_Lang(ir,2) = poweri_Lang(ir,2) + powere_Lang(ir,2)
            power_Lang(ir,3) = poweri_Lang(ir,3) + powere_Lang(ir,3)
         enddo
c
         call readcfdbg(irlim1,irlim2)
c
         if (s23_par_tectrl.eq.0) then
           do ir = irlim1, irlim2
             Te_ctrl(ir,1) = Te_Lang(ir,1)
             Te_ctrl(ir,2) = Te_Lang(ir,2)
           enddo
         endif
c
      endif
      pin_iter_no = pin_count
c
      write(72,*) 'PIN_COUNT:',pin_count,'OUT OF:',nitersol
      write(71,*) 'PIN_COUNT:',pin_count,'OUT OF:',nitersol
c      write(0,*) 'PIN_COUNT:',pin_count,'OUT OF:',nitersol
      write(71,*) 'FILE_CFD:',s23_cfdfile
      write(71,*) 'SEEDPLASMA:', seedplasma
c
      do ir = irlim1, irlim2
         if (ir.eq.14) then
           write(71,*)
           write(71,*) 'PLASMA:'
           write(71,*)
           write(71,'(A10,6A16)') 'ir   ik', 'n [m-3]  ',
     >      'Te [eV]  ','Ti [eV]  ', 'nv [m-2s-1]',
     >      'Qrec [m-3s-1]', 'Qcx [m-3s-1]'
           write(71,*)
         endif
c
         do ik = 1,nks(ir)
           knbs(ik,ir)  = oldknbs(ik,ir)
           ktebs(ik,ir) = oldktebs(ik,ir)
           ktibs(ik,ir) = oldktibs(ik,ir)
           kvhs(ik,ir)  = oldkvhs(ik,ir)
           if (ir.eq.14) then
             write(71,'(2I5,6G16.8)') ir, ik,
     >         knbs(ik,ir), ktebs(ik,ir), ktibs(ik,ir),
     >         kvhs(ik,ir) * knbs(ik,ir),
     >         pinrec(ik,ir), - pinmom_last(ik,ir) /
     >         (crmb * amu * kvhs(ik,ir))
           endif
         enddo
c
c         if (ir.eq.14) then
c           write(71,*)
c           write(71,*) 'PIN DATA:'
c           write(71,*)
c           write(71,'(A10,7A16)') 'ir   ik', 'pinatom  ',
c     >      'pinmol  ','pinena  ','pinion  ', 'pinrec  ',
c     >      'pinqi  ', 'pinqe  '
c           write(71,*)
c
c           do ik = 1,nks(ir)
c             write(71,'(2I5,7G16.8)') ir, ik,
c     >         pinatom(ik,ir), pinmol(ik,ir), pinena(ik,ir),
c     >         pinion(ik,ir), pinrec(ik,ir),
c     >         pinqi(ik,ir), pinqe(ik,ir)
c           enddo
c         endif
c
      enddo
      write(71,*)
c
c     RADIAL GRADIENTS
c
      if (s23_par_drflag.eq.1) then
        i_a = 9
        i_b = 10
        i_c = 11
        i_d = 12
        temp_7 = 1.0E19
        temp_8 = 1.0E23
        temp_9 = 1.0E6
c
        do ik = 1, nks(ir)
           write(71,'(5G16.8)') kss(ik,i_d),
     >    1.0 /  calcwidth2(ik,i_a) ** 2.0,
     >    1.0 /  calcwidth2(ik,i_b) ** 2.0,
     >    1.0 /  calcwidth2(ik,i_c) ** 2.0,
     >    1.0 /  calcwidth2(ik,i_d) ** 2.0
c
c     >         CalcDist(kss(ik,i_a),i_a,knbs,knds,9),
c     >         CalcDist(kss(ik,i_b),i_b,knbs,knds,9),
c     >         CalcDist(kss(ik,i_c),i_c,knbs,knds,9),
c     >         CalcDist(kss(ik,i_d),i_d,knbs,knds,9)
        enddo
        write(71,*)
c
        STOP
c
      endif
c
c     PRINT FILE FOR EXCELL
c
      if (s23_par_dsflag.eq.1) then
c
         call print_sol_excel
c
         STOP
c
      endif
c
c      Te_flag = 0
c      call Te_inv(irsep,irwall-1,midnks,Te_flag)
c      write(71,*) 'Te_flag', Te_flag
c
      if (pin_count.eq.0.and.s23_cfdfile.eq.0) then
c
         write(71,*) 'LOC 1'
c
         call TwoD_zero(pinion_last)
         call TwoD_zero(pinatom_last)
         call TwoD_zero(pinmol_last)
         call TwoD_zero(pinena_last)
         call TwoD_zero(pinqi_last)
         call TwoD_zero(pinqe_last)
         call TwoD_zero(pinrec_last)
         call TwoD_zero(pinmom_last)
         call TwoD_zero(d2ndr_last)
c
         call TwoD_zero(pinion_sm)
         call TwoD_zero(pinatom_sm)
         call TwoD_zero(pinmol_sm)
         call TwoD_zero(pinena_sm)
         call TwoD_zero(pinqi_sm)
         call TwoD_zero(pinqe_sm)
         call TwoD_zero(pinrec_sm)
         call TwoD_zero(pinmom_sm)
         call TwoD_zero(d2ndr_sm)
c
      elseif (pin_count.eq.1.and.s23_cfdfile.eq.0) then
c
         write(71,*) 'LOC 2'
c
         call TwoD_assign(pinion,   pinion_sm, 1.0)
         call TwoD_assign(pinrec,   pinrec_sm, 1.0)
         call TwoD_assign(pinatom,  pinatom_sm, 1.0)
         call TwoD_assign(pinmol,   pinmol_sm, 1.0)
         call TwoD_assign(pinena,   pinena_sm, 1.0)
         call TwoD_zero(pinqi_sm)
         call TwoD_zero(pinqe_sm)
         call TwoD_zero(pinmom_sm)
         call TwoD_zero(d2ndr_sm)
c
      elseif (pin_count.ge.1) then
c
         write(71,*) 'LOC 3'
c
         call TwoD_momloss(irlim1, irlim2)
c
c         do ir = irlim1, irlim2
c           do ik = 1,nks(ir)
c              pinmom(ik,ir) = 1.0 * Q_CX(ik,ir)
c              pinmom(ik,ir) = pinmp(ik,ir)
c           enddo
c         enddo
c
         call perp_diffusion(d2ndr_sm)
         call TwoD_assign(pinqe,  pinqe_sm,  1.0)
         call TwoD_assign(pinion, pinion_sm, 1.0)
         call TwoD_assign(pinatom, pinatom_sm, 1.0)
         call TwoD_assign(pinmol, pinmol_sm, 1.0)
         call TwoD_assign(pinena, pinena_sm, 1.0)
         call TwoD_assign(pinmom, pinmom_sm, 1.0)
         call TwoD_assign(pinrec, pinrec_sm, 1.0)
c
         call pinqi_1234(irlim1, irlim2)
         call TwoD_assign(pinqi,  pinqi_sm,  1.0)
c
        if (s23_par_qrec.eq.0) then
            call TwoD_zero(pinrec)
            call TwoD_zero(pinrec_sm)
            call TwoD_zero(pinrec_last)
        endif
c
      endif
c
      if (pin_count.eq.0.and.s23_cfdfile.eq.0) then
c
         write(71,*) 'LOC 4'
c
         do ir = irlim1, irlim2
           old_te(ir,2)   = kteds(idds(ir,2)) * 2.0
           old_ti(ir,2)   = ktids(idds(ir,2)) * 2.0
           old_ne(ir,2)   = knds(idds(ir,2))
           old_v(ir,2)    = kvds(idds(ir,2))
           old_mach(ir,2) = cmachno(ir,2)
c
           old_te(ir,1)   = kteds(idds(ir,1)) * 2.0
           old_ti(ir,1)   = ktids(idds(ir,1)) * 2.0
           old_ne(ir,1)   = knds(idds(ir,1))
           old_v(ir,1)    = kvds(idds(ir,1))
           old_mach(ir,1) = cmachno(ir,1)
         enddo
c
      else
c
         write(71,*) 'LOC 5'
c
         do ir = irlim1, irlim2
           kteds(idds(ir,1)) = old_te(ir,1)
           ktids(idds(ir,1)) = old_ti(ir,1)
           knds(idds(ir,1))  = old_ne(ir,1)
           kvds(idds(ir,1))  = old_v(ir,1)
           cmachno(ir,1)     = old_mach(ir,1)
c
           kteds(idds(ir,2)) = old_te(ir,2)
           ktids(idds(ir,2)) = old_ti(ir,2)
           knds(idds(ir,2))  = old_ne(ir,2)
           kvds(idds(ir,2))  = old_v(ir,2)
           cmachno(ir,2)     = old_mach(ir,2)
         enddo
c
         write(71,*)
         write(71,*) 'INNER TARGET'
         write(71,*)
c
         if (s23_2T.eq.0) then
           write(71,'(7A16)') 'Te_0    ',
     >     'ne_0   ','M_0     ', 'G_0/G_Lang  ',
     >     ' P_0/P_Lang   ', 'G_0      ', 'P_0 [MW/m2] '
           write(71,*)
           do ir = irlim1, irlim2
              call calctargfluxes(ir)
              write(71,'(7G16.8)') old_te(ir,1),
     >         old_ne(ir,1), old_mach(ir,1),
     >         gtarg(2) / gamma_Lang(ir,2),
     >        (ionptarg(2)+elecptarg(2))/power_Lang(ir,2),
     >         gamma_targ(ir,2), power_targ(ir,2)/1.0E6
           enddo
         elseif (s23_2T.eq.1) then
           write(71,'(9A12)') 'Te_0     ', 'Ti_0    ',
     >     'ne_0    ','M_0    ', 'G_0/G_Lg  ',
     >     ' Pi_0/Pi_Lg   ', ' Pe_0/Pe_Lg   ',
     >     ' Pi_0   ',  'Pe_0   '
           write(71,*)
           do ir = irlim1, irlim2
              call calctargfluxes(ir)
              write(71,'(9G12.6)') old_te(ir,1),
     >         old_ti(ir,1),
     >         old_ne(ir,1), old_mach(ir,1),
     >         gtarg(2) / gamma_Lang(ir,2),
     >         ionptarg(2) / poweri_Lang(ir,2),
     >         elecptarg(2) / powere_Lang(ir,2),
     >         ionptarg(2)/1.0E6,
     >         elecptarg(2)/1.0E6
           enddo
         endif
         write(71,*)
c
         write(71,*)
         write(71,*) 'OUTER TARGET'
         write(71,*)
c
         if (s23_2T.eq.0) then
           write(71,'(7A16)') 'Te_0    ',
     >     'ne_0   ','M_0     ', 'G_0/G_Lang  ',
     >     ' P_0/P_Lang   ', 'G_0      ', 'P_0 [MW/m2] '
           write(71,*)
           do ir = irlim1, irlim2
              call calctargfluxes(ir)
              write(71,'(7G16.8)') old_te(ir,2),
     >         old_ne(ir,2), old_mach(ir,2),
     >         gtarg(1) / gamma_Lang(ir,1),
     >        (ionptarg(1)+elecptarg(1))/power_Lang(ir,1),
     >         gamma_targ(ir,1), power_targ(ir,1)/ 1.0E6
           enddo
         elseif (s23_2T.eq.1) then
           write(71,'(9A12)') 'Te_0     ', 'Ti_0    ',
     >     'ne_0    ','M_0    ', 'G_0/G_Lg  ',
     >     ' Pi_0/Pi_Lg   ', ' Pe_0/Pe_Lg   ',
     >     ' Pi_0   ',  'Pe_0   '
           write(71,*)
           do ir = irlim1, irlim2
              call calctargfluxes(ir)
              write(71,'(9G12.6)') old_te(ir,2),
     >         old_ti(ir,2),
     >         old_ne(ir,2), old_mach(ir,2),
     >         gtarg(1) / gamma_Lang(ir,1),
     >         ionptarg(1) / poweri_Lang(ir,1),
     >         elecptarg(1) / powere_Lang(ir,1),
     >         ionptarg(1)/1.0E6,
     >         elecptarg(1)/1.0E6
           enddo
         endif
         write(71,*)
c
         write(71,*)
         write(71,*) 'INNER MIDPOINT'
         write(71,*)
         write(71,'(5A16)') 'Ti_mid  ',  'Te_mid  ',
     >      'ne_mid   ','P_e_in [MW/m2]  ', 'P_i_in [MW/m2] '                                                                                            '
         write(71,*)
c
         do ir = irlim1, irlim2
           npts = nks(ir)
           smax = ksmaxs(ir)
           do ik = 1, npts
             spts(ik)    = kss(ik,ir)
             sbnds(ik)   = ksb(ik,ir)
             if (ik.ne.1.and.
     >         (spts(ik-1).le.(smax/2.0).and.
     >          spts(ik)  .gt.(smax/2.0)))
     >          midnks = ik-1
                npts_mid = midnks
           enddo
           write(71,'(5G16.8)') ktibs(midnks+1,ir),
     >         ktebs(midnks+1,ir),
     >         knbs(midnks+1,ir),
     >         (power_flux(ir,2)+power_flux(ir,1))/2.0E6,
     > (power_flux(ir,1)+power_flux(ir,2))/2.0E6 * s23_par_ptipte
         enddo
         write(71,*)
c
         write(71,*)
         write(71,*) 'OUTER MIDPOINT'
         write(71,*)
         write(71,'(5A16)') 'Ti_mid  ',  'Te_mid  ',
     >      'ne_mid   ','P_e_in [MW/m2]  ', 'P_i_in [MW/m2] '
         write(71,*)
c
         do ir = irlim1, irlim2
           npts = nks(ir)
           smax = ksmaxs(ir)
           do ik = 1, npts
             spts(ik)    = kss(ik,ir)
             sbnds(ik)   = ksb(ik,ir)
             if (ik.ne.1.and.
     >         (spts(ik-1).le.(smax/2.0).and.
     >          spts(ik)  .gt.(smax/2.0)))
     >          midnks = ik-1
                npts_mid = midnks
           enddo
           write(71,'(5G16.8)') ktibs(midnks,ir),
     >         ktebs(midnks,ir),
     >         knbs(midnks,ir),
     >         (power_flux(ir,1)+power_flux(ir,2))/2.0E6,
     >(power_flux(ir,1)+power_flux(ir,2))/2.0E6 * s23_par_ptipte
         enddo
         write(71,*)
c
c         temp_1 = 0.1
c         temp_2 = 0.0
c
c         call TwoD_smooth(d2ndr_sm,   temp_1, temp_2)
c         call TwoD_smooth(pinion_sm,  temp_1, temp_2)
c         call TwoD_smooth(pinatom_sm, temp_1, temp_2)
c         call TwoD_smooth(pinmol_sm,  temp_1, temp_2)
c         call TwoD_smooth(pinena_sm,  temp_1, temp_2)
c         call TwoD_smooth(pinqi_sm,   temp_1, temp_2)
c         call TwoD_smooth(pinqe_sm,   temp_1, temp_2)
c         call TwoD_smooth(pinrec_sm,  temp_1, temp_2)
c         call TwoD_smooth(pinmom_sm,  temp_1, temp_2)
c
c         ir = 8
c         do ik = 1,nks(ir)
c           if (ir.eq.8) then
c             write(71,'(A10,2I5,3G12.6)') 'mom/mp', ir, ik,
c     >            pinmom_last(ik,ir), pinmom_sm(ik,ir),
c     >            pinmom_sm(ik,ir) * s23_relax +
c     >            (1.0 - s23_relax) * pinmom_last(ik,ir)
c           endif
c         enddo
c
      endif
c
      if (pin_count.le.1.and.s23_cfdfile.eq.0) then
c
         write(71,*) 'LOC 6'
c
        do ir = irlim1, irlim2
           fract1_new(ir)      = 0.0
           fract2_new(ir)      = 0.0
           fract3_new(ir)      = 0.0
           fract4_new(ir)      = 0.0
        enddo
c
        call TwoD_assign(d2ndr_sm,  d2ndr_last,  1.0)
        call TwoD_assign(pinion_sm, pinion_last, 1.0)
        call TwoD_assign(pinatom_sm, pinatom_last, 1.0)
        call TwoD_assign(pinmol_sm, pinmol_last, 1.0)
        call TwoD_assign(pinena_sm, pinena_last, 1.0)
        call TwoD_assign(pinqi_sm,  pinqi_last,  1.0)
        call TwoD_assign(pinqe_sm,  pinqe_last,  1.0)
        call TwoD_assign(pinmom_sm, pinmom_last, 1.0)
        call TwoD_assign(pinrec_sm, pinrec_last, 1.0)
c
      elseif (pin_count.ge.1) then
c
        write(71,*) 'LOC 7'
        zhi_1 = 0.0
        zhi_2 = 0.0
        zhi_3 = 0.0
        zhi_4 = 0.0
        zhi_n = 0.0
        power_ave = 0.0
c
        if (s23_2T.eq.0) then
          write(71,*)
          write(71,*) 'CONVERGENCE'
          write(71,*)
          write(71,'(A6,6A16)') 'ir  ', 'f(1)      ',
     >     'f(2)      ','f(3)      ', 'zh(1)  ',
     >     'zh(2)   ', 'zh(3)  '
          write(71,*)
        elseif (s23_2T.eq.1) then
          write(71,*)
          write(71,*) 'CONVERGENCE'
          write(71,*)
          write(71,'(A6,8A12)') 'ir  ', 'f(1)  ',
     >     'f(2)  ','f(3)  ', 'f(4)  ', 'zh(1)  ',
     >     'zh(2)   ', 'zh(3)  ', 'zh(4)  '
          write(71,*)
        endif
c
        do ir = irlim1, irlim2
c
           npts = nks(ir)
           smax = ksmaxs(ir)
           do ik = 1, npts
             spts(ik)  = kss(ik,ir)
             sbnds(ik) = ksb(ik,ir)
c
             ionpow_src(ik)  = pinqi_sm(ik,ir)
             elecpow_src(ik) = pinqe_sm(ik,ir) / s23_zhi0
             mom_src(ik)     = pinmom_sm(ik,ir)
             part_src(ik)    = pinion_sm(ik,ir)
     >                         - pinrec_sm(ik,ir)
           enddo
c
           sum1a = src_int_ab (part_src, 0.0D0,
     >              smax/2.0D0,spts,sbnds,npts,intopt)
           sum1b = src_int_ab (part_src, smax/2.0D0,
     >              smax,spts,sbnds,npts,intopt)
           sum1c = src_int_ab (part_src, 0.0D0,
     >              smax,spts,sbnds,npts,intopt)
c
           sum2a = src_int_ab (mom_src, 0.0D0,
     >              smax/2.0D0,spts,sbnds,npts,intopt)
           sum2b = src_int_ab (mom_src, smax/2.0D0,
     >              smax,spts,sbnds,npts,intopt)
           sum2c = src_int_ab (mom_src, 0.0D0,
     >              smax,spts,sbnds,npts,intopt)
c
           if (s23_2T.eq.0) then
             sum3a = src_int_ab (ionpow_src, 0.0D0,
     >              smax/2.0D0,spts,sbnds,npts,intopt)
     >          + src_int_ab (elecpow_src, 0.0D0,
     >              smax/2.0D0,spts,sbnds,npts,intopt)
             sum3b = src_int_ab (ionpow_src, smax/2.0D0,
     >              smax,spts,sbnds,npts,intopt)
     >          + src_int_ab (elecpow_src, smax/2.0D0,
     >              smax,spts,sbnds,npts,intopt)
             sum3c = src_int_ab (ionpow_src, 0.0D0,
     >              smax,spts,sbnds,npts,intopt)
     >          + src_int_ab (elecpow_src, 0.0D0,
     >              smax,spts,sbnds,npts,intopt)
           elseif (s23_2T.eq.1) then
             sum3a = src_int_ab (ionpow_src, 0.0D0,
     >              smax/2.0D0,spts,sbnds,npts,intopt)
             sum3b = src_int_ab (ionpow_src, smax/2.0D0,
     >              smax,spts,sbnds,npts,intopt)
             sum3c = src_int_ab (ionpow_src, 0.0D0,
     >              smax,spts,sbnds,npts,intopt)
             sum4a = src_int_ab (elecpow_src, 0.0D0,
     >              smax/2.0D0,spts,sbnds,npts,intopt)
             sum4b = src_int_ab (elecpow_src, smax/2.0D0,
     >              smax,spts,sbnds,npts,intopt)
             sum4c = src_int_ab (elecpow_src, 0.0D0,
     >              smax,spts,sbnds,npts,intopt)
           endif
c
           call calctargfluxes(ir)
c
           targ1a = abs(gtarg(1))
           targ1b = abs(gtarg(2))
           targ1c = abs(gtarg(3))
           targ2a = momtarg(1)
           targ2b = momtarg(2)
           targ2c = momtarg(3)
           if (s23_2T.eq.0) then
             targ3a = ionptarg(1) + elecptarg(1)
             targ3b = ionptarg(2) + elecptarg(2)
             targ3c = ionptarg(3) + elecptarg(3)
             core3a = power_flux(ir,1) / 2.0
             core3b = power_flux(ir,2) / 2.0
             core3c = core3a + core3b
           elseif (s23_2T.eq.1) then
             targ3a = ionptarg(1)
             targ3b = ionptarg(2)
             targ3c = ionptarg(3)
             targ4a = elecptarg(1)
             targ4b = elecptarg(2)
             targ4c = elecptarg(3)
             core3a = power_flux(ir,1) / 4.0
             core3b = power_flux(ir,2) / 4.0
             core3c = core3a + core3b
             core4a = power_flux(ir,1) / 4.0
             core4b = power_flux(ir,2) / 4.0
             core4c = core4a + core4b
           endif
c
           if (s23_symm.eq.1) then
c
c             fract1_pin = max(abs(sum1a / targ1a),
c     >                      abs(sum1b / targ1b))
c             fract2_pin = max(abs(sum2a / targ2a),
c     >                      abs(sum2b / targ2b))
c             fract3_pin = max(abs(sum3a / core3a),
c     >                      abs(sum3b / core3b))
c
c            INNER MASS, MOM, POW
c
             fract1_pin = abs(sum1b / targ1b)
             fract2_pin = abs(sum2b / targ2b)
             if (s23_2T.eq.0) then
               fract3_pin = abs(sum3b / core3b)
             elseif (s23_2T.eq.1) then
               fract3_pin = abs(sum3b / core3b)
               fract4_pin = abs(sum4b / core4b)
             endif
c
           elseif (s23_symm.eq.0) then
c
c            INNER: MOM; TOTAL: MASS, POW
c
             fract1_pin = abs(sum1c / targ1c)
             fract2_pin = abs(sum2b / targ2b)
             if (s23_2T.eq.0) then
               fract3_pin = abs(sum3c / core3c)
             elseif (s23_2T.eq.1) then
               fract3_pin = abs(sum3c / core3c)
               fract4_pin = abs(sum4c / core4c)
             endif
           endif
c
           fract1_new(ir) = s23_relax * fract1_pin +
     >         (1.0 - s23_relax) * fract1_old(ir)
           fract2_new(ir) = s23_relax * fract2_pin +
     >         (1.0 - s23_relax) * fract2_old(ir)
           fract3_new(ir) = s23_relax * fract3_pin +
     >         (1.0 - s23_relax) * fract3_old(ir)
           fract4_new(ir) = s23_relax * fract4_pin +
     >         (1.0 - s23_relax) * fract4_old(ir)
c
           do ik = 1, npts
             spts(ik)    = kss(ik,ir)
             sbnds(ik)   = ksb(ik,ir)
             if (ik.ne.1.and.
     >         (spts(ik-1).le.(smax/2.0).and.
     >          spts(ik)  .gt.(smax/2.0)))
     >          midnks = ik-1
                npts_mid = midnks
           enddo
c
c           temp_7 = knbs(midnks+1,ir) *
c     >              (2.0 * ktebs(midnks+1,ir)) * ech
c           temp_8 = abs(targ2b) * (1.0 + fract2_old(ir))
c
c           temp_5 = src_int_ab (part_src, smax/2.0D0,
c     >              smax,spts,sbnds,npts,intopt)
c           temp_6 = src_int_ab (mom_src, smax/2.0D0,
c     >              smax,spts,sbnds,npts,intopt)
c           temp_7 = temp_6 / temp_5
c
c           do ik = 1, npts
c             spts(ik)    = kss(ik,ir)
c             sbnds(ik)   = ksb(ik,ir)
c             part_src(ik) = pinmom_cx(ik,ir) /
c     >            (knbs(ik,ir) * kvhs(ik,ir))
c             mom_src(ik) = part_src(ik) * kvhs(ik,ir)
c           enddo
c
c           temp_5 = src_int_ab (part_src, smax/2.0D0,
c     >              smax,spts,sbnds,npts,intopt)
c           temp_6 = src_int_ab (mom_src, smax/2.0D0,
c     >              smax,spts,sbnds,npts,intopt)
c           temp_8 = temp_6 / temp_5
c
           zhi_sm(1,ir) = fract1_pin / fract1_old(ir)
           zhi_sm(2,ir) = fract2_pin / fract2_old(ir)
           zhi_sm(3,ir) = fract3_pin / fract3_old(ir)
           zhi_sm(4,ir) = fract4_pin / fract4_old(ir)
c
           if (s23_2T.eq.1.and.fract3_pin.eq.0.D0) then
              zhi_sm(3,ir) = 1.0D0
           endif
c
           if (s23_2T.eq.0) then
             zhi_temp(ir) = zhi_sm(3,ir)
           elseif (s23_2T.eq.1) then
             zhi_temp(ir) = zhi_sm(4,ir)
           endif
c
           delta1 = abs(fract1_new(ir) - fract1_old(ir))
           delta2 = abs(fract2_new(ir) - fract2_old(ir))
           delta3 = abs(fract3_new(ir) - fract3_old(ir))
           delta4 = abs(fract4_new(ir) - fract4_old(ir))
c
           temp_1 = s23_par_limrel1
           temp_2 = s23_par_limrel2
c
           temp_3 = s23_par_limrel3 * min(1.0,
     >               min(old_ti(ir,1),old_ti(ir,2))
     >             / max(ktibs(midnks,ir),ktibs(midnks+1,ir)))
           temp_4 = s23_par_limrel4 * min(1.0,
     >               min(old_te(ir,1),old_te(ir,2))
     >             / max(ktebs(midnks,ir),ktebs(midnks+1,ir)))
c
           relax1(ir) = s23_relax * min(1.0, temp_1/abs(delta1))
           relax2(ir) = s23_relax * min(1.0, temp_2/abs(delta2))
           relax3(ir) = s23_relax * min(1.0, temp_3/abs(delta3))
           relax4(ir) = s23_relax * min(1.0, temp_4/abs(delta4))
c
           if (min(ktebs(midnks,ir),ktebs(midnks+1,ir))
     >           .le.s23_par_Tulimit) then
             finished(ir) = 1
           endif
c
           if (finished(ir).eq.0) then
              rel_min = min(relax1(ir), relax2(ir),
     >                   relax3(ir), relax4(ir))
           else
              rel_min = 0.0D0
           endif
c
           relax1(ir) = rel_min
           relax2(ir) = rel_min
           relax3(ir) = rel_min
           relax4(ir) = rel_min
c
           if (ir.lt.irwall) then
             zhi_n = zhi_n + 1.0
             zhi_1 = zhi_1 + zhi_sm(1,ir)
             zhi_2 = zhi_2 + zhi_sm(2,ir)
             zhi_3 = zhi_3 + zhi_sm(3,ir)
             zhi_4 = zhi_4 + zhi_sm(4,ir)
             power_ave = power_ave +
     >                   (power_flux(ir,1) + power_flux(ir,2)) / 2.0
           endif
c
           if (s23_2T.eq.0) then
             write(71,'(I6,6G16.8)') ir,
     >          fract1_old(ir), fract2_old(ir), fract3_old(ir),
     >          zhi_sm(1,ir), zhi_sm(2,ir), zhi_sm(3,ir)
           elseif (s23_2T.eq.1) then
             write(71,'(I6,8G12.6)') ir,
     >          fract1_old(ir), fract2_old(ir), fract3_old(ir),
     >          fract4_old(ir), zhi_sm(1,ir),
     >          zhi_sm(2,ir), zhi_sm(3,ir), zhi_sm(4,ir)
           endif
        enddo
        write(71,*)
c
        zhi_1 = zhi_1 / zhi_n
        zhi_2 = zhi_2 / zhi_n
        zhi_3 = zhi_3 / zhi_n
        zhi_4 = zhi_4 / zhi_n
        power_ave = power_ave / zhi_n
c
        if (s23_2T.eq.0) then
          write(71,'(A10,3G16.8,A12,G16.8,A12,G16.8)') 'zhi_ave:',
     >      zhi_1, zhi_2 , zhi_3, 'P_SOL [MW]:',s23_pinqe/zhi_3,
     >      'Q(3)|| ave:', power_ave/1.0E6
          write(71,'(A10,3G16.8,A12,G16.8)') 'zhi_inf:'  ,
     >      1.0, 1.0, 1.0,  'n_u_limit   :', sol23_par_nulimit
        elseif (s23_2T.eq.1) then
          write(71,'(A10,4G12.6,A12,G12.6,A12,G12.6)') 'zhi_ave:',
     >      zhi_1, zhi_2 , zhi_3, zhi_4,
     >      'P_SOL [MW]:', s23_pinqe/zhi_4,
     >      'P_in_ave:', power_ave/1.0E6
          write(71,'(A10,4G12.6,A12,G12.6)') 'zhi_inf:'  ,
     >      1.0, 1.0, 1.0, 1.0,     'n_u_limit   :', sol23_par_nulimit
        endif
c
c
c       GLOBAL PARAMETRIZATION
c
c        temp_1 = 1.0
c        temp_2 = 1.0
c        temp_3 = 1.0
c
c        if (zhi_3.le.a_sm_r.and.
c     >      (zhi_1.ge.0.75.and.zhi_1.le.1.25).and.
c     >      (zhi_2.ge.0.75.and.zhi_2.le.sol23_par_nulimit)) then
c           temp_2 = 1.025
c        elseif(zhi_1.ge.3.0) then
c           temp_2 = 0.99
c        elseif (zhi_3.ge.1.0) then
c           do ir = irlim1, irlim2
c             relax1(ir) = 0.0
c             relax2(ir) = 0.0
c             relax3(ir) = 0.0
c           enddo
c        endif
c
c        if (zhi_1.ge.1.0.and.zhi_2.ge.1.0) then
c          temp_3 = 0.99
c        endif
c
c
        if (s23_2T.eq.0) then
          write(71,*)
          write(71,*) 'POWER:'
          write(71,*)
          write(71,'(A5,8A14)') 'ir  ',
     >     'Q(3)||_0','Q(3)||_N',
     >     'Q(3)||_0/P_0','Q(3)||_N/P_N',
     >     'P_0/P_BC_0','P_N/P_BC_N',
     >     'Te_0   ', 'Te_N   '
        elseif (s23_2T.eq.1) then
          write(71,*)
          write(71,*) 'POWER:'
          write(71,*)
          write(71,'(A5,8A14)') 'ir  ',
     >     'Q(4)||_0','Q(4)||_N',
     >     'Qei/Q(4)',
     >     'Te_0','Te_ctrl_0',
     >     'Te_N','Te_ctrl_N','f_rec_N'
        endif
c
        if (s23_2T.eq.0) then
          do ir = irsep, irwall - 1
             if (ir.eq.irsep) then
               zhi_sm(3,ir) = 0.75 * zhi_temp(ir) +
     >           0.25 * zhi_temp(ir+1)
             elseif (ir.eq.irwall-1) then
               zhi_sm(3,ir) = 0.75 * zhi_temp(ir) +
     >           0.25 * zhi_temp(ir-1)
             else
               zhi_sm(3,ir) = 0.5 * zhi_temp(ir) +
     >           0.25 * zhi_temp(ir-1) + 0.25 * zhi_temp(ir+1)
             endif
          enddo
        endif
c
        do ir = irlim1, irlim2
c
c         LOCAL PARAMETRIZATION
c
           npts = nks(ir)
           smax = ksmaxs(ir)
           do ik = 1, npts
             spts(ik)    = kss(ik,ir)
             sbnds(ik)   = ksb(ik,ir)
             if (ik.ne.1.and.
     >         (spts(ik-1).le.(smax/2.0).and.
     >          spts(ik)  .gt.(smax/2.0)))
     >          midnks = ik-1
                npts_mid = midnks
          enddo
c
          do ik = 1, nks(ir)
             part_src(ik) = pinrec_last(ik,ir)
          enddo
c
          temp_1 = src_int_ab (part_src, 0.0D0,
     >              smax/2.0D0,spts,sbnds,npts,intopt)
          temp_3 = src_int_ab (part_src, 0.0D0,
     >              smax,spts,sbnds,npts,intopt)
c
          if (temp_3.ne.0.0) then
            f_rec_0 = temp_1 / temp_3
          else
            f_rec_0 = 0.0
          endif
          f_rec_L = 1.0 - f_rec_0
c
          call calctargfluxes(ir)
          targ1a = abs(gtarg(1))
          targ1b = abs(gtarg(2))
          targ1c = abs(gtarg(3))
          targ2a = momtarg(1)
          targ2b = momtarg(2)
          targ2c = momtarg(3)
          if (s23_2T.eq.0) then
            targ3a = ionptarg(1) + elecptarg(1)
            targ3b = ionptarg(2) + elecptarg(2)
            targ3c = ionptarg(3) + elecptarg(3)
          elseif (s23_2T.eq.1) then
            targ3a = ionptarg(1)
            targ3b = ionptarg(2)
            targ3c = ionptarg(3)
            targ4a = elecptarg(1)
            targ4b = elecptarg(2)
            targ4c = elecptarg(3)
          endif
c
c          temp_1 = 1.0
c          temp_2 = 1.0
c          temp_3 = 1.0
c          temp_4 = 1.0
c
c          temp_7 = max(abs(zhi_sm(1,ir) - 1.0),
c     >                 abs(zhi_sm(2,ir) - 1.0))
c
c          if (s23_2T.eq.0) then
c             temp_4 = min(0.01,
c     >                    3.0 * abs(zhi_sm(3,ir) - a_sm_r)**2.0)
c             temp_5 = zhi_sm(3,ir)
c             temp_6 = s23_pinqe / zhi_3
c          elseif (s23_2T.eq.1) then
c             temp_4 = min(0.01,
c     >                    3.0 * abs(zhi_sm(4,ir) - a_sm_r)**2.0)
c             temp_5 = zhi_sm(4,ir)
c             temp_6 = s23_pinqe / zhi_4
c          endif
c
c          if (temp_5.le.a_sm_r.and.temp_6.le.sol23_par_nulimit) then
c             power_sector = 1
c          elseif (temp_5.ge.a_sm_r.and.temp_6.le.sol23_par_nulimit) then
c             power_sector = 2
c          elseif (temp_5.le.a_sm_r.and.temp_6.ge.sol23_par_nulimit) then
c             power_sector = 3
c          elseif (temp_5.ge.a_sm_r.and.temp_6.ge.sol23_par_nulimit) then
c             power_sector = 4
c          endif
c
c          if (power_sector.eq.1) then
c             temp_4 = 1.0 * temp_4
c          elseif (power_sector.eq.2) then
c             temp_4 = 1.0 * temp_4
c          else
c             temp_6 = min(0.01,
c     >                    0.25 * abs(zhi_sm(2,ir) - 1.0   )**2.0)
c             temp_4 = max(temp_4, temp_6)
c          endif
c
c          if (power_sector.eq.1) then
c
c            Q(3)|| down, T0 const
c
c             temp_3 = 1.0 - temp_4
c             power_flux(ir,1)  = power_flux(ir,1)  * temp_3
c             power_flux(ir,2)  = power_flux(ir,2)  * temp_3
c
c          elseif (power_sector.eq.2) then
c
c            Q(3)|| up,  T0 const
c
c             temp_3 = 1.0 + temp_4
c             power_flux(ir,1)  = power_flux(ir,1)  * temp_3
c             power_flux(ir,2)  = power_flux(ir,2)  * temp_3
c
c          elseif (power_sector.eq.3) then
c
c            Q(3)|| down, T0 down, with power weighting
c
c             if (zhi_sm(2,ir).le.1.0) then
c
c               relax2(ir) = relax2(ir) * 0.5
c
c               temp_5 = (power_flux(ir,1) /
c     >               (power_flux(ir,1) + power_ave)) ** 1.0
c
c               temp_3 = 1.0 - temp_4 * (1.0 - temp_5)
c               power_flux(ir,1)  = power_flux(ir,1)  * temp_3
c
c               temp_3 = 1.0 - temp_4 * temp_5
c               power_targ(ir,1)  = power_targ(ir,1)  * temp_3
c
c               temp_3 = temp_3 + temp_4 * temp_5
c               s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c
c               temp_5 = (power_flux(ir,2) /
c     >               (power_flux(ir,2) + power_ave)) ** 1.0
c
c               temp_3 = 1.0 - temp_4 * (1.0 - temp_5)
c               power_flux(ir,2)  = power_flux(ir,2)  * temp_3
c
c               temp_3 = 1.0 - temp_4 * temp_5
c               power_targ(ir,2)  = power_targ(ir,2)  * temp_3
c
c               temp_3 = temp_3 + temp_4 * temp_5
c               s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c
c             endif
c
c          elseif (power_sector.eq.4) then
c
c            Q(3)|| up,  T0 up, with power weighting
c
c             if (zhi_sm(2,ir).ge.1.0) then
c
c               relax2(ir) = relax2(ir) * 0.1
c
c               temp_5 = (power_flux(ir,1) /
c     >               (power_flux(ir,1) + power_ave)) ** 1.0
c
c               temp_3 = 1.0 + temp_4 * (1.0 - temp_5)
c               power_flux(ir,1)  = power_flux(ir,1)  * temp_3
c
c               temp_3 = 1.0 + temp_4 * temp_5
c               power_targ(ir,1)  = power_targ(ir,1)  * temp_3
c
c               temp_3 = temp_3 - temp_4 * temp_5
c               s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c
c               temp_5 = (power_flux(ir,2) /
c     >               (power_flux(ir,2) + power_ave)) ** 1.0
c
c               temp_3 = 1.0 + temp_4 * (1.0 - temp_5)
c               power_flux(ir,2)  = power_flux(ir,2)  * temp_3
c
c               temp_3 = 1.0 - temp_4 * temp_5
c               power_targ(ir,2)  = power_targ(ir,2)  * temp_3
c
c               temp_3 = temp_3 - temp_4 * temp_5
c               s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c
c             endif
c
c          endif
c
c          if (zhi_sm(3,ir).le.a_sm_r) then
c             temp_3 = 1.0 - temp_4
c             power_flux(ir,1)  = power_flux(ir,1)  * temp_3
c             s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c             power_flux(ir,2)  = power_flux(ir,2)  * temp_3
c             s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c          else
c             temp_3 = 1.0 + temp_4
c             power_flux(ir,1)  = power_flux(ir,1)  * temp_3
c             s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c             power_flux(ir,2)  = power_flux(ir,2)  * temp_3
c             s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c          endif
c
c          if (zhi_sm(2,ir).le.1.0) then
c             temp_3 = 1.0 + temp_6
c             power_flux(ir,1)  = power_flux(ir,1)  * 1.0
c             s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c             power_flux(ir,2)  = power_flux(ir,2)  * 1.0
c             s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c          else
c             temp_3 = 1.0 - temp_6
c             power_flux(ir,1)  = power_flux(ir,1)  * 1.0
c             s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c             power_flux(ir,2)  = power_flux(ir,2)  * 1.0
c             s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c          endif
c
c
c         G_0 ~ G_0_Langmuir
c
          temp_1 = gamma_targ(ir,1)
          temp_2 = gamma_targ(ir,2)
          temp_3 = s23_par_g0relax
c
          gamma_targ(ir,1)  = gamma_targ(ir,1) * (1.0 - temp_3) +
     >                 temp_3 * gamma_Lang(ir,1)
          gamma_targ(ir,2)  = gamma_targ(ir,2) * (1.0 - temp_3) +
     >                 temp_3 * gamma_Lang(ir,2)
c
          temp_1 = gamma_targ(ir,1) / temp_1
          temp_2 = gamma_targ(ir,2) / temp_2
c
          if (abs(1.0 - temp_1).ge.0.01) then
            part_flux(ir,1) = part_flux(ir,1) * temp_1
          endif
c
          if (abs(1.0 - temp_2).ge.0.01) then
            part_flux(ir,2) = part_flux(ir,2) * temp_2
          endif
c
c         P_0 ~ P_0_Langmuir
c
          temp_3 = s23_par_p0relax
          if (s23_2T.eq.0) then
            power_targ(ir,1)  = power_targ(ir,1) * (1.0 - temp_3) +
     >                 temp_3 * power_Lang(ir,1)
            power_targ(ir,2)  = power_targ(ir,2) * (1.0 - temp_3) +
     >                 temp_3 * power_Lang(ir,2)
          elseif (s23_2T.eq.1) then
            poweri_targ(ir,1)  = poweri_targ(ir,1) * (1.0 - temp_3) +
     >                 temp_3 * poweri_Lang(ir,1)
            poweri_targ(ir,2)  = poweri_targ(ir,2) * (1.0 - temp_3) +
     >                 temp_3 * poweri_Lang(ir,2)
            powere_targ(ir,1)  = powere_targ(ir,1) * (1.0 - temp_3) +
     >                 temp_3 * powere_Lang(ir,1)
            powere_targ(ir,2)  = powere_targ(ir,2) * (1.0 - temp_3) +
     >                 temp_3 * powere_Lang(ir,2)
            power_targ(ir,1) = poweri_targ(ir,1) + powere_targ(ir,1)
            power_targ(ir,2) = poweri_targ(ir,2) + powere_targ(ir,2)
          endif
c
c         if (ir.le.irwall-1) then
c
c           P_0 > 1.0 * P_0_Langmuir
c
c            if (targ3a.le.1.0*power_targ(ir,1)) then
c               temp_3 = targ3a /
c     >             (0.9 * targ3a + 0.1 * power_targ(ir,1))
c               temp_3 = 0.99
c               power_flux(ir,1)  = power_flux(ir,1)  * 1.0
c               s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c            endif
c
c            if (targ3b.le.1.0*power_targ(ir,2)) then
c               temp_3 = targ3b /
c     >             (0.9 * targ3b + 0.1 * power_targ(ir,2))
c               temp_3 = 0.99
c               power_flux(ir,2)  = power_flux(ir,2)  * 1.0
c               s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c            endif
c
c           P_0 < 1.5 * P_0_Langmuir
c
c            if (targ3a.ge.1.5*power_targ(ir,1)) then
c               temp_3 = 1.01
c               power_flux(ir,1)  = power_flux(ir,1)  * 1.0
c               s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c            endif
c
c            if (targ3b.ge.1.5*power_targ(ir,2)) then
c               temp_3 = 1.01
c               power_flux(ir,2)  = power_flux(ir,2)  * 1.0
c               s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c            endif
c
c           Q(3)|| > 1.25 * targ_power
c
c            if (s23_powlim(ir,1).le.1.25) then
c               temp_3 = 1.0 + 0.01
c               power_flux(ir,1)  = power_flux(ir,1)  * temp_3
c               s23_powlim(ir,1)  = s23_powlim(ir,1)  * temp_3
c            endif
c
c            if (s23_powlim(ir,2).le.1.25) then
c               temp_3 = 1.0 + 0.01
c               power_flux(ir,2)  = power_flux(ir,2)  * temp_3
c               s23_powlim(ir,2)  = s23_powlim(ir,2)  * temp_3
c            endif
c
c          endif
c
          if (s23_par_dvmrel.eq.1) then
            temp_1 = 0.0D0
            temp_2 = 0.1D0
            temp_3 = 0.3D0
          elseif (s23_par_dvmrel.eq.2) then
            temp_1 = 0.0D0
            temp_2 = 0.03D0
            temp_3 = 0.1D0
          elseif (s23_par_dvmrel.eq.3) then
            temp_1 = 0.0D0
            temp_2 = 0.01D0
            temp_3 = 0.03D0
          elseif (s23_par_dvmrel.eq.4) then
            temp_1 = 0.0D0
            temp_3 = min(0.1, 1.0 * min(old_te(ir,1),old_te(ir,2))
     >                 / max(ktebs(midnks,ir),ktebs(midnks+1,ir)))
            temp_2 = 0.1D0
          elseif (s23_par_dvmrel.eq.5) then
            temp_1 = 0.0D0
            temp_3 = min(0.1, 1.0 * min(old_te(ir,1),old_te(ir,2))
     >                 / max(ktebs(midnks,ir),ktebs(midnks+1,ir)))
            temp_2 = temp_3
          elseif (s23_par_dvmrel.eq.6) then
            temp_1 = 0.0D0
            temp_3 = min(0.03, 0.3 * min(old_te(ir,1),old_te(ir,2))
     >                 / max(ktebs(midnks,ir),ktebs(midnks+1,ir)))
            temp_2 = temp_3
          elseif (s23_par_dvmrel.eq.7) then
            temp_1 = 0.0D0
            temp_3 = min(0.03, 0.3 * min(old_te(ir,1),old_te(ir,2))
     >                 / max(ktebs(midnks,ir),ktebs(midnks+1,ir)))
            temp_2 = 0.1D0
          else
             temp_1 = 0.0D0
             temp_2 = 0.0D0
             temp_3 = 0.0D0
          endif
c
          if (finished(ir).eq.1.or.pin_count.le.1) then
             temp_1 = 0.0D0
             temp_2 = 0.0D0
             temp_3 = 0.0D0
          endif
c
          power_flux(ir,1) = max(0.0D0, power_flux(ir,1) +
     >                               temp_3 * del_power_flux(ir,1))
c
          power_flux(ir,2) = max(0.0D0, power_flux(ir,2) +
     >                                temp_3 * del_power_flux(ir,2))
c
          mom_flux(ir,1) = mom_flux(ir,1) + temp_2 *
     >                       del_mom_flux(ir,1)
c
          mom_flux(ir,2) = mom_flux(ir,2) + temp_2 *
     >                       del_mom_flux(ir,2)
c
          part_flux(ir,1) = part_flux(ir,1) + temp_1 *
     >                       del_part_flux(ir,1)
c
          part_flux(ir,2) = part_flux(ir,2) + temp_1 *
     >                       del_part_flux(ir,2)
c
c           d/dr (Pin||)0+N > 0 and Te,Lang > Temin eV
c
          if (s23_par_dpdrflag.eq.1) then
           if (ir.gt.irsep.and.ir.le.irwall-1.and.
     >       (power_flux(ir-1,1)+power_flux(ir-1,2)).le.
     >       (power_flux(ir,1)+power_flux(ir,2))) then
            if (Te_Lang(ir,1).ge.s23_par_dpdrtemin) then
              power_flux(ir,1) = max(0.0D0, power_flux(ir,1)
     >         - (1.0 + s23_par_dpdrstep) * temp_3
     >           * abs(del_power_flux(ir,1)))
            endif
            if (Te_Lang(ir,2).ge.s23_par_dpdrtemin) then
              power_flux(ir,2) = max(0.0D0, power_flux(ir,2)
     >         - (1.0 + s23_par_dpdrstep) * temp_3
     >           * abs(del_power_flux(ir,2)))
            endif
           endif
          endif
c
c           UPSTREAM DENSITY CONTROL VIA TE_CTRL
c
          if (s23_par_nuflag.eq.1) then
            if (knbs(midnks+1,irsep).le.sol23_par_nulimit) then
              if (zhi_sm(2,irsep).ge.1.0) then
                temp_7 = (knbs(midnks+1,irsep) -
     >                     sol23_par_nulimit) / sol23_par_nulimit
              else
                temp_7 = (knbs(midnks+1,irsep) -
     >                     sol23_par_nulimit) / sol23_par_nulimit
              endif
            else
              if (zhi_sm(2,irsep).ge.1.0) then
                temp_7 = (knbs(midnks+1,irsep) -
     >                     sol23_par_nulimit) / sol23_par_nulimit
                relax2(ir) = relax2(ir) * 0.0
              else
                temp_7 = 0.0D0
              endif
            endif
          else
            temp_7 = 0.0
          endif
c
          temp_1 = temp_7 * f_rec_0
          temp_2 = temp_7 * f_rec_L
c
c          write(71,*) 'Te_CTRL', Te_ctrl(ir,1), Te_ctrl(ir,2)
c          write(71,*) 'Te_LANG', Te_Lang(ir,1), Te_Lang(ir,2)
c
          Te_ctrl(ir,1) = max(Te_Lang(ir,1),
     >                   Te_ctrl(ir,1) * (1.0 + temp_1))
          Te_ctrl(ir,2) = max(Te_Lang(ir,2),
     >                   Te_ctrl(ir,2) * (1.0 + temp_2))
c
c          write(71,*) 'Te_CTRL', Te_ctrl(ir,1), Te_ctrl(ir,2)
c          write(71,*) 'Te_LANG', Te_Lang(ir,1), Te_Lang(ir,2)
c
          power_flux(ir,3) = power_flux(ir,1) + power_flux(ir,2)
          mom_flux(ir,3)   = mom_flux(ir,1)   + mom_flux(ir,2)
          part_flux(ir,3)  = part_flux(ir,1)  + part_flux(ir,2)
c
          if (s23_2T.eq.0) then
            write(71,'(I6,8G14.6)') ir,
     >        power_flux(ir,1)/1.0E6,  power_flux(ir,2)/1.0E6,
     >        power_flux(ir,1)/targ3a, power_flux(ir,2)/targ3b,
     >        targ3a/power_targ(ir,1), targ3b/power_targ(ir,2),
     >        old_te(ir,2), old_te(ir,1)
          elseif (s23_2T.eq.1) then
            write(71,'(I6,8G14.6)') ir,
     >        power_flux(ir,1)/2.0E6,  power_flux(ir,2)/2.0E6,
     >        power_ei(ir) * 0.5D0 /
     >          (power_flux(ir,1) + power_flux(ir,2)),
     >        old_te(ir,2), Te_ctrl(ir,1),
     >        old_te(ir,1), Te_ctrl(ir,2), f_rec_L
          endif
        enddo
        write(71,*)
c
        write(71,*)
        write(71,*) 'MOMENTUM:'
        write(71,*)
        write(71,'(A10,5A16)') 'ir  ',
     >      ' MOM_IN_0   ', ' MOM_IN_N ',
     >      ' DEL_MOM_0  ', ' DEL_MOM_N'
        write(71,*)
c
        do ir = irlim1, irlim2
          write(71,'(I6,4G14.6)') ir,
     >        mom_flux(ir,1) / mom_Lang(ir,3),
     >        mom_flux(ir,2) / mom_Lang(ir,3),
     >        del_mom_flux(ir,1) / mom_Lang(ir,3),
     >        del_mom_flux(ir,2) / mom_Lang(ir,3)
        enddo
        write(71,*)
c
c        write(71,*)
c        write(71,*) 'MASS:'
c        write(71,*)
c        write(71,'(A10,5A16)') 'ir  ',
c     >      ' MASS_IN_0   ', ' MASS_IN_N ',
c     >      ' DEL_MASS_0  ', ' DEL_MASS_N'
c        write(71,*)
c
c        do ir = irlim1, irlim2
c          write(71,'(I6,4G14.6)') ir,
c     >        part_flux(ir,1) / gamma_Lang(ir,3),
c     >        part_flux(ir,2) / gamma_Lang(ir,3),
c     >        del_part_flux(ir,1) / gamma_Lang(ir,3),
c     >        del_part_flux(ir,2) / gamma_Lang(ir,3)
c        enddo
c        write(71,*)
c
c      ir = 8
c         npts = nks(ir)
c         smax = ksmaxs(ir)
c         if (ir.eq.8) then
c           write(71,*)
c           write(71,*) 'MOMENTUM:'
c           write(71,*)
c           write(71,'(A10,5A16)') 'ir   ik',
c     >      ' s [m]   ',  'p + nmv^2 [Pa]',
c     >      'Int[Qcx]  ','Int[Qrec]  ', 'sum  '
c           write(71,*)
c         endif
c
c         do ik = 1, npts
c             spts(ik)    = kss(ik,ir)
c             sbnds(ik)   = ksb(ik,ir)
c             part_src(ik) = - pinmom_cx (ik,ir) /
c     >             zhi_sm(2,ir)
c             mom_src(ik)  = - pinmom_rec(ik,ir) /
c     >             zhi_sm(2,ir)
c             if (ik.ne.1.and.
c     >         (spts(ik-1).le.(smax/2.0).and.
c     >          spts(ik)  .gt.(smax/2.0)))
c     >          midnks = ik-1
c                npts_mid = midnks
c         enddo
c
c         do ik = midnks+1,npts
c           temp_1 = knbs(ik,ir) * (ktebs(ik,ir) + ktibs(ik,ir))
c     >              * ech + crmb * amu * knbs(ik,ir) *
c     >              kvhs(ik,ir) ** 2.0
c           temp_2 = src_int_ab (part_src, smax/2.0D0,
c     >              spts(ik),spts,sbnds,npts,intopt)
c           temp_3 = src_int_ab (mom_src, smax/2.0D0,
c     >              spts(ik),spts,sbnds,npts,intopt)
c           temp_4 = temp_1 + temp_2 + temp_3
c           if (ir.eq.8) then
c             write(71,'(2I5,5G16.8)') ir, ik, spts(ik),
c     >         temp_1, temp_2, temp_3, temp_4
c           endif
c         enddo
c
        write(71,*)
        write(71,*) 'RELAX:'
        write(71,*)
        write(71,'(A6,4A14)') 'ir  ', 'relax1  ',
     >   'relax2  ','relax3  ', 'relax4  '
        write(71,*)
c
        do ir = irlim1, irlim2
           write(71,'(I6,4G14.8)') ir,
     >        relax1(ir), relax2(ir), relax3(ir) , relax4(ir)
        enddo
        write(71,*)
c
c        do ir = irlim1, irlim2
c           delta1 = abs(fract1_new(ir) - fract1_old(ir))
c           write(71,'(A10,I5,4G16.8)') 'dts1',ir,
c     >        fract1_old(ir), fract1_new(ir), delta1, relax1(ir)
c
c        enddo
c        write(71,*)
c
c
c        write(71,*) 'ION:'
c        call TwoD_compare(pinion_sm, pinion_last)
c        write(71,*) 'REC:'
c        call TwoD_compare(pinrec_sm, pinrec_last)
c        write(71,*) 'MOM:'
c        call TwoD_compare(pinmom_sm, pinmom_last)
c        write(71,*) 'QE:'
c        call TwoD_compare(pinqe_sm,  pinqe_last)
c
        call TwoD_relax(d2ndr_sm,  d2ndr_last,  relax1)
        call TwoD_relax(pinion_sm, pinion_last, relax1)
        call TwoD_relax(pinatom_sm, pinatom_last, relax1)
        call TwoD_relax(pinmol_sm, pinmol_last, relax1)
        call TwoD_relax(pinena_sm, pinena_last, relax1)
        call TwoD_relax(pinrec_sm, pinrec_last, relax1)
        call TwoD_relax(pinmom_sm, pinmom_last, relax2)
        call TwoD_relax(pinqi_sm,  pinqi_last,  relax4)
        call TwoD_relax(pinqe_sm,  pinqe_last,  relax3)
c
      endif
c
c     Loop over rings
c
 799  do ir = irlim1,irlim2
c
c        Do not try solving on the grid boundary rings
c
         if (ir.eq.irwall.or.ir.eq.irtrap) goto 1000
c
c        Mark type of ring
c
         if (ir.gt.irwall) then
c
c           PP ring
c
            ring_type = 1
c
         else
c
c           Main Sol ring
c
            ring_type = 1
c
         endif
c
c        Set Ring values
c
         npts = nks(ir)
         smax = ksmaxs(ir)
c
c        Copy grid
c
         sbnds(0) = ksb(0,ir)
c
         do ik = 1,nks(ir)
c
            spts(ik) = kss(ik,ir)
            sbnds(ik) = ksb(ik,ir)
c
            if (ik.ne.1.and.
     >         (spts(ik-1).le.(smax/2.0).and.
     >          spts(ik)  .gt.(smax/2.0)))
     >          midnks = ik-1
                npts_mid = midnks
c
c            write (6,*) 'KSS:',ksb(ik,ir),kss(ik,ir),ksb(ik-1,ir)
c            write (6,*) 'SPT:',sbnds(ik),spts(ik),sbnds(ik-1)
c
c
         end do
c
c         write(6,*) 'smax:',ir,smax,ksmaxs(ir),ksmaxs2(ir),
c     >              sbnds(npts)
c
c
c        Load boundary condtions
c
c        KVDS should contain the sound speed at the target if this
c        is the first iteration.
c
c         te_0 = kteds(idds(ir,2))
c         ti_0 = ktids(idds(ir,2))
c         ne_0 = knds(idds(ir,2))
c         isat_0 = knds(idds(ir,2)) * abs(kvds(idds(ir,2)))
c         mach_0 = 1.0
c
c         te_smax = kteds(idds(ir,1))
c         ti_smax = ktids(idds(ir,1))
c         ne_smax = knds(idds(ir,1))
c         isat_smax = knds(idds(ir,1)) * abs(kvds(idds(ir,1)))
c         mach_smax = 1.0
c
c        Perform setup operations if seeplasma is OFF
c
c
c        Setup the grid - load ring geometry information
c
c         write (6,*) 'Testing Integration:'
c
         do ik = 1,nks(ir)
c
            part_src(ik) = ik
c
         end do
c
         srcsum1 = src_int_ab(part_src,0.0d0,smax,
     >                           spts,sbnds,npts,0)
         srcsum2 = src_int_ab(part_src,0.0d0,smax,
     >                           spts,sbnds,npts,1)
c
c         write (6,*) 'Whole:',srcsum1,srcsum2
c
         srcsum1 = src_int_ab(part_src,0.25d0*smax,0.75d0*smax,
     >                           spts,sbnds,npts,0)
         srcsum2 = src_int_ab(part_src,0.25d0*smax,0.75d0*smax,
     >                           spts,sbnds,npts,1)
c         write (6,*) 'Half :',srcsum1,srcsum2
c
c
         srcsum1 = src_int_ab(part_src,0.0d0,0.5d0*smax,
     >                           spts,sbnds,npts,0)
         srcsum2 = src_int_ab(part_src,0.0d0,0.5d0*smax,
     >                           spts,sbnds,npts,1)
c         write (6,*) 'First Half :',srcsum1,srcsum2
c
         srcsum1 = src_int_ab(part_src,0.5d0*smax,smax,
     >                           spts,sbnds,npts,0)
         srcsum2 = src_int_ab(part_src,0.5d0*smax,smax,
     >                           spts,sbnds,npts,1)
c         write (6,*) 'Second Half :',srcsum1,srcsum2
c
c
         do ik = 1,nks(ir)
c
            srcsum1 = src_int_ab(part_src,sbnds(ik-1),sbnds(ik),
     >                           spts,sbnds,npts,0)
c
            srcsum2 = src_int_ab(part_src,sbnds(ik-1),sbnds(ik),
     >                           spts,sbnds,npts,1)
c
c            write (6,*) 'Cell :',ik,srcsum1,srcsum2
c
         end do
c
         call dzero(part_src,maxpts)
c
         if (pin_count.ge.1) then
c
c            write(71,*) 'LOC 8'
c
c
c           Copy in the initial sources - these will be adjusted
c           so that all quantities are conserved.
c
            do ik = 1, nks(ir)
               part_src(ik) = pinion_last(ik,ir)
            enddo
            srcsum1 = src_int(part_src)
c
            do ik = 1, nks(ir)
               part_src(ik) = pinrec_last(ik,ir)
            enddo
            srcsum2 = src_int(part_src)
c
            f_recomb(ir) = srcsum2 / srcsum1
c
            write(71,'(A16,G16.8)') 'RECOMB/IONIZ', srcsum2 / srcsum1
c            write(72,'(A16,G16.8)') 'RECOMB/IONIZ', srcsum2 / srcsum1
c
            do ik = 1,nks(ir)
c
              part_src(ik) = pinion_last(ik,ir) - pinrec_last(ik,ir)
c
c               ionpow_src(ik)  = min(0.0, pinqi_last(ik,ir))
c               ionpow_src(ik)  = pinqi_last(ik,ir) *
c     >            min(1.0, pinrec_last(ik,ir)/pinion_last(ik,ir))
c
               ionpow_src(ik)  = pinqi_last(ik,ir)
               elecpow_src(ik) = pinqe_last(ik,ir) / s23_zhi0
c
               mom_src(ik)     = pinmom_last(ik,ir)
               perp_src(ik)    = 1.0D0 / calcwidth2(ik,ir) ** 2.0D0
               magn_field(ik)  = abs(bts(ik,ir))
c
               h_neut(ik)      = pinatom(ik,ir)
               h2_neut(ik)     = pinmol(ik,ir)
               th_neut(ik)     = pinena(ik,ir) * 0.66667
c
            enddo
c
            call calctargfluxes(ir)
c
c            if (pin_count.eq.1) then
c                do ik = 1, nks(ir)
c                   part_src(ik) = pinion_last(ik,ir)
c                enddo
c                srcsum1 = src_int(part_src)
c
c               do ik = 1, nks(ir)
c                  part_src(ik) = part_src(ik) * dabs(gtarg(3)/srcsum1)
c                  pinion_last(ik,ir) = part_src(ik)
c               enddo
c            endif
c
c           Ion power Source
c
c            call adjust_src(ionpow_src,ionptarg,opt)
c
c           Electron power Source
c
c            call adjust_src(elecpow_src,elecptarg,opt)
c
c           On first iteration calculate some multiple of
c           the target power fluxes on each ring and use
c           this to fix the power flow on iteration - for testing ONLY
c
         elseif (pin_count.eq.0.and.seedplasma.ne.0) then
c
            write(71,*) 'LOC 9'
c
c
c           Seed plasma parameters - to be moved to INPUT file
c
            S23_izlen    = sol23_izlen  * smax
            S23_izlam    = sol23_izlam  * smax
            S23_izoffset = sol23_izoffset * smax
            S23_momlen   = sol23_momlen * smax
c
c           This assumes sonic flow at the target
c
c            call calctargfluxes(ir)
c
            write(71,*) 'DEFINE POWER_FLUX'
c
c            power_flux(ir,1) = 2.0D0 * (ionptarg(1) + elecptarg(1))
c            power_flux(ir,2) = 2.0D0 * (ionptarg(2) + elecptarg(2))
c
c            power_flux(ir,1) = s23_maxpow_0 * power_flux(ir,1)
c            power_flux(ir,2) = s23_maxpow_N * power_flux(ir,2)
c
c            gamma_targ(ir,1) = abs(gtarg(1))
c            gamma_targ(ir,2) = abs(gtarg(2))
c
c            power_targ(ir,1) =  (ionptarg(1) + elecptarg(1))
c            power_targ(ir,2) =  (ionptarg(2) + elecptarg(2))
c
            gamma_targ(ir,1) = gamma_Lang(ir,1)
            gamma_targ(ir,2) = gamma_Lang(ir,2)
c
            if (s23_2T.eq.0) then
              power_targ(ir,1)  = power_Lang(ir,1)
              power_targ(ir,2)  = power_Lang(ir,2)
            elseif (s23_2T.eq.1) then
              poweri_targ(ir,1) = poweri_Lang(ir,1)
              poweri_targ(ir,2) = poweri_Lang(ir,2)
              powere_targ(ir,1) = powere_Lang(ir,1)
              powere_targ(ir,2) = powere_Lang(ir,2)
              power_targ(ir,1) = poweri_targ(ir,1) + powere_targ(ir,1)
              power_targ(ir,2) = poweri_targ(ir,2) + powere_targ(ir,2)
            endif
c
            power_flux(ir,1) = 2.0D0 * (poweri_Lang(ir,1) +
     >                         powere_Lang(ir,1))
            power_flux(ir,2) = 2.0D0 * (poweri_Lang(ir,2) +
     >                         powere_Lang(ir,2))
c
            power_flux(ir,1) = 1.0 * power_flux(ir,1)
            power_flux(ir,2) = 1.0 * power_flux(ir,2)
c
c            write(71,*) 'gamma_targ', gamma_targ(ir,1),
c     >                  gamma_targ(ir,2)
c            write(71,*) 'power_targ', power_targ(ir,1),
c     >                  power_targ(ir,2)
c            write(71,*) 'power_flux', power_flux(ir,1),
c     >                  power_flux(ir,2)
c
         endif
c
c        SYMMETRIC SOLUTION
c
         if (s23_symm.eq.1) then
c
c          Invoke Solver for OUTER half ring, s = 0
c
           solve_half = 1
c
           write(72,*) 'RING: ', ir, ' [0,N/2] '
           write(71,*) 'RING: ', ir, ' [0,N/2] '
c           write(0,*)  'RING: ', ir, ' [0,N/2] '
c
           pow_in_0     = power_flux(ir,1) / 2.0
           pow_in_N     = power_flux(ir,1) / 2.0
c
           if (s23_2T.eq.0) then
             pow_tg_0    = power_targ(ir,1)
             pow_tg_N    = power_targ(ir,1)
           elseif (s23_2T.eq.1) then
             pow_tg_i0    = poweri_targ(ir,1)
             pow_tg_iN    = poweri_targ(ir,1)
             pow_tg_e0    = powere_targ(ir,1)
             pow_tg_eN    = powere_targ(ir,1)
           endif
c
           mom_in_0     = mom_flux(ir,1)
           mom_in_N     = mom_flux(ir,1)
c
           mom_tg_0     = mom_targ(ir,1)
           mom_tg_N     = mom_targ(ir,1)
c
           isat_0       = gamma_targ(ir,1)
           isat_smax    = gamma_targ(ir,1)
c
           ir_cfd       = ir
           targ_cfd     = 1
c
           te_0   = old_te(ir,2)
           ti_0   = old_ti(ir,2)
           ne_0   = old_ne(ir,2)
           mach_0 = old_mach(ir,2)
c
           te_smax   = old_te(ir,2)
           ti_smax   = old_ti(ir,2)
           ne_smax   = old_ne(ir,2)
           mach_smax = old_mach(ir,2)
c
           do ik = 1,npts
              ne(ik) = oldknbs(ik,ir)
              te(ik) = oldktebs(ik,ir)
              ti(ik) = oldktibs(ik,ir)
              vb(ik) = oldkvhs(ik,ir)
c             write(71,'(I5,3G16.8)') ik, ne(ik), vb(ik), te(ik)
           end do
c
c           call calctargfluxes(ir)
c           write(71,*) 'BEFORE CFD_OSM'
c           write(71,*) 'g_Lang, gtarg', gamma_Lang(ir,1),
c     >                  gtarg(1), 1.0 - gtarg(1)/gamma_Lang(ir,1)
c           write(71,*) 'pow_Lang, pow_targ', power_Lang(ir,1),
c     >                  ionptarg(1) + elecptarg(1),  1.0 -
c     >         (ionptarg(1) + elecptarg(1)) / power_Lang(ir,1)
c
c           s23_powlim(ir,1) = 0.9 * s23_powlim(ir,1) +
c     >         0.1 * power_flux(ir,1) / (2.0 * power_targ(ir,1))
c
c           s23_powlim(ir,1) = power_flux(ir,1)
c     >                        / (2.0 * power_targ(ir,1))
c           s23_maxpow_0 = s23_powlim(ir,1)
c
c           write(71,*) 'a_maxpow, f_rad', s23_maxpow_0,
c     >        1.0 - 1.0 / s23_maxpow_0
c           write(71,*) 'P_in,targ,Lang', power_flux(ir,1),
c     >        power_targ(ir,1), power_Lang(ir,1)
c
c
           call cfd_osm(ir,err_fl_a)
c
           if (err_fl_a.eq.0) then
c
             do ik = 1,midnks
                knbs(ik,ir) = ne(ik)
                ktebs(ik,ir) = te(ik)
                ktibs(ik,ir) = ti(ik)
                kvhs(ik,ir)  = vb(ik)
             end do
c
             old_te(ir,2) = te_0
             old_ti(ir,2) = ti_0
             old_ne(ir,2)  = ne_0
             old_v(ir,2)  = - abs(isat_0)/ne_0
             old_mach(ir,2) = abs(mach_0)
c
             kteds(idds(ir,2)) = te_0
             ktids(idds(ir,2)) = ti_0
             knds(idds(ir,2))  = ne_0
             kvds(idds(ir,2))  = - abs(isat_0)/ne_0
             cmachno(ir,2)     = abs(mach_0)
c
c             write(71,'(A10,5G16.8)') 'TARG0',
c     >         kteds(idds(ir,2)),
c     >         ktids(idds(ir,2)), knds(idds(ir,2)),
c     >         kvds(idds(ir,2)), cmachno(ir,2)
c             write(71,*) 'nTM_0_SOL',ne_0, ti_0, mach_0
c
c             call calctargfluxes(ir)
c             write(71,*) 'AFTER CFD_OSM'
c             write(71,*) 'g_Lang, gtarg', gamma_Lang(ir,1),
c     >                  gtarg(1), 1.0 - gtarg(1)/gamma_Lang(ir,1)
c             write(71,*) 'pow_targ, pow_targ', power_Lang(ir,1),
c     >         ionptarg(1) + elecptarg(1), 1.0 -
c     >         (ionptarg(1) + elecptarg(1)) / power_Lang(ir,1)
c
           endif
c
c          Invoke Solver for INNER half ring, s = s_max
c
           solve_half = 2
c
           write(72,*) 'RING: ', ir, ' [N/2,N] '
           write(71,*) 'RING: ', ir, ' [N/2,N] '
c           write(0,*)  'RING: ', ir, ' [N/2,N] '
c
           pow_in_0     = power_flux(ir,2) / 2.0
           pow_in_N     = power_flux(ir,2) / 2.0
c
           if (s23_2T.eq.0) then
             pow_tg_0    = power_targ(ir,2)
             pow_tg_N    = power_targ(ir,2)
           elseif (s23_2T.eq.1) then
             pow_tg_i0    = poweri_targ(ir,2)
             pow_tg_iN    = poweri_targ(ir,2)
             pow_tg_e0    = powere_targ(ir,2)
             pow_tg_eN    = powere_targ(ir,2)
           endif
c
           mom_in_0     = mom_flux(ir,2)
           mom_in_N     = mom_flux(ir,2)
c
           mom_tg_0     = mom_targ(ir,2)
           mom_tg_N     = mom_targ(ir,2)
c
           isat_0       = gamma_targ(ir,2)
           isat_smax    = gamma_targ(ir,2)
c
           ir_cfd       = ir
           targ_cfd     = 2
c
           te_0   = old_te(ir,1)
           ti_0   = old_ti(ir,1)
           ne_0   = old_ne(ir,1)
           mach_0 = old_mach(ir,1)
c
           te_smax   = old_te(ir,1)
           ti_smax   = old_ti(ir,1)
           ne_smax   = old_ne(ir,1)
           mach_smax = old_mach(ir,1)
c
           do ik = 1,npts
              ne(ik) = oldknbs(ik,ir)
              te(ik) = oldktebs(ik,ir)
              ti(ik) = oldktibs(ik,ir)
              vb(ik) = oldkvhs(ik,ir)
           end do
c
c           call calctargfluxes(ir)
c           write(71,*) 'BEFORE CFD_OSM'
c           write(71,*) 'g_targ, gtarg', gamma_targ(ir,2),
c     >                  gtarg(2), 1.0 - gtarg(2)/gamma_targ(ir,2)
c           write(71,*) 'pow_targ, pow_targ', power_targ(ir,2),
c     >                  ionptarg(2) + elecptarg(2),  1.0 -
c     >         (ionptarg(2) + elecptarg(2)) / power_targ(ir,2)
c
c           s23_powlim(ir,2) = 0.9 * s23_powlim(ir,2) +
c     >         0.1 * power_flux(ir,2) / (2.0 * power_targ(ir,2))
c
c           s23_powlim(ir,2) = power_flux(ir,2)
c     >                        / (2.0 * power_targ(ir,2))
c           s23_maxpow_N = s23_powlim(ir,2)
c
c           write(71,*) 'a_maxpow, f_rad', s23_maxpow_N,
c     >        1.0 - 1.0 / s23_maxpow_N
c           write(71,*) 'P_in, P_targ', power_flux(ir,2),
c     >        power_targ(ir,2)
c
           call cfd_osm(ir,err_fl_b)
c
           if (err_fl_b.eq.0) then
c
             do ik = midnks+1,nks(ir)
                knbs(ik,ir) = ne(ik)
                ktebs(ik,ir) = te(ik)
                ktibs(ik,ir) = ti(ik)
                kvhs(ik,ir)  = vb(ik)
             end do
c
             old_te(ir,1)   = te_smax
             old_ti(ir,1)   = ti_smax
             old_ne(ir,1)   = ne_smax
             old_v(ir,1)    = abs(isat_smax)/ne_smax
             old_mach(ir,1) = abs(mach_smax)
c
             kteds(idds(ir,1)) = te_smax
             ktids(idds(ir,1)) = ti_smax
             knds(idds(ir,1))  = ne_smax
             kvds(idds(ir,1))  = abs(isat_smax)/ne_smax
             cmachno(ir,1)     = abs(mach_smax)
c
c             write(71,'(A10,5G16.8)') 'TARGN',
c     >         kteds(idds(ir,1)),
c     >         ktids(idds(ir,1)), knds(idds(ir,1)),
c     >         kvds(idds(ir,1)), cmachno(ir,1)
c             write(71,*) 'nTM_s_SOL',ne_smax,ti_smax,mach_smax
c
c             call calctargfluxes(ir)
c             write(71,*) 'AFTER CFD_OSM'
c             write(71,*) 'g_targ, gtarg', gamma_targ(ir,2),
c     >                  gtarg(2), 1.0 - gtarg(2)/gamma_targ(ir,2)
c             write(71,*) 'pow_targ, pow_targ', power_targ(ir,2),
c     >                  ionptarg(2) + elecptarg(2), 1.0 -
c     >         (ionptarg(2) + elecptarg(2)) / power_targ(ir,2)
c
           endif
c
c          REPLACE BACKGROUND PLASMA WITH NEW VALUES
c
           if (err_fl_a.eq.0.and.err_fl_b.eq.0) then
             do ik = 1,nks(ir)
                oldknbs(ik,ir)  = knbs(ik,ir)
                oldktebs(ik,ir) = ktebs(ik,ir)
                oldktibs(ik,ir) = ktibs(ik,ir)
                oldkvhs(ik,ir)  = kvhs(ik,ir)
             enddo
c
c            INNER ONLY
c
             fract1_old(ir) = abs(s23_fract1_N)
             fract2_old(ir) = abs(s23_fract2_N)
             fract3_old(ir) = abs(s23_fract3_N)
             fract4_old(ir) = abs(s23_fract4_N)
           else
             err_dummy = 1
           endif
c
c          NON-SYMMETRIC SOLUTION
c
         elseif (s23_symm.eq.0) then
c
c          Invoke Solver for the entire ring
c
           solve_half = 0
c
           write(72,*) 'RING: ', ir, ' [0,N] '
           write(71,*) 'RING: ', ir, ' [0,N] '
c           write(0,*)  'RING: ', ir, ' [0,N] '
c
           pow_in_0     = power_flux(ir,1) / 2.0
           pow_in_N     = power_flux(ir,2) / 2.0
c
           if (s23_2T.eq.0) then
             pow_tg_0    = power_targ(ir,1)
             pow_tg_N    = power_targ(ir,2)
           elseif (s23_2T.eq.1) then
             pow_tg_i0    = poweri_targ(ir,1)
             pow_tg_iN    = poweri_targ(ir,2)
             pow_tg_e0    = powere_targ(ir,1)
             pow_tg_eN    = powere_targ(ir,2)
           endif
c
           part_in_0    = part_flux(ir,1)
           part_in_N    = part_flux(ir,2)
c
           mom_in_0     = mom_flux(ir,1)
           mom_in_N     = mom_flux(ir,2)
c
           mom_tg_0     = mom_targ(ir,1)
           mom_tg_N     = mom_targ(ir,2)
c
           isat_0       = gamma_targ(ir,1)
           isat_smax    = gamma_targ(ir,2)
c
           ir_cfd     = ir
           targ_cfd   = 1
c
           te_0   = old_te(ir,2)
           ti_0   = old_ti(ir,2)
           ne_0   = old_ne(ir,2)
           mach_0 = old_mach(ir,2)
c
           te_smax   = old_te(ir,1)
           ti_smax   = old_ti(ir,1)
           ne_smax   = old_ne(ir,1)
           mach_smax = old_mach(ir,1)
c
           do ik = 1, nks(ir)
              ne(ik) = oldknbs(ik,ir)
              te(ik) = oldktebs(ik,ir)
              ti(ik) = oldktibs(ik,ir)
              vb(ik) = oldkvhs(ik,ir)
c             write(71,'(I5,3G16.8)') ik, ne(ik), vb(ik), te(ik)
           end do
c
           err_fl_b = 0
           call cfd_osm(ir,err_fl_a)
c
c           write(71,*) 'Del, Pin0', del_power_flux(ir,1), pow_in_0
c           write(71,*) 'Del, PinN', del_power_flux(ir,2), pow_in_N
c
           if (err_fl_a.eq.0) then
c
             do ik = 1, nks(ir)
                knbs(ik,ir) = ne(ik)
                ktebs(ik,ir) = te(ik)
                ktibs(ik,ir) = ti(ik)
                kvhs(ik,ir)  = vb(ik)
             end do
c
c             call calctargfluxes(ir)
c             write(71,*) 'BEFORE ', gtarg(1),ionptarg(1),
c     >                    elecptarg(1)
c
             old_te(ir,2) = te_0
             old_ti(ir,2) = ti_0
             old_ne(ir,2)  = ne_0
             old_v(ir,2)  = - abs(isat_0)/ne_0
             old_mach(ir,2) = abs(mach_0)
c
             kteds(idds(ir,2)) = te_0
             ktids(idds(ir,2)) = ti_0
             knds(idds(ir,2))  = ne_0
             kvds(idds(ir,2))  = - abs(isat_0)/ne_0
             cmachno(ir,2)     = abs(mach_0)
c
             old_te(ir,1)   = te_smax
             old_ti(ir,1)   = ti_smax
             old_ne(ir,1)   = ne_smax
             old_v(ir,1)    = abs(isat_smax)/ne_smax
             old_mach(ir,1) = abs(mach_smax)
c
             kteds(idds(ir,1)) = te_smax
             ktids(idds(ir,1)) = ti_smax
             knds(idds(ir,1))  = ne_smax
             kvds(idds(ir,1))  = abs(isat_smax)/ne_smax
             cmachno(ir,1)     = abs(mach_smax)
c
c             write(71,'(A10,5G16.8)') 'TARG0',
c     >         kteds(idds(ir,2)),
c     >         ktids(idds(ir,2)), knds(idds(ir,2)),
c     >         kvds(idds(ir,2)), cmachno(ir,2)
c             write(71,'(A10,5G16.8)') 'TARGN',
c     >         kteds(idds(ir,1)),
c     >         ktids(idds(ir,1)), knds(idds(ir,1)),
c     >         kvds(idds(ir,1)), cmachno(ir,1)
c
c             call calctargfluxes(ir)
c             write(71,*) 'AFTER  ', gtarg(1),ionptarg(1),
c     >                             elecptarg(1)
c
           endif
c
c
c          REPLACE BACKGROUND PLASMA WITH NEW VALUES
c
           if (err_fl_a.eq.0.and.err_fl_b.eq.0) then
c
             do ik = 1,nks(ir)
                oldknbs(ik,ir)  = knbs(ik,ir)
                oldktebs(ik,ir) = ktebs(ik,ir)
                oldktibs(ik,ir) = ktibs(ik,ir)
                oldkvhs(ik,ir)  = kvhs(ik,ir)
             enddo
c
c            INNER: MOM; TOTAL: MASS, POW
c
             fract1_old(ir) = abs(s23_fract1)
             fract2_old(ir) = abs(s23_fract2)
             fract3_old(ir) = abs(s23_fract3)
             fract4_old(ir) = abs(s23_fract4)
           else
             err_dummy = 1
           endif
c
         endif
c
C
C       CALCULATE ELECTRIC FIELD
C
C
C       IN THE FOLLOWING EQUATIONS THE FACTOR E CANCELS WITH THE
C       SAME FACTOR USED IN CONVERTING T IN EV TO KT.
C
c       If OFIELD is turned ON then set the electric field to
c       zero for this case. Note: the electric field is
c       initialized to zero in the plasma.d3a module.
c
c       For OFIELD ... 0=off   1=on
c
        if (ofield.eq.0) then
c
c
          if (kss2(1,ir).eq.0.0) then
            DS1 = KSS2(2,IR) - KSS2(1,IR)
            DP1 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
            DT1 = (KTEBS(2,IR)-KTEBS(1,IR))
            NB1 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
C
            KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
c
            KEDS(idds(ir,2))= kes(1,ir)
          else
            DS2 = KSS2(2,IR) - KSS2(1,IR)
            DP2 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
            DT2 = (KTEBS(2,IR)-KTEBS(1,IR))
            NB2 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
c
            DS1 = KSS2(1,IR)
            DP1 = KNBS(1,IR)*KTEBS(1,IR)-
     >              KNDS(idds(ir,2))*KTEDS(idds(ir,2))
            DT1 = KTEBS(1,IR)-KTEDS(idds(ir,2))
            NB1 = 0.5*(KNBS(1,IR)+KNDS(idds(ir,2)))
c
            KES(1,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
c
            KEDS(idds(ir,2)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
          endif
C
c
          if (kss2(nks(ir),ir).eq.ksmaxs2(ir)) then
            DS1 = KSS2(NKS(IR),IR) - KSS2(NKS(IR)-1,IR)
            DP1 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >         -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
            DT1 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
            NB1 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
C
            KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
c
            KEDS(idds(ir,1))= kes(nks(ir),ir)
          else
c
            DS2 = KSS2(NKS(IR),IR) - KSS2(NKS(IR)-1,IR)
            DP2 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >         -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
            DT2 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
            NB2 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
c
            DS1 = ksmaxs2(ir) - KSS2(nks(ir),IR)
            DP1 = KNDS(idds(ir,1))*KTEDS(idds(ir,1))-
     >             KNBS(nks(ir),IR)*KTEBS(nks(ir),IR)
            DT1 = KTEDS(idds(ir,1))-KTEBS(nks(ir),IR)
            NB1 = 0.5*(KNBS(nks(ir),IR)+KNDS(idds(ir,1)))
c
            KES(nks(ir),IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
c
            KEDS(idds(ir,1)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
          endif
C
C        WRITE(6,*) 'KES:IR:',KES(1,IR),KES(NKS(IR),IR)
C
          DO 500 IK = 2,NKS(IR)-1

            DS1 = KSS2(IK,IR) - KSS2(IK-1,IR)
            DP1 = KNBS(IK,IR)*KTEBS(IK,IR)-KNBS(IK-1,IR)*KTEBS(IK-1,IR)
            DT1 = (KTEBS(IK,IR)-KTEBS(IK-1,IR))
            NB1 = 0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
            DS2 = KSS2(IK+1,IR) - KSS2(IK,IR)
            DP2 = KNBS(IK+1,IR)*KTEBS(IK+1,IR)-KNBS(IK,IR)*KTEBS(IK,IR)
            DT2 = (KTEBS(IK+1,IR)-KTEBS(IK,IR))
            NB2 = 0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))
            KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
C
C            WRITE(6,*) 'KES:',IK,IR,KES(IK,IR)
C
 500      CONTINUE
c
c         End of test on OFIELD
c
        endif
c
c
c
c
c     Loop over more rings
c
 1000 end do
c
c
      end_time = za02as(1)
      elapsed_time = end_time - start_time
c      write(0,*)  'TIME SPENT IN CFD (min)', elapsed_time/60.0
      write(71,*) 'TIME SPENT IN CFD (min)', elapsed_time/60.0
      write(72,*) 'TIME SPENT IN CFD (min)', elapsed_time/60.0
c
c      if (pin_count.eq.nitersol) call wrtcfdbg(irlim1,irlim2)
c
      call wrtcfdbg(irlim1,irlim2)
c
c      call Te_inv(irsep,irwall-1,midnks,Te_flag)
c      write(71,*) 'Te_flag', Te_flag
c      write(71,*) 'PINTEMP:',pin_temp
c
c      if (Te_flag.eq.0.and.pin_temp.ge.5) then
c           pin_temp = 0
c
c           write(71,*) 'f_Q3 up 0.5%'
c           do ir = irlim1, irlim2
c              fract_cfd(ir,1,3) = 1.005 * fract_cfd(ir,1,3)
c              fract_cfd(ir,2,3) = 1.005 * fract_cfd(ir,2,3)
c              write(71,*) ir, fract_cfd(ir,1,3), fract_cfd(ir,2,3)
c           enddo
c
c           do ir = irlim1, irlim2
c              s23_powlim(ir) = s23_powlim(ir) * 1.01
c              write(71,*) 'POWER limit up', s23_powlim(ir)
c           enddo
c
c          do ir = irlim1, irlim2
c             power_flux(ir,1) = power_flux(ir,1) * 1.0139595
c             power_flux(ir,2) = power_flux(ir,2) * 1.0139595
c          enddo
c
c          sol23_sm_s = sol23_sm_s + 0.25
c          write(71,*) 'f_Z up'
c
c      elseif (Te_flag.eq.1.and.pin_temp.ge.5) then
c           pin_temp = 0
c
c          sol23_sm_s = max(0.0, sol23_sm_s - 0.25)
c          write(71,*) 'f_Z down'
c      endif
c
c      write(71,*) 'PINQE limit:', sol23_sm_r
c      write(71,*) 'power_flux:', power_flux(irlim1,1), sol23_maxpow
c
      if (err_dummy.eq.1) then
        write(71,*) 'ERROR IN CFD SOLN'
        STOP
      endif
c
      return
      end
c
c
      subroutine echosol23
      use mod_params
      use mod_comtor
      use mod_sol23_input
      implicit none
c
c     include 'params'
c     include 'comtor'
c     include 'sol23_input'
c
c     ECHOSOL23:
c
c     This routine prints out the input option used for SOL 23 to
c     the output data file.
c


      CALL PRC ('  SOL OPTION      23  : CFD ring by ring plasma solver'
     >)
c
      call prb
      call prc ('              SOL23 SUB-OPTION LIST:')
      call prb
c
c     Integration option
c
      if (sol23_intopt.eq.0) then
c
        call prc ('     INTEGRATION OPT  0: Integrations are performed')
        call prc ('                         over step-function sources')
c
      elseif (sol23_intopt.eq.1) then
c
        call prc ('     INTEGRATION OPT  1: Integrations are performed')
        call prc ('                         by linearly interpolating')
        call prc ('                         between source points.')
c
      endif
c
c     Seed Plasma Option
c
      if (sol23_seed.eq.1) then
c
        call prc ('     SEED PLASMA OPT  1: Exponential'//
     >                                 'Decay Ionization Source')
        call prc ('                         Note: *smax for ring')
        call prr ('                         Ionization Length = ',
     >                                        sol23_izlen)
        call prr ('                         Ionization Decay  = ',
     >                                        sol23_izlam)
        call prr ('                         Momentum Length   = ',
     >                                        sol23_momlen)
c
      elseif (sol23_seed.eq.2) then
c
        call prc ('     SEED PLASMA OPT  2: s5-Gaussian '//
     >                                     'Ionization Source')
        call prc ('                         Note: *smax for ring')
        call prr ('                         Ionization Length = ',
     >                                        sol23_izlen)
        call prr ('                         Ionization Decay  = ',
     >                                        sol23_izlam)
        call prr ('                         Ionization Offset = ',
     >                                        sol23_izoffset)
        call prr ('                         Momentum Length   = ',
     >                                        sol23_momlen)
      endif
c
c
      call prb
c
      return
      end

c
c
c
      subroutine adjust_src(quant,outflow,opt)
      use mod_params
      implicit none
c
c     include 'params'
c
      integer opt
      real*8 quant(maxpts),outflow(3)
c
c     ADJUST_SRC:
c
c     This subroutine adjusts the given source so that it's integral
c     will be equal to the given outflow.
c
c
      integer  ik
      real*8     srcsum
      real*8     extra_src
      real*8     src_int
      external src_int
c
c     The different values of srcsum represent two different methods
c     of integrating the two sources.
c
      srcsum = src_int(quant)
c
      extra_src =  outflow(3) - srcsum
c
c      write (6,*) 'Adjusting Source:',extra_src,outflow(3),
c     >                  srcsum
c
      call srcupdate(quant,extra_src,opt)
c
c     For testing - double check result
c
c      srcsum =  src_int(quant)
c      extra_src =  outflow(3) - srcsum
c
c      write (6,*) 'Checking Source :',extra_src,outflow(3),
c     >                  srcsum
c
      return
      end
c
c
c
      subroutine calctargfluxes(ir)
      use mod_params
      use mod_comtor
      use mod_sol23_input
      use mod_cgeom
      use mod_sol23_com
      implicit none
c
c     include 'params'
c     include 'comtor'
c     include 'sol23_input'
c     include 'cgeom'
c
c     include 'sol23_com'
c
      integer ir
c
c
c
c     This subroutine calculates the target partcle and
c     power fluxes based on the target data.
c
c     Local Variables
c
      real v0
      real mach0
      real gae,gaio,gaii
c
c     Calculate target fluxes (both particles and heat)
c
c     OUTER
c
      mach0 = cmachno(ir,2)
      gae = gae_const
      gaio = gai_const +  0.5* mach0**2
     >         * (1.0 + kteds(idds(ir,2))/ktids(idds(ir,2)))
c
      v0 = sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb
     >          * emi) * mach0
      gtarg(1)= knds(idds(ir,2)) * v0
      momtarg(1) = knds(idds(ir,2)) *
     >   (kteds(idds(ir,2)) + ktids(idds(ir,2))) * ech *
     >   (1.0 + mach0 * mach0)
      elecptarg(1)= gae*kteds(idds(ir,2))*ech*gtarg(1)
      ionptarg(1) = gaio*ktids(idds(ir,2))*ech*gtarg(1)
c
c     INNER
c
      mach0 = cmachno(ir,1)
      gae = gae_const
      gaii = gai_const +  0.5* mach0**2
     >         * (1.0 + kteds(idds(ir,1))/ktids(idds(ir,1)))
c
      v0 = sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))/crmb
     >             * emi) * mach0
      gtarg(2)= knds(idds(ir,1)) * v0
      momtarg(2) = knds(idds(ir,1)) *
     >   (kteds(idds(ir,1)) + ktids(idds(ir,1))) * ech *
     >   (1.0 + mach0 * mach0)
      elecptarg(2)=gae*kteds(idds(ir,1))*ech*gtarg(2)
      ionptarg(2)=gaii*ktids(idds(ir,1))*ech*gtarg(2)
c
c     Total
c
      gtarg(3) = gtarg(1) + gtarg(2)
      momtarg(3) = momtarg(2) - momtarg(1)
      elecptarg(3)= elecptarg(1) + elecptarg(2)
      ionptarg(3) = ionptarg(1) + ionptarg(2)
c
c
      return
      end
c
c
c
      real*8 function src_int(quant)
      use mod_params
      use mod_sol23_com
      implicit none
c
c     include 'params'
c     include 'sol23_com'
c
      real*8 quant(maxpts)
c
c
c     SRC_INT:
c
c     This routine calculates the total integral of the passed in
c     quantity on teh particular ring. In order
c     to evaluate the net total excess or deficit of flux
c     onto the ring. These values are then used to balance the
c     sources on the ring
c
c     Local variables
c
      real*8 srcsum
      integer ik
c
      srcsum  = 0.0
c
c     Calculate the integral in different ways depending on options
c
      if (intopt.eq.0) then
c
         do ik = 1,npts
c
            srcsum = srcsum + quant(ik)*(sbnds(ik)-sbnds(ik-1))
c
         end do
c
      elseif (intopt.eq.1) then
c
         do ik = 1,npts
c
            if (ik.eq.1) then

               srcsum = srcsum + quant(ik)*(spts(ik)-sbnds(ik-1))

            elseif (ik.eq.npts) then

               srcsum = srcsum
     >                   + (quant(ik)+quant(ik-1))/2.0
     >                     * (spts(ik)-spts(ik-1))

               srcsum = srcsum + quant(ik)*(sbnds(ik)-spts(ik))

            else

               srcsum = srcsum
     >                   + (quant(ik)+quant(ik-1))/2.0
     >                     * (spts(ik)-spts(ik-1))

            endif
c
         end do
c
      endif
c
      src_int = srcsum
c
      return
      end
c
c
c
      subroutine srcupdate(quant,extra_src,opt)
      use mod_params
      use mod_sol23_com
      implicit none
c
c     include 'params'
c     include 'sol23_com'
c
      real*8    quant(maxpts)
      real*8    extra_src
      integer opt
c
c
c     SRCUPDATE:
c
c     This routine updates the source terms by including the extra_src
c     in the source. The value of extra_src may be positive or negative.
c
c     There will eventually be various ways of allocating this source
c     to the ring. For now - all that is supported is a uniform disposition.
c
c     Local variables
c
      integer ik
      real*8    extra_src_density
c
      extra_src_density = extra_src/smax
c
      do ik = 1,npts
c
         if (opt.eq.1) then
            quant(ik) = quant(ik) + extra_src_density
         endif
c
      end do
c
c
      return
      end
c
c
c
      integer function src_index(s,spts,sbnds,npts)
      use mod_params
      implicit none
c     include 'params'
c
      integer npts
      real*8 s
      real*8 spts(npts)
      real*8 sbnds(0:npts)
c
c     SRC_INDEX:
c
c     This routine returns an ik index into the S-distances array
c     passed in as an argument.
c
c
      integer iklast,ik
      data iklast /0/
c
c     If S is too small
c
      if (s.le.sbnds(0)) then
c
         src_index = 1
c
c     IF S is too large
c
      elseif (s.ge.sbnds(npts)) then

         src_index = npts
c
c     Check if still in last cell called
c
      elseif (iklast.ge.1.and.iklast.le.npts) then
c
c        Still in last cell
c
         if (s.ge.sbnds(iklast-1).and.s.le.sbnds(iklast)) then
            src_index = iklast
c
c        Start scanning - down
c
         elseif (s.lt.sbnds(iklast-1)) then
c
            do ik = iklast-1,1,-1
               if (s.gt.sbnds(ik-1)) then
                  src_index = ik
                  goto 10
               endif
            end do
 10         continue
c
c        Start scanning - up
c
         elseif (s.gt.sbnds(iklast)) then
c
            do ik = iklast+1,npts
               if (s.lt.sbnds(ik)) then
                  src_index = ik
                  goto 20
               endif
            end do
 20         continue
c
         endif
c
c     Start brand new scan
c
      else
c
         do ik = 1,npts
            if (s.lt.sbnds(ik)) then
               src_index = ik
               goto 30
            endif
         end do
 30      continue
      endif
c
c     Set iklast value and return
c
      iklast = src_index
c
      return
      end
c
c
c
      real*8 function src_val(quant,s,spts,sbnds,npts,intopt)
      use mod_params
      implicit none
c     include 'params'
c
      integer npts,intopt
      real*8   s,quant(npts)
      real*8   sbnds(0:npts)
      real*8   spts(npts)
c
c     SRC_VAL:
c
c     This routine returns a value of the discrete function quant at
c     the value S. Depending on the value of opt - the return quantity
c     will be either the value of the function in the cell OR a value
c     interpolated between adjacent cell centres.
c
c     The value of opt should be consistent with the usage in the
c     integration routine - thus option 0 uses the value in cell
c     result and option 1 returns an interpolated value.
c
c
      integer ik
      integer src_index
      external src_index
c
c     If S is too small
c
      if (s.lt.sbnds(0)) then
c
         src_val = 0.0
c
c     IF S is too large
c
      elseif (s.gt.sbnds(npts)) then
c
         src_val = 0.0
c
c     S within proper range - find corresponding quant value
c
      else
c
         ik = src_index(s,spts,sbnds,npts)
c
c        Cell value
c
         if (intopt.eq.0) then
c
            src_val = quant(ik)
c
c        Interpolated value
c
         elseif (intopt.eq.1) then
c
c           Greater than cell centre - interpolate up
c
            if (s.gt.spts(ik)) then
c
               if (ik.eq.npts) then
c
                  src_val = quant(ik)
c
               else
c
                  src_val = quant(ik) +
     >                 (s-spts(ik))/(spts(ik+1)-spts(ik)) *
     >                 (quant(ik+1)-quant(ik))
c
               endif
c
c           Less than cell centre - interpolate down
c
            elseif (s.le.spts(ik)) then
c
               if (ik.eq.1) then
c
                  src_val = quant(ik)
c
               else
c
                  src_val = quant(ik-1) +
     >                 (s-spts(ik-1))/(spts(ik)-spts(ik-1)) *
     >                 (quant(ik)-quant(ik-1))
c
               endif
c
            end if
c
         end if
c
      endif
c
      return
      end
c
c
c
      real*8 function src_val_ik(quant,s,ik,spts,sbnds,npts,intopt)
      use mod_params
      implicit none
c     include 'params'
c
      integer npts,intopt
      real*8   s,quant(npts)
      real*8   spts(npts)
      real*8   sbnds(0:npts)
      integer ik
c
c     SRC_VAL_IK:
c
c     This routine returns a value of the discrete function quant at
c     the value S. Depending on the value of opt - the return quantity
c     will be either the value of the function in the cell OR a value
c     interpolated between adjacent cell centres.
c
c     The value of opt should be consistent with the usage in the
c     integration routine - thus option 0 uses the value in cell
c     result and option 1 returns an interpolated value.
c
c
      integer src_index,iktmp,ikold
      external src_index
c
      iktmp = ik
c
c     If S is too small
c
      if (s.lt.sbnds(0)) then
c
         src_val_ik = 0.0
c
c     IF S is too large
c
      elseif (s.gt.sbnds(npts)) then
c
         src_val_ik = 0.0
c
c     S within proper range - find corresponding quant value
c
      else
c
c        Verify correct cell - otherwise perform lookup
c
         if (s.lt.sbnds(iktmp-1).or.s.gt.sbnds(iktmp)) then
c
            ikold = iktmp

c
            iktmp = src_index(s,spts,sbnds,npts)
c
c            write (6,*) 'Changing ik:',iktmp,ikold,
c     >                 sbnds(ikold-1),s,sbnds(ikold)
c
c
         endif
c
c
c        Cell value
c
         if (intopt.eq.0) then
c
            src_val_ik = quant(iktmp)
c
c        Interpolated value
c
         elseif (intopt.eq.1) then
c
c           Greater than cell centre - interpolate up
c
            if (s.gt.spts(iktmp)) then
c
               if (iktmp.eq.npts) then
c
                  src_val_ik = quant(iktmp)
c
               else
c
                  src_val_ik = quant(iktmp) +
     >                 (s-spts(iktmp))/
     >                 (spts(iktmp+1)-spts(iktmp)) *
     >                 (quant(iktmp+1)-quant(iktmp))
c
               endif
c
c           Less than cell centre - interpolate down
c
            elseif (s.le.spts(iktmp)) then
c
               if (iktmp.eq.1) then
c
                  src_val_ik = quant(iktmp)
c
               else
c
                  src_val_ik = quant(iktmp-1) +
     >                 (s-spts(iktmp-1))/
     >                 (spts(iktmp)-spts(iktmp-1)) *
     >                 (quant(iktmp)-quant(iktmp-1))
c
               endif
c
            end if
c
         end if
c
      endif
c
c      write (6,'(a,i4,8g12.5)') 'src:',iktmp,
c     >               src_val_ik,s,sbnds(iktmp-1),
c     >                  spts(iktmp),sbnds(iktmp),
c     >                   quant(iktmp+1),quant(iktmp),
c     >                   quant(iktmp-1)
c
c
      return
      end
c
c
c
      real*8 function src_int_ab(quant,s1,s2,spts,sbnds,npts,intopt)
      use mod_params
      implicit none
c
c     include 'params'
c
      integer npts,intopt
      real*8 quant(npts)
      real*8 sbnds(0:npts)
      real*8 spts(npts)
      real*8 s1,s2
c
c
c     SRC_INT_AB:
c
c     This routine calculates the total integral of the passed in
c     quantity on teh particular ring. In order
c     to evaluate the net total excess or deficit of flux
c     onto the ring. These values are then used to balance the
c     sources on the ring
c
c     Local variables
c
      integer ik,ik1,ik2
      real*8 srcsum,sval1,sval2
c
      real*8 src_val_ik
      external src_val_ik
c
      integer src_index
      external src_index
c
c      if (intopt.eq.7) then
c        write(71,*) 'SRC_INT_AB:', npts, intopt, sbnds(0)
c        do ik = 1, npts
c           write(71,*) ik, quant(ik), spts(ik), sbnds(ik)
c           write(0,*) ik, quant(ik), spts(ik), sbnds(ik)
c        enddo
c      endif
c
      srcsum  = 0.0
c
c     Calculate the integral in different ways depending on options
c
      ik1 = src_index(s1,spts,sbnds,npts)
      ik2 = src_index(s2,spts,sbnds,npts)
c
c     Integral over step-wise source function
c
      if (intopt.eq.0) then
c
c
c        Two points in same cell
c
         if (ik1.eq.ik2) then
c
            srcsum = quant(ik1) * (s2-s1)
c
         else
c
            do ik = ik1,ik2
c
               if (ik.eq.ik1) then
c
                  srcsum = srcsum + quant(ik) * (sbnds(ik) - s1)
c
               elseif (ik.eq.ik2) then
c
                  srcsum=srcsum + quant(ik) * (s2 - sbnds(ik-1))
c
               else
c
                  srcsum=srcsum+quant(ik)*(sbnds(ik)-sbnds(ik-1))
c
               endif
c
            end do
c
         endif
c
c     Integration over Interpolated source function
c
      elseif (intopt.eq.1) then
c
c        Two points in same cell
c

c         write (6,*) 'ik:',ik1,ik2,s1,s2
c
         if (ik1.eq.ik2) then
c
c           Both in first or second half
c
            if ((s1.le.spts(ik1).and.s2.le.spts(ik1)).or.
     >          (s1.ge.spts(ik1).and.s2.ge.spts(ik1))) then
c
               sval1 = (s1+s2)/2.0d0
               srcsum = (s2-s1) * src_val_ik(quant,sval1,ik1,
     >                                       spts,sbnds,npts,intopt)
c
c           S1 in first half and S2 in second.
c
            else
c
               sval1 = ((s2+spts(ik1))/2.0d0)
               sval2 = ((s1+spts(ik1))/2.0d0)
c
               srcsum = (s2-spts(ik1))
     >            * src_val_ik(quant,sval1,ik1,
     >                         spts,sbnds,npts,intopt)
     >                +  (spts(ik1)-s1)
     >            * src_val_ik(quant,sval2,ik1,
     >                         spts,sbnds,npts,intopt)
c
            endif
c
c        Points not in same cell
c
         else
c
            do ik = ik1,ik2
c
               if (ik.eq.ik1) then
c
                  if (s1.le.spts(ik)) then
c
                      sval1 = ((s1+spts(ik))/2.0d0)
                      sval2 = (sbnds(ik)+spts(ik))/2.0d0
c
                      srcsum = srcsum + (spts(ik)-s1)
     >                   * src_val_ik(quant,sval1,ik,
     >                                spts,sbnds,npts,intopt)
     >                  + (sbnds(ik)-spts(ik)) *
     >                src_val_ik(quant,sval2,ik,
     >                           spts,sbnds,npts,intopt)
c
c                      write (6,*) '1a: ',ik,srcsum
c                      write(6,*)  'num:',s1,sval1,spts(ik),
c     >                           sval2,sbnds(ik)
c
                  else
c
                      sval1 =  (sbnds(ik)+s1)/2.0d0
c
                      srcsum = srcsum
     >                  + (sbnds(ik)-s1) *
     >                src_val_ik(quant,sval1,ik,
     >                           spts,sbnds,npts,intopt)
c
c                      write (6,*) '1b:',ik,srcsum,sval1
c                      write(6,*)  'num:',spts(ik),s1,
c     >                           sval1,sbnds(ik)
c
                  endif
c
               elseif (ik.eq.ik2) then
c
                  if (s2.le.spts(ik)) then
c
                      sval1 = (sbnds(ik-1)+s2)/2.0d0
c
                      srcsum = srcsum
     >                  + (s2-sbnds(ik-1)) *
     >                src_val_ik(quant,sval1,ik,
     >                           spts,sbnds,npts,intopt)
c
c                      write (6,*) 'Na:',ik,srcsum,sval1
c                      write(6,*)  'num:',sbnds(ik-1),sval1,s2,
c     >                           spts(ik)
c
                  else
c
                      sval1 = ((sbnds(ik-1)+spts(ik))/2.0d0)
                      sval2 = (s2+spts(ik))/2.0d0
c
                      srcsum = srcsum + (spts(ik)-sbnds(ik-1))
     >                     * src_val_ik(quant,
     >                             sval1,ik,
     >                              spts,sbnds,npts,intopt)
     >                  + (s2-spts(ik)) *
     >                       src_val_ik(quant,sval2,ik,
     >                                  spts,sbnds,npts,intopt)
c
c                      write (6,*) 'Nb:',ik,srcsum,sval1,sval2
c                      write(6,*)  'num:',sbnds(ik-1),sval1,spts(ik),
c     >                           sval2,s2
c
                  endif
c
               else
c
                  sval1 = (sbnds(ik-1)+spts(ik))/2.0d0
                  sval2 = (sbnds(ik)+spts(ik))/2.0d0
c
                  srcsum = srcsum
     >                  + (spts(ik)-sbnds(ik-1)) *
     >                src_val_ik(quant,sval1,ik,
     >                           spts,sbnds,npts,intopt)
     >                  + (sbnds(ik)-spts(ik)) *
     >                src_val_ik(quant,sval2,ik,
     >                           spts,sbnds,npts,intopt)
c
c                  write (6,*) 'ik:',ik,srcsum,sval1,sval2
c                      write(6,*)  'num:',sbnds(ik-1),sval1,spts(ik),
c     >                           sval2,sbnds(ik)
c
               endif
c
            end do
c
         endif
c
      endif
c
c

c
      src_int_ab = srcsum
c
c      write (6,'(a,2i4,5g12.4)') 'srcab:',ik1,ik2,src_int_ab
c
c
c
      return
      end
c
c
c
      subroutine TwoD_smooth(Q_old, eps_s, eps_r)
      use mod_params
      use mod_cgeom
      implicit none
c
c     include 'params'
c
      real eps_s, eps_r
      real Q_old(maxnks,maxnrs), Q_new(maxnks,maxnrs)
c
c     include 'cgeom'
      real Q_new_s, Q_new_r, A_s, A_r
      integer ia,ib,ic,ja,jb,jc,ik,ir
c
      do ir = irsep, irwall-1
        do ik = 1, nks(ir)
          if (ik.ne.1) then
            ia = ik-1
          else
            ia = ik
          endif
          ja = ir
          ib = ik
          jb = ir
          if (ik.ne.nks(ir)) then
            ic = ik+1
          else
            ic = ik
          endif
          jc = ir
c
          Q_new_s =  Q_old(ia,ja) - 2.0 * Q_old(ib,jb)
     >            +  Q_old(ic,jc)
c
c          Q_new_s = ( kareas(ia,ja) * Q_old(ia,ja)
c     >              - kareas(ib,jb) * Q_old(ib,jb) * 2.0
c     >              + kareas(ic,jc) * Q_old(ic,jc) )
c          Q_new_s = Q_new_s * 4.0 / (kareas(ia,ja) +
c     >                2.0 * kareas(ib,jb) + kareas(ic,jc))
c
          if (ir.ne.irwall-1) then
            ia = ikouts(ik,ir)
            ja = irouts(ik,ir)
          else
            ia = ik
            ja = ir
          endif
          ib = ik
          jb = ir
          if (ir.ne.irsep) then
            ic = ikins(ik,ir)
            jc = irins(ik,ir)
          else
            ic = ik
            jc = ir
          endif
          Q_new_r = Q_old(ia,ja) - 2.0 * Q_old(ib,jb)
     >            + Q_old(ic,jc)
c
c          Q_new_r = ( kareas(ia,ja) * Q_old(ia,ja)
c     >              - kareas(ib,jb) * Q_old(ib,jb) * 2.0
c     >              + kareas(ic,jc) * Q_old(ic,jc) )
c          Q_new_r = Q_new_r * 4.0 / (kareas(ia,ja) +
c     >                2.0 * kareas(ib,jb) + kareas(ic,jc))
c
          if (eps_r.ge.0.0) then
             Q_new(ik,ir) = Q_old(ik,ir) + 0.5D0 *
     >                   (eps_r * Q_new_r  + eps_s * Q_new_s)
          else
             Q_new(ik,ir) = Q_old(ik,ir) + eps_s * Q_new_s
          endif
        enddo
      enddo
c
      call TwoD_assign(Q_new, Q_old, 1.0)
c
      return
      end
c
c
      subroutine TwoD_assign(Q_1, Q_2, const)
      use mod_params
      use mod_cgeom
      implicit none
c
c     Assign Q_1 to Q_2
c
c     include 'params'
      real Q_1(maxnks,maxnrs), Q_2(maxnks,maxnrs),
     >     const
c
c     include 'cgeom'
      integer ik,ir
c
c      write(71,*) 'ass.const=',const
c
      do ir = 1, nrs
        do ik = 1, nks(ir)
           Q_2(ik,ir) = const * Q_1(ik,ir)
        enddo
      enddo
c
      return
      end
c
c
      subroutine TwoD_relax(Q_new, Q_old, delta)
      use mod_params
      use mod_cgeom
      implicit none
c
c     relax Q_old and Q_new into Q_old
c
c     include 'params'
      real Q_old(maxnks,maxnrs), Q_new(maxnks,maxnrs)
c
      real delta(maxnrs)
      real Q_temp(maxnks,maxnrs)
c
c
c     include 'cgeom'
      integer ik,ir
c
      do ir = 1, nrs
        do ik = 1, nks(ir)
           Q_temp(ik,ir) = Q_old(ik,ir) * (1.0 - delta(ir)) +
     >                     Q_new(ik,ir) * delta(ir)
        enddo
      enddo
c
      call TwoD_assign(Q_temp, Q_old, 1.0)
c
      return
      end
c
c
      subroutine TwoD_compare(Q_new, Q_old)
      use mod_params
      use mod_cgeom
      implicit none
c
c     relax Q_old and Q_new into Q_old
c
c     include 'params'
      real Q_old(maxnks,maxnrs), Q_new(maxnks,maxnrs)
c
      real temp1, temp2
c
c     include 'cgeom'
      integer ik,ir,ik_temp
c
      do ir = irsep, irwall
        temp1 = 0.0
        temp2 = 0.0
        do ik = 1, nks(ir)
           temp1 = abs((Q_old(ik,ir) - Q_new(ik,ir))
     >                   / Q_old(ik,ir))
           if (temp1.ge.temp2) then
              temp2 = temp1
              ik_temp = ik
           endif
        enddo
        write(71,*) 'comp:', ir, ik_temp, temp2
      enddo
c
      return
      end
c
c
      subroutine TwoD_zero(Q_1)
      use mod_params
      use mod_cgeom
      implicit none
c
c     include 'params'
      real Q_1(maxnks,maxnrs)
c
c     include 'cgeom'
      integer ik,ir
c
      do ir = 1, nrs
        do ik = 1, nks(ir)
           Q_1(ik,ir) = 0.0
        enddo
      enddo
c
      return
      end
c
c
      subroutine ups_zero(Q_1)
      use mod_params
      use mod_cgeom
      implicit none
c
c     include 'params'
      real Q_1(maxnks,maxnrs)
c
c     include 'cgeom'
      integer ik,ir
c
      do ir = 1, nrs
        do ik = 1, nks(ir)
           if (kss(ik,ir).ge.ksmaxs(ir)*0.25.and.
     >         kss(ik,ir).le.ksmaxs(ir)*0.75) then
               Q_1(ik,ir) = 0.0
           endif
        enddo
      enddo
c
      return
      end
c
c
      subroutine TwoD_momloss(irlim1, irlim2)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_sol23_input
      use mod_pindata
      use mod_pin_cfd
      use mod_sol23_com
      use mod_cx
      implicit none
c
      integer irlim1, irlim2, ik,ir
c     sazmod - can;t find so commenting out, i.e. breaking sol23
c      real get_bg_mass
      real*8 mass_tmp 

c      external get_bg_mass

c
c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c     include 'sol23_input'
c     include 'pindata'
c     include 'pin_cfd'
c     include 'sol23_com'
c
      real*8  sigmavcx, sigma, Ti_tmp,
     >    vp_tmp, Eh_tmp, nh_tmp, ga_tmp, wn_iso
      real*8   Q_CX
c
c      wn_iso = 1.0/sqrt(3.0d0)
c
      wn_iso = 0.57735d0
c
      do ir = irlim1, irlim2
        do ik = 1,nks(ir)
c
        Eh_tmp = pinena(ik,ir)
        nh_tmp = pinatom(ik,ir)
c
        ga_tmp = - knbs(ik,ir) * kvhs(ik,ir)
        vp_tmp = abs(kvhs(ik,ir))
        Ti_tmp = ktibs(ik,ir)
c
c        mass_tmp = get_bg_mass()
c
        call cxsig(1.0d0,mass_tmp, Eh_tmp, wn_iso, wn_iso,
     >    wn_iso, Ti_tmp, vp_tmp, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1,
     >    sigma, sigmavcx)
c
        if (ga_tmp.eq.0.0.or.nh_tmp.eq.0.0) then
           Q_CX = 0.0d0
        else
           Q_CX = crmb * amu * ga_tmp * nh_tmp * sigmavcx
        endif
c
        pinmom_cx (ik,ir) = Q_CX
        pinmom_rec(ik,ir) = - crmb * amu * kvhs(ik,ir)
     >                             * pinrec_sm(ik,ir)
        if (s23_par_pinmom.eq.0) then
          pinmom(ik,ir) = pinmom_cx(ik,ir) + pinmom_rec(ik,ir)
          pinmom(ik,ir) = s23_par_vnmult * pinmom(ik,ir)
        else
          pinmom(ik,ir) = pinmp(ik,ir)
        endif
c
        enddo
      enddo
c
      return
      end
c
c
      subroutine perp_diffusion(d2ndr)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_sol23_input
      implicit none
c
c     include 'params'
      real d2ndr(maxnks,maxnrs)
c
      integer ik, ir, ikout, irout, ikin, irin
      real dr_a, dr_b, dr_ab
c
c     include 'cgeom'
c     include 'comtor'
c     include 'sol23_input'
c
      do ir = irsep,irwall-1
        do ik = 1,nks(ir)
c
          if (ir.eq.irsep) then
            dr_b  = koutds(ik,ir)
            ikout = ikouts(ik,ir)
            irout = irouts(ik,ir)
c
            if (dr_b.ne.0.0) then
              d2ndr(ik,ir)  =
     >         (knbs(ikout,irout) - knbs(ik,ir)) / dr_b
            else
              d2ndr(ik,ir)  = 0.0
            endif
c
          elseif (ir.eq.irwall-1) then
            dr_a  = kinds(ik,ir)
            ikin  = ikins(ik,ir)
            irin  = irins(ik,ir)
c
            if (dr_a.ne.0.0) then
              d2ndr(ik,ir)  =
     >          (knbs(ikin,irin) - knbs(ik,ir)) / dr_a
            else
              d2ndr(ik,ir)  = 0.0
            endif
          else
            dr_a  = kinds(ik,ir)
            ikin  = ikins(ik,ir)
            irin  = irins(ik,ir)
c
            dr_b  = koutds(ik,ir)
            ikout = ikouts(ik,ir)
            irout = irouts(ik,ir)
            dr_ab = (dr_a + dr_b) / 2.0
c
            if (dr_a.ne.0.0.and.dr_b.ne.0.0) then
              d2ndr(ik,ir)  =
     >         ( (knbs(ikout,irout) - knbs(ik,ir)) / dr_b
     >         - (knbs(ikin,irin) - knbs(ik,ir)) / dr_a ) / dr_ab
            else
              d2ndr(ik,ir)  = 0.0
            endif
          endif
c
        enddo
      enddo
c
      do ir = irtrap + 1, nrs
        do ik = 1, nks(ir)
          d2ndr(ik,ir) = 1.0
        enddo
      enddo
c
      return
      end
c
c
      subroutine pinqi_1234(ir1,ir2)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_sol23_input
      use mod_pindata
      use mod_pin_cfd
      use mod_sol23_com
      use mod_cx
      implicit none
c
      integer ir1,ir2
      real get_bg_mass
      real*8 mass_tmp 

c      external get_bg_mass

c
c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c     include 'sol23_input'
c     include 'pindata'
c     include 'pin_cfd'
c     include 'sol23_com'
c
      integer ik,ir
      real*8  sigmavcx, sigma, T_ion, n_atom,
     > v_ion, E_ion, E_atom, E_molec, n_ion, wn_iso
      real  f_molec, f_atom,  sum_e,
     >      sum_i_new, sum_i_old, sum_1, sum_2, sum_3,
     >      sum_4, sum_5, sum_6, temp_9
c
      sum_i_new = 0.0
      sum_i_old = 0.0
      sum_e     = 0.0
      sum_1     = 0.0
      sum_2     = 0.0
      sum_3     = 0.0
      sum_4     = 0.0
      sum_5     = 0.0
      sum_6     = 0.0
c
      do ir = ir1, ir2
        do ik = 1, nks(ir)
c
c       wn_iso = 1.0/sqrt(3.0d0)
c
        wn_iso = 0.57735d0
        E_atom = pinena(ik,ir)
        E_molec= s23_par_Emolec * 1.0d0
        n_atom = pinatom(ik,ir)
c
        n_ion  = knbs(ik,ir)
        v_ion  = abs(kvhs(ik,ir))
        T_ion  = ktibs(ik,ir)
        E_ion  = 1.5 * T_ion
        if ((pinatom(ik,ir) +
     >       pinmol(ik,ir)).ne.0.0) then
           f_atom = pinatom(ik,ir) /
     >          (pinatom(ik,ir) + pinmol(ik,ir))
        else
           f_atom = 0.5
        endif
        f_molec = 1.0 - f_atom
c
c       THERMAL LOSS DUE TO ATOMIC IONIZ.
c
        pinqi_1(ik,ir) = f_atom * E_atom * ech
     >                          * pinion_sm(ik,ir)
c
c       THERMAL LOSS DUE TO MOLECULAR IONIZ.
c
        pinqi_2(ik,ir) = f_molec * E_molec * ech
     >                           * pinion_sm(ik,ir)
c
c       THERMAL LOSS DUE TO RECOMBINATION
c
        pinqi_3(ik,ir) = - E_ion * ech
     >                           * pinrec_sm(ik,ir)
c
c       THERMAL LOSS DUE TO CHARGE EXCHANGE (CX EQUIPARTITION)
c
c        mass_tmp=get_bg_mass()
        call cxsig(1.0d0,mass_tmp, E_atom, wn_iso, wn_iso,
     >    wn_iso, T_ion, v_ion, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1,
     >    sigma, sigmavcx)
c
c       if (ir.eq.14) then
c          write(71,'(2I5,5G12.6)') ir, ik, sigmavcx,
c     >              E_atom, E_ion, n_atom, n_ion
c       endif
c
        if (n_atom.eq.0.0.or.n_ion.eq.0.0) then
          pinqi_4(ik,ir) = 0.0
        else
          pinqi_4(ik,ir) = (E_atom - E_ion) * ech *
     >                           n_ion * n_atom * sigmavcx
        endif
c
c       CONVECTIVE LOSS DUE TO RECOMBINATION
c
        pinqi_5(ik,ir) = - 0.5 * crmb * amu *
     >                   kvhs(ik,ir) ** 2.0 * pinrec_sm(ik,ir)
c
c       CONVECTIVE LOSS DUE TO CHARGE EXCHANGE
c
        if (n_atom.eq.0.0.or.n_ion.eq.0.0) then
          pinqi_6(ik,ir) = 0.0
        else
          pinqi_6(ik,ir) = - 0.5 * crmb * amu * kvhs(ik,ir) ** 2.0 *
     >                n_ion * n_atom * sigmavcx
        endif
c
c       OLD PINQI SUM [NIMBUS]
c
        sum_i_old  = sum_i_old  + pinqi(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
c
        if (s23_par_pinqiflag.eq.0) then
          pinqi(ik,ir) =  pinqi_1(ik,ir) + pinqi_2(ik,ir) +
     >      pinqi_3(ik,ir) + pinqi_4(ik,ir) + pinqi_5(ik,ir) +
     >      pinqi_6(ik,ir)
        endif
c
c        pinqi(ik,ir) =  pinqi_1(ik,ir) + pinqi_2(ik,ir)
c
c      ADD RECOMBINATION HEATING
c
        pinqe_sm(ik,ir) = pinqe_sm(ik,ir)
     >                    + s23_par_rec_heat * ech * pinrec_sm(ik,ir)
c
c       NEW PINQI SUM [1 to 6]
c
        sum_i_new = sum_i_new + pinqi(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
        sum_1    = sum_1    + pinqi_1(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
        sum_2    = sum_2    + pinqi_2(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
        sum_3    = sum_3    + pinqi_3(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
        sum_4    = sum_4    + pinqi_4(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
        sum_5    = sum_5    + pinqi_5(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
        sum_6    = sum_6    + pinqi_6(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
        sum_e    = sum_e    + pinqe_sm(ik,ir) * karea2(ik,ir) *
     >                krb(ik,ir) * 2.0 * 3.14 * 1.0E-6
c
        pinqi(ik,ir) = s23_par_pinqimult * pinqi(ik,ir)
c
        enddo
      enddo
c
      s23_pinqe = abs(sum_e)
c
      write(71,'(A30,G16.8)') '     TOTAL PINQE        [MW]:', sum_e
      write(71,'(A30,G16.8)') '     TOTAL PINQI NIMBUS [MW]:', sum_i_old
      write(71,'(A30,G16.8)') '     TOTAL PINQI 123456 [MW]:', sum_i_new
      write(71,'(A25,G16.8)') '      atomic iz. [MW]:', sum_1
      write(71,'(A25,G16.8)') '   molecular iz. [MW]:', sum_2
      write(71,'(A25,G16.8)') '  thermal recomb.[MW]:', sum_3
      write(71,'(A25,G16.8)') '  convect recomb.[MW]:', sum_5
      write(71,'(A25,G16.8)') '  thermal CX     [MW]:', sum_4
      write(71,'(A25,G16.8)') '  convect CX     [MW]:', sum_6
      write(71,'(A30,G16.8)') '  PINQI MULTIPLIER          :', temp_9
      write(71,'(A30,G16.8)') '  PINQE MULTIPLIER          :',
     >                                              1.0 / s23_zhi0
c
c      write(71,'(A10,5G16.8)') 'sum:i/e', sum_1234, sum_e,
c     >    sum_core, sum_1234/sum_e, (sum_1234+sum_e)/sum_core
c      write(71,*) 'PINQI mult:', 0.0D0
c
c      if (abs(sum_1234/sum_e).ge.1.0.and.sum_1234.ge.0.0) then
c        do ir = ir1, ir2
c          do ik = 1, nks(ir)
c            pinqi(ik,ir) = pinqi(ik,ir)
c     >                      * 1.0 / abs(sum_1234/sum_e)
c          enddo
c        enddo
c        write(71,*) 'PINQI >= PINQE, clipped'
c      endif
c
      return
      end
c
      REAL FUNCTION Calcwidth2(ik,ir)
      use mod_params
      use mod_cgeom
      use mod_slcom
      IMPLICIT none
c
c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'slcom'
c
      INTEGER ik,ir,mode
c
      REAL    s,r,z,r1,z1,r2,z2,c1,c2,wid(3)
      REAL    z3,z4,r3,r4
      REAL    widtop,widbot,cbot1,ctop1,ctop2,cbot2
      REAL    deltar,deltaz
      REAL    rside1,zside1,rside2,zside2
      INTEGER id
c
      id = korpg(ik,ir)
c
      IF (ir.EQ.1.OR.ir.EQ.irwall.OR.ir.EQ.irtrap) THEN
        WRITE(0,*) 'ERROR (Calcwidth2): Virtual ring (IR = ',ir,')'
        STOP
      ENDIF
c
      IF (nvertp(id).NE.4) THEN
        write(0,*) 'ERROR (Calcwidth2): Invalid cell',ik,ir,id,
     >   nvertp(id)
        STOP 'ERROR (Calcwidth2): Invalid cell'
      ENDIF
c
      r = rs(ik,ir)
      z = zs(ik,ir)
c
c       Perpendicular distance from point on center line to
c       each side:
c
c       Side 1-4:
c
      deltar = rvertp(4,id) - rvertp(1,id)
      deltaz = zvertp(4,id) - zvertp(1,id)
c
      c1 = ((r - rvertp(1,id)) * deltar +
     >        (z - zvertp(1,id)) * deltaz) / (deltar**2 + deltaz**2)
c
      r1 = rvertp(1,id) + c1 * deltar
      z1 = zvertp(1,id) + c1 * deltaz
c
c       Side 2-3:
c
      deltar = rvertp(3,id) - rvertp(2,id)
      deltaz = zvertp(3,id) - zvertp(2,id)
c
      c2 = ((r - rvertp(2,id)) * deltar +
     >        (z - zvertp(2,id)) * deltaz) / (deltar**2 + deltaz**2)
c
      r2 = rvertp(2,id) + c2 * deltar
      z2 = zvertp(2,id) + c2 * deltaz
c
c       Calculate distances:
c
      Calcwidth2 = SQRT((r - r1)**2.0 + (z - z1)**2.0) +
     >             SQRT((r - r2)**2.0 + (z - z2)**2.0)
c
      RETURN
99    STOP
      END
c
c
c
      subroutine print_sol_excel
      use mod_params
      use mod_sol23_com
      use mod_cgeom
      use mod_comtor
      implicit none
c     include 'params'
c     include 'sol23_com'
c     include 'cgeom'
c     include 'comtor'
c
c     print_sol_excel: Print out SOL information from SOL22 or 23
c                      in a format that Wojciech prefers for
c                      creating Excel spreadsheets.
c
c
      integer i_a,i_b,i_c,i_d,ik,ir
      real temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7,
     >     temp_8,temp_9
      real Mach_a,Mach_b, Mach_c, Mach_d
c
        i_a = s23_par_prring0
        i_b = s23_par_prring0 + s23_par_prring1
c
c       The following two lines look to be incorrect. 
c       Set i_c and i_d to the values listed. 
c
c       i_b = s23_par_prring0 + 2 * s23_par_prring1
c       i_b = s23_par_prring0 + 3 * s23_par_prring1
c
c        Replace with i_c and i_d as Steve mentions 

        i_c = s23_par_prring0 + 2 * s23_par_prring1
        i_d = s23_par_prring0 + 3 * s23_par_prring1
c slmod begin - new
c...BUG? Should the above code read?
c        i_c = s23_par_prring0 + 2 * s23_par_prring1
c        i_d = s23_par_prring0 + 3 * s23_par_prring1
c...ARRAY BOUNDS:
        i_a = MAX(1,MIN(nrs,i_a))
        i_b = MAX(1,MIN(nrs,i_b))
        i_c = MAX(1,MIN(nrs,i_c))
        i_d = MAX(1,MIN(nrs,i_d))
c slmod end
c
        temp_7 = 1.0E20
        temp_8 = 1.0E23
        temp_9 = 1.0E6
c
c       kps(ik) poloidal s
c       kss(ik)  along-B s||
c
c        write(0,*) 'KSMAX', ksmaxs(i_a),  ksmaxs(i_b),
c     >                      ksmaxs(i_c),  ksmaxs(i_d)
c
           do ik = 1, nks(i_a)
             write(73,'(8G16.8)')   kss(ik,i_a) / ksmaxs(i_a),
     >         knbs(ik,i_a)/temp_7, kss(ik,i_b) / ksmaxs(i_b),
     >         knbs(ik,i_b)/temp_7, kss(ik,i_c) / ksmaxs(i_c),
     >         knbs(ik,i_c)/temp_7, kss(ik,i_d) / ksmaxs(i_d),
     >         knbs(ik,i_d)/temp_7
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
            write(73,'(8G16.8)')  kss(ik,i_a) / ksmaxs(i_a),
     >         knbs(ik,i_a)*kvhs(ik,i_a)/temp_8,
     >         kss(ik,i_b) / ksmaxs(i_b),
     >         knbs(ik,i_b)*kvhs(ik,i_b)/temp_8,
     >         kss(ik,i_c) / ksmaxs(i_c),
     >         knbs(ik,i_c)*kvhs(ik,i_c)/temp_8,
     >         kss(ik,i_d) / ksmaxs(i_d),
     >         knbs(ik,i_d)*kvhs(ik,i_d)/temp_8
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
             write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         ktibs(ik,i_a),     kss(ik,i_b) / ksmaxs(i_b),
     >         ktibs(ik,i_b),     kss(ik,i_c) / ksmaxs(i_c),
     >         ktibs(ik,i_c),     kss(ik,i_d) / ksmaxs(i_d),
     >         ktibs(ik,i_d)
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
             write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         ktebs(ik,i_a),     kss(ik,i_b) / ksmaxs(i_b),
     >         ktebs(ik,i_b),     kss(ik,i_c) / ksmaxs(i_c),
     >         ktebs(ik,i_c),     kss(ik,i_d) / ksmaxs(i_d),
     >         ktebs(ik,i_d)
           enddo
           write(73,*)
c
          do ik = 1, nks(i_a)
c
           Mach_a = kvhs(ik,i_a) /
     >           sqrt((ktibs(ik,i_a) + ktebs(ik,i_a))/crmb*emi)
           Mach_b = kvhs(ik,i_b) /
     >           sqrt((ktibs(ik,i_b) + ktebs(ik,i_b))/crmb*emi)
           Mach_c = kvhs(ik,i_c) /
     >           sqrt((ktibs(ik,i_c) + ktebs(ik,i_c))/crmb*emi)
           Mach_d = kvhs(ik,i_d) /
     >           sqrt((ktibs(ik,i_d) + ktebs(ik,i_d))/crmb*emi)
c
           temp_1 = knbs(ik,i_a) * (ktibs(ik,i_a) + ktebs(ik,i_a)) * ech
           temp_3 = temp_1 * (1.0 + Mach_a ** 2.0)
c
           temp_1 = knbs(ik,i_b) * (ktibs(ik,i_b) + ktebs(ik,i_b)) * ech
           temp_4 = temp_1 * (1.0 + Mach_b ** 2.0)
c
           temp_1 = knbs(ik,i_c) * (ktibs(ik,i_c) + ktebs(ik,i_c)) * ech
           temp_5 = temp_1 * (1.0 + Mach_c ** 2.0)
c
           temp_1 = knbs(ik,i_d) * (ktibs(ik,i_d) + ktebs(ik,i_d)) * ech
           temp_6 = temp_1 * (1.0 + Mach_d ** 2.0)
c
           write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         temp_3,     kss(ik,i_b) / ksmaxs(i_b),
     >         temp_4,     kss(ik,i_c) / ksmaxs(i_c),
     >         temp_5,     kss(ik,i_d) / ksmaxs(i_d),
     >         temp_6
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
             Mach_a = kvhs(ik,i_a) /
     >           sqrt((ktibs(ik,i_a) + ktebs(ik,i_a))/crmb*emi)
             Mach_b = kvhs(ik,i_b) /
     >           sqrt((ktibs(ik,i_b) + ktebs(ik,i_b))/crmb*emi)
             Mach_c = kvhs(ik,i_c) /
     >           sqrt((ktibs(ik,i_c) + ktebs(ik,i_c))/crmb*emi)
             Mach_d = kvhs(ik,i_d) /
     >           sqrt((ktibs(ik,i_d) + ktebs(ik,i_d))/crmb*emi)
             write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         Mach_a,     kss(ik,i_b) / ksmaxs(i_b),
     >         Mach_b,     kss(ik,i_c) / ksmaxs(i_c),
     >         Mach_c,     kss(ik,i_d) / ksmaxs(i_d),
     >         Mach_d
           enddo
           write(73,*)
c
c       call rinout ('W E2D N ',e2dnbs,maxnks*maxnrs)
c       call rinout ('W E2D TE',e2dtebs,maxnks*maxnrs)
c       call rinout ('W E2D TI',e2dtibs,maxnks*maxnrs)
c       call rinout ('W E2D VB',e2dvhs,maxnks*maxnrs)
c       call rinout ('W E2D E ',e2des,maxnks*maxnrs)
c       call rinout ('W E2D I ',e2dion,maxnks*maxnrs)
c       call rinout ('W E2D A ',e2datom,maxnks*maxnrs)
c       call rinout ('W E2D TA',e2dtarg,maxnrs*8*2)
c       call rinout ('W E2D GP',e2dgpara,maxnks*maxnrs)
c       call rinout ('W E2D GD',e2dgdown,maxnks*maxnrs)
c       call rinout ('W E2D G ',e2dflux,(maxnks+1)*maxnrs)
c       call rinout ('W E2D VE',e2dbvel,(maxnks+1)*maxnrs)
c
c           i_a = 48
c           i_b = 38
c           i_c = 36
c           i_d = 27
c
c           do ir = 1, 7
c              write(73,'(I5,4G16.8)') ir, knbs(20,ir), ktebs(20,ir),
c     >           ktibs(20,ir), kvhs(20,ir)
c           enddo
c
c           do ir = 8, 16
c             write(73,'(I5,2G16.8)') ir,
c     >         power_flux(ir,2)/temp_9,  powere_targ(ir,2)/temp_9
c          enddo
c           write(73,*)
c
c           do ir = 8, 16
c             write(73,'(I5,4G16.8)') ir,
c     >         knbs(i_a,ir)/temp_7,  knbs(i_b,ir)/temp_7,
c     >         knbs(i_c,ir)/temp_7,  knbs(i_d,ir)/temp_7
c           enddo
c           write(73,*)
c
c           do ir = 8, 16
c             write(73,'(I5,4G16.8)') ir,
c     >         ktebs(i_a,ir),  ktebs(i_b,ir),
c     >         ktebs(i_c,ir),  ktebs(i_d,ir)
c           enddo
c           write(73,*)
c
c           do ir = 8, 16
c           temp_1 = knbs(i_a,ir) * (2.0 * ktebs(i_a,ir)) * ech
c           temp_2 = kvhs(i_a,ir) / sqrt(2.0*ktebs(i_a,ir)/crmb*emi)
c           temp_3 = temp_1 * (1.0 + temp_2 ** 2.0)
c
c           temp_1 = knbs(i_b,ir) * (2.0 * ktebs(i_b,ir)) * ech
c           temp_2 = kvhs(i_b,ir) / sqrt(2.0*ktebs(i_b,ir)/crmb*emi)
c           temp_4 = temp_1 * (1.0 + temp_2 ** 2.0)
c
c           temp_1 = knbs(i_c,ir) * (2.0 * ktebs(i_c,ir)) * ech
c           temp_2 = kvhs(i_c,ir) / sqrt(2.0*ktebs(i_c,ir)/crmb*emi)
c           temp_5 = temp_1 * (1.0 + temp_2 ** 2.0)
c
c           temp_1 = knbs(i_d,ir) * (2.0 * ktebs(i_d,ir)) * ech
c           temp_2 = kvhs(i_d,ir) / sqrt(2.0*ktebs(i_d,ir)/crmb*emi)
c           temp_6 = temp_1 * (1.0 + temp_2 ** 2.0)
c           write(73,'(I5,4G16.8)') ir, temp_3,
c     >         temp_4, temp_5, temp_6
c           enddo
c
c
         return
c
      end
c
c
      subroutine print_e2d_excel
      use mod_params
      use mod_sol23_com
      use mod_cgeom
      use mod_comtor
      use mod_cedge2d
      implicit none
c     include 'params'
c     include 'sol23_com'
c     include 'cgeom'
c     include 'comtor'
c     include 'cedge2d'
c
c     print_sol_excel: Print out SOL information from SOL22 or 23
c                      in a format that Wojciech prefers for
c                      creating Excel spreadsheets.
c
c
      integer i_a,i_b,i_c,i_d,ik,ir
      real temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7,
     >     temp_8,temp_9
      real Mach_a,Mach_b, Mach_c, Mach_d
c
        i_a = s23_par_prring0
        i_b = s23_par_prring0 + s23_par_prring1
c
c       The following two lines look to be incorrect. 
c       Set i_c and i_d to the values listed. 
c
c       i_b = s23_par_prring0 + 2 * s23_par_prring1
c       i_b = s23_par_prring0 + 3 * s23_par_prring1
c
c       Replace with i_c and i_d as Steve mention in print_sol_excel

        i_c = s23_par_prring0 + 2 * s23_par_prring1
        i_d = s23_par_prring0 + 3 * s23_par_prring1
c
c       ARRAY BOUNDS CHECK:
c
        i_a = MAX(1,MIN(nrs,i_a))
        i_b = MAX(1,MIN(nrs,i_b))
        i_c = MAX(1,MIN(nrs,i_c))
        i_d = MAX(1,MIN(nrs,i_d))
c 
c
        temp_7 = 1.0E20
        temp_8 = 1.0E23
        temp_9 = 1.0E6
c
c       kps(ik) poloidal s
c       kss(ik)  along-B s||
c
c        write(0,*) 'KSMAX', ksmaxs(i_a),  ksmaxs(i_b),
c     >                      ksmaxs(i_c),  ksmaxs(i_d)
c
           do ik = 1, nks(i_a)
             write(73,'(8G16.8)')   kss(ik,i_a) / ksmaxs(i_a),
     >         e2dnbs(ik,i_a)/temp_7, kss(ik,i_b) / ksmaxs(i_b),
     >         e2dnbs(ik,i_b)/temp_7, kss(ik,i_c) / ksmaxs(i_c),
     >         e2dnbs(ik,i_c)/temp_7, kss(ik,i_d) / ksmaxs(i_d),
     >         e2dnbs(ik,i_d)/temp_7
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
            write(73,'(8G16.8)')  kss(ik,i_a) / ksmaxs(i_a),
     >         e2dnbs(ik,i_a)*e2dvhs(ik,i_a)/temp_8,
     >         kss(ik,i_b) / ksmaxs(i_b),
     >         e2dnbs(ik,i_b)*e2dvhs(ik,i_b)/temp_8,
     >         kss(ik,i_c) / ksmaxs(i_c),
     >         e2dnbs(ik,i_c)*e2dvhs(ik,i_c)/temp_8,
     >         kss(ik,i_d) / ksmaxs(i_d),
     >         e2dnbs(ik,i_d)*e2dvhs(ik,i_d)/temp_8
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
             write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         e2dtibs(ik,i_a),     kss(ik,i_b) / ksmaxs(i_b),
     >         e2dtibs(ik,i_b),     kss(ik,i_c) / ksmaxs(i_c),
     >         e2dtibs(ik,i_c),     kss(ik,i_d) / ksmaxs(i_d),
     >         e2dtibs(ik,i_d)
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
             write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         e2dtebs(ik,i_a),     kss(ik,i_b) / ksmaxs(i_b),
     >         e2dtebs(ik,i_b),     kss(ik,i_c) / ksmaxs(i_c),
     >         e2dtebs(ik,i_c),     kss(ik,i_d) / ksmaxs(i_d),
     >         e2dtebs(ik,i_d)
           enddo
           write(73,*)
c
          do ik = 1, nks(i_a)
c
           Mach_a = e2dvhs(ik,i_a) /
     >           sqrt((e2dtibs(ik,i_a) + e2dtebs(ik,i_a))/crmb*emi)
           Mach_b = e2dvhs(ik,i_b) /
     >           sqrt((e2dtibs(ik,i_b) + e2dtebs(ik,i_b))/crmb*emi)
           Mach_c = e2dvhs(ik,i_c) /
     >           sqrt((e2dtibs(ik,i_c) + e2dtebs(ik,i_c))/crmb*emi)
           Mach_d = e2dvhs(ik,i_d) /
     >           sqrt((e2dtibs(ik,i_d) + e2dtebs(ik,i_d))/crmb*emi)
c
           temp_1 = e2dnbs(ik,i_a) * (e2dtibs(ik,i_a) + e2dtebs(ik,i_a)) * ech
           temp_3 = temp_1 * (1.0 + Mach_a ** 2.0)
c
           temp_1 = e2dnbs(ik,i_b) * (e2dtibs(ik,i_b) + e2dtebs(ik,i_b)) * ech
           temp_4 = temp_1 * (1.0 + Mach_b ** 2.0)
c
           temp_1 = e2dnbs(ik,i_c) * (e2dtibs(ik,i_c) + e2dtebs(ik,i_c)) * ech
           temp_5 = temp_1 * (1.0 + Mach_c ** 2.0)
c
           temp_1 = e2dnbs(ik,i_d) * (e2dtibs(ik,i_d) + e2dtebs(ik,i_d)) * ech
           temp_6 = temp_1 * (1.0 + Mach_d ** 2.0)
c
           write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         temp_3,     kss(ik,i_b) / ksmaxs(i_b),
     >         temp_4,     kss(ik,i_c) / ksmaxs(i_c),
     >         temp_5,     kss(ik,i_d) / ksmaxs(i_d),
     >         temp_6
           enddo
           write(73,*)
c
           do ik = 1, nks(i_a)
           Mach_a = e2dvhs(ik,i_a) /
     >           sqrt((e2dtibs(ik,i_a) + e2dtebs(ik,i_a))/crmb*emi)
           Mach_b = e2dvhs(ik,i_b) /
     >           sqrt((e2dtibs(ik,i_b) + e2dtebs(ik,i_b))/crmb*emi)
           Mach_c = e2dvhs(ik,i_c) /
     >           sqrt((e2dtibs(ik,i_c) + e2dtebs(ik,i_c))/crmb*emi)
           Mach_d = e2dvhs(ik,i_d) /
     >           sqrt((e2dtibs(ik,i_d) + e2dtebs(ik,i_d))/crmb*emi)
             write(73,'(8G16.8)') kss(ik,i_a) / ksmaxs(i_a),
     >         Mach_a,     kss(ik,i_b) / ksmaxs(i_b),
     >         Mach_b,     kss(ik,i_c) / ksmaxs(i_c),
     >         Mach_c,     kss(ik,i_d) / ksmaxs(i_d),
     >         Mach_d
           enddo
           write(73,*)
c
c       call rinout ('W E2D N ',e2dnbs,maxnks*maxnrs)
c       call rinout ('W E2D TE',e2dtebs,maxnks*maxnrs)
c       call rinout ('W E2D TI',e2dtibs,maxnks*maxnrs)
c       call rinout ('W E2D VB',e2dvhs,maxnks*maxnrs)
c       call rinout ('W E2D E ',e2des,maxnks*maxnrs)
c       call rinout ('W E2D I ',e2dion,maxnks*maxnrs)
c       call rinout ('W E2D A ',e2datom,maxnks*maxnrs)
c       call rinout ('W E2D TA',e2dtarg,maxnrs*8*2)
c       call rinout ('W E2D GP',e2dgpara,maxnks*maxnrs)
c       call rinout ('W E2D GD',e2dgdown,maxnks*maxnrs)
c       call rinout ('W E2D G ',e2dflux,(maxnks+1)*maxnrs)
c       call rinout ('W E2D VE',e2dbvel,(maxnks+1)*maxnrs)
c
c           i_a = 48
c           i_b = 38
c           i_c = 36
c           i_d = 27
c
c           do ir = 1, 7
c              write(73,'(I5,4G16.8)') ir, knbs(20,ir), ktebs(20,ir),
c     >           ktibs(20,ir), kvhs(20,ir)
c           enddo
c
c           do ir = 8, 16
c             write(73,'(I5,2G16.8)') ir,
c     >         power_flux(ir,2)/temp_9,  powere_targ(ir,2)/temp_9
c          enddo
c           write(73,*)
c
c           do ir = 8, 16
c             write(73,'(I5,4G16.8)') ir,
c     >         knbs(i_a,ir)/temp_7,  knbs(i_b,ir)/temp_7,
c     >         knbs(i_c,ir)/temp_7,  knbs(i_d,ir)/temp_7
c           enddo
c           write(73,*)
c
c           do ir = 8, 16
c             write(73,'(I5,4G16.8)') ir,
c     >         ktebs(i_a,ir),  ktebs(i_b,ir),
c     >         ktebs(i_c,ir),  ktebs(i_d,ir)
c           enddo
c           write(73,*)
c
c           do ir = 8, 16
c           temp_1 = knbs(i_a,ir) * (2.0 * ktebs(i_a,ir)) * ech
c           temp_2 = kvhs(i_a,ir) / sqrt(2.0*ktebs(i_a,ir)/crmb*emi)
c           temp_3 = temp_1 * (1.0 + temp_2 ** 2.0)
c
c           temp_1 = knbs(i_b,ir) * (2.0 * ktebs(i_b,ir)) * ech
c           temp_2 = kvhs(i_b,ir) / sqrt(2.0*ktebs(i_b,ir)/crmb*emi)
c           temp_4 = temp_1 * (1.0 + temp_2 ** 2.0)
c
c           temp_1 = knbs(i_c,ir) * (2.0 * ktebs(i_c,ir)) * ech
c           temp_2 = kvhs(i_c,ir) / sqrt(2.0*ktebs(i_c,ir)/crmb*emi)
c           temp_5 = temp_1 * (1.0 + temp_2 ** 2.0)
c
c           temp_1 = knbs(i_d,ir) * (2.0 * ktebs(i_d,ir)) * ech
c           temp_2 = kvhs(i_d,ir) / sqrt(2.0*ktebs(i_d,ir)/crmb*emi)
c           temp_6 = temp_1 * (1.0 + temp_2 ** 2.0)
c           write(73,'(I5,4G16.8)') ir, temp_3,
c     >         temp_4, temp_5, temp_6
c           enddo
c
c
         return
c
      end
c
c
c
      subroutine read_sol23_params(ierr)
      use mod_params
      use mod_sol23_input
      implicit none
      integer ierr
c     include 'params'
c     include 'sol23_input'
c
c     READ_SOL23_PARAMS: This routine reads the SOL23 parameters
c     from the input file - it is invoked from IODIV. These have been
c     grouped together in the sol23 module in order to minimize
c     the impact of changes made to SOL 23 on teh rest of the code.
c
      call rdr(sol23_par_ptipte,.true.,0.0,.true.,1.0E8,
     >                                  'par_ptipte       ',ierr)
      call rdi(sol23_par_adaptnr,.true.,0  ,.true.,10   ,
     >                                  'par_adaptnr      ',ierr)
      call rdi(sol23_par_debugflag,.true.,0  ,.true.,10   ,
     >                                  'par_debugflag    ',ierr)
      call rdi(sol23_par_debugnr,.true.,0  ,.true.,10   ,
     >                                  'par_debugnr      ',ierr)
      call rdi(sol23_par_refresh,.true.,0  ,.true.,10000,
     >                                  'par_refresh      ',ierr)
      call rdr(sol23_par_artvisc2,.true.,0.0,.true.,1.0E8,
     >                                  'par_artvisc2     ',ierr)
      call rdr(sol23_par_artvisc4,.true.,0.0,.true.,1.0E8,
     >                                  'par_artvisc4     ',ierr)
      call rdr(sol23_par_dtf,.true.,0.0,.true.,1.0E8,
     >                                  'par_dtf          ',ierr)
      call rdr(sol23_par_dtg,.true.,0.0,.true.,1.0E36,
     >                                  'par_dtg          ',ierr)
      call rdr(sol23_par_grid0,.true.,0.0,.true.,1.0E8,
     >                                  'par_grid0        ',ierr)
      call rdr(sol23_par_gridexp,.true.,0.0,.true.,1.0E8,
     >                                  'par_gridexp      ',ierr)
c
      call rdi(sol23_par_itermax,.true.,0  ,.true.,10000,
     >                                  'par_itermax      ',ierr)
      call rdr(sol23_par_ga1      ,.true.,0.0,.true.,1.0E8,
     >                                  'par_ga1          ',ierr)
      call rdr(sol23_par_ga2    ,.true.,0.0,.true.,1.0E8,
     >                                  'par_ga2          ',ierr)
      call rdr(sol23_par_ga3    ,.true.,0.0,.true.,1.0E8,
     >                                  'par_ga3          ',ierr)
      call rdr(sol23_par_ga4     ,.true.,0.0,.true.,1.0E8,
     >                                  'par_ga4          ',ierr)
      call rdi(sol23_par_updtqpit,.true.,0  ,.true.,10   ,
     >                                  'par_updtqpit     ',ierr)
      call rdr(sol23_par_updtqp2 ,.true.,0.0,.true.,1.0E8,
     >                                  'par_updtqp2      ',ierr)
      call rdr(sol23_par_updtqp3,.true.,0.0,.true.,1.0E8,
     >                                  'par_updtqp3      ',ierr)
      call rdr(sol23_par_updtqpte,.true.,0.0,.true.,1.0E8,
     >                                  'par_updtqpte     ',ierr)
      call rdr(sol23_par_garelax,.true.,0.0,.true.,1.0E8,
     >                                  'par_garelax      ',ierr)
c
      call rdi(sol23_par_gaiter ,.true.,0  ,.true.,10000,
     >                                  'par_gaiter       ',ierr)
      call rdi(sol23_par_limitte  ,.true.,0  ,.true.,10   ,
     >                                  'par_limitte      ',ierr)
      call rdr(sol23_par_celldte,.true.,0.0,.true.,1.0E8,
     >                                  'par_celldte      ',ierr)
      call rdr(sol23_par_updtbcte,.true.,0.0,.true.,1.0E8,
     >                                  'par_updtbcte     ',ierr)
      call rdr(sol23_par_updtdel0,.true.,0.0,.true.,1.0E8,
     >                                  'par_updtdel0     ',ierr)
      call rdr(sol23_par_updtdel1,.true.,0.0,.true.,1.0E8,
     >                                  'par_updtdel1     ',ierr)
      call rdr(sol23_par_updtdelm,.true.,0.0,.true.,1.0E8,
     >                                  'par_updtdelm     ',ierr)
      call rdr(sol23_par_qbrelax,.true.,0.0,.true.,1.0E8,
     >                                  'par_qbrelax      ',ierr)
      call rdr(sol23_par_gridg   ,.true.,0.0,.true.,1.0E8,
     >                                  'par_gridg        ',ierr)
      call rdr(sol23_par_grid_dx0,.true.,0.0,.true.,1.0E8,
     >                                  'par_grid_dx0     ',ierr)
c
      call rdr(sol23_par_gae    ,.true.,0.0,.true.,1.0E8,
     >                                  'par_gae          ',ierr)
      call rdr(sol23_par_gai      ,.true.,0.0,.true.,1.0E8,
     >                                  'par_gai          ',ierr)
      call rdi(sol23_par_tectrl ,.true.,0  ,.true.,10   ,
     >                                  'par_tectrl       ',ierr)
      call rdi(sol23_par_drflag ,.true.,0  ,.true.,10   ,
     >                                  'par_drflag       ',ierr)
      call rdi(sol23_par_dsflag  ,.true.,0  ,.true.,10   ,
     >                                  'par_dsflag       ',ierr)
      call rdr(sol23_par_limrel1 ,.true.,0.0,.true.,1.0E8,
     >                                  'par_limrel1      ',ierr)
      call rdr(sol23_par_limrel2 ,.true.,0.0,.true.,1.0E8,
     >                                  'par_limrel2      ',ierr)
      call rdr(sol23_par_limrel3 ,.true.,0.0,.true.,1.0E8,
     >                                  'par_limrel3      ',ierr)
      call rdr(sol23_par_limrel4 ,.true.,0.0,.true.,1.0E8,
     >                                  'par_limrel4      ',ierr)
      call rdr(sol23_par_tulimit,.true.,0.0,.true.,1.0E8,
     >                                  'par_tulimit      ',ierr)
c
      call rdr(sol23_par_g0relax,.true.,0.0,.true.,1.0E8,
     >                                  'par_g0relax      ',ierr)
      call rdr(sol23_par_p0relax  ,.true.,0.0,.true.,1.0E8,
     >                                  'par_p0relax      ',ierr)
      call rdi(sol23_par_dpdrflag,.true.,0  ,.true.,10   ,
     >                                  'par_dpdrflag     ',ierr)
      call rdr(sol23_par_dpdrtemin,.true.,0.0,.true.,1.0E8,
     >                                  'par_dpdrtemin    ',ierr)
      call rdr(sol23_par_dpdrstep,.true.,0.0,.true.,1.0E8,
     >                                  'par_dpdrstep     ',ierr)
      call rdi(sol23_par_nuflag  ,.true.,0  ,.true.,10   ,
     >                                  'par_nuflag       ',ierr)
      call rdr(sol23_par_nulimit,.true.,0.0,.true.,1.0E24,
     >                                  'par_nulimit      ',ierr)
      call rdr(sol23_par_vnmult ,.true.,0.0,.true.,1.0E8,
     >                                  'par_vnmult       ',ierr)
      call rdi(sol23_par_pinmom  ,.true.,0  ,.true.,10   ,
     >                                  'par_pinmom       ',ierr)
      call rdr(sol23_par_emolec ,.true.,0.0,.true.,1.0E8,
     >                                  'par_emolec       ',ierr)
c
      call rdr(sol23_par_rec_heat,.true.,0.0,.true.,1.0E8,
     >                                  'par_rec_heat     ',ierr)
      call rdr(sol23_par_pinqimult,.true.,0.0,.true.,1.0E8,
     >                                  'par_pinqimult    ',ierr)
      call rdi(sol23_par_pinqiflag,.true.,0  ,.true.,10   ,
     >                                  'par_pinqiflag    ',ierr)
      call rdi(sol23_par_prring0,.true.,0  ,.true.,10   ,
     >                                  'par_prring0      ',ierr)
      call rdi(sol23_par_prring1 ,.true.,0  ,.true.,10   ,
     >                                  'par_prring1      ',ierr)
      call rdi(sol23_par_qperp34 ,.true.,0  ,.true.,10   ,
     >                                  'par_qperp34      ',ierr)
      call rdi(sol23_par_qeiflag ,.true.,0  ,.true.,10   ,
     >                                  'par_qeiflag      ',ierr)
      call rdi(sol23_par_chie   ,.true.,0  ,.true.,10   ,
     >                                  'par_chie         ',ierr)
      call rdi(sol23_par_joule   ,.true.,0  ,.true.,10   ,
     >                                  'par_joule        ',ierr)
      call rdi(sol23_par_fluxexp,.true.,0  ,.true.,10   ,
     >                                  'par_fluxexp      ',ierr)
      call rdi(sol23_par_qrec   ,.true.,0  ,.true.,10   ,
     >                                  'par_qrec         ',ierr)
      call rdr(sol23_par_fzrad,.true.,0.0,.true.,1.0E8,
     >                                  'par_fzrad        ',ierr)
      call rdi(sol23_par_dvmrel,.true.,1,.true.,10,
     >                                  'par_dvmrel       ',ierr)
c
      call rdi(sol23_intopt,.TRUE.,0 ,.true.,1,
     >                                  'intopt',IERR)
      call rdi(sol23_adaptgrid,.TRUE.,0 ,.true.,5,
     >                                  'adaptgrid',IERR)
      call rdi(sol23_seed,.TRUE.,1 ,.true.,2,
     >                                  'seed',IERR)
      call rdr(sol23_izlen,.true.,0.0,.false.,0.5,
     >                                  'izlen',ierr)
      call rdr(sol23_izlam,.true.,0.0,.false.,0.5,
     >                                  'izlam',ierr)
      call rdr(sol23_izoffset,.true.,0.0,.false.,0.5,
     >                                  'izoffset',ierr)
      call rdr(sol23_momlen,.true.,0.0,.false.,0.0,
     >                                  'momlen',ierr)
c
      call rdr(sol23_relax,.true.,0.0,.true.,1.0,
     >                                  'relax',ierr)
      call rdr(sol23_maxtol,.true.,0.0,.true.,1.0,
     >                             'maxtol',ierr)
      call rdr(sol23_rmstol,.true.,0.0,.true.,1.0,
     >                             'rmstol',ierr)
      call rdi(sol23_perp,.true.,0 ,.true.,10 ,
     >                                  'perp',ierr)
c
      return
      end
