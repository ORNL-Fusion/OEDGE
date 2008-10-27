c
c    CFD METHODS FOR ONION SKIN MODELING (OSM)
c
c
c      program cfd_osm
c      implicit none
c
c      call cfd_osm
c
c      return
c      end
c
c
      subroutine cfd_osm (ring,err_flag)
      implicit none
      integer ring, err_flag
c
      include 'params'
      include 'cgeom'
      include 'sol23_com'
      include 'cfd_osm_com'
      integer jj,in,ik
      real start_time, end_time, elapsed_time, test1, delta1
      real     za02as, calcwidth
      external za02as
c
      err_flag     = 0
      start_time   = za02as(1)
c
      call setparams
c
      do ik = 1, ikmax+1
         do in = 1, 4
            RE(in,ik) = 0.0D0
         enddo
      enddo
c
      write(71,*) 'STATUS:', a_error, cfd_failed(ring), finished(ring)
c
      if (ring.ge.irwall.and.(a_pin.ge.0.or.a_cfdfile.eq.1)) then
        call readring(0)
        call flatten_U
        call writering(0)
        goto 720
      endif
c
      if (a_pin.eq.0.and.a_cfdfile.eq.1) then
        call readring(0)
        call writering(0)
        goto 720
      endif
c
      if (finished(ring).eq.1) then
        call readring(0)
        call writering(0)
        goto 720
      endif
c
      if (cfd_failed(ring).eq.0) then
c
c       ACTIVE ADJUSTMENT OF SOURCES
c
        call readring(0)
  650   call solvering(2)
        call updateQperp
c
         write(71,*) 'updateQperp'
c
        if (a_Qup.eq.1.and.a_pin.ge.1
     >               .and.s23_par_updtQpit.eq.1) goto 650
c
c       GRID ADAPTATION OPTIONS
c
        if (s23_adaptgrid.eq.0) then
c
c          G.A. OFF
c
           call solvering(0)
           call print_error
           call check_error
c
        elseif (s23_adaptgrid.eq.2) then
c
c          G.A. ON
c
           if (a_error.eq.0.and.a_pin.ge.1) then
              call readring(1)
              call solvering(0)
              call print_error
              call check_error
           else
              call solvering(0)
              call print_error
              call check_error
           endif
c
        elseif (s23_adaptgrid.eq.1) then
c
c          G.A. ON for ir < irsep + ...
c
           if (a_error.eq.0.and.a_pin.ge.1.and.ring.le.irsep+
     >         s23_par_adaptnr) then
              call readring(1)
              call solvering(0)
              call print_error
              call check_error
           else
              call solvering(0)
              call print_error
              call check_error
           endif
c
        endif
c
        if (a_error.ge.1.and.a_pin.ge.1) then
           cfd_failed(ring) = 1
           finished(ring)   = 1
           write(71,*) 'CFD FAILED'
        else
           call writering(0)
        endif
c
      endif
c
      write(71,*) 'STATUS:', a_error, cfd_failed(ring), finished(ring)
c
c      write(71,*) 'x_xpt0,N', ir_cfd, x_xpt_0/x_max,
c     >                                x_xpt_N/x_max
c
c
      if (s23_par_debugflag.eq.1.and.ring.eq.irsep+s23_par_debugnr) then
        call debugring
        STOP
      endif
c
 720  end_time     = za02as(1)
      elapsed_time = end_time - start_time
c
       write(71,'(A30,G16.8)')  ' ELAPSED TIME  : ', elapsed_time
c
 999  return
      end
c
c
      subroutine readring(n_flag)
      implicit none
      integer n_flag
c
      include 'cfd_osm_com'
c
      call setparams
      call setUBC
      call setQBC
c
      if (a_pin.eq.0) then
        if (a_cfdfile.eq.1) then
           call read_cfd
        else
           call initgrid
           call initQ(n_flag)
           call initU(n_flag)
        endif
      else
        if (n_flag.eq.0) then
           call initgrid
           call initU(n_flag)
           call initQ(n_flag)
        elseif (n_flag.eq.1) then
           call writeU(1)
           call solvegrid
           call initU(n_flag)
           call initQ(n_flag)
        else
           write(71,*) 'READRING ERROR'
           STOP
        endif
      endif
c
      call WfromU
      call mkfromW
      call zeroR
      call writeU(1)
c
c      call print_U
c      call print_Q
c
      return
      end
c
c
      subroutine solvering(n_flag)
      implicit none
      integer n_flag
c
      integer ik,in,n_iter,n_temp,n_tempQ,n_loop
      real*8 c_s_0, g_net_0, htr_err,
     >        old_maxres, old_rmsres,
     >       Mach_N, Ti_N, Te_N, c_s_N
c
c      Mach_0,Te_0 and Ti_0 are declared in the SOL23 common block
c
c      real*8 Mach_0, Ti_0, Te_0, c_s_0, g_net_0, htr_err,
c     >        old_maxres, old_rmsres,
c     >       Mach_N, Ti_N, Te_N, c_s_N
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
 800  n_temp = 0
      n_tempq = 0
      n_iter = 0
c
      do n_iter = 1, itermax
c
         old_maxres = max_res
         old_rmsres = rms_res
         a_error    = 1
c
         call sweepring
         call check_error
c
         n_temp = n_temp + 1
         n_tempq = n_tempq + 1
c
         if (a_2T.eq.0) then
           c_s_0 = dsqrt(W(3,0) / U(1,0))
           c_s_N = dsqrt(W(3,ikmax+1) / U(1,ikmax+1))
         else
           c_s_0 = dsqrt((W(3,0)+W(4,0)) / U(1,0))
           c_s_N = dsqrt((W(3,ikmax+1) + W(4,ikmax+1))
     >                 / U(1,ikmax+1))
         endif
c
         if (n_tempq.ge.s23_par_refresh) then
c
           n_tempq = 0
c
           call initQ(n_flag)
c
         endif
c
         if (n_temp.ge.100.or.n_iter.eq.1) then
c
           n_temp = 0
c
           Mach_0 = abs(W(2,0)) / c_s_0
           Ti_0 = W(3,0) / (a_Temp * W(1,0) * eV)
           Te_0 = W(4,0) / (a_Temp * W(1,0) * eV)
c
           Mach_N = abs(W(2,ikmax+1)) / c_s_N
           Ti_N = W(3,ikmax+1) / (a_Temp * W(1,ikmax+1) * eV)
           Te_N = W(4,ikmax+1) / (a_Temp * W(1,ikmax+1) * eV)
c
           write(71,'(I5,8G12.6)') n_iter, max_res,
     >       rms_res, Mach_0, Mach_N, Ti_0, Te_0,
     >       Ti_N, Te_N
c
         endif
c
         if (n_flag.eq.1.and.n_iter.ge.10*ikmax_0) then
             return
         endif
c
         if (a_error.eq.0.and.n_iter.ge.100) then
c
           Mach_0 = abs(W(2,0)) / c_s_0
           Ti_0 = W(3,0) / (a_Temp * W(1,0) * eV)
           Te_0 = W(4,0) / (a_Temp * W(1,0) * eV)
c
           Mach_N = abs(W(2,ikmax+1)) / c_s_N
           Ti_N = W(3,ikmax+1) / (a_Temp * W(1,ikmax+1) * eV)
           Te_N = W(4,ikmax+1) / (a_Temp * W(1,ikmax+1) * eV)
c
           write(71,'(I5,8G12.6)') n_iter, max_res,
     >       rms_res, Mach_0, Mach_N, Ti_0, Te_0,
     >       Ti_N, Te_N
c
            return
c
         endif
c
      enddo
c
      return
      end
c
c
      subroutine writering(n_flag)
      implicit none
      integer n_flag
c
      include 'cfd_osm_com'
c
      real*8 Mach_0, T_0, c_s_0, g_net_0, c_s_N, g_net_N
c
      if (n_flag.eq.0) then
         call writeU(0)
      endif
c
      write(71,'(A30,2G16.8)') ' MAX/RMS RES.  : ', max_res, rms_res
c
      if (a_2T.eq.0) then
        c_s_0 = dsqrt(W(3,0) / U(1,0))
        g_net_0 = gamma_sum + 0.5D0 * (W(2,0)/c_s_0)**2.0D0
        c_s_N = dsqrt( W(3,ikmax+1) / U(1,ikmax+1))
        g_net_N = gamma_sum + 0.5D0 * (W(2,ikmax+1)/c_s_N)**2.0D0
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW 0:  A/C, B/C, A/B',
     >     (F(3,0) + G(3,0)) / POW_BC_0,
     >     (g_net_0 * W(3,0) * W(2,0)) / POW_BC_0,
     >     (F(3,0) + G(3,0)) / (g_net_0 * W(3,0) * W(2,0))
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW N:  A/C, B/C, A/B',
     >     (F(3,ikmax+2) + G(3,ikmax+2)) / POW_BC_N,
     >     (g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
     >      / POW_BC_N,  (F(3,ikmax+2) + G(3,ikmax+2))
     >      / (g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
      endif
c
      return
      end
c
c
      subroutine sweepring
      implicit none
c
      include 'cfd_osm_com'
      integer in
      real*8 c_s_0, g_net_0
c
      call U0fromU
      call dtfromU
      call fluxes
      call residual
      call form_tri_mtrx
      call updateU
      call updateBC
      call WfromU
      call mkfromW
      call max_res_R
c
      return
      end
c
c
      subroutine setparams
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer in,im
c
      a_PIN       = pin_iter_no
      a_cfdfile   = s23_cfdfile
      a_seed      = seedplasma
      a_symm      = s23_symm
      a_2T        = s23_2T
c
      if (a_symm.eq.0) then
         a_half   = 0
      else
         a_half = solve_half
      endif
c
      A_i         = mb
      Z_i         = cb
c
      if (a_2T.eq.0) then
         a_Temp      = 1.0D0 + Z_i
      else
         a_Temp      = 1.0D0
      endif
c
      mass_i      = A_i * mass_p
      mass        = mass_i + Z_i * mass_e
      c_v         = c_v_0 / mass
      c_p         = c_p_0 / mass
      gamma_i     = gai_const
      gamma_e     = gae_const
      gamma_sum   = (gamma_i + gamma_e) / 2.0D0
      x_max       = smax
      L_iz        = s23_izlam
      L_iz_0      = L_iz
      L_iz_N      = L_iz
      L_izoffset  = s23_izoffset
      L_mom       = s23_momlen
      L_pow       = 10.0D0 * x_max
c
      if (a_2T.eq.0) then
c
        POW_core_0  = pow_in_0
        POW_core_N  = pow_in_N
        POW_core    = POW_core_0 + POW_core_N
c
        POW_BC_0    = pow_tg_0
        POW_BC_N    = pow_tg_N
        POW_BC      = POW_BC_0 + POW_BC_N
c
      else
c
        POW_core_i0  = pow_in_0 * s23_par_ptipte / 2.0D0
        POW_core_iN  = pow_in_N * s23_par_ptipte / 2.0D0
        POW_core_i   = POW_core_i0 + POW_core_iN
c
        POW_core_e0  = pow_in_0 / 2.0D0
        POW_core_eN  = pow_in_N / 2.0D0
        POW_core_e   = POW_core_e0 + POW_core_eN
c
        POW_core_0  = POW_core_i0 + POW_core_e0
        POW_core_N  = POW_core_iN + POW_core_eN
        POW_core    = POW_core_0  + POW_core_N
c
        POW_BC_i0    = pow_tg_i0
        POW_BC_iN    = pow_tg_iN
        POW_BC_e0    = pow_tg_e0
        POW_BC_eN    = pow_tg_eN
c
        POW_BC_i     = POW_BC_i0 + POW_BC_iN
        POW_BC_e     = POW_BC_e0 + POW_BC_eN
c
        POW_BC_0     = POW_BC_i0 + POW_BC_e0
        POW_BC_N     = POW_BC_iN + POW_BC_eN
        POW_BC       = POW_BC_0  + POW_BC_N
c
      endif
c
      MASS_core_0  = mass * part_in_0
      MASS_core_N  = mass * part_in_N
      MASS_core    = MASS_core_0 + MASS_core_N
c
      MOM_core_0  = mom_in_0
      MOM_core_N  = mom_in_N
      MOM_core    = MOM_core_0 + MOM_core_N
c
      MOM_BC_0    = mom_tg_0
      MOM_BC_N    = mom_tg_N
      MOM_BC      = MOM_BC_N + MOM_BC_0
c
      MASS_BC_0    = isat_0 * mass
      MASS_BC_N    = isat_smax * mass
      MASS_BC      = MASS_BC_0 + MASS_BC_N
c
      n0_BC       = ne_0
      M0_BC       = - mach_0
      Ti0_BC      = ti_0 * eV
      Te0_BC      = te_0 * eV
c
      nN_BC       = ne_smax
      MN_BC       = mach_smax
      TiN_BC      = ti_smax * eV
      TeN_BC      = te_smax * eV
c
      a_subM      = 0
      a_error     = 1
      a_resid     = 1
      a_Q2        = 0.0D0
      a_Q3        = 0.0D0
      a_artvisc2  = s23_par_artvisc2
      a_artvisc4  = s23_par_artvisc4
      a_dt_F      = s23_par_dtf
      a_dt_G      = s23_par_dtg
      a_perp      = s23_perp
      a_grid0     = s23_par_grid0
      a_gridexp   = s23_par_gridexp
      a_flag      = 0
c
      ikmax_0     = NNN - 2
      itermax     = s23_par_itermax * ikmax_0
      max_tol     = s23_maxtol
      rms_tol     = s23_rmstol
c
      do in = 1, 4
        do im = 1, 4
          unit(in,im) = 0.0D0
          if (im.eq.in) unit(in,im) = 1.0D0
        enddo
      enddo
c
      return
      end
c
c
      subroutine setUBC
      implicit none
c
      include 'cfd_osm_com'
c
      if (a_half.eq.1) then
         M0_BC       = - abs(M0_BC)
         nN_BC       = n0_BC
         MN_BC       = - M0_BC
         TiN_BC      = Ti0_BC
         TeN_BC      = Te0_BC
         MASS_BC_N   = MASS_BC_0
c         write(71,*) 'nTM0BC',n0_BC, Ti0_BC/eV, M0_BC
      elseif (a_half.eq.2) then
         MN_BC       = abs(MN_BC)
         n0_BC       = nN_BC
         M0_BC       = - MN_BC
         Ti0_BC      = TiN_BC
         Te0_BC      = TeN_BC
         MASS_BC_0   = MASS_BC_N
c         write(71,*) 'nTMNBC',nN_BC, TiN_BC/eV, MN_BC
      endif
c
c      write(71,*) 'nTM0BC',n0_BC, Ti0_BC/eV, M0_BC
c      write(71,*) 'nTMNBC',nN_BC, TiN_BC/eV, MN_BC
c
      if (a_2T.eq.0) then
c
        T0_grid = Ti0_BC
        TN_grid = TiN_BC
c
        if (MASS_BC_0.ne.0.0D0) then
          v_0 = M0_BC * dsqrt(a_Temp*Ti0_BC/mass)
          n0_BC = abs(MASS_BC_0 / v_0) / mass
        endif
        U(1,0) = mass * n0_BC
        W(3,0) = n0_BC * Ti0_BC * a_Temp
        p_0 = W(3,0)
        v_0 = M0_BC * dsqrt( p_0 / U(1,0))
        U(3,0) = p_0 / (g_pv - 1.0D0)
     >         + 0.5D0 * mass * n0_BC * v_0 * v_0
c
        if (MASS_BC_N.ne.0.0D0) then
          v_N = MN_BC * dsqrt(a_Temp*TiN_BC/mass)
          nN_BC = abs(MASS_BC_N / v_N) / mass
        endif
        U(1,NNN) = nN_BC * mass
        W(3,NNN) = nN_BC * TiN_BC * a_Temp
        p_N = W(3,NNN)
        v_N = MN_BC * dsqrt( p_N / U(1,NNN))
        U(3,NNN) = p_N / (g_pv - 1.0D0)
     >         + 0.5D0 * mass * nN_BC * v_N * v_N
c
      else
c
        T0_grid = Te0_BC
        TN_grid = TeN_BC
c
        if (MASS_BC_0.ne.0.0D0) then
          v_0 = M0_BC * dsqrt((Ti0_BC + Te0_BC) / mass)
          n0_BC = abs(MASS_BC_0 / v_0) / mass
        endif
        U(1,0) = mass * n0_BC
        W(3,0) = n0_BC * Ti0_BC
        W(4,0) = n0_BC * Te0_BC
        v_0 = M0_BC * dsqrt( (W(3,0) + W(4,0)) / U(1,0))
        p_0 = W(3,0) + W(4,0)
        U(3,0) = W(3,0) / (g_pv - 1.0D0)
     >         + 0.5D0 * mass * n0_BC * v_0 * v_0
        U(4,0) = W(4,0) / (g_pv - 1.0D0)
c
        if (MASS_BC_N.ne.0.0D0) then
          v_N = MN_BC * dsqrt((TiN_BC + TeN_BC) / mass)
          nN_BC = abs(MASS_BC_N / v_N) / mass
        endif
        U(1,NNN) = mass * nN_BC
        W(3,NNN) = nN_BC * TiN_BC
        W(4,NNN) = nN_BC * TeN_BC
        v_N = MN_BC * dsqrt( (W(3,NNN) + W(4,NNN)) / U(1,NNN))
        p_N = W(3,NNN) + W(4,NNN)
        U(3,NNN) = W(3,NNN) / (g_pv - 1.0D0)
     >         + 0.5D0 * mass * nN_BC * v_N * v_N
        U(4,NNN) = W(4,NNN) / (g_pv - 1.0D0)
c
      endif
c
      U(2,0)   = mass * n0_BC * v_0
      U(2,NNN) = mass * nN_BC * v_N
c
c      write(71,*) 'U0', U(1,0), U(2,0), U(3,0), U(4,0)
c      write(71,*) 'UN', U(1,NNN), U(2,NNN), U(3,NNN), U(4,NNN)
c
c      write(71,*) 'nv = ', n0_BC * v_0
c      write(71,*) 'nv = ', nN_BC * v_N
c
      return
      end
c
c
      subroutine setQBC
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      real*8 src_int_ab
      external src_int_ab
c
      MASS_0 = dabs(U(2,0))
      MASS_N = dabs(U(2,NNN))
      MASS_PIN = 0.0D0
c
c      write(71,*) 'MASS_0,N', MASS_0, MASS_N, MASS_0 / mass
c
      MOM_0   = p_0 + dabs(U(1,0)   * v_0 * v_0)
      MOM_N   = p_N + dabs(U(1,NNN) * v_N * v_N)
      MOM_PIN = 0.0D0
c
c      write(71,*) 'MOM_0,N', MOM_0, MOM_N
c
      if (a_2T.eq.0) then
c
        POWi_0 = POW_core_0
        POWi_N = POW_core_N
        POWe_0 = POWi_0
        POWe_N = POWe_N
c
      else
c
        POWi_0 = POW_core_i0
        POWi_N = POW_core_iN
        POWe_0 = POW_core_e0
        POWe_N = POW_core_eN
        POW_ei = power_ei(ir_cfd)
        POW_ei_0 = POW_ei / 2.0D0
        POW_ei_N = POW_ei / 2.0D0
c
c        write(71,*) 'POW_ei', POW_ei
c
      endif
      POW_PIN  = 0.0D0
c
      Q1 = (MASS_0 + MASS_N) / x_max
      Q2 = (MOM_N  - MOM_0 ) / x_max
c
      if (a_2T.eq.0) then
        Q3 = POW_core / x_max
      else
        Q3 = POW_core_i / x_max
        Q4 = POW_core_e / x_max
      endif
c
c      write(71,*) 'Q34', Q3, Q4
c
      if (a_PIN.ge.1) then
          if (a_2T.eq.0.or.a_2T.eq.1) then
            POW_PIN_0 = (src_int_ab (ionpow_src, 0.0D0,
     >              x_max/2.0D0,spts,sbnds,npts,intopt)
     >          + src_int_ab (elecpow_src, 0.0D0,
     >              x_max/2.0D0,spts,sbnds,npts,intopt))
            POW_PIN_N = (src_int_ab (ionpow_src,
     >            x_max/2.0D0,x_max,spts,sbnds,npts,intopt)
     >          + src_int_ab (elecpow_src, x_max/2.0D0,
     >              x_max,spts,sbnds,npts,intopt))
            POW_PIN_i = src_int_ab (ionpow_src, 0.0D0,
     >                   x_max,spts,sbnds,npts,intopt)
            POW_PIN_e = src_int_ab (elecpow_src, 0.0D0,
     >                   x_max,spts,sbnds,npts,intopt)
            POW_PIN = src_int_ab (ionpow_src, 0.0D0,
     >                   x_max,spts,sbnds,npts,intopt)
     >            + src_int_ab (elecpow_src, 0.0D0,
     >                   x_max,spts,sbnds,npts,intopt)
          endif
          MOM_PIN_0 = abs(src_int_ab (mom_src, 0.0D0,
     >              x_max/2.0D0,spts,sbnds,npts,intopt))
          MOM_PIN_N = abs(src_int_ab (mom_src, x_max/2.0D0,
     >              x_max,spts,sbnds,npts,intopt))
          MOM_PIN = src_int_ab (mom_src, 0.0D0,
     >                   x_max,spts,sbnds,npts,intopt)
          MASS_PIN_0 = mass * src_int_ab (part_src, 0.0D0,
     >              x_max/2.0D0,spts,sbnds,npts,intopt)
          MASS_PIN_N = mass * src_int_ab (part_src,
     >              x_max/2.0D0,x_max,spts,sbnds,npts,intopt)
          MASS_PIN = mass * src_int_ab (part_src, 0.0D0,
     >                   x_max,spts,sbnds,npts,intopt)
          L_iz_0 = MASS_0 / abs(mass * part_src(1))
          L_iz_N = MASS_N / abs(mass * part_src(npts))
c
          if (a_half.eq.1) then
            POW_PIN_N  = POW_PIN_0
            MOM_PIN_N  = - MOM_PIN_0
            MASS_PIN_N = MASS_PIN_0
            POW_PIN    = POW_PIN_0  + POW_PIN_N
            MOM_PIN    = MOM_PIN_0  + MOM_PIN_N
            MASS_PIN   = MASS_PIN_0 + MASS_PIN_N
            L_iz = L_iz_0
          elseif (a_half.eq.2) then
            POW_PIN_0  = POW_PIN_N
            MOM_PIN_0  = - MOM_PIN_N
            MASS_PIN_0 = MASS_PIN_N
            POW_PIN    = POW_PIN_0  + POW_PIN_N
            MOM_PIN    = MOM_PIN_0  + MOM_PIN_N
            MASS_PIN   = MASS_PIN_0 + MASS_PIN_N
            L_iz = L_iz_N
          endif
c
c          POW_PIN  = POW_PIN_0  + POW_PIN_N
c          MOM_PIN  = MOM_PIN_0  + MOM_PIN_N
c          MASS_PIN = MASS_PIN_0 + MASS_PIN_N
c
c        if (POW_PIN.le.0.0D0) then
c
c            POW_RAD_0 = - abs(POW_core_0) * (1.0D0 - 1.0D0/a_maxpow_0)
c            POW_RAD_N = - abs(POW_core_N) * (1.0D0 - 1.0D0/a_maxpow_N)
c            POW_RAD   = - abs(POW_core)   * (1.0D0 - 1.0D0/a_maxpow)
c
c            POW_RAD_0 = - (POW_core_0 - POW_BC_0)
c            POW_RAD_N = - (POW_core_N - POW_BC_N)
c            POW_RAD   = - (POW_core   - POW_BC)
c
c        else
c
c            POW_RAD_0 =  min(abs(POW_PIN_0), abs(POW_core_0) *
c     >                       (1.0D0 - 1.0D0/a_maxpow_0))
c            POW_RAD_N =  min(abs(POW_PIN_N), abs(POW_core_N) *
c     >                       (1.0D0 - 1.0D0/a_maxpow_N))
c            POW_RAD   =  min(abs(POW_PIN), abs(POW_core) *
c     >                       (1.0D0 - 1.0D0/a_maxpow))
c
c            POW_RAD_0 = - (POW_core_0 - POW_BC_0)
c            POW_RAD_N = - (POW_core_N - POW_BC_N)
c            POW_RAD   = - (POW_core   - POW_BC)
c
c        endif
c
c       RADIATED POWER PARAMETRIZED
c
        if (a_2T.eq.0) then
          POW_RAD_0 = - (POW_core_0 - POW_BC_0)
          POW_RAD_N = - (POW_core_N - POW_BC_N)
          POW_RAD   = - (POW_core   - POW_BC)
        elseif (a_2T.eq.1) then
c
           POW_RAD_i = POW_PIN_i
           POW_RAD_e = POW_PIN_e
c
c          POW_RAD_i0 =   0.0D0
c          POW_RAD_iN =   0.0D0
c          POW_RAD_i  =   0.0D0
c          POW_RAD_e0 = - (POW_core_e0 - POW_BC_e0 - POW_ei_0)
c          POW_RAD_eN = - (POW_core_eN - POW_BC_eN - POW_ei_N)
c          POW_RAD_e  =    POW_RAD_e0  + POW_RAD_eN
        endif
c
c       CORE POWER PARAMETRIZED
c
c        POW_BC_0 = POW_core_0 + POW_RAD_0
c        POW_BC_N = POW_core_N + POW_RAD_N
c        POW_BC   = POW_BC_0   + POW_BC_N
c
c       TARGET POWER PARAMETRIZED
c
c        POW_BC = POW_core
c        Q3 =  (POW_core - POW_PIN) / x_max
c
        if (a_2T.eq.0) then
          T0_grid = (0.5D0 / a_Temp) * POW_BC_0 /
     >         ((gamma_sum + 0.5D0 * M0_BC ** 2.0D0) *
     >          (MASS_0 / mass))
c
          TN_grid = (0.5D0 / a_Temp) * POW_BC_N /
     >         ((gamma_sum + 0.5D0 * MN_BC ** 2.0D0) *
     >          (MASS_N / mass))
        elseif (a_2T.eq.1) then
          T0_grid = POW_BC_e0 / (gamma_e * (MASS_0 / mass))
          TN_grid = POW_BC_eN / (gamma_e * (MASS_N / mass))
        endif
c
c      write(71,*) 'MASS_PIN', MASS_PIN_0, MASS_PIN_N
c      write(71,*) 'MOM_PIN',  MOM_PIN_0, MOM_PIN_N
c      write(71,*) 'MOM_core', MOM_core_0, MOM_core_N
c      write(71,*) 'MOM_BC',   MOM_BC_0, MOM_BC_N
c      write(71,*) 'POW_PIN',  POW_PIN_0, POW_PIN_N
c      write(71,*) 'POW_RAD',  POW_RAD_0, POW_RAD_N
c      write(71,*) 'POW_core', POW_core_0, POW_core_N
c      write(71,*) 'POW_BC',   POW_BC_0, POW_BC_N
c      write(71,*) 'POW_core/_BC', POW_core_0 / POW_BC_0,
c     >                            POW_core_N / POW_BC_N
c      write(71,*) 'POW_RAD/_core', POW_RAD_0 / POW_core_0,
c     >                             POW_RAD_N / POW_core_N
c
c        v_0 = M0_BC * dsqrt(a_Temp * T0_grid / mass)
c        U(1,0) = U(2,0) / v_0
c        p_0 =  (U(1,0)/mass) * T0_grid * a_Temp
c        MOM_0 = p_0 + dabs(U(1,0)   * v_0 * v_0)
c        U(3,0) = p_0 / (g_pv - 1.0D0)
c     >         + 0.5D0 * U(1,0) * v_0 * v_0
c        U(1,NNN) = U(1,0)
c        MOM_N = MOM_0
c        U(3,NNN) = U(3,0)
c
      endif
c
      return
      end
c
c
      subroutine grid_error
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer j, in, ik
      real*8 x_a, src_val, temp1, temp2, temp3, f_max, f_min,
     >    dt_grid_0
      external src_val
c
      j = ikmax_cfd(ir_cfd, targ_cfd)
      bbb(0)   = 0.0
      bbb(j+2) = 0.0
c
      do ik = 0, j+2
           uuu(ik) = xc_cfd(ir_cfd, targ_cfd, ik)
           fff(ik) = xf_cfd(ir_cfd, targ_cfd, ik)
c
           if (ik.ne.0.and.ik.ne.j+2) then
             aaa(ik) = abs( (U(1,ik+1) - U(1,ik-1))
     >                       / U(1,0) ) ** 3.0
             bbb(ik) = abs( (U(2,ik+1) - U(2,ik-1))
     >                       / U(2,0) )
             ccc(ik) = abs( (U(3,ik+1) - U(3,ik-1))
     >                       / U(3,0) )
             ddd(ik) = abs( (U(4,ik+1) - U(4,ik-1))
     >                       / U(4,0) )
           endif
c
           if (a_2T.eq.0) then
             ccc(ik) = (aaa(ik) * bbb(ik) * ccc(ik)) ** 0.2
           else
             ccc(ik) = (aaa(ik) * bbb(ik) * ddd(ik)) ** 0.2
           endif
      enddo
c
c      write(0,*) 'OK1'
      f_max = 0.0D0
      f_min = 1.0D99
      do ik = 1, ikmax
         x_a = xc(ik)
         do in = 1, 3
            R(1,ik) = abs(src_val(ccc,x_a,uuu,fff,j+2,intopt))
            if (R(1,ik).ge.f_max) f_max = R(1,ik)
            if (R(1,ik).le.f_min) f_min = R(1,ik)
         enddo
      enddo
c      write(0,*) 'OK2'
c
      R(1,0)       = f_min
      R(1,ikmax+1) = f_min
      do ik = 0, ikmax+1
         W(1,ik) = 1.0D0 + s23_par_ga1 * (R(1,ik) - f_min) /
     >                              (f_max - f_min)
      enddo
c      write(0,*) 'OK3'
c
      do ik = 1, ikmax
         temp1 = (W(1,ik-1) + W(1,ik) + W(1,ik+1)) / 3.0D0
         if (W(1,ik).ge.2.0D0*temp1) then
              W(2,ik) = temp1
         else
              W(2,ik) = W(1,ik)
         endif
      enddo
c      write(0,*) 'OK3'
c
      do ik = 1, ikmax
         W(1,ik) = W(2,ik)
c         write(71,*) ik, W(1,ik)
      enddo
c
      do ik = 1, ikmax
         W(2,ik) = (W(1,ik-1) + 2.0D0*W(1,ik) + W(1,ik+1))/4.0D0
      enddo
c      write(0,*) 'OK4'
c
      do ik = 1, ikmax
         W(1,ik) = W(2,ik)
c         write(71,*) ik, W(1,ik)
      enddo
c
      do ik = 1, ikmax
         R(1,ik) = (W(1,ik+1) - W(1,ik-1)) / (2.0D0 * W(1,ik))
c         write(71,*) ik, R(1,ik)
      enddo
      R(1,0)       = 0.0D0
      R(1,ikmid)   = 0.0D0
      R(1,ikmax+1) = 0.0D0
c      write(0,*) 'OK5'
c
      R(2,0) = R(2,0) * (1.0D0 - s23_par_ga2) -
     >         s23_par_ga2 * (xc(1) - dx0_BC) / dx0_BC
      R(2,ikmax+1) = R(2,ikmax+1) * (1.0D0 - s23_par_ga2) +
     >         s23_par_ga2 * ((xc(ikmax+1) - xc(ikmax)) - dxN_BC)
     >                   / dxN_BC
c
c      R(2,ikmax+1) = - R(2,0)
c
c      write(71,*) 'R20', R(2,0), xc(1), dx0_BC
c
      do ik = 1, ikmid
         R(2,ik)         = R(2,0)       * exp(-s23_par_ga3*(ik-1))
         R(2,ikmax+1-ik) = R(2,ikmax+1) * exp(-s23_par_ga3*(ik-1))
c
c         R(2,ikmax+1-ik) = - R(2,ik)
c
      enddo
c
      if (a_symm.eq.1) then
        R(2,ikmid) = 0.0D0
      endif
c
c      write(0,*) 'OK6'
c
c      write(71,*) 'fcontrol'
      do ik = 1, ikmax
        R(3,ik) = 1.0D0 * R(1,ik) + s23_par_ga4 * R(2,ik)
c        write(71,'(I5,3G16.8)') ik, R(1,ik), R(2,ik), R(3,ik)
      enddo
c
      dt_grid_0 = 1.0D0
c      write(0,*) 'OK7'
c
      do ik = 1, ikmax
           aaa(ik) = (1.0D0 - 0.5D0 * R(3,ik)) * dt_grid_0
           bbb(ik) = - 1.0D0 - 2.0D0 * dt_grid_0
           ccc(ik) = (1.0D0 + 0.5D0 * R(3,ik)) * dt_grid_0
           fff(ik) = - xc(ik)
      enddo
      fff(1)     = fff(1)     - 0.0D0 * aaa(1) * dt_grid_0
      fff(ikmax) = fff(ikmax) - x_max * ccc(ikmax) * dt_grid_0
c
c      write(0,*) 'OK8'
c
c      do ik = 1, ikmax
c         write(71,*) ik, aaa(ik), bbb(ik), ccc(ik), fff(ik)
c      enddo
c
      call tridiag(aaa,bbb,ccc,fff,uuu,ikmax)
c
c      write(0,*) 'OK9'
c
      do ik = 1, ikmax
c         write(71,'(I5,3G16.8)') ik, uuu(ik), uuu(ik) - xc(ik),
c     >                   uuu(ik) - uuu(ik-1)
         xc(ik) = uuu(ik)
      enddo
c
      do ik = 0, ikmid
         xc(ikmax+1-ik) = 2.0D0 * xc(ikmid) - xc(ik)
      enddo
c
      do ik = 1, ikmax+1
         xc(ik) = x_max * (xc(ik) / xc(ikmax+1))
      enddo
c
c      write(71,'(I5,2G16.8)') ikmax+1, x_max, x_max - uuu(ikmax)
c
c      write(0,*) 'GRID OK'
c      STOP
c
      return
      end
c
c
      real*8 function f_grid(ik)
      implicit none
c
      integer ik
      include 'cfd_osm_com'
      real*8 temp1
c
      temp1 = a_gridexp
c      if (xc(ik-1).ge.L_izoffset.and.
c     >    xc(ik-1).le.L_izoffset+L_iz) then
c         write(0,*) '(---)', ik, xc(ik-1)
c         temp1 = 1.0D0
c      endif
      f_grid =  exp(log(temp1) * ik)
c
c      ik_off = max(0.0D0, xc(ik) - L_izoffset)
c      ik_off = abs(ik - 3)
c      f_grid = exp(log(a_gridexp) * ik * ik) * ik ** 5.0D0
c
      return
      end
c
c
      real*8 function d1dx1( A, B, dx_ab)
      implicit none
      real*8 A, B, dx_ab
c
c     d/dx based on two points
c     2nd order accurate on an uneven grid
c
      if (dx_ab.eq.0.0D0) then
         write(71,*) 'ZERO INTERVAL dx in d1dx1'
         STOP
      endif
c
      d1dx1 = (B - A) / dx_ab
c
      return
      end
c
c
      real*8 function d2dx2( A, B, C, dx_ab, dx_bc)
      implicit none
      real*8 A, B, C, dx_ab, dx_bc
c
c     d2/dx2 based on three points
c     2nd order accurate on an uneven grid
c
      if (dx_ab.eq.0.0D0.or.dx_bc.eq.0.0D0) then
         write(71,*) 'ZERO INTERVAL dx in D2DX2'
         STOP
      endif
c
      d2dx2 = ((C - B)/dx_bc - (B - A)/dx_ab) *
     >         2.0D0 / (dx_ab + dx_bc)
c
      return
      end
c
c
      real*8 function d3dx3( A, B, C, D, dx_ab, dx_bc, dx_cd)
      implicit none
      real*8 A, B, C, D, dx_ab, dx_bc, dx_cd
      real*8 divB, divC, d2dx2, d1dx1
      external d2dx2, d1dx1
c
c     d3/dx3 based on four points
c     2nd order accurate on an uneven grid
c
      if (dx_ab.eq.0.0D0.or.dx_bc.eq.0.0D0.or.
     >                      dx_cd.eq.0.0D0) then
         write(71,*) 'ZERO INTERVAL dx in D3DX3'
         STOP
      endif
c
      divB = d2dx2( A, B, C, dx_ab, dx_bc)
      divC = d2dx2( B, C, D, dx_bc, dx_cd)
      d3dx3 = d1dx1( divB, divC, dx_bc)
c
      return
      end
c
c
      real*8 function interpol2p ( A, B, xa, xb, xi)
      implicit none
      real*8 A, B, xa, xb, xi
c
c     two point linear interpolation
c     2nd order on an uneven grid
c
      interpol2p = A + (B - A)*(xi - xa) / (xb - xa)
c
      return
      end
c
c
      real*8 function extrapol2p ( A, B, xa, xb, xi)
      implicit none
      real*8 A, B, xa, xb, xi
c
c     two point linear extrapolation
c     2nd order on an uneven grid
c
      extrapol2p = A + (B - A)*(xi - xa) / (xb - xa)
c
      return
      end
c
c
      real*8 function ave2p ( A, B, dx_ab, dx_bc)
      implicit none
      real*8 A, B, dx_ab, dx_bc
c
c     two point weighted average
c     2nd order on an uneven grid
c
      if (dx_ab.eq.0.0D0.or.dx_bc.eq.0.0D0) then
         write(71,*) 'ZERO INTERVAL dx in AVE2P'
         STOP
      endif
c
      ave2p = (A * dx_bc + B * dx_ab) / (dx_ab + dx_bc)
c
      return
      end
c
c
      real*8 function ave3p ( A, B, C, dx_ab, dx_bc)
      implicit none
      real*8 A, B, C, dx_ab, dx_bc
c
c     three point weighted average
c     2nd order on an uneven grid
c
      if (dx_ab.eq.0.0D0.or.dx_bc.eq.0.0D0) then
         write(71,*) 'ZERO INTERVAL dx in AVE3P'
         STOP
      endif
c
      ave3p = ((C + B)/dx_bc + (B + A)/dx_ab) *
     >         2.0D0 / (dx_ab + dx_bc)
c
      return
      end
c
c
      subroutine tridiag(a,b,c,f,u,n)
      implicit none
c
c     solves tridiag linear system, see p.43
c     Numerical recipes
c
      integer n,nmax
      real*8 a(0:n),b(0:n),c(0:n),f(0:n),u(0:n)
      parameter (nmax=1000)
c
      integer j
      real*8 beta(0:nmax),gamma(0:nmax)
c
c      do j = 0, n
c        write(71,*) 'abcf',j,a(j),b(j),c(j),f(j)
c      enddo
c      write(71,*)
c
      if (b(1).eq.0.0D0) then
         write(71,*) 'ERROR in tridaig,beta=0'
         STOP
      endif
      beta(1) = b(1)
      gamma(1) = f(1) / beta(1)
      do j = 2, n
        beta(j) = b(j) - (a(j) * c(j-1)) / beta(j-1)
        gamma(j) = (f(j) - a(j) * gamma(j-1)) / beta(j)
        if (beta(j).eq.0.0D0) then
            write(71,*) 'ERROR in tridaig,beta=0'
            STOP
        endif
      enddo
c
      u(n) = gamma(n)
      do j = n-1, 1, -1
         u(j) = gamma(j) - u(j+1) * c(j) / beta(j)
      enddo
c
c      do j = 1, n
c        write(71,*) 'bgu',j,beta(j),gamma(j),u(j)/1.6D-19
c      enddo
c      write(71,*)
c      STOP
c
      return
      end
c
c
      subroutine ludcmp(a,n,np,indx,d)
      implicit none
c
      integer n, np ,indx(n), nmax
      real*8 d, a(np,np), tiny
      parameter (nmax=500, tiny=1.0D-20)
      integer i,imax,j,k
      real*8 aamax, dum, sum, vv(nmax)
c
      d = 1.0D0
      do i = 1, n
        aamax = 0.0D0
        do j = 1, n
           if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
        enddo
        if (aamax.eq.0.0D0)  then
          write(71,*) 'AAMAX = 0 in LUDCMP'
c          write(0,*) 'AAMAX = 0 in LUDCMP'
          STOP
c
c          call print_R
c          call print_Q
c          call print_U
c          call print_W
c          call print_symmetry
c          STOP
c
        endif
        vv(i) = 1.0D0 / aamax
      enddo
c
      do j = 1, n
        do i = 1, j-1
          sum = a(i,j)
          do k = 1, i-1
             sum = sum - a(i,k) * a(k,j)
          enddo
          a(i,j) = sum
        enddo
        aamax = 0.0D0
        do i = j, n
          sum = a(i,j)
          do k = 1, j-1
             sum = sum - a(i,k) * a(k,j)
          enddo
          a(i,j) = sum
          dum = vv(i) * dabs(sum)
          if (dum.ge.aamax) then
             imax = i
             aamax = dum
          endif
        enddo
        if (j.ne.imax) then
          do k = 1, n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          enddo
          d = -d
          vv(imax) = vv(j)
        endif
        indx(j) = imax
        if (a(j,j).eq.0.0D0) a(j,j) = tiny
        if (j.ne.n) then
           dum = 1.0D0 / a(j,j)
           do i = j+1, n
              a(i,j) = a(i,j) * dum
           enddo
        endif
      enddo
c
      return
      end
c
c
      subroutine lubksb(a,n,np,indx,b)
      implicit none
c
      integer n, np, indx(n)
      real*8 a(np,np), b(n)
      integer i,ii,j,ll
      real*8 sum
c
      ii = 0
      do i = 1, n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if (ii.ne.0) then
           do j = ii, i-1
             sum = sum - a(i,j) * b(j)
           enddo
        else if (sum.ne.0.0D0) then
           ii = i
        endif
        b(i) = sum
      enddo
      do i = n, 1, -1
         sum = b(i)
         do j = i+1, n
            sum = sum - a(i,j) * b(j)
         enddo
         b(i) = sum / a(i,i)
      enddo
c
      return
      end
c
c
      subroutine blocktridiag(A,B,C,F,X,N,M)
      implicit none
c
      integer N,M
      real*8 A(M,M,0:N), B(M,M,0:N), C(M,M,0:N), beta(M,M,0:N),
     >       betainv(M,M,0:N), F(M,0:N), X(M,0:N), gamma(M,0:N)
      integer k
c
      call findbeta(1,A,B,C,beta,betainv,N,M)
      call invertbeta(1,beta,betainv,N,M)
      call findgamma(1,A,F,betainv,gamma,N,M)
c
      do k = 2, N
        call findbeta(k,A,B,C,beta,betainv,N,M)
        call invertbeta(k,beta,betainv,N,M)
        call findgamma(k,A,F,betainv,gamma,N,M)
      enddo
c
      call findX(N,C,betainv,gamma,X,N,M)
      do k = N-1, 1, -1
        call findX(k,C,betainv,gamma,X,N,M)
      enddo
c
      return
      end
c
c
      subroutine findbeta(k,A,B,C,beta,betainv,N,M)
      implicit none
c
      integer k,N,M
      real*8 A(M,M,0:N), B(M,M,0:N), C(M,M,0:N), beta(M,M,0:N),
     >       betainv(M,M,0:N)
      integer i,j,l
      real*8 sum, Y(M,M)
c
      if (k.eq.1) then
        do i = 1, M
          do j = 1, M
            beta(i,j,k) = B(i,j,k)
          enddo
        enddo
      else
        do i = 1, M
          do j = 1, M
            sum = 0.0D0
            do l = 1, M
              sum = sum + betainv(i,l,k-1) * C(l,j,k-1)
            enddo
            Y(i,j) = sum
          enddo
        enddo
c
        do i = 1, M
          do j = 1, M
            sum = 0.0D0
            do l = 1, M
              sum = sum + A(i,l,k) * Y(l,j)
            enddo
            beta(i,j,k) = B(i,j,k) - sum
          enddo
        enddo
      endif
c
      return
      end
c
c
      subroutine invertbeta(k,beta,betainv,N,M)
      implicit none
c
      integer k,N,M
      real*8 beta(M,M,0:N), betainv(M,M,0:N)
      integer i,j,l, indx(M)
      real*8 sum, Y(M,M), one(M,M), dperm
c
      do i = 1, M
        do j = 1, M
          Y(i,j)   = beta(i,j,k)
          one(i,j) = 0.0D0
        enddo
        one(i,i) = 1.0D0
      enddo
c
      call ludcmp(Y,M,M,indx,dperm)
c
      do i = 1, M
        call lubksb(Y,M,M,indx,one(1,i))
      enddo
c
      do i = 1, M
        do j = 1, M
          betainv(i,j,k) = one(i,j)
        enddo
c        write(71,*) one(i,1), one(i,2), one(i,3)
      enddo
c      write(71,*)
c
        do i = 1, M
          do j = 1, M
            sum = 0.0D0
            do l = 1, M
              sum = sum +
     >              beta(i,l,k) * betainv(l,j,k)
            enddo
            Y(i,j) = sum
          enddo
c          write(71,*) Y(i,1), Y(i,2), Y(i,3)
        enddo
c      write(71,*)
c
      return
      end
c
c
      subroutine findgamma(k,A,F,betainv,gamma,N,M)
      implicit none
c
      integer k,N,M
      real*8 A(M,M,0:N), betainv(M,M,0:N), F(M,0:N), gamma(M,0:N)
      integer i,j,l
      real*8 sum, Y(M,M), Z(M)
c
      if (k.eq.1) then
        do i = 1, M
          sum = 0.0D0
          do j = 1, M
            sum = sum + betainv(i,j,k) * F(j,k)
          enddo
          gamma(i,k) = sum
        enddo
      else
        do i = 1, M
          sum = 0.0D0
          do j = 1, M
            sum = sum + A(i,j,k) * gamma(j,k-1)
          enddo
          Z(i) = F(i,k) - sum
        enddo
c
        do i = 1, M
          sum = 0.0D0
          do j = 1, M
            sum = sum + betainv(i,j,k) * Z(j)
          enddo
          gamma(i,k) = sum
        enddo
      endif
c
      return
      end
c
c
      subroutine findX(k,C,betainv,gamma,X,N,M)
      implicit none
c
      integer k,N,M
      real*8 C(M,M,0:N), betainv(M,M,0:N), gamma(M,0:N), X(M,0:N)
      integer i,j,l
      real*8 sum, Z(M)
c
      if (k.eq.N) then
        do i = 1, M
          X(i,k) = gamma(i,k)
        enddo
      else
        do i = 1, M
          sum = 0.0D0
          do j = 1, M
            sum = sum + C(i,j,k) * X(j,k+1)
          enddo
          Z(i) = sum
        enddo
c
        do i = 1, M
          sum = 0.0D0
          do j = 1, M
            sum = sum + betainv(i,j,k) * Z(j)
          enddo
          X(i,k) = gamma(i,k) - sum
        enddo
      endif
c
      return
      end
c
c
c      subroutine NBTRIP(A,B,C,D,IL,IU,ORDER)
c      implicit none
c
c      DIMENSION A(1), B(1), C(1), D(1)
c      real*8 A, B, C, D, SUM
c      integer IL, IU, ORDER, ORDSQ, I, J, IOMAT, IOVEC,
c     >        IOMATJ, I1MAT, I1MATJ, I1VEC
c
c      ORDSQ = ORDER ** 2
c
c      I = IL
c      IOMAT = 1 + (I-1) * ORDSQ
c      IOVEC = 1 + (I-1) * ORDER
c      CALL LUDECO( B(IOMAT), ORDER )
c      CALL LUSOLV( B(IOMAT), D(IOVEC), D(IOVEC), ORDER )
c      DO J = 1, ORDER
c        IOMATJ = IOMAT + (J-1) * ORDER
c        CALL LUSOLV( B(IOMAT), C(IOMATJ), C(IOMATJ), ORDER )
c      ENDDO
c
c 200  I = I + 1
c      IOMAT = 1 + (I-1) * ORDSQ
c      IOVEC = 1 + (I-1) * ORDER
c      I1MAT = IOMAT - ORDSQ
c      I1VEC = IOVEC - ORDER
c      CALL MULPUT( A(IOMAT), D(I1VEC), D(IOVEC), ORDER )
c      DO J = 1, ORDER
c        IOMATJ = IOMAT + (J-1) * ORDER
c        I1MATJ = I1MAT + (J-1) * ORDER
c        CALL MULPUT( A(IOMAT), C(I1MATJ), B(IOMATJ), ORDER )
c      ENDDO
c      CALL LUDECO( B(IOMAT), ORDER )
c      CALL LUSOLV( B(IOMAT), D(IOVEC), D(IOVEC), ORDER )
c      IF (I.EQ.IU) GOTO 500
c      DO J = 1, ORDER
c        IOMATJ = IOMAT + (J-1) * ORDER
c        CALL LUSOLV( B(IOMAT), C(IOMATJ), C(IOMATJ), ORDER )
c      ENDDO
c      GOTO 200
c 500  CONTINUE
c
c      I = IU
c 600  I = I - 1
c      IOMAT = 1 + (I-1) * ORDSQ
c      IOVEC = 1 + (I-1) * ORDER
c      I1VEC = IOVEC + ORDER
c      CALL MULPUT( C(IOMAT), D(I1VEC), D(IOVEC), ORDER )
c      IF (I.GT.IL) GOTO 600
c
c      return
c      end
c
c
c      subroutine LUDECO( A, ORDER )
c      implicit none
c
c      dimension A(ORDER, 1)
c      real*8  A, SUM
c      integer ORDER, JC, JR, JM, JRJC, JRJCM1, JRJCP1
c
c      DO JC = 2, ORDER
c        A(1,JC) = A(1,JC) / A(1,1)
c      ENDDO
c      JRJC = 1
c
c 10   JRJC   = JRJC + 1
c      JRJCM1 = JRJC - 1
c      JRJCP1 = JRJC + 1
c      DO JR = JRJC, ORDER
c        SUM = A(JR, JRJC)
c        DO JM = 1, JRJCM1
c          SUM = SUM - A(JR, JM) * A(JM, JRJC)
c        ENDDO
c        A(JR, JRJC) = SUM
c      ENDDO
c      IF (JRJC.EQ.ORDER) RETURN
c      DO JC = JRJCP1, ORDER
c        SUM = A(JRJC, JC)
c        DO JM = 1, JRJCM1
c          SUM = SUM - A(JRJC, JM) * A(JM, JC)
c        ENDDO
c        A(JRJC, JC) = SUM / A(JRJC, JRJC)
c      ENDDO
c      GOTO 10
c
c      return
c      end
c
c
c      subroutine LUSOLV(A,B,C,ORDER)
c      implicit none
c
c      dimension A(ORDER,1), B(1), C(1)
c      real*8  A, B, C, SUM
c      integer ORDER, JR, JM, JRM1, JRP1, JMJM, JRJR
c
c      C(1) = C(1) / A(1,1)
c      DO JR = 2, ORDER
c        JRM1 = JR - 1
c        SUM = B(JR)
c        DO JM = 1, JRM1
c          SUM = SUM - A(JR, JM) * C(JM)
c        ENDDO
c        C(JR) = SUM / A(JR,JR)
c      ENDDO
c
c      DO JRJR = 2, ORDER
c        JR = ORDER - JRJR + 1
c        JRP1 = JR + 1
c        SUM = C(JR)
c        DO JMJM = JRP1, ORDER
c          JM = ORDER - JMJM + JRP1
c          SUM = SUM - A(JR,JM) * C(JM)
c        ENDDO
c        C(JR) = SUM
c      ENDDO
c
c      return
c      end
c
c
c      subroutine MULPUT(A,B,C,ORDER)
c      implicit none
c
c      dimension A(1), B(1), C(1)
c      real*8 A, B, C, SUM
c      integer ORDER, JR, JC, IA, ik
c
c      do ik = 1, 3
c        write(0,'(5G16.8)') A(ik), B(ik), C(ik)
c      enddo
c
c      DO JR = 1, ORDER
c        SUM = 0.0D0
c        DO JC = 1, ORDER
c          IA = JR + (JC-1) * ORDER
c          SUM = SUM + A(IA) * B(JC)
c        ENDDO
c        C(JR) = C(JR) - SUM
c      ENDDO
c
c      return
c      end
c
c
c      subroutine MULPUT_0(A,B,C,ORDER)
c      implicit none
c
c      real*8 A,B,C
c      dimension A(1)
c      integer order
c
c      call MULPUT(A,B,C,ORDER)
c
c      return
c      end
c
c
      subroutine debugring
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
      real*8  Mach_n0, T_0
c
      write(71,*) 'DEBUG:'
c
      call print_U
      call print_W
      call print_Q
      call print_R
c
      call print_mass
      call print_pressure
      call print_power
c
c      call print_symmetry
c
      return
      end
c
c
      subroutine U0fromU
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
c
      do ik = 0, ikmax + 1
         if (U(1,ik).le.0.0D0.or.U(3,ik).le.0.0D0) then
           write (71,*) 'U1,U3<0', ik, U(1,ik), U(3,ik), U(4,ik)
         endif
c
         U0(1,ik) = dabs(U(1,ik))
         U0(2,ik) = U(2,ik)
         U0(3,ik) = dabs(U(3,ik))
         U0(4,ik) = dabs(U(4,ik))
      enddo
c
      return
      end
c
c
      subroutine Q0fromQ
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
c
      do ik = 0, ikmax + 1
         do in = 1, 4
             Q0(in,ik) = Q(in,ik)
         enddo
      enddo
c
      return
      end
c
c
      subroutine WfromU
      implicit none
c
c     n = rho / m
c     v = qm / rho
c     p = 2/3 * (Ei - 1/2 * m * n * v^2)
c
      include 'cfd_osm_com'
      integer ik
      real*8 vel,interpol2p, dx_ab, dx_bc, netv, d2dx2,
     >      p_a, p_b, p_c, p_div, p_ave, ave3p
c
      W(2,0) = U(2,0) / U(1,0)
      W(2,ikmax+1) = U(2,ikmax+1) / U(1,ikmax+1)
c
      do ik = 1, ikmax
         W(1,ik) = U(1,ik) / mass
         W(2,ik) = U(2,ik) / U(1,ik)
         W(3,ik) = (g_pv - 1.0D0) *
     >        (U(3,ik) - 0.5D0 * U(1,ik)* W(2,ik) **2.0D0)
         W(4,ik) = (g_pv - 1.0D0) * U(4,ik)
         W(1,ik) = dabs(W(1,ik))
         W(3,ik) = dabs(W(3,ik))
         W(4,ik) = dabs(W(4,ik))
         if (a_2T.eq.0) then
           c_s_j(ik) = dsqrt(g_pv * W(3,ik) / U(1,ik))
           W(4,ik) = W(3,ik)
         else
           c_s_j(ik) = dsqrt(g_pv * (W(3,ik) + W(4,ik)) / U(1,ik))
         endif
      enddo
c
      W(1,0) = U(1,0) / mass
      W(2,0) = U(2,0) / U(1,0)
      W(3,0) = (g_pv - 1.0D0) * (U(3,0) - 0.5D0*U(1,0)*W(2,0)**2.0D0)
      W(4,0) = (g_pv - 1.0D0) *  U(4,0)
      W(1,0) = dabs(W(1,0))
      W(3,0) = dabs(W(3,0))
      W(4,0) = dabs(W(4,0))
      if (a_2T.eq.0) then
        c_s_j(0) = dsqrt(g_pv * W(3,0) / U(1,0))
        W(4,0) = W(3,0)
      else
        c_s_j(0) = dsqrt(g_pv * (W(3,0) + W(4,0)) / U(1,0))
      endif
c
      W(1,ikmax+1) = U(1,ikmax+1) / mass
      W(2,ikmax+1) = U(2,ikmax+1) / U(1,ikmax+1)
      W(3,ikmax+1) = (g_pv - 1.0D0) * (U(3,ikmax+1) -
     >       0.5D0*U(1,ikmax+1)* W(2,ikmax+1) **2.0D0)
      W(4,ikmax+1) = (g_pv - 1.0D0) * U(4,ikmax+1)
      W(1,ikmax+1) = dabs(W(1,ikmax+1))
      W(3,ikmax+1) = dabs(W(3,ikmax+1))
      W(4,ikmax+1) = dabs(W(4,ikmax+1))
      if (a_2T.eq.0) then
        c_s_j(ikmax+1) = dsqrt(g_pv * W(3,ikmax+1) / U(1,ikmax+1))
        W(4,ikmax+1) = W(3,ikmax+1)
      else
        c_s_j(ikmax+1) = dsqrt(g_pv * (W(3,ikmax+1) + W(4,ikmax+1))
     >                    / U(1,ikmax+1))
      endif
c
      do ik = 1, ikmax
         dx_ab = xc(ik) - xc(ik-1)
         dx_bc = xc(ik+1) - xc(ik)
         if (a_2T.eq.0) then
           p_a = W(3,ik-1)
           p_b = W(3,ik)
           p_c = W(3,ik+1)
         else
           p_a = W(3,ik-1) + W(4,ik-1)
           p_b = W(3,ik)   + W(4,ik)
           p_c = W(3,ik+1) + W(4,ik+1)
         endif
         p_div = dabs( d2dx2(p_a, p_b, p_c, dx_ab, dx_bc) )
         p_ave = dabs( ave3p(p_a, p_b, p_c, dx_ab, dx_bc) )
         psi_p(ik) = p_div / p_ave
      enddo
      psi_p(0) = psi_p(1)
      psi_p(ikmax+1) = psi_p(ikmax)
c
      do ik = 1, ikmax
         netv = c_s_j(ik) + dabs(W(2,ik))
         ep2(ik) = a_artvisc2 * netv *
     >             max(psi_p(ik-1), psi_p(ik), psi_p(ik+1))
         ep4(ik) = max(0.0D0, a_artvisc4 * netv - ep2(ik))
      enddo
      ep2(0) = ep2(1)
      ep2(ikmax+1) = ep2(ikmax)
      ep4(0) = ep4(1)
      ep4(ikmax+1) = ep4(ikmax)
c
      return
      end
c
c
      subroutine initU(n_flag)
      implicit none
c
      include 'cfd_osm_com'
c
      integer ik,in,j,n_flag
      real*8  source,flux,sum,ratio,x_c,vel,press,
     >        interpol2p
c
      U(1,0) = U(1,0)
      U(2,0) = U(2,0)
      U(3,0) = U(3,0)
      U(4,0) = U(4,0)
c
      U(1,ikmax+1) = U(1,NNN)
      U(2,ikmax+1) = U(2,NNN)
      U(3,ikmax+1) = U(3,NNN)
      U(4,ikmax+1) = U(4,NNN)
c
      if (a_PIN.ge.1) then
c
        call readU(n_flag)
c
      else
        source = 0.0D0
        do ik = 1, ikmax
           source = source + 0.5D0 * (Q(1,ik) + Q(1,ik-1))
     >                             * (xc(ik) - xc(ik-1))
c
           U(1,ik) = U(1,0) +
     >       (U(1,ikmax+1) - U(1,0)) * xc(ik) / xc(ikmax+1)
           U(2,ik) = U(2,0) + source
           U(4,ik) = U(4,0) +
     >       (U(4,ikmax+1) - U(4,0)) * xc(ik) / xc(ikmax+1)
           U(3,ik) = U(3,0) +
     >       (U(3,ikmax+1) - U(3,0)) * xc(ik) / xc(ikmax+1)
c
        enddo
      endif
c
c      U(1,0) = U(1,0)
c      U(2,0) = U(2,0)
c      U(3,0) = U(3,0)
c
c      U(1,ikmax+1) = U(1,NNN)
c      U(2,ikmax+1) = U(2,NNN)
c      U(3,ikmax+1) = U(3,NNN)
c
c
c      call updateBC
c
c      Ti_ev(0) = Ti0_BC / eV
c      Ti_eV(ikmax+1) = TiN_BC / eV
c      Te_eV(0) = Ti_eV(0)
c      Te_eV(ikmax+1) = Ti_eV(ikmax+1)
c
c      do ik = 1, ikmax
c         x_c = xc(ik)
c         flux = U(2,ik)
c         ratio = dabs(flux / U(2,0))
c         source = x_c * ((POWi_0) - (Q3) * x_c / 2.0D0)
c         sum = (TI0_BC/eV) ** 3.5D0 + (3.5D0 / ka0_ev) * source
c         Ti_eV(ik) = (sum ** (2.0D0/7.0D0)) * (1.0D0 - ratio) +
c     >          (TI0_BC/eV) * ratio
c         Te_eV(ik) = Ti_eV(ik)
c         ratio = a_Temp * (Ti_eV(ik)) * eV / mass
c         source = MOM_0 * MOM_0 - 4.0D0 * ratio * flux ** 2.0D0
c         if (source.le.0.0D0) then
c             source = 0.0D0
c             write(71,*) 'PRESSURE TOO LOW at ik,x ',ik,xc(ik)
c         endif
c         U(1,ik) = (MOM_0 + dsqrt(source)) / (2.0D0 * ratio)
c         vel = flux / U(1,ik)
c         press = U(1,ik) * ratio
c         U(3,ik) = press / (g_pv - 1.0D0)
c     >            + 0.5D0 * U(1,ik) * vel * vel
c      enddo
c
c      source = 0.0D0
c      do ik = 1, ikmax+1
c         flux = U(2,ik)
c         ratio = 0.5D0 * (Ti_eV(ik-1) + Ti_eV(ik)) * eV / mass
c         source = MOM_0 * MOM_0 - 4.0D0 * ratio * flux ** 2.0D0
c         W(2,ik) = (MOM_0 - dsqrt(source)) / (2.0D0 * flux)
c         write(71,'(I5,3G16.8)') ik, ratio, source,
c     >        W(2,ik)
c      enddo
c
c      do ik = 1, ikmax
c         sum = interpol2p(U(2,ik)/W(2,ik), U(2,ik+1)/W(2,ik+1),
c     >             xf(ik), xf(ik+1), xc(ik))
c            write(71,'(I5,1G16.8)') ik, sum / U(1,ik)
c      enddo
c      STOP
c
c      do ik = 1, ikmax+1
c            write(71,'(I5,6G16.8)') ik,
c     >       U(1,ik), U(2,ik), U(3,ik), Ti_eV(ik)
c      enddo
c
c      do ik = 1, ikmax
c        if (ik.le.ikmid) then
c           U(1,ik) = U(1,0) * (1.0D0 - 0.5D0 *
c     >            (xc(ik)/xc(ikmid)) ** 0.5D0)
c           U(3,ik) = U(3,0) * (1.0D0 + 4.0D0 *
c     >             (xc(ik)/xc(ikmid)) ** 0.5D0)
c        else
c           U(1,ik) = U(1,ikmax+1 - ik)
c           U(3,ik) = U(3,ikmax+1 - ik)
c        endif
c      enddo
c
c      do ik = 1, ikmax+1
c        U(2,ik) = U(2,0) +
c     >   (U(2,ikmax+2) - U(2,0)) * xf(ik) / xc(ikmax+1)
c        F(1,ik) = U(2,ik)
c      enddo
c
      return
      end
c
c
      subroutine initQ(n_flag)
      implicit none
c
      include 'cfd_osm_com'
c
      integer ik, in, n_flag
      real*8 sum1,sum2,sum2_0,sum0,sum3,sum4,x_off
c
      do ik = 0, ikmax + 1
        if (xc(ik).le.x_xpt_0) then
          Q_uni(ik)  = 0.0D0
        elseif (xc(ik).ge.x_xpt_N) then
          Q_uni(ik)  = 0.0D0
        else
          Q_uni(ik)  = 1.0D0
        endif
      enddo
c
      if (a_PIN.ge.1) then
c
        call readQ(n_flag)
c
      elseif (a_seed.eq.1) then
c
        do ik = 0, ikmax + 1
c
          if (ik.le.ikmid) then
            Q(1,ik) =  MASS_0 * exp(- xc(ik)/L_iz)
            Q(2,ik) =  MOM_0  * exp(- xc(ik)/L_mom)
            Q(3,ik) =  POWi_0 * Q_uni(ik)
            Q(4,ik) =  POWe_0 * Q_uni(ik)
          else
            Q(1,ik) =  MASS_N * exp(-xc(ikmax+1-ik)/L_iz)
            Q(2,ik) =  - MOM_N * exp(- xc(ikmax+1-ik) / L_mom)
            Q(3,ik) =  POWi_N * Q_uni(ik)
            Q(4,ik) =  POWe_N * Q_uni(ik)
          endif
c
        enddo
      elseif (a_seed.eq.2) then
        do ik = 0, ikmax+1
          if (ik.le.ikmid) then
            x_off = max(0.0D0, xc(ik) - L_izoffset)
            Q(1,ik) =  MASS_0 *
     >               exp(- (x_off/L_iz)**2.0D0) * x_off ** 5.0D0
            Q(2,ik) =  MOM_0 *
     >               exp(- (x_off/L_mom)**2.0D0) * x_off ** 5.0D0
            Q(3,ik) =  POWi_0
            Q(4,ik) =  POWe_0
          else
            x_off = max(0.0D0, xc(ikmax+1-ik) - L_izoffset)
            Q(1,ik) =  MASS_N * exp(- (x_off/L_iz)**2.0D0) *
     >                 x_off ** 5.0D0
            Q(2,ik) =  - MOM_N * exp(- (x_off/L_mom)**2.0D0) *
     >                 x_off ** 5.0D0
            Q(3,ik) =  POWi_N
            Q(4,ik) =  POWe_N
          endif
        enddo
      else
c        write(0,*)  'SOURCES NOT AVAILABLE'
        write(71,*) 'SOURCES NOT AVAILABLE'
        STOP
      endif
c
      if (a_PIN.eq.0) then
c
        sum1 = 0.0D0
        sum2 = 0.0D0
        sum2_0 = 0.0D0
        sum3 = 0.0D0
        sum4 = 0.0D0
c
        do ik = 0, ikmax + 1
           sum1 = sum1 + Q(1,ik) * (xf(ik+1) - xf(ik))
           sum2 = sum2 + Q(2,ik) * (xf(ik+1) - xf(ik))
           sum3 = sum3 + Q(3,ik) * (xf(ik+1) - xf(ik))
           sum4 = sum4 + Q(4,ik) * (xf(ik+1) - xf(ik))
        enddo
c
        do ik = 0, ikmid
           sum2_0 = sum2_0 + Q(2,ik) * (xf(ik+1) - xf(ik))
        enddo
c
c        write(71,'(A14,3G16.8)') '<sum>:1',
c     >     sum1 / MASS_PIN,   sum1 / (MASS_0 + MASS_N)
c        write(71,'(A14,3G16.8)') '<sum>:2',
c     >     sum2,  MOM_PIN,    MOM_N - MOM_0
c        write(71,'(A14,3G16.8)') '<sum>:3',
c     >     sum3 / (POW_core + POW_RAD),  sum3 / POW_BC
c
c      call print_Q
c      call print_symmetry
c      call check_symmetry
c      write(71,*) 'sum123', sum1 / (MASS_0 + MASS_N)
c      write(71,*) 'Q2,sum2',           Q2, sum2
c      write(71,*) 'sum3/POW_PIN',      sum3/  POW_BC
c
        do ik = 0, ikmax + 1
           do in = 1, 4
             Q_B(in,ik) = 0.0D0
           enddo
           L_B(ik) = 1.0D0
           Q(1,ik) = Q(1,ik) * (Q1 * x_max / sum1)
           Q(3,ik) = Q(3,ik) * (Q3 * x_max / sum3)
           Q(4,ik) = Q(4,ik) * (Q4 * x_max / sum4)
           if (abs(sum2_0).gt.0.0D0) then
              Q(2,ik) = a_Q2 * Q(2,ik) * (MOM_0 / sum2_0)
           else
              Q(2,ik) = 0.0D0
           endif
        enddo
        Q(2,ikmid) = 0.0D0
c
      endif
c
      return
      end
c
c
      subroutine updateQperp
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer ik, in
      real*8 temp_2, temp_3, temp_4, temp_5
c
      call writeU(1)
c
      a_Qup  = 0
      temp_2 = s23_par_updtQp2
      temp_3 = s23_par_updtQp3
c
      if (Te_ctrl(ir_cfd,1) - W(4,0)/(W(1,0)*eV).ge.s23_par_updtqpte
     >      .and.min(W(4,0)/(W(1,0)*eV),Te_ctrl(ir_cfd,1)).lt.1.0) then
           power_flux(ir_cfd,1) = max(0.0D0, power_flux(ir_cfd,1) +
     >                               temp_3 * del_power_flux(ir_cfd,1))
           mom_flux(ir_cfd,1) = mom_flux(ir_cfd,1) + temp_2 *
     >                       del_mom_flux(ir_cfd,1)
           pow_in_0     = power_flux(ir_cfd,1) / 2.0
           mom_in_0     = mom_flux(ir_cfd,1)
           POW_core_i0  = pow_in_0 / 2.0D0
           POW_core_e0  = pow_in_0 / 2.0D0
           MOM_core_0   = mom_in_0
           a_Qup = 1
      endif
c
      if (Te_ctrl(ir_cfd,2) - W(4,ikmax+1)/(W(1,ikmax+1)*eV).ge.
     >      s23_par_updtqpte.and.min(Te_ctrl(ir_cfd,2),
     >      W(4,ikmax+1)/(W(1,ikmax+1)*eV)).lt.1.0) then
           power_flux(ir_cfd,2) = max(0.0D0, power_flux(ir_cfd,2) +
     >                               temp_3 * del_power_flux(ir_cfd,2))
           mom_flux(ir_cfd,2) = mom_flux(ir_cfd,2) + temp_2 *
     >                       del_mom_flux(ir_cfd,2)
           pow_in_N     = power_flux(ir_cfd,2) / 2.0
           mom_in_N     = mom_flux(ir_cfd,2)
           POW_core_iN  = pow_in_N / 2.0D0
           POW_core_eN  = pow_in_N / 2.0D0
           MOM_core_N   = mom_in_N
           a_Qup = 1
      endif
c
      if (a_Qup.eq.1) then
        POW_core_i   = POW_core_i0 + POW_core_iN
        POW_core_e   = POW_core_e0 + POW_core_eN
        MOM_core     = MOM_core_0  + MOM_core_N
        call initQ(0)
      endif
c
      return
      end
c
c
      subroutine readgrid
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer in,ik
c
      write(71,*) 'READ GRID'
c
      ikmax = ikmax_cfd(ir_cfd, targ_cfd)
      ikmid = (ikmax + 1) / 2
c
      do ik = 0, ikmax+1
         xc(ik) = xc_cfd(ir_cfd, targ_cfd, ik)
      enddo
c
      return
      end
c
c
      subroutine read_cfd
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer in,ik
c
      write(71,*) 'READ_CFD', targ_cfd, ikmax_cfd(ir_cfd,3)
c
      ikmax = ikmax_cfd(ir_cfd, targ_cfd)
      ikmid = (ikmax + 1) / 2
c
      do ik = 0, ikmax+1
         xc(ik) = xc_cfd(ir_cfd, targ_cfd, ik)
      enddo
c
      xf(0) = xc(0)
      xf(ikmax+2) = xc(ikmax+1)
      do ik = 1, ikmax+1
         xf(ik) = (xc(ik) + xc(ik-1)) / 2.0D0
      enddo
c
      s23_fract1 = fract_cfd(ir_cfd, targ_cfd, 1)
      s23_fract2 = fract_cfd(ir_cfd, targ_cfd, 2)
      s23_fract3 = fract_cfd(ir_cfd, targ_cfd, 3)
      s23_fract4 = fract_cfd(ir_cfd, targ_cfd, 4)
      s23_fract1_0 = s23_fract1
      s23_fract2_0 = s23_fract2
      s23_fract3_0 = s23_fract3
      s23_fract4_0 = s23_fract4
      s23_fract1_N = s23_fract1
      s23_fract2_N = s23_fract2
      s23_fract3_N = s23_fract3
      s23_fract4_0 = s23_fract4
c
      do ik = 0, ikmax+2
         do in = 1, 4
            U(in,ik) = U_cfd(ir_cfd, targ_cfd, in, ik)
            Q(in,ik) = Q_cfd(ir_cfd, targ_cfd, in, ik)
         enddo
      enddo
c
      return
      end
c
c
      subroutine mix_Q(delta)
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer in,ik
      real delta
c
      do ik = 0, ikmax+2
         do in = 1, 4
            Q(in,ik) = Q(in,ik) * delta +
     >        (1.0 - delta) * Q_cfd(ir_cfd, targ_cfd, in, ik)
         enddo
      enddo
c
      return
      end
c
c
      subroutine setQ_new
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer in,ik
c
      do ik = 0, ikmax+2
         do in = 1, 4
            Q_new(in,ik) = Q(in,ik)
         enddo
      enddo
c
      return
      end
c
c
      subroutine solvegrid
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer in,ik,ik_temp
      real*8  temp1, temp2, temp3
c
      write(71,*) 'SOLVEGRID'
c
      temp3 = s23_par_garelax
c
      call relax_RE(temp3)
c
      R(2,0) = 0.0D0
      do in = 1, ikmax_0 * s23_par_gaiter
        call grid_error
      enddo
c
      temp2 = 0.0D0
      do ik = 1, ikmax
        temp1 = abs(xc_cfd(ir_cfd, targ_cfd, ik)
     >                           - xc(ik)) / xc(ik)
        if (temp1.ge.temp2) then
            temp2   = temp1
            ik_temp = ik
        endif
      enddo
c
      write(71,*) 'GRID CHANGE', ik_temp, temp2
c
      temp1 = temp3
c
      do ik = 1, ikmax
         xc(ik) = xc_cfd(ir_cfd, targ_cfd, ik) * (1.0 - temp1) +
     >            temp1 * xc(ik)
      enddo
c
      xf(0) = xc(0)
      xf(ikmax+2) = xc(ikmax+1)
      do ik = 1, ikmax+1
         xf(ik) = (xc(ik) + xc(ik-1)) / 2.0D0
      enddo
c
      write(71,*) 'xc(1):', xc(1)
c
      return
      end
c
c
      subroutine readU(n_flag)
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      real*8 src_int_ab, x_a, x_b ,dx_ab, sum3, p0
      external src_int_ab
      integer in,ik,j,n_flag
c
c      write(71,*) 'READ_U'
c
c      do ik = 0, ikmax+1
c        x_a = xf(ik)
c        x_b = xf(ik+1)
c        dx_ab = x_b - x_a
c        W(1,ik) = src_int_ab(ne, x_a, x_b,
c     >              spts,sbnds,npts,intopt) / dx_ab
c        W(2,ik) = src_int_ab(vb, x_a, x_b,
c     >              spts,sbnds,npts,intopt) / dx_ab
c        W(3,ik) = src_int_ab (te, x_a, x_b,
c     >              spts,sbnds,npts,intopt) / dx_ab
c     >              * eV * W(1,ik) * a_Temp
c      enddo
c
c      if (a_half.eq.1) then
c        do in = 1, 3
c          do ik = 0, ikmid
c            if (in.eq.2) then
c              W(in,ikmax+1-ik) = - W(in,ik)
c            else
c              W(in,ikmax+1-ik) = W(in,ik)
c            endif
c          enddo
c        enddo
c      elseif(a_half.eq.2) then
c        do in = 1, 3
c          do ik = 0, ikmid
c            if (in.eq.2) then
c              W(in,ik) = - W(in,ikmax+1-ik)
c            else
c              W(in,ik) = W(in,ikmax+1-ik)
c            endif
c          enddo
c        enddo
c      endif
c
c      do in = 1, 3
c        W(in,ikmid) = (W(in,ikmid-1) + W(in,ikmid+1)) / 2.0D0
c      enddo
c
c      write(71,*) 'AFTER INVERSION'
c      call print_W
c
c      do j = 1, 5
c        call smoothU(0.1D0)
c      enddo
c
c      j = ikmax_cfd(ir_cfd, targ_cfd)
c      do ik = 0, j+2
c           aaa(ik) = U_cfd(ir_cfd, targ_cfd, 1, ik)
c           bbb(ik) = U_cfd(ir_cfd, targ_cfd, 2, ik)
c           ccc(ik) = U_cfd(ir_cfd, targ_cfd, 3, ik)
c           ddd(ik) = U_cfd(ir_cfd, targ_cfd, 4, ik)
c           uuu(ik) = xc_cfd(ir_cfd, targ_cfd, ik)
c           fff(ik) = xf_cfd(ir_cfd, targ_cfd, ik)
c      enddo
c
c      do ik = 0, ikmax+1
c         x_a = xf(ik)
c         x_b = xf(ik+1)
c         dx_ab = x_b - x_a
c         U(1,ik) = src_int_ab(aaa, x_a, x_b,
c     >                    uuu,fff,j+2,intopt) / dx_ab
c         U(2,ik) = src_int_ab(bbb, x_a, x_b,
c     >                    uuu,fff,j+2,intopt) / dx_ab
c         U(3,ik) = src_int_ab(ccc, x_a, x_b,
c     >                    uuu,fff,j+2,intopt) / dx_ab
c         U(4,ik) = src_int_ab(ddd, x_a, x_b,
c     >                    uuu,fff,j+2,intopt) / dx_ab
c      enddo
c
       if (n_flag.eq.0) then
         do ik = 0, ikmax+2
           xc(ik) = xc_cfd(ir_cfd, targ_cfd, ik)
           xf(ik) = xf_cfd(ir_cfd, targ_cfd, ik)
           do in = 1, 4
             U(in,ik) = U_cfd(ir_cfd, targ_cfd, in, ik)
           enddo
         enddo
       else
         j = ikmax_cfd(ir_cfd, targ_cfd)
         do ik = 0, j+2
           aaa(ik) = U_cfd(ir_cfd, targ_cfd, 1, ik)
           bbb(ik) = U_cfd(ir_cfd, targ_cfd, 2, ik)
           ccc(ik) = U_cfd(ir_cfd, targ_cfd, 3, ik)
           ddd(ik) = U_cfd(ir_cfd, targ_cfd, 4, ik)
           uuu(ik) = xc_cfd(ir_cfd, targ_cfd, ik)
           fff(ik) = xf_cfd(ir_cfd, targ_cfd, ik)
         enddo
c
         do ik = 0, ikmax+1
            x_a = xf(ik)
            x_b = xf(ik+1)
            dx_ab = x_b - x_a
            U(1,ik) = src_int_ab(aaa, x_a, x_b,
     >                    uuu,fff,j+2,intopt) / dx_ab
            U(2,ik) = src_int_ab(bbb, x_a, x_b,
     >                    uuu,fff,j+2,intopt) / dx_ab
            U(3,ik) = src_int_ab(ccc, x_a, x_b,
     >                    uuu,fff,j+2,intopt) / dx_ab
            U(4,ik) = src_int_ab(ddd, x_a, x_b,
     >                    uuu,fff,j+2,intopt) / dx_ab
         enddo
       endif
c
      return
      end
c
c
      subroutine readQ(n_flag)
      implicit none
      integer n_flag
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      real*8 src_int_ab, x_a, x_b ,dx_ab,
     >  sum1, sum2, sum2_0, sum2_N, sum3, sum3_0, sum3_N,
     >  sum4, sum4_0, sum4_N, sum_perp4,
     >  sum_perp1, sum_perp2, sum_perp3,
     >  temp1, temp2, temp3, p0
      real calcwidth
      external src_int_ab, calcwidth
      integer in,ik,j,ik1_max,ik2_max
c
      do ik = 0, ikmax+1
        x_a = xf(ik)
        x_b = xf(ik+1)
        dx_ab = x_b - x_a
c
        B_field(ik) = src_int_ab(magn_field, x_a, x_b,
     >                  spts,sbnds,npts,intopt) / dx_ab
        Q_grad(ik) = src_int_ab(perp_src, x_a, x_b,
     >                  spts,sbnds,npts,intopt) / dx_ab
        Q(1,ik) = mass * src_int_ab(part_src, x_a, x_b,
     >                  spts,sbnds,npts,intopt) / dx_ab
        Q(2,ik) = src_int_ab(mom_src, x_a, x_b,
     >                  spts,sbnds,npts,intopt) / dx_ab
        if (a_2T.eq.0) then
          Q(3,ik) = (src_int_ab (ionpow_src, x_a, x_b,
     >                        spts,sbnds,npts,intopt)
     >           + src_int_ab (elecpow_src, x_a, x_b,
     >             spts,sbnds,npts,intopt))/ dx_ab
        elseif (a_2T.eq.1) then
          Q(3,ik) = src_int_ab (ionpow_src, x_a, x_b,
     >                spts,sbnds,npts,intopt) / dx_ab
          Q(4,ik) = src_int_ab (elecpow_src, x_a, x_b,
     >                 spts,sbnds,npts,intopt) / dx_ab
        endif
      enddo
c
      Q(1,0) = 0.0D0
      Q(2,0) = 0.0D0
      Q(3,0) = 0.0D0
      Q(4,0) = 0.0D0
      Q(1,ikmax+1) = 0.0D0
      Q(2,ikmax+1) = 0.0D0
      Q(3,ikmax+1) = 0.0D0
      Q(4,ikmax+1) = 0.0D0
c
      do ik = 1, ikmax
        L_B(ik) = B_field(ik) * (xc(ik+1) - xc(ik-1)) /
     >                  (B_field(ik+1) - B_field(ik-1))
        if (s23_par_qperp34.eq.0) then
          Q_grad(ik) = 1.0D0
        elseif (s23_par_qperp34.eq.1) then
          Q_grad(ik) = (W(3,ik) + W(4,ik))
        elseif (s23_par_qperp34.eq.2) then
          Q_grad(ik) = B_field(ik)
        endif
        if (xc(ik).le.x_xpt_0) then
          Q_grad(ik) = 0.0D0
        elseif (xc(ik).ge.x_xpt_N) then
          Q_grad(ik) = 0.0D0
        endif
      enddo
c
      L_B(0)          = L_B(1)
      L_B(ikmax+1)    = L_B(ikmax)
      Q_grad(0)       = 0.0D0
      Q_grad(ikmax+1) = 0.0D0
c
      do ik = 0, ikmax+1
        if (a_perp.eq.0) then
          Q_perp(1,ik) = 1.0D0
          Q_perp(2,ik) = abs(U(2,ik))
        elseif (a_perp.eq.1) then
          Q_perp(1,ik) = abs(U(1,ik))
          Q_perp(2,ik) = abs(U(2,ik))
        elseif (a_perp.eq.2) then
          Q_perp(1,ik) = abs(U(1,ik))
          Q_perp(2,ik) = 1.0D0
        elseif (a_perp.eq.3) then
          Q_perp(1,ik) = abs(U(1,ik)) ** 0.5
          Q_perp(2,ik) = abs(U(2,ik))
        elseif (a_perp.eq.4) then
          Q_perp(1,ik) = abs(U(1,ik))
          Q_perp(2,ik) = abs(U(2,ik)) ** 0.5
        endif
      enddo
c
      do ik = 0, ikmax+1
        if (a_2T.eq.0) then
          if (ik.le.ikmid) then
           Q_perp(3,ik) = Q_uni(ik) * POW_core_0
          else
           Q_perp(3,ik) = Q_uni(ik) * POW_core_N
          endif
        elseif (a_2T.eq.1) then
          if (ik.le.ikmid) then
c
           Q_perp(3,ik) = Q_uni(ik) * POW_core_i0
           Q_perp(4,ik) = Q_uni(ik) * POW_core_e0
c
           if (s23_par_qperp34.gt.0) then
             Q_perp(3,ik) = Q_grad(ik) * POW_core_i0
             Q_perp(4,ik) = Q_grad(ik) * POW_core_e0
           endif
c
          else
c
           Q_perp(3,ik) = Q_uni(ik) * POW_core_iN
           Q_perp(4,ik) = Q_uni(ik) * POW_core_eN
c
           if (s23_par_qperp34.gt.0) then
             Q_perp(3,ik) = Q_grad(ik) * POW_core_iN
             Q_perp(4,ik) = Q_grad(ik) * POW_core_eN
           endif
c
          endif
        endif
      enddo
c
      if (a_symm.eq.1) then
c
        Q(2,ikmid) = 0.0D0
        do ik = 0, ikmax+1
          Q_perp(2,ik) = 0.0D0
        enddo
c
        if (a_half.eq.1) then
          do ik = 0, ikmid
            Q(1,ikmax+1-ik) = Q(1,ik)
            Q(2,ikmax+1-ik) = - Q(2,ik)
            Q(3,ikmax+1-ik) = Q(3,ik)
            Q(4,ikmax+1-ik) = Q(4,ik)
            Q_perp(1,ikmax+1-ik) = Q_perp(1,ik)
            Q_perp(3,ikmax+1-ik) = Q_perp(3,ik)
            Q_perp(4,ikmax+1-ik) = Q_perp(4,ik)
          enddo
         elseif(a_half.eq.2) then
          do ik = 0, ikmid
            Q(1,ik) = Q(1,ikmax+1-ik)
            Q(2,ik) = - Q(2,ikmax+1-ik)
            Q(3,ik) = Q(3,ikmax+1-ik)
            Q(4,ik) = Q(4,ikmax+1-ik)
            Q_perp(1,ik) = Q_perp(1,ikmax+1-ik)
            Q_perp(3,ik) = Q_perp(3,ikmax+1-ik)
            Q_perp(4,ik) = Q_perp(4,ikmax+1-ik)
          enddo
        endif
c
      endif
c
      Q_perp(1,0) = 0.0D0
      Q_perp(2,0) = 0.0D0
      Q_perp(3,0) = 0.0D0
      Q_perp(4,0) = 0.0D0
      Q_perp(1,ikmax+1) = 0.0D0
      Q_perp(2,ikmax+1) = 0.0D0
      Q_perp(3,ikmax+1) = 0.0D0
      Q_perp(4,ikmax+1) = 0.0D0
c
      sum1   = 0.0D0
      sum2   = 0.0D0
      sum2_0 = 0.0D0
      sum2_N = 0.0D0
      sum3   = 0.0D0
      sum3_0 = 0.0D0
      sum3_N = 0.0D0
      sum4   = 0.0D0
      sum4_0 = 0.0D0
      sum4_N = 0.0D0
      sum_perp1 = 0.0D0
      sum_perp2 = 0.0D0
      sum_perp3 = 0.0D0
      sum_perp4 = 0.0D0
c
      do ik = 0, ikmax + 1
         sum1 = sum1 + Q(1,ik) * (xf(ik+1) - xf(ik))
         sum2 = sum2 + Q(2,ik) * (xf(ik+1) - xf(ik))
         sum3 = sum3 + Q(3,ik) * (xf(ik+1) - xf(ik))
         sum4 = sum4 + Q(4,ik) * (xf(ik+1) - xf(ik))
         sum_perp1 = sum_perp1 + Q_perp(1,ik)
     >                      * (xf(ik+1) - xf(ik))
         sum_perp2 = sum_perp2 + Q_perp(2,ik)
     >                      * (xf(ik+1) - xf(ik))
         sum_perp3 = sum_perp3 + Q_perp(3,ik)
     >                      * (xf(ik+1) - xf(ik))
         sum_perp4 = sum_perp4 + Q_perp(4,ik)
     >                      * (xf(ik+1) - xf(ik))
      enddo
c
      do ik = 0, ikmid
         sum2_0 = sum2_0 + Q(2,ik) * (xf(ik+1) - xf(ik))
         sum3_0 = sum3_0 + Q(3,ik) * (xf(ik+1) - xf(ik))
         sum4_0 = sum4_0 + Q(4,ik) * (xf(ik+1) - xf(ik))
      enddo
c
      do ik = ikmid+1, ikmax+1
         sum2_N = sum2_N + Q(2,ik) * (xf(ik+1) - xf(ik))
         sum3_N = sum3_N + Q(3,ik) * (xf(ik+1) - xf(ik))
         sum4_N = sum4_N + Q(4,ik) * (xf(ik+1) - xf(ik))
      enddo
c
c      if (n_flag.eq.0) then
c        write(71,'(A14,3G16.8)') '<QPIN>:1',
c     >     MASS_PIN_0 / MASS_0,
c     >     MASS_PIN_N / MASS_N,
c     >     MASS_PIN / (MASS_0 + MASS_N)
c        write(71,'(A14,3G16.8)') '<QPIN>:2',
c     >     MOM_PIN_0 / MOM_0,
c     >     MOM_PIN_N / MOM_N,
c     >     MOM_PIN / (MOM_0 + MOM_N)
c        write(71,'(A14,3G16.8)') '<QPIN>:3',
c     >     POW_PIN_0 / POW_core_0,
c     >     POW_PIN_N / POW_core_N,
c     >     POW_PIN   / POW_core
c        write(71,'(A14,3G16.8)') '<QRAD>:3',
c     >     POW_RAD_0 / POW_core_0,
c     >     POW_RAD_N / POW_core_N,
c     >     POW_RAD   / POW_core
c        write(71,'(A14,3G16.8)') '<QRAD>:3',
c     >     POW_core_0 / POW_BC_0,
c     >     POW_core_N / POW_BC_N,
c     >     POW_core   / POW_BC
c        write(71,'(A14,3G16.8)') '<sum>:1',
c     >     sum1, MASS_PIN, sum_perp1
c        write(71,'(A14,4G16.8)') '<sum>:2',
c     >     sum2, MOM_PIN,  sum_perp2, MOM_core
c        write(71,'(A14,3G16.8)') '<sum>:3',
c     >     sum3, POW_PIN,  sum_perp3
c      endif
c
      s23_fract1 = (MASS_PIN  / (MASS_0 + MASS_N))
      s23_fract2 = (MOM_PIN_N / MOM_N)
c
      s23_fract1_0 = (MASS_PIN_0 / MASS_0)
      s23_fract2_0 = (MOM_PIN_0  / MOM_0)
c
      s23_fract1_N = (MASS_PIN_N / MASS_N)
      s23_fract2_N = (MOM_PIN_N  / MOM_N)
c
      if (a_2T.eq.0) then
        s23_fract3   = (POW_RAD    / POW_core)
        s23_fract3_0 = (POW_RAD_0  / POW_core_0)
        s23_fract3_N = (POW_RAD_N  / POW_core_N)
      elseif (a_2T.eq.1) then
        s23_fract3   = (POW_RAD_i   / POW_core_i)
        s23_fract3_0 = (POW_RAD_i0  / POW_core_i0)
        s23_fract3_N = (POW_RAD_iN  / POW_core_iN)
        s23_fract4   = (POW_RAD_e   / POW_core_e)
        s23_fract4_0 = (POW_RAD_e0  / POW_core_e0)
        s23_fract4_N = (POW_RAD_eN  / POW_core_eN)
      endif
c
c      if (n_flag.eq.0) then
c        write(71,'(A14,3G16.8)') 'f1',
c     >     s23_fract1, s23_fract1_0, s23_fract1_N
c        write(71,'(A14,3G16.8)') 'f2',
c     >     s23_fract2, s23_fract2_0, s23_fract2_N
c        write(71,'(A14,3G16.8)') 'f3',
c     >     s23_fract3, s23_fract3_0, s23_fract3_N
c        write(71,'(A14,3G16.8)') 'f4',
c     >     s23_fract4, s23_fract4_0, s23_fract4_N
c      endif
c
c     ########### CHANGED MASS_0 to MASS_BC_0 ###############
c
      if (sum_perp1.ne.0.0D0) then
        do ik = 0, ikmax+1
           Q_perp(1,ik) = Q_perp(1,ik) *
     >             (MASS_0 + MASS_N - MASS_PIN) / sum_perp1
        enddo
      endif
c
      if (sum_perp2.ne.0.0D0) then
        do ik = 0, ikmax+1
           Q_perp(2,ik) = Q_perp(2,ik) *
     >              (MOM_core - MOM_PIN) / sum_perp2
        enddo
      endif
c
      if (a_2T.eq.0) then
        if (sum_perp3.ne.0.0D0) then
          do ik = 0, ikmax+1
            Q_perp(3,ik) = Q_perp(3,ik) *
     >                               POW_core / sum_perp3
          enddo
        endif
      elseif (a_2T.eq.1) then
        if (sum_perp3.ne.0.0D0) then
          do ik = 0, ikmax+1
            Q_perp(3,ik) = Q_perp(3,ik) *
     >                               POW_core_i / sum_perp3
          enddo
        endif
c
        if (sum_perp4.ne.0.0D0) then
          do ik = 0, ikmax+1
            Q_perp(4,ik) = Q_perp(4,ik) *
     >                               POW_core_e / sum_perp4
          enddo
        endif
      endif
c
      if (sum1.ne.0.0D0) then
        do ik = 0, ikmax+1
           Q(1,ik) = Q(1,ik) * MASS_PIN / sum1
        enddo
      endif
c
c      if (sum2.ne.0.0D0) then
c        do ik = 0, ikmax+1
c           Q(2,ik) = Q(2,ik) * MOM_PIN / sum2
c        enddo
c      endif
c
      if (a_2T.eq.0) then
        if (sum3.ne.0.0D0) then
          do ik = 0, ikmax+1
             Q(3,ik) = Q(3,ik) * POW_RAD / sum3
          enddo
        endif
      elseif (a_2T.eq.1) then
        if (sum3.ne.0.0D0) then
          do ik = 0, ikmax+1
c             Q(3,ik) = Q(3,ik) * POW_RAD_i / sum3
          enddo
        endif
c
        if (sum4.ne.0.0D0) then
          do ik = 0, ikmax+1
c             Q(4,ik) = Q(4,ik) * POW_RAD_e / sum4
          enddo
        endif
      endif
c
      do ik = 0, ikmax+1
         Q(1,ik) = Q(1,ik) + Q_perp(1,ik)
         Q(2,ik) = Q(2,ik) + Q_perp(2,ik)
         Q(3,ik) = Q(3,ik) + Q_perp(3,ik)
         Q(4,ik) = Q(4,ik) + Q_perp(4,ik)
      enddo
c
c      if (sum3.ne.0.0D0) then
c        do ik = 0, ikmax + 1
c           Q(3,ik) = Q3 * x_max * Q_uni(ik) / sum_uni +
c     >               Q(3,ik) * abs(POW_RAD / sum3)
c        enddo
c      else
c        do ik = 0, ikmax + 1
c           Q(3,ik) = Q3
c        enddo
c      endif
c
c      write(71,*) 'sumuni,x_max,Q3',sum_uni,x_max*Q3
c
c      temp1 = 0.0D0
c      do ik = 0, ikmid
c         temp2 = Q(1,ikmid) + const1 * (W(1,ikmid) / sum_n)
c         if (temp2.ge.temp1) then
c            temp1 = temp2
c            ik2_max = ik
c         endif
c      enddo
c
c      temp1 = Q(1,ikmid) + const1 * (W(1,ikmid) / sum_n)
c      temp2 = Q(1,ikmid-2) + const1 * (W(1,ikmid-2) / sum_n)
c      do ik = 0, ikmax + 1
c        if (xc(ik1_max).ge.x_max/4.0D0) then
c          write(71,*) 'Qopt1'
c          Q(1,ik) = Q(1,ik) + (const1 / x_max)
c
c        if (a_perp.eq.0) then
c          write(71,*) 'Qopt2'
c          Q(1,ik) = Q(1,ik) + (const1 / x_max)
c        elseif (a_perp.eq.1.or.a_perp.eq.5) then
c          write(71,*) 'Qopt3'
c          Q(1,ik) = Q(1,ik) + const1 * (W(1,ik) / sum_n)
c        elseif (a_perp.eq.2) then
c          write(71,*) 'Qopt4'
c          Q(1,ik) = Q(1,ik) + const1 *
c     >                          (dsqrt(W(1,ik)) / sum_nh)
c        elseif (a_perp.eq.3) then
c          write(71,*) 'Qopt5'
c          Q(1,ik) = Q(1,ik)
c        elseif(a_perp.eq.4) then
c          write(71,*) 'Qopt6'
c          if ((sum_perp/const1).le.0.0D0) then
c             Q(1,ik) = Q(1,ik) + const1 *
c     >                     (Q(4,ik) / sum_perp)
c          else
c            Q(1,ik) = Q(1,ik) + const1 * W(1,ik) / sum_n
c          endif
c        endif
c      enddo
c
      return
      end
c
c
      subroutine writeU(n_flag)
      implicit none
c
      integer n_flag
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
      include 'pin_cfd'
c
      real*8 ne1,vb1,te1, c_s_0, c_s_smax, temp_4, temp_5
      real*8 src_int_ab, x_a, x_b , dx_ab, sum1
      external src_int_ab
      integer in,ik
c
      ikmax_cfd(ir_cfd, targ_cfd) = ikmax
c
      fract_cfd(ir_cfd, targ_cfd, 1) = s23_fract1
      fract_cfd(ir_cfd, targ_cfd, 2) = s23_fract2
      fract_cfd(ir_cfd, targ_cfd, 3) = s23_fract3
      fract_cfd(ir_cfd, targ_cfd, 4) = s23_fract4
c
      do ik = 0, ikmax+2
         xc_cfd(ir_cfd, targ_cfd, ik) = xc(ik)
         xf_cfd(ir_cfd, targ_cfd, ik) = xf(ik)
         do in = 1, 4
            U_cfd(ir_cfd, targ_cfd, in, ik) = U(in,ik)
         enddo
      enddo
c
      if (a_pin.ge.1) then
        if (a_2T.eq.0) then
          del_power_flux(ir_cfd,1) = POW_core *
     >     (POW_BC_0 - abs(F(3,0) + G(3,0))) / POW_BC
c
          del_power_flux(ir_cfd,2) = POW_core *
     >     (POW_BC_N - abs(F(3,ikmax+2) + G(3,ikmax+2)))
     >      / POW_BC
c
        elseif (a_2T.eq.1) then
c
c          temp_4 = (POW_BC_e0 - abs(F(4,0) + G(4,0))) / POW_BC_e
c          temp_5 = (POW_BC_eN - abs(F(4,ikmax+2) + G(4,ikmax+2)))
c     >      / POW_BC_e
c
           temp_4 = (Te_ctrl(ir_cfd,1) - W(4,0)/(W(1,0)*eV))
     >        / max(Te_ctrl(ir_cfd,1), W(4,0)/(W(1,0)*eV))
           temp_5 = (Te_ctrl(ir_cfd,2) - W(4,ikmax+1)/(W(1,ikmax+1)*eV))
     >        / max(Te_ctrl(ir_cfd,2), W(4,ikmax+1)/(W(1,ikmax+1)*eV))
c
           if (cfd_failed(ir_cfd).eq.0) then
             if (temp_4.ge.0.0.and.temp_5.le.0.0) then
                del_power_flux(ir_cfd,1) = POW_core_e * temp_4
                del_power_flux(ir_cfd,2) = POW_core_e * temp_5
             elseif (temp_4.le.0.0.and.temp_5.ge.0.0) then
                del_power_flux(ir_cfd,1) = POW_core_e * temp_4
                del_power_flux(ir_cfd,2) = POW_core_e * temp_5
             elseif (temp_4.le.0.0.and.temp_5.le.0.0) then
                del_power_flux(ir_cfd,1) = POW_core_e * temp_4
                del_power_flux(ir_cfd,2) = POW_core_e * temp_5
             else
                del_power_flux(ir_cfd,1) = POW_core_e * temp_4
                del_power_flux(ir_cfd,2) = POW_core_e * temp_5
             endif
           else
             del_power_flux(ir_cfd,1) = 0.0D0
             del_power_flux(ir_cfd,2) = 0.0D0
           endif
        endif
c
        mom_targ(ir_cfd,1) = abs(F(2,0) + G(2,0))
        mom_targ(ir_cfd,2) = abs(F(2,ikmax+1) + G(2,ikmax+1))
c
        if (cfd_failed(ir_cfd).eq.0) then
          del_mom_flux(ir_cfd,1)   = - abs(MOM_BC_0) *
     >     (MASS_BC_0 - abs(F(1,0))) / MASS_BC
          del_mom_flux(ir_cfd,2)   =   abs(MOM_BC_N) *
     >     (MASS_BC_N - abs(F(1,ikmax+2))) / MASS_BC
        else
          del_mom_flux(ir_cfd,1) = 0.0D0
          del_mom_flux(ir_cfd,2) = 0.0D0
        endif
c
        del_part_flux(ir_cfd,1)   = abs(MASS_BC_0/mass) *
     >     (MASS_BC_0 - abs(F(1,0))) / MASS_BC
c
        del_part_flux(ir_cfd,2)   = abs(MASS_BC_N/mass) *
     >     (MASS_BC_N - abs(F(1,ikmax+2))) / MASS_BC
c
c        write(71,'(A30,2G16.8)') 'DEL_POW_IN 0,N  ',
c     >     del_power_flux(ir_cfd,1) / POW_core,
c     >     del_power_flux(ir_cfd,2) / POW_core
c        write(71,'(A30,2G16.8)') 'DEL_MOM_IN 0,N  ',
c     >     del_mom_flux(ir_cfd,1) / MOM_BC,
c     >     del_mom_flux(ir_cfd,2) / MOM_BC
c        write(71,'(A30,2G16.8)') 'DEL_MASS_IN 0,N ',
c     >     del_part_flux(ir_cfd,1) / (MASS_BC/mass),
c     >     del_part_flux(ir_cfd,2) / (MASS_BC/mass)
c
      endif
c
      if (n_flag.eq.0) then
c
        sum1 = 0.0D0
        do ik = 0, ikmax + 1
c           sum1 = sum1 + (Q_ei(ik) + Q_joule(ik))
c     >                    * (xf(ik+1) - xf(ik))
           sum1 = sum1 + Q_ei(ik) * (xf(ik+1) - xf(ik))
        enddo
        power_ei(ir_cfd) = sum1
c
        do ik = 0, ikmax+2
           do in = 1, 4
             Q_cfd(ir_cfd, targ_cfd, in, ik) = Q(in,ik)
           enddo
        enddo
c
        if (a_pin.eq.0.and.a_cfdfile.eq.0) then
          do ik = 0, ikmax+1
            aaa(ik) = Q(1,ik) / mass
          enddo
c
          do ik = 1, npts
             dx_ab = sbnds(ik) - sbnds(ik-1)
             if (dx_ab.ne.0.0) then
                x_a = sbnds(ik-1)
                x_b = sbnds(ik)
                dx_ab = x_b - x_a
                pinion_last(ik,ir_cfd) =
     >                    src_int_ab(aaa, x_a, x_b,
     >                    xc,xf, ikmax+2, intopt) / dx_ab
             endif
          enddo
        endif
c
        do ik = 0, ikmax+1
          aaa(ik) = W(1,ik)
          bbb(ik) = W(2,ik)
          ccc(ik) = W(3,ik)/ (W(1,ik) * a_temp * eV)
          ddd(ik) = W(4,ik)/ (W(1,ik) * a_temp * eV)
        enddo
c
        do ik = 1, npts
           dx_ab = sbnds(ik) - sbnds(ik-1)
           if (dx_ab.ne.0.0) then
              x_a = sbnds(ik-1)
              x_b = sbnds(ik)
              dx_ab = x_b - x_a
              ne(ik) = src_int_ab(aaa, x_a, x_b,
     >                    xc,xf, ikmax+2, intopt) / dx_ab
              vb(ik) = src_int_ab(bbb, x_a, x_b,
     >                    xc,xf, ikmax+2, intopt) / dx_ab
              ti(ik) = src_int_ab(ccc, x_a, x_b,
     >                    xc, xf, ikmax+2, intopt) / dx_ab
              if (a_2T.eq.0) then
                 te(ik) = ti(ik)
              else
                 te(ik) = src_int_ab(ddd, x_a, x_b,
     >                    xc, xf, ikmax+2, intopt) / dx_ab
              endif
           endif
        enddo
c
c       NEARLY MONOTONIC TEMPERATURE RISE
c
        if (s23_par_limitte.eq.1) then
c
          do ik = 2, npts_mid
            if (te(ik).le.te(ik-1)*s23_par_celldte) then
               te(ik) = te(ik-1) * s23_par_celldte
            endif
            if (ti(ik).le.ti(ik-1)*s23_par_celldte) then
               ti(ik) = ti(ik-1) * s23_par_celldte
            endif
          enddo
c
          do ik = npts - 1, npts_mid + 1, -1
            if (te(ik).le.te(ik+1)*s23_par_celldte) then
               te(ik) = te(ik+1) * s23_par_celldte
            endif
            if (ti(ik).le.ti(ik+1)*s23_par_celldte) then
               ti(ik) = ti(ik+1) * s23_par_celldte
            endif
          enddo
c
        endif
c
        if (a_half.eq.1) then
c
          ne_0   = W(1,0)
          isat_0 = U(2,0) / mass
          ti_0   = W(3,0) / (W(1,0) * a_Temp * eV)
          if (a_2T.eq.0) then
            te_0  = ti_0
            c_s_0 = dsqrt(W(3,0) / U(1,0))
          else
            te_0  = W(4,0) / (W(1,0) * a_Temp * eV)
            c_s_0 = dsqrt((W(3,0) + W(4,0)) / U(1,0))
          endif
          Mach_0  = abs(W(2,0) / c_s_0)
c
        elseif (a_half.eq.2) then
c
          ne_smax   = W(1,ikmax+1)
          isat_smax = U(2,ikmax+1) / mass
          ti_smax   = W(3,ikmax+1) / (W(1,ikmax+1) * a_Temp * eV)
          if (a_2T.eq.0) then
            te_smax  = ti_smax
            c_s_smax = dsqrt(W(3,ikmax+1) / U(1,ikmax+1))
          else
            te_smax  = W(4,ikmax+1) / (W(1,ikmax+1) * a_Temp * eV)
            c_s_smax = dsqrt((W(3,ikmax+1) + W(4,ikmax+1))
     >                  / U(1,ikmax+1))
          endif
          Mach_smax  = abs(W(2,ikmax+1) / c_s_smax)
c
        elseif (a_half.eq.0) then
c
          ne_0   = W(1,0)
          isat_0 = U(2,0) / mass
          ti_0   = W(3,0) / (W(1,0) * a_Temp * eV)
          if (a_2T.eq.0) then
            te_0  = ti_0
            c_s_0 = dsqrt(W(3,0) / U(1,0))
          else
            te_0  = W(4,0) / (W(1,0) * a_Temp * eV)
            c_s_0 = dsqrt((W(3,0) + W(4,0)) / U(1,0))
          endif
          Mach_0  = abs(W(2,0) / c_s_0)
c
          ne_smax   = W(1,ikmax+1)
          isat_smax = U(2,ikmax+1) / mass
          ti_smax   = W(3,ikmax+1) / (W(1,ikmax+1) * a_Temp * eV)
          if (a_2T.eq.0) then
            te_smax  = ti_smax
            c_s_smax = dsqrt(W(3,ikmax+1) / U(1,ikmax+1))
          else
            te_smax  = W(4,ikmax+1) / (W(1,ikmax+1) * a_Temp * eV)
            c_s_smax = dsqrt((W(3,ikmax+1) + W(4,ikmax+1))
     >                  / U(1,ikmax+1))
          endif
          Mach_smax  = abs(W(2,ikmax+1) / c_s_smax)
c
        endif
c
      endif
c
      return
      end
c
c
      subroutine mkfromW
      implicit none
c
c     calculates the dynamic viscosity mu at each GRID pt.
c     ie.
c        mu = 0.96 * pi * tauii
c        ka = 2.6 * cv * pi * tauii
c     where
c        pi = n*k*Ti, ion pressure
c        tuaii, i-i collision time,
c        tauii = Ti^3/2 * Ai^1/2 / (ne * Lambda)
c        Ti = ion temp in eV, and Ai = ion mass # = mi/mp
c     Note that the units of mu are S.I., Pa*s
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
      integer ik
      real*8  n_i, n_e, p_i, p_e, Ti_J, Te_J,
     >        a_2, a_3, chi_e, chi_i,
     >        eta_e, eta_i, lambda_C, G2, dvdx, dx,
     >        dTdx, TJ_a, TJ_b, TJ_c, c_e
c
      do ik = 0, ikmax+1
c
        if (a_2T.eq.0) then
c
          n_i = W(1,ik)
          n_e = Z_i * n_i
          p_i = W(3,ik)
          p_e = Z_i * p_i
          Ti_eV(ik) = p_i / (a_Temp * n_i * eV)
          Te_eV(ik) = Ti_eV(ik)
          Ti_J  = p_i / (a_Temp * n_i)
          Te_J = Ti_J
          lambda_C = 7.31D0 + max(0.0D0,
     >     log((Te_eV(ik) ** 1.5D0) / dsqrt(n_e / 1.0D20)))
          tau_i(ik) = 2.09D13 * (Ti_eV(ik) ** 1.5D0) * sqrt(A_i)
     >                 / (Z_i ** 4.0D0 * n_i * lambda_C)
          tau_e(ik) = 3.44D11 * (Te_eV(ik) ** 1.5D0)
     >                 / (n_e * lambda_C)
          chi_i = 3.9D0  * p_i * tau_i(ik) / mass_i
          chi_e = 3.16D0 * p_e * tau_e(ik) / mass_e
          eta_i = 1.11D0 * p_i * tau_i(ik)
          eta_e = 0.73D0 * p_e * tau_e(ik)
c
          Q_ei(ik) = 0.0D0
c
          mu(ik) = eta_i + eta_e
          ka(ik) = chi_i + chi_e
c
          a_2 = mu(ik) / U(1,ik)
          a_3 = ka(ik) / (W(1,ik) * c_v_0)
          nu_max(ik) = max(a_2, a_3)
c
        else
c
          n_i = W(1,ik)
          n_e = Z_i * W(1,ik)
          p_i = W(3,ik)
          p_e = W(4,ik)
          Ti_eV(ik) = p_i / (n_i * eV)
          Te_eV(ik) = p_e / (n_e * eV)
          Ti_J      = p_i / n_i
          Te_J      = p_e / n_e
          lambda_C = 7.31D0 + max(0.0D0,
     >     log((Te_eV(ik) ** 1.5D0) / dsqrt(n_e / 1.0D20)))
          tau_i(ik) = 2.09D13 * (Ti_eV(ik) ** 1.5D0) * sqrt(A_i)
     >                 / (Z_i ** 4.0D0 * n_i * lambda_C)
          tau_e(ik) = 3.44D11 * (Te_eV(ik) ** 1.5D0)
     >                 / (n_e * lambda_C)
          tau_ei(ik) = tau_e(ik) * (mass_i / (mass_e * 3.0D0))
          chi_i = 3.9D0  * p_i * tau_i(ik) / mass_i
          chi_e = 3.16D0 * p_e * tau_e(ik) / mass_e
          eta_i = 1.11D0 * p_i * tau_i(ik)
          eta_e = 0.73D0 * p_e * tau_e(ik)
c
          if (s23_par_qeiflag.eq.0) then
             tau_ei(ik) = tau_ei(ik) * 1.0D3
          endif
c
          Q_ei(ik) = (W(4,ik) - W(3,ik)) / tau_ei(ik)
c
          ka_i(ik) = chi_i
          ka_e(ik) = chi_e
          mu(ik)   = eta_i
c
          a_2 = mu(ik) / U(1,ik)
          a_3 = ka(ik) / (W(1,ik) * c_v_0)
          nu_max(ik) = max(a_2, a_3)
c
        endif
      enddo
c
      do ik = 1, ikmax
         dx   = xc(ik+1) - xc(ik-1)
         dvdx = abs((W(2,ik+1) - W(2,ik-1)) / dx)
         TJ_a = W(4,ik-1) / (W(1,ik-1) * a_Temp)
         TJ_b = W(4,ik)   / (W(1,ik)   * a_Temp)
         TJ_c = W(4,ik+1) / (W(1,ik+1) * a_Temp)
         dTdx = abs((TJ_c - TJ_a) / dx)
         c_e  = dsqrt(8.0D0 * TJ_b / (3.14 * mass_e))
c
         mu(ik)   = mu(ik)   /
     >              (1.0D0 + 1.95D0 * tau_i(ik) * dvdx)
c
         if (s23_par_chie.eq.1) then
           ka_e(ik) = ka_e(ik) /
     >         (1.0D0 + 3.16D0 * tau_e(ik) * dTdx /
     >                      (0.2D0 * c_e * mass_e))
         elseif (s23_par_chie.eq.2) then
           ka_e(ik) = ka_e(ik) * 3.2 / 3.16
         endif
c
         if (a_2T.eq.0) then
           Q_joule(ik) = 0.0D0
         else
           Q_joule(ik) = - W(2,ik) * (W(4,ik+1) - W(4,ik-1)) / dx
         endif
c
         if (s23_par_joule.eq.0) then
           Q_joule(ik) = 0.0D0
         endif
c
      enddo
c
      mu(0)             = mu(1)
      mu(ikmax+1)       = mu(ikmax)
      ka_e(0)           = ka_e(1)
      ka_e(ikmax+1)     = ka_e(ikmax)
c
      if (a_2T.eq.0) then
        Q_joule(0)        = 0.0D0
        Q_joule(ikmax+1)  = 0.0D0
      else
        Q_joule(0)       = 0.0D0
        Q_joule(ikmax+1) = 0.0D0
        Q_ei(0)          = 0.0D0
        Q_ei(ikmax+1)    = 0.0D0
      endif
c
      return
      end
c
c
      subroutine dtfromU
      implicit none
c
      include 'cfd_osm_com'
      integer ik
      real*8 c_s, netv, dt_F, dt_G, dx_min,
     >   dt_min, interpol2p, a_2, a_3, beta, nu_min,
     >   Pe_inv
c
      dt_min = 10.0D0
      do ik = 0, ikmax+1
c
         dx_min = xf(ik+1) - xf(ik)
         netv = c_s_j(ik) + dabs(W(2,ik))
         a_2 = mu(ik) / U(1,ik)
         a_3 = ka(ik) / (W(1,ik) * c_v_0)
         beta = netv * dx_min
         nu_min = min(a_2, a_3)
         Pe_inv = nu_max(ik) / beta
c
         dt_F = dx_min / netv
c
         if (a_2T.eq.0) then
           dtc(ik)  = a_dt_F * dt_F
         else
           dtc(ik)  = a_dt_F * dt_F
         endif
c
         if (dtc(ik).lt.dt_min) then
           dt_min = dtc(ik)
         endif
      enddo
c
      do ik = 1, ikmax+1
         dtf(ik) = (dtc(ik-1) +  dtc(ik)) / 2.0D0
      enddo
c
      return
      end
c
c
      subroutine updateU
      implicit none
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
      real*8 source, ratio, flux, interpol2p, U_1, dummy
      integer ik,in,im

      if (a_flag.eq.1.and.a_2T.eq.0) then
        do ik = 1, ikmax
            do in = 1, 3
              ffff(in,ik) = 0.0D0
              do im = 1, 3
                ffff(in,ik) = ffff(in,ik) +  tinv_mat(in,im,ik) *
     >                        (R(im,ik)) * dtc(ik)
              enddo
            enddo
        enddo
c
        do in = 1, 3
           do ik = 1, ikmax
             aaa(ik) = aaaa(in,ik)
             bbb(ik) = bbbb(in,ik)
             ccc(ik) = cccc(in,ik)
             fff(ik) = ffff(in,ik)
             dummy =  (dabs(aaa(ik))+dabs(ccc(ik))) /
     >                 bbb(ik)
          enddo
c
          call tridiag(aaa,bbb,ccc,fff,uuu,ikmax)
c
          do ik = 1, ikmax
             ffff(in,ik) = uuu(ik)
          enddo
        enddo
c
        do in = 1, 3
          do ik = 1, ikmax
            ffff_new(in,ik) = 0.0D0
            do im = 1, 3
                ffff_new(in,ik) = ffff_new(in,ik) +
     >             t_mat(in,im,ik) * ffff(im,ik)
            enddo
          enddo
        enddo
c
        do ik = 1, ikmax
            do in = 1, 3
               U(in,ik) = U0(in,ik) + ffff_new(in,ik)
            enddo
        enddo
c
      elseif (a_flag.eq.0.and.a_2T.eq.0) then
c
        call blocktridiag(aaaa3,bbbb3,cccc3,ffff3,uuuu3,ikmax,3)
c
        do in = 1, 3
          do ik = 1, ikmax
             U(in,ik) = U0(in,ik) + uuuu3(in,ik)
          enddo
        enddo
c
      elseif (a_flag.eq.0.and.a_2T.eq.1) then
c
c        do ik = 1, ikmax
c          write(71,*) 'A(ik) , ik = ', ik
c          write(71,*) aaaa4(1,1,ik), aaaa4(1,2,ik),
c     >                aaaa4(1,3,ik), aaaa4(1,4,ik)
c          write(71,*) aaaa4(2,1,ik), aaaa4(2,2,ik),
c     >                aaaa4(2,3,ik), aaaa4(2,4,ik)
c          write(71,*) aaaa4(3,1,ik), aaaa4(3,2,ik),
c     >                aaaa4(3,3,ik), aaaa4(3,4,ik)
c          write(71,*) aaaa4(4,1,ik), aaaa4(4,2,ik),
c     >                aaaa4(4,3,ik), aaaa4(4,4,ik)
c        enddo
c
c        do ik = 1, ikmax
c          write(71,*) 'B(ik) , ik = ', ik
c          write(71,*) bbbb4(1,1,ik), bbbb4(1,2,ik),
c     >                bbbb4(1,3,ik), bbbb4(1,4,ik)
c          write(71,*) bbbb4(2,1,ik), bbbb4(2,2,ik),
c     >                bbbb4(2,3,ik), bbbb4(2,4,ik)
c          write(71,*) bbbb4(3,1,ik), bbbb4(3,2,ik),
c     >                bbbb4(3,3,ik), bbbb4(3,4,ik)
c          write(71,*) bbbb4(4,1,ik), bbbb4(4,2,ik),
c     >                bbbb4(4,3,ik), bbbb4(4,4,ik)
c        enddo
c
c        do ik = 1, ikmax
c          write(71,*) 'C(ik) , ik = ', ik
c          write(71,*) cccc4(1,1,ik), cccc4(1,2,ik),
c     >                cccc4(1,3,ik), cccc4(1,4,ik)
c          write(71,*) cccc4(2,1,ik), cccc4(2,2,ik),
c     >                cccc4(2,3,ik), cccc4(2,4,ik)
c          write(71,*) cccc4(3,1,ik), cccc4(3,2,ik),
c     >                cccc4(3,3,ik), cccc4(3,4,ik)
c          write(71,*) cccc4(4,1,ik), cccc4(4,2,ik),
c     >                cccc4(4,3,ik), cccc4(4,4,ik)
c        enddo
c
        call blocktridiag(aaaa4,bbbb4,cccc4,ffff4,uuuu4,ikmax,4)
c
        do in = 1, 4
          do ik = 1, ikmax
               U(in,ik) = U0(in,ik) + uuuu4(in,ik)
          enddo
        enddo
c
      else
        write(71,*) 'INCORRECT UPDATE'
        STOP
      endif
c
      do ik = 0, ikmax + 1
         if (U(1,ik).le.0.0D0) then
           U(1,ik) = abs(U0(1,ik))
         endif
c
         if (U(3,ik).le.0.0D0) then
           U(3,ik) = abs(U0(3,ik))
         endif
c
         if (U(4,ik).le.0.0D0) then
           U(4,ik) = abs(U0(4,ik))
         endif
      enddo
c
c      do ik = 0, ikmax+1
c         W(1,ik) = U(1,ik) / mass
c         W(2,ik) = U(2,ik) / U(1,ik)
c         W(3,ik) = (g_pv - 1.0D0) *
c     >        (U(3,ik) - 0.5D0 * U(1,ik)* W(2,ik) **2.0D0)
c         dummy = W(3,ik)/(a_Temp*W(1,ik)*eV)
c         if (dummy.le.0.1D0) then
c            U(3,ik) = U(3,ik) * (0.1D0/dummy)
c         endif
c      enddo
c
      return
      end
c
c
      subroutine updateBC
      implicit none
      real*8 a_last, a_orig, interpol2p, extrapol2p,
     >   W2_0,  T_0, T_1, T_2, ka_ave, delta_0, delta_N, vel1,
     >   flux1, flux2, source, zhi, c_s_0, U1_0, M_0, U3_0,
     >   c_s_N, U4_0, U1_N, U2_N, U3_N, T_N, T_N1, v_1,
     >   dx_N, dx_0, d0_1, d0_2, d0_3, dN_1, dN_2, dN_3,
     >   aa, bb, cc, dd, U2_0, flux0, g_net_0, c_a_0, c_a_N,
     >   mach0, machN, del_m, dadM, dbdM, dcdM, dddM,
     >   del_q, dadq, dbdq, dcdq, g_net,
     >   del_v, dadv, dbdv, dcdv, dN_4, U4_N, d0_4,
     >   del_U1, del_U2, del_U3, del_U4,
     >   temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, d1dx1
      integer in, ik, j, j_max
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      if (a_2T.eq.0) then
c
        c_s_0 = dsqrt( W(3,0) / U(1,0))
        c_a_0 = dsqrt(g_pv * W(3,0) / U(1,0))
c
        T_0 = W(3,0) / (W(1,0) * a_Temp)
        T_1 = W(3,1) / (W(1,1) * a_Temp)
        ka_ave = (ka(0) + ka(1)) / 2.0D0
        v_1 = abs((W(2,1) + W(2,0)) / 2.0D0)
        v_0 = W(2,0)
        dx_0 = xc(1)
        zhi  = (gamma_sum + 0.5D0 * (v_0/c_s_0)**2.0D0) *
     >          dabs(U(2,0)/mass) * a_Temp
        source = - Q(3,0) * dx_0 / 2.0D0
        flux1 = ka_ave / dx_0 - 0.5D0* v_1 * 2.5D0 * W(1,0) *
     >         a_Temp + zhi
        flux2 = - ka_ave / dx_0 * T_1 -
     >   0.5D0* v_1 * (0.5D0 * U(1,0) * v_0 * v_0 +
     >                   W(3,1) + U(3,1))
c
        T_0 = - (source + flux2) / flux1
        p_0 = W(1,0) * T_0 * a_Temp
c
        U3_0 = p_0 / (g_pv - 1.0D0) +
     >              0.5D0 * U(1,0) * v_0 * v_0
c
c      if (dabs(W(2,0)).le.c_s_0) then
c        U1_0 = dabs(U(2,0) / c_s_0)
c      else
c        U1_0 = extrapol2p(U(1,1), U(1,2), xc(1), xc(2), xc(0))
c      endif
c
        dadv = 0.0D0
        dbdv = 1.0D0
        dcdv = 0.0D0
        del_v = - MASS_0 - U(2,0)
c
        c_s_0 = dsqrt( W(3,0) / U(1,0))
        c_a_0 = dsqrt(g_pv * W(3,0) / U(1,0))
        mach0 = dabs(W(2,0) / c_s_0)
c
c        dadM = - 1.0D0 / ((g_pv - 1.0D0) * U(3,0) * mach0 ** 3.0D0 /
c     >          (2.0D0 * U(2,0) ** 2.0D0))
c        dbdM =   1.0D0 / ((g_pv - 1.0D0) * U(1,0) * U(3,0)
c     >          * mach0 ** 3.0D0 / (U(2,0) ** 3.0D0))
c        dcdM = - 1.0D0 / ((g_pv - 1.0D0) * U(1,0) * mach0 ** 3.0D0 /
c     >          (2.0D0 * U(2,0) ** 2.0D0))
c
        dadM = - 1.0D0 / ((g_pv - 1.0D0) * mach0 ** 3.0D0 *
     >           (1.0D0 / U(1,0) - U(3,0) / U(2,0) ** 2.0D0) +
     >           mach0 / U(1,0))
        dbdM = 1.0D0 / ((g_pv - 1.0D0) * mach0 ** 3.0D0 / U(2,0) +
     >                  mach0 / U(2,0))
        dcdM = - 1.0D0 / ((g_pv - 1.0D0) * mach0 ** 3.0D0 * U(1,0) /
     >                   U(2,0) ** 2.0D0)
c
        del_m = max(0.0D0, 1.0D0 - mach0)
c
        v_0    = W(2,0)
        T_0    = W(3,0) / (W(1,0) * a_Temp)
        g_net  = gamma_sum + 0.5D0 * (v_0/c_s_0)**2.0D0
        dadq = 1.0D0 / (g_net * (g_pv - 1.0D0) *
     >       ( (U(2,0)/U(1,0)) ** 3.0 -
     >         (U(2,0)*U(3,0)) / U(1,0) ** 2.0 ) )
        dbdq = 1.0D0 / (g_net * (g_pv - 1.0D0) *
     >       ( (U(3,0)/U(1,0)) -
     >         1.5 * (U(2,0) / U(1,0)) ** 2.0 ) )
        dcdq = 1.0D0 / (g_net * (g_pv - 1.0D0) *
     >         (U(2,0) / U(1,0)) )
c
        del_q = - (POW_BC_0 -
     >           dabs(g_net * W(3,0) * W(2,0)))
c
        delta_0 = 1.0D0
        d0_1 = delta_0
        d0_2 = delta_0
        d0_3 = delta_0
c
        if (mach0.lt.1.0) then
          delta_0 = delta_0 / 10.0
        endif
c
        U1_0 = extrapol2p(U(1,1), U(1,2), xc(1), xc(2), xc(0))
        U2_0 = extrapol2p(U(2,1), U(2,2), xc(1), xc(2), xc(0))
c
        del_U1 = del_m * dadm
        del_U2 = del_m * dbdm
        del_U3 = del_m * 0.0D0
c
        U(1,0) = U0(1,0) + del_U1 * delta_0
     >         + (U1_0 - U0(1,0)) * d0_1
        U(2,0) = U0(2,0) + del_U2 * delta_0
     >         + (U2_0 - U0(2,0)) * d0_2
        U(3,0) = U0(3,0) + del_U3 * delta_0
     >         + (U3_0 - U0(3,0)) * d0_3
c
        if (U(1,0).le.0.0D0)  U(1,0) = abs(U0(1,0))
        if (U(3,0).le.0.0D0)  U(3,0) = abs(U0(3,0))
c
c     SYMMETRIC BC
c
c      if (a_symm.eq.1) then
c
c        U(1,ikmax+1) =   U(1,0)
c        U(2,ikmax+1) = - U(2,0)
c        U(3,ikmax+1) =   U(3,0)
c
c     NOT SYMMETRIC BC
c
c      elseif (a_symm.eq.0) then
c
        c_s_N = dsqrt( W(3,ikmax+1) / U(1,ikmax+1))
        c_a_N = dsqrt(g_pv * W(3,ikmax+1) / U(1,ikmax+1))
c
        T_N = W(3,ikmax+1) / (W(1,ikmax+1) * a_Temp)
        T_1 = W(3,ikmax) / (W(1,ikmax) * a_Temp)
        ka_ave = (ka(ikmax+1) + ka(ikmax)) / 2.0D0
        v_1 = abs((W(2,ikmax+1) + W(2,ikmax)) / 2.0D0)
        v_N = W(2,ikmax+1)
        dx_N = x_max - xc(ikmax)
        zhi  = (gamma_sum + 0.5D0 * (v_N/c_s_N)**2.0D0) *
     >          dabs(U(2,ikmax+1)/mass) * a_Temp
        source = - Q(3,ikmax+1) * dx_N / 2.0D0
        flux1 = ka_ave / dx_N - 0.5D0* v_1 * 2.5D0* W(1,ikmax+1) *
     >         a_Temp + zhi
        flux2 = - ka_ave / dx_N * T_1 -
     >   0.5D0* v_1 * (0.5D0 * U(1,ikmax+1) * v_N * v_N +
     >                   W(3,ikmax) + U(3,ikmax))
        T_N = - (source + flux2) / flux1
        p_N = W(1,ikmax+1) * T_N * a_Temp
c
        U3_N = p_N / (g_pv - 1.0D0) +
     >               0.5D0 * U(1,ikmax+1) * v_N * v_N
c
        dadv = 0.0D0
        dbdv = 1.0D0
        dcdv = 0.0D0
        del_v = MASS_N - U(2,ikmax+1)
c
        c_s_N = dsqrt( W(3,ikmax+1) / U(1,ikmax+1))
        c_a_N = dsqrt(g_pv * W(3,ikmax+1) / U(1,ikmax+1))
        machN = dabs(W(2,ikmax+1) / c_s_N)
c
c        dadM = - 1.0D0 / ((g_pv - 1.0D0) * U(3,ikmax+1) *
c     >       machN ** 3.0D0 / (2.0D0 * U(2,ikmax+1) ** 2.0D0))
c        dbdM =   1.0D0 / ((g_pv - 1.0D0) * U(1,ikmax+1) *
c     >    U(3,ikmax+1) * machN ** 3.0D0 / (U(2,ikmax+1) ** 3.0D0))
c        dcdM = - 1.0D0 / ((g_pv - 1.0D0) * U(1,ikmax+1) *
c     >       machN ** 3.0D0 / (2.0D0 * U(2,ikmax+1) ** 2.0D0))
c
        dadM = - 1.0D0 / ((g_pv - 1.0D0) * machN ** 3.0D0 *
     >           (1.0D0 / U(1,ikmax+1) - U(3,ikmax+1) / U(2,ikmax+1)
     >            ** 2.0D0) + machN / U(1,ikmax+1))
        dbdM = 1.0D0 / ((g_pv - 1.0D0) * machN ** 3.0D0 / U(2,ikmax+1)
     >                  + machN / U(2,ikmax+1))
        dcdM = - 1.0D0 / ((g_pv - 1.0D0) * machN ** 3.0D0 *
     >                   U(1,ikmax+1) / U(2,ikmax+1) ** 2.0D0)
        del_m = max(0.0D0, 1.0D0 - machN)
c
        v_N    = W(2,ikmax+1)
        T_N    = W(3,ikmax+1) / (W(1,ikmax+1) * a_Temp)
        g_net  = gamma_sum + 0.5D0 * (v_N / c_s_N)**2.0D0
        dadq =  - 1.0D0 / (g_net * (g_pv - 1.0D0) *
     >       ( (U(2,ikmax+1)/U(1,ikmax+1)) ** 3.0 -
     >         (U(2,ikmax+1)*U(3,ikmax+1))
     >          / U(1,ikmax+1) ** 2.0 ) )
        dbdq = - 1.0D0 / (g_net * (g_pv - 1.0D0) *
     >       ( (U(3,ikmax+1)/U(1,ikmax+1)) -
     >         1.5 * (U(2,ikmax+1) / U(1,ikmax+1)) ** 2.0 ) )
        dcdq =  - 1.0D0 / (g_net * (g_pv - 1.0D0) *
     >         (U(2,ikmax+1) / U(1,ikmax+1)) )
c
        del_q = - (POW_BC_N -
     >           dabs(g_net * W(3,ikmax+1) * W(2,ikmax+1)))
c
        delta_N = 1.0D0
        dN_1 = delta_N
        dN_2 = delta_N
        dN_3 = delta_N
c
        if (machN.lt.1.0) then
          delta_N = delta_N / 10.0
        endif
c
        U1_N = extrapol2p(U(1,ikmax-1), U(1,ikmax), xc(ikmax-1),
     >                    xc(ikmax), xc(ikmax+1))
        U2_N = extrapol2p(U(2,ikmax-1), U(2,ikmax), xc(ikmax-1),
     >                    xc(ikmax), xc(ikmax+1))
c
        del_U1 = del_m * dadm
        del_U2 = del_m * dbdm
        del_U3 = del_m * 0.0D0
c
        U(1,ikmax+1) = U0(1,ikmax+1) + del_U1 * delta_N
     >                             + dN_1 * (U1_N - U0(1,ikmax+1))
        U(2,ikmax+1) = U0(2,ikmax+1) + del_U2 * delta_N
     >                             + dN_2 * (U2_N - U0(2,ikmax+1))
        U(3,ikmax+1) = U0(3,ikmax+1) + del_U3 * delta_N
     >                             + dN_3 * (U3_N - U0(3,ikmax+1))
c
        if (U(1,ikmax+1).le.0.0D0)  U(1,ikmax+1) = abs(U0(1,ikmax+1))
        if (U(3,ikmax+1).le.0.0D0)  U(3,ikmax+1) = abs(U0(3,ikmax+1))
c
      elseif (a_2T.eq.1) then
c
        c_s_0 = dsqrt( (W(3,0) + W(4,0)) / U(1,0))
        c_a_0 = dsqrt(g_pv * (W(3,0)+W(4,0)) / U(1,0))
        mach0 = dabs(W(2,0) / c_s_0)
c
        T_0 = W(3,0) / W(1,0)
        T_1 = W(3,1) / W(1,1)
        ka_ave = (ka_i(0) + ka_i(1)) / 2.0D0
        v_1 = (W(2,1) + W(2,0)) / 2.0D0
        v_0 = W(2,0)
        dx_0 = xc(1)
        zhi  = (gamma_i + 0.5D0 * (v_0/c_s_0) ** 2.0D0 *
     >          (1.0D0 + W(4,0)/W(3,0)) ) * (U(2,0)/mass)
        source = - (Q(3,0) + Q_B(3,0) + Q_Joule(0) + Q_ei(0))
     >           * dx_0 / 2.0D0
        flux1 = - ka_ave / dx_0 - 0.5D0 * v_1 * W(1,0) *
     >         (1.0D0 + 1.0D0 / (g_pv - 1.0D0)) + zhi
        flux2 = - ka_ave / dx_0 * T_1 +
     >   0.5D0 * v_1 * (0.5D0 * U(1,0) * v_0 * v_0 +
     >      W(3,1) + U(3,1) + 2.0D0 * G(2,1))
c
        T_0 = (source + flux2) / flux1
        p_0 = W(1,0) * T_0
        U3_0 = p_0 / (g_pv - 1.0D0) +
     >              0.5D0 * U(1,0) * v_0 * v_0
c
        T_0 = W(4,0) / W(1,0)
        T_1 = W(4,1) / W(1,1)
        T_2 = W(4,2) / W(1,2)
        ka_ave = (ka_e(0) + ka_e(1)) / 2.0D0
        v_1 = (W(2,1) + W(2,0)) / 2.0D0
        v_0 = W(2,0)
        dx_0 = xc(1)
        zhi  = gamma_e * (U(2,0)/mass)
        source = - (Q(4,0) + Q_B(4,0) - Q_Joule(0) - Q_ei(0))
     >           * dx_0 / 2.0D0
c
        if (T_0/eV.le.3.0D0) then
          j_max = 1
        elseif (T_0/eV.le.10.0D0) then
          j_max = 10
        elseif (T_0/eV.le.30.0D0) then
          j_max = 30
        elseif (T_0/eV.le.60.0D0) then
          j_max = 100
        elseif (T_0/eV.le.100.0D0) then
          j_max = 300
        else
          j_max = 1000
        endif
c
        if (T_0/eV.le.s23_par_updtbcte) then
c
          T_1    = W(4,1) / W(1,1)
          flux1 = - ka_ave / dx_0 - 0.5D0 * v_1 * 2.5D0 * W(1,0) + zhi
          flux2 = - ka_ave / dx_0 * T_1 +
     >     0.5D0 * v_1 * (2.5D0 * W(4,1))
c
          T_0 = (source + flux2) / flux1
          W(4,0) = W(1,0) * T_0
          U4_0 = W(4,0) / (g_pv - 1.0D0)
c
        else
c
          do j = 1, j_max
c
            T_1    = W(4,1) / W(1,1)
            flux1 = - ka_ave / dx_0 - 0.5D0 * v_1 * 2.5D0 * W(1,0) + zhi
            flux2 = - ka_ave / dx_0 * T_1 +
     >       0.5D0 * v_1 * (2.5D0 * W(4,1))
c
            T_0 = (source + flux2) / flux1
            W(4,0) = W(1,0) * T_0
            U4_0 = W(4,0) / (g_pv - 1.0D0)
c
            do ik = 1, 10
c
              temp_1 = 0.5 * (W(2,ik) + W(2,ik-1)) * 0.5 *
     >             2.5 * (W(4,ik) + W(4,ik-1))
              temp_2 = 0.5 * (W(2,ik+1) + W(2,ik)) * 0.5 *
     >             2.5 * (W(4,ik+1) + W(4,ik))
              temp_3 = - 0.5 * (ka_e(ik) + ka_e(ik-1)) *
     >            (W(4,ik)/W(1,ik) - W(4,ik-1)/W(1,ik-1))
     >            / (xc(ik) - xc(ik-1))
              temp_4 = - 0.5 * (ka_e(ik+1) + ka_e(ik)) *
     >            (W(4,ik+1)/W(1,ik+1) - W(4,ik)/W(1,ik))
     >            / (xc(ik+1) - xc(ik))
              temp_5 = RH(4,ik) + Q(4,ik) + Q_B(4,ik) -
     >               Q_ei(ik) - Q_joule(ik)
              temp_6 = 0.5 * (xf(ik+1) - xf(ik)) ** 2.0
     >               / (2.0 * ka_e(ik) / (W(1,ik) * c_v_0))
              R(4,ik) = temp_5
     >                - d1dx1(temp_3, temp_4, (xf(ik+1) - xf(ik)))
     >                - d1dx1(temp_1, temp_2, (xf(ik+1) - xf(ik)))
c
              U(4,ik) = U(4,ik) + R(4,ik) * temp_6
              W(4,ik) = (g_pv - 1.0D0) * U(4,ik)
c
            enddo
c
          enddo
c
        endif
c
        dadM = - 1.0D0 / ((g_pv - 1.0D0) * 0.5D0 * mach0 ** 3.0D0 *
     >           (1.0D0 / U(1,0) - (U(3,0) + U(4,0)) /
     >           U(2,0) ** 2.0D0) + mach0 / U(1,0))
        dbdM = 1.0D0 / ((g_pv - 1.0D0) * 0.5D0 *
     >           mach0 ** 3.0D0 / U(2,0) + mach0 / U(2,0))
        dcdM = - 1.0D0 / ((g_pv - 1.0D0) * 0.5D0 *
     >           mach0 ** 3.0D0 * U(1,0) / U(2,0) ** 2.0D0)
        dddM = dcdM
c
        del_m = max(0.0D0, 1.0D0 - mach0)
c
        delta_0 = s23_par_updtdel0
c
        d0_1 = s23_par_updtdel1
        d0_2 = s23_par_updtdel1
        d0_3 = s23_par_updtdel1
        d0_4 = s23_par_updtdel1
c
        if (mach0.lt.1.0) then
          delta_0 = delta_0 * s23_par_updtdelm
        endif
c
        U1_0 = extrapol2p(U(1,1), U(1,2), xc(1), xc(2), xc(0))
        U2_0 = extrapol2p(U(2,1), U(2,2), xc(1), xc(2), xc(0))
c
        del_U1 = del_m * dadm
        del_U2 = del_m * dbdm
        del_U3 = del_m * 0.0D0
        del_U4 = del_m * 0.0D0
c
        U(1,0) = U0(1,0) + del_U1 * delta_0
     >         + (U1_0 - U0(1,0)) * d0_1
        U(2,0) = U0(2,0) + del_U2 * delta_0
     >         + (U2_0 - U0(2,0)) * d0_2
        U(3,0) = U0(3,0) + del_U3 * delta_0
     >         + (U3_0 - U0(3,0)) * d0_3
        U(4,0) = U0(4,0) + del_U4 * delta_0
     >         + (U4_0 - U0(4,0)) * d0_4
c
        if (U(1,0).le.0.0D0)  U(1,0) = abs(U0(1,0))
        if (U(3,0).le.0.0D0)  U(3,0) = abs(U0(3,0))
        if (U(4,0).le.0.0D0)  U(4,0) = abs(U0(4,0))
c
c
        c_s_N = dsqrt( (W(3,ikmax+1) + W(4,ikmax+1))
     >          / U(1,ikmax+1))
        c_a_N = dsqrt(g_pv * (W(3,ikmax+1) + W(4,ikmax+1))
     >          / U(1,ikmax+1))
        machN = dabs(W(2,ikmax+1) / c_s_N)
c
        T_N = W(3,ikmax+1) / W(1,ikmax+1)
        T_1 = W(3,ikmax) / W(1,ikmax)
        ka_ave = (ka_i(ikmax+1) + ka_i(ikmax)) / 2.0D0
        v_1 = (W(2,ikmax+1) + W(2,ikmax)) / 2.0D0
        v_N = W(2,ikmax+1)
        dx_N = x_max - xc(ikmax)
        zhi  = (gamma_i + 0.5D0 * (v_N/c_s_N)**2.0D0 *
     >           (1.0D0 + W(4,ikmax+1) / W(3,ikmax+1)) )
     >         * dabs(U(2,ikmax+1)/mass)
        source = (Q(3,ikmax+1) + Q_B(3,ikmax+1) +
     >           Q_Joule(ikmax+1) + Q_ei(ikmax+1)) * dx_N / 2.0D0
        flux1 = ka_ave / dx_N - 0.5D0* v_1 * W(1,ikmax+1) *
     >         (1.0D0 + 1.0D0 / (g_pv - 1.0D0)) + zhi
        flux2 = ka_ave / dx_N * T_1 +
     >   0.5D0 * v_1 * (0.5D0 * U(1,ikmax+1) * v_N * v_N +
     >     W(3,ikmax) + U(3,ikmax) + 2.0D0 * G(2,ikmax+1))
        T_N = (source + flux2) / flux1
        p_N = W(1,ikmax+1) * T_N
        U3_N = p_N / (g_pv - 1.0D0) +
     >               0.5D0 * U(1,ikmax+1) * v_N * v_N
c
        T_N = W(4,ikmax+1) / W(1,ikmax+1)
        T_1 = W(4,ikmax) / W(1,ikmax)
        ka_ave = (ka_e(ikmax+1) + ka_e(ikmax)) / 2.0D0
        v_1 = (W(2,ikmax+1) + W(2,ikmax)) / 2.0D0
        v_N = W(2,ikmax+1)
        dx_N = x_max - xc(ikmax)
        zhi  = gamma_e * (U(2,ikmax+1)/mass)
        source = (Q(4,ikmax+1) + Q_B(4,ikmax+1) -
     >           Q_Joule(ikmax+1) - Q_ei(ikmax+1)) * dx_N / 2.0D0
c
        if (T_N/eV.le.3.0D0) then
          j_max = 1
        elseif (T_N/eV.le.10.0D0) then
          j_max = 10
        elseif (T_N/eV.le.30.0D0) then
          j_max = 30
        elseif (T_N/eV.le.60.0D0) then
          j_max = 100
        elseif (T_N/eV.le.100.0D0) then
          j_max = 300
        else
          j_max = 1000
        endif
c
        if (T_N/eV.le.s23_par_updtbcte) then
c
          T_1 = W(4,ikmax) / W(1,ikmax)
          flux1 = ka_ave / dx_N
     >          - 0.5D0* v_1 * 2.5 * W(1,ikmax+1) + zhi
          flux2 = ka_ave / dx_N * T_1 +
     >            0.5D0 * v_1 * 2.5 * W(4,ikmax)
          T_N = (source + flux2) / flux1
          W(4,ikmax+1) = W(1,ikmax+1) * T_N
          U4_N = W(4,ikmax+1) / (g_pv - 1.0D0)
c
        else
c
          do j = 1, j_max
c
            T_1 = W(4,ikmax) / W(1,ikmax)
            flux1 = ka_ave / dx_N
     >          - 0.5D0* v_1 * 2.5 * W(1,ikmax+1) + zhi
            flux2 = ka_ave / dx_N * T_1 +
     >            0.5D0 * v_1 * 2.5 * W(4,ikmax)
            T_N = (source + flux2) / flux1
            W(4,ikmax+1) = W(1,ikmax+1) * T_N
            U4_N = W(4,ikmax+1) / (g_pv - 1.0D0)
c
            do ik = ikmax, ikmax - 10, -1
c
              temp_1 = 0.5 * (W(2,ik) + W(2,ik-1)) * 0.5 *
     >             2.5 * (W(4,ik) + W(4,ik-1))
              temp_2 = 0.5 * (W(2,ik+1) + W(2,ik)) * 0.5 *
     >             2.5 * (W(4,ik+1) + W(4,ik))
              temp_3 = - 0.5 * (ka_e(ik) + ka_e(ik-1)) *
     >            (W(4,ik)/W(1,ik) - W(4,ik-1)/W(1,ik-1))
     >            / (xc(ik) - xc(ik-1))
              temp_4 = - 0.5 * (ka_e(ik+1) + ka_e(ik)) *
     >            (W(4,ik+1)/W(1,ik+1) - W(4,ik)/W(1,ik))
     >            / (xc(ik+1) - xc(ik))
              temp_5 = RH(4,ik) + Q(4,ik) + Q_B(4,ik) -
     >               Q_ei(ik) - Q_joule(ik)
              temp_6 = 0.5 * (xf(ik+1) - xf(ik)) ** 2.0
     >               / (2.0 * ka_e(ik) / (W(1,ik) * c_v_0))
              R(4,ik) = temp_5
     >                - d1dx1(temp_3, temp_4, (xf(ik+1) - xf(ik)))
     >                - d1dx1(temp_1, temp_2, (xf(ik+1) - xf(ik)))
c
              U(4,ik) = U(4,ik) + R(4,ik) * temp_6
              W(4,ik) = (g_pv - 1.0D0) * U(4,ik)
c
            enddo
c
          enddo
c
        endif
c
        dadM = - 1.0D0 / ((g_pv - 1.0D0) * machN ** 3.0D0 *
     >           (1.0D0 / U(1,ikmax+1) - U(3,ikmax+1) / U(2,ikmax+1)
     >            ** 2.0D0) + machN / U(1,ikmax+1))
        dbdM = 1.0D0 / ((g_pv - 1.0D0) * machN ** 3.0D0 / U(2,ikmax+1)
     >                  + machN / U(2,ikmax+1))
        dcdM = - 1.0D0 / ((g_pv - 1.0D0) * machN ** 3.0D0 *
     >                   U(1,ikmax+1) / U(2,ikmax+1) ** 2.0D0)
        dddM = dcdM
        del_m = max(0.0D0, 1.0D0 - machN)
c
        delta_N = s23_par_updtdel0
c
        dN_1 = s23_par_updtdel1
        dN_2 = s23_par_updtdel1
        dN_3 = s23_par_updtdel1
        dN_4 = s23_par_updtdel1
c
        if (machN.lt.1.0) then
          delta_N = delta_N * s23_par_updtdelm
        endif
c
        U1_N = extrapol2p(U(1,ikmax-1), U(1,ikmax), xc(ikmax-1),
     >                    xc(ikmax), xc(ikmax+1))
        U2_N = extrapol2p(U(2,ikmax-1), U(2,ikmax), xc(ikmax-1),
     >                    xc(ikmax), xc(ikmax+1))
c
        del_U1 = del_m * dadm
        del_U2 = del_m * dbdm
        del_U3 = del_m * 0.0D0
        del_U4 = del_m * 0.0D0
c
        U(1,ikmax+1) = U0(1,ikmax+1) + del_U1 * delta_N
     >                             + dN_1 * (U1_N - U0(1,ikmax+1))
        U(2,ikmax+1) = U0(2,ikmax+1) + del_U2 * delta_N
     >                             + dN_2 * (U2_N - U0(2,ikmax+1))
        U(3,ikmax+1) = U0(3,ikmax+1) + del_U3 * delta_N
     >                             + dN_3 * (U3_N - U0(3,ikmax+1))
        U(4,ikmax+1) = U0(4,ikmax+1) + del_U4 * delta_N
     >                             + dN_4 * (U4_N - U0(4,ikmax+1))
c
        if (U(1,ikmax+1).le.0.0D0)  U(1,ikmax+1) = abs(U0(1,ikmax+1))
        if (U(3,ikmax+1).le.0.0D0)  U(3,ikmax+1) = abs(U0(3,ikmax+1))
        if (U(4,ikmax+1).le.0.0D0)  U(4,ikmax+1) = abs(U0(4,ikmax+1))
c
      else
        write(71,*) 'UPDATEBC INCORRECT'
        STOP
      endif
c
      return
      end
c
c
      subroutine zeroR
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
c
      do ik = 1, ikmax+1
         do in = 1, 4
            RF(in,ik) = 0.0D0
            RG(in,ik) = 0.0D0
            RH(in,ik) = 0.0D0
            R(in,ik)  = 0.0D0
         enddo
      enddo
      max_res = 0.0D0
      rms_res = 0.0D0
c
      return
      end
c
c
      subroutine fluxes
      implicit none
c
c               Hj-1          Hj         Hj+1
c               Gj-1          Gj         Gj+1
c               Fj-1          Fj         Fj+1
c
c         o      |     o      |     o     |     o
c        j-2          j-1           j          j+1
c
c         A            B            C           D
c
      include 'params'
      include 'sol23_com'
      include 'cfd_osm_com'
      integer in,ik
      real*8 dx_bc, rho_b, rho_c, rho_ave,
     >       qm_b, qm_c, qm_ave,
     >       v_b, v_c, v_ave, TJ_b, TJ_c,
     >       p_b, p_c, p_ave, T_ave, c_s,
     >       qE_b, qE_c, qE_ave, ka_ave, mu_ave,
     >       c_s_b, c_s_c, c_s_ave, F_b, F_c,
     >       netv_b, netv_c, eps4, eps2, netv_ave,
     >       eps2_b, eps2_c, aa, bb, cc, dd, Q_temp
      real*8 ave3p, d3dx3, d2dx2, d1dx1, extrapol2p
      external ave3p, d3dx3, d2dx2, d1dx1, extrapol2p
c
      do ik = 1, ikmax + 1
         dx_bc = xc(ik) - xc(ik-1)
         F(1,ik) = 0.5D0 *(U(2,ik-1) + U(2,ik))
c
c         F(2,ik) = (2.0D0 * U(3,ik) + cons_gpv * W(3,ik) +
c     >       2.0D0 * U(3,ik-1) + cons_gpv * W(3,ik-1)) / 2.0D0
c
         if (a_2T.eq.0) then
           F_b = W(3,ik-1) + U(2,ik-1) ** 2.0D0 / U(1,ik-1)
           F_c = W(3,ik)   + U(2,ik) ** 2.0D0 / U(1,ik)
           F(2,ik) = (F_b + F_c) / 2.0D0
         else
           F_b = W(3,ik-1) + W(4,ik-1) +
     >                       U(2,ik-1) ** 2.0D0 / U(1,ik-1)
           F_c = W(3,ik) + W(4,ik) + U(2,ik) ** 2.0D0 / U(1,ik)
           F(2,ik) = (F_b + F_c) / 2.0D0
         endif
c
         F(3,ik) = 0.5D0 * (W(2,ik-1) + W(2,ik)) * 0.5D0 *
     >      ((W(3,ik) + U(3,ik)) + (W(3,ik-1) + U(3,ik-1)))
         F(4,ik) = 0.5D0 * (W(2,ik-1) + W(2,ik)) * 0.5D0 *
     >      ((W(4,ik) + U(4,ik)) + (W(4,ik-1) + U(4,ik-1)))
c
         G(1,ik) = 0.0D0
         mu_ave = (mu(ik-1) + mu(ik)) / 2.0D0
         G(2,ik) = - mu_ave * d1dx1(W(2,ik-1), W(2,ik), dx_bc)
c
         if (a_2T.eq.0) then
           TJ_b = W(3,ik-1) / (W(1,ik-1) * a_Temp)
           TJ_c = W(3,ik)   / (W(1,ik) * a_Temp)
           ka_ave = (ka(ik-1) + ka(ik)) / 2.0D0
           G(3,ik) = - ka_ave * d1dx1(TJ_b, TJ_c, dx_bc)
     >               + 0.5D0 * (W(2,ik-1) + W(2,ik)) * G(2,ik)
           TJ_b = W(4,ik-1) / (W(1,ik-1) * a_Temp)
           TJ_c = W(4,ik)   / (W(1,ik) * a_Temp)
           ka_ave = (ka(ik-1) + ka(ik)) / 2.0D0
           G(4,ik) = - ka_ave * d1dx1(TJ_b, TJ_c, dx_bc)
         else
           TJ_b = W(3,ik-1) / (W(1,ik-1) * a_Temp)
           TJ_c = W(3,ik)   / (W(1,ik) * a_Temp)
           ka_ave = (ka_i(ik-1) + ka_i(ik)) / 2.0D0
           G(3,ik) = - ka_ave * d1dx1(TJ_b, TJ_c, dx_bc)
     >               + 0.5D0 * (W(2,ik-1) + W(2,ik)) * G(2,ik)
           TJ_b = W(4,ik-1) / (W(1,ik-1) * a_Temp)
           TJ_c = W(4,ik)   / (W(1,ik) * a_Temp)
           ka_ave = (ka_e(ik-1) + ka_e(ik)) / 2.0D0
           G(4,ik) = - ka_ave * d1dx1(TJ_b, TJ_c, dx_bc)
         endif
c
c         F(1,ik) = (U(2,ik) + U(2,ik-1)) / 2.0D0
c         F_b = (g_pv - 1.0D0) * U(3,ik-1) +
c     >    0.5D0 * (3.0D0 - g_pv) * U(2,ik-1) ** 2.0D0 / U(1,ik-1)
c         F_c = (g_pv - 1.0D0) * U(3,ik) +
c     >    0.5D0 * (3.0D0 - g_pv) * U(2,ik) ** 2.0D0 / U(1,ik)
c         F(2,ik) = (F_b + F_c) / 2.0D0
c         F_b = (g_pv * U(3,ik-1) - 0.5D0 *
c     >         (g_pv - 1.0D0) * U(2,ik-1) ** 2.0D0 / U(1,ik-1))
c     >        * U(2,ik-1) / U(1,ik-1)
c         F_c = (g_pv * U(3,ik) - 0.5D0 *
c     >         (g_pv - 1.0D0) * U(2,ik) ** 2.0D0 / U(1,ik))
c     >        * U(2,ik) / U(1,ik)
c         F(3,ik) = (F_b + F_c) / 2.0D0
c
c         c_s_b = dsqrt(g_pv * W(3,ik-1) / U(1,ik-1))
c         c_s_c = dsqrt(g_pv * W(3,ik) / U(1,ik))
c         c_s_b = c_s_j(ik-1)
c         c_s_c = c_s_j(ik)
c         netv_ave = 0.5D0*(c_s_b + c_s_c) + dabs(W(2,ik))
c         eps4 = a_artvisc4 * netv_ave
c
         eps2 = (ep2(ik) + ep2(ik-1)) / 2.0D0
         eps4 = (ep4(ik) + ep4(ik-1)) / 2.0D0
c
         do in = 1, 4
             if (ik.eq.1) then
               bb = 1.0D0
               cc = (xc(ik+1) - xc(ik)) / (xc(ik) - xc(ik-1))
               H(in,ik) = - eps4 * d2dx2(U(in,ik-1), U(in,ik),
     >                     U(in,ik+1), bb, cc)
     >                  + eps2 * d1dx1(U(in,ik-1), U(in,ik), bb)
             elseif (ik.eq.ikmax+1) then
               aa = (xc(ik-1) - xc(ik-2)) / (xc(ik) - xc(ik-1))
               bb = 1.0D0
               H(in,ik) = eps4 * d2dx2(U(in,ik-2), U(in,ik-1),
     >                      U(in,ik), aa, bb)
     >                  + eps2 * d1dx1(U(in,ik-1), U(in,ik), bb)
             else
               aa = (xc(ik-1) - xc(ik-2)) / (xc(ik) - xc(ik-1))
               bb = 1.0D0
               cc = (xc(ik+1) - xc(ik)) / (xc(ik) - xc(ik-1))
               H(in,ik) = - eps4 * d3dx3(U(in,ik-2), U(in,ik-1),
     >         U(in,ik), U(in,ik+1), aa, bb, cc)
     >                  + eps2 * d1dx1(U(in,ik-1), U(in,ik), bb)
             endif
         enddo
c
      enddo
c
c      H(4,1)     = 0.0D0
c      H(4,ikmax) = 0.0D0
c
      F(1,0) = U(2,0)
      if (a_2T.eq.0) then
        F(2,0) = W(3,0) + U(2,0) ** 2.0D0 / U(1,0)
      else
        F(2,0) = W(3,0) + W(4,0) + U(2,0) ** 2.0D0 / U(1,0)
      endif
c
c      F(2,0) = 2.0D0 * U(3,0) + cons_gpv * W(3,0)
c
      F(3,0) = W(2,0) * (W(3,0) + U(3,0))
      F(4,0) = W(2,0) * (W(4,0) + U(4,0))
      G(2,0) = F(2,1) + G(2,1) - F(2,0) - Q(2,0) * xf(1)
c
c      G(3,0) = F(3,1) + G(3,1) - F(3,0) - Q(3,0) * xf(1)
c      G(4,0) = F(4,1) + G(4,1) - F(4,0) - Q(4,0) * xf(1)
c
c
      G(3,0) = F(3,1) + G(3,1) - F(3,0)
     >         - (Q(3,0) + Q_B(3,0) + Q_Joule(0) + Q_ei(0) ) * xf(1)
      G(4,0) = F(4,1) + G(4,1) - F(4,0)
     >         - (Q(4,0) + Q_B(4,0) - Q_Joule(0) - Q_ei(0) ) * xf(1)
c
c      write(71,*) F(3,0), F(3,1), G(3,0), G(3,1)
c      STOP
c
      F(1,ikmax+2) = U(2,ikmax+1)
      if (a_2T.eq.0) then
        F(2,ikmax+2) = W(3,ikmax+1) +
     >                 U(2,ikmax+1) ** 2.0D0 / U(1,ikmax+1)
      else
        F(2,ikmax+2) = W(3,ikmax+1) + W(4,ikmax+1) +
     >                 U(2,ikmax+1) ** 2.0D0 / U(1,ikmax+1)
      endif
c
c      F(2,ikmax+2) = 2.0D0 * U(3,ikmax+1)
c     >               + cons_gpv * W(3,ikmax+1)
c
      F(3,ikmax+2) = W(2,ikmax+1) *
     >               (W(3,ikmax+1) + U(3,ikmax+1))
      F(4,ikmax+2) = W(2,ikmax+1) *
     >               (W(4,ikmax+1) + U(4,ikmax+1))
c
      G(2,ikmax+2) = F(2,ikmax+1) + G(2,ikmax+1)
     >                - F(2,ikmax+2) + Q(2,ikmax+1)
     >      * (xf(ikmax+2) - xf(ikmax+1))
      Q_temp = Q(3,ikmax+1) + Q_B(3,ikmax+1) + Q_Joule(ikmax+1)
     >         + Q_ei(ikmax+1)
      G(3,ikmax+2) = F(3,ikmax+1) + G(3,ikmax+1)
     >                - F(3,ikmax+2) + Q_temp
     >      * (xf(ikmax+2) - xf(ikmax+1))
      Q_temp = Q(4,ikmax+1) + Q_B(4,ikmax+1) - Q_Joule(ikmax+1)
     >         - Q_ei(ikmax+1)
      G(4,ikmax+2) = F(4,ikmax+1) + G(4,ikmax+1)
     >                - F(4,ikmax+2) + Q_temp
     >      * (xf(ikmax+2) - xf(ikmax+1))
c
c      F(1,ikmax+2) = - F(1,0)
c      F(2,ikmax+2) =   F(2,0)
c      F(3,ikmax+2) = - F(3,0)
c      G(2,ikmax+2) =   G(2,0)
c      G(3,ikmax+2) = - G(3,0)
c
      if (a_flag.eq.0.and.a_2T.eq.0) then
        do ik = 0, ikmax+1
          jac_A3(1,1,ik) = 0.0D0
          jac_A3(1,2,ik) = 1.0D0
          jac_A3(1,3,ik) = 0.0D0
          jac_A3(2,1,ik) = (0.5D0 * (g_pv - 1.0D0) - 1.0D0) *
     >                 W(2,ik) ** 2.0D0
          jac_A3(2,2,ik) = (3.0D0 - g_pv) * W(2,ik)
          jac_A3(2,3,ik) = g_pv - 1.0D0
          jac_A3(3,1,ik) = W(2,ik) * ((g_pv - 1.0D0) * W(2,ik) ** 2.0D0
     >                 - g_pv * U(3,ik) / U(1,ik))
          jac_A3(3,2,ik) = g_pv * U(3,ik) / U(1,ik)
     >                 - 1.5D0 * (g_pv - 1.0D0) * W(2,ik) ** 2.0D0
          jac_A3(3,3,ik) = g_pv * W(2,ik)
c
          jac_B3(1,1,ik) = 0.0D0
          jac_B3(1,2,ik) = 0.0D0
          jac_B3(1,3,ik) = 0.0D0
          jac_B3(2,1,ik) = - (U(2,ik) / (U(1,ik)**2.0D0))
          jac_B3(2,2,ik) = 1.0D0 / U(1,ik)
          jac_B3(2,3,ik) = 0.0D0
          jac_B3(3,1,ik) = (- U(3,ik) + U(2,ik)**2.0D0/U(1,ik) )
     >                    / (U(1,ik)**2.0D0 * c_v)
          jac_B3(3,2,ik) = - U(2,ik) / (U(1,ik)**2.0D0 * c_v)
          jac_B3(3,3,ik) = 1.0D0 / (U(1,ik) * c_v)
        enddo
c
      elseif (a_flag.eq.1.and.a_2T.eq.0) then
c
        do ik = 0, ikmax+1
         c_s = c_s_j(ik)
         lam_mat(1,ik) = W(2,ik)
         lam_mat(2,ik) = W(2,ik) + c_s_j(ik)
         lam_mat(3,ik) = W(2,ik) - c_s_j(ik)
c
         t_mat(1,1,ik) = 1.0D0
         t_mat(1,2,ik) = U(1,ik) / (2.0D0*c_s)
         t_mat(1,3,ik) = - U(1,ik) / (2.0D0*c_s)
         t_mat(2,1,ik) = W(2,ik)
         t_mat(2,2,ik) = (W(2,ik) + c_s) * U(1,ik) / (2.0D0*c_s)
         t_mat(2,3,ik) = - (W(2,ik) - c_s) * U(1,ik) / (2.0D0*c_s)
         t_mat(3,1,ik) = W(2,ik) ** 2.0D0 / 2.0D0
         t_mat(3,2,ik) = (W(2,ik) ** 2.0D0 / 2.0D0 + W(2,ik) * c_s +
     >       c_s * c_s / (g_pv - 1.0D0)) * U(1,ik) / (2.0D0*c_s)
         t_mat(3,3,ik) = - (W(2,ik) ** 2.0D0 / 2.0D0 - W(2,ik) * c_s +
     >       c_s * c_s / (g_pv - 1.0D0)) * U(1,ik) / (2.0D0*c_s)
c
         tinv_mat(1,1,ik) = 1.0D0 -
     >            0.5D0 * (g_pv - 1.0D0) * (W(2,ik)/c_s) ** 2.0D0
         tinv_mat(1,2,ik) = (g_pv - 1.0D0) * W(2,ik) / (c_s * c_s)
         tinv_mat(1,3,ik) = - (g_pv - 1.0D0) / (c_s * c_s)
         tinv_mat(2,1,ik) = (0.5D0 * (g_pv - 1.0D0) * W(2,ik) - c_s) *
     >       W(2,ik) / (U(1,ik) * c_s)
         tinv_mat(2,2,ik) = (c_s - (g_pv - 1.0D0) * W(2,ik)) /
     >              (U(1,ik)*c_s)
         tinv_mat(2,3,ik) = (g_pv - 1.0D0) / (U(1,ik) * c_s)
         tinv_mat(3,1,ik) = - (0.5D0 * (g_pv - 1.0D0) * W(2,ik) + c_s) *
     >       W(2,ik) / (U(1,ik) * c_s)
         tinv_mat(3,2,ik) = (c_s + (g_pv - 1.0D0) * W(2,ik)) /
     >              (U(1,ik)*c_s)
         tinv_mat(3,3,ik) = - (g_pv - 1.0D0) / (U(1,ik) * c_s)
        enddo
c
      elseif (a_flag.eq.0.and.a_2T.eq.1) then
c
        do ik = 0, ikmax+1
          jac_A4(1,1,ik) = 0.0D0
          jac_A4(1,2,ik) = 1.0D0
          jac_A4(1,3,ik) = 0.0D0
          jac_A4(1,4,ik) = 0.0D0
          jac_A4(2,1,ik) = (0.5D0 * (g_pv - 1.0D0) - 1.0D0) *
     >                 W(2,ik) ** 2.0D0
          jac_A4(2,2,ik) = (3.0D0 - g_pv) * W(2,ik)
          jac_A4(2,3,ik) = g_pv - 1.0D0
          jac_A4(2,4,ik) = g_pv - 1.0D0
          jac_A4(3,1,ik) = W(2,ik) * ((g_pv - 1.0D0) * W(2,ik) ** 2.0D0
     >                 - g_pv * U(3,ik) / U(1,ik))
          jac_A4(3,2,ik) = g_pv * U(3,ik) / U(1,ik)
     >                 - 1.5D0 * (g_pv - 1.0D0) * W(2,ik) ** 2.0D0
          jac_A4(3,3,ik) = g_pv * W(2,ik)
          jac_A4(3,4,ik) = 0.0D0
          jac_A4(4,1,ik) = - g_pv * W(2,ik) * U(4,ik) / U(1,ik)
          jac_A4(4,2,ik) = g_pv * U(4,ik) / U(1,ik)
          jac_A4(4,3,ik) = 0.0D0
          jac_A4(4,4,ik) = g_pv * W(2,ik)
c
          jac_B4(1,1,ik) = 0.0D0
          jac_B4(1,2,ik) = 0.0D0
          jac_B4(1,3,ik) = 0.0D0
          jac_B4(1,4,ik) = 0.0D0
          jac_B4(2,1,ik) = - (U(2,ik) / (U(1,ik)**2.0D0))
          jac_B4(2,2,ik) = 1.0D0 / U(1,ik)
          jac_B4(2,3,ik) = 0.0D0
          jac_B4(2,4,ik) = 0.0D0
          jac_B4(3,1,ik) = (- U(3,ik) + U(2,ik)**2.0D0/U(1,ik) )
     >                    / (U(1,ik)**2.0D0 * c_v)
          jac_B4(3,2,ik) = - U(2,ik) / (U(1,ik)**2.0D0 * c_v)
          jac_B4(3,3,ik) = 1.0D0 / (U(1,ik) * c_v)
          jac_B4(3,4,ik) = 0.0D0
          jac_B4(4,1,ik) = - U(4,ik) / (U(1,ik)**2.0D0 * c_v)
          jac_B4(4,2,ik) = 0.0D0
          jac_B4(4,3,ik) = 0.0D0
          jac_B4(4,4,ik) = 1.0D0 / (U(1,ik) * c_v)
c
          jac_Q4(1,1,ik) = 0.0D0
          jac_Q4(1,2,ik) = 0.0D0
          jac_Q4(1,3,ik) = 0.0D0
          jac_Q4(1,4,ik) = 0.0D0
          jac_Q4(2,1,ik) = 0.0D0
          jac_Q4(2,2,ik) = 0.0D0
          jac_Q4(2,3,ik) = 0.0D0
          jac_Q4(2,4,ik) = 0.0D0
          jac_Q4(3,1,ik) = ((g_pv - 1.0D0) / tau_ei(ik)) *
     >                 ( - 0.5D0 * (U(2,ik) / U(1,ik)) ** 2.0D0 )
     >                     + Q_ei(ik) * 2.5D0 / U(1,ik)
          jac_Q4(3,2,ik) =  ((g_pv - 1.0D0) / tau_ei(ik)) *
     >                     ( U(2,ik) / U(1,ik) )
          jac_Q4(3,3,ik) = - (g_pv - 1.0D0) / tau_ei(ik)
          jac_Q4(3,4,ik) =   (g_pv - 1.0D0) / tau_ei(ik)
     >                     - Q_ei(ik) * 1.5D0 / U(4,ik)
          jac_Q4(4,1,ik) = - jac_Q4(3,1,ik)
          jac_Q4(4,2,ik) = - jac_Q4(3,2,ik)
          jac_Q4(4,3,ik) = - jac_Q4(3,3,ik)
          jac_Q4(4,4,ik) = - jac_Q4(3,4,ik)
c
          jac_Q4(3,1,ik) = jac_Q4(3,1,ik) - Q_joule(ik) / U(1,ik)
          jac_Q4(4,1,ik) = jac_Q4(4,1,ik) + Q_joule(ik) / U(1,ik)
c
          if (U(2,ik).ne.0.0D0) then
            jac_Q4(3,2,ik) = jac_Q4(3,2,ik) + Q_joule(ik) / U(2,ik)
            jac_Q4(4,2,ik) = jac_Q4(4,2,ik) - Q_joule(ik) / U(2,ik)
          endif
c
c          Q_B(1,ik)      = (F(1,ik) + G(1,ik) + F(1,ik+1) + G(1,ik+1))
c     >                       / (2.0D0 * L_B(ik))
c          Q_B(2,ik)      = (F(2,ik) + G(2,ik) + F(2,ik+1) + G(2,ik+1))
c     >                       / (2.0D0 * L_B(ik)) +
c     >                     (0.5D0 * (G(2,ik) + G(2,ik+1)) -
c     >                       W(3,ik) - W(4,ik)) / (2.0D0 * L_B(ik))
c          Q_B(3,ik)      = (F(3,ik) + G(3,ik) * 0.01 +
c     >                      F(3,ik+1) + G(3,ik+1))
c     >                       / (2.0D0 * L_B(ik))
c          Q_B(4,ik)      = (F(4,ik) + G(4,ik) + F(4,ik+1) + G(4,ik+1))
c     >                       / (2.0D0 * L_B(ik))
c
          Q_B(1,ik)      = U(2,ik) / L_B(ik)
          Q_B(2,ik)         = (U(2,ik) ** 2.0 / U(1,ik)) / L_B(ik)
     >             + 1.5D0 * (G(2,ik) + G(2,ik+1)) * 0.5 / L_B(ik)
          Q_B(3,ik)      = W(2,ik) * (W(3,ik) + U(3,ik)) / L_B(ik)
     >                     + (G(3,ik) + G(3,ik+1)) * 0.5 / L_B(ik)
          Q_B(4,ik)      = W(2,ik) * (W(4,ik) + U(4,ik)) / L_B(ik)
     >                     + (G(4,ik) + G(4,ik+1)) * 0.5 / L_B(ik)
c
         if (a_cfdfile.eq.0) then
           Q_B(1,ik) = Q_B(1,ik) * min(1.0D0, a_pin*s23_par_Qbrelax)
           Q_B(2,ik) = Q_B(2,ik) * min(1.0D0, a_pin*s23_par_Qbrelax)
           Q_B(3,ik) = Q_B(3,ik) * min(1.0D0, a_pin*s23_par_Qbrelax)
           Q_B(4,ik) = Q_B(4,ik) * min(1.0D0, a_pin*s23_par_Qbrelax)
         endif
c
          if (s23_par_fluxexp.eq.0) then
            Q_B(1,ik) = 0.0D0
            Q_B(2,ik) = 0.0D0
            Q_B(3,ik) = 0.0D0
            Q_B(4,ik) = 0.0D0
          endif
c
        enddo
c
        Q_B(1,0) = 0.0D0
        Q_B(2,0) = 0.0D0
        Q_B(3,0) = 0.0D0
        Q_B(4,0) = 0.0D0
        Q_B(1,ikmax+1) = 0.0D0
        Q_B(2,ikmax+1) = 0.0D0
        Q_B(3,ikmax+1) = 0.0D0
        Q_B(4,ikmax+1) = 0.0D0
c
      else
c
        write(71,*) 'INCORRECT FLUXES'
        STOP
c
      endif
c
      return
      end
c
c
      subroutine residual
      implicit none
c
c
c           Hj-1          Hj         Hj+1         Hj+2
c           Fj-1          Fj         Fj+1         Fj+2
c                 Rj-1         Rj         Rj+1
c
c            |      x      |     x     |     x      |
c                 j-1           j          j+1
c
      include 'cfd_osm_com'
      integer in,ik
      real*8 dx_ab, d1dx1, rho_a, TJ_a, RG_impl,
     >      T_c, v_c, R1_c, R2_c, rg_expl,
     >      interpol2p, dum1
      external d1dx1
c
      do ik = 1, ikmax
        dx_ab = xf(ik+1) - xf(ik)
        do in = 1, 4
          RF(in,ik) =  d1dx1(H(in,ik), H(in,ik+1), dx_ab)
     >               - d1dx1(F(in,ik), F(in,ik+1), dx_ab)
     >               +  Q(in,ik) + Q_B(in,ik)
          if (in.eq.3) then
            RF(3,ik) = RF(3,ik) + Q_ei(ik) + Q_joule(ik)
          elseif (in.eq.4) then
            RF(4,ik) = RF(4,ik) - Q_ei(ik) - Q_joule(ik)
          endif
          RH(in,ik) = d1dx1(H(in,ik), H(in,ik+1), dx_ab)
          RG(in,ik) = - d1dx1(G(in,ik), G(in,ik+1), dx_ab)
          R(in,ik) = RF(in,ik) + RG(in,ik)
          if (abs(R(in,ik)).lt.1.0D99.and.
     >         (R(in,ik).le.1.0D0.or.R(in,ik).ge.1.0D0)) then
              dum1 = 1.0D0
          else
              write(71,*) 'R_NANQ', ik, RF(in,ik), RG(in,ik),
     >        W(1,ik), W(2,ik), W(3,ik)
              write(71,*) U(1,ik), U(2,ik), U(3,ik), U(4,ik)
              write(71,*) F(1,ik), F(2,ik), F(3,ik), F(4,ik)
              write(71,*) G(1,ik), G(2,ik), G(3,ik), G(4,ik)
              write(71,*) H(1,ik), H(2,ik), H(3,ik), H(4,ik)
              write(71,*) Q(1,ik), Q(2,ik), Q(3,ik), Q(4,ik)
              write(71,*)  ka(ik), mu(ik), G(in,ik), G(in,ik+1),
     >        dx_ab
              call debugring
              STOP
c
c             R(in,ik) = 0.0D0
c
          endif
        enddo
      enddo
c
c      write(71,'(A20,6G14.8)') 'U,RF,RG,RH,R0:',
c     >   U(4,0), RF(4,0), RG(4,0),
c     >   RH(4,0), R(4,0), dtc(0)
c      write(71,'(A20,6G14.8)') 'U,RF,RG,RH,R1:',
c     >   U(4,1), RF(4,1), RG(4,1),
c     >   RH(4,1), R(4,1), dtc(1)
c      write(71,'(A20,6G14.8)') 'U,RF,RG,RH,R1:',
c     >   U(4,2), RF(4,2), RG(4,2),
c     >   RH(4,2), R(4,2), dtc(2)
c
c      do ik = 1, ikmid
c          R(1,ikmax+1-ik) = R(1,ik)
c          R(2,ikmax+1-ik) = - R(2,ik)
c          R(3,ikmax+1-ik) = R(3,ik)
c      enddo
c
c      do in = 1, 3
c        R(in,ikmid) = (R(in,ikmid-1) + R(in,ikmid+1)) / 2.0D0
c      enddo
c
      return
      end
c
c
      subroutine form_tri_mtrx
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in,im
      real*8 dx_j, dx_ab, dx_bc, zhi, ep_2a, ep_2b, ep_2c,
     >       dummy, al_2a, al_2b, al_2c
c
      if (a_flag.eq.1.and.a_2T.eq.0) then
        do in = 1, 3
          do ik = 1, ikmax
            dx_ab = xc(ik) - xc(ik-1)
            dx_bc = xc(ik+1) - xc(ik)
            dx_j = (dx_ab + dx_bc) / 2.0D0
            zhi = dtc(ik) / (dx_j * 2.0D0)
            ep_2a =  - ep2(ik-1) * dx_j / dx_ab
            ep_2c =  - ep2(ik+1) * dx_j / dx_bc
            ep_2b =    ep2(ik) *
     >                   dx_j * (1.0D0/dx_ab + 1.0D0/dx_bc)
            if (ik.eq.1) then
              ep_2b =   ep2(ik) * dx_j / dx_bc
            elseif (ik.eq.ikmax) then
              ep_2b =  ep2(ik) * dx_j / dx_ab
            endif
            bbbb(in,ik) = 1.0D0 + zhi * ep_2b
            aaaa(in,ik) = zhi * ( - lam_mat(in,ik-1) + ep_2a)
            cccc(in,ik) = zhi * (   lam_mat(in,ik+1) + ep_2c)
            dummy =  (dabs(aaaa(in,ik))+dabs(cccc(in,ik))) /
     >                 bbbb(in,ik)
          enddo
        enddo
      elseif (a_flag.eq.0.and.a_2T.eq.0) then
        do in = 1, 3
          do im = 1, 3
            do ik = 1, ikmax
             dx_ab = xc(ik) - xc(ik-1)
             dx_bc = xc(ik+1) - xc(ik)
             dx_j = (dx_ab + dx_bc) / 2.0D0
             zhi = dtc(ik) / dx_j
c
             ep_2a =  - ep2(ik-1) * dx_j / dx_ab
             ep_2c =  - ep2(ik+1) * dx_j / dx_bc
             ep_2b =    ep2(ik) *
     >                   dx_j * (1.0D0/dx_ab + 1.0D0/dx_bc)
c
             if (ik.eq.1) then
               ep_2b =   ep2(ik) * dx_j / dx_bc
             elseif (ik.eq.ikmax) then
               ep_2b =   ep2(ik) * dx_j / dx_ab
             endif
c
             aaaa3(in,im,ik) = 0.5D0 * zhi *
     >         ( -jac_A3(in,im,ik-1) + ep_2a * unit(in,im) )
             bbbb3(in,im,ik) = (1.0D0 + 0.5D0 * zhi * ep_2b)
     >                         * unit(in,im)
             cccc3(in,im,ik) = 0.5D0 * zhi *
     >         (  jac_A3(in,im,ik+1) + ep_2c * unit(in,im) )
             ffff3(in,ik)    =   R(in,ik) * dtc(ik)
c
             if (in.eq.1) then
               al_2a = 0.0D0
               al_2c = 0.0D0
             elseif (in.eq.2) then
               al_2a =  - 0.5D0 * (mu(ik-1) + mu(ik)) / dx_ab
               al_2c =  - 0.5D0 * (mu(ik+1) + mu(ik)) / dx_bc
             elseif (in.eq.3) then
               al_2a =  - 0.5D0 * (ka(ik-1) + ka(ik)) / dx_ab
               al_2c =  - 0.5D0 * (ka(ik+1) + ka(ik)) / dx_bc
             endif
             al_2b =  - (al_2a + al_2c)
c
             aaaa3(in,im,ik) = aaaa3(in,im,ik) +
     >                         zhi * al_2a * jac_B3(in,im,ik-1)
             bbbb3(in,im,ik) = bbbb3(in,im,ik) +
     >                         zhi * al_2b * jac_B3(in,im,ik)
             cccc3(in,im,ik) = cccc3(in,im,ik) +
     >                         zhi * al_2c * jac_B3(in,im,ik+1)
c
             if (in.eq.1) then
               al_2a = 0.0D0
               al_2c = 0.0D0
             elseif (in.eq.2) then
               al_2a = - 0.5D0 * ( mu(ik-1) / (W(3,ik-1)
     >              / (W(1,ik-1) * a_Temp)) * jac_B3(3,im,ik-1)
     >              +  mu(ik) / (W(3,ik)
     >              / (W(1,ik) * a_Temp)) * jac_B3(3,im,ik) )
     >                   / dx_ab
               al_2c = - 0.5D0 * ( mu(ik+1) / (W(3,ik+1)
     >              / (W(1,ik+1) * a_Temp)) * jac_B3(3,im,ik+1)
     >              +  mu(ik) / (W(3,ik)
     >              / (W(1,ik) * a_Temp)) * jac_B3(3,im,ik) )
     >                   / dx_bc
             elseif (in.eq.3) then
               al_2a = - 0.5D0 * ( ka(ik-1) / (W(3,ik-1)
     >              / (W(1,ik-1) * a_Temp)) * jac_B3(3,im,ik-1)
     >              +  ka(ik) / (W(3,ik)
     >              / (W(1,ik) * a_Temp)) * jac_B3(3,im,ik) )
     >                   / dx_ab
               al_2c = - 0.5D0 * ( ka(ik+1) / (W(3,ik+1)
     >              / (W(1,ik+1) * a_Temp)) * jac_B3(3,im,ik+1)
     >              +  ka(ik) / (W(3,ik)
     >              / (W(1,ik) * a_Temp)) * jac_B3(3,im,ik) )
     >                   / dx_bc
             endif
             al_2b =  - (al_2a + al_2c)
c
c             if (in.eq.2) then
c               aaaa3(in,im,ik) = aaaa3(in,im,ik) +
c     >                         zhi * al_2a * W(2,ik-1)
c               bbbb3(in,im,ik) = bbbb3(in,im,ik) +
c     >                         zhi * al_2b * W(2,ik)
c               cccc3(in,im,ik) = cccc3(in,im,ik) +
c     >                         zhi * al_2c * W(2,ik+1)
c             elseif (in.eq.3) then
c               aaaa3(in,im,ik) = aaaa3(in,im,ik) +
c     >                         zhi * al_2a * (W(3,ik-1)
c     >                         / (W(1,ik-1) * a_Temp))
c               bbbb3(in,im,ik) = bbbb3(in,im,ik) +
c     >                         zhi * al_2b * (W(3,ik)
c     >                         / (W(1,ik) * a_Temp))
c               cccc3(in,im,ik) = cccc3(in,im,ik) +
c     >                         zhi * al_2c * (W(3,ik+1)
c     >                         / (W(1,ik+1) * a_Temp))
c             endif
c
            enddo
          enddo
        enddo
c
      elseif (a_flag.eq.0.and.a_2T.eq.1) then
c
        do in = 1, 4
          do im = 1, 4
            do ik = 1, ikmax
             dx_ab = xc(ik) - xc(ik-1)
             dx_bc = xc(ik+1) - xc(ik)
             dx_j = (dx_ab + dx_bc) / 2.0D0
             zhi = dtc(ik) / dx_j
c
             ep_2a =  - ep2(ik-1) * dx_j / dx_ab
             ep_2c =  - ep2(ik+1) * dx_j / dx_bc
             ep_2b =    ep2(ik) *
     >                   dx_j * (1.0D0/dx_ab + 1.0D0/dx_bc)
c
             if (ik.eq.1) then
               ep_2b =   ep2(ik) * dx_j / dx_bc
             elseif (ik.eq.ikmax) then
               ep_2b =   ep2(ik) * dx_j / dx_ab
             endif
c
             aaaa4(in,im,ik) = 0.5D0 * zhi *
     >         ( -jac_A4(in,im,ik-1) + ep_2a * unit(in,im) )
             bbbb4(in,im,ik) = (1.0D0 + 0.5D0 * zhi * ep_2b)
     >                          * unit(in,im)
             cccc4(in,im,ik) = 0.5D0 * zhi *
     >         (  jac_A4(in,im,ik+1) + ep_2c * unit(in,im) )
             ffff4(in,ik)    =   R(in,ik) * dtc(ik)
c
             if (in.eq.1) then
               al_2a = 0.0D0
               al_2c = 0.0D0
             elseif (in.eq.2) then
               al_2a =  - 0.5D0 * (mu(ik-1) + mu(ik)) / dx_ab
               al_2c =  - 0.5D0 * (mu(ik+1) + mu(ik)) / dx_bc
             elseif (in.eq.3) then
               al_2a =  - 0.5D0 * (ka_i(ik-1) + ka_i(ik)) / dx_ab
               al_2c =  - 0.5D0 * (ka_i(ik+1) + ka_i(ik)) / dx_bc
             elseif (in.eq.4) then
               al_2a =  - 0.5D0 * (ka_e(ik-1) + ka_e(ik)) / dx_ab
               al_2c =  - 0.5D0 * (ka_e(ik+1) + ka_e(ik)) / dx_bc
             endif
             al_2b =  - (al_2a + al_2c)
c
             aaaa4(in,im,ik) = aaaa4(in,im,ik) +
     >                         zhi * al_2a * jac_B4(in,im,ik-1)
             bbbb4(in,im,ik) = bbbb4(in,im,ik) +
     >                         zhi * al_2b * jac_B4(in,im,ik)
             cccc4(in,im,ik) = cccc4(in,im,ik) +
     >                         zhi * al_2c * jac_B4(in,im,ik+1)
c
             bbbb4(in,im,ik) = bbbb4(in,im,ik) - dtc(ik) *
     >                 jac_Q4(in,im,ik)
c
c             bbbb4(in,im,ik) = bbbb4(in,im,ik) - dtc(ik) *
c     >                jac_A4(in,im,ik) + jac_A4(in,im,ik+1) +
c     >                 jac_B4(in,im,ik) + jac_B4(in,im,ik+1))
c     >                / (2.0D0 * L_B(ik))
c
c             if (in.eq.1) then
c               al_2a = 0.0D0
c               al_2c = 0.0D0
c             elseif (in.eq.2) then
c               al_2a = - 0.5D0 * ( mu(ik-1) / (W(3,ik-1)
c     >              / (W(1,ik-1) * a_Temp)) * jac_B4(3,im,ik-1)
c     >              +  mu(ik) / (W(3,ik)
c     >              / (W(1,ik) * a_Temp)) * jac_B4(3,im,ik) )
c     >                   / dx_ab
c               al_2c = - 0.5D0 * ( mu(ik+1) / (W(3,ik+1)
c     >              / (W(1,ik+1) * a_Temp)) * jac_B4(3,im,ik+1)
c     >              +  mu(ik) / (W(3,ik)
c     >              / (W(1,ik) * a_Temp)) * jac_B4(3,im,ik) )
c     >                   / dx_bc
c             elseif (in.eq.3) then
c               al_2a = - 0.5D0 * ( ka_i(ik-1) / (W(3,ik-1)
c     >              / W(1,ik-1)) * jac_B4(3,im,ik-1)
c     >              +  ka_i(ik) / (W(3,ik)
c     >              / W(1,ik)) * jac_B4(3,im,ik) )
c     >                   / dx_ab
c               al_2c = - 0.5D0 * ( ka_i(ik+1) / (W(3,ik+1)
c     >              / W(1,ik+1)) * jac_B4(3,im,ik+1)
c     >              +  ka_i(ik) / (W(3,ik)
c     >              / W(1,ik)) * jac_B4(3,im,ik) )
c     >                   / dx_bc
c             elseif (in.eq.4) then
c               al_2a = - 0.5D0 * ( ka_e(ik-1) / (W(4,ik-1)
c     >              / W(1,ik-1)) * jac_B4(4,im,ik-1)
c     >              +  ka_e(ik) / (W(4,ik)
c     >              / W(1,ik)) * jac_B4(4,im,ik) )
c     >                   / dx_ab
c               al_2c = - 0.5D0 * ( ka_e(ik+1) / (W(4,ik+1)
c     >              / W(1,ik+1)) * jac_B4(4,im,ik+1)
c     >              +  ka_e(ik) / (W(4,ik)
c     >              / W(1,ik)) * jac_B4(4,im,ik) )
c     >                   / dx_bc
c             endif
c             al_2b =  - (al_2a + al_2c)
c
            enddo
          enddo
        enddo
c
      else
c
        write(71,*) 'INCORRECT TRIMATRIX'
        STOP
c
      endif
c
      return
      end
c
c
      subroutine max_res_R
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in,jj
      real*8 R_max, res_sum, err_sum,aa,bb,dx
c
      max_res = 0.0D0
      max_err = 0.0D0
      res_sum = 0.0D0
      err_sum = 0.0D0
c
      if (a_2T.eq.0) then
        jj = 3
      else
        jj = 4
      endif
c
      do ik = 1, ikmax
        do in = 1, jj
           dx = (xc(ik+1) - xc(ik-1)) / 2.0D0
           aa = abs(R(in,ik) * dx / F(in,0))
           bb = abs(R(in,ik) * dtc(ik)/U(in,ik))
           res_sum = res_sum + aa ** 2.0D0
           err_sum = err_sum + bb ** 2.0D0
           if (aa.gt.max_res) then
              max_res = aa
           endif
           if (bb.gt.max_err) then
              max_err = bb
           endif
        enddo
      enddo
c
      rms_res = dsqrt(res_sum) / (3.0D0 * ikmax)
      rms_err = dsqrt(err_sum) / (3.0D0 * ikmax)
c
c      write(71,*) 'res:', max_res, rms_res
c
      return
      end
c
c
      subroutine relax_RE(relax_var)
      implicit none
      real*8 relax_var
c
      include 'cfd_osm_com'
      integer in,ik
c
      do in = 1, 4
        do ik = 1, ikmax
          RE(in,ik) = RE(in,ik) * (1.0D0 - relax_var)
     >             + relax_var * RH(in,ik)
        enddo
      enddo
c
      return
      end
c
c
      subroutine smoothU(eps_s)
      implicit none
      real*8 eps_s
c
      include 'cfd_osm_com'
      real*8 U_new_s, A_s, dx_ab, dx_bc, dx_b
      integer in,ik
c
      call U0fromU
c
      do in = 1, 4
        do ik = 1, ikmax
          dx_bc = xc(ik+1) - xc(ik)
          dx_ab = xc(ik) - xc(ik-1)
          dx_b  = (xc(ik+1) - xc(ik-1)) / 2.0D0
          U_new_s =    U0(in,ik-1) * dx_ab
     >              -  U0(in,ik)   * 2.0D0 * dx_b
     >              +  U0(in,ik+1) * dx_bc
          U(in,ik) = U0(in,ik) + eps_s * U_new_s
        enddo
      enddo
c
      return
      end
c
c
      subroutine smoothQ(eps_s)
      implicit none
      real*8 eps_s
c
      include 'cfd_osm_com'
      real*8 Q_new_s, A_s, dx_ab, dx_bc, dx_b
      integer in,ik
c
      call Q0fromQ
c
      do in = 1, 4
        do ik = 1, ikmax
          dx_bc = xc(ik+1) - xc(ik)
          dx_ab = xc(ik) - xc(ik-1)
          dx_b  = (xc(ik+1) - xc(ik-1)) / 2.0D0
          Q_new_s =    Q0(in,ik-1) * dx_ab
     >              -  Q0(in,ik)   * 2.0D0 * dx_b
     >              +  Q0(in,ik+1) * dx_bc
          Q(in,ik) = Q0(in,ik) + eps_s * Q_new_s
        enddo
      enddo
c
      return
      end
c
c
      subroutine initgrid
      implicit none
c
      include 'params'
      include 'cgeom'
      include 'sol23_com'
      include 'cfd_osm_com'
c
      integer ik,in,im, ik_xpt
      real*8 grid_exp, beta, c_s, a_3, T_eV, f_grid,
     >      dx_T, dx_G, dx_c, dx_max, P_cond, taue, chi_e,
     >      dx_0, n_0, T_0, M_0,
     >      dx_N, n_N, T_N, M_N
      external f_grid
c
      if (a_2T.eq.0) then
        T_0 = T0_grid
        M_0 = M0_BC
        v_0 = M_0 * dsqrt(a_Temp * T_0 / mass)
        n_0 = MASS_0 / abs(v_0 * mass)
c
        T_N = TN_grid
        M_N = MN_BC
        v_N = M_N * dsqrt(a_Temp * T_N / mass)
        n_N = MASS_N / abs(v_N * mass)
c
        T_eV = T_0 / eV
        taue = 3.44D11 * (T_eV ** 1.5D0) / (n_0 * 15.0D0)
        chi_e = 3.16D0 * n_0 * T_0 * taue / mass_e
        dx_T   = a_grid0 * chi_e /
     >       ((gamma_sum - 2.5D0) * abs(n_0 * v_0))
        dx_G   = abs(a_grid0 * 0.03D0 * L_iz_0)
        dx_0  = min(dx_T , dx_G)
        dx0_BC = dx_0
c
        T_eV = T_N / eV
        taue = 3.44D11 * (T_eV ** 1.5D0) / (n_N * 15.0D0)
        chi_e = 3.16D0 * n_N * T_N * taue / mass_e
        dx_T   = a_grid0 * chi_e /
     >       ((gamma_sum - 2.5D0) * abs(n_N * v_N))
        dx_G   = abs(a_grid0 * 0.03D0 * L_iz_N)
        dx_N  = min(dx_T , dx_G)
        dxN_BC = dx_N
c
      elseif (a_2T.eq.1) then
c
        M_0 = M0_BC
        v_0 = M_0 * dsqrt((Ti0_BC + Te0_BC) / mass)
        n_0 = MASS_0 / abs(v_0 * mass)
c
        M_N = MN_BC
        v_N = M_N * dsqrt((TiN_BC + TeN_BC) / mass)
        n_N = MASS_N / abs(v_N * mass)
c
        T_eV = Te0_BC / eV
        taue = 3.44D11 * (T_eV ** 1.5D0) / (n_0 * 15.0D0)
        chi_e = 3.16D0 * n_0 * Te0_BC * taue / mass_e
        dx_T  = a_grid0 * chi_e /
     >       (gamma_e * abs(n_0 * v_0))
        dx_G  = abs(a_grid0 * s23_par_gridg * L_iz_0)
        dx_0  = min(dx_T , dx_G)
        dx0_BC = dx_0
c
        T_eV = TeN_BC / eV
        taue = 3.44D11 * (T_eV ** 1.5D0) / (n_N * 15.0D0)
        chi_e = 3.16D0 * n_N * TeN_BC * taue / mass_e
        dx_T   = a_grid0 * chi_e /
     >       (gamma_e * abs(n_N * v_N))
        dx_G   = abs(a_grid0 * s23_par_gridg * L_iz_N)
        dx_N  = min(dx_T , dx_G)
        dxN_BC = dx_N
c
        if (a_PIN.ge.1) then
          dx0_BC = s23_par_grid_dx0
          dxN_BC = s23_par_grid_dx0
        endif
c
c        write(71,'(A16,2G16.8)') 'dx0_G,T', dx_G, dx_T
c
      endif
c
c      write(71,'(A16,2G16.8)') 'dx0N_BC', dx0_BC, dxN_BC
c
      if (a_PIN.ge.1) then
c
        call readgrid
c
      else
c
        xc(1) = dx0_BC
        do ik = 2, ikmax_0
           if (xc(ik-1).ge.L_izoffset*0.5D0.and.
     >         xc(ik-1).le.L_izoffset+L_iz*2.0D0) then
              xc(ik) = xc(ik-1) + (xc(ik-1) - xc(ik-2)) *
     >                 1.0D0
           else
              xc(ik) = xc(ik-1) + (xc(ik-1) - xc(ik-2)) *
     >                 a_gridexp
           endif
           if (xc(ik).ge.x_max/2.0D0) then
              ikmid = ik
              goto 300
           endif
        enddo
c
 300    ikmid = ik
        ikmax = 2 * ikmid - 1
c
c        write(0,*) 'ikmid,max', ikmid, ikmax
c
        do ik = 0, ikmid
           xc(ikmax+1-ik) = 2.0D0 * xc(ikmid) - xc(ik)
        enddo
c
      endif
c
c      xc(0) = 0.0D0
c      xc(1) = 1.0D0
c      do ik = 2, ikmax+1
c         if(ik.le.ikmid) then
c           xc(ik) = xc(ik-1) + (xc(ik-1) - xc(ik-2)) *
c     >              grid_exp
c         else
c           xc(ik) = 2.0D0 * xc(ikmid) - xc(ikmax+1-ik)
c         endif
c      enddo
c
      do ik = 1, ikmax+1
         xc(ik) = x_max * (xc(ik) / xc(ikmax+1))
      enddo
c
      xf(0) = xc(0)
      xf(ikmax+2) = xc(ikmax+1)
      do ik = 1, ikmax+1
         xf(ik) = (xc(ik) + xc(ik-1)) / 2.0D0
      enddo
c
c      write(71,*) '   ik        xc[m]          xf[m]'
c      do ik = 1, ikmax+1
c         write(71,'(I5,2G16.8)') ik , xc(ik), xf(ik)
c      enddo
c      write(71,*)
c
      ik_xpt = 1
      do ik  = 1, nks(ir_cfd) / 2
        if (zs(ik,ir_cfd).ge.zxp) then
          ik_xpt = ik
        endif
c        write(71,*) ik, ik_xpt,
c     >              zs(ik,ir_cfd), zxp, kss(ik,ir_cfd)
      enddo
c
      if (ik_xpt.eq.nks(ir_cfd)/2) then
        x_xpt_0 = 0.0D0
      else
        x_xpt_0 = kss(ik_xpt + 1,ir_cfd)
      endif
c
      ik_xpt = nks(ir_cfd)
      do ik  = nks(ir_cfd), nks(ir_cfd) / 2, -1
        if (zs(ik,ir_cfd).ge.zxp) then
          ik_xpt = ik
        endif
c        write(71,*) ik, ik_xpt,
c     >              zs(ik,ir_cfd), zxp, kss(ik,ir_cfd)
      enddo
c
      if (ik_xpt.eq.nks(ir_cfd)/2) then
        x_xpt_N = x_max
      else
        x_xpt_N = kss(ik_xpt - 1,ir_cfd)
      endif
c
c      write(71,*) 'x_xpt0,N', ir_cfd, x_xpt_0, x_xpt_N
c
      return
      end
c
c
      subroutine check_symmetry
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
      real*8 qdel1, qdel2, qdel3, udel1, udel2, udel3,
     >  qdel, udel, x_symm
c
      qdel_max = 0.0D0
      udel_max = 0.0D0
c
      do ik = 0, ikmid-1
          qdel1 = Q(1,ikmax+1-ik) - Q(1,ik)
          qdel2 = Q(2,ikmax+1-ik) + Q(2,ik)
          qdel3 = Q(3,ikmax+1-ik) - Q(3,ik)
          udel1 = U(1,ikmax+1-ik) - U(1,ik)
          udel2 = U(2,ikmax+1-ik) + U(2,ik)
          udel3 = U(3,ikmax+1-ik) - U(3,ik)
          if (Q(1,ik).ne.0.0D0) qdel1 = qdel1 / Q(1,ik)
          if (Q(2,ik).ne.0.0D0) qdel2 = qdel2 / Q(2,ik)
          if (Q(3,ik).ne.0.0D0) qdel3 = qdel3 / Q(3,ik)
          if (U(1,ik).ne.0.0D0) udel1 = udel1 / U(1,ik)
          if (U(2,ik).ne.0.0D0) udel2 = udel2 / U(2,ik)
          if (U(3,ik).ne.0.0D0) udel3 = udel3 / U(3,ik)
          qdel = max(abs(qdel1),abs(qdel2), abs(qdel3))
          udel = max(abs(udel1),abs(udel2), abs(udel3))
          if (qdel.ge.Qdel_max) Qdel_max = qdel
          if (udel.ge.Udel_max)  then
             Udel_max = udel
             x_symm = xc(ik)
          endif
      enddo
c      write(0,'(A30,2G16.8)') ' MAX. ASYMMETRY (Q,U) ',
c     >        qdel_max, udel_max
      write(71,'(A30,2G16.8)') ' MAX. ASYMMETRY (Q,U) ',
     >        qdel_max, udel_max
c
      return
      end
c
c
      subroutine print_symmetry
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
      real*8 qdel1, qdel2, qdel3, udel1, udel2, udel3,
     >      rdel1, rdel2, rdel3
c
      do ik = 0, ikmid-1
          qdel1 = Q(1,ikmax+1-ik) - Q(1,ik)
          qdel2 = Q(2,ikmax+1-ik) + Q(2,ik)
          qdel3 = Q(3,ikmax+1-ik) - Q(3,ik)
          udel1 = U(1,ikmax+1-ik) - U(1,ik)
          udel2 = U(2,ikmax+1-ik) + U(2,ik)
          udel3 = U(3,ikmax+1-ik) - U(3,ik)
          rdel1 = R(1,ikmax+1-ik) - R(1,ik)
          rdel2 = R(2,ikmax+1-ik) + R(2,ik)
          rdel3 = R(3,ikmax+1-ik) - R(3,ik)
          if (Q(1,ik).ne.0.0D0) qdel1 = qdel1 / Q(1,ik)
          if (Q(2,ik).ne.0.0D0) qdel2 = qdel2 / Q(2,ik)
          if (Q(3,ik).ne.0.0D0) qdel3 = qdel3 / Q(3,ik)
          if (U(1,ik).ne.0.0D0) udel1 = udel1 / U(1,ik)
          if (U(2,ik).ne.0.0D0) udel2 = udel2 / U(2,ik)
          if (U(3,ik).ne.0.0D0) udel3 = udel3 / U(3,ik)
          if (R(1,ik).ne.0.0D0) rdel1 = rdel1 / R(1,ik)
          if (R(2,ik).ne.0.0D0) rdel2 = rdel2 / R(2,ik)
          if (R(3,ik).ne.0.0D0) rdel3 = rdel3 / R(3,ik)
          write(71,'(I5,9G12.6)')  ik,
     >     qdel1, qdel2, qdel3, udel1, udel2, udel3,
     >     rdel1, rdel2, rdel3
      enddo
      write(71,*)
c
      return
      end
c
c
      subroutine flatten_U
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
c
      call setUBC
      U(1,ikmax+1) = U(1,NNN)
      U(2,ikmax+1) = U(2,NNN)
      U(3,ikmax+1) = U(3,NNN)
      U(4,ikmax+1) = U(4,NNN)
c
      do ik = 0, ikmax + 1
           U(1,ik) = U(1,0) +
     >       (U(1,ikmax+1) - U(1,0)) * xc(ik) / xc(ikmax+1)
           U(2,ik) = U(2,0) +
     >       (U(2,ikmax+1) - U(2,0)) * xc(ik) / xc(ikmax+1)
           U(3,ik) = U(3,0) +
     >       (U(3,ikmax+1) - U(3,0)) * xc(ik) / xc(ikmax+1)
           U(4,ik) = U(4,0) +
     >       (U(4,ikmax+1) - U(4,0)) * xc(ik) / xc(ikmax+1)
      enddo
c
      call WfromU
c
      return
      end
c
c
      subroutine print_Q
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
      real*8 sum_1, sum_2, sum_3, sum_4, sum_5, sum_6,
     >       sum_7, sum_8, sum_9, sum_10, temp_1, temp_3
c
      sum_1 = 0.0D0
      sum_2 = 0.0D0
      sum_3 = 0.0D0
      sum_4 = 0.0D0
      sum_5 = 0.0D0
      sum_6 = 0.0D0
      sum_7 = 0.0D0
      sum_8 = 0.0D0
      sum_9 = 0.0D0
      sum_10= 0.0D0
c
      write(71,'(A16,5A16)') 'x [m]   ','Q1 [SI]    '
     >    ,'Q2 [SI]    ',
     >    'Q3 [SI]    ', 'Q4 [SI]    '
      do ik = 0, ikmax + 1
         write(71,'(5G16.8)')  xc(ik), Q(1,ik) + Q_B(1,ik),
     >   Q(2,ik) + Q_B(2,ik),
     >   Q(3,ik) + Q_ei(ik) + Q_joule(ik) + Q_B(3,ik),
     >   Q(4,ik) - Q_ei(ik) - Q_joule(ik) + Q_B(4,ik)
      enddo
      write(71,*)
c
      write(71,'(A16,5A16)') 'x [m]   ','Q3_in [SI]    '
     >    ,'Q4_in [SI]    ',
     >    'Q_ei [SI]    ', 'Q_j [SI]    '
      do ik = 0, ikmax + 1
         write(71,'(5G16.8)')  xc(ik), Q(3,ik),
     >     Q(4,ik), Q_ei(ik), Q_joule(ik)
      enddo
      write(71,*)
c
      write(71,'(A16,5A16)') 'x [m]   ','Q1_B   '
     >    ,'Q2_B   ',
     >    'Q3_B    ', 'Q4_B  '
      do ik = 0, ikmax + 1
         write(71,'(5G16.8)')  xc(ik), Q_B(1,ik),
     >     Q_B(2,ik), Q_B(3,ik), Q_B(4,ik)
      enddo
      write(71,*)
c
c
      do ik = 1, ikmax+1
         sum_1 = sum_1 + Q_ei(ik-1) * (xf(ik) - xf(ik-1))
         sum_2 = sum_2 + Q_joule(ik-1) * (xf(ik) - xf(ik-1))
         sum_3 = sum_3 + Q_B(3,ik-1) * (xf(ik) - xf(ik-1))
         sum_4 = sum_4 + Q_B(4,ik-1) * (xf(ik) - xf(ik-1))
         sum_5 = sum_5 + Q_perp(3,ik-1) * (xf(ik) - xf(ik-1))
         sum_6 = sum_6 + Q(3,ik-1) * (xf(ik) - xf(ik-1))
         sum_7 = sum_7 + (Q(3,ik) + Q_ei(ik) + Q_joule(ik)
     >                    + Q_B(3,ik)) * (xf(ik) - xf(ik-1))
         sum_8 = sum_8 + Q_perp(4,ik-1) * (xf(ik) - xf(ik-1))
         sum_9 = sum_9 + Q(4,ik-1) * (xf(ik) - xf(ik-1))
         sum_10= sum_10 + (Q(4,ik) - Q_ei(ik) - Q_joule(ik)
     >                    + Q_B(4,ik)) * (xf(ik) - xf(ik-1))
      enddo
c
      write(71,*) 'P_ei      [MW] ', sum_1/ 1.0E6
      write(71,*) 'P_joule   [MW] ', sum_2/ 1.0E6
      write(71,*) 'P_ion_B   [MW] ', sum_3/ 1.0E6
      write(71,*) 'P_elec_B  [MW] ', sum_4/ 1.0E6
      write(71,*)
      write(71,*) 'P_i_in    [MW] ', sum_5 / 1.0E6
      write(71,*) 'P_i_in+rad[MW] ', sum_6 / 1.0E6
      write(71,*) 'P_i_net   [MW] ', sum_7 / 1.0E6
      write(71,*) 'P_e_in    [MW] ', sum_8 / 1.0E6
      write(71,*) 'P_e_in+rad[MW] ', sum_9 / 1.0E6
      write(71,*) 'P_e_net   [MW] ', sum_10/ 1.0E6
      write(71,*)
c
      write(71,'(A16,4A16)') 'x [m]   ','Qp1 [SI]    '
     >    ,'Qp2 [SI]    ',
     >    'Qp3 [SI]    ', 'Qp4 [SI]   '
      do ik = 0, ikmax + 1
         write(71,'(5G16.8)')  xc(ik), Q_perp(1,ik),
     >   Q_perp(2,ik), Q_perp(3,ik), Q_perp(4,ik)
      enddo
      write(71,*)
c
      temp_1 = 1.0D23 * mass
      temp_3 = 1.0D6
c
      do ik = 0, ikmax + 1
         write(73,'(5G16.8)')  xc(ik)/x_max,
     >   (Q(1,ik) - Q_perp(1,ik))/temp_1, Q_perp(1,ik)/temp_1,
     >   Q_B(1,ik)/temp_1 , (Q(1,ik) + Q_B(1,ik))/temp_1
      enddo
      write(73,*)
c
      do ik = 0, ikmax + 1
         write(73,'(5G16.8)')  xc(ik)/x_max,
     >   Q(2,ik) - Q_perp(2,ik), Q_perp(2,ik), Q_B(2,ik),
     >   Q(2,ik) + Q_B(2,ik)
      enddo
      write(73,*)
c
      do ik = 0, ikmax + 1
         write(73,'(6G16.8)')  xc(ik)/x_max,
     >   (Q(3,ik) - Q_perp(3,ik))/temp_3, Q_perp(3,ik)/temp_3,
     >   Q_B(3,ik)/temp_3, (Q_ei(ik) + Q_joule(ik))/temp_3,
     >   (Q(3,ik) + Q_B(3,ik) + Q_ei(ik) + Q_joule(ik))/temp_3
      enddo
      write(73,*)
c
      do ik = 0, ikmax + 1
         write(73,'(6G16.8)')  xc(ik)/x_max,
     >   (Q(4,ik) - Q_perp(4,ik))/temp_3, Q_perp(4,ik)/temp_3,
     >    Q_B(4,ik)/temp_3, (- Q_ei(ik) - Q_joule(ik))/temp_3,
     >    (Q(4,ik) + Q_B(4,ik) - Q_ei(ik) - Q_joule(ik))/temp_3
      enddo
      write(73,*)
c
      return
      end
c
c
      subroutine print_U
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
      real*8 dum1, dum2, dum3
c
      write(71,'(A16,4A16)') 'x [m]   ','mn [SI]    '
     >    ,'mnv [SI]    ',
     >    'Ti [eV]    ','Te [eV]    '
      do ik = 0, ikmax + 1
         write(71,'(5G16.8)')  xc(ik), U(1,ik),
     >   U(2,ik), W(3,ik) / (a_Temp * W(1,ik) * eV),
     >   W(4,ik) / (a_Temp * W(1,ik) * eV)
         write(73,'(7G16.8)')  xc(ik)/x_max, W(1,ik)/1.0D20,
     >   U(2,ik)/(mass*1.0D23), W(3,ik) / (W(1,ik) * eV),
     >   W(4,ik) / (W(1,ik) * eV), W(3,ik) + W(4,ik) +
     >   U(2,ik) ** 2.0 / U(1,ik),
     >   W(2,ik) / dsqrt((W(3,ik) + W(4,ik)) / U(1,ik))
      enddo
      write(71,*)
      write(73,*)
c
      return
      end
c
c
      subroutine print_W
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
      real*8  a_2, a_3, beta, Re_inv, c_s,
     >  Pe_inv, Pr, T_J, dx_min, Mach1, Mach2, interpol2p
c
      write(71,'(A16,5A16)') 'x [m]   ','n [m^-3]    '
     >    ,'pi [Pa]    ', 'pe [Pa]  ', 'Mach_s [-]  '
     >    ,' v [m/s]   '
c
c      write(71,'(A5,4A16)') 'ik','1/Re [-]  '
c     >    ,'1/Pe [-]  ',
c     >    'Mach_s [-]    ','Mach_a [-]    '
c
      do ik = 0, ikmax+1
         dx_min = xf(ik+1) - xf(ik)
         if (a_2T.eq.0) then
           c_s = dsqrt( W(3,ik) / U(1,ik))
         else
           c_s = dsqrt((W(3,ik) + W(4,ik)) / U(1,ik))
         endif
         Mach1 = W(2,ik) / c_s
         Mach2 = W(2,ik) / c_s_j(ik)
         a_3 = ka(ik) / (W(1,ik)*c_v_0)
         a_2 = mu(ik) / U(1,ik)
         beta = dx_min * (c_s_j(ik) + dabs(W(2,ik)))
         Re_inv = a_2 / beta
         Pe_inv = a_3 / beta
         Pr = a_2 / a_3
         write(71,'(6G16.8)') xc(ik), W(1,ik), W(3,ik),
     >     W(4,ik), Mach1, W(2,ik)
c
      enddo
      write(71,*)
c
      return
      end
c
c
      subroutine print_R
      implicit none
c
      include 'cfd_osm_com'
      real*8  T_c, RG_impl, RG_qE, v_c, R2_c, fff_expl,
     >        R1_c,  R3_T, R3_nv, R3_c, R3_imp, R3_exp,
     >        rho_c, R2_v, dx
      integer ik,in
      real*8 interpol2p
c
c      do ik = 0, ikmax+1
c         write(71,'(A6,I5,3G16.8)') 'RF/F0',
c     >            ik, RF(1,ik)/F(1,0),
c     >            RF(2,ik)/F(2,0), RF(3,ik)/F(3,0)
c      enddo
c      write(71,*)
c
c      do ik = 0, ikmax+1
c         write(71,'(A6,I5,3G16.8)') 'RG/F0',
c     >            ik, RG(1,ik)/F(1,0),
c     >            RG(2,ik)/F(2,0), RG(3,ik)/F(3,0)
c      enddo
c      write(71,*)
c
c      do ik = 0, ikmax+1
c         write(71,'(A6,I5,3G16.8)') 'RG/RF',
c     >            ik, RG(1,ik)/RF(1,ik),
c     >         RG(2,ik)/RF(2,ik), RG(3,ik)/RF(3,ik)
c      enddo
c      write(71,*)
c
      write(71,'(A5,4A16)') 'ik','R1/F1 [SI]  '
     >    ,'R2/F2 [SI]  ', 'R3/F3 [SI]  '
     >    ,'R4/F4 [SI]  '
      do ik = 1, ikmax
         dx = (xc(ik+1) - xc(ik-1)) / 2.0D0
         write(71,'(I5,4G16.8)') ik,
     >   (R(1,ik) * dx) / F(1,0),
     >   (R(2,ik) * dx) / F(2,0),
     >   (R(3,ik) * dx) / F(3,0),
     >   (R(4,ik) * dx) / F(4,0)
      enddo
      write(71,*)
c
c      write(71,'(A5,3A16)') 'ik','RH1/F1 [SI]  '
c     >    ,'RH2/F2 [SI]  ', 'RH3/F3 [SI]  '
c      do ik = 1, ikmax
c         write(71,'(I5,3G16.8)') ik,
c     >   (RH(1,ik)) / F(1,0),
c     >   (RH(2,ik)) / F(2,0),
c     >   (RH(3,ik)) / F(3,0)
c      enddo
c      write(71,*)
c
c      do ik = 1, ikmax
c         T_c = W(3,ik)/(W(1,ik)*eV)
c         R3_T = 1.5D0 * W(1,ik) * (Ti_eV(ik) - T_c)
c     >             / dtc(ik)
c         R1_c = RF(1,ik)
c         v_c = interpol2p(W(2,ik), W(2,ik+1),
c     >                    xf(ik), xf(ik+1), xc(ik))
c         R2_c = interpol2p(RF(2,ik)+RG(2,ik), RF(2,ik+1)
c     >          + RG(2,ik+1), xf(ik), xf(ik+1), xc(ik))
c         R3_nv = ((1.5D0*T_c/mass - 0.5D0*v_c*v_c) *
c     >              R1_c + v_c * R2_c)
c         R3_imp = R3_T + R3_nv
c         write(71,*) 'R3,R3T,R3nv', ik, RG(3,ik),
c     >      R3_T, R3_nv
c      enddo
c      write(71,*)
c
c      do ik = 1,ikmax
c         T_c = W(3,ik)/(W(1,ik)*eV)
c         R1_c = RF(1,ik)
c         v_c = interpol2p(W(2,ik), W(2,ik+1),
c     >                    xf(ik), xf(ik+1), xc(ik))
c         R2_c = interpol2p(RF(2,ik)+RG(2,ik), RF(2,ik+1)
c     >          + RG(2,ik+1), xf(ik), xf(ik+1), xc(ik))
c         fff_expl = ((1.5D0*T_c/mass - 0.5D0*v_c*v_c) *
c     >              R1_c + v_c * R2_c) * (dtc(ik) /
c     >            (1.5D0 * W(1,ik)))
c         write(71,*) 'f,g,f/g', T_c, fff_expl, fff_expl/T_c
c      enddo
c      write(71,*)
c
c      do ik = 1, ikmax+1
c         R1_c = (RF(1,ik-1) + RF(1,ik)) / 2.0D0
c         rho_c = (U(1,ik-1) + U(1,ik)) / 2.0D0
c         v_c = W(2,ik)
c         R2_v = rho_c * (v(ik) - v_c) / dtf(ik)
c         write(71,*) 'R2G,R2v,R2n,ex-im', ik, RG(2,ik),
c     >      R2_v, R1_c * v_c, RG(2,ik) - R2_v
c      enddo
c      write(71,*)
c
c
c      write(71,'(A5,3A16)') 'ik','F1 [SI]  '
c     >    ,'F2 [SI]  ', 'F3 [SI]  '
c      do ik = 1, ikmax+1
c         write(71,'(I5,3G16.8)') ik,
c     >   (F(1,ik)),
c     >   (F(2,ik)),
c     >   (F(3,ik))
c      enddo
c      write(71,*)
c
c      write(71,'(A5,3A16)') 'ik','G1 [SI]  '
c     >    ,'G2 [SI]  ', 'G3 [SI]  '
c      do ik = 1, ikmax+1
c         write(71,'(I5,3G16.8)') ik,
c     >   (G(1,ik)),
c     >   (G(2,ik)),
c     >   (G(3,ik))
c      enddo
c      write(71,*)
c
c      write(71,'(A5,3A16)') 'ik','H1 [SI]  '
c     >    ,'H2 [SI]  ', 'H3 [SI]  '
c      do ik = 1, ikmax+1
c         write(71,'(I5,3G16.8)') ik,
c     >   (H(1,ik)),
c     >   (H(2,ik)),
c     >   (H(3,ik))
c      enddo
c      write(71,*)
c
c      write(71,'(A5,3A16)') 'ik','J1 [SI]  '
c     >    ,'J2 [SI]  ', 'J3 [SI]  '
c      do ik = 1, ikmax+1
c         write(71,'(I5,3G16.8)') ik,
c     >   (jac(1,1,ik)),
c     >   (jac(1,2,ik)),
c     >   (jac(1,3,ik))
c         write(71,'(I5,3G16.8)') ik,
c     >   (jac(2,1,ik)),
c     >   (jac(2,2,ik)),
c     >   (jac(2,3,ik))
c         write(71,'(I5,3G16.8)') ik,
c     >   (jac(3,1,ik)),
c     >   (jac(3,2,ik)),
c     >   (jac(3,3,ik))
c         write(71,*)
c      enddo
c      write(71,*)
c
c
c      write(71,*) ' ----- A : ------'
c      do ik = 1, ikmax
c         do in = 1, 3
c         write(71,'(I5,4G14.6)') ik,
c     >   aaaa3(in,1,ik), aaaa3(in,2,ik), aaaa3(in,3,ik), ffff3(1,ik)
c         enddo
c      write(71,*)
c      enddo
c      write(71,*)
c
c      write(71,*) ' ----- B: ------'
c      do ik = 1, ikmax
c         do in = 1, 3
c         write(71,'(I5,4G14.6)') ik,
c     >   bbbb3(in,1,ik), bbbb3(in,2,ik), bbbb3(in,3,ik), ffff3(1,ik)
c         enddo
c         write(71,*)
c      enddo
c      write(71,*)
c
c      write(71,*) ' ----- C: ------'
c      do ik = 1, ikmax
c         do in = 1, 3
c         write(71,'(I5,4G14.6)') ik,
c     >   cccc3(in,1,ik), cccc3(in,2,ik), cccc3(in,3,ik), ffff3(1,ik)
c         enddo
c         write(71,*)
c      enddo
c      write(71,*)
c
c      write(71,'(A5,3A16)') 'ik','M1 [SI]  '
c     >    ,'M2 [SI]  ', 'M3 [SI]  '
c      do ik = 1, 4
c         write(71,'(I5,7G14.6)') ik,
c     >   aaaa(1,ik), bbbb(1,ik), cccc(1,ik), ffff(1,ik)
c         write(71,'(I5,7G14.6)') ik,
c     >   aaaa(2,ik), bbbb(2,ik), cccc(2,ik), ffff(2,ik)
c         write(71,'(I5,7G14.6)') ik,
c     >   aaaa(3,ik), bbbb(3,ik), cccc(3,ik), ffff(3,ik)
c      enddo
c      write(71,*)
c
      return
      end
c
c
      subroutine print_mass
      implicit none
c
      include 'cfd_osm_com'
      integer ik,jj,in
      real*8  source,flux1,flux2,sum0,sum,ave2p,temp_1
c
      write(71,*) '   MASS BALANCE : '
      write(71,*)
      write(71,'(A5,5A16)') 'ik','Int(Q1) [SI] '
     >    ,'- F1 [SI]    ','- G1 [SI]    '
     >    ,'Sum [SI] ', 'Error [-]'
      temp_1 = mass * 1.0D23
      sum0 = - F(1,0)
      jj = 0
      source = 0.0D0
      flux1 = - F(1,0)
      flux2 = - G(1,0)
      sum = source + flux1
      write(71,'(I5,5G16.8)') 0, source,
     >      flux1, flux2, sum, (sum-sum0)/sum0
      write(73,'(5G16.8)') xc(0), source,
     >      flux1, flux2, sum
      do ik = 1, ikmax + 1
         source = source + (Q(1,ik-1) + Q_B(1,ik-1))
     >             * (xf(ik) - xf(ik-1))
         flux1 = - F(1,ik)
         flux2 = - G(1,ik)
         sum = source + flux1
         write(71,'(I5,5G16.8)') ik, source,
     >            flux1, flux2, sum, (sum-sum0)/sum0
         write(73,'(5G16.8)') xc(ik)/x_max, source/temp_1,
     >      flux1/temp_1, flux2/temp_1, sum/temp_1
      enddo
c
      write(71,*)
      write(73,*)
      write(71,*) '   GLOBAL MASS FLOW ERROR : ', 1.0D0 -
     >   (abs(F(1,0)) + abs(F(1,ikmax+2))) / source
      write(71,*) '   TARGET MASS FLOW ERROR : ', 1.0D0 -
     >   abs(F(1,0)) / MASS_BC_0,
     >   1.0D0 - abs(F(1,ikmax+2)) / MASS_BC_N
      write(71,*) '   TRAGET MASS FLOW : ', abs(F(1,0))/mass,
     >   MASS_BC_0 / mass,
     >   abs(F(1,ikmax+2)) / mass, MASS_BC_N / mass
c
      return
      end
c
c
      subroutine print_pressure
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
      real*8 source,sum0,flux1,flux2,sum
c
      write(71,*) '   MOMENTUM BALANCE : '
      write(71,*)
      write(71,'(A5,5A16)') 'ik',' Int(Q2) [SI] '
     >    ,'F2 [SI]    ','G2 [SI]    ',
     >     'Sum [SI] ', 'Error [-]'
      sum0 = F(2,0) + G(2,0)
      source = 0.0D0
      flux1 = F(2,0)
      flux2 = G(2,0)
      sum = source + flux1 + flux2
      write(71,'(I5,5G16.8)') 0, -source,
     >   flux1, flux2, sum, (sum-sum0)/sum0
      write(73,'(5G16.8)') xc(0), source,
     >      flux1, flux2, sum
      do ik = 1, ikmax + 1
         source = source - (Q(2,ik-1) + Q_B(2,ik-1))
     >                     * (xf(ik) - xf(ik-1))
         flux1 = F(2,ik)
         flux2 = G(2,ik)
         sum = source + flux1 + flux2
         write(71,'(I5,5G16.8)') ik, - source,
     >   flux1, flux2, sum, (sum-sum0)/sum0
         write(73,'(5G16.8)') xc(ik)/x_max, source,
     >            flux1, flux2, sum
      enddo
c
      write(71,*)
      write(73,*)
      write(71,*) '   GLOBAL MOMENTUM FLOW ERROR = ',
     >     (abs(F(2,0) + G(2,0)) - source -
     >      abs(F(2,ikmax+2) + G(2,ikmax+2))) / MOM_0
      write(71,*) '  MOMENTUM 0,N : ', abs(F(2,0) + G(2,0)),
     >      abs(F(2,ikmax+2) + G(2,ikmax+2))
      write(71,*)
c
      return
      end
c
c
      subroutine print_power
      implicit none
c
      include 'cfd_osm_com'
      integer ik,in
      real*8 source,flux1,flux2,sum,T_c,T_b,ka_ave,T_0,
     >      zhi, sum0, g_net_0, g_net_N, c_s_0, c_s_N,
     >      sum1, sum2, temp_3
c
      temp_3 = 1.0D6
c
      if (a_2T.eq.0) then
c
        write(71,*) '   ENERGY BALANCE :'
        write(71,*)
        write(71,'(A5,5A16)') 'ik','Int(Q3) [SI] '
     >    ,' -F3 [SI]    ','-G3 [SI]    ',
     >     'Sum [SI] ', 'Error [-]'
        c_s_0 = dsqrt(W(3,0) / U(1,0))
        g_net_0 = gamma_sum + 0.5D0 * (W(2,0)/c_s_0)**2.0D0
        c_s_N = dsqrt( W(3,ikmax+1) / U(1,ikmax+1))
        g_net_N = gamma_sum + 0.5D0 * (W(2,ikmax+1)/c_s_N)**2.0D0
        sum0 = - (F(3,0) + G(3,0))
        source = 0.0D0
        flux1 = -F(3,0)
        flux2 = -G(3,0)
        sum = source + flux1 + flux2
        write(71,'(I5,5G16.8)') 0, source,
     >     flux1, flux2, sum, (sum-sum0)/sum0
        write(73,'(5G16.8)') xc(0)/x_max, source,
     >            flux1, flux2, sum
        do ik = 1, ikmax + 1
           source = source + (Q(3,ik-1) + Q_B(3,ik-1))
     >                       * (xf(ik) - xf(ik-1))
           flux1 = - F(3,ik)
           flux2 = - G(3,ik)
           sum = source + flux1 + flux2
           write(71,'(I5,5G16.8)') ik, source,
     >       flux1, flux2, sum, (sum-sum0)/sum0
           write(73,'(5G16.8)') xc(ik)/x_max, source,
     >            flux1, flux2, sum
        enddo
c
        write(71,*)
        write(73,*)
        write(71,*) '   GLOBAL ENERGY FLOW ERROR : ',
     >    1.0D0 - abs(F(3,0) + G(3,0) - F(3,ikmax+2) -
     >            G(3,ikmax+2)) / source
        write(71,*) '   TARGET ENERGY FLOW ERROR : ',
     >    1.0D0 - abs(F(3,0) + G(3,0)) / POW_BC_0,
     >    1.0D0 - abs(F(3,ikmax+2) + G(3,ikmax+2))
     >      / POW_BC_N
        write(71,*) ' HEAT TRANSMISSION BC ERROR : ',
     >    1.0D0 - abs(F(3,0) + G(3,0)) /
     >            dabs(g_net_0 * W(3,0) * W(2,0)) ,
     >    1.0D0 - abs(F(3,ikmax+2) + G(3,ikmax+2)) /
     >        dabs(g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW 0:  A/C, B/C, A/B',
     >     (F(3,0) + G(3,0)) / POW_BC_0,
     >     (g_net_0 * W(3,0) * W(2,0)) / POW_BC_0,
     >     (F(3,0) + G(3,0)) / (g_net_0 * W(3,0) * W(2,0))
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW N:  A/C, B/C, A/B',
     >     (F(3,ikmax+2) + G(3,ikmax+2)) / POW_BC_N,
     >     (g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
     >      / POW_BC_N,  (F(3,ikmax+2) + G(3,ikmax+2))
     >      / (g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
        write(71,*)
c
       elseif (a_2T.eq.1) then
c
        write(71,*) '   ION ENERGY BALANCE :'
        write(71,*)
        write(71,'(A5,5A16)') 'ik','Int(Q3) [SI] '
     >    ,'  F3 [SI]    ',' G3 [SI]    ',
     >     'Sum [SI] ', 'Error [-]'
        c_s_0 = dsqrt((W(3,0) + W(4,0)) / U(1,0))
        g_net_0 = gamma_i + 0.5D0 * (W(2,0)/c_s_0)**2.0D0 *
     >             (1.0D0 + W(4,0) / W(3,0))
        c_s_N = dsqrt((W(3,ikmax+1) + W(4,ikmax+1)) / U(1,ikmax+1))
        g_net_N = gamma_i + 0.5D0 * (W(2,ikmax+1)/c_s_N)**2.0D0 *
     >             (1.0D0 + W(4,ikmax+1) / W(3,ikmax+1))
        sum0 = - (F(3,0) + G(3,0))
        source = 0.0D0
        flux1 = -F(3,0)
        flux2 = -G(3,0)
        sum = source + flux1 + flux2
        write(71,'(I5,5G16.8)') 0, source,
     >     flux1, flux2, sum, (sum-sum0)/sum0
        write(73,'(5G16.8)') xc(0)/x_max, source/temp_3,
     >            flux1/temp_3, flux2/temp_3, sum/temp_3
        do ik = 1, ikmax + 1
           source = source + (Q(3,ik-1) + Q_ei(ik-1)
     >      + Q_joule(ik-1) + Q_B(3,ik-1)) * (xf(ik) - xf(ik-1))
           flux1 = - F(3,ik)
           flux2 = - G(3,ik)
           sum = source + flux1 + flux2
           write(71,'(I5,5G16.8)') ik, source,
     >       F(3,ik), G(3,ik), sum, (sum-sum0)/sum0
           write(73,'(5G16.8)') xc(ik)/x_max, source/temp_3,
     >            flux1/temp_3, flux2/temp_3, sum/temp_3
        enddo
c
        write(71,*)
        write(73,*)
        write(71,*) '   GLOBAL ION ENERGY FLOW ERROR : ',
     >    1.0D0 - abs(F(3,0) + G(3,0) - F(3,ikmax+2) -
     >            G(3,ikmax+2)) / source
        write(71,*) '   TARGET ENERGY FLOW ERROR : ',
     >    1.0D0 - abs(F(3,0) + G(3,0)) / POW_BC_i0,
     >    1.0D0 - abs(F(3,ikmax+2) + G(3,ikmax+2))
     >      / POW_BC_iN
        write(71,*) ' HEAT TRANSMISSION BC ERROR : ',
     >    1.0D0 - abs(F(3,0) + G(3,0)) /
     >            dabs(g_net_0 * W(3,0) * W(2,0)),
     >    1.0D0 - abs(F(3,ikmax+2) + G(3,ikmax+2)) /
     >        dabs(g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW 0:  A/C, B/C, A/B',
     >     (F(3,0) + G(3,0)) / POW_BC_i0,
     >     (g_net_0 * W(3,0) * W(2,0)) / POW_BC_i0,
     >     (F(3,0) + G(3,0)) / (g_net_0 * W(3,0) * W(2,0))
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW N:  A/C, B/C, A/B',
     >     (F(3,ikmax+2) + G(3,ikmax+2)) / POW_BC_iN,
     >     (g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
     >      / POW_BC_iN,  (F(3,ikmax+2) + G(3,ikmax+2))
     >      / (g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
        write(71,*)
c
        write(71,*) '   ELECTRON ENERGY BALANCE :'
        write(71,*)
        write(71,'(A5,5A16)') 'ik','Int(Q4) [SI] '
     >    ,'  F4 [SI]    ',' G4 [SI]    ',
     >     'Sum [SI] ', 'Error [-]'
        c_s_0 = dsqrt((W(3,0) + W(4,0)) / U(1,0))
        g_net_0 = gamma_e
        c_s_N = dsqrt((W(3,ikmax+1) + W(4,ikmax+1)) / U(1,ikmax+1))
        g_net_N = gamma_e
        sum0 = - (F(4,0) + G(4,0))
        flux1 = - F(4,0)
        flux2 = - G(4,0)
        source = 0.0D0
        sum1 = 0.0
        sum2 = 0.0
        sum = source + flux1 + flux2
        write(71,'(I5,5G16.8)') 0, source,
     >     flux1, flux2, sum, (sum-sum0)/sum0
        write(73,'(5G16.8)') xc(0)/x_max, source/temp_3,
     >            flux1/temp_3, flux2/temp_3, sum/temp_3
        do ik = 1, ikmax + 1
           source = source + (Q(4,ik-1) - Q_ei(ik-1) -
     >       Q_joule(ik-1) + Q_B(4,ik-1)) * (xf(ik) - xf(ik-1))
           flux1 = - F(4,ik)
           flux2 = - G(4,ik)
           sum = source + flux1 + flux2
           sum1 = sum1 +
     >            (abs(F(4,ik)) + abs(F(3,ik)) + abs(G(3,ik)))
     >                   * (xf(ik) - xf(ik-1))
           sum2 = sum2 + abs(G(4,ik))
     >                   * (xf(ik) - xf(ik-1))
           write(71,'(I5,5G16.8)') ik, source,
     >       F(4,ik), G(4,ik), sum, (sum-sum0)/sum0
           write(73,'(5G16.8)') xc(ik)/x_max, source/temp_3,
     >            flux1/temp_3, flux2/temp_3, sum/temp_3
        enddo
c
        write(71,*)
        write(73,*)
        write(71,*) 'f_cond:', sum2 / (sum1 + sum2)
        write(71,*) '   GLOBAL ELECTRON ENERGY FLOW ERROR : ',
     >    1.0D0 - abs(F(4,0) + G(4,0) - F(4,ikmax+2) -
     >            G(4,ikmax+2)) / source
        write(71,*) '   TARGET ENERGY FLOW ERROR : ',
     >    1.0D0 - abs(F(4,0) + G(4,0)) / POW_BC_e0,
     >    1.0D0 - abs(F(4,ikmax+2) + G(4,ikmax+2))
     >      / POW_BC_eN
        write(71,*) ' HEAT TRANSMISSION BC ERROR : ',
     >    1.0D0 - abs(F(4,0) + G(4,0)) /
     >            dabs(g_net_0 * W(4,0) * W(2,0)),
     >    1.0D0 - abs(F(4,ikmax+2) + G(4,ikmax+2)) /
     >        dabs(g_net_N * W(4,ikmax+1) * W(2,ikmax+1))
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW 0:  A/C, B/C, A/B',
     >     (F(4,0) + G(4,0)) / POW_BC_e0,
     >     (g_net_0 * W(4,0) * W(2,0)) / POW_BC_e0,
     >     (F(4,0) + G(4,0)) / (g_net_0 * W(4,0) * W(2,0))
        write(71,'(A30,3G16.8)')
     >     '  HEAT FLOW N:  A/C, B/C, A/B',
     >     (F(4,ikmax+2) + G(4,ikmax+2)) / POW_BC_eN,
     >     (g_net_N * W(4,ikmax+1) * W(2,ikmax+1))
     >      / POW_BC_eN,  (F(4,ikmax+2) + G(4,ikmax+2))
     >      / (g_net_N * W(4,ikmax+1) * W(2,ikmax+1))
        write(71,*)
c
      endif
c
      return
      end
c
c
      subroutine check_error
      implicit none
c
      include 'cfd_osm_com'
      integer ik,jj,in, a_er_1, a_er_2, a_er_3
      real*8  source,flux,sum0,sum,error,max_error,
     >   flux1, flux2, x_error, T_c, T_b, ka_ave, T_0,
     >      zhi, g_net_0, g_net_N, c_s_0, c_s_N,
     >      err_1a, err_1b, err_3a, err_3b
c
      if (rms_res.le.rms_tol.and.max_res.le.max_tol) then
        a_resid = 0
      elseif (rms_res.ge.1.0.or.max_res.ge.1.0) then
        a_resid = 2
      else
        a_resid = 1
      endif
c
c      c_s_0 = dsqrt(W(3,0) / U(1,0))
c      g_net_0 = gamma_sum + 0.5D0 * (W(2,0)/c_s_0)**2.0D0
c      c_s_N = dsqrt( W(3,ikmax+1) / U(1,ikmax+1))
c      g_net_N = gamma_sum + 0.5D0 * (W(2,ikmax+1)/c_s_N)**2.0D0
c
c      err_3a = 1.0 - (F(3,0) + G(3,0)) /
c     >               (g_net_0 * W(3,0) * W(2,0))
c      err_3b = 1.0 - (F(3,ikmax+2) + G(3,ikmax+2))
c     >      / (g_net_N * W(3,ikmax+1) * W(2,ikmax+1))
c
c      if (abs(err_3a).le.max_tol.and.abs(err_3b).le.max_tol) then
c        a_er_3 = 0
c      else
c        a_er_3 = 1
c      endif
c
      if (a_resid.eq.0) then
        a_error = 0
      else
        a_error = 1
      endif
c
      return
      end
c
c
      subroutine print_error
      implicit none
c
      include 'cfd_osm_com'
      integer ik,jj,in
      real*8  source,flux,sum0,sum,error,max_error,
     >   flux1, flux2, x_error, T_c, T_b, ka_ave, T_0,
     >      zhi, g_net_0, g_net_N, c_s_0, c_s_N
c
c      if (rms_res.le.rms_tol.and.max_res.le.max_tol) then
c        a_resid = 0
c      else
c        a_resid = 1
c      endif
c
      max_error = 0.0D0
      sum0 = - F(1,0)
      source = 0.0D0
      flux = - F(1,0)
      sum = source + flux
c
      do ik = 1, ikmid
         source = source + Q(1,ik-1) * (xf(ik) - xf(ik-1))
         flux = - F(1,ik)
         sum = source + flux
         error = abs((sum - sum0) / sum0)
         if (error.ge.max_error) then
            max_error = error
            x_error   = xc(ik)
         endif
         RH(1,ik) = abs(error)
         RH(1,ikmax+1-ik) = RH(1,ik)
      enddo
      RH(1,0) = 0.0
      RH(1,ikmax+2) = 0.0
      write(71,'(A30,2G16.8)') 'MAX MASS ERROR, x[m] ',
     >     max_error
c      write(0,'(A30,2G16.8)') 'MAX MASS ERROR, x[m] ',
c     >     max_error
c
c      if (max_error.ge.0.1D0) then
c         a_error = 2
c      elseif (max_error.ge.max_tol) then
c         a_error = 1
c      else
c         a_error = 0
c      endif
c
      max_error = 0.0D0
      sum0 = F(2,0) + G(2,0)
      source = 0.0D0
      flux1 = F(2,0)
      flux2 = G(2,0)
      sum = source + flux1 + flux2
      do ik = 1, ikmid
         source = source - Q(2,ik-1) * (xf(ik) - xf(ik-1))
         flux1 = F(2,ik)
         flux2 = G(2,ik)
         sum = source + flux1 + flux2
         error = abs((sum - sum0) / sum0)
         if (error.ge.max_error) then
            max_error = error
            x_error   = xc(ik)
         endif
         RH(2,ik) = abs(error)
         RH(2,ikmax+1-ik) = RH(2,ik)
      enddo
      RH(2,0) = 0.0
      RH(2,ikmax+2) = 0.0
      write(71,'(A30,2G16.8)') 'MAX MOM. ERROR, x[m] ',
     >     max_error
c      write(0,'(A30,2G16.8)') 'MAX MOM. ERROR, x[m] ',
c     >     max_error
c
c      if (max_error.ge.0.1D0) then
c         a_error = 2
c      elseif (max_error.ge.max_tol) then
c         a_error = 1
c      else
c         a_error = 0
c      endif
c
      if (a_2T.eq.0) then
        max_error = 0.0D0
        c_s_0 = dsqrt(W(3,0) / U(1,0))
        g_net_0 = gamma_sum + 0.5D0 * (W(2,0)/c_s_0)**2.0D0
        c_s_N = dsqrt( W(3,ikmax+1) / U(1,ikmax+1))
        g_net_N = gamma_sum + 0.5D0 * (W(2,ikmax+1)/c_s_N)**2.0D0
        sum0 = - (F(3,0) + G(3,0))
        source = 0.0D0
        flux1 = -F(3,0)
        flux2 = -G(3,0)
        sum = source + flux1 + flux2
        do ik = 1, ikmid
           source = source + Q(3,ik-1) * (xf(ik) - xf(ik-1))
           flux1 = - F(3,ik)
           flux2 = - G(3,ik)
           sum = source + flux1 + flux2
           error = abs((sum - sum0) / sum0)
           if (error.ge.max_error) then
              max_error = error
              x_error   = xc(ik)
           endif
           RH(3,ik) = abs(error)
           RH(3,ikmax+1-ik) = RH(3,ik)
        enddo
        RH(3,0) = 0.0
        RH(3,ikmax+2) = 0.0
        write(71,'(A30,2G16.8)') 'MAX POW. ERROR, x[m] ',
     >     max_error
c        write(0,'(A30,2G16.8)') 'MAX POW. ERROR, x[m] ',
c     >     max_error
      endif
c
c      if (max_error.ge.0.1D0) then
c         a_error = 2
c      elseif (max_error.ge.max_tol) then
c         a_error = 1
c      else
c         a_error = 0
c      endif
c
      return
      end
c
c
      subroutine wrtcfdbg(irlim1,irlim2)
      implicit none
c
      integer irlim1, irlim2
      include 'params'
      include 'sol23_com'
      include 'cgeom'
      include 'pin_cfd'
c
      integer of
      parameter(of=75)
      integer ik,ir,in,id
c
      rewind(of)
c
      write (of,10) 'CFD BACKGROUND PLASMA:'
      write (of,20) irlim1, irlim2
      write (of,10) 'INDEX:'
      write (of,400) ((ikmax_cfd(ir,id),ir=irlim1,irlim2),id=1,2)
c
      write (of,10) 'FINISHED:'
      write (of,400) (finished(ir),ir=irlim1,irlim2)
c
      write (of,10) 'MASS:'
      write (of,500) ((part_flux(ir,id),ir=irlim1,irlim2),id=1,2)
c
      write (of,10) 'MOMENTUM:'
      write (of,500) ((mom_flux(ir,id),ir=irlim1,irlim2),id=1,2)
c
      write (of,10) 'POWER:'
      write (of,500) ((power_flux(ir,id),ir=irlim1,irlim2),id=1,2)
c
      write (of,10) 'TARG MASS:'
      write (of,500) ((gamma_targ(ir,id),ir=irlim1,irlim2),id=1,2)
c
      write (of,10) 'TARG MOM:'
      write (of,500) ((mom_targ(ir,id),ir=irlim1,irlim2),id=1,2)
c
      if (s23_2T.eq.0) then
        write (of,10) 'TARG POW:'
        write (of,500) ((power_targ(ir,id),ir=irlim1,irlim2),id=1,2)
      else
        write (of,10) 'TARG POWi:'
        write (of,500) ((poweri_targ(ir,id),ir=irlim1,irlim2),id=1,2)
        write (of,10) 'TARG POWe:'
        write (of,500) ((powere_targ(ir,id),ir=irlim1,irlim2),id=1,2)
        write (of,10) 'POW ei:'
        write (of,500) (power_ei(ir),ir=irlim1,irlim2)
      endif
c
      write (of,10) 'POW LIM:'
      write (of,500) ((s23_powlim(ir,id),ir=irlim1,irlim2),id=1,2)
c
      write (of,10) 'MOM LIM:'
      write (of,500) ((s23_momlim(ir,id),ir=irlim1,irlim2),id=1,2)
c
      write (of,10) 'MASS LIM:'
      write (of,500) ((s23_masslim(ir,id),ir=irlim1,irlim2),id=1,2)
c
      if (s23_2T.eq.0) then
        write (of,10) 'FRACT:'
        write (of,500) (((fract_cfd(ir,id,in),ir=irlim1,irlim2)
     >                             ,id=1,2),in=1,3)
      else
        write (of,10) 'FRACT:'
        write (of,500) (((fract_cfd(ir,id,in),ir=irlim1,irlim2)
     >                             ,id=1,2),in=1,4)
      endif
c
      write (of,10)  'GRID:'
      write (of,500) (((xc_cfd(ir,id,ik),ik=0,ikmax_cfd(ir,id)+1),
     >                             ir=irlim1,irlim2),id=1,2)
c
      if (s23_2T.eq.0) then
        write (of,10)  'U_CFD:'
        write (of,500) ((((U_cfd(ir,id,in,ik),ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,3)
        write (of,10)  'Q_CFD:'
        write (of,500) ((((Q_cfd(ir,id,in,ik),ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,3)
      else
        write (of,10)  'U_CFD:'
        write (of,500) ((((U_cfd(ir,id,in,ik),ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,4)
        write (of,10)  'Q_CFD:'
        write (of,500) ((((Q_cfd(ir,id,in,ik),ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,4)
      endif
c
      write (of,10) 'PINION:'
      write (of,500) ((pinion_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'PINREC:'
      write (of,500) ((pinrec_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'PINMOM:'
      write (of,500) ((pinmom_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'PINQI:'
      write (of,500) ((pinqi_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'PINQE:'
      write (of,500) ((pinqe_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'PINATOM:'
      write (of,500) ((pinatom_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'PINMOL:'
      write (of,500) ((pinmol_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'PINENA:'
      write (of,500) ((pinena_last(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'KNBS:'
      write (of,500) ((knbs(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'KTIBS:'
      write (of,500) ((ktibs(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'KTEBS:'
      write (of,500) ((ktebs(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'KVHS:'
      write (of,500) ((kvhs(ik,ir),ik=1,nks(ir)),
     >                   ir=irlim1,irlim2)
c
      write (of,10) 'TE_CTRL:'
      write (of,500) ((Te_ctrl(ir,id),ir=irlim1,irlim2),
     >                  id=1,2)
c
      write (of,10) 'NE_TARG:'
      write (of,500) ((old_ne(ir,id),ir=irlim1,irlim2),
     >                  id=1,2)
c
      write (of,10) 'TE_TARG:'
      write (of,500) ((old_te(ir,id),ir=irlim1,irlim2),
     >                  id=1,2)
c
      write (of,10) 'TI_TARG:'
      write (of,500) ((old_ti(ir,id),ir=irlim1,irlim2),
     >                  id=1,2)
c
      write (of,10) 'V_TARG:'
      write (of,500) ((old_v(ir,id),ir=irlim1,irlim2),
     >                  id=1,2)
c
      write (of,10) 'MACH_TARG:'
      write (of,500) ((old_mach(ir,id),ir=irlim1,irlim2),
     >                  id=1,2)
c
      return
c
  10  format(a)
  20  format(2i6)
 400  format(12i6)
 500  format(6d18.10)
c
      end
c
c
      subroutine readcfdbg(irlim1,irlim2)
      implicit none
c
      integer irlim1, irlim2, ir1, ir2
      include 'params'
      include 'sol23_com'
      include 'cgeom'
      include 'comtor'
      include 'sol23_input'
      include 'pindata'
      include 'pin_cfd'
c
      integer infile
      parameter(infile=74)
c
      character*120 buffer
      integer ik,ir,id,in
      integer tmpnrs,tmpnds,tmpirsep
      integer tmpnks(maxnrs)
c
      rewind(infile)
      read(infile,10,end=1500,err=1500) buffer
c
      rewind(infile)
      read(infile,10,end=1000,err=2000) buffer
c
      if (buffer(1:8).ne.'CFD BACK') then
         write (71,*) 'ERROR: NOT A CFD PLASMA FILE'
         write (71,*) buffer
         STOP
      endif
c
      read(infile,20,end=1000,err=2000) ir1, ir2
c
      if (ir1.ne.irlim1.or.ir2.ne.irlim2) then
         write (71,*) 'ERROR: IR.ne.IRLIM'
         write (71,*) ir1,irlim1,ir2,irlim2
      endif
c
 50   read(infile,10,end=1000,err=2000) buffer
c
      if (buffer(1:6).eq.'INDEX:') then
c
         read (infile,400) ((ikmax_cfd(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:9).eq.'FINISHED:') then
c
         read (infile,400) (finished(ir),ir=irlim1,irlim2)
c
      elseif (buffer(1:5).eq.'MASS:') then
c
         read (infile,500) ((part_flux(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:9).eq.'MOMENTUM:') then
c
         read (infile,500) ((mom_flux(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:6).eq.'POWER:') then
c
         read (infile,500) ((power_flux(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:10).eq.'TARG MASS:') then
c
         read (infile,500) ((gamma_targ(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:9).eq.'TARG MOM:') then
c
         read (infile,500) ((mom_targ(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:9).eq.'TARG POW:') then
c
         read (infile,500) ((power_targ(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:10).eq.'TARG POWi:') then
c
         read (infile,500) ((poweri_targ(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:10).eq.'TARG POWe:') then
c
         read (infile,500) ((powere_targ(ir,id),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:7).eq.'POW ei:') then
c
         read (infile,500) (power_ei(ir),
     >                           ir=irlim1,irlim2)
c
      elseif (buffer(1:8).eq.'POW LIM:') then
c
         read (infile,500) ((s23_powlim(ir,id),ir=irlim1,irlim2)
     >                           ,id=1,2)
c
c      elseif (buffer(1:8).eq.'MOM LIM:') then
c
c         write(0,*) 'MOM'
c         read (infile,500) ((s23_momlim(ir,id),ir=irlim1,irlim2)
c     >                            ,id=1,2)
c
c      elseif (buffer(1:9).eq.'MASS LIM:') then
c
c         write(0,*) 'MASS'
c         read (infile,500) ((s23_masslim(ir,id),ir=irlim1,irlim2)
c     >                            ,id=1,2)
c
      elseif (buffer(1:6).eq.'FRACT:') then
c
         if (s23_2T.eq.0) then
           read (infile,500) (((fract_cfd(ir,id,in),
     >                   ir=irlim1,irlim2),id=1,2),in=1,3)
         else
           read (infile,500) (((fract_cfd(ir,id,in),
     >                   ir=irlim1,irlim2),id=1,2),in=1,4)
         endif
c
      elseif (buffer(1:5).eq.'GRID:') then
c
         read (infile,500) (((xc_cfd(ir,id,ik),
     >                           ik=0,ikmax_cfd(ir,id)+1),
     >                           ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:6).eq.'U_CFD:') then
c
         if (s23_2T.eq.0) then
           read (infile,500) ((((U_cfd(ir,id,in,ik),
     >                     ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,3)
         else
           read (infile,500) ((((U_cfd(ir,id,in,ik),
     >                     ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,4)
         endif
c
      elseif (buffer(1:6).eq.'Q_CFD:') then
c
         if (s23_2T.eq.0) then
           read (infile,500) ((((Q_cfd(ir,id,in,ik),
     >                     ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,3)
         else
           read (infile,500) ((((Q_cfd(ir,id,in,ik),
     >                     ik=0,ikmax_cfd(ir,id)+1),
     >                     ir=irlim1,irlim2),id=1,2),in=1,4)
         endif
c
      elseif (buffer(1:7).eq.'PINION:') then
c
         read (infile,500) ((pinion_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinion', pinion_last(1,8)
c
      elseif (buffer(1:7).eq.'PINREC:') then
c
         read (infile,500) ((pinrec_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinrec', pinrec_last(1,8)
c
      elseif (buffer(1:7).eq.'PINMOM:') then
c
         read (infile,500) ((pinmom_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinmom', pinmom_last(1,8)
c
      elseif (buffer(1:6).eq.'PINQI:') then
c
         read (infile,500) ((pinqi_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinqi', pinqi_last(1,8)
c
      elseif (buffer(1:6).eq.'PINQE:') then
c
         read (infile,500) ((pinqe_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinqe', pinqe_last(1,8)
c
      elseif (buffer(1:8).eq.'PINATOM:') then
c
         read (infile,500) ((pinatom_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinatom', pinatom_last(1,8)
c
      elseif (buffer(1:7).eq.'PINMOL:') then
c
         read (infile,500) ((pinmol_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinmol', pinmol_last(1,8)
c
      elseif (buffer(1:7).eq.'PINENA:') then
c
         read (infile,500) ((pinena_last(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c         write(0,*) 'pinena', pinena_last(1,8)
c
      elseif (buffer(1:5).eq.'KNBS:') then
c
         read (infile,500) ((oldknbs(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c
      elseif (buffer(1:7).eq.'KTIBS:') then
c
         read (infile,500) ((oldktibs(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c
      elseif (buffer(1:6).eq.'KTEBS:') then
c
         read (infile,500) ((oldktebs(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c
      elseif (buffer(1:5).eq.'KVHS:') then
c
         read (infile,500) ((oldkvhs(ik,ir),
     >                  ik=1,nks(ir)),ir=irlim1,irlim2)
c
      elseif (buffer(1:8).eq.'TE_CTRL:') then
c
         read (infile,500) ((Te_ctrl(ir,id),
     >                  ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:8).eq.'NE_TARG:') then
c
         read (infile,500) ((old_ne(ir,id),
     >                  ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:8).eq.'TE_TARG:') then
c
         read (infile,500) ((old_te(ir,id),
     >                  ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:8).eq.'TI_TARG:') then
c
         read (infile,500) ((old_ti(ir,id),
     >                  ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:7).eq.'V_TARG:') then
c
         read (infile,500) ((old_v(ir,id),
     >                  ir=irlim1,irlim2),id=1,2)
c
      elseif (buffer(1:10).eq.'MACH_TARG:') then
c
         read (infile,500) ((old_mach(ir,id),
     >                  ir=irlim1,irlim2),id=1,2)
c
      endif
c
      goto 50
c
 1000 continue
      s23_cfdfile = 1
      seedplasma  = 0
c      write(0,*) 'READ CFD FILE'
      write(71,*) 'READ CFD FILE'
c
      do ir = irlim1, irwall-1
        do ik = 1, nks(ir)
          knbs(ik,ir) = oldknbs(ik,ir)
          ktibs(ik,ir) = oldktibs(ik,ir)
          ktebs(ik,ir) = oldktebs(ik,ir)
          kvhs(ik,ir) = oldkvhs(ik,ir)
          pinion(ik,ir) = pinion_last(ik,ir)
          pinrec(ik,ir) = pinrec_last(ik,ir)
          pinmom(ik,ir) = pinmom_last(ik,ir)
          pinqi(ik,ir) = pinqi_last(ik,ir)
          pinqe(ik,ir) = pinqe_last(ik,ir)
        enddo
        kteds(idds(ir,1)) = old_te(ir,1)
        ktids(idds(ir,1)) = old_ti(ir,1)
        knds(idds(ir,1))  = old_ne(ir,1)
        kvds(idds(ir,1))  = old_v(ir,1)
        cmachno(ir,1)     = old_mach(ir,1)
        kteds(idds(ir,2)) = old_te(ir,2)
        ktids(idds(ir,2)) = old_ti(ir,2)
        knds(idds(ir,2))  = old_ne(ir,2)
        kvds(idds(ir,2))  = old_v(ir,2)
        cmachno(ir,2)     = old_mach(ir,2)
      enddo
c
      do ir = irwall, irlim2
        old_te(ir,1)   = kteds(idds(ir,1))
        old_ti(ir,1)   = ktids(idds(ir,1))
        old_ne(ir,1)   = knds(idds(ir,1))
        old_v(ir,1)    = kvds(idds(ir,1))
        old_mach(ir,1) = cmachno(ir,1)
c
        old_te(ir,2)   = kteds(idds(ir,2))
        old_ti(ir,2)   = ktids(idds(ir,2))
        old_ne(ir,2)   = knds(idds(ir,2))
        old_v(ir,2)    = kvds(idds(ir,2))
        old_mach(ir,2) = cmachno(ir,2)
      enddo
c
c      do ir = irlim1, irlim2
c        do ik = 0, ikmax_cfd(ir,1) + 1
c          if (U_cfd(ir,1,4,ik).eq.0.0D0) then
c             U_cfd(ir,1,3,ik) = U_cfd(ir,1,3,ik) / 2.0D0
c             U_cfd(ir,1,4,ik) = U_cfd(ir,1,3,ik)
c             Q_cfd(ir,1,3,ik) = Q_cfd(ir,1,3,ik)
c             Q_cfd(ir,1,4,ik) = Q_cfd(ir,1,3,ik)
c          endif
c        enddo
c      enddo
c
      do ir = irlim1, irlim2
        do id = 1, 2
          del_power_flux(ir,id) = 0.0D0
          del_mom_flux(ir,id)   = 0.0D0
c
          if (gamma_targ(ir,id).eq.0.0) then
            gamma_targ(ir,id) = gamma_Lang(ir,id)
          endif
c
          if (mom_targ(ir,id).eq.0.0) then
            mom_targ(ir,id) = mom_Lang(ir,id)
          endif
c
          if (poweri_targ(ir,id).eq.0.0) then
            poweri_targ(ir,id) = poweri_Lang(ir,id)
          endif
c
          if (powere_targ(ir,id).eq.0.0) then
            powere_targ(ir,id) = powere_Lang(ir,id)
          endif
c
          if (power_targ(ir,id).eq.0.0) then
            power_targ(ir,id) = poweri_targ(ir,id) +
     >                          powere_targ(ir,id)
          endif
        enddo
      enddo
c
      return
c
 1500 continue
      s23_cfdfile = 0
c      write(0,*) 'NO CFD FILE'
      return
c
 2000 continue
      write (71,*) 'ERROR READING IN CFD PLASMA BACKGROUND:'
      STOP
c
  10  format(a)
  20  format(2i6)
 400  format(12i6)
 500  format(6d18.10)
c
      end
c
c
      subroutine Te_inv(irlim1,irlim2,ik_mid,flag)
      implicit none
c
      integer irlim1, irlim2, ik_mid, flag, ir, ik
      include 'params'
      include 'sol23_com'
      include 'cgeom'
      include 'comtor'
      include 'sol23_input'
c
      flag = 0
c
      do ir = irlim1, irlim2
        do ik = 1, ik_mid
          if (ktebs(ik,ir).lt.ktebs(1,ir)) then
            flag = 1
            write(71,*) 'Te_flag=1',
     >             ir,ik,ktebs(ik,ir),ktebs(1,ir)
          endif
        enddo
c
        do ik = ik_mid + 1, nks(ir)
          if (ktebs(ik,ir).lt.ktebs(nks(ir),ir)) then
            flag = 1
            write(71,*) 'Te_flag=1',
     >             ir,ik,ktebs(ik,ir),ktebs(nks(ir),ir)
          endif
        enddo
c
c        write(0,*) 'ir,f,Te:', ir, flag, ktebs(1,ir)
c
      enddo
c
      return
      end
c
c
c      real*8 function phi_lim(r)
c      implicit none
c      real*8 r, b
c
c      b = 1.5D0
c      phi_lim = max(0.0D0, min(b*r, 1.0D0), min(r, b))
c
c      return
c      end
c
c
c      subroutine U_Riemann
c      implicit none
c      real*8 root_ff, ff, dt
c
c      include 'cfd_ps1_com'
c      include 'cfd_ps2_com'
c
c      gamma = 1.4D0
c      alpha = (gamma + 1.0D0) / (gamma - 1.0D0)
c      p_L   = WL(3)
c      p_R   = WR(3)
c      rho_L = WL(1)
c      rho_R = WR(1)
c      u_L   = WL(2)
c      u_R   = WR(2)
c
c      p_ratio = root_ff(0.0D0, 10.0D0, tol)
c
c      u_cont = u_L + 2.0D0 / (gamma - 1.0D0) * c_L * (1.0D0 -
c     >  (p_R * p_ratio / p_L) ** ((gamma - 1.0D0)/(2.0D0*gamma)))
c
c      WL(3) = p_R * p_ratio
c      WL(1) = rho_L * (WL(3) / p_L) ** (1.0D0 / gamma)
c      WL(2) = u_cont
c
c      write(6,'(3G14.6)') WL(1), WL(2), WL(3)
c      write(0,'(3G14.6)') WL(1), WL(2), WL(3)
c
c      return
c      end
c
c
c      subroutine WfromULR
c      implicit none
c
c      include 'cfd_ps2_com'
c
c      WL(1) = UL(1)
c      WL(2) = UL(2) / UL(1)
c      WL(3) = (g_pv - 1.0D0) * (UL(3)
c     >               - 0.5D0 * UL(2)** 2.0D0/ UL(1))
c      WR(1) = UR(1)
c      WR(2) = UR(2) / UR(1)
c      WR(3) = (g_pv - 1.0D0) * (UR(3)
c     >               - 0.5D0 * UR(2)** 2.0D0/ UR(1))
c
c      return
c      end
c
c
c      subroutine ULfromWL
c      implicit none
c
c      include 'cfd_ps2_com'
c
c      UL(1) = WL(1)
c      UL(2) = WL(2) * WL(1)
c      UL(3) = WL(3) / (g_pv - 1.0D0)
c     >           + 0.5D0 * WL(2)** 2.0D0 * WL(1)
c      return
c      end
c
c
c      real*8 function root_ff(x_a,x_b,tol)
c      implicit none
c      real*8 x_a, x_b, tol
c
c      real*8 x, ff_x, ff_a, ff_b, ff
c      integer nr_it
c
c      nr_it = 0
c
c 100  x = (x_a + x_b) / 2.0D0
c      ff_a = ff(x_a)
c      ff_b = ff(x_b)
c      ff_x = ff(x)
c
c      if  ((ff_a.le.0.0D0.and.ff_x.le.0.0D0).or.
c     >     (ff_a.ge.0.0D0.and.ff_x.ge.0.0D0)) then
c         x_a = x
c      else
c         x_b = x
c      endif
c
c      if (dabs(ff_x).ge.tol) then
c         write(0,'(A10,4G16.8)') 'xa,x,xb,ff', x_a,x,x_b,ff_x
c         nr_it = nr_it + 1
c         if (nr_it.gt.100) then
c            write(0,*) 'NO SOLN in ROOTFF'
c            STOP
c         endif
c         goto 100
c      endif
c      write(6,*)
c
c      root_ff = x
c
c      return
c      end
c
