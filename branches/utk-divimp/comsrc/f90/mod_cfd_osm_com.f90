module mod_cfd_osm_com
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  real*8,public :: mass_p,mass_e,ev,c_v_0,c_p_0,g_pv,cons_gpv
  !
  !      parameter (mass_p = 1.0d0, mass_e = 0.0d0, ev = 1.0d0, c_v_0 = 1.5d0,
  !     >    c_p_0 = 2.5d0, g_pv = 1.4d0)
  !
  integer,public :: nnn
  !
  parameter (mass_p = 1.6726d-27, ev = 1.6022d-19,mass_e = 9.1095d-31, c_p_0 = 2.5d0,&
        c_v_0 = 1.5d0,g_pv = c_p_0 / c_v_0, nnn = 300,cons_gpv =  1.0d0 - 2.0d0 / (g_pv - 1.0d0))
  ! common /cfd_osm_com/ unit,u,u0,q0,w,f,g,h,q,q_new,q_perp, tau_i, tau_e, q_b, b_field,&
  !      l_b, q_grad, q_uni,q_ei, q_joule, tau_ei, lam_mat, t_mat, tinv_mat,ep2, ep4,&
  !      psi_p, aaaa, bbbb, cccc, ffff, ffff_new,rf,rg,rh,r, re, xc, xf, dxc, dtc, dtf,&
  !      a_q2, a_q3,ka, ka_i, ka_e, mu, v_j, c_s_j, nu_max,aaa, bbb, ccc, ddd, fff, uuu,&
  !     aaaa3, bbbb3, cccc3, ffff3, uuuu3, jac_a3, jac_b3,aaaa4, bbbb4, cccc4, ffff4, uuuu4,&
  !      jac_a4, jac_b4, jac_q4,v, ti_ev, te_ev, a_rk, a_artvisc4, a_artvisc2, a_temp,&
  !     a_grid0, a_gridexp, max_err, rms_err, err_tol, l_pow,l_iz, l_iz_0, l_iz_n, l_mom,&
  !      l_en,mu0, ka0_ev, a_ex, dx_ave, dx0_bc, dxn_bc, x_max,max_res, rms_res, q1,&
  !      q2, q3, q4, mass, mass_i, l_izoffset,a_i, z_i, c_v, c_p, a_dt_f, a_dt_g, n0_bc,&
  !      nn_bc,ti0_bc, tin_bc, te0_bc, ten_bc, m0_bc, mn_bc, f0_bc, fn_bc,p_0, p_n, v_0,&
  !      v_n, x_xpt_0, x_xpt_n,powi_0, powi_n, powe_0, powe_n,pow_ei, pow_ei_0, pow_ei_n,&
  !     mass_bc_0, mass_bc_n, mass_bc,mom_0, mom_n, mass_0, mass_n,mass_pin, mass_pin_0,&
  !      mass_pin_n,mom_pin,  mom_pin_0,  mom_pin_n,pow_pin,  pow_pin_0,  pow_pin_n,pow_pin_i,&
  !       pow_pin_e,pow_bc,   pow_bc_0,   pow_bc_n,pow_bc_i, pow_bc_i0,  pow_bc_in,&
  !     pow_bc_e, pow_bc_e0,  pow_bc_en,mom_bc,   mom_bc_0,   mom_bc_n,pow_rad,  pow_rad_0,&
  !       pow_rad_n,pow_rad_i,  pow_rad_i0,  pow_rad_in,pow_rad_e,  pow_rad_e0,  pow_rad_en,&
  !     pow_core, pow_core_0, pow_core_n,pow_core_i, pow_core_i0, pow_core_in,&
  !     pow_core_e, pow_core_e0, pow_core_en,mom_core, mom_core_0, mom_core_n,mass_core,&
  !      mass_core_0, mass_core_n,gamma_i, gamma_e, t0_min, t0_grid, tn_grid,gamma_sum, gamma_net,&
  !      rms_tol, a_maxmom,a_maxpow, a_maxmass, qdel_max, udel_max, max_tol,const_ei,&
  !     ikmax, ikmax_f, ik_xpt_0, ik_xpt_n,ikmid, ikmax_0, stepik, itermax, a_seed,&
  !      a_subm, a_flag,a_bc, a_half, a_pin, a_cfdfile, a_perp, a_error, a_resid,a_symm,&
  !      a_2t, cfd_failed, a_qup
  !
  ! save /cfd_osm_com/
  real*8,public :: a_q2,a_q3,a_artvisc4,a_artvisc2,a_temp,a_grid0,a_gridexp,l_pow,l_iz,&
       l_iz_0,l_iz_n,l_mom,l_en,a_dt_f,a_dt_g,mu0,ka0_ev,a_ex,dx_ave,dx0_bc,dxn_bc,&
       x_max,max_res,rms_res,max_err,rms_err,err_tol,q1,q2,q3,q4,mass,mass_i,l_izoffset,&
       a_i,z_i,c_v,c_p,n0_bc,nn_bc,ti0_bc,tin_bc,te0_bc,ten_bc,m0_bc,mn_bc,f0_bc,fn_bc,&
       p_0,p_n,v_0,v_n,x_xpt_0,x_xpt_n,powi_0,powi_n,powe_0,powe_n,pow_ei,pow_ei_0,pow_ei_n,&
       mass_bc_0,mass_bc_n,mass_bc,mom_0,mom_n,mass_0,mass_n,mass_pin,mass_pin_0,&
       mass_pin_n,mom_pin,mom_pin_0,mom_pin_n,pow_pin,pow_pin_0,pow_pin_n,pow_pin_i,pow_pin_e,&
       pow_bc,pow_bc_0,pow_bc_n,pow_bc_i,pow_bc_i0,pow_bc_in,pow_bc_e,pow_bc_e0,&
       pow_bc_en,mom_bc,mom_bc_0,mom_bc_n,pow_rad,pow_rad_0,pow_rad_n,pow_rad_i,pow_rad_i0,&
       pow_rad_in,pow_rad_e,pow_rad_e0,pow_rad_en,pow_core,pow_core_0,pow_core_n,pow_core_i,&
       pow_core_i0,pow_core_in,pow_core_e,pow_core_e0,pow_core_en,mom_core,mom_core_0,&
       mom_core_n,mass_core,mass_core_0,mass_core_n,gamma_i,gamma_e,t0_grid,tn_grid,&
       t0_min,gamma_sum,gamma_net,max_tol,rms_tol,a_maxmom,a_maxpow,a_maxmass,qdel_max,&
       udel_max,const_ei
  real*8,public,allocatable :: unit(:,:),u(:,:),u0(:,:),w(:,:),f(:,:),g(:,:),h(:,:),&
       q(:,:),q0(:,:),q_new(:,:),tau_i(:),tau_e(:),q_b(:,:),b_field(:),l_b(:),q_grad(:),&
       q_uni(:),q_perp(:,:),q_ei(:),q_joule(:),tau_ei(:),lam_mat(:,:),t_mat(:,:,:),tinv_mat(:,:,:),&
       aaaa3(:,:,:),bbbb3(:,:,:),cccc3(:,:,:),ffff3(:,:),uuuu3(:,:),jac_a3(:,:,:),&
       jac_b3(:,:,:),aaaa4(:,:,:),bbbb4(:,:,:),cccc4(:,:,:),ffff4(:,:),uuuu4(:,:),&
       jac_a4(:,:,:),jac_b4(:,:,:),jac_q4(:,:,:),ep2(:),ep4(:),psi_p(:),aaaa(:,:),bbbb(:,:),&
       cccc(:,:),dddd(:,:),ffff(:,:),ffff_new(:,:),rf(:,:),rg(:,:),rh(:,:),r(:,:),&
       re(:,:),dxc(:),xc(:),xf(:),dtc(:),dtf(:),ka(:),ka_i(:),ka_e(:),mu(:),v_j(:),&
       c_s_j(:),nu_max(:),ddd(:),aaa(:),bbb(:),ccc(:),fff(:),uuu(:),v(:),ti_ev(:),te_ev(:),&
       a_rk(:)
  
  
  
  
  
  
  
  
  integer,public :: ikmax,ikmax_f,ik_xpt_0,ik_xpt_n,ikmid,ikmax_0,stepik,itermax,a_seed,&
       a_subm,a_flag,a_bc,a_half,a_pin,a_cfdfile,a_perp,a_error,a_resid,a_symm,a_2t,&
       a_qup
  integer,public,allocatable :: cfd_failed(:)

  public :: allocate_mod_cfd_osm_com,deallocate_mod_cfd_osm_com

contains

  subroutine allocate_mod_cfd_osm_com
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cfd_osm_com','ALLOCATE')

    call allocate_array(unit,4,4,'unit',ierr)
    call allocate_array(u,1,4,0,nnn,'u',ierr)
    call allocate_array(u0,1,4,0,nnn,'u0',ierr)
    call allocate_array(w,1,4,0,nnn,'w',ierr)
    call allocate_array(f,1,4,0,nnn,'f',ierr)
    call allocate_array(g,1,4,0,nnn,'g',ierr)
    call allocate_array(h,1,4,0,nnn,'h',ierr)
    call allocate_array(q,1,4,0,nnn,'q',ierr)
    call allocate_array(q0,1,4,0,nnn,'q0',ierr)
    call allocate_array(q_new,1,4,0,nnn,'q_new',ierr)
    call allocate_array(tau_i,0,'tau_i',nnn,ierr)
    call allocate_array(tau_e,0,'tau_e',nnn,ierr)
    call allocate_array(q_b,1,4,0,nnn,'q_b',ierr)
    call allocate_array(b_field,0,'b_field',nnn,ierr)
    call allocate_array(l_b,0,'l_b',nnn,ierr)
    call allocate_array(q_grad,0,'q_grad',nnn,ierr)
    call allocate_array(q_uni,0,'q_uni',nnn,ierr)
    call allocate_array(q_perp,1,4,0,nnn,'q_perp',ierr)
    call allocate_array(q_ei,0,'q_ei',nnn,ierr)
    call allocate_array(q_joule,0,'q_joule',nnn,ierr)
    call allocate_array(tau_ei,0,'tau_ei',nnn,ierr)
    call allocate_array(lam_mat,1,3,0,nnn,'lam_mat',ierr)
    call allocate_array(t_mat,1,3,1,3,0,nnn,'t_mat',ierr)
    call allocate_array(tinv_mat,1,3,1,3,0,nnn,'tinv_mat',ierr)
    call allocate_array(aaaa3,1,3,1,3,0,nnn,'aaaa3',ierr)
    call allocate_array(bbbb3,1,3,1,3,0,nnn,'bbbb3',ierr)
    call allocate_array(cccc3,1,3,1,3,0,nnn,'cccc3',ierr)
    call allocate_array(ffff3,1,3,0,nnn,'ffff3',ierr)
    call allocate_array(uuuu3,1,3,0,nnn,'uuuu3',ierr)
    call allocate_array(jac_a3,1,3,1,3,0,nnn,'jac_a3',ierr)
    call allocate_array(jac_b3,1,3,1,3,0,nnn,'jac_b3',ierr)
    call allocate_array(aaaa4,1,4,1,4,0,nnn,'aaaa4',ierr)
    call allocate_array(bbbb4,1,4,1,4,0,nnn,'bbbb4',ierr)
    call allocate_array(cccc4,1,4,1,4,0,nnn,'cccc4',ierr)
    call allocate_array(ffff4,1,4,0,nnn,'ffff4',ierr)
    call allocate_array(uuuu4,1,4,0,nnn,'uuuu4',ierr)
    call allocate_array(jac_a4,1,4,1,4,0,nnn,'jac_a4',ierr)
    call allocate_array(jac_b4,1,4,1,4,0,nnn,'jac_b4',ierr)
    call allocate_array(jac_q4,1,4,1,4,0,nnn,'jac_q4',ierr)
    call allocate_array(ep2,0,'ep2',nnn,ierr)
    call allocate_array(ep4,0,'ep4',nnn,ierr)
    call allocate_array(psi_p,0,'psi_p',nnn,ierr)
    call allocate_array(aaaa,3,nnn,'aaaa',ierr)
    call allocate_array(bbbb,3,nnn,'bbbb',ierr)
    call allocate_array(cccc,3,nnn,'cccc',ierr)
    call allocate_array(dddd,3,nnn,'dddd',ierr)
    call allocate_array(ffff,3,nnn,'ffff',ierr)
    call allocate_array(ffff_new,3,nnn,'ffff_new',ierr)
    call allocate_array(rf,1,4,0,nnn,'rf',ierr)
    call allocate_array(rg,1,4,0,nnn,'rg',ierr)
    call allocate_array(rh,1,4,0,nnn,'rh',ierr)
    call allocate_array(r,1,4,0,nnn,'r',ierr)
    call allocate_array(re,1,4,0,nnn,'re',ierr)
    call allocate_array(dxc,0,'dxc',nnn,ierr)
    call allocate_array(xc,0,'xc',nnn,ierr)
    call allocate_array(xf,0,'xf',nnn,ierr)
    call allocate_array(dtc,0,'dtc',nnn,ierr)
    call allocate_array(dtf,0,'dtf',nnn,ierr)
    call allocate_array(ka,0,'ka',nnn,ierr)
    call allocate_array(ka_i,0,'ka_i',nnn,ierr)
    call allocate_array(ka_e,0,'ka_e',nnn,ierr)
    call allocate_array(mu,0,'mu',nnn,ierr)
    call allocate_array(v_j,0,'v_j',nnn,ierr)
    call allocate_array(c_s_j,0,'c_s_j',nnn,ierr)
    call allocate_array(nu_max,0,'nu_max',nnn,ierr)
    call allocate_array(ddd,0,'ddd',nnn,ierr)
    call allocate_array(aaa,0,'aaa',nnn,ierr)
    call allocate_array(bbb,0,'bbb',nnn,ierr)
    call allocate_array(ccc,0,'ccc',nnn,ierr)
    call allocate_array(fff,0,'fff',nnn,ierr)
    call allocate_array(uuu,0,'uuu',nnn,ierr)
    call allocate_array(v,0,'v',nnn,ierr)
    call allocate_array(ti_ev,0,'ti_ev',nnn,ierr)
    call allocate_array(te_ev,0,'te_ev',nnn,ierr)
    call allocate_array(a_rk,5,'a_rk',ierr)
    call allocate_array(cfd_failed,nnn,'cfd_failed',ierr)

  end subroutine allocate_mod_cfd_osm_com


  subroutine deallocate_mod_cfd_osm_com
    implicit none

    call pr_trace('mod_cfd_osm_com','DEALLOCATE')

    if (allocated(unit)) deallocate(unit)
    if (allocated(u)) deallocate(u)
    if (allocated(u0)) deallocate(u0)
    if (allocated(w)) deallocate(w)
    if (allocated(f)) deallocate(f)
    if (allocated(g)) deallocate(g)
    if (allocated(h)) deallocate(h)
    if (allocated(q)) deallocate(q)
    if (allocated(q0)) deallocate(q0)
    if (allocated(q_new)) deallocate(q_new)
    if (allocated(tau_i)) deallocate(tau_i)
    if (allocated(tau_e)) deallocate(tau_e)
    if (allocated(q_b)) deallocate(q_b)
    if (allocated(b_field)) deallocate(b_field)
    if (allocated(l_b)) deallocate(l_b)
    if (allocated(q_grad)) deallocate(q_grad)
    if (allocated(q_uni)) deallocate(q_uni)
    if (allocated(q_perp)) deallocate(q_perp)
    if (allocated(q_ei)) deallocate(q_ei)
    if (allocated(q_joule)) deallocate(q_joule)
    if (allocated(tau_ei)) deallocate(tau_ei)
    if (allocated(lam_mat)) deallocate(lam_mat)
    if (allocated(t_mat)) deallocate(t_mat)
    if (allocated(tinv_mat)) deallocate(tinv_mat)
    if (allocated(aaaa3)) deallocate(aaaa3)
    if (allocated(bbbb3)) deallocate(bbbb3)
    if (allocated(cccc3)) deallocate(cccc3)
    if (allocated(ffff3)) deallocate(ffff3)
    if (allocated(uuuu3)) deallocate(uuuu3)
    if (allocated(jac_a3)) deallocate(jac_a3)
    if (allocated(jac_b3)) deallocate(jac_b3)
    if (allocated(aaaa4)) deallocate(aaaa4)
    if (allocated(bbbb4)) deallocate(bbbb4)
    if (allocated(cccc4)) deallocate(cccc4)
    if (allocated(ffff4)) deallocate(ffff4)
    if (allocated(uuuu4)) deallocate(uuuu4)
    if (allocated(jac_a4)) deallocate(jac_a4)
    if (allocated(jac_b4)) deallocate(jac_b4)
    if (allocated(jac_q4)) deallocate(jac_q4)
    if (allocated(ep2)) deallocate(ep2)
    if (allocated(ep4)) deallocate(ep4)
    if (allocated(psi_p)) deallocate(psi_p)
    if (allocated(aaaa)) deallocate(aaaa)
    if (allocated(bbbb)) deallocate(bbbb)
    if (allocated(cccc)) deallocate(cccc)
    if (allocated(dddd)) deallocate(dddd)
    if (allocated(ffff)) deallocate(ffff)
    if (allocated(ffff_new)) deallocate(ffff_new)
    if (allocated(rf)) deallocate(rf)
    if (allocated(rg)) deallocate(rg)
    if (allocated(rh)) deallocate(rh)
    if (allocated(r)) deallocate(r)
    if (allocated(re)) deallocate(re)
    if (allocated(dxc)) deallocate(dxc)
    if (allocated(xc)) deallocate(xc)
    if (allocated(xf)) deallocate(xf)
    if (allocated(dtc)) deallocate(dtc)
    if (allocated(dtf)) deallocate(dtf)
    if (allocated(ka)) deallocate(ka)
    if (allocated(ka_i)) deallocate(ka_i)
    if (allocated(ka_e)) deallocate(ka_e)
    if (allocated(mu)) deallocate(mu)
    if (allocated(v_j)) deallocate(v_j)
    if (allocated(c_s_j)) deallocate(c_s_j)
    if (allocated(nu_max)) deallocate(nu_max)
    if (allocated(ddd)) deallocate(ddd)
    if (allocated(aaa)) deallocate(aaa)
    if (allocated(bbb)) deallocate(bbb)
    if (allocated(ccc)) deallocate(ccc)
    if (allocated(fff)) deallocate(fff)
    if (allocated(uuu)) deallocate(uuu)
    if (allocated(v)) deallocate(v)
    if (allocated(ti_ev)) deallocate(ti_ev)
    if (allocated(te_ev)) deallocate(te_ev)
    if (allocated(a_rk)) deallocate(a_rk)
    if (allocated(cfd_failed)) deallocate(cfd_failed)

  end subroutine deallocate_mod_cfd_osm_com

end module mod_cfd_osm_com