module mod_cgeom
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! slmod begin - new
  !...    had to reorder the target plasma arrays in order to use them in
  !       an equivalence statement in the calcdist routine:
  !
  !     >  crun,pncnt,pshift,actnds,kdmap,rp,zp,idds,kteds,kfeds,kfids,
  !     >  ktids, knds,tagdv,iking,ikoutg,kvds,keds,kti3ls,kpredbar,ktinj,
  ! slmod end
  ! common /cgeom/ nks,nrs,irsep,irwall,rmin,rmax,zmin,zmax,ircent,dr,nxs,nys,dz,nds,&
  !     ndsin,ishot,irtrap,ikt,tslice,r0,z0,rxp,zxp,ikref,cvmf,kbfst,ikto,ikti,rs,zs,kareas,&
  !     ktotas,ktotvs,kpmaxs,kfegs,kfigs,kalphs,kbetas,ksmaxs,kss,kbfs,kes,kvhs,kks,&
  !     kps,kcurvs,kinds,koutds,knorms,ktibs,ktebs,knbs,kvols,kins,ikins,irins,ikouts,&
  !     irouts,kperps,kbacds,kfords,ikds,irds,dds,bts,thetas,sepdist,sepdist2,rspdist,zspdist,&
  !     knds,kteds,ktids,kvds,crun,pncnt,pshift,actnds,kdmap,rp,zp,idds,kfeds,kfids,&
  !     tagdv,iking,ikoutg,keds,kti3ls,kpredbar,ktinj,karea2,kvol2,ktota2,ktotv2,dds2,&
  !     thetas2,thetag,thetat,dthetg,kory,npolyp,korpg,nvertp,rvertp,zvertp,refct,ksmaxs2,&
  !     kpmaxs2,irsep2,irwall2,irtrap2,nrs2,ndsin2,ndsin3,ikt1,ikt2,cosali,cosalo,costet,&
  !     bratio,rhog,hro,hteta,ksb,kss2,kps2,krb,oktebs,oktibs,oknbs,okvhs,okes,nves,rves,&
  !     zves,areap,kpb,kzb,alph,finds,foutds,ikrefsol,ikrefcore,ikrefpp,tempds,nfla,&
  !     distin,distout,tdistin,tdistout,ikmids,n_ik_offsets,ik_offset_data,sol22_halfringlen_opt,&
  !     rcouter,rcinner,zcouter,zcinner,middist,separatrix_dist,target_orth_angle,&
  !     psifl,knes,kpinchs,kpinchs_para,psitarg,ikmidout_sep,ikmidin_sep,asep,asep_eff,&
  !     asep_tor,asep_tor_eff,targfluxdata,wallfluxdata,wallprad,nopriv,delta_psin_core
  
  !
  ! save /cgeom/
  ! common /sreflect/  s_reflect_opt,s_reflect
  
  ! save /sreflect/
  ! common /geo_flags/ inner_targid,outer_targid,xpoint_up,add_virtual_cells_to_grid
  !
  ! save /geo_flags/
  integer,public :: inner_targid,outer_targid,s_reflect_opt,add_virtual_cells_to_grid
  !
  real,public,allocatable :: s_reflect(:,:)
  ! common /sheath/ nsheath_vali,sheath_vali,nsheath_valo,sheath_valo
  !
  ! save /sheath/
  integer,public :: nsheath_vali,nsheath_valo
  !
  real,public,allocatable :: sheath_vali(:,:),sheath_valo(:,:)
  integer,public :: nrs,irsep,irwall,ishot,irtrap,ikt,ircent,nds,ndsin,nxs,nys,ikref,&
       pncnt,pshift,actnds,npolyp,refct,ikti,ikto,irsep2,irwall2,irtrap2,nrs2,ndsin2,&
       ikrefsol,ikrefcore,ikrefpp,ndsin3,ikt1,ikt2,nves,nfla,ikmidout_sep,ikmidin_sep,n_ik_offsets,&
       sol22_halfringlen_opt
  integer,public,allocatable :: nks(:),ikins(:,:),irins(:,:),ikds(:),ikouts(:,:),irouts(:,:),&
       irds(:),kdmap(:),idds(:,:),kory(:,:),nvertp(:),korpg(:,:),tagdv(:,:),iking(:,:),&
       ikoutg(:,:),ikmids(:)
  !
  
  real, target ,public :: rmin,rmax,zmin,zmax,r0,z0,rxp,zxp,tslice,dr,dz,dthetg
  real, target ,public,allocatable :: rs(:,:),zs(:,:),rves(:),zves(:),bts(:,:),knbs(:,:),&
       kbfst(:,:),ktebs(:,:),ktibs(:,:),ksmaxs(:),kss(:,:),kks(:),kpmaxs(:),dds(:),&
       kbfs(:,:),kes(:,:),kvhs(:,:),kvols(:,:),kins(:,:),kcurvs(:,:),kinds(:,:),koutds(:,:),&
       thetas(:),knorms(:,:),kperps(:,:),ktotas(:),kbacds(:,:),kfords(:,:),ktotvs(:),&
       kareas(:,:),kps(:,:),kalphs(:),kfegs(:,:),kfigs(:,:),kbetas(:),cvmf(:,:),ktota2(:),&
       ktotv2(:),karea2(:,:),kvol2(:,:),rp(:),zp(:),kteds(:),ktids(:),knds(:),kvds(:),&
       keds(:),kti3ls(:),kpredbar(:,:,:),kfeds(:),kfids(:),ktinj(:),sepdist(:),sepdist2(:),&
       rspdist(:),zspdist(:),dds2(:),thetas2(:),bratio(:,:),rvertp(:,:),zvertp(:,:),&
       cosali(:,:),cosalo(:,:),costet(:),distin(:,:),distout(:,:),tdistin(:,:),&
       tdistout(:,:),target_orth_angle(:),psifl(:,:),knes(:,:),kpinchs(:,:),kpinchs_para(:,:),&
       psitarg(:,:)
  !
  real,public :: asep,asep_eff,asep_tor,asep_tor_eff,delta_psin_core
  real,public,allocatable :: thetag(:,:),thetat(:),ksb(:,:),rhog(:,:),hro(:,:),hteta(:,:),&
       oktebs(:,:),oktibs(:,:),oknbs(:,:),okvhs(:,:),okes(:,:),kss2(:,:),kps2(:,:),&
       kpmaxs2(:),ksmaxs2(:),krb(:,:),areap(:),finds(:,:),foutds(:,:),kpb(:,:),kzb(:,:),&
       tempds(:),rcouter(:),rcinner(:),zcouter(:),zcinner(:),middist(:,:),separatrix_dist(:,:),&
       targfluxdata(:,:,:),wallfluxdata(:,:,:),wallprad(:,:),ik_offset_data(:,:)
  !
  ! sl begin
  !
  ! slnote - not all of these arrays are necessary - remove asap
  !
  logical,public :: xpoint_up,nopriv
  ! sl end
  real,public,allocatable :: alph(:,:)
  character,public :: crun*40
  

  public :: allocate_mod_cgeom,deallocate_mod_cgeom,allocate_mod_cgeom_input

contains

  subroutine allocate_mod_cgeom
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cgeom','ALLOCATE')

    call allocate_array(s_reflect,maxnrs,2,'s_reflect',ierr)
    call allocate_array(nks,maxnrs,'nks',ierr)
    call allocate_array(ikins,maxnks,maxnrs,'ikins',ierr)
    call allocate_array(irins,maxnks,maxnrs,'irins',ierr)
    call allocate_array(ikds,maxnds,'ikds',ierr)
    call allocate_array(ikouts,maxnks,maxnrs,'ikouts',ierr)
    call allocate_array(irouts,maxnks,maxnrs,'irouts',ierr)
    call allocate_array(irds,maxnds,'irds',ierr)
    call allocate_array(kdmap,maxnds,'kdmap',ierr)
    call allocate_array(idds,maxnrs,2,'idds',ierr)
    call allocate_array(kory,maxnrs,maxnks,'kory',ierr)
    call allocate_array(nvertp,maxnks*maxnrs,'nvertp',ierr)
    call allocate_array(korpg,maxnks,maxnrs,'korpg',ierr)
    call allocate_array(tagdv,maxnks,maxnrs,'tagdv',ierr)
    call allocate_array(iking,maxnks,maxnrs,'iking',ierr)
    call allocate_array(ikoutg,maxnks,maxnrs,'ikoutg',ierr)
    call allocate_array(ikmids,maxnrs,'ikmids',ierr)
    call allocate_array(rs,maxnks,maxnrs,'rs',ierr)
    call allocate_array(zs,maxnks,maxnrs,'zs',ierr)
    call allocate_array(rves,mves,'rves',ierr)
    call allocate_array(zves,mves,'zves',ierr)
    call allocate_array(bts,maxnks,maxnrs,'bts',ierr)
    call allocate_array(knbs,maxnks,maxnrs,'knbs',ierr)
    call allocate_array(kbfst,maxnrs,2,'kbfst',ierr)
    call allocate_array(ktebs,maxnks,maxnrs,'ktebs',ierr)
    call allocate_array(ktibs,maxnks,maxnrs,'ktibs',ierr)
    call allocate_array(ksmaxs,maxnrs,'ksmaxs',ierr)
    call allocate_array(kss,maxnks,maxnrs,'kss',ierr)
    call allocate_array(kks,maxnrs,'kks',ierr)
    call allocate_array(kpmaxs,maxnrs,'kpmaxs',ierr)
    call allocate_array(dds,maxnds,'dds',ierr)
    call allocate_array(kbfs,maxnks,maxnrs,'kbfs',ierr)
    call allocate_array(kes,maxnks,maxnrs,'kes',ierr)
    call allocate_array(kvhs,maxnks,maxnrs,'kvhs',ierr)
    call allocate_array(kvols,maxnks,maxnrs,'kvols',ierr)
    call allocate_array(kins,maxnks,maxnrs,'kins',ierr)
    call allocate_array(kcurvs,maxnks,maxnrs,'kcurvs',ierr)
    call allocate_array(kinds,maxnks,maxnrs,'kinds',ierr)
    call allocate_array(koutds,maxnks,maxnrs,'koutds',ierr)
    call allocate_array(thetas,maxnds,'thetas',ierr)
    call allocate_array(knorms,maxnks,maxnrs,'knorms',ierr)
    call allocate_array(kperps,maxnks,maxnrs,'kperps',ierr)
    call allocate_array(ktotas,maxnrs,'ktotas',ierr)
    call allocate_array(kbacds,maxnks,maxnrs,'kbacds',ierr)
    call allocate_array(kfords,maxnks,maxnrs,'kfords',ierr)
    call allocate_array(ktotvs,maxnrs,'ktotvs',ierr)
    call allocate_array(kareas,maxnks,maxnrs,'kareas',ierr)
    call allocate_array(kps,maxnks,maxnrs,'kps',ierr)
    call allocate_array(kalphs,maxizs,'kalphs',ierr)
    call allocate_array(kfegs,maxnks,maxnrs,'kfegs',ierr)
    call allocate_array(kfigs,maxnks,maxnrs,'kfigs',ierr)
    call allocate_array(kbetas,maxizs,'kbetas',ierr)
    call allocate_array(cvmf,maxnks,maxnrs,'cvmf',ierr)
    call allocate_array(ktota2,maxnrs,'ktota2',ierr)
    call allocate_array(ktotv2,maxnrs,'ktotv2',ierr)
    call allocate_array(karea2,maxnks,maxnrs,'karea2',ierr)
    call allocate_array(kvol2,maxnks,maxnrs,'kvol2',ierr)
    call allocate_array(rp,maxnds,'rp',ierr)
    call allocate_array(zp,maxnds,'zp',ierr)
    call allocate_array(kteds,maxnds,'kteds',ierr)
    call allocate_array(ktids,maxnds,'ktids',ierr)
    call allocate_array(knds,maxnds,'knds',ierr)
    call allocate_array(kvds,maxnds,'kvds',ierr)
    call allocate_array(keds,maxnds,'keds',ierr)
    call allocate_array(kti3ls,maxnds,'kti3ls',ierr)
    call allocate_array(kpredbar,maxnds,3,2,'kpredbar',ierr)
    call allocate_array(kfeds,maxnds,'kfeds',ierr)
    call allocate_array(kfids,maxnds,'kfids',ierr)
    call allocate_array(ktinj,maxnds,'ktinj',ierr)
    call allocate_array(sepdist,maxnds,'sepdist',ierr)
    call allocate_array(sepdist2,maxnds,'sepdist2',ierr)
    call allocate_array(rspdist,maxnds,'rspdist',ierr)
    call allocate_array(zspdist,maxnds,'zspdist',ierr)
    call allocate_array(dds2,maxnds,'dds2',ierr)
    call allocate_array(thetas2,maxnds,'thetas2',ierr)
    call allocate_array(bratio,maxnks,maxnrs,'bratio',ierr)
    call allocate_array(rvertp,5,maxnks*maxnrs,'rvertp',ierr)
    call allocate_array(zvertp,5,maxnks*maxnrs,'zvertp',ierr)
    call allocate_array(cosali,maxnks,maxnrs,'cosali',ierr)
    call allocate_array(cosalo,maxnks,maxnrs,'cosalo',ierr)
    call allocate_array(costet,maxnds,'costet',ierr)
    call allocate_array(distin,maxnks,maxnrs,'distin',ierr)
    call allocate_array(distout,maxnks,maxnrs,'distout',ierr)
    call allocate_array(tdistin,maxnks,maxnrs,'tdistin',ierr)
    call allocate_array(tdistout,maxnks,maxnrs,'tdistout',ierr)
    call allocate_array(target_orth_angle,maxnds,'target_orth_angle',ierr)
    call allocate_array(psifl,maxnks,maxnrs,'psifl',ierr)
    call allocate_array(knes,maxnks,maxnrs,'knes',ierr)
    call allocate_array(kpinchs,maxnks,maxnrs,'kpinchs',ierr)
    call allocate_array(kpinchs_para,maxnks,maxnrs,'kpinchs_para',ierr)
    call allocate_array(psitarg,maxnrs,2,'psitarg',ierr)
    call allocate_array(thetag,maxnks,maxnrs,'thetag',ierr)
    call allocate_array(thetat,maxnds,'thetat',ierr)
    call allocate_array(ksb,0,maxnks,1,maxnrs,'ksb',ierr)
    call allocate_array(rhog,maxnks,maxnrs,'rhog',ierr)
    call allocate_array(hro,maxnks,maxnrs,'hro',ierr)
    call allocate_array(hteta,maxnks,maxnrs,'hteta',ierr)
    call allocate_array(oktebs,maxnks,maxnrs,'oktebs',ierr)
    call allocate_array(oktibs,maxnks,maxnrs,'oktibs',ierr)
    call allocate_array(oknbs,maxnks,maxnrs,'oknbs',ierr)
    call allocate_array(okvhs,maxnks,maxnrs,'okvhs',ierr)
    call allocate_array(okes,maxnks,maxnrs,'okes',ierr)
    call allocate_array(kss2,maxnks,maxnrs,'kss2',ierr)
    call allocate_array(kps2,maxnks,maxnrs,'kps2',ierr)
    call allocate_array(kpmaxs2,maxnrs,'kpmaxs2',ierr)
    call allocate_array(ksmaxs2,maxnrs,'ksmaxs2',ierr)
    call allocate_array(krb,0,maxnks,1,maxnrs,'krb',ierr)
    call allocate_array(areap,maxnks*maxnrs,'areap',ierr)
    call allocate_array(finds,maxnks,maxnrs,'finds',ierr)
    call allocate_array(foutds,maxnks,maxnrs,'foutds',ierr)
    call allocate_array(kpb,0,maxnks,1,maxnrs,'kpb',ierr)
    call allocate_array(kzb,0,maxnks,1,maxnrs,'kzb',ierr)
    call allocate_array(tempds,maxnds,'tempds',ierr)
    call allocate_array(rcouter,maxnrs,'rcouter',ierr)
    call allocate_array(rcinner,maxnrs,'rcinner',ierr)
    call allocate_array(zcouter,maxnrs,'zcouter',ierr)
    call allocate_array(zcinner,maxnrs,'zcinner',ierr)
    call allocate_array(middist,maxnrs,2,'middist',ierr)
    call allocate_array(separatrix_dist,maxnks,maxnrs,'separatrix_dist',ierr)
    call allocate_array(targfluxdata,maxnds+3,4,4,'targfluxdata',ierr)
    call allocate_array(wallfluxdata,maxpts+5,4,4,'wallfluxdata',ierr)
    call allocate_array(wallprad,maxpts+6,3,'wallprad',ierr)
    call allocate_array(alph,maxnks,maxnrs,'alph',ierr)

  end subroutine allocate_mod_cgeom


  subroutine deallocate_mod_cgeom
    implicit none

    call pr_trace('mod_cgeom','DEALLOCATE')

    if (allocated(s_reflect)) deallocate(s_reflect)
    if (allocated(sheath_vali)) deallocate(sheath_vali)
    if (allocated(sheath_valo)) deallocate(sheath_valo)
    if (allocated(nks)) deallocate(nks)
    if (allocated(ikins)) deallocate(ikins)
    if (allocated(irins)) deallocate(irins)
    if (allocated(ikds)) deallocate(ikds)
    if (allocated(ikouts)) deallocate(ikouts)
    if (allocated(irouts)) deallocate(irouts)
    if (allocated(irds)) deallocate(irds)
    if (allocated(kdmap)) deallocate(kdmap)
    if (allocated(idds)) deallocate(idds)
    if (allocated(kory)) deallocate(kory)
    if (allocated(nvertp)) deallocate(nvertp)
    if (allocated(korpg)) deallocate(korpg)
    if (allocated(tagdv)) deallocate(tagdv)
    if (allocated(iking)) deallocate(iking)
    if (allocated(ikoutg)) deallocate(ikoutg)
    if (allocated(ikmids)) deallocate(ikmids)
    if (allocated(rs)) deallocate(rs)
    if (allocated(zs)) deallocate(zs)
    if (allocated(rves)) deallocate(rves)
    if (allocated(zves)) deallocate(zves)
    if (allocated(bts)) deallocate(bts)
    if (allocated(knbs)) deallocate(knbs)
    if (allocated(kbfst)) deallocate(kbfst)
    if (allocated(ktebs)) deallocate(ktebs)
    if (allocated(ktibs)) deallocate(ktibs)
    if (allocated(ksmaxs)) deallocate(ksmaxs)
    if (allocated(kss)) deallocate(kss)
    if (allocated(kks)) deallocate(kks)
    if (allocated(kpmaxs)) deallocate(kpmaxs)
    if (allocated(dds)) deallocate(dds)
    if (allocated(kbfs)) deallocate(kbfs)
    if (allocated(kes)) deallocate(kes)
    if (allocated(kvhs)) deallocate(kvhs)
    if (allocated(kvols)) deallocate(kvols)
    if (allocated(kins)) deallocate(kins)
    if (allocated(kcurvs)) deallocate(kcurvs)
    if (allocated(kinds)) deallocate(kinds)
    if (allocated(koutds)) deallocate(koutds)
    if (allocated(thetas)) deallocate(thetas)
    if (allocated(knorms)) deallocate(knorms)
    if (allocated(kperps)) deallocate(kperps)
    if (allocated(ktotas)) deallocate(ktotas)
    if (allocated(kbacds)) deallocate(kbacds)
    if (allocated(kfords)) deallocate(kfords)
    if (allocated(ktotvs)) deallocate(ktotvs)
    if (allocated(kareas)) deallocate(kareas)
    if (allocated(kps)) deallocate(kps)
    if (allocated(kalphs)) deallocate(kalphs)
    if (allocated(kfegs)) deallocate(kfegs)
    if (allocated(kfigs)) deallocate(kfigs)
    if (allocated(kbetas)) deallocate(kbetas)
    if (allocated(cvmf)) deallocate(cvmf)
    if (allocated(ktota2)) deallocate(ktota2)
    if (allocated(ktotv2)) deallocate(ktotv2)
    if (allocated(karea2)) deallocate(karea2)
    if (allocated(kvol2)) deallocate(kvol2)
    if (allocated(rp)) deallocate(rp)
    if (allocated(zp)) deallocate(zp)
    if (allocated(kteds)) deallocate(kteds)
    if (allocated(ktids)) deallocate(ktids)
    if (allocated(knds)) deallocate(knds)
    if (allocated(kvds)) deallocate(kvds)
    if (allocated(keds)) deallocate(keds)
    if (allocated(kti3ls)) deallocate(kti3ls)
    if (allocated(kpredbar)) deallocate(kpredbar)
    if (allocated(kfeds)) deallocate(kfeds)
    if (allocated(kfids)) deallocate(kfids)
    if (allocated(ktinj)) deallocate(ktinj)
    if (allocated(sepdist)) deallocate(sepdist)
    if (allocated(sepdist2)) deallocate(sepdist2)
    if (allocated(rspdist)) deallocate(rspdist)
    if (allocated(zspdist)) deallocate(zspdist)
    if (allocated(dds2)) deallocate(dds2)
    if (allocated(thetas2)) deallocate(thetas2)
    if (allocated(bratio)) deallocate(bratio)
    if (allocated(rvertp)) deallocate(rvertp)
    if (allocated(zvertp)) deallocate(zvertp)
    if (allocated(cosali)) deallocate(cosali)
    if (allocated(cosalo)) deallocate(cosalo)
    if (allocated(costet)) deallocate(costet)
    if (allocated(distin)) deallocate(distin)
    if (allocated(distout)) deallocate(distout)
    if (allocated(tdistin)) deallocate(tdistin)
    if (allocated(tdistout)) deallocate(tdistout)
    if (allocated(target_orth_angle)) deallocate(target_orth_angle)
    if (allocated(psifl)) deallocate(psifl)
    if (allocated(knes)) deallocate(knes)
    if (allocated(kpinchs)) deallocate(kpinchs)
    if (allocated(kpinchs_para)) deallocate(kpinchs_para)
    if (allocated(psitarg)) deallocate(psitarg)
    if (allocated(thetag)) deallocate(thetag)
    if (allocated(thetat)) deallocate(thetat)
    if (allocated(ksb)) deallocate(ksb)
    if (allocated(rhog)) deallocate(rhog)
    if (allocated(hro)) deallocate(hro)
    if (allocated(hteta)) deallocate(hteta)
    if (allocated(oktebs)) deallocate(oktebs)
    if (allocated(oktibs)) deallocate(oktibs)
    if (allocated(oknbs)) deallocate(oknbs)
    if (allocated(okvhs)) deallocate(okvhs)
    if (allocated(okes)) deallocate(okes)
    if (allocated(kss2)) deallocate(kss2)
    if (allocated(kps2)) deallocate(kps2)
    if (allocated(kpmaxs2)) deallocate(kpmaxs2)
    if (allocated(ksmaxs2)) deallocate(ksmaxs2)
    if (allocated(krb)) deallocate(krb)
    if (allocated(areap)) deallocate(areap)
    if (allocated(finds)) deallocate(finds)
    if (allocated(foutds)) deallocate(foutds)
    if (allocated(kpb)) deallocate(kpb)
    if (allocated(kzb)) deallocate(kzb)
    if (allocated(tempds)) deallocate(tempds)
    if (allocated(rcouter)) deallocate(rcouter)
    if (allocated(rcinner)) deallocate(rcinner)
    if (allocated(zcouter)) deallocate(zcouter)
    if (allocated(zcinner)) deallocate(zcinner)
    if (allocated(middist)) deallocate(middist)
    if (allocated(separatrix_dist)) deallocate(separatrix_dist)
    if (allocated(targfluxdata)) deallocate(targfluxdata)
    if (allocated(wallfluxdata)) deallocate(wallfluxdata)
    if (allocated(wallprad)) deallocate(wallprad)
    if (allocated(ik_offset_data)) deallocate(ik_offset_data)
    if (allocated(alph)) deallocate(alph)

  end subroutine deallocate_mod_cgeom

  subroutine allocate_mod_cgeom_input
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cgeom','ALLOCATE INPUT')

    call allocate_array(sheath_vali,maxnrs,2,'sheath_vali',ierr)
    call allocate_array(sheath_valo,maxnrs,2,'sheath_valo',ierr)
    call allocate_array(ik_offset_data,maxnrs,4,'ik_offset_data',ierr)

    
  end subroutine allocate_mod_cgeom_input
  
end module mod_cgeom
