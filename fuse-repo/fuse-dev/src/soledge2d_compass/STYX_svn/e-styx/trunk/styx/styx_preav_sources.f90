subroutine styx_preav_sources
  USE EIRMOD_PRECISION
  USE EIRMOD_PARMMOD
  USE EIRMOD_COMUSR
  USE EIRMOD_COMPRT, ONLY: IUNOUT
  USE EIRMOD_CCONA
  USE EIRMOD_COMXS
  use all_variables, only : global_parameters
  use eirmod_precision
  use styx2eirene

  implicit none

  real(dp), allocatable :: nu_reac(:)
  real(dp), allocatable :: Sn(:,:)
  real(dp), allocatable :: Sm(:,:)
  real(dp), allocatable :: SE(:,:)
  real(dp), allocatable :: Srad(:,:)

  real(dp), allocatable :: ESIGEI_styx(:,:,:)  

  integer :: IATM,IMOL,IION,ITYP
  integer :: IREI,NRC,KK
  integer :: IRCX,iacx,imcx
  integer :: ipls,ipls1,irc,irrc,irel,i
  integer :: eirene_idez

  real(dp) :: tst2,tend2,omp_get_wtime

  Sn_tot=0._dp
  Sm_tot=0._dp
  SE_tot=0._dp

  Sn_at=0._dp
  Sm_at=0._dp
  SE_at=0._dp

!  SE_i_io=0._dp
!  SE_i_cx=0._dp

  Sn_mol=0._dp
  Sm_mol=0._dp
  SE_mol=0._dp

  Sn_tion=0._dp
  Sm_tion=0._dp
  SE_tion=0._dp
 
  Sn_pls=0._dp
  Sm_pls=0._dp
  SE_pls=0._dp
  
  Srad_at=0._dp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  allocate(ESIGEI_styx(NRDS,4:5,Neir_cells))
  allocate(nu_reac(1:Neir_cells))
  allocate(Sn(Neir_cells,0:nplsi),Sm(Neir_cells,1:nplsi))
  allocate(SE(Neir_cells,0:nplsi),Srad(Neir_cells,1:natmi))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!        ATOMS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     ATOMS      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    ATOMS      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IREI=0

  do iatm=1,NATMI
    if (NRCA(IATM).EQ.0.AND.NCHARA(IATM).LE.2) then
!
!  DEFAULT H,D,T OR HE ELEC. IMP. IONIZATION MODEL
!
      do IPLS=1,NPLSI
        if (NCHARP(IPLS).EQ.NCHARA(IATM).AND.NMASSP(IPLS).EQ.NMASSA(IATM).AND. NCHRGP(IPLS).EQ.1) then

          Sn=0._dp
          Sm=0._dp
          SE=0._dp

          IREI=IREI+1
          ESIGEI_styx(IREI,4,1:Neir_cells)=EPLDS(IREI,1)*E0_at(iatm,1:Neir_cells)+EPLDS(IREI,2)*EHVDS1(IREI,1:Neir_cells)*elcha
          ESIGEI_styx(IREI,5,1:Neir_cells)=EELDS1(IREI,1:Neir_cells)*elcha
  						
  	  nu_reac(1:Neir_cells)=atom_density(iatm,1:Neir_cells)*TABDS1(IREI,1:Neir_cells)

	  Sn(1:Neir_cells,0) = nu_reac(1:Neir_cells)*pelds(IREI)
          Sn(1:Neir_cells,ipls) = Sn(1:Neir_cells,0)
          Sm(1:Neir_cells,ipls) = RMASSA(iatm)*amuakg*v0par_at(iatm,1:Neir_cells)*nu_reac(1:Neir_cells)*pplds(IREI,1)
          SE(1:Neir_cells,0) = ESIGEI_styx(IREI,5,1:Neir_cells)*nu_reac(1:Neir_cells)
          SE(1:Neir_cells,ipls) = ESIGEI_styx(IREI,4,1:Neir_cells)*nu_reac(1:Neir_cells)
               
  	  Sn_at = Sn_at + Sn
          Sm_at = Sm_at + Sm						
	  SE_at = SE_at + SE 
          !SE_i_at = SE_i_at + SE_i
          !SE_i_io = SE_i                                     
        endif
      enddo
!
!  NON DEFAULT ELEC. IMP. COLLISION MODEL,
!
    elseif (NRCA(iatm)>0) then
      do NRC=1,NRCA(IATM)
        KK=IREACA(IATM,NRC)
          if (ISWR(KK) == 1) then

            Sn=0._dp
            Sm=0._dp
            SE=0._dp
            Srad=0._dp

            IREI=IREI+1

            ! index of ion created
            ITYP=EIRENE_IDEZ(ISCD1A(iatm,NRC),1,3)
            if (ityp /= 4) then
              write(*,*) 'problem with ISCD1A ....'
              stop
            endif
            ipls=EIRENE_IDEZ(ISCD1A(iatm,NRC),3,3)
              
            ! the energy loss per ionization for electrons in EELDS1 (eV/event)
            ! the energy gain for ions in EHVDS1s (eV/event)
            
            ESIGEI_styx(IREI,4,1:Neir_cells)=EPLDS(IREI,1)*E0_at(iatm,1:Neir_cells)+EPLDS(IREI,2)*EHVDS1(IREI,1:Neir_cells)*elcha
            ESIGEI_styx(IREI,5,1:Neir_cells)=EELDS1(IREI,1:Neir_cells)*elcha

            nu_reac(1:Neir_cells)=atom_density(iatm,1:Neir_cells)*TABDS1(IREI,1:Neir_cells)
            Sn(1:Neir_cells,0) = nu_reac(1:Neir_cells)*pelds(IREI)
            Sn(1:Neir_cells,ipls) = Sn(1:Neir_cells,0)
            Sm(1:Neir_cells,ipls) = RMASSA(iatm)*amuakg*v0par_at(iatm,1:Neir_cells)*nu_reac(1:Neir_cells)*pplds(IREI,ipls)
            SE(1:Neir_cells,0) = ESIGEI_styx(IREI,5,1:Neir_cells)*nu_reac(1:Neir_cells)
            SE(1:Neir_cells,ipls) = ESIGEI_styx(IREI,4,1:Neir_cells)*nu_reac(1:Neir_cells)

            ! radiated power [plt], resolved per atomic species
            if (iatm == 1) then
              Srad(1:Neir_cells,iatm)=SE(1:Neir_cells,0) - nu_reac(1:Neir_cells)*13.6_dp*elcha
            elseif (iatm < natmi-ispc_add+1) then
              Srad(1:Neir_cells,iatm)=SE(1:Neir_cells,0) - nu_reac(1:Neir_cells)* &
                     global_parameters%element_list(iatm)%amdatas(1)%ionization_potential*elcha 
            endif
            Sn_at = Sn_at + Sn
            Sm_at = Sm_at + Sm
            SE_at = SE_at + SE
            Srad_at = Srad_at + Srad
    !SE_i_io = SE_i
          endif
        enddo
      endif  
    enddo

! then charge-exchange

    IRCX=0
    do iatm=1,NATMI
      iacx=0
      if (NRCA(iatm) == 0) then

! default CXmodel : not implemented
	
      elseif (NRCA(iatm) > 0) then

! non default model

      do NRC=1,NRCA(IATM)
        KK=IREACA(IATM,NRC)
        if (ISWR(KK) == 3) then

          Sn=0._dp
          Sm=0._dp
          SE=0._dp

          IRCX=IRCX+1
          
          iacx=iacx+1
          if (lgacx(iatm,iacx,0) /= ircx) then
            write(*,*) 'issue in cx reaction indexing, from styx_preav_sources'
          endif
          ! index of ion consumed in the reaction          
          ipls1=lgacx(iatm,iacx,1)
    
          ! index of ion created (try first secondary, then second if not bulk ion)
           ! this assumes reaction between one neutral and one ion
          ITYP=EIRENE_IDEZ(ISCD1A(iatm,NRC),1,3)
          if (ityp /= 4) then
            ITYP=EIRENE_IDEZ(ISCD2A(iatm,NRC),1,3)
            if (ityp /= 4) then
              write(*,*) 'problem with ISCD1A ....'
              stop
            else
              ipls=EIRENE_IDEZ(ISCD2A(iatm,NRC),3,3)
            endif
          else
            ipls=EIRENE_IDEZ(ISCD1A(iatm,NRC),3,3)
          endif
          ! check that TABCR3 is calculated with DIIN !

          if (ipls == ipls1) then
            Sm(1:Neir_cells,ipls) =  atom_density(iatm,1:Neir_cells)*TABCX3(IRCX,1:Neir_cells,1)*&
                                                     (RMASSA(iatm)*amuakg*v0par_at(iatm,1:Neir_cells)-&
                                                    RMASSP(ipls)*amuakg*vpar_tri(1:Neir_cells,ipls))
            SE(1:Neir_cells,ipls) =  atom_density(iatm,1:Neir_cells)*TABCX3(IRCX,1:Neir_cells,1)*&
                                                    (E0_at(iatm,1:Neir_cells)- EPLCX3(IRCX,1:Neir_cells,1)*elcha)
          else
            Sn(1:Neir_cells,ipls1)= - atom_density(iatm,1:Neir_cells)*TABCX3(IRCX,1:Neir_cells,1)
            Sn(1:Neir_cells,ipls) = -Sn(1:Neir_cells,ipls1)
        
            Sm(1:Neir_cells,ipls1) = Sn(1:Neir_cells,ipls1)*RMASSP(ipls1)*amuakg*vpar_tri(1:Neir_cells,ipls1)
            Sm(1:Neir_cells,ipls)  = Sn(1:Neir_cells,ipls)*RMASSA(iatm)*amuakg*v0par_at(iatm,1:Neir_cells)

            SE(1:Neir_cells,ipls1) = Sn(1:Neir_cells,ipls1)*EPLCX3(IRCX,1:Neir_cells,1)*elcha
            SE(1:Neir_cells,ipls)  = Sn(1:Neir_cells,ipls)*E0_at(iatm,1:Neir_cells)
          endif
          Sn_at = Sn_at + Sn
          Sm_at = Sm_at + Sm 
          SE_at = SE_at +SE
        endif
      enddo
    endif
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!        MOLECULES       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     MOLECULES      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    MOLECULES      !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do imol=1,NMOLI
    if (NRCM(IMOL).EQ.0.AND.NCHARM(IMOL).EQ.2) then
    !  APPLY THE DEFAULT MODEL FOR H2 DISSOCIATION AND IONIZATION
    !  TO ALL HYDROGENIC MOLECULES again 3 reactions
      do i=1,3

        Sn=0._dp
        Sm=0._dp
        SE=0._dp

        IREI=IREI+1

        ! index of ion created (try first secondary, then second if not bulk ion)
        ITYP=EIRENE_IDEZ(ISCD1M(imol,NRC),1,3)
        if (ityp /= 4) then
          ITYP=EIRENE_IDEZ(ISCD2M(imol,NRC),1,3)
          if (ityp /= 4) then
            !write(*,*) 'no bulk ion formed for reaction ',NRC,' imol = ',imol
            ipls=0
          elseif (ityp /=4) then
            write(*,*) 'problem with ISCD1M, reaction ',NRC,' imol = ',imol
            stop
          else
            ipls=EIRENE_IDEZ(ISCD2M(imol,NRC),3,3) 
          endif
        else
          ipls=EIRENE_IDEZ(ISCD1M(imol,NRC),3,3)
        endif

        ESIGEI_styx(IREI,4,1:Neir_cells)=EPLDS(IREI,1)*E0_mol(imol,1:Neir_cells)+EPLDS(IREI,2)*EHVDS1(IREI,1:Neir_cells)*elcha
        ESIGEI_styx(IREI,5,1:Neir_cells)=EELDS1(IREI,1:Neir_cells)*elcha

        nu_reac(1:Neir_cells)=mol_density(imol,1:Neir_cells)*TABDS1(IREI,1:Neir_cells)
        Sn(1:Neir_cells,0) = nu_reac(1:Neir_cells)*pelds(IREI)
        SE(1:Neir_cells,0) = ESIGEI_styx(IREI,5,1:Neir_cells)*nu_reac(1:Neir_cells)

        if (ipls /= 0) then
          Sn(1:Neir_cells,ipls) = nu_reac(1:Neir_cells)*pplds(IREI,ipls)
          Sm(1:Neir_cells,ipls) = RMASSM(imol)*amuakg*v0par_mol(imol,1:Neir_cells)*nu_reac(1:Neir_cells)*pplds(IREI,ipls)
          SE(1:Neir_cells,ipls) = ESIGEI_styx(IREI,4,1:Neir_cells)*nu_reac(1:Neir_cells)
        endif
 
        Sn_mol = Sn_mol + Sn
        Sm_mol = Sm_mol + Sm
        SE_mol = SE_mol + SE 
      enddo
!
!  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
!
    elseif (NRCM(IMOL).GT.0) then
      do NRC=1,NRCM(IMOL)
        KK=IREACM(IMOL,NRC)
        if (ISWR(KK) == 1) then
 
          Sn=0._dp
          Sm=0._dp
          SE=0._dp

          IREI=IREI+1
! index of ion created (try first secondary, then second if not bulk ion)
          ITYP=EIRENE_IDEZ(ISCD1M(imol,NRC),1,3)
          if (ityp /= 4) then
            ITYP=EIRENE_IDEZ(ISCD2M(imol,NRC),1,3)
            if (ityp /= 4) then
             ! write(*,*) 'no bulk ion formed for reaction ',NRC,' imol = ',imol
              ipls=0
            elseif (ityp /=4) then
              write(*,*) 'problem with ISCD1M, reaction ',NRC,' imol = ',imol
              stop
            else
              ipls=EIRENE_IDEZ(ISCD2M(imol,NRC),3,3) 
            endif
          else
            ipls=EIRENE_IDEZ(ISCD1M(imol,NRC),3,3)
          endif

          ESIGEI_styx(IREI,4,1:Neir_cells)=EPLDS(IREI,1)*E0_mol(imol,1:Neir_cells)+EPLDS(IREI,2)*EHVDS1(IREI,1:Neir_cells)*elcha
          ESIGEI_styx(IREI,5,1:Neir_cells)=EELDS1(IREI,1:Neir_cells)*elcha

          nu_reac(1:Neir_cells)= mol_density(imol,1:Neir_cells)*TABDS1(IREI,1:Neir_cells)

          Sn(1:Neir_cells,0)= nu_reac(1:Neir_cells)*pelds(IREI)
          SE(1:Neir_cells,0) = ESIGEI_styx(IREI,5,1:Neir_cells)*nu_reac(1:Neir_cells)

          if (ipls /= 0) then
            Sn(1:Neir_cells,ipls)= nu_reac(1:Neir_cells)*pplds(IREI,ipls)
  	    Sm(1:Neir_cells,ipls) = RMASSM(imol)*amuakg*v0par_mol(imol,1:Neir_cells)*nu_reac(1:Neir_cells)*pplds(IREI,ipls)
	    SE(1:Neir_cells,ipls) = ESIGEI_styx(IREI,4,1:Neir_cells)*nu_reac(1:Neir_cells)
          endif
 
          Sn_mol = Sn_mol + Sn
          Sm_mol = Sm_mol + Sm
          SE_mol = SE_mol + SE 
        endif
      enddo
    endif
  enddo

!   now charge exchange
  do IMOL=1,NMOLI
    imcx=0
    if (NRCM(IMOL).EQ.0) then
!  default model : no cx

!  non default cx model :
    elseif (NRCM(IMOL).GT.0) then
      do NRC=1,NRCM(IMOL)
        KK=IREACM(IMOL,NRC)
        if (ISWR(KK) == 3) then
          
          Sn=0._dp
          Sm=0._dp
          SE=0._dp

          IRCX=IRCX+1
          imcx=imcx+1
          if (lgmcx(imol,imcx,0) /= ircx) then
            write(*,*) 'issue in cx reaction indexing, from styx_preav_sources: molecules'
          endif
          ! index of ion consumed in the reaction          
          ipls1=lgmcx(imol,imcx,1)

          ! index of ion reacting
          ITYP=EIRENE_IDEZ(IBULKM(imol,NRC),1,3)
          if (ityp /= 4) then
            write(*,*) 'problem with IBULKM, cx ....'
            stop
          endif
          ipls=EIRENE_IDEZ(IBULKM(imol,NRC),3,3)

! not finished, one the ions may be a test ion: to be included propely.


! check here that D+ ion is not a secondary (only losses from the impinging D+ are included here)
          Sm(1:Neir_cells,ipls) = -1._dp *mol_density(imol,1:Neir_cells)*TABCX3(IRCX,1:Neir_cells,1)* &
                                     RMASSP(ipls)*amuakg*vpar_tri(1:Neir_cells,ipls)
          SE(1:Neir_cells,ipls) = -1._dp * mol_density(imol,1:Neir_cells)*TABCX3(IRCX,1:Neir_cells,1)* &
                                     EPLCX3(IRCX,1:Neir_cells,1)*elcha
          Sm_mol = Sm_mol + Sm 
          SE_mol = SE_mol + SE
        endif
      enddo
    endif
   enddo

! now elastic collisions

   IREL=0

   do imol=1,nmoli        	
!
!  DEFAULT EL MODEL: NOT AVAILABLE
!
     IF (NRCM(IMOL).EQ.0) THEN
          		!NMELI(IMOL)=0
!
!  NON DEFAULT EL MODEL
!
     ELSEIF (NRCM(IMOL).GT.0) THEN
       do NRC=1,NRCM(IMOL)
         KK=IREACM(IMOL,NRC)
           if (ISWR(KK)== 5) then

             Sn=0._dp
             Sm=0._dp
             SE=0._dp

             IREL=IREL+1
             ! index of ion reacting
             ITYP=EIRENE_IDEZ(IBULKM(imol,NRC),1,3)
             if (ityp /= 4) then
               write(*,*) 'problem with IBULKM, el ....'
               stop
             endif
             ipls=EIRENE_IDEZ(IBULKM(imol,NRC),3,3)
          
!  BGK SELF AND CROSS COLLISIONS? nothing to do here if neutral/neutral collision
!  reactions keep the same index as in EIRENE
             if (IBGKM(IMOL,NRC) == 0) then
               Sm(1:Neir_cells,ipls) = mol_density(imol,1:Neir_cells)*TABEL3(IREL,1:Neir_cells,1)* &
                                 (RMASSM(iatm)*amuakg*v0par_mol(imol,1:Neir_cells)-RMASSP(ipls)*amuakg*vpar_tri(1:Neir_cells,ipls))
               SE(1:Neir_cells,ipls)= -1._dp*mol_density(imol,1:Neir_cells)*eiv_i(ipls)%dens(1:Neir_cells)*EPLEL3(IREL,1:Neir_cells,1)
               
               Sm_mol=Sm_mol+Sm
               SE_mol=SE_mol+SE
             endif
           endif
         enddo
      endif
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!        TEST IONS       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      TEST IONS      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    TEST IONS      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do iion=1,NIONI
       if (NRCI(IION).EQ.0.AND.NCHARI(IION).EQ.2.AND.NPRT(NSPAM+IION).GT.1) then
!  APPLY THE DEFAULT MODEL FOR H2+ DISSOCIATION
         do i=1,3

           Sn=0._dp
           Sm=0._dp
           SE=0._dp

           IREI=IREI+1

           ! index of ion created (try first secondary, then second if not bulk ion)
           ITYP=EIRENE_IDEZ(ISCD1I(iion,NRC),1,3)
           if (ityp /= 4) then
             ITYP=EIRENE_IDEZ(ISCD2I(iion,NRC),1,3)
             if (ityp /= 4) then
               !write(*,*) 'no bulk ion formed for reaction ',NRC,' iion = ',iion
               ipls=0
             elseif (ityp /=4) then
               write(*,*) 'problem with ISCD1M, reaction ',NRC,' iion = ',iion
               stop
             else
               ipls=EIRENE_IDEZ(ISCD2I(iion,NRC),3,3) 
             endif
           else
             ipls=EIRENE_IDEZ(ISCD1M(iion,NRC),3,3)
           endif

           ESIGEI_styx(IREI,4,1:Neir_cells)=EPLDS(IREI,1)*E0_tion(iion,1:Neir_cells)+EPLDS(IREI,2)*EHVDS1(IREI,1:Neir_cells)*elcha
           ESIGEI_styx(IREI,5,1:Neir_cells)=EELDS1(IREI,1:Neir_cells)*elcha

           nu_reac(1:Neir_cells) = tion_density(iion,1:Neir_cells)*TABDS1(IREI,1:Neir_cells)

           Sn(1:Neir_cells,0) = nu_reac(1:Neir_cells)*pelds(IREI)
           SE(1:Neir_cells,0) = ESIGEI_styx(IREI,5,1:Neir_cells)*nu_reac(1:Neir_cells)

           if (ipls /=0) then
             Sn(1:Neir_cells,ipls)=  nu_reac(1:Neir_cells)*pplds(IREI,ipls)
             Sm(1:Neir_cells,ipls) = RMASSI(iion)*amuakg*v0par_tion(iion,1:Neir_cells)*nu_reac(1:Neir_cells)*pplds(IREI,ipls)           
             SE(1:Neir_cells,ipls) = ESIGEI_styx(IREI,4,1:Neir_cells)*nu_reac(1:Neir_cells)
           endif
 
           Sn_tion = Sn_tion + Sn
           Sm_tion = Sm_tion + Sm
           SE_tion = SE_tion + SE 
         enddo		
       elseif (NRCI(IION).GT.0) then
         do NRC=1,NRCI(IION)
           KK=IREACI(IION,NRC)
           if (ISWR(KK) == 1) then

             Sn=0._dp
             Sm=0._dp
             SE=0._dp

             IREI=IREI+1
             ! index of ion created (try first secondary, then second if not bulk ion)
             ITYP=EIRENE_IDEZ(ISCD1I(iion,NRC),1,3)
             if (ityp /= 4) then
               ITYP=EIRENE_IDEZ(ISCD2I(iion,NRC),1,3)
               if (ityp /= 4) then
                 !write(*,*) 'no bulk ion formed for reaction ',NRC,' iion = ',iion
                 ipls=0
               elseif (ityp /=4) then
                 write(*,*) 'problem with ISCD1M, reaction ',NRC,' iion = ',iion
                 stop
               else
                 ipls=EIRENE_IDEZ(ISCD2I(iion,NRC),3,3) 
               endif
             else
               ipls=EIRENE_IDEZ(ISCD1M(iion,NRC),3,3)
             endif


             ESIGEI_styx(IREI,4,1:Neir_cells)=EPLDS(IREI,1)*E0_tion(iion,1:Neir_cells)+EPLDS(IREI,2)*EHVDS1(IREI,1:Neir_cells)*elcha
             ESIGEI_styx(IREI,5,1:Neir_cells)=EELDS1(IREI,1:Neir_cells)*elcha

             nu_reac(1:Neir_cells) = tion_density(iion,1:Neir_cells)*TABDS1(IREI,1:Neir_cells)

             Sn(1:Neir_cells,0) = nu_reac(1:Neir_cells)*pelds(IREI)
             SE(1:Neir_cells,0) = ESIGEI_styx(IREI,5,1:Neir_cells)*nu_reac(1:Neir_cells)

             if (ipls /=0) then
               Sn(1:Neir_cells,ipls)=  nu_reac(1:Neir_cells)*pplds(IREI,ipls)
               Sm(1:Neir_cells,ipls) = RMASSI(iion)*amuakg*v0par_tion(iion,1:Neir_cells)*nu_reac(1:Neir_cells)*pplds(IREI,ipls)
               SE(1:Neir_cells,ipls) = ESIGEI_styx(IREI,4,1:Neir_cells)*nu_reac(1:Neir_cells)
             endif
 
             Sn_tion = Sn_tion + Sn
	     Sm_tion = Sm_tion + Sm
             SE_tion = SE_tion + SE 

           endif
         enddo	
       endif
     enddo

! Now recombination

! conversion from eV.s-1 to W
     EELRC1=EELRC1*elcha

! as many singly charged ions as species
     do i=1,global_parameters%n_species
       ipls=singly_charged_ions(i)
       ! corresponding atomic species in EIRENE
       iatm=i
       IRRC=0
       do irc=1,NRCP(ipls)
  
         Sn=0._dp
         Sm=0._dp
         SE=0._dp
         Srad=0._dp

         IRRC=IRRC+1

         Sn(1:Neir_cells,0)=-1._dp*eiv_i(ipls)%dens(1:Neir_cells)*TABRC1(IRRC,1:Neir_cells)
         SE(1:Neir_cells,0)= eiv_i(ipls)%dens(1:Neir_cells)*EELRC1(IRRC,1:Neir_cells)

         Sn(1:Neir_cells,ipls)=Sn(1:Neir_cells,0)
         Sm(1:Neir_cells,ipls)=(RMASSP(ipls)*amuakg)*vpar_tri(1:Neir_cells,ipls)*Sn(1:Neir_cells,ipls)
         SE(1:Neir_cells,ipls)=-1._dp*(1.5_dp*eiv_i(ipls)%T(1:Neir_cells)+Edrift(ipls,1:Neir_cells))*elcha*Sn(1:Neir_cells,ipls)
 
         ! radiated power [prb], resolved per atomic species; Sn<0 so gain of electrons
         if (iatm == 1) then
           Srad(1:Neir_cells,iatm)=SE(1:Neir_cells,0) - Sn(1:Neir_cells,0)*13.6_dp*elcha
         else
           Srad(1:Neir_cells,iatm)=SE(1:Neir_cells,0) - Sn(1:Neir_cells,0)* &
                 global_parameters%element_list(i)%amdatas(1)%ionization_potential*elcha
         endif

         Sn_pls = Sn_pls + Sn
         Sm_pls = Sm_pls + Sm
         SE_pls = SE_pls + SE
         Srad_at = Srad_at + Srad

       enddo
     enddo

     deallocate(nu_reac,ESIGEI_styx)
     deallocate(Sn,Sm,SE,Srad)

! implement tweaks on radiated power (working with negative numbers here !)

     if (tweak_chemistry) then
       write(*,*) 'warning, radiation power tricks to be re-implemented !!!'
!       if (.not.allocated(Srad_h)) allocate(Srad_h(Neir_cells))
!         Srad_h=SE_at(:,1)+13.6_dp*elcha*Sn_at(:,1)
!         SE_at(:,1)=-13.6_dp*elcha*Sn_at(:,1)+faktE(1)*Srad_h
     endif

! remove mass from momentum sources

     do ipls=1,nplsi
       Sm_at(:,ipls)  = Sm_at(:,ipls)/(RMASSP(ipls)*amuakg)
       Sm_mol(:,ipls) = Sm_mol(:,ipls)/(RMASSP(ipls)*amuakg)
       Sm_tion(:,ipls)= Sm_tion(:,ipls)/(RMASSP(ipls)*amuakg)
       Sm_pls(:,ipls) = Sm_pls(:,ipls)/(RMASSP(ipls)*amuakg)
     enddo
! calculate total sources (sum over species)

     Sn_tot=Sn_at+Sn_mol+Sn_tion+Sn_pls
     Sm_tot=Sm_at+Sm_mol+Sm_tion+Sm_pls
     SE_tot=SE_at+SE_mol+SE_tion+SE_pls

! calculate volume integrals to check balances (only after calling EIRENE)

     if (.not.short_cycle) then

       tst2 = omp_get_wtime()
 
       do ipls=0,nplsi
         call styx_total_sources(Sn_at(:,ipls),Sn_at_intg(ipls))
         call styx_total_sources(Sn_mol(:,ipls),Sn_mol_intg(ipls))
         call styx_total_sources(Sn_tion(:,ipls),Sn_tion_intg(ipls))
         call styx_total_sources(Sn_pls(:,ipls),Sn_pls_intg(ipls))
         call styx_total_sources(SE_at(:,ipls),SE_at_intg(ipls))
         call styx_total_sources(SE_mol(:,ipls),SE_mol_intg(ipls))
         call styx_total_sources(SE_tion(:,ipls),SE_tion_intg(ipls))
         call styx_total_sources(SE_pls(:,ipls),SE_pls_intg(ipls))
       enddo

       do ipls=1,nplsi
         call styx_total_sources(Sm_at(:,ipls),Sm_at_intg(ipls))
         call styx_total_sources(Sm_mol(:,ipls),Sm_mol_intg(ipls))
         call styx_total_sources(Sm_tion(:,ipls),Sm_tion_intg(ipls))
         call styx_total_sources(Sm_pls(:,ipls),Sm_pls_intg(ipls))
       enddo
    

       !call styx_total_sources(SE_i_io,SE_i_io_intg)
       !call styx_total_sources(SE_i_cx,SE_i_cx_intg)	

       tend2 = omp_get_wtime()
       write(*,*) ' time for integrations in styx_preav_sources (ms) = ',(tend2-tst2)*1e3_dp
    endif


end subroutine

