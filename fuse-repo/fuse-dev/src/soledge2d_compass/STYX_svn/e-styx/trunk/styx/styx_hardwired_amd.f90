subroutine styx_hardwired_amd
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use eirmod_comxs
  use eirmod_ccona
  use styx2eirene
  implicit none

  real(dp) :: TeLo,pls,coua1,coum1,coum2,coum3,coui1,coui2,coui3,ecoui3
  real(dp) :: DeL,ecoua1,coup1,ecoup1,erate
  real(dp) :: eirene_rate_coeff,eirene_energy_rate_coeff
  integer :: itri

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(DeL,pls,Telo,coua1,coum1,coum2,coum3,coui1,coui2,coui3,ecoua1,coup1,ecoup1,ecoui3,erate) 
	do itri=1,Neir_cells
    		if (.not.lgvac(itri,npls)) then
    			DeL=log(DEIN(itri))
    			pls=max(log(1e8_dp),DeL)	
    			TeLo=log(TEIN(itri))

    			! Electron Impact reactions - rate coefficients
    			coua1 = eirene_rate_coeff(1,TeLo,pls,.false.,1,erate)
    			coum1 = eirene_rate_coeff(6,TeLo,pls,.false.,1,erate)
    			coum2 = eirene_rate_coeff(7,TeLo,pls,.false.,1,erate)
    			coum3 = eirene_rate_coeff(8,TeLo,pls,.false.,1,erate)
    			coui1 = eirene_rate_coeff(9,TeLo,pls,.false.,1,erate)
    			coui2 = eirene_rate_coeff(10,TeLo,pls,.false.,1,erate)
        		coui3 = eirene_rate_coeff(11,TeLo,pls,.false.,1,erate)
    		
    			coua1 = coua1 + DeL
   			coum1 = coum1 + DeL
    			coum2 = coum2 + DeL
    			coum3 = coum3 + DeL
			coui1 = coui1 + DeL
    			coui2 = coui2 + DeL
    			coui3 = coui3 + DeL

        		TABDS1(1,itri)  = fa1*exp(max(-100._dp,coua1))
    			TABDS1(2,itri)  = fm1*exp(max(-100._dp,coum1))
   			TABDS1(3,itri)  = fm2*exp(max(-100._dp,coum2))
    			TABDS1(4,itri)  = fm3*exp(max(-100._dp,coum3))
			TABDS1(5,itri)  = fi1*exp(max(-100._dp,coui1))
   			TABDS1(6,itri)  = fi2*exp(max(-100._dp,coui2))
    			TABDS1(7,itri)  = fi3*exp(max(-100._dp,coui3))

			! Electron Impact reactions - energy losses

   			ecoua1 = eirene_energy_rate_coeff(4,TeLo,pls,.false.,1)
   			ecoua1 = ecoua1 + Del
     			EELDS1(1,itri) = - fa1*exp(max(-100._dp,ecoua1))/(TABDS1(1,itri)+eps60)
  		
			! Values for processes 6 and 7 already set in initialisation phase
  			! They do not depend on plasma parameters
    			!EELDS1(4,itri)=-eionh2
    			! idem for processes 9 and 10
                        ecoui3 = eirene_energy_rate_coeff(12,TeLo,pls,.false.,1)
                        ecoui3 = ecoui3 + Del

    			EELDS1(7,itri)= -fi3*exp(max(-100._dp,ecoui3))/(TABDS1(7,itri)+eps60)
     			EHVDS1(7,itri)= - EELDS1(7,itri)

    			! Recombination (if active)
   			coup1 = eirene_rate_coeff(3,TeLo,pls,.true.,1,erate)
   			TABRC1(1,itri)=coup1 *fp1*DEIN(itri)

  			ecoup1 = eirene_energy_rate_coeff(5,TeLo,pls,.false.,1)
  			ecoup1 = ecoup1 + Del
   			EELRC1(1,itri)= fp1*(-exp(max(-100._dp,ecoup1))+delpot(5)*TABRC1(1,itri))
                        

    			! Charge exchange reactions
   			pls=log(TIIN(1,itri))+addtls
    			coua1 = eirene_rate_coeff(2,pls,0._dp,.true.,0,erate)
    			TABCX3(1,itri,1) = coua1*DIIN(1,itri)*fa2
    			EPLCX3(1,itri,1) = 1.5_dp*TIIN(1,itri)+EDRIFT(1,itri)

			
     		endif	

    	enddo
!$OMP END PARALLEL DO

end subroutine styx_hardwired_amd

