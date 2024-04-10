subroutine extra_chis()
  use all_variables, only : global_parameters, zones, megazones
  implicit none
  integer*4,allocatable:: bchi(:,:)
  integer*4,allocatable:: bchis(:,:,:)
  real*8, allocatable:: npen_exp(:,:)
  integer*4 bNx,bNz
  integer*4,allocatable:: nbeg(:),nend(:),begins(:,:),ends(:,:)
  integer*4 :: n_mega,size
  integer*4 :: k,i,j
  do n_mega=1,global_parameters%N_megazones
     size=megazones(n_mega)%size
     bNz=0
     do k=1,size
        bNz=bNz+zones(megazones(n_mega)%zone_number(k))%mesh%Nz
     end do
     bNx=zones(megazones(n_mega)%zone_number(1))%mesh%Nx
     allocate(bchis(1:4,1:bNx,1:bNz))
     bchis=0
     allocate(bchi(1:bNx,1:bNz))
     bchi=0
     !concatenation of the chis
     call concat_chi(n_mega,bchi,bNx,bNz)     
     !now the big chi matrix is defined. Let us determine the extra chis from it
     !let us count the mask beginings on each line (cf doc)
     allocate(nbeg(1:bNx),nend(1:bNx),begins(1:bNx,1:bNz),ends(1:bNx,1:bNz))
     nbeg=0
     nend=0
     do i=1,bNx
        do j=2,bNz
           if((bchi(i,j-1).eq.0).and.(bchi(i,j).eq.1)) then
              nbeg(i)=nbeg(i)+1
              begins(i,nbeg(i))=j
           end if
        end do
     end do
     !let us count the mask endings
     do i=1,bNx
        do j=1,bNz-1
           if((bchi(i,j).eq.1).and.(bchi(i,j+1).eq.0)) then
              nend(i)=nend(i)+1
              ends(i,nend(i))=j
           end if
        end do
     end do
     !4 cases 
     do i=1,bNx
        if(sum(bchi(i,1:bNz)).eq.bNz) then !case XXXXXXXX
           bchis(4,i,1:bNz)=1.D0
           bchis(2,i,1:bNz)=1.D0
        else
           if(sum(bchi(i,1:bNz)).ne.0) then
              if(nend(i).gt.nbeg(i)) then !case  XX\__/XX\___
                 call shape1chi(bchis,bNx,bNz,begins,ends,nbeg,i) 
              else
                 if(nend(i).lt.nbeg(i)) then !case __/XX\__/XX
                    call shape2chi(bchis,bNx,bNz,begins,ends,nbeg,nend,i)
                 else
                    if(bchi(i,1).eq.1) then !case XX\__/XX\__/XX
                       call shape3chi(bchis,bNx,bNz,begins,ends,nbeg,nend,i)
                    else !case __/XX\__/XX\__
                       call shape4chi(bchis,bNx,bNz,begins,ends,nend,i)
                    end if
                 end if
              end if
           end if
        end if
     end do
     call unconcat_chi(n_mega,bchis,bNx,bNz)
     deallocate(bchis,bchi)
     deallocate(begins,ends,nbeg,nend)
  end do
end subroutine extra_chis
