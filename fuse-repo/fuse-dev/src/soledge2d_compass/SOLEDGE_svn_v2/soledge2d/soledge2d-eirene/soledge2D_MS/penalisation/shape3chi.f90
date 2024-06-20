subroutine shape3chi(bchis,bNx,bNz,begins,ends,nbeg,nend,ind) ! XX\__/XX\__/XX
  implicit none
  integer*4 bNx,bNz
  integer*4,dimension(1:4,1:bNx,1:bNz):: bchis
  integer*4,dimension(1:bNx)::nbeg,nend
  integer*4,dimension(1:bNx,1:bNz):: begins,ends
  integer*4 ind
  integer*4 j,k
  integer*4 length,l1,l2
  !masks in the middle of zone
  do k=1,nend(ind)-1
     length=ends(ind,k+1)-begins(ind,k)+1
     if(length.ge.3) then
        l1=min(1,length/3)
        l2=length-2*l1
        do j=begins(ind,k),begins(ind,k)+l1-1
           bchis(1,ind,j)=1.D0
        end do
        do j=begins(ind,k)+l1,begins(ind,k)+l1+l2-1
           bchis(4,ind,j)=1.D0
        end do
        do j=begins(ind,k)+l1+l2,begins(ind,k)+2*l1+l2-1
           bchis(3,ind,j)=1.D0
        end do
        do j=begins(ind,k),begins(ind,k)+length-1
           bchis(2,ind,j)=1.D0
        end do
     else
        if(length.eq.2) then
           bchis(1,ind,begins(ind,k))=1.d0
           bchis(3,ind,ends(ind,k+1))=1.D0
           bchis(2,ind,begins(ind,k))=1.d0
           bchis(2,ind,ends(ind,k+1))=1.D0
        end if
     end if
  end do
  !mask at zone end and beginning
  length=bNz-begins(ind,nbeg(ind))+1+ends(ind,1)
  if(length.ge.3) then
     l1=min(1,length/3)
     l2=length-2*l1
     do j=begins(ind,nbeg(ind)),min(begins(ind,nbeg(ind))+l1-1,bNz)
        bchis(1,ind,j)=1.D0
     end do
     do j=bNz+1,begins(ind,nbeg(ind))+l1-1
        bchis(1,ind,j-bNz)=1.D0
     end do
     do j=ends(ind,1)-l1+1,0
        bchis(3,ind,j+bNz)=1.D0
     end do
     do j=max(1,ends(ind,1)-l1+1),ends(ind,1)
        bchis(3,ind,j)=1.D0
     end do
     do j=begins(ind,nbeg(ind))+l1,min(bNz,begins(ind,nbeg(ind))+l1+l2)
        bchis(4,ind,j)=1.D0
     end do
     do j=max(1,begins(ind,nbeg(ind))+l1-bNz),ends(ind,1)-l1
        bchis(4,ind,j)=1.D0
     end do
     do j=begins(ind,nbeg(ind)),bNz
        bchis(2,ind,j)=1.D0
     end do
     do j=1,ends(ind,1)
        bchis(2,ind,j)=1.D0
     end do
  end if
end subroutine shape3chi
