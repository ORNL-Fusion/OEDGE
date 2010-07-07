      subroutine xyz_to_rst (icell, x1, y1, z1, x2, y2, z2, x3, y3, z3, 
     .                       x4, y4, z4, xp, yp, zp, r, s, t, u)

      use precision
      use parmmod
      use cgrid
      use ccona
      use comprt, only: iunout

      implicit none

      integer, intent(in) :: icell
      real(dp), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, 
     .                        x4, y4, z4, xp, yp, zp
      real(dp), intent(out) :: r, s, t, u

      real(dp), allocatable, save :: xb(:), yb(:), xcx(:), ycx(:), 
     .                               xce(:), yce(:), a(:), j1(:), j2(:), 
     .                               x0(:), y0(:)
      real(dp), allocatable, save :: za23(:), za31(:), za12(:), 
     .                               y23(:), y31(:), y12(:),   
     .                               x32(:), x13(:), x21(:)   
      real(dp), allocatable, save :: am1(:,:,:) 
      real(dp) :: xp0, yp0, b_xi, b_eta, c_xi, c_eta, xip, etap, twoai,
     .            deter4x4, det, deti, root1, root2
      real(dp) :: b(4,4), ad(4,4)
      logical, allocatable, save :: visited(:)

      real(dp) :: dummy 


      if ((levgeo == 2) .or. (levgeo == 3)) then

! isoparametric quadrilateral (see chapter 23, AFEM)
        if (.not.allocated(visited)) then
          allocate (xb(0:nrad))
          allocate (yb(0:nrad))
          allocate (xcx(0:nrad))
          allocate (ycx(0:nrad))
          allocate (xce(0:nrad))
          allocate (yce(0:nrad))
          allocate (a(0:nrad))
          allocate (j1(0:nrad))
          allocate (j2(0:nrad))
          allocate (x0(0:nrad))
          allocate (y0(0:nrad))
          allocate (visited(0:nrad))
          visited = .false.
         end if

        if (.not.visited(icell)) then
          xb(icell) = x1 - x2 + x3 - x4
          yb(icell) = y1 - y2 + y3 - y4

          xcx(icell) = x1 + x2 - x3 - x4
          ycx(icell) = y1 + y2 - y3 - y4
          
          xce(icell) = x1 - x2 - x3 + x4
          yce(icell) = y1 - y2 - y3 + y4
          
          a(icell) = 0.5_dp * ((x3-x1)*(y4-y2) - (x4-x2)*(y3-y1))
        
          j1(icell) = (x3-x4)*(y1-y2) - (x1-x2)*(y3-y4)
          j2(icell) = (x2-x3)*(y1-y4) - (x1-x4)*(y2-y3)
        
          x0(icell) = 0.25_dp * (x1+x2+x3+x4)
          y0(icell) = 0.25_dp * (y1+y2+y3+y4)

          visited(icell) = .true.
          visited(0) = .false.  ! reset cell 0 for cell outside mesh
        end if

        xp0 = xp - x0(icell)
        yp0 = yp - y0(icell)

        b_xi  =  a(icell) - xp0*yb(icell) + yp0*xb(icell)
        b_eta = -a(icell) - xp0*yb(icell) + yp0*xb(icell)

        c_xi  = xp0*ycx(icell) - yp0*xcx(icell)
        c_eta = xp0*yce(icell) - yp0*xce(icell)

!pb        xip = 2._dp*c_xi / 
!pb     .        (-sqrt(b_xi*b_xi - 2._dp*j1(icell)*c_xi) - b_xi)
!pb        etap = 2._dp*c_eta / 
!pb     .        ( sqrt(b_eta*b_eta + 2._dp*j2(icell)*c_eta) - b_eta)

        root1 =  b_xi*b_xi - 2._dp*j1(icell)*c_xi         
        if ( root1 < 0 ) then
	        write(*,*) '!---------------------------------!'
		write(*,*) 'WARNING NEGATIVE ROOT IN XYZ_TO_RST'
		write(*,*) '!---------------------------------!'
	        root1 = 0._DP
        endif    
	root2 =  b_eta*b_eta + 2._dp*j2(icell)*c_eta 
        if ( root2 < 0 ) then
        	write(*,*) '!---------------------------------!'
		write(*,*) 'WARNING NEGATIVE ROOT IN XYZ_TO_RST'
		write(*,*) '!---------------------------------!'
	        root2 = 0._DP
        endif  
		
	dummy = (-sqrt(root1) - b_xi)
	if ( abs(dummy) > EPS30 ) then	
          xip = 2._dp*c_xi / 
     .        (-sqrt(root1) - b_xi)
        else  
          xip = 0._DP
        endif
! keine Ahnung warum das hier 0 wird, aber ich fangs einfach mal ab !     
     	dummy = ( sqrt(root2) - b_eta)
        if ( abs(dummy) > EPS30 ) then        
           etap = 2._dp*c_eta /dummy 
	else
	   etap = 0._DP
	endif

        r = xip
        s = etap
        t = 0._dp
        u = 0._dp

      else if (levgeo == 4) then

! triangle (see chapter 15, IFEM)
        if (.not.allocated(visited)) then
          allocate (za23(0:nrad))
          allocate (za31(0:nrad))
          allocate (za12(0:nrad))
          allocate (y23(0:nrad))
          allocate (y31(0:nrad))
          allocate (y12(0:nrad))
          allocate (x32(0:nrad))
          allocate (x13(0:nrad))
          allocate (x21(0:nrad))
          allocate (visited(0:nrad))
          visited = .false.
        end if

        if (.not.visited(icell)) then

          y23(icell) = y2-y3
          y31(icell) = y3-y1
          y12(icell) = y1-y2
          
          twoai=1._dp / (x1*y23(icell) + x2*y31(icell) + x3*y12(icell))

          y23(icell) = y23(icell) * twoai
          y31(icell) = y31(icell) * twoai
          y12(icell) = y12(icell) * twoai

          x32(icell) = (x3-x2) * twoai
          x13(icell) = (x1-x3) * twoai
          x21(icell) = (x2-x1) * twoai

          za23(icell) = (x2*y3 - x3*y2) * twoai
          za31(icell) = (x3*y1 - x1*y3) * twoai
          za12(icell) = (x1*y2 - x2*y1) * twoai

          visited(icell) = .true.
          visited(0) = .false.  ! reset cell 0 for cell outside mesh
        end if

        r = za23(icell) + y23(icell)*xp + x32(icell)*yp
        s = za31(icell) + y31(icell)*xp + x13(icell)*yp
        t = za12(icell) + y12(icell)*xp + x21(icell)*yp
        u = 0._dp

      else if (levgeo == 5) then
! tetrahedron (see chapter 15, AFEM)

        if (.not.allocated(visited)) then
          allocate (am1(4,4,0:nrad))
          allocate (visited(0:nrad))
          visited = .false.
        end if

        if (.not.visited(icell)) then

          b(1,1:4) = (/ 1._dp, 1._dp, 1._dp, 1._dp /)
          b(2,1:4) = (/ x1, x2, x3, x4 /)
          b(3,1:4) = (/ y1, y2, y3, y4 /)
          b(4,1:4) = (/ z1, z2, z3, z4 /)
          det = deter4x4(b)
          deti = 1._dp / det
        
          ad(1,1) = x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2)
          ad(1,2) = x1*(y4*z3-y3*z4)+x3*(y1*z4-y4*z1)+x4*(y3*z1-y1*z3)
          ad(1,3) = x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4)+x4*(y1*z2-y2*z1)
          ad(1,4) = x1*(y3*z2-y2*z3)+x2*(y1*z3-y3*z1)+x3*(y2*z1-y1*z2)
          
          ad(2,1) = y2*(z4-z3) + y3*(z2-z4) + y4*(z3-z2)
          ad(2,2) = y1*(z3-z4) + y3*(z4-z1) + y4*(z1-z3)
          ad(2,3) = y1*(z4-z2) + y2*(z1-z4) + y4*(z2-z1)
          ad(2,4) = y1*(z2-z3) + y2*(z3-z1) + y4*(z1-z2)
          
          ad(3,1) = x2*(z3-z4) + x3*(z4-z2) + x4*(z2-z3)
          ad(3,2) = x1*(z4-z3) + x3*(z1-z4) + x4*(z3-z1)
          ad(3,3) = x1*(z2-z4) + x2*(z4-z1) + x4*(z1-z2)
          ad(3,4) = x1*(z3-z2) + x2*(z1-z3) + x3*(z2-z1)
          
          ad(4,1) = x2*(y4-y3) + x3*(y2-y4) + x4*(y3-y2)
          ad(4,2) = x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3)
          ad(4,3) = x1*(y4-y2) + x2*(y1-y4) + x4*(y2-y1)
          ad(4,4) = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
          
          am1(1:4,1:4,icell) = transpose(ad) * deti

          visited(icell) = .true.
          visited(0) = .false.  ! reset cell 0 for cell outside mesh
        end if
        
        r = am1(1,1,icell) + am1(1,2,icell)*xp + am1(1,3,icell)*yp
     .    + am1(1,4,icell)*zp
        s = am1(2,1,icell) + am1(2,2,icell)*xp + am1(2,3,icell)*yp
     .    + am1(2,4,icell)*zp
        t = am1(3,1,icell) + am1(3,2,icell)*xp + am1(3,3,icell)*yp
     .    + am1(3,4,icell)*zp
        u = am1(4,1,icell) + am1(4,2,icell)*xp + am1(4,3,icell)*yp
     .    + am1(4,4,icell)*zp
               
      end if           

      return
      end subroutine xyz_to_rst

      
