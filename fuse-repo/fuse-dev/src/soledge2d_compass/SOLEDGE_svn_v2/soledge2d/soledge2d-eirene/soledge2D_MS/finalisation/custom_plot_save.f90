subroutine custom_plot_save(n_ite)
  use Mphysics
  use all_variables, only : global_variables, global_parameters, zones, reference_parameters, custom_plots, N_custom_plots
  implicit none
  integer*4,intent(in) :: n_ite
  integer*4 :: i,j,k,n,k1,Nx,Nz
  character(50) :: filename
  real*8 :: Mtheta,Mdrift
  real*8 :: R0,c0,rs0
  real*8 :: upara,G,cs,uE,ctt
  real*8 :: Total_p, Static_p,A
  R0 = reference_parameters%geometry%R0
  c0 = reference_parameters%fields%c0
  rs0 = reference_parameters%geometry%rs0
  if(modulo(dble(n_ite),dble(global_parameters%N_iterations)/1000.D0).eq.0) then
     do n=1,N_custom_plots
        write(filename,"(A7,I0)") "custom_",n
        if(custom_plots(n)%type.le.2) then
           open(unit=10,file=trim(filename),status='unknown')
        else
           open(unit=10,file=trim(filename),status='unknown',access='append')
        end if
        if(custom_plots(n)%type.eq.1) then ! parallel
           k1=1
           k=custom_plots(n)%zones(k1)
           Nx=zones(k)%mesh%Nx
           Nz=zones(k)%mesh%Nz
           i=custom_plots(n)%coord
           if(custom_plots(n)%nzones.eq.1) then
              do j=0,Nz+1
                 upara=zones(k)%species(1)%var(1)%Gamma(i,j)/zones(k)%species(1)%var(1)%density(i,j)
                 uE=zones(k)%species(1)%drifts%uEt(i,j)+zones(k)%species(1)%drifts%uBt(i,j)
                 G=zones(k)%metric_coefficients%G(i,j)
                 ctt=zones(k)%metric_coefficients%ctt(i,j)
                 cs=sqrt((zones(k)%species(1)%var(1)%temperature(i,j)+zones(k)%species(0)%var(1)%temperature(i,j))&
                      /zones(k)%species(1)%element%mass)
                 Mtheta=(upara*G+uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 Mdrift=(uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 A=G/sqrt(ctt)*rs0/(2.D0*pi*R0)
                 Total_p=zones(k)%species(1)%var(1)%density(i,j)*((zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A+&
                      upara*upara*A+upara*uE)
                 Static_p=zones(k)%species(1)%var(1)%density(i,j)*(zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A
                 if(zones(k)%masks%chi2(i,j).eq.0) then
                    write(10,100) zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0,& !n
                         zones(k)%species(1)%var(1)%Mach(i,j),&                                            !M
                         zones(k)%species(0)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Te
                         zones(k)%species(1)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Ti
                         zones(k)%electric_fields(1)%phi(i,j)*reference_parameters%fields%phi0,&           !phi
                         zones(k)%electric_fields(1)%vorticity(i,j)*reference_parameters%fields%W0,&       !W
                         zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j),&
                         Mtheta,&                                                                          !Mtheta
                         uE*c0,&                                                                           !u drift poloidal
                         Mdrift,&
                         uE*c0+upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                         Total_p, Static_p
                 else
                    write(10,100) 0., 0., 0., 0., 0., 0., zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j), 0., 0., 0., 0., 0., 0.
                 end if
              end do
           else
              do j=0,Nz
                 upara=zones(k)%species(1)%var(1)%Gamma(i,j)/zones(k)%species(1)%var(1)%density(i,j)
                 uE=zones(k)%species(1)%drifts%uEt(i,j)+zones(k)%species(1)%drifts%uBt(i,j)
                 G=zones(k)%metric_coefficients%G(i,j)
                 ctt=zones(k)%metric_coefficients%ctt(i,j)
                 cs=sqrt((zones(k)%species(1)%var(1)%temperature(i,j)+zones(k)%species(0)%var(1)%temperature(i,j))&
                      /zones(k)%species(1)%element%mass)
                 Mtheta=(upara*G+uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 Mdrift=(uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 A=G/sqrt(ctt)*rs0/(2.D0*pi*R0)
                 Total_p=zones(k)%species(1)%var(1)%density(i,j)*((zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A+&
                      upara*upara*A+upara*uE)
                 Static_p=zones(k)%species(1)%var(1)%density(i,j)*(zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A
                 if(zones(k)%masks%chi2(i,j).eq.0) then
                    write(10,100) zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0,& !n
                         zones(k)%species(1)%var(1)%Mach(i,j),&                                            !M
                         zones(k)%species(0)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Te
                         zones(k)%species(1)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Ti
                         zones(k)%electric_fields(1)%phi(i,j)*reference_parameters%fields%phi0,&           !phi
                         zones(k)%electric_fields(1)%vorticity(i,j)*reference_parameters%fields%W0,&       !W
                         zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j),&
                         Mtheta,&                                                                          !Mtheta
                         uE*c0,&                                                                           !u drift poloidal
                         Mdrift,&
                         uE*c0+upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                         uE*c0,&
                         upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                         Total_p, Static_p
                 else
                    write(10,100) 0., 0., 0., 0., 0., 0., zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j), 0., 0., 0., 0., 0., 0.
                 end if
              end do
           end if
           do k1=2,custom_plots(n)%nzones-1
              k=custom_plots(n)%zones(k1)
              Nx=zones(k)%mesh%Nx
              Nz=zones(k)%mesh%Nz
              i=custom_plots(n)%coord
              do j=1,Nz
                 upara=zones(k)%species(1)%var(1)%Gamma(i,j)/zones(k)%species(1)%var(1)%density(i,j)
                 uE=zones(k)%species(1)%drifts%uEt(i,j)+zones(k)%species(1)%drifts%uBt(i,j)
                 G=zones(k)%metric_coefficients%G(i,j)
                 ctt=zones(k)%metric_coefficients%ctt(i,j)
                 cs=sqrt((zones(k)%species(1)%var(1)%temperature(i,j)+zones(k)%species(0)%var(1)%temperature(i,j))&
                      /zones(k)%species(1)%element%mass)
                 Mtheta=(upara*G+uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 Mdrift=(uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 A=G/sqrt(ctt)*rs0/(2.D0*pi*R0)
                 Total_p=zones(k)%species(1)%var(1)%density(i,j)*((zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A+&
                      upara*upara*A+upara*uE)
                 Static_p=zones(k)%species(1)%var(1)%density(i,j)*(zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A
                 if(zones(k)%masks%chi2(i,j).eq.0) then
                    write(10,100) zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0,& !n
                         zones(k)%species(1)%var(1)%Mach(i,j),&                                            !M
                         zones(k)%species(0)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Te
                         zones(k)%species(1)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Ti
                         zones(k)%electric_fields(1)%phi(i,j)*reference_parameters%fields%phi0,&           !phi
                         zones(k)%electric_fields(1)%vorticity(i,j)*reference_parameters%fields%W0,&       !W
                         zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j),&
                         Mtheta,&                                                                          !Mtheta
                         uE*c0,&                                                                           !u drift poloidal
                         Mdrift,&
                         uE*c0+upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                         uE*c0,&
                         upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                         Total_p, Static_p                                       !Mtheta
                 else
                    write(10,100) 0., 0., 0., 0., 0., 0., zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j), 0., 0., 0., 0., 0., 0.
                 end if
              end do
           end do
           if(custom_plots(n)%nzones.gt.1) then
              k1=custom_plots(n)%nzones
              k=custom_plots(n)%zones(k1)
              Nx=zones(k)%mesh%Nx
              Nz=zones(k)%mesh%Nz
              i=custom_plots(n)%coord
              do j=1,Nz+1
                 upara=zones(k)%species(1)%var(1)%Gamma(i,j)/zones(k)%species(1)%var(1)%density(i,j)
                 uE=zones(k)%species(1)%drifts%uEt(i,j)+zones(k)%species(1)%drifts%uBt(i,j)
                 G=zones(k)%metric_coefficients%G(i,j)
                 ctt=zones(k)%metric_coefficients%ctt(i,j)
                 cs=sqrt((zones(k)%species(1)%var(1)%temperature(i,j)+zones(k)%species(0)%var(1)%temperature(i,j))&
                      /zones(k)%species(1)%element%mass)
                 Mtheta=(upara*G+uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 Mdrift=(uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                 A=G/sqrt(ctt)*rs0/(2.D0*pi*R0)
                 Total_p=zones(k)%species(1)%var(1)%density(i,j)*((zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A+&
                      upara*upara*A+upara*uE)
                 Static_p=zones(k)%species(1)%var(1)%density(i,j)*(zones(k)%species(1)%var(1)%temperature(i,j)&
                      +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A
                 if(zones(k)%masks%chi2(i,j).eq.0) then
                    write(10,100) zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0,& !n
                         zones(k)%species(1)%var(1)%Mach(i,j),&                                            !M
                         zones(k)%species(0)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Te
                         zones(k)%species(1)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Ti
                         zones(k)%electric_fields(1)%phi(i,j)*reference_parameters%fields%phi0,&           !phi
                         zones(k)%electric_fields(1)%vorticity(i,j)*reference_parameters%fields%W0,&       !W
                         zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j),&
                         Mtheta,&                                                                          !Mtheta
                         uE*c0,&                                                                           !u drift poloidal
                         Mdrift,&
                         uE*c0+upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                         uE*c0,&
                         upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                         Total_p, Static_p                                                                            !Mthet
                 else
                    write(10,100) 0., 0., 0., 0., 0., 0., zones(k)%mesh%Rgeom(i,j),&
                         zones(k)%mesh%Zgeom(i,j), 0., 0., 0., 0., 0., 0.
                 end if
              end do
           end if
        else
           if(custom_plots(n)%type.eq.2) then ! perp
              do k1=1,custom_plots(n)%nzones
                 k=custom_plots(n)%zones(k1)
                 Nx=zones(k)%mesh%Nx
                 Nz=zones(k)%mesh%Nz
                 j=custom_plots(n)%coord
                 do i=1,Nx
                    upara=zones(k)%species(1)%var(1)%Gamma(i,j)/zones(k)%species(1)%var(1)%density(i,j)
                    uE=zones(k)%species(1)%drifts%uEt(i,j)+zones(k)%species(1)%drifts%uBt(i,j)
                    G=zones(k)%metric_coefficients%G(i,j)
                    ctt=zones(k)%metric_coefficients%ctt(i,j)
                    cs=sqrt((zones(k)%species(1)%var(1)%temperature(i,j)+zones(k)%species(0)%var(1)%temperature(i,j))&
                         /zones(k)%species(1)%element%mass)
                    Mtheta=(upara*G+uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                    Mdrift=(uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
                    A=G/sqrt(ctt)*rs0/(2.D0*pi*R0)
                    Total_p=zones(k)%species(1)%var(1)%density(i,j)*((zones(k)%species(1)%var(1)%temperature(i,j)&
                         +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A+&
                         upara*upara*A+upara*uE)
                    Static_p=zones(k)%species(1)%var(1)%density(i,j)*(zones(k)%species(1)%var(1)%temperature(i,j)&
                         +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A
                    if(zones(k)%masks%chi2(i,j).eq.0) then
                       write(10,100) zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0,& !n
                            zones(k)%species(1)%var(1)%Mach(i,j),&                                            !M
                            zones(k)%species(0)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Te
                            zones(k)%species(1)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Ti
                            zones(k)%electric_fields(1)%phi(i,j)*reference_parameters%fields%phi0,&           !phi
                            zones(k)%electric_fields(1)%vorticity(i,j)*reference_parameters%fields%W0,&       !W
                            zones(k)%mesh%Rgeom(i,j),&
                            zones(k)%mesh%Zgeom(i,j),&
                            Mtheta,&                                                                          !Mtheta
                            uE*c0,&                                                                           !u drift poloidal
                            Mdrift,&
                            uE*c0+upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                            uE*c0,&
                            upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                            Total_p, Static_p                                          !Mtheta
                    else
                       write(10,100) 0., 0., 0., 0., 0., 0., zones(k)%mesh%Rgeom(i,j),&
                            zones(k)%mesh%Zgeom(i,j), 0., 0., 0., 0., 0., 0.
                    end if
                 end do
              end do
           else
              !temporal
              k=custom_plots(n)%zones(1)
              i=custom_plots(n)%coord
              j=custom_plots(n)%coord2
              upara=zones(k)%species(1)%var(1)%Gamma(i,j)/zones(k)%species(1)%var(1)%density(i,j)
              uE=zones(k)%species(1)%drifts%uEt(i,j)+zones(k)%species(1)%drifts%uBt(i,j)
              G=zones(k)%metric_coefficients%G(i,j)
              ctt=zones(k)%metric_coefficients%ctt(i,j)
              cs=sqrt((zones(k)%species(1)%var(1)%temperature(i,j)+zones(k)%species(0)%var(1)%temperature(i,j))&
                   /zones(k)%species(1)%element%mass)
              Mtheta=(upara*G+uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
              Mdrift=(uE*sqrt(ctt)*(2.d0*pi*R0/rs0))/(Cs*G)
              A=G/sqrt(ctt)*rs0/(2.D0*pi*R0)
              Total_p=zones(k)%species(1)%var(1)%density(i,j)*((zones(k)%species(1)%var(1)%temperature(i,j)&
                   +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A+&
                   upara*upara*A+upara*uE)
              Static_p=zones(k)%species(1)%var(1)%density(i,j)*(zones(k)%species(1)%var(1)%temperature(i,j)&
                   +zones(k)%species(0)%var(1)%temperature(i,j))/zones(k)%species(1)%element%mass*A
              if(zones(k)%masks%chi2(i,j).eq.0) then
                 write(10,100) zones(k)%species(1)%var(1)%density(i,j)*reference_parameters%fields%n0,& !n
                      zones(k)%species(1)%var(1)%Mach(i,j),&                                            !M
                      zones(k)%species(0)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Te
                      zones(k)%species(1)%var(1)%temperature(i,j)*reference_parameters%fields%T0eV,&    !Ti
                      zones(k)%electric_fields(1)%phi(i,j)*reference_parameters%fields%phi0,&           !phi
                      zones(k)%electric_fields(1)%vorticity(i,j)*reference_parameters%fields%W0,&       !W
                      zones(k)%mesh%Rgeom(i,j),&
                      zones(k)%mesh%Zgeom(i,j),&
                      Mtheta,&                                                                          !Mtheta
                      uE*c0,&                                                                           !u drift poloidal
                      Mdrift,&
                      uE*c0+upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                      uE*c0,&
                      upara*G/sqrt(ctt)*c0*rs0/(2.d0*pi*R0),&
                      Total_p, Static_p                                                                            !Mtheta
              else
                 write(10,100) 0., 0., 0., 0., 0., 0., zones(k)%mesh%Rgeom(i,j),&
                      zones(k)%mesh%Zgeom(i,j), 0., 0., 0., 0., 0., 0.
              end if
           end if
        end if
        close(10)
     end do
100  format(512es15.7)
  end if
end subroutine custom_plot_save
