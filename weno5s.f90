                 subroutine weno5_solid(nvar2,u,uiphm,uimhp)

           implicit none

           integer :: nvar2  
           real(kind=8) :: u(1:5,1:nvar2), uiphm(1:nvar2), uimhp(1:nvar2)           
           real(kind=8) :: d0, d1, d2, dt2, dt1, dt0
           real(kind=8) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, & 
                               omega0, omega1, omega2
           real(kind=8) :: alphat0, alphat1, alphat2
           real(kind=8) :: sum1, p0, p1, p2, eps
           integer :: i, l

           real(kind=8) :: phi, phi1, phi2, ke
           real(kind=8) :: rho2min, rho2max, rr
           real(kind=8), parameter :: rho2cut = 1.0d-3
           integer :: monotonic

                 

                  monotonic = 1





                 rho2min = minval(u(:,1))
                 rho2max = maxval(u(:,1))
                 rr = rho2max/(rho2min+1.0d-6)


                 if(u(1,1).le.u(2,1).and.u(2,1).le.u(3,1).and.u(3,1).le.u(4,1).and.u(4,1).le.u(5,1)) then
                   monotonic = 1 
                 else if(u(1,1).ge.u(2,1).and.u(2,1).ge.u(3,1).and.u(3,1).ge.u(4,1).and.u(4,1).ge.u(5,1)) then
                   monotonic = 1
                 else
                   monotonic = 0 
                 endif



               
              !if(rr.ge.2.5d0.or.rho2min.le.rho2cut.or.monotonic.eq.0) then
              !if(rr.ge.4.0d0.or.rho2min.le.rho2cut) then
              if(rho2min.le.rho2cut) then
               ! DANGER
                uiphm = u(3,:)
                uimhp = u(3,:)
                return
              endif  

             


                  i = 3


                  d0 = 3.0d0/10.0d0
                  d1 = 6.0d0/10.0d0
                  d2 = 1.0d0/10.0d0
                  dt2 = 3.0d0/10.0d0
                  dt1 = 6.0d0/10.0d0
                  dt0 = 1.0d0/10.0d0

                  eps = 1.0d-12


                  uiphm = 0.0d0
                  uimhp = 0.0d0


      do l = 1, nvar2

              ! compute u(i+1/2)-

           beta0 = 13.0d0/12.0d0*(u(i,l) - 2.0d0*u(i+1,l) + &
            u(i+2,l))**2.0d0 + 1.0d0/4.0d0*(3.0d0*u(i,l) - 4.0d0*u(i+1,l) + u(i+2,l))**2.0d0

           beta1 = 13.0d0/12.0d0*(u(i-1,l) - 2.0d0*u(i,l) + &
                u(i+1,l))**2.0d0 + 1.0d0/4.0d0*(u(i-1,l) - u(i+1,l))**2.0d0

           beta2 = 13.0d0/12.0d0*(u(i-2,l) - 2.0d0*u(i-1,l) + u(i,l))**2.0d0 + &
                 1.0d0/4.0d0*(u(i-2,l) - 4.0d0*u(i-1,l) + 3.0d0*u(i,l))**2.0d0

           alpha0 = d0/(eps + beta0)**2.0d0
           alpha1 = d1/(eps + beta1)**2.0d0
           alpha2 = d2/(eps + beta2)**2.0d0

           sum1 = alpha0 + alpha1 + alpha2

           omega0 = alpha0/sum1
           omega1 = alpha1/sum1
           omega2 = alpha2/sum1


           p0 = 1.0d0/3.0d0*u(i,l) + 5.0d0/6.0d0*u(i+1,l) - 1.0d0/6.0d0*u(i+2,l)
           p1 = -1.0d0/6.0d0*u(i-1,l) + 5.0d0/6.0d0*u(i,l) + 1.0d0/3.0d0*u(i+1,l)
           p2 = 1.0d0/3.0d0*u(i-2,l) - 7.0d0/6.0d0*u(i-1,l) + 11.0d0/6.0d0*u(i,l)


           uiphm(l) = omega0*p0 + omega1*p1 + omega2*p2

                  ! compute u(i-1/2)+

           alphat0 = dt0/(eps + beta0)**2.0d0
           alphat1 = dt1/(eps + beta1)**2.0d0
           alphat2 = dt2/(eps + beta2)**2.0d0

           sum1 = alphat0 + alphat1 + alphat2

           omega0 = alphat0/sum1
           omega1 = alphat1/sum1
           omega2 = alphat2/sum1

           p0 = 11.0d0/6.0d0*u(i,l) - 7.0d0/6.0d0*u(i+1,l) + 1.0d0/3.0d0*u(i+2,l)
           p1 = 1.0d0/3.0d0*u(i-1,l) + 5.0d0/6.0d0*u(i,l) - 1.0d0/6.0d0*u(i+1,l)
           p2 = -1.0d0/6.0d0*u(i-2,l) + 5.0d0/6.0d0*u(i-1,l) + 1.0d0/3.0d0*u(i,l)


           uimhp(l) = omega0*p0 + omega1*p1 + omega2*p2






            
                  
            if(1.eq.2) then

                      ! slope limiter
                 
              if(u(i,l).ne.u(i-1,l)) then
               phi1 = (u(i+1,l)-u(i,l))/(u(i,l)-u(i-1,l))
               phi2 = (uiphm(l)-u(i,l))/(u(i,l)-u(i-1,l))           
               phi = max(0.0d0,min(2.0d0,2.0d0*phi1,2.0d0*phi2))
                                       
               uiphm(l) = u(i,l) + 0.5d0*(u(i,l) - u(i-1,l))*phi
              else
               uiphm(l) = u(i,l)
              endif 



              if(u(i,l).ne.u(i+1,l)) then   
               phi1 = (u(i-1,l)-u(i,l))/(u(i,l)-u(i+1,l))
               phi2 = (uimhp(l)-u(i,l))/(u(i,l)-u(i+1,l))           
               phi = max(0.0d0,min(2.0d0,2.0d0*phi1,2.0d0*phi2))
                        
               uimhp(l) = u(i,l) - 0.5d0*(u(i+1,l) - u(i,l))*phi
              else
               uimhp(l) = u(i,l)
              endif

            endif        



                      

      enddo


            
                  
             !  NUMBER DENSITY FIRST ORDER
              uiphm(6) = u(3,6)
              uimhp(6) = u(3,6)

      


           return
           end subroutine

!---------------------------------------------------------------------------


