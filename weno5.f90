                 subroutine weno5_gas(nprim,u,uiphm,uimhp)

           implicit none

           integer :: nprim  
           real(kind=8) :: u(1:5,1:nprim), uiphm(1:nprim), uimhp(1:nprim)          
           real(kind=8) :: d0, d1, d2, dt2, dt1, dt0
           real(kind=8) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, & 
                               omega0, omega1, omega2
           real(kind=8) :: alphat0, alphat1, alphat2
           real(kind=8) :: sum1, p0, p1, p2, eps
           integer :: i, l, ns

           real(kind=8) :: phi, phi1, phi2, ke, sumyk
           real(kind=8), dimension(1:2) :: w0_rho, w1_rho, w2_rho


               
                   
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


!------------------------------------------------------------------------------------------------
      l = 1

319           continue


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


                  if(l.eq.1) then
                   w0_rho(1) = omega0
                   w1_rho(1) = omega1
                   w2_rho(1) = omega2 
                  endif


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


                 if(l.eq.1) then
                   w0_rho(2) = omega0
                   w1_rho(2) = omega1
                   w2_rho(2) = omega2 
                  endif


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




          if(l.lt.3) then
            l = l + 1
            goto 319 
          endif 

          if(l.eq.3) then
            l = 5
            goto 319
          endif


!---------------------------------------------------------------------------------- 
            ! MASS FRACTIONS


                  
          do l = 6, 10


           p0 = 1.0d0/3.0d0*u(i,l) + 5.0d0/6.0d0*u(i+1,l) - 1.0d0/6.0d0*u(i+2,l)
           p1 = -1.0d0/6.0d0*u(i-1,l) + 5.0d0/6.0d0*u(i,l) + 1.0d0/3.0d0*u(i+1,l)
           p2 = 1.0d0/3.0d0*u(i-2,l) - 7.0d0/6.0d0*u(i-1,l) + 11.0d0/6.0d0*u(i,l)


           uiphm(l) = w0_rho(1)*p0 + w1_rho(1)*p1 + w2_rho(1)*p2


           p0 = 11.0d0/6.0d0*u(i,l) - 7.0d0/6.0d0*u(i+1,l) + 1.0d0/3.0d0*u(i+2,l)
           p1 = 1.0d0/3.0d0*u(i-1,l) + 5.0d0/6.0d0*u(i,l) - 1.0d0/6.0d0*u(i+1,l)
           p2 = -1.0d0/6.0d0*u(i-2,l) + 5.0d0/6.0d0*u(i-1,l) + 1.0d0/3.0d0*u(i,l)


           uimhp(l) = w0_rho(2)*p0 + w1_rho(2)*p1 + w2_rho(2)*p2
            

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


          enddo 


           sumyk = sum(uiphm(6:10))
           uiphm(6:10) = uiphm(6:10)/sumyk 

           sumyk = sum(uimhp(6:10))
           uimhp(6:10) = uimhp(6:10)/sumyk 


!---------------------------------------------------------------------------------- 
          ! INTERNAL ENERGY

                    

             call compute_ie(uiphm(1), uiphm(5), uiphm(6:10), uiphm(4)) 
             call compute_ie(uimhp(1), uimhp(5), uimhp(6:10), uimhp(4)) 


                    
!---------------------------------------------------------------------------------- 


           return
           end subroutine

!---------------------------------------------------------------------------


