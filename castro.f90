
         subroutine castro_2d(priml,primr,sweep,flux)

              use thermo

      implicit none   

        real(kind=8) :: rhol, rhor, pl, pr, ul, ur, vl, vr, iel,  &
                   ier, ysl(1:5), ysr(1:5), priml(1:10), primr(1:10)
        real(kind=8) :: term1, term2, pstar, ustar
        real(kind=8), dimension(1:10) :: ustarl, ustarr, left, right, godunov
        real(kind=8) :: flux(1:9)
        integer :: sweep, ns   
        real(kind=8) :: ql, qr, cl, cr, gaml, gamr, sl, sr, hl, hr
        real(kind=8) :: laml, lamr, lamlstar, lamrstar, qn, nx, ny, energy, ke
        real(kind=8) :: clstar, crstar, sigma   

               

               rhol = priml(1) 
               rhor = primr(1)
               ul = priml(2)
               ur = primr(2)
               vl = priml(3)
               vr = primr(3)
               iel = priml(4)
               ier = primr(4)
               pl = priml(5)
               pr = primr(5)
               ysl = priml(6:10)
               ysr = primr(6:10)  
                              
         
         
            


          call compute_sound(pl,rhol,ysl,cl) 
          call compute_sound(pr,rhor,ysr,cr) 



               if(abs(1.0d0-sum(ysl(:))).gt.1.0d-6) then
                 print*, 'ysl not adding to one ', sum(ysl), ysl
                 stop
               endif
               if(abs(1.0d0-sum(ysr(:))).gt.1.0d-6) then
                 print*, 'ysr not adding to one ', sum(ysr), ysr
                 stop
               endif
                   



             
            if(sweep.eq.1) then
              ql = ul
              qr = ur
              nx = 1.0d0
              ny = 0.0d0 
            else if(sweep.eq.2) then
              ql = vl
              qr = vr
              nx = 0.0d0              
              ny = 1.0d0
            else
             print*, 'bug in sweep ', sweep
             stop
            endif

 
             gaml = cl*cl*rhol/pl
             gamr = cr*cr*rhor/pr

             sl = sqrt(gaml*pl*rhol)
             sr = sqrt(gamr*pr*rhor)     



             term1 = sl*pr + sr*pl + sl*sr*(ql-qr)
             pstar = term1/(sl+sr)


             term2 = sl*ql + sr*qr + pl-pr 
             ustar = term2/(sl+sr)


             ustarl(1) = rhol + (pstar-pl)/cl/cl
             ustarr(1) = rhor + (pstar-pr)/cr/cr

             ustarl(5) = pstar
             ustarr(5) = pstar
 
             ustarl(6:10) = ysl
             ustarr(6:10) = ysr


            if(sweep.eq.1) then  
             ustarl(2) = ustar
             ustarr(2) = ustar
             ustarl(3) = vl
             ustarr(3) = vr     
            else if(sweep.eq.2) then
             ustarl(2) = ul 
             ustarr(2) = ur  
             ustarl(3) = ustar
             ustarr(3) = ustar
            else
              print*, 'bug in sweep ', sweep
              stop 
            endif 


            ! 4 is actually rho*e
            
            hl = iel + pl/rhol
            ustarl(4) = rhol*iel + (pstar-pl)*hl/cl/cl 
 
            hr = ier + pr/rhor
            ustarr(4) = rhor*ier + (pstar-pr)*hr/cr/cr 



            clstar = sqrt(gaml*pstar/ustarl(1))
            crstar = sqrt(gamr*pstar/ustarr(1)) 


            left = priml
            left(4) = left(4)*left(1)
 
            right = primr  
            right(4) = right(4)*right(1)          



            if(ustar.gt.0.0d0) then

              laml = ql - cl
              lamlstar = ustar - clstar

              if(pstar.ge.pl) then
                ! shock to the left
   
                sigma = 0.5d0*(laml + lamlstar) 

                if(sigma.ge.0.0d0) then
                 godunov = left                 
                else
                 godunov = ustarl
                endif 
  
              else
                ! rarefaction to the left 

                 if(laml.le.0.0d0.and.lamlstar.le.0.0d0) then
                   godunov = ustarl
                 else if(laml.ge.0.0d0.and.lamlstar.ge.0.0d0) then
                   godunov = left                   
                 else 
                   ! transonic rarefaction
                   godunov = left*(lamlstar/(lamlstar-laml)) - ustarl*(laml/(lamlstar-laml))
                 endif

              endif 
 
            else if(ustar.lt.0.0d0) then

              lamr = qr + cr
              lamrstar = ustar + crstar

              if(pstar.ge.pr) then
                ! shock to the right
   
                sigma = 0.5d0*(lamr + lamrstar) 

                if(sigma.gt.0.0d0) then
                 godunov = ustarr                 
                else
                 godunov = right
                endif 
  
              else
                ! rarefaction to the right 

                 if(lamr.ge.0.0d0.and.lamrstar.ge.0.0d0) then
                   godunov = ustarr
                 else if(lamr.le.0.0d0.and.lamrstar.le.0.0d0) then
                   godunov = right                   
                 else 
                   ! transonic rarefaction
                   godunov = right*(lamrstar/(lamrstar-lamr)) - ustarr*(lamr/(lamrstar-lamr))
                 endif

              endif 


            else if(ustar.eq.0.0d0) then
              godunov = 0.5d0*(ustarl+ustarr)
            endif



         
           qn = godunov(2)*nx + godunov(3)*ny

           flux(1) = godunov(1)*qn
           flux(2) = godunov(1)*qn*godunov(2) + godunov(5)*nx
           flux(3) = godunov(1)*qn*godunov(3) + godunov(5)*ny

           ke = 0.5d0*(godunov(2)**2.0d0 + godunov(3)**2.0d0)
           energy = godunov(4) + godunov(1)*ke
           flux(4) = qn*(energy + godunov(5))   

           flux(5:9) = godunov(1)*godunov(6:10)



          
               return
               end subroutine
!----------------------------------------------------------------------------------




 
