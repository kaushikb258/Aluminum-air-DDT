
         subroutine hllc_2d(priml,primr,sweep,flux)

              use thermo

      implicit none   

        real(kind=8) :: rhol, rhor, pl, pr, ul, ur, vl, vr, iel,  &
                   ier, ysl(1:5), ysr(1:5), priml(1:10), primr(1:10)
        real(kind=8) :: spl, spr, term1, term2, sm, pstar
        real(kind=8) :: ustarl(1:9), ustarr(1:9), matl(1:9), matr(1:9)
        real(kind=8) :: flux(1:9), energyl, energyr, omegal, omegar
        integer :: sweep, ns   
        real(kind=8) :: ql, qr, cl, cr
        real(kind=8) :: nx, ny
        real(kind=8), dimension(1:9) :: fl, fr, consl, consr
       

! HLLC FLUX FOR THE GAS


               

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

            SPL = ql - cl
            SPR = qr + cr
            TERM1 = RHOR*qr*(SPR-qr) - RHOL*ql*(SPL-ql)+pl-pr
            TERM2 = RHOR*(SPR-qr) - RHOL*(SPL-ql)
            SM = TERM1/TERM2

           PSTAR = RHOL*(ql-SPL)*(ql-SM) + pl


           omegal = rhol*(spl-ql)/(spl-sm)
           omegar = rhor*(spr-qr)/(spr-sm)



           matl(1) = 1.0d0
           matr(1) = 1.0d0
           if(sweep.eq.1) then
            matl(2) = sm
            matl(3) = vl
            matr(2) = sm
            matr(3) = vr 
           else if(sweep.eq.2) then
            matl(2) = ul
            matl(3) = sm
            matr(2) = ur
            matr(3) = sm 
           else
            print*, 'hllc sweep ', sweep
            stop 
           endif  


           energyl = iel + 0.5d0*(ul*ul + vl*vl)
           energyr = ier + 0.5d0*(ur*ur + vr*vr)


           matl(4) = energyl + (sm-ql)*(sm+pl/rhol/(spl-ql)) 
           matr(4) = energyr + (sm-qr)*(sm+pr/rhor/(spr-qr)) 

           do ns = 1, nspeci
             matl(4+ns) = ysl(ns)
             matr(4+ns) = ysr(ns)
           enddo 



           ustarl(:) = omegal*matl(:) 
           ustarr(:) = omegar*matr(:)

            

           fl(1) = rhol*ql
           fr(1) = rhor*qr
           fl(2) = rhol*ql*ul + pl*nx
           fr(2) = rhor*qr*ur + pr*nx
           fl(3) = rhol*ql*vl + pl*ny
           fr(3) = rhor*qr*vr + pr*ny

           fl(4) = rhol*ql*energyl + pl*ql
           fr(4) = rhor*qr*energyr + pr*qr 

           do ns = 1, nspeci
            fl(4+ns) = rhol*ql*ysl(ns) 
            fr(4+ns) = rhor*qr*ysr(ns) 
           enddo 


           consl(1) = rhol
           consr(1) = rhor
           consl(2) = rhol*ul
           consr(2) = rhor*ur
           consl(3) = rhol*vl
           consr(3) = rhor*vr
           consl(4) = rhol*energyl 
           consr(4) = rhor*energyr

           do ns = 1, nspeci
            consl(4+ns) = rhol*ysl(ns)
            consr(4+ns) = rhor*ysr(ns) 
           enddo


          
           if(spl.ge.0.0d0) then
             flux = fl
           else if(spl.le.0.0d0.and.sm.ge.0.0d0) then
             flux = fl + spl*(ustarl-consl)
           else if(sm.le.0.0d0.and.spr.ge.0.0d0) then
             flux = fr + spr*(ustarr-consr) 
           else if(spr.le.0.0d0) then
             flux = fr
           else
            print*, 'hllc bug '
            stop
           endif   
          



          
               return
               end subroutine
!----------------------------------------------------------------------------------




 
