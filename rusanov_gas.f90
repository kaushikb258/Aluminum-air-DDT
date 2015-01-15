
         subroutine rusanov_gas_2d(priml,primr,sweep,flux)

              use thermo

      implicit none   

        real(kind=8) :: rhol, rhor, pl, pr, ul, ur, vl, vr, iel,  &
                   ier, ysl(1:5), ysr(1:5), priml(1:10), primr(1:10)
        real(kind=8) :: spl, spr, term1, term2, sm, pstar
        real(kind=8) :: rhostarl, rhostarr, ustarl(1:9), ustarr(1:9)
        real(kind=8) :: flux(1:9), energyl, energyr
        integer :: sweep, ns   
        real(kind=8) :: ql, qr, qtl, qtr, omegal, omegar, cl, cr
        real(kind=8), dimension(1:9) :: consl, consr, fl, fr   
        real(kind=8) :: lamtil, nx, ny, ke   
        

! RUSANOV FLUX FOR THE GAS


               

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

            

           consl(1) = rhol
           consr(1) = rhor
           consl(2) = rhol*ul
           consr(2) = rhor*ur
           consl(3) = rhol*vl
           consr(3) = rhor*vr
           
           ke = 0.5d0*(ul*ul + vl*vl)
           consl(4) = rhol*(iel + ke) 
           ke = 0.5d0*(ur*ur + vr*vr)
           consr(4) = rhor*(ier + ke)     

           do ns = 1, nspeci
            consl(4+ns) = rhol*ysl(ns)
            consr(4+ns) = rhor*ysr(ns) 
           enddo


           fl(1) = rhol*ql
           fr(1) = rhor*qr
           fl(2) = rhol*ql*ul + pl*nx
           fr(2) = rhor*qr*ur + pr*nx 
           fl(3) = rhol*ql*vl + pl*ny
           fr(3) = rhor*qr*vr + pr*ny  

           ke = 0.5d0*(ul*ul + vl*vl)
           fl(4) = rhol*ql*(iel + ke) + pl*ql 
           ke = 0.5d0*(ur*ur + vr*vr)
           fr(4) = rhor*qr*(ier + ke) + pr*qr

           do ns = 1, nspeci
            fl(4+ns) = rhol*ql*ysl(ns) 
            fr(4+ns) = rhor*qr*ysr(ns) 
           enddo   



            lamtil = max(abs(ql)+cl,abs(qr)+cr)


            do ns = 1, 4+nspeci
             flux(ns) = 0.5d0*(fl(ns)+fr(ns)) - 0.5d0*lamtil*(consr(ns)-consl(ns)) 
            enddo



               return
               end subroutine
!----------------------------------------------------------------------------------




 
