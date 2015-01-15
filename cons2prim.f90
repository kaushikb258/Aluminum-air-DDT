

                subroutine cons_2_prim(rho2cut,ncons,nvar2,nspec,cons,cons2,prim,prim2)

          use thermo


         implicit none

            
 
             integer :: ncons, nvar2,nspec, ns
             real(kind=8) :: rho2cut, cons(1:ncons), cons2(1:nvar2), prim(1:ncons+1), prim2(1:nvar2)
             real(kind=8) :: ke, sumyk


                 prim(1) = cons(1)
                 prim(2) = cons(2)/cons(1)
                 prim(3) = cons(3)/cons(1)


                 ke = prim(1)*(prim(2)*prim(2) + prim(3)*prim(3))/2.0d0
                 prim(4) = (cons(4)-ke)/cons(1)

              

                 ! species
            do ns = 1, nspec
             prim(5+ns) = cons(4+ns)/cons(1)
             prim(5+ns) = max(min(prim(5+ns),1.0d0),0.0d0)          
            enddo

            sumyk = sum(prim(6:5+nspec)) 
            prim(6:5+nspec) = prim(6:5+nspec)/sumyk


                 ! pressure              
         call compute_pres(prim(1),prim(6:5+nspec),prim(4),prim(5)) 
              


                    ! SOLID


                    if(cons2(1).gt.rho2cut) then
                 prim2(1) = cons2(1)
                 prim2(2) = cons2(2)/cons2(1)
                 prim2(3) = cons2(3)/cons2(1)

                 ke = prim2(1)*(prim2(2)*prim2(2) + prim2(3)*prim2(3))/2.0d0
                 prim2(4) = (cons2(4)-ke)/cons2(1)  
                 prim2(5) = cons2(5)/cons2(1)
                 prim2(6) = cons2(6)  

                    if(prim2(6).le.0.0d0) then
                     !print*, 'num < 0; rho2 > 0 ', prim2(6), prim2(1)
                     prim2(6) = 0.0d0
                     !stop
                    endif

                    else
                 cons2(1:nvar2) = 0.0d0
                 prim2(1:nvar2) = 0.0d0
                    endif




          return
        end subroutine

!-------------------------------------------------------------- 
