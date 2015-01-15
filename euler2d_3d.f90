                               program euler2d

         

           USE THERMO


         implicit none 


              integer, parameter :: imax = 4001 
              integer, parameter :: jmax = 501
              integer, parameter :: nmax =  15000 ! number of time steps
              integer, parameter :: ofile = 200 ! how often to output vtk file             

              real(kind=8) :: dx, dy, x(1:imax), y(1:jmax), xlen, ylen
              real(kind=8) :: time, dt, cfl
              integer :: i, j, n, ns, irestart, l

              integer, parameter :: nghost = 5 ! number of ghost cells
              integer, parameter :: nspec = 5 ! number of species
              integer, parameter :: ncons = 4+nspec ! number of conservative variables
              integer, parameter :: nprim = 5+nspec ! number of primitive variables
              integer, parameter :: nvar2 = 6 ! number of solid variables

              real(kind=8) :: cons(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: prim(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: xflux(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: yflux(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: xflux2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: yflux2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)

              real(kind=8) :: cons2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)

              real(kind=8) :: augmented_prim(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim+nvar2+3)

              real(kind=8) :: prim_intxl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: prim_intxr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: cons_intxl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: cons_intxr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: prim_intyl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: prim_intyr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: cons_intyl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: cons_intyr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)

              real(kind=8) :: prim_intxl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim_intxr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intxl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intxr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim_intyl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim_intyr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intyl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intyr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8), dimension(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost) :: Jd, Jk, Jdk
             
              real(kind=8) :: d_ip12, d_im12, d_im32, eps1, kappa, phi1, phi2, d1, d2
              real(kind=8) :: spl, spr, wavespeed, terml, termr, termrl, etotl, etotr

              real(kind=8) :: gam, vel, cc, ke, xc, yc, rc, etot       
              real(kind=8) :: rhocut 
              real(kind=8) :: temp, molwt(nspec), molmix
              real(kind=8) :: Rgas
              real(kind=8) :: rhomin, rhomax, rhomid, pres
              integer :: iter
              real(kind=8) :: dist
              character ( len = 100 ) filename
              real(kind=8) :: gravity, source(1:ncons), omegadot(1:nspec)
              real(kind=8) :: maxvel, maxrho, maxp, minvel, minrho, minp
              real(kind=8) :: hfk(1:nspec), enth, ie, conc(1:nspec)   

              integer :: rkstep 
              integer, parameter :: rk = 3
              real(kind=8) :: k_rk(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons,1:rk)
              real(kind=8) :: cons0(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
            
              real(kind=8) :: k_rk2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2,1:rk)
              real(kind=8) :: cons02(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)

              real(kind=8) :: dummy, T1, hfs
 
                                                      
              integer, parameter :: space_scheme = 2 ! enter 1/2/3  

              ! boundary conditions (1: periodic; 2: outflow; 3: slip wall)
              integer, parameter :: xbc = 2
              integer, parameter :: ybc = 3

              !--------------------------------------
              ! Variables for paraview 
              integer, parameter :: output_unit = 29
              character (len=100) title
              real(kind=8) :: xyz(1:imax,1:jmax,1:3)   
              !--------------------------------------

              
              real(kind=8) :: dia, rho2cut, rhosol
              real(kind=8) :: num, dia1, cd, nuss, sigmadot, T2, reyn, mu 
              real(kind=8) :: fdrag(2), qheat, cpmix, cpk(1:nspec), lambda, prandtl  
              real(kind=8) :: Jkin, Jdiff, Zhyb, Ea_R, taup, muc, source_fign, tign
              real(kind=8) :: sumyk, yk(1:nspec), mach 
 
                 call thermo_data()


              ! conservative variables
              ! 1: rho
              ! 2: rho*u
              ! 3: rho*v 
              ! 4: rho*(e + u^2/2 + v^2/2)
              ! 5: rho*Y1
              ! 6: rho*Y2

              ! primitive variables
              ! 1: rho
              ! 2: u
              ! 3: v
              ! 4: e
              ! 5: p
              ! 6: Y1
              ! 7: Y2

              !----------------------

              ! Set values here
              xlen = 1.0d0    
              ylen = 0.125d0    
              cfl = 0.2d0 
              rhocut = 1.0d-4 ! cutoff rho
              temp = 293.0d0 ! temperature
              gravity = 0.0d0 !-1.0d0 !-9.8d0 ! gravity

              molwt = MW

              !----------------------
              ! MUSCL PARAMETERS
              
               eps1 = 1.0d0 
               kappa = 1.0d0/3.0d0 
              !----------------------

              ! SOLID PARAMETERS

              dia = 1.0d-6
              rho2cut = 1.0d-4
              rhosol = 2700.0d0               
              prandtl = 0.71d0  

              Zhyb = 7.5d4
              Ea_R = 8500.0d0 

              !----------------------

 
              dx = xlen/dble(imax-1)
              dy = ylen/dble(jmax-1)
              print*, 'dx, dy: ', dx, dy 
              

              do i = 1, imax
               x(i) = dx*dble(i-1)
              enddo
              do j = 1, jmax
               y(j) = dy*dble(j-1)
              enddo

              ! initialize arrays
              cons = 0.0d0
              prim = 0.0d0

 

               Jd = 0.0d0
               Jk = 0.0d0
               Jdk = 0.0d0
 

               open(4,file='restart.data',form='formatted')
                read(4,*) irestart
               close(4)


                   n = 0
                   time = 0.0d0

                 if(irestart.eq.0) then

             do i = 1, imax-1 !1-nghost, imax-1+nghost
              do j = 1, jmax-1 !1-nghost, jmax-1+nghost

                 xc = 0.5d0*(x(i+1)+x(i)) !- xlen/2.0d0
                 yc = 0.5d0*(y(j+1)+y(j)) !- ylen/2.0d0

                 !rc = sqrt(xc**2.0d0 + yc**2.0d0)

               
                 
               prim2(i,j,1) = 330.0d-3

               if(xc.ge.0.05d0.and.xc.le.0.1d0) then
               temp = 3500.0d0 ! SURROUNDING TEMPERATURE
               prim(i,j,5) = 25.0d5 ! pressure
               else
               temp = 300.0d0  
               prim(i,j,5) = 1.0d5 ! pressure
               endif 


               ! PERTURBATION
               rc = sqrt((xc-0.13d0)**2.0d0 + (yc-0.04d0)**2.0d0)
               if(rc.le.0.01d0) then 
                temp = 3500.0d0 ! SURROUNDING TEMPERATURE
                prim(i,j,5) = 25.0d5 ! pressure
               endif
                
              
                




               prim(i,j,2) = 0.0d0 ! x-velocity
               prim(i,j,3) = 0.0d0 ! y-velocity
               

               ! species mass fraction
               prim(i,j,6) = 0.0d0 ! Al
               prim(i,j,7) = 0.21d0 ! O2
               prim(i,j,8:9) = 0.0d0 ! AlO, Al2O3
               prim(i,j,10) = 0.79d0 ! N2


               ! mol wt of mixture
               molmix = 0.0d0
               do ns = 1, nspec
                molmix = molmix + prim(i,j,5+ns)/molwt(ns)
               enddo 
               molmix = 1.0d0/molmix
               Rgas = Runiv/molmix 


               ! density
               prim(i,j,1) = prim(i,j,5)/Rgas/temp 

               ! internal energy
               call compute_ie(prim(i,j,1),prim(i,j,5),prim(i,j,6:10),prim(i,j,4)) 
               

              
               prim2(i,j,2) = 0.0d0 ! x-velocity
               prim2(i,j,3) = 0.0d0 ! y-velocity
               if(prim2(i,j,1).gt.rho2cut) then
                prim2(i,j,4) = 2.666d5
                num = prim2(i,j,1)/rhosol*6.0d0/3.14159265d0/(dia**3.0d0)
                prim2(i,j,6) = num ! number density
               else
                prim2(i,j,4) = 0.0d0 ! e2
                prim2(i,j,6) = 0.0d0 ! number density
               endif
               prim2(i,j,5) = 0.0d0 ! fign
 
             enddo
            enddo

                   else

              write(filename,'("output/output_",I5.5,".dat")'),irestart
              open(34,file=filename,form='formatted') 
               read(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1 
                 read(34,*) prim(i,j,1), prim(i,j,2), prim(i,j,3), T1, prim(i,j,5), prim(i,j,6), & 
                   prim(i,j,7), prim(i,j,8), prim(i,j,9), prim(i,j,10), prim2(i,j,1), prim2(i,j,2), & 
                   prim2(i,j,3), prim2(i,j,4), prim2(i,j,5), prim2(i,j,6), cc 
                   

               ! internal energy
               call COMPUTE_HFK (T1, HFK) 
               enth = 0.0d0
               do ns = 1, nspec
                 enth = enth + prim(i,j,5+ns)*hfk(ns)
               enddo
               ie = enth - prim(i,j,5)/prim(i,j,1)
               prim(i,j,4) = ie

                enddo 
               enddo 
              close(34)


                   endif


                  !---------------------------------------
                  !---------------------------------------
                  ! set primitive variables in ghost cells

                      if(xbc.eq.1) then
                 ! BC in x-direction (periodic)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(imax-i,j,1:nprim)
                 prim(imax-1+i,j,1:nprim) = prim(i,j,1:nprim)

                 prim2(1-i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)
                 prim2(imax-1+i,j,1:nvar2) = prim2(i,j,1:nvar2)
                enddo
               enddo
                        else if(xbc.eq.2) then
                 ! BC in x-direction (outflow)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(1,j,1:nprim)
                 prim(imax-1+i,j,1:nprim) = prim(imax-1,j,1:nprim)

                 prim2(1-i,j,1:nvar2) = prim2(1,j,1:nvar2)
                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-1,j,1:nvar2)
                enddo
               enddo
                        else if(xbc.eq.3) then
                 ! BC in x-direction (slip wall)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(i,j,1:nprim)
                 prim(imax-1+i,j,1:nprim) = prim(imax-i,j,1:nprim)

                 prim2(1-i,j,1:nvar2) = prim2(i,j,1:nvar2)
                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)

                 ! u velocity should be opposite
                 prim(1-i,j,2) = -prim(i,j,2)
                 prim(imax-1+i,j,2) = -prim(imax-i,j,2)

                 prim2(1-i,j,2) = -prim2(i,j,2)
                 prim2(imax-1+i,j,2) = -prim2(imax-i,j,2)
                enddo
               enddo
                        endif

                !-----------------------------------------------------

                        if(ybc.eq.1) then
                 ! BC in y-direction (periodic)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,1-j,1:nprim) = prim(i,jmax-j,1:nprim)
                 prim(i,jmax-1+j,1:nprim) = prim(i,j,1:nprim)

                 prim2(i,1-j,1:nvar2) = prim2(i,jmax-j,1:nvar2)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,j,1:nvar2)
                enddo
               enddo
                        else if(ybc.eq.2) then
                 ! BC in y-direction (outflow)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,1-j,1:nprim) = prim(i,1,1:nprim)
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-1,1:nprim)

                 prim2(i,1-j,1:nvar2) = prim2(i,1,1:nvar2)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-1,1:nvar2)
                enddo
               enddo
                        else if(ybc.eq.3) then
                 ! BC in y-direction (slip wall)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,1-j,1:nprim) = prim(i,j,1:nprim)
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-j,1:nprim)

                 prim2(i,1-j,1:nvar2) = prim2(i,j,1:nvar2)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-j,1:nvar2)

                 ! v velocity should be opposite
                 prim(i,1-j,3) = -prim(i,j,3)
                 prim(i,jmax-1+j,3) = -prim(i,jmax-j,3)

                 prim2(i,1-j,3) = -prim2(i,j,3)
                 prim2(i,jmax-1+j,3) = -prim2(i,jmax-j,3)
                enddo
               enddo
                        endif
                      

                  !---------------------------------------
                  !---------------------------------------


             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
               cons(i,j,1) = prim(i,j,1)
               cons(i,j,2) = prim(i,j,1)*prim(i,j,2)
               cons(i,j,3) = prim(i,j,1)*prim(i,j,3)
               cons(i,j,4) = prim(i,j,1)*prim(i,j,4) + prim(i,j,1)*(prim(i,j,2)*prim(i,j,2) &
               + prim(i,j,3)*prim(i,j,3))/2.0d0
                 
                 do ns = 1, nspec 
               cons(i,j,4+ns) = prim(i,j,1)*prim(i,j,5+ns)
                 enddo


               cons2(i,j,1) = prim2(i,j,1)
               cons2(i,j,2) = prim2(i,j,1)*prim2(i,j,2)
               cons2(i,j,3) = prim2(i,j,1)*prim2(i,j,3)
               cons2(i,j,4) = prim2(i,j,1)*prim2(i,j,4) + prim2(i,j,1)*(prim2(i,j,2)*prim2(i,j,2) &
               + prim2(i,j,3)*prim2(i,j,3))/2.0d0

               cons2(i,j,5) = prim2(i,j,1)*prim2(i,j,5)
               cons2(i,j,6) = prim2(i,j,1)*prim2(i,j,6)

              enddo
             enddo
             

                       if(irestart.eq.0) then
              print*, 'writing vtk file for paraview '
              ! write paraview vtk file
              title = 'vtk_initial'       
              n = 0
              write(filename,'("vtk/output_",I5.5,".vtk")'),n   



             do i = 1, imax
              xyz(i,:,1) = x(i) + dx 
             enddo
             do j = 1, jmax
              xyz(:,j,2) = y(j) + dy
             enddo
             xyz(:,:,3) = 0.0d0 ! no z-direction

             augmented_prim(:,:,1:nprim) = prim(:,:,1:nprim) 
             augmented_prim(:,:,nprim+1:nprim+nvar2) = prim2(:,:,1:nvar2) 
             augmented_prim(:,:,nprim+nvar2+1) = Jd(:,:)
             augmented_prim(:,:,nprim+nvar2+2) = Jk(:,:)
             augmented_prim(:,:,nprim+nvar2+3) = Jdk(:,:)
             

              call vtk_write(output_unit,filename,title,imax,jmax,xyz,nprim+nvar2+3    & 
                    ,augmented_prim(1:imax,1:jmax,1:nprim+nvar2+3))

              n = 0
              time = 0.0d0
              write(filename,'("output/output_",I5.5,".dat")'),n
              open(34,file=filename,form='formatted') 
               write(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1 

                 call compute_sound(prim(i,j,5),prim(i,j,1),prim(i,j,6:5+nspec),cc)  
                 call compute_temperature(prim(i,j,1), prim(i,j,5), prim(i,j,6:5+nspec), T1)

                 write(34,*) prim(i,j,1), prim(i,j,2), prim(i,j,3), T1, prim(i,j,5), prim(i,j,6), & 
                   prim(i,j,7), prim(i,j,8), prim(i,j,9), prim(i,j,10), prim2(i,j,1), prim2(i,j,2), & 
                   prim2(i,j,3), prim2(i,j,4), prim2(i,j,5), prim2(i,j,6), cc 

                enddo 
               enddo 
              close(34)

                      endif

!------------------------------------------------------------
              print*, 'beginning time loop '
 
              ! begin time loop
              do n = irestart, nmax
               
               vel = 0.0d0
               do i = 1, imax-1
                do j = 1, jmax-1
                       if(prim(i,j,5).le.0.0d0.or.prim(i,j,1).le.0.0d0) then
                         print*, 'something wrong: p, rho: ', prim(i,j,5), prim(i,j,1), i, j
                         stop
                       endif
                 call compute_sound(prim(i,j,5),prim(i,j,1),prim(i,j,6:5+nspec),cc) 
                 vel = max(vel,sqrt(prim(i,j,2)**2.0d0 +prim(i,j,3)**2.0d0)+cc) 
                 vel = max(vel,sqrt(prim2(i,j,2)**2.0d0 + prim2(i,j,3)**2.0d0))
                enddo 
               enddo 
                dt = cfl*min(dx,dy)/vel
                time = time + dt
                print*, ' '
                print*, 'iteration, dt, time: ', n, dt, time
                call flush() 


                ! compute fluxes at cell interfaces
                ! here, flux(i) refers to flux(i-1/2) 

                k_rk = 0.0d0                  
                cons0 = cons ! cons0 is the conservative variables (gas) at time step n (beginning of RK)   
 
                k_rk2 = 0.0d0                  
                cons02 = cons2 ! cons02 is the conservative variables (solid) at time step n (beginning of RK)   

         do rkstep = 1, rk

                !-----------------------------------------------------
                ! apply boundary conditions 

                        if(xbc.eq.1) then
                 ! BC in x-direction (periodic)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost 
                 prim(1-i,j,1:nprim) = prim(imax-i,j,1:nprim)
                 prim(imax-1+i,j,1:nprim) = prim(i,j,1:nprim)        

                 prim2(1-i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)
                 prim2(imax-1+i,j,1:nvar2) = prim2(i,j,1:nvar2)        

                 cons(1-i,j,1:ncons) = cons(imax-i,j,1:ncons)          
                 cons(imax-1+i,j,1:ncons) = cons(i,j,1:ncons)        

                 cons2(1-i,j,1:nvar2) = cons2(imax-i,j,1:nvar2)        
                 cons2(imax-1+i,j,1:nvar2) = cons2(i,j,1:nvar2)        

                 k_rk(1-i,j,1:ncons,:) = k_rk(imax-i,j,1:ncons,:)       
                 k_rk(imax-1+i,j,1:ncons,:) = k_rk(i,j,1:ncons,:)    

                 k_rk2(1-i,j,1:nvar2,:) = k_rk2(imax-i,j,1:nvar2,:)    
                 k_rk2(imax-1+i,j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)    
                enddo
               enddo
                        else if(xbc.eq.2) then
                 ! BC in x-direction (outflow)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost 
                 prim(1-i,j,1:nprim) = prim(1,j,1:nprim)
                 prim(imax-1+i,j,1:nprim) = prim(imax-1,j,1:nprim)     

                 prim2(1-i,j,1:nvar2) = prim2(1,j,1:nvar2)
                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-1,j,1:nvar2)   

                 cons(1-i,j,1:ncons) = cons(1,j,1:ncons)          
                 cons(imax-1+i,j,1:ncons) = cons(imax-1,j,1:ncons)     

                 cons2(1-i,j,1:nvar2) = cons2(1,j,1:nvar2)          
                 cons2(imax-1+i,j,1:nvar2) = cons2(imax-1,j,1:nvar2)   

                 k_rk(1-i,j,1:ncons,:) = k_rk(1,j,1:ncons,:)       
                 k_rk(imax-1+i,j,1:ncons,:) = k_rk(imax-1,j,1:ncons,:) 

                 k_rk2(1-i,j,1:nvar2,:) = k_rk2(1,j,1:nvar2,:)       
                 k_rk2(imax-1+i,j,1:nvar2,:) = k_rk2(imax-1,j,1:nvar2,:)
                enddo
               enddo
                        else if(xbc.eq.3) then
                 ! BC in x-direction (slip wall)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(i,j,1:nprim)
                 prim(imax-1+i,j,1:nprim) = prim(imax-i,j,1:nprim)

                 prim2(1-i,j,1:nvar2) = prim2(i,j,1:nvar2)
                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)

                 cons(1-i,j,1:ncons) = cons(i,j,1:ncons)          
                 cons(imax-1+i,j,1:ncons) = cons(imax-i,j,1:ncons)     

                 cons2(1-i,j,1:nvar2) = cons2(i,j,1:nvar2)          
                 cons2(imax-1+i,j,1:nvar2) = cons2(imax-i,j,1:nvar2)   

                 k_rk(1-i,j,1:ncons,:) = k_rk(i,j,1:ncons,:)       
                 k_rk(imax-1+i,j,1:ncons,:) = k_rk(imax-i,j,1:ncons,:) 

                 k_rk2(1-i,j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)       
                 k_rk2(imax-1+i,j,1:nvar2,:) = k_rk2(imax-i,j,1:nvar2,:)

                 ! u velocity should be opposite
                 prim(1-i,j,2) = -prim(i,j,2)
                 prim(imax-1+i,j,2) = -prim(imax-i,j,2)

                 prim2(1-i,j,2) = -prim2(i,j,2)
                 prim2(imax-1+i,j,2) = -prim2(imax-i,j,2)

                 cons(1-i,j,2) = -cons(i,j,2)          
                 cons(imax-1+i,j,2) = -cons(imax-i,j,2)     

                 cons2(1-i,j,2) = -cons2(i,j,2)          
                 cons2(imax-1+i,j,2) = -cons2(imax-i,j,2)     

                 k_rk(1-i,j,2,:) = -k_rk(i,j,2,:)       
                 k_rk(imax-1+i,j,2,:) = -k_rk(imax-i,j,2,:) 

                 k_rk2(1-i,j,2,:) = -k_rk2(i,j,2,:)       
                 k_rk2(imax-1+i,j,2,:) = -k_rk2(imax-i,j,2,:) 

                enddo
               enddo
                        endif

                !-----------------------------------------------------

                        if(ybc.eq.1) then
                 ! BC in y-direction (periodic)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,1-j,1:nprim) = prim(i,jmax-j,1:nprim)
                 prim(i,jmax-1+j,1:nprim) = prim(i,j,1:nprim)        

                 prim2(i,1-j,1:nvar2) = prim2(i,jmax-j,1:nvar2)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,j,1:nvar2)        

                 cons(i,1-j,1:ncons) = cons(i,jmax-j,1:ncons)          
                 cons(i,jmax-1+j,1:ncons) = cons(i,j,1:ncons)        

                 cons2(i,1-j,1:nvar2) = cons2(i,jmax-j,1:nvar2)        
                 cons2(i,jmax-1+j,1:nvar2) = cons2(i,j,1:nvar2)        

                 k_rk(i,1-j,1:ncons,:) = k_rk(i,jmax-j,1:ncons,:)       
                 k_rk(i,jmax-1+j,1:ncons,:) = k_rk(i,j,1:ncons,:)    

                 k_rk2(i,1-j,1:nvar2,:) = k_rk2(i,jmax-j,1:nvar2,:)     
                 k_rk2(i,jmax-1+j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)    
                enddo
               enddo
                        else if(ybc.eq.2) then
                 ! BC in y-direction (outflow)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,1-j,1:nprim) = prim(i,1,1:nprim)
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-1,1:nprim)     

                 prim2(i,1-j,1:nvar2) = prim2(i,1,1:nvar2)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-1,1:nvar2)    

                 cons(i,1-j,1:ncons) = cons(i,1,1:ncons)          
                 cons(i,jmax-1+j,1:ncons) = cons(i,jmax-1,1:ncons)     

                 cons2(i,1-j,1:nvar2) = cons2(i,1,1:nvar2)          
                 cons2(i,jmax-1+j,1:nvar2) = cons2(i,jmax-1,1:nvar2)   

                 k_rk(i,1-j,1:ncons,:) = k_rk(i,1,1:ncons,:)       
                 k_rk(i,jmax-1+j,1:ncons,:) = k_rk(i,jmax-1,1:ncons,:)  

                 k_rk2(i,1-j,1:nvar2,:) = k_rk2(i,1,1:nvar2,:)       
                 k_rk2(i,jmax-1+j,1:nvar2,:) = k_rk2(i,jmax-1,1:nvar2,:)
                enddo
               enddo
                        else if(ybc.eq.3) then
                 ! BC in y-direction (slip wall)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,1-j,1:nprim) = prim(i,j,1:nprim)
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-j,1:nprim)     

                 prim2(i,1-j,1:nvar2) = prim2(i,j,1:nvar2)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-j,1:nvar2)     

                 cons(i,1-j,1:ncons) = cons(i,j,1:ncons)          
                 cons(i,jmax-1+j,1:ncons) = cons(i,jmax-j,1:ncons)     

                 cons2(i,1-j,1:nvar2) = cons2(i,j,1:nvar2)          
                 cons2(i,jmax-1+j,1:nvar2) = cons2(i,jmax-j,1:nvar2)     

                 k_rk(i,1-j,1:ncons,:) = k_rk(i,j,1:ncons,:)       
                 k_rk(i,jmax-1+j,1:ncons,:) = k_rk(i,jmax-j,1:ncons,:)  

                 k_rk2(i,1-j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)       
                 k_rk2(i,jmax-1+j,1:nvar2,:) = k_rk2(i,jmax-j,1:nvar2,:)  

                 ! v velocity should be opposite
                 prim(i,1-j,3) = -prim(i,j,3)
                 prim(i,jmax-1+j,3) = -prim(i,jmax-j,3)     

                 prim2(i,1-j,3) = -prim2(i,j,3)
                 prim2(i,jmax-1+j,3) = -prim2(i,jmax-j,3)     

                 cons(i,1-j,3) = -cons(i,j,3)          
                 cons(i,jmax-1+j,3) = -cons(i,jmax-j,3)     

                 cons2(i,1-j,3) = -cons2(i,j,3)          
                 cons2(i,jmax-1+j,3) = -cons2(i,jmax-j,3)     

                 k_rk(i,1-j,3,:) = -k_rk(i,j,3,:)       
                 k_rk(i,jmax-1+j,3,:) = -k_rk(i,jmax-j,3,:)  

                 k_rk2(i,1-j,3,:) = -k_rk2(i,j,3,:)       
                 k_rk2(i,jmax-1+j,3,:) = -k_rk2(i,jmax-j,3,:)  

                enddo
               enddo
                        endif


                !-----------------------------------------------------

                  !print*, 'interpol x '

                ! Interpolate in x-direction 
                do i = 0, imax+1 
                 do j = 0, jmax+1 

                 ! compute primitive and conservative variables at the cell interfaces
                 ! i refers to i-1/2 for interface variables


             if(space_scheme.eq.1) then

               prim_intxl(i+1,j,:) = prim(i,j,:)
               prim_intxr(i,j,:) = prim(i,j,:)

               prim_intxl2(i+1,j,:) = prim2(i,j,:)
               prim_intxr2(i,j,:) = prim2(i,j,:)

             else if(space_scheme.eq.2) then

               
              
              call weno5_gas(nprim,prim(i-2:i+2,j,1:nprim),prim_intxl(i+1,j,1:nprim),prim_intxr(i,j,1:nprim))

              call weno5_solid(nvar2,prim2(i-2:i+2,j,1:nvar2),prim_intxl2(i+1,j,1:nvar2),prim_intxr2(i,j,1:nvar2))


             else if(space_scheme.eq.3) then

               do ns = 1, ncons             
               d_im12 = cons(i,j,ns) - cons(i-1,j,ns)
               d_ip12 = cons(i+1,j,ns) - cons(i,j,ns)

               cons_intxl(i+1,j,ns) = cons(i,j,ns)  +&
                  eps1/4.0d0*((1.0d0-kappa)*d_im12+(1.0d0+kappa)*d_ip12)
               cons_intxr(i,j,ns) = cons(i,j,ns) -&
                  eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo

              do ns = 1, nvar2
               d_im12 = cons2(i,j,ns) - cons2(i-1,j,ns)
               d_ip12 = cons2(i+1,j,ns) - cons2(i,j,ns)

               cons_intxl2(i+1,j,ns) = cons2(i,j,ns) +&
                   eps1/4.0d0*((1.0d0-kappa)*d_im12+(1.0d0+kappa)*d_ip12)
               cons_intxr2(i,j,ns) = cons2(i,j,ns) -&
                   eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo

                     
              call cons_2_prim(rho2cut,ncons,nvar2,nspec,cons_intxl(i+1,j,:), &
                   cons_intxl2(i+1,j,:),prim_intxl(i+1,j,:),prim_intxl2(i+1,j,:))
              call cons_2_prim(rho2cut,ncons,nvar2,nspec,cons_intxr(i,j,:), &
                   cons_intxr2(i,j,:),prim_intxr(i,j,:),prim_intxr2(i,j,:))


             else
               
               print*, 'wrong entry for space_scheme ', space_scheme

             endif

              

                

              ! internal energy
              call compute_ie(prim_intxl(i+1,j,1),prim_intxl(i+1,j,5),prim_intxl(i+1,j,6:10),prim_intxl(i+1,j,4)) 
              call compute_ie(prim_intxr(i,j,1),prim_intxr(i,j,5),prim_intxr(i,j,6:10),prim_intxr(i,j,4))             

                 enddo
                enddo  

                !-----------------------------------------------------

                   !print*, 'interpol y'

                ! Interpolate in y-direction 
                do j = 0, jmax+1 
                 do i = 0, imax+1 

                 ! compute primitive and conservative variables at the cell interfaces
                 ! i refers to i-1/2 for interface variables


                  if(space_scheme.eq.1) then

               prim_intyl(i,j+1,:) = prim(i,j,:)
               prim_intyr(i,j,:) = prim(i,j,:)

               prim_intyl2(i,j+1,:) = prim2(i,j,:)
               prim_intyr2(i,j,:) = prim2(i,j,:)

             else if(space_scheme.eq.2) then

               
         call weno5_gas(nprim,prim(i,j-2:j+2,1:nprim),prim_intyl(i,j+1,1:nprim),prim_intyr(i,j,1:nprim))

         call weno5_solid(nvar2,prim2(i,j-2:j+2,1:nvar2),prim_intyl2(i,j+1,1:nvar2),prim_intyr2(i,j,1:nvar2))



              else if(space_scheme.eq.3) then

               do ns = 1, ncons
               d_im12 = cons(i,j,ns) - cons(i,j-1,ns)
               d_ip12 = cons(i,j+1,ns) - cons(i,j,ns)

               cons_intyl(i,j+1,ns) = cons(i,j,ns)   &
                + eps1/4.0d0*((1.0d0-kappa)*d_im12+(1.0d0+kappa)*d_ip12)
               cons_intyr(i,j,ns) = cons(i,j,ns)     &
                - eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo

              do ns = 1, nvar2
               d_im12 = cons2(i,j,ns) - cons2(i,j-1,ns)
               d_ip12 = cons2(i,j+1,ns) - cons2(i,j,ns)

               cons_intyl2(i,j+1,ns) = cons2(i,j,ns)   &
                + eps1/4.0d0*((1.0d0-kappa)*d_im12+(1.0d0+kappa)*d_ip12)
               cons_intyr2(i,j,ns) = cons2(i,j,ns)     &
                - eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo

              call cons_2_prim(rho2cut,ncons,nvar2,nspec,cons_intyl(i,j+1,:), &
                   cons_intyl2(i,j+1,:),prim_intyl(i,j+1,:),prim_intyl2(i,j+1,:))
              call cons_2_prim(rho2cut,ncons,nvar2,nspec,cons_intyr(i,j,:), &
                   cons_intyr2(i,j,:),prim_intyr(i,j,:),prim_intyr2(i,j,:))


             else
               
               print*, 'wrong entry for space_scheme ', space_scheme

             endif


              



               ! internal energy
               call compute_ie(prim_intyl(i,j+1,1),prim_intyl(i,j+1,5),prim_intyl(i,j+1,6:10),prim_intyl(i,j+1,4)) 
               call compute_ie(prim_intyr(i,j,1),prim_intyr(i,j,5),prim_intyr(i,j,6:10),prim_intyr(i,j,4))             
   
                 enddo
                enddo  

                !-----------------------------------------------------
          
                    !print*, 'flux computation '

                          


 
           xflux = 0.0d0
           yflux = 0.0d0
           xflux2 = 0.0d0
           yflux2 = 0.0d0

 
                 ! compute x-fluxes
           do i = 1, imax !0, imax
            do j = 1, jmax !0, jmax



                 
                call hllc_2d(prim_intxl(i,j,1:nprim),prim_intxr(i,j,1:nprim),1,xflux(i,j,1:ncons))
                !call rusanov_gas_2d(prim_intxl(i,j,1:nprim),prim_intxr(i,j,1:nprim),1,xflux(i,j,1:ncons))
                !call castro_2d(prim_intxl(i,j,1:nprim),prim_intxr(i,j,1:nprim),1,xflux(i,j,1:ncons))

                 
                            ! SOLID 

              spl = abs(prim_intxl2(i,j,2))  
              spr = abs(prim_intxr2(i,j,2))  
              wavespeed = max(spl,spr)
         
                 ! continuity
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)
                 termrl = prim_intxr2(i,j,1) - prim_intxl2(i,j,1)
                 xflux2(i,j,1) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! x-momentum
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)*prim_intxl2(i,j,2)  
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)*prim_intxr2(i,j,2)  
                 termrl = prim_intxr2(i,j,1)*prim_intxr2(i,j,2) - prim_intxl2(i,j,1)*prim_intxl2(i,j,2)
                 xflux2(i,j,2) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! y-momentum
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)*prim_intxl2(i,j,3)  
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)*prim_intxr2(i,j,3)  
                 termrl = prim_intxr2(i,j,1)*prim_intxr2(i,j,3) - prim_intxl2(i,j,1)*prim_intxl2(i,j,3)
                 xflux2(i,j,3) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! energy
                 ke = 0.5d0*prim_intxl2(i,j,1)*(prim_intxl2(i,j,2)**2.0d0 + prim_intxl2(i,j,3)**2.0d0)
                 etotl = prim_intxl2(i,j,1)*prim_intxl2(i,j,4) + ke
                 terml = etotl*prim_intxl2(i,j,2)
                 ke = 0.5d0*prim_intxr2(i,j,1)*(prim_intxr2(i,j,2)**2.0d0 + prim_intxr2(i,j,3)**2.0d0)
                 etotr = prim_intxr2(i,j,1)*prim_intxr2(i,j,4) + ke
                 termr = etotr*prim_intxr2(i,j,2)
                 termrl = etotr - etotl
                 xflux2(i,j,4) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl 

                 ! fign
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)*prim_intxl2(i,j,5)  
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)*prim_intxr2(i,j,5)  
                 termrl = prim_intxr2(i,j,1)*prim_intxr2(i,j,5) - prim_intxl2(i,j,1)*prim_intxl2(i,j,5)
                 xflux2(i,j,5) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! Num
                 terml = prim_intxl2(i,j,2)*prim_intxl2(i,j,6)  
                 termr = prim_intxr2(i,j,2)*prim_intxr2(i,j,6)  
                 termrl = prim_intxr2(i,j,6) - prim_intxl2(i,j,6)
                 xflux2(i,j,6) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl


            enddo
           enddo

                          !--------!


                 ! compute y-fluxes
           do i = 1, imax !0, imax
            do j = 1, jmax !0, jmax


                
                 call hllc_2d(prim_intyl(i,j,1:nprim),prim_intyr(i,j,1:nprim),2,yflux(i,j,1:ncons))
                 !call rusanov_gas_2d(prim_intyl(i,j,1:nprim),prim_intyr(i,j,1:nprim),2,yflux(i,j,1:ncons))
                 !call castro_2d(prim_intyl(i,j,1:nprim),prim_intyr(i,j,1:nprim),2,yflux(i,j,1:ncons))


                            ! SOLID


              spl = abs(prim_intyl2(i,j,3))
              spr = abs(prim_intyr2(i,j,3)) 
              wavespeed = max(spl,spr)

                 ! continuity
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)
                 termrl = prim_intyr2(i,j,1) - prim_intyl2(i,j,1)
                 yflux2(i,j,1) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! x-momentum
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)*prim_intyl2(i,j,2)  
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)*prim_intyr2(i,j,2)  
                 termrl = prim_intyr2(i,j,1)*prim_intyr2(i,j,2) - prim_intyl2(i,j,1)*prim_intyl2(i,j,2)
                 yflux2(i,j,2) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! y-momentum
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)*prim_intyl2(i,j,3)  
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)*prim_intyr2(i,j,3)  
                 termrl = prim_intyr2(i,j,1)*prim_intyr2(i,j,3) - prim_intyl2(i,j,1)*prim_intyl2(i,j,3)
                 yflux2(i,j,3) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! energy
                 ke = 0.5d0*prim_intyl2(i,j,1)*(prim_intyl2(i,j,2)**2.0d0 + prim_intyl2(i,j,3)**2.0d0)
                 etotl = prim_intyl2(i,j,1)*prim_intyl2(i,j,4) + ke
                 terml = etotl*prim_intyl2(i,j,3)
                 ke = 0.5d0*prim_intyr2(i,j,1)*(prim_intyr2(i,j,2)**2.0d0 + prim_intyr2(i,j,3)**2.0d0)
                 etotr = prim_intyr2(i,j,1)*prim_intyr2(i,j,4) + ke
                 termr = etotr*prim_intyr2(i,j,3)
                 termrl = etotr - etotl
                 yflux2(i,j,4) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! fign
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)*prim_intyl2(i,j,5)  
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)*prim_intyr2(i,j,5)  
                 termrl = prim_intyr2(i,j,1)*prim_intyr2(i,j,5) - prim_intyl2(i,j,1)*prim_intyl2(i,j,5)
                 yflux2(i,j,5) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl 
                
                 ! Num
                 terml = prim_intyl2(i,j,3)*prim_intyl2(i,j,6)  
                 termr = prim_intyr2(i,j,3)*prim_intyr2(i,j,6)  
                 termrl = prim_intyr2(i,j,6) - prim_intyl2(i,j,6)
                 yflux2(i,j,6) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl 
                

            enddo
           enddo


                !-----------------------------------------------------



                ! Runge-Kutta increment
                do i = 1, imax-1 
                 do j = 1, jmax-1 
             k_rk(i,j,1:ncons,rkstep) = - (xflux(i+1,j,:)-xflux(i,j,:))*dt/dx   &
                     - (yflux(i,j+1,:)-yflux(i,j,:))*dt/dy
             k_rk2(i,j,1:nvar2,rkstep) = - (xflux2(i+1,j,:)-xflux2(i,j,:))*dt/dx   &
                     - (yflux2(i,j+1,:)-yflux2(i,j,:))*dt/dy
                 enddo
                enddo


                ! update the conserved variables
                ! assume no source terms for now
            do i = 1, imax-1 
             do j = 1, jmax-1 
              if(rkstep.eq.1) then
                   cons(i,j,:) = cons0(i,j,:) + k_rk(i,j,:,1)/2.0d0
                   cons2(i,j,:) = cons02(i,j,:) + k_rk2(i,j,:,1)/2.0d0
              endif 
              if(rkstep.eq.2) then
                   cons(i,j,:) = cons0(i,j,:) - k_rk(i,j,:,1) + 2.0d0*k_rk(i,j,:,2)
                   cons2(i,j,:) = cons02(i,j,:) - k_rk2(i,j,:,1) + 2.0d0*k_rk2(i,j,:,2)
              endif 
             enddo
            enddo

               if(rkstep.eq.rk) then
                do i = 1, imax-1
                 do j = 1, jmax-1
         cons(i,j,:) = cons0(i,j,:) + 1.0d0/6.0d0*(1.0d0*k_rk(i,j,:,1)    &
         + 4.0d0*k_rk(i,j,:,2) + 1.0d0*k_rk(i,j,:,3))
         cons2(i,j,:) = cons02(i,j,:) + 1.0d0/6.0d0*(1.0d0*k_rk2(i,j,:,1)    &
         + 4.0d0*k_rk2(i,j,:,2) + 1.0d0*k_rk2(i,j,:,3))
                 enddo
                enddo
               endif

     
!-------------------------------------------------------------------
           
                ! recompute primitive variables
                do i = 1, imax-1
                 do j = 1, jmax-1

                  if(cons(i,j,1).le.rhocut) then
                    print*, 'too small a rho ', i,j, cons(i,j,:),    &
         xflux(i+1,j,1), xflux(i,j,1), yflux(i,j+1,1), yflux(i,j,1)
                    stop
                  endif

        call cons_2_prim(rho2cut,ncons,nvar2,nspec,cons(i,j,:),cons2(i,j,:),prim(i,j,:),prim2(i,j,:))

                 enddo ! do j             
                enddo ! do i             

 
!---------------------------------------------------              



         ! close rk loop
         enddo   


                          

!----------------------------------------------------------------
               ! add source terms (gravity)
             

                  !print*, 'source terms '

              
             !do i = 1, imax-1 !1-nghost, imax-1+nghost
             ! do j = 1, jmax-1 !1-nghost, jmax-1+nghost
             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
               source(1:ncons) = 0.0d0
               source(3) = gravity*dt*cons(i,j,1)
               source(4) = gravity*dt*cons(i,j,3)

               cons(i,j,1:ncons) = cons(i,j,1:ncons) + source(1:ncons)
              enddo
             enddo


                       ! chemistry

             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
               source(1:ncons) = 0.0d0

                  if(cons(i,j,5).gt.0.0d0.and.cons(i,j,6).gt.0.0d0) then
               molmix = 0.0d0
               do ns = 1, nspec
                molmix = molmix + prim(i,j,5+ns)/molwt(ns)
               enddo 
               molmix = 1.0d0/molmix
               Rgas = Runiv/molmix 

                     temp = prim(i,j,5)/Rgas/prim(i,j,1)

              do ns = 1, nspec
               !conc(ns) = (cons(i,j,4+ns)/1000.0d0)/MW(ns)
               conc(ns) = cons(i,j,4+ns)
              enddo


              omegadot = 0.0d0
              call combustion(temp,conc,omegadot)
              omegadot = omegadot/dt

              do ns = 1, nspec
               source(4+ns) = omegadot(ns)
              enddo

                  endif

               cons(i,j,1:ncons) = cons(i,j,1:ncons) + source(1:ncons)*dt

                     do ns = 1, nspec
                      cons(i,j,4+ns) = max(cons(i,j,4+ns),0.0d0)
                     enddo


              enddo
             enddo
              

 
                     Jd = 0.0d0 
                     Jk = 0.0d0
                     Jdk = 0.0d0
 

 
                        ! TWO-PHASE SOURCE TERMS
             do i = 1, imax-1
              do j = 1, jmax-1
                source(1:ncons) = 0.0d0
                source_fign = 0.0d0
                hfs = 0.0d0
                sigmadot = 0.0d0

                if(prim2(i,j,1).gt.rho2cut.and.prim2(i,j,6).gt.0.0d0) then

               ! find gas temperature 
               molmix = 0.0d0
               do ns = 1, nspec
                molmix = molmix + prim(i,j,5+ns)/molwt(ns)
               enddo 
               molmix = 1.0d0/molmix
               Rgas = Runiv/molmix 
               T1 = prim(i,j,5)/Rgas/prim(i,j,1)
                if(T1.le.100.0d0) then
                  print*, 'T1 very low ', T1
                  stop
                endif
             

                 ! find particle diameter
                 num = prim2(i,j,6)
                 dia1 = (6.0d0/3.14159265d0*prim2(i,j,1)/rhosol/num)**(1.0d0/3.0d0)  
                 dia1 = max(min(dia1,dia),0.1d-6) 
 

                 ! find reynolds number
                 vel = sqrt((prim(i,j,2)-prim2(i,j,2))**2.0d0 + (prim(i,j,3)-prim2(i,j,3))**2.0d0)
                 mu = 4.0d-5*sqrt(T1/300.0d0) 
                 reyn = dia1*vel*prim(i,j,1)/mu   

                 ! find Mach number
                 call compute_sound(prim(i,j,5),prim(i,j,1),prim(i,j,6:5+nspec),cc) 
                 mach = vel/cc

                 ! find drag coeff.
                 call drag(reyn,mach,cd)      
                 

 
                 ! find nusselt number 
                 nuss = 2.0d0 + 0.6d0*(prandtl**0.33d0)*sqrt(reyn)    
                 !nuss = min(nuss,500.0d0)


                 ! find cpmix
                 call COMPUTE_CPK (T1, CPK)
                 cpmix = 0.0d0
                 do ns = 1, nspeci
                  cpmix = cpmix + prim(i,j,5+ns)*cpk(ns)
                 enddo

                 ! find lambda
                 lambda = mu*cpmix/prandtl

                 ! T2
                 call from_es_to_Ts(prim2(i,j,4),T2)
                 


                          

                 ! find fdrag and qheat
                 fdrag(1) = 3.0d0/4.0d0*prim(i,j,1)/rhosol*prim2(i,j,1)/dia1*cd*vel*(prim(i,j,2)-prim2(i,j,2)) 
                 fdrag(2) = 3.0d0/4.0d0*prim(i,j,1)/rhosol*prim2(i,j,1)/dia1*cd*vel*(prim(i,j,3)-prim2(i,j,3)) 
                 qheat = 6.0d0*prim2(i,j,1)/rhosol/dia1*nuss*lambda*(T1-T2)/dia1


                               
                               

                 Tign = 933.0d0
                 ! mass transfer
                 
                 if(prim2(i,j,1).gt.0.0d0.and.prim2(i,j,6).ge.0.0d0.and.prim(i,j,7).gt.0.0d0.and.T1.ge.500.0d0) then 
                  Jkin = 6.0d0*prim2(i,j,1)/rhosol/dia1*Zhyb*exp(-Ea_R/T2) 
                  taup = 1.5d6*(dia*dia)/(prim(i,j,7)**0.9d0) 
                  Jdiff = 3.0d0*prim2(i,j,1)*(1.0d0 + 0.276d0*sqrt(reyn))/taup
                  muc = 1.0d0/(1.0d0 + exp((130.0d-3-prim2(i,j,1))/(20.0d-3)))

                  sigmadot = 1.0d0/(1.0d0/Jkin + 1.0d0/Jdiff)*muc

                  Jd(i,j) = Jdiff
                  Jk(i,j) = Jkin  

                 else
                  sigmadot = 0.0d0  
                  Jd(i,j) = 0.0d0
                  Jk(i,j) = 0.0d0                  
                 endif             


                   Jdk(i,j) = Jd(i,j)/(Jd(i,j)+Jk(i,j)+1.0d-12)
                   if(Jdk(i,j).gt.1.0d0) then
                    print*, 'Jdk > 1 ', i, j, Jd(i,j), Jk(i,j), Jdk(i,j)
                    stop
                   endif


                  if(sigmadot.lt.0.0d0) then
                   print*, 'bug in sigmadot: ', sigmadot, Jkin, Jdiff
                   call flush()
                   stop
                  endif

                  sigmadot = min(sigmadot,cons2(i,j,1)/dt)
                  sigmadot = max(sigmadot,0.0d0)



                  call HTOTAL(T2,hfs)

                              
                              

                             
                  source(1) = sigmadot
                  source(2) = sigmadot*prim2(i,j,2) - fdrag(1)
                  source(3) = sigmadot*prim2(i,j,3) - fdrag(2)
                  source(4) = sigmadot*(prim2(i,j,4) + 0.5d0*(prim2(i,j,2)**2.0d0 + prim2(i,j,3)**2.0d0)) &
                               - fdrag(1)*prim2(i,j,2) - fdrag(2)*prim2(i,j,3) - qheat 
                  source(5) = sigmadot
                  source(6:4+nspec) = 0.0d0

                  ! KUHL-BOIKO CURVE FIT
                  tign = 1.0d0/(6.25d10)*exp(30000/T1) 
                  source_fign = prim2(i,j,1)/tign

                else
                 source = 0.0d0 
                 source_fign = 0.0d0
                endif


 
                                
                   if(prim2(i,j,1).le.rho2cut.or.prim2(i,j,6).le.0.0d0) then
                       source = 0.0d0
                   endif
                           


               cons(i,j,1:ncons) = cons(i,j,1:ncons) + source(1:ncons)*dt               
               cons(i,j,4) = cons(i,j,4) 
                                 

               cons2(i,j,1:4) = cons2(i,j,1:4) - source(1:4)*dt  
               cons2(i,j,5) = cons2(i,j,5) + source_fign*dt 

               if(cons2(i,j,1).le.rho2cut) then
                 prim2(i,j,1:nvar2) = 0.0d0 
                 cons2(i,j,1:nvar2) = 0.0d0
               endif

              enddo  
             enddo



                       
!----------------------------------------------------------------

                    ! upper bound on fign
                 
             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
                cons2(i,j,5) = min(cons2(i,j,5),1.5d0*cons2(i,j,1))
              enddo
             enddo

!----------------------------------------------------------------

                    ! recompute primitive variables
             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost

                  if(cons(i,j,1).le.rhocut) then
                    print*, 'too small a rho ', i,j, cons(i,j,:),    &
         xflux(i+1,j,1), xflux(i,j,1), yflux(i,j+1,1), yflux(i,j,1)
                    stop
                  endif

         call cons_2_prim(rho2cut,ncons,nvar2,nspec,cons(i,j,:),cons2(i,j,:),prim(i,j,:),prim2(i,j,:))

             enddo ! do j             
            enddo ! do i




               

!----------------------------------------------------------------
                  maxvel = 0.0d0
                  maxp = 0.0d0
                  maxrho = 0.0d0
                  minvel = 1.0d10
                  minp = 1.0d10
                  minrho = 1.0d10

                  do i = 1, imax-1 
                   do j = 1, jmax-1
                      maxp = max(maxp,prim(i,j,5))
                      maxrho = max(maxrho,prim(i,j,1))
                      maxvel = max(maxvel,abs(prim(i,j,2)),abs(prim(i,j,3)))

                      minp = min(minp,prim(i,j,5))
                      minrho = min(minrho,prim(i,j,1))
                      minvel = min(minvel,abs(prim(i,j,2)),abs(prim(i,j,3)))
                   enddo
                  enddo  
                 print*, 'max (pressure,density,velocity): ', maxp, maxrho, maxvel 
                 print*, 'min (pressure,density,velocity): ', minp, minrho, minvel 

             augmented_prim(:,:,1:nprim) = prim(:,:,1:nprim) 
             augmented_prim(:,:,nprim+1:nprim+nvar2) = prim2(:,:,1:nvar2) 
             augmented_prim(:,:,nprim+nvar2+1) = Jd(:,:)
             augmented_prim(:,:,nprim+nvar2+2) = Jk(:,:)
             augmented_prim(:,:,nprim+nvar2+3) = Jdk(:,:)


                 ! write paraview vtk file 
                 if(mod(n,ofile).eq.0) then

                         


              print*, ' '
              print*, 'writing vtk file for paraview '
              title = 'vtk/vtk_output'       

             do i = 1, imax
              xyz(i,:,1) = x(i) + dx 
             enddo
             do j = 1, jmax
              xyz(:,j,2) = y(j) + dy
             enddo
             xyz(:,:,3) = 0.0d0 ! no z-direction

              write(filename,'("vtk/output_",I5.5,".vtk")'),n   
              call vtk_write(output_unit,filename,title,imax,jmax,xyz,nprim+nvar2    & 
                       ,augmented_prim(1:imax,1:jmax,1:nprim+nvar2+3))

                   print*, 'writing output file '

              write(filename,'("output/output_",I5.5,".dat")'),n
              open(34,file=filename,form='formatted') 
               write(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1 

                 call compute_temperature(prim(i,j,1), prim(i,j,5), prim(i,j,6:5+nspec), T1)
                 call compute_sound(prim(i,j,5),prim(i,j,1),prim(i,j,6:5+nspec),cc) 

                 write(34,*) prim(i,j,1), prim(i,j,2), prim(i,j,3), T1, prim(i,j,5), prim(i,j,6), & 
                   prim(i,j,7), prim(i,j,8), prim(i,j,9), prim(i,j,10), prim2(i,j,1), prim2(i,j,2), & 
                   prim2(i,j,3), prim2(i,j,4), prim2(i,j,5), prim2(i,j,6), cc
 
                enddo 
               enddo 
              close(34)


              write(filename,'("combustion/J_",I5.5,".dat")'),n
              open(34,file=filename,form='formatted') 
               write(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1                  

                 write(34,*) Jd(i,j), Jk(i,j), Jdk(i,j)
                enddo 
               enddo 
              close(34)
              


                 endif 


                ! close time loop (do n)
              enddo
              

 
              print*, 'done ' 


                               end program

!--------------------------------------------------------------

         subroutine compute_sound(p,rho,Yk,c)

         use thermo

          real(kind=8) :: p, rho, Yk(1:nspeci), c, temp, molwt, Rgas
          real(kind=8) :: cpk(1:nspeci), cvk(1:nspeci), gam
          integer :: ns
          real(kind=8) :: cpmix, cvmix

          molwt = 0.0d0
          do ns = 1, nspeci
           molwt = molwt + Yk(ns)/MW(ns) 
          enddo  
          molwt = 1.0d0/molwt
          Rgas = runiv/molwt
          temp = p/rho/Rgas 


          call COMPUTE_CPK (TEMP, CPK)
          do ns = 1, nspeci
           cvk(ns) = cpk(ns) - runiv/MW(ns)  
          enddo 
          
          cpmix = 0.0d0
          cvmix = 0.0d0
          do ns = 1, nspeci
           cpmix = cpmix + Yk(ns)*cpk(ns)
           cvmix = cvmix + Yk(ns)*cvk(ns)
          enddo
           
          gam = cpmix/cvmix


          if(p.le.0.0d0.or.rho.le.0.0d0) then
            print*, 'sound speed not physical ', p, rho, Yk 
            stop
          endif
 
          c = sqrt(gam*p/rho)


!          if(c.lt.50.0d0.or.c.gt.1000.0d0) then
!           print*, 'something awry with c ', c, p, rho
!           stop
!          endif 

         end subroutine

!--------------------------------------------------------------

        subroutine compute_pres(rho, yk, ie, pres) 

        USE THERMO 
 
        implicit none

        real(kind=8) :: rho, yk(1:nspeci), ie, pres
        real(kind=8) :: temp, molwt, Rgas  
        integer :: ns

        call FROM_IE_TO_T(ie, YK, TEMP)
       
            if(temp.le.0.0d0) then
             print*, 'temp < 0 ', temp, yk, sum(yk), ie, rho
             stop 
            endif
 
        molwt = 0.0d0
        do ns = 1, nspeci
          molwt = molwt + yk(ns)/MW(ns)
        enddo 
        molwt = 1.0d0/molwt
        Rgas = runiv/molwt 

        pres = rho*Rgas*temp


        end subroutine

!--------------------------------------------------------------

        subroutine compute_ie(rho, pres, yk, ie) 

        USE THERMO 
 
        implicit none

        real(kind=8) :: rho, yk(1:nspeci), ie, pres
        real(kind=8) :: temp, molwt, Rgas, hfk(1:nspeci)  
        integer :: ns
        real(kind=8) :: enth  

        
        call compute_temperature(rho,pres,yk,TEMP)
        call COMPUTE_HFK (TEMP, HFK) 
        enth = 0.0d0
        do ns = 1, nspeci
         enth = enth + yk(ns)*hfk(ns)
        enddo
        ie = enth - pres/rho

               
        return
        end subroutine

!--------------------------------------------------------------

        subroutine compute_temperature(rho, pres, yk, T)

        USE THERMO

        implicit none

        real(kind=8) :: rho, pres, T, yk(1:nspeci)
        real(kind=8) :: molwt, Rgas
        integer :: ns

        molwt = 0.0d0
        do ns = 1, nspeci
          molwt = molwt + yk(ns)/MW(ns)
        enddo
        molwt = 1.0d0/molwt
        Rgas = runiv/molwt

        T= pres/rho/Rgas

           if(T.le.100.0d0) then
            print*, 'T very low ', T, rho, pres, yk 
            T = 200.0d0  
            !stop
           endif   


        end subroutine
  


!--------------------------------------------------------------

        subroutine digit_to_ch ( digit, ch )

        implicit none
        character ch
        integer ( kind = 4 ) digit

        if ( 0 <= digit .and. digit <= 9 ) then
          ch = achar ( digit + 48 )
        else
         ch = '*'
         end if

        return
        end

!*****************************************************************************
       subroutine i4_to_s_left ( i4, s )

        implicit none

        character c
        integer ( kind = 4 ) i
        integer ( kind = 4 ) i4
        integer ( kind = 4 ) idig
        integer ( kind = 4 ) ihi
        integer ( kind = 4 ) ilo
        integer ( kind = 4 ) ipos
        integer ( kind = 4 ) ival
        character ( len = * ) s

         s = ' '
         ilo = 1
         ihi = len ( s )
         if ( ihi <= 0 ) then
          return
         end if

         ival = i4

          if ( ival < 0 ) then
           if ( ihi <= 1 ) then
        s(1:1) = '*'
        return
           end if

         ival = -ival
          s(1:1) = '-'
          ilo = 2

         end if

         ipos = ihi

        do

         idig = mod ( ival, 10 )
         ival = ival / 10

           if ( ipos < ilo ) then
       do i = 1, ihi
        s(i:i) = '*'
       end do
       return
        end if

         call digit_to_ch ( idig, c )

         s(ipos:ipos) = c
         ipos = ipos - 1

        if ( ival == 0 ) then
       exit
        end if

         end do

         s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
         s(ilo+ihi-ipos:ihi) = ' '

         return
        end

!*****************************************************************************
        subroutine vtk_write(output_unit,filename,title,imax,jmax,xyz,nprim,prim)

        USE THERMO 

        implicit none

         integer :: imax, jmax

         character ( len = 20 ) cell_size_string
         integer :: node
         character ( len = 20 ) node_num_string
         integer :: output_unit
         character ( len = 20 ) s_node_num, s_imax, s_jmax, s_cells
         character ( len = 100 ) title
         character ( len = 100 ) filename
         real(kind=8) :: xyz(1:imax,1:jmax,1:3)
         integer :: nprim 
         real(kind=8) :: prim(1:imax,1:jmax,1:nprim)
         real(kind=8) :: uvw(1:imax,1:jmax,1:3), uvw2(1:imax,1:jmax,1:3)
         integer :: i, j, k
         real(kind=8) :: ie, yk(1:5), temp(1:imax,1:jmax) 


         call i4_to_s_left ( imax*jmax, s_node_num )
         call i4_to_s_left ( (imax-1)*(jmax-1), s_cells )  ! don't use kmax for 2D
         call i4_to_s_left ( imax, s_imax )
         call i4_to_s_left ( jmax, s_jmax )

         uvw(1:imax,1:jmax,1) = prim(1:imax,1:jmax,2)
         uvw(1:imax,1:jmax,2) = prim(1:imax,1:jmax,3)
         uvw(1:imax,1:jmax,3) = 0.0d0
 
         uvw2(1:imax,1:jmax,1) = prim(1:imax,1:jmax,12)
         uvw2(1:imax,1:jmax,2) = prim(1:imax,1:jmax,13)
         uvw2(1:imax,1:jmax,3) = 0.0d0

 

         open(output_unit,file=filename,form='formatted')

         write ( output_unit, '(a)' ) '# vtk DataFile Version 2.0'
         write ( output_unit, '(a)' ) title
         write ( output_unit, '(a)' ) 'ASCII'
         write ( output_unit, '(a)' ) 'DATASET STRUCTURED_GRID'
       write ( output_unit, '(a)' ) 'DIMENSIONS ' // (s_imax) //        &
                        (s_jmax) // '1'  ! kmax = 1 for now
       write (output_unit, '(a)' ) 'POINTS ' // (s_node_num) // 'double'


         do j = 1, jmax
          do i = 1, imax
         write ( output_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xyz(i,j,1:3)
          end do
         end do

        write ( output_unit, '(a)' ) 'CELL_DATA ' // (s_cells)
        write ( output_unit, '(a)' ) 'POINT_DATA ' // (s_node_num)

        write ( output_unit, '(a)' ) 'SCALARS pressure double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,5),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'VECTORS velocity double'
        write ( output_unit, *) ((uvw(i,j,1:3),i=1,imax),j=1,jmax)

        write ( output_unit, '(a)' ) 'SCALARS density double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,1),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS internal_energy double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,4),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Y1 double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,6),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Y2 double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,7),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Y3 double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,8),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Y4 double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,9),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Y5 double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,10),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS solid_density double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,11),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'VECTORS solid_velocity double'
        write ( output_unit, *) ((uvw2(i,j,1:3),i=1,imax),j=1,jmax)

        write ( output_unit, '(a)' ) 'SCALARS solid_energy double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,14),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS fign double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,15),i=1,imax),j=1,jmax)  

        
        ! COMPUTE TEMPERATURE   
        do i = 1, imax
         do j = 1, jmax
          ie = prim(i,j,4)
          Yk(1:5) = prim(i,j,6:10)   
        call FROM_IE_TO_T(ie, YK, temp(i,j))
         enddo 
        enddo    

        write ( output_unit, '(a)' ) 'SCALARS temperature double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((temp(i,j),i=1,imax),j=1,jmax)  


        write ( output_unit, '(a)' ) 'SCALARS Jd double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,17),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Jk double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,18),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Jdk double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,19),i=1,imax),j=1,jmax)  


        call flush()
        close(output_unit)

        return
        end subroutine

!---------------------------------------------------------------------------------

!-------------------------------------------------------

        SUBROUTINE HTOTAL(TEMP,HT)

      IMPLICIT NONE
      REAL(KIND=8) :: TEMP, TT, A1, A2, A3, A4, A5, A6, A7, HT, RMW
      REAL(KIND=8) :: HTOINT 

       RMW = 26.98D0
       TT = MIN(TEMP,3000.0D0)

         IF(TT.GT.600.D0) THEN
            A1 =  0.02559589D+02
            A2 = -0.10632239D-03
            A3 =  0.07202828D-06
            A4 = -0.02121105D-09
            A5 =  0.02289429D-13
            A6 =  0.03890214D+06
            A7 =  0.05234522D+02
         ELSE
            A1 =  0.02736825D+02
            A2 = -0.05912374D-02
            A3 = -0.04033937D-05
            A4 =  0.02322343D-07
            A5 = -0.01705599D-10
            A6 =  0.03886794D+06
            A7 =  0.04363879D+02
         END IF

       HTOINT =  A5 * TT**5 / 5.D0 +   &
                 A4 * TT**4 / 4.D0 +   &
                 A3 * TT**3 / 3.D0 +   &
                 A2 * TT**2 / 2.D0 +   &
                 A1 * TT           +   &
                 A6

       HT = HTOINT * 8314.15D0 / RMW

      RETURN
      END SUBROUTINE

!-------------------------------------------------------

 
