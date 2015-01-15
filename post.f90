                                     program post

                 integer, parameter :: imax = 8001
                 integer, parameter :: jmax = 2
                 character ( len = 100 ) filename
                 integer :: n, i, j, l, fno
                 real(kind=8) :: prim(1:imax,1:jmax,1:16)
                 real(kind=8) :: dx, dy
                 real(kind=8) :: rho2, time, Alox(1:2), rho20
                 LOGICAL :: file_exists
                 real(kind=8), dimension(1:500) :: xp, t, maxp, maxsigma


                 xp = 0.0d0
                 t = 0.0d0
                 fno = 0
                 maxp = 0.0d0
                 maxsigma = 0.0d0


                 dx = 4.0d0/(imax-1)
                 dy = dx


                 open(17,file='massAl.txt',form='formatted')
                 open(18,file='Aloxides.txt',form='formatted')

                 do n = 0, 50000, 250    
                  write(filename,'("output/output_",I5.5,".dat")'),n
                      INQUIRE(FILE=filename, EXIST=file_exists)


                

                       if(file_exists) then

                  fno = fno + 1

                  open(34,file=filename,form='formatted')
                  read(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1
                 read(34,*) (prim(i,j,l),l=1,16)
                enddo
               enddo
                  close(34)
                  rho2 = 0.0d0
                  do i = 1, imax-1
                   do j = 1, jmax-1
                     rho2 = rho2 + prim(i,j,11)*dx*dy*1.0d0
                   enddo 
                  enddo

                      if(n.eq.0) then
                       rho20 = rho2
                      endif

                  Alox = 0.0d0
                  do i = 1, imax-1
                   do j = 1, jmax-1
                     Alox(1) = Alox(1) + prim(i,j,1)*prim(i,j,8)*dx*dy*1.0d0
                     Alox(2) = Alox(2) + prim(i,j,1)*prim(i,j,9)*dx*dy*1.0d0
                   enddo 
                  enddo


                  write(17,*) time, rho2/rho20
                  write(18,*) time, Alox(1), Alox(2)
                  print*, 'time, rho2: ', time, rho2, Alox
                  call flush()


                  maxp(fno) = maxval(prim(:,1,5))
                  maxsigma(fno) = maxval(prim(:,1,11))
                  xp(fno) = 0.0d0
                  do i = imax-1, 1, -1
                    if(prim(i,1,5).ge.1.1d5) then
                      xp(fno) = (dble(i)-0.5d0)*dx
                      goto 223  
                    endif
                  enddo

223               continue

                  t(fno) = time


                      ! file_exists
                    endif

                 enddo 


                    ! open(3,file='output_1D',form='formatted')
                    !   do i = 1, imax-1
                    !     write(3,*) (dble(i)-0.5d0)*dx, prim(i,3,1:2), prim(i,3,4:12)
                    !   enddo
                    ! close(3) 


                      open(45,file='xt',form='formatted')
                      do i = 1, fno 
                       write(45,*) t(i)*1.0d3, xp(i), maxp(i)/1.0d5, maxsigma(i)*1.0d3  
                      enddo
                      close(45)

                 close(17)
                 close(18)
                  

                               end program
