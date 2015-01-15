                                     program post

                 integer, parameter :: imax = 8001
                 integer, parameter :: jmax = 2
                 character ( len = 100 ) filename
                 integer :: n, i, j, l
                 real(kind=8) :: prim(1:imax,1:jmax,1:16)
                 real(kind=8) :: dx, dy, time
                 real(kind=8), dimension(1:imax,1:jmax) :: Jd, Jk, Jdk



                 dx = 4.0d0/(imax-1)
                 dy = dx


                 

                 n = 44940
                  write(filename,'("output/output_",I5.5,".dat")'),n
                    

                  open(34,file=filename,form='formatted')
                  read(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1
                 read(34,*) (prim(i,j,l),l=1,16)
                enddo
               enddo
                  close(34)
                  

                  
                  write(filename,'("combustion/J_",I5.5,".dat")'),n
                    

                  open(34,file=filename,form='formatted')
                  read(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1
                 read(34,*) Jd(i,j), Jk(i,j), Jdk(i,j)
                enddo
               enddo
                 close(34)
                  



                        print*, 'n, time ', n, time
              

                     open(3,file='output_1D',form='formatted')
                       do i = 1, imax-1
                         write(3,*) (dble(i)-0.5d0)*dx, prim(i,1,1:2), prim(i,1,4:12)
                       enddo
                     close(3) 


                       open(3,file='J_1D',form='formatted')
                       do i = 1, imax-1
                         write(3,*) (dble(i)-0.5d0)*dx, Jd(i,1), Jk(i,1), Jdk(i,1)
                       enddo
                       close(3) 

                      
                  

                               end program
