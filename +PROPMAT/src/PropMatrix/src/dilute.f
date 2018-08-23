
      character*80 input_file,output_file
      dimension c1(6,6),c2(6,6),c3(6,6),c(6,6)

                             
      print *,'input elasticity file = ?'
      read(*,'(a)')input_file
      open(1,file=input_file)

      print *,'output new elasticity file = ?'
      read(*,'(a)')output_file
      open(20,file=output_file,status='new')

      print *,'give weights of the 3 models:'
      read(5,*)w1,w2,w3


      read(1,*)
      read(1,*)ifiso,rho1, dm
      read(1,*)c1(1,1),c1(1,2),c1(1,3),c1(1,4),c1(1,5),c1(1,6),
     &                 c1(2,2),c1(2,3),c1(2,4),c1(2,5),c1(2,6),
     &                         c1(3,3),c1(3,4),c1(3,5),c1(3,6),
     &                                 c1(4,4),c1(4,5),c1(4,6),
     &                                         c1(5,5),c1(5,6),
     &                                                 c1(6,6)

	write(*,*)  ifiso, rho1, c1(1,1)	
      read(1,*)
      read(1,*)ifiso,rho2 ,dm
      read(1,*)c2(1,1),c2(1,2),c2(1,3),c2(1,4),c2(1,5),c2(1,6),
     &                 c2(2,2),c2(2,3),c2(2,4),c2(2,5),c2(2,6),
     &                         c2(3,3),c2(3,4),c2(3,5),c2(3,6),
     &                                 c2(4,4),c2(4,5),c2(4,6),
     &                                         c2(5,5),c2(5,6),
     &                                                 c2(6,6)

	write(*,*) ifiso, rho2, c2(1,1)
      read(1,*)
      read(1,*)ifiso,rho3, dm
      read(1,*)c3(1,1),c3(1,2),c3(1,3),c3(1,4),c3(1,5),c3(1,6),
     &                 c3(2,2),c3(2,3),c3(2,4),c3(2,5),c3(2,6),
     &                         c3(3,3),c3(3,4),c3(3,5),c3(3,6),
     &                                 c3(4,4),c3(4,5),c3(4,6),
     &                                         c3(5,5),c3(5,6),
     &                                                 c3(6,6)
	write(*,*)ifiso, rho3, c3(6,6)
                     
               
      w=w1+w2+w3
      w1=w1/w
      w2=w2/w          
      w3=w3/w
                         
      rho=w1*rho1+w2*rho2+w3*rho3

      do 100 i=1,6
      do 100 j=1,6
      c(j,i)=w1*c1(j,i)+w2*c2(j,i)+w3*c3(j,i)
  100 continue


      write(20,'(3f6.2)')w1,w2,w3
      write(20,'(i3,f7.3)')ifiso,rho
      write(20,1001)c(1,1),c(1,2),c(1,3),c(1,4),c(1,5),c(1,6)
      write(20,1002)       c(2,2),c(2,3),c(2,4),c(2,5),c(2,6)
      write(20,1003)              c(3,3),c(3,4),c(3,5),c(3,6)
      write(20,1004)                     c(4,4),c(4,5),c(4,6)
      write(20,1005)                            c(5,5),c(5,6)
      write(20,1006)                                   c(6,6)

 1001 format(    6f9.3)
 1002 format(9x, 5f9.3)
 1003 format(18x,4f9.3)
 1004 format(27x,3f9.3)
 1005 format(36x,2f9.3)
 1006 format(45x, f9.3)
  
      stop
      end




