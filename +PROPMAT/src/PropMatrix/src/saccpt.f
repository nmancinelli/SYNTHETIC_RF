c------------------------------------------------------------
c     convert output from BYK's programs to SAC format
c------------------------------------------------------------
      dimension  y(5000), nn(10), x(1500)
      character*32 infil, outfil1, outfil2, outfil3
	
      print *, ' input data file =?'
      read(*,*)infil
      iblank=mblank(infil)
      write(*,*)iblank, mblank(infil)
      outfil1=infil(1:iblank)//'.r'
      outfil2=infil(1:iblank)//'.t'
      outfil3=infil(1:iblank)//'.z'      
      ncpt=3
	
      open(1, file=infil)
      mcpt=0                   
      nn(0)=0            
      do 100 icpt=1, ncpt
      read(1,*)
      read(1,*)nn(icpt),dx
      write(*,*) nn(icpt), dx
        mcpt=mcpt+nn(icpt-1)
      write(*,*)mcpt
      read(1,*)(y(i+mcpt),i=1,nn(icpt))
  100 continue
	

      print *,'# of data for each file ='
      print *,(nn(i),i=1,ncpt)
      n=9999999
      do icpt=1,ncpt
      if(nn(icpt).lt.n)n=nn(icpt)
      end do                  
  
      print *,'updated n=',n
      mcpt=0
      m1=0
      do i=2, ncpt
       mcpt=(i-1)*n
       m1=m1+nn(i-1)
       do j=1, n
        y(j+mcpt)=y(j+m1)
       end do	
      end do

c----------------------------------------------------------
c       set scaling parameter for all the data in y
c---------------------------------------------------------
      amax=0
      do i=1,n*ncpt
      if(abs(y(i)).gt.amax)amax=abs(y(i))
      end do

      print *,'calculated maximum=',amax
      print *,'input maximum (0=calculate) :'
      read(*,*)amax1

      if(amax1.eq.0.)then
c--------------------------- use the calculated max.
      scale=amax/10.           
      else
      scale=amax1/10.
      end if





      npps=int(1./dx+0.000001)

      do i=1,n
      x(i)=(i-1)*dx
      end do

c==================================================================
c   now, y contains ncpt components with n data points for each
c   components and with resolution dx.
c   scale y as 0-10, calculate x, and output.
c==================================================================


                               
      do 800 icpt=1,ncpt
                
      mcpt=(icpt-1)*n

      do i=1,n
      y(i+mcpt)=y(i+mcpt)/scale
      end do

      if (icpt.eq.1) then
       call wsac1(outfil1, y, n, 0., dx, neer) 
      end if
c
      if (icpt.eq.2) then
       do j=1, n
        y(j)=y(j+n)
       end do
       call wsac1(outfil2, y, n, 0., dx, neer) 
      end if
c
      if (icpt.eq.3) then
       do j=1, n
        y(j)=y(j+2*n)
       end do
       call wsac1(outfil3, y, n, 0., dx, neer) 
      end if


  800 continue
      
      stop
      end



      function mblank(file)
      character file*32
      do 1 i=1,32
      if(file(i:i).ne.' ') goto 1
      mblank=i-1
      goto 2
 1    continue
 2    end        
