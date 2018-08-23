c======================================================================
c                     
c    plane, body-wave synthetic seismograms for anisotropic medium.
c                        part 1 : "synth" :
c          generating transfer function for stacked layers.

c                       b.y.k. 5/13/87

c                     (part 2 is "sourc1")
c
c       for an upcoming incident wave through layered isotropic or
c       anisotropic media.  seismograms correspond to the 3-component
c       displacements at z=(1) surface, &/or (2) bottom interface.
c
c         algorithm is based on keith & crampin (1977 I,II,III)
c
c         vertical slowness   "q" = sqrt [ 1/(v*v)-p*p ]   is real.
c
c         v : p or s velocity
c         p : horizontal slowness = sin(i)/v = ray parameter.
c
c======================================================================
c      <q3>  : vertical slowness; real.
c      [f]   : f matrix; determinant=0; obtained from "frank".
c      <a>   : polarization vectors; obtained from "plriz".
c      [e]   : e matrix; relating <v> and <f>; obtained from "conti".
c      [aa]  : an operating matrix for phase calculation later.
c      [al]  : transfer function for each layer; phase considered.
c      [als] : stacked [al]; not include the bottom layer.
c      [ac]  : [einv][als], matrix relating bottom excitation factor
c                           with surface displacement-stress.
c      [ar]  : reduced [ac] (3 by 3).
c      [ari] : inverse [ar], the transfer function to obtain surface
c              displacements in the subsequent program "sour1".
c      [arb] : transfer function to obtain bottom displacements in
c              the subsequent program "sour1".
c====================================================================== 
      implicit real*8 (a-h,o-z)                                         
      real*4 arg
      character*60 elast_file, name1, synthin                             
      character*60 synthout0, synthout1, synthout2                                     
      dimension  q3(6),f(3,3,6),a(3,6),e(6,6),einv(6,6)
      dimension  p(3,3,0:2),polyq(0:6),abair(13)
      dimension  solu(2,6),rr(3),ss(3)
      complex  ari(3,3),arb(3,3),ar(3,3),carg
c     SIZE OF al DETERMINES MAX NUMBER OF LAYERS
      complex  phase(6),al(100,6,6),als(6,6),ac(6,6)
c------------- everything involved with sub "elast" is real*4.
c     SIZE OF dm DETERMINES MAX NUMBER OF LAYERS
      real*4  c(3,3,3,3),dm(100),rho
c------------- all involved with complex variables are real*4.
      real*4  ph,realpt,imagpt,aa(6,6,6),reinv(6,6)
c  added at Ved's recommendation
      real*8  q1, theta

      
 
      pi=3.1415926536d0                                                 
      rad=pi/180.d0                                                     

      print *,'input elastic constants:'
      read(*,*)name1 
      open(1,file=name1)



      print *,'synth.out filename:'
      read(*,*)synthout0 
      print *,'synthout1 filename:'
      read(*,*)synthout1 
      print *,'synthout2 filename:'
      read(*,*)synthout2 
      
      open(21,file=synthout0)
      open(30,file=synthout1,form='unformatted')
      open(31,file=synthout2,form='unformatted')

      open(40,status='scratch',form='unformatted')
                                                                        
c--------------------------------------------------------------------   
c             input related information                                 
c--------------------------------------------------------------------   
	  print *,'input parms filename:'
      read(*,*)synthin 
      open(2,file=synthin)

      print *,'     velocity of incident wave = ?'                      
      read(2,*)vi                                                       
      print *,'     incidence angle ?'                           
      read(*,*)theta                                                    
      print *,'     how many layers (including half space) ?'           
      read(2,*)nl                                                       
      print *,'     length of the time series you will have later ?'    
      print *,'           = pts of "fft" '                              
      read(2,*)lx                                                       
      print *,'     "npps" points per sec,  "npps"=?'                   
      read(2,*)npps                                                     
      tsec=real(lx)/real(npps)                                          
      write(*,'(1x,''======>'',f7.2,''sec record'')')tsec               
      write(*,'(1x,''cutting freq for transfer function calculation'')')
      write(*,'(1x,''= nc/'',i4,''  nc=?'')')lx
      read(2,*)nc
      print *,' (1) surface motion? (2)motion at depth? (3)both ?'
      read(2,*)motion
              
      
      hv=vi/dsin(theta*rad)
      q1=1.d0/hv


      write(30)lx,npps,nc,motion
      write(30)q1,theta

      write(21,'(1x,''============================================='')')
      write(21,1001)vi,theta,hv
 1001 format(1x,'incident v=',f5.2,3x,'i. angle=',f4.1,3x,'app. v=',
     &f14.8)
      write(21,1002)dcos(theta*rad)/vi
 1002 format(1x,'   for isotropic medium : q3 for i. wave =',e14.6)


c********************************************************************
      do 300 l=1,nl
      write(21,'(//''******** layer '',i2)')l
      print *,'layer',l     
c--------------------------------------------------------------------
c          input elastic moduli;  read file 1 in subroutine.
c--------------------------------------------------------------------
      call elast  (ifiso,rho,dm(l),ifsym,azim,c)
      if (ifiso.eq.2) then
      write(*,*) c(1,1,1,1), c(2,2,2,2), c(3,3,3,3)
      end if
c--------------------------------------------------------------------
c          create 6-degree polynomial for slowness q3 
c                   by setting det [f] = 0
c--------------------------------------------------------------------
      call polyf (c,hv,rho,p,polyq)

      npoly=6
      do 10 i=3,9
   10 abair(i)=polyq(9-i)

c--------------------------------------------------------------------
c           solve the 6-degree polynomial for for "q3"
c--------------------------------------------------------------------
      call bairs (npoly,abair,solu)

      do 90 n=1,6
      q3(n)=solu(1,n)
        if(dabs(solu(2,n)).gt.1.d-6)then
      print *,'========== imaginary vertical slowness ========='
      print *,' inhomogeneous waves have been excited ! stop'
c        stop
        end if
      abair(n)=polyq(0)
      xp=1.
      do 89 i=1,6
      xp=xp*solu(1,n)
   89 abair(n)=abair(n)+polyq(i)*xp
   90 continue

      write(21,'(''residual for each solution'')')
      write(21,'(6e12.4)')(abair(n),n=1,6)   


c---------------------------------------------------------------------
c          calculate matrix [f] ; plugging "q3"
c                    cols correspond to <u,w,v>     (displ.)
c                    rows correspond to <a1,a3,a2>  (polarz.)
c            <q3> will be rearranged in "frank"
c---------------------------------------------------------------------  
      call frank (l,p,q3,f,ifiso,ifsym)

c------------------------------------------------------
c  calculate incidence angle for each wave in each layer
c  and write on 'info.out' for reference.
c  *** after q3 have been rearranged
c------------------------------------------------------
      call angle(l,q1,q3)

c---------------------------------------------------------------------  
c          calculate polarization vector <a>
c---------------------------------------------------------------------  
      call plriz (l,f,a,ifiso,ifsym)

c---------------------------------------------------------------------  
c          calculate [e] & [d] & [al]=[d]*inv[e] for each layer
c      but output [aa] to which the phase of each layer can be
c          easily incorperated in stacking all [al] to [ac]
c---------------------------------------------------------------------  
      call conti (nl,l,a,c,hv,q3,e,einv,aa)
c     has the call to gauss inside it, makes the matrix "a"

      write(40)(q3(n),n=1,6)
      do k=1,6
      do n=1,6
      write(40)(aa(i,n,k),i=1,6)
      end do
      end do


  300 continue
c=====================================================================  
c this is the end of the loop of layers 1 to nl
c=====================================================================  

c---------------------- write down the [e] for the bottom layer.
      do 294 j=1,6
  294 write(30)(e(i,j),i=1,6)


      do 295 j=1,6
      do 295 i=1,6
  295 reinv(i,j)=einv(i,j)

c=====================================================================  
c
c                   go through "frequency loop"
c
c      [ aa(i,n,k) ] < phase(l,n,freq) > = [ al(l,i,k,freq) ]
c             where [al] = [d] [einv] for each layer
c
c  note : imaginary vertical slowness, which correponds to
c         imhomogeneous waves, should be avoided.
c=====================================================================  
      wp=2.*pi*real(npps)/real(lx)

      do 800 ifre=1,nc

      w1=wp*real(ifre-1)

      rewind 40
      do 370 l=1,nl-1


      read(40)(q3(n),n=1,6)
      do k=1,6
      do n=1,6
      read(40)(aa(i,n,k),i=1,6)
      end do
      end do


      c1=w1*dm(l)
      do 348 n=1,6
      arg=c1*q3(n)
      carg=cmplx(0.,1.)*arg
  348 phase(n)=cexp(carg)
      do 352 k=1,6
      do 352 i=1,6
      al(l,i,k)=cmplx(0.,0.)
      do 350 n=1,6
  350 al(l,i,k)=al(l,i,k)+aa(i,n,k)*phase(n)
  352 continue

  	  
  370 continue
c this is the end of the loop of layers 1 to nl-1

c---------------------------------------------------------------------  
c            "stack" the transfer functions between
c       the free surface (0) and the lowest interface (nl-1)
c---------------------------------------------------------------------

c--------stacking matrices from 1st layer
      do 280 j=1,6
      do 280 i=1,6
  280 ac(i,j)=al(1,i,j)
     

      do 297 l=2,nl-1

      do 290 j=1,6
      do 290 i=1,6        
c-- [als] works as an operating matrix here
      als(i,j)=cmplx(0.,0.)
      do 285 k=1,6
  285 als(i,j)=als(i,j)+al(l,i,k)*ac(k,j)
  290 continue
      do 292 j=1,6
      do 292 i=1,6
  292 ac(i,j)=als(i,j)

  297 continue  
      
c---------- return the stacked transfer function to [als]
      do 299 j=1,6
      do 299 i=1,6
      als(i,j)=ac(i,j)
  299 continue
      
c--------------------------- get rid of possible "underflow" later.     
      do 315 i=1,6                                                      
      do 315 j=1,6                                                      
        realpt=real(als(i,j))                                              
        imagpt=aimag(als(i,j))                                              
        if(abs(realpt).lt.1.e-30)realpt=0.
        if(abs(imagpt).lt.1.e-30)imagpt=0.
        als(i,j)=cmplx(realpt,imagpt)                                          
  315 continue                                                          
c-----------------------------------------------------------------------
c        link to the bottom "excitation factor" by multiplying [einv]
c  then, the transfer function of the stacked layers is [ac]
c-----------------------------------------------------------------------
      do 320 i=1,6                                                      
      do 320 j=1,6                                                      
      ac(i,j)=cmplx(0.,0.)                                         
      do 310 k=1,6                                                      
  310 ac(i,j)=ac(i,j)+reinv(i,k)*als(k,j)                                
        realpt=real(ac(i,j))                                               
        imagpt=aimag(ac(i,j))                                               
        if(abs(realpt).lt.1.e-15)realpt=0.
        if(abs(imagpt).lt.1.e-15)imagpt=0.
        ac(i,j)=cmplx(realpt,imagpt)                                           
  320 continue                                                          
                                                                        
c---------------------------------------------------------------------- 
c      solve the 6 by 6 equation  [ac]<v> = <ef> =====>                 
c      by introducing the free surface (3 stress cpts = 0) and          
c      the known upgoing excitation factors (e.g.  1. 0. 0.),            
c      we construct a "reduced"  "3 by 3"  equation [ar]<v>=<fr>.       
c          inverse [ar] = [ari] will be written on file 30.             

c   *** 1,2,5 columns are about displacement.
c   *** 2,4,6 rows are about upward propagation.
c---------------------------------------------------------------------- 
                                                                        
                                                                        
      ar(1,1)=ac(2,1)
      ar(2,1)=ac(4,1)
      ar(3,1)=ac(6,1)
      ar(1,2)=ac(2,2)
      ar(2,2)=ac(4,2)
      ar(3,2)=ac(6,2)
      ar(1,3)=ac(2,5)
      ar(2,3)=ac(4,5)
      ar(3,3)=ac(6,5)
                                 
      do 2 j=1,3
      do 1 i=1,3
    1 ari(i,j)=cmplx(0.,0.)
    2 ari(j,j)=cmplx(1.,0.)

      call cgauss (ar,ari,3,3,3)
                   
c-------- main frame IMSL subroutine :
c      call leq2c (ar, 3, 3, ari, 3, 3,   0,wa,wk,ier)
	  
      write(30)ari         

c      write(21,'(''ari='')')
c      do i=1,3
c      write(21,'(6e12.4)')(ari(i,j),j=1,3)
c      end do

      if(motion.eq.1)go to 800

c--------------------------------------------------------------------   
c        for the motion in the bottom layer, save part of [ac]          
c              [ar] will be written on file 31                          
c--------------------------------------------------------------------   
      arb(1,1)=ac(1,1)
      arb(2,1)=ac(3,1)
      arb(3,1)=ac(5,1)
      arb(1,2)=ac(1,2)
      arb(2,2)=ac(3,2)
      arb(3,2)=ac(5,2)
      arb(1,3)=ac(1,5)
      arb(2,3)=ac(3,5)
      arb(3,3)=ac(5,5)
                      
      write(31)arb    
                      
  800 continue        
                      
      stop            
      end             
