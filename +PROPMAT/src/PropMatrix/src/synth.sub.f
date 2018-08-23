


c************* "synth.sub.ftn" ---> binded with "synth.ftn" *********



c====================================================================
c      calculate generalized elastic moduli  (called by "synth.ftn")
c       after a rotation about symmetry axes;
c   note that this subroutine cannot calculate rotation about
c   an axis other than symmetry axes.
c
c *** The 3 axes of crystal are x1 (N), x2 (E), and x3 (Down) in
c *** the right-hand system. The crystal is then rotated clockwisely
c *** by the input angle. 
c *** 
c
c====================================================================
      subroutine elast  (ifiso,rho,dm,ifsym,azim,c)

c      implicit real*8 (a-h,o-z)
      dimension c(3,3,3,3),d(3,3,3,3)
      real*4  ce(6,6)

      pi=3.1415926536d0
      rad=pi/180.d0

c-----------------------------------------------------------------
c   set initial state :
c   horizontal (a-b) plane, propagation in symmetry plane.
c      ma=1, ifsym=1
c-----------------------------------------------------------------
      ma=1
      ifsym=1
      do 1 n=1,3
      do 1 m=1,3
      do 1 k=1,3
      do 1 j=1,3
    1 c(j,k,m,n)=0.d0
      do 2 j=1,6
      do 2 i=1,6
    2 ce(i,j)=0.
c--------------------------------------------------------------------
c        read in elastic moduli in (gpa = 10**9 pa);
c     =====>1 mbar = 100 gpa
c--------------------------------------------------------------------
      read(1,*)
      read(1,*)ifiso,rho,dm
      read(1,*)ce(1,1),ce(1,2),ce(1,3),ce(1,4),ce(1,5),ce(1,6),
     &                 ce(2,2),ce(2,3),ce(2,4),ce(2,5),ce(2,6),
     &                         ce(3,3),ce(3,4),ce(3,5),ce(3,6),
     &                                 ce(4,4),ce(4,5),ce(4,6),
     &                                         ce(5,5),ce(5,6),
     &                                                 ce(6,6)

      write(*,*) ce(1,1), ce(2,2), ce(3,3), ce(4,4), ce(5,5), ce(6,6)
      
c------------------------------------------------------------------     
c       dump ce(m,n) to c(i,j,k,l)
c------------------------------------------------------------------     
      call dump(ce,c)

c-------------------------------------------------------------------    
c        generate constants to a total of 81
c-------------------------------------------------------------------    
      do 90 i=1,3
      do 90 j=i,3
      do 90 m=i,3
      do 90 n=m,3
      c(i,j,n,m)=c(i,j,m,n)
      c(j,i,m,n)=c(i,j,m,n)
      c(j,i,n,m)=c(i,j,m,n)

      c(m,n,i,j)=c(i,j,m,n)
      c(n,m,i,j)=c(i,j,n,m)
      c(m,n,j,i)=c(j,i,m,n)
   90 c(n,m,j,i)=c(j,i,n,m)

      if(ifiso.eq.1)return

c--------------------------------------------------------------------
c         *** this layer is anisotropic ***
c     choose the desired horizontal symmetry plane :
c           ma=1  test (a-b) plane; skip rotation.
c           ma=2  test (b-c) plane
c           ma=3  test (c-a) plane
c--------------------------------------------------------------------
      print *,' azimuthal variation in (1) a-b (2) b-c (3) c-a ?'
      read(5,*)ma

      if(ma.eq.2)then
c-------to form b-c, form c-b then rotate 90 degrees.
      lx=2    
      th0=90.      
      call elast_rotate (lx,th0,c,d)
      lx=3
      th0=-90.                     
      call elast_rotate (lx,th0,c,d)
      
      else if(ma.eq.3)then

c-------to form c-a, form a-c then rotate 90 degrees.
      lx=1    
      th0=90. 
      call elast_rotate (lx,th0,c,d)
      lx=3
      th0=90.
      call elast_rotate (lx,th0,c,d)

      end if                       

c---------------------------------------------------------------
c          global dip of the original a-axis
c  we do this prior to the azimuthal rotation, so the dip of the
c  original a-axis is fixed relative to the original coordinates.
c  in other words, for propagation in absolute N, old a-axis is 
c  dipping down; for propagation in E, old a-axis is dipping to 
c  the left. the later azimuthal rotation is done about new x3.
c---------------------------------------------------------------
      print *,'dip of a-axis = ?'
      read(*,*)dip
      if(dip.ne.0.)then
      lx=2
      th0=-dip
      call elast_rotate (lx,th0,c,d)
      end if


c---------------------------------------------------------------
c   now, rotate (1,2,3) about lx=3 to (by) the desired azimuth
c        and generate new c(3,3,3,3); lx=3 is + down.
c---------------------------------------------------------------
      print *,'=======> azimuth from x1 (or a-axis) = ?'
      read(*,*)azim

      if(azim.eq.0.)return

      if(azim.ne.0. .and. azim.ne.90.)ifsym=2
                           
      lx=3
      call elast_rotate (lx,azim,c,d)
      return
      end



      subroutine dump(ce,c)
c      implicit real*8 (a-h,o-z)
      dimension  c(3,3,3,3)
      real*4  ce(6,6)
      c(1,1,1,1)=ce(1,1)
      c(1,1,2,2)=ce(1,2)
      c(1,1,3,3)=ce(1,3)
      c(1,1,2,3)=ce(1,4)
      c(1,1,1,3)=ce(1,5)
      c(1,1,1,2)=ce(1,6)

      c(2,2,2,2)=ce(2,2)
      c(2,2,3,3)=ce(2,3)
      c(2,2,2,3)=ce(2,4)
      c(1,3,2,2)=ce(2,5)
      c(1,2,2,2)=ce(2,6)

      c(3,3,3,3)=ce(3,3)
      c(2,3,3,3)=ce(3,4)
      c(1,3,3,3)=ce(3,5)
      c(1,2,3,3)=ce(3,6)

      c(2,3,2,3)=ce(4,4)
      c(1,3,2,3)=ce(4,5)
      c(1,2,2,3)=ce(4,6)

      c(1,3,1,3)=ce(5,5)
      c(1,2,1,3)=ce(5,6)

      c(1,2,1,2)=ce(6,6)

      return
      end



      subroutine elast_rotate (lx,th0,c,d)     
c---------------------------------------------------------------

c   rotate the coordinates (1,2,3) about rotation axis lx by an
c   angle th0 to a new coordinates (1',2',3'),
c   and convert c(i,j,k,l) to c'(i,j,k,l) in the new coordinates.

c      input  : c(3,3,3,3), rotation axis lx, rotation angle th0.
c      output : c(3,3,3,3) at new rotated coordinate.


c    in 2-D system, this is just rotation by a 2-D tensor :

c                 x1'= cos(1',1)*x1+cos(1',2)*x2
c                 x2'= cos(2',1)*x1+cos(2',2)*x2

c    in 4-index, 3-D system, this is :

c                    c'(m1,m2,m3,m4) = 

c    SUM  cos(m1,k1)*cos(m2,k2)*cos(m3,k3)*cos(m4,k4)*c(k1,k2,k3,k4)  
c   k1=1,3
c   k2=1,3
c   k3=1,3
c   k4=1,3
                                       
c      where m refers to axis no. in prime system 
c        and k refers to axis no. in original system.
                   
c            for lx=3

c     if m=k=3,      cos(m,k)=1             e.g., cos(3',3)
c     if m=k.ne.3,   cos(m,k)=cos(theta)    e.g., cos(1',1)
c                                                 cos(2',2)
c     if m=lperm(3)  cos(m,k)=sin(theta),   e.g., cos(1',2)
c     else           cos(m,k)=-sin(theta)   e.g., cos(2',1)

c--------------------------------------------------------------------   

      dimension c(3,3,3,3),d(3,3,3,3),lperm(3)
      rad=3.1415926536d0/180.d0
      lperm(1)=2
      lperm(2)=3
      lperm(3)=1

      co=cos(th0*rad)
      si=sin(th0*rad)

      do 150 m1=1,3
      do 150 m2=1,3
      do 150 m3=1,3
      do 150 m4=1,3

      d(m1,m2,m3,m4)=0.
                       
      do 140 k1=1,3    
                       
       if(k1.eq.m1)then
                       
          if(k1.eq.lx)then
          d1=1.  
          else   
          d1=co
          end if 
                 
       else if (m1.eq.lx.or.k1.eq.lx)then
          d1=0.                    
       else if (m1.eq.lperm(lx))then
          d1=si
       else      
          d1=-si
       end if     
                  
      do 140 k2=1,3
                   
       if(k2.eq.m2)then

          if(k2.eq.lx)then
          d2=1.    
          else     
          d2=co
          end if   
                                                                        

       else if (m2.eq.lx.or.k2.eq.lx)then
          d2=0.              
       else if (m2.eq.lperm(lx))then      
          d2=si
       else                              
          d2=-si
       end if       
                    
      do 140 k3=1,3 
                    
       if(k3.eq.m3)then

          if(k3.eq.lx)then
          d3=1.   
          else    
          d3=co
          end if  
                  
       else if (m3.eq.lx.or.k3.eq.lx)then
          d3=0.                          
       else if (m3.eq.lperm(lx))then
          d3=si
       else       
          d3=-si
       end if     
                  
      do 140 k4=1,3
                   
       if(k4.eq.m4)then
                       
          if(k4.eq.lx)then
          d4=1.           
          else            
          d4=co
          end if          
                          
       else if (m4.eq.lx.or.k4.eq.lx)then
          d4=0.                          
       else if (m4.eq.lperm(lx))then
          d4=si
       else
          d4=-si
       end if

      d(m1,m2,m3,m4)=d(m1,m2,m3,m4)+d1*d2*d3*d4*c(k1,k2,k3,k4)

  140 continue

  150 continue

c--------------------- new elastic moduli -------------------------
      do 160 k1=1,3
      do 160 k2=1,3
      do 160 k3=1,3
      do 160 k4=1,3
  160 c(k1,k2,k3,k4)=d(k1,k2,k3,k4)

      return
      end





c---------------------------------------------------------------------  
c       called by "synth" to create a 6-degree polynomial for q3        
c                q3 : vertical slowness                                 
c---------------------------------------------------------------------  
      subroutine polyf (c,av,rho,p,polyq)
      implicit real*8 (a-h,o-z)
      real*8  neg1,neg2,neg3,p(3,3,0:2),polyq(0:6)
      real*4  c(3,3,3,3),rho
      do 10 n=0,6
   10 polyq(n)=0.d0

      q1=1.d0/av
      q3=1.d0
      do 100 j=1,3
      do 100 m=1,3
      p(j,m,0)=c(j,1,m,1)*q1*q1
      p(j,m,1)=( c(j,1,m,3)+c(j,3,m,1) )*q1*q3
      p(j,m,2)=c(j,3,m,3)*q3*q3
  100 continue
      do 101 j=1,3
  101 p(j,j,0)=p(j,j,0)-rho

      do 140 i=0,2
      do 140 j=0,2
      do 140 k=0,2
      m=i+j+k
      pos1=p(1,1,i)*p(2,2,j)*p(3,3,k)
  140 polyq(m)=polyq(m)+pos1
      do 150 i=0,2
      do 150 j=0,2
      do 150 k=0,2
      m=i+j+k
      pos2=p(2,1,i)*p(3,2,j)*p(1,3,k)
  150 polyq(m)=polyq(m)+pos2
      do 160 i=0,2
      do 160 j=0,2
      do 160 k=0,2
      m=i+j+k
      pos3=p(1,2,i)*p(2,3,j)*p(3,1,k)
  160 polyq(m)=polyq(m)+pos3
      do 170 i=0,2
      do 170 j=0,2
      do 170 k=0,2
      m=i+j+k
      neg1=p(1,3,i)*p(2,2,j)*p(3,1,k)
  170 polyq(m)=polyq(m)-neg1
      do 180 i=0,2
      do 180 j=0,2
      do 180 k=0,2
      m=i+j+k
      neg2=p(2,3,i)*p(3,2,j)*p(1,1,k)
  180 polyq(m)=polyq(m)-neg2
      do 190 i=0,2
      do 190 j=0,2
      do 190 k=0,2
      m=i+j+k
      neg3=p(1,2,i)*p(2,1,j)*p(3,3,k)
  190 polyq(m)=polyq(m)-neg3
      return
      end







c--------------------------------------------------------------------   
c          bairstow's algarithm : called by "synth"                     
c--------------------------------------------------------------------   
      subroutine bairs (n,a,solu)                                    
      implicit real*8 (a-h,o-z)
      real*8 a(13),b(13),c(13),denom,o,r,s,delr,dels,tol,tolp,
     -       xr1,xi1,xr2,xi2,solu(2,13)
      complex solu4

      np3=n+3       

c--------- set initial guess r=0, s=0;  tol0=1.d-7  tolp0=1.d-7% -----  
      tol0=1.d-20
      tolp0=1.d-20
      b(1)=0.
      b(2)=0.
      c(1)=0.
      c(2)=0.
      kr=0                                                              
      r=0.
      s=0.
   80 write(21,'( /''tolp='',e7.2,''%'',3x,''**new solutions**'')')  
     &tolp0*100.
      tol=tol0                                                          
      tolp=tolp0                                                        
   81 write(21,'(1x,''iterations :'')')                                 
                                                                        
      do 20 i=1,80                                                      
      do 10 j=3,np3                                                     
      b(j)=a(j)+r*b(j-1)+s*b(j-2)                                       
   10 c(j)=b(j)+r*c(j-1)+s*c(j-2)                                       
      denom=c(n+1)*c(n+1)-c(n+2)*c(n)                                   
      if(denom.ne.0.d0)go to 12                                         
      r=r+1.d0                                                            
      s=s+1.d0                                                            
      go to 81                                                          
   12 delr=(-b(n+2)*c(n+1)+b(n+3)*c(n))/denom                           
      dels=(-c(n+1)*b(n+3)+c(n+2)*b(n+2))/denom                         
      r=r+delr                                                          
      s=s+dels                                                          
      if(dabs(r).lt.1.d-40)r=0.d0                                         
      if(dabs(s).lt.1.d-40)s=0.d0                                         
                                                                        
      mr=0                                                              
      if(r.eq.0.d0)then                                                 
      if(dabs(delr).lt.tol)mr=1                                         
      else                                                              
      if(dabs(delr/r).lt.tolp)mr=1                                                
      end if                                                            
                                                                        
      ms=0                                                              
      if(s.eq.0.d0)then                                                 
      if(dabs(dels).lt.tol)ms=1                                         
      else                                                              
      if(dabs(dels/s).lt.tolp)ms=1                                                
      end if                                                            
                                                                        
      if(mr.eq.1.and.ms.eq.1)go to 21                                   
   20 continue                                                          
      print *,'accuracy not met in 80 iterations --> relax'            
      tolp=tolp*10.
      write(21,'(1x,''accuracy not met in 80 steps==> relax'')')        
      write(21,'(1x,''tolp='',e7.2,''%'')')tolp*100.
      go to 81                                                          

c===================================================================
c  now, the 1st quradratic equation has been factored from
c  the original polynomial : find the roots (may be complex)
c===================================================================    
   21 write(21,'(15x,''satistied at step '',i2)')i                      
      write(21,'(''quadratic form that was factored out='')')                              
      o=1.d0
      write(21,'(e14.6,'' xx +'',e14.6,'' x +'',e14.6)')o,-r,-s      
      call root(o,-r,-s,xr1,xi1,xr2,xi2)                                
      kr=kr+1                                                           
      solu(1,kr)=xr1
      solu(2,kr)=xi1            

      kr=kr+1                                                           
      solu(1,kr)=xr2
      solu(2,kr)=xi2
      write(21,1001)xr1,xi1,xr2,xi2
 1001 format('roots=','(',2e23.15,')',/6x,'(',2e23.15,')')
                    
c===================================================================
c  here, the reduced equation is checked against 3 possibilities:
c  (1)degree=1 (2)degree=2 (3)degree>2.
c===================================================================    
      n=n-2                                                             
      if(n-2)22,23,24                                                   
c(1)--------- reduced equation is of degree 1                              
   22 continue                                                          
      write(21,'(''reduced linear form='')')                         
      write(21,'(e14.6,'' xx +'',e14.6,'' x +'',e14.6)')b(n+2),      
     &b(n+3)                                                            
      xr1=-b(n+3)/b(n+2)                                                  
      kr=kr+1                                                           
      solu(1,kr)=xr1
      solu(2,kr)=0.
      write(21,*)solu(1,kr),solu(2,kr)                                               

      return
            
                                                            
c(2)--------- reduced equation is of degree 2                              
   23 continue                                                          
      write(21,'( /''reduced quadratic form='')')                    
      write(21,'(e14.6,'' xx +'',e14.6,'' x +'',e14.6)')b(n+1),      
     &  b(n+2),b(n+3)                                                     
      call root(b(n+1),b(n+2),b(n+3),xr1,xi1,xr2,xi2)                   
      kr=kr+1                                                           
      solu(1,kr)=xr1
      solu(2,kr)=xi1
      write(21,1001)xr1,xi1,xr2,xi2

      kr=kr+1                                                           
      solu(1,kr)=xr2
      solu(2,kr)=xi2

      return

                                                                        
c(3)--------- degree > 2; set coeffs. into "a" and get next factor
   24 np3=n+3                                                           
      do 30 i=3,np3                                                     
   30 a(i)=b(i)                                                         
      write(21,'( /''reduced equation='')')                          
      write(21,'(5e14.5)')(a(i),i=3,np3)                                



                                                                        
c====================================================================   
c       this is the case : a(3)*(xxxx) + a(5)*(xx) + a(7) :
c       degree 4 equation for x but quadratic for x*x.
c---------------------------------------------------------------------  


      if(n.ne.4)go to 80                                                
      if(dabs(a(4)).gt.1.d-30.and.dabs(a(6)).gt.1.d-30)go to 80
        call root(a(3),a(5),a(7),xr1,xi1,xr2,xi2)                       
        kr=kr+1  
      write(21,1002)xr1,xi1,xr2,xi2                        
 1002 format('x**2= ','(',2e23.15,')',/6x,'(',2e23.15,')')

c--------------------------------------------------------------------
c  now, xr1,xi1; xr2,xi2 are solution to the 4th degree equation.
c  before take the square-root of the solutions, check if they
c  are complex. if they are, make them complex of real*4 and use
c  "csqrt" to obtain the solution.
c--------------------------------------------------------------------
c------------ 1st solution
      if(dabs(xi1).gt.1.d-7)then
        solu4=cmplx(xr1,xi1)
        solu4=csqrt(solu4)
        solu(1,kr)=real(solu4)
        solu(2,kr)=aimag(solu4)
      else if(xr1.gt.0.)then
        solu(1,kr)=dsqrt(xr1)
        solu(2,kr)=0.
      else
        solu(1,kr)=0.
        solu(2,kr)=dsqrt(-xr1)
      end if

        kr=kr+1                                                         
        solu(1,kr)=-solu(1,kr-1)                                            
        solu(2,kr)=-solu(2,kr-1)
                          
c------------ 2nd solution                                              
        kr=kr+1                                                         
      if(dabs(xi2).gt.1.d-7)then
        solu4=cmplx(xr2,xi2)
        solu4=csqrt(solu4)
        solu(1,kr)=real(solu4)
        solu(2,kr)=aimag(solu4)
      else if(xr2.gt.0.)then
        solu(1,kr)=dsqrt(xr2)
        solu(2,kr)=0.
      else
        solu(1,kr)=0.
        solu(2,kr)=dsqrt(-xr2)
      end if

        kr=kr+1                                                         
        solu(1,kr)=-solu(1,kr-1)                                            
        solu(2,kr)=-solu(2,kr-1)

        return

      end        

                                                       
      subroutine root(a,b,c,xr1,xi1,xr2,xi2)
c--------------------------------------------------
c        root of a(x*x) + b(x) + c = 0.
c--------------------------------------------------
      implicit real*8 (a-h,o-z)                                         
c------this criterion is set particularly for calculating "q3"          
      tol=1.d-20
      t=b*b-4.d0*a*c   
                  
      abst=dabs(t)
      sqabst=dsqrt(abst)
      if(b.eq.0.d0)then
      write(21,'(''b='',e23.15,'' sqrt(bb-4ac)='',e23.15)')b,sqabst
      ratio=1
      else
      ratio=dabs(sqabst/b)
      write(21,'('' sqrt(bb-4ac)/b = '',e23.15)')ratio
      end if

        if(ratio.gt.1.d-7)then                                             
        tt=sqabst
        else                                                            
        tt=0.d0                                                         
        end if                                                          
                                                                        
        if(t.lt.0.d0)then     
c--------------------------- complex roots
        xr1=-b/(2.d0*a)
        xi1=tt/(2.d0*a)
        xr2=xr1
        xi2=-xi1
        else                              
c--------------------------- real roots
        xr1=(-b+tt)/(2.*a)
        xr2=(-b-tt)/(2.*a)
        xi1=0.d0
        xi2=0.d0
        end if                                                          
                                                                        
      return                                                            
      end                                                               
           





      subroutine frank (l,p,q3,f,ifiso,ifsym)
c--------------------------------------------------------
c              called by "synth"
c--------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer max,min                                                   
      dimension p(3,3,0:2),q3(6),f(3,3,6)                             
      do 100 n=1,6                                                      
      q6=q3(n)*q3(n)                                                
      do 90  m=1,3                                                      
      do 90  j=1,3                                                      
   90 f(j,m,n)=p(j,m,0)+p(j,m,1)*q3(n)+p(j,m,2)*q6                    
  100 continue                                                          
                                                                        
c---------------------------------------------------------------------  
c     to separate "sh" in isotropic medium in order to get [f] :        
c                           a b 0                                       
c                  [f] :    c d 0                                       
c                           0 0 e                                       
c----------------------------------------------------------------       
      do 200 n=1,6                                                      
                                                                        
      do 180 i=1,3                                                      
      temp=f(i,2,n)                                                     
      f(i,2,n)=f(i,3,n)                                                 
  180 f(i,3,n)=temp                                                     
      do 190 j=1,3                                                      
      temp=f(2,j,n)                                                     
      f(2,j,n)=f(3,j,n)                                                 
  190 f(3,j,n)=temp                                                     
                                                                        
  200 continue                                                          
                                                                        
      if(ifsym.eq.2)go to 888                                           
                                                                        
c---------------------------------------------------------------------  
c        for iso-  and  aniso-in-symmetry layers                        
c        vertical slownesses are arranged in the order of               
c              p+  p-  sv+  sv-  sh+  sh-                               
c   but,  "q3" may not be in this order when sorted out from "bairs"    
c               ======> check and rearrange [f] and [q3]                
c     q3(n) which makes f(3,3,n) non-zero correspond to "p"           
c     q3(n) which makes f(3,3,n) zero correspond to "sh"              
c---------------------------------------------------------------------  
      tmax= dabs(f(3,3,1))                                              
      max=1                                                             
      do 20 n=3,5,2                                                     
      if( dabs(f(3,3,n)).gt.tmax)then                                   
      tmax= dabs(f(3,3,n))                                              
      max=n                                                             
      end if                                                            
   20 continue                                                          
      if(max.eq.1)go to 29                                              
      do 14 i=1,3                                                       
      do 14 j=1,3                                                       
      temp=f(i,j,1)                                                     
      f(i,j,1)=f(i,j,max)                                               
      f(i,j,max)=temp                                                   
      temp=f(i,j,2)                                                     
      f(i,j,2)=f(i,j,max+1)                                             
   14 f(i,j,max+1)=temp                                                 
      temp=q3(1)                                                      
      q3(1)=q3(max)                                                 
      q3(max)=temp                                                    
      temp=q3(2)                                                      
      q3(2)=q3(max+1)                                               
      q3(max+1)=temp                                                  
c---------------------------------------------------------------------  
c     now, max f(3,3,n) in n=1,2; search for the min and put in n=5,6   
c---------------------------------------------------------------------  
   29 if( dabs(f(3,3,3)).lt. dabs(f(3,3,5)))then                        
                                                                        
      do 17 i=1,3                                                       
      do 17 j=1,3                                                       
      temp=f(i,j,3)                                                     
      f(i,j,3)=f(i,j,5)                                                 
      f(i,j,5)=temp                                                     
      temp=f(i,j,4)                                                     
      f(i,j,4)=f(i,j,6)                                                 
   17 f(i,j,6)=temp                                                     
      temp=q3(3)                                                      
      q3(3)=q3(5)                                                   
      q3(5)=temp                                                      
      temp=q3(4)                                                      
      q3(4)=q3(6)                                                   
      q3(6)=temp                                                      
                                                                        
      end if                                                            

  888 write(21,'(1x,'' "q3" : rearranged vertical slowness='')')        
      do 101 n=1,6                                                      
  101 write(21,'(1x,f20.8)')q3(n)                                     

 
c---------------- write up f matrix ?                                                                       

c      write(21,'(1x,''f='')')                                     
c      do 95 n=1,6                                                       
c      write(21,*)n                                                      
c      do 95 i1=1,3                                                      
c   95 write(21,'(1x,3e13.5)')(f(i1,j1,n)   ,j1=1,3)
                                                                        
      return                                                            
      end                                                               







      subroutine angle(l,q1,q3)
c----------------------------------------------------------
c    "angle" called by synth to calculate the incidence angle
c            for each waves in each layer, knowing q1 and q3.
c     first determine velocity by  sqrt(1/v*v-q1*q1)=q3;
c     then use q3=cos(i)/v
c-----------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension q3(6)  

      rad=180.d0/3.1415926536d0
             
      sqq1=q1*q1
      vp=1.d0/dsqrt(sqq1+q3(1)*q3(1))
      vs1=1.d0/dsqrt(sqq1+q3(3)*q3(3))
      vs2=1.d0/dsqrt(sqq1+q3(5)*q3(5))

      tp=dacos(vp*q3(1))
      ts1=dacos(vs1*q3(3))
      ts2=dacos(vs2*q3(5))

      write(21,'(''incidence angle of P, S1, S2 :'')')
      write(21,'(3f7.3)')tp*rad,ts1*rad,ts2*rad

      return
      end





      subroutine plriz (l,f,a,ifiso,ifsym)
c====================================================================   
c                       called by "synth"                               
c       calculate "polarization vector" <a> for corresponding "q3"      
c             <a>=<x1,x3,x2>=<x,z,y>==<u,w,v>                           
c====================================================================   
c                  for incident "p" wave :                              
c     in isotropic layer or in symmetry plane of anisotropic layer,     
c     f(1,3,n)=f(2,3,n)=f(3,1,n)=f(3,2,n)=0,        so                  
c                        "p","sv" are separateed from "sh"              
c     in isotropic layer, only q3 for n=1,2 lead f(3,3,n) non-zero;     
c     therefore <a> for "p" = <a1,a3,0>                                 
c     q3 for n=3,4,5,6 all lead f(3,3,n)=0.  but we shold determine     
c     <a> for "sv" from the upperleft 2 by 2 submatrix, which will also 
c     have a form of <a1,a3,0>, and <a> for "sh" is ought to be <0,0,1>.
c                        however,                                       
c     in symmetry plane of aniso-layer, n=1,2 plus n=3,4 "or" n=5,6     
c     lead f(3,3,n) non-zero, so <a> for "p" & "sv" can be determined   
c     as <a1,a3,0>;    "sh", corresponding to those q3(n) that make     
c     f(3,3,n) zero, has <a>=<0,0,1>                                    
c                                                                       
c     in off-symmetry plane of aniso-layer,  <a> for one "p" wave and   
c     "two shear waves" are determined through the whole 3 by 3 [f]     
c--------------------------------------------------------------------   
      implicit real*8 (a-h,o-z)
      dimension  f(3,3,6),a(3,6),d(6)


      if(ifsym.ne.1)go to 700                                           
c--------------------------------------------------------------------   
c             in symmetry plane, "sh" needs special treatment           
c--------------------------------------------------------------------   
      do 10 n=1,6                                                       
   10 d(n)= dsqrt(f(1,1,n)*f(1,1,n)+f(1,2,n)*f(1,2,n))                  
                                                                        
      do 13 n=1,4                                                       
      a(1,n)= f(1,2,n)/d(n)                                           
      a(3,n)=-f(1,1,n)/d(n)                                           
   13 a(2,n)=0.d0                                                     
                                                                        
c------------------- in "frank", n=5,6 has been arranged for "sh"       
      do 15 n=5,6                                                       
      a(1,n)=0.
      a(3,n)=0.
   15 a(2,n)=1.
                                                                        
      go to 801                                                         
c-------------------------------------------------------------------    
c        anisotropic medium                                             
c-------------------------------------------------------------------    
  700 do 800 n=1,6                                                      
      d(1)=f(2,2,n)*f(3,3,n)-f(2,3,n)*f(3,2,n)                          
      d(3)=f(2,3,n)*f(3,1,n)-f(2,1,n)*f(3,3,n)                          
      d(2)=f(2,1,n)*f(3,2,n)-f(2,2,n)*f(3,1,n)                          
      d0= dsqrt(d(1)*d(1)+d(3)*d(3)+d(2)*d(2))                          
      a(1,n)=d(1)/d0                                                  
      a(3,n)=d(3)/d0                                                  
  800 a(2,n)=d(2)/d0                                                  
                                                                        
  801 write(21,*)                                                       
      write(21,'(''check if the polarization gives "0"'')')          
      mmm=0                                                             
      do 400 n=1,6                                                      
      do 300 i=1,3                                                      
      d(i)=0.                                                           
      d(i)=d(i)+f(i,1,n)*a(1,n)+f(i,2,n)*a(3,n)+f(i,3,n)*a(2,n)   
      if( dabs(d(i)).gt.1.d-5)then                                      
      print *,dabs(d(i)),'==>polarization vectors are not correct'        
      mmm=1                                                             
      end if                                                            
  300 continue                                                          
  400 write(21,'(3e15.6)')( dabs(d(i)) , i=1,3)                         
                                                                        
      write(21,'(''===(x,y,z) of polarization vectors==='')')  
      do 200 n=1,6                                                      
  200 write(21,'(3f20.8)')(a(jj,n),jj=1,3)                            
c      if(mmm.eq.1)stop                                                  
                                                                        
      return                                                            
      end                                                               






      subroutine conti (nl,l,a,c,av,q3,e,einv,aa)
c--------------------------------------------------------------------
c                called by "synth" 
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension  q3(6),e(6,6),a(3,6),einv(6,6),e1(6,6)           
c           [c] is real*4 generated by sub "elast"
c           [aa] will be multipled by complex variables
      real*4  c(3,3,3,3),aa(6,6,6)                          



c-------------------------------------------
c                     calculate [e]
c----------------------------------------u1
      do 101 n=1,6
  101 e(1,n)=a(1,n)
c----------------------------------------u3
      do 102 n=1,6
  102 e(2,n)=a(3,n)
c----------------------------------------t33
      do 103 n=1,6
      e(3,n)=0.
      do 103 m=1,3
  103 e(3,n)=e(3,n)+a(m,n)*(c(3,3,m,1)+c(3,3,m,3)*av*q3(n))
c----------------------------------------t31
      do 104 n=1,6
      e(4,n)=0.
      do 104 m=1,3
  104 e(4,n)=e(4,n)+a(m,n)*(c(3,1,m,1)+c(3,1,m,3)*av*q3(n))
c----------------------------------------u2
      do 105 n=1,6
  105 e(5,n)=a(2,n)
c----------------------------------------t32
      do 106 n=1,6
      e(6,n)=0.
      do 106 m=1,3
  106 e(6,n)=e(6,n)+a(m,n)*(c(3,2,m,1)+c(3,2,m,3)*av*q3(n))

c=====================================================================
c                     calculate "inverse-e"
c               save 'e' before calling 'gauss'
c=====================================================================
      write(21,'('' e-matrix='')')
      do 116 i=1,6
      write(21,'(6f12.5)')(e(i,j),j=1,6)
      do 116 j=1,6
  116 e1(i,j)=e(i,j)

      do 118 j=1,6
      do 117 i=1,6
  117 einv(i,j)=0.
  118 einv(j,j)=1.

      call gauss (e1,einv,6,6,6)

c--------- IMSL subroutine
c     call linv1f (e1,6,6,einv,9,wk,ier)
c----------IMSL subroutine for complex [e].
c     call leq2c (e, 6, 6, einv, 6, 6,   0,wa,wk,ier)

      do 206 i=1,6
      do 206 j=1,6
       if(dabs( einv(i,j) ) .lt. 1.d-30) einv(i,j)=0.d0
  206 continue

c*****************************************************************
      if(l.eq.nl)return
c=====================================================================
c        construct a matrix [aa] for the phase effects between
c        interfaces.  no phase change in the bottom half space.

c     [ aa(i,j,k) ] < phase(l,n,freq) > = [ al(l,i,k,freq) ]
c=====================================================================


      do 555 k=1,6
      do 555 i=1,6
      do 555 n=1,6
      aa(i,n,k)=e(i,n)*einv(n,k)
  555 continue
  
      return
      end




      subroutine gauss (a,b,n,nr,ndim)
c----------------------------------------------

c      gaussian elimination method with pivoting

c                   solve [a][x]=[b] ;  solution will be in [b]
c          n=number of equations=number of unknowns
c          nr=number of rhs vectors=number of columns of [b]
c          ndim=first dimension of [a] and [b] claimed 
c               in calling program   

c----------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(ndim,n),b(ndim,nr)

      
      np1=n+1
c-------------- now pivoting and decomposition
      do 200 i=1,n-1                
     
      ipvt=i
      ip1=i+1
      do 10 j=ip1,n
      if( dabs(a(ipvt,i)) .lt. dabs(a(j,i)) )ipvt=j
   10 continue

      if( dabs(a(ipvt,i)) .lt. 1.d-6) go to 99
      if(ipvt.eq.i)go to 50
                       
c------------------------- switch
      do 20 jcol=i,n
      save=a(i,jcol)
      a(i,jcol)=a(ipvt,jcol)
   20 a(ipvt,jcol)=save
      do 21 jcol=1,nr
      save=b(i,jcol)
      b(i,jcol)=b(ipvt,jcol)
   21 b(ipvt,jcol)=save   



   50 do 100 jrow=ip1,n
      if(a(jrow,i).eq.0.)go to 100
      ratio=a(jrow,i)/a(i,i)
      do 30 kcol=ip1,n
   30 a(jrow,kcol)=a(jrow,kcol)-ratio*a(i,kcol)
      do 31 kcol=1,nr                             
   31 b(jrow,kcol)=b(jrow,kcol)-ratio*b(i,kcol)
  100 continue

  200 continue  
                                      

      if(dabs(a(n,n)) .lt. 1.d-6)go to 99
c--------------- now, backsubstitution : solve for each rhs vector.

      do 500 kcol=1,nr

      b(n,kcol)=b(n,kcol)/a(n,n)
      do 400 j=2,n
      nvbl=np1-j
      l=nvbl+1
      value=b(nvbl,kcol)
      do 300 k=l,n
      value=value-a(nvbl,k)*b(k,kcol)
  300 continue
      b(nvbl,kcol)=value/a(nvbl,nvbl)
  400 continue

  500 continue
      return

   99 print *,'pivot=',dabs(a(n,n)),'n=',n,'program stops!!!'
      return
      end



      subroutine cgauss (a,b,n,nr,ndim)
c----------------------------------------------
c                 complex, single precision
c      gaussian elimination method with pivoting
c                   solve [a][x]=[b] ;  solution will be in [b]
c          n=number of equations=number of unknowns
c          nr=number of rhs vectors=number of columns of [b]
c          ndim=first dimension of [a] and [b] claimed 
c               in calling program   
c----------------------------------------------
      complex a(ndim,n),b(ndim,nr),save,ratio,value
       
      np1=n+1
c-------------- now pivoting and decomposition
      do 200 i=1,n-1                
     
      ipvt=i
      ip1=i+1
      do 10 j=ip1,n
      if( cabs(a(ipvt,i)) .lt. cabs(a(j,i)) )ipvt=j
   10 continue

      if( cabs(a(ipvt,i)) .lt. 1.e-6) go to 99
      if(ipvt.eq.i)go to 50
                       
c------------------------- switch
      do 20 jcol=i,n
      save=a(i,jcol)
      a(i,jcol)=a(ipvt,jcol)
   20 a(ipvt,jcol)=save
      do 21 jcol=1,nr
      save=b(i,jcol)
      b(i,jcol)=b(ipvt,jcol)
   21 b(ipvt,jcol)=save   



   50 do 100 jrow=ip1,n
      if(a(jrow,i).eq.0.)go to 100
      ratio=a(jrow,i)/a(i,i)
      do 30 kcol=ip1,n
   30 a(jrow,kcol)=a(jrow,kcol)-ratio*a(i,kcol)
      do 31 kcol=1,nr                             
   31 b(jrow,kcol)=b(jrow,kcol)-ratio*b(i,kcol)
  100 continue

  200 continue  
                                      

      if(cabs(a(n,n)) .lt. 1.e-6)go to 99
c--------------- now, backsubstitution : solve for each rhs vector.

      do 500 kcol=1,nr

      b(n,kcol)=b(n,kcol)/a(n,n)
      do 400 j=2,n
      nvbl=np1-j
      l=nvbl+1
      value=b(nvbl,kcol)
      do 300 k=l,n
      value=value-a(nvbl,k)*b(k,kcol)
  300 continue
      b(nvbl,kcol)=value/a(nvbl,nvbl)
  400 continue

  500 continue
      return

   99 print *,'pivot=',cabs(a(n,n)),'cn=',n,'program stops!!!!'

      return
      end





