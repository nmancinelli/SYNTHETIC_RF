
c                 fork.ftn & fold.ftn

c===================================================================
c            1-d fft of the function cx(j),j=1,lx
c                     lx
c   cx(k)= SUM (cx(j)* exp {2*pi*signi*i*(j-1)*(k-1)/lx}
c                    j=1

c        point corresponding to Nyquist frequency is lx/2+1

c===================================================================
      subroutine fork (lx,isigni,cx)
      complex cx(lx),carg,cexp,cw,ctemp
      j=1
c------------- in inverse fft, a factor of (1/lx) is applied.
      sc=1.
      if(isigni.lt.0)sc=1./real(lx)
      do 30 i=1,lx
      if(i.gt.j)go to 10
      ctemp=cx(j)*sc
      cx(j)=cx(i)*sc
      cx(i)=ctemp
   10 m=lx/2
   20 if(j.le.m)go to 30
      j=j-m
      m=m/2
      if(m.ge.1)go to 20
   30 j=j+m
      l=1
   40 istep=2*l
      do 50 m=1,l
      carg=cmplx(0.,1.)*(3.141592654*isigni*(m-1))/real(l)
      cw=cexp(carg)
      do 50 i=m,lx,istep
      ctemp=cw*cx(i+l)
      cx(i+l)=cx(i)-ctemp
   50 cx(i)=cx(i)+ctemp
      l=istep
      if(l.lt.lx)go to 40
      return
      end

c========================================================
c         Nyquist frequency is at nq = lx/2
c       folding across the symmetry axis nq+1=nq1
c   first variable (k=1) folds nowhere, starts from k=2
c========================================================
      subroutine fold  (n,nq1,u)
      complex u(1)   
      np2=n+2
      do 820 k=2,nq1
  820 u(np2-k)=conjg(u(k))
      return
      end

