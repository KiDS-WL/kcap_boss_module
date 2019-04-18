!este modulo calcula las funciones de distancia para una dada cosmologia

  module numrec
    implicit none
!******************************************************************
!* NUMERICAL RECIPIES
!******************************************************************
    contains
                                                                                
    subroutine pinter(npts,xx,yy,xv,yv)
      implicit none
      integer,intent(in)                       :: npts
      integer                                  :: klo, khi, kp
      real(kind=8),dimension(npts),intent(in) :: xx, yy
      real(kind=8),intent(in)                 :: xv
      real(kind=8),intent(out)                :: yv
      real(kind=8)                            :: dx
      if(xx(1) > xv) then
     !  extrapolate solution
        klo = 1
        khi = 2
      elseif(xx(npts) < xv)then
     !  extrapolate solution
        klo = npts-1
        khi = npts
      else
        klo = 1
        khi = npts
        do while (khi-klo > 1)
          kp = (khi+klo)/2
          if(xx(kp).gt.xv) then
            khi = kp
          else
            klo = kp
          end if
        end do
      end if
      dx = (xv - xx(klo))/(xx(khi)-xx(klo))
      yv = dx*yy(khi) + (1. - dx)*yy(klo)
      return
    end subroutine pinter



      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      DOUBLE PRECISION ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
!     PARAMETER (MBIG=4000000.d0,MSEED=1618033.d0,MZ=0.d0,FAC=1.d0/MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END function ran3

      FUNCTION gasdev(idum)
      INTEGER idum
      DOUBLE PRECISION gasdev
!U    USES ran1
      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.d0*ran3(idum)-1.d0
        v2=2.d0*ran3(idum)-1.d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.d0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END FUNCTION gasdev


      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-5, JMAX=22, JMAXP=JMAX+1, K=5, KM=K-1)
!U    USES polint,trapzd
      INTEGER j
      DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      write(*,*)'too many steps in qromb'
      END subroutine qromb

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
    END subroutine trapzd


    SUBROUTINE qromo(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
!     USES polint
      INTEGER j
      DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do 11 j=1,JMAX
        call midpnt(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.d0
11    continue
      write(*,*) 'too many steps in qromo'
    END subroutine qromo


      SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      DOUBLE PRECISION dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=20,MMAX=20)
!CU    USES polint
      INTEGER j,k
      DOUBLE PRECISION ymtmp(MMAX),yntmp(NMAX)
    
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
          !write(*,*)x2a(k),ya(j,k)
11      continue
        write(*,*)'interpolo 1',n,x2
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
        !write(*,*)x2,ymtmp(j)
        !read(*,*)
12    continue
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      END subroutine polin2
!C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

SUBROUTINE polint(xa,ya,n,x,y,dy)
  INTEGER n,NMAX
  DOUBLE PRECISION dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER i,m,ns
  DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=abs(x-xa(1))
  do 11 i=1,n
    dift=abs(x-xa(i))
    if (dift.lt.dif) then
      ns=i
      dif=dift
    endif
    c(i)=ya(i)
    d(i)=ya(i)
11  continue
    y=ya(ns)
    ns=ns-1
    do 13 m=1,n-1
      do 12 i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.d0)then
          write(*,*)'failure in polint'
          read(*,*)
        end if
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
12    continue
      if (2*ns.lt.n-m)then
        dy=c(ns+1)
      else
        dy=d(ns)
        ns=ns-1
      endif
      y=y+dy
13  continue
  return
END subroutine polint
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

!**************************************************************************

      SUBROUTINE midpnt(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5d0*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.d0*tnm)
        ddel=del+del
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.d0
      endif
      return
      END subroutine midpnt
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

!************************************************************************

      FUNCTION gammln(xx)
      double precision gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
     &24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
     &-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      end do
      gammln=tmp+log(stp*ser/x)
      return
      END function gammln

!************************************************************************

      FUNCTION gammq(a,x)
      DOUBLE PRECISION a,gammq,x
!U    USES gcf,gser
      DOUBLE PRECISION gammcf,gamser,gln
      if(x.lt.0.d0.or.a.le.0.d0)then
        write(*,*) 'bad arguments in gammq'
        read(*,*)
      end if
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.d0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END function gammq
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

!************************************************************************
      SUBROUTINE gser(gamser,a,x,gln)
      implicit none
      INTEGER ITMAX
      DOUBLE PRECISION a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
!U    USES gammln
      INTEGER n
      DOUBLE PRECISION ap,del,sum
      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0)then
          write(*,*) 'x < 0 in gser'
          read(*,*)
        end if
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      write(*,*)'a too large, ITMAX too small in gser'
      read(*,*)
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END subroutine gser
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.
                                                                    
!************************************************************************
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      DOUBLE PRECISION a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
!U    USES gammln
      INTEGER i
      DOUBLE PRECISION an,b,c,d,del,h   
      gln=gammln(a)
      b=x+1.d0-a
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(abs(d) < FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c) < FPMIN)c=FPMIN
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0)< EPS)goto 1
11    continue
      write(*,*)'a too large, ITMAX too small in gcf'
      read(*,*)
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END subroutine gcf
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.
  
    
      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      DOUBLE PRECISION x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
!U    USES spline
      INTEGER j,k
      DOUBLE PRECISION y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      END SUBROUTINE splie2
!C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.


      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=3000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+&
     &1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*&
     &u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      END SUBROUTINE spline

      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      DOUBLE PRECISION x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
!CU    USES spline,splint
      INTEGER j,k
      DOUBLE PRECISION y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      END SUBROUTINE SPLIN2
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.


      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) then
        write(*,*) 'bad xa input in splint'
        read(*,*)
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**&
     &2)/6.d0
      return
      END SUBROUTINE splint
!C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

      
      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          write(*,*)'NSTACK too small in sort'
          read(*,*)
        end if
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END subroutine sort
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

      SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      DOUBLE PRECISION arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)then
          write(*,*) 'NSTACK too small in sort2'
          read(*,*)
        end if
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END subroutine sort2

    SUBROUTINE gaussj(a,n,np,b,m,mp)
      implicit none
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION a(np,np),b(np,mp)
      PARAMETER (NMAX=150)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      DOUBLE PRECISION big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(*,*) 'singular matrix in gaussj'
                read(*,*)
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) then
          write(*,*)'singular matrix in gaussj'
          read(*,*)
        end if
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
    END subroutine gaussj
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

      SUBROUTINE savgol(c,np,nl,nr,ld,m)
      INTEGER ld,m,nl,np,nr,MMAX
      DOUBLE PRECISION c(np)
      PARAMETER (MMAX=6)
!U    USES lubksb,ludcmp
      INTEGER imj,ipj,j,k,kk,mm,indx(MMAX+1)
      DOUBLE PRECISION d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
      if(np.lt.nl+nr+&
     &1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m)&
     &then
        write(*,*)'bad args in savgol'
        read(*,*)
      end if
      do 14 ipj=0,2*m
        sum=0.d0
        if(ipj.eq.0)sum=1.d0
        do 11 k=1,nr
          sum=sum+ dble(k)**ipj
11      continue
        do 12 k=1,nl
          sum=sum+ dble(-k)**ipj
12      continue
        mm=min(ipj,2*m-ipj)
        do 13 imj=-mm,mm,2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
13      continue
14    continue
      call ludcmp(a,m+1,MMAX+1,indx,d)
      do 15 j=1,m+1
        b(j)=0.d0
15    continue
      b(ld+1)=1.d0
      call lubksb(a,m+1,MMAX+1,indx,b)
      do 16 kk=1,np
        c(kk)=0.d0
16    continue
      do 18 k=-nl,nr
        sum=b(1)
        fac=1.d0
        do 17 mm=1,m
          fac=fac*k
          sum=sum+b(mm+1)*fac
17      continue
        kk=mod(np-k,np)+1
        c(kk)=sum
18    continue
      return
      END SUBROUTINE savgol

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END SUBROUTINE lubksb

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0)then
          write(*,*)'singular matrix in ludcmp'
          read(*,*)
        end if
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END SUBROUTINE ludcmp


      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      DOUBLE PRECISION ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=5000)
      EXTERNAL funk
!CU    USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      DOUBLE PRECISION rtol,sum,swap,ysave,ytry,psum(NMAX)
      iter=0
1     do 12 n=1,ndim
        sum=0.d0
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
      if (iter.ge.ITMAX) then
        write(*,*)'ITMAX exceeded in amoeba'
        read(*,*)
      end if
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0d0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END subroutine amoeba
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.


      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      DOUBLE PRECISION amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
!CU    USES funk
      INTEGER j
      DOUBLE PRECISION fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.d0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END function amotry
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

     SUBROUTINE gauleg(x1,x2,x,w,n) 
      INTEGER n 
      DOUBLE PRECISION x1,x2,x(n),w(n) 
      DOUBLE PRECISION EPS 
      PARAMETER (EPS=3.d-14) 
      INTEGER i,j,m 
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1 
      m=(n+1)/2 
      xm=0.5d0*(x2+x1) 
      xl=0.5d0*(x2-x1) 
      do 12 i=1,m 
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0)) 
1       continue 
          p1=1.d0 
          p2=0.d0 
          do 11 j=1,n 
            p3=p2 
            p2=p1 
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j 
11        continue 
          pp=n*(z*p1-p2)/(z*z-1.d0) 
          z1=z 
          z=z1-p1/pp 
        if(abs(z-z1).gt.EPS)goto 1 
        x(i)=xm-xl*z 
        x(n+1-i)=xm+xl*z 
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp) 
        w(n+1-i)=w(i) 
12    continue 
      return 
      END SUBROUTINE gauleg 
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.


      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      DOUBLE PRECISION a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      DOUBLE PRECISION c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.d0
11      continue
        v(ip,ip)=1.d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.d0
23      continue
24    continue
      write(*,*)'too many iterations in jacobi'
      read(*,*)
      return
      END subroutine jacobi
!  (C) Copr. 1986-92 Numerical Recipes Software 63$&.

      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      INTEGER n
      DOUBLE PRECISION sj,sjp,sy,syp,x
!CU    USES bessjy
      DOUBLE PRECISION factor,order,rj,rjp,ry,ryp,RTPIO2
      PARAMETER (RTPIO2=1.2533141d0)
      if(n.lt.0.or.x.le.0.d0)then
        write(*,*)'bad arguments in sphbes'
        read(*,*)
      end if
      order=n+0.5d0
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.d0*x)
      syp=factor*ryp-sy/(2.d0*x)
      return
      END SUBROUTINE  sphbes

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      INTEGER MAXIT
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-10,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,PI=3.141592653589793d0)
!CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,&
     &f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,&
     &r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,&
     &temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0.d0.or.xnu.lt.0.d0)then
        write(*,*)'bad arguments in bessjy'
        read(*,*)
      end if
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      write(*,*)'x too large in bessjy; try asymptotic expansion'
      read(*,*)
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        write(*,*)'bessy series failed to converge'
        read(*,*)
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        write(*,*) 'cf2 failed in bessjy'
        read(*,*)
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END SUBROUTINE bessjy

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=5,NUSE2=5)
!     USES chebev
      DOUBLE PRECISION xx,c1(7),c2(8)
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,&
     &-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,&
     &-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END SUBROUTINE beschb

      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      DOUBLE PRECISION chebev,a,b,x,c(m)
      INTEGER j
      DOUBLE PRECISION d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.d0)then
        write(*,*)'x not in range in chebev'
        read(*,*)
      end if
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END FUNCTION chebev


      FUNCTION plgndr(l,m,x)
      INTEGER l,m
      DOUBLE PRECISION plgndr,x
      INTEGER i,ll
      DOUBLE PRECISION fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0)then
        write(*,*)'bad arguments in plgndr'
        read(*,*)
      end if
      pmm=1.d0
      if(m.gt.0) then
        somx2=sqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END FUNCTION plgndr


  end module numrec







