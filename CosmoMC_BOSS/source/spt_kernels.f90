! ****************************************************************
!     kernels.f
!     Author: Rob Speare, 2013-14
!     Purpose: Subroutine/function library written for Bispectrum
!     and Power Spectrum Calculation.
!
! ****************************************************************

! ******************************************************
!     Angle Averaged F3 kernel for the b321_I diagram
! Convention: k3 is the total incoming momentum. p=k1-q. k2 is the
!           lower "leg" of the diagram, which has it's own P(k2)
!           F3(k3, k1-q, q)
! ******************************************************

  module spt_kernels
    use settings, only : dp => mcp
    implicit none
    integer                     :: cyc   ! Triangle Wave Vector Block
    !real(kind=dp)               :: k1, k2, k3
    !real(kind=dp),dimension(3)  :: k1vec, k2vec, k3vec

    contains

    function F3_2D(k3,k2,k1,p,q) ! symmetric in p,q
      implicit none
      real(kind=dp) F3_2D
      real(kind=dp) k1,k2,k3,p,q
      real(kind=dp) t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
      real(kind=dp) t21,t22,t23,t24,t25,t26,t27,t28,t29,t30
      
      t1=1._dp
      t2=2._dp

        t1=1.240079365079365d-4/(k1**4*k2**2)
        t2=p**(-2)
        t3=q**(-2)
        t4=-1.4d1*k1**10
        t5=-2.1d1*k1**8*(2.d0*k3**2+5.d0*(p**2+q**2))
        t6=-3.d0*(k2**2-1.d0*k3**2)**2*(p**2-1.d0*q**2)**2*(2.8d1*k2**2+1.d1*k3**2+7.d0*(p**2+q**2))
        t7=2.d0*k1**6*(-7.d0*k2**4+2.3d1*k3**4+5.6d1*(p**2-1.d0*q**2)**2+5.3d1*k3**2*(p**2+q**2)+k2**2* &
            (6.8d1*k3**2+1.61d2*(p**2+q**2)))
        t8=k1**4
        t9=2.8d1*k2**6+1.d1*k3**6-5.d0*k3**4*(p**2+q**2)
        t10=7.d0*(p**2+q**2)*(p**4-1.d1*p**2*q**2+q**4)-2.d0*k3**2*(3.1d1*p**4-7.4d1*p**2*q**2+3.1d1*q**4)
        t11=-1.d0*k2**4*(3.d1*k3**2+3.29d2*(p**2+q**2))
        t12=k2**2*(2.4d1*k3**4+2.2d1*k3**2*(p**2+q**2)-5.6d1*(7.d0*p**4-1.7d1*p**2*q**2+7.d0*q**4))
        t13=k1**2
        t14=2.8d1*k2**6*(p**2+q**2)-2.d0*k2**4*(-9.1d1*p**4+1.68d2*p**2*q**2-9.1d1*q**4+2.3d1*k3**2*(p**2+q**2))
        t15=k3**2*(1.d1*k3**4*(p**2+q**2)+7.d0*(p**2-1.d0*q**2)**2*(p**2+q**2)-2.d0*k3**2*(5.d0*p**4-2.4d1*p**2*q**2+5.d0*q**4))
        t16=k2**2*(8.d0*k3**4*(p**2+q**2)+7.d0*(p**2-1.d0*q**2)**2*(p**2+q**2)-8.d0*k3**2*(1.2d1*p**4-1.7d1*p**2*q**2+1.2d1*q**4))
        t17=(-9.920634920634921d-4*k1*(k3-1.d0*p))/(k2**2*p**2)
        t18=(k3+p)*(2.d0*k3**2+7.d0*p**2)*(k2-1.d0*q)**2
        t19=q**(-2)
        t20=(k2+q)**2
        t21=k2**4*p**2+k1**2*(k3**2-1.d0*p**2)*(k2**2-1.d0*q**2)
        t22=k3**2*q**2*(k3**2-1.d0*p**2+q**2)-1.d0*k2**2*(-1.d0*p**4+p**2*q**2+k3**2*(p**2+q**2))
        t23=(-9.920634920634921d-4*k1*(k2-1.d0*p)**2)/(k2**2*p**2)
        t24=((k2+p)**2*(k3-1.d0*q))/q**2
        t25=k3+q
        t26=2.d0*k3**2+7.d0*q**2
        t27=k2**4*q**2+k1**2*(k2**2-1.d0*p**2)*(k3**2-1.d0*q**2)
        t28=k3**2*p**2*(k3**2+p**2-1.d0*q**2)-1.d0*k2**2*(q**2*(p**2-1.d0*q**2)+k3**2*(p**2+q**2))
        t29=t1*t2*t3*(t4+t5+t6+t7+t8*(t9+t10+t11+t12)+2.d0*t13*(t14+t15+t16))
        t30=dsqrt(1/(t21+t22))*t17*t18*t19*t20+dsqrt(1/(t27+t28))*t23*t24*t25*t26
      F3_2D=t29+t30
    end function F3_2D

!    function F3_2D(k3,k2,k1,p,q) ! symmetric in p,q
!      implicit none
!      real(kind=dp) F3_2D
!      real(kind=dp) k1,k2,k3,p,q
!      real(kind=dp) t1,t2
!      t1 = 1._dp
!      t2 = 2._dp
!      F3_2D=t1+t2
!    end function F3_2D

! ******************************************************
! The F2 kernel, called with three scalar momenta.
! k is the total incoming momenta. q and p are the "legs"
! This kernel is symmetric in p and q, and one can also
! swap the sign of p and q, yieldling the same result.
! ******************************************************
    function F2hat(k,p,q) ! symmetric in p,q!
      implicit none
      real(kind=dp) F2hat
      real(kind=dp) k,p,q,cos
      cos=(k**2-q**2-p**2)/(2._dp*p*q)
      F2hat=5._dp/7._dp+0.5_dp*cos*(p/q+q/p)+2._dp/7._dp*cos**2
    end function F2hat

! ******************************************************
! The symmetric F2 and G2 kernels, called with two wave vectors
! ******************************************************
    function F2s(v1,v2) ! checked
      implicit none
      real(kind=dp) F2s
      real(kind=dp) v1(3),v2(3),x,av1,av2
!      F2s=(F2(v1,v2)+F2(v2,v1))/2._dp
      av1=av(v1)
      av2=av(v2)
      x=dot_product(v1,v2)/av1/av2
      F2s=5._dp/7._dp+0.5_dp*x*(av1/av2+av2/av1)+2._dp/7._dp*x**2
!      F2s=(5._dp/7._dp)+0.5_dp*dot_product(v1,v2)*(dot_product(v1,v1)+dot_product(v2,v2))
!     &/(dot_product(v1,v1)*dot_product(v2,v2))+(2._dp/7._dp)*dot_product(v1,v2)*dot_product(v1,v2)
!     &/(dot_product(v1,v1)*dot_product(v2,v2))
    end function F2s

! ******************************************************
    function G2s(v1,v2) ! checked
      implicit double precision (A-Z)
      real(kind=dp)       :: G2s
      real(kind=dp)       :: v1(3),v2(3)
      av1=av(v1)
      av2=av(v2)
      x=dot_product(v1,v2)/av1/av2
      G2s=3._dp/7._dp+0.5_dp*x*(av1/av2+av2/av1)+4._dp/7._dp*x**2
    end function G2s
   
! ***************************************************
    function av(wav) ! checked
      implicit none
      real(kind=dp) av 
      real(kind=dp) wav(3),temp
      temp=wav(1)*wav(1)+wav(2)*wav(2)+wav(3)*wav(3)
      av = dsqrt(temp)
    end function av

! ***************************************************
    subroutine kvectors(k1,k2,k3,k1vec,k2vec,k3vec) !checked
      implicit none
      real(kind=dp) k1,k2,k3,k1vec(3),k2vec(3),k3vec(3)
      real(kind=dp) cosi,sine,zero
      zero = 0._dp
      if(abs(k1-k2-k3) .lt. k3*1.d-9) then
        cosi=-1._dp
        sine=0._dp
      else
        cosi=(k3**2-k2**2-k1**2)/(2._dp*k1*k2)
        sine=dsqrt(max(0._dp,1._dp-cosi**2)) ! safety
      endif
      k1vec(1)=0._dp
      k1vec(2)=0._dp
      k1vec(3)=k1
      k2vec(1)=k2*sine
      k2vec(2)=0
      k2vec(3)=k2*cosi
      k3vec(1)=-k2vec(1)!-k1vec(1)
      k3vec(2)=0._dp
      k3vec(3)=-k2vec(3)-k1vec(3)
    end subroutine kvectors

! ***************************************************
    subroutine permute(k1, k2, k3, k1vec, k2vec, k3vec) !checked
      implicit none
      real(kind=dp)               :: k1, k2, k3
      real(kind=dp),dimension(3)  :: k1vec, k2vec, k3vec
      real(kind=dp)               :: aw1, aw2, aw3
      real(kind=dp),dimension(3)  :: w1, w2, w3
      aw1 = k1
      aw2 = k2
      aw3 = k3
      w1 = k1vec
      w2 = k2vec
      w3 = k3vec
      k2 = aw1
      k3 = aw2
      k1 = aw3
      k2vec = w1
      k3vec = w2
      k1vec = w3
    end subroutine permute

! ***************************************************
    function WTH(x) !
      real(kind=dp) WTH
      real(kind=dp) x
      WTH=3._dp*(sin(x)-x*cos(x))/x/x/x
    end function WTH

! ***************************************************

    function number_of_triangles(nkmin,nkmax,step)
      implicit none
      integer number_of_triangles
      integer i,j,l,nkmin,nkmax,step,nT
      nT=0
      do i=nkmin,nkmax,step
        do j=nkmin,i,step
          do l=max(nkmin,abs(i-j)),j,step
            nT=nT+1
          enddo
        enddo
      enddo
      number_of_triangles=nT
    end function number_of_triangles

    subroutine diagnose(k1, k2, k3, tclass)
      implicit none
      real(kind=dp),intent(in)  :: k1, k2, k3
      integer, intent(out)      :: tclass
      real(kind=dp)             :: min,  max, eps

      min = dmin1(k1,k2,k3)
      max = dmax1(k1,k2,k3)
      eps = min*10.e-6_dp
      if (k1 == k2) then
        if(k2 == k3) then
          tclass = 4
!         write(*,*) k1,k2,k3,'equilateral'
        else
          tclass = 3
!         write(*,*) k1,k2,k3,'isosceles'
        endif
      else
        if(k2 == k3) then
          tclass = 1
!         write(*,*) k1,k2,k3,'isosceles'
        else
          if (k1 == k3) then
!           write(*,*) k1,k2,k3,'isosceles'
            tclass = 2
          else
!           write(*,*) k1,k2,k3,'undef'
            tclass = -1
          endif
        endif
      endif

      if(dabs(k1-k2-k3) < eps .or. dabs(k2-k1-k3) < eps .or. dabs(k3-k2-k1) < eps) then
        tclass = 0
!         write(*,*) k1,k2,k3,'colinear'
      endif

    end subroutine diagnose
   
  end module spt_kernels




