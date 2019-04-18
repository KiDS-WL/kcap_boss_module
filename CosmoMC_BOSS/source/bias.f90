!! scaling version !!
module bias
  use settings, only : dp => mcp
  use quadpack 
  use cubaint
  use NLRSD_settings, only : fac_norm, matterpowerat  
  implicit none
  real(kind=dp),  private, parameter :: kmin = 0.0001_dp, kmax = 10._dp
  real(kind=dp),  private, parameter :: relacc = 0.01_dp, absacc = 0.1_dp
  real(kind=dp),  private            :: k

  contains

    !! b1b2 !!

   function Pk_b1b2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b1b2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer  :: method
     k = arg_k
     epsabs    = absacc/10._dp
     epsrel    = relacc/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(b1b2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_b1b2   = ans(1)/fac_norm
   end function Pk_b1b2

   subroutine b1b2_int(ndim,x,ncomp,f)
     implicit   none
     integer    :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,mu12,powq,powp,f2

     q    = kmax*x(1)+kmin*(1._dp-x(1))
     pmin = max(q,dabs(k-q))
     pmax = min(kmax,k+q)

     if (pmin.lt.pmax) then
       p    = pmax*x(2)+pmin*(1._dp-x(2))
       mu12 = (k**2-q**2-p**2)/(2._dp*p*q)
       powq = matterpowerat(q)
       powp = matterpowerat(p)        
       f2   = 5._dp/7._dp+0.5_dp*mu12*(p/q+q/p)+2._dp/7._dp*mu12**2
       f(1) = 2._dp*(4._dp*pi/k)*q*p*powq*powp*f2*(kmax-kmin)*(pmax-pmin)
     else
       f(1)=0._dp
     endif

   end subroutine b1b2_int

   !! b1g2 !!

   function Pk_b1g2_mc(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b1g2_mc,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer         :: method

     k = arg_k

     epsabs    = absacc/10._dp
     epsrel    = relacc/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(b1g2mc_int,method,epsrel,epsabs,ans,error,prob)
     Pk_b1g2_mc = ans(1)/fac_norm

   end function Pk_b1g2_mc

   subroutine b1g2mc_int(ndim,x,ncomp,f)
     implicit   none
     integer           :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,mu12,kernel,powq,powp,f2

     q    = kmax*x(1)+kmin*(1._dp-x(1))
     pmin = max(q,dabs(k-q))
     pmax = min(kmax,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        mu12 = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q)
        powp = matterpowerat(p)        
        f2   = 5._dp/7._dp+0.5_dp*mu12*(p/q+q/p)+2._dp/7._dp*mu12**2  
        kernel = (mu12**2._dp-1._dp)*f2
        f(1)   = 4._dp*(4._dp*pi/k)*q*p*powq*powp*kernel*(kmax-kmin)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine b1g2mc_int

   function Pk_b1g2_prop(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b1g2_prop,arg_k,epsabs,epsrel,a,b
     integer, parameter :: limit=1000, key = 2 !! local integration rule !!
     real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
     integer  :: neval,ier,last,iord(limit)
     real(kind=dp)   :: ans

     k = arg_k
     a = kmin      ! lower limit of integration
     b = kmax      ! upper limit of integration
     epsrel = relacc/10._dp
     epsabs = absacc/10._dp
     ans = 0._dp

     call qage (b1g2prop_dq, a, b, epsabs, epsrel, key, limit, ans, &
     abserr, neval, ier, alist, blist, rlist, elist, iord, last)
     Pk_b1g2_prop = -4._dp*pi*ans/fac_norm !remove *matterpowerat(k) to include Pdv later
     if(ier /= 0)  call warning('Pk_b1g2_prop',ier)  

   end function Pk_b1g2_prop

   function b1g2prop_dq(q)        
     implicit none
     real(kind=dp)   :: b1g2prop_dq,q,kernel

     kernel = ((k**2 + q**2)*(33*k**4 + 14*k**2*q**2 + 33*q**4))/(42.*k**2*q**2) +               &
            ((k**2 - q**2)**2*(11*k**4 + 34*k**2*q**2 + 11*q**4)*log((k - q)**2/(k + q)**2))/  &
            (56.*k**3*q**3)
     b1g2prop_dq = kernel*matterpowerat(q)

   end function b1g2prop_dq

!! b2b2 !!

   function Pk_b2b2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b2b2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer  :: method

     k = arg_k

     epsabs  = absacc/10._dp
     epsrel  = relacc/10._dp
     method  = 3 !4:for high lambda, divonne is more accurate than cuhre
     call    integrate2D(b2b2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_b2b2 = ans(1)/fac_norm

   end function Pk_b2b2

   subroutine b2b2_int(ndim,x,ncomp,f)
     implicit   none
     integer    :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,kernel,powq,powp

     q    = kmax*x(1)+kmin*(1._dp-x(1))
     pmin = max(kmin,dabs(k-q))
     pmax = min(kmax,k+q)

     if (pmin.lt.pmax) then
       p    = pmax*x(2)+pmin*(1._dp-x(2))
       powq = matterpowerat(q)
       powp = matterpowerat(p)          
       kernel = (1._dp-powq/powp)       
       f(1)   = (1._dp/2._dp)*(2._dp*pi/k)*q*p*powq*powp*kernel*(kmax-kmin)*(pmax-pmin)
     else
       f(1)=0._dp
     endif

   end subroutine b2b2_int

!! b2g2 !!

   function Pk_b2g2(arg_k)
     implicit none
     real(kind=dp)   :: Pk_b2g2,arg_k,epsabs,epsrel
     real(kind=dp)   :: ans(1),error(1),prob(1)
     integer  :: method

     k = arg_k

     epsabs  = absacc/10._dp
     epsrel  = relacc/10._dp
     method  = 3 !4:for high lambda, divonne is more accurate than cuhre
     call    integrate2D(g2g2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_b2g2 = ans(1)/fac_norm

   end function Pk_b2g2

   subroutine b2g2_int(ndim,x,ncomp,f)
     implicit   none
     integer    :: ndim,ncomp
     real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,kernel,mu12,powq,powp

     q    = kmax*x(1)+kmin*(1._dp-x(1))
     pmin = max(q,dabs(k-q))
     pmax = min(kmax,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        mu12 = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q)
        powp = matterpowerat(p)          
        kernel = (mu12**2._dp-1._dp)   
        f(1)   = 2._dp*(4._dp*pi/k)*q*p*powq*powp*kernel*(kmax-kmin)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine b2g2_int

!! g2g2 !!

   function Pk_g2g2(arg_k)
   implicit none
   real(kind=dp)   :: Pk_g2g2,arg_k,epsabs,epsrel
   real(kind=dp)   :: ans(1),error(1),prob(1)
   integer  :: method

     k = arg_k

     epsabs  = absacc/10._dp
     epsrel  = relacc/10._dp
     method  = 3 !4:for high lambda, divonne is more accurate than cuhre
     call    integrate2D(g2g2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_g2g2 = ans(1)/fac_norm

   end function Pk_g2g2

   subroutine g2g2_int(ndim,x,ncomp,f)
   implicit   none
   integer    :: ndim,ncomp
   real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,mu12,kernel,powq,powp

     q    = kmax*x(1)+kmin*(1._dp-x(1))
     pmin = max(q,dabs(k-q))
     pmax = min(kmax,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        mu12 = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q)
        powp = matterpowerat(p)          
        kernel = (mu12**2._dp-1._dp)**2._dp     
        f(1)   = 2._dp*(4._dp*pi/k)*q*p*powq*powp*kernel*(kmax-kmin)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine g2g2_int

!! b1g3m !!

   function Pk_b1g3m_prop(arg_k)
   implicit none
   real(kind=dp)   :: Pk_b1g3m_prop,arg_k,epsabs,epsrel,a,b
   integer, parameter :: limit=1000, key = 2 !! local integration rule !!
   real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
   integer  :: neval,ier,last,iord(limit)
   real(kind=dp)   :: ans

     k = arg_k
     a = kmin      ! lower limit of integration
     b = kmax      ! upper limit of integration
     epsrel = relacc
     epsabs = absacc
     ans = 0._dp

     call qage (b1g3m_dq,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)

     Pk_b1g3m_prop = 4._dp*pi*ans/fac_norm !remove *matterpowerat(k) to include Pdv later
     if(ier /= 0) call warning('Pk_b1g3m_prop',ier)  

   end function Pk_b1g3m_prop

   function b1g3m_dq(q)        
   implicit none
   real(kind=dp)   :: b1g3m_dq,q,kernel

   kernel  = ((k**2 + q**2)*(3*k**4 - 14*k**2*q**2 + 3*q**4))/(21.*k**2*q**2) + &
             ((k**2 - q**2)**4*log((k - q)**2/(k + q)**2))/(28.*k**3*q**3)
   b1g3m_dq = kernel*matterpowerat(q)

   end function b1g3m_dq

   !! b2 !!

   function Pk_b2(arg_k)
   implicit none
   real(kind=dp)   :: Pk_b2,arg_k,epsabs,epsrel
   real(kind=dp)   :: ans(1),error(1),prob(1)
   integer  :: method

     k = arg_k

     epsabs    = absacc/10._dp
     epsrel    = relacc/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(b2_int,method,epsrel,epsabs,ans,error,prob)
     Pk_b2   = ans(1)/fac_norm

   end function Pk_b2

   subroutine b2_int(ndim,x,ncomp,f)
   implicit   none
   integer    :: ndim,ncomp
   real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,mu12,powq,powp,g2

     q    = kmax*x(1)+kmin*(1._dp-x(1))
     pmin = max(q,dabs(k-q))
     pmax = min(kmax,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        mu12 = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q)
        powp = matterpowerat(p)        
        g2   = 3._dp/7._dp+0.5_dp*mu12*(p/q+q/p)+4._dp/7._dp*mu12**2
        f(1) = (4._dp*pi/k)*q*p*powq*powp*g2*(kmax-kmin)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine b2_int

   !! g2 !!

   function Pk_g2_mc(arg_k)
   implicit none
   real(kind=dp)   :: Pk_g2_mc,arg_k,epsabs,epsrel
   real(kind=dp)   :: ans(1),error(1),prob(1)
   integer  :: method

     k = arg_k

     epsabs    = absacc/10._dp
     epsrel    = relacc/10._dp
     method    = 3 !4:for high lambda, divonne is more accurate than cuhre
     call      integrate2D(g2mc_int,method,epsrel,epsabs,ans,error,prob)
     Pk_g2_mc = ans(1)/fac_norm

   end function Pk_g2_mc

   subroutine g2mc_int(ndim,x,ncomp,f)
   implicit   none
   integer    :: ndim,ncomp
   real(kind=dp)     :: x(ndim),f(ncomp),q,p,pmin,pmax,mu12,kernel,powq,powp,g2

     q    = kmax*x(1)+kmin*(1._dp-x(1))
     pmin = max(q,dabs(k-q))
     pmax = min(kmax,k+q)

     if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        mu12 = (k**2-q**2-p**2)/(2._dp*p*q)
        powq = matterpowerat(q)
        powp = matterpowerat(p)        
        g2   = 3._dp/7._dp+0.5_dp*mu12*(p/q+q/p)+4._dp/7._dp*mu12**2  
        kernel = (mu12**2._dp-1._dp)*g2
        f(1)   = 2._dp*(4._dp*pi/k)*q*p*powq*powp*kernel*(kmax-kmin)*(pmax-pmin)
     else
        f(1)=0._dp
     endif

   end subroutine g2mc_int


   function Pk_g2_prop(arg_k)
     implicit none
     real(kind=dp)   :: Pk_g2_prop,arg_k

     Pk_g2_prop = Pk_b1g2_prop(arg_k)/2._dp

   end function Pk_g2_prop

!! g3m !!

! has been moved to model_power.f90 Feb 5 2015

!! shot-noise renormalization - prop to b2^2       !!
!! i.e. the low-k noise is N --> N + b^2 Pk_noise  !!

   function Pk_noise(arg_k)
   implicit none
   real(kind=dp)   :: Pk_noise,arg_k,epsabs,epsrel,a,b
   integer, parameter :: limit=1000, key = 2 !! local integration rule !!
   real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
   integer  :: neval,ier,last,iord(limit)
   real(kind=dp)   :: ans

     k = arg_k
     a = kmin      ! lower limit of integration
     b = kmax      ! upper limit of integration
     epsrel = relacc
     epsabs = absacc
     ans = 0._dp

     call qage (noise_dq,a,b,epsabs,epsrel,key,limit,ans, &
     abserr,neval,ier,alist,blist,rlist,elist,iord,last)

     Pk_noise = (1._dp/2._dp)*4._dp*pi*ans/fac_norm
     if(ier /= 0) call warning('Pk_noise',ier)  

   end function Pk_noise

   function noise_dq(q)        
     implicit none
     real(kind=dp)   :: noise_dq, q
     noise_dq = matterpowerat(q)**2*q**2
   end function noise_dq

!! warning for 1D integrals !!

   subroutine warning(name,ier)
   implicit   none
   integer    :: ier
   character  :: name*(*)

     if (ier .eq. 1) then
       write(*,*) trim(name)//' achieved Maximum number of subdivisions'
     elseif (ier .eq. 2) then
       write(*,*) trim(name)//' Round off error detected'
     elseif (ier .eq. 3) then
       write(*,*) trim(name)//' Extremely bad integrand behaviour.'
     elseif (ier .eq. 4) then
       write(*,*) trim(name)//' algorithm does not converge. Accur. imposs'
     elseif (ier .eq. 5) then
       write(*,*) trim(name)//' Integrand is probably divergent.'
     elseif (ier .eq. 6) then
       write(*,*) trim(name)//' input is invalid'
     elseif (ier .eq. 7) then
       write(*,*) trim(name)//' bad integrand behaviour occurs in one cycle'
     endif

   end subroutine warning

end module bias
