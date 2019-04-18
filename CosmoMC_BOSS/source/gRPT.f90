! scaling version!
module nonlinear_power
  use settings, only : dp => mcp, feedback
  !use linear_power, only : MatterPowerAt !, sigv, get_sigv
  use numrec
  use quadpack 
  use cubaint
  use rootfinding
  use NLRSD_settings

  implicit none
  logical                             :: init_za = .true.
  real(kind=dp),  private, parameter  :: kmin = 0.0001_dp, kmax = 10._dp
  real(kind=dp),  private, parameter  :: kcut = 100._dp !for I0,2 integration
  real(kind=dp),  private, parameter  :: rcut = 10000._dp !see compute_I0_I2
  real(kind=dp),  private, parameter  :: kroot = 0.1_dp !for gRPT roots: PT or full ZA
  real(kind=dp),  private, parameter  :: relacc = 0.01_dp, absacc = 0.1_dp
!  real(kind=dp),  private, parameter  :: relacc = 0.0001_dp, absacc = 0.0001_dp !for better accuracy
  real(kind=dp),  private, parameter  :: rmin = 0.01_dp, rmax = 20000._dp
  real(kind=dp),  private             :: k
  integer, parameter                  :: nr   = 701
  real(kind=dp)                       :: rvec(nr)
  real(kind=dp),  private             :: I0vec(nr), I0vecd(nr), I2vec(nr), I2vecd(nr)
  real(kind=dp),  private             :: e3sigv2, r
  real(kind=dp)                       :: sigv
  real(kind=dp),  private             :: root_store, a, b, c
  
  contains

    subroutine gRPTpower(arg_k, x, dP13dd, dP13dv, dP13vv, P22inv_dd, P22inv_dv, P22inv_vv)
      implicit none  
      real(kind=dp)   :: arg_k, x, gRPTdd, gRPTdv, gRPTvv
      real(kind=dp)   :: P0, P22inv_dd, P22inv_dv, P22inv_vv, dP13dd, dP13dv, dP13vv

      x      = root(arg_k) !! This call also evaluates sigv !!

      P0     = MatterPowerAt(arg_k)/fac_norm

      P22inv_dd = p22dd(arg_k) - (arg_k*sigv)**2*P0
      dP13dd    = p13dd(arg_k)/P0 + (arg_k*sigv)**2            
      
      P22inv_dv = p22dv(arg_k) - (arg_k*sigv)**2*P0
      dP13dv    = p13dv(arg_k)/P0 + (arg_k*sigv)**2            

      P22inv_vv = p22vv(arg_k) - (arg_k*sigv)**2*P0
      dP13vv    = p13vv(arg_k)/P0 + (arg_k*sigv)**2            

    end subroutine gRPTpower

    subroutine gRPTpower_22inv(arg_k, x, P22inv_dd, P22inv_dv, P22inv_vv)
      implicit none  
      real(kind=dp)   :: arg_k, x, gRPTdd, gRPTdv, gRPTvv
      real(kind=dp)   :: P0, P22inv_dd, P22inv_dv, P22inv_vv

      x      = root(arg_k) !! This call also evaluates sigv !!
      P0     = MatterPowerAt(arg_k)/fac_norm
      P22inv_dd = p22dd(arg_k) - (arg_k*sigv)**2*P0
      P22inv_dv = p22dv(arg_k) - (arg_k*sigv)**2*P0
      P22inv_vv = p22vv(arg_k) - (arg_k*sigv)**2*P0

    end subroutine gRPTpower_22inv

    subroutine gRPTpower_13(arg_k, dP13dd, dP13dv, dP13vv)
      implicit none  
      real(kind=dp)   :: arg_k, x, gRPTdd, gRPTdv, gRPTvv
      real(kind=dp)   :: P0, dP13dd, dP13dv, dP13vv

      P0     = MatterPowerAt(arg_k)/fac_norm
      dP13dd    = p13dd(arg_k)/P0 + (arg_k*sigv)**2            
      dP13dv    = p13dv(arg_k)/P0 + (arg_k*sigv)**2            
      dP13vv    = p13vv(arg_k)/P0 + (arg_k*sigv)**2            

    end subroutine gRPTpower_13


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    Boost                !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function root(arg_k)
      implicit none
      real(kind=dp)   :: root,arg_k,disc
      real(kind=dp)   :: p0L,p1L,p2L,pza
      real(kind=dp)   :: xinit
      real(kind=dp)   :: tol = 1.0e-6_dp
      integer         :: maxiter = 200
      logical         :: success

      call  pk_za(arg_k,p0L,p1L,p2L,pza)

      if(arg_k <= kroot) then
        b  = -2._dp*p1L/p0L
        c  =  2._dp*p2L/p0L
        disc = (b/2._dp)**2 - c
        if(disc < 0._dp) disc = 0._dp
        root = -b/2._dp - dsqrt(disc)
      else
        a     = pza
        b     = p0L 
        c     = p1L
        xinit = root_store 
        !xinit = -0.5_dp
        call find_root( gRPTza_1L , xinit, tol, maxiter, root, success )
      endif
      
      root_store = root 

    end function root


    function gRPTza_1L( x )
      implicit none
      real(kind=dp)   :: gRPTza_1L
      real(kind=dp),  intent(in) :: x
   
      gRPTza_1L = exp(x)*(b*(1._dp-x)+c)-a      
            
    end function gRPTza_1L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    Zeldovich Power Spectrum     !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine pk_za(ak, p0L, p1L, p2L, pza)
      implicit   none
      real(kind=dp)     :: ak, p0L, p1L, p2L, pza
      real(kind=dp)     :: I1k, I2k, I3k, x

      if(init_za)then
        init_za = .false.
        call compute_I0_I2()
      end if
       
      I1k = MatterPowerAt(ak)/fac_norm
      I2k = p1loop(ak)    !! Mode-coupling loops !!
      I3k = p2loops(ak)

      x  = (ak*sigv)**2

      p0L = I1k           !! Loop orders !! 
      p1L = I2k - I1k*x
      p2L = I3k - I2k*x + I1k*x**2/2._dp

      pza = pza_nloops(ak)
           
    end subroutine pk_za

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Zeldovich Generating Functions !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_sigv(sigmav)
        implicit none
        integer        :: iwr
        real(kind=dp)  :: tol, sigmav
        tol     = 1.e-4_dp                                                    
        iwr     = 0          
        call gabq(MatterPowerAt, kmin, kmax, sigmav, tol, iwr)
        sigmav   = sqrt(4._dp*pi*sigmav/3._dp/fac_norm)
      end subroutine get_sigv

    subroutine compute_I0_I2()
      implicit   none
      integer           :: ir, iwr
      logical           :: log_binning
      real(kind=dp)     :: eps, ar, tol
 
      eps = kmin
 
      !sigv is computed in module linear_power
      call get_sigv(sigv)
      if(feedback > 1)write(*,'(A,F10.5)')'sigv = ',sigv
      e3sigv2 = 3._dp*sigv**2
 
 
      log_binning = .true.
 
      if(log_binning)then
        do ir = 1,nr
          rvec(ir) = 10._dp**(log10(rmax/rmin)*dble(ir-1)/dble(nr-1)+log10(rmin))
        enddo
      else
        do ir = 1,nr
          rvec(ir) = (rmax-rmin)*dble(ir-1)/dble(nr-1)+rmin
        enddo
      endif
 
      do ir = 1,nr
        ar = rvec(ir)
        if(feedback > 1 .and. mod(ir,100) == 0) write(*,*)ir,ar
        !revisar esta parte!!!! ARIEL
        !revisar esta parte!!!! ARIEL
        if(ar >= rcut) then 
          I0vec(ir) = 0._dp
          I2vec(ir) = 0._dp
        else
          I0vec(ir) = I0fun(ar,eps)
          I2vec(ir) = Ithfun(ar,eps) + I0vec(ir)
        endif
      enddo
 
      call spline(rvec,I0vec,nr,3d30,3d30,I0vecd)
      call spline(rvec,I2vec,nr,3d30,3d30,I2vecd)
      
      if(feedback > 1)write(*,'(A,2E16.5)')'> Done with I0 and I2, limits of r variable',rmin,rmax
 
    end subroutine compute_I0_I2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!       Compute and evaluate I0              !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function I0(RR)   
      implicit none
      real(kind=dp)   :: I0,RR

      if (RR >= rmin .and. RR <= rmax) then 
        call splint(rvec,I0vec,I0vecd,nr,RR,I0)
      elseif (RR < rmin) then
        I0 = e3sigv2
      elseif (RR > rmax) then
        I0 = 0._dp
      endif

    end function I0

    function I0fun(arg_r, eps)   
      implicit none
      real(kind=dp)                    :: I0fun, arg_r, eps
      integer, parameter               :: limit = 10000, limlst = 5000, maxp1 = 5000
      integer                          :: neval, ier
      integer, dimension(limit)        :: iord, nnlog
      integer, dimension(limlst)       :: ierlst
      real(kind=dp), dimension(limit)  :: alist, blist, elist, rlist
      real(kind=dp), dimension(limlst) :: rslst, erlst
      real(kind=dp)                    :: epsabs, abserr, ans
      integer                          :: integr, lst
      real(kind=dp)                    :: chebmo(maxp1,25)

      r      = arg_r
      epsabs = 1.e-7_dp         !! Relative Error !!
      integr = 2
      
      call qawfe(I0q,eps,r,integr,epsabs,limlst,limit,maxp1,      &
      ans,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist,    &
      rlist,elist,iord,nnlog,chebmo)     
      I0fun = (twopi/r)*ans/fac_norm
      !!if(ier.gt.0_dp) print*,'error in', r, abserr/ans

    end function I0fun
      
    function I0q(q)   
      implicit none
      real(kind=dp)   :: I0q, q
      I0q = 2._dp*MatterPowerAt(q)/q
    end function I0q

    function I0fun_slow(arg_r,eps)   
      implicit none
      real(kind=dp)      :: I0fun_slow,arg_r,eps
      integer, parameter :: limit = 5000, key = 3
      integer            :: neval,ier,iord(limit),last
      real(kind=dp)      :: alist(limit), blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans

      r      = arg_r
      epsabs = 0._dp
      epsrel = 0.1e-7_dp         !! Relative Error !!

      call qage(I0q_slow, eps, kcut, epsabs, epsrel, key, limit,   &
         ans, abserr, neval, ier, alist, blist, rlist, elist, iord, last)    
      I0fun_slow = twopi*ans/fac_norm

    end function I0fun_slow

    function I0q_slow(q)   
      implicit none
      real(kind=dp)   :: I0q_slow,q, a, ff

      a   = q*r

      if(a >= 0.001_dp .and. a <= 1000._dp)then
        f = 2._dp*sin(a)/a
      elseif(a < 0.001_dp)then
        f = 2._dp - a**2/3._dp
      elseif(a > 1000._dp)then
        f = 2._dp*sin(a)/a
      endif
     
      I0q_slow = MatterPowerAt(q)*f

    end function I0q_slow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!       Compute and evaluate I2              !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function I2(RR)   
      implicit none
      real(kind=dp)   :: I2, RR
    
      if (RR >= rmin .and. RR <= rmax) then 
        call splint(rvec,I2vec,I2vecd,nr,RR,I2)
      else
        I2 = 0._dp
      endif

    end function I2

    function I2fun(arg_r,eps)      
      implicit none
      real(kind=dp)      :: I2fun,arg_r,eps
      integer, parameter :: limit=5000, key = 3
      integer            :: neval,ier,iord(limit),last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs,epsrel,abserr,ans

      r      = arg_r
      epsabs = 0._dp
      epsrel = 1.e-7_dp         !! Relative Error !!
 
      call qage(I2q, eps, kcut, epsabs, epsrel, key, limit, &
      ans, abserr, neval, ier, alist, blist, rlist, elist, iord, last)    
      I2fun = twopi*ans/fac_norm
      !! if(ier.gt.0_dp) print*,'error in', r, abserr/ans
 
    end function I2fun
      
    function I2q(q)   
      implicit none
      real(kind=dp)   :: I2q, q, a, ff

      a   = q*r
      
      if(a >= 0.001_dp .and. a <= 1000._dp)then
        ff = 6._dp*cos(a)/a**2._dp - 6._dp*sin(a)/a**3._dp + 2._dp*sin(a)/a
      elseif(a < 0.001_dp)then
        ff = -2._dp*a**2/15._dp+a**4/105._dp-a**6/3780._dp
      elseif(a > 1000._dp)then
        ff = (2._dp*sin(a))/a
      endif

      I2q = MatterPowerAt(q)*ff

    end function I2q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!       Compute and evaluate I top-hat       !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function Ithfun(arg_r,eps)   
      implicit none
      real(kind=dp)   :: Ithfun,arg_r,eps
      integer, parameter :: limit=10000, key = 3
      integer  :: neval,ier,iord(limit),last
      real(kind=dp)   :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)   :: epsabs,epsrel,abserr,ans
 
      r      = arg_r
      epsabs = 0._dp
      epsrel = 1.e-7_dp         !! Relative Error !!
 
      call qage(Ithq,eps,kcut,epsabs,epsrel, key, limit,             &
      ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      Ithfun = twopi*ans/fac_norm
      !!if(ier.gt.0_dp) print*,'error in', r, abserr/ans
    end function Ithfun
      
    function Ithq(q)   
      implicit none
      real(kind=dp)   :: Ithq, q, a, ff
      a   = q*r
      if(a >= 0.001_dp .and. a <= 1000._dp)then
        ff = 6._dp*cos(a)/a**2 - 6._dp*sin(a)/a**3 
      elseif(a < 0.001_dp)then
        ff = -2._dp + a**2/5._dp - a**4/140._dp
      elseif(a > 1000._dp)then
        ff = (2._dp*cos(a))/a**2
      endif
      Ithq = MatterPowerAt(q)*ff

    end function Ithq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Mode-Coupling Loops of the Zeldovich Power Spectrum !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function p1loop(arg_k)   
      implicit none
      real(kind=dp)   :: p1loop,arg_k 
      integer, parameter :: limit = 5000  
      integer  :: neval,ier,iord(limit),last
      real(kind=dp)   :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)   :: epsabs,epsrel,abserr,ans
       
      k      = arg_k
      epsabs = 0._dp
      epsrel = 0.00001_dp         !! Relative Error !!
 
      call qage(p1loop_r,rmin,rmax,epsabs,epsrel,3,limit,    &
      ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
 
      p1loop = (1._dp/2._dp)*ans*(k**2/3._dp)**2/(twopi)**2
 
    end function p1loop
      
    function p1loop_r(r)   
      implicit none
      real(kind=dp)   :: p1loop_r, r, ans, I2r, I0r, kr
  
      I2r = I2(r)
      I0r = I0(r)
      kr = k*r
  
      ans = (24._dp*I2r*kr*(-18._dp*I2r + (I0r + 2._dp*I2r)*kr**2   &
          )*dcos(kr) + 2._dp*(216._dp*I2r**2 - 12._dp*I2r*(I0r + 8._dp*I2r   &
          )*kr**2 + (I0r + 2._dp*I2r)**2*kr**4)*dsin(kr))  &
          /(kr**5)

      p1loop_r = ans*r**2
   
    end function p1loop_r
   
    function p2loops(arg_k)   
      implicit none
      real(kind=dp)      :: p2loops,arg_k 
      integer, parameter :: limit=1000
      integer            :: neval,ier,iord(limit),last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs,epsrel,abserr,ans
   
      k      = arg_k
      epsabs = 0._dp
      epsrel = 0.00001_dp         !! Relative Error !!
  
      call  qage(p2loops_r,rmin,rmax,epsabs,epsrel,3,limit,   &
      ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
  
      p2loops = (1._dp/6._dp)*ans*(k**2/3._dp)**3/(twopi)**2
  
    end function p2loops
       
    function p2loops_r(r)   
      implicit none
      real(kind=dp)   :: p2loops_r,r,ans,I2r,I0r
  
      I2r = I2(r)
      I0r = I0(r)
  
      ans = (36._dp*I2r*k*r*(1080._dp*I2r**2._dp - 36._dp*I2r*(I0r +       &
          4._dp*I2r)*k**2._dp*r**2._dp + (I0r + 2._dp*I2r)**2._dp*k**4._dp    &
          *r**4._dp)*dcos(k*r) + 2._dp*(-19440._dp*I2r**3._dp + 648._dp*I2r  &
          **2._dp*(I0r + 14._dp*I2r)*k**2._dp*r**2._dp - 18._dp*I2r*(I0r     &  
          + 2._dp*I2r)*(I0r + 14._dp*I2r)*k**4._dp*r**4._dp + (I0r + 2._dp   &
          *I2r)**3._dp*k**6._dp*r**6._dp)*dsin(k*r))/(k**7._dp*r**7._dp)

!     check this!!!!! ARIEL
!<<<<<<< HEAD
      p2loops_r = ans*r**2
!=======
!     p1loop = (1._dp/2._dp)*ans*(k**2/3._dp)**2/(twopi)**2
!>>>>>>> master

    end function p2loops_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ZA nonlinear power spectrum non-perturbatively !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

    function pza_nloops(arg_k)   
      implicit none
      real(kind=dp)      :: pza_nloops, arg_k
      integer, parameter :: limit=100000, key=3
      integer            :: neval, ier, iord(limit), last
      real(kind=dp)      :: alist(limit),blist(limit),elist(limit),rlist(limit)
      real(kind=dp)      :: epsabs, epsrel, abserr, ans
 
      k      = arg_k
      epsabs = 0_dp
      epsrel = 0.000001_dp         !! Relative Error !!
 
      call qage(pnloops_r,rmin,rmax,epsabs,epsrel,key,limit,    &
      ans,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        
      pza_nloops = exp(-(k*sigv)**2)*ans/twopisqr
 
    end function pza_nloops
    
    function pnloops_r(arg_r)   
      implicit none
      real(kind=dp)   :: pnloops_r, arg_r
      real(kind=dp)   :: ans, arg, pre, I2r, I0r, x
 
      I2r = I2(arg_r)
      I0r = I0(arg_r)
        
      x   = -2._dp*k*I2r/arg_r
      pre = exp(k**2*(I0r + 2._dp*I2r)/3._dp)
 
      arg = k*arg_r
 
      ans = x*sphbe1(arg) + x**2*sphbe2(arg) + x**3*sphbe3(arg) + &
            x**4*sphbe4(arg) + x**5*sphbe5(arg) + x**6*sphbe6(arg) + &
            x**7*sphbe7(arg) + x**8*sphbe8(arg)
      pnloops_r = ans

      ans = (pre-1._dp)*sphbe0(arg) + pre*ans
 
      pnloops_r = ans*arg_r**2
 
    end function pnloops_r

    function sphbe0(x) !added small arg expansion
      implicit none
      real(kind=dp)   :: sphbe0,x
      if (x.gt.0.001_dp) then
         sphbe0 = dsin(x)/x 
      else
         sphbe0 = 1._dp - x**2/6._dp !O(1.e-12)
      endif   
    end function sphbe0

    function sphbe1(x) !added small arg expansion
      implicit none
      real(kind=dp)   :: sphbe1,x
         if (x.lt.0.001_dp) then
            sphbe1 = x/3._dp - x**3/30._dp !O(1.e-15)
         else
            sphbe1 = (-(x*dcos(x)) + dsin(x))/x**2
         endif   
    end function sphbe1

    function sphbe2(x) !added small arg expansion
      implicit none
      real(kind=dp)   :: sphbe2,x
         if (x.lt.0.01_dp) then
            sphbe2 = x**2/15._dp - x**4/210._dp !O(1.e-12)
         else
            sphbe2 = -((3._dp*x*dcos(x) + (-3._dp + x**2)*dsin(x))/x**3)
         endif   
    end function sphbe2

    function sphbe3(x) !added small arg expansion
      implicit none
      real(kind=dp)   :: sphbe3,x
         if (x.lt.0.01_dp) then
            sphbe3 = x**3/105._dp - x**5/1890._dp !O(1.e-14)
         else      
            sphbe3 = (x*(-15._dp + x**2)*dcos(x) + 3._dp*(5._dp - 2._dp*x**2 &
                     )*dsin(x))/x**4
         endif  
    end function sphbe3

    function sphbe4(x) !added small arg expansion
      implicit none
      real(kind=dp)   :: sphbe4,x
         if (x.lt.0.01_dp) then
             sphbe4 = x**4._dp/945._dp !O(1.e-12)
          else
             sphbe4 = (5._dp*x*(-21._dp + 2._dp*x**2)*dcos(x) + (105._dp - 45._dp & 
                      *x**2 + x**4)*dsin(x))/x**5
          endif 
    end function sphbe4

    function sphbe5(x)
      implicit none
      real(kind=dp)   :: sphbe5,x
      sphbe5 = (-(x*(945._dp - 105._dp*x**2 + x**4._dp)*dcos(x)) + 15._dp*( &
           63._dp - 28._dp*x**2 + x**4)*dsin(x))/x**6            
    end function sphbe5

    function sphbe6(x)
      implicit none
      real(kind=dp)   :: sphbe6,x
      sphbe6 = -((21._dp*x*(495._dp - 60._dp*x**2 + x**4)*dcos(x) +   &
           (-10395._dp + 4725._dp*x**2 - 210._dp*x**4 + x**6)*   &
           dsin(x))/x**7)
    end function sphbe6

    function sphbe7(x)
      implicit none
      real(kind=dp)   :: sphbe7,x
      sphbe7 = (x*(-135135._dp + 17325._dp*x**2 - 378._dp*x**4 +      &
           x**6)*dcos(x) - 7._dp*(-19305._dp + 8910._dp*x**2 - 450._dp* &  
           x**4 + 4._dp*x**6)*dsin(x))/x**8
    end function sphbe7

    function sphbe8(x)
      implicit none
      real(kind=dp)   :: sphbe8,x
      sphbe8 = (9._dp*x*(-225225._dp + 30030._dp*x**2 - 770._dp*x**4 + &
             4._dp*x**6)*dcos(x) + (2027025._dp - 945945._dp*x**2 +   &
             51975._dp*x**4 - 630._dp*x**6 + x**8)*dsin(x))       &
             /x**9
    end function sphbe8

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Standard Pertubation Theory !! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! P22 density-density !!   

    function p22dd(arg_k)
      implicit none
      real(kind=dp)   :: p22dd,arg_k,epsabs,epsrel
      real(kind=dp)   :: ans(1),error(1),prob(1)
      integer  :: method
 
      k = arg_k
 
      epsabs = absacc/100._dp
      epsrel = relacc/100._dp
      method = 3 !4:for high lambda, divonne is more accurate than cuhre
      call   integrate2D(p22dd_int,method,epsrel,epsabs,ans,error,prob)
      p22dd  = ans(1)/fac_norm**2
 
    end function p22dd
 
    subroutine p22dd_int(ndim,x,ncomp,f)
      implicit   none
      integer    :: ndim,ncomp
      real(kind=dp)     :: x(ndim),f(ncomp),q,pmin,pmax,y,f2,cos,p
 
      q    = kmax*x(1)+kmin*(1._dp-x(1))
      pmin = max(q,dabs(k-q))
      pmax = min(kmax,k+q)
 
      if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        cos  = (k**2-q**2-p**2)/(2._dp*p*q)
        f2   = 5._dp/7._dp+0.5_dp*cos*(p/q+q/p)+2._dp/7._dp*cos**2
        y    = MatterPowerAt(p)*MatterPowerAt(q)
        f(1) = (4._dp*pi/k)*q*p*y*2._dp*f2**2*(kmax-kmin)*(pmax-pmin)
      else
        f(1) = 0._dp
      endif
 
    end subroutine p22dd_int
 
    !!P22 density-velocity !!   
 
    function p22dv(arg_k)
      implicit none
      real(kind=dp)   :: p22dv,arg_k,epsabs,epsrel
      real(kind=dp)   :: ans(1),error(1),prob(1)
      integer  :: method
 
      k = arg_k
 
      epsabs = absacc/100._dp
      epsrel = relacc/100._dp
      method = 3 !4:for high lambda, divonne is more accurate than cuhre
      call   integrate2D(p22dv_int, method, epsrel, epsabs, ans, error, prob)
      p22dv  = ans(1)/fac_norm**2
 
    end function p22dv
 
    subroutine p22dv_int(ndim,x,ncomp,f)
      implicit none
      integer           :: ndim,ncomp
      real(kind=dp)     :: q, pmin, pmax, y, f2, g2, vcos, p
      real(kind=dp)     :: term, term2
      real(kind=dp)     :: x(ndim), f(ncomp)
 
      q    = kmax*x(1)+ kmin*(1._dp - x(1))
      pmin = max(q,dabs(k-q))
      pmax = min(kmax,k+q)
 
      if (pmin < pmax) then
        p     = pmax*x(2)+pmin*(1._dp-x(2))
        vcos  = (k**2-q**2-p**2)/(2._dp*p*q)
        term  = 0.5_dp*vcos*(p/q+q/p) 
        term2 = 2._dp*vcos**2/7._dp
        f2    = 5._dp/7._dp + term + term2
        g2    = 3._dp/7._dp + term + 2._dp*term2
        y     = MatterPowerAt(p)*MatterPowerAt(q)
        f(1)  = (4._dp*pi/k)*q*p*y*2._dp*f2*g2*(kmax-kmin)*(pmax-pmin)
      else
        f(1) = 0._dp
      endif
 
    end subroutine p22dv_int
 
    !! P22 velocity-velocity !!   
 
    function p22vv(arg_k)
      implicit none
      real(kind=dp)   :: p22vv,arg_k,epsabs,epsrel
      real(kind=dp)   :: ans(1),error(1),prob(1)
      integer  :: method
 
      k = arg_k
 
      epsabs = absacc/100._dp
      epsrel = relacc/100._dp
      method = 3 !4:for high lambda, divonne is more accurate than cuhre
      call   integrate2D(p22vv_int,method,epsrel,epsabs,ans,error,prob)
      p22vv  = ans(1)/fac_norm**2
 
    end function p22vv
 
    subroutine p22vv_int(ndim,x,ncomp,f)
      implicit   none
      integer           :: ndim,ncomp
      real(kind=dp)     :: x(ndim),f(ncomp),q,pmin,pmax,y,g2,vcos,p
 
      q    = kmax*x(1)+kmin*(1._dp-x(1))
      pmin = max(q,dabs(k-q))
      pmax = min(kmax,k+q)
 
      if (pmin.lt.pmax) then
        p    = pmax*x(2)+pmin*(1._dp-x(2))
        vcos  = (k**2 - q**2 - p**2)/(2._dp*p*q)
        g2   = 3._dp/7._dp + 0.5_dp*vcos*(p/q+q/p)+4._dp/7._dp*vcos**2
        y    = MatterPowerAt(p)*MatterPowerAt(q)
        f(1) = (4._dp*pi/k)*q*p*y*2._dp*g2**2*(kmax-kmin)*(pmax-pmin)
      else
        f(1) = 0._dp
      endif
 
    end subroutine p22vv_int
 
    !! P13 density-density !!   
 
    function p13dd(arg_k)
      implicit none
      real(kind=dp)      :: p13dd, arg_k, epsabs, epsrel, aa, bb
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)      :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer            :: neval, ier, last, iord(limit)
      real(kind=dp)      :: ans
 
      k  = arg_k
      aa = kmin      ! lower limit of integration
      bb = kmax      ! upper limit of integration
      epsrel = relacc/10._dp
      epsabs = absacc/10_dp
      ans = 0._dp
 
      call qage (ps13dd_dq, aa, bb, epsabs, epsrel, key, limit, ans, &
                 abserr, neval, ier, alist, blist, rlist, elist, iord, last)
      p13dd = ans
      if(ier /= 0) call warning('p13dd',ier)
 
 
    end function p13dd
 
    function ps13dd_dq(q)        
      implicit none
      real(kind=dp)   :: ps13dd_dq, q, asum
 
      asum = (6._dp*k**6 - 79._dp*k**4*q**2 + 50._dp*k**2*q**4 - 21._dp*q**6)/    & 
            (63._dp*k**2*q**2) +(q**2 - k**2)**3*(7._dp*q**2 + 2._dp*k**2)/   &
            (42._dp*k**3*q**3)*log((k+q)/abs(k-q))
 
      ps13dd_dq = pi*asum*MatterPowerAt(q)*MatterPowerAt(k)/fac_norm**2
 
    end function ps13dd_dq
 
!!! P13 density - velocity !!   
 
    function p13dv(arg_k)
      implicit none
      real(kind=dp)      :: p13dv, arg_k, epsabs, epsrel, aa, bb
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)      :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer            :: neval,ier,last,iord(limit)
      real(kind=dp)      :: ans
 
      k  = arg_k
      aa = kmin      ! lower limit of integration
      bb = kmax      ! upper limit of integration
      epsrel = relacc/10._dp
      epsabs = absacc/10._dp
      ans = 0._dp
 
      call qage (ps13dv_dq, aa, bb, epsabs, epsrel, key, limit, ans, &
              abserr, neval, ier, alist, blist, rlist, elist, iord, last)
      p13dv = ans
      if(ier /= 0) call warning('p13dv',ier)
 
    end function p13dv
 
    function ps13dv_dq(q)        
      implicit none
      real(kind=dp)   :: ps13dv_dq, q, asum
 
      asum  = (12._dp*k**6 - 101._dp*k**4*q**2 + 28._dp*k**2*q**4 - 15._dp*q**6)/ &
           (63._dp*k**2*q**2) +(q**2 - k**2)**3*(5._dp*q**2 + 4._dp*k**2)/    &
           (42._dp*k**3*q**3)*log((k+q)/abs(k-q))
 
      ps13dv_dq = pi*asum*MatterPowerAt(q)*MatterPowerAt(k)/fac_norm**2
 
    end function ps13dv_dq

   !! P13 velocity - velocity !!   
 
    function p13vv(arg_k)
      implicit none
      real(kind=dp)   :: p13vv, arg_k, epsabs, epsrel, aa, bb
      integer, parameter :: limit=1000, key = 2 !! local integration rule !!
      real(kind=dp)   :: abserr,alist(limit),blist(limit),rlist(limit),elist(limit)
      integer  :: neval,ier,last,iord(limit)
      real(kind=dp)   :: ans
 
      k = arg_k
      aa = kmin      ! lower limit of integration
      bb = kmax      ! upper limit of integration
      epsrel = relacc/10._dp
      epsabs = absacc/10._dp
      ans = 0._dp
 
      call qage (ps13vv_dq, aa, bb, epsabs, epsrel, key, limit, ans, &
                abserr, neval, ier, alist, blist, rlist, elist, iord, last)
      p13vv = ans
      if(ier /= 0) call warning('p13vv',ier)
 
    end function p13vv
 
    function ps13vv_dq(q)        
      implicit none
      real(kind=dp)   :: ps13vv_dq,q,asum
 
      asum  = (6._dp*k**6 - 41._dp*k**4*q**2 + 2._dp*k**2*q**4 - 3._dp*q**6)/    &
           (21._dp*k**2*q**2) + (q**2 - k**2)**3*(1._dp*q**2 + 2._dp*k**2)/   &
           (14._dp*k**3*q**3)*log((k + q)/abs(k - q))
 
      ps13vv_dq = pi*asum*MatterPowerAt(q)*MatterPowerAt(k)/fac_norm**2
 
    end function ps13vv_dq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Standard Pertubation Theory !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function ps13_dq(x)        
      implicit none
      real(kind=dp)   :: ps13_dq, x, y

      y = (pi*(4._dp*(6._dp*k**7*x - 79._dp*k**5*x**3 + 50._dp*k**3*x**5 - 21._dp*k*x**7) +  &
          3._dp*(k**2 - x**2)**3*(2._dp*k**2 + 7._dp*x**2)*(log((k - x)**2) -        &
          log((k + x)**2))))/(1512._dp*k**3*x**3)
 
      ps13_dq = 6._dp*y*MatterPowerAt(x)*MatterPowerAt(k)

    end function ps13_dq
    !! warning for 1D integrations !!
 
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
       write(*,*) trim(name)//' bad integrand occurs in one cycle/abnormal termination'
     elseif (ier .eq. 8) then
       write(*,*) trim(name)//' maximum number of cycles allowed has been achieved'
     elseif (ier .eq. 9) then
       write(*,*) trim(name)//' the extrapol table for convergence does not converge'
     endif

   end subroutine warning

end module nonlinear_power
