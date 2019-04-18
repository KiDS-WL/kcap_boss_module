module rootfinding
  use settings, only: dp => mcp 
  implicit none  

  contains

  subroutine find_root( f, xinit, tol, maxiter, result, success )

    interface
       function f(x)
         real(kind=8)             :: f
         real(kind=8), intent(in) :: x
       end function f
    end interface

    real(kind=dp),  intent(in)  :: xinit, tol
    integer, intent(in)         :: maxiter
    real(kind=dp),  intent(out) :: result
    logical, intent(out)        :: success

    real(kind=dp), parameter    :: eps = 1.0e-3_dp
    real(kind=dp)               :: fx1,fx2,fprime,x,xnew
    integer                     :: i

    result  = 0._dp
    success = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!!    Newton-Raphson method     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    x = xinit
    do i = 1,max(1,maxiter)

      fx1    = f(x)
      fx2    = f(x+eps)
      fprime = (fx2 - fx1) / eps
      xnew   = x - fx1 / fprime

      if ( abs(xnew-x) <= tol ) then
          success = .true.
          result  = xnew
          !write(*,*) i, f(xnew)
          exit
      endif

      x = xnew
      
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!!      Secant method        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$    x0 = xinit + 0.03d0
!!$    x1 = xinit - 0.03d0
!!$    do i = 1,max(1,maxiter)
!!$
!!$        DX   = x1-x0
!!$        DF   = f(x1) - f(x0)
!!$        x2   = x1 - f(x1)*DX/DF
!!$
!!$        if ( abs(x2-x1) <= tol ) then
!!$            success = .true.
!!$            result  = x2
!!$            !write(*,*) i, f(x2)
!!$            exit
!!$        endif
!!$
!!$        x0   = x1
!!$        x1   = x2
!!$        
!!$     enddo

  end subroutine find_root

end module
