module m_func
! Supports two functions
!  y = x + beta*x^3 +q
!  y = 1 + sin(x*pi)
   use mod_inistat
   implicit none
   real, parameter :: pi=3.141592654
   real, save ::  beta               ! model parameter
   integer, save ::  funcmode        ! which function to use

contains

real function func(x,q)
   real, intent(in) :: x
   real, intent(in) :: q
   select case (funcmode)
   case(0)
      func=x+beta*x**3 + q
   case(1)
      func=1.0+sin(x*pi) + q
   case default
      stop 'func: invalid mode for function func'
   end select
end function func

function dfunc(x,q)
! Analytic gradients used for adjoint sensitivities 
   real dfunc(2)
   real, intent(in) :: x
   real, intent(in) :: q
   select case (funcmode)
   case(0)
      dfunc(1)=1.0+3.0*beta*x**2  ! derivative wrt x
      dfunc(2)=1.0                ! derivative wrt q
   case(1)
      dfunc(1)=pi*cos(x*pi)       ! derivarive wrt x
      dfunc(2)=1.0                ! derivative wrt q
   case default
      stop 'dfunc: invalid mode for function dfunc'
   end select
end function dfunc

function ddfunc(x,q)
! Analytic Hessian 
! Only computes ddf/ddx since ddf/ddq = 0
   real ddfunc(2,2)
   real, intent(in) :: x
   real, intent(in) :: q
   select case (funcmode)
   case(0)
      ddfunc(1,1)=6.0*beta*x
      ddfunc(2,2)=0.0
      ddfunc(1,2)=0.0
      ddfunc(2,1)=0.0
   case(1)
      ddfunc=-pi**2*sin(x*pi)
      ddfunc(2,2)=0.0
      ddfunc(1,2)=0.0
      ddfunc(2,1)=0.0
   case default
      stop 'ddfunc: invalid mode for function ddfunc'
   end select
end function ddfunc
end module m_func
