module m_func
   use mod_inistat
   implicit none
   real, parameter :: pi=3.141592654
   real, save ::  beta    ! model parameter
   integer, save ::  funcmode ! which function to use

contains

real function func(x,q)
   real, intent(in) :: x
   real, intent(in) :: q
   select case (funcmode)
   case(0)
      func=x+beta*x**3 + q 
   case(1)
      func=x + beta*x**3 + q + beta*q**3
!      func=func-0.5*(func-d-q)
   case(2)
      func=1.0+sin(x*pi) + q
   case default
      stop 'func: invalid mode for function func'
   end select
end function func

real function dfunc(x,q)
   real, intent(in) :: x
   real, intent(in) :: q
   select case (funcmode)
   case(0)
      dfunc=1.0+3.0*beta*x**2
   case(1)
      dfunc=beta*(3.0*x**2-1.0)
   case(2)
      dfunc=pi*cos(x*pi)
   case default
      stop 'dfunc: invalid mode for function dfunc'
   end select
end function dfunc

real function ddfunc(x,q)
   real, intent(in) :: x
   real, intent(in) :: q
   select case (funcmode)
   case(0)
      ddfunc=6.0*beta*x
   case(1)
      ddfunc=beta*(6.0*x)
   case(2)
      ddfunc=-pi**2*sin(x*pi)
   case default
      stop 'ddfunc: invalid mode for function ddfunc'
   end select
end function ddfunc
end module m_func
