module m_costf
use mod_inistat
use m_func
contains

real function jxfunc(x)
   implicit none
   real, intent(in) :: x
   jxfunc=(x-x0)**2/siga**2 + (func(x)-d)**2/sigo**2
end function jxfunc

real function djxfunc(x)
   implicit none
   real, intent(in) :: x
   djxfunc=sigo**2*(x-x0) + siga**2*(func(x)-d)*dfunc(x)
end function djxfunc

real function ddjxfunc(x)
   implicit none
   real, intent(in) :: x
   ddjxfunc=sigo**2 + siga**2*dfunc(x)**2
end function ddjxfunc

end module m_costf
