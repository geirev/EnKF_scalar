module m_iescostf
use mod_inistat
use m_func
contains

real function jxfunc(x,xf,dpert)
   implicit none
   real, intent(in) :: x
   real, intent(in) :: xf
   real, intent(in) :: dpert
   jxfunc=0.5*(x-xf)**2/siga**2 + 0.5*(func(x,0.0)-dpert)**2/sigo**2
end function jxfunc

real function djxfunc(x,xf,dpert)
   implicit none
   real, intent(in) :: x
   real, intent(in) :: xf
   real, intent(in) :: dpert
   djxfunc=(x-xf)/siga**2 + (func(x,0.0)-dpert)*dfunc(x,0.0)/sigo**2
   djxfunc=cdd*(x-xf) + (func(x,0.0)-dpert)*dfunc(x,0.0)*siga**2
end function djxfunc

real function ddjxfunc(x,xf,dpert)
   implicit none
   real, intent(in) :: x
   real, intent(in) :: xf
   real, intent(in) :: dpert
   ddjxfunc=1.0/siga**2 + dfunc(x,0.0)**2/sigo**2
end function ddjxfunc

end module m_iescostf
