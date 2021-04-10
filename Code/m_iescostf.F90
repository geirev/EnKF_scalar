module m_iescostf
use mod_inistat
use m_func
contains

real function costf(x,xf,q,dpert)
   implicit none
   real, intent(in) :: x
   real, intent(in) :: q
   real, intent(in) :: xf
   real, intent(in) :: dpert
   costf=0.5*(x-xf)**2/siga**2 + 0.5*q**2/sigw**2 + 0.5*(func(x,q)-dpert)**2/sigo**2
end function costf

function dcost(x,xf,q,dpert)
   implicit none
   real dcost(2)
   real, intent(in) :: x
   real, intent(in) :: q
   real, intent(in) :: xf
   real, intent(in) :: dpert
   real             :: dg(2)
   dg=dfunc(x,q)
   dcost(1)=(x-xf)/siga**2 + (func(x,q)-dpert)*dg(1)/sigo**2
   dcost(2)=q/sigw**2 + (func(x,q)-dpert)*dg(2)/sigo**2  
end function dcost

function ddjfunc(x,xf,q,dpert)
   implicit none
   real ddjfunc(2,2)
   real, intent(in) :: x
   real, intent(in) :: q
   real, intent(in) :: xf
   real, intent(in) :: dpert
   ddjfunc=0.0
 !  ddjxfunc=1.0/siga**2 + dfunc(x,0.0)**2/sigo**2
end function ddjfunc

end module m_iescostf
